/*
    Tcl initialization and main program
    Copyright (C) 1999 Mads Bak, 2000-2005 Thomas Vosegaard
                  2010 Zdenek Tosner, Rasmus Andersen, Thomas Vosegaard

    This file is part of the SIMPSON General NMR Simulation Package

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
    
    This is the main program that initates the Tcl commands declared
    in various source code files via the 'tclcmd' commands.
    Evaluates the 'main.tcl' Tcl code that loads other statically linked
    Tcl code and evaluates the input file. Next step in program
    flow is main.tcl.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <tcl.h>
#include <pthread.h>

#ifdef WIN32
#include <windows.h>
#include <winbase.h>
#endif

#include "defs.h"
#include "tclutil.h"
#include "OCroutines.h"
#include "rfshapes.h"
#include "cryst.h"

#include "pthread_barrier_mac.h"

int verbose = 0;
int various = 0;
pthread_barrier_t simpson_b_start, simpson_b_end;
glob_info_type glob_info;
int MAXFULLDIM;
int MAXDIMDIAGONALIZE;
double SPARSE_TOL;
double SPARSITY;

extern void simpson_mpi_slave(Tcl_Interp *interp);
extern void simpson_thread_slave( void *thr_id);
//extern int simpson_nfft_test(void);
//extern void simpson_fftw_test(void);

void noprintf(const char* format, ...)
{
   return;
}

int main (int argc,char *argv[])
{
	int process_id = 0;
	int num_processes = 1;
	int par_num_cores, num_threads;
	char *env_np;
	int *thread_id;
	pthread_t *threads;
	pthread_attr_t threadattr;
	Tcl_Interp *interp;
	int i;
	char buf[32];

	OCpar.isinit = 0;
	RFshapes_init();

#ifdef MPI
   MPI_Init(NULL, NULL);
   //--> might need changes to this:
   //int mpi_thread_support;
   //MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &mpi_thread_support);
   //printf("asked for %d, got %d\n",MPI_THREAD_MULTIPLE,mpi_thread_support);
   //if (mpi_thread_support < MPI_THREAD_FUNNELED) {
   //   fprintf(stderr,"Error: MPI lib doesn't support threads\n");
   //   exit(1);
   //}
   //--> end of changes
   //--> NOTE: simpson theoretically needs just MPI_THREAD_FUNNELED level of support
   MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
   MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   int name_len;
   MPI_Get_processor_name(processor_name, &name_len);
   if (process_id == 0) printf("MPI activated\n-------------\n");
   printf("I am process %i of %i running on %s\n", process_id, num_processes, processor_name);
#endif


   /* test input arguments */
   if (argc >= 2) {
	   glob_info.mpi_size = num_processes;
	   glob_info.mpi_rank = process_id;
	   glob_info.cont_thread = 0;
	   glob_info.inputfile = strdup(argv[1]);

	   /* create and initialize master Tcl interpreter */
	   DEBUGPRINT("Main: creating Tcl interpreter\n");
	   interp = Tcl_CreateInterp();
	   if (Tcl_Init(interp) == TCL_ERROR) {
		   fprintf(stderr,"%s is unable to initialize Tcl interpreter. Is init.tcl on your path?\n",PACKAGE);
		   fprintf(stderr,"Error: %s\n",Tcl_GetStringResult(interp));
		   exit(1);
	   }
	   TclSimpsonInterpInit(interp);
	   TclSetSignalHandler(interp,"signalhandler");
	   /* mark variables that will not be transferred to slaves */
	   DEBUGPRINT("Main: executing 'markvars'\n");
	   if (Tcl_Eval(interp, "markvars") != TCL_OK) {
	   	   fprintf(stderr,"Error when executing 'markvars':\n%s\n",Tcl_GetStringResult(interp));
	   	   exit(1);
	   }
	   /* pass input arguments to Tcl interpreter as globals */
	   DEBUGPRINT("Main: passing input arguments to Tcl interpreter\n");
	   sprintf(buf,"%d",argc-1);
	   Tcl_SetVar(interp,"argc",buf,TCL_GLOBAL_ONLY);
	   Tcl_SetVar(interp,"argv0",argv[1],TCL_GLOBAL_ONLY);
	   for (i=1; i<argc; i++) {
		   Tcl_SetVar(interp,"argv",argv[i],TCL_GLOBAL_ONLY|TCL_LIST_ELEMENT|TCL_APPEND_VALUE);
	   }
       if (NULL == Tcl_SetVar(interp,"simpson_version",VERSION,
                              TCL_GLOBAL_ONLY | TCL_LEAVE_ERR_MSG)) {
           fprintf(stderr,"error: %s\n",Tcl_GetStringResult(interp));
           return TCL_ERROR;
       }
	   /* evaluate input file */
	   DEBUGPRINT("Main: evaluating input file '%s'\n\n",argv[1]);
	   if (Tcl_EvalFile(interp, argv[1]) != TCL_OK) {
		   fprintf(stderr,"Error when evaluating input script %s:\n%s\n",argv[1],Tcl_GetStringResult(interp));
		   exit(1);
	   }

	   /* set par(name) variable */
	   if (Tcl_Eval(interp, "set par(name) [file rootname [lindex [file split $argv0] end]]") != TCL_OK) {
	   	   fprintf(stderr,"Error when setting 'par(name)' variable from argument '%s':\n%s\n",argv[1],Tcl_GetStringResult(interp));
	   }
#ifdef MPI
	   /* set par(MPI_rank) and par(MPI_size) */
	   Tcl_Obj *obj = Tcl_NewIntObj(process_id);
	   obj = Tcl_SetVar2Ex(interp,"par","MPI_rank",obj,TCL_GLOBAL_ONLY);
	   if (obj == NULL) {
		   fprintf(stderr,"Error when setting 'par(MPI_rank)':\n%s\n",Tcl_GetStringResult(interp));
		   exit(1);
	   }
	   obj = Tcl_NewIntObj(num_processes);
	   obj = Tcl_SetVar2Ex(interp,"par","MPI_size",obj,TCL_GLOBAL_ONLY);
	   if (obj == NULL) {
		   fprintf(stderr,"Error when setting 'par(MPI_size)':\n%s\n",Tcl_GetStringResult(interp));
		   exit(1);
	   }
#endif
	   /* set up some globals */
	   MAXFULLDIM = TclGetInt(interp,"par","maxfulldim",0,10);
	   MAXDIMDIAGONALIZE = TclGetInt(interp,"par","maxdimdiagonalize",0,4096);
	   SPARSITY = TclGetDouble(interp,"par","sparsity",0,0.8);
	   SPARSE_TOL = TclGetDouble(interp,"par","sparse_tol",0,1.0e-6);
	   /* read par(num_cores) and start threads */
	   par_num_cores = TclGetInt(interp,"par","num_cores",0,0);
#ifdef WIN32
		SYSTEM_INFO sysinfo;
		GetSystemInfo( &sysinfo );
		num_threads = sysinfo.dwNumberOfProcessors;
#else
		num_threads = sysconf (_SC_NPROCESSORS_CONF);
#endif
		if ((env_np = getenv("SIMPSON_NUM_CORES")) != NULL){
			char *endptr;
			num_threads = strtol(env_np, &endptr, 10);
		}
		if (par_num_cores != 0) num_threads = par_num_cores;
		glob_info.num_threads = num_threads;
		threads = (pthread_t *) malloc(sizeof(pthread_t) * num_threads);
		thread_id = (int *) malloc(sizeof(int) * num_threads);
		for (i=0; i<num_threads; i++) thread_id[i] = i;
		pthread_barrier_init(&simpson_b_start, NULL, num_threads+1);
		pthread_barrier_init(&simpson_b_end, NULL, num_threads+1);
		pthread_attr_init(&threadattr);
		pthread_attr_setdetachstate(&threadattr, PTHREAD_CREATE_JOINABLE);
		for (i=0; i<num_threads; i++) {
			pthread_create(&threads[i], &threadattr, (void *)simpson_thread_slave, (void *)(thread_id+i));
		}

	   if (process_id == 0) {
		   /* MASTER does Tcl main evaluations */

		   /* execute the 'main' section from input file */
		   DEBUGPRINT("\n\nMain: evaluating 'main' section\n\n");
		   if (Tcl_Eval(interp, "main") != TCL_OK) {
			   fprintf(stderr,"Error when evaluating 'main' section of input script %s:\n%s\n",argv[1],Tcl_GetStringResult(interp));
		   }

		    //simpson_fftw_test();
			//simpson_nfft_test();

		   // send message to MPI slaves to terminate (tag = 1)
#ifdef MPI
		   int terminate = 1;
		   for (i=1; i< num_processes; i++) {
			   MPI_Send(&terminate, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
		   }
#endif
		   Tcl_DeleteInterp(interp);
		   DEBUGPRINT("\n\nSuccessful end of %s program!\n",PACKAGE);

	   } else {
#ifdef MPI

		   /* slaves jump to simpson() and wait for master */
		   simpson_mpi_slave(interp);

#endif
	   }

	   /* clean up and terminate threads */
	   glob_info.cont_thread = 0;
	   pthread_barrier_wait(&simpson_b_start);
	   for (i=0; i<num_threads; i++) {
		   pthread_join(threads[i], NULL);
	   }
	   DEBUGPRINT("threads successfully terminated.\n");
	   pthread_attr_destroy(&threadattr);
	   pthread_barrier_destroy(&simpson_b_start);
	   pthread_barrier_destroy(&simpson_b_end);
	   free(threads);
	   free(thread_id);



   } else {
	   /* not enough arguments */
	   if (process_id == 0) {
		   fprintf(stderr,"%s version %s, Copyright (C)\n",PACKAGE,VERSION);
		   fprintf(stderr,"1999-2000 Mads Bak and Jimmy T. Rasmussen\n"
	                  "2001 Mads Bak and Thomas Vosegaard\n"
	                  "2002-2007 Thomas Vosegaard\n"
	                  "2008-2009 Zdenek Tosner, Thomas Vosegaard, and Niels Chr. Nielsen\n"
			          "2009 plus Rasmus Andersen\n"
	                  "2010-2014 Zdenek Tosner, Niels Chr. Nielsen, Thomas Vosegaard\n");
		   fprintf(stderr,"%s comes with ABSOLUTELY NO WARRANTY, for details\n",PACKAGE);
		   fprintf(stderr,"read the COPYING file included in this distribution.\n"
	                  "This is free software, and you are welcome to redistribute\n"
	                  "it under certain conditions according to the GNU General Public License.\n"
	                  "\nPlease specify an inputfile, optionally with other arguments.\n");
	   }
   }


   /* clean up and terminate */
#ifdef MPI
   MPI_Finalize();
#endif

   return 0;
}

