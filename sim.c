/*
    Simulation setup and calculation
    Copyright (C) 1999 Mads Bak, Jimmy T. Rasmussen
                  2010 Zdenek Tosner

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
    
    Makes available a setup routine, and a routine for calculation of
    a single crystallite.
    Uses the 'readsys' routine to create the spin-system. Uses the
    'pulse_propagate' routine that performs evaluation of the
    pulses. Uses the functions in 'fidcalc' to perform the 
    evolution of the fid in case of a smart gamma averaged method.
    
    Called by simpson.c that makes the crystallite powder averaging.
*/

#include <stdlib.h>
#include <math.h>
#include <tcl.h>
#include <string.h>
#include "matrix.h"
#include "tclutil.h"
#include "sim.h"
#include "ham.h"
#include "fidcalc.h"
#include "relax.h"
#include "defs.h"
#include "pulse.h"
#include "spinsys.h"
#include "cm.h"
#include "blockdiag.h"
#include "OCroutines.h"
#include "rfprof.h"

#define SPINRATE_SMALL 0.001
#define NOTSET -123.456

extern glob_info_type glob_info;

void readsys(Tcl_Interp* interp,Sim_info * s);

/* finds the number of the observable nucleus given the detect operator*/
int obs_nuc(SpinSys* ss, char* det,int* was_all)
{
  char buf[256];
  int i,nuc=0,l,j,nuc2,iso2,iso=0;

  *was_all=0;
  l=strlen(det);
  for (i=0;i<l;i++) {
    if (det[i] == 'I') {
      j=0; i++;
      if (det[i] == 'n') {
        *was_all=1;
      } else {
        while (det[i] >= '0' && det[i] <= '9') { 
          buf[j++]=det[i++];
        }
        buf[j]=0;
        nuc2=strtol(buf,(char**)NULL, 10);
        if (nuc2 == 0) {
          fprintf(stderr,"error: unable to convert '%s' in '%s' to an integer value\n",buf,det);
          exit(1);
        }
        iso2 = ss_isotope(ss,nuc2)->number;
        if (nuc == 0) {
          nuc=nuc2;
          iso=iso2;
        } else if (iso != iso2)
          return 0;
      }      
    }
  }
  if (*was_all && nuc == 0) return 1;
  return nuc;
}

Sim_info * sim_initialize(Tcl_Interp* interp)
{
  char detectop[128],startop[256],pulseq[256],buf[256];
  double f;
  Sim_info * s;
  Tcl_Obj *obj, **objv;
  int i, objc;

  s = (Sim_info*)malloc(sizeof(Sim_info));
  if (!s) {
	  fprintf(stderr,"error: can not allocate Sim_info structure\n");
	  exit(1);
  }

  s->ss=(SpinSys*)malloc(sizeof(SpinSys));
  if (!s->ss) {
    fprintf(stderr,"error: allocation failure in sim_initialize()\n");
    exit(1);
  } 

  f=TclGetDouble(interp,"par","proton_frequency",0,400e6);  
  if (f < 0.0) {
     fprintf(stderr,"error: proton_frequency must be positive\n");
     exit(1);
  }
  if (f != 0 && f <= 10000) {
     fprintf(stderr,"error: proton_frequency must be given in Hz\n");
     exit(1);
  }
  /* Internal in the program we use omega= -gamma B0, that is the
     proton resonance frequency is negative because gamma for at
     proton is positive. */
  s->specfreq = -f;
    
  s->wr=TclGetDouble(interp,"par","spin_rate",0,0);
  if (s->wr < 0.0) {
     fprintf(stderr,"error: spin_rate cannot be negative\n");
     exit(1);
  }
  s->sw=TclGetDouble(interp,"par","sw",0,s->wr);
  s->sw1=TclGetDouble(interp,"par","sw1",0,0);
  s->ni=TclGetInt(interp,"par","ni",0,0);

  // read and decompose par(method)
  //  DEFAULTS:
  s->imethod = M_DIRECT_TIME;
  s->sparse = 0;
  s->block_diag = 0;
  s->propmethod = 0; // via diagonalization
  s->interpolation = INTERPOL_NOT_USED;
  s->domain = 0; // time domain simulation
  s->labframe = 0; // rotating frame is default
  obj = Tcl_GetVar2Ex(interp,"par","method",0);
  if (obj != NULL) {
	  if (Tcl_ListObjGetElements(interp,obj,&objc,&objv) != TCL_OK) {
		  fprintf(stderr,"Error: problem reading method definition in par:\n%s\n",Tcl_GetStringResult(interp));
		  exit(1);
	  }
	  for (i=0; i<objc; i++) {
		  strcpy(buf,Tcl_GetString(objv[i]));
		  if (!strcmp(buf,"gcompute")){
		    s->imethod=M_GCOMPUTE_TIME;
		  } else if ( !strcmp(buf,"direct")){
		    s->imethod=M_DIRECT_TIME;
		  } else if ( !strcmp(buf,"idirect")){
		    s->imethod=M_DIRECT_TIME;
		  } else if ( !strcmp(buf,"igcompute")){
		    s->imethod=M_GCOMPUTE_TIME;
		  } else if ( !strcmp(buf,"spinach")){
		    s->imethod=M_SPINACH;
		  } else if (!strncmp(buf,"diag",4)) {
			  s->propmethod = 0; // this uses dsyevr that might not be thread safe
		  } else if (!strncmp(buf,"pade",4)) {
			  s->propmethod = 1;
		  } else if (!strncmp(buf,"cheby1",6)) {
			  s->propmethod = 2; // with Scaling & Squaring
		  } else if (!strncmp(buf,"cheby2",6)) {
			  s->propmethod = 6; // with Shifting & Scaling
		  } else if (!strncmp(buf,"taylor",6)) {
			  s->propmethod = 3;
		  } else if (!strncmp(buf,"lanczos",7)) {
			  s->propmethod = 4;
		  } else if (!strncmp(buf,"dsyev",5)) {
			  s->propmethod = 5; // diag with dsyev found always thread safe but slower
		  } else if (!strncmp(buf,"sparse",6)) {
			  s->sparse = 1;
		  } else if (!strncmp(buf,"block_diag",10)) {
			  s->block_diag = 1;
		  } else if (!strncmp(buf,"FWTinterpolation",9)) {
			  s->interpolation = INTERPOL_FWT_ALL;
		  } else if (!strncmp(buf,"FWT2interpolation",10)) {
			  s->interpolation = INTERPOL_FWT_LAM;
		  } else if (!strncmp(buf,"ASGinterpolation",9)) {
			  s->interpolation = INTERPOL_ASG;
		  } else if (!strncmp(buf,"FWTASGinterpolation",12)) {
			  s->interpolation = INTERPOL_FWTASG_ALL;
		  } else if (!strncmp(buf,"FWT2ASGinterpolation",12)) {
			  s->interpolation = INTERPOL_FWTASG_LAM;
		  } else if (!strncmp(buf,"time",4)) {
			  s->domain = 0;
		  } else if (!strncmp(buf,"frequency",4)) {
			  s->domain = 1; // calculations will be in frequency domain
			  s->imethod++;
		  } else if (!strncmp(buf,"labframe",8)) {
			  s->labframe = 1;
		  } else {
			  fprintf(stderr,"Error: method option '%s' not known.\n",buf);
			  fprintf(stderr,"Must be one of : direct/idirect/gcompute/igcompute,\n");
			  fprintf(stderr,"               : diag/dsyev/pade/cheby1/cheby2/taylor, \n");
			  fprintf(stderr,"               : sparse, block_diag, \n");
			  fprintf(stderr,"               : time/frequency, \n");
			  fprintf(stderr,"               : FWTinterpolation/ASGinterpolation/FWTASGinterpolation\n");
			  exit(1);

		  }
	  }  // for loop over objc
  }
  // par(method) parsed.
#ifdef NO_NFFT
  if ( (s->interpolation != INTERPOL_NOT_USED) && (s->interpolation != INTERPOL_ASG) ) {
	  fprintf(stderr,"Error: simpson compiled without NFFT library - FWT interpolation not possible!\n");
	  exit(1);
  }
#endif

  s->points_per_cycle = TclGetInt(interp,"par","points_per_cycle",0,8);

  TclGetString(interp,buf,"par","name",1,"");
  strcpy(s->parname,buf);
  
  s->dor=TclGetInt(interp,"par","dor",0,0);
  if (s->dor == 1) {
    s->brl1=TclGetDouble(interp,"par","inner_rotor_angle",1,0);
    s->brl2=TclGetDouble(interp,"par","outer_rotor_angle",1,0);
    s->wr1=TclGetDouble(interp,"par","inner_spin_rate",1,0);
    s->wr2=TclGetDouble(interp,"par","outer_spin_rate",1,0);
    s->wr=s->wr1;
  } else {
    s->brl=TclGetDouble(interp,"par","rotor_angle",0,NOTSET);
    if (s->brl == NOTSET) {    
      s->brl=(s->wr == 0.0 ? 0.0 : MAGIC_ANGLE);
      if (verbose & VERBOSE_SIMINFO) {
        printf("The rotor_angle is set to %g\n",s->brl);      
      }
    } 
    s->gamma_zero=TclGetDouble(interp,"par","gamma_zero",0,0);
  }
  s->ngamma=TclGetInt(interp,"par","gamma_angles",0,1);

  s->acq_adjoint = TclGetInt(interp,"par","acq_adjoint",0,0);

  TclGetString(interp,buf,"par","dipole_check",0,"true");
  if (!strcmp(buf,"true")) {
    s->dipole_check=1;
  } else if (!strcmp(buf,"false")) {
    s->dipole_check=0;
  } else {
    fprintf(stderr,"error: 'dipole_check' must be 'false' or 'true'\n");
    exit(1);
  }

  s->np=TclGetInt(interp,"par","np",0,1);
  s->ntot=s->np*(s->ni > 1 ? s->ni : 1);

  s->realspec=TclGetInt(interp,"par","real_spec",0,0);

  //TclGetString(interp,s->crystfile,"par","crystal_file",0,"alpha0beta90");
  obj = Tcl_GetVar2Ex(interp,"par","crystal_file",0);
  if (obj == NULL) {
	  strcpy(s->crystfile,"alpha0beta90");
	  s->crystfile_from = 1;
	  s->crystfile_to = 1;
	  strcpy(s->targetcrystfile,"");
  } else {
	  if (Tcl_ListObjGetElements(interp,obj,&objc,&objv) != TCL_OK) {
		  fprintf(stderr,"Error: problem reading crystal_file definition in par:\n%s\n",Tcl_GetStringResult(interp));
		  exit(1);
	  }
	  switch (objc) {
	  case 1:
		  strcpy(s->crystfile,Tcl_GetString(objv[0]));
		  s->crystfile_from = 1;
		  s->crystfile_to = -1;
		  strcpy(s->targetcrystfile,"");
		  break;
	  case 2:
		  strcpy(s->crystfile,Tcl_GetString(objv[0]));
		  s->crystfile_from = 1;
		  s->crystfile_to = -1;
		  strcpy(s->targetcrystfile,Tcl_GetString(objv[1]));
		  break;
	  case 3:
		  strcpy(s->crystfile,Tcl_GetString(objv[0]));
		  if (Tcl_GetIntFromObj(interp,objv[1],&(s->crystfile_from)) != TCL_OK) {
			  fprintf(stderr,"Error: problem reading crystal_file definition in par (from):\n%s\n",Tcl_GetStringResult(interp));
			  exit(1);
		  }
		  if (Tcl_GetIntFromObj(interp,objv[2],&(s->crystfile_to)) != TCL_OK) {
			  fprintf(stderr,"Error: problem reading crystal_file definition in par (to):\n%s\n",Tcl_GetStringResult(interp));
			  exit(1);
		  }
		  if (s->crystfile_to < s->crystfile_from) {
			  fprintf(stderr,"Error: problem reading crystal_file definition in par:\n 'to' is smaller than 'from'\n");
			  exit(1);
		  }
		  strcpy(s->targetcrystfile,"");
		  break;
	  default:
		  fprintf(stderr,"Error: problem reading crystal_file definition in par:\nlist elements count mismatch\n");
		  fprintf(stderr,"Allowed formats: 'rep10', 'rep100 1 25','LEB217 ROSELEB11617'\n");
		  exit(1);
	  }
  }
  //printf("\nsim_init test cryst_file: '%s' (%d, %d), '%s'\n",s->crystfile,s->crystfile_from,s->crystfile_to,s->targetcrystfile);
  s->crdata = read_crystfile(s->crystfile, s->crystfile_from, s->crystfile_to);
  if (s->interpolation == INTERPOL_FWT_ALL || s->interpolation == INTERPOL_FWT_LAM || s->interpolation == INTERPOL_FWTASG_ALL || s->interpolation == INTERPOL_FWTASG_LAM) {
	  if (strlen(s->targetcrystfile) == 0) {
		  fprintf(stderr,"Error: no target crystallite set defined for FWT interpolation.\n");
		  exit(1);
	  }
	  s->targetcrdata = read_crystfile(s->targetcrystfile, s->crystfile_from, s->crystfile_to);
	  // WARNING: this holds ONLY for LEBh/ROSELEBh sets!!!
	  if (strncmp(s->targetcrystfile,"LEBh",4) && strncmp(s->targetcrystfile,"ROSELEBh",8)) {
		  fprintf(stderr,"Error: at the moment, only LEBhXXX/ROSELEBhXXX source orientation sets are allowed in FWT interpolation!\n");
		  exit(1);
	  }
	  s->Jinterpol[0] = ( (int)sqrt(6*(LEN(s->crdata)-1)) - 1) / 2;
	  s->Jinterpol[1] = 2 << (int)(floor(log2((double)s->Jinterpol[0])));
  } else {
	  s->targetcrdata = NULL;
  }
  if (s->interpolation == INTERPOL_FWT_LAM || s->interpolation == INTERPOL_FWTASG_LAM) {
	  // load, or create and load, map of nearest crystallites target->source
	  s->crmap = read_cryst_map(s->crystfile, s->crdata, s->targetcrystfile, s->targetcrdata);
  } else {
	  s->crmap = NULL;
  }
  if (s->interpolation == INTERPOL_ASG || s->interpolation == INTERPOL_FWTASG_ALL || s->interpolation == INTERPOL_FWTASG_LAM) {
	  // test for triangle data
	  // at the moment they are read in mpi_ASG_interpol ...
  }

  TclGetString(interp,detectop,"par","detect_operator",0,"Inp");
  TclGetString(interp,s->rfproffile,"par","rfprof_file",0,"none");
  TclGetString(interp,startop,"par","start_operator",0,"Inz");
  TclGetString(interp,pulseq,"par","pulse_sequence",0,"pulseq");
  strcpy(s->pulsename, pulseq);

  obj = Tcl_GetVar2Ex(interp,"par","averaging_file",0);
  if (obj == NULL) {
	  s->averaging = 0;
  } else {
	  s->averaging = 1;
  }

  readsys(interp,s);

  s->matdim=s->ss->matdim;
  s->fdetect = ss_readoper(s,detectop);
  if (s->fdetect == NULL) {
     fprintf(stderr,"error: unable to interpret detect operator '%s'\n",detectop);
     exit(1);
  }
  if (!strncmp(startop,"equil",5)) {
     s->fstart = ss_eq_sigma(s);
  } else {
     s->fstart = ss_readoper(s,startop);
  }
  if (s->fstart == NULL) {
     fprintf(stderr,"error: unable to interpret start operator '%s'\n",startop);
     exit(1);
  }
  /* conjugate_fid: This parameter is used overrule the automatic detection of the magnetogyric
    ratio for the observable nucleus. It can be set to 'auto', 'true' or 'false.
  */
  TclGetString(interp,buf,"par","conjugate_fid",0,"auto");
  if (!strcmp(buf,"auto")) {
    int was_all=0;
    s->obs_nuc = obs_nuc(s->ss,detectop,&was_all);
    if ( (s->obs_nuc < 1) || (was_all && s->ss->nchan > 1) ) {
      fprintf(stderr,"error: in detect-operator '%s'. More than one type of nucleus\n"
                     "was detected. Set 'conjugate_fid' to 'true' or 'false' to turn off\n"
                     "this sanity check. 'conjugate_fid' is 'auto' per default which conjugates\n"
                     "the fid if the gyromagnetic ratio is larger than zero for the observable\n" 
                     "nucleus. That corrects for the standard axis convention of plotting spectra\n"
                     "that ignores the sign of gamma.\n"
                     ,detectop);
      exit(1);
    }
    s->conjugate_fid=  (ss_gamma(s->ss,s->obs_nuc) > 0 ? 1 : 0);
  } else if (!strcmp(buf,"true")) {
    s->conjugate_fid = 1;
  } else if (!strcmp(buf,"false")) {
    s->conjugate_fid = 0;
  } else {
    fprintf(stderr,"error: 'conjugate_fid' must be 'true', 'false' or 'auto'");
    exit(1);
  }

  if (verbose & VERBOSE_OPER) {
    printf("Acquisition data will%s be complex conjugated.\n",
       (s->conjugate_fid ? "" : " not"));
  }

  // Excitation-Detection symmetry check
  //  1 ... assume there is ED symmetry
  //  0 ... force internal check on ED symmetry
  // -1 ... no ED symmetry, no internal check
  s->EDsymmetry = TclGetInt(interp,"par","ED_symmetry",0,-1);
  if (s->EDsymmetry == 0) s->EDsymmetry = is_rhosymmetry(s->fstart,s->fdetect);
  if (abs(s->EDsymmetry)>1) {
	  fprintf(stderr,"Error: ED_symmetry parameter has wrong value (can be -1, 0, 1)\n");
	  exit(1);
  }
  //printf("ED symmetry is %d\n",s->EDsymmetry);

  if (verbose & VERBOSE_OPER) {
     cm_print(s->fstart,"Start operator");
     cm_print(s->fdetect,"Detect operator");
     printf("Excitation-Detection symmetry: ");
     if (s->EDsymmetry < 0)
    	 printf("not assumed\n");
   	 else
   		 printf("is assumed\n");
  }

  if ( (s->imethod == M_DIRECT_TIME) || (s->imethod == M_DIRECT_FREQ) || (s->imethod == M_SPINACH)) {
    if (s->wr == 0.0) {
      if (s->ngamma != 1) {
        fprintf(stderr,"error: 'gamma_angles' must be one when calculating static spectra\n");
        exit(1);
      }
    }
  } else  {
    if (s->dor==1) {
      fprintf(stderr,"error: '(i)gcompute' method not compatible with calculating DOR spectra\n");
      exit(1);
    }
    if (s->wr < SPINRATE_SMALL) {
      fprintf(stderr,"error: must use the 'direct' or 'idirect' method when calculating static spectra (spin-rate=0).\n");
      exit(1);
    }
  }
  s->wr *= M_PI*2.0;
  if (s->dor) {
     s->wr1 *= M_PI*2.0;
     s->wr2 *= M_PI*2.0;
  }

  /* ZT: relaxation setup */
  char rxbuf[32];
  TclGetString(interp,rxbuf,"par","relax",0,"off");
  /* printf("RX: relax is %s\n",rxbuf); */
  if (!strncmp(rxbuf,"on",2)) {
     s->relax = 1;
     read_relax(interp,s);
  } else {
     s->relax = 0;
  }

  /* other parameters */
  //s->dw = 1.0e6/s->sw; // dwelltime in us
  if (s->wr > TINY) {
	  s->taur = 2.0e6*M_PI/s->wr;  // rotor period in us
  } else {
	  s->taur = -1;
  }

  /* spinach interfers with conjugate_fid strangely */
  //if (s->imethod == M_SPINACH) s->conjugate_fid = !(s->conjugate_fid);

  s->rfdata = read_rfproffile(s->rfproffile,s->ss->nchan);

  s->tridata = NULL;
  s->ASG_freq = NULL;
  s->ASG_ampl = NULL;
  s->FWTASG_nnz = NULL;
  s->FWT_frs = NULL;
  s->FWT_lam = NULL;
  s->FWT_icol = NULL;
  s->FWT_irow = NULL;


  /* hack for FFTW thread safety - create plan only once and remember it */
  if (s->domain == 1) { /* calculations in frequency domain */
	  if (s->imethod == M_DIRECT_FREQ) {
		  i = s->points_per_cycle;
	  } else {
		  i = s->ngamma;
	  }
	  fftw_complex *fftin1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * i);
	  fftw_complex *fftout1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * i);
	  s->fftw_plans = (fftw_plan*)malloc(glob_info.num_threads*sizeof(fftw_plan));
	  int kk;
	  //printf("\nCreating %d FFTW plans with size of %d \n",glob_info.num_threads,i);
	  for (kk=0; kk<glob_info.num_threads; kk++) s->fftw_plans[kk] = fftw_plan_dft_1d(i, fftin1, fftout1, FFTW_FORWARD, FFTW_MEASURE);
	  fftw_free(fftin1);
	  fftw_free(fftout1);
  } else {
	  s->fftw_plans = NULL;
  }

  s->do_avg = TclGetInt(interp,"par","do_avg",0,0);

  return s;
}

void sim_destroy(Sim_info* s, int this_is_copy)
{
  int i, ii;

  /* ZT: relaxation setup */
  if (s->relax) destroy_Relax(s->ss->nspins);

  DEBUGPRINT("SIM_DESTROY 1\n");
  
  free(s->ss);
  DEBUGPRINT("SIM_DESTROY 2\n");
  for (i=0; i<s->nCS; i++) {
	  Shift *ptr = s->CS[i];
	  if (ptr->T) free_double_matrix(ptr->T);
	  free_complx_vector(ptr->Rmol);
	  for (ii=0; ii<5; ii++) {
		  if (ptr->T2q[ii] != NULL) free_double_matrix(ptr->T2q[ii]);
	  }
	  free((char*)(ptr));
  }
  free(s->CS);
  DEBUGPRINT("SIM_DESTROY 3\n");
  for (i=0; i<s->nJ; i++) {
	  Jcoupling *ptr = s->J[i];
	  if (ptr->blk_T) free_blk_mat_double(ptr->blk_T);
	  if (ptr->blk_Tiso) free_blk_mat_double(ptr->blk_Tiso);
	  free_complx_vector(ptr->Rmol);
	  for (ii=0; ii<5; ii++) {
		  if (ptr->T2q[ii] != NULL) free_double_matrix(ptr->T2q[ii]);
	  }
 	  free((char*)(ptr));
 }
  free(s->J);
  DEBUGPRINT("SIM_DESTROY 4\n");
  for (i=0; i<s->nDD; i++) {
	  Dipole *ptr = s->DD[i];
	  if (ptr->blk_T) free_blk_mat_double(ptr->blk_T);
	  free_complx_vector(ptr->Rmol);
	  for (ii=0; ii<5; ii++) {
		  if (ptr->T2q[ii] != NULL) free_double_matrix(ptr->T2q[ii]);
	  }
	  free((char*)(ptr));
  }
  free(s->DD);
  DEBUGPRINT("SIM_DESTROY 5\n");
  for (i=0; i<s->nQ; i++) {
	  Quadrupole *ptr = s->Q[i];
	  if (ptr->T) free_double_matrix(ptr->T);
	  if (ptr->Ta) free_double_matrix(ptr->Ta);
	  if (ptr->Tb) free_double_matrix(ptr->Tb);
	  if (ptr->T3a) free_double_matrix(ptr->T3a);
	  if (ptr->T3b) free_double_matrix(ptr->T3b);
	  if (ptr->T3c) free_double_matrix(ptr->T3c);
	  free_complx_vector(ptr->Rmol);
	  for (ii=0; ii<5; ii++) {
		  if (ptr->T2q[ii] != NULL) free_double_matrix(ptr->T2q[ii]);
	  }
 	  free((char*)(ptr));
 }
  free(s->Q);
  DEBUGPRINT("SIM_DESTROY 6\n");
  for (i=0; i<s->nMIX; i++) {
	  Mixing *ptr = s->MIX[i];
	  if (ptr->T) free_double_matrix(ptr->T);
	  if (ptr->Ta) free_blk_mat_double(ptr->Ta);
	  if (ptr->Tb) free_blk_mat_double(ptr->Tb);
	  free((char*)(ptr));
  }
  free(s->MIX);
  DEBUGPRINT("SIM_DESTROY 7\n");
  if (s->fstart) free_complx_matrix(s->fstart);
  DEBUGPRINT("SIM_DESTROY 8\n");
  if (s->fdetect) free_complx_matrix(s->fdetect);
  DEBUGPRINT("SIM_DESTROY 9\n");

  /* new Hamiltonian assembling */
  if (s->Hiso != NULL) free_blk_mat_double(s->Hiso);
  DEBUGPRINT("SIM_DESTROY 10\n");
  for (i=0; i<5; i++) {
	  if (s->HQ[i] != NULL) free_blk_mat_double(s->HQ[i]);
  }
  DEBUGPRINT("SIM_DESTROY 11\n");

  if (s->perm_table != NULL) {
	  for (i=0; i< (s->Nbasis*s->Nbasis); i++) {
		  if (s->perm_table[i]) free_int_vector(s->perm_table[i]);
	  }
	  free(s->perm_table);
	  s->perm_table = NULL;
  }
  DEBUGPRINT("SIM_DESTROY 12\n");
  if (s->dims_table != NULL) {
	  for (i=0; i<s->Nbasis; i++) {
		  if (s->dims_table[i]) free_int_vector(s->dims_table[i]);
	  }
	  free(s->dims_table);
	  s->dims_table = NULL;
  }
  DEBUGPRINT("SIM_DESTROY 13\n");

  if (!this_is_copy) {
	  rfprof_free(s->rfdata);
	  cryst_free(s->crdata);
	  if (s->targetcrdata != NULL) {
		  cryst_free(s->targetcrdata);
		  s->targetcrdata = NULL;
	  }
  }

  if (s->tridata != NULL) {
	  triangle_free(s->tridata);
	  s->tridata = NULL;
  }
  if (s->ASG_ampl != NULL) {
	  free(s->ASG_ampl);
	  s->ASG_ampl = NULL;
  }
  if (s->ASG_freq != NULL) {
	  free(s->ASG_freq);
	  s->ASG_freq = NULL;
  }
  if (s->FWT_lam != NULL) {
	  free(s->FWT_lam);
	  s->FWT_lam = NULL;
  }
  if (s->FWT_frs != NULL) {
	  free(s->FWT_frs);
	  s->FWT_frs = NULL;
  }
  if (s->FWT_irow != NULL) {
	  free(s->FWT_irow);
	  s->FWT_irow = NULL;
  }
  if (s->FWT_icol != NULL) {
	  free(s->FWT_icol);
	  s->FWT_icol = NULL;
  }
  if (s->crmap != NULL) free(s->crmap);

  if (s->fftw_plans != NULL) {
	  for (i=0; i<glob_info.num_threads; i++) fftw_destroy_plan(s->fftw_plans[i]);
	  free(s->fftw_plans);
	  s->fftw_plans = NULL;
  }

  if (s->FWTASG_nnz != NULL) free_int_vector(s->FWTASG_nnz);
  if (s->FWT_lam != NULL) {
	  free(s->FWT_lam);
	  s->FWT_lam = NULL;
  }
  if (s->FWT_frs != NULL) {
	  free(s->FWT_frs);
	  s->FWT_frs = NULL;
  }
  if (s->FWT_irow != NULL) {
	  free(s->FWT_irow);
	  s->FWT_irow = NULL;
  }
  if (s->FWT_icol != NULL) {
	  free(s->FWT_icol);
	  s->FWT_icol = NULL;
  }
  if (s->ASG_ampl != NULL) {
	  free(s->ASG_ampl);
	  s->ASG_ampl = NULL;
  }
  if (s->ASG_freq != NULL) {
	  free(s->ASG_freq);
	  s->ASG_freq = NULL;
  }

}

/****
 * decided NOT to use this anymore - NOT UPDATED!!!!
 *
Sim_info * sim_duplicate(Sim_info *sim)
{
	int i, j;

	Sim_info * sim2 = (Sim_info*)malloc(sizeof(Sim_info));
	if (sim2 == NULL) {
		fprintf(stderr,"Error: sim_duplicate - can not allocate Sim_info structure\n");
		exit(1);
	}
	sim2->nCS = sim->nCS;
	sim2->nDD = sim->nDD;
	sim2->nJ = sim->nJ;
	sim2->nQ = sim->nQ;
	sim2->nMIX = sim->nMIX;
	sim2->specfreq = sim->specfreq;
	sim2->wr = sim->wr;
	sim2->sw = sim->sw;
	sim2->sw1 = sim->sw1;
	sim2->brl = sim->brl;
	sim2->gamma_zero = sim->gamma_zero;
	sim2->np = sim->np;
	sim2->ni = sim->ni;
	sim2->ntot = sim->ntot;
	sim2->ngamma = sim->ngamma;
	sim2->matdim = sim->matdim;
	sim2->obs_nuc = sim->obs_nuc;
	sim2->imethod = sim->imethod;
	sim2->crystfile_from = sim->crystfile_from;
	sim2->crystfile_to = sim->crystfile_to;
	sim2->acq_adjoint = sim->acq_adjoint;
	sim2->dipole_check = sim->dipole_check;
	sim2->realspec = sim->realspec;
	sim2->sparse = sim->sparse;
	sim2->conjugate_fid = sim->conjugate_fid;
	sim2->EDsymmetry = sim->EDsymmetry;
	sim2->Hint_isdiag = sim->Hint_isdiag;
	sim2->averaging = sim->averaging;
	sim2->block_diag = sim->block_diag;
	sim2->dor = sim->dor;
	sim2->brl1 = sim->brl1;
	sim2->brl2 = sim->brl2;
	sim2->wr1 = sim->wr1;
	sim2->wr2 = sim->wr2;
	sim2->relax = sim->relax;
	sim2->propmethod = sim->propmethod;
	sim2->interpolation = sim->interpolation;
	//sim2->dw = sim->dw;
	sim2->dw1 = sim->dw1;
	sim2->taur = sim->taur;
	sim2->basis = sim->basis;
	sim2->Nbasis = sim->Nbasis;
	sim2->Nmz = sim->Nmz;
	sim2->Hassembly = sim->Hassembly;
	if (sim->fstart != NULL) {
		sim2->fstart = cm_dup(sim->fstart);
	} else {
		sim2->fstart = NULL;
	}
	if (sim->fdetect != NULL) {
		sim2->fdetect = cm_dup(sim->fdetect);
	} else {
		sim2->fdetect = NULL;
	}
	if (sim->Hiso != NULL) {
		sim2->Hiso = blk_dm_dup(sim->Hiso);
	} else {
		sim2->Hiso = NULL;
	}
	for (i=0; i<5; i++) {
		if (sim->HQ[i] != NULL) {
			sim2->HQ[i] = blk_dm_dup(sim->HQ[i]);
		} else {
			sim2->HQ[i] = NULL;
		}
	}
	strcpy(sim2->crystfile,sim->crystfile);
	strcpy(sim2->targetcrystfile,sim->targetcrystfile);
	strcpy(sim2->rfproffile,sim->rfproffile);
	strcpy(sim2->pulsename,sim->pulsename);
	strcpy(sim2->parname,sim->parname);
	if (sim->perm_table != NULL) {
		sim2->perm_table = (int**)malloc(sim->Nbasis*sim->Nbasis*sizeof(int*));
		for (i=0; i< (sim->Nbasis*sim->Nbasis); i++) {
			if (sim->perm_table[i] != NULL) {
				sim2->perm_table[i] = iv_dup(sim->perm_table[i]);
			} else {
				sim2->perm_table[i] = NULL;
			}
		}
	} else {
		sim2->perm_table = NULL;
	}
	if (sim->dims_table != NULL) {
		sim2->dims_table = (int**)malloc(sim->Nbasis*sizeof(int*));
		for (i=0; i<sim->Nbasis; i++) {
			if (sim->dims_table[i] != NULL) {
				sim2->dims_table[i] = iv_dup(sim->dims_table[i]);
			} else {
				sim2->dims_table[i] = NULL;
			}
		}
	} else {
		sim2->dims_table = NULL;
	}

	sim2->CS = (Shift**)malloc(sim->nCS*sizeof(Shift*));
	for (i=0; i<sim2->nCS; i++) {
		Shift *ptr = sim2->CS[i] = (Shift*)malloc(sizeof(Shift));
		ptr->nuc = sim->CS[i]->nuc;
		ptr->iso = sim->CS[i]->iso;
		ptr->delta = sim->CS[i]->delta;
		ptr->eta = sim->CS[i]->eta;
		ptr->pas[0] = sim->CS[i]->pas[0];
		ptr->pas[1] = sim->CS[i]->pas[1];
		ptr->pas[2] = sim->CS[i]->pas[2];
		if (sim->CS[i]->T != NULL) {
			ptr->T = dm_dup(sim->CS[i]->T);
		} else {
			ptr->T = NULL;
		}
		if (sim->CS[i]->Rmol != NULL) {
			ptr->Rmol = cv_dup(sim->CS[i]->Rmol);
		} else {
			ptr->Rmol = NULL;
		}
	}

	sim2->J = (Jcoupling**)malloc(sim->nJ*sizeof(Jcoupling*));
	for (i=0; i<sim2->nJ; i++) {
		Jcoupling *ptr = sim2->J[i] = (Jcoupling*)malloc(sizeof(Jcoupling));
		ptr->nuc[0] = sim->J[i]->nuc[0];
		ptr->nuc[1] = sim->J[i]->nuc[1];
		ptr->iso = sim->J[i]->iso;
		ptr->delta = sim->J[i]->delta;
		ptr->eta = sim->J[i]->eta;
		ptr->pas[0] = sim->J[i]->pas[0];
		ptr->pas[1] = sim->J[i]->pas[0];
		ptr->pas[2] = sim->J[i]->pas[0];
		if (sim->J[i]->blk_T != NULL){
			ptr->blk_T = blk_dm_dup(sim->J[i]->blk_T);
		} else {
			ptr->blk_T = NULL;
		}
		if (sim->J[i]->blk_Tiso != NULL) {
			ptr->blk_Tiso = blk_dm_dup(sim->J[i]->blk_Tiso);
		} else {
			ptr->blk_Tiso = NULL;
		}
		if (sim->J[i]->Rmol != NULL) {
			ptr->Rmol = cv_dup(sim->J[i]->Rmol);
		} else {
			ptr->Rmol = NULL;
		}
	}

	sim2->DD = (Dipole**)malloc(sim->nDD*sizeof(Dipole*));
	for (i=0; i<sim2->nDD; i++) {
		Dipole *ptr = sim2->DD[i] = (Dipole*)malloc(sizeof(Dipole));
		ptr->nuc[0] = sim->DD[i]->nuc[0];
		ptr->nuc[1] = sim->DD[i]->nuc[1];
		ptr->delta = sim->DD[i]->delta;
		ptr->eta = sim->DD[i]->eta;
		ptr->pas[0] = sim->DD[i]->pas[0];
		ptr->pas[1] = sim->DD[i]->pas[1];
		ptr->pas[2] = sim->DD[i]->pas[2];
		if (sim->DD[i]->blk_T != NULL) {
			ptr->blk_T = blk_dm_dup(sim->DD[i]->blk_T);
		} else {
			ptr->blk_T = NULL;
		}
		ptr->Rmol = cv_dup(sim->DD[i]->Rmol);
	}

	sim2->Q = (Quadrupole**)malloc(sim->nQ*sizeof(Quadrupole*));
	for (i=0; i<sim2->nQ; i++) {
		Quadrupole *ptr = sim2->Q[i] = (Quadrupole*)malloc(sizeof(Quadrupole));
		ptr->nuc = sim->Q[i]->nuc;
		ptr->order = sim->Q[i]->order;
		ptr->w0 = sim->Q[i]->w0;
		ptr->wq = sim->Q[i]->wq;
		ptr->eta = sim->Q[i]->eta;
		ptr->pas[0] = sim->Q[i]->pas[0];
		ptr->pas[1] = sim->Q[i]->pas[1];
		ptr->pas[2] = sim->Q[i]->pas[2];
		if (sim->Q[i]->T != NULL) {
			ptr->T = dm_dup(sim->Q[i]->T);
		} else {
			ptr->T = NULL;
		}
		if (sim->Q[i]->Ta != NULL) {
			ptr->Ta = dm_dup(sim->Q[i]->Ta);
		} else {
			ptr->Ta = NULL;
		}
		if (sim->Q[i]->Tb != NULL) {
			ptr->Tb = dm_dup(sim->Q[i]->Tb);
		} else {
			ptr->Tb = NULL;
		}
		ptr->Rmol = cv_dup(sim->Q[i]->Rmol);
	}

	sim2->MIX = (Mixing**)malloc(sim->nMIX*sizeof(Mixing*));
	for (i=0; i<sim2->nMIX; i++) {
		Mixing *ptr = sim2->MIX[i] = (Mixing*)malloc(sizeof(Mixing));
		ptr->type = sim->MIX[i]->type;
		ptr->couple[0] = sim->MIX[i]->couple[0];
		ptr->couple[1] = sim->MIX[i]->couple[1];
		ptr->idx = sim->MIX[i]->idx;
		ptr->qidx = sim->MIX[i]->qidx;
		ptr->T = dm_dup(sim->MIX[i]->T);
		if (sim->MIX[i]->Ta != NULL) {
			ptr->Ta = blk_dm_dup(sim->MIX[i]->Ta);
		} else {
			ptr->Ta = NULL;
		}
		if (sim->MIX[i]->Tb != NULL) {
			ptr->Tb = blk_dm_dup(sim->MIX[i]->Tb);
		} else {
			ptr->Tb = NULL;
		}
	}

	sim2->ss=(SpinSys*)malloc(sizeof(SpinSys));
	if (sim2->ss == NULL) {
		fprintf(stderr,"Error: sim_duplicate - can not allocate SpinSys\n");
		exit(1);
	}
	sim2->ss->nchan = sim->ss->nchan;
	sim2->ss->matdim = sim->ss->matdim;
	sim2->ss->nspins = sim->ss->nspins;
	for (i=1; i<=sim2->ss->nspins; i++) sim2->ss->iso[i] = sim->ss->iso[i];
	for (i=1; i<=sim2->ss->nchan; i++) {
		sim2->ss->nchanelem[i] = sim->ss->nchanelem[i];
		strcpy(sim2->ss->channames[i],sim->ss->channames[i]);
		for (j=1; j<=sim2->ss->nchanelem[i]; j++) sim2->ss->chan[i][j] = sim->ss->chan[i][j];
	}

	return sim2;
}
***/

Sim_wsp * wsp_initialize(Sim_info *s)
{
	Sim_wsp * wsp;
	int i, j, k, nchan, nspins, matdim = s->matdim;
	mat_double *dmx;

	DEBUGPRINT("wsp_initialize start\n");
	wsp = (Sim_wsp*)malloc(sizeof(Sim_wsp));
	if (!wsp) {
		fprintf(stderr,"error: unable to allocate workspace structure\n");
		exit(1);
	}
	/* fid */
	wsp->fid = complx_vector(s->ntot);
	cv_zero(wsp->fid);
	wsp->curr_nsig = 0;
	/* interaction data */
	wsp->CS_Rlab = (complx**)malloc(sizeof(complx*)*s->nCS);
	wsp->CS_Rrot = (complx**)malloc(sizeof(complx*)*s->nCS);
	for (i=0; i<s->nCS; i++) {
		if (s->CS[i]->Rmol == NULL) {
			wsp->CS_Rrot[i] = NULL;
			wsp->CS_Rlab[i] = NULL;
		} else {
			wsp->CS_Rrot[i] = complx_vector(5);
			wsp->CS_Rlab[i] = complx_vector(5);
		}
	}
	wsp->DD_Rlab = (complx**)malloc(sizeof(complx*)*s->nDD);
	wsp->DD_Rrot = (complx**)malloc(sizeof(complx*)*s->nDD);
	for (i=0; i<s->nDD; i++) {
		wsp->DD_Rrot[i] = complx_vector(5);
		wsp->DD_Rlab[i] = complx_vector(5);
	}
	wsp->Q_Rlab = (complx**)malloc(sizeof(complx*)*s->nQ);
	wsp->Q_Rrot = (complx**)malloc(sizeof(complx*)*s->nQ);
	for (i=0; i<s->nQ; i++) {
		wsp->Q_Rrot[i] = complx_vector(5);
		wsp->Q_Rlab[i] = complx_vector(5);
	}
	wsp->J_Rlab = (complx**)malloc(sizeof(complx*)*s->nJ);
	wsp->J_Rrot = (complx**)malloc(sizeof(complx*)*s->nJ);
	for (i=0; i<s->nJ; i++) {
		if (s->J[i]->Rmol == NULL) {
			wsp->J_Rrot[i] = NULL;
			wsp->J_Rlab[i] = NULL;
		} else {
			wsp->J_Rrot[i] = complx_vector(5);
			wsp->J_Rlab[i] = complx_vector(5);
		}
	}
	wsp->CS_used = (int*)malloc(sizeof(int)*s->nCS);
	wsp->DD_used = (int*)malloc(sizeof(int)*s->nDD);
	wsp->J_used = (int*)malloc(sizeof(int)*s->nJ);
	wsp->Q_used = (int*)malloc(sizeof(int)*s->nQ);
	wsp->spinused = NULL;

	/* rf channels and offsets */
	nchan = s->ss->nchan;
	nspins = s->ss->nspins;
	wsp->rfv = double_vector(nchan);
	wsp->phv = double_vector(nchan);
	wsp->inhom_offset = NULL;
	wsp->offset = NULL;
	wsp->chan_Ix = (mat_double**)malloc(sizeof(mat_double*)*(nchan+1));
	wsp->chan_Iy = (mat_double**)malloc(sizeof(mat_double*)*(nchan+1));
	wsp->chan_Iz = (mat_double**)malloc(sizeof(mat_double*)*(nchan+1));
	wsp->Ix = (mat_double**)malloc(sizeof(mat_double*)*(nspins+1));
	wsp->Iz = (mat_double**)malloc(sizeof(mat_double*)*(nspins+1));
	for (i=1; i<=nchan; i++) {
		if (s->sparse) {
			wsp->chan_Ix[i] = double_matrix(matdim,matdim,MAT_SPARSE,1,0);
			wsp->chan_Iy[i] = double_matrix(matdim,matdim,MAT_SPARSE,1,0);
			wsp->chan_Iz[i] = double_matrix(matdim,matdim,MAT_SPARSE_DIAG,1,s->basis);
		} else {
			wsp->chan_Ix[i] = double_matrix(matdim,matdim,MAT_DENSE,0,0);
			wsp->chan_Iy[i] = double_matrix(matdim,matdim,MAT_DENSE,0,0);
			wsp->chan_Iz[i] = double_matrix(matdim,matdim,MAT_DENSE_DIAG,1,s->basis);
		}
		dm_zero(wsp->chan_Ix[i]); dm_zero(wsp->chan_Iy[i]); dm_zero(wsp->chan_Iz[i]);
		for (j=1; j<=s->ss->nchanelem[i]; j++) {
			k = s->ss->chan[i][j];
			dmx = Ix_real(s,k); dm_addto(wsp->chan_Ix[i],dmx); free_double_matrix(dmx);
			dmx = Iz_ham(s,k); dm_addto(wsp->chan_Iz[i],dmx); free_double_matrix(dmx);
			dmx = Iy_real(s,k); dm_addto(wsp->chan_Iy[i],dmx); free_double_matrix(dmx);
			wsp->Ix[k] = NULL;
			wsp->Iz[k] = NULL; // will be filled in 'select' command, on demand
		}
		if (s->basis != 0) {
			dm_permute(wsp->chan_Ix[i],s->perm_table[0+s->basis*s->Nbasis]);
			dm_permute(wsp->chan_Iy[i],s->perm_table[0+s->basis*s->Nbasis]);
			wsp->chan_Ix[i]->basis = wsp->chan_Iy[i]->basis = s->basis;
			/* chan_Iz is already permuted */
		}
	}

	int dum_type;
	if (s->sparse) dum_type = MAT_SPARSE; else dum_type = MAT_DENSE;
	if (s->Hint_isdiag) dum_type = MAT_DENSE_DIAG;
	if (s->basis != 0) {
		wsp->ham_blk = create_blk_mat_double(matdim,LEN(s->dims_table[s->basis]),s->dims_table[s->basis],dum_type,s->basis);
	} else {
		wsp->ham_blk = create_blk_mat_double(matdim,1,NULL,dum_type,s->basis);
	}

	wsp->sumUph = double_matrix(matdim,matdim,MAT_DENSE_DIAG,1,s->basis);

	if (s->block_diag) {
		wsp->U = NULL;
		wsp->dU = NULL;
		wsp->tmpU = NULL;
		wsp->sumHrf = NULL;
	} else {
		//wsp->U = create_blk_mat_complx(matdim,1,NULL,dum_type,0);
		//wsp->dU = create_blk_mat_complx(matdim,1,NULL,dum_type,0);
		//wsp->tmpU = create_blk_mat_complx(matdim,1,NULL,dum_type,0);
		wsp->U = create_blk_mat_complx(matdim,1,NULL,MAT_DENSE_DIAG,0);
		wsp->dU = create_blk_mat_complx(matdim,1,NULL,MAT_DENSE_DIAG,0);
		wsp->tmpU = NULL;
		wsp->sumHrf = create_blk_mat_double(matdim,1,NULL,(s->sparse) ? MAT_SPARSE : MAT_DENSE,0);
	}

    /* other */
    wsp->t = wsp->tstart = 0.0;
    if (s->wr == 0.0)
      wsp->dtmax=1e15; /* dtmax is infinity in the static case. */
    else
      wsp->dtmax= 1.0; /* dt max is 1 usec per default in case of spinning*/
    for (i=0;i<ACQBLOCK_STO_END; i++) {
    	wsp->STO[i] = NULL;
    	wsp->matrix[i] = NULL;
    }
    wsp->fstart = s->fstart;
    wsp->fdetect = s->fdetect;
    wsp->sigma = NULL;
    wsp->Nacq = s->np;
    wsp->acqphase = 0;
	wsp->evalmode = EM_NORMAL;

	/* Optimal control related */
	if (OCpar.isinit) {
		wsp->OC_propstatus = 0;
		wsp->OC_mxpos = 0;
		for (i=0; i<MAXOCPROPS; i++) {
			wsp->OC_dens[i] = NULL;
			wsp->OC_props[i] = NULL;
			wsp->OC_mxcode[i] = NULL;
			for (j=0; j<MAXOCSHAPES*2; j++) {
				wsp->OC_deriv[i][j] = NULL;
			}
		}
	}

	/* parametrs for spinach */
    wsp->spinach_delay_Nsub = 0;
    wsp->spinach_delay_sub_idx = NULL;
    wsp->spinach_pulse_Nsub = 0;
    wsp->spinach_pulse_sub_idx = NULL;

    wsp->Dmol_rot = NULL;
    wsp->Hiso_off = NULL;
    for (i=0; i<5; i++) wsp->HQ_off[i] = NULL;
    wsp->Nint_off = 0;

    /* when block-diagonalization will be used */
    wsp->Hiso = s->Hiso;
    for (i=0; i<5; i++) wsp->HQ[i] = s->HQ[i];
    if (s->nQ != 0) {
    	wsp->QTa = (mat_double**)malloc(s->nQ*sizeof(mat_double*));
    	wsp->QTb = (mat_double**)malloc(s->nQ*sizeof(mat_double*));
    	wsp->QT3a = (mat_double**)malloc(s->nQ*sizeof(mat_double*));
    	wsp->QT3b = (mat_double**)malloc(s->nQ*sizeof(mat_double*));
    	wsp->QT3c = (mat_double**)malloc(s->nQ*sizeof(mat_double*));
    	for (i=0; i<s->nQ; i++) {
    		wsp->QTa[i] = s->Q[i]->Ta;
    		wsp->QTb[i] = s->Q[i]->Tb;
    		wsp->QT3a[i] = s->Q[i]->T3a;
    		wsp->QT3b[i] = s->Q[i]->T3b;
    		wsp->QT3c[i] = s->Q[i]->T3c;
    	}
    } else {
    	wsp->QTa = NULL;
    	wsp->QTb = NULL;
    	wsp->QT3a = NULL;
    	wsp->QT3b = NULL;
    	wsp->QT3c = NULL;
    }
    if (s->nMIX != 0) {
    	wsp->MT  = (mat_double**)malloc(s->nMIX*sizeof(mat_double*));
    	wsp->MTa = (blk_mat_double**)malloc(s->nMIX*sizeof(blk_mat_double*));
    	wsp->MTb = (blk_mat_double**)malloc(s->nMIX*sizeof(blk_mat_double*));
    	for (i=0; i<s->nMIX; i++) {
    		wsp->MT[i]  = s->MIX[i]->T;
    		wsp->MTa[i] = s->MIX[i]->Ta;
    		wsp->MTb[i] = s->MIX[i]->Tb;
    	}
    } else {
    	wsp->MT  = NULL;
    	wsp->MTa = NULL;
    	wsp->MTb = NULL;
    }

    /* when averaging_file was added */
    wsp->CS = (Shift **)malloc(s->nCS*sizeof(Shift*));
    for (i=0; i<s->nCS; i++) wsp->CS[i] = s->CS[i];
    wsp->DD = (Dipole **)malloc(s->nDD*sizeof(Dipole*));
    for (i=0; i<s->nDD; i++) wsp->DD[i] = s->DD[i];
    wsp->J = (Jcoupling **)malloc(s->nJ*sizeof(Jcoupling*));
    for (i=0; i<s->nJ; i++) wsp->J[i] = s->J[i];
    wsp->Q = (Quadrupole **)malloc(s->nQ*sizeof(Quadrupole*));
    for (i=0; i<s->nQ; i++) wsp->Q[i] = s->Q[i];

    wsp->dw = 1.0e6/s->sw;

    wsp->FWTASG_irow = NULL;
    wsp->FWTASG_icol = NULL;

    // labframe
    if (s->labframe == 1) {
    	wsp->Hlab = complx_matrix(s->matdim,s->matdim,MAT_DENSE,0,s->basis);
    	wsp->Hrflab = complx_matrix(s->matdim,s->matdim,MAT_DENSE,0,s->basis);
    } else {
    	wsp->Hlab = NULL;
    	wsp->Hrflab = NULL;
    }

	DEBUGPRINT("wsp_initialize end\n");

	return wsp;
}


void wsp_destroy(Sim_info *s, Sim_wsp *wsp)
{
	int i, j;

    if (wsp->fdetect != s->fdetect) free_complx_matrix(wsp->fdetect);
    if (wsp->fstart != s->fstart) free_complx_matrix(wsp->fstart);
    if (wsp->sigma) free_complx_matrix(wsp->sigma);
    free_complx_vector(wsp->fid);

	for (i=0; i<s->nCS; i++) {
		free_complx_vector(wsp->CS_Rlab[i]);
		free_complx_vector(wsp->CS_Rrot[i]);
	}
	free((char*)(wsp->CS_Rlab));
	free((char*)(wsp->CS_Rrot));
	free((char*)(wsp->CS_used));
	for (i=0; i<s->nDD; i++) {
		free_complx_vector(wsp->DD_Rlab[i]);
		free_complx_vector(wsp->DD_Rrot[i]);
	}
	free((char*)(wsp->DD_Rlab));
	free((char*)(wsp->DD_Rrot));
	free((char*)(wsp->DD_used));
	for (i=0; i<s->nQ; i++) {
		free_complx_vector(wsp->Q_Rlab[i]);
		free_complx_vector(wsp->Q_Rrot[i]);
	}
	free((char*)(wsp->Q_Rlab));
	free((char*)(wsp->Q_Rrot));
	free((char*)(wsp->Q_used));
	for (i=0; i<s->nJ; i++) {
		free_complx_vector(wsp->J_Rlab[i]);
		free_complx_vector(wsp->J_Rrot[i]);
	}
	free((char*)(wsp->J_Rlab));
	free((char*)(wsp->J_Rrot));
	free((char*)(wsp->J_used));

    if (wsp->inhom_offset) free_double_vector(wsp->inhom_offset);
    if (wsp->offset) free_double_vector(wsp->offset);
    free_double_vector(wsp->rfv);
    free_double_vector(wsp->phv);
    for (i=1; i<=s->ss->nchan; i++) {
    	free_double_matrix(wsp->chan_Ix[i]);
    	free_double_matrix(wsp->chan_Iy[i]);
    	free_double_matrix(wsp->chan_Iz[i]);
    	for (j=1; j<=s->ss->nchanelem[i]; j++) {
    		if (wsp->Ix[s->ss->chan[i][j]] != NULL) free_double_matrix(wsp->Ix[s->ss->chan[i][j]]);
    		if (wsp->Iz[s->ss->chan[i][j]] != NULL) free_double_matrix(wsp->Iz[s->ss->chan[i][j]]);
    	}
    }
    free(wsp->chan_Ix);
    free(wsp->chan_Iy);
    free(wsp->chan_Iz);
    free(wsp->Ix);
    free(wsp->Iz);
   	if (wsp->ham_blk) free_blk_mat_double(wsp->ham_blk);
   	if (wsp->U) free_blk_mat_complx(wsp->U);
   	if (wsp->dU) free_blk_mat_complx(wsp->dU);
   	if (wsp->tmpU) free_blk_mat_complx(wsp->tmpU);
   	if (wsp->sumHrf) free_blk_mat_double(wsp->sumHrf);
   	if (wsp->sumUph) free_double_matrix(wsp->sumUph);

    for (i=0;i<ACQBLOCK_STO_END; i++) {
    	if (wsp->STO[i]) free_blk_mat_complx(wsp->STO[i]);
    	if (wsp->matrix[i]) free_complx_matrix(wsp->matrix[i]);
    }
    if (wsp->spinused) free_int_vector(wsp->spinused);

	/* Optimal control related */
    if (OCpar.isinit) {
    	for (i=0; i<MAXOCPROPS; i++) {
    		if (wsp->OC_dens[i]) free_complx_matrix(wsp->OC_dens[i]);
    		if (wsp->OC_props[i]) free_blk_mat_complx(wsp->OC_props[i]);
    		if (wsp->OC_mxcode[i]) free(wsp->OC_mxcode[i]);
			for (j=0; j<MAXOCSHAPES*2; j++) {
				if (wsp->OC_deriv[i][j] != NULL) free_complx_matrix(wsp->OC_deriv[i][j]);
			}

    	}
    }
	/* parametrs for spinach */
    if (wsp->spinach_delay_Nsub > 0) {
    	for (i=0; i<wsp->spinach_delay_Nsub; i++) free_int_vector(wsp->spinach_delay_sub_idx[i]);
    }
    if (wsp->spinach_pulse_Nsub > 0) {
    	for (i=0; i<wsp->spinach_pulse_Nsub; i++) free_int_vector(wsp->spinach_pulse_sub_idx[i]);
    }
    free(wsp->spinach_delay_sub_idx);
    free(wsp->spinach_pulse_sub_idx);

    if (wsp->Dmol_rot) free_complx_matrix(wsp->Dmol_rot);
    if (wsp->Hiso_off) free_blk_mat_double(wsp->Hiso_off);
    for (i=0; i<5; i++) {
    	if (wsp->HQ_off[i]) free_blk_mat_double(wsp->HQ_off[i]);
    }

    /* when block-diagonalization will be used */
    if (wsp->Hiso != s->Hiso) free_blk_mat_double(wsp->Hiso);
    for (i=0; i<5; i++) {
    	if (wsp->HQ[i] != s->HQ[i]) free_blk_mat_double(wsp->HQ[i]);
    }
    if (s->nQ != 0) {
    	for (i=0; i<s->nQ; i++) {
    		if (wsp->QTa[i] != s->Q[i]->Ta) free_double_matrix(wsp->QTa[i]);
    		if (wsp->QTb[i] != s->Q[i]->Tb) free_double_matrix(wsp->QTb[i]);
    		if (wsp->QT3a[i] != s->Q[i]->T3a) free_double_matrix(wsp->QT3a[i]);
    		if (wsp->QT3b[i] != s->Q[i]->T3b) free_double_matrix(wsp->QT3b[i]);
    		if (wsp->QT3c[i] != s->Q[i]->T3c) free_double_matrix(wsp->QT3c[i]);
    	}
    	free(wsp->QTa);
    	free(wsp->QTb);
    	free(wsp->QT3a);
    	free(wsp->QT3b);
    	free(wsp->QT3c);
    }
    if (s->nMIX != 0) {
    	for (i=0; i<s->nMIX; i++) {
    		if (wsp->MT[i]  != s->MIX[i]->T)  free_double_matrix(wsp->MT[i]);
    		if (wsp->MTa[i] != s->MIX[i]->Ta) free_blk_mat_double(wsp->MTa[i]);
    		if (wsp->MTb[i] != s->MIX[i]->Tb) free_blk_mat_double(wsp->MTb[i]);
    	}
    	free(wsp->MT);
    	free(wsp->MTa);
    	free(wsp->MTb);
    }

    /* when averaging_file was added */
    for (i=0; i<s->nCS; i++) {
    	if (wsp->CS[i] != s->CS[i]) {
    		free_complx_vector(wsp->CS[i]->Rmol);
    		if (wsp->CS[i]->T != s->CS[i]->T) free_double_matrix(wsp->CS[i]->T);
    		free((char *)(wsp->CS[i]));
    	}
    }
    free((char *)(wsp->CS)); wsp->CS = NULL;
    for (i=0; i<s->nDD; i++) {
    	if (wsp->DD[i] != s->DD[i]) {
    		free_complx_vector(wsp->DD[i]->Rmol);
    		if (wsp->DD[i]->blk_T != s->DD[i]->blk_T) free_blk_mat_double(wsp->DD[i]->blk_T);
    		free((char *)(wsp->DD[i]));
    	}
    }
    free((char *)(wsp->DD)); wsp->DD = NULL;
    for (i=0; i<s->nJ; i++) {
    	if (wsp->J[i] != s->J[i]) {
    		free_complx_vector(wsp->J[i]->Rmol);
    		if (wsp->J[i]->blk_Tiso != s->J[i]->blk_Tiso) free_blk_mat_double(wsp->J[i]->blk_Tiso);
    		if (wsp->J[i]->blk_T != s->J[i]->blk_T) free_blk_mat_double(wsp->J[i]->blk_T);
    		free((char *)(wsp->J[i]));
    	}
    }
    free((char *)(wsp->J)); wsp->J = NULL;
    for (i=0; i<s->nQ; i++) {
    	if (wsp->Q[i] != s->Q[i]) {
    		free_complx_vector(wsp->Q[i]->Rmol);
    		if (wsp->Q[i]->T != s->Q[i]->T) free_double_matrix(wsp->Q[i]->T);
    		free((char *)(wsp->Q[i]));
    	}
    }
    free((char *)(wsp->Q)); wsp->Q = NULL;

    if (wsp->FWTASG_irow != NULL) {
    	free(wsp->FWTASG_irow);
    	wsp->FWTASG_irow = NULL;
    }
    if (wsp->FWTASG_icol != NULL) {
    	free(wsp->FWTASG_icol);
    	wsp->FWTASG_icol = NULL;
    }

    // labframe
    if (wsp->Hlab != NULL) {
    	free_complx_matrix(wsp->Hlab);
    	wsp->Hlab = NULL;
    }
    if (wsp->Hrflab != NULL) {
    	free_complx_matrix(wsp->Hrflab);
    	wsp->Hrflab = NULL;
    }

}

void store_sim_pointers(Tcl_Interp* interp, Sim_info* sim, Sim_wsp * wsp)
{
	//char buf1[2*sizeof(Sim_info*)+3];
	//char buf2[2*sizeof(Sim_wsp*)+3];
	//sprintf(buf1,"%p",sim);
	//sprintf(buf2,"%p",wsp);
	//Tcl_SetVar(interp,"__sim_info_ptr",buf1,TCL_GLOBAL_ONLY);
	//Tcl_SetVar(interp,"__sim_wsp_ptr",buf2,TCL_GLOBAL_ONLY);
	//DEBUGPRINT("store_sim_pointers done (sim=%p, wsp=%p)\n",sim,wsp);

	char *buf = (char*)malloc(2*(sizeof(Sim_info)+sizeof(Sim_wsp))+3);
	sprintf(buf,"%p",sim);
	Tcl_SetVar(interp,"__sim_info_ptr",buf,TCL_GLOBAL_ONLY);
	sprintf(buf,"%p",wsp);
	Tcl_SetVar(interp,"__sim_wsp_ptr",buf,TCL_GLOBAL_ONLY);
	free(buf);
}

void read_sim_pointers(Tcl_Interp* interp, Sim_info **sim, Sim_wsp **wsp)
{
	char *buf;

	buf = Tcl_GetVar(interp,"__sim_info_ptr",TCL_GLOBAL_ONLY);
	if (buf == NULL) {
		fprintf(stderr,"Error: can not read __sim_info_ptr\n");
		exit(1);
	}
	sscanf(buf,"%p",sim);
	buf = Tcl_GetVar(interp,"__sim_wsp_ptr",TCL_GLOBAL_ONLY);
	if (buf == NULL) {
		fprintf(stderr,"Error: can not read __sim_wsp_ptr\n");
		exit(1);
	}
	sscanf(buf,"%p",wsp);
	//DEBUGPRINT("read_sim_pointers done (sim=%p, wsp=%p)\n",*sim,*wsp);
}

/****
 *
 */
void sim_prepare_interpol(Sim_info *sim)
{
	int i;
	switch (sim->interpolation) {
	case INTERPOL_NOT_USED:
		break;
	case INTERPOL_FWT_ALL:
	case INTERPOL_FWT_LAM:
		sim->FWTASG_nnz = int_vector(LEN(sim->crdata));
		sim->FWT_lam = (complx*)malloc(LEN(sim->crdata)*sim->matdim*sizeof(complx));
		sim->FWT_frs = (complx**)malloc(LEN(sim->crdata)*sizeof(complx*));
		sim->FWT_irow = (int**)malloc(LEN(sim->crdata)*sizeof(int*));
		sim->FWT_icol = (int**)malloc(LEN(sim->crdata)*sizeof(int*));
		if (sim->FWT_lam == NULL || sim->FWT_frs==NULL || sim->FWT_irow==NULL || sim->FWT_icol==NULL) {
			fprintf(stderr,"Error: no more memory for sim->FWT_xxx structures\n");
			exit(1);
		}
		for (i=0; i<LEN(sim->crdata); i++) {
			sim->FWT_frs[i] = NULL;
			sim->FWT_irow[i] = NULL;
			sim->FWT_icol[i] = NULL;
		}
		break;
	case INTERPOL_ASG:
		sim->FWTASG_nnz = int_vector(LEN(sim->crdata));
		sim->ASG_freq = (double**)malloc(LEN(sim->crdata)*sizeof(double*));
		sim->ASG_ampl = (complx**)malloc(LEN(sim->crdata)*sizeof(complx*));
		if (sim->ASG_freq == NULL || sim->ASG_ampl == NULL) {
			fprintf(stderr,"Error: no more memory for sim->ASG_xxx structures\n");
			exit(1);
		}
		for (i=0; i<LEN(sim->crdata); i++) {
			sim->ASG_freq[i] = NULL;
			sim->ASG_ampl[i] = NULL;
		}
		break;
	case INTERPOL_FWTASG_ALL:
	case INTERPOL_FWTASG_LAM:
		sim->FWTASG_nnz = int_vector(LEN(sim->crdata));
		sim->FWT_lam = (complx*)malloc(LEN(sim->crdata)*sim->matdim*sizeof(complx));
		sim->FWT_frs = (complx**)malloc(LEN(sim->crdata)*sizeof(complx*));
		sim->FWT_irow = (int**)malloc(LEN(sim->crdata)*sizeof(int*));
		sim->FWT_icol = (int**)malloc(LEN(sim->crdata)*sizeof(int*));
		if (sim->FWT_lam == NULL || sim->FWT_frs==NULL || sim->FWT_irow==NULL || sim->FWT_icol==NULL) {
			fprintf(stderr,"Error: no more memory for sim->FWT_xxx structures\n");
			exit(1);
		}
		for (i=0; i<LEN(sim->crdata); i++) {
			sim->FWT_frs[i] = NULL;
			sim->FWT_irow[i] = NULL;
			sim->FWT_icol[i] = NULL;
		}
		sim->ASG_freq = (double**)malloc(LEN(sim->targetcrdata)*sizeof(double*));
		sim->ASG_ampl = (complx**)malloc(LEN(sim->targetcrdata)*sizeof(complx*));
		if (sim->ASG_freq == NULL || sim->ASG_ampl == NULL) {
			fprintf(stderr,"Error: no more memory for sim->ASG_xxx structures\n");
			exit(1);
		}
		for (i=0; i<LEN(sim->targetcrdata); i++) {
			sim->ASG_freq[i] = NULL;
			sim->ASG_ampl[i] = NULL;
		}
		break;
	}
}

void sim_preempty_interpol(Sim_info *sim)
{
	int i;
	complx *zp;

	switch (sim->interpolation) {
	case INTERPOL_NOT_USED:
		break;
	case INTERPOL_FWT_ALL:
		if (LEN(sim->FWTASG_nnz) != LEN(sim->crdata)) {
			free_int_vector(sim->FWTASG_nnz);
			sim->FWTASG_nnz = int_vector(LEN(sim->crdata));
		}
		zp = (complx*) realloc(sim->FWT_lam, LEN(sim->crdata)*sim->matdim*sizeof(complx));
		if (zp == NULL) {
			fprintf(stderr,"Error: sim_preempty_interpol - reallocation failure\n");
			exit(1);
		}
		sim->FWT_lam = zp;
		for (i=0; i<LEN(sim->targetcrdata); i++) {
			if (sim->FWT_frs[i] != NULL) {
				free(sim->FWT_frs[i]);
				sim->FWT_frs[i] = NULL;
			}
		}
		for (i=0; i<LEN(sim->crdata); i++) { // _ALL done for Udiag only -> irow and icol are not changed
			if (sim->FWT_irow[i] != NULL) {
				free(sim->FWT_irow[i]);
				sim->FWT_irow[i] = NULL;
			}
			if (sim->FWT_icol[i] != NULL) {
				free(sim->FWT_icol[i]);
				sim->FWT_icol[i] = NULL;
			}
		}
		break;
	case INTERPOL_FWT_LAM:
		if (LEN(sim->FWTASG_nnz) != LEN(sim->crdata)) {
			free_int_vector(sim->FWTASG_nnz);
			sim->FWTASG_nnz = int_vector(LEN(sim->crdata));
		}
		zp = (complx*) realloc(sim->FWT_lam, LEN(sim->crdata)*sim->matdim*sizeof(complx));
		if (zp == NULL) {
			fprintf(stderr,"Error: sim_preempty_interpol - reallocation failure\n");
			exit(1);
		}
		sim->FWT_lam = zp;
		for (i=0; i<LEN(sim->crdata); i++) {
			if (sim->FWT_frs[i] != NULL) {
				free(sim->FWT_frs[i]);
				sim->FWT_frs[i] = NULL;
			}
			if (sim->FWT_irow[i] != NULL) {
				free(sim->FWT_irow[i]);
				sim->FWT_irow[i] = NULL;
			}
			if (sim->FWT_icol[i] != NULL) {
				free(sim->FWT_icol[i]);
				sim->FWT_icol[i] = NULL;
			}
		}
		break;
	case INTERPOL_ASG:
		for (i=0; i<LEN(sim->crdata); i++) {
			if (sim->ASG_freq[i] != NULL) {
				free(sim->ASG_freq[i]);
				sim->ASG_freq[i] = NULL;
			}
			if (sim->ASG_ampl[i] != NULL) {
				free(sim->ASG_ampl[i]);
				sim->ASG_ampl[i] = NULL;
			}
		}
		break;
	case INTERPOL_FWTASG_ALL:
		for (i=LEN(sim->crdata); i<LEN(sim->targetcrdata); i++) {
			if (sim->FWT_frs[i] != NULL) {
				free(sim->FWT_frs[i]);
				sim->FWT_frs[i] = NULL;
			}
		}
	case INTERPOL_FWTASG_LAM:
		if (LEN(sim->FWTASG_nnz) != LEN(sim->crdata)) {
			free_int_vector(sim->FWTASG_nnz);
			sim->FWTASG_nnz = int_vector(LEN(sim->crdata));
		}
		zp = (complx*) realloc(sim->FWT_lam, LEN(sim->crdata)*sim->matdim*sizeof(complx));
		if (zp == NULL) {
			fprintf(stderr,"Error: sim_preempty_interpol - reallocation failure (FWTASG)\n");
			exit(1);
		}
		sim->FWT_lam = zp;
		for (i=0; i<LEN(sim->crdata); i++) {
			if (sim->FWT_frs[i] != NULL) {
				free(sim->FWT_frs[i]);
				sim->FWT_frs[i] = NULL;
			}
			if (sim->FWT_irow[i] != NULL) {
				free(sim->FWT_irow[i]);
				sim->FWT_irow[i] = NULL;
			}
			if (sim->FWT_icol[i] != NULL) {
				free(sim->FWT_icol[i]);
				sim->FWT_icol[i] = NULL;
			}
		}
		for (i=0; i<LEN(sim->targetcrdata); i++) {
			if (sim->ASG_freq[i] != NULL) {
				free(sim->ASG_freq[i]);
				sim->ASG_freq[i] = NULL;
			}
			if (sim->ASG_ampl[i] != NULL) {
				free(sim->ASG_ampl[i]);
				sim->ASG_ampl[i] = NULL;
			}
		}
		break;
	}
}

/****
 * calculates pulse sequence response without any interpolations
 */
int sim_calcfid(Sim_info *s, Sim_wsp *wsp)
{
  int ig;

  wsp->cryst.gamma += s->gamma_zero;
  ham_rotate(s,wsp);
  cv_zero(wsp->fid);

  switch (s->imethod) {
  case M_GCOMPUTE_TIME:
      if ( s->relax ) {
         fprintf(stderr,"error: methods gcompute and igcompute not allowed with relaxation\n");
         exit(1);
      }
	  wsp->dt_gcompute = s->taur/(double)s->ngamma;
      gcompute_fid(s,wsp);
      if (s->conjugate_fid) cv_conj(wsp->fid);
	  break;
  case M_GCOMPUTE_FREQ:
      if ( s->relax ) {
         fprintf(stderr,"error: methods gcompute and igcompute not allowed with relaxation\n");
         exit(1);
      }
	  wsp->dt_gcompute = s->taur/(double)s->ngamma;
      gcompute_spc(s,wsp);
	  break;
  case M_DIRECT_TIME:
	  wsp->dt_gcompute = 1e99;
      for (ig=0;ig<s->ngamma;ig++) {
        if (s->wr == 0.0) {
          wsp->tstart = 0.0;
        } else {
	      wsp->tstart = ig/(double)s->ngamma*2.0e6*M_PI/s->wr;
        }
        wsp->ig = ig;
        direct_propagate(s,wsp);
        wsp->cryst.gamma += 360.0/(double)s->ngamma;
      }
      cv_muld(wsp->fid, 1.0/(double)(s->ngamma));
      if (s->conjugate_fid) cv_conj(wsp->fid);
      break;
  case M_DIRECT_FREQ:
	  wsp->dt_gcompute = 1e99;
      for (ig=0;ig<s->ngamma;ig++) {
        if (s->wr == 0.0) {
          wsp->tstart = 0.0;
        } else {
	      wsp->tstart = ig/(double)s->ngamma*2.0e6*M_PI/s->wr;
        }
        wsp->ig = ig;
        direct_propagate(s,wsp);
        wsp->cryst.gamma += 360.0/(double)s->ngamma;
      }
      cv_muld(wsp->fid, 1.0/(double)(s->ngamma));
      break;
  case M_SPINACH:
	  wsp->dt_gcompute = 1e99;
      for (ig=0;ig<s->ngamma;ig++) {
        if (s->wr == 0.0) {
          wsp->tstart = 0.0;
        } else {
	      wsp->tstart = ig/(double)s->ngamma*2.0e6*M_PI/s->wr;
        }
        wsp->ig = ig;
        direct_propagate(s,wsp);
        wsp->cryst.gamma += 360.0/(double)s->ngamma;
      }
      cv_muld(wsp->fid, 1.0/(double)(s->ngamma));
      if (s->conjugate_fid) cv_conj(wsp->fid);
      break;
  default:
	  fprintf(stderr,"Error: sim_calcfid - imethod (%d) not recognized\n",s->imethod);
	  exit(1);
  }

  return 0;
}

/****
 * calculates pulse sequence response data for interpolations
 */
int sim_calcfid_interpol(Sim_info *sim, Sim_wsp *wsp)
{
  if ( sim->relax ) {
     fprintf(stderr,"Error: relaxation not compatible with interpolation\n");
     exit(1);
  }

  wsp->cryst.gamma += sim->gamma_zero;
  ham_rotate(sim,wsp);

  switch (sim->imethod) {
  case M_GCOMPUTE_TIME:
      if ( !(sim->interpolation == INTERPOL_FWT_ALL || sim->interpolation == INTERPOL_FWT_LAM) ) {
    	  fprintf(stderr,"Error: incompatible interpolation method for time domain. Use FWT.\n");
    	  exit(1);
      }
	  wsp->dt_gcompute = sim->taur/(double)sim->ngamma;
      gcompute_FWTdata(sim,wsp);
	  break;
  case M_GCOMPUTE_FREQ:
	  wsp->dt_gcompute = sim->taur/(double)sim->ngamma;
      if (sim->interpolation == INTERPOL_FWT_ALL || sim->interpolation == INTERPOL_FWT_LAM || sim->interpolation == INTERPOL_FWTASG_ALL || sim->interpolation == INTERPOL_FWTASG_LAM) {
    	  gcompute_FWTdata(sim,wsp);
      } else if (sim->interpolation == INTERPOL_ASG) {
    	  gcompute_ASGdata(sim,wsp);
      } else {
    	  fprintf(stderr,"Error: incompatible interpolation method for frequency domain. Try FWT or ASG.\n");
    	  exit(1);
      }
	  break;
  case M_DIRECT_TIME:
	  if ( !(sim->interpolation == INTERPOL_FWT_ALL || sim->interpolation == INTERPOL_FWT_LAM) ) {
		  fprintf(stderr,"Error: incompatible interpolation method for time domain. Use FWT.\n");
		  exit(1);
	  }
	  wsp->dt_gcompute = 1e99;
	  fprintf(stderr,"\n\n...NOTHING DONE HERE!!!...\n\n");
	  exit(1);
      break;
  case M_DIRECT_FREQ:
	  wsp->dt_gcompute = 1e99;
	  // works only for static samples
	  if (sim->wr > SPINRATE_SMALL) {
		  fprintf(stderr,"Error: interpolation implemented only for static samples\n");
		  exit(1);
	  }
	  // static samples do not depend on gamma angle
	  direct_propagate(sim,wsp);
      break;
  case M_SPINACH:
	  fprintf(stderr,"Error: the method is not compatible with interpolations\n");
	  exit(1);
      break;
  default:
	  fprintf(stderr,"Error: sim_calcfid_interpol - imethod (%d) not recognized\n",sim->imethod);
	  exit(1);
  }

  return 0;
}

