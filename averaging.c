/*
 * averaging.c
 *
 *  Created on: 14.9.2011
 *      Author: Zdenek Tosner
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <tcl.h>
#include "averaging.h"
#include "defs.h"
#include "sim.h"
#include "spinsys.h"
#include "ham.h"
#include "blockdiag.h"
#include "wigner.h"
#include "cm.h"

#define MAXLINELEN 1024

void set_ave_val_type(Ave_elem *ave, int i, char *str)
{
	if (strcmp(str,"iso") == 0) {
		ave[i].val_type = 1;
		return;
	}
	if (strcmp(str,"aniso") == 0) {
		ave[i].val_type = 2;
		return;
	}
	if (strcmp(str,"eta") == 0) {
		ave[i].val_type = 3;
		return;
	}
	if (strcmp(str,"alpha") == 0) {
		ave[i].val_type = 4;
		return;
	}
	if (strcmp(str,"beta") == 0) {
		ave[i].val_type = 5;
		return;
	}
	if (strcmp(str,"gamma") == 0) {
		ave[i].val_type = 6;
		return;
	}

	ave[i].val_type = 0;
	return;
}

void read_averaging_file(Tcl_Interp *interp, Sim_info *sim, Ave_elem **aveptr, double **weightptr, int *Npar, int *Nval)
{
	FILE *fp;
	char fname[256], name[256], line[MAXLINELEN], **parsv, buf[8], *sptr;
	Tcl_Obj *obj, **objv;
	int i, j, k, N, N_from, N_to, Nlines;
	Ave_elem *ave;
	double dval, *ave_weight;
	const int verb = verbose & VERBOSE_AVERAGING_FILE;


	obj = Tcl_GetVar2Ex(interp,"par","averaging_file",TCL_GLOBAL_ONLY);
	if (obj == NULL) {
		/* NO averaging */
		*aveptr = NULL;
		ave_weight = (double*)malloc(sizeof(double));
		ave_weight[0] = 1.0;
		*weightptr = ave_weight;
		*Npar = 0;
		*Nval = 1;
		return;
	}
	if (Tcl_ListObjGetElements(interp,obj,&N,&objv) != TCL_OK) {
		fprintf(stderr,"Error: problem reading averaging_file definition in par:\n%s\n",Tcl_GetStringResult(interp));
		exit(1);
	}
	if (N == 1) {
		strcpy(name,Tcl_GetString(objv[0]));
		N_from = 1;
		N_to = -1;
	} else if (N == 3) {
		strcpy(name,Tcl_GetString(objv[0]));
		if (Tcl_GetIntFromObj(interp,objv[1],&N_from) != TCL_OK) {
			fprintf(stderr,"Error: problem reading averaging_file definition in par (from):\n%s\n",Tcl_GetStringResult(interp));
			exit(1);
		}
		if (Tcl_GetIntFromObj(interp,objv[2],&N_to) != TCL_OK) {
			fprintf(stderr,"Error: problem reading averaging_file definition in par (to):\n%s\n",Tcl_GetStringResult(interp));
			exit(1);
		}
		if (N_to < N_from) {
			fprintf(stderr,"Error: problem reading averaging_file definition in par:\n 'to' is smaller than 'from'\n");
			exit(1);
		}
	} else {
		fprintf(stderr,"Error: problem reading averaging_file definition in par:\nlist elements count mismatch\n");
		exit(1);
	}

	strcpy(fname,name);
#ifdef UNIX
   if (name[0] == '~') {
	   char* p=getenv("HOME");
	   if (p != NULL) {
		   strcpy(fname,p);
		   strcat(fname,&name[1]);
	   }
   }
#endif
   if (!(fp=fopen(fname,"r"))) {
	   strcat(fname,".ave");
	   fp=fopen(fname,"r");
   }
   if (!fp) {
	   fprintf(stderr,"error: unable to open averaging_file '%s'\n",fname);
	   exit(1);
   }

   /* scan file for number of lines with values */
   Nlines = -1;
   while ( fgets(line, MAXLINELEN, fp) ) {
      Nlines++;
   }

   /* first line defines parameter names and weight */
   fseek(fp, 0, SEEK_SET);
   fgets(line, MAXLINELEN, fp);
   if (Tcl_SplitList(interp,line,&N,&parsv) != TCL_OK) {
	   fprintf(stderr,"Error: can not decompose first line of averaging_file:\n%s\n",Tcl_GetStringResult(interp));
	   exit(1);
   }
   if (strcmp(parsv[N-1],"weight") != 0 ) {
	   fprintf(stderr,"Error: last parameter in averaging_file must be 'weight'\n");
	   exit(1);
   }
   /* debug output start */
   if (verb) {
	   printf("averaging_file number of value lines = %d\n",Nlines);
	   printf("               number of parameters  = %d\n",N);
	   for (i=0; i<N; i++) {
		   printf("                 %d) %s\n",i+1,parsv[i]);
	   }
   }
   /* debug output end */

   /* set range */
   if (N_from > Nlines) {
	   fprintf(stderr,"Error in averaging_file range: from=%d is too large (max=%d)\n",N_from,Nlines);
	   exit(1);
   }
   if (N_to > Nlines || N_to < 0) N_to = Nlines;
   Nlines = N_to - N_from + 1;

   /* analyze first line and allocate ave structure */
   ave = (Ave_elem *)malloc(N*sizeof(Ave_elem));
   ave_weight = (double *)malloc(Nlines*sizeof(double));
   for (i=0; i<N; i++) {
	   if (sscanf(parsv[i],"shift_%d_%s",&j,buf) == 2) {
		   /* shift must exist in sim!!! */
		   if (shift_exist(sim,j) < 0) {
			   fprintf(stderr,"Error reading averaging_file: evaluating '%s' but interaction 'shift %d' not properly defined in spinsys\n",parsv[i],j);
			   exit(1);
		   }
		   ave[i].nuc1 = j;
		   //printf("shift found: %d '%s'\n",ave[i].nuc1,buf);
		   ave[i].par_type = 1;
		   ave[i].txt = NULL;
		   set_ave_val_type(ave,i,buf);
		   if (ave[i].val_type == 0) {
			   fprintf(stderr,"Error in reading averaging_file: incorrect parameter '%s'\n",parsv[i]);
			   exit(1);
		   }
		   ave[i].data = (double*)malloc(Nlines*sizeof(double));
		   continue;
	   }
	   if (sscanf(parsv[i],"dipole_%d_%d_%s",&j,&k,buf) == 3) {
		   /* dipole must exist in sim!!! */
		   if (dipole_exist(sim,j,k) < 0) {
			   fprintf(stderr,"Error reading averaging_file: evaluating '%s' but interaction 'dipole %d %d' not properly defined in spinsys\n",parsv[i],j,k);
			   exit(1);
		   }
		   if (k<j) {
			   ave[i].nuc1 = k;
			   ave[i].nuc2 = j;
		   } else {
			   ave[i].nuc1 = j;
			   ave[i].nuc2 = k;
		   }
		   //printf("dipole found: %d %d '%s'\n",ave[i].nuc1,ave[i].nuc2,buf);
		   ave[i].par_type = 2;
		   ave[i].txt = NULL;
		   set_ave_val_type(ave,i,buf);
		   if (ave[i].val_type == 0 || ave[i].val_type == 1) {
			   fprintf(stderr,"Error in reading averaging_file: incorrect parameter '%s'\n",parsv[i]);
			   exit(1);
		   }
		   ave[i].data = (double*)malloc(Nlines*sizeof(double));
		   continue;
	   }
	   if (sscanf(parsv[i],"jcoupling_%d_%d_%s",&j,&k,buf) == 3) {
		   /* jcoupling must exist in sim!!! */
		   if (jcoupling_exist(sim,j,k) < 0) {
			   fprintf(stderr,"Error reading averaging_file: evaluating '%s' but interaction 'jcoupling %d %d' not properly defined in spinsys\n",parsv[i],j,k);
			   exit(1);
		   }
		   if (k<j) {
			   ave[i].nuc1 = k;
			   ave[i].nuc2 = j;
		   } else {
			   ave[i].nuc1 = j;
			   ave[i].nuc2 = k;
		   }
		   //printf("jcoupling found: %d %d '%s'\n",ave[i].nuc1,ave[i].nuc2,buf);
		   ave[i].par_type = 3;
		   ave[i].txt = NULL;
		   set_ave_val_type(ave,i,buf);
		   if (ave[i].val_type == 0) {
			   fprintf(stderr,"Error in reading averaging_file: incorrect parameter '%s'\n",parsv[i]);
			   exit(1);
		   }
		   ave[i].data = (double*)malloc(Nlines*sizeof(double));
		   continue;
	   }
	   if (sscanf(parsv[i],"quadrupole_%d_%s",&j,buf) == 2) {
		   /* quadrupole must exist in sim!!! */
		   if (quadrupole_exist(sim,j) < 0) {
			   fprintf(stderr,"Error reading averaging_file: evaluating '%s' but interaction 'quadrupole %d' not properly defined in spinsys\n",parsv[i],j);
			   exit(1);
		   }
		   ave[i].nuc1 = j;
		   //printf("quadrupole found: %d '%s'\n",ave[i].nuc1,buf);
		   ave[i].par_type = 4;
		   ave[i].txt = NULL;
		   set_ave_val_type(ave,i,buf);
		   if (ave[i].val_type == 0 || ave[i].val_type == 1) {
			   fprintf(stderr,"Error in reading averaging_file: incorrect parameter '%s'\n",parsv[i]);
			   exit(1);
		   }
		   ave[i].data = (double*)malloc(Nlines*sizeof(double));
		   continue;
	   }
	   if (sscanf(parsv[i],"dipole_ave_%d_%d_%s",&j,&k,buf) == 3) {
		   /* dipole must exist in sim!!! */
		   if (dipole_exist(sim,j,k) < 0) {
			   fprintf(stderr,"Error reading averaging_file: evaluating '%s' but interaction 'dipole %d %d' not properly defined in spinsys\n",parsv[i],j,k);
			   exit(1);
		   }
		   if (k<j) {
			   ave[i].nuc1 = k;
			   ave[i].nuc2 = j;
		   } else {
			   ave[i].nuc1 = j;
			   ave[i].nuc2 = k;
		   }
		   //printf("dipole found: %d %d '%s'\n",ave[i].nuc1,ave[i].nuc2,buf);
		   ave[i].par_type = 2;
		   ave[i].txt = NULL;
		   set_ave_val_type(ave,i,buf);
		   if (ave[i].val_type == 0 || ave[i].val_type == 1) {
			   fprintf(stderr,"Error in reading averaging_file: incorrect parameter '%s'\n",parsv[i]);
			   exit(1);
		   }
		   ave[i].data = (double*)malloc(Nlines*sizeof(double));
		   continue;
	   }
	   /* otherwise it is par(*) */
	   ave[i].par_type = 0;
	   ave[i].txt = strdup(parsv[i]);
	   ave[i].data = (double*)malloc(Nlines*sizeof(double));
   }
   Tcl_Free((char *)parsv);

   /* roll on to the first relevant line */
   for (i=1; i<N_from; i++) fgets(line, MAXLINELEN, fp);
   /* read in values */
   for (i=0; i<Nlines; i++) {
	   fgets(line, MAXLINELEN, fp);
	   if (Tcl_SplitList(interp,line,&k,&parsv) != TCL_OK) {
		   fprintf(stderr,"Error: can not decompose line %d of averaging_file:\n%s\n",i+N_from+1,Tcl_GetStringResult(interp));
		   exit(1);
	   }
	   if (k != N) {
		   fprintf(stderr,"Error in reading averaging_file: line %d has incorrect number of elements (%d, should be %d)\n",i+N_from+1,k,N);
		   exit(1);
	   }
	   for (j=0; j<N-1; j++) {
		   dval = strtod(parsv[j],&sptr);
		   if (sptr[0] == 'p') {
			   if (ave[j].par_type != 1 || ave[j].val_type > 2) {
				   fprintf(stderr,"Error reading averaging_file (line %d, col %d): 'p' can be used only for shift iso/aniso\n",i+N_from+1,j+1);
				   exit(1);
			   }
			   ave[j].data[i] = dval*fabs(ss_gamma(sim->ss,ave[j].nuc1)/ss_gamma1H()*sim->specfreq*1e-6);
			   if (verb) printf("\t file line %d, el %d : read %lg ppm, translated to %lg Hz\n",i+N_from+1,j+1,dval,ave[j].data[i]);
		   } else if (sptr[0] == '\0') {
			   if (verb) printf("\t file line %d, el %d : read %lg\n",i+N_from+1,j+1,dval);
			   ave[j].data[i] = dval;
		   } else {
			   fprintf(stderr,"Error reading averaging_file (line %d, col %d): conversion not valid for '%s'\n",i+N_from+1,j+1,parsv[j]);
			   exit(1);
		   }
	   }
	   /* the last value is weight */
	   if (sscanf(parsv[N-1],"%lg",&(ave_weight[i])) != 1) {
		   fprintf(stderr,"Error reading averaging_file: can not read weight on line %d\n",i+N_from+1);
		   exit(1);
	   }
	   if (verb) printf("\t file line %d, el %d : read weight %lg\n",i+N_from+1,N,ave_weight[i]);

	   Tcl_Free((char *)parsv);
   }
    fclose(fp);
   *Npar = N;
   *Nval = Nlines;
   *aveptr = ave;
   *weightptr = ave_weight;
}


void free_averaging_data(Ave_elem **avestruct,int Navepar)
{
	int i;
	Ave_elem *ave = *avestruct;
	if (ave != NULL) {
		for (i=0; i<Navepar; i++) free(ave[i].data);
		free(ave);
		*avestruct = NULL;
	}
}

double signum(double x) {
	if (x > 0) return 1.0;
	if (x < 0) return -1.0;
	return 0.0;
}

void set_averaging_parameters(Sim_info *sim,Sim_wsp *wsp,Ave_elem *ave_struct,int Navepar,int idx)
{
	int i, j, Naniso=0, Niso=0;

	if (Navepar == 0) return;

	for (i=0; i<Navepar; i++) {
		switch (ave_struct[i].par_type) {
		case 0: { // par(*)
			Tcl_Obj *obj = Tcl_NewDoubleObj(ave_struct[i].data[idx]);
			if (Tcl_SetVar2Ex(wsp->interp,"par",ave_struct[i].txt,obj,TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG) == NULL) {
				fprintf(stderr,"Error in setting parameter par(%s) from averaging_file:\n%s\n",ave_struct[i].txt,Tcl_GetStringResult(wsp->interp));
				exit(1);
			}
			break; }
		case 1: { // shift
			j = shift_exist(sim,ave_struct[i].nuc1);
			assert( j >= 0 );
			if (wsp->CS[j] == sim->CS[j]) {
				wsp->CS[j] = (Shift *)malloc(sizeof(Shift));
				wsp->CS[j]->nuc = sim->CS[j]->nuc;
				wsp->CS[j]->iso = sim->CS[j]->iso;
				wsp->CS[j]->delta = sim->CS[j]->delta;
				wsp->CS[j]->eta = sim->CS[j]->eta;
				wsp->CS[j]->T = sim->CS[j]->T;
				wsp->CS[j]->pas[0] = sim->CS[j]->pas[0];
				wsp->CS[j]->pas[1] = sim->CS[j]->pas[1];
				wsp->CS[j]->pas[2] = sim->CS[j]->pas[2];
			}
			switch (ave_struct[i].val_type) {
			case 1: Niso++; wsp->CS[j]->iso = -signum(ss_gamma(sim->ss,ave_struct[i].nuc1))*2.0*M_PI*ave_struct[i].data[idx]; break;
			case 2: Naniso++; wsp->CS[j]->delta = -signum(ss_gamma(sim->ss,ave_struct[i].nuc1))*2.0*M_PI*ave_struct[i].data[idx]; break;
			case 3: Naniso++; wsp->CS[j]->eta = ave_struct[i].data[idx]; break;
			case 4: wsp->CS[j]->pas[0] = ave_struct[i].data[idx]; break;
			case 5: wsp->CS[j]->pas[1] = ave_struct[i].data[idx]; break;
			case 6: wsp->CS[j]->pas[2] = ave_struct[i].data[idx]; break;
			}
			break; }
		case 2: { // dipole
			j = dipole_exist(sim,ave_struct[i].nuc1,ave_struct[i].nuc2);
			assert( j >= 0 );
			if (wsp->DD[j] == sim->DD[j]) {
				wsp->DD[j] = (Dipole *)malloc(sizeof(Dipole));
				wsp->DD[j]->nuc[0] = sim->DD[j]->nuc[0];
				wsp->DD[j]->nuc[1] = sim->DD[j]->nuc[1];
				wsp->DD[j]->delta = sim->DD[j]->delta;
				wsp->DD[j]->eta = sim->DD[j]->eta;
				wsp->DD[j]->blk_T = sim->DD[j]->blk_T;
				wsp->DD[j]->pas[0] = sim->DD[j]->pas[0];
				wsp->DD[j]->pas[1] = sim->DD[j]->pas[1];
				wsp->DD[j]->pas[2] = sim->DD[j]->pas[2];
			}
			Naniso++;
			switch (ave_struct[i].val_type) {
			case 2: wsp->DD[j]->delta = ave_struct[i].data[idx]; break;
			case 3: wsp->DD[j]->eta = ave_struct[i].data[idx]; break;
			case 4: wsp->DD[j]->pas[0] = ave_struct[i].data[idx]; break;
			case 5: wsp->DD[j]->pas[1] = ave_struct[i].data[idx]; break;
			case 6: wsp->DD[j]->pas[2] = ave_struct[i].data[idx]; break;
			}
			break; }
		case 3: { // jcoupling
			j = jcoupling_exist(sim,ave_struct[i].nuc1,ave_struct[i].nuc2);
			assert( j >= 0 );
			if (wsp->J[j] == sim->J[j]) {
				wsp->J[j] = (Jcoupling *)malloc(sizeof(Jcoupling));
				wsp->J[j]->nuc[0] = sim->J[j]->nuc[0];
				wsp->J[j]->nuc[1] = sim->J[j]->nuc[1];
				wsp->J[j]->iso = sim->J[j]->iso;
				wsp->J[j]->delta = sim->J[j]->delta;
				wsp->J[j]->eta = sim->J[j]->eta;
				wsp->J[j]->blk_Tiso = sim->J[j]->blk_Tiso;
				wsp->J[j]->blk_T = sim->J[j]->blk_T;
				wsp->J[j]->pas[0] = sim->J[j]->pas[0];
				wsp->J[j]->pas[1] = sim->J[j]->pas[1];
				wsp->J[j]->pas[2] = sim->J[j]->pas[2];
			}
			switch (ave_struct[i].val_type) {
			case 1: Niso++; wsp->J[j]->iso = 2.0*M_PI*ave_struct[i].data[idx]; break;
			case 2: Naniso++; wsp->J[j]->delta = 2.0*M_PI*ave_struct[i].data[idx]; break;
			case 3: Naniso++; wsp->J[j]->eta = ave_struct[i].data[idx]; break;
			case 4: wsp->J[j]->pas[0] = ave_struct[i].data[idx]; break;
			case 5: wsp->J[j]->pas[1] = ave_struct[i].data[idx]; break;
			case 6: wsp->J[j]->pas[2] = ave_struct[i].data[idx]; break;
			}
			break; }
		case 4: { // quadrupole
			j = quadrupole_exist(sim,ave_struct[i].nuc1);
			assert( j >= 0 );
			if (wsp->Q[j] == sim->Q[j]) {
				wsp->Q[j] = (Quadrupole *)malloc(sizeof(Quadrupole));
				wsp->Q[j]->nuc = sim->Q[j]->nuc;
				wsp->Q[j]->order = sim->Q[j]->order;
				wsp->Q[j]->wq = sim->Q[j]->wq;
				wsp->Q[j]->eta = sim->Q[j]->eta;
				wsp->Q[j]->w0 = sim->Q[j]->w0;
				wsp->Q[j]->T = sim->Q[j]->T;
				wsp->Q[j]->Ta = sim->Q[j]->Ta;
				wsp->Q[j]->Tb = sim->Q[j]->Tb;
				wsp->Q[j]->pas[0] = sim->Q[j]->pas[0];
				wsp->Q[j]->pas[1] = sim->Q[j]->pas[1];
				wsp->Q[j]->pas[2] = sim->Q[j]->pas[2];
			}
			Naniso++;
			switch (ave_struct[i].val_type) {
			case 2: wsp->Q[j]->wq = ave_struct[i].data[idx]*2*M_PI/(4.0*ss_qn(sim->ss,ave_struct[i].nuc1)*(2.0*ss_qn(sim->ss,ave_struct[i].nuc1)-1)); break;
			case 3: wsp->Q[j]->eta = ave_struct[i].data[idx]; break;
			case 4: wsp->Q[j]->pas[0] = ave_struct[i].data[idx]; break;
			case 5: wsp->Q[j]->pas[1] = ave_struct[i].data[idx]; break;
			case 6: wsp->Q[j]->pas[2] = ave_struct[i].data[idx]; break;
			}
			break; }
		}
	}

	/* now complete definitions with Rmol vectors and create new Hamiltonian matrices */
	if (wsp->Hiso != sim->Hiso) {
		free_blk_mat_double(wsp->Hiso);
	}
	wsp->Hiso = create_blk_mat_double_copy(sim->Hiso);
	blk_dm_copy(wsp->Hiso,sim->Hiso);
	if (sim->HQ[0] != NULL) {
		for (i=0; i<5; i++) {
			if (wsp->HQ[i] != sim->HQ[i]) free_blk_mat_double(wsp->HQ[i]);
			wsp->HQ[i] = create_blk_mat_double_copy(sim->HQ[i]);
			blk_dm_copy(wsp->HQ[i],sim->HQ[i]);
		}
	} else if (Naniso>0 && sim->block_diag) {
		sim->Hassembly = 1;
		for (i=0; i<5; i++) {
			if (wsp->HQ[i] != NULL) free_blk_mat_double(wsp->HQ[i]);
			wsp->HQ[i] = create_blk_mat_double_copy(sim->Hiso);
			blk_dm_zero(wsp->HQ[i]);
		}
	}
    for (i=0; i<sim->nCS; i++) {
    	if (wsp->CS[i] != sim->CS[i]) {
    		Shift *csptr = wsp->CS[i];
    		if (fabs(csptr->delta) > TINY ) {
    			csptr->Rmol = Dtensor2(SQRT2BY3*csptr->delta, csptr->eta);
    			csptr->Rmol = wig2roti(csptr->Rmol, csptr->pas[0], csptr->pas[1], csptr->pas[2]);
    			if (wsp->CS_Rrot[i] == NULL) {
    				wsp->CS_Rrot[i] = complx_vector(5);
    				wsp->CS_Rlab[i] = complx_vector(5);
    			}
    		} else {
    			csptr->Rmol = NULL;
    			if (wsp->CS_Rrot[i] != NULL) {
    				free_complx_vector(wsp->CS_Rrot[i]);
    				wsp->CS_Rrot[i] = NULL;
    				free_complx_vector(wsp->CS_Rlab[i]);
    				wsp->CS_Rlab[i] = NULL;
    			}
    		}
    		if (csptr->T == NULL) csptr->T = Iz_ham(sim,csptr->nuc);
    		blk_dm_multod_diag(wsp->Hiso,csptr->T,csptr->iso - sim->CS[i]->iso);
    		if (sim->Hassembly) {
    			complx *zdum = complx_vector(5);
    			cv_zero(zdum);
    			if (sim->CS[i]->Rmol != NULL) cv_multod(zdum,sim->CS[i]->Rmol,-1.0);
    			if (csptr->Rmol != NULL) cv_multod(zdum,csptr->Rmol,1.0);
    			blk_dm_multod_diag(wsp->HQ[0],csptr->T,zdum[3].re);
				blk_dm_multod_diag(wsp->HQ[1],csptr->T,zdum[4].re);
				blk_dm_multod_diag(wsp->HQ[2],csptr->T,zdum[4].im);
				blk_dm_multod_diag(wsp->HQ[3],csptr->T,zdum[5].re);
				blk_dm_multod_diag(wsp->HQ[4],csptr->T,zdum[5].im);
				free_complx_vector(zdum);
    		}
    	}
    }
    for (i=0; i<sim->nDD; i++) {
    	if (wsp->DD[i] != sim->DD[i]) {
    		Dipole * ddptr = wsp->DD[i];
        	ddptr->Rmol = Dtensor2(4.0*M_PI*ddptr->delta, ddptr->eta);
        	ddptr->Rmol = wig2roti(ddptr->Rmol, ddptr->pas[0], ddptr->pas[1], ddptr->pas[2]);
        	if (ddptr->blk_T == NULL) {
        		if (ss_issame(sim->ss,ddptr->nuc[0],ddptr->nuc[1])) {
        			ddptr->blk_T = T20(sim,ddptr->nuc[0],ddptr->nuc[1]);
        		} else {
        			ddptr->blk_T = IzIz_sqrt2by3(sim,ddptr->nuc[0],ddptr->nuc[1]);
        		}
        	}
        	if (sim->Hassembly) {
        		blk_dm_multod(wsp->HQ[0],ddptr->blk_T,ddptr->Rmol[3].re - sim->DD[i]->Rmol[3].re);
        		blk_dm_multod(wsp->HQ[1],ddptr->blk_T,ddptr->Rmol[4].re - sim->DD[i]->Rmol[4].re);
        		blk_dm_multod(wsp->HQ[2],ddptr->blk_T,ddptr->Rmol[4].im - sim->DD[i]->Rmol[4].im);
        		blk_dm_multod(wsp->HQ[3],ddptr->blk_T,ddptr->Rmol[5].re - sim->DD[i]->Rmol[5].re);
        		blk_dm_multod(wsp->HQ[4],ddptr->blk_T,ddptr->Rmol[5].im - sim->DD[i]->Rmol[5].im);
        	}
    	}
    }
    for (i=0; i<sim->nJ; i++) {
    	if (wsp->J[i] != sim->J[i]) {
    		Jcoupling *jptr = wsp->J[i];
    		if (fabs(jptr->delta) > TINY) {
    			jptr->Rmol = Dtensor2(jptr->delta*2.0, jptr->eta);
    			jptr->Rmol = wig2roti(jptr->Rmol,jptr->pas[0],jptr->pas[1],jptr->pas[2]);
    			if (wsp->J_Rrot[i] == NULL) {
    				wsp->J_Rrot[i] = complx_vector(5);
    				wsp->J_Rlab[i] = complx_vector(5);
    			}
    		} else {
    			jptr->Rmol = NULL;
    			if (wsp->J_Rrot[i] != NULL) {
    				free_complx_vector(wsp->J_Rrot[i]);
    				wsp->J_Rrot[i] = NULL;
    				free_complx_vector(wsp->J_Rlab[i]);
    				wsp->J_Rlab[i] = NULL;
    			}
    		}
    		if (jptr->blk_Tiso == NULL) {
    			if (ss_issame(sim->ss,jptr->nuc[0],jptr->nuc[1])) {
    				jptr->blk_Tiso = II(sim,jptr->nuc[0],jptr->nuc[1]);
    			} else {
    				jptr->blk_Tiso = IzIz_sqrt2by3(sim,jptr->nuc[0],jptr->nuc[1]);
    			}
    		}
    		double corr_factor;
			if (ss_issame(sim->ss,jptr->nuc[0],jptr->nuc[1])) {
				corr_factor = 1.0;
			} else {
				corr_factor = 1.0/SQRT2BY3;
			}
    		blk_dm_multod(wsp->Hiso,jptr->blk_Tiso,(jptr->iso - sim->J[i]->iso)*corr_factor);
    		if (jptr->blk_T == NULL) {
    			if (ss_issame(sim->ss,jptr->nuc[0],jptr->nuc[1])) {
    				jptr->blk_T = T20(sim,jptr->nuc[0],jptr->nuc[1]);
    			} else {
    				jptr->blk_T = IzIz_sqrt2by3(sim,jptr->nuc[0],jptr->nuc[1]);
    			}
    		}
    		if (sim->Hassembly) {
    			complx *zdum = complx_vector(5);
    		    cv_zero(zdum);
    		    if (sim->J[i]->Rmol != NULL) cv_multod(zdum,sim->J[i]->Rmol,-1.0);
    		    if (jptr->Rmol != NULL) cv_multod(zdum,jptr->Rmol,1.0);
    		    blk_dm_multod(wsp->HQ[0],jptr->blk_T,zdum[3].re);
    		    blk_dm_multod(wsp->HQ[1],jptr->blk_T,zdum[4].re);
    		    blk_dm_multod(wsp->HQ[2],jptr->blk_T,zdum[4].im);
    		    blk_dm_multod(wsp->HQ[3],jptr->blk_T,zdum[5].re);
    		    blk_dm_multod(wsp->HQ[4],jptr->blk_T,zdum[5].im);
    		    free_complx_vector(zdum);
    		}
    	}
    }
    for (i=0; i<sim->nQ; i++) {
    	if (wsp->Q[i] != sim->Q[i]) {
    		Quadrupole * qptr = wsp->Q[i];
        	qptr->Rmol = Dtensor2(2.0*qptr->wq,qptr->eta);
        	qptr->Rmol = wig2roti(qptr->Rmol,qptr->pas[0],qptr->pas[1],qptr->pas[2]);
        	if (qptr->T == NULL) qptr->T =  T20II(sim,qptr->nuc);
        	if (sim->Hassembly) {
        		blk_dm_multod_diag(wsp->HQ[0],qptr->T,qptr->Rmol[3].re - sim->Q[i]->Rmol[3].re);
        		blk_dm_multod_diag(wsp->HQ[1],qptr->T,qptr->Rmol[4].re - sim->Q[i]->Rmol[4].re);
        		blk_dm_multod_diag(wsp->HQ[2],qptr->T,qptr->Rmol[4].im - sim->Q[i]->Rmol[4].im);
        		blk_dm_multod_diag(wsp->HQ[3],qptr->T,qptr->Rmol[5].re - sim->Q[i]->Rmol[5].re);
        		blk_dm_multod_diag(wsp->HQ[4],qptr->T,qptr->Rmol[5].im - sim->Q[i]->Rmol[5].im);
        	}
    	}
    }

}
