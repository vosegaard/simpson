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
        if (det[i] == 'z' || det[i] == 'a' || det[i] == 'b') continue;
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
  int i, NN, objc;

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
  f = TclGetDouble(interp,"par","field",0,-100001);
  if (f > -100000) { //i.e. it actually was defined by user
	  if ( f < TINY ) {
		  fprintf(stderr,"Error: field must be positive and nonzero\n");
		  exit(1);
	  }
	  s->Bzero = f;
	  s->specfreq = -ss_gamma1H()*f/(2*M_PI);
  } else {
	  s->Bzero = -s->specfreq*2*M_PI/ss_gamma1H();
  }

  s->wr=TclGetDouble(interp,"par","spin_rate",0,0);
  if (s->wr < 0.0) {
     fprintf(stderr,"error: spin_rate cannot be negative\n");
     exit(1);
  }
  s->sw=TclGetDouble(interp,"par","sw",0,s->wr);
  s->sw1=TclGetDouble(interp,"par","sw1",0,0);
  s->np=TclGetInt(interp,"par","np",0,1);
  s->ni=TclGetInt(interp,"par","ni",0,0);

  // read and decompose par(method)
  //  DEFAULTS:
  s->imethod = M_DIRECT_TIME;
  s->sparse = 0;
  s->block_diag = 0;
  s->propmethod = 0; // via diagonalization
  s->interpolation = INTERPOL_NOT_USED;
  s->domain = 0; // time domain simulation
  s->frame = ROTFRAME; // rotating frame is default
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
		  } else if (!strncmp(buf,"mcn",3)) {
			  s->propmethod = 7;
			  s->do_avg = 1;
		  } else if (!strncmp(buf,"sparse",6)) {
			  s->sparse = 1;
		  } else if (!strncmp(buf,"block_diag",10)) {
			  s->block_diag = 1;
		  } else if (!strncmp(buf,"FWTinterpolation",9)) {
			  s->interpolation = INTERPOL_FWT_LAM;
		  } else if (!strncmp(buf,"FWT2interpolation",10)) {
			  s->interpolation = INTERPOL_FWT_ALL;
		  } else if (!strncmp(buf,"ASGinterpolation",9)) {
			  s->interpolation = INTERPOL_ASG;
		  } else if (!strncmp(buf,"FWTASGinterpolation",12)) {
			  s->interpolation = INTERPOL_FWTASG_LAM;
		  } else if (!strncmp(buf,"FWT2ASGinterpolation",12)) {
			  s->interpolation = INTERPOL_FWTASG_ALL;
		  } else if (!strncmp(buf,"time",4)) {
			  s->domain = 0;
		  } else if (!strncmp(buf,"frequency",4)) {
			  s->domain = 1; // calculations will be in frequency domain
			  s->imethod++;
		  } else if (!strncmp(buf,"labframe",8)) {
			  s->frame = LABFRAME;
		  } else {
			  fprintf(stderr,"Error: method option '%s' not known.\n",buf);
			  fprintf(stderr,"Must be one of : direct/idirect/gcompute/igcompute,\n");
			  fprintf(stderr,"               : diag/dsyev/pade/cheby1/cheby2/taylor/mcn, \n");
			  fprintf(stderr,"               : sparse, block_diag, \n");
			  fprintf(stderr,"               : time/frequency, \n");
			  fprintf(stderr,"               : FWTinterpolation/ASGinterpolation/FWTASGinterpolation\n");
			  exit(1);

		  }
	  }  // for loop over objc
  }
  // par(method) parsed.
    
//    printf("Interpolation: %d\n", s->interpolation);
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
//    printf("Here 1\n");
  if (s->interpolation == INTERPOL_FWT_LAM || s->interpolation == INTERPOL_FWTASG_LAM || s->interpolation == INTERPOL_FWT_ALL || s->interpolation == INTERPOL_FWTASG_ALL) {
	  // load, or create and load, map of nearest crystallites target->source
//      printf("Crystal info: %s %s %ld %ld\n", s->crystfile, s->targetcrystfile, s->crdata, s->targetcrdata);
	  s->crmap = read_cryst_map(s->crystfile, s->crdata, s->targetcrystfile, s->targetcrdata);
  } else {
	  s->crmap = NULL;
  }
//    printf("Here 2\n");
  if (s->interpolation == INTERPOL_ASG || s->interpolation == INTERPOL_FWTASG_ALL || s->interpolation == INTERPOL_FWTASG_LAM) {
	  // test for triangle data
	  // at the moment they are read in mpi_ASG_interpol ...
  }

  TclGetString(interp,s->rfproffile,"par","rfprof_file",0,"none");
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

//  TclGetString(interp,startop,"par","start_operator",0,"Inz");
  obj = Tcl_GetVar2Ex(interp,"par","start_operator",TCL_GLOBAL_ONLY);
  if (obj == NULL) {  // default value is Inz
	  obj = Tcl_NewStringObj("Inz",3);
	  obj = Tcl_SetVar2Ex(interp,"par","start_operator",obj,TCL_GLOBAL_ONLY);
	  if (obj == NULL) {
		  fprintf(stderr,"Error when setting 'par(start_operator)':\n%s\n",Tcl_GetStringResult(interp));
		  exit(1);
	  }
  }
  if (Tcl_ListObjGetElements(interp, obj, &objc, &objv) != TCL_OK) {
	  fprintf(stderr,"Error: sim_initialize cannot decompose list par(start_operator):\n%s\n",Tcl_GetStringResult(interp));
	  exit(1);
  }
  if (objc < 1) {
	  fprintf(stderr,"Error: par(start_operator) must not be empty\n");
	  exit(1);
  }
  s->Nfstart = objc;
  s->fstart = (mat_complx**)malloc(s->Nfstart*sizeof(mat_complx*));
  for (i=0; i<objc; i++) {
	  strcpy(startop,Tcl_GetString(objv[i]));
	  if (!strncmp(startop,"equil",5)) {
		  s->fstart[i] = ss_eq_sigma(s);
	  } else {
		  s->fstart[i] = ss_readoper(s,startop);
	  }
	  if (s->fstart[i] == NULL) {
		  fprintf(stderr,"error: unable to interpret start operator '%s' (element %d)\n",startop,i+1);
		  exit(1);
	  }
	  if (verbose & VERBOSE_OPER) {
		  if (objc == 1) {
			  sprintf(buf,"Start operator '%s'",startop);
		  } else {
			  sprintf(buf,"Start operator element %d: '%s'",i+1,startop);
		  }
		  cm_print(s->fstart[i],buf);
	  }
  }

  //  TclGetString(interp,detectop,"par","detect_operator",0,"Inp");
    obj = Tcl_GetVar2Ex(interp,"par","detect_operator",TCL_GLOBAL_ONLY);
    if (obj == NULL) {  // default value is Inp
  	  obj = Tcl_NewStringObj("Inp",3);
  	  obj = Tcl_SetVar2Ex(interp,"par","detect_operator",obj,TCL_GLOBAL_ONLY);
  	  if (obj == NULL) {
  		  fprintf(stderr,"Error when setting 'par(detect_operator)':\n%s\n",Tcl_GetStringResult(interp));
  		  exit(1);
  	  }
    }
    if (Tcl_ListObjGetElements(interp, obj, &objc, &objv) != TCL_OK) {
    	fprintf(stderr,"Error: sim_initialize cannot decompose list par(detect_operator):\n%s\n",Tcl_GetStringResult(interp));
    	exit(1);
    }
    if (objc < 1) {
  	  fprintf(stderr,"Error: par(detect_operator) must not be empty\n");
  	  exit(1);
    }
    s->Nfdetect = objc;
    s->fdetect = (mat_complx**)malloc(s->Nfdetect*sizeof(mat_complx*));
    s->obs_nuc = (int*)malloc(s->Nfdetect*sizeof(int));
    s->conjugate_fid = (int*)malloc(s->Nfdetect*sizeof(int));
    /* conjugate_fid: This parameter is used overrule the automatic detection of the magnetogyric
      ratio for the observable nucleus. It can be set to 'auto', 'true' or 'false.
    */
    TclGetString(interp,buf,"par","conjugate_fid",0,"auto");
    if (!strcmp(buf,"auto")) {
    	NN = -1;
    } else if (!strcmp(buf,"true")) {
        NN = 1;
    } else if (!strcmp(buf,"false")) {
        NN = 0;
    } else {
        fprintf(stderr,"error: 'conjugate_fid' must be 'true', 'false' or 'auto'");
        exit(1);
    }
    for (i=0; i<s->Nfdetect; i++) s->conjugate_fid[i] = NN;
    for (i=0; i<objc; i++) {
    	strcpy(detectop,Tcl_GetString(objv[i]));
    	s->fdetect[i] = ss_readoper(s,detectop);
    	if (s->fdetect[i] == NULL) {
    		fprintf(stderr,"error: unable to interpret detect operator '%s' (element %d)\n",detectop,i+1);
    		exit(1);
    	}
    	NN = 0;
    	if (s->conjugate_fid[i] == -1) {
    		s->obs_nuc[i] = obs_nuc(s->ss,detectop,&NN);
            //printf("Observe nucleus: %d %d %d\n", s->obs_nuc[i], NN, s->ss->nchan);
    		if ( (s->obs_nuc[i] < 1) || (NN && s->ss->nchan > 1) ) {
    			fprintf(stderr,"error: in detect-operator '%s'. \nNone or more than one type of nucleus was detected. \n"
    					"Set 'conjugate_fid' to 'true' or 'false' to turn off\n"
    					"this sanity check. 'conjugate_fid' is 'auto' per default which conjugates\n"
    					"the fid if the gyromagnetic ratio is larger than zero for the observable\n"
    					"nucleus. That corrects for the standard axis convention of plotting spectra\n"
    					"that ignores the sign of gamma.\n",detectop);
    			exit(1);
    		}
    		s->conjugate_fid[i] = (ss_gamma(s->ss,s->obs_nuc[i]) > 0 ? 1 : 0);
    	} else {
    		s->obs_nuc[i] = 0;
    	}
    	if (s->frame == LABFRAME) { // obs_nuc is needed for LABFRAME
    		s->obs_nuc[i] = obs_nuc(s->ss,detectop,&NN);
    		if ( (s->obs_nuc[i] < 1) || (NN && s->ss->nchan > 1) ) {
    			fprintf(stderr,"error: in LABFRAME detect-operator '%s'. None or more than one type of nucleus was detected. \n" ,detectop);
    			exit(1);
    		}
    	}

    	//if (OCpar.isinit != 1) {
    	//	s->obs_nuc[i] = obs_nuc(s->ss,detectop,&NN);
    	//	printf("obs_nuc result from '%s': %d\n", detectop, s->obs_nuc[i]); exit(1);
    	//} // this is omitted for Optimal control, anything is allowed as detectop
    	//if ( (s->obs_nuc[i] < 1) || (NN && s->ss->nchan > 1) ) {
    	//	fprintf(stderr,"error: in detect-operator '%s'. More than one type of nucleus was detected. \n" ,detectop);
    	//	exit(1);
    	//}
    	//if (s->conjugate_fid[i] == -1) s->conjugate_fid[i]=(ss_gamma(s->ss,s->obs_nuc[i]) > 0 ? 1 : 0);

    	if (verbose & VERBOSE_OPER) {
    		if (objc == 1) {
    			sprintf(buf,"Detect operator '%s'",detectop);
    		} else {
    			sprintf(buf,"Detect operator element %d: '%s'",i+1,detectop);
    		}
    		cm_print(s->fdetect[i],buf);
    		printf("Acquisition data will%s be complex conjugated.\n",(s->conjugate_fid ? "" : " not"));
    	}
    }

    /* check numbers of start/detect operators */
    if (s->Nfstart != 1 && s->Nfdetect != 1) {
    	if (s->Nfstart != s->Nfdetect) {
    		fprintf(stderr,"Error: unsupported numbers of start/detect operators (%d/%d)\n"
    				       "       can be (1/x) or (x/1) or (x/x)\n",s->Nfstart,s->Nfdetect);
    		exit(1);
    	}
    }
    /* Total length of FID is np*ni*max(Nfstart,Nfdetect)
     * => multiple start/detect are concatenated in FID
     * but ntot is the length of individual FID
     */
    s->ntot=s->np*(s->ni > 1 ? s->ni : 1);

    // issue error when interpolation and multiple fstart/fdetect
    if (s->interpolation != INTERPOL_NOT_USED) {
  	  if (s->Nfstart != 1 || s->Nfdetect != 1) {
  		  fprintf(stderr,"Error: Interpolation can not be used when using multiple start/detect operators. Feature not implemented\n");
  		  exit(1);
  	  }
    }

  // Excitation-Detection symmetry check
  //  1 ... assume there is ED symmetry
  //  0 ... force internal check on ED symmetry
  // -1 ... no ED symmetry, no internal check
  s->EDsymmetry = TclGetInt(interp,"par","ED_symmetry",0,-1);
  if (abs(s->EDsymmetry)>1) {
	  fprintf(stderr,"Error: ED_symmetry parameter has wrong value (can be -1, 0, 1)\n");
	  exit(1);
  }
  if (s->EDsymmetry != -1) {
	  if (s->Nfstart != 1 || s->Nfdetect != 1) {
		  fprintf(stderr,"Error: ED_symmetry can not be used when using multiple start/detect operators. Feature not implemented\n");
		  exit(1);
	  }
  }
  if (s->EDsymmetry == 0) {
	  s->EDsymmetry = is_rhosymmetry(s->fstart[0],s->fdetect[0]);
  }
  //printf("ED symmetry is %d\n",s->EDsymmetry);
  if (verbose & VERBOSE_OPER) {
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

  if (s->propmethod != 7) s->do_avg = TclGetInt(interp,"par","do_avg",0,0);

  return s;
}

//#define DEBUGPRINT printf

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

  for (i=0; i<s->nG; i++) {
	  Gtensor *ptr = s->G[i];
	  if (ptr->T) free_double_matrix(ptr->T);
	  free_complx_vector(ptr->Rmol);
	  for (ii=0; ii<5; ii++) {
		  if (ptr->T2q[ii] != NULL) free_double_matrix(ptr->T2q[ii]);
	  }
	  free((char*)(ptr));
  }
  free(s->G);

  for (i=0; i<s->nHF; i++) {
	  Hyperfine *ptr = s->HF[i];
	  if (ptr->blk_Tiso) free_blk_mat_double(ptr->blk_Tiso);
	  if (ptr->blk_T) free_blk_mat_double(ptr->blk_T);
	  if (ptr->blk_Ta) free_blk_mat_double(ptr->blk_Ta);
	  if (ptr->blk_Tb) free_blk_mat_double(ptr->blk_Tb);
	  free_complx_vector(ptr->Rmol);
	  for (ii=0; ii<5; ii++) {
		  if (ptr->T2q[ii] != NULL) free_double_matrix(ptr->T2q[ii]);
	  }
	  free((char*)(ptr));
  }
  free(s->HF);

  DEBUGPRINT("SIM_DESTROY 7\n");
  for (i=0; i<s->Nfstart; i++) {
	  if (s->fstart[i]) free_complx_matrix(s->fstart[i]);
  }
  free(s->fstart);
  DEBUGPRINT("SIM_DESTROY 8\n");
  for (i=0; i<s->Nfdetect; i++) {
	  if (s->fdetect[i]) free_complx_matrix(s->fdetect[i]);
  }
  free(s->fdetect);
  free(s->obs_nuc);
  free(s->conjugate_fid);
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

Sim_wsp * wsp_initialize(Sim_info *s)
{
	Sim_wsp * wsp;
	int i, j, k, nchan, nspins, NN, matdim = s->matdim;
	mat_double *dmx;

	DEBUGPRINT("wsp_initialize start\n");
	wsp = (Sim_wsp*)malloc(sizeof(Sim_wsp));
	if (!wsp) {
		fprintf(stderr,"error: unable to allocate workspace structure\n");
		exit(1);
	}
	/* fid */
	NN = (s->Nfstart > s->Nfdetect ? s->Nfstart : s->Nfdetect);
	wsp->fid = complx_vector(s->ntot*NN);
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
	wsp->G_Rlab = (complx**)malloc(sizeof(complx*)*s->nG);
	wsp->G_Rrot = (complx**)malloc(sizeof(complx*)*s->nG);
	for (i=0; i<s->nG; i++) {
		if (s->G[i]->Rmol == NULL) {
			wsp->G_Rrot[i] = NULL;
			wsp->G_Rlab[i] = NULL;
		} else {
			wsp->G_Rrot[i] = complx_vector(5);
			wsp->G_Rlab[i] = complx_vector(5);
		}
	}
	wsp->HF_Rlab = (complx**)malloc(sizeof(complx*)*s->nHF);
	wsp->HF_Rrot = (complx**)malloc(sizeof(complx*)*s->nHF);
	for (i=0; i<s->nHF; i++) {
		if (s->HF[i]->Rmol == NULL) {
			wsp->HF_Rrot[i] = NULL;
			wsp->HF_Rlab[i] = NULL;
		} else {
			wsp->HF_Rrot[i] = complx_vector(5);
			wsp->HF_Rlab[i] = complx_vector(5);
		}
	}

	wsp->CS_used = (int*)malloc(sizeof(int)*s->nCS);
	wsp->DD_used = (int*)malloc(sizeof(int)*s->nDD);
	wsp->J_used = (int*)malloc(sizeof(int)*s->nJ);
	wsp->Q_used = (int*)malloc(sizeof(int)*s->nQ);
	wsp->G_used = (int*)malloc(sizeof(int)*s->nG);
	wsp->HF_used = (int*)malloc(sizeof(int)*s->nHF);
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
    //wsp->fstart = s->fstart;
    wsp->fstart = (mat_complx**)malloc(s->Nfstart*sizeof(mat_complx*));
    wsp->sigma  = (mat_complx**)malloc(s->Nfstart*sizeof(mat_complx*));
    for (i=0; i<s->Nfstart; i++) {
    	wsp->fstart[i] = s->fstart[i];
    	wsp->sigma[i] = NULL;
    }
    //wsp->fdetect = s->fdetect;
    wsp->fdetect = (mat_complx**)malloc(s->Nfdetect*sizeof(mat_complx*));
    for (i=0; i<s->Nfdetect; i++) wsp->fdetect[i] = s->fdetect[i];
    wsp->Nacq = s->np;
    wsp->acqphase = 0;
	wsp->evalmode = EM_NORMAL;

	/* Optimal control related */
	if (OCpar.isinit) {
		wsp->OC_propstatus = 0;
		wsp->OC_mxpos = 0;
		wsp->OC_dens = (mat_complx**)malloc(MAXOCPROPS*s->Nfstart*sizeof(mat_complx*));
		for (i=0; i<MAXOCPROPS; i++) {
			//wsp->OC_dens[i] = NULL;
			wsp->OC_props[i] = NULL;
			wsp->OC_mxcode[i] = NULL;
			for (j=0; j<MAXOCSHAPES*2; j++) {
				wsp->OC_deriv[i][j] = NULL;
			}
		}
		for (i=0; i<s->Nfstart*MAXOCPROPS; i++) {
			wsp->OC_dens[i] = NULL;
		}
		wsp->OC_phivals = complx_vector(NN);
		cv_zero(wsp->OC_phivals);
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
    wsp->G = (Gtensor **)malloc(s->nG*sizeof(Gtensor*));
    for (i=0; i<s->nG; i++) wsp->G[i] = s->G[i];
    wsp->HF = (Hyperfine **)malloc(s->nHF*sizeof(Hyperfine*));
    for (i=0; i<s->nHF; i++) wsp->HF[i] = s->HF[i];

    wsp->dw = 1.0e6/s->sw;

    wsp->FWTASG_irow = NULL;
    wsp->FWTASG_icol = NULL;

    // labframe
    if (s->frame == LABFRAME) {
    	wsp->Hlab = complx_matrix(s->matdim,s->matdim,MAT_DENSE,0,s->basis);
    	wsp->Hrflab = complx_matrix(s->matdim,s->matdim,MAT_DENSE,0,s->basis);
    } else {
    	wsp->Hlab = NULL;
    	wsp->Hrflab = NULL;
    }

    // complex Hamiltonian
    wsp->Hcplx = NULL;
    wsp->Hrf_blk = NULL;

	DEBUGPRINT("wsp_initialize end\n");

	return wsp;
}


void wsp_destroy(Sim_info *s, Sim_wsp *wsp)
{
	int i, j;

	DEBUGPRINT("wsp_destroy start\n");

    //if (wsp->fdetect != s->fdetect) free_complx_matrix(wsp->fdetect);
    //if (wsp->fstart != s->fstart) free_complx_matrix(wsp->fstart);
	for (i=0; i<s->Nfstart; i++) {
		if (wsp->fstart[i] != s->fstart[i]) free_complx_matrix(wsp->fstart[i]);
		if (wsp->sigma[i] != NULL) free_complx_matrix(wsp->sigma[i]);
	}
	free(wsp->fstart);
    free(wsp->sigma);
	for (i=0; i<s->Nfdetect; i++) {
		if (wsp->fdetect[i] != s->fdetect[i]) free_complx_matrix(wsp->fdetect[i]);
	}
	free(wsp->fdetect);
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
	for (i=0; i<s->nG; i++) {
		free_complx_vector(wsp->G_Rlab[i]);
		free_complx_vector(wsp->G_Rrot[i]);
	}
	free((char*)(wsp->G_Rlab));
	free((char*)(wsp->G_Rrot));
	free((char*)(wsp->G_used));
	for (i=0; i<s->nHF; i++) {
		free_complx_vector(wsp->HF_Rlab[i]);
		free_complx_vector(wsp->HF_Rrot[i]);
	}
	free((char*)(wsp->HF_Rlab));
	free((char*)(wsp->HF_Rrot));
	free((char*)(wsp->HF_used));

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
    		//if (wsp->OC_dens[i]) free_complx_matrix(wsp->OC_dens[i]);
    		if (wsp->OC_props[i]) free_blk_mat_complx(wsp->OC_props[i]);
    		if (wsp->OC_mxcode[i]) free(wsp->OC_mxcode[i]);
			for (j=0; j<MAXOCSHAPES*2; j++) {
				if (wsp->OC_deriv[i][j] != NULL) free_complx_matrix(wsp->OC_deriv[i][j]);
			}

    	}
    	for (i=0; i<s->Nfstart*MAXOCPROPS; i++) {
    		if (wsp->OC_dens[i]) free_complx_matrix(wsp->OC_dens[i]);
    	}
    	free(wsp->OC_dens);
    	free_complx_vector(wsp->OC_phivals);
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
    for (i=0; i<s->nG; i++) {
    	if (wsp->G[i] != s->G[i]) {
    		free_complx_vector(wsp->G[i]->Rmol);
    		if (wsp->G[i]->T != s->G[i]->T) free_double_matrix(wsp->G[i]->T);
    		free((char *)(wsp->G[i]));
    	}
    }
    free((char *)(wsp->G)); wsp->G = NULL;
    for (i=0; i<s->nHF; i++) {
    	if (wsp->HF[i] != s->HF[i]) {
    		free_complx_vector(wsp->HF[i]->Rmol);
    		if (wsp->HF[i]->blk_Tiso != s->HF[i]->blk_Tiso) free_blk_mat_double(wsp->HF[i]->blk_Tiso);
    		if (wsp->HF[i]->blk_T != s->HF[i]->blk_T) free_blk_mat_double(wsp->HF[i]->blk_T);
    		if (wsp->HF[i]->blk_Ta != s->HF[i]->blk_Ta) free_blk_mat_double(wsp->HF[i]->blk_Ta);
    		if (wsp->HF[i]->blk_Tb != s->HF[i]->blk_Tb) free_blk_mat_double(wsp->HF[i]->blk_Tb);
    		free((char *)(wsp->HF[i]));
    	}
    }
    free((char *)(wsp->HF)); wsp->HF = NULL;

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

    // complex Hamiltonian
    if (wsp->Hcplx != NULL) {
    	free_blk_mat_complx(wsp->Hcplx);
    	wsp->Hcplx = NULL;
    }
    if (wsp->Hrf_blk != NULL) {
    	free_blk_mat_complx(wsp->Hrf_blk);
    	wsp->Hrf_blk = NULL;
    }

	DEBUGPRINT("wsp_destroy end\n");

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

/* fid is with 1-based indexing!!! */
void sim_conjfid(complx *fid, Sim_info* sim)
{
	int i, N, NN;

	N = (sim->Nfstart > sim->Nfdetect ? sim->Nfstart : sim->Nfdetect);
	NN = sim->ntot;
	for (i=0; i<N; i++) {
		if (sim->conjugate_fid[i]) {
			 cv_conj2(fid+1+i*NN,NN);
		}
	}
}


/****
 * calculates pulse sequence response without any interpolations
 */
int sim_calcfid(Sim_info *s, Sim_wsp *wsp)
{
  int ig;

    //printf("Angles: %g %g %g %g\n", wsp->cryst.alpha, wsp->cryst.beta, wsp->cryst.gamma, wsp->cryst.weight);
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
      //if (s->conjugate_fid) cv_conj(wsp->fid);
      sim_conjfid(wsp->fid,s);
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
      //if (s->conjugate_fid) cv_conj(wsp->fid);
      sim_conjfid(wsp->fid,s);
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
      //if (s->conjugate_fid) cv_conj(wsp->fid);
      sim_conjfid(wsp->fid,s);
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

