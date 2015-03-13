/*
    Routines connected with relaxation
    Copyright (C) 2008 Zdenek Tosner

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
    
    Reads the spin system from the 'spinsys' section in the input
    file and creates the Hamiltonian for the spin system.
    
     Plenty os functions called from different parts of the code :-)
*/

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include "cm.h"
#include "spinsys.h"
#include "ham.h"
#include "sim.h"
#include "pulse.h"
#include "defs.h"
#include "relax.h"

/****  Defining global variable Relax ****/
RXstruct Relax;


/****
 * Encoding spectral density of particular relaxation mechanism
 ****/
int get_Jtype(char* glm,char *jlocname)
{
   int type;
   
   if (!strncmp(glm,"non",3)) {
      type = 0;
   } else if (!strncmp(glm,"sph",3)) {
      type = 1;
   } else if (!strncmp(glm,"sym",3)) {
      type = 2;
   } else if (!strncmp(glm,"asy",3)) {
      type = 3;
   } else {
      fprintf(stderr,"get_Jtype error: unknown type of glob. motion '%s'\n",glm);
      exit(1);
   }
   if (!strncmp(jlocname,"rig",3)) {
      type += 0;
   } else if (!strncmp(jlocname,"model_free_ext",14)) {
      type += 20;
   } else if (!strncmp(jlocname,"model_free",10)) {
      type += 10;
   } else if (!strncmp(jlocname,"diffusion_on",12)) {
      type += 30;
   } else if (!strncmp(jlocname,"diffusion_in",12)) {
      type += 40;
   } else if (!strncmp(jlocname,"3_sites",7)) {
      type += 50;
   } else {
      fprintf(stderr,"get_Jtype error: unknown type of local motion '%s'\n",jlocname);
      exit(1);
   }
   if (type==0) {
      fprintf(stderr,"get_Jtype error: combination of global an local motion is not allowed\n");
      exit(1);
   }
   /* LIMIT TO LIPARI-SZABO DURING DEVELOPMENT */
   if (type != 11) {
      fprintf(stderr,"get_Jtype error: USE L-S SETTINGS DURING DEVELOPMENT! \n");
      exit(1);
   }
   return type;   
}



/****
 * Reading relax Tcl variable and initial set up of C global Relax structure
 ****/
void read_relax(Tcl_Interp *interp, Sim_info *s)
{
   int N, i, j, k, n1, n2, lam, matdim, Nj, Nglm, jpars, index;
   char **elname, *jlocname, *glmname;
   char *bufstr;
   Tcl_Obj *o1, *glm, *jlist;
   Tcl_Obj **jelptr, **glmptr;
   double om1, om2, field, cseta, dum;
   SpinSys *ss;
   SecCouple *scdum;
   mat_complx *mp1, *mp2, *mm1, *mm2, *mz1, *mz2;
   
   printf("read_relax  starts now\n");
   
   ss = s->ss;
   matdim = ss->matdim;
   field = fabs((s->specfreq)/ss_gamma1H())*2.0*M_PI; /* this is in 10^-7 Tesla */
   printf("   Magnetic field is %f\n",field*1.0e-7);
   Relax.eqnorm = 1.0/(2.0*M_PI*fabs(s->specfreq));
   printf("   Relax.eqnorm is %g\n",Relax.eqnorm);
   
   /* get type of global motion into Tcl_Obj pointer */
   if ((glm=Tcl_GetVar2Ex(interp,"relax","global_motion",TCL_GLOBAL_ONLY)) == NULL) {
      fprintf(stderr,"read_relax: error reading relax(global_motion)");
      exit(1);
   }
   if (Tcl_ListObjGetElements(interp, glm, &Nglm, &glmptr) != TCL_OK) {
      fprintf(stderr,"read_relax error: cannot decompose parameters from a list for global_motion\n");
      exit(1);
   }
   glmname = Tcl_GetString(glmptr[0]);
   if (!glmname) {
      fprintf(stderr,"read_relax error: cannot convert glm to string\n");
      exit(1);
   }
   printf("   global motion is %s\n",glmname);
   
   /* get number of relaxation mechanisms */
   if ((o1=Tcl_GetVar2Ex(interp,"relax","auto",TCL_GLOBAL_ONLY)) == NULL) {
      fprintf(stderr,"read_relax: error reading relax(auto)");
      exit(1);
   }
   if (Tcl_GetIntFromObj(interp,o1,&N) != TCL_OK) {
      fprintf(stderr,"read_relax: cannot get int from relax(auto)");
      exit(1);
   }
   Relax.Nauto = N;
   printf("   Number of relaxation mechanisms is %d\n",N);
   
   if ((o1=Tcl_GetVar2Ex(interp,"relax","cross",TCL_GLOBAL_ONLY)) == NULL) {
      fprintf(stderr,"read_relax: error reading relax(cross)");
      exit(1);
   }
   if (Tcl_GetIntFromObj(interp,o1,&N) != TCL_OK) {
      fprintf(stderr,"read_relax: cannot get int from relax(cross)");
      exit(1);
   }
   Relax.Ncross = N;
   printf("   Number of cross correlations is %d\n",N);

   /* allocate rows in Q array */
   Relax.Q = (mat_complx***)malloc(Relax.Nauto*sizeof(mat_complx**));
   if (!Relax.Q) {
      fprintf(stderr,"read_relax error: cannot allocate rows in Q\n");
      exit(1);
   }
   /* allocate rows in Y array */
   Relax.Y = (mat_complx***)malloc(Relax.Nauto*sizeof(mat_complx**));
   if (!Relax.Y) {
      fprintf(stderr,"read_relax error: cannot allocate rows in Y\n");
      exit(1);
   }
   /* allocate rows in NQ vector */
   Relax.NQ = (int*)malloc(Relax.Nauto*sizeof(int));
   if (!Relax.NQ) {
      fprintf(stderr,"read_relax error: cannot allocate rows in NQ\n");
      exit(1);
   }
   /* allocate rows in omega array */
   Relax.omega = (double**)malloc(Relax.Nauto*sizeof(double*));
   if (!Relax.omega) {
      fprintf(stderr,"read_relax error: cannot allocate rows in omega\n");
      exit(1);
   }
   /* allocate rows in Jcode array */
   Relax.Jcode = (SpecDens*)malloc(Relax.Nauto*sizeof(SpecDens));
   if (!Relax.Jcode) {
      fprintf(stderr,"read_relax error: cannot allocate rows in Jcode\n");
      exit(1);
   }
   /* allocate rows in secular array */
   Relax.secular = (SecCouple**)malloc(Relax.Nauto*sizeof(SecCouple*));
   if (!Relax.secular) {
      fprintf(stderr,"read_relax error: cannot allocate rows in secular\n");
      exit(1);
   }

   /* get elements in relax array */
   if (Tcl_Eval(interp,"array names relax") != TCL_OK) {
      fprintf(stderr,"read_relax: %s\n",Tcl_GetStringResult(interp));
      exit(1);
   }
   bufstr = Tcl_GetStringResult(interp);

   if (Tcl_SplitList(interp, bufstr, &N, &elname) != TCL_OK) {
      fprintf(stderr,"relax_read error: cannot decompose array names from a list\n");
      exit(1);
   }

   /* go through the elements and create all necessary info */
   lam = -1;
   for (i=0; i<N; i++) {
      printf("   inspecting element %s\n",elname[i]);
      
      if (!strncmp(elname[i],"auto",4)) continue;
      if (!strncmp(elname[i],"cross",5)) continue;
      if (!strncmp(elname[i],"global_motion",13)) continue;
      
      /* dipole dipole relaxation */
      if (!strncmp(elname[i],"dipole",6)) {
	 if (sscanf(elname[i],"dipole_%d_%d",&n1,&n2) != 2) {
	    fprintf(stderr,"read_relax error: cannot get nuclei numbers from %s\n",elname[i]);
	    exit(1);
	 }
         printf("read_relax found dipole between nuclei %d and %d\n", n1, n2);
	 if (n1>n2) {
	    int dum;
	    dum = n2;
	    n2 = n1;
	    n1 = dum;
	 }
	 index = dipole_exist(s,n1,n2);
	 if ( index < 0 ) {
	    fprintf(stderr,"read_relax error: dipole %d %d not defined in spinsys\n",n1, n2);
	    exit(1);
	 }
	 lam++;
	 printf("   it has index lam = %d\n",lam);
	 /* allocate 9 matrix representations of Q, fill them */
	 Relax.Q[lam] = (mat_complx**)malloc(9*sizeof(mat_complx*));
	 mp1 = Ip(s,n1);
	 mp2 = Ip(s,n2);
	 mm1 = Im(s,n1);
	 mm2 = Im(s,n2);
	 mz1 = Iz(s,n1);
	 mz2 = Iz(s,n2);
	 Relax.Q[lam][0] = cm_mul(mm1,mm2); cm_muld(Relax.Q[lam][0],0.5);
	 Relax.Q[lam][1] = cm_mul(mm1,mz2); cm_muld(Relax.Q[lam][1],0.5);
	 Relax.Q[lam][2] = cm_mul(mz1,mm2); cm_muld(Relax.Q[lam][2],0.5);
	 Relax.Q[lam][3] = cm_mul(mz1,mz2); cm_muld(Relax.Q[lam][3],2.0/sqrt(6.0));
	 Relax.Q[lam][4] = cm_mul(mp1,mm2); cm_muld(Relax.Q[lam][4],-0.5/sqrt(6.0));
	 Relax.Q[lam][5] = cm_mul(mm1,mp2); cm_muld(Relax.Q[lam][5],-0.5/sqrt(6.0));
	 Relax.Q[lam][6] = cm_mul(mp1,mz2); cm_muld(Relax.Q[lam][6],-0.5);
	 Relax.Q[lam][7] = cm_mul(mz1,mp2); cm_muld(Relax.Q[lam][7],-0.5);
	 Relax.Q[lam][8] = cm_mul(mp1,mp2); cm_muld(Relax.Q[lam][8],0.5);
	 free_complx_matrix(mp1);
	 free_complx_matrix(mp2);
	 free_complx_matrix(mm1);
	 free_complx_matrix(mm2);
	 free_complx_matrix(mz1);
	 free_complx_matrix(mz2);
	 printf("   Relax.Q[%d][0-8] filled.\n",lam);
	 /* allocate space for 9 Y matrices */
	 Relax.Y[lam] = (mat_complx**)malloc(9*sizeof(mat_complx*));
	 for (j=0; j<9; j++) {
            Relax.Y[lam][j] = NULL;
         }
	 /* set number of allocated Q matrices */
	 Relax.NQ[lam] = 9;
         /* allocate and fill frequencies */
	 Relax.omega[lam] = (double*)malloc(9*sizeof(double));
	 om1 = ss_gamma(s->ss,n1)*field; /* this is now in rad.s-1 */
	 om2 = ss_gamma(s->ss,n2)*field;
	 Relax.omega[lam][0] = -om1-om2;
	 Relax.omega[lam][1] = -om1;
	 Relax.omega[lam][2] = -om2;
	 Relax.omega[lam][3] = 0.0;
	 Relax.omega[lam][4] = om1-om2;
	 Relax.omega[lam][5] = -om1+om2;
	 Relax.omega[lam][6] = om1;
	 Relax.omega[lam][7] = om2;
	 Relax.omega[lam][8] = om1+om2;
	 printf("   Relax.omega[%d][0-8] filled.\n",lam);
	 /* allocate and fill secular couples */
	 if (ss_issame(ss,n1,n2)) {
	    /* homonuclear */
	    scdum = (SecCouple*)malloc(20*sizeof(SecCouple));
	    *(int*)(scdum) = 19;
	    for (k=1; k<=9; k++) {
	       scdum[k].q = k-1;
	       scdum[k].y = k-1;
	    }
            scdum[10].q = 1; scdum[10].y = 2;
            scdum[11].q = 2; scdum[11].y = 1;
            scdum[12].q = 3; scdum[12].y = 4;
            scdum[13].q = 4; scdum[13].y = 3;
            scdum[14].q = 3; scdum[14].y = 5;
            scdum[15].q = 5; scdum[15].y = 3;
            scdum[16].q = 4; scdum[16].y = 5;
            scdum[17].q = 5; scdum[17].y = 4;
            scdum[18].q = 6; scdum[18].y = 7;
            scdum[19].q = 7; scdum[19].y = 6;
	    printf("   homonuclear scdum[1-19] filled");
	 } else {
	    /* heteronuclear */
            scdum = (SecCouple*)malloc(10*sizeof(SecCouple));
	    *(int*)scdum = 9;
	    for (k=1; k<=9; k++) {
	       scdum[k].q = k-1;
	       scdum[k].y = k-1;
	    }
	    printf("   heteronuclear scdum[1-9] filled");
	 }
	 Relax.secular[lam] = scdum;
	 printf(" and finaly assigned to Relax.secular[%d].\n",lam);
	 printf(" Cross check: %d is %d\n",*(int*)scdum,*(int*)(Relax.secular[lam]));
	 /* decide on spectral density */
	 if ((jlist=Tcl_GetVar2Ex(interp,"relax",elname[i],TCL_GLOBAL_ONLY)) == NULL) {
            fprintf(stderr,"read_relax: error reading relax(%s)",elname[i]);
            exit(1);
         }
	 if (Tcl_ListObjGetElements(interp, jlist, &Nj, &jelptr) != TCL_OK) {
            fprintf(stderr,"read_relax error: cannot decompose parameters from a list for %s\n",elname[i]);
            exit(1);
         }
         jlocname = Tcl_GetString(jelptr[0]);
         if (!jlocname) {
            fprintf(stderr,"relax_read error: cannot convert local motion to string for %s\n", elname[i]);
            exit(1);
         }
 	 Relax.Jcode[lam].type = get_Jtype(glmname,jlocname);
     printf("   Relax.Jcode[%d] asigned to %d",lam,Relax.Jcode[lam].type);
	 jpars = Nj+Nglm-2+1;
	 printf(", pars will have %d elems:\n",jpars);
	 Relax.Jcode[lam].par = double_vector(jpars);
	 /* set interaction constant */
	 dum = s->DD[index]->delta*2.0*M_PI;
	 Relax.Jcode[lam].par[1] = 0.5*dum*dum*6.0;
     printf("   [%g;",Relax.Jcode[lam].par[1]);
	 /* set other parameters */
	 j = 2;
	 for (k=1; k<Nglm; k++) {
	    if (Tcl_GetDoubleFromObj(interp,glmptr[k],&(Relax.Jcode[lam].par[j])) != TCL_OK) {
	       fprintf(stderr,"read_relax error: (glm) cannot get J parameter as double (%s)\n",elname[i]);
	       exit(1);
	    }
	    printf(" %g;",Relax.Jcode[lam].par[j]);
	    j++;
	 }
	 for (k=1; k<Nj; k++) {
	    if (Tcl_GetDoubleFromObj(interp,jelptr[k],&(Relax.Jcode[lam].par[j])) != TCL_OK) {
	       fprintf(stderr,"read_relax error: (loc) cannot get J parameter as double (%s)\n",elname[i]);
	       exit(1);
	    }
	    printf(" %g;",Relax.Jcode[lam].par[j]);
	    j++;
	 }
	 printf("\n");
	 Tcl_Free((char *) jelptr);
	 continue;
      }
      
      /* chemical shift anisotropy relaxation */
      if (!strncmp(elname[i],"shift",5)) {
	 if (sscanf(elname[i],"shift_%d",&n1) != 1) {
	    fprintf(stderr,"read_relax error: cannot get nucleus number from %s\n",elname[i]);
	    exit(1);
	 }
     printf("read_relax found shift for nucleus %d\n",n1);
     index = shift_exist(s,n1);
	 if ( index < 0 ) {
	    fprintf(stderr,"read_relax error: shift %d not defined in spinsys\n",n1);
	    exit(1);
	 }
	 lam++;
	 printf("   it has index lam = %d\n",lam);
	 /* allocate 3 matrix representations of Q, fill them */
	 Relax.Q[lam] = (mat_complx**)malloc(3*sizeof(mat_complx*));
	 Relax.Q[lam][0] = Im(s,n1); cm_muld(Relax.Q[lam][0],0.5);
	 Relax.Q[lam][1] = Iz(s,n1); cm_muld(Relax.Q[lam][1],2.0/sqrt(6.0));
	 Relax.Q[lam][2] = Ip(s,n1); cm_muld(Relax.Q[lam][2],-0.5);
	 printf("   Relax.Q[%d][0-2] filled.\n",lam);
	 /* allocate space for 3 Y matrices */
	 Relax.Y[lam] = (mat_complx**)malloc(3*sizeof(mat_complx*));
	 for (j=0; j<3; j++) {
            Relax.Y[lam][j] = NULL;
     }
	 /* set number of allocated Q matrices */
	 Relax.NQ[lam] = 3;
         /* allocate and fill frequencies */
	 Relax.omega[lam] = (double*)malloc(3*sizeof(double));
	 om1 = ss_gamma(s->ss,n1)*field; /* this is in rad.s-1 */
	 Relax.omega[lam][0] = -om1;
	 Relax.omega[lam][1] = 0.0;
	 Relax.omega[lam][2] = om1;
	 printf("   Relax.omega[%d][0-2] filled.\n",lam);
	 /* allocate and fill secular couples */
	 Relax.secular[lam] = (SecCouple*)malloc(4*sizeof(SecCouple));
	 *(int*)(Relax.secular[lam]) = 3;
	 Relax.secular[lam][1].q = 0;
	 Relax.secular[lam][1].y = 0;
	 Relax.secular[lam][2].q = 1;
	 Relax.secular[lam][2].y = 1;
	 Relax.secular[lam][3].q = 2;
	 Relax.secular[lam][3].y = 2;
	 printf("   Relax.secular[%d][1-3] filled\n",lam);
	 /* decide on spectral density */
	 if ((jlist=Tcl_GetVar2Ex(interp,"relax",elname[i],TCL_GLOBAL_ONLY)) == NULL) {
            fprintf(stderr,"read_relax: error reading relax(%s)",elname[i]);
            exit(1);
         }
	 if (Tcl_ListObjGetElements(interp, jlist, &Nj, &jelptr) != TCL_OK) {
            fprintf(stderr,"read_relax error: cannot decompose parameters from a list for %s\n",elname[i]);
            exit(1);
         }
         jlocname = Tcl_GetString(jelptr[0]);
         if (!jlocname) {
            fprintf(stderr,"get_Jtype error: cannot convert local motion to string for %s\n", elname[i]);
            exit(1);
         }
 	 Relax.Jcode[lam].type = get_Jtype(glmname,jlocname);
     printf("   Relax.Jcode[%d] asigned to %d",lam,Relax.Jcode[lam].type);
	 jpars = Nj+Nglm-2+1;
	 printf(", pars will have %d elems:\n",jpars);
	 Relax.Jcode[lam].par = double_vector(jpars);
	 /* set interaction constant */
	 dum = s->CS[index]->delta*2.0*M_PI;
	 cseta = s->CS[index]->eta;
	 Relax.Jcode[lam].par[1] = 0.5*dum*dum*(1.0+cseta*cseta/3.0);
     printf("   [%g;",Relax.Jcode[lam].par[1]);
	 /* set other parameters */
	 j = 2;
	 for (k=1; k<Nglm; k++) {
	    if (Tcl_GetDoubleFromObj(interp,glmptr[k],&(Relax.Jcode[lam].par[j])) != TCL_OK) {
	       fprintf(stderr,"read_relax error: (glm) cannot get J parameter as double (%s)\n",elname[i]);
	       exit(1);
	    }
	    printf(" %g;",Relax.Jcode[lam].par[j]);
	    j++;
	 }
	 for (k=1; k<Nj; k++) {
	    if (Tcl_GetDoubleFromObj(interp,jelptr[k],&(Relax.Jcode[lam].par[j])) != TCL_OK) {
	       fprintf(stderr,"read_relax error: (loc) cannot get J parameter as double (%s)\n",elname[i]);
	       exit(1);
	    }
	    printf(" %g;",Relax.Jcode[lam].par[j]);
	    j++;
	 }
	 printf("\n");
	 Tcl_Free((char *) jelptr);
	 continue;
      }
      
      /* quadrupolar relaxation */
      if (!strncmp(elname[i],"quadrupole",10)) {
	 if (sscanf(elname[i],"quadrupole_%d",&n1) != 1) {
	    fprintf(stderr,"read_relax error: cannot get nucleus number from %s\n",elname[i]);
	    exit(1);
	 }
         printf("read_relax found quadrupole for nucleus %d\n", n1);
	 /* NOT IMPLEMENTED FURTHER !!! */
	 fprintf(stderr,"read_relax error: quadrupole not implemented yet\n");
	 exit(1);
	 lam++;
      }
      
      /* random field relaxation */
      if (!strncmp(elname[i],"random_field",12)) {
	 if (sscanf(elname[i],"random_field_%d",&n1) != 1) {
	    fprintf(stderr,"read_relax error: cannot get nucleus number from %s\n",elname[i]);
	    exit(1);
	 }
         printf("read_relax found random_field for nucleus %d\n", n1);
	 /* NOT IMPLEMENTED FURTHER !!! */
	 fprintf(stderr,"read_relax error: random field not implemented yet\n");
	 exit(1);
	 lam++;
      }
      
      /* cross correlations not implemented yet */
      if (!strncmp(elname[i],"cross_correlation",6)) {
	 /* NOT IMPLEMENTED FURTHER !!! */
         printf("read_relax found cross_correlation\n");
	 fprintf(stderr,"read_relax: cross correlations not implemented yet\n");
	 exit(1);
      }

   }
   printf("read_relax is DONE.\n");
}


void destroy_Relax(int Nnuc)
{
   int i, j;
   int N, NN, qmx, ymx;
   
   printf("Destroy Relax\n=============\n");
   printf("Nauto = %d, Ncross = %d\n",Relax.Nauto, Relax.Ncross);
   if (!(Relax.secular)) {
     printf("  problem, Relax.secular does not exist\n");
     exit(1);
   } else {
     printf("  OK 1\n");
   }
      
   for (i=0; i<Relax.Nauto; i++) {
      if (!(Relax.secular[i])) {
         printf("  problem, Relax.secular[%d] does not exist\n",i);
	 exit(1);
      } else {
         printf("  OK 2\n");
      }
      N = *(int*)(Relax.secular[i]);
      printf("  OK 3, N = %d\n",N);
      qmx = 0;
      ymx = 0;
      for (j=1; j<=N; j++) {
         if (Relax.secular[i][j].q > qmx) qmx = Relax.secular[i][j].q;
	 if (Relax.secular[i][j].y > ymx) ymx = Relax.secular[i][j].y;
	 /* printf("  OKej %d\n",j); */
      }

      printf("Mech. %d : %d secular terms, Qmax = %d, Ymax = %d\n",i,N,qmx,ymx);
      printf("           free Q[%d]:",i);
      
      for (j=0; j<=qmx; j++) {
         if (Relax.Q[i][j]) {
            free_complx_matrix(Relax.Q[i][j]);
            printf("(%d,%d)",i,j);
         }
      }
      free((char *)(Relax.Q[i]));
      printf(" done.\n");
      
      if (Relax.Y[i]) {
         printf("           free Y[%d]:",i);
         for (j=0; j<=ymx; j++) {
            if (Relax.Y[i][j]) {
               free_complx_matrix(Relax.Y[i][j]);
               printf(" (%d,%d)",i,j);
            }
         }
         free((char *)(Relax.Y[i]));
         printf(" done.\n");
      } else {
         printf("           Y[%d] not present\n",i);
      }
      
      printf("           Jtype = %d, par = ",Relax.Jcode[i].type);
      NN = LEN(Relax.Jcode[i].par);
      for (j=1; j<=NN; j++) {
         printf(" %g;",Relax.Jcode[i].par[j]);
      }
      printf("\n           freeing Jcode[%d].par",i);
      free_double_vector(Relax.Jcode[i].par);
      printf(" done.\n");
      
      printf("           secular =");
      for (j=1; j<=N; j++) {
         printf(" (%d,%d)",Relax.secular[i][j].q, Relax.secular[i][j].y);
      }
      printf("\n           freeing secular[%d]",i);
      free((char *)(Relax.secular[i]));
      printf(" done.\n");

      printf("           omega =");
      for (j=0; j<=qmx; j++) {
         printf(" %g;",Relax.omega[i][j]);
      }
      printf("\n           freeing omega[%d]",i);
      free((char *)(Relax.omega[i]));
      printf(" done.\n");
      
   }
   printf("Freeing Jcode");
   free((char *)(Relax.Jcode));
   printf(" done.\n");
   
   printf("Freeing secular");
   free((char *)(Relax.secular));
   printf(" done.\n");
   
   printf("Freeing omega");
   free((char *)(Relax.omega));
   printf(" done.\n");
   
   printf("Freeing Q");
   free((char *)(Relax.Q));
   printf(" done.\n");
   
   if (Relax.Y)
   printf("Freeing omega[%d]",i);
   free((char *)(Relax.Y));
   printf(" done.\n");
   
}

/****
 * Spectral density functions, based on key Relax.Jcode
 ****/
double spectral_density(int lam, double om)
{
   double res, kk;
   int N, type;
   
   if (!(Relax.Jcode)) {
      fprintf(stderr,"spectral_density error: Relax.Jcode does not exist.\n");
      exit(1);
   } 
   if (!(Relax.Jcode[lam].par)) {
      fprintf(stderr,"spectral_density error: Relax.Jcode[%d].par does not exist.\n",lam);
      exit(1);
   }
   
   type = Relax.Jcode[lam].type;
   N = LEN(Relax.Jcode[lam].par);
   
   switch (type)
     {
      case 11:
       { /* Lipari and Szabo */
	 double tauM, taue, SS;
	 if (N != 4) {
	    fprintf(stderr,"spectral_density error: not enough parameters in Relax.Jcode[%d].par for type %d (%d, not 4)\n",lam,type,N);
	    exit(1);
	 }
	 kk = Relax.Jcode[lam].par[1];
	 tauM = 1.0/(6.0*Relax.Jcode[lam].par[2]);
	 taue = 1.0/tauM + 1.0/(Relax.Jcode[lam].par[4]);
	 taue = 1.0/taue;
	 SS = Relax.Jcode[lam].par[3];
	 res = (SS*tauM)/(1.0+om*om*tauM*tauM)+((1.0-SS)*taue)/(1+om*om*taue*taue);
	 res = 2.0/5.0*res*kk;
         break;
	}
      default:
         {
	  fprintf(stderr,"spectral_density error: unknown type %d\n",type);
	  exit(1);
	 }
     }
     
   return res;
}


/****
 * Here we go with relaxation evolution...
 ****/
void _delay_relax(Sim_info *sim, Sim_wsp *wsp, double duration)
{
  /*
  int i,j,k,l,n,q,y, Ndim, Nsec, isHdiag;
  double dt, omq, om, omsc, J;
  mat_complx *Q, *Y, *mx, *mx1, *mx2;
  mat_double *Ham;
  */

  if (sim->wr != 0) {
    fprintf(stderr,"error: relaxation not supported for sample rotation\n");
    exit(1);
  }

  
  fprintf(stderr,"Error: relaxation not implemented...\n");
  exit(1);

}

