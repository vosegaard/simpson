/*
    Reads B0 along z axis inhomogeneity file and interpolates values 
    used in loop over z coordinate.
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
    
    The B0z inhomogeneity file must have the following data format:

       -1.0 <delta B0z>
        <z> <delta B0z>
        <z> <delta B0z>
       ...
        1.0 <delta B0z>

    This means that z coordinate spans exactly interval <-1,1>, with end points 
    strictly defined in the file. <delta B0z> is given in Hz, proton frequency.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <tcl.h>
#include "matrix.h"
#include "defs.h"
#include "cryst.h"
#include "tclutil.h"
#include "sim.h"
#include "B0inhom.h"
#include "spinsys.h"

/****
 * ZT: generates values for z coordinate (from -1 to +1)
 ****/
double* get_zlims(int Nz) 
{
   double *zval, s;
   int i;
   
   if (Nz <=1) {
      zval=double_vector(1);
      zval[1] = 0.0;
   } else {
      zval = double_vector(Nz);
      s = 2.0/((double)Nz-1.0);
      for (i=1; i<=Nz; i++) {
         zval[i] = -1.0+((double)i-1.0)*s;
      }
   }
   return zval;
}



/****
 * ZT: reads in information for averaging over z coordinate
 ****/
void prepare_zaveraging(Tcl_Interp* interp,double** zv,double** zoffs)
{
  char profilename[256], fname[256];
  FILE* fp;
  int nz,i,ver,NZ;
  double zs, zn, Bs, Bn;
  double *zvals, *zoffsetvals;

  TclGetString(interp,profilename,"par","zprofile",0,"none"); 
  NZ=TclGetInt(interp,"par","zvals",0,1);

  zvals = get_zlims(NZ); /* zvals are allocated */
  nz = LEN(zvals);
  zoffsetvals = double_vector(nz);
  
  if (NZ < 1) {
     zoffsetvals[1] = 0.0;
     *zv = zvals;
     *zoffs = zoffsetvals;
     return;
  }
  
  /* now read in field inhomogeneity profile */
  if ( !strcmp(profilename,"none") ) {
     for (i=1; i<=nz; i++) {
        zoffsetvals[i] = 0.0;
     }
     *zv = zvals;
     *zoffs = zoffsetvals;
     return;
  }
  /* now realy read the file in */
  ver= verbose & VERBOSE_RFPROF;
  strcpy(fname,profilename);
#ifdef UNIX
  if (profilename[0] == '~') {
    char* p=getenv("HOME");
    if (p != NULL) {
      strcpy(fname,p);
      strcat(fname,&profilename[1]);
    }
  }
#endif
  if (!(fp=fopen(fname,"r"))) {
    strcat(fname,".B0z");
    fp=fopen(fname,"r");
  }
  
  if (!fp) {
    fprintf(stderr,"error: unable to open file '%s'\n",fname);
    fprintf(stderr,"\n");
    exit(1);
  } 
  
  if (ver) printf("loading external B0 field inhom. profile file '%s'\n",fname);

  if (fscanf(fp,"%lg%lg",&zn,&Bn) != 2) {
    fprintf(stderr,"unable to read first line from file '%s'\n",fname);
    exit(1);
  }
  if (NZ > 1) {
    if ( fabs(zn-zvals[1]) > TINY ) {
      fprintf(stderr,"error: Reading zprofile, first line must have z=%f, not %f\n",zvals[1],zn);
      exit(1);
    }
    zs = zn;
    Bs = Bn;
    zoffsetvals[1] = Bn;
    if (ver) printf("z = %f, offset = %f\n",zvals[1],zoffsetvals[1]);
    i = 2;
  } else {
    zs = zn;
    Bs = Bn;
    i = 1;
  } 
  while (fscanf(fp,"%lg%lg",&zn,&Bn) == 2) {
     if (zn >= (zvals[i])) {
        zoffsetvals[i] = (Bn-Bs)/(zn-zs)*(zvals[i]-zs)+Bs;
	if (ver) printf("z = %f, offset = %f\n",zvals[i],zoffsetvals[i]);
	i++;
	if (i > nz) break;
     }
     zs = zn;
     Bs = Bn;
  }
  if (i <= nz) {
     fprintf(stderr,"error reading zprofile, insufficient data in file %s (less entries than par(zvals)\n",fname);
     exit(1);
  }

  fclose(fp);
  
 if (ver) printf("Reading B0 field inhom. profile '%s' done!'\n",profilename);

 *zv = zvals;
 *zoffs = zoffsetvals;

}

/****
 * ZT: helper function calculating effects of offset on different channels
 ****/
void get_chan_offset_ratios(SpinSys* ss, double zoffnominal, double* ovals)
{
  int Nchan,i;
  Nchan = LEN(ovals);
  for (i=1;i<=Nchan;i++) {
        ovals[i] = zoffnominal/ss_gamma(ss,ss->chan[i][1])*ss_gamma1H();
  }
}  

/****
 *  ZT: sets proper offsets to wsp structure
 ****/
void set_inhom_offsets(Sim_info* s, Sim_wsp* wsp,double zoffnominal)
{
  int i, Nchan;
  double chanoffset;
  double *inh_offset;
  
  inh_offset = wsp->inhom_offset;
  Nchan = s->ss->nchan;
  
  if (zoffnominal == 0.0) {
     if (inh_offset) free_double_vector(inh_offset);
     inh_offset = NULL;
  } else {
     if (!inh_offset) inh_offset = double_vector(Nchan);
     for (i=1;i<=Nchan;i++) {
        chanoffset = zoffnominal/ss_gamma(s->ss,s->ss->chan[i][1])*ss_gamma1H();
	    inh_offset[i] = chanoffset*M_PI*2.0;
     }
  }
  wsp->inhom_offset = inh_offset;
}

/*********** z gradient stuff ********/

 /* define global array of pointers to rf shapes */
 double *ZgradShapes[MAXZGRADSHAPES];
 
 /****
  * provide free slot in ZgradShapes
  ****/
 int ZgradShapes_slot() {
   int a;
   
   for (a=0; a<MAXZGRADSHAPES; a++) {
      if (!ZgradShapes[a]) {
         break;
      }
   }
   if (a >= MAXZGRADSHAPES) {
      fprintf(stderr,"ZgradShapes error: no more free slots available, maximum is %d\n",MAXZGRADSHAPES);
      a=-1;
   }
   return a;
 }
 
 
 
/****
 * allocation function for z grad shapes 
 ****/
 double* ZgradShapes_alloc(int len) {
   double* v;
   
   
   v = (double*)malloc((len+1)*sizeof(double));
   if (!v) {
      fprintf(stderr,"error: unable to alocate ZgradShapes");
      exit(-1);
   }
   
   /* store its length to the first element */
   *(int*)v=len;
   return v;
 }
 
 
 /****
  * freeing the shape from memory 
  ****/
 void free_ZgradShapes(int a) {
 
    free((char *)ZgradShapes[a]);
    ZgradShapes[a] = NULL;
 }

/****
 * free all remaining slots in RFshapes
 ****/
 void ZgradShapes_reset() {
   int a;
   
   for (a=0; a<MAXZGRADSHAPES; a++) {
     if (ZgradShapes[a] != NULL) {
       free_ZgradShapes(a);
     }
   }
 }
 
 
 
/****
 * length of shape in a given slot
 ****/
 int ZgradShapes_len(int slot) {
 
    if (!ZgradShapes[slot]) {
      fprintf(stderr,"ZgradShapes error: testing length of empty slot\n");
      exit(1);
    }
    return *(int*)ZgradShapes[slot];
 }

/****
 * implementation of Tcl zgrad_shape_create routine 
 ****/
int tclZgradCreate(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
   int slot, i, j, len;
   char expr[2048];
   double a;
   Tcl_Obj *tclres;

  if ( (argc < 2) || (argc > 3) )
    return TclError(interp,"usage: zgrad_shape_create <num of lems>  ?<ampl expr>? ");
  if (Tcl_GetInt(interp,argv[1],&len) == TCL_ERROR)
    return TclError(interp,"zgrad_shape_create: argument 1 must be integer <num of elems>");

  /* get a new slot and allocate */
  slot = ZgradShapes_slot();
  if (slot == -1) {
     return TclError(interp,"zgrad_shape_create error: no more free slots available, free some shape first!");
  }
  ZgradShapes[slot] = ZgradShapes_alloc(len);

  for (i=1; i<=len; i++) {
     ZgradShapes[slot][i] = 0.0;
  }

  if (argc == 3) {
     /* evaluate expression in Tcl */
     for (j=1; j<=len; j++) {
        sprintf(expr,"\n set i %f\n expr %s\n", (double)j, argv[2]);
	if ( Tcl_EvalEx(interp, expr, -1,TCL_EVAL_DIRECT) != TCL_OK ) 
           return TclError(interp,"error in zgrad_shape_create: can not evaluate %s for index %d",expr, j);
        tclres = Tcl_GetObjResult(interp);
        if ( Tcl_GetDoubleFromObj(interp,tclres,&a) != TCL_OK )
           return TclError(interp,"error in zgrad_shape_create: can not get amplitude result for index %d",j);
        ZgradShapes[slot][j] = a;
     }
  }   

  return TclSetResult(interp,"%d",slot);
}

/****
 * implementation of Tcl free_zgrad routine 
 ****/
int tclZgradFree(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int slot;
  
  if ( argc != 2 )
    return TclError(interp,"usage: free_zgrad <zgrad shape>");
  if (Tcl_GetInt(interp,argv[1],&slot) == TCL_ERROR)
    return TclError(interp,"free_zgrad: argument 1 must be integer <zgrad shape>");
  
  free_ZgradShapes(slot);
  
  return TCL_OK;
}

/****
 * implementation of Tcl load_zgrad routine 
 ****/
int tclZgradLoad(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int slot;
   FILE* fp;
   char fname[256], dum[256],name[256];
   int Nelem, i;
   double am;

   if (argc != 2) 
     return TclError(interp,"usage: <zgrad shape> load_zgrad <name of file>");

   /* decide which slot to use */
   slot = ZgradShapes_slot();
   strcpy(name,argv[1]);
    
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
   fp=fopen(fname,"r");
   if (!fp) {
      fprintf(stderr,"load_zgrad error: unable to open file %s\n\n",fname);
      exit(1);
   }
   
   /* scan file for number of lines, e.g. number of elements in zgrad shape */
   Nelem = 0;
   while ( fgets(dum, 256, fp) ) {
      Nelem++;
   }
   fseek(fp, 0, SEEK_SET);
   /* printf("Number of lines = %d\n",Nelem); */
 
   ZgradShapes[slot]=ZgradShapes_alloc(Nelem);
   
   for (i=1; i<=Nelem; i++) {
      fgets(dum, 256, fp);
      if ( sscanf(dum,"%lg",&am) != 1 ) {
         fprintf(stderr,"load_zgrad error: unable to read line %d in %s\n",i,fname);
	 exit(1);
      }
      ZgradShapes[slot][i] = am;
   }
   fclose(fp);
   
   return TclSetResult(interp,"%d",slot);
}

/****
 * implementation of Tcl zgrad2list routine 
 ****/
int tclZgradList(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int slot, i;
  Tcl_Obj *lptr1;
  Tcl_Obj *elemptr;
  
  if ( argc != 2 )
    return TclError(interp,"usage: zgrad2list <zgrad shape>");
  if (Tcl_GetInt(interp,argv[1],&slot) == TCL_ERROR)
    return TclError(interp,"zgrad2list: argument 1 must be integer <zgrad shape>");
  if (!ZgradShapes[slot])
    return TclError(interp,"zgrad2list: zgrad shape doesn't seem to exist!");
    
  /* create list objects */
  lptr1 = Tcl_NewListObj(0,NULL);
  if (!lptr1) return TclError(interp,"zgrad2list unable to create outer list");

  for (i=1; i<=ZgradShapes_len(slot); i++) {
     elemptr = Tcl_NewDoubleObj(ZgradShapes[slot][i]);
     if (!elemptr) {
	return TclError(interp,"zgrad2list unable to create double from zgrad shape element %d",i);
     }
     if ( Tcl_ListObjAppendElement(interp,lptr1,elemptr) != TCL_OK ) {
        /* Tcl_Free(lptr2);
	Tcl_Free(lptr1); */
	return TclError(interp,"zgrad2list unable to append element %d to the list",i);
     }
  }
  Tcl_SetObjResult(interp,lptr1);

  return TCL_OK;
}









void tclcmd_inhom(Tcl_Interp* interp) {
  
Tcl_CreateCommand(interp,"zgrad_shape_create",(Tcl_CmdProc *)tclZgradCreate,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"free_zgrad",(Tcl_CmdProc *)tclZgradFree,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"load_zgrad",(Tcl_CmdProc *)tclZgradLoad,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateCommand(interp,"zgrad2list",(Tcl_CmdProc *)tclZgradList,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
}
