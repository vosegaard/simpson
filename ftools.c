/*
    Tcl data manipulation routines
    Copyright (C) 1999 Mads Bak, 2001,2002 Thomas Vosegaard
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
    
    
    Various routines that manipulate the SIMPSON data format.
    Read about them in the manual. They operate on data sets
    that are numbers referring to a structure containing the data.
    ZT: changing duplicate type double2 for complx
*/

#ifdef __APPLE__
#  include <malloc/malloc.h>
#else
#  include <malloc.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <tcl.h>
#include "iodata.h"
#include "tclutil.h"
#include <unistd.h>
#include "defs.h"

#define UNITS_HZ	(1<<1)
#define UNITS_PPM	(1<<2)
#define UNITS_KHZ	(1<<3)
#define UNITS_S		(1<<4)
#define UNITS_MS	(1<<5)
#define UNITS_US	(1<<6)



int TclGetDoubleWithUnits(Tcl_Interp* interp,char* value,
                     FD* f, double* returnvalue, double scale, int units)
{
  double val,factor;
  char ivalue[128], *p;
  
  strcpy(ivalue, value);
  ivalue[strlen(value)] = 0;
  p = ivalue;
  while (!isalpha(*p) && *p) p++;
  if (*p == 0) {
    factor = 1;
  } else {
    if (!strcasecmp(p, "hz") && units & UNITS_HZ) {
      factor = scale*f->np/f->sw;
    } else if (!strcasecmp(p, "khz") && units & UNITS_KHZ) {
      factor = scale*1000*f->np/f->sw;
    } else if (!strcasecmp(p, "ppm") && units & UNITS_PPM) {
      factor = scale*f->sfrq*f->np/f->sw;
    } else if (!strcasecmp(p, "s") && units & UNITS_S) {
      factor = f->sw*scale;
    } else if (!strcasecmp(p, "ms") && units & UNITS_MS) {
      factor = f->sw*scale/1.0e3;
    } else if (!strcasecmp(p, "us") && units & UNITS_US) {
      factor = f->sw*scale/1.0e6;
    } else {
      return TCL_ERROR;
    }
    *p = 0;
  }
  if (Tcl_GetDouble(interp,ivalue,&val) == TCL_ERROR) {
    fprintf(stderr, "unable to convert '%s' to a double value",ivalue);
    exit(-1);
  }
  val *= factor;
  *returnvalue = val;
  return TCL_OK;
}


/* vec and vec2 can be the same data. */
void dscalephase(complx* vec,complx* vec2,int nvec,double scale,double offset,
                 double rp,double lp)
{
  double w,phase;
  double cosph,sinph,re,im;
  int i;

  for (i=1;i<=nvec;i++) {
    w = (double)(i-1)/(double)nvec;
    phase= ( rp + lp * w ) * DEG2RAD;
    cosph=cos(phase);
    sinph=sin(phase);
    re=vec[i].re;
    im=vec[i].im;
    vec2[i].re=scale*(re*cosph-im*sinph)+offset;
    vec2[i].im=scale*(re*sinph+im*cosph)+offset;
  }
}

/* vec and vec2 can be the same data. */
void dphase(complx* vec,complx* vec2,int nvec,double rp,double lp)
{
  double w,phase;
  double cosph,sinph,re,im;
  int i;

  for (i=1;i<=nvec;i++) {
    w = (double)(i-1)/(double)nvec;
    phase= ( rp + lp * w ) * DEG2RAD;
    cosph=cos(phase);
    sinph=sin(phase);
    re=vec[i].re;
    im=vec[i].im;
    vec2[i].re=(re*cosph-im*sinph);
    vec2[i].im=(re*sinph+im*cosph);
  }
}

/* vec and vec2 can be the same data. */
void dscale(complx* vec,complx* vec2,int nvec,double scale,double offset)
{
  int i;

  for (i=1;i<=nvec;i++) {
    vec2[i].re=scale*vec[i].re+offset;
    vec2[i].im=scale*vec[i].im+offset;
  }
}

void daddlb(complx* vec,int nvec,double sw,double lineb,double
  glfrac,int sym,double top)
{
  double tmp,lb,gb,dt,fac,fag,x;
  int i;

  lb=lineb*(1.0-glfrac); /* lorentz line broadening */
  gb=lineb*glfrac;       /* gauss line broadening */
  dt=1.0/sw;
  fac= -M_PI*lb*dt;
  fag= 2.0*M_PI*gb*dt/(4.0*sqrt(log(2.0)));
  if (sym == 1) {
    int np2 = nvec/2, j;
    for (i=0,j=nvec;i<np2;i++,j--) {
      x=fag*i;
      tmp=exp(fac*i-x*x);
      vec[i+1].re *= tmp;
      vec[i+1].im *= tmp;
      vec[j].re *= tmp;
      vec[j].im *= tmp;
    }  
  } else {
    for (i=0;i<(int)top;i++) {
      x=fag*(top-i);
      tmp=exp(fac*(top-i)-x*x);
      vec[i+1].re *= tmp;
      vec[i+1].im *= tmp;
    }  
    for (i=(int)top;i<nvec;i++) {
      x=fag*(i-top);
      tmp=exp(fac*(i-top)-x*x);
      vec[i+1].re *= tmp;
      vec[i+1].im *= tmp;
    }  
  }
}

void daddlb2d(complx* vec,int np,int ni,double sw,double lineb,double glfrac,
              double sw1,double lineb1,double glfrac1,
	      int phsens,double _top)
{
  double tmp,lb,gb,dt,fac,fag,fac1,fag1,f,x,top;
  int i,j,k,l;

  lb=lineb*(1.0-glfrac); /* lorentz line broadening */
  gb=lineb*glfrac;       /* gauss line broadening */
  dt=1.0/sw;
  fac= -M_PI*lb*dt;
  fag= 2.0*M_PI*gb*dt/(4.0*sqrt(log(2.0)));

  lb=lineb1*(1.0-glfrac1); /* lorentz line broadening */
  gb=lineb1*glfrac1;       /* gauss line broadening */
  dt=1.0/sw1;
  fac1= -M_PI*lb*dt;
  fag1= 2.0*M_PI*gb*dt/(4.0*sqrt(log(2.0)));

  top=_top;
  if (top < 1) top = 0;
  if (top > np) top = np;

  if (top == 0) {
    for (j=0;j<ni;j++) {
      if (phsens == 0)
        f=fac1*j-fag1*fag1*j*j;
      else {
        l=j/2;
        f=fac1*l-fag1*fag1*l*l;
      }
      k=np*j+1;
      for (i=0;i<np;i++) {
        x=fag*i;
        tmp=exp(fac*i-x*x+f);
        vec[k+i].re *= tmp;
        vec[k+i].im *= tmp;
      }  
    }
  } else {
    for (j=0;j<ni;j++) {
      if (phsens == 0)
        f=fac1*j-fag1*fag1*j*j;
      else {
        l=j/2;
        f=fac1*l-fag1*fag1*l*l;
      }
      k=np*j+1;
      for (i=0;i<(int)top;i++) {
        x=fag*(top-i);
        tmp=exp(fac*i-x*x+f);
        vec[k+i].re *= tmp;
        vec[k+i].im *= tmp;
      }  
      for (i=(int)top;i<np;i++) {
        x=fag*(i-top);
        tmp=exp(fac*i-x*x+f);
        vec[k+i].re *= tmp;
        vec[k+i].im *= tmp;
      }  
    }
  }
}

char ferrormsg[256];

#define F_ALLOC_INIT 3

int nfd=0;
int nallocfd=0;
FD** fd;

int fnew(FD* f)
{
  FD **fdold;
  int i,empty;

  if (!nallocfd || (nallocfd < nfd + 1)) {
    fdold=NULL;
    if (!nallocfd) {
      nallocfd = F_ALLOC_INIT;
    } else {
      fdold=fd;
      nallocfd *= 2;
    }
    fd=(FD**)malloc(sizeof(FD*)*(nallocfd+1));
    if (!fd) {
      sprintf(ferrormsg,"fload: unable to allocate %d FD structures\n",nallocfd+1);
      return -1;
    }
    if (fdold != NULL) {
      for (i=1;i<=nfd;i++) {
        fd[i]=fdold[i];
      }
      free((char*)fdold);    
    }
    for (i=nfd+1;i<=nallocfd;i++) {
      fd[i]=NULL;
    }
  }
  empty=0;
  for (i=1;i<=nfd;i++) {
    if (fd[i] == NULL) {
      empty=i;
      break;
    }
  }
  if (empty == 0) {
    empty = ++nfd;
  }
  fd[empty]=f;
  return empty;
}

int fload(char *fname)
{
  FD* f;
  f=FD_read(fname);
  return fnew(f); 
}

int floaddata(char *data)
{
  FD* f;
  f=FD_readstr(data);
  return fnew(f); 
}

int floadnmrpipe(char *fname)
{
  FD* f;
  f=FD_read_nmrpipe(fname);
  return fnew(f); 
}

int funload(int fidN)
{
  int i;

  if (fidN < 0 || fidN > nfd) {
    fprintf(stderr,"funload: argument out of range"
                   " (%d not inside 1 to %d)\n",fidN,nfd);
    return -1;
  }

  if (fidN == 0) {
     for (i=1;i<=nallocfd;i++) {
       if (fd[i] != NULL) {
          FD_free(fd[i]);
       }
       fd[i]=NULL;
     }
     free((char*)fd);
     nfd=0;
     nallocfd=0;
  } else {
    if (fd[fidN] != NULL) {
      FD_free(fd[fidN]);
      fd[fidN]=NULL;
    }
  }
  return 0;
}


int tclFLoad(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN;
  char buf1[16], buf2[128];

  if (argc != 2 && argc != 3)
     return TclError(interp,"Usage: fload <filename> ?-nmrpipe?");

  if (argc == 3) {
    if (strcmp(argv[2], "-nmrpipe"))
      return TclError(interp,"Error: fload: unknown option '%s'", argv[2]);
    fidN=floadnmrpipe(argv[1]);
  } else {
    fidN=fload(argv[1]);
  }
  if (fidN == -1)
     return TclError(interp,ferrormsg);     
  sprintf(buf1, "%d", fidN);
  sprintf(buf2, "%ld", (long)fd[fidN]);
  Tcl_Eval(interp,"namespace eval FD {variable f}");
  Tcl_SetVar2(interp,"FD::f",buf1,buf2,TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG);
  Tcl_SetVar2(interp,"FD_Internal",buf1,buf2,TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG);
  TclSetResult(interp,"%d",fidN);
  
  return TCL_OK;
}  

int tclFLoaddata(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN;
  char buf1[16], buf2[128];

  if (argc != 2)
    return TclError(interp,"Usage: floaddata <data>");

  fidN=floaddata(argv[1]);
  if (fidN == -1)
     return TclError(interp,ferrormsg);     
  sprintf(buf1, "%d", fidN);
  sprintf(buf2, "%ld", (long)fd[fidN]);
  Tcl_Eval(interp,"namespace eval FD {variable f}");
  Tcl_SetVar2(interp,"FD::f",buf1,buf2,TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG);
  Tcl_SetVar2(interp,"FD_Internal",buf1,buf2,TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG);
  TclSetResult(interp,"%d",fidN);
  
  return TCL_OK;
}  


int tclFUnload(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN;
  char buf1[16];
  
  if (argc != 1 && argc != 2)
     return TclError(interp,"Usage: funload ?<data set>?");

  fidN=0;
  if (argc == 2) {
    if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
      return TclError(interp,"funload: argument 1 must be integer <data set>");
  }     
  if (funload(fidN) == -1)
    return TclError(interp,"%s",ferrormsg);
  sprintf(buf1,"%d",fidN);
  Tcl_UnsetVar2(interp,"FD_Internal::f",buf1,0);
  return TCL_OK;
}  


int fsave(FD* f,char* fname,int format,int prec)
{
  FD_write(fname,f,format,prec);
  return 0;
}

int fsave_xreim(FD* f,char* fname)
{
  FILE* fp;
  int i,np;
  complx* data;
  
  fp=fopen(fname,"wt");
  if (!fp) {
    sprintf(ferrormsg,"fsave: unable to create file '%s'\n",fname);
    return -1;  
  }
  np=f->np*(f->ni > 1 ? f->ni : 1);
  data=f->data;
  if (f->type == FD_TYPE_SPE) {
    for (i=1;i<=np;i++) {
      fprintf(fp,"%g %g %g\n",FD_FREQ(f,i),data[i].re,data[i].im);
    }
  } else {
    for (i=1;i<=np;i++) {
      fprintf(fp,"%g %g %g\n",FD_TIME(f,i),data[i].re,data[i].im);
    }  
  }
  fclose(fp);
  return 0;
}
/*
#define FREQ(i,np,sw,ref) (sw*(((i)-1)/(double)(np)-0.5)+ref)
#define TIME(i,sw) ((double)(i-1)/(sw))
*/

int fsave_xyreim(FD* f,char* fname)
{
  FILE* fp;
  int i,j,np,ni;
  complx* data;
  double sw1,ref1;
  
  fp=fopen(fname,"w");
  if (!fp) {
    sprintf(ferrormsg,"fsave: unable to create file '%s'\n",fname);
    return -1;  
  }
  np=f->np;
  ni=f->ni;
  data=f->data;
  sw1=f->sw1;
  ref1=f->ref1;
  
  if (f->type == FD_TYPE_SPE) {
    for (j=1;j<=ni;j++) {
      for (i=1;i<=np;i++) {
         complx* p = &data[(j-1)*np+i];
         fprintf(fp,"%g %g %g %g\n",FREQ(j,ni,sw1,ref1),FD_FREQ(f,i),p->re,p->im);
      }
      fprintf(fp,"\n");
    }    
  } else {
    for (j=1;j<=ni;j++) {
      for (i=1;i<=np;i++) {
         complx* p = &data[(j-1)*np+i];
         fprintf(fp,"%g %g %g %g\n",TIME(j,sw1),FD_TIME(f,i),p->re,p->im);
      }
      fprintf(fp,"\n");
    }
  }
  fclose(fp);
  return 0;
}

/*    
 Use 
   fsave $f name.gnudata -gnu2d -binary

 and a gnuplot file like this:
 
  set param
  set view 0,0,1
  set cntrparam bspline
  set cntrparam levels 10
  set nosurface
  set xlabel 'Xlabel'
  set ylabel 'Ylabel'
  set contour
  set term x11
  splot 'name.gnudata' binary w l
  pause -1

 From the GNUPLOT help file:
 
     Single precision floats are stored in a binary file as follows:

       <N+1>  <y0>   <y1>   <y2>  ...  <yN>
        <x0> <z0,0> <z0,1> <z0,2> ... <z0,N>
        <x1> <z1,0> <z1,1> <z1,2> ... <z1,N>
         :      :      :      :   ...    :

     which are converted into triplets:
       <x0> <y0> <z0,0>
       <x0> <y1> <z0,1>
       <x0> <y2> <z0,2>
        :    :     :
       <x0> <yN> <z0,N>

       <x1> <y0> <z1,0>
       <x1> <y1> <z1,1>
       :    :     :

     These triplets are then converted into `gnuplot` iso-curves and then
     `gnuplot` proceeds in the usual manner to do the rest of the plotting.
 
 
*/    

int fsave_gnu2d(FD* f,char* fname,int is_binary)
{
  if (!is_binary) {
    fsave_xyreim(f,fname);
    return 0;
  } else {
  
    FILE* fp;
    int i,j,np,ni;
    complx* data;
    double sw1,ref1;
    float* v;
    fp=fopen(fname,"w");
    if (!fp) {
      sprintf(ferrormsg,"fsave: unable to create file '%s'\n",fname);
      return -1;  
    }
    np=f->np;
    ni=f->ni;
    data=f->data;
    sw1=f->sw1;
    ref1=f->ref1;
    v=(float*)malloc(sizeof(float)*(np+1));
    if (!v) {
      fprintf(stderr,"error: unable to allocate %d floats\n",np+1);
      exit(1);
    }
    
    if (f->type == FD_TYPE_SPE) {
      v[0]=(float)np;
      for (i=1;i<=np;i++) {
         v[i]= FD_FREQ(f,i);
      }
      if (fwrite(v,sizeof(float),np+1,fp) != np+1) {
         fprintf(stderr,"error: unable to write %d floats to file '%s'\n",np+1,fname);
         exit(1);
      }
      for (j=1;j<=ni;j++) {
        v[0] = FREQ(j,ni,sw1,ref1);
        for (i=1;i<=np;i++) {
           v[i]= data[(j-1)*np+i].re;
        }
        if (fwrite(v,sizeof(float),np+1,fp) != np+1) {
           fprintf(stderr,"error: unable to write %d floats to file '%s'\n",np+1,fname);
           exit(1);
        }
      }    
    } else {
      v[0]=(float)np;
      for (i=1;i<=np;i++) {
         v[i]= FD_TIME(f,i);
      }
      if (fwrite(v,sizeof(float),np+1,fp) != np+1) {
         fprintf(stderr,"error: unable to write %d floats to file '%s'\n",np+1,fname);
         exit(1);
      }
      for (j=1;j<=ni;j++) {
        v[0] = TIME(j,sw1);
        for (i=1;i<=np;i++) {
           v[i]= data[(j-1)*np+i].re;
        }
        if (fwrite(v,sizeof(float),np+1,fp) != np+1) {
           fprintf(stderr,"error: unable to write %d floats to file '%s'\n",np+1,fname);
           exit(1);
        }
      }    
    }
    free(v);
    fclose(fp);
  }
  return 0;
}


int tclFSave(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  FD* f;
  int i,fidN,is_xreim=0,is_binary=0,is_double=0,is_gnu2d=0,is_xyreim=0,rmn=0,pipe=0,phsens=0;
  int is_raw_bin=0;
  /* ZT: -raw_bin format is just binary file without any headers and trailers */

  if (argc < 3)
    return TclError(interp,"Usage: fsave <data set> <filename> ?-xreim? ?-xyreim? ?-gnu2d? ?-binary? ?-double? ?-rmn? ?-nmrpipe? ?-phsens? ?-raw_bin?");
/* xyreim */
  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"fsave: argument 1 must be integer <data set>");

  for (i=3;i<argc;i++) {
    if (!strcmp(argv[i],"-xreim")) {
      is_xreim=1;
    } else if (!strcmp(argv[i],"-binary")) {
      is_binary=1;
    } else if (!strcmp(argv[i],"-double")) {
      is_double=1;
    } else if (!strcmp(argv[i],"-gnu2d")) {
      is_gnu2d=1;
    } else if (!strcmp(argv[i],"-xyreim")) {
      is_xyreim=1;
    } else if (!strcmp(argv[i],"-rmn")) {
      rmn=1;
    } else if (!strcmp(argv[i],"-nmrpipe")) {
      pipe=1;
    } else if (!strcmp(argv[i],"-magnitude")) {
      phsens=0;
    } else if (!strcmp(argv[i],"-states")) {
      phsens=1;
    } else if (!strcmp(argv[i],"-tppi")) {
      phsens=2;
    } else if (!strcmp(argv[i],"-rance-kay")) {
      phsens=3;
    } else if (!strcmp(argv[i],"-raw_bin")) {
      is_raw_bin=1;
    } else
      return TclError(interp,"fsave: argument %d must be -xreim, -xyreim, -gnu2d, -rmn, -binary and/or -double, -nmrpipe, -raw_bin\n",i);
  }

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fsave: data set %d was not previously loaded\n",fidN);

  f=fd[fidN];
  /* RA: worker process, ie. processes with rank>0, do not write file */
  //if (f->rank > 0){
  //  return 0;
  //}
  // ZT: obsolete in new master-slaves format

  if (rmn) {
    FD_write_rmn(argv[2], f);
  } else if (pipe) {
    FD_write_nmrpipe(argv[2], f, phsens);
  } else if (is_raw_bin) {
    FD_write_raw_bin(argv[2], f);
  } else if (is_xreim) {
    if (is_binary || is_double)
      return TclError(interp,"fsave: -xreim cannot be combined with -binary or -double\n");

    if (fsave_xreim(f,argv[2]) == -1)
      return TclError(interp,ferrormsg);

  } else if (is_xyreim) {

    if (is_double)
      return TclError(interp,"fsave: -double cannot be combined with -xyreim\n");
    if (is_binary)
      return TclError(interp,"fsave: -binary cannot be combined with -xyreim\n");

    if (f->sw1 <= 0.0 || f->ni < 1)
      return TclError(interp,"fsave: it must be 2d data to use -xyreim (use 'fset')\n");

    fsave_xyreim(f,argv[2]);
  
  } else if (is_gnu2d) {

    if (is_double)
      return TclError(interp,"fsave: -double cannot be combined with -gnu2d\n");

    if (f->sw1 <= 0.0 || f->ni < 1)
      return TclError(interp,"fsave: it must be 2d data to use -gnu2d (use 'fset')\n");

    fsave_gnu2d(f,argv[2],is_binary);
  } else {
    if (is_double && !is_binary)
      return TclError(interp,"fsave: -double cannot be used without specifying -single\n");
    if (fsave(f,argv[2],is_binary,is_double) == -1)
      return TclError(interp,"%s",ferrormsg);

  }
  return TCL_OK;
} 

int fdup(int fidN)
{
  FD* f;
  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL) {
    sprintf(ferrormsg,"fdup: data set %d was not previously loaded\n",fidN);
    return -1;
  }
  f=FD_dup(fd[fidN]);
  if (!f) return -1;
  return fnew(f);
}

int fdupzero(int fidN)
{
  FD* f;
  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL) {
    sprintf(ferrormsg,"fdupzero: data set %d was not previously loaded\n",fidN);
    return -1;
  }
  f=FD_dup(fd[fidN]);
  if (!f) return -1;
  memset(f->data,0,sizeof(complx)*(f->np*(f->ni > 1 ? f->ni : 1)+1));
  return fnew(f);
}

int check_2n(int np)
{
   int n,n2;
   n=np;
   for (;;) {
      n2=n/2;
      if (2*n2 != n) return 0;
      if (n2 == 1) break;
      n=n2;
   }
   return 1;
}

int round_2n(int np)
{
   int n,n2,nn;
   n=np;
   nn=1;
   for (;;) {
      nn *= 2;
      n2=n/2;
      if (n2 == 1) break;
      n=n2;
   }
   if (abs(nn*2-np) < abs(nn-np)) nn *= 2; else
   if (abs(nn/2-np) < abs(nn-np)) nn /= 2;
   return nn;
}


int tclFZerofill(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN;
  int npz;
  complx *d;
  FD *f;

  if (argc != 3 && argc != 4)
    return TclError(interp,"Usage: fzerofill <data set> <np/zerofill upto> ?<ni/zerofill upto>?");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"fzerofill: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fzerofill: data set %d was not previously loaded\n",fidN);

  if (Tcl_GetInt(interp,argv[2],&npz) == TCL_ERROR)
    return TclError(interp,"fzerofill: argument 2 must be integer <np/zerofill upto>");

  f=fd[fidN];
  npz=round_2n(npz);

  if (argc == 4) {
    int i,np,niz;
    
    if (Tcl_GetInt(interp,argv[3],&niz) == TCL_ERROR)
      return TclError(interp,"fzerofill: argument 3 must be integer <ni/zerofill upto>");
    if (f->ni <= 1)
      return TclError(interp,"fzerofill: 2D zerofilling requires 'ni' to be set (use 'fset')");

    np=f->np;
    niz=round_2n(niz);


    if (npz < np) npz=np;
    if (niz < f->ni) niz=f->ni;

    if (npz > np || niz > f->ni) {
      d=(complx*)malloc(sizeof(complx)*(niz*npz+1));
      if (!d)
        return TclError(interp,"fzerofill: unable to allocate %d complex numbers\n",niz*npz);
      memset(d,0,sizeof(complx)*(niz*npz+1));
      for (i=0;i<f->ni;i++) {
        memcpy(&d[i*npz+1],&f->data[i*np+1],sizeof(complx)*np);
      }
      free(f->data);
      f->ni=niz;
      f->np=npz;
      f->data=d;      
    }
  } else {
    if (f->ni > 1)
        return TclError(interp,"fzerofill: 1D zerofilling requires 'ni' to be zero or one\n");

    if (npz > f->np) {
      d=(complx*)malloc(sizeof(complx)*(npz+1));
      if (!d)
        return TclError(interp,"fzerofill: unable to allocate %d complex numbers\n",npz);

      memset(d,0,sizeof(complx)*(npz+1));
      memcpy(&d[1],&f->data[1],sizeof(complx)*f->np);
      free(f->data);
      f->np=npz;
      f->data=d;      
    }
  }
  return TCL_OK;
}


int tclFAddlb(ClientData data,Tcl_Interp* interp,int _argc, char *argv[])
{
  int i,fidN,phsens=0,sym=0,argc = _argc;
  double lineb,glfrac,top;
  FD *f;

  if (argc < 4) goto usage;
  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"faddlb: argument 1 must be integer <data set>");
  
  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
   return TclError(interp,"faddlb: data set %d was not previously loaded\n",fidN);
  f = fd[fidN];
  
  for (i=4; i<argc;i++) {
    if (!strcmp(argv[i],"-symmetric")) {
      sym = 1;
      argc--;
      for (; i<argc; i++) argv[i] = argv[i+1];
    }
  }

  top = 0.0;
  for (i=4; i<argc-1;i++) {
    if (!strcmp(argv[i],"-top")) {
      if (TclGetDoubleWithUnits(interp,argv[i+1],fd[fidN],&top,1,UNITS_S|UNITS_MS|UNITS_US) == TCL_ERROR)
       return TclError(interp,"faddlb: argument %d must be double (with units)",
         i);
      argc-=2;
      for (; i<argc; i++) argv[i] = argv[i+2];
    }
  }
  if (top < 1) top = 0;
  if (top > f->np) top = f->np;

  for (i=4; i<argc;i++) {
    if (!strcmp(argv[i],"-phsens")) {
      phsens=1;
      argc-=1;
      for (; i<argc; i++) argv[i] = argv[i+1];
    }
  }
  if (argc != 4 && argc != 6 && argc != 7 && argc != 8) {
usage:
     return TclError(interp,"Usage: faddlb <data set> <linebroadening> <gauss/lorentz fraction> "
      " ?<linebroadening1> <gauss/lorentz fraction1>? ?-phsens? ?-symmetric? ?-top <toppoint>?");
  }
  
  if (Tcl_GetDouble(interp,argv[2],&lineb) == TCL_ERROR)
    return TclError(interp,"faddlb: argument 2 must be double <linebroadening>");

  if (Tcl_GetDouble(interp,argv[3],&glfrac) == TCL_ERROR)
   return TclError(interp,"faddlb: argument 3 must be double <gauss/lorentz fraction>");

  if (glfrac < 0.0 || glfrac > 1.0)
    return TclError(interp,"faddlb: gauss/lorentz fraction out of range:"
                           " %g (must be between 0 and 1)",glfrac);
  f=fd[fidN];

  if (f->sw <= 0.0)
    return TclError(interp,"faddlb: sw must be different from zero");

  if (f->ni > 1) {
    double lineb1,glfrac1;

/*
    if (f->ni <= 1)
      return TclError(interp,"faddlb: 2D linebroadening requires 'ni' to be set (use 'fset')");
*/    
    if (Tcl_GetDouble(interp,argv[4],&lineb1) == TCL_ERROR)
      return TclError(interp,"faddlb: argument 4 must be double <linebroadening1>");

    if (Tcl_GetDouble(interp,argv[5],&glfrac1) == TCL_ERROR) 
      return TclError(interp,"faddlb: argument 5 must be double <gauss/lorentz fraction1>");

    if (f->sw1 <= 0.0)
      return TclError(interp,"faddlb: sw1 must be different from zero when adding"
      " linebroadening to 2D data");

/*
    phsens=0;
    if (argc == 7) {
      phsens=1;
      if (strcmp(argv[6],"-phsens"))
        return TclError(interp,"faddlb: argument 6 must be -phsens or nothing\n");
    }
*/
    daddlb2d(f->data,f->np,f->ni,f->sw,lineb,glfrac,f->sw1,lineb1,glfrac1,phsens,top);
  } else {
    if (f->ni > 1)
        return TclError(interp,"faddlb: 1D linebroadining requires 'ni' to be zero or one\n");

    daddlb((complx*)f->data,f->np,f->sw,lineb,glfrac,sym,top);
  }
  return TCL_OK;
}


int tclFDup(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN,fidN2;
  char buf1[16], buf2[128];

  if (argc != 2)
     return TclError(interp,"Usage: <data set> fdup <data set>");


  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"fdup: argument 1 must be integer <data set>");

  if ((fidN2=fdup(fidN)) == -1) 
     return TclError(interp,"%s",ferrormsg);

  sprintf(buf1, "%d", fidN);
  sprintf(buf2, "%ld", (long)fd[fidN]);
  Tcl_Eval(interp,"namespace eval FD_Internal {variable f}");
  Tcl_SetVar2(interp,"FD_Internal::f",buf1,buf2,TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG);
  TclSetResult(interp,"%d",fidN2);
  return TCL_OK;
} 


int tclFDupZero(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN,fidN2;
  char buf1[16], buf2[128];

  if (argc != 2)
    return TclError(interp,"Usage: <data set> fdupzero <data set>");


  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"fdupzero: argument 1 must be integer <data set>");


  if ((fidN2=fdupzero(fidN)) == -1)
    return TclError(interp,"%s",ferrormsg);
  sprintf(buf1, "%d", fidN);
  sprintf(buf2, "%ld", (long)fd[fidN]);
  Tcl_Eval(interp,"namespace eval FD_Internal {variable f}");
  Tcl_SetVar2(interp,"FD_Internal::f",buf1,buf2,TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG);
  TclSetResult(interp,"%d",fidN2);
  TclSetResult(interp,"%d",fidN2);
  return TCL_OK;
} 

int fsamesize(FD* fd1,FD* fd2)
{
  if (fd1->np != fd2->np || fd1->ni != fd2->ni) {
     strcpy(ferrormsg,"the data sets have different number of data points");
     return 0;
  }
  if (fabs(fd1->sw-fd2->sw) > 1.0e-5*fabs(fd1->sw)) {
     sprintf(ferrormsg,"the data sets have different spectral-widths (%g) and (%g)",
     fd1->sw, fd2->sw);
     return 0;
  }
  if (fabs(fd1->sw1-fd2->sw1) > 1.0e-5*fabs(fd1->sw1)) {
     sprintf(ferrormsg,"the data sets have different indirect spectral-widths (%g) and (%g)",
     fd1->sw1, fd2->sw1);
     return 0;
  }
  if (fabs(fd1->ref-fd2->ref) > 1.0e-7*fabs(fd1->sw)) {
     sprintf(ferrormsg,"the data sets have different reference values (%g) and (%g)",
     fd1->ref, fd2->ref);
     return 0;
  }
  if (fabs(fd1->ref1-fd2->ref1) > 1.0e-5*fabs(fd1->sw1)) {
     sprintf(ferrormsg,"the data sets have different indirect reference values (%g) and (%g)",
     fd1->ref1, fd2->ref1);
     return 0;
  }
  return 1;
}

int tclFCopy(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN,fidN2;
  
  if (argc != 3)
    return TclError(interp,"Usage: fcopy <desc-to> <desc-from>");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"fcopy: argument 1 must be integer <desc-to>");

  if (Tcl_GetInt(interp,argv[2],&fidN2) == TCL_ERROR) 
    return TclError(interp,"fcopy: argument 1 must be integer <desc-from>");

  if (!fsamesize(fd[fidN],fd[fidN2]))
    return TclError(interp,"fcopy: %s\n",ferrormsg);

  memcpy(&(fd[fidN]->data[1]),&(fd[fidN2]->data[1]),sizeof(complx)*(fd[fidN]->np)*(fd[fidN]->ni > 1 ? fd[fidN]->np : 1));
  return TCL_OK;
} 

double amin(double a, double b, double c)
{
  double min = a<b ? a:b;
  return min<c ? min:c;
}
double amid(double a, double b, double c)
{
  double ab = a>b ? a:b;
  double ac = a>c ? a:c;
  double bc = b>c ? b:c;

  return amin(ab, ac, bc);
}
double amax(double a, double b, double c)
{
  double max = a>b ? a:b;
  return max>c ? max:c;
}

int freq2i(double f, double fmin, double df)
{
  return (int)((f-fmin)/df)+1;
}

double i2freq(int i, double fmin, double df)
{
  return df*(i-1)+fmin;
}

void addpoint(complx* vec, int i, int np, complx w)
{
  int ii = i;
  while (ii < 1) ii+= np;
  while (ii > np) ii-= np;
  vec[ii].re += w.re;
  vec[ii].im += w.im;
}

int tclFAddtriangle(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN, np, i, imin, imid, imax;
  double f1, f2, f3, xmin, xmid, xmax, ixmin, ixmid, ixmax;
  double fmin, sw, df;
  complx a, y, weight, area, hmin, hmax, sum;
  complx *vec;

  if (argc < 7)
    return TclError(interp,"Usage: faddtriangle <desc> freq1 freq2 freq3 weight_re weight_im");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"faddtriangle: argument 1 must be integer <data set from>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
      return TclError(interp,"faddtriangle: data set %d was not previously loaded\n",fidN);

  if (fd[fidN]->ni > 1)
    return TclError(interp,"faddtriangle: command cannot be used for 2D data.");

  if (fd[fidN]->type != FD_TYPE_SPE)
    return TclError(interp,"faddtriangle: command only operates on descriptors of type 'spe'.");

  if (Tcl_GetDouble(interp,argv[2],&f1) == TCL_ERROR) 
    return TclError(interp, "faddtriangle: cannot interpret argument 3 as double");

  if (Tcl_GetDouble(interp,argv[3],&f2) == TCL_ERROR) 
    return TclError(interp, "faddtriangle: cannot interpret argument 4 as double");

  if (Tcl_GetDouble(interp,argv[4],&f3) == TCL_ERROR) 
    return TclError(interp, "faddtriangle: cannot interpret argument 5 as double");

  if (Tcl_GetDouble(interp,argv[5],&weight.re) == TCL_ERROR) 
    return TclError(interp, "faddtriangle: cannot interpret argument 6 as double");

  if (Tcl_GetDouble(interp,argv[6],&weight.im) == TCL_ERROR) 
    return TclError(interp, "faddtriangle: cannot interpret argument 7 as double");

  vec = fd[fidN]->data;
  np = fd[fidN]->np;
  sw = fd[fidN]->sw;
  df = sw/np;
  fmin = fd[fidN]->ref-sw/2.0;

  xmin = amin(f1,f2,f3);
  xmid = amid(f1,f2,f3);
  xmax = amax(f1,f2,f3);
  imin = freq2i(xmin, fmin, df);
  imid = freq2i(xmid, fmin, df);
  imax = freq2i(xmax, fmin, df);

  sum.re = 0;
  sum.im = 0;

  if (imin == imax) {
    return TCL_OK;
    addpoint(vec, imin, np, weight);
    return TCL_OK;
  }
  
  y = CRmul(weight, 2.0/(xmax-xmin)*df);
  ixmin = i2freq(imin, fmin, df);
  ixmid = i2freq(imid, fmin, df);
  ixmax = i2freq(imax, fmin, df);

/* Uphill */
  if (imid > imin) {
    a = CRmul(y, 1.0/(xmid-xmin));
    
    /* First point */
    hmax = CRmul(a, ixmin+df-xmin);
    area = CRmul(hmax, 0.5*(ixmin+df-xmin)/df);
    addpoint(vec, imin, np, area);
    
    /* Other points */
    for (i=2; i<= imid-imin; i++) {
      hmin = hmax;
      hmax = Cadd(hmin, CRmul(a, df));
      area = CRmul(Cadd(hmax,hmin), 0.5);
      addpoint(vec, imin+i-1, np, area);
    }
    /* Final point */
    area = CRmul(Cadd(hmax,y), 0.5*(xmid-ixmid)/df);
    addpoint(vec, imid, np, area);
  } else {
    area = CRmul(y, 0.5*(xmid-xmin)/df);
    addpoint(vec, imin, np, area);
  }
  
/* Downhill */
  if (imax > imid) {
    a = CRmul(y, 1.0/(xmax-xmid));
    
    /* First point */
    hmax = Csub(y, CRmul(a, ixmid+df-xmid));
    area = CRmul(Cadd(y,hmax), 0.5*(ixmid+df-xmid)/df);
    addpoint(vec, imid, np, area);
    /* Other points */
    for (i=2; i<= imax-imid; i++) {
      hmin = Csub(hmax, CRmul(a, df));
      area = CRmul(Cadd(hmin,hmax), 0.5);
      addpoint(vec, imid+i-1, np, area);
      hmax = hmin;
    }
    /* Final point */
    area = CRmul(hmax, 0.5*(xmax-ixmax)/df);
    addpoint(vec, imax, np, area);
  } else {
    area = CRmul(y, 0.5*(xmax-xmid)/df);
    addpoint(vec, imax, np, area);
  }

  return TCL_OK;
}


int tclFPhase(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN;
  int i,j;
  int is_rp=0,is_lp=0,is_scale=0,is_offset=0,is_dlp=0;
  double v,rp=0.0,lp=0.0,scale=1.0,offset=0.0,dlp=0.0;
  int nvec;
  complx *vec;

  if (argc < 4 || (argc % 2))
    return TclError(interp,"Usage: fphase <desc> [-rp <v>] [-lp <v>] [-scale <v>] [-offset <v>]");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"fphase: argument 1 must be integer <data set from>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
      return TclError(interp,"fphase: data set %d was not previously loaded\n",fidN);

/*
  if (fd[fidN]->ni > 1)
    return TclError(interp,"fphase: cannot phase 2D data. That is done by the 2D FFT.");
*/
  for (i=2;i<argc-1;i += 2) {
    if (fd[fidN]->type == FD_TYPE_SPE) {
      if (TclGetDoubleWithUnits(interp,argv[i+1],fd[fidN],&v,360,
         UNITS_S|UNITS_MS|UNITS_US) != TCL_OK)
	 return TclError(interp,"fphase: unable to convert '%s' to a double value",
         argv[i+1]);
    } else {
      if (TclGetDoubleWithUnits(interp,argv[i+1],fd[fidN],&v,360,
         UNITS_HZ|UNITS_KHZ|UNITS_PPM) != TCL_OK) 
	 return TclError(interp,"fphase: unable to convert '%s' to a double value",
         argv[i+1]);
    }

    if (!strcmp(argv[i],"-rp")) {
      is_rp=1;
      rp=v;
    } else if (!strcmp(argv[i],"-lp")) {
      is_lp=1;
      lp=v;
    } else if (!strcmp(argv[i],"-dlp")) {
      is_dlp=1;
      dlp=v;
    } else if (!strcmp(argv[i],"-scale")) {
      is_scale=1;
      scale=v;
    } else if (!strcmp(argv[i],"-offset")) {
      is_offset=1;
      offset=v;
    } else 
      return TclError(interp,"fphase: unknown option '%s'",argv[i]);
  }



  nvec = fd[fidN]->np;
  if (is_lp == 0 && is_rp == 0) {
    if (fd[fidN]->ni > 1) {
      for (j=0; j<fd[fidN]->ni; j++) {
        vec = (complx*) &fd[fidN]->data[j*fd[fidN]->np];
        dscale(vec,vec,nvec,scale,offset);
      }
    } else {
      vec = (complx*) fd[fidN]->data;
      dscale(vec,vec,nvec,scale,offset);
    }
  } else {
    if (fd[fidN]->ni > 1) {
      for (j=0; j<fd[fidN]->ni; j++) {
        vec = (complx*) &fd[fidN]->data[j*fd[fidN]->np];
        if (is_scale == 0 && is_offset == 0)
          dphase(vec,vec,nvec,rp,lp+(j*dlp));
        else
          dscalephase(vec,vec,nvec,scale,offset,rp,lp+(j*dlp));
      }
    } else {
      vec = (complx*) fd[fidN]->data;
      if (is_scale == 0 && is_offset == 0)
        dphase(vec,vec,nvec,rp,lp);
      else
        dscalephase(vec,vec,nvec,scale,offset,rp,lp);
    }
  }
  return TCL_OK;
} 


int tclFExtract(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  FD *f;
  int fidN,np,newnp;  
  complx *newdat,*dat;
  double xfrom,xto,newsw,dx;
  int i,j,ifrom,ito;

  if (argc != 4)
    return TclError(interp,"Usage: fextract <data set> <frq-from> <frq-to>");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"fextract: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fextract: data set %d was not previously loaded\n",fidN);

  if (fd[fidN]->type != FD_TYPE_SPE)
    return TclError(interp,"fextract: the data data set must contain spectrum type data.");

  if (fd[fidN]->ni > 1)
    return TclError(interp,"fextract: cannot extract from 2D data.");
 
  if (Tcl_GetDouble(interp,argv[2],&xfrom) == TCL_ERROR)
    return TclError(interp,"fextract: argument 2 must be double <frq-from>");

  if (Tcl_GetDouble(interp,argv[3],&xto) == TCL_ERROR)
    return TclError(interp,"fextract: argument 3 must be double <frq-to>");

  f=fd[fidN];
  np = f->np;
  dat = (complx*) f->data;

  if (xfrom > xto) {
    dx=xfrom;
    xfrom=xto;
    xto=dx;
  }
  ifrom=FD_INDEX(f,xfrom);
  if (ifrom < 1) ifrom=1;
  else if (ifrom > np) ifrom=np;
  ito=FD_INDEX(f,xto);
  if (ito < 1) ito=1;
  else if (ito > np) ito=np;
  if (ifrom >= ito) {
    /* extracting no points, but does'nt change the data.*/
    return TCL_OK;
  }
  if (ifrom == 1 && ito == np) {
    /* extracting all points*/
    return TCL_OK;
  }
  xfrom = FD_FREQ(f,ifrom);
  xto   = FD_FREQ(f,ito);
  newnp= ito - ifrom;
  newdat=(complx*)malloc(sizeof(complx)*(newnp+1));
  if (!newdat)
    return TclError(interp,"fextract: unable to allocate %d times 2 doubles\n",newnp+1);

  dx = FD_DELTAFREQ(f);
  newsw = newnp * dx;

  for (j=1,i=ifrom;i<ito;i++,j++) {
     newdat[j]=dat[i];
  }
  f->ref = (xto+xfrom)/2.0;
  free(f->data);
  f->data=newdat;
  f->np=newnp;
  f->sw=newsw;  
  return TCL_OK;
}  


#define PART_COMPLEX 0
#define PART_RE 1
#define PART_IM 2

void daddrms(complx* vec, complx* vec2,int from,int to,int part,double* sumrms,double* sumint)
{
  int i;
  double sre,sim,dre,dim;

  if (part == PART_COMPLEX) {
    for (i=from;i<=to;i++) {
      sre=vec[i].re+vec2[i].re;
      sim=vec[i].im+vec2[i].im;
      dre=vec[i].re-vec2[i].re;
      dim=vec[i].im-vec2[i].im;
      *sumint += sre*sre+sim*sim;
      *sumrms += dre*dre+dim*dim;
    }
  } else if (part == PART_RE) {
    for (i=from;i<=to;i++) {
      sre=vec[i].re+vec2[i].re;
      dre=vec[i].re-vec2[i].re;
      *sumint += sre*sre;
      *sumrms += dre*dre;
    }
  } else if (part == PART_IM) {
    for (i=from;i<=to;i++) {
      sim=vec[i].im+vec2[i].im;
      dim=vec[i].im-vec2[i].im;
      *sumint += sim*sim;
      *sumrms += dim*dim;
    }
  }
}

int tclFRealrms(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  FD* f;
  int fidN,fidN2;
  int i,nvec;
  complx *vec,*vec2;
  double v1,sumrms;

  if (argc != 3)
    return TclError(interp,"Usage: frealrms <descr 1> <descr 2>");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"frms: argument 1 must be integer <data set from>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"frms: data set %d was not previously loaded\n",fidN);

  if (Tcl_GetInt(interp,argv[2],&fidN2) == TCL_ERROR)
    return TclError(interp,"frms: argument 2 must be integer <data set to>");

  if (fidN2 < 1 || fidN2 > nfd || fd[fidN2] == NULL)
    return TclError(interp,"frms: data set %d was not previously loaded\n",fidN);

  if (!fsamesize(fd[fidN],fd[fidN2]))
    return TclError(interp,"frms: %s\n",ferrormsg);

  f=fd[fidN];
  vec = (complx*) f->data;
  nvec = f->np*(f->ni > 1 ? f->ni : 1);
  vec2 = (complx*) fd[fidN2]->data;

  sumrms = 0;
  for (i=1; i<=nvec; i++) {
    v1 = vec[i].re - vec2[i].re;
    sumrms += v1*v1;
  }

  TclSetResult(interp,"%g",sqrt(sumrms/(double)nvec));
  return TCL_OK;
} 


int tclFRms(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  FD* f;
  int fidN,fidN2;
  int nvec;
  complx *vec,*vec2;
  double v1,v2;
  int part;
  char **par,**par2,*range;
  int i,i1,i2,npar,npar2;
  double sumrms,sumint;

  if (argc < 3 || argc > 6)
    return TclError(interp,"Usage: frms <descr 1> <descr 2> ?-re | -im? ?{{from to} {from to} ..}?");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"frms: argument 1 must be integer <data set from>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"frms: data set %d was not previously loaded\n",fidN);

  if (Tcl_GetInt(interp,argv[2],&fidN2) == TCL_ERROR)
    return TclError(interp,"frms: argument 2 must be integer <data set to>");

  if (fidN2 < 1 || fidN2 > nfd || fd[fidN2] == NULL)
    return TclError(interp,"frms: data set %d was not previously loaded\n",fidN);


  if (!fsamesize(fd[fidN],fd[fidN2]))
    return TclError(interp,"frms: %s\n",ferrormsg);

  f=fd[fidN];
  vec = (complx*) f->data;
  nvec = f->np*(f->ni > 1 ? f->ni : 1);
  vec2 = (complx*) fd[fidN2]->data;

  part=PART_COMPLEX;
  range=NULL;
  for (i=3;i<argc;i++) {
    if (!strcmp(argv[i],"-re")) {
      part=PART_RE;
    } else if (!strcmp(argv[i],"-im")) {
      part=PART_IM;
    } else {
      range=argv[i];
    } 
  }
  sumint=0;
  sumrms=0;
  if (range != NULL) {
    if (Tcl_SplitList(interp,range,&npar,&par) != TCL_OK)
      return TclError(interp,"frms: list is not formed correctly\n");

    for (i=0;i<npar;i++) {
      if (Tcl_SplitList(interp,par[i],&npar2,&par2) != TCL_OK)
        return TclError(interp,"frms: list element number %d is not formed correctly\n",i+1);
  
      if (npar2 != 2)
        return TclError(interp,"frms: list element number %d must contain two values, not %d\n",i+1,npar2);

      if (Tcl_GetDouble(interp,par2[0],&v1) != TCL_OK)
        return TclError(interp,"frms: unable to convert '%s' to a value in list element number %d\n",par2[0],i+1);

      if (Tcl_GetDouble(interp,par2[1],&v2) != TCL_OK)
        return TclError(interp,"frms: unable to convert '%s' to a value in list element number %d\n",par2[1],i+1);

      if (v1 >= v2)
        return TclError(interp,"frms: value 2 must be larger than value 1 in list element number %d\n",i+1);

      i1=FD_INDEX(f,v1);
      if (i1 < 1) i1=1; else if (i1 > nvec) i1=nvec;
      i2=FD_INDEX(f,v2);
      if (i2 < 1) i2=1; else if (i2 > nvec) i2=nvec;

      if (i1 < i2) {
         daddrms(vec,vec2,i1,i2,part,&sumrms,&sumint);
      }
      Tcl_Free((char*)par2);
    }
    Tcl_Free((char*)par);
  } else {
    daddrms(vec,vec2,1,nvec,part,&sumrms,&sumint);
  }
  TclSetResult(interp,"%g",1000.0*sumrms/sumint);
  return TCL_OK;
} 


void daddscale(complx* vec, complx* vec2,int from,int to,int part,double* sumxx,double* sumxy)
{
  int i;

  if (part == PART_COMPLEX) {
    for (i=from;i<=to;i++) {
      *sumxx += vec[i].re*vec[i].re;
      *sumxx += vec[i].im*vec[i].im;
      *sumxy += vec[i].re*vec2[i].re;
      *sumxy += vec[i].im*vec2[i].im;
    }
  } else if (part == PART_RE) {
    for (i=from;i<=to;i++) {
      *sumxx += vec[i].re*vec[i].re;
      *sumxy += vec[i].re*vec2[i].re;
    }
  } else if (part == PART_IM) {
    for (i=from;i<=to;i++) {
      *sumxx += vec[i].im*vec[i].im;
      *sumxy += vec[i].im*vec2[i].im;
    }
  }
}


int tclFAutoscale(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  FD* f;
  int fidN,fidN2;
  int nvec;
  complx *vec,*vec2;
  double v1,v2;
  int part;
  char **par,**par2,*range;
  int i,i1,i2,npar,npar2;
  double sumxx,sumxy,scale;

  if (argc < 3 || argc > 6)
    return TclError(interp,"Usage: fautoscale <descr 1> <descr 2> ?-re | -im? ?{{from to} {from to} ..}?");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"fautoscale: argument 1 must be integer <data set from>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fautoscale: data set %d was not previously loaded\n",fidN);

  if (Tcl_GetInt(interp,argv[2],&fidN2) == TCL_ERROR)
    return TclError(interp,"fautoscale: argument 2 must be integer <data set to>");

  if (fidN2 < 1 || fidN2 > nfd || fd[fidN2] == NULL)
    return TclError(interp,"fautoscale: data set %d was not previously loaded\n",fidN);


  if (!fsamesize(fd[fidN],fd[fidN2]))
    return TclError(interp,"fautoscale: %s\n",ferrormsg);

  f=fd[fidN];
  vec = (complx*) f->data;
  nvec = f->np*(f->ni > 1 ? f->ni : 1);
  vec2 = (complx*) fd[fidN2]->data;

  part=PART_COMPLEX;
  range=NULL;
  for (i=3;i<argc;i++) {
    if (!strcmp(argv[i],"-re")) {
      part=PART_RE;
    } else if (!strcmp(argv[i],"-im")) {
      part=PART_IM;
    } else {
      range=argv[i];
    } 
  }
  sumxx=0;
  sumxy=0;
  if (range != NULL) {
    if (Tcl_SplitList(interp,range,&npar,&par) != TCL_OK)
      return TclError(interp,"fautoscale: list is not formed correctly\n");

    for (i=0;i<npar;i++) {
      if (Tcl_SplitList(interp,par[i],&npar2,&par2) != TCL_OK)
        return TclError(interp,"fautoscale: list element number %d is not formed correctly\n",i+1);
  
      if (npar2 != 2)
        return TclError(interp,"fautoscale: list element number %d must contain two values, not %d\n",i+1,npar2);

      if (Tcl_GetDouble(interp,par2[0],&v1) != TCL_OK)
        return TclError(interp,"fautoscale: unable to convert '%s' to a value in list element number %d\n",par2[0],i+1);

      if (Tcl_GetDouble(interp,par2[1],&v2) != TCL_OK)
        return TclError(interp,"fautoscale: unable to convert '%s' to a value in list element number %d\n",par2[1],i+1);

      if (v1 >= v2)
        return TclError(interp,"fautoscale: value 2 must be larger than value 1 in list element number %d\n",i+1);

      i1=FD_INDEX(f,v1);
      if (i1 < 1) i1=1; else if (i1 > nvec) i1=nvec;
      i2=FD_INDEX(f,v2);
      if (i2 < 1) i2=1; else if (i2 > nvec) i2=nvec;

      if (i1 < i2) {
         daddscale(vec,vec2,i1,i2,part,&sumxx,&sumxy);
      }
      Tcl_Free((char *)par2);
    }
    Tcl_Free((char *)par);
  } else {
    daddscale(vec,vec2,1,nvec,part,&sumxx,&sumxy);
  }
  scale = sumxy/sumxx;

  for (i=1; i<=nvec;i++) {
    vec[i].re *= scale;
    vec[i].im *= scale;
  }
  TclSetResult(interp,"%g",scale);
  return TCL_OK;
} 

int tclFInt(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  FD* f;
  int fidN;
  int nvec;
  complx *vec;
  double v1,v2;
  char **par,**par2,buf[256];
  int i,j,i1,i2,npar,npar2;
  double area;

  if (argc != 3)
    return TclError(interp,"Usage: {area1 area2 ..} fint <descr> {{from1 to1} {from2 to2} ..}");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"fint: argument 1 must be integer <data set from>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fint: data set %d was not previously loaded\n",fidN);

  if (fd[fidN]->ni > 1)
    return TclError(interp,"fint: cannot integrate 2D data.");

  f=fd[fidN];
  vec = (complx*) f->data;
  nvec = f->np;

  Tcl_ResetResult(interp);
  if (Tcl_SplitList(interp,argv[2],&npar,&par) != TCL_OK)
    return TclError(interp,"fint: list is not formed correctly\n");

  for (i=0;i<npar;i++) {
    if (Tcl_SplitList(interp,par[i],&npar2,&par2) != TCL_OK)
      return TclError(interp,"fint: list element number %d is not formed correctly\n",i+1);
   
    if (npar2 != 2)
      return TclError(interp,"fint: list element number %d must contain two values, not %d\n",i+1,npar2);

    if (Tcl_GetDouble(interp,par2[0],&v1) != TCL_OK)
      return TclError(interp,"fint: unable to convert '%s' to a value in list element number %d\n",par2[0],i+1);

    if (Tcl_GetDouble(interp,par2[1],&v2) != TCL_OK)
      return TclError(interp,"fint: unable to convert '%s' to a value in list element number %d\n",par2[1],i+1);

    if (v1 >= v2)
      return TclError(interp,"fint: value 2 must be larger than value 1 in list element number %d\n",i+1);

    i1=FD_INDEX(f,v1);
    if (i1 < 1) i1=1; else if (i1 > nvec) i1=nvec;
    i2=FD_INDEX(f,v2);
    if (i2 < 1) i2=1; else if (i2 > nvec) i2=nvec;

    area=0;
    if (i1 < i2) {
       for (j=i1;j<=i2;j++) {
         area += vec[j].re;
       }
       area *=  (v2-v1)/(double)(i2-i1+1);
    }
    sprintf(buf,"%g",area);
    Tcl_AppendElement(interp,buf);
    Tcl_Free((char *)par2);
  }
  Tcl_Free(par);
  return TCL_OK;
} 


int tclFRev(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  FD *f;
  int fidN,np,np2;  
  complx *dat,tmp;
  int i,j;

  if (argc != 2)
    return TclError(interp,"Usage: frev <data set>");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"frev: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"frev: data set %d was not previously loaded\n",fidN);

  f=fd[fidN];
  np = f->np;
  dat = (complx*) f->data;
  f->ref = -f->ref;
  np2=np/2;
  if (fd[fidN]->ni > 1) {
    int k;
    for (k=0;k<fd[fidN]->ni;k++) {
      for (i=1;i<=np2;i++) {
        j=np-i+1;
        tmp=dat[j+k*np];
        dat[j+k*np]=dat[i+k*np];
        dat[i+k*np]=tmp;
      }
    }
  } else {
    for (i=1;i<=np2;i++) {
      j=np-i+1;
      tmp=dat[j];
      dat[j]=dat[i];
      dat[i]=tmp;
    }
  }
  return TCL_OK;
}  
    
int tclFZero(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  FD *f;
  int fidN,nvec;  
  complx* vec;
  double v1,v2;
  char **par,**par2;
  int i,i1,i2,npar,npar2;

  if (argc != 3 && argc != 2)
    return TclError(interp,"Usage: fzero <data set> ?{{from to} {from to} ..}?");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"fzero: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fzero: data set %d was not previously loaded\n",fidN);

  f=fd[fidN];

  if (argc == 2) {
    FD_zero(f);
    return TCL_OK;
  }
  if (f->ni > 1)
    return TclError(interp,"fzero: cannot zero ranges for a 2D data set\n");

  nvec = f->np;
  vec = (complx*) f->data;

  
  if (Tcl_SplitList(interp,argv[2],&npar,&par) != TCL_OK)
    return TclError(interp,"fzero: list is not formed correctly\n");

  for (i=0;i<npar;i++) {
    if (Tcl_SplitList(interp,par[i],&npar2,&par2) != TCL_OK)
      return TclError(interp,"fzero: list element number %d is not formed correctly\n",i+1);

    if (npar2 != 2)
      return TclError(interp,"fzero: list element number %d must contain two values, not %d\n",i+1,npar2);

    if (Tcl_GetDouble(interp,par2[0],&v1) != TCL_OK)
      return TclError(interp,"fzero: unable to convert '%s' to a value in list element number %d\n",par2[0],i+1);

    if (Tcl_GetDouble(interp,par2[1],&v2) != TCL_OK)
      return TclError(interp,"fzero: unable to convert '%s' to a value in list element number %d\n",par2[1],i+1);
      
    if (v1 >= v2)
      return TclError(interp,"fzero: value 2 must be larger than value 1 in list element number %d\n",i+1);

    i1=FD_INDEX(f,v1);
    if (i1 < 1) i1=1; else if (i1 > nvec) i1=nvec;
    i2=FD_INDEX(f,v2);
    if (i2 < 1) i2=1; else if (i2 > nvec) i2=nvec;
    if (i1 < i2) {
      memset(&vec[i1],0,sizeof(complx)*(i2-i1+1));
    }
    free(par2);
  }
  free(par);
  return TCL_OK;
}
   
/* Given four data points y1,y2,y3,y4 with x-values starting from x1
   and incremented with dx1, return the function value
   at x-value x using a polynomial of third degree as the interpolating
   function. Calculated using Lagrange's classical formula e.g. decribed in
   Numerical Recipes. */

complx pol3int(double x,double x1,double dx,
            complx y1, complx y2, complx y3, complx y4)
{
   double f1,f2,f3,f4,t;
   double t1,t2,t3,t4,dx3;
   complx r;
   
   t1=x-x1;
   t2=t1-dx;
   t3=t2-dx;
   t4=t3-dx;

   dx3=dx*dx*dx;

   t=t2*t3/(6*dx3);
   f1= -t4*t;
   f4=t1*t;
   t=t1*t4/(2*dx3);
   f2=t3*t;
   f3= -t2*t;

   r.re=f1*y1.re+f2*y2.re+f3*y3.re+f4*y4.re;
   r.im=f1*y1.im+f2*y2.im+f3*y3.im+f4*y4.im;
   return r;   
}
   
int tclFNewnp(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int i,j,np,newnp,np2,fidN;
  double dx1,x1,dx2,x2,x,r,newsw=0;
  complx* newdata,*dat;
  FD *f;

  if (argc < 3)
    return TclError(interp,"Usage: fnewnp <data set> <points>");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"fnewnp: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fnewnp: data set %d was not previously loaded\n",fidN);

/*
  if (fd[fidN]->ni > 1)
    return TclError(interp,"fnewnp: cannot operate on 2D data.");
*/

  if (Tcl_GetInt(interp,argv[2],&newnp) == TCL_ERROR)
    return TclError(interp,"fnewnp: argument 2 must be integer <points>");

  f=fd[fidN];
  np=f->np;
  dat=(complx*)f->data;
  
  if ( newnp == np )
    return TCL_OK;
  
  if (f->ni > 1) {
    newdata=(complx*)malloc(sizeof(complx)*(newnp+1)*f->ni);
  } else {
    newdata=(complx*)malloc(sizeof(complx)*(newnp+1));
  }
  if (!newdata)
    return TclError(interp,"fnewnp: unable to allocate %d times 2 doubles\n",newnp+1);

  if (f->type == FD_TYPE_FID) {
    newsw= f->sw * (double)(newnp-1)/(double)(np-1);
    x1=FD_TIME(f,1);
    dx1=FD_TIME(f,2)-x1;
    x2=TIME(1,newsw);
    dx2=TIME(2,newsw)-x2;
    np2=np+99;
  } else {
    x1=FD_FREQ(f,1);
    dx1=FD_FREQ(f,2)-x1;
    x2=FREQ(1,newnp,f->sw,f->ref);
    dx2=FREQ(2,newnp,f->sw,f->ref)-x2;
    np2=np-2;
  }

  r=dx2/dx1;

  if (f->ni > 1) {
    int k;
    for (k=0; k<f->ni;k++) {
      for (j=1;j<=newnp; j++) {
        i= (int)((double)(j-1)*r)+1;
        x= (j-1)*dx2+x2;
        if (i<2) {
          newdata[j+k*newnp] = pol3int(x, x1, dx1, dat[1+k*np],dat[2+k*np], dat[3+k*np], dat[4+k*np]);
        } else if (i>np2) {
          newdata[j+k*newnp] = pol3int(x, x1+np*dx1-dx1*3, dx1, dat[np-3+k*np],dat[np-2+k*np], dat[np-1+k*np], dat[np+k*np]);
        } else {
          newdata[j+k*newnp] = pol3int(x, (i-2)*dx1+x1, dx1, dat[i-1+k*np],dat[i+k*np], dat[i+1+k*np], dat[i+2+k*np]);
        }
      }
    }    
  } else {
    for (j=1;j<=newnp; j++) {
      i= (int)((double)(j-1)*r)+1;
      x= (j-1)*dx2+x2;
      if (i<2) {
        newdata[j] = pol3int(x, x1, dx1, dat[1],dat[2], dat[3], dat[4]);
      } else if (i>np2) {
        newdata[j] = pol3int(x, x1+np*dx1-dx1*3, dx1, dat[np-3],dat[np-2], dat[np-1], dat[np]);
      } else {
        newdata[j] = pol3int(x, (i-2)*dx1+x1, dx1, dat[i-1],dat[i], dat[i+1], dat[i+2]);
      }
    }    
  }
  free(f->data);
  f->data=newdata;
  f->np=newnp;
  if (f->type == FD_TYPE_FID) {
    f->sw=newsw;
  }
  return TCL_OK;
}


int tclFCreate(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  FD* f;
  int ac,is_sw=0,is_np=0,fidN;
  char buf1[16], buf2[128];
  
  if ((argc % 2) == 0)
    return TclError(interp,"Usage: <data set> fcreate -np v -sw v ?-type (spe|fid)? ?-ref v? ?-ni v? ?-sw1 v? ?-ref1 v?" );


  f=FD_alloc();
  if (!f)
    return TclError(interp,"fcreate: unable to allocate fid data set data structure\n");
  
  for (ac=1;ac<argc;ac += 2) {
    if (!strcmp(argv[ac],"-ref")) {
      if (Tcl_GetDouble(interp,argv[ac+1],&(f->ref)) == TCL_ERROR)
         return TclError(interp,"fcreate: argument %d must be double",ac+1);
    } else if (!strcmp(argv[ac],"-sw")) {
      if (Tcl_GetDouble(interp,argv[ac+1],&(f->sw)) == TCL_ERROR)
         return TclError(interp,"fcreate: argument %d must be double",ac+1);
      is_sw=1;
    } else if (!strcmp(argv[ac],"-sw1")) {
      if (Tcl_GetDouble(interp,argv[ac+1],&(f->sw1)) == TCL_ERROR)
        return TclError(interp,"fcreate: argument %d must be double",ac+1);
    } else if (!strcmp(argv[ac],"-ref")) {
      if (Tcl_GetDouble(interp,argv[ac+1],&(f->ref)) == TCL_ERROR)
        return TclError(interp,"fcreate: argument %d must be double",ac+1);
    } else if (!strcmp(argv[ac],"-ref1")) {
      if (Tcl_GetDouble(interp,argv[ac+1],&(f->ref1)) == TCL_ERROR)
        return TclError(interp,"fcreate: argument %d must be double",ac+1);
    } else if (!strcmp(argv[ac],"-np")) {
      if (Tcl_GetInt(interp,argv[ac+1],&(f->np)) == TCL_ERROR)
         return TclError(interp,"fcreate: argument %d must be integer",ac+1);
      is_np=1;
    } else if (!strcmp(argv[ac],"-ni")) {
      if (Tcl_GetInt(interp,argv[ac+1],&(f->ni)) == TCL_ERROR)
         return TclError(interp,"fcreate: argument %d must be integer",ac+1);
    } else if (!strcmp(argv[ac],"-type")) {
      if (!strcmp(argv[ac+1],"fid")) {
         f->type=FD_TYPE_FID;
      } else if (!strcmp(argv[ac+1],"spe")) {
         f->type=FD_TYPE_SPE;
      } else 
         return TclError(interp,"fcreate: argument to -type must be fid or spe\n");
    } else
      return TclError(interp,"fcreate: unknown option '%s'\n",argv[ac]);
  }
  if (!is_np)
    return TclError(interp,"fcreate: the -np option is required\n");
  if (f->np < 1)
    return TclError(interp,"fcreate: the argument to -np must be larger than zero\n");
  if (f->ni < 0)
    return TclError(interp,"fcreate: the argument to -ni must be larger than zero\n");
  if (!is_sw)
    return TclError(interp,"fcreate: the -sw option is required\n");
  if (f->sw <= 0.0)
    return TclError(interp,"fcreate: the argument to -sw must be larger than zero\n");
  if (f->sw1 < 0.0)
    return TclError(interp,"fcreate: the argument to -sw1 must be larger than zero\n");
  if (!FD_alloc1ddata(f))
    return TclError(interp,"fcreate: unable to allocate %d complex data points\n",f->np);
  FD_zero(f);
  fidN = fnew(f);
  sprintf(buf1, "%d", fidN);
  sprintf(buf2, "%ld", (long)fd[fidN]);
  Tcl_Eval(interp,"namespace eval FD {variable f}");
  Tcl_SetVar2(interp,"FD::f",buf1,buf2,TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG);
  Tcl_SetVar2(interp,"FD_Internal",buf1,buf2,TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG);
  TclSetResult(interp,"%d",fidN);
  return TCL_OK;
} 

int tclFSet(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN;
  int ac;
  FD* f;
  
  if (argc < 4 || (argc % 2) != 0) 
    return TclError(interp,"Usage: fset <data set> ?-ref v?  ?-sw v? ?-type (spe|fid)? ?-ni v? ?-ref1 v?  ?-sw1 v?");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"fset: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fset: data set %d was not previously loaded\n",fidN);

  f=fd[fidN];
  
  for (ac=2;ac<argc;ac += 2) {
    if (!strcmp(argv[ac],"-ref")) {
      if (Tcl_GetDouble(interp,argv[ac+1],&(f->ref)) == TCL_ERROR)
        return TclError(interp,"fset: argument %d must be double",ac+1);

    } else if (!strcmp(argv[ac],"-sw")) {
      if (Tcl_GetDouble(interp,argv[ac+1],&(f->sw)) == TCL_ERROR)
        return TclError(interp,"fset: argument %d must be double",ac+1);

    } else if (!strcmp(argv[ac],"-sw1")) {
      if (Tcl_GetDouble(interp,argv[ac+1],&(f->sw1)) == TCL_ERROR)
        return TclError(interp,"fset: argument %d must be double",ac+1);

    } else if (!strcmp(argv[ac],"-ref1")) {
      if (Tcl_GetDouble(interp,argv[ac+1],&(f->ref1)) == TCL_ERROR)
        return TclError(interp,"fset: argument %d must be double",ac+1);

    } else if (!strcmp(argv[ac],"-ni")) {
      if (Tcl_GetInt(interp,argv[ac+1],&(f->ni)) == TCL_ERROR)
        return TclError(interp,"fset: argument %d must be integer",ac+1);
      
    } else if (!strcmp(argv[ac],"-np")) {
      int np;
      if (Tcl_GetInt(interp,argv[ac+1],&np) == TCL_ERROR)
        return TclError(interp,"fset: argument %d must be integer",ac+1);
      if (np <= f->np && f->ni <= 1) {
	f->np = np;
      } else {
	return TclError(interp,"fset: check that you operate on 1D data and that the new np is smaller than the old");
      }
    } else if (!strcmp(argv[ac],"-type")) {
      if (!strcmp(argv[ac+1],"fid")) {
         f->type=FD_TYPE_FID;
      } else if (!strcmp(argv[ac+1],"spe")) {
         f->type=FD_TYPE_SPE;
      } else 
        return TclError(interp,"fset: argument to -type must be fid or spe\n");

    } else
      return TclError(interp,"fset: unknown option '%s'\n",argv[ac]);
  }
  return TCL_OK;
} 


int tclFGet(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN;

  if (argc != 3)
    return TclError(interp,"Usage: <result> fget <data set> [-ref | -sw | -np | -ni | -sw1 | -ref1 | -type]");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"fget: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fget: data set %d was not previously loaded\n",fidN);

  if (!strcmp(argv[2],"-ref")) {
    TclSetResult(interp,"%e",fd[fidN]->ref);
  } else if (!strcmp(argv[2],"-sw")) {
    TclSetResult(interp,"%e",fd[fidN]->sw);
  } else if (!strcmp(argv[2],"-ref1")) {
    TclSetResult(interp,"%e",fd[fidN]->ref1);
  } else if (!strcmp(argv[2],"-sw1")) {
    TclSetResult(interp,"%e",fd[fidN]->sw1);
  } else if (!strcmp(argv[2],"-np")) {
    TclSetResult(interp,"%d",fd[fidN]->np);
  } else if (!strcmp(argv[2],"-ni")) {
    TclSetResult(interp,"%d",fd[fidN]->ni);
  } else if (!strcmp(argv[2],"-type")) {
     if (fd[fidN]->type == FD_TYPE_SPE)      
       TclSetResult(interp,"spe");
     else
       TclSetResult(interp,"fid");
  } else 
    return TclError(interp,"fget: unknown option '%s'\n",argv[2]);

  return TCL_OK;
} 

int tclFIndex(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int i;
  int fidN;
  int ii = 0, ii2 = 0, indx, indx2, option = 0;
  int np;
  complx val;

  if (argc != 3 && argc != 4 && argc != 5)
    return TclError(interp,"Usage: <value> findex <data set> <point-number> ?-re|-im?");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"findex: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL) 
    return TclError(interp,"findex: data set %d was not previously loaded\n",fidN);

  for (i = 2; i<argc; i++) {
    if (*argv[i] == '-') {
      option = i;
    } else {
      if (ii == 0) ii = i;
      else ii2 = i;
    }
  }

  if (ii == 0) {
    return TclError(interp, "Usage: <value> findex <data set> <point-number> ?-re|-im?");
  }
  if (Tcl_GetInt(interp,argv[ii],&indx) == TCL_ERROR) 
    return TclError(interp,"findex: argument %d must be integer <point-number>",
    ii);
  if (ii2 != 0) {
    if (Tcl_GetInt(interp,argv[ii2],&indx2) == TCL_ERROR)
    return TclError(interp,"findex: argument %d must be integer <point-number>",
    ii2);
  }

  np=fd[fidN]->np;
  if (ii2==0 && fd[fidN]->ni>1) np *= fd[fidN]->ni;

  if (indx < 1 || indx > np)
    return TclError(interp,"findex: argument %d out of range (%d is not between 1 and %d)",
      ii,indx,np);
  if ((ii2 > 0) && (indx2 < 1 || indx2 > fd[fidN]->ni))
    return TclError(interp,"findex: argument %d out of range (%d is not between 1 and %d)",
      ii2,indx2,fd[fidN]->ni);
  
  if (ii2 > 0) indx += (indx2-1)*fd[fidN]->np;
  val=((complx*)(fd[fidN]->data))[indx];

  /* ZT: originally was like this:
  if (option > 0) {
     if (!strcmp(argv[option],"-re"))
        TclSetResult(interp,"%g",val.re);
     else if (!strcmp(argv[option],"-im"))
        TclSetResult(interp,"%g",val.im);
     else
        TclError(interp,"findex: argument %d must be either '-re' or '-im'",
	option);
  } else
    TclSetResult(interp,"%g %g",val.re,val.im);
  */
  /* ZT: here I change just format of output */
  if (option > 0) {
     if (!strcmp(argv[option],"-re"))
        TclSetResult(interp,"%.15g",val.re);
     else if (!strcmp(argv[option],"-im"))
        TclSetResult(interp,"%.15g",val.im);
     else
        TclError(interp,"findex: argument %d must be either '-re' or '-im'",
	option);
  } else
    TclSetResult(interp,"%.15g %.15g",val.re,val.im);
  
  return TCL_OK;
} 

int tclFX(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  FD* f;
  int fidN;
  int indx;
  int np;

  if (argc != 3)
    return TclError(interp,"Usage: <value> fx <data set> <point-number>");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"fx: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fx: data set %d was not previously loaded\n",fidN);

  if (fd[fidN]->ni > 1)
    return TclError(interp,"fx: cannot operate on 2D data.");

  if (Tcl_GetInt(interp,argv[2],&indx) == TCL_ERROR) 
    return TclError(interp,"fx: argument 2 must be integer <point-number>");

  f=fd[fidN];
  np=f->np;
  if (indx < 1 || indx > np) 
    return TclError(interp,"fx: argument 2 out of range (%d is not between 1 and %d)",indx,np);
    
  if (f->type == FD_TYPE_SPE)
    TclSetResult(interp,"%g",FD_FREQ(f,indx));
  else
    TclSetResult(interp,"%g",FD_TIME(f,indx));
  return TCL_OK;
} 



int tclFSetIndex(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN;
  int indx, indx2 = 1;
  int np, icor = 0;
  complx* val;
  double re,im;
  
  if (argc != 5 && argc != 6)
    return TclError(interp,"Usage: fsetindex <data set> <point-number> <real-value> <imag-value>");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"fsetindex: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fsetindex: data set %d was not previously loaded\n",fidN);

  if (Tcl_GetInt(interp,argv[2],&indx) == TCL_ERROR) 
    return TclError(interp,"fsetindex: argument 2 must be integer <point-number>");

  if (fd[fidN]->ni > 1) {
    if (argc == 6) {
      if (Tcl_GetInt(interp, argv[3],&indx2) == TCL_ERROR)
        return TclError(interp,"fsetindex: argument 3 must be integer <point-number>");
      icor = 1;
    }
  }

  if (Tcl_GetDouble(interp,argv[3+icor],&re) == TCL_ERROR) 
    return TclError(interp,"fsetindex: argument %d must be double <real-value>",
    3+icor);

  if (Tcl_GetDouble(interp,argv[4+icor],&im) == TCL_ERROR) 
    return TclError(interp,"fsetindex: argument %d must be double <imag-value>",
    4+icor);

  np=fd[fidN]->np;
  if (argc == 5 && fd[fidN]->ni > 1) np *= fd[fidN]->ni;

  if (indx < 1 || indx > np) 
    return TclError(interp,"fsetindex: argument 2 out of range (%d is not between 1 and %d)",indx,np);
  if (indx2 < 0 || (indx2 > fd[fidN]->ni && fd[fidN]->ni > 1)) 
    return TclError(interp,"fsetindex: argument 3 out of range (%d is not between 1 and %d)",
    indx2,fd[fidN]->ni);

  indx += (indx2-1)*fd[fidN]->np;

  val= &((complx*)(fd[fidN]->data))[indx];

  val->re=re;
  val->im=im;
  return TCL_OK;
} 



void fft2d(complx* data,int np,int ni,double sw,double sw1,double
     rp,double lp,double rp1,double lp1,double dlp,int phsens)
{
  void fft(complx* data,int np,int is);
  void dphase(complx* vec,complx* vec2,int nvec,double rp,double lp);
  int i,j,ni2;
  complx *d;

  for (i=0;i<ni;i++) {   
    d = &data[i*np];
    fft(d,np,1);
    if (rp != 0 || lp != 0.0 || dlp != 0) dphase(d,d,np,rp,(lp+dlp*i)); 
  }

  d=(complx*)malloc(sizeof(complx)*(ni+1));
  if (!d) {
    fprintf(stderr,"error: fft2d: unable to allocate vector\n");
    exit(1);
  }
  ni2=ni/2;

  for (i=1;i<=np;i++) {
    if (phsens) {
      for (j=0;j<ni2;j++) {
         d[j+1].re = data[2*j*np+i].re;
         d[j+1].im = data[(2*j+1)*np+i].re;
      }
      for (j=ni2;j<ni;j++) {
         d[j+1].re = 0.0;
         d[j+1].im = 0.0;
      }    
    } else {
      for (j=0;j<ni;j++) {
         d[j+1]= data[j*np+i];
         d[j+1].im=0;
      }
    }
    fft(d,ni,1);
    if (rp1 != 0 || lp1 != 0.0) dphase(d,d,ni,rp1,lp1); 
    for (j=0;j<ni;j++) {
       data[j*np+i] = d[j+1];
    }
  }
  free(d);
}    


int tclFFt(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  void fft(complx* data,int np,int is);
  int fidN;
  FD* f;

  if (argc != 3 && argc != 2 && argc != 6 && argc != 7 && argc != 8 &&
    argc != 9)
    return TclError(interp,
      "Usage:\n"
      "  1D FFT:  fft <data set> ?-inv?\n"
      "  2D FFT:  fft <data set> <rp> <lp1> <rp1> <lp1> ?-phsens? ?-dlp <dlp>?\n");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"fft: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fft: data set %d was not previously loaded\n",fidN);
    
  f=fd[fidN];

  if (argc == 3 || argc == 2) { 
    if (!check_2n(f->np))
      return TclError(interp,"fft: size of data to be fourier transformed"
                             " must be a fourier number 2^n but was %d",f->np);

    if (argc == 3) {
      if (strcmp(argv[2],"-inv"))
        return TclError(interp,"fft: unknown option '%s' (must be -inv or noting)\n",argv[2]);
      fft((complx*)f->data,f->np,-1);
    } else {
      fft((complx*)f->data,f->np,1);
    }
  } else {
    double rp,lp,rp1,lp1, dlp = 0;
    int isphsens = 0;

    if (f->ni <= 1)
      return TclError(interp,"fft: data must contain 'ni' to be fouriertranformed as 2D data");

    if (f->sw1 <= 0.0)
      return TclError(interp,"fft: data must contain 'sw1' to be fouriertranformed as 2D data");

    if (!check_2n(f->np))
      return TclError(interp,"fft: size of data to be fourier transformed"
                             " must be a fourier number 2^n but was %d",f->np);
    if (!check_2n(f->ni))
      return TclError(interp,"fft: size of data to be fourier transformed"
                             " must be a fourier number 2^n but was %d",f->ni);

    if (Tcl_GetDouble(interp,argv[2],&rp) == TCL_ERROR)
      return TclError(interp,"fft: argument 2 must be double <rp>");
    if (Tcl_GetDouble(interp,argv[3],&lp) == TCL_ERROR)
      return TclError(interp,"fft: argument 3 must be double <lp>");
    if (Tcl_GetDouble(interp,argv[4],&rp1) == TCL_ERROR)
      return TclError(interp,"fft: argument 4 must be double <rp1>");
    if (Tcl_GetDouble(interp,argv[5],&lp1) == TCL_ERROR)
      return TclError(interp,"fft: argument 5 must be double <lp1>");

    if (argc >= 7) {
      int i=6;
      while (i < argc) {
        if (!strcmp(argv[i], "-phsens")) isphsens = 1;
	else if (!strcmp(argv[i], "-dlp")) {
          if (i == argc-1) 
	    return TclError(interp,"fft: specify value to option -dlp");
          i++;
          if (TclGetDoubleWithUnits(interp,argv[i], fd[fidN],
	  &dlp,360,UNITS_S|UNITS_MS|UNITS_US) == TCL_ERROR)
            return TclError(interp,"fft: unable to convert argument %d ('%s') to a double", 
	      i, argv[i]);
	} else {
          return TclError(interp,"fft: argument %d ('%s') is not valid", 
	    i, argv[i]);
	}
	i++;
      }
    }
    fft2d(f->data,f->np,f->ni,f->sw,f->sw1,rp,lp,rp1,lp1,dlp,isphsens);
  }
  f->type= (f->type == FD_TYPE_FID ? FD_TYPE_SPE : FD_TYPE_FID);
  return TCL_OK;
} 



int tclFFt1d(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  void fft(complx* data,int np,int is);
  int fidN, i;
  FD* f;

  if (argc != 2 && argc != 3)
    return TclError(interp,
      "Usage:\n"
      "  fft1d <data set> ?-inv?\n");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"fft1d: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fft1d: data set %d was not previously loaded\n",fidN);
    
  f=fd[fidN];

  if (!check_2n(f->np))
    return TclError(interp,"fft1d: size of data to be fourier transformed"
                             " must be a fourier number 2^n but was %d",f->np);

  if (f->ni < 2) 
    return TclError(interp,"fft1d: data set must be two-dimensional\n");
    
  if (argc == 3) {
    if (strcmp(argv[2],"-inv"))
      return TclError(interp,"fft: unknown option '%s' (must be -inv or noting)\n",argv[2]);
    for (i=0; i<f->ni; i++) {
      fft((complx*)&f->data[i*f->np],f->np,-1);
    }
  } else {
    for (i=0; i<f->ni; i++) {
      fft((complx*)&f->data[i*f->np],f->np,1);
    }
  }
  f->type = (f->type == FD_TYPE_FID ? FD_TYPE_SPE : FD_TYPE_FID);
  return TCL_OK;
} 

#define SWAPI(a, b) {i = (a); (a) = (b); (b) = i;}
#define SWAPD(a, b) {dd = (a); (a) = (b); (b) = dd;}

int tclFTranspose(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN, i, j, ii, jj;
  FD* f;
  complx* d;
  double dd;
  
  if (argc != 2)
    return TclError(interp,
      "Usage:\n"
      "  ftranspose <data set>\n");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"ftranspose: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"ftranspose: data set %d was not previously loaded\n",fidN);
    
  f=fd[fidN];

  if (f->ni < 2) 
    return TclError(interp,"ftranspose: data set must be two-dimensional\n");
    
  d = (complx*)malloc(sizeof(complx)*(f->np*f->ni+1));
  if (!d) {
    fprintf(stderr,"error: ftranspose: unable to allocate data\n");
    exit(1);
  }
  for (i=0; i<f->np; i++) {
    for (j=0; j<f->ni; j++) {
      ii=i*f->ni+j+1;
      jj=i+j*f->np+1;
      d[ii].re = f->data[jj].re;
      d[ii].im = f->data[jj].im;
    }
  }
  free((char*)f->data);
  f->data = d;
  SWAPI(f->np,f->ni);
  SWAPI(f->type,f->type1);
  SWAPD(f->ref,f->ref1);
  SWAPD(f->sw,f->sw1);
  SWAPD(f->sfrq,f->sfrq1);
  return TCL_OK;
} 

#define FAC_GAUSS   (8.0*log(2.0))
#define FAC_LORENTZ (4.0)

double fgausslorentz(double x,double a0,double a1,double a2,double a3)
{
  double a,b;

  a0=fabs(a0);
  a2=fabs(a2);
  a3=fabs(a3);
  if (a3 > 1.0) a3=1.0;
  a= (x-a1)/a2;
  a *= a;
 
  if (a > 11) { /* exponential less than 5.68e-14 */
   /* lorentz is called a lot more times than gauss because the footer */
    b = a0*((1.0-a3)/(1.0+a*FAC_LORENTZ));
  } else {
    b = a0*(a3*exp(-0.5*FAC_GAUSS*a)+(1.0-a3)/(1.0+a*FAC_LORENTZ));
  }
  return b;
}

double fgausslorentz_area(double a0,double a2,double a3)
{
  a0=fabs(a0);
  a2=fabs(a2);
  a3=fabs(a3);
  if (a3 > 1.0) a3=1.0;
  return a2*a0*(M_PI*(1.0-a3)+sqrt(2.0*M_PI)*a3);
}

int tclFAddpeaks2D(ClientData data,Tcl_Interp* interp,int argc, char *argv[]);


int tclFAddpeaks(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  FD *f;
  int fidN;  
  complx* dat;
  double *x;
  double frq,in,lb,lg,cutoff;
  char **par,**par2,buf[128];
  int i,i0,k,npar,npar2,np;
  double re,relast;

  if (argc != 4)
    return TclError(interp,"Usage: faddpeaks <data set> <cutoff> "
                     "{{freq int linebroad gausslorentzratio} {...} ...}");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"faddpeaks: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL) 
    return TclError(interp,"faddpeaks: data set %d was not previously loaded\n",fidN);

  if (fd[fidN]->ni > 1)
    return tclFAddpeaks2D(data,interp,argc,argv);

  if (Tcl_GetDouble(interp,argv[2],&cutoff) != TCL_OK)  
    return TclError(interp,"faddpeaks: unable to convert '%s' in argument 2 to a cutoff value\n",argv[2]);

  if (Tcl_SplitList(interp,argv[3],&npar,&par) != TCL_OK) 
    return TclError(interp,"faddpeaks: list is not formed correctly\n");

  f=fd[fidN];
  np = f->np;
  dat = (complx*) f->data;

  x=(double*)malloc(sizeof(double)*(np+1));
  if (!x) 
    return TclError(interp,"faddpeaks: unable to allocate %d doubles\n",np+1);


  for (i=1;i<=np;i++) {
    x[i]=FD_FREQ(f,i);
    dat[i].im=0.0;
  }
  Tcl_ResetResult(interp);
  for (k=0;k<npar;k++) {
    if (Tcl_SplitList(interp,par[k],&npar2,&par2) != TCL_OK)
      return TclError(interp,"faddpeaks: list element number %d is not formed correctly\n",k+1);

    if (npar2 != 4)
      return TclError(interp,"faddpeaks: list element number %d must contain four values\n"
                             "{frequency intensity linebroadening gauss/lorentz-fraction}, not %d "
                             "values\n",k+1,npar2);

    if (Tcl_GetDouble(interp,par2[0],&frq) != TCL_OK)
      return TclError(interp,"faddpeaks: unable to convert '%s' to a value in list element"
                             " number %d\n",par2[0],k+1);

    if (Tcl_GetDouble(interp,par2[1],&in) != TCL_OK)
      return TclError(interp,"faddpeaks: unable to convert '%s' to a value in list element"
                             " number %d\n",par2[1],k+1);

    if (Tcl_GetDouble(interp,par2[2],&lb) != TCL_OK)
      return TclError(interp,"faddpeaks: unable to convert '%s' to a value in list element"
                             " number %d\n",par2[2],k+1);

    if (Tcl_GetDouble(interp,par2[3],&lg) != TCL_OK)
      return TclError(interp,"faddpeaks: unable to convert '%s' to a value in list element"
                             " number %d\n",par2[3],k+1);

    i0=FD_INDEX(f,frq);
    if (i0 < 4 || i0 > np) continue;
    
    for (i=i0-3;i<=i0+3;i++) {
        re = fgausslorentz(x[i],in,frq,lb,lg);
       dat[i].re += re;  
    }
    relast = fgausslorentz(x[i0+4],in,frq,lb,lg);
    dat[i0+4].re += relast;
    for (i=i0+6;i<=np;i += 2) {
       re = fgausslorentz(x[i],in,frq,lb,lg);
       dat[i].re += re;
       dat[i-1].re += (re+relast)*0.5;
       relast=re;        
       if (re < cutoff) break;
    }
    relast = fgausslorentz(x[i0-4],in,frq,lb,lg);
    dat[i0-4].re += relast;
    for (i=i0-6;i>0;i -= 2) {
       re = fgausslorentz(x[i],in,frq,lb,lg);
       dat[i].re += re;
       dat[i+1].re += (re+relast)*0.5;
       relast=re;
       if (re < cutoff) break;
    }
    sprintf(buf,"%g",fgausslorentz_area(in,lb,lg));
    Tcl_AppendElement(interp,buf);
    Tcl_Free((char*)par2);
  }
  Tcl_Free((char*)par);
  free((char*)x);
  return TCL_OK;
}

double fgausslorentz2D(double x,double y,double a0,double a3)
{
  double a,b,c;

  a0=fabs(a0);
  a3=fabs(a3);
  if (a3 > 1.0) a3=1.0;
  a= x*x;
  b= y*y;
 
  if (a+b > 11) { /* exponential less than 5.68e-14 */
   /* lorentz is called a lot more times than gauss because the footer */
    c = a0*((1.0-a3)/(1.0+a*FAC_LORENTZ)/(1.0+b*FAC_LORENTZ));
  } else {
    c = a0*(a3*exp(-0.5*FAC_GAUSS*(a+b))+(1.0-a3)/(1.0+a*FAC_LORENTZ)/(1.0+b*FAC_LORENTZ));
  }
  return c;
}




int tclFAddpeaks2D(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  FD *f;
  int fidN;  
  complx* dat;
  double *x, *y;
  double frq,in,lb,lg,cutoff,frq1,lb1;
  char **par,**par2;
  int i,i0,j0,k,npar,npar2,np,j,imin,imax,jmin,jmax,ni;
  
  if (argc != 4)
    return TclError(interp,"Usage: faddpeaks <data set> <cutoff> "
                     "{{freq freq1 int linebroad gausslorentzratio linebroad1} {...} ...}");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"faddpeaks: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL) 
    return TclError(interp,"faddpeaks: data set %d was not previously loaded\n",fidN);

  if (Tcl_GetDouble(interp,argv[2],&cutoff) != TCL_OK)  
    return TclError(interp,"faddpeaks: unable to convert '%s' in argument 2 to a cutoff value\n",argv[2]);

  if (Tcl_SplitList(interp,argv[3],&npar,&par) != TCL_OK) 
    return TclError(interp,"faddpeaks: list is not formed correctly\n");

  f=fd[fidN];
  np = f->np;
  ni = f->ni;
  dat = (complx*) f->data;

  x=(double*)malloc(sizeof(double)*(np+1));
  if (!x) 
    return TclError(interp,"faddpeaks: unable to allocate %d doubles\n",np+1);
  y=(double*)malloc(sizeof(double)*(ni+1));
  if (!y) 
    return TclError(interp,"faddpeaks: unable to allocate %d doubles\n",ni+1);


  for (i=1;i<=np;i++) {
    x[i]=FD_FREQ(f,i);
  }
  for (i=1;i<=ni;i++) {
    y[i]=FD_FREQ1(f,i);
  }
  for (i=1;i<=np*ni;i++) dat[i].im = 0;


  Tcl_ResetResult(interp);
  for (k=0;k<npar;k++) {
    if (Tcl_SplitList(interp,par[k],&npar2,&par2) != TCL_OK)
      return TclError(interp,"faddpeaks: list element number %d is not formed correctly\n",k+1);

    if (npar2 != 6)
      return TclError(interp,"faddpeaks: list element number %d must contain six values\n"
                             "{frequency frequency1 intensity lb lb1 gauss/lorentz-fraction}, not %d "
                             "values\n",k+1,npar2);

    if (Tcl_GetDouble(interp,par2[0],&frq) != TCL_OK)
      return TclError(interp,"faddpeaks: unable to convert '%s' to a value in list element"
                             " number %d\n",par2[0],k+1);

    if (Tcl_GetDouble(interp,par2[1],&frq1) != TCL_OK)
      return TclError(interp,"faddpeaks: unable to convert '%s' to a value in list element"
                             " number %d\n",par2[1],k+1);

    if (Tcl_GetDouble(interp,par2[2],&in) != TCL_OK)
      return TclError(interp,"faddpeaks: unable to convert '%s' to a value in list element"
                             " number %d\n",par2[2],k+1);

    if (Tcl_GetDouble(interp,par2[3],&lb) != TCL_OK)
      return TclError(interp,"faddpeaks: unable to convert '%s' to a value in list element"
                             " number %d\n",par2[3],k+1);

    if (Tcl_GetDouble(interp,par2[4],&lb1) != TCL_OK)
      return TclError(interp,"faddpeaks: unable to convert '%s' to a value in list element"
                             " number %d\n",par2[4],k+1);

    if (Tcl_GetDouble(interp,par2[5],&lg) != TCL_OK)
      return TclError(interp,"faddpeaks: unable to convert '%s' to a value in list element"
                             " number %d\n",par2[5],k+1);


    i0=FD_INDEX(f,frq);
    if (i0 < 1 || i0 > np) continue;
    
    j0=FD_INDEX1(f,frq1);
    if (j0 < 1 || j0 > ni) continue;

    imin = 1; jmin = 1;
    imax = np; jmax = ni;
    
    if (cutoff > 0) {
      for (i=1; i<=np; i*=2) {
        if (fgausslorentz2D(f->sw*i/np/lb,0,in,lg) < cutoff) break;
      }
      imin = i0-i-1; imax = i0+i+1;
      if (imin < 1) imin = 1;
      if (imax > np) imax = np;
      for (j=1; j<=ni; j*=2) {
        if (fgausslorentz2D(0,f->sw1*j/ni/lb1,in,lg) < cutoff) break;
      }
      jmin = j0-j-1; jmax = j0+j+1;
      if (jmin < 1) jmin = 1;
      if (jmax > ni) jmax = ni;
    }

    for (j=jmin-1;j<=jmax-1;j++) {
      for (i=imin;i<=imax;i++) {
        dat[j*np+i].re +=
	fgausslorentz2D((frq-x[i])/lb,(frq1-y[j+1])/lb1,in,lg);
      }
    }
    Tcl_Free((char *)par2);
  }
  Tcl_Free((char *)par);
  free((char*)x);
  free((char*)y);
  return TCL_OK;
}


int tclFFindPeaks(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN;
  int np,i1,i2;
  double sw,dx,sw2,x1= -1e10,x2 = 1e10;
  complx* dat;
  double dy,lastdy;
  int i,up,down,found,lastdown,lastup,lastswap,sensitivity;
  double x,threshold,maxy,y;
  char buf[256];

  if (argc != 6 && argc != 4) 
    return TclError(interp,"Usage: ffindpeaks <desc> <threshold 0-1> <sensitivity 1-> ?<from-freq> <to-freq>? ");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"ffindpeaks: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL) 
    return TclError(interp,"ffindpeaks: data set %d was not previously loaded\n",fidN);

  if (fd[fidN]->ni > 1)
    return TclError(interp,"ffindpeaks: cannot operate on 2D data.");

  np=fd[fidN]->np;
  sw=fd[fidN]->sw;
  dat=(complx*)(fd[fidN]->data);

  if (Tcl_GetDouble(interp,argv[2],&threshold) != TCL_OK)
    return TCL_ERROR;
  if (Tcl_GetInt(interp,argv[3],&sensitivity) != TCL_OK)
    return TCL_ERROR;
  if (sensitivity < 1) 
    return TclError(interp,"ffindpeaks: sensitivity parameter must be larger than zero");

  dx=sw/(double)np;
  sw2=sw/2.0;

  if (argc == 6) {
    if (Tcl_GetDouble(interp,argv[4],&x1) != TCL_OK)
      return TclError(interp,"ffindpeaks: unable to get double <from-freq>");

    if (Tcl_GetDouble(interp,argv[5],&x2) != TCL_OK)
      return TclError(interp,"ffindpeaks: unable to get double <to-freq>");

    i1= (x1+sw2)/dx+1;
    i2= (x2+sw2)/dx+1;
    if (i1 > np - 1) 
      return TclError(interp,"i1 is too large");

    if (i1 < 2) i1=2;
    if (i2 < 3) 
      return TclError(interp,"i2 is too small");

    if (i2 > np) i2=np;
  } else {
    i1=1;
    i2=np;
  }

  lastdy=dat[i1-1].re - dat[i1].re;
  maxy= -1e99;
  for (i=1;i<=np;i++) {
    if (dat[i].re > maxy) maxy=dat[i].re;
  }
  threshold *= maxy;
  found=0;
  up=0;
  down=0;
  lastup=0;
  lastdown=0;
  lastswap=0;
  Tcl_ResetResult(interp);
  for (i=i1;i<=i2;i++) {
    x = (i-1)*dx-sw2;
    dy = dat[i].re - dat[i-1].re;
    lastdy=dy;
    if (dy > 0.0) {
      if (down != 0) {
        lastswap=i-1;
        lastdown=down;
      }
      down=0;        
      up++;
      lastup=up;
    } else {
      if (up != 0) {
        lastswap=i-1;
        lastup=up;
      }
      up=0;
      down++;
      lastdown=down;
    }
    if (lastup >= sensitivity && lastdown >= sensitivity) {       
      if (lastswap) {
        if (down) {
          x = (lastswap-1)*dx-sw2;
          y = dat[lastswap].re;
          if (dat[i].re > threshold) {
            sprintf(buf,"%.20g %.20g",x,y);
            Tcl_AppendElement(interp,buf);
          }
        }
        lastswap=0;
      }
    }     
  }
  return TCL_OK;
}

int tclFMaxheight(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  char buf[256];
  int fidN,np,ni,i,i1,k,n;
  double maxy,t,t1;
  complx* dat;
  int abs=0;
  FD* f;
  
  if (argc != 2 && argc != 3) 
    return TclError(interp,"Usage: {height t/freq ?t1/freq1?} fmaxheight <desc> ?-abs?");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"fmaxheight: argument 1 must be integer <data set>");

  if (argc == 3) {
    if (strcmp(argv[2],"-abs"))
      return TclError(interp,"fmaxheight: argument 2 must be '-abs' or nothing");
    abs=1;
  }
  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL) 
    return TclError(interp,"fmaxheight: data set %d was not previously loaded\n",fidN);

  f= fd[fidN];
  np=f->np;
  ni=f->ni;
  dat=(complx*)(f->data);

  n=np*(ni > 1 ? ni : 1);
  if (abs) {
    maxy=0;
    k=0;
    for (i=1;i<=n;i++) {
      if (fabs(dat[i].re) > maxy) {
        maxy=fabs(dat[i].re);
        k=i;
      }
    }
  } else {
    maxy=-1e30;
    k=0;
    for (i=1;i<=n;i++) {
      if (dat[i].re > maxy) {
        maxy=dat[i].re;
        k=i;
      }
    }
  }
  if (k == 0)
    return TclError(interp,"fmaxheight: data set %d contains erroneous data\n",fidN);

  if (ni > 1) {
    i = ((k-1) % np) + 1;
    i1 = k / np + 1 ;
    if (fd[fidN]->type == FD_TYPE_FID) {
      t1=FD_TIME1(f,i1);
      t=FD_TIME(f,i);
    } else {
      t1=FD_FREQ1(f,i1);
      t=FD_FREQ(f,i);
    }
    sprintf(buf,"%g %g %g",maxy,t,t1);
  } else {
    if (fd[fidN]->type == FD_TYPE_FID) {
      t=FD_TIME(f,k);
    } else {
      t=FD_FREQ(f,k);
    }
    sprintf(buf,"%g %g",maxy,t);
  }
  Tcl_SetResult(interp,buf,TCL_VOLATILE);  
  return TCL_OK;
}


int tclFAdd(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN,fidN2;
  int i,np;
  complx *v1,*v2;
  
  if (argc != 3) 
    return TclError(interp,"Usage: fadd <data set 1/result> <data set 2>");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"fadd: argument 1 must be integer <data set>");

  if (Tcl_GetInt(interp,argv[2],&fidN2) == TCL_ERROR) 
    return TclError(interp,"fadd: argument 2 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL) 
    return TclError(interp,"fadd: data set %d was not previously loaded\n",fidN);

  if (fidN2 < 1 || fidN2 > nfd || fd[fidN2] == NULL) 
    return TclError(interp,"fadd: data set %d was not previously loaded\n",fidN2);

  if (!fsamesize(fd[fidN],fd[fidN2])) 
    return TclError(interp,"fadd: %s\n",ferrormsg);

  np=fd[fidN]->np*(fd[fidN]->ni > 1 ? fd[fidN]->ni : 1);
  v1=(complx*)(fd[fidN]->data);
  v2=(complx*)(fd[fidN2]->data);
  for (i=1;i<=np;i++) {
     v1[i].re += v2[i].re;
     v1[i].im += v2[i].im;
  }
  return TCL_OK;
} 

int tclFSub(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN,fidN2;
  int i,np;
  complx *v1,*v2;
  
  if (argc != 3) 
    return TclError(interp,"Usage: fsub <data set 1/result> <data set 2>");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"fsub: argument 1 must be integer <data set>");

  if (Tcl_GetInt(interp,argv[2],&fidN2) == TCL_ERROR) 
    return TclError(interp,"fsub: argument 2 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL) 
    return TclError(interp,"fsub: data set %d was not previously loaded\n",fidN);


  if (fidN2 < 1 || fidN2 > nfd || fd[fidN2] == NULL) 
    return TclError(interp,"fsub: data set %d was not previously loaded\n",fidN2);

  if (!fsamesize(fd[fidN],fd[fidN2])) 
    return TclError(interp,"fsub: %s\n",ferrormsg);

  np=fd[fidN]->np*(fd[fidN]->ni > 1 ? fd[fidN]->ni : 1);
  v1=(complx*)(fd[fidN]->data);
  v2=(complx*)(fd[fidN2]->data);
  for (i=1;i<=np;i++) {
     v1[i].re -= v2[i].re;
     v1[i].im -= v2[i].im;
  }
  return TCL_OK;
} 


typedef struct _RANGE {
  int i1,i2;
} RANGE;


/* Currently implemented with Numerical Recipes. Removed due to Numrep license.

void fbc(complx* src,complx* dst,int np,int order,RANGE* range,int nrange,int skip)
{
  void lfit(double x[], double y[], double sig[], int ndat, double a[], int ia[],
        int ma, double **covar, double *chisq, void (*funcs)(double, double [], int));

  void fpoly(double x, double p[], int np);
  double *vector(long nl, long nh);
  int *ivector(long nl, long nh);
  double **matrix(long nrl, long nrh, long ncl, long nch);        
  void free_vector(double *v, long nl, long nh);
  void free_ivector(int *v, long nl, long nh);
  void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);

  double sum,xx,ff;
  double *x,*y,*sig,*a,chisq,**covar;
  int *ia;
  int n,i,j,i1,i2;

  x=vector(1,np);
  y=vector(1,np);

  n=0;

#define ADDPT(i) {n++; x[n]=(double)(i); y[n]=src[(i)].re;}

  for (i=1;i<=nrange;i++) {
    i1=range[i].i1;
    i2=range[i].i2;
    if (i1 < i2) {
      for (j=i1;j<i2;j += skip) {
        ADDPT(j);
      }
    }
    ADDPT(i2);
  }
  sig=vector(1,n);
  for (i=1;i<=n;i++) {
    sig[i]=1.0;
  }

  a=vector(1,order);
  ia=ivector(1,order);
  covar=matrix(1,order,1,order);
  for (i=1;i<=order;i++) {
    ia[i]=1;
  }

  lfit(x,y,sig,n,a,ia,order,covar,&chisq,fpoly);

  for (i=1;i<=np;i++) {
     sum=0.0;
     xx=1.0;
     sum=a[1];
     ff=(double)i;
     for (j=2;j<=order;j++) {
        xx *= ff;
        sum += a[j]*xx;        
     }
     dst[i].re -= sum;
     dst[i].im = 0.0;
  }
  free_vector(x,1,np);
  free_vector(y,1,np);
  free_vector(sig,1,n);
  free_vector(a,1,order);
  free_ivector(ia,1,order);
  free_matrix(covar,1,order,1,order);
}

int tclFBc(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  FD *f;
  int fidN,np;
  complx *dat;
  int nrange;
  RANGE* range;  
  double v1,v2;
  char **par,**par2;
  int i,i1,i2,npar,npar2;
  int order,skip=1;
  
  if (argc != 4 && argc != 5) 
    return TclError(interp,"Usage: fbc <descr> <order> {{from to} {from to} ..} ?<skip=1>?");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"fbc: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL) 
    return TclError(interp,"fbc: data set %d was not previously loaded\n",fidN);

  if (Tcl_GetInt(interp,argv[2],&order) == TCL_ERROR) 
    return TclError(interp,"fbc: argument 2 must be integer <order of polynomial>");

  if (argc == 5) {
    if (Tcl_GetInt(interp,argv[4],&skip) == TCL_ERROR) 
      return TclError(interp,"fbc: argument 4 must be integer <skip (default=1)>");
  }
  if (Tcl_SplitList(interp,argv[3],&npar,&par) != TCL_OK) 
    return TclError(interp,"fbc: list is not formed correctly\n");

  f=fd[fidN];
  np=f->np;
  dat=(complx*)f->data;
  
  nrange=npar;
  range=(RANGE*)malloc(sizeof(RANGE)*(nrange+1));
  if (!range) 
    return TclError(interp,"fbc: unable to allocate 2 x %d integers\n",nrange);


  for (i=0;i<npar;i++) {
    if (Tcl_SplitList(interp,par[i],&npar2,&par2) != TCL_OK)
      return TclError(interp,"fbc: list element number %d is not formed correctly\n",i+1);

    if (npar2 != 2)
      return TclError(interp,"fbc: list element number %d must contain two values, not %d\n",i+1,npar2);

    if (Tcl_GetDouble(interp,par2[0],&v1) != TCL_OK)
      return TclError(interp,"fbc: unable to convert '%s' to a value in list element number %d\n",par2[0],i+1);

    if (Tcl_GetDouble(interp,par2[1],&v2) != TCL_OK)
      return TclError(interp,"fbc: unable to convert '%s' to a value in list element number %d\n",par2[1],i+1);

    if (v1 >= v2)
      return TclError(interp,"fbc: value 2 must be larger than value 1 in list element number %d\n",i+1);

    i1=FD_INDEX(f,v1);
    if (i1 < 1) i1=1; else if (i1 > np) i1=np;
    i2=FD_INDEX(f,v2);
    if (i2 < 1) i2=1; else if (i2 > np) i2=np;
    range[i+1].i1=i1;
    range[i+1].i2=i2;
    free(par2);
  }
  free(par);

  fbc(dat,dat,np,order,range,nrange,skip);

  free(range);
  return TCL_OK;
}

int tclFSmooth(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  void savgol(double c[], int np, int nl, int nr, int ld, int m);
  void convlv(double data[], unsigned long n, double respns[], unsigned long m, int isign, double ans[]);
  int i,fidN;
  int pnt,ord,pnt2;
  int np;
  complx* dat;
  double *vec,*vec2,*vec3;

  if (argc != 4)
    return TclError(interp,"Usage: fsmooth <data set> <points [16]> <order [4]>");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"fsmooth: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fsmooth: data set %d was not previously loaded\n",fidN);

  if (fd[fidN]->ni > 1)
    return TclError(interp,"fsmooth: cannot operate on 2D data.");

  if (Tcl_GetInt(interp,argv[2],&pnt) == TCL_ERROR) 
    return TclError(interp,"fsmooth: argument 2 must be integer <points [16]>");

  if (Tcl_GetInt(interp,argv[3],&ord) == TCL_ERROR) 
    return TclError(interp,"fsmooth: argument 3 must be integer <order [4]>");


  np=fd[fidN]->np;
  dat=(complx*)(fd[fidN]->data);

  vec=(double*)malloc(sizeof(double)*(np+1));
  vec2=(double*)malloc(sizeof(double)*(np+1));
  vec3=(double*)malloc(sizeof(double)*2*(np+1));
  if (!vec || !vec2 || !vec3)
    return TclError(interp,"fsmooth: unable to allocate 4 times %d doubles\n",np+1);

  pnt2=pnt+pnt+1;
  savgol(vec2,pnt2,pnt,pnt,0,ord);
  for (i=1;i<=np;i++) {
     vec[i]=dat[i].re;
  }
  convlv(vec,np,vec2,pnt2,1,vec3);
  for (i=1;i<=np;i++) {
     dat[i].re=vec3[i];
  }

  savgol(vec2,pnt2,pnt,pnt,0,ord);
  for (i=1;i<=np;i++) {
     vec[i]=dat[i].im;
  }
  convlv(vec,np,vec2,pnt2,1,vec3);
  for (i=1;i<=np;i++) {
     dat[i].im=vec3[i];
  }
  free(vec);
  free(vec2);
  free(vec3);
  return TCL_OK;
} 
*/

int tclFBc(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  return TclError(interp,"fbc: Not implemented yet. Must be rewritten without"
                         " Numerical Recipes dependency (see source code)\n");
}

int tclFSmooth(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  return TclError(interp,"fsmooth: Not implemented yet. Must be rewritten without"
                         " Numerical Recipes dependency (see source code)\n");
}



int tclFSsbint(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN,i,j;
  FD *fn,*f;
  double srate,shift,peakwidth,x;
  
  if (argc != 5) 
    return TclError(interp,"Usage: <data set> fssbint <desc> <srate> <shift> <peakwidth>");

  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR) 
    return TclError(interp,"fssbint: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fssbint: data set %d was not previously loaded\n",fidN);


  f=fd[fidN];
  
  if (Tcl_GetDouble(interp,argv[2],&srate) == TCL_ERROR) 
    return TclError(interp,"fssbint: argument 2 must be double <spinrate>");

  if (Tcl_GetDouble(interp,argv[3],&shift) == TCL_ERROR) 
    return TclError(interp,"fssbint: argument 3 must be double <shift>");

  if (Tcl_GetDouble(interp,argv[4],&peakwidth) == TCL_ERROR) 
    return TclError(interp,"fssbint: argument 4 must be double <peakwidth>");
    
  if (peakwidth != srate ||  shift != 0.0)
    return TclError(interp,"fssbint: this part is not implemented yet. Until then\n"
                           "<peakwidth> must be <spinrate> and <shift> must be 0");

  fn=FD_alloc();
  if (!fn)
     return TclError(interp,"fssbint: unable to allocate fid data set data structure\n");
  
  fn->sw=f->sw;
  fn->np=f->sw/srate;
  fn->type=FD_TYPE_SPE;
  if (!FD_alloc1ddata(fn)) 
    return TclError(interp,"fssbint: unable to allocate %d complex data points\n",fn->np);

  memset(fn->data,0,sizeof(complx)*(fn->np+1));
  
  for (i=1;i<=f->np;i++) {

     x=FD_FREQ(f,i);
     j=FD_INDEX(fn,x);
     /* fprintf(stderr,"x: %g, j: %d, i: %d, newx: %g\n",x,j,i,FD_FREQ(fn,j)); */
     if (j < 1 || j > fn->np) continue;  
     ((complx*)fn->data)[j].re +=((complx*)f->data)[i].re;
     ((complx*)fn->data)[j].im +=((complx*)f->data)[i].im;
  }
  TclSetResult(interp,"%d",fnew(fn));
  return TCL_OK;
} 


void dsave2d_reim(complx* data,int ni,int np,char* name,double scale,double autoscale)
{
  FILE* fp;
  int i,j,k;
  
  fp=fopen(name,"w");
  if (!fp) {
    fprintf(stderr,"error: unable to create file '%s'\n",name);
    exit(1);
  }
  for (j=1,k=1;j<=ni;j++) {
    for (i=1;i<=np;i++,k++) {
      fprintf(fp,"%g %g\n",data[k].re*scale,data[k].im*scale);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}

const char hex_table[] = {'0','1','2','3','4','5','6','7','8','9',
'a','b','c','d','e','f'};

void hex_print(FILE* fp,unsigned char byte)
{ 
  putc(hex_table[byte / 16],fp);
  putc(hex_table[byte % 16],fp);
}

void dsave2d_ps(complx* data,int ni,int np,char* name,double scale,int autoscale)
{
  FILE* fp;
  int i,j,k,l,nmax,npr,nir;
  unsigned char r,g,b;
  double max,min,x;
  
  fp=fopen(name,"w");
  if (!fp) {
    fprintf(stderr,"error: unable to create file '%s'\n",name);
    exit(1);
  }
  nmax=(ni > np ? ni : np);
  npr= (512*np)/nmax;
  nir= (512*ni)/nmax;
  fprintf(fp,
"%%!PS-Adobe-3.0\n"
"%%%%Creator: SIMPSON\n"
"%%%%Title: \n"
"%%%%DocumentData: Clean7Bit\n"
"%%%%Pages: 1\n"
"%%%%BoundingBox: 0 0 %d %d\n"
"%%%%EndComments\n"
"%%%%BeginProlog\n"
"5 dict begin\n"
"%%%%EndProlog\n"
"%%%%Page: 1 1\n"
"%d 0 translate 90 rotate\n"
"%d %d scale\n"
"/scanline %d 3 mul string def\n"
"%d %d 8\n"
"[ %d 0 0 %d 0 0 ]\n"
"{ currentfile scanline readhexstring pop } false 3\n"
"colorimage\n",nir,npr,nir,npr,nir,ni,ni,np,ni,np);

  if (autoscale) {
    max=-1e99;
    min=1e99;
    for (k=1,j=1;j<=ni;j++) {
      for (i=0;i<np;i++) {
        x=data[k++].re;
        if (x > max) max=x; 
        else if (x < min) min=x;
      }
    }
    if (-min > max) max=-min;
    scale *= 65536.0/sqrt(max);
  }
  for (k=1,j=1;j<=np;j++) {

    for (l=0,i=0;i<ni;i++) {
      r=0;
      g=0;
      b=0;
      x= data[i*np+j].re*scale;
      if (x < 0.0) {
         if (x < -65535.0) b=255;
         else b = (unsigned char)(sqrt(-x));
      } else {
         if (x > 65535.0) r=255;
         else r = (unsigned char)(sqrt(x));
      }      
      hex_print(fp,r);
      hex_print(fp,g);
      hex_print(fp,b);
      if ((k++ % 60) == 0) 
        fprintf(fp,"\n");
    }
  }
  fprintf(fp,"\nshowpage\n%%%%Trailer\n"
"end\n"
"%%%%EOF\n");
  fclose(fp);
}

int tclFRead(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  
  if (argc != 3 && argc != 4) {
    return TclError(interp,"usage: fread <handle> <type> ?-swap?\n"
     "\n"
     "type may be short, int, long, float, double, char, bits");
  }
  
  return TCL_OK;
}



void dsave2d_ppm(complx* data,int ni,int np,char* name,double scale,int autoscale)
{
  FILE* fp;
  int i,j,k,l;
  unsigned char *chr,r,g,b;
  double max,min,x;
  
  fp=fopen(name,"w");
  if (!fp) {
    fprintf(stderr,"error: unable to create file '%s'\n",name);
    exit(1);
  }
  chr=(unsigned char*)malloc(ni*3);
  if (!chr) {
    fprintf(stderr,"error: unable to allocate bytes\n");
    exit(1);
  }
  if (autoscale) {
    max=-1e99;
    min=1e99;
    for (k=1,j=1;j<=ni;j++) {
      for (i=0;i<np;i++) {
        x=data[k++].re;
        if (x > max) max=x; 
        else if (x < min) min=x;
      }
    }
    if (-min > max) max=-min;
    scale *= 65536.0/sqrt(max);
  }
  fprintf(fp,"P6\n%d %d\n255\n",ni,np);
  for (j=1;j<=np;j++) {

    for (l=0,i=0;i<ni;i++) {
      r=0;
      g=0;
      b=0;

      x= data[i*np+j].re*scale;
      if (x < 0.0) {
         if (x < -65535.0) b=255;
         else b = (unsigned char)(sqrt(-x));
      } else {
         if (x > 65535.0) r=255;
         else r = (unsigned char)(sqrt(x));
      }      
      chr[l++]=r;
      chr[l++]=g;
      chr[l++]=b;
    }
    fwrite((char*)chr,3,ni,fp);
  }
  free(chr);
  fclose(fp);
}

int tclFPlot2d(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN,autoscale;
  double scale;
  FD* f;
  
  if (argc != 4 && argc != 5) 
    return TclError(interp,"usage: fplot2d <data set> <output file> (-ps | -ppm) ?<scale>?");
  
  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"fplot2d: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fplot2d: data set %d was not previously loaded\n",fidN);
    
  f=fd[fidN];

  if (f->ni <= 0)
    return TclError(interp,"fplot2d: ni must be different from zero when saving 2d data");

  if (argc == 5) {
    autoscale=1;
    if (Tcl_GetDouble(interp,argv[4],&scale) == TCL_ERROR)
      return TclError(interp,"fplot2d: argument 4 must be double <scale>");
  } else {
    autoscale=1;
    scale=1;
  }
  if (!strcmp(argv[3],"-ps")) {
    dsave2d_ps(f->data,f->ni,f->np,argv[2],scale,autoscale);
  } else if (!strcmp(argv[3],"-ppm")) {
    dsave2d_ppm(f->data,f->ni,f->np,argv[2],scale,autoscale);
  } else 
    return TclError(interp,"fplot2d: argument 3 must be either -ps or -ppm");
  
  return TCL_OK;
} 


#define R_HORZ 0
#define R_VERT 1
#define R_DIAG 2
#define R_ANTI 3

int tclFReconstruct(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN, fidN2, i,j, mode=-1,jj,ii;
  FD *f, *f2;
  char buf1[64], buf2[64];
  
  if (argc != 3) 
    return TclError(interp,"usage: <desc> %s <data set> (-horizontal | -diagonal | -antidiagonal | -vertical)", 
      argv[0]);
  
  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"%s: argument 1 must be integer <data set>",
      argv[0]);

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"%s: data set %d was not previously loaded\n",
      argv[0],fidN);
  f=fd[fidN];

  if (f->ni > 1)
    return TclError(interp,"%s: input data must be 1d",
      argv[0]);

  if (!strcmp(argv[2],"-horizontal")) mode=R_HORZ;
  else if (!strcmp(argv[2],"-vertical")) mode=R_VERT;
  else if (!strcmp(argv[2],"-diagonal")) mode=R_DIAG;
  else if (!strcmp(argv[2],"-antidiagonal")) mode=R_ANTI;

  if (mode == -1)
    return TclError(interp,"%s: '%s' is not a valid option. Must be one of -horizontal, -diagonal, -antidiagonal, or -vertical)\n",
      argv[0],fidN);

  f2 = FD_alloc();
  f2->np=f->np;
  f2->ni=f->np;
  f2->sw=f->sw;
  f2->sw1=f->sw;
  f2->ref=f->ref;
  f2->ref1=f->ref;
  f2->type=f->type;
  f2->prec=f->prec;
  FD_alloc1ddata(f2);
  
  switch (mode) {
    case R_ANTI:
      for (j=1; j<= f2->np; j++) {
	for (i=1; i<= f2->np; i++) {
          jj=((i-j+f2->np/2)%f2->np);
          while (jj < 0) jj+= f2->np;
	  while (jj > f2->np) jj-= f2->np;
	  f2->data[jj*f2->np + j].re = f->data[i].re;
	  f2->data[jj*f2->np + j].im = f->data[i].im;
	}
      }
      break;
    case R_DIAG:
      for (j=1; j<= f2->np; j++) {
	for (i=1; i<= f2->np; i++) {
          jj=((j-i+f2->np/2)%f2->np);
          while (jj < 0) jj+= f2->np;
	  while (jj > f2->np) jj-= f2->np;
	  f2->data[jj*f2->np + j].re = f->data[i].re;
	  f2->data[jj*f2->np + j].im = f->data[i].im;
	}
      }
      break;
    case R_VERT:
      for (j=0; j<f2->ni; j++) {
        jj = j*f2->np;
	for (i=1; i<= f2->np; i++) {
	  ii = jj+i;
	  f2->data[ii].re = f->data[i].re;
	  f2->data[ii].im = f->data[i].im;
	}
      }
      break;
    case R_HORZ:
      for (j=0; j<f2->ni; j++) {
	for (i=1; i<= f2->np; i++) {
	  ii = j+(i-1)*f2->np+1;
	  f2->data[ii].re = f->data[i].re;
	  f2->data[ii].im = f->data[i].im;
	}
      }
      break;
  }

  fidN2 = fnew(f2);
  sprintf(buf1, "%d", fidN2);
  sprintf(buf2, "%ld", (long)fd[fidN2]);
  Tcl_Eval(interp,"namespace eval FD_Internal {variable f}");
  Tcl_SetVar2(interp,"FD_Internal::f",buf1,buf2,TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG);
  TclSetResult(interp,"%d",fidN2);
  return TCL_OK;
} 

double frx(FD *f, int i) 
{
  return f->sw*(i-1)/(double)f->np-f->sw/2+f->ref;
}

double fry(FD *f, int j) 
{
  return f->sw1*(j-1)/(double)f->ni-f->sw1/2+f->ref1;
}

int fint(FD *f, double x, double y, complx* v)
{
  int xi, yi;
  double wx, wy;
  complx a, b, c, d;

  xi = (int)(f->np*(x-f->ref+f->sw/2)/f->sw)+1;
  yi = (int)(f->ni*(y-f->ref1+f->sw1/2)/f->sw1)+1;
  
  if (xi < 1 || xi >= f->np || yi < 1 || yi >= f->ni) {
    return 0;
  }
  wx = (x-frx(f,xi))/(frx(f,xi+1)-frx(f,xi));
  wy = (y-fry(f,yi))/(fry(f,yi+1)-fry(f,yi));

  a = f->data[(yi-1)*f->np+xi];
  b = f->data[(yi-1)*f->np+xi+1];
  c = f->data[yi*f->np+xi];
  d = f->data[yi*f->np+xi+1];

  v->re = a.re*(1-wx)*(1-wy) + b.re*wx*(1-wy) + c.re*(1-wx)*wy + d.re*wx*wy;
  v->im = a.im*(1-wx)*(1-wy) + b.im*wx*(1-wy) + c.im*(1-wx)*wy + d.im*wx*wy;
  return 1;
}

int tclFCompare(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN1, fidN2, i, ndata,j;
  FD *f1, *f2;
  complx c;
  
  if (argc != 4) 
    return TclError(interp,"usage: <desc> %s <data set1> <data set2> (-min | -max)", 
      argv[0]);
  
  if (Tcl_GetInt(interp,argv[1],&fidN1) == TCL_ERROR)
    return TclError(interp,"%s: argument 1 must be integer <data set>",
      argv[0]);

  if (fidN1 < 1 || fidN1 > nfd || fd[fidN1] == NULL)
    return TclError(interp,"%s: data set %d was not previously loaded\n",
      argv[0],fidN1);
  f1=fd[fidN1];
  ndata = (f1->ni>1) ? f1->ni*f1->np:f1->np;
  
  if (Tcl_GetInt(interp,argv[2],&fidN2) == TCL_ERROR)
    return TclError(interp,"%s: argument 2 must be integer <data set>",
      argv[0]);

  if (fidN1 < 1 || fidN2 > nfd || fd[fidN2] == NULL)
    return TclError(interp,"%s: data set %d was not previously loaded\n",
      argv[0],fidN2);
  f2=fd[fidN2];
  i = (f1->ni>1) ? f1->ni*f1->np:f1->np;
  
  if (i != ndata)
    return TclError(interp,"%s: data sets must have same number of points\n",
      argv[0]);


  if (f1->sw1 == f2->sw1 && f1->ref1 == f2->ref1) {
    if (!strcmp(argv[3],"-max")) {
      for (i=1; i<=ndata; i++) {
        f1->data[i].re = f1->data[i].re > f2->data[i].re ? f1->data[i].re:f2->data[i].re;
        f1->data[i].im = f1->data[i].im > f2->data[i].im ? f1->data[i].im:f2->data[i].im;
      }
    } else if (!strcmp(argv[3],"-min")) {
      for (i=1; i<=ndata; i++) {
        f1->data[i].re = f1->data[i].re < f2->data[i].re ? f1->data[i].re:f2->data[i].re;
        f1->data[i].im = f1->data[i].im < f2->data[i].im ? f1->data[i].im:f2->data[i].im;
      }
    } else if (!strcmp(argv[3],"-add")) {
      for (i=1; i<=ndata; i++) {
        f1->data[i].re += f2->data[i].re;
        f1->data[i].im += f2->data[i].im;
      }
    }
  } else {
    if (!strcmp(argv[3],"-max")) {
      for (j=1; j<= f1->ni; j++) {
        for (i=1; i<=f1->np; i++) {
          if (fint(f2, FD_FREQ(f1,i), FD_FREQ1(f1,j), &c) == 1) {
  	    ndata = (j-1)*f1->np+i;
	    f1->data[ndata].re = f1->data[ndata].re < c.re ? c.re:f1->data[ndata].re;
	    f1->data[ndata].im = f1->data[ndata].im < c.im ? c.im:f1->data[ndata].im;
	  }
        }
      }
    } else if (!strcmp(argv[3],"-min")) {
      for (j=1; j<= f1->ni; j++) {
        for (i=1; i<=f1->np; i++) {
          if (fint(f2, FD_FREQ(f1,i), FD_FREQ1(f1,j), &c) == 1) {
  	    ndata = (j-1)*f1->np+i;
	    f1->data[ndata].re = f1->data[ndata].re > c.re ? c.re:f1->data[ndata].re;
	    f1->data[ndata].im = f1->data[ndata].im > c.im ? c.im:f1->data[ndata].im;
	  }
        }
      }
    } else if (!strcmp(argv[3],"-add")) {
      for (j=1; j<= f1->ni; j++) {
        for (i=1; i<=f1->np; i++) {
          if (fint(f2, FD_FREQ(f1,i), FD_FREQ1(f1,j), &c) == 1) {
  	    ndata = (j-1)*f1->np+i;
	    f1->data[ndata].re += c.re;
	    f1->data[ndata].im += c.im;
	  }
        }
      }
    }
  }
  return TCL_OK;
} 

int tclFCovariance(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN1, fidN2, i, ndata,j, fidN;
  FD *f1, *f2, *f;
  complx a1, a2, c;
  char buf1[128], buf2[128];
  
  if (argc != 4) 
    return TclError(interp,"usage: <desc> %s <data set1> <data set2> (-h | -v)", 
      argv[0]);
  
  if (Tcl_GetInt(interp,argv[1],&fidN1) == TCL_ERROR)
    return TclError(interp,"%s: argument 1 must be integer <data set>",
      argv[0]);

  if (fidN1 < 1 || fidN1 > nfd || fd[fidN1] == NULL)
    return TclError(interp,"%s: data set %d was not previously loaded\n",
      argv[0],fidN1);
  f1=fd[fidN1];
  ndata = (f1->ni>1) ? f1->ni*f1->np:f1->np;
  
  if (Tcl_GetInt(interp,argv[2],&fidN2) == TCL_ERROR)
    return TclError(interp,"%s: argument 2 must be integer <data set>",
      argv[0]);

  if (fidN1 < 1 || fidN2 > nfd || fd[fidN2] == NULL)
    return TclError(interp,"%s: data set %d was not previously loaded\n",
      argv[0],fidN2);
  f2=fd[fidN2];
  i = (f1->ni>1) ? f1->ni*f1->np:f1->np;
  
  if (i != ndata)
    return TclError(interp,"%s: data sets must have same number of points\n",
      argv[0]);

  if (f1->sw1 != f2->sw1 && f1->ref1 != f2->ref1) {
    return TclError(interp, "%s: data sets must have same sw and ref values\n",
      argv[0]);
  }

  if (!strncmp(argv[3], "-v", 2)) {
    /* vertical projection */
    f=FD_alloc();
    if (!f)
      return TclError(interp,"fcreate: unable to allocate fid data set data structure\n");

    f->np = f1->np;
    f->sw = f1->sw;
    f->ref = f1->ref;
    f->ni = 0;
    if (!FD_alloc1ddata(f))
      return TclError(interp,"fcreate: unable to allocate %d complex data points\n",f->np);
    FD_zero(f);
    for (i=1; i<= f->np; i++) {
      a1.re = a1.im = a2.re = a2.im = c.re = c.im = 0;
      for (j = 1; j<=f1->ni; j++) {
        a1.re += f1->data[(j-1)*f1->np+i].re/f1->ni;
        a1.im += f1->data[(j-1)*f1->np+i].im/f1->ni;
        a2.re += f2->data[(j-1)*f2->np+i].re/f1->ni;
        a2.im += f2->data[(j-1)*f2->np+i].im/f1->ni;
      }
      for (j = 1; j<=f1->ni; j++) {
        c.re += (f1->data[(j-1)*f1->np+i].re-a1.re)*(f2->data[(j-1)*f2->np+i].re-a2.re)/f1->ni;
        c.im += (f1->data[(j-1)*f1->np+i].im-a1.im)*(f2->data[(j-1)*f2->np+i].im-a2.im)/f1->ni;
      }
      f->data[i].re = c.re;
      f->data[i].im = c.im;
    }
  } else {
    /* horizontal projection */
    f->np = f1->ni;
    f->sw = f1->sw1;
    f->ref = f1->ref1;
    f->ni = 0;
    if (!FD_alloc1ddata(f))
      return TclError(interp,"fcreate: unable to allocate %d complex data points\n",f->np);
    for (i=1; i<= f->np; i++) {
      a1.re = a1.im = a2.re = a2.im = c.re = c.im = 0;
      for (j = 1; j<=f1->np; j++) {
        a1.re += f1->data[(i-1)*f1->np+j].re/f1->ni;
        a1.im += f1->data[(i-1)*f1->np+j].im/f1->ni;
        a2.re += f2->data[(i-1)*f2->np+j].re/f1->ni;
        a2.im += f2->data[(i-1)*f2->np+j].im/f1->ni;
      }
      for (j = 1; j<=f1->np; j++) {
        c.re += (f1->data[(i-1)*f1->np+j].re-a1.re)*(f2->data[(i-1)*f2->np+j].re-a2.re)/f1->ni;
        c.im += (f1->data[(i-1)*f1->np+j].im-a1.im)*(f2->data[(i-1)*f2->np+j].im-a2.im)/f1->ni;
      }
      f->data[i].re = c.re;
      f->data[i].im = c.im;
    }
  }

  fidN = fnew(f);
  sprintf(buf1, "%d", fidN);
  sprintf(buf2, "%ld", (long)fd[fidN]);
  Tcl_Eval(interp,"namespace eval FD {variable f}");
  Tcl_SetVar2(interp,"FD::f",buf1,buf2,TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG);
  Tcl_SetVar2(interp,"FD_Internal",buf1,buf2,TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG);
  TclSetResult(interp,"%d",fidN);
  return TCL_OK;
} 



int tclFAbs(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int fidN, i, n;
  FD* f;
  
  if (argc != 2) 
    return TclError(interp,"usage: fabs <data set>");
  
  if (Tcl_GetInt(interp,argv[1],&fidN) == TCL_ERROR)
    return TclError(interp,"fabs: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
    return TclError(interp,"fplot2d: data set %d was not previously loaded\n",fidN);
    
  f=fd[fidN];
  
  n = f->np * (f->ni > 0 ? f->ni:1);
  for (i=1;i<=n;i++) {f->data[i].re = sqrt(f->data[i].re*f->data[i].re+f->data[i].im*f->data[i].im); f->data[i].im = 0;}
  return TCL_OK;
}







void tclcmd_ftools(Tcl_Interp* interp)
{
  Tcl_CreateCommand(interp,"fabs",(Tcl_CmdProc *)tclFAbs,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fload",(Tcl_CmdProc *)tclFLoad,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"floaddata",(Tcl_CmdProc *)tclFLoaddata,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"funload",(Tcl_CmdProc *)tclFUnload,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fsave",(Tcl_CmdProc *)tclFSave,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fphase",(Tcl_CmdProc *)tclFPhase,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fscale",(Tcl_CmdProc *)tclFPhase,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fdup",(Tcl_CmdProc *)tclFDup,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fdupzero",(Tcl_CmdProc *)tclFDupZero,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"frms",(Tcl_CmdProc *)tclFRms,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"frealrms",(Tcl_CmdProc *)tclFRealrms,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fautoscale",(Tcl_CmdProc *)tclFAutoscale,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fextract",(Tcl_CmdProc *)tclFExtract,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fzero",(Tcl_CmdProc *)tclFZero,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fzerofill",(Tcl_CmdProc *)tclFZerofill,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"frev",(Tcl_CmdProc *)tclFRev,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fget",(Tcl_CmdProc *)tclFGet,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fset",(Tcl_CmdProc *)tclFSet,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fcreate",(Tcl_CmdProc *)tclFCreate,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fft",(Tcl_CmdProc *)tclFFt,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"faddpeaks",(Tcl_CmdProc *)tclFAddpeaks,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fint",(Tcl_CmdProc *)tclFInt,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fssbint",(Tcl_CmdProc *)tclFSsbint,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fadd",(Tcl_CmdProc *)tclFAdd,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fsub",(Tcl_CmdProc *)tclFSub,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fcopy",(Tcl_CmdProc *)tclFCopy,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"faddlb",(Tcl_CmdProc *)tclFAddlb,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"findex",(Tcl_CmdProc *)tclFIndex,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fx",(Tcl_CmdProc *)tclFX,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"ffindpeaks",(Tcl_CmdProc *)tclFFindPeaks,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fmaxheight",(Tcl_CmdProc *)tclFMaxheight,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fsetindex",(Tcl_CmdProc *)tclFSetIndex,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fplot2d",(Tcl_CmdProc *)tclFPlot2d,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"faddtriangle",(Tcl_CmdProc *)tclFAddtriangle,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"freconstruct",(Tcl_CmdProc *)tclFReconstruct,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fcompare",(Tcl_CmdProc *)tclFCompare,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

  Tcl_CreateCommand(interp,"fnewnp",(Tcl_CmdProc *)tclFNewnp,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

  Tcl_CreateCommand(interp,"fread",(Tcl_CmdProc *)tclFRead,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

  Tcl_CreateCommand(interp,"fft1d",(Tcl_CmdProc *)tclFFt1d,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"ftranspose",(Tcl_CmdProc *)tclFTranspose,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);


/* Must find replacement for numerical recipes code to work */
  Tcl_CreateCommand(interp,"fbc",(Tcl_CmdProc *)tclFBc,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"fsmooth",(Tcl_CmdProc *)tclFSmooth,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
}


