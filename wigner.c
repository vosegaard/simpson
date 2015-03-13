/*
    Wigner rotation routines
    Copyright (C) 1999 Mads Bak
    
    2009: Z.T. changes for new structures of matrices

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
    
    Creating of Wigner roatation matrices -- either the full matrix
    or the zero'th column.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complx.h"
#include "defs.h"
#include "matrix.h"
#include "cm.h"

#ifdef INTEL_MKL
#include "mkl.h"
#elif defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#elif defined(GSL)
#include <gsl/gsl_cblas.h>
#else
#include "cblas.h"
#endif

complx* Dtensor2(double deltazz,double eta)
{
  complx* v;

  v=complx_vector(5);
  v[1]=v[5]= Complx(-eta*deltazz/2.0,0.0);
  v[2]=v[4]= Complx(0.0,0.0);
  v[3]=      Complx(deltazz*sqrt(3.0/2.0),0.0);
  return v;
}


mat_complx * wigner2(double alpha,double beta,double gamma)
{
  double cosb,sinb,cos2b,sin2b,cplus,cminus,gamma2,alpha2;
  double cplus2,cminus2,sinof2b,cplussinb,cminussinb,SQRT3BY8sinof2b,SQRT3BY8sin2b;
  complx em2am2g,em2amg,em2a,em2apg,em2ap2g,emam2g,emamg,ema;
  complx emapg,emap2g,em2g,emg,epg,ep2g,epam2g,epamg,epa;    
  complx epapg,epap2g,ep2am2g,ep2amg,ep2a,ep2apg,ep2ap2g;
  mat_complx * res;
  
  alpha *= DEG2RAD;
  beta *= DEG2RAD;
  gamma *= DEG2RAD;
  cosb=cos(beta);
  sinb=sin(beta);
  cos2b=cosb*cosb;
  sin2b=sinb*sinb;
  cplus=(1.0+cosb)*0.5;
  cminus=(1.0-cosb)*0.5;
  cplus2=cplus*cplus;
  cminus2=cminus*cminus;
  sinof2b=sin(2.0*beta);

  alpha2=-2.0*alpha;
  gamma2=2.0*gamma;
  
  cplussinb=cplus*sinb;
  cminussinb=cminus*sinb;
  SQRT3BY8sinof2b=SQRT3BY8*sinof2b;
  SQRT3BY8sin2b=SQRT3BY8*sin2b;
  
  epa.re = cos(alpha); epa.im = sin(alpha);
  ema.re = epa.re; ema.im = -epa.im;
  ep2a.re = epa.re*epa.re - epa.im*epa.im;
  ep2a.im = epa.re*epa.im*2;
  em2a.re = ep2a.re; em2a.im = -ep2a.im;

  epg.re = cos(gamma); epg.im = sin(gamma);
  emg.re = epg.re; emg.im = -epg.im;
  ep2g.re = epg.re*epg.re - epg.im*epg.im;
  ep2g.im = epg.re*epg.im*2;
  em2g.re = ep2g.re; em2g.im = -ep2g.im;


  em2am2g.re = em2a.re*em2g.re - em2a.im*em2g.im; em2am2g.im = em2a.im*em2g.re + em2a.re*em2g.im;
  em2amg.re = em2a.re*emg.re - em2a.im*emg.im; em2amg.im = em2a.im*emg.re + em2a.re*emg.im;
  //em2a
  em2apg.re = em2a.re*epg.re - em2a.im*epg.im; em2apg.im = em2a.im*epg.re + em2a.re*epg.im;
  em2ap2g.re = em2a.re*ep2g.re - em2a.im*ep2g.im; em2ap2g.im = em2a.im*ep2g.re + em2a.re*ep2g.im;

  emam2g.re = ema.re*em2g.re - ema.im*em2g.im; emam2g.im = ema.im*em2g.re + ema.re*em2g.im;
  emamg.re  = ema.re*emg.re - ema.im*emg.im; emamg.im = ema.im*emg.re + ema.re*emg.im;
  //ema
  emapg.re  = ema.re*epg.re - ema.im*epg.im; emapg.im = ema.im*epg.re + ema.re*epg.im;
  emap2g.re = ema.re*ep2g.re - ema.im*ep2g.im; emap2g.im = ema.im*ep2g.re + ema.re*ep2g.im;

  //em2g
  //emg
  //epg
  //ep2g

  epam2g.re = emap2g.re; epam2g.im = -emap2g.im;
  epamg.re  = emapg.re;  epamg.im  = -emapg.im;
  //epa
  epapg.re  = emamg.re; epapg.im  = -emamg.im;
  epap2g.re = emam2g.re; epap2g.im = -emam2g.im;

  ep2am2g.re = em2ap2g.re; ep2am2g.im = -em2ap2g.im;
  ep2amg.re = em2apg.re; ep2amg.im = -em2apg.im;
  //ep2a
  ep2apg.re = em2amg.re; ep2apg.im = -em2amg.im;
  ep2ap2g.re = em2am2g.re; ep2ap2g.im = -em2am2g.im;

  res = complx_matrix(5,5,MAT_DENSE,0,0);
  complx *d2 = res->data;
  /* first column D_i-2 */
  d2[ 0].re = cplus2*ep2ap2g.re; d2[ 0].im = cplus2*ep2ap2g.im;
  d2[ 1].re = -cplussinb*epap2g.re; d2[ 1].im = -cplussinb*epap2g.im;
  d2[ 2].re = SQRT3BY8sin2b*ep2g.re; d2[ 2].im = SQRT3BY8sin2b*ep2g.im;
  d2[ 3].re = -cminussinb*emap2g.re; d2[ 3].im = -cminussinb*emap2g.im;
  d2[ 4].re = cminus2*em2ap2g.re; d2[ 4].im = cminus2*em2ap2g.im;
  /* second column D_i-1 */
  d2[ 5].re = cplussinb*ep2apg.re; d2[ 5].im = cplussinb*ep2apg.im;
  d2[ 6].re = (cos2b-cminus)*epapg.re; d2[ 6].im = (cos2b-cminus)*epapg.im;
  d2[ 7].re = -SQRT3BY8sinof2b*epg.re; d2[ 7].im = -SQRT3BY8sinof2b*epg.im;
  d2[ 8].re = (cplus-cos2b)*emapg.re; d2[ 8].im = (cplus-cos2b)*emapg.im;
  d2[ 9].re = -cminussinb*em2apg.re; d2[ 9].im = -cminussinb*em2apg.im;
  /* third column D_i0 */
  d2[10].re = SQRT3BY8sin2b *ep2a.re; d2[10].im = SQRT3BY8sin2b *ep2a.im;
  d2[11].re = SQRT3BY8sinof2b*epa.re; d2[11].im = SQRT3BY8sinof2b*epa.im;
  d2[12].re = 0.5*(3.0*cos2b-1.0); d2[12].im = 0.0;
  d2[13].re = -SQRT3BY8sinof2b*ema.re; d2[13].im = -SQRT3BY8sinof2b*ema.im;
  d2[14].re = SQRT3BY8sin2b*em2a.re; d2[14].im = SQRT3BY8sin2b*em2a.im;
  /* fourth column D_i+1 */
  d2[15].re = cminussinb*ep2amg.re; d2[15].im = cminussinb*ep2amg.im;
  d2[16].re = (cplus-cos2b)*epamg.re; d2[16].im = (cplus-cos2b)*epamg.im;
  d2[17].re = SQRT3BY8sinof2b*emg.re; d2[17].im = SQRT3BY8sinof2b*emg.im;
  d2[18].re = (cos2b-cminus)*emamg.re; d2[18].im = (cos2b-cminus)*emamg.im;
  d2[19].re = -cplussinb*em2amg.re; d2[19].im = -cplussinb*em2amg.im;
  /* fifth column D_i+2 */
  d2[20].re = cminus2*ep2am2g.re; d2[20].im = cminus2*ep2am2g.im;
  d2[21].re = cminussinb*epam2g.re; d2[21].im = cminussinb*epam2g.im;
  d2[22].re = SQRT3BY8sin2b*em2g.re; d2[22].im = SQRT3BY8sin2b*em2g.im;
  d2[23].re = cplussinb*emam2g.re; d2[23].im = cplussinb*emam2g.im;
  d2[24].re = cplus2*em2am2g.re; d2[24].im = cplus2*em2am2g.im;

  return res;
}


complx * wigner20(double alpha,double beta)
{
  double cosb,sinb,Ksin2b,Ksinof2b;
  complx em2a,ema,epa,ep2a;
  complx *res;

  alpha *= DEG2RAD;
  beta *= DEG2RAD;

  cosb   =cos(beta);
  sinb   =sin(beta);
  Ksin2b  =sinb*sinb*SQRT3BY8;
  Ksinof2b=sin(2.0*beta)*SQRT3BY8;

  epa.re =cos(alpha); epa.im = sin(alpha);
  ema.re = epa.re; ema.im = -epa.im;
  ep2a.re = epa.re*epa.re - epa.im*epa.im; ep2a.im = epa.re*epa.im*2.0;
  em2a.re = ep2a.re; em2a.im = -ep2a.im;

  res = complx_vector(5);
  res[5].re = Ksin2b*em2a.re; res[5].im = Ksin2b*em2a.im;
  res[4].re = -Ksinof2b*ema.re; res[4].im = -Ksinof2b*ema.im;
  res[3].re = 0.5*(3.0*cosb*cosb-1.0); res[3].im = 0.0;
  res[2].re =  Ksinof2b*epa.re; res[2].im =  Ksinof2b*epa.im;
  res[1].re =  Ksin2b*ep2a.re; res[1].im =  Ksin2b*ep2a.im;
  return res;
}

void wig2rot(complx* res, complx* vec, mat_complx *d2)
{
   const int N=5;
   complx zone = {1.0,0.0};

   //cblas_zgemv(CblasColMajor,CblasTrans,N,N,&zone,d2->data,N,&vec[1],1,&Cnull,&res[1],1);
   cblas_zgemv(CblasColMajor,CblasTrans,N,N,&zone,d2->data,N,vec+1,1,&Cnull,res+1,1);
}

void wig2rot_t(complx* res, complx* vec, mat_complx *d2)
{
   const int N=5;
   complx zone = {1.0,0.0};

   cblas_zgemv(CblasColMajor,CblasNoTrans,N,N,&zone,d2->data,N,vec+1,1,&Cnull,res+1,1);
}

double wig20rot(complx *vec, complx *d20)
{
   complx res;
   const int len=5;
   //cblas_zdotu_sub(len,&vec[1],1,&d20[1],1,&res);
   cblas_zdotu_sub(len,vec+1,1,d20+1,1,&res);
   if ( fabs(res.im) > 1e-8 ) {
      fprintf(stderr,"wig20rot error: result is not pure real\n");
      exit(1);
   }
   return res.re; 
}

complx * wig2roti(complx *vec, double a, double b, double c)
{
	mat_complx *w;
	complx *v;

	w = wigner2(a,b,c);
	v = complx_vector(5);
	wig2rot(v,vec,w);
	free_complx_matrix(w);
	free_complx_vector(vec);
	return v;
}
