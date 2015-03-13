/*
    Complex numbers
    Copyright (C) 1999 Mads Bak, Jimmy T. Rasmussen

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
    
*/

#include <math.h>
#include "complx.h"

complx Cnull={0.0,0.0};
complx Cunit={1.0,0.0};

inline complx Complx(double re,double im)
{
  complx c;

  c.re=re;
  c.im=im;
  return c;
}

inline complx Cmul(complx a,complx b)
{
  complx c;

  c.re=a.re*b.re-a.im*b.im;
  c.im=a.im*b.re+a.re*b.im;
  return c;
}

inline complx Cadd(complx a,complx b)
{
  complx c;
  c.re=a.re+b.re;
  c.im=a.im+b.im;
  return c;
}

inline complx Conj(complx z)
{
  z.im= -z.im;
  return z;
}

complx Cneg(complx a)
{
  complx c;
  c.re= -a.re;
  c.im= -a.im;
  return c;
}

complx Csub(complx a,complx b)
{
  complx c;

  c.re=a.re-b.re;
  c.im=a.im-b.im;
  return c;
}




complx Cdiv(complx a,complx b)
{
   complx c;
   double r,den;

   if (fabs(b.re) >= fabs(b.im)) {
      r=b.im/b.re;
      den=b.re+r*b.im;
      c.re=(a.re+r*a.im)/den;
      c.im=(a.im-r*a.re)/den;
   } else {
      r=b.re/b.im;
      den=b.im+r*b.re;
      c.re=(a.re*r+a.im)/den;
      c.im=(a.im*r-a.re)/den;
   }
   return c;
}

double Cabs(complx z)
{
  double x,y,ans,temp;
  x=fabs(z.re);
  y=fabs(z.im);
  if (x == 0.0)
    ans=y;
  else if (y == 0.0)
    ans=x;
  else if (x > y) {
    temp=y/x;
    ans=x*sqrt(1.0+temp*temp);
  } else {
    temp=x/y;
    ans=y*sqrt(1.0+temp*temp);
  }  
  return ans;
}

complx Csqrt(complx z)
{
  complx c;
  double x,y,w,r;

  if ((z.re == 0.0) && (z.im == 0.0)) 
    return z;
  else {
    x=fabs(z.re);
    y=fabs(z.im);
    if (x >= 0.0) {
      r=y/x;
      w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
    } else {
      r=x/y;
      w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));    
    }
    if (z.re >= 0.0) {
      c.re=w;
      c.im=z.im/(2.0*w);
    } else {
      c.im=(z.im >= 0.0) ? w : -w;
      c.re=z.im/(2.0*c.im);
    }
    return c;
  }
}

complx RCmul(double x,complx a)
{
  a.re *= x;
  a.im *= x;
  return a;
}

inline complx CRmul(complx a,double x)
{
  a.re *= x;
  a.im *= x;
  return a;
}

complx Cexp(complx z)
{
  double f;

  f=exp(z.re);
  z.re=f*cos(z.im);
  z.im=f*sin(z.im);
  return z;
}

inline complx Cexpi(double im)
{
  complx z;
  z.re=cos(im);
  z.im=sin(im);
  return z;
}

complx CRexp(double re,double im)
{
  complx z;
  double f;

  f=exp(re);
  z.re=f*cos(im);
  z.im=f*sin(im);
  return z;
}


complx Cacos(complx z)
{
  complx c;
  double phi,rp;

  c.re = 1-z.re*z.re+z.im*z.im;
  c.im = -2*z.re*z.im;
  phi = Carg(c)/2.0;
  rp = sqrt(Cabs(c));
  c.re=z.re-rp*sin(phi);
  c.im=z.im+rp*cos(phi);
  rp=Cabs(c);
  c.re=Carg(c);
  c.im= -log(rp);
  return c;
}


complx Csqr(complx z)
{
  complx c;

  c.re=z.re*z.re-z.im*z.im;
  c.im=2.0*z.re*z.im;
  return c;
}


complx Clog(complx z)
{
  complx c;

  c.re=log(Cabs(z));
  c.im=Carg(z);
  return c;
}


complx Clog10(complx z)
{
  complx c;

  c.re=log10(Cabs(z));
  c.im=Carg(z)*LOG10E;
  return c;
}


complx Cpolar(double mag,double angle)
{
  complx c;

  c.re=mag*cos(angle);
  c.im=mag*sin(angle);
  return c;
}


complx CRpow(complx base,double expon)
{
  double mag,angle,ans,x,y,temp;
  complx c;

  if (base.im == 0.0 && base.re == 0.0 && expon > 0.0) return Cnull;
  x=fabs(base.re);
  y=fabs(base.im);
  if (x > y) {
    temp=y/x;
    ans=x*sqrt(1.0+temp*temp);
  } else {
    temp=x/y;
    ans=y*sqrt(1.0+temp*temp);
  }
  mag=pow(ans,expon);
  angle=expon*Carg(base);
  c.re=mag*cos(angle);
  c.im=mag*sin(angle);
  return c;
}

complx RCpow(double base,complx expon)
{
  complx c;
  double lnx,f;

  if (base == 0.0 && expon.re > 0.0) return Cnull;
  lnx=log(fabs(base));
  if (base > 0.0) {
    c.re=expon.re * lnx;
    c.im=expon.im * lnx;
  } else {
    c.re=expon.re*lnx-expon.im*M_PI;
    c.im=expon.im*lnx+expon.re*M_PI;
  }
  f=exp(c.re);
  c.re=f*cos(c.im);
  c.im=f*sin(c.im);
  return c;
}

complx Cpow(complx base,complx expon)
{
  complx c;
  double x,y,temp,ans,f;

  if (base.re == 0.0 && base.im == 0.0 && expon.re > 0.0) return Cnull;
  x=fabs(base.re);
  y=fabs(base.im);
  if (x > y) {
    temp=y/x;
    ans=x*sqrt(1.0+temp*temp);
  } else {
    temp=x/y;
    ans=y*sqrt(1.0+temp*temp);
  }
  c.re=log(ans);
  c.im=Carg(base);
  base.re=expon.re*c.re-expon.im*c.im;
  base.im=expon.im*c.re+expon.re*c.im;
  f=exp(base.re);
  c.re=f*cos(base.im);
  c.im=f*sin(base.im);
  return c;
}

complx Casin(complx z)
{
  complx c;
  double phi,rp;

  c.re = 1-z.re*z.re+z.im*z.im;
  c.im = -2*z.re*z.im;
  phi = Carg(c)/2.0;
  rp = sqrt(Cabs(c));
  c.re= -z.im+rp*cos(phi);
  c.im= z.re+rp*sin(phi);
  rp=Cabs(c);
  c.re=Carg(c);
  c.im= -log(rp);
  return c;
}

complx Catan(complx z)
{
  double opb,a2,den;
  complx c;

  opb=1+z.im;
  a2=z.re*z.re;
  den=opb*opb+a2;
  c.re=((1-z.im)*opb-a2)/den;
  c.im= 2*z.re/den;
  den=Cabs(c);
  c.re=Carg(c)/2.0;
  c.im= -log(den)/2.0;
  return c;
}

complx Ccos(complx z)
{
  complx c;
  double eb,emb;

  eb=exp(z.im);
  emb=1.0/eb;
  c.re=cos(z.re)*(emb+eb)/2.0;
  c.im=sin(z.re)*(emb-eb)/2.0;
  return c;
}

complx Ccosh(complx z)
{
  complx c;
  double ea,eainv;

  ea=exp(z.re);
  eainv=1.0/ea;
  c.re=cos(z.im)*(ea + eainv)/2.0;
  c.im=sin(z.im)*(ea - eainv)/2.0;
  return c;
}

complx Csin(complx z)
{
  complx c;
  double eb,emb;

  eb=exp(z.im);
  emb=1.0/eb;
  c.re=sin(z.re)*(emb+eb)/2.0;
  c.im= -0.5*cos(z.re)*(emb-eb);
  return c;
}

complx Csinh(complx z)
{
  complx c;
  double ea,eainv;

  ea=exp(z.re);
  eainv=1.0/ea;
  c.re=cos(z.im)*(ea-eainv)/2.0;
  c.im=sin(z.im)*(ea+eainv)/2.0;
  return c;
}


complx Ctan(complx z)
{
  complx c;
  double sina,cosa,emin,eplus,temp,temp2;
  double eb,emb;

  sina=sin(z.re);
  cosa=cos(z.re);
  eb=exp(z.im);
  emb=1.0/eb;
  emin=emb-eb;
  eplus=emb+eb;
  temp=cosa*eplus;
  temp2=sina*emin;
  temp=temp*temp+temp2*temp2;
  c.re= 4*sina*cosa/temp;
  c.im= -emin*eplus/temp;
  return c;
}

complx Ctanh(complx z)
{
  complx c;
  double sinb,cosb,eamin,eapls,temp,temp2;
  double ea,eainv;

  sinb = sin(z.im);
  cosb = cos(z.im);
  ea = exp(z.re);
  eainv = 1 / ea;
  eamin = ea - eainv;
  eapls = ea + eainv;
  temp=cosb*eapls;
  temp2=sinb*eamin;
  temp=temp*temp+temp2*temp2;
  c.re= eamin*eapls/temp;
  c.im= 4*sinb*cosb/temp;
  return c;
}
