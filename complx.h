/*
    Complex number
    Copyright (C) 1999 Mads Bak

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

#ifndef __COMPLX_H
#define __COMPLX_H

#ifndef __MATH_H
#include <math.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif
 
typedef struct __complx { double re,im; } complx;

#define Cnorm(z)(z.re*z.re+z.im*z.im)
#define Carg(z) (z.re == 0.0 && z.im == 0.0 ? 0.0 : atan2(z.im,z.re))
#define Czero(c) {(c).re=0.0;(c).im=0.0;}
#define LOG10E 0.434294481903251827651

  
#ifdef __APPLE__
    complx Complx(double re,double im);
    complx Cadd(complx a,complx b);
    complx Cmul(complx a,complx b);
    complx Conj(complx z);
    complx Cexpi(double im);
#else
    inline complx Complx(double re,double im);
    inline complx Cadd(complx a,complx b);
    inline complx Cmul(complx a,complx b);
    inline complx Conj(complx z);
    inline complx Cexpi(double im);
#endif


complx Cneg(complx a);
complx Csub(complx a,complx b);

complx Cdiv(complx a,complx b);
double  Cabs(complx z);
complx Csqrt(complx z);
inline complx RCmul(double x,complx a);
complx CRmul(complx a,double x);
complx Cexp(complx z);
complx CRexp(double re,double im);
complx Cacos(complx z);
complx Csqr(complx z);
complx Clog(complx z);
complx Clog10(complx z);
complx Cpolar(double mag,double angle);
complx CRpow(complx base,double expon);
complx RCpow(double base,complx expon);
complx Cpow(complx base,complx expon);
complx Casin(complx z);
complx Catan(complx z);
complx Ccos(complx z);
complx Ccosh(complx z);
complx Csin(complx z);
complx Csinh(complx z);
complx Ctan(complx z);
complx Ctanh(complx z);

extern complx Cnull;
extern complx Cunit;

#ifdef __cplusplus
}
#endif

#endif


