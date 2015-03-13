/*
    A fast fourier transformation routine translated from some fortran code

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

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif
   

const int POW2[] = {
 1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,
 32768,65536,131072,262144,524288,1048576
};

/* data[1..np]   fourier transform if is=1, invers if is=-1 */

void fft(complx* data,int np,int is)
{
    int ira[19],nr[19],span,step,step2,span2;
    int i,j,j1,j2,kr,ki,jr,j2r,ji,m,n,n2,iff;
    int ir,ifreal,ifimag,irreal,irimag;
    double wreal,wimag,tempr,tempi;
    double dreal,dimag,x;
    double* a;
    complx c;

    if (is > 0) {
      data[1].re *= 0.5;
      data[1].im *= 0.5;
    }
    a = &data[1].re-1;

    m=(int)(log((double)np)/log(2.0)+0.5);
    n=POW2[m];
  
    if (is <= 0) {
      n2=n/2;
      for (i=1,j=n2+1;i<=n2;i++,j++) {   
        c=data[i];
        data[i]=data[j];
        data[j]=c;
      }
      for (i=1;i<=n;i++) {
        data[i].re /= (double)n;
        data[i].im /= (double)n;
      }
    }
  
    n2=n+n-1;
    for (j=1;j<=m;j++) {
      ira[j]=0;
      nr[j]=POW2[j-1];
    }
    iff=1;
L2: ir = ira[m]+1;
    if (ir <= iff) goto L3;
    ifreal=iff+iff-1;
    ifimag=iff+iff;
    irreal=ir+ir-1;
    irimag=ir+ir;
    tempr=a[ifreal];
    tempi=a[ifimag];
    a[ifreal]=a[irreal];
    a[ifimag]=a[irimag];
    a[irreal]=tempr;
    a[irimag]=tempi;
    
L3: iff++;
    if (iff > n) goto L7;
    j=m;
L4: if (ira[j] < nr[j]) goto L5;
    j--;
    goto L4;
L5: ira[j] += nr[j];
L6: if (j == m) goto L2;
    ira[j+1]=ira[j];
    j++;
    goto L6;
L7: for (j1=1;j1<=m;j1++) {
      span=POW2[j1-1];
      step=span+span;
      span2=step;
      step2=step+step;
      wreal=1.0;
      wimag=0.0;
      x=M_PI/(double)span;
      dreal=cos(x);
      dimag=sin(x);
      if (is > 0) dimag = -dimag;
      for (j2=1;j2<=span;j2++) {
        j2r=j2+j2-1;
        for (jr=j2r;jr<=n2;jr += step2) {
          ji=jr+1;
          kr=jr+span2;
          ki=kr+1;
          tempr=a[kr]*wreal-a[ki]*wimag;
          tempi=a[ki]*wreal+a[kr]*wimag;
          a[kr]=a[jr]-tempr;
          a[ki]=a[ji]-tempi;
          a[jr] += tempr;
          a[ji] += tempi;
        }
        tempr=wreal;
        tempi=wimag;
        wreal=tempr*dreal-tempi*dimag;
        wimag=tempi*dreal+tempr*dimag;
      }
    }
  if (is > 0) {
    n2=n/2;
    for (i=1,j=n2+1;i<=n2;i++,j++) {   
      c=data[i];
      data[i]=data[j];
      data[j]=c;
    }
  }
  
  if (is <= 0) {
    data[1].re *= 2;
    data[1].im *= 2;
  }
}
