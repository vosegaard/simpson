/*
    Data format declaration and manipulation
    Copyright (C) 1999 Mads Bak, 2001 Thomas Vosegaard

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

#ifndef __IODATA_H
#define __IODATA_H

#ifdef __cplusplus
extern "C" {
#endif

#include "complx.h"

/*
  Definition of fid:
    1..............np

    dt = 1/sw
    t(1)=0
    t(2)=dt
    t(np) = (np-1)*dt;

  Definition of spectrum:
    1..............np
   freq(1) = -sw/2 + ref
   freq(np) = sw/2 - dx  + ref
   dx = sw/np;

   FD_INDEX returns index
      1 if  frequency is freq(1) + 0.5*dx
      0 if  frequency is freq(1) - 0.5*dx
    that is the index rounds down.
*/
typedef struct _FD {
  complx* data;
  int np,ni,type,type1,format,prec;
  double ref,ref1,sw,sw1;
  double sfrq,sfrq1;
  int rank;
} FD;


#define FREQ(i,np,sw,ref) (sw*(((i)-1)/(double)(np)-0.5)+ref)
#define TIME(i,sw) ((double)(i-1)/(sw))

#define FD_FREQ(fd,i) ((fd)->sw*(((i)-1)/(double)((fd)->np)-0.5)+(fd)->ref)
#define FD_TIME(fd,i) ((double)(i-1)/((fd)->sw))

#define FD_FREQ1(fd,i) ((fd)->sw1*(((i)-1)/(double)((fd)->ni)-0.5)+(fd)->ref1)
#define FD_TIME1(fd,i) ((double)(i-1)/((fd)->sw1))

#define FD_INDEX(fd,freq) ((int)((fd)->np*(((freq)-(fd)->ref)/(fd)->sw+0.5)+1.5))
#define FD_DELTAFREQ(fd) ((fd)->sw/(double)((fd)->np))

#define FD_INDEX1(fd,freq) ((int)((fd)->ni*(((freq)-(fd)->ref1)/(fd)->sw1+0.5)+1.5))
#define FD_DELTAFREQ1(fd) ((fd)->sw1/(double)((fd)->ni))

#define FD_TYPE_FID 0
#define FD_TYPE_SPE 1

#define FD_FORMAT_TEXT 0
#define FD_FORMAT_BINARY 1
#define FD_FORMAT_BZLIB 2

#define FD_PREC_SINGLE 0
#define FD_PREC_DOUBLE 1

void FD_write(char* fname,FD* fd,int format,int prec);
void FD_free(FD* fd);
FD* FD_read(char* fname);
FD* FD_readstr(char* data);
FD* FD_alloc();
FD* FD_dup(FD* fd);
void FD_zero(FD* fd);
int FD_alloc1ddata(FD* fd);

FD* FD_data2fd(char* file,complx* data,int np,int ni,double sw,double sw1);

void float_to_bits(float f,char* bits);
float bits_to_float(char* bits);

void FD_write_rmn(char* fname, FD* fd);
void FD_write_nmrpipe(char* fname, FD* fd, int phsens);
FD* FD_read_nmrpipe(char* fname);
void swap(char *, int);
/****
 * ZT: write in raw binary format, without any header and trailer
 ****/
void FD_write_raw_bin(char* fname,FD* fd);
#ifdef __cplusplus
}
#endif

#endif /* __IODATA_H */
