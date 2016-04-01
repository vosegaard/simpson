/*
 * spinsys_struct.h
 *
 *  Created on: Jun 24, 2010
 *      Author: zdenek
 */

#ifndef SPINSYS_STRUCT_H_
#define SPINSYS_STRUCT_H_

#define MAXCHAN  31
#define MAXSPINS 127
#define MAXCHANSPINS 127

typedef struct _ISOTOPE {
  int  number;
  char name[8];
  double spin;
  double gamma;
} ISOTOPE;


typedef struct _SpinSys {
  ISOTOPE* iso[MAXSPINS+1];
  int nspins;
  int matdim;

  int chan[MAXCHAN+1][MAXCHANSPINS+1];
  int nchanelem[MAXCHAN+1];
  char channames[MAXCHAN+1][8];
  int nchan;
} SpinSys;

#endif /* SPINSYS_STRUCT_H_ */
