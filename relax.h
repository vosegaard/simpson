#ifndef _RELAX_H
#define _RELAX_H

#include "matrix.h"
#include "sim.h"

#define MAXRXTERMS 32

typedef struct _SecCouple {
   int q;
   int y;
} SecCouple;

typedef struct _SecQuartet {
   int l1, l2;
   int q;
   int y;
} SecQuartet;
typedef struct _SpecDens {
   int type;
   double *par;
} SpecDens;

typedef struct _RXstruct {
   int Nauto, Ncross;
   SecCouple** secular;
   SecQuartet** ccsecular; 
   mat_complx ***Q, ***Y, ***Ycc;
   double **omega;
   SpecDens *Jcode;
   int *NQ;
   double eqnorm;
} RXstruct;
/* this structure is missing code for cross corelated spectral densities */
/* rows in Q, omega, Y, Ycc and Jcode run from 0 to Nauto-1     */
/* rows in ccsecular run from 0 to Ncross-1                              */

extern RXstruct Relax;

void read_relax(Tcl_Interp *interp, Sim_info* s);
void destroy_Relax(int Nnuc);
double spectral_density(int lam, double om);

void _delay_relax(Sim_info *sim, Sim_wsp *wsp, double dur);

#endif
