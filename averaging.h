/*
 * averaging.h
 *
 *  Created on: 14.9.2011
 *      Author: Zdenek Tosner
 */

#ifndef AVERAGING_H_
#define AVERAGING_H_

#include <tcl.h>
#include "sim.h"

typedef struct _Ave_elem {
	int par_type, nuc1, nuc2, val_type;
	char *txt;
	double *data;
} Ave_elem;

/*
 * par_type key-code:
 * ------------------
 *  0 ... par(*)
 *  1 ... shift
 *  2 ... dipole
 *  3 ... jcoupling
 *  4 ... quadrupole
 *
 *  val_type key-code:
 *  ------------------
 *  0 ... not specified
 *  1 ... iso
 *  2 ... aniso
 *  3 ... eta
 *  4 ... alpha
 *  5 ... beta
 *  6 ... gamma
 */


void read_averaging_file(Tcl_Interp *interp, Sim_info *sim, Ave_elem **aveptr, double **weightprt, int *Npar, int *Nval);
void free_averaging_data(Ave_elem **avestruct,int Navepar);
void set_averaging_parameters(Sim_info *sim,Sim_wsp *wsp,Ave_elem *ave_struct,int Navepar,int ave_start);

#endif /* AVERAGING_H_ */
