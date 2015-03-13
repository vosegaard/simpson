#ifndef _OCroutines_H
#define _OCroutines_H

#include "complx.h"
#include "sim.h"

/* ZT: structure to hold and pass OC parameters */
typedef struct _OCoptPars {
	/* all are read only */
	double eps, tol, mnbkstep, cut, stepmin, grad_level, fvals[2];
	int isinit, ndim, method, nIterations, nreport, ncut, max_brack_eval, max_brent_eval, verb, lbfgs_m, ispenalty;
	char VarSaveProc[64];
	int gradmode, gradmodeprop, *var_shapes, *grad_shapes, *var_shapes_penalty_order;
	double *var_shapes_min, *var_shapes_max, *var_shapes_rmsmax, *var_shapes_penalty_factor;
} OCoptPars;
  
/* global variable holding all OC parameters */
extern OCoptPars OCpar;


void store_OCprop(Sim_wsp *wsp);
void store_OCdens(Sim_wsp *wsp);
void set_OCmx_code(Sim_wsp *wsp, char *c);
void incr_OCmx_pos(Sim_wsp *wsp);
void _filterOC(Sim_info *sim, Sim_wsp *wsp,int num);
void _pulse_shapedOC(Sim_info *sim, Sim_wsp *wsp, char *code, int Nch, int Nelem, int *mask, double steptime);
void _pulse_shapedOCprops(Sim_info *sim, Sim_wsp *wsp, char *code, int Nch, int Nelem, int *mask, double steptime);
void _pulse_and_zgrad_shapedOC(Sim_info *sim, Sim_wsp *wsp, char *code, int Nch, int Nelem, int *mask, int zgrs, double steptime);
void _pulse_and_zgrad_shapedOCprops(Sim_info *sim, Sim_wsp *wsp, char *code, int Nch, int Nelem, int *mask, int zgrs, double steptime);
void _pulse_shapedOC_2(Sim_info *sim, Sim_wsp *wsp, int Nelem, int *OCchanmap, int *mask, double duration);
void _pulse_shapedOCprops_2(Sim_info *sim, Sim_wsp *wsp, int Nelem, int *OCchanmap, int *mask, double duration);

void test_print_codes(Sim_wsp *wsp);

#endif
