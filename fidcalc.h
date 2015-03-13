/*
 * fidcalc.h
 *
 *  Created on: Jun 24, 2010
 *      Author: zdenek
 */

#ifndef FIDCALC_H_
#define FIDCALC_H_

#include "tcl.h"
#include "sim.h"

typedef struct _acq_data {
	complx *fid, *lams, *frs;
	int *irow, *icol;
	double cosph, sinph, weight, taur, wr, sw;
	int Nacq, Ng, Np, m, dim, isphase, ncp;
} acq_data;

int is_rhosymmetry(mat_complx *fstart, mat_complx *fdetect);
void acqblock_disable_command(Sim_wsp *wsp, char *cmd);
void acqblock_sto_incr(Sim_wsp *wsp);
void direct_acqblock(Tcl_Interp *interp,Tcl_Obj *obj,Sim_info *sim,Sim_wsp *wsp);
void direct_acqblock_freq(Tcl_Interp *interp,Tcl_Obj *obj,Sim_info *sim,Sim_wsp *wsp);
void direct_acqblock_time_FWT(Tcl_Interp *interp,Tcl_Obj *obj,Sim_info *sim,Sim_wsp *wsp);
void gcompute_acqblock(Tcl_Interp *interp,Tcl_Obj *obj,Sim_info *sim,Sim_wsp *wsp);

//void new_gcompute(Sim_info *s, Sim_wsp *wsp);
//void new_gcompute_time(Sim_info *s, Sim_wsp *wsp);
//void new_gcompute_freq(Sim_info *s, Sim_wsp *wsp);
//void collect_fid_interpol_direct(int icr, Sim_info *sim, complx *fid);
//void collect_spc_interpol_direct(int icr, Sim_info *sim, complx *fid);

void gcompute_fid(Sim_info *sim, Sim_wsp *wsp);
void gcompute_spc(Sim_info *sim, Sim_wsp *wsp);
void gcompute_ASGdata(Sim_info *sim, Sim_wsp *wsp);
void gcompute_FWTdata(Sim_info *sim, Sim_wsp *wsp);
void collect_fid_interpol_all(int icr, Sim_info *sim, complx *fid);
void collect_fid_interpol_lam(int icr, Sim_info *sim, complx *fid);
void collect_spc_interpol_all(int icr, Sim_info *sim, complx *fid, int thrd_id);
void collect_spc_interpol_lam(int icr, Sim_info *sim, complx *fid, int thrd_id);
void convert_FWTtoASG_gcompute(Sim_info *sim, int icr, int thrd_id);

void direct_acqblock_ASGdata(Tcl_Interp *interp,Tcl_Obj *obj,Sim_info *sim,Sim_wsp *wsp);
void direct_acqblock_FWTdata(Tcl_Interp *interp,Tcl_Obj *obj,Sim_info *sim,Sim_wsp *wsp);
void collect_spc_direct_interpol_all(int icr, Sim_info *sim, complx *fid, int thrd_id);
void collect_spc_direct_interpol_lam(int icr, Sim_info *sim, complx *fid, int thrd_id);
void convert_FWTtoASG_direct(Sim_info *sim, int icr, int thrd_id);


#endif /* FIDCALC_H_ */
