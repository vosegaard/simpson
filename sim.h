/*
    Simulation setup and calculation
    Copyright (C) 1999 Mads Bak, Jimmy T. Rasmussen
                  2010 Zdenek Tosner

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

#ifndef __SIM_H
#define __SIM_H

#include <tcl.h>
#include <stdio.h>
#include "cryst.h"
#include "spinsys_struct.h"
#include "matrix.h"
#include "fftw3.h"

typedef struct _Shift {
	int nuc;
	double iso, delta, eta, pas[3];
	mat_double *T;
	complx *Rmol;
	mat_double *T2q[5];
} Shift;

typedef struct _Dipole {
	int nuc[2];
	double delta, eta, pas[3];
	blk_mat_double *blk_T;
	complx *Rmol;
	mat_double *T2q[5];
} Dipole;

typedef struct _Jcoupling {
	int nuc[2];
	double iso, delta, eta, pas[3];
	blk_mat_double *blk_Tiso, *blk_T;
	complx *Rmol;
	mat_double *T2q[5];
} Jcoupling;

typedef struct _Quadrupole {
	int nuc, order;
	double wq, eta, pas[3], w0;
	mat_double *T, *Ta, *Tb, *T3a, *T3b, *T3c;
	complx *Rmol;
	mat_double *T2q[5];
} Quadrupole;

/*****************************
 * type = 1  quadrupole-CSA
 * type = 2  quadrupole-DD
 *****************************/
typedef struct _Mixing {
	int type, couple[2], qidx, idx;
	mat_double *T;
	blk_mat_double *Ta, *Tb;
} Mixing;

//typedef struct _Perms {
//	int *permvec;
//	int *dims;
//} Perms;

typedef struct _Sim_info {
  /* remains from simpson 2 */
  SpinSys* ss;
  mat_complx *fstart, *fdetect;
  /* spinsys interactions */
  Shift **CS;
  Dipole **DD;
  Jcoupling **J;
  Quadrupole **Q;
  Mixing **MIX;
  int nCS, nDD, nJ, nQ, nMIX;
  double specfreq, wr, sw, sw1, brl, gamma_zero;
  int np, ni, ntot, ngamma, matdim, obs_nuc, imethod, crystfile_from, crystfile_to;
  char crystfile[256],rfproffile[256],pulsename[64],parname[256],targetcrystfile[256];
  /* switches */
  int acq_adjoint, dipole_check, realspec, sparse, conjugate_fid, EDsymmetry, Hint_isdiag;
  int block_diag, averaging, interpolation;
  /* DOR related */
  int dor;
  double brl1, brl2, wr1, wr2;
  /* ZT add-ons */
  int relax, propmethod, domain, points_per_cycle;
  /* ZT: tricks for new direct method */
  double dw1, taur;
  /* Hamiltonian block structure info */
  int basis, Nbasis, Nmz;
  int **perm_table, **dims_table;
  /* new Hamiltonian assembling */
  int Hassembly;
  blk_mat_double *Hiso, *HQ[5];
  Cryst *crdata, *targetcrdata;
  double **rfdata;
  int Jinterpol[2], *crmap;
  TRIANGLE *tridata;
  double **ASG_freq, ASG_period;
  complx **ASG_ampl, *FWT_lam, **FWT_frs;
  int *FWTASG_nnz, **FWT_irow, **FWT_icol;
  fftw_plan *fftw_plans;
  // testing to take average value instead of initial point Hamiltonian
  int do_avg;
  int labframe;
} Sim_info;

#define MAXSTO 1000
#define ACQBLOCK_STO_INI 1001
#define ACQBLOCK_STO_END 15000
#define MAXOCPROPS 5120
#define MAXOCSHAPES 4

/* all times are stored in microseconds */
typedef struct _Sim_wsp {
	double zcoor, tstart, t, dtmax, tpropstart, brl, gamma_add, dt_gcompute;
	double *phv, *rfv, *rffac, *inhom_offset, *offset;
	Cryst cryst;
	blk_mat_double *ham_blk, *sumHrf;
	mat_double *sumUph;
	mat_double **chan_Ix, **chan_Iy, **chan_Iz, **Ix, **Iz;
	mat_complx *fstart, *fdetect, *sigma;
	blk_mat_complx *U, *dU, *tmpU;
	complx * fid;
	int curr_nsig;
    complx **CS_Rrot, **CS_Rlab;
    complx **DD_Rrot, **DD_Rlab;
    complx **Q_Rrot, **Q_Rlab;
    complx **J_Rrot, **J_Rlab;
    int *CS_used, *DD_used, *Q_used, *J_used, *spinused;
    blk_mat_complx *STO[ACQBLOCK_STO_END];
    mat_complx *matrix[ACQBLOCK_STO_END];
    double   STO_tproplength_usec[ACQBLOCK_STO_END], STO_tpropstart_usec[ACQBLOCK_STO_END], STO_brl[ACQBLOCK_STO_END];
    int STO_waspulse[ACQBLOCK_STO_END];
    int isselectivepulse,cannotbestored,hamchanged,Uisunit,waspulse;
    Tcl_Interp *interp;
    /* parameters for acq_block feature */
    int evalmode, Nacq, acqblock_sto;
    double acqphase, acqblock_t0;
    /* parameters for Optimal Control */
    int OC_propstatus, OC_mxpos;
    mat_complx *OC_dens[MAXOCPROPS], *OC_deriv[MAXOCPROPS][MAXOCSHAPES*2];
    blk_mat_complx *OC_props[MAXOCPROPS];
    char *OC_mxcode[MAXOCPROPS];
    /* parametrs for spinach */
    int spinach_delay_Nsub, **spinach_delay_sub_idx, spinach_pulse_Nsub, **spinach_pulse_sub_idx;
    /* new Hamiltonian assembling */
    mat_complx *Dmol_rot;
    blk_mat_double *Hiso, *Hiso_off, *HQ[5], *HQ_off[5];
    int Nint_off;
    mat_double **QTa, **QTb, **QT3a, **QT3b, **QT3c, **MT;
    blk_mat_double **MTa, **MTb;
    /* when averaging_file was added */
    Shift **CS;
    Dipole **DD;
    Jcoupling **J;
    Quadrupole **Q;
    int thread_id, cryst_idx, ig;
    int *FWTASG_irow, *FWTASG_icol;
    double dw;
    // labframe
    mat_complx *Hlab, *Hrflab;
} Sim_wsp;



Sim_info * sim_initialize(Tcl_Interp* interp);
Sim_info * sim_duplicate(Sim_info *sim);
Sim_wsp * wsp_initialize(Sim_info * sim);
void sim_destroy(Sim_info* sim, int this_is_copy);
void wsp_destroy(Sim_info* sim, Sim_wsp * wsp);
int sim_calcfid(Sim_info* sim,Sim_wsp * wsp);
int sim_calcfid_interpol(Sim_info *s, Sim_wsp *wsp);
void store_sim_pointers(Tcl_Interp* interp, Sim_info* sim, Sim_wsp * wsp);
void read_sim_pointers(Tcl_Interp* interp, Sim_info **sim, Sim_wsp **wsp);
void sim_prepare_interpol(Sim_info *sim);
void sim_preempty_interpol(Sim_info *sim);

// important: _FREQ must be _TIME + 1
#define    M_GCOMPUTE_TIME   2001
#define    M_GCOMPUTE_FREQ   2002
#define    M_DIRECT_TIME     2003
#define    M_DIRECT_FREQ     2004
#define    M_SPINACH         2005

#define    EM_NORMAL    101
#define    EM_ACQBLOCK  201
#define    EM_MEASURE   301

#define    INTERPOL_NOT_USED    0
#define    INTERPOL_FWT_ALL     1
#define    INTERPOL_FWT_LAM     2
#define    INTERPOL_ASG         3
#define    INTERPOL_FWTASG_ALL  4
#define    INTERPOL_FWTASG_LAM  5

#endif /* __SIM_H */
