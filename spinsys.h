/*
    Spinsystem declaration and routines
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

#ifndef __SPINSYS_H
#define __SPINSYS_H

#include "complx.h"
#include "matrix.h"
#include "spinsys_struct.h"
#include "sim.h"

extern ISOTOPE isotopes[];

void ss_showspins(SpinSys* S);
ISOTOPE* ss_findisotope(char* name);
double ss_qn(SpinSys* S,int spin);
double ss_gamma1H();
double ss_gamma(SpinSys* S,int spin);
ISOTOPE* ss_isotope(SpinSys* S,int spin);
int ss_matdim(SpinSys* S);
void ss_addspin(SpinSys* S,char* name);
void ss_initialize(SpinSys* S);
int ss_issame(SpinSys* S,int spin1,int spin2);

mat_complx * ss_qdiag(Sim_info* sim);
mat_complx * Iqdelta(Sim_info* sim);
mat_complx * Icoherence(Sim_info* sim, double* coh);

/**** tip: this could be re-coded to real matrices  *****/
//mat_complx * Ie(Sim_info* sim);
//mat_complx * Ic(Sim_info* sim,int spin);
mat_complx * Ip(Sim_info* sim,int spin);
mat_complx * Im(Sim_info* sim,int spin);
mat_complx * Ix(Sim_info* sim,int spin);
mat_complx * Iy(Sim_info* sim,int spin);
mat_complx * Iz(Sim_info* sim,int spin);

/* tip: if the above re-coded, this still returns complex matrices */
//mat_complx * ss_oper(Sim_info* sim,char* name);
mat_complx * ss_readoper(Sim_info *sim,char* sop);
mat_complx * ss_eq_sigma(Sim_info* sim);

/* returns result as real matrix, NOT permuted! */
mat_double * Ip_real(Sim_info* sim,int spin);
mat_double * Im_real(Sim_info* sim,int spin);
mat_double * Ix_real(Sim_info* sim,int spin);
mat_double * Iy_real(Sim_info* sim,int spin);

/*** for block diagonalized Hamiltonians ***/
mat_double * Iz_ham(Sim_info* sim, int nuc);
blk_mat_double * IzIz_sqrt2by3(Sim_info *sim, int n1, int n2);
blk_mat_double * II(Sim_info* sim, int n1, int n2);
blk_mat_double * T20(Sim_info* sim, int n1, int n2);
mat_double * T20II(Sim_info *sim, int nuc);
void fill_Tquad_2(Sim_info *sim, int nuc, Quadrupole *qptr);
void fill_Tquad_3(Sim_info *sim, Quadrupole *qptr);
void fill_Tmix_dipole(Sim_info *sim, Mixing *mptr);
void fill_Tabmix_dipole(Sim_info *sim, Mixing *mptr);


#endif


