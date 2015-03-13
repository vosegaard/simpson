/*
    Hamiltonian calculation
    Copyright (C) 1999 Mads Bak, Jimmy T. Rasmussen
                  2002 Thomas Vosegaard
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

#ifndef __HAM_H
#define __HAM_H

#include "matrix.h"
#include "sim.h"

void ham_print(Sim_info *s, Sim_wsp *wsp);
void ham_set_offset(Sim_info *s, Sim_wsp *wsp, double* offset,int used);

void ham_turnoff(Sim_info *s, Sim_wsp *wsp,char* name);
void ham_turnon(Sim_info *s, Sim_wsp *wsp,char* name);

int ham_ischanged(Sim_info *s, Sim_wsp *wsp);
int shift_exist(Sim_info* s,int n);
int dipole_exist(Sim_info* s,int n1, int n2);
int jcoupling_exist(Sim_info* s,int n1, int n2);
int quadrupole_exist(Sim_info* s,int n);

void ham_rotate(Sim_info *s, Sim_wsp *wsp);

void ham_hamilton(Sim_info *s, Sim_wsp *wsp);
void ham_hamilton_integrate(Sim_info *s, Sim_wsp *wsp, double dur);
blk_mat_complx * ham_rf(Sim_wsp *wsp);

/* DOR is not really finished implementation */

 
#endif
