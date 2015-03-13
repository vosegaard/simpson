/*
    Pulse propagation
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

#ifndef __PULSE_H
#define __PULSE_H

#include "sim.h"

void direct_propagate(Sim_info *sim, Sim_wsp *wsp);

/********* ZT: few extra definitions for OC to work ************/
void _evolve_with_prop(Sim_info *sim, Sim_wsp *wsp);
void _reset_prop(Sim_info *sim, Sim_wsp *wsp);
void _ph(Sim_wsp *wsp, int channel,double phase);
void _rf(Sim_wsp *wsp, int channel,double rffield);
void _pulse(Sim_info *sim, Sim_wsp *wsp, double duration);
void _delay(Sim_info *sim, Sim_wsp *wsp, double duration);
void _pulseid(Sim_info *sim,Sim_wsp *wsp,double duration);
int _setrfprop(Sim_info *sim, Sim_wsp *wsp);

#endif /* __Pulse_H */
