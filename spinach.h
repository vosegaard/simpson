/*
 * spinach.h
 *
 *  Created on: Aug 30, 2010
 *      Author: zdenek
 */

#ifndef SPINACH_H_
#define SPINACH_H_

#include "tcl.h"
#include "sim.h"

void spinach_disable_command(Sim_info *sim, char *cmd);
void spinach_pulseid(Sim_info *sim,Sim_wsp *wsp,double duration);
void spinach_pulse(Sim_info *sim,Sim_wsp *wsp,double duration);
void spinach_delay(Sim_info *sim,Sim_wsp *wsp,double duration);
void spinach_acqblock(Tcl_Interp *interp, Tcl_Obj *obj, Sim_info *sim, Sim_wsp *wsp);

#endif /* SPINACH_H_ */
