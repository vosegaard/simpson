#ifndef _B0zprofile_H
#define _B0zprofile_H

#include "sim.h" 

#define MAXZGRADSHAPES 32

extern double* ZgradShapes[];

void prepare_zaveraging(Tcl_Interp* interp,double** zvals,double** zoffsetvals);
void set_inhom_offsets(Sim_info* s, Sim_wsp* wsp, double zoffnominal);
void get_chan_offset_ratios(SpinSys* ss, double zoffnominal, double* ovals);

int ZgradShapes_slot();
double* ZgradShapes_alloc(int len);
void free_ZgradShapes(int a);
int ZgradShapes_len(int slot);

/*
double* get_zlims(int Nz) 
*/

#endif
