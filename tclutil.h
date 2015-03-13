/*
    Tcl/C utillity routines
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


#ifndef __TCLUTIL_H
#define __TCLUTIL_H
  
#include <tcl.h>
#include "complx.h"
#include "matrix.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct TCLCODE_ { char name[32]; char* code; } TCLCODE;

extern TCLCODE tclcode_pointers[];

int TclAppendResult(Tcl_Interp* interp,const char* format, ...);
int TclError(Tcl_Interp* interp,const char* format, ...);
int TclSetResult(Tcl_Interp* interp,const char* format, ...);

//int TclAppendRealMatrix(Tcl_Interp* interp,mv_double * m);
int TclAppendMatrix(Tcl_Interp* interp,mat_complx * m);
int TclGetInt(Tcl_Interp* interp,char *aryname,char* varname,
                     int mustexist,int defval);
double TclGetDouble(Tcl_Interp* interp,char *aryname,char* varname,
                     int mustexist,double defval);

char* TclGetString(Tcl_Interp* interp,char *dst,char* aryname,char* varname,
                     int mustexist,char* defval);
                     
double* TclGetVector(Tcl_Interp* interp,char* aryname,char* varname,
                     int mustexist,double* defval);

void TclSetSignalHandler(Tcl_Interp* interp,char* function);

void TclSimpsonInterpInit(Tcl_Interp* interp);
void TclSlaveInterpInit(Tcl_Interp* interp);

#ifdef __cplusplus
}
#endif


#endif /* __TCLUTIL_H */
