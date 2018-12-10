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
    
    Makes available som Tcl function to 
      - evaluate statically linkes Tcl code
      - get parameters from the 'par' array.
      - signal handling so Ctrl-C is handled properly
*/


#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <signal.h>
#include <stdarg.h>
#include "defs.h"
#include "cm.h"
#include "tclutil.h"

	/* functions that define new Tcl commands from C level */
	void tclcmd_rfshape(Tcl_Interp* interp);
	void tclcmd_simpson(Tcl_Interp* interp);
	void tclcmd_ftools(Tcl_Interp* interp);
	void tclcmd_spinsys(Tcl_Interp* interp);
	void tclcmd_pulse(Tcl_Interp* interp);
	void tclcmd_OCroutines(Tcl_Interp* interp);
	void tclcmd_spinach(Tcl_Interp* interp);
	void tclcmd_distortions(Tcl_Interp* interp);

void TclSimpsonInterpInit(Tcl_Interp* interp)
{
	int i = 0;

	while ( strlen(tclcode_pointers[i].name)) {
		/* execute all except 'slave' */
		if ( strncmp(tclcode_pointers[i].name,"slave",5) ) {
			DEBUGPRINT("Initialize: evaluating code from '%s'\n",tclcode_pointers[i].name);
			if (Tcl_Eval(interp,tclcode_pointers[i].code) != TCL_OK) {
				fprintf(stderr,"Error when evaluating initialization script from '%s':\n%s\n",tclcode_pointers[i].name,Tcl_GetStringResult(interp));
				exit(1);
			}
		}
		i++;
	}

	DEBUGPRINT("Initialize: internal code all done, now adding C commands.\n");
	tclcmd_rfshape(interp);
	tclcmd_simpson(interp);
	tclcmd_ftools(interp);
	tclcmd_spinsys(interp);
	tclcmd_pulse(interp);
	tclcmd_OCroutines(interp);
	tclcmd_spinach(interp);
	tclcmd_distortions(interp);

	DEBUGPRINT("Initialize: all done.\n");
}

void TclSlaveInterpInit(Tcl_Interp* interp)
{
	int i = 0;

	while ( strlen(tclcode_pointers[i].name)) {
		/* execute only 'slave' */
		if ( !strncmp(tclcode_pointers[i].name,"slave",5) ) {
			DEBUGPRINT("Initialize slave: evaluating code from '%s'\n",tclcode_pointers[i].name);
			if (Tcl_Eval(interp,tclcode_pointers[i].code) != TCL_OK) {
				fprintf(stderr,"Error when evaluating initialization script from '%s':\n%s\n",tclcode_pointers[i].name,Tcl_GetStringResult(interp));
				exit(1);
			}
			break;
		}
		i++;
	}

	DEBUGPRINT("Initialize slave: internal code all done, now adding C commands.\n");
	tclcmd_rfshape(interp);
	tclcmd_simpson(interp);
	tclcmd_ftools(interp);
	tclcmd_spinsys(interp);
	tclcmd_pulse(interp);
	tclcmd_OCroutines(interp);
	tclcmd_spinach(interp);
	tclcmd_distortions(interp);

	DEBUGPRINT("Initialize slave: all done.\n");

}



int TclAppendResult(Tcl_Interp* interp,const char* format, ...)
{
   char buffer [512];
   va_list argptr;
   va_start(argptr, format);
   vsprintf(buffer, format, argptr);
   va_end(argptr);
   Tcl_AppendElement(interp,buffer);
   return TCL_OK;
}

int TclSetResult(Tcl_Interp* interp,const char* format, ...)
{
   char buffer [512];
   va_list argptr;
   va_start(argptr, format);
   vsprintf(buffer, format, argptr);
   va_end(argptr);
   Tcl_SetResult(interp,buffer,TCL_VOLATILE);
   return TCL_OK;
}

int TclError(Tcl_Interp* interp,const char* format, ...)
{
   char buffer [512];
   va_list argptr;
   va_start(argptr, format);
   vsprintf(buffer, format, argptr);
   va_end(argptr);
   Tcl_SetResult(interp,buffer,TCL_VOLATILE);   
   return TCL_ERROR;
}


int TclGetInt(Tcl_Interp* interp,char *aryname,char* varname,
                     int mustexist,int defval)
{
  int  val;
  Tcl_Obj* src;

  src = Tcl_GetVar2Ex(interp,aryname,varname,0);
  if ( src == NULL) {
     if (mustexist) {
       fprintf(stderr,"error: could not read integer variable %s(%s)\n",aryname,varname);
       exit(1);
     }
     if (verbose & VERBOSE_PAR)
       printf("integer variable %s in array %s is set to default value %d\n",varname,aryname,defval);

     return defval;
  }

  if (Tcl_GetIntFromObj(interp,src,&val) != TCL_OK) return TclError(interp,"GetInt(2)");
 
  if (verbose & VERBOSE_PAR)
    printf("integer variable %s in array %s is set to %d\n",varname,aryname,val);
  return val;
}

double TclGetDouble(Tcl_Interp* interp,char *aryname,char* varname,
                     int mustexist,double defval)
{
  double val;
  Tcl_Obj* src;
  
  src = Tcl_GetVar2Ex(interp,aryname,varname,0);
  if (src == NULL) {
     if (mustexist) {
       fprintf(stderr,"error: could not read double variable %s(%s)\n",aryname,varname);
       exit(-1);
     }
     if (verbose & VERBOSE_PAR)
       printf("double variable %s in array %s is set to default value %g\n",varname,aryname,defval);
     return defval;
  }
  if (Tcl_GetDoubleFromObj(interp,src,&val) != TCL_OK) return TclError(interp,"GetDouble(2)");
 
  if (verbose & VERBOSE_PAR)
    printf("double variable %s in array %s is set to %g\n",varname,aryname,val);
  return val;
}

char* TclGetString(Tcl_Interp* interp,char *dst,char* aryname,char* varname,
                     int mustexist,char* defval)
{
  Tcl_Obj* src;

  src = Tcl_GetVar2Ex(interp,aryname,varname,0);
  if (src == NULL) {
     if (mustexist) {
       fprintf(stderr,"error: could not read string variable %s(%s)\n",aryname,varname);
       exit(-1);
     }
     if (verbose & VERBOSE_PAR)
       printf("string variable %s in array %s is set to default value %s\n",varname,aryname,defval);
     strcpy(dst,defval);
     return dst;
  }
  strcpy(dst,Tcl_GetString(src));

  if (verbose & VERBOSE_PAR) {
    printf("string variable %s in array %s is set to %s\n",varname,aryname,dst);
  }
  return dst;
}


char* sigtxt[] = {
"-", /* 0 */
"-", /* 1 */
"Interrupt from keyboard", /* 2 */
"Quit from keyboard", /* 3 */
"-", /* 4 */
"-", /* 5 */
"Abort", /* 6 */
"-", /* 7 */
"-", /* 8 */
"Kill signal", /* 9 */
"-", /* 10 */
"-", /* 11 */
"-", /* 12 */
"-", /* 13 */
"-", /* 14 */
"Termination signal", /* 15 */
};

int lastsig;

int TclSignalHandler(ClientData clientData,Tcl_Interp *interp,int code)
{
   char buf[256];
   if (interp == NULL) {
     fprintf(stderr,"signalhandler was called with code %d "
                    "when no Tcl interpreter was available\n",code);
     return TCL_OK;
   }
   sprintf(buf,"signalhandler %d \"%s\"",lastsig,sigtxt[lastsig]);
   return Tcl_Eval(interp,buf);
}           

Tcl_AsyncHandler asynchandler;

void signal_handler(int sig)
{
  lastsig=sig;
  Tcl_AsyncMark(asynchandler);
}

void TclSetSignalHandler(Tcl_Interp* interp,char* function)
{
  asynchandler=Tcl_AsyncCreate(TclSignalHandler,(ClientData)NULL);

#ifdef SIGINT
  signal(SIGINT,signal_handler);
#else
  signal(SIGBREAK,signal_handler);
#endif
}

#define BUFLEN 16384

int TclAppendMatrix(Tcl_Interp* interp,mat_complx * m)
{
   int i,j,r,c;
   char num[64];
   char buf[BUFLEN];
   int nbuf;
   complx z;
   
   r = m->row;
   c = m->col;
   for (i=0;i<r;i++) {
     nbuf=0;
     buf[0]=0;
     for (j=0;j<c;j++) {
        z = cm_getelem(m,i+1,j+1);
        sprintf(num,"{%g %g}",z.re,z.im);
        nbuf += strlen(num);
        if (nbuf >= BUFLEN) {
          Tcl_SetResult(interp,"getmatrix: internal buffer overflow\n",TCL_STATIC);
          return TCL_ERROR;          
        }
        if (j != 0) strcat(buf," ");
        strcat(buf,num);
     } 
     Tcl_AppendElement(interp,buf);
  }
  // DEBUGPRINT("TclAppendMatrix created interp result:\n%s\n",Tcl_GetStringResult(interp));
  return TCL_OK;
}
