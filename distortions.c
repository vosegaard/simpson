/*
 * distortions.c - utilities  for RF shape distortions characterized by
 *                 Impulse Response Functions (IRF)
 *
 *  Created on: 15. 10. 2015
 *      Author: Zdenìk Tošner
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include <tcl.h>
 #include <math.h>
 #include <ctype.h>
 #include <unistd.h>
 #include <string.h>
 #include <time.h>
 #include <assert.h>
 #include "defs.h"
 #include "tclutil.h"
 #include "rfshapes.h"
 #include "spinsys_struct.h"
 #include "cm.h"
#ifdef INTEL_MKL
#include "mkl.h"
#include "mkl_spblas.h"
#elif defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#elif defined(GSL)
#include <gsl/gsl_cblas.h>
#else
#include "cblas.h"
#endif
#include "iodata.h"
#include "OCroutines.h"

 /* define global array of matrices used to store distortion operators */
 mat_complx *DOPs[MAXCHAN];

 /****
  * DOPs array initialization
  ****/
 void DOPs_init() {
    int a;

    for (a=0; a<MAXCHAN; a++) {
       DOPs[a] = NULL;
    }
 }

 /****
  * provide free slot in DOPs array
  ****/
 int DOPs_slot() {
    int a;

    for (a=0; a<MAXCHAN; a++) {
       if (!DOPs[a]) {
          break;
       }
    }
    if (a >= MAXCHAN) {
       fprintf(stderr,"DOPs error: no more free slots available\n");
       a=-1;
    }
    return a;
 }



 /****
  * implementation of Tcl create_distortion_operator routine
  ****/
int tclCreateDOP(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
	int slot, i, j, k, M, N, ratio, Nirf;
	double dt_small, dt_big, re, im;
	char *fname, name[256], dum[256];
	complx *irf, z;
	//Tcl_Obj *tclres;
	FILE* fp;

	if ( argc != 5)
		return TclError(interp,"usage: create_distortion_operator <IRF-filename> <sampling step> <# of input elements> <input sampling step>");

	fname = Tcl_GetString(argv[1]);
	if (Tcl_GetDoubleFromObj(interp,argv[2],&dt_small) == TCL_ERROR)
		return TclError(interp,"create_distortion_operator: argument 2 must be double <IRF sampling step>");
	if (Tcl_GetIntFromObj(interp,argv[3],&N) == TCL_ERROR)
	    return TclError(interp,"create_distortion_operator: argument 3 must be integer <num of input elements>");
	if (Tcl_GetDoubleFromObj(interp,argv[4],&dt_big) == TCL_ERROR)
		return TclError(interp,"create_distortion_operator: argument 4 must be double <input sampling step>");
	// debug
	//printf("\ncreate_distortion_operator %s %g %d %g\n",fname,dt_small,N,dt_big);

	// check time step parameters
	ratio = (int)floor(dt_big/dt_small);
	if (fabs(ratio*dt_small-dt_big) > 1e-10) ratio++;
	if (fabs(ratio*dt_small-dt_big) > 1e-10)
		return TclError(interp,"create_distortion_operator: input sampling (%g) has to be multiple of IRF sampling (%g)",dt_big, dt_small);
	// debug
	//printf("ratio = %d\n",ratio);

	slot = DOPs_slot();
	// load the IRF into memory
	strcpy(name,fname);
#ifdef UNIX
	if (fname[0] == '~') {
		char* p=getenv("HOME");
		if (p != NULL) {
			strcpy(fname,p);
			strcat(name,&fname[1]);
		}
	}
#endif
	fp=fopen(name,"r");
	if (!fp) {
		fprintf(stderr,"create_distortion_operator error: unable to open file %s\n\n",name);
		exit(1);
	}
	/* scan file for number of lines, i.e. number of elements in rf shape */
	Nirf = 0;
	while ( fgets(dum, 256, fp) ) {
		Nirf++;
	}
	fseek(fp, 0, SEEK_SET);
	//printf("Number of lines = %d\n",Nirf);
	irf = complx_vector(Nirf);
	z.re = 0; z.im = 0;
	for (i=1; i<=Nirf; i++) {
		fgets(dum, 256, fp);
		if ( sscanf(dum,"%lg%lg",&re,&im) != 2 ) {
			fprintf(stderr,"create_distortion_operator error: unable to read line %d in %s\n",i,name);
			exit(1);
		}
		irf[i].re = re;
		irf[i].im = im;
		z.re += re;
		z.im += im;
		//printf("row %d: %g %g\n",i,irf[i].re,irf[i].im);
	}
	fclose(fp);
	// scale the IRF properly
	z = Cdiv(Cunit,z);
	cv_mulc(irf,z);

	// allocate distortion operator
	M = ratio*N;
	DOPs[slot] = complx_matrix(M,N,MAT_DENSE,1,0);
	cm_zero(DOPs[slot]);
	for (i=0; i<M; i++) {
		for (j=0; j<N; j++) {
			z = Cnull;
			for (k=j*ratio; k<(j+1)*ratio; k++) {
				int idx = i-k+1;
				if (idx > 0 && idx <= Nirf) {
					z.re += irf[idx].re;
					z.im += irf[idx].im;
				}
			}
			DOPs[slot]->data[i+j*M] = z;
		}
	}
	// debug
	//cm_print(DOPs[slot],"DOP");

	free_complx_vector(irf);

	//printf("...\n");
	return TclSetResult(interp,"%d",slot);
}

int tclDistortShape(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
	int dslot, shape_in, shape_out, M, N, i;
	complx *input, *output;

	if ( argc != 4)
		return TclError(interp,"usage: distort_shape <distortion operator> <input shape> <output shape>");
	if (Tcl_GetIntFromObj(interp,argv[1],&dslot) == TCL_ERROR)
	    return TclError(interp,"distort_shape: argument 1 must be integer <distortion operator>");
	if (Tcl_GetIntFromObj(interp,argv[2],&shape_in) == TCL_ERROR)
	    return TclError(interp,"distort_shape: argument 2 must be integer <input shape>");
	if (Tcl_GetIntFromObj(interp,argv[3],&shape_out) == TCL_ERROR)
	    return TclError(interp,"distort_shape: argument 3 must be integer <output shape>");
	// parameter checks
	if (DOPs[dslot] == NULL)
		return TclError(interp,"distort_shape: wrong reference to distortion operator");
	if (RFshapes[shape_in] == NULL)
		return TclError(interp,"distort_shape: wrong reference to input shape");
	if (RFshapes[shape_out] == NULL)
		return TclError(interp,"distort_shape: wrong reference to output shape");
	M = DOPs[dslot]->row;
	N = DOPs[dslot]->col;
	//printf("DOP = [%d, %d], in = %d, out = %d\n",M,N,RFshapes_len(shape_in),RFshapes_len(shape_out));
	if (RFshapes_len(shape_in) != N)
		return TclError(interp,"distort_shape: distortion operator size not compatible with input shape length");
	if (RFshapes_len(shape_out) != M)
		return TclError(interp,"distort_shape: distortion operator size not compatible with output shape length");

	// convert input shape (ampl,phase) into (re,im)
	input = complx_vector(N);
	output = complx_vector(M);
	for (i=1; i<=N; i++) {
		input[i].re =RFshapes[shape_in][i].ampl*cos(RFshapes[shape_in][i].phase*DEG2RAD);
		input[i].im =RFshapes[shape_in][i].ampl*sin(RFshapes[shape_in][i].phase*DEG2RAD);
	}
	// distort and resample
	cv_matmul(output,DOPs[dslot],input);
	// convert output (re,im) into (ampl,phase)
	for (i=1; i<=M; i++) {
		RFshapes[shape_out][i].ampl = sqrt(output[i].re*output[i].re+output[i].im*output[i].im);
		RFshapes[shape_out][i].phase = RAD2DEG*atan2(output[i].im,output[i].re);
	}

	free_complx_vector(input);
	free_complx_vector(output);
	//printf("distort_shape done\n");
	return TCL_OK;
}

int tclReconstructGradient(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
	int fidNin, fidNout, Nsh, *dslots, M, N, i, j, Ngrout, Ngrin;
	complx *input, *output;
	extern FD **fd;
	extern int nfd;

	if (!OCpar.grad_shapes)
		return TclError(interp,"Error: recontruct_gradient can be used only in OC gradient procedure");
	if ( argc < 4)
		return TclError(interp,"usage: recontruct_gradient <fid grads distorted> <fid grads input> <distorsion op shape 1> ?<distorsion op shape 2>? ...");
	if (Tcl_GetIntFromObj(interp,argv[1],&fidNout) == TCL_ERROR)
	    return TclError(interp,"recontruct_gradient: argument 1 must be integer <fid grads distorted>");
	if (fidNout < 1 || fidNout > nfd || fd[fidNout] == NULL)
	    return TclError(interp,"recontruct_gradient: distorted gradient data set %d was not previously loaded\n",fidNout);
	if (Tcl_GetIntFromObj(interp,argv[2],&fidNin) == TCL_ERROR)
	    return TclError(interp,"recontruct_gradient: argument 2 must be integer <fid grads input>");
	if (fidNin < 1 || fidNin > nfd || fd[fidNin] == NULL)
	    return TclError(interp,"recontruct_gradient: input gradient data set %d was not previously loaded\n",fidNin);
	Nsh = argc - 3;
	if (Nsh != LEN(OCpar.grad_shapes))
		return TclError(interp,"reconstruct_gradient: number of distortion operators does not match number of grad shapes");
	dslots = int_vector(Nsh);
	Ngrin = 0;
	Ngrout = 0;
	for (i=1; i<=Nsh; i++) {
		if (Tcl_GetIntFromObj(interp,argv[i+2],&j) == TCL_ERROR)
		    return TclError(interp,"recontruct_gradient: argument %d must be integer <distortion op shape %d>",i+2,i);
		if (DOPs[j] == NULL)
			return TclError(interp,"recontruct_gradient: distortion operator %d indicated for shape %d does not exist>",j,i);
		if (DOPs[j]->row != RFshapes_len(OCpar.grad_shapes[i]))
			return TclError(interp,"reconstruc_gradient: distortion operator %d rows does not match OC grad shape %d",j,i);
		dslots[i] = j;
		Ngrout += DOPs[j]->row;
		Ngrin += DOPs[j]->col;
	}
	// check lengths of gradient fids
	if (fd[fidNout]->np != Ngrout)
		return TclError(interp,"reconstruc_gradient: length of distorted grads (%d) does not match sum of distorted shape lengths (%d)",fd[fidNout]->np,Ngrout);
	if (fd[fidNin]->np != Ngrin)
		return TclError(interp,"reconstruc_gradient: length of normal grads (%d) does not match sum of original shape lengths (%d)",fd[fidNin]->np,Ngrin);
	// do the transformation
	input = fd[fidNin]->data+1;
	output = fd[fidNout]->data+1;
	for (i=1; i<=Nsh; i++) {
		M = DOPs[dslots[i]]->row;
		N = DOPs[dslots[i]]->col;
		assert(DOPs[dslots[i]]->type == MAT_DENSE);
		cblas_zgemv(CblasColMajor,CblasTrans,M,N,&Cunit,DOPs[dslots[i]]->data,M,output,1,&Cnull,input,1);
		input += N;
		output += M;
	}
	free_int_vector(dslots);

	return TCL_OK;
}

 void tclcmd_distortions(Tcl_Interp* interp) {

  Tcl_CreateObjCommand(interp,"create_distortion_operator",(Tcl_ObjCmdProc *)tclCreateDOP,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"distort_shape",(Tcl_ObjCmdProc *)tclDistortShape,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"reconstruct_gradient",(Tcl_ObjCmdProc *)tclReconstructGradient,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

 }
