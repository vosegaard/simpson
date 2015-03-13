/*
    Tcl/C code performing the simulation
    Copyright (C) 1999 Mads Bak. Modifications by Thomas Vosegaard, 2001.
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
    
    Makes available the 'internalsimpson' command used by the
    'fsimpson' Tcl function located in 'simpson.tcl'
    Uses the commands given in sim.c that performs the simulation.
*/

#include <stdio.h>
#include <string.h>
#include <tcl.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
#include <limits.h>
#include <pthread.h>
#ifdef WIN32
//#include <windows.h>
//#include <winbase.h>
#endif

#ifndef NO_NFFT
//#include <complex.h>
#include "nfft3util.h"
#include "nfft3.h"
#endif

#include "defs.h"
#include "tclutil.h"
#include "iodata.h"
#include "cryst.h"
#include "rfprof.h"
#include "matrix.h"
#include "cm.h"
#include "sim.h"
#include "rfshapes.h"
#include "OCroutines.h"
#include "B0inhom.h"
#include "averaging.h"
#include "rfshapes.h"
#include "OCroutines.h"
#include "fidcalc.h"
#include "pthread_barrier_mac.h"

//extern void simpson_nfft_test(void);

extern pthread_barrier_t simpson_b_start, simpson_b_end;
extern glob_info_type glob_info;

typedef struct _td{
  int Ntot, ncr, nrf, nz, Nacq;
  Sim_info *sim;
  const char *state;
  double *zvals;
  double *zoffsetvals;
  Ave_elem *ave_struct;
  int Navepar, Naveval;
  double *ave_weight;
  complx **fids;
  double phase;
} thread_data;

thread_data thrd;

int getbits(int* e,char *bits)
{
  int i;
  char* p;

  *e=0;
  for (p=bits,i=0;*p != 0 && i < sizeof(int)*8;p++,i++) {
    if (*p == '1')
      *e |= 1 << i; 
    else if (*p != '0')
      return 0;
  }
  return 1;
}

void putbits(char *bits,int e)
{
  int i;
  char* p;

  for (p=bits,i=0; i < sizeof(int)*8;p++,i++) {
    if (e & (1 << i))
      *p = '1';
    else
      *p = '0';
  }
  *p=0;
}

void thread_work(int thread_id, Tcl_Interp *interp){
  Sim_wsp * wsp;
  complx *fidsum;
  double weight;
  int i, Np_start, ncr_start, icr, nrf_start, irf, nz_start, iz;
  int ave_start;

  Sim_info *sim = thrd.sim;
  int Naveval = thrd.Naveval;

  fidsum = thrd.fids[thread_id];
  cv_zero(fidsum);
  wsp = wsp_initialize(sim);
  /* store sim and wsp memory addresses to Tcl interp */
  store_sim_pointers(interp, sim, wsp);
  wsp->interp = interp;
  wsp->thread_id = thread_id;

  i=0;
  while (1) {
	  // WARNING!!! This works ONLY IF all MPI slaves have the same number of threads!!!
	  Np_start = glob_info.mpi_rank + thread_id*glob_info.mpi_size + glob_info.mpi_size*glob_info.num_threads*i;
	  if (Np_start >= thrd.Ntot) break;
	  ncr_start = Np_start/thrd.nrf/thrd.nz/Naveval;
	  nrf_start = (Np_start-ncr_start*thrd.nrf*thrd.nz*Naveval)/thrd.nz/Naveval;
	  nz_start = (Np_start-ncr_start*thrd.nrf*thrd.nz*Naveval-nrf_start*thrd.nz*Naveval)/Naveval;
	  ave_start = Np_start - nz_start*Naveval - nrf_start*Naveval*thrd.nz - ncr_start*Naveval*thrd.nrf*thrd.nz;
	  icr = ncr_start+1;
	  irf = nrf_start+1;
	  iz = nz_start+1;
	  weight = sim->crdata[icr].weight * sim->rfdata[irf][sim->ss->nchan+1] * thrd.ave_weight[ave_start];
	  wsp->cryst =sim->crdata[icr];
	  wsp->cryst_idx = icr;
	  wsp->rffac = sim->rfdata[irf];
	  wsp->zcoor = thrd.zvals[iz];
	  set_inhom_offsets(sim,wsp,thrd.zoffsetvals[iz]);
	  set_averaging_parameters(sim,wsp,thrd.ave_struct,thrd.Navepar,ave_start);
	  sim_calcfid(sim, wsp);
	  cv_multod(fidsum, wsp->fid, weight);
	  //printf("process (%i/%i): fidsum[end].re = %f\n", process_id, thread_id, fidsum[sim->ntot].re);
	  if (verbose & VERBOSE_PROGRESS) {
		  printf("Worker %i/%i: [cryst %d/%d]",thread_id+1,glob_info.mpi_rank+1,icr,thrd.ncr);
		  if (thrd.nrf > 1) printf(" [rfsc %d/%d]",irf,thrd.nrf);
		  if (thrd.nz > 1) printf(" [z-coor %d/%d]",iz,thrd.nz);
		  if (thrd.Navepar > 0) printf(" [ave %d/%d]",ave_start+1,Naveval);
		  printf("\n");
		  //printf("Worker %i/%i: [cryst %d/%d] [rfsc %d/%d] [z-coor %d/%d]\n", thread_id+1,process_id+1,icr,ncr,irf,nrf,iz,nz);
		  fflush(stdout);
	  }
	  i++;
  }

  wsp_destroy(sim, wsp);
  free(wsp);
}

void thread_work_interpol_source(int thread_id, Tcl_Interp *interp) {
	Sim_wsp * wsp;
	int i, icr;

	Sim_info *sim = thrd.sim;

	wsp = wsp_initialize(sim);
	/* store sim and wsp memory addresses to Tcl interp */
	store_sim_pointers(interp, sim, wsp);
	wsp->interp = interp;
	wsp->thread_id = thread_id;

	set_averaging_parameters(sim,wsp,thrd.ave_struct,thrd.Navepar,thrd.Naveval);
	wsp->rffac = sim->rfdata[thrd.nrf];
	wsp->zcoor = thrd.zvals[thrd.nz];
	set_inhom_offsets(sim,wsp,thrd.zoffsetvals[thrd.nz]);

	i = -1;
	while (1) {
		i++;
		// WARNING!!! This works ONLY IF all MPI slaves have the same number of threads!!!
		icr = 1 + thread_id + glob_info.num_threads*glob_info.mpi_rank + glob_info.mpi_size*glob_info.num_threads*i;
		if (icr > thrd.ncr) break;
		//if (icr == sim->icr_done) continue;
		wsp->cryst = sim->crdata[icr];
		wsp->cryst_idx = icr;
		sim_calcfid_interpol(sim, wsp);
		if (verbose & VERBOSE_PROGRESS) {
			printf("Worker %i/%i: [cryst %d/%d] for interpolation\n",thread_id+1,glob_info.mpi_rank+1,icr,thrd.ncr);
			fflush(stdout);
		}
	}

	wsp_destroy(sim, wsp);
	free(wsp);
}

/******** NFFT library does not work on win with threads...
void thread_work_interpol_enlarge(int thread_id) {
}
***************/

/****
 * FWT interpolation on the level of master thread
 ****
void mpi_work_interpolate(Sim_info *sim)
{
	int ncr = LEN(sim->crdata);
	int ntcr = LEN(sim->targetcrdata);
	int idata, i, j, k, icr, nnztmp, Npoints, Ndata;
	char buf[256];
	FILE *fh, *fout;
	long fpos, flen, header;
	int *intdata = NULL;
	int info[2];
	double d, period;
	complx z;
	nfsft_plan sourceplan, targetplan;

	// open output file
	sprintf(buf,"%s_out_P%d.bin",sim->parname,glob_info.mpi_rank);
	fout = fopen(buf,"wb");
	if (!fout) {
		fprintf(stderr,"Error: worker (%d/master) can not create temporary file %s\n",glob_info.mpi_rank+1,buf);
		exit(1);
	}
	//printf("Worker %d/master opened file %s\n",glob_info.mpi_rank+1, buf);

	sprintf(buf,"%s_P%d_T%d.bin",sim->parname,0,0); // the first source file
	fh = fopen(buf,"rb");
	if (!fh) {
		fprintf(stderr,"Error: worker (%d/master) can not access temporary file %s\n",glob_info.mpi_rank+1,buf);
		exit(1);
	}
	//printf("Worker %d/master opened file %s\n",glob_info.mpi_rank+1, buf);

	fread(&nnztmp,sizeof(int),1,fh); // analyze the header
	fwrite(&nnztmp,sizeof(int),1,fout); // and copy header info
	switch (sim->imethod) {
	case M_GCOMPUTE_TIME:
	case M_GCOMPUTE_FREQ:
		Ndata = sim->matdim + nnztmp*sim->ngamma*2;
		Npoints = 1;
		header = (sim->matdim+nnztmp+2)*sizeof(int);
		break;
	case M_DIRECT_TIME:
	case M_DIRECT_FREQ:
		fread(&Npoints,sizeof(int),1,fh);
		fwrite(&Npoints,sizeof(int),1,fout);
		fread(&period,sizeof(int),1,fh);
		fwrite(&period,sizeof(int),1,fout);
		Ndata = (sim->matdim + nnztmp*Npoints)*sim->ngamma;
		header = (sim->matdim+nnztmp+3)*sizeof(int)+sizeof(double);
	}
	k = sim->matdim+1+nnztmp; // length of irow and icol
	intdata = (int*)malloc(k*sizeof(int));
	fread(intdata,sizeof(int),k,fh);
	//for (i=0; i<k; i++) printf(" %d ",intdata[i]);
	fwrite(intdata,sizeof(int),k,fout); // this is sparse info
	//printf(" ... stored\n");
	free(intdata);
	fclose(fh);

	nfsft_precompute(sim->Jinterpol[1],1000.0,0U,0U);
	//printf("Jinterpol = %d, %d\n",sim->Jinterpol[0],sim->Jinterpol[1]);

	// initialize nfft plans (nfsft_precompute was called from master thread)
	nfsft_init_guru(&sourceplan, sim->Jinterpol[1], ncr, NFSFT_NO_DIRECT_ALGORITHM | NFSFT_MALLOC_X | NFSFT_MALLOC_F |
			NFSFT_MALLOC_F_HAT | NFSFT_DESTROY_X | NFSFT_DESTROY_F | NFSFT_DESTROY_F_HAT,
			PRE_PHI_HUT | PRE_PSI | FFTW_INIT | FFT_OUT_OF_PLACE, 6); // NFSFT_NORMALIZED shall not be included
	//printf("sourceplan init done\n");
	for (i=1; i<=ncr; i++) { // load source orientations
		d = sim->crdata[i].alpha;
		if (d < 180.0 ) d /= 360.0; else d = d/360.0 - 1.0;
		sourceplan.x[2*(i-1)] = d;
		sourceplan.x[2*(i-1)+1] = sim->crdata[i].beta/360.0;
	}
	//printf("sourceplan.x filled \n");
	nfsft_precompute_x(&sourceplan); // Do pre-computation for nodes
	//printf("sourceplan.x precomputed \n");
	// used for initializing the odd rank terms
	//memset(sourceplan.f_hat, 0, sourceplan.N_total * sizeof(complx));
	//printf("sourceplan.f_hat nulled \n");
	// and for target
	nfsft_init_guru(&targetplan, sim->Jinterpol[1], ntcr, NFSFT_NO_DIRECT_ALGORITHM | NFSFT_MALLOC_X | NFSFT_MALLOC_F |
			NFSFT_MALLOC_F_HAT | NFSFT_DESTROY_X | NFSFT_DESTROY_F | NFSFT_DESTROY_F_HAT,
			PRE_PHI_HUT | PRE_PSI | FFTW_INIT | FFT_OUT_OF_PLACE, 6); // NFSFT_NORMALIZED shall not be included
	for (i=1; i<=ntcr; i++) { // load source orientations
		d = sim->targetcrdata[i].alpha;
		if (d < 180.0 ) d /= 360.0; else d = d/360.0 - 1.0;
		targetplan.x[2*(i-1)] = d;
		targetplan.x[2*(i-1)+1] = sim->targetcrdata[i].beta/360.0;
	}
	nfsft_precompute_x(&targetplan); // Do pre-computation for nodes
	// used for initializing the odd rank terms
	//memset(targetplan.f_hat, 0, targetplan.N_total * sizeof(complx));
	//printf("The same done for targetplan, N_total = %d, M_total = %d\n",targetplan.N_total, targetplan.M_total);


	i = 0;
	while (1) {
		idata = glob_info.mpi_rank + glob_info.mpi_size*i;
		if (idata >= Ndata) break;
		for (j=0; j<glob_info.mpi_size; j++) {
			for (k=0; k<glob_info.num_threads; k++) {
				sprintf(buf,"%s_P%d_T%d.bin",sim->parname,j,k);
				fh = fopen(buf,"rb");
				if (!fh) {
					fprintf(stderr,"Error: worker (%d/master) can not access temporary file %s\n",glob_info.mpi_rank+1,buf);
					exit(1);
				}
				fseek(fh,0,SEEK_END);
				flen = ftell(fh);
				fpos = header;
				//printf("file %s, header = %ld, idata = %d\n",buf,header,idata);
				while (fpos < flen) {
					fseek(fh,fpos,SEEK_SET);
					fread(&icr,sizeof(int),1,fh);
					//printf("icr = %d, fpos = %ld\n",icr,fpos);
					fseek(fh,idata*sizeof(complx),SEEK_CUR);
					//fread(sourcedata+icr,sizeof(complx),1,fh);
					//sourcedata[icr].re *= sim->crdata[icr].weight; // normalization for interpolation
					//sourcedata[icr].im *= sim->crdata[icr].weight;
					fread(&z,sizeof(complx),1,fh);
					sourceplan.f[icr-1][0] = sim->crdata[icr].weight * z.re;
					sourceplan.f[icr-1][1] = sim->crdata[icr].weight * z.im;
					fpos += sizeof(int)+Ndata*sizeof(complx);
				}
				fclose(fh);
			}
		}
		// set element info
		switch (sim->imethod) {
		case M_GCOMPUTE_TIME:
		case M_GCOMPUTE_FREQ:
			if (idata < sim->matdim) {
				info[0] = 0;         // matrix type
				info[1] = idata;     // element number
			} else {
					info[0] = ((idata-sim->matdim)%2)*sim->ngamma+1 + ((idata-sim->matdim)/2)%sim->ngamma; // which gamma index (1..Ng)
					info[1] = (idata - sim->matdim) / (2*sim->ngamma) ;     // element number (0..nnz-1) in amplitude matrix
			}
			break;
		case M_DIRECT_TIME:
		case M_DIRECT_FREQ:
			k = idata / (sim->matdim+nnztmp*Npoints);
			j = idata % (sim->matdim+nnztmp*Npoints);
			if (j < sim->matdim) {
				info[0] = k;
				info[1] = j;
			} else {
				j -= sim->matdim;
				info[0] = sim->ngamma + k*Npoints + (j % Npoints);
				info[1] = j / Npoints;
			}
			break;
		}
		fwrite(info,sizeof(int),2,fout);
		// sourcedata should be ready
		//sprintf(buf,"Sourcedata index %d",idata);
		//nfft_vpr_complex(sourceplan.f,8,buf); // sourceplan.M_total
		//int xyz;
		//for (xyz=0; xyz<8; xyz++)
		//	printf("xyz = %d, %g ... %g\n",xyz,sourceplan.f[xyz][0]/sim->crdata[xyz+1].weight,sourceplan.f[xyz][1]/sim->crdata[xyz+1].weight);
		// do the adjoint transform on source
		//printf("now nfsft_adjoint\n");
		nfsft_adjoint(&sourceplan);
		//printf("sourceplan transformed \n");
		// Copy only even rank terms, the rest shall remain zero
		memset(targetplan.f_hat, 0, targetplan.N_total * sizeof(complx));
		for( j = 0; j <= sim->Jinterpol[0]; j += 2) {
			double normalization = 2.0 * ((double) j) + 1.0;
			for( k = -j; k <= j; k++) {
				targetplan.f_hat[NFSFT_INDEX(j, k, &targetplan)][0] = sourceplan.f_hat[NFSFT_INDEX(j, k, &sourceplan)][0] * normalization;
				targetplan.f_hat[NFSFT_INDEX(j, k, &targetplan)][1] = sourceplan.f_hat[NFSFT_INDEX(j, k, &sourceplan)][1] * normalization;
				//targetplan.f_hat[NFSFT_INDEX(j, k, &targetplan)] = sourceplan.f_hat[NFSFT_INDEX(j, k, &sourceplan)] * normalization;
			}
		}
		//printf("targetplan.f_hat filled \n");
		// do the transformation on target
		nfsft_trafo(&targetplan);
		//printf("targetplan transformed \n");
		// result is in targetplan.f
		fwrite(targetplan.f,sizeof(complx),ntcr,fout);
		//sprintf(buf,"Index %d, interpolated target vector f(1:8)",idata);
		//nfft_vpr_complex(targetplan.f,8,buf);
		i++;
	}
	fclose(fout);
	//free_complx_vector(sourcedata);
	//free_complx_vector(targetdata);
	nfsft_finalize(&sourceplan);
	nfsft_finalize(&targetplan);
}
****/

/******
 * FWT interpolation of data stored in memory (option to interpolate only lams)
 */
void master_FWTinterpolate(Sim_info *sim)
{
#ifdef NO_NFFT
	fprintf(stderr,"Error: simpson compiled without NFFT library so no FWT interpolation possible!\n");
	exit(1);
#else
	int ncr = LEN(sim->crdata);
	int ntcr = LEN(sim->targetcrdata);
	int i, j, k, icr;
	double d;
	complx *zz;
	nfsft_plan sourceplan, targetplan;

	nfsft_precompute(sim->Jinterpol[1],1000.0,0U,0U);
	//printf("Jinterpol = %d, %d\n",sim->Jinterpol[0],sim->Jinterpol[1]);

	// initialize nfft plans
	nfsft_init_guru(&sourceplan, sim->Jinterpol[1], ncr, NFSFT_NO_DIRECT_ALGORITHM | NFSFT_MALLOC_X | NFSFT_MALLOC_F |
			NFSFT_MALLOC_F_HAT | NFSFT_DESTROY_X | NFSFT_DESTROY_F | NFSFT_DESTROY_F_HAT,
			PRE_PHI_HUT | PRE_PSI | FFTW_INIT | FFT_OUT_OF_PLACE, 6); // NFSFT_NORMALIZED shall not be included
	for (i=1; i<=ncr; i++) { // load source orientations
		d = sim->crdata[i].alpha;
		if (d < 180.0 ) d /= 360.0; else d = d/360.0 - 1.0;
		sourceplan.x[2*(i-1)] = d;
		sourceplan.x[2*(i-1)+1] = sim->crdata[i].beta/360.0;
	}
	//printf("sourceplan.x filled \n");
	nfsft_precompute_x(&sourceplan); // Do pre-computation for nodes
	//printf("sourceplan.x precomputed \n");

	nfsft_init_guru(&targetplan, sim->Jinterpol[1], ntcr, NFSFT_NO_DIRECT_ALGORITHM | NFSFT_MALLOC_X | NFSFT_MALLOC_F |
			NFSFT_MALLOC_F_HAT | NFSFT_DESTROY_X | NFSFT_DESTROY_F | NFSFT_DESTROY_F_HAT,
			PRE_PHI_HUT | PRE_PSI | FFTW_INIT | FFT_OUT_OF_PLACE, 6); // NFSFT_NORMALIZED shall not be included
	for (i=1; i<=ntcr; i++) { // load source orientations
		d = sim->targetcrdata[i].alpha;
		if (d < 180.0 ) d /= 360.0; else d = d/360.0 - 1.0;
		targetplan.x[2*(i-1)] = d;
		targetplan.x[2*(i-1)+1] = sim->targetcrdata[i].beta/360.0;
	}
	nfsft_precompute_x(&targetplan); // Do pre-computation for nodes
	//printf("The same done for targetplan, N_total = %d, M_total = %d\n",targetplan.N_total, targetplan.M_total);

	// interpolate lams
	complx *new_lam = (complx*)malloc(sim->matdim*ntcr*sizeof(complx));
	if (new_lam == NULL) {
		fprintf(stderr,"Error: no more memory for new lam vector (FWT)\n");
		exit(1);
	}
	for (i=0; i<sim->matdim; i++) {
		zz = &sim->FWT_lam[i*ncr];
		for (icr=1; icr<=ncr; icr++) {
			sourceplan.f[icr-1][0] = sim->crdata[icr].weight * zz->re;
			sourceplan.f[icr-1][1] = sim->crdata[icr].weight * zz->im;
			//printf("\t(%d)\t %d: %g %g\n",i,icr,zz->re,zz->im);
			zz++;
		}
		//printf("\n\t(%d) ...",i);
		// sourcedata should be ready, do the adjoint transform on source
		nfsft_adjoint(&sourceplan);
		// Copy only even rank terms, the rest shall remain zero
		//printf(" done ...");
		memset(targetplan.f_hat, 0, targetplan.N_total * sizeof(complx));
		//printf(" done ...");
		for( j = 0; j <= sim->Jinterpol[0]; j += 2) {
			double normalization = 2.0 * ((double) j) + 1.0;
			for( k = -j; k <= j; k++) {
				targetplan.f_hat[NFSFT_INDEX(j, k, &targetplan)][0] = sourceplan.f_hat[NFSFT_INDEX(j, k, &sourceplan)][0] * normalization;
				targetplan.f_hat[NFSFT_INDEX(j, k, &targetplan)][1] = sourceplan.f_hat[NFSFT_INDEX(j, k, &sourceplan)][1] * normalization;
			}
		}
		//printf(" done ...");
		// do the transformation on target
		nfsft_trafo(&targetplan);
		//printf(" done ...");
		// result is in targetplan.f
		memcpy(&new_lam[i*ntcr],targetplan.f,ntcr*sizeof(complx));
		//printf(" done !\n");
		//for (icr=0; icr<ntcr; icr++) {
		//	printf("(%d)\t %d: %g %g\n",i,icr,new_lam[i*ntcr+icr].re,new_lam[i*ntcr+icr].im);
		//}
		//if (i == 2) {
		//	printf("STOP\n");
		//	exit(1);
		//}
	}
	free(sim->FWT_lam);
	sim->FWT_lam = new_lam;
	//printf("lams\n");
	//for (i=0; i<ntcr*sim->matdim; i++) printf("(%d) %g, %g\n",i,sim->FWT_lam[i].re,sim->FWT_lam[i].im);

	if ( sim->interpolation == INTERPOL_FWT_ALL || sim->interpolation == INTERPOL_FWTASG_ALL) { // interpolate also amplitudes data
		if (sim->FWTASG_nnz[1] > 0) {
			fprintf(stderr,"Error: Can not interpolate frs for non-diagonal cases, use FWT2interpolation\n");
			exit(1);
		}
		int Ndata = ((sim->EDsymmetry == 1) ? 1 : 2)*(-sim->FWTASG_nnz[1])*sim->points_per_cycle;
		assert(Ndata > 0);
			//printf("%d, %d, %d -> %d\n",sim->EDsymmetry,sim->FWTASG_nnz,sim->points_per_cycle,Ndata);
			//complx *new_frs = (complx *)malloc(Ndata*ntcr*sizeof(complx));
			//if (new_frs == NULL) {
			//	fprintf(stderr,"Error: no more memory for new_frs (FWT)\n");
			//	exit(1);
			//}
		complx **frs_pool = (complx**)malloc(ntcr*sizeof(complx*));
		if (frs_pool == NULL) {
			fprintf(stderr,"Error: master_interpolateFWT - no memory for frs_pool\n");
			exit(1);
		}
		for (i=0; i<ntcr; i++) {
			frs_pool[i] = (complx*)malloc(Ndata*sizeof(complx));
			if (frs_pool[i] == NULL) {
				fprintf(stderr,"Error: master_interpolateFWT - no memory for frs_pool element %d\n",i);
				exit(1);
			}
		}
		for (i=0; i<Ndata; i++) {
			//zz = &sim->FWT_frs[i*ncr];
			//printf("\t(%d)\t",i);
			for (icr=1; icr<=ncr; icr++) {
				zz = &(sim->FWT_frs[icr-1][i]);
				sourceplan.f[icr-1][0] = sim->crdata[icr].weight * zz->re;
				sourceplan.f[icr-1][1] = sim->crdata[icr].weight * zz->im;
				//zz++;
			}
			//printf("...done");
			// sourcedata should be ready, do the adjoint transform on source
			nfsft_adjoint(&sourceplan);
			//printf("...done");
			// Copy only even rank terms, the rest shall remain zero
			memset(targetplan.f_hat, 0, targetplan.N_total * sizeof(complx));
			//printf("...done");
			for( j = 0; j <= sim->Jinterpol[0]; j += 2) {
				double normalization = 2.0 * ((double) j) + 1.0;
				for( k = -j; k <= j; k++) {
					targetplan.f_hat[NFSFT_INDEX(j, k, &targetplan)][0] = sourceplan.f_hat[NFSFT_INDEX(j, k, &sourceplan)][0] * normalization;
					targetplan.f_hat[NFSFT_INDEX(j, k, &targetplan)][1] = sourceplan.f_hat[NFSFT_INDEX(j, k, &sourceplan)][1] * normalization;
				}
			}
			//printf("...done");
			// do the transformation on target
			nfsft_trafo(&targetplan);
			//printf("...done");
			// result is in targetplan.f
			//memcpy(&new_frs[i*ntcr],targetplan.f,ntcr*sizeof(complx));
			for (icr=0; icr<ntcr; icr++) {
				frs_pool[icr][i].re = targetplan.f[icr][0];
				frs_pool[icr][i].im = targetplan.f[icr][1];
			}
			//printf("...done!\n");
		}
		//free(sim->FWT_frs);
		//sim->FWT_frs = new_frs;
		for (i=0; i<ncr; i++) free(sim->FWT_frs[i]);
		free(sim->FWT_frs);
		sim->FWT_frs = frs_pool;
	}

	nfsft_finalize(&sourceplan);
	nfsft_finalize(&targetplan);
#endif
}

/****
 * this should somehow collect all ASG data from all workers.
 * NOT implemented yet but not really needed without MPI
 */
void mpi_work_ASGread(Sim_info *sim)
{
	//int i, j, k, l, icr, npts = 0, ncr, nnz = 0;
	//char buf[256];
	//FILE *fh;
	//long fpos, flen;

	assert(sim->interpolation == INTERPOL_ASG || sim->interpolation == INTERPOL_FWTASG_ALL || sim->interpolation == INTERPOL_FWTASG_LAM);

	if (sim->interpolation == INTERPOL_ASG) {
		sim->tridata = read_triangle_file(sim->crystfile);
	} else {
		sim->tridata = read_triangle_file(sim->targetcrystfile);
	}
/*	ncr = LEN(sim->crdata);

	if (sim->imethod == M_DIRECT_FREQ) {
		for (i=0; i<glob_info.mpi_size; i++) {
			for (j=0; j<glob_info.num_threads; j++) {
				sprintf(buf,"%s_P%d_T%d.bin",sim->parname,i,j);
				fh = fopen(buf,"rb");
				if (nnz == 0) {
					fread(&nnz,sizeof(int),1,fh);
					fread(&npts,sizeof(int),1,fh);
					fread(&(sim->ASG_period),sizeof(double),1,fh);
					sim->ASG_ampl = (complx*)malloc(nnz*npts*ncr*sim->ngamma*sizeof(complx));
					sim->ASG_freq = (double*)malloc(nnz*ncr*sim->ngamma*sizeof(double));
				}
				fseek(fh,0,SEEK_END);
				flen = ftell(fh);
				fpos = 2*sizeof(int)+sizeof(double);
				fseek(fh,fpos,SEEK_SET);
				//printf("ASG reads file %s, nnz = %d\n",buf, nnz);
				while (fpos < flen) {
					fread(&icr,sizeof(int),1,fh);
					//printf("icr = %d, fpos = %ld\n",icr,fpos);
					for (k=0; k<nnz*sim->ngamma; k++) {
						fread(sim->ASG_freq+k+nnz*sim->ngamma*(icr-1),sizeof(double),1,fh);
						fread(sim->ASG_ampl+k*npts+(icr-1)*npts*nnz*sim->ngamma,sizeof(complx),npts,fh);
					}
					fpos = ftell(fh);
				}
				fclose(fh);
			}
		}
		sim->ASG_nnz = nnz;
		// END of DIRECT freq option
	} else if (sim->imethod == M_GCOMPUTE_FREQ) {
		npts = sim->ngamma;
		for (i=0; i<glob_info.mpi_size; i++) {
			for (j=0; j<glob_info.num_threads; j++) {
				sprintf(buf,"%s_P%d_T%d.bin",sim->parname,i,j);
				fh = fopen(buf,"rb");
				if (nnz == 0) {
					fread(&nnz,sizeof(int),1,fh);
					sim->ASG_ampl = (complx*)malloc(nnz*npts*ncr*sizeof(complx));
					sim->ASG_freq = (double*)malloc(nnz*ncr*sizeof(double));
				}
				fseek(fh,0,SEEK_END);
				flen = ftell(fh);
				fpos = sizeof(int);
				fseek(fh,fpos,SEEK_SET);
				//printf("ASG reads file %s, nnz = %d\n",buf, nnz);
				while (fpos < flen) {
					fread(&icr,sizeof(int),1,fh);
					//printf("icr = %d, fpos = %ld\n",icr,fpos);
					for (k=0; k<nnz; k++) {
						fread(sim->ASG_freq+k+nnz*(icr-1),sizeof(double),1,fh);
						fread(sim->ASG_ampl+k*npts+(icr-1)*npts*nnz,sizeof(complx),npts,fh);
					}
					fpos = ftell(fh);
				}
				fclose(fh);
			}
		}
		sim->ASG_nnz = nnz;
		// END of gCOMPUTE freq option
	} else {
		fprintf(stderr,"Error: ASGread - wrong imethod\n");
		exit(1);
	}
*/
	//printf("ASG data read into memory\n");
}


void thread_work_interpol_calcfid(int thread_id)
{
	int i, icr;
	Sim_info *sim = thrd.sim;
	int ncr = LEN(sim->targetcrdata);
	complx *fid = thrd.fids[thread_id];

	assert(sim->interpolation == INTERPOL_FWT_ALL || sim->interpolation == INTERPOL_FWT_LAM);

	cv_zero(fid);
	i=0;
	while (1) {
		// WARNING!!! This works ONLY IF all MPI slaves have the same number of threads!!!
		icr = 1+ thread_id + glob_info.num_threads*glob_info.mpi_rank + glob_info.mpi_size*glob_info.num_threads*i;
		if (icr > ncr) break;
		switch (sim->imethod) {
		case M_GCOMPUTE_TIME:
			if (sim->interpolation == INTERPOL_FWT_ALL) {
				collect_fid_interpol_all(icr, sim, fid);
			} else {
				collect_fid_interpol_lam(icr, sim, fid);
			}
			break;
		case M_GCOMPUTE_FREQ:
			if (sim->interpolation == INTERPOL_FWT_ALL) {
				collect_spc_interpol_all(icr, sim, fid,thread_id);
			} else {
				collect_spc_interpol_lam(icr, sim, fid,thread_id);
			}
			break;
		case M_DIRECT_TIME:
			//collect_fid_interpol_direct(icr, sim, fid);
			break;
		case M_DIRECT_FREQ:
			if (sim->interpolation == INTERPOL_FWT_ALL) {
				collect_spc_direct_interpol_all(icr, sim, fid,thread_id);
			} else {
				collect_spc_direct_interpol_lam(icr, sim, fid,thread_id);
			}
			break;
		}
		i++;
		if (verbose & VERBOSE_PROGRESS) {
			printf("Worker %i/%i: [cryst %d/%d] \n",thread_id+1,glob_info.mpi_rank+1,icr,ncr);
			fflush(stdout);
		}
	}
}

void thread_work_FWTtoASG(int thread_id)
{
	int i, icr;
	Sim_info *sim = thrd.sim;
	int ncr = LEN(sim->targetcrdata);

	assert(sim->interpolation == INTERPOL_FWTASG_ALL || sim->interpolation == INTERPOL_FWTASG_LAM);

	i=0;
	while (1) {
		// WARNING!!! This works ONLY IF all MPI slaves have the same number of threads!!!
		icr = 1+ thread_id + glob_info.num_threads*glob_info.mpi_rank + glob_info.mpi_size*glob_info.num_threads*i;
		if (icr > ncr) break;
		switch (sim->imethod) {
		case M_GCOMPUTE_TIME:
			fprintf(stderr,"Error: thread_work_FWTtoASG - wrong imethod %d\n",sim->imethod);
			exit(1);
			break;
		case M_GCOMPUTE_FREQ:
			convert_FWTtoASG_gcompute(thrd.sim,icr,thread_id);
			break;
		case M_DIRECT_TIME:
			fprintf(stderr,"Error: thread_work_FWTtoASG - wrong imethod %d\n",sim->imethod);
			exit(1);
			break;
		case M_DIRECT_FREQ:
			convert_FWTtoASG_direct(thrd.sim,icr,thread_id);
			break;
		}
		i++;
	}
}

void sort_triangle_data(double *frq, TRIANGLE *tri, int *unfold)
{
	double dum;
	int qum;

	if (frq[0] > frq[1]) {
		dum = frq[1]; frq[1] = frq[0]; frq[0] = dum;
		qum = tri->b; tri->b = tri->a; tri->a = qum;
		qum = unfold[1]; unfold[1] = unfold[0]; unfold[0] = qum;
	}
	if (frq[1] > frq[2]) {
		dum = frq[2]; frq[2] = frq[1]; frq[1] = dum;
		qum = tri->c; tri->c = tri->b; tri->b = qum;
		qum = unfold[2]; unfold[2] = unfold[1]; unfold[1] = qum;
	}
	if (frq[0] > frq[1]) {
		dum = frq[1]; frq[1] = frq[0]; frq[0] = dum;
		qum = tri->b; tri->b = tri->a; tri->a = qum;
		qum = unfold[1]; unfold[1] = unfold[0]; unfold[0] = qum;
	}

}

int myround(double d)
{
	return ( d >= 0 ) ? (int)(d+0.5) : (int)(d-0.5);
}

void thread_work_ASG_interpol(int thread_id)
{
	int i, j, k, Ntri, itr, nnz, npts, unfold[3];
	int ia, ib, ic, bin;
	complx *fid, zw, area, prev, height;
	double frq[3], w[3], dum, wbin, fa, fb, fc, freqT;
	TRIANGLE tri;
	Sim_info *sim = thrd.sim;
	Cryst *crdata = NULL;

	if (sim->interpolation == INTERPOL_ASG) { //pure ASG
		crdata = sim->crdata;
	} else { // FWTASG interpolations
		crdata = sim->targetcrdata;
	}
	fid = thrd.fids[thread_id];
	Ntri = LEN(sim->tridata);
	//nnz = sim->FWTASG_nnz;
	if (sim->imethod == M_DIRECT_FREQ) {
		npts = sim->points_per_cycle;
		freqT = 2.0*M_PI*1e6/sim->ASG_period;
		//fprintf(stderr,"Error: ASG thread work - imethod M_DIRECT_FREQ not correctly implemented\n");
		//exit(1);
	} else if (sim->imethod == M_GCOMPUTE_FREQ) {
		npts = sim->ngamma;
		freqT = sim->wr;
	} else {
		fprintf(stderr,"Error: ASG thread work - wrong imethod\n");
		exit(1);
	}
	wbin = sim->sw*2.0*M_PI/sim->np;

	cv_zero(fid);
	i=0;
	while (1) {
		// WARNING!!! This works ONLY IF all MPI slaves have the same number of threads!!!
		itr = thread_id + glob_info.num_threads*glob_info.mpi_rank + glob_info.mpi_size*glob_info.num_threads*i;
		if (itr >= Ntri) break;
		tri = sim->tridata[itr + 1];
		if (sim->interpolation == INTERPOL_ASG) {
			nnz = abs(sim->FWTASG_nnz[tri.a]);
			if (nnz != abs(sim->FWTASG_nnz[tri.b]) || nnz != abs(sim->FWTASG_nnz[tri.c])) {
				fprintf(stderr,"Error: oops, non-equal nnz for triangle %d\nGiving up ASG...\n",itr+1);
				exit(1);
			}
		} else { // FWTASG variants
			nnz = abs(sim->FWTASG_nnz[sim->crmap[tri.a-1]]);
			if (nnz != abs(sim->FWTASG_nnz[sim->crmap[tri.b-1]]) || nnz != abs(sim->FWTASG_nnz[sim->crmap[tri.c-1]])) {
				fprintf(stderr,"Error: oops, non-equal nnz for triangle %d\nGiving up FWTASG...\n",itr+1);
				exit(1);
			}
		}
		for (j=0; j<nnz; j++) {
			//frq[0] = sim->ASG_freq[j+nnz*(tri.a - 1)];
			frq[0] = sim->ASG_freq[tri.a-1][j];
			//frq[1] = sim->ASG_freq[j+nnz*(tri.b - 1)];
			frq[1] = sim->ASG_freq[tri.b - 1][j];
			//frq[2] = sim->ASG_freq[j+nnz*(tri.c - 1)];
			frq[2] = sim->ASG_freq[tri.c - 1][j];
			// resolve possible fold of frequencies
			unfold[0] = unfold[1] = unfold[2] = 0;
			//unfold[1] = myround( (frq[0]-frq[1])/sim->wr );
			//unfold[2] = myround( (frq[0]-frq[2])/sim->wr );
			unfold[1] = myround( (frq[0]-frq[1])/freqT );
			unfold[2] = myround( (frq[0]-frq[2])/freqT );
			//frq[1] += unfold[1]*sim->wr;
			//frq[2] += unfold[2]*sim->wr;
			frq[1] += unfold[1]*freqT;
			frq[2] += unfold[2]*freqT;
			sort_triangle_data(frq,&tri,unfold);
			w[0] = crdata[tri.a].weight;
			w[1] = crdata[tri.b].weight;
			w[2] = crdata[tri.c].weight;
			for (k=0; k<npts; k++) {
				int k0 = k+unfold[0]; while (k0 < 0) k0 += npts; while (k0 >= npts) k0 -= npts;
				int k1 = k+unfold[1]; while (k1 < 0) k1 += npts; while (k1 >= npts) k1 -= npts;
				int k2 = k+unfold[2]; while (k2 < 0) k2 += npts; while (k2 >= npts) k2 -= npts;
				//dum = sim->wr*(k-npts/2+1);
				dum = freqT*(k-npts/2+1);
				//zw.re = w[0]*sim->ASG_ampl[k0+npts*j+(tri.a-1)*npts*nnz].re;
				zw.re = w[0]*sim->ASG_ampl[tri.a-1][k0+npts*j].re;
				//zw.re += w[1]*sim->ASG_ampl[k1+npts*j+(tri.b-1)*npts*nnz].re;
				zw.re += w[1]*sim->ASG_ampl[tri.b-1][k1+npts*j].re;
				//zw.re += w[2]*sim->ASG_ampl[k2+npts*j+(tri.c-1)*npts*nnz].re;
				zw.re += w[2]*sim->ASG_ampl[tri.c-1][k2+npts*j].re;
				//zw.im = w[0]*sim->ASG_ampl[k0+npts*j+(tri.a-1)*npts*nnz].im;
				zw.im = w[0]*sim->ASG_ampl[tri.a-1][k0+npts*j].im;
				//zw.im += w[1]*sim->ASG_ampl[k1+npts*j+(tri.b-1)*npts*nnz].im;
				zw.im += w[1]*sim->ASG_ampl[tri.b-1][k1+npts*j].im;
				//zw.im += w[2]*sim->ASG_ampl[k2+npts*j+(tri.c-1)*npts*nnz].im;
				zw.im += w[2]*sim->ASG_ampl[tri.c-1][k2+npts*j].im;
				zw.re *= 1.0/3.0;
				zw.im *= 1.0/3.0;
				fa = frq[0] + dum; ia = (int)((fa+sim->sw*M_PI)/wbin) + 1;
				fb = frq[1] + dum; ib = (int)((fb+sim->sw*M_PI)/wbin) + 1;
				fc = frq[2] + dum; ic = (int)((fc+sim->sw*M_PI)/wbin) + 1;
				height.re = 2.0*zw.re/(fc-fa); height.im = 2.0*zw.im/(fc-fa);
				//printf("\tk=%i, indices (%i %i %i), weight (%g, %g), height (%g, %g)\n",k,ia,ib,ic,zw.re,zw.im,height.re, height.im);
				// uphill
				prev.re = prev.im = 0.0;
				while (ia<ib) {
					area.re = height.re/(fb-fa)*(-sim->sw*M_PI+ia*wbin-fa)*(-sim->sw*M_PI+ia*wbin-fa)/2.0;
					area.im = height.im/(fb-fa)*(-sim->sw*M_PI+ia*wbin-fa)*(-sim->sw*M_PI+ia*wbin-fa)/2.0;
					bin = ia;
					while (bin < 1) bin += sim->np;
					while (bin > sim->np) bin -= sim->np;
					//printf("uphill first contrib bin %d (%g, %g)\n",bin,area.re,area.im);
					//bin = sim->np - bin + 1;
					bin = sim->conjugate_fid ? (sim->np-bin+1) : (bin);
					fid[bin].re += area.re - prev.re;
					fid[bin].im += area.im - prev.im;
					ia++;
					prev.re = area.re;
					prev.im = area.im;
				}
				bin = ia;
				while (bin < 1) bin += sim->np;
				while (bin > sim->np) bin -= sim->np;
				//printf("uphill last contrib bin %d \n",bin);
				//bin = sim->np - bin + 1;
				bin = sim->conjugate_fid ? (sim->np-bin+1) : (bin);
				fid[bin].re += height.re*(fb-fa)/2.0 - prev.re;
				fid[bin].im += height.im*(fb-fa)/2.0 - prev.im;
				// downhill but from the other side, so it is also uphill...
				prev.re = prev.im = 0.0;
				while (ic > ib) {
					area.re = height.re/(fc-fb)*(fc-(-sim->sw*M_PI+(ic-1)*wbin))*(fc-(-sim->sw*M_PI+(ic-1)*wbin))/2.0;
					area.im = height.im/(fc-fb)*(fc-(-sim->sw*M_PI+(ic-1)*wbin))*(fc-(-sim->sw*M_PI+(ic-1)*wbin))/2.0;
					bin = ic;
					while (bin < 1) bin += sim->np;
					while (bin > sim->np) bin -= sim->np;
					//printf("downhill first contrib bin %d (%g, %g)\n",bin,area.re,area.im);
					//bin = sim->np - bin + 1;
					bin = sim->conjugate_fid ? (sim->np-bin+1) : (bin);
					fid[bin].re += area.re - prev.re;
					fid[bin].im += area.im - prev.im;
					ic--;
					prev.re = area.re;
					prev.im = area.im;
				}
				bin = ic;
				while (bin < 1) bin += sim->np;
				while (bin > sim->np) bin -= sim->np;
				//printf("downhill last contrib bin %d \n",bin);
				//bin = sim->np - bin + 1;
				bin = sim->conjugate_fid ? (sim->np-bin+1) : (bin);
				fid[bin].re += height.re*(fc-fb)/2.0 - prev.re;
				fid[bin].im += height.im*(fc-fb)/2.0 - prev.im;
			}
		}
		i++;
		if (verbose & VERBOSE_PROGRESS) {
			printf("Worker %i/%i: [triangle %d/%d] \n",thread_id+1,glob_info.mpi_rank+1,i,Ntri);
			fflush(stdout);
		}
	}

}


void simpson_thread_slave(void *thr_id)
{
	int thread_id = *(int*)thr_id;

	DEBUGPRINT("Starting simpson_thread_slave %d\n",thread_id);
	//printf("DBGinfo: starting simpson_thread_slave rank %d, thread %d\n",glob_info.mpi_rank, thread_id);

	/* prepare Tcl interpreter */
	Tcl_Interp *interp = Tcl_CreateInterp();
	if (Tcl_Init(interp) == TCL_ERROR) {
		fprintf(stderr,"simpson_thread_slave is unable to initialize Tcl interpreter for process %d, thread %d. Is init.tcl on your path?\n", glob_info.mpi_rank, thread_id);
		fprintf(stderr,"Error: %s\n",Tcl_GetStringResult(interp));
		exit(1);
	}
	TclSlaveInterpInit(interp);
	TclSetSignalHandler(interp,"signalhandler");
	if (Tcl_EvalFile(interp, glob_info.inputfile) != TCL_OK) {
		fprintf(stderr,"Error when evaluating input script '%s' in slave interpreter for process %d, thread %d:\n%s\n",glob_info.inputfile,glob_info.mpi_rank,thread_id,Tcl_GetStringResult(interp));
		exit(1);
	}
	Tcl_Obj *obj = Tcl_NewIntObj(thread_id);
	obj = Tcl_SetVar2Ex(interp,"par","thread_id",obj,TCL_GLOBAL_ONLY);
	if (obj == NULL) {
		fprintf(stderr,"Error in simpson_thread_slave (rank %d, thread %d) when setting 'par(thread_id)':\n%s\n",glob_info.mpi_rank,thread_id,Tcl_GetStringResult(interp));
		exit(1);
	}

	/* wait for job and do it */
	while (1) {
		pthread_barrier_wait(&simpson_b_start);
		/* should thread die? */
		if (glob_info.cont_thread == 0) {
			DEBUGPRINT("simpson_thread_slave %d terminates\n",thread_id);
			break;
		}
		if (glob_info.cont_thread < 3) {
			/* set state of the Tcl interpreter */
			if (Tcl_Eval(interp, thrd.state) != TCL_OK) {
				fprintf(stderr,"Error: can not set slave interpreter state for process %d, thread %d:\n%s\n",glob_info.mpi_rank,thread_id,Tcl_GetStringResult(interp));
				exit(1);
			}
			/* correct for actual MPI rank (state reflects only master) */
			obj = Tcl_NewIntObj(glob_info.mpi_rank);
			obj = Tcl_SetVar2Ex(interp,"par","MPI_rank",obj,TCL_GLOBAL_ONLY);
			if (obj == NULL) {
				fprintf(stderr,"Error in simpson_thread_slave (rank %d, thread %d) when setting 'par(MPI_rank)':\n%s\n",glob_info.mpi_rank,thread_id,Tcl_GetStringResult(interp));
				exit(1);
			}
		}
		switch (glob_info.cont_thread) {
		case 1:
			/* do the ordinary work */
			thread_work(thread_id, interp);
			break;
		case 2:
			// FWTinterpolation - evaluate source orientation set
			thread_work_interpol_source(thread_id, interp);
			break;
		case 3:
			// FWTinterpolation - read data files and interpolate them
			//thread_work_interpol_enlarge(thread_id);
			fprintf(stderr,"Error: glob_info.cont_thread = 3 reserved for NFFT with threads which is not implemented.\n");
			exit(1);
			break;
		case 4:
			// FWTinterpolation - enlarged dataset is ready, construct FID
			thread_work_interpol_calcfid(thread_id);
			break;
		case 5:
			thread_work_ASG_interpol(thread_id);
			break;
		case 6:
			thread_work_FWTtoASG(thread_id);
		}
		/* inform master that job is done by reaching this barrier */
		pthread_barrier_wait(&simpson_b_end);
	}

	/* clean up Tcl interpreter and terminate */
	Tcl_DeleteInterp(interp);
	pthread_exit(NULL);
}
/*****
 * this is called once in order to allocate interpolation data structures in sim
 * ONLY thread MASTER can call it !!!
 * NOT compatible with MPI !!!
 *
void initial_interpol_run(Tcl_Interp *interp, Sim_info *sim) {
	Sim_wsp * wsp;

	wsp = wsp_initialize(sim);
	store_sim_pointers(interp, sim, wsp);
	wsp->interp = interp;
	wsp->thread_id = 0;
	set_averaging_parameters(sim,wsp,thrd.ave_struct,thrd.Navepar,thrd.Naveval);
	wsp->rffac = sim->rfdata[thrd.nrf];
	wsp->zcoor = thrd.zvals[thrd.nz];
	set_inhom_offsets(sim,wsp,thrd.zoffsetvals[thrd.nz]);
	wsp->cryst = sim->crdata[sim->icr_done];
	wsp->cryst_idx = sim->icr_done;
	sim_calcfid_interpol(sim, wsp);
	if (verbose & VERBOSE_PROGRESS) {
		printf("Worker %i/%i: [cryst %d/%d] for interpolation\n",1,glob_info.mpi_rank+1,sim->icr_done,thrd.ncr);
		fflush(stdout);
	}
	wsp_destroy(sim, wsp);
	free(wsp);
}
****/

complx * simpson_common_work(Tcl_Interp *interp, char *state, int *ints_out, double *dbls_out)
{
	int i, ncr, nrf, nz, Navepar, Naveval, Ntot, num_threads;
	double *zvals=NULL, *zoffsetvals=NULL, *aveweight=NULL;
	Ave_elem *avestruct=NULL;
	Sim_info *sim;
	complx *fidsum=NULL;

	/* read global simulation information */
	sim = sim_initialize(interp);
	/* read all averaging profiles (crystallites, rf inhom, z-averaging)*/
	prepare_zaveraging(interp,&zvals,&zoffsetvals);
	read_averaging_file(interp, sim, &avestruct, &aveweight, &Navepar, &Naveval);
	//printf("\n\ntest ave: %d, %lg\n\n",avestruct[1].par_type,aveweight[2]);
	ncr = LEN(sim->crdata);
	nrf = rfprof_len(sim->rfdata);
	nz = LEN(zvals);
	Ntot = ncr*nrf*nz*Naveval; /* Total number of jobs */
	num_threads = glob_info.num_threads;

	thrd.fids = (complx**)malloc(num_threads*sizeof(complx*));
	if (thrd.fids == NULL) {
		fprintf(stderr,"Error: ups, not even started and already out of memory...\n");
		exit(1);
	}
	for (i=0; i<num_threads; i++) thrd.fids[i] = complx_vector(sim->ntot);

	switch (sim->interpolation) {
	case INTERPOL_NOT_USED: // no interpolation
		/* prepare input data for threads */
		thrd.Ntot = Ntot;
		thrd.ncr = ncr;
		thrd.nrf = nrf;
		thrd.nz = nz;
		thrd.Nacq = 0;
		thrd.phase = 0;
		thrd.zvals = zvals;
		thrd.zoffsetvals = zoffsetvals;
		thrd.state = state;
		thrd.Navepar = Navepar;
		thrd.Naveval = Naveval;
		thrd.ave_struct = avestruct;
		thrd.ave_weight = aveweight;
		thrd.sim = sim;
		/* start calculation in threads */
		glob_info.cont_thread = 1;
		pthread_barrier_wait(&simpson_b_start);
		/* harvest results from threads */
		pthread_barrier_wait(&simpson_b_end);
		fidsum = complx_vector(sim->ntot);
		cv_zero(fidsum);
		for(i=0; i<num_threads; i++) {
			cv_multod(fidsum, thrd.fids[i], 1.0);
			//free_complx_vector(thrd.fids[i]);
		}
		break;
	case INTERPOL_FWT_ALL:
	case INTERPOL_FWT_LAM: { // FWTinterpolation
		int iz, iave, irf;
		double weight;
		fidsum = complx_vector(sim->ntot);
		cv_zero(fidsum);
		sim_prepare_interpol(sim);
		// all averaging except powder needs to be serial
		for (iz=1; iz<=nz; iz++) {
			for (iave=0; iave<Naveval; iave++) {
				for (irf=1; irf<=nrf; irf++) {
					/* prepare input data for threads */
					thrd.Ntot = Ntot;
					thrd.ncr = ncr;
					thrd.nrf = irf;
					thrd.nz = iz;
					thrd.Nacq = 0;
					thrd.phase = 0;
					thrd.zvals = zvals;
					thrd.zoffsetvals = zoffsetvals;
					thrd.state = state;
					thrd.Navepar = Navepar;
					thrd.Naveval = iave;
					thrd.ave_struct = avestruct;
					thrd.ave_weight = aveweight;
					thrd.sim = sim;
					/* start calculation in threads */
					glob_info.cont_thread = 2;
					pthread_barrier_wait(&simpson_b_start);
					/* wait for completion */
					pthread_barrier_wait(&simpson_b_end);
#ifdef MPI
					// wait somehow for all MPI slaves to finish, only then continue
					MPI_Barrier(MPI_COMM_WORLD);
#endif
					// read data files and interpolate, use ONLY MPI workers
					// NFFT library does not work properly on windows with threads
					printf("\nInterpolating... \n");
					master_FWTinterpolate(sim);
#ifdef MPI
					// wait somehow for all MPI slaves to finish, only then continue
					MPI_Barrier(MPI_COMM_WORLD);
#endif
					printf(" done.\n");
					// calculate fid
					glob_info.cont_thread = 4;
					pthread_barrier_wait(&simpson_b_start);
					/* wait for completion */
					pthread_barrier_wait(&simpson_b_end);
					printf("\nfids generated ... ");
					weight = sim->rfdata[irf][sim->ss->nchan+1] * aveweight[iave];
					for(i=0; i<num_threads; i++){
						cv_multod(fidsum, thrd.fids[i], weight);
						//free_complx_vector(thrd.fids[i]);
					}
					printf("and collected.\n");
					sim_preempty_interpol(sim);
					// repeat for averaging
				}
			}
		}
		//nfsft_forget(); // free NFFT pre-calculated globals
		break; }
	case INTERPOL_ASG : { // ASGinterpolation
		int iz, iave, irf;
		double weight;
		fidsum = complx_vector(sim->ntot);
		cv_zero(fidsum);
		sim_prepare_interpol(sim);
		// initialize nfft library
		//nfsft_precompute(sim->Jinterpol[1],1000.0,0U,0U);
		// all averaging except powder needs to be serial
		for (iz=1; iz<=nz; iz++) {
			for (iave=0; iave<Naveval; iave++) {
				for (irf=1; irf<=nrf; irf++) {
					/* prepare input data for threads */
					thrd.Ntot = Ntot;
					thrd.ncr = ncr;
					thrd.nrf = irf;
					thrd.nz = iz;
					thrd.Nacq = 0;
					thrd.phase = 0;
					thrd.zvals = zvals;
					thrd.zoffsetvals = zoffsetvals;
					thrd.state = state;
					thrd.Navepar = Navepar;
					thrd.Naveval = iave;
					thrd.ave_struct = avestruct;
					thrd.ave_weight = aveweight;
					thrd.sim = sim;
					/* start calculation in threads */
					glob_info.cont_thread = 2;
					pthread_barrier_wait(&simpson_b_start);
					/* wait for completion */
					pthread_barrier_wait(&simpson_b_end);
#ifdef MPI
					// wait somehow for all MPI slaves to finish, only then continue
					MPI_Barrier(MPI_COMM_WORLD);
#endif
					// read all data to memory
					mpi_work_ASGread(sim);
					// calculate subspectra in threads
					glob_info.cont_thread = 5;
					pthread_barrier_wait(&simpson_b_start);
					/* wait for completion */
					pthread_barrier_wait(&simpson_b_end);
					//printf("\nASG subspectra generated ... ");
					weight = sim->rfdata[irf][sim->ss->nchan+1] * aveweight[iave];
					for(i=0; i<num_threads; i++){
						cv_multod(fidsum, thrd.fids[i], weight);
						//free_complx_vector(thrd.fids[i]);
					}
					//printf("and collected.\n");
					sim_preempty_interpol(sim);
					// repeat for averaging
				}
			}
		}
		break; }
	case INTERPOL_FWTASG_ALL:
	case INTERPOL_FWTASG_LAM: {// FWTASGinterpolation
		int iz, iave, irf;
		double weight;
		fidsum = complx_vector(sim->ntot);
		cv_zero(fidsum);
		sim_prepare_interpol(sim);
		// all averaging except powder needs to be serial
		for (iz=1; iz<=nz; iz++) {
			for (iave=0; iave<Naveval; iave++) {
				for (irf=1; irf<=nrf; irf++) {
					/* prepare input data for threads */
					thrd.Ntot = Ntot;
					thrd.ncr = ncr;
					thrd.nrf = irf;
					thrd.nz = iz;
					thrd.Nacq = 0;
					thrd.phase = 0;
					thrd.zvals = zvals;
					thrd.zoffsetvals = zoffsetvals;
					thrd.state = state;
					thrd.Navepar = Navepar;
					thrd.Naveval = iave;
					thrd.ave_struct = avestruct;
					thrd.ave_weight = aveweight;
					thrd.sim = sim;
					/* start FWT calculation in threads first */
					glob_info.cont_thread = 2;
					pthread_barrier_wait(&simpson_b_start);
					/* wait for completion */
					pthread_barrier_wait(&simpson_b_end);
#ifdef MPI
					// wait somehow for all MPI slaves to finish, only then continue
					MPI_Barrier(MPI_COMM_WORLD);
#endif
					// read data files and interpolate, use ONLY MPI workers
					// NFFT library does not work properly on windows with threads
//					printf("\nInterpolating with FWT ");
					master_FWTinterpolate(sim);
#ifdef MPI
					// wait somehow for all MPI slaves to finish, only then continue
					MPI_Barrier(MPI_COMM_WORLD);
#endif
//					printf(" done\n");
					// the commented below should be done in FWTinterpolate()!!!
					// prepare data for ASG interpolation
					//if (sim->ASG_freq == NULL) {
					//	sim->ASG_freq = (double*)malloc(LEN(sim->targetcrdata)*sim->FWTASG_nnz*sizeof(double));
					//	sim->ASG_ampl = (complx*)malloc(LEN(sim->targetcrdata)*sim->FWTASG_nnz*sim->points_per_cycle*sizeof(complx));
					//	if ( sim->ASG_freq == NULL || sim->ASG_ampl == NULL ) {
					//		fprintf(stderr,"Error: thread master - no more memory for ASG data");
					//		exit(1);
					//	}
					//}
					//assert(sim->ASG_freq != NULL && sim->ASG_ampl != NULL);
					glob_info.cont_thread = 6;
					pthread_barrier_wait(&simpson_b_start);
					/* wait for completion */
					pthread_barrier_wait(&simpson_b_end);
//					printf("\nASG data generated \n");
					// read all data to memory
					mpi_work_ASGread(sim);
					// calculate fid
					glob_info.cont_thread = 5;
					pthread_barrier_wait(&simpson_b_start);
					/* wait for completion */
					pthread_barrier_wait(&simpson_b_end);
//					printf("\nfids generated ... ");
					weight = sim->rfdata[irf][sim->ss->nchan+1] * aveweight[iave];
					for(i=0; i<num_threads; i++){
						cv_multod(fidsum, thrd.fids[i], weight);
						//free_complx_vector(thrd.fids[i]);
					}
//					printf("and collected.\n");
					sim_preempty_interpol(sim);
					// repeat for averaging
				}
			}
		}
		break; }
	default:
		fprintf(stderr,"Error: invalid value for sim->interpolation (%d)\n",sim->interpolation);
		exit(1);
	} // end switch (interpolation)

	/* set output parameters */
	if (ints_out != NULL) {
		ints_out[0] = sim->np;
		ints_out[1] = sim->ni;
		ints_out[2] = sim->domain;
	}
	if (dbls_out != NULL) {
		dbls_out[0] = sim->sw;
		dbls_out[1] = sim->sw1;
		dbls_out[2] = (double)nz * rfprof_sumweight(sim->rfdata);
	}
	sim_destroy(sim,0);
	free(sim);
	free_double_vector(zvals);
	free_double_vector(zoffsetvals);
	free_averaging_data(&avestruct,Navepar);
	free(aveweight); aveweight = NULL;
	for (i=0; i<num_threads; i++) free_complx_vector(thrd.fids[i]);
	free(thrd.fids);

	return fidsum;
}


void simpson_mpi_slave(Tcl_Interp *interp)
{
#ifdef MPI
	// MPI slave waits here for all possible calls of fsimpson
	MPI_Status status;
	int i, j, k, Nints, Ndbls, Nstate, Nsh, terminate = 0;
	int *intmessage = NULL;
	char *state = NULL;
	double *dblmessage = NULL;
	complx *fidsum;
	int mpirank = glob_info.mpi_rank;
	int mpisize = glob_info.mpi_size;

	while (1) {
		// wait for message whether to terminate or continue (tag = 1)
		MPI_Recv(&terminate, 1, MPI_INT, 0, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		if (terminate) {
			break;
		}

		// wait for message from master - integer data (tag = 10)
		// Probe for an incoming message from process zero (master)
		MPI_Probe(0, 10, MPI_COMM_WORLD, &status);
		// When probe returns, the status object has the size and other
		// attributes of the incoming message. Get the size of the message
		MPI_Get_count(&status, MPI_INT, &Nints);
		// Allocate a buffer just big enough to hold the incoming data
		intmessage = (int*)malloc(sizeof(int)*Nints);
		MPI_Recv(intmessage, Nints, MPI_INT, 0, 10, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		assert(Nints >= 4);
		// process this message
		if (OCpar.isinit) OCpar_destroy();
		OCpar.isinit = intmessage[0];
		j = 1;
		if (OCpar.isinit != 0) {
			OCpar.gradmode = intmessage[j++];
			OCpar.gradmodeprop = intmessage[j++];
			if (intmessage[j] > 0) {
				//if (OCpar.grad_shapes != NULL) free_int_vector(OCpar.grad_shapes);
				OCpar.grad_shapes = int_vector(intmessage[j++]);
				for (i=1; i<=OCpar.grad_shapes[0]; i++) OCpar.grad_shapes[i] = intmessage[j++];
			} else {
				//if (OCpar.grad_shapes != NULL) free_int_vector(OCpar.grad_shapes);
				//OCpar.grad_shapes = NULL;
				j++;
			}
		}
		Nsh = intmessage[j++];
		Nstate = intmessage[Nints-2];
		Ndbls = intmessage[Nints-1];
		// Now receive the state message (tag = 20);
		state = (char*)malloc(sizeof(char)*Nstate);
		MPI_Recv(state, Nstate, MPI_CHAR, 0, 20, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		// receive eventual dbl message
		RFshapes_reset();
		if (Ndbls > 0) {
			dblmessage = (double*)malloc(sizeof(double)*Ndbls);
			MPI_Recv(dblmessage, Ndbls, MPI_DOUBLE, 0, 30, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			int slot, len, pos = 0;
			for (i=0; i<Nsh; i++) {
				slot = intmessage[j++];
				len = intmessage[j++];
				RFshapes[slot] = RFshapes_alloc(len);
				for (k=1; k<=len; k++) {
					RFshapes[slot][k].ampl  = dblmessage[pos++];
					RFshapes[slot][k].phase = dblmessage[pos++];
				}
			}
			free(dblmessage);
			dblmessage = NULL;
		}
		free(intmessage);

		if (Tcl_Eval(interp, state) != TCL_OK) {
			fprintf(stderr,"Error: can not set MPI slave interpreter state for process %d:\n%s\n",glob_info.mpi_rank,Tcl_GetStringResult(interp));
			exit(1);
		}
		char buf[256];
		TclGetString(interp,buf,"par","verbose",0,"0");
		getbits(&verbose,buf);

		// do the job
	//printf("DBGinfo: MPIslave %d state:\n%s\n---------\n",glob_info.mpi_rank,state);
		fidsum = simpson_common_work(interp, state, NULL, NULL);
		// send data to master
		MPI_Send(fidsum, 2*(LEN(fidsum)+1), MPI_DOUBLE, 0, 0,MPI_COMM_WORLD);
		free_complx_vector(fidsum);
		free(state);
	}
#endif
}

FD* simpson(Tcl_Interp* interp, int mpirank, int mpisize)
{
	// only MPI master should get here
	assert(mpirank == 0);

	int npnidm[3];
	double swsinten[3];
	complx *fidsum;
	FD *fd;

	Tcl_ResetResult(interp);
	if (Tcl_Eval(interp,"savestate") != TCL_OK) {
		fprintf(stderr,"Error: can not get main interpreter state:\n%s\n",Tcl_GetStringResult(interp));
		exit(1);
	}
	char *state = strdup(Tcl_GetStringResult(interp));
	//printf("\n\nSAVESTATE result:\n%s\n\n",state);
#ifdef MPI
	/* prepare messages for slaves */
	int i, j, k, Nints = 1, Ndbls = 0, Nc = 0, terminate = 0;
	int *intmessage;
	double *dblmessage;
	if (OCpar.isinit != 0) {
		Nints += 3;
		if (OCpar.grad_shapes != NULL) Nints += OCpar.grad_shapes[0];
	}
	for (i=0; i<MAXRFSHAPES; i++) {
		if (RFshapes[i] != NULL) Nc++;
	}
	Nints += 1+Nc*2 + 2;
	intmessage = (int*)malloc(Nints*sizeof(int));
	intmessage[0] = OCpar.isinit;
	j = 1;
	if (OCpar.isinit != 0) {
		intmessage[j++] = OCpar.gradmode;
		intmessage[j++] = OCpar.gradmodeprop;
		if (OCpar.grad_shapes != NULL) {
			for (i=0; i<= OCpar.grad_shapes[0]; i++) intmessage[j++] = OCpar.grad_shapes[i];
		} else {
			intmessage[j++] = 0;
		}
	}
	intmessage[j++] = Nc;
	Ndbls = 0;
	if (Nc != 0) {
		for (i=0; i<MAXRFSHAPES; i++) {
			if (RFshapes[i] != NULL) {
				intmessage[j++] = i;
				intmessage[j++] = RFshapes_len(i);
				Ndbls += 2*RFshapes_len(i);
			}
		}
	}
	intmessage[j++] = strlen(state)+1; // size of state message
	intmessage[j++] = Ndbls; // size of dblmessage
	assert(j == Nints);
	// intmessage complete, now create dblmessage
	if (Ndbls != 0) {
		dblmessage = (double*)malloc(Ndbls*sizeof(double));
		int pos = 0;
		for (i=0; i<MAXRFSHAPES; i++) {
			if (RFshapes[i] != NULL) {
				for (k=1; k<=RFshapes_len(i); k++) {
					dblmessage[pos++] = RFshapes[i][k].ampl;
					dblmessage[pos++] = RFshapes[i][k].phase;
				}
			}
		}
		assert(pos == Ndbls);
	} else {
		dblmessage = NULL;
	}
	Nc = strlen(state)+1;
	/* send messages out */
	for (i=1; i<mpisize; i++) {
		MPI_Send(&terminate,1,MPI_INT,i,1,MPI_COMM_WORLD);  // inform slave we will work
		MPI_Send(intmessage,Nints,MPI_INT,i,10,MPI_COMM_WORLD); // give him intmessage (tag = 10);
		MPI_Send(state,Nc,MPI_CHAR,i,20,MPI_COMM_WORLD); // give him state (tag = 20)
		if (Ndbls != 0) MPI_Send(dblmessage,Ndbls,MPI_DOUBLE,i,30,MPI_COMM_WORLD); // give him dblmessage if any (tag = 30)
	}
	free(intmessage);
	free(dblmessage);
	// state will be used later
#endif

	fidsum = simpson_common_work(interp, state, npnidm, swsinten);

	/* wait for MPI slaves to send their results */
#ifdef MPI
	complx *buffer = complx_vector(LEN(fidsum));
	for (i=1; i<mpisize; i++) {
		MPI_Recv(buffer, 2*(LEN(fidsum)+1), MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		/* update global vector with results from each local vector */
		cv_multod(fidsum, buffer, 1.0);
	}
	free_complx_vector(buffer);
#endif

	/* Master: global vector is now complete */
	//totalintensity = cryst_sumweight(crdata)*rfprof_sumweight(rfdata)*(double)nz;
	cv_muld(fidsum,1.0/swsinten[2]);
	if (verbose & VERBOSE_PROGRESS) printf("All done!\n");

	fd = FD_data2fd(NULL,fidsum,npnidm[0],npnidm[1],swsinten[0],swsinten[1]);
	if (npnidm[2] == 1) fd->type = FD_TYPE_SPE;

	free(state);
	free_complx_vector(fidsum);

	return fd;
}

int tclCrystallites(ClientData data,Tcl_Interp* interp,
      int argc, Tcl_Obj *argv[])
{
  int i, *map;
  char name1[256], name2[256];
  Cryst *cr1=NULL, *cr2=NULL;
  TRIANGLE *tria;
  int savebin = 0;
  int what = 0;

  for (i=2; i<argc; i++) {
	  strcpy(name1, Tcl_GetString(argv[i]));
	  if (!strncmp(name1,"-savebin",8)) {
		  savebin = 1;
	  } else if (!strncmp(name1,"-tri",4)) {
		  what = 1;
	  } else if (!strncmp(name1,"-map",4)) {
		  what = 2;
	  } else { // this seems to be a crystallite set
		  if (cr2 != NULL) {
			  Tcl_SetResult(interp,"Error in crystallites: confusion in parameters",TCL_VOLATILE);
			  return TCL_ERROR;
		  }
		  cr2 = read_crystfile(name1, 1, -1);
		  strcpy(name2, name1);
	  }
  }

  strcpy(name1, Tcl_GetString(argv[1]));
  //printf("what = %d, name1 = '%s', savebin = %d\n",what,name1,savebin);
  switch (what) {
  case 0: // work on crystallites
	  cr1 = read_crystfile(name1, 1, -1);
	  if (savebin == 0) {
		  Tcl_ResetResult(interp);
		  for (i=1;i<=LEN(cr1);i++) {
			  TclAppendResult(interp,"%g %g %g",cr1[i].alpha,cr1[i].beta,cr1[i].weight);
		  }
	  } else {
		  save_bin_crystfile(name1,cr1);
	  }
	  cryst_free(cr1);
	  break;
  case 1: // work on triangles
	  tria = read_triangle_file(name1);
	  if (savebin == 0) {
		  for (i=1;i<=LEN(tria);i++) {
			  TclAppendResult(interp,"%d %d %d",tria[i].a,tria[i].b,tria[i].c);
		  }
	  } else {
		  save_bin_triangle_file(name1,tria);
	  }
	  triangle_free(tria);
	  break;
  case 2: // work on crystallite sets map
	  cr1 = read_crystfile(name1, 1, -1);
	  if (cr1 == NULL || cr2 == NULL) {
		  Tcl_SetResult(interp,"Error in crystallites: need two crystallite sets when -map set",TCL_VOLATILE);
		  return TCL_ERROR;
	  }
	  map = read_cryst_map(name1, cr1, name2, cr2);
	  if (savebin == 0) {
		  for (i=1;i<=LEN(cr2);i++) {
			  TclAppendResult(interp,"%d %d",i,map[i-1]);
		  }
	  } else {
		  save_bin_cryst_map(name1, name2, map, LEN(cr2));
	  }
	  cryst_free(cr1);
	  cryst_free(cr2);
	  break;
  }
  return TCL_OK;

	/***
	   sprintf(name, "%s_cryst", Tcl_GetString(argv[1]));
	  while (strlen(cryst_names[n])) {
		  if (!strcmp(cryst_names[n],name)) {
			  c = cryst_pointers[n];
			  for (i=0;i<cryst_numbers[n];i++) {
				  TclAppendResult(interp,"%g %g %g",c[i].alpha,c[i].beta,c[i].weight);
			  }
			  return TCL_OK;
		  }
		  n++;
	  }
	  return TclError(interp,"crystallites: unable to find internal crystal file '%s'", argv[1]);
	***/
}

int tclInternalSimpson(ClientData data,Tcl_Interp* interp,
      int argc, Tcl_Obj *argv[]){
  char buf[256],buf1[16];
  int fnew(FD* f);
  int fidN, mpirank, mpisize;
  extern FD **fd;
  FD* f;

  if (argc != 1) {
    return TclError(interp, "Usage: <desc> internalsimpson");
  }

  TclGetString(interp,buf,"par","verbose",0,"0");
  if (!getbits(&verbose,buf))
    return TclError(interp,"error: 'verbose' parameter contains illegal characters '%s' (must be a string of '0' and '1')",buf);

  TclGetString(interp,buf,"par","various",0,"0");
  if (!getbits(&various,buf))
    return TclError(interp,"error: 'various' parameter contains illegal characters '%s' (must be a string of '0' and '1')",buf);
  mpirank = TclGetInt(interp,"par","MPI_rank",0,0);
  mpisize = TclGetInt(interp,"par","MPI_size",0,1);
  assert(mpirank == 0);

  f = simpson(interp,mpirank,mpisize);
  
  fidN = fnew(f);
  /*  fprintf(stderr,"fidN = %d\n", fidN); */
  sprintf(buf1,"%d", fidN);
  sprintf(buf,"%ld", (long)fd[fidN]);
  Tcl_Eval(interp,"namespace eval FD {variable f}");
  Tcl_SetVar2(interp,"FD::f",buf1,buf,TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG);
  Tcl_SetVar2(interp,"FD_Internal",buf1,buf,TCL_GLOBAL_ONLY);

  return TclSetResult(interp,"%d",fidN);
}

/****** remove this ********
int tclnfft_test(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[]) {

	printf("\nKUKUKUKUKUKIU\n\n");

	simpson_nfft_test();

	return TCL_OK;
}

extern void simpson_fftw_test(void);

int tclfftw_test(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[]) {


	printf("\nExternal NFFT test\n++++++++++++++++++++++++\n\n");
	system("D:\\Tosner\\workspace\\NFFTW_test\\Debug\\NFFTW_test.exe");


	printf("\nINTERNAL  FFTW TEST\n\n");

	simpson_fftw_test();

	return TCL_OK;
}
*************/


void tclcmd_simpson(Tcl_Interp* interp)
{
  Tcl_CreateObjCommand(interp,"internalsimpson",(Tcl_ObjCmdProc *)tclInternalSimpson,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"crystallites",(Tcl_ObjCmdProc *)tclCrystallites,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateObjCommand(interp,"nfft_test",(Tcl_ObjCmdProc *)tclnfft_test,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateObjCommand(interp,"fftw_test",(Tcl_ObjCmdProc *)tclfftw_test,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

} 
