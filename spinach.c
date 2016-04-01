/*
 * spinach.c
 *
 *  Created on: Aug 30, 2010
 *      Author: zdenek
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <tcl.h>
#include "complx.h"
#include "matrix.h"
#include "cm.h"
#include "blockdiag.h"
#include "ham.h"
#include "tclutil.h"
#include "pulse.h"
#include "spinach.h"
#include "defs.h"
#include "sim.h"

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



/* works on parts of sparse matrix representation: n=A->row/col, rp=A->irow, ci=A->icol */
int * scomponents(int n, int *rp, int *ci)
{
	int *sci, *root, *dt, *cs, *rs, cn, t, css, rss;
	int sv, v, ri, w;

	sci = int_vector(n); iv_zero(sci);
	cn = 1;
	root = int_vector(n); iv_zero(root);
	dt = int_vector(n); iv_zero(dt);
	t = 0;
	cs = int_vector(n); iv_zero(cs);
	css = 0;
	rs = int_vector(2*n); iv_zero(rs);
	rss = 0;

	for (sv=1; sv<=n; sv++) {
		v = sv;
		if (root[v]>0) continue;
		rss++;
		rs[2*rss-1] = v;
		rs[2*rss] = rp[v-1];
		root[v] = v;
		sci[v] = -1;
		dt[v] = t;
		t++;
		css++;
		cs[css] = v;
		while (rss>0) {
			v = rs[2*rss-1];
			ri = rs[2*rss];
			rss--;
			while (ri<rp[v]) {
				w = ci[ri-1];
				ri++;
				if (root[w]==0) {
					root[w] = w;
					sci[w] = -1;
					dt[w] = t;
					t++;
					css++;
					cs[css] = w;
					rss++;
					rs[2*rss-1] = v;
					rs[2*rss] = ri;
					v = w;
					ri = rp[w-1];
					continue;
				}
			}
			for (ri=rp[v-1]; ri<=rp[v]-1; ri++) {
				w = ci[ri-1];
				if (sci[w]==-1) {
					if (dt[root[v]]>dt[root[w]]) {
						root[v] = root[w];
					}
				}
			}
			if (root[v] == v) {
				while (css>0) {
					w = cs[css];
					css--;
					sci[w] = cn;
					if (w==v) break;
				}
				cn++;
			}
		}
	}
	free_int_vector(root);
	free_int_vector(dt);
	free_int_vector(cs);
	free_int_vector(rs);
	return sci;
}


mat_double * liou_from_ham(blk_mat_double *ham)
{
	int i, j, k, r, c, N, NN, Nmax;
	double val;
	mat_double *liou;

	N = ham->dim;
	NN = N*N;
	Nmax = 0;
	for (i=0; i<ham->Nblocks; i++) {
		if (ham->blk_dims[i] == 1) {
			if (fabs(ham->m[i].data[0]) > SPARSE_TOL) Nmax++;
		} else {
			Nmax += dm_nnz(ham->m+i);
		}
	}
	Nmax *= 2*N;
	liou = double_matrix(NN,NN,MAT_SPARSE,Nmax,0);

	liou->irow[0] = 1;
	k = 0;
	for (i=1; i<=NN; i++) {
		liou->irow[i] = liou->irow[i-1];
		for (j=1; j<=NN; j++) {
			val = 0.0;
			if (i%N == j%N) {
				r = (j-1)/N+1;
				c = (i-1)/N+1;
				val += blk_dm_getelem(ham,r,c);
			}
			if ((i-1)/N == (j-1)/N) {
				r = i%N; if (r == 0) r = N;
				c = j%N; if (c == 0) c = N;
				val -= blk_dm_getelem(ham,r,c);
			}
			if (fabs(val) >= SPARSE_TOL) {
				(liou->irow[i])++;
				liou->icol[k] = j;
				liou->data[k] = val;
				if (++k > Nmax) {
					fprintf(stderr,"Error: liou_from_ham - Nmax=%d overflow after (%d,%d)\n",Nmax,i,j);
					exit(1);
				}
			}
		}
	}
	DEBUGPRINT("liou_from_ham: nnz=%d (=%d), initially allocated %d\n",liou->irow[NN]-1,k,Nmax);
	if (Nmax > k) {
		Nmax = k;
		double *new_data = (double*)realloc(liou->data, Nmax*sizeof(double));
		MKL_INT *new_icol = (MKL_INT*)realloc(liou->icol, Nmax*sizeof(MKL_INT));
		assert( (new_data != NULL) && (new_icol != NULL) );
		liou->data = new_data;
		liou->icol = new_icol;
	}

	return liou;
}

complx * get_sigma(mat_complx *mat, int adjoint)
{
	complx *vec, *vp, *mp;
	int i, j, *ic;

	vec = complx_vector(mat->row*mat->col);

	switch (mat->type) {
	case MAT_DENSE:
		if (adjoint) {
			cv_zero(vec);
			for (i=0; i<mat->row; i++) {
				cblas_dcopy(mat->col,(double*)(mat->data+i),2*mat->row,(double*)(vec+1+i*mat->col),2);
				cblas_daxpy(mat->col,-1.0,(double*)(mat->data+i)+1,2*mat->row,(double*)(vec+1+i*mat->col)+1,2);
			}
		} else {
			memcpy(&(vec[1]),mat->data,mat->row*mat->col*sizeof(complx));
		}
		break;
	case MAT_DENSE_DIAG:
		cv_zero(vec);
		vp = vec+1; mp = mat->data;
		for (i=0; i<mat->row; i++) {
			if (adjoint) {
				vp->re = mp->re;
				vp->im = -mp->im;
			} else {
				*vp = *mp;
			}
			vp += (mat->row + 1);
			mp++;
		}
		break;
	case MAT_SPARSE:
	case MAT_SPARSE_DIAG:
		cv_zero(vec);
		ic = mat->icol;
		mp = mat->data;
		for (i=1; i<=mat->row; i++) {
			for (j=0; j<mat->irow[i]-mat->irow[i-1]; j++) {
				if (adjoint) {
					vec[*ic-1+(i-1)*mat->row].re = mp->re;
					vec[*ic-1+(i-1)*mat->row].im = -mp->im;
				} else {
					vec[i-1+(*ic-1)*mat->row] = *mp;
				}
				ic++;
				mp++;
			}
		}
		break;
	default :
		fprintf(stderr,"Error: get_sigma - unknown matrix type\n");
		exit(1);
	}
	return vec;
}

complx * get_rho(mat_complx *mat, int *idx, int adjoint)
{
	complx *rho, z;
	int i, row, col, dim, len;

	len = LEN(idx);
	dim = mat->row;
	rho = complx_vector(len);
	for (i=1; i<=len; i++) {
		row = idx[i] % dim; if (row == 0) row = dim;
		col = (idx[i] - 1)/dim + 1;
		if (adjoint) {
			z = cm_getelem(mat,col,row);
			rho[i].re = z.re;
			rho[i].im = -z.im;
		} else {
			rho[i] = cm_getelem(mat,row,col);
		}
	}

	return rho;
}

/* creates restricted current Liouvillian as dense matrix */
void ham_liou_restr(blk_mat_double *ham, int *idx, mat_double *Lr)
{
	int i, j, r, c, k, l, Nr, N;
	double val;

	Nr = Lr->row;
	N = ham->dim;
	dm_zero(Lr);
	for (i=1; i<=Nr; i++) {
		r = idx[i];
		for (j=1; j<=Nr; j++) {
			c = idx[j];
			//printf("index %d, %d ",r,c);
			if (r%N == c%N) {
				k = (r-1)/N +1;
				l = (c-1)/N +1;
				val = blk_dm_getelem(ham,k,l);
				//printf(" (1) n=%d, k=%d\n",k,l);
			} else {
				val = 0.0;
			}
			if ( (r-1)/N == (c-1)/N) {
				k = r%N; if (k==0) k=N;
				l = c%N; if (l==0) l=N;
				val -= blk_dm_getelem(ham,k,l);
				//printf(" (2) n=%d, k=%d\n",k,l);
			}
			Lr->data[i-1 + (j-1)*Nr] = val;
		}
	}
}

/* this returns Liouvillian as sparse matrix from complex Hamiltonian */
mat_complx * create_liou(blk_mat_complx *ham, int *idx)
{
	mat_complx *L;
	int i, j, k, r, c, N, NN, Nmax, ii, jj;
	complx z, zval;

	N = ham->dim;
	Nmax = 0;
	for (i=0; i<ham->Nblocks; i++) {
		if (ham->blk_dims[i] == 1) {
			zval = ham->m[i].data[0];
			if (fabs(zval.re)>SPARSE_TOL || fabs(zval.im)>SPARSE_TOL) Nmax++;
		} else {
			Nmax += cm_nnz(ham->m+i);
		}
	}
	DEBUGPRINT("create_liou: hamiltonian non-zeros = %d\n",Nmax);
	Nmax *= 2*N;

	if (idx == NULL) {
		NN = N*N;
		L = complx_matrix(NN,NN,MAT_SPARSE,Nmax,ham->basis);

		L->irow[0] = 1;
		k = 0;
		for (i=1; i<=NN; i++) {
			L->irow[i] = L->irow[i-1];
			for (j=1; j<=NN; j++) {
				zval = Cnull;
				if (i%N == j%N) {
					r = (j-1)/N+1;
					c = (i-1)/N+1;
					z = blk_cm_getelem(ham,r,c);
					zval.re += z.re;
					zval.im += z.im;
				}
				if ((i-1)/N == (j-1)/N) {
					r = i%N; if (r == 0) r = N;
					c = j%N; if (c == 0) c = N;
					z = blk_cm_getelem(ham,r,c);
					zval.re -= z.re;
					zval.im -= z.im;
				}
				if (Cnorm(zval) >= SPARSE_TOL) {
					(L->irow[i])++;
					L->icol[k] = j;
					L->data[k] = zval;
					if (++k > Nmax) {
						fprintf(stderr,"Error: create_liou - Nmax=%d overflow after (%d,%d)\n",Nmax,i,j);
						exit(1);
					}
				}
			}
		}
		DEBUGPRINT("create_liou: nnz=%d (=%d), initially allocated %d\n",L->irow[NN]-1,k,Nmax);
		if (Nmax > k) {
			Nmax = k;
			complx *new_data = (complx*)realloc(L->data, Nmax*sizeof(complx));
			MKL_INT *new_icol = (MKL_INT*)realloc(L->icol, Nmax*sizeof(MKL_INT));
			assert( (new_data != NULL) && (new_icol != NULL) );
			L->data = new_data;
			L->icol = new_icol;
		}

	} else {
		NN = LEN(idx);
		Nmax = (Nmax > NN*NN) ? NN*NN : Nmax;
		L = complx_matrix(NN,NN,MAT_SPARSE,Nmax,ham->basis);
		L->irow[0] = 1;
		k = 0;
		for (ii=1; ii<=NN; ii++) {
			L->irow[ii] = L->irow[ii-1];
			i = idx[ii];
			for (jj=1; jj<=NN; jj++) {
				j = idx[jj];
				zval = Cnull;
				if (i%N == j%N) {
					r = (j-1)/N+1;
					c = (i-1)/N+1;
					z = blk_cm_getelem(ham,r,c);
					zval.re += z.re;
					zval.im += z.im;
				}
				if ((i-1)/N == (j-1)/N) {
					r = i%N; if (r == 0) r = N;
					c = j%N; if (c == 0) c = N;
					z = blk_cm_getelem(ham,r,c);
					zval.re -= z.re;
					zval.im -= z.im;
				}
				if (Cnorm(zval) >= SPARSE_TOL) {
					(L->irow[ii])++;
					L->icol[k] = jj;
					L->data[k] = zval;
					if (++k > Nmax) {
						fprintf(stderr,"Error: create_liou - Nmax=%d overflow after (%d,%d)\n",Nmax,i,j);
						exit(1);
					}
				}
			}
		}
		DEBUGPRINT("create_liou: nnz=%d (=%d), initially allocated %d\n",L->irow[NN]-1,k,Nmax);
		if (Nmax > k) {
			Nmax = k;
			complx *new_data = (complx*)realloc(L->data, Nmax*sizeof(complx));
			MKL_INT *new_icol = (MKL_INT*)realloc(L->icol, Nmax*sizeof(MKL_INT));
			assert( (new_data != NULL) && (new_icol != NULL) );
			L->data = new_data;
			L->icol = new_icol;
		}
	}

	return L;
}

void sigma_update(mat_complx *sigma, int Nsub, int **sub_idx, complx** rhos)
{
	int i, j, row, col, dim = sigma->row;

	switch (sigma->type) {
	case MAT_DENSE_DIAG:
		free(sigma->data);
		sigma->data = (complx*)malloc(dim*dim*sizeof(complx));
		sigma->type = MAT_DENSE;
	case MAT_DENSE:
		cm_zero(sigma);
		for (i=1; i<=Nsub; i++) {
			for (j=1; j<=LEN(sub_idx[i]); j++) {
				row = sub_idx[i][j] % dim; if (row == 0) row = dim;
				col = (sub_idx[i][j] - 1)/dim;
				sigma->data[row-1+col*dim] = rhos[i][j];
			}
		}
		break;
	case MAT_SPARSE_DIAG:
		sigma->type = MAT_SPARSE;
	case MAT_SPARSE: {
		int len = 0;
		int *dumrow = (int*) malloc((dim+1)*sizeof(int));
		memset(dumrow,0,(dim+1)*sizeof(int));
		for (i=1; i<=Nsub; i++) {
			len += LEN(sub_idx[i]);
			for (j=1; j<=LEN(sub_idx[i]); j++) {
				row = sub_idx[i][j] % dim; if (row == 0) row = dim;
				(dumrow[row])++;
			}
		}
		cm_change_nnz(sigma,len);
		sigma->irow[0]=1;
		for (i=1; i<=dim; i++) sigma->irow[i] = sigma->irow[i-1]+dumrow[i];
		int pos;
		memset(dumrow,0,(dim+1)*sizeof(int));
		for (i=1; i<=Nsub; i++) {
			for (j=1; j<=LEN(sub_idx[i]); j++) {
				row = sub_idx[i][j] % dim; if (row == 0) row = dim;
				col = (sub_idx[i][j] - 1)/dim + 1;
				pos = sigma->irow[row-1] + dumrow[row] - 1;
				sigma->icol[pos] = col;
				sigma->data[pos] = rhos[i][j];
				dumrow[row]++;
			}
		}
		free(dumrow);
		break;}
	default:
		fprintf(stderr,"Error: wsp_sigma_update - invalid sigma matrix type\n");
		exit(1);
	}
}

void sigma_update_partial(mat_complx *sigma, int *idx, complx *rho)
{
	int i, j, row, col, dim;
	dim = sigma->row;

	switch (sigma->type) {
	case MAT_DENSE_DIAG:
		free(sigma->data);
		sigma->data = (complx*)malloc(dim*dim*sizeof(complx));
		sigma->type = MAT_DENSE;
	case MAT_DENSE:
		for (j=1; j<=LEN(idx); j++) {
			row = idx[j] % dim; if (row == 0) row = dim;
			col = (idx[j] - 1)/dim;
			sigma->data[row-1+col*dim] = rho[j];
		}
		break;
	case MAT_SPARSE_DIAG:
	case MAT_SPARSE: {
		int len = LEN(idx);
		int *dumrow = (int*)malloc((dim+1)*sizeof(int));
		memset(dumrow,0,(dim+1)*sizeof(int));
		mat_complx *dum = complx_matrix(dim,dim,MAT_SPARSE,len,sigma->basis);
		for (j=1; j<=len; j++) {
			row = idx[j] % dim; if (row == 0) row = dim;
			(dumrow[row])++;
		}
		dum->irow[0]=1;
		for (i=1; i<=dim; i++) dum->irow[i] = dum->irow[i-1]+dumrow[i];
		int pos;
		memset(dumrow,0,(dim+1)*sizeof(int));
		for (j=1; j<=len; j++) {
			row = idx[j] % dim; if (row == 0) row = dim;
			col = (idx[j] - 1)/dim + 1;
			pos = dum->irow[row-1] + dumrow[row] - 1;
			dum->icol[pos] = col;
			dum->data[pos] = rho[j];
			dumrow[row]++;
		}
		cm_addto(sigma,dum);
		free_complx_matrix(dum);
		free(dumrow);
		break;	}
	default:
		fprintf(stderr,"Error: sigma_update_partial - invalid sigma matrix type\n");
		exit(1);
	}
}

void sigma_update_full(mat_complx *sigma, complx *rho)
{
	int i, j, row, col, dim;
	dim = sigma->row;

	switch (sigma->type) {
	case MAT_DENSE_DIAG:
		free(sigma->data);
		sigma->data = (complx*)malloc(dim*dim*sizeof(complx));
		sigma->type = MAT_DENSE;
	case MAT_DENSE:
		memcpy(sigma->data,&(rho[1]),LEN(rho)*sizeof(complx));
		break;
	case MAT_SPARSE_DIAG:
		sigma->type = MAT_SPARSE;
	case MAT_SPARSE: {
		int len = LEN(rho);
		int *dumrow = (int*)malloc((dim+1)*sizeof(int));
		memset(dumrow,0,(dim+1)*sizeof(int));
		int nnz = 0;
		for (j=1; j<=len; j++) {
			if (Cnorm(rho[j])>SPARSE_TOL) {
				row = j % dim; if (row == 0) row = dim;
				(dumrow[row])++;
				nnz++;
			}
		}
		cm_change_nnz(sigma,nnz);
		sigma->irow[0]=1;
		for (i=1; i<=dim; i++) sigma->irow[i] = sigma->irow[i-1]+dumrow[i];
		int pos;
		memset(dumrow,0,(dim+1)*sizeof(int));
		for (j=1; j<=len; j++) {
			if (Cnorm(rho[j])>SPARSE_TOL) {
				row = j % dim; if (row == 0) row = dim;
				col = (j - 1)/dim + 1;
				pos = sigma->irow[row-1] + dumrow[row] - 1;
				sigma->icol[pos] = col;
				sigma->data[pos] = rho[j];
				dumrow[row]++;
			}
		}
		free(dumrow);
		break;	}
	default:
		fprintf(stderr,"Error: sigma_update_full - invalid sigma matrix type\n");
		exit(1);
	}
}


void spinach_point_acq_partial(Sim_info *sim, Sim_wsp *wsp, complx *rho, complx *detect)
{
	complx z, *ptr;

	ptr = &(wsp->fid[++(wsp->curr_nsig)]);
	z = cv_dotc(detect,rho);
	if (wsp->acqphase != 0.0) {
		double cosph,sinph,phase;
		phase = wsp->acqphase*DEG2RAD;
		cosph=cos(phase);
		sinph=sin(phase);
		ptr->re += cosph*z.re+sinph*z.im;
		ptr->im += -sinph*z.re+cosph*z.im;
	} else {
		ptr->re += z.re;
		ptr->im += z.im;
	}
}

/* exit program if command not supported for spinach */
void spinach_disable_command(Sim_info *sim, char *cmd)
{
	if (sim->imethod == M_SPINACH) {
		fprintf(stderr,"Error: command '%s' disabled for spinach method\n",cmd);
		exit(1);
	}
}

mat_double * pick_liou_ints_mask(Sim_info *sim,Sim_wsp *wsp)
{
	mat_double *Lmask, *m1;
	double t_tmp;

	/* pick Liouvillian-mask at 3 "random" rotor positions */
	ham_hamilton(sim,wsp);
	//dm_print(wsp->ham,"Hamiltonian");
	Lmask = liou_from_ham(wsp->ham_blk);
	//dm_print(L,"Liouvillian");
	dm_nnz2one(Lmask);
	if (sim->taur > 0) {
		t_tmp = wsp->t;
		wsp->t = t_tmp + sim->taur*0.369;
		ham_hamilton(sim,wsp);
		m1 = liou_from_ham(wsp->ham_blk);
		dm_nnz2one(m1);
		dm_addto(Lmask,m1);
		free_double_matrix(m1);
		wsp->t = t_tmp + sim->taur*0.747;
		ham_hamilton(sim,wsp);
		m1 = liou_from_ham(wsp->ham_blk);
		dm_nnz2one(m1);
		dm_addto(Lmask,m1);
		free_double_matrix(m1);
		wsp->t = t_tmp;
	}
	return Lmask;
}

/* Warning: assumes _setrfprop has been called!!! */
mat_double * pick_liou_pulse_mask(Sim_info *sim,Sim_wsp *wsp)
{
	mat_double *Lmask, *m;

	Lmask = pick_liou_ints_mask(sim,wsp);
	m = liou_from_ham(wsp->sumHrf);
	dm_addto(Lmask,m);
	free_double_matrix(m);
	return Lmask;
}

/* works on parts of sparse matrix representation: n=A->row/col, rp=A->irow, ci=A->icol */
void path_tracing(int n, int *rp, int *ci, int *Nsub_ptr, int ***sub_idx_ptr)
{
	int *member_nodes, i, j, Nsub, *sub_dim, **sub_idx;

	member_nodes = scomponents(n, rp, ci);
		/* test printout */
		//printf("Path tracing - member_nodes:\n");
		//for (i=1; i<=LEN(member_nodes); i++) printf("\t%d\n",member_nodes[i]);
	Nsub = iv_max(member_nodes);
	sub_dim = int_vector(Nsub);
	iv_zero(sub_dim);
	for (i=1; i<=LEN(member_nodes); i++) {
		(sub_dim[member_nodes[i]])++;
	}
	sub_idx = (int**)malloc((Nsub+1)*sizeof(int*));
	for (i=1; i<=Nsub; i++) {
		sub_idx[i] = int_vector(sub_dim[i]);
		sub_dim[i] = 0;
		for (j=1; j<=LEN(member_nodes); j++) {
			if (member_nodes[j] == i) {
				(sub_dim[i])++;
				sub_idx[i][sub_dim[i]] = j;
			}
		}
	}
	free_int_vector(sub_dim);
	free_int_vector(member_nodes);
	*Nsub_ptr = Nsub;
	*sub_idx_ptr = sub_idx;
}

void spinach_pulseid(Sim_info *sim,Sim_wsp *wsp,double duration)
{
	_pulseid(sim,wsp,duration);
	_evolve_with_prop(sim,wsp);
	_reset_prop(sim,wsp);
}

void spinach_pulse_simple(Sim_info *sim, Sim_wsp *wsp,double duration,complx *rho,int *idx, blk_mat_complx *Hrf)
{
	int i, j, n;
	double dt, dt_us;
	mat_complx *Ltot;
	blk_mat_complx *Htot = NULL;

	n = (int)ceil(duration/wsp->dtmax);
	if (n < 1) n = 1;
	dt_us = duration/(double)n;
	dt = dt_us*1.0e-6;
	DEBUGPRINT("spinach_pulse_simple duration of %f us split into %d steps of %f us\n",duration,n,dt_us);
	for (i=1;i<=n;i++) {
		ham_hamilton(sim,wsp);
		if (Htot == NULL) Htot = create_blk_mat_complx_copy2(wsp->ham_blk);
		assert(wsp->ham_blk->Nblocks == Hrf->Nblocks);
		for (j=0; j<Htot->Nblocks; j++) {
			assert(wsp->ham_blk->blk_dims[j] == Hrf->blk_dims[j]);
			if (Htot->blk_dims[j] == 1) {
				Htot->m[j].data[0] = Complx(wsp->ham_blk->m[j].data[0]+Hrf->m[j].data[0].re,Hrf->m[j].data[0].im);
			} else {
				dm_copy2cm(wsp->ham_blk->m + j, Htot->m + j);
				cm_addto(Htot->m + j, Hrf->m + j);
			}
		}
		//dm_copy2cm(wsp->ham,Htot);
		//cm_addto(Htot,Hrf);
		Ltot = create_liou(Htot,idx);
		prop_krylov(-dt, Ltot, rho);
		wsp->t += dt_us;
		free_complx_matrix(Ltot);
	}
	free_blk_mat_complx(Htot);
	DEBUGPRINT("spinach_pulse_simple: ---> DONE <---\n");
}

/* there were serious trials to use zero track elimination and path tracing
 * also with rf pulse evolutions but without general success. Space reduction
 * and subspace division do not lead neither to dimensions smaller to Hilbert
 * space nor dimensions where it pays back by matrix.vector operations instead
 * of matrix.matrix.matrix. Final decision is to make pulse propagation via
 * Krylov subspaces and avoiding matrix exponentialization. No propagators are
 * created...
 */
void spinach_pulse(Sim_info *sim,Sim_wsp *wsp,double duration)
{
	int Nchan_active, Nsub, **sub_idx, temp_curr_nsig, temp_Nacq, n, i, j;
	mat_double *Lmask;
	mat_complx *temp_sigma;
	blk_mat_complx *Hrf;
	double temp_curr_time, temp_curr_acqblock_t0, dt, t_tmp;
	complx *rho, *detect;
	const int verb = verbose & VERBOSE_PROGRESS;

	if (wsp->evalmode == EM_ACQBLOCK) {
		if (wsp->Nacq == 0) return;
	}

	/* set up sumHrf matrix and see if there is any RF */
	Nchan_active = _setrfprop(sim,wsp);
	if (Nchan_active == 0) {
		spinach_delay(sim,wsp,duration);
		return;
	}

	/* decide on using subspaces */
	if (wsp->spinach_pulse_Nsub == 0 ) {
		if (Nchan_active < sim->ss->nchan) {
			/* there is good chance to get smaller matrices using path tracing */
			Lmask = pick_liou_pulse_mask(sim,wsp);
			path_tracing(Lmask->row,Lmask->irow,Lmask->icol,&Nsub,&sub_idx);
			free_double_matrix(Lmask);
			if (verb) printf("spinach_pulse: fresh new path tracing of pulse evolution <---\n");
		} else {
			Nsub = 0;
			sub_idx = NULL;
		}
	} else {
		Nsub = wsp->spinach_pulse_Nsub;
		sub_idx = wsp->spinach_pulse_sub_idx;
	}

		/* test print out *
		for (i=1; i<=Nsub; i++) {
			printf("\nSubspace %d/%d, dim=%d\n\telements:",i,Nsub,LEN(sub_idx[i]));
			for (j=1; j<=LEN(sub_idx[i]); j++) {
				printf(" %d",sub_idx[i][j]);
			}
		}
		printf("\n");
		**************/

	/* generate complex RF hamiltonian (includes pulse phase) */
	Hrf = ham_rf(wsp);

	if (Nsub > 0) {
		temp_curr_nsig = wsp->curr_nsig;
		temp_Nacq = wsp->Nacq;
		temp_curr_time = wsp->t;
		temp_curr_acqblock_t0 = wsp->acqblock_t0;
		temp_sigma = cm_dup(wsp->sigma);
		cm_zero(wsp->sigma);
		for (i=1; i<=Nsub; i++) {
			rho = get_rho(temp_sigma,sub_idx[i],0);
			if (cv_asum(rho) < TINY) {
				if (verb) printf("\tSubspace %d  (dim=%d) will be skipped\n",i,LEN(rho));
				free_complx_vector(rho);
				continue;
			}
			if (verb) printf("spinach_pulse working in subspace %d (dim=%d)\n--------------\n",i,LEN(rho));
			if (wsp->evalmode == EM_ACQBLOCK) {
				detect = get_rho(wsp->fdetect,sub_idx[i],!(sim->acq_adjoint));
				wsp->t = temp_curr_time;
				wsp->curr_nsig = temp_curr_nsig;
				wsp->Nacq = temp_Nacq;
				wsp->acqblock_t0 = temp_curr_acqblock_t0;
				if ( fabs(wsp->t - wsp->acqblock_t0) > TINY ) {
					/* there was some previous event not synchronized with dw */
					dt = wsp->dw - (wsp->t - wsp->acqblock_t0);
					if ( duration < dt) {
						spinach_pulse_simple(sim,wsp,duration,rho,sub_idx[i],Hrf);
						t_tmp = 0;
					} else {
						spinach_pulse_simple(sim,wsp,dt,rho,sub_idx[i],Hrf);
						t_tmp = duration - dt;
					}
					/* did it fill time up to dw? */
					if (fabs(wsp->t - wsp->acqblock_t0 - wsp->dw) < TINY ) {
						spinach_point_acq_partial(sim,wsp,rho,detect);
						wsp->acqblock_t0 = wsp->t;
						(wsp->Nacq)--;
						DEBUGPRINT("\tpulse\t wsp->Nacq = %d\n",wsp->Nacq);
						if (wsp->Nacq == 0) {
							free_complx_vector(detect);
							free_complx_vector(rho);
							continue;
						}
					}
				} else {
					t_tmp = duration;
				}
				if ( fabs(t_tmp) > TINY) {
					n = (int)floor(t_tmp/wsp->dw+1e-6);
					dt = t_tmp - wsp->dw*(double)n;
					for (j=1; j<=n; j++) {
						spinach_pulse_simple(sim,wsp,wsp->dw,rho,sub_idx[i],Hrf);
						spinach_point_acq_partial(sim,wsp,rho,detect);
						wsp->acqblock_t0 = wsp->t;
						(wsp->Nacq)--;
						DEBUGPRINT("\tpulse\t wsp->Nacq = %d\n",wsp->Nacq);
						if (wsp->Nacq == 0) break;
					}
					if (dt > TINY && wsp->Nacq > 0) {
						/* remaining time not synchronized with dw */
						spinach_pulse_simple(sim,wsp,dt,rho,sub_idx[i],Hrf);
					}
				}
				free_complx_vector(detect);
			} else {
				/* just evolve */
				wsp->t = temp_curr_time;
				spinach_pulse_simple(sim,wsp,duration,rho,sub_idx[i],Hrf);
			}
			/* update resulting density matrix wsp->sigma */
			sigma_update_partial(wsp->sigma,sub_idx[i],rho);
			free_complx_vector(rho);
		}
		free_complx_matrix(temp_sigma);
	} else {
		/* need to do full space evolution */
		rho = get_sigma(wsp->sigma,0);
		if (wsp->evalmode == EM_ACQBLOCK) {
			detect = get_sigma(wsp->fdetect,!(sim->acq_adjoint));
			if ( fabs(wsp->t - wsp->acqblock_t0) > TINY ) {
				/* there was some previous event not synchronized with dw */
				dt = wsp->dw - (wsp->t - wsp->acqblock_t0);
				if ( duration < dt) {
					spinach_pulse_simple(sim,wsp,duration,rho,NULL,Hrf);
					t_tmp = 0;
				} else {
					spinach_pulse_simple(sim,wsp,dt,rho,NULL,Hrf);
					t_tmp = duration - dt;
				}
				/* did it fill time up to dw? */
				if (fabs(wsp->t - wsp->acqblock_t0 - wsp->dw) < TINY ) {
					spinach_point_acq_partial(sim,wsp,rho,detect);
					wsp->acqblock_t0 = wsp->t;
					(wsp->Nacq)--;
					DEBUGPRINT("\t pulse\t wsp->Nacq = %d\n",wsp->Nacq);
					if (wsp->Nacq == 0) {
						free_complx_vector(detect);
						free_complx_vector(rho);
						free_blk_mat_complx(Hrf);
						return;
					}
				}
			} else {
				t_tmp = duration;
			}
			if ( fabs(t_tmp) > TINY) {
				n = (int)floor(t_tmp/wsp->dw+1e-6);
				dt = t_tmp - wsp->dw*(double)n;
				for (i=1; i<=n; i++) {
					spinach_pulse_simple(sim,wsp,wsp->dw,rho,NULL,Hrf);
					spinach_point_acq_partial(sim,wsp,rho,detect);
					wsp->acqblock_t0 = wsp->t;
					(wsp->Nacq)--;
					DEBUGPRINT("\t pulse\t wsp->Nacq = %d\n",wsp->Nacq);
					if (wsp->Nacq == 0) {
						free_complx_vector(detect);
						free_complx_vector(rho);
						free_blk_mat_complx(Hrf);
						return;
					}
				}
				if (dt > TINY) {
					/* remaining time not synchronized with dw */
					spinach_pulse_simple(sim,wsp,dt,rho,NULL,Hrf);
				}
			}
			free_complx_vector(detect);
		} else {
			spinach_pulse_simple(sim,wsp,duration,rho,NULL,Hrf);
		}
		/* update wsp->sigma */
		sigma_update_full(wsp->sigma,rho);
		free_complx_vector(rho);
	}

	free_blk_mat_complx(Hrf);
	/* do not remember subspaces if they were not known previously */
	if (wsp->spinach_pulse_Nsub != Nsub ) {
		for (i=1; i<=Nsub; i++) free_int_vector(sub_idx[i]);
		free(sub_idx);
	}
	if (verb) printf("spinach_pulse: ---> DONE <---\n");
}

void spinach_delay_simple(Sim_info *sim,Sim_wsp *wsp,double duration, complx *rho, mat_double *Liou, mat_complx *Prop,int *idx)
{
	int Ndt, i;
	double dt, dt_us;

	Ndt = (int)ceil(duration/wsp->dtmax);
	if (Ndt < 1) Ndt = 1;
	dt_us = duration/(double)Ndt;
	dt = dt_us*1.0e-6;
	DEBUGPRINT("\t spinach_delay_simple %f us split into %d steps of %f us length\n",duration,Ndt,dt_us);
	for (i=0; i<Ndt; i++) {
		/* prepare current hamiltonian and get Liou */
		ham_hamilton(sim,wsp);
		ham_liou_restr(wsp->ham_blk,idx,Liou);
		prop_real(Prop,Liou,-dt,sim->propmethod);
		/* apply prop to rho */
		cv_matmulto(rho,Prop);
		wsp->t += dt_us;
	}
}

void spinach_delay(Sim_info *sim,Sim_wsp *wsp,double duration)
{
	mat_double *L, *Liou;
	int i, j, n, Nsub, **sub_idx, temp_curr_nsig, temp_Nacq, subdim, analyze_fdetect;
	complx *rho, *detect;
	double t_tmp, dt, temp_curr_time,temp_curr_acqblock_t0;
	mat_complx *Prop, *temp_sigma;
	const int verb = verbose & VERBOSE_PROGRESS;

	if (wsp->evalmode == EM_ACQBLOCK) {
		if (wsp->Nacq == 0) return;
	}

	if (sim->Hint_isdiag) {
		if (verb) printf("spinach_delay: diagonal Hamiltonian, evolution in Hilbert space\n");
		_delay(sim,wsp,duration);
		_evolve_with_prop(sim,wsp);
		_reset_prop(sim,wsp);
		return;
	}

	if (wsp->spinach_delay_Nsub == 0) {
		L = pick_liou_ints_mask(sim,wsp);
		//dm_print(L,"final Liou-mask");
		/* do path tracing */
		path_tracing(L->row,L->irow,L->icol,&(wsp->spinach_delay_Nsub),&(wsp->spinach_delay_sub_idx));
		free_double_matrix(L);
		if (verb) printf("spinach_delay: fresh new path tracing of delay evolution <---\n");
	}
	Nsub = wsp->spinach_delay_Nsub;
	sub_idx = wsp->spinach_delay_sub_idx;

		/* test print out */
	/*	for (i=1; i<=Nsub; i++) {
			printf("\nSubspace %d, dim=%d\n\telements:",i,LEN(sub_idx[i]));
			for (j=1; j<=LEN(sub_idx[i]); j++) {
				printf(" %d",sub_idx[i][j]);
			}
		}
		printf("\n");
	*/
		/* test print out in nice form for analysis */
	/*	printf("======= Subspace report - all dims =======\n");
		for (i=1; i<=Nsub; i++) {
			printf("%d\n",LEN(sub_idx[i]));
		}
	*/

	temp_curr_nsig = wsp->curr_nsig;
	temp_Nacq = wsp->Nacq;
	temp_curr_time = wsp->t;
	temp_curr_acqblock_t0 = wsp->acqblock_t0;
	temp_sigma = cm_dup(wsp->sigma);
	cm_zero(wsp->sigma);

	/* do calculation only in relevant subspaces (= those with non-zero state vector) */
	analyze_fdetect = (wsp->dw*wsp->Nacq <= duration) ? 1 : 0;
	//printf("\n++++  %f vs %f ++++\n",sim->dw*wsp->Nacq,duration);
	if (verb) printf("spinach_delay: fdetect will%s be analyzed",analyze_fdetect ? "" : " not");
	for (i=1; i<=Nsub; i++) {
		rho = get_rho(temp_sigma,sub_idx[i],0);
		subdim = LEN(rho);
		if (cv_asum(rho) < TINY) {
			if (verb) printf("\tSubspace %d  (dim=%d) will be skipped\n",i,subdim);
			free_complx_vector(rho);
			continue;
		}
		if (wsp->evalmode == EM_ACQBLOCK) {
			detect = get_rho(wsp->fdetect,sub_idx[i],!(sim->acq_adjoint));
			if (analyze_fdetect) {
				if (cv_asum(detect) < TINY) {
					if (verb) printf("\tSubspace %d (dim=%d) will be skipped - no detection\n",i,LEN(detect));
					free_complx_vector(rho);
					free_complx_vector(detect);
					continue;
				}
			}
			if (verb) printf("Subspace %d (dim=%d) will be retained for acquisition\n",i,subdim);
			wsp->t = temp_curr_time;
			wsp->curr_nsig = temp_curr_nsig;
			wsp->Nacq = temp_Nacq;
			wsp->acqblock_t0 = temp_curr_acqblock_t0;
			Prop = complx_matrix(subdim,subdim,MAT_DENSE,0,0);
			Liou = double_matrix(subdim,subdim,MAT_DENSE,0,0);
			if ( fabs(wsp->t - wsp->acqblock_t0) > TINY ) {
				/* there was some previous event not synchronized with dw */
				dt = wsp->dw - (wsp->t - wsp->acqblock_t0);
				if ( duration < dt) {
					spinach_delay_simple(sim,wsp,duration,rho,Liou,Prop,sub_idx[i]);
					t_tmp = 0;
				} else {
					spinach_delay_simple(sim,wsp,dt,rho,Liou,Prop,sub_idx[i]);
					t_tmp = duration - dt;
				}
				/* did it fill time up to dw? */
				if (fabs(wsp->t - wsp->acqblock_t0 - wsp->dw) < TINY ) {
					spinach_point_acq_partial(sim,wsp,rho,detect);
					(wsp->Nacq)--;
					DEBUGPRINT("\t\t wsp->Nacq = %d\n",wsp->Nacq);
					if (wsp->Nacq == 0) {
						free_double_matrix(Liou);
						free_complx_matrix(Prop);
						free_complx_vector(detect);
						free_complx_vector(rho);
						continue;}
					wsp->acqblock_t0 = wsp->t;
				}
			} else {
				t_tmp = duration;
			}
			if ( fabs(t_tmp) > TINY) {
				n = (int)floor(t_tmp/wsp->dw+1e-6);
				dt = t_tmp - wsp->dw*(double)n;
				if (sim->wr <TINY) {
					/* static case, can re-use dwelltime propagator */
					if (n > 0) {
						spinach_delay_simple(sim,wsp,wsp->dw,rho,Liou,Prop,sub_idx[i]);
						spinach_point_acq_partial(sim,wsp,rho,detect);
						wsp->acqblock_t0 = wsp->t;
						(wsp->Nacq)--;
						DEBUGPRINT("\t\t wsp->Nacq = %d\n",wsp->Nacq);
						if (wsp->Nacq == 0) {
							free_double_matrix(Liou);
							free_complx_matrix(Prop);
							free_complx_vector(detect);
							free_complx_vector(rho);
							continue;
						}
						for (j=2; j<=n; j++) {
							/* apply prop to rho */
							cv_matmulto(rho,Prop);
							wsp->t += wsp->dw;;
							spinach_point_acq_partial(sim,wsp,rho,detect);
							wsp->acqblock_t0 = wsp->t;
							(wsp->Nacq)--;
							DEBUGPRINT("\t\t wsp->Nacq = %d\n",wsp->Nacq);
							if (wsp->Nacq == 0) break;
						}
					}
				} else {
					for (j=1; j<=n; j++) {
						spinach_delay_simple(sim,wsp,wsp->dw,rho,Liou,Prop,sub_idx[i]);
						spinach_point_acq_partial(sim,wsp,rho,detect);
						wsp->acqblock_t0 = wsp->t;
						(wsp->Nacq)--;
						DEBUGPRINT("\t\t wsp->Nacq = %d\n",wsp->Nacq);
						if (wsp->Nacq == 0) break;
					}
				}
				if (dt > TINY && wsp->Nacq > 0) {
					/* remaining time not synchronized with dw */
					spinach_delay_simple(sim,wsp,dt,rho,Liou,Prop,sub_idx[i]);
				}
			}
			free_double_matrix(Liou);
			free_complx_matrix(Prop);
			free_complx_vector(detect);
		} else {
			if (verb) printf("Subspace %d (dim=%d) will be retained\n",i,subdim);
			/* just evolve */
			wsp->t = temp_curr_time;
			Prop = complx_matrix(subdim,subdim,MAT_DENSE,0,0);
			Liou = double_matrix(subdim,subdim,MAT_DENSE,0,0);
			spinach_delay_simple(sim,wsp,duration,rho,Liou,Prop,sub_idx[i]);
			free_double_matrix(Liou);
			free_complx_matrix(Prop);
		}
		/* update resulting density matrix wsp->sigma */
		sigma_update_partial(wsp->sigma,sub_idx[i],rho);
		free_complx_vector(rho);
	} //end of loop over subspaces

	free_complx_matrix(temp_sigma);
	wsp->Uisunit = 0;
}



void spinach_acqblock(Tcl_Interp *interp, Tcl_Obj *obj, Sim_info *sim, Sim_wsp *wsp)
{
	complx z, *ptr;
	double t_tmp, t_acq;
	int i, Nrep;

	t_tmp = wsp->acqblock_t0 = wsp->t;
	t_acq = wsp->dw*(wsp->Nacq-1);

	/* acquire first point */
	ptr = &(wsp->fid[++(wsp->curr_nsig)]);
	if (sim->acq_adjoint == 0) {
		z = cm_trace(wsp->fdetect,wsp->sigma);
	} else {
		z = cm_trace_adjoint(wsp->fdetect,wsp->sigma);
	}
	if (fabs(wsp->acqphase) > TINY) {
		double cosph,sinph;
		cosph=cos(wsp->acqphase*DEG2RAD);
		sinph=sin(wsp->acqphase*DEG2RAD);
		ptr->re += cosph*z.re+sinph*z.im;
		ptr->im += -sinph*z.re+cosph*z.im;
	} else {
		ptr->re += z.re;
		ptr->im += z.im;
	}
	(wsp->Nacq)--;
	if (Tcl_EvalObjEx(interp,obj,0) != TCL_OK) {
		fprintf(stderr,"acq_block error: (1) can not execute block:\n'%s'\n",Tcl_GetString(obj));
		exit(1);
	}
	if (fabs(t_tmp-wsp->t) < TINY) {
		fprintf(stderr,"Error: acq_block of zero duration\n");
		exit(1);
	}
	Nrep = (int)ceil(t_acq/(wsp->t-t_tmp));
	DEBUGPRINT("spinach_acq_block: elapsed time %g us, acq time %g us, Nrep %d\n",wsp->t-t_tmp,t_acq,Nrep);
	for (i=1; i<Nrep; i++) {
		if (Tcl_EvalObjEx(interp,obj,0) != TCL_OK) {
			fprintf(stderr,"acq_block error: (%d) can not execute block:\n'%s'\n",i+1,Tcl_GetString(obj));
			exit(1);
		}
	}
}



void tclcmd_spinach(Tcl_Interp* interp) {


}
