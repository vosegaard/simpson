/*
 * blockdiag.c
 *
 *  Created on: 12.3.2011
 *      Author: zdenek
 */

/***********************************************************
 *   B L O C K   D I A G O N A L   M A T R I C E S
 *********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <inttypes.h>
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
#include "cm.h"
#include "blockdiag.h"
#include "defs.h"

	/* for acurate timings on windows */
//#include <windows.h>
//#include <winbase.h>

void blk_dm_zero(blk_mat_double *obj)
{
	int i;

	for (i=0; i<obj->Nblocks; i++) {
		dm_zero(&(obj->m[i]));
	}
}

void blk_cm_zero(blk_mat_complx *obj)
{
	int i;

	for (i=0; i<obj->Nblocks; i++) {
		cm_zero(&(obj->m[i]));
	}
}

void blk_cm_unit(blk_mat_complx *obj)
{
	int i;

	for (i=0; i<obj->Nblocks; i++) {
		cm_unit(&(obj->m[i]));
	}
}

void blk_dm_multod_diag(blk_mat_double *mblk, mat_double *mdiag, double d)
{
	assert(mdiag->type == MAT_DENSE_DIAG);
	assert(mdiag->row == mblk->dim);
	int i;
	double *orig, *curr;

	orig = curr = mdiag->data;
	for (i=0; i<mblk->Nblocks; i++) {
		if (mblk->blk_dims[i] == 1) {
			mblk->m[i].data[0] += (*curr)*d;
			curr++;
		} else {
			mdiag->data = curr;
			mdiag->row = mdiag->col = mblk->blk_dims[i];
			dm_multod(&(mblk->m[i]),mdiag,d);
			curr += mblk->blk_dims[i];
		}
	}
	mdiag->data = orig;
	mdiag->row = mdiag->col = mblk->dim;
}

void blk_dm_multod(blk_mat_double *res, blk_mat_double *m, double d)
{
	assert(res->dim == m->dim);
	assert(res->Nblocks == m->Nblocks);
	int i;

	for (i=0; i<res->Nblocks; i++) {
		assert(res->blk_dims[i] == m->blk_dims[i]);
		if (m->blk_dims[i] == 1) {
			res->m[i].data[0] += d*(m->m[i].data[0]);
		} else {
			dm_multod(&(res->m[i]),&(m->m[i]),d);
		}
	}
}


/****
 * extract from double matrix a diagonal block of size dim defined by
 * permutation vector permvec (starting from position sft)
 ****/
mat_double * dm_get_diagblock_permute(mat_double *a, int *permvec, int sft, int dim)
{
	mat_double *res;
	int i,j;

	switch (a->type) {
		case MAT_DENSE:
			res = double_matrix(dim,dim,MAT_DENSE,0,a->basis);
			for (i=0; i<dim; i++) {
				for (j=0; j<dim; j++) {
					res->data[i+j*dim] = a->data[permvec[sft+i]-1+(permvec[sft+j]-1)*a->row];
				}
			}
			break;
		case MAT_SPARSE: {
			uint64_t NNzmax = (uint64_t)floor((uint64_t)dim*dim*(1.0-SPARSITY));
			if (dim < MAXFULLDIM) {
				res = double_matrix(dim,dim,MAT_DENSE,0,a->basis);
				for (i=0; i<dim; i++) {
					for (j=0; j<dim; j++) {
						res->data[i+j*dim] = dm_getelem(a,permvec[sft+i],permvec[sft+j]);
					}
				}
			} else {
				//printf("in DGDP\n");
				int r, nnz = 0;
				double val;
				// count non-zeros
				for (i=0; i<dim; i++) {
					r = permvec[sft+i];
					if (a->irow[r] - a->irow[r-1] == 0) continue;
					for (j=0; j<dim; j++) {
						val = dm_getelem(a,r,permvec[sft+j]);
						if (fabs(val) < SPARSE_TOL) continue;
						nnz++;
					}
				}
				// allocate and fill
				res = double_matrix(dim,dim,MAT_SPARSE,nnz,a->basis);
				res->irow[0] = 1;
				nnz = 0;
				for (i=0; i<dim; i++) {
					res->irow[i+1] = res->irow[i];
					r = permvec[sft+i];
					if (a->irow[r] - a->irow[r-1] == 0) continue;
					for (j=0; j<dim; j++) {
						val = dm_getelem(a,r,permvec[sft+j]);
						if (fabs(val) < SPARSE_TOL) continue;
						res->icol[nnz] = j+1;
						res->data[nnz] = val;
						res->irow[i+1]++;
						nnz++;
					}
				}

				if (nnz > NNzmax) {
					dm_dense(res);
				} else {
					dm_change_nnz(res,nnz);
				}
			}
			break; }
		case MAT_DENSE_DIAG:
			res = double_matrix(dim,dim,MAT_DENSE_DIAG,0,a->basis);
			for (i=0; i<dim; i++) {
				res->data[i] = a->data[permvec[sft+i]-1];
			}
			break;
		case MAT_SPARSE_DIAG:
			res = double_matrix(dim,dim,MAT_DENSE_DIAG,dim,a->basis);
			for (i=0; i<dim; i++) {
				j = permvec[sft+i];
				res->data[i] = dm_getelem(a,j,j);
			}
			break;
		default:
			fprintf(stderr,"Error: dm_get_diagblock_permute - unknown matrix type %d\n",a->type);
			exit(1);
	}
	return res;
}

blk_mat_double * blk_dm_dup(blk_mat_double *m)
{
	// nezapomenout kopirovat i info o bazi
	blk_mat_double *res;
	mat_double *rmx, *mx;
	int i, n;

	res = (blk_mat_double*)malloc(sizeof(blk_mat_double));
	res->dim = m->dim;
	res->Nblocks = m->Nblocks;
	res->basis = m->basis;
	res->m = (mat_double*)malloc(res->Nblocks*sizeof(mat_double));
	if (res->m == NULL) {
		fprintf(stderr,"Error: blk_dm_dup can not allocate block-matrix pointers\n");
		exit(1);
	}
	res->blk_dims = (int*)malloc(res->Nblocks*sizeof(int));
	if (res->blk_dims == NULL) {
		fprintf(stderr,"Error: blk_dm_dup can not allocate blk_dims\n");
		exit(1);
	}
	for (i=0; i<res->Nblocks; i++) {
		res->blk_dims[i] = n = m->blk_dims[i];
		rmx = res->m + i;
		mx = m->m + i;
		if (n == 1) {
			create_sqmat_double(rmx,1,MAT_DENSE,1,m->basis);
			rmx->data[0] = mx->data[0];
		} else {
			switch (mx->type) {
			case MAT_DENSE:
				create_sqmat_double(rmx,n,MAT_DENSE,1,m->basis);
				memcpy(rmx->data,mx->data,n*n*sizeof(double));
				break;
			case MAT_DENSE_DIAG:
				create_sqmat_double(rmx,n,MAT_DENSE_DIAG,1,m->basis);
				memcpy(rmx->data,mx->data,n*sizeof(double));
				break;
			case MAT_SPARSE: {
				int len = mx->irow[n] -1;
				create_sqmat_double(rmx,n,MAT_SPARSE,len,m->basis);
				memcpy(rmx->data, mx->data, len*sizeof(double));
				memcpy(rmx->icol, mx->icol, len*sizeof(int));
				memcpy(rmx->irow, mx->irow, (n+1)*sizeof(int));
				break; }
			case MAT_SPARSE_DIAG: {
				int len = mx->irow[n] -1;
				create_sqmat_double(rmx,n,MAT_SPARSE_DIAG,len,m->basis);
				memcpy(rmx->data, mx->data, len*sizeof(double));
				memcpy(rmx->icol, mx->icol, len*sizeof(int));
				memcpy(rmx->irow, mx->irow, (n+1)*sizeof(int));
				break; }
			default:
				fprintf(stderr,"Error: blk_dm_dup - unknown input matrix %d type %d\n",i,mx->type);
				exit(1);
			}
		}
	}
	return res;
}


void blk_dm_copy(blk_mat_double *res, blk_mat_double *m)
{
	int i;

	assert(res->Nblocks == m->Nblocks);
	assert(res->dim == m->dim);
	for (i=0; i<res->Nblocks; i++) {
		assert(res->blk_dims[i] == m->blk_dims[i]);
		dm_copy(&(res->m[i]),&(m->m[i]));
	}
}

void blk_dm_muld(blk_mat_double *m, double d)
{
	int i;

	for (i=0; i<m->Nblocks; i++) {
		dm_muld(&(m->m[i]),d);
	}
}

void blk_dm_print(blk_mat_double *obj, char *title)
{
	int i;
	char cc[16];

	printf("\n\n%s : block diagonal matrix of dim %d in basis %d\n",title,obj->dim, obj->basis);
	for (i=0; i<obj->Nblocks; i++) {
		sprintf(cc,"BLOCK %d",i);
		dm_print(&(obj->m[i]),cc);
	}
	printf("\n");
}

void blk_cm_print(blk_mat_complx *obj, char *title)
{
	int i;
	char cc[16];

	printf("\n\n%s : block diagonal matrix of dim %d in basis %d\n",title,obj->dim, obj->basis);
	for (i=0; i<obj->Nblocks; i++) {
		sprintf(cc,"BLOCK %d",i);
		cm_print(&(obj->m[i]),cc);
	}
	printf("\n");
}





mat_double * dm_change_basis_2(mat_double *m, int basis, Sim_info *sim)
{
	// creates new matrix from m by changing its basis
	mat_double *res;

	if (m->basis == basis) {
		res = dm_dup(m);
	} else {
		int *permvec = sim->perm_table[m->basis + basis*sim->Nbasis];
		int dim = m->row;
		assert( (LEN(permvec) == dim) && (dim == m->col) );
		switch (m->type) {
		case MAT_DENSE: {
			int i, j;
			res = double_matrix(dim, dim, MAT_DENSE, 0, basis);
			for (i=0; i<dim; i++) {
				for (j=0; j<dim; j++) {
					res->data[i+j*dim] = m->data[permvec[i+1]-1 +(permvec[j+1]-1)*dim];
				}
			}
			break;
		}
		case MAT_DENSE_DIAG: {
			int i;
			res = double_matrix(dim,dim,MAT_DENSE_DIAG,0, basis);
			for(i=0; i<dim; i++) res->data[i] = m->data[permvec[i+1]-1];
			break;
		}
		case MAT_SPARSE: {
			int i, j, jj, c_orig, k, r, nc, nnz = m->irow[dim]-1;
			res = double_matrix(dim,dim,MAT_SPARSE,nnz,basis);
			res->irow[0] = 1;
			int pos = 0;
			for (i=1; i<=dim; i++) {
				r = permvec[i];
				nc = m->irow[r] - m->irow[r-1];
				for (j=0; j<nc; j++) {
					c_orig = m->icol[m->irow[r-1]-1+j];
					for (k=1; k<=dim; k++) {
						if (permvec[k] == c_orig) break;
					}
					if (j == 0) {
						res->icol[pos] = k;
						res->data[pos] = m->data[m->irow[r-1]-1+j];
					} else {
						for (jj=0; jj<j; jj++) {
							if (k < res->icol[pos-jj-1]) {
								res->icol[pos-jj] = res->icol[pos-jj-1];
								res->data[pos-jj] = res->data[pos-jj-1];
							} else {
								break;
							}
						}
						res->icol[pos-jj] = k;
						res->data[pos-jj] = m->data[m->irow[r-1]-1+j];
					}
					pos++;
				}
				res->irow[i] = pos+1;
			}
			break;
		}
		case MAT_SPARSE_DIAG: {
			int i,j, nnz = m->irow[dim]-1;
			res = double_matrix(dim,dim,MAT_SPARSE_DIAG,nnz,basis);
			res->irow[0] = 1;
			j = 0;
			double val;
			for (i=1; i<=dim; i++) {
				res->irow[i] = res->irow[i-1];
				val = dm_getelem(m,permvec[i],permvec[i]);
				if (fabs(val)<SPARSE_TOL ) continue;
				res->icol[j] = i;
				res->data[j] = val;
				(res->irow[i])++;
				j++;
			}
			break;
		}
		default:
			fprintf(stderr,"Error: dm_change_basis_2 - unknown matrix type %d\n",m->type);
			exit(1);
		}
	}

	return res;
}

mat_complx * cm_change_basis_2(mat_complx *m, int basis, Sim_info *sim)
{
	// create new matrix with change of basis
	mat_complx *res;

	if (m->basis == basis) {
		res = cm_dup(m);
	} else {
		int *permvec = sim->perm_table[m->basis + basis*sim->Nbasis];
		int dim = m->row;
		assert( (LEN(permvec) == dim) && (dim == m->col) );
		switch (m->type) {
		case MAT_DENSE: {
			int i, j;
			res = complx_matrix(dim, dim, MAT_DENSE, 0,basis);
			for (i=0; i<dim; i++) {
				for (j=0; j<dim; j++) {
					res->data[i+j*dim] = m->data[permvec[i+1]-1 +(permvec[j+1]-1)*dim];
				}
			}
			break;
		}
		case MAT_DENSE_DIAG: {
			int i;
			res = complx_matrix(dim,dim,MAT_DENSE_DIAG,0,basis);
			for(i=0; i<dim; i++) res->data[i] = m->data[permvec[i+1]-1];
			break;
		}
		case MAT_SPARSE: {
			int i, j, jj, c_orig, k, r, nc, nnz = m->irow[dim]-1;
			res = complx_matrix(dim,dim,MAT_SPARSE,nnz,basis);
			res->irow[0] = 1;
			int pos = 0;
			for (i=1; i<=dim; i++) {
				r = permvec[i];
				nc = m->irow[r] - m->irow[r-1];
				for (j=0; j<nc; j++) {
					c_orig = m->icol[m->irow[r-1]-1+j];
					for (k=1; k<=dim; k++) {
						if (permvec[k] == c_orig) break;
					}
					if (j == 0) {
						res->icol[pos] = k;
						res->data[pos] = m->data[m->irow[r-1]-1+j];
					} else {
						for (jj=0; jj<j; jj++) {
							if (k < res->icol[pos-jj-1]) {
								res->icol[pos-jj] = res->icol[pos-jj-1];
								res->data[pos-jj] = res->data[pos-jj-1];
							} else {
								break;
							}
						}
						res->icol[pos-jj] = k;
						res->data[pos-jj] = m->data[m->irow[r-1]-1+j];
					}
					pos++;
				}
				res->irow[i] = pos+1;
			}
			break;
		}
		case MAT_SPARSE_DIAG: {
			int i,j, nnz = m->irow[dim]-1;
			res = complx_matrix(dim,dim,MAT_SPARSE_DIAG,nnz,basis);
			res->irow[0] = 1;
			j = 0;
			complx zval;
			for (i=1; i<=dim; i++) {
				res->irow[i] = res->irow[i-1];
				zval = cm_getelem(m,permvec[i],permvec[i]);
				if (fabs(zval.re)<SPARSE_TOL && fabs(zval.im)<SPARSE_TOL) continue;
				res->icol[j] = i;
				res->data[j] = zval;
				(res->irow[i])++;
				j++;
			}
			break;
		}
		default:
			fprintf(stderr,"Error: cm_change_basis_2 - unknown matrix type %d\n",m->type);
			exit(1);
		}
	}

	return res;
}


double blk_dm_getelem(blk_mat_double *obj, int row, int col)
{
	// row and col are 1-based
	int rm=0, cm=0;
	double val = 0.0;

	assert(row <= obj->dim && col <= obj->dim);
	while (row > obj->blk_dims[rm]) {
		row -= obj->blk_dims[rm];
		rm++;
	}
	while (col > obj->blk_dims[cm]) {
		col -= obj->blk_dims[cm];
		cm++;
	}

	if (rm == cm) {
		mat_double dm = obj->m[rm];
		switch (dm.type) {
		case MAT_DENSE:
			val = dm.data[row-1 + (col-1)*dm.row];
			break;
		case MAT_DENSE_DIAG:
			if (row == col) val = dm.data[row-1];
			break;
		case MAT_SPARSE: {
			int i;
			for (i=dm.irow[row-1]; i<dm.irow[row]; i++) {
				if (dm.icol[i-1] == col) {
					val = dm.data[i-1];
					break;
				}
			}
			break;
		}
		case MAT_SPARSE_DIAG:
			if (row == col) {
				int cc = dm.irow[row-1];
				if (dm.irow[row]-cc != 0) val = dm.data[cc-1];
			}
			break;
		default:
			fprintf(stderr,"Error: blk_dm_getelem - unknown matrix %d type %d\n",rm,dm.type);
			exit(1);
		}
	}

	return val;
}

complx blk_cm_getelem(blk_mat_complx *obj, int row, int col)
{
	// row and col are 1-based
	int rm=0, cm=0;
	complx val = Cnull;

	assert(row <= obj->dim && col <= obj->dim);
	while (row > obj->blk_dims[rm]) {
		row -= obj->blk_dims[rm];
		rm++;
	}
	while (col > obj->blk_dims[cm]) {
		col -= obj->blk_dims[cm];
		cm++;
	}

	if (rm == cm) {
		mat_complx dm = obj->m[rm];
		switch (dm.type) {
		case MAT_DENSE:
			val = dm.data[row-1 + (col-1)*dm.row];
			break;
		case MAT_DENSE_DIAG:
			if (row == col) val = dm.data[row-1];
			break;
		case MAT_SPARSE: {
			int i;
			for (i=dm.irow[row-1]; i<dm.irow[row]; i++) {
				if (dm.icol[i-1] == col) {
					val = dm.data[i-1];
					break;
				}
			}
			break;
		}
		case MAT_SPARSE_DIAG:
			if (row == col) {
				int cc = dm.irow[row-1];
				if (dm.irow[row]-cc != 0) val = dm.data[cc-1];
			}
			break;
		default:
			fprintf(stderr,"Error: blk_cm_getelem - unknown matrix %d type %d\n",rm,dm.type);
			exit(1);
		}
	}

	return val;
}


void blk_dm_change_basis(blk_mat_double *dest, blk_mat_double *from, Sim_info *sim)
{
	// base do ktere se transformuje je uz definovana v cilove matici
	// zdrojova matice ma bud spravnou blokovou strukturu, nebo jeden blok
	// cilova matice zrovna tak...

	int i;

	if (dest->basis == from->basis) {
		if (dest->Nblocks == from->Nblocks) {
			for (i=0; i<dest->Nblocks; i++) {
				assert(dest->blk_dims[i] == from->blk_dims[i]);
				dm_copy(&(dest->m[i]),&(from->m[i]));
			}
		} else if (dest->Nblocks == 1) {
			blk_dm_copy_2(dest->m,from);
		} else {
			assert(from->Nblocks == 1);
			//blk_dm_copy_3(dest,from->m);
			blk_dm_zero(dest);
			blk_dm_multod_extract(dest,from->m,1);
		}
	} else {
		int *permvec = sim->perm_table[from->basis + dest->basis*sim->Nbasis];
		assert(LEN(permvec) == sim->matdim);
		int dim, sft = 1;
		mat_double *dm;
		for (i=0; i<dest->Nblocks; i++) {
			dim = dest->blk_dims[i];
			dm = dest->m+i;
			assert(dim == sim->dims_table[dest->basis][i+1] || dim == sim->matdim);
			assert(dim == dm->row);
			switch (dm->type) {
			case MAT_DENSE: {
				int r, c, ro, co;
				for (r=0; r<dim; r++) {
					ro = permvec[sft+r];
					for (c=0; c<dim; c++) {
						co = permvec[sft+c];
						dm->data[r+c*dim] = blk_dm_getelem(from,ro,co);
					}
				}
				break;
			}
			case MAT_DENSE_DIAG: {
				int r, ro;
				for (r=0; r<dim; r++) {
					ro = permvec[sft+r];
					dm->data[r] = blk_dm_getelem(from,ro,ro);
				}
				break;
			}
			case MAT_SPARSE: {
				int k, l, ro, co, nnz = 0;
				const uint64_t nnzmax = (uint64_t)((uint64_t)dim*dim*(1.0-SPARSITY));
				dm_change_nnz(dm,nnzmax);
				dm->irow[0] = 1;
				for (k=0; k<dim; k++) {
					dm->irow[k+1] = dm->irow[k];
					ro = permvec[sft+k];
					for (l=0; l<dim; l++) {
						co = permvec[sft+l];
						double val = blk_dm_getelem(from,ro,co);
						if (fabs(val) < SPARSE_TOL) continue;
						dm->icol[nnz] = l+1;
						dm->data[nnz] = val;
						dm->irow[k+1]++;
						nnz++;
						if (nnz >= nnzmax) break;
					}
					if (l<dim) break;
				}
				if (nnz >= nnzmax) {
					dm_dense(dm);
					int kk, ll;
					for (kk=k; kk<dim; kk++) {
						ro = permvec[sft+kk];
						for (ll=l; ll<dim; ll++) dm->data[kk+ll*dim] = blk_dm_getelem(from,ro,permvec[sft+ll]);
						l = 0;
					}
				} else {
					assert( k==dim && l==dim );
					dm_change_nnz(dm,dm->irow[dim]-1);
				}
				break;
			}
			case MAT_SPARSE_DIAG: {
				int r, ro, cnt=0;
				double val;
				dm_change_nnz(dm,dim);
				dm->irow[0] = 1;
				for (r=0; r<dim; r++) {
					dm->irow[r+1] = dm->irow[r];
					ro = permvec[sft+r];
					val = blk_dm_getelem(from,ro,ro);
					if (fabs(val) < SPARSE_TOL) continue;
					dm->data[cnt] = val;
					dm->icol[cnt] = r+1;
					(dm->irow[r+1])++;
					cnt++;
				}
				dm_change_nnz(dm,cnt);
				break;
			}
			default:
				fprintf(stderr,"Error: blk_dm_change_basis - unknown dest matrix %d type %d\n",i,dm->type);
				exit(1);
			}
			sft += dim;
		}
	}
	//blk_dm_print(dest,"blk_dm_change_basis DEST");
}

void blk_cm_change_basis(blk_mat_complx *dest, blk_mat_complx *from, Sim_info *sim)
{
	// base do ktere se transformuje je uz definovana v cilove matici
	// zdrojova matice ma bud spravnou blokovou strukturu, nebo jeden blok
	// cilova matice zrovna tak...

	int i;

	if (dest->basis == from->basis) {
		if (dest->Nblocks == from->Nblocks) {
			for (i=0; i<dest->Nblocks; i++) {
				assert(dest->blk_dims[i] == from->blk_dims[i]);
				cm_copy(&(dest->m[i]),&(from->m[i]));
			}
		} else if (dest->Nblocks == 1) {
			blk_cm_copy_2(dest->m,from);
		} else {
			assert(from->Nblocks == 1);
			//blk_dm_copy_3(dest,from->m);
			blk_cm_zero(dest);
			blk_cm_multod_extract(dest,from->m,1);
		}
	} else {
		int *permvec = sim->perm_table[from->basis + dest->basis*sim->Nbasis];
		assert(LEN(permvec) == sim->matdim);
		int dim, sft = 1;
		mat_complx dm;
		for (i=0; i<dest->Nblocks; i++) {
			dim = dest->blk_dims[i];
			dm = dest->m[i];
			assert(dim == sim->dims_table[dest->basis][i+1] || dim == sim->matdim);
			assert(dim == dm.row);
			switch (dm.type) {
			case MAT_DENSE: {
				int r, c, ro, co;
				for (r=0; r<dim; r++) {
					ro = permvec[sft+r];
					for (c=0; c<dim; c++) {
						co = permvec[sft+c];
						dm.data[r+c*dim] = blk_cm_getelem(from,ro,co);
					}
				}
				break;
			}
			case MAT_DENSE_DIAG: {
				int r, ro;
				for (r=0; r<dim; r++) {
					ro = permvec[sft+r];
					dm.data[r] = blk_cm_getelem(from,ro,ro);
				}
				break;
			}
			case MAT_SPARSE: {
				int k, l, ro, co, nnz = 0;
				const uint64_t nnzmax = (uint64_t)((uint64_t)dim*dim*(1.0-SPARSITY));
				cm_change_nnz(&dm,nnzmax);
				dm.irow[0] = 1;
				for (k=0; k<dim; k++) {
					dm.irow[k+1] = dm.irow[k];
					ro = permvec[sft+k];
					for (l=0; l<dim; l++) {
						co = permvec[sft+l];
						complx val = blk_cm_getelem(from,ro,co);
						if (fabs(val.re) < SPARSE_TOL && fabs(val.im)< SPARSE_TOL) continue;
						dm.icol[nnz] = l+1;
						dm.data[nnz] = val;
						dm.irow[k+1]++;
						nnz++;
						if (nnz > nnzmax) break;
					}
					if (l<dim) break;
				}
				if (nnz > nnzmax) {
					cm_dense(&dm);
					int kk, ll;
					for (kk=k; kk<dim; kk++) {
						ro = permvec[sft+kk];
						for (ll=l; ll<dim; ll++) dm.data[kk+ll*dim] = blk_cm_getelem(from,ro,permvec[sft+ll]);
						l = 0;
					}
				} else {
					assert( k==dim && l==dim );
					cm_change_nnz(&dm,dm.irow[dim]-1);
				}
				break;
			}
			case MAT_SPARSE_DIAG: {
				int r, ro, cnt=0;
				cm_change_nnz(&dm,dim);
				dm.irow[0] = 1;
				for (r=0; r<dim; r++) {
					dm.irow[r+1] = dm.irow[r];
					ro = permvec[sft+r];
					complx val = blk_cm_getelem(from,ro,ro);
					if (fabs(val.re) < SPARSE_TOL && fabs(val.im)<SPARSE_TOL) continue;
					dm.data[cnt] = val;
					dm.icol[cnt] = r+1;
					(dm.irow[r+1])++;
					cnt++;
				}
				cm_change_nnz(&dm,cnt);
				break;
			}
			default:
				fprintf(stderr,"Error: blk_cm_change_basis - unknown dest matrix %d type %d\n",i,dm.type);
				exit(1);
			}
			sft += dim;
		}
	}
}

void blk_simtrans(mat_complx *sg, blk_mat_complx *U, Sim_info *sim)
{
	mat_complx *sigma;

	//LARGE_INTEGER tv1, tv2, tickpsec, tv3, tv4;
	//QueryPerformanceFrequency(&tickpsec);
	//QueryPerformanceCounter(&tv1);

	if (sg->basis != U->basis) {
		DEBUGPRINT("\nblk_simtrans: changing basis of sigma\n");
		sigma = cm_change_basis_2(sg,U->basis,sim);
		//cm_print(sigma,"blk_simtrans - sigma");
	} else {
		sigma = sg;
	}

 //cm_print(sigma,"blk_simtrans sigma");
 //blk_cm_print(U,"blk_simtrans prop");


 if (U->Nblocks == 1) {
	simtrans(sigma,U->m);
 } else {
	switch (sigma->type) {
	case MAT_DENSE: {
		int i, j, Nr, Nc, pr = 0, pc;
		mat_complx *cm1, *cm2, *dum = NULL;
		for (i=0; i<U->Nblocks; i++) {
			Nr = U->blk_dims[i];
			cm1 = U->m + i;
			// diagonal blocks; if Nr==1, no change to sigma submatrix
			if (Nr != 1) {
				dum = cm_extract_block(sigma,pr,pr,Nr,Nr);
				if (dum != NULL) {
					switch (cm1->type) {
					case MAT_DENSE:
						//dum = complx_matrix(Nr,Nr,MAT_DENSE,0,sigma->basis);
						assert(dum->type == MAT_DENSE);
						cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,Nr,Nr,Nr,&Cunit,sigma->data+pr+pr*sigma->row,sigma->row,cm1->data,Nr,&Cnull,dum->data,Nr);
						cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,Nr,Nr,Nr,&Cunit,cm1->data,Nr,dum->data,Nr,&Cnull,sigma->data+pr+pr*sigma->row,sigma->row);
						//free_complx_matrix(dum);
						break;
					case MAT_DENSE_DIAG: {
						int i2,j2;
						complx *dA, *dA2, *dU, *dU2;
						dU = cm1->data;
						for (i2=0; i2<Nr; i2++) { /* columns */
							dA = dA2 = sigma->data+pr+sigma->row*pr + (sigma->row+1)*i2; // sit on the diagonal element
							dU2 = cm1->data + i2;
							for (j2=i2+1; j2<Nr; j++) { /* rows, lower triangle */
								dU2++;
								dA++;
								dA2 += sigma->row;
								double re = dA->re*dU->re + dA->im*dU->im;
								double im = dA->im*dU->re - dA->re*dU->im;
								dA->re = re*dU2->re - im*dU2->im;
								dA->im = im*dU2->re + re*dU2->im;
								re = dA2->re*dU->re - dA2->im*dU->im;
								im = dA2->im*dU->re + dA2->re*dU->im;
								dA2->re = re*dU2->re + im*dU2->im;
								dA2->im = im*dU2->re - re*dU2->im;
							}
							dU++;
						}
						break; }
					case MAT_SPARSE: {
						char trans='N', matdescra[]="GNNFUU";
						mat_complx *dum2;
						//dum = complx_matrix(Nr,Nr,MAT_DENSE,0,sigma->basis);
						assert(dum->type == MAT_DENSE);
#ifdef INTEL_MKL
						mkl_zcsrmm(&trans, &Nr, &Nr, &Nr, (MKL_Complex16*)&Cunit, matdescra, (MKL_Complex16*)(cm1->data), cm1->icol, cm1->irow, &(cm1->irow[1]), (MKL_Complex16*)(sigma->data+pr), &(sigma->row), (MKL_Complex16*)&Cnull, (MKL_Complex16*)(dum->data), &Nr);
#else
						fprintf(stderr,"Error: blk_simtrans - sparse blocks without MKL support\n");
						exit(1);
#endif
						dum2 = cm_dup(cm1);
						cm_dense(dum2);
						cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,Nr,Nr,Nr,&Cunit,dum->data,Nr,dum2->data,Nr,&Cnull,sigma->data+pr,sigma->row);
						//free_complx_matrix(dum);
						free_complx_matrix(dum2);
						break; }
					case MAT_SPARSE_DIAG:
						fprintf(stderr,"Error: blk_simtrans - diag blocks - prop matrix %d of type MAT_SPARSE_DIAG\n",i);
						exit(1);
					default :
						fprintf(stderr,"Error: blk_simtrans - diag blocks - unknown matrix %d type %d\n",i, cm1->type);
						exit(1);
					}
					free_complx_matrix(dum);
					dum = NULL;
				}
			}
			// off-diagonal blocks, first above and then below diagonal
			pc = pr + Nr;
			for (j=i+1; j<U->Nblocks; j++) {
				Nc = U->blk_dims[j];
				cm2 = U->m + j;
				if (Nr == 1) {
					if (Nc == 1) {
						double re1, im1;
						complx *z1 = cm1->data;
						complx *z2 = cm2->data;
						complx *z3 = sigma->data+pr+pc*sigma->row;
						re1 = z3->re*z1->re - z3->im*z1->im;
						im1 = z3->re*z1->im + z3->im*z1->re;
						z3->re = re1*z2->re + im1*z2->im;
						z3->im = im1*z2->re - re1*z2->im;
						z3 = sigma->data+pc+pr*sigma->row;
						re1 = z3->re*z1->re + z3->im*z1->im;
						im1 = z3->im*z1->re - z3->re*z1->im;
						z3->re = re1*z2->re - im1*z2->im;
						z3->im = im1*z2->re + re1*z2->im;
					} else {
						switch (cm2->type) {
						case MAT_DENSE: {
							complx *zz = (complx*)malloc(2*Nc*sizeof(complx));
							int i2;
							for (i2=0; i2<Nc; i2++) {
								cblas_zdotc_sub(Nc,cm2->data+i2,Nc,sigma->data+pr+pc*sigma->row,sigma->row,zz+i2);
								cblas_zdotu_sub(Nc,sigma->data+pc+pr*sigma->row,1,cm2->data+i2,Nc,zz+i2+Nc);
							}
							for (i2=0; i2<Nc; i2++) {
								sigma->data[pr+(pc+i2)*sigma->row] = Cmul(zz[i2],cm1->data[0]);
								sigma->data[pc+i2+pr*sigma->row] = Cmul(zz[i2+Nc],Conj(cm1->data[0]));
							}
							free(zz);
							break;	}
						case MAT_DENSE_DIAG: {
							complx *z1 = cm1->data;
							complx *z2 = cm2->data;
							complx *z3 = sigma->data+pr+pc*sigma->row;
							complx *z4 = sigma->data+pc+pr*sigma->row;
							int i2;
							double re1, im1, re2;
							for (i2=0; i2<Nc; i2++) {
								re1 = z1->re*z2->re + z1->im*z2->im;
								im1 = z1->im*z2->re - z1->re*z2->im;
								re2 = z3->re;
								z3->re = re1*z3->re - im1*z3->im;
								z3->im = re1*z3->im + im1*re2;
								re2 = z4->re;
								z4->re = re1*z4->re + im1*z4->im;
								z4->im = re1*z4->im - im1*re2;
								z2++;
								z3 += sigma->row;
								z4++;
							}
							break; }
						case MAT_SPARSE: {
							complx *zz = (complx*)calloc(2*Nc,sizeof(complx));
							int i2, j2, n2, c2;
							complx *z3 = sigma->data+pr+pc*sigma->row;
							complx *z4 = sigma->data+pc+pr*sigma->row;
							complx zn;
							for (i2=0; i2<Nc; i2++) {
								n2 = cm2->irow[i2+1] - cm2->irow[i2];
								for (j2=0; j2<n2; j2++) {
									c2 = cm2->icol[cm2->irow[i2]-1+j2] - 1;
									zn = cm2->data[cm2->irow[i2]-1+j2];
									zz[c2].re += z3->re*zn.re + z3->im*zn.im;
									zz[c2].im += z3->im*zn.re - z3->re*zn.im;
									zz[c2+Nc].re += z4->re*zn.re - z4->im*zn.im;
									zz[c2+Nc].im += z4->im*zn.re + z4->re*zn.im;
								}
								z3 += sigma->row;
								z4++;
							}
							for (i2=0; i2<Nc; i2++) {
								sigma->data[pr+(pc+i2)*sigma->row] = Cmul(zz[i2],cm1->data[0]);
								sigma->data[pc+i2+pr*sigma->row] = Cmul(zz[i2+Nc],Conj(cm1->data[0]));
							}
							free(zz);
							break; }
						case MAT_SPARSE_DIAG:
						default:
							fprintf(stderr,"Error: blk_simtrans - unsupported combination of matrix types, place 1\n");
							exit(1);
						}
					}
				} else {
					// Nr != 1
					if (Nc == 1) {
						switch (cm1->type) {
						case MAT_DENSE: {
							complx *zz = (complx*)malloc(2*Nr*sizeof(complx));
							int i2;
							for (i2=0; i2<Nr; i2++) {
								cblas_zdotu_sub(Nr,sigma->data+pr+pc*sigma->row,1,cm1->data+i2,Nr,zz+i2);
								cblas_zdotc_sub(Nr,cm1->data+i2,Nr,sigma->data+pc+pr*sigma->row,sigma->row,zz+i2+Nr);
							}
							for (i2=0; i2<Nr; i2++) {
								sigma->data[pr+i2+pc*sigma->row] = Cmul(zz[i2],Conj(cm2->data[0]));
								sigma->data[pc+(pr+i2)*sigma->row] = Cmul(zz[i2+Nr],cm2->data[0]);
							}
							free(zz);
							break;	}
						case MAT_DENSE_DIAG: {
							complx *z1 = cm1->data;
							complx *z2 = cm2->data;
							complx *z3 = sigma->data+pr+pc*sigma->row;
							complx *z4 = sigma->data+pc+pr*sigma->row;
							int i2;
							double re1, im1, re2;
							for (i2=0; i2<Nr; i2++) {
								re1 = z1->re*z2->re + z1->im*z2->im;
								im1 = z1->im*z2->re - z1->re*z2->im;
								re2 = z3->re;
								z3->re = re1*z3->re - im1*z3->im;
								z3->im = re1*z3->im + im1*re2;
								re2 = z4->re;
								z4->re = re1*z4->re + im1*z4->im;
								z4->im = re1*z4->im - im1*re2;
								z1++;
								z3++;
								z4 += sigma->row;
							}
							break; }
						case MAT_SPARSE: {
							complx *zz = (complx*)calloc(2*Nr,sizeof(complx));
							int i2, j2, n2, c2;
							complx *z3 = sigma->data+pr+pc*sigma->row;
							complx *z4 = sigma->data+pc+pr*sigma->row;
							complx zn;
							for (i2=0; i2<Nr; i2++) {
								n2 = cm1->irow[i2+1] - cm1->irow[i2];
								for (j2=0; j2<n2; j2++) {
									c2 = cm1->icol[cm1->irow[i2]-1+j2] - 1;
									zn = cm1->data[cm1->irow[i2]-1+j2];
									zz[c2].re += z3->re*zn.re - z3->im*zn.im;
									zz[c2].im += z3->im*zn.re + z3->re*zn.im;
									zz[c2+Nc].re += z4->re*zn.re + z4->im*zn.im;
									zz[c2+Nc].im += z4->im*zn.re - z4->re*zn.im;
								}
								z3++;
								z4 += sigma->row;
							}
							for (i2=0; i2<Nc; i2++) {
								sigma->data[pr+i2+pc*sigma->row] = Cmul(zz[i2],Conj(cm2->data[0]));
								sigma->data[pc+(pr+i2)*sigma->row] = Cmul(zz[i2+Nc],cm2->data[0]);
							}
							free(zz);
							break; }
						case MAT_SPARSE_DIAG:
						default:
							fprintf(stderr,"Error: blk_simtrans - unsupported combination of matrix types, place 2\n");
							exit(1);
						}
					} else {

						//above diagonal
						// IMPROVE: calculation could be skipped if block in sigma is empty!!!
						//          (in the same way as for diagonal blocks of sigma)
						dum = complx_matrix(Nr,Nc,MAT_DENSE,0,sigma->basis);
						switch (cm1->type) {
						case MAT_DENSE:
							cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,Nr,Nc,Nr,&Cunit,cm1->data,Nr,sigma->data+pr+pc*sigma->row,sigma->row,&Cnull,dum->data,Nr);
							break;
						case MAT_DENSE_DIAG: {
							int i2, j2;
							complx *z2 = dum->data;
							complx *z3 = sigma->data+pr+pc*sigma->row;
							for (i2=0; i2<Nc; i2++) {
								complx *z1 = cm1->data;
								for (j2=0; j2<Nr; j2++) {
									z2->re = z1->re*z3->re - z1->im*z3->im;
									z2->im = z1->im*z3->re + z1->re*z3->im;
									z1++;
									z2++;
									z3++;
								}
								z3 += sigma->row - Nr;
							}
							break; }
						case MAT_SPARSE: {
#ifdef INTEL_MKL
							char trans='N', matdescra[]="GNNFUU";
							mkl_zcsrmm(&trans, &Nr, &Nc, &Nr, (MKL_Complex16*)&Cunit, matdescra, (MKL_Complex16*)(cm1->data), cm1->icol, cm1->irow, &(cm1->irow[1]), (MKL_Complex16*)(sigma->data+pr+pc*sigma->row), &(sigma->row), (MKL_Complex16*)&Cnull, (MKL_Complex16*)(dum->data), &Nr);
#else
							fprintf(stderr,"Error: blk_simtrans - unsupported sparse operation, place 7A\n");
							exit(1);
#endif
							break; }
						default:
							fprintf(stderr,"Error: blk_simtrans - unsupported combination of matrix types, place 3A\n");
							exit(1);
						}
						switch (cm2->type) {
						case MAT_DENSE:
							cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,Nr,Nc,Nc,&Cunit,dum->data,Nr,cm2->data,Nc,&Cnull,sigma->data+pr+pc*sigma->row,sigma->row);
							break;
						case MAT_DENSE_DIAG: {
							complx *z1 = dum->data;
							complx *z2 = cm2->data;
							complx *z3 = sigma->data+pr+pc*sigma->row;
							int i2, j2;
							for (i2=0; i2<Nc; i2++) {
								for (j2=0; j2<Nr; j2++) {
									z3->re = z1->re*z2->re + z1->im*z2->im;
									z3->im = z1->im*z2->re - z1->re*z2->im;
									z1++;
									z3++;
								}
								z2++;
								z3 += sigma->row - Nr;
							}
							break; }
						case MAT_SPARSE: {
							mat_complx *dum2;
							dum2 = cm_dup(cm2);
							cm_dense(dum2);
							cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,Nr,Nc,Nc,&Cunit,dum->data,Nr,dum2->data,Nc,&Cnull,sigma->data+pr+pc*sigma->row,sigma->row);
							free_complx_matrix(dum2);
							break; }
						default:
							fprintf(stderr,"Error: blk_simtrans - unsupported combination of matrix types, place 4\n");
							exit(1);
						}
						// below diagonal
						dum->row = Nc;
						dum->col = Nr;
						switch (cm1->type) {
						case MAT_DENSE:
							cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,Nc,Nr,Nr,&Cunit,sigma->data+pc+pr*sigma->row,sigma->row,cm1->data,Nr,&Cnull,dum->data,Nc);
							break;
						case MAT_DENSE_DIAG: {
							int i2, j2;
							complx *z1 = cm1->data;
							complx *z2 = dum->data;
							complx *z3 = sigma->data+pc+pr*sigma->row;
							for (i2=0; i2<Nr; i2++) {
								for (j2=0; j2<Nc; j2++) {
									z2->re =  z1->re*z3->re + z1->im*z3->im;
									z2->im = -z1->im*z3->re + z1->re*z3->im;
									z2++;
									z3++;
								}
								z1++;
								z3 += sigma->row - Nc;
							}
							break; }
						case MAT_SPARSE: {
							mat_complx *dum2 = cm_dup(cm1);
							cm_dense(dum2);
							cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,Nc,Nr,Nr,&Cunit,sigma->data+pc+pr*sigma->row,sigma->row,dum2->data,Nr,&Cnull,dum->data,Nc);
							free_complx_matrix(dum2);
							break; }
						default:
							fprintf(stderr,"Error: blk_simtrans - unsupported combination of matrix types, place 3B\n");
							exit(1);
						}
						switch (cm2->type) {
						case MAT_DENSE:
							cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,Nc,Nr,Nc,&Cunit,cm2->data,Nc,dum->data,Nc,&Cnull,sigma->data+pc+pr*sigma->row,sigma->row);
							break;
						case MAT_DENSE_DIAG: {
							complx *z1 = dum->data;
							complx *z3 = sigma->data+pc+pr*sigma->row;
							int i2, j2;
							for (i2=0; i2<Nr; i2++) {
								complx *z2 = cm2->data;
								for (j2=0; j2<Nc; j2++) {
									z3->re = z1->re*z2->re - z1->im*z2->im;
									z3->im = z1->im*z2->re + z1->re*z2->im;
									z1++;
									z2++;
									z3++;
								}
								z3 += sigma->row - Nc;
							}
							break; }
						case MAT_SPARSE: {
#ifdef INTEL_MKL
							char trans='N', matdescra[]="GNNFUU";
							mkl_zcsrmm(&trans, &Nc, &Nr, &Nc, (MKL_Complex16*)&Cunit, matdescra, (MKL_Complex16*)(cm2->data), cm2->icol, cm2->irow, &(cm2->irow[1]), (MKL_Complex16*)(dum->data), &Nc, (MKL_Complex16*)&Cnull, (MKL_Complex16*)(sigma->data+pc+pr*sigma->row), &(sigma->row));
#else
							fprintf(stderr,"Error: blk_simtrans - unsupported sparse operation, place 6\n");
							exit(1);
#endif
							break; }
						default:
							fprintf(stderr,"Error: blk_simtrans - unsupported combination of matrix types, place 5\n");
							exit(1);
						}
						free_complx_matrix(dum);
					} // end else from if Nc == 1, with Nr != 1
				} // end else from if Nr==1
				pc += Nc;
			} // end for-loop over off-diagonal blocks
			pr += Nr;
		} //
		break;}
	case MAT_DENSE_DIAG:
		// sigma is diagonal; if sparse regime, fall through this to sparse case
	  if (!sim->sparse)	{
		int i;
		for (i=0; i<U->Nblocks; i++) {
			if (U->blk_dims[i] != 1 &&  U->m[i].type != MAT_DENSE_DIAG && U->m[i].type != MAT_SPARSE_DIAG ) break;
		}
		if (i == U->Nblocks) break; // if U is all diagonal, then there is no change to sigma
		int N, pos = 0;
		mat_complx *dum, *cm, *sigmanew;
		sigmanew = complx_matrix(sigma->row, sigma->col, MAT_DENSE, 0,sigma->basis);
		cm_zero(sigmanew);
		for (i=0; i<U->Nblocks; i++) {
			N = U->blk_dims[i];
			if (N == 1) {
				sigmanew->data[pos+pos*sigma->row] = sigma->data[pos];
			} else {
				cm = U->m + i;
				switch (cm->type) {
				case MAT_DENSE: {
					int i2, j2;
					dum = complx_matrix(N,N,MAT_DENSE,0,sigma->basis);
					complx *z1 = cm->data;
					for (i2=0; i2<N; i2++) {
						complx *z2 = dum->data + i2;
						complx *z3 = sigma->data + pos+i2;
						for (j2=0; j2<N; j2++) {
							z2->re = z1->re*z3->re + z1->im*z3->im;
							z2->im = z1->re*z3->im - z1->im*z3->re;
							z1++;
							z2 += N;
						}
					}
					cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,N,N,N,&Cunit,cm->data,N,dum->data,N,&Cnull,sigmanew->data+pos+pos*sigma->row,sigmanew->row);
					free_complx_matrix(dum);
					break; }
				case MAT_DENSE_DIAG: {
					// no change to sigma, just copy
					int i2;
					complx *z1 = sigma->data + pos;
					complx *z2 = sigmanew->data + pos + pos*sigmanew->row;
					for (i2=0; i2<N; i2++) {
						*z2 = *z1;
						z1++;
						z2 += sigmanew->row + 1;
					}
					break; }
				case MAT_SPARSE: {
					dum = cm_adjoint(cm);
					int i2, j2, nc2;
					complx *z1 = dum->data;
					complx *z2 = sigma->data;
					double re1;
					for (i2=0; i2<N; i2++) {
						nc2 = dum->irow[i2+1] - dum->irow[i2];
						for (j2=0; j2<nc2; j2++) {
							re1 = z1->re;
							z1->re = re1*z2->re - z1->im*z2->im;
							z1->im = re1*z2->im + z1->im*z2->re;
							z1++;
						}
						z2++;
					}
					mat_complx *dum2 = simpson_zcsrmultcsr(cm,dum);
					z1 = dum2->data;
					z2 = sigmanew->data + pos + pos*sigmanew->row;
					int *ic1 = dum2->icol;
					for (i2=0; i2<N; i2++) {
						nc2 = dum2->irow[i2+1] - dum2->irow[i2];
						for (j2=0; j2<nc2; j2++) {
							z2[ (*ic1 - 1)*sigmanew->row ] = *z1;
							z1++;
							ic1++;
						}
						z2++;
					}
					free_complx_matrix(dum2);
					free_complx_matrix(dum);
					break; }
				default:
					fprintf(stderr,"Error: blk_simtrans - unsupported combination of matrix types, place 11\n");
					exit(1);
				}
			}
			pos += N;
		}
		cm_swap_innards_and_destroy(sigma,sigmanew);
		break;}
	case MAT_SPARSE:
	case MAT_SPARSE_DIAG: {
		// sigma is sparse
		int i, j, Nr, Nc, pr, pc, N1 = 0;
		mat_complx *cm1, *cm2, *dum = NULL, *sigmanew = NULL;
		for (i=0; i<U->Nblocks; i++) {
			if (U->blk_dims[i] == 1) N1++;
		}
		if (N1>0) {
			//diagonal blocks: if Nr==1, no change to sigma submatrix, copy it directly
			sigmanew = complx_matrix(sigma->row,sigma->col,MAT_SPARSE,N1,sigma->basis);
			sigmanew->irow[0] = 1;
			pr = 0;
			pc = 0;
			for (i=0; i<U->Nblocks; i++) {
				Nr = U->blk_dims[i];
				if (Nr == 1) {
					sigmanew->irow[pr+1] = sigmanew->irow[pr];
					pr++;
					sigmanew->irow[pr]++;
					sigmanew->data[pc] = cm_getelem(sigma,pr,pr);
					sigmanew->icol[pc] = pr;
					pc++;
				} else {
					for (j=0; j<Nr; j++) {
						sigmanew->irow[pr+1] = sigmanew->irow[pr];
						pr++;
					}
				}
			}
			//cm_print(sigmanew,"DEBUG blk_simtrans test: is this diagonal sigmanew?");
		} else {
			sigmanew = complx_matrix(sigma->row,sigma->col,MAT_SPARSE,1,sigma->basis);
			cm_zero(sigmanew);
		}
		assert(sigmanew != NULL);
		pr = 0;
		for (i=0; i<U->Nblocks; i++) {
			Nr = U->blk_dims[i];
			cm1 = U->m + i;
			// diagonal blocks
			if (Nr != 1)  {
				dum = cm_extract_block(sigma, pr, pr, Nr, Nr);
				if (dum != NULL) {
					//cm_print(dum,"blk_simtrans extract diag block");
					simtrans(dum,cm1);
					cm_addto_block(sigmanew,dum,pr,pr);
					free_complx_matrix(dum);
					dum = NULL;
				}
			}
			// off-diagonal blocks
			pc = pr + Nr;
			for (j=i+1; j<U->Nblocks; j++) {
				Nc = U->blk_dims[j];
				cm2 = U->m + j;
				// above diagonal
				dum = cm_extract_block(sigma,pr,pc,Nr,Nc);
				if (dum !=  NULL) {
					// transform, resize, add to sigmanew
					//printf("blk_simtrans extract block from (%i,%i)",pr,pc);
					//cm_print(dum," up");
					//cm_print(cm1,"cm1");
					//cm_print(cm2,"cm2");
					if (Nr == 1) {
						cm_mulc(dum,cm1->data[0]);
					} else {
						cm_multo_rev(dum,cm1);
					}
					if (Nc == 1) {
						cm_mulc(dum,Conj(cm2->data[0]));
					} else {
						mat_complx *tmp = cm_adjoint(cm2);
						cm_multo(dum,tmp);
						free_complx_matrix(tmp);
					}
					cm_addto_block(sigmanew,dum,pr,pc);
					//cm_print(dum,"transformed");
					free_complx_matrix(dum);
					dum = NULL;
				}
				// below diagonal
				dum = cm_extract_block(sigma,pc, pr, Nc, Nr);
				if (dum != NULL) {
					//printf("blk_simtrans extract block from (%i,%i)",pc,pr);
					//cm_print(dum," down");
					// transform, resize, add to sigmanew
					if (Nc == 1) {
						cm_mulc(dum,cm2->data[0]);
					} else {
						cm_multo_rev(dum,cm2);
					}
					if (Nr == 1) {
						cm_mulc(dum,Conj(cm1->data[0]));
					} else {
						mat_complx *tmp = cm_adjoint(cm1);
						cm_multo(dum,tmp);
						free_complx_matrix(tmp);
					}
					cm_addto_block(sigmanew,dum,pc,pr);
					//cm_print(dum,"transformed");
					free_complx_matrix(dum);
					dum = NULL;
				}
				pc += Nc;
			}
			pr += Nr;
		}	// end for loop over block rows
		cm_swap_innards_and_destroy(sigma,sigmanew);
		break; }
	default:
		fprintf(stderr,"ERROR: blk_simtrans - unknown sigma matrix type %d\n",sigma->type);
		exit(1);
	}
 }

 //cm_print(sigma,"blk_simtrans sigma RESULT");
 //exit(1);

	if (sg != sigma) {
		int base = sg->basis;
		//cm_print(sigma,"blk_simtrans - sigma");
		DEBUGPRINT("blk_simtrans: sigma: change basis back to %d\n",base);
		cm_swap_innards_and_destroy(sg,sigma);
		sigma = cm_change_basis_2(sg,base,sim);
		cm_swap_innards_and_destroy(sg,sigma);
	}

	//QueryPerformanceCounter(&tv2);
	//printf("timing blk_simtrans time: %.9f\n",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);

}

/***
 * extract block from matrix
 *   r0, c0 are zero-based!!!
 *   if empty block, then NULL returned
 ***/
mat_complx * cm_extract_block(mat_complx *cm, int r0, int c0, int Nr, int Nc)
{
	mat_complx *res = NULL;
	int i, j;

	switch (cm->type) {
	case MAT_DENSE: {
		double dbl;
		// test NNZ
		if (Nr == 1) {
			dbl = cblas_dzasum(Nc,cm->data+r0+c0*cm->row,cm->row);
		} else {
			dbl = 0;
			for (i=0; i<Nc; i++) {
				dbl += cblas_dzasum(Nr,cm->data+r0+(c0+i)*cm->row,1);
			}
		}
		if (dbl > TINY) {
			// extract
			res = complx_matrix(Nr,Nc,MAT_DENSE,0,cm->basis);
			if (Nr == 1) {
				cblas_zcopy(Nc,cm->data+r0+c0*cm->row,cm->row,res->data,1);
			} else {
				for (i=0; i<Nc; i++) {
					cblas_zcopy(Nr,cm->data+r0+(c0+i)*cm->row,1,res->data+i*Nr,1);
				}
			}
		}
		break;	}
	case MAT_DENSE_DIAG: {
		double dbl;
		if (r0 == c0) {
			assert(Nr == Nc);
			dbl = cblas_dzasum(Nr,cm->data+r0,1);
			if (dbl > TINY) {
				// extract
				res = complx_matrix(Nr,Nc,MAT_DENSE_DIAG,0,cm->basis);
				memcpy(res->data,cm->data+r0,Nr*sizeof(complx));
			}
		}
		break;	}
	case MAT_SPARSE: {
		// test NNZ
		int nnz = 0;
		for (i=0; i<Nr; i++) {
			for (j=cm->irow[r0+i]; j<cm->irow[r0+i+1]; j++) {
				if (cm->icol[j-1]>c0 && c0+Nc-cm->icol[j-1]>=0) nnz++;
			}
		}
		if (nnz != 0) {
			// extract
			res = complx_matrix(Nr,Nc,MAT_SPARSE,nnz,cm->basis);
			complx *z = res->data;
			int *ic = res->icol;
			res->irow[0] = 1;
			for (i=0; i<Nr; i++) {
				res->irow[i+1] = res->irow[i];
				for (j=cm->irow[r0+i]; j<cm->irow[r0+i+1]; j++) {
					if (cm->icol[j-1]>c0 && c0+Nc-cm->icol[j-1]>=0) {
						res->irow[i+1]++;
						*z = cm->data[j-1];
						*ic = cm->icol[j-1] - c0;
						z++;
						ic++;
					}
				}
			}
		}
		break;	}
	case MAT_SPARSE_DIAG: {
		if (r0 == c0) {
			assert(Nr == Nc);
			for (i=0; i<Nr; i++) {
				if (cm->irow[r0+i+1] - cm->irow[r0+i] != 0) break;
			}
			if (i < Nr) {
				// extract
				res = complx_matrix(Nr,Nc,MAT_DENSE_DIAG,0,cm->basis);
				for (i=0; i<Nr; i++) {
					j = cm->irow[r0+i];
					if (cm->irow[r0+i+1] - j != 0) {
						res->data[i] = cm->data[j-1];
					} else {
						res->data[i] = Cnull;
					}
				}
			}
		}
		break;	}
	default:
		fprintf(stderr,"Error: cm_extract_block - unknown matrix type %d\n",cm->type);
		exit(1);
	}

	return res;
}

/***
 * extract block from matrix
 *   r0, c0 are zero-based!!!
 *   if empty block, then NULL returned
 ***/
mat_double * dm_extract_block(mat_double *cm, int r0, int c0, int Nr, int Nc)
{
	mat_double *res = NULL;
	int i, j;

	switch (cm->type) {
	case MAT_DENSE: {
		double dbl;
		// test NNZ
		if (Nr == 1) {
			dbl = cblas_dasum(Nc,cm->data+r0+c0*cm->row,cm->row);
		} else {
			dbl = 0;
			for (i=0; i<Nc; i++) {
				dbl += cblas_dasum(Nr,cm->data+r0+(c0+i)*cm->row,1);
			}
		}
		if (dbl > TINY) {
			// extract
			res = double_matrix(Nr,Nc,MAT_DENSE,0,cm->basis);
			if (Nr == 1) {
				cblas_dcopy(Nc,cm->data+r0+c0*cm->row,cm->row,res->data,1);
			} else {
				for (i=0; i<Nc; i++) {
					cblas_dcopy(Nr,cm->data+r0+(c0+i)*cm->row,1,res->data+i*Nr,1);
				}
			}
		}
		break;	}
	case MAT_DENSE_DIAG: {
		double dbl;
		if (r0 == c0) {
			assert(Nr == Nc);
			dbl = cblas_dasum(Nr,cm->data+r0,1);
			if (dbl > TINY) {
				// extract
				res = double_matrix(Nr,Nc,MAT_DENSE_DIAG,0,cm->basis);
				memcpy(res->data,cm->data+r0,Nr*sizeof(double));
			}
		}
		break;	}
	case MAT_SPARSE: {
		// test NNZ
		int nnz = 0;
		for (i=0; i<Nr; i++) {
			for (j=cm->irow[r0+i]; j<cm->irow[r0+i+1]; j++) {
				if (cm->icol[j-1]>c0 && c0+Nc-cm->icol[j-1]>=0) nnz++;
			}
		}
		if (nnz != 0) {
			// extract
			res = double_matrix(Nr,Nc,MAT_SPARSE,nnz,cm->basis);
			double *z = res->data;
			int *ic = res->icol;
			res->irow[0] = 1;
			for (i=0; i<Nr; i++) {
				res->irow[i+1] = res->irow[i];
				for (j=cm->irow[r0+i]; j<cm->irow[r0+i+1]; j++) {
					if (cm->icol[j-1]>c0 && c0+Nc-cm->icol[j-1]>=0) {
						res->irow[i+1]++;
						*z = cm->data[j-1];
						*ic = cm->icol[j-1] - c0;
						z++;
						ic++;
					}
				}
			}
		}
		break;	}
	case MAT_SPARSE_DIAG: {
		if (r0 == c0) {
			assert(Nr == Nc);
			for (i=0; i<Nr; i++) {
				if (cm->irow[r0+i+1] - cm->irow[r0+i] != 0) break;
			}
			if (i < Nr) {
				// extract
				res = double_matrix(Nr,Nc,MAT_DENSE_DIAG,0,cm->basis);
				for (i=0; i<Nr; i++) {
					j = cm->irow[r0+i];
					if (cm->irow[r0+i+1] - j != 0) {
						res->data[i] = cm->data[j-1];
					} else {
						res->data[i] = 0.0;
					}
				}
			}
		}
		break;	}
	default:
		fprintf(stderr,"Error: dm_extract_block - unknown matrix type %d\n",cm->type);
		exit(1);
	}

	return res;
}

/****
 * adds smaller sub-block matrix to full sized matrix at position r0,c0 (zero based)
 */
void cm_addto_block(mat_complx *cm, mat_complx *bm, int r0, int c0)
{
	// limit usage to this:
	//assert(cm->type == MAT_SPARSE || cm->type == MAT_DENSE);

	int i, j;

	if (cm->type == MAT_SPARSE) {
		switch (bm->type) {
		case MAT_DENSE: {
			int nnz = cm_nnz(bm);
			mat_complx *dum = complx_matrix(cm->row,cm->col,MAT_SPARSE,nnz,cm->basis);
			for (i=0; i<=r0; i++) dum->irow[i]=1;
			complx *z = dum->data, *zz;
			int *ic = dum->icol;
			for (i=0; i<bm->row; i++) {
				dum->irow[r0+i+1] = dum->irow[r0+i];
				for (j=0; j<bm->col; j++) {
					zz = bm->data + i + j*bm->row;
					if (zz->re*zz->re + zz->im*zz->im > TINY) {
						*z = *zz;
						*ic = c0 + j + 1;
						dum->irow[r0+i+1]++;
						z++;
						ic++;
					}
				}
			}
			for (i=r0+bm->row+1; i<=dum->row; i++) dum->irow[i] = dum->irow[r0+bm->row];
			cm_addto(cm,dum);
			free_complx_matrix(dum);
			break; }
		case MAT_DENSE_DIAG: {
			int nnz = cm_nnz(bm);
			mat_complx *dum = complx_matrix(cm->row,cm->col,MAT_SPARSE,nnz,cm->basis);
			for (i=0; i<=r0; i++) dum->irow[i]=1;
			complx *z = dum->data, *zz = bm->data;
			int *ic = dum->icol;
			for (i=0; i<bm->row; i++) {
				dum->irow[r0+i+1] = dum->irow[r0+i];
				if (zz->re*zz->re + zz->im*zz->im > TINY) {
					*z = *zz;
					*ic = c0 + i + 1;
					dum->irow[r0+i+1]++;
					z++;
					ic++;
				}
				zz++;
			}
			for (i=r0+bm->row+1; i<=dum->row; i++) dum->irow[i] = dum->irow[r0+bm->row];
			cm_addto(cm,dum);
			free_complx_matrix(dum);
			break; }
		case MAT_SPARSE:
		case MAT_SPARSE_DIAG: {
			int *irow = (int*)malloc((cm->row+1)*sizeof(int));
			for (i=0; i<=r0; i++) irow[i] = 1;
			for (i=1; i<=bm->row; i++) {
				irow[r0+i] =  bm->irow[i];
			}
			for (i=r0+bm->row+1; i<=cm->row; i++) irow[i] = irow[r0+bm->row];
			for (i=0; i<bm->irow[bm->row]-1; i++) bm->icol[i] += c0;
			free(bm->irow);
			bm->irow = irow;
			bm->row = cm->row;
			bm->col = cm->col;
			cm_addto(cm,bm);
			break; }
		default:
			fprintf(stderr,"Error: cm_addto_block - unknown sub-matrix type %d\n",bm->type);
			exit(1);
		}
	} else if (cm->type == MAT_DENSE) {
		switch (bm->type) {
		case MAT_DENSE:
			for (i=0; i<bm->col; i++) cblas_zaxpy(bm->row,&Cunit,&(bm->data[i*bm->row]),1,&(cm->data[r0+(c0+i)*cm->row]),1);
			break;
		case MAT_DENSE_DIAG:
			for (i=0; i<bm->row; i++) {
				cm->data[r0+i+(c0+i)*cm->row].re += bm->data[i].re;
				cm->data[r0+i+(c0+i)*cm->row].im += bm->data[i].im;
			}
			break;
		case MAT_SPARSE:
		case MAT_SPARSE_DIAG: {
			int *ic = bm->icol;
			complx *z = bm->data;
			for (i=0; i<bm->row; i++) {
				int n = bm->irow[i+1] - bm->irow[i];
				for (j=0; j<n; j++) {
					cm->data[r0+i+(c0+(*ic)-1)*cm->row].re += z->re;
					cm->data[r0+i+(c0+(*ic)-1)*cm->row].im += z->im;
					z++;
					ic++;
				}
			}
			break; }
		default:
			fprintf(stderr,"Error: cm_addto_block - unknown sub-matrix type %d\n",bm->type);
			exit(1);
		}
	} else {
		fprintf(stderr,"Error: cm_addto_block - wrong cm matrix type (%s)\n",matrix_type(cm->type));
		exit(1);
	}
}

/* IMPROVE: now it creates adjoint block matrix and uses blk_simtrans
 *         more efficient would be directly coded routine
 */
void blk_simtrans_adj(mat_complx *sg, blk_mat_complx *U, Sim_info *sim)
{
	mat_complx *sigma;

	if (sg->basis != U->basis) {
		DEBUGPRINT("blk_simtrans_adj: changing basis of sigma\n");
		sigma = cm_change_basis_2(sg,U->basis,sim);
	} else {
		sigma = sg;
	}

	if (U->Nblocks == 1) {
		simtrans_adj(sigma,U->m);
	} else {
		//DEBUGPRINT("blk_simtrans_adj point HA\n");
		blk_mat_complx *blkm = blk_cm_adjoint(U);
		//DEBUGPRINT("blk_simtrans_adj point HUHU\n");
		blk_simtrans(sigma,blkm,sim);
		free_blk_mat_complx(blkm);
	}
	if (sg != sigma) {
		int base = sg->basis;
		cm_swap_innards_and_destroy(sg,sigma);
		sigma = cm_change_basis_2(sg,base,sim);
		cm_swap_innards_and_destroy(sg,sigma);
	}

}

blk_mat_complx * blk_cm_adjoint(blk_mat_complx *m)
{
	blk_mat_complx *res;
	mat_complx *rmx, *mx;
	int i, n;

	res = (blk_mat_complx*)malloc(sizeof(blk_mat_complx));
	res->dim = m->dim;
	res->Nblocks = m->Nblocks;
	res->basis = m->basis;
	res->m = (mat_complx*)malloc(res->Nblocks*sizeof(mat_complx));
	if (res->m == NULL) {
		fprintf(stderr,"Error: blk_cm_adjoint can not allocate block-matrix pointers\n");
		exit(1);
	}
	res->blk_dims = (int*)malloc(res->Nblocks*sizeof(int));
	if (res->blk_dims == NULL) {
		fprintf(stderr,"Error: blk_cm_adjoint can not allocate blk_dims\n");
		exit(1);
	}
	for (i=0; i<res->Nblocks; i++) {

		res->blk_dims[i] = n = m->blk_dims[i];
		rmx = res->m + i;
		mx = m->m + i;
		//DEBUGPRINT("blk_cm_adjoint block %d ...",i);
		if (n == 1) {
			//DEBUGPRINT("(dim=%d)",n);
			create_sqmat_complx(rmx,1,MAT_DENSE,1,m->basis);
			rmx->data[0].re = mx->data[0].re;
			rmx->data[0].im = -mx->data[0].im;
		} else {
				//DEBUGPRINT("(dim=%d)",n);
			switch (mx->type) {
			case MAT_DENSE: {
				int j;
				create_sqmat_complx(rmx,n,MAT_DENSE,1,m->basis);
				for (j=0; j<n; j++) {
					cblas_zcopy(n, &(mx->data[j]),n,&(rmx->data[j*n]),1);
					cblas_dscal(n,-1.0,(double*)(&(rmx->data[j*n]))+1,2);
				}
				break; }
			case MAT_DENSE_DIAG:
				create_sqmat_complx(rmx,n,MAT_DENSE_DIAG,1,m->basis);
				memcpy(rmx->data,mx->data,n*sizeof(complx));
				cblas_dscal(n,-1.0,(double*)(rmx->data)+1,2);
				break;
			case MAT_SPARSE: {
				int len = mx->irow[n] - 1;
				create_sqmat_complx(rmx,n,MAT_SPARSE,len,m->basis);
				int *dumrow, *dumcol, *ic = mx->icol;
				dumrow = (int*)malloc((n+1+len)*sizeof(int));
				dumcol = dumrow+n+1;
				int ii, jj, nn=0;
				complx *z;
				memset(dumrow,0,(mx->col+1)*sizeof(int));
				for (ii=0; ii<n; ii++) {
					int nnn = mx->irow[ii+1] - mx->irow[ii];
					for (jj=0; jj<nnn; jj++) {
						dumrow[*ic]++;
						dumcol[nn] = ii+1;
						ic++;
						nn++;
					}
				}
				rmx->irow[0] = dumrow[0] = 1;
				for (ii=1; ii<=n; ii++) {
					rmx->irow[ii] = rmx->irow[ii-1] + dumrow[ii];
					dumrow[ii] = rmx->irow[ii];
				}
				z = mx->data;
				for (ii=0; ii<len; ii++) {
					jj = dumrow[mx->icol[ii]-1]-1;
					rmx->data[jj].re = z->re;
					rmx->data[jj].im = -z->im;
					rmx->icol[jj] = dumcol[ii];
					z++;
					dumrow[mx->icol[ii]-1]++;
				}
				free(dumrow);
				break;}
			case MAT_SPARSE_DIAG: {
				int nnz = mx->irow[n] - 1;
				create_sqmat_complx(rmx,n,MAT_SPARSE_DIAG,nnz,m->basis);
				memcpy(rmx->data, mx->data, nnz*sizeof(complx));
				memcpy(rmx->icol, mx->icol, nnz*sizeof(MKL_INT));
				memcpy(rmx->irow, mx->irow, (n+1)*sizeof(MKL_INT));
				cblas_dscal(n,-1.0,(double*)(rmx->data)+1,2);
				break;}
			default:
				fprintf(stderr,"Error: blk_cm_adjoint - unknown input matrix %d type %d\n",i,mx->type);
				exit(1);
			}
		}
		//DEBUGPRINT(" ... done!\n");
	}

	return res;
}

blk_mat_complx * blk_cm_dup(blk_mat_complx *m)
{
	// nezapomenout kopirovat i info o bazi
	blk_mat_complx *res;
	mat_complx *rmx, *mx;
	int i, n;

	res = (blk_mat_complx*)malloc(sizeof(blk_mat_complx));
	res->dim = m->dim;
	res->Nblocks = m->Nblocks;
	res->basis = m->basis;
	res->m = (mat_complx*)malloc(res->Nblocks*sizeof(mat_complx));
	if (res->m == NULL) {
		fprintf(stderr,"Error: blk_cm_dup can not allocate block-matrix pointers\n");
		exit(1);
	}
	res->blk_dims = (int*)malloc(res->Nblocks*sizeof(int));
	if (res->blk_dims == NULL) {
		fprintf(stderr,"Error: blk_cm_dup can not allocate blk_dims\n");
		exit(1);
	}
	for (i=0; i<res->Nblocks; i++) {
		res->blk_dims[i] = n = m->blk_dims[i];
		rmx = res->m + i;
		mx = m->m + i;
		if (n == 1) {
			create_sqmat_complx(rmx,1,MAT_DENSE,1,m->basis);
			rmx->data[0].re = mx->data[0].re;
			rmx->data[0].im = mx->data[0].im;
		} else {
			switch (mx->type) {
			case MAT_DENSE:
				create_sqmat_complx(rmx,n,MAT_DENSE,1,m->basis);
				memcpy(rmx->data,mx->data,n*n*sizeof(complx));
				break;
			case MAT_DENSE_DIAG:
				create_sqmat_complx(rmx,n,MAT_DENSE_DIAG,1,m->basis);
				memcpy(rmx->data,mx->data,n*sizeof(complx));
				break;
			case MAT_SPARSE: {
				int len = mx->irow[n] -1;
				create_sqmat_complx(rmx,n,MAT_SPARSE,len,m->basis);
				memcpy(rmx->data, mx->data, len*sizeof(complx));
				memcpy(rmx->icol, mx->icol, len*sizeof(int));
				memcpy(rmx->irow, mx->irow, (n+1)*sizeof(int));
				break; }
			case MAT_SPARSE_DIAG: {
				int len = mx->irow[n] -1;
				create_sqmat_complx(rmx,n,MAT_SPARSE_DIAG,len,m->basis);
				memcpy(rmx->data, mx->data, len*sizeof(complx));
				memcpy(rmx->icol, mx->icol, len*sizeof(int));
				memcpy(rmx->irow, mx->irow, (n+1)*sizeof(int));
				break; }
			default:
				fprintf(stderr,"Error: blk_cm_dup - unknown input matrix %d type %d\n",i,mx->type);
				exit(1);
			}
		}
	}
	return res;
}

void blk_cm_copy(blk_mat_complx *dest, blk_mat_complx *from)
{
	assert(dest != NULL);
	assert(dest->dim == from->dim);

	int i;

	//if (dest->basis == from->basis) {
	//	assert(dest->Nblocks == from->Nblocks);
	//	for (i=0; i<from->Nblocks; i++) {
	//		cm_copy(dest->m+i,from->m+i);
	//	}
	//} else
	if (dest->Nblocks == from->Nblocks) {
		dest->basis = from->basis;
		for (i=0; i<dest->Nblocks; i++) {
			dest->blk_dims[i] = from->blk_dims[i];
			cm_copy(dest->m+i,from->m+i);
		}
	} else {
		mat_complx *mx, *rmx;
		int n, *ddd;
		for (i=0; i<dest->Nblocks; i++) {
			mx = dest->m + i;
			if (mx->data != NULL) free(mx->data);
			if (mx->irow != NULL) free(mx->irow);
			if (mx->icol != NULL) free(mx->icol);
		}
		mx = (mat_complx*)realloc(dest->m, from->Nblocks*sizeof(mat_complx));
		if (mx == NULL) {
			fprintf(stderr, "Error: blk_cm_copy - problem when reallocating m\n");
			exit(1);
		}
		dest->m = mx;
		ddd = (int*)realloc(dest->blk_dims, from->Nblocks*sizeof(int));
		if (ddd == NULL) {
			fprintf(stderr, "Error: blk_cm_copy - problem when reallocating blk_dims\n");
			exit(1);
		}
		dest->blk_dims = ddd;
		dest->Nblocks = from->Nblocks;
		dest->basis = from->basis;
		for (i=0; i<dest->Nblocks; i++) {
			rmx = dest->m + i;
			mx = from->m + i;
			n = dest->blk_dims[i] = from->blk_dims[i];
			switch (mx->type) {
			case MAT_DENSE:
				create_sqmat_complx(rmx,n,MAT_DENSE,1,mx->basis);
				memcpy(rmx->data,mx->data,n*n*sizeof(complx));
				break;
			case MAT_DENSE_DIAG:
				create_sqmat_complx(rmx,n,MAT_DENSE_DIAG,1,mx->basis);
				memcpy(rmx->data,mx->data,n*sizeof(complx));
				break;
			case MAT_SPARSE: {
				int len = mx->irow[n] -1;
				create_sqmat_complx(rmx,n,MAT_SPARSE,len,mx->basis);
				memcpy(rmx->data,mx->data,len*sizeof(complx));
				memcpy(rmx->icol,mx->icol,len*sizeof(int));
				memcpy(rmx->irow,mx->irow,(n+1)*sizeof(int));
				break; }
			case MAT_SPARSE_DIAG: {
				int len = mx->irow[n] -1;
				create_sqmat_complx(rmx,n,MAT_SPARSE_DIAG,len,mx->basis);
				memcpy(rmx->data,mx->data,len*sizeof(complx));
				memcpy(rmx->icol,mx->icol,len*sizeof(int));
				memcpy(rmx->irow,mx->irow,(n+1)*sizeof(int));
				break; }
			default:
				fprintf(stderr,"Error: blk_cm_copy - unknown input matrix %d type %d\n",i,mx->type);
				exit(1);
			}
		}
	}
}

/* NOT complete; intended to convert block matrix to sparse matrix
mat_complx * blk_cm_sparse(blk_mat_complx *blkm, double tol)
{
	mat_complx *res, *bm;
	int i, j, NN, dim = blkm->dim;

	for (i=0; i<blkm->Nblocks; i++) {
		if (blkm->blk_dims[i] > 1 ) {
			if (blkm->m[i].type == MAT_DENSE || blkm->m[i].type == MAT_SPARSE) break;
		}
	}
	if (i == blkm->Nblocks) {
		// block matrix is diagonal
		res = complx_matrix(dim,dim,MAT_DENSE_DIAG,0,blkm->basis);
		complx *zres = res->data;
		for (i=0; i<blkm->Nblocks; i++) {
			NN = blkm->blk_dims[i];
			assert(NN==1 || blkm->m[i].type == MAT_DENSE_DIAG);
			memcpy(zres,blkm->m[i].data,NN*sizeof(complx));
			zres += NN;
		}
	} else {
		int nnz = 0, pr = 0;
		const int nnzmax = (int)(dim*dim*(1.0-SPARSITY));
		res = complx_matrix(dim,dim,MAT_SPARSE,nnzmax,blkm->basis);
		int *ir = res->irow;
		int *ic = res->icol;
		complx *zz = res->data;
		ir[0] = 1;
		for (i=0; i<blkm->Nblocks; i++) {
			NN = blkm->blk_dims[i];
			bm = blkm->m + i;
			ir[pr+1] = ir[pr];
			if (NN == 1) {
				complx *z = bm->data;
				if (fabs(z->re)>tol || fabs(z->im)>tol) {
					ir[pr+1]++;
					*ic = pr+1;
					*zz = *z;
					nnz++;
					ic++;
					zz++;
				}
				pr++;
			} else {
				// analyze and copy this block, check if nnz < nnzmax or switch to dense matrix
			}
		}
	}
}
*******************************/


blk_mat_complx * blk_cm_mul(blk_mat_complx *blkm1, blk_mat_complx *blkm2)
{
	blk_mat_complx *res;
	mat_complx *rm, *mm;
	int i, NN;

	// currently implemented only with these restrictions
	assert(blkm1->basis == blkm2->basis);
	assert(blkm1->Nblocks == blkm2->Nblocks);

	res = (blk_mat_complx *)malloc(sizeof(blk_mat_complx));
	res->basis = blkm1->basis;
	res->Nblocks = blkm1->Nblocks;
	res->dim = blkm1->dim;
	res->blk_dims = (int*)malloc(res->Nblocks*sizeof(int));
	res->m = (mat_complx*)malloc(res->Nblocks*sizeof(mat_complx));
	if (res->blk_dims == NULL || res->m == NULL) {
		fprintf(stderr,"Error: blk_cm_mul can not allocate m/blk_dims\n");
		exit(1);
	}
	for (i=0; i<blkm1->Nblocks; i++) {
		NN = res->blk_dims[i] = blkm1->blk_dims[i];
		rm = res->m + i;
		if (NN == 1) {
			create_sqmat_complx(rm,1,MAT_DENSE_DIAG,0,blkm1->basis);
			rm->data[0] = Cmul(blkm1->m[i].data[0],blkm2->m[i].data[0]);
		} else {
			mm = cm_mul(blkm1->m + i, blkm2->m + i);
			rm->row = mm->row;
			rm->col = mm->col;
			rm->basis = mm->basis;
			rm->type = mm->type;
			rm->irow = mm->irow;
			rm->icol = mm->icol;
			rm->data = mm->data;
			free(mm);
		}
	}

	//blk_cm_print(blkm1,"blk_cm_print blkm1");
	//blk_cm_print(blkm2,"blk_cm_print blkm2");
	//blk_cm_print(res,"blk_cm_print res");
	//exit(1);
	return res;
}

mat_complx * blk_cm_mul_v1(blk_mat_complx *blkm, mat_complx *cm)
{
	assert(cm->basis == blkm->basis);
	assert(cm->row == blkm->dim && cm->col == blkm->dim);
	mat_complx *res, *bm;
	int i, j;
	const int dim = blkm->dim;

	switch (cm->type) {
	case MAT_DENSE: {
		int NN, pr=0, pc;
		res = complx_matrix(dim,dim,MAT_DENSE,0,blkm->basis);
		for (i=0; i<blkm->Nblocks; i++) {
			bm = blkm->m+i;
			NN = blkm->blk_dims[i];
			if (NN == 1) {
				complx *z = bm->data;
				complx *zz = cm->data + pr;
				complx *zres = res->data + pr;
				for (j=0; j<dim; j++) {
					zres->re = z->re*zz->re - z->im*zz->im;
					zres->im = z->re*zz->im + z->im*zz->re;
					zres += dim;
					zz += dim;
				}
			} else {
				pc = 0;
				int NNc;
				switch (bm->type) {
				case MAT_DENSE:
					for (j=0; j<blkm->Nblocks; j++) {
						NNc = blkm->blk_dims[j];
						if (NNc == 1) {
							cblas_zgemv(CblasColMajor,CblasNoTrans,NN,NN,&Cunit,bm->data,NN,cm->data+pr+pc*dim,1,&Cnull,res->data+pr+pc*dim,1);
						} else {
							cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,NN,NNc,NN,&Cunit,bm->data,NN,cm->data+pr+pc*dim,dim,&Cnull,res->data+pr+pc*dim,dim);
						}
						pc += NNc;
					}
					break;
				case MAT_DENSE_DIAG: {
					int j2;
					for (j=0; j<NN; j++) {
						complx *z = bm->data + j;
						complx *zz = cm->data + pr + j;
						complx *zres = res->data + pr + j;
						for (j2=0; j2<dim; j2++) {
							zres->re = z->re*zz->re - z->im*zz->im;
							zres->im = z->re*zz->im + z->im*zz->re;
							zres += dim;
							zz += dim;
						}
					}
					break; }
				case MAT_SPARSE: {
#ifdef INTEL_MKL
					char trans='N', matdescra[]="GNNFUU";
					for (j=0; j<blkm->Nblocks; j++) {
						NNc = blkm->blk_dims[j];
						if (NNc == 1) {
							mkl_zcsrmv(&trans, &NN, &NN, (MKL_Complex16*)&Cunit, matdescra, (MKL_Complex16*)(bm->data), bm->icol, bm->irow, &(bm->irow[1]), (MKL_Complex16*)(cm->data+pr+pc*dim), (MKL_Complex16*)&Cnull, (MKL_Complex16*)(res->data+pr+pc*dim));
						} else {
							mkl_zcsrmm(&trans, &NN, &NNc, &NN, (MKL_Complex16*)&Cunit, matdescra, (MKL_Complex16*)(bm->data), bm->icol, bm->irow, &(bm->irow[1]), (MKL_Complex16*)(cm->data+pr+pc*dim), &dim, (MKL_Complex16*)&Cnull, (MKL_Complex16*)(res->data+pr+pc*dim), &dim);
						}
					}
#else
					fprintf(stderr,"Error: blk_cm_mul_v1 - sparse blocks without MKL support\n");
					exit(1);
#endif
					break; }
				case MAT_SPARSE_DIAG:
					fprintf(stderr,"Error: blk_cm_mul_v1 - submatrix %d has unsupported type %d (place AAA)\n",i,bm->type);
					break;
				default:
					fprintf(stderr,"Error: blk_cm_mul_v1 - undefined matrix %d type %d\n",i,bm->type);
					exit(1);
				}
			}
			pr += NN;
		}
		break; }
	case MAT_DENSE_DIAG: {
		int NN;
		for (i=0; i<blkm->Nblocks; i++) {
			if (blkm->blk_dims[i] > 1 ) {
				if (blkm->m[i].type == MAT_DENSE || blkm->m[i].type == MAT_SPARSE) break;
			}
		}
		if (i == blkm->Nblocks) {
			res = complx_matrix(dim,dim,MAT_DENSE_DIAG,0,blkm->basis);
			complx *z = cm->data;
			complx *zres = res->data;
			for (i=0; i<blkm->Nblocks; i++) {
				NN = blkm->blk_dims[i];
				complx *zz = blkm->m[i].data;
				if (NN == 1) {
					zres->re = z->re*zz->re - z->im*zz->im;
					zres->im = z->re*zz->im + z->im*zz->re;
					zres++;
					z++;
				} else {
					for (j=0; j<NN; j++) {
						zres->re = z->re*zz->re - z->im*zz->im;
						zres->im = z->re*zz->im + z->im*zz->re;
						zres++;
						zz++;
						z++;
					}
				}
			}
		} else {
			res = complx_matrix(dim,dim,MAT_DENSE,0,blkm->basis);
			cm_zero(res);
			complx *z = cm->data;
			int pr = 0;
			for (i=0; i<blkm->Nblocks; i++) {
				NN = blkm->blk_dims[i];
				bm = blkm->m + i;
				if (NN == 1) {
					complx *zres = res->data + pr+pr*dim;
					complx *zz = bm->data;
					zres->re = z->re*zz->re - z->im*zz->im;
					zres->im = z->re*zz->im + z->im*zz->re;
					z++;
				} else {
					for (j=0; j<NN; j++) {
						cblas_zaxpy(NN,z,bm->data+j*NN,1,res->data+pr+(pr+j)*dim,1);
						z++;
					}
				}
				pr += NN;
			}
		}
		break; }
	case MAT_SPARSE:
		fprintf(stderr,"Error: blk_cm_mul_v1 - sparse algebra not implemented yet (place uu)\n");
		exit(1);
		//res = blk_cm_sparse(blkm,SPARSE_TOL);
		//cm_multo(res,cm);
		break;
	case MAT_SPARSE_DIAG:
		fprintf(stderr,"Error: blk_cm_mul_v1 - MAT_SPARSE_DIAG not implemented, place rrr\n");
		exit(1);
		break;
	default:
		fprintf(stderr,"Error: blk_cm_mul_v1 - undefined matrix type (cm -> %d)\n",cm->type);
		exit(1);
	}

	return res;
}

mat_complx * blk_cm_mul_v2(mat_complx *cm, blk_mat_complx *blkm)
{
	assert(cm->basis == blkm->basis);
	assert(cm->row == blkm->dim && cm->col == blkm->dim);
	mat_complx *res, *bm;
	int i, j;
	const int dim = blkm->dim;

	switch (cm->type) {
	case MAT_DENSE: {
		int NN, pr, pc=0;
		res = complx_matrix(dim,dim,MAT_DENSE,0,blkm->basis);
		for (i=0; i<blkm->Nblocks; i++) {
			bm = blkm->m+i;
			NN = blkm->blk_dims[i];
			if (NN == 1) {
				complx *z = bm->data;
				complx *zz = cm->data + pc*dim;
				complx *zres = res->data + pc*dim;
				for (j=0; j<dim; j++) {
					zres->re = z->re*zz->re - z->im*zz->im;
					zres->im = z->re*zz->im + z->im*zz->re;
					zres++;
					zz++;
				}
			} else {
				pr = 0;
				int NNr;
				switch (bm->type) {
				case MAT_DENSE:
					for (j=0; j<blkm->Nblocks; j++) {
						NNr = blkm->blk_dims[j];
						// kdyz NNr==1 tak to je radka*blk_matice=radka, stejne jako matice*blk_matice=matice
						cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,NNr,NN,NNr,&Cunit,cm->data+pr+pc*dim,dim,bm->data,NN,&Cnull,res->data+pr+pc*dim,dim);
						pr += NNr;
					}
					break;
				case MAT_DENSE_DIAG: {
					int j2;
					for (j=0; j<NN; j++) {
						complx *z = bm->data + j;
						complx *zz = cm->data + (pc+j)*dim;
						complx *zres = res->data + (pc+j)*dim;
						for (j2=0; j2<dim; j2++) {
							zres->re = z->re*zz->re - z->im*zz->im;
							zres->im = z->re*zz->im + z->im*zz->re;
							zres++;
							zz++;
						}
					}
					break; }
				case MAT_SPARSE: {
					fprintf(stderr,"Error: blk_mcm_mul_v2 - sparse submatrix %d, unfinished code, place HUHA\n",i);
					exit(1);
					break; }
				case MAT_SPARSE_DIAG:
					fprintf(stderr,"Error: blk_cm_mul_v2 - submatrix %d has unsupported type %d (place AAA)\n",i,bm->type);
					break;
				default:
					fprintf(stderr,"Error: blk_cm_mul_v2 - undefined matrix %d type %d\n",i,bm->type);
					exit(1);
				}
			}
			pc += NN;
		}
		break; }
	case MAT_DENSE_DIAG: {
		int NN;
		for (i=0; i<blkm->Nblocks; i++) {
			if (blkm->blk_dims[i] > 1 ) {
				if (blkm->m[i].type == MAT_DENSE || blkm->m[i].type == MAT_SPARSE) break;
			}
		}
		if (i == blkm->Nblocks) {
			res = complx_matrix(dim,dim,MAT_DENSE_DIAG,0,blkm->basis);
			complx *z = cm->data;
			complx *zres = res->data;
			for (i=0; i<blkm->Nblocks; i++) {
				NN = blkm->blk_dims[i];
				assert(blkm->m[i].type != MAT_SPARSE_DIAG);
				complx *zz = blkm->m[i].data;
				if (NN == 1) {
					zres->re = z->re*zz->re - z->im*zz->im;
					zres->im = z->re*zz->im + z->im*zz->re;
					zres++;
					z++;
				} else {
					for (j=0; j<NN; j++) {
						zres->re = z->re*zz->re - z->im*zz->im;
						zres->im = z->re*zz->im + z->im*zz->re;
						zres++;
						zz++;
						z++;
					}
				}
			}
		} else {
			res = complx_matrix(dim,dim,MAT_DENSE,0,blkm->basis);
			cm_zero(res);
			complx *z = cm->data;
			int pr = 0;
			for (i=0; i<blkm->Nblocks; i++) {
				NN = blkm->blk_dims[i];
				bm = blkm->m + i;
				if (NN == 1) {
					complx *zres = res->data + pr+pr*dim;
					complx *zz = bm->data;
					zres->re = z->re*zz->re - z->im*zz->im;
					zres->im = z->re*zz->im + z->im*zz->re;
					z++;
				} else {
					assert(bm->type == MAT_DENSE);  // bm can be sparse or diag, not implemented
					for (j=0; j<NN; j++) {
						cblas_zaxpy(NN,z,bm->data+j,NN,res->data+pr+j+pr*dim,dim);
						z++;
					}
				}
				pr += NN;
			}
		}
		break; }
	case MAT_SPARSE:
		fprintf(stderr,"Error: blk_cm_mul_v2 - sparse algebra not implemented yet (place uu)\n");
		exit(1);
		//res = blk_cm_sparse(blkm,SPARSE_TOL);
		//cm_multo(res,cm);
		break;
	case MAT_SPARSE_DIAG:
		fprintf(stderr,"Error: blk_cm_mul_v2 - MAT_SPARSE_DIAG not implemented, place rrr\n");
		exit(1);
		break;
	default:
		fprintf(stderr,"Error: blk_cm_mul_v2 - undefined matrix type (cm -> %d)\n",cm->type);
		exit(1);
	}

	return res;
}

void update_propagator(blk_mat_complx *U, blk_mat_complx *dU, Sim_info *sim, Sim_wsp *wsp)
{
	// vychazi z cm_multo_rev, ale musi dat pozor na baze
	// kdyz wsp == NULL, pak se neohlizi na wsp->Uisunit (vzdy prinasobit), jinak jo a meni ho

	int i;
	mat_complx *mU, *mdU;


	if (wsp == NULL) {
		// this is used only within _pulse_simple  and _delay_simple with no base change
		// always update U with dU (no check for unity of U)
		assert(U->basis == dU->basis && U->Nblocks == dU->Nblocks);
		DEBUGPRINT("update_propagator 1\n");
		for (i=0; i<U->Nblocks; i++) {
			mU = U->m + i;
			mdU = dU->m + i;
			if (U->blk_dims[i] == 1) {
				mU->data[0] = Cmul(mU->data[0],mdU->data[0]);
			} else {
				cm_multo_rev(mU,mdU);
			}
		}
	} else {
		if (wsp->Uisunit) {
			DEBUGPRINT("update_propagator 2A\n");
			// just copy
			blk_change_structure_complx2(U,dU);
			blk_cm_copy(U,dU);
			wsp->Uisunit = 0;
		} else {
			if (U->basis == dU->basis) {
				DEBUGPRINT("update_propagator 2B1\n");
				if (U->Nblocks == dU->Nblocks) {
					DEBUGPRINT("update_propagator 2B1a\n");
					for (i=0; i<U->Nblocks; i++) {
						if (U->blk_dims[i] == 1) {
							U->m[i].data[0] = Cmul(U->m[i].data[0],dU->m[i].data[0]);
						} else {
							cm_multo_rev(U->m+i, dU->m+i);
						}
					}
				} else if (U->Nblocks == 1) {
					DEBUGPRINT("update_propagator 2B1b\n");
					// U bude mit jeden blok
					assert(dU->Nblocks > 1);
					mU = blk_cm_mul_v1(dU,U->m); // block matrix times full matrix
					cm_swap_innards_and_destroy(U->m,mU);
				} else {
					DEBUGPRINT("update_propagator 2B1c\n");
					assert(dU->Nblocks == 1);
					// U bude mit jeden blok
					mU = blk_cm_mul_v2(dU->m,U);  // full matrix times block matrix
					for (i=0; i<U->Nblocks; i++) {
						mdU = U->m + i;
						if (mdU->data != NULL) free(mdU->data);
						if (mdU->irow != NULL) free(mdU->irow);
						if (mdU->icol != NULL) free(mdU->icol);
					}
					free(U->m);
					U->m = mU;
					int *ddd = (int*)realloc(U->blk_dims, sizeof(int));
					if (ddd == NULL) {
						fprintf(stderr, "Error: update_propagator - problem when reallocating blk_dims, place A\n");
						exit(1);
					}
					U->blk_dims = ddd;
					U->blk_dims[0] = U->dim;
					U->Nblocks = 1;
				}
			} else {
				// ruzne baze
				DEBUGPRINT("update_propagator 2B2\n");
				if (U->Nblocks == 1) {
					assert(dU->Nblocks != 1);
					if (U->basis == 0 && U->m->type == MAT_DENSE_DIAG) {
						// predelat U do baze dU a prenasobit tou diag. matici
						mU = cm_change_basis_2(U->m,dU->basis,sim);
						assert(mU->type == MAT_DENSE_DIAG);
						int NN, pp = 0;
						complx *zstart = mU->data;
						for (i=0; i<dU->Nblocks; i++) {
							NN = dU->blk_dims[i];
							mdU = dU->m + i;
							if (NN == 1) {
								mdU->data[0] = Cmul(mdU->data[0],zstart[pp]);
							} else {
								mU->col = mU->row = NN;
								mU->data = zstart + pp;
								cm_multo(mdU,mU);
							}
							pp += NN;
						}
						mU->data = zstart;
						free_complx_matrix(mU);
						// vysledek je ted v dU, je treba U<->dU a zrejme i zmenit strukturu dU (aby byla porad stejna)
						blk_change_structure_complx2(U,dU);
						mU = U->m; U->m = dU->m; dU->m = mU;
					} else {
						// predelat dU do base U (bud 0 nebo 11111) a jednim blokem
						blk_mat_complx *blkdum = create_blk_mat_complx_copy(U);
						blk_cm_change_basis(blkdum,dU,sim);
						assert(blkdum->Nblocks == 1);
						cm_multo_rev(U->m,blkdum->m);
						free_blk_mat_complx(blkdum);
					}
				} else if (dU->Nblocks == 1) {
					if (dU->basis == 0 && dU->m->type == MAT_DENSE_DIAG) {
						//predelat dU do baze U a prenasobit tou diag. matici
						mdU = cm_change_basis_2(dU->m,U->basis,sim);
						assert(mdU->type == MAT_DENSE_DIAG);
						int NN, pp = 0;
						complx *zstart = mdU->data;
						for (i=0; i<U->Nblocks; i++) {
							NN = U->blk_dims[i];
							mU = U->m + i;
							if (NN == 1) {
								mU->data[0] = Cmul(mU->data[0],zstart[pp]);
							} else {
								mdU->col = mdU->row = NN;
								mdU->data = zstart + pp;
								cm_multo_rev(mU,mdU);
							}
							pp += NN;
						}
						mdU->data = zstart;
						free_complx_matrix(mdU);
					} else {
						// predelat U do base dU (bud 0 nebo 11111) a jednim blokem
						blk_mat_complx *blkdum = create_blk_mat_complx_copy(dU);
						blk_cm_change_basis(blkdum,U,sim);
						assert(blkdum->Nblocks == 1);
						cm_multo_rev(blkdum->m,dU->m);
						// swap U with blkdum and destroy work space stuff
						mU = U->m;
						for (i=0; i<U->Nblocks; i++) {
							destroy_sqmat_complx(mU);
							mU++;
						}
						free(U->m);
						U->m = blkdum->m;
						free(U->blk_dims);
						U->blk_dims = blkdum->blk_dims;
						U->Nblocks = blkdum->Nblocks;
						U->basis = blkdum->basis;
						free(blkdum);
					}
				} else {
					// oba maji vic bloku, zjistovat kompatibilni bazi
					assert(U->basis != 0 && dU->basis != 0);
					int new_basis = U->basis & dU->basis; // bitwise AND
					if (new_basis == 0) {
						// there is no possible common blocking
						// switch to basis 0 or 11111 with one block
						new_basis = sim->Hiso->basis;
						if (U->basis == new_basis) {
							// U is in good basis but with blocks
							blk_mat_complx *blkdum_dU = create_blk_mat_complx(dU->dim,1,NULL,(sim->sparse == 2) ? MAT_SPARSE : MAT_DENSE,U->basis);
							blk_cm_change_basis(blkdum_dU,dU,sim);
							mU = blk_cm_mul_v2(blkdum_dU->m,U);  // full matrix times block matrix
							free_blk_mat_complx(blkdum_dU);
							for (i=0; i<U->Nblocks; i++) {
								mdU = U->m + i;
								if (mdU->data != NULL) free(mdU->data);
								if (mdU->irow != NULL) free(mdU->irow);
								if (mdU->icol != NULL) free(mdU->icol);
							}
							free(U->m);
							U->m = mU;
							int *ddd = (int*)realloc(U->blk_dims, sizeof(int));
							if (ddd == NULL) {
								fprintf(stderr, "Error: update_propagator - problem when reallocating blk_dims, place X\n");
								exit(1);
							}
							U->blk_dims = ddd;
							U->blk_dims[0] = U->dim;
							U->Nblocks = 1;
						} else if (dU->basis == new_basis) {
							// dU is in good basis but with blocks
							blk_mat_complx *blkdum_U = create_blk_mat_complx(U->dim,1,NULL,(sim->sparse == 2) ? MAT_SPARSE : MAT_DENSE,dU->basis);
							blk_cm_change_basis(blkdum_U,U,sim);
							mU = blk_cm_mul_v1(dU,blkdum_U->m); // block matrix times full matrix
							free_blk_mat_complx(blkdum_U);
							for (i=0; i<U->Nblocks; i++) {
								mdU = U->m + i;
								if (mdU->data != NULL) free(mdU->data);
								if (mdU->irow != NULL) free(mdU->irow);
								if (mdU->icol != NULL) free(mdU->icol);
							}
							free(U->m);
							U->m = mU;
							int *ddd = (int*)realloc(U->blk_dims, sizeof(int));
							if (ddd == NULL) {
								fprintf(stderr, "Error: update_propagator - problem when reallocating blk_dims, place XX\n");
								exit(1);
							}
							U->blk_dims = ddd;
							U->blk_dims[0] = U->dim;
							U->Nblocks = 1;
						} else {
							blk_mat_complx *blkdum_U = create_blk_mat_complx(U->dim,1,NULL,(sim->sparse == 2) ? MAT_SPARSE : MAT_DENSE,new_basis);
							blk_cm_change_basis(blkdum_U,U,sim);
							blk_mat_complx *blkdum_dU = create_blk_mat_complx_copy(blkdum_U);
							blk_cm_change_basis(blkdum_dU,dU,sim);
							cm_multo_rev(blkdum_U->m, blkdum_dU->m);
							free_blk_mat_complx(blkdum_dU);
							// swap U with blkdum and destroy work space stuff
							mU = U->m;
							for (i=0; i<U->Nblocks; i++) {
								destroy_sqmat_complx(mU);
								mU++;
							}
							free(U->m);
							U->m = blkdum_U->m;
							free(U->blk_dims);
							U->blk_dims = blkdum_U->blk_dims;
							U->Nblocks = blkdum_U->Nblocks;
							U->basis = blkdum_U->basis;
							free(blkdum_U);
						}
					} else {
						// switch to new_basis with some blocking
						if (U->basis == new_basis) {
							assert(dU->basis != new_basis);
							blk_mat_complx *blkdum = create_blk_mat_complx_copy(U);
							blk_cm_change_basis(blkdum,dU,sim);
							for (i=0; i<U->Nblocks; i++) {
								if (U->blk_dims[i] == 1) {
									U->m[i].data[0] = Cmul(U->m[i].data[0],blkdum->m[i].data[0]);
								} else {
									cm_multo_rev(U->m+i, blkdum->m+i);
								}
							}
							free_blk_mat_complx(blkdum);
						} else if (dU->basis == new_basis) {
							blk_mat_complx *blkdum = create_blk_mat_complx_copy(dU);
							blk_cm_change_basis(blkdum,U,sim);
							for (i=0; i<dU->Nblocks; i++) {
								if (dU->blk_dims[i] == 1) {
									blkdum->m[i].data[0] = Cmul(dU->m[i].data[0],blkdum->m[i].data[0]);
								} else {
									cm_multo_rev(blkdum->m+i, dU->m+i);
								}
							}
							// swap U with blkdum and destroy work space stuff
							mU = U->m;
							for (i=0; i<U->Nblocks; i++) {
								destroy_sqmat_complx(mU);
								mU++;
							}
							free(U->m);
							U->m = blkdum->m;
							free(U->blk_dims);
							U->blk_dims = blkdum->blk_dims;
							U->Nblocks = blkdum->Nblocks;
							U->basis = blkdum->basis;
							free(blkdum);
						} else {
							// zmenit oba dva
							blk_mat_complx *blkdum_U = create_blk_mat_complx(U->dim, LEN(sim->dims_table[new_basis]),sim->dims_table[new_basis],(sim->sparse == 2) ? MAT_SPARSE : MAT_DENSE,new_basis);
							blk_cm_change_basis(blkdum_U,U,sim);
							blk_mat_complx *blkdum_dU = create_blk_mat_complx_copy(blkdum_U);
							blk_cm_change_basis(blkdum_dU,dU,sim);
							for (i=0; i<blkdum_U->Nblocks; i++) {
								if (blkdum_U->blk_dims[i] == 1) {
									blkdum_U->m[i].data[0] = Cmul(blkdum_U->m[i].data[0],blkdum_dU->m[i].data[0]);
								} else {
									cm_multo_rev(blkdum_U->m+i, blkdum_dU->m+i);
								}
							}
							free_blk_mat_complx(blkdum_dU);
							// swap U with blkdum and destroy work space stuff
							mU = U->m;
							for (i=0; i<U->Nblocks; i++) {
								destroy_sqmat_complx(mU);
								mU++;
							}
							free(U->m);
							U->m = blkdum_U->m;
							free(U->blk_dims);
							U->blk_dims = blkdum_U->blk_dims;
							U->Nblocks = blkdum_U->Nblocks;
							U->basis = blkdum_U->basis;
							free(blkdum_U);
						}
					}
				}
			}
		}
	}
	//blk_cm_print(wsp->U,"updatePropagator NEW PROP");
}

void blk_prop_real(blk_mat_complx *U, blk_mat_double *ham, double duration, Sim_info *sim)
{
	int i, N;
	mat_complx *dU;
	mat_double *dH;

	//blk_dm_print(ham,"blk_prop_real Hamiltonian in");
	//printf("blk_prop_real: Ham NNZ = %i, norm = %g, ",blk_dm_nnz(ham),dm_normest(ham->m));

	assert(U != NULL);
	DEBUGPRINT("blk_prop_real: U->basis = %d, ham->basis = %d\n",U->basis,ham->basis);
	assert(U->basis == ham->basis); // assume U already has proper block-structure (done in change_basis)
	assert(U->Nblocks == ham->Nblocks);

	for (i=0; i<ham->Nblocks; i++) {
		N = ham->blk_dims[i];
		assert(U->blk_dims[i] == N);
		dU = U->m + i;
		dH = ham->m + i;
		if (N==1) {
			double c = cos(dH->data[0]*duration);
			double s = -sin(dH->data[0]*duration);
			double r = dU->data[0].re;
			dU->data[0].re = r*c - dU->data[0].im*s;
			dU->data[0].im = r*s + dU->data[0].im*c;
		} else {
			//printf("\tblk_prop_real: block %i is %s with %d nnz (dim = %d)\n",i,matrix_type(dH->type),dm_nnz(dH), N);
			//printf("\t               dU was %s \n",matrix_type(dU->type));
			if (sim->sparse) {
				if (dH->type == MAT_DENSE) {
					prop_real(dU,dH,duration,0); // via dsyevr
				} else {
					prop_real(dU,dH,duration,sim->propmethod);
				}
			} else {
				prop_real(dU,dH,duration,sim->propmethod);
			}
			//printf(" blk_prop_real: END dU is %s with %d nnz\n",matrix_type(dU->type),cm_nnz(dU));
		}
	}
	//blk_cm_print(U,"blk_cm_prop_real U");
	//exit(1);
	//printf("prop NNz = %i\n",blk_cm_nnz(U));
}

void blk_prop_complx(blk_mat_complx *U, mat_complx *ham, double dur, Sim_info *sim)
{
	assert(sim->labframe == 1);
	assert(U->Nblocks == 1);
	mat_complx *dum = complx_matrix(ham->row,ham->col,MAT_DENSE_DIAG, 0, ham->basis);

	prop_complx(dum,ham,dur,sim->propmethod);
	cm_multo_rev(U->m,dum);
	free_complx_matrix(dum);

}

void blk_prop_complx_2(blk_mat_complx *U, mat_complx *mtx, double dur, int propmethod)
{
	assert(U != NULL);

	if (U->Nblocks != 1) {
		// squeeze to 1 block
		int i;
		mat_complx *tmp_m;
		for (i=0; i<U->Nblocks; i++) {
			tmp_m = U->m + i;
			if (tmp_m->data != NULL) free(tmp_m->data);
			if (tmp_m->irow != NULL) free(tmp_m->irow);
			if (tmp_m->icol != NULL) free(tmp_m->icol);
		}
		tmp_m = (mat_complx*)realloc(U->m,sizeof(mat_complx));
		if (tmp_m != NULL) {
			U->m = tmp_m;
		} else {
			fprintf(stderr,"Error: blk_prop_complx_2 - can not reallocate U->m\n");
			exit(1);
		}
		int *tmp_dims = (int*)realloc(U->blk_dims,sizeof(int));
		if (tmp_dims != NULL) {
			U->blk_dims = tmp_dims;
		} else {
			fprintf(stderr,"Error: blk_prop_complx_2 - can not reallocate U->blk_dims\n");
			exit(1);
		}
		U->Nblocks = 1;
		U->blk_dims[0] = U->dim = mtx->row;
		create_sqmat_complx(U->m,U->dim,(mtx->type == MAT_SPARSE) ? MAT_SPARSE : MAT_DENSE_DIAG,1,mtx->basis);
	}

	U->basis = mtx->basis;
	assert(U->dim == mtx->row);
	prop_complx(U->m,mtx,dur,propmethod);

}

void blk_dm_multod_extract(blk_mat_double *blkm, mat_double *dm, double dval)
{
	// z mat_double si vytahne prislusne diagonalni bloky a ty pricte do blkm

	assert(blkm->basis == dm->basis);

	int i, j, N, pos = 0;
	mat_double *bm;

	for (i=0; i<blkm->Nblocks; i++) {
		N = blkm->blk_dims[i];
		bm = blkm->m + i;
		if (N == 1) {
			pos++;
			bm->data[0] += dval*dm_getelem(dm,pos,pos);
			continue;
		}
		//DEBUGPRINT("\t blk_multod_extract: (%i) N = %d, %s\n",i,N,matrix_type(bm->type));
		switch (bm->type) {
		case MAT_DENSE:
			switch (dm->type) {
			case MAT_DENSE:
				for (j=0; j<N; j++) {
					cblas_daxpy(N,dval,dm->data+pos+(j+pos)*dm->row,1,bm->data+j*N,1);
				}
				break;
			case MAT_DENSE_DIAG:
				cblas_daxpy(N,dval,dm->data+pos,1,bm->data,N+1);
				break;
			case MAT_SPARSE: {
				int r, c, nc, cs;
				for (r=0; r<N; r++) {
					cs = dm->irow[r+pos];
					nc = dm->irow[r+pos+1] - cs;
					for (j=0; j<nc; j++) {
						c = dm->icol[cs+j-1] - pos - 1;
						if (c>=0 && c<N) bm->data[r+c*N] += dval*dm->data[cs+j-1];
					}
				}
				break; }
			case MAT_SPARSE_DIAG: {
				int cs;
				for (j=0; j<N; j++) {
					cs = dm->irow[j+pos];
					if (dm->irow[j+pos+1] - cs == 1) bm->data[j+j*N] += dval*dm->data[cs-1];
				}
				break; }
			default:
				fprintf(stderr,"Error: blk_dm_multod_extract - unknown from-matrix type %d\n",dm->type);
				exit(1);
			}
			break;
		case MAT_DENSE_DIAG:
			fprintf(stderr,"Error: blk_dm_multod_extract - matrix %d type %d not implemented\n",i,bm->type);
			exit(1);
			break;
		case MAT_SPARSE: {
			mat_double *dmx = dm_extract_block(dm,pos,pos,N,N);
			if (dmx != NULL) {
				dm_multod(bm,dmx,dval);
				free_double_matrix(dmx);
			}
			/***** original coding - START *****
			switch (dm->type) {
			case MAT_DENSE: {
				mat_double *dmx = double_matrix(N,N,MAT_DENSE,1,blkm->basis);
				for (j=0; j<N; j++) {
					cblas_dcopy(N,dm->data+pos+(j+pos)*dm->row,1,dmx->data+j*N,1);
				}
				dm_sparse(dmx,SPARSE_TOL);
				dm_multod(bm,dmx,dval);
				free_double_matrix(dmx);
				break; }
			case MAT_DENSE_DIAG: {
				mat_double *dmx = double_matrix(N,N,MAT_SPARSE_DIAG,N,blkm->basis);
				for (j=0; j<N; j++) {
					dmx->irow[j] = j+1;
					dmx->icol[j] = j+1;
					dmx->data[j] = dm->data[j+pos];
				}
				dmx->irow[N] = N + 1;
				dm_multod(bm,dmx,dval);
				free_double_matrix(dmx);
				break; }
			case MAT_SPARSE: {
				int nzmax = (int)(N*N*(1.0-SPARSITY));
				int r, c, nc, cs, pp = 0;
				double val;
				mat_double *dmx = double_matrix(N,N,MAT_SPARSE,nzmax,blkm->basis);
				dmx->irow[0] = 1;
				for (r=0; r<N; r++) {
					cs = dm->irow[r+pos];
					nc = dm->irow[r+pos+1] - cs;
					dmx->irow[r+1] = dmx->irow[r];
					for (j=0; j<nc; j++) {
						c = dm->icol[cs-1+j] - pos - 1;
						if (c>=0 && c<N) val = dm->data[cs-1+j]; else val=0;
						if (fabs(val)<SPARSE_TOL) continue;
						dmx->irow[r+1]++;
						dmx->icol[pp] = c+1;
						dmx->data[pp] = val;
						pp++;
						if (pp>=nzmax) {
							fprintf(stderr,"Error: blk_dm_multod_extract - in SPARSE SPARSE case, nzmax exceeded\n");
							exit(1);
						}
					}
				}
				dm_multod(bm,dmx,dval);
				free_double_matrix(dmx);
				break; }
			case MAT_SPARSE_DIAG: {
				mat_double *dmx = double_matrix(N,N,MAT_SPARSE_DIAG,N,blkm->basis);
				int pp = 0;
				double val;
				dmx->irow[0] = 1;
				for (j=0; j<N; j++) {
					dmx->irow[j+1] = dmx->irow[j];
					if (dm->irow[j+pos+1] - dm->irow[j+pos] == 1)
						val = dm->data[dm->irow[j+pos]-1];
					else
						val = 0;
					if (fabs(val)<SPARSE_TOL) continue;
					dmx->icol[pp] = j+1;
					dmx->data[pp] = val;
					dmx->irow[j+1]++;
					pp++;
				}
				dm_multod(bm,dmx,dval);
				free_double_matrix(dmx);
				break; }
			default:
				fprintf(stderr,"Error: blk_dm_multod_extract - unknown from-matrix type %d\n",dm->type);
				exit(1);
			}
			*********** original coding - END *****/
			break; }
		case MAT_SPARSE_DIAG: {
			mat_double *dmx = dm_extract_block(dm,pos,pos,N,N);
			if (dmx != NULL) {
				dm_multod(bm,dmx,dval);
				free_double_matrix(dmx);
			}
			break; }
		default:
			fprintf(stderr,"Error: blk_dm_multod_extract - unknown matrix %d type %d\n",i,bm->type);
			exit(1);
		}
		pos += N;
	}
}

void blk_cm_multod_extract(blk_mat_complx *blkm, mat_complx *dm, double dval)
{
	// z mat_double si vytahne prislusne diagonalni bloky a ty pricte do blkm

	assert(blkm->basis == dm->basis);

	int i, j, N, pos = 0;
	mat_complx *bm;
	complx z;

	for (i=0; i<blkm->Nblocks; i++) {
		N = blkm->blk_dims[i];
		bm = blkm->m + i;
		if (N == 1) {
			pos++;
			z = cm_getelem(dm,pos,pos);
			bm->data[0].re += dval*z.re;
			bm->data[0].im += dval*z.im;
			continue;
		}
		switch (bm->type) {
		case MAT_DENSE:
			switch (dm->type) {
			case MAT_DENSE:
				z = Complx(dval,0.0);
				for (j=0; j<N; j++) {
					cblas_zaxpy(N,&z,dm->data+pos+(j+pos)*dm->row,1,bm->data+j*N,1);
				}
				break;
			case MAT_DENSE_DIAG:
				z =Complx(dval,0.0);
				cblas_zaxpy(N,&z,dm->data+pos,1,bm->data,N+1);
				break;
			case MAT_SPARSE: {
				int r, c, nc, cs;
				for (r=0; r<N; r++) {
					cs = dm->irow[r+pos];
					nc = dm->irow[r+pos+1] - cs;
					for (j=0; j<nc; j++) {
						c = dm->icol[cs+j-1] - pos - 1;
						if (c>=0 && c<N) {
							complx *zz = bm->data + r+c*N;
							z = dm->data[cs+j-1];
							zz->re += dval*z.re;
							zz->im += dval*z.im;
						}
					}
				}
				break; }
			case MAT_SPARSE_DIAG: {
				int cs;
				for (j=0; j<N; j++) {
					cs = dm->irow[j+pos];
					if (dm->irow[j+pos+1] - cs == 1) {
						complx *zz = bm->data + j+j*N;
						z = dm->data[cs-1];
						zz->re += dval*z.re;
						zz->im += dval*z.im;
					}
				}
				break; }
			default:
				fprintf(stderr,"Error: blk_cm_multod_extract - unknown from-matrix type %d\n",dm->type);
				exit(1);
			}
			break;
		case MAT_DENSE_DIAG:
			fprintf(stderr,"Error: blk_cm_multod_extract - matrix %d type %d not implemented\n",i,bm->type);
			exit(1);
			break;
		case MAT_SPARSE:
			switch (dm->type) {
			case MAT_DENSE: {
				mat_complx *dmx = complx_matrix(N,N,MAT_DENSE,1,blkm->basis);
				for (j=0; j<N; j++) {
					cblas_zcopy(N,dmx->data+j*N,1,dm->data+pos+(j+pos)*dm->row,1);
				}
				cm_sparse(dmx,SPARSE_TOL);
				cm_multod(bm,dmx,dval);
				free_complx_matrix(dmx);
				break; }
			case MAT_DENSE_DIAG: {
				mat_complx *dmx = complx_matrix(N,N,MAT_SPARSE_DIAG,N,blkm->basis);
				for (j=0; j<N; j++) {
					dmx->irow[j] = j+1;
					dmx->icol[j] = j+1;
					dmx->data[j] = dm->data[j+pos];
				}
				dmx->irow[N] = N + 1;
				cm_multod(bm,dmx,dval);
				free_complx_matrix(dmx);
				break; }
			case MAT_SPARSE: {
				uint64_t nzmax = (uint64_t)((uint64_t)N*N*(1.0-SPARSITY));
				int r, c, nc, cs, pp = 0;
				complx val;
				mat_complx *dmx = complx_matrix(N,N,MAT_SPARSE,nzmax,blkm->basis);
				dmx->irow[0] = 1;
				for (r=0; r<N; r++) {
					cs = dm->irow[r+pos];
					nc = dm->irow[r+pos+1] - cs;
					dmx->irow[r+1] = dmx->irow[r];
					for (j=0; j<nc; j++) {
						c = dm->icol[cs-1] - pos - 1;
						if (c>=0 && c<N) val = dm->data[cs-1]; else val=Cnull;
						if (fabs(val.re)<SPARSE_TOL && fabs(val.im)<SPARSE_TOL) continue;
						dmx->irow[r+1]++;
						dmx->icol[pp] = c+1;
						dmx->data[pp] = val;
						pp++;
						if (pp>=nzmax) {
							fprintf(stderr,"Error: blk_cm_multod_extract - in SPARSE SPARSE case, nzmax exceeded\n");
							exit(1);
						}
					}
				}
				cm_multod(bm,dmx,dval);
				free_complx_matrix(dmx);
				break; }
			case MAT_SPARSE_DIAG: {
				mat_complx *dmx = complx_matrix(N,N,MAT_SPARSE_DIAG,N,blkm->basis);
				int pp = 0;
				complx val;
				dmx->irow[0] = 1;
				for (j=0; j<N; j++) {
					dmx->irow[j+1] = dmx->irow[j];
					if (dm->irow[j+pos+1] - dm->irow[j+pos] == 1)
						val = dm->data[dm->irow[j+pos]-1];
					else
						val = Cnull;
					if (fabs(val.re)<SPARSE_TOL && fabs(val.im)<SPARSE_TOL) continue;
					dmx->icol[pp] = j+1;
					dmx->data[pp] = val;
					dmx->irow[j+1]++;
					pp++;
				}
				cm_multod(bm,dmx,dval);
				free_complx_matrix(dmx);
				break; }
			default:
				fprintf(stderr,"Error: blk_cm_multod_extract - unknown from-matrix type %d\n",dm->type);
				exit(1);
			}
			break;
		case MAT_SPARSE_DIAG:
			fprintf(stderr,"Error: blk_cm_multod_extract - matrix %d type %d not implemented\n",i,bm->type);
			exit(1);
			break;
		default:
			fprintf(stderr,"Error: blk_cm_multod_extract - unknown matrix %d type %d\n",i,bm->type);
			exit(1);
		}
		pos += N;
	}
}

void blk_simtrans_zrot2(blk_mat_complx *U, mat_double *sumUph)
{
	assert(sumUph->type == MAT_DENSE_DIAG || sumUph->type == MAT_SPARSE_DIAG);
	assert(sumUph->basis == U->basis);

	int i, N, pos = 0;
	mat_complx *cm;
	complx expT[U->dim];

	complx *z = expT;
	double *dd = sumUph->data;
	if (sumUph->type == MAT_DENSE_DIAG) {
		for(i=0; i<U->dim; i++) {
			z->re = cos(*dd);
			z->im = sin(*dd);
			z++;
			dd++;
		}
	} else {
		/* sumUph->type == MAT_SPARSE_DIAG */
		/* strictly assumes only diagonal elements are present! */
		for (i=0; i<U->dim; i++) {
			int nc = sumUph->irow[i+1] - sumUph->irow[i];
			if (nc != 0) {
				z->re = cos(*dd);
				z->im = sin(*dd);
				dd++;
			} else {
				z->re = 1.0;
				z->im = 0.0;
			}
			z++;
		}
	}


	for (i=0; i<U->Nblocks; i++) {
		N = U->blk_dims[i];
		cm = U->m + i;
		if (N == 1) {
			// no change of cm
			pos++;
			continue;
		}
		switch (cm->type) {
		case MAT_DENSE: {
			int i2,j2;
			complx *zz, *dU, *dU2;
			double rr, ii, c, s;
			z = expT+pos;
			for (i2=0; i2<N; i2++) { /* columns */
				dU = dU2 = cm->data + (N+1)*i2; // sit on the diagonal
				zz = expT+pos + i2; //dP2 = T->data + i;
				for (j2=i2+1; j2<N; j2++) { /* rows, lower triangle */
					dU++;
					dU2 += N;
					zz++;
					c = z->re*zz->re + z->im*zz->im;
					s = z->im*zz->re - z->re*zz->im;
					rr = dU->re;
					ii = dU->im;
					dU->re = rr*c - ii*s;
					dU->im = ii*c + rr*s;
					rr = dU2->re;
					ii = dU2->im;
					dU2->re = rr*c + ii*s;
					dU2->im = ii*c - rr*s;
				}
				z++; //dP++;
			}
			break; }
		case MAT_DENSE_DIAG:
			// no change of cm
			break;
		case MAT_SPARSE: {
			int i2, j2, nc, *ic = cm->icol;
			double re, s, c;
			complx *z1, *z2;
			z = cm->data;
			for (i2=0; i2<N; i2++) {
				nc = cm->irow[i2+1] - cm->irow[i2];
				z1 = expT+pos+i2;
				for (j2=0; j2<nc; j2++) {
					z2 = expT+pos+(*ic)-1;
					c = z1->re*z2->re + z1->im*z2->im;
					s = z1->im*z2->re - z1->re*z2->im;
					re = z->re;
					z->re = re*c + z->im*s;
					z->im = z->im*c - re*s;
					z++;
					ic++;
				}
			}
			break; }
		case MAT_SPARSE_DIAG:
			// no change of cm
			break;
		default:
			fprintf(stderr,"Error: blk_simtrans_zrot2 - unknown matrix %d type %d\n",i,cm->type);
			exit(1);
		}
		pos += N;
	}

}

/* change structure of blkm according to structure of ham ( = hamiltonian) */
void blk_change_structure_double(blk_mat_double *blkm, blk_mat_double *ham)
{
	if ( (blkm->basis != ham->basis) || (blkm->Nblocks != ham->Nblocks) ) {
		const int NN = ham->Nblocks;
		int i;
		mat_double *mx;
		for (i=0; i<blkm->Nblocks; i++) {
			mx = blkm->m + i;
			if (mx->data != NULL) free(mx->data);
			if (mx->irow != NULL) free(mx->irow);
			if (mx->icol != NULL) free(mx->icol);
		}
		if (blkm->Nblocks != NN) {
			mx = (mat_double*)realloc(blkm->m, NN*sizeof(mat_double));
			if (mx == NULL) {
				fprintf(stderr, "Error: blk_change_structure_double - problem when reallocating m\n");
				exit(1);
			}
			blkm->m = mx;
			int *ddd = (int*)realloc(blkm->blk_dims, NN*sizeof(int));
			if (ddd == NULL) {
				fprintf(stderr, "Error: blk_change_structure_double - problem when reallocating blk_dims\n");
				exit(1);
			}
			blkm->blk_dims = ddd;
			blkm->Nblocks = NN;
		}
		blkm->basis = ham->basis;
		for (i=0; i<NN; i++) {
			mx = blkm->m+i;
			create_sqmat_double(mx,ham->blk_dims[i],ham->m[i].type,1,ham->basis);
			blkm->blk_dims[i] = ham->blk_dims[i];
		}
	}

}

void blk_change_structure_double_nondiag(blk_mat_double *blkm, blk_mat_double *ham, int issparse)
{
	if ( (blkm->basis != ham->basis) || (blkm->Nblocks != ham->Nblocks) ) {
		const int NN = ham->Nblocks;
		int i, mtype;
		mat_double *mx;
		for (i=0; i<blkm->Nblocks; i++) {
			mx = blkm->m + i;
			if (mx->data != NULL) free(mx->data);
			if (mx->irow != NULL) free(mx->irow);
			if (mx->icol != NULL) free(mx->icol);
		}
		if (blkm->Nblocks != NN) {
			mx = (mat_double*)realloc(blkm->m, NN*sizeof(mat_double));
			if (mx == NULL) {
				fprintf(stderr, "Error: blk_change_structure_double - problem when reallocating m\n");
				exit(1);
			}
			blkm->m = mx;
			int *ddd = (int*)realloc(blkm->blk_dims, NN*sizeof(int));
			if (ddd == NULL) {
				fprintf(stderr, "Error: blk_change_structure_double - problem when reallocating blk_dims\n");
				exit(1);
			}
			blkm->blk_dims = ddd;
			blkm->Nblocks = NN;
		}
		blkm->basis = ham->basis;
		for (i=0; i<NN; i++) {
			mx = blkm->m+i;
			if (issparse && ham->blk_dims[i]>MAXFULLDIM) {
				mtype = MAT_SPARSE;
			} else {
				mtype = MAT_DENSE;
			}
			create_sqmat_double(mx,ham->blk_dims[i],mtype,1,ham->basis);
			blkm->blk_dims[i] = ham->blk_dims[i];
		}
	}

}

void blk_change_structure_complx(blk_mat_complx *blkm, blk_mat_double *ham)
{
	if ( (blkm->basis != ham->basis) || (blkm->Nblocks != ham->Nblocks) ) {
		const int NN = ham->Nblocks;
		int i;
		mat_complx *mx;
		for (i=0; i<blkm->Nblocks; i++) {
			mx = blkm->m + i;
			if (mx->data != NULL) free(mx->data);
			if (mx->irow != NULL) free(mx->irow);
			if (mx->icol != NULL) free(mx->icol);
		}
		if (blkm->Nblocks != NN) {
			mx = (mat_complx*)realloc(blkm->m, NN*sizeof(mat_complx));
			if (mx == NULL) {
				fprintf(stderr, "Error: blk_change_structure_complx - problem when reallocating m\n");
				exit(1);
			}
			blkm->m = mx;
			int *ddd = (int*)realloc(blkm->blk_dims, NN*sizeof(int));
			if (ddd == NULL) {
				fprintf(stderr, "Error: blk_change_structure_complx - problem when reallocating blk_dims\n");
				exit(1);
			}
			blkm->blk_dims = ddd;
			blkm->Nblocks = NN;
		}
		blkm->basis = ham->basis;
		for (i=0; i<NN; i++) {
			mx = blkm->m+i;
			create_sqmat_complx(mx,ham->blk_dims[i],ham->m[i].type,1,ham->basis);
			blkm->blk_dims[i] = ham->blk_dims[i];
		}
	}

}


void blk_change_structure_complx2(blk_mat_complx *blkm, blk_mat_complx *cmx)
{
	if ( (blkm->basis != cmx->basis) || (blkm->Nblocks != cmx->Nblocks) ) {
		const int NN = cmx->Nblocks;
		int i;
		mat_complx *mx;
		for (i=0; i<blkm->Nblocks; i++) {
			mx = blkm->m + i;
			if (mx->data != NULL) free(mx->data);
			if (mx->irow != NULL) free(mx->irow);
			if (mx->icol != NULL) free(mx->icol);
		}
		if (blkm->Nblocks != NN) {
			mx = (mat_complx*)realloc(blkm->m, NN*sizeof(mat_complx));
			if (mx == NULL) {
				fprintf(stderr, "Error: blk_change_structure_complx2 - problem when reallocating m\n");
				exit(1);
			}
			blkm->m = mx;
			int *ddd = (int*)realloc(blkm->blk_dims, NN*sizeof(int));
			if (ddd == NULL) {
				fprintf(stderr, "Error: blk_change_structure_complx2 - problem when reallocating blk_dims\n");
				exit(1);
			}
			blkm->blk_dims = ddd;
			blkm->Nblocks = NN;
		}
		blkm->basis = cmx->basis;
		for (i=0; i<NN; i++) {
			mx = blkm->m+i;
			create_sqmat_complx(mx,cmx->blk_dims[i],cmx->m[i].type,1,cmx->basis);
			blkm->blk_dims[i] = cmx->blk_dims[i];
		}
	}
}


/* copy block matrix to full matrix prepared in advance */
void blk_dm_copy_2(mat_double *dest, blk_mat_double *blkm) {
	int dim = blkm->dim;

	assert(dest->row == dim && dest->col == dim);

	switch (dest->type ) {
	case MAT_DENSE: {
		int i, NN, pos = 0;
		mat_double *smx;
		dm_zero(dest);
		for (i=0; i<blkm->Nblocks; i++) {
			smx = blkm->m + i;
			NN = blkm->blk_dims[i];
			if (NN == 1) {
				dest->data[pos+pos*dim] = smx->data[0];
			} else {
				switch (smx->type) {
				case MAT_DENSE: {
					int j2;
					for (j2=0; j2<NN; j2++) {
						memcpy(dest->data+pos+(pos+j2)*dim,smx->data+j2*NN,NN*sizeof(double));
					}
					break; }
				case MAT_DENSE_DIAG:
					cblas_dcopy(NN,smx->data,1,dest->data+pos+pos*dim,dim+1);
					break;
				case MAT_SPARSE:
					fprintf(stderr,"...ups, this should not happen, unfinished sparse code\n");
				default:
					fprintf(stderr,"Error: blk_dm_copy_2 - unsupported matrix type (%d), place A\n",smx->type);
					exit(1);
				}
			}
			pos += NN;
		}
		break; }
	case MAT_DENSE_DIAG: {
		int i, NN, pos = 0;
		mat_double *mx;
		for (i=0; i<blkm->Nblocks; i++) {
			mx = blkm->m + i;
			NN = blkm->blk_dims[i];
			if (NN == 1) {
				dest->data[pos] = mx->data[0];
			} else {
				assert(mx->type == MAT_DENSE_DIAG);
				memcpy(dest->data+pos,mx->data,NN*sizeof(double));
			}
			pos += NN;
		}
		break; }
	case MAT_SPARSE: {
		// cil je ridka, ale zdroj muze byt husty blok s neznamym poctem nnz
		const uint64_t nnzmax = (uint64_t)floor((uint64_t)dim*dim*(1.0-SPARSITY));
		int nnz = blk_dm_nnz(blkm);
		int i, NN, pos = 0;
		if (nnz > nnzmax) {
			// cil bude husta
			//printf("blk_dm_copy_2 entered not tested code!!!\n");
			//blk_dm_print(blkm,"blk_dm_copy_2 initial");
			double *data = (double*)calloc(dim*dim,sizeof(double));
			for (i=0; i<blkm->Nblocks; i++) {
				mat_double *bmx = blkm->m+i;
				NN = blkm->blk_dims[i];
				if (NN == 1) {
					data[pos+pos*dim] = bmx->data[0];
				} else {
					switch (bmx->type) {
					case MAT_DENSE: {
						int i2;
						for (i2=0; i2<NN; i2++) {
							memcpy(data+pos+(pos+i2)*dim,bmx->data+i2*NN,NN*sizeof(double));
						}
						break; }
					case MAT_DENSE_DIAG:
						cblas_dcopy(NN, bmx->data,1,data+pos+pos*dim,dim);
						break;
					case MAT_SPARSE:
					case MAT_SPARSE_DIAG: {
						int r, c, *ic = bmx->icol;
						double *z = bmx->data;
						for (r=0; r<NN; r++) {
							for (c=bmx->irow[r]; c<bmx->irow[r+1]; c++) {
								data[pos+r+(pos+(*ic)-1)*dim] = *z;
								z++;
								ic++;
							}
						}
						break; }
					default:
						fprintf(stderr,"Error: blk_dm_copy_2 - unknown submatrix type (%d), place XUX\n",bmx->type);
						exit(1);
					}
				}
				pos += NN;
			}
			free(dest->irow);
			dest->irow = NULL;
			free(dest->icol);
			dest->icol = NULL;
			dest->type = MAT_DENSE;
			free(dest->data);
			dest->data = data;
			//dm_print(dest,"blk_dm_copy_2 result");
		} else {
			dm_change_nnz(dest,nnz);
			double *dd = dest->data;
			int *ic = dest->icol, *ir = dest->irow;
			ir[0] = 1;
			for (i=0; i<blkm->Nblocks; i++) {
				NN = blkm->blk_dims[i];
				mat_double *bmx = blkm->m+i;
				if (NN == 1) {
					*dd = bmx->data[0];
					*ic = pos+1;
					ir[pos+1] = ir[pos] + 1;
					dd++;
					ic++;
				} else {
					switch (bmx->type) {
					case MAT_DENSE: {
						int rr, cc;
						for (rr=0; rr<NN; rr++) {
							ir[pos+1+rr] = ir[pos+rr];
							for (cc=0; cc<NN; cc++) {
								double bdd = bmx->data[rr+cc*NN];
								if (fabs(bdd) > SPARSE_TOL) {
									*dd = bdd;
									*ic = pos+1+cc;
									ir[pos+1+rr]++;
									dd++;
									ic++;
								}
							}
						}
						break; }
					case MAT_DENSE_DIAG: {
						int rr;
						for (rr=0; rr<NN; rr++) {
							ir[pos+1+rr] = ir[pos+rr];
							double bdd = bmx->data[rr];
							if (fabs(bdd) > SPARSE_TOL) {
								*dd = bdd;
								*ic = pos+1+rr;
								ir[pos+1+rr]++;
								dd++;
								ic++;
							}
						}
						break; }
					case MAT_SPARSE:
					case MAT_SPARSE_DIAG: {
						int rr, cc, *icb=bmx->icol;
						double *bdd = bmx->data;
						for (rr=0; rr<NN; rr++) {
							ir[pos+1+rr] = ir[pos+rr];
							for (cc=bmx->irow[rr+1]; cc<bmx->irow[rr]; cc++) {
								if (fabs(*bdd) > SPARSE_TOL) {
									*dd = *bdd;
									*ic = pos + *icb;
									dd++;
									ic++;
									ir[pos+1+rr]++;
								}
								bdd++;
								icb++;
							}
						}
						break; }
					default:
						fprintf(stderr,"Errror: blk_dm_copy_2- unknown submatrix type (%d), place YHJ\n",bmx->type);
						exit(1);
					}
				}
				pos += NN;
			}
		}
		break; }
	default :
		fprintf(stderr,"Error: blk_dm_copy_2 - unsupported matrix type %d, place G\n",dest->type);
		exit(1);
	}

	dest->basis = blkm->basis;
}

void blk_cm_copy_2(mat_complx *dest, blk_mat_complx *blkm) {
	int dim = blkm->dim;

	assert(dest->row == dim && dest->col == dim);

	switch (dest->type ) {
	case MAT_DENSE: {
		int i, NN, pos = 0;
		mat_complx *smx;
		cm_zero(dest);
		for (i=0; i<blkm->Nblocks; i++) {
			smx = blkm->m + i;
			NN = blkm->blk_dims[i];
			if (NN == 1) {
				dest->data[pos+pos*dim] = smx->data[0];
			} else {
				switch (smx->type) {
				case MAT_DENSE: {
					int j2;
					for (j2=0; j2<NN; j2++) {
						memcpy(dest->data+pos+(pos+j2)*dim,smx->data+j2*NN,NN*sizeof(complx));
					}
					break; }
				case MAT_DENSE_DIAG:
					cblas_zcopy(NN,smx->data,1,dest->data+pos+pos*dim,dim+1);
					break;
				case MAT_SPARSE: {
					int i2, j2, *ic = smx->icol;
					complx *iz = smx->data;
					for (i2=0;i2<NN;i2++) {
						int nc = smx->irow[i2+1] - smx->irow[i2];
						for (j2=0; j2<nc; j2++) {
							dest->data[pos+i2 + (pos+*ic-1)*dim] = *iz;
							iz++;
							ic++;
						}
					}
					break; }
				default:
					fprintf(stderr,"Error: blk_cm_copy_2 - unsupported matrix type (%d), place A\n",smx->type);
					exit(1);
				}
			}
			pos += NN;
		}
		break; }
	case MAT_DENSE_DIAG: {
		int i, NN, pos = 0;
		mat_complx *mx;
		for (i=0; i<blkm->Nblocks; i++) {
			mx = blkm->m + i;
			NN = blkm->blk_dims[i];
			if (NN == 1) {
				dest->data[pos] = mx->data[0];
			} else {
				assert(mx->type == MAT_DENSE_DIAG);
				memcpy(dest->data+pos,mx->data,NN*sizeof(complx));
			}
			pos += NN;
		}
		break; }
	case MAT_SPARSE:
		fprintf(stderr,"Error: blk_cm_copy_2 - unfinished code dealing with sparse matrices, place B\n");
		exit(1);
		// cil je ridka, ale zdroj muze byt husty blok s neznamym poctem nnz
	default :
		fprintf(stderr,"Error: blk_cm_copy_2 - unsupported matrix type %d, place G\n",dest->type);
		exit(1);
	}

	dest->basis = blkm->basis;

}


/****** copy full matrix into predefined block matrix, skiping nondiagonal blocks
void blk_dm_copy_3(blk_mat_double *dest, mat_double *m)
{
	assert(dest->basis == m->basis);
	assert(dest->dim == m->row && dest->dim == m->col);

	int i, NN, pos = 0;
	mat_double *dm;
	for (i=0; i<dest->Nblocks; i++) {
		dm = dest->m + i;
		NN = dest->blk_dims[i];
		if (NN == 1) {
			dm->data[0] = dm_getelem(m, pos+1, pos+1);
		} else {
			// tady podobnej kod jako v blk_dm_multod_extract
			fprintf(stderr,"Error: unfinished code in blk_dm_copy_3\n");
			exit(1);
		}
		pos += NN;
	}
}
*********************/


void change_basis_pulse(Sim_info *sim, Sim_wsp *wsp, int basis)
{
	// zmeni bazi pro hamiltonian i pro rf, vcetne ham_blk a sumHrf, wsp->dU a tmpU (pokud != NULL)
	// basis = 0 means there is no channel useful for blocking, make one whole block
	int i;

	//assert(sim->Hassembly == 1);

	// change basis of wsp->Hiso, this serves also as structure template
	if (sim->averaging == 0) {
		if (basis == 0) {
			// make one whole block
			if (sim->Hiso->basis == 0) {
				assert(sim->Hiso->Nblocks == 1); // also, sim->Hiso->m->type == MAT_DENSE_DIAG
				if (wsp->Hiso != sim->Hiso) {
					DEBUGPRINT("change_basis_pulse: no blocking and diag Ham\n");
					free_blk_mat_double(wsp->Hiso);
					wsp->Hiso = sim->Hiso;
				}
			} else {
				assert(sim->Hiso->basis == sim->Nbasis-1);
				assert(sim->Hint_isdiag != 1);
				DEBUGPRINT("change_basis_pulse: no blocking and stay in basis %d\n",sim->Hiso->basis);
				// no need to change basis back to 0, just make it one block
				mat_double *wmx;
				if (wsp->Hiso != sim->Hiso) {
					if (wsp->Hiso->Nblocks != 1) {
						// shrink to one block
						for (i=0; i<wsp->Hiso->Nblocks; i++) {
							wmx = wsp->Hiso->m + i;
							if (wmx->data != NULL) free(wmx->data);
							if (wmx->irow != NULL) free(wmx->irow);
							if (wmx->icol != NULL) free(wmx->icol);
						}
						wmx = (mat_double*)realloc(wsp->Hiso->m, sizeof(mat_double));
						if (wmx == NULL) {
							fprintf(stderr, "Error: change_basis_pulse - problem when reallocating m, place 1\n");
							exit(1);
						}
						wsp->Hiso->m = wmx;
						int *ddd = (int*)realloc(wsp->Hiso->blk_dims, sizeof(int));
						if (ddd == NULL) {
							fprintf(stderr, "Error: change_basis_pulse - problem when reallocating blk_dims, place 1\n");
							exit(1);
						}
						wsp->Hiso->blk_dims = ddd;
						wsp->Hiso->Nblocks = 1;
						create_sqmat_double(wsp->Hiso->m, sim->matdim, sim->sparse ? MAT_SPARSE : MAT_DENSE, 1,sim->Hiso->basis);
					}
				} else {
					// create wsp->Hiso
					wsp->Hiso = create_blk_mat_double(sim->matdim,1,NULL,sim->sparse ? MAT_SPARSE : MAT_DENSE,sim->Hiso->basis);
				}
				// fill wsp->Hiso
				blk_dm_copy_2(wsp->Hiso->m,sim->Hiso);
				//wsp->Hiso->basis = sim->Nbasis - 1;
			}
		} else {
			DEBUGPRINT("change_basis_pulse: blocking in basis %d and change there from basis %d\n",basis,sim->Hiso->basis);
			if (wsp->Hiso->basis != basis) {
				if (wsp->Hiso != sim->Hiso) free_blk_mat_double(wsp->Hiso);
				if (sim->Hint_isdiag) {
					wsp->Hiso = create_blk_mat_double(sim->matdim,LEN(sim->dims_table[basis]),sim->dims_table[basis],MAT_DENSE_DIAG, basis);
				} else {
					wsp->Hiso = create_blk_mat_double(sim->matdim,LEN(sim->dims_table[basis]),sim->dims_table[basis],sim->sparse ? MAT_SPARSE : MAT_DENSE, basis);
				}
				blk_dm_change_basis(wsp->Hiso,sim->Hiso,sim);
			}
		}
		// change basis of wsp->HQ[]
		if (wsp->HQ[0] != NULL) {
			int j;
			if (basis == 0) {
				if (sim->Hiso->basis == 0) {
					for (j=0; j<5; j++) {
						if (wsp->HQ[j] != sim->HQ[j]) free_blk_mat_double(wsp->HQ[j]);
						wsp->HQ[j] = sim->HQ[j];
					}
				} else {
					for (j=0; j<5; j++) {
						if (wsp->HQ[j]->basis == sim->Hiso->basis && wsp->HQ[j]->Nblocks == 1) continue;
						if (wsp->HQ[j] != sim->HQ[j]) free_blk_mat_double(wsp->HQ[j]);
						wsp->HQ[j] = create_blk_mat_double(sim->matdim,1,NULL,sim->sparse ? MAT_SPARSE : MAT_DENSE,sim->HQ[j]->basis);
						blk_dm_copy_2(wsp->HQ[j]->m,sim->HQ[j]);
					}
				}
			} else {
				for (j=0; j<5; j++) {
					if (wsp->HQ[j]->basis == basis) continue;
					if (wsp->HQ[j] != sim->HQ[j]) free_blk_mat_double(wsp->HQ[j]);
					wsp->HQ[j] = create_blk_mat_double_copy(wsp->Hiso);
					blk_dm_change_basis(wsp->HQ[j],sim->HQ[j],sim);
				}
			}
		}
	} else {
		assert(sim->averaging == 1);
		blk_mat_double *dum;
		if (basis == 0) {
			// basis remains that defined in sim->Hiso, but we need single block
			if (sim->Hiso->basis == 0) {
				// hamiltonian is diagonal
				if (wsp->Hiso->basis != 0 || wsp->Hiso->Nblocks != 1) {
					dum = create_blk_mat_double(sim->matdim,1,NULL,MAT_DENSE_DIAG, 0);
					blk_dm_change_basis(dum,wsp->Hiso,sim);
					free_blk_mat_double(wsp->Hiso);
					wsp->Hiso = dum;
				}
			} else {
				if (wsp->Hiso->basis != sim->Nbasis-1 || wsp->Hiso->Nblocks != 1) {
					dum = create_blk_mat_double(sim->matdim,1,NULL,sim->sparse ? MAT_SPARSE : MAT_DENSE, sim->Nbasis-1);
					blk_dm_change_basis(dum,wsp->Hiso,sim);
					free_blk_mat_double(wsp->Hiso);
					wsp->Hiso = dum;
				}
			}
		} else {
			if (wsp->Hiso->basis != basis) {
				if (sim->Hint_isdiag) {
					dum = create_blk_mat_double(sim->matdim,LEN(sim->dims_table[basis]),sim->dims_table[basis],MAT_DENSE_DIAG, basis);
				} else {
					dum = create_blk_mat_double(sim->matdim,LEN(sim->dims_table[basis]),sim->dims_table[basis],sim->sparse ? MAT_SPARSE : MAT_DENSE, basis);
				}
				blk_dm_change_basis(dum,wsp->Hiso,sim);
				free_blk_mat_double(wsp->Hiso);
				wsp->Hiso = dum;
			}
		}
		// change basis of wsp->HQ[]
		if (wsp->HQ[0] != NULL) {
			int j;
			for (j=0; j<5; j++) {
				if (wsp->HQ[j]->basis != wsp->Hiso->basis || wsp->HQ[j]->Nblocks != wsp->Hiso->Nblocks) {
					dum = create_blk_mat_double_copy(wsp->Hiso);
					blk_dm_change_basis(dum, wsp->HQ[j], sim);
					free_blk_mat_double(wsp->HQ[j]);
					wsp->HQ[j] = dum;
				}
			}
		}
	}


	// change basis of wsp->Hiso_off
	if (wsp->Hiso_off != NULL) {
		if ( (wsp->Hiso_off->basis != wsp->Hiso->basis) || (wsp->Hiso_off->Nblocks != wsp->Hiso->Nblocks) ) {
			blk_mat_double *dum = create_blk_mat_double_copy(wsp->Hiso);
			blk_dm_change_basis(dum,wsp->Hiso_off,sim);
			free_blk_mat_double(wsp->Hiso_off);
			wsp->Hiso_off = dum;
		}
	}

	// change basis of wsp->HQ_off[]
	if (wsp->HQ_off[0] != NULL) {
		if ( (wsp->HQ_off[0]->basis != wsp->Hiso->basis) || (wsp->HQ_off[0]->Nblocks != wsp->Hiso->Nblocks) ) {
			int i;
			for (i=0; i<5; i++) {
				blk_mat_double *dum = create_blk_mat_double_copy(wsp->Hiso);
				blk_dm_change_basis(dum,wsp->HQ_off[i],sim);
				free_blk_mat_double(wsp->HQ_off[i]);
				wsp->HQ_off[i] = dum;
			}
		}
	}

	// change basis for 2nd and 3rd order quadrupoles
	if (sim->nQ != 0) {
		for (i=0; i<sim->nQ; i++) {
			if (sim->Q[i]->order == 1) continue;
			assert( (wsp->QTa[i] != NULL) && (wsp->QTb[i] != NULL) );
			if (wsp->QTa[i]->basis != wsp->Hiso->basis) {
				if (wsp->QTa[i] != sim->Q[i]->Ta) free_double_matrix(wsp->QTa[i]);
				wsp->QTa[i] = dm_change_basis_2(sim->Q[i]->Ta,wsp->Hiso->basis,sim);
				if (wsp->QTb[i] != sim->Q[i]->Tb) free_double_matrix(wsp->QTb[i]);
				wsp->QTb[i] = dm_change_basis_2(sim->Q[i]->Tb,wsp->Hiso->basis,sim);
			}
			if (sim->Q[i]->order == 3) {
				assert( (wsp->QT3a[i] != NULL) && (wsp->QT3b[i] != NULL) && (wsp->QT3c[i] != NULL) );
				if (wsp->QT3a[i]->basis != wsp->Hiso->basis) {
					if (wsp->QT3a[i] != sim->Q[i]->T3a) free_double_matrix(wsp->QT3a[i]);
					wsp->QT3a[i] = dm_change_basis_2(sim->Q[i]->T3a,wsp->Hiso->basis,sim);
					if (wsp->QT3b[i] != sim->Q[i]->T3b) free_double_matrix(wsp->QT3b[i]);
					wsp->QT3b[i] = dm_change_basis_2(sim->Q[i]->T3b,wsp->Hiso->basis,sim);
					if (wsp->QT3c[i] != sim->Q[i]->T3c) free_double_matrix(wsp->QT3c[i]);
					wsp->QT3c[i] = dm_change_basis_2(sim->Q[i]->T3c,wsp->Hiso->basis,sim);
				}
			}
		}
	}

	// change basis for mixing terms
	if (sim->nMIX != 0) {
		for (i=0; i<sim->nMIX; i++) {
			assert(wsp->MT[i] != NULL);
			if (wsp->MT[i]->basis != wsp->Hiso->basis) {
				if (wsp->MT[i] != sim->MIX[i]->T) free_double_matrix(wsp->MT[i]);
				wsp->MT[i] = dm_change_basis_2(sim->MIX[i]->T,wsp->Hiso->basis,sim);
			}
			if (wsp->MTa[i] != NULL) {
				if ((wsp->MTa[i]->basis != wsp->Hiso->basis) || (wsp->MTa[i]->Nblocks != wsp->Hiso->Nblocks)) {
					blk_mat_double *dum = create_blk_mat_double_copy(wsp->Hiso);
					blk_dm_change_basis(dum,wsp->MTa[i],sim);
					if (wsp->MTa[i] != sim->MIX[i]->Ta) free_blk_mat_double(wsp->MTa[i]);
					wsp->MTa[i] = dum;
				}
				if ((wsp->MTb[i]->basis != wsp->Hiso->basis) || (wsp->MTb[i]->Nblocks != wsp->Hiso->Nblocks)) {
					blk_mat_double *dum = create_blk_mat_double_copy(wsp->Hiso);
					blk_dm_change_basis(dum,wsp->MTb[i],sim);
					if (wsp->MTb[i] != sim->MIX[i]->Tb) free_blk_mat_double(wsp->MTb[i]);
					wsp->MTb[i] = dum;
				}
			}
		}
	}

	// change basis for rf pulse channels
	for (i=1; i<=sim->ss->nchan; i++) {
		if (wsp->chan_Ix[i]->basis == wsp->Hiso->basis) continue;
		mat_double *dum = dm_change_basis_2(wsp->chan_Ix[i],wsp->Hiso->basis,sim);
		free_double_matrix(wsp->chan_Ix[i]);
		wsp->chan_Ix[i] = dum;
		assert(wsp->chan_Iz[i]->basis != wsp->Hiso->basis);
		dum = dm_change_basis_2(wsp->chan_Iz[i],wsp->Hiso->basis,sim);
		free_double_matrix(wsp->chan_Iz[i]);
		wsp->chan_Iz[i] = dum;
		// chan_Iy is not needed for pulses, only for OC gradients
	}
	wsp->sumUph->basis = wsp->Hiso->basis;
	if (wsp->isselectivepulse) {
		for (i=1; i<=LEN(wsp->spinused); i++) {
			if (wsp->spinused[i] == 0) continue;
			assert(wsp->Ix[i] != NULL);
			if (wsp->Ix[i]->basis == wsp->Hiso->basis) continue;
			mat_double *dum = dm_change_basis_2(wsp->Ix[i],wsp->Hiso->basis,sim);
			free_double_matrix(wsp->Ix[i]);
			wsp->Ix[i] = dum;
			assert(wsp->Iz[i]->basis != wsp->Hiso->basis);
			dum = dm_change_basis_2(wsp->Iz[i],wsp->Hiso->basis,sim);
			free_double_matrix(wsp->Iz[i]);
			wsp->Iz[i] = dum;
		}
	}

	// change structure of ham_blk, sumHrf, dU, tmpU
	blk_change_structure_double(wsp->ham_blk,wsp->Hiso);
	if (wsp->sumHrf == NULL) {
		wsp->sumHrf = create_blk_mat_double(sim->matdim,wsp->Hiso->Nblocks,sim->dims_table[wsp->Hiso->basis],sim->sparse ? MAT_SPARSE : MAT_DENSE,wsp->Hiso->basis);
	} else {
		blk_change_structure_double_nondiag(wsp->sumHrf,wsp->Hiso,sim->sparse);
	}
	if (wsp->dU != NULL) {
		blk_change_structure_complx(wsp->dU,wsp->Hiso);
	}
	if (wsp->tmpU != NULL) {
		blk_change_structure_complx(wsp->tmpU,wsp->Hiso);
	}

}

void change_basis_delay(Sim_info *sim, Sim_wsp *wsp)
{
	// basis change before 'delay' period
	// zmeni bazi pro hamiltonian pri delay, vcetne ham_blk, wsp->dU a tmpU (pokud != NULL)
	// (ucel je rusit zmenu baze pred pulsama)
	int i;

	//assert(sim->Hassembly == 1);

	// point wsp->Hiso, HQ, QTa, QTb back to sim
	// wsp->Hiso serves as structure template
	if (sim->averaging == 0) {
		if (wsp->Hiso != sim->Hiso) {
			free_blk_mat_double(wsp->Hiso);
			wsp->Hiso = sim->Hiso;
			DEBUGPRINT("change_basis_delay: wsp->Hiso reset\n");
			if (wsp->HQ[0] != NULL) {
				for (i=0; i<5; i++) {
					assert(wsp->HQ[i] != sim->HQ[i]);
					free_blk_mat_double(wsp->HQ[i]);
					wsp->HQ[i] = sim->HQ[i];
				}
			}
			if (sim->nQ != 0) {
				for (i=0; i<sim->nQ; i++) {
					if (sim->Q[i]->order == 1) continue;
					assert(wsp->QTa[i] != sim->Q[i]->Ta && wsp->QTa[i] != NULL);
					free_double_matrix(wsp->QTa[i]);
					wsp->QTa[i] = sim->Q[i]->Ta;
					assert(wsp->QTb[i] != sim->Q[i]->Tb && wsp->QTb[i] != NULL);
					free_double_matrix(wsp->QTb[i]);
					wsp->QTb[i] = sim->Q[i]->Tb;
					if (sim->Q[i]->order == 3) {
						assert(wsp->QT3a[i] != sim->Q[i]->T3a && wsp->QT3a[i] != NULL);
						free_double_matrix(wsp->QT3a[i]);
						wsp->QT3a[i] = sim->Q[i]->T3a;
						assert(wsp->QT3b[i] != sim->Q[i]->T3b && wsp->QT3b[i] != NULL);
						free_double_matrix(wsp->QT3b[i]);
						wsp->QT3b[i] = sim->Q[i]->T3b;
						assert(wsp->QT3c[i] != sim->Q[i]->T3c && wsp->QT3c[i] != NULL);
						free_double_matrix(wsp->QT3c[i]);
						wsp->QT3c[i] = sim->Q[i]->T3c;
					}
				}
			}
			if (sim->nMIX != 0) {
				for (i=0; i<sim->nMIX; i++) {
					assert(wsp->MT[i] != sim->MIX[i]->T);
					free_double_matrix(wsp->MT[i]);
					wsp->MT[i] = sim->MIX[i]->T;
					if (wsp->MTa[i] != NULL) {
						assert(wsp->MTa[i] != sim->MIX[i]->Ta);
						free_blk_mat_double(wsp->MTa[i]);
						wsp->MTa[i] = sim->MIX[i]->Ta;
						assert(wsp->MTb[i] != sim->MIX[i]->Tb);
						free_blk_mat_double(wsp->MTb[i]);
						wsp->MTb[i] = sim->MIX[i]->Tb;
					}
				}
			}
		}
	} else {
		assert(sim->averaging == 1);
		blk_mat_double *blkm;
		if (wsp->Hiso->basis != sim->Hiso->basis) {
			blkm = create_blk_mat_double_copy(sim->Hiso);
			blk_dm_change_basis(blkm, wsp->Hiso,sim);
			free_blk_mat_double(wsp->Hiso);
			wsp->Hiso = blkm;
		}
		if (sim->Hassembly) {
			for (i=0; i<5; i++) {
				assert(wsp->HQ[i] != sim->HQ[i]);
				if (wsp->HQ[i]->basis != wsp->Hiso->basis) {
					blkm = create_blk_mat_double_copy(wsp->Hiso);
					blk_dm_change_basis(blkm, wsp->HQ[i],sim);
					free_blk_mat_double(wsp->HQ[i]);
					wsp->HQ[i] = blkm;
				}
			}
		}
		for (i=0; i<sim->nQ; i++) {
			if (sim->Q[i]->order == 1) continue;
			if (wsp->QTa[i] != sim->Q[i]->Ta) {
				free_double_matrix(wsp->QTa[i]);
				wsp->QTa[i] = sim->Q[i]->Ta;
			}
			if (wsp->QTb[i] != sim->Q[i]->Tb) {
				free_double_matrix(wsp->QTb[i]);
				wsp->QTb[i] = sim->Q[i]->Tb;
			}
			if (sim->Q[i]->order == 3) {
				if (wsp->QT3a[i] != sim->Q[i]->T3a) {
					free_double_matrix(wsp->QT3a[i]);
					wsp->QT3a[i] = sim->Q[i]->T3a;
				}
				if (wsp->QT3b[i] != sim->Q[i]->T3b) {
					free_double_matrix(wsp->QT3b[i]);
					wsp->QT3b[i] = sim->Q[i]->T3b;
				}
				if (wsp->QT3c[i] != sim->Q[i]->T3c) {
					free_double_matrix(wsp->QT3c[i]);
					wsp->QT3c[i] = sim->Q[i]->T3c;
				}
			}
		}
		for (i=0; i<sim->nMIX; i++) {
			if (wsp->MT[i] != sim->MIX[i]->T) {
				free_double_matrix(wsp->MT[i]);
				wsp->MT[i] = sim->MIX[i]->T;
			}
			if (wsp->MTa[i] != sim->MIX[i]->Ta) {
				free_blk_mat_double(wsp->MTa[i]);
				wsp->MTa[i] = sim->MIX[i]->Ta;
			}
			if (wsp->MTb[i] != sim->MIX[i]->Tb) {
				free_blk_mat_double(wsp->MTb[i]);
				wsp->MTb[i] = sim->MIX[i]->Tb;
			}
		}
	}

	// change basis of Hiso_off and HQ_off
	if (wsp->Hiso_off != NULL) {
		if (wsp->Hiso_off->basis != wsp->Hiso->basis) {
			blk_mat_double *blkm = create_blk_mat_double_copy(wsp->Hiso);
			blk_dm_change_basis(blkm, wsp->Hiso_off,sim);
			free_blk_mat_double(wsp->Hiso_off);
			wsp->Hiso_off = blkm;
		}
	}
	if (wsp->HQ_off[0] != NULL) {
		if (wsp->HQ_off[0]->basis != wsp->Hiso->basis) {
			for (i=0; i<5; i++) {
				assert(wsp->HQ_off[i] != NULL);
				assert(wsp->HQ_off[i]->basis != wsp->Hiso->basis);
				blk_mat_double *blkm = create_blk_mat_double_copy(wsp->Hiso);
				blk_dm_change_basis(blkm, wsp->HQ_off[i],sim);
				free_blk_mat_double(wsp->HQ_off[i]);
				wsp->HQ_off[i] = blkm;
			}
		}
	}
	// change structure of ham_blk, dU, tmpU
	blk_change_structure_double(wsp->ham_blk,wsp->Hiso);
	if (wsp->dU != NULL) {
		blk_change_structure_complx(wsp->dU,wsp->Hiso);
	}
	if (wsp->tmpU != NULL) {
		blk_change_structure_complx(wsp->tmpU,wsp->Hiso);
	}
	DEBUGPRINT("change_basis_delay: done.\n");

}

int blk_cm_isdiag(blk_mat_complx *blkm)
{
	int i;

	for (i=0; i<blkm->Nblocks; i++) {
		if (blkm->blk_dims[i] > 1 ) {
			if (blkm->m[i].type == MAT_DENSE || blkm->m[i].type == MAT_SPARSE) break;
		}
	}

	return ( (i==blkm->Nblocks) ? 1 : 0 );
}

/*****
 *  diagonalize block matrix; on output Ud is replaced with its diagonal form
 * and transformation block matrix is returned separately
 *****/
blk_mat_complx * blk_cm_diag(blk_mat_complx *Ud)
{
	blk_mat_complx *T;
	mat_complx *bm, *tm;
	int i, NN, lwsp = 0;

	// determine maximal workspace needed for zgeev
	for (i=0; i<Ud->Nblocks; i++) {
		NN = Ud->blk_dims[i];
		if (NN > lwsp) lwsp = NN;
	}
	lwsp *= 3;
	int info=0;
	complx *wsp = (complx*)malloc(lwsp*sizeof(complx));
	double *rwork = (double*)malloc(lwsp*sizeof(double));
	const int ione = 1;

	mat_complx *dm = complx_matrix(Ud->dim,Ud->dim,MAT_DENSE_DIAG,0,Ud->basis);
	complx *dmpos = dm->data;
	T = (blk_mat_complx*)malloc(sizeof(blk_mat_complx));
	T->Nblocks = Ud->Nblocks;
	T->dim = Ud->dim;
	T->basis = Ud->basis;
	T->blk_dims = (int*)malloc(T->Nblocks*sizeof(int));
	T->m = (mat_complx*)malloc(T->Nblocks*sizeof(mat_complx));
	if (T->blk_dims == NULL || T->m == NULL) {
		fprintf(stderr,"Error: blk_cm_diag can not allocate T->m/blk_dims\n");
		exit(1);
	}

	for (i=0; i<Ud->Nblocks; i++) {
		NN = T->blk_dims[i] = Ud->blk_dims[i];
		bm = Ud->m + i;
		tm = T->m + i;
		if (NN == 1) {
			*dmpos = bm->data[0];
			create_sqmat_complx(tm,1,MAT_DENSE_DIAG,0,Ud->basis);
			tm->data[0] = Cunit;
		} else {
			switch (bm->type) {
			case MAT_DENSE:
				create_sqmat_complx(tm,NN,MAT_DENSE,0,Ud->basis);
				zgeev_("N","V",&NN,bm->data,&NN,dmpos,NULL,&ione,tm->data,&NN,wsp,&lwsp,rwork,&info);
				if ( info != 0) {
					fprintf(stderr,"blk_cm_diag error: diagonalization failed for submatrix %d with the code '%d'\n", i, info);
				    exit(1);
				}
				break;
			case MAT_DENSE_DIAG:
				create_sqmat_complx(tm,NN,MAT_DENSE_DIAG,0,Ud->basis);
				cm_unit(tm);
				memcpy(dmpos,bm->data,NN*sizeof(complx));
				break;
			case MAT_SPARSE:
				printf("WARNING!!! In order to diagonalize the submatrix %d (dim = %d) I need to \n"
						"convert it from sparse to full format. The whole operation may take long time...\n",i,NN);
				cm_dense(bm);
				create_sqmat_complx(tm,NN,MAT_DENSE,0,Ud->basis);
				zgeev_("N","V",&NN,bm->data,&NN,dmpos,NULL,&ione,tm->data,&NN,wsp,&lwsp,rwork,&info);
				if ( info != 0) {
					fprintf(stderr,"blk_cm_diag error: diagonalization failed for submatrix %d with the code '%d'\n", i, info);
				    exit(1);
				}
				free(bm->irow);
				free(bm->icol);
				break;
			case MAT_SPARSE_DIAG:
				fprintf(stderr,"Error: blk_cm_diag - block %d of MAT_SPASE_DIAG not implemented\n",i);
				exit(1);
			default:
				fprintf(stderr,"Error: blk_cm_diag - unknown submatrix %d type %d\n",i,bm->type);
				exit(1);
			}
		}
		dmpos += NN;

		free(bm->data);
	}

	free(wsp);
	free(rwork);
	free(Ud->m);
	Ud->m = dm;
	Ud->Nblocks = 1;
	int *ddd = realloc(Ud->blk_dims,sizeof(int));
	if (ddd == NULL) {
		fprintf(stderr,"Error: blk_cm_diag - can not reallocate Ud->blk_dims\n");
		exit(1);
	}
	Ud->blk_dims = ddd;
	Ud->blk_dims[0] = Ud->dim;

	return T;
}

/*****
 *  diagonalize block matrix; on output Ud is replaced with its diagonal form
 * and transformation block matrix is returned separately
 * ONLY reasonably small blocks will be diagonalized
 *****/
blk_mat_complx * blk_cm_diag_partly(blk_mat_complx *Ud)
{
	blk_mat_complx *T;
	mat_complx *bm, *tm;
	int i, NN, lwsp = 0;

	// determine maximal workspace needed for zgeev
	for (i=0; i<Ud->Nblocks; i++) {
		NN = Ud->blk_dims[i];
		if (NN > lwsp) lwsp = NN;
	}
	if (lwsp > MAXDIMDIAGONALIZE) lwsp = MAXDIMDIAGONALIZE;
	lwsp *= 3;
	int info=0;
	complx *wsp = (complx*)malloc(lwsp*sizeof(complx));
	double *rwork = (double*)malloc(lwsp*sizeof(double));
	if (wsp == NULL || rwork == NULL) {
		fprintf(stderr,"Error: blk_cm_diag_partly can not allocate workspace\n");
		exit(1);
	}
	const int ione = 1;

	T = (blk_mat_complx*)malloc(sizeof(blk_mat_complx));
	if (T == NULL) {
		fprintf(stderr,"Error: blk_cm_diag_partly can not allocate T\n");
		exit(1);
	}
	T->Nblocks = Ud->Nblocks;
	T->dim = Ud->dim;
	T->basis = Ud->basis;
	T->blk_dims = (int*)malloc(T->Nblocks*sizeof(int));
	T->m = (mat_complx*)malloc(T->Nblocks*sizeof(mat_complx));
	if (T->blk_dims == NULL || T->m == NULL) {
		fprintf(stderr,"Error: blk_cm_diag_partly can not allocate T->m/blk_dims\n");
		exit(1);
	}

	for (i=0; i<Ud->Nblocks; i++) {
		NN = T->blk_dims[i] = Ud->blk_dims[i];
		bm = Ud->m + i;
		tm = T->m + i;
		if (NN == 1) {
			create_sqmat_complx(tm,1,MAT_DENSE_DIAG,0,Ud->basis);
			tm->data[0] = Cunit;
		} else if (NN > MAXDIMDIAGONALIZE) {
			create_sqmat_complx(tm,NN,MAT_DENSE_DIAG,0,Ud->basis);
			cm_unit(tm);
		} else {
			//printf("diagonalizing block %d of dim %d\n",i+1,NN);
			switch (bm->type) {
			case MAT_DENSE: {
				create_sqmat_complx(tm,NN,MAT_DENSE,0,Ud->basis);
				mat_complx *ud = complx_matrix(NN,NN,MAT_DENSE_DIAG,0,Ud->basis);
				zgeev_("N","V",&NN,bm->data,&NN,ud->data,NULL,&ione,tm->data,&NN,wsp,&lwsp,rwork,&info);
				if ( info != 0) {
					fprintf(stderr,"blk_cm_diag error: diagonalization failed for submatrix %d with the code '%d'\n", i, info);
				    exit(1);
				}
				cm_swap_innards_and_destroy(bm,ud);
				break; }
			case MAT_DENSE_DIAG:
				create_sqmat_complx(tm,NN,MAT_DENSE_DIAG,0,Ud->basis);
				cm_unit(tm);
				break;
			case MAT_SPARSE: {
				printf("WARNING!!! In order to diagonalize the submatrix %d (dim = %d) I need to \n"
						"convert it from sparse to full format. The whole operation may take long time...\n",i,NN);
				cm_dense(bm);
				create_sqmat_complx(tm,NN,MAT_DENSE,0,Ud->basis);
				mat_complx *ud = complx_matrix(NN,NN,MAT_DENSE_DIAG,0,Ud->basis);
				zgeev_("N","V",&NN,bm->data,&NN,ud->data,NULL,&ione,tm->data,&NN,wsp,&lwsp,rwork,&info);
				if ( info != 0) {
					fprintf(stderr,"blk_cm_diag error: diagonalization failed for submatrix %d with the code '%d'\n", i, info);
				    exit(1);
				}
				cm_swap_innards_and_destroy(bm,ud);
				break; }
			case MAT_SPARSE_DIAG:
				cm_dense(bm);
				create_sqmat_complx(tm,NN,MAT_DENSE_DIAG,0,Ud->basis);
				cm_unit(tm);
			default:
				fprintf(stderr,"Error: blk_cm_diag - unknown submatrix %d type %d\n",i,bm->type);
				exit(1);
			}
		}
	}
	free(wsp);
	free(rwork);

	return T;
}

/*****
 * sorts N complex numbers according to their phase. No change done on eigs,
 * just coefficient map is filled.
 */
void eigsort(complx *eig, double *ph, int N, int *map)
{
	int i, j, k, l, ir, indxt, itemp, istack[64], jstack=0;
	double a;

	for (i=0; i<N; i++) {
		map[i] = i;
		ph[i] = Carg(eig[i]);
		//printf("%d: %g\n",i,ph[i]);
	}

	l = 0;
	ir = N-1;
	for (;;) {
		if (ir-l < 10) {
			/* insertion sort when subarray small */
			for (j=l+1; j<=ir; j++) {
				indxt = map[j];
				a = ph[indxt];
				for (i=j-1; i>=l; i--) {
					if (ph[map[i]] <= a) break;
					map[i+1] = map[i];
				}
				map[i+1] = indxt;
			}
			if (jstack == 0) break;
			ir = istack[jstack--];
			l = istack[jstack--];
		} else {
			k = (l +ir) >> 1;
			itemp = map[k]; map[k] = map[l+1]; map[l+1] = itemp;
			if (ph[map[l]] > ph[map[ir]]) {
				itemp = map[l]; map[l] = map[ir]; map[ir] = itemp;
			}
			if (ph[map[l+1]] > ph[map[ir]]) {
				itemp = map[l+1]; map[l+1] = map[ir]; map[ir] = itemp;
			}
			if (ph[map[l]] > ph[map[l+1]]) {
				itemp = map[l]; map[l] = map[l+1]; map[l+1] = itemp;
			}
			i = l+1;
			j = ir;
			indxt = map[l+1];
			a = ph[indxt];
			for (;;) {
				do i++; while (ph[map[i]] < a);
				do j--; while (ph[map[j]] > a);
				if (j < i) break;
				itemp = map[i]; map[i] = map[j]; map[j] = itemp;
			}
			map[l+1] = map[j];
			map[j] = indxt;
			jstack += 2;
			if (jstack > 64) {
				fprintf(stderr,"Error in eigsort: istack too small\n");
				exit(1);
			}
			if (ir-i+1 >= j-l) {
				istack[jstack] = ir;
				istack[jstack-1] = i;
				ir = j-1;
			} else {
				istack[jstack] = j-1;
				istack[jstack-1] = l;
				l = i;
			}
		}
	}
	//for (i=0; i<N; i++) {
	//	printf("%d: %g\n",i,ph[map[i]]);
	//}
}

/*****
 *  diagonalize block matrix; on output Ud is replaced with its diagonal form
 * and transformation block matrix is returned separately
 * ADDITIONALLY, eigenvalues are sorted according to their phase
 *              (does not check for unit eigenvalue amplitudes)
 *****/
blk_mat_complx * blk_cm_diag_sort(blk_mat_complx *Ud)
{
	blk_mat_complx *T;
	mat_complx *bm, *tm;
	int i, j, NN, lwsp = 0;

	// determine maximal workspace needed for zgeev
	for (i=0; i<Ud->Nblocks; i++) {
		NN = Ud->blk_dims[i];
		if (NN > lwsp) lwsp = NN;
	}
	int info=0;
	complx *wsp = (complx*)malloc((lwsp*4+lwsp*lwsp)*sizeof(complx));
	complx *eigs = wsp+3*lwsp;
	complx *eigvecs = eigs+lwsp;
	int *map = (int*)malloc(lwsp*sizeof(int));
	lwsp *= 3;
	double *rwork = (double*)malloc(lwsp*sizeof(double));
	const int ione = 1;
	if (!wsp || !map || !rwork) {
		fprintf(stderr,"Error: blk_cm_diag_sort - out of memory\n");
		exit(1);
	}

	mat_complx *dm = complx_matrix(Ud->dim,Ud->dim,MAT_DENSE_DIAG,0,Ud->basis);
	complx *dmpos = dm->data;
	T = (blk_mat_complx*)malloc(sizeof(blk_mat_complx));
	T->Nblocks = Ud->Nblocks;
	T->dim = Ud->dim;
	T->basis = Ud->basis;
	T->blk_dims = (int*)malloc(T->Nblocks*sizeof(int));
	T->m = (mat_complx*)malloc(T->Nblocks*sizeof(mat_complx));
	if (T->blk_dims == NULL || T->m == NULL) {
		fprintf(stderr,"Error: blk_cm_diag can not allocate T->m/blk_dims\n");
		exit(1);
	}

	for (i=0; i<Ud->Nblocks; i++) {
		NN = T->blk_dims[i] = Ud->blk_dims[i];
		bm = Ud->m + i;
		tm = T->m + i;
		if (NN == 1) {
			*dmpos = bm->data[0];
			create_sqmat_complx(tm,1,MAT_DENSE_DIAG,0,Ud->basis);
			tm->data[0] = Cunit;
		} else {
			switch (bm->type) {
			case MAT_DENSE:
				//cm_print(bm,"Propagator before");
				create_sqmat_complx(tm,NN,MAT_DENSE,0,Ud->basis);
				//zgeev_("N","V",&NN,bm->data,&NN,dmpos,NULL,&ione,tm->data,&NN,wsp,&lwsp,rwork,&info);
				zgeev_("N","V",&NN,bm->data,&NN,eigs,NULL,&ione,eigvecs,&NN,wsp,&lwsp,rwork,&info);
				if ( info != 0) {
					fprintf(stderr,"blk_cm_diag error: diagonalization failed for submatrix %d with the code '%d'\n", i, info);
				    exit(1);
				}
				//for (j=0; j<NN; j++) printf(" (%g %g)",eigs[j].re,eigs[j].im);
				//printf("\n---\n");
				eigsort(eigs,rwork,NN,map);
				// print the maping
				//printf("\nmap:");
				//for (j=0; j<NN; j++) printf(" %d",map[j]);
				//printf("\n");
				for (j=0; j<NN; j++) {
					dmpos[j] = eigs[map[j]];
					memcpy(tm->data+j*NN,eigvecs+map[j]*NN,NN*sizeof(complx));
				}
				//for (j=0; j<NN; j++) printf(" (%g %g)",dmpos[j].re,dmpos[j].im);
				//printf("\n+++++\n");
				//cm_print(tm,"Trasf. matrix");
				break;
			case MAT_DENSE_DIAG:
				create_sqmat_complx(tm,NN,MAT_DENSE_DIAG,0,Ud->basis);
				cm_unit(tm);
				memcpy(dmpos,bm->data,NN*sizeof(complx));
				break;
			case MAT_SPARSE:
				printf("WARNING!!! In order to diagonalize the submatrix %d (dim = %d) I need to \n"
						"convert it from sparse to full format. The whole operation may take long time...\n",i,NN);
				cm_dense(bm);
				create_sqmat_complx(tm,NN,MAT_DENSE,0,Ud->basis);
				//zgeev_("N","V",&NN,bm->data,&NN,dmpos,NULL,&ione,tm->data,&NN,wsp,&lwsp,rwork,&info);
				zgeev_("N","V",&NN,bm->data,&NN,eigs,NULL,&ione,eigvecs,&NN,wsp,&lwsp,rwork,&info);
				if ( info != 0) {
					fprintf(stderr,"blk_cm_diag error: diagonalization failed for submatrix %d with the code '%d'\n", i, info);
				    exit(1);
				}
				eigsort(eigs,rwork,NN,map);
				for (j=0; j<NN; j++) {
					dmpos[j] = eigs[map[j]];
					memcpy(tm->data+j*NN,eigvecs+map[j]*NN,NN*sizeof(complx));
				}
				//free(bm->irow);
				//free(bm->icol);
				break;
			case MAT_SPARSE_DIAG:
				fprintf(stderr,"Error: blk_cm_diag - block %d of MAT_SPASE_DIAG not implemented\n",i);
				exit(1);
			default:
				fprintf(stderr,"Error: blk_cm_diag - unknown submatrix %d type %d\n",i,bm->type);
				exit(1);
			}
		}
		dmpos += NN;

		free(bm->data);
	}

	free(wsp);
	free(rwork);
	free(map);
	free(Ud->m);
	Ud->m = dm;
	Ud->Nblocks = 1;
	int *ddd = realloc(Ud->blk_dims,sizeof(int));
	if (ddd == NULL) {
		fprintf(stderr,"Error: blk_cm_diag - can not reallocate Ud->blk_dims\n");
		exit(1);
	}
	Ud->blk_dims = ddd;
	Ud->blk_dims[0] = Ud->dim;

	return T;
}

blk_mat_complx * blk_cm_power(blk_mat_complx *obj, int n)
{
	int i, NN;
	blk_mat_complx *res = (blk_mat_complx*)malloc(sizeof(blk_mat_complx));
	mat_complx *m, *mm;

	res->Nblocks = obj->Nblocks;
	res->dim = obj->dim;
	res->basis = obj->basis;
	res->blk_dims = (int*)malloc(res->Nblocks*sizeof(int));
	res->m = (mat_complx*)malloc(res->Nblocks*sizeof(mat_complx));
	if (res->blk_dims == NULL || res->m == NULL) {
		fprintf(stderr,"Error: blk_cm_power can not allocate m/blk_dims\n");
		exit(1);
	}

	for (i=0; i<res->Nblocks; i++) {
		res->blk_dims[i] = NN = obj->blk_dims[i];
		m = res->m + i;
		if (NN == 1) {
			create_sqmat_complx(m,1,MAT_DENSE_DIAG,1,res->basis);
			m->data[0] = CRpow(obj->m[i].data[0],(double)n);
		} else {
			mm = cm_power(obj->m+i,n);
			m->data = mm->data;
			m->row = mm->row;
			m->col = mm->col;
			m->irow = mm->irow;
			m->icol = mm->icol;
			m->type = mm->type;
			m->basis = res->basis;
		}
	}

	return res;
}

mat_complx * cm_get_diagblock(mat_complx *cm, int sft,int dim)
{
	mat_complx *res;
	int i;

	switch (cm->type) {
	case MAT_DENSE:
		res = complx_matrix(dim,dim,MAT_DENSE,0,cm->basis);
		for (i=0; i<dim; i++) {
			memcpy(res->data+i*dim,cm->data+sft+sft*cm->row + i*cm->row,dim);
		}
		break;
	case MAT_DENSE_DIAG:
		res = complx_matrix(dim,dim,MAT_DENSE_DIAG,0,cm->basis);
		memcpy(res->data,cm->data+sft,dim);
		break;
	case MAT_SPARSE: {
		uint64_t nnzmax = (uint64_t)floor((uint64_t)dim*dim*(1.0-SPARSITY));
		int nnz = 0;
		if (dim<MAXFULLDIM) {
			res = complx_matrix(dim,dim,MAT_DENSE,0,cm->basis);
			cm_zero(res);
		} else {
			res = complx_matrix(dim,dim,MAT_SPARSE,nnzmax,cm->basis);
			res->irow[0] = 1;
		}
		for (i=0; i<dim; i++) {
			if (res->type == MAT_SPARSE) res->irow[i+1] = res->irow[i];
			int j;
			for (j=cm->irow[sft+i]; j<cm->irow[sft+i+1]; j++) {
				int c = cm->icol[j-1];
				if (c>sft && c<=sft+dim) {
					if (res->type == MAT_SPARSE) {
						res->icol[nnz] = c-sft;
						res->data[nnz] = cm->data[j-1];
						res->irow[i+1]++;
						nnz++;
						if (nnz > nnzmax) cm_dense(res);
					} else {
						res->data[i+(c-sft-1)*dim] = cm->data[j-1];
					}
				}
			}
		}
		if (res->type == MAT_SPARSE) cm_change_nnz(res,nnz);
		break; }
	case MAT_SPARSE_DIAG:
		res = complx_matrix(dim,dim,MAT_DENSE_DIAG,0,cm->basis);
		for (i=0; i<dim; i++) {
			res->data[i] = cm_getelem(cm,sft+i+1,sft+i+1);
		}
		break;
	default:
		fprintf(stderr,"Error: cm_get_diagblock - unknown matrix type\n");
		exit(1);
	}
	return res;
}

complx blk_cm_trace_adjoint(blk_mat_complx *A, blk_mat_complx *B, Sim_info *sim)
{
	complx res;

	if (A->basis == B->basis) {
		// same basis
		if (A->Nblocks == 1 && B->Nblocks == 1) {
			res = cm_trace_adjoint(A->m, B->m);
		} else if (A->Nblocks == 1) {
			int i, dim, sft=0;
			Czero(res);
			for (i=0; i<B->Nblocks; i++) {
				dim = B->blk_dims[i];
				mat_complx *mx = cm_get_diagblock(A->m,sft,dim);
				complx z = cm_trace_adjoint(mx,B->m+i);
				free_complx_matrix(mx);
				res.re += z.re;
				res.im += z.im;
				sft += dim;
			}
		} else if (B->Nblocks == 1) {
			int i, dim, sft=0;
			Czero(res);
			for (i=0; i<A->Nblocks; i++) {
				dim = A->blk_dims[i];
				mat_complx *mx = cm_get_diagblock(B->m,sft,dim);
				complx z = cm_trace_adjoint(A->m+i,mx);
				free_complx_matrix(mx);
				res.re += z.re;
				res.im += z.im;
				sft += dim;
			}
		} else {
			assert(A->Nblocks == B->Nblocks);
			int i;
			Czero(res);
			for (i=0; i<A->Nblocks; i++) {
				complx z = cm_trace_adjoint(A->m+i, B->m+i);
				res.re += z.re;
				res.im += z.im;
			}
		}
	} else {
		// different basis
		Czero(res);
		if (A->Nblocks >= B->Nblocks) {
			// do the job in basis of A
			int *permvec = sim->perm_table[B->basis + A->basis*sim->Nbasis];
			assert(LEN(permvec) == sim->matdim);
			int i, sft=0;
			for (i=0; i<A->Nblocks; i++) {
				int dim = A->blk_dims[i];
				mat_complx *mx = A->m+i;
				switch (mx->type) {
				case MAT_DENSE: {
					int i2, j2;
					complx *zz = mx->data;
					for (i2=0; i2<dim; i2++) {
						for (j2=0; j2<dim; j2++) { // row in mx
							complx z = blk_cm_getelem(B,permvec[sft+j2+1],permvec[sft+i2+1]);
							res.re += z.re*zz->re + z.im*zz->im; // complex conjugate of mx
							res.im += z.im*zz->re - z.re*zz->im;
							zz++;
						}
					}
					break; }
				case MAT_DENSE_DIAG: {
					complx *zz = mx->data;
					int i2;
					for (i2=0; i2<dim; i2++) {
						complx z = blk_cm_getelem(B,permvec[sft+i2+1],permvec[sft+i2+1]);
						res.re += z.re*zz->re + z.im*zz->im; // complex conjugate of mx
						res.im += z.im*zz->re - z.re*zz->im;
						zz++;
					}
					break; }
				case MAT_SPARSE:
				case MAT_SPARSE_DIAG: {
					complx *zz = mx->data;
					int *ic = mx->icol;
					int i2, j2;
					for (i2=0; i2<dim; i2++) {
						for (j2=mx->irow[i2]; j2<mx->irow[i2+1]; j2++) {
							complx z = blk_cm_getelem(B,permvec[sft+i2+1],permvec[sft+(*ic)]);
							res.re += z.re*zz->re + z.im*zz->im; // complex conjugate of mx
							res.im += z.im*zz->re - z.re*zz->im;
							zz++;
							ic++;
						}
					}
					break; }
				default:
					fprintf(stderr,"Error: blk_cm_trace_adjoint - unknown matrix A type\n");
					exit(1);
				}
				sft += dim;
			}
		} else {
			// do the job in basis of B
			int *permvec = sim->perm_table[A->basis + B->basis*sim->Nbasis];
			assert(LEN(permvec) == sim->matdim);
			int i, sft=0;
			for (i=0; i<B->Nblocks; i++) {
				int dim = B->blk_dims[i];
				mat_complx *mx = B->m+i;
				switch (mx->type) {
				case MAT_DENSE: {
					int i2, j2;
					complx *zz = mx->data;
					for (i2=0; i2<dim; i2++) {
						for (j2=0; j2<dim; j2++) { // row in mx
							complx z = blk_cm_getelem(A,permvec[sft+j2+1],permvec[sft+i2+1]);
							res.re += z.re*zz->re + z.im*zz->im; // complex conjugate of A element
							res.im += -z.im*zz->re + z.re*zz->im;
							zz++;
						}
					}
					break; }
				case MAT_DENSE_DIAG: {
					complx *zz = mx->data;
					int i2;
					for (i2=0; i2<dim; i2++) {
						complx z = blk_cm_getelem(A,permvec[sft+i2+1],permvec[sft+i2+1]);
						res.re += z.re*zz->re + z.im*zz->im; // complex conjugate of A element
						res.im += -z.im*zz->re + z.re*zz->im;
						zz++;
					}
					break; }
				case MAT_SPARSE:
				case MAT_SPARSE_DIAG: {
					complx *zz = mx->data;
					int *ic = mx->icol;
					int i2, j2;
					for (i2=0; i2<dim; i2++) {
						for (j2=mx->irow[i2]; j2<mx->irow[i2+1]; j2++) {
							complx z = blk_cm_getelem(A,permvec[sft+i2+1],permvec[sft+(*ic)]);
							res.re += z.re*zz->re + z.im*zz->im; // complex conjugate of A element
							res.im += -z.im*zz->re + z.re*zz->im;
							zz++;
							ic++;
						}
					}
					break; }
				default:
					fprintf(stderr,"Error: blk_cm_trace_adjoint - unknown matrix B type\n");
					exit(1);
				}
				sft += dim;
			}
		}
	}

	return res;
}

int blk_dm_nnz(blk_mat_double *blkm)
{
	int i, nnz = 0;

	for (i=0; i<blkm->Nblocks; i++) {
		if (blkm->blk_dims[i] == 1) {
			if (fabs(blkm->m[i].data[0]) > TINY) nnz++;
		} else {
			nnz += dm_nnz(blkm->m+i);
		}
	}
	return nnz;
}

int blk_cm_nnz(blk_mat_complx *blkm)
{
	int i, nnz = 0;

	for (i=0; i<blkm->Nblocks; i++) {
		if (blkm->blk_dims[i] == 1) {
			if (fabs(blkm->m[i].data[0].re) > TINY && fabs(blkm->m[i].data[0].im) > TINY) nnz++;
		} else {
			nnz += cm_nnz(blkm->m+i);
		}
	}
	return nnz;
}

blk_mat_complx * blk_cm_ln(blk_mat_complx *blkm)
{
	int i, N;
	mat_complx *m, *mln;
	blk_mat_complx *res = create_blk_mat_complx_copy(blkm);

	for (i=0; i<blkm->Nblocks; i++) {
		N = blkm->blk_dims[i];
		m = blkm->m + i;
		mln = res->m + i;
		if (N==1) {
			mln->data[0] = Clog(m->data[0]);
		} else {
			mat_complx *dum = cm_ln(m);
			cm_swap_innards_and_destroy(mln,dum);
		}
	}
	return res;
}

