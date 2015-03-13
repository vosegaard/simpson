/*
 * cm.c
 *
 *  Created on: Jun 11, 2010
 *      Author: zdenek
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>
#include "matrix.h"
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
#if !defined(INTEL_MKL) && !defined(__APPLE__)
	void zgeev_(const char *jobl, const char *jobr, const int *d1, void *data,const int *d2,
			void *eigs, void *vl, const int *dvl,void *vr, const int *dvr, void *wsp,
			const int *lwsp, double *rwork, int *info);
	void dsyev_(const char *job,const char *uplo, const int *K, double * dT, const int *L, double *eigs, double *dwsp, const int *ldwsp, int *info);
	void zlarcm_(const int *K,const int *L,double *a,const int *M,void *cm,const int *ldcm,void *res,const int *ldres, double *dwsp);
#endif
#include "defs.h"
#include "auxmath.h"

	/* for acurate timings on windows */
//#define TIMING
#include "timing.h"

/* vector operations */
void dv_muld(double *vec, double d)
{
	//cblas_dscal(LEN(vec),d,&vec[1],1);
	cblas_dscal(LEN(vec),d,vec+1,1);
}
void dv_multod(double *v1, double *v2, double d)
{
	int len = LEN(v1);
	if (len != LEN(v2)) {
		fprintf(stderr,"Error in dv_multod: dimension mismatch\n");
		exit(1);
	}
	//cblas_daxpy(len,d,&v2[1],1,&v1[1],1);
	cblas_daxpy(len,d,v2+1,1,v1+1,1);
}

void dv_zero(double *v)
{
	//memset(&v[1],0,LEN(v)*sizeof(double));
	memset(v+1,0,LEN(v)*sizeof(double));
}

void cv_muld(complx *vec, double d)
{
	//cblas_zdscal(LEN(vec),d,&vec[1],1);
	cblas_zdscal(LEN(vec),d,vec+1,1);
}

void cv_mulc(complx *vec, complx z)
{
	cblas_zscal(LEN(vec),&z,vec+1,1);
}

void cv_multod(complx *v1, complx *v2, double d)
{
	int len = LEN(v1);
	if (len != LEN(v2)) {
		fprintf(stderr,"Error in cv_multod: dimension mismatch\n");
		exit(1);
	}
	//cblas_daxpy(len*2,d,(double*)(&v2[1]),1,(double*)(&v1[1]),1);
	cblas_daxpy(len*2,d,(double*)(v2+1),1,(double*)(v1+1),1);
}

void cv_multoc(complx *v1, complx *v2, complx z)
{
	int len = LEN(v1);
	if (len != LEN(v2)) {
		fprintf(stderr,"Error in cv_multoc: dimension mismatch\n");
		exit(1);
	}
	cblas_zaxpy(len,&z,&(v2[1]),1,&(v1[1]),1);
}

void cv_conj(complx *vec)
{
	//double *ptr = (double*)(&vec[1]) +1;
	double *ptr = (double*)(vec+1) +1;
	cblas_dscal(LEN(vec),-1.0,ptr,2);
}

double cv_asum(complx *vec)
{
	double res;

	res = cblas_dzasum(LEN(vec),&(vec[1]),1);
	return res;
}

int * iv_dup(int *vec)
{
	int *res;
	int len = LEN(vec);

	res = int_vector(len);
	memcpy(&(res[1]),&(vec[1]),len*sizeof(int));

	return res;
}


int iv_max(int *vec)
{
	int i, res = -1000000;

	for (i=1; i<=LEN(vec); i++) {
		if (vec[i] > res) res = vec[i];
	}
	return res;
}

/* remains to be tested!!! */
void cv_matmulto(complx *vec, mat_complx *A)
{
	int dim = LEN(vec);
	if ( dim != A->col) {
		fprintf(stderr,"Error: cv_matmulto - dimension mismatch A(%d,%d).v(%d)\n",A->row,A->col,dim);
		exit(1);
	}

	switch (A->type) {
	case MAT_DENSE: {
		complx *res = malloc(dim*sizeof(complx));
		cblas_zgemv(CblasColMajor,CblasNoTrans,A->row,A->col,&Cunit,A->data,A->row,&(vec[1]),1,&Cnull,res,1);
		memcpy(&(vec[1]),res,dim*sizeof(complx));
		free(res);
		break; }
	case MAT_DENSE_DIAG: {
		int i;
		complx *vval=vec+1, *Aval=A->data;
		for (i=0; i<dim; i++) {
			double re = vval->re;
			vval->re = re*Aval->re - vval->im*Aval->im;
			vval->im = re*Aval->im + vval->im*Aval->re;
			vval++;
			Aval++;
		}
		break;}
	case MAT_SPARSE_DIAG: {
		int i;
		complx *vval=vec+1, Aval;
		for (i=1; i<=dim; i++) {
			double re = vval->re;
			Aval = cm_getelem(A,i,i);
			vval->re = re*Aval.re - vval->im*Aval.im;
			vval->im = re*Aval.im + vval->im*Aval.re;
			vval++;
		}
		break;}
	case MAT_SPARSE: {
#ifdef INTEL_MKL
		complx *res = (complx*)malloc(dim*sizeof(complx));
		if (A->row == A->col) {
			mkl_zcsrgemv("N",&dim,A->data,A->irow,A->icol,&(vec[1]),res);
		} else {
			mkl_zcsrmv("N",&(A->row),&(A->col),&Cunit,"GNNFUU",A->data,A->icol,A->irow,&(A->irow[1]),&(vec[1]),&Cnull,res);
		}
		memcpy(&(vec[1]),res,dim*sizeof(complx));
		free(res);
#else
		fprintf(stderr,"Error: cv_matmulto - sparse algebra not compiled\n");
		exit(1);
#endif
		break;}
	default:
		fprintf(stderr,"Error: cv_matmulto - invalid matrix type\n");
		exit(1);
	}
}

/* remains to be tested!!! res = A*vec */
void cv_matmul(complx *res, mat_complx *A, complx *vec)
{
	int dim = LEN(vec);
	if ( dim != A->col || dim != LEN(res)) {
		fprintf(stderr,"Error: cv_matmulto - dimension mismatch res(%d)=A(%d,%d).v(%d)\n",LEN(res),A->row,A->col,dim);
		exit(1);
	}

	switch (A->type) {
	case MAT_DENSE: {
		cblas_zgemv(CblasColMajor,CblasNoTrans,A->row,A->col,&Cunit,A->data,A->row,&(vec[1]),1,&Cnull,&(res[1]),1);
		break; }
	case MAT_DENSE_DIAG: {
		int i;
		complx *vval=vec+1, *resval=res+1, *Aval=A->data;
		for (i=0; i<dim; i++) {
			resval->re = vval->re*Aval->re - vval->im*Aval->im;
			resval->im = vval->re*Aval->im + vval->im*Aval->re;
			resval++;
			vval++;
			Aval++;
		}
		break;}
	case MAT_SPARSE_DIAG: {
		int i;
		complx *vval=vec+1, *resval=res+1, Aval;
		for (i=1; i<=dim; i++) {
			Aval = cm_getelem(A,i,i);
			resval->re = vval->re*Aval.re - vval->im*Aval.im;
			resval->im = vval->re*Aval.im + vval->im*Aval.re;
			resval++;
			vval++;
		}
		break;}
	case MAT_SPARSE: {
#ifdef INTEL_MKL
		if (A->row == A->col) {
			mkl_zcsrgemv("N",&dim,A->data,A->irow,A->icol,&(vec[1]),&(res[1]));
		} else {
			mkl_zcsrmv("N",&(A->row),&(A->col),&Cunit,"GNNFUU",A->data,A->icol,A->irow,&(A->irow[1]),&(vec[1]),&Cnull,&(res[1]));
		}
#else
		fprintf(stderr,"Error: cv_matmulto - sparse algebra not compiled\n");
		exit(1);
#endif
		break;}
	default:
		fprintf(stderr,"Error: cv_matmul - invalid matrix type\n");
		exit(1);
	}
}

complx cv_dotu(complx *a, complx *b)
{
	complx res;
	cblas_zdotu_sub(LEN(a),&(a[1]),1,&(b[1]),1,&res);
	return res;
}

complx cv_dotc(complx *a, complx *b)
{
	complx res;
	cblas_zdotc_sub(LEN(a),&(a[1]),1,&(b[1]),1,&res);
	return res;
}

double cv_norm(complx *vec)
{
	complx res;

	cblas_zdotc_sub(LEN(vec),&(vec[1]),1,&(vec[1]),1,&res);
	return sqrt(res.re);
}

complx * cv_dup(complx *vec)
{
	complx *res;
	int len = LEN(vec);

	res = complx_vector(len);
	memcpy(&(res[1]),&(vec[1]),len*sizeof(complx));

	return res;
}

void cv_copy(complx *dest, complx *vec)
{
	int len = LEN(vec);

	if (len != LEN(dest)) {
		fprintf(stderr,"Error: cv_copy - dimension mismatch dest(%d), orig(%d)\n",LEN(dest),LEN(vec));
		exit(1);
	}
	memcpy(&(dest[1]),&(vec[1]),len*sizeof(complx));
}

void cv_print(complx *vec, char *title)
{
	int i;
	printf("\n%s: complex vector (%d)\n",title,LEN(vec));
	for (i=1; i<=LEN(vec); i++) printf("(%9.6g , %9.6g )\n",vec[i].re,vec[i].im);
	printf("\n");
}

/****************************************************/

/* double matrix operations */

void dm_swap_innards_and_destroy(mat_double *dest,mat_double *from)
{
	free((char*)(dest->data));
	dest->data = from->data;
	dest->type = from->type;
	dest->row = from->row;
	dest->col = from->col;
	dest->basis = from->basis;
	if (dest->irow != NULL) free((char*)(dest->irow));
	dest->irow = from->irow;
	if (dest->icol != NULL) free((char*)(dest->icol));
	dest->icol = from->icol;
	free((char*)from);
}

void dm_zero(mat_double *m)
{
	switch (m->type) {
	case MAT_DENSE :
		memset(m->data, 0, (m->row * m->col) * sizeof(double));
		break;
	case MAT_DENSE_DIAG :
		memset(m->data, 0, m->row * sizeof(double));
		break;
	case MAT_SPARSE :
	case MAT_SPARSE_DIAG : {
		//DEBUGPRINT("dm_zero sparse:");
		if (m->irow[m->row] != 0) {
			//DEBUGPRINT(" reallocation and");
			double * new_data = (double*)realloc(m->data, sizeof(double));
			MKL_INT * new_icol = (MKL_INT*)realloc(m->icol, sizeof(MKL_INT));
			assert( (new_data != NULL) && (new_icol != NULL) );
			m->data = new_data;
			m->icol = new_icol;
		}
		//DEBUGPRINT(" filling\n");
		m->data[0] = 0.0;
		m->irow[0] = 1;
		int i; for (i=1; i<=m->row; i++) m->irow[i] = 2;
		m->icol[0] = 1;
		break; }
   default :
	   fprintf(stderr,"Error: dm_zero - unknown type '%d'\n",m->type);
	   exit(1);
	}
}

mat_double * dm_dup(mat_double *m)
{
	int len;
	mat_double *mm;


	switch (m->type) {
	case MAT_DENSE :
		mm = double_matrix(m->row,m->col,m->type,0,m->basis);
		memcpy(mm->data,m->data,(m->row*m->col)*sizeof(double));
		break;
	case MAT_DENSE_DIAG :
		mm = double_matrix(m->row,m->col,m->type,0,m->basis);
		memcpy(mm->data,m->data,(m->row)*sizeof(double));
		break;
	case MAT_SPARSE :
	case MAT_SPARSE_DIAG :
		len = m->irow[m->row] -1;
		mm = double_matrix(m->row,m->col,m->type,len,m->basis);
		memcpy(mm->data,m->data,len*sizeof(double));
		memcpy(mm->icol, m->icol, len*sizeof(MKL_INT));
		memcpy(mm->irow, m->irow, (m->row+1)*sizeof(MKL_INT));
		break;
   default :
	   fprintf(stderr,"Error: dm_dup - unknown type '%d'\n",m->type);
	   exit(1);
	}

	return mm;
}

/* result is always full dense matrix */
mat_double * dm_dup2(mat_double *m)
{
	int i, j, N, c, *ic;
	double *ddata, *dd;
	int len = m->row;
	mat_double *mm = double_matrix(m->row,m->col,MAT_DENSE,0,m->basis);

	switch (m->type) {
		case MAT_DENSE :
			memcpy(mm->data,m->data,(m->row*m->col)*sizeof(double));
			break;
		case MAT_DENSE_DIAG :
			dm_zero(mm);
			cblas_dcopy(len,m->data,1,mm->data,len+1);
			break;
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG :
			dm_zero(mm);
			ddata = mm->data;
			ic = m->icol;
			dd = m->data;
			for (i=0; i<m->row; i++) {
				N = m->irow[i+1] - m->irow[i];
				for (j=0;j<N; j++) {
					c = (*ic) -1;
					ic++;
					ddata[i+c*m->row] = *dd;
					dd++;
				}
			}
			break;
	   default :
		   fprintf(stderr,"Error: dm_dup2 - unknown type '%d'\n",m->type);
		   exit(1);
		}

	return mm;
}



void dm_copy(mat_double *m1, mat_double *m2)
{
	if (m1 == NULL) {
		fprintf(stderr,"Error: dm_copy - destination is NULL\n");
		exit(1);
	}
	if ( (m1->type != m2->type) || (m1->row != m2->row) || (m1->col != m2->col) ) {
		DEBUGPRINT("dm_copy: different matrix types or dimensions (%d,%d)<-(%d,%d)\n",m1->row,m1->col,m2->row,m2->col);
		mat_double *dum = dm_dup(m2);
		dm_swap_innards_and_destroy(m1,dum);
		return;
	}
	switch (m1->type) {
		case MAT_DENSE :
			memcpy(m1->data,m2->data,m2->row*m2->col*sizeof(double));
			break;
		case MAT_DENSE_DIAG :
			memcpy(m1->data,m2->data,m2->row*sizeof(double));
			break;
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG : {
			int nnz = m2->irow[m2->row] - 1;
			dm_change_nnz(m1,nnz);
			memcpy(m1->data,m2->data,nnz*sizeof(double));
			memcpy(m1->icol,m2->icol,nnz*sizeof(MKL_INT));
			memcpy(m1->irow,m2->irow,(m2->row+1)*sizeof(MKL_INT));
			break; }
	   default :
		   fprintf(stderr,"Error: dm_copy - unknown type '%d'\n",m1->type);
		   exit(1);
	}
	m1->basis = m2->basis;
}

void dm_dense(mat_double *m)
{
	int i, j, N, c, *ic;
	double *ddata, *dd;

	switch (m->type) {
		case MAT_DENSE :
		case MAT_DENSE_DIAG :
			break;
		case MAT_SPARSE :
			ddata = (double*)calloc(m->row*m->col,sizeof(double));
			ic = m->icol;
			dd = m->data;
			for (i=0; i<m->row; i++) {
				N = m->irow[i+1] - m->irow[i];
				for (j=0;j<N; j++) {
					c = (*ic) -1;
					ic++;
					ddata[i+c*m->row] = *dd;
					dd++;
				}
			}
			free(m->icol); m->icol = NULL;
			free(m->irow); m->irow = NULL;
			free(m->data); m->data = ddata;
			m->type = MAT_DENSE;
			break;
		case MAT_SPARSE_DIAG :
			ddata = (double*)calloc(m->row,sizeof(double));
			N = m->irow[m->row] - 1;
			ic = m->icol;
			dd = m->data;
			for (i=0; i<N; i++) {
				ddata[*ic -1] = *dd;
				ic++;
				dd++;
			}
			free(m->icol); m->icol = NULL;
			free(m->irow); m->irow = NULL;
			free(m->data); m->data = ddata;
			m->type = MAT_DENSE_DIAG;
			break;
	   default :
		   fprintf(stderr,"Error: dm_dense - unknown type '%d'\n",m->type);
		   exit(1);
		}
}

void dm_sparse(mat_double *m, double tol)
{
	int realloc_size;
	double *new_data;
	MKL_INT *new_icol;
	uint64_t Nmax = (uint64_t)floor((uint64_t)m->row*m->col*(1.0-SPARSITY));

	//printf("dm_sparse: %d x %d = %" PRId64 ", sp = %g, Nmax = %" PRId64 "\n",m->row,m->col,(uint64_t)m->row*m->col,SPARSITY,Nmax);
	//printf("dm_sparse will allow %d NNZ\n",Nmax);
	switch (m->type) {
		case MAT_DENSE : {
			int i, j, n=0, rc=1;
			double *d1, *d2;
			MKL_INT *ir, *ic;

			m->irow = ir = (MKL_INT*)malloc((m->row+1)*sizeof(MKL_INT));
			m->icol = ic = (MKL_INT*)malloc(Nmax*sizeof(MKL_INT));
			new_data = d2 = (double*)malloc(Nmax*sizeof(double));
		    *ir = rc;
		    for (i=0; i<m->row; i++) {
		    	d1 = m->data + i;
			    for (j=0; j<m->col; j++) {
			    	if ( fabs(*d1) >= tol ) {
			    		if (++n <= Nmax) {
			    			*d2 = *d1;
			    			d2++;
			    			*ic = j+1;
			    			ic++;
			    			rc++;
			    		}
			    	} else {
			    		*d1 = 0.0;
			    	}
			    	d1 += m->row;
		    	}
		    	ir++;
		    	*ir = rc;
		    }
		    if (n>Nmax) {
		    	DEBUGPRINT("dm_sparse: dense matrix has %d NNZ (%d allowed) - remains dense\n",n,Nmax);
		    	free(new_data);
		    	free(m->irow); m->irow = NULL;
		    	free(m->icol); m->icol = NULL;
		    	return;
		    }
		    free(m->data);
		    m->data = new_data;
		    new_data = NULL;
		    m->type = MAT_SPARSE;
			break; }
		case MAT_DENSE_DIAG : {
			int r, rr=1;
			double *d1, *d2;
			MKL_INT *ic;
			m->type = MAT_SPARSE_DIAG;
			m->icol = ic = (MKL_INT*)malloc(m->row*sizeof(MKL_INT));
			m->irow = (MKL_INT*)malloc((m->row+1)*sizeof(MKL_INT));
			m->irow[0] = rr;
			d1 = d2 = m->data;
			for (r=1; r<=m->row; r++) {
				if ( fabs(*d2) < tol) {
					d2++;
				} else {
					if (d1 != d2) {
						*d1 = *d2;
					}
					*ic = r;
					d1++;
					d2++;
					rr++;
					ic++;
				}
				m->irow[r] = rr;
				//m->icol[r-1] = r;
			}
			break; }
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG : {
			double *d1, *d2;
			MKL_INT *c1, *c2, *r1, *r2;
			int nc, ncnew, rr, ic, ir;

			r1 = m->irow;
			r2 = r1+1;
			rr = *r1;
			d1 = d2 = m->data;
			c1 = c2 = m->icol;
			for (ir=0; ir<m->row; ir++) {
				nc = ncnew = *r2 - rr;
				for (ic=0; ic<nc; ic++) {
					if (fabs(*d2) < tol) {
						d2++;
						c2++;
						ncnew--;
					} else {
						if (d1 != d2) {
							*d1 = *d2;
							*c1 = *c2;
						}
						d1++;
						d2++;
						c1++;
						c2++;
					}
				}
				rr = *r2;
				*r2 = *r1 + ncnew;
				r1++;
				r2++;
			}
			break; }
	   default :
		   fprintf(stderr,"Error: dm_sparse - unknown type '%d'\n",m->type);
		   exit(1);
	}
	realloc_size = m->irow[m->row]-1;
	assert(realloc_size >= 0);
	if (realloc_size > Nmax) {
		DEBUGPRINT("dm_sparse: %d NNZ was detected but %d is allowed - converting to dense matrix\n",realloc_size,Nmax);
		dm_dense(m);
		return;
	}
	//DEBUGPRINT("dm_sparse nnz = %d\n",realloc_size);
	if (realloc_size == 0) {
		realloc_size = 1;
		m->data[0] = 0.0;
		m->icol[0] = m->col;
		m->irow[m->row] = 2;
	}
	new_data = (double*)realloc(m->data, realloc_size*sizeof(double));
	new_icol = (MKL_INT*)realloc(m->icol, realloc_size*sizeof(MKL_INT));
	assert( (new_data != NULL) && (new_icol != NULL) );
	m->data = new_data;
	m->icol = new_icol;
}

mat_double * dm_creatediag(const char c, int dim, double d, int basis)
{
	mat_double * res;
	int i;

	if ( (c == 'd') || (c == 'D') ) {
		res = double_matrix(dim,dim,MAT_DENSE_DIAG,dim,basis);
		for (i=0; i<dim; i++) res->data[i] = d;
		return res;
	}
	if ( (c == 's') || (c == 'S') ) {
		res = double_matrix(dim,dim,MAT_SPARSE_DIAG,dim,basis);
		for (i=0; i<dim; i++) {
			res->data[i] = d;
			res->icol[i] = res->irow[i] = i+1;
		}
		res->irow[dim] = dim+1;
		return res;
	}
    fprintf(stderr,"Error: dm_creatediag called with wrong code '%c'\n",c);
    exit(1);
    return NULL;
}

mat_complx * dm_complx(mat_double *m)
{
	mat_complx * res;
	int len=0;

	switch (m->type) {
	case MAT_DENSE :
		res = complx_matrix(m->row, m->col, m->type, 0,m->basis);
		len = m->row * m->col;
		break;
	case MAT_DENSE_DIAG :
		res = complx_matrix(m->row,m->col,m->type,0,m->basis);
		len = m->row;
		break;
	case MAT_SPARSE :
	case MAT_SPARSE_DIAG :
		len = m->irow[m->row] -1;
		res = complx_matrix(m->row,m->col,m->type,len,m->basis);
		memcpy(res->icol, m->icol, len*sizeof(MKL_INT));
		memcpy(res->irow, m->irow, (m->row+1)*sizeof(MKL_INT));
		break;
   default :
	   fprintf(stderr,"Error: dm_complx - unknown type '%d'\n",m->type);
	   exit(1);
	}
	memset(res->data,0,len*sizeof(complx));
	cblas_dcopy(len,m->data,1,(double*)(res->data),2);
	return res;
}

void cm_swap_innards_and_destroy(mat_complx *dest, mat_complx *from)
{
	free((char*)(dest->data));
	dest->data = from->data;
	dest->type = from->type;
	dest->row = from->row;
	dest->col = from->col;
	dest->basis = from->basis;
	if (dest->irow != NULL) free((char*)(dest->irow));
	dest->irow = from->irow;
	if (dest->icol != NULL) free((char*)(dest->icol));
	dest->icol = from->icol;
	free((char*)from);
}

void dm_copy2cm(mat_double *dm, mat_complx *cm)
{
	int dlen;

	if (cm == NULL) {
		fprintf(stderr,"Error: dm_copy2cm - destination is NULL\n");
		exit(1);
	}

	if ( (dm->type != cm->type) || (dm->row != cm->row) || (dm->col != cm->col) ) {
		DEBUGPRINT("dm_copy2cm: different dimensions or types, reallocating cm\n");
		mat_complx *dum = dm_complx(dm);
		cm_swap_innards_and_destroy(cm,dum);
		return;
	}

	switch (dm->type) {
	case MAT_DENSE :
		dlen = cm->row*cm->col;
		memset(cm->data,0,dlen*sizeof(complx));
		break;
	case MAT_DENSE_DIAG :
		dlen = cm->row;
		memset(cm->data,0,dlen*sizeof(complx));
		break;
	case MAT_SPARSE :
	case MAT_SPARSE_DIAG :
		dlen = dm->irow[dm->row] -1;
		cm_change_nnz(cm,dlen);
		memset(cm->data,0,dlen*sizeof(complx));
		memcpy(cm->icol, dm->icol, dlen*sizeof(MKL_INT));
		memcpy(cm->irow, dm->irow, (dm->row+1)*sizeof(MKL_INT));
		break;
   default :
	   fprintf(stderr,"Error: dm_copy2cm - unknown type '%d'\n",dm->type);
	   exit(1);
	}
	cblas_dcopy(dlen,dm->data,1,(double*)(cm->data),2);
	cm->basis = dm->basis;
}

mat_complx * dm_imag(mat_double *m)
{
	mat_complx * res;
	int len;

	switch (m->type) {
	case MAT_DENSE :
		res = complx_matrix(m->row,m->col,m->type,0,m->basis);
		len = m->row * m->col;
		break;
	case MAT_DENSE_DIAG :
		res = complx_matrix(m->row,m->col,m->type,0,m->basis);
		len = m->row;
		break;
	case MAT_SPARSE :
	case MAT_SPARSE_DIAG :
		len = m->irow[m->row] -1;
		res = complx_matrix(m->row,m->col,m->type,len,m->basis);
		memcpy(res->icol, m->icol, len*sizeof(MKL_INT));
		memcpy(res->irow, m->irow, (m->row+1)*sizeof(MKL_INT));
		break;
   default :
	   fprintf(stderr,"Error: dm_complx - unknown type '%d'\n",m->type);
	   exit(1);
	}

	memset(res->data,0,len*sizeof(complx));
	cblas_dcopy(len,m->data,1,(double*)(res->data)+1,2);
	return res;
}

/* create MAT_DANSE with imaginary part from m */
mat_complx * dm_imag2(mat_double *m)
{
	int len;
	mat_complx * res = complx_matrix(m->row,m->col,MAT_DENSE,0,m->basis);
	cm_zero(res);

	switch (m->type) {
	case MAT_DENSE :
		len = m->row * m->col;
		cblas_dcopy(len,m->data,1,(double*)(res->data)+1,2);
		break;
	case MAT_DENSE_DIAG :
		len = m->row;
		cblas_dcopy(len,m->data,1,(double*)(res->data)+1,2*(len+1));
		break;
	case MAT_SPARSE :
	case MAT_SPARSE_DIAG : {
		int i, j;
		int *ic = m->icol;
		double *dd = m->data;
		for (i=0; i<m->row; i++) {
			len = m->irow[i+1] - m->irow[i];
			for (j=0;j<len; j++) {
				res->data[i+(*ic - 1)*m->row].im = *dd;
				ic++;
				dd++;
			}
		}
		break; }
   default :
	   fprintf(stderr,"Error: dm_imag2 - unknown type '%d'\n",m->type);
	   exit(1);
	}
	return res;
}

void dm_copy2cm_imag(mat_double *dm, mat_complx *cm)
{
	int dlen;

	if (cm == NULL) {
		fprintf(stderr,"Error: dm_copy2cm - destination is NULL\n");
		exit(1);
	}
	if ( (dm->type != cm->type) || (dm->row != cm->row) || (dm->col != cm->col) ) {
		DEBUGPRINT("dm_copy2cm_imag: different dimensions or types, reallocating cm\n");
		mat_complx *dum = dm_imag(dm);
		cm_swap_innards_and_destroy(cm,dum);
		return;
	}
	switch (dm->type) {
	case MAT_DENSE :
		dlen = cm->row*cm->col;
		memset(cm->data,0,dlen*sizeof(complx));
		break;
	case MAT_DENSE_DIAG :
		dlen = dm->row;
		memset(cm->data,0,dlen*sizeof(complx));
		break;
	case MAT_SPARSE :
	case MAT_SPARSE_DIAG :
		dlen = dm->irow[dm->row] -1;
		cm_change_nnz(cm,dlen);
		memset(cm->data,0,dlen*sizeof(complx));
		memcpy(cm->icol, dm->icol, dlen*sizeof(MKL_INT));
		memcpy(cm->irow, dm->irow, (dm->row+1)*sizeof(MKL_INT));
		break;
   default :
	   fprintf(stderr,"Error: dm_copy2cm - unknown type '%d'\n",dm->type);
	   exit(1);
	}
	cblas_dcopy(dlen,dm->data,1,(double*)(cm->data)+1,2);
	cm->basis = dm->basis;
}

void dm_muld(mat_double *m, double d)
{
	int len=0;

	switch (m->type) {
	case MAT_DENSE  : len = m->row * m->col; break;
	case MAT_DENSE_DIAG  : len = m->row; break;
	case MAT_SPARSE :
	case MAT_SPARSE_DIAG : len = m->irow[m->row]-1;
	}
	cblas_dscal(len,d,m->data,1);
}

/* sparse res = a + d*b */
mat_double * simpson_dcsradd(mat_double *a, double d, mat_double *b)
{
	assert(a->basis == b->basis);
#ifdef INTEL_MKL
		char trans='N';
		int sort=0;
		int request=1, info=0, nzmax=a->row*a->col;
		int type;
		mat_double *res;

		if ( (a->type == MAT_SPARSE) || (b->type == MAT_SPARSE) ) {
			type = MAT_SPARSE;
		} else {
			type = MAT_SPARSE_DIAG;
		}
		//DEBUGPRINT("   dcsradd types: a = %d, b = %d, res = %d\n",a->type,b->type,type);
		res = double_matrix(a->row,a->col,type,0,a->basis);
		mkl_dcsradd(&trans, &request, &sort, &(a->row), &(a->col), a->data, a->icol, a->irow, &d, b->data, b->icol, b->irow, res->data, res->icol, res->irow, &nzmax, &info);
		dm_change_nnz(res,res->irow[res->row]-1);
		request = 2;
		mkl_dcsradd(&trans, &request, &sort, &(a->row), &(a->col), a->data, a->icol, a->irow, &d, b->data, b->icol, b->irow, res->data, res->icol, res->irow, &nzmax, &info);
		if (info != 0) {
			fprintf(stderr,"Error: dcsradd failed with the code '%d'\n",info);
			exit(1);
		}
		return res;
#else
		fprintf(stderr,"Error: sparse matrix routine dcsradd not compiled\n");
		exit(1);
		return NULL;
#endif
}

void dm_multod(mat_double *m1, mat_double *m2, double d)
{
	if ( (m2->row != m1->row) || (m2->col != m1->col) ) {
		fprintf(stderr,"Error: dm_multod - dimension mismatch (%d,%d) and (%d,%d)\n",m1->row,m1->col,m2->row,m2->col);
		exit(1);
	}
	assert(m1->basis == m2->basis);

	switch (m1->type + m2->type) {
	case MAT_DENSE+MAT_DENSE :
		//DEBUGPRINT("dm_multod: F-F\n");
		cblas_daxpy(m1->row*m1->col,d,m2->data,1,m1->data,1);
		break;
	case MAT_DENSE_DIAG+MAT_DENSE_DIAG :
		//DEBUGPRINT("dm_multod: D-D\n");
		cblas_daxpy(m1->row,d,m2->data,1,m1->data,1);
		break;
	case MAT_DENSE+MAT_DENSE_DIAG : { /* dense, diagonal + full */
		mat_double *res, *dadd;
		double *ds1, *ds2, dd;
		int i;
		if (m1->type == MAT_DENSE) {
			//DEBUGPRINT("dm_multod: F-D\n");
			res = m1;
			dadd = m2;
			dd = d;
		} else {
			//DEBUGPRINT("dm_multod: D-F\n");
			res = dm_dup(m2);
			dm_muld(res,d);
			dadd = m1;
			dd = 1.0;
		}
		ds1 = res->data;
		ds2 = dadd->data;
		for (i=0; i<m1->row; i++) {
		      *ds1 += (*ds2)*dd;
		      ds1 += (m1->row+1);
		      ds2++;
		}
		//DEBUGPRINT("m1=%p, m2=%p, res=%p, dadd=%p\n",m1,m2,res,dadd);
		if (res != m1) {
			dm_swap_innards_and_destroy(m1,res);
		}
		break; }
	case MAT_DENSE+MAT_SPARSE :
	case MAT_DENSE+MAT_SPARSE_DIAG : { /* sparse plus dense full */
		mat_double *res, *dadd;
		double *ds, dd;
		int i, r, c, n;
		MKL_INT *ic;
		if (m1->type == MAT_DENSE) {
			res = m1;
			dadd = m2;
			dd = d;
		} else {
			res = dm_dup(m2);
			dm_muld(res,d);
			dadd = m1;
			dd = 1.0;
		}
		ic = dadd->icol;
		ds = dadd->data;
		for (r=0; r<dadd->row; r++) {
			n = dadd->irow[r+1] - dadd->irow[r];
			for (i=0; i<n; i++) {
				c = (*ic) - 1;
				ic++;
				res->data[r+c*res->row] += dd*(*ds);
				ds++;
			}
		}
		if (res != m1) {
			dm_swap_innards_and_destroy(m1,res);
		}
		break; }
	case MAT_SPARSE+MAT_SPARSE:
	case MAT_SPARSE+MAT_SPARSE_DIAG:
	case MAT_SPARSE_DIAG+MAT_SPARSE_DIAG: { /* both sparse */
		//DEBUGPRINT("  multo TYPES 1 --- Iz = %d, chan Iz = %d\n",m2->type,m1->type);
		mat_double *res;
		//dm_print(m1,"dm_multo m1");
		//printf("dm_multo d = %f\n",d);
		//dm_print(m2,"dm_multo m2");
		res = simpson_dcsradd(m1,d,m2);
		//DEBUGPRINT("  multo TYPES 2 --- Iz = %d, chan Iz = %d, res = %d\n",m2->type,m1->type,res->type);
		dm_swap_innards_and_destroy(m1,res);
		//DEBUGPRINT("  multo TYPES 3 --- Iz = %d, chan Iz = %d\n",m2->type,m1->type);
		break; }
	case MAT_DENSE_DIAG+MAT_SPARSE :
	case MAT_DENSE_DIAG+MAT_SPARSE_DIAG : {
		mat_double *res, *dum;
		if (m1->type == MAT_DENSE_DIAG) {
			dum = m2;
			dm_sparse(m1,SPARSE_TOL);
		} else {
			dum = dm_dup(m2);
			dm_sparse(dum,SPARSE_TOL);
		}
		res = simpson_dcsradd(m1,d,dum);
		if (dum != m2) free_double_matrix(dum);
		dm_swap_innards_and_destroy(m1,res);
		break; }
	default :
		fprintf(stderr,"Error: dm_multod - unknown matrix types '%d', '%d'\n",m1->type,m2->type);
		exit(1);
	}
	//dm_print(m1,"dm_multo result");
	if (m1->type == MAT_SPARSE || m1->type == MAT_SPARSE_DIAG) dm_sparse(m1,SPARSE_TOL);
}

/* to a double matrix, add real part of z*mat_complx */
void dm_multocc(mat_double *res, mat_complx *m, complx z)
{
	if (res->row != m->row || res->col != m->col) {
		fprintf(stderr,"Error: dm_multocc - dimension mismatch (%d,%d) + (%d,%d)\n",res->row,res->col, m->row, m->col);
		exit(1);
	}
	assert(res->basis == m->basis);

	switch (100*res->type + m->type) {
	case 100*MAT_DENSE+MAT_DENSE : {
		double *dv = (double*)(m->data);
		int len = res->row*res->col;
		cblas_daxpy(len,z.re,dv,2,res->data,1);
		cblas_daxpy(len,-z.im,dv+1,2,res->data,1);
		break; }
	case 100*MAT_DENSE_DIAG+MAT_DENSE_DIAG : {
		double *dv = (double*)(m->data);
		int len = res->row;
		cblas_daxpy(len,z.re,dv,2,res->data,1);
		cblas_daxpy(len,-z.im,dv+1,2,res->data,1);
		break; }
	case 100*MAT_DENSE+MAT_DENSE_DIAG : { /* dense, diagonal + full */
		double *ds1 = res->data;
		complx *ds2 = m->data;
		int i, step = res->row + 1;
		for (i=0; i<m->row; i++) {
			*ds1 += ds2->re*z.re-ds2->im*z.im;
			ds1 += step;
			ds2++;
		}
		break;}
	case 100*MAT_DENSE_DIAG+MAT_DENSE : {
		int len = m->row*m->col, step = m->row+1;
		complx *zs = m->data;
		double *dum, *ds, *dstop, *dss;
		dum = ds = (double*)malloc(len*sizeof(double));
		dstop = ds + len;
		do {
			*ds = zs->re*z.re-zs->im*z.im;
			ds++;
			zs++;
		} while (ds < dstop);
		ds = dum;
		dum = dss = res->data;
		res->data = ds;
		res->type = MAT_DENSE;
		do {
			*ds += *dss;
			ds += step;
			dss++;
		} while (ds < dstop);
		free(dum);
		break; }
	case 100*MAT_DENSE+MAT_SPARSE :
	case 100*MAT_DENSE+MAT_SPARSE_DIAG : { /* sparse plus dense full */
		MKL_INT *ic = m->icol;
		complx *zs = m->data;
		int i, r, c, n;
		for (r=0; r<m->row; r++) {
			n = m->irow[r+1] - m->irow[r];
			for (i=0; i<n; i++) {
				c = (*ic) - 1;
				ic++;
				res->data[r+c*m->row] += z.re*zs->re-z.im*zs->im;
				zs++;
			}
		}
		break;}
	case 100*MAT_SPARSE+MAT_DENSE:
	case 100*MAT_SPARSE_DIAG+MAT_DENSE: {
		int len = m->row*m->col;
		complx *zs = m->data;
		double *dum, *ds, *dstop, *dss;
		dum = ds = (double*)malloc(len*sizeof(double));
		dstop = ds + len;
		do {
			*ds = zs->re*z.re-zs->im*z.im;
			ds++;
			zs++;
		} while (ds < dstop);
		int i, r, c, n;
		MKL_INT *ic = res->icol;
		ds = dum;
		dum = dss = res->data;
		res->data = ds;
		for (r=0; r<res->row; r++) {
			n = res->irow[r+1] - res->irow[r];
			for (i=0; i<n; i++) {
				c = (*ic) - 1;
				ic++;
				res->data[r+c*res->row] += *dss;
				dss++;
			}
		}
		free(dum);
		free(res->irow); res->irow = NULL;
		free(res->icol); res->icol = NULL;
		res->type = MAT_DENSE;
		break; }
	case 100*MAT_SPARSE+MAT_SPARSE:
	case 100*MAT_SPARSE+MAT_SPARSE_DIAG:
	case 100*MAT_SPARSE_DIAG + MAT_SPARSE:
	case 100*MAT_SPARSE_DIAG+MAT_SPARSE_DIAG: { /* both sparse */
		mat_double *tmp, *dum;
		int len = m->irow[m->row]-1;
		dum = (mat_double*)malloc(sizeof(mat_double));
		dum->row = m->row; dum->col = m->col;
		dum->irow = m->irow;
		dum->icol = m->icol;
		dum->type = m->type;
		dum->basis = m->basis;
		dum->data = (double*)malloc(len*sizeof(double));
		double *ds = dum->data, *dstop = dum->data + len;
		complx *zs = m->data;
		do {
			*ds = zs->re*z.re-zs->im*z.im;
			ds++;
			zs++;
		} while (ds<dstop);
		tmp = simpson_dcsradd(res,1,dum);
		free(dum->data);
		free(dum);
		dm_swap_innards_and_destroy(res,tmp);
		break; }
	case 100*MAT_DENSE_DIAG+MAT_SPARSE :
	case 100*MAT_DENSE_DIAG+MAT_SPARSE_DIAG : {
		mat_double *tmp, *dum;
		int len = m->irow[m->row]-1;
		dum = (mat_double*)malloc(sizeof(mat_double));
		dum->row = m->row; dum->col = m->col;
		dum->irow = m->irow;
		dum->icol = m->icol;
		dum->type = m->type;
		dum->basis = m->basis;
		dum->data = (double*)malloc(len*sizeof(double));
		double *ds = dum->data, *dstop = dum->data + len;
		complx *zs = m->data;
		do {
			*ds = zs->re*z.re-zs->im*z.im;
			ds++;
			zs++;
		} while (ds<dstop);
		dm_sparse(res,SPARSE_TOL);
		tmp = simpson_dcsradd(res,1,dum);
		free(dum->data);
		free(dum);
		dm_swap_innards_and_destroy(res,tmp);
		break; }
	case 100*MAT_SPARSE + MAT_DENSE_DIAG:
	case 100*MAT_SPARSE_DIAG + MAT_DENSE_DIAG: {
		mat_double *dum, *tmp;
		int i, len = m->row;
		dum = double_matrix(len,len,MAT_SPARSE_DIAG,len,m->basis);
		complx *zs = m->data;
		for (i=0; i<len; i++){
			dum->data[i] = zs->re*z.re-zs->im*z.im;
			zs++;
			dum->icol[i] = dum->irow[i] = i+1;
		}
		dum->irow[len] = len + 1;
		tmp = simpson_dcsradd(res,1,dum);
		free_double_matrix(dum);
		dm_swap_innards_and_destroy(res,tmp);
		break;}
	default :
		fprintf(stderr,"Error: dm_multocc - unknown matrix types '%d', '%d'\n",res->type,m->type);
		exit(1);
	}
	if (res->type == MAT_SPARSE || res->type == MAT_SPARSE_DIAG) dm_sparse(res,SPARSE_TOL);

}

void dm_addto(mat_double *m1, mat_double *m2)
{
	//DEBUGPRINT(" addto TYPES 1 --- Iz = %d, chan Iz = %d\n",m2->type,m1->type);
	dm_multod(m1, m2, 1.0);
	//DEBUGPRINT(" addto TYPES 2 --- Iz = %d, chan Iz = %d\n",m2->type,m1->type);
}

void dm_addtodiag(mat_double *m, double d)
{
	if (m->row != m->col) {
		fprintf(stderr,"dm_addtodiag error: matrix is not square (%d,%d)\n",m->row,m->col);
		exit(1);
	}

	switch (m->type) {
		case MAT_DENSE : {
			int i, N = m->row;
			double *dd = m->data;
			for (i=0; i<N; i++) {
				*dd += d;
				dd += (N+1);
			}
			break; }
		case MAT_DENSE_DIAG : {
			int i, N = m->row;
			double *dd = m->data;
			for (i=0; i<N; i++) {
				*dd += d;
				dd++;
			}
			break; }
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG : {
			mat_double *dum, *res;
			dum = dm_creatediag('S',m->row,d,m->basis);
			res = simpson_dcsradd(m,1.0,dum);
			free_double_matrix(dum);
			dm_swap_innards_and_destroy(m,res);
			dm_sparse(m,SPARSE_TOL);
			break; }
		default :
		   fprintf(stderr,"Error: dm_addtodiag - unknown type '%d'\n",m->type);
		   exit(1);
		}
}

void dm_shrinktodiag(mat_double *m)
{
	int i, dim;

	if (m->row != m->col) {
		fprintf(stderr,"Error: dm_shrinktodiag - not a square matrix (%d,%d)\n",m->row,m->col);
		exit(1);
	}
	dim = m->row;

	switch (m->type) {
		case MAT_DENSE_DIAG :
		case MAT_SPARSE_DIAG :
			break;
		case MAT_DENSE : {
			double *dg, *z;
			z = dg = m->data;
			for (i=1; i<dim; i++) {
				dg++;
				z += (dim+1);
				*dg = *z;
			}
			m->data = (double*)realloc(m->data,dim*sizeof(double));
			assert(m->data != NULL);
			m->type = MAT_DENSE_DIAG;
			break; }
		case MAT_SPARSE : {
			mat_double *dum;
			dum = double_matrix(dim, dim, MAT_SPARSE_DIAG,dim,m->basis);
			for (i=0; i<dim; i++) {
				dum->data[i] = dm_getelem(m,i+1,i+1);
				dum->icol[i] = dum->irow[i] = i+1;
			}
			dum->irow[dim] = dim+1;
			dm_sparse(dum,SPARSE_TOL);
			dm_swap_innards_and_destroy(m,dum);
			break; }
	   default :
		   fprintf(stderr,"Error: dm_shrinktodiag - unknown type '%d'\n",m->type);
		   exit(1);
		}
}


/* return matrix element, 1-based indexing */
double dm_getelem(mat_double *m, int r, int c)
{
	double res;

	if ( (r>m->row) || (c>m->col) ) {
		fprintf(stderr,"Error: dm_getelem - index exceeds matrix dimension\n");
		exit(1);
	}

	switch (m->type) {
		case MAT_DENSE :
			res = m->data[(r-1) + (c-1)*m->row];
			break;
		case MAT_DENSE_DIAG :
			if (r != c)
				res = 0.0;
			else
				res = m->data[r-1];
			break;
		case MAT_SPARSE : {
			int i;
			res = 0.0;
			for (i=m->irow[r-1]; i<m->irow[r]; i++) {
				if (m->icol[i-1] == c) {
					res = m->data[i-1];
					break;
				}
			}
			break; }
		case MAT_SPARSE_DIAG :
			if (r != c) {
				res = 0.0;
			} else {
				int cc = m->irow[r-1];
				if (m->irow[r]-cc != 0)
					res = m->data[cc-1];
				else
					res = 0.0;
			}
			break;
	   default :
		   fprintf(stderr,"Error: dm_getelem - unknown type '%d'\n",m->type);
		   exit(1);
		}
		return res;
}

/* sparse res = a * b  , both sparse */
mat_double * simpson_dcsrmultcsr(mat_double *a, mat_double *b)
{
	assert(a->basis == b->basis);
#ifdef INTEL_MKL
		char trans='N';
		int sort=0;
		int request=1, info=0, nnzmax=a->row*a->col;
		int type;
		mat_double *res;

		if ( (a->type == MAT_SPARSE) || (b->type == MAT_SPARSE) ) {
			type = MAT_SPARSE;
		} else {
			type = MAT_SPARSE_DIAG;
		}
		res = double_matrix(a->row,b->col,type,0,a->basis);
		mkl_dcsrmultcsr(&trans, &request, &sort, &(a->row), &(a->col), &(b->col), a->data, a->icol, a->irow, b->data, b->icol, b->irow, res->data, res->icol, res->irow, &nnzmax, &info);
		dm_change_nnz(res,res->irow[res->row]-1);
		request = 2;
		mkl_dcsrmultcsr(&trans, &request, &sort, &(a->row), &(a->col), &(b->col), a->data, a->icol, a->irow, b->data, b->icol, b->irow, res->data, res->icol, res->irow, &nnzmax, &info);
		if (info != 0) {
			fprintf(stderr,"Error: dcsrmultcsr failed with the code '%d'\n",info);
			exit(1);
		}
		return res;
#else
		fprintf(stderr,"Error: sparse matrix routine dcsrmultcsr not compiled\n");
		exit(1);
		return NULL;
#endif
}

/* dense res = sparse a * dense b  */
mat_double * simpson_dcsrmm(mat_double *a, mat_double *b)
{
	assert(a->basis == b->basis);
#ifdef INTEL_MKL
		char trans='N', matdescra[]="GNNFUU";
		double done=1.0, dzero=0.0;
		mat_double *res;

		//DEBUGPRINT("Warning!!! Please check sparse * dense multiplication result!\n");
		// checked 8.9.2010 Zd.
		res = double_matrix(a->row,b->col,MAT_DENSE,0,a->basis);
		//dm_zero(res);
		mkl_dcsrmm(&trans, &(a->row), &(res->col), &(a->col), &done, matdescra, a->data, a->icol, a->irow, &(a->irow[1]), b->data, &(b->row), &dzero, res->data, &(res->row));
		return res;
#else
		fprintf(stderr,"Error: sparse matrix routine dcsrmmm not compiled\n");
		exit(1);
		return NULL;
#endif
}

/* dense res =alpha * sparse a * dense b + beta *dense c */
void simpson_dcsrmm2(double alpha, mat_double *a, mat_double *b, double beta, mat_double *c)
{
	assert(a->basis == b->basis && c->basis == a->basis);
#ifdef INTEL_MKL
		char trans='N', matdescra[]="GNNFUU";
		mkl_dcsrmm(&trans, &(a->row), &(c->col), &(a->col), &alpha, matdescra, a->data, a->icol, a->irow, &(a->irow[1]), b->data, &(b->row), &beta, c->data, &(c->row));
#else
		fprintf(stderr,"Error: sparse matrix routine dcsrmm2 not compiled\n");
		exit(1);
		return;
#endif
}

mat_double * dm_mul(mat_double *m1, mat_double *m2)
{
	mat_double *res;
	int i;

	if (m1->col != m2->row) {
		fprintf(stderr,"Error: dm_mul - dimension mismatch (%d,%d)*(%d,%d)\n",m1->row,m1->col,m2->row,m2->col);
		exit(1);
	}
	assert(m1->basis == m2->basis);

	switch (m1->type*100+m2->type) {
		case 100*MAT_DENSE+MAT_DENSE :
			res = double_matrix(m1->row,m2->col,MAT_DENSE,0,m1->basis);
			cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m1->row,res->col,m1->col,1.0,m1->data,m1->row,m2->data,m2->row,0.0,res->data,res->row);
			break;
		case 100*MAT_DENSE+MAT_DENSE_DIAG :
			assert(m2->row == m2->col);
			res = dm_dup(m1);
			for (i=0; i<m1->col; i++) {
				cblas_dscal(res->row,m2->data[i],&(res->data[i*res->row]),1);
			}
			break;
		case 100*MAT_DENSE+MAT_SPARSE : {
			DEBUGPRINT("warning from dm_mul: converting sparse to dense\n");
			mat_double *dum;
			dum = dm_dup(m2);
			dm_dense(dum);
			res = double_matrix(m1->row,m2->col,MAT_DENSE,0,m1->basis);
			cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m1->row,res->col,m1->col,1.0,m1->data,m1->row,dum->data,dum->row,0.0,res->data,res->row);
			free_double_matrix(dum);
			break; }
		case 100*MAT_DENSE+MAT_SPARSE_DIAG : {
			double d;
			res = dm_dup(m1);
			assert(m2->row == m2->col);
			for (i=0; i<res->col; i++) {
				d = dm_getelem(m2,i+1,i+1);
				cblas_dscal(res->row,d,&(res->data[i*res->row]),1);
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_DENSE :
			assert(m1->row == m1->col);
			res = dm_dup(m2);
			for (i=0; i<res->row; i++) {
				cblas_dscal(res->col,m1->data[i],&(res->data[i]),res->row);
			}
			break;
		case 100*MAT_DENSE_DIAG+MAT_DENSE_DIAG : {
			assert(m1->row == m1->col && m2->row == m2->col);
			res = dm_dup(m1);
			double *id1 = res->data, *id2 = m2->data;
			for (i=0; i<res->row; i++) {
				*id1 *= *id2;
				id1++;
				id2++;
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_SPARSE :
			dm_sparse(m1,SPARSE_TOL);
			res = simpson_dcsrmultcsr(m1,m2);
			break;
		case 100*MAT_DENSE_DIAG+MAT_SPARSE_DIAG :
			assert(m1->row == m1->col && m2->row == m2->col);
			res = dm_dup(m1);
			for (i=0; i<m1->row; i++) {
				res->data[i] *= dm_getelem(m2,i+1,i+1);
			}
			break;
		case 100*MAT_SPARSE+MAT_DENSE :
			res = simpson_dcsrmm(m1,m2);
			break;
		case 100*MAT_SPARSE+MAT_DENSE_DIAG : {
			mat_double *dum;
			dum = dm_dup(m2);
			dm_sparse(dum, SPARSE_TOL);
			res = simpson_dcsrmultcsr(m1,dum);
			free_double_matrix(dum);
			break; }
		case 100*MAT_SPARSE+MAT_SPARSE :
		case 100*MAT_SPARSE+MAT_SPARSE_DIAG :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE_DIAG :
			res = simpson_dcsrmultcsr(m1,m2);
			break;
		case 100*MAT_SPARSE_DIAG+MAT_DENSE : {
			assert(m1->row == m1->col);
			double d;
			res = dm_dup(m2);
			for (i=0; i<res->row; i++) {
				d = dm_getelem(m1,i+1,i+1);
				cblas_dscal(res->col,d,&(res->data[i]),res->row);
			}
			break; }
		case 100*MAT_SPARSE_DIAG+MAT_DENSE_DIAG :
			assert(m1->row == m1->col && m2->row == m2->col);
			res = dm_dup(m1);
			for (i=0; i<m1->irow[m1->row]-1; i++) {
				res->data[i] *= m2->data[res->icol[i]-1];
			}
			break;
	   default :
		   fprintf(stderr,"Error: dm_mul - unknown types '%d' and '%d'\n",m1->type, m2->type);
		   exit(1);
	}
	if (m1->type == MAT_SPARSE || m1->type == MAT_SPARSE_DIAG) dm_sparse(res,SPARSE_TOL);

	return res;
}

void dm_multo(mat_double *m1, mat_double *m2)
{
	mat_double *res;
	int i;

	if (m1->col != m2->row) {
		fprintf(stderr,"Error: dm_multo - dimension mismatch (%d,%d)*(%d,%d)\n",m1->row,m1->col,m2->row,m2->col);
		exit(1);
	}
	assert(m1->basis == m2->basis);

	switch (m1->type*100+m2->type) {
		case 100*MAT_DENSE+MAT_DENSE :
			res = double_matrix(m1->row,m2->col,MAT_DENSE,0,m1->basis);
			cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m1->row,res->col,m1->col,1.0,m1->data,m1->row,m2->data,m2->row,0.0,res->data,res->row);
			dm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_DENSE+MAT_DENSE_DIAG :
			assert(m2->row == m2->col);
			for (i=0; i<m1->col; i++) {
				cblas_dscal(m1->row,m2->data[i],&(m1->data[i*m1->row]),1);
			}
			break;
		case 100*MAT_DENSE+MAT_SPARSE : {
			DEBUGPRINT("warning from dm_multo: converting sparse to dense\n");
			mat_double *dum;
			dum = dm_dup(m2);
			dm_dense(dum);
			res = double_matrix(m1->row,m2->col,MAT_DENSE,0,m1->basis);
			cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m1->row,res->col,m1->col,1.0,m1->data,m1->row,dum->data,dum->row,0.0,res->data,res->row);
			free_double_matrix(dum);
			dm_swap_innards_and_destroy(m1,res);
			break; }
		case 100*MAT_DENSE+MAT_SPARSE_DIAG : {
			double d;
			assert(m2->row == m2->col);
			for (i=0; i<m1->col; i++) {
				d = dm_getelem(m2,i+1,i+1);
				cblas_dscal(m1->row,d,&(m1->data[i*m1->row]),1);
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_DENSE :
			assert(m1->row == m1->col);
			res = dm_dup(m2);
			for (i=0; i<res->row; i++) {
				cblas_dscal(res->col,m1->data[i],&(res->data[i]),res->row);
			}
			dm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_DENSE_DIAG+MAT_DENSE_DIAG : {
			assert(m1->row == m1->col && m2->row == m2->col);
			double *id1 = m1->data, *id2 = m2->data;
			for (i=0; i<m1->row; i++) {
				*id1 *= *id2;
				id1++;
				id2++;
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_SPARSE :
			dm_sparse(m1,SPARSE_TOL);
			res = simpson_dcsrmultcsr(m1,m2);
			dm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_DENSE_DIAG+MAT_SPARSE_DIAG :
			assert(m1->row == m1->col && m2->row == m2->col);
			for (i=0; i<m1->row; i++) {
				m1->data[i] *= dm_getelem(m2,i+1,i+1);
			}
			break;
		case 100*MAT_SPARSE+MAT_DENSE :
			res = simpson_dcsrmm(m1,m2);
			dm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_SPARSE+MAT_DENSE_DIAG : {
			mat_double *dum;
			dum = dm_dup(m2);
			dm_sparse(dum, SPARSE_TOL);
			res = simpson_dcsrmultcsr(m1,dum);
			free_double_matrix(dum);
			dm_swap_innards_and_destroy(m1,res);
			break; }
		case 100*MAT_SPARSE+MAT_SPARSE :
		case 100*MAT_SPARSE+MAT_SPARSE_DIAG :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE_DIAG :
			res = simpson_dcsrmultcsr(m1,m2);
			dm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_SPARSE_DIAG+MAT_DENSE : {
			assert(m1->row == m1->col);
			double d;
			res = dm_dup(m2);
			for (i=0; i<res->row; i++) {
				d = dm_getelem(m1,i+1,i+1);
				cblas_dscal(res->col,d,&(res->data[i]),res->row);
			}
			dm_swap_innards_and_destroy(m1,res);
			break; }
		case 100*MAT_SPARSE_DIAG+MAT_DENSE_DIAG :
			assert(m1->row == m1->col && m2->row == m2->col);
			for (i=0; i<m1->irow[m1->row]-1; i++) {
				m1->data[i] *= m2->data[m1->icol[i]-1];
			}
			break;
	   default :
		   fprintf(stderr,"Error: dm_multo - unknown types '%d' and '%d'\n",m1->type, m2->type);
		   exit(1);
	}
	if (m1->type == MAT_SPARSE || m1->type == MAT_SPARSE_DIAG) dm_sparse(m1,SPARSE_TOL);
}

void dm_multo_rev(mat_double *m1, mat_double *m2)
{
	mat_double *res;
	int i;

	if (m2->col != m1->row) {
		fprintf(stderr,"Error: dm_multo_rev - dimension mismatch (%d,%d)*(%d,%d)\n",m2->row,m2->col,m1->row,m1->col);
		exit(1);
	}
	assert(m1->basis == m2->basis);

	switch (m1->type*100+m2->type) {
		case 100*MAT_DENSE+MAT_DENSE :
			res = double_matrix(m2->row,m1->col,MAT_DENSE,0,m2->basis);
			cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m2->row,res->col,m2->col,1.0,m2->data,m2->row,m1->data,m1->row,0.0,res->data,res->row);
			dm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_DENSE+MAT_DENSE_DIAG :
			assert(m2->row == m2->col);
			for (i=0; i<m1->row; i++) {
				cblas_dscal(m1->col,m2->data[i],&(m1->data[i]),m1->row);
			}
			break;
		case 100*MAT_DENSE+MAT_SPARSE :
			res = simpson_dcsrmm(m2,m1);
			dm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_DENSE+MAT_SPARSE_DIAG : {
			double d;
			assert(m2->row == m2->col);
			for (i=0; i<m1->row; i++) {
				d = dm_getelem(m2,i+1,i+1);
				cblas_dscal(m1->col,d,&(m1->data[i]),m1->row);
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_DENSE :
			assert(m1->row == m1->col);
			res = dm_dup(m2);
			for (i=0; i<res->col; i++) {
				cblas_dscal(res->row,m1->data[i],&(res->data[i*res->row]),1);
			}
			dm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_DENSE_DIAG+MAT_DENSE_DIAG : {
			assert(m1->row == m1->col && m2->row == m2->col);
			double *id1 = m1->data, *id2 = m2->data;
			for (i=0; i<m1->row; i++) {
				*id1 *= *id2;
				id1++;
				id2++;
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_SPARSE :
			dm_sparse(m1,SPARSE_TOL);
			res = simpson_dcsrmultcsr(m2,m1);
			dm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_DENSE_DIAG+MAT_SPARSE_DIAG :
			assert(m1->row == m1->col && m2->row == m2->col);
			for (i=0; i<m1->row; i++) {
				m1->data[i] *= dm_getelem(m2,i+1,i+1);
			}
			break;
		case 100*MAT_SPARSE+MAT_DENSE :
			DEBUGPRINT("warning from dm_multo_rev: converting sparse to dense\n");
			dm_dense(m1);
			res = double_matrix(m2->row,m1->col,MAT_DENSE,0,m2->basis);
			cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m2->row,res->col,m2->col,1.0,m2->data,m2->row,m1->data,m1->row,0.0,res->data,res->row);
			dm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_SPARSE+MAT_DENSE_DIAG : {
			mat_double *dum;
			dum = dm_dup(m2);
			dm_sparse(dum, SPARSE_TOL);
			res = simpson_dcsrmultcsr(dum,m1);
			free_double_matrix(dum);
			dm_swap_innards_and_destroy(m1,res);
			break; }
		case 100*MAT_SPARSE+MAT_SPARSE :
		case 100*MAT_SPARSE+MAT_SPARSE_DIAG :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE_DIAG :
			res = simpson_dcsrmultcsr(m2,m1);
			dm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_SPARSE_DIAG+MAT_DENSE : {
			assert(m1->row == m1->col);
			double d;
			res = dm_dup(m2);
			for (i=0; i<res->col; i++) {
				d = dm_getelem(m1,i+1,i+1);
				cblas_dscal(res->row,d,&(res->data[i*res->row]),1);
			}
			dm_swap_innards_and_destroy(m1,res);
			break; }
		case 100*MAT_SPARSE_DIAG+MAT_DENSE_DIAG :
			assert(m1->row == m1->col && m2->row == m2->col);
			for (i=0; i<m1->irow[m1->row]-1; i++) {
				m1->data[i] *= m2->data[m1->icol[i]-1];
			}
			break;
	   default :
		   fprintf(stderr,"Error: dm_multo_rev - unknown types '%d' and '%d'\n",m1->type, m2->type);
		   exit(1);
	}
	if (m1->type == MAT_SPARSE || m1->type == MAT_SPARSE_DIAG) dm_sparse(m1,SPARSE_TOL);
}

/****
 *  C = alpha*A*B + beta*C
 ****/
void dm_mm(double alpha, mat_double *A, mat_double *B, double beta, mat_double *C)
{

	long int key = A->type*10000+B->type*100+C->type;

	//printf("\t\t\tdm_mm types: %d, %d, %d\n",A->type,B->type,C->type);

	if (key == 10000*MAT_DENSE+100*MAT_DENSE+MAT_DENSE) {
		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,A->row,C->col,A->col,alpha,A->data,A->row,B->data,B->row,beta,C->data,C->row);
	} else if (key == 10000*MAT_SPARSE+100*MAT_DENSE+MAT_DENSE) {
		simpson_dcsrmm2(alpha,A,B,beta,C);
	} else if (key == 10000*MAT_SPARSE_DIAG+100*MAT_DENSE+MAT_DENSE) {
		simpson_dcsrmm2(alpha,A,B,beta,C);
	} else {
		mat_double *dum = dm_mul(A,B);
		dm_muld(C,beta);
		dm_multod(C,dum,alpha);
		free_double_matrix(dum);
	}
}

void dm_print(mat_double *m, char * title)
{
	int i, j, Nr = m->row, Nc = m->col;

	switch (m->type) {
		case MAT_DENSE :
			printf("%s : full dense matrix(%i,%i) in basis %d\n",title,Nr,Nc,m->basis);
			for (i=0; i<Nr; i++) {
				printf(" ");
			    for (j=0; j<Nc; j++) {
			    	printf("%9.6g ",m->data[i+j*Nr]);
			    }
			    printf("\n");
			}
			break;
		case MAT_DENSE_DIAG :
			printf("%s : diagonal dense matrix(%i,%i) in basis %d\n ",title,Nr,Nc,m->basis);
			for (i=0; i<Nr; i++) {
				printf("%9.6g ",m->data[i]);
			}
			printf("\n");
			break;
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG : {
			int nn, *ic = m->icol, len = m->irow[Nr] - 1;
			double *z = m->data;
			printf("%s : %s matrix(%i,%i) in basis %d, %d elements\n",title,matrix_type(m->type),Nr,Nc,m->basis,len);
			for (i=0; i<Nr; i++) {
				nn = m->irow[i+1] - m->irow[i];
				for (j=0; j<nn; j++) {
					printf(" [%d, %d] = %9.6g\n",i+1,*ic,*z);
					ic++;
					z++;
				}
			}
			break; }
	   default :
		   fprintf(stderr,"Error: dm_print - unknown type '%d'\n",m->type);
		   exit(1);
	}
}

/* norm estimation taken as largest sum of absolute values on a row */
/* i.e. infinity norm */
double dm_normest(mat_double *m)
{
	double res=-10.0, dum;
	int i,j;

	switch (m->type) {
	case MAT_DENSE_DIAG:
		for (i=0; i<m->row; i++) {
			dum = fabs(m->data[i]);
			if (dum > res) res = dum;
		}
		break;
	case MAT_SPARSE_DIAG:
		for (i=0; i<m->irow[m->row]-1; i++) {
			dum = fabs(m->data[i]);
			if (dum > res) res = dum;
		}
		break;
	case MAT_DENSE:
		for (i=0; i<m->row; i++) {
			dum = cblas_dasum(m->col,&(m->data[i]),m->row);
			if (dum > res) res = dum;
		}
		break;
	case MAT_SPARSE: {
		int N;
		double *dval = m->data;
		for (i=0; i<m->row; i++) {
			N = m->irow[i+1] - m->irow[i];
			dum = 0;
			for (j=0; j<N; j++ ) {
				dum += fabs(*dval);
				dval++;
			}
			if (dum > res) res = dum;
		}
		break; }
	default:
		fprintf(stderr,"Error: dm_normest - invalid matrix type\n");
		exit(1);
	}

	return res;
}

double dm_dm_normest(mat_double *re, mat_double *im)
{
	double res=-10.0, dum, *r1, *r2;
	int i,j;

	assert(re->type == im->type);
	switch (re->type) {
	case MAT_DENSE_DIAG:
		r1 = re->data;
		r2 = im->data;
		for (i=0; i<re->row; i++) {
			dum = sqrt((*r1)*(*r1)+(*r2)*(*r2));
			if (dum > res) res = dum;
			r1++;
			r2++;
		}
		break;
	case MAT_SPARSE_DIAG:
		fprintf(stderr,"Error: dm_dm_normest - MAT_SPARSE_DIAG\n");
		break;
	case MAT_DENSE:
		r1 = re->data;
		r2 = im->data;
		for (i=0; i<re->row; i++) {
			dum = 0;
			for (j=0; j<re->col; j++) {
				dum += sqrt((*r1)*(*r1)+(*r2)*(*r2));
				r1++;
				r2++;
			}
			if (dum > res) res = dum;
		}
		break;
	case MAT_SPARSE: {
		fprintf(stderr,"Error: dm_dm_normest - MAT_SPARSE\n");
		break; }
	default:
		fprintf(stderr,"Error: dm_dm_normest - invalid matrix type\n");
		exit(1);
	}

	return res;
}

double dm_max(mat_double *m)
{
	int len=0;

	switch (m->type) {
	case MAT_DENSE  : len = m->row * m->col; break;
	case MAT_DENSE_DIAG  : len = m->row; break;
	case MAT_SPARSE :
	case MAT_SPARSE_DIAG : len = m->irow[m->row]-1;
	}
	return m->data[cblas_idamax(len,m->data,1)];
}

int dm_nnz(mat_double *m)
{
	int nnz, i;
	double *d;

	switch (m->type) {
	case MAT_DENSE :
		nnz = 0;
		d = m->data;
		for (i=0; i<m->row*m->col; i++) {
			if (fabs(*d) > TINY) nnz++;
			d++;
		}
		break;
	case MAT_DENSE_DIAG:
		nnz = 0;
		d = m->data;
		for (i=0; i<m->row; i++) {
			if (fabs(*d) > TINY) nnz++;
			d++;
		}
		break;
	case MAT_SPARSE:
	case MAT_SPARSE_DIAG:
		nnz = m->irow[m->row] -1;
		break;
	default:
		fprintf(stderr,"Error: dm_nnz - invalid matrix type\n");
		exit(1);
	}

	return nnz;
}

void dm_nnz2one(mat_double *m)
{
	int i;
	double *d = m->data;

	switch (m->type) {
		case MAT_DENSE :
			for (i=0; i<m->row*m->col; i++) {
				if ( fabs(*d) > TINY )
					*d = 1;
				else
					*d = 0;
				d++;
			}
			break;
		case MAT_DENSE_DIAG :
			for (i=0; i<m->row; i++) {
				if ( fabs(*d) > TINY )
					*d = 1;
				else
					*d = 0;
				d++;
			}
			break;
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG : {
			int dim = m->irow[m->row] - 1;
			for (i=0; i<dim; i++) {
				if ( fabs(*d) > TINY )
					*d = 1;
				else
					*d = 0;
				d++;
			}
			dm_sparse(m,SPARSE_TOL);
			break;}
	   default :
		   fprintf(stderr,"Error: dm_nnz2one - unknown type '%d'\n",m->type);
		   exit(1);
	}
}

/********************************************************************
 *           C O M P L E X    M A T R I C E S                       *
 ********************************************************************/
mat_complx * cm_creatediag(const char c, int dim, complx z, int basis)
{
	mat_complx * res;
	int i;

	if ( (c == 'd') || (c == 'D') ) {
		res = complx_matrix(dim,dim,MAT_DENSE_DIAG,dim,basis);
		for (i=0; i<dim; i++) res->data[i] = z;
		return res;
	}
	if ( (c == 's') || (c == 'S') ) {
		res = complx_matrix(dim,dim,MAT_SPARSE_DIAG,dim,basis);
		for (i=0; i<dim; i++) {
			res->data[i] = z;
			res->icol[i] = res->irow[i] = i+1;
		}
		res->irow[dim] = dim+1;
		return res;
	}
    fprintf(stderr,"Error: cm_creatediag called with wrong code '%c'\n",c);
    exit(1);
    return NULL;
}

mat_complx *cm_dup(mat_complx *m)
{
	int len;
	mat_complx *mm;


	switch (m->type) {
	case MAT_DENSE :
		mm = complx_matrix(m->row,m->col,m->type,0,m->basis);
		memcpy(mm->data,m->data,(m->row*m->col)*sizeof(complx));
		break;
	case MAT_DENSE_DIAG :
		assert(m->row == m->col);
		mm = complx_matrix(m->row,m->col,m->type,0,m->basis);
		memcpy(mm->data,m->data,(m->row)*sizeof(complx));
		break;
	case MAT_SPARSE :
	case MAT_SPARSE_DIAG :
		len = m->irow[m->row] -1;
		mm = complx_matrix(m->row,m->col,m->type,len,m->basis);
		memcpy(mm->data,m->data,len*sizeof(complx));
		memcpy(mm->icol, m->icol, len*sizeof(MKL_INT));
		memcpy(mm->irow, m->irow, (m->row+1)*sizeof(MKL_INT));
		break;
   default :
	   fprintf(stderr,"Error: cm_dup - unknown type '%d'\n",m->type);
	   exit(1);
	}
	return mm;
}

void cm_copy(mat_complx *m1, mat_complx *m2)
{
	if (m1 == NULL) {
		fprintf(stderr,"Error: cm_copy - destination is NULL\n");
		exit(1);
	}
	if ( (m1->type != m2->type) || (m1->row != m2->row) || (m1->col != m2->col) ) {
		DEBUGPRINT("cm_copy: different matrix types or dimensions (%d,%d)/%d <- (%d,%d)/%d\n",m1->row,m1->col,m1->type,m2->row,m2->col,m2->type);
		mat_complx *dum = cm_dup(m2);
		cm_swap_innards_and_destroy(m1,dum);
		return;
	}
	switch (m1->type) {
		case MAT_DENSE :
			memcpy(m1->data,m2->data,m2->row*m2->col*sizeof(complx));
			break;
		case MAT_DENSE_DIAG :
			memcpy(m1->data,m2->data,m2->row*sizeof(complx));
			break;
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG : {
			int nnz = m2->irow[m2->row] - 1;
			cm_change_nnz(m1,nnz);
			memcpy(m1->data,m2->data,nnz*sizeof(complx));
			memcpy(m1->icol,m2->icol,nnz*sizeof(MKL_INT));
			memcpy(m1->irow,m2->irow,(m2->row+1)*sizeof(MKL_INT));
			break; }
	   default :
		   fprintf(stderr,"Error: cm_copy - unknown type '%d'\n",m1->type);
		   exit(1);
	}
	m1->basis = m2->basis;
}

void cm_unit(mat_complx *m)
{
	int i, dim;
	complx *z;

	if (m->row != m->col) {
		fprintf(stderr,"cm_unit error: not square matrix (%d,%d)\n",m->row,m->col);
		exit(1);
	}
	dim = m->row;

	/*
	switch (m->type) {
		case MAT_DENSE :
			free(m->data);
			m->data = (complx*)malloc(dim*sizeof(complx));
		case MAT_DENSE_DIAG :
			//for (i=0; i<dim; i++) m->data[i] = Complx(1.0,0.0);
			z = m->data;
			for (i=0; i<dim; i++) {
				z->re = 1.0;
				z->im = 0.0;
				z++;
			}
			m->type = MAT_DENSE_DIAG;
			break;
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG :
			cm_change_nnz(m,dim);
			z = m->data;
			for (i=0; i<dim; i++) {
				z->re = 1.0;
				z->im = 0.0;
				z++;
				m->irow[i] = m->icol[i] = i+1;
			}
			m->irow[dim] = dim+1;
			m->type = MAT_SPARSE_DIAG;
			break;
	   default :
		   fprintf(stderr,"Error: cm_unit - unknown type '%d'\n",m->type);
		   exit(1);
		}
		*/
	switch (m->type) {
		case MAT_DENSE :
			free(m->data);
			m->data = (complx*)malloc(dim*sizeof(complx));
			break;
		case MAT_DENSE_DIAG :
			break;
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG :
			free(m->data);
			free(m->irow); m->irow = NULL;
			free(m->icol); m->icol = NULL;
			m->data = (complx*)malloc(dim*sizeof(complx));
			break;
	   default :
		   fprintf(stderr,"Error: cm_unit - unknown type '%d'\n",m->type);
		   exit(1);
		}
	z = m->data;
	for (i=0; i<dim; i++) {
		z->re = 1.0;
		z->im = 0.0;
		z++;
	}
	m->type = MAT_DENSE_DIAG;
}

void cm_zero(mat_complx *m)
{
	switch (m->type) {
	case MAT_DENSE :
		memset(m->data, 0, (m->row * m->col) * sizeof(complx));
		break;
	case MAT_DENSE_DIAG :
		memset(m->data, 0, m->row * sizeof(complx));
		break;
	case MAT_SPARSE :
	case MAT_SPARSE_DIAG : {
		DEBUGPRINT("cm_zero sparse:");
		if (m->irow[m->row] != 0) {
			DEBUGPRINT(" reallocation and");
			complx * new_data = (complx*)realloc(m->data, sizeof(complx));
			MKL_INT * new_icol = (MKL_INT*)realloc(m->icol, sizeof(MKL_INT));
			assert( (new_data != NULL) && (new_icol != NULL) );
			m->data = new_data;
			m->icol = new_icol;
		}
		DEBUGPRINT(" filling\n");
		m->data[0] = Cnull;
		m->irow[0] = 1;
		int i; for (i=1; i<=m->row; i++) m->irow[i] = 2;
		m->icol[0] = 1;
		break; }
   default :
	   fprintf(stderr,"Error: cm_zero - unknown type '%d'\n",m->type);
	   exit(1);
	}
}

void cm_dense(mat_complx *m)
{
	int i, j, N, c, *ic;
	complx *ddata, *dd;

	switch (m->type) {
		case MAT_DENSE :
		case MAT_DENSE_DIAG :
			break;
		case MAT_SPARSE :
			ddata = (complx*)calloc(m->row*m->col,sizeof(complx));
			ic = m->icol;
			dd = m->data;
			for (i=0; i<m->row; i++) {
				N = m->irow[i+1] - m->irow[i];
				for (j=0;j<N; j++) {
					c = (*ic) -1;
					ic++;
					ddata[i+c*m->row] = *dd;
					dd++;
				}
			}
			free(m->icol); m->icol = NULL;
			free(m->irow); m->irow = NULL;
			free(m->data); m->data = ddata;
			m->type = MAT_DENSE;
			break;
		case MAT_SPARSE_DIAG :
			assert(m->row == m->col);
			ddata = (complx*)calloc(m->row,sizeof(complx));
			N = m->irow[m->row] - 1;
			ic = m->icol;
			dd = m->data;
			for (i=0; i<N; i++) {
				ddata[*ic -1] = *dd;
				ic++;
				dd++;
			}
			free(m->icol); m->icol = NULL;
			free(m->irow); m->irow = NULL;
			free(m->data); m->data = ddata;
			m->type = MAT_DENSE_DIAG;
			break;
	   default :
		   fprintf(stderr,"Error: cm_dense - unknown type '%d'\n",m->type);
		   exit(1);
		}
}

/* creates full dense matrix even from diagonal */
void cm_dense_full(mat_complx *m)
{
	int i, j, N, c, *ic;
	complx *ddata, *dd;

	switch (m->type) {
		case MAT_DENSE :
			break;
		case MAT_DENSE_DIAG :
			ddata = (complx*)calloc(m->row*m->col,sizeof(complx));
			cblas_zcopy(m->col,m->data,1,ddata,m->row+1);
			free(m->data); m->data = ddata;
			m->type = MAT_DENSE;
			break;
		case MAT_SPARSE :
			ddata = (complx*)calloc(m->row*m->col,sizeof(complx));
			ic = m->icol;
			dd = m->data;
			for (i=0; i<m->row; i++) {
				N = m->irow[i+1] - m->irow[i];
				for (j=0;j<N; j++) {
					c = (*ic) -1;
					ic++;
					ddata[i+c*m->row] = *dd;
					dd++;
				}
			}
			free(m->icol); m->icol = NULL;
			free(m->irow); m->irow = NULL;
			free(m->data); m->data = ddata;
			m->type = MAT_DENSE;
			break;
		case MAT_SPARSE_DIAG :
			assert(m->row == m->col);
			ddata = (complx*)calloc(m->row*m->col,sizeof(complx));
			N = m->irow[m->row] - 1;
			ic = m->icol;
			dd = m->data;
			for (i=0; i<N; i++) {
				ddata[(*ic -1)*m->row + *ic -1] = *dd;
				ic++;
				dd++;
			}
			free(m->icol); m->icol = NULL;
			free(m->irow); m->irow = NULL;
			free(m->data); m->data = ddata;
			m->type = MAT_DENSE;
			break;
	   default :
		   fprintf(stderr,"Error: cm_dense_full - unknown type '%d'\n",m->type);
		   exit(1);
		}
}

void cm_sparse(mat_complx *m, double tol)
{
	int realloc_size;
	uint64_t Nmax = (uint64_t)floor((uint64_t)m->row*m->col*(1.0-SPARSITY));
	complx *new_data;
	MKL_INT *new_icol;

	//DEBUGPRINT("cm_sparse will allow %d NNZ\n",Nmax);
	//cm_print(m,"cm_sparse IN");
	switch (m->type) {
		case MAT_DENSE : {
			int i, j, n=0, rc=1;
			complx *d1, *d2;
			MKL_INT *ir, *ic;

			m->irow = ir = (MKL_INT*)malloc((m->row+1)*sizeof(MKL_INT));
			m->icol = ic = (MKL_INT*)malloc(Nmax*sizeof(MKL_INT));
			new_data = d2 = (complx*)malloc(Nmax*sizeof(complx));
		    *ir = rc;
		    for (i=0; i<m->row; i++) {
		    	d1 = m->data + i;
			    for (j=0; j<m->col; j++) {
			    	if ( (fabs(d1->re) >= tol) || (fabs(d1->im) >= tol) ) {
			    		if (++n <= Nmax) {
			    			*d2 = *d1;
			    			d2++;
			    			*ic = j+1;
			    			ic++;
			    			rc++;
			    		}
			    	} else {
			    		d1->re = 0.0;
			    		d1->im = 0.0;
			    	}
			    	d1 += m->row;
		    	}
		    	ir++;
		    	*ir = rc;
		    }
		    if (n>Nmax) {
		    	DEBUGPRINT("cm_sparse: dense matrix has %d NNZ (%d allowed) - remains dense\n",n,Nmax);
		    	free(m->irow); m->irow = NULL;
		    	free(m->icol); m->icol = NULL;
		    	free(new_data);
		    	return;
		    }
		    free(m->data);
		    m->data = new_data;
		    m->type = MAT_SPARSE;
		    new_data = NULL;
			break; }
		case MAT_DENSE_DIAG : {
			int r, rr=1;
			complx *d1, *d2;
			MKL_INT *ic;
			m->type = MAT_SPARSE_DIAG;
			m->icol = ic = (MKL_INT*)malloc(m->row*sizeof(MKL_INT));
			m->irow = (MKL_INT*)malloc((m->row+1)*sizeof(MKL_INT));
			m->irow[0] = rr;
			d1 = d2 = m->data;
			for (r=1; r<=m->row; r++) {
				if ( (fabs(d2->re)<tol) && (fabs(d2->im)<tol) ) {
					d2++;
				} else {
					if (d1 != d2) {
						*d1 = *d2;
					}
					*ic = r;
					d1++;
					d2++;
					rr++;
					ic++;
				}
				m->irow[r] = rr;
			}
			break; }
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG : {
			complx *d1, *d2;
			MKL_INT *c1, *c2, *r1, *r2;
			int nc, ncnew, rr, ic, ir;

			r1 = m->irow;
			r2 = r1+1;
			rr = *r1;
			d1 = d2 = m->data;
			c1 = c2 = m->icol;
			for (ir=0; ir<m->row; ir++) {
				nc = ncnew = *r2 - rr;
				for (ic=0; ic<nc; ic++) {
					if ( (fabs(d2->re)<tol) && (fabs(d2->im)<tol) ) {
						d2++;
						c2++;
						ncnew--;
					} else {
						if (d1 != d2) {
							*d1 = *d2;
							*c1 = *c2;
						}
						d1++;
						d2++;
						c1++;
						c2++;
					}
				}
				rr = *r2;
				*r2 = *r1 + ncnew;
				r1++;
				r2++;
			}
			break; }
	   default :
		   fprintf(stderr,"Error: cm_sparse - unknown type '%d'\n",m->type);
		   exit(1);
	}
	realloc_size = m->irow[m->row]-1;
	assert(realloc_size >=0);
	if (realloc_size > Nmax) {
		DEBUGPRINT("cm_sparse: %d NNZ was detected but %d is allowed - converting to dense matrix\n",realloc_size,Nmax);
		cm_dense(m);
		return;
	}
	//DEBUGPRINT("cm_sparse restricted to %d NNZ\n",realloc_size);
	if (realloc_size == 0) {
			realloc_size = 1;
			m->data[0] = Cnull;
			m->icol[0] = m->col;
			m->irow[m->row] = 2;
	}
	new_data = (complx*)realloc(m->data, realloc_size*sizeof(complx));
	new_icol = (MKL_INT*)realloc(m->icol, realloc_size*sizeof(MKL_INT));
	assert( (new_data != NULL) && (new_icol != NULL) );
	m->data = new_data;
	m->icol = new_icol;
	/*cm_print(m,"cm_sparse OUT");*/
	/*exit(1); */
}

complx cm_getelem(mat_complx *m, int r, int c)
{
	complx res;

	if ( (r>m->row) || (c>m->col) ) {
		fprintf(stderr,"Error: cm_getelem - index exceeds matrix dimension (%d,%d)>(%d,%d)\n",r,c,m->row,m->col);
		exit(1);
	}

	switch (m->type) {
		case MAT_DENSE :
			res = m->data[(r-1) + (c-1)*m->row];
			break;
		case MAT_DENSE_DIAG :
			if (r != c)
				res = Cnull;
			else
				res = m->data[r-1];
			break;
		case MAT_SPARSE : {
			int i;
			res = Cnull;
			for (i=m->irow[r-1]; i<m->irow[r]; i++) {
				if (m->icol[i-1] == c) {
					res = m->data[i-1];
					break;
				}
			}
			break; }
		case MAT_SPARSE_DIAG :
			if (r != c) {
				res = Cnull;
			} else {
				int cc = m->irow[r-1];
				if (m->irow[r]-cc != 0)
					res = m->data[cc-1];
				else
					res = Cnull;
			}
			break;
	   default :
		   fprintf(stderr,"Error: cm_getelem - unknown type '%d'\n",m->type);
		   exit(1);
	}
	return res;
}

void cm_shrinktodiag(mat_complx *m)
{
	int i, dim;

	if (m->row != m->col) {
		fprintf(stderr,"Error: cm_shrinktodiag - not a square matrix (%d,%d)\n",m->row,m->col);
		exit(1);
	}
	dim = m->row;

	switch (m->type) {
		case MAT_DENSE_DIAG :
		case MAT_SPARSE_DIAG :
			break;
		case MAT_DENSE : {
			complx *cdata, *z;
			cdata = (complx*)malloc(dim*sizeof(complx));
			z = m->data;
			for (i=0; i<dim; i++) {
				cdata[i] = *z;
				z += (dim+1);
			}
			free(m->data);
			m->data = cdata;
			m->type = MAT_DENSE_DIAG;
			break; }
		case MAT_SPARSE : {
			mat_complx *dum;
			dum = complx_matrix(dim, dim, MAT_SPARSE_DIAG,dim,m->basis);
			dum->basis = m->basis;
			for (i=0; i<dim; i++) {
				dum->data[i] = cm_getelem(m,i+1,i+1);
				dum->icol[i] = dum->irow[i] = i+1;
			}
			dum->irow[dim] = dim+1;
			cm_sparse(dum,SPARSE_TOL);
			cm_swap_innards_and_destroy(m,dum);
			break; }
	   default :
		   fprintf(stderr,"Error: cm_shrinktodiag - unknown type '%d'\n",m->type);
		   exit(1);
		}
}

mat_double * cm_real(mat_complx *m)
{
	mat_double *res;
	int len;

	switch (m->type) {
		case MAT_DENSE :
			len = m->row*m->col;
			res = double_matrix(m->row, m->col, m->type, 0,m->basis);
			break;
		case MAT_DENSE_DIAG :
			len = m->row;
			res = double_matrix(len, len, m->type, 0, m->basis);
			break;
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG :
			len = m->irow[m->row] - 1;
			res = double_matrix(m->row, m->col, m->type, len, m->basis);
			memcpy(res->icol,m->icol,len*sizeof(MKL_INT));
			memcpy(res->irow,m->irow,(m->row+1)*sizeof(MKL_INT));
			break;
	   default :
		   fprintf(stderr,"Error: cm_real - unknown type '%d'\n",m->type);
		   exit(1);
	}
	cblas_dcopy(len,(double*)(m->data),2,res->data,1);
	return res;
}

/* will return just diagonal no matter what m is */
mat_double *cm_realdiag(mat_complx *m)
{
	mat_double *res;
	int i, dim;

	if (m->row != m->col) {
		fprintf(stderr,"Error: cm_realdiag - not a square matrix (%d,%d)\n",m->row,m->col);
		exit(1);
	}
	dim = m->row;

	switch (m->type) {
		case MAT_DENSE : {
			complx *z = m->data;
			res = double_matrix(dim, dim, MAT_DENSE_DIAG, 0, m->basis);
			for (i=0; i<dim; i++) {
				res->data[i] = z->re;
				z += (dim+1);
			}
			break; }
		case MAT_DENSE_DIAG :
			res = double_matrix(dim, dim, MAT_DENSE_DIAG, 0, m->basis);
			cblas_dcopy(dim,(double*)(m->data),2,res->data,1);
			break;
		case MAT_SPARSE : {
			complx z;
			res = double_matrix(dim, dim, MAT_SPARSE_DIAG, dim, m->basis);
			for (i=0; i<dim; i++) {
				z = cm_getelem(m,i+1,i+1);
				res->data[i] = z.re;
				res->icol[i] = res->irow[i] = i+1;
			}
			res->irow[dim] = dim+1;
			dm_sparse(res,SPARSE_TOL);
			break; }
		case MAT_SPARSE_DIAG : {
			int len = m->irow[dim] - 1;
			res = double_matrix(dim, dim, MAT_SPARSE_DIAG, dim, m->basis);
			cblas_dcopy(len,(double*)(m->data),2,res->data,1);
			memcpy(res->icol, m->icol, len*sizeof(MKL_INT));
			memcpy(res->irow, m->irow, (dim+1)*sizeof(MKL_INT));
			break; }
	   default :
		   fprintf(stderr,"Error: matrix_type - unknown type '%d'\n",m->type);
		   exit(1);
	}
	return res;
}

mat_double * cm_imag(mat_complx *m)
{
	mat_double *res;
	int len;

	switch (m->type) {
		case MAT_DENSE :
			len = m->row*m->col;
			res = double_matrix(m->row, m->col, m->type, 0, m->basis);
			break;
		case MAT_DENSE_DIAG :
			len = m->row;
			res = double_matrix(len,len,m->type,0, m->basis);
			break;
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG :
			len = m->irow[m->row] - 1;
			res = double_matrix(m->row, m->col, m->type, len, m->basis);
			memcpy(res->icol,m->icol,len*sizeof(MKL_INT));
			memcpy(res->irow,m->irow,(m->row+1)*sizeof(MKL_INT));
			break;
	   default :
		   fprintf(stderr,"Error: cm_real - unknown type '%d'\n",m->type);
		   exit(1);
	}
	cblas_dcopy(len,(double*)(m->data)+1,2,res->data,1);
	return res;
}


int cm_ishermit(mat_complx *m)
{
	int i, j, len, dim, res = 0;
	double d;

	if (m->row != m->col) {
		return 0;
	}
	dim = m->row;

	switch (m->type) {
		case MAT_DENSE_DIAG :
			d = cblas_dasum(dim,(double*)(m->data)+1,2);
			if (d < TINY) res = 1;
			break;
		case MAT_SPARSE_DIAG :
			len = m->irow[dim] - 1;
			d = cblas_dasum(len,(double*)(m->data)+1,2);
			if (d < TINY) res = 1;
			break;
		case MAT_DENSE : {
			complx v1, v2;
			d = cblas_dasum(dim,(double*)(m->data)+1,2*(dim+1));
			if (d > TINY) return 0;
			for (i=0; i<dim; i++) {
				for (j=i+1; j<dim; j++) {
			         v1 = m->data[i+dim*j];
			         v2 = m->data[j+dim*i];
			         if (fabs(v1.re - v2.re) > TINY) return 0;
			         if (fabs(v1.im + v2.im) > TINY) return 0;
			    }
			}
			res = 1;
			break; }
		case MAT_SPARSE : {
			int *ic = m->icol;;
			complx val, *z = m->data;
			for (i=0; i<dim; i++) {
				len = m->irow[i+1] - m->irow[i];
				for (j=0; j<len; j++) {
					if (*ic == i+1) {
						if (fabs(z->im) > TINY) return 0;
					} else {
						val = cm_getelem(m,*ic,i+1);
						if (fabs(val.re - z->re) > TINY) return 0;
						if (fabs(val.im + z->im) > TINY) return 0;
					}
					ic++;
					z++;
				}
			}
			res = 1;
			break; }
	   default :
		   fprintf(stderr,"Error: matrix_type - unknown type '%d'\n",m->type);
		   exit(1);
	}
	return res;
}

void cm_print(mat_complx *m, char * title)
{
	int i, j;

	switch (m->type) {
		case MAT_DENSE :
			printf("%s : full dense matrix(%i,%i) in basis %d\nreal part:\n",title,m->row,m->col,m->basis);
			for (i=0; i<m->row; i++) {
				printf(" ");
			    for (j=0; j<m->col; j++) {
			    	printf("%9.6g ",m->data[i+j*m->row].re);
			    }
			    printf("\n");
			}
			printf("imag part:\n");
			for (i=0; i<m->row; i++) {
			    printf(" ");
			    for (j=0; j<m->col; j++) {
			    	printf("%9.6g ",m->data[i+j*m->row].im);
			    }
			    printf("\n");
			}
			break;
		case MAT_DENSE_DIAG :
			printf("%s : diagonal dense matrix(%i,%i) in basis %d\nreal part:\n ",title,m->row,m->col,m->basis);
			for (i=0; i<m->row; i++) {
				printf("%9.6g ",m->data[i].re);
			}
			printf("\nimag part:\n ");
			for (i=0; i<m->row; i++) {
				printf("%9.6g ",m->data[i].im);
			}
			printf("\n");
			break;
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG : {
			int nn, *ic = m->icol, len = m->irow[m->row] - 1;
			complx *z = m->data;
			printf("%s : %s matrix(%i,%i) in basis %d, %d elements\n",title,matrix_type(m->type),m->row,m->col,m->basis,len);
			for (i=0; i<m->row; i++) {
				nn = m->irow[i+1] - m->irow[i];
				for (j=0; j<nn; j++) {
					printf(" [%d, %d] = (%9.6g, %9.6g)\n",i+1,*ic,z->re, z->im);
					ic++;
					z++;
				}
			}
			break; }
	   default :
		   fprintf(stderr,"Error: cm_print - unknown type '%d'\n",m->type);
		   exit(1);
	}
}

mat_complx * cm_adjoint(mat_complx *m)
{
	mat_complx *res;

	switch (m->type) {
		case MAT_DENSE_DIAG :
			assert(m->row == m->col);
			res = cm_dup(m);
			cblas_dscal(m->row,-1.0,(double*)(res->data)+1,2);
			break;
		case MAT_SPARSE_DIAG : {
			int len;
			assert(m->row == m->col);
			res = cm_dup(m);
			len = res->irow[res->row] - 1;
			cblas_dscal(len,-1.0,(double*)(res->data)+1,2);
			break; }
		case MAT_DENSE : {
			int i;
			res = complx_matrix(m->col, m->row, MAT_DENSE,0,m->basis);
			for (i=0; i<m->row; i++) {
				cblas_zcopy(m->col, &(m->data[i]),m->row,&res->data[i*res->row],1);
				cblas_dscal(res->row,-1.0,(double*)(&(res->data[i*res->row]))+1,2);
			}
			break; }
		case MAT_SPARSE : {
			int len = m->irow[m->row] - 1;
			MKL_INT *dumrow, *dumcol, *ic = m->icol;
			dumrow = (MKL_INT*)malloc((m->col+1+len)*sizeof(MKL_INT));
			dumcol = dumrow+m->col+1;
			int i, j, nn=0;
			complx *z;
			res = complx_matrix(m->col, m->row, MAT_SPARSE, len, m->basis);
			memset(dumrow,0,(m->col+1)*sizeof(MKL_INT));
			for (i=0; i<m->row; i++) {
				int n = m->irow[i+1] - m->irow[i];
				for (j=0; j<n; j++) {
					dumrow[*ic]++;
					dumcol[nn] = i+1;
					ic++;
					nn++;
				}
			}
			res->irow[0] = dumrow[0] = 1;
			for (i=1; i<=res->row; i++) {
				res->irow[i] = res->irow[i-1] + dumrow[i];
				//DEBUGPRINT("   cm_adjoint res->irow[%d]=%d\n",i,res->irow[i]);
				dumrow[i] = res->irow[i];
			}
			z = m->data;
			for (i=0; i<len; i++) {
				j = dumrow[m->icol[i]-1]-1;
				res->data[j].re = z->re;
				res->data[j].im = -z->im;
				res->icol[j] = dumcol[i];
				z++;
				dumrow[m->icol[i]-1]++;
				//DEBUGPRINT("   cm_adjoint res->icol[%d] = %d\n",j,res->icol[j]);
			}
			free(dumrow);
			break; }
	   default :
		   fprintf(stderr,"Error: cm_adjoint - unknown type '%d'\n",m->type);
		   exit(1);
	}
	return res;
}

void cm_adjointi(mat_complx *m)
{
	switch (m->type) {
		case MAT_DENSE_DIAG :
			assert(m->row == m->col);
			cblas_dscal(m->row,-1.0,(double*)(m->data)+1,2);
			break;
		case MAT_SPARSE_DIAG :
			assert(m->row == m->col);
			cblas_dscal(m->irow[m->row]-1,-1.0,(double*)(m->data)+1,2);
			break;
		case MAT_DENSE :
		case MAT_SPARSE : {
			mat_complx *dum;
			dum = cm_adjoint(m);
			cm_swap_innards_and_destroy(m,dum);
			break; }
	   default :
		   fprintf(stderr,"Error: cm_adjointi - unknown type '%d'\n",m->type);
		   exit(1);
	}
}

mat_complx * cm_transpose(mat_complx *m)
{
	mat_complx *res;

	switch (m->type) {
		case MAT_DENSE_DIAG :
			res = cm_dup(m);
			break;
		case MAT_SPARSE_DIAG : {
			res = cm_dup(m);
			break; }
		case MAT_DENSE : {
			int i;
			res = complx_matrix(m->col, m->row, MAT_DENSE,0,m->basis);
			for (i=0; i<m->row; i++) {
				cblas_zcopy(m->col, &(m->data[i]),m->row,&res->data[i*res->row],1);
			}
			break; }
		case MAT_SPARSE : {
			int len = m->irow[m->row] - 1;
			MKL_INT *dumrow, *dumcol, *ic = m->icol;
			dumrow = (MKL_INT*)malloc((m->col+1+len)*sizeof(MKL_INT));
			dumcol = dumrow + m->col+1;
			int i, j, nn=0;
			complx *z;
			res = complx_matrix(m->col, m->row, MAT_SPARSE, len,m->basis);
			memset(dumrow,0,(m->col+1)*sizeof(MKL_INT));
			for (i=0; i<m->row; i++) {
				int n = m->irow[i+1] - m->irow[i];
				for (j=0; j<n; j++) {
					dumrow[*ic]++;
					dumcol[nn] = i+1;
					ic++;
					nn++;
				}
			}
			res->irow[0] = dumrow[0] = 1;
			for (i=1; i<=res->row; i++) {
				res->irow[i] = res->irow[i-1] + dumrow[i];
				//DEBUGPRINT("   cm_transpose res->irow[%d]=%d\n",i,res->irow[i]);
				dumrow[i] = res->irow[i];
			}
			z = m->data;
			for (i=0; i<len; i++) {
				j = dumrow[m->icol[i]-1]-1;
				res->data[j].re = z->re;
				res->data[j].im = z->im;
				res->icol[j] = dumcol[i];
				z++;
				dumrow[m->icol[i]-1]++;
				//DEBUGPRINT("   cm_transpose res->icol[%d] = %d\n",j,res->icol[j]);
			}
			free(dumrow);
			break; }
	   default :
		   fprintf(stderr,"Error: cm_transpose - unknown type '%d'\n",m->type);
		   exit(1);
	}
	return res;
}



void cm_filter(mat_complx *m, mat_complx *mask)
{
	int i, j, n;
	complx zz, *z = m->data;

	assert(m->basis == mask->basis);
	//cm_print(m,"cm_filter input matrix");
	//cm_print(mask,"cm_filter mask matrix");
	switch (m->type) {
		case MAT_DENSE :
			for (i=1; i<=m->row; i++) {
				for (j=1; j<=m->col; j++) {
					zz = cm_getelem(mask,j,i);
					if (zz.re*zz.re+zz.im*zz.im < TINY) *z = Cnull;
					z++;
				}
			}
			break;
		case MAT_DENSE_DIAG :
			for (i=1; i<=m->row; i++) {
				zz = cm_getelem(mask,i,i);
				if (zz.re*zz.re+zz.im*zz.im < TINY) *z = Cnull;
				z++;
			}
			break;
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG : {
			int *ic = m->icol;
			for (i=0; i<m->row; i++) {
				n = m->irow[i+1] - m->irow[i];
				for (j=0; j<n; j++) {
					zz = cm_getelem(mask,i+1,*ic);
					if (zz.re*zz.re+zz.im*zz.im < TINY) *z = Cnull;
					z++;
					ic++;
				}
			}
			cm_sparse(m,SPARSE_TOL);
			break; }
	   default :
		   fprintf(stderr,"Error: cm_filter - unknown type '%d'\n",m->type);
		   exit(1);
	}
	//cm_print(m,"cm_filter result");
}

/* this replaces original use of cm_or. That only converted nonzero real values to ones */
void cm_nnz2one(mat_complx *m)
{
	int i;
	complx *z = m->data;

	switch (m->type) {
		case MAT_DENSE :
			for (i=0; i<m->row*m->col; i++) {
				if ( fabs(z->re) > TINY )
					*z = Complx(1.0,0.0);
				else
					*z = Cnull;
				z++;
			}
			break;
		case MAT_DENSE_DIAG :
			for (i=0; i<m->row; i++) {
				if ( fabs(z->re) > TINY )
					*z = Complx(1.0,0.0);
				else
					*z = Cnull;
				z++;
			}
			break;
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG : {
			int dim = m->irow[m->row] - 1;
			for (i=0; i<dim; i++) {
				if ( fabs(z->re) > TINY )
					*z = Complx(1.0,0.0);
				else
					*z = Cnull;
				z++;
			}
			cm_sparse(m,SPARSE_TOL);
			break;}
	   default :
		   fprintf(stderr,"Error: cm_nnz2one - unknown type '%d'\n",m->type);
		   exit(1);
	}
}

void cm_muld(mat_complx *m, double d)
{
	switch (m->type) {
		case MAT_DENSE :
			cblas_zdscal(m->row*m->col,d,m->data,1);
			break;
		case MAT_DENSE_DIAG :
			cblas_zdscal(m->row,d,m->data,1);
			break;
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG :
			cblas_zdscal(m->irow[m->row]-1,d,m->data,1);
			break;
	   default :
		   fprintf(stderr,"Error: cm_muld - unknown type '%d'\n",m->type);
		   exit(1);
	}
}

void cm_mulc(mat_complx *m, complx z)
{
	switch (m->type) {
		case MAT_DENSE :
			cblas_zscal(m->row*m->col,&z,m->data,1);
			break;
		case MAT_DENSE_DIAG :
			cblas_zscal(m->row,&z,m->data,1);
			break;
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG :
			cblas_zscal(m->irow[m->row]-1,&z,m->data,1);
			break;
	   default :
		   fprintf(stderr,"Error: cm_mulc - unknown type '%d'\n",m->type);
		   exit(1);
	}
}

/* sparse complex res = a + z*b */
mat_complx * simpson_zcsradd(mat_complx *a, complx z, mat_complx *b)
{
	assert(a->basis == b->basis);
#ifdef INTEL_MKL
		char trans='N';
		int sort=0;
		int request=1, info=0, nzmax=a->row*a->col;
		int type;
		mat_complx *res;

		if ( (a->type == MAT_SPARSE) || (b->type == MAT_SPARSE) ) {
			type = MAT_SPARSE;
		} else {
			type = MAT_SPARSE_DIAG;
		}
		res = complx_matrix(a->row,a->col,type,0,a->basis);
		mkl_zcsradd(&trans, &request, &sort, &(a->row), &(a->col), a->data, a->icol, a->irow, &z, b->data, b->icol, b->irow, res->data, res->icol, res->irow, &nzmax, &info);
		cm_change_nnz(res,res->irow[res->row]-1);
		request = 2;
		mkl_zcsradd(&trans, &request, &sort, &(a->row), &(a->col), a->data, a->icol, a->irow, &z, b->data, b->icol, b->irow, res->data, res->icol, res->irow, &nzmax, &info);
		if (info != 0) {
			fprintf(stderr,"Error: zcsradd failed with the code '%d'\n",info);
			exit(1);
		}
		return res;
#else
		fprintf(stderr,"Error: sparse matrix routine dcsradd not compiled\n");
		exit(1);
		return NULL;
#endif
}

void cm_multod(mat_complx *m1, mat_complx *m2, double d)
{
	complx z;

	if (m1->row != m2->row || m1->col != m2->col) {
		fprintf(stderr,"Error: cm_multod - dimension mismatch (%d,%d) + (%d,%d)\n",m1->row,m1->col, m2->row, m2->col);
		exit(1);
	}
	assert(m1->basis == m2->basis);

	switch (m1->type + m2->type) {
	case MAT_DENSE+MAT_DENSE :
		z = Complx(d,0.0);
		cblas_zaxpy(m1->row*m1->col,&z,m2->data,1,m1->data,1);
		break;
	case MAT_DENSE_DIAG+MAT_DENSE_DIAG :
		z = Complx(d,0.0);
		cblas_zaxpy(m1->row,&z,m2->data,1,m1->data,1);
		break;
	case MAT_DENSE+MAT_DENSE_DIAG : { /* dense, diagonal + full */
		mat_complx *res, *dadd;
		complx *ds1, *ds2;
		double dd;
		int i;
		if (m1->type == MAT_DENSE) {
			res = m1;
			dadd = m2;
			dd = d;
		} else {
			res = cm_dup(m2);
			cm_muld(res,d);
			dadd = m1;
			dd = 1.0;
		}
		ds1 = res->data;
		ds2 = dadd->data;
		for (i=0; i<res->row; i++) {
		      ds1->re += (ds2->re)*dd;
		      ds1->im += (ds2->im)*dd;
		      ds1 += (res->row+1);
		      ds2++;
		}
		if (res != m1) {
			cm_swap_innards_and_destroy(m1,res);
		}
		break; }
	case MAT_DENSE+MAT_SPARSE :
	case MAT_DENSE+MAT_SPARSE_DIAG : { /* sparse plus dense full */
		mat_complx *res, *dadd;
		complx *ds, *zs;
		double dd;
		int i, r, c, n;
		MKL_INT *ic;
		if (m1->type == MAT_DENSE) {
			res = m1;
			dadd = m2;
			dd = d;
		} else {
			res = cm_dup(m2);
			cm_muld(res,d);
			dadd = m1;
			dd = 1.0;
		}
		ic = dadd->icol;
		ds = dadd->data;
		for (r=0; r<dadd->row; r++) {
			n = dadd->irow[r+1] - dadd->irow[r];
			for (i=0; i<n; i++) {
				c = (*ic) - 1;
				ic++;
				zs = res->data + (r+c*res->row);
				zs->re += dd*(ds->re);
				zs->im += dd*(ds->im);
				ds++;
			}
		}
		if (res != m1) {
			cm_swap_innards_and_destroy(m1,res);
		}
		break; }
	case MAT_SPARSE+MAT_SPARSE:
	case MAT_SPARSE+MAT_SPARSE_DIAG:
	case MAT_SPARSE_DIAG+MAT_SPARSE_DIAG: { /* both sparse */
		mat_complx *res;
		res = simpson_zcsradd(m1,Complx(d,0.0),m2);
		cm_swap_innards_and_destroy(m1,res);
		break; }
	case MAT_DENSE_DIAG+MAT_SPARSE :
	case MAT_DENSE_DIAG+MAT_SPARSE_DIAG : {
		mat_complx *res, *dum;
		if (m1->type == MAT_DENSE_DIAG) {
			dum = m2;
			cm_sparse(m1,SPARSE_TOL);
		} else {
			dum = cm_dup(m2);
			cm_sparse(dum,SPARSE_TOL);
		}
		res = simpson_zcsradd(m1,Complx(d,0.0),dum);
		if (dum != m2) free_complx_matrix(dum);
		cm_swap_innards_and_destroy(m1,res);
		break; }
	default :
		fprintf(stderr,"Error: cm_multod - unknown matrix types '%d', '%d'\n",m1->type,m2->type);
		exit(1);
	}
	if (m1->type == MAT_SPARSE || m1->type == MAT_SPARSE_DIAG) cm_sparse(m1,SPARSE_TOL);
}

void cm_multoc(mat_complx *m1, mat_complx *m2, complx z)
{
	if (m1->row != m2->row || m1->col != m2->col) {
		fprintf(stderr,"Error: cm_multoc - dimension mismatch (%d,%d) + (%d,%d)\n",m1->row,m1->col, m2->row, m2->col);
		exit(1);
	}
	assert(m1->basis == m2->basis);

	switch (100*m1->type + m2->type) {
	case 100*MAT_DENSE+MAT_DENSE :
		cblas_zaxpy(m1->row*m1->col,&z,m2->data,1,m1->data,1);
		break;
	case 100*MAT_DENSE_DIAG+MAT_DENSE_DIAG :
		cblas_zaxpy(m1->row,&z,m2->data,1,m1->data,1);
		break;
	case 100*MAT_DENSE+MAT_DENSE_DIAG : { /* dense, diagonal + full */
		complx *ds1 = m1->data;
		complx *ds2 = m2->data;
		int i;
		for (i=0; i<m1->row; i++) {
			double rre = ds2->re;
			double iim = ds2->im;
			ds1->re += rre*z.re - iim*z.im;
			ds1->im += rre*z.im + iim*z.re;
			ds1 += (m1->row+1);
			ds2++;
		}
		break;}
	case 100*MAT_DENSE_DIAG+MAT_DENSE : {
		mat_complx *res;
		complx *ds1, *ds2;
		int i;
		res = cm_dup(m2);
		cm_mulc(res,z);
		ds1 = res->data;
		ds2 = m1->data;
		for (i=0; i<res->row; i++) {
			ds1->re += ds2->re;
			ds1->im += ds2->im;
			ds1 += (res->row+1);
			ds2++;
		}
		cm_swap_innards_and_destroy(m1,res);
		break; }
	case 100*MAT_DENSE+MAT_SPARSE :
	case 100*MAT_DENSE+MAT_SPARSE_DIAG : { /* sparse plus dense full */
		MKL_INT *ic = m2->icol;
		complx *ds = m2->data;
		int i, r, c, n;
		for (r=0; r<m2->row; r++) {
			n = m2->irow[r+1] - m2->irow[r];
			for (i=0; i<n; i++) {
				c = (*ic) - 1;
				ic++;
				complx *zs = m1->data + (r+c*m1->row);
				zs->re += z.re*ds->re - z.im*ds->im;
				zs->im += z.re*ds->im + z.im*ds->re;
				ds++;
			}
		}
		break;}
	case 100*MAT_SPARSE+MAT_DENSE:
	case 100*MAT_SPARSE_DIAG+MAT_DENSE: {
		mat_complx *res;
		complx *ds, *zs;
		int i, r, c, n;
		MKL_INT *ic;
		res = cm_dup(m2);
		cm_mulc(res,z);
		ic = m1->icol;
		ds = m1->data;
		for (r=0; r<m1->row; r++) {
			n = m1->irow[r+1] - m1->irow[r];
			for (i=0; i<n; i++) {
				c = (*ic) - 1;
				ic++;
				zs = res->data + (r+c*res->row);
				zs->re += ds->re;
				zs->im += ds->im;
				ds++;
			}
		}
		cm_swap_innards_and_destroy(m1,res);
		break; }
	case 100*MAT_SPARSE+MAT_SPARSE:
	case 100*MAT_SPARSE+MAT_SPARSE_DIAG:
	case 100*MAT_SPARSE_DIAG + MAT_SPARSE:
	case 100*MAT_SPARSE_DIAG+MAT_SPARSE_DIAG: { /* both sparse */
		mat_complx *res;
		res = simpson_zcsradd(m1,z,m2);
		cm_swap_innards_and_destroy(m1,res);
		break; }
	case 100*MAT_DENSE_DIAG+MAT_SPARSE :
	case 100*MAT_DENSE_DIAG+MAT_SPARSE_DIAG : {
		mat_complx *res;
		cm_sparse(m1,SPARSE_TOL);
		res = simpson_zcsradd(m1,z,m2);
		cm_swap_innards_and_destroy(m1,res);
		break; }
	case 100*MAT_SPARSE + MAT_DENSE_DIAG:
	case 100*MAT_SPARSE_DIAG + MAT_DENSE_DIAG: {
		mat_complx *dum, *res;
		dum = cm_dup(m2);
		cm_sparse(dum,SPARSE_TOL);
		res = simpson_zcsradd(m1,z,dum);
		free_complx_matrix(dum);
		cm_swap_innards_and_destroy(m1,res);
		break;}
	default :
		fprintf(stderr,"Error: cm_multoc - unknown matrix types '%d', '%d'\n",m1->type,m2->type);
		exit(1);
	}
	if (m1->type == MAT_SPARSE || m1->type == MAT_SPARSE_DIAG) cm_sparse(m1,SPARSE_TOL);
}

/* to a complex matrix, add z*real_matrix */
void cm_multocr(mat_complx *m1, mat_double *m2, complx z)
{
	if (m1->row != m2->row || m1->col != m2->col) {
		fprintf(stderr,"Error: cm_multocr - dimension mismatch (%d,%d) + (%d,%d)\n",m1->row,m1->col, m2->row, m2->col);
		exit(1);
	}
	assert(m1->basis == m2->basis);

	switch (100*m1->type + m2->type) {
	case 100*MAT_DENSE+MAT_DENSE : {
		double *dv = (double*)(m1->data);
		cblas_daxpy(m1->row*m1->col,z.re,m2->data,1,dv,2);
		cblas_daxpy(m1->row*m1->col,z.im,m2->data,1,dv+1,2);
		break; }
	case 100*MAT_DENSE_DIAG+MAT_DENSE_DIAG : {
		double *dv = (double*)(m1->data);
		cblas_daxpy(m1->row,z.re,m2->data,1,dv,2);
		cblas_daxpy(m1->row,z.im,m2->data,1,dv+1,2);
		break; }
	case 100*MAT_DENSE+MAT_DENSE_DIAG : { /* dense, diagonal + full */
		complx *zs = m1->data;
		double *ds = m2->data;
		int i, step = m1->row+1;
		for (i=0; i<m1->row; i++) {
			zs->re += (*ds)*z.re;
			zs->im += (*ds)*z.im;
			zs += step;
			ds++;
		}
		break;}
	case 100*MAT_DENSE_DIAG+MAT_DENSE : {
		mat_complx *res;
		complx *ds1, *ds2;
		int i, step = m2->row+1;
		res = dm_complx(m2);
		cm_mulc(res,z);
		ds1 = res->data;
		ds2 = m1->data;
		for (i=0; i<res->row; i++) {
			ds1->re += ds2->re;
			ds1->im += ds2->im;
			ds1 += step;
			ds2++;
		}
		cm_swap_innards_and_destroy(m1,res);
		break; }
	case 100*MAT_DENSE+MAT_SPARSE :
	case 100*MAT_DENSE+MAT_SPARSE_DIAG : { /* sparse plus dense full */
		MKL_INT *ic = m2->icol;
		double *ds = m2->data;
		int i, r, c, n;
		for (r=0; r<m2->row; r++) {
			n = m2->irow[r+1] - m2->irow[r];
			for (i=0; i<n; i++) {
				c = (*ic) - 1;
				ic++;
				complx *zs = m1->data + (r+c*m1->row);
				zs->re += z.re*(*ds);
				zs->im += z.im*(*ds);
				ds++;
			}
		}
		break;}
	case 100*MAT_SPARSE+MAT_DENSE:
	case 100*MAT_SPARSE_DIAG+MAT_DENSE: {
		mat_complx *res;
		complx *ds, *zs;
		int i, r, c, n;
		MKL_INT *ic;
		res = dm_complx(m2);
		cm_mulc(res,z);
		ic = m1->icol;
		ds = m1->data;
		for (r=0; r<m1->row; r++) {
			n = m1->irow[r+1] - m1->irow[r];
			for (i=0; i<n; i++) {
				c = (*ic) - 1;
				ic++;
				zs = res->data + (r+c*res->row);
				zs->re += ds->re;
				zs->im += ds->im;
				ds++;
			}
		}
		cm_swap_innards_and_destroy(m1,res);
		break; }
	case 100*MAT_SPARSE+MAT_SPARSE:
	case 100*MAT_SPARSE+MAT_SPARSE_DIAG:
	case 100*MAT_SPARSE_DIAG + MAT_SPARSE:
	case 100*MAT_SPARSE_DIAG+MAT_SPARSE_DIAG: { /* both sparse */
		mat_complx *res, *dum;
		dum = dm_complx(m2);
		res = simpson_zcsradd(m1,z,dum);
		free_complx_matrix(dum);
		cm_swap_innards_and_destroy(m1,res);
		break; }
	case 100*MAT_DENSE_DIAG+MAT_SPARSE :
	case 100*MAT_DENSE_DIAG+MAT_SPARSE_DIAG : {
		mat_complx *res, *dum;
		cm_sparse(m1,SPARSE_TOL);
		dum = dm_complx(m2);
		res = simpson_zcsradd(m1,z,dum);
		free_complx_matrix(dum);
		cm_swap_innards_and_destroy(m1,res);
		break; }
	case 100*MAT_SPARSE + MAT_DENSE_DIAG:
	case 100*MAT_SPARSE_DIAG + MAT_DENSE_DIAG: {
		mat_complx *dum, *res;
		dum = dm_complx(m2);
		cm_sparse(dum,SPARSE_TOL);
		res = simpson_zcsradd(m1,z,dum);
		free_complx_matrix(dum);
		cm_swap_innards_and_destroy(m1,res);
		break;}
	default :
		fprintf(stderr,"Error: cm_multocr - unknown matrix types '%d', '%d'\n",m1->type,m2->type);
		exit(1);
	}
	if (m1->type == MAT_SPARSE || m1->type == MAT_SPARSE_DIAG) cm_sparse(m1,SPARSE_TOL);
}


void cm_addto(mat_complx *res, mat_complx *m)
{
	cm_multod(res,m,1.0);
}

void cm_subfrom(mat_complx *res, mat_complx *m)
{
	cm_multod(res,m,-1.0);
}

/* complex sparse res = a * b  , both sparse */
mat_complx * simpson_zcsrmultcsr(mat_complx *a, mat_complx *b)
{
	assert(a->basis == b->basis);

#ifdef INTEL_MKL
		const char trans='N';
		const int sort=0;
		int request=1, info=0, nnzmax=a->row*a->col;
		int type;
		mat_complx *res;

		if ( (a->type == MAT_SPARSE) || (b->type == MAT_SPARSE) ) {
			type = MAT_SPARSE;
		} else {
			type = MAT_SPARSE_DIAG;
		}
		res = complx_matrix(a->row, b->col, type, 0, a->basis);
		mkl_zcsrmultcsr(&trans, &request, &sort, &(a->row), &(a->col), &(b->col), a->data, a->icol, a->irow, b->data, b->icol, b->irow, res->data, res->icol, res->irow, &nnzmax, &info);
		cm_change_nnz(res,res->irow[res->row]-1);
		request = 2;
		mkl_zcsrmultcsr(&trans, &request, &sort, &(a->row), &(a->col), &(b->col), a->data, a->icol, a->irow, b->data, b->icol, b->irow, res->data, res->icol, res->irow, &nnzmax, &info);
		if (info != 0) {
			fprintf(stderr,"Error: zcsrmultcsr failed with the code '%d'\n",info);
			exit(1);
		}
		return res;
#else
		fprintf(stderr,"Error: sparse matrix routine zcsrmultcsr not compiled\n");
		exit(1);
		return NULL;
#endif
}

/* complex dense res = sparse a * dense b  */
mat_complx * simpson_zcsrmm(mat_complx *a, mat_complx *b)
{
	assert(a->basis == b->basis);
#ifdef INTEL_MKL
		char trans='N', matdescra[]="GNNFUU";
		complx zone, zzero;
		mat_complx *res;

		//DEBUGPRINT("Warning!!! Please check sparse * dense multiplication result!\n");
		// checked 8.9.2010 Zd.
		zone = Complx(1.0,0.0); zzero = Cnull;
		res = complx_matrix(a->row,b->col,MAT_DENSE,0,a->basis);
		//dm_zero(res);
		mkl_zcsrmm(&trans, &(a->row), &(res->col), &(a->col), (MKL_Complex16*)&zone, matdescra, (MKL_Complex16*)(a->data), a->icol, a->irow, &(a->irow[1]), (MKL_Complex16*)(b->data), &(b->row), (MKL_Complex16*)&zzero, (MKL_Complex16*)(res->data), &(res->row));
		//cm_print(a,"A");
		//cm_print(b,"B");
		//cm_print(res,"res = A*B");
		//exit(1);
		res->basis = a->basis;
		return res;
#else
		fprintf(stderr,"Error: sparse matrix routine zcsrmmm not compiled\n");
		exit(1);
		return NULL;
#endif
}


mat_complx * cm_mul(mat_complx *m1, mat_complx *m2)
{
	mat_complx *res;
	int i;
	complx zone, zzero;

	if (m1->col != m2->row) {
		fprintf(stderr,"Error: cm_mul - dimension mismatch (%d,%d)x(%d,%d)\n",m1->row,m1->col,m2->row,m2->col);
		exit(1);
	}
	assert(m1->basis == m2->basis);

	zone = Complx(1.0,0.0); zzero = Cnull;

	switch (m1->type*100+m2->type) {
		case 100*MAT_DENSE+MAT_DENSE :
			res = complx_matrix(m1->row,m2->col,MAT_DENSE,0,m1->basis);
			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m1->row,res->col,m1->col,&zone,m1->data,m1->row,m2->data,m2->row,&zzero,res->data,res->row);
			break;
		case 100*MAT_DENSE+MAT_DENSE_DIAG :
			assert(m2->row == m2->col);
			res = cm_dup(m1);
			for (i=0; i<res->col; i++) {
				cblas_zscal(res->row,&(m2->data[i]),&(res->data[i*res->row]),1);
			}
			break;
		case 100*MAT_DENSE+MAT_SPARSE : {
			DEBUGPRINT("warning from cm_mul: converting sparse to dense\n");
			mat_complx *dum;
			dum = cm_dup(m2);
			cm_dense(dum);
			res = complx_matrix(m1->row,dum->col,MAT_DENSE,0,m1->basis);
			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m1->row,res->col,m1->col,&zone,m1->data,m1->row,dum->data,dum->row,&zzero,res->data,res->row);
			free_complx_matrix(dum);
			break; }
		case 100*MAT_DENSE+MAT_SPARSE_DIAG : {
			complx z;
			assert(m2->row == m2->col);
			res = cm_dup(m1);
			for (i=0; i<res->col; i++) {
				z = cm_getelem(m2,i+1,i+1);
				cblas_zscal(res->row,&z,&(res->data[i*res->row]),1);
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_DENSE :
			assert(m1->row == m1->col);
			res = cm_dup(m2);
			for (i=0; i<res->row; i++) {
				cblas_zscal(res->col,&(m1->data[i]),&(res->data[i]),res->row);
			}
			break;
		case 100*MAT_DENSE_DIAG+MAT_DENSE_DIAG : {
			assert( m1->row == m1->col && m2->row == m2->col);
			res = complx_matrix(m1->row,m2->col,MAT_DENSE_DIAG,0,m1->basis);
			complx *id1 = m1->data, *id2 = m2->data, *id3 = res->data;
			for (i=0; i<m1->row; i++) {
				//*id3 = Cmul(*id1,*id2);
				id3->re = id1->re*id2->re - id1->im*id2->im;
				id3->im = id1->re*id2->im + id1->im*id2->re;
				id1++;
				id2++;
				id3++;
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_SPARSE : {
			mat_complx *dum;
			assert(m1->row == m1->col);
			dum = cm_dup(m1);
			cm_sparse(dum,SPARSE_TOL);
			res = simpson_zcsrmultcsr(dum,m2);
			free_complx_matrix(dum);
			break; }
		case 100*MAT_DENSE_DIAG+MAT_SPARSE_DIAG :
			assert(m1->row == m1->col && m2->row == m2->col);
			res = complx_matrix(m1->row, m2->col,MAT_DENSE_DIAG,0,m1->basis);
			for (i=0; i<m1->row; i++) {
				res->data[i] = Cmul(m1->data[i],cm_getelem(m2,i+1,i+1));
			}
			break;
		case 100*MAT_SPARSE+MAT_DENSE :
			res = simpson_zcsrmm(m1,m2);
			break;
		case 100*MAT_SPARSE+MAT_DENSE_DIAG : {
			assert(m2->row == m2->col);
			mat_complx *dum;
			dum = cm_dup(m2);
			cm_sparse(dum, SPARSE_TOL);
			res = simpson_zcsrmultcsr(m1,dum);
			free_complx_matrix(dum);
			break;}
		case 100*MAT_SPARSE+MAT_SPARSE :
		case 100*MAT_SPARSE+MAT_SPARSE_DIAG :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE_DIAG :
			res = simpson_zcsrmultcsr(m1,m2);
			break;
		case 100*MAT_SPARSE_DIAG+MAT_DENSE : {
			complx z;
			assert(m1->row == m1->col);
			res = cm_dup(m2);
			for (i=0; i<res->row; i++) {
				z = cm_getelem(m1,i+1,i+1);
				cblas_zscal(res->col,&z,&(res->data[i]),res->row);
			}
			break; }
		case 100*MAT_SPARSE_DIAG+MAT_DENSE_DIAG : {
			assert(m1->row == m1->col && m2->row == m2->col);
			int len = m1->irow[m1->row]-1;
			res = complx_matrix(m1->row,m2->col,MAT_SPARSE_DIAG,len,m1->basis);
			memcpy(res->icol,m1->icol,len*sizeof(MKL_INT));
			memcpy(res->irow,m1->irow,(m1->row+1)*sizeof(MKL_INT));
			for (i=0; i<len; i++) {
				res->data[i] = Cmul(m1->data[i],m2->data[m1->icol[i]-1]);
			}
			break; }
	   default :
		   fprintf(stderr,"Error: cm_mul - unknown types '%d' and '%d'\n",m1->type, m2->type);
		   exit(1);
	}
	if (res->type == MAT_SPARSE || res->type == MAT_SPARSE_DIAG) cm_sparse(res,SPARSE_TOL);
	return res;
}

void cm_multo(mat_complx *m1, mat_complx *m2)
{
	mat_complx *res;
	int i;
	complx zone, zzero;

	if (m1->col != m2->row) {
		fprintf(stderr,"Error: cm_multo - dimension mismatch (%d,%d)x(%d,%d)\n",m1->row,m1->col,m2->row, m2->col);
		exit(1);
	}
	assert(m1->basis == m2->basis);

	zone = Complx(1.0,0.0); zzero = Cnull;
	switch (m1->type*100+m2->type) {
		case 100*MAT_DENSE+MAT_DENSE :
			res = complx_matrix(m1->row,m2->col,MAT_DENSE,0,m1->basis);
			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m1->row,res->col,m1->col,&zone,m1->data,m1->row,m2->data,m2->row,&zzero,res->data,res->row);
			cm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_DENSE+MAT_DENSE_DIAG :
			assert(m2->row == m2->col);
			for (i=0; i<m1->col; i++) {
				cblas_zscal(m1->row,&(m2->data[i]),&(m1->data[i*m1->row]),1);
			}
			break;
		case 100*MAT_DENSE+MAT_SPARSE : {
			DEBUGPRINT("warning from cm_multo: converting sparse to dense\n");
			mat_complx *dum;
			dum = cm_dup(m2);
			cm_dense(dum);
			res = complx_matrix(m1->row,dum->col,MAT_DENSE,0,m1->basis);
			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m1->row,res->row,m1->col,&zone,m1->data,m1->row,dum->data,dum->row,&zzero,res->data,res->row);
			free_complx_matrix(dum);
			cm_swap_innards_and_destroy(m1,res);
			break; }
		case 100*MAT_DENSE+MAT_SPARSE_DIAG : {
			complx z;
			assert(m2->row == m2->col);
			for (i=0; i<m1->col; i++) {
				z = cm_getelem(m2,i+1,i+1);
				cblas_zscal(m1->row,&z,&(m1->data[i*m1->row]),1);
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_DENSE :
			assert(m1->row == m1->col);
			res = cm_dup(m2);
			for (i=0; i<res->row; i++) {
				cblas_zscal(res->col,&(m1->data[i]),&(res->data[i]),res->row);
			}
			cm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_DENSE_DIAG+MAT_DENSE_DIAG : {
			assert(m1->row == m1->col && m2->row == m2->col);
			complx *id1 = m1->data, *id2 = m2->data;
			for (i=0; i<m1->row; i++) {
				double rre = id1->re;
				//*id1 = Cmul(*id1,*id2);
				id1->re = rre*id2->re - id1->im*id2->im;
				id1->im = rre*id2->im + id1->im*id2->re;
				id1++;
				id2++;
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_SPARSE :
			cm_sparse(m1,SPARSE_TOL);
			res = simpson_zcsrmultcsr(m1,m2);
			cm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_DENSE_DIAG+MAT_SPARSE_DIAG :
			assert(m1->row == m1->col && m2->row == m2->col);
			for (i=0; i<m1->row; i++) {
				m1->data[i] = Cmul(m1->data[i],cm_getelem(m2,i+1,i+1));
			}
			break;
		case 100*MAT_SPARSE+MAT_DENSE :
			res = simpson_zcsrmm(m1,m2);
			cm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_SPARSE+MAT_DENSE_DIAG : {
			mat_complx *dum;
			dum = cm_dup(m2);
			cm_sparse(dum, SPARSE_TOL);
			res = simpson_zcsrmultcsr(m1,dum);
			free_complx_matrix(dum);
			cm_swap_innards_and_destroy(m1,res);
			break; }
		case 100*MAT_SPARSE+MAT_SPARSE :
		case 100*MAT_SPARSE+MAT_SPARSE_DIAG :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE_DIAG :
			res = simpson_zcsrmultcsr(m1,m2);
			cm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_SPARSE_DIAG+MAT_DENSE : {
			complx z;
			assert(m1->row == m1->col);
			res = cm_dup(m2);
			for (i=0; i<res->row; i++) {
				z = cm_getelem(m1,i+1,i+1);
				cblas_zscal(res->col,&z,&(res->data[i]),res->row);
			}
			cm_swap_innards_and_destroy(m1,res);
			break; }
		case 100*MAT_SPARSE_DIAG+MAT_DENSE_DIAG :
			assert(m1->row == m1->col && m2->row == m2->col);
			for (i=0; i<m1->irow[m1->row]-1; i++) {
				m1->data[i] = Cmul(m1->data[i],m2->data[m1->icol[i]-1]);
			}
			break;
	   default :
		   fprintf(stderr,"Error: cm_multo - unknown types '%d' and '%d'\n",m1->type, m2->type);
		   exit(1);
	}
	if (m1->type == MAT_SPARSE || m1->type == MAT_SPARSE_DIAG) cm_sparse(m1,SPARSE_TOL);
}

void cm_multo_rev(mat_complx *m1, mat_complx *m2)
{
	mat_complx *res;
	int i;
	complx zone, zzero;

	//DEBUGPRINT("cm_multo_rev: m1 = %s, m2 = %s\n",matrix_type(m1->type),matrix_type(m2->type));
	if (m2->col != m1->row) {
		fprintf(stderr,"Error: cm_multo_rev - dimension mismatch m1=m2xm1, i.e. (%d,%d)x(%d,%d)\n",m2->row,m2->col, m1->row, m1->col);
		exit(1);
	}
	assert(m1->basis == m2->basis);

	zone = Complx(1.0,0.0); zzero = Cnull;
	switch (m1->type*100+m2->type) {
		case 100*MAT_DENSE+MAT_DENSE :
			res = complx_matrix(m2->row,m1->col,MAT_DENSE,0,m1->basis);
			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m2->row,res->col,m2->col,&zone,m2->data,m2->row,m1->data,m1->row,&zzero,res->data,res->row);
			cm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_DENSE+MAT_DENSE_DIAG :
			assert(m2->row == m2->col);
			for (i=0; i<m1->row; i++) {
				cblas_zscal(m1->col,&(m2->data[i]),&(m1->data[i]),m1->row);
			}
			break;
		case 100*MAT_DENSE+MAT_SPARSE :
			res = simpson_zcsrmm(m2,m1);
			cm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_DENSE+MAT_SPARSE_DIAG : {
			complx z;
			assert(m2->row == m2->col);
			for (i=0; i<m1->row; i++) {
				z = cm_getelem(m2,i+1,i+1);
				cblas_zscal(m1->col,&z,&(m1->data[i]),m1->row);
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_DENSE :
			assert(m1->row == m1->col);
			res = cm_dup(m2);
			for (i=0; i<res->col; i++) {
				cblas_zscal(res->row,&(m1->data[i]),&(res->data[i*res->row]),1);
			}
			cm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_DENSE_DIAG+MAT_DENSE_DIAG : {
			assert(m1->row == m1->col && m2->row == m2->col);
			complx *id1 = m1->data, *id2 = m2->data;
			for (i=0; i<m1->row; i++) {
				//*id1 = Cmul(*id1,*id2);
				double rre = id1->re;
				id1->re = rre*id2->re - id1->im*id2->im;
				id1->im = rre*id2->im + id1->im*id2->re;
				id1++;
				id2++;
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_SPARSE :
			cm_sparse(m1,SPARSE_TOL);
			res = simpson_zcsrmultcsr(m2,m1);
			cm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_DENSE_DIAG+MAT_SPARSE_DIAG :
			assert(m1->row == m1->col && m2->row == m2->col);
			for (i=0; i<m1->row; i++) {
				m1->data[i] = Cmul(m1->data[i],cm_getelem(m2,i+1,i+1));
			}
			break;
		case 100*MAT_SPARSE+MAT_DENSE : {
			DEBUGPRINT("warning from cm_multo_rev: converting sparse to dense\n");
			mat_complx *dum;
			dum = cm_dup(m1);
			cm_dense(dum);
			res = complx_matrix(m2->row,m1->col,MAT_DENSE,0,m1->basis);
			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m2->row,res->col,m2->col,&zone,m2->data,m2->row,dum->data,dum->row,&zzero,res->data,res->row);
			free_complx_matrix(dum);
			cm_swap_innards_and_destroy(m1,res);
			break; }
		case 100*MAT_SPARSE+MAT_DENSE_DIAG : {
			mat_complx *dum;
			dum = cm_dup(m2);
			cm_sparse(dum, SPARSE_TOL);
			res = simpson_zcsrmultcsr(dum,m1);
			free_complx_matrix(dum);
			cm_swap_innards_and_destroy(m1,res);
			break; }
		case 100*MAT_SPARSE+MAT_SPARSE :
		case 100*MAT_SPARSE+MAT_SPARSE_DIAG :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE_DIAG :
			res = simpson_zcsrmultcsr(m2,m1);
			cm_swap_innards_and_destroy(m1,res);
			break;
		case 100*MAT_SPARSE_DIAG+MAT_DENSE : {
			complx z;
			assert(m1->row == m1->col);
			res = cm_dup(m2);
			for (i=0; i<res->col; i++) {
				z = cm_getelem(m1,i+1,i+1);
				cblas_zscal(res->row,&z,&(res->data[i*res->row]),1);
			}
			cm_swap_innards_and_destroy(m1,res);
			break; }
		case 100*MAT_SPARSE_DIAG+MAT_DENSE_DIAG :
			assert(m1->row == m1->col && m2->row == m2->col);
			for (i=0; i<m1->irow[m1->row]-1; i++) {
				m1->data[i] = Cmul(m1->data[i],m2->data[m1->icol[i]-1]);
			}
			break;
	   default :
		   fprintf(stderr,"Error: cm_multo_rev - unknown types '%d' and '%d'\n",m1->type, m2->type);
		   exit(1);
	}
	if (m1->type == MAT_SPARSE || m1->type == MAT_SPARSE_DIAG) cm_sparse(m1,SPARSE_TOL);
}

mat_complx * cm_commutator(mat_complx *a, mat_complx *b)
{
	/* this is not the best code for dense matrices... */
	mat_complx *dum, *res;

	assert(a->basis == b->basis);

	dum = cm_mul(b,a);
	res = cm_mul(a,b);
	cm_multod(res,dum,-1.0);
	free_complx_matrix(dum);
	if (res->type == MAT_SPARSE || res->type == MAT_SPARSE_DIAG) cm_sparse(res, SPARSE_TOL);
	return res;
}

complx cm_trace(mat_complx *a, mat_complx *b)
{
	int i, j, dim = a->row;
	complx res;

	if ( (a->row != a->col) || (b->row != b->col) || (a->row != b->col) ) {
		fprintf(stderr,"Error: cm_trace - non-square matrix or dimension mismatch Tr{(%d,%d)x(%d,%d)}\n",a->row,a->col,b->row,b->col);
		exit(1);
	}

	//DEBUGPRINT("cm_trace: a is %s in basis %d, b is %s in basis %d\n",matrix_type(a->type),a->basis,matrix_type(b->type),b->basis);
	assert(a->basis == b->basis);

	switch (a->type*100+b->type) {
		case 100*MAT_DENSE+MAT_DENSE : {
			complx z;
			res = Cnull;
			for (i=0; i<dim; i++) {
				cblas_zdotu_sub(dim, &(a->data[i*dim]),1,&(b->data[i]),dim,&z);
				//res = Cadd(res, z);
				res.re += z.re;
				res.im += z.im;
			}
			break; }
		case 100*MAT_DENSE+MAT_DENSE_DIAG : {
		    complx *va = a->data, *vb = b->data;
			res = Cnull;
			for (i=0; i<dim; i++) {
				res.re += va->re*vb->re-va->im*vb->im;
				res.im += va->im*vb->re+va->re*vb->im;
				vb++;
				va += (dim+1);
			}
			break; }
		case 100*MAT_DENSE+MAT_SPARSE :
		case 100*MAT_DENSE+MAT_SPARSE_DIAG : {
			int n, *ic = b->icol;
			complx *va, *vb = b->data;
			res = Cnull;
			for (i=0; i<dim; i++) {
				n = b->irow[i+1] - b->irow[i];
				for (j=0; j<n; j++) {
					va = a->data + (*ic-1) + i*dim;
					res.re += va->re*vb->re-va->im*vb->im;
					res.im += va->im*vb->re+va->re*vb->im;
					ic++;
					vb++;
				}
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_DENSE : {
			complx *va = a->data, *vb = b->data;
			res = Cnull;
			for (i=0; i<dim; i++) {
				res.re += va->re*vb->re-va->im*vb->im;
				res.im += va->im*vb->re+va->re*vb->im;
				va++;
				vb += (dim+1);
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_DENSE_DIAG :
			cblas_zdotu_sub(dim, a->data,1,b->data,1,&res);
			break;
		case 100*MAT_DENSE_DIAG+MAT_SPARSE : {
			complx *va = a->data, vb;
			res = Cnull;
			for (i=1; i<=dim; i++) {
				vb = cm_getelem(b,i,i);
				res.re += va->re*vb.re-va->im*vb.im;
				res.im += va->im*vb.re+va->re*vb.im;
				va++;
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_SPARSE_DIAG : {
			complx *vb = b->data, *va;
			int *ic = b->icol;
			res = Cnull;
			j = b->irow[dim]-1;
			for (i=0; i<j; i++) {
				va = a->data + *ic - 1;
				res.re += va->re*vb->re-va->im*vb->im;
				res.im += va->im*vb->re+va->re*vb->im;
				vb++;
				ic++;
			}
			break; }
		case 100*MAT_SPARSE+MAT_DENSE : {
			int n, *ic = a->icol;
			complx *vb, *va = a->data;
			res = Cnull;
			for (i=0; i<dim; i++) {
				n = a->irow[i+1] - a->irow[i];
				for (j=0; j<n; j++) {
					vb = b->data + (*ic-1) + i*dim;
					res.re += va->re*vb->re-va->im*vb->im;
					res.im += va->im*vb->re+va->re*vb->im;
					ic++;
					va++;
				}
			}
			break; }
		case 100*MAT_SPARSE_DIAG+MAT_DENSE : {
			complx *va = a->data, *vb;
			int *ic = a->icol;
			res = Cnull;
			j = a->irow[dim]-1;
			for (i=0; i<j; i++) {
				vb = b->data + (*ic-1) + (*ic-1)*dim;
				res.re += va->re*vb->re-va->im*vb->im;
				res.im += va->im*vb->re+va->re*vb->im;
				va++;
				ic++;
			}
			break; }
		case 100*MAT_SPARSE+MAT_DENSE_DIAG : {
			complx *vb = b->data, va;
			res = Cnull;
			for (i=1; i<=dim; i++) {
				va = cm_getelem(a,i,i);
				res.re += va.re*vb->re-va.im*vb->im;
				res.im += va.im*vb->re+va.re*vb->im;
				vb++;
			}
			break; }
		case 100*MAT_SPARSE+MAT_SPARSE :
		case 100*MAT_SPARSE+MAT_SPARSE_DIAG :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE_DIAG : {
			int la, lb, len, *ic, *ir;
			complx *va, vb;
			mat_complx *dum;
			la = a->irow[dim]-1;
			lb = b->irow[dim]-1;
			if (la < lb) {
				va = a->data;
				ic = a->icol;
				ir = a->irow;
				len = la;
				dum = b;
			} else {
				va = b->data;
				ic = b->icol;
				ir = b->irow;
				len = lb;
				dum = a;
			}
			res = Cnull;
			for (i=0; i<dim; i++) {
				la = ir[i+1] - ir[i];
				for (j=0; j<la; j++) {
					vb = cm_getelem(dum,*ic,i+1);
					res.re += va->re*vb.re-va->im*vb.im;
					res.im += va->im*vb.re+va->re*vb.im;
					va++;
					ic++;
				}
			}
			break; }
		case 100*MAT_SPARSE_DIAG+MAT_DENSE_DIAG : {
			complx *va = a->data, *vb;
			int *ic = a->icol;
			res = Cnull;
			j = a->irow[dim]-1;
			for (i=0; i<j; i++) {
				vb = b->data + *ic - 1;
				res.re += va->re*vb->re-va->im*vb->im;
				res.im += va->im*vb->re+va->re*vb->im;
				va++;
				ic++;
			}
			break; }
	   default :
		   fprintf(stderr,"Error: cm_trace - unknown types '%d' and '%d'\n",a->type, b->type);
		   exit(1);
	}
	return res;

}

/* Tr { A+  B  }   */
complx cm_trace_adjoint(mat_complx *a, mat_complx *b)
{
	int i, j, dim = a->row;
	complx res;

	if ( (a->row != a->col) || (b->row != b->col) || (a->row != b->col) ) {
		fprintf(stderr,"Error: cm_trace_adjoint - non-square matrix or dimension mismatch Tr{(%d,%d)+ x (%d,%d)}\n",a->row,a->col,b->row,b->col);
		exit(1);
	}

	assert(a->basis == b->basis);

	switch (a->type*100+b->type) {
		case 100*MAT_DENSE+MAT_DENSE :
			cblas_zdotc_sub(dim*dim, a->data, 1, b->data, 1, &res);
			break;
		case 100*MAT_DENSE+MAT_DENSE_DIAG : {
		    complx *va = a->data, *vb = b->data;
			res = Cnull;
			for (i=0; i<dim; i++) {
				res.re += va->re*vb->re+va->im*vb->im;
				res.im += -va->im*vb->re+va->re*vb->im;
				vb++;
				va += (dim+1);
			}
			break; }
		case 100*MAT_DENSE+MAT_SPARSE :
		case 100*MAT_DENSE+MAT_SPARSE_DIAG : {
			int n, *ic = b->icol;
			complx *va, *vb = b->data;
			res = Cnull;
			for (i=0; i<dim; i++) {
				n = b->irow[i+1] - b->irow[i];
				for (j=0; j<n; j++) {
					va = a->data + i + (*ic-1)*dim;
					res.re += va->re*vb->re+va->im*vb->im;
					res.im += -va->im*vb->re+va->re*vb->im;
					ic++;
					vb++;
				}
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_DENSE : {
			complx *va = a->data, *vb = b->data;
			res = Cnull;
			for (i=0; i<dim; i++) {
				res.re += va->re*vb->re+va->im*vb->im;
				res.im += -va->im*vb->re+va->re*vb->im;
				va++;
				vb += (dim+1);
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_DENSE_DIAG :
			cblas_zdotc_sub(dim, a->data,1,b->data,1,&res);
			break;
		case 100*MAT_DENSE_DIAG+MAT_SPARSE : {
			complx *va = a->data, vb;
			res = Cnull;
			for (i=1; i<=dim; i++) {
				vb = cm_getelem(b,i,i);
				res.re += va->re*vb.re+va->im*vb.im;
				res.im += -va->im*vb.re+va->re*vb.im;
				va++;
			}
			break; }
		case 100*MAT_DENSE_DIAG+MAT_SPARSE_DIAG : {
			complx *vb = b->data, *va;
			int *ic = b->icol;
			res = Cnull;
			j = b->irow[dim]-1;
			for (i=0; i<j; i++) {
				va = a->data + *ic - 1;
				res.re += va->re*vb->re+va->im*vb->im;
				res.im += -va->im*vb->re+va->re*vb->im;
				vb++;
				ic++;
			}
			break; }
		case 100*MAT_SPARSE+MAT_DENSE : {
			int n, *ic = a->icol;
			complx *vb, *va = a->data;
			res = Cnull;
			for (i=0; i<dim; i++) {
				n = a->irow[i+1] - a->irow[i];
				for (j=0; j<n; j++) {
					vb = b->data + i + (*ic-1)*dim;
					res.re += va->re*vb->re+va->im*vb->im;
					res.im += -va->im*vb->re+va->re*vb->im;
					ic++;
					va++;
				}
			}
			break; }
		case 100*MAT_SPARSE_DIAG+MAT_DENSE : {
			complx *va = a->data, *vb;
			int *ic = a->icol;
			res = Cnull;
			j = a->irow[dim]-1;
			for (i=0; i<j; i++) {
				vb = b->data + (*ic-1) + (*ic-1)*dim;
				res.re += va->re*vb->re+va->im*vb->im;
				res.im += -va->im*vb->re+va->re*vb->im;
				va++;
				ic++;
			}
			break; }
		case 100*MAT_SPARSE+MAT_DENSE_DIAG : {
			complx *vb = b->data, va;
			res = Cnull;
			for (i=1; i<=dim; i++) {
				va = cm_getelem(a,i,i);
				res.re += va.re*vb->re+va.im*vb->im;
				res.im += -va.im*vb->re+va.re*vb->im;
				vb++;
			}
			break; }
		case 100*MAT_SPARSE+MAT_SPARSE :
		case 100*MAT_SPARSE+MAT_SPARSE_DIAG :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE_DIAG : {
			int la, lb, len, *ic, *ir, kk;
			complx *va, vb;
			mat_complx *dum;
			la = a->irow[dim]-1;
			lb = b->irow[dim]-1;
			if (la < lb) {
				va = a->data;
				ic = a->icol;
				ir = a->irow;
				len = la;
				dum = b;
				kk = 0;
			} else {
				va = b->data;
				ic = b->icol;
				ir = b->irow;
				len = lb;
				dum = a;
				kk = 1;
			}
			res = Cnull;
			for (i=0; i<dim; i++) {
				la = ir[i+1] - ir[i];
				for (j=0; j<la; j++) {
					vb = cm_getelem(dum,i+1,*ic);
					res.re += va->re*vb.re+va->im*vb.im;
					res.im += -va->im*vb.re+va->re*vb.im;
					va++;
					ic++;
				}
			}
			if (kk) res.im *= -1.0;
			break; }
		case 100*MAT_SPARSE_DIAG+MAT_DENSE_DIAG : {
			complx *va = a->data, *vb;
			int *ic = a->icol;
			res = Cnull;
			j = a->irow[dim]-1;
			for (i=0; i<j; i++) {
				vb = b->data + *ic - 1;
				res.re += va->re*vb->re+va->im*vb->im;
				res.im += -va->im*vb->re+va->re*vb->im;
				va++;
				ic++;
			}
			break; }
	   default :
		   fprintf(stderr,"Error: cm_trace_adjoint - unknown types '%d' and '%d'\n",a->type, b->type);
		   exit(1);
	}
	return res;

}

void cm_diag(mat_complx *m, complx *eigs, mat_complx *T)
{

	/* this is called only with non-diagonal matrices */
	assert( (m->type != MAT_DENSE_DIAG) && (m->type != MAT_SPARSE_DIAG));

	if (m->row != m->col) {
		fprintf(stderr,"Error: cm_diag - non-square matrix (%d,%d)\n",m->row,m->col);
		exit(1);
	}

	mat_complx *dum;
	const int ione = 1, dim = m->row, lwsp = 3*dim;
	int info=0;
	complx *wsp = (complx*)malloc(lwsp*sizeof(complx));
	double *rwork = (double*)malloc(lwsp*sizeof(double));

	TIMING_INIT;
	TIMING_INIT_VAR(tv1);
	TIMING_INIT_VAR(tv2);
	TIMING_TIC(tv1);

	if (m->type == MAT_SPARSE) {
		printf("WARNING!!! In order to diagonalize the matrix (dim = %d) I need to \n"
			   "convert it from sparse to full format. The whole operation may take long time...\n",dim);
		dum = cm_dup(m);
		cm_dense(dum);
	} else {
		dum = cm_dup(m);
	}

	if (LEN(eigs) != dim) {
		fprintf(stderr,"Error: cm_diag - wrong dimension of eigs\n");
		exit(1);
	}
	if (T->row != dim || T->col != dim || T->type != MAT_DENSE) {
		fprintf(stderr,"Error: cm_diag - wrong dimension or type of T (%d,%d)\n",T->row,T->col);
		exit(1);
	}
	T->basis = m->basis;

	zgeev_("N","V",&dim,dum->data,&dim,&(eigs[1]),NULL,&ione,T->data,&dim,wsp,&lwsp,rwork,&info);

	if ( info != 0) {
		fprintf(stderr,"cm_diag error: diagonalization failed with the code '%d'\n", info);
	    exit(1);
	}
	free_complx_matrix(dum);
	free(wsp);
	free(rwork);

	TIMING_TOC(tv1,tv2,"cm_diag");
}

double cm_sumnorm1(mat_complx *m)
{
	double res;
	int len;

	switch (m->type) {
		case MAT_DENSE :
			len = m->row*m->col;
			break;
		case MAT_DENSE_DIAG :
			len = m->row;
			break;
		case MAT_SPARSE :
		case MAT_SPARSE_DIAG :
			len = m->irow[m->row] - 1;
			break;
	   default :
		   fprintf(stderr,"Error: cm_sumnorm1 - unknown type '%d'\n",m->type);
		   exit(1);
	}

	res = cblas_dzasum(len,m->data,1);
	return res;
}

/* norm estimation taken as largest sum of absolute values on a row */
double cm_normest(mat_complx *m)
{
	double res=-10.0, dum;
	int i,j,N;
	complx *z = m->data;

	switch (m->type) {
	case MAT_DENSE_DIAG:
		for (i=0; i<m->row; i++) {
			dum = sqrt(z->re*z->re+z->im*z->im);
			if (dum > res) res = dum;
			z++;
		}
		break;
	case MAT_SPARSE_DIAG:
		for (i=0; i<m->irow[m->row]-1; i++) {
			dum = sqrt(z->re*z->re+z->im*z->im);
			if (dum > res) res = dum;
			z++;
		}
		break;
	case MAT_DENSE:
		for (i=0; i<m->row; i++) {
			dum = 0;
			for (j=0; j<m->col; j++) {
				dum += sqrt(z->re*z->re+z->im*z->im);
				z++;
			}
			if (dum > res) res = dum;
		}
		break;
	case MAT_SPARSE:
		for (i=0; i<m->row; i++) {
			N = m->irow[i+1] - m->irow[i];
			dum = 0;
			for (j=0; j<N; j++ ) {
				dum += sqrt(z->re*z->re+z->im*z->im);
				z++;
			}
			if (dum > res) res = dum;
		}
		break;
	default:
		fprintf(stderr,"Error: cm_normest - invalid matrix type\n");
		exit(1);
	}

	return res;
}


int cm_nnz(mat_complx *m)
{
	int nnz, i;
	complx *d;

	switch (m->type) {
	case MAT_DENSE :
		nnz = 0;
		d = m->data;
		for (i=0; i<m->row*m->col; i++) {
			if ( d->re*d->re+d->im*d->im > TINY) nnz++;
			d++;
		}
		break;
	case MAT_DENSE_DIAG:
		nnz = 0;
		d = m->data;
		for (i=0; i<m->row; i++) {
			if (d->re*d->re+d->im*d->im > TINY) nnz++;
			d++;
		}
		break;
	case MAT_SPARSE:
	case MAT_SPARSE_DIAG:
		nnz = m->irow[m->row] -1;
		break;
	default:
		fprintf(stderr,"Error: cm_nnz - invalid matrix type\n");
		exit(1);
	}

	return nnz;
}

void cm_vec_on_diag(mat_complx *m, complx *vec, int basis)
{
	int dim = LEN(vec);
	mat_complx *res;

	res = complx_matrix(dim,dim,MAT_DENSE_DIAG,0,basis);
	res->basis = m->basis;
	memcpy(res->data,&(vec[1]),dim*sizeof(complx));
	cm_swap_innards_and_destroy(m,res);

}

/* only for full dense matrices */
mat_complx * cm_direct(mat_complx *a, mat_complx *b)
{
	mat_complx *res;
	int Nar, Nac, Nbr, Nbc, i, j, k, s1, s2;
	complx * tmp, *da, *db, *dr;

	if (a->type != MAT_DENSE || b->type != MAT_DENSE) {
		fprintf(stderr,"Error: cm_direct does not support other matrix types than MAT_DENSE\n");
		exit(1);
	}
	assert(a->basis == b->basis);

	  Nar = a->row;
	  Nac = a->col;
	  Nbr = b->row;
	  Nbc = b->col;
	  s1 = Nar*Nbr;
	  s2 = s1*Nbc;

	  res = complx_matrix(s1, Nac*Nbc, MAT_DENSE, 0, a->basis);
	  cm_zero(res);
	  da = a->data;
	  db = b->data;
	  dr = res->data;
	  for (i=0; i<Nbc; i++) {
	     tmp = &(db[i*Nbr]); /* this is i-column of b */
	     for (j=0; j<Nar; j++) {
	        for (k=0; k<Nac; k++) {
	        	cblas_zaxpy(Nbr, &(da[j+k*Nar]), tmp, 1, &(dr[j*Nbr+i*s1+k*s2]), 1);
	        }
	     }
	  }
	res->basis = a->basis;
	return res;
}

/****
 * sort of direct sum of two complex matrices
 *     A (+) B = ( a_11+B a_12+B ...  )
 *               ( a_21+B a_22+B ...  )
 *               (  ...               )
 *  note: this is different from direct sum and Kronecker sum in Mathematica...
 *  only for full dense matrices
 ****/
mat_complx * cm_directadd(mat_complx *a, mat_complx *b)
{
	mat_complx *tmp;
	int Nar, Nac, Nbr, Nbc, i, j, r, c, s1, ptmp;
	complx *da, *db, *dr, aval;

	if (a->type != MAT_DENSE || b->type != MAT_DENSE) {
		fprintf(stderr,"Error: cm_directadd does not support other matrix types than MAT_DENSE\n");
		exit(1);
	}
	assert(a->basis == b->basis);

	Nar = a->row;
	Nac = a->col;
	Nbr = b->row;
	Nbc = b->col;
	s1 = Nar*Nbr;

	tmp = complx_matrix(s1, Nac*Nbc, MAT_DENSE, 0, a->basis);
	cm_zero(tmp);
	da = a->data;
	db = b->data;
	dr = tmp->data;

	for (r=0; r<Nar; r++) {
		for (c=0; c<Nac; c++) {
			aval = da[r+c*Nar];
			for (i=0; i<Nbr; i++) {
				for (j=0; j<Nbc; j++) {
					ptmp = r*Nbr+c*s1*Nbc+i+j*s1;
					dr[ptmp] = Cadd(aval,db[i+j*Nbr]);
				}
			}
		}
	}
	return tmp;
}




/*********************** matrix transformations *******************/

/* Note: assumes unitary T, T-1 = T+ */
void simtrans(mat_complx *m, mat_complx *T)
{
	/* m = T m T+ */
	int dim = m->row;
	int dosparse = 0;

	//LARGE_INTEGER tv1, tv2, tickpsec;
	//QueryPerformanceFrequency(&tickpsec);
	//QueryPerformanceCounter(&tv1);

	if (dim != m->col || T->row != T->col || dim != T->row) {
		fprintf(stderr,"Error: simtrans - non-square matrix or dimension mismatch in T m T+, m=(%d,%d), T=(%d,%d)\n",m->row,m->col,T->row,T->col);
		exit(1);
	}
	assert(m->basis == T->basis);

	//DEBUGPRINT("simtrans gets types m='%s', T='%s'\n",matrix_type(m->type),matrix_type(T->type));
	switch (m->type*100+T->type) {
		case 100*MAT_DENSE+MAT_DENSE : {
			mat_complx *tmp;
			complx zone = Complx(1.0,0.0);
			tmp = complx_matrix(dim, dim, MAT_DENSE, 0,m->basis);
			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,dim,dim,dim,&zone,m->data,dim,T->data,dim,&Cnull,tmp->data,dim);
			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,&zone,T->data,dim,tmp->data,dim,&Cnull,m->data,dim);
			free_complx_matrix(tmp);
			break; }
		case 100*MAT_DENSE_DIAG+MAT_DENSE_DIAG :
		case 100*MAT_DENSE_DIAG+MAT_SPARSE_DIAG :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE_DIAG :
		case 100*MAT_SPARSE_DIAG+MAT_DENSE_DIAG :
			/* T commutes with m => m does not change */
			break;
		case 100*MAT_DENSE+MAT_DENSE_DIAG : {
			int i,j;
			complx *dA, *dA2, *dU, *dU2;
			dU = T->data;
			for (i=0; i<dim; i++) { /* columns */
				dA = dA2 = m->data + (dim+1)*i; // sit on the diagonal element
				dU2 = T->data + i;
				for (j=i+1; j<dim; j++) { /* rows, lower triangle */
					dU2++;
					dA++;
					dA2 += dim;
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
		case 100*MAT_DENSE+MAT_SPARSE_DIAG : {
			int i,j;
			complx *dA, *dA2, dU, dU2;
			for (i=0; i<dim; i++) { /* columns */
				dU = cm_getelem(T,i+1,i+1);
				dA = dA2 = m->data + (dim+1)*i; // sit on the diagonal element
				for (j=i+1; j<dim; j++) { /* rows, lower triangle */
					dU2 = cm_getelem(T,j+1,j+1);
					dA++;
					dA2 += dim;
					double re = dA->re*dU.re + dA->im*dU.im;
					double im = dA->im*dU.re - dA->re*dU.im;
					dA->re = re*dU2.re - im*dU2.im;
					dA->im = im*dU2.re + re*dU2.im;
					re = dA2->re*dU.re - dA2->im*dU.im;
					im = dA2->im*dU.re + dA2->re*dU.im;
					dA2->re = re*dU2.re + im*dU2.im;
					dA2->im = im*dU2.re - re*dU2.im;
				}
			}
			break; }
		case 100*MAT_SPARSE+MAT_DENSE_DIAG :
		case 100*MAT_SPARSE+MAT_SPARSE_DIAG : {
			int i, j, n, *ic = m->icol;
			complx dU, dU2, *dA = m->data;
			for (i=0; i<dim; i++) {
				dU = cm_getelem(T,i+1,i+1);
				n = m->irow[i+1] - m->irow[i];
				for (j=0; j<n; j++) {
					dU2 = cm_getelem(T, *ic, *ic);
					double re = dA->re*dU.re - dA->im*dU.im;
					double im = dA->im*dU.re + dA->re*dU.im;
					dA->re = re*dU2.re + im*dU2.im;
					dA->im = im*dU2.re - re*dU2.im;
					dA++;
					ic++;
				}
			}
			dosparse = 1;
			break; }
		case 100*MAT_SPARSE+MAT_DENSE : {
			/*** this seems slower...
			mat_complx *tmp = cm_adjoint(T);
			mat_complx *dum = simpson_zcsrmm(m,tmp);
			assert(dum->type == MAT_DENSE);
			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,&Cunit,T->data,dim,dum->data,dim,&Cnull,tmp->data,dim);
			free_complx_matrix(dum);
			cm_swap_innards_and_destroy(m,tmp);
			*** original code faster ***/
			mat_complx *tmp, *dum;
			complx zone = Complx(1.0,0.0);
			tmp = cm_mul(T,m); assert(tmp->type == MAT_DENSE);
			dum = complx_matrix(dim, dim, MAT_DENSE,0,m->basis);
			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,dim,dim,dim,&zone,tmp->data,dim,T->data,dim,&Cnull,dum->data,dim);
			free_complx_matrix(tmp);
			cm_swap_innards_and_destroy(m,dum);
			/************/
			break; }
		case 100*MAT_DENSE+MAT_SPARSE : {
			mat_complx *tmp, *tmpA, *dum;
			tmp = simpson_zcsrmm(T,m); assert(tmp->type == MAT_DENSE);
			tmpA = cm_adjoint(T);
			dum = cm_mul(tmp,tmpA);
			free_complx_matrix(tmp);
			free_complx_matrix(tmpA);
			cm_swap_innards_and_destroy(m,dum);
			break; }
		case 100*MAT_SPARSE_DIAG+MAT_DENSE :
		case 100*MAT_DENSE_DIAG+MAT_DENSE : {
			mat_complx *tmp, *tmp2;
			complx zone = Complx(1.0,0.0);
			tmp = cm_mul(T,m); assert(tmp->type == MAT_DENSE);
			tmp2 = complx_matrix(dim,dim,MAT_DENSE,0,m->basis);
			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,dim,dim,dim,&zone,tmp->data,dim,T->data,dim,&Cnull,tmp2->data,dim);
			free_complx_matrix(tmp);
			cm_swap_innards_and_destroy(m,tmp2);
			break; }
		case 100*MAT_DENSE_DIAG+MAT_SPARSE :
		case 100*MAT_SPARSE+MAT_SPARSE :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE : {
			mat_complx *tmp;
			//DEBUGPRINT("   simtrans sparse - point 1\n");
			tmp = cm_adjoint(T);
			//DEBUGPRINT("   simtrans sparse - point 2\n");
			//cm_print(T,"matrix T"); cm_print(tmp,"adjoint T");
			cm_multo(m,tmp);
			//DEBUGPRINT("   simtrans sparse - point 3\n");
			cm_multo_rev(m,T);
			//DEBUGPRINT("   simtrans sparse - point 4\n");
			free_complx_matrix(tmp);
			dosparse = 0;
			break; }
	   default :
			   fprintf(stderr,"Error: simtrans - unknown types '%d' and '%d'\n",m->type, T->type);
			   exit(1);
	}
	if (dosparse) cm_sparse(m,SPARSE_TOL);

	//QueryPerformanceCounter(&tv2);
	//printf("timing simtrans: %.9f\n",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);

	return;
}

/* Note: assumes unitary T, T-1 = T+ */
void simtrans_adj(mat_complx *m, mat_complx *T)
{
	/* m = T+ m T */
	//printf("simtrans_adj - matrix types '%d' and '%d'\n",m->type, T->type);
	int dim = m->row;
	int dosparse = 0;

	//LARGE_INTEGER tv1, tv2, tickpsec, tv3, tv4;
	//QueryPerformanceFrequency(&tickpsec);
	//QueryPerformanceCounter(&tv1);

	if (dim != m->col || T->row != T->col || dim != T->row) {
		fprintf(stderr,"Error: simtrans_adj - non-square matrix or dimension mismatch in T+ m T, m=(%d,%d), T=(%d,%d)\n",m->row,m->col,T->row,T->col);
		exit(1);
	}
	assert(m->basis == T->basis);

	switch (m->type*100+T->type) {
		case 100*MAT_DENSE+MAT_DENSE : {
			mat_complx *tmp;
			//complx zone = Complx(1.0,0.0);
			complx zone = {1.0, 0.0};
			tmp = complx_matrix(dim, dim, MAT_DENSE, 0,m->basis);
			//QueryPerformanceCounter(&tv3);
			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,&zone,m->data,dim,T->data,dim,&Cnull,tmp->data,dim);
			//QueryPerformanceCounter(&tv4);
			//printf("\t\ttiming zgemm 1 in simtrans_adj: %.9f\n",((float)(tv4.QuadPart-tv3.QuadPart))/(float)tickpsec.QuadPart);
			cblas_zgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,dim,dim,dim,&zone,T->data,dim,tmp->data,dim,&Cnull,m->data,dim);
			//QueryPerformanceCounter(&tv3);
			//printf("timing simtrans_adj time: %.9f\n",((float)(tv3.QuadPart-tv4.QuadPart))/(float)tickpsec.QuadPart);
			free_complx_matrix(tmp);
			break; }
		case 100*MAT_DENSE_DIAG+MAT_DENSE_DIAG :
		case 100*MAT_DENSE_DIAG+MAT_SPARSE_DIAG :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE_DIAG :
		case 100*MAT_SPARSE_DIAG+MAT_DENSE_DIAG :
			/* T commutes with m => m does not change */
			break;
		case 100*MAT_DENSE+MAT_DENSE_DIAG : {
			int i,j;
			complx *dA, *dA2, *dU, *dU2;
			dU = T->data;
			for (i=0; i<dim; i++) { /* columns */
				dA = dA2 = m->data + (dim+1)*i; // sit on the diagonal element
				dU2 = T->data + i;
				for (j=i+1; j<dim; j++) { /* rows, lower triangle */
					dU2++;
					dA++;
					dA2 += dim;
					double re = dA->re*dU->re - dA->im*dU->im;
					double im = dA->im*dU->re + dA->re*dU->im;
					dA->re = re*dU2->re + im*dU2->im;
					dA->im = im*dU2->re - re*dU2->im;
					re = dA2->re*dU->re + dA2->im*dU->im;
					im = dA2->im*dU->re - dA2->re*dU->im;
					dA2->re = re*dU2->re - im*dU2->im;
					dA2->im = im*dU2->re + re*dU2->im;
				}
				dU++;
			}
			break; }
		case 100*MAT_DENSE+MAT_SPARSE_DIAG : {
			int i,j;
			complx *dA, *dA2, dU, dU2;
			for (i=0; i<dim; i++) { /* columns */
				dU = cm_getelem(T,i+1,i+1);
				dA = dA2 = m->data + (dim+1)*i; // sit on the diagonal element
				for (j=i+1; j<dim; j++) { /* rows, lower triangle */
					dU2 = cm_getelem(T,j+1,j+1);
					dA++;
					dA2 += dim;
					double re = dA->re*dU.re - dA->im*dU.im;
					double im = dA->im*dU.re + dA->re*dU.im;
					dA->re = re*dU2.re + im*dU2.im;
					dA->im = im*dU2.re - re*dU2.im;
					re = dA2->re*dU.re + dA2->im*dU.im;
					im = dA2->im*dU.re - dA2->re*dU.im;
					dA2->re = re*dU2.re - im*dU2.im;
					dA2->im = im*dU2.re + re*dU2.im;
				}
			}
			break; }
		case 100*MAT_SPARSE+MAT_DENSE_DIAG :
		case 100*MAT_SPARSE+MAT_SPARSE_DIAG : {
			int i, j, n, *ic = m->icol;
			complx dU, dU2, *dA = m->data;
			for (i=0; i<dim; i++) {
				dU = cm_getelem(T,i+1,i+1);
				n = m->irow[i+1] - m->irow[i];
				for (j=0; j<n; j++) {
					dU2 = cm_getelem(T, *ic, *ic);
					double re = dA->re*dU.re + dA->im*dU.im;
					double im = dA->im*dU.re - dA->re*dU.im;
					dA->re = re*dU2.re - im*dU2.im;
					dA->im = im*dU2.re + re*dU2.im;
					dA++;
					ic++;
				}
			}
			dosparse = 1;
			break; }
		case 100*MAT_SPARSE+MAT_DENSE : {
			mat_complx *tmp, *dum;
			complx zone = {1.0,0.0};
			tmp = simpson_zcsrmm(m,T); assert(tmp->type == MAT_DENSE);
			dum = complx_matrix(dim, dim, MAT_DENSE,0,m->basis);
			cblas_zgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,dim,dim,dim,&zone,T->data,dim,tmp->data,dim,&Cnull,dum->data,dim);
			free_complx_matrix(tmp);
			cm_swap_innards_and_destroy(m,dum);
			break; }
		case 100*MAT_DENSE+MAT_SPARSE : {
			mat_complx *tmp, *tmpA, *dum;
			tmpA = cm_adjoint(T);
			tmp = simpson_zcsrmm(tmpA,m);
			dum = cm_mul(tmp,T);
			free_complx_matrix(tmp);
			free_complx_matrix(tmpA);
			cm_swap_innards_and_destroy(m,dum);
			break; }
		case 100*MAT_SPARSE_DIAG+MAT_DENSE :
		case 100*MAT_DENSE_DIAG+MAT_DENSE : {
			mat_complx *tmp, *dum;
			complx zone = Complx(1.0,0.0);
			tmp = cm_mul(m,T); assert(tmp->type == MAT_DENSE);
			dum = complx_matrix(dim, dim, MAT_DENSE,0,m->basis);
			cblas_zgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,dim,dim,dim,&zone,T->data,dim,tmp->data,dim,&Cnull,dum->data,dim);
			free_complx_matrix(tmp);
			cm_swap_innards_and_destroy(m,dum);
			break; }
		case 100*MAT_DENSE_DIAG+MAT_SPARSE :
		case 100*MAT_SPARSE+MAT_SPARSE :
		case 100*MAT_SPARSE_DIAG+MAT_SPARSE : {
			mat_complx *tmp;
			tmp = cm_adjoint(T);
			cm_multo(m,T);
			cm_multo_rev(m,tmp);
			free_complx_matrix(tmp);
			dosparse = 1;
			break; }
		default :
			fprintf(stderr,"Error: simtrans_adj - unknown types '%d' and '%d'\n",m->type, T->type);
			exit(1);
	}
	if (dosparse) cm_sparse(m,SPARSE_TOL);

	//QueryPerformanceCounter(&tv2);
	//printf("timing simtrans_adj time: %.9f\n",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);

	return;
}

void simtrans_zrot2(mat_complx *m, mat_double *T)
{
	/* exp(-i T) m exp(i T)  */
	assert( (T->type == MAT_DENSE_DIAG) || (T->type == MAT_SPARSE_DIAG));
	assert( T->row == m->row );
	if (m->row != m->col) {
		fprintf(stderr,"Error: simtrans_zrot2 - non-square matrix m (%d,%d)\n",m->row, m->col);
		exit(1);
	}
	assert(m->basis == T->basis);

// tip for optimization: take sin and cos functions out of the double loop,
//	multiply complex numbers there
	switch (m->type) {
		case MAT_DENSE_DIAG :
		case MAT_SPARSE_DIAG :
			/* no change here */
			break;
		case MAT_DENSE : {
			int i,j,N=m->row;
			complx *expT, *z, *zz, *dU, *dU2;
			double *dd, rr, ii, c, s;
			z = expT = (complx*)malloc(N*sizeof(complx));
			dd = T->data;
			if (T->type == MAT_DENSE_DIAG) {
				for(i=0; i<N; i++) {
					z->re = cos(*dd);
					z->im = sin(*dd);
					z++;
					dd++;
				}
			} else {
				/* T->type == MAT_SPARSE_DIAG */
				/* strictly assumes only diagonal elements are present! */
				for (i=1; i<=N; i++) {
					int nc = T->irow[i+1] - T->irow[i];
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
			z = expT;
			for (i=0; i<N; i++) { /* columns */
				dU = dU2 = m->data + (N+1)*i; // sit on the diagonal
				zz = expT + i; //dP2 = T->data + i;
				for (j=i+1; j<N; j++) { /* rows, lower triangle */
					dU++;
					dU2 += N;
					zz++; //dP2++;
					//ph = *dP - *dP2;
					c = z->re*zz->re + z->im*zz->im;
					s = z->im*zz->re - z->re*zz->im;
					//if (fabs(ph)<TINY) continue;
					//c = cos(ph);
					//s = sin(ph);
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
			free(expT);
			break; }
		case MAT_SPARSE : {
			int i, j, n, *ic=m->icol;
			complx *z=m->data;
			double ph, re, s, c;
			for (i=0; i<m->row; i++) {
				n = m->irow[i+1] - m->irow[i];
				for (j=0; j<n; j++) {
					ph = dm_getelem(T,i+1,i+1) - dm_getelem(T,*ic,*ic);
					re = z->re;
					s = sin(ph);
					c = cos(ph);
					z->re = re*c + z->im*s;
					z->im = z->im*c - re*s;
					z++;
					ic++;
				}
			}
			break; }
		default :
			fprintf(stderr,"Error: simtrans_zrot2 - unknown type '%d'\n",m->type);
			exit(1);
	}
}

/******************** propagators and exponentials *****************/

void prop_pade_real(mat_complx *prop, mat_double *ham, double dt)
{
	//int m = 0;
	int nsq = 0;
	double scl = 1.0, norm;
	mat_double *U, *V;
	const double b[] = {64764752532480000, 32382376266240000, 7771770303897600,
						1187353796428800, 129060195264000, 10559470521600,
						670442572800, 33522128640, 1323241920, 40840800, 960960,
						16380, 182, 1};
	int dim = ham->row;

	TIMING_INIT_VAR2(tv1,tv2);
	TIMING_TIC(tv1);

	norm = dm_normest(ham)*dt; // NOTE: here 1-norm must be used!!!
	// I use infinity-norm since it is easier to code for sparse matrices
	// ham is real symmetric and therefore the result is identical
	printf("prop_pade_real norm = %f\n",norm);
	if (norm < 1.495585217958292e-2) {
		//m = 3;
		U = dm_creatediag('d',dim,b[1],ham->basis);
		V = dm_creatediag('d',dim,b[0],ham->basis);
		dm_muld(ham, dt);
		mat_double *H2 = dm_mul(ham,ham);
		dm_multod(U,H2,-b[3]);
		dm_multod(V,H2,-b[2]);
		free_double_matrix(H2);
	} else if (norm < 0.253939833006323) {
		//m = 5;
		U = dm_creatediag('d',dim,b[1],ham->basis);
		V = dm_creatediag('d',dim,b[0],ham->basis);
		dm_muld(ham, dt);
		mat_double *H2 = dm_mul(ham,ham);
		dm_multod(U,H2,-b[3]);
		dm_multod(V,H2,-b[2]);
		mat_double *H4 = dm_mul(H2,H2);
		free_double_matrix(H2);
		dm_multod(U,H4,b[5]);
		dm_multod(V, H4,b[4]);
		free_double_matrix(H4);
	} else if (norm < 0.9504178996162932) {
		//m = 7;
		U = dm_creatediag('d',dim,b[1],ham->basis);
		V = dm_creatediag('d',dim,b[0],ham->basis);
		dm_muld(ham, dt);
		mat_double *H2 = dm_mul(ham,ham);
		dm_multod(U,H2,-b[3]);
		dm_multod(V,H2,-b[2]);
		mat_double *H4 = dm_mul(H2,H2);
		mat_double *H6 = dm_mul(H2,H4);
		free_double_matrix(H2);
		dm_multod(U,H4,b[5]);
		dm_multod(V, H4,b[4]);
		free_double_matrix(H4);
		dm_multod(U,H6,-b[7]);
		dm_multod(V, H6,-b[6]);
		free_double_matrix(H6);
	} else if (norm < 2.097847961257068) {
		//m = 9;
		U = dm_creatediag('d',dim,b[1],ham->basis);
		V = dm_creatediag('d',dim,b[0],ham->basis);
		dm_muld(ham, dt);
		mat_double *H2 = dm_mul(ham,ham);
		dm_multod(U,H2,-b[3]);
		dm_multod(V,H2,-b[2]);
		mat_double *H4 = dm_mul(H2,H2);
		dm_multod(U,H4,b[5]);
		dm_multod(V,H4,b[4]);
		mat_double *H6 = dm_mul(H2,H4);
		free_double_matrix(H4);
		dm_multod(U,H6,-b[7]);
		dm_multod(V,H6,-b[6]);
		mat_double *H8 = dm_mul(H2,H6);
		free_double_matrix(H6);
		free_double_matrix(H2);
		dm_multod(U,H8,b[9]);
		dm_multod(V,H8,b[8]);
		free_double_matrix(H8);
	} else {
		//m = 13;
		if (norm > 5.371920351148152) {
			// need scaling
			nsq = (int)ceil(log(norm/5.371920351148152)/log(2.0));
			scl = pow(2,nsq);
		}
		printf("prop_pade_real nsq = %d, scl = %f\n",nsq,scl);
		dm_muld(ham, dt/scl);
		mat_double *H2 = dm_mul(ham, ham);
		mat_double *H4 = dm_mul(H2, H2);
		mat_double *H6 = dm_mul(H2, H4);
		U = double_matrix(dim,dim,MAT_DENSE,0,ham->basis);
		dm_zero(U);
		dm_multod(U,H2, b[9]);
		dm_multod(U,H4,-b[11]);
		dm_multod(U,H6, b[13]);
		dm_multo_rev(U,H6);
		dm_multod(U,H6,-b[7]);
		dm_multod(U,H4, b[5]);
		dm_multod(U,H2,-b[3]);
		dm_addtodiag(U,b[1]);
		V = double_matrix(dim,dim,MAT_DENSE,0,ham->basis);
		dm_zero(V);
		dm_multod(V,H2, b[8]);
		dm_multod(V,H4,-b[10]);
		dm_multod(V,H6, b[12]);
		dm_multo_rev(V,H6);
		dm_multod(V,H6,-b[6]);
		dm_multod(V,H4, b[4]);
		dm_multod(V,H2,-b[2]);
		dm_addtodiag(V,b[0]);
		free_double_matrix(H2);
		free_double_matrix(H4);
		free_double_matrix(H6);
	}
	dm_multo_rev(U,ham);
	mat_complx *lhs = dm_complx(V);
	//dm_copy2cm(V,prop);
	free_double_matrix(V);
	cm_dense(lhs);
	//cm_dense(prop);
	mat_complx *rhs = cm_dup(lhs);
	//assert( (lhs->type == MAT_DENSE) && (prop->type == MAT_DENSE) );
	dm_dense(U);
	assert(U->type == MAT_DENSE);
	cblas_daxpy(dim*dim,1.0,U->data,1,(double*)(lhs->data)+1,2);
	//cblas_daxpy(dim*dim,-1.0,U->data,1,(double*)(prop->data)+1,2);
	cblas_daxpy(dim*dim,-1.0,U->data,1,(double*)(rhs->data)+1,2);
	TIMING_TOC(tv1,tv2,"Pade after polynomials");
	int *pvec = (int*)malloc((dim+1)*sizeof(int));
	int info = 0;
	//zgesv_(&dim,&dim,lhs->data,&dim,pvec,prop->data,&dim,&info);
	zgesv_(&dim,&dim,lhs->data,&dim,pvec,rhs->data,&dim,&info);
	// tried zsysv solver for symmetric matrix but it is again slower...
	if (info != 0) {
		fprintf(stderr,"prop_pade error: zgesv failed with info=%d\n",info);
		exit(1);
	}
    free(pvec);
    free_complx_matrix(lhs);
	TIMING_TOC(tv1,tv2,"Pade after ratio solver");
	/* squaring */
	if (nsq) {
		int k;
		for (k=0; k<nsq; k++) {
			//cm_multo(prop,prop);
			cm_multo(rhs,rhs);
		}
	}
	TIMING_TOC(tv1,tv2,"Pade after squarings");
	// update propagator
	cm_multo_rev(prop,rhs);
	free_complx_matrix(rhs);
	TIMING_TOC(tv1,tv2,"Pade after prop update");

}

void dm_lams(mat_double *ham, double *lmax, double *lmin)
{
	int dim = ham->row;
	int i;
	double l1, l2;
	double *dd = ham->data;
	*lmax = -1.0e99;
	*lmin =  1.0e99;

	switch (ham->type) {
	case MAT_DENSE: {
		for(i=0; i<dim; i++) {
			double d2 = cblas_dasum(dim,dd,1); // adds numbers in a column
			double d1 = dd[i];
			if (d1<0) {
				l1 = 2.0*d1+d2;
				l2 = -d2;
			} else {
				l1 = d2;
				l2 = 2.0*d1-d2;
			}
			if (l1 > *lmax) *lmax = l1;
			if (l2 < *lmin) *lmin = l2;
			dd += dim;
		}
		break; }
	case MAT_SPARSE: {
		for (i=0; i< dim; i++) {
			int len = ham->irow[i+1] - ham->irow[i];
			if (len == 0) {
				l1 = l2 = 0.0;
			} else {
				double d2 = cblas_dasum(len,dd,1); // adds numbers in a row
				double d1 = dm_getelem(ham,i+1,i+1);
				if (d1<0) {
					l1 = 2.0*d1+d2;
					l2 = -d2;
				} else {
					l1 = d2;
					l2 = 2.0*d1-d2;
				}
			}
			if (l1 > *lmax) *lmax = l1;
			if (l2 < *lmin) *lmin = l2;
			dd += len;
		}
		break; }
	default:
		fprintf(stderr,"dm_lams error: matrix type %d not supported\n",ham->type);
		exit(1);
	}

	assert(*lmax > *lmin);
}

/* prop is updated with exp(-i*ham*dt), i.e. prop = exp(-i*ham*dt)*prop */
/*      (ham gets overwritten)   */
void prop_cheby_real(mat_complx *prop, mat_double *ham, double dt)
{
	double TOL = 1e-6;
	double lmax = 0.0;
	double lmin = 0.0;

	//printf("\tprop_cheby IN: ham is %s, prop is %s\n",matrix_type(ham->type),matrix_type(prop->type));
	dm_muld(ham,-dt);
	//dm_print(ham,"HAMILTONIAN");
	//printf("\t--> ham norm is %g",dm_normest(ham));
	dm_lams(ham,&lmax, &lmin);
	double dw = 0.5*(lmax-lmin);
	double wb = 0.5*(lmax+lmin);
	dm_addtodiag(ham,-wb);
	dm_muld(ham,1.0/dw);
	//printf(", after scaling is %g (%g, %g)\n",dm_normest(ham),wb,dw);
	cm_mulc(prop,Cexpi(wb));

	mat_double *Ure = cm_real(prop);
	mat_double *Uim = cm_imag(prop);
	double aa = bessj(0,dw);
	dm_muld(Ure,aa);
	dm_muld(Uim,aa);

	mat_double *T0re = cm_real(prop);
	mat_double *T0im = cm_imag(prop);
	mat_double *T1re = dm_mul(ham,T0re);
	mat_double *T1im = dm_mul(ham,T0im);
	aa = 2.0*bessj(1,dw);
	dm_multod(Uim,T1re,aa);
	dm_multod(Ure,T1im,-aa);
	double norm = (fabs(dm_normest(T1re)) + fabs(dm_normest(T1im)))*fabs(aa);
	//double norm = dm_dm_normest(T1re,T1im)*fabs(aa);
	int iter = 1;
	mat_double *p0re = T0re;
	mat_double *p0im = T0im;
	mat_double *p1re = T1re;
	mat_double *p1im = T1im;
	mat_double *pdum;
	// fast convergence guaranteed for iter > dw
	while ( (iter < dw) || (norm > TOL) ) {
		iter++;
		aa = 2.0*bessj(iter,dw);
		//dm_dense(p1re); dm_dense(p1im);
		dm_mm(2.0,ham,p1re,-1.0,p0re);
		dm_mm(2.0,ham,p1im,-1.0,p0im);
		pdum = p0re; p0re = p1re; p1re = pdum;
		pdum = p0im; p0im = p1im; p1im = pdum;
		switch (iter % 4) {
		case 0:
			dm_multod(Ure,p1re,aa);
			dm_multod(Uim,p1im,aa);
			break;
		case 1:
			dm_multod(Uim,p1re,aa);
			dm_multod(Ure,p1im,-aa);
			break;
		case 2:
			dm_multod(Ure,p1re,-aa);
			dm_multod(Uim,p1im,-aa);
			break;
		case 3:
			dm_multod(Uim,p1re,-aa);
			dm_multod(Ure,p1im, aa);
		}
		norm = (fabs(dm_normest(p1re)) + fabs(dm_normest(p1im)))*fabs(aa);
		//norm = dm_dm_normest(p1re,p1im)*fabs(aa);
		//printf("\t iter %d.) T_norm = %g, Ji = %g, norm = %g\n",iter,fabs(dm_normest(p1re)) + fabs(dm_normest(p1im)),fabs(aa),norm);
		//printf("Real part maximum: %g\n",p1re->data[cblas_idamax(p1re->row*p1re->col,p1re->data,1)]);
		//printf("Imag part maximum: %g\n",p1im->data[cblas_idamax(p1im->row*p1im->col,p1im->data,1)]);
	}
	//printf("prop_cheby_real done in %d iterations\n",iter);

	free_double_matrix(T0re);
	free_double_matrix(T0im);
	free_double_matrix(T1im);
	free_double_matrix(T1re);

	dm_copy2cm(Ure,prop);
	free_double_matrix(Ure);
	cm_multocr(prop,Uim,Complx(0,1));
	free_double_matrix(Uim);

	//printf("\tprop_cheby OUT: ham is %s, prop is %s\n",matrix_type(ham->type),matrix_type(prop->type));

}

/* Chebyshev expm with Ham possibly sparse and other dense, with Scaling & Squaring */
void prop_cheby2_real(mat_complx *prop, mat_double *ham, double dt)
{
	double TOL = 1e-6;
	double dw=0.5;
	int nsq = 0;
	int dim = ham->row;
	int len = dim*dim;
	int i;

	//printf("ham is %s, nnz = %d\n",matrix_type(ham->type),dm_nnz(ham));

	TIMING_INIT_VAR2(tv1,tv2);
	TIMING_TIC(tv1);

	//dm_print(ham,"HAMILTONIAN");
	/* scaling step */
	double scaling = dm_normest(ham)*dt;
	//printf("prop_cheby2_real: Hamiltonian norm estimate = %g",scaling);
	if (scaling > 1) {
		nsq = (int)ceil(log(scaling)/log(2.0));
		scaling = pow(2,nsq);
	} else {
		scaling = 1.0;
	}
	//printf(" --> will do %d squarings\n",nsq);
	dm_muld(ham,-dt/scaling/dw);


	mat_double *T0 = double_matrix(dim,dim,MAT_DENSE,0,ham->basis);
	dm_zero(T0);
	mat_double *T1 = dm_dup2(ham);
	mat_complx *U = complx_matrix(dim,dim,MAT_DENSE,0,ham->basis);
	cm_zero(U);
	double aa = bessj(0,dw);
	for (i=0; i<len; i+=dim+1) {
		T0->data[i] = 1.0;
		U->data[i].re = aa;
	}

	aa = 2.0*bessj(1,dw);
	cblas_daxpy(len,aa,T1->data,1,(double*)(U->data)+1,2);

	double norm = fabs(dm_max(T1)*aa);
	int iter = 1;
	mat_double *p0 = T0;
	mat_double *p1 = T1;
	mat_double *pdum;
	while ( norm > TOL ) {
		iter++;
		if (iter > 20) {
			fprintf(stderr,"ERROR: prop_cheby2_real - exceeded max iter 20\n");
			exit(1);
		}
		aa = 2.0*bessj(iter,dw);
		dm_mm(2.0,ham,p1,-1.0,p0); // this is faster!!!
		//cblas_dsymm(CblasColMajor,CblasLeft,CblasUpper,dim,dim,2.0,ham->data,dim,p1->data,dim,-1.0,p0->data,dim);
		pdum = p0; p0 = p1; p1 = pdum;
		switch (iter % 4) {
		case 0:
			cblas_daxpy(len,aa,p1->data,1,(double*)(U->data),2);
			break;
		case 1:
			cblas_daxpy(len,aa,p1->data,1,(double*)(U->data)+1,2);
			break;
		case 2:
			cblas_daxpy(len,-aa,p1->data,1,(double*)(U->data),2);
			break;
		case 3:
			cblas_daxpy(len,-aa,p1->data,1,(double*)(U->data)+1,2);
		}
		norm = fabs(dm_max(p1)*aa);
		//printf("\t iter %d.) T_norm = %g, Ji = %g, norm = %g\n",iter,fabs(dm_max(p1)),aa,norm);
	}
	//printf("prop_cheby_real done in %d iterations\n",iter);

	free_double_matrix(T0);
	free_double_matrix(T1);

	TIMING_TOC(tv1,tv2,"chebyshev after iterations");

	// squaring
	if (nsq) {
		for (i=0; i<nsq; i++) {
			cm_multo(U,U);
		}
	}

	TIMING_TOC(tv1,tv2,"chebyshev after squarings");

	// update propagator
	cm_multo_rev(prop,U);
	free_complx_matrix(U);
	TIMING_TOC(tv1,tv2,"cheby after prop update");
}

/* Chebyshev expm with Ham possibly sparse and other dense,  shift&scale (no squaring!) */
void prop_cheby3_real(mat_complx *prop, mat_double *ham, double dt)
{
	double TOL = 1e-6;
	double lmax = 0.0;
	double lmin = 0.0;
	int dim = ham->row;
	int len = dim*dim;
	int i;

	//printf("ham is %s, nnz = %d\n",matrix_type(ham->type),dm_nnz(ham));

	TIMING_INIT_VAR2(tv1,tv2);
	TIMING_TIC(tv1);

	dm_muld(ham,-dt);
	//dm_print(ham,"HAMILTONIAN");
	//printf("\t--> ham norm is %g",dm_normest(ham));
	dm_lams(ham,&lmax, &lmin);
	double dw = 0.5*(lmax-lmin);
	double wb = 0.5*(lmax+lmin);
	dm_addtodiag(ham,-wb);
	dm_muld(ham,1.0/dw);
	//printf(", after scaling is %g (%g, %g)\n",dm_normest(ham),wb,dw);
	cm_mulc(prop,Cexpi(wb));

	TIMING_TOC(tv1,tv2,"chebyshev3 after preparations");

	mat_double *T0 = double_matrix(dim,dim,MAT_DENSE,0,ham->basis);
	dm_zero(T0);
	mat_double *T1 = dm_dup2(ham);
	mat_complx *U = complx_matrix(dim,dim,MAT_DENSE,0,ham->basis);
	cm_zero(U);
	double aa = bessj(0,dw);
	for (i=0; i<len; i+=dim+1) {
		T0->data[i] = 1.0;
		U->data[i].re = aa;
	}

	aa = 2.0*bessj(1,dw);
	cblas_daxpy(len,aa,T1->data,1,(double*)(U->data)+1,2);

	double norm = fabs(dm_max(T1)*aa);
	int iter = 1;
	mat_double *p0 = T0;
	mat_double *p1 = T1;
	mat_double *pdum;
	while ( norm > TOL ) {
		iter++;
		if (iter > 20) {
			fprintf(stderr,"ERROR: prop_cheby3_real - exceeded max iter 20\n");
			exit(1);
		}
		aa = 2.0*bessj(iter,dw);
		dm_mm(2.0,ham,p1,-1.0,p0); // this is faster!!!
		//cblas_dsymm(CblasColMajor,CblasLeft,CblasUpper,dim,dim,2.0,ham->data,dim,p1->data,dim,-1.0,p0->data,dim);
		pdum = p0; p0 = p1; p1 = pdum;
		switch (iter % 4) {
		case 0:
			cblas_daxpy(len,aa,p1->data,1,(double*)(U->data),2);
			break;
		case 1:
			cblas_daxpy(len,aa,p1->data,1,(double*)(U->data)+1,2);
			break;
		case 2:
			cblas_daxpy(len,-aa,p1->data,1,(double*)(U->data),2);
			break;
		case 3:
			cblas_daxpy(len,-aa,p1->data,1,(double*)(U->data)+1,2);
		}
		norm = fabs(dm_max(p1)*aa);
		//printf("\t iter %d.) T_norm = %g, Ji = %g, norm = %g\n",iter,fabs(dm_max(p1)),aa,norm);
	}
	//printf("prop_cheby_real done in %d iterations\n",iter);

	free_double_matrix(T0);
	free_double_matrix(T1);

	TIMING_TOC(tv1,tv2,"chebyshev3 after iterations");

	// update propagator
	cm_multo_rev(prop,U);
	free_complx_matrix(U);
	TIMING_TOC(tv1,tv2,"cheby3 after prop update");
}

/* Chebyshev expm with all matrices possibly sparse, no scaling&squaring */
void prop_cheby4_real(mat_complx *prop, mat_double *ham, double dt)
{
	double TOL = 1e-6;
	double lmax = 0.0;
	double lmin = 0.0;
	int dim = ham->row;
	//int len = dim*dim;
	//int i;

	TIMING_INIT_VAR2(tv1,tv2);
	TIMING_TIC(tv1);

	//printf("\tprop_cheby IN: ham is %s, prop is %s\n",matrix_type(ham->type),matrix_type(prop->type));
	dm_muld(ham,-dt);
	//dm_print(ham,"HAMILTONIAN");
	//printf("\t--> ham norm is %g",dm_normest(ham));
	dm_lams(ham,&lmax, &lmin);
	double dw = 0.5*(lmax-lmin);
	double wb = 0.5*(lmax+lmin);
	dm_addtodiag(ham,-wb);
	dm_muld(ham,1.0/dw);
	//printf(", after scaling is %g (%g, %g)\n",dm_normest(ham),wb,dw);
	cm_mulc(prop,Cexpi(wb));

	TIMING_TOC(tv1,tv2,"chebyshev3 preparations");

	mat_double *T0 = dm_creatediag('s',dim,1.0,ham->basis);
	mat_double *T1 = dm_dup(ham);
	double aa = bessj(0,dw);
	mat_double *ure = dm_creatediag('s',dim,aa,ham->basis);
	mat_double *uim = dm_dup(T1);
	aa = 2.0*bessj(1,dw);
	dm_muld(uim,aa);

	double norm = fabs(dm_max(T1)*aa);
	int iter = 1;
	mat_double *p0 = T0;
	mat_double *p1 = T1;
	mat_double *pdum;
	while ( norm > TOL ) {
		iter++;
		if (iter > 20) {
			fprintf(stderr,"ERROR: prop_cheby3_real - exceeded max iter 20\n");
			exit(1);
		}
		aa = 2.0*bessj(iter,dw);
		dm_mm(2.0,ham,p1,-1.0,p0); // this is faster!!!
		pdum = p0; p0 = p1; p1 = pdum;
		switch (iter % 4) {
		case 0:
			//cblas_daxpy(len,aa,p1->data,1,(double*)(U->data),2);
			dm_multod(ure,p1,aa);
			break;
		case 1:
			//cblas_daxpy(len,aa,p1->data,1,(double*)(U->data)+1,2);
			dm_multod(uim,p1,aa);
			break;
		case 2:
			//cblas_daxpy(len,-aa,p1->data,1,(double*)(U->data),2);
			dm_multod(ure,p1,-aa);
			break;
		case 3:
			//cblas_daxpy(len,-aa,p1->data,1,(double*)(U->data)+1,2);
			dm_multod(uim,p1,-aa);
		}
		norm = fabs(dm_max(p1)*aa);
		//printf("\t iter %d.) T_norm = %g, Ji = %g, norm = %g\n",iter,fabs(dm_max(p1)),aa,norm);
	}
	//printf("prop_cheby_real done in %d iterations\n",iter);

	free_double_matrix(T0);
	free_double_matrix(T1);

	TIMING_TOC(tv2,tv1,"chebyshev4 iterations");

	// update propagator
	mat_complx *U = dm_complx(ure);
	free_double_matrix(ure);
	mat_complx *dumc = dm_imag(uim);
	free_double_matrix(uim);
	cm_addto(U,dumc);
	free_complx_matrix(dumc);
	cm_multo_rev(prop,U);
	free_complx_matrix(U);
}

/* exp via Taylor expansion, Ham possibly sparse and other dense, Scaling & Squaring */
void prop_taylor_real(mat_complx *prop, mat_double *ham, double dt)
{
	const double norm_tol = 1.0e-6;
	const double dim = ham->row;
	mat_double *next_term;
	double next_term_norm, fac, scaling;
	int nsq, i;
	const int len = dim*dim;

	//printf("ham is %s, nnz = %d\n",matrix_type(ham->type),dm_nnz(ham));

	TIMING_INIT_VAR2(tv1,tv2);
	TIMING_TIC(tv1);

	/* scaling step */
	scaling = dm_normest(ham)*dt;
	//printf("prop_taylor_real: Hamiltonian norm estimate = %g",scaling);
	if (scaling > 1) {
		nsq = (int)ceil(log(scaling)/log(2.0));
		scaling = pow(2,nsq);
	} else {
		nsq = 0;
		scaling = 1;
	}
	//printf(", will do %d squarings\n",nsq);
	dm_muld(ham,-dt/scaling);
	mat_complx *prop_dum = dm_imag2(ham); // result is always MAT_DENSE
	for (i=0; i<len; i+=dim+1) {
		prop_dum->data[i].re += 1;
	}
	// now prop = 1 - i H dt / sc
	next_term = dm_dup2(ham);
	//cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,1.0,ham->data,dim,ham->data,dim,0.0,next_term->data,dim);
	dm_multo_rev(next_term,ham);
	fac = 0.5;
	//next_term_norm = dm_normest(next_term)*fac;
	next_term_norm = fabs(dm_max(next_term)*fac);
	int n = 2;
	while (next_term_norm > norm_tol) {
		assert(next_term->type == MAT_DENSE);
		switch (n%4) {
		case 0:
			//cblas_daxpy(len,fac,next_term->data,1,(double*)(prop->data),2);
			cblas_daxpy(len,fac,next_term->data,1,(double*)(prop_dum->data),2);
			break;
		case 1:
			//cblas_daxpy(len,fac,next_term->data,1,(double*)(prop->data)+1,2);
			cblas_daxpy(len,fac,next_term->data,1,(double*)(prop_dum->data)+1,2);
			break;
		case 2:
			//cblas_daxpy(len,-fac,next_term->data,1,(double*)(prop->data),2);
			cblas_daxpy(len,-fac,next_term->data,1,(double*)(prop_dum->data),2);
			break;
		case 3:
			//cblas_daxpy(len,-fac,next_term->data,1,(double*)(prop->data)+1,2);
			cblas_daxpy(len,-fac,next_term->data,1,(double*)(prop_dum->data)+1,2);
			break;
		}
		fac /= (double)(++n);
		dm_multo_rev(next_term,ham);
		//next_term_norm = dm_normest(next_term)*fac;
		next_term_norm = fabs(dm_max(next_term)*fac);
	}
	free_double_matrix(next_term);
	//printf("prop_real dense taylor: %d Taylor steps\n",n-1);
	TIMING_TOC(tv1,tv2,"Taylor after iterations");

	if (nsq) {
		for (n=0; n<nsq; n++) {
			//cm_multo(prop,prop);
			cm_multo(prop_dum,prop_dum);
			DEBUGPRINT("prop_real dense taylor: SQUARING %d done\n",n);
		}
	}
	TIMING_TOC(tv1,tv2,"Taylor after squarings");

	cm_multo_rev(prop,prop_dum);
	free_complx_matrix(prop_dum);
	TIMING_TOC(tv1,tv2,"Taylor after prop update");
}

void prop_diag1_real(mat_complx *prop, mat_double *ham, double dt)
{
	int i, j, info = 0, ldwsp;
	double *dwsp, *eigs, *dT = ham->data;
	int dim = ham->row;

	TIMING_INIT_VAR2(tv1,tv2);
	TIMING_TIC(tv1);

	dm_muld(ham, dt);
	dwsp = (double*)malloc((2*dim*dim)*sizeof(double));
	ldwsp = 2*dim*dim-dim;
	eigs = dwsp+ldwsp;

	dsyev_("V","U",&dim,dT,&dim, eigs, dwsp, &ldwsp, &info);
	if (info != 0) {
		fprintf(stderr,"Error: prop_diag1_real - hamiltonian diagonalization failed (info = %d)\n",info);
		exit(1);
	}
	mat_complx *cm = complx_matrix(dim,dim,MAT_DENSE,0,ham->basis);
	complx *dV;
	double cs, ss;
	for (j=0; j<dim; j++) {
		dV = cm->data + j;
		cs = cos(eigs[j]);
		ss = -sin(eigs[j]);
		for (i=0; i<dim; i++) {
			dV->re = (*dT)*cs;
			dV->im = (*dT)*ss;
			dT++;
			dV += dim;
		}
	}
	/* result of this(^) should be cm = exp(-i*eigs)*Tansf.Mat.+ */
	cm_multo(cm,prop); // in order to update existing prop
	if (prop->type != MAT_DENSE) {
		mat_complx *dum = complx_matrix(dim,dim,MAT_DENSE, dim, ham->basis);
		cm_swap_innards_and_destroy(prop,dum);
	}
	zlarcm_(&dim,&dim,ham->data,&dim,cm->data,&dim,prop->data,&dim,dwsp);
	free_complx_matrix(cm);
	free(dwsp);
	TIMING_TOC(tv1,tv2,"DSYEV done");
}

void prop_diag2_real(mat_complx *prop, mat_double *ham, double dt)
{
	int i, j, info = 0;
	int dim = ham->row;

	TIMING_INIT_VAR2(tv1,tv2);
	TIMING_TIC(tv1);

	dm_muld(ham, dt);
	double *eigs = (double*)malloc(dim*sizeof(double));

	i = 1;
	double abstol = -1.0;
	mat_double *T = double_matrix(dim,dim,MAT_DENSE,0,ham->basis);
	int *isupp = (int *)malloc(2*dim*sizeof(int));
	double wkopt;
	int lwork = -1;
	int iwkopt;
	int liwork = -1;
	dsyevr_("V","I","U",&dim,ham->data,&dim,NULL,NULL,&i,&dim,&abstol,&j,
			eigs,T->data,&dim,isupp,&wkopt,&lwork,&iwkopt,&liwork, &info);
	// j is total number of eigval found
	lwork = (int)wkopt;
	double *work = (double*)malloc( lwork*sizeof(double) );
	liwork = iwkopt;
	int *iwork = (int*)malloc( liwork*sizeof(int) );
	// using wkopt as dummy variable to prevent segfault on some systems
	dsyevr_("V","I","U",&dim,ham->data,&dim,&wkopt,&wkopt,&i,&dim,&abstol,&j,
			eigs,T->data,&dim,isupp,work,&lwork,iwork,&liwork, &info);
	// j is total number of eigval found
	if (info != 0) {
		fprintf(stderr,"Error: prop_diag2_real - hamiltonian diagonalization failed (info = %d)\n",info);
		exit(1);
	}
	//dm_print(T,"Transf. matrix");
	free(work);
	free(iwork);
	free(isupp);
	mat_complx *cm = complx_matrix(dim,dim,MAT_DENSE,0,ham->basis);
	complx *dV;
	double cs, ss, *dT = T->data;
	for (j=0; j<dim; j++) {
		dV = cm->data + j;
		cs = cos(eigs[j]);
		ss = -sin(eigs[j]);
		for (i=0; i<dim; i++) {
			dV->re = (*dT)*cs;
			dV->im = (*dT)*ss;
			dT++;
			dV += dim;
		}
	}
	/* result of this(^) should be cm = exp(-i*eigs)*Tansf.Mat.+ */
	free(eigs);
	cm_multo(cm,prop); // in order to update existing prop
	if (prop->type != MAT_DENSE) {
		mat_complx *dum = complx_matrix(dim,dim,MAT_DENSE, dim, ham->basis);
		cm_swap_innards_and_destroy(prop,dum);
	}
	work = (double *)malloc(2*dim*dim*sizeof(double));
	zlarcm_(&dim,&dim,T->data,&dim,cm->data,&dim,prop->data,&dim,work);
	free_complx_matrix(cm);
	free_double_matrix(T);
	free(work);
	TIMING_TOC(tv1,tv2,"DSYEVR done");
}

/****
 *  Note: calculates prop = exp(-i*H*dt)*prop
 *        matrix 'ham' gets overwritten with intermediate results!!!
 ****/
void prop_real(mat_complx *prop, mat_double *ham, double dt, int method)
{

	int dim = ham->row;

	assert(dim == ham->col && dim == prop->row && dim == prop->col);
	assert(prop->basis == ham->basis);
	switch (ham->type) {
		case MAT_DENSE_DIAG : {
			assert(prop->type == MAT_DENSE_DIAG || prop->type == MAT_SPARSE_DIAG);
			if (prop->type == MAT_SPARSE_DIAG) cm_dense(prop);
			int i;
			for (i=0; i<dim; i++) {
				double re = cos(ham->data[i]*dt);
				double im = -sin(ham->data[i]*dt);
				double dum = prop->data[i].re;
				prop->data[i].re = dum*re - prop->data[i].im*im;
				prop->data[i].im = dum*im + prop->data[i].im*re;
			}
			break; }
		case MAT_SPARSE_DIAG : {
			assert(prop->type == MAT_DENSE_DIAG || prop->type == MAT_SPARSE_DIAG);
			if (prop->type == MAT_SPARSE_DIAG) cm_dense(prop);
			int i;
			for (i=0; i<dim; i++) {
				double dum = dm_getelem(ham,i+1,i+1);
				double re = cos(dum*dt);
				double im = -sin(dum*dt);
				dum = prop->data[i].re;
				prop->data[i].re = dum*re - prop->data[i].im*im;
				prop->data[i].im = dum*im + prop->data[i].im*re;
			}
			break; }
		case MAT_DENSE : {
			switch (method) {
			case 0: {// via diagonalization using dsyevr (might not be thread safe with some LAPACK libs)
				prop_diag2_real(prop,ham,dt);
				break; }
			case 1: {// via scaling & squaring with Pade approx
				//printf("\nHURA, Padde!\n\n");
				prop_pade_real(prop,ham,dt);
				break; }
			case 2: {// via Chebyshev with Scaling & Squaring
#ifdef INTEL_MKL
				dm_sparse(ham,SPARSE_TOL);
#endif
				prop_cheby2_real(prop,ham,dt);
				break; }
			case 3: { // via Taylor
#ifdef INTEL_MKL
				dm_sparse(ham,SPARSE_TOL);
#endif
				prop_taylor_real(prop, ham, dt);
				break; }
			case 4: // via Lanczos
				fprintf(stderr,"prop_real error: lanczos method not implemented yet\n");
				exit(-1);
			case 5: {// via diagonalization using dsyev (found thread safe but slower)
				prop_diag1_real(prop,ham,dt);
				break; }
			case 6: {// via Chebyshev with Shifting & Scaling
#ifdef INTEL_MKL
				dm_sparse(ham,SPARSE_TOL);
#endif
				prop_cheby3_real(prop,ham,dt);
				break; }
			default:
				fprintf(stderr,"prop_real error: unknown calculation method (%d), use diagonalization, pade, chebyshev, taylor or lanczos\n",method);
				exit(-1);
			}
			break; }
		case MAT_SPARSE : {
			switch (method) {
			case 0:
				fprintf(stderr,"prop_real error: diagonalization method not allowed for sparse matrices, use Chebyshev or Taylor\n");
				exit(-1);
			case 1:
				fprintf(stderr,"prop_real error: Pade method not allowed for sparse matrices, use Chebyshev or Taylor\n");
				exit(-1);
			case 6:  // Chebyshev
				prop_cheby4_real(prop,ham,dt);
				break;
			case 3: { // Taylor
				const double norm_tol = 1.0e-6;
				mat_double *next_term, *ure, *uim;
				double next_term_norm, fac, scaling;
				int nsq;
				/* scaling step */
				scaling = dm_normest(ham)*dt;
				DEBUGPRINT("prop_real: sparse variant, Hamiltonian largest eigval estimate = %g\n",scaling);
				if (scaling > 1) {
					nsq = (int)ceil(log(scaling)/log(2.0));
					scaling = pow(2,nsq);
				} else {
					nsq = 0;
					scaling = 1;
				}
				//printf("prop_real: sparse variant, will do %d squarings\n",nsq);
				dm_muld(ham,-dt/scaling);
				ure = dm_creatediag('s',dim,1,ham->basis);	// U_re = 1
				uim = dm_dup(ham);							// U_im = - H dt / sc
				next_term = simpson_dcsrmultcsr(ham,ham);
				fac = 0.5;
				dm_sparse(next_term,SPARSE_TOL*fac);
				next_term_norm = dm_normest(next_term)*fac;
				int n = 2;
				while (next_term_norm > norm_tol) {
					//printf("prop_real: ure is %s (%d), uim is %s (%d)\n",matrix_type(ure->type),dm_nnz(ure),matrix_type(uim->type),dm_nnz(uim));
					//printf("         : next_term is %s (%d)\n",matrix_type(next_term->type),dm_nnz(next_term));
					switch (n%4) {
					case 0:
						// add to real part
						dm_multod(ure,next_term,fac);
						break;
					case 1:
						// add to imaginary part
						dm_multod(uim,next_term,fac);
						break;
					case 2:
						// subtract from real part
						dm_multod(ure,next_term,-fac);
						break;
					case 3:
						// subtract from imaginary part
						dm_multod(uim,next_term,-fac);
						break;
					}
					fac /= (double)(++n);
					dm_multo_rev(next_term,ham); // cleaning done, ham is sparse, next_term can be dense
					next_term_norm = dm_normest(next_term)*fac;
					DEBUGPRINT("next term norm = %f\n",next_term_norm);
				}
				free_double_matrix(next_term);
				//printf("prop_real sparse: %d Taylor steps\n",n-1);
				// assemble complex propagator
				//if (prop->type != MAT_SPARSE) {
				//	mat_complx *dum = dm_complx(ure);
				//	cm_swap_innards_and_destroy(prop,dum);
				//} else {
				//	dm_copy2cm(ure,prop);
				//}
				mat_complx *prop_dum = dm_complx(ure);
				mat_complx *dumim = dm_imag(uim);
				//cm_addto(prop,dumim);
				cm_addto(prop_dum,dumim);
				free_complx_matrix(dumim);
				free_double_matrix(ure);
				free_double_matrix(uim);
				if (nsq) {
					for (n=0; n<nsq; n++) {
						//cm_multo(prop,prop);
						cm_multo(prop_dum,prop_dum);
						DEBUGPRINT("prop_real sparse: SQUARING %d done\n",n);
					}
				}
				cm_multo_rev(prop,prop_dum);
				cm_sparse(prop,SPARSE_TOL);
				free_complx_matrix(prop_dum);
				break; }
			default:
				fprintf(stderr,"prop_real error: unknown calculation method %d. Use cheby2 or taylor\n",method);
				exit(-1);
			}
			break; }
		default :
			fprintf(stderr,"Error: prop_real - unknown hamiltonian type '%d'\n",ham->type);
			exit(1);
		}
	//cm_print(prop,"PROPAGATOR");
}

/****
 * Note: this does NOT updates prop, it REPLACES prop with exp(-i*H*dt)
 ****/
void prop_complx(mat_complx *prop, mat_complx *ham, double dt, int method)
{
	/* assumes hermitian Hamiltonian, but only for diagonal types... */
	int dim = ham->row;

	assert(dim == ham->col && dim == prop->row && dim == prop->col);
	switch (ham->type) {
		case MAT_DENSE_DIAG : {
			if (prop->type != MAT_DENSE_DIAG) {
				mat_complx *dum = complx_matrix(dim,dim,MAT_DENSE_DIAG, dim, ham->basis);
				cm_swap_innards_and_destroy(prop,dum);
			}
			prop->basis = ham->basis;
			complx *dres = prop->data, *ddv = ham->data, *stop = ddv + dim;
			double arg;
			do {
				arg = (ddv->re)*dt;
				dres->re = cos(arg);
				dres->im = -sin(arg);
				ddv++;
				dres++;
			} while(ddv != stop);
			break; }
		case MAT_SPARSE_DIAG : {
			if (prop->type != MAT_SPARSE_DIAG) {
				mat_complx *dum = complx_matrix(dim,dim,MAT_SPARSE_DIAG, dim, ham->basis);
				cm_swap_innards_and_destroy(prop,dum);
			} else {
				cm_change_nnz(prop, dim);
			}
			prop->basis = ham->basis;
			int i;
			complx z;
			for (i=0; i<dim; i++) {
				z = cm_getelem(ham,i+1,i+1);
				prop->data[i] = Cexpi(-dt*z.re);
				prop->icol[i] = prop->irow[i] = i+1;
			}
			prop->irow[dim] = dim + 1;
			break; }
		case MAT_DENSE : {
			/* this is just diag */
			complx *eigs=complx_vector(dim);
			mat_complx *T=complx_matrix(dim,dim,MAT_DENSE,0,ham->basis);
			int i;
			cm_muld(ham,dt);
		//cm_print(ham,"------> H*0.01");
			cm_diag(ham,eigs,T);
		//for (i=1; i<=dim; i++) printf("   eig[%i] = (%g,%g)\n",i,eigs[i].re,eigs[i].im);
		//cm_print(T,"====> T");
			if (prop->type != MAT_DENSE_DIAG) {
				mat_complx *dum = complx_matrix(dim,dim,MAT_DENSE_DIAG, dim, ham->basis);
				cm_swap_innards_and_destroy(prop,dum);
			}
			prop->basis = ham->basis;
			for (i=0; i<dim; i++) {
				prop->data[i] = Cexpi(-eigs[i+1].re);
			}
		//cm_print(prop,"+++++>prop");
			simtrans(prop,T);
			free_complx_vector(eigs);
			free_complx_matrix(T);
		//cm_print(prop,"~~~~~> prop");
			break; }
		case MAT_SPARSE : {
			/* this is just Taylor */
			const double norm_tol = 1.0e-3; // machine precision is reached with 1e-4
			mat_complx *next_term;
			double next_term_norm, fac, scaling;
			int nsq;
			/* scaling step */
			scaling = cm_normest(ham)*dt;
			DEBUGPRINT("prop_complx: sparse variant, Hamiltonian largest eigval estimate = %g\n",scaling);
			if (scaling > 1) {
				nsq = (int)ceil(log(scaling)/log(2.0));
				scaling = pow(2,nsq);
			} else {
				nsq = 0;
				scaling = 1;
			}
			DEBUGPRINT("prop_complx: sparse variant, will do %d squarings\n",nsq);
			cm_mulc(ham,Complx(0.0,-dt/scaling));
			if (prop->type != MAT_SPARSE) {
				mat_complx *dum = complx_matrix(dim,dim,MAT_SPARSE, dim, ham->basis);
				cm_swap_innards_and_destroy(prop,dum);
			}
			prop->basis = ham->basis;
			cm_unit(prop);
			next_term = complx_matrix(dim, dim, MAT_SPARSE, dim, ham->basis);
			cm_unit(next_term);
			fac = 1.0;
			int n = 1;
			do {
				cm_multo(next_term, ham);
				cm_multod(prop,next_term,fac);
				// cleaning is already done in cm_multo
				//cm_sparse(next_term,SPARSE_TOL);
				next_term_norm = cblas_dzasum(next_term->irow[dim]-1,next_term->data,1)*fac;
				fac /= (double)(++n);
			} while (next_term_norm > norm_tol);
			free_complx_matrix(next_term);
			cm_sparse(prop,SPARSE_TOL);
			DEBUGPRINT("prop_real sparse: %d Taylor steps\n",n-1);
			if (nsq) {
				for (n=0; n<nsq; n++) {
					cm_multo(prop,prop);
					DEBUGPRINT("prop_complx sparse: SQUARING %d done\n",n);
				}
			}
			break; }
		default :
			fprintf(stderr,"Error: prop_complx - unknown hamiltonian type '%d'\n",ham->type);
			exit(1);
	}

}


mat_complx * cm_ln(mat_complx *m)
{
	int i, dim = m->row;
	mat_complx *res;

	if (dim != m->col) {
		fprintf(stderr,"Error: cm_ln - non-square matrix (%d,%d)\n",m->row,m->col);
		exit(1);
	}

	switch (m->type) {
		case MAT_DENSE_DIAG :
			res = complx_matrix(dim, dim, MAT_DENSE_DIAG,0,m->basis);
			for (i=0; i<dim; i++) {
				res->data[i] = Clog(m->data[i]);
			}
			break;
		case MAT_SPARSE_DIAG : {
			complx z;
			res = complx_matrix(dim, dim, MAT_SPARSE_DIAG,dim,m->basis);
			for (i=0; i<dim; i++) {
				z = cm_getelem(m,i+1,i+1);
				res->data[i] = Clog(z);
				res->irow[i] = res->icol[i] = i+1;
			}
			res->irow[dim] = dim+1;
			break; }
		case MAT_DENSE :
		case MAT_SPARSE : {
			complx *eigs=complx_vector(dim);
			mat_complx *T=complx_matrix(dim,dim,MAT_DENSE,0,m->basis);
			cm_diag(m,eigs,T);
			res = complx_matrix(dim,dim,MAT_DENSE_DIAG,0,m->basis);
			for (i=0; i<dim; i++) {
				res->data[i] = Clog(eigs[i+1]);
			}
			simtrans(res,T);
			free_complx_vector(eigs);
			free_complx_matrix(T);
			break; }
	   default :
		   fprintf(stderr,"Error: cm_ln - unknown type '%d'\n",m->type);
		   exit(1);
	}
	return res;
}

/* matrix exponential using scaling and squaring with Pade approximation */
mat_complx * cm_expm_pade(mat_complx *m)
{
	//double c[] ={1.0, 0.5, 5.0/44.0, 1.0/66.0, 1.0/792.0, 1.0/15840.0, 1.0/665280.0};
	double c, scaling;
	int nsq, dim, p, q, k, info=1;
	mat_complx *res, *mX, *mD, *mA;

	if (m->type != MAT_DENSE) {
		fprintf(stderr,"Error: cm_expm_pade can be used only with dense matrices\n");
		exit(1);
	}
	dim = m->row;
	if (dim != m->col) {
		fprintf(stderr,"Error: cm_expm_pade - matrix is not square (%d,%d)\n",m->row,m->col);
		exit(1);
	}

	/* scaling step */
	scaling = cm_normest(m);
	DEBUGPRINT("cm_expm_pade: matrix largest eigval estimate = %g\n",scaling);
	if (scaling > 1) {
		nsq = (int)ceil(log(scaling)/log(2.0))+1;
		scaling = pow(2,-nsq);
	} else {
		nsq = 0;
		scaling = 1;
	}
	DEBUGPRINT("cm_expm_pade: will do %d squarings\n",nsq);
	mA = cm_dup(m);
	cm_muld(mA,scaling);

	/* algorithm from matlab's expm1 */
	c = 0.5;
	mX = cm_dup(mA);
	res = cm_creatediag('D', dim, Complx(1.0,0.0),m->basis);
	cm_multod(res,mA,c);
	mD = cm_creatediag('D', dim, Complx(1.0,0.0),m->basis);
	cm_multod(mD,mA,-c);
	p = 1; q = 6;
	for (k=2; k<=q; k++) {
		c *= (double)(q-k+1) / (double)(k*(2*q-k+1));
		cm_multo_rev(mX,mA);
		cm_multod(res,mX,c);
		if (p) {
			cm_multod(mD,mX,c);
		} else {
			cm_multod(mD,mX,-c);
		}
		p = !p;
	}
	int *pvec = (int*)malloc((dim+1)*sizeof(int));
	zgesv_(&dim,&dim,mD->data,&dim,pvec,res->data,&dim,&info);
	if (info != 0) {
		fprintf(stderr,"cm_expm_pade error: zgesv failed with info=%d\n",info);
		exit(1);
	}
    free(pvec);

	/* squaring */
	if (nsq) {
		for (k=0; k<nsq; k++) {
			cm_multo(res,res);
		}
	}
	free_complx_matrix(mA);
	free_complx_matrix(mD);
	free_complx_matrix(mX);
	return res;
}

/****
 * this uses expv from expokit rewritten from matlab code
 * should make use of Krylov subspaces technique
 * calculates vec = exp(-i*t*A)*vec
 * AND          A = -i*A !!!
 */
void prop_krylov(double t, mat_complx *A, complx *vec)
{
	int n,m,mb,nstep,cc,rr,dimH,i,j,k1,mx,ireject,mxrej;
	double tol,anorm,avnorm,normv,btol,beta,gamma,delta,t_out,err_loc,fact,phi1,phi2;
	double rndoff, s_error,s,sgn,t_new,t_now,t_step,xm;
	complx *w, **Vecs, *H, *dumvec,*p,zval;
	mat_complx *dummat,*F;

	n = A->row;
	if (n != A->col) {
		fprintf(stderr,"Error: prop_krylov - got non-square matrix\n");
		exit(1);
	}
	if (n != LEN(vec)) {
		fprintf(stderr,"Error: prop_krylov - dimension mismatch A=(%d,%d), vec=(%d)\n", A->row, A->col, LEN(vec));
		exit(1);
	}
	if (fabs(t) < TINY) {
		fprintf(stderr,"Warning: prop_krylov got tiny time-step, does NOTHING\n");
		return;
	}

	cm_mulc(A,Complx(0.0,-1.0));
	tol = 1e-7;
	if (n<30) m = n; else m = 30;
	anorm = cm_normest(A);
	mxrej = 10;
	btol = 1e-7;
	gamma = 0.9;
	delta = 1.2;
	mb = m;
	t_out = fabs(t);
	nstep = 0;
	//t_new = 0;
	t_now = 0;
	s_error = 0;
	rndoff = anorm*1e-15;
	k1 = 2;
	xm = 1.0/(double)(m);
	normv = cv_norm(vec);
	beta = normv;
	fact = pow((m+1)/2.71828182845904523536,m+1)*sqrt(2*M_PI*(m+1));
	t_new = 1.0/anorm*pow((fact*tol)/(4*beta*anorm),xm);
	s = pow(10,floor(log(t_new)/log(10.0))-1);
	t_new = ceil(t_new/s)*s;
	if (t>0) sgn = 1; else sgn = -1;

	w = cv_dup(vec);
	Vecs = (complx**)malloc((m+1)*sizeof(complx*));
	for (i=0; i<=m; i++) {
		Vecs[i] = complx_vector(n);
	}
	dimH = m+2;
	H = (complx*)malloc(dimH*dimH*sizeof(complx));
	p = complx_vector(n); cv_zero(p);
	dumvec = complx_vector(n);
	while (t_now < t_out) {
		for (i=0; i<=m; i++) {
			cv_zero(Vecs[i]);
		}
		memset(H,0,dimH*dimH*sizeof(complx));
		nstep++;
		t_step = (t_out-t_now < t_new) ? t_out-t_now : t_new;
		cv_copy(Vecs[0],w); cv_muld(Vecs[0],1.0/beta);
		for (j=1; j<=m; j++) {
			cv_matmul(p,A,Vecs[j-1]);
			for (i=1; i<=j; i++) {
				H[i-1+dimH*(j-1)] = zval = cv_dotc(Vecs[i-1],p);
				zval.re = -zval.re; zval.im = -zval.im;
				cv_multoc(p,Vecs[i-1],zval);
			}
			s = cv_norm(p);
			if (s < btol) {
				k1 = 0;
				mb = j;
				t_step = t_out - t_now;
				break;
			}
			H[j+dimH*(j-1)] = Complx(s,0.0);
			cv_copy(Vecs[j],p); cv_muld(Vecs[j],1.0/s);
		}
		if (k1 != 0) {
			H[m+1 + dimH*m] = Complx(1,0);
			cv_matmul(dumvec,A,Vecs[m]);
			avnorm = cv_norm(dumvec);
		}
		ireject = 0;
		while (ireject <= mxrej) {
			mx = mb + k1;
			if (mx == 1) {
				F = complx_matrix(1,1,MAT_DENSE,0,A->basis);
				F->data[0] = Cexp( CRmul(H[0],sgn*t_step) );
			} else {
				dummat = complx_matrix(mx,mx,MAT_DENSE,0,A->basis);
				for (rr=0; rr<mx; rr++) {
					for (cc=0; cc<mx; cc++) {
						dummat->data[rr+mx*cc] = CRmul(H[rr+dimH*cc],sgn*t_step);
					}
				}
				F = cm_expm_pade(dummat);
				free_complx_matrix(dummat);
			}
			if (k1 == 0) {
				err_loc = btol;
				break;
			} else {
				zval = F->data[m];
				phi1 = fabs( beta*sqrt(zval.re*zval.re+zval.im*zval.im));
				zval = F->data[m+1];
				phi2 = fabs( beta*sqrt(zval.re*zval.re+zval.im*zval.im)*avnorm);
				if (phi1 > 10*phi2) {
					err_loc = phi2;
					xm = 1.0/(double)m;
				} else if (phi1>phi2) {
					err_loc = phi1*phi2/(phi1-phi2);
					xm = 1.0/(double)m;
				} else {
					err_loc = phi1;
					xm = 1.0/(double)(m-1);
				}
			}
			if (err_loc <= delta*t_step*tol) {
				break;
			} else {
				t_step = gamma*t_step*pow(t_step*tol/err_loc,xm);
				s = pow(10,floor(log(t_step)/log(10.0))-1);
				t_step = ceil(t_step/s)*s;
				if (ireject == mxrej) {
					fprintf(stderr,"Error: prop_krylov - tolerance is too high\n");
					exit(1);
				}
				ireject++;
			}
		}
		mx = mb + ( (k1-1 > 0) ? k1-1 : 0 );
		for (cc=1; cc<=n; cc++) {
			w[cc] = Cnull;
			for (rr=0; rr<mx; rr++) {
				zval = Cmul(Vecs[rr][cc],F->data[rr]);
				w[cc].re += beta*zval.re;
				w[cc].im += beta*zval.im;
			}
		}
		free_complx_matrix(F);
		beta = cv_norm(w);

		t_now += t_step;
		t_new = gamma*t_step*pow(t_step*tol/err_loc,xm);
		s = pow(10,floor(log(t_new)/log(10.0))-1);
		t_new = ceil(t_new/s)*s;

		err_loc = (err_loc > rndoff) ? err_loc : rndoff;
	}

	/* clean up */
	free(H);
	for (i=0; i<=m; i++) free_complx_vector(Vecs[i]);
	free(Vecs);
	free_complx_vector(p);
	free_complx_vector(dumvec);

	cv_copy(vec,w);
	free_complx_vector(w);
	DEBUGPRINT("prop_krylov: ---> DONE <---\n\n\n");

}

mat_complx * cm_power(mat_complx *mx, int n)
{
	unsigned int cnt;
	mat_complx *res=NULL, *dum=NULL;

	assert(mx->row == mx->col);
	if ( n <= 0 ) {
		fprintf(stderr,"cm_power error: matrix power has to be nonzero positive (not %d)\n",n);
		exit(1);
	}

	cnt = (unsigned int)n;
	dum = cm_dup(mx);
	while (1) {
		if (cnt & (1<<0)) {
			if (res != NULL) {
				cm_multo(res,dum);
			} else {
				res = cm_dup(dum);
			}
		}
		cnt >>= 1;
		if (cnt == 0) break;
		cm_multo(dum,dum);
	}
	free_complx_matrix(dum);
	assert(res != NULL);
	return res;
}


/****
 * transpose double matrix in-place
 ****/
void dm_transposei(mat_double *m)
{
	switch (m->type) {
		case MAT_DENSE_DIAG :
		case MAT_SPARSE_DIAG :
			break;
		case MAT_DENSE : {
			int i;
			mat_double *res = double_matrix(m->col, m->row, MAT_DENSE,0,m->basis);
			for (i=0; i<m->row; i++) {
				cblas_dcopy(m->col, &(m->data[i]),m->row,&res->data[i*res->row],1);
			}
			dm_swap_innards_and_destroy(m,res);
			break; }
		case MAT_SPARSE : {
			int len = m->irow[m->row] - 1;
			MKL_INT *dumrow, *dumcol, *ic = m->icol;
			dumrow = (MKL_INT*)malloc((m->col+1+len)*sizeof(MKL_INT));
			dumcol = dumrow + m->col+1;
			int i, j, nn=0;
			double *z;
			mat_double *res = double_matrix(m->col, m->row, MAT_SPARSE, len, m->basis);
			memset(dumrow,0,(m->col+1)*sizeof(MKL_INT));
			for (i=0; i<m->row; i++) {
				int n = m->irow[i+1] - m->irow[i];
				for (j=0; j<n; j++) {
					dumrow[*ic]++;
					dumcol[nn] = i+1;
					ic++;
					nn++;
				}
			}
			res->irow[0] = dumrow[0] = 1;
			for (i=1; i<=res->row; i++) {
				res->irow[i] = res->irow[i-1] + dumrow[i];
				//DEBUGPRINT("   dm_transposei res->irow[%d]=%d\n",i,res->irow[i]);
				dumrow[i] = res->irow[i];
			}
			z = m->data;
			for (i=0; i<len; i++) {
				j = dumrow[m->icol[i]-1]-1;
				res->data[j] = *z;
				res->icol[j] = dumcol[i];
				z++;
				dumrow[m->icol[i]-1]++;
				//DEBUGPRINT("   dm_transposei res->icol[%d] = %d\n",j,res->icol[j]);
			}
			dm_swap_innards_and_destroy(m,res);
			break; }
	   default :
		   fprintf(stderr,"Error: dm_transposei - unknown type '%d'\n",m->type);
		   exit(1);
	}

}


void cm_permute(mat_complx *m, int *permvec)
{
	assert(m->row == m->col);
	assert(m->row == LEN(permvec));
	int dim = m->row;

	switch (m->type) {
	case MAT_DENSE: {
		int i, j;
		complx *zv = (complx*)malloc(dim*dim*sizeof(complx));
		for (i=0; i<dim; i++) {
			for (j=0; j<dim; j++) {
				zv[i+j*dim] = m->data[(permvec[i+1]-1) + (permvec[j+1]-1)*dim];
			}
		}
		free(m->data);
		m->data = zv;
		break; }
	case MAT_DENSE_DIAG: {
		complx *zv = (complx*)malloc(dim*sizeof(complx));
		int i;
		for (i=0; i<dim; i++) zv[i] = m->data[permvec[i+1]-1];
		free(m->data);
		m->data = zv;
		break; }
	case MAT_SPARSE: {
		//cm_print(m,"cm_permute stara matice");
		int i, j, jj, c_orig, k, r, nc, nnz = m->irow[dim]-1;
		mat_complx *mp = complx_matrix(dim,dim,MAT_SPARSE,nnz,m->basis);
		mp->irow[0] = 1;
		int pos = 0;
		for (i=1; i<=dim; i++) {
			r = permvec[i];
			nc = m->irow[r] - m->irow[r-1];
			//printf("novy radek %i je stary %i, ktery ma %i clenu\n",i,r,nc);
			for (j=0; j<nc; j++) {
				c_orig = m->icol[m->irow[r-1]-1+j];
				for (k=1; k<=dim; k++) {
					if (permvec[k] == c_orig) break;
				}
				//printf("\tstary sl %i patri na novy sl %i\n",c_orig,k);
				if (j == 0) {
					mp->icol[pos] = k;
					mp->data[pos] = m->data[m->irow[r-1]-1+j];
				} else {
					for (jj=0; jj<j; jj++) {
						if (k < mp->icol[pos-jj-1]) {
							mp->icol[pos-jj] = mp->icol[pos-jj-1];
							mp->data[pos-jj] = mp->data[pos-jj-1];
						} else {
							break;
						}
					}
					mp->icol[pos-jj] = k;
					mp->data[pos-jj] = m->data[m->irow[r-1]-1+j];
				}
				pos++;
			}
			mp->irow[i] = pos+1;
		}
		cm_swap_innards_and_destroy(m,mp);
		//cm_print(m,"cm_permute nova matice");
		//exit(1);
		break; }
	case MAT_SPARSE_DIAG: {
		int i,j, nnz;
		nnz = m->irow[dim]-1;
		mat_complx *mp = complx_matrix(dim,dim,MAT_SPARSE_DIAG,nnz,m->basis);
		mp->irow[0] = 1;
		j = 0;
		complx zval;
		for (i=1; i<=dim; i++) {
			mp->irow[i] = mp->irow[i-1];
			zval = cm_getelem(m,permvec[i],permvec[i]);
			if (fabs(zval.re)<SPARSE_TOL && fabs(zval.im)<SPARSE_TOL) continue;
			mp->icol[j] = i;
			mp->data[j] = zval;
			(mp->irow[i])++;
			j++;
		}
		cm_swap_innards_and_destroy(m,mp);
		break; }
	default :
		fprintf(stderr,"Error: cm_permute - unknown matrix type %d\n",m->type);
		exit(1);
	}
}

void dm_permute(mat_double *m, int *permvec)
{
	assert(m->row == m->col);
	assert(m->row == LEN(permvec));
	int dim = m->row;

	switch (m->type) {
	case MAT_DENSE: {
		int i, j;
		double *zv = (double*)malloc(dim*dim*sizeof(double));
		for (i=0; i<dim; i++) {
			for (j=0; j<dim; j++) {
				zv[i+j*dim] = m->data[(permvec[i+1]-1) + (permvec[j+1]-1)*dim];
			}
		}
		free(m->data);
		m->data = zv;
		break; }
	case MAT_DENSE_DIAG: {
		double *zv = (double*)malloc(dim*sizeof(double));
		int i;
		for (i=0; i<dim; i++) zv[i] = m->data[permvec[i+1]-1];
		free(m->data);
		m->data = zv;
		break; }
	case MAT_SPARSE: {
		int i, j, jj, c_orig, k, r, nc, nnz = m->irow[dim]-1;
		mat_double *mp = double_matrix(dim,dim,MAT_SPARSE,nnz,m->basis);
		mp->irow[0] = 1;
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
					mp->icol[pos] = k;
					mp->data[pos] = m->data[m->irow[r-1]-1+j];
				} else {
					for (jj=0; jj<j; jj++) {
						if (k < mp->icol[pos-jj-1]) {
							mp->icol[pos-jj] = mp->icol[pos-jj-1];
							mp->data[pos-jj] = mp->data[pos-jj-1];
						} else {
							break;
						}
					}
					mp->icol[pos-jj] = k;
					mp->data[pos-jj] = m->data[m->irow[r-1]-1+j];
				}
				pos++;
			}
			mp->irow[i] = pos+1;
		}
		dm_swap_innards_and_destroy(m,mp);
		break; }
	case MAT_SPARSE_DIAG: {
		int i,j, nnz;
		nnz = m->irow[dim]-1;
		mat_double *mp = double_matrix(dim,dim,MAT_SPARSE_DIAG,nnz,m->basis);
		mp->irow[0] = 1;
		j = 0;
		double zval;
		for (i=1; i<=dim; i++) {
			mp->irow[i] = mp->irow[i-1];
			zval = dm_getelem(m,permvec[i],permvec[i]);
			if (fabs(zval)<SPARSE_TOL ) continue;
			mp->icol[j] = i;
			mp->data[j] = zval;
			(mp->irow[i])++;
			j++;
		}
		dm_swap_innards_and_destroy(m,mp);
		break; }
	default :
		fprintf(stderr,"Error: dm_permute - unknown matrix type %d\n",m->type);
		exit(1);
	}
}

complx cm_true_trace(mat_complx *m)
{
	complx *z, res = Cnull;
	int i;

	assert(m->row == m->col);

	switch (m->type) {
		case MAT_DENSE :
			z = m->data;
			for (i=0; i<m->row; i++) {
				res.re += z->re;
				res.im += z->im;
				z += m->row+1;
			}
			break;
		case MAT_DENSE_DIAG :
			z = m->data;
			for (i=0; i<m->row; i++) {
				res.re += z->re;
				res.im += z->im;
				z++;
			}
			break;
		case MAT_SPARSE :
			for (i=1; i<=m->row; i++) {
				res = Cadd(res,cm_getelem(m,i,i));
			}
			break;
		case MAT_SPARSE_DIAG : {
			int len = m->irow[m->row] - 1;
			z = m->data;
			for (i=0; i<len; i++) {
				res.re += z->re;
				res.im += z->im;
				z++;
			}
			break; }
	   default :
		   fprintf(stderr,"Error: cm_true_trace - unknown type '%d'\n",m->type);
		   exit(1);
	}

	return res;
}
