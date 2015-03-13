/*
 * cm.h
 *
 *  Created on: Jun 11, 2010
 *      Author: zdenek
 */

#ifndef CM_H_
#define CM_H_

#include "complx.h"
#include "matrix.h"

/* vector operations */
void dv_muld(double *vec, double d);
void dv_multod(double *v1, double *v2, double d);
void dv_zero(double *v);
complx * cv_dup(complx *vec);
void cv_copy(complx *dest, complx *vec);
void cv_muld(complx *vec, double d);
void cv_multod(complx *v1, complx *v2, double d);
void cv_mulc(complx *vec, complx z);
void cv_multoc(complx *v1, complx *v2, complx z);
void cv_conj(complx *vec);
double cv_asum(complx *vec);
void cv_matmulto(complx *vec, mat_complx *A);
void cv_matmul(complx *res, mat_complx *A, complx *vec);
complx cv_dotu(complx *a, complx *b);
complx cv_dotc(complx *a, complx *b);
double cv_norm(complx *vec);
int iv_max(int *vec);
int * iv_dup(int *vec);
void cv_print(complx *vec, char *title);

/* double matrix operations */
void dm_zero(mat_double *m);
mat_double * dm_dup(mat_double *m);
void dm_copy(mat_double *dest, mat_double *from);
void dm_sparse(mat_double *m, double tol);
void dm_dense(mat_double *m);
mat_double * dm_creatediag(const char c, int dim, double d,int basis);
mat_complx * dm_complx(mat_double *m);
mat_complx * dm_imag(mat_double *m);
void dm_copy2cm(mat_double *dm, mat_complx *cm);
void dm_copy2cm_imag(mat_double *dm, mat_complx *cm);
void dm_muld(mat_double *m, double d);
void dm_multod(mat_double *m1, mat_double *m2, double d);
void dm_multocc(mat_double *res, mat_complx *m, complx z);
mat_double * dm_mul(mat_double *m1, mat_double *m2);
void dm_multo(mat_double *m1, mat_double *m2);
void dm_multo_rev(mat_double *m1, mat_double *m2);
void dm_addto(mat_double *m1, mat_double *m2);
void dm_addtodiag(mat_double *m, double d);
void dm_shrinktodiag(mat_double *m);
/* return matrix element, 1-based indexing */
double dm_getelem(mat_double *m, int r, int c);
void dm_print(mat_double *m, char *title);
double dm_normest(mat_double *m);
int dm_nnz(mat_double *m);
void dm_nnz2one(mat_double *m);
void dm_mm(double alpha, mat_double *A, mat_double *B, double beta, mat_double *C);

/* complex martix operations */
mat_complx * cm_creatediag(const char c, int dim, complx d, int basis);
mat_complx *cm_dup(mat_complx *m);
void cm_copy(mat_complx *m1, mat_complx *m2);
void cm_unit(mat_complx *m);
void cm_zero(mat_complx *m);
void cm_filter(mat_complx *m, mat_complx *mask);
mat_complx * cm_adjoint(mat_complx *m);
void cm_adjointi(mat_complx *m);
mat_complx * cm_transpose(mat_complx *m);
int cm_ishermit(mat_complx *m);
void cm_dense(mat_complx *m);
void cm_dense_full(mat_complx *m);
void cm_sparse(mat_complx *m, double tol);
void cm_shrinktodiag(mat_complx *m);
complx cm_getelem(mat_complx *m, int r, int c);
mat_double * cm_real(mat_complx *m);
mat_double *cm_realdiag(mat_complx *m); // will return just diagonal no matter what m is
mat_double * cm_imag(mat_complx *m);
void cm_print(mat_complx *m, char * title);
void cm_nnz2one(mat_complx *m);
void cm_muld(mat_complx *m, double d);
void cm_mulc(mat_complx *m, complx z);
void cm_multod(mat_complx *res, mat_complx *m, double d);
void cm_multoc(mat_complx *res, mat_complx *m, complx z);
void cm_multocr(mat_complx *res, mat_double *m, complx z);
void cm_addto(mat_complx *res, mat_complx *m);
void cm_subfrom(mat_complx *res, mat_complx *m);
mat_complx * cm_mul(mat_complx *m1, mat_complx *m2);
void cm_multo(mat_complx *m1, mat_complx *m2);
void cm_multo_rev(mat_complx *m1, mat_complx *m2);
mat_complx * cm_commutator(mat_complx *a, mat_complx *b);
complx cm_trace(mat_complx *a, mat_complx *b);
complx cm_trace_adjoint(mat_complx *a, mat_complx *b);
mat_complx * cm_ln(mat_complx *m);
void cm_diag(mat_complx *m, complx *eigs, mat_complx *T);
mat_complx * cm_direct(mat_complx *m1, mat_complx *m2);
mat_complx * cm_directadd(mat_complx *m1, mat_complx *m2);
double cm_sumnorm1(mat_complx *m);
double cm_normest(mat_complx *m);
int cm_nnz(mat_complx *m);
void cm_vec_on_diag(mat_complx *m, complx *vec, int basis);
mat_complx * cm_expm_pade(mat_complx *m);
mat_complx * cm_power(mat_complx *m, int n);
mat_complx * simpson_zcsrmultcsr(mat_complx *a, mat_complx *b);
void cm_swap_innards_and_destroy(mat_complx *dest, mat_complx *from);

/* transformations and propagators */
void simtrans(mat_complx *m, mat_complx *T);
void simtrans_adj(mat_complx *m, mat_complx *T);
void simtrans_zrot2(mat_complx *m, mat_double *T);
void prop_real(mat_complx *prop, mat_double *ham, double dt, int method);
void prop_complx(mat_complx *prop, mat_complx *ham, double dt, int method);
void prop_krylov(double dt, mat_complx *A, complx *vec);

void dm_transposei(mat_double *m);
void cm_permute(mat_complx *m, int *permvec);
void dm_permute(mat_double *m, int *permvec);
complx cm_true_trace(mat_complx *m);

#endif /* CM_H_ */
