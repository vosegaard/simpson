/*
 * blockdiag.h
 *
 *  Created on: 12.3.2011
 *      Author: zdenek
 */

#ifndef BLOCKDIAG_H_
#define BLOCKDIAG_H_

/* block diagonal matrix routines and auxiliary */

#include "complx.h"
#include "matrix.h"
#include "sim.h"


void blk_dm_zero(blk_mat_double *obj);
void blk_cm_zero(blk_mat_complx *obj);
void blk_cm_unit(blk_mat_complx *obj);
void blk_dm_multod_diag(blk_mat_double *mblk, mat_double *mdiag, double d);
void blk_dm_multod(blk_mat_double *res, blk_mat_double *m, double d);
mat_double * dm_get_diagblock_permute(mat_double *a, int *permvec, int sft, int dim);
blk_mat_double * blk_dm_dup(blk_mat_double *m);
void blk_dm_copy(blk_mat_double *res, blk_mat_double *m);
void blk_dm_muld(blk_mat_double *m, double d);
void blk_dm_print(blk_mat_double *obj, char *title);
void blk_cm_print(blk_mat_complx *obj, char *title);
mat_double * dm_change_basis_2(mat_double *m, int basis, Sim_info *sim);
mat_complx * cm_change_basis_2(mat_complx *m, int basis, Sim_info *sim);
double blk_dm_getelem(blk_mat_double *obj, int row, int col);
complx blk_cm_getelem(blk_mat_complx *obj, int row, int col);
void blk_dm_change_basis(blk_mat_double *dest, blk_mat_double *from, Sim_info *sim);
void blk_cm_change_basis(blk_mat_complx *dest, blk_mat_complx *from, Sim_info *sim);
void blk_simtrans(mat_complx *sigma, blk_mat_complx *U, Sim_info *sim);
void blk_simtrans_adj(mat_complx *sigma, blk_mat_complx *U, Sim_info *sim);
blk_mat_complx * blk_cm_adjoint(blk_mat_complx *m);
blk_mat_complx * blk_cm_dup(blk_mat_complx *m);
void blk_cm_copy(blk_mat_complx *dest, blk_mat_complx *from);
blk_mat_complx * blk_cm_mul(blk_mat_complx *blkm1, blk_mat_complx *blkm2);
mat_complx * blk_cm_mul_v1(blk_mat_complx *blkm, mat_complx *cm);
mat_complx * blk_cm_mul_v2(mat_complx *cm, blk_mat_complx *blkm);
void update_propagator(blk_mat_complx *U, blk_mat_complx *dU, Sim_info *sim, Sim_wsp *wsp);
void blk_prop_real(blk_mat_complx *U, blk_mat_double *ham, double duration, Sim_info *sim);
void blk_prop_complx(blk_mat_complx *U, mat_complx *ham, double duration, Sim_info *sim);
void blk_prop_complx_2(blk_mat_complx *U, mat_complx *mtx, double dur, int propmethod);
void blk_dm_multod_extract(blk_mat_double *blkm, mat_double *dm, double dval);
void blk_cm_multod_extract(blk_mat_complx *blkm, mat_complx *dm, double dval);
void blk_simtrans_zrot2(blk_mat_complx *U, mat_double *sumUph);
void blk_change_structure_double(blk_mat_double *blkm, blk_mat_double *ham);
void blk_change_structure_double_nondiag(blk_mat_double *blkm, blk_mat_double *ham, int issparse);
void blk_change_structure_complx(blk_mat_complx *blkm, blk_mat_double *ham);
void blk_change_structure_complx2(blk_mat_complx *blkm, blk_mat_complx *cmx);
void blk_dm_copy_2(mat_double *dest, blk_mat_double *blkm);
void blk_cm_copy_2(mat_complx *dest, blk_mat_complx *blkm);
void change_basis_pulse(Sim_info *sim, Sim_wsp *wsp, int basis);
void change_basis_delay(Sim_info *sim, Sim_wsp *wsp);
int blk_cm_isdiag(blk_mat_complx *blkm);
blk_mat_complx * blk_cm_diag(blk_mat_complx *Ud);
blk_mat_complx * blk_cm_diag_sort(blk_mat_complx *Ud);
blk_mat_complx * blk_cm_diag_partly(blk_mat_complx *Ud);
blk_mat_complx * blk_cm_power(blk_mat_complx *blkm, int n);
complx blk_cm_trace_adjoint(blk_mat_complx *A, blk_mat_complx *B, Sim_info *sim);
int blk_dm_nnz(blk_mat_double *blkm);
int blk_cm_nnz(blk_mat_complx *blkm);
mat_complx * cm_extract_block(mat_complx *cm, int r0, int c0, int Nr, int Nc);
mat_double * dm_extract_block(mat_double *cm, int r0, int c0, int Nr, int Nc);
void cm_addto_block(mat_complx *cm, mat_complx *bm, int r0, int c0);
blk_mat_complx * blk_cm_ln(blk_mat_complx *blkm);

#endif /* BLOCKDIAG_H_ */
