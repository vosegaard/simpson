/*
    Matrix and vector definition
    Copyright (C) 2009-2011 Zdenek Tosner

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

    Matrix definition and implementation.
*/

#ifndef __MATRIX_H
#define __MATRIX_H

#include "complx.h"


/* VECTORS are defined as array with N+1 elements, 1-based, element 0 contains number of elements (i.e. N)*/
#define LEN(v)   (*(int*)v)
int * int_vector(int len);
double * double_vector(int len);
complx * complx_vector(int len);
void free_int_vector(int *);
void free_double_vector(double *);
void free_complx_vector(complx *);
void cv_zero(complx * v);
void iv_zero(int * v);

/* MATRICES are structures that contain two representations:
 *    dense as 0-based column major
 *    sparse CSR format of MKL library
 *    Key for type:
 */
#define MAT_DENSE           1
#define MAT_DENSE_DIAG      2
#define MAT_SPARSE         10
#define MAT_SPARSE_DIAG    20

typedef struct _mat_complx {
       int type, row, col, basis;
       int *irow, *icol;
       complx * data;
} mat_complx;

typedef struct _mat_double {
       int type, row, col, basis;
       int *irow, *icol;
       double * data;
} mat_double;

/* Hamiltonians and propagators are square and potentially block-diagonal */
typedef struct _blk_mat_double {
	int dim, Nblocks, basis;
	int *blk_dims;
	mat_double *m;
} blk_mat_double;

typedef struct _blk_mat_complx {
	int dim, Nblocks, basis;
	int *blk_dims;
	mat_complx *m;
} blk_mat_complx;


/***
 * complex
 ***/
blk_mat_complx * create_blk_mat_complx(int matdim, int Nblocks, int *blk_dims, int type, int basis);
blk_mat_complx * create_blk_mat_complx_copy(blk_mat_complx *obj);
blk_mat_complx * create_blk_mat_complx_copy2(blk_mat_double *obj);
void free_blk_mat_complx(blk_mat_complx * obj);
mat_complx * complx_matrix(int row, int col, int type, int nnz, int basis);
void free_complx_matrix(mat_complx * obj);
void cm_change_nnz(mat_complx *cm, int nnz);
void create_sqmat_complx(mat_complx *obj, int dim, int type, int nnz, int basis);
void destroy_sqmat_complx(mat_complx * obj);

/***
 * double
 ***/
blk_mat_double * create_blk_mat_double(int matdim, int Nblocks, int *blk_dims, int type, int basis);
blk_mat_double * create_blk_mat_double_copy(blk_mat_double *obj);
void free_blk_mat_double(blk_mat_double * obj);
mat_double * double_matrix(int row, int col, int type, int nnz, int basis);
void free_double_matrix(mat_double * obj);
void dm_change_nnz(mat_double *cm, int nnz);
void create_sqmat_double(mat_double *obj, int dim, int type, int nnz, int basis);

char * matrix_type(int type);

#ifndef INTEL_MKL
#define MKL_INT int
#endif


#endif /* __MATRIX_H */

