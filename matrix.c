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

    Matrix definition and implementation. Data structure according to BLAS,
    i.e. Column major. Vector is just a single column.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "defs.h"
#include "complx.h"
#include "matrix.h"
#ifdef INTEL_MKL
#include "mkl.h"
#endif

/* vectors, 1-based indexing */
int * int_vector(int len)
{
   int * v;

   v = (int*)malloc((len+1)*sizeof(int));
   if (!v) {
      fprintf(stderr,"Can not allocate int vector (%i)\n",len);
      exit(1);
   }
   *(int*)v = len;
   return v;
}

double * double_vector(int len)
{
   double * v;

   v = (double*)malloc((len+1)*sizeof(double));
   if (!v) {
      fprintf(stderr,"Can not allocate double vector (%i)\n",len);
      exit(1);
   }
   *(int*)v = len;
   return v;
}

complx * complx_vector(int len)
{
   complx * v;

   v = (complx*)malloc((len+1)*sizeof(complx));
   if (!v) {
      fprintf(stderr,"Can not allocate complx vector (%i)\n",len);
      exit(1);
   }
   *(int*)v = len;
   return v;
}

void free_int_vector(int *v)
{
   free((char*)v);
}

void free_double_vector(double *v)
{
   free((char*)v);
}

void free_complx_vector(complx *v)
{
   free((char*)v);
}

void cv_zero(complx * v)
{
  //memset(&v[1],0,LEN(v)*sizeof(complx));
	memset(v+1,0,LEN(v)*sizeof(complx));
}

void iv_zero(int * v)
{
  //memset(&v[1],0,LEN(v)*sizeof(int));
	memset(v+1,0,LEN(v)*sizeof(int));
}


/****************************************
 *          C O M P L E X
 ****************************************/

/***
 *  nnz is used only for sparse matrices
 *  when nnz = 0, only irow is allocated
 ***/
void create_sqmat_complx(mat_complx *obj, int dim, int type, int nnz, int basis)
{
   assert(obj != NULL);
   if (dim == 1) type = MAT_DENSE_DIAG;
   obj->row = dim;
   obj->col = dim;
   obj->icol = NULL;
   obj->irow = NULL;
   obj->type = type;
   obj->basis = basis;
   switch (type) {
   case MAT_DENSE :
	   nnz = dim*dim;
	   obj->data = (complx*)malloc(nnz*sizeof(complx));
	   if (obj->data == NULL) {
		   fprintf(stderr,"Error: out of memory (create_sqmat_complx data)\n");
		   exit(1);
	   }
	   break;
   case MAT_DENSE_DIAG :
	   nnz = dim;
	   obj->data = (complx*)malloc(nnz*sizeof(complx));
	   if (obj->data == NULL) {
		   fprintf(stderr,"Error: out of memory (create_sqmat_complx data)\n");
		   exit(1);
	   }
	   break;
   case MAT_SPARSE :
   case MAT_SPARSE_DIAG :
		if (nnz != 0) {
			obj->data = (complx*)malloc(nnz*sizeof(complx));
			if (obj->data == NULL) {
				fprintf(stderr,"Error: out of memory (create_sqmat_complx data)\n");
				exit(1);
			}
			obj->icol = (MKL_INT*)calloc(nnz,sizeof(MKL_INT));
			if (obj->icol == NULL) {
				fprintf(stderr,"Error: out of memory (create_sqmat_complx icol)\n");
				exit(1);
			}
		} else {
			obj->data = NULL;
			obj->icol = NULL;
		}
		obj->irow = (MKL_INT*)calloc((dim+1),sizeof(MKL_INT));
		if (obj->irow == NULL) {
			fprintf(stderr,"Error: out of memory (create_sqmat_complx irow)\n");
			exit(1);
		}
		break;
   default :
	   fprintf(stderr,"Error: create_sqmat_complx - unknown type '%d'\n",type);
	   exit(1);
   }
}

void destroy_sqmat_complx(mat_complx * obj)
{
	assert(obj != NULL);
	free((char*)(obj->data)); obj->data = NULL;
	if (obj->irow != NULL) {
	   free((char*)(obj->irow));
	   obj->irow = NULL;
	}
	if (obj->icol != NULL) {
	   free((char*)(obj->icol));
	   obj->icol = NULL;
	}
}

void cm_change_nnz(mat_complx *cm, int nnz)
{
	assert(cm != NULL);
	if (cm->data != NULL) {
		if (cm->irow[cm->row]-1 != nnz) {
			complx * new_data = (complx*)realloc(cm->data, nnz*sizeof(complx));
			MKL_INT * new_icol = (MKL_INT*)realloc(cm->icol, nnz*sizeof(MKL_INT));
			if ( (new_data == NULL) || (new_icol == NULL) ) {
				fprintf(stderr,"Error: cm_change_nnz can not reallocate to given size\n");
				exit(1);
			}
			cm->icol = new_icol;
			cm->data = new_data;
		}
	} else {
		cm->icol = (MKL_INT *)calloc(nnz,sizeof(MKL_INT));
		cm->data = (complx*)malloc(nnz*sizeof(complx));
		if ( (cm->data == NULL) || (cm->icol == NULL) ) {
			fprintf(stderr,"Error: cm_change_nnz can not allocate given size\n");
			exit(1);
		}
	}
}

blk_mat_complx * create_blk_mat_complx(int matdim, int Nblocks, int *blk_dims, int type, int basis)
{
	blk_mat_complx *obj;
	int i;

	obj = (blk_mat_complx*)malloc(sizeof(blk_mat_complx));
	if (obj == NULL) {
		fprintf(stderr,"Error: out of memory (create_blk_mat_complx obj)\n");
		exit(1);
	}
	obj->dim = matdim;
	obj->Nblocks = Nblocks;
	obj->basis = basis;
	obj->m = (mat_complx*)malloc(Nblocks*sizeof(mat_complx));
	if (obj->m == NULL) {
		fprintf(stderr,"Error: out of memory (create_blk_mat_complx m)\n");
		exit(1);
	}
	if (Nblocks > 1) {
		obj->blk_dims = (int*)malloc(Nblocks*sizeof(int));
		if (obj->blk_dims == NULL) {
			fprintf(stderr,"Error: out of memory (create_blk_mat_complx blk_dims)\n");
			exit(1);
		}
		assert(Nblocks == LEN(blk_dims));
		const int dum_type = ( (type == MAT_DENSE_DIAG) || (type == MAT_SPARSE_DIAG) ) ? MAT_DENSE_DIAG : MAT_DENSE;
		for (i=0; i<Nblocks; i++) {
			obj->blk_dims[i] = blk_dims[i+1];
			if (obj->blk_dims[i] < MAXFULLDIM) {
				create_sqmat_complx(&(obj->m[i]),obj->blk_dims[i],dum_type,0,basis);
			} else {
				create_sqmat_complx(&(obj->m[i]),obj->blk_dims[i],type,1,basis);
			}
		}
	} else {
		obj->blk_dims = (int*)malloc(sizeof(int));
		if (obj->blk_dims == NULL) {
			fprintf(stderr,"Error: out of memory (create_blk_mat_complx blk_dims)\n");
			exit(1);
		}
		obj->blk_dims[0] = matdim;
		create_sqmat_complx(&(obj->m[0]),matdim,type,1,basis);
	}

	return obj;
}

blk_mat_complx * create_blk_mat_complx_copy(blk_mat_complx *obj)
{
	int i;
	blk_mat_complx *res = (blk_mat_complx*)malloc(sizeof(blk_mat_complx));
	if (res == NULL) {
		fprintf(stderr,"Error: out of memory (create_blk_mat_complx_copy res)\n");
		exit(1);
	}

	res->Nblocks = obj->Nblocks;
	res->dim = obj->dim;
	res->basis = obj->basis;
	res->blk_dims = (int*)malloc(res->Nblocks*sizeof(int));
	if (res->blk_dims == NULL) {
		fprintf(stderr,"Error: out of memory (create_blk_mat_complx_copy blk_dims)\n");
		exit(1);
	}
	res->m = (mat_complx*)malloc(res->Nblocks*sizeof(mat_complx));
	if (res->m == NULL) {
		fprintf(stderr,"Error: out of memory (create_blk_mat_complx_copy m)\n");
		exit(1);
	}
	for (i=0; i<res->Nblocks; i++) {
		res->blk_dims[i] = obj->blk_dims[i];
		create_sqmat_complx(&(res->m[i]),res->blk_dims[i],obj->m[i].type,1,res->basis);
	}

	return res;
}



blk_mat_complx * create_blk_mat_complx_copy2(blk_mat_double *obj)
{
	int i;
	blk_mat_complx *res = (blk_mat_complx*)malloc(sizeof(blk_mat_complx));
	if (res == NULL) {
		fprintf(stderr,"Error: out of memory (create_blk_mat_complx_copy2 res)\n");
		exit(1);
	}

	res->Nblocks = obj->Nblocks;
	res->dim = obj->dim;
	res->basis = obj->basis;
	res->blk_dims = (int*)malloc(res->Nblocks*sizeof(int));
	if (res->blk_dims == NULL) {
		fprintf(stderr,"Error: out of memory (create_blk_mat_complx_copy2 blk_dims)\n");
		exit(1);
	}
	res->m = (mat_complx*)malloc(res->Nblocks*sizeof(mat_complx));
	if (res->m == NULL) {
		fprintf(stderr,"Error: out of memory (create_blk_mat_complx_copy2 m)\n");
		exit(1);
	}
	for (i=0; i<res->Nblocks; i++) {
		res->blk_dims[i] = obj->blk_dims[i];
		create_sqmat_complx(&(res->m[i]),res->blk_dims[i],obj->m[i].type,1,res->basis);
	}

	return res;
}

void free_blk_mat_complx(blk_mat_complx * obj)
{
	assert(obj != NULL);
	int i;
	mat_complx *mx = obj->m;

	for (i=0; i<obj->Nblocks; i++) {
		destroy_sqmat_complx(mx);
		mx++;
	}
	free(obj->m);
	if (obj->blk_dims != NULL) free(obj->blk_dims);
	free(obj);
}




/***
 * double
 ***/

/***
 *  nnz is used only for sparse matrices
 *  when nnz = 0, only irow is allocated
 ***/
void create_sqmat_double(mat_double *obj, int dim, int type, int nnz, int basis)
{
   assert(obj != NULL);
   if (dim == 1) type = MAT_DENSE_DIAG;
   obj->row = dim;
   obj->col = dim;
   obj->icol = NULL;
   obj->irow = NULL;
   obj->type = type;
   obj->basis = basis;
   switch (type) {
   case MAT_DENSE :
	   nnz = dim*dim;
	   obj->data = (double*)malloc(nnz*sizeof(double));
	   if (obj->data == NULL) {
		   fprintf(stderr,"Error: out of memory (create_sqmat_double data)\n");
		   exit(1);
	   }
	   break;
   case MAT_DENSE_DIAG :
	   nnz = dim;
	   obj->data = (double*)malloc(nnz*sizeof(double));
	   if (obj->data == NULL) {
		   fprintf(stderr,"Error: out of memory (create_sqmat_double data)\n");
		   exit(1);
	   }
	   break;
   case MAT_SPARSE :
   case MAT_SPARSE_DIAG :
		if (nnz != 0) {
			obj->data = (double*)malloc(nnz*sizeof(double));
			if (obj->data == NULL) {
				fprintf(stderr,"Error: out of memory (create_sqmat_double data)\n");
				exit(1);
			}
			obj->icol = (MKL_INT*)calloc(nnz,sizeof(MKL_INT));
			if (obj->icol == NULL) {
				fprintf(stderr,"Error: out of memory (create_sqmat_double icol)\n");
				exit(1);
			}
		} else {
			obj->data = NULL;
			obj->icol = NULL;
		}
		obj->irow = (MKL_INT*)calloc((dim+1),sizeof(MKL_INT));
		if (obj->irow == NULL) {
			fprintf(stderr,"Error: out of memory (create_sqmat_double irow)\n");
			exit(1);
		}
		break;
   default :
	   fprintf(stderr,"Error: create_sqmat_double - unknown type '%d'\n",type);
	   exit(1);
   }
}

void destroy_sqmat_double(mat_double * obj)
{
	assert(obj != NULL);
	free((char*)(obj->data)); obj->data = NULL;
	if (obj->irow != NULL) {
	   free((char*)(obj->irow));
	   obj->irow = NULL;
	}
	if (obj->icol != NULL) {
	   free((char*)(obj->icol));
	   obj->icol = NULL;
	}
}


void dm_change_nnz(mat_double *cm, int nnz)
{
	assert(cm != NULL);
	if (cm->data != NULL) {
		if (cm->irow[cm->row]-1 != nnz) {
			double * new_data = (double*)realloc(cm->data, nnz*sizeof(double));
			MKL_INT * new_icol = (MKL_INT*)realloc(cm->icol, nnz*sizeof(MKL_INT));
			if ( (new_data == NULL) || (new_icol == NULL) ) {
				fprintf(stderr,"Error: dm_change_nnz can not reallocate to given size\n");
				exit(1);
			}
			cm->icol = new_icol;
			cm->data = new_data;
		}
	} else {
		cm->icol = (MKL_INT *)calloc(nnz,sizeof(MKL_INT));
		cm->data = (double*)malloc(nnz*sizeof(double));
		if ( (cm->data == NULL) || (cm->icol == NULL) ) {
			fprintf(stderr,"Error: dm_change_nnz can not allocate given size (%d)\n",nnz);
			exit(1);
		}
	}
}

blk_mat_double * create_blk_mat_double(int matdim, int Nblocks, int *blk_dims, int type, int basis)
{
	blk_mat_double *obj;
	int i;

	obj = (blk_mat_double*)malloc(sizeof(blk_mat_double));
	if (obj == NULL) {
		fprintf(stderr,"Error: out of memory (create_blk_mat_double obj)\n");
	}
	obj->dim = matdim;
	obj->Nblocks = Nblocks;
	obj->basis = basis;
	obj->m = (mat_double*)malloc(Nblocks*sizeof(mat_double));
	if (obj->m == NULL) {
		fprintf(stderr,"Error: out of memory (create_blk_mat_double m)\n");
		exit(1);
	}
	if (Nblocks > 1) {
		obj->blk_dims = (int*)malloc(Nblocks*sizeof(int));
		if (obj->blk_dims == NULL) {
			fprintf(stderr,"Error: out of memory (create_blk_mat_double blk_dims)\n");
			exit(1);
		}
		assert(Nblocks == LEN(blk_dims));
		const int dum_type = ( (type == MAT_DENSE_DIAG) || (type == MAT_SPARSE_DIAG) ) ? MAT_DENSE_DIAG : MAT_DENSE;
		for (i=0; i<Nblocks; i++) {
			obj->blk_dims[i] = blk_dims[i+1];
			if (obj->blk_dims[i] < MAXFULLDIM) {
				create_sqmat_double(&(obj->m[i]),obj->blk_dims[i],dum_type,0,basis);
			} else {
				create_sqmat_double(&(obj->m[i]),obj->blk_dims[i],type,1,basis);
			}
		}
	} else {
		obj->blk_dims = (int*)malloc(sizeof(int));
		if (obj->blk_dims == NULL) {
			fprintf(stderr,"Error: out of memory (create_blk_mat_double blk_dims)\n");
			exit(1);
		}
		obj->blk_dims[0] = matdim;
		create_sqmat_double(&(obj->m[0]),matdim,type,1,basis);
	}

	return obj;
}

blk_mat_double * create_blk_mat_double_copy(blk_mat_double *obj)
{
	int i;
	blk_mat_double *res = (blk_mat_double*)malloc(sizeof(blk_mat_double));
	if (res == NULL) {
		fprintf(stderr,"Error: out of memory (create_blk_mat_double_copy res)\n");
		exit(1);
	}

	res->Nblocks = obj->Nblocks;
	res->dim = obj->dim;
	res->basis = obj->basis;
	res->blk_dims = (int*)malloc(res->Nblocks*sizeof(int));
	if (res->blk_dims == NULL) {
		fprintf(stderr,"Error: out of memory (create_blk_mat_double_copy blk_dims)\n");
		exit(1);
	}
	res->m = (mat_double*)malloc(res->Nblocks*sizeof(mat_double));
	if (res->m == NULL) {
		fprintf(stderr,"Error: out of memory (create_blk_mat_double_copy m)\n");
		exit(1);
	}
	for (i=0; i<res->Nblocks; i++) {
		res->blk_dims[i] = obj->blk_dims[i];
		create_sqmat_double(&(res->m[i]),res->blk_dims[i],obj->m[i].type,1,res->basis);
	}

	return res;
}


void free_blk_mat_double(blk_mat_double * obj)
{
	assert(obj != NULL);
	int i;
	mat_double *mx = obj->m;

	for (i=0; i<obj->Nblocks; i++) {
		destroy_sqmat_double(mx);
		mx++;
	}
	free(obj->m);
	if (obj->blk_dims != NULL) free(obj->blk_dims);
	free(obj);
}


/**********************************/

char * matrix_type(int type)
{
	char *res = (char*)malloc(16*sizeof(char));

	switch (type) {
	case MAT_DENSE :
		sprintf(res,"%s","dense general");
		break;
	case MAT_DENSE_DIAG :
		sprintf(res,"%s","dense diagonal");
		break;
	case MAT_SPARSE :
		sprintf(res,"%s","sparse general");
		break;
	case MAT_SPARSE_DIAG :
		sprintf(res,"%s","sparse diagonal");
		break;
   default :
	   fprintf(stderr,"Error: matrix_type - unknown type '%d'\n",type);
	   exit(1);
	}
	return res;
}


/* back-compatibility functions during development */
mat_complx * complx_matrix(int row, int col, int type, int nnz, int basis)
{
	   mat_complx * obj;

	   if ( (type == MAT_DENSE_DIAG || type == MAT_SPARSE_DIAG) && (row != col) ) {
		   fprintf(stderr,"Error: complx_matrix - request for non-square diagonal matrix (%d,%d)\n",row,col);
		   exit(1);
	   }
	   if (row == 1 && col == 1) type = MAT_DENSE_DIAG;

	   obj = (mat_complx*)(malloc(sizeof(mat_complx)));
		if (obj == NULL) {
			fprintf(stderr,"Error: out of memory (complx_matrix obj)\n");
			exit(1);
		}
	   obj->row = row;
	   obj->col = col;
	   obj->icol = NULL;
	   obj->irow = NULL;
	   obj->type = type;
	   obj->basis = basis;
	   switch (type) {
	   case MAT_DENSE :
		   nnz = row*col;
		   obj->data = (complx*)malloc(nnz*sizeof(complx));
			if (obj->data == NULL) {
				fprintf(stderr,"Error: out of memory (complx_matrix data)\n");
				exit(1);
			}
		   break;
	   case MAT_DENSE_DIAG :
		   nnz = row;
		   obj->data = (complx*)malloc(nnz*sizeof(complx));
			if (obj->data == NULL) {
				fprintf(stderr,"Error: out of memory (complx_matrix data)\n");
				exit(1);
			}
		   break;
	   case MAT_SPARSE :
	   case MAT_SPARSE_DIAG :
			if (nnz) {
				obj->data = (complx*)malloc(nnz*sizeof(complx));
				if (obj->data == NULL) {
					fprintf(stderr,"Error: out of memory (complx_matrix data)\n");
					exit(1);
				}
				obj->icol = (MKL_INT*)calloc(nnz,sizeof(MKL_INT));
				if (obj->icol == NULL) {
					fprintf(stderr,"Error: out of memory (complx_matrix icol)\n");
					exit(1);
				}
			} else {
				obj->data = NULL;
				obj->icol = NULL;
			}
			obj->irow = (MKL_INT*)calloc((row+1),sizeof(MKL_INT));
			if (obj->irow == NULL) {
				fprintf(stderr,"Error: out of memory (complx_matrix irow)\n");
				exit(1);
			}
			break;
	   default :
		   fprintf(stderr,"Error: complx_matrix - unknown type '%d'\n",type);
		   exit(1);
	   }
	   return obj;
}

void free_complx_matrix(mat_complx * obj)
{
   free((char*)(obj->data)); obj->data = NULL;
   if (obj->irow) {
	   free((char*)(obj->irow));
	   obj->irow = NULL;
   }
   if (obj->icol) {
	   free((char*)(obj->icol));
	   obj->icol = NULL;
   }
   obj->row = obj->col = 0;
   obj->type = -1;
   free((char*)obj);
}

mat_double * double_matrix(int row, int col, int type, int nnz, int basis)
{
   mat_double * obj;

   if ( (type == MAT_DENSE_DIAG || type == MAT_SPARSE_DIAG) && (row != col) ) {
	   fprintf(stderr,"Error: double_matrix - request for non-square diagonal matrix (%d,%d)\n",row,col);
	   exit(1);
   }
   if (row == 1 && col == 1) type = MAT_DENSE_DIAG;

   obj = (mat_double*)(malloc(sizeof(mat_double)));
	if (obj == NULL) {
		fprintf(stderr,"Error: out of memory (double_matrix obj)\n");
		exit(1);
	}
   obj->row = row;
   obj->col = col;
   obj->icol = NULL;
   obj->irow = NULL;
   obj->type = type;
   obj->basis = basis;
   switch (type) {
   case MAT_DENSE :
	   nnz = row*col;
	   obj->data = (double*)malloc(nnz*sizeof(double));
		if (obj->data == NULL) {
			fprintf(stderr,"Error: out of memory (double_matrix data)\n");
			exit(1);
		}
	   break;
   case MAT_DENSE_DIAG :
	   nnz = row;
	   obj->data = (double*)malloc(nnz*sizeof(double));
		if (obj->data == NULL) {
			fprintf(stderr,"Error: out of memory (double_matrix data)\n");
			exit(1);
		}
	   break;
   case MAT_SPARSE :
   case MAT_SPARSE_DIAG :
		obj->irow = (MKL_INT*)calloc((row+1),sizeof(MKL_INT));
		if (obj->irow == NULL) {
			fprintf(stderr,"Error: out of memory (double_matrix irow)\n");
			exit(1);
		}
		if (nnz) {
			obj->data = (double*)malloc(nnz*sizeof(double));
			if (obj->data == NULL) {
				fprintf(stderr,"Error: out of memory (double_matrix data)\n");
				exit(1);
			}
			obj->icol = (MKL_INT*)calloc(nnz,sizeof(MKL_INT));
			if (obj->icol == NULL) {
				fprintf(stderr,"Error: out of memory (double_matrix icol)\n");
				exit(1);
			}
		} else {
			obj->data = NULL;
			obj->icol = NULL;
		}
		break;
   default :
	   fprintf(stderr,"Error: double_matrix - unknown type '%d'\n",type);
	   exit(1);
   }

   return obj;
}

void free_double_matrix(mat_double * obj)
{
   free((char*)(obj->data)); obj->data = NULL;
   if (obj->irow) {
	   free((char*)(obj->irow));
	   obj->irow = NULL;
   }
   if (obj->icol) {
	   free((char*)(obj->icol));
	   obj->icol = NULL;
   }
   obj->row = obj->col = 0;
   obj->type = -1;
   free((char*)obj);
}

