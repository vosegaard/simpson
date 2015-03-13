/*
    Spinsystem routines
    Copyright (C) 1999 Mads Bak
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
    
    Routines for calculation of spin-tensors for the spin-system.
    Used by readsys.c that creates the spinsystem, and 
    sim.c to create the start and detect operators.
*/

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <math.h>
#include <assert.h>
#include "matrix.h"
#include "cm.h"
#include "blockdiag.h"
#include "spinsys.h"
#include "tclutil.h"
#include "ham.h"
#include "defs.h"

double ss_gamma1H()
{
   if (isotopes[0].number != 1) {
     fprintf(stderr,"error: isotope list does not have 1H as first entry\n");
     exit(1);
   }
   return isotopes[0].gamma;
}

void ss_showspins(SpinSys* S)
{
  int i;
  
  for (i=1;i<=S->nspins;i++) {
    printf("Spin[%d] : I = %g\n",i,S->iso[i]->spin);
  }
}

ISOTOPE* ss_findisotope(char* name)
{
  char buf[256],*src,*dst;
  int number,niso, gotnumber=0;
  
  if (!name) {
    fprintf(stderr,"ss_findisotope: argument must not be NULL\n");
    exit(1);
  }
  src=name;
  dst=buf;
    
  while (isdigit(*src) && *src != 0) {
    *dst++ = *src++;
  }

  if (src == name) {
    fprintf(stderr,"error: name of isotope '%s' is wrong, '23Na' is a correct example\n",name);
    exit(1);
  }

  *dst=0;
  number=atoi(buf);

  niso=0;
  while (isotopes[niso].number) {
    if (isotopes[niso].number == number) {
      gotnumber=1;
      sprintf(buf,"%d%s",isotopes[niso].number,isotopes[niso].name);
      if (!strcmp(name,buf)) {
/*
        fprintf(stderr,"error: name of isotope '%s' is wrong, maybe '%s' is correct ?\n",name,buf);
        exit(1);
*/
        return &isotopes[niso];   
      }
    }
    niso++;
  }
  fprintf(stderr,"error: unable to find name '%s' in the internal isotope table.\n", name);
  if (gotnumber==1) {
    niso=0;
    fprintf(stderr, "       My suggestions are ");
    while (isotopes[niso].number) {
      if (isotopes[niso].number == number)
        fprintf(stderr, "'%d%s' ", number, isotopes[niso].name);
      if (!strcmp(isotopes[niso].name, src))
        fprintf(stderr, "'%d%s' ", isotopes[niso].number, isotopes[niso].name);
      niso++;
    }
    fprintf(stderr, "\n");
  }
  exit(1);
}

double ss_qn(SpinSys* S,int spin)
{ 
  if (spin < 1 || spin > S->nspins) {
    fprintf(stderr,"error: invalid range of spin (was %d)\n",spin);
    exit(1);
  }
  return S->iso[spin]->spin;
}

double ss_gamma(SpinSys* S,int spin)
{ 
  if (spin < 1 || spin > S->nspins) {
    fprintf(stderr,"error: invalid range of spin (was %d)\n",spin);
    exit(1);
  }
  return S->iso[spin]->gamma;
}

ISOTOPE* ss_isotope(SpinSys* S,int spin)
{ 
  if (spin < 1 || spin > S->nspins) {
    fprintf(stderr,"error: invalid range of spin (was %d)\n",spin);
    exit(1);
  }
  return S->iso[spin];
}

int ss_matdim(SpinSys* S)
{
  int i;
  double N;

  if (S->nspins < 1) return 0;
  N = 1.0;
  for (i=1;i<=S->nspins;i++) 
    N *= (double)(2.0*(S->iso[i]->spin)+1.0);
  return (int)N;
}


void ss_addspin(SpinSys* S,char* name)
{
  S->nspins++;

  if (S->nspins >= MAXSPINS) {
    fprintf(stderr,"error: max number of spins reached, increase MAXSPINS\n");
    exit(1);
  }
  S->iso[S->nspins]=ss_findisotope(name);  
  S->matdim=ss_matdim(S);
}

void ss_initialize(SpinSys* S)
{
  S->nspins=0;
  S->nchan=0;
}

int ss_issame(SpinSys* S,int spin1,int spin2)
{
  if (spin1 < 1 || spin1 > S->nspins) {
    fprintf(stderr,"error: invalid range of spin (was %d)\n",spin1);
    exit(1);
  }
  if (spin2 < 1 || spin2 > S->nspins) {
    fprintf(stderr,"error: invalid range of spin (was %d)\n",spin2);
    exit(1);
  }
  return (S->iso[spin1] == S->iso[spin2]);
}

mat_complx * _Ip(double I)
{
  int i,j,S;  
  double mm,m;
  mat_complx * M;
  
  S= (int)(2*I+1);
  M=complx_matrix(S,S,MAT_DENSE,0,0);
  for (i=0;i<S;i++) {
    mm= I-i;
    for (j=0;j<S;j++) {
      m= I-j;
      if (mm == m + 1)
        M->data[i+j*S]=Complx(sqrt(I*(I+1.0)-m*(m+1.0)),0.0);
      else
        M->data[i+j*S]=Complx(0.0,0.0);
    }
  }
  return M;
}

mat_complx * _Im(double I)
{
  int i,j,S;  
  double mm,m;
  mat_complx * M;
  
  S= (int)(2*I+1);
  M=complx_matrix(S,S,MAT_DENSE,0,0);
  for (i=0;i<S;i++) {
    mm= I-i;
    for (j=0;j<S;j++) {
      m= I-j;
      if (mm == m - 1)
        M->data[i+j*S]=Complx(sqrt(I*(I+1.0)-m*(m-1.0)),0.0);
      else
        M->data[i+j*S]=Complx(0.0,0.0);
    }
  }
  return M;
}

/***** moved to Iz directly *****
mat_complx * _Iz(double I)
{
  int i,j,S;  
  double mm,m;
  mat_complx * M;
  
  S= (int)(2*I+1);
  M=complx_matrix(S,S,MAT_DENSE,0);
  for (i=0;i<S;i++) {
    mm= I-i;
    for (j=0;j<S;j++) {
      m= I-j;
      if (mm == m)
        M->data[i+j*S]=Complx(m,0.0);
      else
        M->data[i+j*S]=Complx(0.0,0.0);
    }
  }
  return M;
}
*/

mat_complx * _Iq(double I)
{
  int i,j,S;  
  double mm,m;
  mat_complx * M;
  
  S= (int)(2*I+1);
  M=complx_matrix(S,S,MAT_DENSE,0,0);
  for (i=0;i<S;i++) {
    mm= I-i;
    for (j=0;j<S;j++) {
      m= I-j;
      M->data[i+j*S]=Complx(m,mm);
    }
  }
  return M;
}

/* this probably wont work in multi-spin cases. */
mat_complx * _Ic(double I)
{
  int S;
  mat_complx * M;
  /* ZT: adding test for I>1 and half-integer spin */
  if ( (I<1.4) || ( fabs(I-floor(I)-0.5)>0.01 ) ) {
     fprintf(stderr,"error in central transition operator, spin=%f is not allowed\n",I);
     exit(1);
  }

  S= (int)(2*I+1);
  M = complx_matrix(S,S,MAT_DENSE,0,0);
  cm_zero(M);
  M->data[S/2-1+(S/2)*S]=Complx(sqrt(I*(I+1)+0.25),0.0);
  return M;
}

mat_complx * _Ix(double I)
{
/* return (Ip(I)+Im(I))*complx(0.5,0); */
  mat_complx *ip,*im;

  ip=_Ip(I);
  im=_Im(I);
  cm_addto(ip,im);
  cm_muld(ip,0.5);
  free_complx_matrix(im);
  return ip;
}

mat_complx * _Iy(double I)
{
/* return (Ip(I)-Im(I))*complx(0.0,-5.0); */
  mat_complx *ip,*im;

  ip=_Ip(I);
  im=_Im(I);
  cm_subfrom(ip,im);
  cm_mulc(ip,Complx(0.0,-0.5));
  free_complx_matrix(im);
  return ip;
}

/* creates the direct product, for example creates
  for a four spin-system with spin=3 the matrix:
    1 (x) 1 (x) m (x) 1
 ****
 *  ZT: new and faster version of In (not using direct product)
 ****/
mat_complx * In(Sim_info* sim, mat_complx * m,int spin)
{
	  int k,l,i,j,r,c,K,L,N,matdim;
	  complx z;
	  SpinSys * S;
	  mat_complx *res;

	  assert(m->row == m->col);
	  S = sim->ss;
	  /*printf("spin = %d\n",spin);
	  cm_print(m,"matrix is");*/
	  N = m->row;
	  matdim = S->matdim;
	  K = 1;
	  for (i=1;i<spin;i++) {
	     K *= (int)(2.0*ss_qn(S,i)+1.0);
	  }
	  L = 1;
	  for (i=spin+1; i<=S->nspins; i++) {
	     L *= (int)(2.0*ss_qn(S,i)+1.0);
	  }
	  /*printf("sp_In report: K = %d; N = %d; L = %d; matdim = %d\n",K,N,L,matdim);*/
	  if (K*N*L != matdim) {
	     fprintf(stderr,"In error: mismatch in dimensions (%d * %d * %d != %d)\n",K,N,L,matdim);
	     exit(1);
	  }

	  if (sim->sparse) {
		  const int nnz = cm_nnz(m)*K*L;
		  res = complx_matrix(matdim,matdim,MAT_SPARSE,nnz,0);
		  complx *zz = res->data;
		  int *ic = res->icol;
		  res->irow[0] = 1;
		  r = 0;
		  for (k=0; k<K; k++) {
			  for (i=0; i<N; i++) {
				  for (l=0; l<L; l++) {
					  res->irow[r+1] = res->irow[r];
					  for (j=0; j<N; j++) {
						  z = m->data[i+j*N];
						  if ( fabs(z.re) < SPARSE_TOL && fabs(z.im) < SPARSE_TOL) continue;
						  *zz = z;
						  *ic = j*L + l + 1 + k*N*L;
						  (res->irow[r+1])++;
						  zz++;
						  ic++;
					  }
					  r++;
				  }
			  }
		  }
		  assert(res->irow[matdim] == nnz+1);
	  } else {
		  res = complx_matrix(matdim,matdim,MAT_DENSE,0,0);
		  cm_zero(res);
		  for (i=0; i<N; i++) {
			  for (j=0; j<N; j++) {
				  z = m->data[i+j*N];
				  if ( fabs(z.re) < SPARSE_TOL && fabs(z.im) < SPARSE_TOL) continue;
				  for (k=0; k<K; k++) {
					  for (l=0; l<L; l++) {
						  r = i*L +l + k*L*N;
						  c = j*L +l + k*L*N;
						  res->data[r+c*matdim] = z;
					  }
				  }
			  }
		  }
	  }

	  return res;
}

mat_double * In_real(Sim_info* sim, mat_double * m,int spin)
{
	  int k,l,i,j,r,c,K,L,N,matdim;
	  double z;
	  SpinSys * S;
	  mat_double *res;

	  assert(m->row == m->col);
	  S = sim->ss;
	  /*printf("spin = %d\n",spin);
	  cm_print(m,"matrix is");*/
	  N = m->row;
	  matdim = S->matdim;
	  K = 1;
	  for (i=1;i<spin;i++) {
	     K *= (int)(2.0*ss_qn(S,i)+1.0);
	  }
	  L = 1;
	  for (i=spin+1; i<=S->nspins; i++) {
	     L *= (int)(2.0*ss_qn(S,i)+1.0);
	  }
	  /*printf("sp_In report: K = %d; N = %d; L = %d; matdim = %d\n",K,N,L,matdim);*/
	  if (K*N*L != matdim) {
	     fprintf(stderr,"In_real error: mismatch in dimensions (%d * %d * %d != %d)\n",K,N,L,matdim);
	     exit(1);
	  }

	  if (sim->sparse) {
		  const int nnz = dm_nnz(m)*K*L;
		  res = double_matrix(matdim,matdim,MAT_SPARSE,nnz,0);
		  double *zz = res->data;
		  int *ic = res->icol;
		  res->irow[0] = 1;
		  r = 0;
		  for (k=0; k<K; k++) {
			  for (i=0; i<N; i++) {
				  for (l=0; l<L; l++) {
					  res->irow[r+1] = res->irow[r];
					  for (j=0; j<N; j++) {
						  z = m->data[i+j*N];
						  if ( fabs(z) < SPARSE_TOL ) continue;
						  *zz = z;
						  *ic = j*L + l + 1 + k*N*L;
						  (res->irow[r+1])++;
						  zz++;
						  ic++;
					  }
					  r++;
				  }
			  }
		  }
		  assert(res->irow[matdim] == nnz+1);
	  } else {
		  res = double_matrix(matdim,matdim,MAT_DENSE,0,0);
		  dm_zero(res);
		  for (i=0; i<N; i++) {
			  for (j=0; j<N; j++) {
				  z = m->data[i+j*N];
				  if ( fabs(z) < SPARSE_TOL ) continue;
				  for (k=0; k<K; k++) {
					  for (l=0; l<L; l++) {
						  r = i*L +l + k*L*N;
						  c = j*L +l + k*L*N;
						  res->data[r+c*matdim] = z;
					  }
				  }
			  }
		  }
	  }

	  return res;
}


mat_complx * Icoherence(Sim_info* sim, double* coh)
{
  SpinSys *S=sim->ss;
  int i,j,k,N;  
  mat_complx *m=NULL,*curr=NULL,*tmp=NULL;

  for (k=1;k <= S->nspins;k++) {

    m = _Iq(ss_qn(S,k));
    N = m->row;
    for (i=0;i<N;i++) {    
      for (j=0;j<N;j++) {
        double ok = 0;
	    int pp = i+j*N;
        if ( (m->data[pp].im - m->data[pp].re) == coh[k] ) ok=1;
        m->data[pp]=Complx(ok,0.0);
        //DEBUGPRINT("Icoherence: k=%d: element [%d, %d] = %f\n",k,i,j,ok);
      }
    }
    if (k == 1) {
      curr=m; 
    } else {
      tmp=cm_direct(curr,m);
      free_complx_matrix(curr);
      free_complx_matrix(m);
      curr=tmp;      
    }
  }
  if (sim->sparse) cm_sparse(curr, SPARSE_TOL);
  return curr;
}

/* this returns full matrix */
mat_complx * Iq(SpinSys *S)
{
  int i;
  mat_complx *curr=NULL,*tmp,*tmp2;

  for (i=1;i <= S->nspins;i++) {
    if (i == 1) {
      curr=_Iq(ss_qn(S,i)); 
    } else {
      tmp2=_Iq(ss_qn(S,i));
      tmp=cm_directadd(curr,tmp2);
      free_complx_matrix(curr);
      free_complx_matrix(tmp2);
      curr=tmp;      
    }
  }
  return curr;
}



mat_complx* Iqdelta(Sim_info* sim)
{
  SpinSys* S = sim->ss;
  int i,j,N;
  mat_complx * M;
  
  N=S->matdim;

  M=Iq(S);

  for (i=0;i<N;i++) {
    for (j=0;j<N;j++) {
      int pp = i+j*N;
      M->data[pp] = Complx( M->data[pp].im - M->data[pp].re, 0.0);
    }
  }  
  if (sim->sparse) cm_sparse(M, SPARSE_TOL);
  return M;
}

mat_complx * Ie(Sim_info* sim)
{
  mat_complx *a;

  if (sim->sparse)
	  a = complx_matrix(sim->matdim,sim->matdim,MAT_SPARSE_DIAG,sim->matdim,0);
  else
	  a = complx_matrix(sim->matdim,sim->matdim,MAT_DENSE_DIAG,1,0);
  cm_unit(a);
  return a;
}

mat_complx * _Ia(double I)
{
  mat_complx * M;

  if ( fabs(I-0.5)>0.01 ) {
     fprintf(stderr,"error in alpha-state operator, it is undefined for spin=%f\n",I);
     exit(1);
  }
  M = complx_matrix(2,2,MAT_DENSE,0,0);
  cm_zero(M);
  M->data[0].re = 1.0;
  return M;
}

mat_complx * _Ib(double I)
{
  mat_complx * M;

  if ( fabs(I-0.5)>0.01 ) {
     fprintf(stderr,"error in beta-state operator, it is undefined for spin=%f\n",I);
     exit(1);
  }
  M = complx_matrix(2,2,MAT_DENSE,0,0);
  cm_zero(M);
  M->data[3].re = 1.0;
  return M;
}

/****
 * ZT: adding I_alpha operator for spin=1/2
 ****/
mat_complx * Ia(Sim_info* sim, int spin) {
  SpinSys* S = sim->ss;
  mat_complx *a, *b;

  a = _Ia( ss_qn(S,spin) );
  b = In(sim,a,spin);
  free_complx_matrix(a);
  return b;
}

/****
 * ZT: adding I_beta operator for spin=1/2
 ****/
mat_complx * Ib(Sim_info *sim, int spin) {
  mat_complx *a, *b;

  a = _Ib(ss_qn(sim->ss,spin));
  b = In(sim,a,spin);
  free_complx_matrix(a);
  return b;
}

#define DECLARE_SIFUNC(TYPE) \
mat_complx * TYPE(Sim_info* sim,int spin)\
{\
  mat_complx *a,*b;\
  a= _##TYPE(ss_qn(sim->ss,spin));\
  b=In(sim,a,spin);\
  free_complx_matrix(a);\
  return b;\
}

DECLARE_SIFUNC(Ic)
DECLARE_SIFUNC(Ip)
DECLARE_SIFUNC(Im)
DECLARE_SIFUNC(Ix)
DECLARE_SIFUNC(Iy)

/* Iz is diagonal matrix */
mat_complx * Iz(Sim_info* sim, int nuc)
{
	mat_complx *res;
	int i, ii, iii, dim, K, L, N, cnt;
	double I;

	I = ss_qn(sim->ss,nuc);
	N = (int)(2*I+1);
	dim = sim->matdim;
	K = 1;
	for (i=1;i<nuc;i++) {
		K *= (int)(2.0*ss_qn(sim->ss,i)+1.0);
	}
	L = 1;
	for (i=nuc+1; i<=sim->ss->nspins; i++) {
		L *= (int)(2.0*ss_qn(sim->ss,i)+1.0);
	}
	assert(K*N*L == dim);

	res = complx_matrix(dim,dim,MAT_DENSE_DIAG,0,0);
	cnt = 0;
	for (i=0; i<K; i++) {
		for (ii=0; ii<N; ii++) {
			for (iii=0; iii<L; iii++) {
				res->data[cnt] = Complx(I - ii,0.0);
				cnt++;
			}
		}
	}
	return res;
}

/* moved here from readsys */

mat_double * Iz_ham(Sim_info* sim, int nuc)
{
	mat_double *res;
	int i, ii, iii, dim, K, L, N, cnt;
	double I;

	I = ss_qn(sim->ss,nuc);
	N = (int)(2*I+1);
	dim = sim->matdim;
	K = 1;
	for (i=1;i<nuc;i++) {
		K *= (int)(2.0*ss_qn(sim->ss,i)+1.0);
	}
	L = 1;
	for (i=nuc+1; i<=sim->ss->nspins; i++) {
		L *= (int)(2.0*ss_qn(sim->ss,i)+1.0);
	}
	assert(K*N*L == dim);

	res = double_matrix(dim,dim,MAT_DENSE_DIAG,0,sim->basis);
	cnt = 0;
	for (i=0; i<K; i++) {
		for (ii=0; ii<N; ii++) {
			for (iii=0; iii<L; iii++) {
				res->data[cnt] = I - ii;
				cnt++;
			}
		}
	}
	if (sim->basis != 0) {
		double* ddum = (double*)malloc(dim*sizeof(double));
		int *permvec = sim->perm_table[0+sim->basis*sim->Nbasis];
		for (i=0; i<dim; i++) ddum[i] = res->data[permvec[i+1]-1];
		free(res->data);
		res->data = ddum;
	}

	return res;
}

blk_mat_double * IzIz_sqrt2by3(Sim_info *sim, int n1, int n2)
{
	blk_mat_double *res;
	mat_double *Iz1, *Iz2;
	int i, j;
	double *d1, *d2, *d3;

	Iz1 = Iz_ham(sim,n1);
	Iz2 = Iz_ham(sim,n2);
	if (sim->basis == 0)
		res = create_blk_mat_double(sim->matdim,1,NULL,MAT_DENSE_DIAG,sim->basis);
	else
		res = create_blk_mat_double(sim->matdim,LEN(sim->dims_table[sim->basis]),sim->dims_table[sim->basis],MAT_DENSE_DIAG,sim->basis);
	d1 = Iz1->data;
	d2 = Iz2->data;
	for (i=0; i<res->Nblocks; i++) {
		d3 = res->m[i].data;
		assert(res->m[i].type == MAT_DENSE_DIAG);
		for (j=0; j<res->blk_dims[i]; j++) {
			*d3 = (*d1)*(*d2)*SQRT2BY3;
			d1++;
			d2++;
			d3++;
		}
	}
	free_double_matrix(Iz1);
	free_double_matrix(Iz2);
	return res;
}

mat_double * Ip_real(Sim_info *sim, int nuc)
{
	mat_double *res, *pom;
	int i,j,S;
	double I, mm,m;

	I = ss_qn(sim->ss,nuc);
	S= (int)(2*I+1);
	pom = double_matrix(S,S,MAT_DENSE,0,0);
	for (i=0;i<S;i++) {
		mm = I-i;
		for (j=0;j<S;j++) {
			m = I-j;
			if (mm == m + 1)
				pom->data[i+j*S] = sqrt(I*(I+1.0)-m*(m+1.0));
			else
				pom->data[i+j*S] = 0.0;
		}
	}
	res = In_real(sim,pom,nuc);
	free_double_matrix(pom);
	return res;
}

mat_double * Im_real(Sim_info *sim, int nuc)
{
	mat_double *res, *pom;
	int i,j,S;
	double I, mm,m;

	I = ss_qn(sim->ss,nuc);
	S= (int)(2*I+1);
	pom = double_matrix(S,S,MAT_DENSE,0,0);
	for (i=0;i<S;i++) {
		mm = I-i;
		for (j=0;j<S;j++) {
			m = I-j;
			if (mm == m - 1)
				pom->data[i+j*S] = sqrt(I*(I+1.0)-m*(m-1.0));
			else
				pom->data[i+j*S] = 0.0;
		}
	}
	res = In_real(sim,pom,nuc);
	free_double_matrix(pom);
	return res;
}

mat_double * Ix_real(Sim_info *sim, int nuc)
{
	mat_double *res, *pom;
	int i,j,S;
	double I, mm,m;

	I = ss_qn(sim->ss,nuc);
	S= (int)(2*I+1);
	pom = double_matrix(S,S,MAT_DENSE,0,0);
	dm_zero(pom);
	for (i=0;i<S;i++) {
		mm = I-i;
		for (j=0;j<S;j++) {
			m = I-j;
			/* this is for Im */
			if (mm == m - 1) pom->data[i+j*S] += 0.5*sqrt(I*(I+1.0)-m*(m-1.0));
			/* this is for Ip */
			if (mm == m + 1) pom->data[i+j*S] += 0.5*sqrt(I*(I+1.0)-m*(m+1.0));
		}
	}
	res = In_real(sim,pom,nuc);
	free_double_matrix(pom);
	return res;
}

mat_double * Iy_real(Sim_info *sim, int nuc)
{
	mat_double *res, *pom;
	int i,j,S;
	double I, mm,m;

	I = ss_qn(sim->ss,nuc);
	S= (int)(2*I+1);
	pom = double_matrix(S,S,MAT_DENSE,0,0);
	dm_zero(pom);
	for (i=0;i<S;i++) {
		mm = I-i;
		for (j=0;j<S;j++) {
			m = I-j;
			/* this is for Im */
			if (mm == m - 1) pom->data[i+j*S] += 0.5*sqrt(I*(I+1.0)-m*(m-1.0));
			/* this is for Ip */
			if (mm == m + 1) pom->data[i+j*S] -= 0.5*sqrt(I*(I+1.0)-m*(m+1.0));
		}
	}
	res = In_real(sim,pom,nuc);
	free_double_matrix(pom);
	return res;
}

/* ( I+ S-  +  I- S+ )/2  =  IxSx + IySy   */
blk_mat_double * zq_ham(Sim_info *sim, int n1, int n2)
{
	blk_mat_double *res;
	mat_double *a, *b;
	int nb,dim,i,j;

	assert(n1 != n2); // this code hold for this ONLY!!!
	a = Ip_real(sim,n1);
	b = Im_real(sim,n2);
	dm_multo(a,b);
	free_double_matrix(b);

	if (sim->basis == 0) {
		res = create_blk_mat_double(sim->matdim,1,NULL,sim->sparse ? MAT_SPARSE : MAT_DENSE,sim->basis);
		dm_zero(res->m);
		dm_multod(res->m,a,0.5);
		dm_transposei(a);
		dm_multod(res->m,a,0.5);
		free_double_matrix(a);
		return res;
	}

	// when blocking is active:
	//printf("a is %s\n",matrix_type(a->type));
	res = create_blk_mat_double(sim->matdim,LEN(sim->dims_table[sim->basis]),sim->dims_table[sim->basis],sim->sparse ? MAT_SPARSE : MAT_DENSE,sim->basis);
	int sft = 1;
	int *permvec = sim->perm_table[0+sim->basis*sim->Nbasis];
	for (nb=0; nb<res->Nblocks; nb++) {
		dim = res->blk_dims[nb];
		//printf("nb = %d, dim = %d, %s\n",nb,dim,matrix_type(res->m[nb].type));
		if (dim == 1) {
			res->m[nb].data[0] = dm_getelem(a,permvec[sft],permvec[sft]);
		} else if (res->m[nb].type == MAT_DENSE) {
			for (i=0; i<dim; i++) {
				res->m[nb].data[i+i*dim] = dm_getelem(a,permvec[i+sft],permvec[i+sft]);
				for (j=i+1; j<dim; j++) {
					res->m[nb].data[i+j*dim] = 0.5*dm_getelem(a,permvec[i+sft],permvec[j+sft]);
					res->m[nb].data[i+j*dim] += 0.5*dm_getelem(a,permvec[j+sft],permvec[i+sft]);
					res->m[nb].data[j+i*dim] = res->m[nb].data[i+j*dim];
				}
			}
		} else {
			assert(res->m[nb].type == MAT_SPARSE);
			b = dm_get_diagblock_permute(a,permvec,sft,dim);
			b->basis = sim->basis;
			dm_zero(&(res->m[nb]));
			dm_multod(&(res->m[nb]),b,0.5);
			dm_transposei(b);
			dm_multod(&(res->m[nb]),b,0.5);
			free_double_matrix(b);
		}
		sft += dim;
	}
	free_double_matrix(a);
	return res;
}

blk_mat_double * zq_ham_permuted(Sim_info *sim, int n1, int n2)
{
	blk_mat_double *res;
	mat_double *a, *b;
	int nb,dim,i,j;

	assert(n1 != n2); // this code hold for this ONLY!!!

	if (sim->basis == 0) {
		a = Ip_real(sim,n1);
		b = Im_real(sim,n2);
		dm_multo(a,b);
		free_double_matrix(b);
		res = create_blk_mat_double(sim->matdim,1,NULL,sim->sparse ? MAT_SPARSE : MAT_DENSE,sim->basis);
		dm_zero(res->m);
		dm_multod(res->m,a,0.5);
		dm_transposei(a);
		dm_multod(res->m,a,0.5);
		free_double_matrix(a);
	} else {
		// when blocking is active:
		a = Ip_real(sim,n1);
		b = Im_real(sim,n2);
		dm_multo(a,b);
		free_double_matrix(b);
		dm_permute(a,sim->perm_table[0+sim->basis*sim->Nbasis]);
		res = create_blk_mat_double(sim->matdim,LEN(sim->dims_table[sim->basis]),sim->dims_table[sim->basis],sim->sparse ? MAT_SPARSE : MAT_DENSE,sim->basis);
		int sft = 0;
		for (nb=0; nb<res->Nblocks; nb++) {
			dim = res->blk_dims[nb];
			if (dim == 1) {
				res->m[nb].data[0] = dm_getelem(a,sft+1,sft+1);
			} else {
				b = dm_extract_block(a,sft,sft,dim,dim);
				dm_zero(&(res->m[nb]));
				if (b != NULL) {
					b->basis = sim->basis;
					dm_multod(&(res->m[nb]),b,0.5);
					dm_transposei(b);
					dm_multod(&(res->m[nb]),b,0.5);
					free_double_matrix(b);
				}
			}
			sft += dim;
		}
		free_double_matrix(a);
	}
	return res;
}


blk_mat_double * II(Sim_info* sim, int n1, int n2)
{
	blk_mat_double *a, *b;

	a = IzIz_sqrt2by3(sim,n1,n2);
	b = zq_ham(sim,n1,n2);
	blk_dm_multod(b,a,1.0/SQRT2BY3);
	free_blk_mat_double(a);

	return b;
}

blk_mat_double * T20(Sim_info* sim, int n1, int n2)
{
	blk_mat_double *a, *b;

	a = IzIz_sqrt2by3(sim,n1,n2);
	//b = zq_ham(sim,n1,n2);
	b = zq_ham_permuted(sim,n1,n2);
	//blk_dm_print(b,"zq_ham");
	blk_dm_multod(a,b,-1.0/sqrt(6.0));
	free_blk_mat_double(b);

	return a;
}

mat_double * T20II(Sim_info *sim, int nuc)
{
	mat_double *res;
	int i, ii, iii, dim, K, L, N, cnt;
	double I, val;

	I = ss_qn(sim->ss,nuc);
	N = (int)(2*I+1);
	dim = sim->matdim;
	K = 1;
	for (i=1;i<nuc;i++) {
		K *= (int)(2.0*ss_qn(sim->ss,i)+1.0);
	}
	L = 1;
	for (i=nuc+1; i<=sim->ss->nspins; i++) {
		L *= (int)(2.0*ss_qn(sim->ss,i)+1.0);
	}
	assert(K*N*L == dim);

	res = double_matrix(dim,dim,MAT_DENSE_DIAG,0,sim->basis);
	cnt = 0;
	for (i=0; i<K; i++) {
		for (ii=0; ii<N; ii++) {
			val = (I - ii)*(I - ii)*SQRT2BY3;
			if (ii < N-1) val -= (I*(I+1)-(I-ii-1)*(I-ii))*(0.5/sqrt(6.0));
			if (ii>0) val -= (I*(I+1)-(I-ii+1)*(I-ii))*(0.5/sqrt(6.0));
			for (iii=0; iii<L; iii++) {
				res->data[cnt] = val;
				cnt++;
			}
		}
	}
	if (sim->basis != 0) {
		double* ddum = (double*)malloc(dim*sizeof(double));
		int *permvec = sim->perm_table[0+sim->basis*sim->Nbasis];
		for (i=0; i<dim; i++) ddum[i] = res->data[permvec[i+1]-1];
		free(res->data);
		res->data = ddum;
	}

	return res;
}

void fill_Tquad_2(Sim_info *sim, int nuc, Quadrupole *qptr)
{
	mat_double *Ta, *Tb;
	int i, ii, iii, dim, K, L, N, cnt;
	double I, val1, val2;

	I = ss_qn(sim->ss,nuc);
	N = (int)(2*I+1);
	dim = sim->matdim;
	K = 1;
	for (i=1;i<nuc;i++) {
		K *= (int)(2.0*ss_qn(sim->ss,i)+1.0);
	}
	L = 1;
	for (i=nuc+1; i<=sim->ss->nspins; i++) {
		L *= (int)(2.0*ss_qn(sim->ss,i)+1.0);
	}
	assert(K*N*L == dim);

	Ta = double_matrix(dim,dim,MAT_DENSE_DIAG,0,sim->basis);
	Tb = double_matrix(dim,dim,MAT_DENSE_DIAG,0,sim->basis);
	cnt = 0;
	for (i=0; i<K; i++) {
		for (ii=0; ii<N; ii++) {
			val1 = (I-ii)*( 2*I*(I+1) -2*(I-ii)*(I-ii) -1 );
			val2 = (I-ii)*( 4*I*(I+1) -8*(I-ii)*(I-ii) -1 );
			for (iii=0; iii<L; iii++) {
				Ta->data[cnt] = val1;
				Tb->data[cnt] = val2;
				cnt++;
			}
		}
	}

	if (sim->basis != 0) {
		double *pom, *ddum = (double*)malloc(dim*sizeof(double));
		int *permvec = sim->perm_table[0+sim->basis*sim->Nbasis];
		for (i=0; i<dim; i++) ddum[i] = Ta->data[permvec[i+1]-1];
		pom = Ta->data; Ta->data = ddum;
		for (i=0; i<dim; i++) pom[i] = Tb->data[permvec[i+1]-1];
		free(Tb->data);
		Tb->data = pom;
	}

	qptr->Ta = Ta;
	qptr->Tb = Tb;
}

void fill_Tquad_3(Sim_info *sim, Quadrupole *qptr)
{
	int i, ii, iii, K, L;
	double I = ss_qn(sim->ss,qptr->nuc);
	int dim = (int)(2*I+1);

	qptr->T3a = double_matrix(sim->matdim,sim->matdim,MAT_DENSE_DIAG,0,sim->basis);
	qptr->T3b = double_matrix(sim->matdim,sim->matdim,MAT_DENSE_DIAG,0,sim->basis);
	qptr->T3c = double_matrix(sim->matdim,sim->matdim,MAT_DENSE_DIAG,0,sim->basis);
	K = 1;
	for (i=1;i<qptr->nuc;i++) {
		K *= (int)(2.0*ss_qn(sim->ss,i)+1.0);
	}
	L = 1;
	for (i=(qptr->nuc)+1; i<=sim->ss->nspins; i++) {
		L *= (int)(2.0*ss_qn(sim->ss,i)+1.0);
	}
	assert(K*dim*L == sim->matdim);
	int cnt = 0;
	for (i=0; i<K; i++) {
		for (ii=0; ii<dim; ii++) {
			double m = I - ii;
			double mm = m*m;
			for (iii=0; iii<L; iii++) {
				qptr->T3a->data[cnt] = sqrt(6.0)/4.0*(20*mm*mm+7*mm-(12*mm+1)*I*(I+1));
				qptr->T3b->data[cnt] = sqrt(6.0)/4.0*(5*mm*mm+7*mm-(6*mm+2)*I*(I+1)+I*I*(I+1)*(I+1));
				qptr->T3c->data[cnt] = 0.25*(30*mm*mm+15*mm-(24*mm+3)*I*(I+1)+2*I*I*(I+1)*(I+1));
				cnt++;
			}
		}
	}

	if (sim->basis != 0) {
		double* ddum1 = (double*)malloc(sim->matdim*sizeof(double));
		double* ddum2 = (double*)malloc(sim->matdim*sizeof(double));
		double* ddum3 = (double*)malloc(sim->matdim*sizeof(double));
		int *permvec = sim->perm_table[0+sim->basis*sim->Nbasis];
		for (i=0; i<sim->matdim; i++) {
			ddum1[i] = qptr->T3a->data[permvec[i+1]-1];
			ddum2[i] = qptr->T3b->data[permvec[i+1]-1];
			ddum3[i] = qptr->T3c->data[permvec[i+1]-1];
		}
		free(qptr->T3a->data); qptr->T3a->data = ddum1;
		free(qptr->T3b->data); qptr->T3a->data = ddum2;
		free(qptr->T3c->data); qptr->T3a->data = ddum3;
	}

}

void fill_Tmix_dipole(Sim_info *sim, Mixing *mptr)
{
	/* operator T, common for both hetero and homonuclear cases */
	mat_double *mx1, *mx2;
	int i;

	mx1 = T20II(sim,mptr->couple[0]);
	mx2 = Iz_ham(sim,mptr->couple[1]);
	assert(mx1->type == MAT_DENSE_DIAG);
	assert(mx2->type == MAT_DENSE_DIAG);
	for (i=0; i<mx1->row; i++) {
		mx1->data[i] *= mx2->data[i]*sqrt(6.0);
	}
	free_double_matrix(mx2);
	mptr->T = mx1;
}

void fill_Tabmix_dipole(Sim_info *sim, Mixing *mptr)
{
	mat_double *a, *b, *c;

	a = Ip_real(sim,mptr->couple[0]);
	b = Im_real(sim,mptr->couple[1]);
	dm_multo(a,b);
	free_double_matrix(b);
	b = Im_real(sim,mptr->couple[0]);
	c = Ip_real(sim,mptr->couple[1]);
	dm_multo(b,c);
	free_double_matrix(c);
	a->basis = b->basis = 0;

	int basis = sim->basis;
	sim->basis = 0;
	c = Iz_ham(sim, mptr->couple[0]);
	sim->basis = basis;
	dm_muld(c,2.0);
	dm_addtodiag(c,1.0);
	dm_multo(a,c);
	dm_addtodiag(c,-2.0);
	dm_multo(b,c);
	free_double_matrix(c);

	if (sim->basis == 0) {
		mptr->Ta = create_blk_mat_double(sim->matdim,1,NULL,sim->sparse ? MAT_SPARSE : MAT_DENSE,sim->basis);
		dm_copy(mptr->Ta->m,a);
		mptr->Tb = create_blk_mat_double(sim->matdim,1,NULL,sim->sparse ? MAT_SPARSE : MAT_DENSE,sim->basis);
		dm_copy(mptr->Tb->m,b);
	} else {
		// when blocking is active:
		mptr->Ta = create_blk_mat_double(sim->matdim,LEN(sim->dims_table[sim->basis]),sim->dims_table[sim->basis],sim->sparse ? MAT_SPARSE : MAT_DENSE,sim->basis);
		mptr->Tb = create_blk_mat_double(sim->matdim,LEN(sim->dims_table[sim->basis]),sim->dims_table[sim->basis],sim->sparse ? MAT_SPARSE : MAT_DENSE,sim->basis);
		int nb, sft = 1;
		int *permvec = sim->perm_table[0+sim->basis*sim->Nbasis];
		for (nb=0; nb<mptr->Ta->Nblocks; nb++) {
			int dim = mptr->Ta->blk_dims[nb];
			if (dim == 1) {
				mptr->Ta->m[nb].data[0] = dm_getelem(a,permvec[sft],permvec[sft]);
				mptr->Tb->m[nb].data[0] = dm_getelem(b,permvec[sft],permvec[sft]);
			} else {
				c = dm_get_diagblock_permute(a,permvec,sft,dim);
				c->basis = sim->basis;
				dm_copy(&(mptr->Ta->m[nb]),c);
				free_double_matrix(c);
				c = dm_get_diagblock_permute(b,permvec,sft,dim);
				c->basis = sim->basis;
				dm_copy(&(mptr->Tb->m[nb]),c);
				free_double_matrix(c);
			}
			sft += dim;
		}
	}
	free_double_matrix(a);
	free_double_matrix(b);
}
/* end of readsys stuff */

#ifdef DEBUG
#define DEBUG_PARSER(x) printf(x)
#define DEBUG_PARSER2(x,y) printf(x,y)
#else
#define DEBUG_PARSER(x)
#define DEBUG_PARSER2(x,y)
#endif

mat_complx * ss_oper(Sim_info* sim,char* name)
{
  SpinSys* S = sim->ss;
  int i,spin,n,N;
  char buf[16],*pname,*pbuf;
  mat_complx *sum,*tmp;
  
  pname=name;
  if (*pname != 'I') {
     fprintf(stderr,"operator `%s` not known, must be of type "
             "Ina, where n=1,2,.., or n='n' for sum of all spins and "      
             "a=x,y,z,p,m\n",name);
     exit(1);
  }
  pname++;

  if (*pname == 'n') {
    pname++;

    n=S->nspins;
    N=S->matdim;
    if (*(pname+1) == 'c') {
      /* The syntax I2xc is not properply implemented yet. */
      fprintf(stderr,"oper: unknown operator: %s\n",name);
      exit(1);
      pname++;
    } else {
      switch (*pname) {
      case 'x':
    	  sum = complx_matrix(N,N,(sim->sparse) ? MAT_SPARSE : MAT_DENSE,1,0);
    	  cm_zero(sum);
    	  for (i=1;i<=n;i++) {
    		  tmp=Ix(sim,i);
    		  cm_addto(sum,tmp);
    		  free_complx_matrix(tmp);
    	  }
    	  break;
      case 'y':
    	  sum = complx_matrix(N,N,(sim->sparse) ? MAT_SPARSE : MAT_DENSE,1,0);
    	  cm_zero(sum);
    	  for (i=1;i<=n;i++) {
    		  tmp=Iy(sim,i);
    		  cm_addto(sum,tmp);
    		  free_complx_matrix(tmp);
    	  }
    	  break;
      case 'z':
    	  sum = complx_matrix(N,N,(sim->sparse) ? MAT_SPARSE_DIAG : MAT_DENSE_DIAG,1,0);
    	  cm_zero(sum);
    	  for (i=1;i<=n;i++) {
    		  tmp=Iz(sim,i);
    		  cm_addto(sum,tmp);
    		  free_complx_matrix(tmp);
    	  }
    	  break;
      case 'p':
    	  sum = complx_matrix(N,N,(sim->sparse) ? MAT_SPARSE : MAT_DENSE,1,0);
    	  cm_zero(sum);
    	  for (i=1;i<=n;i++) {
    		  tmp=Ip(sim,i);
    		  cm_addto(sum,tmp);
    		  free_complx_matrix(tmp);
    	  }
    	  break;
      case 'm':
    	  sum = complx_matrix(N,N,(sim->sparse) ? MAT_SPARSE : MAT_DENSE,1,0);
    	  cm_zero(sum);
    	  for (i=1;i<=n;i++) {
    		  tmp=Im(sim,i);
    		  cm_addto(sum,tmp);
    		  free_complx_matrix(tmp);
    	  }
    	  break;
      case 'c':
    	  sum = complx_matrix(N,N,(sim->sparse) ? MAT_SPARSE : MAT_DENSE,1,0);
    	  cm_zero(sum);
    	  for (i=1;i<=n;i++) {
    		  tmp=Ic(sim,i);
    		  cm_addto(sum,tmp);
    		  free_complx_matrix(tmp);
    	  }
    	  break;
      case 'a':
    	  sum = complx_matrix(N,N,(sim->sparse) ? MAT_SPARSE : MAT_DENSE,1,0);
    	  cm_zero(sum);
    	  for (i=1;i<=n;i++) {
    		  tmp=Ia(sim,i);
    		  cm_addto(sum,tmp);
    		  free_complx_matrix(tmp);
    	  }
    	  break;
      case 'b':
    	  sum = complx_matrix(N,N,(sim->sparse) ? MAT_SPARSE : MAT_DENSE,1,0);
    	  cm_zero(sum);
    	  for (i=1;i<=n;i++) {
    		  tmp=Ib(sim,i);
    		  cm_addto(sum,tmp);
    		  free_complx_matrix(tmp);
    	  }
    	  break;
       default:
         fprintf(stderr,"oper: unknown operator: '%s'\n",name);
         exit(1);
      }
    }
    pname++;
    if (*pname != 0) {
       fprintf(stderr,"oper: trailing garbage at operator: '%s'\n",name);
       exit(1);
    }
    return sum;

  } else {
    if (!isdigit(*pname) || *pname == '0') {
       fprintf(stderr,"oper: unknown operator: '%s'\n",name);
       exit(1);
    }
    pbuf=buf;
    *pbuf++ = *pname++;
    while (isdigit(*pname)) {
      *pbuf++ = *pname++;
      if (pname-name > 6) {
         fprintf(stderr,"oper: unknown operator: '%s'\n",name);
         exit(1);
      }
    }
    *pbuf=0;
    spin=atoi(buf);
    
    if (*(pname+1) == 'c') {
      /* The syntax I2xc is not properply implemented yet. */
      fprintf(stderr,"oper: unknown operator: %s\n",name);
      exit(1);
      pname++;
    } else {
      switch (*pname) {
      case 'x': return Ix(sim,spin); break;
      case 'y': return Iy(sim,spin); break;
      case 'z': return Iz(sim,spin); break;
      case 'p': return Ip(sim,spin); break;
      case 'm': return Im(sim,spin); break;
      case 'c': return Ic(sim,spin); break;
      case 'a': return Ia(sim,spin); break;
      case 'b': return Ib(sim,spin); break;
      default:
         fprintf(stderr,"oper: unknown operator: %s\n",name);
         exit(1);
      }
    }
    if (*(pname+1) != 0) {
       fprintf(stderr,"oper: trailing garbage at operator: '%s'\n",name);
       exit(1);
    }

  }
  return NULL;
}

/*
 Description:  A parser with the following grammar

 Grammar:
 ----------------------------
 program:
	expression RETURN
	
 expression:
	expression + term
	expression - term
	term
	
 term:
	term * primary
	primary
	
 primary:
	SPINOPERATOR
	NUMBER
	IMAG
    -primary
	( expression )	
 ----------------------------
*/ 
/* Defines of tokens */
#define	sop_token_SPINOPERATOR 1
#define	sop_token_NUMBER 2
#define	sop_token_IMAG 3
#define	sop_token_LP 4
#define	sop_token_RP 5

#define	sop_token_MINUS 6
#define	sop_token_MULT 7
#define	sop_token_PLUS 8
#define	sop_token_RETURN 9

typedef struct _sop_struct {
	char		initial[256];
	complx		complexval;
	mat_complx*	operator;
	int			current_Token;
	int			matrix_dimension;
	int			nspins;
	char*		current_chrptr;
	Sim_info*	siminfoptr;
} sop_struct;

int sop_get_token(sop_struct *sop)
{
	char *chr,*chr0;
	char buf[100],*bufp;
	
	chr = sop->current_chrptr;
	chr0 = chr;

	while (*chr==' ') {
		chr++;
	}

	switch (*chr)
	{
	case '\0':
		return sop_token_RETURN;
	case '+':
		DEBUG_PARSER("PLUS\n");
		sop->current_chrptr = ++chr;
		return sop_token_PLUS;
	case '*':
		DEBUG_PARSER("MULT\n");
		sop->current_chrptr = ++chr;
		return sop_token_MULT;
	case '-':
		DEBUG_PARSER("MINUS\n");
		sop->current_chrptr = ++chr;
		return sop_token_MINUS;
	case 'i':
		DEBUG_PARSER("IMAG\n");
		sop->current_chrptr = ++chr;
		return sop_token_IMAG;
	case 'I':
		bufp=buf;
		*bufp++ = *chr++;
		if (*chr=='n') {
			*bufp++ = *chr++;
		} else if (isdigit(*chr) && *chr != '0') {
			*bufp++ = *chr++;
			while (isdigit(*chr)) {
				*bufp++ = *chr++;
			}
		} else {
			fprintf(stderr,"Cannot parse spinoperator at character"
					" '%s' in operator '%s'\n",chr,sop->initial);
			exit(-1);
		}
		switch (*chr) {
		case 'x': *bufp ='x';break;
		case 'y': *bufp ='y';break;
		case 'z': *bufp ='z';break;
		case 'p': *bufp ='p';break;
		case 'm': *bufp ='m';break;
		case 'c': *bufp ='c';break;
		case 'a': *bufp ='a';break;
		case 'b': *bufp ='b';break;
		default:
			fprintf(stderr,"Cannot parse spinoperator at character"
					" '%s' in operator '%s'\n",chr,sop->initial);
			exit(-1);
		}
		bufp++;
		if (*(chr+1) == 'c') {
			*bufp++ = *(++chr);
		}
		*bufp=0;
				
		DEBUG_PARSER2("SPINOPERATOR: %s\n",buf);
		sop->current_chrptr = ++chr;
		sop->operator = ss_oper(sop->siminfoptr,buf);
		return sop_token_SPINOPERATOR;
	case '(':
		DEBUG_PARSER2("LP : %c\n",*chr);
		sop->current_chrptr = ++chr;
		return sop_token_LP;
	case ')':
		DEBUG_PARSER2("RP : %c\n",*chr);
		sop->current_chrptr = ++chr;
		return sop_token_RP;
	case '1': case '2': case '3': case '4': case '5':
	case '6': case '7': case '8': case '9': case '0': case '.':
		bufp = buf;
		if (*chr != '.') {
			*bufp++ = *chr++;
			while (isdigit(*chr)) *bufp++ = *chr++;
		}
		if (*chr == '.') {
			*bufp++ = *chr++;
			if (!isdigit(*chr)) {
				fprintf(stderr,"Cannot parse spinoperator at character"
						" '%s' in operator '%s'\n",chr,sop->initial);
				exit(-1);
			}
			while (isdigit(*chr)) *bufp++ = *chr++;
		}
		if (*chr == 'E' || *chr == 'e') {
			*bufp++ = *chr++;
			if (*chr == '+' || *chr == '-') *bufp++ = *chr++;
			if (!isdigit(*chr)) {
				fprintf(stderr,"Cannot parse spinoperator at character"
						" '%s' in operator '%s'\n",chr,sop->initial);
				exit(-1);
			}
			while (isdigit(*chr)) *bufp++ = *chr++;
		}
		if (!isdigit(*chr)) chr--;
		*bufp=0;
		DEBUG_PARSER2("NUMBER : %s\n",buf);
		bufp=buf;
		sop->complexval.re = strtod(buf,&bufp);
		sop->complexval.im = 0.0;
		sop->current_chrptr = ++chr;
		return sop_token_NUMBER;
	default:
		fprintf(stderr,"Cannot parse spinoperator at character"
				" '%s' in operator '%s'\n",chr,sop->initial);
		exit(-1);
	}
}

mat_complx * sop_term(sop_struct *sop, int get);

mat_complx * sop_expression(sop_struct *sop, int get)
{
	mat_complx * M;
	mat_complx * left;

	left = sop_term(sop, get);

	for (;;)
	{
		switch (sop->current_Token)
		{
		case sop_token_PLUS:
			M = sop_term(sop,1);
			cm_addto(left,M);
			free_complx_matrix(M);
			break;
		case sop_token_MINUS:
			M = sop_term(sop,1);
			cm_subfrom(left,M);
			free_complx_matrix(M);
			break;
		default:
			return left;
		}
	}
}

mat_complx * sop_prim(sop_struct *sop, int get);

mat_complx * sop_term(sop_struct *sop, int get)
{
	mat_complx* M;
	mat_complx* M2;
	mat_complx* left;
	
	left = sop_prim(sop,get);

	for (;;)
	{
		switch (sop->current_Token)
		{
		case sop_token_MULT:
			M2 = sop_prim(sop,1);
			M = cm_mul(left,M2);
			free_complx_matrix(M2);
			free_complx_matrix(left);
			left = M;
			break;
		default:
			return left;
		}
	}
}

mat_complx * sop_prim(sop_struct *sop, int get)
{
	mat_complx * M;

	if (get) 
	{
		sop->current_Token = sop_get_token(sop);
	}

	switch (sop->current_Token)
	{
	case sop_token_SPINOPERATOR:
		M = cm_dup(sop->operator);
		free_complx_matrix(sop->operator);
		sop->current_Token = sop_get_token(sop);
		if(sop->current_Token < sop_token_RP)
		{
			fprintf(stderr,"Parse error: primary not expected in operator '%s'\n",sop->initial);
			exit(1);
		}
		return M;

	case sop_token_NUMBER:
		M = cm_creatediag('d',sop->matrix_dimension,sop->complexval,0); //creates diagonal matrix filled with complexval
		sop->current_Token = sop_get_token(sop);
		if(sop->current_Token < sop_token_RP)
		{
			fprintf(stderr,"Parse error: primary not expected in operator '%s'\n",sop->initial);
			exit(1);
		}
		return M;

	case sop_token_IMAG:
		M = cm_creatediag('d',sop->matrix_dimension,Complx(0.0,1.0),0);
		sop->current_Token = sop_get_token(sop);
		if(sop->current_Token < sop_token_RP)
		{
			fprintf(stderr,"Parse error: primary not expected in operator '%s'\n",sop->initial);
			exit(1);
		}
		return M;

	case sop_token_RP:
	  fprintf(stderr,"Parse error: \")\" not expected in operator '%s'\n",sop->initial);
	  exit(1);
    return 0;
    
	case sop_token_MINUS:
		M = sop_prim(sop,1);
		cm_muld(M,-1.0);
		return M;

	case sop_token_LP:
		M = sop_expression(sop,1);
		if (sop->current_Token != sop_token_RP)
		{
		  fprintf(stderr,"Parse error: expected \")\" in operator '%s'\n",sop->initial);
			exit(1);
		}
		sop->current_Token = sop_get_token(sop);
		if(sop->current_Token < 6)
		{
			fprintf(stderr,"Parse error: primary not expected in operator '%s'\n",sop->initial);
			exit(1);
		}
		return M;
	default:
	  fprintf(stderr,"Parse error: primary expected in operator '%s'\n",sop->initial);
		exit(1);
	}
	return 0;
}

mat_complx* ss_readoper(Sim_info * sim,char* psop)
{  
	sop_struct *sop;
	mat_complx *res;

	sop = (sop_struct*)malloc(sizeof(sop_struct));
	sop->siminfoptr = sim;
	sop->matrix_dimension = sim->matdim;
    sop->nspins = sim->ss->nspins;
	sop->current_chrptr = psop;
	strcpy(sop->initial,sop->current_chrptr);
	sop->current_Token = sop_get_token(sop);
	if (sop->current_Token==sop_token_RETURN) {
		fprintf(stderr,"Parse error: operator is empty\n");
		exit(1);
	}
	res = sop_expression(sop,0);
	if (sop->current_Token != sop_token_RETURN)
	{
		fprintf(stderr,"Parse error: invalid syntax in '%s'\n",sop->initial);
		exit(1);
	}
	free(sop);

	if (sim->basis != 0) {
		cm_permute(res, sim->perm_table[0+sim->basis*sim->Nbasis]);
	}
	res->basis = sim->basis;
	return res;
}

/* 
END OF PARSERCODE
  */
  
  
/******
 * ZT: this creates Boltzman equilibrium density matrix for a given spin system
 *     High temperature limit, normalized according to 1H polarization
 *****/
mat_complx * ss_eq_sigma(Sim_info* sim)
{
   int n, i;
   mat_complx *mx, *sum;
   double b;
   SpinSys *S=sim->ss;
   
   n = S->nspins;
   if (sim->sparse)
	   sum = complx_matrix(S->matdim,S->matdim,MAT_SPARSE_DIAG,1,0);
   else
	   sum = complx_matrix(S->matdim,S->matdim,MAT_DENSE_DIAG,1,0);
   cm_zero(sum);
   for (i=1; i<=n; i++) {
      mx = Iz(sim,i);
      b = ss_gamma(S,i)/ss_gamma1H();
      cm_multod(sum,mx,b);
      free_complx_matrix(mx);
   }

   if (sim->basis != 0) {
		cm_permute(sum, sim->perm_table[0+sim->basis*sim->Nbasis]);
   }
   sum->basis = sim->basis;

   return sum;
}

int tclIsotopes(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  char buf[256];
  ISOTOPE* ptr;
  if (argc != 1) {
    //interp->result = "Usage: isotopes  (returns a list of isotopes: number, nuclei, spin, and gyromag ratio)";
    //return TCL_ERROR;
    return TclError(interp,"%s","Usage: isotopes  (returns a list of isotopes: number, nuclei, spin, and gyromag ratio)");
  }

  Tcl_ResetResult(interp);
  ptr=isotopes;
  while (ptr->number != 0) {
     sprintf(buf,"%d %s %g %g",ptr->number,ptr->name,ptr->spin,ptr->gamma);
     Tcl_AppendElement(interp,buf);
     ptr++;
  }
  return TCL_OK;
}    

#define HBAR 1054.59198

int tclDist2Dip(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  char* nuc1,*nuc2;
  double dist,dip,gamma1,gamma2;
  ISOTOPE* isop;

  if (argc != 4) {
    //interp->result = "dist2dip <nuc1> <nuc2> <distance> : calculates dipolar coupling (Hz) from internuclear distance (Aangstrom)";
    //return TCL_ERROR;
    return TclError(interp,"%s","dist2dip <nuc1> <nuc2> <distance> : calculates dipolar coupling (Hz) from internuclear distance (Aangstrom)");
  }
  nuc1=Tcl_GetString(argv[1]);
  nuc2=Tcl_GetString(argv[2]);
  if (Tcl_GetDoubleFromObj(interp,argv[3],&dist) != TCL_OK)
    return TCL_ERROR;   

  isop=ss_findisotope(nuc1);
  gamma1=isop->gamma;
  isop=ss_findisotope(nuc2);
  gamma2=isop->gamma;
 
  dip= -gamma1*gamma2*HBAR/(2.0*M_PI*dist*dist*dist);

  return TclSetResult(interp,"%g",dip);
}

int tclDip2Dist(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  char* nuc1,*nuc2;
  double dist,dip,gamma1,gamma2;
  ISOTOPE* isop;

  if (argc != 4)
    return TclError(interp,"dip2dist <nuc1> <nuc2> <dipolar coupling> :"
    " calculates internuclear distance (Aangstrom) from dipolar coupling (Hz)");

  nuc1=Tcl_GetString(argv[1]);
  nuc2=Tcl_GetString(argv[2]);
  if (Tcl_GetDoubleFromObj(interp,argv[3],&dip) != TCL_OK)
    return TCL_ERROR;   

  isop=ss_findisotope(nuc1);
  gamma1=isop->gamma;
  isop=ss_findisotope(nuc2);
  gamma2=isop->gamma;
  if (-dip/gamma1*gamma2 < 0.0)
     return TclError(interp,"dip2dist: the dipolar coupling for nuclei %s and %s must be %s",nuc1,nuc2,
           (gamma1*gamma2 < 0.0 ? "positive" : "negative"));

  dist=pow(fabs(gamma1*gamma2*HBAR/(M_PI*2.0*dip)),1.0/3.0);  
  return TclSetResult(interp,"%g",dist);
}

int tclGamma(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  ISOTOPE* isop;

  if (argc != 2)
    return TclError(interp,
    "gamma <nuc> : returns the magnetogyric ratio of the nucleus in the unit 10^7 rad/(T s)");
    
  isop=ss_findisotope(argv[1]);
  //sprintf(interp->result,"%g",isop->gamma);
  //return TCL_OK;
  return TclSetResult(interp,"%g",isop->gamma);
}


int tclResfreq(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  double protonfield;
  ISOTOPE* isop;

  if (argc != 2 && argc != 3)
    return TclError(interp,
     "resfreq <nuc> ?<proton_frequency/Hz default= 1e6Hz>?: returns the absolute resonance\n"
     " frequency in Hz for a nucleus optionally at a given proton frequency");
    
  isop=ss_findisotope(Tcl_GetString(argv[1]));
  protonfield=1.0;
  if (argc == 3) {
    if (Tcl_GetDoubleFromObj(interp,argv[2],&protonfield) != TCL_OK)
      return TCL_ERROR;   
    if (protonfield < 10000.0)
      return TclError(interp,"resfreq: illegal value of proton frequency (must be positive and > 10000Hz)\n");
  }
  return TclSetResult(interp, "%g", fabs(isop->gamma/ss_gamma1H()*protonfield));
}

void tclcmd_spinsys(Tcl_Interp* interp)
{

  Tcl_CreateObjCommand(interp,"isotopes",(Tcl_ObjCmdProc *)tclIsotopes,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"gamma",(Tcl_CmdProc *)tclGamma,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"resfreq",(Tcl_ObjCmdProc *)tclResfreq,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"dist2dip",(Tcl_ObjCmdProc *)tclDist2Dip,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"dip2dist",(Tcl_ObjCmdProc *)tclDip2Dist,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

}

              



