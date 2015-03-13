/*
 * auxmath.c
 *
 *  Created on: 9.1.2012
 *      Author: Zdenek
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
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

double bessj0(double x)
{
	double ax,z, xx,y,ans,ans1,ans2;

	if ((ax = fabs(x)) < 8.0) {
		y = x*x;
		ans1 = 57568490574.0+y*(-13362590354.0+y*(651619640.7+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2 = 57568490411.0+y*(1029532985.0+y*(9494680.718+y*(59272.64853+y*(267.8532712+y*1.0))));
		ans = ans1/ans2;
	} else {
		z = 8.0/ax;
		y = z*z;
		xx = ax-0.785398164;
		ans1 = 1.0 + y*(-0.1098628627e-2 + y*(0.2734510407e-4 + y*(-0.2073370639e-5 + y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1 + y*(0.1430488765e-3 + y*(-0.6911147651e-5 + y*(0.7621095161e-6 - y*0.934935152e-7)));
		ans = sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}

double bessj1(double x)
{
	double ax,z,xx,y,ans,ans1,ans2;

	if ((ax = fabs(x)) < 8.0) {
		y = x*x;
		ans1 = x*(72362614232.0+y*(-7895059235.0+y*(242396853.1+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
		ans2 = 144725228442.0+y*(2300535178.0+y*(18583304.74+y*(99447.43394+y*(376.9991397+y*1.0))));
		ans	=	ans1/ans2;
	} else {
		z = 8.0/ax;
		y = z*z;
		xx = ax-2.356194491;
		ans1 = 1.0+y*(0.183105e-2+y*(-0.3516396496e-4+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2 = 0.04687499995+y*(-0.2002690873e-3+y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
		ans = sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}

#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

/****
 * Returns the Bessel function Jn(x) for any real x and n >= 0.
 ****/
double bessj(int n, double x)
{
	double ans;

	switch (n) {
	case 0:
		ans = bessj0(x);
		break;
	case 1:
		ans = bessj1(x);
		break;
	default: {
		int j,jsum,m;
		double ax,bj,bjm,bjp,sum,tox;

		ax = fabs(x);
		if (ax == 0.0) {
			return 0.0;
		} else if (ax > (float) n) { //Upwards recurrence
			tox = 2.0/ax;
			bjm = bessj0(ax);
			bj = bessj1(ax);
			for (j=1;j<n;j++) {
				bjp = j*tox*bj-bjm;
				bjm = bj;
				bj = bjp;
			}
			ans = bj;
		} else { //Downwards recurrence
			tox = 2.0/ax;
			m = 2*((n+(int) sqrt(ACC*n))/2);
			jsum = 0;
			bjp = ans = sum = 0.0;
			bj = 1.0;
			for (j=m;j>0;j--) {
				bjm = j*tox*bj-bjp;
				bjp = bj;
				bj = bjm;
				if (fabs(bj) > BIGNO) {
					bj *= BIGNI;
					bjp *= BIGNI;
					ans *= BIGNI;
					sum *= BIGNI;
				}
				if (jsum) sum += bj;
				jsum = !jsum;
				if (j == n) ans=bjp;
			}
			sum = 2.0*sum-bj;
			ans /= sum;
		}
		if (x < 0.0 && (n & 1) ) ans *= -1;
		}
	}
	return ans;
}

/****
 * Lanczos tridiagonalization with partial orthogonalization.
 * IN: A
 * OUT: dg0 main diagonal, dg1 subdiagonal, Q transformation matrix
 ****/
void slanpro(mat_double *A, double *dg0, double *dg1, mat_double *Q)
{
	int i, j, k, p, second = 0, nVec = 0, doOrtho = 0;
	int dim = A->row;
	double d;
	const double eps = 2.2e-16;
	const double SQRTEPS = sqrt(eps);
	const double MIDEPS = pow(eps,3.0/4.0);

	srand(time(NULL));

	// orthogonality estimates
	double *wOld = (double*)malloc(dim*sizeof(double));
	double *wCur = (double*)malloc(dim*sizeof(double));
	wOld[0] = 1.0;
	// upper and lower bounds for orthogonalization estimates
	int *up = (int*)malloc(dim/2*sizeof(int));
	int *low = (int*)malloc(dim/2*sizeof(int));
	int interNum = 0;

	assert(Q->type == MAT_DENSE && dim == Q->row);
printf("\nslanpro\n");
	dm_zero(Q);
	// starting vector of unit norm
	d = 1.0/sqrt(dim);
	for (i=0; i<dim; i++) {
		Q->data[i] = d;
	}
	double *qCur = Q->data;
	double *qOld = NULL;
	double *r = (double *)malloc(dim*sizeof(double));

	for (i=0; i<dim; i++) {
		// r = A*qCurr;
		if (A->type == MAT_DENSE) {
			cblas_dsymv(CblasColMajor,CblasUpper,dim,1.0,A->data,dim,qCur,1,0.0,r,1);
		} else {
			assert(A->type == MAT_SPARSE);
#ifdef INTEL_MKL
			mkl_dcsrgemv("N",&dim,A->data,A->irow,A->icol,qCur,r);
#else
			fprintf(stderr,"Error: slanpro - sparse algebra not compiled\n");
			exit(1);
#endif
		}
		// a[i] = qCur^+ * r
		dg0[i] = cblas_ddot(dim,qCur,1,r,1);
		if (i == dim-1) break; // stop here in the last run
		if (i==0) {
			// r = r - a[i]*qCur
			cblas_daxpy(dim,-dg0[i],qCur,1,r,1);
		} else {
			// r = r - a[i]*qCur - b[i-1]*qOld
			for (j=0;j<dim; j++) {
				r[j] -= ( dg0[i]*qCur[j] + dg1[i-1]*qOld[j] );
			}
		}
		// b[i] = norm(r);
		dg1[i] = cblas_dnrm2(dim,r,1);
		if (dg1[i] <= SQRTEPS) { // b is small
			// orthogonalize r against all previous q
			interNum = 1;
			up[0] = i;
			low[0] = 0;
			doOrtho = 1;
			// update wOld and wCurr
			if (i > 1) {
				cblas_dcopy(i-1,wCur,1,wOld,1);
				wOld[i] = 1.0;
			}
			wCur[i+1] = 1.0;
		} else {
			if (i > 1) {
				// compute orthogonality estimates
				wOld[0] = ( dg1[0]*wCur[1]+dg0[0]*wCur[0]-dg0[i]*wCur[0]-dg1[i-1]*wOld[0] )/dg1[i];
				for (j=1;j<i;j++) {
					wOld[j] = ( dg1[j]*wCur[j+1] + dg0[j]*wCur[j] - dg0[i]*wCur[j] + dg1[j-1]*wCur[j-1] - dg1[i-1]*wOld[j] )/dg1[i];
					wOld[j] += eps*0.3*( dg1[j] + dg1[i] )*(rand()/(double)RAND_MAX - 0.5);
				}
				// swap wOld and wCurr
				for (j=0;j<i; j++) {
					d = wOld[j];
					wOld[j] = wCur[j];
					wCur[j] = d;
				}
				wOld[i] = 1.0;
			}
			wCur[i] = eps*dim*dg1[0]/dg1[i]*0.6*(rand()/(double)RAND_MAX-0.5);
			wCur[i+1] = 1.0;

			if (second == 0) {
				// not the second time, determine intervals
				doOrtho = 0;
				interNum = 0;
				k = 0;
				while (k<=i) {
					if (fabs(wCur[k]) >= SQRTEPS) {
						// lost orthogonality
						doOrtho = 1;
						// find the upper bound
						p = k+1;
						while ( (p<i+1) && ( fabs(wCur[p]) >= MIDEPS) ) p++;
						up[interNum] = p-1;
						// find the lower bound
						p = k-1;
						while ( (p>=0) && ( fabs(wCur[p]) >= MIDEPS) ) p--;
						low[interNum] = p+1;
						k += up[interNum] + 1;
						interNum++;
					} else {
						k++;
					}
				}
			}
		} // if b small
		if (doOrtho || second) {
			printf("orthog. in iter %d (%d, %d)\n",i,doOrtho,second);
			printf("\t input vector:\n");
			int fff;
			for (fff=0; fff<dim; fff++) printf("%g ",r[fff]);
			printf("\n");
			// do orthogonalization
			for (j=0; j<interNum; j++) { // for each interval
				printf("\t (%d, %d)\n",low[j],up[j]);
				for (p=low[j]; p<=up[j]; p++) {
					// reset ortho estimates
					wCur[p] = eps*1.5*(rand()/(double)RAND_MAX-0.5);
					// orthogonalization
					d = cblas_ddot(dim,Q->data+p*dim,1,r,1);
					printf("\t\t d = %g\n",d);
					cblas_daxpy(dim,-d,Q->data+p*dim,1,r,1);
				}
				// count the number of orthogonalizations performed
				nVec += up[j]-low[j]+1;
				if (second == 0) {
					// adjust intervals for the second time
					if (low[j]-1 >= 0) low[j]--; else low[j]=0;
					if (up[j] <= i) up[j]++; else up[j]=i+1;
				}
			}
			// recalculate b[i]
			dg1[i] = cblas_dnrm2(dim,r,1);
			// set logicals
			if (second == 1) {
				second = 0;
			} else {
				second = 1;
				doOrtho = 0;
			}
			printf("\t output vector of norm = %g:\n",dg1[i]);
			for (fff=0; fff<dim; fff++) printf("%g ",r[fff]);
			printf("\n");
		}
		// store the vector to Q
		qOld = qCur;
		qCur = Q->data+(i+1)*dim;
		cblas_daxpy(dim,1.0/dg1[i],r,1,qCur,1);

	} // main loop (i)

	// cleaning
	free(r);
	free(up);
	free(low);
	free(wOld);
	free(wCur);

	printf("slanpro did %d orthogonalizations\n",nVec);

}

