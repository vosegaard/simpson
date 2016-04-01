/*
 * wigner.h
 *
 *  Created on: Jun 2, 2010
 *      Author: zdenek
 */

#ifndef WIGNER_H_
#define WIGNER_H_

complx* Dtensor2(double deltazz,double eta);
mat_complx * wigner2(double alpha,double beta,double gamma);
complx * wigner20(double alpha,double beta);
void wig2rot(complx* res, complx* vec, mat_complx *d2);
void wig2rot_t(complx* res, complx* vec, mat_complx *d2);
complx * wig2roti(complx *vec, double a, double b, double c);
double wig20rot(complx *vec, complx *d20);

#endif /* WIGNER_H_ */
