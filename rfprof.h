/*
 * rfprof.h
 *
 *  Created on: Jun 1, 2010
 *      Author: zdenek
 */

#ifndef RFPROF_H_
#define RFPROF_H_

double ** rfprof_alloc(int N, int chan);
void rfprof_free(double ** obj);
double ** read_rfproffile(const char * fname, int chan);
int rfprof_len(double ** obj);
double rfprof_sumweight(double **rfdata);

#endif /* RFPROF_H_ */
