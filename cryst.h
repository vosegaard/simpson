#ifndef __CRYST_H
#define __CRYST_H

typedef struct _CRYSTALLITE {
  double alpha,beta,weight;
} CRYSTALLITE;

extern char* cryst_names[];
extern int cryst_numbers[];
extern CRYSTALLITE* cryst_pointers[];

typedef struct _TRIANGLE {
  int a, b, c;
  //double weight;
} TRIANGLE;

typedef struct _Cryst {
   double alpha, beta, gamma, weight;
} Cryst;

Cryst * cryst_alloc(int N);
void cryst_free(Cryst * crdata);
Cryst * read_crystfile(char* crystname, int from, int to);
double cryst_sumweight(Cryst *crdata);
int * read_cryst_map(char *crystfile, Cryst *crdata, char *targetcrystfile, Cryst *targetcrdata);

TRIANGLE *triangle_alloc(int N);
void triangle_free(TRIANGLE * triadata);
TRIANGLE * read_triangle_file(char* trianame);
void save_bin_crystfile(char* filename, Cryst * crdata);
void save_bin_triangle_file(char* filename,TRIANGLE *tria);
void save_bin_cryst_map(char *crystfile, char *targetcrystfile, int *map, int N);

#endif
