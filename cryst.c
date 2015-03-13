/*
    Reads powder averaging/crystallite files
    Copyright (C) 1999 Mads Bak

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
    
    
    The crystallite file must have the following data format:

       <N>
       <alpha1> <beta1> <weight1>
       <alpha2> <beta2> <weight2>
       <alpha3> <beta3> <weight3>
       ...
       <alphaN> <betaN> <weightN>

    where <N> is the number of crystallites, alpha is the 
    alpha Euler angle ranging from 0 to 360 degrees, beta is
    the beta Euler angle ranging from 0 to 180 degrees.
    Weight is the solid angle covered by the crystallite
    and all weights must sum up to 1. Be avare that the number
    of significant digits must large, otherwise errors can be
    introduced.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "cryst.h"
#include "defs.h"
#include "matrix.h"

Cryst * cryst_alloc(int N)
{
	Cryst * crdata;

	crdata = (Cryst*)malloc((N+1)*sizeof(Cryst));
    *(int*)crdata = N;
    return crdata;
}

void cryst_free(Cryst * crdata)
{
	free((char*)crdata);
}

Cryst * read_crystfile(char* crystname, int from, int to)
{
  FILE* fp;
  char fname[256], line[256];
  CRYSTALLITE* c;
  int N,i,j;
  int n=0,ver;
  Cryst *crdata;
  int cnt, cnt_old;
  double da, db, dg, dw;
  int isbin = 0;
  int a = 2;
  
  ver= verbose & VERBOSE_POWDER;
  strcpy(fname,crystname);
  strcat(fname,"_cryst");
  while (strlen(cryst_names[n])) {
    if (!strcmp(cryst_names[n],fname)) {
       if (ver) {
         printf("found internal crystallite file '%s'\n",crystname);
         printf("to overwrite with external file, specify './%s' instead\n",crystname);
       }
       N=cryst_numbers[n];
       if (from > N) {
    	   fprintf(stderr,"Error: bad crystallite range - initial value too large/n");
    	   exit(1);
       }
       if (to > N || to < 0) to = N;
       N = to - from + 1;
       assert( N > 0 );
       if (ver) printf("crystallites (%d - %d): %d\n",from, to, N);
       crdata = cryst_alloc(N);
       c=cryst_pointers[n];
       for (i=from-1,j=1;i<to;i++,j++) {
/* Set the powder angles, gamma is set to zero for the internal powder files */
         crdata[j].alpha = c[i].alpha;
         crdata[j].beta = c[i].beta;
         crdata[j].gamma = 0.0;
         crdata[j].weight = c[i].weight;
         if (c[i].weight == 0.0) {
           fprintf(stderr,"error: crystallite number %d in file '%s' has zero weight\n",i+1,crystname);
           exit(1);
         }
         if (ver) 
      	  printf("%5d %15g %15g %15g %15g\n",j, c[i].alpha,c[i].beta,0.0,c[i].weight);
       }
       return crdata;
    }
    n++;
  }

  strcpy(fname,crystname);

#ifdef UNIX
  if (name[0] == '~') {
    char* p=getenv("HOME");
    if (p != NULL) {
      strcpy(fname,p);
      strcat(fname,&name[1]);
    }
  }
#endif
  fp = fopen(fname,"rb");
  if (!fp) {
	  strcpy(line,fname);
	  strcat(line,".cry");
	  fp = fopen(line,"r");
  }
  if (!fp) {
	  strcpy(line,fname);
	  strcat(line,".crb");
	  fp = fopen(line,"rb");
  }
  if (!fp) {
    int n=0;
    int l,tl=0;

    fprintf(stderr,"error: unable to open file '%s'\n",fname);
    fprintf(stderr,"internal crystallite files are: \n");
    while (strlen(cryst_names[n])) {
       char nam[256];
       strcpy(nam,cryst_names[n++]);
       l=strlen(nam)-6;
       nam[l]=0;
       tl += l;
       fprintf(stderr,"%s ",nam);
       if (tl > 32) {
         fprintf(stderr,"\n");
         tl=0;
       }
    }
    fprintf(stderr,"\n");
    exit(1);
  }
  
  if (ver) printf("loading external crystallite file '%s'\n",fname);

  fgets(line,256,fp);
  if (!strncmp(line,"Orientations:",13)) {
	  isbin = 1;
	  if (sscanf(&line[13],"%d",&N) != 1) {
		  fprintf(stderr,"Error: unable to read number of crystallites from file '%s'\n",fname);
		  exit(1);
	  }
	  //printf("filepos 2: %ld\n",ftell(fp));
	  fgets(line,256,fp);
	  //printf("filepos 2: %ld\n",ftell(fp));
	  if (strncmp(line,"Angles      :",13)) {
		  fprintf(stderr,"Error: wrong format of file '%s'\n",fname);
		  exit(1);
	  }
	  if (sscanf(&line[13],"%d",&a) != 1) {
		  fprintf(stderr,"Error: unable to read number of crystallites from file '%s'\n",fname);
		  exit(1);
	  }
	  if (a != 2 && a != 3) {
		  fprintf(stderr,"Error: wrong number of angles in file '%s'\nShould be 2 or 3 but read %d\n",fname,a);
		  exit(1);
	  }
  } else {
	  isbin = 0;
	  if (sscanf(line,"%d",&N) != 1) {
	    fprintf(stderr,"unable to read number of crystallites from file '%s'\n",fname);
	    exit(1);
	  }
  }
  if (from > N) {
	   fprintf(stderr,"Error: bad crystallite range - initial value too large/n");
	   exit(1);
  }
  if (to > N || to < 0) to = N;
  N = to - from + 1;
  assert( N > 0 );
  if (ver) printf("crystallites (%d - %d): %d\n",from, to, N);
  crdata = cryst_alloc(N);

  if (isbin == 0) {
	  for (i=1;i<from;i++) {
		  if (!fgets(line, 256, fp) ) {
			  fprintf(stderr,"Error: unable to read line %d in file %s\n",i,fname);
			  exit(1);
		  }
	  }
	  for (i=1;i<=N;i++) {
		  /* read in 3 or 4 numbers, just decide on the content of the file */
		  if (!fgets(line, 256, fp) ) {
			  fprintf(stderr,"Error: unable to read line %d in file %s\n",i+from-1,fname);
			  exit(1);
		  }
		  cnt = sscanf(line,"%lg%lg%lg%lg",&da,&db,&dg,&dw);
		  if ( (cnt != 3) && (cnt != 4) ) {
			  fprintf(stderr,"Error: wrong number of parameters on line %d in file %s\n",i,fname);
			  exit(1);
		  }
		  if (i == 1)
			  cnt_old = cnt;
		  else {
			  if (cnt != cnt_old) {
				  fprintf(stderr,"Error: wrong number of parameters on line %d in file %s\n",i,fname);
				  exit(1);
			  }
		  }
		  if (cnt == 3) {
			  dw = dg;
			  dg = 0.0;
		  }
		  if (dw == 0.0) {
			  fprintf(stderr,"error: crystallite number %d in file '%s' has zero weight\n",i,crystname);
			  exit(1);
		  }
		  crdata[i].alpha =da;
		  crdata[i].beta = db;
		  crdata[i].gamma = dg;
		  crdata[i].weight = dw;

		  if (ver)
			  printf("%5d %15g %15g %15g %15g\n",i, da,db,dg,dw);
	  }
  } else { // binary format
	  fseek(fp,(from-1)*(a+1)*sizeof(double),SEEK_CUR);
	  for (i=1; i<=N; i++) {
		  fread(&da,sizeof(double),1,fp);
		  fread(&db,sizeof(double),1,fp);
		  if (a == 3) fread(&dg,sizeof(double),1,fp); else dg = 0.0;
		  fread(&dw,sizeof(double),1,fp);
		  crdata[i].alpha =da;
		  crdata[i].beta = db;
		  crdata[i].gamma = dg;
		  crdata[i].weight = dw;
		  if (ver)
			  printf("%5d %15g %15g %15g %15g\n",i, da,db,dg,dw);
	  }
  }

  fclose(fp);
  return crdata;
}

double cryst_sumweight(Cryst *crdata)
{
	double w=0.0;
	int i,N;

	N = LEN(crdata);
	for (i=1; i<=N; i++)
		w += crdata[i].weight;

	return w;
}



TRIANGLE * triangle_alloc(int N)
{
	TRIANGLE * tria;

	tria = (TRIANGLE*)malloc((N+1)*sizeof(TRIANGLE));
    *(int*)tria = N;
    return tria;
}

void triangle_free(TRIANGLE *tria)
{
	free((char*)tria);
}

TRIANGLE * read_triangle_file(char* name)
{
  FILE* fp;
  char fname[256], line[256];
  int ver, N, i, cnt, isbin;
  TRIANGLE *tria;

  ver= verbose & VERBOSE_POWDER;
  strcpy(fname,name);

#ifdef UNIX
  if (name[0] == '~') {
    char* p=getenv("HOME");
    if (p != NULL) {
      strcpy(fname,p);
      strcat(fname,&name[1]);
    }
  }
#endif
  // assumes fname is without suffix !!!
  strcpy(line,fname);
  strcat(line,".tri");
  isbin = 0;
  fp=fopen(line,"r");
  if (!fp) {
	  strcpy(line,fname);
	  strcat(line,".trb");
	  fp=fopen(line,"rb");
	  isbin = 1;
  }

  if (!fp) {
    fprintf(stderr,"error: unable to open file '%s' (.tri nor .trb)\n",fname);
    exit(1);
  }

  if (ver) printf("loading external crystallite triangles file '%s'\n",fname);

  /* scan file for number of lines, i.e. number of elements in rf shape */
  N = 0;
  if (isbin == 0) { // text file
	  while ( fgets(line, 256, fp) ) {
		  N++;
	  }
  } else { // binary file
	  fseek(fp,0,SEEK_END);
	  cnt = ftell(fp);
	  if (cnt % (3*sizeof(int)) != 0) {
		  fprintf(stderr,"error: triangle file '%s' has wrong size, possibly corrupted...\n",fname);
		  exit(1);
	  }
	  N = cnt / (3*sizeof(int));
  }
  fseek(fp, 0, SEEK_SET);
  if (ver) printf("Number of lines = %d\n",N);

  tria = triangle_alloc(N);
  for (i=1;i<=N;i++) {
	  if (isbin == 0) { // text file
		  // read in just 3 integers from each line
		  if (!fgets(line, 256, fp) ) {
			  fprintf(stderr,"Error: unable to read line %d in file %s\n",i,fname);
			  exit(1);
		  }
		  cnt = sscanf(line,"%d%d%d",&tria[i].a,&tria[i].b,&tria[i].c);
		  if ( cnt != 3 ) {
			  fprintf(stderr,"Error: can not read 3 integers on line %d in file %s\n",i,fname);
			  exit(1);
		  }
	  } else { // binary file
		  fread(&tria[i].a,sizeof(int),1,fp);
		  fread(&tria[i].b,sizeof(int),1,fp);
		  fread(&tria[i].c,sizeof(int),1,fp);
	  }
	  if (ver) printf("%5d %5d %5d %5d\n",i, tria[i].a,tria[i].b,tria[i].c);
  }

  fclose(fp);
  return tria;
}

int * read_cryst_map(char *crystfile, Cryst *crdata, char *targetcrystfile, Cryst *targetcrdata)
{
	FILE* fp;
	char fname[256], line[32];
	int ver, i, j, ncr, ntcr, k, l;
	int *crmap, isbin;
	double d, dmin;

	ver = verbose & VERBOSE_POWDER;
	ncr = LEN(crdata);
	ntcr = LEN(targetcrdata);
	crmap = (int *)malloc(ntcr*sizeof(int));
	if (crmap == NULL) {
		fprintf(stderr,"Error: can not allocate crmap\n");
		exit(1);
	}
	sprintf(fname,"%s_%s.mab",targetcrystfile,crystfile);
    fp = fopen(fname,"rb");
    isbin = 1;
    if (fp != NULL) {
        fseek(fp,0,SEEK_END);
        k = ftell(fp);
        if ( ( k/sizeof(int) != ntcr ) || (k%sizeof(int) != 0) ) {
        	fprintf(stderr,"Error: read_cryst_map: wrong size of the binary file %s",fname);
        	exit(1);
        }
        fseek(fp,0,SEEK_SET);
    } else { // fp == NULL
    	sprintf(fname,"%s_%s.map",targetcrystfile,crystfile);
        fp = fopen(fname,"r");
        isbin = 0;
    }

	if (fp != NULL) { // load from file if the data file exist
		if (ver) printf("Reading crystallites NEAREST map from %s file %s\n",isbin ? "binary" : "text",fname);
		for (i=0; i<ntcr; i++) {
			if (isbin == 0) {
				if (!fgets(line, 32, fp) ) {
					fprintf(stderr,"Error: unable to read line %d in file %s\n",i,fname);
					exit(1);
				}
				j = sscanf(line,"%d%d",&k,&l);
				if ( j != 2 ) {
					fprintf(stderr,"Error: can not read 2 integers on line %d in file %s\n",i,fname);
					exit(1);
				}
				crmap[k-1] = l;
			    if (ver) printf("%5d %5d\n", k, l);
			} else {
				fread(&crmap[i],sizeof(int),1,fp);
			    if (ver) printf("%5d %5d\n", i+1, crmap[i]);
			}
		}
	} else { // create map de novo and store it to file
		sprintf(fname,"%s_%s",targetcrystfile,crystfile);
		if (ver) printf("Creating crystallites NEAREST map and saving to file %s.map/.mab\n",fname);
		strcat(fname,".map");
		fp = fopen(fname,"w");
		for (i=1; i<=ntcr; i++) {
			dmin = 1e99;
			for (j=1; j<=ncr; j++) {
				d = fabs(crdata[j].alpha - targetcrdata[i].alpha) + fabs(crdata[j].beta - targetcrdata[i].beta) + fabs(crdata[j].gamma - targetcrdata[i].gamma);
				if (d < dmin) {
					dmin = d;
					l = j;
				}
			}
			if (dmin>1e90) {
				fprintf(stderr,"Error: read_cryst_map - ups, can not find nearest for %d\n",i);
				exit(1);
			}
			crmap[i-1] = l;
			fprintf(fp,"%d %d\n",i,l);
		    if (ver) printf("%5d %5d\n", i, l);
		}
		fclose(fp);
		sprintf(fname,"%s_%s.mab",targetcrystfile,crystfile);
		fp = fopen(fname,"wb");
		for (i=0; i<ntcr; i++) fwrite(&crmap[i],sizeof(int),1,fp);
		fclose(fp);
	}
	return crmap;
}


void save_bin_crystfile(char* filename, Cryst * crdata)
{
	FILE* fp;
	int i, a;
	int N = LEN(crdata);
	double g;
	char fname[256];

	// decide if it is 2 or 3 angle set
	g = 0.0;
	for (i=1; i<=N; i++) {
		g += fabs(crdata[i].gamma);
	}
	if (g > 1e-6) {
		a = 3;
	} else {
		a = 2;
	}

	strcpy(fname,filename);
	strcat(fname,".crb");
	fp = fopen(fname,"wb");
	if (!fp) {
		fprintf(stderr,"save_bin_crystfile: unable to create file '%s'\n",fname);
	    exit(1);
	}
	fprintf(fp,"Orientations: %d\n",N);
	fprintf(fp,"Angles      : %d\n",a);
	for (i=1; i<=N; i++) {
		fwrite(&crdata[i].alpha,sizeof(double),1,fp);
		fwrite(&crdata[i].beta,sizeof(double),1,fp);
		if (a == 3) fwrite(&crdata[i].gamma,sizeof(double),1,fp);
		fwrite(&crdata[i].weight,sizeof(double),1,fp);
   	  //printf("%5d %15g %15g %15g %15g\n",i, crdata[i].alpha,crdata[i].beta,crdata[i].gamma,crdata[i].weight);
	}
	fclose(fp);
}

void save_bin_triangle_file(char* filename,TRIANGLE *tria)
{
	FILE* fp;
	int i;
	int N = LEN(tria);
	char fname[256];

	strcpy(fname,filename);
	strcat(fname,".trb");
	fp = fopen(fname,"wb");
	if (!fp) {
		fprintf(stderr,"save_bin_triangle_file: unable to create file '%s'\n",fname);
	    exit(1);
	}
	for (i=1; i<=N; i++) {
		fwrite(&tria[i].a,sizeof(int),1,fp);
		fwrite(&tria[i].b,sizeof(int),1,fp);
		fwrite(&tria[i].c,sizeof(int),1,fp);
	}
	fclose(fp);
}

void save_bin_cryst_map(char *crystfile, char *targetcrystfile, int *map, int N)
{
	FILE* fp;
	int i;
	char fname[256];

	sprintf(fname,"%s_%s.mab",targetcrystfile,crystfile);
    fp = fopen(fname,"wb");
	if (!fp) {
		fprintf(stderr,"save_bin_cryst_map: unable to create file '%s'\n",fname);
	    exit(1);
	}
	for (i=0; i<N; i++) {
		fwrite(&map[i],sizeof(int),1,fp);
	}
	fclose(fp);
}
