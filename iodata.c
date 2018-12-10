/*
    Data format manipulation routines
    Copyright (C) 1999 Mads Bak, 2001 Thomas Vosegaard

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
    
    Routines for internal manipulation with the SIMPSON format.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <tcl.h>
#include <sys/types.h>
#include <unistd.h>
#include "iodata.h"

/* NMRpipe specific inclusion */
#include "fdatap.h"

#ifdef __APPLE__
#  ifndef __i386__
#    define REVERSEBYTES 1
#  endif
#endif

typedef union {
   char bits[4];

#ifdef WORDS_BIGENDIAN
   struct {
      unsigned int negative:1;
      unsigned int exponent:8;
      unsigned int mantissa:23;
   } f;
#else
   struct {
      unsigned int mantissa:23;
      unsigned int exponent:8;
      unsigned int negative:1;
   } f;
#endif

} floatbits;

#define FLOATBITS_EXPADD	0x7f
#define FLOATBITS_MANT    23

void double_to_bits(double fd, char *bits) {
   float f = 0;
   
   if ((fd > __FLT_MIN__ && fd < __FLT_MAX__) || (-fd > __FLT_MIN__ && -fd < __FLT_MAX__)) {
      f = (float)fd;
   }
   float_to_bits(f, bits);
}

void float_to_bits(float f,char* bits)
{
   int e;
   double m;
   floatbits* fb=(floatbits*)bits;
   
   m=frexp(f,&e);
   fb->f.negative = (f < 0);
   fb->f.exponent = e + FLOATBITS_EXPADD;
   fb->f.mantissa = fabs(m)*(double)(1<<FLOATBITS_MANT);

#ifdef WORDS_BIGENDIAN
   {
     char tmp;
     tmp=bits[3];
     bits[3]=bits[0];
     bits[0]=tmp;

     tmp=bits[2];
     bits[2]=bits[1];
     bits[1]=tmp;
   }
#endif

}

float bits_to_float(char* bits)
{
   int e;
   double m;
#ifdef WORDS_BIGENDIAN
   char bits_rev[4];
   floatbits* fb=(floatbits*)bits_rev;
   bits_rev[0]=bits[3];
   bits_rev[1]=bits[2];
   bits_rev[2]=bits[1];
   bits_rev[3]=bits[0];
#else
   floatbits* fb=(floatbits*)bits;
#endif

   e = fb->f.exponent - FLOATBITS_EXPADD;
   m = (double)fb->f.mantissa/(double)(1<<FLOATBITS_MANT);
   if (fb->f.negative) {
     return -ldexp(m,e);
   }
   return ldexp(m,e);
}


#define FIRST(f,x) ((x) & ~(~0 << f))
#define LAST(f,x) ((x) & (~0 << (8-f)))
/* remove spaces */
#define BASE 33 

int _pack_newline;
int _pack_nbuf = -1;
unsigned char _pack_buf[4];

void pack_begin(FILE *fp)
{
  if (_pack_nbuf != -1) {
    fprintf(stderr,"error: pack_begin(): previous session was not ended with pack_end()\n");
    exit(-1);
  }
  _pack_nbuf=0;
  _pack_newline=0;
}

int pack_putc(int c,FILE* fp)
{
  _pack_buf[_pack_nbuf++] = c;
  if (_pack_nbuf == 3) {
    putc(FIRST(6,_pack_buf[0]) + BASE,fp);
    putc(((LAST(2,_pack_buf[0]) >> 2) | FIRST(4,_pack_buf[1]))+BASE,fp);
    putc(((LAST(4,_pack_buf[1]) >> 2) | FIRST(2,_pack_buf[2]))+BASE,fp);
    if (putc((LAST(6,_pack_buf[2]) >> 2)+BASE,fp) == EOF) {
      _pack_nbuf= -1;
      return EOF;
    }
    _pack_newline++;
    if (_pack_newline == 16) {
      putc('\n',fp);
      _pack_newline=0;
    }
    _pack_nbuf=0;
  }  
  return 1;  
}


int pack_end_putc(FILE* fp)
{
  if (_pack_nbuf > 0) {
    int i;

    for (i=_pack_nbuf;i<3;i++) {
      if (pack_putc(0,fp) == EOF) return EOF;
    }
    if (_pack_newline != 0)
      putc('\n',fp);    
  }
  _pack_nbuf=-1;
  return 1;
}

int pack_getc(FILE* fp)
{
  if (_pack_nbuf == 0) {
     int c;
     unsigned char c0,c1,c2,c3;

     c=getc(fp);
     if (c == '\n')
       c=getc(fp);
     c0=c-BASE;
     c1=getc(fp)-BASE;
     c2=getc(fp)-BASE;     
     if ((c=getc(fp)) == EOF) {
       _pack_nbuf= -1;
       return EOF;
     }
     c3=c-BASE;
     _pack_buf[0]=FIRST(6,c0) | LAST(2,c1 << 2);
     _pack_buf[1]=FIRST(4,c1) | LAST(4,c2 << 2);
     _pack_buf[2]=FIRST(2,c2) | LAST(6,c3 << 2);
     _pack_nbuf=3;
  }
  return _pack_buf[3-_pack_nbuf--];
}

int pack_end_getc(FILE* fp)
{
  _pack_nbuf=-1;  
  return 1;
}


/*
  Data Format Conventions:

  Arrays are 1 offset
*/

/*
This is a test of
the new dataformat.:

SIMP
NP=4
SW=10000
REF=100
NI=2
SW1=20000
REF1=200
DATA
1 0
2.3 0
3.1 0
4.01 0
1 0
2.3 0
3.1 0
4.01 0
END
*/

void FD_write(char* fname,FD* fd,int format,int prec)
{
  int i,ntot;
  FILE* fp;
  char bits[4];
  
  if (*fname != '-') {
    fp=fopen(fname,"w");
  } else {
    fp=stdout;
  }
  if (!fp) {
    fprintf(stderr,"FD_write: unable to create file '%s'\n",fname);
    exit(1);
  }
  fprintf(fp,"SIMP\n");  
  if (fd->np <= 0) {
    fprintf(stderr,"FD_write: number of points with must be larger"
                    " than zero when writing file '%s'\n",fname);
    exit(1);
  }
  fprintf(fp,"NP=%d\n",fd->np);
  if (fd->sw <= 0.0) {
    fprintf(stderr,"FD_write: spectral-width must be larger"
                   " than zero when writing file '%s'\n",fname);
    exit(1);
  }
  fprintf(fp,"SW=%.12g\n",fd->sw);  
  if (fd->ref != 0.0)
    fprintf(fp,"REF=%.12g\n",fd->ref);

  if (fd->ni != 0) {
    fprintf(fp,"NI=%d\n",fd->ni);
  }
  if (fd->sw1 != 0.0)
    fprintf(fp,"SW1=%.12g\n",fd->sw1);  
  if (fd->ref1 != 0.0)
    fprintf(fp,"REF1=%.12g\n",fd->ref1);  
  if (fd->nelem != 1)
    fprintf(fp,"NELEM=%d\n",fd->nelem);
  if (fd->type == FD_TYPE_FID)
    fprintf(fp,"TYPE=FID\n");
  else
    fprintf(fp,"TYPE=SPE\n");
  fd->format=format;
  if (fd->format == FD_FORMAT_BINARY) {
    fprintf(fp,"FORMAT=BINARY\n");  
    fd->prec=prec;
    if (fd->prec == FD_PREC_DOUBLE)
      fprintf(fp,"PREC=DOUBLE\n");    
  }
  fprintf(fp,"DATA\n");
  ntot=fd->np*(fd->ni > 1 ? fd->ni : 1)*fd->nelem;

  if (fd->format == FD_FORMAT_BINARY) {    
    if (fd->prec == FD_PREC_SINGLE) {
       int i,doubles;
       double* data;

       doubles=ntot*2;
       data = &(fd->data)[1].re;
       pack_begin(fp);
       
       for (i=0;i<doubles;i++) {
         int j;         
         double_to_bits(data[i],bits);
         for (j=0;j<4;j++) {
           if (pack_putc(bits[j],fp) == EOF) {
             fprintf(stderr,"error: FD_write: cannot write float %d\n",i);
             exit(1);  
           }
         }
       }
       if (pack_end_putc(fp) == EOF) {
         fprintf(stderr,"error: FD_write: cannot write last float\n");
         exit(1);  
       }       
    } else {
    }
  } else {
    for (i=1;i<=ntot;i++)
      fprintf(fp,"%.9g %.9g\n",fd->data[i].re,fd->data[i].im);
  }
  fprintf(fp,"END\n");
  if (*fname != '-') {
    fclose(fp);
  }
}



int FD_alloc1ddata(FD* fd)
{
  int ntot;
  ntot=fd->np*(fd->ni > 1 ? fd->ni : 1)*fd->nelem;
  fd->data=(complx*)malloc(sizeof(complx)*(ntot+1));
  if (!fd->data) {
    fprintf(stderr,"FD_alloc1ddata: unable to allocate 2x%d doubles\n",ntot);
    return 0;
  }
  return 1;
}


void FD_defvalues(FD* fd)
{
  if (!fd) {
    fprintf(stderr,"FD_defvalues: called with NULL argument\n");
    exit(1);
  }  
  fd->data=NULL;
  fd->np=0;
  fd->ni=0;
  fd->nelem = 1;
  fd->type=FD_TYPE_FID;
  fd->format=FD_FORMAT_TEXT;
  fd->prec=FD_PREC_SINGLE;
  fd->ref=0.0;
  fd->sw=0.0;
  fd->ref1=0.0;
  fd->sw1=0.0;
}



FD* FD_alloc()
{
  FD* fd;
  fd = (FD*)malloc(sizeof(FD));
  if (!fd) {
    fprintf(stderr,"FD_alloc: unable to allocate FD data structure\n");
    exit(1);
  }  
  FD_defvalues(fd);
  return fd;  
}

void FD_free(FD* fd)
{
  if (!fd) {
    fprintf(stderr,"FD_free: called with NULL argument\n");
    exit(1);
  }  
  if (!fd->data) {
    fprintf(stderr,"FD_free: FD structure contained no data !\n");
    exit(1);
  }
  free(fd->data);
  free(fd);
}

void FD_zero(FD* fd)
{
  if (!fd) {
    fprintf(stderr,"FD_zero: called with NULL argument\n");
    exit(1);
  }  
  if (!fd->data) {
    fprintf(stderr,"FD_zero: FD structure contained no data !\n");
    exit(1);
  }
  memset(fd->data,0,sizeof(complx)*(fd->np*(fd->ni > 1 ? fd->ni : 1)*fd->nelem+1));
}


FD* FD_read(char* fname)
{
  int got_np=0,got_sw=0,got_data=0;
  FILE* fp;
  FD* fd;
  int i,lineno,ntot;
  char buf[256],*p,*q;
  char com[256];
  char val[256];
  char bits[4];

  fd=FD_alloc();
  if (!fd){
    fprintf(stderr,"FD_read: unable to allocate fd datatype\n");
    exit(1);
  }
  if (!fname) {
    fprintf(stderr,"FD_read: filename pointer is null\n");
    exit(1);
  }  
  fp=fopen(fname,"r");
  if (!fp) {
    fprintf(stderr,"FD_read: unable to open file '%s'\n",fname);
    exit(1);
  }
  lineno=1;
  if (!fgets(buf,sizeof(buf),fp)) {
    fprintf(stderr,"FD_read: unable to read line %d from file '%s'\n",lineno,fname);
    exit(1);
  }
  // this should fix mixed line-ending errors
  if (strcmp(buf,"SIMP\n") != 0 && strcmp(buf,"SIMP\r\n") != 0) {
    fprintf(stderr,"FD_read: invalid type of file '%s' (expected 'SIMP'), in file '%s'\n",buf,fname);
    exit(1);
  } 
  while (!feof(fp) && !got_data) {
    lineno++;
    if (!fgets(buf,sizeof(buf),fp)) {
      if (got_sw && got_data && got_np) break;
      buf[0]=0;
      if (!got_sw) {
        strcat(buf," SW=<spectralwidth>,");
      }
      if (!got_data) {
        strcat(buf," datasection,");
      }
      if (!got_np) {
        strcat(buf," NP=<number of points>,");
      }
      fprintf(stderr,"FD_read: missing%s in file '%s'\n",buf,fname);
      exit(1);
    }
    p=buf;
    q=p+strlen(buf)-1;
    if (*q == '\n') *q=0;
    while (isspace(*p)) p++;
    if (*p == 0) continue;
    q=(char*)strchr(p,'=');    
    if (q) {
       *q=0; q++;
       while (isspace(*p)) p++;

       strcpy(com,p);
       strcpy(val,q);

       if (!strcmp(com,"NP")) {
          fd->np=atoi(val);
          if (fd->np <= 0) {
            fprintf(stderr,"FD_read: number of points NP must be larger than zero\n");
            exit(1);
          }
          got_np=1;
       } else if (!strcmp(com,"NI")) {
          fd->ni=atoi(val);
          if (fd->ni <= 0) {
            fprintf(stderr,"FD_read: number of points NI must be larger than zero\n");
            exit(1);
          }
       } else if (!strcmp(com,"NELEM")) {
           fd->nelem=atoi(val);
           if (fd->nelem <= 0) {
             fprintf(stderr,"FD_read: number of elements NELEM must be larger than zero\n");
             exit(1);
           }
       } else if (!strcmp(com,"SW")) {
          if (sscanf(val,"%lg",&fd->sw) != 1) {
            fprintf(stderr,"FD_read: unable to convert '%s' to double in file '%s'\n",val,fname);
            exit(1);
          }
          got_sw=1;
       } else if (!strcmp(com,"SW1")) {
          if (sscanf(val,"%lg",&fd->sw1) != 1) {
            fprintf(stderr,"FD_read: unable to convert '%s' to double in file '%s'\n",val,fname);
            exit(1);
          }
       } else if (!strcmp(com,"TYPE")) {
          if (!strcmp(val,"FID")) {
            fd->type=FD_TYPE_FID;
          } else if (!strcmp(val,"SPE")) {
            fd->type=FD_TYPE_SPE;
          } else {
            fprintf(stderr,"FD_read: unable to convert '%s' to SPE or FID in file '%s'\n",val,fname);
            exit(1);
          }
       } else if (!strcmp(com,"FORMAT")) {
          if (!strcmp(val,"TEXT")) {
            fd->format=FD_FORMAT_TEXT;
          } else if (!strcmp(val,"BINARY")) {
            fd->format=FD_FORMAT_BINARY;
          } else {
            fprintf(stderr,"FD_read: unable to convert '%s' to TEXT or BINARY in file '%s'\n",val,fname);
            exit(1);
          }
       } else if (!strcmp(com,"PREC")) {
          if (!strcmp(val,"DOUBLE")) {
            fd->prec=FD_PREC_DOUBLE;
          } else if (!strcmp(val,"SINGLE")) {
            fd->prec=FD_PREC_SINGLE;
          } else {
            fprintf(stderr,"FD_read: unable to convert '%s' to DOUBLE or SINGLE in file '%s'\n",val,fname);
            exit(1);
          }
       } else if (!strcmp(com,"REF")) {
          if (sscanf(val,"%lg",&fd->ref) != 1) {
            fprintf(stderr,"FD_read: unable to convert '%s' to double in file '%s'\n",val,fname);
            exit(1);
          }
       } else if (!strcmp(com,"REF1")) {
          if (sscanf(val,"%lg",&fd->ref1) != 1) {
            fprintf(stderr,"FD_read: unable to convert '%s' to double in file '%s'\n",val,fname);
            exit(1);
          }
       } else {
/*
         fprintf(stderr,"FD_read: unknown command '%s' in file '%s'\n",com,fname);
         exit(1);
*/
       }
    } else {
       if (strcmp(buf,"DATA")) {
         fprintf(stderr,"FD_read: data section must start with DATA and not '%s' in file '%s'\n",buf,fname);
         exit(1);
       }
       FD_alloc1ddata(fd);
       ntot=fd->np*(fd->ni > 1 ? fd->ni : 1)*fd->nelem;

       if (fd->format == FD_FORMAT_TEXT) {
         for (i=1;i<=ntot;i++) {
           if (fscanf(fp,"%lg%lg",&(fd->data[i].re),&(fd->data[i].im)) != 2) {
              fprintf(stderr,"FD_read: unable to read data point number %d from file '%s'\n",i,fname);
              exit(1);
           }
         }
       } else {
         if (fd->prec == FD_PREC_SINGLE) {
           int i,doubles,c;
           double* data;

           doubles=ntot*2;
           data = &(fd->data)[1].re;
           pack_begin(fp);
           for (i=0;i<doubles;i++) {
             int j;
             for (j=0;j<4;j++) {
               if ((c=pack_getc(fp)) == EOF) {
                 fprintf(stderr,"error: FD_read: cannot read float %d\n",i);
                 exit(1);  
               }
               bits[j]=c;
             }
             data[i]=bits_to_float(bits);
           }
           pack_end_getc(fp);
         } else {
           int nb;
           if ((nb=fread(&(fd->data[1]),sizeof(complx),ntot,fp)) != ntot) {
             fprintf(stderr,"FD_read: read %d but needed %d complex numbers"
                            " from file '%s'\n",nb,ntot,fname);
             exit(1);           
           }
         }
       }
       strcpy(buf,"");
       if (!fscanf(fp,"%s",buf)) {
         fprintf(stderr,"FD_read: must read a string (\"END\") at the end of data\n");
         exit(1);
       }
       if (strcmp(buf,"END")) {
         fprintf(stderr,"FD_read: data section must end with END and not '%s' in file '%s'\n",buf,fname);
	 /*         exit(1); */
       }
       got_data=1;
    }
  }
  fclose(fp);
  return fd;
}

char *sgets(char *str, int n, char *data) {
  int i=0;
  while (data[i] != 0 && data[i] != '\n' && i < n) {
    str[i] = data[i];
    i++;
  }
  if (i < n) str[i] = 0;
  if (i==0) return NULL;
  return &data[i+1];
}

int seof(char *data) {
  if (*data == 0) return 1;
  return 0;
}

FD* FD_readstr(char* d)
{
  int got_np=0,got_sw=0,got_data=0;
  FD* fd;
  int i,lineno,ntot,j;
  char buf[256],*p,*q;
  char com[256];
  char val[256];
  char bits[4];
  char *dd = d;

  fd=FD_alloc();
  if (!fd){
    fprintf(stderr,"FD_readstr: unable to allocate fd datatype\n");
    exit(1);
  }
  if (!dd) {
    fprintf(stderr,"FD_readstr: data pointer is null\n");
    exit(1);
  }  
  lineno=1;
  if ((dd = sgets(buf, sizeof(buf), dd)) == NULL) {
    fprintf(stderr,"FD_readstr: unable to read line %d\n",lineno);
    exit(1);
  }
  if (strncmp(buf,"SIMP",4)) {
    fprintf(stderr,"FD_readstr: invalid type of file '%s' (expected 'SIMP')\n",buf);
    exit(1);
  }
  while (!seof(dd) && !got_data) {
    lineno++;
    if ((dd = sgets(buf,sizeof(buf),dd)) == NULL) {
      if (got_sw && got_data && got_np) break;
      buf[0]=0;
      if (!got_sw) {
        strcat(buf," SW=<spectralwidth>,");
      }
      if (!got_data) {
        strcat(buf," datasection,");
      }
      if (!got_np) {
        strcat(buf," NP=<number of points>,");
      }
      fprintf(stderr,"FD_readstr: missing%s\n",buf);
      exit(1);
    }
    p=buf;
    q=p+strlen(buf)-1;
    if (*q == '\n') *q=0;
    while (isspace(*p)) p++;
    if (*p == 0) continue;
    q=(char*)strchr(p,'=');    
    if (q) {
       *q=0; q++;
       while (isspace(*p)) p++;

       strcpy(com,p);
       strcpy(val,q);

       if (!strcmp(com,"NP")) {
          fd->np=atoi(val);
          if (fd->np <= 0) {
            fprintf(stderr,"FD_readstr: number of points must be larger than zero\n");
            exit(1);
          }
          got_np=1;
       } else if (!strcmp(com,"NI")) {
          fd->ni=atoi(val);
          if (fd->ni <= 0) {
            fprintf(stderr,"FD_readstr: number of points must be larger than zero\n");
            exit(1);
          }
       } else if (!strcmp(com,"NELEM")) {
          fd->nelem=atoi(val);
          if (fd->nelem <= 0) {
            fprintf(stderr,"FD_readstr: number of elements NELEM must be larger than zero\n");
            exit(1);
          }
       } else if (!strcmp(com,"SW")) {
          if (sscanf(val,"%lg",&fd->sw) != 1) {
            fprintf(stderr,"FD_readstr: unable to convert '%s' to double\n",val);
            exit(1);
          }
          got_sw=1;
       } else if (!strcmp(com,"SW1")) {
          if (sscanf(val,"%lg",&fd->sw1) != 1) {
            fprintf(stderr,"FD_readstr: unable to convert '%s' to double\n",val);
            exit(1);
          }
       } else if (!strcmp(com,"TYPE")) {
          if (!strcmp(val,"FID")) {
            fd->type=FD_TYPE_FID;
          } else if (!strcmp(val,"SPE")) {
            fd->type=FD_TYPE_SPE;
          } else {
            fprintf(stderr,"FD_readstr: unable to convert '%s' to SPE or FID\n",val);
            exit(1);
          }
       } else if (!strcmp(com,"FORMAT")) {
          if (!strcmp(val,"TEXT")) {
            fd->format=FD_FORMAT_TEXT;
          } else if (!strcmp(val,"BINARY")) {
            fd->format=FD_FORMAT_BINARY;
          } else {
            fprintf(stderr,"FD_readstr: unable to convert '%s' to TEXT or BINARY\n",val);
            exit(1);
          }
       } else if (!strcmp(com,"PREC")) {
          if (!strcmp(val,"DOUBLE")) {
            fd->prec=FD_PREC_DOUBLE;
          } else if (!strcmp(val,"SINGLE")) {
            fd->prec=FD_PREC_SINGLE;
          } else {
            fprintf(stderr,"FD_readstr: unable to convert '%s' to DOUBLE or SINGLE\n",val);
            exit(1);
          }
       } else if (!strcmp(com,"REF")) {
          if (sscanf(val,"%lg",&fd->ref) != 1) {
            fprintf(stderr,"FD_readstr: unable to convert '%s' to double\n",val);
            exit(1);
          }
       } else if (!strcmp(com,"REF1")) {
          if (sscanf(val,"%lg",&fd->ref1) != 1) {
            fprintf(stderr,"FD_readstr: unable to convert '%s' to double\n",val);
            exit(1);
          }
       }
    } else {
       if (strcmp(buf,"DATA")) {
         fprintf(stderr,"FD_readstr: data section must start with DATA and not '%s'\n",buf);
         exit(1);
       }
       FD_alloc1ddata(fd);
       ntot=fd->np*(fd->ni > 1 ? fd->ni : 1)*fd->nelem;

       if (fd->format == FD_FORMAT_TEXT) {
         for (i=1;i<=ntot;i++) {
           if (sscanf(dd,"%lg%lg%n",&(fd->data[i].re),&(fd->data[i].im),&j) != 2) {
              fprintf(stderr,"FD_read: unable to read data point number %d\n",i);
              exit(1);
           }
	   dd += j;
	   if (*dd == '\n') dd++;
         }
       } else {
         if (fd->prec == FD_PREC_SINGLE) {
           int i,doubles,c;
           double* data;

           doubles=ntot*2;
           data = &(fd->data)[1].re;
           _pack_nbuf=0;
	   _pack_newline=0;
           for (i=0;i<doubles;i++) {
             for (j=0;j<4;j++) {
               if (_pack_nbuf==0) {
	         unsigned char c0,c1,c2,c3;
		 if (*dd == '\n') dd++;
		 c0=(unsigned char)*dd-BASE;dd++;
		 c1=(unsigned char)*dd-BASE;dd++;
		 c2=(unsigned char)*dd-BASE;dd++;
		 c3=(unsigned char)*dd-BASE;dd++;
                 _pack_buf[0]=FIRST(6,c0) | LAST(2,c1 << 2);
                 _pack_buf[1]=FIRST(4,c1) | LAST(4,c2 << 2);
                 _pack_buf[2]=FIRST(2,c2) | LAST(6,c3 << 2);
                 _pack_nbuf=3;
	       }
	       c = (int)_pack_buf[3-_pack_nbuf--];
               bits[j]=c;
             }
             data[i]=bits_to_float(bits);
           }
           _pack_nbuf=-1;  
	   while (*dd == '\n') dd++;
         } else {
           int nb;
	   if (strlen(dd) < ntot*sizeof(complx)) {
             fprintf(stderr,"FD_read: read %d but needed %d complex numbers\n",nb,ntot);
             exit(1);           
           } else {
	     memcpy(&fd->data[1], dd, ntot*sizeof(complx));
	   }
         }
       }
       strcpy(buf,"");
       if ((dd = sgets(buf,sizeof(buf),dd)) == NULL) {
         fprintf(stderr,"FD_readstr: must read a string (\"END\") at the end of data\n");
         exit(1);
       }
       if (strcmp(buf,"END")) {
         fprintf(stderr,"FD_readstr: data section must end with END and not '%s'\n",buf);
	 exit(1);
       }
       got_data=1;
    }
  }
  return fd;
}

FD* FD_dup(FD* fd)
{
  FD* to;

  if (!fd) {
    fprintf(stderr,"FD_dup: called with NULL argument\n");
    exit(1);
  }
  to=FD_alloc();
  to->np=fd->np;
  to->ni=fd->ni;
  to->nelem=fd->nelem;
  to->type=fd->type;
  to->format=fd->format;
  to->prec=fd->prec;
  to->ref=fd->ref;
  to->ref1=fd->ref1;
  to->sw=fd->sw;
  to->sw1=fd->sw1;
  FD_alloc1ddata(to);
  memcpy(to->data,fd->data,sizeof(complx)*(to->np*(to->ni > 1 ? to->ni : 1)*to->nelem+1));
  return to;
}

FD* FD_data2fd(char* file,complx* data,int np,int ni,double sw,double sw1, int nelem)
{
  FD* fd;

  fd=FD_alloc();
  fd->np=np;
  fd->ni=ni;
  fd->nelem=nelem;
  fd->sw=sw;
  fd->sw1=sw1;
  FD_alloc1ddata(fd);
  memcpy(fd->data,data,sizeof(complx)*(fd->np*(fd->ni > 1 ? fd->ni : 1)*fd->nelem+1));
  return fd;
}




/* Routines inserted by Thomas Vosegaard, 2001 */

void mac2unix(char *in)
{
  int i = 0;
  
  while (in[i]) {
    if (in[i] == 13) in[i] = 10;
    i++;
  }
}

void unix2mac(char *in)
{
  int i = 0;
  
  while (in[i]) {
    if (in[i] == 10) in[i] = 13;
    i++;
  }
}

void swap(char *ptr, int size)
{
  int i;
  char *p, *rp, swp;

  for (i = 0; i < size/2; i++) {
    p = ptr + i;
    rp = ptr + (size-1) - i;
    swp = *p;
    *p = *rp;
    *rp = swp;
  }
}


size_t freadreversedbyteorder(void *ptr, size_t itemsize, size_t nitem, FILE *file)
{
  int i, tmp;
  char *cptr = ptr;
  tmp = fread(ptr, itemsize, nitem, file);

  if (itemsize > 1)
    for (i = 0; i < nitem; i++) swap(cptr+i*itemsize, itemsize);

  return (size_t)tmp;
}

size_t fwritereversedbyteorder(const void *ptr, size_t itemsize, size_t nitem, FILE *file)
{
  int i, tmp;
  char *p, *cptr = (char *)ptr;
  
  if (itemsize == 1) {
    tmp = fwrite(ptr, 1, nitem, file);
  } else {
    p = malloc(itemsize*nitem);
    for (i = 0; i < nitem*itemsize; i++) p[i] = cptr[i];
    for (i = 0; i < nitem; i++) swap((char *)(p+i*itemsize), itemsize);
    tmp = fwrite(p, itemsize, nitem, file);
    free(p);
  }

  return (size_t)tmp;
}


#ifdef REVERSEBYTES
#  define freadreversed freadreversedbyteorder
#  define fwritereversed fwritereversedbyteorder
#else
#  define freadreversed fread
#  define fwritereversed fwrite
#endif

void FD_write_rmn(char *fname, FD* fd)
{
  int i,ntot;
  char version, title[512],c;
  double dd;
  float df;
  FILE* fp;

  if (fd->nelem != 1) {
	  fprintf(stderr,"Error: FD_write_rmn detected multiple elements in fid/spc - feature not supported here\n");
	  exit(1);
  }
  fp=fopen(fname,"w");
  if (!fp) {
    fprintf(stderr,"FD_write_rmn: unable to create file '%s'\n",fname);
    exit(1);
  }
  if (fd->ni == 0) {
    /*
    if (fd->type == FD_TYPE_FID) {
      printf("Saving 1D fid\n");
    } else {
      printf("Saving 1D spec\n");
    }
    */
    version=2;
  } else {
    /*
    if (fd->type == FD_TYPE_FID) {
      printf("Saving 2D fid\n");
    } else {
      printf("Saving 2D spec\n");
    }
    */
    version=4;
  }
  
  if (!fwritereversed(&version,1,1,fp)) {
    fprintf(stderr,"FD_write_rmn: unable to write version"
                   " to file '%s'\n",fname);
    exit(1);
  }

  if (fd->np <= 0) {
    fprintf(stderr,"FD_write_rmn: number of points with must be larger"
                    " than zero when writing file '%s'\n",fname);
    exit(1);
  }
  i = fd->np;
  if (!fwritereversed(&i, sizeof(int), 1, fp)) {
    fprintf(stderr,"FD_write_rmn: unable to write number of points"
                   " to file '%s'\n",fname);
    exit(1);
  }

  if (fd->sw <= 0.0) {
    fprintf(stderr,"FD_write: spectral-width must be larger"
                   " than zero when writing file '%s'\n",fname);
    exit(1);
  }
  dd = 1.0/fd->sw;
  if (!fwritereversed(&dd, sizeof(double), 1, fp)) {
    fprintf(stderr,"FD_write_rmn: unable to write dwell-time"
                   " to file '%s'\n",fname);
    exit(1);
  }

  dd = 0;
  if (!fwritereversed(&dd, sizeof(double), 1, fp)) {
    fprintf(stderr,"FD_write_rmn: unable to write starting time"
                   " to file '%s'\n",fname);
    exit(1);
  }
  
  if (!fwritereversed(&dd, sizeof(double), 1, fp)) {
    fprintf(stderr,"FD_write_rmn: unable to write larmor frequency"
                   " to file '%s'\n",fname);
    exit(1);
  }

  dd = fd->ref;
  if (!fwritereversed(&dd, sizeof(double), 1, fp)) {
    fprintf(stderr,"FD_write_rmn: unable to write reference point"
                   " to file '%s'\n",fname);
    exit(1);
  }
  sprintf(title,fname);
  unix2mac(title);
  ntot=fd->np*(fd->ni > 1 ? fd->ni : 1);

  switch((int)version) {
  case 2:
    if (!fwritereversed(&title, 1, 512, fp)) {
      fprintf(stderr,"FD_write_rmn: unable to write file name"
                     " to file '%s'\n",fname);
      exit(1);
    }
    c = fd->type == FD_TYPE_FID ? 'f':'s';
    if (!fwritereversed(&c, 1, 1, fp)) {
      fprintf(stderr,"FD_write_rmn: unable to write spectrum type"
                     " to file '%s'\n",fname);
      exit(1);
    }
    c = 0;
    if (!fwritereversed(&title, 1, 511, fp)) {
      fprintf(stderr,"FD_write_rmn: unable to write 511 dummy bytes"
                     " to file '%s'\n",fname);
      exit(1);
    }
    break;
  case 4:
    i = fd->ni;
    if (!fwritereversed(&i, sizeof(int), 1, fp)) {
      fprintf(stderr,"FD_write_rmn: unable to write number of increments"
                     " to file '%s'\n",fname);
      exit(1);
    }

    if (fd->sw1 <= 0.0) {
      fprintf(stderr,"FD_write_rmn: indirect spectral-width must be larger"
                     " than zero when writing file '%s'\n",fname);
      exit(1);
    }
    dd = 1.0/fd->sw1;
    if (!fwritereversed(&dd, sizeof(double), 1, fp)) {
      fprintf(stderr,"FD_write_rmn: unable to write indirect dwell-time"
                     " to file '%s'\n",fname);
      exit(1);
    }

    dd = 0; 	
    if (!fwritereversed(&dd, sizeof(double), 1, fp)) {
      fprintf(stderr,"FD_write_rmn: unable to write indirect starting time"
                     " to file '%s'\n",fname);
      exit(1);
    }
  
    if (!fwritereversed(&dd, sizeof(double), 1, fp)) {
      fprintf(stderr,"FD_write_rmn: unable to write indirect larmor frequency"
                     " to file '%s'\n",fname);
      exit(1);
    }

    dd = fd->ref1;
    if (!fwritereversed(&dd, sizeof(double), 1, fp)) {
      fprintf(stderr,"FD_write_rmn: unable to write indirect reference point"
                     " to file '%s'\n",fname);
      exit(1);
    }
    if (!fwritereversed(&title, 1, 512, fp)) {
      fprintf(stderr,"FD_write_rmn: unable to write file name"
                     " to file '%s'\n",fname);
      exit(1);
    }
    c = fd->type == FD_TYPE_FID ? 'f':'s';
    if (!fwritereversed(&c, 1, 1, fp)) {
      fprintf(stderr,"FD_write_rmn: unable to write spectrum type"
                     " to file '%s'\n",fname);
      exit(1);
    }
    if (!fwritereversed(&c, 1, 1, fp)) {
      fprintf(stderr,"FD_write_rmn: unable to write indirect spectrum type"
                     " to file '%s'\n",fname);
      exit(1);
    }
    if (!fwritereversed(&title, 1, 510, fp)) {
      fprintf(stderr,"FD_write_rmn: unable to write 510 dummy bytes"
                     " to file '%s'\n",fname);
      exit(1);
    }
    break;
  }

  for (i=1; i<= ntot; i++) {
    df = fd->data[i].re;
    if (!fwritereversed(&df, sizeof(float), 1, fp)) {
      fprintf(stderr,"FD_write_rmn: unable to write real data point %d"
                     " to file '%s'\n",i,fname);
      exit(1);
    }
    df = fd->data[i].im;
    if (!fwritereversed(&df, sizeof(float), 1, fp)) {
      fprintf(stderr,"FD_write_rmn: unable to write imaginary data point %d"
                     " to file '%s'\n",i,fname);
      exit(1);
    }
    if ((i % fd->np) == 0) {
      df = fd->data[i-fd->np+1].re;
      if (!fwritereversed(&df, sizeof(float), 1, fp)) {
        fprintf(stderr,"FD_write_rmn: unable to write real data point **"
                       " to file '%s'\n",fname);
        exit(1);
      }
      df = fd->data[i-fd->np+1].im;
      if (!fwritereversed(&df, sizeof(float), 1, fp)) {
        fprintf(stderr,"FD_write_rmn: unable to write imaginary data point **"
                       " to file '%s'\n",fname);
        exit(1);
      }
    }
  }
  fclose(fp);
}







/*
  Varian's file format


  FORMAT (1D):
  
  data file header
  header for block 1
  data of block 1
  header for block 2
  data of block 2
  ...

  FORMAT (Hypercomplex 2D) - datafile.nbheaders=2:
  
  data file header
  header for block 1
  second header for block 1
  data of block 1
  header for block 2
  second header for block 2
  data of block 2
  ...

*/
struct VNMR_datafilehead {
  long nblocks;		/* number of blocks in file */
  long ntraces;		/* number of traces per block */
  long np; 		/* number of elements per trace */
  long ebytes;		/* number of bytes per elements */
  long tbytes;		/* number of bytes per trace */
  long bbytes;		/* number of bytes per block */
  short vers_id;	/* software version, file_id status bits */
  short status;		/* status of whole file */
  long nbheaders;	/* number of block headers per block */
};

#define S_DATA		(1<<0)
#define S_SPEC		(1<<1)
#define S_32		(1<<2)
#define S_FLOAT		(1<<3)
#define S_COMPLEX	(1<<4)
#define S_HYPERCOMPLEX	(1<<5)
#define S_ACQPAR	(1<<7)
#define S_SECND		(1<<8)
#define S_TRANSF	(1<<9)
#define S_NP		(1<<11)
#define S_NF		(1<<12)
#define S_NI		(1<<13)
#define S_NI2		(1<<14)

int freadreversed_datafilehead(struct VNMR_datafilehead *h, FILE *f)
{
  if (freadreversed(&h->nblocks,   sizeof(long), 1, f) < 0) return -1;
  if (freadreversed(&h->ntraces,   sizeof(long), 1, f) < 0) return -1;
  if (freadreversed(&h->np,        sizeof(long), 1, f) < 0) return -1;
  if (freadreversed(&h->ebytes,    sizeof(long), 1, f) < 0) return -1;
  if (freadreversed(&h->tbytes,    sizeof(long), 1, f) < 0) return -1;
  if (freadreversed(&h->bbytes,    sizeof(long), 1, f) < 0) return -1;
  if (freadreversed(&h->vers_id,   sizeof(short),1, f) < 0) return -1;
  if (freadreversed(&h->status,    sizeof(short),1, f) < 0) return -1;
  if (freadreversed(&h->nbheaders, sizeof(long), 1, f) < 0) return -1;
  return 1;
}

int fwritereversed_datafilehead(struct VNMR_datafilehead *h, FILE *f)
{
  if (fwritereversed(&h->nblocks,   sizeof(long), 1, f) < 0) return -1;
  if (fwritereversed(&h->ntraces,   sizeof(long), 1, f) < 0) return -1;
  if (fwritereversed(&h->np,        sizeof(long), 1, f) < 0) return -1;
  if (fwritereversed(&h->ebytes,    sizeof(long), 1, f) < 0) return -1;
  if (fwritereversed(&h->tbytes,    sizeof(long), 1, f) < 0) return -1;
  if (fwritereversed(&h->bbytes,    sizeof(long), 1, f) < 0) return -1;
  if (fwritereversed(&h->vers_id,   sizeof(short),1, f) < 0) return -1;
  if (fwritereversed(&h->status,    sizeof(short),1, f) < 0) return -1;
  if (fwritereversed(&h->nbheaders, sizeof(long), 1, f) < 0) return -1;
  return 1;
}

int print_datafilehead(struct VNMR_datafilehead *h)
{
  int i;
  if (!h) return -1;

  printf("Number of blocks:                  %6ld\n", h->nblocks);  
  printf("Number of traces:                  %6ld\n", h->ntraces);  
  printf("Number of elements per trace:      %6ld\n", h->np);  
  printf("Number of bytes per element:       %6ld\n", h->ebytes);  
  printf("Number of bytes per trace:         %6ld\n", h->tbytes);  
  printf("Number of bytes per block:         %6ld\n", h->bbytes);  
  printf("Version ID:                        %6d\n", (int)h->vers_id);  
  printf("File status: (%5d)     ", h->status);
  for (i=sizeof(h->status)*8-1; i>=0; i--) 
    printf("%d", (h->status & (1<<i))?1:0);
  printf("\n");
  printf("Number of block headers per block: %6ld\n", h->nbheaders);  

  return 1;
}

struct VNMR_datablockhead {
  short scale;		/* scaling factor */
  short status;		/* status of data in block */
  short index;		/* block index */
  short mode;		/* mode of data in block */
  long ctcount;		/* ct value for FID */
  float lpvalue;	/* f2 (2D-f1) left phase in phasefile */
  float rpvalue;	/* f2 (2D-f1) right phase in phasefile */
  float lvl;		/* level drift correction */
  float tlt;		/* tilt drift correction */
};

#define MORE_BLOCKS	(1<<7)
#define NP_CMPLX	(1<<8)
#define NF_CMPLX	(1<<9)
#define NI_CMPLX	(1<<10)
#define NI2_CMPLX	(1<<11)

#define NP_PHMODE	(1<<0)
#define NP_AVMODE	(1<<1)
#define NP_PWRMODE	(1<<2)
#define NF_PHMODE	(1<<4)
#define NF_AVMODE	(1<<5)
#define NF_PWRMODE	(1<<6)
#define NI_PHMODE	(1<<8)
#define NI_AVMODE	(1<<9)
#define NI_PWRMODE	(1<<10)
#define NI2_PHMODE	(1<<12)
#define NI2_AVMODE	(1<<13)
#define NI2_PWRMODE	(1<<14)

int freadreversed_datablockhead(struct VNMR_datablockhead *h, FILE *f)
{
  if (freadreversed(&h->scale,   sizeof(short), 1, f) < 0) return -1;
  if (freadreversed(&h->status,  sizeof(short), 1, f) < 0) return -1;
  if (freadreversed(&h->index,   sizeof(short), 1, f) < 0) return -1;
  if (freadreversed(&h->mode,    sizeof(short), 1, f) < 0) return -1;
  if (freadreversed(&h->ctcount, sizeof(long),  1, f) < 0) return -1;
  if (freadreversed(&h->lpvalue, sizeof(float), 1, f) < 0) return -1;
  if (freadreversed(&h->rpvalue, sizeof(float), 1, f) < 0) return -1;
  if (freadreversed(&h->lvl,     sizeof(float), 1, f) < 0) return -1;
  if (freadreversed(&h->tlt,     sizeof(float), 1, f) < 0) return -1;
  return 1;
}

int fwritereversed_datablockhead(struct VNMR_datablockhead *h, FILE *f)
{
  if (fwritereversed(&h->scale,   sizeof(short), 1, f) < 0) return -1;
  if (fwritereversed(&h->status,  sizeof(short), 1, f) < 0) return -1;
  if (fwritereversed(&h->index,   sizeof(short), 1, f) < 0) return -1;
  if (fwritereversed(&h->mode,    sizeof(short), 1, f) < 0) return -1;
  if (fwritereversed(&h->ctcount, sizeof(long),  1, f) < 0) return -1;
  if (fwritereversed(&h->lpvalue, sizeof(float), 1, f) < 0) return -1;
  if (fwritereversed(&h->rpvalue, sizeof(float), 1, f) < 0) return -1;
  if (fwritereversed(&h->lvl,     sizeof(float), 1, f) < 0) return -1;
  if (fwritereversed(&h->tlt,     sizeof(float), 1, f) < 0) return -1;
  return 1;
}

int print_datablockhead(struct VNMR_datablockhead *h)
{
  int i;
  if (!h) return -1;

  printf("Scaling factor:                    %6d\n", (int)h->scale);  
  printf("Block status: (%5d)    ", h->status);
  for (i=sizeof(h->status)*8-1; i>=0; i--) 
    printf("%d", (h->status & (1<<i))?1:0);
  printf("\n");
  printf("Block index:                       %6d\n", (int)h->index);  
  printf("Data mode:                         %6d\n", (int)h->mode);  
  printf("ct value:                          %6ld\n", h->ctcount);  
  printf("lp value:                          %6.1f\n", h->lpvalue);  
  printf("rp value:                          %6.1f\n", h->rpvalue);  
  printf("Level drift correction:            %6.1f\n", h->lvl);
  printf("Tilt drift correction:             %6.1f\n", h->tlt);

  return 1;
}

struct VNMR_hypercmplxbhead {
  short s_spare1;	/* short word: spare */
  short status;		/* status word for block header */
  short s_spare2;	/* short word: spare */
  short s_spare3;	/* short word: spare */
  long l_spare1;	/* long word: spare */
  float lpval1;		/* 2D-f2 left phase */
  float rpval1;		/* 2D-f2 right phase */
  float f_spare1;	/* float word: spare */
  float f_spare2;	/* float word: spare */
};

#define U_HYPERCOMPLEX (1<<1)

int freadreversed_hypercmplxbhead(struct VNMR_hypercmplxbhead *h, FILE *f)
{
  if (freadreversed(&h->s_spare1, sizeof(short), 1, f) < 0) return -1;
  if (freadreversed(&h->status,   sizeof(short), 1, f) < 0) return -1;
  if (freadreversed(&h->s_spare2, sizeof(short), 1, f) < 0) return -1;
  if (freadreversed(&h->s_spare3, sizeof(short), 1, f) < 0) return -1;
  if (freadreversed(&h->l_spare1, sizeof(long),  1, f) < 0) return -1;
  if (freadreversed(&h->lpval1,   sizeof(float), 1, f) < 0) return -1;
  if (freadreversed(&h->rpval1,   sizeof(float), 1, f) < 0) return -1;
  if (freadreversed(&h->f_spare1, sizeof(float), 1, f) < 0) return -1;
  if (freadreversed(&h->f_spare2, sizeof(float), 1, f) < 0) return -1;
  return 1;
}

int fwritereversed_hypercmplxbhead(struct VNMR_hypercmplxbhead *h, FILE *f)
{
  if (fwritereversed(&h->s_spare1, sizeof(short), 1, f) < 0) return -1;
  if (fwritereversed(&h->status,   sizeof(short), 1, f) < 0) return -1;
  if (fwritereversed(&h->s_spare2, sizeof(short), 1, f) < 0) return -1;
  if (fwritereversed(&h->s_spare3, sizeof(short), 1, f) < 0) return -1;
  if (fwritereversed(&h->l_spare1, sizeof(long),  1, f) < 0) return -1;
  if (fwritereversed(&h->lpval1,   sizeof(float), 1, f) < 0) return -1;
  if (fwritereversed(&h->rpval1,   sizeof(float), 1, f) < 0) return -1;
  if (fwritereversed(&h->f_spare1, sizeof(float), 1, f) < 0) return -1;
  if (fwritereversed(&h->f_spare2, sizeof(float), 1, f) < 0) return -1;
  return 1;
}

int print_hypercmplxbhead(struct VNMR_hypercmplxbhead *h)
{
  int i;
  if (!h) return -1;

  printf("Short spare 1:       %6d\n", (int)h->s_spare1);  
  printf("Block header status: ");
  for (i=sizeof(h->status)*8-1; i>=0; i--) 
    printf("%d", (h->status & (1<<i))?1:0);
  printf("\n");
  printf("Short spare 2:       %6d\n", (int)h->s_spare2);  
  printf("Short spare 3:       %6d\n", (int)h->s_spare3);  
  printf("Long spare 1:        %6ld\n", h->l_spare1);  
  printf("lp1 value:           %6.1f\n", h->lpval1);  
  printf("rp1 value:           %6.1f\n", h->rpval1);  
  printf("Float spare 1:       %6.1f\n", h->f_spare1);
  printf("Float spare 2:       %6.1f\n", h->f_spare2);

  return 1;
}

float getvnmrval(FILE *f, char *par)
{
  int i;
  float fl = 0;
  char buf[512], p[32];
  
  rewind(f);
  while (fgets(buf, sizeof(buf), f)) {
    i = 0;
    if (isdigit(*buf)) continue;
    while (buf[i] != ' ') p[i] = buf[i++];
    p[i] = 0;
    if (!strcmp(p, par)) {
      fgets(buf, sizeof(buf), f);
      i=0;
      while (buf[i] != ' ') i++;
      sscanf(&buf[i+1], "%g", &fl);
      break;
    }
  }
  return fl;
}
  

void FD_read_vnmr(char *fname, int swap)
{
	printf("Error: FD_read_vnmr not implemented, nothing done.\n");
}



void FD_write_nmrpipe(char* fname,FD* fd, int phsens)
{
#ifndef __MINGW32__
  int i,ni,j, no, ii, jj;
  FILE* fp;
  char tmp[256], cmd[1024];
  float f;
  double ref, ref1;
  char *aq2d;
  char *yMODE;

  if (fd->nelem != 1) {
	  fprintf(stderr,"Error: FD_write_nmrpipe detected multiple fid/spc elements - feature not supported here\n");
	  exit(1);
  }

#ifdef REVERSEBYTES
  char *swap = "-no";
#else
  char *swap = "-no";
#endif

  switch (phsens) {
    case 1:
      aq2d = "States";
      yMODE = "Complex";
      break;
    case 2:
      aq2d = "TPPI";
      yMODE = "Real";
      break;
    case 3:
      aq2d = "States";
      yMODE = "Rance-Kay";
      break;
    default:
      aq2d = "Magnitude";
      yMODE = "Real";
      break;
  }

  sprintf(tmp, "/tmp/%s-XXXXXX", getenv("USER"));
  if ((fp=fdopen(mkstemp(tmp), "w")) == NULL) {
    fprintf(stderr,"FD_write_nmrpipe: unable to connect to 'bin2pipe'\n");
    exit(1);
  }
  if (fd->ni <= 1) {
    ref = fd->ref / (fd->sfrq == 0 ? 1:fd->sfrq);
    sprintf(cmd, "bin2pipe -ndim 1 -in %s -out %s -ov %sswap"
                 " -xSW %g"
		 " -xMODE Complex"
		 " -xN %d"
		 " -xOBS %g"
		 " -xCAR %g"
		 " -xLAB \"Obs\""
		 " -xFT %s",
		 tmp, fname, swap, fd->sw, fd->np*2, fd->sfrq,
		 ref,
		 (fd->type == FD_TYPE_FID) ? "Time":"Freq");
  } else {
    ref = fd->ref / (fd->sfrq == 0 ? 1:fd->sfrq);
    ref1 = fd->ref1 / (fd->sfrq1 == 0 ? 1:fd->sfrq1);
    sprintf(cmd, "bin2pipe -ndim 2 -in %s -out %s -ov %sswap -aq2D %s"
                 " -xSW %g"
		 " -xMODE Complex"
		 " -xN  %d"
		 " -xOBS %g"
		 " -xCAR %g"
		 " -xORIG %g"
		 " -xLAB \"Direct\""
		 " -xFT %s"
                 " -ySW %g"
		 " -yMODE %s"
		 " -yN  %d"
		 " -yOBS %g"
		 " -yCAR %g"
		 " -yORIG %g"
		 " -yLAB \"Indirect\""
		 " -yFT %s",
		 tmp, fname, swap, aq2d,
		 fd->sw, fd->np*2, fd->sfrq, ref, 
		 -fd->sw/2.0,
		 (fd->type == FD_TYPE_FID) ? "Time":"Freq",
                 fd->sw1, yMODE,
		 fd->ni, fd->sfrq1, ref1,
		 -fd->sw1/2.0+fd->sw1/(fd->ni*2),
		 (fd->type == FD_TYPE_FID) ? "Time":"Freq");
  }
  printf("%s\n", cmd);

  ni=(fd->ni > 1 ? fd->ni : 1);

  if (fd->type == FD_TYPE_SPE) {
    if (ni > 1) {
      for (j=ni+1; j>1; j--) { 
        if (j > ni) jj = 1;
	else jj = j;
        no = (jj-1)*fd->np;
        for (i=no+fd->np+1;i>no+1;i--) {
	  if (i > no+fd->np) ii = 1;
	  else ii = i;
          f = fd->data[ii].re;
          fwrite(&f, sizeof(float), 1, fp);
        }
        for (i=no+fd->np;i>no;i--) {
	  if (i > no+fd->np) ii = 1;
	  else ii = i;
          f = fd->data[ii].im;
          fwrite(&f, sizeof(float), 1, fp);
        }
      }
    } else {
      for (i=fd->np;i>=1;i--) {
        f = fd->data[i].re;
        fwrite(&f, sizeof(float), 1, fp);
      }
      for (i=fd->np;i>=1;i--) {
        f = fd->data[i].im;
        fwrite(&f, sizeof(float), 1, fp);
      }
    }
  } else {
    if (ni > 1) {
      for (j=1; j<= ni; j++) { 
        no = (j-1)*fd->np+1;
        for (i=no;i<no+fd->np;i++) {
          f = fd->data[i].re;
          fwrite(&f, sizeof(float), 1, fp);
        }
        for (i=no;i<no+fd->np;i++) {
          f = fd->data[i].im;
          fwrite(&f, sizeof(float), 1, fp);
        }
      }
    } else {
      for (i=1;i<=fd->np;i++) {
        f = fd->data[i].re;
        fwrite(&f, sizeof(float), 1, fp);
      }
      for (i=1;i<=fd->np;i++) {
        f = fd->data[i].im;
        fwrite(&f, sizeof(float), 1, fp);
      }
    }
  }
  fclose(fp);
  system(cmd);
  sprintf(cmd, "rm -f %s", tmp);
  system(cmd);
#endif
}

/****
 * ZT: write in raw binary format, without any header and trailer
 ****/
void FD_write_raw_bin(char* fname,FD* fd)
{
  int i,ni,j, no, ii, jj;
  FILE* fp;
  float f;

  if (fd->nelem != 1) {
	  fprintf(stderr,"Error: FD_write_raw_bin detected multiple fid/spc elements - feature not supported here\n");
	  exit(1);
  }

  if ((fp=fopen(fname, "w")) == NULL) {
    fprintf(stderr,"FD_write_raw_bin: unable to create file %s\n",fname);
    exit(1);
  }

  ni=(fd->ni > 1 ? fd->ni : 1);

  if (fd->type == FD_TYPE_SPE) {
    if (ni > 1) {
      for (j=ni+1; j>1; j--) { 
        if (j > ni) jj = 1;
	else jj = j;
        no = (jj-1)*fd->np;
        for (i=no+fd->np+1;i>no+1;i--) {
	  if (i > no+fd->np) ii = 1;
	  else ii = i;
          f = fd->data[ii].re;
          fwrite(&f, sizeof(float), 1, fp);
        }
        for (i=no+fd->np;i>no;i--) {
	  if (i > no+fd->np) ii = 1;
	  else ii = i;
          f = fd->data[ii].im;
          fwrite(&f, sizeof(float), 1, fp);
        }
      }
    } else {
      for (i=fd->np;i>=1;i--) {
        f = fd->data[i].re;
        fwrite(&f, sizeof(float), 1, fp);
      }
      for (i=fd->np;i>=1;i--) {
        f = fd->data[i].im;
        fwrite(&f, sizeof(float), 1, fp);
      }
    }
  } else {
    if (ni > 1) {
      for (j=1; j<= ni; j++) { 
        no = (j-1)*fd->np+1;
        for (i=no;i<no+fd->np;i++) {
          f = fd->data[i].re;
          fwrite(&f, sizeof(float), 1, fp);
        }
        for (i=no;i<no+fd->np;i++) {
          f = fd->data[i].im;
          fwrite(&f, sizeof(float), 1, fp);
        }
      }
    } else {
      for (i=1;i<=fd->np;i++) {
        f = fd->data[i].re;
        fwrite(&f, sizeof(float), 1, fp);
      }
      for (i=1;i<=fd->np;i++) {
        f = fd->data[i].im;
        fwrite(&f, sizeof(float), 1, fp);
      }
    }
  }
  fclose(fp);

}



double getpipeval(char *fname, char *flag, int dim)
{
  FILE* fp;
  double value;
  char buf[256];
  
  if (dim) {
    sprintf(buf, "getParm -in %s -parm %s -dim %d", fname, flag, dim);
  } else {
    sprintf(buf, "getParm -in %s -parm %s", fname, flag);
  }

  fp=popen(buf,"r");
  if (!fp) {
    fprintf(stderr,"FD_read: unable to open file '%s'\n",fname);
    exit(1);
  }

  fgets(buf, 256, fp);
  if (sscanf(buf, "%lf", &value) != 1) {
    fprintf(stderr,"getParm: unable to retrieve parameter '%s'\n",
    flag);
    exit(1);
  }
  pclose(fp);

  return value;
}



/***********************************************
  testHdr, swapHdr, isNullHdr, isHdrStr
  are part of the NMRPipe processing
  system. Used here with permission
  from F. Delaglio.
************************************************/

int testHdr( fdata )

   float  fdata[FDATASIZE];
{
    int   status;

    status = HDR_OK;

#ifndef CRAY

    if (fdata[FDMAGIC] != 0.0)
       {
        return( HDR_BAD );
       }

    if ((fdata[FDFLTORDER] != 0.0) && (fdata[FDFLTORDER] != (float)FDORDERCONS))
       {
        (void) swapHdr( fdata );

        if (fdata[FDFLTORDER] != (float)FDORDERCONS)
           {
            return( HDR_BAD );
           }

        status = HDR_SWAPPED;
       }

    if (fdata[FDSIZE] <= 0.0 || fdata[FDSIZE] != (int)fdata[FDSIZE])
       {
        return( HDR_BAD );
       }
#endif

    return( status );
}

int swapHdr( fdata )

   float fdata[FDATASIZE];
{
    int   i;
    union fsw { float f; char s[4]; } in, out;

    for( i = 0; i < FDATASIZE; i++ )
       {
        if (!isHdrStr( i ))
           {
            in.f = fdata[i];

            out.s[0] = in.s[3]; 
            out.s[1] = in.s[2]; 
            out.s[2] = in.s[1]; 
            out.s[3] = in.s[0]; 

            fdata[i] = out.f;
           }
       }

    return( 0 );
}

int isNullHdr( fdata )

   float *fdata;
{
    int i, count;

    count = 0;

    for( i = 0; i < FDATASIZE; i++ ) if (fdata[i] == 0.0) count++;

    if (count == FDATASIZE) return( 1 );

    return( 0 );
}


#define FLTSIZE 4
int isHdrStr( parmLoc )

   int parmLoc;
{
    if (parmLoc >= FDF2LABEL &&
          parmLoc < FDF2LABEL + SIZE_NDLABEL/FLTSIZE) return( 1 );

    if (parmLoc >= FDF1LABEL &&
          parmLoc < FDF1LABEL + SIZE_NDLABEL/FLTSIZE) return( 1 );

    if (parmLoc >= FDF3LABEL &&
          parmLoc < FDF3LABEL + SIZE_NDLABEL/FLTSIZE) return( 1 );

    if (parmLoc >= FDF4LABEL &&
          parmLoc < FDF4LABEL + SIZE_NDLABEL/FLTSIZE) return( 1 );

    if (parmLoc >= FDSRCNAME &&
          parmLoc < FDSRCNAME + SIZE_SRCNAME/FLTSIZE) return( 1 );

    if (parmLoc >= FDUSERNAME &&
          parmLoc < FDUSERNAME + SIZE_USERNAME/FLTSIZE) return( 1 );

    if (parmLoc >= FDOPERNAME &&
          parmLoc < FDOPERNAME + SIZE_OPERNAME/FLTSIZE) return( 1 );

    if (parmLoc >= FDTITLE  &&
          parmLoc < FDTITLE + SIZE_TITLE/FLTSIZE) return( 1 );

    if (parmLoc >= FDCOMMENT &&
          parmLoc < FDCOMMENT + SIZE_COMMENT/FLTSIZE) return( 1 );

    if (parmLoc == NDLABEL)  return( 1 );

    if (parmLoc == NDLABEL1) return( 1 );

    if (parmLoc == NDLABEL2) return( 1 );

    return( 0 );
}

/***************************************************/



FD* FD_read_nmrpipe(char* fname)
{
  FILE* fp;
  FD* fd;
  int i,ndim, type,j,no,swap, ii, jj;
  float f, head[FDATASIZE];
  size_t (*freadpipe)(void *, size_t, size_t, FILE*);

  fd=FD_alloc();
  if (!fd){
    fprintf(stderr,"FD_read_nmrpipe: unable to allocate fd datatype\n");
    exit(1);
  }
  if (!fname) {
    fprintf(stderr,"FD_read_nmrpipe: filename pointer is null\n");
    exit(1);
  }

  ndim = (int)getpipeval(fname, "FDDIMCOUNT", 0);
  if (ndim != 1 && ndim != 2) {
    fprintf(stderr, "FD_read_nmrpipe: only support for 1 and 2 dimensional spectra\n");
    exit(1);
  }
  
  fd->np = (int)getpipeval(fname, "NDSIZE", 1);
  fd->sw = getpipeval(fname, "NDSW", 1);
  fd->ref = getpipeval(fname, "NDCAR", 1);
  fd->sfrq = getpipeval(fname, "NDOBS", 1);
  if (fd->sfrq != 0) fd->ref *= fd->sfrq;
  type = (int)getpipeval(fname, "NDFTFLAG", 1);

  if (ndim == 2) {
    fd->ni = (int)getpipeval(fname, "NDSIZE", 2);
    fd->sw1 = getpipeval(fname, "NDSW", 2);
    fd->ref1 = getpipeval(fname, "NDCAR", 2);
    fd->sfrq1 = getpipeval(fname, "NDOBS", 2);
    if (fd->sfrq1 != 0) fd->ref1 *= fd->sfrq1;
    if (type != (int)getpipeval(fname, "NDFTFLAG", 2)) {
      fprintf(stderr, "FD_read_nmrpipe: all dimensions must time or frequency\n");
      exit(0);
    }
  }

  fd->type = type ? FD_TYPE_SPE:FD_TYPE_FID;
  FD_alloc1ddata(fd);
  
  fp = fopen(fname, "r");
  /* Skip header of file */
  fread(head, sizeof(float), FDATASIZE, fp);
  switch (testHdr(head)) {
    case HDR_BAD:
      fprintf(stderr,"FD_read_nmrpipe: Unable to read header of file '%s'\n",
              fname);
      exit (0);
      break;
    case HDR_SWAPPED:
      swap = 1;
      break;
    case HDR_OK:
    default:
      swap = 0;
      break;
  }

  if (swap) {
    freadpipe = freadreversedbyteorder;
  } else {
    freadpipe = fread;
  }

  if (type == FD_TYPE_SPE) {
    if (ndim == 2) {
      for (j=fd->ni+1; j>1; j--) {
        if (j > fd->ni) jj = 1;
	else jj = j;
        no = (jj-1)*fd->np;
        for (i=no+fd->np+1; i>no+1; i--) {
          if (i > no+fd->np) ii = no+1;
	  else ii = i;
	  freadpipe(&f, sizeof(float), 1, fp);
          fd->data[ii].re = f;
        }
        for (i=no+fd->np; i>no; i--) {
          if (i > no+fd->np) ii = no+1;
	  else ii = i;
          freadpipe(&f, sizeof(float), 1, fp);
          fd->data[ii].im = f;
        }
      }
    } else {
      for (i=fd->np+1; i>1; i--) {
        if (i > fd->np) ii = 1;
        else ii = i;
        freadpipe(&f, sizeof(float), 1, fp);
        fd->data[ii].re = f;
      }
      for (i=fd->np; i>=1; i--) {
        if (i > fd->np) ii = 1;
        else ii = i;
        freadpipe(&f, sizeof(float), 1, fp);
        fd->data[ii].im = f;
      }
    }
  } else {
    if (ndim == 2) {
      for (j=1; j<=fd->ni; j++) {
        no = (j-1)*fd->np+1;
        for (i=no; i<no+fd->np; i++) {
          freadpipe(&f, sizeof(float), 1, fp);
          fd->data[i].re = f;
        }
        for (i=no; i<no+fd->np; i++) {
          freadpipe(&f, sizeof(float), 1, fp);
          fd->data[i].im = f;
        }
      }
    } else {
      for (i=1; i<=fd->np; i++) {
        freadpipe(&f, sizeof(float), 1, fp);
        fd->data[i].re = f;
      }
      for (i=1; i<=fd->np; i++) {
        freadpipe(&f, sizeof(float), 1, fp);
        fd->data[i].im = f;
      }
    }
  }
  return fd;
}
#undef FIRST
#undef LAST
