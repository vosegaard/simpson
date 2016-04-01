#ifndef __RFSHAPES_H
#define __RFSHAPES_H


/* maximal number of rf shapes */
#define MAXRFSHAPES 100

 /* Define rf shape element structure */
 typedef struct _RFelem {
     double ampl;
     double phase;
 } RFelem;


extern RFelem* RFshapes[];


void RFshapes_init(void);
void RFshapes_reset(void);
int RFshapes_len(int slot);
RFelem* RFshapes_alloc(int len);











#endif
