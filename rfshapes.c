 #include <stdlib.h>
 #include <stdio.h>
 #include <tcl.h>
 #include <math.h>
 #include <ctype.h>
 #include <unistd.h>
 #include <string.h>
 #include <time.h>
 #include "defs.h"
 #include "tclutil.h"
 #include "rfshapes.h"
 
 /* define global array of pointers to rf shapes */
 RFelem *RFshapes[MAXRFSHAPES];
 
/****
 * initialization of RFshapes
 ****/
  void RFshapes_init() {
     int a;
     
     for (a=0; a<MAXRFSHAPES; a++) {
        RFshapes[a] = NULL;
     }
  }
 
 
 /****
  * provide free slot in RFshapes
  ****/
 int RFshapes_slot() {
   int a;
   
   for (a=0; a<MAXRFSHAPES; a++) {
      if (!RFshapes[a]) {
         break;
      }
   }
   if (a >= MAXRFSHAPES) {
      fprintf(stderr,"RFshapes error: no more free slots available\n");
      a=-1;
   }
   return a;
 }
 
 
 
/****
 * allocation function for rf shapes 
 ****/
 RFelem* RFshapes_alloc(int len) {
   RFelem* v;
   
   
   v = (RFelem*)malloc((len+1)*sizeof(RFelem));
   if (!v) {
      fprintf(stderr,"error: unable to alocate RFshape");
      exit(-1);
   }
   
   /* store its length to the first element */
   *(int*)v=len;
   return v;
 }
 
 
 /****
  * freeing the shape from memory 
  ****/
 void free_RFshapes(int a) {
 
    free((char *)RFshapes[a]);
    RFshapes[a] = NULL;
 }

/****
 * free all remaining slots in RFshapes
 ****/
 void RFshapes_reset() {
   int a;
   
   for (a=0; a<MAXRFSHAPES; a++) {
     if (RFshapes[a] != NULL) {
       free_RFshapes(a);
       RFshapes[a] = NULL;
     }
   }
 }
 
 
 
/****
 * length of shape in a given slot
 ****/
 int RFshapes_len(int slot) {
 
    if (!RFshapes[slot]) {
      fprintf(stderr,"RFshapes error: testing length of empty slot\n");
      exit(1);
    }
    return *(int*)RFshapes[slot];
 }
 
 /****
  * create RFshape from file
  ****/
 int load_RFshape(const char* name){
   FILE* fp;
   char fname[256], dum[256];
   int Nelem, RFidx, i;
   double am, ph;
    
   /* decide which slot to use */
   RFidx = RFshapes_slot();
    
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
   fp=fopen(fname,"r");
   if (!fp) {
      fprintf(stderr,"load_shape error: unable to open file %s\n\n",fname);
      exit(1);
   }
   
   /* scan file for number of lines, i.e. number of elements in rf shape */
   Nelem = 0;
   while ( fgets(dum, 256, fp) ) {
      Nelem++;
   }
   fseek(fp, 0, SEEK_SET);
   /* printf("Number of lines = %d\n",Nelem); */
 
   RFshapes[RFidx]=RFshapes_alloc(Nelem);
   
   for (i=1; i<=Nelem; i++) {
      fgets(dum, 256, fp);
      if ( sscanf(dum,"%lg%lg",&am,&ph) != 2 ) {
         fprintf(stderr,"load_shape error: unable to read line %d in %s\n",i,fname);
	 exit(1);
      }
      RFshapes[RFidx][i].ampl = am;
      RFshapes[RFidx][i].phase = ph;   
   }
   fclose(fp);
   
   return RFidx;
 }
 
 
 /****
  * save slot to file in ASCII
  ****/
  void save_RFshape(int slot, const char* name) {
    FILE* fp;
   char fname[256];
   int  i, Nelem;
   
   /* does the slot hold a shape? */
   if (!RFshapes[slot]) {
      fprintf(stderr,"save_shape error: the shape does not exist\n");
      exit(1);
   }
   Nelem = RFshapes_len(slot); 
    
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
   fp=fopen(fname,"w");
   if (!fp) {
      fprintf(stderr,"save_shape error: unable to create file %s\n\n",fname);
      exit(1);
   }
    
   for (i=1;i<=Nelem; i++) {
     fprintf(fp,"%.15g %.15g\n",RFshapes[slot][i].ampl,RFshapes[slot][i].phase);
   }
   
   fclose(fp); 
 }
 
 
 /*======================================================*/

/****
 * implementation of Tcl free_shape routine 
 ****/
int tclFreeShape(ClientData data,Tcl_Interp* interp,int objc, Tcl_Obj *objv[])
{
   int slot;
   
  if (objc != 2)
    return TclError(interp,"usage: free_shape <RFshape>");
  if (Tcl_GetIntFromObj(interp,objv[1],&slot) == TCL_ERROR)
    return TclError(interp,"free_shape: argument 1 must be integer <RFshape>");

  /* test if slot exists */
  if (slot>MAXRFSHAPES)
     return TclError(interp,"free_shape: provided slot exceeds slots in memory");
  if (RFshapes[slot]==NULL) {
/*     return TclError(interp,"free_shape: shape seems not to exist");  */
     fprintf(stderr,"free_shape warning: shape seems not to exist\n");
  } else {
     free_RFshapes(slot);
  }

  return TCL_OK;
}

/****
 * implementation of Tcl free_all_shapes routine 
 ****/
int tclFreeAllShapes(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
    
  RFshapes_reset();
  
  return TCL_OK;
}

/****
 * implementation of Tcl shape_len routine 
 ****/
int tclShapeLen(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
   int slot;
   
  if (argc != 2) 
    return TclError(interp,"usage: <int> shape_len <RFshape>");
  if (Tcl_GetIntFromObj(interp,argv[1],&slot) == TCL_ERROR)
    return TclError(interp,"shape_len: argument 1 must be integer <RFshape>");
  
  /* test existence of slot */
  if (!RFshapes[slot])
     return TclError(interp,"shape_len: shape seems not to exist");

  return TclSetResult(interp,"%d",RFshapes_len(slot));
}



/****
 * implementation of Tcl load_shape routine 
 ****/
int tclLoadShape(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
   int slot;
   char * fname;
   
  if (argc != 2) 
    return TclError(interp,"usage: <RFshape> load_shape <name of file>");

   fname = Tcl_GetString(argv[1]);
   slot = load_RFshape(fname);

  return TclSetResult(interp,"%d",slot);
}



/****
 * implementation of Tcl save_shape routine 
 ****/
int tclSaveShape(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
   int slot;
   char *fname;
   
  if (argc != 3) 
    return TclError(interp,"usage: save_shape <RFshape> <name of file>");
  if (Tcl_GetIntFromObj(interp,argv[1],&slot) == TCL_ERROR)
    return TclError(interp,"save_shape: argument 1 must be integer <RFshape>");
    
  fname = Tcl_GetString(argv[2]);
  save_RFshape(slot,fname);

  return TCL_OK;
}

 
 /****
 *  routine for generating random rf shape
 *     create, fill the slot, return index of that slot
 ****/
int tclRandShape(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  complx pol3int(double x,double x1,double dx, complx y1, complx y2, complx y3, complx y4);

  int newnp, np, i, j, slot;
  long ltime, ltimebuf1;
  static long ltimebuf2=1;
  static unsigned long cntr=0;
  double ampl;
  complx *dat;

  if (argc != 4) 
    return TclError(interp,"usage: <rfshape> rand_shape <max. amplitude> <# of elements> <# of random elements>");

  if (Tcl_GetDoubleFromObj(interp,argv[1],&ampl) == TCL_ERROR)
    return TclError(interp,"rand_shape: argument 1 must be double <max. amplitude>");
  if (Tcl_GetIntFromObj(interp,argv[2],&newnp) == TCL_ERROR)
    return TclError(interp,"rand_shape: argument 2 must be integer <# of elements>");
  if (Tcl_GetIntFromObj(interp,argv[3],&np) == TCL_ERROR)
    return TclError(interp,"rand_shape: argument 3 must be integer <# of random elements>");
  if (np>newnp)
    return TclError(interp,"rand_shape: number of random elements must not exceed total number of elements");


  // get a new slot and allocate
  slot = RFshapes_slot();
  if (slot == -1) 
     return TclError(interp,"rand_shape error: no more free slots available, free some shape first!");
  RFshapes[slot] = RFshapes_alloc(newnp);

  // generate random points
  dat=(complx*)malloc(sizeof(complx)*(np+1));
  // initialize random sequence
  cntr++;
  ltime = time(NULL);
  ltimebuf1 = ltime;
   if (ltimebuf2 == ltime) {
      ltime /= cntr;
   }
   if ( (ltimebuf2 == 1) || (ltimebuf2 != ltimebuf1) ) {
    ltimebuf2 = ltimebuf1;
    cntr = 1;
   } 
  srand( (unsigned int) ltime);
  for (i=1;i<=np;i++) {
    dat[i].re = ((double)rand()/((double)RAND_MAX+1.0))*ampl;
    dat[i].im = (((double)rand()/((double)RAND_MAX+1.0))-0.5)*360.0;
    // DEBUG output
    //printf("%g %f\n",dat[i].re, dat[i].im);
    //
  }


  if (newnp == np) {
    // don't do spline, just copy
    for (i=1; i<=newnp; i++) {
      if (dat[i].re < 0) {
         dat[i].re = 0.0;
         dat[i].im = 0.0;
      }
      if (dat[i].re > ampl) dat[i].re=ampl;

      RFshapes[slot][i].ampl = dat[i].re;
      RFshapes[slot][i].phase = dat[i].im;
    }
  } else {
     // do spline interpolation
     double dx1, dx2, r, x1, x2, x;
     int np2;
     complx *newdata;
       
     x1 = 1;
     dx1 = 1.0;
     x2 = 1;
     dx2 = (double)(np-1)/(double)(newnp-1)*dx1;
     newdata=(complx*)malloc(sizeof(complx)*(newnp+1));
     np2 = np-2;
     r = dx2/dx1;
  
     for (j=1;j<=newnp; j++) {
       i= (int)((double)(j-1)*r)+1;
       x= (j-1)*dx2+x2;
       if (i<2) {
         newdata[j] = pol3int(x, x1, dx1, dat[1],dat[2], dat[3], dat[4]);
       } else if (i>np2) {
         newdata[j] = pol3int(x, x1+(np-4)*dx1, dx1, dat[np-3],dat[np-2], dat[np-1], dat[np]);
       } else {
         newdata[j] = pol3int(x, (i-2)*dx1+x1, dx1, dat[i-1],dat[i], dat[i+1], dat[i+2]);
       }
       if (newdata[j].re < 0) {
         newdata[j].re = 0.0;
         newdata[j].im = 0.0;
       }
       if (newdata[j].re > ampl) newdata[j].re=ampl;
       if (newdata[j].im > 180.0) newdata[j].im=180.0;
       if (newdata[j].im < -180.0) newdata[j].im=-180.0;
       

       RFshapes[slot][j].ampl = newdata[j].re;
       RFshapes[slot][j].phase = newdata[j].im;
       
     }  
     free(newdata);
   }
   free(dat);
 
  // DEBUG output
  //printf("=========\n");
  //for (i=1;i<=newnp;i++) {
  //  printf("%g %f\n",RFshapes[slot][i].ampl, RFshapes[slot][i].phase);
  //}

  
  // set Tcl result to slot index
  //sprintf(interp->result,"%d",slot);
  //return TCL_OK;
  return TclSetResult(interp,"%d",slot);

}


/****
 * implementation of Tcl shape_index routine 
 ****/
int tclShapeIndex(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
   int slot, idx;
   char *buf;
   
  if (argc != 3 && argc != 4)
    return TclError(interp,"Usage: <value> shape_index <RFshape> <index> ?-ampl|-phase?");

  if (Tcl_GetIntFromObj(interp,argv[1],&slot) == TCL_ERROR)
    return TclError(interp,"shape_index: argument 1 must be integer <RFshape>");
 
  /* check for RFshape existence */
  if (!RFshapes[slot])
     return TclError(interp,"shape_index: trying to acces non-existing RFshape");
     
  /* read second argument */
  if (Tcl_GetIntFromObj(interp,argv[2],&idx) == TCL_ERROR)
    return TclError(interp,"shape_index: argument 2 must be integer <index>");
  if (idx<1 || idx>RFshapes_len(slot))
    return TclError(interp,"shape_index_set: index out of shape size");

  if (argc == 4) {
    buf = Tcl_GetString(argv[3]);
    if (!strcmp(buf,"-ampl") )
       return TclSetResult(interp,"%.20f",RFshapes[slot][idx].ampl);
    else if (!strcmp(buf,"-phase") )
       return TclSetResult(interp,"%.20f",RFshapes[slot][idx].phase);
    else
       return TclError(interp,"shape_index: argument 3 must be either '-ampl' or '-phase'");
   } else
      return TclSetResult(interp,"%.20f %.20f",RFshapes[slot][idx].ampl, RFshapes[slot][idx].phase);
}

/****
 * implementation of Tcl shape_index_set routine 
 ****/
int tclShapeIndexSet(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
   int slot, idx;
   double am, ph;
   char *buf;
   
  if (argc != 5 && argc != 7)
    return TclError(interp,"Usage: shape_index_set <RFshape> <index> ?-ampl <a>? ?-phase <p>?");

  if (Tcl_GetIntFromObj(interp,argv[1],&slot) == TCL_ERROR)
    return TclError(interp,"shape_index_set: argument 1 must be integer <RFshape>");
 
  /* check for RFshape existence */
  if (!RFshapes[slot])
     return TclError(interp,"shape_index_set: trying to acces non-existing RFshape");
     
  /* read second argument */
  if (Tcl_GetIntFromObj(interp,argv[2],&idx) == TCL_ERROR)
    return TclError(interp,"shape_index_set: argument 2 must be integer <index>");
  if (idx<1 || idx>RFshapes_len(slot))
    return TclError(interp,"shape_index_set: index out of shape size");

  buf = Tcl_GetString(argv[3]);
  if (!strcmp(buf,"-ampl") ) {
    if (Tcl_GetDoubleFromObj(interp,argv[4],&am) == TCL_ERROR)
      return TclError(interp,"shape_index_set: -ampl must be double");
    RFshapes[slot][idx].ampl = am;
  } else if (!strcmp(buf,"-phase") ) {
    if (Tcl_GetDoubleFromObj(interp,argv[4],&ph) == TCL_ERROR)
      return TclError(interp,"shape_index_set: -phase must be double");
    RFshapes[slot][idx].phase = ph;
  } else
       return TclError(interp,"shape_index_set: argument 3 must be either '-ampl' or '-phase'");

  if ( argc == 5) return TCL_OK;
  buf = Tcl_GetString(argv[5]);
  if (!strcmp(buf,"-ampl") ) {
    if (Tcl_GetDoubleFromObj(interp,argv[6],&am) == TCL_ERROR)
      return TclError(interp,"shape_index_set: -ampl must be double");
    RFshapes[slot][idx].ampl = am;
  } else if (!strcmp(buf,"-phase") ) {
    if (Tcl_GetDoubleFromObj(interp,argv[6],&ph) == TCL_ERROR)
      return TclError(interp,"shape_index_set: -phase must be double");
    RFshapes[slot][idx].phase = ph;
  } else
       return TclError(interp,"shape_index_set: argument 5 must be either '-ampl' or '-phase'");

  return TCL_OK;
}

/****
 * implementation of list2shape (creates RFshape from a list { {a p} {a p} ... }
 ****/
int tclList2Shape(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  Tcl_Obj **list1, **list2;
  int nlist1, nlist2, i, slot;
  
  if (argc != 2)
    return TclError(interp,"Usage: <RFshape> list2shape { {a1 p1} {a2 p2} ... }");

  if (Tcl_ListObjGetElements(interp,argv[1],&nlist1,&list1) != TCL_OK)
     return TclError(interp,"list2shape: unable to decompose list argument");

  /* get a new slot and allocate */
  slot = RFshapes_slot();
  if (slot == -1) {
     return TclError(interp,"list2shape error: no more free slots available, free some shape first!");
  }
  RFshapes[slot] = RFshapes_alloc(nlist1);


  for (i=0; i<nlist1; i++) {
     if (Tcl_ListObjGetElements(interp,list1[i],&nlist2,&list2) != TCL_OK) {
          return TclError(interp,"list2shape can not read list element %d",i+1);
     }
     if (nlist2 != 2) {
        return TclError(interp,"list2shape: expecting two elements like {amplitude phase} in list");
     }
     if (Tcl_GetDoubleFromObj(interp,list2[0],&RFshapes[slot][i+1].ampl) != TCL_OK) {
        return TclError(interp,"lis2shape cannot interpret amplitude in element %d",i+1);
     }
     if (Tcl_GetDoubleFromObj(interp,list2[1],&RFshapes[slot][i+1].phase) != TCL_OK) {
        return TclError(interp,"lis2shape cannot interpret phase in element %d",i+1);
     }
  }

  return TclSetResult(interp,"%d",slot);
}
 
/****
 * implementation of shape2list (from RFshape creates a list { {a p} {a p} ... }
 ****/
int tclShape2List(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  Tcl_Obj *lptr1, *lptr2;
  Tcl_Obj *elemptr[2];
  int i, slot;

  if (argc != 2)
    return TclError(interp,"Usage: <list> shape2list <RFshape>");

  if (Tcl_GetIntFromObj(interp,argv[1],&slot) == TCL_ERROR)
    return TclError(interp,"shape2list: argument must be integer <RFshape>");
 
  /* check for RFshape existence */
  if (!RFshapes[slot])
     return TclError(interp,"shape2list: trying to acces non-existing RFshape");

  /* create list objects */
  lptr1 = Tcl_NewListObj(0,NULL);
  if (!lptr1) return TclError(interp,"shape2list unable to create outer list");
  
  for (i=1; i<=RFshapes_len(slot); i++) {
     elemptr[0] = Tcl_NewDoubleObj(RFshapes[slot][i].ampl);
     if (!elemptr[0]) {
        return TclError(interp,"shape2list unable to create double from RFshape amplitude element %d",i);
     }
     elemptr[1] = Tcl_NewDoubleObj(RFshapes[slot][i].phase);
     if (!elemptr[1]) {
        return TclError(interp,"shape2list unable to create double from RFshape amplitude element %d",i);
     }
     lptr2 = Tcl_NewListObj(2,elemptr);
     if (!lptr2) return TclError(interp,"shape2list unable to create inner list");

     if ( Tcl_ListObjAppendElement(interp,lptr1,lptr2) != TCL_OK ) {
        return TclError(interp,"shape2list unable to append element %d to oute list",i);
     }
  }
  
  Tcl_SetObjResult(interp,lptr1);
  
  return TCL_OK;
} 

/****
 * implementation of shape_join (creates RFshape by joining all argument shapes}
 ****/
int tclShapeJoin(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  int i, j, ll, slot;
  int *jsh = malloc(argc*sizeof(int));
  
  ll = 0;
  for (i=1; i<argc; i++) {
     if (Tcl_GetIntFromObj(interp,argv[i],&slot) == TCL_ERROR) {
        return TclError(interp,"shape_join: all arguments must be integers <RFshape>");
     }
     /* check for RFshape existence */
     if (!RFshapes[slot]) {
        return TclError(interp,"shape_join: trying to access non-existing RFshape in argument %d",i);
     }
     ll += RFshapes_len(slot);
     jsh[i]=slot;
  }
  
  /* get a new slot and allocate */
  slot = RFshapes_slot();
  if (slot == -1) {
     return TclError(interp,"shape_join error: no more free slots available, free some shape first!");
  }
  RFshapes[slot] = RFshapes_alloc(ll);
  
  /* fill this slot */
  ll=0;
  for (i=1; i<argc; i++) {
     for (j=1; j<=RFshapes_len(jsh[i]); j++) {
        ll++;
	//RFshapes[slot][ll].ampl = RFshapes[jsh[i]][j].ampl;
	//RFshapes[slot][ll].phase = RFshapes[jsh[i]][j].phase;
        RFshapes[slot][ll] = RFshapes[jsh[i]][j];
     }
  }
  
  free(jsh);
  return TclSetResult(interp,"%d",slot);
}

/****
 * implementation of shape_dup (creates new RFshape by adding phase to existent RFshape}
 ****/
int tclShapeDup(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  int slot, newslot, i, N;
  double ph;
  double k=1.0;
  
  if ( (argc < 3) || (argc > 4) )
    return TclError(interp,"Usage: <RFshape> shape_dup <RFshape> <phase> ?<ampl. scale factor>?");

  if (Tcl_GetIntFromObj(interp,argv[1],&slot) == TCL_ERROR)
    return TclError(interp,"shape_dup: first argument must be integer <RFshape>");
 
  /* check for RFshape existence */
  if (!RFshapes[slot])
     return TclError(interp,"shape_dup: trying to acces non-existing RFshape");

  if (Tcl_GetDoubleFromObj(interp,argv[2],&ph) == TCL_ERROR)
    return TclError(interp,"shape_dup: second argument must be double <phase in deg.>");

  if (argc == 4) {
    if (Tcl_GetDoubleFromObj(interp,argv[3],&k) == TCL_ERROR)
      return TclError(interp,"shape_dup: third argument must be double <ampl. scale factor>");
  }

  /* get a new slot and allocate */
  newslot = RFshapes_slot();
  if (newslot == -1) {
     return TclError(interp,"shape_dup error: no more free slots available, free some shape first!");
  }
  N = RFshapes_len(slot);
  RFshapes[newslot] = RFshapes_alloc(N);

  for (i=1; i<=N; i++) {
     RFshapes[newslot][i].ampl = (RFshapes[slot][i].ampl)*k;
     RFshapes[newslot][i].phase = RFshapes[slot][i].phase+ph;
     if ( RFshapes[newslot][i].phase < 0.0 )
        RFshapes[newslot][i].phase += 360.0;
     if ( RFshapes[newslot][i].phase > 360.0 )
        RFshapes[newslot][i].phase -= 360.0;
  }
  
  return TclSetResult(interp,"%d",newslot);
}

/****
 * implementation of shape_ampl (RFshape amplitude statistics)
 ****/
int tclShapeAmpl(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  int slot, i, N;
  double dum, avg, min, max, rms;
   
  if (argc<3)
     return TclError(interp,"Usage: <result> shape_ampl <RFshape> ?[-avg|-min|-max|-rms]?");
      
  if (Tcl_GetIntFromObj(interp,argv[1],&slot) == TCL_ERROR)
    return TclError(interp,"shape_ampl: first argument must be integer <RFshape>");
 
  /* check for RFshape existence */
  if (!RFshapes[slot])
     return TclError(interp,"shape_ampl: trying to access non-existing RFshape");
     
  min = 1.0e6;
  max = -1.0;
  avg = 0.0;
  rms = 0.0;
  N = RFshapes_len(slot);
  for (i=1; i<=N; i++) {
     dum = RFshapes[slot][i].ampl;
     avg += dum;
     rms += dum*dum;
     if (dum < min) min = dum;
     if (dum > max) max = dum;
  }
  avg /= (double)N;
  rms = sqrt(rms/((double)N));
  
  /* create output */
  Tcl_ResetResult(interp);
  for (i=2; i<argc; i++) {
     if (!strcmp(Tcl_GetString(argv[i]),"-avg")) {
        TclAppendResult(interp,"%.15g",avg);
     } else if (!strcmp(Tcl_GetString(argv[i]),"-min")) {
        TclAppendResult(interp,"%.15g",min);
     } else if (!strcmp(Tcl_GetString(argv[i]),"-max")) {
        TclAppendResult(interp,"%.15g",max);
     } else if (!strcmp(Tcl_GetString(argv[i]),"-rms")) {
        TclAppendResult(interp,"%.15g",rms);
     } else {
        return TclError(interp,"shape_ampl: unknown argument '%s'",Tcl_GetString(argv[i]));
     }
  }
  
  return TCL_OK;
}

/****
 * implementation of shape_energy [rf in Hz, energy in (rad.s-1)^2 ]
 ****/
int tclShapeEnergy(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  int slot, i, N;
  double nrg, dur, dum;
   
  if (argc != 3)
     return TclError(interp,"Usage: <result> shape_energy <RFshape> <duration>");
      
  if (Tcl_GetIntFromObj(interp,argv[1],&slot) == TCL_ERROR)
    return TclError(interp,"shape_energy: first argument must be integer <RFshape>");
  /* check for RFshape existence */
  if (!RFshapes[slot])
     return TclError(interp,"shape_energy: trying to acces non-existing RFshape");

  if (Tcl_GetDoubleFromObj(interp,argv[2],&dur) == TCL_ERROR)
    return TclError(interp,"shape_energy: second argument must be double <duration>");
  if (dur <= 0.0)
    return TclError(interp,"shape_energy: duration should be greater than zero");

  nrg = 0.0;
  N = RFshapes_len(slot);
  for (i=1; i<=N; i++) {
     dum = RFshapes[slot][i].ampl;
     nrg += dum*dum;
  }
  nrg = 4.0*M_PI*M_PI*nrg*dur*1e-6/(double)N;

  return TclSetResult(interp,"%.15g",nrg);

}


/****
 * implementation of Tcl shape_manipulate routine 
 ****/
int tclShapeManipulate(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
   int slot, i, j, k, len, shft;
   char expr[2048], *buf;
   double a, p;
   Tcl_Obj *tclres;

  if ( (argc < 3) || (argc > 10) )
    return TclError(interp,"usage: shape_manipulate <RFshape> -ampl <ampl expr> | -phase <phase expr> | -time_reversal | -phase_invert | -cyclic_shift N");
  if (Tcl_GetIntFromObj(interp,argv[1],&slot) == TCL_ERROR)
    return TclError(interp,"save_manipulate: argument 1 must be integer <RFshape>");
  /* check for RFshape existence */
  if (!RFshapes[slot])
     return TclError(interp,"shape_manipulate: trying to access non-existing RFshape");

  len = RFshapes_len(slot);
  for (i=2; i<argc; i++) {
      buf = Tcl_GetString(argv[i]);
	  if (!strcmp(buf,"-time_reversal")) {
        for (j=1; j<= len/2; j++) {
           k = len+1-j;
           a = RFshapes[slot][k].ampl;
           p = RFshapes[slot][k].phase;
           RFshapes[slot][k].ampl = RFshapes[slot][j].ampl;
           RFshapes[slot][k].phase = RFshapes[slot][j].phase;
           RFshapes[slot][j].ampl = a;
           RFshapes[slot][j].phase = p;
        }
     } else if (!strcmp(buf,"-phase_invert")) {
        for (j=1; j<=len; j++) {
           RFshapes[slot][j].phase = -RFshapes[slot][j].phase;
        }
     } else if (!strcmp(buf,"-cyclic_shift")) {
        RFelem* v;

        i++;
        if (Tcl_GetIntFromObj(interp,argv[i],&shft) == TCL_ERROR)
           return TclError(interp,"shape_manipulate: argument following -cyclic_shift must be integer");
        v = (RFelem*)malloc((len+1)*sizeof(RFelem));
        if (!v) {
           fprintf(stderr,"error in shape_manipulate: unable to allocate temporary RFshape");
           exit(-1);
        }
        /* store its length to the first element */
        *(int*)v=len;
        for (j=1; j<=len; j++) {
           v[j].ampl = RFshapes[slot][j].ampl;
           v[j].phase = RFshapes[slot][j].phase;
        }
        shft = shft % len;
        for (j=1; j<= len; j++) {
           k = (len + shft + j) % len;
           if (k <= 0) k += len;
           RFshapes[slot][k].ampl = v[j].ampl;
           RFshapes[slot][k].phase = v[j].phase;
        }
        free((char *)v);
     } else if (!strcmp(buf,"-ampl")) {
        i++;
     /* evaluate expression in Tcl */
        for (j=1; j<=len; j++) {
           sprintf(expr,"\n set i %d\n set ampl %f\n set phase %f\n expr %s\n", j, RFshapes[slot][j].ampl, RFshapes[slot][j].phase, Tcl_GetString(argv[i]));
           if ( Tcl_EvalEx(interp, expr, -1,TCL_EVAL_DIRECT) != TCL_OK ) 
              return TclError(interp,"error in shape_manipulate: can not evaluate %s for index %d",Tcl_GetString(argv[i]), j);
           tclres = Tcl_GetObjResult(interp);
           if ( Tcl_GetDoubleFromObj(interp,tclres,&a) != TCL_OK )
              return TclError(interp,"error in shape_manipulate: can not get amplitude result for index %d",j);
           RFshapes[slot][j].ampl = a;
        }
     } else if (!strcmp(buf,"-phase")) {
        i++;
     /* evaluate expression in Tcl */
        for (j=1; j<=len; j++) {
           sprintf(expr,"\n set i %d\n set ampl %f\n set phase %f\n expr %s\n", j, RFshapes[slot][j].ampl, RFshapes[slot][j].phase, Tcl_GetString(argv[i]));
           if ( Tcl_EvalEx(interp, expr, -1,TCL_EVAL_DIRECT) != TCL_OK ) 
              return TclError(interp,"error in shape_manipulate: can not evaluate %s for index %d",Tcl_GetString(argv[i]), j);
           tclres = Tcl_GetObjResult(interp);
           if ( Tcl_GetDoubleFromObj(interp,tclres,&p) != TCL_OK )
              return TclError(interp,"error in shape_manipulate: can not get amplitude result for index %d",j);
           RFshapes[slot][j].phase = p;
        }
     }
  }   

  return TCL_OK;
}

/****
 * implementation of Tcl shape_create routine 
 ****/
int tclShapeCreate(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
   int slot, i, j, len;
   char expr[2048], *buf;
   double a, p;
   Tcl_Obj *tclres;

  if ( (argc < 2) || (argc > 6) )
    return TclError(interp,"usage: shape_create <num of lems> | -ampl <ampl expr> | -phase <phase expr>");
  if (Tcl_GetIntFromObj(interp,argv[1],&len) == TCL_ERROR)
    return TclError(interp,"shape_create: argument 1 must be integer <num of elems>");

  /* get a new slot and allocate */
  slot = RFshapes_slot();
  if (slot == -1) {
     return TclError(interp,"shape_create error: no more free slots available, free some shape first!");
  }
  RFshapes[slot] = RFshapes_alloc(len);

  for (i=1; i<=len; i++) {
     RFshapes[slot][i].ampl = 0.0;
     RFshapes[slot][i].phase = 0.0;
  }

  for (i=2; i<argc; i++) {
     buf = Tcl_GetString(argv[i]);
	 if (!strcmp(buf,"-ampl")) {
        i++;
        /* evaluate expression in Tcl */
        for (j=1; j<=len; j++) {
           sprintf(expr,"\n set i %d\n expr %s\n", j, Tcl_GetString(argv[i]));
           if ( Tcl_EvalEx(interp, expr, -1,TCL_EVAL_DIRECT) != TCL_OK ) 
              return TclError(interp,"error in shape_create: can not evaluate %s for index %d",Tcl_GetString(argv[i]), j);
           tclres = Tcl_GetObjResult(interp);
           if ( Tcl_GetDoubleFromObj(interp,tclres,&a) != TCL_OK )
              return TclError(interp,"error in shape_create: can not get amplitude result for index %d",j);
           RFshapes[slot][j].ampl = a;
        }
     } else if (!strcmp(buf,"-phase")) {
        i++;
        /* evaluate expression in Tcl */
        for (j=1; j<=len; j++) {
           sprintf(expr,"\n set i %d\n expr %s\n", j, Tcl_GetString(argv[i]));
           if ( Tcl_EvalEx(interp, expr, -1,TCL_EVAL_DIRECT) != TCL_OK ) 
              return TclError(interp,"error in shape_create: can not evaluate %s for index %d",Tcl_GetString(argv[i]), j);
           tclres = Tcl_GetObjResult(interp);
           if ( Tcl_GetDoubleFromObj(interp,tclres,&p) != TCL_OK )
              return TclError(interp,"error in shape_create: can not get amplitude result for index %d",j);
           RFshapes[slot][j].phase = p;
        }
     }
  }   

  return TclSetResult(interp,"%d",slot);
}



/* implement new commands */
void tclcmd_rfshape(Tcl_Interp* interp) {
 
   Tcl_CreateObjCommand(interp,"load_shape",(Tcl_ObjCmdProc *)tclLoadShape,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateObjCommand(interp,"save_shape",(Tcl_ObjCmdProc *)tclSaveShape,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateObjCommand(interp,"rand_shape",(Tcl_ObjCmdProc *)tclRandShape,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateObjCommand(interp,"free_shape",(Tcl_ObjCmdProc *)tclFreeShape,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateObjCommand(interp,"free_all_shapes",(Tcl_ObjCmdProc *)tclFreeAllShapes,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateObjCommand(interp,"shape_len",(Tcl_ObjCmdProc *)tclShapeLen,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateObjCommand(interp,"shape_index",(Tcl_ObjCmdProc *)tclShapeIndex,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateObjCommand(interp,"list2shape",(Tcl_ObjCmdProc *)tclList2Shape,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateObjCommand(interp,"shape2list",(Tcl_ObjCmdProc *)tclShape2List,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateObjCommand(interp,"shape_join",(Tcl_ObjCmdProc *)tclShapeJoin,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateObjCommand(interp,"shape_dup",(Tcl_ObjCmdProc *)tclShapeDup,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateObjCommand(interp,"shape_ampl",(Tcl_ObjCmdProc *)tclShapeAmpl,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateObjCommand(interp,"shape_energy",(Tcl_ObjCmdProc *)tclShapeEnergy,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateObjCommand(interp,"shape_manipulate",(Tcl_ObjCmdProc *)tclShapeManipulate,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateObjCommand(interp,"shape_create",(Tcl_ObjCmdProc *)tclShapeCreate,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
   Tcl_CreateObjCommand(interp,"shape_index_set",(Tcl_ObjCmdProc *)tclShapeIndexSet,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
 }
