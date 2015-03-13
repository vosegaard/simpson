/***
 * namelist2: defines structure for various lists mapping keywords to values.
 ***/

#ifndef __namelist2_h

#define __namelist2_h

int intStrArgD(), fltStrArgD(), getValByName(), getNameByVal();

struct NameVal
   {
    char *name;
    float val;
   };

#endif
