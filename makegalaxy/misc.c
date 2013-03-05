
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"


#include "prototypes.h"
#include "globvars.h"




double splint_xl_yl_D2yl(double t)
{
  double value;

  splint(xl, yl, D2yl, ZSIZE + 1, t, &value);

  return value;
}







/* assuming f(x) is a smooth function of x */
int find_idx(double x, double *fx, int size)
{
  int hi, low, idx, incr, iter= 0;

  if((x<fx[0]) && (x<fx[size])) {
    //assert(0);
    return -1;
  }

  if((x>fx[0]) && (x>fx[size])) {
    //assert(0);
    return -1;
  }

  hi= size;
  low= 0;
  idx= size/2;
  incr=fx[size]>fx[0];
  do {
        if(fx[idx] > x)
          {
	    /* is it positively increasing or not */
	    if(incr)
              hi= idx;
	    else
	      low= idx;
          }
        else
          {
	    if(incr)
              low= idx;
	    else
	      hi= idx;
          }

        idx= (hi-low)/2 + low;
        iter++;

  } while ((hi-low) > 1);

  /* always return low index */
  return low;
}

