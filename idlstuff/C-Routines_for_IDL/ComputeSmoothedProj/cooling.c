#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#include "cooling.h"



/* returns new internal energy per unit mass. 
 * Arguments are passed in code units, density is proper density.
 */     
double DoCooling(double u_old, double rho, double dt, double *ne_guess) 
{

  /*  Here you can put your own implementation of cooling/heating....
   *
   *  Have fun!
   */

  return u_old;
}


/*  this function computes the self-consistent temperature
 *  and abundance ratios 
 */
double AbundanceRatios(double u, double rho, double *ne_guess, double *nH0_pointer)
{
  double temp=0;

  /* Again, you may fill in this yourself!
   */
  
  return temp;
}
 


void IonizeParams(void)
{

}

void InitCool(void)
{
  IonizeParams();
}






















