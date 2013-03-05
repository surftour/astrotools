#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "allvars.h"
#include "proto.h"



/* This routine computes various global properties of the particle
 * distribution and stores the result in the struct `SysState'.
 * Currently, not all the information that's computed here is 
 * actually used (e.g. momentum is not really used anywhere),
 * just the energies are written to a log-file every once in a while.
 */
void compute_global_quantities_of_system(void)
{
  int i, j, n;
  
  for(n=0; n<5; n++)   
    { 
      SysState.MassComp[n]= SysState.EnergyKinComp[n]= SysState.EnergyPotComp[n]= SysState.EnergyIntComp[n]=0; 
      
      for(j=0;j<4;j++) 
	SysState.CenterOfMassComp[n][j]= SysState.MomentumComp[n][j] = SysState.AngMomentumComp[n][j] = 0; 
    } 

  
  for(i=1;i<=NumPart;i++) 
    { 
      SysState.MassComp[P[i].Type] += P[i].Mass; 
      
      SysState.EnergyPotComp[P[i].Type] += 0.5*P[i].Mass*P[i].Potential; 
      

      if(P[i].Type==0)
	SysState.EnergyIntComp[0] += P[i].Mass * SphP[i].EgySpecPred; 

      SysState.EnergyKinComp[P[i].Type] += 0.5*P[i].Mass*(P[i].VelPred[0]*P[i].VelPred[0] +
							  P[i].VelPred[1]*P[i].VelPred[1] + 
							  P[i].VelPred[2]*P[i].VelPred[2]); 

      for(j=0;j<3;j++) 
	{ 
	  SysState.MomentumComp[P[i].Type][j] += P[i].Mass* P[i].VelPred[j]; 
	  SysState.CenterOfMassComp[P[i].Type][j] += P[i].Mass*P[i].PosPred[j]; 
	} 
      
      SysState.AngMomentumComp[P[i].Type][0] += P[i].Mass* ( P[i].PosPred[1]*P[i].VelPred[2] -
							  P[i].PosPred[2]*P[i].VelPred[1] ); 
      SysState.AngMomentumComp[P[i].Type][1] += P[i].Mass* ( P[i].PosPred[2]*P[i].VelPred[0] -
							  P[i].PosPred[0]*P[i].VelPred[2] ); 
      SysState.AngMomentumComp[P[i].Type][2] += P[i].Mass* ( P[i].PosPred[0]*P[i].VelPred[1] -
							  P[i].PosPred[1]*P[i].VelPred[0] ); 
    } 
  


  for(i=0;i<5;i++) 
    { 
      SysState.EnergyTotComp[i] = SysState.EnergyKinComp[i] + 
	                          SysState.EnergyPotComp[i] + 
	                          SysState.EnergyIntComp[i];
    } 
  
  
  SysState.Mass=SysState.EnergyKin=SysState.EnergyPot=SysState.EnergyInt=SysState.EnergyTot=0; 
  for(j=0;j<3;j++) 
    SysState.Momentum[j]=SysState.AngMomentum[j]=SysState.CenterOfMass[j]=0; 
  
  for(i=0;i<5;i++) 
    { 
      SysState.Mass+=SysState.MassComp[i]; 
      SysState.EnergyKin+=SysState.EnergyKinComp[i]; 
      SysState.EnergyPot+=SysState.EnergyPotComp[i]; 
      SysState.EnergyInt+=SysState.EnergyIntComp[i]; 
      SysState.EnergyTot+=SysState.EnergyTotComp[i]; 
      
      for(j=0;j<3;j++) 
	{ 
	  SysState.Momentum[j]+=SysState.MomentumComp[i][j]; 
	  SysState.AngMomentum[j]+=SysState.AngMomentumComp[i][j]; 
	  SysState.CenterOfMass[j]+=SysState.CenterOfMassComp[i][j]; 
	} 
    } 

  for(i=0;i<5;i++) 
    for(j=0;j<3;j++) 
      { 
	if(SysState.MassComp[i]>0) 
	  SysState.CenterOfMassComp[i][j]/=SysState.MassComp[i]; 
      } 
  
  for(j=0;j<3;j++) 
    { 
	  if(SysState.Mass>0) 
	    SysState.CenterOfMass[j]/=SysState.Mass; 
    } 
  
  for(i=0;i<5;i++) 
    { 
      SysState.CenterOfMassComp[i][3]=SysState.MomentumComp[i][3]=SysState.AngMomentumComp[i][3]=0; 
      for(j=0;j<3;j++) 
	{ 
	  SysState.CenterOfMassComp[i][3]+=SysState.CenterOfMassComp[i][j]*SysState.CenterOfMassComp[i][j]; 
	  SysState.MomentumComp[i][3]+=SysState.MomentumComp[i][j]*SysState.MomentumComp[i][j]; 
	  SysState.AngMomentumComp[i][3]+=SysState.AngMomentumComp[i][j]*SysState.AngMomentumComp[i][j]; 
	} 
      SysState.CenterOfMassComp[i][3]=sqrt(SysState.CenterOfMassComp[i][3]); 
      SysState.MomentumComp[i][3]=sqrt(SysState.MomentumComp[i][3]); 
      SysState.AngMomentumComp[i][3]=sqrt(SysState.AngMomentumComp[i][3]); 
    } 
      
  SysState.CenterOfMass[3]=SysState.Momentum[3]=SysState.AngMomentum[3]=0; 
  
  for(j=0;j<3;j++) 
    { 
      SysState.CenterOfMass[3]+=SysState.CenterOfMass[j]*SysState.CenterOfMass[j]; 
      SysState.Momentum[3]+=SysState.Momentum[j]*SysState.Momentum[j]; 
      SysState.AngMomentum[3]+=SysState.AngMomentum[j]*SysState.AngMomentum[j]; 
    } 
  
  SysState.CenterOfMass[3]=sqrt(SysState.CenterOfMass[3]); 
  SysState.Momentum[3]=sqrt(SysState.Momentum[3]); 
  SysState.AngMomentum[3]=sqrt(SysState.AngMomentum[3]); 
}
