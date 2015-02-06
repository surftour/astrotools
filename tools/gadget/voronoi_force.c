#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#ifdef VORONOI
#include "voronoi.h"


void voronoi_hydro_force(void)
{
  int i, j, q1, q2, li, ri;
  double length, pressure1, pressure2, cx, cy, cz, c, fac1, fac2, ex, ey, ez, forcex, forcey, forcez;
  double mass1, mass2, dEdt, fviscx, fviscy, fviscz, w, dens1, dens2, entr1, entr2;
  MyFloat *vel1, *vel2;

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
    if(P[i].Type == 0)
			{
      SphP[i].e.DtEntropy = 0;

      for(j = 0; j < 3; j++)
				SphP[i].a.dHydroAccel[j] = 0;
			}
    }


  voronoi_exchange_ghost_variables();

  for(i = 0; i < Nvf; i++)
    {
      q1 = li = VF[i].p1->index;
      q2 = ri = VF[i].p2->index;

      if(!((li >= 0 && li < N_gas && VF[i].p1->task == ThisTask) ||
	   (ri >= 0 && ri < N_gas && VF[i].p2->task == ThisTask)))
	continue;

      if(li >= N_gas && VF[i].p1->task == ThisTask)
	li -= N_gas;

      if(ri >= N_gas && VF[i].p2->task == ThisTask)
	ri -= N_gas;

      if(VF[i].p1->task == ThisTask)
	{
          pressure1 = SphP[li].Pressure;
          dens1 = SphP[li].d.Density;
          vel1 = SphP[li].VelPred;
          mass1 = P[li].Mass;
          entr1 = SphP[li].Entropy;
	}
      else
	{
	  pressure1 = PrimExch[q1].Pressure;
	  dens1 = PrimExch[q1].Density;
	  vel1 = PrimExch[q1].VelPred;
	  mass1 = PrimExch[q1].Mass;
	  entr1 = PrimExch[q1].Entropy;
	}

      if(VF[i].p2->task == ThisTask)
	{
          pressure2 = SphP[ri].Pressure;
          dens2 = SphP[ri].d.Density;
          vel2 = SphP[ri].VelPred;
          mass2 = P[ri].Mass;
          entr2 = SphP[ri].Entropy;
	}
      else
	{
	  pressure2 = PrimExch[q2].Pressure;
	  dens2 = PrimExch[q2].Density;
	  vel2 = PrimExch[q2].VelPred;
	  mass2 = PrimExch[q2].Mass;
	  entr2 = PrimExch[q2].Entropy;
	}

#ifdef VORONOI_MESHRELAX
     if(VF[i].p1->task == ThisTask)
	pressure1 = 1.0 / (P[li].Mass / All.MeanMass);
      else
	pressure1 =
	  1.0 / (PrimExch[li].Mass / All.MeanMass);

      if(VF[i].p2->task == ThisTask)
	pressure2 = 1.0 / (P[ri].Mass / All.MeanMass);
      else
	pressure2 =
	  1.0 / (PrimExch[ri].Mass / All.MeanMass);

      mass1 = mass2 = 1.0;
#endif

      cx = VF[i].p2->x - VF[i].p1->x;
      cy = VF[i].p2->y - VF[i].p1->y;
      cz = VF[i].p2->z - VF[i].p1->z;

      c = sqrt(cx * cx + cy * cy + cz * cz);	/* distance of the two points */
      length = VF[i].area;	/* length/area of common face */



      fac1 = 0.5 * (pressure2 + pressure1) * length / c;
      fac2 = (pressure2 - pressure1) * length / c;

#ifdef VORONOI_MESHRELAX
      fac2 = 0;
#ifdef TWODIMS
      double dist_min = 0.1 * VF[i].area;
#else
      double dist_min = 0.1 * sqrt(VF[i].area);
#endif
      if(c < dist_min)
	{
	  double x = 1.0 - c / dist_min;
          fac1 *= 10 * x * x;
	}
#endif

      ex = VF[i].cx - 0.5 * (VF[i].p1->x + VF[i].p2->x);
      ey = VF[i].cy - 0.5 * (VF[i].p1->y + VF[i].p2->y);
      ez = VF[i].cz - 0.5 * (VF[i].p1->z + VF[i].p2->z);

      forcex = -fac1 * cx - fac2 * ex;
      forcey = -fac1 * cy - fac2 * ey;
      forcez = -fac1 * cz - fac2 * ez;

      /* calculate viscous force */

      w = cx * (vel2[0] - vel1[0]) + cy * (vel2[1] - vel1[1]) + cz * (vel2[2] - vel1[2]);
      if(w < 0)
	{
	  w /= c;

	  /*
	     fviscx = - 0.5 * All.ArtBulkViscConst * mass1 * mass2 / (mass1 + mass2) * cx * w * w / (c * c);
	     fviscy = - 0.5 * All.ArtBulkViscConst * mass1 * mass2 / (mass1 + mass2) * cy * w * w / (c * c);
	     fviscz = - 0.5 * All.ArtBulkViscConst * mass1 * mass2 / (mass1 + mass2) * cz * w * w / (c * c);
	   */

	  double csound, pvisc;

	  /* calculate viscous force */

	  csound = 0.5 * (sqrt(GAMMA * pressure1 / dens1) + sqrt(GAMMA * pressure2 / dens2));
	  pvisc = 0.5 * All.ArtBulkViscConst * (dens1 + dens2) * (-w * csound + 1.5 * w * w);

	  fviscx = -pvisc * length * cx / c;
	  fviscy = -pvisc * length * cy / c;
	  fviscz = -pvisc * length * cz / c;


	  /* rate at which energy is dissipated */
	  dEdt = fviscx * (vel2[0] - vel1[0]) + fviscy * (vel2[1] - vel1[1]) + fviscz * (vel2[2] - vel1[2]);
	  if(dEdt < 0)
	    endrun(88);
	}
      else
	{
	  fviscx = 0;
	  fviscy = 0;
	  fviscz = 0;

	  dEdt = 0;
	}

#ifdef VORONOI_MESHRELAX
      fviscx = fviscy = fviscz = dEdt = 0;
#endif

      if(VF[i].p1->task == ThisTask && q1 >= 0 && q1 < N_gas)
	{
	  if(TimeBinActive[P[q1].TimeBin])
	    {
	      SphP[q1].a.dHydroAccel[0] += (forcex + fviscx) / mass1;
	      SphP[q1].a.dHydroAccel[1] += (forcey + fviscy) / mass1;
	      SphP[q1].a.dHydroAccel[2] += (forcez + fviscz) / mass1;

	      SphP[q1].e.DtEntropy += 0.5 * dEdt * GAMMA_MINUS1 / pow(dens1, GAMMA_MINUS1) / mass1;
	    }
	}

      if(VF[i].p2->task == ThisTask && q2 >= 0 && q2 < N_gas)
	{
	  if(TimeBinActive[P[q2].TimeBin])
	    {
	      SphP[q2].a.dHydroAccel[0] -= (forcex + fviscx) / mass2;
	      SphP[q2].a.dHydroAccel[1] -= (forcey + fviscy) / mass2;
	      SphP[q2].a.dHydroAccel[2] -= (forcez + fviscz) / mass2;

	      SphP[q2].e.DtEntropy += 0.5 * dEdt * GAMMA_MINUS1 / pow(dens2, GAMMA_MINUS1) / mass2;
	    }
	}
    }


#ifdef VORONOI_MESHRELAX
  voronoi_meshrelax();

  myfree(Grad);
  myfree(GradExch);
#endif

  myfree(PrimExch);

  myfree(List_P);
  myfree(ListExports);

  myfree(DT);
  myfree(DP - 5);
  myfree(VF);			/* free the list of faces */
}

#endif
