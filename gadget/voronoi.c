#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef VORONOI

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

static double AllocFactor_MaxNinlist = 1.5;

int Ndp;			/* number of delaunay points */
int MaxNdp;			/* maximum number of delaunay points */
point *DP;			/* delaunay points */

int Ndt;			/* number of delaunary tetrahedra/triangles */
int MaxNdt;			/* maximum number of delaunary tetrahedra/triangles */
tetra *DT;			/* Delaunay tetrahedra/triangles */


int Nvf;			/* number of Voronoi faces */
int MaxNvf;			/* maximum number of Voronoi faces */
face *VF;			/* Voronoi faces */


#ifdef VORONOI_MESHRELAX
struct grad_data *Grad, *GradExch;
#endif



point *DPinfinity;
double CentralOffsetX, CentralOffsetY, CentralOffsetZ, ConversionFac;

struct list_export_data *ListExports;
struct list_P_data *List_P;
struct primexch *PrimExch;
struct grad_data *Grad, *GradExch;

int CountInSphereTests, CountInSphereTestsExact;
int CountConvexEdgeTest, CountConvexEdgeTestExact;
int Ninlist, MaxNinlist;
int Flag_MaxNinlistReached;

int CountFlips, Count_1_to_3_Flips2d, Count_2_to_4_Flips2d;
int Count_1_to_4_Flips, Count_2_to_3_Flips, Count_3_to_2_Flips, Count_4_to_4_Flips;
int Count_EdgeSplits, Count_FaceSplits;
int Count_InTetra, Count_InTetraExact;
double drand48();

void voronoi_mesh(void)
{
  tetra *tlast;
  int i, n, glob_flag, ntot, skip;
  double timeinsert = 0, tstart, tend;

  if(ThisTask == 0)
    printf("\ncreate delaunay mesh\n");


  do
    {
      for(i = 0; i < N_gas; i++)
	SphP[i].Volume = SphP[i].MaxDelaunayRadius;	/* just for back-up, in case we have to repeat with larger AllocFactor_MaxNinlist */

      initialize_and_create_first_tetra();

      CountInSphereTests = CountInSphereTestsExact = 0;
      CountConvexEdgeTest = CountConvexEdgeTestExact = 0;
      CountFlips = Count_1_to_3_Flips2d = Count_2_to_4_Flips2d = 0;
      Count_1_to_4_Flips = 0;
      Count_2_to_3_Flips = 0;
      Count_3_to_2_Flips = 0;
      Count_4_to_4_Flips = 0;
      Count_EdgeSplits = 0;
      Count_FaceSplits = 0;
      Count_InTetra = Count_InTetraExact = 0;

      Ninlist = 0;
      Flag_MaxNinlistReached = 0;
      MaxNinlist = AllocFactor_MaxNinlist * N_gas;
      ListExports = mymalloc(MaxNinlist * sizeof(struct list_export_data));
      List_P = mymalloc(N_gas * sizeof(struct list_P_data));

      tlast = &DT[0];

      do
	{
	  skip = Ndp;

	  if(Ndp == 0)
	    n = voronoi_get_local_particles();
	  else
	    n = voronoi_get_additional_points();

	  MPI_Allreduce(&Flag_MaxNinlistReached, &glob_flag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	  if(glob_flag)
	    break;

	  MPI_Allreduce(&n, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	  if(ThisTask == 0)
	    printf("have obtained %d additional points (%d on task 0)\n", ntot, n);

	  tstart = second();
	  for(i = 0; i < n; i++)
	    {
	      //      printf("\ninsert=%d\n", i);

	      set_integers_for_point(&DP[skip + i]);

	      tlast = insert_point(&DP[skip + i], tlast);
	    }
	  tend = second();
	  timeinsert += timediff(tstart, tend);

	  if(ThisTask == 0)
	    printf("points inserted.\n");

	  n = compute_max_delaunay_radius();

	  MPI_Allreduce(&n, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	}
      while(ntot > 0);

      if(glob_flag)
	{
	  myfree(List_P);
	  myfree(ListExports);
	  myfree(DT);
	  myfree(DP - 5);
	  myfree(VF);

	  for(i = 0; i < N_gas; i++)
	    SphP[i].MaxDelaunayRadius = SphP[i].Volume;

	  if(AllocFactor_MaxNinlist < 10)
	    {
	      AllocFactor_MaxNinlist *= 1.5;
	      if(ThisTask == 0)
		{
		  printf("\nnot enough space in export list, increasing AllocFactor_MaxNinlist to %g\n\n",
			 AllocFactor_MaxNinlist);
		}
	    }
	  else
	    {
	      printf("\nAllocFactor_MaxNinlist becomes too large, stopping\n");
	      endrun(0);
	    }
	}
    }
  while(glob_flag);


  if(ThisTask == 0)
    {
#ifndef TWODIMS
      printf("Delaunay tetrahedra are calculated, point insertion took=%g sec\n", timeinsert);
      printf("D-Points=%d  D-Tetrahedra=%d  InSphereTests=%d InSphereTestsExact=%d  Flips=%d\n",
	     Ndp, Ndt, CountInSphereTests, CountInSphereTestsExact, CountFlips);
      printf
	("   1_to_4_Flips=%d  2_to_3_Flips=%d  3_to_2_Flips=%d  4_to_4_Flips=%d  FaceSplits=%d  EdgeSplits=%d\n",
	 Count_1_to_4_Flips, Count_2_to_3_Flips, Count_3_to_2_Flips, Count_4_to_4_Flips, Count_FaceSplits,
	 Count_EdgeSplits);
      printf("   InTetra=%d  InTetraExact=%d  ConvexEdgeTest=%d  ConvexEdgeTestExact=%d\n", Count_InTetra,
	     Count_InTetraExact, CountConvexEdgeTest, CountConvexEdgeTestExact);
      printf("\n");
#else
      printf("Delaunay triangles are calculated, point insertion took=%g sec\n", timeinsert);
      printf("D-Points=%d  D-Triangles=%d  InCircleTests=%d InCircleTestsExact=%d  Flips=%d\n",
	     Ndp, Ndt, CountInSphereTests, CountInSphereTestsExact, CountFlips);
      printf("   1_to_3_Flips=%d  2_to_4_Flips=%d  InTriangle=%d  InTriangleExact=%d\n",
	     Count_1_to_3_Flips2d, Count_2_to_4_Flips2d, Count_InTetra, Count_InTetraExact);
      printf("\n");
#endif
    }


  //  dump_points();


  compute_circumcircles();

  for(i = 0; i < N_gas; i++)
    SphP[i].Volume = 0;

  compute_voronoi_faces_and_volumes();

  double vol, voltot;

  for(i = 0, vol = 0; i < N_gas; i++)
    {
      if(SphP[i].Volume > (boxSize_X * boxSize_Y * boxSize_Z))
	printf("i=%d vol=%g  hsml=%g MaxDelaunayRadius=%g\n", i, SphP[i].Volume, SphP[i].Hsml,
	       SphP[i].MaxDelaunayRadius);

      //      printf("i=%d  vol=%g\n", i,  SphP[i].Volume);
      vol += SphP[i].Volume;
    }

  MPI_Allreduce(&vol, &voltot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    printf("total volume = %g\n", voltot);

}



int compute_max_delaunay_radius(void)
{
  int i, j, count;
  point *p, *q;
  double dx, dy, dz, r;

  for(i = 0; i < N_gas; i++)
    SphP[i].MaxDelaunayRadius = 0;

  for(i = 0; i < Ndt; i++)
    {
      if(DT[i].deleted)
	continue;

      for(j = 0; j < (DIMS + 1); j++)
	{
	  p = DT[i].p[j];
	  q = DT[i].t[j]->p[DT[i].s[j]];
	  dx = p->x - q->x;
	  dy = p->y - q->y;
	  dz = p->z - q->z;

	  r = sqrt(dx * dx + dy * dy + dz * dz);

	  if(p->task == ThisTask && p->index < N_gas && p->index >= 0)
	    if(r > SphP[p->index].MaxDelaunayRadius)
	      SphP[p->index].MaxDelaunayRadius = r;

	  if(q->task == ThisTask && q->index < N_gas && q->index >= 0)
	    if(r > SphP[q->index].MaxDelaunayRadius)
	      SphP[q->index].MaxDelaunayRadius = r;
	}
    }

  for(i = 0, count = 0; i < N_gas; i++)
    if(SphP[i].MaxDelaunayRadius > SphP[i].Hsml)
      {
	//  printf("i=%d SphP[i].MaxDelaunayRadius=%g  SphP[i].Hsml=%g\n", i, SphP[i].MaxDelaunayRadius, SphP[i].Hsml);

	count++;

	if(SphP[i].Hsml == boxHalf_X)
	  {
	    printf("too big Hsml: i=%d vol=%g  hsml=%g MaxDelaunayRadius=%g\n", i,
		   SphP[i].Volume, SphP[i].Hsml, SphP[i].MaxDelaunayRadius);
	    endrun(6864132);
	  }

	SphP[i].Hsml *= HSML_INCREASE_FACTOR;

	if(SphP[i].Hsml > boxHalf_X || SphP[i].Hsml > boxHalf_Y || SphP[i].Hsml > boxHalf_Z)
	  {
	    SphP[i].Hsml = boxHalf_X;
	  }

	SphP[i].MaxDelaunayRadius = MAX_REAL_NUMBER;
      }

  return count;
}


void compute_voronoi_faces_and_volumes(void)
{
  int i, bit, nr;

  for(i = 0; i < N_gas; i++)
    {
      SphP[i].Volume = 0;
      SphP[i].Center[0] = 0;
      SphP[i].Center[1] = 0;
      SphP[i].Center[2] = 0;
    }

  for(i = 0; i < Ndt; i++)
    DT[i].egde_visited = 0;


  for(i = 0; i < Ndt; i++)
    {
      if(DT[i].deleted)
	continue;

      bit = 1;
      nr = 0;

      while(DT[i].egde_visited != EDGE_ALL)
	{
	  if((DT[i].egde_visited & bit) == 0)
	    process_edge_faces_and_volumes(&DT[i], nr);

	  bit <<= 1;
	  nr++;
	}
    }

  for(i = 0; i < N_gas; i++)
    {
      if(SphP[i].Volume)
	{
	  SphP[i].Center[0] /= SphP[i].Volume;
	  SphP[i].Center[1] /= SphP[i].Volume;
	  SphP[i].Center[2] /= SphP[i].Volume;
	}
    }
}

void dump_points(void)
{
  FILE *fd;
  int i;
  float xyz[3];

  fd = fopen("points.dat", "w");
  my_fwrite(&Ndp, sizeof(int), 1, fd);
  for(i = 0; i < Ndp; i++)
    {
      xyz[0] = DP[i].x;
      xyz[1] = DP[i].y;
      xyz[2] = DP[i].z;
      my_fwrite(xyz, sizeof(float), 3, fd);
    }
  fclose(fd);
}

#endif
