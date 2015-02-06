#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

#ifdef VORONOI

#include "voronoi.h"


static struct vorodata_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int Origin;
  int NodeList[NODELISTLENGTH];
}
 *VoroDataIn, *VoroDataGet;



struct data_primexch_compare
{
  int rank, task, index;
}
 *SortPrimExch, *SortPrimExch2;



static point *DP_Buffer;
int MaxN_DP_Buffer, N_DP_Buffer;

int NadditionalPoints;


int voronoi_get_additional_points(void)
{
  int i, j, ndone, ndone_flag, dummy;
  int ngrp, sendTask, recvTask, place, nexport, nimport;

  NadditionalPoints = 0;

  /* allocate buffers to arrange communication */

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     2 * sizeof(struct vorodata_in)));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));



  i = 0;

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */

      for(nexport = 0; i < N_gas; i++)
	{

	  /*          if(SphP[i].Hsml > boxHalf_X ||  SphP[i].Hsml > boxHalf_Y || SphP[i].Hsml > boxHalf_Z)
	     {
	     printf("too big Hsml: i=%d vol=%g  hsml=%g MaxDelaunayRadius=%g\n", i, 
	     SphP[i].Volume, SphP[i].Hsml,
	     SphP[i].MaxDelaunayRadius);

	     endrun(68642);
	     }
	   */

	  if(SphP[i].Hsml < SphP[i].MaxDelaunayRadius)
	    {
	      if(voronoi_exchange_evaluate(i, 0, &nexport, Send_count) < 0)
		break;
	    }
	}

#ifdef MYSORT
      mysort_dataindex(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#else
      qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#endif

      MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);


      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
	  nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      VoroDataGet = (struct vorodata_in *) mymalloc(nimport * sizeof(struct vorodata_in));
      VoroDataIn = (struct vorodata_in *) mymalloc(nexport * sizeof(struct vorodata_in));

      /* prepare particle data for export */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  VoroDataIn[j].Pos[0] = P[place].Pos[0];
	  VoroDataIn[j].Pos[1] = P[place].Pos[1];
	  VoroDataIn[j].Pos[2] = P[place].Pos[2];
	  VoroDataIn[j].Origin = ThisTask;
	  VoroDataIn[j].Hsml = PPP[place].Hsml;

	  memcpy(VoroDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	}

      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&VoroDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct vorodata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &VoroDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct vorodata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}

      myfree(VoroDataIn);

      MaxN_DP_Buffer = 4 * N_gas;
      N_DP_Buffer = 0;

      DP_Buffer = (point *) mymalloc(MaxN_DP_Buffer * sizeof(point));

      for(j = 0; j < NTask; j++)
	Send_count[j] = 0;

      /* now do the particles that were sent to us */

      for(j = 0; j < nimport; j++)
	voronoi_exchange_evaluate(j, 1, &dummy, Send_count);


      MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);


      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
	  nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      if(nimport + Ndp > MaxNdp)
	endrun(123);


      /* get the delaunay points */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&DP_Buffer[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(point), MPI_BYTE,
			       recvTask, TAG_DENS_B,
			       &DP[Ndp + Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(point), MPI_BYTE,
			       recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}

      Ndp += nimport;
      NadditionalPoints += nimport;

      myfree(DP_Buffer);
      myfree(VoroDataGet);


      if(i >= N_gas)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }
  while(ndone < NTask);


  myfree(DataNodeList);
  myfree(DataIndexTable);

  return NadditionalPoints;
}



int voronoi_exchange_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int origin;
  int startnode, numngb, numngb_inbox, listindex = 0;
  double h;
  MyDouble *pos;


  if(mode == 0)
    {
      pos = P[target].Pos;
      h = PPP[target].Hsml;
      h *= HSML_INCREASE_FACTOR;
      origin = ThisTask;
    }
  else
    {
      pos = VoroDataGet[target].Pos;
      h = VoroDataGet[target].Hsml;
      origin = VoroDataGet[target].Origin;
    }


  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = VoroDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  numngb = 0;

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb_inbox = ngb_treefind_voronoi(pos, h, target, origin, &startnode, mode, nexport, nsend_local);

	  if(numngb_inbox < 0)
	    return -1;
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = VoroDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  return 0;
}



int ngb_treefind_voronoi(MyDouble searchcenter[3], MyFloat hsml, int target, int origin,
			 int *startnode, int mode, int *nexport, int *nsend_local)
{
  int numngb, no, p, task, nexport_save, ndp_save, nadditionalpoints_save;
  int image_flag;
  struct NODE *current;
  MyDouble dx, dy, dz, dist, dist2;
  struct list_export_data *listp;


  nadditionalpoints_save = NadditionalPoints;
  ndp_save = Ndp;
  nexport_save = *nexport;

  numngb = 0;
  no = *startnode;


  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(P[p].Type > 0)
	    continue;

	  if(P[p].Ti_current != All.Ti_Current)
	    drift_particle(p, All.Ti_Current);

	  dist = hsml;

	  image_flag = 0;	/* for each coordinates there are three possibilities. We
				   encode them to basis three, i.e. x*3^0 + y*3^1 + z*3^2 */

	  dx = P[p].Pos[0] - searchcenter[0];
#if defined(PERIODIC) && !defined(REFLECTIVE_X)
	  if(dx < -boxHalf_X)
	    {
	      dx += boxSize_X;
	      image_flag += 1;
	    }
	  else if(dx > boxHalf_X)
	    {
	      dx -= boxSize_X;
	      image_flag += 2;
	    }
#endif
	  if(fabs(dx) > dist)
	    continue;


	  dy = P[p].Pos[1] - searchcenter[1];
#if defined(PERIODIC) && !defined(REFLECTIVE_Y)
	  if(dy < -boxHalf_Y)
	    {
	      dy += boxSize_Y;
	      image_flag += 1 * 3;
	    }
	  else if(dy > boxHalf_Y)
	    {
	      dy -= boxSize_Y;
	      image_flag += 2 * 3;
	    }
#endif
	  if(fabs(dy) > dist)
	    continue;


	  dz = P[p].Pos[2] - searchcenter[2];
#if defined(PERIODIC) && !defined(REFLECTIVE_Z)
	  if(dz < -boxHalf_Z)
	    {
	      dz += boxSize_Z;
	      image_flag += 1 * 9;
	    }
	  else if(dz > boxHalf_Z)
	    {
	      dz -= boxSize_Z;
	      image_flag += 2 * 9;
	    }
#endif
	  if(fabs(dz) > dist)
	    continue;


	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  /* now we need to check whether this particle has already been sent to 
	     the requesting cpu for this particular image shift */



	  if(List_P[p].firstexport)
	    {
	      if(List_P[p].currentexport->origin != origin)
		{
		  listp = List_P[p].firstexport;
		  while(listp)
		    {
		      if(listp->origin == origin)
			{
			  List_P[p].currentexport = listp;
			  break;
			}
		      if(!listp->nextexport)
			{
			  if(Ninlist < MaxNinlist)
			    {
			      List_P[p].currentexport = &ListExports[Ninlist++];
			      List_P[p].currentexport->image_bits = 0;
			      List_P[p].currentexport->nextexport = NULL;
			      List_P[p].currentexport->origin = origin;
			      List_P[p].currentexport->index = p;
			      listp->nextexport = List_P[p].currentexport;
			      break;
			    }
			  else
			    {
			      Flag_MaxNinlistReached = 1;
			      break;
			    }
			}
		      listp = listp->nextexport;
		    }
		}
	    }
	  else
	    endrun(22);

	  if(Flag_MaxNinlistReached)
	    continue;


	  if(!(List_P[p].currentexport->image_bits & (1 << image_flag)))
	    {
	      List_P[p].currentexport->image_bits |= (1 << image_flag);

	      /* add the particle to the ones that need to be exported */

	      if(origin == ThisTask)
		{
		  if(Ndp >= MaxNdp)
		    endrun(1313);

		  DP[Ndp].x = dx + searchcenter[0];
		  DP[Ndp].y = dy + searchcenter[1];
		  DP[Ndp].z = dz + searchcenter[2];
		  DP[Ndp].task = ThisTask;
		  DP[Ndp].index = p + N_gas;	/* this is a mirrored local point */
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
		  DP[Ndp].image_flags = (1 << image_flag);
#endif
		  Ndp++;
		  NadditionalPoints++;
		}
	      else
		{
		  if(N_DP_Buffer >= MaxN_DP_Buffer)
		    endrun(12319998);

		  DP_Buffer[N_DP_Buffer].x = dx + searchcenter[0];
		  DP_Buffer[N_DP_Buffer].y = dy + searchcenter[1];
		  DP_Buffer[N_DP_Buffer].z = dz + searchcenter[2];
		  DP_Buffer[N_DP_Buffer].task = ThisTask;
		  DP_Buffer[N_DP_Buffer].index = p;
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
		  DP_Buffer[N_DP_Buffer].image_flags = (1 << image_flag);
#endif
		  nsend_local[origin]++;
		  N_DP_Buffer++;
		}
	    }


#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
	  int repx = 0, repy = 0, repz = 0, image_flag_refl;

	  double dxx, dyy, dzz;

#if defined(REFLECTIVE_X)
	  for(repx = -1; repx <= 1; repx++)
#endif
#if defined(REFLECTIVE_Y)
	    for(repy = -1; repy <= 1; repy++)
#endif
#if defined(REFLECTIVE_Z)
	      for(repz = -1; repz <= 1; repz++)
#endif
		{
		  image_flag_refl = image_flag;

		  if(repx == -1)
		    {
		      dxx = -(P[p].Pos[0] + searchcenter[0]);
		      image_flag_refl += 1;
		    }
		  else if(repx == 1)
		    {
		      dxx = 2 * boxSize_X - (P[p].Pos[0] + searchcenter[0]);
		      image_flag_refl += 2;
		    }
		  else
		    dxx = dx;

		  if(repy == -1)
		    {
		      dyy = -(P[p].Pos[1] + searchcenter[1]);
		      image_flag_refl += (1) * 3;
		    }
		  else if(repy == 1)
		    {
		      dyy = 2 * boxSize_Y - (P[p].Pos[1] + searchcenter[1]);
		      image_flag_refl += (2) * 3;
		    }
		  else
		    dyy = dy;

		  if(repz == -1)
		    {
		      dzz = -(P[p].Pos[2] + searchcenter[2]);
		      image_flag_refl += (1) * 9;
		    }
		  else if(repz == 1)
		    {
		      dzz = repz * 2 * boxSize_Z - (P[p].Pos[2] + searchcenter[2]);
		      image_flag_refl += (2) * 9;
		    }
		  else
		    dzz = dz;

		  if(repx != 0 || repy != 0 || repz != 0)
		    if(dxx * dxx + dyy * dyy + dzz * dzz <= dist * dist)
		      {
			if(!(List_P[p].currentexport->image_bits & (1 << image_flag_refl)))
			  {
			    List_P[p].currentexport->image_bits |= (1 << image_flag_refl);

			    /* add the particle to the ones that need to be exported */

			    if(origin == ThisTask)
			      {
				if(Ndp >= MaxNdp)
				  endrun(1313);

				DP[Ndp].x = dxx + searchcenter[0];
				DP[Ndp].y = dyy + searchcenter[1];
				DP[Ndp].z = dzz + searchcenter[2];
				DP[Ndp].task = ThisTask;
				DP[Ndp].index = p + N_gas;	/* this is a mirrored local point */
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
				DP[Ndp].image_flags = (1 << image_flag_refl);
#endif
				Ndp++;
				NadditionalPoints++;
			      }
			    else
			      {
				if(N_DP_Buffer >= MaxN_DP_Buffer)
				  endrun(1255319998);

				DP_Buffer[N_DP_Buffer].x = dxx + searchcenter[0];
				DP_Buffer[N_DP_Buffer].y = dyy + searchcenter[1];
				DP_Buffer[N_DP_Buffer].z = dzz + searchcenter[2];
				DP_Buffer[N_DP_Buffer].task = ThisTask;
				DP_Buffer[N_DP_Buffer].index = p;
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
				DP_Buffer[N_DP_Buffer].image_flags = (1 << image_flag_refl);
#endif
				nsend_local[origin]++;

				N_DP_Buffer++;
			      }
			  }
		      }
		}
#endif
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      if(target >= 0)	/* if no target is given, export will not occur */
		{
		  if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
		    {
		      Exportflag[task] = target;
		      Exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(Exportnodecount[task] == NODELISTLENGTH)
		    {
		      if(*nexport >= All.BunchSize)
			{
			  Ndp = ndp_save;
			  NadditionalPoints = nadditionalpoints_save;
			  *nexport = nexport_save;
			  if(nexport_save == 0)
			    endrun(13004);	/* in this case, the buffer is too small to process even a single particle */
			  for(task = 0; task < NTask; task++)
			    nsend_local[task] = 0;
			  for(no = 0; no < nexport_save; no++)
			    nsend_local[DataIndexTable[no].Task]++;
			  return -1;
			}
		      Exportnodecount[task] = 0;
		      Exportindex[task] = *nexport;
		      DataIndexTable[*nexport].Task = task;
		      DataIndexTable[*nexport].Index = target;
		      DataIndexTable[*nexport].IndexGet = *nexport;
		      *nexport = *nexport + 1;
		      nsend_local[task]++;
		    }

		  DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
		    DomainNodeIndex[no - (All.MaxPart + MaxNodes)];

		  if(Exportnodecount[task] < NODELISTLENGTH)
		    DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
		}

	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
		  return numngb;
		}
	    }

	  if(current->Ti_current != All.Ti_Current)
	    force_drift_node(no, All.Ti_Current);

	  if(!(current->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
	    {
	      if(current->u.d.mass)	/* open cell */
		{
		  no = current->u.d.nextnode;
		  continue;
		}
	    }

	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;

	  dx = current->center[0] - searchcenter[0];
#if defined(PERIODIC) && !defined(REFLECTIVE_X)
	  if(dx < -boxHalf_X)
	    {
	      dx += boxSize_X;
	    }
	  else if(dx > boxHalf_X)
	    {
	      dx -= boxSize_X;
	    }
#endif
	  if((dx = fabs(dx)) > dist)
	    continue;

	  dy = current->center[1] - searchcenter[1];
#if defined(PERIODIC) && !defined(REFLECTIVE_Y)
	  if(dy < -boxHalf_Y)
	    {
	      dy += boxSize_Y;
	    }
	  else if(dy > boxHalf_Y)
	    {
	      dy -= boxSize_Y;
	    }
#endif
	  if((dy = fabs(dy)) > dist)
	    continue;

	  dz = current->center[2] - searchcenter[2];
#if defined(PERIODIC) && !defined(REFLECTIVE_Z)
	  if(dz < -boxHalf_Z)
	    {
	      dz += boxSize_Z;
	    }
	  else if(dz > boxHalf_Z)
	    {
	      dz -= boxSize_Z;
	    }
#endif
	  if((dz = fabs(dz)) > dist)
	    continue;

	  /* now test against the minimal sphere enclosing everything */
	  dist2 = dist + FACT1 * current->len;
	  if(dx * dx + dy * dy + dz * dz > dist2 * dist2)
	    continue;

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
  return numngb;
}




int voronoi_get_local_particles(void)
{
  int p;

  for(p = 0; p < N_gas; p++)
    {
      if(Ninlist >= MaxNinlist)
	endrun(787879);

      List_P[p].currentexport = List_P[p].firstexport = &ListExports[Ninlist++];
      List_P[p].currentexport->image_bits = 1;
      List_P[p].currentexport->nextexport = NULL;
      List_P[p].currentexport->origin = ThisTask;
      List_P[p].currentexport->index = p;

      if(Ndp >= MaxNdp)
	endrun(1313);

      DP[Ndp].x = P[p].Pos[0];
      DP[Ndp].y = P[p].Pos[1];
      DP[Ndp].z = P[p].Pos[2];
      DP[Ndp].task = ThisTask;
      DP[Ndp].index = p;
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
      DP[Ndp].image_flags = 1;
#endif
      Ndp++;
    }

  return N_gas;
}




int voronoi_bitcount(unsigned int x)
{
  int count;

  for(count = 0; x; count++)
    x &= x - 1;
  return count;
}



void voronoi_exchange_ghost_variables(void)
{
  struct list_export_data *listp;
  struct primexch *tmpPrimExch;
  int i, j, p, task, off, count;
  int ngrp, sendTask, recvTask, place, nexport, nimport;

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(p = 0; p < N_gas; p++)
    {
      listp = List_P[p].firstexport;
      while(listp)
	{
	  if(listp->origin != ThisTask)
	    {
	      Send_count[listp->origin]++;
	    }
	  listp = listp->nextexport;
	}
    }

  MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, nexport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
      nimport += Recv_count[j];
      nexport += Send_count[j];

      if(j > 0)
	{
	  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	}
    }

  PrimExch = (struct primexch *) mymalloc(nimport * sizeof(struct primexch));
#ifdef VORONOI_MESHRELAX
  GradExch = (struct grad_data *) mymalloc(nimport * sizeof(struct grad_data));
#endif
  tmpPrimExch = (struct primexch *) mymalloc(nexport * sizeof(struct primexch));



  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(p = 0; p < N_gas; p++)
    {
      listp = List_P[p].firstexport;
      while(listp)
	{
	  if((task = listp->origin) != ThisTask)
	    {
	      place = listp->index;
	      off = Send_offset[task] + Send_count[task]++;

	      tmpPrimExch[off].Pressure = SphP[place].Pressure;
	      tmpPrimExch[off].Mass = P[place].Mass;
	      tmpPrimExch[off].Density = SphP[place].d.Density;
	      tmpPrimExch[off].Entropy = SphP[place].Entropy;

	      for(j = 0; j < 3; j++)
		tmpPrimExch[off].VelPred[j] = SphP[place].VelPred[j];

	      tmpPrimExch[off].task = ThisTask;
	      tmpPrimExch[off].index = place;
	    }
	  listp = listp->nextexport;
	}
    }


  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* get the particles */
	      MPI_Sendrecv(&tmpPrimExch[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct primexch), MPI_BYTE,
			   recvTask, TAG_DENS_A,
			   &PrimExch[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct primexch), MPI_BYTE,
			   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  myfree(tmpPrimExch);



  /* now we need to associate the imported data with the points stored in the DP[] array */


  SortPrimExch = (struct data_primexch_compare *) mymalloc(nimport * sizeof(struct data_primexch_compare));

  for(i = 0; i < nimport; i++)
    {
      SortPrimExch[i].rank = i;
      SortPrimExch[i].task = PrimExch[i].task;
      SortPrimExch[i].index = PrimExch[i].index;
    }

  /* let sort the data according to task and index */
  qsort(SortPrimExch, nimport, sizeof(struct data_primexch_compare), compare_primexch);


  SortPrimExch2 = (struct data_primexch_compare *) mymalloc(Ndp * sizeof(struct data_primexch_compare));

  for(i = 0, count = 0; i < Ndp; i++)
    {
      if(DP[i].task != ThisTask)
	{
	  SortPrimExch2[count].rank = i;
	  SortPrimExch2[count].task = DP[i].task;
	  SortPrimExch2[count].index = DP[i].index;
	  count++;
	}
    }

  /* let sort according to task and index */
  qsort(SortPrimExch2, count, sizeof(struct data_primexch_compare), compare_primexch);


  for(i = 0, j = 0; i < count; i++)
    {
      if(SortPrimExch2[i].task != SortPrimExch[j].task || SortPrimExch2[i].index != SortPrimExch[j].index)
	j++;

      if(j >= nimport)
	endrun(888);

      DP[SortPrimExch2[i].rank].index = SortPrimExch[j].rank;	/* note: this change is now permanent and available for next exchange */
    }

  myfree(SortPrimExch2);
  myfree(SortPrimExch);
}


#ifdef VORONOI_MESHRELAX

void voronoi_exchange_gradients(void)
{
  struct list_export_data *listp;
  struct grad_data *tmpGradExch;
  struct primexch *tmpPrimExch;
  int j, p, task, off;
  int ngrp, sendTask, recvTask, place, nexport, nimport;

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(p = 0; p < N_gas; p++)
    {
      listp = List_P[p].firstexport;
      while(listp)
	{
	  if(listp->origin != ThisTask)
	    {
	      Send_count[listp->origin]++;
	    }
	  listp = listp->nextexport;
	}
    }

  MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, nexport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
      nimport += Recv_count[j];
      nexport += Send_count[j];

      if(j > 0)
	{
	  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	}
    }

  tmpPrimExch = (struct primexch *) mymalloc(nexport * sizeof(struct primexch));
  tmpGradExch = (struct grad_data *) mymalloc(nexport * sizeof(struct grad_data));


  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(p = 0; p < N_gas; p++)
    {
      listp = List_P[p].firstexport;
      while(listp)
	{
	  if((task = listp->origin) != ThisTask)
	    {
	      place = listp->index;
	      off = Send_offset[task] + Send_count[task]++;

	      tmpPrimExch[off].Pressure = SphP[place].Pressure;
	      tmpPrimExch[off].Mass = P[place].Mass;
	      tmpPrimExch[off].Density = SphP[place].d.Density;
	      tmpPrimExch[off].Entropy = SphP[place].Entropy;

	      for(j = 0; j < 3; j++)
		{
		  tmpPrimExch[off].VelPred[j] = SphP[place].VelPred[j];
		  tmpPrimExch[off].HydroAccel[j] = SphP[place].a.HydroAccel[j];
		  tmpPrimExch[off].Center[j] = SphP[place].Center[j];
		}

	      tmpPrimExch[off].task = ThisTask;
	      tmpPrimExch[off].index = place;

	      tmpGradExch[off] = Grad[place];
	    }
	  listp = listp->nextexport;
	}
    }


  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* exchange the data */
	      MPI_Sendrecv(&tmpPrimExch[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct primexch), MPI_BYTE,
			   recvTask, TAG_DENS_A,
			   &PrimExch[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct primexch), MPI_BYTE,
			   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	      MPI_Sendrecv(&tmpGradExch[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct grad_data), MPI_BYTE,
			   recvTask, TAG_HYDRO_A,
			   &GradExch[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct grad_data), MPI_BYTE,
			   recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }


  myfree(tmpGradExch);
  myfree(tmpPrimExch);

  /* note: because the sequence is the same as before, we don't have to do the sorts again */
}
#endif






int compare_primexch(const void *a, const void *b)
{
  if(((struct data_primexch_compare *) a)->task < ((struct data_primexch_compare *) b)->task)
    return -1;

  if(((struct data_primexch_compare *) a)->task > ((struct data_primexch_compare *) b)->task)
    return +1;

  if(((struct data_primexch_compare *) a)->index < ((struct data_primexch_compare *) b)->index)
    return -1;

  if(((struct data_primexch_compare *) a)->index > ((struct data_primexch_compare *) b)->index)
    return +1;

  return 0;
}

#endif
