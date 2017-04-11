#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef VORONOI
#include <gmp.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"




#ifndef TWODIMS			/* will only be compiled in 3D case */


#define INSIDE_EPS   1.0e-7
#define USEDBITS 31
#define DOUBLE_to_VORONOIINT(y)       ((int)(((*((long long *) &y)) & 0xFFFFFFFFFFFFFllu) >> (52 - USEDBITS)))


const int access_triangles[4][3] = {
  {1, 3, 2},
  {0, 2, 3},
  {0, 3, 1},
  {0, 1, 2}
};

const int edge_start[6] = { 0, 0, 0, 1, 1, 2 };
const int edge_end[6] = { 1, 2, 3, 2, 3, 3 };
const int edge_opposite[6] = { 3, 1, 2, 3, 0, 1 };
const int edge_nexttetra[6] = { 2, 3, 1, 0, 2, 0 };





void initialize_and_create_first_tetra(void)
{
  point *p;
  int i, n;


  for(i = 0; i < N_gas; i++)
    {
      SphP[i].Hsml = SphP[i].MaxDelaunayRadius;
      SphP[i].MaxDelaunayRadius = MAX_REAL_NUMBER;
    }

  MaxNdp = 2.0 * N_gas;
  MaxNdt = MaxNdp * 20;
  MaxNvf = MaxNdt;

  Ndp = 0;
  Ndt = 0;
  Nvf = 0;


  VF = mymalloc(MaxNvf * sizeof(face));

  DP = mymalloc((MaxNdp + 5) * sizeof(point));
  DP += 5;

  DT = mymalloc(MaxNdt * sizeof(tetra));
 
 if(ThisTask == 0)
    printf("space for mesh allocated  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

 
  /* construct all encompassing huge tetrahedron */

  double box, tetra_incircle, tetra_sidelength, tetra_height, tetra_face_height;

  box = boxSize_X;
  if(box < boxSize_Y)
    box = boxSize_Y;
  if(box < boxSize_Z)
    box = boxSize_Z;

  tetra_incircle = 1.5 * box;
  tetra_sidelength = tetra_incircle * sqrt(24);
  tetra_height = sqrt(2.0 / 3) * tetra_sidelength;
  tetra_face_height = sqrt(3.0) / 2.0 * tetra_sidelength;

  if(ThisTask == 0)
    printf("side-length of enclosing tetraeder=%g\n", tetra_sidelength);

  /* first, let's make the points */
  DP[-4].x = 0.5 * tetra_sidelength;
  DP[-4].y = -1.0 / 3 * tetra_face_height;
  DP[-4].z = -0.25 * tetra_height;

  DP[-3].x = 0;
  DP[-3].y = 2.0 / 3 * tetra_face_height;
  DP[-3].z = -0.25 * tetra_height;

  DP[-2].x = -0.5 * tetra_sidelength;
  DP[-2].y = -1.0 / 3 * tetra_face_height;
  DP[-2].z = -0.25 * tetra_height;

  DP[-1].x = 0;
  DP[-1].y = 0;
  DP[-1].z = 0.75 * tetra_height;

  for(i = -4; i <= -1; i++)
    {
      DP[i].x += 0.5 * box;
      DP[i].y += 0.5 * box;
      DP[i].z += 0.5 * box;
    }

  for(i = -4, p = &DP[-4]; i < 0; i++, p++)
    {
      p->index = -1;
      p->task = ThisTask;
    }

  /* we also define a neutral element at infinity */
  DPinfinity = &DP[-5];

  DPinfinity->x = 0;
  DPinfinity->y = 0;
  DPinfinity->z = 0;
  DPinfinity->index = -1;
  DPinfinity->task = ThisTask;

  /* now let's make the big tetrahedron */
  DT[0].p[0] = DP - 4;
  DT[0].p[1] = DP - 3;
  DT[0].p[2] = DP - 2;
  DT[0].p[3] = DP - 1;

  /* On the outer faces, we attach tetrahedra with the neutral element as tip.
   * This way we will be able to navigate nicely within the tesselation, 
   * and all tetrahedra have defined neighbouring tetrahedra.
   */

  for(i = 0; i < 4; i++)
    {
      n = i + 1;		/* tetra index */

      DT[0].t[i] = &DT[n];
      DT[0].s[i] = 3;

      DT[n].t[3] = &DT[0];
      DT[n].s[3] = i;
      DT[n].p[3] = DPinfinity;
    }


  DT[1].p[0] = DT[0].p[1];
  DT[1].p[1] = DT[0].p[2];
  DT[1].p[2] = DT[0].p[3];

  DT[2].p[0] = DT[0].p[0];
  DT[2].p[1] = DT[0].p[3];
  DT[2].p[2] = DT[0].p[2];

  DT[3].p[0] = DT[0].p[0];
  DT[3].p[1] = DT[0].p[1];
  DT[3].p[2] = DT[0].p[3];

  DT[4].p[0] = DT[0].p[0];
  DT[4].p[1] = DT[0].p[2];
  DT[4].p[2] = DT[0].p[1];


  DT[1].t[0] = &DT[2];
  DT[2].t[0] = &DT[1];
  DT[1].s[0] = 0;
  DT[2].s[0] = 0;

  DT[1].t[1] = &DT[3];
  DT[3].t[0] = &DT[1];
  DT[1].s[1] = 0;
  DT[3].s[0] = 1;

  DT[1].t[2] = &DT[4];
  DT[4].t[0] = &DT[1];
  DT[1].s[2] = 0;
  DT[4].s[0] = 2;

  DT[2].t[2] = &DT[3];
  DT[3].t[1] = &DT[2];
  DT[2].s[2] = 1;
  DT[3].s[1] = 2;

  DT[2].t[1] = &DT[4];
  DT[4].t[2] = &DT[2];
  DT[2].s[1] = 2;
  DT[4].s[2] = 1;

  DT[3].t[2] = &DT[4];
  DT[4].t[1] = &DT[3];
  DT[3].s[2] = 1;
  DT[4].s[1] = 2;

  Ndt = 5;			/* we'll start out with 5 tetras */

  for(i = 0; i < Ndt; i++)
    DT[i].deleted = 0;

  CentralOffsetX = 0.5 * box - 0.5000001 * tetra_sidelength;
  CentralOffsetY = 0.5 * box - (1.0000001 / 3) * tetra_face_height;
  CentralOffsetZ = 0.5 * box - 0.25000001 * tetra_height;

  ConversionFac = 1.0 / (1.001 * tetra_sidelength);

  for(i = -4; i < 0; i++)
    set_integers_for_point(&DP[i]);
}





void process_edge_faces_and_volumes(tetra * t, int nr)
{
  int i, j, k, l, m, ii, jj, kk, ll, count, nr_next;
  face *f;
  tetra *prev, *next;
  double ax, ay, az;
  double bx, by, bz;
  double cx, cy, cz;
  double nx, ny, nz;
  double sx, sy, sz;
  double hx, hy, hz;
  double darea, dvol, h;

  i = edge_start[nr];
  j = edge_end[nr];
  k = edge_opposite[nr];
  l = edge_nexttetra[nr];

  t->egde_visited |= (1 << nr);

  if(Nvf + 1 >= MaxNvf)
    endrun(1123);

  f = &VF[Nvf++];

  f->area = 0;
  f->p1 = t->p[i];
  f->p2 = t->p[j];

  f->cx = 0;
  f->cy = 0;
  f->cz = 0;


  cx = t->c.x;
  cy = t->c.y;
  cz = t->c.z;

  count = 0;

  prev = t;
  do
    {
      next = prev->t[l];

      if(prev != t && next != t)
	{
	  ax = prev->c.x - cx;
	  ay = prev->c.y - cy;
	  az = prev->c.z - cz;

	  bx = next->c.x - cx;
	  by = next->c.y - cy;
	  bz = next->c.z - cz;

	  nx = ay * bz - az * by;
	  ny = az * bx - ax * bz;
	  nz = ax * by - ay * bx;

	  sx = next->c.x + prev->c.x + cx;
	  sy = next->c.y + prev->c.y + cy;
	  sz = next->c.z + prev->c.z + cz;

	  darea = 0.5 * sqrt(nx * nx + ny * ny + nz * nz);
	  f->area += darea;

	  darea *= (1.0 / 3);

	  f->cx += darea * sx;
	  f->cy += darea * sy;
	  f->cz += darea * sz;
	}

      for(m = 0, ll = ii = jj = -1; m < 4; m++)
	{
	  if(next->p[m] == prev->p[k])
	    ll = m;
	  if(next->p[m] == prev->p[i])
	    ii = m;
	  if(next->p[m] == prev->p[j])
	    jj = m;
	}

      if(ll < 0 || ii < 0 || jj < 0)
	endrun(1231);

      kk = 6 - (ll + ii + jj);


      /* need to determine the edge number to be able to flag it */

      for(nr_next = 0; nr_next < 6; nr_next++)
	if((edge_start[nr_next] == ii && edge_end[nr_next] == jj) ||
	   (edge_start[nr_next] == jj && edge_end[nr_next] == ii))
	  {
	    if((next->egde_visited & (1 << nr_next)) && next != t)
	      endrun(1312377);

	    next->egde_visited |= (1 << nr_next);
	    break;
	  }

      prev = next;
      i = ii;
      l = ll;
      j = jj;
      k = kk;

      count++;

      if(count > 1000)
	endrun(1231);
    }
  while(next != t);

  i = edge_start[nr];
  j = edge_end[nr];

  if(f->area)
    {
      f->cx /= f->area;
      f->cy /= f->area;
      f->cz /= f->area;
    }

  hx = 0.5 * (t->p[i]->x - t->p[j]->x);
  hy = 0.5 * (t->p[i]->y - t->p[j]->y);
  hz = 0.5 * (t->p[i]->z - t->p[j]->z);

  h = sqrt(hx * hx + hy * hy + hz * hz);
  dvol = (1.0 / 3) * f->area * h;

  if(t->p[i]->task == ThisTask && t->p[i]->index >= 0 && t->p[i]->index < N_gas)
    {
      SphP[t->p[i]->index].Volume += dvol;

      /* let's now compute the center-of-mass of the pyramid at the bottom top */
      sx = 0.75 * f->cx + 0.25 * t->p[i]->x;
      sy = 0.75 * f->cy + 0.25 * t->p[i]->y;
      sz = 0.75 * f->cz + 0.25 * t->p[i]->z;

      SphP[t->p[i]->index].Center[0] += dvol * sx;
      SphP[t->p[i]->index].Center[1] += dvol * sy;
      SphP[t->p[i]->index].Center[2] += dvol * sz;
    }


  if(t->p[j]->task == ThisTask && t->p[j]->index >= 0 && t->p[j]->index < N_gas)
    {
      SphP[t->p[j]->index].Volume += dvol;

      /* let's now compute the center-of-mass of the pyramid on top */
      sx = 0.75 * f->cx + 0.25 * t->p[j]->x;
      sy = 0.75 * f->cy + 0.25 * t->p[j]->y;
      sz = 0.75 * f->cz + 0.25 * t->p[j]->z;

      SphP[t->p[j]->index].Center[0] += dvol * sx;
      SphP[t->p[j]->index].Center[1] += dvol * sy;
      SphP[t->p[j]->index].Center[2] += dvol * sz;
    }
}


tetra *insert_point(point * p, tetra * tstart)	/* returns a tetra that (currently) contains the point p */
{
  tetra *t0, *t1, *t2, *t3, *t4;
  tetra *tetra_with_p, *t, *to_check[STACKSIZE_TETRA], *freestack[STACKSIZE_TETRA];
  int n_faces_to_check = 0, nfree_on_stack = 0, moves;
  int tip_index, flag, edgeface_nr;
  int non_convex, convex_edge = 0, i, j;

  /* first, need to do a point location */
  t0 = get_tetra(p, &moves, tstart, &flag, &edgeface_nr);

  tetra_with_p = t0;


  if(flag == 1)			/* that's the normal split of a tetrahedron into 4 */
    {
      if(n_faces_to_check >= STACKSIZE_TETRA - 4)
	endrun(897);

      /* we now need to split this tetrahedron into four  */
      if(nfree_on_stack)
	t1 = freestack[--nfree_on_stack];
      else
	t1 = &DT[Ndt++];

      if(nfree_on_stack)
	t2 = freestack[--nfree_on_stack];
      else
	t2 = &DT[Ndt++];

      if(nfree_on_stack)
	t3 = freestack[--nfree_on_stack];
      else
	t3 = &DT[Ndt++];

      if(Ndt > MaxNdt)
	endrun(12199);

      make_a_1_to_4_flip(p, t0, t1, t2, t3);

      /* now we have a triangulation again - need to check whether there are 
         facets that are not Delaunay */

      /* let's initialize a stack with the facets that we need to check */

      n_faces_to_check = 0;

      to_check[n_faces_to_check++] = t0;
      to_check[n_faces_to_check++] = t1;
      to_check[n_faces_to_check++] = t2;
      to_check[n_faces_to_check++] = t3;
    }


  if(flag == 2)
    {
      //      printf("need to split a face (nr=%d)\n", edgeface_nr);

      /* create four new tetra  */
      if(nfree_on_stack)
	t1 = freestack[--nfree_on_stack];
      else
	t1 = &DT[Ndt++];

      if(nfree_on_stack)
	t2 = freestack[--nfree_on_stack];
      else
	t2 = &DT[Ndt++];

      if(nfree_on_stack)
	t3 = freestack[--nfree_on_stack];
      else
	t3 = &DT[Ndt++];

      if(nfree_on_stack)
	t4 = freestack[--nfree_on_stack];
      else
	t4 = &DT[Ndt++];

      n_faces_to_check = 0;

      to_check[n_faces_to_check++] = t0;
      to_check[n_faces_to_check++] = t0->t[edgeface_nr];
      to_check[n_faces_to_check++] = t1;
      to_check[n_faces_to_check++] = t2;
      to_check[n_faces_to_check++] = t3;
      to_check[n_faces_to_check++] = t4;

      make_a_face_split(t0, edgeface_nr, p, t1, t2, t3, t4);
    }

  if(flag == 3)			/* here we need to split an edge */
    {
      int i, j, k, l, ii, jj, kk, ll, m, count;
      tetra *prev, *next;

      //      printf("need to split an edge (nr=%d)\n", edgeface_nr);

      /* count how many triangles share the edge */
      i = edge_start[edgeface_nr];
      j = edge_end[edgeface_nr];
      k = edge_opposite[edgeface_nr];
      l = edge_nexttetra[edgeface_nr];

      count = 0;
      n_faces_to_check = 0;

      prev = t0;
      do
	{
	  to_check[n_faces_to_check++] = prev;

	  next = prev->t[l];

	  for(m = 0, ll = ii = jj = -1; m < 4; m++)
	    {
	      if(next->p[m] == prev->p[k])
		ll = m;
	      if(next->p[m] == prev->p[i])
		ii = m;
	      if(next->p[m] == prev->p[j])
		jj = m;
	    }

	  if(ll < 0 || ii < 0 || jj < 0)
	    endrun(122231);

	  kk = 6 - (ll + ii + jj);

	  prev = next;
	  i = ii;
	  l = ll;
	  j = jj;
	  k = kk;

	  count++;

	  if(count > 1000)
	    endrun(1231);
	}
      while(next != t0);

      //      printf("there are %d tetras around the edge\n", count);

      tetra **tlist = mymalloc(count * sizeof(tetra *));

      for(i = 0; i < count; i++)
	{
	  if(nfree_on_stack)
	    tlist[i] = freestack[--nfree_on_stack];
	  else
	    {
	      tlist[i] = &DT[Ndt++];
	      if(Ndt > MaxNdt)
		endrun(12191219);
	    }

	  to_check[n_faces_to_check++] = tlist[i];
	}

      make_an_edge_split(t0, edgeface_nr, count, p, tlist);

      myfree(tlist);
    }

  int iter = 0;

  while(n_faces_to_check)
    {
      iter++;
      if(iter > 200000)
	endrun(112312);

      t = to_check[--n_faces_to_check];	/* this is the current tetra to look at. 
					   The facet in question lies opposite to q */
      if(t->deleted)
	continue;

      for(tip_index = 0; tip_index < 4; tip_index++)
	if(t->p[tip_index] == p)
	  break;

      if(tip_index < 4)		/* otherwise the facet has been removed in a 3-2 flip */
	{
	  tetra *q;
	  point *pp;

	  q = t->t[tip_index];	/* tetrahedron that's opposite of ours and shares the facet */
	  pp = q->p[t->s[tip_index]];	/* point that's opposite of the facet in the other tetrahedron */


	  int ret, ret_exact;

	  ret = InSphere_Errorbound(q->p[0], q->p[1], q->p[2], q->p[3], p);
	  CountInSphereTests++;

	  if(ret != 0)
	    ret_exact = ret;
	  else
	    {
	      ret_exact = InSphere_Exact(q->p[0], q->p[1], q->p[2], q->p[3], p);	// let's decide with exact integer arithmetic 
	      CountInSphereTestsExact++;
	    }

	  /*     printf
	     ("ret = %g  t-tetra=%d q-tetra=%d n_faces_to_check=%d  points q=%d %d %d %d  insert point %d   points t=%d %d %d %d\n",
	     ret, t-DT, q - DT, n_faces_to_check, 
	     q->p[0] - DP, q->p[1] - DP, q->p[2] - DP, q->p[3] - DP, p - DP,
	     t->p[0] - DP, t->p[1] - DP, t->p[2] - DP, t->p[3] - DP);
	   */

	  if(ret_exact > 0)	/* facet is illegal, because point lies inside */
	    {
	      /* let's see whether the point lies in the triangle, or on a side, or opposite of one convex edge */

	      non_convex = convex_edge_test(t, tip_index, &convex_edge);

	      //   printf("non_convex=%d convex_edge=%d\n", non_convex,  convex_edge);

	      if(non_convex == 0)	/* we can make a 2-3 flip */
		{
		  tetra *w;

		  if(nfree_on_stack)
		    w = freestack[--nfree_on_stack];
		  else
		    w = &DT[Ndt++];

		  if(Ndt > MaxNdt)
		    endrun(1217799);

		  if(n_faces_to_check >= STACKSIZE_TETRA - 3)
		    endrun(896);

		  make_a_2_to_3_flip(t, tip_index, q, t->s[tip_index], pp, w);

		  to_check[n_faces_to_check++] = t;
		  to_check[n_faces_to_check++] = q;
		  to_check[n_faces_to_check++] = w;
		}
	      else if(non_convex == 1)	/* we might be able to make a 3-2 flip, or we deal with a convex edge on the outer hull */
		{
		  /* test whether the reflex edge is surrounded by exactly three tetrahedra */

		  i = convex_edge + 2;
		  if(i >= 3)
		    i -= 3;
		  i = access_triangles[tip_index][i];

		  for(j = 0; j < 4; j++)
		    if(t->p[i] == q->p[j])
		      break;

		  if(j >= 4)
		    {
		      endrun(1313);
		    }


		  if(t->t[i] == q->t[j])	/* this means there is exactly one tetrahedron between them, i.e. we have found the
						   third partner for the flip */
		    {
		      tetra *w;

		      w = t->t[i];

		      make_a_3_to_2_flip(t, q, w, tip_index, convex_edge, t->s[tip_index]);

		      w->deleted = 1;

		      if(nfree_on_stack < STACKSIZE_TETRA)
			freestack[nfree_on_stack++] = w;
		      else
			endrun(1212999);


		      tetra_with_p = t;
		      if(n_faces_to_check >= STACKSIZE_TETRA - 2)
			endrun(894);

		      to_check[n_faces_to_check++] = t;
		      to_check[n_faces_to_check++] = q;
		    }
		  else
		    {
		      if(t->t[i]->p[t->s[i]] == DPinfinity && q->t[j]->p[q->s[j]] == DPinfinity)
			{
			  printf("convex edge between points=%d %d on outer hull found\n",
				 (int) (t->p[access_triangles[tip_index][convex_edge]] - DP),
				 (int) (t->
					p[access_triangles[tip_index][convex_edge < 2 ? convex_edge + 1 : 0]]
					- DP));

			  endrun(1231);	/* this should not occur since we have embedded the points into a convex big triangle */
			}
		    }
		}
	      else if(non_convex == 2)	/* we might be able to make a 4-4 flip */
		{
		  i = convex_edge + 2;
		  if(i >= 3)
		    i -= 3;
		  i = access_triangles[tip_index][i];	/* this is the point opposite of edge (but not tip) */

		  for(j = 0; j < 4; j++)
		    if(t->p[i] == q->p[j])
		      break;



		  /* count how many tetras there are around the edge */
		  /*
		     {
		     int i, j, k, l, m, count;
		     int ii, jj, kk, ll;
		     tetra *prev, *next;

		     i = access_triangles[tip_index][convex_edge];
		     j = convex_edge+1;
		     if(j>=3)
		     j-=3;
		     j = access_triangles[tip_index][j];
		     l = convex_edge+2;
		     if(l>=3)
		     l-=3;
		     l = access_triangles[tip_index][l];

		     k = tip_index;

		     count = 0;

		     prev = t;
		     do
		     {
		     printf("tetra=%d points=%d %d %d %d\n",
		     prev-DT, prev->p[0]-DP, prev->p[1]-DP, prev->p[2]-DP, prev->p[3]-DP);
		     next = prev->t[l];

		     for(m = 0, ll = ii = jj = -1; m < 4; m++)
		     {
		     if(next->p[m] == prev->p[k])
		     ll = m;
		     if(next->p[m] == prev->p[i])
		     ii = m;
		     if(next->p[m] == prev->p[j])
		     jj = m;
		     }

		     if(ll < 0 || ii < 0 || jj < 0)
		     endrun(122231);

		     kk = 6 - (ll + ii + jj);

		     prev = next;
		     i = ii;
		     l = ll;
		     j = jj;
		     k = kk;

		     count++;

		     if(count > 1000)
		     endrun(1231);
		     }
		     while(next != t);

		     printf("Count around edge=%d\n", count);
		     }
		   */


		  if(t->t[i]->p[t->s[i]] == q->t[j]->p[q->s[j]])
		    {
		      /* ok, so we really have 4 tetra. The opposite points match up */
		      //                      printf("we should do a 4-4 flip\n");

		      to_check[n_faces_to_check++] = t;
		      to_check[n_faces_to_check++] = q;
		      to_check[n_faces_to_check++] = t->t[i];
		      to_check[n_faces_to_check++] = q->t[j];

		      make_a_4_to_4_flip(t, tip_index, convex_edge);
		    }
		  else
		    {
		      /*
		         printf("t->t[i]->p[t->s[i]]=%d   q->t[j]->p[q->s[j]] =%d i=%d j=%d\n",
		         t->t[i]->p[t->s[i]] -DP,  q->t[j]->p[q->s[j]]-DP, i, j);


		         printf("can't make a 4-to-4 flip\n");
		       */
		    }
		}
	    }
	  else
	    tetra_with_p = t;
	}
    }

  return tetra_with_p;
}


int convex_edge_test(tetra * t, int tip, int *edgenr)
{
  int i0, i1, i2, i3;
  point *p0, *p1, *p2, *p3, *p4;
  int vol, flag0, flag1, flag2;
  int count_zeros = 0;

  i0 = access_triangles[tip][0];
  i1 = access_triangles[tip][1];
  i2 = access_triangles[tip][2];
  i3 = tip;

  p0 = t->p[i0];
  p1 = t->p[i1];
  p2 = t->p[i2];
  p3 = t->p[i3];
  p4 = t->t[i3]->p[t->s[i3]];


  CountConvexEdgeTest++;

  double ax = p1->xx - p0->xx;
  double ay = p1->yy - p0->yy;
  double az = p1->zz - p0->zz;

  double bx = p2->xx - p0->xx;
  double by = p2->yy - p0->yy;
  double bz = p2->zz - p0->zz;

  double cx = p3->xx - p0->xx;
  double cy = p3->yy - p0->yy;
  double cz = p3->zz - p0->zz;

  double qx = p4->xx - p0->xx;
  double qy = p4->yy - p0->yy;
  double qz = p4->zz - p0->zz;


  double mv_data[] = { ax, bx, cx, qx, ay, by, cy, qy, az, bz, cz, qz };
  double x[3];



  int status;

  status = solve_linear_equations(mv_data, x);

  if(status < 0)
    {
      /*
         vol = Orient3d_Exact(t->p[0], t->p[1], t->p[2], t->p[3]);

         if(vol <= 0)
         {
         printf("vol = %d\n", vol);
         printf("a= %g %g %g\n", ax, ay, az);
         printf("b= %g %g %g\n", bx, by, bz);
         printf("c= %g %g %g\n", cx, cy, cz);
         printf("q= %g %g %g\n", qx, qy, qz);
         printf("(axb)*c) = %g\n",
         (ay * bz - az * by) * cx + (az * bx - ax * bz) * cy + (ax * by - ay * bx) * cz);

         endrun(1312);
         }
       */
    }


  /* x now contains the coordinates of the point p4 expanded in the basis (a,b,c) */
  /* the coordinates of point 3 in this basis are (0,0,1) */


  if(status >= 0)
    {
      if(fabs(1.0 - x[2]) < INSIDE_EPS)
	endrun(1312376543);

      double u, v, w;

      w = 1.0 / (1.0 - x[2]);

      u = w * x[0];
      v = w * x[1];

      //  printf("u=%g  v=%g   u+v=%g  1-(u+v)=%g\n", u, v, u + v, 1 - (u+v));


      if(u > INSIDE_EPS && v > INSIDE_EPS && (1 - (u + v)) > INSIDE_EPS)
	{
	  /* we have a point safely in the triangle: 2-3 flip should be fine */
	  return 0;
	}

      if(u > INSIDE_EPS && v < -INSIDE_EPS && (1 - (u + v)) > INSIDE_EPS)
	{
	  /* edge 0 is clearly reflect,  3-2 flip allowed around edge 0 */
	  *edgenr = 0;
	  return 1;
	}

      if(u > INSIDE_EPS && v > INSIDE_EPS && (1 - (u + v)) < -INSIDE_EPS)
	{
	  // printf("3-2 flip allowed since edge 1 is reflex\n");
	  *edgenr = 1;
	  return 1;
	}

      if(u < -INSIDE_EPS && v > INSIDE_EPS && (1 - (u + v)) > INSIDE_EPS)
	{
	  // printf("3-2 flip allowed since edge 2 is reflex\n");
	  *edgenr = 2;
	  return 1;
	}

      if(u < -INSIDE_EPS && v < -INSIDE_EPS && (1 - (u + v)) > INSIDE_EPS)
	return -1;		/* two reflex edges */

      if(u < -INSIDE_EPS && v > INSIDE_EPS && (1 - (u + v)) < -INSIDE_EPS)
	return -1;		/* two reflex edges */

      if(u > INSIDE_EPS && v < -INSIDE_EPS && (1 - (u + v)) < -INSIDE_EPS)
	return -1;		/* two reflex edges */

    }

  CountConvexEdgeTestExact++;

  /* Now we need to test in more detail if we are on one of the edges */

  vol = Orient3d_Exact(p0, p1, p2, p3);

  if(vol <= 0)
    {
      printf("flat or negatively tetrahedron found (vol=%d)\n", vol);
      {
	dump_points();
	endrun(13777123);
      }
    }

  flag0 = Orient3d_Exact(p1, p3, p2, p4);
  flag1 = Orient3d_Exact(p0, p2, p3, p4);
  flag2 = Orient3d_Exact(p0, p3, p1, p4);

  if(flag0 == 0)
    count_zeros++;

  if(flag1 == 0)
    count_zeros++;

  if(flag2 == 0)
    count_zeros++;

  if(flag0 >= 0 && flag1 >= 0 && flag2 < 0)
    {
      //  printf("3-2 flip allowed since edge 0 is reflex\n");
      *edgenr = 0;
      return 1;
    }

  if(flag0 < 0 && flag1 >= 0 && flag2 >= 0)
    {
      // printf("3-2 flip allowed since edge 1 is reflex\n");
      *edgenr = 1;
      return 1;
    }

  if(flag0 >= 0 && flag1 < 0 && flag2 >= 0)
    {
      // printf("3-2 flip allowed since edge 2 is reflex\n");
      *edgenr = 2;
      return 1;
    }


  if(flag0 >= 0 && flag1 >= 0 && flag2 == 0)
    {
      // printf("4-4 flip around edge 0 may be possible\n");
      *edgenr = 0;
      return 2;
    }

  if(flag0 >= 0 && flag1 == 0 && flag2 >= 0)
    {
      // printf("4-4 flip around edge 2 may be possible\n");
      *edgenr = 2;
      return 2;
    }

  if(flag0 == 0 && flag1 >= 0 && flag2 >= 0)
    {
      // printf("4-4 flip around edge 1 may be possible\n");
      *edgenr = 1;
      return 2;
    }

  if(flag0 >= 0 && flag1 >= 0 && flag2 >= 0)
    {
      /* we seem to have a point in the triangle: 2-3 flip should be fine */
      return 0;
    }

  return -1;
}


void make_a_face_split(tetra * t0, int face_nr, point * p, tetra * t1, tetra * t2, tetra * q1, tetra * q2)
{
  tetra *q0;
  int m, i0 = -1, i1 = -1, i2 = -1, i3 = -1, j0 = -1, j1 = -1, j2 = -1, j3 = -1;

  Count_FaceSplits++;
  CountFlips++;

  //  printf("face split of face_nr=%d\n", face_nr);

  q0 = t0->t[face_nr];



  *t1 = *t0;
  *t2 = *t0;

  *q1 = *q0;
  *q2 = *q0;

  t1->deleted = 0;
  t2->deleted = 0;
  q1->deleted = 0;
  q2->deleted = 0;


  i3 = face_nr;
  j3 = t0->s[face_nr];

  switch (i3)
    {
    case 3:
      i0 = 0;
      i1 = 1;
      i2 = 2;
      break;
    case 2:
      i0 = 0;
      i1 = 3;
      i2 = 1;
      break;
    case 1:
      i0 = 0;
      i1 = 2;
      i2 = 3;
      break;
    case 0:
      i0 = 1;
      i1 = 3;
      i2 = 2;
      break;
    }

  for(m = 0; m < 4; m++)
    {
      if(q0->p[m] == t0->p[i0])
	j0 = m;
      if(q0->p[m] == t0->p[i1])
	j2 = m;
      if(q0->p[m] == t0->p[i2])
	j1 = m;
    }

  if(i0 < 0 || i1 < 0 || i2 < 0 || i3 < 0 || j0 < 0 || j1 < 0 || j2 < 0 || j3 < 0)
    endrun(12121999);


  t0->p[i2] = p;
  t1->p[i0] = p;
  t2->p[i1] = p;

  q0->p[j1] = p;
  q1->p[j0] = p;
  q2->p[j2] = p;


  t0->t[i0] = t1;
  t1->t[i2] = t0;
  t0->s[i0] = i2;
  t1->s[i2] = i0;

  t1->t[i1] = t2;
  t2->t[i0] = t1;
  t1->s[i1] = i0;
  t2->s[i0] = i1;

  t2->t[i2] = t0;
  t0->t[i1] = t2;
  t2->s[i2] = i1;
  t0->s[i1] = i2;


  q0->t[j0] = q1;
  q1->t[j1] = q0;
  q0->s[j0] = j1;
  q1->s[j1] = j0;

  q1->t[j2] = q2;
  q2->t[j0] = q1;
  q1->s[j2] = j0;
  q2->s[j0] = j2;

  q2->t[j1] = q0;
  q0->t[j2] = q2;
  q2->s[j1] = j2;
  q0->s[j2] = j1;


  t0->t[i3] = q0;
  q0->t[j3] = t0;
  t0->s[i3] = j3;
  q0->s[j3] = i3;

  t1->t[i3] = q1;
  q1->t[j3] = t1;
  t1->s[i3] = j3;
  q1->s[j3] = i3;

  t2->t[i3] = q2;
  q2->t[j3] = t2;
  t2->s[i3] = j3;
  q2->s[j3] = i3;


  t0->t[i2]->t[t0->s[i2]] = t0;
  t1->t[i0]->t[t1->s[i0]] = t1;
  t2->t[i1]->t[t2->s[i1]] = t2;

  q0->t[j1]->t[q0->s[j1]] = q0;
  q1->t[j0]->t[q1->s[j0]] = q1;
  q2->t[j2]->t[q2->s[j2]] = q2;

}


void make_an_edge_split(tetra * t0, int edge_nr, int count, point * p, tetra ** tlist)
{
  int i, j, k, l, ii, jj, kk, ll, m, nr, nrm, nrp;
  tetra *prev, *next;
  tetra **t_orig_list;
  int *i_list, *j_list, *k_list, *l_list;


  Count_EdgeSplits++;
  CountFlips++;


  //  printf("edge split\n");

  t_orig_list = mymalloc(sizeof(tetra *) * count);
  i_list = mymalloc(sizeof(int) * count);
  j_list = mymalloc(sizeof(int) * count);
  k_list = mymalloc(sizeof(int) * count);
  l_list = mymalloc(sizeof(int) * count);


  i = edge_start[edge_nr];
  j = edge_end[edge_nr];
  k = edge_opposite[edge_nr];
  l = edge_nexttetra[edge_nr];

  nr = 0;
  prev = t0;
  do
    {
      t_orig_list[nr] = prev;
      i_list[nr] = i;
      j_list[nr] = j;
      k_list[nr] = k;
      l_list[nr] = l;

      next = prev->t[l];

      for(m = 0, ll = ii = jj = -1; m < 4; m++)
	{
	  if(next->p[m] == prev->p[k])
	    ll = m;
	  if(next->p[m] == prev->p[i])
	    ii = m;
	  if(next->p[m] == prev->p[j])
	    jj = m;
	}

      if(ll < 0 || ii < 0 || jj < 0)
	endrun(122231);

      kk = 6 - (ll + ii + jj);

      prev = next;
      i = ii;
      l = ll;
      j = jj;
      k = kk;

      nr++;
    }
  while(next != t0);


  for(nr = 0; nr < count; nr++)
    {
      *tlist[nr] = *t_orig_list[nr];

      tlist[nr]->deleted = 0;

      t_orig_list[nr]->p[j_list[nr]] = p;
      tlist[nr]->p[i_list[nr]] = p;


      t_orig_list[nr]->t[i_list[nr]] = tlist[nr];
      tlist[nr]->t[j_list[nr]] = t_orig_list[nr];

      t_orig_list[nr]->s[i_list[nr]] = j_list[nr];
      tlist[nr]->s[j_list[nr]] = i_list[nr];

      tlist[nr]->t[i_list[nr]]->t[tlist[nr]->s[i_list[nr]]] = tlist[nr];

      nrp = nr + 1;
      if(nrp >= count)
	nrp -= count;

      nrm = nr - 1;
      if(nrm < 0)
	nrm += count;

      tlist[nr]->t[l_list[nr]] = tlist[nrp];
      tlist[nr]->s[l_list[nr]] = k_list[nrp];

      tlist[nr]->t[k_list[nr]] = tlist[nrm];
      tlist[nr]->s[k_list[nr]] = l_list[nrm];
    }

  myfree(l_list);
  myfree(k_list);
  myfree(j_list);
  myfree(i_list);

  myfree(t_orig_list);

}

void make_a_4_to_4_flip(tetra * t, int tip_index, int edge_nr)
{
  //  printf("4-to-4 flip\n");

  int i0, i1, i2, j;

  tetra *w, *q, *u;
  tetra *t_top[4], *t_bottom[4];
  int s_top[4], s_bottom[4];
  point *p[6];


  Count_4_to_4_Flips++;
  CountFlips++;


  u = NULL;

  for(j = 0; j < 4; j++)
    {
      t_top[j] = NULL;
      t_bottom[j] = NULL;
      s_top[j] = -1;
      s_bottom[j] = -1;
    }

  for(j = 0; j < 6; j++)
    {
      p[j] = NULL;
    }

  //  printf("edge_nr=%d\n", edge_nr);


  i0 = access_triangles[tip_index][edge_nr];
  edge_nr += 1;
  if(edge_nr >= 3)
    edge_nr -= 3;
  i1 = access_triangles[tip_index][edge_nr];
  edge_nr += 1;
  if(edge_nr >= 3)
    edge_nr -= 3;
  i2 = access_triangles[tip_index][edge_nr];

  t_top[0] = t->t[i0];
  s_top[0] = t->s[i0];

  t_top[1] = t->t[i1];
  s_top[1] = t->s[i1];

  w = t->t[i2];
  q = t->t[tip_index];

  //  printf("i0=%d i1=%d i2=%d tip_index=%d   points=%d %d %d\n", i0, i1, i2, tip_index, t->p[i0] - DP, t->p[i1] - DP, t->p[i2] - DP);

  for(j = 0; j < 4; j++)
    {
      if(w->p[j] == t->p[i0])
	{
	  t_top[3] = w->t[j];
	  s_top[3] = w->s[j];
	}

      if(w->p[j] == t->p[i1])
	{
	  t_top[2] = w->t[j];
	  s_top[2] = w->s[j];
	}

      if(w->p[j] == t->p[tip_index])
	{
	  u = w->t[j];
	}
    }

  //  printf("w(%d)-tetra=%d %d %d %d\n", w - DT, w->t[0] - DT, w->t[1] - DT, w->t[2] - DT, w->t[3] - DT);
  // printf("t(%d)-tetra=%d %d %d %d\n", t - DT, t->t[0] - DT, t->t[1] - DT, t->t[2] - DT, t->t[3] - DT);
  // printf("q(%d)-tetra=%d %d %d %d\n", q - DT, q->t[0] - DT, q->t[1] - DT, q->t[2] - DT, q->t[3] - DT);
  // printf("u(%d)-tetra=%d %d %d %d\n", u - DT, u->t[0] - DT, u->t[1] - DT, u->t[2] - DT, u->t[3] - DT);


  for(j = 0; j < 4; j++)
    {
      if(u->p[j] == t->p[i0])
	{
	  t_bottom[3] = u->t[j];
	  s_bottom[3] = u->s[j];
	}

      if(u->p[j] == t->p[i1])
	{
	  t_bottom[2] = u->t[j];
	  s_bottom[2] = u->s[j];
	}

      if(q->p[j] == t->p[i0])
	{
	  t_bottom[0] = q->t[j];
	  s_bottom[0] = q->s[j];
	}

      if(q->p[j] == t->p[i1])
	{
	  t_bottom[1] = q->t[j];
	  s_bottom[1] = q->s[j];
	}
    }

  //  printf("t=%d w=%d q=%d u=%d \n", t - DT, w - DT, q - DT, u - DT);




  p[0] = t->p[i1];
  p[1] = t->p[i2];
  p[2] = t->p[i0];
  p[3] = t->t[i2]->p[t->s[i2]];
  p[4] = t->p[tip_index];
  p[5] = t->t[tip_index]->p[t->s[tip_index]];


  for(j = 0; j < 4; j++)
    {
      if(t_top[j] == NULL || t_bottom[j] == NULL)
	{
	  printf("bad!\n");
	  endrun(121887123);
	}

      // printf("top=%3d bottom=%3d\n", t_top[j] - DT, t_bottom[j] - DT);
    }


  for(j = 0; j < 4; j++)
    {
      if(t_top[j] == NULL || t_bottom[j] == NULL)
	{
	  printf("bad!\n");
	  endrun(1312398541);
	}

      // printf("stop=%3d sbottom=%3d\n", s_top[j], s_bottom[j]);
    }



  for(j = 0; j < 6; j++)
    {
      if(p[j] == NULL)
	{
	  printf("wrong!\n");
	  endrun(1312377);
	}

      // printf("point=%d\n", p[j] - DP);
    }



  t->p[0] = p[0];
  t->p[1] = p[1];
  t->p[2] = p[5];
  t->p[3] = p[4];

  q->p[0] = p[1];
  q->p[1] = p[2];
  q->p[2] = p[5];
  q->p[3] = p[4];

  u->p[0] = p[2];
  u->p[1] = p[3];
  u->p[2] = p[5];
  u->p[3] = p[4];

  w->p[0] = p[3];
  w->p[1] = p[0];
  w->p[2] = p[5];
  w->p[3] = p[4];


  t->t[0] = q;
  q->t[1] = t;
  t->s[0] = 1;
  q->s[1] = 0;

  q->t[0] = u;
  u->t[1] = q;
  q->s[0] = 1;
  u->s[1] = 0;

  u->t[0] = w;
  w->t[1] = u;
  u->s[0] = 1;
  w->s[1] = 0;

  w->t[0] = t;
  t->t[1] = w;
  w->s[0] = 1;
  t->s[1] = 0;


  t->t[2] = t_top[0];
  t->s[2] = s_top[0];
  t->t[2]->t[t->s[2]] = t;
  t->t[2]->s[t->s[2]] = 2;

  t->t[3] = t_bottom[0];
  t->s[3] = s_bottom[0];
  t->t[3]->t[t->s[3]] = t;
  t->t[3]->s[t->s[3]] = 3;


  q->t[2] = t_top[1];
  q->s[2] = s_top[1];
  q->t[2]->t[q->s[2]] = q;
  q->t[2]->s[q->s[2]] = 2;

  q->t[3] = t_bottom[1];
  q->s[3] = s_bottom[1];
  q->t[3]->t[q->s[3]] = q;
  q->t[3]->s[q->s[3]] = 3;


  u->t[2] = t_top[2];
  u->s[2] = s_top[2];
  u->t[2]->t[u->s[2]] = u;
  u->t[2]->s[u->s[2]] = 2;

  u->t[3] = t_bottom[2];
  u->s[3] = s_bottom[2];
  u->t[3]->t[u->s[3]] = u;
  u->t[3]->s[u->s[3]] = 3;


  w->t[2] = t_top[3];
  w->s[2] = s_top[3];
  w->t[2]->t[w->s[2]] = w;
  w->t[2]->s[w->s[2]] = 2;

  w->t[3] = t_bottom[3];
  w->s[3] = s_bottom[3];
  w->t[3]->t[w->s[3]] = w;
  w->t[3]->s[w->s[3]] = 3;
}



void make_a_1_to_4_flip(point * p, tetra * t0, tetra * t1, tetra * t2, tetra * t3)
{
  //  printf("1-to-4 flip\n");

  Count_1_to_4_Flips++;
  CountFlips++;

  t1->deleted = 0;
  t2->deleted = 0;
  t3->deleted = 0;

  *t1 = *t0;
  *t2 = *t0;
  *t3 = *t0;

  t0->p[0] = p;
  t1->p[1] = p;
  t2->p[2] = p;
  t3->p[3] = p;

  t0->t[1] = t1;
  t1->t[0] = t0;
  t0->s[1] = 0;
  t1->s[0] = 1;

  t1->t[2] = t2;
  t2->t[1] = t1;
  t1->s[2] = 1;
  t2->s[1] = 2;

  t2->t[0] = t0;
  t0->t[2] = t2;
  t2->s[0] = 2;
  t0->s[2] = 0;

  t0->t[3] = t3;
  t3->t[0] = t0;
  t0->s[3] = 0;
  t3->s[0] = 3;

  t1->t[3] = t3;
  t3->t[1] = t1;
  t1->s[3] = 1;
  t3->s[1] = 3;

  t2->t[3] = t3;
  t3->t[2] = t2;
  t2->s[3] = 2;
  t3->s[2] = 3;

  t0->t[0]->t[t0->s[0]] = t0;
  t1->t[1]->t[t1->s[1]] = t1;
  t2->t[2]->t[t2->s[2]] = t2;
  t3->t[3]->t[t3->s[3]] = t3;

}



void make_a_3_to_2_flip(tetra * t0, tetra * t1, tetra * t2, int tip, int edge, int bottom)
{
  int i, j, k, ii, jj, iii, jjj;
  tetra qbak, tbak, wbak;

  Count_3_to_2_Flips++;
  CountFlips++;

  //  printf("3-to-2 flip\n");


  tbak = *t0;
  qbak = *t1;
  wbak = *t2;

  i = edge;
  j = i + 1;
  k = i + 2;
  if(j >= 3)
    j -= 3;
  if(k >= 3)
    k -= 3;

  i = access_triangles[tip][i];
  j = access_triangles[tip][j];
  k = access_triangles[tip][k];

  for(ii = 0; ii < 4; ii++)
    if(tbak.p[i] == qbak.p[ii])
      break;

  for(iii = 0; iii < 4; iii++)
    if(tbak.p[i] == wbak.p[iii])
      break;

  for(jj = 0; jj < 4; jj++)
    if(tbak.p[j] == qbak.p[jj])
      break;

  for(jjj = 0; jjj < 4; jjj++)
    if(tbak.p[j] == wbak.p[jjj])
      break;

  t0->p[0] = qbak.p[bottom];
  t0->p[1] = tbak.p[k];
  t0->p[2] = tbak.p[i];
  t0->p[3] = tbak.p[tip];

  t1->p[0] = qbak.p[bottom];
  t1->p[1] = tbak.p[j];
  t1->p[2] = tbak.p[k];
  t1->p[3] = tbak.p[tip];

  t0->t[2] = t1;
  t1->t[1] = t0;
  t0->s[2] = 1;
  t1->s[1] = 2;

  t0->t[0] = tbak.t[j];
  t0->s[0] = tbak.s[j];
  t0->t[0]->s[t0->s[0]] = 0;
  t0->t[0]->t[t0->s[0]] = t0;

  t0->t[3] = qbak.t[jj];
  t0->s[3] = qbak.s[jj];
  t0->t[3]->s[t0->s[3]] = 3;
  t0->t[3]->t[t0->s[3]] = t0;

  t0->t[1] = wbak.t[jjj];
  t0->s[1] = wbak.s[jjj];
  t0->t[1]->s[t0->s[1]] = 1;
  t0->t[1]->t[t0->s[1]] = t0;


  t1->t[0] = tbak.t[i];
  t1->s[0] = tbak.s[i];
  t1->t[0]->s[t1->s[0]] = 0;
  t1->t[0]->t[t1->s[0]] = t1;

  t1->t[3] = qbak.t[ii];
  t1->s[3] = qbak.s[ii];
  t1->t[3]->s[t1->s[3]] = 3;
  t1->t[3]->t[t1->s[3]] = t1;

  t1->t[2] = wbak.t[iii];
  t1->s[2] = wbak.s[iii];
  t1->t[2]->s[t1->s[2]] = 2;
  t1->t[2]->t[t1->s[2]] = t1;

  CountFlips++;
}


void make_a_2_to_3_flip(tetra * t0, int tip, tetra * t1, int bottom, point * qq, tetra * t2)
{
  tetra qbak, tbak;
  int k;

  Count_2_to_3_Flips++;

  //  printf("2-to-3 flip\n");


  t2->deleted = 0;

  tbak = *t0;
  qbak = *t1;			/* to save info */

  *t1 = *t0;
  *t2 = *t0;

  /* redefine points */
  t0->p[access_triangles[tip][0]] = qq;
  t1->p[access_triangles[tip][1]] = qq;
  t2->p[access_triangles[tip][2]] = qq;

  /* make neighbour connections */
  t0->t[access_triangles[tip][1]] = t1;
  t1->t[access_triangles[tip][0]] = t0;
  t0->s[access_triangles[tip][1]] = access_triangles[tip][0];
  t1->s[access_triangles[tip][0]] = access_triangles[tip][1];

  t0->t[access_triangles[tip][2]] = t2;
  t2->t[access_triangles[tip][0]] = t0;
  t0->s[access_triangles[tip][2]] = access_triangles[tip][0];
  t2->s[access_triangles[tip][0]] = access_triangles[tip][2];

  t1->t[access_triangles[tip][2]] = t2;
  t2->t[access_triangles[tip][1]] = t1;
  t1->s[access_triangles[tip][2]] = access_triangles[tip][1];
  t2->s[access_triangles[tip][1]] = access_triangles[tip][2];

  /* these are the ones on the top */
  t0->t[access_triangles[tip][0]]->t[t0->s[access_triangles[tip][0]]] = t0;
  t1->t[access_triangles[tip][1]]->t[t1->s[access_triangles[tip][1]]] = t1;
  t2->t[access_triangles[tip][2]]->t[t2->s[access_triangles[tip][2]]] = t2;

  /* now the one at the bottom */

  if(qbak.p[access_triangles[bottom][0]] == tbak.p[access_triangles[tip][0]])
    k = 0;
  else if(qbak.p[access_triangles[bottom][1]] == tbak.p[access_triangles[tip][0]])
    k = 1;
  else
    k = 2;

  t0->t[tip] = qbak.t[access_triangles[bottom][k]];
  t0->s[tip] = qbak.s[access_triangles[bottom][k]];
  t0->t[tip]->t[t0->s[tip]] = t0;
  t0->t[tip]->s[t0->s[tip]] = tip;


  if(qbak.p[access_triangles[bottom][0]] == tbak.p[access_triangles[tip][1]])
    k = 0;
  else if(qbak.p[access_triangles[bottom][1]] == tbak.p[access_triangles[tip][1]])
    k = 1;
  else
    k = 2;

  t1->t[tip] = qbak.t[access_triangles[bottom][k]];
  t1->s[tip] = qbak.s[access_triangles[bottom][k]];
  t1->t[tip]->t[t1->s[tip]] = t1;
  t1->t[tip]->s[t1->s[tip]] = tip;


  if(qbak.p[access_triangles[bottom][0]] == tbak.p[access_triangles[tip][2]])
    k = 0;
  else if(qbak.p[access_triangles[bottom][1]] == tbak.p[access_triangles[tip][2]])
    k = 1;
  else
    k = 2;

  t2->t[tip] = qbak.t[access_triangles[bottom][k]];
  t2->s[tip] = qbak.s[access_triangles[bottom][k]];
  t2->t[tip]->t[t2->s[tip]] = t2;
  t2->t[tip]->s[t2->s[tip]] = tip;

}




static int ErrorFlag = 0;


tetra *get_tetra(point * p, int *moves, tetra * tstart, int *flag, int *edgeface_nr)
{
  int ret, count_moves = 0;
  tetra *t, *next_tetra;

  t = tstart;

#define MAX_COUNT_MOVES 1000000

  while((ret = InTetra(t, p, edgeface_nr, &next_tetra)) == 0)
    {
      count_moves++;

      if(count_moves > MAX_COUNT_MOVES)
	{
	  ErrorFlag = 1;

	  if(count_moves > MAX_COUNT_MOVES + 10)
	    endrun(11993123);
	}

      t = next_tetra;
    }

  *moves = count_moves;
  *flag = ret;

  return t;
}




int InTetra(tetra * t, point * p, int *edgeface_nr, tetra ** nexttetra)	/* tests whether point p lies in the tetraeder */
{
  point *p0, *p1, *p2, *p3;

  p0 = t->p[0];
  p1 = t->p[1];
  p2 = t->p[2];
  p3 = t->p[3];

  if(p0 == DPinfinity || p1 == DPinfinity || p2 == DPinfinity || p3 == DPinfinity)
    {
      printf("we are in a tetraeder with in infinity point. tetra=%d\n", (int) (t - DT));
      endrun(8765);
    }

  Count_InTetra++;

  double ax = p1->xx - p0->xx;
  double ay = p1->yy - p0->yy;
  double az = p1->zz - p0->zz;

  double bx = p2->xx - p0->xx;
  double by = p2->yy - p0->yy;
  double bz = p2->zz - p0->zz;

  double cx = p3->xx - p0->xx;
  double cy = p3->yy - p0->yy;
  double cz = p3->zz - p0->zz;

  double qx = p->xx - p0->xx;
  double qy = p->yy - p0->yy;
  double qz = p->zz - p0->zz;

  double mv_data[] = { ax, bx, cx, qx, ay, by, cy, qy, az, bz, cz, qz };
  double x[3];

  int ivol, flag3, flag2, flag1, flag0;
  int count_zeros = 0;


  int status;

  status = solve_linear_equations(mv_data, x);

  if(status < 0)
    {
      ivol = Orient3d_Exact(t->p[0], t->p[1], t->p[2], t->p[3]);
      if(ivol <= 0)
	{
	  printf("flat or negatively tetrahedron found (ivol=%d)\n", ivol);
	  endrun(1121312);
	}
    }

  /* x now contains the coordinates of the point p expanded in the basis (a,b,c) */

  if(ErrorFlag)
    {
      ivol = Orient3d_Exact(p0, p1, p2, p3);
      flag3 = Orient3d_Exact(p0, p1, p2, p);
      flag2 = Orient3d_Exact(p0, p3, p1, p);
      flag1 = Orient3d_Exact(p0, p2, p3, p);
      flag0 = Orient3d_Exact(p1, p3, p2, p);

      printf("\n\nTetra=%d\n", (int) (t - DT));
      printf("ivol=%d  flag0=%d %d %d %d\n", ivol, flag0, flag1, flag2, flag3);
      printf("xx = %g %g %g   1-sum=%g\n", x[0], x[1], x[2], 1 - (x[0] + x[1] + x[2]));
      printf("a= %g %g %g\n", ax, ay, az);
      printf("b= %g %g %g\n", bx, by, bz);
      printf("c= %g %g %g\n", cx, cy, cz);
      printf("q= %g %g %g\n", qx, qy, qz);
      printf("(axb)*c) = %g\n",
	     (ay * bz - az * by) * cx + (az * bx - ax * bz) * cy + (ax * by - ay * bx) * cz);
      printf("next tetras=%d %d %d %d\n", (int) (t->t[0] - DT), (int) (t->t[1] - DT), (int) (t->t[2] - DT),
	     (int) (t->t[3] - DT));
    }

  if(status >= 0)
    {

      if(x[0] > INSIDE_EPS && x[1] > INSIDE_EPS && x[2] > INSIDE_EPS
	 && (1 - (x[0] + x[1] + x[2])) > INSIDE_EPS)
	{
	  /* looks like we are safely inside the tetrahedron */

	  return 1;		/* our point is really nicely inside the tetrahedron */
	}

      if(x[0] < -INSIDE_EPS || x[1] < -INSIDE_EPS || x[2] < -INSIDE_EPS
	 || (1 - (x[0] + x[1] + x[2])) < -INSIDE_EPS)
	{
	  /* looks like we are clearly outside the tetrahedron. 
	     Let's look for a good neighbouring tetrahedron to continue the search */

	  /* note: in the (a,b,c) basis, the center-of-mass has coordinates (1/4, 1/4, 1/4) */

	  double w, u, v;

	  if(ErrorFlag)
	    {
	      w = 0.25 / (0.25 - x[2]);
	      u = 0.25 + w * (x[0] - 0.25);
	      v = 0.25 + w * (x[1] - 0.25);
	      printf("[3] w=%g u=%g v=%g    fabs(x[2] - 0.25)=%g\n", w, u, v, fabs(x[2] - 0.25));



	      w = 0.25 / (0.25 - x[1]);
	      u = 0.25 + w * (x[0] - 0.25);
	      v = 0.25 + w * (x[2] - 0.25);
	      printf("[3] w=%g u=%g v=%g    fabs(x[1] - 0.25)=%g\n", w, u, v, fabs(x[1] - 0.25));


	      w = 0.25 / (0.25 - x[0]);
	      u = 0.25 + w * (x[1] - 0.25);
	      v = 0.25 + w * (x[2] - 0.25);
	      printf("[3] w=%g u=%g v=%g    fabs(x[0] - 0.25)=%g\n", w, u, v, fabs(x[0] - 0.25));

	    }


	  if(fabs(x[2] - 0.25) > INSIDE_EPS)
	    {
	      w = 0.25 / (0.25 - x[2]);
	      if(w > 0)
		{
		  u = 0.25 + w * (x[0] - 0.25);
		  v = 0.25 + w * (x[1] - 0.25);
		  if(u > -INSIDE_EPS && v > -INSIDE_EPS && (1 - (u + v) > -INSIDE_EPS))
		    {
		      //              printf("[3] u=%g v=%g \n", u, v);

		      *nexttetra = t->t[3];
		      return 0;
		    }
		}
	    }

	  if(fabs(x[1] - 0.25) > INSIDE_EPS)
	    {
	      w = 0.25 / (0.25 - x[1]);
	      if(w > 0)
		{
		  u = 0.25 + w * (x[0] - 0.25);
		  v = 0.25 + w * (x[2] - 0.25);
		  if(u > -INSIDE_EPS && v > -INSIDE_EPS && (1 - (u + v) > -INSIDE_EPS))
		    {
		      //   printf("[2] u=%g v=%g \n", u, v);

		      *nexttetra = t->t[2];
		      return 0;
		    }
		}
	    }


	  if(fabs(x[0] - 0.25) > INSIDE_EPS)
	    {
	      w = 0.25 / (0.25 - x[0]);
	      if(w > 0)
		{
		  u = 0.25 + w * (x[1] - 0.25);
		  v = 0.25 + w * (x[2] - 0.25);
		  if(u > -INSIDE_EPS && v > -INSIDE_EPS && (1 - (u + v) > -INSIDE_EPS))
		    {
		      // printf("[1] u=%g v=%g \n", u, v);

		      *nexttetra = t->t[1];
		      return 0;
		    }
		}
	    }

	  // printf("[1]  xx = %g %g %g\n", x[0], x[1], x[2]);

	  *nexttetra = t->t[0];
	  return 0;
	}
    }


  /* here we need to decide whether we have a degenerate case, i.e.
     whether we think the point lies on a face or an edge of the tetrahedron */

  if(ErrorFlag)
    {
      printf("doing exact test for tetra=%d\n", (int) (t - DT));
    }

  Count_InTetraExact++;

  if((ivol = Orient3d_Exact(p0, p1, p2, p3)) <= 0)
    {
      printf("flat or negatively oriented tetrahedron found (vol=%d)\n", ivol);
      endrun(13123);
    }

  flag3 = Orient3d_Exact(p0, p1, p2, p);
  flag2 = Orient3d_Exact(p0, p3, p1, p);
  flag1 = Orient3d_Exact(p0, p2, p3, p);
  flag0 = Orient3d_Exact(p1, p3, p2, p);

  if(flag0 == 0)
    count_zeros++;

  if(flag1 == 0)
    count_zeros++;

  if(flag2 == 0)
    count_zeros++;

  if(flag3 == 0)
    count_zeros++;

  if(count_zeros > 2)
    {
      printf("flags=%d %d %d %d\n", flag0, flag1, flag2, flag3);
      printf("(axb)*c) = %g\n",
	     (ay * bz - az * by) * cx + (az * bx - ax * bz) * cy + (ax * by - ay * bx) * cz);
      endrun(1312399812);
    }


  if(flag0 >= 0 && flag1 >= 0 && flag2 >= 0 && flag3 >= 0)
    {
      /* we have a point inside the tetra, but it may still be on one of the edges */

      if(count_zeros == 0)
	{
	  /* ok, let's split the tetra in 4, we are apparently well enough inside */
	  return 1;
	}

      if(count_zeros == 1)	/* we lie on a face */
	{
	  if(flag0 == 0)
	    {
	      *edgeface_nr = 0;
	      return 2;
	    }

	  if(flag1 == 0)
	    {
	      *edgeface_nr = 1;
	      return 2;
	    }

	  if(flag2 == 0)
	    {
	      *edgeface_nr = 2;
	      return 2;
	    }

	  if(flag3 == 0)
	    {
	      *edgeface_nr = 3;
	      return 2;
	    }
	}

      if(count_zeros == 2)	/* we lie on an edge */
	{
	  if(flag0 == 0 && flag1 == 0)
	    {
	      *edgeface_nr = 5;
	      return 3;
	    }

	  if(flag0 == 0 && flag2 == 0)
	    {
	      *edgeface_nr = 4;
	      return 3;
	    }

	  if(flag0 == 0 && flag3 == 0)
	    {
	      *edgeface_nr = 3;
	      return 3;
	    }

	  if(flag1 == 0 && flag2 == 0)
	    {
	      *edgeface_nr = 2;
	      return 3;
	    }

	  if(flag1 == 0 && flag3 == 0)
	    {
	      *edgeface_nr = 1;
	      return 3;
	    }

	  if(flag2 == 0 && flag3 == 0)
	    {
	      *edgeface_nr = 0;
	      return 3;
	    }
	}
    }


  /* we seem to be lying clearly outside the tetrahedron */
  /* Let's determine a suitable neighbour */

  /* if there is a single negative value, let's pick this side */

  if(flag0 < 0 && flag1 >= 0 && flag2 >= 0 && flag3 >= 0)
    {
      *nexttetra = t->t[0];
      return 0;
    }

  if(flag0 >= 0 && flag1 < 0 && flag2 >= 0 && flag3 >= 0)
    {
      *nexttetra = t->t[1];
      return 0;
    }

  if(flag0 >= 0 && flag1 >= 0 && flag2 < 0 && flag3 >= 0)
    {
      *nexttetra = t->t[2];
      return 0;
    }
  if(flag0 >= 0 && flag1 >= 0 && flag2 >= 0 && flag3 < 0)
    {
      *nexttetra = t->t[3];
      return 0;
    }

  /* there are at least two negative values. Let's pick a random one */

  int ind = -1;

  if(flag0 < 0)
    {
      if(ind < 0)
	ind = 0;
      else
	{
	  if(drand48() < 0.5)
	    ind = 0;
	}
    }

  if(flag1 < 0)
    {
      if(ind < 0)
	ind = 1;
      else
	{
	  if(drand48() < 0.5)
	    ind = 1;
	}
    }

  if(flag2 < 0)
    {
      if(ind < 0)
	ind = 2;
      else
	{
	  if(drand48() < 0.5)
	    ind = 2;
	}

    }

  if(flag3 < 0)
    {
      if(ind < 0)
	ind = 3;
      else
	{
	  if(drand48() < 0.5)
	    ind = 3;
	}
    }

  *nexttetra = t->t[ind];
  return 0;
}





void compute_circumcircles(void)
{
  int i;

  for(i = 0; i < Ndt; i++)
    {
      if(DT[i].deleted)
	continue;

      if(DT[i].p[0] == DPinfinity)
	continue;
      if(DT[i].p[1] == DPinfinity)
	continue;
      if(DT[i].p[2] == DPinfinity)
	continue;
      if(DT[i].p[3] == DPinfinity)
	continue;

      update_circumcircle(&DT[i]);
    }
}



void calc_mpz_determinant(mpz_t det, mpz_t ax, mpz_t ay, mpz_t az, mpz_t bx, mpz_t by, mpz_t bz, mpz_t cx,
			  mpz_t cy, mpz_t cz)
{
  mpz_t bz_cy, by_cz, cz_ay, cy_az, az_by, ay_bz;

  mpz_init(bz_cy);
  mpz_mul(bz_cy, bz, cy);

  mpz_init(by_cz);
  mpz_mul(by_cz, by, cz);

  mpz_init(cz_ay);
  mpz_mul(cz_ay, cz, ay);

  mpz_init(cy_az);
  mpz_mul(cy_az, cy, az);

  mpz_init(az_by);
  mpz_mul(az_by, az, by);

  mpz_init(ay_bz);
  mpz_mul(ay_bz, ay, bz);

  mpz_t bzcy_bycz, czay_cyaz, azby_aybz;

  mpz_init(bzcy_bycz);
  mpz_init(czay_cyaz);
  mpz_init(azby_aybz);

  mpz_sub(bzcy_bycz, bz_cy, by_cz);
  mpz_sub(czay_cyaz, cz_ay, cy_az);
  mpz_sub(azby_aybz, az_by, ay_bz);

  mpz_t a, b, c, ab;

  mpz_init(a);
  mpz_init(b);
  mpz_init(c);

  mpz_mul(a, bzcy_bycz, ax);
  mpz_mul(b, czay_cyaz, bx);
  mpz_mul(c, azby_aybz, cx);

  mpz_init(ab);

  mpz_add(ab, a, b);
  mpz_add(det, ab, c);

  mpz_clear(ab);
  mpz_clear(c);
  mpz_clear(b);
  mpz_clear(a);
  mpz_clear(azby_aybz);
  mpz_clear(czay_cyaz);
  mpz_clear(bzcy_bycz);
  mpz_clear(ay_bz);
  mpz_clear(az_by);
  mpz_clear(cy_az);
  mpz_clear(cz_ay);
  mpz_clear(by_cz);
  mpz_clear(bz_cy);
}


void get_circumcircle_exact(tetra * t, double *x, double *y, double *z)
{
  point *p0, *p1, *p2, *p3;

  p0 = t->p[0];
  p1 = t->p[1];
  p2 = t->p[2];
  p3 = t->p[3];

  mpz_t det, detA, detB, detC;
  mpz_t qx, qy, qz;
  mpz_t a2, b2, c2, tmp, AA, BB, CC;
  mpz_t ax, ay, az, bx, by, bz, cx, cy, cz;

  mpz_init(det);
  mpz_init(detA);
  mpz_init(detB);
  mpz_init(detC);
  mpz_init(qx);
  mpz_init(qy);
  mpz_init(qz);

  mpz_init(a2);
  mpz_init(b2);
  mpz_init(c2);
  mpz_init(tmp);
  mpz_init(AA);
  mpz_init(BB);
  mpz_init(CC);

  mpz_init(ax);
  mpz_init(ay);
  mpz_init(az);
  mpz_init(bx);
  mpz_init(by);
  mpz_init(bz);
  mpz_init(cx);
  mpz_init(cy);
  mpz_init(cz);


  mpz_set_si(tmp, p1->ix);
  mpz_sub_ui(ax, tmp, p0->ix);
  mpz_set_si(tmp, p1->iy);
  mpz_sub_ui(ay, tmp, p0->iy);
  mpz_set_si(tmp, p1->iz);
  mpz_sub_ui(az, tmp, p0->iz);

  mpz_set_si(tmp, p2->ix);
  mpz_sub_ui(bx, tmp, p0->ix);
  mpz_set_si(tmp, p2->iy);
  mpz_sub_ui(by, tmp, p0->iy);
  mpz_set_si(tmp, p2->iz);
  mpz_sub_ui(bz, tmp, p0->iz);

  mpz_set_si(tmp, p3->ix);
  mpz_sub_ui(cx, tmp, p0->ix);
  mpz_set_si(tmp, p3->iy);
  mpz_sub_ui(cy, tmp, p0->iy);
  mpz_set_si(tmp, p3->iz);
  mpz_sub_ui(cz, tmp, p0->iz);


  mpz_set(tmp, ax);
  mpz_mul(AA, tmp, ax);
  mpz_set(tmp, ay);
  mpz_mul(BB, tmp, ay);
  mpz_set(tmp, az);
  mpz_mul(CC, tmp, az);
  mpz_add(tmp, AA, BB);
  mpz_add(AA, tmp, CC);
  mpz_tdiv_q_2exp(a2, AA, 1);


  mpz_set(tmp, bx);
  mpz_mul(AA, tmp, bx);
  mpz_set(tmp, by);
  mpz_mul(BB, tmp, by);
  mpz_set(tmp, bz);
  mpz_mul(CC, tmp, bz);
  mpz_add(tmp, AA, BB);
  mpz_add(AA, tmp, CC);
  mpz_tdiv_q_2exp(b2, AA, 1);

  mpz_set(tmp, cx);
  mpz_mul(AA, tmp, cx);
  mpz_set(tmp, cy);
  mpz_mul(BB, tmp, cy);
  mpz_set(tmp, cz);
  mpz_mul(CC, tmp, cz);
  mpz_add(tmp, AA, BB);
  mpz_add(AA, tmp, CC);
  mpz_tdiv_q_2exp(c2, AA, 1);

  calc_mpz_determinant(det, ax, ay, az, bx, by, bz, cx, cy, cz);
  calc_mpz_determinant(detA, a2, ay, az, b2, by, bz, c2, cy, cz);
  calc_mpz_determinant(detB, ax, a2, az, bx, b2, bz, cx, c2, cz);
  calc_mpz_determinant(detC, ax, ay, a2, bx, by, b2, cx, cy, c2);

  mpz_cdiv_q(qx, detA, det);
  mpz_cdiv_q(qy, detB, det);
  mpz_cdiv_q(qz, detC, det);

  mpz_add_ui(AA, qx, p0->ix);
  mpz_add_ui(BB, qy, p0->iy);
  mpz_add_ui(CC, qz, p0->iz);

  double xx, yy, zz;

  xx = mpz_get_d(AA);
  yy = mpz_get_d(BB);
  zz = mpz_get_d(CC);

  xx /= (1LLu << USEDBITS);
  yy /= (1LLu << USEDBITS);
  zz /= (1LLu << USEDBITS);


  xx = xx / ConversionFac + CentralOffsetX;
  yy = yy / ConversionFac + CentralOffsetY;
  zz = zz / ConversionFac + CentralOffsetZ;


  *x = xx;
  *y = yy;
  *z = zz;


  mpz_clear(det);
  mpz_clear(detA);
  mpz_clear(detB);
  mpz_clear(detC);
  mpz_clear(qx);
  mpz_clear(qy);
  mpz_clear(qz);

  mpz_clear(a2);
  mpz_clear(b2);
  mpz_clear(c2);
  mpz_clear(tmp);
  mpz_clear(AA);
  mpz_clear(BB);
  mpz_clear(CC);

  mpz_clear(ax);
  mpz_clear(ay);
  mpz_clear(az);
  mpz_clear(bx);
  mpz_clear(by);
  mpz_clear(bz);
  mpz_clear(cx);
  mpz_clear(cy);
  mpz_clear(cz);
}


void update_circumcircle(tetra * t)
{
  point *p0, *p1, *p2, *p3;

  if(t->deleted)
    return;

  p0 = t->p[0];
  p1 = t->p[1];
  p2 = t->p[2];
  p3 = t->p[3];

  if(t->p[0] == DPinfinity)
    return;
  if(t->p[1] == DPinfinity)
    return;
  if(t->p[2] == DPinfinity)
    return;
  if(t->p[3] == DPinfinity)
    return;

  double ax = p1->xx - p0->xx;
  double ay = p1->yy - p0->yy;
  double az = p1->zz - p0->zz;

  double bx = p2->xx - p0->xx;
  double by = p2->yy - p0->yy;
  double bz = p2->zz - p0->zz;

  double cx = p3->xx - p0->xx;
  double cy = p3->yy - p0->yy;
  double cz = p3->zz - p0->zz;

  double aa = 0.5 * (ax * ax + ay * ay + az * az);
  double bb = 0.5 * (bx * bx + by * by + bz * bz);
  double cc = 0.5 * (cx * cx + cy * cy + cz * cz);

  double mv_data[] = { ax, ay, az, aa, bx, by, bz, bb, cx, cy, cz, cc };
  double x[3];

  int status = solve_linear_equations(mv_data, x);

  if(status < 0)
    {
      printf("resort to exact circum-sphere calculation\n");

      if(Orient3d_Exact(p0, p1, p2, p3) != 1)
	{
	  printf("p0 = %g %g %g\n", p0->x, p0->y, p0->z);
	  printf("p1 = %g %g %g\n", p1->x, p1->y, p1->z);
	  printf("p2 = %g %g %g\n", p2->x, p2->y, p2->z);
	  printf("p3 = %g %g %g\n", p3->x, p3->y, p3->z);

	  printf("Orient-Test=%d\n", Orient3d_Exact(p0, p1, p2, p3));
	  printf("tetra-volume=%g  tetra=%d\n", calculate_tetra_volume(p0, p1, p2, p3), (int) (t - DT));

	  return;
	}

      double xc, yc, zc;

      get_circumcircle_exact(t, &xc, &yc, &zc);

      t->c.x = xc;
      t->c.y = yc;
      t->c.z = zc;
    }
  else
    {
      x[0] += p0->xx;
      x[1] += p0->yy;
      x[2] += p0->zz;

      t->c.x = (x[0] - 1.0) / ConversionFac + CentralOffsetX;
      t->c.y = (x[1] - 1.0) / ConversionFac + CentralOffsetY;
      t->c.z = (x[2] - 1.0) / ConversionFac + CentralOffsetZ;

    }
}




int test_tetra_orientation(point * p0, point * p1, point * p2, point * p3)
{
  double nx, ny, nz;

  if(p0 == DPinfinity || p1 == DPinfinity || p2 == DPinfinity || p3 == DPinfinity)
    return +1;

  nx = (p1->yy - p0->yy) * (p2->zz - p0->zz) - (p1->zz - p0->zz) * (p2->yy - p0->yy);
  ny = (p1->zz - p0->zz) * (p2->xx - p0->xx) - (p1->xx - p0->xx) * (p2->zz - p0->zz);
  nz = (p1->xx - p0->xx) * (p2->yy - p0->yy) - (p1->yy - p0->yy) * (p2->xx - p0->xx);

  if(nx * (p3->xx - p0->xx) + ny * (p3->yy - p0->yy) + nz * (p3->zz - p0->zz) >= 0)
    return +1;
  else
    return -1;
}


double calculate_tetra_volume(point * p0, point * p1, point * p2, point * p3)
{
  double nx, ny, nz;

  if(p0 == DPinfinity || p1 == DPinfinity || p2 == DPinfinity || p3 == DPinfinity)
    return +1;

  nx = (p1->y - p0->y) * (p2->z - p0->z) - (p1->z - p0->z) * (p2->y - p0->y);
  ny = (p1->z - p0->z) * (p2->x - p0->x) - (p1->x - p0->x) * (p2->z - p0->z);
  nz = (p1->x - p0->x) * (p2->y - p0->y) - (p1->y - p0->y) * (p2->x - p0->x);

  return nx * (p3->x - p0->x) + ny * (p3->y - p0->y) + nz * (p3->z - p0->z);
}



inline void add_row(double *m, int r1, int r2, double fac)
{
  int i;

  for(i = 0; i < 4; i++)
    m[r1 * 4 + i] += fac * m[r2 * 4 + i];
}


int solve_linear_equations(double *m, double *res)
{
  int ix, iy, iz, itmp;

  if(fabs(m[4]) > fabs(m[0]))
    {
      ix = 1;
      iy = 0;
      iz = 2;
    }
  else
    {
      ix = 0;
      iy = 1;
      iz = 2;
    }

  if(fabs(m[8]) > fabs(m[ix * 4]))
    {
      ix = 2;
      iy = 0;
      iz = 1;
    }

  add_row(m, iy, ix, -m[iy * 4] / m[ix * 4]);
  add_row(m, iz, ix, -m[iz * 4] / m[ix * 4]);

  if(fabs(m[iz * 4 + 1]) > fabs(m[iy * 4 + 1]))
    {
      /* swap iy/iz */
      itmp = iy;
      iy = iz;
      iz = itmp;
    }

  if(fabs(m[iy * 4 + 1]) < 1.0e-12)
    return -1;

  add_row(m, iz, iy, -m[iz * 4 + 1] / m[iy * 4 + 1]);

  if(fabs(m[iz * 4 + 2]) < 1.0e-12)
    {
      /*
         printf("\n");
         printf("m = %15g %15g %15g    %15g\n", m[ix * 4 + 0], m[ix * 4 + 1], m[ix * 4 + 2], m[ix * 4 + 3]);
         printf("m = %15g %15g %15g    %15g\n", m[iy * 4 + 0], m[iy * 4 + 1], m[iy * 4 + 2], m[iy * 4 + 3]);
         printf("m = %15g %15g %15g    %15g\n", m[iz * 4 + 0], m[iz * 4 + 1], m[iz * 4 + 2], m[iz * 4 + 3]);
         printf("\n");
       */
      return -1;
    }
  if(fabs(m[iy * 4 + 1]) < 1.0e-12)
    {
      /*
         printf("\n");
         printf("m = %15g %15g %15g    %15g\n", m[ix * 4 + 0], m[ix * 4 + 1], m[ix * 4 + 2], m[ix * 4 + 3]);
         printf("m = %15g %15g %15g    %15g\n", m[iy * 4 + 0], m[iy * 4 + 1], m[iy * 4 + 2], m[iy * 4 + 3]);
         printf("m = %15g %15g %15g    %15g\n", m[iz * 4 + 0], m[iz * 4 + 1], m[iz * 4 + 2], m[iz * 4 + 3]);
         printf("\n");
       */
      return -2;
    }
  if(fabs(m[ix * 4]) < 1.0e-12)
    {
      /*
         printf("\n");
         printf("m = %15g %15g %15g    %15g\n", m[ix * 4 + 0], m[ix * 4 + 1], m[ix * 4 + 2], m[ix * 4 + 3]);
         printf("m = %15g %15g %15g    %15g\n", m[iy * 4 + 0], m[iy * 4 + 1], m[iy * 4 + 2], m[iy * 4 + 3]);
         printf("m = %15g %15g %15g    %15g\n", m[iz * 4 + 0], m[iz * 4 + 1], m[iz * 4 + 2], m[iz * 4 + 3]);
         printf("\n");
       */
      return -3;
    }

  res[2] = m[iz * 4 + 3] / m[iz * 4 + 2];
  res[1] = (m[iy * 4 + 3] - res[2] * m[iy * 4 + 2]) / m[iy * 4 + 1];
  res[0] = (m[ix * 4 + 3] - res[2] * m[ix * 4 + 2] - res[1] * m[ix * 4 + 1]) / m[ix * 4];

  return 0;
}



void set_integers_for_point(point * p)
{
  p->xx = (p->x - CentralOffsetX) * ConversionFac + 1.0;
  p->yy = (p->y - CentralOffsetY) * ConversionFac + 1.0;
  p->zz = (p->z - CentralOffsetZ) * ConversionFac + 1.0;

  if(p->xx < 1.0 || p->xx >= 2.0 || p->yy < 1.0 || p->yy >= 2.0 || p->zz < 1.0 || p->zz >= 2.0)
    {
      printf("%g %g %g\n", p->xx, p->yy, p->zz);
      endrun(13123);
    }

  p->ix = DOUBLE_to_VORONOIINT(p->xx);
  p->iy = DOUBLE_to_VORONOIINT(p->yy);
  p->iz = DOUBLE_to_VORONOIINT(p->zz);

  /*
     printf("%g %g %g   %d %d %d     (%g %g %g)\n", 
     p->xx, p->yy, p->zz, p->ix, p->iy, p->iz, p->x, p->y, p->z);
   */

  unsigned long long *x;

  x = (unsigned long long *) &p->xx;
  *x = *x & (~((1llu << (52 - USEDBITS)) - 1));

  x = (unsigned long long *) &p->yy;
  *x = *x & (~((1llu << (52 - USEDBITS)) - 1));

  x = (unsigned long long *) &p->zz;
  *x = *x & (~((1llu << (52 - USEDBITS)) - 1));

}




int InSphere_Exact(point * p0, point * p1, point * p2, point * p3, point * p)
{
  int ax, bx, cx, dx;
  int ay, by, cy, dy;
  int az, bz, cz, dz;

  if(p0 == DPinfinity || p1 == DPinfinity || p2 == DPinfinity || p3 == DPinfinity)
    return -1;

  ax = p0->ix - p->ix;
  ay = p0->iy - p->iy;
  az = p0->iz - p->iz;

  bx = p1->ix - p->ix;
  by = p1->iy - p->iy;
  bz = p1->iz - p->iz;

  cx = p2->ix - p->ix;
  cy = p2->iy - p->iy;
  cz = p2->iz - p->iz;

  dx = p3->ix - p->ix;
  dy = p3->iy - p->iy;
  dz = p3->iz - p->iz;


  mpz_t ab, bc, cd, da, ac, bd;

  mpz_init(ab);
  mpz_init(bc);
  mpz_init(cd);
  mpz_init(da);
  mpz_init(ac);
  mpz_init(bd);


  mpz_t tmp, AA, BB, CC;

  mpz_init(tmp);
  mpz_init(AA);
  mpz_init(BB);
  mpz_init(CC);


  mpz_set_si(tmp, ax);
  mpz_mul_si(AA, tmp, by);
  mpz_set_si(tmp, bx);
  mpz_mul_si(BB, tmp, ay);
  mpz_sub(ab, AA, BB);

  mpz_set_si(tmp, bx);
  mpz_mul_si(AA, tmp, cy);
  mpz_set_si(tmp, cx);
  mpz_mul_si(BB, tmp, by);
  mpz_sub(bc, AA, BB);

  mpz_set_si(tmp, cx);
  mpz_mul_si(AA, tmp, dy);
  mpz_set_si(tmp, dx);
  mpz_mul_si(BB, tmp, cy);
  mpz_sub(cd, AA, BB);

  mpz_set_si(tmp, dx);
  mpz_mul_si(AA, tmp, ay);
  mpz_set_si(tmp, ax);
  mpz_mul_si(BB, tmp, dy);
  mpz_sub(da, AA, BB);

  mpz_set_si(tmp, ax);
  mpz_mul_si(AA, tmp, cy);
  mpz_set_si(tmp, cx);
  mpz_mul_si(BB, tmp, ay);
  mpz_sub(ac, AA, BB);

  mpz_set_si(tmp, bx);
  mpz_mul_si(AA, tmp, dy);
  mpz_set_si(tmp, dx);
  mpz_mul_si(BB, tmp, by);
  mpz_sub(bd, AA, BB);


  mpz_t abc, bcd, cda, dab;

  mpz_init(abc);
  mpz_init(bcd);
  mpz_init(cda);
  mpz_init(dab);

  mpz_mul_si(AA, bc, az);
  mpz_mul_si(BB, ac, -bz);
  mpz_mul_si(CC, ab, cz);
  mpz_add(tmp, AA, BB);
  mpz_add(abc, tmp, CC);

  mpz_mul_si(AA, cd, bz);
  mpz_mul_si(BB, bd, -cz);
  mpz_mul_si(CC, bc, dz);
  mpz_add(tmp, AA, BB);
  mpz_add(bcd, tmp, CC);

  mpz_mul_si(AA, da, cz);
  mpz_mul_si(BB, ac, dz);
  mpz_mul_si(CC, cd, az);
  mpz_add(tmp, AA, BB);
  mpz_add(cda, tmp, CC);

  mpz_mul_si(AA, ab, dz);
  mpz_mul_si(BB, bd, az);
  mpz_mul_si(CC, da, bz);
  mpz_add(tmp, AA, BB);
  mpz_add(dab, tmp, CC);


  mpz_t a2, b2, c2, d2;

  mpz_init(a2);
  mpz_init(b2);
  mpz_init(c2);
  mpz_init(d2);

  mpz_set_si(tmp, ax);
  mpz_mul_si(AA, tmp, ax);
  mpz_set_si(tmp, ay);
  mpz_mul_si(BB, tmp, ay);
  mpz_set_si(tmp, az);
  mpz_mul_si(CC, tmp, az);
  mpz_add(tmp, AA, BB);
  mpz_add(a2, tmp, CC);

  mpz_set_si(tmp, bx);
  mpz_mul_si(AA, tmp, bx);
  mpz_set_si(tmp, by);
  mpz_mul_si(BB, tmp, by);
  mpz_set_si(tmp, bz);
  mpz_mul_si(CC, tmp, bz);
  mpz_add(tmp, AA, BB);
  mpz_add(b2, tmp, CC);

  mpz_set_si(tmp, cx);
  mpz_mul_si(AA, tmp, cx);
  mpz_set_si(tmp, cy);
  mpz_mul_si(BB, tmp, cy);
  mpz_set_si(tmp, cz);
  mpz_mul_si(CC, tmp, cz);
  mpz_add(tmp, AA, BB);
  mpz_add(c2, tmp, CC);

  mpz_set_si(tmp, dx);
  mpz_mul_si(AA, tmp, dx);
  mpz_set_si(tmp, dy);
  mpz_mul_si(BB, tmp, dy);
  mpz_set_si(tmp, dz);
  mpz_mul_si(CC, tmp, dz);
  mpz_add(tmp, AA, BB);
  mpz_add(d2, tmp, CC);

  /* now calculate final result */

  mpz_mul(AA, c2, dab);
  mpz_mul(BB, d2, abc);
  mpz_sub(tmp, AA, BB);

  mpz_mul(AA, a2, bcd);
  mpz_mul(BB, b2, cda);
  mpz_sub(CC, AA, BB);

  mpz_add(AA, tmp, CC);

  /* AA now contains the result */

  int sign = mpz_sgn(AA);

  mpz_clear(d2);
  mpz_clear(c2);
  mpz_clear(b2);
  mpz_clear(a2);
  mpz_clear(dab);
  mpz_clear(cda);
  mpz_clear(bcd);
  mpz_clear(abc);
  mpz_clear(CC);
  mpz_clear(BB);
  mpz_clear(AA);
  mpz_clear(tmp);
  mpz_clear(bd);
  mpz_clear(ac);
  mpz_clear(da);
  mpz_clear(cd);
  mpz_clear(bc);
  mpz_clear(ab);

  return sign;
}


int InSphere_Quick(point * p0, point * p1, point * p2, point * p3, point * p)
{
  double ax, bx, cx, dx;
  double ay, by, cy, dy;
  double az, bz, cz, dz;
  double a2, b2, c2, d2;
  double ab, bc, cd, da, ac, bd;
  double abc, bcd, cda, dab;
  double x;

  if(p0 == DPinfinity || p1 == DPinfinity || p2 == DPinfinity || p3 == DPinfinity)
    return -1;

  ax = p0->xx - p->xx;
  ay = p0->yy - p->yy;
  az = p0->zz - p->zz;

  bx = p1->xx - p->xx;
  by = p1->yy - p->yy;
  bz = p1->zz - p->zz;

  cx = p2->xx - p->xx;
  cy = p2->yy - p->yy;
  cz = p2->zz - p->zz;

  dx = p3->xx - p->xx;
  dy = p3->yy - p->yy;
  dz = p3->zz - p->zz;

  ab = ax * by - bx * ay;
  bc = bx * cy - cx * by;
  cd = cx * dy - dx * cy;
  da = dx * ay - ax * dy;
  ac = ax * cy - cx * ay;
  bd = bx * dy - dx * by;

  abc = az * bc - bz * ac + cz * ab;
  bcd = bz * cd - cz * bd + dz * bc;
  cda = cz * da + dz * ac + az * cd;
  dab = dz * ab + az * bd + bz * da;

  a2 = ax * ax + ay * ay + az * az;
  b2 = bx * bx + by * by + bz * bz;
  c2 = cx * cx + cy * cy + cz * cz;
  d2 = dx * dx + dy * dy + dz * dz;

  x = ((c2 * dab - d2 * abc) + (a2 * bcd - b2 * cda));

  if(x < 0)
    return -1;
  if(x > 0)
    return +1;

  return 0;
}




int InSphere_Errorbound(point * p0, point * p1, point * p2, point * p3, point * p)
{
  double ax, bx, cx, dx;
  double ay, by, cy, dy;
  double az, bz, cz, dz;
  double a2, b2, c2, d2;
  double ab, bc, cd, da, ac, bd;
  double abc, bcd, cda, dab;
  double x;

  if(p0 == DPinfinity || p1 == DPinfinity || p2 == DPinfinity || p3 == DPinfinity)
    return -1;

  ax = p0->xx - p->xx;
  ay = p0->yy - p->yy;
  az = p0->zz - p->zz;

  bx = p1->xx - p->xx;
  by = p1->yy - p->yy;
  bz = p1->zz - p->zz;

  cx = p2->xx - p->xx;
  cy = p2->yy - p->yy;
  cz = p2->zz - p->zz;

  dx = p3->xx - p->xx;
  dy = p3->yy - p->yy;
  dz = p3->zz - p->zz;

  double axby = ax * by;
  double bxay = bx * ay;
  double bxcy = bx * cy;
  double cxby = cx * by;
  double cxdy = cx * dy;
  double dxcy = dx * cy;
  double dxay = dx * ay;
  double axdy = ax * dy;
  double axcy = ax * cy;
  double cxay = cx * ay;
  double bxdy = bx * dy;
  double dxby = dx * by;

  ab = axby - bxay;
  bc = bxcy - cxby;
  cd = cxdy - dxcy;
  da = dxay - axdy;
  ac = axcy - cxay;
  bd = bxdy - dxby;

  abc = az * bc - bz * ac + cz * ab;
  bcd = bz * cd - cz * bd + dz * bc;
  cda = cz * da + dz * ac + az * cd;
  dab = dz * ab + az * bd + bz * da;

  a2 = ax * ax + ay * ay + az * az;
  b2 = bx * bx + by * by + bz * bz;
  c2 = cx * cx + cy * cy + cz * cz;
  d2 = dx * dx + dy * dy + dz * dz;

  x = ((c2 * dab - d2 * abc) + (a2 * bcd - b2 * cda));


  /* calculate absolute maximum size */

  ab = fabs(axby) + fabs(bxay);
  bc = fabs(bxcy) + fabs(cxby);
  cd = fabs(cxdy) + fabs(dxcy);
  da = fabs(dxay) + fabs(axdy);
  ac = fabs(axcy) + fabs(cxay);
  bd = fabs(bxdy) + fabs(dxby);

  az = fabs(az);
  bz = fabs(bz);
  cz = fabs(cz);
  dz = fabs(dz);

  abc = az * bc + bz * ac + cz * ab;
  bcd = bz * cd + cz * bd + dz * bc;
  cda = cz * da + dz * ac + az * cd;
  dab = dz * ab + az * bd + bz * da;

  double sizelimit = ((c2 * dab + d2 * abc) + (a2 * bcd + b2 * cda));

  double errbound = 1.0e-14 * sizelimit;

  if(x < -errbound)
    return -1;
  else if(x > errbound)
    return +1;

  return 0;
}




int InSphere_Gauss(point * p0, point * p1, point * p2, point * p3, point * p)
{
  /* according to my benchmarks, this test is a bit slower than InSphere_Quick() */

  if(p0 == DPinfinity || p1 == DPinfinity || p2 == DPinfinity || p3 == DPinfinity)
    return -1;

  double ax = p1->xx - p0->xx;
  double ay = p1->yy - p0->yy;
  double az = p1->zz - p0->zz;

  double bx = p2->xx - p0->xx;
  double by = p2->yy - p0->yy;
  double bz = p2->zz - p0->zz;

  double cx = p3->xx - p0->xx;
  double cy = p3->yy - p0->yy;
  double cz = p3->zz - p0->zz;

  double aa = 0.5 * (ax * ax + ay * ay + az * az);
  double bb = 0.5 * (bx * bx + by * by + bz * bz);
  double cc = 0.5 * (cx * cx + cy * cy + cz * cz);

  double mv_data[] = { ax, ay, az, aa, bx, by, bz, bb, cx, cy, cz, cc };
  double x[3];

  int status = solve_linear_equations(mv_data, x);

  if(status < 0)
    {
      return 0;
      printf("trouble in in-circle test\n");
      endrun(123121388);
    }

  double r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];

  double dx = p->xx - (p0->xx + x[0]);
  double dy = p->yy - (p0->yy + x[1]);
  double dz = p->zz - (p0->zz + x[2]);

  double d2 = dx * dx + dy * dy + dz * dz;

  if((r2 - d2) > 1.0e-5 * r2)
    return 1;
  else if((r2 - d2) < -1.0e-5 * r2)
    return -1;

  return 0;
}



int Orient3d_Exact(point * p0, point * p1, point * p2, point * p3)
{
  int ax, bx, cx;
  int ay, by, cy;
  int az, bz, cz;

  ax = p0->ix - p3->ix;
  ay = p0->iy - p3->iy;
  az = p0->iz - p3->iz;

  bx = p1->ix - p3->ix;
  by = p1->iy - p3->iy;
  bz = p1->iz - p3->iz;

  cx = p2->ix - p3->ix;
  cy = p2->iy - p3->iy;
  cz = p2->iz - p3->iz;

  mpz_t bz_cy, by_cz, cz_ay, cy_az, az_by, ay_bz;
  mpz_t bz2, by2, cz2, cy2, az2, ay2;

  mpz_init(bz_cy);
  mpz_init_set_si(bz2, bz);
  mpz_mul_si(bz_cy, bz2, cy);

  mpz_init(by_cz);
  mpz_init_set_si(by2, by);
  mpz_mul_si(by_cz, by2, cz);

  mpz_init(cz_ay);
  mpz_init_set_si(cz2, cz);
  mpz_mul_si(cz_ay, cz2, ay);

  mpz_init(cy_az);
  mpz_init_set_si(cy2, cy);
  mpz_mul_si(cy_az, cy2, az);

  mpz_init(az_by);
  mpz_init_set_si(az2, az);
  mpz_mul_si(az_by, az2, by);

  mpz_init(ay_bz);
  mpz_init_set_si(ay2, ay);
  mpz_mul_si(ay_bz, ay2, bz);

  mpz_t bzcy_bycz, czay_cyaz, azby_aybz;

  mpz_init(bzcy_bycz);
  mpz_init(czay_cyaz);
  mpz_init(azby_aybz);

  mpz_sub(bzcy_bycz, bz_cy, by_cz);
  mpz_sub(czay_cyaz, cz_ay, cy_az);
  mpz_sub(azby_aybz, az_by, ay_bz);

  mpz_t a, b, c, ab, res;

  mpz_init(a);
  mpz_init(b);
  mpz_init(c);

  mpz_mul_si(a, bzcy_bycz, ax);
  mpz_mul_si(b, czay_cyaz, bx);
  mpz_mul_si(c, azby_aybz, cx);

  mpz_init(ab);
  mpz_init(res);

  mpz_add(ab, a, b);
  mpz_add(res, ab, c);

  int sign = mpz_sgn(res);

  mpz_clear(res);
  mpz_clear(ab);
  mpz_clear(c);
  mpz_clear(b);
  mpz_clear(a);
  mpz_clear(azby_aybz);
  mpz_clear(czay_cyaz);
  mpz_clear(bzcy_bycz);
  mpz_clear(ay2);
  mpz_clear(ay_bz);
  mpz_clear(az2);
  mpz_clear(az_by);
  mpz_clear(cy2);
  mpz_clear(cy_az);
  mpz_clear(cz2);
  mpz_clear(cz_ay);
  mpz_clear(by2);
  mpz_clear(by_cz);
  mpz_clear(bz2);
  mpz_clear(bz_cy);

  return sign;
}



int Orient3d_Quick(point * p0, point * p1, point * p2, point * p3)
{
  double ax, bx, cx;
  double ay, by, cy;
  double az, bz, cz;

  ax = p0->xx - p3->xx;
  ay = p0->yy - p3->yy;
  az = p0->zz - p3->zz;

  bx = p1->xx - p3->xx;
  by = p1->yy - p3->yy;
  bz = p1->zz - p3->zz;

  cx = p2->xx - p3->xx;
  cy = p2->yy - p3->yy;
  cz = p2->zz - p3->zz;

  double x = (ax * (bz * cy - by * cz) + bx * (cz * ay - cy * az) + cx * (az * by - ay * bz));

  if(x < 0)
    return -1;
  else if(x > 0)
    return +1;

  return 0;
}

#endif

#endif
