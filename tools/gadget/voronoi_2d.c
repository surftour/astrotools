#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef VORONOI

#include <gmp.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"


#ifdef TWODIMS			/* will only be compiled in 2D case */

#define INSIDE_EPS   1.0e-8
#define USEDBITS 31
#define DOUBLE_to_VORONOIINT(y)       ((int)(((*((long long *) &y)) & 0xFFFFFFFFFFFFFllu) >> (52 - USEDBITS)))



void initialize_and_create_first_tetra(void)
{
  point *p;
  int i, n;

  for(i = 0; i < N_gas; i++)
    {
      SphP[i].Hsml = SphP[i].MaxDelaunayRadius;
      SphP[i].MaxDelaunayRadius = MAX_REAL_NUMBER;
    }


  MaxNdp = 10 * N_gas;
  MaxNdt = MaxNdp * 16;
  MaxNvf = MaxNdt;

  Ndp = 0;
  Nvf = 0;
  Ndt = 0;


  VF = mymalloc(MaxNvf * sizeof(face));

  DP = mymalloc((MaxNdp + 5) * sizeof(point));
  DP += 5;

  DT = mymalloc(MaxNdt * sizeof(tetra));


  /* construct all encompassing huge triangle */

  double box, tetra_incircle, tetra_sidelength, tetra_height;

  box = boxSize_X;
  if(box < boxSize_Y)
    box = boxSize_Y;

  tetra_incircle = 2.001 * (1 + sqrt(3)) / 3.0 * box;	/* to give room for ghost particles needed for periodic/reflective
							   boundary conditions, the incircle is twice as large, i.e.
							   [-0.5*box, 1.5*box,-0.5*box, 1.5*box] should be inside triangle */
  tetra_sidelength = tetra_incircle * sqrt(12);
  tetra_height = sqrt(3.0) / 2 * tetra_sidelength;

  printf("side-length of enclosing triangle=%g tetra_height=%g box=%g\n", tetra_sidelength, tetra_height,
	 box);

  /* first, let's make the points */
  DP[-3].x = 0.5 * tetra_sidelength;
  DP[-3].y = -1.0 / 3 * tetra_height;
  DP[-3].z = 0;

  DP[-2].x = 0;
  DP[-2].y = 2.0 / 3 * tetra_height;
  DP[-2].z = 0;

  DP[-1].x = -0.5 * tetra_sidelength;
  DP[-1].y = -1.0 / 3 * tetra_height;
  DP[-1].z = 0;

  for(i = -3; i <= -1; i++)
    {
      DP[i].x += 0.5 * box;
      DP[i].y += 1.0 / 3 * tetra_height - 0.5 * box;
    }

  for(i = -3, p = &DP[-3]; i < 0; i++, p++)
    {
      p->index = -1;
      p->task = ThisTask;
    }

  /* we also define a neutral element at infinity */
  DPinfinity = &DP[-4];

  DPinfinity->x = 0;
  DPinfinity->y = 0;
  DPinfinity->z = 0;
  DPinfinity->index = -1;
  DPinfinity->task = ThisTask;

  /* now let's make the big triangle */
  DT[0].p[0] = DP - 3;
  DT[0].p[1] = DP - 2;
  DT[0].p[2] = DP - 1;


  /* On the outer faces, we attach tetrahedra with the neutral element as tip.
   * This way we will be able to navigate nicely within the tesselation, 
   * and all tetrahedra have defined neighbouring tetrahedra.
   */

  for(i = 0; i < 3; i++)
    {
      n = i + 1;		/* tetra index */

      DT[0].t[i] = &DT[n];
      DT[0].s[i] = 2;

      DT[n].t[2] = &DT[0];
      DT[n].s[2] = i;
      DT[n].p[2] = DPinfinity;
    }


  DT[1].p[0] = DT[0].p[2];
  DT[1].p[1] = DT[0].p[1];

  DT[2].p[0] = DT[0].p[0];
  DT[2].p[1] = DT[0].p[2];

  DT[3].p[0] = DT[0].p[1];
  DT[3].p[1] = DT[0].p[0];

  DT[1].t[0] = &DT[3];
  DT[3].t[1] = &DT[1];
  DT[1].s[0] = 1;
  DT[3].s[1] = 0;

  DT[1].t[1] = &DT[2];
  DT[2].t[0] = &DT[1];
  DT[1].s[1] = 0;
  DT[2].s[0] = 1;

  DT[2].t[1] = &DT[3];
  DT[3].t[0] = &DT[2];
  DT[2].s[1] = 0;
  DT[3].s[0] = 1;


  Ndt = 4;			/* we'll start out with 4 triangles */

  for(i = 0; i < Ndt; i++)
    DT[i].deleted = 0;

  CentralOffsetX = 0.5 * box - 0.5000001 * tetra_sidelength;
  CentralOffsetY = -0.5000001 * box;

  ConversionFac = 1.0 / (1.001 * tetra_sidelength);

  for(i = -3; i < 0; i++)
    set_integers_for_point(&DP[i]);
}




tetra *insert_point(point * p, tetra * tstart)	/* returns a triangle that (currently) contains the point p */
{
  tetra *t0, *t1, *t2, *t3;
  tetra *tetra_with_p;
  int moves, degenerate_flag;


  /* first, need to do a point location */
  t0 = get_triangle(p, &moves, &degenerate_flag, tstart);

  tetra_with_p = t0;

  if(degenerate_flag == 1)	/* that's the normal split of a triangle into 3 */
    {
      /* we now need to split this triangle into three  */
      t1 = &DT[Ndt++];
      t2 = &DT[Ndt++];
      *t1 = *t0;
      *t2 = *t0;

      if(Ndt > MaxNdt)
	endrun(121997);

      make_a_1_to_3_flip(p, t0, t1, t2);

      check_edge_and_flip_if_needed(p, t0);
      check_edge_and_flip_if_needed(p, t1);
      check_edge_and_flip_if_needed(p, t2);
    }
  else
    {
      degenerate_flag -= 10;

      t1 = t0->t[degenerate_flag];

      /* we now need to split this into two triangles */
      t2 = &DT[Ndt++];
      t3 = &DT[Ndt++];
      *t2 = *t0;
      *t3 = *t1;

      if(Ndt > MaxNdt)
	endrun(121999);

      make_a_2_to_4_flip(p, t0, t1, t2, t3, degenerate_flag, t0->s[degenerate_flag]);

      check_edge_and_flip_if_needed(p, t0);
      check_edge_and_flip_if_needed(p, t1);
      check_edge_and_flip_if_needed(p, t2);
      check_edge_and_flip_if_needed(p, t3);
    }

  return tetra_with_p;
}



void make_a_2_to_4_flip(point * p, tetra * t0, tetra * t1, tetra * t2, tetra * t3, int i0, int j0)
{
  int i1, i2, j1, j2;

  CountFlips++;
  Count_2_to_4_Flips2d++;


  i1 = i0 + 1;
  i2 = i0 + 2;
  j1 = j0 + 1;
  j2 = j0 + 2;

  if(i1 > 2)
    i1 -= 3;
  if(i2 > 2)
    i2 -= 3;

  if(j1 > 2)
    j1 -= 3;
  if(j2 > 2)
    j2 -= 3;

  t0->p[i1] = p;
  t1->p[j2] = p;
  t2->p[i2] = p;
  t3->p[j1] = p;

  t0->t[i0] = t1;
  t1->t[j0] = t0;
  t0->s[i0] = j0;
  t1->s[j0] = i0;

  t1->t[j1] = t3;
  t3->t[j2] = t1;
  t1->s[j1] = j2;
  t3->s[j2] = j1;

  t2->t[i1] = t0;
  t0->t[i2] = t2;
  t2->s[i1] = i2;
  t0->s[i2] = i1;

  t2->t[i0] = t3;
  t3->t[j0] = t2;
  t2->s[i0] = j0;
  t3->s[j0] = i0;

  t0->t[i1]->t[t0->s[i1]] = t0;
  t1->t[j2]->t[t1->s[j2]] = t1;
  t2->t[i2]->t[t2->s[i2]] = t2;
  t3->t[j1]->t[t3->s[j1]] = t3;
}


void make_a_1_to_3_flip(point * p, tetra * t0, tetra * t1, tetra * t2)
{
  CountFlips++;
  Count_1_to_3_Flips2d++;

  t0->p[0] = p;
  t1->p[1] = p;
  t2->p[2] = p;


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

  t0->t[0]->t[t0->s[0]] = t0;
  t1->t[1]->t[t1->s[1]] = t1;
  t2->t[2]->t[t2->s[2]] = t2;
}


void check_edge_and_flip_if_needed(point * p, tetra * t)
{
  tetra *tt, *t0, *t2;
  point *pp;
  int pi, pi1, pi2;
  int ni, ni1, ni2;
  int st2, st0;

  if(t->p[0] == p)
    pi = 0;
  else if(t->p[1] == p)
    pi = 1;
  else
    pi = 2;

  /* get the point that lies accross the edge to obtain the quadriliteral */

  tt = t->t[pi];
  ni = t->s[pi];

  pp = tt->p[ni];


  int ret, ret_exact;

  ret = InCircle_Errorbound(t->p[0], t->p[1], t->p[2], pp);
  CountInSphereTests++;

  if(ret != 0)
    ret_exact = ret;
  else
    {
      ret_exact = InCircle_Exact(t->p[0], t->p[1], t->p[2], pp);
      CountInSphereTestsExact++;
    }


  if(ret_exact > 0)
    {
      /* pp lies in the triangle, the edge is not Delaunay. Need to do a flip */

      CountFlips++;

      ni1 = ni + 1;
      if(ni1 > 2)
	ni1 -= 3;
      ni2 = ni + 2;
      if(ni2 > 2)
	ni2 -= 3;

      pi1 = pi + 1;
      if(pi1 > 2)
	pi1 -= 3;
      pi2 = pi + 2;
      if(pi2 > 2)
	pi2 -= 3;


      t0 = tt->t[ni1];
      t2 = t->t[pi1];

      st0 = tt->s[ni1];
      st2 = t->s[pi1];

      /* change the points of the triangles */
      t->p[pi2] = pp;
      tt->p[ni2] = p;

      /* change the pointers to the neighbouring triangles, and fix
         the adjency relations */

      t->t[pi1] = tt;
      tt->t[ni1] = t;
      t->s[pi1] = ni1;
      tt->s[ni1] = pi1;


      t->t[pi] = t0;
      t0->t[st0] = t;
      t->s[pi] = st0;
      t0->s[st0] = pi;


      tt->t[ni] = t2;
      t2->t[st2] = tt;
      tt->s[ni] = st2;
      t2->s[st2] = ni;

      /* now we need to test also the two sides opposite of p */

      check_edge_and_flip_if_needed(p, t);
      check_edge_and_flip_if_needed(p, tt);
    }
}




tetra *get_triangle(point * p, int *moves, int *degenerate_flag, tetra * tstart)
{
  int count_moves = 0;
  int ret;
  tetra *t, *next_tetra;

  t = tstart;

#define MAX_COUNT_MOVES 1000000

  while((ret = FindTriangle(t, p, degenerate_flag, &next_tetra)) == 0)
    {
      /* we need to see in which of the three possible neighbouring triangles
         we should walk. We'll choose the one which lies along the face that
         is traversed by a line from the cm of the triangle to the point in
         question.
       */
      count_moves++;

      if(count_moves > MAX_COUNT_MOVES)
	{
	  printf("ta=%d triangle=%d\n", ThisTask, (int) (t - DT));

	  if(count_moves > MAX_COUNT_MOVES + 10)
	    endrun(113123);
	}

      t = next_tetra;
    }

  *moves = count_moves;

  return t;
}




inline void add_row_2d(double *m, int r1, int r2, double fac)
{
  int i;

  for(i = 0; i < 3; i++)
    m[r1 * 3 + i] += fac * m[r2 * 3 + i];
}


int solve_linear_equations_2d(double *m, double *res)
{
  int ix, iy;

  if(fabs(m[0]) > fabs(m[3]))
    {
      ix = 0;
      iy = 1;
    }
  else
    {
      ix = 1;
      iy = 0;
    }

  if(fabs(m[ix * 3]) < 1.0e-12)
    {
      return -1;
    }

  add_row_2d(m, iy, ix, -m[iy * 3] / m[ix * 3]);


  res[1] = m[iy * 3 + 2] / m[iy * 3 + 1];
  res[0] = (m[ix * 3 + 2] - res[1] * m[ix * 3 + 1]) / m[ix * 3];

  return 0;
}



/* tests whether point p lies in the triangle, on an edge, or outside. In the latter case, a nighbouring triangle is returned */
int FindTriangle(tetra * t, point * p, int *degnerate_flag, tetra ** nexttetra)
{
  point *p0, *p1, *p2;

  p0 = t->p[0];
  p1 = t->p[1];
  p2 = t->p[2];

  if(p0 == DPinfinity || p1 == DPinfinity || p2 == DPinfinity)
    {
      printf("we are in a triangle with in infinity point. tetra=%d\n", (int) (t - DT));
      endrun(87658);
    }

  Count_InTetra++;

  double ax = p1->xx - p0->xx;
  double ay = p1->yy - p0->yy;

  double bx = p2->xx - p0->xx;
  double by = p2->yy - p0->yy;

  double qx = p->xx - p0->xx;
  double qy = p->yy - p0->yy;

  double mv_data[] = { ax, bx, qx, ay, by, qy };
  double x[2];

  int ivol, flag2, flag1, flag0;
  int count_zeros = 0;


  int status;

  status = solve_linear_equations_2d(mv_data, x);


  if(status < 0)
    {
      ivol = Orient2d_Exact(t->p[0], t->p[1], t->p[2]);
      if(ivol <= 0)
	{
	  printf("flat or negatively triangle found (ivol=%d)\n", ivol);
	  endrun(11213192);
	}
    }

  if(status >= 0)
    {
      if(x[0] > INSIDE_EPS && x[1] > INSIDE_EPS && (1 - (x[0] + x[1])) > INSIDE_EPS)
	{
	  /* looks like we are safely inside the triangle */

	  *degnerate_flag = 1;
	  return 1;
	}


      if(x[0] < -INSIDE_EPS || x[1] < -INSIDE_EPS || (1 - (x[0] + x[1])) < -INSIDE_EPS)
	{
	  /* looks like we are clearly outside the triangle. 
	     Let's look for a good neighbouring triangle to continue the search */

	  /* note: in the (a,b) basis, the center-of-mass has coordinates (1/3, 1/3) */

	  double w, u;

	  if(fabs(x[1] - (1.0 / 3)) > INSIDE_EPS)
	    {
	      w = (1.0 / 3) / ((1.0 / 3) - x[1]);
	      if(w > 0)
		{
		  u = (1.0 / 3) + w * (x[0] - (1.0 / 3));
		  if(u > -INSIDE_EPS && (1 - u) > -INSIDE_EPS)
		    {
		      *nexttetra = t->t[2];
		      return 0;
		    }
		}
	    }


	  if(fabs(x[0] - (1.0 / 3)) > INSIDE_EPS)
	    {
	      w = (1.0 / 3) / ((1.0 / 3) - x[0]);
	      if(w > 0)
		{
		  u = (1.0 / 3) + w * (x[1] - (1.0 / 3));
		  if(u > -INSIDE_EPS && (1 - u) > -INSIDE_EPS)
		    {
		      *nexttetra = t->t[1];
		      return 0;
		    }
		}
	    }

	  *nexttetra = t->t[0];
	  return 0;
	}
    }

  /* here we need to decide whether we have a degenerate case, i.e.
     whether we think the point lies on an edge of the triangle */

  Count_InTetraExact++;

  ivol = Orient2d_Exact(t->p[0], t->p[1], t->p[2]);

  if(ivol <= 0)
    {
      printf("flat or negatively triangle found (ivol=%d)\n", ivol);
      endrun(1128813192);
    }

  flag0 = Orient2d_Exact(p1, p2, p);
  flag1 = Orient2d_Exact(p2, p0, p);
  flag2 = Orient2d_Exact(p0, p1, p);

  if(flag0 == 0)
    count_zeros++;

  if(flag1 == 0)
    count_zeros++;

  if(flag2 == 0)
    count_zeros++;

  if(count_zeros >= 2)
    {
      printf("flags=%d %d %d\n", flag0, flag1, flag2);
      endrun(1312399812);
    }

  if(flag0 >= 0 && flag1 >= 0 && flag2 >= 0)
    {

      /* we have a point inside the triangle, but it may still be on one of the edges */

      if(count_zeros == 0)
	{
	  /* ok, we are inside */
	  *degnerate_flag = 1;
	  return 1;
	}

      if(count_zeros == 1)	/* we lie on a face */
	{
	  if(flag2 == 0)
	    {
	      *degnerate_flag = 12;
	      return 12;	/* point lies on side A */
	    }
	  if(flag1 == 0)
	    {
	      *degnerate_flag = 11;
	      return 11;	/* point lies on side C */
	    }

	  if(flag0 == 0)
	    {
	      *degnerate_flag = 10;
	      return 10;	/* point lies on side B */
	    }
	}
    }

  /* we are clearly outside, let's select the suitable neighbour */

  if(flag0 < 0 && flag1 >= 0 && flag2 >= 0)
    {
      *nexttetra = t->t[0];
      return 0;
    }

  if(flag0 >= 0 && flag1 < 0 && flag2 >= 0)
    {
      *nexttetra = t->t[1];
      return 0;
    }

  if(flag0 >= 0 && flag1 >= 0 && flag2 < 0)
    {
      *nexttetra = t->t[2];
      return 0;
    }

  /* there are apparently two negative values. Let's pick a random one */

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
	ind = 0;
      else
	{
	  if(drand48() < 0.5)
	    ind = 0;
	}
    }

  if(flag2 < 0)
    {
      if(ind < 0)
	ind = 0;
      else
	{
	  if(drand48() < 0.5)
	    ind = 0;
	}
    }

  *nexttetra = t->t[ind];
  return 0;
}


/* tests whether point p lies in the circumcircle around triangle p0,p1,p3 */

int InCircle_Quick(point * p0, point * p1, point * p2, point * p)
{
  double ax, ay, bx, by, cx, cy;
  double ab, bc, ca, a2, b2, c2, x;

  if(p0 == DPinfinity || p1 == DPinfinity || p2 == DPinfinity || p == DPinfinity)
    return -1;

  ax = p0->xx - p->xx;
  ay = p0->yy - p->yy;
  bx = p1->xx - p->xx;
  by = p1->yy - p->yy;
  cx = p2->xx - p->xx;
  cy = p2->yy - p->yy;

  ab = ax * by - bx * ay;
  bc = bx * cy - cx * by;
  ca = cx * ay - ax * cy;

  a2 = ax * ax + ay * ay;
  b2 = bx * bx + by * by;
  c2 = cx * cx + cy * cy;

  x = a2 * bc + b2 * ca + c2 * ab;

  if(x < 0)
    return -1;
  if(x > 0)
    return +1;

  return 0;
}


int InCircle_Errorbound(point * p0, point * p1, point * p2, point * p)
{
  double ax, ay, bx, by, cx, cy;
  double ab, bc, ca, a2, b2, c2, x;
  double axby, bxay, bxcy, cxby, cxay, axcy;

  if(p0 == DPinfinity || p1 == DPinfinity || p2 == DPinfinity || p == DPinfinity)
    return -1;

  ax = p0->xx - p->xx;
  ay = p0->yy - p->yy;
  bx = p1->xx - p->xx;
  by = p1->yy - p->yy;
  cx = p2->xx - p->xx;
  cy = p2->yy - p->yy;

  axby = ax * by;
  bxay = bx * ay;
  bxcy = bx * cy;
  cxby = cx * by;
  cxay = cx * ay;
  axcy = ax * cy;

  ca = cxay - axcy;
  ab = axby - bxay;
  bc = bxcy - cxby;

  a2 = ax * ax + ay * ay;
  b2 = bx * bx + by * by;
  c2 = cx * cx + cy * cy;

  x = a2 * bc + b2 * ca + c2 * ab;

  /* calculate absolute maximum size */

  double sizelimit =
    a2 * (fabs(bxcy) + fabs(cxby)) + b2 * (fabs(cxay) + fabs(axcy)) + c2 * (fabs(axby) + fabs(bxay));

  double errbound = 1.0e-14 * sizelimit;

  if(x < -errbound)
    return -1;
  else if(x > errbound)
    return +1;

  return 0;
}


int InCircle_Exact(point * p0, point * p1, point * p2, point * p)
{
  int ax, ay, bx, by, cx, cy;

  if(p0 == DPinfinity || p1 == DPinfinity || p2 == DPinfinity || p == DPinfinity)
    return -1;

  ax = p0->ix - p->ix;
  ay = p0->iy - p->iy;
  bx = p1->ix - p->ix;
  by = p1->iy - p->iy;
  cx = p2->ix - p->ix;
  cy = p2->iy - p->iy;

  mpz_t axby, bxay, bxcy, cxby, cxay, axcy, tmp;

  mpz_init(tmp);

  mpz_init(axby);
  mpz_set_si(tmp, ax);
  mpz_mul_si(axby, tmp, by);
  mpz_init(bxay);
  mpz_set_si(tmp, bx);
  mpz_mul_si(bxay, tmp, ay);
  mpz_init(bxcy);
  mpz_set_si(tmp, bx);
  mpz_mul_si(bxcy, tmp, cy);
  mpz_init(cxby);
  mpz_set_si(tmp, cx);
  mpz_mul_si(cxby, tmp, by);
  mpz_init(cxay);
  mpz_set_si(tmp, cx);
  mpz_mul_si(cxay, tmp, ay);
  mpz_init(axcy);
  mpz_set_si(tmp, ax);
  mpz_mul_si(axcy, tmp, cy);

  mpz_t ca, ab, bc;

  mpz_init(ca);
  mpz_init(ab);
  mpz_init(bc);

  mpz_sub(ca, cxay, axcy);
  mpz_sub(ab, axby, bxay);
  mpz_sub(bc, bxcy, cxby);


  mpz_t AA, BB, a2, b2, c2;

  mpz_init(AA);
  mpz_init(BB);
  mpz_init(a2);
  mpz_init(b2);
  mpz_init(c2);

  mpz_set_si(tmp, ax);
  mpz_mul_si(AA, tmp, ax);
  mpz_set_si(tmp, ay);
  mpz_mul_si(BB, tmp, ay);
  mpz_add(a2, AA, BB);

  mpz_set_si(tmp, bx);
  mpz_mul_si(AA, tmp, bx);
  mpz_set_si(tmp, by);
  mpz_mul_si(BB, tmp, by);
  mpz_add(b2, AA, BB);

  mpz_set_si(tmp, cx);
  mpz_mul_si(AA, tmp, cx);
  mpz_set_si(tmp, cy);
  mpz_mul_si(BB, tmp, cy);
  mpz_add(c2, AA, BB);

  /* now calculate the final result */

  mpz_mul(AA, a2, bc);
  mpz_mul(BB, b2, ca);
  mpz_add(tmp, AA, BB);
  mpz_mul(BB, c2, ab);
  mpz_add(AA, BB, tmp);

  int sign = mpz_sgn(AA);

  mpz_clear(c2);
  mpz_clear(b2);
  mpz_clear(a2);
  mpz_clear(BB);
  mpz_clear(AA);
  mpz_clear(bc);
  mpz_clear(ab);
  mpz_clear(ca);
  mpz_clear(axcy);
  mpz_clear(cxay);
  mpz_clear(cxby);
  mpz_clear(bxcy);
  mpz_clear(bxay);
  mpz_clear(axby);
  mpz_clear(tmp);

  return sign;
}




double test_triangle_orientation(point * p0, point * p1, point * p2)
{
  return (p1->x - p0->x) * (p2->y - p0->y) - (p1->y - p0->y) * (p2->x - p0->x);
}


int Orient2d_Quick(point * p0, point * p1, point * p2)
{
  double x;

  x = (p1->xx - p0->xx) * (p2->yy - p0->yy) - (p1->yy - p0->yy) * (p2->xx - p0->xx);

  if(x < 0)
    return -1;
  if(x > 0)
    return +1;
  return 0;
}

int Orient2d_Exact(point * p0, point * p1, point * p2)
{
  signed long long dx1, dy1, dx2, dy2, x;

  dx1 = (p1->ix - p0->ix);
  dy1 = (p1->iy - p0->iy);
  dx2 = (p2->ix - p0->ix);
  dy2 = (p2->iy - p0->iy);

  x = dx1 * dy2 - dy1 * dx2;

  if(x < 0)
    return -1;
  if(x > 0)
    return +1;
  return 0;
}





const int edge_start[3] = { 1, 2, 0 };
const int edge_end[3] = { 2, 0, 1 };


void process_edge_faces_and_volumes(tetra * t, int nr)
{
  int i, j;
  face *f;
  tetra *q;
  double nx, ny;
  double sx, sy;
  double hx, hy;
  double dvol, h;

  i = edge_start[nr];
  j = edge_end[nr];

  q = t->t[nr];

  t->egde_visited |= (1 << nr);
  q->egde_visited |= (1 << (t->s[nr]));


  if(Nvf + 1 >= MaxNvf)
    endrun(1123);

  f = &VF[Nvf++];

  f->p1 = t->p[i];
  f->p2 = t->p[j];

  f->cx = 0.5 * (t->c.x + q->c.x);
  f->cy = 0.5 * (t->c.y + q->c.y);
  f->cz = 0;


  nx = t->c.x - q->c.x;
  ny = t->c.y - q->c.y;

  f->area = sqrt(nx * nx + ny * ny);


  hx = 0.5 * (t->p[i]->x - t->p[j]->x);
  hy = 0.5 * (t->p[i]->y - t->p[j]->y);

  h = sqrt(hx * hx + hy * hy);
  dvol = 0.5 * f->area * h;

  if(t->p[i]->task == ThisTask && t->p[i]->index >= 0 && t->p[i]->index < N_gas)
    {
      SphP[t->p[i]->index].Volume += dvol;

      /* let's now compute the center-of-mass of the pyramid at the bottom top */
      sx = (2.0 / 3) * f->cx + (1.0 / 3) * t->p[i]->x;
      sy = (2.0 / 3) * f->cy + (1.0 / 3) * t->p[i]->y;

      SphP[t->p[i]->index].Center[0] += dvol * sx;
      SphP[t->p[i]->index].Center[1] += dvol * sy;
    }


  if(t->p[j]->task == ThisTask && t->p[j]->index >= 0 && t->p[j]->index < N_gas)
    {
      SphP[t->p[j]->index].Volume += dvol;

      /* let's now compute the center-of-mass of the pyramid on top */
      sx = (2.0 / 3) * f->cx + (1.0 / 3) * t->p[j]->x;
      sy = (2.0 / 3) * f->cy + (1.0 / 3) * t->p[j]->y;

      SphP[t->p[j]->index].Center[0] += dvol * sx;
      SphP[t->p[j]->index].Center[1] += dvol * sy;
    }
}



void compute_circumcircles(void)
{
  int i;

  for(i = 0; i < Ndt; i++)
    {
      if(DT[i].p[0] == DPinfinity)
	continue;
      if(DT[i].p[1] == DPinfinity)
	continue;
      if(DT[i].p[2] == DPinfinity)
	continue;

      update_circumcircle(&DT[i]);
    }
}


void update_circumcircle(tetra * t)
{
  point *p0, *p1, *p2;

  p0 = t->p[0];
  p1 = t->p[1];
  p2 = t->p[2];

  if(t->p[0] == DPinfinity)
    return;
  if(t->p[1] == DPinfinity)
    return;
  if(t->p[2] == DPinfinity)
    return;

  double ax = p1->x - p0->x;
  double ay = p1->y - p0->y;

  double bx = p2->x - p0->x;
  double by = p2->y - p0->y;

  double aa = 0.5 * (ax * ax + ay * ay);
  double bb = 0.5 * (bx * bx + by * by);

  double mv_data[] = { ax, ay, aa, bx, by, bb };
  double x[2];

  int status = solve_linear_equations_2d(mv_data, x);

  if(status < 0)
    {
      printf("trouble in circum-circle calculation\n");
      endrun(17231388);
    }

  t->c.x = x[0] + p0->x;
  t->c.y = x[1] + p0->y;
  t->c.z = 0;
}




void get_circle_center(point * a, point * b, point * c, double *x, double *y)
{
  double s;

  s =
    0.5 * ((b->x - c->x) * (a->x - b->x) + (b->y - c->y) * (a->y - b->y))
    / ((a->y - c->y) * (a->x - b->x) - (a->x - c->x) * (a->y - b->y));

  *x = 0.5 * (a->x + c->x) + s * (a->y - c->y);
  *y = 0.5 * (a->y + c->y) - s * (a->x - c->x);
}



void set_integers_for_point(point * p)
{
  p->xx = (p->x - CentralOffsetX) * ConversionFac + 1.0;
  p->yy = (p->y - CentralOffsetY) * ConversionFac + 1.0;

  if(p->xx < 1.0 || p->xx >= 2.0 || p->yy < 1.0 || p->yy >= 2.0)
    {
      printf("%g %g %g\n", p->xx, p->yy, p->zz);
      endrun(13123);
    }

  p->ix = DOUBLE_to_VORONOIINT(p->xx);
  p->iy = DOUBLE_to_VORONOIINT(p->yy);


  unsigned long long *x;

  x = (unsigned long long *) &p->xx;
  *x = *x & (~((1llu << (52 - USEDBITS)) - 1));

  x = (unsigned long long *) &p->yy;
  *x = *x & (~((1llu << (52 - USEDBITS)) - 1));
}



void write_voronoi_mesh(char *fname, int writeTask, int lastTask)
{
  FILE *fd;
  MPI_Status status;
  int i, j, k, MaxNel, Nel;
  int ngas_tot, nel_tot, ndt_tot, nel_before, ndt_before, task;
  int *EdgeList, *Nedges, *NedgesOffset, *whichtetra;
  int *ngas_list, *nel_list, *ndt_list, *tmp;
  float *xyz_edges;
  tetra *q, *qstart;

  MaxNel = 10 * N_gas;		/* max edge list */
  Nel = 0;			/* length of edge list */

  EdgeList = mymalloc(MaxNel * sizeof(int));
  Nedges = mymalloc(N_gas * sizeof(int));
  NedgesOffset = mymalloc(N_gas * sizeof(int));
  whichtetra = mymalloc(N_gas * sizeof(int));
  xyz_edges = mymalloc(Ndt * DIMS * sizeof(float));
  ngas_list = mymalloc(sizeof(int) * NTask);
  nel_list = mymalloc(sizeof(int) * NTask);
  ndt_list = mymalloc(sizeof(int) * NTask);

  for(i = 0; i < Ndt; i++)
    {
      xyz_edges[i * DIMS + 0] = DT[i].c.x;
      xyz_edges[i * DIMS + 1] = DT[i].c.y;
    }

  for(i = 0; i < N_gas; i++)
    {
      Nedges[i] = 0;
      whichtetra[i] = -1;
    }

  for(i = 0; i < Ndt; i++)
    {
      for(j = 0; j < DIMS + 1; j++)
	if(DT[i].p[j]->task == ThisTask && DT[i].p[j]->index >= 0 && DT[i].p[j]->index < N_gas)
	  whichtetra[DT[i].p[j]->index] = i;
    }

  for(i = 0; i < N_gas; i++)
    {
      if(whichtetra[i] < 0)
	endrun(1212);

      qstart = q = &DT[whichtetra[i]];

      do
	{
	  Nedges[i]++;

	  if(Nel >= MaxNel)
	    endrun(7654);

	  EdgeList[Nel++] = q - DT;

	  for(j = 0; j < 3; j++)
	    if(q->p[j]->task == ThisTask && q->p[j]->index == i)
	      break;

	  k = j + 1;
	  if(k >= 3)
	    k -= 3;

	  q = q->t[k];
	}
      while(q != qstart);
    }

  for(i = 1, NedgesOffset[0] = 0; i < N_gas; i++)
    NedgesOffset[i] = NedgesOffset[i - 1] + Nedges[i - 1];


  /* determine particle numbers and number of edges in file */

  if(ThisTask == writeTask)
    {
      ngas_tot = N_gas;
      nel_tot = Nel;
      ndt_tot = Ndt;

      for(task = writeTask + 1; task <= lastTask; task++)
	{
	  MPI_Recv(&ngas_list[task], 1, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
	  MPI_Recv(&nel_list[task], 1, MPI_INT, task, TAG_LOCALN + 1, MPI_COMM_WORLD, &status);
	  MPI_Recv(&ndt_list[task], 1, MPI_INT, task, TAG_LOCALN + 2, MPI_COMM_WORLD, &status);

	  MPI_Send(&nel_tot, 1, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
	  MPI_Send(&ndt_tot, 1, MPI_INT, task, TAG_N + 1, MPI_COMM_WORLD);

	  ngas_tot += ngas_list[task];
	  nel_tot += nel_list[task];
	  ndt_tot += ndt_list[task];
	}

      if(!(fd = fopen(fname, "w")))
	{
	  printf("can't open file `%s' for writing snapshot.\n", fname);
	  endrun(123);
	}

      my_fwrite(&ngas_tot, sizeof(int), 1, fd);
      my_fwrite(&nel_tot, sizeof(int), 1, fd);
      my_fwrite(&ndt_tot, sizeof(int), 1, fd);

      my_fwrite(Nedges, sizeof(int), N_gas, fd);
      for(task = writeTask + 1; task <= lastTask; task++)
	{
	  tmp = mymalloc(sizeof(int) * ngas_list[task]);
	  MPI_Recv(tmp, ngas_list[task], MPI_INT, task, TAG_N + 2, MPI_COMM_WORLD, &status);
	  my_fwrite(tmp, sizeof(int), ngas_list[task], fd);
	  myfree(tmp);
	}

      my_fwrite(NedgesOffset, sizeof(int), N_gas, fd);
      for(task = writeTask + 1; task <= lastTask; task++)
	{
	  tmp = mymalloc(sizeof(int) * ngas_list[task]);
	  MPI_Recv(tmp, ngas_list[task], MPI_INT, task, TAG_N + 3, MPI_COMM_WORLD, &status);
	  my_fwrite(tmp, sizeof(int), ngas_list[task], fd);
	  myfree(tmp);
	}

      my_fwrite(EdgeList, sizeof(int), Nel, fd);
      for(task = writeTask + 1; task <= lastTask; task++)
	{
	  tmp = mymalloc(sizeof(int) * nel_list[task]);
	  MPI_Recv(tmp, nel_list[task], MPI_INT, task, TAG_N + 4, MPI_COMM_WORLD, &status);
	  my_fwrite(tmp, sizeof(int), nel_list[task], fd);
	  myfree(tmp);
	}

      my_fwrite(xyz_edges, sizeof(float), Ndt * DIMS, fd);
      for(task = writeTask + 1; task <= lastTask; task++)
	{
	  tmp = mymalloc(sizeof(float) * DIMS * ndt_list[task]);
	  MPI_Recv(tmp, sizeof(float) * DIMS * ndt_list[task], MPI_BYTE, task, TAG_N + 5, MPI_COMM_WORLD,
		   &status);
	  my_fwrite(tmp, sizeof(float), DIMS * ndt_list[task], fd);
	  myfree(tmp);
	}

      fclose(fd);
    }
  else
    {
      MPI_Send(&N_gas, 1, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
      MPI_Send(&Nel, 1, MPI_INT, writeTask, TAG_LOCALN + 1, MPI_COMM_WORLD);
      MPI_Send(&Ndt, 1, MPI_INT, writeTask, TAG_LOCALN + 2, MPI_COMM_WORLD);

      MPI_Recv(&nel_before, 1, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
      MPI_Recv(&ndt_before, 1, MPI_INT, writeTask, TAG_N + 1, MPI_COMM_WORLD, &status);

      for(i = 0; i < N_gas; i++)
	NedgesOffset[i] += nel_before;
      for(i = 0; i < Nel; i++)
	EdgeList[i] += ndt_before;

      MPI_Send(Nedges, N_gas, MPI_INT, writeTask, TAG_N + 2, MPI_COMM_WORLD);
      MPI_Send(NedgesOffset, N_gas, MPI_INT, writeTask, TAG_N + 3, MPI_COMM_WORLD);
      MPI_Send(EdgeList, Nel, MPI_INT, writeTask, TAG_N + 4, MPI_COMM_WORLD);
      MPI_Send(xyz_edges, sizeof(float) * DIMS * Ndt, MPI_BYTE, writeTask, TAG_N + 5, MPI_COMM_WORLD);
    }

  myfree(ndt_list);
  myfree(nel_list);
  myfree(ngas_list);
  myfree(xyz_edges);
  myfree(whichtetra);
  myfree(NedgesOffset);
  myfree(Nedges);
  myfree(EdgeList);

  if(ThisTask == 0)
    printf("wrote Voronoi mesh to file\n");
}


#endif

#endif
