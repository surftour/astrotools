/* interpolation/accel.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Author:  G. Jungman
 */
/*
#include <config.h>
*/
#include <stdlib.h>
#include <stdio.h>
#include "manual_gsl_interp.h"

gsl_interp_accel *
gsl_interp_accel_alloc (void)
{
  gsl_interp_accel *a = (gsl_interp_accel *) malloc (sizeof (gsl_interp_accel));
  if (a == 0)
    {
      printf("could not allocate space for gsl_interp_accel\n");
    }

  a->cache = 0;
  a->hit_count = 0;
  a->miss_count = 0;

  return a;
}

int
gsl_interp_accel_reset (gsl_interp_accel * a)
{
  a->cache = 0;
  a->hit_count = 0;
  a->miss_count = 0;

  return 0;
}

#ifndef HIDE_INLINE_STATIC
int
gsl_interp_accel_find (gsl_interp_accel * a, const double xa[], int len, double x)
{
  int x_index = a->cache;

  if (x < xa[x_index])
    {
      a->miss_count++;
      a->cache = gsl_interp_bsearch (xa, x, 0, x_index);
    }
  else if (x > xa[x_index + 1])
    {
      a->miss_count++;
      a->cache = gsl_interp_bsearch (xa, x, x_index, len - 1);
    }
  else
    {
      a->hit_count++;
    }

  return a->cache;
}
#endif

void
gsl_interp_accel_free (gsl_interp_accel * a)
{
  free (a);
}
