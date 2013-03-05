#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "proto.h"
#include "allvars.h"


int findHsmlAndProject(int argc, void *argv[])
{
  if(argc != 21)
    {
      fprintf(stderr, "\n\nwrong number of arguments ! (found %d) \n\n", argc);
      exit(0);
    }

  NumPart = *(int *) argv[0];
  P = (struct particle_data *) argv[1];
  Hsml = (float *) argv[2];
  Mass = (float *) argv[3];
  Quantity = (float *) argv[4];


  Xmin = *(float *) argv[5];
  Xmax = *(float *) argv[6];
  Ymin = *(float *) argv[7];
  Ymax = *(float *) argv[8];
  Zmin = *(float *) argv[9];
  Zmax = *(float *) argv[10];

  Xpixels = *(int *) argv[11];
  Ypixels = *(int *) argv[12];

  DesDensNgb =  *(int *) argv[13];

  Axis1 = *(int *) argv[14];
  Axis2 = *(int *) argv[15];
  Axis3 = *(int *) argv[16];

  Hmax = *(float *) argv[17];

  BoxSize = *(double *) argv[18];

  Value =         (float *) argv[19];
  ValueQuantity = (float *) argv[20];

  printf("N=%d\n", NumPart);

  
  printf("peano-hilbert order...\n");
  peano_hilbert_order();
  printf("done");

  tree_treeallocate(2.0 * NumPart, NumPart);
  Softening = (Xmax-Xmin)/Xpixels/100;

  printf("build tree...\n");
  tree_treebuild();
  printf("done.\n");

  printf("finding neighbours...\n");
  determine_hsml();
  printf("done.\n");

  tree_treefree();

  printf("projecting\n");
  make_map();
  printf("done\n");

  return 0;
}

void determine_hsml(void)
{
  int i, signal;
  double h;

  for(i = 0, signal = 0, h = 0; i < NumPart; i++)
    {
      if(i > (signal / 100.0) * NumPart)
        {
          printf("x");
          fflush(stdout);
          signal++;
        }

      if(Hsml[i] == 0)
        Hsml[i] = h  = ngb_treefind(P[i].Pos, DesDensNgb, h * 1.1);
    }

  printf("\n");
}


#ifdef PERIODIC
#define NGB_PERIODIC(x) (((x)>BoxHalf)?((x)-BoxSize):(((x)<-BoxHalf)?((x)+BoxSize):(x)))
#else
#define NGB_PERIODIC(x) (x)
#endif



void make_map(void)
{
  int i, j, n;
  int dx, dy, nx, ny;
  double h, r, u, wk;
  double pos[2];
  double LengthX;
  double LengthY;
  double r2, h2;
  double sum, hmin, hmax, x, y, xx, yy, xxx, yyy;
  double pixelsizeX, pixelsizeY;

  BoxHalf = 0.5 * BoxSize;

  for(i = 0; i < Xpixels; i++)
    for(j = 0; j < Ypixels; j++)
      {
        Value[i * Ypixels + j] = 0;
        ValueQuantity[i * Ypixels + j] = 0;
      }
 

  LengthX = Xmax-Xmin;
  LengthY = Ymax-Ymin;

  Xc = 0.5 * (Xmax + Xmin);
  Yc = 0.5 * (Ymax + Ymin);

  pixelsizeX = LengthX / Xpixels;
  pixelsizeY = LengthY / Ypixels;

  if(pixelsizeX < pixelsizeY)
    hmin = 1.001 * pixelsizeX / 2;
  else
    hmin = 1.001 * pixelsizeY / 2;


  hmax = Hmax;


  for(n = 0; n < NumPart; n++)
    {
      if((n % (NumPart / 100)) == 0)
	{
	  printf(".");
	  fflush(stdout);
	}

      if(P[n].Pos[Axis3]< Zmin || P[n].Pos[Axis3] > Zmax)
        continue;

      pos[0]= P[n].Pos[Axis1]-Xmin;
      pos[1]= P[n].Pos[Axis2]-Ymin;


      h = Hsml[n];
      
      if(h < hmin)
        h = hmin;

      if(h > hmax)
        h = hmax;

#ifdef PERIODIC
      if((NGB_PERIODIC(Xc - P[n].Pos[Axis1]) + 0.5 * (Xmax - Xmin)) < -Hsml[n])
        continue;
      if((NGB_PERIODIC(Xc - P[n].Pos[Axis1]) - 0.5 * (Xmax - Xmin)) > Hsml[n])
        continue;
      
      if((NGB_PERIODIC(Yc - P[n].Pos[Axis2]) + 0.5 * (Ymax - Ymin)) < -Hsml[n])
        continue;
      if((NGB_PERIODIC(Yc - P[n].Pos[Axis2]) - 0.5 * (Ymax - Ymin)) > Hsml[n])
        continue;
#else
      if(pos[0] + h < 0 || pos[0] - h >  LengthX
         || pos[1] + h < 0 || pos[1] - h > LengthY)
        continue;
#endif
 
      h2 = h * h;

      nx = h / pixelsizeX + 1;
      ny = h / pixelsizeY + 1;

      /* x,y central pixel of region covered by the particle on the mesh */
      
      x = (floor(pos[0] / pixelsizeX) + 0.5) * pixelsizeX;
      y = (floor(pos[1] / pixelsizeY) + 0.5) * pixelsizeY;

      /* determine kernel normalizaton */

      
      sum = 0;

      
      for(dx = -nx; dx <= nx; dx++)
        for(dy = -ny; dy <= ny; dy++)
          {
            xx = x + dx * pixelsizeX - pos[0];
            yy = y + dy * pixelsizeY - pos[1];
            r2 = xx * xx + yy * yy;
            
            if(r2 < h2)
              {
                r = sqrt(r2);
                u = r / h;
                
                if(u < 0.5)
                  wk = (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
                else
                  wk = 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                
                sum += wk;
              }
          }
      
      if(sum < 1.0e-10)
        continue;

      for(dx = -nx; dx <= nx; dx++)
        for(dy = -ny; dy <= ny; dy++)
          {
            xxx = x + dx * pixelsizeX;
            yyy = y + dy * pixelsizeY;

#ifdef PERIODIC
            xxx = NGB_PERIODIC(xxx + Xmin - Xc) + Xc - Xmin;
            yyy = NGB_PERIODIC(yyy + Ymin - Yc) + Yc - Ymin;
#endif
           
            if(xxx >= 0 && yyy >= 0)
              {
                i = xxx / pixelsizeX;
                j = yyy / pixelsizeY;
                
                if(i >= 0 && i < Xpixels)
                  if(j >= 0 && j < Ypixels)
                    {
                      xx = x + dx * pixelsizeX - pos[0];
                      yy = y + dy * pixelsizeY - pos[1];
                      r2 = xx * xx + yy * yy;
                      
                      if(r2 < h2)
                        {
                          r = sqrt(r2);
                          u = r / h;
                          
                          if(u < 0.5)
                            wk = (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
                          else
                            wk = 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                          
                          Value[i * Ypixels + j] += Mass[n] * wk / sum;
                          ValueQuantity[i * Ypixels + j] += Mass[n]*Quantity[n]*wk / sum;
                        }
                    }
              }
          }
    }


  for(i = 0; i < Xpixels; i++)
    for(j = 0; j < Ypixels; j++)
      if(Value[i * Ypixels + j]>0)
        ValueQuantity[i * Ypixels + j] /= Value[i * Ypixels + j];


  printf("\n");
}

