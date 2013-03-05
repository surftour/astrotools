#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "proto.h"


int NumPart;

struct particle
{
  float Pos[2];
} *P;


float *Hsml;
float *Quantity;
float *Mass;
float *Rho;

float Xmin, Ymin, Xmax, Ymax;

float Xc, Yc;

int Xpixels, Ypixels;

float *Value, *ValueQuantity;


int slice(int argc, void *argv[])
{
  if(argc != 13)
    {
      fprintf(stderr, "\n\nwrong number of arguments ! (found %d) \n\n", argc);
      exit(0);
    }

  NumPart = *(int *) argv[0];
  P = (struct particle *) argv[1];
  Hsml = (float *) argv[2];
  Mass = (float *) argv[3];
  Quantity = (float *) argv[4];


  Xmin = *(float *) argv[5];
  Xmax = *(float *) argv[6];
  Ymin = *(float *) argv[7];
  Ymax = *(float *) argv[8];

  Xpixels = *(int *) argv[9];
  Ypixels = *(int *) argv[10];

  Value =         (float *) argv[11];
  ValueQuantity = (float *) argv[12];

  printf("N=%d\n", NumPart);

  make_map();

  return 0;
}



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



  for(i = 0; i < Xpixels; i++)
    for(j = 0; j < Ypixels; j++)
      {
        Value[i * Ypixels + j] = 0;
        ValueQuantity[i * Ypixels + j] = 0;
      }


  LengthX = Xmax-Xmin;
  LengthY = Ymax-Ymin;

  pixelsizeX = LengthX / Xpixels;
  pixelsizeY = LengthY / Ypixels;

  if(pixelsizeX < pixelsizeY)
    hmin = 1.001 * pixelsizeX / 2;
  else
    hmin = 1.001 * pixelsizeY / 2;

  if(pixelsizeX < pixelsizeY)
    hmax = 64 * pixelsizeX;
  else
    hmax = 64 * pixelsizeY;



  for(n = 0; n < NumPart; n++)
    {
      if((n % (NumPart / 100)) == 0)
	{
	  printf(".");
	  fflush(stdout);
	}

      pos[0]= P[n].Pos[0]-Xmin;
      pos[1]= P[n].Pos[1]-Ymin;


      h = Hsml[n];
      
      if(h < hmin)
        h = hmin;

      if(h > hmax)
        h = hmax;

      if(pos[0] + h < 0 || pos[0] - h >  LengthX
         || pos[1] + h < 0 || pos[1] - h > LengthY)
        continue;
 
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

  printf("\n");
}

