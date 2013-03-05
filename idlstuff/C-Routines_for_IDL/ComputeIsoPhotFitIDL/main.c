#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "proto.h"
#include "nrsrc/nrutil.h"


#define PI  3.141592
#define N   360
#define MAXITER 100



int Size;
char fname_input[100];

double Th, Scale;

float *Xfit, *Yfit;
float *Xcon, *Ycon;
float *Result_Phi, *Result_A, *Result_B, *Result_X0, *Result_Y0, *Result_a4;
float *Rawdata;

double **value;

double Phi, A, B, X0, Y0, a4;
double *t, *ts, *tc, *phi, *xx, *yy, *xi, *yi;

double dx, dy;

static int flag;

int fit_isophot(int argc, char **argv)
{
  int i, j;
  float *fp;
  float min, max;
  
  flag=0;

  if(argc != 14)
    {
      fprintf(stderr, "\n\nwrong number of arguments ! (found %d) \n\n", argc);
      exit(0);
    }

  Size = *(int *) argv[0];
  Rawdata = (float *) argv[1];
  Th = *(float *) argv[2];
  Scale = *(float *) argv[3];

  Xfit = (float *) argv[4];	/* must give room for a result of 72 points */
  Yfit = (float *) argv[5];

  Result_Phi = (float *) argv[6];
  Result_A = (float *) argv[7];
  Result_B = (float *) argv[8];
  Result_X0 = (float *) argv[9];
  Result_Y0 = (float *) argv[10];
  Result_a4 = (float *) argv[11];

  Xcon = (float *) argv[12];
  Ycon = (float *) argv[13];


  for(i = 0, min = max = *Rawdata, fp = Rawdata; i < Size; i++)
    for(j = 0; j < Size; j++)
      {
	if(*fp < min)
	  min = *fp;
	if(*fp > max)
	  max = *fp;
	fp++;
      }
  printf("min= %g\n", min);
  printf("max= %g\n", max); fflush(stdout);

  value = dmatrix(-1, Size, -1, Size);

  for(i = -1; i <= Size; i++)
    for(j = -1; j <= Size; j++)
      {
	value[j][i] = min;
      }

  for(i = 0, fp = Rawdata; i < Size; i++)
    for(j = 0; j < Size; j++)
      {
	value[j][i] = *fp;
	fp++;
      }



  fit_ellipse();

  free_dmatrix(value, -1, Size, -1, Size - 1);

  return 0;
}









void fit_ellipse(void)
{
  int i;
  double ps, pc;
  double sum1, sum2, sum3, sum4, chi2, chi2old, dx, dy;
  int iter = 0;

  t = dvector(1, N);
  ts = dvector(1, N);
  tc = dvector(1, N);
  phi = dvector(1, N);
  xx = dvector(1, N);
  yy = dvector(1, N);
  xi = dvector(1, N);
  yi = dvector(1, N);

  for(i = 1; i <= N; i++)
    {
      t[i] = (2 * PI / N) * (i - 1);
      ts[i] = sin(t[i]);
      tc[i] = cos(t[i]);
    }


  /* set initial guess */

  Phi = 0.0;
  A = B = 1.0;
  X0 = Y0 = 0;

  chi2old = 1e30;


  do
    {
      iter++;

      if(flag==0)
	compute_isophotal_points();

      compute_ideal_points();

      for(i = 1, X0 = 0; i <= N; i++)
	X0 += xx[i];
      for(i = 1, Y0 = 0; i <= N; i++)
	Y0 += yy[i];
      X0 /= N;
      Y0 /= N;


      ps = sin(Phi);
      pc = cos(Phi);


      for(i = 1, sum1 = sum2 = 0; i <= N; i++)
	{
	  sum1 += xx[i] * ts[i] + yy[i] * tc[i];
	  sum2 += xx[i] * tc[i] - yy[i] * ts[i];
	}

      Phi = 0.5 * atan(sum1 / sum2);

      for(i = 1, sum1 = sum2 = sum3 = sum4 = 0; i <= N; i++)
	{
	  sum1 +=
	    (xx[i] + B * sin(t[i] - Phi) * ps - X0) * cos(t[i] - Phi) * pc + (yy[i] -
									      B * sin(t[i] - Phi) * pc -
									      Y0) * cos(t[i] - Phi) * ps;
	  sum3 += cos(t[i] - Phi) * cos(t[i] - Phi);

	  sum2 +=
	    -(xx[i] - A * cos(t[i] - Phi) * pc - X0) * sin(t[i] - Phi) * ps + (yy[i] -
									       A * cos(t[i] - Phi) * ps -
									       Y0) * sin(t[i] - Phi) * pc;
	  sum4 += sin(t[i] - Phi) * sin(t[i] - Phi);
	}

      A = sum1 / sum3;
      B = sum2 / sum4;

      compute_ideal_points();

      for(i = 1, chi2 = 0; i <= N; i++)
	{
	  dx = xx[i] - xi[i];
	  dy = yy[i] - yi[i];

	  chi2 += dx * dx + dy * dy;
	}

      if(fabs(chi2 - chi2old) / chi2 < 1e-8)
	break;

      if(iter > MAXITER)
	{
	  printf("maximum iteration count reached\n");
	  break;
	}

      chi2old = chi2;
    }
  while(1);

  if(flag==0)
    compute_isophotal_points();

  compute_ideal_points();

  for(i = 1; i <= N; i++)
    {
      Xfit[i - 1] = xi[i];
      Yfit[i - 1] = yi[i];

      Xcon[i - 1] = xx[i];
      Ycon[i - 1] = yy[i];
    }


  /* write(); */

  a4 = compute_a(4);

  *Result_Phi = Phi;
  *Result_A = fabs(A);
  *Result_B = fabs(B);
  *Result_X0 = X0;
  *Result_Y0 = Y0;
  *Result_a4 = a4;


  printf("A=%g  B=%g  X0=%g  Y0=%g  Phi=%g  a4/A=%g  Chi2=%g\n\n", A, B, X0, Y0, Phi, a4/A, chi2);
  fflush(stdout);


  free_dvector(t, 1, N);
  free_dvector(ts, 1, N);
  free_dvector(tc, 1, N);
  free_dvector(phi, 1, N);
  free_dvector(xx, 1, N);
  free_dvector(yy, 1, N);
  free_dvector(xi, 1, N);
  free_dvector(yi, 1, N);
}



double compute_a(int k)
{
  int i;
  double sum, dr;

  for(i = 1, sum = 0; i <= N; i++)
    {

      dr = sqrt((xx[i] - X0) * (xx[i] - X0) + (yy[i] - Y0) * (yy[i] - Y0))
	- sqrt((xi[i] - X0) * (xi[i] - X0) + (yi[i] - Y0) * (yi[i] - Y0));

      sum += dr * cos(k * (t[i] - Phi));
    }

  return sum / (N / 2);
}



void compute_ideal_points(void)
{
  int i;

  for(i = 1; i <= N; i++)
    {
      xi[i] = X0 + A * cos(t[i] - Phi) * cos(Phi) - B * sin(t[i] - Phi) * sin(Phi);
      yi[i] = Y0 + A * cos(t[i] - Phi) * sin(Phi) + B * sin(t[i] - Phi) * cos(Phi);
    }
}


void compute_isophotal_points(void)
{
  int i, j, errflag;
  double tt, tmax, tb;
  float x, y, rm;

  /*** below a test suite ***/

  /*
     float AFit, BFit, XFit, YFit, PhiFit;
     float dr, r; 

     AFit=1.0;
     BFit=0.5;
     XFit=0.0;
     YFit=0.0;
     PhiFit=PI/6;

     a4=0.04;

     for(i=1;i<=N;i++)
     {
     xx[i]=XFit+AFit*cos(t[i]-PhiFit)*cos(PhiFit) - BFit*sin(t[i]-PhiFit)*sin(PhiFit);
     yy[i]=YFit+AFit*cos(t[i]-PhiFit)*sin(PhiFit) + BFit*sin(t[i]-PhiFit)*cos(PhiFit);

     r=  sqrt(xx[i]*xx[i] + yy[i]*yy[i]);

     dr = a4*cos((t[i]-PhiFit)*4);

     xx[i]+=  xx[i]/r * dr;
     yy[i]+=  yy[i]/r * dr;
     }

     return 0; 
   */


  for(i = 1; i < N; i++)
    {
      x = A * cos(t[i] - Phi);
      y = B * sin(t[i] - Phi);

      phi[i] = arg(x / sqrt(x * x + y * y), y / sqrt(x * x + y * y)) + Phi;

      /*      printf("%g %g\n",phi[i],t[i]); */
    }



  if(fabs(X0) > fabs(Y0))
    rm = fabs(X0);
  else
    rm = fabs(Y0);


  for(i = 1; i <= N; i++)
    {
      dx = cos(phi[i]);
      dy = sin(phi[i]);

      /*
         tmax=Scale/2-0.6*Scale/Size-rm;
       */

      tmax = Scale / 2;

      for(j = 1; j <= 100; j++)
	{
	  tb = j * tmax / 100.0;
	  if(ray(tb) < 0)
	    break;
	}

      if(j>100)
	{
	  tt = 1.5*Scale;

	  for(i = 1; i <= N; i++)
	    {
	      dx = cos(phi[i]);
	      dy = sin(phi[i]);
	      
	      xx[i] = X0 + dx * tt;
	      yy[i] = Y0 + dy * tt;
	    }

	  printf("outside!\n");
	  flag = 1;
	  return;
	}
      else
	{
	  tt = zriddr(ray, 0, tb, 1e-8, &errflag);
	}

      if(errflag)
	{
	  printf("warning.\n");
	  exit(3);
	}

      xx[i] = X0 + dx * tt;
      yy[i] = Y0 + dy * tt;

      /*      printf("i:%d      x: %g   y: %g\n",i,xx[i],yy[i]); */

    }
}





void write(void)
{
  FILE *fd;
  int i;

  fd = fopen("test.txt", "w");
  for(i = 1; i <= N; i++)
    fprintf(fd, "%g %g\n", xx[i], yy[i]);
  for(i = 1; i <= N; i++)
    fprintf(fd, "%g %g\n", xi[i], yi[i]);
  fclose(fd);

}





double arg(double co, double si)
{
  if(si > 0)
    return (acos(co));
  else
    return (2 * PI - acos(co));
}



double ray(double t)
{
  double x, y;

  x = X0 + dx * t;
  y = Y0 + dy * t;

  return sb(x, y) - Th;
}




double sb(double x, double y)	/* returns surface brightness */
{
  int i, j;
  double u, v;
  double xold, yold;

  xold = x;
  yold = y;



  x += Scale / 2 - 0.5 * Scale / Size;
  y += Scale / 2 - 0.5 * Scale / Size;

  x *= Size / Scale;
  y *= Size / Scale;

  x++;
  y++;

  i = (int) x;
  j = (int) y;

  x--;
  y--;
  i--;
  j--;

  u = x - i;
  v = y - j;



  if(x >= (Size) || x < -1)
    {
      if(x >= (Size))
	{
	  i = (Size - 1);
	  u = x - i;
	}
      if(x < -1)
	{
	  i = -1;
	  u = x - i;
	}
    }
  if(y >= (Size) || y < -1)
    {
      if(y >= (Size))
	{
	  j = (Size - 1);
	  v = y - j;
	}
      if(y < -1)
	{
	  j = 0;
	  v = y - j;
	}
    }



  return ((1 - u) * (1 - v) * value[i][j] +
	  (u) * (1 - v) * value[i + 1][j] +
	  (1 - u) * (v) * value[i][j + 1] + (u) * (v) * value[i + 1][j + 1]);
}
