/********************************************************
 *  initialparms--compute initial values for alpha, beta, 
 *  and rho for the NCTR and RaiVR models.
 *
 *  Input:
 *    nobs = number of observations (int)
 *    min = the lower bound on rho (double)
 *    max = the upper bound on rho (double)
 *    X[] = array of independent variables (double)
 *    Y[] = array of dependent variables (double)
 *    Parms[] = array of alpha, beta, and rho (double)
 *    EPS = convergence critria
 *    
 *  By: Micheal Ferree
 *  Date: May 23, 2001
 *
 ********************************************************/
#define ASSIST
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <benchmark.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include <ERRORPRT.h>
#include <allo_memo.h>
#include <matrix_agb.h>
#include <specialfun.h>
#include <computation.h>
#include <in_outfun.h>

int    Nobs;
double *X;
double *Y;
double *p;
/* SumSquares(c): a function that returns the sum of squares */
/* for the power model: Y = a + b*X^c                        */
double SumSquares (double *c);

void initialparms(int nobs, double min, double max, double tmX[],
		 double tmY[], double Parms[], double gtol)
{

  int i;
  double junk;

  Nobs = nobs;
  X = DVECTOR(1, nobs);
  Y = DVECTOR(1, nobs);
  p = DVECTOR(1, 5);

  for (i = 1; i <= nobs; i++)
    {
      X[i] = tmX[i];
      Y[i] = tmY[i];
    }
  for (i = 1; i <= 5; i++)
    p[i] = Parms[i];

  Parms[5] = fmin_(&min, &max, SumSquares, &gtol);
  junk = SumSquares(&Parms[5]);
  Parms[1] = p[1];
  Parms[2] = p[2];

  FREE_DVECTOR(X, 1, nobs);
  FREE_DVECTOR(Y, 1, nobs);
  FREE_DVECTOR(p, 1, 5);

  return;
}

/*  This function takes a given power(c) and finds the regression coeffients */
/*  for a and b in the model: Y = a + b*X^c and returns the SSE              */ 

double SumSquares(double *c)
{
  int i;
  double a, b;
  double *newX;
  double meanX, meanY, sumXY, sumX2;
  double Yhat, SSE;

  /* This part finds the least squares line using rho=c */
  newX = DVECTOR(1, Nobs);
  meanX = meanY = sumXY = sumX2 = 0.0;
  for (i=1; i<=Nobs; i++)
    {
      newX[i] = pow(X[i], *c);
      meanX += newX[i];
      meanY += Y[i];
      sumXY += newX[i]*Y[i];
      sumX2 += pow(newX[i], 2);
    }
  meanX = meanX/Nobs;
  meanY = meanY/Nobs;
  b = (sumXY - Nobs*meanX*meanY)/(sumX2 - Nobs*pow(meanX,2));
  a = meanY - b*meanX;
  p[1] = a;
  p[2] = b;

  /* This part finds the Sum of Squares Error using the parms *
   * found in the previous part.                              */
  SSE = 0;
  for (i = 1; i <= Nobs; i++)
    {
      Yhat = a + b*newX[i];
      SSE += (Y[i]-Yhat)*(Y[i]-Yhat);
    }

  FREE_DVECTOR(newX, 1, Nobs);
  return SSE;
}
  
