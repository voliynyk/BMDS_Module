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
#include <benchmark.h>
#include <ERRORPRT.h>
#include <allo_memo.h>

int    N;      
double *int_X;        /* internal independent var. array */
double *int_Y;        /* internal dependent var. array   */
double *int_p;        /* internal parmeter array         */

/* SumSquares(c): a function that returns the sum of squares */
/* for the power model: Y = a + b*X^c                        */
double SumSquares (double *c);

void initialparms(int nobs, double min, double max, double tmX[],
		 double tmY[], double P[], double gtol)
{

  int i;
  double junk;

  N = nobs;
  int_X = DVECTOR(1, nobs);
  int_Y = DVECTOR(1, nobs);
  int_p = DVECTOR(1, 5);

  for (i = 1; i <= nobs; i++)
    {
      int_X[i] = tmX[i];
      int_Y[i] = tmY[i];
    }
  for (i = 1; i <= 5; i++)
    int_p[i] = P[i];


  P[5] = fmin_(&min, &max, SumSquares, &gtol);

  junk = SumSquares(&P[5]);


  P[1] = int_p[1];
  P[2] = int_p[2];



  FREE_DVECTOR(int_X, 1, nobs);
  FREE_DVECTOR(int_Y, 1, nobs);
  FREE_DVECTOR(int_p, 1, 5);

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
  printf("N=%d",N);

  newX = DVECTOR(1, N);
  meanX = meanY = sumXY = sumX2 = 0.0;
  for (i=1; i<=N; i++)
    {
      newX[i] = pow(int_X[i], *c);
      meanX += newX[i];
      meanY += int_Y[i];
      sumXY += newX[i]*int_Y[i];
      sumX2 += pow(newX[i], 2);
    }
  meanX = meanX/N;
  meanY = meanY/N;
  b = (sumXY - N*meanX*meanY)/(sumX2 - N*pow(meanX,2));
  a = meanY - b*meanX;
  int_p[1] = a;
  int_p[2] = b;

  /* This part finds the Sum of Squares Error using the parms *
   * found in the previous part.                              */
  SSE = 0;
  for (i = 1; i <= N; i++)
    {
      Yhat = a + b*newX[i];
      SSE += (int_Y[i]-Yhat)*(int_Y[i]-Yhat);
    }

  FREE_DVECTOR(newX, 1, N);
  return SSE;
}
  
