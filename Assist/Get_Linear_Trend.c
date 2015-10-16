#include "benchmark.h"

int 
  Get_Linear_Trend (int N, double *doses, double *means, int *numi)
/***********************************************************************
* 
* This function will determine the general linear trend of the data.
*	
*	INPUT:	N - the number of dose levels
*			doses[i] - the ith dose
*			means[i] - the sample mean at the ith dose
*			numi[i] - the number individuals in the ith dose group
*
*	Returns:	 1 if the trend is positive or there is no trend
*				-1 if the trend is negative
*
************************************************************************/
{
  int default_dir, i, NTotal;
  double Xbar, Ybar, temp, temp2, slope;

  Xbar = 0.0;
  temp = 0.0;
  temp2 = 0.0;
  Ybar = 0.0;
  NTotal = 0;

  for (i = 1; i <= N; i++)
    {
      Xbar += doses[i] * numi[i];
      Ybar += means[i] * numi[i];
      NTotal += numi[i];
    }				/* end for */

  Xbar = Xbar / NTotal;
  Ybar = Ybar / NTotal;

  for (i = 1; i <= N; i++)
    {
      temp += numi[i] * (means[i] - Ybar) * (doses[i] - Xbar);
      temp2 += numi[i] * (doses[i] - Xbar) * (doses[i] - Xbar);
    }				/* end for */

  slope = temp / temp2;

  if (slope >= 0)
    {
      default_dir = 1;
    }
  else
    {
      default_dir = -1;
    }				/* end if */

  return default_dir;

}				/* End Get_Linear_Trend() */

