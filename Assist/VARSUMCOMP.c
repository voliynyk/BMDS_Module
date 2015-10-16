#include <math.h>
#include "benchmark.h"


/**********************************************************************
** VARSUMCOMP--subfunction used to compute summary data for a variable.
***********************************************************************/
void 
  VARSUMCOMP (int n, int Nobs, double v[], VarList varsum[])
{
  int i;

  varsum[n].flag = n;
  varsum[n].S = 0;
  varsum[n].SS = 0;

  for (i = 1; i <= Nobs; i++)
    {
      varsum[n].S += v[i];
      varsum[n].SS += v[i] * v[i];
    }

  varsum[n].M = varsum[n].S / Nobs;
  varsum[n].CSS = varsum[n].SS - varsum[n].M * varsum[n].S;

}				/*end of VARSUMCOMP subfunction */


