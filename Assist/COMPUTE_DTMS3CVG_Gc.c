#include "benchmark.h"


/**************************************************************
**  COMPUTE_DTMS3CVG_Gc--compute gradient convergence variable
*                        Gc for 3-parameters dichotomous model.
************************************************************/
void 
  COMPUTE_DTMS3CVG_Gc (int nparm, double **y, double **d)
{
  int i;

  CVG.Gc = 0.0;
  for (i = 1; i <= nparm; i++)
    CVG.Gc += y[i][1] * d[i][1];
}


