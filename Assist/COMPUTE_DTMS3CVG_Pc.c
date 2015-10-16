#include <math.h>
#include "benchmark.h"


/**************************************************************
**  COMPUTE_DTMS3CVG_Pc--compute parameter convergence variable
*                        Pc for 3-parameters dichotomous model.
*************************************************************/
void 
  COMPUTE_DTMS3CVG_Pc (int nparm, double **bb, double bp[])
{
  int i;
  double x;

  CVG.Pc = 0;
  for (i = 1; i <= nparm; i++)
    {
      x = fabs (bb[i][1] - bp[i]) / (fabs (bb[i][1]) + 0.00001);
      if (x > CVG.Pc)
	CVG.Pc = x;
      bb[i][1] = bp[i];
    }
}

