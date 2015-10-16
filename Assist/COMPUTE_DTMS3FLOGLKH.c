#include <math.h>
#include "benchmark.h"


/**************************************************************
**  COMPUTE_DTMS3FLOGLKH--compute log likelihood xlk for fitted
*                         dichotomous model.
***************************************************************/
void 
  COMPUTE_DTMS3FLOGLKH (double *xlk, double ex, double Yp, double Yn)
{
  /* ex and w have to be used same prob_control */
  /* inorder to get same magnitude for xlk, xlkf and xlkr */
  /* if (ex > 0) *xlk += Yp*log(ex); */
  /* if (ex < 1) *xlk += Yn*log(1-ex); */

  if (ex <= 0.0)
    *xlk -= Yp * 1.0e+10;
  else if (ex >= 1.0)
    *xlk -= Yn * 1.0e+10;
  else
    *xlk += Yp * log (ex) + Yn * log (1.0 - ex);


}

