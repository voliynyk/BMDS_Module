#include <math.h>
#include "benchmark.h"


/*******************************************************************
**  COMPUTE_DTMS3LOGLKH--compute log likelihood for full and reduced
*                        model xlkf and xlkr for dichotomous models.
****************************************************************/
void 
COMPUTE_DTMS3LOGLKH (int Nobs, double *xlkf, double *xlkr, VarList varsum[],
		     double Yp[], double Yn[])
{
  int i;
  double W;

  *xlkf = 0.0;
  varsum[1].S = 0;
  varsum[2].S = 0;
  for (i = 1; i <= Nobs; i++)
    {
      varsum[1].S += Yp[i];
      varsum[2].S += Yn[i];
      W = Yp[i] / (Yp[i] + Yn[i]);
      if (W > 0)
	*xlkf += Yp[i] * log (W);
      if (W < 1)
	*xlkf += Yn[i] * log (1 - W);
      /* W=0 <=> Yp[i]=0 <=> log(1-W)=0 <=> xlkf=0 */
      /* W=1 <=> Yn[i]=0 <=> log(W)=0 <=> xlkf=0 */
      /* PROBABILITY_W(&W); */
      /* *xlkf += Yp[i] * log(W) + Yn[i] * log(1- W); */
    }
  W = varsum[1].S / (varsum[1].S + varsum[2].S);
  *xlkr = varsum[1].S * log (W) + varsum[2].S * log (1 - W);
}


