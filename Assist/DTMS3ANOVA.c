#include "benchmark.h"
#include "ERRORPRT.h"
#include "specialfun.h"
#include "computation.h"



/********************************************************
**  DTMS3ANOVA--compute element values for ANOVA table
*               of a dichotomous model with 3 parameters.
*				If the Chi-Square test is invalid, then
*				anasum[i].TEST < 0 upon return
*********************************************************/

extern int Model_DF(int []); 
void 
  DTMS3ANOVA (int nparm, int Nobs, int Spec[], double xlkf, double xlk, double xlkr, AnaList anasum[], int bounded[])
{
  /*compute ANOVA table:
     anasum[1]--Full model
     anasum[2]--Fitted model
     anasum[3]--Reduced model */

  double pv, dev;
  int df;

  anasum[1].DF = COMPUTEDF (nparm, Spec) + 1;
  anasum[1].SS = xlkf;

  dev = 2 * (xlkf - xlk);
  df = Nobs - Model_DF(bounded);
  
  if (dev < 0.0)
    {
      Warning (" Warning: Likelihood for the fitted model larger than the Likelihood for the full model.");
      pv = -1;
    }
  else
    {
      if (df > 0)
	{
	  pv = CHISQ (dev, df);
	}
      else
	{
	  pv = -1;
	}
    }
  anasum[2].SS = xlk;
  anasum[2].MSE = dev;
  anasum[2].DF = df;

  anasum[2].TEST = pv;

  if (Nobs > 1)
    {
      dev = 2 * (xlkf - xlkr);
      if (dev <= 0)
	{
	  Warning ("\n Warning: Likelihood for the full model smaller than the Likelihood for the reduced model.");
	  pv = -1;
	}
      else
	{
	  pv = CHISQ (dev, Nobs - 1);
	}
      anasum[3].SS = xlkr;
      anasum[3].MSE = dev;
      anasum[3].DF = Nobs - 1;
    }
  else
    {
      pv = -9999;
    }
  anasum[3].TEST = pv;

}
