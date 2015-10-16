/******************************************************************* */
/* Goodness.c -- Quantal Goodness-of-fit function */
/* Feb. 7, 2000 */
/******************************************************************* */

/******************************************************************* */
/* Quantal_Goodness.c  Compute goodness-of-fit statistics for */
/*                     Quantal models. */
/*               Global (declared in benchmark.h) : */
/*               Xi, Yn */
/*		input: */
/*			nparm is the number of parameters in the model */
/*			bounded is the vector giving which parameters */
/*                       are fixed or bounded */
/*			Parms[] is the parameter vector */
/*				of length nparm */
/*                       Xi[] vector of doses */
/*                       Yp[] vector of positive responses */
/*                       Yn[] vector of negative responses */
/* */
/************************************************************/
#include <stdlib.h>
#include <math.h>
#include "benchmark.h"
#include "specialfun.h"
#include "allo_memo.h"
#include "in_outfun.h"

extern int Model_DF (int []);
extern void Predict (double [], int, double [], double []);

void Quantal_Goodness(int nparm, int bounded[], double Parms[], int Nobs,
		      double Xi[], double Yp[], double Yn[], double scale)
{
  double chs, pv, x, t, *Ep, *Ypp, resid;
  int  i, df;

  Ep = DVECTOR(1, Nobs);
  Ypp = DVECTOR(1, Nobs);
  Predict(Xi, Nobs, Parms, Ep);
  for (i=1;i<=Nobs;i++)
    {
      Ypp[i] = Ep[i]*(Yn[i]+Yp[i]);
    }

  /*Output Goodness of fit table*/
  OUTPUT_TEXT("\n\n                                  Goodness  of  Fit ");
  OUTPUT_TEXT(    "                                                                 Scaled");
  OUTPUT_TEXT(    "     Dose     Est._Prob.    Expected    Observed     Size       Residual"); 
  OUTPUT_TEXT(    "  ------------------------------------------------------------------------");
  chs = 0.0;   

  for (i=1;i<=Nobs;i++)
    {               
      t=Yp[i]+Yn[i];
      resid = (Yp[i] - Ypp[i]);
      if (Ep[i] > 0 && Ep[i] < 1)
	{
	  resid = resid / sqrt(t * Ep[i] * (1 - Ep[i]));
	}
      else
	{
	  if (resid > 0.0 && (Ep[i] == 0.0 || Ep[i] == 1.0))
	    resid = Max_double;
	}
      chs += resid * resid;
      x=Xi[i]*scale;          /*back to original scale.*/

#ifdef MISC_OUT
      printf("%10.4f %11.4f %11.0f %13.3f %9.3f %12.4g\n", 
	     x, Ep[i], Ypp[i], Yp[i], t, resid); 
#endif

      fprintf(fp_out,
#ifndef RBMDS
	      "%10.4f %10.4f %13.3f %9.3f %11.3f %#12.3f\n"
#else
	      "%10.4f %10.4f %30.22g %9.3f %11.3f %30.22g\n"
#endif
	      , x,    Ep[i], Ypp[i], Yp[i], t,   resid); 
    }
  
  df = Nobs - Model_DF(bounded);

  if (df > 0)
    pv = CHISQ(chs, df);   /*p-value for chi-square test*/ 
  else
    pv = 2;

  if (pv != 2)
    {

#ifdef MISC_OUT
      printf("\n Chi-square = %10.2f     DF = %d        P-value = %6.4f\n",
	     chs,df,pv);
#endif

      fprintf(fp_out,
#ifndef RBMDS
	      "\n Chi^2 = %-5.2f     d.f. = %d        P-value = %6.4f\n"
#else
	      "\n Chi^2 = %30.22g     d.f. = %d        P-value = %30.22g\n"
#endif
	      , chs,df,pv);
    }
  else
    {

#ifdef MISC_OUT
      printf("\n Chi-square = %10.2f     DF = %d        P-value =     NA\n",
	     chs,df);
#endif

      fprintf(fp_out,
#ifndef RBMDS
	      "\n Chi^2 = %-5.2f     d.f. = %d        P-value =     NA\n"
#else
	      "\n Chi^2 = %30.22g     d.f. = %d        P-value =     NA\n"
#endif
	      , chs,df);
    }
  FREE_DVECTOR(Ep,1,Nobs);
  FREE_DVECTOR(Ypp,1,Nobs);

}	/*end of Goodness*/

