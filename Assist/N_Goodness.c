#include <stdlib.h>
#include <math.h>
#include "benchmark.h"
#include "specialfun.h"
#include "allo_memo.h"
#include "in_outfun.h"
#include "gsl/gsl_randist.h"

extern int Model_DF (int []);
extern void Predict(double [], double [], int, double [], double []);


//#define BSDEBUG 1000		/* Bootstrapping debug code (N>0 outputs every N BSDebug Iterations) */

/***********************************************************
*N_Goodness -- Used to test the Goodness of Fit
*
************************************************************/
void
N_Goodness (int ngrp, int nparm, double Parms[], int bounded[],
	    int Nobs, double Xi[], double Yp[], double Yn[], double Ls[],
	    int Xg[], double SR[] )
{
  double *GYpp;			/*predicted positive depdendent values */
  double *GEp;			/*estimated probability */
  double *P;
  double *Gvar, *GYsum, *Gcount;
  int *GrpSize;			/*numbers of obs. in each group. */
  int i, j, k, l;
  double gfit;
    

  /* Note: have to redefine those var. because the data must be */
  /* changed in order to compute goodness-of-fit statistic. */
  GYpp = DVECTOR (1, Nobs);
  GEp = DVECTOR (1, Nobs);
  GYsum = DVECTOR(1, Nobs);
  Gcount = DVECTOR(1, Nobs);
  Gvar = DVECTOR(1, Nobs);
  GrpSize = IVECTOR (1, ngrp);
  P = DVECTOR (1, nparm);



/**** compute the GrpSize[] -- # of obs. in each group***/
  for (i = 1; i <= ngrp; i++)
    GrpSize[i] = 0;

  for (j = 1; j <= Nobs; j++)
    GrpSize[Xg[j]] += 1;


  for (j = 6; j <= nparm; j++)
    P[j - 5] = Parms[j]; 

/** compute the estimated prob. and expected obs.  *******/

  Predict(Xi, Ls, Nobs, Parms, GEp);
  SortByLs (Nobs, ngrp, GrpSize, Ls, Yp, Yn, GYpp, GEp);

  

  for (i = 1; i <= Nobs; i++)
    {
      GYsum[i] = Yp[i] + Yn[i];
      GYpp[i] = GEp[i] * GYsum[i];
      Gvar[i] = GYsum[i] * GEp[i] * (1 - GEp[i]) * (1.0 + (GYsum[i] - 1.0) * P[Xg[i]]);
    }


  /*** compute the goodness-of-fit test for litter data. *********/
  gfit = 0.0;

  for (k = 1; k <= Nobs; k++)
    {
    if (Gvar[k] > 0)
      {
      SR[k] = Gvar[k] > 0 ? (Yp[k] - GYpp[k])/sqrt(Gvar[k]) : 0;
      gfit += SR[k]*SR[k];
      }

    }


 /*** declare random array size for debug output */


  OUTPUT_TEXT ("\n\n                               Litter Data");
  OUTPUT_TEXT ("\n\n           Lit.-Spec.              Litter                          Scaled");
  OUTPUT_TEXT (    "   Dose       Cov.     Est._Prob.   Size    Expected   Observed   Residual");
  OUTPUT_TEXT (    "--------------------------------------------------------------------------");
/*                           11111111112222222222333333333344444444445555555555666666666677777777778 */
/*                  12345678901234567890123456789012345678901234567890123456789012345678901234567890 */
/*                  XXXXXXXXX XXXXXXXXX     XXXXXX     XXXXX     XXXXXXX     XXXXX   XXXXXXXXX */

/*  GXi[0] = 1000000.0;*/
  for (i = 1; i <= Nobs; i++)
    {
      if (i > 1 && Xi[i] > Xi[i - 1])
	{
	  fprintf (fp_out, "\n");
	}
      fprintf (fp_out,
#ifndef RBMDS
	       "%9.4f %9.4f     %6.3f     %5.0f     %7.3f     %5.0f   %9.4f\n"
#else
	       "%9.4f %9.4f     %30.22g %5.0f     %30.22g     %5.0f   %9.4f\n"
#endif
	       , Xi[i], Ls[i], GEp[i], GYsum[i], GYpp[i], Yp[i],
	       (Gvar[i] > 0 ? (Yp[i] - GYpp[i])/sqrt(Gvar[i]) : 0));
    }




  FREE_DVECTOR (GYpp, 1, Nobs);
  FREE_DVECTOR (GEp, 1, Nobs);
  FREE_DVECTOR (GYsum, 1, Nobs);
  FREE_DVECTOR (Gcount, 1, Nobs);
  FREE_DVECTOR (Gvar, 1, Nobs);
  FREE_IVECTOR (GrpSize, 1, ngrp);
  FREE_IVECTOR (P, 1, nparm);
}
