#include "benchmark.h"
#include "ERRORPRT.h"
#include "specialfun.h"

/********************************************************
 **  DTMS3ANOVAC--compute element values for ANOVA table
 *                of a continuous model.  anasum.TEST[i]
 *				 will be equal to -1 if the ith p-value
 *				 is undefined.
 *
 * Modified by: GLN
 * Date: 03/15/2005
 * Purpose: Correct the degree of freedom computation
 *
 * Modified By: Micheal Ferree
 * Date: 25MAY05
 * Purpose: Change to output all tests regardless of type
 *          value.  Fix DF if rho or alpha is specified.
 *          Reorganized to make easier to read.
 *
 *********************************************************/
void DTMS3ANOVAC (int nparm, int Nobs, int Spec[], double A3, double xlk,
		  double A2, double A1, double R, AnaList anasum[], int type,
		  int bounded[])
{

  double pv, dev;
  int parm_known = 0, i;

  /* MJF 25MAY05, changed nparm to 2 because we only care about
     rho and alpha being specified. */
  for (i = 1; i <= 2; i++)
    if (Spec[i] == 1) parm_known++;

  /*********************************************/
  /* Setup model DF and Log Likelihood values. */
  /*********************************************/
  anasum[1].DF = Nobs + 1;	/* DF for model A1.               */
  anasum[1].SS = A1;		    /* Likelihood value for model A1. */

  anasum[2].DF = 2 * Nobs;	/* DF for model A2.               */
  anasum[2].SS = A2;          /* Likelihood value for model A2. */

  /* GLN - 03/18/2005 */
  anasum[3].DF = Nobs + 2 - parm_known;	/* DF for model A3.               */
  anasum[3].SS = A3;          /* Likelihood value for model A3. */

  anasum[4].DF = 2;           /* DF for Reduced model.               */
  anasum[4].SS = R;		    /* Likelihood value for Reduced model. */

  anasum[5].DF = nparm;       /* DF for Fitted model.               */
  for (i = 1; i <= nparm; i++) {
    if (bounded[i]==1 || Spec[i]==1) {
      anasum[5].DF = anasum[5].DF - 1;
    }
  } /* could also just have added up the values in bounded[]. RWS 10/19/2005 */
  anasum[5].SS = xlk;         /* Likelihood value for Fitted model. */

  /***************************************/
  /* Setup values for Test 2 (A1 vs A2). */
  /***************************************/
  dev = 2 * (A2 - A1);		/* Test 2 MSE: Use this to test whether the */
  /* variances are homogeneous (A1 vs. A2)    */

  if (dev < 0.0){
    Warning(" Warning: Likelihood for model A1 larger than the Likelihood for model A2.\n");
    pv = -1;
  }
  else {
    if (anasum[2].DF-anasum[1].DF <= 0) {
      Warning ("Degrees of freedom for Test A1 vs A2 <= 0\n");
      pv = -1;
    }
    else
      /* MJF 25MAY05, changed to use difference for DF. */
      pv = CHISQ (dev, anasum[2].DF-anasum[1].DF /* Nobs - 1 */);	/* P-value for Test 2. */
  }

  anasum[1].MSE = dev;		/* Likelihood ratio for Test 2. */
  anasum[1].TEST = pv;        /* P-value for Test 2.          */

  /***************************************/
  /* Setup values for Test 3 (A2 vs A3). */
  /***************************************/
  dev = 2 * (A2 - A3);	/* Test 3 MSE: Use this to test whether or not
			   the model adequately describes the variances
			   (A2 vs. A3) */
  if (dev < 0.0) {
    Warning(" Warning: Likelihood for model A3 larger than the Likelihood for model A2.\n");
    pv = -1;
  }
  else {
    if (anasum[2].DF-anasum[3].DF <= 0) {
      Warning ("Degrees of freedom for Test A2 vs A3 <= 0\n");
      pv = -1;
    }
    else
      pv = CHISQ (dev, anasum[2].DF-anasum[3].DF /* Nobs - 2 */); /* P-value for Test 3. */
  }

  anasum[2].MSE = dev;	/* Likelihood ratio for Test 3. */
  anasum[2].TEST = pv;	/* P-value for Test 3. */

  /*******************************************/
  /* Setup values for Test 4 (A3 vs fitted). */
  /*******************************************/
  dev = 2 * (A3 - xlk);	/* Test 4 MSE: Use this to test whether the
			   mean model fits (A3 vs. fitted model) */
  if (dev < 0.0) {
    Warning(" Warning: Likelihood for fitted model larger than the Likelihood for model A3.\n");
    pv = -1;
  }
  else {
    if (anasum[3].DF-anasum[5].DF <= 0) {
      Warning ("Degrees of freedom for Test A3 vs fitted <= 0\n");
      pv = -1;
    }
    else
      pv = CHISQ (dev, anasum[3].DF-anasum[5].DF /* Nobs + 2 - nparm + parm_known */); /* P-value for Test 4. */
  }

  anasum[3].MSE = dev;	/* Likelihood ratio for Test 4. */
  anasum[3].TEST = pv;	/* P-value for Test 4. */

  /*******************************************/
  /* Setup values for Test 1 (A2 vs reduced). */
  /*******************************************/
  dev = 2 * (A2 - R);	  /* Test 1 MSE: Use this to test whether there
			     really is a relationship between dose and
			     response (A2 vs reduced model) */
  if (dev < 0.0) {
    Warning(" Warning: Likelihood for model R larger than the Likelihood for model A2.\n");
    pv = -1;
  }
  else {
    if (anasum[2].DF - anasum[4].DF <= 0) {
      Warning ("Degrees of freedom for Test A2 vs R <= 0\n");
      pv = -1;
    }
    else
      pv = CHISQ (dev, anasum[2].DF - anasum[4].DF /* 2 * Nobs - 2 */);
  }

  anasum[4].MSE = dev;	/* Likelihood ratio for Test 4. */
  anasum[4].TEST = pv;	/* P-value for Test 4. */


}

