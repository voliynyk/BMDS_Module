#include <stdlib.h>
#include "benchmark.h"
#include "allo_memo.h"

/* NESTED_CI.c -- compute confidence intervals for nested data using */
/*                the Rao-Scott transformation */

/******************************************************************* */
/* raoscott -- computes the rao-scott transformed summary for litter */
/*             data from multiple dose groups. */
/*    inputs: */
/*      ngrp    number of groups */
/*      nobs    total number of observations */
/*      xij     1-based vector of number of affected littermates */
/*              each element corresponds to a single litter */
/*      nij     1-based vector of unaffected littermates; */
/*              each element corresponds to a single litter */
/*      Xg      group indicator */
/*    outputs: */
/*      num     1-based vector containing the numerator of result */
/*      den     1-based vector containing the denominator of result */
/******************************************************************* */

void raoscott(int ngrp, int nobs, double xij[], double nij[],
	      int Xg[], double num[], double den[])
{
  double phat, rij, vi, di, sumnum, sumden, sumsq;
  int i, j, grpsize;

  for (j = 1; j <= ngrp; j++)
    {
      sumnum = sumden = 0.0;
      for (i = 1; i <= nobs; i++)
	{
	  if (j == Xg[i])
	    {
	      sumnum += xij[i];
	      sumden += (xij[i] + nij[i]);
	    }
	}
      phat = sumnum/sumden;
      sumsq = 0;
      grpsize = 0;
      for (i = 1; i <= nobs; i++)
	{
	  if (j == Xg[i])
	    {
	      grpsize++;
	      rij = xij[i] - phat * (xij[i] + nij[i]);
	      sumsq += rij * rij;
	    }
	}
      vi = grpsize * sumsq/((grpsize - 1) * sumden * sumden);
      if (phat > 0 && phat < 1)
	{
	  di = sumden * vi/(phat * (1.0 - phat));
	}
      else
	{
	  di = 1.0;
	}
      num[j] = sumnum/di;
      den[j] = sumden/di;
    }

}
  
/******************************************************************* */
/* NESTED_CI  -- use the rao-scott transformation to compute confidence */
/*               limits for nested (litter) quantal data */
/*   inputs: */
/*      ngrp  number of dose groups */
/*      nobs  total number of observations */
/*      Yp    the ith element contains the number of affected pups in the */
/*            ith litter */
/*      Yn    the ith element contains the number of non-affected pups */
/*            in the ith litter */
/*      Xg    indicator that shows which dose group the ith litter */
/*            belongs to */
/*      size  size of desired confidence limits (ie, result will be */
/*            approximate size*100% confidence limits) */
/*   outputs: */
/*      LL     vector of lower confidence limits */
/*      phat   estimate of proportion affected for each dose group */
/*      UL     vector of upper confidence limits */
/******************************************************************* */


void Nested_CI(int ngrp, int nobs, double Yp[], double Yn[], int Xg[],
	       double size, double LL[], double phat[], double UL[])
{
  double *num, *den;
  int i;

  num = DVECTOR(1, ngrp);
  den = DVECTOR(1, ngrp);

  raoscott(ngrp,nobs,Yp, Yn, Xg, num, den);
  /* Quantal_CI requires number pos and number neg */
  /* after this, den contains the number unaffected (transformed) */
  for (i = 1; i <= ngrp; i++) den[i] -= num[i];
  Quantal_CI (ngrp, num, den, size, LL, phat, UL);

  FREE_DVECTOR(num, 1, ngrp);
  FREE_DVECTOR(den, 1, ngrp);
}

  
