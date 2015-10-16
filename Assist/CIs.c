/* CIs.c -- compute confidence intervals for data */

#include <stdlib.h>
#include <math.h>
#include "specialfun.h"

/* ************************************************************************/
/* Quantal_CI:  compute confidence intervals (as in Fleiss) for quantal */
/*              data expressed as number positive and number negative */
/* */
/*              All input and output vectors are based at 1, that is, the */
/*              elements go from 1..Nobs */
/*    Input: */
/*       Nobs   number of data items */
/*       Yp     vector of number of affected (starts with 1) */
/*       Yn     vector of number NOT affected (total is Yp + Yn) */
/*       conf   confidence required (e.g., if 95% confidence intervals, conf */
/*              is 0.95 */
/*    Output: */
/*       LL     vector if lower confidence limits */
/*       estp   estimated propotions */
/*       UL     vector of upper confidence limits */

void Quantal_CI (int Nobs, double Yp[], double Yn[], double conf, double LL[],
		 double estp[], double UL[])
{
  double cc, cccc, n, phat;
  int i;

  cc = RNORM(1.0 - (1.0 - conf)/2.0);
  cccc = cc * cc;

  for (i = 1; i <= Nobs; i++)
    {
      n = Yp[i] + Yn[i];
      estp[i] = phat = Yp[i]/n;
      LL[i] = (2 * Yp[i] + cccc - 1.0) -
	cc*sqrt(cccc - (2.0 + 1.0/n) + 4.0 * phat * (Yn[i] + 1.0));
      LL[i] = LL[i] / (2.0*(n + cccc));
      UL[i] = (2 * Yp[i] + cccc + 1.0) +
	cc*sqrt(cccc + (2.0 - 1.0/n) + 4.0 * phat * (Yn[i] - 1.0));
      UL[i] = UL[i] / (2.0*(n + cccc));
    }
}
  
