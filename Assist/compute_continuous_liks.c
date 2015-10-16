/* Given vectors of data, sample sizes, and standard deviations,
   calculate the likelihoods for models A1, A2, and R for continuous
   data. 
*/
#include <math.h>

void compute_continuous_liks(int Nobs, int *Ni, double *Ym, double *Yd,
			     double *lk1, double *lk2, double *lkR) {
  int Ntot, i;
  double sigma2, ybar;

  /* total sample size: used in all the models */
  Ntot = Ni[1];
  for (i = 2; i <= Nobs; i++) Ntot += Ni[i];

  /* first model: separate mean for each dose group, common estimate
     of the variance
  */
  sigma2 = Yd[1] * (Ni[1] - 1);
  for (i = 2; i <= Nobs; i++)
    sigma2 += Yd[i] * (Ni[i] - 1);
  /* normal MLE for variance divides by N, not N-1 */
  sigma2 = sigma2/Ntot;
  *lk1 = - Ntot*(1.0 + log(sigma2))/2.0;

  /* Model 2: separate mean and variance for each dose group.  The MLE
     variance is (N - 1)/N times the sample variance.
  */

  *lk2 = - Ni[1]*log(Yd[1]*(Ni[1] - 1)/Ni[1])/2.0 - Ntot/2.0;
  for (i = 2; i <= Nobs; i++) 
    *lk2 -= Ni[i]*log(Yd[i]*(Ni[i] - 1)/Ni[i])/2.0;

  /* Model R: common mean and variance. */
  ybar = Ym[1] * Ni[1];
  for (i = 2; i <= Nobs; i++) ybar += Ym[i] * Ni[i];
  ybar = ybar / Ntot;
  sigma2 = Yd[1] * (Ni[1] - 1) + Ni[1]*(Ym[1] - ybar)*(Ym[1] - ybar);
  for (i = 2; i <= Nobs; i++)
    sigma2 += Yd[i] * (Ni[i] - 1) + Ni[i]*(Ym[i] - ybar)*(Ym[i] - ybar);
  sigma2 = sigma2 / Ntot;
  *lkR = -Ntot * (1.0 + log(sigma2))/2.0;
}
