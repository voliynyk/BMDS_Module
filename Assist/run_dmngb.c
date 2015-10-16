/***************************************************************
 * run_dmngb -- sets up auxiliary arrays to run dmngb, and tries
 * repeatedly (by tweaking the initial parameter values) to find
 * an optimum.
 * input:
 *       nparm: the number of parameters to be estimated
 *       start: array of initial values
 *       lower: array of lower bounds
 *       upper: array of upper bounds
 *       Dscale: array of scales for parameters
 *       maxloglik: -log-likelihood at full model (for absolute function
 *                  convergence)
 *       relconv: relative function convergence limit
 *       parmconv: parameter convergence limit
 *       itmax: maximum number of iterations allowed dmngb
 *       trymax: maximum number of retries allowed when dmngb does not
 *               converge
 *       func: function to minimize
 *       grad: gradient of func
 *       uiparm: array of user-specified long int values for func and grad
 *       urparm: array of user-specified double values for func and grad
 *       ufparm: user specified function argument for func and grad
 *       debug: if 0, no debugging output; if 1, debugging output to unit 16
 * output:
 *       start: final parameter vector returned by dmngb
 *       fval: value of the function at the optimum
 * return value:
 *      -2: unable to allocate required memory
 *      -1: problem too big (MAXPARMS exceeded)
 *       0: for any of the values returned by dmngb that indicate convergence
 *          (including singular convergence): iv[0] in 1..7.
 *       other: the value returned in iv[0] on return from dmngb the last try.
 ******************************************************************/
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

#define MAXPARMS 30

#define fmin(x,y) ((x) < (y) ? (x) : (y))
#define fmax(x,y) ((x) < (y) ? (y) : (x))

extern int  
dmngb_(long int *nparm, double *d, double *x, 
       double b[][2], 
       void (*Logist_lk)(long int *nvar, double *x, long int *nf, double *f,
			 long int *uiparm, double *urparm, void (*ufparm)()),
       void (*Logist_g)(long int *nvar, double *x, long int *nf, double *g,
			long int *uiparm, double *urparm, void (*ufparm)()),  
       long int *iv, long int *liv, 
       long int *lv, double *v, long int *uiparm, double *urparm,
       void (*ufparm)()); 
extern void divset_ (long int *arg, long int *iv, long int *liv,
                     long int *lv, double *v);


int run_dmngb(int nparm, double start[], double lower[], double upper[],
	      double maxloglik,
	      double relconv, double parmconv, int itmax, int trymax,
	      void (*func)(), void (*grad)(),
	      long int *uiparm, double *urparm, void (*ufparm)(),
	      int debug, double *fval)
{
  long int nvar, *iv, liv, lv, divsetarg;
  double *v, bounds[MAXPARMS][2], *x, Dscale[MAXPARMS];
  int i, j;
  long int retval;

  nvar = nparm;
  /* make sure we're not going to overrun our bounds matrix */
  if (nvar > MAXPARMS) return(-1);
  /* set up the bounds matrix */
  for (i = 0; i < nvar; i++)
    {
      bounds[i][0] = lower[i];
      bounds[i][1] = upper[i];
    }
  x = (double *) malloc((size_t) nvar * sizeof(double));
  for (i = 0; i < nvar; i++)
    {
      x[i] = start[i];
    }
  liv = 100;
  iv = (long int *) malloc((size_t) liv * sizeof(long int));
  lv = 171 + nvar*(nvar+15)/2;
  v = (double *) malloc((size_t) lv*sizeof(double));
  if (iv == (long int *) NULL || v == (double *) NULL || x == (double *) NULL)
    return -2;
  divsetarg = 2;
  /* make sure we've allocated enough memory */
  j = 0;
  do
    {
      iv[0] = 13;
      divset_(&divsetarg, iv, &liv, &lv, v);
      j++;
      retval = iv[0];
      if (retval != 12 && j < 5)
	{
	  liv = iv[43];
	  lv = iv[44];
	  free(iv);
	  free(v);
	  iv = (long int *) malloc((size_t) liv * sizeof(long int));
	  v = (double *) malloc((size_t) lv*sizeof(double));
	}
    } while (retval != 12 && j < 5);
  if (retval != 12) return -2; /* unable to allocate memory */
  i = 0;
  do
    {
      /* Now set up the problem properly */
      iv[0] = 0;
      divset_(&divsetarg, iv, &liv, &lv, v);
      if (!debug)
	iv[20]=0;
      else
	iv[20]=16;
      v[30] = maxloglik;
      v[31] = relconv;
      v[32] = parmconv;
      iv[17] = itmax;
  /* scale by reciprocal of starting values for DMNGB */
      for (j = 0; j < nvar; j++) 
	{
	  if (fabs(x[j]) > 0) 
	    {
	      Dscale[j] = 1/fabs(x[j]);
	    }
	  else 
	    {
	      Dscale[j] = 1.0;
	      
	    } /* end if */
	  
	} /* end for */


      /*dmngb_ is the Fortran optimization routine*/
      /*x will have the maximum likelihood estimate when finished*/
      dmngb_(&nvar, Dscale, x, bounds, func, grad, iv, &liv, &lv,
	     v, uiparm, urparm, ufparm);
      retval = iv[0];
      /* an experiment: if we don't converge, restart where we */
      /* left off, but rescale things. */
      if (i == 0 && retval > 6)
	{
	  i++;
	  continue;
	}
      i++;
      /* if we did not converge, perturb the initial point and try again */
      if (retval > 6 && i < trymax)
	{
	  double U, U2, span;
	  for (j = 0; j < nvar; j++)
	    {
	      U = rand()/32768.0;
	      U2 = rand()/32768.0;
	      if (U2 < 0.5)
		{
		  if (lower[j] == -DBL_MAX)
		    {
		      x[j] = start[j] * U;
		    }
		  else
		    {
		      span = start[j] - lower[j];
		      x[j] = start[j]  - 0.1*U * span;
		    }
		}
	      else
		{
		  if (upper[j] == DBL_MAX)
		    {
		      x[j] = start[j] * (1.0 + U);
		    }
		  else
		    {
		      span = upper[j] - start[j];
		      x[j] = start[j] + 0.1* U * span;
		    }
		}
	    }
	}
    } while (retval > 6 && i < trymax);
  *fval = v[9];
  free(iv);
  free(v);
  for (i = 0; i < nvar; i++) start[i] = x[i];
  free(x);
  if (retval <= 6)
    return(0);
  else
    return ((int) retval);
}

