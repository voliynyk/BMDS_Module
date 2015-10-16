#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"




/***************************************************************
 * Slog -- Safe log: a "log" function which can be evaluated at
 *         0 and negative numbers.
 *
 *         if x >= 1e-16, Slog(x) = log(x) (natural log)
 *         if x < 1e-16, Slog(x) is a cubic polynomial in x
 *         such that:  Slog(1e-16) = log(1e-16)
 *                     Slog(0)     = -710 (~log(.Machine$double.min))
 *                     dSlog/dx|x=1e-16 = dlog/dx|1e-16
 *                     d^2Slog/dx^2|x=1e-16 = d^2log/dx^2|x=1e-16
 * The idea is to replace log in binomial log-likelihoods with
 * Slog, and stop testing for special values of x
 ****************************************************************/
static double _Slogcoefs[4] = {6.7165863851209542e+50,
			       -2.0154759155362862e+35,
			       2.0169759155362859e+19,
			       -710};

double Slog(double x)
{
  int i;
  double v;
  if (x >= 1e-16) v = log(x);
  else
    {
      v = 0.0;
      for (i = 0; i < 4; i++)
	v = x * v + _Slogcoefs[i];
    }
  return v;
}
		    
