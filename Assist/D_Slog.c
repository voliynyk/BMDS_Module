#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"




/**************************************************************** */
/* D_Slog -- derivative of Slog wrt its argument */
/**************************************************************** */

static double _Slogcoefs[4] = {6.7165863851209542e+50,
			       -2.0154759155362862e+35,
			       2.0169759155362859e+19,
			       -710};

double D_Slog(double x)
{
  int i;
  double v;
  if (x >= 1e-16) v = 1/x;
  else
    {
      v = 0.0;
      for (i = 0; i < 3; i++)
	v = (3 - i) * _Slogcoefs[i] + x * v;
    }
  return v;
}
