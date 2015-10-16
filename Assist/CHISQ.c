 #include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"




/******************************************************
**  CHISQ--compute cumulative Chi-square for the model.
*******************************************************/
double 
  CHISQ (double x, int m)
{
  int which, status;
  double p, q, df, bound;

  which = 1;
  df = (double) m;
  cdfchi(&which, &p, &q, &x, &df, &status, &bound);

  if (status == 0) return q;
  else 
    {
      Warning("Error in computing chi-square; returning 2");
      return 2.0;
    }
      
}
