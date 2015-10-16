#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"




/****************************************************************** */
/* QCHISQ -- inverse qui-square function */
/****************************************************************** */

double QCHISQ(double p, int m)
{
  int which, status;
  double chisq, q, df, bound;

  which = 2;
  df = (double) m;
  q = 1.0 - p;
  cdfchi(&which, &p, &q, &chisq, &df, &status, &bound);

  if (status == 0) return chisq;
  else 
    {
      Warning("Error in computing chi-square; returning -1");
      return -1;
    }
}
