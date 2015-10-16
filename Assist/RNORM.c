#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"




/********************************************************
**  RNORM--subfunction used to inverse cumulative normal.
*********************************************************/
double 
  RNORM (double p)
{
  double q, x, mean, sd, bound;
  int which, status;
  which = 2;
  mean = 0.0;
  sd = 1.0;
  q = 1.0 - p;
  cdfnor(&which,&p, &q, &x, &mean, &sd, &status, &bound);
  return x;
}
