#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"



/********************************************************************
** CNORM--subfunction used to compute cumulative normal distribution.

*********************************************************************/
double 
  CNORM (double x)
{
  int which, status;
  double p, q, mean, sd, bound;
  which = 1;
  mean = 0.0;
  sd = 1.0;
  cdfnor(&which, &p, &q, &x, &mean, &sd, &status, &bound);
  return(p);
}
