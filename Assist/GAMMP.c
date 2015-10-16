#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"


/*****************************************************
** GAMMP--return the incomplete gamma function P(a,x).
*         Uses cdflib.
******************************************************/
double 
  GAMMP (double shape, double x)
{
  double p,q, scale, bound;
  int which, status;

  which = 1;
  scale = 1.0;
  /* bound = 0.0; */
  /*status = 0; */
  /*p = q = 0.0;*/
  cdfgam(&which, &p, &q, &x, &shape, &scale, &status, &bound);
  if (status != 0)
    fprintf(stderr,
	   "GAMMP:  status from cdfgam is %d. x is %14.4g and shape is %14.4g\n", 
	   status, x, shape);
  return p;
}

