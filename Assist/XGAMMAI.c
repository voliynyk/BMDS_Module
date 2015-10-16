#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"



/**********************************************
**  XGAMMAI--inverse incomplete Gamma function.
***********************************************/
double 
  XGAMMAI (double p, double shape)
{
  double q, scale, bound, x;
  int which, status;

  which = 2;
  /*bound = x = 0.0;*/ /* returnvals; to avoid error in valgrind */
  scale = 1.0;
  q = 1.0 - p;
  /*status = 0;*/ /* return val; to avoid error in valgrind */
  cdfgam(&which, &p, &q, &x, &shape, &scale, &status, &bound);
  /*  if (status != 0)
    fprintf(stderr,
	   "XGAMMAI:  status from cdfgam is %d. p is %14.4g and shape is %14.4g\n", 
	   status, p, shape);
  */
  return x;
}
