#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"



/****************************************************** */
/* XGAMMAI_A -- return shape parameter for incomplete */
/*              Gamma function given p and x */
/****************************************************** */
double XGAMMAI_A(double p, double x)
{
  double q, scale, bound, shape;
  int which, status;

  which = 3;
  scale = 1.0;
  q = 1.0 - p;
  cdfgam(&which, &p, &q, &x, &shape, &scale, &status, &bound);
  return shape;
}
