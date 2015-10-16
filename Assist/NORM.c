#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"



/***************************************************
** NORM--subfunction used to compute normal density.
****************************************************/
double 
  NORM (double x)
{
  double y;

  y = OneUponSqrt2Pi * exp (-1 * x * x * 0.5);
  return y;
}
