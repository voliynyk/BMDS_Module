#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"




/*****************************************************************
**  PROBABILITY_INRANGE--control P[response] value to be in (0,1).
******************************************************************/
void 
  PROBABILITY_INRANGE (double *ex)
{
  if (*ex <= 0.00000001)
    *ex = 1.0e-7;
  if (*ex >= 0.9999999)
    *ex = 0.9999999;
}
