#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"




/****************************************************************
**  PROBABILITY_W--control positive probability W value in (0,1).
*****************************************************************/
void 
  PROBABILITY_W (double *W)
{
  if (*W <= 0.01)
    *W = 0.01;
  if (*W >= 0.99)
    *W = 0.99;
}
