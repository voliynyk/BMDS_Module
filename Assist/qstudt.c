#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"




/**************************************************************** */
/* qstudt -- student-t function */
/**************************************************************** */
double qstudt(double p, int m)
{
    int which, status;
    double q, df, bound, t;
    q = 1.0 - p;
    df = (double) m;
    which = 2;

    cdft (&which, &p, &q, &t, &df, &status, &bound);

    if (status == 0) return t;
    else
      {
	Warning("Error in computing t; returning -1");
	return -1;
      }
}
