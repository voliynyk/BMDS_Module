#include <math.h>
#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"

#define LITTLE 1.0e-20;



/******************************************************************************
**  MATMPYV2--subfunction used to multiply two matrices m and v, and store in c.
*******************************************************************************/                                
void MATMPYV2 (int n, int k, double **m, double v[],double *c)
{
	int    i,j;
	double t;
	
	for (i=1;i<=n;i++)
	 {
      t = 0.0;
      for (j=1;j<=k;j++)
        t += m[i][j] * v[j];
      c[i] = t;
     }
} /*end of MATMPYV subfunction*/
