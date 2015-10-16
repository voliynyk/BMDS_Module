#include <math.h>
#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"

#define LITTLE 1.0e-20;



/******************************************************************************
**  MATMPYV--subfunction used to multiply two matrices m and v, and store in c.
*******************************************************************************/                                
void MATMPYV (int n,double **m,double v[],double *c)
{
	int    i,j;
	double t;
	
	for (i=1;i<=n;i++)
	 {
      t = 0.0;
      for (j=1;j<=n;j++)
        t += m[i][j] * v[j];
      c[i] = t;
     }
} /*end of MATMPYV subfunction*/

