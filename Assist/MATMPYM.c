#include <math.h>
#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"

#define LITTLE 1.0e-20;



/******************************************************************************
**  MATMPYM--subfunction used to multiply two matrices m and v, and store in c.
*******************************************************************************/
void MATMPYM (double **A,double **B,double **C,int n,int m)
{
	int i,j,k;           
	
	for (i=1;i<=n;i++)
	 {
	  for (j=1;j<=m;j++)
       {
        C[i][j] = 0.0;
        for (k=1;k<=n;k++)
          C[i][j] += A[i][k]*B[k][j];
       }
     }
} 
