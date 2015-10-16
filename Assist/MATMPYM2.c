#include <math.h>
#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"

#define LITTLE 1.0e-20;



/******************************************************************************
**  MATMPYM2--subfunction used to multiply two matrices A and B, and store in c.
**  n is the number of rows in A, m is the number of columns in A and number of 
**  rows in B, and s is the number of columns in B
*******************************************************************************/
void MATMPYM2 (double **A,double **B,double **C,int n,int m, int s)
{
	int i,j,k;           
	
	for (i=1;i<=n;i++)
	 {
	  for (j=1;j<=s;j++)
       {
        C[i][j] = 0.0;
        for (k=1;k<=m;k++)
          C[i][j] += A[i][k]*B[k][j];
       }
     }
} 
