#include <math.h>
#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"

#define LITTLE 1.0e-20;



/********************************************************************
**  FILL_LDMATRIXBYUD--subfunction used to fill low-diagonal elements
*                      by up-diagonal elements.
*********************************************************************/ 
void FILL_LDMATRIXBYUD(double **m, int r, int c)
{
	int i,j;
	
	for (i=1;i<=r;i++)
	 {
	   for (j=i+1;j<=c;j++)
	     m[j][i]=m[i][j];
	 }
}	 
