#include <math.h>
#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"

#define LITTLE 1.0e-20;



/*********************************************************************
** INITIALIZE_DMATRIX--subfunction used to initialize a double matrix. 
**********************************************************************/    
void INITIALIZE_DMATRIX(double **m, int r, int c)
{       
	int i,j;

	for (i=1;i<=r;i++)
	 {
	  for (j=1;j<=c;j++)
	   m[i][j]=0;
	 }  
}                                          
