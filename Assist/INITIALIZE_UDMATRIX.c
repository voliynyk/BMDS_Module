#include <math.h>
#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"

#define LITTLE 1.0e-20;



/****************************************************************
** INITIALIZE_UDMATRIX--subfunction used to initialize updiagonal 
*                       elements of a double matrix. 
*****************************************************************/    
void INITIALIZE_UDMATRIX(double **m, int r, int c)
{       
	int i,j;

	for (i=1;i<=r;i++)
	 {
	  for (j=i;j<=c;j++)
	   m[i][j]=0;
	 }  
}
