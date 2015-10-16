#include <math.h>
#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"

#define LITTLE 1.0e-20;



/*********************************************************************
** INITIALIZE_DVECTOR--subfunction used to initialize a double vector.
**********************************************************************/    
void  INITIALIZE_DVECTOR(double v[], int n)
{               
	int i;   
	
	for (i=1;i<=n;i++)
	  v[i]=0;
}
