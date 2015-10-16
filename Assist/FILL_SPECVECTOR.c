#include <math.h>
#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"

#define LITTLE 1.0e-20;



/**************************************************************
** FILL_SPECVECTOR--fill an int vector based on another vector.
***************************************************************/  
void FILL_SPECVECTOR(int nparm,double Parms[], int Spec[])
{
	int i;
	double dfabs;
	/*fill in Spec[i]=1 if Parms[i] is unknown */ 
	for (i=1;i<=nparm;i++){
		/* if (Parms[i] != MISSING) */
		dfabs = fabs(Parms[i] - MISSING);
		if( dfabs <= 1e-10 )
			Spec[i] = 0;
		else
			Spec[i] = 1;
	}
}
