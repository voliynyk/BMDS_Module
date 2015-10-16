#include <math.h>
#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"

#define LITTLE 1.0e-20;



/***************************************
**  GETKNOWNPARMS--get known parameters.
****************************************/ 
void GETKNOWNPARMS(int nparm, int Spec[], double **m, double Parms[])
{
	int   i;
	
	/*assign known parameters*/
	for (i=1;i<=nparm;i++)
	 {
	  if (Spec[i]==1)   
	    m[i][1]=Parms[i];
	 }       
}

