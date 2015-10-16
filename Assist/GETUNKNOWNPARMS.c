#include <math.h>
#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"

#define LITTLE 1.0e-20;



/***************************************************************
**  GETUNKNOWNPARMS--get computed values for unknown parameters.
****************************************************************/ 
void GETUNKNOWNPARMS(int nparm, int Spec[], double Parms[], double d[])
{
	int   i;
	
	/*assign computed unknown parameters*/
	for (i=1;i<=nparm;i++)
	 {
	  if (Spec[i]==0)   
	    Parms[i]=d[i];
	 }       
}
