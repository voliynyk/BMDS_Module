#include <math.h>
#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"

#define LITTLE 1.0e-20;



/**********************************************************
**COUNT_SPECVECTOR--- count the number of non-zero  elemants
**                    in the Specvector
***********************************************************/
int COUNT_SPECVECTOR(int nparm, int Spec[])
{
  int    nparm_known;         /* # of specified parameters */
  int    i;
                                
    nparm_known = 0;
    for (i=1; i<=nparm; i++)
	  {
		if(Spec[i] > 0) nparm_known += 1;
      }
		return nparm_known ;
}
