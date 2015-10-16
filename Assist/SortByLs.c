#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"




/*************************************************************
*SortByLs --- uesd to sort the data by litter size in each dose
*             group (also rearrange Yn, Yn and Ypp), for 
*             computing goodnees-of-fit statistic.
*
**************************************************************/
void 
SortByLs (int nobs, int ngrp, int GrpSize[], double GLs[],
	  double GYp[], double GYn[], double GYpp[],
	  double GEp[])
     /* GrpSize[1]+ .... + Grpsize[ngrp] = nobs . */
{
  int i, j, k;
  double a1, a2, a3, a4, a5;
  int *Cgs;			/* cumulation GrpSize; */

  Cgs = IVECTOR (1, ngrp);
  for (i = 1; i <= ngrp; i++)	/* counting the cumulating GS. */

    {
      Cgs[i] = 0;
      for (j = 1; j <= i; j++)
	Cgs[i] += GrpSize[j];
    }

  Cgs[0] = 0;
  for (k = 1; k <= ngrp; k++)
    {
      for (j = Cgs[k - 1] + 2; j <= Cgs[k]; j++)
	{
	  a1 = GLs[j];
	  a2 = GYp[j];
	  a3 = GYn[j];
	  a4 = GYpp[j];
	  a5 = GEp[j];
	  i = j - 1;
	  while (i > Cgs[k - 1] && GLs[i] > a1)
	    {
	      GLs[i + 1] = GLs[i];
	      GYp[i + 1] = GYp[i];
	      GYn[i + 1] = GYn[i];
	      GYpp[i + 1] = GYpp[i];
	      GEp[i + 1] = GEp[i];
	      i--;
	    }
	  GLs[i + 1] = a1;
	  GYp[i + 1] = a2;
	  GYn[i + 1] = a3;
	  GYpp[i + 1] = a4;
	  GEp[i + 1] = a5;
	}
    }
  FREE_IVECTOR (Cgs, 1, ngrp);
}
