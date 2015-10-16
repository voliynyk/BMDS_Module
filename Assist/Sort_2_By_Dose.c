#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"




/***************************************************************
* Sort_2_By_Dose() sorts an array xxi in ascending order, while
* simultaneously moving the elements of yyi of the same index as 
* the xxi elements into the same indices as the sorted xxi.  A 
* Heap Sort is used to perform this task.  Added to Press pg. 337
****************************************************************/

void 
  Sort_2_By_Dose (int size, double xxi[], double yyi[])
{
  int i, sze, jj, k;
  double temp1, temp2;

  if (size < 2)
    return;

  k = (size >> 1) + 1;
  sze = size;

  for (;;)
    {
      if (k > 1)
	{
	  temp1 = xxi[--k];
	  temp2 = yyi[k];
	}
      else
	{
	  temp1 = xxi[sze];
	  temp2 = yyi[sze];
	  xxi[sze] = xxi[1];
	  yyi[sze] = yyi[1];

	  if (--sze == 1)
	    {
	      xxi[1] = temp1;
	      yyi[1] = temp2;
	      break;
	    }

	}

      i = k;
      jj = k + k;

      while (jj <= sze)
	{
	  if (jj < sze && xxi[jj] < xxi[jj + 1])
	    jj++;

	  if (temp1 < xxi[jj])
	    {
	      xxi[i] = xxi[jj];
	      yyi[i] = yyi[jj];
	      i = jj;
	      jj <<= 1;
	    }
	  else
	    jj = sze + 1;
	}
      xxi[i] = temp1;
      yyi[i] = temp2;
    }
}
