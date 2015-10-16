#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"




/******************************************************************
* This function will sort a double array in ascending order while
* repositioning the elements of 3 integer arrays accordingly.
* Press, Pg. 337
*******************************************************************/
void 
  Sort_4_By_Dose (int size, double dos[], double tot[], double res[], double litter[])
{
  int i, sze, jj, k;
  double temp1, temp2, temp3, temp4;

  if (size < 2)
    return;

  k = (size >> 1) + 1;
  sze = size;

  for (;;)
    {
      if (k > 1)
	{
	  temp1 = dos[--k];
	  temp2 = tot[k];
	  temp3 = res[k];
	  temp4 = litter[k];
	}
      else
	{
	  temp1 = dos[sze];
	  temp2 = tot[sze];
	  temp3 = res[sze];
	  temp4 = litter[sze];
	  dos[sze] = dos[1];
	  tot[sze] = tot[1];
	  res[sze] = res[1];
	  litter[sze] = litter[1];

	  if (--sze == 1)
	    {
	      dos[1] = temp1;
	      tot[1] = temp2;
	      res[1] = temp3;
	      litter[1] = temp4;
	      break;
	    }

	}

      i = k;
      jj = k + k;

      while (jj <= sze)
	{
	  if (jj < sze && dos[jj] < dos[jj + 1])
	    jj++;

	  if (temp1 < dos[jj])
	    {
	      dos[i] = dos[jj];
	      tot[i] = tot[jj];
	      res[i] = res[jj];
	      litter[i] = litter[jj];
	      i = jj;
	      jj <<= 1;
	    }
	  else
	    jj = sze + 1;
	}
      dos[i] = temp1;
      tot[i] = temp2;
      res[i] = temp3;
      litter[i] = temp4;
    }
}
