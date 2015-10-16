#include "benchmark.h"


/******************************************************
**  COMPUTEDF--compute degree of freedom for the model.
*******************************************************/
int 
  COMPUTEDF (int nparm, int Spec[])
{
  int i, df;

  df = nparm - 1;
  for (i = 1; i <= nparm; i++)
    df = df - Spec[i];

  return df;
}


