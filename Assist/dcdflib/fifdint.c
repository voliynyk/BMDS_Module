#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cdflib.h"

/************************************************************************
FIFDINT:
Truncates a double precision number to an integer and returns the
value in a double.
************************************************************************/
double fifdint(double a)
/* a     -     number to be truncated */
{
  long temp;
  temp = (long)(a);
  return (double)(temp);
}
