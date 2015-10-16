#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cdflib.h"

/************************************************************************
FIFIDINT:
Truncates a double precision number to a long integer
************************************************************************/
long fifidint(double a)
/* a - number to be truncated */
{
  return (long)(a);
}
