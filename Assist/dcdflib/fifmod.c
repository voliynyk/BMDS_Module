#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cdflib.h"

/************************************************************************
FIFMOD:
returns the modulo of a and b
************************************************************************/
long fifmod(long a,long b)
/* a - numerator */
/* b - denominator */
{
  return a % b;
}
