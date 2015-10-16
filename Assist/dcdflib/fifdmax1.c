#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cdflib.h"

/************************************************************************
FIFDMAX1:
returns the maximum of two numbers a and b
************************************************************************/
double fifdmax1(double a,double b)
/* a     -      first number */
/* b     -      second number */
{
  if (a < b) return b;
  else return a;
}
