#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cdflib.h"

/************************************************************************
FIFDMIN1:
returns the minimum of two numbers a and b
************************************************************************/
double fifdmin1(double a,double b)
/* a     -     first number */
/* b     -     second number */
{
  if (a < b) return a;
  else return b;
}
