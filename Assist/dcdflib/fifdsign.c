#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cdflib.h"

/************************************************************************
FIFDSIGN:
transfers the sign of the variable "sign" to the variable "mag"
************************************************************************/
double fifdsign(double mag,double sign)
/* mag     -     magnitude */
/* sign    -     sign to be transfered */
{
  if (mag < 0) mag = -mag;
  if (sign < 0) mag = -mag;
  return mag;

}
