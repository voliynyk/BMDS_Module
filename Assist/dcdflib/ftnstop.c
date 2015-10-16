#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cdflib.h"

/************************************************************************
FTNSTOP:
Prints msg to standard error and then exits
************************************************************************/
void ftnstop(char* msg)
/* msg - error message */
{
  if (msg != NULL) fprintf(stderr,"%s\n",msg);
  exit(EXIT_FAILURE); /* EXIT_FAILURE from stdlib.h, or use an int */
}
