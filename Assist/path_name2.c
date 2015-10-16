#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"




void path_name2(int argc, char** argv, char* outstr)
{
  int i;
  char *space = " ";

  strcpy(outstr, argv[1]);

  if (argc > 3)
    {
      for (i = 2; i < argc-1; i++) {
	strcat(outstr, space);
	strcat(outstr, argv[i]);
      }
    }
  strcat(outstr, space);
  strcat(outstr, argv[argc-1]);
}
