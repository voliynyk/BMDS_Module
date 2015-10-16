#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"




/**********************************************************/
void path_name(char *input_name)
{
  int i, file_flag;

  /* Allow spaces in path name */
  file_flag = 0;

  for (i=0;i<=255;i++)
    {
      if (input_name[i] == '.')
	file_flag = 1;

      if ((input_name[i] == '\0') && (file_flag == 0))
	input_name[i] = ' ';
          
    } /* end for */

} /* end path_name */


