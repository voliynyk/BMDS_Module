#include <float.h>
#include <limits.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"

/****************************************************
 **  CLOSE_FILES--used to close input and output files.
 *****************************************************/                                         
void CLOSE_FILES (void)
{
  if (fclose (fp_in) != 0 || fclose (fp_out) != 0 
#ifndef RBMDS
      || fclose (fp_out2) != 0
#endif
      )
    ERRORPRT ("Error in closing opened files.");
}
