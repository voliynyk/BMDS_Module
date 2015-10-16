// Nov. 22 1996

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "benchmark.h"

/******************************************
**  ERRORPRT--used to handle error message.
*******************************************/  
void ERRORPRT (char error_text[])
{
	fprintf(stderr,"\n%s\n",error_text);
	fprintf(fp_out," \n%s",error_text);
	
	CLOSE_FILES();
	exit (1);                                
}
/******************************************
**  Warning -- used to handle warning message.
*******************************************/  
void Warning (char error_text[])
{
	fprintf(fp_out," \n%s",error_text);
}
