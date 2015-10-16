#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "specialfun.h"

/************************************************************************
* {QH 2004/01/14 PR# }
* Added to show version number if executed from command line with -v
* show_version can be modified to take more options. Right now it 
* simply prints out the version number if the parameter is not .(d) file
************************************************************************/

void show_version(char argv[], char version[])
{
  int i, file_flag;

  
  file_flag = 0;

  for (i=0;i<=255;i++)
    {
      if (argv[i] == '.' && argv[i+1] == '(' && argv[i+2] == 'd' && argv[i+3] == ')' && argv[i+4] == NULL)
	  {
	     file_flag = 1;
		 break;
	  }
	  else if (argv[i] == NULL)
			  break;

	  if (i == 255)
	  {
		  fprintf(stderr, "\n\n***File name too long!\n");
		  exit(1);
	  }
    } /* end for */

  if (file_flag == 0)
  {
	fprintf(stderr, "\n%s\n", version);
	exit(0);
  }

} /* end show_version */


