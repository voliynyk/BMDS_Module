#include <stdio.h>
#include <math.h>
#include <string.h>
#include "benchmark.h"
#include "allo_memo.h"
#include "ERRORPRT.h"
#include "dcdflib/cdflib.h"




/***********************************************************************************************
*
* Get_Names()
*	
*    This function will create the files filename.out, filename.002, and filename.plt
*    (filenames right now are restricted to be 255 characters long)
*			
*    INPUT:
*	 input_name - the name of the filename.(d) file
*
*    OUTPUT:
*	 fileout_name - will contain filename.out
*	 file002_name - will contain filename.002
*	 plot_file_name - will contain filename.plt
*
************************************************************************************************/

void 
  Get_Names (char *input_name, char *fileout_name, char *file002_name, char *plot_file_name)
{
  int i, start_of_ext;

  /* Find where the .(d) extension starts.  If it isn't found, make sure we know (notfound = -9999) */

  start_of_ext = -9999;

  for (i = 0; i <= 255; i++)
    {
      if (input_name[i] == '.')
	{
	  if (input_name[i + 1] == '(' && (input_name[i + 2] == 'd' || input_name[i + 2] == 'D') && input_name[i + 3] == ')')
	    {
	      start_of_ext = i + 1;
	      break;

	    }			/* end if input_name[i+1] ... */

	}			/* end if input_name[i] */

    }				/* end for */


  /* put the filenames together */

  if (start_of_ext < 0)
    {
      ERRORPRT ("ERROR:  Invalid input file");

    }				/* end if */
  else
    {
      for (i = 0; i < start_of_ext; i++)
	{
	  fileout_name[i] = input_name[i];
	  file002_name[i] = input_name[i];
	  plot_file_name[i] = input_name[i];

	}			/* end for */

      fileout_name[start_of_ext] = 'o';
      fileout_name[start_of_ext + 1] = 'u';
      fileout_name[start_of_ext + 2] = 't';
      fileout_name[start_of_ext + 3] = '\0';
      file002_name[start_of_ext] = '0';
      file002_name[start_of_ext + 1] = '0';
      file002_name[start_of_ext + 2] = '2';
      file002_name[start_of_ext + 3] = '\0';
      plot_file_name[start_of_ext] = 'p';
      plot_file_name[start_of_ext + 1] = 'l';
      plot_file_name[start_of_ext + 2] = 't';
      plot_file_name[start_of_ext + 3] = '\0';

    }				/* end if */

}				/* end Get_Names() */
