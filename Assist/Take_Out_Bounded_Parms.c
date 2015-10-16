#include "benchmark.h"
#include "allo_memo.h"

int 
  Take_Out_Bounded_Parms (int nparm, int *bounded, double **vcv_matrix)
/***********************************************************************
*  Take_Out_Bounded_Parms() will take a square matrix and remove the
*  row and column associated with any parameter that is at a bound.
*
*	INPUT:	nparm - the number of rows/columns on vcv_matrix
*			bounded - bounded[i] = 0 if the parameter is not on a bound
*					  bounded[i] = 1 if the parameter is on a bound
*
*	INPUT/OUTPUT:	vcv_matrix is the matrix that is to be adjusted.
*					Upon output, vcv_matrix will be a square matrix
*					of smaller or equal size as it was upon input, with
*					0's occupying the last num_deleted rows and columns
*					of the matrix, where num_deleted is the number of
*					rows/columns deleted
*
*	RETURN VALUE = the number of rows in the new, adjusted vcv_matrix
*
**********************************************************************/
{
  int i, j, num_bounded, counti, countj;
  double **temp_vcv;

  num_bounded = 0;

  for (i = 1; i <= nparm; i++)
    {
      if (bounded[i] == 1)
	{
	  num_bounded++;	/* Finds the number of bounded
				   bounded parameters */
	}			/* end if */
    }				/* end for */

  /* Create a temporary matrix to hold the new matrix */

  temp_vcv = DMATRIX (1, nparm - num_bounded, 1, nparm - num_bounded);

  /* Take out the appropriate rows/columns */

  counti = 0;

  for (i = 1; i <= nparm; i++)
    {
      if (bounded[i] == 0)
	{
	  counti++;		/* So far, this element stays */
	  countj = 0;

	  /* Now check the column */

	  for (j = 1; j <= nparm; j++)
	    {
	      if (bounded[j] == 0)
		{
		  /* This element stays */

		  countj++;
		  temp_vcv[counti][countj] = vcv_matrix[i][j];

		}		/* end if(bounded[j]) */

	    }			/* end for j */

	}			/* end if(bounded[i]) */

    }				/* end for i */

  /* Put it back into vcv_matrix */

  for (i = 1; i <= nparm - num_bounded; i++)
    {
      for (j = 1; j <= nparm - num_bounded; j++)
	{
	  vcv_matrix[i][j] = temp_vcv[i][j];

	}			/* end for j */

    }				/* end for i */

  /*Fill the rest of it up with 0's */

  for (i = nparm - num_bounded + 1; i <= nparm; i++)
    {
      for (j = 1; j <= nparm; j++)
	{
	  vcv_matrix[i][j] = 0.0;
	  vcv_matrix[j][i] = 0.0;

	}			/* end for j */

    }				/* end for i */

  FREE_DMATRIX(temp_vcv, 1, nparm-num_bounded,1, nparm - num_bounded);
  return nparm - num_bounded;

}				/* end Take_Out_Bounded_Parms() */

