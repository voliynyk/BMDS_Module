#include <math.h>
#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"

#define LITTLE 1.0e-20;



/*****************************************************************
*  INVMAT will replace the matrix mat with its inverse, n 
*  gives the number of rows/columns of mat
*  Uses lapack and blas routines
*****************************************************************/
void dgesv_(int *N, int *NRHS, double *A, int *LDA,
	    int *IPIV,
	    double *B, int *LDB, int *INFO);

int INVMAT(double **mat, int n)
{
  int N, NRHS, *IPIV, LDB, INFO, LDA;
  double *A, *B;
  int i, j, k;

  LDA = LDB = N = NRHS = n;
  A = (double *) malloc((size_t) n * n * sizeof(double));
  if (A == (double *) NULL) ERRORPRT("Out of memory\n");
  B = (double *) malloc((size_t) n * n * sizeof(double));
  if (B == (double *) NULL) ERRORPRT("Out of memory\n");
  IPIV = (int *) malloc((size_t) n * sizeof(long int));
  if (IPIV == (int *) NULL) ERRORPRT("Out of memory\n");
  k = 0;
  for (j = 1; j <= n; j++)
    for (i = 1; i <= n; i++)
      {
	A[k] = mat[i][j];
	if (i == j)
	  {
	    B[k] = 1.0;
	  }
	else
	  {
	    B[k] = 0.0;
	  }
	k++;
      }
  dgesv_(&N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO);
  k = 0;
  for (j = 1; j <= n; j++)
    for (i = 1; i <= n; i++)
      {
	mat[i][j] = B[k];
	k++;
      }
  free(A);
  free(B);
  free(IPIV);
  return ((int) INFO);
}
