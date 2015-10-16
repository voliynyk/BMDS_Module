#include <math.h>
#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"

#define LITTLE 1.0e-20;



/*****************************************************************
*  Given a matrix A, and the dimensions of A, m being the number
*  of rows, and n the number of columns, TRANSPOSE will place the 
*  transpose of A into a matrix B
*****************************************************************/

void TRANSPOSE(double **A, double **B, int m, int n)
{
	int i, j;

	for(i = 1; i <= m; i++)
		for(j = 1; j <= n; j++)
			B[j][i] = A[i][j];
}
