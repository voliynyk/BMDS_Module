#include <math.h>
#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"

#define LITTLE 1.0e-20;



/***************************************************
** SWAP_IVECTOR--swap two elements in an int vector.
****************************************************/                                                             
void SWAP_IVECTOR (int v[], int a, int b)
{
	int temp;
	
	temp = v[a];
	v[a] = v[b];
	v[b] = temp;
}

