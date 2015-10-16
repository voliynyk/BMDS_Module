//
//allo_memo.c  
//Nov. 22 1996
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "benchmark.h"
#include "ERRORPRT.h"

/******************************************************************************
**  VLVECTOR--allocate a VarList vector with subscript range [n1,n2].
*******************************************************************************/ 
VarList *VLVECTOR(int n1, int n2)
{
	VarList *v;
	
	v=(VarList *) malloc((size_t) ((n2-n1+1+NR_END)*sizeof(VarList))); 
	 if (!v) ERRORPRT ("Memory allocation failed in VLVECTOR()");
	return v-n1+NR_END; 
}                                         
                                          

/******************************************************************************
**  ALVECTOR--allocate a AnaList vector with subscript range [n1,n2].
*******************************************************************************/ 
AnaList *ALVECTOR(int n1, int n2)
{
	AnaList *v;
	
	v=(AnaList *) malloc((size_t) ((n2-n1+1+NR_END)*sizeof(AnaList))); 
	 if (!v) ERRORPRT ("Memory allocation failed in ALVECTOR()");
	return v-n1+NR_END; 
}                                           


/******************************************************************************
**  IVECTOR--allocate an int vector with subscript range [n1,n2].
*******************************************************************************/ 
int *IVECTOR(int n1, int n2)
{
	int *v;
	
	v=(int *) malloc((size_t) ((n2-n1+1+NR_END)*sizeof(int))); 
	 if (!v) ERRORPRT ("Memory allocation failed in IVECTOR()");
	return v-n1+NR_END; 
}                         
             
                         
/******************************************************************************
**  DVECTOR--allocate a double vector with subscript range [n1,n2].
*******************************************************************************/ 
double *DVECTOR(int n1, int n2)
{
	double *v;
	v=(double *) malloc((size_t) ((n2-n1+1+NR_END)*sizeof(double))); 
	 if (!v) ERRORPRT ("Memory allocation failed in DVECTOR()");
	return v-n1+NR_END; 
}


/******************************************************************************
**  IMATRIX--allocate an int matrix with subscript range [r1,r2][c1,c2].
*******************************************************************************/ 
int **IMATRIX(int r1, int r2, int c1, int c2)
{                   
    int i,nrow,ncol;
	int **m;
	                                  
    nrow = r2-r1+1;
    ncol = c2-c1+1;	   
                          
    /*allocate pointers to rows*/ 
	m=(int **) malloc((size_t) ((nrow+NR_END)*sizeof(int*))); 
	 if (!m) ERRORPRT ("Memory allocation failed in IMATRIX()");
	m += NR_END;                                            
	m -= r1;
	
	/*allocate rows and set pointers to them*/
	m[r1]=(int *) malloc((size_t) ((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[r1]) ERRORPRT ("Memory allocation failed in IMATRIX()");
	m[r1] += NR_END;
	m[r1] -= c1;
	
	for (i=r1+1;i<=r2;i++)
	  m[i]=m[i-1]+ncol;
	  
	/*return pointer to array of pointers to rows*/  
	return m;
}


/******************************************************************************
**  DMATRIX--allocate a double matrix with subscript range [r1,r2][c1,c2].
*******************************************************************************/ 
double **DMATRIX(int r1, int r2, int c1, int c2)
{                   
    int i,nrow,ncol;
	double **m;
	   
    nrow = r2-r1+1;
    ncol = c2-c1+1;	                         
    
   	/*allocate pointers to rows*/     
	m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*))); 
	 if (!m) ERRORPRT ("Memory allocation failed in DMATRIX()");
	m += NR_END;                                            
	m -= r1;
	
	/*allocate rows and set pointers to them*/
	m[r1]=(double *) malloc((size_t) ((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[r1]) ERRORPRT ("Memory allocation failed in DMATRIX()");
	m[r1] += NR_END;
	m[r1] -= c1;
	
	for (i=r1+1;i<=r2;i++)
	  m[i]=m[i-1]+ncol;
	  
	/*return pointer to array of pointers to rows*/  
	return m;
}


/*********************************************************************
**  FREE_VLVECTOR--free a VarList vector with subscript range [n1,n2].
**********************************************************************/ 
void FREE_VLVECTOR (VarList *v, int n1, int n2)
{
	free ((FREE_ARG) (v+n1-NR_END));
}                                             


/*********************************************************************
**  FREE_ALVECTOR--free a AnaList vector with subscript range [n1,n2].
**********************************************************************/ 
void FREE_ALVECTOR (AnaList *v, int n1, int n2)
{
	free ((FREE_ARG) (v+n1-NR_END));
}

 
/******************************************************************************
**  FREE_IVECTOR--free an int vector with subscript range [n1,n2].
*******************************************************************************/ 
void  FREE_IVECTOR (int *v, int n1, int n2)
{
	free ((FREE_ARG) (v+n1-NR_END));
} 
                       

/******************************************************************************
**  FREE_DVECTOR--free a double vector with subscript range [n1,n2].
*******************************************************************************/ 
void FREE_DVECTOR (double *v, int n1, int n2)
{
	free ((FREE_ARG) (v+n1-NR_END));
} 
                              

/******************************************************************************
**  FREE_IMATRIX--free an int matrix with subscript range [r1,r2][c1,c2].
*******************************************************************************/ 
void FREE_IMATRIX (int **m, int r1, int r2, int c1, int c2)
{
	free ((FREE_ARG) (m[r1]+c1-NR_END));
	free ((FREE_ARG) (m+r1-NR_END));
} 


/******************************************************************************
**  FREE_DMATRIX--free a double matrix with subscript range [r1,r2][c1,c2].
*******************************************************************************/ 
void FREE_DMATRIX (double **m, int r1, int r2, int c1, int c2)
{
	free ((FREE_ARG) (m[r1]+c1-NR_END));
	free ((FREE_ARG) (m+r1-NR_END));
}

/******************************************************************************
**  LIVECTOR--allocate a long int vector with subscript range [n1,n2].
*******************************************************************************/
long int *
LIVECTOR (int n1, int n2)
{
  long int *v;

  v = (long int *) malloc ((size_t) ((n2 - n1 + 1 + NR_END) * sizeof (long int)));
  if (!v)
    ERRORPRT ("Memory allocation failed in LIVECTOR()");
  return v - n1 + NR_END;
}


/******************************************************************************
**  FREE_IVECTOR--free an int vector with subscript range [n1,n2].
*******************************************************************************/
void
FREE_LIVECTOR (long int *v, int n1, int n2)
{
  free ((FREE_ARG) (v + n1 - NR_END));
}
