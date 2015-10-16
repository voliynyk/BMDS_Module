#include <math.h>
#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"
#include "matrix_agb.h"

#define LITTLE 1.0e-20;



/**********************************************                                                                                 
** MATINVS--invert a double covariance matrix.                                                            
***********************************************/
int MATINVS (int n, double **m, int v[])
{
	int     *cols;          /*temp int vector*/
	double  *scale;         /*temp double vector*/
	int     i,j,k;          /*iteration control variables*/
    double  pivot;          /*control variable*/ 
    int     col;            /*control variable*/
    double  temp;      
    
    /*allocate memory for vectors*/  
    cols = IVECTOR(1,n);
    scale = DVECTOR (1,n);

	for (i=1;i<=n;i++)
	 {
	  for (j=1;j<=n;j++)
	   {
	    if (i==j)
	     {
	      if (v[i]==1)
	        m[i][j] = 1;      
	     }
	    else if (v[i]==1 || v[j]==1)
	     m[i][j] = 0; 
	   }
	 } 
   
	for (i=1;i<=n;i++)
     {
      cols[i] = i;
      scale[i] = 0.0;
      for (j=1;j<=n;j++)
        scale[i] += fabs(m[i][j]);
     } 
  
	for (i=1;i<=n;i++)
	 {               
	  pivot = m[cols[i]][cols[i]];
	  for (j=i;j<=n;j++)
	   {
	    if (scale[cols[i]]==0.0)
          return 0;
        if (fabs(pivot/scale[cols[i]]) < fabs(m[cols[j]][cols[j]]) / scale[cols[i]]) 
         {
          pivot = m[cols[j]][cols[j]];
          SWAP_IVECTOR(cols,j,i);
         }
        col = cols[i];
       }  
       
      if (pivot == 0)      
        return 0;
        
      m[col][col] = 1;
      for (j=1;j<=n;j++)
        m[col][j] = m[col][j] / pivot;

      for (j=1;j<=n;j++)
       {
        if (j!=col)
         { 
          temp = m[j][col];
          m[j][col] = 0;
          for (k=1;k<=n;k++)
            m[j][k] = m[j][k] - m[col][k] * temp;
         }

       }
    
      if (i < (n-1))
       {
        for (j=i+1;j<=n;j++)
         {
          scale[cols[j]] = 0;
          for (k=i+1;k<=n;k++)
            scale[cols[j]] += fabs(m[cols[j]][cols[k]]);
         }
       } 
	 }  
	 
	for(i=1;i<=n;i++)
	 {
      if (v[i] == 1)
       m[i][i] = 0;
     }   
     
    /*free allocated memory*/  
    FREE_IVECTOR (cols,1,n);
    FREE_DVECTOR (scale,1,n);
    return 1;
	     
} /*end of MATINVS*/

