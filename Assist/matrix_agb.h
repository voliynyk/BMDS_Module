/****************
*  matrix_agb.h
*  Nov. 22 1996
*****************/


void  SWAP_IVECTOR();
void  FILL_SPECVECTOR();
int   COUNT_SPECVECTOR(); 
void  GETUNKNOWNPARMS();
void  GETKNOWNPARMS();
void  FILL_LDMATRIXBYUD();
void  INITIALIZE_DVECTOR();
void  INITIALIZE_DMATRIX();
void  INITIALIZE_UDMATRIX();
void  MATMPYV ();
void  MATMPYM ();
void MATMPYV2 (int n, int k, double **m,double v[],double *c);
void TRANSPOSE(double **A, double **B, int m, int n);
int INVMAT(double **mat, int n);
void MATMPYM2 (double **A,double **B,double **C,int n,int m, int s);
