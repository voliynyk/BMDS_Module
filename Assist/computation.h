/******************
*  computation.h
*  Nov. 22 1996
******************/
void INITIALIZE_CVG3(void);
void COMPUTE_DTMS3LOGLKH(int Nobs, double *xlkf,double *xlkr,
                         VarList varsum[],double Yp[],double Yn[]);
void COMPUTE_DTMS3FLOGLKH(double *xlk, double ex, double Yp, double Yn);
void COMPUTE_DTMS3CVG_Gc(int nparm,double **y, double **d);
void COMPUTE_DTMS3CVG_Pc(int nparm,double **bb,double bp[]);
int  COMPUTEDF(int nparm,int Spec[]);
void DTMS3ANOVA (int nparm,int Nobs,int Spec[],double xlkf,double xlk,
                 double xlkr,AnaList anasum[], int []);
void VARSUMCOMP(int n, int Nobs, double v[], VarList varsum[]);
void DTMS3ANOVAC (int nparm,int Nobs,int Spec[],double A3, double xlk, double A2, double A1, double R, AnaList anasum[], int type, int *bounded);
int Take_Out_Bounded_Parms(int nparm, int *bounded, double **vcv_matrix);
int Get_Linear_Trend(int N, double *doses, double *means, int *numi);

