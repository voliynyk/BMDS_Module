/****************
*  benchmark.h
*  Nov. 22 1996
*****************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h> 
 
#define  MAXDEG 100
#define  MISSING -9999
#define  NR_END 1    
#define  FREE_ARG void *
#define  Yes 1
#define  No 0
#define  ADDED 1
#define  EXTRA 0
#define  true 1
#define  false 0

#define  Log_zero -22				/* e^(-22) ~ 2.7e-10 */
#define  CloseToZero 0.0000001		/* 1.0e-7 */
#define  Pi  3.1415926535897932
#define  OneUponSqrt2Pi .39894228040143267794
#define  TwoPi 6.283185307179587
#define  LnSqrt2Pi -0.9189385332046727417803296 
#define  SQRT2 1.414213562373095049
#define  SQRTPI 1.772453850905516027

/***********************************************
 ** Define indexes for arrays iv and v for dmngb
 ***********************************************/

#define MXITER 17
#define RFCTOL 31
#define XCTOL 32
#define LMAX0 34
#define SCTOL 36

/********************************************
**  Define Data Structures used in BIOSUBCC.C
*********************************************/
/*data type for a variable summary*/
typedef struct var_struct {
	int  flag;
	double S;
	double SS;  
	double CSS;
	double M;
} VarList;
/* VarList varsum; */

/*data type for ANOVA summary*/              
typedef struct ana_struct {
	double SS;
	double MSE;
	int    DF;
	double TEST;
} AnaList;
/* AnaList anasum; */

/*data type for three convergence criteria*/
typedef struct convg_struct {
    double MinPc;      /*parameter criterion*/
    double MinOc;      /*objective criterion*/
    double MinGc;      /*gradient criterion*/
	double Pc;         /*parameter convergence variable*/
	double Oc;         /*objective convergence variable*/
	double Gc;         /*gradient convergence variable*/
} ConvgList; 

/*data type for benchmark dose parameters*/
typedef struct bmd_struct {
	double effect;
	int    risk;
	double level;
	int    *Xg;     /*dose group ID for nested binary model*/
} BMDList;
/* BMDList bmdparm; */
        
/*data type for model IDs in benchmark dose computation*/
enum ModelIDs {
	Probit_linear,            /*Probit_linear = 0*/
	Probit_log,               /*Probit_log = 1*/
	Weibull,                  /*Weibull = 2*/
	Logist_linear,            /*Logist_linear = 3*/
	Logist_log,               /*Logist_log = 4*/
	Gammhit,                  /*Gammhit = 5*/
	Multist,                  /*Multist = 6*/
	Nestlog,                  /*Nestlog = 7*/
	Raivanr,                  /*Raivanr = 8*/
	Nestnctr                  /*Nestnctr = 9*/
};	


/*Define global variables*/
FILE *fp_log;          /*  log file if requested */
FILE         *fp_in, *fp_out, *fp_out2;  /*file pointers*/  
ConvgList    CVG;  /*convergence  criterion structure*/
BMDList      bmdparm; /* parameters for benchmark dose */
int          ITMAX;  /*Maximum allowed  number of iterations*/
#ifndef ASSIST
double Min_increment = DBL_EPSILON, Max_double = DBL_MAX, Rel_Conv, Parm_Conv;
#else
extern double Min_increment, Max_double, Rel_Conv, Parm_Conv;
#endif

/*GLN [10/27/2009]-standardize the file name, model name, user notes, and column name length */
#define  FLENGTH 256	/*File name length */
#define  MNLENGTH 81	/*Model name length */
#define  UNLENGTH 125	/*User note's length */
#define  CNLENGTH 51	/*Column name length */
/*End of GLN [10/27/2009] */

/* Some useful macros */
#define VERYCLOSE(x,y) (fabs((x) - (y)) < Min_increment)



 /*function prototypes*/
void CLOSE_FILES (void);

int run_dmngb(int nparm, double start[], double lower[], double upper[],
	      double maxloglik,
	      double relconv, double parmconv, int itmax, int trymax,
	      void (*func)(), void (*grad)(),
	      long int *uiparm, double *urparm, void (*ufparm)(),
	      int debug, double *fval);

double zeroin(double, double, double,
	      double (*f)(int, double [], double, double),
	      int, double [], double);
void Quantal_CI (int, double [], double [], double, double [], double [],
		 double []);
void Quantal_Goodness (int, int [], double [], int, double [], double [],
		       double [], double);
void Nested_CI(int, int, double [], double [], int [],
	       double, double [], double [], double []);
void N_Goodness (int, int, double [], int [],
		 int, double [], double [], double [], double [],
		 int [], double []);
void N_Bootstrap (int, int, double [], int [],
		 int, double [], double [], double [], double [],
		 int [], double [], int, long);
void SRoI (int, int, double [], double [], int [], double [], double [], double, double);
void PrintData(double [], double [], double [], double [],
	       double [], int [], int);
void Goodness(int nparm, int nparm_known, double Parms[], int type,
	      AnaList anasum[]);
void compute_continuous_liks(int Nobs, int *Ni, double *Ym, double *Yd,
			     double *lk1, double *lk2, double *lkR);
#ifdef WIN32
void BMDS_Win_fpset(void);
#endif
