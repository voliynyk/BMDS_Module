/****************************************************************
*
* IMPORTANT NOTE:  The following variable is the version number for
*                  the current model.  THIS MUST BE CHANGED as
*				   important changes are made to the models.
*
*****************************************************************/
char Version_no[]="Exponential Model. (Version: 1.10;  Date: 01/12/2015)";
/*
char Version_no[]="Exponential Model. (Version: 1.8;  Date: 11/12/2014)";
char Version_no[]="Exponential Model. (Version: 1.7;  Date: 12/10/2009)";
char Version_no[]="Exponential Model. (Version: 1.61;  Date: 7/24/2009)";
char Version_no[]="Exponential Model. (Version: 1.6;  Date: 7/20/2009)";
char Version_no[]="Exponential Model. (Version: 1.5;  Date: 4/23/2009)";
char Version_no[]="Exponential Model. (Version: 1.42;  Date: 12/10/2008)";
char Version_no[]="Exponential Model. (Version: 1.32;  Date: 11/24/2008)";
char Version_no[]="Exponential Model. (Version: 1.31;  Date: 11/19/2008)";
char Version_no[]="Exponential Model. (Version: 1.2;  Date: 10/23/2008)";
char Version_no[]="Exponential Model. (Version: 1.1;  Date: 08/15/2007)";
char Version_no[]="Exponential Model. (Version: 1.0;  Date: 3/30/2007)";
char Version_no[]="Exponential Model. (Version: 0.9;  Date: 2/20/2006)";
*/
/****************************************************************
* Exponential.C - a ANSI C program for Exponential model fitting.
*
* Date: Feb 28, 2008
* 
*****************************************************************/
/********************************************************************
* Modification Log:
*
* Version Number: 1.0
* Modified By: Geoffrey Nonato
* Date: 3/31/2007
* Reason:	1. Changed most variables to a naming convention (i.e, prefix "g" for global scope, 
*			   followed with data type, then variable name, example: gpiSpecVector, for global 
*			   pointer to integer SpecVector).
*			2. Reformat/reorganize source code
*			3. Added clean-up code (free-up allocated memories).
*			4. Added printing routine as specifed (to conform with current BMDS OUT file format).
*			5. Deleted variables that are declared/allocated memory, but never use.
*			6. Restructure all exponential variables by grouping them in one structure.
*			7. Added functions and structures.
*			8. Fixed AThree_Fit.
*			9. Fixed Exponential_Fit.
*
* Version Number: 1.1
* Modified By: Geoffrey Nonato
* Date: 8/15/2007
* Reason:	1. Added BMDL computation.
*
* Version Number: 1.2
* Modified By: Geoffrey Nonato
* Date: 10/23/2008
* Reason: Modify code to use donlp2 fortran instead of dmngb fortran.
*
* Version Number: 1.31
* Modified By: Geoffrey Nonato
* Date: 11/19/2008
* Reason: Modify code that writes to "out" file.
*
* Version Number: 1.32
* Modified By: Geoffrey Nonato
* Date: 11/19/2008
* Reason: Added code to print .002 file.
*
* Version Number: 1.42
* Modified By: Geoffrey Nonato
* Date: 12/10/2008
* Reason:	Added A3 computation when constant variance is "0".
*
* Version Number: 1.5
* Modified By: R. Woodrow Setzer
* Date: 4/23/2009
* Reason:
*   Relative BMR was giving BMD == BMDL.  Reformatted FORTRAN codes
*   to remove all tab characters and keep all statements within standard
*   F77 columns: 7 - 72.  The resulting code now gives BMDL < BMD
*
* Version Number: 1.6
* Modified By: Geoffrey Nonato
* Date: 7/20/2009
* Reason:
*   Fix Observation # < parameter # for Exponential model, to allow 
*   continuation of processing with models who are valid
*   
*
* Version Number: 1.61
* Modified By: Geoffrey Nonato
* Date: 7/24/2009
* Reason:
*   Fix DoInitParam() and DoParamEstimates() to display the right values
*   
* Version Number: 1.7
* Modified By: G. Nonato
* Modification Date: 12/10/2009
* Reason:
*      To be able to process files/folders with spaces (PR 257)
*      Fix program freeze due to long variable names (PR 278)
*      Process long files up to 256 characters (PR 303 and 308)
*      Modify code for easy maintenance, all lengths for file name,
*        model name, and column names are placed in benchmark.h
*
* Version Number: 1.8
* Modified By: Cody Simmons
* Modification Date: 11/12/2014
* Reason:
*	Added standard error calculations and output.  Updated 
*	DoParamEstimates() to simplify formatting. 
*
* Version Number: 1.10
* Modified By: Cody Simmons
* Modification Date: 01/12/2015
* Reason:
*     Version number updated for consistency with BMDS website
*********************************************************************/

#include <float.h>
#include <limits.h>
#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>

#include <benchmark.h>
#include <ERRORPRT.h>
#include <allo_memo.h>
#include <matrix_agb.h>
#include <specialfun.h>
#include <computation.h>
#include <in_outfun.h>

extern void getmle_(long int *ndoses, double doses[], double means[],
					long int nanimals[], double svar[], long int *nparm,
					double parms[], long int fixed[], double fixedval[],
					long int *restrict,long int *adverse, double parms2[],double *ll,
					long int *optite, long int *nresm, long int bind[],
					long int *model_type, long int *flag, long int *lognormal);

extern void getcl_(long int *which, long int *ndoses, double doses[], double means[],
				   long int nanimals[], double svar[], long int *nparm, double *bmr,
				   double *bmd, double *target, double parms[],
				   long int fixed[], double fixedval[], long int *risktype,
				   long int *restrict, double *bmdl, double parms2[],
				   long int *optite, long int *nresm, long int bind[], long int *adv,
				   long int *model_type, long int *flag, long int *lognormal);

extern void getmlea3_(long int *ndoses, double doses[], double means[],
					  long int nanimals[], double svar[], long int *nparm,
					  double parms[], long int fixed[], double fixedval[],
					  long int *restrict, double parms2[],double *ll,
					  long int *optite, long int *nresm, long int bind[]);

void Exponential_fit(int nModel);
void GetMLEParms(double *p, int size, int nModel);
void GetMoreParms(double *p, int size, int nModel);
void GetNewParms(double *p, int size, int nModel);
void AThree_Fit(int nparm, double p[], double gtol, int *iter, int nModel);

#define EPS 3.0e-8
//#define DO_LOG		//Uncomment to generate log file code
#define NBR_OF_PARMS 6
#define NBR_OF_MODELS 4

#define dmax(x,y) ((x) < (y) ? (y) : (x))
#define dmin(x,y) ((x) > (y) ? (y) : (x))

/* control models */
/* it may be replaced by gpiSpecVector[i] later */
/* For DMNGB */
int giErrorFlag; /* Error States from DMNGB */
int giErrorFlag_A3;

int    girestrict;    /* flag for restricting slope>=1 */

char    gacFileOut[FLENGTH];  /*output temp file*/
char    gacFileOut2[FLENGTH];
char	gacPltFileName[FLENGTH];  /* file to pass to GnuPlot */


/*************************************************************
*   giDo_Log = true -> log file is made
*   giDo_Log = false -> no log
*************************************************************/
int     giDo_Log = false;			/*  Set true, to activate switch for log file   */
int     giDo_PrintLog = false;		/*  option to print info to log file   */

char   *gaParm_Name[]={"lnalpha", "rho", "a", "b", "c", "d"};
char   gacFileName2[FLENGTH], gacLogFile[FLENGTH], *gcDot2;
int    giNbrParm_A3;

int	   gpimodtype;
double *gdBMD;
double *gdBMDL;
int *giRun;
int *giBmdRun;

//int    *gpiSpecVector;	/* vector used to identify user input parm.
//						will be 1 if user specified that parameter,
//						0 otherwise */
//int    *gpiIniSpVector;			/* user initial specified parameter vector */
//double *gpdIniParaVector;	    /* user intialized parameter vector */

/* New Matrices for each model*/
int		**gppiSpecPara;	/* int vector used to identify user specified parm vector, will be 1 if user specified that parameter.*/ 
int		**gppiInitPara;	/* int vector used to identify user initialized parm vector.*/ 
double	**gppdSpecPara; /* vector used to identify user specified parm values.*/ 
double	**gppdInitPara; /* vector used to identify user initialized parm values.*/ 
int		*gpiInitialized; /*did the user choose to initialize*/
int		*gpiExp_Known; /* # of known parameters for each model*/
double	**gppdPredict; /* double matrix for prediction values*/

int    *gpiNi;
double *gpdYm;
double *gpdYd;
double *gpdStDev;
double *gpdScXi;		/* new one, added 01/28/08 */
double *gpdXi;		/* independent variable data array */
double *gpdxxi;     /*dose values when not divided to dose groups*/
double *gpdyyi;     /*response values when not divided to dose groups*/
double *gpdYsum;    /*sum of responses within a group*/
double *gpdScYm;    /*scaled mean response data */
double gdmaxYm;     /*max mean response */
double *gpdScYd;    /*scaled sample variances of responses */
double gdmaxYd;     /*max sample variance */

double gdxmax, gdxmin, gdyymin, gdyymax;


double *mg;         /* first derivative of mean wrt parameters */
double **mg2;        /* second derivative of mean wrt parameters */

double **gppdVCV2;      /* variance/covariance (actually hessian) matrix for Exp2 model */
double **gppdVCV3;      /* variance/covariance (actually hessian) matrix for Exp3 model */
double **gppdVCV4;      /* variance/covariance (actually hessian) matrix for Exp4 model */
double **gppdVCV5;      /* variance/covariance (actually hessian) matrix for Exp5 model */
double **gppdVCV[NBR_OF_MODELS];     /* Array of pointers to variance/covariance (actually hessian) matrices for each model */

double **gppdVCV2_adj;  /* variance/covariance matrix for Exp2 model adjusted for bounded/specified parameters */
double **gppdVCV3_adj;  /* variance/covariance matrix for Exp3 model adjusted for bounded/specified parameters */
double **gppdVCV4_adj;  /* variance/covariance matrix for Exp4 model adjusted for bounded/specified parameters */
double **gppdVCV5_adj;  /* variance/covariance matrix for Exp5 model adjusted for bounded/specified parameters */
double **gppdVCV_adj[NBR_OF_MODELS];    /* Array of pointers to adjusted variance/covariance matrices for each model */
      
int *adj_vcv_rows;
int **bounded;       /* 1 if a parameter is estimated at a boundary, 0 otherwise */ 

double **parmSE;     /* contains standard error estimates for each parameter in each model */



int giPrintLL = 0; /* if 1, print a line after the anodev table with the */
/* log-likelihood to high precision.  Triggered by a */
/* negative value for gsExpoVars.iMaxIter*/

/* changing variable */
int giReplace, giBrat;
int gpiIter;


//Miscellaneous functions and structure added by Geoffrey, 03/24/2007

int giNmiss = 0;
void HeaderAndReadVar_2Out(char *dFileName, char* clocktime, int iModelNbr);	//Header and Read variables to out file
typedef struct exponent 
{
	int iSelect;
	char caUser_note[UNLENGTH];
	int iIn_Type;
	int iNbrObs_Total;	/* number of observations when type=1 or sample size when type=0 */
	int iIsInitial;
	int iSign;
	int iLogNormal;
	int iExactSolution;
	int iMaxIter;
	double dRel_Conv;
	double dParm_Conv;
	int iBmdlCurve;		/* flag for bmdl curve option */
	int iBmdose;
	int iAppendIt;		/* flag: 0 for append output, 1 for overwrite */
	int iSmooth;		/* flag: 0 for unique bmdl curve, 1 for C-spline */
	int iBmr_Type;
	double dBmdEffect;	//bmdparm.effect
	int iCons_Var;		/* 1: constant, 0: power of mean */
	int iBmdRisk;			//bmdparm.risk
	double dBmdConfi_Level;	//bmdparm.level
	char caDoseName[CNLENGTH];
	char caNiName[CNLENGTH];
	char caMeanName[CNLENGTH];
	char caStdevName[CNLENGTH];
	char caResponseName[CNLENGTH];
	char caModelName[MNLENGTH];
} exponentialVars;
exponentialVars gsExpoVars;

typedef struct LikeInterest {
	char caModel[8];
	double dLogLikelihood;
	int    iDF;
	double dAIC_PValue;
} LikeInterestList;

//typedef struct BenchMark {
//	int    iModelNbr;
//	double dBMR;
//	double dBMD;
//	double dBMDL;
//} BMarkList;
//
//#define IBMarkList 4		//Nbr of Models in a BenchMark List
//BMarkList *aBMarks;

double **gppdMLEs;
double *gpdLikelihoods;
double *gpdA3s;
double gdAdverseBMDEffect;

typedef enum {
	eM2=1, eM3, eM4, eM5
} eMs;

typedef enum {
	eAlpha=1, eRho, ea, eb, ec, ed
} eParms;

int DoInitParam(int iModelNbr);
int DoParamEstimates(int iModelNbr);
void DoDataEstimateInt(int iModelNbr);
int DoLikehoods(double A1, double A2, double A3, double R, int iModelNbr);
void DoBenchMark(int iGrouped, int iModelNbr);

void DoInitValuesM2(double p[]);
void DoInitValuesM4(double p[]);
double getBMD23(double dBMD, double da, double db, double dd);
double getBMD45(double dBMD, double da, double db, double dc, double dd);
double BMDL_func(int nModel);

//int Get_Linear_Trend2 (int N, double *doses, double *means, int *numi);
double LogLike(int nmodel, double p[]);
int doDot002(int iStart, int iEnd, int iGrouped);

int READ_OBSDATA4V(int iNbr_Obs,double gpdXi[],int gpiNi[],double gpdYm[],double gpdYd[]);
int READ_OBSDATA2V(int iNTotal, double pdxxi[], double pdyyi[]);



void MeanPart (int obs, double *mg, int mod);
void Mean2Part (int obs, double **mg2, int mod);
void VarPart (int obs, double *mg, double *vg, int mod);
void Var2Part (int obs, double *mg, double **mg2, double **vg2, int mod);
void F1iDoublePart( double **Fn1i, int obs, int mod);
void F2iDoublePart (double **Fn2i, int obs, int mod);
void F3iDoublePart (double **Fn3i, int obs, int mod);
void Exp_vcv (double ***gppdVCV, int mod);
void Get_DTMSVCV (int mod);
void Calc_ParmSE (int mod);



void FreeUp_mem()
{
	FREE_DVECTOR (gpdScYm, 1, gsExpoVars.iNbrObs_Total);
	FREE_DVECTOR (gpdScYd, 1, gsExpoVars.iNbrObs_Total);

	if (gsExpoVars.iIn_Type==1) 
	{
		FREE_DVECTOR (gpdYd, 1, gsExpoVars.iNbrObs_Total);
		FREE_DVECTOR (gpdYm, 1, gsExpoVars.iNbrObs_Total);
		FREE_IVECTOR (gpiNi, 1, gsExpoVars.iNbrObs_Total);
		FREE_DVECTOR (gpdXi, 1, gsExpoVars.iNbrObs_Total);
		FREE_DVECTOR (gpdStDev, 1, gsExpoVars.iNbrObs_Total);
	}
	else
	{
		FREE_IVECTOR (gpiNi, 1, gsExpoVars.iNbrObs_Total);
		FREE_DVECTOR (gpdYm, 1, gsExpoVars.iNbrObs_Total);
		FREE_DVECTOR (gpdYd, 1, gsExpoVars.iNbrObs_Total);
		FREE_DVECTOR (gpdXi, 1, gsExpoVars.iNbrObs_Total);
		FREE_DVECTOR (gpdxxi, 1, gsExpoVars.iNbrObs_Total);
		FREE_DVECTOR (gpdyyi, 1, gsExpoVars.iNbrObs_Total);
		FREE_DVECTOR (gpdYsum, 1, gsExpoVars.iNbrObs_Total);
	}

	FREE_IMATRIX (gppiSpecPara, 1, NBR_OF_MODELS, 1, NBR_OF_PARMS);
	FREE_IMATRIX (gppiInitPara, 1, NBR_OF_MODELS, 1, NBR_OF_PARMS);

	FREE_DMATRIX (gppdSpecPara, 1, NBR_OF_MODELS, 1, NBR_OF_PARMS);
	FREE_DMATRIX (gppdInitPara, 1, NBR_OF_MODELS, 1, NBR_OF_PARMS);
	FREE_DMATRIX (gppdMLEs, 1, NBR_OF_MODELS, 1, NBR_OF_PARMS);
	FREE_DVECTOR (gpdLikelihoods, 1, NBR_OF_MODELS);

	FREE_DMATRIX (gppdPredict, 1, NBR_OF_MODELS, 1, gsExpoVars.iNbrObs_Total);

	FREE_IVECTOR (gpiInitialized, 1, NBR_OF_MODELS);
	FREE_IVECTOR (gpiExp_Known, 1, NBR_OF_MODELS);
	FREE_IVECTOR (giRun, 1, NBR_OF_MODELS);
	FREE_IVECTOR (giBmdRun, 1, NBR_OF_MODELS);

	FREE_DVECTOR (gdBMD, 1, NBR_OF_MODELS);
	FREE_DVECTOR (gdBMDL, 1, NBR_OF_MODELS);
	FREE_DVECTOR (gpdA3s, 1, NBR_OF_MODELS);

	if (fp_log != (FILE *) NULL)
		fclose(fp_log);

	return;
}

/****************************************************************
** main--main function used to call Exponential mode fitting program.
Includes: biosubcc.c--common subfunction C program.	int iNtot=0;
*****************************************************************/
int main (int argc, char *argv[])
{
	int		i, j, jj, junk;	     /* iteration variable */
	int     iGroup;
	double  dlikeA1, dlikeA2, dlikeA3, dlikeR;     /* log likelihoods */
	//double  *pdParms;	     /* parameter array */
	//double  *pdInitParms;	 /* Initial values parameter array */

	//AnaList *psAnasum;       /*information for ANONA analysis*/  
	char	junkname[FLENGTH];
	char	caLong_path_name[FLENGTH];
	char	caModelRuns[6]; /* models that are being run */
	char	caModelSpan[4];
	int		iModelSpan;
	int		iStart;
	int		iEnd;

	double dvv, dbmr_root;
	gpiIter = 0;

	int iTemp_sign;
	time_t  ltime;
	time( &ltime );

	if(argc == 2) show_version(argv[1], Version_no);

	if(argc < 2){
		fprintf(stderr, "ERROR:  Requires two arguments\nUsage:  %s <file.(d)>\n", argv[0]);
		fprintf (stderr, "   or:  %s -v for version number.\n", argv[0]);
		exit (1);
	}

	if (argc > 2){
		path_name2(argc, argv, caLong_path_name);
		argv[1] = caLong_path_name;
	}

	/* open the log file if giDo_Log = 1 */
#ifdef DO_LOG
	if (giDo_Log)
	{
		strcpy(gacLogFile,argv[1]);
		gcDot2 = strchr(gacLogFile, (int) '.');
		(*gcDot2) = (char) 0;
		strcpy(gacFileName2,gacLogFile);
		strcat(gacLogFile,"-Exp.log");
		fp_log = fopen(gacLogFile, "w");

		if (fp_log == (FILE *) NULL)
			ERRORPRT("Unable to open log for Expoential.C.");
		fprintf(fp_log,"  (d)File = %s\n", argv[1]);
		fprintf(fp_log,"caLong_path_name = %s\n", caLong_path_name);
		fprintf(fp_log,"  argc = %d\n", argc);
	}
#endif

	fp_in=fopen(argv[1], "r");

	/* Check if input file is open, if not, print error message and exit */
	if (fp_in==NULL){
		fprintf(stderr,"Error in opening input file.\n");
		fprintf (stderr,"...Exited to system!\n");
		exit (1);
	}

	/* begin reading input file from batch file (.(d) ext.) */
	fscanf(fp_in, "%s",gsExpoVars.caModelName );

	/* select is the Exponential model number */
	//fscanf(fp_in, "%d", &gsExpoVars.iSelect);
	gsExpoVars.iSelect = 2;
	fscanf(fp_in, "%[ ^\n]", gsExpoVars.caUser_note);
	fscanf(fp_in, "%[^\n]", gsExpoVars.caUser_note);
	fscanf(fp_in, "%s", junkname);
	fscanf(fp_in, "%s", junkname);

	/* gsExpoVars.iIn_Type=1 if input format is gpdXi, gpiNi, gpdYm, gpdYd. */
	/* gsExpoVars.iIn_Type=0 if input format is giNTotal, gpdXi, Y_ij. */

	//if (gsExpoVars.iIn_Type==1)
	//	fscanf(fp_in, "%d", &gsExpoVars.iNbrObs_Total);
	//else
	fscanf(fp_in, "%d", &gsExpoVars.iIn_Type);

	fscanf(fp_in, "%d", &gsExpoVars.iNbrObs_Total);

	fscanf(fp_in, "%d", &gsExpoVars.iSign);

	fscanf(fp_in, "%s", caModelRuns);

	fscanf(fp_in, "%s", caModelSpan);

	fscanf(fp_in, "%d", &gsExpoVars.iLogNormal); // = 0 for normal; = 1 for lognormal

	fscanf(fp_in, "%d", &gsExpoVars.iExactSolution); // = 1 for exact solution; 0 for approx solution (in lognormal case)

	// always exact solution when normal distribution
	if(gsExpoVars.iLogNormal == 0) gsExpoVars.iExactSolution = 1;

	// never exact solution when lognormal distribution and summary data
	if((gsExpoVars.iLogNormal == 1) && (gsExpoVars.iIn_Type==1)) gsExpoVars.iExactSolution = 0;

	i = sscanf(caModelSpan, "%d", &iModelSpan);
	iStart = iModelSpan/10;
	iEnd = iModelSpan%10;

	iModelSpan = strlen(caModelSpan)>2?1:0;	//grouped

	//if ((gsExpoVars.iSign != 1) || (gsExpoVars.iSign != -1))
	//	gsExpoVars.iSign = 0;

	fscanf(fp_in,"%d%lf%lf%d%d%d%d", &gsExpoVars.iMaxIter, &gsExpoVars.dRel_Conv, &gsExpoVars.dParm_Conv,
		&gsExpoVars.iBmdlCurve, &gsExpoVars.iBmdose, &gsExpoVars.iAppendIt, &gsExpoVars.iSmooth);

	fscanf(fp_in,"%d%lf%d%lf",&gsExpoVars.iBmr_Type,&gsExpoVars.dBmdEffect,&gsExpoVars.iCons_Var,&gsExpoVars.dBmdConfi_Level);

	if(gsExpoVars.iLogNormal == 1) gsExpoVars.iCons_Var = 1;  // always constant variance for lognormal distrib

#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"   \n\niIn_Type = %d\n", gsExpoVars.iIn_Type);
		fprintf(fp_log,"   iNbrObs_Total = %d\n", gsExpoVars.iNbrObs_Total);
		fprintf(fp_log,"   iSign = %d\n", gsExpoVars.iSign);

		fprintf(fp_log,"   iMaxIter = %d\n", gsExpoVars.iMaxIter);
		fprintf(fp_log,"   dRel_Conv = %g\n", gsExpoVars.dRel_Conv);
		fprintf(fp_log,"   dParm_Conv = %g\n", gsExpoVars.dParm_Conv);

		fprintf(fp_log,"   iBmr_Type = %d\n", gsExpoVars.iBmr_Type);
		fprintf(fp_log,"   dBmdEffect = %g\n", gsExpoVars.dBmdEffect);
		fprintf(fp_log,"   iCons_Var = %d\n", gsExpoVars.iCons_Var);
		fprintf(fp_log,"   dBmdConfi_Level = %g\n", gsExpoVars.dBmdConfi_Level);
		fprintf(fp_log,"   iStart = %d\n", iStart);
		fprintf(fp_log,"   iEnd = %d\n", iEnd);
		fprintf(fp_log,"   iModelSpan = %d\n", iModelSpan);
		fprintf(fp_log,"   LogNormal = %d\n", gsExpoVars.iLogNormal);
		fprintf(fp_log,"   ExactSolution = %d\n", gsExpoVars.iExactSolution);
	}
#endif

	gsExpoVars.iBmdRisk = 1;
	junk = 0;    /* used to see if an extension was added to output file name */


	if (gsExpoVars.iMaxIter < 0)
	{
		gsExpoVars.iMaxIter = -1 * gsExpoVars.iMaxIter;
		giPrintLL = 1;
	}
	/* gsExpoVars.dBmdEffect is the benchmark risk level */
	/* gsExpoVars.iBmdRisk is 1 is extra risk, 0 if added risk */
	/* gsExpoVars.dBmdConfi_Level is the confidence level */


	Get_Names(argv[1], gacFileOut, gacFileOut2, gacPltFileName);

	//if(gsExpoVars.iAppendIt==Yes)
	//	fp_out=fopen(gacFileOut,"a");
	//else
	fp_out=fopen(gacFileOut,"w");

#ifndef RBMDS
	fp_out2=fopen(gacFileOut2,"w");
#endif
	if (fp_out==NULL 
#ifndef RBMDS
		||fp_out2==NULL
#endif
		)
	{
		printf("Error in opening  output files.\n");
		printf ("...now exiting to system...\n");

		fprintf(fp_out,"Error in opening output files.\n");
		fprintf (fp_out,"...Exited to system!\n");
#ifdef DO_LOG
		if (giDo_Log) 
		{
			fflush(fp_log);
			fclose(fp_log);
		}
#endif
		CLOSE_FILES();
		exit (1);
	}
	if (gsExpoVars.iBmdose < 0 || gsExpoVars.iBmdose > 1)
		ERRORPRT("Error in choosing benchmark dose computation.");

	/*Allocate memories for the model matrices*/
	gppiSpecPara = IMATRIX (1, NBR_OF_MODELS, 1, NBR_OF_PARMS);
	gppiInitPara = IMATRIX (1, NBR_OF_MODELS, 1, NBR_OF_PARMS);
	gppdSpecPara = DMATRIX (1, NBR_OF_MODELS, 1, NBR_OF_PARMS);
	gppdInitPara = DMATRIX (1, NBR_OF_MODELS, 1, NBR_OF_PARMS);
	gpiInitialized = IVECTOR (1, NBR_OF_MODELS);
	gpiExp_Known = IVECTOR (1, NBR_OF_MODELS);

	gppdPredict = DMATRIX (1, NBR_OF_MODELS, 1, gsExpoVars.iNbrObs_Total);

	gdBMD = DVECTOR (1, NBR_OF_MODELS);
	gdBMDL = DVECTOR (1, NBR_OF_MODELS);
	giRun = IVECTOR (1, NBR_OF_MODELS);
	giBmdRun = IVECTOR(1, NBR_OF_MODELS);
	gpdA3s = DVECTOR(1, NBR_OF_MODELS);


	mg = DVECTOR(1, NBR_OF_PARMS);
        mg2 = DMATRIX(1, NBR_OF_PARMS, 1, NBR_OF_PARMS);
        gppdVCV2 = DMATRIX(1, NBR_OF_PARMS, 1, NBR_OF_PARMS);
        gppdVCV3 = DMATRIX(1, NBR_OF_PARMS, 1, NBR_OF_PARMS);
        gppdVCV4 = DMATRIX(1, NBR_OF_PARMS, 1, NBR_OF_PARMS);
        gppdVCV5 = DMATRIX(1, NBR_OF_PARMS, 1, NBR_OF_PARMS);

        gppdVCV[0] = gppdVCV2;
        gppdVCV[1] = gppdVCV3;
        gppdVCV[2] = gppdVCV4;
        gppdVCV[3] = gppdVCV5;

	bounded = IMATRIX(1, NBR_OF_MODELS, 1, NBR_OF_PARMS);
	adj_vcv_rows = IVECTOR(1, NBR_OF_MODELS);

	parmSE = DMATRIX(1, NBR_OF_MODELS, 1, NBR_OF_PARMS);




	/*Check if memory allocation is a success*/
	if (!gppiSpecPara || !gppiInitPara || !gppdSpecPara || !gppdInitPara || !gpiInitialized || !gpiExp_Known || !gppdPredict) 
	{
#ifdef DO_LOG
		if (giDo_Log) 
		{
			fflush(fp_log);
			fclose(fp_log);
		}
#endif
		ERRORPRT ("Memory allocation failed in parameter matrices");
	}

	/*Initialize MATRICES*/
	for(i=1; i<=NBR_OF_MODELS; i++)
	{
		for(j=1; j<=NBR_OF_PARMS; j++)
		{
			gppiSpecPara[i][j] = 0;
			gppiInitPara[i][j] = 0;
			gppdSpecPara[i][j] = 0.0;
			gppdInitPara[i][j] = 0.0;
			bounded[i][j] = 0;
			parmSE[i][j] = 0.0;
		}
		gpiInitialized[i]=0;
		gpiExp_Known[i] = 0;
		giRun[i] = 1;
		giBmdRun[i] = 1;
		gdBMD[i] = 0.0;
		gdBMDL[i] = 0.0;
	}




//CWS   CV
/*	bounded[1][2] = 1;
	bounded[1][5] = 1;
	bounded[1][6] = 1;
	bounded[2][2] = 1;
	bounded[2][5] = 1;
	bounded[3][2] = 1;
	bounded[3][6] = 1;
	bounded[4][2] = 1;
*/
// NCV
/*	bounded[1][5] = 1;
	bounded[1][6] = 1;
	bounded[2][5] = 1;
	bounded[2][6] = 1;
        bounded[3][6] = 1;
*/	



	/*Read parameters of all the models here,
	then populate specified and initialized vectors*/
	for(i=1; i <= NBR_OF_MODELS; i++)
	{
		READ_PARAMETERS(NBR_OF_PARMS, gppdSpecPara[i]);  //Populate specified para

		FILL_SPECVECTOR(NBR_OF_PARMS, gppdSpecPara[i], gppiSpecPara[i]);

		fscanf(fp_in,"%d", &gpiInitialized[i]);

		READ_PARAMETERS(NBR_OF_PARMS, gppdInitPara[i]);  //Populate initialized para

		FILL_SPECVECTOR(NBR_OF_PARMS, gppdInitPara[i], gppiInitPara[i]);

		if(gsExpoVars.iCons_Var==1)  
		{
			gppdSpecPara[i][(int)eRho]=0.0; //this is the pdParms before, set to 0, no matter what's specified
			gppiSpecPara[i][(int)eRho]=1; //rho specified
			//Bruce, the rest of Instruction 5 should be here, but needs more info
			/* GLN TODO (Done {02/19/08}): {02/18/08}  I think the following 3 lines are what the rest of instruction 5 should be -- uncomment if looks correct */
			gppdInitPara[i][(int)eRho]=0.0; //this is the pdParms before, set to 0, no matter what's specified
			gppiInitPara[i][(int)eRho]=1; //rho initialized
			//AnaList *psAnasum = ALVECTOR(1, 5);
		}

		//Bruce, as per your instruction 3, do the following assignments ALWAYS
		//eAlpha=1, eRho, ea, eb, ec, ed
		switch(i)
		{
		case (int)eM2:
			gppdSpecPara[i][(int)ec] = 0;
			gppdSpecPara[i][(int)ed] = 1;
			gppiSpecPara[i][(int)ec] = 1;
			gppiSpecPara[i][(int)ed] = 1;
			break;
		case (int)eM3:
			gppdSpecPara[i][(int)ec] = 0;
			gppiSpecPara[i][(int)ec] = 1;
			break;
		case (int)eM4:
			gppdSpecPara[i][(int)ed] = 1;
			gppiSpecPara[i][(int)ed] = 1;
			break;
		default:
			break;
		}
		/* GLN TODO (Done 02/19/08): {02/18/08}  move this for loop AFTER the following switch*/
		for(j=((int)ea); j<=NBR_OF_PARMS; ++j)
		{
			if(gppiSpecPara[i][j]==1)
				gpiExp_Known[i]+=1;
		}

	} /*End of reading parameters for all the models*/

	/* Print values to log for investigation */
#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		for(i=1; i<= NBR_OF_MODELS; i++)
		{
			fprintf(fp_log,"\nModel %d Specified Values: ", i+1);
			for(j=1; j <= NBR_OF_PARMS; j++)
			{
				fprintf(fp_log," %g", gppdSpecPara[i][j]);
			}

			fprintf(fp_log,"\nModel %d Specified Vectors: ", i+1);
			for(j=1; j <= NBR_OF_PARMS; j++)
			{
				fprintf(fp_log," %d", gppiSpecPara[i][j]);
			}

			fprintf(fp_log,"\n\nModel %d Initialize?: %d\n", i+1, gpiInitialized[i]);

			fprintf(fp_log,"\nModel %d Initialized Values: ", i+1);
			for(j=1; j <= NBR_OF_PARMS; j++)
			{
				fprintf(fp_log," %g", gppdInitPara[i][j]);
			}
			fprintf(fp_log,"\nModel %d Initialized Vectors: ", i+1);
			for(j=1; j <= NBR_OF_PARMS; j++)
			{
				fprintf(fp_log," %d", gppiInitPara[i][j]);
			}

			fprintf(fp_log,"\n\nModel %d gpiExp_Known: %d\n\n", i+1, gpiExp_Known[i]);
		}
	}
#endif
	for(i = 0; i < 4; i++)
	{
		if(caModelRuns[i]=='1')
			giRun[i+1] = 1;
		else
			giRun[i+1] = 0;
	}

	/*Test if # of observation < # of parameters, then exit (instruction 7) */
	for(i=1; i<=NBR_OF_MODELS; i++)
	{
		if (giRun[i] == 1 && gsExpoVars.iNbrObs_Total < (4 - gpiExp_Known[i]))
		{
			//GLN TODO 2/19/08: have program print message in output where BMD would be that says:
			//"BMR choice is outside the predicted range of mean values; no BMD estimate is possible."
			// GLN TODO initialize and use vector giRun[1,4] starting with all elements = 1 (will run).
			// when giRun[i] = 0 then no further processing on that model.
			// this needs to be done before in main: where it says "instruction 7" where the decision about
			// too many parameters is done model by model: each time there are too many parameters, for 
			// that model (i), set giRun[i]=0 and do not do further processing on that model (e.g., do not
			// run exponential_fit) -- this means setting up if statements in some places conditioning on 
			// giRun[i] being equal to 1.
			//ERRORPRT("Observation # < parameter # for Exponential model.");
			fprintf(stderr,"\n%s %d.\n","Observation # < parameter # for Exponential model", i+1);
			giRun[i] = 0;
		}
	}

	/*obtain observation data into gpdYp, gpdYn, gpdXi, Ls, Xg vectors*/
	if (gsExpoVars.iIn_Type == 1)
		fscanf(fp_in,"%s%s%s%s", gsExpoVars.caDoseName, gsExpoVars.caNiName, gsExpoVars.caMeanName, gsExpoVars.caStdevName);
	else
		fscanf(fp_in,"%s%s", gsExpoVars.caDoseName, gsExpoVars.caResponseName);

	//Allocate memories for sMLEs + likelihood
	gppdMLEs = DMATRIX(1, NBR_OF_MODELS, 1, NBR_OF_PARMS);
	gpdLikelihoods = DVECTOR(1, NBR_OF_MODELS);

	if (gsExpoVars.iIn_Type==1) 
	{
		gpdYm = DVECTOR (1, gsExpoVars.iNbrObs_Total);
		gpdYd = DVECTOR (1, gsExpoVars.iNbrObs_Total);
		gpdXi = DVECTOR (1, gsExpoVars.iNbrObs_Total);
		gpiNi = IVECTOR (1, gsExpoVars.iNbrObs_Total);
		gpdStDev = DVECTOR (1, gsExpoVars.iNbrObs_Total);
		for(i =1; i <= gsExpoVars.iNbrObs_Total; i++)
			gpdStDev[i] = 0.0;

		giNmiss = READ_OBSDATA4V(gsExpoVars.iNbrObs_Total, gpdXi, gpiNi, gpdYm, gpdStDev);

		for(i =1; i <= gsExpoVars.iNbrObs_Total; i++)
			gpdYd[i] = gpdStDev[i]*gpdStDev[i];

		/* extern variable gsExpoVars.iNbrObs_Total has been changed. */
		gsExpoVars.iNbrObs_Total -= giNmiss;

		// working with lognormal and summary data – get approx log-scale sample means and variances
		if(gsExpoVars.iLogNormal == 1)
		{
			for(i =1; i <= gsExpoVars.iNbrObs_Total; i++)
			{
				/* need to add check: if any gpdYm[i] is <=0 
				print error message on screen and in output file: 
				Data set includes negative responses; no negative 
				values are allowed when a lognormal distribution is assumed.”  
				Then stop processing */
				if(gpdYm[i] <= 0) //Exit
				{
					fprintf(fp_out, "\n\n\nData set includes negative responses:\nno negative values are allowed when a lognormal distribution is assumed."); 
					FreeUp_mem();
					ERRORPRT("Data set includes negative responses:\nno negative values are allowed when a lognormal distribution is assumed.");
				}
				gpdYd[i] = log(1.0 + gpdYd[i] / (gpdYm[i] * gpdYm[i]));
				gpdYm[i] = log(gpdYm[i]) - gpdYd[i]/2.0;
				gpdStDev[i] = sqrt(gpdYd[i]);
			}
		}
	}
	else 
	{
		gpdxxi = DVECTOR (1, gsExpoVars.iNbrObs_Total);
		gpdyyi = DVECTOR (1, gsExpoVars.iNbrObs_Total);
		gpdYm = DVECTOR (1, gsExpoVars.iNbrObs_Total);
		gpdYd = DVECTOR (1, gsExpoVars.iNbrObs_Total);
		gpdXi = DVECTOR (1, gsExpoVars.iNbrObs_Total);
		gpiNi = IVECTOR (1, gsExpoVars.iNbrObs_Total);
		gpdYsum = DVECTOR (1, gsExpoVars.iNbrObs_Total);

		double *xxi, *yyi;
		xxi = DVECTOR (1, gsExpoVars.iNbrObs_Total);
		yyi = DVECTOR (1, gsExpoVars.iNbrObs_Total);

		giNmiss = READ_OBSDATA2V(gsExpoVars.iNbrObs_Total, xxi, yyi);

#ifdef DO_LOG
		if (giDo_Log)	// Print values to log for investigation
		{
			fprintf(fp_log,"\n\nxxi[i]\n");
			for (i=1; i<=gsExpoVars.iNbrObs_Total; i++)
			{
				fprintf(fp_log,"xxi[%d]=%g\n",i,xxi[i]);
			}
			fprintf(fp_log,"\n\nyyi[i]\n");
			for (i=1; i<=gsExpoVars.iNbrObs_Total; i++)
			{
				fprintf(fp_log,"yyi[%d]=%g\n",i,yyi[i]);
			}
		}
#endif

		gsExpoVars.iNbrObs_Total -= giNmiss;  /* extern variable giNbr_Obs has been changed. */

		for (i=1; i<=gsExpoVars.iNbrObs_Total; i++) 
		{
			gpdXi[i] = gpdYsum[i] = gpdYm[i] = gpdYd[i] = 0;
			gpiNi[i] = 0;
		}

		Sort_2_By_Dose(gsExpoVars.iNbrObs_Total, xxi, yyi);  /*Sort gpdxxi and make appropriate changes to gpdyyi */

		for(i=1; i<=gsExpoVars.iNbrObs_Total; i++)
		{
			gpdxxi[i] = xxi[i];
			gpdyyi[i] = yyi[i];
			if ((gsExpoVars.iLogNormal == 1) && (gsExpoVars.iExactSolution == 1)) 
			{
				/* exact soln with individ data – need to add check: 
				if any gpdyyi is <=0 print error message on screen 
				and in output file: “Data set includes negative responses; 
				no negative values are allowed when a lognormal distribution is assumed.”  
				Then stop processing */
				if(gpdyyi[i] <= 0) //Exit
				{
					fprintf(fp_out, "\n\n\nData set includes negative responses:\nno negative values are allowed when a lognormal distribution is assumed."); 
					FreeUp_mem();
					ERRORPRT("Data set includes negative responses:\nno negative values are allowed when a lognormal distribution is assumed.");
				}
				gpdyyi[i] = log(gpdyyi[i]);
			}
		}

		//iGroup = 0;
		iGroup = 1;

		for (i=1; i<=gsExpoVars.iNbrObs_Total; i++) {
			gpdXi[iGroup]=gpdxxi[i];
			gpiNi[iGroup] += 1;
			gpdYsum[iGroup] += gpdyyi[i];
			if(i < gsExpoVars.iNbrObs_Total) {
				if(gpdxxi[i] != gpdxxi[i + 1])
					iGroup +=1;
			}
		}
		gsExpoVars.iNbrObs_Total=iGroup;

		jj=1;

		for (i=1; i<=gsExpoVars.iNbrObs_Total; i++) {
			gpdYm[i] = gpdYsum[i]/gpiNi[i];
			for (j=1; j<=gpiNi[i]; j++)
			{
				gpdYd[i] += (gpdyyi[jj]-gpdYm[i])*(gpdyyi[jj]-gpdYm[i])/(gpiNi[i]-1);
				jj += 1;
			}
			if ((gsExpoVars.iLogNormal == 1) && (gsExpoVars.iExactSolution == 0)) // lognormal, approx soln
			{
				/* need to add check: if any gpdYm[i] is <=0 
				print error message on screen and in output file: 
				Data set includes negative responses; 
				no negative values are allowed when a lognormal distribution is assumed.”  
				Then stop processing */
				if(gpdYm[i] <= 0) //Exit
				{
					fprintf(fp_out, "\n\n\nData set includes negative responses:\nno negative values are allowed when a lognormal distribution is assumed."); 
					FreeUp_mem();
					ERRORPRT("Data set includes negative responses:\nno negative values are allowed when a lognormal distribution is assumed.");
				}
				gpdYd[i] = log(1.0+ gpdYd[i] / (gpdYm[i] * gpdYm[i]));
				gpdYm[i] = log(gpdYm[i]) - gpdYd[i]/2.0;
				//gpdStDev[i] = sqrt(gpdYd[i]);
			}
		}
		FREE_DVECTOR (xxi, 1, gsExpoVars.iNbrObs_Total);
		FREE_DVECTOR (yyi, 1, gsExpoVars.iNbrObs_Total);

	} // end else (gsExpoVars.iIn_Type==1)

#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\ngpdYd[i]\n");
		for (i=1; i<=gsExpoVars.iNbrObs_Total; i++)
		{
			fprintf(fp_log,"gpdYd[%d]=%g\n",i,gpdYd[i]);
		}
		fprintf(fp_log,"\n\ngpdYm[i]\n");
		for (i=1; i<=gsExpoVars.iNbrObs_Total; i++)
		{
			fprintf(fp_log,"gpdYm[%d]=%g\n",i,gpdYm[i]);
		}
		if (gsExpoVars.iIn_Type==1)
		{
			fprintf(fp_log,"\n\ngpdStDev[i]\n");
			for (i=1; i<=gsExpoVars.iNbrObs_Total; i++)
			{
				fprintf(fp_log,"gpdStDev[%d]=%g\n",i,gpdStDev[i]);
			}
		}
	}
#endif

	gpdScXi = DVECTOR(1, gsExpoVars.iNbrObs_Total);
	gpdScYd = DVECTOR(1, gsExpoVars.iNbrObs_Total);
	gpdScYm = DVECTOR(1, gsExpoVars.iNbrObs_Total);
	gdmaxYd = gpdYd[1]; 
	gdmaxYm = gpdYm[1];

	for (i=1; i<=gsExpoVars.iNbrObs_Total; i++)
	{
		if(gpdYd[i] > gdmaxYd)
			gdmaxYd = gpdYd[i];
		if(gpdYm[i] > gdmaxYm)
			gdmaxYm = gpdYm[i];
	}

	for (i=1; i<=gsExpoVars.iNbrObs_Total; i++){
		gpdScYd[i] = gpdYd[i]/gdmaxYd;
		gpdScYm[i] = gpdYm[i]/gdmaxYm;
	}

	if (gsExpoVars.iIn_Type == 1) 
	{
		for (i=1; i<=gsExpoVars.iNbrObs_Total; i++) {
			if ( (gpdXi[i] < 0) || (gpiNi[i] < 0) || (gpdYd[i] < 0) )
				ERRORPRT("Values of dose, group size, and sample variance should be positive...");
		}
	}
	else {
		for (i=1; i<=gsExpoVars.iNbrObs_Total; i++) {
			if (gpdXi[i] < 0)
				ERRORPRT("Dose value should be positive ...");
		}
	}

	gdxmin = gdyymin = Max_double;
	gdxmax = gdyymax = 0.0;

	for (i=1;i<=gsExpoVars.iNbrObs_Total;i++){
		if (gpdXi[i] < gdxmin) {
			gdxmin = gpdXi[i];
			dbmr_root = sqrt(gpdYd[i]);
		}
		if (gpdXi[i] > gdxmax)
			gdxmax = gpdXi[i];
		if (gpdYm[i] > gdyymax)
			gdyymax = gpdYm[i];
		if (gpdYm[i] < gdyymin)
			gdyymin = gpdYm[i];
	}

	if ((gsExpoVars.iSign != 1) && (gsExpoVars.iSign != -1))	//Change || to && Bruce, 01/28/08 and moved from L817
	{
		gsExpoVars.iSign = Get_Linear_Trend(gsExpoVars.iNbrObs_Total, gpdXi, gpdYm, gpiNi);
		iTemp_sign = 0;
	}
	else 
	{
		iTemp_sign = 9999;
	}/* end if ((gsExpoVars.iSign != 1) || (gsExpoVars.iSign != -1)) */

	gdAdverseBMDEffect = gsExpoVars.dBmdEffect * gsExpoVars.iSign;

	// Added 1/28/2008, with Bruce
	for (i=1; i<=gsExpoVars.iNbrObs_Total; i++)
	{
		gpdScXi[i] = gpdXi[i]/gdxmax;
		gpdScYm[i] = gpdYm[i]/gdmaxYm;
		gpdScYd[i] = gpdYd[i]/(gdmaxYm*gdmaxYm);
	}

	//Compute for likelihoods
	dlikeA1=dlikeA2=dlikeA3=dvv=dlikeR=0.0;

	double dNtot=0;
	double dtemp=0.0;

	// Likelihood A1: Yij = Mu(i) + e(ij), Var{e(ij)} = Sigma^2
	for (i=1; i<=gsExpoVars.iNbrObs_Total; i++){
		dNtot += gpiNi[i];
	}

	for(i =1; i <= gsExpoVars.iNbrObs_Total; i++)
		dvv += ((gpiNi[i] - 1)*gpdYd[i]);

	dvv = dvv/dNtot;

	/* Compute Likelihood for model A1: Yij = Mu(i) + e(ij)
	Var(e(ij)) = Sigma^2 */

	dlikeA1 = - log(dvv)*dNtot/2 - dNtot/2;

	/* Compute Likelihood for model A2: Yij = Mu(i) + e(ij)
	Var(e(ij)) = Sigma(i)^2*/
	for(i =1; i <= gsExpoVars.iNbrObs_Total; i++)
	{
		dtemp = (gpiNi[i] - 1)*gpdYd[i];
		dtemp = dtemp/gpiNi[i];
		if (dtemp < 1E-20) dtemp = 0.000000001;
		dtemp = log(dtemp) + 1;
		dlikeA2 -= .5*gpiNi[i]*dtemp;
	}

	/* Compute Likelihood for model A3: Yij = Mu(i) + e(ij)
	Var(e(ij)) = k*(m(gpdXi))^rho*/

	/* Calculate A3 likelihood regardless of constant variance choice. */
	// Yij = Mu(i) + e(ij), Var{e(ij)} = exp(lalpha + log(mean(i)) * rho)
#ifdef DO_LOG
	if (giDo_Log && giDo_PrintLog)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\niCons_Var before AThree_fit\n");
		fprintf(fp_log,"   iCons_Var = %d\n", gsExpoVars.iCons_Var);
	}
#endif

	giNbrParm_A3=gsExpoVars.iNbrObs_Total+2;
	dlikeA3 = 0.0;
	for(i=1; i < 5; i++)
		gpdA3s[i] = 0.0;

	if (gsExpoVars.iCons_Var == 0 ) { // Parameters for fitting the model above
#ifdef DO_LOG
		if (giDo_Log && giDo_PrintLog)	// Print values to log for investigation
		{
			fprintf(fp_log,"\n\nList of Data before AThree_fit\n");
			for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
			{
				fprintf(fp_log,"     %5.4g%7d%14.3g%14.3g\n", gpdXi[i], gpiNi[i], gpdYm[i], gpdStDev[i]);
			}
		}
#endif
		double *pdLKParms;
		// Test if all models have the same specs
		//int piSum[] = {0,0,0,0};
		int iVar, iSum, iSame;
		iSame = iSum = iVar = 0;
		for(i = 1; i < 5; i++)
			iSum += (gppiSpecPara[i][1] + gppiSpecPara[i][2]);

		//for(i = 1; i < 5; i++)
		//	iSum += piSum[i];

		if(iSum == 0 || iSum == 8) //the same specs for 4 models
			iSame = 1;
		else
		{
			if(iSum == 4)
			{
				for(i = 1; i < 5; i++)
					iVar = gppiSpecPara[i][1];

				if(iSum == iVar || iVar == 0) //Alphas are set to 1 OR Rhos are 0 and all models have same specs
					iSame = 1;
			}
		}
		if(iSame == 0)
			iSame = 4;

		pdLKParms = DVECTOR(1, giNbrParm_A3);

		/******* Fit the A3 model above *********/
		for(i = 1; i<= iSame; i++)
		{
			for(j=1; j <= giNbrParm_A3; j++)
				pdLKParms[j]=0.0;

			AThree_Fit(giNbrParm_A3, pdLKParms, EPS, &gpiIter, i);
		}
		dlikeA3 = gpdA3s[1];
#ifdef DO_LOG
		if (giDo_Log && giDo_PrintLog)	// Print values to log for investigation
		{
			for (i=1; i<=giNbrParm_A3; i++) {
				fprintf(fp_log, "LKParms[%d] = %.8f\n", i, pdLKParms[i]);
			}
			fprintf(fp_log,"\n\nList of Data after AThree_fit\n");
			for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
			{
				fprintf(fp_log,"     %5.4g%7d%14.3g%14.3g\n", gpdXi[i], gpiNi[i], gpdYm[i], gpdStDev[i]);
			}
		}
#endif
		//int iter;
		//AThree_Fit(giNbrParm_A3, pdLKParms, EPS, &iter, &lkA3, nModel);

		FREE_DVECTOR (pdLKParms, 1, giNbrParm_A3);
	}
	/* If constant variance model then just set equal to model A1. */
	else {
		dlikeA3 = dlikeA1;
	}

	///* Calculate "reduced" model likelihood R:  Yi = Mu + e(i)
	//Var{e(i)} = Sigma^2 */
	dlikeR = 0.0;
	double dYbar = 0.0;

	/* Model R: common mean and variance. */
	dYbar = gpdYm[1] * gpiNi[1];

	for (i = 2; i <= gsExpoVars.iNbrObs_Total; i++)
		dYbar += gpdYm[i] * gpiNi[i];

	dYbar = dYbar / dNtot;

	dvv = gpdYd[1] * (gpiNi[1] - 1) + gpiNi[1]*(gpdYm[1] - dYbar)*(gpdYm[1] - dYbar);

	for (i = 2; i <= gsExpoVars.iNbrObs_Total; i++)
		dvv += gpdYd[i] * (gpiNi[i] - 1) + gpiNi[i]*(gpdYm[i] - dYbar)*(gpdYm[i] - dYbar);

	dvv = dvv / dNtot;

	dlikeR = -dNtot * (1.0 + log(dvv))/2.0;

	if (!gppdMLEs || !gpdLikelihoods) 
	{
#ifdef DO_LOG
		if (giDo_Log) 
		{
			fflush(fp_log);
			fclose(fp_log);
		}
#endif
		ERRORPRT ("Memory allocation failed in sMLEs");
	}

	//INITIALIZE_DMATRIX(gppdMLEs, NBR_OF_MODELS, NBR_OF_PARMS);
	for(i=1; i<= NBR_OF_MODELS; i++)
	{
		for(j=1; j<=NBR_OF_PARMS; j++)
			gppdMLEs[i][j] = 0.0;

		for(j=1; j<=gsExpoVars.iNbrObs_Total ; j++)
			gppdPredict[i][j] = 0.0;

		gpdLikelihoods[i] = -Max_double;
	}

	/* Main loop to process Exponential Models*/
#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\nMAIN LOOP\n");
	}
#endif

	// NOTE: This computation does take into account user parameter input yet.
	for( i=iStart; i <= iEnd; i++)
	{
#ifdef DO_LOG
		if (giDo_Log)	// Print values to log for investigation
		{
			fprintf(fp_log,"\n\nModel Run #%d BEFORE Processing, giRun[%d]=%d: ", i+1, i, giRun[i]);
		}
#endif

		if(giRun[i] == 0)
			continue;

		for(j=1; j<=NBR_OF_PARMS; j++)
		{
			//if specified, get the specified values to gppdInitPara array, 
			//else if initialized just use what were read into the gppdInitPara
			if(gppiSpecPara[i][j] == 1)
				gppdInitPara[i][j] = gppdSpecPara[i][j];
		}

#ifdef DO_LOG
		if (giDo_Log)	// Print values to log for investigation
		{
			fprintf(fp_log,"\n\nBEFORE CASE, Model %d Initial Values: ", i+1);
			for(j=1; j <= NBR_OF_PARMS; j++)
			{
				fprintf(fp_log," %g", gppdInitPara[i][j]);
			}
			fprintf(fp_log,"\n\nModel %d MLEs: ", i+1);
			for(j=1; j <= NBR_OF_PARMS; j++)
			{
				fprintf(fp_log," %g", gppdMLEs[i][j]);
			}
		}
#endif

		if(giRun[i] == 0)
			goto DONE;

		switch(i)
		{
		case (int)eM2:	//Model 2

			gsExpoVars.iSelect = i+1;

			if(gpiInitialized[i] == 0) //Not initialized generate initial value for Model 2
				DoInitValuesM2(gppdInitPara[i]);

			for(j=1; j<= NBR_OF_PARMS; j++)
				gppdMLEs[i][j] = gppdInitPara[i][j];

			gpdLikelihoods[i] = LogLike(i, gppdInitPara[i]); 
#ifdef DO_LOG
			if (giDo_Log)	// Print values to log for investigation
			{
				fprintf(fp_log,"\n\nBEFORE Exponential_fit CALL, gpdLikelihoods[i] = %g", gpdLikelihoods[i]);
			}
#endif

			gpimodtype = 3;

			Exponential_fit(i);
			gpdLikelihoods[i] = LogLike(i, gppdMLEs[i]);
		
			Exp_vcv(gppdVCV, i);
			Get_DTMSVCV(i);
			Calc_ParmSE(i);



#ifdef DO_LOG
			if (giDo_Log)	// Print values to log for investigation
			{
				fprintf(fp_log,"\n\nAFTER Exponential_fit CALL, gpdLikelihoods[i] = %g", gpdLikelihoods[i]);
			}
#endif

			break;

		case (int)eM3:	//Model 3

			gsExpoVars.iSelect = i+1;

			if(gpiInitialized[i] == 0) //Not initialized generate initial value for Model 2
				DoInitValuesM2(gppdInitPara[i]);

			for(j=1; j<= NBR_OF_PARMS; j++)
				gppdMLEs[i][j] = gppdInitPara[i][j];

			gpdLikelihoods[i] = LogLike(i, gppdInitPara[i]);
#ifdef DO_LOG
			if (giDo_Log)	// Print values to log for investigation
			{
				fprintf(fp_log,"\n\nBEFORE Exponential_fit CALL, gpdLikelihoods[i] = %g", gpdLikelihoods[i]);
			}
#endif
			gpimodtype = 4;

			Exponential_fit(i);
			gpdLikelihoods[i] = LogLike(i, gppdMLEs[i]);

			Exp_vcv(gppdVCV, i);
			Get_DTMSVCV(i);
			Calc_ParmSE(i);
			



#ifdef DO_LOG
			if (giDo_Log)	// Print values to log for investigation
			{
				fprintf(fp_log,"\n\nAFTER Exponential_fit CALL, gpdLikelihoods[i] = %g", gpdLikelihoods[i]);
			}
#endif
			break;

		case (int)eM4:	//Model 4

			gsExpoVars.iSelect = i+1;

			if(gpiInitialized[i] == 0)
				DoInitValuesM4(gppdInitPara[i]);


			for(j=1; j<= NBR_OF_PARMS; j++)
				gppdMLEs[i][j] = gppdInitPara[i][j];

			gpdLikelihoods[i] = LogLike(i, gppdInitPara[i]);
#ifdef DO_LOG
			if (giDo_Log)	// Print values to log for investigation
			{
				fprintf(fp_log,"\n\nBEFORE Exponential_fit CALL, gpdLikelihoods[i] = %g", gpdLikelihoods[i]);
			}
#endif
			gpimodtype = 5;

			Exponential_fit(i);
			gpdLikelihoods[i] = LogLike(i, gppdMLEs[i]);

			Exp_vcv(gppdVCV, i);
			Get_DTMSVCV(i);
			Calc_ParmSE(i);



#ifdef DO_LOG
			if (giDo_Log)	// Print values to log for investigation
			{
				fprintf(fp_log,"\n\nAFTER Exponential_fit CALL, gpdLikelihoods[i] = %g", gpdLikelihoods[i]);
			}
#endif
			break;

		case (int)eM5:	//Model 5

			gsExpoVars.iSelect = i+1;

			if(gpiInitialized[i] == 0)
				DoInitValuesM4(gppdInitPara[i]);


			for(j=1; j<= NBR_OF_PARMS; j++)
				gppdMLEs[i][j] = gppdInitPara[i][j];

			gpdLikelihoods[i] = LogLike(i, gppdInitPara[i]);
#ifdef DO_LOG
			if (giDo_Log)	// Print values to log for investigation
			{
				fprintf(fp_log,"\n\nBEFORE Exponential_fit CALL, gpdLikelihoods[i] = %g", gpdLikelihoods[i]);
			}
#endif
			gpimodtype = 6;

			Exponential_fit(i);
			gpdLikelihoods[i] = LogLike(i, gppdMLEs[i]);

			Exp_vcv(gppdVCV, i);
			Get_DTMSVCV(i);
			Calc_ParmSE(i);




#ifdef DO_LOG
			if (giDo_Log)	// Print values to log for investigation
			{
				fprintf(fp_log,"\n\nAFTER Exponential_fit CALL, gpdLikelihoods[i] = %g", gpdLikelihoods[i]);
			}
#endif
			break;
		}

#ifdef DO_LOG
		if (giDo_Log)	// Print values to log for investigation
		{
			fprintf(fp_log,"\n\nLine 1038, Model Run #%d AFTER Processing, giRun[%d]=%d: \n\n", gsExpoVars.iSelect, i, giRun[i]);

			fprintf(fp_log,"\n\nAFTER CASE, Model %d Initial Values:  ", i+1);
			for(j=1; j <= NBR_OF_PARMS; j++)
			{
				fprintf(fp_log,"  %g", gppdInitPara[i][j]);
			}
			fprintf(fp_log,"\nModel %d MLEs: ", i+1);
			for(j=1; j <= NBR_OF_PARMS; j++)
			{
				fprintf(fp_log,"\ngppdMLEs[%d][%d] = %g", i, j, gppdMLEs[i][j]);
			}
			fflush(fp_log);
		}
#endif

	} //End of for( i=iStart; i <= iEnd; i++), main loop

	/* Print model and file information on output page */
#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\nModel Runs\n");
		for(i=1; i<=NBR_OF_MODELS; i++)
		{
			fprintf(fp_log,"   \nRun[%d] = %d\n", i+1, giRun[i]);
		}
	}
#endif


	if(iModelSpan == 1)
	{
		HeaderAndReadVar_2Out(argv[1], ctime(&ltime), 0);
		DoInitParam(0);
		DoParamEstimates(0);
		DoDataEstimateInt(0);
		DoLikehoods(dlikeA1, dlikeA2, dlikeA3, dlikeR, 0);
		DoBenchMark(iModelSpan, 0);
	}
	else
	{
		int i, j, nLastBS, nSrcLen;
		nLastBS = 0;
		char sPath[255];
		char modelName[255];
		char outFile[255];
		strcpy(sPath,gacFileOut);
		nSrcLen = strlen(gacFileOut);
#ifdef DO_LOG
		if (giDo_Log == true)
		{
			fprintf(fp_log,"\n\n gagFileOut=%s, sPath=%s, nSrcLen=%d", gacFileOut, sPath, nSrcLen);
		}
#endif

		for(i = 0; i < nSrcLen; i++) //Test for backslash
		{
			if(gacFileOut[i] == '\\')
				if(i > nLastBS)
					nLastBS = i;
		}

		if(nLastBS == 0) //Test for forward slash
		{
			for(i = 0; i < nSrcLen; i++)
			{
				if(gacFileOut[i] == '/')
					if(i > nLastBS)
						nLastBS = i;
			}
		}
//#ifdef DO_LOG
//		if (giDo_Log == true)
//		{
//			fprintf(fp_log,"\n\n nLastBS = %d", nLastBS);
//		}
//#endif

		if(nLastBS > 0 || nSrcLen > 0)
		{
			sPath[nLastBS+(nLastBS>0?1:0)] = '\0';
			modelName[0] = 'M';
			modelName[1] = '1';
			i = nLastBS;
			j = 2;
			for(i = nLastBS+(nLastBS>0?1:0); i < nSrcLen; i++)
			{
				if(gacFileOut[i] == '\0')
					break;
				modelName[j] = gacFileOut[i];
				j++;
			}
			modelName[j] = '\0';
		}
		else
			goto DONE;

		if(fp_out != NULL)
			fclose(fp_out);
#ifdef DO_LOG
		if (giDo_Log == true)
		{
			fprintf(fp_log,"\n\niStart = %d, iEnd = %d\n", iStart, iEnd);
		}
#endif

		for( i=iStart; i <= iEnd; i++)  //Start debug here
		{
			if(giRun[i] == 0)
				continue;

			if(fp_out != NULL)
			{
				fflush(fp_out);
				fclose(fp_out);
			}
			if(i == 1)
				modelName[1] = '2';
			else if(i == 2)
				modelName[1] = '3';
			else if (i == 3)
				modelName[1] = '4';
			else if(i == 4)
				modelName[1] = '5';

			strcpy(outFile, sPath);
			strcat(outFile, modelName);
			fp_out=fopen(outFile,"w");

			if(fp_out == NULL)
				continue;

			HeaderAndReadVar_2Out(argv[1], ctime(&ltime), i);
			DoInitParam(i);
			DoParamEstimates(i);
			DoDataEstimateInt(i);
			DoLikehoods(dlikeA1, dlikeA2, dlikeA3, dlikeR, i);
			DoBenchMark(iModelSpan, i);
		}
	}

	//Print the .002 files
	i = 0;
	if(doDot002(iStart, iEnd, iModelSpan)==0)
		i = 1;

DONE:
	FreeUp_mem();

#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\n\nEND OF MAIN LOOP");
		fflush(fp_log);
		fclose(fp_log);
	}
#endif

	CLOSE_FILES ();
	return(i);
}

/**********************************************************
** READ_OBSDATA2V--used to read 2 column data in 2 vectors.
***********************************************************/
int READ_OBSDATA2V(int iNTotal, double pdxxi[], double pdyyi[])

{
	int     iNmiss;          /*number of records with missing values*/
	int     i,j,n,m;        /*count and iteration control variables*/
	double  value;          /*temp variable*/

	iNmiss = 0;
	for(i=1; i<=iNTotal; i++){
		n = i-iNmiss;
		m = 0;
		for (j=1;j<=2;j++){
			fscanf(fp_in,"%lf",&value);
			if (value != MISSING){
				if (j==1)     pdxxi[n]=value;
				if (j==2)     pdyyi[n]=value;
			} 
			else m++;
		}
		if (m != 0)
			iNmiss++;
		else if(pdxxi[n] < 0)
			iNmiss++;
	}
	return iNmiss;
}

void HeaderAndReadVar_2Out(char *dFileName, char* clocktime, int iModelNbr)
{
	strcpy(gacPltFileName,"");
	Output_Header(Version_no, dFileName, gacPltFileName, clocktime, gsExpoVars.caUser_note);

	/* output title and summary of input data */
	OUTPUT_TEXT("\n   The form of the response function by Model: ");
	OUTPUT_TEXT("      Model 2:     Y[dose] = a * exp{sign * b * dose}");  //took out - sign on b, 01/28/08
	OUTPUT_TEXT("      Model 3:     Y[dose] = a * exp{sign * (b * dose)^d}"); //took out - sign on b, 01/28/08
	OUTPUT_TEXT("      Model 4:     Y[dose] = a * [c-(c-1) * exp{-b * dose}]");
	OUTPUT_TEXT("      Model 5:     Y[dose] = a * [c-(c-1) * exp{-(b * dose)^d}]");

	OUTPUT_TEXT("\n    Note: Y[dose] is the median response for exposure = dose;");
	OUTPUT_TEXT("          sign = +1 for increasing trend in data;");
	OUTPUT_TEXT("          sign = -1 for decreasing trend.");

	OUTPUT_TEXT("\n      Model 2 is nested within Models 3 and 4.");
	OUTPUT_TEXT("      Model 3 is nested within Model 5.");
	OUTPUT_TEXT("      Model 4 is nested within Model 5.");

	if (gsExpoVars.iIn_Type == 1)
		fprintf(fp_out,"\n\n   Dependent variable = %s", gsExpoVars.iLogNormal==1?"Calculated Median":gsExpoVars.caMeanName);
	else
		fprintf(fp_out,"\n\n   Dependent variable = %s", gsExpoVars.iLogNormal==1?"Calculated Median":gsExpoVars.caResponseName);

	fprintf(fp_out,"\n   Independent variable = %s", gsExpoVars.caDoseName);
	fprintf(fp_out,"\n   %s", gsExpoVars.iLogNormal==1?"Data are assumed to be distributed: lognormally ":"Data are assumed to be distributed: normally");
	fprintf(fp_out,"\n   Variance Model: %s", gsExpoVars.iLogNormal==0?"exp(lnalpha +rho *ln(Y[dose]))":"Log-scale variance = exp(lnalpha)");

	if(gsExpoVars.iCons_Var == 0)
		fprintf(fp_out,"\n   The variance is to be modeled as Var(i) = exp(lalpha + log(mean(i)) * rho)");
	else
		fprintf(fp_out,"\n   rho is set to 0.\n   A constant%s variance model is fit.", gsExpoVars.iLogNormal==1?" log-scale":"");

	fprintf(fp_out,"\n\n   Total number of dose groups = %d", gsExpoVars.iIn_Type==1 ? gsExpoVars.iNbrObs_Total-giNmiss : gsExpoVars.iNbrObs_Total-giNmiss);
	fprintf(fp_out,"\n   Total number of records with missing values = %d", giNmiss);
	fprintf(fp_out,"\n   Maximum number of iterations = %d", gsExpoVars.iMaxIter);
	fprintf(fp_out,"\n   Relative Function Convergence has been set to: %g", gsExpoVars.dRel_Conv);
	fprintf(fp_out,"\n   Parameter Convergence has been set to: %g\n", gsExpoVars.dParm_Conv);
	fprintf(fp_out,"\n   MLE solution provided: %s\n", gsExpoVars.iExactSolution==0?"Approximate":"Exact");
}


/************************************************/
//  DoInitValuesM2(double p[]) -- modified to include power so can be used to search for good start for model 3 also
/***********************************************/

void DoInitValuesM2(double p[])

{
	double dNtot=0.0;
	//, Nd=0, Nc=0;
	double y_bar=0., x_bar=0., stot=0.;
	double ssxy=0., ssx=0.;
	double sigma_bar=0., Ym_bar=0.;
	double ss_Ym=0., ss_YmSigma=0.;
	double ll=-Max_double, maxll= -Max_double;
	int i, j;
	double *ptemp;
	ptemp = DVECTOR(1, NBR_OF_PARMS);
	int n = 0;

#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\nList of Parms inside DoInitValuesM2\n");
		for(i=1; i <= NBR_OF_PARMS; i++)
		{
			fprintf(fp_log,"p[%1d] = %4.4g\n", i, p[i]);
		}
	}
#endif

	for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i){

		x_bar+=(double)gpiNi[i]*gpdXi[i]; //Changed BACK TO gpdXi, Bruce 08/11/08

		if(gsExpoVars.iLogNormal == 0)
			y_bar+=(double)gpiNi[i]*Slog(gpdYm[i]);
		else
			y_bar+=(double)gpiNi[i]*gpdYm[i];

		dNtot+=gpiNi[i];
		stot+=((double)gpiNi[i]-1.)*gpdYd[i];
	}

	y_bar=y_bar/dNtot;
	x_bar=x_bar/dNtot;

	// initial values for alpha and rho  
	if(gsExpoVars.iCons_Var==1)
	{
		p[(int)eAlpha]=(double)log(stot/dNtot); //added log, 01/28/08
		p[(int)eRho]=0.;
	} 
	else 
	{

		for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
		{
			sigma_bar+=Slog(gpdYd[i]);
			Ym_bar+=Slog(gpdYm[i]);
		}

		sigma_bar=sigma_bar/gsExpoVars.iNbrObs_Total;
		Ym_bar=Ym_bar/gsExpoVars.iNbrObs_Total;

		for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
		{
			ss_Ym +=(Slog(gpdYm[i])-Ym_bar)*(Slog(gpdYm[i])-Ym_bar);
			ss_YmSigma+=(Slog(gpdYm[i])-Ym_bar)*Slog(gpdYd[i]);
		}

		p[(int)eRho]=ss_YmSigma/ss_Ym;		//1
		p[(int)eAlpha]=(sigma_bar-p[(int)eRho]*Ym_bar);	//0, del exp, 01/28/08
	}

	// initial values for a, b


	for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i){

		ssx +=(double)gpiNi[i]*(gpdXi[i]-x_bar)*(gpdXi[i]-x_bar); //Changed back to gpdXi, Bruce 08/11/08

		if(gsExpoVars.iLogNormal == 0)
			ssxy+=(double)gpiNi[i]*(gpdXi[i]-x_bar)*Slog(gpdYm[i]); //Changed back to gpdXi, Bruce 08/11/08
		else
			ssxy+=(double)gpiNi[i]*(gpdXi[i]-x_bar)*gpdYm[i];
	}



	//Alpha, Rho, a, b, c, d
	//    0,   1, 2, 3, 4, 5
	p[(int)eb]=(double)gsExpoVars.iSign*ssxy/ssx;

	p[(int)ea]=exp(y_bar-p[(int)eb]*x_bar);  //a

	p[(int)ec]= 0.0;
	p[(int)ed]= 1.0;

	maxll = LogLike(1, p);

	// for model 3
#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n(inside DoInitValuesM2)gsExpoVars.iSelect = %d\n", gsExpoVars.iSelect);
	}
#endif
	if (gsExpoVars.iSelect != 2)
	{
		ptemp[(int)eAlpha] = p[(int)eAlpha];
		ptemp[(int)eRho] = p[(int)eRho];

		for (j=2; j<=18; j++)
		{
			ptemp[(int)ed]= (double)j;

			for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
			{
				x_bar+=(double)gpiNi[i]*pow(gpdXi[i],ptemp[(int)ed]); 

				if(gsExpoVars.iLogNormal == 0)
					y_bar+=(double)gpiNi[i]*Slog(gpdYm[i]);
				else
					y_bar+=(double)gpiNi[i]*gpdYm[i];

				dNtot+=gpiNi[i];
				stot+=((double)gpiNi[i]-1.)*gpdYd[i];
			}

			y_bar=y_bar/dNtot;
			x_bar=x_bar/dNtot;

			for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
			{
				ssx +=(double)gpiNi[i]*(pow(gpdXi[i],ptemp[(int)ed])-x_bar)*(pow(gpdXi[i],ptemp[(int)ed])-x_bar); 

				if(gsExpoVars.iLogNormal == 0)
					ssxy+=(double)gpiNi[i]*(pow(gpdXi[i],ptemp[(int)ed])-x_bar)*Slog(gpdYm[i]);
				else
					ssxy+=(double)gpiNi[i]*(pow(gpdXi[i],ptemp[(int)ed])-x_bar)*gpdYm[i];
			}

			ptemp[(int)eb]=(double)pow((gsExpoVars.iSign*ssxy/ssx), (1/p[(int)ed]));  //b

			ptemp[(int)ea]=exp(y_bar-ptemp[(int)eb]*x_bar);  //a

			ptemp[(int)ec]= 0.0;

			ll = LogLike(2, ptemp);

			if (ll > maxll)
			{
				maxll = ll;
				for (i=(int)eAlpha; i <=(int)ed; i++)
				{
					p[i] = ptemp[i];
				}
			}

#ifdef DO_LOG
			if (giDo_Log)	// Print values to log for investigation
			{
				fprintf(fp_log,"\nloop j=%d, ll=%g, maxll=%g\n", j, ll, maxll);
			}
#endif
			for (n = 1; n <= NBR_OF_PARMS; n++)
			{
#ifdef DO_LOG
				if (giDo_Log)	// Print values to log for investigation
				{
					fprintf(fp_log,"ptemp[%d] = %g\n", n, ptemp[n]);
				}
#endif
			}
		}
	}
	FREE_DVECTOR(ptemp, 1, NBR_OF_PARMS);
}


/*****************************************************/
//DoInitValuesM4(double p[])
/*****************************************************/
void DoInitValuesM4(double p[])
{
	int i, j, iLowest=0, iGreatest=0;
	double *Pred;
	Pred = DVECTOR(1, gsExpoVars.iNbrObs_Total);

	double dLowest, dGreatest, dmG, dma, dmc, dmb, dmult, dlltemp, dll, dxlk;

	dLowest=dGreatest=dmG=dma=dmc=dmb=dmult=dlltemp=dxlk=0.0;

	dll=-Max_double;

	dLowest = Max_double;
	dGreatest = -1;
	for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
	{
		if(dLowest >= gpdYm[i])
		{
			dLowest = gpdYm[i];
			iLowest = i;
		}
		if(dGreatest <= gpdYm[i])
		{
			dGreatest = gpdYm[i];
			iGreatest = i;
		}
	}

#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\nBEFORE CHANGE IN DoInitValuesM4(double p[])\n");
		for(i=(int)eAlpha; i <= NBR_OF_PARMS; i++)
		{
			fprintf(fp_log,"p[%d] = %g,  ", i, p[i]);
		}
	}
#endif

	//Code from DoInitValuesM2()
	double dNtot=0.0;
	double y_bar=0., x_bar=0., stot=0.;
	double ssxy=0., ssx=0.;
	double sigma_bar=0., Ym_bar=0.;
	double ss_Ym=0., ss_YmSigma=0.;

	for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
	{

		x_bar+=(double)gpiNi[i]*gpdXi[i]; //Changed to gpdXi 

		if(gsExpoVars.iLogNormal == 0)
			y_bar+=(double)gpiNi[i]*Slog(gpdYm[i]);
		else
			y_bar+=(double)gpiNi[i]*gpdYm[i];

		dNtot+=gpiNi[i];
		stot+=((double)gpiNi[i]-1.)*gpdYd[i];
	}

	y_bar=y_bar/dNtot;
	x_bar=x_bar/dNtot;

	for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
	{
		ssx +=(double)gpiNi[i]*(gpdXi[i]-x_bar)*(gpdXi[i]-x_bar); //Changed to gpdXi 
		if(gsExpoVars.iLogNormal == 0)
			ssxy+=(double)gpiNi[i]*(gpdXi[i]-x_bar)*Slog(gpdYm[i]); //Changed to gpdXi
		else
			ssxy+=(double)gpiNi[i]*(gpdXi[i]-x_bar)*gpdYm[i];
	}

	// initial values for alpha and rho  
	if(gsExpoVars.iCons_Var==1)
	{
		p[(int)eAlpha]=(double)log(stot/dNtot); //added log, 01/28/08
		p[(int)eRho]=0.;
	} 
	else 
	{

		for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
		{
			sigma_bar+=Slog(gpdYd[i]);
			Ym_bar+=Slog(gpdYm[i]);
		}

		sigma_bar=sigma_bar/gsExpoVars.iNbrObs_Total;
		Ym_bar=Ym_bar/gsExpoVars.iNbrObs_Total;

		for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
		{
			ss_Ym +=(Slog(gpdYm[i])-Ym_bar)*(Slog(gpdYm[i])-Ym_bar);
			ss_YmSigma+=(Slog(gpdYm[i])-Ym_bar)*Slog(gpdYd[i]);
		}

		p[(int)eRho]=ss_YmSigma/ss_Ym;		//1
		p[(int)eAlpha]=(sigma_bar-p[(int)eRho]*Ym_bar);	//0, del exp, 01/28/08
	}

	//End of Code from DoInitValuesM2()


	//for model 4
	if (gsExpoVars.iSign == 1)
	{
#ifdef DO_LOG
		if (giDo_Log)	// Print values to log for investigation
		{
			fprintf(fp_log,"\nInside DoInitValuesM4() loop (gsExpoVars.iSign == 1)");
		}
#endif

		for (j = 1; j <= 10; j++)
		{
			dxlk = 0.0;
			dmult = (j == 1 ? 1.05 : (j == 2 ? 2 : (j == 3 ? 5 : (j == 4 ? 10 : (j == 5 ? 20 : (j == 6 ? 50 : (j == 7 ? 100 : (j == 8 ? 200 : (j == 9 ? 500 : 1000)))))))));

			if(gsExpoVars.iLogNormal == 0)
				dma = 0.95*gpdYm[iLowest];
			else
				dma = 0.95*exp(gpdYm[iLowest]);

			if(gsExpoVars.iLogNormal == 0)
				dmG = dmult*gpdYm[iGreatest];
			else
				dmG = dmult*exp(gpdYm[iGreatest]);


			//Compute for b:
			dGreatest	= 0.0;	//Borrowed as the numerator
			dLowest		= 0.0;	//Borrowed as the denominator

			for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
			{
				if(gsExpoVars.iLogNormal == 0)
					dGreatest += (double)( gpdXi[i] * Slog( (dmG - gpdYm[i])/(dmG - dma) ) );
				else
					dGreatest += (double) (gpdXi[i] * Slog( (dmG - exp(gpdYm[i]))/(dmG-dma)));

				dLowest += (double)(gpdXi[i] * gpdXi[i]);

			}

			dmb = -1.0*(double)(dGreatest/dLowest);
			dmc = dmG/dma;

			for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
			{
				Pred[i] = dma*(dmc-(dmc-1)*exp(-1.0*(dmb*gpdXi[i])));
				if(gsExpoVars.iLogNormal == 1)
					Pred[i] = Slog(Pred[i]);

				dxlk+= 0.5*(double)gpiNi[i]*(p[(int)eAlpha]+p[(int)eRho]*Slog(Pred[i])) +0.5*(((double)gpiNi[i]-1)*gpdYd[i]
				+ (double)gpiNi[i]*(gpdYm[i]-Pred[i])*(gpdYm[i]-Pred[i]))/exp(p[(int)eAlpha]+p[(int)eRho]*Slog(Pred[i]));
			}

			dlltemp=-dxlk;
			if (dlltemp > dll)
			{
				dll = dlltemp;
				p[(int)ea] = dma;			// a
				p[(int)ec] = dmc;			// c
				p[(int)eb] = dmb;			// b
				p[(int)ed] = 1.0;			// d
			}
#ifdef DO_LOG
			if (giDo_Log)	// Print values to log for investigation
			{
				fprintf(fp_log,"\ndma=%g, dmult=%g, dmb=%g, dmc=%g, dxlk=%g, dll=%g", dma, dmult, dmb, dmc, dxlk, dll);
			}
#endif
		}
	}
	else
	{
#ifdef DO_LOG
		if (giDo_Log)	// Print values to log for investigation
		{
			fprintf(fp_log,"\nInside DoInitValuesM4() Line 1702 loop (gsExpoVars.iSign != 1)");
		}
#endif
		for (j=1; j <= 10; j++)
		{
			dxlk = 0.0;
			dmult = (j == 1 ? 1.05 : (j == 2 ? 2 : (j == 3 ? 5 : (j == 4 ? 10 : (j == 5 ? 20 : (j == 6 ? 50 : (j == 7 ? 100 : (j == 8 ? 200 : (j == 9 ? 500 : 1000)))))))));

			if(gsExpoVars.iLogNormal == 0)
				dmG = gpdYm[iLowest]/dmult;
			else
				dmG = exp(gpdYm[iLowest])/dmult;

			if(gsExpoVars.iLogNormal == 0)
				dma = 1.05 * gpdYm[iGreatest];
			else
				dma = 1.05 * exp(gpdYm[iGreatest]);

			//Compute for b:
			dGreatest	= 0.0;	//Borrowed as the numerator
			dLowest		= 0.0;	//Borrowed as the denominator

			for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
			{
				if(gsExpoVars.iLogNormal == 0)
					dGreatest += (double)( gpdXi[i] * Slog( (dmG - gpdYm[i])/(dmG - dma) ) );
				else
					dGreatest += (double)( gpdXi[i] * Slog( (dmG - exp(gpdYm[i]))/(dmG -dma)));

				dLowest += (double)(gpdXi[i] * gpdXi[i]);

			}

			dmb = -1.0*(double)(dGreatest/dLowest);
			dmc = dmG/dma;

			for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
			{
				Pred[i] = dma*(dmc-(dmc-1)*exp(-1.0*(dmb*gpdXi[i])));
				if(gsExpoVars.iLogNormal == 1)
					Pred[i] = Slog(Pred[i]);

				dxlk+= 0.5*(double)gpiNi[i]*(p[(int)eAlpha]+p[(int)eRho]*Slog(Pred[i])) +0.5*(((double)gpiNi[i]-1)*gpdYd[i]
				+ (double)gpiNi[i]*(gpdYm[i]-Pred[i])*(gpdYm[i]-Pred[i]))/exp(p[(int)eAlpha]+p[(int)eRho]*Slog(Pred[i]));
			}

			dlltemp=-dxlk;
			if (dlltemp > dll)
			{
				dll = dlltemp;
				p[(int)ea] = dma;			// a
				p[(int)ec] = dmc;			// c
				p[(int)eb] = dmb;			// b
				p[(int)ed] = 1.0;			// d
			}
#ifdef DO_LOG
			if (giDo_Log)	// Print values to log for investigation
			{
				fprintf(fp_log,"\ndma=%g, dmult=%g, dmb=%g, dmc=%g, dxlk=%g, dll=%g", dma, dmult, dmb, dmc, dxlk, dll);
			}
#endif
		}	
	}
#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\nAFTER CHANGE IN DoInitValuesM4(double p[])\n");
		for(i=(int)eAlpha; i <= NBR_OF_PARMS; i++)
		{
			fprintf(fp_log,"p[%d] = %g,  ", i, p[i]);
		}
	}
#endif

	// for Model 5
	int k = 0;
	double dmd = 0.0;
	if (gsExpoVars.iSelect != 4)
	{
		for (k = 2; k<=18; k++)
		{
			dmd = (double)k;
			if (gsExpoVars.iSign == 1)
			{
#ifdef DO_LOG
				if (giDo_Log)	// Print values to log for investigation
				{
					fprintf(fp_log,"\nInside DoInitValuesM4() loop (gsExpoVars.iSign == 1)");
				}
#endif

				for (j = 1; j <= 10; j++)
				{
					dxlk = 0.0;
					dmult = (j == 1 ? 1.05 : (j == 2 ? 2 : (j == 3 ? 5 : (j == 4 ? 10 : (j == 5 ? 20 : (j == 6 ? 50 : (j == 7 ? 100 : (j == 8 ? 200 : (j == 9 ? 500 : 1000)))))))));

					if(gsExpoVars.iLogNormal == 0)
						dma = 0.95*gpdYm[iLowest];
					else
						dma = 0.95*exp(gpdYm[iLowest]);

					if(gsExpoVars.iLogNormal == 0)
						dmG = dmult*gpdYm[iGreatest];
					else
						dmG = dmult*exp(gpdYm[iGreatest]);

					//Compute for b:
					dGreatest	= 0.0;	//Borrowed as the numerator
					dLowest		= 0.0;	//Borrowed as the denominator

					for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
					{
						if(gsExpoVars.iLogNormal == 0)
							dGreatest += (double)pow((-1.0*(gpdXi[i] * Slog( (dmG - gpdYm[i])/(dmG - dma)))), (1.0/dmd));
						else
							dGreatest += (double)pow((-1.0*(gpdXi[i] * Slog( (dmG - exp(gpdYm[i]))/(dmG - dma)))), (1.0/dmd));

						dLowest += (double)(gpdXi[i] * gpdXi[i]);
					}

					dmb = (double)(dGreatest/dLowest);
					dmc = dmG/dma;

					for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
					{
						Pred[i] = dma*(dmc-(dmc-1)*exp(-1.0*pow((dmb*gpdXi[i]), dmd)));
						if(gsExpoVars.iLogNormal == 1)
							Pred[i] = Slog(Pred[i]);

						dxlk+= 0.5*(double)gpiNi[i]*(p[(int)eAlpha]+p[(int)eRho]*Slog(Pred[i])) +0.5*(((double)gpiNi[i]-1)*gpdYd[i]
						+ (double)gpiNi[i]*(gpdYm[i]-Pred[i])*(gpdYm[i]-Pred[i]))/exp(p[(int)eAlpha]+p[(int)eRho]*Slog(Pred[i]));
					}

					dlltemp=-dxlk;
					if (dlltemp > dll)
					{
						dll = dlltemp;
						p[(int)ea] = dma;			// a
						p[(int)ec] = dmc;			// c
						p[(int)eb] = dmb;			// b
						p[(int)ed] = dmd;			// d
					}
#ifdef DO_LOG
					if (giDo_Log)	// Print values to log for investigation
					{
						fprintf(fp_log,"\ndma=%g, dmult=%g, dmb=%g, dmc=%g, dxlk=%g, dll=%g", dma, dmult, dmb, dmc, dxlk, dll);
					}
#endif
				}
			}
			else
			{
#ifdef DO_LOG
				if (giDo_Log)	// Print values to log for investigation
				{
					fprintf(fp_log,"\nInside DoInitValuesM4() Line 1859 loop (gsExpoVars.iSign != 1)");
				}
#endif
				for (j=1; j <= 10; j++)
				{
					dxlk = 0.0;
					dmult = (j == 1 ? 1.05 : (j == 2 ? 2 : (j == 3 ? 5 : (j == 4 ? 10 : (j == 5 ? 20 : (j == 6 ? 50 : (j == 7 ? 100 : (j == 8 ? 200 : (j == 9 ? 500 : 1000)))))))));

					if(gsExpoVars.iLogNormal == 0)
						dmG = gpdYm[iLowest]/dmult;
					else
						dmG = exp(gpdYm[iLowest])/dmult;

					if(gsExpoVars.iLogNormal == 0)
						dma = 1.05 * gpdYm[iGreatest];
					else
						dma = 1.05 * exp(gpdYm[iGreatest]);

					//Compute for b:
					dGreatest	= 0.0;	//Borrowed as the numerator
					dLowest		= 0.0;	//Borrowed as the denominator

					for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
					{
						if(gsExpoVars.iLogNormal == 0)
							dGreatest += (double)pow((-1.0*(gpdXi[i] * Slog( (dmG - gpdYm[i])/(dmG - dma)))), (1.0/dmd));
						else
							dGreatest += (double)pow((-1.0*(gpdXi[i] * Slog( (dmG - exp(gpdYm[i]))/(dmG - dma)))), (1.0/dmd));

						dLowest += (double)(gpdXi[i] * gpdXi[i]);
					}

					dmb = (double)(dGreatest/dLowest);
					dmc = dmG/dma;

					for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
					{
						Pred[i] = dma*(dmc-(dmc-1)*exp(-1.0*pow((dmb*gpdXi[i]), dmd)));
						if(gsExpoVars.iLogNormal == 1)
							Pred[i] = Slog(Pred[i]);

						dxlk+= 0.5*(double)gpiNi[i]*(p[(int)eAlpha]+p[(int)eRho]*Slog(Pred[i])) +0.5*(((double)gpiNi[i]-1)*gpdYd[i] + (double)gpiNi[i]*(gpdYm[i]-Pred[i])*(gpdYm[i]-Pred[i]))/exp(p[(int)eAlpha]+p[(int)eRho]*Slog(Pred[i]));
					}

					dlltemp=-dxlk;
					if (dlltemp > dll)
					{
						dll = dlltemp;
						p[(int)ea] = dma;			// a
						p[(int)ec] = dmc;			// c
						p[(int)eb] = dmb;			// b
						p[(int)ed] = dmd;			// d
					}
#ifdef DO_LOG
					if (giDo_Log)	// Print values to log for investigation
					{
						fprintf(fp_log,"\ndma=%g, dmult=%g, dmb=%g, dmc=%g, dxlk=%g, dll=%g", dma, dmult, dmb, dmc, dxlk, dll);
					}
#endif
				}
			}
		}
	}

#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\nAFTER CHANGE IN DoInitValuesM4(double p[])\n");
		for(i=(int)eAlpha; i <= NBR_OF_PARMS; i++)
		{
			fprintf(fp_log,"p[%d] = %g,  ", i, p[i]);
		}
	}
#endif
	FREE_DVECTOR (Pred, 1, gsExpoVars.iNbrObs_Total);
}

/**********************************************************
** READ_OBSDATA4V--used to read 4 column data in 4 vectors.
***********************************************************/
int READ_OBSDATA4V(int nObs,double gpdXi[],int gpiNi[],double gpdYm[],double gpdYd[])
{
	int     iNmiss;          /*number of records with missing values*/
	int     i,j,n,m;        /*count and iteration control variables*/
	double  value;          /*temp variable*/

	iNmiss = 0;
	for(i=1;i<=nObs;i++){
		n = i-iNmiss;
		m = 0;
		for (j=1;j<=4;j++)
		{
			fscanf(fp_in,"%lf",&value);
			if (value != MISSING){
				if (j==1)     gpdXi[n]=value;
				if (j==2)     gpiNi[n]=(int)value;
				if (j==3)     gpdYm[n]=value;
				if (j==4)     gpdYd[n]=value;
			}
			else m++;
		}
		if (m != 0)	   iNmiss++;
		else if (gpdXi[n] < 0)      iNmiss++;
	}
	return iNmiss;
}

void Exponential_fit(int nModel)
{
	double *fitparms, *fitparms2, *fitparms3, *beginp, ll, ll2, ll3;
	double *doses, *means, *svar;
	long int *nAnim, *Spec2, *bind, *bind2, *bind3, optite, optite2, optite3;
	long int nvar, signs, nparms, restr;
	long int nresm, model_type, flag, optiteflag, lognormal;
	int gmcnt, n, j, i, k;
	double *dMLEtemp, dDivisor;
	double dLikelihoodtemp = 0.0;

	nvar = gsExpoVars.iNbrObs_Total;
	signs = gsExpoVars.iSign;
	nparms = NBR_OF_PARMS;
	restr = 0;
	nresm = 0;
	model_type = gpimodtype;
	lognormal = (long int)gsExpoVars.iLogNormal;
	ll = 0.0; ll2 = 0.0; ll3 = 0.0;

	fitparms = DVECTOR (0, nparms - 1);
	fitparms2 = DVECTOR (0, nparms - 1);
	fitparms3 = DVECTOR (0, nparms - 1);
	beginp = DVECTOR (0, nparms - 1);

	doses = DVECTOR(0, nvar - 1);
	means = DVECTOR(0, nvar - 1);
	svar =  DVECTOR(0, nvar - 1);
	nAnim = LIVECTOR (0, nvar - 1);

	Spec2 = LIVECTOR (0, nparms - 1);
	bind = LIVECTOR (0, nparms - 1);
	bind2 = LIVECTOR (0, nparms - 1);
	bind3 = LIVECTOR (0, nparms - 1);
	dMLEtemp = DVECTOR(0, nparms - 1);

	for(i=0; i < nparms; i++)
	{
		fitparms[i] = 0.00;
		fitparms2[i] = 0.00;
		fitparms3[i] = 0.00;
		beginp[i] = 0.00;
		Spec2[i] = 0.00;
		bind[i] = 0;
		bind2[i] = 0;
		bind3[i] = 0;
		dMLEtemp[i] = 0.00;
	}
	for(i=0; i < nvar; i++)
	{
		doses[i] = 0.00;
		means[i] = 0.00;
		svar[i] = 0.00;
		nAnim[i] = 0.00;
	}

	optiteflag = 0;

	//HERE WE WANT A LOOP; k = 1 TO 4  starting here and ending at spot shown below
	//1: NO SCALING
	//2: USE SCALED MEANS ONLY
	//3: USE SCALED DOSES ONLY
	//4: USE SCALED MEANS AND DOSES
	for (k = 1; k < 5; k++)
	{
		if((k != 4) && (lognormal == 0))
			continue;

		if((k != 3) && (lognormal == 1))
			continue;

		optite = optite2 = optite3 = -5;
		ll = ll2 = ll3 = -Max_double;
		for(n = 0; n < nparms; n++)
		{
			beginp[n] =  gppdInitPara[nModel][n+1];
			Spec2[n] = gppiSpecPara[nModel][n+1];
#ifdef DO_LOG
			if (giDo_Log)	// Print values to log for investigation
			{
				if(n == 0)
					fprintf(fp_log,"\n\ngppdInitPara[%d][j], for Model #%d\n",nModel,nModel+1);
				fprintf(fp_log,"gppdInitPara[%d][%d] = %g\n", nModel, n+1, gppdInitPara[nModel][n+1]);
			}
#endif
		}
		// START OF STUFF ADDED BY BRUCE 8/11/08 -- this first block needs to be formatted correctly.
		if (k == 2)
		{
			beginp[0] = beginp[0] - Slog(gdmaxYm*gdmaxYm) + beginp[1]*Slog(gdmaxYm);
			beginp[2] = beginp[2]/gdmaxYm;
		}
		if (k == 3)
			beginp[3] = beginp[3]*gdxmax;

		if (k == 4)
		{
			beginp[0] = beginp[0] - Slog(gdmaxYm*gdmaxYm) + beginp[1]*Slog(gdmaxYm);
			beginp[2] = beginp[2]/gdmaxYm;
			beginp[3] = beginp[3]*gdxmax;
		}

		// MORE STUFF ADDED BY BRUCE

		for(n=0; n < gsExpoVars.iNbrObs_Total; n++)
		{
			nAnim[n] = gpiNi[n+1];

			if (k == 3 || k == 4) 					//BCA 8/11/08
				doses[n] = gpdScXi[n+1]; 
			else 
				doses[n] = gpdXi[n+1];

			if (k == 2 || k == 4)					//BCA 8/11/08		
				means[n] = gpdScYm[n+1];	
			else
				means[n] = gpdYm[n+1];

			if (k == 2 || k == 4)					//BCA 8/11/08		
				svar[n] = gpdScYd[n+1];	
			else
				svar[n] = gpdYd[n+1];
		}

#ifdef DO_LOG
		if (giDo_Log == true)
		{
			fprintf(fp_log,"\n*******Spec Init Values, Model #%d, k=%d **********\n", gsExpoVars.iSelect, k);
			for(n = 0; n < nparms; n++)
				fprintf(fp_log,"Spec2[%d] = %ld\n", n, Spec2[n]);

			fprintf(fp_log,"\n*******nAnim Values, Model #%d, k=%d **********\n", gsExpoVars.iSelect, k);
			for(n = 0; n < gsExpoVars.iNbrObs_Total; n++)
				fprintf(fp_log,"nAnim[%d] = %ld\n", n, nAnim[n]);

			fprintf(fp_log,"\n*******doses Values, Model #%d, k=%d **********\n", gsExpoVars.iSelect, k);
			for(n = 0; n < gsExpoVars.iNbrObs_Total; n++)
				fprintf(fp_log,"doses[%d] = %g\n", n, doses[n]);

			fprintf(fp_log,"\n*******means Values, Model #%d, k=%d **********\n", gsExpoVars.iSelect, k);
			for(n = 0; n < gsExpoVars.iNbrObs_Total; n++)
				fprintf(fp_log,"means[%d] = %g\n", n, means[n]);

			fprintf(fp_log,"\n*******svar Values, Model #%d, k=%d **********\n", gsExpoVars.iSelect, k);
			for(n = 0; n < gsExpoVars.iNbrObs_Total; n++)
				fprintf(fp_log,"svar[%d] = %g\n", n, svar[n]);
			fflush(fp_log);
		}
#endif

		gmcnt = 1;
		flag = -1;

		/* This is the first call to getmle.  The parameters will be scaled */
		/* internally by donlp2 by their starting values.                   */
		while((optite2 < 0) || (optite2 > 2))
		{
			if(gmcnt > 1 && gmcnt < 6)
				GetNewParms(beginp, nparms, nModel);

			if(gmcnt > 5 && gmcnt < 10)
			{
#ifdef DO_LOG
				if (giDo_Log == true)
				{
					fprintf(fp_log,"\n\nbeginp before GetMoreParms. Call #%d\n", gmcnt);
					for(i = 0; i < nparms; i++)
					{
						fprintf(fp_log,"beginp[%d] = %12.5g  (%s)\n", i, beginp[i], gaParm_Name[i]);
					}
				}
#endif

				GetMoreParms(beginp, nparms, nModel);

#ifdef DO_LOG
				if (giDo_Log == true)
				{
					fprintf(fp_log,"\n\nbeginp after GetMoreParms. Call #%d\n", gmcnt);

					for(i = 0; i < nparms; i++)
					{
						fprintf(fp_log,"beginp[%d] = %12.5g  (%s)\n", i, beginp[i], gaParm_Name[i]);
					}
					fflush(fp_log);
				}
#endif
			}

			if(gmcnt > 9 && gmcnt < 14)
				GetMLEParms(beginp, nparms, nModel);

#ifdef DO_LOG
			if (giDo_Log == true)
			{
				fprintf(fp_log,"\n***********Call #%d to getmle. k=%d, Model #%d optite2****************\n",gmcnt, k, gsExpoVars.iSelect);
				fprintf(fp_log,"(this call scales the parameters by their starting values in donlp2)\n");
				fprintf(fp_log,"nparms = %ld\n",nparms);
				fprintf(fp_log,"flag = %ld\n",flag);
				fprintf(fp_log,"Model = %d\n",nModel+1);
				fprintf(fp_log,"bmrtype = %d\n",gsExpoVars.iBmr_Type);
				fprintf(fp_log,"These are the parameters going into getmle.  The slope\n");
				fprintf(fp_log,"has been scaled by maxdose^power.\n");

				for(i = 0; i < nparms; i++)
				{
					fprintf(fp_log,"beginp[%d] = %12.5g  (%s)\n", i, beginp[i], gaParm_Name[i]);
				}
				fprintf(fp_log,"************************************************\n");
				fflush(fp_log);
			}
#endif

			getmle_(&nvar, doses, means, nAnim, svar, &nparms, 
				beginp, Spec2, beginp, &restr, &signs, 
				fitparms2, &ll2, &optite2, &nresm, bind2, &model_type, &flag, &lognormal);

#ifdef DO_LOG
			if (giDo_Log == true)
			{
				fprintf(fp_log,"\n*******After Call #%d to getmle. k=%d, Model #%d optite2**********\n",gmcnt, k, gsExpoVars.iSelect);
				fprintf(fp_log,"optite2 = %ld    (good optimum if 0<=optite2<=2)\n",optite2);
				fprintf(fp_log,"flag = %ld\n",flag);
				fprintf(fp_log,"ll2 = %10.5g   (likelihood)\n",ll2);
				fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
				fprintf(fp_log,"These are the parameters coming out of the first run\n");
				fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");

				for(i = 0; i < nparms; i++)
				{
					fprintf(fp_log,"fitparms2[%d] = %12.5g  (%s)\n", i, fitparms2[i], gaParm_Name[i]);
				}
				fprintf(fp_log,"*************************************************\n");
				fflush(fp_log);
			}
#endif

			gmcnt++;

			if(gmcnt == 14)
				break;
		}  //End of while((optite2 < 0) || (optite2 > 2)), first run.

		flag = 0;
		/* This is at most the 14th call to getmle.  The initial parameters will be */
		/* the fit parameters from the scaled (previous) call                       */
#ifdef DO_LOG
		if (giDo_Log == true)
		{
			fprintf(fp_log,"\n***********Call #%d to getmle. k=%d, Model #%d optite3****************\n",gmcnt, k, gsExpoVars.iSelect);
			fprintf(fp_log,"(This call doesn't scale the parameters by their starting values in\n");
			fprintf(fp_log," donlp2.  This call is using the parms coming out of the previous calls.)\n");
			fprintf(fp_log,"nparms = %ld\n",nparms);
			fprintf(fp_log,"flag = %ld\n",flag);
			fprintf(fp_log,"Model = %d\n",nModel+1);
			fprintf(fp_log,"bmrtype = %d\n",gsExpoVars.iBmr_Type);
			fprintf(fp_log,"These are the parameters going into getmle.  The slope\n");
			fprintf(fp_log,"has been scaled by maxdose^power.\n");

			for(i = 0; i < nparms; i++)
			{
				fprintf(fp_log,"fitparms2[%d] = %12.5g  (%s)\n", i, fitparms2[i], gaParm_Name[i]);
			}
			fprintf(fp_log,"************************************************\n");
			fflush(fp_log);
		}
#endif

		getmle_(&nvar, doses, means, nAnim, svar, &nparms, 
			fitparms2, Spec2, beginp, &restr, &signs, 
			fitparms3, &ll3, &optite3, &nresm, bind3, &model_type, &flag, &lognormal);

#ifdef DO_LOG
		if (giDo_Log == true)
		{
			fprintf(fp_log,"\n*******After Call #%d to getmle. k=%d, Model #%d optite3****************\n",gmcnt, k, gsExpoVars.iSelect);
			fprintf(fp_log,"optite3 = %ld    (good optimum if 0<=optite3<=2)\n",optite3);
			fprintf(fp_log,"flag = %ld\n",flag);
			fprintf(fp_log,"ll3 = %10.5g   (likelihood)\n",ll3);
			fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
			fprintf(fp_log,"These are the parameters coming out of the first run\n");
			fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");

			for(i = 0; i < nparms; i++)
			{
				fprintf(fp_log,"fitparms3[%d] = %12.5g  (%s)\n", i, fitparms3[i], gaParm_Name[i]);
			}
			fprintf(fp_log,"*************************************************\n");
			fflush(fp_log);
		}
#endif
		gmcnt++;

		/* This starts the loops without scaling in donlp2 */
		n = -1;
		for(flag=0; flag < 2; flag++)
		{
			while((optite < 0) || (optite > 2))
			{
				if(n < 30 && n > -1)
				{
					if (optite != 3)
						GetNewParms(beginp, nparms, nModel);	/* Get a new starting point */
					else
					{
						for (j = 0; j < nparms; j++)
							beginp[j] = fitparms[j];
					}
				}

				if(n > 29 && n < 60)
				{
					if (optite != 3)
						GetMoreParms(beginp, nparms, nModel);		/* Get a new starting point */
					else
					{
						for (j = 0; j < nparms; j++)
							beginp[j] = fitparms[j];
					}
				}

#ifdef DO_LOG
				if (giDo_Log == true)
				{
					fprintf(fp_log,"\n***********Call #%d to getmle. k=%d, Model #%d optite****************\n",gmcnt, k, gsExpoVars.iSelect);
					fprintf(fp_log,"optite = %ld    (good optimum if 0<=optite<=2)\n",optite);
					fprintf(fp_log,"flag = %ld,  signs = %ld\n",flag, signs);
					fprintf(fp_log,"ll = %10.5g   (likelihood)\n",ll);
					fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
					fprintf(fp_log,"These are the parameters coming out of the first run\n");
					fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");

					for(i = 0; i < nparms; i++)
					{
						fprintf(fp_log,"beginp[%d] = %12.5g  (%s)\n", i, beginp[i], gaParm_Name[i]);
					}
					fprintf(fp_log,"*************************************************\n");
					fflush(fp_log);
				}
#endif

				getmle_(&nvar, doses, means, nAnim, svar, &nparms, 
					beginp, Spec2, beginp, &restr, &signs, 
					fitparms, &ll, &optite, &nresm, bind, &model_type, &flag, &lognormal);

#ifdef DO_LOG
				if (giDo_Log == true)
				{
					fprintf(fp_log,"\n*******After Call #%d to getmle. k=%d, Model #%d optite****************\n",gmcnt, k, gsExpoVars.iSelect);
					fprintf(fp_log,"optite = %ld    (good optimum if 0<=optite<=2)\n",optite);
					fprintf(fp_log,"flag = %ld,  signs = %ld\n",flag, signs);
					fprintf(fp_log,"ll = %10.5g   (likelihood)\n",ll);
					fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
					fprintf(fp_log,"These are the parameters coming out of the first run\n");
					fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");

					for(i = 0; i < nparms; i++)
					{
						fprintf(fp_log,"fitparms[%d] = %12.5g  (%s)\n", i, fitparms[i], gaParm_Name[i]);
					}
					fprintf(fp_log,"*************************************************\n");
					fflush(fp_log);
				}
#endif

				//#ifdef DO_LOG
				//if (giDo_Log == true)
				//{
				//}
				//#endif
				gmcnt++;
				n++;
#ifdef DO_LOG
				if (giDo_Log == true)
				{
					fprintf(fp_log,"\n*******Ln 2040. Model #%d  after optite ****************\n", gsExpoVars.iSelect);
					fflush(fp_log);
				}
#endif

				if(n > 60)
				{
					n = -1;
					break;
				}

			} //End of while((optite < 0) || (optite > 2)) 3rd run.
#ifdef DO_LOG
			if (giDo_Log == true)
			{
				fprintf(fp_log,"\n*******Ln 2055. End of while((optite < 0) || (optite > 2))************\n");
				fflush(fp_log);
			}
#endif

			if ((optite >= 0) && (optite <= 2))
			{
				flag = 2;
				break;
			}

		}	/* end  for (flag=0; flag<2; flag++) */

		/* This decides if the scaling model is better than the unscaled model */
		/* or not. */
#ifdef DO_LOG
		if (giDo_Log == true)
		{
			fprintf(fp_log, "\n\nOptite Values:\n");
			fprintf(fp_log, "\noptite = %ld", optite);
			fprintf(fp_log, "\noptite2 = %ld", optite2);
			fprintf(fp_log, "\noptite3 = %ld", optite3);
		}
#endif

		if ((optite2 >= 0) && (optite2 <= 2))
		{
#ifdef DO_LOG
			if (giDo_Log == true)
			{
				fprintf(fp_log, "\nINSIDE if ((optite2 >= 0) && (optite2 <= 2))\n");
			}
#endif
			if ((ll2 > ll) || ((optite < 0) || (optite > 2)))
			{
#ifdef DO_LOG
				if (giDo_Log == true)
				{
					fprintf(fp_log, "\nINSIDE if ((ll2 > ll) || ((optite < 0) || (optite > 2)))\n");
				}
#endif
				for (j = 0; j < nparms; j++)
				{
					fitparms[j] = fitparms2[j];
					bind[j] = bind2[j];

#ifdef DO_LOG
					if (giDo_Log == true)
					{
						//fprintf(fp_log, "fitparms2[%d] (%12.5g)\n", nModel, j+1, j, fitparms2[j]);
						fprintf(fp_log, "fitparms2[%d] (%12.5g)\n", j, fitparms2[j]);
						fflush(fp_log);
					}
#endif
				}
				optite = optite2;
				ll = ll2;
#ifdef DO_LOG
				if (giDo_Log == true)
				{
					fprintf(fp_log, "likelihood (ll2) = %12.5g\n", ll);
					fflush(fp_log);
				}
#endif
			}
		}

		if ((optite3 >= 0) && (optite3 <= 2))
		{
#ifdef DO_LOG
			if (giDo_Log == true)
			{
				fprintf(fp_log, "\nINSIDE if ((optite3 >= 0) && (optite3 <= 2))\n");
			}
#endif
			if ((ll3 > ll) || ((optite < 0) || (optite > 2)))
			{
#ifdef DO_LOG
				if (giDo_Log == true)
				{
					fprintf(fp_log, "\nINSIDE if ((ll3 > ll) || ((optite < 0) || (optite > 2)))\n");
				}
#endif
				for (j = 0; j < nparms; j++)
				{
					fitparms[j] = fitparms3[j];
					bind[j] = bind3[j];

#ifdef DO_LOG
					if (giDo_Log == true)
					{
						fprintf(fp_log, "fitparms3[%d] (%12.5g)\n", j, fitparms3[j]);
						fflush(fp_log);
					}
#endif
				}
				optite = optite3;
				ll = ll3;

#ifdef DO_LOG
				if (giDo_Log == true)
				{
					fprintf(fp_log, "\nlikelihood (ll3) = %12.5g\n", ll);
					fflush(fp_log);
				}
#endif
			}
		}

		/* Let user know if no optimum was found */
		if ((optite < 0) || (optite > 2))
		{
			fprintf(fp_out, "\n\n!!! Warning:  optimum may not have been found for Model %d          !!!", nModel+1);
			fprintf(fp_out, "\n!!! Bad completion code in maximum likelihood optimization routine  !!!");
			fprintf(fp_out, "\n!!! Try choosing different initial values                           !!!\n\n");

#ifdef DO_LOG
			if (giDo_Log == true)
			{
				fprintf(fp_log, "\n\nWARNING: Optimum may not have been found for Model %d\n\n", nModel + 1);
				fflush(fp_log);
			}
#endif
		}

		if (optite >= 0 && optite <=2)
		{
#ifdef DO_LOG
			if (giDo_Log == true)
			{
				fprintf(fp_log, "\n\nINSIDE if (optite >= 0 && optite <=2)\n");
			}
#endif
			optiteflag = 1;
			if (k == 1)
			{
				for (j = 0; j < nparms; j++)
				{
					dMLEtemp[j] = fitparms[j];
				}

			}
			else if(k == 2)
			{
				for (j = 0; j < nparms; j++)
				{
					dMLEtemp[j] = fitparms[j];
				}
				dMLEtemp[0] = fitparms[0] + Slog(gdmaxYm*gdmaxYm) - fitparms[1]*Slog(gdmaxYm);
				dMLEtemp[2] = fitparms[2]*gdmaxYm;


			}
			else if (k == 3)
			{
				for (j = 0; j < nparms; j++)
				{
					dMLEtemp[j] = fitparms[j];
				}
				dMLEtemp[3] = fitparms[3]/gdxmax;


			}
			else if (k == 4)
			{
				for (j = 0; j < nparms; j++)
				{
					dMLEtemp[j] = fitparms[j];
				}
				dMLEtemp[0] = fitparms[0] + Slog(gdmaxYm*gdmaxYm) - fitparms[1]*Slog(gdmaxYm);
				dMLEtemp[2] = fitparms[2]*gdmaxYm;
				dMLEtemp[3] = fitparms[3]/gdxmax;
#ifdef DO_LOG
				if (giDo_Log)	// Print values to log for investigation
				{
					fprintf(fp_log,"\n\nLine 2518,else if (k == 4), Model %d: gdxmax=%g, gdmaxYm=%g", nModel+1, gdxmax, gdmaxYm);
					for(j = 0; j < nparms; j++)
					{
						fprintf(fp_log,"\ndMLEtemp[%d]=%g: ", j, dMLEtemp[j]);
					}
				}
#endif
			}

			for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
			{
				if (nModel == 3 || nModel == 4)
				{
					gppdPredict[nModel][i] = dMLEtemp[2]*(dMLEtemp[4]-(dMLEtemp[4]-1)*exp(-1.0*pow((dMLEtemp[3]*gpdXi[i]),dMLEtemp[5])));
					if(gsExpoVars.iLogNormal == 1)
						gppdPredict[nModel][i] = Slog(gppdPredict[nModel][i]);

					dDivisor = exp(dMLEtemp[0]+dMLEtemp[1]*Slog(gppdPredict[nModel][i]));
					if(dDivisor == 0)
					{
#ifdef DO_LOG
						if (giDo_Log)	// Print values to log for investigation
						{
							fprintf(fp_log,"\n\n,Line 2542 INSIDE if(dDivisor == 0), giRun[%d]=%d: ", nModel, giRun[nModel]);
						}
#endif
						giRun[nModel] = 0;
						return;
					}
					dLikelihoodtemp += 0.5*(double)gpiNi[i]*(dMLEtemp[0]+dMLEtemp[1]*Slog(gppdPredict[nModel][i]))  
						+ 0.5*(((double)gpiNi[i]-1)*gpdYd[i] + (double)gpiNi[i]*(gpdYm[i]-gppdPredict[nModel][i])*(gpdYm[i]-gppdPredict[nModel][i]))/
						dDivisor;
				}
				else
				{
					gppdPredict[nModel][i] = dMLEtemp[2]*exp(gsExpoVars.iSign*pow((dMLEtemp[3]*gpdXi[i]),dMLEtemp[5]));
					if(gsExpoVars.iLogNormal == 1)
						gppdPredict[nModel][i] = Slog(gppdPredict[nModel][i]);


					dDivisor = exp(dMLEtemp[0]+dMLEtemp[1]*Slog(gppdPredict[nModel][i]));
					if(dDivisor == 0)
					{
#ifdef DO_LOG
						if (giDo_Log)	// Print values to log for investigation
						{
							fprintf(fp_log,"\n\n,Line 2566 INSIDE if(dDivisor == 0), giRun[%d]=%d: ", nModel, giRun[nModel]);
						}
#endif
						giRun[nModel] = 0;
						return;
					}
					dLikelihoodtemp+= 0.5*(double)gpiNi[i]*(dMLEtemp[0]+dMLEtemp[1]*Slog(gppdPredict[nModel][i])) +0.5*(((double)gpiNi[i]-1)*gpdYd[i]
					+ (double)gpiNi[i]*(gpdYm[i]-gppdPredict[nModel][i])*(gpdYm[i]-gppdPredict[nModel][i]))/dDivisor;
				}
			}
			dLikelihoodtemp = -dLikelihoodtemp; 
#ifdef DO_LOG
			if (giDo_Log)	// Print values to log for investigation
			{
				fprintf(fp_log,"\n\nLine 2364, if (nModel == 3 || nModel == 4): dLikelihoodtemp=%g, gpdLikelihoods[%d]=%g", dLikelihoodtemp, nModel, gpdLikelihoods[nModel]);
				for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
				{
					fprintf(fp_log,"\ngppdPredict[%d][%d]=%g", nModel, i, gppdPredict[nModel][i]);
				}
			}
#endif

			if (dLikelihoodtemp > gpdLikelihoods[nModel])
			{
				gpdLikelihoods[nModel] = dLikelihoodtemp;
				for (j = 1; j <= nparms; j++)
				{
					gppdMLEs[nModel][j] = dMLEtemp[j-1];
					
				}
			}
		}

		for (j = 1; j <= nparms; j++)
		{
			bounded[nModel][j] = bind[j-1]; /* put back into the program */
			//if(bounded[nModel][j] != 1 && bounded[nModel][j] != 0)
			//{
			//	bounded[nModel][j] = 1; /* Hopefully, we never get here; originally, this was 0, but 1 looks
			//					like a better warning flag */
			//}
		}
	
	} //End of for(k = 1; k < 5; k++)

	// Here at the end of the loop over k if optiteflag=1 then there was a run that converged
	// If optiteflag = 0, no convergence, and we want to print a warning (take out printing of warnings above)
	// if optiteflag = 1, then we want to print MLEs, maxLL, calc and print BMD and go to BMDL calculation
	//  This conditional processing needs to be added still (BCA 8/14/08)
	if(optiteflag==1)
	{
#ifdef DO_LOG
		if (giDo_Log)	// Print values to log for investigation
		{
			fprintf(fp_log,"\n\nLine 2357 INSIDE if(optiteflag==1), giRun[%d]=%d: ", nModel, giRun[nModel]);
		}
#endif
		giRun[nModel]=1;
	}

	if (gsExpoVars.iBmr_Type == 0) // Absolute difference in means only need LogNorm if here for BMRtype = 1

		gdBMD[nModel] = gppdMLEs[nModel][(int)ea] + gdAdverseBMDEffect;

	else if (gsExpoVars.iBmr_Type == 1) // Standard deviation
	{
		gdBMD[nModel] = gppdMLEs[nModel][(int)ea] + gdAdverseBMDEffect*sqrt(exp(gppdMLEs[nModel][1]+gppdMLEs[nModel][(int)eRho]*log(gppdMLEs[nModel][(int)ea])));
		if(gsExpoVars.iLogNormal == 1)
			gdBMD[nModel] = exp(log(gppdMLEs[nModel][(int)ea]) + gdAdverseBMDEffect*sqrt(exp(gppdMLEs[nModel][1])));
	}
	else if (gsExpoVars.iBmr_Type == 2) // Relative deviation

		gdBMD[nModel] = gppdMLEs[nModel][(int)ea] + gdAdverseBMDEffect*gppdMLEs[nModel][(int)ea];

	else if (gsExpoVars.iBmr_Type == 3) // Point

		gdBMD[nModel] = gsExpoVars.dBmdEffect;

	else if (gsExpoVars.iBmr_Type == 4) // Extra
	{
		gdBMD[nModel] = gppdMLEs[nModel][(int)ea] + gsExpoVars.dBmdEffect*(gppdMLEs[nModel][(int)ea]*gppdMLEs[nModel][(int)ec] - gppdMLEs[nModel][(int)ea]);
#ifdef DO_LOG
		if (giDo_Log)	// Print values to log for investigation
		{
			fprintf(fp_log,"\n\nLine 2436: gdBMD[%d]=%g, gsExpoVars.dBmdEffect=%g", nModel, gdBMD[nModel], gsExpoVars.dBmdEffect);
			for (j = 1; j <= nparms; j++)
			{
				fprintf(fp_log,"\ngppdMLEs[%d][%d]=%g: ", nModel, j, gppdMLEs[nModel][j]);
			}
		}
#endif
	}

#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\nLine 2448: gdBMD[%d]=%g, giBmdRun[%d]=%d, giRun[%d]=%d", nModel, gdBMD[nModel], nModel, giBmdRun[nModel], nModel, giRun[nModel]);
	}
#endif

	if (gdBMD[nModel] < 0) // no Models can have negative responses as BMR

		giBmdRun[nModel] = 0;

	else if ((nModel < 3) && (gsExpoVars.iBmr_Type == 4)) // Models 2 and 3 cannot take BMR type 4 (extra risk)

		giBmdRun[nModel] = 0;

	else if (gsExpoVars.iSign*(gdBMD[nModel] - gppdMLEs[nModel][(int)ea]) < 0) // BMR not in right direction

		giBmdRun[nModel] = 0;

	else if ((nModel > 2) && ((gdBMD[nModel] - gppdMLEs[nModel][(int)ea])/(gppdMLEs[nModel][(int)ea]*gppdMLEs[nModel][(int)ec] - gppdMLEs[nModel][(int)ea])) >= 1.0)  // BMR outside range from a to a*c 

		giBmdRun[nModel] = 0;

	else if (nModel < 3)  // everything looks good for computing BMD for models 2 and 3

		gdBMD[nModel] = getBMD23(gdBMD[nModel], gppdMLEs[nModel][(int)ea], gppdMLEs[nModel][(int)eb], gppdMLEs[nModel][(int)ed]);

	else if (nModel > 2)  // everything looks good for computing BMD for models 4 and 5

		gdBMD[nModel] = getBMD45(gdBMD[nModel], gppdMLEs[nModel][(int)ea], gppdMLEs[nModel][(int)eb], gppdMLEs[nModel][(int)ec], gppdMLEs[nModel][(int)ed]);

#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\nLine 2479: gdBMD[%d]=%g, giBmdRun[%d]=%d, giRun[%d]=%d\n\n", nModel, gdBMD[nModel], nModel, giBmdRun[nModel], nModel, giRun[nModel]);
	}
#endif

	if(giBmdRun[nModel]==1)  //Compute BMDL
		gdBMDL[nModel] = BMDL_func(nModel);

	FREE_DVECTOR (dMLEtemp, 0, nparms - 1);

	FREE_LIVECTOR (bind3, 0, nparms - 1);
	FREE_LIVECTOR (bind2, 0, nparms - 1);
	FREE_LIVECTOR (bind, 0, nparms - 1);
	FREE_LIVECTOR (Spec2, 0, nparms - 1);

	FREE_LIVECTOR (nAnim, 0, nvar - 1);
	FREE_DVECTOR (svar, 0, nvar - 1);
	FREE_DVECTOR (means, 0, nvar - 1);
	FREE_DVECTOR (doses, 0, nvar - 1);

	FREE_DVECTOR (beginp, 0, nparms - 1);
	FREE_DVECTOR (fitparms3, 0, nparms - 1);
	FREE_DVECTOR (fitparms2, 0, nparms - 1);
	FREE_DVECTOR (fitparms, 0, nparms - 1);
}

double getBMD23(double dBMD, double da, double db, double dd)
{
	double dBMD23 = gsExpoVars.iSign*log(dBMD/da);
	dBMD23 = pow(dBMD23, (1/dd));
	return (dBMD23/db);
}

double getBMD45(double dBMD, double da, double db, double dc, double dd)
{
	double dBMD45 = (dBMD - da*dc)/(da - da*dc);
	dBMD45 = pow((-log(dBMD45)), (1/dd));
	return (dBMD45/db);
}

double BMDL_func(int nModel)
{
	double *fitparms, *beginp, *doses, *means, *svar, *parms;
	double bmdl, target, dLR, dBMR, dDose;
	long int *nAnim, *Spec2, *bind;
	long int nvar, signs, nparms, restr, bmrtype, which;
	long int nresm, model_type, flag, optite, lognormal;
	int gccnt, j, i;

	/* Initialize local variables */
	bmdl = target = dLR = dBMR = 0;
	model_type = (long int) gpimodtype;
	dBMR = gsExpoVars.dBmdEffect;
	lognormal = (long int)gsExpoVars.iLogNormal;

	nvar = (long int)gsExpoVars.iNbrObs_Total;
	nparms = (long int)NBR_OF_PARMS;
	flag = restr = nresm = 0;
	optite = -5;
	signs = (long int)gsExpoVars.iSign;
	bmrtype = (long int)gsExpoVars.iBmr_Type;
	which = 1; //Want a lower confidense limit

	fitparms = DVECTOR (0, nparms - 1);
	beginp = DVECTOR (0, nparms - 1);
	parms = DVECTOR(0, nparms - 1);
	Spec2 = LIVECTOR (0, nparms - 1);
	bind = LIVECTOR (0, nparms - 1);

	doses = DVECTOR(0, nvar - 1);
	means = DVECTOR(0, nvar - 1);
	svar =  DVECTOR(0, nvar - 1);
	nAnim = LIVECTOR (0, nvar - 1);
	bmdl = gdBMD[nModel];

	for(i=0; i < nparms; i++)
	{
		fitparms[i] = 0.00;
		beginp[i] = 0.00;
		parms[i] = 0.00;
		Spec2[i] = 0;
		bind[i] = 0;
	}
	for(i=0; i < nparms; i++)
	{
		beginp[i] = gppdMLEs[nModel][i+1];
		parms[i] = gppdSpecPara[nModel][i+1];
		Spec2[i] = gppiSpecPara[nModel][i+1];
	}
	for(i=0; i < nvar; i++)
	{
		doses[i] = gpdXi[i+1];
		means[i] = gpdYm[i+1];
		svar[i] = gpdYd[i+1];
		nAnim[i] = gpiNi[i+1];
	}
	/* End of Initialize local variables */

	if(gsExpoVars.dBmdConfi_Level < 0.5)
		dLR = 0.5*QCHISQ(1.0 - 2.0 * gsExpoVars.dBmdConfi_Level, 1);
	else
		dLR = 0.5*QCHISQ(2.0 * gsExpoVars.dBmdConfi_Level - 1.0, 1);
	target = gpdLikelihoods[nModel] - dLR;

	// scaling by max dose:
	for(i=0; i < nvar; i++)
	{
		doses[i] = doses[i] / gdxmax;
	}
	/*Need help on computation of dDose from Power.c, Bruce verify if correct.
	dDose also known as "xb", "thedose" in Pow_BMD() and "Dose" in BMDL_func(), 
	computed starting from line 3735 to line 3865 and passed to 
	BMDL_func(nparm, BMD_lk, xb, p, tol)
	*/
	dDose = gdBMD[nModel]/gdxmax;

	//rescaling parameter b (index 3 when starting from 0) because of scaling of dose:
	parms[3] = parms[3]*gdxmax;
	beginp[3] = beginp[3]*gdxmax;

	for(gccnt=1; gccnt < 34; gccnt++)
	{
		/*Next 10, 2 to 11*/
		if((gccnt > 1 && gccnt < 12))
		{
			for(i=0; i < nparms; i++)
				beginp[i] = fitparms[i];
		}// End if(gccnt > 1 && gccnt < 12), 2 to 11

		if(gccnt == 12 || gccnt == 28)
		{
			for (j = 1; j <= nparms; j++)
				beginp[j - 1] = gppdMLEs[nModel][j];
			beginp[3] = beginp[3]*gdxmax;
		}

		if((gccnt > 12 && gccnt < 18)) //13 to 17
		{
			if(optite != 4)
				GetNewParms(beginp, nparms, nModel);
			else
			{
				for(i=0; i < nparms; i++)
					beginp[i] = fitparms[i];
			}
		}// End if(gccnt > 12 && gccnt < 18), 13 to 17

		if(gccnt != 28 && gccnt > 17) //18 to 27 && 28 to 38
		{
			if(optite != 4)
				GetMoreParms(beginp, nparms, nModel);
			else
			{
				for(i=0; i < nparms; i++)
					beginp[i] = fitparms[i];
			}
		}// End if(gccnt != 28 && gccnt > 17), 18 to 27 && 28 to 38

#ifdef DO_LOG
		if (giDo_Log == true)
		{
			fprintf(fp_log,"\n***********Call #%d to getcl.****************\n",gccnt);
			fprintf(fp_log,"(this call scales the parameters by their starting values in donlp2)\n");
			fprintf(fp_log,"nparms = %ld\n",nparms);
			fprintf(fp_log,"flag = %ld\n",flag);
			fprintf(fp_log,"BMR = %10.5g\n",dBMR);
			fprintf(fp_log,"bmrtype = %ld\n",bmrtype);
			fprintf(fp_log,"target = %10.5g\n",target);
			fprintf(fp_log,"bmdl = %10.5g  (This should be the BMD)\n",bmdl);
			fprintf(fp_log,"These are the parameters going into getcl.  The slope\n");
			fprintf(fp_log,"has been scaled by maxdose^power.\n");
			for(i = 0; i < nparms; i++)
			{
				fprintf(fp_log,"beginp[%d] = %12.5g  (%s)\n", i, beginp[i], gaParm_Name[i]);
			}
			fprintf(fp_log,"************************************************\n");
			fflush(fp_log);
		}
#endif

		//if(gccnt==1 || gccnt==17)
		//{
		getcl_(&which, &nvar, doses, means, nAnim, svar, &nparms, &dBMR,
			&dDose, &target, beginp, Spec2, parms, &bmrtype, &restr,
			&bmdl, fitparms, &optite, &nresm, bind, &signs, &model_type, &flag, &lognormal);
		//}
		//else
		//{
		//	getcl_(&which, &nvar, doses, means, nAnim, svar, &nparms, &dBMR,
		//		&dDose, &target, beginp, Spec2, beginp, &bmrtype, &restr,
		//		&bmdl, fitparms, &optite, &nresm, bind, &signs, &model_type, &flag);
		//}// End if(gccnt==1 || gccnt==17)

#ifdef DO_LOG
		if (giDo_Log == true)
		{
			fprintf(fp_log,"\n*******After Call #%d to getcl.**********\n",gccnt);
			fprintf(fp_log,"nparms = %ld\n",nparms);
			fprintf(fp_log,"flag = %ld\n",flag);
			fprintf(fp_log,"bmdl = %10.5g   (BMDL)\n",bmdl);
			fprintf(fp_log,"optite = %ld    (good optimum if 0<=optite<=2)\n",optite);
			fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
			if(gccnt > 17)
				fprintf(fp_log,"These are the parameters coming out of the second run\n");
			else
				fprintf(fp_log,"These are the parameters coming out of the first run\n");
			fprintf(fp_log,"of getcl.  The slope is still in scaled form.\n");
			for(i = 0; i < nparms; i++)
			{
				fprintf(fp_log,"fitparms[%d] = %12.5g  (%s)\n", i, fitparms[i], gaParm_Name[i]);
			}
			fprintf(fp_log,"************************************************\n");
			fflush(fp_log);
		}
#endif
		if(optite>=0 && optite <=3)
			break;

	}//End for(gccnt=0; gccnt < 35; gccnt++)
	/* Let user know if no optimum was found */
	if ((optite < 0) || (optite > 3))
	{
		fprintf(fp_out, "Warning:  optimum may not have been found from model %d (BMDL).  Bad process completion in Optimization routine.\n", nModel+1);
		bmdl = -1;
	}
	else
	{
		bmdl = bmdl * gdxmax;
	}

	FREE_LIVECTOR (nAnim, 0, nvar - 1);
	FREE_DVECTOR (svar, 0, nvar - 1);
	FREE_DVECTOR (means, 0, nvar - 1);
	FREE_DVECTOR (doses, 0, nvar - 1);

	FREE_LIVECTOR (bind, 0, nparms - 1);
	FREE_LIVECTOR (Spec2, 0, nparms - 1);
	FREE_DVECTOR (parms, 0, nparms - 1);
	FREE_DVECTOR (beginp, 0, nparms - 1);
	FREE_DVECTOR (fitparms, 0, nparms - 1);

	return bmdl;
}

/***********************************************************
*	Given a vector of parameter values, and the number of
*	parameters in that vector, this function will return three
*	new parameter values to restart the optimization if a "bad"
*	completion code is returned from GETCL(), using a uniform
*	random number centered at p[i]
***********************************************************/
void GetMoreParms(double *p, int size, int nModel)
{
	int i;
	/* Find parameters by randomly selecting new parameters in
	/  a uniform interval of p[i] +/- .25*p[i] */
	int *Spec = gppiSpecPara[nModel];

	for (i = 0; i < size; i++)
	{
		if (Spec[i + 1] != 1)
		{
			p[i] = ((p[i]*0.5) * (double)rand() / 32767.0) + p[i] * .75;
		}
	}

	/* If parameters are to be restricted, make sure restrictions
	/  are not violated */
	if(nModel == 3 || nModel == 4) //Model 4 && 5
	{
		if(gsExpoVars.iSign == 1 && p[4] < 1)
			p[4] = 1.00001;
		else if(gsExpoVars.iSign == -1 && p[4] > 1)
			p[4] = .9999;
		else if(gsExpoVars.iSign == -1 && p[4] < 0)
			p[4] = -p[4];
	}
	if(nModel == 2 || nModel == 4)
	{
		if(p[5] <= 1)
			p[5] = 1.00001;
	}
	if(p[3] < 0)
		p[3] = -p[3];
	if(p[2] < 0)
		p[2] = -p[2];

#ifdef DO_LOG
	if (giDo_Log == true)
	{
		fprintf(fp_log,"\n*******New Values after GetMoreParms() Function, Model #%d**********\n", nModel + 1);
		for(i = 0; i < size; i++)
			fprintf(fp_log,"p[%d] = %g\n", i, p[i]);

		fflush(fp_log);
	}
#endif


}	/* end GetMoreParms */

/***********************************************************
*	Given a vector of parameter values, and the number of
*	parameters in that vector, this function will return three
*	new parameter values to restart the optimization if a "bad"
*	completion code is returned from GETCL(), using a uniform
*	random number centered at p[i]
***********************************************************/
void GetMLEParms(double *p, int size, int nModel)
{
	int i;
	/* Find parameters by randomly selecting new parameters in
	/  a uniform interval of p[i] +/- .5*p[i] */
	int *Spec = gppiSpecPara[nModel];

	for (i = 0; i < size; i++)
	{
		if (Spec[i + 1] != 1)
		{
			p[i] = (p[i] * (double)rand() / 32767.0) + p[i] * .5;
		}
	}

	/* If parameters are to be restricted, make sure restrictions
	/  are not violated */
	if(nModel == 3 || nModel == 4) //Model 4 && 5
	{
		if(gsExpoVars.iSign == 1 && p[4] < 1)
			p[4] = 1.00001;
		else if(gsExpoVars.iSign == -1 && p[4] > 1)
			p[4] = .9999;
		else if(gsExpoVars.iSign == -1 && p[4] < 0)
			p[4] = -p[4];
	}
	if(nModel == 2 || nModel == 4)
	{
		if(p[5] <= 1)
			p[5] = 1.00001;
	}
	if(p[3] < 0)
		p[3] = -p[3];
	if(p[2] < 0)
		p[2] = -p[2];

#ifdef DO_LOG
	if (giDo_Log == true)
	{
		fprintf(fp_log,"\n*******New Values after GetMLEParms() Function, Model #%d**********\n", nModel + 1);
		for(i = 0; i < size; i++)
			fprintf(fp_log,"p[%d] = %g\n", i, p[i]);

		fflush(fp_log);
	}
#endif

}	/* end GetMLEParms */

/***********************************************************
*	Given a vector of parameter values, and the number of
*	parameters in that vector, this function will return three
*	new parameter values to restart the optimization if a "bad"
*	completion code is returned from GETCL(), using a uniform
*	random number centered at p[i]
***********************************************************/
void GetNewParms(double *p, int size, int nModel)
{
	int i;
	/* Find parameters by randomly selecting new parameters in
	/  a uniform interval of p[i] +/- .005*p[i] */
	int *Spec = gppiSpecPara[nModel];

	for (i = 0; i < size; i++)
	{
		if (Spec[i + 1] != 1)
		{
			p[i] = ((p[i] * .010) * (double)rand() / 32767.0) + p[i] * .995;
		}
	}

	/* If parameters are to be restricted, make sure restrictions
	/  are not violated */
	if(nModel == 3 || nModel == 4) //Model 4 && 5
	{
		if(gsExpoVars.iSign == 1 && p[4] < 1)
			p[4] = 1.00001;
		else if(gsExpoVars.iSign == -1 && p[4] > 1)
			p[4] = .9999;
		else if(gsExpoVars.iSign == -1 && p[4] < 0)
			p[4] = -p[4];
	}
	if(nModel == 2 || nModel == 4)
	{
		if(p[5] <= 1)
			p[5] = 1.00001;
	}
	if(p[3] < 0)
		p[3] = -p[3];
	if(p[2] < 0)
		p[2] = -p[2];

#ifdef DO_LOG
	if (giDo_Log == true)
	{
		fprintf(fp_log,"\n*******New Values after GetNewParms() Function, Model #%d**********\n", nModel + 1);
		for(i = 0; i < size; i++)
			fprintf(fp_log,"p[%d] = %g\n", i, p[i]);

		fflush(fp_log);
	}
#endif
}	/* end GetNewParms */

double LogLike(int nmodel, double p[])
{
	//int iNtot=0;
	double *Pred;
	Pred = DVECTOR(1, gsExpoVars.iNbrObs_Total);
	//double *p;
	//p = DVECTOR(0, NBR_OF_PARMS-1);
	//double ll=-Max_double;
	double dLikelihoodtemp=0.0;
	int i;
	for(i = 1; i <= NBR_OF_PARMS; i++)
	{
#ifdef DO_LOG
		if (giDo_Log)	// Print values to log for investigation
		{
			fprintf(fp_log,"\n(inside LogLike, model %d)p[%d] = %g\n", nmodel+1, i, p[i]);
		}
#endif
	}
	for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
	{
		if ((nmodel == 1) || (nmodel == 2))    // model numbering from 1 to 4
		{
			Pred[i] = p[3]*exp(gsExpoVars.iSign*pow((p[4]*gpdXi[i]),p[6]));
			if(gsExpoVars.iLogNormal == 1)
				Pred[i] = Slog(Pred[i]);

		}
		else	
		{
			Pred[i] = p[3]*(p[5]-(p[5]-1)*exp(-1.0*pow((p[4]*gpdXi[i]),p[6])));
			if(gsExpoVars.iLogNormal == 1)
				Pred[i] = Slog(Pred[i]);
		}
		dLikelihoodtemp += 0.5*(double)gpiNi[i]*(p[1]+p[2]*Slog(Pred[i])) + 0.5*(((double)gpiNi[i]-1)*gpdYd[i] + 
			(double)gpiNi[i]*(gpdYm[i]-Pred[i])*(gpdYm[i]-Pred[i]))/exp(p[1]+p[2]*Slog(Pred[i]));
	}

#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n(inside LogLike)dLikelihoodtemp = %g\n", -dLikelihoodtemp);
	}
#endif

	return -dLikelihoodtemp;
}

int DoInitParam(int iModelNbr)
{
	int i, j, n;
	//char gaPName[6];
	//if(gsExpoVars.iCons_Var == 0)
	//{
	//	strcpy(gaPName, "l");
	//	strcat(gaPName, gaParm_Name[0]);
	//}
	//else
	//	strcpy(gaPName, gaParm_Name[0]);

	char nc[]="NC";
	char model[24];
	char model1[24];
	char model2[24];
	char model3[24];
	char model4[24]; 
	char specified[] = " *";
	char specblank[] = "  ";
	char specified2[] = " Specified";
	char specblank2[] = "          ";

	if(iModelNbr > 0)
	{
		OUTPUT_TEXT("\n\n                  Initial Parameter Values");
		fprintf(fp_out,"\n                  Variable          Model %d\n", iModelNbr+1);
		fprintf(fp_out,"                  --------          --------\n");
	}
	else
	{
		OUTPUT_TEXT("\n\n                                 Initial Parameter Values");
		OUTPUT_TEXT("\n     Variable          Model 2             Model 3             Model 4             Model 5");
		OUTPUT_TEXT("     --------          -------             -------             -------             -------");
	}

#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\nINSIDE DoInitParam(), Model %d Initial Values:\n", iModelNbr);
		if(iModelNbr == 0)
		{
			for(i=1; i <= NBR_OF_MODELS; i++)
			{
				fprintf(fp_log,"Model %d\n", i+1);
				for(j=1; j <= NBR_OF_PARMS; j++)
				{
					fprintf(fp_log,"  gppdInitPara[%d][%d] = %g\n", i, j, gppdInitPara[i][j]);
				}
			}
		}
		else
		{
			fprintf(fp_log,"Model %d\n", iModelNbr+1);
			for(j=1; j <= NBR_OF_PARMS; j++)
			{
				fprintf(fp_log,"  gppdInitPara[%d][%d] = %g\n", iModelNbr, j, gppdInitPara[iModelNbr][j]);
			}
		}
		fflush(fp_log);
	}
#endif

	for(i=1; i <= NBR_OF_PARMS; i++)
	{
		if(iModelNbr > 0)
		{
			j = sprintf(model, "%18.6g", (gppdInitPara[iModelNbr][i]));
		}
		else
		{
			for(n = 1; n <= NBR_OF_MODELS; n++)
			{
				if(giRun[n] == 0)
				{
					switch(n)
					{
					case 1:
						j = sprintf(model1, "             %s", nc);
						break;
					case 2:
						j = sprintf(model2, "             %s", nc);
						break;
					case 3:
						j = sprintf(model3, "             %s  ", nc);
						break;
					default:
						j = sprintf(model4, "             %s", nc);
						break;
					}
				}
				else
				{
					switch(n)
					{
					case 1:
						if(i == 5)
							j = sprintf(model1, "%s", "                 0");
						else if (i == 6)
							j = sprintf(model1, "%s", "                 1");
						else
							j = sprintf(model1, "%18.6g", gppdInitPara[n][i]);
						break;
					case 2:
						if(i == 5)
							j = sprintf(model2, "%s", "                 0");
						else
							j = sprintf(model2, "%18.6g", gppdInitPara[n][i]);
						break;
					case 3:
						if(i == 6)
							j = sprintf(model3, "%s", "                 1");
						else
							j = sprintf(model3, "%18.6g", gppdInitPara[n][i]);
						break;
					default:
						j = sprintf(model4, "%18.6g", gppdInitPara[n][i]);
						break;
					}
				}
			}
		}

		if(iModelNbr > 0)
		{
			fprintf(fp_out,"                   %8s%s%2s\n", gaParm_Name[i-1], model, (gppiSpecPara[iModelNbr][i] == 1 ? specified2 : specblank2));
		}
		else
		{
			fprintf(fp_out,"    %8s%s%2s%s%2s%s%2s%s%2s\n", gaParm_Name[i-1], model1, ((gppiSpecPara[1][i] == 1 && giRun[1]==1) ? specified : specblank), model2, ((gppiSpecPara[2][i] == 1 && giRun[2]==1) ? specified : specblank), 
				model3, ((gppiSpecPara[3][i] == 1 && giRun[3]==1) ? specified : specblank), model4, ((gppiSpecPara[4][i] == 1 && giRun[4]==1) ? specified : specblank) );
		}
	}

	if (iModelNbr == 0)
		fprintf(fp_out,"\n     * Indicates that this parameter has been specified\n");

	return 1;
}

int DoParamEstimates(int iModelNbr)
{
	int i, j, k, n, checkval;
	//char gaPName[6];
	char nc[]="NC";
	char model[24];
	char model1[24];
	char model2[24];
	char model3[24];
	char model4[24];
	char modelSE[24];

	char specified[] = " *";
	char specblank[] = "  ";
	char dashdash[] = "              --  ";

	checkval = 0;
	//if(gsExpoVars.iCons_Var == 0)
	//{
	//	strcpy(gaPName, "l");
	//	strcat(gaPName, gaParm_Name[0]);
	//}
	//else
	//	strcpy(gaPName, gaParm_Name[0]);

#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\nINSIDE DoParamEstimates(), Model %d Values:\n", iModelNbr);
		if(iModelNbr == 0)
		{
			for(i=1; i <= NBR_OF_MODELS; i++)
			{
				fprintf(fp_log,"Model %d\n", i+1);
				for(j=1; j <= NBR_OF_PARMS; j++)
				{
					fprintf(fp_log,"  gppdMLEs[%d][%d] = %g\n", i, j, gppdMLEs[i][j]);
				}
			}
		}
		else
		{
			fprintf(fp_log,"Model %d\n", iModelNbr+1);
			for(j=1; j <= NBR_OF_PARMS; j++)
			{
				fprintf(fp_log,"  gppdMLEs[%d][%d] = %g\n", iModelNbr, j, gppdMLEs[iModelNbr][j]);
			}
		}
		fflush(fp_log);
	}
#endif

	if(iModelNbr > 0)
	{
		OUTPUT_TEXT("\n\n\n                     Parameter Estimates\n");  /*Print computed MLEs here*/
		fprintf(fp_out,"                   %s          Model %d          Std. Err.\n", "Variable", iModelNbr+1);
		OUTPUT_TEXT("                   --------          -------          ---------");
	}
	else
	{
		OUTPUT_TEXT("\n\n\n                               Parameter Estimates by Model");  /*Print computed MLEs here*/
		OUTPUT_TEXT("\n     Variable          Model 2             Model 3            Model 4             Model 5");
		OUTPUT_TEXT("     --------          -------             -------            -------             -------");
	}
	for(i=1; i <= NBR_OF_PARMS; i++)
	{

		if(iModelNbr > 0)
		{
			if(iModelNbr == 3)
			{
				j = sprintf(model, "%20.6g", gppdMLEs[iModelNbr][i]);
				if((int)parmSE[iModelNbr][i] == -9999)
				{
					k = sprintf(modelSE, "%s", "             NA");
					checkval = 1;
				}
				else
					k = sprintf(modelSE, "%20.6g", parmSE[iModelNbr][i]);
			}
			else
			{
				j = sprintf(model, "%18.6g", gppdMLEs[iModelNbr][i]);
				if((int)parmSE[iModelNbr][i] == -9999)
				{
					k = sprintf(modelSE, "%s", "             NA");
					checkval = 1;
				}
				else
					k = sprintf(modelSE, "%18.6g", parmSE[iModelNbr][i]);
			}
		}
		else
		{
			for(n = 1; n <= NBR_OF_MODELS; n++)
			{
				if(giRun[n] == 0)
				{
					switch(n)
					{
					case 1:
						j = sprintf(model1, "              %s  ", nc);
						break;
					case 2:
						j = sprintf(model2, "              %s  ", nc);
						break;
					case 3:
						j = sprintf(model3, "              %s  ", nc);
						break;
					default:
						j = sprintf(model4, "              %s  ", nc);
						break;
					}
				}
				else
				{
					switch(n)
					{
					case 1:
						if(i == 5 || i == 6)
							j = sprintf(model1, "%s", dashdash);
						else
							j = sprintf(model1, "%18.6g", (gppdMLEs[n][i]));
						break;
					case 2:
						if(i == 5)
							j = sprintf(model2, "%s", "              --  ");
						else
							j = sprintf(model2, "%18.6g", (gppdMLEs[n][i]));
						break;
					case 3:
						if(i == 6)
							j = sprintf(model3, "%s", "              --  ");
						else
							j = sprintf(model3, "%18.6g", (gppdMLEs[n][i]));
						break;
					default:
						j = sprintf(model4, "%18.6g", (gppdMLEs[n][i]));
						break;
					}
				}
			}
		}

		if(iModelNbr > 0)
		{
			if( gppiSpecPara[iModelNbr][i] != 1)
			{
				fprintf(fp_out,"                   %8s%s%s\n", gaParm_Name[i-1], model, modelSE);
			}
		}
		else
		{
			fprintf(fp_out,"    %8s%s%2s%s%2s%s%2s%s%2s\n", gaParm_Name[i-1], model1, ((gppiSpecPara[1][i] == 1 && i==2 && giRun[1]==1) ? specified : specblank), model2, ((gppiSpecPara[2][i] == 1 && i==2 && giRun[2]==1) ? specified : specblank), model3, ((gppiSpecPara[3][i] == 1 && i==2 && giRun[3]==1) ? specified : specblank), model4, ((gppiSpecPara[4][i] == 1 && i==2 && giRun[4]==1) ? specified : specblank) );
		}
	}
	for(i = 1; i <= 4; i++)
	{
		if(giRun[i] == 0 || giBmdRun[i] == 0)
		{
			fprintf(fp_out,"\n     NC = No Convergence\n");
			break;
		}
	}

	if (iModelNbr == 0)
		fprintf(fp_out, "\n    -- Indicates that this parameter does not appear in model");

	if(gsExpoVars.iCons_Var != 0 && iModelNbr == 0)
		fprintf(fp_out,"\n     * Indicates that this parameter has been specified\n");


	// Output SE estimates for grouped output
	if(iModelNbr == 0)
	{
		OUTPUT_TEXT("\n\n\n                               Std. Err. Estimates by Model");  /*Print computed MLEs here*/
		OUTPUT_TEXT("\n     Variable          Model 2           Model 3           Model 4           Model 5");
		OUTPUT_TEXT("     --------          -------           -------           -------           -------");

	
		for(i=1; i <= NBR_OF_PARMS; i++)
		{
			if(i == 1)
			{
				for(n = 1; n <= NBR_OF_MODELS; n++)
				{
					if(giRun[n] == 0)
					{
						switch(n)
						{
						case 1:
							j = sprintf(model1, "              %s  ", nc);
							break;
						case 2:
							j = sprintf(model2, "              %s  ", nc);
							break;
						case 3:
							j = sprintf(model3, "              %s  ", nc);
							break;
						default:
							j = sprintf(model4, "              %s  ", nc);
							break;
						}
					}
					else
					{
						switch(n)
						{
						case 1:
							if ((int)parmSE[n][i] == -9999)
							{
								j = sprintf(model1, "%s", "              NA  ");
								checkval = 1;
							}
							else
								j = sprintf(model1, "%18.6g", parmSE[n][i]);
							break;
						case 2:
							if ((int)parmSE[n][i] == -9999)
							{
								j = sprintf(model2, "%s", "              NA  ");
								checkval = 1;
							}
							else
								j = sprintf(model2, "%18.6g", parmSE[n][i]);
							break;
						case 3:
							if ((int)parmSE[n][i] == -9999)
							{
								j = sprintf(model3, "%s", "              NA  ");
								checkval = 1;
							}
							else
								j = sprintf(model3, "%18.6g", parmSE[n][i]);
							break;
						case 4:
							if ((int)parmSE[n][i] == -9999)
							{
								j = sprintf(model4, "%s", "              NA  ");
								checkval = 1;
							}
							else
								j = sprintf(model4, "%18.6g", parmSE[n][i]);
							break;
						}
					}
				}
			}
			else
			{
				for(n = 1; n <= NBR_OF_MODELS; n++)
				{
					if(giRun[n] == 0)
					{
						switch(n)
						{
						case 1:
							j = sprintf(model1, "              %s  ", nc);
							break;
						case 2:
							j = sprintf(model2, "              %s  ", nc);
							break;
						case 3:
							j = sprintf(model3, "              %s  ", nc);
							break;
						default:
							j = sprintf(model4, "              %s  ", nc);
							break;
						}
					}
					else
					{
						switch(n)
						{
						case 1:
							if(i == 5 || i == 6 || (int)parmSE[n][i] == -9999 )
							{
								j = sprintf(model1, "%s", "              NA  ");
								checkval = 1;
							}
							else
								j = sprintf(model1, "%18.6g", parmSE[n][i]);
							break;
						case 2:
							if(i == 5  || (int)parmSE[n][i] == -9999)
							{
								j = sprintf(model2, "%s", "              NA  ");
								checkval = 1;
							}
							else
								j = sprintf(model2, "%18.6g", parmSE[n][i]);
								break;
						case 3:
							if(i == 6 || (int)parmSE[n][i] == -9999 )
							{
								j = sprintf(model3, "%s", "              NA  ");
								checkval = 1;
							}
							else
								j = sprintf(model3, "%18.6g", parmSE[n][i]);
								break;
						default:
							if ((int)parmSE[n][i] == -9999)
							{
								j = sprintf(model4, "%s", "              NA  ");
								checkval = 1;
							}
							else
								j = sprintf(model4, "%18.6g", parmSE[n][i]);
							break;
						}
					}
				}
			}

			fprintf(fp_out,"    %8s%s%s%s%s\n", gaParm_Name[i-1], model1, model2, model3, model4 );
		}

	}

	if ((checkval = 1 && iModelNbr == 0))
	{
		fprintf(fp_out, "\nNA - Indicates that this parameter was specified (by the user or because of the model form)");
		fprintf(fp_out, "\n     or has hit a bound implied by some inequality constraint and thus has no standard error.");
	}	



	return i;
}

void DoDataEstimateInt(int iModelNbr)
{
	int i, n;
	OUTPUT_TEXT("\n\n            Table of Stats From Input Data");
	if (gsExpoVars.iIn_Type==1)
	{
		if(gsExpoVars.iLogNormal==0)
			OUTPUT_TEXT("\n     Dose      N         Obs Mean     Obs Std Dev");
		else
			OUTPUT_TEXT("\n     Dose      N     Calc'd Median   Calc'd GSD");
		OUTPUT_TEXT("     -----    ---       ----------   -------------");	
		for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
		{
			if(gsExpoVars.iLogNormal==0)
				fprintf(fp_out,"     %5.4g%7d%13.4g%13.4g\n", gpdXi[i], gpiNi[i], gpdYm[i], gpdStDev[i]);
			else
				fprintf(fp_out,"     %5.4g%7d%13.4g%13.4g\n", gpdXi[i], gpiNi[i], exp(gpdYm[i]), exp(gpdStDev[i]));
		}
	}
	else
	{
		OUTPUT_TEXT("\n     Dose      Obs Mean");
		OUTPUT_TEXT("     ------       ----------");	
		for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
		{
			if(gsExpoVars.iLogNormal==0)
				fprintf(fp_out,"     %5.4g%10.4g\n", gpdXi[i], gpdYm[i]);
			else
				fprintf(fp_out,"     %5.4g%10.4g\n", gpdXi[i], exp(gpdYm[i]));
		}
	}
	char caModel[1];
	caModel[0] = ' ';
	double dalpha=0.0;
	double drho=0.0;
	double da=0.0;
	double db=0.0;
	double dc=0.0;
	double dd=0.0;

	double dEstMean=0.0;
	double dEstStd=0.0;
	double dChi2=0.0;
	double dCalcMean = 0.0;
	if(iModelNbr > 0)
		OUTPUT_TEXT("\n\n                  Estimated Values of Interest");
	else
		OUTPUT_TEXT("\n\n                      Estimated Values of Interest");

	if(gsExpoVars.iLogNormal==0)
	{
		if(iModelNbr > 0)
			OUTPUT_TEXT("\n      Dose      Est Mean      Est Std     Scaled Residual");  /* Compute Values here*/
		else
			OUTPUT_TEXT("\n      Model      Dose      Est Mean      Est Std     Scaled Residual");  /* Compute Values here*/
	}
	else
	{
		if(iModelNbr > 0)
			OUTPUT_TEXT("\n      Dose     Est Median     Est GSD      Scaled Residual");  /* Compute Values here*/
		else
			OUTPUT_TEXT("\n      Model      Dose     Est Median     Est GSD      Scaled Residual");  /* Compute Values here*/
	}

	if(iModelNbr > 0)
		OUTPUT_TEXT("    ------    ----------    ---------    ----------------");
	else
		OUTPUT_TEXT("     -------    ------    ----------    ---------    ----------------");
	for(i=1; i <= NBR_OF_MODELS; i++)
	{
		if(giRun[i] == 0)
			continue;

		if(iModelNbr > 0 && i != iModelNbr)
			continue;

		for(n = 1; n <= gsExpoVars.iNbrObs_Total; n++)
		{
			//dalpha = log(gppdMLEs[i][(int) eAlpha]);
			dalpha = gppdMLEs[i][(int) eAlpha];
			drho = gppdMLEs[i][(int) eRho];
			da = gppdMLEs[i][(int) ea];
			db = gppdMLEs[i][(int) eb];
			dc = gppdMLEs[i][(int) ec];
			dd = gppdMLEs[i][(int) ed];
			if(i == 1 && n == 1)
			{
				caModel[0]='2';
			}
			if(i == 2 && n == 1)
			{
				caModel[0]='3';
			}
			if(i == 3 && n == 1)
			{
				caModel[0]='4';
			}
			if(i == 4 && n == 1)
			{
				caModel[0]='5';
			}

			//dEstMean = da * (dc - (dc - 1) * exp(db * pow(gpdXi[n], dd)) );
			dEstMean = gppdPredict[i][n];
			dCalcMean = gpdYm[n];
			if(gsExpoVars.iLogNormal==1)
				dEstMean = exp(dEstMean);
			//dEstStd =  sqrt(exp(dalpha + log(gpdYm[n]) * drho));
			dEstStd = sqrt(exp(gppdMLEs[i][1] + gppdMLEs[i][2] * Slog(dEstMean)));
			if(gsExpoVars.iLogNormal==1)
			{
				dEstStd = exp(dEstStd);
				dCalcMean = exp(gpdYm[n]);
			}

			//dChi2 = (gpdYm[n] - dEstMean)/(dEstStd/sqrt(gpiNi[n]));
			dChi2 = (dCalcMean - dEstMean)/(dEstStd/sqrt(gpiNi[n]));

//#ifdef DO_LOG
//			if (giDo_Log)	// Print values to log for investigation
//			{
//				fprintf(fp_log,"\n\nEstimated Values of Interest (gsExpoVars.iNbrObs_Total=%d): i=%d, n=%d;\ndalpha=%g\ndrho=%g\nda=%g\ndb=%g\ndc=%g\ndd=%g",
//					gsExpoVars.iNbrObs_Total, i,n,dalpha,drho,da,db,dc,dd);
//				fprintf(fp_log,"\n\ngppdMLEs[%d][eAlpha]=%g, dEstMean=%g, dEstStd=%g, dChi2=%g", i, gppdMLEs[i][(int) eAlpha], dEstMean, dEstStd, dChi2);
//			}
//#endif

			if(iModelNbr > 0)
				fprintf(fp_out,"%10.4g%14.4g%13.4g%17.4g\n", 
				gpdXi[n],
				dEstMean,
				dEstStd,
				dChi2);
			else
				fprintf(fp_out,"     %6c%10.4g%14.4g%13.4g%17.4g\n", 
				n==1?caModel[0]:' ',
				gpdXi[n],
				dEstMean,
				dEstStd,
				dChi2);
		}
	}
}

int DoLikehoods(double A1, double A2, double A3, double R, int iModelNbr)
{
	fprintf(fp_out,"\n\n\n   Other models for which likelihoods are calculated:");
	fprintf(fp_out,"\n\n     Model A1:        Yij = Mu(i) + e(ij)");
	fprintf(fp_out,"\n               Var{e(ij)} = Sigma^2");

	fprintf(fp_out,"\n\n     Model A2:        Yij = Mu(i) + e(ij)");
	fprintf(fp_out,"\n               Var{e(ij)} = Sigma(i)^2");

	fprintf(fp_out,"\n\n     Model A3:        Yij = Mu(i) + e(ij)");
	fprintf(fp_out,"\n               Var{e(ij)} = exp(lalpha + log(mean(i)) * rho)");

	fprintf(fp_out,"\n\n     Model  R:        Yij = Mu + e(i)");
	fprintf(fp_out,"\n               Var{e(ij)} = Sigma^2");

	typedef enum {
		eA1=1, eA2, eA3, eMR, eM2, eM3, eM4, eM5
	} eModels;

	int iLikeIntList = 8;		//Nbr of elements in a Likelihood Interest List
	int iTstOfIntLst = 11;		//Nbr of elements in a Tests Interest List
	int i;
	LikeInterestList *aLikeIntList;
	LikeInterestList *aTstOfIntLst;
	aLikeIntList = (LikeInterestList *) malloc((size_t) ((iLikeIntList-1+1+NR_END)*sizeof(LikeInterestList))); 
	if (!aLikeIntList) 
	{
		ERRORPRT ("Memory allocation failed in aLikeIntList");
#ifdef DO_LOG
		if (giDo_Log) 
		{
			if(fp_log != NULL)
			{
				fflush(fp_log);
				fclose(fp_log);
			}
		}
#endif

		CLOSE_FILES ();
		exit(1);
	}

	aTstOfIntLst = (LikeInterestList *) malloc((size_t) ((iTstOfIntLst-1+1+NR_END)*sizeof(LikeInterestList))); 
	if (!aTstOfIntLst) 
	{
		ERRORPRT ("Memory allocation failed in aTstOfIntLst");
#ifdef DO_LOG
		if (giDo_Log) 
		{
			if(fp_log != NULL)
			{
				fflush(fp_log);
				fclose(fp_log);
			}
		}
#endif

		CLOSE_FILES ();
		exit(1);
	}

	strcpy(aLikeIntList[(int) eA1].caModel, "A1");
	strcpy(aLikeIntList[(int) eA2].caModel, "A2");
	strcpy(aLikeIntList[(int) eA3].caModel, "A3");
	strcpy(aLikeIntList[(int) eMR].caModel, "R");
	strcpy(aLikeIntList[(int) eM2].caModel, "2");
	strcpy(aLikeIntList[(int) eM3].caModel, "3");
	strcpy(aLikeIntList[(int) eM4].caModel, "4");
	strcpy(aLikeIntList[(int) eM5].caModel, "5");

	aLikeIntList[(int) eA1].dLogLikelihood = A1;
	aLikeIntList[(int) eA2].dLogLikelihood = A2;
	aLikeIntList[(int) eA3].dLogLikelihood = A3;
	aLikeIntList[(int) eMR].dLogLikelihood = R;
	aLikeIntList[(int) eM2].dLogLikelihood = gpdLikelihoods[1];
	aLikeIntList[(int) eM3].dLogLikelihood = gpdLikelihoods[2];
	aLikeIntList[(int) eM4].dLogLikelihood = gpdLikelihoods[3];
	aLikeIntList[(int) eM5].dLogLikelihood = gpdLikelihoods[4];

	//Computation for DF and AIC
	aLikeIntList[1].iDF = gsExpoVars.iNbrObs_Total + 1;
	aLikeIntList[2].iDF = (gsExpoVars.iNbrObs_Total*2);
	aLikeIntList[3].iDF = gsExpoVars.iNbrObs_Total + (gsExpoVars.iCons_Var==0?2:1);
	aLikeIntList[4].iDF = 2;

	aLikeIntList[5].iDF = 2 + (gsExpoVars.iCons_Var==0?2:1) - 
		(gppdMLEs[1][(int) ea]==0?1:0);

	aLikeIntList[6].iDF = 3 + (gsExpoVars.iCons_Var==0?2:1) - 
		(gppdMLEs[2][(int) ea]==0?1:0 + gppdMLEs[2][(int) ed]==1?1:0);

	aLikeIntList[7].iDF = 3 + (gsExpoVars.iCons_Var==0?2:1) - 
		(gppdMLEs[3][(int) ea]==0?1:0 + gppdMLEs[3][(int) eb]==0?1:0 + gppdMLEs[3][(int) ec]==0?1:0);

	aLikeIntList[8].iDF = 4 + (gsExpoVars.iCons_Var==0?2:1) -
		(gppdMLEs[4][(int) ea]==0?1:0 + gppdMLEs[4][(int) eb]==0?1:0 + gppdMLEs[4][(int) ec]==0?1:0 + gppdMLEs[4][(int) ed]==1?1:0);

	for(i = 1; i <= iLikeIntList; i++)
	{
		aLikeIntList[i].dAIC_PValue = (-2*aLikeIntList[i].dLogLikelihood) + (2 * aLikeIntList[i].iDF);
	}

	OUTPUT_TEXT("\n\n\n                                Likelihoods of Interest");
	OUTPUT_TEXT("\n                     Model      Log(likelihood)      DF         AIC");
	fprintf(fp_out,"                    -------    -----------------    ----   ------------");	
	for(i = 1; i <= iLikeIntList; i++)
	{
		if(i > 4 && giRun[i-4]==0)
			continue;
		if(i > 4 && iModelNbr > 0 && i-4 != iModelNbr)
			continue;
		fprintf(fp_out,"\n                  %8s  %14.7g  %11d  %12.7g",
			aLikeIntList[i].caModel, 
			aLikeIntList[i].dLogLikelihood, 
			aLikeIntList[i].iDF,
			aLikeIntList[i].dAIC_PValue);	//Changed significant digits to 7 as per Bruce's "Output file changes for exponential" 08/09/07 email
	}

	double dX, dSumN;
	dX = dSumN = 0.0;
	for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
		dSumN += gpiNi[i];
	dX = -(dSumN * log(TwoPi)/2);
	fprintf(fp_out,"\n\n\n   Additive constant for all log-likelihoods = %10.4g.  This constant added to the", dX);
	fprintf(fp_out,"\n   above values gives the log-likelihood including the term that does not");
	fprintf(fp_out,"\n   depend on the model parameters.");

	OUTPUT_TEXT("\n\n\n                                 Explanation of Tests");
	fprintf(fp_out,"\n   Test 1:  Does response and/or variances differ among Dose levels? (A2 vs. R)");

	fprintf(fp_out,"\n   Test 2:  Are Variances Homogeneous? (A2 vs. A1)");

	fprintf(fp_out,"\n   Test 3:  Are variances adequately modeled? (A2 vs. A3)");

	if(giRun[1]==1 && (iModelNbr == 0 || iModelNbr == 1))
		fprintf(fp_out,"\n   Test 4:  Does Model 2 fit the data? (A3 vs. 2)");

	if(giRun[2]==1 && (iModelNbr == 0 || iModelNbr == 2))
		fprintf(fp_out,"\n\n   Test 5a: Does Model 3 fit the data? (A3 vs 3)");

	if(iModelNbr == 0)
	{
		if(giRun[1]==1 && giRun[2]==1)
			fprintf(fp_out,"\n   Test 5b: Is Model 3 better than Model 2? (3 vs. 2)");
	}

	if(giRun[3]==1 && (iModelNbr == 0 || iModelNbr == 3))
		fprintf(fp_out,"\n\n   Test 6a: Does Model 4 fit the data? (A3 vs 4)");

	if(iModelNbr == 0)
	{
		if(giRun[1]==1 && giRun[3]==1)
			fprintf(fp_out,"\n   Test 6b: Is Model 4 better than Model 2? (4 vs. 2)");
	}

	if(giRun[4]==1 && (iModelNbr == 0 || iModelNbr == 4))
		fprintf(fp_out,"\n\n   Test 7a: Does Model 5 fit the data? (A3 vs 5)");

	if(iModelNbr == 0)
	{
		if(giRun[4]==1 && giRun[2]==1)
			fprintf(fp_out,"\n   Test 7b: Is Model 5 better than Model 3? (5 vs. 3)");
	}

	if(iModelNbr == 0)
	{
		if(giRun[4]==1 && giRun[3]==1)
			fprintf(fp_out,"\n   Test 7c: Is Model 5 better than Model 4? (5 vs. 4)");
	}

	//Initialized the Test of Interest List
	strcpy(aTstOfIntLst[1].caModel, "Test 1");
	strcpy(aTstOfIntLst[2].caModel, "Test 2");
	strcpy(aTstOfIntLst[3].caModel, "Test 3");
	strcpy(aTstOfIntLst[4].caModel, "Test 4");
	strcpy(aTstOfIntLst[5].caModel, "Test 5a");
	strcpy(aTstOfIntLst[6].caModel, "Test 5b");
	strcpy(aTstOfIntLst[7].caModel, "Test 6a");
	strcpy(aTstOfIntLst[8].caModel, "Test 6b");
	strcpy(aTstOfIntLst[9].caModel, "Test 7a");
	strcpy(aTstOfIntLst[10].caModel,"Test 7b");
	strcpy(aTstOfIntLst[11].caModel,"Test 7c");

	for(i = 1; i <= iTstOfIntLst; i++)
	{
		aTstOfIntLst[i].dLogLikelihood = 0.0;
		aTstOfIntLst[i].iDF = 0;
		aTstOfIntLst[i].dAIC_PValue = 0.0;
	}

	// Computation here
	for(i = 1; i <= iTstOfIntLst; i++)
	{
		//Test 1:  Does response and/or variances differ among Dose levels? (A2 vs. R)
		if( strcmp(aTstOfIntLst[i].caModel, "Test 1") == 0)
		{
			aTstOfIntLst[i].dLogLikelihood = 2 * (aLikeIntList[(int) eA2].dLogLikelihood - aLikeIntList[(int) eMR].dLogLikelihood);	//Likelihood ratio
			aTstOfIntLst[i].iDF = aLikeIntList[(int) eA2].iDF - aLikeIntList[(int) eMR].iDF;						// Difference DF
		}

		//Test 2:  Are Variances Homogeneous? (A1 vs. A2)
		if( strcmp(aTstOfIntLst[i].caModel, "Test 2") == 0)
		{
			aTstOfIntLst[i].dLogLikelihood = 2 * (aLikeIntList[(int) eA2].dLogLikelihood - aLikeIntList[(int) eA1].dLogLikelihood);	//Likelihood ratio
			aTstOfIntLst[i].iDF = aLikeIntList[(int) eA2].iDF - aLikeIntList[(int) eA1].iDF;						// Difference DF
		}

		//Test 3:  Are variances adequately modeled? (A2 vs. A3)
		if( strcmp(aTstOfIntLst[i].caModel, "Test 3") == 0)
		{
			aTstOfIntLst[i].dLogLikelihood = 2 * (aLikeIntList[(int) eA2].dLogLikelihood - aLikeIntList[(int) eA3].dLogLikelihood);	//Likelihood ratio
			aTstOfIntLst[i].iDF = aLikeIntList[(int) eA2].iDF - aLikeIntList[(int) eA3].iDF;						// Difference DF
		}

		//Test 4:  Does Model 2 fit the data? (2 vs A3)
		if( strcmp(aTstOfIntLst[i].caModel, "Test 4") == 0)
		{
			aTstOfIntLst[i].dLogLikelihood = 2 * (aLikeIntList[(int) eA3].dLogLikelihood - aLikeIntList[(int) eM2].dLogLikelihood);	//Likelihood ratio
			aTstOfIntLst[i].iDF = aLikeIntList[(int) eA3].iDF - aLikeIntList[(int) eM2].iDF;						// Difference DF
		}

		//Test 5a: Does Model 3 fit the data? (3 vs A3)
		if( strcmp(aTstOfIntLst[i].caModel, "Test 5a") == 0)
		{
			aTstOfIntLst[i].dLogLikelihood = 2 * (aLikeIntList[(int) eA3].dLogLikelihood - aLikeIntList[(int) eM3].dLogLikelihood);	//Likelihood ratio
			aTstOfIntLst[i].iDF = aLikeIntList[(int) eA3].iDF - aLikeIntList[(int) eM3].iDF;						// Difference DF
		}

		//Test 5b: Is Model 3 better than Model 2? (3 vs. 2)
		if( strcmp(aTstOfIntLst[i].caModel, "Test 5b") == 0)
		{
			aTstOfIntLst[i].dLogLikelihood = 2 * (aLikeIntList[(int) eM3].dLogLikelihood - aLikeIntList[(int) eM2].dLogLikelihood);	//Likelihood ratio
			aTstOfIntLst[i].iDF = aLikeIntList[(int) eM3].iDF - aLikeIntList[(int) eM2].iDF;						// Difference DF
		}

		//Test 6a: Does Model 4 fit the data? (4 vs A3)
		if( strcmp(aTstOfIntLst[i].caModel, "Test 6a") == 0)
		{
			aTstOfIntLst[i].dLogLikelihood = 2 * (aLikeIntList[(int) eA3].dLogLikelihood - aLikeIntList[(int) eM4].dLogLikelihood);	//Likelihood ratio
			aTstOfIntLst[i].iDF = aLikeIntList[(int) eA3].iDF - aLikeIntList[(int) eM4].iDF;						// Difference DF
		}

		//Test 6b: Is Model 4 better than Model 2? (4 vs. 2)
		if( strcmp(aTstOfIntLst[i].caModel, "Test 6b") == 0)
		{
			aTstOfIntLst[i].dLogLikelihood = 2 * (aLikeIntList[(int) eM4].dLogLikelihood - aLikeIntList[(int) eM2].dLogLikelihood);	//Likelihood ratio
			aTstOfIntLst[i].iDF = aLikeIntList[(int) eM4].iDF - aLikeIntList[(int) eM2].iDF;						// Difference DF
		}

		//Test 7a: Does Model 5 fit the data? (5 vs A3)
		if( strcmp(aTstOfIntLst[i].caModel, "Test 7a") == 0)
		{
			aTstOfIntLst[i].dLogLikelihood = 2 * (aLikeIntList[(int) eA3].dLogLikelihood - aLikeIntList[(int) eM5].dLogLikelihood);	//Likelihood ratio
			aTstOfIntLst[i].iDF = aLikeIntList[(int) eA3].iDF - aLikeIntList[(int) eM5].iDF;						// Difference DF
		}

		//Test 7b: Is Model 5 better than Model 3? (5 vs. 3)
		if( strcmp(aTstOfIntLst[i].caModel,"Test 7b") == 0)
		{
			aTstOfIntLst[i].dLogLikelihood = 2 * (aLikeIntList[(int) eM5].dLogLikelihood - aLikeIntList[(int) eM3].dLogLikelihood);	//Likelihood ratio
			aTstOfIntLst[i].iDF = aLikeIntList[(int) eM5].iDF - aLikeIntList[(int) eM3].iDF;						// Difference DF
		}

		//Test 7c: Is Model 5 better than Model 4? (5 vs. 4)
		if( strcmp(aTstOfIntLst[i].caModel,"Test 7c") == 0)
		{
			aTstOfIntLst[i].dLogLikelihood = 2 * (aLikeIntList[(int) eM5].dLogLikelihood - aLikeIntList[(int) eM4].dLogLikelihood);	//Likelihood ratio
			aTstOfIntLst[i].iDF = aLikeIntList[(int) eM5].iDF - aLikeIntList[(int) eM4].iDF;						// Difference DF
		}
		//aTstOfIntLst[i].dLogLikelihood = fabs(aTstOfIntLst[i].dLogLikelihood);
		aTstOfIntLst[i].dLogLikelihood = aTstOfIntLst[i].dLogLikelihood;
		//aTstOfIntLst[i].iDF = fabs(aTstOfIntLst[i].iDF);
		if( (aTstOfIntLst[i].dLogLikelihood < 0.0) || (aTstOfIntLst[i].iDF <= 0) )
			aTstOfIntLst[i].dAIC_PValue = -1;
		else
			aTstOfIntLst[i].dAIC_PValue = CHISQ(aTstOfIntLst[i].dLogLikelihood, aTstOfIntLst[i].iDF);	// chi-square
	}

	OUTPUT_TEXT("\n\n\n                            Tests of Interest");
	OUTPUT_TEXT("\n     Test          -2*log(Likelihood Ratio)       D. F.         p-value");
	fprintf(fp_out, "   --------        ------------------------      ------     --------------");
	for(i = 1; i <= iTstOfIntLst; i++)
	{
		if(i > 3)
		{
			if((i==4 && giRun[1]==0) || (i==5 && giRun[2]==0) || (i==6 && (giRun[1]==0 || giRun[2]==0)) ||(i==7 && giRun[3]==0) || 
				(i==8 && (giRun[1]==0 || giRun[3]==0)) || (i==9 && giRun[4]==0) || (i==10 && (giRun[4]==0 || giRun[2]==0)) || (i==11 && (giRun[4]==0 || giRun[3]==0)) )
				continue;

			if(iModelNbr > 0)
			{
				if( strcmp(aTstOfIntLst[i].caModel, "Test 4") == 0 && iModelNbr != 1)
					continue;
				if( strcmp(aTstOfIntLst[i].caModel, "Test 5a") == 0 && iModelNbr != 2)
					continue;
				if( strcmp(aTstOfIntLst[i].caModel, "Test 5b") == 0)
					continue;
				if( strcmp(aTstOfIntLst[i].caModel, "Test 6a") == 0 && iModelNbr != 3)
					continue;
				if( strcmp(aTstOfIntLst[i].caModel, "Test 6b") == 0)
					continue;
				if( strcmp(aTstOfIntLst[i].caModel, "Test 7a") == 0 && iModelNbr != 4)
					continue;
				if( strcmp(aTstOfIntLst[i].caModel, "Test 7b") == 0)
					continue;
				if( strcmp(aTstOfIntLst[i].caModel, "Test 7c") == 0)
					continue;
			}
		}
		if(aTstOfIntLst[i].dAIC_PValue == -1)
			fprintf(fp_out, "\n   %8s        %22.4g      %6d                 N/A",
			aTstOfIntLst[i].caModel,
			aTstOfIntLst[i].dLogLikelihood,
			aTstOfIntLst[i].iDF);
		else if(aTstOfIntLst[i].dAIC_PValue < 0.0001)
			fprintf(fp_out, "\n   %8s        %22.4g      %6d            < 0.0001",
			aTstOfIntLst[i].caModel,
			aTstOfIntLst[i].dLogLikelihood,
			aTstOfIntLst[i].iDF);
		else
			fprintf(fp_out, "\n   %8s        %22.4g      %6d        %12.4g",
			aTstOfIntLst[i].caModel,
			aTstOfIntLst[i].dLogLikelihood,
			aTstOfIntLst[i].iDF,
			aTstOfIntLst[i].dAIC_PValue );

	}

	//The following are the output conclusions.
	//Test 1:  Does response and/or variances differ among Dose levels? (A2 vs. R)"
	if(aTstOfIntLst[1].dAIC_PValue < 0.05)
	{
		fprintf(fp_out, "\n\n\n     The p-value for Test 1 is less than .05.  There appears to be a");
		fprintf(fp_out, "\n     difference between response and/or variances among the dose");
		fprintf(fp_out, "\n     levels, it seems appropriate to model the data.");
	}
	else
	{
		fprintf(fp_out, "\n\n\n     The p-value for Test 1 is greater than .05.  There may not be a");
		fprintf(fp_out, "\n     diffence between responses and/or variances among the dose levels");
		fprintf(fp_out, "\n     Modelling the data with a dose/response curve may not be appropriate.");
	}

	//Test 2:  Are Variances Homogeneous? (A1 vs. A2)
	if(gsExpoVars.iCons_Var == 1)
	{
		if(aTstOfIntLst[2].iDF <= 0)
		{
			fprintf(fp_out, "\n\n     Degrees of freedom for Test 2 are less than or equal to 0.");
			fprintf(fp_out, "\n     The Chi-Square test for fit is not valid.");
		}
		else if(aTstOfIntLst[2].dAIC_PValue < 0.1)
		{
			fprintf(fp_out, "\n\n     The p-value for Test 2 is less than .1.  Consider running");
			fprintf(fp_out, "\n     a non-homogeneous variance model.");
		}
		else
		{
			fprintf(fp_out, "\n\n     The p-value for Test 2 is greater than .1.  A homogeneous");
			fprintf(fp_out, "\n     variance model appears to be appropriate here.");
		}
	}
	else
	{
		if(aTstOfIntLst[2].iDF <= 0)
		{
			fprintf(fp_out, "\n\n     Degrees of freedom for Test 2 are less than or equal to 0.");
			fprintf(fp_out, "\n     The Chi-Square test for fit is not valid.");
		}
		else if(aTstOfIntLst[2].dAIC_PValue < 0.1)
		{
			fprintf(fp_out, "\n\n     The p-value for Test 2 is less than .1.  A non-homogeneous");
			fprintf(fp_out, "\n     variance model appears to be appropriate.");
		}
		else
		{
			fprintf(fp_out, "\n\n     The p-value for Test 2 is greater than .1.  Consider");
			fprintf(fp_out, "\n     running a homogeneous model.");
		}
	}

	//Test 3:  Are variances adequately modeled? (A2 vs. A3)
	if(aTstOfIntLst[3].iDF <= 0)
	{
		fprintf(fp_out, "\n\n     Degrees of freedom for Test 3 are less than or equal to 0.");
		fprintf(fp_out, "\n     The Chi-Square test for fit is not valid.");
	}
	else if(aTstOfIntLst[3].dAIC_PValue < 0.1)
	{
		fprintf(fp_out, "\n\n     The p-value for Test 3 is less than .1.  You may want to");
		fprintf(fp_out, "\n     consider a different variance model.");
	}
	else
	{
		fprintf(fp_out, "\n\n     The p-value for Test 3 is greater than .1.  The modeled");
		fprintf(fp_out, "\n     variance appears to be appropriate here.");
	}

	//Test 4:  Does Model 2 fit the data? (2 vs A3)
	if(giRun[1]==1)
	{
		if(iModelNbr > 0 && iModelNbr != 1)
			goto Test5;

		if(aTstOfIntLst[4].iDF <= 0)
		{
			fprintf(fp_out, "\n\n     Degrees of freedom for Test 4 are less than or equal to 0.");
			fprintf(fp_out, "\n     The Chi-Square test for fit is not valid.");
		}
		else if(aTstOfIntLst[4].dAIC_PValue < 0.1)
		{
			fprintf(fp_out, "\n\n     The p-value for Test 4 is less than .1.  Model 2 may not adequately");
			fprintf(fp_out, "\n     describe the data; you may want to consider another model.");
		}
		else
		{
			fprintf(fp_out, "\n\n     The p-value for Test 4 is greater than .1.  Model 2 seems");
			fprintf(fp_out, "\n     to adequately describe the data.");
		}
	}

	//Test 5a: Does Model 3 fit the data? (3 vs A3)
Test5:
	if(giRun[2]==1)
	{
		if(iModelNbr > 0 && iModelNbr != 2)
			goto Test5b;

		if(aTstOfIntLst[5].iDF <= 0)
		{
			fprintf(fp_out, "\n\n     Degrees of freedom for Test 5a are less than or equal to 0.");
			fprintf(fp_out, "\n     The Chi-Square test for fit is not valid.");
		}
		else if(aTstOfIntLst[5].dAIC_PValue < 0.1)
		{
			fprintf(fp_out, "\n\n     The p-value for Test 5a is less than .1.  Model 3 may not adequately");
			fprintf(fp_out, "\n     describe the data; you may want to consider another model.");
		}
		else
		{
			fprintf(fp_out, "\n\n     The p-value for Test 5a is greater than .1.  Model 3 seems");
			fprintf(fp_out, "\n     to adequately describe the data.");
		}
	}

	//Test 5b: Is Model 3 better than Model 2? (3 vs. 2)
Test5b:
	if(iModelNbr > 0)
		goto Test6a;
	if(giRun[1]==1 && giRun[2]==1)
	{
		if(aTstOfIntLst[6].iDF <= 0)
		{
			fprintf(fp_out, "\n\n     Degrees of freedom for Test 5b are less than or equal to 0.");
			fprintf(fp_out, "\n     The Chi-Square test for fit is not valid.");
		}
		else
		{
			if(aTstOfIntLst[6].dAIC_PValue < 0.05)
			{
				fprintf(fp_out, "\n\n     The p-value for Test 5b is less than .05.  Model 3 appears");
				fprintf(fp_out, "\n     to fit the data better than Model 2.");
			}
			else
			{
				fprintf(fp_out, "\n\n     The p-value for Test 5b is greater than .05.  Model 3 does");
				fprintf(fp_out, "\n     not seem to fit the data better than Model 2.");
			}
		}
	}

	//Test 6a: Does Model 4 fit the data? (4 vs A3)
Test6a:
	if(giRun[3]==1)
	{
		if(iModelNbr > 0 && iModelNbr != 3)
			goto Test6b;

		if(aTstOfIntLst[7].iDF <= 0)
		{
			fprintf(fp_out, "\n\n     Degrees of freedom for Test 6a are less than or equal to 0.");
			fprintf(fp_out, "\n     The Chi-Square test for fit is not valid.");
		}
		else if(aTstOfIntLst[7].dAIC_PValue < 0.1)
		{
			fprintf(fp_out, "\n\n     The p-value for Test 6a is less than .1.  Model 4 may not adequately");
			fprintf(fp_out, "\n     describe the data; you may want to consider another model.");
		}
		else
		{
			fprintf(fp_out, "\n\n     The p-value for Test 6a is greater than .1.  Model 4 seems");
			fprintf(fp_out, "\n     to adequately describe the data.");
		}
	}

	//Test 6b: Is Model 4 better than Model 2? (4 vs. 2)
Test6b:
	if(iModelNbr > 0)
		goto Test7a;
	if(giRun[1]==1 && giRun[3]==1)
	{
		if(aTstOfIntLst[8].iDF <= 0)
		{
			fprintf(fp_out, "\n\n     Degrees of freedom for Test 6b are less than or equal to 0.");
			fprintf(fp_out, "\n     The Chi-Square test for fit is not valid.");
		}
		else
		{
			if(aTstOfIntLst[8].dAIC_PValue < 0.05)
			{
				fprintf(fp_out, "\n\n     The p-value for Test 6b is less than .05.  Model 4 appears");
				fprintf(fp_out, "\n     to fit the data better than Model 2.");
			}
			else
			{
				fprintf(fp_out, "\n\n     The p-value for Test 6b is greater than .05.  Model 4 does");
				fprintf(fp_out, "\n     not seem to fit the data better than Model 2.");
			}
		}
	}

	//Test 7a: Does Model 5 fit the data? (5 vs A3)
Test7a:
	if(giRun[4]==1)
	{
		if(iModelNbr > 0 && iModelNbr != 4)
			goto Test7b;

		if(aTstOfIntLst[9].iDF <= 0)
		{
			fprintf(fp_out, "\n\n     Degrees of freedom for Test 7a are less than or equal to 0.");
			fprintf(fp_out, "\n     The Chi-Square test for fit is not valid.");
		}
		else if(aTstOfIntLst[9].dAIC_PValue < 0.1)
		{
			fprintf(fp_out, "\n\n     The p-value for Test 7a is less than .1.  Model 5 may not adequately");
			fprintf(fp_out, "\n     describe the data; you may want to consider another model.");
		}
		else
		{
			fprintf(fp_out, "\n\n     The p-value for Test 7a is greater than .1.  Model 5 seems");
			fprintf(fp_out, "\n     to adequately describe the data.");
		}
	}

	//Test 7b: Is Model 5 better than Model 3? (5 vs. 3)
Test7b:
	if(iModelNbr > 0)
		goto Test7c;
	if(giRun[4]==1 && giRun[2]==1)
	{
		if(aTstOfIntLst[10].iDF <= 0)
		{
			fprintf(fp_out, "\n\n     Degrees of freedom for Test 7b are less than or equal to 0.");
			fprintf(fp_out, "\n     The Chi-Square test for fit is not valid.");
		}
		else
		{
			if(aTstOfIntLst[10].dAIC_PValue < 0.05)
			{
				fprintf(fp_out, "\n\n     The p-value for Test 7b is less than .05.  Model 5 appears");
				fprintf(fp_out, "\n     to fit the data better than Model 3.");
			}
			else
			{
				fprintf(fp_out, "\n\n     The p-value for Test 7b is greater than .05.  Model 5 does");
				fprintf(fp_out, "\n     not seem to fit the data better than Model 3.");
			}
		}
	}

	//Test 7c: Is Model 5 better than Model 4? (5 vs. 4)
Test7c:
	if(iModelNbr > 0)
		goto TestEnd;
	if(giRun[4]==1 && giRun[3]==1)
	{
		if(aTstOfIntLst[11].iDF <= 0)
		{
			fprintf(fp_out, "\n\n     Degrees of freedom for Test 7c are less than or equal to 0.");
			fprintf(fp_out, "\n     The Chi-Square test for fit is not valid.");
		}
		else
		{
			if(aTstOfIntLst[11].dAIC_PValue < 0.05)
			{
				fprintf(fp_out, "\n\n     The p-value for Test 7c is less than .05.  Model 5 appears");
				fprintf(fp_out, "\n     to fit the data better than Model 4.");
			}
			else
			{
				fprintf(fp_out, "\n\n     The p-value for Test 7c is greater than .05.  Model 5 does");
				fprintf(fp_out, "\n     not seem to fit the data better than Model 4.");
			}
		}
	}
TestEnd:

	free ((FREE_ARG) (aLikeIntList));
	free ((FREE_ARG) (aTstOfIntLst));
	return 1;
}

void DoBenchMark(int iGrouped, int iModelNbr)
{
	int n = 0;
	char acBMDL[15];
	fprintf(fp_out, "\n\n\n   Benchmark Dose Computations:");
	fprintf(fp_out, "\n\n     Specified Effect = %f", gsExpoVars.dBmdEffect);
	fflush(fp_out);

	if(gsExpoVars.iBmr_Type == 4 && (giRun[1]==1 || giRun[2]==1))
	{
		fprintf(fp_out, "\n\n            Risk Type = Extra risk; not applicable for models 2 and 3," );
		fprintf(fp_out, "\n                        but calculated for models 4 and 5." );
	}
	else
	{
		if(gsExpoVars.iBmr_Type == 0)
			fprintf(fp_out, "\n\n            Risk Type = Absolute deviation" );
		else if(gsExpoVars.iBmr_Type == 1)
		{
			if(gsExpoVars.iLogNormal == 0)
				fprintf(fp_out, "\n\n            Risk Type = Estimated standard deviations from control" );
			else
				fprintf(fp_out, "\n\n            Risk Type = Log-scale standard deviations from log control median" );
		}
		else if(gsExpoVars.iBmr_Type == 2)
			fprintf(fp_out, "\n\n            Risk Type = Relative deviation" );
		else if(gsExpoVars.iBmr_Type == 3)
			fprintf(fp_out, "\n\n            Risk Type = Point estimate" );
		else if(gsExpoVars.iBmr_Type == 4)
		{
			if(giRun[1]==1 || giRun[2]==1)
				fprintf(fp_out, "\n\n            Risk Type = Extra risk, not applicable for models 2 & 3" );
			else
				fprintf(fp_out, "\n\n            Risk Type = Extra risk" );
		}
	}

	fprintf(fp_out, "\n\n     Confidence Level = %f\n", gsExpoVars.dBmdConfi_Level);

	if(iGrouped == 0)
	{
		if(giRun[iModelNbr] == 0 || giBmdRun[iModelNbr] == 0)
		{
			if(gsExpoVars.iBmr_Type == 4 && n < 3)
				fprintf(fp_out, "\n                  BMD = %s", "Not_Applicable");
			else
				fprintf(fp_out, "\n                  BMD = %s", "Not_Computed");
		}
		else
		{
			fprintf(fp_out, "\n                  BMD = %12.6g", gdBMD[iModelNbr]);
		}
		if(giRun[iModelNbr] == 0 || gdBMDL[iModelNbr] == -1)
			if(giRun[iModelNbr] == 0)
				strcpy(acBMDL, "Not_Computed");
			else
				strcpy(acBMDL, "Bad_Completion");
		else
			sprintf(acBMDL, "%12.6g", gdBMDL[iModelNbr]);
		fprintf(fp_out, "\n\n                 BMDL = %s\n", acBMDL);
		fflush(fp_out);
	}

	if(iGrouped == 1)
	{
		OUTPUT_TEXT( "\n\n                BMD and BMDL by Model");
		OUTPUT_TEXT( "\n      Model             BMD                BMDL");
		fprintf(fp_out, "     -------        ------------        ----------");
		fflush(fp_out);
		for(n = 1; n <= NBR_OF_MODELS; n++)
		{
			if(giRun[n] == 0 || giBmdRun[n] == 0)
			{
#ifdef DO_LOG
				if (giDo_Log == true)
				{
					fprintf(fp_log,"\n*******Line 4077: INSIDE if(giRun[n] == 0 || giBmdRun[1] == 0)**********\n");
				}
#endif
				if(gsExpoVars.iBmr_Type == 4 && n < 3)
					fprintf(fp_out, "\n       %2d        %s       %s       %s", n+1, "        NA ", "            ", "Not Applicable");
				else
					fprintf(fp_out, "\n       %2d        %12.6g       %12.6g       %s", n+1, -0.00, -0.00, "Not computed");
			}
			else
			{
#ifdef DO_LOG
				if (giDo_Log == true)
				{
					fprintf(fp_log,"\n*******Line 4090: INSIDE else if(giRun[n] == 0 || giBmdRun[1] == 0)**********\n");
				}
#endif
				//fprintf(fp_out, "\n       %2d        %12.6g       %12.6g", n+1, gdBMD[n], gdBMDL[n]);
				if(gdBMDL[n] == -1)
					strcpy(acBMDL, "Bad completion");
				else
					sprintf(acBMDL, "%12.6g", gdBMDL[n]);
				fprintf(fp_out, "\n       %2d        %12.6g       %s", n+1, gdBMD[n], acBMDL);
			}
		}
		fflush(fp_out);
	}
}

int doDot002(int iStart, int iEnd, int iGrouped)
{
	//FILE *fp_002;
	int i, j, jj, nLastBS, nSrcLen;
	double lep, upep, react;
	nLastBS = 0;
	react = 0.0;
	char sPath[255];
	char modelName[255];
	char File002[255];
	strcpy(sPath,gacFileOut2);
	nSrcLen = strlen(gacFileOut2);
#ifdef DO_LOG
	if (giDo_Log == true)
	{
		fprintf(fp_log,"\n\n gagFileOut2=%s, sPath=%s, nSrcLen=%d", gacFileOut2, sPath, nSrcLen);
	}
#endif

	for(i = 0; i < nSrcLen; i++) //Test for backslash
	{
		if(gacFileOut2[i] == '\\')
			if(i > nLastBS)
				nLastBS = i;
	}

	if(nLastBS == 0) //Test for forward slash
	{
		for(i = 0; i < nSrcLen; i++)
		{
			if(gacFileOut2[i] == '/')
				if(i > nLastBS)
					nLastBS = i;
		}
	}
#ifdef DO_LOG
	if (giDo_Log == true)
	{
		fprintf(fp_log,"\n\n nLastBS = %d", nLastBS);
	}
#endif

	if(nLastBS > 0 || nSrcLen > 0)
	{
		sPath[nLastBS+(nLastBS>0?1:0)] = '\0';
		modelName[0] = 'M';
		modelName[1] = '1';
		i = nLastBS;
		j = 2;
		for(i = nLastBS+(nLastBS>0?1:0); i < nSrcLen; i++)
		{
			if(gacFileOut2[i] == '\0')
				break;
			modelName[j] = gacFileOut2[i];
			j++;
		}
		modelName[j] = '\0';
	}
	else
		goto Ending;

#ifdef DO_LOG
	if (giDo_Log == true)
	{
		fprintf(fp_log,"\n\n After BS");
	}
#endif

	if(fp_out2 != NULL)
		fclose(fp_out2);

	for( i=iStart; i <= iEnd; i++)
	{
		if(giRun[i] == 0)
			continue;

		if(fp_out2 != NULL)
		{
#ifdef DO_LOG
			if (giDo_Log == true)
			{
				fprintf(fp_log,"\n\n Before fclose");
			}
#endif
			fflush(fp_out2);
			fclose(fp_out2);
#ifdef DO_LOG
			if (giDo_Log == true)
			{
				fprintf(fp_log,"\n\n After fclose");
			}
#endif
		}

		if(i == 1)
			modelName[1] = '2';
		else if(i == 2)
			modelName[1] = '3';
		else if (i == 3)
			modelName[1] = '4';
		else if(i == 4)
			modelName[1] = '5';

#ifdef DO_LOG
		if (giDo_Log)	// Print values to log for investigation
		{
			fprintf(fp_log,"\n\n Before concat and fopen: modelName = %s", modelName);
		}
#endif
		strcpy(File002, sPath);
		strcat(File002, modelName);
		fp_out2=fopen(File002,"w");

		if(fp_out2 == NULL)
			continue;

#ifdef DO_LOG
		if (giDo_Log)	// Print values to log for investigation
		{
			fprintf(fp_log,"\n\n After concat and fopen: modelName = %s, gacFileOut2=%s", modelName, gacFileOut2);
		}
#endif

		fprintf (fp_out2, "ModelNbr \t %d\n BMD_flag \t %d \n Nobs \t%d \n nparm \t%d \n sign \t%d \n lognorm \t%d", i+1, 
			gsExpoVars.iBmdose, gsExpoVars.iNbrObs_Total, NBR_OF_PARMS, gsExpoVars.iSign, gsExpoVars.iLogNormal);
		fprintf (fp_out2, "\n Con_lev \t%3.3g ", gsExpoVars.dBmdConfi_Level);
		fprintf (fp_out2, "\n BMRT \t%d ", gsExpoVars.iBmr_Type);
		fprintf (fp_out2, "\n BMRF \t%3.3g ", gsExpoVars.dBmdEffect);

		for(j = 1; j <= NBR_OF_PARMS; j++)
			fprintf (fp_out2, "\n %s \t %5.5g", gaParm_Name[j-1], gppdMLEs[i][j]);

		fprintf (fp_out2,"\n\n Data");
		for(j = 1; j <= gsExpoVars.iNbrObs_Total; j++)
		{
			jj = gpiNi[j];
			lep = gpdYm[j] + qstudt(0.025, jj - 1.) * sqrt(gpdYd[j]/jj);
			upep= gpdYm[j] + qstudt(0.975, jj - 1.) * sqrt(gpdYd[j]/jj);
			fprintf (fp_out2,"\n %f %f %f %f", gpdXi[j], gpdYm[j], lep, upep);
		}
		fprintf (fp_out2,"\n Max_Min_dose \n  %f %f ", gdxmax, gdxmin);
		if(i == 1)
		{
			react = gppdMLEs[i][(int)ea] * exp(gsExpoVars.iSign * gppdMLEs[i][(int)eb] * gdBMD[i]);
		}
		else if(i == 2)
		{
			react = gppdMLEs[i][(int)ea] * exp(gsExpoVars.iSign * pow((gppdMLEs[i][(int)eb] * gdBMD[i]), gppdMLEs[i][(int)ed]));
		}
		else if (i == 3)
		{
			react = gppdMLEs[i][(int)ea] * (gppdMLEs[i][(int)ec] - ((gppdMLEs[i][(int)ec] - 1) * exp(-1.0 * gppdMLEs[i][(int)eb] * gdBMD[i])));
		}
		else if(i == 4)
		{
			react = gppdMLEs[i][(int)ea] * (gppdMLEs[i][(int)ec] - ((gppdMLEs[i][(int)ec] - 1) * exp(-1.0 * pow((gppdMLEs[i][(int)eb]*gdBMD[i]), gppdMLEs[i][(int)ed]))));
		}
		if(giBmdRun[i] == 0)
			react = 0;

		fprintf (fp_out2, "\n  RSL \t%f",react);
		fprintf (fp_out2, "\n  BMD \t%f", giBmdRun[i]==0?0:gdBMD[i]);

		fprintf (fp_out2,"\n\n BMD_line");
		fprintf (fp_out2,"\n %f %f", -1.0, react);
		fprintf (fp_out2,"\n %f %f", giBmdRun[i]==0?0:gdBMD[i], giBmdRun[i]==0?0:react);
		fprintf (fp_out2,"\n %f %f", giBmdRun[i]==0?0:gdBMD[i], -0.1);

		if(gdBMD[i] > 1000*gdxmax || giBmdRun[i] == 0)
			fprintf (fp_out2, "\n\n BMDL_comput_ind %d", No);  //computation failed.
		else
		{
			fprintf (fp_out2, "\n\n BMDL_comput_ind %d", Yes);
			fprintf (fp_out2, "\n\n  BMDL \t%f",gdBMDL[i]);

			fprintf (fp_out2,"\n\n BMDL_line");
			fprintf (fp_out2,"\n %f %f", gdBMDL[i], -0.1);
			fprintf (fp_out2,"\n %f %f", gdBMDL[i], react);

			fprintf (fp_out2, "\n\n BMDL_Curve_flag \t %d  \n smooth_opt  %d", gsExpoVars.iBmdlCurve, gsExpoVars.iSmooth);
		}

	} // End for(i = 1; i <= NBR_OF_MODELS; i++)

Ending:
#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\n After 002");
	}	
#endif
	return 0;
}

/*******************************************************************
*	AThree_Fit fits a likelihood to a "full" model of the form
*	Yij = Mu(i) + e(ij)  Var(eij) = k*Mu(i)^b.  The parameters Mu(i)
*	i = 1,2,...,Nobs, k, and p are estimated using the likelihood
*	maximization procedure, and the lkA3 = log-likelihood value
*	upon return.
*******************************************************************/

void AThree_Fit(int nparm, double p[], double gtol, int *iter, int nModel)
{
	int i;
	double *parms, *fitparms, *doses, ll, *svar, *means;
	double *bsv;
	double **X, **XP, **XPX, *XPY, *Y;
	long int *Spec2, *bind, *nanim, nresm, optite;
	long int restr, nparms, nvar;

	switch (gppiSpecPara[nModel][1] * 10 + gppiSpecPara[nModel][2]) 
	{
	case 11: /* Both alpha and rho fixed */
		break;
	case 1: /* Estimating alpha; rho fixed */
		/* estimate alpha as mean of Y/pow(mean, rho) */
		p[1] = 0.0;
		for (i = 1; i <= gsExpoVars.iNbrObs_Total; i++)
			p[1] += gpdYd[i] / pow(gpdYm[i], p[2]);
		p[1] = log(p[1] / (double) gsExpoVars.iNbrObs_Total);
		break;
	case 10: /* Estimate rho, alpha fixed */
		/* estimate rho as mean of (log(Var) - log(alpha)) / log(mu) */
		p[2] = 0.0;
		for (i = 1; i <= gsExpoVars.iNbrObs_Total; i++)
			p[2] += (Slog(gpdYd[i]) - Slog(p[1])) / Slog(gpdYm[i]);
		p[2] = p[2] / (double) gsExpoVars.iNbrObs_Total;
		break;
	default:  /* estimating both variance parameters */
		/* Allocate memory for linear regression */
		X = DMATRIX(1, gsExpoVars.iNbrObs_Total, 1, 2);
		XP = DMATRIX(1, 2, 1, gsExpoVars.iNbrObs_Total);
		XPX = DMATRIX(1, 2, 1, 2);
		XPY = DVECTOR(1, 2);
		bsv = DVECTOR(1, 2);
		Y = DVECTOR(1, gsExpoVars.iNbrObs_Total);

		/* set up regression problem; regress log(variance)
		on log means */
		for (i = 1; i <= gsExpoVars.iNbrObs_Total; i++)
		{
			X[i][1] = 1.0;
			X[i][2] = Slog(gpdYm[i]);
			Y[i] = Slog(gpdYd[i]);
		}
		/* linear regression */
		TRANSPOSE(X, XP, gsExpoVars.iNbrObs_Total, 2);
		MATMPYM2(XP, X, XPX, 2, gsExpoVars.iNbrObs_Total, 2);
		INVMAT(XPX, 2);
		MATMPYV2(2, gsExpoVars.iNbrObs_Total, XP, Y, XPY);
		MATMPYV2(2, 2, XPX, XPY, bsv);
		/* transfer the results to p[1:2] */
		p[1] = bsv[1];
		p[2] = bsv[2];
		/* Free memory */
		FREE_DMATRIX(X, 1, gsExpoVars.iNbrObs_Total, 1, 2);
		FREE_DMATRIX(XP, 1, 2, 1, gsExpoVars.iNbrObs_Total);
		FREE_DMATRIX(XPX, 1, 2, 1, 2);
		FREE_DVECTOR(XPY, 1, 2);
		FREE_DVECTOR(bsv, 1, 2);
		FREE_DVECTOR(Y, 1, gsExpoVars.iNbrObs_Total);
		break;
	}

	for (i = 1; i <= nparm-2; i++)
	{
		p[i + 2] = gpdYm[i];
	}

	nvar = gsExpoVars.iNbrObs_Total;
	restr = 1;
	nparms = nparm;

	doses = DVECTOR(0, nvar-1);
	means = DVECTOR(0, nvar-1);
	svar =  DVECTOR(0, nvar-1);
	nanim = LIVECTOR(0, nvar-1);
	parms = DVECTOR(0, nparm-1);
	fitparms = DVECTOR(0, nparm-1);
	Spec2 = LIVECTOR(0, nparm-1);
	bind = LIVECTOR(0, nparm-1);

	for (i = 1; i <= nvar; i++)
	{
		nanim[i - 1] = (long int) gpiNi[i];
		doses[i - 1] = gpdXi[i];
		means[i - 1] = gpdYm[i];
		svar[i - 1] = gpdYd[i];
	}

	for (i = 1; i <= nparm; i++)
	{
		Spec2[i - 1] = gppiSpecPara[nModel][1];
		parms[i - 1] = p[i];
	}

	getmlea3_(&nvar, doses, means, nanim, svar, &nparms, parms,
		Spec2, parms, &restr, fitparms, &ll, &optite,
		&nresm, bind);

	for (i = 1; i <= nparm; i++)
	{
		p[i] = fitparms[i - 1];
	}

	//*fret = ll;
	gpdA3s[nModel] = ll;

	FREE_DVECTOR(doses, 0, gsExpoVars.iNbrObs_Total-1);
	FREE_DVECTOR(means, 0, gsExpoVars.iNbrObs_Total-1);
	FREE_DVECTOR(svar, 0, gsExpoVars.iNbrObs_Total-1);
	FREE_LIVECTOR(nanim, 0, gsExpoVars.iNbrObs_Total-1);
	FREE_DVECTOR(parms, 0, nparm-1);
	FREE_DVECTOR(fitparms, 0, nparm-1);
	FREE_LIVECTOR(Spec2, 0, nparm-1);
	FREE_LIVECTOR(bind, 0, nparm-1);

}	/* end AThree_Fit */


/*******************************************************************
 *	First partial derivatives of the mean function for model mod
 *      at observation obs, contained in mg[] at exit.
 *******************************************************************/
void
MeanPart (int obs, double *mg, int mod)
{
  int i;
  double *p;
  double sign, dose;

  p = DVECTOR(1, NBR_OF_PARMS);

  sign = gsExpoVars.iSign;
  dose = gpdXi[obs];
  // copy parms to temp array
  for (i=1; i<=NBR_OF_PARMS; i++)
      p[i] = gppdMLEs[mod][i];

  mg[1] = mg[2] = 0.0;		/* Variance parameters not in mean function */

  switch(mod)
      {
        case (int)eM2:
	    mg[3] = exp(sign*p[4]*dose);
	    mg[4] = p[3]*dose*sign*exp(sign*p[4]*dose);
	    mg[5] = 0;
	    mg[6] = 0;
	    break;
	case (int)eM3:
  	    mg[3] = exp(sign*pow(p[4]*dose, p[6]));
	    mg[4] = p[3]*p[6]*pow(p[4]*dose,p[6]-1)*dose*sign*exp(sign*pow(p[4]*dose,p[6]));
	    mg[5] = 0;
            if (dose==0)
              mg[6] = 0;
            else
	      mg[6] = p[3]*pow(p[4]*dose,p[6])*sign*log(p[4]*dose)*exp(pow(p[4]*dose,p[6]));
            break;
	case (int)eM4:
	    mg[3] = p[5] - (p[5]-1)*exp(-1*p[4]*dose);
	    mg[4] = p[3]*(p[5]-1)*dose*exp(-1*p[4]*dose);
	    mg[5] = p[3]*(1-exp(-1*p[4]*dose));
	    mg[6] = 0;
            break;
	case (int)eM5:
	    mg[3] = p[5] - (p[5]-1)*exp(-1* pow(p[4]*dose, p[6]));
	    mg[4] = p[3]*(p[5]-1)*p[6]*dose*pow(p[4]*dose,p[6]-1)*exp(-1*pow(p[4]*dose,p[6]));
	    mg[5] = p[3]*(1-exp(-1*pow(p[4]*dose,p[6])));
            if (dose == 0)
              mg[6] = 0;
            else
	      mg[6] = p[3]*(p[5]-1)*pow(p[4]*dose,p[6])*log(p[4]*dose)*exp(-1*pow(p[4]*dose,p[6]));
            break;
      }

}

/********************************************************************
 *	Second partial derivates with respect to the mean function for model mod.
 *	mg2[][] contains all second partials upon exit.
 ********************************************************************/
void Mean2Part (int obs, double **mg2, int mod)
{
  int i, j;
  double *p;
  double sign, dose, bdd;

  p = DVECTOR(1, NBR_OF_PARMS);

  sign = gsExpoVars.iSign;
  dose = gpdXi[obs];
  // copy parms to temp array
  for (i=1; i<=NBR_OF_PARMS; i++)
      p[i] = gppdMLEs[mod][i];

  for (i=1; i<= NBR_OF_PARMS; i++)
      for (j=1; j<= NBR_OF_PARMS; j++)
        mg2[i][j] = 0;

  bdd = pow(p[4]*dose,p[6]);

  switch(mod)
      {
        case (int)eM2:
	    mg2[4][4] = p[3]*dose*dose*sign*exp(p[4]*dose*sign);
	    mg2[3][4] = dose*sign*exp(p[4]*dose*sign);
            break;
	case (int)eM3:
	    mg2[4][4] = sign*p[3]*p[6]*bdd*exp(sign*bdd)*(-1+p[6]+p[6]*bdd*sign)/pow(p[4],2);
	    mg2[3][4] = sign*p[6]*dose*pow(p[4]*dose,p[6]-1)*exp(sign*bdd);

            if (dose != 0.0)
              {
                mg2[6][6] = sign*p[3]*bdd*exp(sign*bdd)*(1+sign*bdd)*log(pow(p[4]*dose,2));
		mg2[3][6] = sign*bdd*exp(sign*bdd)*log(p[4]*dose);
		mg2[4][6] = sign*p[3]/p[4]*bdd*exp(sign*bdd)*(1+p[6]*(1+bdd*sign)*log(p[4]*dose));
              }
            break;
	case (int)eM4:
	    mg2[4][4] = -p[3]*(p[5]-1)*pow(dose,2)*exp(-1*p[4]*dose);
	    mg2[3][4] = (p[5]-1)*dose*exp(-1*p[4]*dose);
	    mg2[3][5] = 1-exp(-1*p[4]*dose);
            mg2[4][5] = p[3]*dose*exp(-1*p[4]*dose);
            break;
	case (int)eM5:
	    mg2[4][4] = -p[3]/pow(p[4],2)*(p[5]-1)*p[6]*bdd*(1+p[6]*(-1+bdd))*exp(-1*bdd);
	    mg2[3][4] = (p[5]-1)*p[6]*dose*pow(p[4]*dose,p[6]-1)*exp(-1*bdd);
	    mg2[3][5] = 1-exp(-1*bdd);
            mg2[4][5] = p[3]*p[6]*dose*pow(p[4]*dose,p[6]-1)*exp(-1*bdd);
            if (dose != 0.0)
              {
		mg2[6][6] = -p[3]*(p[5]-1)*bdd*(bdd-1)*exp(-1*bdd)*pow(log(p[4]*dose),2);
		mg2[3][6] = (p[5]-1)*bdd*exp(-1*bdd)*log(p[4]*dose);
		mg2[4][6] = -p[3]/p[4]*(p[5]-1)*bdd*exp(-1*bdd)*(-1+p[6]*(bdd-1)*log(p[4]*dose));
            	mg2[5][6] = p[3]*bdd*exp(-1*bdd)*log(p[4]*dose);
              }
            break;
      }



  /* make mg2 symmetric */
  for (i=1; i<=NBR_OF_PARMS; i++)
    for (j=i; j<=NBR_OF_PARMS; j++)
      mg2[j][i] = mg2[i][j];


}


/*********************************************************************
 *	First partial derivatives of the variance function.  Vi and meani
 *	are the estimated mean and variance respectively, and should be
 *	passed by the calling function.  const_var = 1 if variance is
 *	constant, = 0 if not.  mg[] are the first partials of the mean
 *	function and should be passed by the calling function.  Partials
 *	stored in vg[] upon exit.
 *********************************************************************/
void VarPart (int obs, double *mg, double *vg, int mod)
{

  double estVar, meani, dose;
  double *p;
  int i, sign;

  p = DVECTOR(1, NBR_OF_PARMS);

  sign = gsExpoVars.iSign;
  dose = gpdXi[obs];
  // copy parms to temp array
  for (i=1; i<=NBR_OF_PARMS; i++)
      p[i] = gppdMLEs[mod][i];

  switch(mod)
      {
        case 1:
		meani = p[3]*exp(sign*p[4]*dose);
		break;
	case (int)eM3:
		meani = p[3]*exp(sign*pow(p[4]*dose,p[6])); 
		break;
	case (int)eM4:
		meani = p[3]*(p[5]-(p[5]-1)*exp(-1*p[4]*dose)); 
		break;
	case (int)eM5:
		meani = p[3]*(p[5]-(p[5]-1)*exp(-1*pow(p[4]*dose,p[6])));
		break;
      }


  if (gsExpoVars.iCons_Var == 1)
      estVar = exp(p[1]);
  else
      estVar = exp(p[1]+log(fabs(meani))*p[2]);

  for (i=1; i<=NBR_OF_PARMS; i++)
      vg[i] = 0;

  if (gsExpoVars.iCons_Var == 1)  /* constant variance */
    {
      vg[1] = 1.0;
      vg[2] = 0.0;
    }
  else  /* non-constant variance */
    {
      vg[1] = estVar;
      if (meani == 0)
        vg[2] = 0;
      else 
	vg[2] = estVar*log(fabs(meani));
    }

  for (i=3; i<=NBR_OF_PARMS; i++)
    {
    if (fabs (meani) > 1e-20)
	vg[i] = p[2] * estVar * mg[i] / fabs (meani);
      else
	vg[i] = 0.0;
    }


}




/*********************************************************************
 *	Second partial derivatives of the variance function.  Vi and meani
 *	are the estimated mean and variance respectively, and should be
 *	passed by the calling function.  const_var = 1 if variance is
 *	constant, = 0 if not.  mg[] are the first partials of the mean
 *	function, mg2[][] are the second partials of the mean function.
 *	Both should be passed by the calling function.  Partials
 *	stored in vg2[][] upon exit.
 *********************************************************************/
void Var2Part (int obs, double *mg, double **mg2, double **vg2, int mod)
{
  double logam, abmn, temp, meani, dose, estVar;
  double *p;
  int i, j, k, Sign, sign;

  p = DVECTOR(1, NBR_OF_PARMS);

  sign = gsExpoVars.iSign;
  dose = gpdXi[obs];
  // copy parms to temp array
  for (i=1; i<=NBR_OF_PARMS; i++)
      p[i] = gppdMLEs[mod][i];


  switch(mod)
      {
        case (int)eM2:
		meani = p[3]*exp(sign*p[4]*dose);
		break;
	case (int)eM3:
		meani = p[3]*exp(sign*pow(p[4]*dose,p[6])); 
		break;
	case (int)eM4:
		meani = p[3]*(p[5]-(p[5]-1)*exp(-1*p[4]*dose)); 
		break;
	case (int)eM5:
		meani = p[3]*(p[5]-(p[5]-1)*exp(-1*pow(p[4]*dose,p[6]))); 
		break;
      }

  abmn = fabs (meani);
  logam = Slog (abmn);

  if (gsExpoVars.iCons_Var == 1)
      estVar = exp(p[1]);
  else
      estVar = exp(p[1]+log(fabs(meani))*p[2]);


  if (gsExpoVars.iCons_Var == 1)
    {
      /* constant variance.  Vi = alpha.  All second partials = 0 */
      for (j = 1; j <= NBR_OF_PARMS; j++)
	{
	  for (k = 1; k <= NBR_OF_PARMS; k++)
	    {
	      vg2[j][k] = 0.0;
	    }
	}
    }
  else
    {
      /* non constant variance.  estVar = alpha*|Meani|**rho */

      if (meani < 0)
	{
	  Sign = -1;
	}
      else
	{
	  Sign = 1;
	}

      vg2[1][1] = estVar;
      vg2[1][2] = estVar * logam;
      vg2[2][1] = vg2[1][2];

      for (j = 3; j <= NBR_OF_PARMS; j++)
	{
	  vg2[1][j] = Sign * p[2] * estVar * mg[j] / abmn;
	  vg2[j][1] = vg2[1][j];
	}
      vg2[2][2] = estVar * logam * logam;

      for (j = 3; j <= NBR_OF_PARMS; j++)
	{
	  vg2[2][j] = Sign * estVar * mg[j] * ((p[2] * logam) + 1) / abmn;
	  vg2[j][2] = vg2[2][j];
	}

      for (j = 3; j <= NBR_OF_PARMS; j++)
	{
	  for (k = j; k <= NBR_OF_PARMS; k++)
	    {
	      temp = ((p[2] - 1) * mg[j] * mg[k] / abmn) + (Sign * mg2[j][k]);
	      vg2[j][k] = p[2] * estVar * temp / abmn;
	      vg2[k][j] = vg2[j][k];
	    }
	}
    }				/* end if (const_var == 1) */

    /* make mg2 symmetric */
  for (i=1; i<=NBR_OF_PARMS; i++)
    for (j=i; j<=NBR_OF_PARMS; j++)
      vg2[j][i] = vg2[i][j];

}				/* end Var2Part */



void F1iDoublePart( double **Fn1i, int obs, int mod)
/***********************************************************
Calculates the second partial derivatives of the function
F = ln(Vi) at dose[obs] where Vi is the estimated variance.  Fn1i[j][k]
contains the second partial with respect to parameters j and
k.
**********************************************************/
{
	double *mg, *vg, **mg2, **vg2, *p, estVar;
	double meani, temp, dose;
	int i, j, k, sign;


	for(j = 1; j <= NBR_OF_PARMS; j++)
		for(k = 1; k <= NBR_OF_PARMS; k++)
			Fn1i[j][k] = 0.0;

	mg = DVECTOR(1,NBR_OF_PARMS);
	vg = DVECTOR(1,NBR_OF_PARMS);
	mg2 = DMATRIX(1,NBR_OF_PARMS,1,NBR_OF_PARMS);
	vg2 = DMATRIX(1,NBR_OF_PARMS,1,NBR_OF_PARMS);
        p = DVECTOR(1, NBR_OF_PARMS);

  	sign = gsExpoVars.iSign;
  	dose = gpdXi[obs];
  	// copy parms to temp array
  	for (i=1; i<=NBR_OF_PARMS; i++)
      	  p[i] = gppdMLEs[mod][i];

	/* Compute the estimated mean */

	switch(mod)
        {
          case (int)eM2:
		meani = p[3]*exp(sign*p[4]*dose);
		break;
	  case (int)eM3:
		meani = p[3]*exp(sign*pow(p[4]*dose,p[6])); 
		break;
	  case (int)eM4:
		meani = p[3]*(p[5]-(p[5]-1)*exp(-1*p[4]*dose)); 
		break;
	  case (int)eM5:
		meani = p[3]*(p[5]-(p[5]-1)*exp(-1*pow(p[4]*dose,p[6]))); 
		break;
        }

        if (gsExpoVars.iCons_Var == 1)
       	    estVar = exp(p[1]);
  	else
      	    estVar = exp(p[1]+log(fabs(meani))*p[2]);
		     



	if (estVar == 0) estVar = .000001;
	/* Get the partial derivatives of the mean function at dose[obs] */
	MeanPart(obs, mg, mod);
	/* Get the partial derivatives of the variance function */
	VarPart(obs, mg, vg, mod);
	/* Get second partials of the mean function */
	Mean2Part(obs, mg2, mod); 
	/* Get second partials of the variance function */
	Var2Part(obs, mg, mg2, vg2, mod);




	/* Calculate partial derivative at dose[obs] */

	for(j = 1; j <= NBR_OF_PARMS; j++)
	{
		for(k = j; k <= NBR_OF_PARMS; k++)
		{
			temp = vg2[j][k] - vg[j]*vg[k]/estVar ;
			Fn1i[j][k] = temp/estVar;
			Fn1i[k][j] = Fn1i[j][k];
		}
	}

	FREE_DVECTOR(mg, 1, NBR_OF_PARMS);
	FREE_DVECTOR(vg, 1, NBR_OF_PARMS);
	FREE_DMATRIX(mg2,1,NBR_OF_PARMS,1,NBR_OF_PARMS);
	FREE_DMATRIX(vg2,1,NBR_OF_PARMS,1,NBR_OF_PARMS);

}


/***********************************************************
 *	Calculates the second partial derivatives of the function
 *	F = 1/Vi at dose[obs] where Vi is the estimated variance.  Fn2i[j][k]
 *	contains the second partial with respect to parameters j and
 *	k.
 ***********************************************************/
void F2iDoublePart (double **Fn2i, int obs, int mod)
{
  double *mg, *vg, **mg2, **vg2, *p, estVar;
  double meani, temp, dose;
  int i, j, k, sign;


  for (j = 1; j <= NBR_OF_PARMS; j++)
    {
      for (k = 1; k <= NBR_OF_PARMS; k++)
	{
	  Fn2i[j][k] = 0.0;
	}
    }

  mg = DVECTOR (1, NBR_OF_PARMS);
  vg = DVECTOR (1, NBR_OF_PARMS);
  mg2 = DMATRIX (1, NBR_OF_PARMS, 1, NBR_OF_PARMS);
  vg2 = DMATRIX (1, NBR_OF_PARMS, 1, NBR_OF_PARMS);


  /* Compute the estimated mean */

  p = DVECTOR(1, NBR_OF_PARMS);

  sign = gsExpoVars.iSign;
  dose = gpdXi[obs];
  // copy parms to temp array
  for (i=1; i<=NBR_OF_PARMS; i++)
    p[i] = gppdMLEs[mod][i];

  /* Compute the estimated mean */

  switch(mod)
  {
    case (int)eM2:
	meani = p[3]*exp(sign*p[4]*dose);
	break;
    case (int)eM3:
	meani = p[3]*exp(sign*pow(p[4]*dose,p[6])); 
	break;
    case (int)eM4:
	meani = p[3]*(p[5]-(p[5]-1)*exp(-1*p[4]*dose)); 
	break;
    case (int)eM5:
	meani = p[3]*(p[5]-(p[5]-1)*exp(-1*pow(p[4]*dose,p[6]))); 
	break;
  }

  if (gsExpoVars.iCons_Var == 1)
        estVar = exp(p[1]);
  else
        estVar = exp(p[1]+log(fabs(meani))*p[2]);

  /* Get the partial derivatives of the mean function at dose[obs] */
  /* Get the partial derivatives of the mean function at dose[obs] */
  MeanPart(obs, mg, mod);
  /* Get the partial derivatives of the variance function */
  VarPart(obs, mg, vg, mod);
  /* Get second partials of the mean function */
  Mean2Part(obs, mg2, mod); 
  /* Get second partials of the variance function */
  Var2Part(obs, mg, mg2, vg2, mod);


  /* Calculate partial derivative at dose[obs] */

  for (j = 1; j <= NBR_OF_PARMS; j++)
    {
      for (k = j; k <= NBR_OF_PARMS; k++)
	{
	  temp = (2 * vg[j] * vg[k] / estVar) - vg2[j][k];
	  Fn2i[j][k] = temp / (estVar * estVar);
	  Fn2i[k][j] = Fn2i[j][k];

	}
    }


  FREE_DVECTOR (mg, 1, NBR_OF_PARMS);
  FREE_DVECTOR (vg, 1, NBR_OF_PARMS);
  FREE_DMATRIX (mg2, 1, NBR_OF_PARMS, 1, NBR_OF_PARMS);
  FREE_DMATRIX (vg2, 1, NBR_OF_PARMS, 1, NBR_OF_PARMS);


}				/* end F2iDoublePart */


/***********************************************************
 *	Calculates the second partial derivatives of the function
 *	F = (Ybar - Mi)**2/Vi at dose[obs] where Vi is the estimated variance
 *	Ybar is the sample mean, and Mi is the estimated mean.
 *	Fn1i[j][k]
 *	contains the second partial with respect to parameters j and
 *	k.
 ***********************************************************/
void F3iDoublePart (double **Fn3i, int obs, int mod)
{
  double *mg, *vg, **mg2, **vg2, *p, estVar;
  double Devi, meani, temp, temp2, temp3, dose;
  int i, j, k, sign;


  for (j = 1; j <= NBR_OF_PARMS; j++)
    {
      for (k = 1; k <= NBR_OF_PARMS; k++)
	{
	  Fn3i[j][k] = 0.0;
	}
    }

  mg = DVECTOR (1, NBR_OF_PARMS);
  vg = DVECTOR (1, NBR_OF_PARMS);
  mg2 = DMATRIX (1, NBR_OF_PARMS, 1, NBR_OF_PARMS);
  vg2 = DMATRIX (1, NBR_OF_PARMS, 1, NBR_OF_PARMS);


  /* Compute the estimated mean */

  p = DVECTOR(1, NBR_OF_PARMS);

  sign = gsExpoVars.iSign;
  dose = gpdXi[obs];
  // copy parms to temp array
  for (i=1; i<=NBR_OF_PARMS; i++)
    p[i] = gppdMLEs[mod][i];

  /* Compute the estimated mean */

  switch(mod)
  {
    case (int)eM2:
	meani = p[3]*exp(sign*p[4]*dose);
        break;
    case (int)eM3:
	meani = p[3]*exp(sign*pow(p[4]*dose,p[6])); 
        break;
    case (int)eM4:
	meani = p[3]*(p[5]-(p[5]-1)*exp(-1*p[4]*dose)); 
        break;
    case (int)eM5:
	meani = p[3]*(p[5]-(p[5]-1)*exp(-1*pow(p[4]*dose,p[6]))); 
        break;
  }

  Devi = gpdYm[obs] - meani;


  if (gsExpoVars.iCons_Var == 1)
        estVar = exp(p[1]);
  else
        estVar = exp(p[1]+log(fabs(meani))*p[2]);

  /* Get the partial derivatives of the mean function at dose[obs] */
  /* Get the partial derivatives of the mean function at dose[obs] */
  MeanPart(obs, mg, mod);
  /* Get the partial derivatives of the variance function */
  VarPart(obs, mg, vg, mod);
  /* Get second partials of the mean function */
  Mean2Part(obs, mg2, mod); 
  /* Get second partials of the variance function */
  Var2Part(obs, mg, mg2, vg2, mod);


  /* Calculate partial derivative at dose[obs] */


  for (j = 1; j <= NBR_OF_PARMS; j++)
    {
      for (k = j; k <= NBR_OF_PARMS; k++)
	{
	  temp = 2 * estVar * estVar * (mg[j] * mg[k] - Devi * mg2[j][k]);
	  temp2 = 2 * Devi * estVar * (vg[j] * mg[k] + mg[j] * vg[k]);
	  temp3 = Devi * Devi * (2 * vg[j] * vg[k] - estVar * vg2[j][k]);
	  Fn3i[j][k] = (temp + temp2 + temp3) / (estVar * estVar * estVar);
	  Fn3i[k][j] = Fn3i[j][k];
	}
    }


  FREE_DVECTOR (mg, 1, NBR_OF_PARMS);
  FREE_DVECTOR (vg, 1, NBR_OF_PARMS);
  FREE_DMATRIX (mg2, 1, NBR_OF_PARMS, 1, NBR_OF_PARMS);
  FREE_DMATRIX (vg2, 1, NBR_OF_PARMS, 1, NBR_OF_PARMS);
  
}




/******************************************************************
 *	Exp_vcv -- used to compute the vcv for Exponential model.
 *		Extern var.: smean, smax, Nobs, Xi, Yp, Yn, Ls, Xg.
 *
 ******************************************************************/
void Exp_vcv (double ***gppdVCV, int mod)
{
  void F1iDoublePart (double **Fn1i, int obs, int mod);

  void F2iDoublePart (double **Fn2i, int obs, int mod);

  void F3iDoublePart (double **Fn3i, int obs, int mod);

  void MeanPart (int obs, double *mg, int mod);
  void VarPart (int obs, double *mg, double *vg, int mod);
  void Mean2Part (int obs, double **mg2, int mod);
  void Var2Part (int obs, double *mg, double **mg2, double **vg2, int mod);

  double **Fn1i, **Fn2i, **Fn3i, *p, numi;
  int i, j, k;

  Fn1i = DMATRIX (1, NBR_OF_PARMS, 1, NBR_OF_PARMS);
  Fn2i = DMATRIX (1, NBR_OF_PARMS, 1, NBR_OF_PARMS);
  Fn3i = DMATRIX (1, NBR_OF_PARMS, 1, NBR_OF_PARMS);
  p = DVECTOR(1, NBR_OF_PARMS);

  /* Compute partials at parameter estimates */

  // copy parms to temp array
  for (i=1; i<=NBR_OF_PARMS; i++)
    p[i] = gppdMLEs[mod][i];

  for (j = 1; j <= NBR_OF_PARMS; j++)
    {
      for (k = 1; k <= NBR_OF_PARMS; k++)
        {
          gppdVCV[mod-1][j][k] = 0.0;
        }
    }



  for (i = 1; i <= gsExpoVars.iNbrObs_Total; i++)
    {
      numi = gpiNi[i];
      for (j = 1; j <= NBR_OF_PARMS; j++)
	{
	  for (k = 1; k <= NBR_OF_PARMS; k++)
	    {
	      F1iDoublePart (Fn1i, i, mod);
	      F2iDoublePart (Fn2i, i, mod);
	      F3iDoublePart (Fn3i, i, mod);

              gppdVCV[mod-1][j][k] = gppdVCV[mod-1][j][k] + (numi * Fn1i[j][k] / 2);
	      gppdVCV[mod-1][j][k] += (numi - 1) * gpdYd[i] * Fn2i[j][k] / 2;
	      gppdVCV[mod-1][j][k] += numi * Fn3i[j][k] / 2;
	    }
	}
    }


  FREE_DMATRIX (Fn1i, 1, NBR_OF_PARMS, 1, NBR_OF_PARMS);
  FREE_DMATRIX (Fn2i, 1, NBR_OF_PARMS, 1, NBR_OF_PARMS);
  FREE_DMATRIX (Fn3i, 1, NBR_OF_PARMS, 1, NBR_OF_PARMS);
  FREE_DVECTOR (p, 1, NBR_OF_PARMS);

}				/* end Exp_vcv */





/*********************************************************************************
*
*  GET_DTMSVCV--This function will take the matrix of Likelihood
*		    second derivatives and remove the columns/rows
*		    that correspond to any bounded or fixed parameters,
*		    invert this modified matrix.
*********************************************************************************/
void Get_DTMSVCV (int mod)
{
  int adj_rows, i, nparm, nboundparm;
  double **temp_vcv; 


  /* Create an adjusted variance-covariance matrix that does not
  /  have correlations/variances for parameters that are estimated
  /  at a boundary point */

  nparm = NBR_OF_PARMS;
  temp_vcv = DMATRIX (1, nparm, 1, nparm);
  temp_vcv = gppdVCV[mod-1];

  nboundparm = 0;
  for (i = 1; i <= nparm; i++)
    {
      if (bounded[mod][i] == 1)
	{
	  nboundparm++;

	}			/* end if */

    }				/* end for */


  /* Adjust matrix to only contain non-fixed or bounded parameters */
  adj_rows = Take_Out_Bounded_Parms (nparm, bounded[mod], temp_vcv);
  gppdVCV_adj[mod-1] = temp_vcv;

  /* Invert to get the covariance matrix */
  INVMAT (gppdVCV_adj[mod-1], adj_rows);

  FREE_DMATRIX (temp_vcv, 1, nparm, 1, nparm);
    
}






/*********************************************************
** Calc_ParmsE--output model fitting parameter SEs 
*********************************************************/
 void Calc_ParmSE (int mod)
{
	int i, num_bound, pos_def;


	// Check to make sure variances are all positive; if not, only output estimates, with a stern warning 
	pos_def = 1;
	num_bound = 0;
	for (i = 1; i <= NBR_OF_PARMS; i++) 
	{
		if (bounded[mod][i] == 0) 
		{
			if (gppdVCV_adj[mod-1][i-num_bound][i-num_bound] <= 0) 
			{
				pos_def = 0;
				break;
			}
		} 
		else 
		{
			num_bound++;
		}
	}

	

	if (pos_def == 1)
	{
		num_bound = 0;
        	for (i=1; i<=NBR_OF_PARMS; i++)
		{
			
            		if (bounded[mod][i] == 0)
			{
				parmSE[mod][i] = sqrt(fabs(gppdVCV_adj[mod-1][i - num_bound][i - num_bound] ));
			}
			else
			{
				if (gppiSpecPara[mod][i] == 0)
				{
					parmSE[mod][i] = -9999; // indicates that NA should be printed
				}
				else
				{
					parmSE[mod][i] = -9999; // indicates that NA should be printed
					num_bound++;
				}
			}
		}
	}
	else
	{
		for (i=1; i<=NBR_OF_PARMS; i++)
		{
			parmSE[mod][i] = -9999;	// indicates that NA should be printed
		}
	}        
        
}

