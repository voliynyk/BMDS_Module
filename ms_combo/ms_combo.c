/******************************************************************
* IMPORTANT NOTE:  The following variable is the version number for
*                  the current model.  THIS MUST BE CHANGED as
*		   important changes are made to the models.
*
*****************************************************************/
char Version_no[]="MS_COMBO. (Version: 1.9;  Date: 05/20/2014)";

/****************************************************************
** MS_COMBO is built on Multistage version 2.8
* 
* MS_COMBO - a ANSI C program for fitting the Multistage model 
*            fitting with/without a natural background rate in Benchmark Dose
*            to two seperate endpoints and then doing a combination of the
*            endpoints to get a BMDL
* Start Date: July 6, 2007
*
********************************************************************
* Modification Log for MS_comboo.c:
*
* Version Number: 1.0
* Modified By: Cynthia Van Landingham 
* Date: July - Aug 2007
* Reason: Implementaion of new method for combining tumors
******************************************************************
* Version Number: 1.1
* Modified By: Geoffrey Nonato 
* Date: 02/12/2009
* Reason: Implementation of multiple tumors (greater than 2)
******************************************************************
* Version Number: 1.2
* Modified By: Geoffrey Nonato 
* Date: 04/30/2009
* Reason: Recompiled fortran code for format correction
******************************************************************
* Version Number: 1.3
* Modified By: Louis Olszyk
* Date: 07/22/2010
* Reason: Resolved issues with >2 tumors and added debug logging,
*         which can be enabled at compile time.
******************************************************************
* Version Number: 1.4
* Modified By: Louis Olszyk
* Date: 10/20/2010
* Reason: Resolved problem calculating combo BMDL when the maximum
*         dose differs across all tumors.
******************************************************************
* Version Number: 1.5
* Modified By: Geoffrey Nonato
* Date: 01/25/2011
* Reason: Changed version to "Beta"
******************************************************************
* Version Number: 1.6 Beta
* Modified By: Louis Olszyk
* Date: 02/28/2013
* Reason: - PR 460: Removed unnecessary if statement that prevent the Cancer
*           Slope Factor from being printed, and added CSL for the combo bmdl.
*         - PR 444: Fixed plot titles
*         - PR 463: Display dataset name for each dose in output
*         - PR 470: Resolve combo bmdl issue that occurs from inconsistent
*                   scaling when the first dataset does not contain the
*                   maximum dose across all of the datasets.
******************************************************************
* Version Number: 1.7 Beta
* Modified By: Louis Olszyk
* Date: 03/25/2014
* Reason: - PR 494: Fixed memory corruption resulting from a call to
*           strncpy() that used the size of the source buffer (512)
*           instead of the target buffer (64).
*                   
******************************************************************
* Version Number: 1.8 Beta
* Modified By: Louis Olszyk
* Modification Date: 04/30/2014
* Reason: PRs 439,486 - Allow non-integer values for N and #affected
******************************************************************
* Version Number: 1.9
* Modified By: Louis Olszyk
* Modification Date: 05/20/2014
* Reason: PRs 355, 357 - Fixed convergence problem when betas are specified.
******************************************************************/

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



extern void getmle_(double parms[], long int fixed[], double fixedval[],
					double parms2[], double *ll,
					long int *optite, long int *nresm, long int bind[]);


extern void getcl_(long int *which, double *bmr, double *bmd,
				   double *target, double parms[], long int fixed[],
				   double fixedval[], long int *risktype, double *bmdl,
				   double parms2[], long int *optite, long int *nresm,
				   long int bind[], int *);

//extern void getcl2_(long int *which, double *bmr, double *bmd,
//					double *target, double parms[], double aparms[], long int fixed[],
//					double fixedval[],  long int afixed[],
//					double afixedval[], long int *risktype, double *bmdl,
//					double parms2[], long int *optite, long int *nresm,
//					long int bind[], int *, long int, double *,double *,
//					double *,double *);

extern void loadcommbloc_(long int *ndoses, double *xmax,
			  double nanimals[],double doses[],
			  double affect[], long int *polyord,
			  long int *restrict, double *, double *, int *n);
/** cvl added 8/2007 - start block **/
//extern void loadcommbl2_(long int *which);
/** cvl added 8/2007 - end block **/


#define EPS 3.0e-8
#define float double

void Multistage_fit(int nparm, double p[], double gtol,
					int *iter, double *fret);
void Multistage_BMD(int nparm, double p[], double gtol, int *iter, double xlk,
					double Rlevel[], double Bmdl[],double Bmdu[],double *BMD);
void Multistage_vcv(int nparm, int Spec[], double p[], double **vcv);
void GetNewParms(double *p, int size);
void GetMoreParms(double *p, int size);
void Copy_dvec(double *a,double *b, int size);
void GetMLEParms(double *p, int size);

/** cvl added 8/2007 - start block **/
void Multistage_ComboBMD (double gtol, int *iter, double xlk,
						  double Rlevel[], double Bmdl[], double Bmdu[], double *BMD);
double LogLik_Constant(int Nobs,  double Yp[], double Yn[]);
double    ComboMaxLike(int , double, double *);
/** cvl added 8/2007 - end block **/

double P_RESPONSE(double q, double p[]);
int Model_DF (int []);

/*** Define input and output files's name  *********************/
char     fin[128];  /* input temp file */
char    fout[128];  /* output temp file */
char    fout2[128]; /* output file for plot */
char    fout3[128]; /* output file for logging */

/* Store Logging/debug level */
#define LOG_NONE 0
#define LOG_DEVELOPER  0x00008000
/* We currently enable/disable debug logging at compile time. */
#ifdef MISC_OUT
#define LOG_LEVEL LOG_DEVELOPER
#else
#define LOG_LEVEL LOG_NONE
#endif /* MISC_OUT */
/* Use a structure so that it maps to a fortran common block */
struct {
  int geloglevel;
} _loginfo = {0x00008000};	/* value was LOG_LEVEL */

char *Parm_name[] ={"Background", "Beta(1)", "Beta(2)", "Beta(3)",
"Beta(4)",  "Beta(5)",  "Beta(6)",  "Beta(7)",
"Beta(8)",  "Beta(9)",  "Beta(10)", "Beta(11)",
"Beta(12)", "Beta(13)", "Beta(14)", "Beta(15)",
"Beta(16)", "Beta(17)", "Beta(18)", "Beta(19)",
"Beta(20)", "Beta(21)", "Beta(22)", "Beta(23)",
"Beta(24)", "Beta(25)", "Beta(26)", "Beta(27)",
"Beta(28)", "Beta(29)", "Beta(30)", "Beta(31)",
"Beta(32)", "Beta(33)", "Beta(34)", "Beta(35)",
"Beta(36)", "Beta(37)", "Beta(38)", "Beta(39)",
"Beta(40)", "Beta(41)", "Beta(42)"};


char plotfilename[128];	/* plot file name */
char *anatxt[]={"Full model", "Fitted model", "Reduced model"};


/*** variables will not be changed except Spec  *******/
/********************************************************************
* CVL - 8/3/2007
* these varibles will be use once for each model and then 
*  with a combination style for the final BMD, BMDL estimates 
*********************************************************************/
int ITMAX;         /* maximum # of iterations */
int    *Spec;      /* vector used to identify user input parm. */
/* Spec[i]=1 if parameter i was specified */
int    *IniSp;     /* vector used to identify initialized parameters */
/* IniSp[i]=1 if parameter i was initialized */
double *Yp;        /* positive dependent variable data array */
double *Yn;        /* negative dependent variable data array */
double *Xi;        /* independent variable data array */
double *Ypp;       /* predicted positive dependent values */
double *Ep;        /* estimated probability at each dose */
double *Rlevel;    /* response levels, .05, .1, .2 .3 */
double *Bmdl;      /* Bmdl values at the response levels */
double *Bmdu;      /* Bmdu values at the response levels */
double *IniP;      /* initial parameter values */
long int Nobs;     /* # of observations */
/** cvl added 8/8/2007 - start block **/
//long int Anobs;     /* Model A number of obs */
//double *AYp;        /* model A positive dependent variable data array */
//double *AYn;        /* model A negative dependent variable data array */
//double *AXi;        /* model A independent variable data array */
//int    *ASpec;      /* vector used to identify user input parm. */
long int which;     /* indicator flag for second load common function */
int flag;     /* indicator flag for ComboMaxLike function */
double con, Acon;  /* storage for constants of log-likelihood */
/** cvl added 8/8/2007 - end block **/
int nparm;         /* # of parameters */
long int restrict;      /* flag for restricting the betas */
/* restrict=1 if the betas >= 0 */
int initial;       /* flag for initialized parameters, initialized=1 */
int appendix;      /* flag for append or overwrite */
/* append = 1, overwrite = 0 */
int smooth;        /* flag for smooth option, 0=unique, 1=C-spline */
int	bmdlCurve;     /* flag for BMDL curve calculation, 1=yes */
int	bmduCurve;     /* flag for BMDU curve calculation, 1=yes */
double xmax=0, xmin, cxmax; /* max and min dose levels */
double scale=0;      /* used to scale dose */
int *boundary;     /* same as bind array, =1 if a boundary was hit */

/* This structure holds details that are not already stored elsewhere */
/* for each tumor. They are used mostly in outputting summary info. */
typedef struct {
  int initial;			/* initialized parms */
  int Nmiss;			/* # incomplete data rows */
  int bmdose;			/* BMD dose calc flag */
  int bmdlCurve;
  int itmax;
  double Rel_Conv;
  double Parm_Conv;
  char DoseHdr[64];
  char PosRespHdr[64];
  char NegRespHdr[64];
  char acDataName[64];
} OtherInfo_t;
OtherInfo_t *pzOtherInfo;

/** changing variable **/
int brat;      /* flag used in BMD calculation */
/* brat = no if background is specified and small */
double BMR;    /* benchmark response */
double BMD_lk; /* log-likelihood  */
double LR;     /* likelihood ratio used in call to Multistage_BMD */
double ck;     /* ck = log(1-A), where A is the added or extra risk */
double back;   /* variable for background parameter */
double back1;  /* variable for background parameter */
int bmdl_bmr_flag; /* flag to see if bmdl converged */
int bmdu_bmr_flag; /* flag to see if bmdu converged */
int PlotLinearExtrapolation;

/* By GeoffreyLNonato-01/14/2009 */
int nT; //Number of tumors
int nt;  /* increment for for loop */

typedef struct ObsData {
	int nObs;
	double *pdYp;
	double *pdYn;
	double *pdXi;
} ObsDataList;
ObsDataList *aDataList;

typedef struct parmVar {
	int nParms;
	long int lnRestrict;
	int nRisk;
	double dEffect;
	double dBMD;
	double dBMDL;
	double dBMDU;
	double dCon;
	double dLogLikelihood;
	double *pdParms;
	int	   *pnSpec;
	int	   *pnIniSp;
	double *pdIniP;
} ParmListVar;
ParmListVar *aParmList;

void GetNewParms2(double *p, int nMaxSize);
void GetMoreParms2(double *p, int nMaxSize);
double ComboMaxLike2(int , double, double *, double **);

extern void getclmt_(long int *lnProbType, long int *nMaxParms, double *dbmr, double *dbmd,
					 double *dtarget, double *pdParms, int *pnFixed,
					 double *pdFixedVal, long int *lnRiskType, double *bmdl,
					 double *pdParms2, long int *optite, long int *nresm,
					 long int *pplnBind, int *);

extern void multiloadcommbloc_(long int *lnT, long int *lnObs, long int *lnMaxNParms, 
		   long int anObs[], long int anPoly[], long int anRest[], long int anParms[]);

extern void parmsloadcommbloc_(long int *, double xParms[]);

extern void getparmscommbloc_(long int *, double xParms[]);

/* End of By GeoffreyLNonato-01/14/2009 */


/****************************************************************
** main--main function used to call Multistage model fitting program.

*****************************************************************/
int main (int argc, char *argv[])
{
	int      iter,i, j;        /* iteration variable */
	int      junk;             /* argument passed in call to Multistage_BMD */
	int      bmdose;           /* flag for computing benchmark dose, 1=yes */
	int      Nmiss;            /* number of records with missing values */
	long int ndegree;          /* degree of polynomial */
	int      nparm_known;      /* number of specified parameters */
	double   lkf;              /* log likelihood for full model */
	double   lkr;              /* log likelihood for reduced model */
	double   xlk;              /* log likelihood for fitted model */
	double   W;                /* proportion of response */
	double   BMD;              /* benchmark dose */
	double   Rel_Conv;         /* relative function convergence */
	double   Parm_Conv;        /* parameter convergence */
	double   ymin, x;
	VarList  *varsum;          /* info for variables--p. dep.,n. dep., indep. */
	AnaList  *anasum;          /* information for ANOVA analysis */
	double   **vcv;            /* variance and covariance matrix */
	char     model_name[82], user_note[82];
	char     dose_name[64], posi_name[64], nega_name[64], junkname[82];
	char     acNameBuffer[512]; /* Hold temp path, possibly long string */
	int      *bounded;         /* bounded[i] = 1 if ith parm hits a bound  */
	double   **vcv_adj;        /* adjusted vcv matrix */
	int      adj_vcv_rows;     /* number of rows in adj. vcv matrix */
	double   *parameters;      /* small parameters are set to 0 */
	char long_path_name[300];
	/********************************************************************
	* CVL 8/2007 - Start block
	* additonal parameters needed for new method of combining tumors
	*********************************************************************/
	int    inm;      /* increment for for loop */
	//int    Anparm;  /* Model A number of parms */
	//double *AParms;     /* model A parameter array */
	double CLL;          /* combined log likelihood of Models A and B */
	/*  Added by CVL 8/2007 - End Block */


	time_t   ltime;

	/* Set time zone from TZ environment variable. If TZ is not set,
	* the operating system is queried to obtain the default value
	* for the variable.
	*/
	/*_tzset();*/
	time( &ltime );


	/********************************************************************
	* {QH 2004/01/14 PR# }
	* Added to show version number if executed from command line with -v
	*********************************************************************/
	if(argc == 2)
		show_version(argv[1], Version_no);

	if(argc < 2)
	{
		fprintf(stderr, "ERROR:  Requires two arguments\nUsage:  %s <file.(d)>\n", argv[0]);
		fprintf (stderr, "   or:  %s -v for version number.\n", argv[0]);
		exit (1);
	} /* end if */


	if (argc > 2)
	{
		path_name2(argc, argv, long_path_name);
		argv[1] = long_path_name;
	}

	/* open the input file */
	fp_in=fopen(argv[1], "r");
	if (fp_in==NULL)
	{
		fprintf(stderr,"Error in opening input  file.\n");
		fprintf (stderr,"...now exiting to system...\n");

		/*	    fprintf(fp_out,"Error in opening input file.\n");*/
		/*	   fprintf (fp_out,"...Exited to system!\n");*/
		exit (1);
	} /* end if */

	/* read in data from input file */
	fscanf(fp_in, "%s", model_name);
	fscanf(fp_in, "%[ ^\n]", user_note);
	fscanf(fp_in, "%[^\n]", user_note);
	fscanf(fp_in, "%s", junkname);
	/* GLN-01/14/2009 */
	//fscanf(fp_in, "%s", junkname);
	fscanf(fp_in, "%d", &nT);
	/* End GLN-01/14/2009 */
	/********************************************************************
	* {CVL - 7/2007
	* Set up a for loop to do through the model fitting etc for two sets
	*********************************************************************/

	/*  Moved by CVL 7/2007 - Start Block 
	*  Code below moved outside of loop    */

	/* get filenames */
	Get_Names(argv[1], fout, fout2, plotfilename);
	/*  Added by CVL 7/2007 - Start Block - get log name */
	for (inm=0;inm<128 && fout[inm] != '.'; inm++)
		fout3[inm] = fout[inm]; /* find period */
	fout3[inm] = '.';
	fout3[inm+1] = 'l';
	fout3[inm+2] = 'o';
	fout3[inm+3] = 'g';
	/*  Added by CVL 7/2007 - End Block */
	if(appendix==Yes)   /* if append was selected, then append output */
		fp_out=fopen(fout,"a");
	else                /* otherwise, overwrite the output */
		fp_out=fopen(fout,"w");

#ifndef RBMDS
	/* open the output file used for the plot */
	fp_out2=fopen(fout2,"w");
	/*  Added by CVL 7/2007 - Start Block */
	/* open the log output file used for tracking progress */
	fp_log=fopen(fout3,"w");
	/*  Added by CVL 7/2007 - End Block */
#endif
	if (fp_out==NULL
#ifndef RBMDS
		||fp_out2==NULL
#endif
		)
	{
#ifdef MISC_OUT
		printf("Error in opening  output files.\n");
		printf ("...now exiting to system...\n");
#endif
		fprintf(fp_out,"Error in opening output files.\n");
		fprintf (fp_out,"...Exited to system!\n");
		exit (1);
	} /* end if */

	/*  Added by CVL 7/2007 - Start Block*/
	if (fp_log==NULL
#ifndef RBMDS
		||fp_log==NULL
#endif
		)
	{
#ifdef MISC_OUT
		printf("Error in opening log output file.\n");
		printf ("...now exiting to system...\n");
#endif
		fprintf(fp_out,"Error in opening log output file.\n");
		fprintf (fp_out,"...Exited to system!\n");
		exit (1);
	} /* end if */
	CLL = 0;
	cxmax = 0;

	/* Allocate memory for meta tumor information.
	 * NOTE: We reserve space for an extra element because the code
	 * expects the first element at index 1 instead of 0.
	 */
	aDataList = (ObsDataList *) malloc((size_t) ((nT+1)*sizeof(ObsDataList)));
	aParmList = (ParmListVar *) malloc((size_t) ((nT+1)*sizeof(ParmListVar)));
	pzOtherInfo = (OtherInfo_t *) malloc((size_t)((nT+1)*sizeof(OtherInfo_t)));
	int nParmMax = 0;

	/* Allocate reusable data structures. These had been incorrectly
	 * allocated within the tumor loop. */
	varsum = VLVECTOR(1, 3);
	anasum = ALVECTOR(1, 3);
	Rlevel = DVECTOR(1, 5);
	Bmdl = DVECTOR(1, 5);
	Bmdu = DVECTOR(1, 5);

	/* The tumor dose groups were previously read in - and the
	 * corresponding single tumor analysis performed - as part
	 * of the same loop. This causes a problem with the combined
	 * BMDL when the maximum dose varies across the groups and
	 * the group with the highest max dose is not first. In this
	 * case, one or more dose groups will be scaled by a different
	 * value.
	 */
	/* This loop reads the data and finds the max dose for all tumors. */
	for (nt = 1; nt <= nT; nt++) {
	  fscanf(fp_in, "%ld%ld%s\n",&Nobs,&ndegree,acNameBuffer);
	  strncpy(pzOtherInfo[nt].acDataName, acNameBuffer, 64);

	  /* GLN-01/14/2009 {Allocate memory and initialize} */
	  aDataList[nt].nObs = Nobs;
	  aDataList[nt].pdYp = DVECTOR(1, Nobs);
	  aDataList[nt].pdYn = DVECTOR(1, Nobs);
	  aDataList[nt].pdXi = DVECTOR(1, Nobs);
	  for(i = 1; i <= Nobs; i++)
	    {
	      aDataList[nt].pdYp[i] = 0;
	      aDataList[nt].pdYn[i] = 0;
	      aDataList[nt].pdXi[i] = 0.0;
	    }

	  /* assign number of parameters */
	  nparm = 1+ndegree;

	  /* GLN-01/14/2009 {Allocate memory and initialized} */
	  if(nparm > nParmMax) nParmMax = nparm;

	  aParmList[nt].nParms = nparm;
	  aParmList[nt].lnRestrict = 0L;
	  aParmList[nt].nRisk = 0;
	  aParmList[nt].dEffect = 0.0;
	  aParmList[nt].dBMD = 0.0;
	  aParmList[nt].dBMDL = 0.0;
	  aParmList[nt].dBMDU = 0.0;
	  aParmList[nt].dCon = 0.0;
	  aParmList[nt].dLogLikelihood = 0.0;
	  aParmList[nt].pdParms = DVECTOR(1, nparm);
	  aParmList[nt].pdIniP = DVECTOR(1, nparm);
	  aParmList[nt].pnSpec = IVECTOR(1, nparm);
	  aParmList[nt].pnIniSp = IVECTOR(1, nparm);
	  for(i = 1; i <= nparm; i++)
	    {
	      aParmList[nt].pdParms[i] = 0.0;
	      aParmList[nt].pdIniP[i] = 0.0;
	      aParmList[nt].pnSpec[i] = 0;
	      aParmList[nt].pnIniSp[i] = 0;
	    }

	  /* read in user input info */
	  fscanf(fp_in,"%d%lf%lf%d%ld%d%d%d", &ITMAX, &Rel_Conv, &Parm_Conv, &bmdlCurve, &restrict, &bmdose, &appendix, &smooth);
	  fscanf(fp_in,"%lf%d%lf",&bmdparm.effect,&bmdparm.risk,&bmdparm.level);
	  /* Save these options for use later */
	  pzOtherInfo[nt].bmdose = bmdose;
	  pzOtherInfo[nt].bmdlCurve = bmdlCurve;
	  pzOtherInfo[nt].itmax = ITMAX;
	  pzOtherInfo[nt].Rel_Conv = Rel_Conv;
	  pzOtherInfo[nt].Parm_Conv = Parm_Conv;

	  /* For now, set these options here in code.  Control
	     whether to plot the linear extrapolation */
	  PlotLinearExtrapolation = 0;

	  junk = 0;

	  /* GLN-01/23/2009 {Assign values} */
	  aParmList[nt].nRisk = bmdparm.risk;
	  aParmList[nt].lnRestrict = restrict;
	  aParmList[nt].dEffect = bmdparm.effect;

	  if (bmdose < 0 || bmdose > 1)
	    ERRORPRT("Error in choosing benchmark dose computation.");

	  /* obtain user input parameters, -9999 is the default value */
	  READ_PARAMETERS(nparm,aParmList[nt].pdParms);

	  /* determine if parameters are to be initialized */
	  fscanf(fp_in,"%d", &initial);
	  pzOtherInfo[nt].initial = initial;

	  /* obtain the initial parameter values */
	  READ_PARAMETERS(nparm,aParmList[nt].pdIniP);
	  /* define vector of flags aParmList[nt].pnIniSp[], 0=initialized */
	  FILL_SPECVECTOR(nparm,aParmList[nt].pdIniP,aParmList[nt].pnIniSp);
	  /* define the specified initial parameter values to be 1 */
	  for(i = 1; i <= nparm; i++)
	    {
	      if(aParmList[nt].pnSpec[i] == 1)
		aParmList[nt].pdIniP[i] = 1.0;
	    } /* end for */

	  /* Read dose group header names. NOTE: We should check */
	  /* to ensure that the string arrays are not overrun!!!!  */
	  fscanf(fp_in,"%s%s%s", dose_name, posi_name, nega_name);
	  /* determine the number of records with missing values */
	  strncpy(pzOtherInfo[nt].DoseHdr, dose_name, 64);
	  strncpy(pzOtherInfo[nt].PosRespHdr, posi_name, 64);
	  strncpy(pzOtherInfo[nt].NegRespHdr, nega_name, 64);

	  Nmiss = READ_OBSDATA3V(Nobs,3,2,3,1, aDataList[nt].pdYp,aDataList[nt].pdYn, aDataList[nt].pdXi);
	  pzOtherInfo[nt].Nmiss = Nmiss;

	  /* aDataList[nt].pdXi[]=dose levels, aDataList[nt].pdYp[]= # of positive responses at the ith dose level */
	  /* aDataList[nt].pdYn[]= total minus aDataList[nt].pdYp[] */

	  Nobs -= Nmiss;       /* extern variable Nobs has been changed */
	  /* Nobs is now the # of obs w/o missing values */

	  /* GLN-01/14/2009 {Allocate memory and initialized} */
	  aDataList[nt].nObs = Nobs;

	  /* print warning if there are too few parameters */
	  if (Nobs < nparm)
	    Warning("Observation # < parameter # for Multistage model.");

	  /*********** end of input data *************/

	  /* Loop across all tumors to ensure correct maximum dose,
	   * which will used to scale all doses: 0 <= Dose <= 1
	   */

	  ymin = aDataList[nt].pdXi[1]; /* not used ?????? */
	  xmin = ymin;
	  /* xmax and scale are global variables that contain the max
	   * dose across all tumors. Do not reinitialize them!!
	   */
	  for (i=1;i<=Nobs;i++) /* obtain min and max dose levels */
	    {
	      x=aDataList[nt].pdXi[i];
	      if (x < xmin) xmin = x;
	      if (x > xmax) xmax = x;
	    } /* end for */
	  scale = xmax;
	} /* end for loop to read data and get xmax */

	/* This loop displays output and analyzes each tumor group */
	for (nt = 1; nt <= nT; nt++) {
	  /* First, populate local variables for the current tumor */
	  Nobs = aDataList[nt].nObs;
	  nparm = aParmList[nt].nParms;
	  ndegree = nparm - 1;
	  restrict = aParmList[nt].lnRestrict;
	  initial = pzOtherInfo[nt].initial;
	  Nmiss = pzOtherInfo[nt].Nmiss;
	  bmdose = pzOtherInfo[nt].bmdose;
	  bmdlCurve = pzOtherInfo[nt].bmdlCurve;
	  ITMAX = pzOtherInfo[nt].itmax;
	  Rel_Conv = pzOtherInfo[nt].Rel_Conv;
	  Parm_Conv = pzOtherInfo[nt].Parm_Conv;
	  strncpy(acNameBuffer, pzOtherInfo[nt].acDataName, 64);
	  strncpy(dose_name, pzOtherInfo[nt].DoseHdr, 64);
	  strncpy(posi_name, pzOtherInfo[nt].PosRespHdr, 64);
	  strncpy(nega_name, pzOtherInfo[nt].NegRespHdr, 64);

	  if (aParmList[nt].pnSpec[1]==1 && aParmList[nt].pdParms[1]<EPS) {
	    /* if background is spec. and small */
	    brat=No;
	  } else {
	    brat=Yes;
	  }

	  vcv = DMATRIX (1,nparm,1,nparm);

	  /* Set bmduCurve equal to bmdlCurve for now. MJF 07AUG05. */
	  bmduCurve = bmdlCurve;

	  /* Define vector of flags for specified parms: */
	  /* 0=unknown, 1=specified. Then count # specified parms. */
	  FILL_SPECVECTOR(nparm,aParmList[nt].pdParms,aParmList[nt].pnSpec);
	  nparm_known = COUNT_SPECVECTOR(nparm,aParmList[nt].pnSpec);

	  /* Print model and file information on output page */
	  Output_Header(Version_no, argv[1], plotfilename, ctime(&ltime), user_note);
	  fflush(fp_out);
	  OUTPUT_TEXT("\n   The form of the probability function is: ");
	  OUTPUT_TEXT("\n   P[response] = background + (1-background)*[1-EXP(");
	  fprintf (fp_out,"                 -beta1*dose^1");
	  for (i=2;i<=ndegree;i++)
	    fprintf (fp_out,"-beta%d*dose^%d",i,i);
	  fprintf (fp_out,")]");

	  if (restrict==Yes) /* if the option to restrict the betas is chosen */
	    OUTPUT_TEXT("\n\n   The parameter betas are restricted to be positive");
	  else
	    OUTPUT_TEXT("\n\n   The parameter betas are not restricted");
	  fprintf(fp_out,"\n\n   Dependent variable = %s", posi_name);
	  fprintf(fp_out,"\n   Independent variable = %s", dose_name);
	  fprintf(fp_out,"\n   Data file name = %s", acNameBuffer);
	  if (aParmList[nt].pnSpec[1]==Yes)  /* if the background parameter was specified */
	    {
	      fprintf(fp_out,"\n");
	      if (aParmList[nt].pdParms[1] <= 0.0000001) /* if the background is small */
		fprintf(fp_out,"\n   Background parameter is set to zero");
	      else
		fprintf(fp_out,"\n   Background parameter is set to %g", aParmList[nt].pdParms[1]);
	      if ((aParmList[nt].pdParms[1] < 1.0e-20) && (aDataList[nt].pdXi[1] < 1.0e-20) && (aDataList[nt].pdYp[1] > 0))
		ERRORPRT("ERROR: Background parameter specified as 0, but there were responses\n     at the control dose");
	    } /* end if */

	  if(nparm_known > 0+aParmList[nt].pnSpec[1])
	    {
	      fprintf (fp_out, "\n\n   User specifies the following parameters:");
	      for (i=2; i<= nparm; i++)
		{
		  if(aParmList[nt].pnSpec[i] == 1)
		    fprintf (fp_out, "\n %15s = %10.5g", Parm_name[i-1], aParmList[nt].pdParms[i]);
		}  /* end for */
	      fprintf (fp_out, "\n");
	    } /* end if */
	  /* Why is this being done again??? It looks the same as before. */
	  /* determine the # of specified parameters */
	  nparm_known = COUNT_SPECVECTOR(nparm, aParmList[nt].pnSpec);
	  /* print out information */
	  fprintf (fp_out, "\n\n Total number of observations = %ld",Nobs+Nmiss);
	  fprintf (fp_out, "\n Total number of records with missing values = %d",Nmiss);
	  fprintf (fp_out, "\n Total number of parameters in model = %d", nparm);
	  fprintf (fp_out, "\n Total number of specified parameters = %d\n", nparm_known);
	  fprintf(fp_out, " Degree of polynomial = %ld\n\n", ndegree);
	  fprintf(fp_out, "\n Maximum number of iterations = %d\n", ITMAX);
	  fprintf(fp_out, " Relative Function Convergence has been set to: %g\n", Rel_Conv);
	  fprintf(fp_out, " Parameter Convergence has been set to: %g\n\n", Parm_Conv);


	  if(Rel_Conv != 1.0e-8 || Parm_Conv != 1.0e-8)
	    {
	      fprintf(fp_out, "****  We are sorry but Relative Function and Parameter Convergence    ****\n");
	      fprintf(fp_out, "****  are currently unavailable in this model.  Please keep checking  ****\n");
	      fprintf(fp_out, "****  the web sight for model updates which will eventually           ****\n");
	      fprintf(fp_out, "****  incorporate these convergence criterion.  Default values used.  ****\n\n");
	    } /* end if */

	  if(initial==Yes)  /* if initial option is chosen */
	    {
	      OUTPUT_TEXT("\n\n                 User Inputs Initial Parameter Values  ");
	      OUTPUT_Init(nparm, aParmList[nt].pnSpec, aParmList[nt].pdIniP, Parm_name);

	      for (i=1; i<=nparm; i++)
		{
		  if(aParmList[nt].pnIniSp[i]==1)    /* if the parameter was initialized */
		    {
		      if(aParmList[nt].pnSpec[i]==1 )  /* check to see if the parameter is fixed */
			Warning("The initial value for the fixed parameter is ignored.");
		    } /* end if */
		  else
		    {
		      /* check if all the unspecified parms were initialized */
		      if (aParmList[nt].pnSpec[i]==0)
			ERRORPRT("When the initial option is chosen, one has to initial ALL unspecified parameters.");
		    }  /* end else */
		} /* end for */
	      i=0;
	      if(restrict==Yes)  /* if the restrict betas option was chosen */
		{
		  for(j=2;j<=nparm;j++)
		    {
		      if(aParmList[nt].pdIniP[j]<0) i++;  /* if initial value is negative */
		    } /* end for */
		} /* end if */
	      if (aParmList[nt].pdIniP[1] < 0 || aParmList[nt].pdIniP[1]>1 || i>0)
		ERRORPRT("The initial values have to be: 1 > Bg >= 0 and Betas >= 0 ). ");

	    }  /* end if (initial==Yes) */

	  /* compute init_lkf for full model and init_lkr for reduced model */
	  lkf = 0.0;
	  varsum[1].S = 0;
	  varsum[2].S = 0;
	  for (i=1;i<=Nobs;i++)
	    {
	      varsum[1].S += aDataList[nt].pdYp[i];
	      varsum[2].S += aDataList[nt].pdYn[i];
	      W = aDataList[nt].pdYp[i] / (aDataList[nt].pdYp[i]+aDataList[nt].pdYn[i]); /* proportion of response */
	      if (W > 0)   lkf += aDataList[nt].pdYp[i] * log(W);
	      if (W < 1)   lkf += aDataList[nt].pdYn[i] * log(1- W);
	    } /* end for */
	  W = varsum[1].S / (varsum[1].S + varsum[2].S);
	  lkr = varsum[1].S * log(W) + varsum[2].S * log(1- W);

	  /* fitting Multistage model and output parameter estimators */
	  fflush(fp_out);

	  if(nparm_known<nparm)
	    Multistage_fit(nparm, aParmList[nt].pdParms, EPS, &iter, &xlk);

	  /* aParmList[nt].pdParms[] is now the fitted parameters */
	  /* xlk is the log-likelihood            */
	  fflush(fp_out);
	  /* initialize vcv (varaince-covariance matrix) */
	  INITIALIZE_DMATRIX(vcv, nparm, nparm);
	  /* compute the approximate variance-covariance matrix */
	  /* variance-covariance matrix is the inverse of the info. matrix */
	  Multistage_vcv(nparm,aParmList[nt].pnSpec,aParmList[nt].pdParms,vcv);

	  /* define bounded vector */
	  bounded = IVECTOR(1, nparm);
	  for (i=1; i<=nparm; i++)
	    {
	      if(nparm_known<nparm)
		bounded[i] = boundary[i];
	      else
		bounded[i]=0;
	    }
	  if (aParmList[nt].pdParms[1] == 0)
	    bounded[1] = 1;

	  adj_vcv_rows = 0;
	  for (i=1; i<=nparm; i++)
	    {
	      if (bounded[i] == 0)
		{
		  adj_vcv_rows ++;
		} /* end if */
	    } /* end for */

	  vcv_adj = DMATRIX(1, adj_vcv_rows, 1, adj_vcv_rows);

	  /* output correlation matrix */
	  fflush(fp_out);
	  if (nparm-nparm_known>0)
	    Get_and_OUTPUT_DTMSVCV(nparm,aParmList[nt].pnSpec,Parm_name,vcv,vcv_adj,bounded);

	  fflush(fp_out);
	  parameters = DVECTOR(1, nparm);
	  for (i=1; i<=nparm; i++)
	    {
	      if (bounded[i] == 1)
		parameters[i] = 0;
	      else
		parameters[i] = aParmList[nt].pdParms[i];
	    }

	  /* output model fitting parameters and standard errors */
	  OP_ParmsE(nparm,aParmList[nt].pnSpec,parameters,Parm_name,vcv_adj, bounded, bmdparm.level, 0);
	  fflush(fp_out);
	  /* free memory */
	  /*  FREE_IVECTOR(bounded, 1, nparm);   need bounded for DTMS3ANOVA */
	  FREE_DVECTOR(parameters, 1, nparm);

	  /* compute ANOVA table elements */
	  DTMS3ANOVA (nparm,Nobs,aParmList[nt].pnSpec,lkf,xlk,lkr,anasum, bounded);
	  fflush(fp_out);

	  anasum[2].TEST = CHISQ(anasum[2].MSE, anasum[2].DF);
	  fflush(fp_out);
	  /* output ANOVA table */
	  OUTPUT_DTMS3ANOVA(anatxt,anasum);
	  fflush(fp_out);

	  con = LogLik_Constant(Nobs, aDataList[nt].pdYp, aDataList[nt].pdYn);
	  aParmList[nt].dCon = con;
	  fflush(fp_out);

	  /* print a goodness of fit table */
	  Quantal_Goodness(nparm, bounded, aParmList[nt].pdParms, Nobs, aDataList[nt].pdXi, aDataList[nt].pdYp, aDataList[nt].pdYn, scale);
	  fflush(fp_out);

	  /******************* compute benchmark dose ***********************/

#ifndef RBMDS
	  /* print info to file used for plot */
	  fprintf (fp_out2, "\n BMD_flag \t %d \n Nobs \t%ld \n nparm \t%d",  bmdose, Nobs, nparm );
	  fprintf (fp_out2, "\n  Con_lev \t%3.3g ", bmdparm.level);
	  fprintf (fp_out2, "\n  RiskType \t%d ", bmdparm.risk);
	  fprintf (fp_out2, "\n  Effect \t%3.3g ", bmdparm.effect);
	  for (i=1;i<=nparm; i++) fprintf (fp_out2, "\n %s \t %5.3g", Parm_name[i-1], aParmList[nt].pdParms[i]);

		/* Calculate 95% CIs at each dose level for graphical output */
		{
			double *LL, *UL, *estp;

			LL = DVECTOR(1, Nobs);
			UL = DVECTOR(1, Nobs);
			estp = DVECTOR(1, Nobs);
			Quantal_CI(Nobs, aDataList[nt].pdYp, aDataList[nt].pdYn, 0.95, LL, estp, UL);
			fprintf (fp_out2,"\n\n Data");
			for (i=1;i<=Nobs;i++)
			{
				fprintf (fp_out2,"\n %f %f %f %f", aDataList[nt].pdXi[i], estp[i], LL[i], UL[i]);
			} /* end for */
			FREE_DVECTOR(LL, 1, Nobs);
			FREE_DVECTOR(UL, 1, Nobs);
			FREE_DVECTOR(estp, 1, Nobs);
		}

		fprintf (fp_out2,"\n Max_Min_dose \n  %f %f ", xmax, xmin);
#endif
		if (bmdose==Yes) /* if BMD calculation was selected */
		{
			back=0.0;
			/* brat = no if aParmList[nt].pnSpec[1] = 1 and Parms[1] < EPS */
			if(brat==1) back=aParmList[nt].pdParms[1];
			back1=1-back;
			if (bmdparm.risk==1) /* if risk type is added */
				back1=1;

			bmdl_bmr_flag = 0;
			bmdu_bmr_flag = 0;
			/* obtain the BMD and BMDL's */
			Multistage_BMD (nparm, aParmList[nt].pdParms, EPS, &junk, xlk, Rlevel, Bmdl, Bmdu, &BMD);

		} /* end: if (bmdose==Yes) */

		/*  Added by CVL 8/2007 - Start Block */
		CLL = CLL + xlk;
		if (cxmax < xmax) cxmax = xmax;
		/* the Model B fit  */
	} /* end of loop for doing two sets of data */
	/* Start of the additional code for the combined log-likelihood and 
	the combined bmd and bmdl - CVL 8/2007*/

	scale = 1.0;

	fprintf(fp_out, "\n\n**** Start of combined BMD and BMDL Calculations.****\n\n");
	fprintf(fp_log, "\n\n**** Start of combined BMD and BMDL Calculations.****\n\n");

	con = 0.0;
	for (nt = 1; nt <= nT; nt++)
	{
		con = con + aParmList[nt].dCon;
	}

	fprintf(fp_out, "  Combined Log-Likelihood          %30.22g \n",CLL);
	fprintf(fp_log, "  Combined Log-Likelihood          %30.22g \n",CLL);
	fprintf(fp_out,"\n  Combined Log-likelihood Constant %30.22g \n",con);
	fprintf(fp_log,"\n  Combined Log-likelihood Constant %30.22g \n",con);

	/* GLN-02/03/2009 Allocate memories to for loadblock */
	long int lnT = (long int)nT;
	long int lnObs, lnMaxNParms;
	long int *anObs, *anPoly, *anRest, *anParms;

	anObs = (long int *) malloc((size_t)10*sizeof(long int));
	anPoly = (long int *) malloc((size_t)10*sizeof(long int));
	anRest = (long int *) malloc((size_t)10*sizeof(long int));
	anParms = (long int *) malloc((size_t)10*sizeof(long int));

	lnObs = lnMaxNParms = 0;
	for(nt = 1; nt<= nT; nt++)
	{
		if((long int)aDataList[nt].nObs > lnObs)
			lnObs = (long int)aDataList[nt].nObs;

		if((long int)aParmList[nt].nParms > lnMaxNParms)
			lnMaxNParms = (long int)aParmList[nt].nParms;
	}

	for(nt = 0; nt < 10; nt++)
	{
		anObs[nt] = 0L;
		anPoly[nt] = 0L;
		anRest[nt] = 0L;
		anParms[nt] = 0L;
	}

	for(nt = 1; nt <= nT; nt++)
	{
		anObs[nt-1] = (long int)aDataList[nt].nObs;
		anPoly[nt-1] = (long int)aParmList[nt].nParms -1;
		anRest[nt-1] = aParmList[nt].lnRestrict;
		anParms[nt-1] = (long int)aParmList[nt].nParms;
	}
	fprintf(fp_log,"\n\nObs         Poly         nRest         nParms\n");
	for(nt = 1; nt <= nT; nt++)
	{
		fprintf(fp_log,"%ld          %ld          %ld         %ld\n", anObs[nt-1], anPoly[nt-1], anRest[nt-1], anParms[nt-1]);
	}

	multiloadcommbloc_(&lnT, &lnObs, &lnMaxNParms, 
		anObs, anPoly, anRest, anParms);

	Multistage_ComboBMD (1e-12, &junk, CLL, Rlevel, Bmdl, Bmdu, &BMD);

	/* End GLN-02/03/2009 Allocate memories to for loadblock */

	/* GLN-01/14/2009 {Release allocated memories} */
	int nO;
	for(nt = 1; nt <= nT; nt++)
	{
		Nobs = aDataList[nt].nObs;
		nparm = aParmList[nt].nParms;
		fprintf(fp_log,"\n\nTumor %d\n  Observation = %d\n", nt, aDataList[nt].nObs);
		fprintf(fp_log,"\n  pdYp -> {");
		for(nO = 1; nO <= aDataList[nt].nObs; nO++)
		{
			fprintf(fp_log,"%g, ", aDataList[nt].pdYp[nO]);
		}
		fprintf(fp_log,"}\n");

		fprintf(fp_log,"\n  pdYn -> {");
		for(nO = 1; nO <= aDataList[nt].nObs; nO++)
		{
			fprintf(fp_log,"%g, ", aDataList[nt].pdYn[nO]);
		}
		fprintf(fp_log,"}\n");

		fprintf(fp_log,"\n  pdXi -> {");
		for(nO = 1; nO <= aDataList[nt].nObs; nO++)
		{
			fprintf(fp_log,"%g, ", aDataList[nt].pdXi[nO]);
		}
		fprintf(fp_log,"}\n");

		FREE_DVECTOR (aDataList[nt].pdYp, 1, Nobs);
		FREE_DVECTOR (aDataList[nt].pdYn, 1, Nobs);
		FREE_DVECTOR (aDataList[nt].pdXi, 1, Nobs);

		FREE_DVECTOR (aParmList[nt].pdParms, 1, nparm);
		FREE_DVECTOR (aParmList[nt].pdIniP, 1, nparm);
		FREE_IVECTOR (aParmList[nt].pnSpec, 1, nparm);
		FREE_IVECTOR (aParmList[nt].pnIniSp, 1, nparm);
	}
	free ((FREE_ARG) (aDataList));
	free ((FREE_ARG) (aParmList));
	free (pzOtherInfo);

	//FREE_DVECTOR(mdDoses, 0, nT-1, 0, (int)(lnObs-1));
	//FREE_IMATRIX(mnAnim, 0, nT-1, 0, (int)(lnObs-1));
	//FREE_IMATRIX(mnAffect, 0, nT-1, 0, (int)(lnObs-1));
	free(anObs);
	free(anPoly);
	free(anRest);
	free(anParms);
	//free(adDoses);
	//free(anAnim);
	//free(anAffect);

	//FREE_DVECTOR(adXmax, 0, nT-1);

	/* End of GLN-01/14/2009 {Release allocated memories} */

	FREE_DVECTOR (Ep,1,Nobs);
	FREE_IVECTOR (bounded, 1, nparm);
	FREE_VLVECTOR(varsum,1,3);
	FREE_ALVECTOR(anasum,1,3);
	FREE_DVECTOR (Rlevel,1,5);
	FREE_DVECTOR (Bmdl, 1, 5);
	FREE_DVECTOR (Bmdu, 1, 5);
	FREE_DMATRIX(vcv,1,nparm,1,nparm);
	FREE_IVECTOR(boundary, 1, nparm);
	FREE_DMATRIX(vcv_adj,1,adj_vcv_rows,1,adj_vcv_rows);
	CLOSE_FILES ();
	/*  Added by CVL 7/2007 - Start Block */
#ifndef RBMDS
	/* close the log output file used for tracking progress */
	fclose(fp_log); 
#endif
	/*  Added by CVL 7/2007 - End Block */

	return(0);     /* This returns an int = 0 until completion of debugging */

} /* end of main */


/***********************************************************
* Predict -- returns predicted values for given parameter
*            values.
*
*            input:
*                   doses is the vector of doses (1 - based)
*                   ndoses is the number of doses (1 - based)
*                   Parms is the vector of parameters (1 - based)
*            output:
*                   P is the vector (of length ndoses) of predicted
*                     proportions
***********************************************************/

void Predict(double doses[], int ndoses, double Parms[], double P[])
{
	int i, j;
	double x, ex1;
	for (i=1;i<=Nobs;i++)
	{
		x = aDataList[nt].pdXi[i];   /* ith dose level */
		/* calculate the probability of response */
		ex1=Parms[nparm];
		for (j=nparm-1; j>=2; j--) ex1=ex1*x+Parms[j];
		ex1=ex1*x;
		P[i] = Parms[1]+(1-Parms[1])*(1-exp(-ex1)); /* prob. of response */
	} /* end for */
}


/******************************************************************
**Multistage_vcv -- used to compute the info matrix for Multistage model.
*  External:  Nobs, Xi[], Yp[], Yn[]
*  input:
*   nparm is the number of parameters
*   Spec[] is a vector of flags for specified parameters
*   p[] is the vector of fitted parameters
*   vcv is the variance-covariance matrix
*  output: vcv
******************************************************************/
void Multistage_vcv(int nparm, int Spec[], double p[], double **vcv)
{
	double  x, ex, ploy, ex2, lk1, lk2;
	double  *d1;
	int i,j,k;

	d1=DVECTOR(1,nparm); /* vector for first derivatives */

	for (i=1;i<=Nobs;i++)
	{ /* compute the probability of response at the ith dose */
		x = aDataList[nt].pdXi[i];
		ploy=p[nparm];
		for (j=nparm-1; j>=2; j--) ploy=ploy*x+p[j];
		ploy=ploy*x;
		ex = p[1]+(1-p[1])*(1-exp(-ploy)); /* probability of response */
		PROBABILITY_INRANGE(&ex); /* make sure prob is between 0 and 1 */
		ex2=exp(-ploy);
		lk1 = (aDataList[nt].pdYp[i]-ex*aDataList[nt].pdYn[i]/(1-ex));
		lk2 = -aDataList[nt].pdYp[i]/ex - ex*aDataList[nt].pdYn[i]/(1-ex)/(1-ex);

		d1[1]= ex2;
		for (j=2;j<=nparm;j++)
			d1[j]= (1-p[1])*ex2*pow(x, j-1); /* first der. of ex wrt beta_j */


		vcv[1][1] -= lk2*d1[1]*d1[1];
		for (j=2;j<=nparm;j++)
			vcv[1][j] -= lk2*d1[1]*d1[j] - lk1*ex2*pow(x, j-1);

		for (j=2;j<=nparm;j++)
		{
			if (Spec[j]==0) /* if p[j] is not specified */
			{
				for(k=j;k<=nparm; k++)
					vcv[j][k] -= lk2*d1[j]*d1[k] - lk1*(1-p[1])*ex2*pow(x, k+j-2);

			} /* end if */
		} /* end for */

	} /* end for */
	for (j=1;j<=nparm;j++)
	{
		for(k=j;k<=nparm; k++)
			vcv[k][j] =vcv[j][k]; /* vcv matrix is symmetric */
	} /* end for */
	FREE_DVECTOR(d1, 1,nparm); /* free memory */

	return;

} /* end: Multistage_vcv */

/*********************************************************************
*Multistage_fit -- Used to "prepare" the data for further computation,
*            i.e. compute the extern variables, give the initial
*            parameters, etc. THEN fit the Multistage model.
*            (In fact, these jobs could be done in main().)
* external: Nobs, scale, Xi[], Yp[], Yn[]
* input:
*  nparm is the number of parameters
*  p[] is the vector of parameters
*  gtol is a small positive number
*  iter is the number of iterations
*  fret is a log-likelihood
* output: prints initial parameter values
**********************************************************************/
void Multistage_fit(int nparm, double p[], double gtol,
					int *iter, double *fret)
{
	long int *Spec2, *bind, ndegree, optite, nresm;
	double *affect, *nanim;
	int    *SpBak, i, j, k, nparm2=0, l;
	int     ii, ntNow;
	double  Yt, N, ll;
	double *pBak, *tmy, *t, **tmv, **X, **XPX, *Y, **XP, *stParm, *XPY, **XT;
	double *fitparms, *parms, *doses, lminbmd, lmaxbmd;


	/* allocate memory for arrays */
	pBak=DVECTOR(1, nparm);
	tmy=DVECTOR(1, nparm);
	t=DVECTOR(1, nparm);
	tmv=DMATRIX(1,nparm,1,nparm);
	SpBak=IVECTOR(1, nparm);
	X = DMATRIX(1, Nobs, 1, nparm);
	XP = DMATRIX(1, nparm, 1, Nobs);
	XPX = DMATRIX(1, nparm, 1, nparm);
	Y = DVECTOR(1, Nobs);
	stParm = DVECTOR(1, nparm);
	XPY = DVECTOR(1, nparm);

	ndegree = nparm-1;

	for (j=1; j<=nparm; j++)
		pBak[j]=p[j];             /* save the input p[] */

	/* allocate memory */
	doses = (double *) malloc((size_t)(Nobs)*sizeof(double));
	affect = (double *) malloc((size_t)(Nobs)*sizeof(double));
	nanim = (double *) malloc((size_t)(Nobs)*sizeof(double));
	parms = (double *) malloc((size_t)(nparm)*sizeof(double));
	fitparms = (double *) malloc((size_t)(nparm)*sizeof(double));
	Spec2 = (long int *) malloc((size_t)(nparm)*sizeof(long int));
	bind = (long int *) malloc((size_t)(nparm)*sizeof(long int));

	for(i = 1; i <= Nobs; i++)
	{
		nanim[i-1] = aDataList[nt].pdYp[i] + aDataList[nt].pdYn[i]; /* # of animals at ith dose level */
		doses[i-1] = aDataList[nt].pdXi[i]/xmax;                 /* dose levels; rescaled to: 0 <= Dose <= 1  */
		affect[i-1] = aDataList[nt].pdYp[i];        /* # of affected animals at ith dose level */
	} /* end for */

	/* calculate the values that will correspond to '0' and 'Infinity'
	in BMD calculations */
	lminbmd = log(DBL_MIN) - log(xmax);
	lmaxbmd = log(DBL_MAX) - log(xmax);
	/* load values for use in getmle, getcl, and getprofile */
	/* into the common block */
	ntNow = nt-1;
	loadcommbloc_(&Nobs,&xmax,nanim,doses,affect,&ndegree,&restrict,
		&lminbmd,&lmaxbmd,&ntNow);


	/* scale parameters */
	for(i = 1; i<=nparm; i++)
	{
		p[i] = p[i]*(pow(xmax,(i-1)));
		if(initial == Yes) aParmList[nt].pdIniP[i]= aParmList[nt].pdIniP[i]*(pow(xmax,(i-1)));
	}


	/****** Obtain initial estimations for p[] ******/
	if(initial==Yes)
	{
		for(j=1; j<=nparm; j++)
			p[j]=aParmList[nt].pdIniP[j];
	} /* end if */
	else
	{
		/* Get initial values for the parameters by regressing
		log(1-P) on log(1 - p[1]) - Sum(p[i]*X^i), if the
		parameters are restricted to be greater than 0, then if a
		negative parameter value occurs, that parameter is set
		to zero, and all the other parameters are re-estimated
		by regressing the same model above, minus the negative
		parameter  (Here, P is the sample proportion at each
		dose */

		for(i = 1; i <= Nobs; i++)	/* Setting up appropriate  */
		{							/* X and Y matrices for    */
			X[i][1] = 1.0;			/* regression              */
			if(aDataList[nt].pdYn[i] > 0)
				Y[i] = log( 1 - aDataList[nt].pdYp[i]/(aDataList[nt].pdYp[i]+aDataList[nt].pdYn[i]) );
			else
				Y[i] = -1.0e20;		/* To "quantify" log(0)    */

			for(j = 2; j <= nparm; j++)
			{
				if(doses[i-1] == 0.0)
					X[i][j] = 0.0;
				else
					X[i][j] = -1*pow(doses[i-1], j-1);
			} /* end for */
		} /* end for */

		/* Perform Multiple Linear Regression to get starting point */

		TRANSPOSE(X, XP, Nobs, nparm); /* XP is the transpose of X */

		MATMPYM2(XP, X, XPX, nparm, Nobs, nparm); /* XPX = XP*X */

		INVMAT(XPX, nparm); /* XPX is now the inverse of XPX */

		MATMPYV2(nparm, Nobs, XP, Y, XPY); /* XPY = XP*Y */

		MATMPYV2(nparm, nparm, XPX, XPY, stParm); /* stParm is the vector of reg. parms */

		/* Need to prevent values of p[1] that are too large */
		if(stParm[1] >= 0)
			p[1] = 0.0;
		else
			p[1] = 1.0 - exp(stParm[1]);

		for(i = 2; i <= nparm; i++)
			p[i] = stParm[i]; /* p[] are initial parameter estimates */

		nparm2 = nparm;

		if (restrict==Yes ) /* if the betas are being restricted */
		{

			k = 0;
			XT = DMATRIX(1, Nobs, 1, nparm);

			do{

				for (j=nparm;j>=2;j--)
				{
					if(p[j] < 0) /* if there is a negative estimate */
					{
						k = j;  /* beta(k) will be left out of the next regression */
						nparm2 = nparm2-1;
						break;
					} /* end if */
					else
						k = 0;
				} /* end for */

				if(nparm2 == 1) /* if there is only 1 parm in the model */
				{
					N = Yt = 0;
					for(i = 1; i <= Nobs; i++)
					{
						Yt += aDataList[nt].pdYp[i];
						N += (aDataList[nt].pdYp[i]+aDataList[nt].pdYn[i]);
					} /* end for */
					p[1] = Yt/N; /* define the parm to be a proportion */
					k = 0;
				} /* end if */

				if(k != 0)
				{
					for(i = 1; i <= Nobs; i++)
					{
						l = 0;
						for(j = 1; j <= nparm2+1; j++)
						{
							if(j != k)
							{
								XT[i][j-l] = X[i][j]; /* XT is the design  matrix X */
							}  /* end if */           /* with the kth column removed */
							else
								l = 1;
						} /* end for */
					} /* end for */
					/* define matrices for another regression with one less parameter */
					TRANSPOSE(XT, XP, Nobs, nparm2); /* XP is the transpose of XT */

					MATMPYM2(XP, XT, XPX, nparm2, Nobs, nparm2); /* XPX = XP*XT */

					INVMAT(XPX, nparm2); /* XPX is now the inverse of XPX */

					MATMPYV2(nparm2, Nobs, XP, Y, XPY); /* XPY = XP*Y */

					MATMPYV2(nparm2, nparm2, XPX, XPY, stParm); /* stParm is the new vector */
					/* of initial estimates     */
					if(stParm[1] >= 0)
						p[1] = 0.0; /* prevent p[1] from being too large */
					else
						p[1] = 1.0 - exp(stParm[1]);

					j = 0;

					for(i = 2; i <= nparm; i++)
					{
						if(i == k || p[i] == 0.0)
						{
							p[i] = 0.0;
							j++;
						} /* end if */
						else
							p[i] = stParm[i-j];
					} /* end for */
					/* free memory */
					FREE_DMATRIX(X, 1, Nobs, 1, nparm2+1);
					/* alocate memory */
					X = DMATRIX(1, Nobs, 1, nparm2);

					for(i = 1; i <= Nobs; i++)
						for(j= 1; j <= nparm2; j++)
							X[i][j] = XT[i][j];

				} /* end if (k != 0) */

			} /* end do */
			while(k != 0);

			FREE_DMATRIX(XT, 1, Nobs, 1, nparm2); /* free up some memory */

		} /* end if (restrict==Yes) */
		if(p[1] < 0)
			p[1] = 0;

		/* unscale for output line */
		for(i = 1; i<=nparm; i++)
		{
			p[i] = p[i]/(pow(xmax,(i-1)));
			if(initial == Yes) aParmList[nt].pdIniP[i]= aParmList[nt].pdIniP[i]/(pow(xmax,(i-1)));
		}

		/* output initial parameter values */
		OUTPUT_TEXT("\n\n                  Default Initial Parameter Values  ");
		OUTPUT_Init(nparm, aParmList[nt].pnSpec, p, Parm_name);

		/* scale back after output */
		for(i = 1; i<=nparm; i++)
		{
			p[i] = p[i]*(pow(xmax,(i-1)));
			if(initial == Yes) aParmList[nt].pdIniP[i]= aParmList[nt].pdIniP[i]*(pow(xmax,(i-1)));
		}
	} /* end else */


	/* free up memory */
	FREE_DMATRIX(X, 1, Nobs, 1, nparm2);
	FREE_DMATRIX(XP, 1, nparm2, 1, Nobs);
	FREE_DMATRIX(XPX, 1, nparm2, 1, nparm2);
	FREE_DVECTOR(Y, 1, Nobs);
	FREE_DVECTOR(stParm, 1, nparm2);
	FREE_DVECTOR(XPY, 1, nparm2);

	/** get specified parameters **/
	for (j=1;j<=nparm;j++)
		if (aParmList[nt].pnSpec[j]==Yes) p[j]=pBak[j];
	/* p[] is now a vector of specified parms and regressed parms */


	for(i = 1; i <= nparm; i++)
	{
		if(initial == Yes && Spec[i] == 0)  /* if the initial option is chosen and */
		{                                   /* if the parameter is not specified   */
			p[i] = aParmList[nt].pdIniP[i];
			Spec2[i] = 0;
		} /* end if */
		else
			Spec2[i-1] = aParmList[nt].pnSpec[i];
		parms[i-1] = p[i]; /* shift index */
	}   /* end for */
	/* parms[] are the starting values to get the MLE's */

	/* need background between 0 and 1 */
	if ((parms[0] < 0) || (parms[0] >= 1))
		parms[0] = 0;

	parms[0] = -log(1.0 - parms[0]);

	/* prevent any parameter from being too large */
	for (i=0; i<=nparm-1; i++)
	{
		if (fabs(parms[i]) > 999999999)
			parms[i] = (double) rand()/ (double) RAND_MAX;
	}

	/* obtain the mle's */
	getmle_(parms, Spec2, parms, fitparms, &ll, &optite,
		&nresm, bind);
	/* note: bind[]=1 if a boundary was hit */
	/* nresm is the # of equality and inequality constraints */

	if((optite < 0 ) || (optite > 2))
	{
#ifdef MISC_OUT
		/* Warn user */
		fprintf(fp_out, "**** WARNING:  Completion code = %ld.  Optimum not found. Trying new starting pont****\n\n", optite);
#endif
		/* Try up to 10 times if needed */
		for(ii = 0; ii < 10; ii++)
		{
			/* again, reparameterize p[0] */
			parms[0] = -log(1-p[1]);

			for(j = 2; j <= nparm; j++)  /* Get original values */
				parms[j-1] = p[j];

			GetNewParms(parms, nparm);  /* Get a new starting point */
			/* parms[] are now new starting values */

			/* Try again */
			getmle_(parms, Spec2, parms, fitparms, &ll, &optite,
				&nresm, bind);

			/* if optite >= 0, it is successful, and we can stop */
			if ((optite >= 0) && (optite <= 2))
				break;
#ifdef MISC_OUT
			/* otherwise, issues another warning, and continue trying */
			else
				fprintf(fp_out, "**** WARNING %d:  Completion code = %ld trying new start****\n\n", ii, optite);
#endif
		} /* end: for (ii=0; ii < 10; ii++) */
	} /* end: if (optite < 0) */




	/***************************************************************/


	/* try perturbing the starting values a little more */
	if((optite < 0) || (optite > 2))
	{
#ifdef MISC_OUT
		/* Warn user */
		printf("**** WARNING:  Completion code = %ld.  Optimum not found. Trying new starting point****\n\n", optite);
#endif
		/* Try up to 10 times if needed */
		for(ii = 0; ii < 10; ii++)
		{
			/* again, reparameterize p[0] */
			parms[0] = -log(1-p[1]);

			for(j = 2; j <= nparm; j++)  /* Get original values */
				parms[j-1] = p[j];

			GetMLEParms(parms, nparm);  /* Get a new starting point */
			/* parms[] are now new starting values */

			/* Try again */
			getmle_(parms, Spec2, parms, fitparms, &ll, &optite,
				&nresm, bind);

			/* if optite >= 0, it is successful, and we can stop */
			if ((optite >= 0) && (optite <= 2))
				break;

			/* otherwise, issues another warning, and continue trying */
#ifdef MISC_OUT
			else
				printf("**** WARNING %d:  Completion code = %ld trying new start****\n\n", ii, optite);
#endif

		} /* end: for (ii=0; ii < 10; ii++) */
	} /* end: if (optite < 0) */

	/*try random numbers as starting values */




	/****************************************************************/


	if((optite < 0) || (optite > 2))
	{
#ifdef MISC_OUT
		/* Warn user */
		fprintf(fp_out, "**** WARNING:  Completion code = %ld.  Optimum not found. Trying new starting point****\n\n", optite);
#endif
		/* Try up to 10 times if needed */
		for(ii = 0; ii < 10; ii++)
		{
			/* again, reparameterize p[0] */
			parms[0] = -log(1-p[1]);

			for(j = 2; j <= nparm; j++)  /* Get original values */
				parms[j-1] = p[j];

			GetMoreParms(parms, nparm);  /* Get a new starting point */
			/* parms[] are now new starting values */

			/* Try again */
			getmle_(parms, Spec2, parms, fitparms, &ll, &optite,
				&nresm, bind);

			/* if optite >= 0, it is successful, and we can stop */
			if ((optite >= 0) && (optite <= 2))
				break;
#ifdef MISC_OUT
			/* otherwise, issues another warning, and continue trying */
			else
				fprintf(fp_out, "**** WARNING %d:  Completion code = %ld trying new start****\n\n", ii, optite);
#endif
		} /* end: for (ii=0; ii < 10; ii++) */
	} /* end: if (optite < 0) */

	if ((optite < 0) || (optite > 2))
		ERRORPRT("Fitted parameters did not converge");

	p[1] = 1.0 - exp(-fitparms[0]);

	for(i = 2; i <= nparm; i++)
		p[i] = fitparms[i-1];
	/* p[] are now the fitted values */

	*fret = ll; /* log-likelihood */

	boundary = IVECTOR(1, nparm);
	for (i=1; i<=nparm; i++)
		boundary[i] = bind[i-1];

	/* remove scaling from doses */
	scale=1;
	for (i=1; i<=Nobs; i++)
	{
		/*     Xi[i]=Xi[i]*xmax;  */
		doses[i-1]=aDataList[nt].pdXi[i];
	}
	/* unscale parameters */
	for(i = 1; i<=nparm; i++)
	{
		p[i] = p[i]/(pow(xmax,(i-1)));
		if(initial == Yes) aParmList[nt].pdIniP[i]= aParmList[nt].pdIniP[i]/(pow(xmax,(i-1)));
	}
	/* free malloc'd memory */
	free(doses);
	free(affect);
	free(nanim);
	free(parms);
	free(fitparms);
	free(Spec2);
	free(bind);
	FREE_DVECTOR(pBak, 1, nparm);
	FREE_DVECTOR(tmy, 1, nparm);
	FREE_DVECTOR(t, 1, nparm);
	FREE_DMATRIX(tmv, 1, nparm, 1, nparm);
	FREE_IVECTOR(SpBak, 1, nparm);

}  /* end: Multistage_fit */

/***************************************************************************
* Multistage_BMD -- Used to calculate the BMD and BMDL for Multistage model.
*  external: bmdparm
*  input:
*   nparm is the number of parameters
*   p[] is the vector of fitted parameters
*   gtol is a small positive number
*   iter is not used ????????
*   xlk is the log-likelihood for the fitted model
*   Rlevel[] is the vector of BMR's
*   Bmdl[] is the vector of BMDL's for the BMR's
*   Bmdu[] is the vector of BMDU's for the BMR's
*   BMD is the benchmark dose
*  output: BMD, Bmdl[], prints BMDL
****************************************************************************/
void Multistage_BMD (int nparm, double p[], double gtol, int *iter, double xlk,
					 double Rlevel[], double Bmdl[], double Bmdu[], double *BMD)
{
	double BMD_func(int nparm, double p[], double x, double gtol);
	double BMDL_func(int nparm, double xlk, double D, double p[], double x, 
		double gtol, int *is_zero);
	double BMDU_func(int nparm, double xlk, double D, double p[], double x, 
		double gtol, int *is_inf);
	double MaxLike(int nparm, double p[]);

	double   tol;
	double   xa,xb,fa,fb;
	double   D, Drange, poly;
	double   *pBak;
	int      i, j, k, kk, is_zero, is_inf;

	is_zero = is_inf = 0;
	pBak=DVECTOR(1, nparm);  /* allocate memory for vector */

	for(j=1; j<=nparm; j++)
		pBak[j]= p[j];          /* save the parameters */


	/**************** compute Chi-squared value  **************/
	/* IF ML is the value of the maximized log-likelihood,
	then ML - LR is the value
	log-likelihood at the BMDL or BMDU */
	if (bmdparm.level<0.5)  /* if confidence level < 0.5 */
		LR = 0.5*QCHISQ(1.0 - 2.0 * bmdparm.level, 1);
	else
		LR = 0.5*QCHISQ(2.0 * bmdparm.level-1.0,1);

	Rlevel[1] = BMR = bmdparm.effect;

	if (bmdparm.risk==1) /* if risk type is added */
	{
		if (1-BMR/(1-p[1]) <= 0)
		{ ERRORPRT("\nBMD computation failed, error: 1 - background <= BMR");}
		else
		{ck = -log( 1-BMR/(1-p[1]) ); }  /* Added risk */
	} /* end if */
	else
	{
		if (1-BMR <= 0)
		{ ERRORPRT("\nBMD computation failed, error: BMR >= 1"); }
		else
		{ck = -log(1-BMR); }              /* Extra risk */
	} /* end else */

	/*************** solve the BMD *********************/
	xa = D = 0.0;
	fa = -ck;  /* Note: ck > 0.0 */
	fb = fa;
	Drange = xmax;   /* xmax for new scale */
	k=1;
	while(k<300 && fb<0)
	{ /* see if BMD is larger than 3 times max dose level */
		fa=fb;
		xa=D;
		D=Drange*k/100.0;
		poly=p[nparm];
		for (j=nparm-1; j>=2; j--) poly = poly*D+p[j];
		poly=poly*D;
		fb= poly - ck; /* fb is the polynomial for computing BMD */
		k++;
	} /* end while */


	if(fb<0) ERRORPRT("BMD computation failed. BMD is larger than three times maximum input doses.");
	xb=D;
	tol=1.0E-6;

	/* compute the BMD */
	/* BMD_func works on log scale, so convert xa and xb to logs. */
	/* IF xa == 0, set xa = log(1e-300) */
	if (xa == 0.0) xa = -690.7755;
	else xa = log(xa);
	xb = log(xb);
	xb=zeroin(xa,xb,tol,BMD_func,nparm,p,tol);
	/* Now convert back to arithmetic scale */
	xa = exp(xa);
	*BMD = xb = exp(xb);

	/* print Benchmark Dose Computation */
	OUTPUT_BENCHMD(1, *BMD*scale);
	fflush(fp_out);
#ifndef RBMDS
	/* print the info that the plot file will use              */
	/* this info goes to the file with the 002 extension       */
	fprintf (fp_out2, "\n  RSL \t%f",bmdparm.effect*back1+back);
	fprintf (fp_out2, "\n  BMD \t%f",*BMD*scale);

	fprintf (fp_out2,"\n\n BMD_line");
	fprintf (fp_out2,"\n %f %f", (xmin-xmax/100), bmdparm.effect*back1+back );
	fprintf (fp_out2,"\n %f %f", *BMD*scale, bmdparm.effect*back1+back );
	fprintf (fp_out2,"\n %f %f", *BMD*scale, -0.1);
#endif

	BMD_lk = xlk;  /* get the lk at BMD */
	fb = -LR;
	/* calculate the BMDL */
	fa = BMDL_func(nparm, BMD_lk, xb, p, xa, tol, &is_zero);

	if (bmdl_bmr_flag == 1) {
		fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  No); /* computation failed */
		ERRORPRT("BMDL Calculation failed");
	}
#ifndef RBMDS
	fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  Yes); /* computation will succeed */
#endif
	Bmdl[1] = fa;
#ifdef MISC_OUT
	if (is_zero != 1)
		printf("           BMDL = %14.6g\n\n", (Bmdl[1]));
	else
		printf("           BMDL = 0.0\n\n");
#endif
	if (is_zero != 1)
		fprintf(fp_out, 
#ifndef RBMDS
		"            BMDL = %14.6g\n\n"
#else
		"            BMDL = %30.22g\n\n"
#endif
		, Bmdl[1]);
	else
		fprintf(fp_out,
		"            BMDL = 0.0\n\n");

#ifndef RBMDS
	fprintf(fp_out2,"\n BMDL \t%f", Bmdl[1]);
	fprintf (fp_out2,"\n\n BMDL_line");
	fprintf (fp_out2,"\n %f %f", Bmdl[1], -0.1);
	fprintf (fp_out2,"\n %f %f", Bmdl[1], bmdparm.effect*back1+back );
#endif
	/* calculate the BMDU */
	fa = BMDU_func(nparm, BMD_lk, xb, p, xa, tol, &is_inf);
	if (bmdu_bmr_flag == 1) {
#ifndef RBMDS
		fprintf (fp_out2, "\n\n BMDU_comput_ind %d",  No); /* computation failed */
#endif
		fprintf(fp_out, "BMDU calculation failed\n");
		fflush(fp_out);
	}
	Bmdu[1] = fa;
	if (is_inf != 1) {
		fprintf(fp_out, 
#ifndef RBMDS
			"            BMDU = %14.6g\n\n"
#else
			"            BMDU = %30.22g\n\n"
#endif
			, Bmdu[1]);
	} else {
		fprintf(fp_out,
			"            BMDU = Inf\n\n");
	}
	fflush(fp_out);
	switch ((bmdu_bmr_flag + is_inf)*10 + (bmdl_bmr_flag + is_zero)) {
  case 0: /* BMDL > 0, BMDU < Inf */
	  fprintf(fp_out,
		  "Taken together, (%-7.6g, %-7.6g) is a %-7.6g%% two-sided confidence\ninterval for the BMD\n\n",
		  Bmdl[1],Bmdu[1], 100*(2.0*bmdparm.level - 1.0));
  case 10: /* BMDL > 0, BMDU == Inf */
    fprintf(fp_out, 
#ifndef RBMDS
	    "Multistage Cancer Slope Factor =%14.6g\n\n"
#else
	    "Multistage Cancer Slope Factor =%30.22g\n\n"
#endif
	    , bmdparm.effect/Bmdl[1]);

	  break;
  case 1: /* BMDL == 0, BMDU < Inf */
  case 11: /* BMDL == 0, BMDU == Inf */
  default:
	  break;
	}
	if(bmdlCurve==Yes)
	{
		/****** calculate  Bmdl[] for each BMR ***********/
		for (k=2; k<=5;k++)
		{
			for(j=1; j<=nparm; j++)
				p[j]= pBak[j];          /* get the "old" p[] */
			/* recall: p[] can change in GETCL */
			if (k==2)
				Rlevel[k]=BMR=0.01;
			else if (k==3)
				Rlevel[k]=BMR=0.05;
			else
				Rlevel[k]= BMR = (k-2)*0.1;

			if (bmdparm.risk==1)
				ck = -log( 1-BMR/(1-p[1]) );   /* Added risk */
			else
			{
				ck = -log(1-BMR);                 /* Extra risk */
			} /* end else */

			/**** solve the BMD[] ********************************/
			xa = D = 0.0;
			fa = -ck;  /* Note: ck > 0.0 */
			fb = fa;
			Drange = xmax;   /* xmax for new scale */
			kk=1;
			while(kk<300 && fb<0)
			{ /* see if BMDL is out of range */
				fa=fb;
				xa=D;
				D=Drange*kk/100;
				poly=p[nparm];
				for (j=nparm-1; j>=2; j--) poly = poly*D+p[j];
				poly=poly*D;
				fb= poly - ck;
				kk++;
			} /* end while */

			if(fb<0)
			{
				Warning("Warning: BMDL is out of the three times range of dose for some BMR in BMDL curve computation.");
				Bmdl[k]=-1;
				continue;
			} /* end if */
			if (xa == 0.0) xa = -690.7755;
			else xa = log(xa);
			xb=log(D);
			/* calculate the BMD at the kth response level */
			xb=zeroin(xa,xb,tol,BMD_func,nparm,p,tol);
			xb = exp(xb);
			xa = exp(xa);
			/********* search for BMDL[] **************************/

			/* get the lk at BMD */
			/******** BMD_lk = MaxLike(nparm, p); **********/

			fb = -LR;

			/* compute the BMDL at the kth response level */
			/* note: before 8/18/99, xlk below was BMD_lk.
			This was causing problems with some data sets */
			fa = BMDL_func(nparm, xlk, xb, p, xa, tol, &is_zero);

			if (bmdl_bmr_flag == 1)
			{
				Bmdl[k]=-1;
				fprintf(fp_out, "\n BMDL curve computation failed for BMR = %f . \n The BMDL curve appearing in the graph may not be accurate.", BMR);
			}
			else
				Bmdl[k] = fa;

		}  /* end: for (k=2; k<=5; k++) */

		/* write BMDL curve info. into plot file */
#ifndef RBMDS
		fprintf (fp_out2, "\n\n BMDL_Curve_flag \t %d  \n smooth_opt  %d", bmdlCurve, smooth);
		fprintf (fp_out2,"\n\n BMDL_curve");
		fprintf (fp_out2,"\n 0.00000 %f", back);
		for (i=1;i<=5;i++)
			fprintf (fp_out2,"\n %f %f", Bmdl[i], Rlevel[i]*back1+back);
#endif
	}   /* end: if (bmdlCurve==Yes) */
	else for (k=2; k<=5;k++) /* if bmdlCurve = no */
		Bmdl[k] = Rlevel[k]= -1;
	FREE_DVECTOR(pBak, 1, nparm);
}  /* end: Multistage_BMD */

/*****************************************************************
* BMD_func -- used to compute the values of functions BMD_f at
*            the point x, given the parm p[] and number of parm.
*            (ck is another parameter).
*            This routine is called by Binary_root().
*  external: ck
*  input: n is the number of parameters
*         p[] is a vector of parameters
*         x is the natural log of a  dose level
*         gtol is useless here
*  output: value of the function
*****************************************************************/
double BMD_func(int n, double p[], double x, double gtol)
{
	double  poly, fx, D;
	int j;

	D = exp(x);
	poly=p[n];
	for (j=n-1; j>=2; j--) poly=poly*D+p[j];
	poly = poly*D;
	fx = poly - ck; /* ck = log(1-A) */
	return fx;
}

/*****************************************************************
* BMDL_func -- returns the lower confidence limit, BMDL.
*  external: Spec[]
*  input:
*   nparm is the number of parameters
*   xlk is the log-likelihood of the fitted model
*   Dose is the BMD or upper dose limit
*   pBak[] is the vector of fitted parameters
*   D is a lower dose limit
*   gtol is a small positive number (tolerance)
*  output: lower confidence limit
*****************************************************************/
double BMDL_func(int nparm, double xlk, double Dose, double pBak[],
				 double D, double gtol, int *is_zero)
{ 	/* ck , BMD_lk and LR are calculated in Multistage_BMD() */

	long int *Spec2, optite, nresm, *bind;
	long int which, polyord, temprisk;
	double fD, bmdl,  target, *parms, *parms2, xmin, x;
	int i, j, ii = 0;

	/* GETCL risks are switched as opposed to bmdparm.risk      */
	/* in Multistage.c  Make sure right risk is going to GETCL: */

	temprisk = bmdparm.risk + 1;

	/* Get the degree of polynomial */
	polyord = nparm - 1;

	/* allocate memory for vectors to be passed to FORTRAN code */
	/*    doses = (double *) malloc((size_t)(Nobs)*sizeof(double)); */
	/*    affect = (long int *) malloc((size_t)(Nobs)*sizeof(long int)); */
	/*    nanim = (long int *) malloc((size_t)(Nobs)*sizeof(long int)); */
	parms = (double *) malloc((size_t)(nparm)*sizeof(double));
	parms2 = (double *) malloc((size_t)(nparm)*sizeof(double));
	Spec2 = (long int *) malloc((size_t)(nparm)*sizeof(long int));
	bind = (long int *) malloc((size_t)(nparm)*sizeof(long int));

	xmin = aDataList[nt].pdXi[1];
	/* xmax is a global variable, which is initialized to 0. We do not
	 * re-initialize it to 0 so that it contains the max dose across
	 * all tumors.
	 */
	for (i=1;i<=Nobs;i++) /* obtain min and max dose levels */
	{
		x=aDataList[nt].pdXi[i];
		if (x < xmin) xmin = x;
		if (x > xmax) xmax = x;
	} /* end for */

	/** rescale Dose to be: 0 <= Dose <= 1 **/
	scale=xmax;
	Dose = Dose/xmax;

	which = 1;          /* Want a lower confidence limit */

	target = (xlk - LR);  /* The value we want the likelihood */
	/* at the BMDL to match             */

	/* Get appropriate values, offset by one */
	for (j=1;j<=nparm;j++)
	{
		parms[j-1]=pBak[j]*(pow(xmax,(j-1)));  /* get the "old" p[] */
		Spec2[j-1] = aParmList[nt].pnSpec[j];
	} /* end for */

	parms[0] = -log(1.0 - parms[0]);   /* FORTRAN code uses
									   this as the value
									   of p[0], so we must transform it*/


	getcl_(&which, &BMR, &Dose, &target, parms, Spec2, parms,
		&temprisk, &bmdl, parms2, &optite, &nresm, bind, is_zero);

	/* optite is a value that is passed back from GETCL which         */
	/* determines whether the optimization was completed successfully */
	/* If optite is less than 0, then it did not, and we want         */
	/* to try a different starting point and recompute                */

	if(optite < 0 )
	{
#ifdef MISC_OUT
		/* Warn user */
		fprintf(fp_out, "**** WARNING:  Completion code = %ld.  Optimum not found. Trying new starting point****\n\n", optite);
#endif
		/* Try up to 10 times if needed */
		for(ii = 0; ii < 10; ii++)
		{
			/* again, reparameterize p[0] */
			parms[0] = -log(1-pBak[1]);

			for(j = 2; j <= nparm; j++)  /* Get original values */
				parms[j-1] = pBak[j]*(pow(xmax,(j-1)));

			GetNewParms(parms, nparm);  /* Get a new starting point */

			/* Try again */
			getcl_(&which, &BMR, &Dose, &target, parms, Spec2, parms,
				&temprisk, &bmdl, parms2, &optite, &nresm, bind, is_zero);

			/* if optite >= 0, it is successful, and we can stop */
			if(optite >= 0)
				break;
#ifdef MISC_OUT
			/* otherwise, issues another warning, and continue trying */
			else
				fprintf(fp_out, "**** WARNING %d:  Completion code = %ld trying new start****\n\n", ii, optite);
#endif
		} /* end for */

	} /* end: if (optite < 0) */

	if(optite < 0 )
	{
#ifdef MISC_OUT
		/* Warn user */
		fprintf(fp_out, "**** WARNING:  Completion code = %ld.  Optimum not found. Trying new starting point****\n\n", optite);
#endif
		/* Try up to 10 times if needed */
		for(ii = 0; ii < 10; ii++)
		{
			/* again, reparameterize p[0] */
			parms[0] = -log(1-pBak[1]);

			for(j = 2; j <= nparm; j++)  /* Get original values */
				parms[j-1] = pBak[j]*(pow(xmax,(j-1)));

			GetMoreParms(parms, nparm);  /* Get a new starting point */

			/* Try again */

			getcl_(&which, &BMR, &Dose, &target, parms, Spec2, parms,
				&temprisk, &bmdl, parms2, &optite, &nresm, bind, is_zero);

			/* if optite >= 0, it is successful, and we can stop */
			if(optite >= 0)
				break;
#ifdef MISC_OUT
			/* otherwise, issues another warning, and continue trying */
			else
				fprintf(fp_out, "**** WARNING %d:  Completion code = %ld trying new start****\n\n", ii, optite);
#endif
		} /* end for */

	} /* end: if (optite < 0) */


	/* Let user know if no optimum was found */
	if(ii == 10)
	{
#ifdef MISC_OUT
		fprintf(fp_out, "\nWarning:  completion code still negative");
#endif
		fprintf(fp_out, "\nBMDL did not converge for BMR = %f\n", BMR);
		bmdl_bmr_flag = 1;
	} /* end if */

	/* unscale parameters */
	for(i = 0; i<=nparm-1; i++)
	{
		parms2[i] = parms2[i]/(pow(xmax,(i)));
#ifdef MISC_OUT
		printf("parms2[%d]=%f\n",i,parms2[i]);
#endif
	}
	bmdl = bmdl*xmax;

	fD = bmdl;      /* Get bmdl and...                   */
	/* free malloc'ed memory */
	free(parms);
	free(parms2);
	free(Spec2);
	free(bind);
	return fD;      /* return it to the calling function */
}   /* end: BMDL_func */

/*****************************************************************
* BMDU_func -- returns the upper confidence limit, BMDU.
*  external: Spec[]
*  input:
*   nparm is the number of parameters
*   xlk is the log-likelihood of the fitted model
*   Dose is the BMD or upper dose limit
*   pBak[] is the vector of fitted parameters
*   D is a lower dose limit
*   gtol is a small positive number (tolerance)
*  output: lower confidence limit
*****************************************************************/
double BMDU_func(int nparm, double xlk, double Dose, double pBak[],
				 double D, double gtol, int *is_inf)
{ 	/* ck , BMD_lk and LR are calculated in Multistage_BMD() */

	long int *Spec2, optite, nresm, *bind;
	long int which, polyord, temprisk;
	double fD, bmdu,  target, *parms, *parms2, xmin, x;
	int i, j, ii = 0;

	/* GETCL risks are switched as opposed to bmdparm.risk      */
	/* in Multistage.c  Make sure right risk is going to GETCL: */

	temprisk = bmdparm.risk + 1;

	/* Get the degree of polynomial */
	polyord = nparm - 1;

	/* allocate memory for vectors to be passed to FORTRAN code */
	parms = (double *) malloc((size_t)(nparm)*sizeof(double));
	parms2 = (double *) malloc((size_t)(nparm)*sizeof(double));
	Spec2 = (long int *) malloc((size_t)(nparm)*sizeof(long int));
	bind = (long int *) malloc((size_t)(nparm)*sizeof(long int));

	xmin = aDataList[nt].pdXi[1];
	/* xmax is a global variable, which is initialized to 0. We do not
	 * re-initialize it to 0 so that it contains the max dose across
	 * all tumors.
	 */
	for (i=1;i<=Nobs;i++) /* obtain min and max dose levels */
	{
		x=aDataList[nt].pdXi[i];
		if (x < xmin) xmin = x;
		if (x > xmax) xmax = x;
	} /* end for */

	/** rescale Dose to be: 0 <= Dose <= 1 **/
	scale=xmax;
	Dose = Dose/xmax;
	which = 2;          /* Want a upper confidence limit */

	target = (xlk - LR);  /* The value we want the likelihood */
	/* at the BMDU to match             */

	/* Get appropriate values, offset by one */
	for (j=1;j<=nparm;j++)
	{
		parms[j-1]=pBak[j]*(pow(xmax,(j-1)));  /* get the "old" p[] */
		Spec2[j-1] = aParmList[nt].pnSpec[j];
	} /* end for */

	parms[0] = -log(1.0 - parms[0]);   /* FORTRAN code uses
									   this as the value
									   of p[0], so we must transform it*/


	/* Call subroutine that calculates BMDL */

	getcl_(&which, &BMR, &Dose, &target, parms, Spec2, parms,
		&temprisk, &bmdu, parms2, &optite, &nresm, bind, is_inf);

	/* optite is a value that is passed back from GETCL which         */
	/* determines whether the optimization was completed successfully */
	/* If optite is less than 0, then it did not, and we want         */
	/* to try a different starting point and recompute                */
	if(optite < 0 )
	{
#ifdef MISC_OUT
		/* Warn user */
		fprintf(fp_out, "**** WARNING:  Completion code = %ld.  Optimum not found. Trying new starting point****\n\n", optite);
#endif
		/* Try up to 10 times if needed */
		for(ii = 0; ii < 10; ii++)
		{
			/* again, reparameterize p[0] */
			parms[0] = -log(1-pBak[1]);

			for(j = 2; j <= nparm; j++)  /* Get original values */
				parms[j-1] = pBak[j]*(pow(xmax,(j-1)));

			GetNewParms(parms, nparm);  /* Get a new starting point */

			/* Try again */
			getcl_(&which, &BMR, &Dose, &target, parms, Spec2, parms,
				&temprisk, &bmdu, parms2, &optite, &nresm, bind, is_inf);

			/* if optite >= 0, it is successful, and we can stop */
			if(optite >= 0)
				break;
#ifdef MISC_OUT
			/* otherwise, issues another warning, and continue trying */
			else
				fprintf(fp_out, "**** WARNING %d:  Completion code = %ld trying new start****\n\n", ii, optite);
#endif
		} /* end for */

	} /* end: if (optite < 0) */

	if(optite < 0 )
	{

		/* Warn user */
#ifdef MISC_OUT
		fprintf(fp_out, "**** WARNING:  Completion code = %ld.  Optimum not found. Trying new starting point****\n\n", optite);
#endif
		/* Try up to 10 times if needed */
		for(ii = 0; ii < 10; ii++)
		{
			/* again, reparameterize p[0] */
			parms[0] = -log(1-pBak[1]);

			for(j = 2; j <= nparm; j++)  /* Get original values */
				parms[j-1] = pBak[j]*(pow(xmax,(j-1)));

			GetMoreParms(parms, nparm);  /* Get a new starting point */

			/* Try again */

			getcl_(&which, &BMR, &Dose, &target, parms, Spec2, parms,
				&temprisk, &bmdu, parms2, &optite, &nresm, bind, is_inf);

			/* if optite >= 0, it is successful, and we can stop */
			if(optite >= 0)
				break;

#ifdef MISC_OUT
			/* otherwise, issues another warning, and continue trying */
			else
				fprintf(fp_out, "**** WARNING %d:  Completion code = %ld trying new start****\n\n", ii, optite);
#endif
		} /* end for */

	} /* end: if (optite < 0) */


	/* Let user know if no optimum was found */
	if(ii == 10)
	{
#ifdef MISC_OUT
		fprintf(fp_out, "\nWarning:  completion code still negative");
#endif
		fprintf(fp_out, "\nBMDU did not converge for BMR = %f\n", BMR);
		bmdu_bmr_flag = 1;
	} /* end if */

	/* unscale parameters */
	for(i = 0; i<=nparm-1; i++)
	{
		parms2[i] = parms2[i]/(pow(xmax,(i)));
#ifdef MISC_OUT
		printf("parms2[%d]=%f\n",i,parms2[i]);
#endif
	}
	bmdu = bmdu*xmax;
	fD = bmdu;      /* Get bmdl and...                   */
	/* free allocated memory */
	free(parms);
	free(parms2);
	free(Spec2);
	free(bind);
	return fD;      /* return it to the calling function */
}   /* end: BMDU_func */


void GetNewParms(double *p, int size)
/***********************************************************
* Given a vector of parameter values, and the number of
* parameters in that vector, this function will return
* new parameter values to restart the optimization if a "bad"
* completion code is returned from GETCL() or GETMLE, using
* a uniform random number centered at p[i]
*  external: restrict
*  input:
*   p[] is the vector of initial parameter estimates
*   size is the number of parameters
*  output: returns new starting values of the parameters
***********************************************************/
{
	int i;

	/* Find parameters by randomly selecting new parameters in */
	/* a uniform interval of p[i] +/- .0005*p[i]               */

	for(i = 0; i < size; i++)
	{
		p[i] = (p[i]*.001)*(rand()/32768.0) + p[i] - p[i]*.0005;

		/* If parameters are to be restricted to >= 0, check       */
		/* to make sure randomization did not force any below zero */
		/* if so, just try 0                                       */

		if(restrict == Yes && p[i] < 0)
			p[i] = 0;
	} /* end for */
}  /* end: GetNewParms */

/********************************************
* MaxLike computes the maximum log-likelihood
*  external: Nobs, Xi[], Yn[], Yp[]
*  input:
*   nparm is the number of parameters
*   p[] is the vector of fitted parameters
*  output: returns log-likelihood
*********************************************/
double MaxLike(int nparm, double p[])
{
	int i, j;
	double like, prob;

	printf("inside max like");
	like = 0.0;

	for(i = 1; i <= Nobs; i++)
	{
		prob = p[nparm];
		for(j = nparm-1; j >= 2; j--)
		{
			prob = aDataList[nt].pdXi[i]*prob + p[j];
		} /* end for */

		prob = aDataList[nt].pdXi[i]*prob;
		prob = p[1] + (1.0 - p[1])*(1 - exp(-prob));
		printf("prob= %le aDataList[nt].pdYp[%d] = %le \n",prob,i,aDataList[nt].pdYp[i]);
		if ((prob == 0.0) || (prob == 1.0))
		{
			if (aDataList[nt].pdYp[i] <= 0)
				like += 0;
			else
				like += -1.0e20;
		} /*end if */
		else
		{
			like += aDataList[nt].pdYp[i]*log(prob) + Yn[i]*log(1.0 - prob);
		} /* end else */
	} /* end for */

	return like;
} /* end: MaxLike */

/**********************************************************
*  P_RESPONSE - computes the probability of response
*   external: nparm
*   input:
*    q is a dose level
*    p[] is the vector of fitted parameters
*   output: probability of response
***********************************************************/
double P_RESPONSE(double q, double p[])
{
	double response;
	int k;

	response = p[nparm];
	for (k=nparm-1; k>=2; k--)
	{
		response = q*response + p[k];
	} /* end for */
	response = response*q;
	response = p[1] + (1-p[1])*(1-exp(-response));

	return response;

} /* end: P_RESPONSE */

void GetMoreParms(double *p, int size)
/***********************************************************
* Given a vector of parameter values, and the number of
* parameters in that vector, this function will return
* new parameter values to restart the optimization if a "bad"
* completion code is returned from GETCL() or GETMLE, using
* a uniform random number centered at p[i]
*  external: restrict
*  input:
*   p[] is the vector of initial parameter estimates
*   size is the number of parameters
*  output: returns new starting values of the parameters
***********************************************************/
{
	int i;

	/* Find parameters by randomly selecting new parameters in */
	/* a uniform interval of (0,1)                             */
	/* background is chosen randomly from (0, 0.01)            */

	p[0] = (.01)*rand()/32768.0;
	for(i = 1; i < size; i++)
	{
		p[i] = -1 + 2*rand()/32768.0;

		/* If parameters are to be restricted to >= 0, check       */
		/* to make sure randomization did not force any below zero */

		if(restrict == Yes && p[i] < 0)
			p[i] = -p[i];
	} /* end for */
}  /* end: GetMoreParms */

/************************************************************************ */
/* Model_DF -- Decrements the number of model parameters by the number of */
/*             parameters prespecified or "stuck" on a boundary of */
/*             parameter space. */
/* */
/*             Global: */
/*               nparm number of parameters */
/* */
/*             input: */
/*               bounded 1-based vector st bounded[i] means that the ith */
/*                       parameter is either fixed or bounded*/
/* */
/*             returns: */
/*               Model degrees of freedom - # fixed - #stuck or fixed */
/* */
/*             NOTE: Subtracting the number "stuck" on a boundary is */
/*                   incorrect, and needs to be revisited.*/
/************************************************************************ */

int Model_DF(int bounded[])
{
	int i, df;

	df = nparm;
	for (i = 1; i <= nparm; i++) df -=bounded[i];
	return df;
}


void GetMLEParms(double *p, int size)
/***********************************************************
* Given a vector of parameter values, and the number of
* parameters in that vector, this function will return
* new parameter values to restart the optimization if a "bad"
* completion code is returned from GETCL() or GETMLE, using
* a uniform random number centered at p[i]
*  external: restrict
*  input:
*   p[] is the vector of initial parameter estimates
*   size is the number of parameters
*  output: returns new starting values of the parameters
***********************************************************/
{
	int i;

	/* Find parameters by randomly selecting new parameters in */
	/* a uniform interval of (0,1)                             */
	/* background is chosen randomly from (0, 0.01)            */

	for(i = 0; i < size; i++)
	{
		p[i] = p[i] + .1*rand()/32768.0;

		/* If parameters are to be restricted to >= 0, check       */
		/* to make sure randomization did not force any below zero */

		if(restrict == Yes && p[i] < 0)
			p[i] = -p[i];
	} /* end for */
}  /* end: GetMLEParms */


void Copy_dvec(double *a, double *b, int size)
/***********************************************************
* Routine added by CVL 8/2007
* Given two double vectors (a and b) of size size,
* this routine copies a onto b
*   a[] is the vector to be copied from
*   b[] is the vector to be copied into
*   size is the number of elements in each of the vectors
*  output: returns copied vector
***********************************************************/
{
	int i;

	for(i = 0; i <= size; i++)
	{
		b[i] = a[i];
	}
}  /* end: Copy_dvec */


/***************************************************************************
* Multistage_ComboBMD -- Used to calculate the BMD and BMDL for combined 
*                         Multistage models A and B
*  external: bmdparm
*  input:
*   nparm is the number of parameters
*   p[] is the vector of fitted parameters
*   gtol is a small positive number
*   iter is not used ????????
*   xlk is the sum of log-likelihood for the fitted models (A + B)
*   Rlevel[] is the vector of BMR's
*   Bmdl[] is the vector of BMDL's for the BMR's
*   Bmdu[] is the vector of BMDU's for the BMR's
*   BMD is the benchmark dose
*  output: BMD, Bmdl[], prints BMDL
****************************************************************************/
void Multistage_ComboBMD (double gtol, int *iter, double xlk,
						  double Rlevel[], double Bmdl[], double Bmdu[], double *BMD)
{
	double BMD_func(int nparm, double p[], double x, double gtol);
	double BMDL_combofunc(double xlk, double Dose, double D, double gtol, int *is_zero);

	double   tol;
	double   xa,xb,fa,fb;
	double   D, Drange, poly;
	double   *cp, xlk2, crisk;
	int      j, k, is_zero, is_inf, cnparm, flag;

	is_zero = is_inf = 0;
	bmdl_bmr_flag = 0;
	/* GLN-01/30/2009 */
	//Get the greatest nparm
	cnparm = 0;
	for(j=1; j<=nT; j++)
	{
		if(aParmList[j].nParms > cnparm)
			cnparm = aParmList[j].nParms;
	}
	/* End of GLN-01/30/2009 */

	cp=DVECTOR(1, cnparm);    /* allocate memory for combined vector */
	for(j=1; j<=cnparm; j++)  cp[j] = 0;  /* zero out combined p vector*/

	/* GLN-01/30/2009 */
	//* Add all Model parameters to combined p
	for(j=1; j<=nT; j++)
	{
		cp[1] = cp[1] - log(1.0 - aParmList[j].pdParms[1]);
/*  3/18/09, BCA added above line and changed following for loop to start at 2: because cp needs to have beta[0] rather than "background" */
		for(k=2; k<=aParmList[j].nParms; k++)
		{
			cp[k] = cp[k] + aParmList[j].pdParms[k];
		}
	}
	/* End of GLN-01/30/2009 */

	/**************** compute Chi-squared value  **************/
	/* IF ML is the value of the maximized log-likelihood,
	then ML - LR is the value
	log-likelihood at the BMDL or BMDU */
	if (bmdparm.level<0.5)  /* if confidence level < 0.5 */
		LR = 0.5*QCHISQ(1.0 - 2.0 * bmdparm.level, 1);
	else
		LR = 0.5*QCHISQ(2.0 * bmdparm.level-1.0,1);

	Rlevel[1] = BMR = bmdparm.effect;

	bmdparm.risk=0;
	ck = -log(1-BMR);               /* Extra risk */

	/*************** solve the BMD *********************/
	xa = D = 0.0;
	fa = -ck;  /* Note: ck > 0.0 */
	fb = fa;
	Drange = cxmax;   /* combined xmax for new scale */
	k=1;
	while(k<300 && fb<0)
	{ /* see if BMD is larger than 3 times max dose level */
		fa=fb;
		xa=D;
		D=Drange*k/100.0;
		poly=cp[cnparm];
		for (j=cnparm-1; j>=2; j--) poly = poly*D+cp[j];
		poly=poly*D;
		fb= poly - ck; /* fb is the polynomial for computing BMD */
		k++;
	} /* end while */

	if(fb<0) ERRORPRT("BMD computation failed. BMD is larger than three times maximum input doses.");
	xb=D;
	tol=1.0E-15;

	/* compute the BMD */
	/* BMD_func works on log scale, so convert xa and xb to logs. */
	/* IF xa == 0, set xa = log(1e-300) */
	if (xa == 0.0) xa = -690.7755;
	else xa = log(xa);
	xb = log(xb);

	fprintf(fp_log,"\n\nIn Multistage_ComboBMD, Before zeroin() Line 2423, printing variables\n");
	fprintf(fp_log,"xa=%g, xb=%g, cxmax=%g\n\nPrinting values in cp:\n", xa, xb, cxmax);
	for(j=1; j<= cnparm; j++) fprintf(fp_log,"cp[%d]=%g\n", j, cp[j]);

	xb=zeroin(xa,xb,tol,BMD_func,cnparm,cp,tol);
	/* Now convert back to arithmetic scale */
	xa = exp(xa);
	*BMD = xb = exp(xb);

	/* print Benchmark Dose Computation */
	OUTPUT_BENCHMD(1, *BMD*scale);

	/* Convert to Exponential */
	for(j=1; j<=nT; j++)
	{
		aParmList[j].pdParms[1] = -log(1.0 - aParmList[j].pdParms[1]);
	}
	flag = 0; 
	xlk2 = ComboMaxLike(flag,xb,&crisk);
	fprintf(fp_log,"\n\n  Computed Combined Log-Likelihood %30.22g\n",xlk2);
	fprintf(fp_log,"\n  Combined BMD                     %24.10g  \n",xb);
	fprintf(fp_log,"\n  Combined Risk at BMD             %16.3g\n",crisk);

	/* convert from Exponential form */
	for(j=1; j<=nT; j++)
	{
		aParmList[j].pdParms[1] = 1.0 - exp(-aParmList[j].pdParms[1]); 
	}

	/*************** solve for the BMDL *********************/
	BMD_lk = xlk;  /* get the lk at BMD */
	fb = -LR;

	fa = BMDL_combofunc(BMD_lk, xb, xa, tol, &is_zero);

	if (bmdl_bmr_flag == 1) {
		fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  No); /* computation failed */
		ERRORPRT("BMDL Calculation failed");
	}

#ifndef RBMDS
	fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  Yes); /* computation will succeed */
#endif
	Bmdl[1] = fa;
#ifdef MISC_OUT
	if (is_zero != 1)
		printf("           BMDL = %14.6g\n\n", (Bmdl[1]));
	else
		printf("           BMDL = 0.0\n\n");
#endif
	if (is_zero != 1) {
	  fprintf(fp_out, 
#ifndef RBMDS
		  "            BMDL = %14.6g\n\n"
#else
		  "            BMDL = %30.22g\n\n"
#endif
		  , Bmdl[1]);
	  fprintf(fp_out, 
#ifndef RBMDS
		  "Multistage Cancer Slope Factor =%14.6g\n\n"
#else
		  "Multistage Cancer Slope Factor =%30.22g\n\n"
#endif
		  , bmdparm.effect/Bmdl[1]);
	}
	else {
	  fprintf(fp_out,"            BMDL = 0.0\n\n");
	  fprintf(fp_out,"Multistage Cancer Slope Factor = 0.0\n\n");
	}

	fflush(fp_out);
	/*  combo BMDL computation goes here */
	FREE_DVECTOR(cp, 1, cnparm);
}  /* end: Multistage_ComboBMD */



/*****************************************************************
*  Added by CVL 7/2007 - start block
* BMDL_combofunc -- returns the lower confidence limit, BMDL.
*  external: Spec[]
*  input:
*   nparm is the number of parameters for Model B
*   Anparm is the number of parameters for Model A
*   xlk is the log-likelihood of the fitted model
*   Dose is the BMD or upper dose limit
*   pBak[] is the vector of fitted parameters for Model B
*   pABak[] is the vector of fitted parameters for Model A
*   D is a lower dose limit
*   gtol is a small positive number (tolerance)
*  output: lower confidence limit
*****************************************************************/
double BMDL_combofunc(double xlk, double Dose, double D, double gtol, int *is_zero)
{ 	/* ck and LR are calculated in Multistage_ComboBMD() */

	long int optite, nresm, *bind, CBnparm;
	long int which, temprisk, lnParmMax;
	double fD, bmdl,  target, xmax, xlk2, xlk3, crisk;
	int i, j, nCall, k, nParmMax, ii = 0;
	
	/* GLN-01/15/2009 */
	double	*pdParms, *pdParms2, *pdVals, *adxmax;
	double *pdParmsBak;	/* LCO 03/2010 */
	int *piSpec2, nParms;
	fD = bmdl =  target = xmax = xlk2 = xlk3 = crisk = 0.0;
	optite = -5;
	nCall = 1;
	/*End GLN-01/15/2009 */

	/* GETCL2 risks are switched as opposed to bmdparm.risk      */
	/* in Multistage.c  Make sure right risk is going to GETCL: */

	temprisk = bmdparm.risk + 1;

	/* Get the degree of polynomial */
	//CBnparm = nparm + Anparm + 1;
	/* GLN-01/15/2009 also, get xmax, populate, print data to log */
	adxmax = (double *) malloc((size_t)nT*sizeof(double));
	for(i=0;i<nT;i++) adxmax[i] = 0.0;

	CBnparm = lnParmMax = 0;
	for(i = 1; i <= nT; i++)
	{
		xmax = aDataList[i].pdXi[1];
		CBnparm = CBnparm + aParmList[i].nParms;
		
		if(aParmList[i].nParms > lnParmMax)
			lnParmMax = aParmList[i].nParms;

		fprintf(fp_log,"\n\nIn BMDL_combofunc, Tumor %d data\n", i);
		fprintf(fp_log,"       DOSE      Inc    N \n");
		for(j = 1; j <= aDataList[i].nObs; j++)
		{
			if(aDataList[i].pdXi[j] > xmax)
				xmax = aDataList[i].pdXi[j];

			fprintf(fp_log,"%10.5g    %5.0f %5.0f \n",aDataList[i].pdXi[j], 
				aDataList[i].pdYp[j], aDataList[i].pdYp[j]+aDataList[i].pdYn[j]);
		}
		adxmax[i-1] = xmax;
	}
	CBnparm = CBnparm + 1;
	bind = (long int *) malloc((size_t)(CBnparm)*sizeof(long int));

	/** rescale all doses to be: 0 <= Dose <= 1 **/
	xmax = adxmax[0];
	for(i=1; i<nT; i++)
	{
	  if(adxmax[i] > xmax) xmax = adxmax[i];
	}
	scale = xmax;

	nParmMax = (int)lnParmMax;
	nParms = nT*nParmMax;
	pdParms = (double *) malloc((size_t)(nParms)*sizeof(double));
	pdParmsBak = (double *) malloc((size_t)(nParms)*sizeof(double));
	pdParms2 = (double *) malloc((size_t)(nParms)*sizeof(double));
	pdVals = (double *) malloc((size_t)(nParms)*sizeof(double));
	piSpec2 = (int *) malloc((size_t)(nParms)*sizeof(int));
	//Initialized parameters
	for(i= 0; i < nParms; i++)
	{
		pdParms[i] = 0.0;
		pdParms2[i] = 0.0;
		pdVals[i] = 0.0;
		piSpec2[i] = 0L;
	}

	fprintf(fp_log,"\n\nIn BMDL_combofunc, aParmList[i+1].pdParms[j+1](MLEs)\n");
	for(i = 1; i <= nT; i++)
	{
		fprintf(fp_log,"Tumor %d=>", i);
		/* LCO 03/29/2010 - Use actual # parms for tumor */
		for(j = 1; j <= aParmList[i].nParms; j++)
		{
		  fprintf(fp_log,"%10.5g\t", aParmList[i].pdParms[j]);
		}
		fprintf(fp_log,"\n");
	}

	k = -1;
	for(j = 1; j <= nParmMax; j++)
	{
	  /* Changed by lco 2009/12/09 so that order of parameters matches
	   * the order of data values.
		for(i = 1; i <= nT; i++)
	   */
		for(i = nT; i > 0; i--)
		{
			k++;
			if(j <= aParmList[i].nParms)
			{
				pdParms[k] = aParmList[i].pdParms[j];
				piSpec2[k] = aParmList[i].pnSpec[j];
			}
		}
	}

	fprintf(fp_log,"\n\nIn BMDL_combofunc, pdParms Values (MLEs, k=%d, nParms=%d)\n", k, nParms);
	i = 0;
	/* fprintf(fp_log,"\n"); */
	for(j = 1; j<=nT; j++)
		fprintf(fp_log,"    Tumor %d\t", j);
	fprintf(fp_log,"\n");

	for(k = 0; k < nParmMax; k++)
	{
		for(j = 1; j <= nT; j++)
			fprintf(fp_log,"%10.5g\t", pdParms[i++]);
		fprintf(fp_log,"\n");
	}
	/*End GLN-01/15/2009, get xmax, populate, print data to log  */

	j=0;
	fprintf(fp_log,"\nIn BMDL_combofunc, Tumor Starting Values\n");
	for(i = 0; i < nT; i++)
	{
		fprintf(fp_log,"Tumor %d => %10.5g\n", i+1, pdParms[i]);
	}

	fprintf(fp_log,"\nMaximum Dose  = %12.5g \n",xmax); 

	Dose = Dose/scale;
	fprintf(fp_log,"Scale  = %12.5g \n",scale); 

	which = 4;          /* Want a combined  lower confidence limit */

	target = (xlk - LR);  /* The value we want the likelihood */
	/* at the BMDL to match             */
	fprintf(fp_log,"Combined Loglikelihood         %30.22g \n",xlk); 
	fprintf(fp_log,"Target                         %30.22g \n",target); 

	/* LCO 03/29/2010 - Rework the following loop because it references
	 * invalid parameter values when the tumors have different
	 * polynomial degrees.
	 */
	k = -1;
	for(j = 1; j <= nParmMax; j++)
	{
	  /* Changed by lco 2009/12/09 so that order of parameters matches
	   * the order of data values.
		for(i = 1; i <= nT; i++)
	   */
		for(i = nT; i > 0; i--)
		{
		  int iParms = aParmList[i].nParms;

		  k++;
		  if (j <= iParms) {
		    pdParmsBak[k] =  aParmList[i].pdParms[j];
		    pdParms[k] = aParmList[i].pdParms[j] * (pow(scale,(j-1)));
		  } else {
		    pdParmsBak[k] = pdParms[k] = -99999;
		  }
		}
	}

	/* One more step for the background terms */
	for(i = 0; i < nT; i++) {
	  pdParms[i] = -log(1.0 - pdParms[i]);
	}
	/* GLN-01/15/2009 */
	fprintf(fp_log,"\n\n\nValues BEFORE call %d to getclmt_()", nCall);
	fprintf(fp_log,"BMR=%10.5g target=%10.5g\n",BMR, target);
	fprintf(fp_log,"bmdl=%10.5g optite=%ld", bmdl, optite);
	i = 0;
	fprintf(fp_log,"\n");
	for(j = 1; j<=nT; j++)
		fprintf(fp_log,"    Tumor %d\t\t\t", j);
	fprintf(fp_log,"\n");
	for(j = 1; j<=nT; j++)
		fprintf(fp_log,"Scaled | Unscaled\t\t");
	fprintf(fp_log,"\n");
	for(k = 0; k < nParmMax; k++)
	{
	  for(j = 1; j <= nT; j++) {
	    fprintf(fp_log,"%10.5g | %10.5g\t\t", pdParms[i], pdParmsBak[i]);
	    i++;
	  }
	  fprintf(fp_log,"\n");
	}
	fflush(fp_log);

	getclmt_(&which, &lnParmMax, &BMR, &Dose,
			 &target, pdParms, piSpec2,
			 pdParms, &temprisk, &bmdl,
			 pdParms2, &optite, &nresm,
			 bind, is_zero);

	fprintf(fp_log,"\n\nValues AFTER call %d to getclmt_()", nCall);
	fprintf(fp_log,"BMR=%10.5g target=%10.5g\n",BMR, target);
	fprintf(fp_log,"bmdl=%10.5g optite=%ld", bmdl, optite);
	i = 0;
	fprintf(fp_log,"\n          ");
	for(j = 1; j<=nT; j++)
		fprintf(fp_log,"    Tumor %d\t", j);
	fprintf(fp_log,"\n");

	for(k = 0; k < nParmMax; k++)
	{
		for(j = 1; j <= nT; j++)
			fprintf(fp_log,"%10.5g\t", pdParms[i++]);
		fprintf(fp_log,"\n");
	}
	fflush(fp_log);
	nCall++;

	/* optite is a value that is passed back from GETCL which         */
	/* determines whether the optimization was completed successfully */
	/* If optite is less than 0, then it did not, and we want         */
	/* to try a different starting point and recompute                */

	if(optite < 0 )
	{
#ifdef MISC_OUT
		/* Warn user */
		fprintf(fp_out, "**** WARNING:  Completion code = %ld.  Optimum not found. Trying new starting point****\n\n", optite);
#endif
		/* Try up to 10 times if needed */
		for(ii = 0; ii < 10; ii++)
		{
#if !0
		  GetNewParms2(pdParms, nParmMax);  /* Get a new starting point */
#else
			/* Get original values */
			k = -1;
			for(i = 1; i <= nT; i++)
			{
				for(j = 1; j <= nParmMax; j++)
				{
					k++;
					pdParms[k] = aParmList[i].pdParms[j] * (pow(scale,(j-1)));;
				}
			}
			GetNewParms2(pdParms, nParmMax);  /* Get a new starting point */

			/* again, reparameterize p[0] */
			for(i = 0; i < nT; i++)
			{
				pdParms[i] = -log(1-aParmList[i+1].pdParms[1]);
			}

#endif
			/* Try again */
	fprintf(fp_log,"\n\n\nValues BEFORE call %d to getclmt_()\n", nCall);
	fprintf(fp_log,"BMR=%10.5g target=%10.5g\n",BMR, target);
	fprintf(fp_log,"bmdl=%10.5g optite=%ld", bmdl, optite);
	i = 0;
	fprintf(fp_log,"\n");
	for(j = 1; j<=nT; j++)
		fprintf(fp_log,"    Tumor %d\t", j);
	fprintf(fp_log,"\n");
	for(k = 0; k < nParmMax; k++)
	{
		for(j = 1; j <= nT; j++)
			fprintf(fp_log,"%10.5g\t", pdParms[i++]);
		fprintf(fp_log,"\n");
	}
	fflush(fp_log);

			getclmt_(&which, &lnParmMax, &BMR, &Dose,
					 &target, pdParms, piSpec2,
					 pdParms, &temprisk, &bmdl,
					 pdParms2, &optite, &nresm,
					 bind, is_zero);

	fprintf(fp_log,"\n\nValues AFTER call %d to getclmt_()", nCall);
	fprintf(fp_log,"BMR=%10.5g target=%10.5g\n",BMR, target);
	fprintf(fp_log,"bmdl=%10.5g optite=%ld", bmdl, optite);
	i = 0;
	fprintf(fp_log,"\n");
	for(j = 1; j<=nT; j++)
		fprintf(fp_log,"    Tumor %d\t", j);
	fprintf(fp_log,"\n");

	for(k = 0; k < nParmMax; k++)
	{
		for(j = 1; j <= nT; j++)
			fprintf(fp_log,"%10.5g\t", pdParms[i++]);
		fprintf(fp_log,"\n");
	}
	nCall++;
	fflush(fp_log);

			/* if optite >= 0, it is successful, and we can stop */
			if(optite >= 0)
				break;
#ifdef MISC_OUT
			/* otherwise, issues another warning, and continue trying */
			else
				fprintf(fp_out, "**** WARNING %d:  Completion code = %ld trying new start****\n\n", ii, optite);
#endif
		} /* end for */

	} /* end: if (optite < 0) */

	if(optite < 0 )
	{
#ifdef MISC_OUT
		/* Warn user */
		fprintf(fp_out, "**** WARNING:  Completion code = %ld.  Optimum not found. Trying new starting point****\n\n", optite);
#endif
		/* Try up to 10 times if needed */
		for(ii = 0; ii < 10; ii++)
		{
			/* Get original values */
			k = -1;
			for(i = 1; i <= nT; i++)
			{
				for(j = 1; j <= nParmMax; j++)
				{
					k++;
					pdParms[k] = aParmList[i].pdParms[j] * (pow(scale,(j-1)));;
				}
			}
			GetMoreParms2(pdParms, nParmMax);  /* Get a new starting point */

			/* again, reparameterize p[0] */
			for(i = 0; i < nT; i++)
			{
				pdParms[i] = -log(1-aParmList[i+1].pdParms[1]);
			}

			/* Try again */
	fprintf(fp_log,"\n\n\nValues BEFORE call %d to getclmt_()", nCall);
	fprintf(fp_log,"BMR=%10.5g target=%10.5g\n",BMR, target);
	fprintf(fp_log,"bmdl=%10.5g optite=%ld", bmdl, optite);
	fprintf(fp_log,"\n");
	for(j = 1; j<=nT; j++)
		fprintf(fp_log,"    Tumor %d\t", j);
	fprintf(fp_log,"\n");
	i = 0;
	for(k = 0; k < nParmMax; k++)
	{
		for(j = 1; j <= nT; j++)
			fprintf(fp_log,"%10.5g\t", pdParms[i++]);
		fprintf(fp_log,"\n");
	}
	fflush(fp_log);

			getclmt_(&which, &lnParmMax, &BMR, &Dose,
					 &target, pdParms, piSpec2,
					 pdParms, &temprisk, &bmdl,
					 pdParms2, &optite, &nresm,
					 bind, is_zero);

	fprintf(fp_log,"\n\nValues AFTER call %d to getclmt_()", nCall);
	fprintf(fp_log,"BMR=%10.5g target=%10.5g\n",BMR, target);
	fprintf(fp_log,"bmdl=%10.5g optite=%ld", bmdl, optite);
	fprintf(fp_log,"\n");
	for(j = 1; j<=nT; j++)
		fprintf(fp_log,"    Tumor %d\t", j);
	fprintf(fp_log,"\n");
	i = 0;
	for(k = 0; k < nParmMax; k++)
	{
		for(j = 1; j <= nT; j++)
			fprintf(fp_log,"%10.5g\t", pdParms[i++]);
		fprintf(fp_log,"\n");
	}
	nCall++;
	fflush(fp_log);

			/* if optite >= 0, it is successful, and we can stop */
			if(optite >= 0)
				break;
#ifdef MISC_OUT
			/* otherwise, issues another warning, and continue trying */
			else
				fprintf(fp_out, "**** WARNING %d:  Completion code = %ld trying new start****\n\n", ii, optite);
#endif
		} /* end for */

	} /* end: if (optite < 0) */


	/* Let user know if no optimum was found */
	if(ii == 10)
	{
#ifdef MISC_OUT
		fprintf(fp_out, "\nWarning:  completion code still negative");
#endif
		fprintf(fp_out, "\nBMDL did not converge for BMR = %f\n", BMR);
		bmdl_bmr_flag = 1;
	} /* end if */

	double **ppdParms;
	ppdParms = DMATRIX (1, nT, 1, nParmMax);

	fprintf(fp_log,"******* pdParms2 Values ********\n");
	for(j = 1; j<=nT; j++)
		fprintf(fp_log,"    Tumor %d\t", j);
	fprintf(fp_log,"\n");

	k = -1;
	for(j = 1; j <= nParmMax; j++)
	{
		for(i = 1; i <= nT; i++)
		{
			k++;
			fprintf(fp_log,"%10.5g\t", pdParms2[k]);
			if(j <= aParmList[i].nParms)
			{
				 ppdParms[i][j] = pdParms2[k];
				 if(k < nT)
				 {
					ppdParms[i][j] = 1-exp(-pdParms2[k]);
				 }
			}
		}
		fprintf(fp_log,"\n");
	}

	/* unscale parameters  Note for use with ComboMaxLike use 1-nparm+1 locations
	in parm and Aparm here instead of 0-nparm             */
	//for(i = 1; i<=nT; i++)
	//{
	//	for(j = 2; j<= aParmList[i].nParms; j++)
	//		ppdParms[i-1][j-1] = ppdParms2[i-1][j-1]/(pow(scale,(j-1)));
	//}
	flag = 1; 
	bmdl = bmdl*scale;
	xlk3 = ComboMaxLike2(flag,bmdl,&crisk, ppdParms);
	/*  convert background parameters to non-exponential form for printing */
	//for(i = 1; i<= nT; i++)
	//{
	//	ppdParms[i-1][1] = 1-exp(-ppdParms[i-1][1]);
	//	ppdParms2[i-1][0] = 1-exp(-ppdParms2[i-1][0]);
	//	for(j = 1; j < aParmList[i].nParms; j++)
	//		ppdParms2[i-1][j] = 1-exp(-ppdParms2[i-1][j]);
	//}
	fprintf(fp_log,"\n\n Combined Log-Likelihood at BMDL (getcl) %30.22g\n",xlk2);
	fprintf(fp_log,"\n\n Combined Log-Likelihood at BMDL (combomaxlike) %30.22g\n",xlk3);
	fprintf(fp_log,"\n Combined BMDL                        %15.7g  \n",bmdl);
	fprintf(fp_log,"\n Combined Risk at BMDL                     %10.3g\n",crisk);

	fD = bmdl;      /* Get bmdl and...                   */
	/* free malloc'ed memory */
	free(pdParms);
	free(pdParms2);
	free(pdVals);
	free(piSpec2);
	free(bind);
	free(adxmax);
	return fD;      /* return it to the calling function */
}  /* end: BMDL_combofunc */

/***  Added by CVL 7/2007 - end block ***/

/** Added by CVL  8/2007 - start block **/
/*****************************************************************
* LogLik_Constant - computes and prints the constant for the log-likelihood
*  input:
*   Nobs is the number of observation
*  Yp[] is the positive dependent variable data array 
*  Yn[] is the negative dependent variable data array 
*****************************************************************/
double LogLik_Constant(int Nobs, double Yp[], double Yn[])
{
	double X, N, NX, con, xf, nf, nxf, val;
	double DLgama(double);
	int i;
	con = 0;
	for (i=1;i<=Nobs;i++)
	{
		X = Yp[i]+1;
		N = Yp[i] + Yn[i] + 1;
		NX = Yn[i] +1;
		xf = exp(DLgama(X));
		nf = exp(DLgama(N));
		nxf = exp(DLgama(NX));
		val = log(nf/(xf*nxf));
		con = con + val;
	} 
	fprintf(fp_out,"\n Log-likelihood Constant %30.22g \n",con);
	return con;
}


/*****************************************************************
* DLgama - double log gamma
*  input:
*   X - double value
*  Returns log(z-1)!
*****************************************************************/
double DLgama(double x)
{
	double r, az, z, y, dl;
	if (x == 0) return (1);
	z = x;
	r = 12.0;
	az = 1.0;
	for (;z-r < 0; z++) az = az * z;
	y = z * z;
	dl=0.91893853320467274-z+(z-0.5)*log(z)+z*
		(0.08333333333333333*y+0.0648322851140734)/(y*(y+
		0.811320754825416)+0.01752021563249146);
	dl = dl - log(az);
	return dl ;
}

/***********************************************************************
* ComboMaxLike computes the maximum log-likelihood for the combined 
*  Model A and Model B sets 
*  external: Nobs, Xi[], Yn[], Yp[]
*  external: Anobs, AXi[], AYn[], AYp[]
*  input:
*   flag - indicator of calc of log like at MLE or at lower bound 
*   dose - dose level for computing extra risk
*   crisk - computed extra risk
*   nparm is the number of parameters in Model B
*   Anparm is the number of parameters in Model A
*   p[] is the vector of fitted parameters for Model B
*   Ap[] is the vector of fitted parameters for Model A
*  output: returns log-likelihood
*  Note that this version of hte log-likelihood uses the expoential form 
*     or the background parameters
*************************************************************************/
double ComboMaxLike(int flag, double dose, double *crisk)
{
	int i, j, n, nt;
	double like, prob, bkg, pr, dSumParm1;

	prob = like = dSumParm1 = 0.0;

	for(n = 1; n <= nT; n++)
	{
		dSumParm1 = dSumParm1 + aParmList[n].pdParms[1];
		for(i = 1; i <= aDataList[n].nObs; i++)
		{
			prob = aParmList[n].pdParms[aParmList[n].nParms];

			if (n == 1 && flag==1 && aParmList[n].nParms==2) 
			{
				prob = aParmList[1].pdParms[2];
				for(nt = 2; nt <= nT; nt++)
					prob = prob - aParmList[nt].pdParms[2];
			}

			for(j = aParmList[n].nParms; j >= 1; j--)
			{
				if (n == 1 && flag==1 && j==2) 
				{
					pr = aParmList[1].pdParms[2];
					
					for(nt = 2; nt <= nT; nt++)
						pr = pr - aParmList[nt].pdParms[2];

					prob = aDataList[nt].pdXi[i]*prob + pr;				
				}
				else 
					prob = aDataList[n].pdXi[i]*prob + aParmList[n].pdParms[j];
			}
			prob = (1 - exp(-prob));
			if ((prob == 0) || (prob == 1))
			{
				if (aDataList[n].pdYp[i] <= 0 || aDataList[n].pdYn[i] <=0 ) 
					like += 0;
				else
				{
					if (prob == 1) 
						like += aDataList[n].pdYn[i]*(-710);
					else
						like += aDataList[n].pdYp[i]*(-710);
				}
			}
			else
			{
				like += aDataList[n].pdYp[i]*log(prob) + aDataList[n].pdYn[i]*log(1 - prob);
			} //if ((prob == 0) || (prob == 1))

		} //End for(i = 1; i <= aDataList[n].nObs; i++)
	} //End for(n = 1; n <= nT; n++)

	bkg = 1.0 - exp(-(dSumParm1));

	for(n = 1; n <= nT; n++)
	{
		prob = aParmList[n].pdParms[aParmList[n].nParms];
		if (n == 1 && flag==1 && aParmList[n].nParms==2) 
		{
			prob = aParmList[1].pdParms[2];

			for(nt = 2; nt <= nT; nt++)
				prob = prob - aParmList[nt].pdParms[2];
		}
		for(j = aParmList[n].nParms-1; j >= 1; j--)
		{
			if (n == 1 && flag==1 && j==2) 
			{
				pr = aParmList[1].pdParms[2];
				
				for(nt = 2; nt <= nT; nt++)
					pr = pr - aParmList[nt].pdParms[2];

				prob = dose*prob + pr;				
			}
			else
				prob = dose*prob + aParmList[n].pdParms[j];
		}
	}

	if (bkg == 1.0) 
		*crisk = 0.0;
	else
		*crisk = ((1.0 - exp(-(prob))) - bkg)/(1.0-bkg);
	return like;
} /* end: Combo MaxLike */
/** Added by CVL  8/2007 - end block **/

/** Added by GLN  02/11/2009 **/
double ComboMaxLike2(int flag, double dose, double *crisk, double **p)
{
	int i, j, n, nt, nObs, nParms;
	double like, prob, bkg, pr, dSumParm1;

	prob = like = dSumParm1 = 0.0;

	for(n = 1; n <= nT; n++)
	{
		dSumParm1 = dSumParm1 + p[n][1];
		nObs = aDataList[n].nObs;
		nParms = aParmList[n].nParms;
		for(i = 1; i <= nObs; i++)
		{
			prob = p[n][nParms];

			if (n == 1 && flag==1 && nParms==2) 
			{
				prob = p[1][2];
				for(nt = 2; nt <= nT; nt++)
					prob = prob - p[nt][2];
			}

			for(j = nParms; j >= 1; j--)
			{
				if (n == 1 && flag==1 && j==2) 
				{
					pr = p[1][2];
					
					for(nt = 2; nt <= nT; nt++)
						pr = pr - p[nt][2];

					prob = aDataList[n].pdXi[i]*prob + pr;				
				}
				else 
					prob = aDataList[n].pdXi[i]*prob + p[n][j];
			}
			prob = (1 - exp(-prob));
			if ((prob == 0) || (prob == 1))
			{
				if (aDataList[n].pdYp[i] <= 0 || aDataList[n].pdYn[i] <=0 ) 
					like += 0;
				else
				{
					if (prob == 1) 
						like += aDataList[n].pdYn[i]*(-710);
					else
						like += aDataList[n].pdYp[i]*(-710);
				}
			}
			else
			{
				like += aDataList[n].pdYp[i]*log(prob) + aDataList[n].pdYn[i]*log(1 - prob);
			} //if ((prob == 0) || (prob == 1))

		} //End for(i = 1; i <= nObs; i++)
	} //End for(n = 1; n <= nT; n++)

	bkg = 1.0 - exp(-(dSumParm1));

	for(n = 1; n <= nT; n++)
	{
		prob = p[n][aParmList[n].nParms];
		if (n == 1 && flag==1 && aParmList[n].nParms==2) 
		{
			prob = p[1][2];

			for(nt = 2; nt <= nT; nt++)
				prob = prob - p[nt][2];
		}
		for(j = aParmList[n].nParms; j >= 1; j--)
		{
			if (n == 1 && flag==1 && j==2) 
			{
				pr = p[1][2];
				
				for(nt = 2; nt <= nT; nt++)
					pr = pr - p[nt][2];

				prob = dose*prob + pr;				
			}
			else
				prob = dose*prob + p[n][j];
		}
	}

	if (bkg == 1.0) 
		*crisk = 0.0;
	else
		*crisk = ((1.0 - exp(-(prob))) - bkg)/(1.0-bkg);
	return like;
} /* end: Combo MaxLike 2*/


void GetNewParms2(double *p, int nMaxSize)
/***********************************************************
* Given a vector of parameter values, and the number of
* parameters in that vector, this function will return
* new parameter values to restart the optimization if a "bad"
* completion code is returned from GETCL() or GETMLE, using
* a uniform random number centered at p[i]
*  external: restrict
*  input:
*   p[] is the vector of initial parameter estimates
*   size is the number of parameters
*  output: returns new starting values of the parameters
***********************************************************/
{
	int i, k, j;

	/* Find parameters by randomly selecting new parameters in */
	/* a uniform interval of p[i] +/- .0005*p[i]               */
	k = -1;
	i = 0;
	for(j = 1; j <= nMaxSize; j++)
	{
		for(i = 1; i <= nT; i++)
		{
			k++;
			p[k] = (p[k]*.001)*(rand()/32768.0) + p[k] - p[k]*.0005;
			/* If parameters are to be restricted to >= 0, check       */
			/* to make sure randomization did not force any below zero */
			/* if so, just try 0                                       */
			if(restrict == Yes && p[k] < 0)
				p[k] = 0;
		}
	}

	//for(i = nMaxSize * nTum; i < nMaxSize; i++)
	//{
	//	if(i > ((nMaxSize * nTum)+aParmList[nTum].nParms)-1)
	//		break;
	//	p[i] = (p[i]*.001)*(rand()/32768.0) + p[i] - p[i]*.0005;
	//	/* If parameters are to be restricted to >= 0, check       */
	//	/* to make sure randomization did not force any below zero */
	//	/* if so, just try 0                                       */

	//	if(restrict == Yes && p[i] < 0)
	//		p[i] = 0;
	//} /* end for */
}  /* end: GetNewParms2 */

void GetMoreParms2(double *p, int nMaxSize)
/***********************************************************
* Given a vector of parameter values, and the number of
* parameters in that vector, this function will return
* new parameter values to restart the optimization if a "bad"
* completion code is returned from GETCL() or GETMLE, using
* a uniform random number centered at p[i]
*  external: restrict
*  input:
*   p[] is the vector of initial parameter estimates
*   size is the number of parameters
*  output: returns new starting values of the parameters
***********************************************************/
{
	int i, k, j;

	/* Find parameters by randomly selecting new parameters in */
	/* a uniform interval of (0,1)                             */
	/* background is chosen randomly from (0, 0.01)            */
	k = -1;
	for(j = 1; j <= nMaxSize; j++)
	{
		for(i = 1; i <= nT; i++)
		{
			k++;
			p[k] = -1 + 2*rand()/32768.0;
			/* If parameters are to be restricted to >= 0, check       */
			/* to make sure randomization did not force any below zero */
			if(restrict == Yes && p[i] < 0)
				p[k] = -p[k];
		}
	}
	for(i = 0; i < nT; i++)
	{
		p[i] = (.01)*rand()/32768.0;
	}
}  /* end: GetMoreParms2 */


/** End of Added by GLN  02/11/2009 **/



