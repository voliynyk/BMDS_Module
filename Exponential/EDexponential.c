/****************************************************************
*
* IMPORTANT NOTE:  The following variable is the version number for
*                  the current model.  THIS MUST BE CHANGED as
*				   important changes are made to the models.
*
*****************************************************************/
char Version_no[]="Exponential Model. (Version: 1.1;  Date: 08/15/2007)";
/*
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
********************************************************************/

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

extern void getmle_(int *ndoses, double doses[], double means[],
		    int nanimals[], double svar[], int *nparm,
		    double parms[], int fixed[], double fixedval[],
		    int *restrict,int *adverse, double parms2[],double *ll,
		    int *optite, int *nresm, int bind[],
		    int *model_type, int *flag);

//double BMDL_func(int nparm, double xlk, double Dose, double pBak[], double gtol);
void Exponential_fit(int i);
void GetMLEParms(double *p, int size);
void GetNewParms(double *p, int size);
void GetMoreParms(double *p, int size);

#define EPS 3.0e-8
#define DO_LOG		//Uncomment to generate log file code
#define NBR_OF_PARMS 6
#define NBR_OF_MODELS 4
#define dmax(x,y) ((x) < (y) ? (y) : (x))
#define dmin(x,y) ((x) > (y) ? (y) : (x))

/* control models */
/* it may be replaced by gpiSpecVector[i] later */
/* For DMNGB */
double gdMaxloglik=.0;
double gdMaxloglik_A3=.0;

int giErrorFlag; /* Error States from DMNGB */
int giErrorFlag_A3;

int    girestrict;    /* flag for restricting slope>=1 */

char    gacFileOut[128];  /*output temp file*/
char    gacFileOut2[128];
char	gacPltFileName[128];  /* file to pass to GnuPlot */

/*************************************************************
*   giDo_Log = true -> log file is made
*   giDo_Log = false -> no log
*************************************************************/
int     giDo_Log = true;     /*  Uncomment to activate switch for log file   */
int     giDo_PrintLog = false;     /*  option to print info to log file   */

char   *gaParm_Name[]={"alpha", "rho", "a", "b", "c", "d"};
char   gacFileName2[186], gacLogFile[186], *gcDot2;
int    giNbrParm_A3;

int	   gpimodtype;
double *gdBMD;
int *giRun;

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

int    *gpiNi;
double *gpdYm;
double *gpdYd;
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

int giPrintLL = 0; /* if 1, print a line after the anodev table with the */
/* log-likelihood to high precision.  Triggered by a */
/* negative value for gsExpoVars.iMaxIter*/

/* changing variable */
int giReplace, giBrat;
int gpiIter;


//Miscellaneous functions and structure added by Geoffrey, 03/24/2007

int giNmiss = 0;
void HeaderAndReadVar_2Out(char *dFileName, char* clocktime);	//Header and Read variables to out file
typedef struct exponent 
{
	int iSelect;
	char caUser_note[82];
	int iIn_Type;
	int iNbrObs_Total;	/* number of observations when type=1 or sample size when type=0 */
	int iIsInitial;
	int iSign;
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
	char caDoseName[20];
	char caNiName[20];
	char caMeanName[20];
	char caStdevName[20];
	char caResponseName[20];
	char caModelName[82];
	int iNTry;
} exponentialVars;
exponentialVars gsExpoVars;

typedef struct LikeInterest {
	char caModel[8];
	double dLogLikelihood;
	int    iDF;
	double dAIC_PValue;
} LikeInterestList;

double gdBMR;
double gdBMDL_r;
double gdb;

typedef struct BenchMark {
	int    iModelNbr;
	double dBMR;
	double dBMD;
	double dBMDL;
} BMarkList;

#define IBMarkList 4		//Nbr of Models in a BenchMark List
BMarkList *aBMarks;

double **gppdMLEs;
double *gpdLikelihoods;
double gdAdverseBMDEffect;

typedef enum {
	eM2=1, eM3, eM4, eM5
} eMs;

typedef enum {
	eAlpha=1, eRho, ea, eb, ec, ed
} eParms;

void DoInitValuesM2(double p[]);
void DoInitValuesM4(double p[]);
double getBMD23(double dBMD, double da, double db, double dd);
double getBMD45(double dBMD, double da, double db, double dc, double dd);
int Get_Linear_Trend2 (int N, double *doses, double *means, int *numi);

int READ_OBSDATA4V(int iNbr_Obs,double gpdXi[],int gpiNi[],double gpdYm[],double gpdYd[]);
int READ_OBSDATA2V(int iNTotal, double gpdxxi[], double gpdyyi[]);

/****************************************************************
** main--main function used to call Exponential mode fitting program.
Includes: biosubcc.c--common subfunction C program.	int iNtot=0;
*****************************************************************/
int main (int argc, char *argv[])
{
	int		i, j, jj, junk;	     /* iteration variable */
	int     iGroup;
	double  dlikeA1, dlikeA2, dlikeA3, dlikeR;     /* log likelihoods */
	double  *pdStDev;
	//double  *pdParms;	     /* parameter array */
	//double  *pdInitParms;	 /* Initial values parameter array */

	//AnaList *psAnasum;       /*information for ANONA analysis*/  
	char	junkname[82];
	char	caLong_path_name[300];

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

	fp_in=fopen(argv[1], "r");

	/* Check if input file is open, if not, print error message and exit */
	if (fp_in==NULL){
		fprintf(stderr,"Error in opening input file.\n");
		fprintf (stderr,"...Exited to system!\n");
		exit (1);
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
	}
	#endif

	/* begin reading input file from batch file (.(d) ext.) */
	fscanf(fp_in, "%s",gsExpoVars.caModelName );

	/* select is the Exponential model number */
	//fscanf(fp_in, "%d", &gsExpoVars.iSelect);
	gsExpoVars.iSelect = 2;
	fscanf(fp_in, "%[ ^\n]", gsExpoVars.caUser_note);
	fscanf(fp_in, "%[^\n]", gsExpoVars.caUser_note);
	fscanf(fp_in, "%s", junkname);
	fscanf(fp_in, "%s", junkname);
	fscanf(fp_in, "%d", &gsExpoVars.iIn_Type);

	/* gsExpoVars.iIn_Type=1 if input format is gpdXi, gpiNi, gpdYm, gpdYd. */
	/* gsExpoVars.iIn_Type=0 if input format is giNTotal, gpdXi, Y_ij. */

	if (gsExpoVars.iIn_Type==1)
		fscanf(fp_in, "%d", &gsExpoVars.iNbrObs_Total);
	else
		fscanf(fp_in, "%d", &gsExpoVars.iNbrObs_Total);

	fscanf(fp_in, "%d", &gsExpoVars.iSign);

	//if ((gsExpoVars.iSign != 1) || (gsExpoVars.iSign != -1))
	//	gsExpoVars.iSign = 0;

	giNbrParm_A3=gsExpoVars.iNbrObs_Total+2;

	fscanf(fp_in,"%d%lf%lf%d%d%d%d", &gsExpoVars.iMaxIter, &gsExpoVars.dRel_Conv, &gsExpoVars.dParm_Conv,
		&gsExpoVars.iBmdlCurve, &gsExpoVars.iBmdose, &gsExpoVars.iAppendIt, &gsExpoVars.iSmooth);

	fscanf(fp_in,"%d%lf%d%lf",&gsExpoVars.iBmr_Type,&gsExpoVars.dBmdEffect,&gsExpoVars.iCons_Var,&gsExpoVars.dBmdConfi_Level);

	#ifdef DO_LOG
		fscanf(fp_in,"%d", &gsExpoVars.iNTry);
	#endif

	#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"   \n\niIn_Type = %d\n", gsExpoVars.iIn_Type);
		fprintf(fp_log,"   iNbrObs_Total = %d\n", gsExpoVars.iNbrObs_Total);
		fprintf(fp_log,"   iSign = %d\n", gsExpoVars.iBmr_Type);

		fprintf(fp_log,"   iMaxIter = %d\n", gsExpoVars.iMaxIter);
		fprintf(fp_log,"   dRel_Conv = %g\n", gsExpoVars.dRel_Conv);
		fprintf(fp_log,"   dParm_Conv = %g\n", gsExpoVars.dParm_Conv);

		fprintf(fp_log,"   iBmr_Type = %d\n", gsExpoVars.iBmr_Type);
		fprintf(fp_log,"   dBmdEffect = %g\n", gsExpoVars.dBmdEffect);
		fprintf(fp_log,"   iCons_Var = %d\n", gsExpoVars.iCons_Var);
		fprintf(fp_log,"   dBmdConfi_Level = %g\n\n", gsExpoVars.dBmdConfi_Level);
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

	gdBMD = DVECTOR (1, NBR_OF_MODELS);
	giRun = IVECTOR (1, NBR_OF_MODELS);

	/*Check if memory allocation is a success*/
	if (!gppiSpecPara || !gppiInitPara || !gppdSpecPara || !gppdInitPara || !gpiInitialized || !gpiExp_Known) 
	{
		ERRORPRT ("Memory allocation failed in parameter matrices");
		#ifdef DO_LOG
		if (giDo_Log) 
		{
			fflush(fp_log);
			fclose(fp_log);
		}
		#endif

		CLOSE_FILES ();
		exit(1);
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
		}
		gpiInitialized[i]=0;
		gpiExp_Known[i] = 0;
		giRun[i] = 1;
		gdBMD[i] = 0.0;
	}

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

	/*Test if # of observation < # of parameters, then exit (instruction 7) */
	for(i=1; i<=NBR_OF_MODELS; i++)
	{
		if (gsExpoVars.iNbrObs_Total < (4 - gpiExp_Known[i]))
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
			ERRORPRT("Observation # < parameter # for Exponential model.");
			giRun[i] = 0;
			//#ifdef DO_LOG
			//if (giDo_Log) 
			//{
			//	fflush(fp_log);
			//	fclose(fp_log);
			//}
			//#endif
			//CLOSE_FILES();
			//exit (2);
		}
	}

	//for(i = 0; i < NBR_OF_PARMS; i++){
	//	if(gpiSpecVector[i] == 1)
	//		gpdIniParaVector[i] = 1;
	//}

	/*obtain observation data into gpdYp, gpdYn, gpdXi, Ls, Xg vectors*/
	if (gsExpoVars.iIn_Type == 1)
		fscanf(fp_in,"%s%s%s%s", gsExpoVars.caDoseName, gsExpoVars.caNiName, gsExpoVars.caMeanName, gsExpoVars.caStdevName);
	else
		fscanf(fp_in,"%s%s", gsExpoVars.caDoseName, gsExpoVars.caResponseName);

	/* The following code was commented by GLN 02/14/08 */
	//Change index to 2 in comparison, Bruce 01/28/2008
	//if(gpiSpecVector[2] == 1){ /* Determine whether it a heterogeneous or homogeneous
	//						   model for likelihood tests */
	//	if(pdParms[2] == 0) { //Change index to 2 in comparison, Bruce 01/28/2008
	//		//      gsExpoVars.iCons_Var = 1;
	//		psAnasum = ALVECTOR(1, 4);
	//	}
	//	else{
	//		//      gsExpoVars.iCons_Var = 0;
	//		psAnasum = ALVECTOR(1, 5);
	//	}
	//}
	//else {
	//	psAnasum = ALVECTOR(1, 5);
	//	//    gsExpoVars.iCons_Var = 0;
	//}
	/* End of commented code by GLN 02/14/08 */

	if (gsExpoVars.iIn_Type==1) 
	{
		gpdYm = DVECTOR (1, gsExpoVars.iNbrObs_Total);
		gpdYd = DVECTOR (1, gsExpoVars.iNbrObs_Total);
		gpdXi = DVECTOR (1, gsExpoVars.iNbrObs_Total);
		gpiNi = IVECTOR (1, gsExpoVars.iNbrObs_Total);
		pdStDev = DVECTOR (1, gsExpoVars.iNbrObs_Total);
		for(i =1; i <= gsExpoVars.iNbrObs_Total; i++)
			pdStDev[i] = 0.0;

		giNmiss = READ_OBSDATA4V(gsExpoVars.iNbrObs_Total, gpdXi, gpiNi, gpdYm, pdStDev);

		for(i =1; i <= gsExpoVars.iNbrObs_Total; i++)
			gpdYd[i] = pdStDev[i]*pdStDev[i];

		/* extern variable gsExpoVars.iNbrObs_Total has been changed. */
		gsExpoVars.iNbrObs_Total -= giNmiss;
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

		gsExpoVars.iNbrObs_Total -= giNmiss;  /* extern variable giNbr_Obs has been changed. */

		for (i=1; i<=gsExpoVars.iNbrObs_Total; i++) 
		{
			gpdXi[i] = gpdYsum[i] = gpdYm[i] = gpdYd[i] = 0;
			gpiNi[i] = 0;
		}

		Sort_2_By_Dose(gsExpoVars.iNbrObs_Total, xxi, yyi);  /*Sort gpdxxi and make appropriate changes to gpdyyi */

		for(i=1; i<=gsExpoVars.iNbrObs_Total; i++)
		{
			gpdxxi[i] = xxi[i+1];
			gpdyyi[i] = yyi[i+1];
		}

		iGroup = 0;

		for (i=1; i<=gsExpoVars.iNbrObs_Total; i++) {
			gpdXi[iGroup]=gpdxxi[i];
			gpiNi[iGroup] += 1;
			gpdYsum[iGroup] += gpdyyi[i];
			if(i < gsExpoVars.iNbrObs_Total) {
				if(gpdxxi[i] != gpdxxi[i + 1])
					iGroup +=1;
			}
		}
		gsExpoVars.iNbrObs_Total=iGroup+1;

		jj=1;

		for (i=1; i<=gsExpoVars.iNbrObs_Total; i++) {
			gpdYm[i] = gpdYsum[i]/gpiNi[i];
			for (j=1; j<=gpiNi[i]; j++)
			{
				gpdYd[i] += (gpdyyi[jj]-gpdYm[i])*(gpdyyi[jj]-gpdYm[i])/(gpiNi[i]-1);
				jj += 1;
			}
		}

		FREE_DVECTOR (xxi, 1, gsExpoVars.iNbrObs_Total);
		FREE_DVECTOR (yyi, 1, gsExpoVars.iNbrObs_Total);

	} // end else (gsExpoVars.iIn_Type==1)

	gpdScXi = DVECTOR(1, gsExpoVars.iNbrObs_Total);
	gpdScYd = DVECTOR(1, gsExpoVars.iNbrObs_Total);
	gpdScYm = DVECTOR(1, gsExpoVars.iNbrObs_Total);
	gdmaxYd = gpdYd[0]; 
	gdmaxYm = gpdYm[0];

	for (i=1; i<=gsExpoVars.iNbrObs_Total; i++){
		if(gpdYd[i] > gdmaxYd)
			gdmaxYd = gpdYd[i];
		if(gpdYm[i] > gdmaxYm)
			gdmaxYm = gpdYm[i];
	}

	for (i=1; i<=gsExpoVars.iNbrObs_Total; i++){
		gpdScYd[i] = gpdYd[i]/gdmaxYd;
		gpdScYm[i] = gpdYm[i]/gdmaxYm;
	}

	if (gsExpoVars.iIn_Type == 1) {
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
	}

	/* Print model and file information on output page */
	HeaderAndReadVar_2Out(argv[1], ctime(&ltime));

	//Compute for likelihoods
	dlikeA1=dlikeA2=dlikeA3=dvv=dlikeR=0.0;

	int iNtot=0;
	double dtemp=0.0;

	// Likelihood A1: Yij = Mu(i) + e(ij), Var{e(ij)} = Sigma^2
	for (i=1; i<=gsExpoVars.iNbrObs_Total; i++){
		iNtot += gpiNi[i];
	}

	for(i =1; i <= gsExpoVars.iNbrObs_Total; i++)
		dvv += ((gpiNi[i] - 1)*gpdYd[i]);

	dvv = dvv/iNtot;

	/* Compute Likelihood for model A1: Yij = Mu(i) + e(ij)
	Var(e(ij)) = Sigma^2 */

	dlikeA1 = - log(dvv)*iNtot/2 - iNtot/2;


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
	if (gsExpoVars.iCons_Var == 0) { // Parameters for fitting the model above
		#ifdef DO_LOG
		if (giDo_Log && giDo_PrintLog)	// Print values to log for investigation
		{
			fprintf(fp_log,"\n\nList of Data before AThree_fit\n");
			for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
			{
				fprintf(fp_log,"     %5.4g%7d%14.3g%14.3g\n", gpdXi[i], gpiNi[i], gpdYm[i], pdStDev[i]);
			}
		}
		#endif
		double *pdLKParms;

		pdLKParms = (double *) malloc((size_t) giNbrParm_A3*sizeof(double));
		for(i=1; i <= giNbrParm_A3; i++)
			pdLKParms[i]=0.0;

		/******* Fit the A3 model above *********/
		//AThree_fit(giNbrParm_A3, pdLKParms, EPS, &gpiIter, &dlikeA3);
		#ifdef DO_LOG
		if (giDo_Log && giDo_PrintLog)	// Print values to log for investigation
		{
			for (i=1; i<=giNbrParm_A3; i++) {
				fprintf(fp_log, "LKParms[%d] = %.8f\n", i, pdLKParms[i]);
			}
			fprintf(fp_log,"\n\nList of Data after AThree_fit\n");
			for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
			{
				fprintf(fp_log,"     %5.4g%7d%14.3g%14.3g\n", gpdXi[i], gpiNi[i], gpdYm[i], pdStDev[i]);
			}
		}
		#endif

		free(pdLKParms);
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
	dYbar = gpdYm[0] * gpiNi[0];

	for (i = 1; i <= gsExpoVars.iNbrObs_Total; i++)
		dYbar += gpdYm[i] * gpiNi[i];

	dYbar = dYbar / iNtot;

	dvv = gpdYd[0] * (gpiNi[0] - 1) + gpiNi[0]*(gpdYm[0] - dYbar)*(gpdYm[0] - dYbar);

	for (i = 1; i <= gsExpoVars.iNbrObs_Total; i++)
		dvv += gpdYd[i] * (gpiNi[i] - 1) + gpiNi[i]*(gpdYm[i] - dYbar)*(gpdYm[i] - dYbar);

	dvv = dvv / iNtot;

	dlikeR = -iNtot * (1.0 + log(dvv))/2.0;

	//Allocate memories for sMLEs + likelihood
	gppdMLEs = DMATRIX(1, NBR_OF_MODELS, 1, NBR_OF_PARMS);
	gpdLikelihoods = DVECTOR(1, NBR_OF_MODELS);

	if (!gppdMLEs || !gpdLikelihoods) 
	{
		ERRORPRT ("Memory allocation failed in sMLEs");
		#ifdef DO_LOG
		if (giDo_Log) 
		{
			fflush(fp_log);
			fclose(fp_log);
		}
		#endif

		CLOSE_FILES ();
		exit(1);
	}

	//INITIALIZE_DMATRIX(gppdMLEs, NBR_OF_MODELS, NBR_OF_PARMS);
	for(i=1; i<= NBR_OF_MODELS; i++)
	{
		for(j=1; j<=NBR_OF_PARMS; j++)
			gppdMLEs[i][j] = 0.0;

		gpdLikelihoods[i]=0.0;
	}


	/* Main loop to process Exponential Models*/
	// NOTE: This computation does take into account user parameter input yet.
	for( i=1; i <= NBR_OF_MODELS; i++)
	{
		#ifdef DO_LOG
		if (giDo_Log)	// Print values to log for investigation
		{
			fprintf(fp_log,"\n\nModel Run #%d, giRun[%d]=%d: ", i+1, i, giRun[i]);
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
		switch(i)
		{
		case (int)eM2:	//Model 2

			gsExpoVars.iSelect = i+1;

			/* GLN TODO: 2/18/08 -- change gppdSpecPara to gppdInitPara */
			if(gpiInitialized[i] == 0) //Not initialized generate initial value for Model 2
				DoInitValuesM2(gppdInitPara[i]);

			for(j=1; j<= NBR_OF_PARMS; j++)
				gppdMLEs[i][j] = gppdInitPara[i][j];

			gpimodtype = 3;
			/* TODO: */
			Exponential_fit(i);
			break;

		case (int)eM3:	//Model 3

			gsExpoVars.iSelect = i+1;

			/* GLN TODO: 2/18/08 -- change gppdSpecPara to gppdInitPara */
			for(j=1; j<= NBR_OF_PARMS; j++)
			{
				if(gpiInitialized[i] == 0)
				{
					gppdInitPara[i][j] = gppdMLEs[i-1][j];  //Initial parameter for Model 3, from Model 2
					gppdMLEs[i][j] = gppdMLEs[i-1][j];
				}
				else
					gppdMLEs[i][j] = gppdInitPara[i][j];
			}
			gpimodtype = 3;
			/* TODO: */
			Exponential_fit(i);
			break;

		case (int)eM4:	//Model 4

			gsExpoVars.iSelect = i+1;

			if(gpiInitialized[i] == 0)
				DoInitValuesM4(gppdInitPara[i]);

			/* GLN TODO: 2/18/08 -- change gppdSpecPara to gppdInitPara */
			for(j=1; j<= NBR_OF_PARMS; j++)
				gppdMLEs[i][j] = gppdInitPara[i][j];

			gpimodtype = 4;
			/* TODO: */
			Exponential_fit(i);
			break;

		default:		//Model 5

			gsExpoVars.iSelect = i+1;

			/* GLN TODO: 2/18/08 -- change gppdSpecPara to gppdInitPara */
			if(gpiInitialized[i] == 0)
			{
				for(j=1; j<= NBR_OF_PARMS; j++)
					gppdInitPara[i][j] = gppdMLEs[i-1][j];  //Initial parameter for Model 5, from Model 4
			}

			/* GLN TODO: 2/18/08 -- change gppdSpecPara to gppdInitPara */
			for(j=1; j<=NBR_OF_PARMS; j++)
				gppdMLEs[i][j] = gppdInitPara[i][j];

			gpimodtype = 4;
			/* TODO: */
			Exponential_fit(i);
			break;
		}

		#ifdef DO_LOG
		if (giDo_Log)	// Print values to log for investigation
		{
			fprintf(fp_log,"\n\nModel %d Initial Values: ", i+2);
			for(j=1; j <= NBR_OF_PARMS; j++)
			{
				fprintf(fp_log," %g", gppdInitPara[i][j]);
			}
			fprintf(fp_log,"\nModel %d MLEs: ", i+2);
			for(j=1; j <= NBR_OF_PARMS; j++)
			{
				fprintf(fp_log," %g", gppdMLEs[i][j]);
			}
		}
		#endif

		/* We need to change the OUT file format to output individual model. */
		//Print Initial Parameters and Estimates by Model
		//DoParamEstimates(gppdMLEs[i]);

		//DoDataEstimateInt(pdStDev);

		//DoLikehoods(dlikeA1, dlikeA2, dlikeA3, dlikeR);

		//DoBenchMark();

	}

	FREE_DVECTOR (gpdScYm, 1, gsExpoVars.iNbrObs_Total);
	FREE_DVECTOR (gpdScYd, 1, gsExpoVars.iNbrObs_Total);

	if (gsExpoVars.iIn_Type==1) 
	{
		FREE_DVECTOR (gpdYd, 1, gsExpoVars.iNbrObs_Total);
		FREE_DVECTOR (gpdYm, 1, gsExpoVars.iNbrObs_Total);
		FREE_IVECTOR (gpiNi, 1, gsExpoVars.iNbrObs_Total);
		FREE_DVECTOR (gpdXi, 1, gsExpoVars.iNbrObs_Total);
		FREE_DVECTOR (pdStDev, 1, gsExpoVars.iNbrObs_Total);
		FREE_IVECTOR (gpiNi, 1, gsExpoVars.iNbrObs_Total);
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

	FREE_IVECTOR (gpiInitialized, 1, NBR_OF_MODELS);
	FREE_IVECTOR (gpiExp_Known, 1, NBR_OF_MODELS);
	FREE_IVECTOR (giRun, 1, NBR_OF_MODELS);

#ifdef DO_LOG
	if (giDo_Log) 
	{
		fflush(fp_log);
		fclose(fp_log);
	}
#endif

	CLOSE_FILES ();

	return(0);
}

/**********************************************************
** READ_OBSDATA2V--used to read 2 column data in 2 vectors.
***********************************************************/
int READ_OBSDATA2V(int iNTotal, double gpdxxi[], double gpdyyi[])

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
				if (j==1)     gpdxxi[n]=value;
				if (j==2)     gpdyyi[n]=value;
			} 
			else m++;
		}
		if (m != 0)	   iNmiss++;
		else if (gpdxxi[n] < 0)      iNmiss++;
	}
	return iNmiss;
}
/**********************************************************
** READ_OBSDATA4V--used to read 4 column data in 4 vectors.
***********************************************************/
int READ_OBSDATA4V(int giNbr_Obs,double gpdXi[],int gpiNi[],double gpdYm[],double gpdYd[])
{
	int     iNmiss;          /*number of records with missing values*/
	int     i,j,n,m;        /*count and iteration control variables*/
	double  value;          /*temp variable*/

	iNmiss = 0;
	for(i=1;i<=giNbr_Obs;i++){
		n = i-iNmiss;
		m = 0;
		for (j=1;j<=4;j++){
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

void HeaderAndReadVar_2Out(char *dFileName, char* clocktime)
{
	strcpy(gacPltFileName,"");
	Output_Header(Version_no, dFileName, gacPltFileName, clocktime, gsExpoVars.caUser_note);

	/* output title and summary of input data */
	OUTPUT_TEXT("\n   The form of the response function by Model: ");
	OUTPUT_TEXT("      Model 2:     Y[dose] = a * exp{b * dose}");  //took out - sign on b, 01/28/08
	OUTPUT_TEXT("      Model 3:     Y[dose] = a * exp{b * dose^d}"); //took out - sign on b, 01/28/08
	OUTPUT_TEXT("      Model 4:     Y[dose] = a * [c – (c – 1)*exp{-b * dose}]");
	OUTPUT_TEXT("      Model 5:     Y[dose] = a * [c – (c – 1)*exp{-b * dose^d}]");

	OUTPUT_TEXT("\n      Model 2 is nested within Models 3 and 4.");
	OUTPUT_TEXT("      Model 3 is nested within Model 5.");
	OUTPUT_TEXT("      Model 4 is nested within Model 5.");

	if (gsExpoVars.iIn_Type == 1)
		fprintf(fp_out,"\n\n   Dependent variable = %s", gsExpoVars.caMeanName);
	else
		fprintf(fp_out,"\n\n   Dependent variable = %s", gsExpoVars.caResponseName);

	fprintf(fp_out,"\n   Independent variable = %s", gsExpoVars.caDoseName);
	if(gsExpoVars.iCons_Var == 0)
		fprintf(fp_out,"\n   The variance is to be modeled as Var(i) = exp(lalpha + log(mean(i)) * rho)");
	else
		fprintf(fp_out,"\n   rho is set to 0.\n   A constant variance model is fit.");

	fprintf(fp_out,"\n\n   Total number of dose groups = %d", gsExpoVars.iIn_Type==1 ? gsExpoVars.iNbrObs_Total-giNmiss : gsExpoVars.iNbrObs_Total-giNmiss);
	fprintf(fp_out,"\n   Total number of records with missing values = %d", giNmiss);
	fprintf(fp_out,"\n   Maximum number of iterations = %d", gsExpoVars.iMaxIter);
	fprintf(fp_out,"\n   Relative Function Convergence has been set to: %g", gsExpoVars.dRel_Conv);
	fprintf(fp_out,"\n   Parameter Convergence has been set to: %g\n", gsExpoVars.dParm_Conv);
}

void DoInitValuesM2(double p[])
{
	int iNtot=0;
	//, Nd=0, Nc=0;
	double y_bar=0., x_bar=0., stot=0.;
	double ssxy=0., ssx=0.;
	double sigma_bar=0., Ym_bar=0.;
	double ss_Ym=0., ss_YmSigma=0.;
	int i;

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

	// initial values for a, b
	for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i){

		x_bar+=(double)gpiNi[i]*gpdScXi[i]; //Changed gpdXi to gpdScXi, Bruce 01/28/08

		y_bar+=(double)gpiNi[i]*Slog(gpdYm[i]);
		iNtot+=gpiNi[i];
		stot+=((double)gpiNi[i]-1.)*gpdYd[i];
	}

	y_bar=y_bar/iNtot;
	x_bar=x_bar/iNtot;

	for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i){

		ssx +=(double)gpiNi[i]*(gpdScXi[i]-x_bar)*(gpdScXi[i]-x_bar); //Changed gpdXi to gpdScXi, Bruce 01/28/08

		ssxy+=(double)gpiNi[i]*(gpdScXi[i]-x_bar)*Slog(gpdYm[i]); //Changed gpdXi to gpdScXi, Bruce 01/28/08
	}

	// initial values for alpha and rho  
	if(gsExpoVars.iCons_Var==1)
	{
		p[(int)eAlpha]=(double)log(stot/iNtot); //added log, 01/28/08
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

	//Alpha, Rho, a, b, c, d
	//    0,   1, 2, 3, 4, 5
	p[(int)eb]=(double)gsExpoVars.iSign*ssxy/ssx;

	p[(int)ea]=exp(y_bar-p[(int)eb]*x_bar);  //a

	p[(int)ec]= 0.0;
	p[(int)ed]= 1.0;
}

void DoInitValuesM4(double p[])
{
	int i, j, iLowest=0, iGreatest=0;
	double *Pred;
	Pred = DVECTOR(1, gsExpoVars.iNbrObs_Total);

	double dLowest=0., dGreatest=0., dmG=0.,dma=0., dmc=0., dmb=0., dmult=0., dlltemp=0., dll=-Max_double, dxlk=0.0;
	dLowest = Max_double;
	dGreatest = -1;
	for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
	{
		if(dLowest >= gpdScXi[i])
		{
			dLowest = gpdScXi[i];
			iLowest = i;
		}
		if(dGreatest <= gpdScXi[i])
		{
			dGreatest = gpdScXi[i];
			iGreatest = i;
		}
	}

	if (gsExpoVars.iSign == 1)
	{
		for (j = 1; j <= 10; j++)
		{
			dmult = (j == 1 ? 1.05 : (j == 2 ? 2 : (j == 3 ? 5 : (j == 4 ? 10 : (j == 5 ? 20 : (j == 6 ? 50 : (j == 7 ? 100 : (j == 8 ? 200 : (j == 9 ? 500 : 1000)))))))));
			dma = 0.95*gpdYm[iLowest];
			dmG = dmult*gpdYm[iGreatest];

			//Compute for b:
			dGreatest	= 0.0;	//Borrowed as the numerator
			dLowest		= 0.0;	//Borrowed as the denominator

			for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
			{

				dGreatest += (double)( gpdScXi[i] * Slog( (dmG - gpdYm[i])/(dmG - dma) ) );
				dLowest += (double)(gpdScXi[i] * gpdScXi[i]);

			}

			dmb = -1.0*(double)(dGreatest/dLowest);
			dmc = dmG/dma;

			for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
			{
				Pred[i] = dma*(dmc-(dmc-1)*exp(-1.0*(dmb*gpdScXi[i])));
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
			}
		}
	}
	else
	{
		for (j=1; j <= 10; j++)
		{
			dmult = (j == 1 ? 1.05 : (j == 2 ? 2 : (j == 3 ? 5 : (j == 4 ? 10 : (j == 5 ? 20 : (j == 6 ? 50 : (j == 7 ? 100 : (j == 8 ? 200 : (j == 9 ? 500 : 1000)))))))));
			dmG = gpdYm[iLowest]/dmult;
			dma = 1.05 * gpdYm[iGreatest];

			//Compute for b:
			dGreatest	= 0.0;	//Borrowed as the numerator
			dLowest		= 0.0;	//Borrowed as the denominator

			for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
			{
				dGreatest += (double)( gpdScXi[i] * Slog( (dmG - gpdYm[i])/(dmG - dma) ) );
				dLowest += (double)(gpdScXi[i] * gpdScXi[i]);

			}

			dmb = -1.0*(double)(dGreatest/dLowest);
			dmc = dmG/dma;

			for(i=1; i<=gsExpoVars.iNbrObs_Total; ++i)
			{
				Pred[i] = dma*(dmc-(dmc-1)*exp(-1.0*(dmb*gpdScXi[i])));
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
			}
		}	
	}
	FREE_DVECTOR (Pred, 1, gsExpoVars.iNbrObs_Total);
}

void Exponential_fit(int i)
{
	double *fitparms, *fitparms2, *fitparms3, *beginp, ll, ll2, ll3;
	int *nAnim, *Spec2, *bind, *bind2, *bind3, optite, optite2, optite3;
	int nvar, signs, nparms, restr;
	int nresm, model_type, flag;
	int gmcnt, n, j;

	nvar = gsExpoVars.iNbrObs_Total;
	signs = gsExpoVars.iSign;
	nparms = NBR_OF_PARMS;
	restr = 0;
	nresm = 0;
	model_type = gpimodtype;

	fitparms = DVECTOR (1, nparms);
	fitparms2 = DVECTOR (1, nparms);
	fitparms3 = DVECTOR (1, nparms);
	beginp = DVECTOR (1, nparms);

	Spec2 = IVECTOR (1, nparms);
	nAnim = IVECTOR (1, nparms);
	bind = IVECTOR (1, nparms);
	bind2 = IVECTOR (1, nparms);
	bind3 = IVECTOR (1, nparms);

	optite = optite2 = optite3 = -5;

	ll = ll2 = ll3 = 0.0;

	for(n =1; n <= nparms; n++)
	{
		beginp[n] =  gppdInitPara[i][n];
		fitparms[n] = 0.0;
		fitparms2[n] = 0.0;
		fitparms3[n] = 0.0;
		Spec2[n] = gppiSpecPara[i][n];
	}
	
	for(n=1; n <= gsExpoVars.iNbrObs_Total; n++)
		nAnim[n] = gpiNi[n];

	fprintf(fp_log,"\n*******Spec Init Values, Model #%d **********\n", i);
	for(n =1; n <= nparms; n++)
		fprintf(fp_log,"Spec2[%d] = %d\n", n, Spec2[n]);

	//exp_getmle_(gsExpoVars.iNbrObs_Total, gpdScXi, gpdYm, gpiNi, gpdYd, NBR_OF_PARMS, 
	//	gppdInitPara[i], gppiSpecPara[i], gppdSpecPara[i], 0, gsExpoVars.iSign, 
	//	gppdMLEs[i], &gpdLikelihoods[i], &optite2, &nresm, bind2, gpimodtype, &flag);
	gmcnt = 1;
	flag = -1;

	/* This is the first call to getmle.  The parameters will be scaled */
	/* internally by donlp2 by their starting values.                   */
	while((optite2 < 0) || (optite2 > 2))
	{
		if(gmcnt > 1 && gmcnt < 6)
			GetNewParms(beginp+1, nparms);

		if(gmcnt > 5 && gmcnt < 10)
		{
			if (giDo_Log == true)
			{
				fprintf(fp_log,"\n\nbeginp before GetMoreParms. Call #%d\n", gmcnt);
				fprintf(fp_log,"beginp[1] = %12.5g  (alpha)\n",beginp[1]);
				fprintf(fp_log,"beginp[2] = %12.5g  (rho)\n",beginp[2]);
				fprintf(fp_log,"beginp[3] = %12.5g  (a)\n",beginp[3]);
				fprintf(fp_log,"beginp[4] = %12.5g  (b)\n",beginp[4]);
				fprintf(fp_log,"beginp[5] = %12.5g  (c)\n",beginp[5]);
				fprintf(fp_log,"beginp[6] = %12.5g  (d)\n",beginp[6]);
			}

			GetMoreParms(beginp+1, nparms);

			if (giDo_Log == true)
			{
				fprintf(fp_log,"\n\nbeginp after GetMoreParms. Call #%d\n", gmcnt);
				fprintf(fp_log,"beginp[1] = %12.5g  (alpha)\n",beginp[1]);
				fprintf(fp_log,"beginp[2] = %12.5g  (rho)\n",beginp[2]);
				fprintf(fp_log,"beginp[3] = %12.5g  (a)\n",beginp[3]);
				fprintf(fp_log,"beginp[4] = %12.5g  (b)\n",beginp[4]);
				fprintf(fp_log,"beginp[5] = %12.5g  (c)\n",beginp[5]);
				fprintf(fp_log,"beginp[6] = %12.5g  (d)\n",beginp[6]);
			}
		}

		if(gmcnt > 9 && gmcnt < 14)
			GetMLEParms(beginp+1, nparms);

		if (giDo_Log == true)
		{
			fprintf(fp_log,"\n***********Call #%d to getmle.****************\n",gmcnt);
			fprintf(fp_log,"(this call scales the parameters by their starting values in donlp2)\n");
			fprintf(fp_log,"nparms = %d\n",nparms);
			fprintf(fp_log,"flag = %d\n",flag);
			fprintf(fp_log,"Model = %d\n",i+1);
			fprintf(fp_log,"bmrtype = %d\n",gsExpoVars.iBmr_Type);
			fprintf(fp_log,"These are the parameters going into getmle.  The slope\n");
			fprintf(fp_log,"has been scaled by maxdose^power.\n");
			fprintf(fp_log,"beginp[1] = %12.5g  (alpha)\n",beginp[1]);
			fprintf(fp_log,"beginp[2] = %12.5g  (rho)\n",beginp[2]);
			fprintf(fp_log,"beginp[3] = %12.5g  (a)\n",beginp[3]);
			fprintf(fp_log,"beginp[4] = %12.5g  (b)\n",beginp[4]);
			fprintf(fp_log,"beginp[5] = %12.5g  (c)\n",beginp[5]);
			fprintf(fp_log,"beginp[6] = %12.5g  (d)\n",beginp[6]);
			fprintf(fp_log,"************************************************\n");
			fflush(fp_log);
		}

		getmle_(&nvar, gpdScXi+1, gpdYm+1, nAnim+1, gpdYd+1, &nparms, 
			beginp+1, Spec2+1, beginp+1, &restr, &signs, 
			fitparms2+1, &ll2, &optite2, &nresm, bind2+1, &model_type, &flag);

		if (giDo_Log == true)
		{
			fprintf(fp_log,"\n*******After Call #%d to getmle.**********\n",gmcnt);
			fprintf(fp_log,"optite2 = %d    (good optimum if 0<=optite2<=2)\n",optite2);
			fprintf(fp_log,"flag = %d\n",flag);
			fprintf(fp_log,"ll2 = %10.5g   (likelihood)\n",ll2);
			fprintf(fp_log,"nresm = %d    (no idea what this does)\n",nresm);
			fprintf(fp_log,"These are the parameters coming out of the first run\n");
			fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");
			fprintf(fp_log,"fitparms2[1] = %12.5g  (alpha)\n",fitparms2[1]);
			fprintf(fp_log,"fitparms2[2] = %12.5g  (rho)\n",fitparms2[2]);
			fprintf(fp_log,"fitparms2[3] = %12.5g  (a)\n",fitparms2[3]);
			fprintf(fp_log,"fitparms2[4] = %12.5g  (b)\n",fitparms2[4]);
			fprintf(fp_log,"fitparms2[5] = %12.5g  (c)\n",fitparms2[5]);
			fprintf(fp_log,"fitparms2[6] = %12.5g  (d)\n",fitparms2[6]);
			fprintf(fp_log,"*************************************************\n");
			fflush(fp_log);
		}

		gmcnt++;

		//FROM Pow_fit:
		//getmle_(&nvar, doses, means, nanim, svar, &nparms, 
		//	beginp, Spec2, beginp, &restr, &signs, 
		//	fitparms2, &ll2, &optite2, &nresm, bind2, &model_type, &flag);

		if(gmcnt == 14)
			break;
	}  //End of while((optite2 < 0) || (optite2 > 2)), first run.

	flag = 0;
	/* This is at most the 14th call to getmle.  The initial parameters will be */
	/* the fit parameters from the scaled (previous) call                       */
	if (giDo_Log == true)
	{
		fprintf(fp_log,"\n***********Call #%d to getmle.****************\n",gmcnt);
		fprintf(fp_log,"(This call doesn't scale the parameters by their starting values in\n");
		fprintf(fp_log," donlp2.  This call is using the parms coming out of the previous calls.)\n");
		fprintf(fp_log,"nparms = %d\n",nparms);
		fprintf(fp_log,"flag = %d\n",flag);
		fprintf(fp_log,"Model = %d\n",i+1);
		fprintf(fp_log,"bmrtype = %d\n",gsExpoVars.iBmr_Type);
		fprintf(fp_log,"These are the parameters going into getmle.  The slope\n");
		fprintf(fp_log,"has been scaled by maxdose^power.\n");
		fprintf(fp_log,"fitparms2[1] = %12.5g  (alpha)\n",fitparms2[1]);
		fprintf(fp_log,"fitparms2[2] = %12.5g  (rho)\n",fitparms2[2]);
		fprintf(fp_log,"fitparms2[3] = %12.5g  (a)\n",fitparms2[3]);
		fprintf(fp_log,"fitparms2[4] = %12.5g  (b)\n",fitparms2[4]);
		fprintf(fp_log,"fitparms2[5] = %12.5g  (c)\n",fitparms2[5]);
		fprintf(fp_log,"fitparms2[6] = %12.5g  (d)\n",fitparms2[6]);
		fprintf(fp_log,"************************************************\n");
		fflush(fp_log);
	}

	getmle_(&nvar, gpdScXi+1, gpdYm+1, nAnim+1, gpdYd+1, &nparms, 
		fitparms2+1, Spec2+1, beginp+1, &restr, &signs, 
		fitparms3+1, &ll3, &optite3, &nresm, bind3+1, &model_type, &flag);

	if (giDo_Log == true)
	{
		fprintf(fp_log,"\n*******After Call #%d to getmle.**********\n",gmcnt);
		fprintf(fp_log,"optite3 = %d    (good optimum if 0<=optite3<=2)\n",optite3);
		fprintf(fp_log,"flag = %d\n",flag);
		fprintf(fp_log,"ll3 = %10.5g   (likelihood)\n",ll3);
		fprintf(fp_log,"nresm = %d    (no idea what this does)\n",nresm);
		fprintf(fp_log,"These are the parameters coming out of the first run\n");
		fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");
		fprintf(fp_log,"fitparms3[1] = %12.5g  (alpha)\n",fitparms3[1]);
		fprintf(fp_log,"fitparms3[2] = %12.5g  (rho)\n",fitparms3[2]);
		fprintf(fp_log,"fitparms3[3] = %12.5g  (a)\n",fitparms3[3]);
		fprintf(fp_log,"fitparms3[4] = %12.5g  (b)\n",fitparms3[4]);
		fprintf(fp_log,"fitparms3[5] = %12.5g  (c)\n",fitparms3[5]);
		fprintf(fp_log,"fitparms3[6] = %12.5g  (d)\n",fitparms3[6]);
		fprintf(fp_log,"*************************************************\n");
		fflush(fp_log);
	}

	gmcnt++;

	/* This starts the loops without scaling in donlp2 */
	n = -1;
	for(flag=0; flag <= 1; flag++)
	{
		while((optite < 0) || (optite > 2))
		{
			if(n < 30 && n > -1)
			{
				if (optite != 3)
					GetNewParms(beginp+1, nparms);	/* Get a new starting point */
				else
				{
					for (j =1; j <= nparms; j++)
						beginp[j] = fitparms[j];
				}
			}

			if(n > 29 && n < 60)
			{
				if (optite != 3)
					GetMoreParms(beginp+1, nparms);		/* Get a new starting point */
				else
				{
					for (j =1; j <= nparms; j++)
						beginp[j] = fitparms[j];
				}
			}

			if(n > 59 && n < 91)
			{
				if(n == 60)
				{
					//for (i =1; i <= nparms; i++)
					//{
					//	parms[i] = p[i+1];
					//}
					//parms[3] = parms[3]*pow(maxdose, parms[4]);
					for (j =1; j <= nparms; j++)
					{
						if(j == (int)eb)
							beginp[j] = gppdInitPara[i][j] * pow(gdxmax, gppdInitPara[i][(int)ed]);
						else
							beginp[j] = gppdInitPara[i][j];
					}
				}
				else
				{
					if (optite != 3)
						GetMLEParms(beginp+1, nparms);	/* Get a new starting point */
					else
					{
						for (j =1; j <= nparms; j++)
							beginp[j] = fitparms[j];
					}
				}
			}

			if (giDo_Log == true)
			{
				fprintf(fp_log,"\n***********Call #%d to getmle.****************\n",gmcnt);
				fprintf(fp_log,"(this call doesn't scale the parameters by their starting values in donlp2)\n");
				fprintf(fp_log,"nparms = %d\n",nparms);
				fprintf(fp_log,"flag = %d\n",flag);
				fprintf(fp_log,"Model = %d\n",i+1);
				fprintf(fp_log,"bmrtype = %d\n",gsExpoVars.iBmr_Type);
				fprintf(fp_log,"These are the parameters going into getmle.  The slope\n");
				fprintf(fp_log,"has been scaled by maxdose^power.\n");
				fprintf(fp_log,"beginp[1] = %12.5g  (alpha)\n",beginp[1]);
				fprintf(fp_log,"beginp[2] = %12.5g  (rho)\n",beginp[2]);
				fprintf(fp_log,"beginp[3] = %12.5g  (a)\n",beginp[3]);
				fprintf(fp_log,"beginp[4] = %12.5g  (b)\n",beginp[4]);
				fprintf(fp_log,"beginp[5] = %12.5g  (c)\n",beginp[5]);
				fprintf(fp_log,"beginp[6] = %12.5g  (d)\n",beginp[6]);
				fprintf(fp_log,"************************************************\n");
				fflush(fp_log);
			}

			getmle_(&nvar, gpdScXi+1, gpdYm+1, nAnim+1, gpdYd+1, &nparms, 
				beginp+1, Spec2+1, beginp+1, &restr, &signs, 
				fitparms+1, &ll, &optite, &nresm, bind+1, &model_type, &flag);

			if (giDo_Log == true)
			{
				fprintf(fp_log,"\n*******After Call #%d to getmle.**********\n",gmcnt);
				fprintf(fp_log,"optite = %d    (good optimum if 0<=optite<=2)\n",optite);
				fprintf(fp_log,"flag = %d\n",flag);
				fprintf(fp_log,"ll = %10.5g   (likelihood)\n",ll);
				fprintf(fp_log,"nresm = %d    (no idea what this does)\n",nresm);
				fprintf(fp_log,"These are the parameters coming out of the first run\n");
				fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");
				fprintf(fp_log,"fitparms[1] = %12.5g  (alpha)\n",fitparms[1]);
				fprintf(fp_log,"fitparms[2] = %12.5g  (rho)\n",fitparms[2]);
				fprintf(fp_log,"fitparms[3] = %12.5g  (a)\n",fitparms[3]);
				fprintf(fp_log,"fitparms[4] = %12.5g  (b)\n",fitparms[4]);
				fprintf(fp_log,"fitparms[5] = %12.5g  (c)\n",fitparms[5]);
				fprintf(fp_log,"fitparms[6] = %12.5g  (d)\n",fitparms[6]);
				fprintf(fp_log,"*************************************************\n");
				fflush(fp_log);
			}

			gmcnt++;
			n++;

			if(n > 91)
			{
				n = 0;
				break;
			}

		} //End of while((optite < 0) || (optite > 2)) 3rd run.

		if ((optite >= 0) && (optite <= 2))
		{
			flag = 2;
			break;
		}

	}	/* end  for (flag=0; flag<=1; flag++) */

	/* This decides if the scaling model is better than the unscaled model */
	/* or not. */
	if ((optite2 >= 0) && (optite2 <= 2))
	{
		if ((ll2 > ll) || ((optite < 0) || (optite > 2)))
		{
			for (j =1; j <= nparms; j++)
			{
				fitparms[j] = fitparms2[j];
				bind[j] = bind2[j];
				gppdMLEs[i][j] = fitparms2[j];
			}
			optite = optite2;
			ll = ll2;
		}
	}

	if ((optite3 >= 0) && (optite3 <= 2))
	{
		if ((ll3 > ll) || ((optite < 0) || (optite > 2)))
		{
			for (j =1; j <= nparms; j++)
			{
				fitparms[j] = fitparms3[j];
				bind[j] = bind3[j];
				gppdMLEs[i][j] = fitparms3[j];
			}
			optite = optite3;
			ll = ll3;
		}
	}

	/* Let user know if no optimum was found */
	if ((optite < 0) || (optite > 2))
	{
		fprintf(fp_out, "\n\n!!! Warning:  optimum may not have been found for Model %d          !!!", i+2);
		fprintf(fp_out, "\n!!! Bad completion code in maximum likelihood optimization routine  !!!");
		fprintf(fp_out, "\n!!! Try choosing different initial values                           !!!\n\n");
		//*is_conv = -1;
		//exit(0);
	}
	gpdLikelihoods[i] = ll;

	if (gsExpoVars.iBmr_Type == 0) // Absolute difference in means
	{
		if (gsExpoVars.iSelect <= 3) // Models 2 and 3
		{
			//if ((gppdMLEs[i][3]+gdAdverseBMDEffect) <= 0)
			if ((gppdMLEs[i][(int)ea]+gdAdverseBMDEffect) <= 0)
			{
				giRun[i] = 0;
				//GLN TODO 2/19/08: have program print message in output where BMD would be that says:
				//"BMR choice is outside the predicted range of mean values; no BMD estimate is possible."
				// GLN TODO initialize and use vector giRun[1,4] starting with all elements = 1 (will run).
				// when giRun[i] = 0 then no further processing on that model.
				// this needs to be done before in main: where it says "instruction 7" where the decision about
				// too many parameters is done model by model: each time there are too many parameters, for 
				// that model (i), set giRun[i]=0 and do not do further processing on that model (e.g., do not
				// run exponential_fit) -- this means setting up if statements in some places conditioning on 
				// giRun[i] being equal to 1.
			}
			else // BMR within range
			{
				//gdBMD[i] = gsExpoVars.iSign*Log((gppdMLEs[i][3]+gdAdverseBMDEffect)/gppdMLEs[i][3]);
				//gdBMD[i] = pow(gdBMD[i], (1/gppdMLEs[i][6]));
				//gdBMD[i] = gdBMD[i]/gppdMLEs[i][4];

				//gdBMD[i] = gsExpoVars.iSign*Log((gppdMLEs[i][(int)ea]+gdAdverseBMDEffect)/gppdMLEs[i][(int)ea]);
				//gdBMD[i] = pow(gdBMD[i], (1/gppdMLEs[i][(int)ed]));
				//gdBMD[i] = gdBMD[i]/gppdMLEs[i][(int)eb];
				gdBMD[i] = getBMD23(gdBMD[i], gppdMLEs[i][(int)ea], gppdMLEs[i][(int)eb], gppdMLEs[i][(int)ed]);
			}
		}
		else  //Models 4 and 5
		{
			//if (gsExpoVars.dBmdEffect >= Abs(gppdMLEs[i][3] - gppdMLEs[i][3]*gppdMLEs[i][5]))
			if (gsExpoVars.dBmdEffect >= abs(gppdMLEs[i][(int)ea] - gppdMLEs[i][(int)ea]*gppdMLEs[i][(int)ec]))
			{
				giRun[i] = 0;
				//GLN TODO 2/19/08: have program print message in output where BMD would be that says:
				//"BMR choice is outside the predicted range of mean values; no BMD estimate is possible."
			}
			else // BMR within range
			{
				//gdBMD[i] = gppdMLEs[i][3] - gppdMLEs[i][3]*gppdMLEs[i][5] + gdAdverseBMDEffect;
				//gdBMD[i] = gdBMD[i]/(gppdMLEs[i][3] - gppdMLEs[i][3]*gppdMLEs[i][5]);
				//gdBMD[i] = pow((-Log(gdBMD[i])), (1/gppdMLEs[i][6]));
				//gdBMD[i] = gdBMD[i]/gppdMLEs[i][4];
				//gdBMD[i] = gppdMLEs[i][(int)ea] - gppdMLEs[i][(int)ea]*gppdMLEs[i][(int)ec] + gdAdverseBMDEffect;
				//gdBMD[i] = gdBMD[i]/(gppdMLEs[i][(int)ea] - gppdMLEs[i][(int)ea]*gppdMLEs[i][(int)ec]);
				//gdBMD[i] = pow((-log(gdBMD[i])), (1/gppdMLEs[i][(int)ed]));
				//gdBMD[i] = gdBMD[i]/gppdMLEs[i][(int)eb];
				gdBMD[i] = getBMD45(gdBMD[i], gppdMLEs[i][(int)ea], gppdMLEs[i][(int)eb], gppdMLEs[i][(int)ec], gppdMLEs[i][(int)ed]);
			}	
		}
	} // end of BMR type 0 (absolute difference)
	else if (gsExpoVars.iBmr_Type == 1)  // Standard deviation
	{
		gdBMD[i] = gdAdverseBMDEffect*sqrt(exp(gppdMLEs[i][1]+gppdMLEs[i][(int)eRho]*log(gppdMLEs[i][(int)ea])));
		if (gsExpoVars.iSelect <= 3) // Models 2 and 3
		{

			//if ((gppdMLEs[i][3]+gdBMD[i]) <= 0)
			if ((gppdMLEs[i][(int)ea]+gdBMD[i]) <= 0)
			{
				giRun[i] = 0;
				//GLN TODO 2/19/08: have program print message in output where BMD would be that says:
				//"BMR choice is outside the predicted range of mean values; no BMD estimate is possible."
			}
			else // BMR within range
			{
				//gdBMD[i] = gsExpoVars.iSign*Log((gppdMLEs[i][3]+gdBMD[i])/gppdMLEs[i][3]);
				//gdBMD[i] = pow(gdBMD[i], (1/gppdMLEs[i][6]));
				//gdBMD[i] = gdBMD[i]/gppdMLEs[i][4];
				gdBMD[i] = getBMD23(gdBMD[i], gppdMLEs[i][(int)ea], gppdMLEs[i][(int)eb], gppdMLEs[i][(int)ed]);
			}
		}
		else  //Models 4 and 5
		{
			//if (Abs(gdBMD[i]) >= Abs(gppdMLEs[i][3] - gppdMLEs[i][3]*gppdMLEs[i][5]))
			if (abs(gdBMD[i]) >= abs(gppdMLEs[i][(int)ea] - gppdMLEs[i][(int)ea]*gppdMLEs[i][(int)ec]))
			{
				giRun[i] = 0;
				//GLN TODO 2/19/08: have program print message in output where BMD would be that says:
				//"BMR choice is outside the predicted range of mean values; no BMD estimate is possible."
			}
			else // BMR within range
			{
				//gdBMD[i] = gppdMLEs[i][3] - gppdMLEs[i][3]*gppdMLEs[i][5] + gdBMD[i];
				//gdBMD[i] = gdBMD[i]/(gppdMLEs[i][3] - gppdMLEs[i][3]*gppdMLEs[i][5]);
				//gdBMD[i] = pow((-Log(gdBMD[i])), (1/gppdMLEs[i][6]));
				//gdBMD[i] = gdBMD[i]/gppdMLEs[i][4];
				gdBMD[i] = getBMD45(gdBMD[i], gppdMLEs[i][(int)ea], gppdMLEs[i][(int)eb], gppdMLEs[i][(int)ec], gppdMLEs[i][(int)ed]);
			}	
		}
	} // end of BMR type 1 (Standard deviation)
	else if (gsExpoVars.iBmr_Type == 2)  // Relative deviation
	{
		//eAlpha=0, eRho,	ea,	eb, ec, ed
		//	1		2		3	4	5	6
		//gdBMD[i] = gdAdverseBMDEffect*gppdMLEs[i][3];
		gdBMD[i] = gdAdverseBMDEffect*gppdMLEs[i][(int)ea];
		if (gsExpoVars.iSelect <= 3) // Models 2 and 3
		{

			//if ((gppdMLEs[i][3]+gdBMD[i]) <= 0)
			if ((gppdMLEs[i][(int)ea]+gdBMD[i]) <= 0)
			{
				giRun[i] = 0;
				//GLN TODO 2/19/08: have program print message in output where BMD would be that says:
				//"BMR choice is outside the predicted range of mean values; no BMD estimate is possible."
			}
			else // BMR within range
			{
				//gdBMD[i] = gsExpoVars.iSign*Log((gppdMLEs[i][3]+gdBMD[i])/gppdMLEs[i][3]);
				//gdBMD[i] = pow(gdBMD[i], (1/gppdMLEs[i][6]));
				//gdBMD[i] = gdBMD[i]/gppdMLEs[i][4];
				gdBMD[i] = getBMD23(gdBMD[i], gppdMLEs[i][(int)ea], gppdMLEs[i][(int)eb], gppdMLEs[i][(int)ed]);
			}
		}
		else  //Models 4 and 5
		{
			if (abs(gdBMD[i]) >= abs(gppdMLEs[i][3] - gppdMLEs[i][3]*gppdMLEs[i][5]))
			{
				giRun[i] = 0;
				//GLN TODO 2/19/08: have program print message in output where BMD would be that says:
				//"BMR choice is outside the predicted range of mean values; no BMD estimate is possible."
			}
			else // BMR within range
			{
				//gdBMD[i] = gppdMLEs[i][3] - gppdMLEs[i][3]*gppdMLEs[i][5] + gdBMD[i];
				//gdBMD[i] = gdBMD[i]/(gppdMLEs[i][3] - gppdMLEs[i][3]*gppdMLEs[i][5]);
				//gdBMD[i] = pow((-Log(gdBMD[i])), (1/gppdMLEs[i][6]));
				//gdBMD[i] = gdBMD[i]/gppdMLEs[i][4];
				gdBMD[i] = getBMD45(gdBMD[i], gppdMLEs[i][(int)ea], gppdMLEs[i][(int)eb], gppdMLEs[i][(int)ec], gppdMLEs[i][(int)ed]);
			}	
		}
	} // end of BMR type 2 (Relative deviation)
	else if (gsExpoVars.iBmr_Type == 3)  // Point
	{
		if (gsExpoVars.iSelect <= 3) // Models 2 and 3
		{

			if (gsExpoVars.iSign*(gppdMLEs[i][3]+gsExpoVars.dBmdEffect) > 0 || gsExpoVars.dBmdEffect < 0)
			{
				giRun[i] = 0;
				//GLN TODO 2/19/08: have program print message in output where BMD would be that says:
				//"BMR choice is outside the predicted range of mean values; no BMD estimate is possible."
			}
			else // BMR within range
			{ 
				//gdBMD[i] = gsExpoVars.iSign*Log(gsExpoVars.dBmdEffect/gppdMLEs[i][3]);
				//gdBMD[i] = pow(gdBMD[i], (1/gppdMLEs[i][6]));
				//gdBMD[i] = gdBMD[i]/gppdMLEs[i][4];
				gdBMD[i] = getBMD23(gdBMD[i], gppdMLEs[i][(int)ea], gppdMLEs[i][(int)eb], gppdMLEs[i][(int)ed]);
			}
		}
		else  // Models 4 and 5
		{
			if (gsExpoVars.dBmdEffect >= dmax(gppdMLEs[i][3],gppdMLEs[i][3]*gppdMLEs[i][5]) || gsExpoVars.dBmdEffect <= dmin(gppdMLEs[i][3],gppdMLEs[i][3]*gppdMLEs[i][5]))
			{
				giRun[i] = 0;
				//GLN TODO 2/19/08: have program print message in output where BMD would be that says:
				//"BMR choice is outside the predicted range of mean values; no BMD estimate is possible."
			}
			else // BMR within range
			{
				//gdBMD[i] = gsExpoVars.dBmdEffect - gppdMLEs[i][3]*gppdMLEs[i][5];
				//gdBMD[i] = gdBMD[i]/(gppdMLEs[i][3] - gppdMLEs[i][3]*gppdMLEs[i][5]);
				//gdBMD[i] = pow((-Log(gdBMD[i])), (1/gppdMLEs[i][6]));
				//gdBMD[i] = gdBMD[i]/gppdMLEs[i][4];
				gdBMD[i] = getBMD45(gdBMD[i], gppdMLEs[i][(int)ea], gppdMLEs[i][(int)eb], gppdMLEs[i][(int)ec], gppdMLEs[i][(int)ed]);
			}	
		}
	}  // end of BMR type 3 (point)
	else if (gsExpoVars.iBmr_Type == 4)  // Extra risk
	{
		if (gsExpoVars.iSelect <= 3) // Models 2 and 3
		{
			giRun[i] = 0;
			//GLN TODO 2/19/08: have program print message in output where BMD would be that says:
			//"Extra risk BMR type not applicable for Models 2 and 3; no BMD estimate is possible for those 			// models."

		}
		else  //Models 4 and 5
		{
			gdBMD[i] = gdAdverseBMDEffect*(gppdMLEs[i][3]*gppdMLEs[i][5]-gppdMLEs[i][3]);
			if (abs(gdBMD[i]) >= abs(gppdMLEs[i][3] - gppdMLEs[i][3]*gppdMLEs[i][5]))
			{
				giRun[i] = 0;
				//GLN TODO 2/19/08: have program print message in output where BMD would be that says:
				//"BMR choice is outside the predicted range of mean values; no BMD estimate is possible."
			}
			else // BMR within range
			{
				//gdBMD[i] = gppdMLEs[i][3] - gppdMLEs[i][3]*gppdMLEs[i][5] + gdBMD[i];
				//gdBMD[i] = gdBMD[i]/(gppdMLEs[i][3] - gppdMLEs[i][3]*gppdMLEs[i][5]);
				//gdBMD[i] = pow((-Log(gdBMD[i])), (1/gppdMLEs[i][6]));
				//gdBMD[i] = gdBMD[i]/gppdMLEs[i][4];
				gdBMD[i] = getBMD45(gdBMD[i], gppdMLEs[i][(int)ea], gppdMLEs[i][(int)eb], gppdMLEs[i][(int)ec], gppdMLEs[i][(int)ed]);
			}	
		}
	} // end of BMR type 4 (Extra)
	else  // invalid BMR type
	{
		for (i=1; i<=4; i++)
			giRun[i] = 0;
		// GLN TODO 2/19/08: print to screen and to output file (where BMDs and BMDLs would appear):
		// "An invalid choice has been made for the type of BMR.  Valid choices have integer flags between 0 and 4.
		//  No BMDs or BMDLs have been computed."

	}  //end of BMD calculations

	if (giDo_Log == true)
	{
		fprintf(fp_log,"\n*******End of Exponential_fit, Call #%d to getmle, Model #%d.**********\n",gmcnt, i+1);
		fprintf(fp_log,"giRun[%d] = %d\n", i, giRun[i]);
		fprintf(fp_log,"optite = %d    (good optimum if 0<=optite2<=2)\n",optite2);
		fprintf(fp_log,"flag = %d\n",flag);
		fprintf(fp_log,"ll = %10.5g   (likelihood)\n",ll);
		fprintf(fp_log,"nresm = %d    (no idea what this does)\n",nresm);
		fprintf(fp_log,"These are the parameters coming out of the first run\n");
		fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");
		fprintf(fp_log,"fitparms2[1] = %12.5g  (alpha)\n",fitparms2[1]);
		fprintf(fp_log,"fitparms2[2] = %12.5g  (rho)\n",fitparms2[2]);
		fprintf(fp_log,"fitparms2[3] = %12.5g  (a)\n",fitparms2[3]);
		fprintf(fp_log,"fitparms2[4] = %12.5g  (b)\n",fitparms2[4]);
		fprintf(fp_log,"fitparms2[5] = %12.5g  (c)\n",fitparms2[5]);
		fprintf(fp_log,"fitparms2[6] = %12.5g  (d)\n",fitparms2[6]);
		fprintf(fp_log,"*************************************************\n");
		fflush(fp_log);
	}

	FREE_DVECTOR (fitparms, 1, nparms);
	FREE_DVECTOR (fitparms2, 1, nparms);
	FREE_DVECTOR (fitparms3, 1, nparms);
	FREE_DVECTOR (beginp, 1, nparms);

	FREE_IVECTOR (Spec2, 1, nparms);
	FREE_IVECTOR (nAnim, 1, nparms);
	FREE_IVECTOR (bind, 1, nparms);
	FREE_IVECTOR (bind2, 1, nparms);
	FREE_IVECTOR (bind3, 1, nparms);

}	//end of Exponential_fit

double getBMD23(double dBMD, double da, double db, double dd)
{
	double dBMD23 = gsExpoVars.iSign*log((da+dBMD)/da);
	dBMD23 = pow(dBMD23, (1/dd));
	return (dBMD23/db);
}

double getBMD45(double dBMD, double da, double db, double dc, double dd)
{
	double dBMD45, dAddend, dsub;
	dBMD45 = dAddend = dsub = 0.0;
	switch(gsExpoVars.iBmr_Type)
	{
	case 0:
		dsub = da;
		dAddend = gdAdverseBMDEffect;
		break;
	case 3:
		dsub = gsExpoVars.dBmdEffect;
		dAddend = 0.0;
		break;
	default:	//1,2,4
		dsub = da;
		dAddend = dBMD;
		break;
	}
	dBMD45 = dsub - da*dc + dAddend;
	dBMD45 = dBMD45/(da - da*dc);
	dBMD45 = pow((-log(dBMD45)), (1/dd));
	return (dBMD45/db);
}

/***********************************************************
*	Given a vector of parameter values, and the number of
*	parameters in that vector, this function will return three
*	new parameter values to restart the optimization if a "bad"
*	completion code is returned from GETCL(), using a uniform
*	random number centered at p[i]
***********************************************************/
void GetNewParms(double *p, int size)
{
	int i;
	int restrict = 0;

	/* Find parameters by randomly selecting new parameters in
	/  a uniform interval of p[i] +/- .005*p[i] */
	int *Spec = gppiSpecPara[gsExpoVars.iSelect - 1];

	for (i =1; i <= size; i++)
	{
		if (Spec[i + 1] != 1)
		{
			p[i] = ((p[i] * .010) * (double)rand() / 32767.0) + p[i] * .995;
		}
	}

	/* If parameters are to be restricted, make sure restrictions
	/  are not violated */

	if (p[0] <= 0)
	{
		p[0] = -p[0];
	}

	if ((restrict == 1) && (p[4] <= 1))
	{
		p[4] = 1.00000001;
	}
	if ((restrict == 0) && (p[4] < 0))
	{
		p[4] = -p[4];
	}

}	/* end GetNewParms */

/***********************************************************
*	Given a vector of parameter values, and the number of
*	parameters in that vector, this function will return three
*	new parameter values to restart the optimization if a "bad"
*	completion code is returned from GETCL(), using a uniform
*	random number centered at p[i]
***********************************************************/
void GetMoreParms(double *p, int size)
{
	int restrict, i;
	restrict = 0;

	int *Spec = gppiSpecPara[gsExpoVars.iSelect-1];

	if (giDo_Log == true)
	{
		fprintf(fp_log,"\n*******Spec Values, Model #%d **********\n",gsExpoVars.iSelect);
		for(i =1; i <= size; i++)
			fprintf(fp_log,"gppiSpecPara[gsExpoVars.iSelect-2][%d] = %d\n", i, gppiSpecPara[gsExpoVars.iSelect-1][i]);

		for(i =1; i <= size; i++)
			fprintf(fp_log,"Spec[%d] = %d\n", i, Spec[i]);
		fflush(fp_log);
	}

	if (Spec[1] != 1)
		p[0] = 1*((double)rand() / (double) RAND_MAX);
	if (Spec[2] != 1)
		p[1] = 1*((double)rand() / (double) RAND_MAX);
	if (Spec[3] != 1)
		p[2] = 1*((double)rand() / (double) RAND_MAX);
	if (Spec[4] != 1)
	{
		if (p[3] >= 0)
			p[3] = 1*((double)rand() / (double) RAND_MAX);
		else
			p[3] = -1*((double)rand() / (double) RAND_MAX);
	} /* end if */
	if (Spec[5] != 1)
	{
		if (restrict == 1)
			p[4] = 1 + ((double)rand() / (double) RAND_MAX);
		else
			p[4] = ((double)rand() / (double) RAND_MAX);
	} /* end if */

}	/* end GetMoreParms */

/***********************************************************
*	Given a vector of parameter values, and the number of
*	parameters in that vector, this function will return three
*	new parameter values to restart the optimization if a "bad"
*	completion code is returned from GETCL(), using a uniform
*	random number centered at p[i]
***********************************************************/
void GetMLEParms(double *p, int size)
{
	int restrict = 0;
	
	int *Spec = gppiSpecPara[gsExpoVars.iSelect -1];

	if (Spec[1] != 1)
		p[0] = p[0] + 2*((double)rand() / (double) RAND_MAX);
	if (Spec[2] != 1)
		p[1] = p[1] + 1*((double)rand() / (double) RAND_MAX);
	if (Spec[3] != 1)
		p[2] = p[2] + 2*((double)rand() / (double) RAND_MAX);
	if (Spec[4] != 1)
	{
		if (p[3] >= 0)
			p[3] = p[3] + 2*((double)rand() / (double) RAND_MAX);
		else
			p[3] = p[3] -2*((double)rand() / (double) RAND_MAX);
	} /* end if */
	if (Spec[5] != 1)
	{
		if (restrict == 1)
			p[4] = p[4] + 1 + ((double)rand() / (double) RAND_MAX);
		else
			p[4] = p[4] + ((double)rand() / (double) RAND_MAX);
	} /* end if */

}	/* end GetMLEParms */
