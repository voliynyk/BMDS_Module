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
* Date: Feb 28, 2006
* 
* by C. Ahn
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

extern void getcl_(long int *which, long int *ndoses, double doses[], double means[],
				   long int nanimals[], double svar[], long int *nparm, double *bmr,
				   double *bmd, double *target, double parms[],
				   long int fixed[], double fixedval[], long int *risktype,
				   long int *restrict, double *bmdl, double parms2[],
				   long int *optite, long int *nresm, long int bind[], long int *adv,
				   long int *model_type, long int *flag);

double BMDL_func(int nparm, double xlk, double Dose, double pBak[], double gtol);
void GetOtherParms(double *p, int size);
void GetNewParms(double *p, int size);

#define EPS 3.0e-8
#define DO_LOG		//Uncomment to generate log file code

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

char   *gaParm_Name[]={"alpha", "rho", "a", "b", "c", "d"};
char   gacFileName2[186], gacLogFile[186], *gcDot2;
int    giNbrParm_A3;

int    *gpiSpecVector;	/* vector used to identify user input parm.
						will be 1 if user specified that parameter,
						0 otherwise */
int    *gpiIniSpVector;			/* user initial specified parameter vector */
double *gpdIniParaVector;	    /* user intialized parameter vector */
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

#define IBMarkList 4;		//Nbr of Models in a BenchMark List
BMarkList *aBMarks;

#define NBR_OF_PARMS 6

double *pdLogLikelihood;
double *pdModel2;
double *pdModel3;
double *pdModel4;
double *pdModel5;
double gdAdverseBMDEffect;

typedef enum {
	eLogM2=1, eLogM3, eLogM4, eLogM5
} eLogM;

typedef enum {
	eAlpha=1, eRho, ea, eb, ec, ed
} eParms;

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
	int     iNParm_Known, iExponential_Known;	 /* number of specified parameters */  
	double  dlikeA1, dlikeA2, dlikeA3, dlikeR;     /* log likelihoods */
	double  *pdStDev;
	double  *pdParms;	     /* parameter array */
	double  *pdInitParms;	 /* Initial values parameter array */

	AnaList *psAnasum;       /*information for ANONA analysis*/  
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

	/* Get initial default adverse direction by finding the general
	linear trend of the data */
	if ((gsExpoVars.iSign != 1) && (gsExpoVars.iSign != -1))	//Change || to && Bruce, 01/28/08 and moved from L817
	{
		gsExpoVars.iSign = Get_Linear_Trend(gsExpoVars.iNbrObs_Total, gpdXi, gpdYm, gpiNi);
		iTemp_sign = 0;
	}
	else {
		iTemp_sign = 9999;
	}/* end if ((gsExpoVars.iSign != 1) || (gsExpoVars.iSign != -1)) */


	// number of parameters
	//giNbr_Parm=6;

	giNbrParm_A3=gsExpoVars.iNbrObs_Total+2;

	// allocate momory for arrays
	pdParms = DVECTOR(1, NBR_OF_PARMS);
	if (!pdParms) 
	{
		ERRORPRT ("Memory allocation failed in pd0sParms");
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

	pdInitParms = DVECTOR(1, NBR_OF_PARMS);
	if (!pdInitParms) 
	{
		ERRORPRT ("Memory allocation failed in pdInitParms");
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

	gpiSpecVector = IVECTOR(1, NBR_OF_PARMS);
	gpdIniParaVector = DVECTOR(1, NBR_OF_PARMS);
	gpiIniSpVector = IVECTOR(1, NBR_OF_PARMS);

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

	gdAdverseBMDEffect = gsExpoVars.dBmdEffect * gsExpoVars.iSign;

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

	/* obtain user input parameters */
	READ_PARAMETERS(NBR_OF_PARMS, pdParms);
	#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\nList of Data after\nREAD_PARAMETERS(NBR_OF_PARMS, pdParms)\n");
		for(i=1; i <= NBR_OF_PARMS; i++)
		{
			fprintf(fp_log,"     pdParms[%d] = %g\n", i, pdParms[i]);
		}
	}
	#endif

	FILL_SPECVECTOR(NBR_OF_PARMS, pdParms, gpiSpecVector);
	#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\nList of Data after\nFILL_SPECVECTOR(NBR_OF_PARMS, pdParms, gpiSpecVector)\n");
		for(i=1; i <= NBR_OF_PARMS; i++)
		{
			fprintf(fp_log,"     gpiSpecVector[%d] = %d\n", i, gpiSpecVector[i]);
		}
	}
	#endif

	iNParm_Known = COUNT_SPECVECTOR(NBR_OF_PARMS, gpiSpecVector);
	iExponential_Known = 0;

	for(i=3; i<=6; ++i){
		if(gpiSpecVector[i]==1) ++iExponential_Known;
	}

	//??? 
	giBrat=Yes;

	if (gpiSpecVector[1]==1 && pdParms[1]<EPS) giBrat=No;

	/* obtain user input initial parameters values */
	fscanf(fp_in,"%d", &gsExpoVars.iIsInitial);
	READ_PARAMETERS(NBR_OF_PARMS, gpdIniParaVector);
	FILL_SPECVECTOR(NBR_OF_PARMS, gpdIniParaVector, gpiIniSpVector);
	#ifdef DO_LOG
	if (giDo_Log)	// Print values to log for investigation
	{
		fprintf(fp_log,"\n\nList of Data after\nFILL_SPECVECTOR(NBR_OF_PARMS, gpdIniParaVector, gpiIniSpVector)\n");
		for(i=1; i <= NBR_OF_PARMS; i++)
		{
			fprintf(fp_log,"     gpiIniSpVector[%d] = %d\n", i, gpiIniSpVector[i]);
		}
	}
	#endif

	for(i = 0; i < NBR_OF_PARMS; i++){
		if(gpiSpecVector[i] == 1)
			gpdIniParaVector[i] = 1;
	}

	/*obtain observation data into gpdYp, gpdYn, gpdXi, Ls, Xg vectors*/
	if (gsExpoVars.iIn_Type == 1)
		fscanf(fp_in,"%s%s%s%s", gsExpoVars.caDoseName, gsExpoVars.caNiName, gsExpoVars.caMeanName, gsExpoVars.caStdevName);
	else
		fscanf(fp_in,"%s%s", gsExpoVars.caDoseName, gsExpoVars.caResponseName);

	//Change index to 2 in comparison, Bruce 01/28/2008
	if(gpiSpecVector[2] == 1){ /* Determine whether it a heterogeneous or homogeneous
							   model for likelihood tests */
		if(pdParms[2] == 0) { //Change index to 2 in comparison, Bruce 01/28/2008
			//      gsExpoVars.iCons_Var = 1;
			psAnasum = ALVECTOR(1, 4);
		}
		else{
			//      gsExpoVars.iCons_Var = 0;
			psAnasum = ALVECTOR(1, 5);
		}
	}
	else {
		psAnasum = ALVECTOR(1, 5);
		//    gsExpoVars.iCons_Var = 0;
	}
	
	if (gsExpoVars.iIn_Type==1) {
		gpdYm = DVECTOR(1, gsExpoVars.iNbrObs_Total);
		gpdYd = DVECTOR(1, gsExpoVars.iNbrObs_Total);
		gpdXi = DVECTOR(1, gsExpoVars.iNbrObs_Total);
		gpiNi = IVECTOR(1, gsExpoVars.iNbrObs_Total);
		pdStDev = DVECTOR(1, gsExpoVars.iNbrObs_Total);
		giNmiss = READ_OBSDATA4V(gsExpoVars.iNbrObs_Total, gpdXi, gpiNi, gpdYm, pdStDev);
		for(i = 1; i <= gsExpoVars.iNbrObs_Total; i++)
			gpdYd[i] = pdStDev[i]*pdStDev[i];

		/* extern variable gsExpoVars.iNbrObs_Total has been changed. */
		gsExpoVars.iNbrObs_Total -= giNmiss;
	}
	else {
		gpdxxi = DVECTOR(1, gsExpoVars.iNbrObs_Total);
		gpdyyi = DVECTOR(1, gsExpoVars.iNbrObs_Total);
		gpdYm = DVECTOR(1, gsExpoVars.iNbrObs_Total);
		gpdYd = DVECTOR(1, gsExpoVars.iNbrObs_Total);
		gpdXi = DVECTOR(1, gsExpoVars.iNbrObs_Total);
		gpiNi = IVECTOR(1, gsExpoVars.iNbrObs_Total);
		gpdYsum = DVECTOR(1, gsExpoVars.iNbrObs_Total);
		
		giNmiss = READ_OBSDATA2V(gsExpoVars.iNbrObs_Total, gpdxxi, gpdyyi);

		gsExpoVars.iNbrObs_Total -= giNmiss;  /* extern variable giNbr_Obs has been changed. */

		for (i=1; i<=gsExpoVars.iNbrObs_Total; i++) {
			gpdXi[i] = gpdYsum[i] = gpdYm[i] = gpdYd[i] = 0;
			gpiNi[i] = 0;
		}
		Sort_2_By_Dose(gsExpoVars.iNbrObs_Total, gpdxxi, gpdyyi);  /*Sort gpdxxi and make appropriate changes to gpdyyi */
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
			for (j=1; j<=gpiNi[i]; j++){
				gpdYd[i] += (gpdyyi[jj]-gpdYm[i])*(gpdyyi[jj]-gpdYm[i])/(gpiNi[i]-1);
				jj += 1;
			}
		}
	} // end else (gsExpoVars.iIn_Type==1)
	gpdScXi = DVECTOR(1, gsExpoVars.iNbrObs_Total);
	gpdScYd = DVECTOR(1, gsExpoVars.iNbrObs_Total);
	gpdScYm = DVECTOR(1, gsExpoVars.iNbrObs_Total); 
	gdmaxYd = gpdYd[1]; 
	gdmaxYm = gpdYm[1];

	for (i=2; i<=gsExpoVars.iNbrObs_Total; i++){
		if(gpdYd[i] > gdmaxYd)
			gdmaxYd = gpdYd[i];
		if(gpdYm[i] > gdmaxYm)
			gdmaxYm = gpdYm[i];
	}

	for (i=1; i<=gsExpoVars.iNbrObs_Total; i++){
		gpdScYd[i] = gpdYd[i]/gdmaxYd;
		gpdScYm[i] = gpdYm[i]/gdmaxYm;
	}

	if (gsExpoVars.iNbrObs_Total < 3-iExponential_Known)
		ERRORPRT("Observation # < parameter # for Exponential model.");

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
	
	// Added 1/28/2008, with Bruce
	for (i=1;i<=gsExpoVars.iNbrObs_Total;i++)
	{
		gpdScXi[i] = gpdXi[i]/gdxmax;
	}

	if(gsExpoVars.iCons_Var==1) 
		gpiSpecVector[2]=1; 

	/* Print model and file information on output page */
	//HeaderAndReadVar_2Out(argv[1], ctime(&ltime));

	//Compute for likelihoods
	dlikeA1=dlikeA2=dlikeA3=dvv=dlikeR=0.0;

	int iNtot=0;
	double dtemp=0.0;

	// Likelihood A1: Yij = Mu(i) + e(ij), Var{e(ij)} = Sigma^2
	for (i=1; i<=gsExpoVars.iNbrObs_Total; i++){
		iNtot += gpiNi[i];
	}

	for(i = 1; i <= gsExpoVars.iNbrObs_Total; i++)
		dvv += ((gpiNi[i] - 1)*gpdYd[i]);

	dvv = dvv/iNtot;

	/* Compute Likelihood for model A1: Yij = Mu(i) + e(ij)
	Var(e(ij)) = Sigma^2 */

	dlikeA1 = - log(dvv)*iNtot/2 - iNtot/2;


	/* Compute Likelihood for model A2: Yij = Mu(i) + e(ij)
	Var(e(ij)) = Sigma(i)^2*/
	for(i = 1; i <= gsExpoVars.iNbrObs_Total; i++) {
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
	//#ifdef DO_LOG
	//if (giDo_Log)	// Print values to log for investigation
	//{
	//	fprintf(fp_log,"\n\niCons_Var before AThree_fit\n");
	//	fprintf(fp_log,"   iCons_Var = %d\n", gsExpoVars.iCons_Var);
	//}
	//#endif
	if (gsExpoVars.iCons_Var == 0) { // Parameters for fitting the model above
		//#ifdef DO_LOG
		//if (giDo_Log)	// Print values to log for investigation
		//{
		//	fprintf(fp_log,"\n\nList of Data before AThree_fit\n");
		//	for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
		//	{
		//		fprintf(fp_log,"     %5.4g%7d%14.3g%14.3g\n", gpdXi[i], gpiNi[i], gpdYm[i], pdStDev[i]);
		//	}
		//}
		//#endif
		double *pdLKParms;

		pdLKParms = (double *) malloc((size_t) giNbrParm_A3*sizeof(double));
		for(i=0; i < giNbrParm_A3; i++)
			pdLKParms[i]=0.0;

		/******* Fit the A3 model above *********/
		//AThree_fit(giNbrParm_A3, pdLKParms, EPS, &gpiIter, &dlikeA3);
		//#ifdef DO_LOG
		//if (giDo_Log)	// Print values to log for investigation
		//{
		//	for (i=0; i<giNbrParm_A3; i++) {
		//		fprintf(fp_log, "LKParms[%d] = %.8f\n", i, pdLKParms[i]);
		//	}
		//	fprintf(fp_log,"\n\nList of Data after AThree_fit\n");
		//	for(i=1; i <= gsExpoVars.iNbrObs_Total; i++)
		//	{
		//		fprintf(fp_log,"     %5.4g%7d%14.3g%14.3g\n", gpdXi[i], gpiNi[i], gpdYm[i], pdStDev[i]);
		//	}
		//}
		//#endif

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
	dYbar = gpdYm[1] * gpiNi[1];

	for (i = 2; i <= gsExpoVars.iNbrObs_Total; i++)
		dYbar += gpdYm[i] * gpiNi[i];

	dYbar = dYbar / iNtot;

	dvv = gpdYd[1] * (gpiNi[1] - 1) + gpiNi[1]*(gpdYm[1] - dYbar)*(gpdYm[1] - dYbar);

	for (i = 2; i <= gsExpoVars.iNbrObs_Total; i++)
		dvv += gpdYd[i] * (gpiNi[i] - 1) + gpiNi[i]*(gpdYm[i] - dYbar)*(gpdYm[i] - dYbar);

	dvv = dvv / iNtot;

	dlikeR = -iNtot * (1.0 + log(dvv))/2.0;

	//Allocate memories for sMLEs
	pdLogLikelihood = DVECTOR(1, 4);
	pdModel2 = DVECTOR(1, NBR_OF_PARMS);
	pdModel3 = DVECTOR(1, NBR_OF_PARMS);
	pdModel4 = DVECTOR(1, NBR_OF_PARMS);
	pdModel5 = DVECTOR(1, NBR_OF_PARMS);

	if (!pdLogLikelihood || !pdModel2 || !pdModel3 || !pdModel4 || !pdModel5) 
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

	FREE_DVECTOR (pdLogLikelihood, 1, 4);
	FREE_DVECTOR (pdModel2, 1, NBR_OF_PARMS);
	FREE_DVECTOR (pdModel3, 1, NBR_OF_PARMS);
	FREE_DVECTOR (pdModel4, 1, NBR_OF_PARMS);
	FREE_DVECTOR (pdModel5, 1, NBR_OF_PARMS);

	FREE_DVECTOR (gpdIniParaVector, 1, NBR_OF_PARMS);
	FREE_IVECTOR (gpiSpecVector, 1, NBR_OF_PARMS);
	FREE_IVECTOR (gpiIniSpVector, 1, NBR_OF_PARMS);
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
	free(pdInitParms);
	free(pdLogLikelihood);
	free(pdModel2);
	free(pdModel3);
	free(pdModel4);
	free(pdModel5);

	if(gpiSpecVector[1] == 1 && pdParms[1] == 0)
		FREE_ALVECTOR(psAnasum, 1, 4);
	else
		FREE_ALVECTOR(psAnasum, 1, 5);

	free(pdParms);

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
	for(i=1;i<=iNTotal;i++){
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
