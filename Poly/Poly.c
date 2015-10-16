/*******************************************************************
*
* IMPORTANT NOTE:  The following variable is the version number for
*                  the current model.  THIS MUST BE CHANGED as
*                  important changes are made to the models.
*
********************************************************************/
char Version_no[]="Polynomial Model. (Version: 2.20;  Date: 10/22/2014)";
/*
char Version_no[]="Polynomial Model. (Version: 2.19;  Date: 06/25/2014)";
char Version_no[]="Polynomial Model. (Version: 2.18;  Date: 05/19/2014)";
char Version_no[]="Polynomial Model. (Version: 2.17;  Date: 01/28/2013)";
char Version_no[]="Polynomial Model. (Version: 2.16;  Date: 05/26/2010)";
*/
/*******************************************************************
*
* Poly.C - a ANSI C program for Polynomial model fitting
*          with/without a natural background rate in Benchmark Dose.
*
* Date: Oct, 2000
*
********************************************************************
* Modification Log:
*
* Version Number: 2.3
* Modified By: Q. He
* Modified Date: 9/03/2003
* Reason:
*
* Version Number: 2.4
* Modified By: Geoffrey Nonato
* Modified Date: 3/18/2005
* Reason: Fix problem with specified data not being used properly
*         and degree of freedom.  Changed/beautify the main() and
*         Poly_fit() function to a coding standard.  Changed the
*         variable "silent" to bNo_Log.
*
* Version Number: 2.5
* Modified By: Micheal Ferree
* Modified Date: 6/21/2005
* Reason: Fixed degrees of freedom problems.  Changed output to
*		  display all Tests of interest with note about Test 2 and 3
*         being the same if rho=0.  Changed confidence level from
*         0.05 to 0.1 for tests 2, 3, and 4.  Will try to clean code
*         for version 2.5.1.
* 
* Version Number: 2.6
* Modified By: R. Woodrow Setzer
* Modified Date: 10/14/2005
* Reason: Fixed invalid memory writes, free all allocated memory, 
*         changed access to bind (bind is parallel to X, and has
*         the same length); changed all 'float' declarations to 'double'
*         (this had been done with '#define float double', which seems
*         dangerous to me).
*         - Turn on logging with a compile switch (LOGGING_ON)
*         - Change slog to Slog, and delete local definition of slog
*         - replace calls to local Binary_root with calls to zeroin
*
* Version Number: 2.7
* Modified By: R. Woodrow Setzer
* Modified Date: 12/06/2005
* Reason: Always allocate 5 elements to anasum
*
* Version Number: 2.8
* Modified By: R. Woodrow Setzer
* Modified Date: 03/21/2006
* Reason: moved calculations for likelihoods for models A1, A2, and
*         R to calculate_continuous_liks(), which fixes errors
*         in calculation of likelihood for R.
*
* Version Number: 2.9
* Modified By: R. Woodrow Setzer
* Modified Date: 03/23/2006
* Reason: Allow likelihood for A3 to used user specified variance
*         parameter values.
*         Fixed (?) parameter standard error calculations.
*
* Version Number: 2.10
* Modified By: R. Woodrow Setzer
* Modified Date: 09/8/2006
* Reason: Change variance model, and, in BMDL calculation,
*         change likelihood constraint from equality to inequality
*
* Version Number: 2.11
* Modified By: Geoffrey
* Date: 1/12/2007
* Reason: Incremented version number.
*		 Added last parameter "1" (print SE) in OP_ParmsE().
*
* Version Number: 2.12
* Modified By: Woodrow Setzer
* Date: 2/20/2007
* Reason: Incremented version number to reflect changed compilation options.
*
* Version Number: 2.13
* Modified By: G. Nonato
* Modification Date: 04/08/2008
* Reason: (Per BMDS 2.0: Problem Report 157 & 147)
*       Fix the Observation # < parameter # for Polynomial model problem.
*       Added code to free-up allocated memories before exiting thru ERRORPRT()
*
* Version Number: 2.14
* Modified By: Geoffrey
* Date: 12/09/2008
* Reason: Modified the BMDL_func, by putting the call to getcl_() into a loop
*		 and added printing for logging information to trace/investigate 
*		 the data being passed to the fortran code.  The modification includes printing of
*		 .002 data as well.
*
*
* Version Number: 2.15
* Modified By: G. Nonato
* Modification Date: 10/28/2009
* Reason:
*      To be able to process files/folders with spaces (PR 257)
*      Fix program freeze due to long variable names (PR 278)
*      Process long files up to 256 characters (PR 303 and 308)
*      Modify code for easy maintenance, all lengths for file name,
*        model name, and column names are placed in benchmark.h
*
* Version Number: 2.16
* Modified By: G. Nonato
* Modification Date: 05/26/2010
* Reason: PR 319
*      Change ERRORPRT("Observation # < parameter # for Multistage Cancer model.");
*      to Warning("Observation # < parameter # for Multistage Cancer model.");
*
* Version Number: 2.17
* Modified by: Louis Olszyk
* Modification Date: 01/28/2013
* Reason: - PR461 - replace "Absolute risk", "Relative risk" and
*           "Point risk" in output with "Absolute deviation",
*           "Relative deviation" and "Point estimate", respectively.
*         - PR444 - Improve plot titles
*
* Version Number: 2.18
* Modified by: Louis Olszyk
* Modification Date: 05/19/2014
* Reason: - PR506 - Corrected Fortran function declarations
*
* Version Number: 2.19
* Modified by: Louis Olszyk
* Modification Date: 06/25/2014
* Reason: - PR502 - Fixed gunit() index in PolyCompIneq() fortran routine,
*           which produced bad MLEs for certain datasets with non-constant
*           variance.
*
* Version Number: 2.20
* Modified by: Cody Simmons
* Modification Date: 10/22/2014
* Reason: - PR487 - allow . in path/filename
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

extern void getmle_(int *nparm, double parms[], double parms2[],double *ll,
					int *optite, int *nresm, int bind[], int *flag);


extern void getcl_(int *which, int *ndoses, double doses[], double means[],
				   int nanimals[], double svar[], int *nparm, double *bmr,
				   double *bmd, double *target, double parms[],
				   int fixed[], double fixedval[], int *risktype,
				   int *restrict, double *bmdl, double parms2[],
				   int *optite, int *nresm, int bind[], int *adv,
				   int *model_type, int *flag, double mleparms[]);

extern void getmlea3_(int *ndoses, double doses[], double means[],
					  int nanimals[], double svar[], int *nparm,
					  double parms[],
					  int fixed[], double fixedval[],
					  int *restrict, double parms2[],double *ll,
					  int *optite, int *nresm, int bind[]);


extern void loadcommbloc_(int *ndoses, double doses[], double means[],
						  int nanimals[], double svar[], int *nparm,
						  int fixed[], double fixedval[], int *restrict,
						  int *adverse, int *model_type, double *xmax,
						  double *xmin);

#define EPS 3.0e-8
#define EPSS 3.0e-7
#define FPMIN 1.0e-30
#define MAXIT 100
#define TOLX (10*EPS)
#define STPMX1 10.0
#define ALF 0.000001
#define float double
#define GOLD 1.618034
#define GLIMIT 100
#define TINY 1.0e-20
#define SWAP(a,b, junk) (junk)=(a); (a)=(b); (b)=(junk);
#define SHFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d);
#define ZEPS 1.0e-8
#define MOV3(a,b,c, d,e,f)(a)=(d);(b)=(e);(c)=(f);
#define TOL 2.0e-6

//#define LOGGING_ON 1

int fixedParm(int);
void GetNewParms(double *p, int size);
void GetMoreParms(double *p, int size);
void GetOtherParms(double *p, int size);
double BMDL_func(int nparm, double xlk, double Dose, double pBak[], double gtol);
void PolyMeans(int nobs, double p[], double Doses[], double means[]);
void Get_BMRS(double *p, double Dose, double bmdl1, double *BMRVals, int sign, int
			  bmr_type);
void OUTPUT_BENCHMD2(int pdcol, double BMD);
void Var2Part(int obs, int const_var, double Vi, double meani, double *p,
			  double *mg, double **mg2, double **vg2);
void MeanPart(int obs, double *p, double *mg);
void VarPart(int obs, int const_var, double Vi, double meani, double *p,
			 double *mg, double *vg);
void Mean2Part(int obs, double *p, double **mg2);
int Model_DF(int []);

/*** Define input and output files's name  *********************/
char     fin[FLENGTH];  /*input temp file*/
char    fout[FLENGTH];  /*output temp file*/
char    fout2[FLENGTH];
char  *Parm_name[]={"alpha", "rho", "beta_0", "beta_1", "beta_2", "beta_3", "beta_4",
"beta_5", "beta_6", "beta_7", "beta_8", "beta_9", "beta_10",
"beta_11", "beta_12", "beta_13", "beta_14", "beta_15",
"beta_16", "beta_17", "beta_18", "beta_19", "beta_20"};
char *alt_Parm1_name = "lalpha";
char  *anatxt[]={"Full model", "Fitted model", "Reduced model"};
char  plotfilename[FLENGTH]; /* name of plot file */
char fname2[FLENGTH], logfile[FLENGTH], *dot2;

/*  bNo_Log = true -> no log file
bNo_Log = false -> log file is made */
int bNo_Log = true;     /*  switch for log file   */

/*** variables will not be changed except Spec,xxi,and yyi will be
sorted together *******/
int    *Spec, *SpecBkp;    /*vector used to identify user input parm.*/
int    *IniSp;
int     *Ni;       /*number of animals in the i-th dose group*/
double *Ym;      /*mean response data array*/
double *Yd;      /*sample variance (s_i^2) of response data array*/
double *Xi;      /*independent variable (dose value) data array*/
double *xxi;      /*dose values when not divided to dose groups*/
double *yyi;     /*response values when not divided to dose groups*/
double *Ysum;    /*sum of responses within a group*/
double *Rlevel;
double *Bmdl;
double *IniP;
int   in_type, Nobs, ntotal, nparm, restrict, ndeg, initial, appendix, smooth;
int   bmdlCurve, sign;
double xmax, xmin, scale;
double ck;
double PE_bmdl;
int ProfFlag, PEFlag, EmptyFlag;


/** changing variables **/
int   replace, cons_var, brat, bmr_type;
double tD, BMD_lk, LR, upb=18, BMR, react;
int ErrorFlag;                   /* Error States from DMNGB */

/* GLN - 03/07/2008
*  Free-up allocated memory before exit
*  upon encountering fatal error.
*/
void FreeUp_mem(double *Parms, double *Parmscopy, VarList *varsum, AnaList *anasum, double  **vcv, double **vcv_adj, int var_type, int adj_vcv_rows)
{
	if(Parmscopy)
		FREE_DVECTOR (Parmscopy, 1, nparm);
	FREE_IVECTOR (Ni,1,Nobs);
	FREE_DVECTOR (Xi,1,Nobs);
	FREE_DVECTOR (Ym,1,Nobs);
	FREE_DVECTOR (Yd,1,Nobs);
	FREE_DVECTOR (IniP,1,nparm);
	FREE_IVECTOR (Spec,1,nparm);
	FREE_IVECTOR (IniSp,1,nparm);
	FREE_VLVECTOR(varsum,1,3);

	if(var_type == 0)
		FREE_ALVECTOR(anasum,1,4);
	else
		FREE_ALVECTOR(anasum,1,5);

	FREE_DVECTOR (Rlevel,1,5);
	FREE_DVECTOR (Bmdl, 1, 5);
	FREE_DMATRIX (vcv,1,nparm,1,nparm);

	if(adj_vcv_rows > 0)
		FREE_DMATRIX (vcv_adj, 1, adj_vcv_rows, 1, adj_vcv_rows);

	FREE_DVECTOR (Parms, 1, nparm);

	if (fp_log != (FILE *) NULL)
		fclose(fp_log);

	return;
}

/****************************************************************
** main--main function used to call Poly mode fitting program.
Includes: biosubcc.c--common subfunction C program.
*****************************************************************/
int main (int argc, char *argv[])
{
	void Poly_fit(int nparm, double p[], int *is_conv,int *iter, double *fret, int *bounded);
	void Poly_BMD (int nparm, double p[], double gtol, int *iter, double xlk,
		double Rlevel[], double Bmdl[],double *BMD, int *bmderror);
	void Poly_vcv(int nparm, int Spec[], double p[], double **vcv);
	int READ_OBSDATA4V(int Nobs,double Xi[],int Ni[],double Ym[],double Yd[]);
	int READ_OBSDATA2V(int ntotal, double xxi[], double yyi[]);
	void OUTPUT_DTMSnVCV(int nparm, int nboundparm, int Spec[], char *parmtxt[],double **vcv,int *bounded);
	double betai(double a, double b, double x);

	void AThree_Fit(int nparm, double p[], double gtol, int *iter, double *fret);
	double BMD_func(int nparm, double p[], double x, double gtol);
	double BMD1_func(int nparm, double p[], double x, double gtol);


	int    iter,i, j, jj, junk;          /*iteration variable*/
	int    group;                            /*temp counter eventually equals Nobs*/
	int    conv_check;                       /* = -1 if the optimizer fails */
	int    bmdflag;  /*flag for computing profile likelihood */
	int    bmdose;        /*flag for computing benchmark dose*/
	int    Nmiss;         /*number of records with missing values*/
	int    nparm_known, poly_known;   /*number of specified parameters */
	double  xlk, lkA3, lkA1, lkA2, lkR;   /*log likelihoods */
	double  BMD, lep, upep, frac1;
	double *stdev;         /*Sample Standard deviation*/
	int linflag = 0;       /* ==1 if model is linear, ==0 if polynomial
						   of degree > 1.  Since linear and poly-
						   nomial models are both given as model
						   options, they should bevoid Mean2Part(int obs,
						   double *p, double **mg2) distinguished
						   in output pages*/

	double  *Parms, *LKParms, *Parmscopy;   /*parameter array, parameter array for likelihood
											test*/
	VarList *varsum;      /*info for variables--p. dep.,n. dep., indep.*/
	AnaList *anasum;      /*information for ANONA analysis*/
	double  **vcv;        /*variance and covariance matrix*/
	double **vcv_adj;        /* The vcv matrix with fixed parameters and
							 parameters estimated at a boundary point
							 removed */
	int *bounded;            /* 1 if a parameter is estimated at a boundary
							 0 otherwise */
	int adj_vcv_rows;        /* number of rows/columns in vcv matrix after
							 fixed parameters or parameters estimated
							 at a bound are removed from it */

	char   model_name[MNLENGTH], user_note[UNLENGTH];
	char  dose_name[CNLENGTH], no_name[CNLENGTH], mean_name[CNLENGTH], stdev_name[CNLENGTH], junkname[FLENGTH];
	char  response_name[CNLENGTH];
	char long_path_name[FLENGTH];

	int   tmpi1, tmpi2, Ntot, frac2, var_type;
	double  bmr_root=0, Rel_Conv, Parm_Conv;
	double sp, Nd;
	double *doses, *means, *svar, *parms;  /* These are for loading in the common block */
	int *nanim, *Spec3;
	int adverse;      /* adverse direction flag to go into the common block */
	int model_type;   /* model_type = 0 for polynomial */
	int NobsL, nparmL, restrictL; /* These are just copies to send to FORTRAN */
	/* in int format                       */
	double *mean, *std;

	int temp_sign;

	time_t ltime;
	/* Set time zone from TZ environment variable. If TZ is not set,
	* the operating system is queried to obtain the default value
	* for the variable.
	*/
	/* _tzset();*/

	time( &ltime );

	//Min_increment=1.0e-30;
	//Max_double= 1.0e30;

	/********************************************************************
	* {QH 2004/01/14 PR# }
	* Added to show version number if executed from command line with -v
	*********************************************************************/
	if(argc == 2)
		show_version(argv[1], Version_no);

	if(argc < 2) {
		fprintf(stderr, "ERROR:  Requires two arguments\nUsage:  %s <file.(d)>\n", argv[0]);
		fprintf (stderr, "   or:  %s -v for version number.\n", argv[0]);
		exit (1);
	} /* end if */

	/* commented out for PR487 - does not allow "." in pathname */
	/*open bmdswrk.001 input and Poly.out output temp files*/
	/*for(i = 0; i < FLENGTH; i++){
		if(argv[1][i] == '.') {
			if(argv[1][i+1] != '(' && argv[1][i+2] != 'd' && argv[1][i+3] != ')'){
				argv[1][i+1] = '(';
				argv[1][i+2] = 'd';
				argv[1][i+3] = ')';
				break;
			}
			else
				break;
		}
	} */

	/*if (argc >= 2)
	path_name(argv[1]);*/

	/* allow spaces in the path name */
	if (argc > 2){
		path_name2(argc, argv, long_path_name);
		argv[1] = long_path_name;
	}

	fp_in=fopen(argv[1], "r");

	if (fp_in==NULL) {
		/*printf("Error in opening input  file.\n");
		printf ("...now exiting to system...\n");*/

		fprintf(stderr, "Error in opening input file.\n");
		fprintf (stderr, "...Exited to system!\n");
		exit (1);
	}

	/* open the log file if bNo_Log = 0 */
#ifdef LOGGING_ON
	{
		strcpy(logfile,argv[1]);
		dot2 = strchr(logfile, (int) '.');
		(*dot2) = (char) 0;
		strcpy(fname2,logfile);
		strcat(logfile,"-poly.log");
		fp_log = fopen(logfile, "w");

		if (fp_log == (FILE *) NULL)
			ERRORPRT("Unable to open log for Poly.C.");
	}
#endif

	fscanf(fp_in, "%s", model_name);
	fscanf(fp_in, "%[ ^\n]", user_note);
	fscanf(fp_in, "%[^\n]", user_note);
	/*fscanf(fp_in, "%s", input_name);
	fscanf(fp_in, "%s", output_name);*/
	fscanf(fp_in, "%s", junkname);
	fscanf(fp_in, "%s", junkname);

	fscanf(fp_in, "%d",&ndeg);
	fscanf(fp_in, "%d",&in_type);
	/* in_type=1 if input format is Xi, Ni, Ym, Yd. */
	/* in_type=0 if input format is Ntotal, Xi, Y_ij. */

	if (in_type==1)
		fscanf(fp_in, "%d", &Nobs);
	else
		fscanf(fp_in, "%d", &ntotal);

	fscanf(fp_in, "%d", &sign);

	if(ndeg == 1)
		linflag = 1;
	else
		linflag = 0;

	/*assign number of parameters*/
	nparm = ndeg+1+2;

#ifdef LOGGING_ON
	{
		fprintf(fp_log,"model_name: %s\n\n",model_name);
		fprintf(fp_log,"INPUT VALUES FOR DATA SET\n\n");
		fprintf(fp_log,"ndeg = %d           (Polynomial Degree)\n", ndeg);
		fprintf(fp_log,"in_type = %d        (Input Format  in_type=1 if input format is Xi, Ni, Ym, Yd\n",in_type);
		fprintf(fp_log,"                                  in_type=0 if input format is Ntotal, Xi, Y_ij)\n");
		if (in_type==1)
			fprintf(fp_log,"Nobs = %d           (Number of Observations)\n",Nobs);
		else
			fprintf(fp_log,"ntotal = %d         (Total Number of Subjects)\n",ntotal);
		fprintf(fp_log,"sign = %2d          (Adverse Direction; 0=Automatic, 1=Up, -1=Down)\n",sign);
		fprintf(fp_log,"linflag = %d        (Flag if the Degree is one)\n",linflag);
		fprintf(fp_log,"nparm = %d          (The Number of Parameters in the Model)\n",nparm);
		fflush(fp_log);
	}
#endif

	/*allocate memory for arrays*/
	Parms = DVECTOR(1, nparm);
	Spec = IVECTOR(1, nparm);
	IniP = DVECTOR(1, nparm);
	IniSp = IVECTOR(1, nparm);
	varsum = VLVECTOR(1, 3);
	Rlevel = DVECTOR(1, 5);
	Bmdl = DVECTOR(1, 5);

	fscanf(fp_in,"%d%lf%lf%d%d%d%d%d", &ITMAX, &Rel_Conv, &Parm_Conv, &bmdlCurve, &restrict, &bmdose, &appendix, &smooth);
	/* restrict=0 if no restriction, +1 if beta_j>0, -1 if beta_j<0 */

	fscanf(fp_in,"%d%lf%d%lf",&bmr_type,&bmdparm.effect,&cons_var,&bmdparm.level);
	bmdparm.risk = 1;
	junk = 0; /*Used to see if an extension was added to output file name */
	/*Search the output_name array for any extensions, as we don't
	want to have two extensions attached to the file name */

#ifdef LOGGING_ON
	{
		fprintf(fp_log,"ITMAX = %d                 (Maximum number of iterations)\n",ITMAX);
		fprintf(fp_log,"Rel_Conv = %g     (Rel Function Convergence, default=2.22045e-16)\n",Rel_Conv);
		fprintf(fp_log,"Parm_Conv = %g    (Parameter Convergence)\n",Parm_Conv);
		fprintf(fp_log,"bmdlCurve = %d               (BMDL Curve Calculation; 1=yes, 0=no)\n",bmdlCurve);
		fprintf(fp_log,"restrict = %2d               (Restriction on polynomial coefficients)\n",restrict);
		fprintf(fp_log,"bmdose = %d                  (BMD Calculation; 1=yes, 0=no)\n",bmdose);
		fprintf(fp_log,"appendix = %d                (Append or Overwrite output file; 1=append, 0=overwrite)\n",appendix);
		fprintf(fp_log,"smooth = %2d                 (Smooth Option; 0=unique, 1=c-spline)\n",smooth);
		fprintf(fp_log,"bmr_type = %d                (BMR Type; 0=AbsDev, 1=StdDev, 2=RelDev(default), 3=Point, 4=Extra)\n",bmr_type);
		fprintf(fp_log,"bmdparm.effect = %4g       (BMR Factor, default=0.1)\n",bmdparm.effect);
		fprintf(fp_log,"cons_var = %d                (Constant Variance; 0=no, 1=yes)\n",cons_var);
		fprintf(fp_log,"bmdparm.level = %4g        (Confidence level, default=0.95)\n",bmdparm.level);
	}
#endif

	Get_Names(argv[1], fout, fout2, plotfilename);

	if(appendix==Yes)
		fp_out=fopen(fout,"a");
	else
		fp_out=fopen(fout,"w");
#ifndef RBMDS
	fp_out2=fopen(fout2,"w");
#endif
	if (fp_out==NULL 
#ifndef RBMDS
		||fp_out2==NULL 
#endif
		) {
			printf("Error in opening  output files.\n");
			printf ("...now exiting to system...\n");
			fprintf(fp_out,"Error in opening output files.\n");
			fprintf (fp_out,"...Exited to system!\n");
			exit (1);
	}

	/* Print model and file information on output page */
	Output_Header(Version_no, argv[1], plotfilename, ctime(&ltime), user_note);

	/*obtain user input parameters*/
	READ_PARAMETERS(nparm,Parms);
	FILL_SPECVECTOR(nparm,Parms,Spec);
	nparm_known = COUNT_SPECVECTOR(nparm,Spec);
	if (cons_var == 1) {
		Spec[2] = 1;
		Parms[2] = 0.0;
	} else {
		Parm_name[0] = alt_Parm1_name;
	}


	poly_known=0;
	for (i=3; i<=nparm; i++){
		if (Spec[i] == 1)
			poly_known += 1;
	}


	brat=Yes;

	if (Spec[1]==1 && Parms[1]<EPS)
		brat=No;

	/*obtain user input initial parameters values      *******/
	fscanf(fp_in,"%d", &initial);
	READ_PARAMETERS(nparm,IniP);
	FILL_SPECVECTOR(nparm,IniP,IniSp);

	for(i = 1; i <= nparm; i++){
		if(Spec[i] == 1)
			IniP[i] = 1;
	}

	/* var_type should always be 1 if cons_var is 0, and 0 if cons_var is 1 */
	/* GLN - 03/18/2005, new code base on the commented code below */
	if(cons_var==0) // not checked
		var_type = 1;
	else
		var_type = 0;
	/* - end of modification - GLN - 03/18/2005*/

	/* commented by GLN - 03/18/2005 to reflect the users selection of constant variance in GUI see modified code above
	var_type = 1;
	cons_var = 0;
	- end - GLN - 03/18/2005 */

	if(Spec[2] == 1){
		if(Parms[2] == 0){
			anasum = ALVECTOR(1, 5);
			var_type = 0;
			cons_var = 1;
		}
		else
			anasum = ALVECTOR(1, 5);
	}
	else
		anasum = ALVECTOR(1, 5);

#ifdef LOGGING_ON
	{
		fprintf(fp_log,"\nThe specified or default parameters are:\n");
		fprintf(fp_log,"                                           Spec Value\n");
		fprintf(fp_log,"Parms[1] = %12.6g        (alpha)        %d\n",Parms[1],Spec[1]);
		fprintf(fp_log,"Parms[2] = %12.6g        (rho)          %d\n",Parms[2],Spec[2]);

		for (i=3; i<=nparm; i++)
			fprintf(fp_log,"Parms[%d] = %12.6g        (beta%d)        %d\n",i,Parms[i],i-3,Spec[i]);

		fprintf(fp_log,"\ninitial = %d                (Number of initialized parameters)\n",initial);
		fprintf(fp_log,"The initial parameter values are:\n");
		fprintf(fp_log,"                                           IniSp Value\n");
		fprintf(fp_log,"IniP[1] = %12.6g         (alpha)         %d\n",IniP[1],IniSp[1]);
		fprintf(fp_log,"IniP[2] = %12.6g         (rho)           %d\n",IniP[2],IniSp[2]);

		for (i=3; i<=nparm; i++)
			fprintf(fp_log,"IniP[%d] = %12.6g         (beta%d)         %d\n",i,IniP[i],i-3,IniSp[i]);

		fprintf(fp_log,"The Spec and IniSp vectors are just flags to tell if the user gave\n");
		fprintf(fp_log,"a value for the specific parameter.\n");
		fprintf(fp_log,"\nvar_type = %d     (for some reason this value should be opposite cons_var)\n",var_type);
	}
#endif

	/*obtain observation data into Yp, Yn, Xi, Ls, Xg vectors*/
	if (in_type == 1)
		fscanf(fp_in,"%s%s%s%s", dose_name, no_name, mean_name, stdev_name);
	else
		fscanf(fp_in,"%s%s", dose_name, response_name);

	if (in_type==1) {
		Ym = DVECTOR(1, Nobs);
		Yd = DVECTOR(1, Nobs);
		Xi = DVECTOR(1, Nobs);
		Ni = IVECTOR(1, Nobs);

		/*The following three lines of code were added because
		* the user is asked to enter and standard deviation, but the
		* program treats entry as a variance */

		stdev = DVECTOR(1, Nobs);

		Nmiss = READ_OBSDATA4V(Nobs, Xi, Ni, Ym, stdev);

		for(i = 1; i <= Nobs; i++)
			Yd[i] = stdev[i]*stdev[i];

		Nobs -= Nmiss;   /* extern variable Nobs has been changed. */

		FREE_DVECTOR(stdev,1,Nobs);

	}
	else {
		xxi = DVECTOR(1, ntotal);
		yyi = DVECTOR(1, ntotal);
		Ym = DVECTOR(1, ntotal);
		Yd = DVECTOR(1, ntotal);
		Xi = DVECTOR(1, ntotal);
		Ni = IVECTOR(1, ntotal);
		Ysum = DVECTOR(1, ntotal);
		Nmiss = READ_OBSDATA2V(ntotal, xxi, yyi);
		ntotal -= Nmiss; /* extern variable Nobs has been changed. */

		for (i=1; i<=ntotal; i++){
			Xi[i] = Ysum[i] = Ym[i] = Yd[i] = 0;
			Ni[i] = 0;
		}

		Sort_2_By_Dose(ntotal, xxi, yyi);  /* Sort arrays so that the
										   following code
										   works.  Located
										   specialfun.h */
		group = 1;

		for (i=1; i<=ntotal; i++){
			Xi[group]=xxi[i];
			Ni[group] += 1;
			Ysum[group] += yyi[i];

			if(i < ntotal){
				if(xxi[i] != xxi[i + 1])
					group +=1;
			}
		} //end of for

		Nobs=group;
		jj=1;
		for (i=1; i<=Nobs; i++) {
			Ym[i] = Ysum[i]/Ni[i];
			if (Ni[i] > 1) {
				for (j=1; j<=Ni[i]; j++) {
					Yd[i] += (yyi[jj]-Ym[i])*(yyi[jj]-Ym[i])/(Ni[i]-1);
					jj += 1;
				}
			}
			else
				Yd[i]=0.0;
		} //end of for

		FREE_DVECTOR (xxi,1,ntotal);
		FREE_DVECTOR (yyi,1,ntotal);
		FREE_DVECTOR (Ysum, 1, ntotal);
	}

#ifdef LOGGING_ON
	{
		fprintf(fp_log,"\nDOSE's          N's                MEAN's             VAR's\n");

		for (i=1; i<=Nobs; i++)
			fprintf(fp_log,"Xi[%d]=%4g      Ni[%d]=%4d         Ym[%d]=%4g         Yd[%d]=%4g\n",

			i,Xi[i],i,Ni[i],i,Ym[i],i,Yd[i]);
		fprintf(fp_log,"\nNmiss = %2d             (Number of missing values)\n",Nmiss);
		fprintf(fp_log,"Nobs = %2d              (Number of observations)\n",Nobs);
		fflush(fp_log);
	}
#endif

	if (Nobs < ndeg+1-poly_known && restrict == 0)
	{
		//FreeUp_mem(Parms, Parmscopy, varsum, anasum, vcv, vcv_adj, var_type, adj_vcv_rows);
		if(linflag)
		{
			// Commented the following line of code for PR 319
			//ERRORPRT("Observation # < parameter # for Linear model.");
			Warning("Observation # < parameter # for Linear model.");
		}
		else
		{
			// Commented the following line of code for PR 319
			//ERRORPRT("Observation # < parameter # for Polynomial model.");
			Warning("Observation # < parameter # for Polynomial model.");
		}
	}

	if (in_type == 1){
		for (i=1; i<=Nobs; i++){
			if ( (Xi[i] < 0) || (Ni[i] < 0) || (Yd[i] < 0) )
			{
				FreeUp_mem(Parms, Parmscopy, varsum, anasum, vcv, vcv_adj, var_type, adj_vcv_rows);
				ERRORPRT("Values of dose, group size, and sample variance should be positive...");
			}
		}
	}
	else {
		for (i=1; i<=Nobs; i++){
			if (Xi[i] < 0)
			{
				FreeUp_mem(Parms, Parmscopy, varsum, anasum, vcv, vcv_adj, var_type, adj_vcv_rows);
				ERRORPRT("Dose value should be positive ...");
			}
		}
	}
	/*********** end of input data **********/

	xmin = Max_double;
	xmax = 0.0;
	frac1=0.0;
	frac2=0;

	for (i=1;i<=Nobs;i++){
		if (Xi[i] < xmin){
			xmin = Xi[i];
			bmr_root = sqrt(Yd[i]);
		}

		if (Xi[i] > xmax)
			xmax = Xi[i];

		frac1 += (Ni[i]-1)*Yd[i];
		frac2 += Ni[i]-1;
	}
#ifdef LOGGING_ON
	{
		fprintf(fp_log,"xmin = %4g            (minimum dose level)\n",xmin);
		fprintf(fp_log,"xmax = %4g            (maximum dose level)\n",xmax);
		fprintf(fp_log,"bmr_root = %4g       (sqrt of std for min dose level, not sure purpose)\n",bmr_root);
		fprintf(fp_log,"frac1 = %4g        (no idea)\n",frac1);
		fprintf(fp_log,"frac2 = %4d           (sum of Ni's minus Nobs)\n",frac2);
		fflush(fp_log);
	}
#endif

	/*output title and summary of intput data  ****************************/
	OUTPUT_TEXT("\n   The form of the response function is: ");

	OUTPUT_TEXT("\n   Y[dose] = beta_0 + beta_1*dose + beta_2*dose^2 + ...");

	if (in_type == 1)
		fprintf(fp_out,"\n\n   Dependent variable = %s", mean_name);
	else
		fprintf(fp_out,"\n\n   Dependent variable = %s", response_name);

	fprintf(fp_out,"\n   Independent variable = %s", dose_name);

	for (i=1; i<=nparm; i++){
		if (Spec[i] == 1)
			fprintf(fp_out,"\n   %s is set to %g", Parm_name[i-1],Parms[i]);
	}

	if (restrict==1)
		fprintf(fp_out,"\n   The polynomial coefficients are restricted to be positive");
	else
		if (restrict==-1)
			fprintf(fp_out,"\n   The polynomial coefficients are restricted to be negative");
		else
			if (restrict==0)
				fprintf(fp_out,"\n   Signs of the polynomial coefficients are not restricted");

	if(cons_var == 1)
		fprintf(fp_out, "\n   A constant variance model is fit");
	else
		fprintf(fp_out, 
		"\n   The variance is to be modeled as Var(i) = exp(lalpha + log(mean(i)) * rho)");

	nparm_known = COUNT_SPECVECTOR(nparm, Spec);
	fprintf (fp_out, "\n\n   Total number of dose groups = %d",Nobs+Nmiss);
	fprintf (fp_out, "\n   Total number of records with missing values = %d",Nmiss);

	fprintf(fp_out, "\n   Maximum number of iterations = %d\n", ITMAX);
	fprintf(fp_out, "   Relative Function Convergence has been set to: %g\n", Rel_Conv);
	fprintf(fp_out, "   Parameter Convergence has been set to: %g\n\n", Parm_Conv);


	if ( (fabs(Rel_Conv - 1e-8) > 1e-10) || (fabs(Parm_Conv - 1e-8) > 1e-10) ) {
		fprintf(fp_out, "**** We are sorry but Relative Function and Parameter Convergence   ****\n");
		fprintf(fp_out, "**** are currently unavailable in this model.  Please keep checking ****\n");
		fprintf(fp_out, "**** the web site for model updates which will eventually           ****\n");
		fprintf(fp_out, "**** incorporate these convergence criterion.  Default values used. ****\n\n");
	}

	if( bmr_type > 3 || bmr_type < 0){
		fprintf(fp_out, "****  Invalid BMR Type.  Defaulting to std dev risk, value = 1.00\n");
		bmr_type = 1;
		bmdparm.effect = 1.0;
	}

	if(initial==Yes) 
	{
		OUTPUT_TEXT("\n\n                 User Inputs Initial Parameter Values  ");
		OUTPUT_Init(nparm, Spec, IniP, Parm_name);

		for (i=1; i<=nparm; i++)
			if(IniSp[i]==1) {   /* have been initialized. */
				if(Spec[i]==1 )  /* check if it is for fixed parm. */
					Warning("The initial value for the fixed parameter is ignored.");
			}
			else {
				/* check if all the unspecifird parms were initialized. */
				if (Spec[i]==0)
				{
					FreeUp_mem(Parms, Parmscopy, varsum, anasum, vcv, vcv_adj, var_type, adj_vcv_rows);
					ERRORPRT("You have to initialize either ALL or NONE of the unspecified parameters.");
				}
			}

			if (IniP[1] <= 0)
			{
				FreeUp_mem(Parms, Parmscopy, varsum, anasum, vcv, vcv_adj, var_type, adj_vcv_rows);
				ERRORPRT("The initial value of variance has to be positive.");
			}

			tmpi1=tmpi2=0;

			if (restrict==1) 
			{
				for (i=4; i<=nparm; i++) 
				{
					if (IniP[i] <0)
						tmpi1 += 1;
				}
				if (tmpi1 != 0)
				{
					FreeUp_mem(Parms, Parmscopy, varsum, anasum, vcv, vcv_adj, var_type, adj_vcv_rows);
					ERRORPRT("You have restricted all polynomial coefficients to be positive.");
				}
			}
			else
				if (restrict ==-1)
				{
					for (i=4; i<=nparm; i++)
					{
						if (IniP[i] >0)
							tmpi2 += 1;
					}
					if (tmpi2 != 0)
					{
						FreeUp_mem(Parms, Parmscopy, varsum, anasum, vcv, vcv_adj, var_type, adj_vcv_rows);
						ERRORPRT("You have restricted all polynomial coefficients to be negative.");
					}

				}
	}

	/*compute likelihoods for A1, A2, and R*/
	lkA1 = lkA2 = lkA3 = lkR = 0.0;
	compute_continuous_liks(Nobs, Ni, Ym, Yd, &lkA1, &lkA2, &lkR);

	/* Compute Likelihood for model A3: Yij = Mu(i) + e(ij)
	Var(e(ij)) = k*(m(xi))^rho*/
	if(var_type == 1 && (Spec[2] == 0 || (Spec[2] == 1 && Parms[2] != 0.0))) {
		/* Parameters for fitting the model above */
		LKParms = DVECTOR(1, Nobs+2);
		LKParms[1] = Parms[1];
		LKParms[2] = Parms[2];

		/* Fit the A3 model above */
		nparm = Nobs+2;

		AThree_Fit(nparm, LKParms, EPS, &iter, &lkA3);

		nparm = ndeg+3;

		FREE_DVECTOR(LKParms, 1, Nobs+2);
		/* FREE_IVECTOR(Spec, 1, Nobs+2);

		Spec = IVECTOR(1, nparm); */

		/*for(i = 1; i <= nparm; i++)
		Spec[i] = SpecBkp[i];*/ /* Put Spec back together */

		/* FREE_IVECTOR(SpecBkp, 1, nparm); */
	}
	/* If constant variance model then just set equal to model A1. */
	/* MJF, 25MAY2005. */
	else {
		lkA3 = lkA1;
	}

	/* Get initial default adverse direction by finding the largest
	absolute difference between the observed dose means and the
	observed control mean */

	if (sign == 0) {
		sign = Get_Linear_Trend(Nobs, Xi, Ym, Ni);
		temp_sign = 0;
	}
	else
		temp_sign = 9999;

	/* This is all for loading the common block for donlp2 */
	doses = DVECTOR(0, Nobs-1);
	means = DVECTOR(0, Nobs-1);
	svar = DVECTOR(0, Nobs-1);
	nanim = IVECTOR(0, Nobs-1);
	parms = DVECTOR(0, nparm-1);
	Spec3 = IVECTOR(0, nparm-1);

	for(i = 1; i <= Nobs; i++){
		nanim[i-1] = (int) Ni[i];
		doses[i-1] = Xi[i]/xmax;
		means[i-1] = Ym[i];
		svar[i-1] = Yd[i];
	}

	for(i=0; i<nparm; i++) 
	{
		Spec3[i] = Spec[i+1];
		parms[i] = Parms[i+1];
	}

	NobsL = Nobs;
	nparmL = nparm;
	restrictL = restrict;
	adverse = sign;
	model_type = 0;

	loadcommbloc_(&NobsL, doses, means, nanim, svar, &nparmL, Spec3,
		parms, &restrictL, &adverse, &model_type, &xmax, &xmin);

	FREE_DVECTOR(doses, 0, Nobs-1);
	FREE_DVECTOR(means, 0, Nobs-1);
	FREE_DVECTOR(svar, 0, Nobs-1);
	FREE_IVECTOR(nanim, 0, Nobs-1);
	FREE_DVECTOR(parms, 0, nparm-1);
	FREE_IVECTOR(Spec3, 0, nparm-1);
	/* End of load common block */

#ifdef LOGGING_ON
	{
		fprintf(fp_log,"\nlkA1 = %4g          (Likelihood for model A1)\n", lkA1);
		fprintf(fp_log,"lkA2 = %4g          (Likelihood for model A2)\n", lkA2);
		fprintf(fp_log,"lkA3 = %4g          (Likelihood for model A3)\n", lkA3);
		fprintf(fp_log,"lkR = %4g           (Likelihood for model R)\n",lkR);
		fprintf(fp_log,"temp_sign = %4d     (This is 9999 if automatic is not selected and 0 if auto is selected)\n",temp_sign);
		fprintf(fp_log,"sign = %4d          (Sign is assigned here if auto is chosen)\n",sign);
		fflush(fp_log);
	}
#endif

	bounded = IVECTOR(1, nparm);

	/*fitting Polynomial model and output parameter estimators */

#ifdef LOGGING_ON
	{
		fprintf(fp_log, "\n********** Call to Poly_Fit ****************\n");
		fprintf(fp_log, "Variables going in:\n");
		fprintf(fp_log, "nparm = %2d;  conv_check = %2d;  iter = %4d;  xlk = %6g\n",
			nparm, conv_check, iter, xlk);

		for (i=1; i<=nparm; i++) {
			fprintf(fp_log,"Parms[%d] = %13g          bounded[%d] = %2d\n",
				i,Parms[i],i,bounded[i]);
		}
		fflush(fp_log);
	}
#endif

	Poly_fit(nparm, Parms, &conv_check, &iter, &xlk, bounded);

#ifdef LOGGING_ON
	{
		fprintf(fp_log, "Variables coming out:\n");
		fprintf(fp_log, "nparm = %2d;  conv_check = %2d;  iter = %4d;  xlk = %6g\n",
			nparm, conv_check, iter, xlk);
		for (i=1; i<=nparm; i++) {
			fprintf(fp_log,"Parms[%d] = %13g          bounded[%d] = %2d\n",	i,Parms[i],i,bounded[i]);
		}
		fflush(fp_log);
	}
#endif
	/* If a bad completion code is returned (conv_check == -1), then
	exit the program */

	if(conv_check == -1){
#ifndef RBMDS
		fprintf(fp_out2, "END  %d", -1);
#endif
		/* Let plotting program know that it didn't make it    */
#ifdef LOGGING_ON
		{ 
			fflush(fp_log);
			fclose (fp_log);
		}
#endif
		CLOSE_FILES ();
		exit(0);
	}

	/* Compute the aprox. covariance matrix */

	vcv = DMATRIX (1,nparm,1,nparm);

	INITIALIZE_DMATRIX(vcv, nparm, nparm);

	/* Creates matrix of second partial derivatives of the negative
	log likelihood function */
	Poly_vcv(nparm,Spec,Parms,vcv);

	adj_vcv_rows = 0;

	for (i = 1; i <= nparm; i++)
	{
		if (bounded[i] == 0)
		{
			adj_vcv_rows++;
		}
	}				/* end for (i = 1; i <= nparm; i++) */

	vcv_adj = DMATRIX(1, adj_vcv_rows, 1, adj_vcv_rows);


	/*  Make a copy of mle parameters for later use, (getprofile) */
	Parmscopy  = DVECTOR (1, nparm);

	for (i = 0 ; i < nparm ; i ++){
		Parmscopy[i] = Parms[i+1];
	}


	/*output covariance matrix (correlation matrix) */
	Get_and_OUTPUT_DTMSVCV (nparm, Spec, Parm_name, vcv, vcv_adj, bounded);

	/* Output parameter estimates and standard errors */
	OP_ParmsE(nparm,Spec,Parms,Parm_name,vcv_adj, bounded, bmdparm.level, 1);


	/* Calculate fitted means and standard deviations, and print
	them out in a table with the data, for comparison. */
	mean = DVECTOR (1, Nobs);
	std = DVECTOR (1, Nobs);

	for(i = 1; i <= Nobs; i++)
	{
		mean[i] = Parms[nparm];
		for(j = nparm-1; j >= 3; j--)
		{
			mean[i] = Xi[i]*mean[i] + Parms[j];
		}

		if(Parms[2] == 0)
			std[i] = sqrt(Parms[1]);
		else
			std[i] = sqrt(exp(Parms[1] + log(fabs(mean[i])) * Parms[2]));
	}
	PrintData(mean, std, Xi, Ym, Yd, Ni, Nobs);

	FREE_DVECTOR(mean, 1, Nobs);
	FREE_DVECTOR (std, 1, Nobs);

#ifdef LOGGING_ON
	{
		fprintf(fp_log,"*************** Call to DTMS3ANOVAC ****************\n");
		fprintf(fp_log,"Variables going in:\n");
		fprintf(fp_log,"nparm = %4d;  Nobs = %4d;  var_type = %3d;\n",nparm, Nobs, var_type);
		fprintf(fp_log,"lkA1 = %13g\n",lkA1);
		fprintf(fp_log,"lkA2 = %13g\n",lkA2);
		fprintf(fp_log,"lkA3 = %13g\n",lkA3);
		fprintf(fp_log,"lkR = %13g\n",lkR);
		fprintf(fp_log,"xlk = %13g          (fitted likelihood)\n",xlk);

		for (i=1; i<=nparm; i++) {
			fprintf(fp_log,"Spec[%d] = %4d     bounded[%d] = %4d\n",i,Spec[i],i,bounded[i]);
		}

		if (var_type == 0) {
			for (i=1; i<=4; i++) {
				fprintf(fp_log,"anasum[%d].SS = %13g  .MSE = %13g  .DF = %4d  .TEST = %13g\n",
					i,anasum[i].SS,anasum[i].MSE,anasum[i].DF, anasum[i].TEST);
			}
		}
		else {
			for (i=1; i<=5; i++) {
				fprintf(fp_log,"anasum[%d].SS = %13g  .MSE = %13g  .DF = %4d  .TEST = %13g\n",
					i,anasum[i].SS,anasum[i].MSE,anasum[i].DF, anasum[i].TEST);
			}
		}
	}
#endif

	/* Send 1 for var_type so we always assume modeled variance. */
	/* MJF, 02JUN2005. */
	/*compute and output ANOVA table elements*/
	DTMS3ANOVAC (nparm,Nobs,Spec,lkA3,xlk,lkA2, lkA1, lkR, anasum, 1, bounded);

#ifdef LOGGING_ON
	{
		fprintf(fp_log,"Variables coming out:\n");
		fprintf(fp_log,"nparm = %4d;  Nobs = %4d;  var_type = %3d;\n",nparm, Nobs, var_type);
		fprintf(fp_log,"lkA1 = %13g\n",lkA1);
		fprintf(fp_log,"lkA2 = %13g\n",lkA2);
		fprintf(fp_log,"lkA3 = %13g\n",lkA3);
		fprintf(fp_log,"lkR = %13g\n",lkR);
		fprintf(fp_log,"xlk = %13g          (fitted likelihood)\n",xlk);

		for (i=1; i<=nparm; i++) {
			fprintf(fp_log,"Spec[%d] = %4d     bounded[%d] = %4d\n",i,Spec[i],i,bounded[i]);
		}
		for (i=1; i<=5; i++) {
			fprintf(fp_log,"anasum[%d].SS = %13g  .MSE = %13g  .DF = %4d  .TEST = %13g\n",
				i,anasum[i].SS,anasum[i].MSE,anasum[i].DF, anasum[i].TEST);
		}
	}
#endif

	/* free memory */
	FREE_IVECTOR (bounded,1,nparm);
	/* Send 1 for var_type so we always assume modeled variance. */
	/* MJF, 02JUN2005. */
	/*output ANOVA table*/
	OUTPUT_DTMS3ANOVAC(anatxt,anasum, var_type);

#ifdef LOGGING_ON
	{
		fprintf(fp_log,"\n************ Call to Goodness *************\n");
	}
#endif
	/*print a goodness of fit table*/
	Goodness(nparm, nparm_known, Parms, var_type, anasum);

#ifndef RBMDS
	/*compute benchmark dose*/
	fprintf (fp_out2, "\n BMD_flag \t %d \n Nobs \t%d \n nparm \t%d",  bmdose, Nobs, nparm );
	fprintf (fp_out2, "\n  Con_lev \t%3.3g ", bmdparm.level);
	fprintf (fp_out2, "\n  BMRType \t%d ", bmr_type);
	fprintf (fp_out2, "\n  BMRF \t%f ", bmdparm.effect);
#endif
	sp = 0;

	if(Spec[2] == 1 && Parms[2] == 0){
		for(i = 1; i <= Nobs; i++){
			sp += ( (Ni[i] - 1)*Yd[i]/(Ntot - Nobs) );
		}
	}
#ifndef RBMDS
	for (i=1;i<=nparm; i++)
		fprintf (fp_out2, "\n %s \t %5.3g", Parm_name[i-1], Parms[i]);

	fprintf (fp_out2,"\n\n Data");

	for (i=1;i<=Nobs;i++) {
		Nd = Ni[i];
		lep=Ym[i] + qstudt(0.025, Nd - 1) * sqrt(Yd[i]/Nd);
		upep=Ym[i] + qstudt(0.975, Nd - 1) * sqrt(Yd[i]/Nd);

		fprintf (fp_out2,"\n %f %f %f %f", Xi[i], Ym[i], lep, upep);
	}

	fprintf (fp_out2,"\n Max_Min_dose \n  %f %f ", xmax, xmin);
#endif

#ifdef LOGGING_ON
	{
		fprintf(fp_log,"\n************* Call to Poly_BMD (if bmdose = 1) **************\n");
	}
#endif

	if (bmdose==Yes){
		Poly_BMD (nparm, Parms, EPS, &junk, xlk, Rlevel, Bmdl, &BMD, &bmdflag);
	}

#ifdef LOGGING_ON
	{
		fprintf(fp_log,"Variables coming out(and going into ProfLik):\n");
		fprintf(fp_log,"nparm = %3d\n",nparm);
		for (i=1; i<=nparm; i++) {
			fprintf(fp_log,"Parms[%d] = %13g         Spec[%d] = %d\n",
				i, Parms[i], i, Spec[i]);
		}
	}
#endif

	/*  Compute profile likelihood */

	/* if bmd is set to maxdose*100, then proflik is not called */

	/*    if (bmdflag ==0) */
	/*     */
	/*    HERE IS WHERE THE PROFILE STUFF STARTS. REMEMBER TO UNCOMMENT AFTER THE BMD STARTS TO WORK  */
	/*    if (bmdflag != 0) */
	/*      BMD = xmax; */

	/*    for (i =0 ; i < nparm; i++) */
	/*      { */
	/*        Parms[i+1]= Parmscopy[i]; */
	/*      }  */
	/*              Parms*/
	/*    ProfLik (nparm, IniP, xlk, BMD, Bmdl[1], argv[1], bmdparm.level, bmdparm.effect,&ProfFlag,&PEFlag,&PE_bmdl,&EmptyFlag); */
	/* If ProfLik does not plot to the bmdl then try new initial parameters */
	/* ProfFlag    0 = good plot                              */
	/*             1 = bad plot                               */
	/* PEFlag  0 = got an estimate for the bmdl from profile  */
	/*         1 = no estimate for bmdl from profile          */
	/* EmptyFlag  0 = more than 2 points in plot file         */
	/*            1 = less than 3 points in plot file         */
	/*    if(EmptyFlag == 1) */
	/*      { */
	/*        for(i=1; i<=nparm; i++) */
	/*  	{ */
	/*  	  Parms[i]= 1; */
	/*  	} */
	/*        ProfLik(nparm, Parms, xlk, BMD, Bmdl[1], argv[1],bmdparm.level, bmdparm.effect,  */
	/*  	      &ProfFlag,&PEFlag,&PE_bmdl,&EmptyFlag);  */
	/*      } */
	/*    if(EmptyFlag == 1)  */
	/*      { */
	/*        for(i=1; i<=nparm; i++) */
	/*  	{ */
	/*  	  Parms[i]= -1; */
	/*  	} */
	/*        ProfLik(nparm, Parms, xlk, BMD, Bmdl[1], argv[1],bmdparm.level, bmdparm.effect,  */
	/*  	      &ProfFlag,&PEFlag,&PE_bmdl,&EmptyFlag);  */
	/*      } */
	/*    if (PEFlag == 0) */
	/*      fprintf(fp_out,"\n         Profile = %14.6g\n",PE_bmdl); */
	/*    if (ProfFlag == 1) */
	/*      fprintf(fp_out,"\nThe profile plot or the BMDL may be inaccurate.\n\n"); */
	/*    else */
	/*      ERRORPRT("BMD computation failed. Profile likelihood function will not be calculated."); */

	FREE_DVECTOR (Parmscopy, 1, nparm);
	FREE_IVECTOR (Ni,1,Nobs);
	FREE_DVECTOR (Xi,1,Nobs);
	FREE_DVECTOR (Ym,1,Nobs);
	FREE_DVECTOR (Yd,1,Nobs);
	FREE_DVECTOR (IniP,1,nparm);
	FREE_IVECTOR (Spec,1,nparm);
	FREE_IVECTOR (IniSp,1,nparm);
	FREE_VLVECTOR(varsum,1,3);

	if(var_type == 0)
		FREE_ALVECTOR(anasum,1,4);
	else
		FREE_ALVECTOR(anasum,1,5);

	FREE_DVECTOR (Rlevel,1,5);
	FREE_DVECTOR (Bmdl, 1, 5);
	FREE_DMATRIX (vcv,1,nparm,1,nparm);
	FREE_DMATRIX (vcv_adj, 1, adj_vcv_rows, 1, adj_vcv_rows);
	FREE_DVECTOR (Parms, 1, nparm);
#ifdef LOGGING_ON
	{

		fflush(fp_log);
		fclose (fp_log);
	}
#endif
	CLOSE_FILES ();

	return(0);
} /*end of main*/


/******************************************************************
**Poly_vcv -- used to compute the vcv for Polynomial model.
*           Extern var.: smean, smax, Nobs, Xi, Yp, Yn, Ls, Xg.
*
******************************************************************/
void Poly_vcv(int nparm, int Spec[], double p[], double **vcv)
{
	void F1iDoublePart(int nparm,int const_var,double p[], double **Fn1i,int obs);
	void F2iDoublePart(int nparm,int const_var,double p[], double **Fn2i,int obs);
	void F3iDoublePart(int nparm,int const_var,double p[], double **Fn3i,int obs);
	void MeanPart(int obs, double *p, double *mg);
	void VarPart(int obs,int const_var,double Vi, double meani, double *p,
		double *mg, double *vg);
	void Mean2Part(int obs, double *p, double **mg2);
	void Var2Part(int obs,int const_var,double Vi, double meani, double *p,
		double *mg, double **mg2, double **vg2);


	double  **Fn1i,**Fn2i,**Fn3i,numi;
	int i, j, k, const_var;

	Fn1i = DMATRIX(1,nparm,1,nparm);
	Fn2i = DMATRIX(1,nparm,1,nparm);
	Fn3i = DMATRIX(1,nparm,1,nparm);

	/* Compute partials at parameter estimates */

	if(Spec[2] == 1 && p[2] == 0)
		const_var = 1;
	else
		const_var = 0;

	for(i = 1; i <= Nobs; i++)
	{
		numi = Ni[i];
		for(j = 1; j <= nparm; j++)
		{
			for(k = 1;k <= nparm; k++)
			{
				F1iDoublePart(nparm,const_var,p,Fn1i,i);
				F2iDoublePart(nparm,const_var,p,Fn2i,i);
				F3iDoublePart(nparm,const_var,p,Fn3i,i);

				vcv[j][k] = vcv[j][k] + numi*Fn1i[j][k]/2;
				vcv[j][k] += (numi-1)*Yd[i]*Fn2i[j][k]/2;
				vcv[j][k] += numi*Fn3i[j][k]/2;

			}
		}
	}

	FREE_DMATRIX(Fn1i, 1,nparm,1,nparm);
	FREE_DMATRIX(Fn2i, 1,nparm,1,nparm);
	FREE_DMATRIX(Fn3i, 1,nparm,1,nparm);
}

/**************************************************************
*Poly_fit -- Used to "prepare" the data for further computation,
*            i.e. compute the extern variables, give the initial
*            parameters, etc. THEN fit the Polynomial model.
*            (In fact, these jobs could be done in main().)
*
***************************************************************/
void Poly_fit(int nparm, double p[], int *is_conv,
			  int *iter, double *fret, int *bounded)
{
	int    i, j;
	double *parms, *fitparms, *fitparms2, *fitparms3, *doses, *savep;
	double ll, ll2, ll3;
	double *beginp, temp1, temp2;
	double *svar, *means, maxdose;
	int    *SpBak, inc_flag;
	double *pBak /*, *tmy, *t, **tmv*/;
	double *bsv, *doses1, *means1;
	double dfpe, pe, **X, **XP, **XPX, *XPY;
	int *Spec2, *bind, *bind2, *bind3, *nanim, nresm;
	int optite, optite2, optite3, model_type;
	int nvar, nparms, restr, signs, flag,  save_op;
	double save_ll, new_ll;


	pBak=DVECTOR(1, nparm);
	SpBak=IVECTOR(1, nparm);

	X = DMATRIX(1, Nobs, 1, ndeg+1);
	XP = DMATRIX(1, ndeg + 1, 1, Nobs);
	XPX = DMATRIX(1, ndeg + 1, 1, ndeg + 1);
	XPY = DVECTOR(1, ndeg + 1);
	bsv = DVECTOR(1, ndeg + 1);

#ifdef LOGGING_ON
	{
		fprintf(fp_log, "\n\n!!!!!!!!!!!!!Inside Poly_Fit!!!!!!!!!!!!!\n\n");
		fflush(fp_log);
	}
#endif

	/* Polynomial model */
	model_type = 0;
	dfpe = 0;
	pe = 0;

	/** rescale Dose to be: 0 <= Dose <= 1 **/
	scale=1;


	for (j=1; j<=nparm; j++)
		pBak[j]=p[j];      /* save the input p[]. */

#ifdef LOGGING_ON
	{//log file output
		fprintf(fp_log,"pBak values \n");
		for (i=1; i<=nparm; i++)
			fprintf(fp_log,"pBak[%d] = %12.6g\n",i,pBak[i]);

		fprintf(fp_log, "Starting to find Initial Parameter Estimates.\n");
		fflush(fp_log);
	}
#endif

	/****** Obtain initial estimations for p[] ******/
	if(initial==Yes) {
		for(j=1; j<=nparm; j++)
			p[j]=IniP[j];
	}
	else { /* The following block of code gives initial starting
		   values as the ordinary least squares solution */

		for(i = 1; i <= Nobs; i++){
			X[i][1] = 1.0;
			for(j = 2; j <= ndeg + 1; j++){
				if(Xi[i] != 0)
					X[i][j] = pow(Xi[i], j-1);
				else
					X[i][j] = 0;
			}
		}

		TRANSPOSE(X, XP, Nobs, ndeg + 1);

		MATMPYM2(XP, X, XPX, ndeg + 1, Nobs, ndeg + 1);

		INVMAT(XPX, ndeg + 1);

		MATMPYV2(ndeg + 1, Nobs, XP, Ym, XPY);

		MATMPYV2(ndeg + 1, ndeg + 1, XPX, XPY, bsv);

#ifdef LOGGING_ON
		{//log file output
			fprintf(fp_log,"p value BEFORE for(i = 1; i <= ndeg + 1; i ++).\n");
			for (i=1; i<=nparm; i++)
				fprintf(fp_log,"p[%d] = %12.6g\n",i,p[i]);
		}
#endif

		for(i = 1; i <= ndeg + 1; i ++)
			p[i+2] = bsv[i];

		for (i=4; i<=ndeg+3; i++){
			if (p[i]*restrict < 0)
				p[i]=0.0;//.1*restrict;
		}

		for(i = 1; i <= Nobs; i++) {
			pe += (Ni[i]-1)*Yd[i];
			dfpe += Ni[i]-1;
		}

		/*GLN - 03/18/2005 to fix the specified value not printing see below block of code */
		for (i = 1; i <= nparm; i++) {
			if (Spec[i] == 0){
				if(i==1)
				{
					if (cons_var == 0)
						p[1] = log(pe/dfpe);
					else
						p[1] = pe/dfpe;
				}
				if(i==2)
					p[i]=0;
			}
			else {
				p[i] = pBak[i];
			}
		} /*GLN - 03/18/2005*/

#ifdef LOGGING_ON
		{//log file output
			fprintf(fp_log,"p value AFTER for(i = 1; i <= ndeg + 1; i ++)\n");
			for (i=1; i<=nparm; i++)
				fprintf(fp_log,"p[%d] = %12.6g\n",i,p[i]);
		}
#endif

		/* commented by GLN - 03/18/2005 to fix the specified value not printing see above block of code
		p[1]=pe/dfpe;
		p[2]=0;
		- end commented by GLN - 03/18/2005 */
		if (p[1] > 1000)
			p[1] = 1.0;

#ifdef LOGGING_ON
		{
			fprintf(fp_log, "Getting Better Starting Values.\n");
			fflush(fp_log);
		}
#endif

		/* get better starting values in the case when parameters are
		restricted to be positive and the data is not increasing. */
		if (restrict == 1) {
			inc_flag = 0;
			doses1 = DVECTOR(0, Nobs);
			means1 = DVECTOR(0, Nobs);
			for(i = 1; i <= Nobs; i++)
			{
				doses1[i] = Xi[i];
				means1[i] = Ym[i];
			} /* end for */
			/* rearrange data to be in the form
			doses1[1] <= doses1[2] ... <= doses1[Nobs] */
			for(i = 2; i <= Nobs; i++) {
				if (doses1[i] < doses1[1]) {
					temp1 = doses1[1];
					temp2 = means1[1];
					doses1[1] = doses1[i];
					means1[1] = means1[i];
					doses1[i] = temp1;
					means1[i] = temp2;
				}
			}  /* end for */
			/* determine if data is increasing */
			for(i = 1; i <= Nobs-1; i++) {
				if (means1[i+1] < means1[i]) {
					inc_flag = 1;
					break;
				}
			}
			FREE_DVECTOR(doses1, 0, Nobs);
			FREE_DVECTOR(means1, 0, Nobs);

			if (inc_flag == 1){
				for(i = 4; i <= nparm; i++)
					p[i] = 0.0;
			} /* end if */
		} /* if (restrict == 1) */

		OUTPUT_TEXT("\n\n                  Default Initial Parameter Values  ");
		fflush(fp_out);

		OUTPUT_Init(nparm, Spec, p, Parm_name);

		fflush(fp_out);

	} /* initial==Yes */

	FREE_IVECTOR(SpBak, 1, nparm);
	FREE_DMATRIX(X, 1, Nobs, 1, ndeg+1);
	FREE_DMATRIX(XP, 1, ndeg + 1, 1, Nobs);
	FREE_DMATRIX(XPX, 1, ndeg + 1, 1, ndeg + 1);
	FREE_DVECTOR(XPY, 1, ndeg + 1);
	FREE_DVECTOR(bsv, 1, ndeg + 1);

	/* Scale by the maximum dose */
	maxdose = Xi[1];
	for(i = 2; i <= Nobs; i++){
		if(maxdose < Xi[i])
			maxdose = Xi[i];
	}

	for(i = 1; i <= Nobs; i++){
		Xi[i] = Xi[i]/maxdose;
	}

#ifdef LOGGING_ON
	{
		fprintf(fp_log, "Scaled Dose Levels.\n");
		for (j=1; j<=Nobs; j++){
			fprintf(fp_log,"Xi[%d] = %6.4g\n",j,Xi[j]);
		}
		fflush(fp_log);
	}
#endif
	/* adjust the parameters due to scaling */
	for(i = 4; i <= nparm; i++){
		p[i] = p[i]*pow(maxdose,i-3);
	}

#ifdef LOGGING_ON
	{
		fprintf(fp_log, "Scaled Parameters.\n");
		for (j=4; j<=nparm; j++){
			fprintf(fp_log,"p[%d] = %6.4g\n",j,p[j]);
		}
		fflush(fp_log);
	}
#endif

	/**get specified parameters.**/
	for (j=1;j<=nparm;j++)
		if (Spec[j]==Yes)
			p[j]=pBak[j];

	nvar = Nobs;
	nparms = nparm;
	restr = restrict;
	signs = sign;

#ifdef LOGGING_ON
	{
		fprintf(fp_log, "Starting to allocate memory.\n");
		fflush(fp_log);
	}
#endif

	doses = DVECTOR(0, Nobs-1);
	means = DVECTOR(0, Nobs-1);
	svar = DVECTOR(0, Nobs-1);
	parms = DVECTOR(0, nparm-1);
	beginp = DVECTOR(0, nparm-1);
	fitparms = DVECTOR(0, nparm-1);
	fitparms2 = DVECTOR(0, nparm-1);
	fitparms3 = DVECTOR(0, nparm-1);
	savep = DVECTOR(0, nparm-1);

	bind = IVECTOR(0, nparm-1);
#ifdef LOGGING_ON
	{
		fprintf(fp_log, "After bind.\n");
		fflush(fp_log);
	}
#endif
	bind2 = IVECTOR(0, nparm-1);
#ifdef LOGGING_ON
	{
		fprintf(fp_log, "After bind2.\n");
		fflush(fp_log);
	}
#endif
	bind3 = IVECTOR(0, nparm-1);

#ifdef LOGGING_ON 
	{
		fprintf(fp_log, "After bind3.\n");
		fflush(fp_log);
	}
#endif

	nanim = IVECTOR(0, Nobs-1);
#ifdef LOGGING_ON
	{
		fprintf(fp_log, "After nanim.\n");
		fflush(fp_log);
	}
#endif

	Spec2 = IVECTOR(0, nparm-1);
#ifdef LOGGING_ON
	{
		fprintf(fp_log, "After Spec2.\n");
		fflush(fp_log);
		fprintf(fp_log, "Memory Allocation\n");
		fflush(fp_log);
		fprintf(fp_log, "Getting ready to load common block values.\n");
		fflush(fp_log);
	}
#endif

	for(i = 1; i <= Nobs; i++){
		nanim[i-1] = (int) Ni[i];
		doses[i-1] = Xi[i];
		means[i-1] = Ym[i];
		svar[i-1] = Yd[i];
	}

	for(i = 1; i <= nparm; i++){
		if(initial == Yes && Spec[i] == 0){
			p[i] = IniP[i];
			Spec2[i-1] = 0;
		}
		else
		{
			Spec2[i-1] = Spec[i];
		}
		parms[i-1] = p[i];
	}

	for(i = 0; i <= nparm-1; i++){
		beginp[i] = parms[i];
	}

	flag = 0; /* del0 = 1.0D-4, tau0 = 1.0D0  */
	save_ll = -9e20;
	save_op = -99;

	/* This is the first call to getmle and internally in donlp2   */
	/* the parameters are scaled by the abs of their initial value */
	/* So this will give the next call different parameters to start w/ */
	getmle_(&nparms, parms, fitparms2, &ll2, &optite2, &nresm, bind2, &flag);

	/* This uses the fitparms from the scaled call as the starting values. */
	getmle_(&nparms, fitparms2, fitparms3, &ll3, &optite3, &nresm, bind3, &flag);


#ifdef LOGGING_ON
	{
		fprintf(fp_log,"Before the 1st call to getmle in Poly_Fit\n");
		for (j=0; j<=nparm-1; j++){
			fprintf(fp_log,"parms[%d] = %6.4g\n",j,parms[j]);
		}
		for (j=0; j<=Nobs-1; j++) {
			fprintf(fp_log,"means[%d] = %6.4g     svar[%d] = %6.4g\n",j,means[j],j,svar[j]);
		}
		fprintf(fp_log,"nvar = %d\n",nvar);
		fprintf(fp_log,"nparms = %d\n",nparms);
		fprintf(fp_log,"restr = %d\n",restr);
		fprintf(fp_log,"signs = %d\n",signs);
		fprintf(fp_log,"model_type = %d\n", model_type);
		fflush(fp_log);
	}
#endif
	getmle_(&nparms, parms, fitparms, &ll, &optite, &nresm, bind, &flag);

#ifdef LOGGING_ON
	{
		fprintf(fp_log,"After the 1st call to getmle in Poly_Fit\n");
		for (j=0; j<=nparm-1; j++){
			fprintf(fp_log,"fitparms[%d] = %6.4g  bind[%d] = %d\n",j,
				fitparms[j],j, bind[j]);
		}
		fprintf(fp_log,"nvar = %d\n",nvar);
		fprintf(fp_log,"nparms = %d\n",nparms);
		fprintf(fp_log,"restr = %d\n",restr);
		fprintf(fp_log,"signs = %d\n",signs);
		fprintf(fp_log,"model_type = %d\n", model_type);
		fprintf(fp_log,"optite = %d\n",optite);
		fprintf(fp_log,"ll = %6.4g\n",ll);
		fflush(fp_log);
	}
#endif
	if(optite < 0 || optite > 2){
		for(i = 0; i < 6; i++){
			if(optite != 3)
				GetNewParms(beginp, nparm);  /* Get a new starting point */
			else
			{
				for(j = 0; j <= nparm-1; j++)
					beginp[j] = fitparms[j];
			} /* end else */

			/* Try again */
#ifdef LOGGING_ON
			{
				fprintf(fp_log,"Before the %dth call to getmle in Poly_Fit\n",i+2);
				for (j=0; j<=nparm-1; j++){
					fprintf(fp_log,"beginp[%d] = %6.4g\n",j,beginp[j]);
				}
				fprintf(fp_log,"nvar = %d\n",nvar);
				fprintf(fp_log,"nparms = %d\n",nparms);
				fprintf(fp_log,"restr = %d\n",restr);
				fprintf(fp_log,"signs = %d\n",signs);
				fprintf(fp_log,"model_type = %d\n", model_type);
				fflush(fp_log);
			}
#endif
			getmle_(&nparms, beginp, fitparms, &ll, &optite, &nresm, bind, &flag);

#ifdef LOGGING_ON
			{
				fprintf(fp_log,"After the %dth call to getmle in Poly_Fit\n",i+2);
				for (j=0; j<=nparm-1; j++){
					fprintf(fp_log,"fitparms[%d] = %6.4g\n",j,fitparms[j]);
				}
				fprintf(fp_log,"nvar = %d\n",nvar);
				fprintf(fp_log,"nparms = %d\n",nparms);
				fprintf(fp_log,"restr = %d\n",restr);
				fprintf(fp_log,"signs = %d\n",signs);
				fprintf(fp_log,"model_type = %d\n", model_type);
				fprintf(fp_log,"optite = %d\n",optite);
				fprintf(fp_log,"ll = %6.4g\n",ll);
				fflush(fp_log);
			}
#endif

			new_ll = ll;
			if (new_ll > save_ll){
				save_ll = new_ll;
				for (j=0; j <= nparm-1; j++)
					savep[j] = fitparms[j];
				save_op = optite;
			}


			if(optite >= 0 && optite <= 2)
				break;
		} /* end for */
	} /* end if */

	flag = 1; /* del0 = 1.0D-2, tau0 = 1.0D-6 */
	if(optite < 0 || optite > 2){
		for(j = 0; j <= nparm-1; j++)
			beginp[j] = parms[j];

		for(i = 0; i < 6; i++){
			if(optite != 3)
				GetNewParms(beginp, nparm); /*Get a new starting point*/
			else{
				for(j = 0; j <= nparm-1; j++)
					beginp[j] = fitparms[j];
			} /* end else */

			/* Try again */
#ifdef LOGGING_ON
			{
				fprintf(fp_log,"After the %dth call to getmle in Poly_Fit\n",i+8);
				for (j=0; j<=nparm-1; j++){
					fprintf(fp_log,"beginp[%d] = %6.4g\n",j,beginp[j]);
				}
				fprintf(fp_log,"nvar = %d\n",nvar);
				fprintf(fp_log,"nparms = %d\n",nparms);
				fprintf(fp_log,"restr = %d\n",restr);
				fprintf(fp_log,"signs = %d\n",signs);
				fprintf(fp_log,"model_type = %d\n", model_type);
				fflush(fp_log);
			}
#endif
			getmle_(&nparms, beginp, fitparms, &ll, &optite,&nresm, bind, &flag);

#ifdef LOGGING_ON
			{
				fprintf(fp_log,"After the %dth call to getmle in Poly_Fit\n",i+8);
				for (j=0; j<=nparm-1; j++){
					fprintf(fp_log,"fitparms[%d] = %6.4g\n",j,fitparms[j]);
				}
				fprintf(fp_log,"nvar = %d\n",nvar);
				fprintf(fp_log,"nparms = %d\n",nparms);
				fprintf(fp_log,"restr = %d\n",restr);
				fprintf(fp_log,"signs = %d\n",signs);
				fprintf(fp_log,"model_type = %d\n", model_type);
				fprintf(fp_log,"optite = %d\n",optite);
				fprintf(fp_log,"ll = %6.4g\n",ll);
				fflush(fp_log);
			}
#endif
			new_ll = ll;
			if (new_ll > save_ll){
				save_ll = new_ll;
				for (j=0; j <= nparm-1; j++)
					savep[j] = fitparms[j];
				save_op = optite;
			}

			if(optite >= 0 && optite <= 2)
				break;

		} /* end for */

	} /* end if */

	if(optite < 0 || optite > 2){
		for(j = 0; j <= nparm-1; j++)
			beginp[j] = parms[j];

		for(flag = 0; flag <= 1; flag++){

			for(i = 0; i < 3; i++){
				GetMoreParms(beginp, nparm);

				/* Try again */
#ifdef LOGGING_ON
				{
					fprintf(fp_log,"After the %dth call to getmle in Poly_Fit\n",i+14);
					for (j=0; j<=nparm-1; j++){
						fprintf(fp_log,"beginp[%d] = %6.4g\n",j,beginp[j]);
					}
					fprintf(fp_log,"nvar = %d\n",nvar);
					fprintf(fp_log,"nparms = %d\n",nparms);
					fprintf(fp_log,"restr = %d\n",restr);
					fprintf(fp_log,"signs = %d\n",signs);
					fprintf(fp_log,"model_type = %d\n", model_type);
					fflush(fp_log);
				}
#endif
				getmle_(&nparms, beginp, fitparms, &ll, &optite,&nresm, bind, &flag);

#ifdef LOGGING_ON
				{
					fprintf(fp_log,"After the %dth call to getmle in Poly_Fit\n",i+14);
					for (j=0; j<=nparm-1; j++){
						fprintf(fp_log,"fitparms[%d] = %6.4g\n",j,fitparms[j]);
					}
					fprintf(fp_log,"nvar = %d\n",nvar);
					fprintf(fp_log,"nparms = %d\n",nparms);
					fprintf(fp_log,"restr = %d\n",restr);
					fprintf(fp_log,"signs = %d\n",signs);
					fprintf(fp_log,"model_type = %d\n", model_type);
					fprintf(fp_log,"optite = %d\n",optite);
					fprintf(fp_log,"ll = %6.4g\n",ll);
					fflush(fp_log);
				}
#endif

				new_ll = ll;
				if (new_ll > save_ll){
					save_ll = new_ll;
					for (j=0; j <= nparm-1; j++)
						savep[j] = fitparms[j];
					save_op = optite;
				}

				if(optite >= 0 && optite <= 2){
					flag = 2;
					break;
				}
			} /* end for */
		}
	}

	flag = 0;
	if(optite < 0 || optite > 2){
		for(i = 0; i <= nparm-1; i++){
			beginp[i] = p[i+1];
		}

		for(i = 0; i < 3; i++){
			GetOtherParms(beginp, nparm);  /*Get a new starting point*/

			/* Try again */
#ifdef LOGGING_ON
			{
				fprintf(fp_log,"After the %dth call to getmle in Poly_Fit\n",i+17);
				for (j=0; j<=nparm-1; j++){
					fprintf(fp_log,"beginp[%d] = %6.4g\n",j,beginp[j]);
				}
				fprintf(fp_log,"nvar = %d\n",nvar);
				fprintf(fp_log,"nparms = %d\n",nparms);
				fprintf(fp_log,"restr = %d\n",restr);
				fprintf(fp_log,"signs = %d\n",signs);
				fprintf(fp_log,"model_type = %d\n", model_type);
			}
#endif
			getmle_(&nparms, beginp, fitparms, &ll,&optite, &nresm, bind, &flag);

#ifdef LOGGING_ON
			{
				fprintf(fp_log,"After the %dth call to getmle in Poly_Fit\n",i+17);
				for (j=0; j<=nparm-1; j++){
					fprintf(fp_log,"fitparms[%d] = %6.4g\n",j,fitparms[j]);
				}
				fprintf(fp_log,"nvar = %d\n",nvar);
				fprintf(fp_log,"nparms = %d\n",nparms);
				fprintf(fp_log,"restr = %d\n",restr);
				fprintf(fp_log,"signs = %d\n",signs);
				fprintf(fp_log,"model_type = %d\n", model_type);
				fprintf(fp_log,"optite = %d\n",optite);
				fprintf(fp_log,"ll = %6.4g\n",ll);
			}
#endif
			new_ll = ll;
			if (new_ll > save_ll){
				save_ll = new_ll;
				for (j=0; j < nparm; j++)
					savep[j] = fitparms[j];
				save_op = optite;
			}

			if(optite >= 0 && optite <= 2)
				break;

		} /* end for */

	} /* end if */

	/* This decides if the scaling model is better than the unscaled model */
	/* or not. */
	if ((optite2 >= 0) && (optite2 <= 2)){
		if (ll2 > ll){
			for (j = 0; j < nparm; j++)
			{
				fitparms[j] = fitparms2[j];
				bind[j] = bind2[j];
			}
			optite = optite2;
			ll = ll2;
		}
	}

	if ((optite3 >= 0) && (optite3 <= 2)){
		if (ll3 > ll){
			for (j = 0; j < nparm; j++){
				fitparms[j] = fitparms3[j];
				bind[j] = bind3[j];
			}
			optite = optite3;
			ll = ll3;
		}
	}

	/* Let user know if no optimum was found */
	if(optite < 0 || optite > 2){
		fprintf(fp_out, "\n\n!!! Warning:  optimum may not have been found. !!!");
		/*    fprintf(fp_out, "\n!!! Bad completion code in maximum likelihood optimization routine  !!!");*/
		/*   fprintf(fp_out, "\n!!! Program halted
		!!!\n\n"); */
		fprintf(fp_out, "\n!!! You may want to try choosing different initial values.  !!!");

		/*    *is_conv = -1; */
	} /* end if */
	else
		*is_conv = 1;

	for(i = 1; i <= nparm; i++){
		if (optite < 0 || optite > 2)
			p[i] = savep[i-1];
		else
			p[i] = fitparms[i-1];

		bounded[i] = bind[i-1]; /* put back into the program */
		if(bounded[i] != 1 && bounded[i] != 0)
		{
			bounded[i] = 1; /* Hopefully, we never get here; originally, this was 0, but 1 looks
							like a better warning flag */
		}
	} /* end for */

	for(i = 1; i <= Nobs; i++)
		Xi[i] = Xi[i]*maxdose;

	for(i = 4; i <= nparm; i++)	{
		p[i] = p[i]/pow(maxdose,i-3);
	}

	/* This looks like an attempt to second guess donlp2; possibly triggered by errors in indexing bind()
	in getmle_().  RWS 10-19-2005
	bounded[1] = 0;
	bounded[3] = 0;
	if ((fabs(p[2]-18) < 1e-20) || (fabs(p[2]+18) < 1e-20))
	bounded[2] = 1;
	else
	bounded[2] = 0;
	*/

	if (optite < 0 || optite > 2)
		*fret = save_ll;
	else
		*fret = ll;

	FREE_DVECTOR(pBak, 1, nparm);

	FREE_DVECTOR(doses, 0, Nobs-1);
	FREE_DVECTOR(means, 0, Nobs-1);
	FREE_DVECTOR(svar, 0, Nobs-1);
	FREE_IVECTOR(nanim, 0, Nobs-1);
	FREE_DVECTOR(parms, 0, nparm-1);
	FREE_DVECTOR(fitparms, 0, nparm-1);
	FREE_DVECTOR(fitparms2, 0, nparm-1);
	FREE_DVECTOR(fitparms3, 0, nparm-1);
	FREE_DVECTOR(beginp, 0, nparm-1);
	FREE_IVECTOR(Spec2, 0, nparm-1);
	FREE_IVECTOR(bind, 0, nparm-1);
	FREE_IVECTOR(bind2, 0, nparm-1);
	FREE_IVECTOR(bind3, 0, nparm-1);
	FREE_DVECTOR(savep, 0, nparm-1);

} /* end of Poly_fit */

/************************************************************
* Poly_BMD -- Used to calculate the BMD and BMDL for Polynomial model.
*
*************************************************************/
void Poly_BMD (int nparm, double p[], double gtol, int *iter, double xlk,
			   double Rlevel[], double Bmdl[], double *BMD, int *bmderror)
{
	double BMD_func(int nparm, double p[], double x, double gtol);
	double BMD1_func(int nparm, double p[], double x, double gtol);

	double   tol, Drange, poly, fb1;
	double   xa,xb,fa,fb, bDose;
	double   D, bmrtemp, temp, temp_poly;
	double   *pBak, *BMRVals;
	double   incre, divide, thedose, tempdose[2], bmdlmean[2];
	int      j, k, bmrfncsign, stop_flag, Errorbmdlc ;



	pBak=DVECTOR(1, nparm);
	BMRVals=DVECTOR(1,5);

	for(j=1; j<=nparm; j++)
		pBak[j]= p[j];          /* save the p[]. */

	/**** compute X^2 value  ************************/
	/* IF ML is the value of the maximized log-likelihood, then ML - LR is the value
	log-likelihood at the BMDL or BMDU */
	if (bmdparm.level<0.5)
		LR = 0.5*QCHISQ(1.0 - 2.0 * bmdparm.level, 1);
	else
		LR = 0.5*QCHISQ(2.0 * bmdparm.level - 1.0, 1);

	if(bmr_type != 3 && bmr_type != 4)
		Rlevel[1] = BMR = fabs(bmdparm.effect);
	else
		Rlevel[1] = BMR = bmdparm.effect;

	xa = D = 0.0;
	tol=0.0000001;
	incre=1.0;
	divide=50.0;
	bmrfncsign = 1;
	fa = -1;


	/* Get a BMR value that is equal to just Mu(BMD) regardless
	of the BMR type to simplify the root finding routine */

	if(bmr_type == 0 || bmr_type == 1 || bmr_type == 2)
	{
		if(sign == -1)
			bmrtemp = -BMR;
		else
			bmrtemp = BMR;

		if(bmr_type == 1)
		{
			if(p[2] != 0)
				bmrtemp = bmrtemp*sqrt(exp(p[1] + log(p[3]) * p[2]));
			else
				bmrtemp = bmrtemp*sqrt(p[1]);
		}
		else if(bmr_type == 2)
		{
			bmrtemp = bmrtemp*p[3];
		}
	}
	else
	{
		if(bmr_type == 3)
		{
			bmrtemp = BMR - p[3];
		}
		else
			bmrtemp = p[4]*BMR;
	}
	/* bmrtype is now a Point BMR type regardless of the
	user input type */
	bmrtemp = bmrtemp + p[3];

	temp = BMR;
	BMR = bmrtemp;

	/*************** solve the BMD *********************/

	ck = BMR;
	xa = D = 0.0;
	Drange = xmax;
	stop_flag = 0;
	fb = p[3] - ck;
	*bmderror = 0;
	k=1;
	while ((k<10000) && (stop_flag == 0))
	{
		xa=D;
		D=Drange*k/100;
		poly=p[nparm];
		for (j=nparm-1; j>=3; j--) poly = poly*D+p[j];
		fb1 = poly - ck;
		if ((fb == 0) || (fb1 == 0))
			stop_flag = 1;
		if ( ((fb > 0) && (fb1 < 0)) || ((fb < 0) && (fb1 > 0)) )
			stop_flag = 1;
		fa = fb;
		k++;
		fb = fb1;
	}
	if(k == 10000)
	{
		*bmderror = 1;
		/*  printf("\nBMD computation failed for BMR = %14.6g", BMR); */
		fprintf(fp_out, "\nBMD computation failed for BMR = %14.6g", BMR);
		fprintf(fp_out, "\nSetting BMD = 100*(maximum dose)");
		thedose = 100*xmax;
		// for testing only
		*BMD = -9999/scale;
		OUTPUT_BENCHMD2(1,(*BMD)*scale);
	} /* end if */
	else
	{
		xb = D + 0.01*xmax;
		tol=0.0000001;

		xb = zeroin(xa, xb, tol, BMD1_func, nparm, p, tol);
		thedose = xb;
		// for testing only
		*BMD = thedose;
		OUTPUT_BENCHMD2(1,(*BMD)*scale);
	} /* end else */

	/* for testing only; when testing is done, unblock this section and
	delete the two testing only sections above

	*BMD = thedose;
	OUTPUT_BENCHMD2(1,(*BMD)*scale);
	*/

	/*********************************************/
	bDose = thedose;

	BMR = temp;
	BMRVals[1] = BMR; /* ********** */

	temp_poly = p[nparm];
	for (j=nparm-1; j>=3; j--)
		temp_poly = temp_poly*bDose + p[j];
	react = temp_poly;

#ifndef RBMDS
	fprintf (fp_out2, "\n  RSL \t%f",react);
	fprintf (fp_out2, "\n  BMD \t%f",bDose);
	fprintf (fp_out2,"\n\n BMD_line");
	fprintf (fp_out2,"\n %f %f", -1.0, react);
	fprintf (fp_out2,"\n %f %f", bDose, react);
	fprintf (fp_out2,"\n %f %f", bDose, -0.1);
#endif
	if (bDose <= 0.0)
	{
#ifndef RBMDS
		fprintf (fp_out2, "\n\n BMDL_comput_ind %d", No);  /*computation failed.*/
#endif
		Warning("BMDL computation failed.  Lower limit includes zero.");
		/*        Warning("BMDL is set to be 0.0000001 for lilelihood function calculation."); */
		// BMDL IS SET TO BE 0.0000001 FOR GETPROFILE
		/*        Bmdl[1]=0.0000001; */
	}

	/********* search for BMDL **************************/
	tol = FMAX((*BMD)*0.001, 0.0000001);

	xa = thedose/100.0;
	xb = thedose;
	BMD_lk = xlk;       /* get the lk at BMD. */

	fb = -LR;

	fa = BMDL_func(nparm, BMD_lk, xb, p, tol);

	if (fa<0.0)
	{   /* computation failed. */
#ifndef RBMDS
		Bmdl[1] = -1.0; 
		fprintf (fp_out2, "\n\n BMDL_comput_ind %d", No);
#endif
		Warning("BMDL computation failed.");
		/*        Warning("BMDL is set to be 0.0000001 for lilelihood function calculation."); */

		// BMDL IS SET TO BE 0.0000001 FOR GETPROFILE
		/*        Bmdl[1] = 0.0000001; */
	}
	else
	{
#ifndef RBMDS
		/* computation will succeed. */
		fprintf (fp_out2, "\n\n BMDL_comput_ind %d", Yes);
#endif
		Bmdl[1] = fa;

		tempdose[1] = thedose;

		for(j=1; j<=nparm; j++)
			p[j] = pBak[j];          /* get the "old" p[]. */

		PolyMeans(1, p, tempdose, bmdlmean);

		Rlevel[1]= bmdlmean[1];

		fprintf(fp_out, 
#ifndef RBMDS
			"\n            BMDL = %14.6g\n\n"
#else
			"\n            BMDL = %30.22g\n\n"
#endif
			, Bmdl[1]);
#ifndef RBMDS
		fprintf (fp_out2, "\n  BMDL \t%f",Bmdl[1]);
		fprintf (fp_out2,"\n\n BMDL_line");
		fprintf (fp_out2,"\n %f %f", Bmdl[1], -0.1);
		fprintf (fp_out2,"\n %f %f", Bmdl[1], react);
#endif
		fflush(fp_out);
	} // endif
	fflush(fp_out2);

	// bmdlCurve = No;
	if (bmdlCurve==Yes)
	{

		Get_BMRS(p, bDose, Bmdl[1], BMRVals, sign, bmr_type);


		/* Now BMRVals[2]...BMRVals[5] contain Point BMRs for BMDL curve
		point calculations.  */
		Errorbmdlc = 0;
		for (k=2; k<=5;k++)
		{
			/**** solve the BMD ********************************/
			if (Errorbmdlc == 0)
			{

				for(j=1; j<=nparm; j++)
				{
					p[j]= pBak[j];               /* Get back p[] */
				}

				xa = D = 0.0;

				tol=0.0000001;
				incre=1.0;

				fa = -1;

				/* Get a BMR value that is equal to just Mu(BMD) regardless
				of the BMR type to simplify the root finding routine */

				if(bmr_type == 0 || bmr_type == 1 || bmr_type == 2)
				{
					if(sign == -1)
						bmrtemp = -BMRVals[k];
					else
						bmrtemp = BMRVals[k];

					if(bmr_type == 1)
					{
						if(p[2] != 0)
							bmrtemp = bmrtemp*sqrt(p[1]*pow(p[3],p[2]));
						else
							bmrtemp = bmrtemp*sqrt(p[1]);
					}
					else if(bmr_type == 2)
					{
						bmrtemp = bmrtemp*p[3];
					}
				}
				else
				{
					if(bmr_type == 3)
					{
						bmrtemp = BMRVals[k] - p[3];
					}
					else
						bmrtemp = p[4]*BMRVals[k];
				}

				bmrtemp = bmrtemp + p[3];     /* bmrtype is now a Point BMR type
											  regardless of the user input type */

				BMR = bmrtemp;


				if (fa > 0)
				{
					fa = BMD_func(nparm,p,xa,tol);
					fb = BMD_func(nparm,p,xb,tol);
					thedose = zeroin (xa, xb, tol, BMD_func, nparm, p, tol);
				}

				/********* search for BMDL **************************/
				tol = FMAX((*BMD)*0.001, 0.0000001);

				xa = thedose/100.0;
				xb = thedose;
				BMD_lk = xlk;       /* get the lk at BMD. */

				fb = -LR;

				BMR = BMRVals[k];

				fa = BMDL_func(nparm, BMD_lk, xb, p, tol);

				if (fa<0.0)
				{
					Errorbmdlc=1;
#ifndef RBMDS
					fprintf (fp_out2, "\n\n BMDL_Curve_flag \t %d  \n smooth_opt  %d", 0, smooth);
#endif
					Warning("BMDL computation failed for one or more point on the BMDL curve. \n The BMDL curve will not be plotted\n");
					fflush(fp_out2);
					fflush(fp_out);
				}


				Bmdl[k] = fa;
				tempdose[1] = thedose;

				for(j=1; j<=nparm; j++)
					p[j] = pBak[j];          /* get the "old" p[]. */

				PolyMeans(1, p, tempdose, bmdlmean);
				Rlevel[k]= bmdlmean[1];
			}
		}
	}
	else 
		for (k=2; k<=5;k++)
			Bmdl[k] = Rlevel[k]= -1;

	for(j=1; j<=nparm; j++)
		p[j]= pBak[j];
#ifndef RBMDS
	fprintf (fp_out2, "\n\n BMDL_Curve_flag \t %d  \n smooth_opt  %d", bmdlCurve, smooth);
	fprintf (fp_out2,"\n\n BMDL_curve");
	fprintf (fp_out2,"\n 0.00000 %f", p[3]);
	for (k=1;k<=5;k++)
		fprintf (fp_out2,"\n %f %f", Bmdl[k], Rlevel[k]);

	fflush(fp_out2);
#endif
	FREE_DVECTOR(BMRVals, 1, 5);
	FREE_DVECTOR(pBak, 1, nparm);
}  /*end BENCHMD*/

/******************************************************************
* Poly_func -- used to compute the value of the polynomial function
*              at the point xx, given parm p[] and nparm
*******************************************************************/
double Poly_func (int nparm, double pBak[], double xx, double gtol)
{
	double fP;
	int j;

	fP = pBak[nparm];

	for (j=nparm-1; j>=1; j--) fP=fP*xx+pBak[j];

	return fP;
}

double BMD1_func(int nparm, double pBak[], double D, double gtol)
{
	double fP, fx;
	int j;

	fP = pBak[nparm];
	for (j=nparm-1; j>=3; j--) fP = fP*D+pBak[j];
	fx= fP - ck;
	return fx;
}



double BMD_func(int nparm, double pBak[], double D, double gtol)
{
	double fd, fbmd, tD[2], m[2];
	double *p;
	int    j;

	/*tD = DVECTOR(1, 1);*/
	p = DVECTOR(1, nparm);
	/*m = DVECTOR(1, 1);*/
	for (j=1;j<=nparm;j++)
		p[j]=pBak[j];            /* get the "old" p[]. */

	tD[1] = D;

	PolyMeans(1, p, tD, m);

	fd = m[1];

	fbmd = fd - BMR;

	return fbmd;

	FREE_DVECTOR(p, 1, nparm);


}

/*****************************************************************
* BMDL_func -- used to compute the values of functions BMDL_f (the
*              X^2 value) at the point D, given the parm p[] and
*              number of parm. If GETCL finishes with bad
*              convergence, then the return value is -1
*
*****************************************************************/
double BMDL_func(int nparm, double xlk, double Dose, double pBak[], double gtol)
{    /*  ck , BMD_lk and LR are calculated in Multistage_BMD() */
	double fD, bmdl, *doses, target, *parms, *means, *beginp;
	double *svar, *fitparms,maxdose;
	int j, i, gccnt;
	int *Spec2, flag,nvar, nparms, restr, signs, bmrtype;
	int which,*nanim, *bind, nresm, optite, model_type;


	/* polynomial model */

	model_type = 0;
	/* flag is used to change the del0 and tau0 in donlp2usrfc.f */
	flag = 0;

	nvar = Nobs;
	restr = restrict;
	nparms = nparm;
	signs = sign;
	bmrtype = bmr_type;
	optite = -5;

	doses = DVECTOR(0, Nobs-1);
	means = DVECTOR(0, Nobs-1);
	svar = DVECTOR(0, Nobs-1);
	nanim = IVECTOR(0, Nobs-1);
	parms = DVECTOR(0, nparm-1);
	beginp = DVECTOR(0, nparm-1);
	fitparms = DVECTOR(0, nparm-1);
	Spec2 = IVECTOR(0, nparm-1);
	bind = IVECTOR(0, nparm-1);


	for(i = 1; i <= Nobs; i++)
	{
		nanim[i-1] = (int) Ni[i];
		doses[i-1] = Xi[i];
		means[i-1] = Ym[i];
		svar[i-1] = Yd[i];
	}

	which = 1;          /* Want an lower confidence limit */

	target = xlk - LR;  /* The value we want the likelihood
						at the BMDL to match */

	for (j=1;j<=nparm;j++)
	{
		parms[j-1]= pBak[j];  /* get the "old" p[]. */
		Spec2[j-1] = Spec[j];
		beginp[j-1] = pBak[j];
	}

	/* Set up and call subroutine that calculates BMDL */

	maxdose = doses[0];
	for(i = 1; i < Nobs; i++)
	{
		if(doses[i] > maxdose)
			maxdose = doses[i];
	}

	for(i = 0; i < Nobs; i++)
		doses[i] = doses[i]/maxdose;    /* Scale dose to be between 0
										and 1 */
	Dose = Dose/maxdose;

	/* adjust parameters effected by scaling */

	for(i = 3; i < nparm; i++)
	{
		parms[i] = parms[i]*pow(maxdose,i-2);
		beginp[i] = beginp[i]*pow(maxdose,i-2);
	}

	for(gccnt = 1; gccnt < 49; gccnt++)
	{
		if(gccnt == 12 || gccnt == 22 || gccnt == 25 || gccnt == 36 || gccnt == 46)
		{
			if(gccnt == 36 || gccnt == 46)
			{
				beginp[0] = pBak[1];
				beginp[1] = pBak[2];
				beginp[2] = pBak[3];
			}

			for(i = 3; i < nparm; i++)
				beginp[i] = pBak[i+1]*pow(maxdose,i-2);

			if(gccnt == 25)
			{
				flag = 1;
				beginp[0] = pBak[1];
				beginp[1] = pBak[2];
				beginp[2] = pBak[3];
			}
		}
		if((gccnt > 1 && gccnt < 12) || (gccnt > 25 && gccnt < 36))
		{
			if(optite != 4)
				GetNewParms(beginp, nparm);  /*Get a new starting point */
			else
				for(i = 0; i < nparm; i++)
					beginp[i] = fitparms[i];
		}
		if((gccnt > 11 && gccnt < 22) || (gccnt > 35 && gccnt < 46))
		{
			if(optite != 4)
				GetOtherParms(beginp, nparm);  /* Get a new starting point */
			else
				for(i = 0; i < nparm; i++)
					beginp[i] = fitparms[i];
		}
		if((gccnt > 21 && gccnt < 25) || (gccnt > 45 && gccnt < 49))
		{
			if(optite != 4)
				GetMoreParms(beginp, nparm);  /* Get a new starting point */
			else
				for(i = 0; i < nparm; i++)
					beginp[i] = fitparms[i];
		}
#ifdef LOGGING_ON
		{
			fprintf(fp_log,"\n***********Before Call #%d to getcl.****************\n",gccnt);
			fprintf(fp_log,"(this call scales the parameters by their starting values in donlp2)\n");
			fprintf(fp_log,"optite = %d\n",optite);
			fprintf(fp_log,"nparms = %d\n",nparms);
			fprintf(fp_log,"flag = %d\n",flag);
			fprintf(fp_log,"BMR = %10.5g\n",BMR);
			fprintf(fp_log,"bmrtype = %d\n",bmrtype);
			fprintf(fp_log,"Dose = %10.5g\n",Dose);
			fprintf(fp_log,"target = %10.5g\n",target);
			fprintf(fp_log,"bmdl = %10.5g  (This should be the BMD)\n",bmdl);
			if(gccnt == 1)
			{
				for(i = 0; i < nparm; i++)
					fprintf(fp_log,"parms[%d] = %12.5g\n", i, parms[i]);
			}
			else
			{
				for(i = 0; i < nparm; i++)
					fprintf(fp_log,"beginp[%d] = %12.5g\n", i, beginp[i]);
			}
			fprintf(fp_log,"************************************************\n");
			fflush(fp_log);
		}
#endif

		if(gccnt == 1)
			getcl_(&which,&nvar,doses,means,nanim,svar,&nparms,&BMR,
			&Dose,&target,parms,Spec2,parms,&bmrtype,
			&restr,&bmdl,fitparms,&optite,&nresm,bind,&signs,
			&model_type, &flag, pBak);
		else
			getcl_(&which,&nvar,doses,means,nanim,svar,&nparms,&BMR,
			&Dose,&target,beginp,Spec2,beginp,&bmrtype,
			&restr,&bmdl,fitparms,&optite,&nresm,bind,&signs,
			&model_type,&flag,pBak);

#ifdef LOGGING_ON
		{
			fprintf(fp_log,"\n***********After Call #%d to getcl.****************\n",gccnt);
			fprintf(fp_log,"(this call scales the parameters by their starting values in donlp2)\n");
			fprintf(fp_log,"optite = %d\n",optite);
			fprintf(fp_log,"nparms = %d\n",nparms);
			fprintf(fp_log,"flag = %d\n",flag);
			fprintf(fp_log,"BMR = %10.5g\n",BMR);
			fprintf(fp_log,"bmrtype = %d\n",bmrtype);
			fprintf(fp_log,"Dose = %10.5g\n",Dose);
			fprintf(fp_log,"target = %10.5g\n",target);
			fprintf(fp_log,"bmdl = %10.5g  (This should be the BMD)\n",bmdl);
			if(gccnt == 1)
			{
				for(i = 0; i < nparm; i++)
					fprintf(fp_log,"parms[%d] = %12.5g\n", i, parms[i]);
			}
			else
			{
				for(i = 0; i < nparm; i++)
					fprintf(fp_log,"beginp[%d] = %12.5g\n", i, beginp[i]);
			}
			fprintf(fp_log,"************************************************\n");
			fflush(fp_log);
		}
#endif

		if(optite >= 0 && optite <= 3)
			break;

	}// End of for(gccnt = 1; gccnt < 34; gccnt++)

	/* optite is a value that is passed back from GETCL which
	determines whether the optimization was completed successfully
	If optite is less than 0, then it did not, and we want
	to try a different starting point and recompute */
	if(optite < 0 || optite > 3)
	{
#ifdef LOGGING_ON
		{
			fprintf(fp_log,"\n***********Call #%d inside if(optite < 0 || optite > 3) NOT VALID OPTITE.****************\n",gccnt);
		}
#endif
		fD = -1;
	} /* end if */
	else
	{
		for (j=1;j<=nparm;j++)
		{
			pBak[j] = fitparms[j-1];
		} /* end for */
		for(i = 3; i <= nparm-1; i++)
		{
			parms[i] = parms[i]/pow(maxdose,i-3);
		} /* end for */
		Dose = Dose*maxdose;
		bmdl = bmdl*maxdose;
		/*  printf(" optite = %d\n", optite); */
		fD = bmdl;
#ifdef LOGGING_ON
		{
			fprintf(fp_log,"\n***********Call #%d inside else (optite >= 0 && optite <= 3).****************\n",gccnt);
			fprintf(fp_log,"optite = %d\n",optite);
			fprintf(fp_log,"nparms = %d\n",nparms);
			fprintf(fp_log,"flag = %d\n",flag);
			fprintf(fp_log,"BMR = %10.5g\n",BMR);
			fprintf(fp_log,"Dose = %10.5g\n",Dose);
			fprintf(fp_log,"bmdl = %10.5g  (This should be the BMD)\n",bmdl);
			for (j=1;j<=nparm;j++)
				fprintf(fp_log,"pBak[%d] = %12.5g, fitparms[%d-1]\n", j, pBak[j], j);

			for(i = 3; i <= nparm-1; i++)
				fprintf(fp_log,"parms[%d] = %12.5g\n", i, parms[i]);

			fprintf(fp_log,"************************************************\n");
			fflush(fp_log);
		}
#endif
	} /* end else */

	FREE_DVECTOR(doses, 0, Nobs-1);
	FREE_DVECTOR(means, 0, Nobs-1);
	FREE_DVECTOR(svar, 0, Nobs-1);
	FREE_IVECTOR(nanim, 0, Nobs-1);
	FREE_DVECTOR(parms, 0, nparm-1);
	FREE_DVECTOR(fitparms, 0, nparm-1);
	FREE_DVECTOR(beginp, 0, nparm-1);
	FREE_IVECTOR(Spec2, 0, nparm-1);
	FREE_IVECTOR(bind, 0, nparm-1);
	return fD;

} // End of double BMDL_func(int nparm, double xlk, double Dose, double pBak[], double gtol)

/**********************************************************
** READ_OBSDATA4V--used to read 4 column data in 4 vectors.
***********************************************************/
int READ_OBSDATA4V(int Nobs,double Xi[],int Ni[],double Ym[],double Yd[])

{
	int     Nmiss;          /*number of records with missing values*/
	int     i,j,n,m;        /*count and iteration control variables*/
	double  value;          /*temp variable*/

	Nmiss = 0;
	for(i=1;i<=Nobs;i++)
	{
		n = i-Nmiss;
		m = 0;
		for (j=1;j<=4;j++)
		{
			fscanf(fp_in,"%lf",&value);
			if (value != MISSING)
			{
				if (j==1)     Xi[n]=value;
				if (j==2)     Ni[n]=(int)value;
				if (j==3)     Ym[n]=value;
				if (j==4)     Yd[n]=value;
			}
			else
				m++;
		}
		if (m != 0)     Nmiss++;
		else if (Xi[n] < 0)      Nmiss++;
	}

	return Nmiss;
}


/**********************************************************
** READ_OBSDATA2V--used to read 2 column data in 2 vectors.
***********************************************************/
int READ_OBSDATA2V(int ntotal, double xxi[], double yyi[])

{
	int     Nmiss;          /*number of records with missing values*/
	int     i,j,n,m;        /*count and iteration control variables*/
	double  value;          /*temp variable*/

	Nmiss = 0;
	for(i=1;i<=ntotal;i++)
	{
		n = i-Nmiss;
		m = 0;
		for (j=1;j<=2;j++)
		{
			fscanf(fp_in,"%lf",&value);
			if (value != MISSING)
			{
				if (j==1)     xxi[n]=value;
				if (j==2)     yyi[n]=value;
			}
			else
				m++;
		}
		if (m != 0)     Nmiss++;
		else if (xxi[n] < 0)      Nmiss++;
	}

	return Nmiss;
}

/****************************************************************
**  OUTPUT_DTMSnVCV--output variance and covariance matrix for a
*                    n-th degree polynomial model.  nboundparm is
*                         the number of elements in the bounded[] array
*****************************************************************/
void OUTPUT_DTMSnVCV(int nparm, int nboundparm, int Spec[], char *parmtxt[],
					 double **vcv, int *bounded)
{
	int i,j,checkvar, bound, counter, *notbounded;
	double  cor;

	checkvar = 0;
	/*output to screen and bmdswrk.002 temp file*/
	OUTPUT_TEXT("\n\n           Asymptotic Correlation Matrix of Parameter Estimates\n");


	fprintf (fp_out,"          ");

	notbounded = IVECTOR(1, nparm);
	bound = 0;
	counter = 0;

	for (i=1;i<=nboundparm;i++)
	{
		if(bounded[i]==0)
		{
			if(Spec[i]==0)
			{
				/*  printf ("%15s",parmtxt[i-1]); */
				fprintf (fp_out,"%13s",parmtxt[i-1]);
				counter++;
				notbounded[counter] = i;
			}
		}
		else if(Spec[i] == 0)
			bound++;
	}
	/*  printf ("\n"); */

	counter = 0;

	for (i=1;i<=nparm;i++)
	{
		counter++;
		/*  printf ("\n%10s",parmtxt[notbounded[counter]-1]); */
		fprintf (fp_out,"\n%10s",parmtxt[notbounded[counter]-1]);

		for (j=1;j<=nparm;j++)
		{
			if (vcv[i][i]==0.0 || vcv[j][j]==0.0)
				checkvar = 1;
			else
				checkvar = 0;

			if(checkvar == 0)
			{

				cor =  vcv[i][j]/(sqrt(fabs(vcv[i][i]))*sqrt(fabs(vcv[j][j])));
				fprintf (fp_out,
#ifndef RBMDS
					"%13.2g"
#else
					"%30.22g "
#endif
					,cor);
			}
			else
			{
				fprintf (fp_out," NA");
			}
		}
	}

	if(bound > 0)
	{
		fprintf(fp_out, "\n\nThe following parameter(s) have been estimated at a boundary\n");
		fprintf(fp_out, "point or have been specified.  Correlations are not computed:  \n\n");
	}

	for(i = 1; i <= nboundparm; i++)
	{
		if(bounded[i] == 1 && Spec[i] == 0)
		{
			fprintf(fp_out, "%s  ", parmtxt[i-1]);
		}
	}
	if(checkvar == 1)
	{
		fprintf(fp_out, "\n\nNA - This parameter's variance has been estimated at zero.\n");
	}

}
/**********************************************************
** READ_OBSDATAV--used to read Ni column data in Ni vectors.
***********************************************************/
int READ_OBSDATAV(int Nobs,double Xi[],int Ni[],double Ym[],double Yd[])

{
	int     Nmiss;          /*number of records with missing values*/
	int     i,j,n,m;        /*count and iteration control variables*/
	double  value;          /*temp variable*/
	double * yy;            /*response of each animal in a group*/
	yy = DVECTOR(1,Nobs+2);/*this line was added to define yy array*/

	Nmiss = 0;
	for(i=1;i<=Nobs;i++)
	{
		n = i-Nmiss;
		m = 0;
		for (j=1;j<=2;j++)
		{
			fscanf(fp_in,"%lf",&value);
			if (value != MISSING)
			{
				if (j==1)     Xi[n]=value;
				if (j==2)     Ni[n]=(int)value;
			}
			else
				m++;
		}
		for (j=3; j<=Ni[n]+2; j++)
		{
			fscanf(fp_in,"%lf",&value);
			if (value != MISSING)
			{
				yy[j-2]=value;
				Ym[j-2] += value/Ni[n];
			}
			else  m++;
		}
		for (j=3; j<=Ni[n]+2; j++)
		{
			fscanf(fp_in,"%lf",&value);
			if (value != MISSING)
				Yd[j-2] += (yy[j-2]-Ym[j-2])*(yy[j-2]-Ym[j-2])/(Ni[n]-1.0);
			else m++;
		}

		if (m != 0)                         Nmiss++;
		else if (Xi[n] < 0)      Nmiss++;
	}
	FREE_DVECTOR(yy,1, Nobs+2);   /*this line was added to free yy */
	return Nmiss;
}

/*******************************************************************
**  A3_Fit fits a likelihood to a "full" model of the form
**  Yij = Mu(i) + e(ij)  Var(eij) = k*Mu(i)^b.  The parameters Mu(i)
**  i = 1,2,...,Nobs, k, and p are estimated using the likelihood
**  maximization procedure, and the lkA3 = log-likelihood value
**  upon return.
*******************************************************************/

void AThree_Fit(int nparm, double p[], double gtol, int *iter, double *fret)
{

	int i;
	double *parms, *fitparms, *doses, ll, *svar, *means;
	double /**tmy, *t, **tmv,*/ *bsv;
	double **X, **XP, **XPX, *XPY, *Y;
	int *Spec2, *bind, *nanim, nresm, optite;
	int restr, nparms, nvar;

	switch (Spec[1] * 10 + Spec[2]) {
  case 11: /* Both alpha and rho fixed */
	  break;
  case 1: /* Estimating alpha; rho fixed */
	  /* estimate alpha as mean of Y/pow(mean, rho) */
	  p[1] = 0.0;
	  for (i = 1; i <= Nobs; i++)
		  p[1] += Yd[i] / pow(Ym[i], p[2]);
	  p[1] = log(p[1] / (double) Nobs);
	  break;
  case 10: /* Estimate rho, alpha fixed */
	  /* estimate rho as mean of (log(Var) - log(alpha)) / log(mu) */
	  p[2] = 0.0;
	  for (i = 1; i <= Nobs; i++)
		  p[2] += (Slog(Yd[i]) - Slog(p[1])) / Slog(Ym[i]);
	  p[2] = p[2] / (double) Nobs;
	  break;
  default:  /* estimating both variance parameters */
	  /* Allocate memory for linear regression */
	  X = DMATRIX(1, Nobs, 1, 2);
	  XP = DMATRIX(1, 2, 1, Nobs);
	  XPX = DMATRIX(1, 2, 1, 2);
	  XPY = DVECTOR(1, 2);
	  bsv = DVECTOR(1, 2);
	  Y = DVECTOR(1, Nobs);

	  /* set up regression problem; regress log(variance)
	  on log means */
	  for (i = 1; i <= Nobs; i++)
	  {
		  X[i][1] = 1.0;
		  X[i][2] = Slog(Ym[i]);
		  Y[i] = Slog(Yd[i]);
	  }
	  /* linear regression */
	  TRANSPOSE(X, XP, Nobs, 2);
	  MATMPYM2(XP, X, XPX, 2, Nobs, 2);
	  INVMAT(XPX, 2);
	  MATMPYV2(2, Nobs, XP, Y, XPY);
	  MATMPYV2(2, 2, XPX, XPY, bsv);
	  /* transfer the results to p[1:2] */
	  p[1] = bsv[1];
	  p[2] = bsv[2];
	  /* Free memory */
	  FREE_DMATRIX(X, 1, Nobs, 1, 2);
	  FREE_DMATRIX(XP, 1, 2, 1, Nobs);
	  FREE_DMATRIX(XPX, 1, 2, 1, 2);
	  FREE_DVECTOR(XPY, 1, 2);
	  FREE_DVECTOR(bsv, 1, 2);
	  FREE_DVECTOR(Y, 1, Nobs);
	  break;


	}

	for (i = 1; i <= nparm-2; i++)
	{
		p[i + 2] = Ym[i];
	}

	nvar = Nobs;
	restr = restrict;
	nparms = nparm;

	doses = DVECTOR(0, Nobs-1);
	means = DVECTOR(0, Nobs-1);
	svar =  DVECTOR(0, Nobs-1);
	nanim = IVECTOR(0, Nobs-1);
	parms = DVECTOR(0, nparm-1);
	fitparms = DVECTOR(0, nparm-1);
	Spec2 = IVECTOR(0, nparm-1);
	bind = IVECTOR(0, nparm-1);

	for (i = 1; i <= Nobs; i++)
	{
		nanim[i - 1] = (int) Ni[i];
		doses[i - 1] = Xi[i];
		means[i - 1] = Ym[i];
		svar[i - 1] = Yd[i];
	}

	for (i = 1; i <= nparm; i++)
	{
		Spec2[i - 1] = Spec[i];
		parms[i - 1] = p[i];
	}

	getmlea3_(&nvar, doses, means, nanim, svar, &nparms, parms,
		Spec2, parms, &restr, fitparms, &ll, &optite,
		&nresm, bind);

	for (i = 1; i <= nparm; i++)
	{
		p[i] = fitparms[i - 1];
	}

	*fret = ll;


	FREE_DVECTOR(doses, 0, Nobs-1);
	FREE_DVECTOR(means, 0, Nobs-1);
	FREE_DVECTOR(svar, 0, Nobs-1);
	FREE_IVECTOR(nanim, 0, Nobs-1);
	FREE_DVECTOR(parms, 0, nparm-1);
	FREE_DVECTOR(fitparms, 0, nparm-1);
	FREE_IVECTOR(Spec2, 0, nparm-1);
	FREE_IVECTOR(bind, 0, nparm-1);


}	/* end AThree_Fit */



void GetNewParms(double *p, int size)
/***********************************************************
Given a vector of parameter values, and the number of
parameters in that vector, this function will return three
new parameter values to restart the optimization if a "bad"
completion code is returned from GETCL(), using a uniform
random number centered at p[i]
***********************************************************/
{
	int i;

	/* Find parameters by randomly selecting new parameters in
	a uniform interval of p[i] +/- .005*p[i] */

	for(i = 0; i < size; i++){
		if(Spec[i+1] != 1)
			p[i] = (p[i]+p[i]*.005 - (p[i]-p[i]*.005))*(double)rand()/ (double)	RAND_MAX + p[i] - p[i]*.005;
	}

	/* If parameters are to be restricted, make sure restrictions
	are not violated */
	if(p[0] <= 0)
		p[0] = .00000001;

	for (i=1; i<size; i++){
		if ((p[i] < 0) && (restrict == 1))
			p[i] = -p[i];

		if ((p[i] > 0) && (restrict == -1))
			p[i] = -p[i];

	} /* end if */
}

void PolyMeans(int nobs, double p[], double Doses[], double means[])
/****************************************************
/   Given the number of dose levels, this function
/    calculates the mean at each dose level and stores
/    it in the array means[]
*****************************************************/
{
	int i,j;

	for(i = 1; i <= nobs; i++)
	{
		means[i] = p[nparm];
		for(j = nparm - 1; j >= 3; j--)
		{
			means[i] = means[i]*Doses[i] + p[j];
		}
	}
}

void F1iDoublePart(int nparm, int const_var, double p[], double **Fn1i,int obs)
/***********************************************************
Calculates the second partial derivatives of the function
F = ln(Vi) at dose[obs] where Vi is the estimated variance.  Fn1i[j][k]
contains the second partial with respect to parameters j and
k.
**********************************************************/
{
	double *mg, *vg, **mg2, **vg2, Vi;
	double meani, temp;
	int i, j, k;


	for(j = 1; j <= nparm; j++)
		for(k = 1; k <= nparm; k++)
			Fn1i[j][k] = 0.0;

	mg = DVECTOR(1,nparm);
	vg = DVECTOR(1,nparm);
	mg2 = DMATRIX(1,nparm,1,nparm);
	vg2 = DMATRIX(1,nparm,1,nparm);


	/* Compute the estimated mean */

	meani = p[nparm];
	for(i = nparm-1; i >= 3; i--)
	{
		meani = meani*Xi[obs] + p[i];
	}

	if(const_var == 1)
		Vi = p[1];
	else
		Vi = exp(p[1]  + log(fabs(meani)) * p[2]);
	if (Vi == 0) Vi = .000001;
	/* Get the partial derivatives of the mean function at dose[obs] */
	MeanPart(obs, p, mg);
	/* Get the partial derivatives of the variance function */
	VarPart(obs, const_var, Vi, meani, p, mg, vg);
	/* Get second partials of the mean function */
	Mean2Part(obs, p, mg2);
	/* Get second partials of the variance function */
	Var2Part(obs, const_var, Vi, meani, p, mg, mg2, vg2);


	/* Calculate partial derivative at dose[obs] */

	for(j = 1; j <= nparm; j++)
	{
		for(k = j; k <= nparm; k++)
		{
			temp = vg2[j][k] - vg[j]*vg[k]/Vi ;
			Fn1i[j][k] = temp/Vi;
			Fn1i[k][j] = Fn1i[j][k];
		}
	}


	FREE_DVECTOR(mg, 1, nparm);
	FREE_DVECTOR(vg, 1, nparm);
	FREE_DMATRIX(mg2,1,nparm,1,nparm);
	FREE_DMATRIX(vg2,1,nparm,1,nparm);

}

void F2iDoublePart(int nparm, int const_var, double p[], double **Fn2i,int obs)
/***********************************************************
/  Calculates the second partial derivatives os the function
/  F = 1/Vi at dose[obs] where Vi is the estimated variance.  Fn2i[j][k]
/  contains the second partial with respect to parameters j and
/  k.
**********************************************************/
{
	double *mg, *vg, **mg2, **vg2, Vi;
	double meani, temp;
	int i, j, k;


	for(j = 1; j <= nparm; j++)
		for(k = 1; k <= nparm; k++)
			Fn2i[j][k] = 0.0;

	mg = DVECTOR(1,nparm);
	vg = DVECTOR(1,nparm);
	mg2 = DMATRIX(1,nparm,1,nparm);
	vg2 = DMATRIX(1,nparm,1,nparm);


	/* Compute the estimated mean */

	meani = p[nparm];
	for(i = nparm-1; i >= 3; i--)
	{
		meani = meani*Xi[obs] + p[i];
	}

	if(const_var == 1)
		Vi = p[1];
	else
		Vi = exp(p[1] + log(fabs(meani)) * p[2]);
	if (Vi == 0) Vi = .000001;
	/* Get the partial derivatives of the mean function at dose[obs] */
	MeanPart(obs, p, mg);
	/* Get the partial derivatives of the variance function */
	VarPart(obs, const_var, Vi, meani, p, mg, vg);
	/* Get second partials of the mean function */
	Mean2Part(obs, p, mg2);
	/* Get second partials of the variance function */
	Var2Part(obs, const_var, Vi, meani, p, mg, mg2, vg2);


	/* Calculate partial derivative at dose[obs] */

	for(j = 1; j <= nparm; j++)
	{
		for(k = j; k <= nparm; k++)
		{
			temp = 2*vg[j]*vg[k]/Vi - vg2[j][k];
			Fn2i[j][k] = temp/(Vi*Vi);
			Fn2i[k][j] = Fn2i[j][k];
		}
	}

	FREE_DVECTOR(mg, 1, nparm);
	FREE_DVECTOR(vg, 1, nparm);
	FREE_DMATRIX(mg2,1,nparm,1,nparm);
	FREE_DMATRIX(vg2,1,nparm,1,nparm);

}

void F3iDoublePart(int nparm, int const_var, double p[], double **Fn3i,int obs)
/***********************************************************
/  Calculates the second partial derivatives of the function
/  F = (Ybar - Mi)**2/Vi at dose[obs] where Vi is the estimated variance
/  Ybar is the sample mean, and Mi is the estimated mean.
/  Fn1i[j][k]
/  contains the second partial with respect to parameters j and
/  k.
**********************************************************/
{
	double *mg, *vg, **mg2, **vg2, Vi;
	double Devi, meani, temp, temp2, temp3;
	int i, j, k;


	for(j = 1; j <= nparm; j++)
		for(k = 1; k <= nparm; k++)
			Fn3i[j][k] = 0.0;

	mg = DVECTOR(1,nparm);
	vg = DVECTOR(1,nparm);
	mg2 = DMATRIX(1,nparm,1,nparm);
	vg2 = DMATRIX(1,nparm,1,nparm);


	/* Compute the estimated mean */

	meani = p[nparm];
	for(i = nparm-1; i >= 3; i--)
	{
		meani = meani*Xi[obs] + p[i];
	}

	Devi = Ym[obs] - meani;

	if(const_var == 1)
		Vi = p[1];
	else
		Vi = exp(p[1] + log(fabs(meani)) * p[2]);
	if (Vi == 0) Vi = .000001;
	/* Get the partial derivatives of the mean function at dose[obs] */
	MeanPart(obs, p, mg);
	/* Get the partial derivatives of the variance function */
	VarPart(obs, const_var, Vi, meani, p, mg, vg);
	/* Get second partials of the mean function */
	Mean2Part(obs, p, mg2);
	/* Get second partials of the variance function */
	Var2Part(obs, const_var, Vi, meani, p, mg, mg2, vg2);

	for(j = 1; j <= nparm; j++)
	{
		for(k = j; k <= nparm; k++)
		{
			temp = 2*Vi*Vi*(mg[j]*mg[k] - Devi*mg2[j][k]);
			temp2 = 2*Devi*Vi*(vg[j]*mg[k] + mg[j]*vg[k]);
			temp3 = Devi*Devi*(2*vg[j]*vg[k] - Vi*vg2[j][k]);
			Fn3i[j][k] = (temp + temp2 + temp3)/(Vi*Vi*Vi);
			Fn3i[k][j] = Fn3i[j][k];
		}
	}

	FREE_DVECTOR(mg, 1, nparm);
	FREE_DVECTOR(vg, 1, nparm);
	FREE_DMATRIX(mg2,1,nparm,1,nparm);
	FREE_DMATRIX(vg2,1,nparm,1,nparm);
}

void MeanPart(int obs, double *p, double *mg)
/*******************************************************************
*  First partial derivatives of the mean function at observation obs,
*  contained in mg[] at exit.
******************************************************************/
{
	int i;

	mg[1] = 0.0; /* Variance parameters not in mean function */
	mg[2] = 0.0;
	mg[3] = 1.0;

	if(Xi[obs] != 0)
	{
		for(i = 4; i <= nparm; i++)
		{
			mg[i] = pow(Xi[obs], i-3);
		}
	}
	else
	{
		for(i = 4; i <= nparm; i++)
		{
			mg[i] = 0;
		}
	}
}

void VarPart(int obs, int const_var, double Vi, double meani, double *p,
			 double *mg, double *vg)
			 /*********************************************************************
			 First partial derivatives of the variance function.  Vi and meani
			 are the estimated mean and variance respectively, and should be
			 passed by the calling function.  const_var = 1 if variance is
			 constant, = 0 if not.  mg[] are the first partials of the mean
			 function and should be passed by the calling function.  Partials
			 stored in vg[] upon exit.
			 *********************************************************************/
{
	int j;

	if(const_var == 1) 
	{
		vg[1] = Vi/p[1];
		vg[2] = 0.0;
	}
	else
	{
		vg[1] = Vi;
		if (meani == 0)
			vg[2] = 0;
		else
			vg[2] = Vi*log(fabs(meani));
	} /* end else */
	for(j = 3; j <= nparm; j++)
		/* If constant variance, then
		p[2] = 0, and thus, vg[j] = 0 */
		vg[j] = 0.0;
} /* end varpart */

void Mean2Part(int obs, double *p, double **mg2)
/********************************************************************
*  Second partial derivates with respect to the mean function.
*  mg2[][] contains all second partials upon exit.
*******************************************************************/
{
	int i,j;

	for(i = 1; i <= nparm; i++)
	{
		for(j = 1; j <= nparm; j++)
		{
			mg2[i][j] = 0.0;
		}
	}

}



void Var2Part(int obs, int const_var, double Vi, double meani, double *p,
			  double *mg, double **mg2, double **vg2)
			  /*********************************************************************
			  *  Second partial derivatives of the variance function.  Vi and meani
			  *  are the estimated mean and variance respectively, and should be
			  *  passed by the calling function.  const_var = 1 if variance is
			  *  constant, = 0 if not.  mg[] are the first partials of the mean
			  *  function, mg2[][] are the second partials of the mean function.
			  *  Both should be passed by the calling function.  Partials
			  *  stored in vg2[][] upon exit.
			  *********************************************************************/
{
	double logam, abmn, temp;
	int j, k, Sign;

	abmn = fabs(meani);
	logam = Slog(abmn);

	if(const_var == 1)
	{
		/* constant variance.  Vi = alpha.  All second partials = 0 */

		for(j = 1; j <= nparm; j++)
			for(k = 1; k <= nparm; k++)
				vg2[j][k] = 0.0;
	}
	else
	{
		/* non constant variance.  Vi = exp(lalpha + log(|Meani|)*rho) */

		if(meani < 0)
			Sign = -1;
		else
			Sign = 1;

		vg2[1][1] = Vi;
		if (meani == 0)
			vg2[1][2] = 0.0;
		else
			vg2[1][2] = Vi*logam;
		vg2[2][1] = vg2[1][2];
		for(j = 3; j <= nparm; j++)
		{
			vg2[1][j] = Sign*p[2]*Vi*mg[j]/abmn;
			vg2[j][1] = vg2[1][j];

		}
		if (meani == 0)
			vg2[2][2] = 0.0;
		else
			vg2[2][2] = Vi*logam*logam;
		for(j = 3; j <= nparm; j++)
		{
			vg2[2][j] = 0.0;
			vg2[2][j] = Sign*Vi*mg[j]*(p[2]*logam + 1)/abmn;
			vg2[j][2] = vg2[2][j];
		}

		for(j = 3; j <= nparm; j++)
		{
			for(k = j; k <= nparm; k++)
			{
				temp = (p[2]-1)*mg[j]*mg[k]/abmn + Sign*mg2[j][k];
				vg2[j][k] = p[2]*Vi*temp/abmn;
				vg2[k][j] = vg2[j][k];
			}
		}
	}
}

void Get_BMRS(double *p, double Dose, double bmdl1, double *BMRVals, int sign, int
			  bmr_type)
			  /*********************************************************************
			  *
			  *  Finds four BMR values for computation of BMDL curve.  bmdl1 is
			  *  an input variable. It is the user BMDL for the BMR specified
			  *  by the user. bmdl2 is a second bmdl value corresponding to the
			  *  BMR in BMRVals[2] (both assigned prior to the call to this
			  *  function.  BMRVals[1] should also contain
			  *  the user input BMR value that corresponds to bmdl2
			  *  (as a Point BMR type). when this function is called.  After
			  *  completion, BMRVals[1]...BMRVals[6] will contain five BMR values
			  *  (Point BMR values), including
			  *  the user input BMR in BMRVals[1] and the "next" BMR value in
			  *  BMRVals[2].
			  *
			  **********************************************************************/
{
	double range, mn;
	double ymax, ymin, bmr_percentage;
	double Percent_Change[4] = {.1, .225, .375, .5};
	int i;

	ymax = Ym[1];
	ymin = Ym[1];

	/* Find the maximum and minimum observed mean to get BMR's in terms
	of a percentage of the range of observed data */

	for(i = 2; i <= Nobs; i++)
	{
		if(Ym[i] > ymax)
			ymax = Ym[i];
		if(Ym[i] < ymin)
			ymin = Ym[i];
	}

	range = ymax - ymin;          /* The range of means */
	/* bmrtype as a percentage of the range */
	bmr_percentage = (BMRVals[1]-ymin)/range;


	/* If adverse direction is "down", then make this as a percentage
	of the range, but the pecent decrease over the maximum observed
	mean as opposed to the percent increase over the minimum */

	if(sign == -1)
		bmr_percentage = 1 - bmr_percentage;

	if(sign == 1)
		mn = ymin;
	else
		mn = ymax;

	/* Find a set of  Point type BMRs that are spaced throughout
	the range, and such that they are not too crowded by the already
	computed BMDL */

	if(.05 < bmr_percentage && bmr_percentage < .10)
	{
		Percent_Change[0] = bmr_percentage + .05;
	}
	else
	{
		if(.2 < bmr_percentage && bmr_percentage < .25)
		{
			Percent_Change[1] = bmr_percentage + .05;
		}
		else
		{
			if(.35 < bmr_percentage && bmr_percentage < .4)
			{
				Percent_Change[2] = bmr_percentage + .05;
			}
			else
			{
				if(.475 < bmr_percentage && bmr_percentage < .525)
				{
					Percent_Change[3] = bmr_percentage + .05;
				}
			}
		}
	}

	/* Get BMR values and give them back as the BMR type specified by
	the user */

	for(i = 0; i <= 3; i++)
	{
		BMRVals[i+2] = mn + sign*Percent_Change[i]*range;

		if(bmr_type == 0)   /* Absolute deviation */
		{
			BMRVals[i+2] = fabs(BMRVals[i+2] - p[3]);
		}
		else
		{
			if(bmr_type == 1)        /* Std Dev */
			{
				BMRVals[i+2] = fabs(BMRVals[i+2] - p[3]);
				BMRVals[i+2] = BMRVals[i+2]/(sqrt(exp(p[1] + log(p[3]) * p[2])));
			}
			else
			{
				if(bmr_type == 2)        /* Relative dev */
				{
					BMRVals[i+2] = fabs(BMRVals[i+2] - p[3])/p[3];
				}
				else
				{
					if(bmr_type == 4)   /* Extra risk */
					{
						BMRVals[i+2] = (BMRVals[i+2] - p[3])/p[4];
					}
					else      /* Point deviation */
						BMRVals[i+2] = BMRVals[i+2];
				}
			}
		}
	}

}

/***************************************************
*    OUTPUT_BENCHMD2--output specified benchmark dose.
****************************************************/
void OUTPUT_BENCHMD2(int pdcol, double BMD)
{


	OUTPUT_TEXT(" \n\n             Benchmark Dose Computation");

	/* output to the screen and to bmdswrk.002 temp file */

	/*  printf("Specified effect =%14.6g\n\n", bmdparm.effect); */
	fprintf(fp_out, "\nSpecified effect =%14.6g\n\n", bmdparm.effect);

	if (bmr_type == 0)
	{
		/*  printf("Risk Type        =     Absolute deviation \n\n");*/
		fprintf(fp_out,"Risk Type        =     Absolute deviation \n\n");
	}
	else if (bmr_type == 1)
	{
		/*  printf("Risk Type        =     Estimated standard deviations from the control mean \n\n"); */
		fprintf(fp_out,"Risk Type        =     Estimated standard deviations from the control mean\n\n");
	}
	else if (bmr_type == 2)
	{
		/*  printf("Risk Type        =     Relative deviation \n\n"); */
		fprintf(fp_out,"Risk Type        =     Relative deviation \n\n");
	}
	else if (bmr_type == 3)
	{
		/*  printf("Risk Type        =     Point estimate \n\n"); */
		fprintf(fp_out,"Risk Type        =     Point deviation \n\n");
	}
	else
	{
		/*  printf("Risk Type        =     Extra risk \n\n"); */
		fprintf(fp_out,"Risk Type        =     Extra risk \n\n");
	}

	/*  printf("Confidence level =%14.6g\n\n",bmdparm.level); */
	/*    printf("             BMD =%14.6g\n\n",BMD); */

	fprintf(fp_out, "Confidence level =%14.6g\n\n",bmdparm.level);
	fprintf(fp_out, 
#ifndef RBMDS
		"             BMD = %14.6g\n\n"
#else
		"             BMD = %30.22g\n\n"
#endif
		,BMD);

}    /* end OUTPUT_BENCHMD2 */

void GetMoreParms(double *p, int size)
/***********************************************************
Given a vector of parameter values, and the number of
parameters in that vector, this function will return three
new parameter values to restart the optimization if a "bad"
completion code is returned from GETCL(), using a uniform
random number centered at p[i]
***********************************************************/
{
	int i;

	/* Find parameters by randomly selecting new parameters in
	a uniform interval of (-2, 2) */

	for(i = 0; i < size; i++)
	{
		if(Spec[i+1] != 1)
			p[i] = 2*(double)rand()/ (double) RAND_MAX;
	}
	/* If parameters are to be restricted, make sure restrictions
	are not violated */
	if(p[0] <= 0)
		p[0] = .01;

	for (i=1; i<size; i++)
	{
		if ((p[i] < 0) && (restrict == 1))
			p[i] = -p[i];
		if ((p[i] > 0) && (restrict == -1))
			p[i] = -p[i];
	} /* end if */

}

void GetOtherParms(double *p, int size)
/***********************************************************
Given a vector of parameter values, and the number of
parameters in that vector, this function will return three
new parameter values to restart the optimization if a "bad"
completion code is returned from GETCL(), using a uniform
random number centered at p[i]
***********************************************************/
{
	int i;

	/* Find parameters by randomly selecting new parameters in
	a interval of p[i] +(-2, 2) */

	for(i = 0; i < size; i++)
	{
		if(Spec[i+1] != 1)
			p[i] = p[i] + 2*(double)rand()/ (double) RAND_MAX;
	}
	/*If parameters are to be restricted, make sure restrictions
	are not violated */
	if(p[0] <= 0)
		p[0] = .1;

	for (i=1; i<size; i++)
	{
		if ((p[i] < 0) && (restrict == 1))
			p[i] = -p[i];
		if ((p[i] > 0) && (restrict == -1))
			p[i] = -p[i];
	} /* end if */

}


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

