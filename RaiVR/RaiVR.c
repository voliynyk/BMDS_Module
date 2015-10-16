/****************************************************************
*
* IMPORTANT NOTE:  The following variable is the version number for
*                  the current model.  THIS MUST BE CHANGED as
*				   important changes are made to the models.
*
*****************************************************************/
char Version_no[]="RaiVR Model. (Version: 2.12;  Date: 04/27/2015)";

/****************************************************************
*
* RaiVR.C - a ANSI C program for RaiVR model fitting with/without
*             a natural background rate in Benchmark Dose.
*
* Date: March 25, 1997
*
********************************************************************
* Modification Log:
*
* Version Number: 2.3
* Modified By: Micheal Ferree
* Modified Date: 8/13/2005
* Reason:
*
* Version Number: 2.4
* Modified By: Woodrow Setzer
* Modified Date: 10/26/2005
* Reason:  1) Free all allocated memory
*          2) RaiVR_probs was returning without calculating any probs
*             for some input parameters when calculating BMDLs
*          3) Drop calculation of likelihoods for "full" and "reduced"
*             models; the reduced model likelihood is useless, and the
*             full model likelihood is bogus.
*          4) Fix an error in Assist/N_Goodness.c that hosed the
*             goodness-of-fit table.
*          5) Upper bound for theta2 was set to 0; should be (and now is)
*             Max_double.
*          6) Added conditional compilation flags for RBMDS and logging.
*          7) Removed references to and definition of Binary_root()
* Outstanding Problems: 1) standard errors for parameters reported as "1".
*                       2) BMDL calculation seems unstable
*
* Version Number: 2.5
* Modified By: Geoffrey
* Date: 1/12/2007
* Reason: Incremented version number.
*		 Added last parameter "0" (don't print SE) in OUTPUT_DTMS3PARMS().
*
* Version Number: 2.6
* Modified By: Woodrow Setzer
* Date: 2/20/2007
* Reason: Incremented version number to reflect changed compilation options.
*
* Version Number: 2.7
* Modified By: G. Nonato
* Modification Date: 04/10/2008
* Reason: (Per BMDS 2.0: Problem Report 157 & 147)
*       Fix the Observation # < parameter # for RaiVR model problem.
*       Added code to free-up allocated memories before exiting thru ERRORPRT()
*
* Version Number: 2.8
* Modified By: G. Nonato
* Modification Date: 10/28/2009
* Reason:
*      To be able to process files/folders with spaces (PR 257)
*      Fix program freeze due to long variable names (PR 278)
*      Process long files up to 256 characters (PR 303 and 308)
*      Modify code for easy maintenance, all lengths for file name,
*        model name, and column names are placed in benchmark.h
*
* Version Number: 2.9
* Modified By: Louis Olszyk
* Modification Date: 02/28/2013
* Reason: PR 444 - Fix wording in plot titles
*
* Version Number: 2.10
* Modified By: Cody Simmons
* Modification Date: 09/08/2014
* Reason: PR232 and PR244
*       Added the bootstrap method for goodness of fit calculation
*
* Version Number: 2.11
* Modified By: Cody Simmons
* Modification Date: 09/15/2014
* Reason: PR300
*       Created scaled residual of interest table for nested models.
*
* Version Number: 2.12
* Modified By: Cody Simmons
* Modification Date: 4/27/2015
* Reason: PR545
*       Display the min, max and mean of absolute value of scaled residuals.
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

void RaiVR_fit(int nparm, int ngrp, double p[], double gtol,
			   int *iter, double *fret);
void RaiVR_BMD (int nparm, double p[], double gtol, int *iter, double xlk,
				double Rlevel[], double Bmdl[], double *BMD);
void RaiVR_grad(int nparm, int Spec[], double ptf[], double grad[]);
int RaiVR_vcv(int nparm, int Spec[], double p[], double **vcv);
void Which_Bounded (int Spec[], double Parms[], int bounded[]);
void Predict(double doses[], double Lsc[], int nobs, double Parms[], double P[]);
extern void initialparms(int, double, double, double [], double [],
						 double [], double);
void TEMP_ANOVA_OUTPUT(char *anatxt[], AnaList anasum[]);

#define EPS 3.0e-8
#define TOLX (10*EPS)
#define STPMX1 10.0
#define ALF 0.000001
#define FREEALL FREE_DVECTOR (xi,1,n); FREE_DVECTOR (pnew, 1,n);\
	FREE_DVECTOR(hdg,1,n); FREE_DVECTOR(g,1,n);\
	FREE_DVECTOR(dg,1,n); \
	FREE_DMATRIX(hessin,1,n,1,n);

#define float double
#define GOLD 1.618034
#define GLIMIT 100
#define TINY 1.0e-20
#define SWAP(a,b, junk) (junk)=(a); (a)=(b); (b)=(junk);
#define SHFT(a,b,c,d,junk) (junk)=(a); (a)=(b); (b)=(c); (c)=(d); (d)=(junk);
#define ZEPS 1.0e-8
#define MOV3(a,b,c, d,e,f)(a)=(d);(b)=(e);(c)=(f);
#define TOL 2.0e-6

/*** Define input and output files's name  *********************/
char     fin[FLENGTH];  /*input temp file*/
char    fout[FLENGTH];  /*output temp file*/
char    fout2[FLENGTH];
char	plotfilename[FLENGTH];	/* file to pass to GnuPlot */
char    *Parm_name[]={"alpha", "beta", "theta1", "theta2","rho",
"phi1", "phi2", "phi3", "phi4", "phi5",
"phi6", "phi7", "phi8", "phi9", "phi10"};
char    *anatxt[]={"Full model", "Fitted model", "Reduced model"};
char    fname2[FLENGTH], logfile[FLENGTH], *dot2;
#ifdef LOGGING_ON
FILE    *fp_log;            /*  log file if requested */
#endif
int     silent = true;     /*  switch for log file   */


typedef enum {
	dummy, alpha, beta, theta1, theta2, rho, phi1, phi2, phi3,
	phi4, phi5, phi6, phi7, phi8, phi9, phi10} Model_Parms;

/*** variables will not be changed execpt Spec  *******/
int    *Spec;    /* vector used to identify user input parm. */
int    *IniSp;
double *Yp;      /* positive dependent variable data array */
double *Yn;      /* negative dependent variable data array */
double *Xi;      /* independent variable data array */
double *ScXi;    /* dose scaled by the max dose (for optimizing) */
double *Ypp;     /* predicted positive depdendent values */
double *Ep;      /* estimated probability */
double *Ls;      /* litter specific covariate */
int    *Xg;      /* dose group index */
double *SR; 			/*scaled residual */
double *Rlevel;
double *Bmdl;
double *IniP;
int    Nobs, nparm, ngrp, restrict, initial, appendix, smooth, fixedSize, bmdlCurve;
double maxdose;  /* max dose level */
double xmax, xmin;
double sijfixed;
double Rel_Conv, Parm_Conv, Maxloglik;
double PowerUpperBound = 18.0;
int    MxLkCnt = 0;   /* counts the number of calls to MAX_lk */
/** changing variable **/
int    replace;
double tD, BMD_lk, LR, BMR;
int    DeBuG = 0;
FILE   *fp_dbg;
int    ErrorFlag; /* Error States from DMNGB */

int BSIter;      /*number of iterations for bootstrapping method */     
long BSSeed;     /*seed values for bootstrapping method - 0 defaults to time-generated*/

/*  int tmpcnt = 0; */

/* GLN - 04/10/2008
*  Free-up allocated memory before exit
*  upon encountering fatal error.
*/
void FreeUp_mem(double *Parms, VarList *varsum, AnaList *anasum, double  **vcv, double *GXi, double *GYp, double *GYn, int *bounded)
{
	FREE_DVECTOR (Parms, 1, nparm);
	FREE_IVECTOR (bounded, 1, nparm);
	FREE_DVECTOR (IniP, 1, nparm);
	FREE_DVECTOR (Xi, 1, Nobs);
	FREE_DVECTOR (ScXi, 1, Nobs);
	FREE_DVECTOR (Yp, 1, Nobs);
	FREE_DVECTOR (Yn, 1, Nobs);
	FREE_DVECTOR (Ls, 1, Nobs);
	FREE_IVECTOR (Xg, 1, Nobs);
	FREE_IVECTOR (IniSp, 1, nparm);
	FREE_IVECTOR (Spec, 1, nparm);
	FREE_VLVECTOR(varsum, 1, 3);
	FREE_ALVECTOR(anasum, 1, 3);
	FREE_DVECTOR (Rlevel, 1, 5);
	FREE_DVECTOR (Bmdl, 1, 5);
	FREE_DMATRIX (vcv, 1, nparm, 1, nparm);
	if(GXi)
		FREE_DVECTOR (GXi, 1, ngrp);
	if(GYp)
		FREE_DVECTOR (GYp, 1, ngrp);
	if(GYn)
		FREE_DVECTOR (GYn, 1, ngrp);

	if (fp_log != (FILE *) NULL)
		fclose(fp_log);

	return;
}

/****************************************************************
* main--main function used to call RaiVR mode fitting program.
*****************************************************************/
int main (int argc, char *argv[])
{
	int     iter, i, junk;             /* iteration variable */
	int     bmdose;                       /* flag for computing benchmark dose */
	int     Nmiss;                        /* number of records with missing values */
	int     nparm_known;                  /* number of specified parameters */
	int     *bounded;
	int     nvar;
	double  BMD, lkf, lkr, xlk, W, junk1; /* log likelihoods */
	double  *Parms;                       /* parameter array */
	VarList *varsum;                      /* info for variables--p. dep.,n. dep., indep. */
	AnaList *anasum;                      /* information for ANONA analysis */
	double  **vcv;                        /* variance and covariance matrix */
	double  back, back1;
	double  *GYp, *GYn, *GXi;
	char    model_name[MNLENGTH], user_note[UNLENGTH], junkname[FLENGTH];
	char    dose_name[CNLENGTH], posi_name[CNLENGTH], nega_name[CNLENGTH], junkname1[FLENGTH], junkname2[FLENGTH];
	time_t  ltime;
	char long_path_name[FLENGTH];

	//struct _timeb tstruct;

	/* Set time zone from TZ environment variable. If TZ is not set,
	* the operating system is queried to obtain the default value
	* for the variable.
	*/
	/* _tzset(); */
	time( &ltime );

	//Min_increment = DBL_EPSILON;
	//Max_double    = DBL_MAX;

	/********************************************************************
	* {QH 2004/01/14 PR# }
	* Added to show version number if executed from command line with -v
	*********************************************************************/
	if(argc == 2)
		show_version(argv[1], Version_no);

	if(argc < 2)
	{
		fprintf(stderr, "ERROR:  Requires two arguments\nUsage:  %s <file.(d)>\n", argv[0]);
		exit(1);

	} /* end if */

	if (argc > 2)
	{
		path_name2(argc, argv, long_path_name);
		argv[1] = long_path_name;
	}

	fp_in=fopen(argv[1], "r");
	if (fp_in==NULL)
	{
		/*printf("Error in opening input  file.\n");
		printf ("...now exiting to system...\n");*/

		fprintf(stderr,"Error in opening input file.\n");
		fprintf (stderr,"...Exited to system!\n");
		exit (1);
	}

	/* open the log file if silent = 0 */
#ifdef LOGGING_ON
	strcpy(logfile,argv[1]);
	dot2 = strchr(logfile, (int) '.');
	(*dot2) = (char) 0;
	strcpy(fname2,logfile);
	strcat(logfile,"-rai.log");
	fp_log = fopen(logfile, "w");
	if (fp_log == (FILE *) NULL)
	{
		ERRORPRT("Unable to open log for RaiVR.c");
	}
#endif

	fscanf(fp_in, "%s", model_name);
	fscanf(fp_in, "%[ ^\n]", user_note);
	fscanf(fp_in, "%[^\n]", user_note);
	fscanf(fp_in, "%s", junkname);
	fscanf(fp_in, "%s", junkname);
	fscanf(fp_in, "%d%d",&Nobs, &ngrp);

	/*assign number of parameters*/
	nparm = 5+ngrp;

#ifdef LOGGING_ON
	fprintf(fp_log,"model_name: %s\n\n",model_name);
	fprintf(fp_log,"INPUT VALUES FOR DATA SET\n\n");
	fprintf(fp_log,"Nobs = %d                  (Number of Litters Observed)\n",Nobs);
	fprintf(fp_log,"ngrp = %d                    (Number of Dose Groups)\n",ngrp);
	fprintf(fp_log,"nparm = %d                  (The Number of Parameters in the Model)\n",nparm);
	fflush(fp_log);
#endif
	/*allocate memory for arrays*/
	Parms   = DVECTOR(1, nparm);
	IniP    = DVECTOR(1, nparm);
	IniSp   = IVECTOR(1, nparm);
	Spec    = IVECTOR(1, nparm);
	bounded = IVECTOR(1, nparm);
	Xi      = DVECTOR(1, Nobs);
	ScXi    = DVECTOR(1, Nobs);
	Yp      = DVECTOR(1, Nobs);
	Yn      = DVECTOR(1, Nobs);
	Ls      = DVECTOR(1, Nobs);
	Xg      = IVECTOR(1, Nobs);
        SR      = DVECTOR(1, Nobs);
	varsum  = VLVECTOR(1, 3);
	anasum  = ALVECTOR(1, 3);
	Rlevel  = DVECTOR(1, 5);
	Bmdl    = DVECTOR(1, 5);
	vcv     = DMATRIX (1,nparm,1,nparm);


	junk = 0;    /* Used to see if an extension was added to output file name */

	/* get filenames */
	Get_Names(argv[1], fout, fout2, plotfilename);

	//************ input All the data from  .(d) file
	fscanf(fp_in,"%d%lf%lf%d%d%d%d%d%d",&ITMAX, &Rel_Conv, &Parm_Conv, &bmdlCurve, &restrict, &bmdose,
		&fixedSize,&appendix, &smooth);
	if(appendix==Yes)
	{
		fp_out=fopen(fout,"a");
	}
	else
	{
		fp_out=fopen(fout,"w");
	}
#ifndef RBMDS
	fp_out2=fopen(fout2,"w");
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
		exit (1);
	}

	fscanf(fp_in,"%lf%d%lf%d%ld", &bmdparm.effect, &bmdparm.risk, &bmdparm.level, &BSIter, &BSSeed);

        if (BSSeed == 0)              /* Set seed from time clock if default(BSSeed=0) is specified */
          {
          BSSeed = time (NULL);
          }

#ifdef LOGGING_ON
	fprintf(fp_log,"ITMAX = %d                 (Maximum number of iterations)\n",ITMAX);
	fprintf(fp_log,"Rel_Conv = %g     (Rel Function Convergence, default=1e-8)\n",Rel_Conv);
	fprintf(fp_log,"Parm_Conv = %g    (Parameter Convergence)\n",Parm_Conv);
	fprintf(fp_log,"bmdlCurve = %d               (BMDL Curve Calculation; 1=yes, 0=no)\n",bmdlCurve);
	fprintf(fp_log,"restrict = %d                (Restrict Rho>=1; 1=yes, 0=no)\n",restrict);
	fprintf(fp_log,"bmdose = %d                  (BMD Calculation; 1=yes, 0=no)\n",bmdose);
	fprintf(fp_log,"appendix = %d                (Append or Overwrite output file)\n",appendix);
	fprintf(fp_log,"smooth = %d                  (Smooth Option)\n",smooth);
	fprintf(fp_log,"bmdparm.effect = %4g       (BMR Factor)\n",bmdparm.effect);
	fprintf(fp_log,"bmdparm.risk = %d            (Risk Type; 1=Added, 0=Extra)\n",bmdparm.risk);
	fprintf(fp_log,"bmdparm.level = %4g        (Confidence Level)\n",bmdparm.level);
#endif

	/* Print model and file information on output page */

	Output_Header(Version_no, argv[1], plotfilename, ctime(&ltime), user_note);

	if (bmdose < 0 || bmdose > 1)
	{
		ERRORPRT("Error in choosing benchmark dose computation.");
	}



	/*obtain user input parameters*/

	READ_PARAMETERS(nparm,Parms);
	SHFT(Parms[2], Parms[3], Parms[4], Parms[5], junk1);
	FILL_SPECVECTOR(nparm,Parms,Spec);
	nparm_known = COUNT_SPECVECTOR(nparm,Spec);

	/*obtain user input initial parameters values      *******/

	fscanf(fp_in,"%d", &initial);
	READ_PARAMETERS(nparm,IniP);
	SHFT(IniP[2], IniP[3], IniP[4], IniP[5], junk1);
	FILL_SPECVECTOR(nparm,IniP,IniSp);

	for(i = 1; i <= nparm; i++)
	{
		if(Spec[i] == Yes)
		{
			IniP[i] = 1;
		}
	}

#ifdef LOGGING_ON
	fprintf(fp_log,"\nThe specified or default parameters are:\n");
	fprintf(fp_log,"                                            Spec Value\n");
	fprintf(fp_log,"Parms[ 1] = %12.6g        (alpha)         %d\n",Parms[1],Spec[1]);
	fprintf(fp_log,"Parms[ 2] = %12.6g        (beta)          %d\n",Parms[2],Spec[2]);
	fprintf(fp_log,"Parms[ 3] = %12.6g        (theta1)        %d\n",Parms[3],Spec[3]);
	fprintf(fp_log,"Parms[ 4] = %12.6g        (theta2)        %d\n",Parms[4],Spec[4]);
	fprintf(fp_log,"Parms[ 5] = %12.6g        (rho)           %d\n",Parms[5],Spec[5]);
	for(i=6; i<=nparm; i++)
	{
		fprintf(fp_log,"Parms[%2d] = %12.6g        (Phi%d)          %d\n",i,Parms[i],i-5,Spec[i]);
	}
	fprintf(fp_log,"\ninitial = %d                (Initialized Parameters; 1=yes, 0=no)\n",initial);
	fprintf(fp_log,"The initial parameter values are:\n");
	fprintf(fp_log,"                                           IniSp Value\n");
	fprintf(fp_log,"IniP[ 1] =  %12.6g        (alpha)         %d\n",IniP[1],IniSp[1]);
	fprintf(fp_log,"IniP[ 2] =  %12.6g        (beta)          %d\n",IniP[2],IniSp[2]);
	fprintf(fp_log,"IniP[ 3] =  %12.6g        (theta1)        %d\n",IniP[3],IniSp[3]);
	fprintf(fp_log,"IniP[ 4] =  %12.6g        (theta2)        %d\n",IniP[4],IniSp[4]);
	fprintf(fp_log,"IniP[ 5] =  %12.6g        (rho)           %d\n",IniP[5],IniSp[5]);
	for(i=6; i<=nparm; i++)
	{
		fprintf(fp_log,"IniP[%2d] =  %12.6g        (Phi%d)          %d\n",i,IniP[i],i-5,IniSp[i]);
	}
	fprintf(fp_log,"The Spec and IniSp vectors are just flags to tell if the user gave\n");
	fprintf(fp_log,"a value for the specific parameter.\n");
#endif

	/*obtain observation data into Yp, Yn, Xi, Ls, Xg vectors*/
	fscanf(fp_in,"%s%s%s%s%s", dose_name, posi_name, nega_name, junkname1, junkname2);
	Nmiss = READ_OBSDATA5V(Nobs,Xi,Yp,Yn,Ls,Xg);

	Sort_4_By_Dose(Nobs, Xi, Yn, Yp, Ls);   /* Sort the arrays by dose and move around */
	/* the other arrays appropriately */

	/* Create scaled vector of doses */
	maxdose = Xi[Nobs];
	for(i = 1; i <= Nobs; i++)
	{
		ScXi[i] = Xi[i]/maxdose;
	}


	Xg[1] = 1;
	for(i=1; i<Nobs; i++)				         /* Hopefully someday, this "group" variable */
	{			                                 /* will not have to be entered by the user */
		if(Xi[i+1] == Xi[i])				 /* This little loop gets the "group" values */
		{                                                /* without explicitly entering them at the */
			Xg[i+1] = Xg[i];				 /* spreadsheet level */
		}
		else
		{
			Xg[i+1] = Xg[i] + 1;
		}
	}

	ngrp = Xg[Nobs];

	GXi = DVECTOR(1, ngrp);
	GYp = DVECTOR(1, ngrp);
	GYn = DVECTOR(1, ngrp);


	Nobs -= Nmiss;             /* extern variable Nobs has been changed */

#ifdef LOGGING_ON
	fprintf(fp_log,"\n  Dose             Positive            Negative         Litter Specific         Dose");
	fprintf(fp_log,"\n  Level           Respondents         Respondents          Covariate            Group\n");
	for(i=1; i<=Nobs; i++)
		fprintf(fp_log,"Xi[%3d]=%4g      Yp[%3d]=%3g         Yn[%3d]=%3g         Ls[%3d]=%3g         Xg[%3d]=%3d\n",
		i,Xi[i],i,Yp[i],i,Yn[i],i,Ls[i],i,Xg[i]);
	fprintf(fp_log,"\nNmiss = %3d             (Number of missing values)\n",Nmiss);
	fprintf(fp_log,"Nobs  = %3d             (Number of observations)\n",Nobs);
	fflush(fp_log);
#endif

	nparm_known = COUNT_SPECVECTOR(nparm, Spec);

	if (Nobs < (nparm-nparm_known))
	{
		FreeUp_mem(Parms, varsum, anasum, vcv, GXi, GYp, GYn, bounded);
		ERRORPRT("Observation # < parameter # for RaiVR model.");
	}
	//*********** end of input data.

	/*output title and summary of intput data  ****************************/
	OUTPUT_TEXT("\n The probability function is: ");
	OUTPUT_TEXT("\n\n Prob. = [1-exp(-Alpha-Beta*Dose^Rho)]*exp(-(Th1+Th2*Dose)*Rij),");
	OUTPUT_TEXT("\n           where Rij is the litter specific covariate.");
	/*indicates restrict power*/
	if (restrict==Yes)
	{
		OUTPUT_TEXT("\n Restrict Power rho >= 1. ");
	}

	if (Spec[(int) theta1]==1 && Parms[(int) theta1]<=0.0)
	{
		if (Spec[(int) theta2]==0 || Parms[(int) theta2]<0 )
		{
			Warning ("\n Warning: Theta2 has to be fixed at zero since Theta1=0.");
		}
		Spec[(int) theta2]=1;
		Parms[(int) theta2]=0;
	}

	fprintf (fp_out, "\n\n\n Total number of observations = %d",Nobs+Nmiss);
	fprintf (fp_out, "\n Total number of records with missing values = %d",Nmiss);
	fprintf (fp_out, "\n Total number of parameters in model = %d", nparm);
	fprintf (fp_out, "\n Total number of specified parameters = %d\n\n", nparm_known);

	fprintf(fp_out, "\n Maximum number of iterations = %d\n", ITMAX);
	fprintf(fp_out, " Relative Function Convergence has been set to: %g\n", Rel_Conv);
	fprintf(fp_out, " Parameter Convergence has been set to: %g\n\n", Parm_Conv);
        fprintf (fp_out, " Number of Bootstrap Iterations per run: %d\n", BSIter);
        fprintf (fp_out, " Bootstrap Seed:  %ld\n\n", BSSeed);


	if(nparm_known > 0)
	{
		fprintf (fp_out, " User specifies the following parameters:");
		for (i=1; i<=nparm; i++)
		{
			if(Spec[i] == Yes)
			{
				fprintf (fp_out, "\n %15s = %10.5g", Parm_name[i-1], Parms[i]);
			}
		}
		fprintf (fp_out, "\n\n");
	}

	fflush(fp_out);

	/*compute init_lkf for full model and init_lkr for reduced model*/
	//  COMPUTE_DTMS3LOGLKH(Nobs,&lkf,&lkr,varsum,Yp,Yn);

	lkf = 0.0;
	varsum[1].S = 0;
	varsum[2].S = 0;
	for (i=1; i<=Nobs; i++)
	{
		varsum[1].S += Yp[i];
		varsum[2].S += Yn[i];
		W = Yp[i] / (Yp[i]+Yn[i]);
		if (W > 0)
		{
			lkf += Yp[i] * log(W);
		}
		if (W < 1)
		{
			lkf += Yn[i] * log(1 - W);
		}
	}
	W = varsum[1].S / (varsum[1].S + varsum[2].S);
	lkr = varsum[1].S * log(W) + varsum[2].S * log(1- W);

#ifdef LOGGING_ON
	fprintf(fp_log,"\nlkf = %4g          (Likelihood for full model)\n", lkf);
	fprintf(fp_log,"xlk = %4g          (Likelihood for fitted model)\n", xlk);
	fprintf(fp_log,"lkr = %4g          (Likelihood for reduced model)\n", lkr);
#endif

	/*fitting RaiVR model and output parameter estimators */

#ifdef LOGGING_ON
	fprintf(fp_log, "\n************Before call to RaiVR_fit*****************\n");
	fprintf(fp_log, "\nnparm = %2d;  ngrp = %2d;  iter = %4d;  \nxlk = %6g \n\n", nparm, ngrp, iter, xlk);
	fprintf(fp_log,"Parms[ 1] = %12.6g  (alpha)\n",Parms[1]);
	fprintf(fp_log,"Parms[ 2] = %12.6g  (beta)\n",Parms[2]);
	fprintf(fp_log,"Parms[ 3] = %12.6g  (theta1)\n",Parms[3]);
	fprintf(fp_log,"Parms[ 4] = %12.6g  (theta2)\n",Parms[4]);
	fprintf(fp_log,"Parms[ 5] = %12.6g  (rho)\n",Parms[5]);
	for(i=6; i<=nparm; i++)
	{
		fprintf(fp_log,"Parms[%2d] = %12.6g  (phi%d)\n",i,Parms[i],i-5);
	}
	fprintf(fp_log,"******************************************************\n");
	fflush(fp_log);
#endif

	RaiVR_fit (nparm, ngrp, Parms, EPS, &iter, &xlk);

#ifdef LOGGING_ON
	fprintf(fp_log,"\n************After call to RaiVR_fit*********************\n");
	fprintf(fp_log, "\nnparm = %2d;  ngrp = %2d;  iter = %4d; \nxlk = %6g \n\n", nparm, ngrp, iter, xlk);
	fprintf(fp_log,"Parms[ 1] = %12.6g  (alpha)\n",Parms[1]);
	fprintf(fp_log,"Parms[ 2] = %12.6g  (beta)\n",Parms[2]);
	fprintf(fp_log,"Parms[ 3] = %12.6g  (theta1)\n",Parms[3]);
	fprintf(fp_log,"Parms[ 4] = %12.6g  (theta2)\n",Parms[4]);
	fprintf(fp_log,"Parms[ 5] = %12.6g  (rho)\n",Parms[5]);
	for(i=6; i<=nparm; i++)
	{
		fprintf(fp_log,"Parms[%2d] = %12.6g  (phi%d)\n",i,Parms[i],i-5);
	}
	fprintf(fp_log,"*********************************************************\n");
	fflush(fp_log);
#endif

	/* Which parameters are on their bounds? */
	Which_Bounded (Spec, Parms, bounded);

#ifdef LOGGING_ON
	fprintf(fp_log,"\n************After call to Which_Bounded*********************\n");
	fprintf(fp_log,"Parms[ 1] = %12.6g  (alpha)\n",Parms[1]);
	fprintf(fp_log,"Parms[ 2] = %12.6g  (beta)\n",Parms[2]);
	fprintf(fp_log,"Parms[ 3] = %12.6g  (theta1)\n",Parms[3]);
	fprintf(fp_log,"Parms[ 4] = %12.6g  (theta2)\n",Parms[4]);
	fprintf(fp_log,"Parms[ 5] = %12.6g  (rho)\n",Parms[5]);
	for(i=6; i<=nparm; i++)
	{
		fprintf(fp_log,"Parms[%2d] = %12.6g  (phi%d)\n",i,Parms[i],i-5);
	}
	fprintf(fp_log,"*********************************************************\n");
	fflush(fp_log);
#endif


	/* If We didn't get an acceptable fit, don't compute a bmd and */
	/* don't compute standard errors */
	if (ErrorFlag != 0)
	{
		bmdose = No;
		for (i=1; i<=nparm; i++)
		{
			bounded[i] = 1;
		}
	}
	/* Some things to do if we DID get a fit */
	if (ErrorFlag == 0)
	{
		/* Compute the approximate covariance matrix using the inverse of the */
		/* hessian matrix*/
		/* Initialize the memory */
		INITIALIZE_DMATRIX (vcv, nparm, nparm);
		/* vcv is now an nparm x nparm matrix of 0's */
		/* compute the hessian */
		nvar = RaiVR_vcv (nparm, Spec, Parms, vcv);
		fflush(fp_out);
		/* the first nvar x nvar rows and columns of vcv now contain the */
		/* hessian of the loglikelihood wrt the unfixed parameters */

#ifdef LOGGING_ON

		for (i = 1; i <= nparm; i++)
			fprintf(fp_log, "bounded[%2d] = %2d\n", i, bounded[i]);
		fprintf(fp_log, "***************************************************\n");
		fprintf(fp_log, "                  VCV MATRIX\n   ");
		for (i = 1; i <= nparm; i++)
			fprintf(fp_log,"        (%2d)", i);
		fprintf(fp_log, "\n");
		for (i = 1; i <= nparm; i++)
		{
			fprintf(fp_log, "(%2d) " , i);
			for (j = 1; j <= nparm; j++)
				fprintf(fp_log, "%10.3g  ", vcv[i][j]);
			fprintf(fp_log, "\n");
		}
#endif

		/* Remove the rows and columns of vcv that correspond to bounded */
		/* parameters */
		nvar = Take_Out_Bounded_Parms(nparm, bounded, vcv);

		/* Do the inverse */
		INVMAT (vcv, nvar);
	}


	/*output Parameter estimates and standard errors, if not all the */
	/*parameters were fixed */
	if (ErrorFlag != -1)
	{
		OUTPUT_DTMS3PARMS (nparm, Spec, bounded, Parms, Parm_name, vcv, 1);
		fflush(fp_out);
	}

	/* Now compute and output the analysis of deviance */
	DTMS3ANOVA (nparm, Nobs, Spec, lkf, xlk, lkr, anasum, bounded);
	TEMP_ANOVA_OUTPUT (anatxt, anasum);	/*output ANOVA table */
	/* OUTPUT_DTMS3ANOVA (anatxt, anasum);*/	/*output ANOVA table */
	fflush(fp_out);

	/*print a goodness of fit table if we fit something*/
	if (ErrorFlag == 0)
	{
		N_Goodness (ngrp, nparm, Parms, bounded, Nobs, Xi, Yp, Yn,
                            Ls, Xg, SR);
		fflush(fp_out);
	}

	/* output to **.002 ************/
#ifndef RBMDS
	fprintf (fp_out2, "\n BMD_flag \t %d \n Nosb \t%d \n nparm \t%d",  bmdose, ngrp, 6 );
	fprintf (fp_out2, "\n  Con_lev \t%3.3g ", bmdparm.level);
	fprintf (fp_out2, "\n  RiskType \t%d ", bmdparm.risk);
	fprintf (fp_out2, "\n  Effect \t%3.3g ", bmdparm.effect);
	for (i=1; i<=5; i++)
	{
		fprintf (fp_out2, "\n %s \t %5.3g", Parm_name[i-1], Parms[i]);
	}
	fprintf (fp_out2, "\n fixedSize \t%f ", sijfixed);
	{
		/* Calculation of 95% confidence intervals at each dose level for */
		/* graphical output */
		double *LL, *UL, *phat;

		LL = DVECTOR(1, ngrp);
		UL = DVECTOR(1, ngrp);
		phat = DVECTOR(1, ngrp);

		Nested_CI(ngrp, Nobs, Yp, Yn, Xg, 0.95, LL, phat, UL);

		/* Gets dose level at each dose group */
		GXi[1] =Xi[1];
		for (i=1; i<Nobs; i++)
		{
			if (Xg[i+1] != Xg[i])
			{
				GXi[ Xg[i+1] ] = Xi[i+1];
			}
		}

		fprintf (fp_out2, "\n\n Data");
		for (i = 1; i <= ngrp; i++)
		{
			fprintf (fp_out2, "\n %f %f %f %f", GXi[i], phat[i], LL[i], UL[i]);
		}
		FREE_DVECTOR(LL,1,ngrp);
		FREE_DVECTOR(UL,1,ngrp);
		FREE_DVECTOR(phat,1,ngrp);
	} /* End of CI computation and output */

	fprintf (fp_out2,"\n Max_Min_dose \n  %f %f \n", xmax, xmin);
	fflush(fp_out2);
#endif

	/*compute benchmark dose*/
	if (bmdose==Yes)
	{
		if(Spec[(int) beta]==Yes)
		{
#ifndef RBMDS
			fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  No);
#endif
			fprintf (fp_out,"\n\n  %s parameter is fixed. The likelihood function can not be reparameterized in BMD.",
				Parm_name[1]);
			exit(0);
		}

		back=(1-exp(-Parms[(int) alpha]) )*exp(-Parms[(int) theta1]*sijfixed);
		back1=1-back;
		if (bmdparm.risk==ADDED)
		{
			back1=1;
		}
		RaiVR_BMD (nparm, Parms, EPS, &junk, xlk, Rlevel, Bmdl, &BMD);
                SRoI(ngrp, Nobs, GXi, Xi, Xg, SR, Ls, sijfixed, BMD);
#ifndef RBMDS
		fprintf (fp_out2, "\n  RSL \t%f",bmdparm.effect*back1+back);
		fprintf (fp_out2, "\n  BMD \t%f",BMD);
		fprintf (fp_out2, "\n  BMDL \t%f",Bmdl[1]);

		fprintf (fp_out2,"\n\n BMD_line");
		fprintf (fp_out2,"\n %f %f", (xmin-xmax/100), bmdparm.effect*back1+back );
		fprintf (fp_out2,"\n %f %f", BMD, bmdparm.effect*back1+back );
		fprintf (fp_out2,"\n %f %f", BMD, -0.1);

		fprintf (fp_out2,"\n\n BMDL_line");
		fprintf (fp_out2,"\n %f %f", Bmdl[1], -0.1);
		fprintf (fp_out2,"\n %f %f", Bmdl[1], bmdparm.effect*back1+back );

		fprintf (fp_out2, "\n\n BMDL_Curve_flag \t %d  \n smooth_opt  %d", bmdlCurve, smooth);
		fprintf (fp_out2,"\n\n BMDL_curve");
		fprintf (fp_out2,"\n 0.00000 %f", back);
		for (i=1; i<=5; i++)
		{
			fprintf (fp_out2,"\n %f %f", Bmdl[i], Rlevel[i]*back1+back);
		}
#endif
	}

        N_Bootstrap (ngrp, nparm, Parms, bounded, Nobs, Xi, Yp, Yn,
                            Ls, Xg, SR, BSIter, BSSeed);
        if (bmdose==Yes)
		{
                	if(fixedSize == Yes)
			{
				fprintf(fp_out, "\n\n\nTo calculate the BMD and BMDL, the litter specific covariate is fixed at\n");
				fprintf(fp_out,     "the mean litter specific covariate of control group: %f", sijfixed);
			}
			else
			{
				fprintf(fp_out, "\n\n\nTo calculate the BMD and BMDL, the litter specific covariate is fixed at\n");
				fprintf(fp_out,     "the mean of litter specific covariate: %f", sijfixed);
			}

		OUTPUT_BENCHMD(1,BMD);
#ifndef RBMDS
		fprintf(fp_out, "            BMDL = %14.6g\n\n", Bmdl[1]);
#else
		fprintf(fp_out, "            BMDL = %30.22g\n\n", Bmdl[1]);

#endif
                fflush(fp_out);
		}

#ifndef RBMDS
	fprintf (fp_out2,"\n\n Check_result %d", Yes); //indicate all requried data have been output.
#endif

#ifdef LOGGING_ON
	fflush(fp_log);
	fclose (fp_log);
#endif

	FREE_DVECTOR (Parms,1,nparm);
	FREE_DVECTOR (IniP,1,nparm);
	FREE_DVECTOR (Xi,1,Nobs);
	FREE_DVECTOR (ScXi,1,Nobs);
	FREE_DVECTOR (Yp,1,Nobs);
	FREE_DVECTOR (Yn,1,Nobs);
	FREE_DVECTOR (Ls,1,Nobs);
	FREE_IVECTOR (Xg,1,Nobs);
	FREE_IVECTOR (IniSp,1,nparm);
	FREE_IVECTOR (Spec,1,nparm);
	FREE_IVECTOR (bounded,1,nparm);
	FREE_VLVECTOR(varsum,1,3);
	FREE_ALVECTOR(anasum,1,3);
	FREE_DVECTOR (Rlevel,1,5);
	FREE_DVECTOR (Bmdl, 1, 5);
	FREE_DMATRIX (vcv,1,nparm,1,nparm);
	FREE_DVECTOR (GXi, 1, ngrp);
	FREE_DVECTOR (GYp, 1, ngrp);
	FREE_DVECTOR (GYn, 1, ngrp);

	/*close opened temp files*/
	CLOSE_FILES ();
	return(0);
} /*end of main*/

/*****************************************************************************/
/* Which_Bounded -- Fills the 1-based vector bounded with 1 if the           */
/*                  corresponding parameter is on one of its boundaries or   */
/*                  was fixed, 0 otherwise.                                  */
/*                                                                           */
/*                  Global:                                                  */
/*                    nparm                                                  */
/*                                                                           */
/*                  input:                                                   */
/*                    Spec  1-based vector giving which parameters are fixed */
/*                    Parms 1-based vector giving the full vector of         */
/*                          parameter values                                 */
/*                  output:                                                  */
/*                    bounded 1-based vector whose elements are 1 if         */
/*                            the corresponding parameter is either          */
/*                            fixed                                          */
/*                            (i.e., Spec[i] == 1) or its value is on        */
/*                            a boundary                                     */
/*****************************************************************************/

void Which_Bounded (int Spec[], double Parms[], int bounded[])
{
	int i;
	double maxdose;

	for (i=1; i<=nparm; i++)
	{
		bounded[i] = Spec[i];
	} /* end for */

	maxdose = Xi[1];

	for(i=1; i<=ngrp; i++)
	{
		if(Xi[i] > maxdose)
		{
			maxdose = Xi[i];
		}
	}

	/*    There are many parameter restrictions in the RaiVR model.  The if statements */
	/*    below address the parameter space restrictions, and adjust the degrees of */
	/*    freedom appropriately */

	/*    Theta1 + Theta2*dose >= 0: */

	if((Spec[(int) theta1] == No || Spec[(int) theta2] == No) && VERYCLOSE(Parms[(int) theta1], -Parms[(int) theta2]*maxdose))
	{
		if (Spec[(int) theta1] == Yes) bounded[(int) theta2] = 1;
		if (Spec[(int) theta2] == Yes) bounded[(int) theta1] = 1;
	}


	/*  Rho >= 1 when specified by the user: */

	if(restrict == Yes && Spec[(int) rho] == No && VERYCLOSE(Parms[(int) rho], 1))
	{
		bounded[(int) rho] = 1;
	}

	/*  Rho >= 0: */

	if(Spec[(int) rho] == No && VERYCLOSE(Parms[(int) rho], 0))
	{
		bounded[(int) rho] = 1;
	}

	/*  alpha > = 0: */

	if(Spec[(int) alpha] == No && VERYCLOSE(Parms[(int) alpha], 0))
	{
		bounded[(int) alpha] = 1;
	}

	/*  beta > = 0: */

	if(Spec[(int) beta] == No && VERYCLOSE(Parms[(int) beta], 0))
	{
		bounded[(int) beta] = 1;
	}

	/*  Phi1,..., Phig >= 0: */

	for(i=6; i<=nparm; i++)
	{
		if(Spec[i] == No && VERYCLOSE(Parms[i], 0))
		{
			bounded[i] = 1;
		}
	}
}/* End of Which_Bounded */

/**************************************************************************/
/* Model_DF -- Decrements the number of model parameters by the number of */
/*             parameters prespecified or "stuck" on a boundary of        */
/*             parameter space.                                           */
/*                                                                        */
/*             Global:                                                    */
/*               nparm number of parameters                               */
/*                                                                        */
/*             input:                                                     */
/*               bounded 1-based vector st bounded[i] means that the ith  */
/*                       parameter is either fixed or bounded             */
/*                                                                        */
/*             returns:                                                   */
/*               Model degrees of freedom - # fixed - #stuck or fixed     */
/*                                                                        */
/*             NOTE: Subtracting the number "stuck" on a boundary is      */
/*                   incorrect, and needs to be revisited.                */
/**************************************************************************/

int Model_DF(int bounded[])
{
	int i, df;

	df = nparm;
	for(i=1; i<=nparm; i++)
	{
		df -=bounded[i];
	}
	return df;
}/* End of Model_DF */

/***********************************************************
* Predict -- returns predicted values for given parameter
*            values.
*            For use OUTSIDE of MAX_lk, that is, with
*            UNTRANSFORMED model parameters.
*
*            input:
*                   doses is the vector of doses (1 - based)
*                   Lsc   is the 1-based vector of litter-specific covariates
*                   nobs is the number of observations (1 - based)
*                   Parms is the vector of UNTRANSFORMED parameters (1 - based)
*            output:
*                   P is the vector (of length nobs) of predicted
*                     values
***********************************************************/
void Predict(double doses[], double Lsc[], int nobs, double Parms[],
	double P[])
{
	int i;

	for (i=1; i<=nobs; i++)
	{
		if (doses[i] <= 0.0)
		{
			P[i] = (1.0 - exp(-Parms[(int) alpha])) * exp(-Parms[(int) theta1] * Lsc[i]);
		}
		else
		{
			P[i] = (1.0 - exp(-(Parms[(int) alpha] + Parms[(int) beta] * pow(doses[i], Parms[(int) rho])))) *
				exp(-(Parms[(int) theta1] + Parms[(int) theta2] * doses[i]) * Lsc[i]);
		}
	}
}/* End of Predict */

/*******************************************************************
**RaiVR_probs -- computes litter-specific probabilities of an
*                  adverse response, and, optionally, the gradient
of the probabilities wrt the parameters.
This uses the transformed parameters.
* Inputs: nobs: number of observations (litters)
*         doses: vector of length nobs containing the administered dose
*                for each litter.
*         sij: vector of length nobs containing the litter-specific
*              covariate for each litter.
*         nparm: number of parameters
*         compgrad: compute the gradient (0 = No, 1 = Yes)
* Outputs: probs: vector of length nobs containing the probability of
*                 adverse response for each litter (based at 1)
gradij: matrix that is nobs x 5 of the derivative
of the probability wrt each parameter, for each litter.
Rows correspond to litters, columns correspond to
parameters.
*
* Global variables used: replace, Spec, bmdparm, BMR, tD
*********************************************************************/
void RaiVR_probs (double probs[], double p[], int compgrad, double **gradij)
{
	int i, j;
	double  *pint, ex, ex1, ex2, ex4, pwx;
	double  P0=0, PR=0, Del, lambda=0, C1, C2, C3, C4, S=0;

	pint = DVECTOR(1, nparm);
	for (i = 1; i <= 5; i++)
	{
		pint[i] = p[i];
		for (j = 1; j <= Nobs; j++)
			gradij[j][i] = 0;
	}

	if (replace==Yes)  /* replace beta */
	{
		S=sijfixed;
		P0=(1-exp(-pint[(int) alpha]))*exp(-S*pint[(int) theta1]);               //Prob.(Dose=0)
		lambda=exp( pint[(int) theta1]*(1+pint[(int) theta2]*tD)*S );
		PR = BMR + P0;                                /* Added risk */
		if (bmdparm.risk==EXTRA)
		{
			/* PR=BMR+(1-BMR)*P0; */      /* extra risk: this looks wrong */
			PR = P0 + BMR * (1 - P0);
		}
		if (1-lambda*PR<=0) /* if this is true, the estimate of beta is infinite, so fitted probs are all 1.0.  Return that. */
		{
			pint[(int) beta] = Max_double;
		}
		else
		{
			pint[(int) beta] = ( -log(1 - lambda*PR) - pint[(int) alpha] ) / pow(tD, pint[(int) rho]);
		}
	}

	for (i=1; i<=Nobs; i++)
	{
		if (ScXi[i] > 0)
		{
			pwx = pow(ScXi[i], pint[(int) rho]);
			ex1 = exp( -pint[(int) alpha] - pint[(int) beta]*pwx );
		}
		else
		{
			pwx = 0.0;
			ex1 = exp( -pint[(int) alpha] );
		}
		if (pint[(int) theta1] <=0 )
		{
			ex2=1;
		}
		else
		{
			ex2 = exp( -pint[(int) theta1]*(1+pint[(int) theta2]*ScXi[i])*Ls[i] );
		}
		ex = (1.0 - ex1)*ex2;

		PROBABILITY_INRANGE(&ex);

		probs[i] = ex;

		if (compgrad == 1)  /* Compute partial derivatives of log likelihood function wrt parameters */
		{

			if (replace==Yes)  /* Must compute partial derivatives of beta wrt other parameters */
			{
				Del = lambda/(1-lambda*PR);
				ex4 = exp(-pint[(int) alpha]-pint[(int) theta1]*S);
				if (bmdparm.risk==EXTRA)
				{
					C1 = 1-Del*ex4*(1-BMR);                                  /* Extra Risk */
					C3 = -Del*S*( (1 + pint[(int) theta2]*tD)*PR - (1-BMR)*P0);
				}
				else
				{
					C1 = 1-Del*ex4;                                          /* Added Risk */
					C3 = -Del*S*( (1 + pint[(int) theta2]*tD)*PR - P0);
				}
				C2 = -pow(tD, pint[(int) rho]);
				//	  C3 = -Del*S*BMR*(1+p[(int) theta2]*tD);
				C4 = -Del*S*tD*PR*pint[(int) theta1];

				gradij[i][2] = ex1*ex2*pwx;      //it is Xi[i](i.e pwx), not tD.
				gradij[i][1] = ( ex1*ex2 + gradij[i][2]*(C1/C2) );
				gradij[i][3] = ( (ex1-1)*ex2*Ls[i]*(1+pint[(int) theta2]*ScXi[i]) + gradij[i][2]*(C3/C2) );
				gradij[i][4] = ( (ex1-1)*ex2*Ls[i]*ScXi[i]*pint[(int) theta1] + gradij[i][2]*(C4/C2) );
				gradij[i][5] = -gradij[i][2]*pint[(int) beta]*log(tD);

				if (ScXi[i] > 0.0)
				{
					gradij[i][5] += ex1*ex2*pint[(int) beta]*log(ScXi[i])*pwx;
				}

			}
			else  /* Beta is not defined in terms of the other parameters */
			{
				gradij[i][1] = ex1*ex2;
				gradij[i][2] = ex1*ex2*pwx;
				gradij[i][3] = (ex1-1)*ex2*Ls[i]*(1+pint[(int) theta2]*ScXi[i]);
				gradij[i][4] = (ex1-1)*ex2*Ls[i]*ScXi[i]*pint[(int) theta1];
				if (ScXi[i] <= 0.0)
				{
					gradij[i][5]=0.0;
				}
				else
				{
					gradij[i][5] = ex1*ex2*pint[(int) beta]*log(ScXi[i])*pwx;
				}

			}

		}
	}

	FREE_DVECTOR(pint, 1, nparm);
}

/*******************************************************************
**RaiVR_lk -- used to compute the log likelihood for RaiVR model.
* 		     Extern var.: Nobs, Xi, Yp, Yn, Ls, Xg.
*
*********************************************************************/
void RaiVR_lk(long int *nvar, double *x, long int *nf, double *f,
	long int *uiparm, double *urparm, void (*ufparm)())
{
	double  xlk;          		  // log likelihood.
	int     i, j, k, plus5, jfixed, jvar;
	double  tm1, tm2, tm3, tm;   	          //temp var.
	double  *p;	                          /* for "untransform" parms. */
	double  *probs;                         /* litter-specific probabilities */
	double  **gradij;                       /* not used right here, but needed for function call */
	int     compgrad;                       /* if compgrad=1 RaiVR_probs computes gradient of the litter specific */
	/* probabilities in addition to litter specific probabilities */
#ifdef LOGGING_ON
	int     run=-9999;                      /* number of the run that we wish output in the log file */

	if ((MxLkCnt == run))
	{
		fprintf(fp_log, "\nFunction: RaiVR_lk (Beginning)\n");
		fprintf(fp_log, "The fixed parameters are:\n");
		for (i = 0; i <= (nparm-*nvar-1); i++)
			fprintf(fp_log, "urparm[%2d] = %12.6g   %d\n", i, urparm[i], &urparm[i]);
		fprintf(fp_log, "The varying parameters are:\n");
		for (i = 0; i <= (*nvar-1); i++)
			fprintf(fp_log, "x[%2d] = %12.6g   %d\n", i, x[i], &x[i]);
	}
#endif
	p     = DVECTOR(1, nparm);
	probs = DVECTOR(1, Nobs);
	gradij = DMATRIX(1, Nobs, 1, 5);

	jfixed = jvar = 0;
	for(j=1; j<=nparm; j++)  //reconstruct the parameter vector.
	{
		if(Spec[j]==Yes)
		{
			p[j] = urparm[jfixed];
			jfixed++;
		}
		else
		{
			p[j] = x[jvar];
			jvar++;
		}
	}

	compgrad = 0;

#ifdef LOGGING_ON
	if ((MxLkCnt == run))
	{
		fprintf(fp_log, "\nFunction: RaiVR_lk (Before RaiVR_probs)\n");
		fprintf(fp_log, "The p and Spec array are:\n");
		for (i = 1; i <= nparm; i++)
			fprintf(fp_log, "p[%2d] = %12.6g  %d      Spec[%2d] = %2d  %d\n",
			i, p[i], &p[i], i, Spec[i], &Spec[i]);
		fprintf(fp_log, "The fixed parameters are:\n");
		for (i = 0; i <= (nparm-*nvar-1); i++)
			fprintf(fp_log, "urparm[%2d] = %12.6g   %d\n", i, urparm[i], &urparm[i]);
		fprintf(fp_log, "The varying parameters are:\n");
		for (i = 0; i <= (*nvar-1); i++)
			fprintf(fp_log, "x[%2d] = %12.6g   %d\n", i, x[i], &x[i]);
	}
#endif
	RaiVR_probs(probs, p, compgrad, gradij);
#ifdef LOGGING_ON
	if ((MxLkCnt == run))
	{
		fprintf(fp_log, "\nFunction: RaiVR_lk (After RaiVR_probs)\n");
		fprintf(fp_log, "The fixed parameters are:\n");
		for (i = 0; i <= (nparm-*nvar-1); i++)
			fprintf(fp_log, "urparm[%2d] = %12.6g   %d\n", i, urparm[i], &urparm[i]);
		fprintf(fp_log, "The varying parameters are:\n");
		for (i = 0; i <= (*nvar-1); i++)
			fprintf(fp_log, "x[%2d] = %12.6g   %d\n", i, x[i], &x[i]);
	}
#endif
	xlk = 0.0;
	for (i=1; i<=Nobs; i++)
	{
		tm1 = 0.0;
		tm2 = 0.0;
		tm3 = 0.0;
		plus5=5+Xg[i];
		j = (int)Yp[i];

		if (probs[i] == 0.0 && j > 0)
		{
			tm1 -= 40.0; /* was 400 */
		}
		else
		{
			for (k = 1; k <= j; k++)
			{
				tm = probs[i] + (k - 1) * p[plus5];
				tm1 += log (tm);
			}
		}
		j = (int) Yn[i];
		if (probs[i] >= 1.0 && j > 0)
		{
			tm2 -= 40.0; /* was 400 */
		}
		else
		{
			for (k = 1; k <= j; k++)
			{
				tm = 1.0 - probs[i] + (k - 1) * p[plus5];
				tm2 += log (tm);
			}
		}
		j = (int) (Yn[i] + Yp[i]);
		for (k = 1; k <= j; k++)
		{
			tm = 1.0 + (k - 1) * p[plus5];
			tm3 += log (tm);
		}
		xlk += (tm1 + tm2 - tm3);
	}

	FREE_DVECTOR(p, 1, nparm);
	FREE_DVECTOR(probs, 1, Nobs);
	FREE_DMATRIX (gradij, 1, Nobs, 1, 5);
	*f = -xlk;
#ifdef LOGGING_ON
	if ((MxLkCnt == run))
	{
		fprintf(fp_log, "MxLkCnt = %3d     (Counts Calls to MAX_lk)\n", MxLkCnt);
		fprintf(fp_log, "\nFunction: RaiVR_lk (Exiting)\n");
		fprintf(fp_log, "The fixed parameters are:\n");
		for (i = 0; i <= (nparm-*nvar-1); i++)
			fprintf(fp_log, "urparm[%2d] = %12.6g   %d\n", i, urparm[i], &urparm[i]);
		fprintf(fp_log, "The varying parameters are:\n");
		for (i = 0; i <= (*nvar-1); i++)
			fprintf(fp_log, "x[%2d] = %12.6g   %d\n", i, x[i], &x[i]);
		fprintf(fp_log, "The likelihood value for these is: %12.6g\n", -xlk);
	}
#endif
}/* End of RaiVR_lk */

/*******************************************************************
**RaiVR_g -- used to computer the gradients for RaiVR model.
*		Extern var.: Nobs, Xi, Yp, Yn, Ls, Xg.
*
**********************************************************************/
void RaiVR_g(long int *nvar, double *x, long int *nf, double *g,
	long int *uiparm, double *urparm, void (*ufparm)())
{
	double  tm1, tm2, tm3, tm1a, tm2a, tm3a, tm, tm12;   //temp var.
	double  *dd, *p, *tmp_g, **gradij, *probs;
	int     i, j, k, plus5, jfixed, jvar, compgrad;
#ifdef LOGGING_ON
	int     run=-9999;                          /* number of the run that we wish output in the log file */

	if ((MxLkCnt == run))
	{
		fprintf(fp_log, "\nFunction: RaiVR_g (Beginning)\n");
		fprintf(fp_log, "The fixed parameters are:\n");
		for (i = 0; i <= (nparm-*nvar-1); i++)
			fprintf(fp_log, "urparm[%2d] = %12.6g\n", i, urparm[i]);
		fprintf(fp_log, "The varying parameters are:\n");
		for (i = 0; i <= (*nvar-1); i++)
			fprintf(fp_log, "x[%2d] = %12.6g\n", i, x[i]);
	}
#endif
	dd = DVECTOR(1, nparm);
	p  = DVECTOR(1, nparm);
	probs = DVECTOR(1, Nobs);
	gradij = DMATRIX(1, Nobs, 1, 5);
	tmp_g = DVECTOR(1, nparm);    /* this variable was create to ensure that g wasn't */
	/* overwritten since g is passed from DMNGB with a  */
	/* length of nvar which is usually less than nparm  */

	jfixed=jvar=0;
	for(j=1; j<=nparm; j++)  //reconstruct the parameter vector.
	{
		if(Spec[j]==Yes)
		{
			p[j] = urparm[jfixed];
			jfixed++;
		}
		else
		{
			p[j] = x[jvar];
			jvar++;
		}
	}

	compgrad = 1;
	RaiVR_probs (probs, p, compgrad, gradij);

	/** initial tmp_g[j]'s            **************/
	for (j=1; j<=nparm; j++)
	{
		tmp_g[j]=0.0;
	}
	for (i=1; i<=Nobs; i++)
	{
		/*Compute first partial derivatives*/

		tm1 = 0.0;
		tm2 = 0.0;
		tm3 = 0.0;
		tm1a = 0.0;
		tm2a = 0.0;
		tm3a = 0.0;

		for(j=6; j<=nparm; j++)
		{
			dd[j]=0.0;
		}
		plus5 = 5+Xg[i];

		j = (int) Yp[i];
		if (probs[i] > 0.0)
		{
			for (k=1 ; k<=j; k++)
			{
				tm = probs[i] + (k - 1) * p[plus5];
				tm1 += 1.0 / tm;
				tm1a += (1.0 / tm) * (k-1);
			}
		}

		j = (int)Yn[i];
		if (probs[i] < 1.0)
		{
			for (k = 1; k <= j; k++)
			{
				tm = 1.0 - probs[i] + (k - 1) * p[plus5];
				tm2 += 1.0 / tm;
				tm2a += (1.0 / tm) * (k - 1);
			}
		}

		j = (int) (Yp[i]+Yn[i]);
		for (k = 1; k <= j; k++)
		{
			tm = 1.0 + (k - 1) * p[plus5];
			if (tm == 0.0)
			{
				tm = 0.000001;
			}
			tm3 += 1.0 / tm;
			tm3a += (1.0 / tm) * (k - 1);
		}

		tm12 = (tm1 - tm2);
		for (j=1; j<=5; j++)
		{
			dd[j] = gradij[i][j] * tm12;
		}

		dd[plus5] = (tm1a + tm2a - tm3a);

		for (j=1; j<=nparm; j++)
		{
			tmp_g[j] -= dd[j];

		}
	}
	/* end of first partial deri. */

	jvar=0;
	for(j=1; j<=nparm; j++)  //reconstruct the parameter vector.
	{
		if(Spec[j]==No)
		{
			g[jvar] = tmp_g[j];
			jvar++;
		}
	}
#ifdef LOGGING_ON
	if ((MxLkCnt == run))
	{
		fprintf(fp_log, "\nFunction: RaiVR_g (Exiting)\n");
		fprintf(fp_log, "The fixed parameters are:\n");
		for (i = 0; i <= (nparm-*nvar-1); i++)
			fprintf(fp_log, "urparm[%2d] = %12.6g\n", i, urparm[i]);
		fprintf(fp_log, "The varying parameters are:\n");
		for (i = 0; i <= (*nvar-1); i++)
			fprintf(fp_log, "x[%2d] = %12.6g\n", i, x[i]);
		fprintf(fp_log, "The gradients are:\n");
		for (i = 0; i <= (*nvar-1); i++)
			fprintf(fp_log, "g[%2d] = %12.6g\n", i, g[i]);
	}
#endif
	FREE_DVECTOR(dd, 1, nparm);
	FREE_DVECTOR(p, 1, nparm);
	FREE_DVECTOR(probs, 1, Nobs);
	FREE_DMATRIX (gradij, 1, Nobs, 1, 5);
	FREE_DVECTOR(tmp_g, 1, nparm);
	return;
}/* End of RaiVR_g */

/****************************************************************************
**
* RaiVR_grad -- Computes the gradient of the nested Rai Van Ryzin likelihood
* function with respect to the user form of the parameters.  This is to
* be used in RaiVR_vcv, to compute a finite difference approximation to
* the hessian of the likelihood function
* Input: nparm -- the number of parameters, the dimension of all the following
arrays.
Spec[] -- if the ith parameter is fixed by the user, Spec[i] == 1,
otherwise Spec[i] == 0.
ptf[] -- vector of parameters, in user form (external form),
based at 1.
index -- 1 for first derivatives
2 for second derivatives
* Output: grad[] -- the gradient of the loglikelihood function
(N.B. RaiVR_g returns the gradient of -loglikelihood)
Based at 1.

*****************************************************************************/
void RaiVR_grad(int nparm, int Spec[], double ptf[], double grad[])
{

	int i, j, k, plus5;
	double *p, *dd;
	double x, tm, tm1, tm1a, tm2, tm2a, tm3, tm3a, pwx, ex, ex1, ex2, ex3, tm12; 
	double *outgrad;

	p     = DVECTOR(1, nparm);
	dd = DVECTOR(1, nparm);
	//prob  = DVECTOR(1, Nobs);
	//pgrad = DMATRIX(1, Nobs, 1, 5);
	outgrad  = DVECTOR(1, nparm);

	for (i = 1; i <= nparm; i++)
	{
		outgrad[i] = 0;
		p[i] = ptf[i];
	}

	/* ------- Transform the parameters to the "internal" form -------    */
  	for (j = 6; j <= nparm; j++)
    		p[j] = p[j] / (1 - p[j]);	// Phi --> Psi.

	for (i = 1; i <= Nobs; i++)
	{
		// reinitialize variables
		tm1 = 0.0;
		tm2 = 0.0;
		tm3 = 0.0;
		tm1a = 0.0;
		tm2a = 0.0;
		tm3a = 0.0;
		for (j=6;j<=nparm; j++) dd[j]=0.0;

		x = Xi[i];
		plus5 = 5 + Xg[i];
		pwx = pow(x, p[5]);
		ex1 = exp(-p[1]-p[2]*pwx);
		ex2 = exp(-(p[3]+p[4]*x)*Ls[i]);
		ex = (1 - ex1) * ex2;
		PROBABILITY_INRANGE(&ex);
		j = (int) Yp[i];

		for (k = 1; k <= j; k++)
		{
			tm = ex + (k - 1) * p[plus5];
			tm1 += 1.0 / tm;
			tm1a += (1.0 / tm) * (k-1);
		}

		j = (int) Yn[i];

		for (k = 1; k <= j; k++)
		{
			tm = 1.0 - ex + (k - 1) * p[plus5];
			tm2 += 1.0 / tm;
			tm2a += (1.0 / tm) * (k - 1);
		}


		j = (int) (Yp[i]+Yn[i]);
		for (k = 1; k <= j; k++)
		{
			tm = 1.0 + (k - 1) * p[plus5];
			if (tm == 0.0)
			{
				tm = 1e-8;
			}
			tm3 += 1.0 / tm;
			tm3a += (1.0 / tm) * (k - 1);
		}

		tm12 = (tm1 - tm2);
		ex3 = ex1*ex2*tm12;

		dd[1] = ex1*ex2*tm12;     				//ex3;
		dd[2] = ex1*ex2*pwx*tm12; 				//ex3*pwx;
		dd[3] = -ex*Ls[i]*tm12;  			//-ex2*(1-ex1)*Ls[i]*(tm1-tm2);
		dd[4] = -x*ex*Ls[i]*tm12;  		//-x*ex2*(1-ex1)*Ls[i]*(tm1-tm2);
		if (Xi[i] <= 0.0) dd[5] = 0.0;
		else dd[5] = ex1*ex2*pwx*p[2]*log(Xi[i])*tm12;  	//ex3*pwx*p[2]*log(Xi[i]);

		dd[plus5] = (tm1a + tm2a -tm3a);

		for (j=1; j<=nparm; j++)
		{
			outgrad[j] += dd[j];
		}
	}


	for (i=1; i<=nparm; i++)
		grad[i] = outgrad[i];

	FREE_DVECTOR(p, 1, nparm);
	FREE_DVECTOR(outgrad, 1, nparm);

}

/*******************************************************************
**RaiVR_vcv -- used to compute the vcv for RaiVR model.
*		  Extern var.: Nobs, Xi, Yp, Yn, Ls, Xg.
*
******************************************************************/
int RaiVR_vcv(int nparm, int Spec[], double ptf[], double **vcv)
{
	int i, j, ivar, jvar, nvar;
	double *ptemp, *saveparms, *h, *gradp, *gradm, hrat, tmp;

	/* initialize memory for all the arrays */
	ptemp = DVECTOR(1, nparm);
	saveparms = DVECTOR(1, nparm);
	h = DVECTOR(1, nparm);
	gradp = DVECTOR(1, nparm);
	gradm = DVECTOR(1, nparm);

	for (i = 1; i <= nparm; i++)
    	{
      		ptemp[i] = ptf[i];
    	}

	/* Get a value of h for each parameter */
	hrat = pow(1.0e-16, 0.333333);

	for (i = 1; i <= nparm; i++)
	{
		if (fabs(ptemp[i] > 1.0e-7))
		{
			h[i] = hrat * fabs(ptemp[i]);
			tmp = ptemp[i] + h[i];
			h[i] = tmp - ptemp[i];
		}
		else
			h[i] = hrat;
	}

	
	/* initialize saveparms with the parameter values */
	for (i = 1; i <= nparm; i++)
	{
		saveparms[i] = ptemp[i];
	}

	/* for each unfixed parameter, compute the second derivative wrt each of the others */
	ivar = 0;
	nvar = 0;
	for (i = 1; i <= nparm; i++)
	{
		if (i > 1) saveparms[i-1] = ptemp[i-1];
		if (Spec[i] != Yes)
		{
			ivar++;
			nvar++;
			saveparms[i] = ptemp[i] + h[i];
			RaiVR_grad(nparm, Spec, saveparms, gradp);
			saveparms[i] = ptemp[i] - h[i];
			RaiVR_grad(nparm, Spec, saveparms, gradm);
			/* Now compute the second derivative */
			jvar = 0;
			for (j = 1; j <= nparm; j++)
			{
				if(Spec[j] != Yes)
				{
					jvar++;
					vcv[ivar][jvar] = - (gradp[jvar] - gradm[jvar])/(2.0 * h[i]);
				}
			}
		}
	}


	FREE_DVECTOR(ptemp, 1, nparm);
	FREE_DVECTOR(saveparms, 1, nparm);
	FREE_DVECTOR(h, 1, nparm);
	FREE_DVECTOR(gradp, 1, nparm);
	FREE_DVECTOR(gradm, 1, nparm);
	
	return nvar;

}/* End of RaiVR_vcv */

/**************************************************************
*MAX_lk -- used to obtain the Maxmun log-lilikehood as well as
*           the estimatiors of parameters, given initial p[1..n],
*           object func. , and gradient func. G_func.
*
**************************************************************/
void MAX_lk(int nparm, double p[], double gtol, int *iter, double *fret)
{
	int i, jfixed, jvar;
	long int nvar;
	long int *uiparm;
	double *urparm, *start, *lower, *upper;
	void (*ufparm)();

	/* Set up initial parameter array, start.  All parameters go either
	into start (if they are changing to improve the fit) or urparm
	(if they are fixed). */


	nvar = nparm;
	for (i = 1; i <= nparm; i++)
	{
		nvar = nvar - Spec[i];               /* Count the varying parameters */
	}

	uiparm = LIVECTOR(0,nparm + 2);
	urparm = DVECTOR(0, nparm-nvar);
	if (nvar > 0)
	{
		start  = DVECTOR(0, nvar);
		lower  = DVECTOR(0, nvar);
		upper  = DVECTOR(0, nvar);
	}
	else
	{
		long int dummy_nf;
		double dummy_x;

		for(i=1; i<=nparm; i++)
		{
			urparm[i-1]=p[i];
		}

		RaiVR_lk(&nvar, &dummy_x,  &dummy_nf, fret,
			uiparm, urparm, ufparm);

		*fret= -(*fret);
		ErrorFlag = -1;
		FREE_LIVECTOR(uiparm, 0, nparm+2);
		FREE_DVECTOR(urparm, 0, nparm-nvar);
		return;
	}

	/* seperate the fixed and variable parameters */
	jfixed = 0;
	jvar = 0;
	for (i = 1; i <= nparm; i++)
	{
		if (Spec[i] == Yes)
		{
			urparm[jfixed] = p[i];
			jfixed++;
		}
		else
		{
			start[jvar] = p[i];
			jvar++;
		}
	}

	/* set up the lower and upper variables.  Each parameter is unique, here. */
	jvar = 0;
	if (Spec[(int) alpha] == No)
	{
		lower[jvar] = 0.0;
		upper[jvar] = Max_double ;
		jvar++;
	}
	if (Spec[(int) beta] == No)
	{
		lower[jvar] = 0.0;
		upper[jvar] = Max_double;
		jvar++;
	}
	if (Spec[(int) theta1] == No)
	{
		lower[jvar] = 0.0;
		upper[jvar] = Max_double;
		jvar++;
	}

	if (Spec[(int) theta2] == No)
	{
		lower[jvar] = -1;  /* -1/xmax; */
		upper[jvar] = Max_double;
		jvar++;
	}
	if (Spec[(int) rho] == No)
	{
		if (restrict == Yes)
		{
			lower[jvar] = 1.0;
		}
		else
		{
			lower[jvar] = 0.0;
		}
		upper[jvar] = PowerUpperBound; /* Maybe this will be a variable someday */
		jvar++;
	}
	for(i=6; i<=nparm; i++)
	{
		if (Spec[i] == No)
		{
			lower[jvar] = 0.0;
			upper[jvar] = Max_double;
			jvar++;
		}
	}

	Maxloglik = 0.0; /* we don't have a good value here */

#ifdef LOGGING_ON
	fprintf(fp_log,"\n***************Inside MAX_lk**************************\n");
	fprintf(fp_log,"\n************Before call to run_dmngb******************\n");
	fprintf(fp_log,"nparm = %3d            (The number of parms in the model)\n",nparm);
	fprintf(fp_log,"nvar = %3ld             (The number of varying parameters)\n",nvar);
	fprintf(fp_log,"jfixed = %3d           (The number of fixed parameters)\n",jfixed);
	fprintf(fp_log,"jvar = %3d             (The number of restrictions on the parameters)\n",jvar);
	fprintf(fp_log,"Max_double = %6.3g  (The maximum double)\n",Max_double);
	fprintf(fp_log,"xmax = %6.3g          (The maximum dose)\n", xmax);
	fprintf(fp_log,"Maxloglik = %6.3g\n", Maxloglik);
	fprintf(fp_log,"Rel_Conv = %6.3g\n", Rel_Conv);
	fprintf(fp_log,"Parm_Conv = %6.3g\n", Parm_Conv);
	fprintf(fp_log,"ITMAX = %3d            (The maximum number of iterations)\n",ITMAX);
	fprintf(fp_log,"gtol = %6.3g          (The convergence critria)\n", gtol);
	fprintf(fp_log,"************Starting Vules & Bounds Going Into DMNGB*************\n");
	for (i=0; i<jvar; i++)
		fprintf(fp_log,"start[%2d] = %12.6g    lower[%2d] = %12.6g    upper[%2d] = %12.6g\n",
		i,start[i],i,lower[i],i,upper[i]);
	fprintf(fp_log,"************************************************\n");
	for (i=0; i<nparm; i++)
		fprintf(fp_log,"uiparm[%2d] = %12.6ld\n",i,uiparm[i]);
	fprintf(fp_log,"************************************************\n");
	for (i=0; i<jfixed; i++)
		fprintf(fp_log,"urparm[%2d] = %12.6g\n",i,urparm[i]);
	fprintf(fp_log,"**************************************************************\n");
	fflush(fp_log);
#endif
	/* maximize the log-likelihood.  ErrorFlag gives the return code */
	ErrorFlag = run_dmngb((int) nvar, start, lower, upper, Maxloglik,
		Rel_Conv, Parm_Conv, ITMAX, 10,
		RaiVR_lk, RaiVR_g, uiparm, urparm, ufparm,
		DeBuG, fret);

	/* We actually minimized -log-likelihood.  Want to return log-likelihood */
	*fret = -*fret;

	/* Put the parameter values back in p */
	jfixed=jvar=0;
	for (i = 1; i <= nparm; i++)
	{
		if (Spec[i]==Yes)
		{
			p[i] = urparm[jfixed];
			jfixed++;
		}
		else
		{
			p[i] = start[jvar];
			jvar++;
		}
	}

#ifdef LOGGING_ON
	fprintf(fp_log,"\n***************Inside MAX_lk**************************\n");
	fprintf(fp_log,"\n************After call to run_dmngb*******************\n");
	fprintf(fp_log,"*fret = %12.6g   (Maximum log-likelihood value)\n",*fret);
	fprintf(fp_log,"ErrorFlag = %2d         (Status of dmngb; Good <= 6)\n",ErrorFlag);
	fprintf(fp_log,"nvar = %3ld             (The number of varying parameters)\n",nvar);
	fprintf(fp_log,"jfixed = %3d           (The number of fixed parameters)\n",jfixed);
	fprintf(fp_log,"jvar = %3d             (The number of restrictions on the parameters)\n",jvar);
	fprintf(fp_log,"Maxloglik = %6.3g\n", Maxloglik);
	fprintf(fp_log,"ITMAX = %3d            (The maximum number of iterations)\n",ITMAX);
	fprintf(fp_log,"************Parmeter Estimates & Bounds Coming Out DMNGB*************\n");
	for (i=0; i<jvar; i++)
		fprintf(fp_log,"start[%2d] = %12.6g    lower[%2d] = %12.6g    upper[%2d] = %12.6g\n",
		i,start[i],i,lower[i],i,upper[i]);
	fprintf(fp_log,"************************************************\n");
	for (i=0; i<nparm; i++)
		fprintf(fp_log,"uiparm[%2d] = %12.6ld\n",i,uiparm[i]);
	fprintf(fp_log,"************************************************\n");
	for (i=0; i<jfixed; i++)
		fprintf(fp_log,"urparm[%2d] = %12.6g\n",i,urparm[i]);
	fprintf(fp_log,"**************************************************************\n");
	fflush(fp_log);
#endif
	FREE_LIVECTOR(uiparm, 0, nparm+2);
	FREE_DVECTOR(urparm, 0, nparm-nvar);
	FREE_DVECTOR(start, 0, nvar);
	FREE_DVECTOR(lower, 0, nvar);
	FREE_DVECTOR(upper, 0, nvar);
}/* End of MAX_lk */

/**************************************************************
*RaiVR_fit -- Uesd to "prepare" the data for further computation,
*            i.e. compute the extern variables, give the initial
*            parameters, etc. THEN fit the RaiVR model.
*            (In fact, these jobs could be done in main().)
*
***************************************************************/
void RaiVR_fit(int nparm, int ngrp, double p[], double gtol,
	int *iter, double *fret)
{
	int    *SpBak, i, j, junk, count;
	double smin, smax, sum, nij, ymin, W, tmvcv, tmy, xlk, xlkold, delta;
	double *tmYp, *tmYn, *tmXi, *pBak;
	double *newY;   /* This holds the values for finding starting values
					for theta1 and theta2 */
	double meanX, meanY, sumXY, sumX2;

	tmYn  = DVECTOR(1, ngrp);
	tmYp  = DVECTOR(1, ngrp);
	tmXi  = DVECTOR(1, ngrp);
	pBak  = DVECTOR(1, nparm);
	SpBak = IVECTOR(1, nparm);
	newY  = DVECTOR(1, Nobs);



	replace = No;

	sum  = 0.0;
	nij  = 0.0;
	xmax = Xi[1];
	xmin = Xi[1];
	smax = Ls[1];
	smin = Ls[1];

	for (i=1; i<=Nobs; i++)
	{
		sum += Yp[i] + Yn[i];
		nij  +=  (Yp[i]+Yn[i])/Ls[i];

		if (Ls[i]>smax)
		{
			smax = Ls[i];        /* litter size max */
		}
		if (Ls[i]<smin)
		{
			smin = Ls[i];        /* litter size min */
		}
		if (Xi[i] > xmax)
		{
			xmax = Xi[i];        /* dose level max */
		}
		if (Xi[i] < xmin)
		{
			xmin = Xi[i];        /* dose level min */
		}
		if (Xg[i]==1 && Xg[i+1]==2 && fixedSize==Yes)
		{
			sijfixed = sum/nij;  /* mean litter size across control group */
		}
	}

	if(fixedSize==No)
	{
		sijfixed = sum/nij;      /* mean litter size across all groups */
	}

	if (smax<=1.0 || sum <=0.0)
		ERRORPRT("all litter sizes are zero.");

	if(initial==Yes)
	{
		for (i=1; i<=nparm; i++)
		{
			SpBak[i]= 1 - IniSp[i];
		}

		//Have to do this because the fun OUTPUT_Init.
		OUTPUT_TEXT("\n\n                 User Inputs Initial Parameter Values  ");
		OUTPUT_Init(nparm, IniSp, IniP, Parm_name);
		//obtain user input initial values for unspecified parms.
		for (i=1; i<=nparm; i++)
		{
			if(IniSp[i]==1)    // have been initialized.
			{
				if(Spec[i]==1 )  // check if it is for fixed parm.
				{
					Warning("The initial value for the fixed parameter is ignored.");
				}
				//p[i]=p[i]. no change.
				else
				{
					p[i]=IniP[i];
				}
			}
			else
			{
				//check if all the unspecified parms were initialized.
				if (Spec[i]==0)
				{
					ERRORPRT("When the initial option is chosen, one has to initial ALL unspecified parameters.");
				}
			}
		}

		if (p[(int) alpha] < 0 || p[(int) beta]<0 || p[(int) rho]<restrict)
		{
			ERRORPRT("The initial values have to be: alpha >= 0,  beta >= 0 and rho >= 0 (or 1 when there is restriction on Rho). ");
		}
		if (p[(int) theta1]<0 || p[(int) theta1]+p[(int) theta2]*xmax<0)
		{
			ERRORPRT("The initial values have to be: theta1 >= 0 and theta1+theta2*Dose >= 0 for all dose.");
		}
		for (j=6; j<=nparm; j++)
		{
			if (p[j]<0 || p[j]>=1)
			{
				ERRORPRT("The initial values have to be: 0 <= Phi[j] < 1, for the correlation parameters.");
			}
		}
	}
	else
	{
		for (i=1; i<=nparm; i++)  //save the Spec[] and p[].
		{
			SpBak[i]= Spec[i];
			pBak[i]=p[i];
		}

		/**First, get the initial for RaiVR model.**********/
		i=1;
		for (j=1; j<=ngrp; j++) //converse nest_data to dicho_data.
		{
			tmYn[j]=0; tmYp[j]=0;
			while (i <= Nobs && Xg[i]==j)
			{
				tmYn[j] += Yn[i];
				tmYp[j] += Yp[i];
				tmXi[j]  = Xi[i];
				i++;
			}
			if (i > Nobs) break;
		}

		ymin = 1.0;
		tmvcv = 0.0;
		tmy = 0.0;
		if(Spec[(int) rho]==0)
		{
			p[(int) rho]=1.001;
		}
		for (i=1; i<=ngrp; i++)
		{
			W = tmYp[i]/(tmYp[i]+tmYn[i]);

			if (ymin>W)
			{
				ymin = W;
			}
			PROBABILITY_INRANGE(&W);
			W = -1.0 * log(1.0-W);
			tmvcv += pow(tmXi[i],p[(int) rho]);
			tmy += W;
		}
		p[(int) alpha] = ymin;
		p[(int) beta] = tmy/tmvcv;
		//GETKNOWNPARMS(3,SpBak, p, pBak);
		for (i=1; i<=2; i++)
		{
			if (SpBak[i]==1)
			{
				p[i]=pBak[i];
			}
		}
		if (SpBak[(int) rho]==1)
		{
			p[(int) rho]=pBak[(int) rho];
		}

		Spec[(int) theta1]=Spec[(int) theta2]=1;     //setup Weibull model.
		p[(int) theta1]=p[(int) theta2]=0.0;

		for (i=6; i<=nparm; i++)
		{
			Spec[i]=1;
			p[i]=0.0;
		}

		/* Scale starting values by max dose level */
		p[2] = p[2]*pow(maxdose, p[5]);
		p[4] = p[4]*maxdose;

#ifdef LOGGING_ON
		fprintf(fp_log,"\n***************Inside RaiVR_fit***********************\n");
		fprintf(fp_log,"\n***********Before call to first MAX_lk****************\n");
		fprintf(fp_log,"\nnparm = %d \n\n",nparm);
		fprintf(fp_log,"p[ 1] = %12.6g  (alpha)\n",p[1]);
		fprintf(fp_log,"p[ 2] = %12.6g  (beta)\n",p[2]);
		fprintf(fp_log,"p[ 3] = %12.6g  (theta1)\n",p[3]);
		fprintf(fp_log,"p[ 4] = %12.6g  (theta2)\n",p[4]);
		fprintf(fp_log,"p[ 5] = %12.6g  (rho)\n",p[5]);
		for(i=6; i<=nparm; i++)
		{
			fprintf(fp_log,"p[%2d] = %12.6g  (phi%d)\n",i,p[i],i-5);
		}
		fprintf(fp_log,"******************************************************\n");
		MxLkCnt = 1;
		fflush(fp_log);
#endif
		MAX_lk(nparm, p, gtol, &junk, &xlk);

		/*** end of the initial estimation from RaiVR model. ***/

#ifdef LOGGING_ON
		fprintf(fp_log,"\n****************Inside RaiVR_fit**********************\n");
		fprintf(fp_log,"\n************After call to first MAX_lk****************\n");
		fprintf(fp_log,"\nnparm = %d \nxlk = %12.6g \n\n",nparm,xlk);
		fprintf(fp_log,"p[ 1] = %12.6g  (alpha)\n",p[1]);
		fprintf(fp_log,"p[ 2] = %12.6g  (beta)\n",p[2]);
		fprintf(fp_log,"p[ 3] = %12.6g  (theta1)\n",p[3]);
		fprintf(fp_log,"p[ 4] = %12.6g  (theta2)\n",p[4]);
		fprintf(fp_log,"p[ 5] = %12.6g  (rho)\n",p[5]);
		for(i=6; i<=nparm; i++)
		{
			fprintf(fp_log,"p[%2d] = %12.6g  (phi%d)\n",i,p[i],i-5);
		}
		fprintf(fp_log,"*******************************************************\n");
		fflush(fp_log);
#endif
		/* Uncale parameter estimates by max dose level */
		p[2] = p[2]/pow(maxdose, p[5]);
		p[4] = p[4]/maxdose;

		/*** Second, get initial values for Phi's.      ***********/
		count=0;
		for (i=6; i<=nparm; i++)
		{
			count += SpBak[i];
		}
		if (count < ngrp)
		{
			for(i=6; i<=nparm; i++)         //Spec[(int) theta1], [(int) theta2] remain as 1.
			{
				Spec[i] = SpBak[i];      //set back to original ones.
				p[i] = pBak[i];
				if (SpBak[i]==0)
				{
					p[i]=0.01;    //set init. if unknown.
				}
			}
			for(j=6; j<=nparm; j++)
			{
				p[j]=p[j]/(1-p[j]);     // Phi --> Psi.
			}

			if( p[(int) theta1] > 0)
			{
				p[(int) theta2] = p[(int) theta2]/p[(int) theta1];
			}
			else
			{
				p[(int) theta2] = 0;
			}

			/* Scale starting values by max dose level */
			p[2] = p[2]*pow(maxdose, p[5]);
			p[4] = p[4]*maxdose;

#ifdef LOGGING_ON
			fprintf(fp_log,"\n****************Inside RaiVR_fit***********************\n");
			fprintf(fp_log,"\n***********Before call to second MAX_lk****************\n");
			fprintf(fp_log,"\nnparm = %d\n",nparm);
			fprintf(fp_log,"p[ 1] = %12.6g  (alpha)\n",p[1]);
			fprintf(fp_log,"p[ 2] = %12.6g  (beta)\n",p[2]);
			fprintf(fp_log,"p[ 3] = %12.6g  (theta1)\n",p[3]);
			fprintf(fp_log,"p[ 4] = %12.6g  (theta2)\n",p[4]);
			fprintf(fp_log,"p[ 5] = %12.6g  (rho)\n",p[5]);
			for(i=6; i<=nparm; i++)
			{
				fprintf(fp_log,"p[%2d] = %12.6g  (psi%d)\n",i,p[i],i-5);
			}
			fprintf(fp_log,"*******************************************************\n");
			MxLkCnt = 2;
			fflush(fp_log);
#endif
			MAX_lk(nparm, p, gtol, &junk, &xlk);

#ifdef LOGGING_ON
			fprintf(fp_log,"\n*****************Inside RaiVR_fit**********************\n");
			fprintf(fp_log,"\n************After call to second MAX_lk****************\n");
			fprintf(fp_log,"\n**********Default Initial Parameter Values*************\n");
			fprintf(fp_log,"\nnparm = %d \nxlk = %12.6g \n\n",nparm,xlk);
			fprintf(fp_log,"p[ 1] = %12.6g  (alpha)\n",p[1]);
			fprintf(fp_log,"p[ 2] = %12.6g  (beta)\n",p[2]);
			fprintf(fp_log,"p[ 3] = %12.6g  (theta1)\n",p[3]);
			fprintf(fp_log,"p[ 4] = %12.6g  (theta2)\n",p[4]);
			fprintf(fp_log,"p[ 5] = %12.6g  (rho)\n",p[5]);
			for(i=6; i<=nparm; i++)
			{
				fprintf(fp_log,"p[%2d] = %12.6g  (psi%d)\n",i,p[i],i-5);
			}
			fflush(fp_log);
#endif
			/* Uncale parameter estimates by max dose level */
			p[2] = p[2]/pow(maxdose, p[5]);
			p[4] = p[4]/maxdose;

			for(j=6; j<=nparm; j++)
			{
				p[j]=p[j]/(1+p[j]);     // Psi --> Phi.
			}

			p[(int) theta2]= p[(int) theta2]*p[(int) theta1];

#ifdef LOGGING_ON
			fprintf(fp_log,"\nThese have been transformed back for output\n");
			fprintf(fp_log,"p[ 3] = %12.6g  (theta1)\n",p[3]);
			fprintf(fp_log,"p[ 4] = %12.6g  (theta2)\n",p[4]);
			for(j=6; j<=nparm; j++)
			{
				fprintf(fp_log,"p[%2d] = %12.6g  (phi%d)\n",j,p[j],j-5);
			}
			fprintf(fp_log,"********************************************************\n");
			fflush(fp_log);
#endif
		}

		/* Find starting values for theta1 and theta2 but also make sure they stay within
		the constraint. */
		meanX = meanY = sumXY = sumX2 = 0.0;
		for (i = 1; i <= Nobs; i++)
		{
			if ((Yp[i] > 0) && ((-p[1]-p[2]*Xi[i]) != 0))
				newY[i] = -(1/Ls[i])*log((Yp[i]/(Yp[i]+Yn[i]))/(1-exp(-p[1]-p[2]*Xi[i])));
			else
				if (Yp[i] > 0)
					newY[i] = -(1/Ls[i])*log((Yp[i]/(Yp[i]+Yn[i]))/(1e-8));
				else
					if ((-p[1]-p[2]*Xi[i]) != 0)
						newY[i] = -(1/Ls[i])*log((1e-8/(Yp[i]+Yn[i]))/(1-exp(-p[1]-p[2]*Xi[i])));
					else
						newY[i] = -(1/Ls[i])*log(1e-8);
			meanX += Xi[i];
			meanY += newY[i];
			sumXY += Xi[i]*newY[i];
			sumX2 += pow(Xi[i], 2);
		}
		meanX = meanX/Nobs;
		meanY = meanY/Nobs;
		p[4] = (sumXY - Nobs*meanX*meanY)/(sumX2 - Nobs*pow(meanX,2));
		p[3] = meanY - p[4]*meanX;
		if ((p[3]/(-maxdose)) >= p[4])
			p[4] = p[3]/(-maxdose) + 1e-8;


		/** Finally, get initial for Theta's.       **************/
		for (i=3; i<=4; i++)
		{
			Spec[i] = SpBak[i];
			/*  	  p[i] = pBak[i]; */
			if (SpBak[i]==1)
			{
				p[i]=pBak[i];  // have to put this one.
			}
		}
#ifdef MISC_OUT
		printf("\n\n                      Default Initial Parameter Values\n");
		printf("                          alpha = %12.5g\n", p[1]);
		printf("                           beta = %12.5g\n", p[2]);
		printf("                         theta1 = %12.5g\n", p[3]);
		printf("                         theta2 = %12.5g\n", p[4]);
		printf("                            rho = %12.5g\n", p[5]);
#endif
#ifdef MISC_OUT
		for (i = 1; i <= ngrp; i++)
			printf("                           phi%1d = %12.5g\n", i, p[i+5]);
#endif
		OUTPUT_TEXT("\n\n                  Default Initial Parameter Values  ");
		OUTPUT_Init(nparm, Spec, p, Parm_name);
		fflush(fp_out);
	}

	/*** Fit the model.              ************************/
	for(j=6; j<=nparm; j++)
	{
		p[j]=p[j]/(1-p[j]);     // Phi --> Psi.
	}

	if( p[(int) theta1] > 0)
	{
		p[(int) theta2] = p[(int) theta2]/p[(int) theta1];
	}
	else
	{
		p[(int) theta2] = 0;
	}

	xlkold = -Max_double/1000;
	delta=1.0;
	i=1;

	while(i<=5 && delta > 0.01)
	{

		/* Scale starting values by max dose level */
		p[2] = p[2]*pow(maxdose, p[5]);
		p[4] = p[4]*maxdose;

#ifdef LOGGING_ON

		fprintf(fp_log,"\n******************Inside RaiVR_fit**********************\n");
		fprintf(fp_log,"\n***********Before call %d to third MAX_lk****************\n",i);
		fprintf(fp_log,"\nnparm = %d \nxlk = %12.6g \nxlkold = %12.6g\n\n",nparm,xlk,xlkold);
		fprintf(fp_log,"p[ 1] = %12.6g  (alpha)\n",p[1]);
		fprintf(fp_log,"p[ 2] = %12.6g  (beta)\n",p[2]);
		fprintf(fp_log,"p[ 3] = %12.6g  (theta1)\n",p[3]);
		fprintf(fp_log,"p[ 4] = %12.6g  (theta2)\n",p[4]);
		fprintf(fp_log,"p[ 5] = %12.6g  (rho)\n",p[5]);
		for(j=6; j<=nparm; j++)
		{
			fprintf(fp_log,"p[%2d] = %12.6g  (psi%d)\n",j,p[j],j-5);
		}
		fprintf(fp_log,"********************************************************\n");
		fflush(fp_log);

#endif
		MAX_lk(nparm, p, gtol, &junk, &xlk);//first fit.

#ifdef LOGGING_ON
		fprintf(fp_log,"\n*****************Inside RaiVR_fit***********************\n");
		fprintf(fp_log,"\n************After call %d to third MAX_lk****************\n",i);
		fprintf(fp_log,"\nnparm = %d \nxlk = %12.6g \ndelta = %12.6g \n\n",nparm,xlk,xlk-xlkold);
		fprintf(fp_log,"p[ 1] = %12.6g  (alpha)\n",p[1]);
		fprintf(fp_log,"p[ 2] = %12.6g  (beta)\n",p[2]);
		fprintf(fp_log,"p[ 3] = %12.6g  (theta1)\n",p[3]);
		fprintf(fp_log,"p[ 4] = %12.6g  (theta2)\n",p[4]);
		fprintf(fp_log,"p[ 5] = %12.6g  (rho)\n",p[5]);
		for(j=6; j<=nparm; j++)
		{
			fprintf(fp_log,"p[%2d] = %12.6g  (psi%d)\n",j,p[j],j-5);
		}
		fprintf(fp_log,"********************************************************\n");
		fflush(fp_log);
#endif
		/* Uncale parameter estimates by max dose level */
		p[2] = p[2]/pow(maxdose, p[5]);
		p[4] = p[4]/maxdose;

		if (ErrorFlag == 9)
		{
			Warning("Warning: Maximum iteration may be not large enough. Iterations reaches the maximum.");
		}
		delta = xlk - xlkold;
		if (delta>0.01)
		{
			xlkold = xlk;
			for (j=1; j<=nparm; j++)
			{
				pBak[j]=p[j];  //save the good one.
			}
			//restart at the near point.
			for (j=1; j<=nparm; j++)
			{
				if (SpBak[j] == 0)
				{
					p[j]= p[j] * 1.001;
				}
			}
		}
		// printf("\n Repeating Number = %d      In(L) =%14.6lf  ", i, xlk);
		//fprintf(fp_out, "\n  Repeating Number = %d      In(L) =%14.6lf  ", i, xlk);
		i++;
	}

	do_dmngb_warning(&ErrorFlag);

	*fret = xlk;
	for (j=1; j<=nparm; j++)
	{
		p[j]=pBak[j];        //take the last good one.
	}

	for(j=6; j<=nparm; j++)
	{
		p[j]=p[j]/(1+p[j]);     // Psi --> Phi.
	}

	p[(int) theta2]= p[(int) theta2]*p[(int) theta1];

	FREE_DVECTOR(tmYn, 1, ngrp);
	FREE_DVECTOR(tmYp, 1, ngrp);
	FREE_DVECTOR(tmXi, 1, ngrp);
	FREE_DVECTOR(pBak, 1, nparm);
	FREE_IVECTOR(SpBak, 1, nparm);
	FREE_DVECTOR(newY, 1, Nobs);
}/* End of RaiVR_fit */

/************************************************************
* RaiVR_BMD -- Used to calculate the BMD and BMDL for RaiVR model.
*
*************************************************************/
void RaiVR_BMD (int nparm, double p[], double gtol, int *iter, double xlk,
	double Rlevel[], double Bmdl[], double *BMD)
{
	double BMD_func(int nparm, double p[], double x, double gtol);
	double BMDL_func(int nparm, double p[], double x, double gtol);

	double   tol;
	double   xa, xb, fa, fb, xmax3;
	double   Bmdjunk;
	double   *pBak, *porig;
	int      j, k;

	pBak=DVECTOR(1, nparm);
        porig=DVECTOR(1, nparm);

        for (j=1; j<=nparm; j++) porig[j] = p[j];   /* save orig value of p */

	/**** compute X^2 value  ************************/
	/* IF ML is the value of the maximized log-likelihood, then ML - LR is the value
	log-likelihood at the BMDL or BMDU */
	if (bmdparm.level < 0.5)
	{
		LR = QCHISQ(1.0 - 2 * bmdparm.level,1) / 2.0;
	}
	else
	{
		LR = QCHISQ(2 * bmdparm.level - 1.0,1) / 2.0;
	}

	tol = sqrt(DBL_EPSILON);

	/* Transform parameters into "internal" form */
	for(j=6; j<=nparm; j++)
	{
		p[j] = p[j]/(1-p[j]);     /* Phi --> Psi */
	}

	if( p[(int) theta1] > 0)
	{
		p[(int) theta2] = p[(int) theta2]/p[(int) theta1];
	}
	else
	{
		p[(int) theta2] = 0;
	}

	for(j=1; j<=nparm; j++)
	{
		pBak[j]= p[j];         /* save the p[]. */
	}

	Rlevel[1] = BMR = bmdparm.effect;


	/**** solve the BMD ********************************/
	fa = BMD_func(nparm, pBak, 0, tol);
	fb = BMD_func(nparm, pBak, xmax, tol);
	xmax3 = xmax;

	if(fa*fb>=0)
	{
		fb=BMD_func(nparm, pBak, 3*xmax, tol);
		if(fa*fb>=0)
		{
#ifndef RBMDS
			fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  No); //computation failed
#endif
			ERRORPRT("BMD computation failed: BMD is out of the three times dose boundary.");
		}
		else
		{
			xmax3=3*xmax;
		}
	}

	/* the BMD cannot be solved analytically, so zeroin is used  */
	/* in conjunction with BMD_func to find the BMD numerically  */
	/* i.e. the RaiVR model cannot be solved for dose explicitly */

	xb = zeroin(0, xmax3, 1.0e-10, BMD_func, nparm, p, 1.0e-14);
	*BMD = xb;

	Spec[(int) beta] = Yes;
	replace = Yes;        //replace is extern var., has to be changed now

	/********* search for BMDL **************************/
	xa = xb*0.5;

	BMD_lk = xlk;  //get the lk at BMD.
	fb = -LR;
	fa = BMDL_func(nparm, p, xa, tol);
        
 
	while (fa<0.0 && xa > DBL_EPSILON)
	{
		xa *= 0.5; //prevent that xa=0.1*BMD is not small enough.
		fa = BMDL_func(nparm, p, xa, tol);
	}
	if (fa<0.0)
	{
		ERRORPRT("Benchmark dose computation failed.  Lower limit includes zero.");
	}

#ifndef RBMDS
	fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  Yes); //computation will succeed
#endif
	Bmdl[1] = zeroin(xa, xb, 1.0e-10, BMDL_func, nparm, p, 1.0e-14);

	if(bmdlCurve==Yes)
	{
		/****** calculate Bmd[] and Bmdl[] ***********/
		for (k=2; k<=5; k++)
		{
			for(j=1; j<=nparm; j++)
			{
				p[j]= pBak[j];          //get the "old" p[].
			}

			if (k==2)
			{
				Rlevel[k]=BMR=0.05;
			}
			else
			{
				Rlevel[k]= BMR = (k-2)*0.1;
			}

			/**** solve the BMD []********************************/
			fa=BMD_func(nparm, pBak, 0, tol);
			fb=BMD_func(nparm, pBak, xmax, tol);
			if(fa*fb>=0)
			{
				xb=xmax;
				Bmdjunk = xb;
				fprintf(fp_out, "\n BMDL curve computation failed for BMR = %f . \n The BMDL curve appearing in the graph may not accurate.", BMR);

			}
			else
			{
				xb = zeroin(0, xmax, 1.0e-10, BMD_func, nparm, p, 1.0e-14);
				Bmdjunk = xb;
			}

			/********* search for BMDL []**************************/
			xa = xb*0.1;
			tol = FMAX((Bmdjunk)*0.001, 0.0000001);
			BMD_lk = xlk;  //get the lk at BMD.
			fb = -LR;
			fa = BMDL_func(nparm, p, xa, tol);

			if (fa<0.0)
			{
				xa *= 0.01; //prevent that xa=0.1*BMD is not small enough.
				fa = BMDL_func(nparm, p, xa, tol);
				if (fa<0.0)
				{
					Bmdl[k]=-1;
					fprintf(fp_out, "\n BMDL curve computation failed for BMR = %f . \n The BMDL curve appearing in the graph may not acurate.", BMR);
				}
				else
				{
					Bmdl[k] = zeroin(xa, xb, 1.0e-10, BMDL_func, nparm, p, 1.0e-14);
				}
			}
			else
			{
				Bmdl[k] = zeroin(xa, xb, 1.0e-10, BMDL_func, nparm, p, 1.0e-14);
			}
		}
	}
	else
	{
		for (k=2; k<=5; k++)
		{
			Bmdl[k] = Rlevel[k]= -1;
		}
	}

        for (j=1; j<=nparm; j++) p[j] = porig[j];   /* reset orig value of p */

        FREE_DVECTOR(porig, 1, nparm);
	FREE_DVECTOR(pBak, 1, nparm);
}/* End of RaiVR_BMD */

/*****************************************************************
* BMD_func -- used to compute the values of functions BMD_f at
*            the point D, given the parm p[] and number of parm.
*
*            This routine is called by zeroin().
*
*****************************************************************/
double BMD_func(int n, double p[], double D, double gtol)
{
	double fx, P0, PR, lambda, S;

	S=sijfixed;
	P0=(1-exp(-p[(int) alpha]))*exp(-S*p[(int) theta1]);               /* Prob.(Dose=0) i.e. Background */
	lambda=exp( p[(int) theta1]*(1+p[(int) theta2]*D)*S );
	PR = BMR + P0;                                //Added.

	if (bmdparm.risk==EXTRA)
	{
		PR=BMR+(1-BMR)*P0;      //extra
	}

	fx = p[(int) alpha]+p[(int) beta]*pow(D, p[(int) rho]) + log(1-PR*lambda);

	return fx;
}/* End of BMD_func */

/*****************************************************************
* BMDL_func -- used to compute the values of functions BMDL_f (the
*              X^2 value) at the point D, given the parm p[] and
*              number of parm.
*
*              This routine is called by zeroin().
*
*****************************************************************/
double BMDL_func(int nparm, double pBak[], double D, double gtol)
{ 	/* BMD_lk and LR are calculated in RaiVR_BMD() */

	double fD, xlk, *p;
	int j, junk;
	double P0, PR, lambda, S;
	p=DVECTOR(1,nparm);

	for (j=1;j<=nparm;j++)
	{
		p[j]=pBak[j];  //get the "old" p[].
	}

	S=sijfixed;
	P0=(1-exp(-p[(int) alpha]))*exp(-S*p[(int) theta1]);               //Prob.(Dose=0)
	lambda=exp( p[(int) theta1]*(1+p[(int) theta2]*D)*S );
	PR = BMR + P0;                                //Added.

	if (bmdparm.risk==EXTRA)
	{
		PR=BMR+(1-BMR)*P0;      //extra
	}

	tD = D;   //tD is global var. have to change before call MAX_lk().

	if (1-lambda*PR<=0)
	{
		fD = -Max_double/1000;
		return fD;
	}
	else
	{
		p[(int) beta]= (-log(1-lambda*PR)-p[(int) alpha])/pow(tD, p[(int) rho]);
	}

	MAX_lk( nparm, p, gtol,  &junk,  &xlk);
	fD = BMD_lk - xlk - LR;

	FREE_DVECTOR(p, 1, nparm);
	return fD;
}/* End of BMDL_func */

/******************************************************
**  TEMP_ANOVA_OUTPUT--output ANOVA table values for a
*                      3-parameter dichotomous model.
*******************************************************/
void TEMP_ANOVA_OUTPUT (char *anatxt[], AnaList anasum[])
{
	double AIC;
	AIC = -2*anasum[2].SS + 2*(1.0 + anasum[3].DF - anasum[2].DF);
	/*output to screen */
	/*  OUTPUT_TEXT ("\n\n\n                      Analysis of Deviance Table");*/
	/*  OUTPUT_TEXT ("\n Fitted Model Log(likelihood)"); */

	/*output to *.out file */
	/*  fprintf (fp_out, "                %15s %15.6g\n", anatxt[0], anasum[1].SS); */
#ifndef RBMDS
	fprintf (fp_out, "\n Log-likelihood: %7.6g   AIC: %7.6g\n", anasum[2].SS, AIC);
#else
	fprintf (fp_out, "\n Log-likelihood: %30.22g AIC: %30.22g\n", anasum[2].SS, AIC);
#endif
	/*  fprintf (fp_out, "                %15s %15.6g\n", anatxt[2], anasum[3].SS);*/

#ifndef RBMDS
	/*   fprintf (fp_out,"\n                %15s %15.6g\n", "AIC:", */
	/* 	   -2*anasum[2].SS + 2*(1.0 + anasum[3].DF - anasum[2].DF)); */
#else
	/*   fprintf (fp_out,"\n                %15s %30.22g\n", "AIC:", */
	/* 	   -2*anasum[2].SS + 2*(1.0 + anasum[3].DF - anasum[2].DF)); */
#endif
}
