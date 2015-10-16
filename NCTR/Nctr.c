/*******************************************************************
**
* IMPORTANT NOTE:  The following variable is the version number for
*                  the current model.  THIS MUST BE CHANGED as
*				   important changes are made to the models.
*
********************************************************************/
char Version_no[]="NCTR Model. (Version: 2.13; Date: 04/27/2015)";

/*******************************************************************
*
* NCTR.C - a ANSI C program for NCTR model fitting with/without
*          a natural background rate in Benchmark Dose.
*
* Date: February 16, 2001
*
*
********************************************************************
* Modification Log:
*
* Version Number: 2.3
* Modified By: Micheal Ferree
* Modified Date: 8/07/2005
* Reason: Changed to use Nest_CI function.
*
* Version Number: 2.4
* Modified By: Woodrow Setzer
* Modified Date: 10/27/2005
* Reason:
*   1) Fixed (minor) memory access errors
*   2) Added conditional compilation for RBMDS and logging
*   3) Report only likelihood for fitted model (not reduced and full
*      likelihoods).
*   4) Modified NCTR_BMD so searching starts closer to BMD.
*
* Version Number: 2.5
* Modified By: Woodrow Setzer
* Modified Date: 03/22/2006
* Reason: Replaced INVS_Mat with INVMAT, with code modifications to
*         allow for the difference.
*
* Version Number: 2.6
* Modified By: Geoffrey
* Date: 1/12/2007
* Reason: Incremented version number.
*		 Added last parameter "0" (don't print SE) in OUTPUT_DTMS3PARMS().
*
* Version Number: 2.7
* Modified By: Woodrow Setzer
* Date: 2/20/2007
* Reason: Incremented version number to reflect changed compilation options.
*
* Version Number: 2.8
* Modified By: G. Nonato
* Modification Date: 04/10/2008
* Reason: (Per BMDS 2.0: Problem Report 157 & 147)
*       Fix the Observation # < parameter # for NCTR model problem.
*       Added code to free-up allocated memories before exiting thru ERRORPRT()
*
* Version Number: 2.9
* Modified By: G. Nonato
* Modification Date: 10/28/2009
* Reason:
*      To be able to process files/folders with spaces (PR 257)
*      Fix program freeze due to long variable names (PR 278)
*      Process long files up to 256 characters (PR 303 and 308)
*      Modify code for easy maintenance, all lengths for file name,
*        model name, and column names are placed in benchmark.h
*
* Version Number: 2.10
* Modified By: Louis Olszyk
* Modification Date: 02/28/2013
* Reason: PR 444 - Fix wording in plot titles
*
* Version Number: 2.11
* Modified By: Cody Simmons
* Modification Date: 09/08/2014
* Reason: PR232 and PR244
*       Added the bootstrap method for goodness of fit calculation
*       PR349 - Fix error caused when LSC is much larger than litter size.
*
* Version Number: 2.12
* Modified By: Cody Simmons
* Modification Date: 09/15/2014
* Reason: PR300
*       Created scaled residual of interest table for nested models.
*
* Version Number: 2.13
* Modified By: Cody Simmons
* Modification Date: 4/27/2015
* Reason: PR545
*       Display the min, max and mean of absolute value of scaled residuals.
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

void NCTR_fit(int n, int ngrp, double p[], double gtol,
			  int *iter, double *fret);
void NCTR_BMD (int nparm, double p[], double gtol, int *iter, double xlk,
			   double Rlevel[], double Bmdl[], double *BMD);
void NCTR_grad(int nparm, int Spec[], double ptf[], double grad[]);
int NCTR_vcv(int nparm, int Spec[], double p[], double **vcv);
int Model_DF (int []);
void Which_Bounded (int [], double [], int []);
void Predict(double [], double [], int, double [], double []);
double BMDL_func(int nparm, double p[], double x, double gtol);
void initialparms(int, double, double, double [], double [],
				  double [], double);
void TEMP_ANOVA_OUTPUT(char *anatxt[], AnaList anasum[]);

#define EPS 3.0e-8
#define TOLX (10*EPS)
#define STPMX1 10.0
#define ALF 0.000001
#define FREEALL FREE_DVECTOR (xi,1,n); FREE_DVECTOR (pnew, 1,n);\
	FREE_DVECTOR(hdg,1,n); FREE_DVECTOR(g,1,n);\
	FREE_DVECTOR(dg,1,n);\
	FREE_DMATRIX(hessin,1,n,1,n);
#define float double
#define GOLD 1.618034
#define GLIMIT 100
#define TINY 1.0e-20
#define SWAP(a,b, junk) (junk)=(a); (a)=(b); (b)=(junk);
#define SHFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d);
#define ZEPS 1.0e-8
#define MOV3(a,b,c, d,e,f)(a)=(d);(b)=(e);(c)=(f);
#define TOL 2.0e-6
#define BSDEBUG 1000		/* Bootstrapping debug code (0 = NoDebug, N>0 outputs every N BSDebug Iterations) */


/*** Define input and output file names  *********************/
char   fin[FLENGTH];           /* input temp file  */
char   fout[FLENGTH];          /* output temp file */
char   fout2[FLENGTH];
char   plotfilename[FLENGTH];  /* file to pass to GnuPlot */
char   *Parm_name[]={"alpha", "beta", "theta1", "theta2","rho",
"phi1", "phi2", "phi3", "phi4", "phi5",
"phi6", "phi7", "phi8", "phi9", "phi10"};
char   *anatxt[]={"Full model", "Fitted model", "Reduced model"};
char    fname2[FLENGTH], logfile[FLENGTH], *dot2;

#ifdef LOGGING_ON
FILE    *fp_log;           /*  log file if requested */
#endif

/****** Variables that will not be changed except Spec  *******/
int    *Spec;    /* vector used to identify user input parm. */
int    *IniSp;
double *Yp;      /* positive dependent variable data array   */
double *Yn;      /* negative dependent variable data array   */
double *Xi;      /* independent variable data array          */
double *ScXi;    /* scaled dose levels                       */
double *Ypp;     /* predicted positive dependent values      */
double *Ep;      /* estimated probability                    */
double *Ls;      /* litter-specific covariate                */
int    *Xg;
double *SR; 			/*scaled residual */
double *Rlevel;
double *Bmdl;
double *IniP;
int    Nobs, nparm, ngrp, restrict, initial, appendix, smooth, fixedSize,
bmdlCurve;
double xmax, xmin, scale;
double maxdose;
double Maxloglik;
int    DeBuG = 0;
int    MxLkCnt = 0;
double smean, smax, smin, sijfixed;    /* smean: mean of litter size.
									   smax: maximum litter size-smean.
									   smin: minimun litter size-smean. */
FILE *fp_dbg;
int BSIter;      /*number of iterations for bootstrapping method */     
long BSSeed;     /*seed values for bootstrapping method - 0 defaults to time-generated*/


/* also used in O_func, _g, _vcv. */
double smean1, xmax, gamm0;            /* mean litter size in group */

/** changing variable **/
int    replace, brat;
double tD,  BMD_lk, LR, ck, upb=18, BMR;
int ErrorFlag; /* Error States from DMNGB */


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
** main--main function used to call NCTR mode fitting program.
Includes: biosubcc.c--common subfunction C program.
*****************************************************************/
int main (int argc, char *argv[])
{


	int     iter,i, j, junk;           /*iteration variable*/
	int     nvar;                      /*number of non-fixed parameters*/
	int     bmdose;                    /*flag for computing benchmark dose*/
	int     Nmiss;                     /*number of records with missing values*/
	int     nparm_known;               /*number of specified parameters */
	double  BMD, lkf,lkr,xlk,W, junk1; /*log likelihoods */
	double  *Parms;                    /*parameter array*/
	int     *bounded;                  /*array containing 0 or 1 depending if the parm is hit a bound */
	VarList *varsum;                   /*info for variables--p. dep.,n. dep., indep.*/
	AnaList *anasum;                   /*information for ANONA analysis*/
	double  **vcv;                     /*variance and covariance matrix*/
	double  back, back1;
	double  *GYp, *GYn, *GXi;
	double  maxsize;
	char    model_name[MNLENGTH], user_note[UNLENGTH], junkname[FLENGTH];
	char    dose_name[CNLENGTH], posi_name[CNLENGTH], nega_name[CNLENGTH], junkname1[CNLENGTH], junkname2[FLENGTH];
	time_t  ltime;
	char long_path_name[FLENGTH];

	/* Set time zone from TZ environment variable. If TZ is not set,
	* the operating system is queried to obtain the default value
	* for the variable.
	*/
	/* _tzset(); */
	time( &ltime );

	/* The following are arguments for D1MACH */

	//Min_increment = DBL_EPSILON;
	//Max_double = DBL_MAX;

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

	/*if (argc >= 2)
	path_name(argv[1]);*/

	if (argc > 2)
	{
		path_name2(argc, argv, long_path_name);
		argv[1] = long_path_name;
	}


	fp_in=fopen(argv[1], "r");
	if (fp_in==NULL)
	{
		fprintf(stderr, "Error in opening input file.\n");
		fprintf(stderr, "...Exited to system!\n");
		exit (1);
	}

#ifdef LOGGING_ON
	/* open the log file*/
	strcpy(logfile,argv[1]);
	dot2 = strchr(logfile, (int) '.');
	(*dot2) = (char) 0;
	strcpy(fname2,logfile);
	strcat(logfile,"-nctr.log");
	fp_log = fopen(logfile, "w");
	if (fp_log == (FILE *) NULL)
	{
		ERRORPRT("Unable to open log for Nctr.C.");
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
	fprintf(fp_log,"model_name: %s\n",model_name);
	fprintf(fp_log,"argc = %d\n",argc);
	fprintf(fp_log,"argv[1] = %s\n\n",argv[1]);
	fprintf(fp_log,"INPUT VALUES FOR DATA SET\n\n");
	fprintf(fp_log,"Nobs = %4d                 (Number of Observations)\n",Nobs);
	fprintf(fp_log,"ngrp = %4d                 (Total Number of Groups)\n",ngrp);
	fprintf(fp_log,"nparm = %3d                 (The Number of Parameters in the Model)\n",nparm);
	fflush(fp_log);
#endif

	/*allocate memory for arrays*/
	Parms   = DVECTOR(1, nparm);
	bounded = IVECTOR(1, nparm);
	IniP    = DVECTOR(1, nparm);
	IniSp   = IVECTOR(1, nparm);
	Spec    = IVECTOR(1, nparm);
	Xi      = DVECTOR(1, Nobs);
	ScXi    = DVECTOR(1, Nobs);
	Yp      = DVECTOR(1, Nobs);
	Yn      = DVECTOR(1, Nobs);
	Ls      = DVECTOR(1, Nobs);
        SR      = DVECTOR(1, Nobs);
	Xg      = IVECTOR(1, Nobs);
	varsum  = VLVECTOR(1, 3);
	anasum  = ALVECTOR(1, 3);
	Rlevel  = DVECTOR(1, 5);
	Bmdl    = DVECTOR(1, 5);
	vcv     = DMATRIX (1,nparm,1,nparm);

	junk = 0;    /*Used to see if an extension was added to output file name*/

	/* get filenames */
	Get_Names(argv[1], fout, fout2, plotfilename);

	/************* input All the data from  .(d) file ************/
	fscanf(fp_in,"%d%lf%lf%d%d%d%d%d%d",&ITMAX, &Rel_Conv, &Parm_Conv, &bmdlCurve, &restrict, &bmdose,
		&fixedSize,&appendix, &smooth);

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
		)
	{
		printf("Error in opening  output files.\n");
		printf ("...now exiting to system...\n");
		fprintf(fp_out,"Error in opening output files.\n");
		fprintf (fp_out,"...Exited to system!\n");
		exit (1);
	}

	fscanf(fp_in,"%lf%d%lf%d%ld",&bmdparm.effect,&bmdparm.risk,&bmdparm.level, &BSIter, &BSSeed);

        if (BSSeed == 0)              /* Set seed from time clock if default(BSSeed=0) is specified */
          {
          BSSeed = time (NULL);
          }

#ifdef LOGGING_ON
	fprintf(fp_log,"ITMAX = %4d                (Maximum number of iterations)\n",ITMAX);
	fprintf(fp_log,"Rel_Conv = %g     (Rel Function Convergence, default=2.22045e-16)\n",Rel_Conv);
	fprintf(fp_log,"Parm_Conv = %g    (Parameter Convergence)\n",Parm_Conv);
	fprintf(fp_log,"bmdlCurve = %1d               (BMDL Curve Calculation; 1=yes, 0=no)\n",bmdlCurve);
	fprintf(fp_log,"restrict = %1d                (Restriction on nctr coefficients)\n",restrict);
	fprintf(fp_log,"bmdose = %1d                  (BMD Calculation; 1=yes, 0=no)\n",bmdose);
	fprintf(fp_log,"fixedSize = %1d               (Fixed Litter Size; 1=ctrl group mean, 0=overall mean)\n",fixedSize);
	fprintf(fp_log,"appendix = %1d                (Append or Overwrite output file)\n",appendix);
	fprintf(fp_log,"smooth = %2d                 (Smooth Option)\n",smooth);
	fprintf(fp_log,"bmdparm.effect = %4g       (BMR Factor)\n",bmdparm.effect);
	fprintf(fp_log,"bmdparm.risk = %1d            (BMR Type; 0=extra, 1=added)\n",bmdparm.risk);
	fprintf(fp_log,"bmdparm.level = %4g        (Confidence level)\n",bmdparm.level);
#endif

	/* Print model and file information on output page */
	Output_Header(Version_no, argv[1], plotfilename, ctime(&ltime), user_note);

	if (bmdose < 0 || bmdose > 1)
	{
		FreeUp_mem(Parms, varsum, anasum, vcv, GXi, GYp, GYn, bounded);
		ERRORPRT("Error in choosing bemchmark dose computation.");
	}

	/*obtain user input parameters*/
	READ_PARAMETERS(nparm,Parms);
	SWAP(Parms[4], Parms[5], junk1);
	SWAP(Parms[3], Parms[5], junk1);
	SWAP(Parms[2], Parms[5], junk1);
	FILL_SPECVECTOR(nparm,Parms,Spec);
	nparm_known = COUNT_SPECVECTOR(nparm,Spec);


	/******  obtain user input initial parameters values      *******/
	fscanf(fp_in,"%d", &initial);
	READ_PARAMETERS(nparm,IniP);
	SWAP(IniP[4], IniP[5], junk1);
	SWAP(IniP[3], IniP[5], junk1);
	SWAP(IniP[2], IniP[5], junk1);
	FILL_SPECVECTOR(nparm,IniP,IniSp);
	for(i = 1; i <= nparm; i++)
	{
		if(Spec[i] == 1)
			IniP[i] = 1;
	}
	
#ifdef LOGGING_ON
	fprintf(fp_log,"\nThe specified or default parameters are:\n");
	fprintf(fp_log,"                                           Spec Value\n");
	fprintf(fp_log,"Parms[ 1] = %12.6g        (alpha)        %d\n",Parms[1],Spec[1]);
	fprintf(fp_log,"Parms[ 2] = %12.6g        (rho)          %d\n",Parms[2],Spec[2]);
	fprintf(fp_log,"Parms[ 3] = %12.6g        (theta1)       %d\n",Parms[3],Spec[3]);
	fprintf(fp_log,"Parms[ 4] = %12.6g        (theta2)       %d\n",Parms[4],Spec[4]);
	fprintf(fp_log,"Parms[ 5] = %12.6g        (rho)          %d\n",Parms[5],Spec[5]);
	for (i=1; i<=ngrp; i++)
		fprintf(fp_log,"Parms[%2d] = %12.6g        (phi%d)         %d\n",i+5,Parms[i],i,Spec[i]);
	fprintf(fp_log,"\ninitial = %d     (Number of initialized parameters)\n",initial);
	fprintf(fp_log,"The initial parameter values are:\n");
	fprintf(fp_log,"                                           IniSp Value\n");
	fprintf(fp_log,"IniP[ 1] = %12.6g         (alpha)         %d\n",IniP[1],IniSp[1]);
	fprintf(fp_log,"IniP[ 2] = %12.6g         (rho)           %d\n",IniP[2],IniSp[2]);
	fprintf(fp_log,"IniP[ 3] = %12.6g         (theat1)        %d\n",IniP[3],IniSp[3]);
	fprintf(fp_log,"IniP[ 4] = %12.6g         (theta2)        %d\n",IniP[4],IniSp[4]);
	fprintf(fp_log,"IniP[ 5] = %12.6g         (rho)           %d\n",IniP[5],IniSp[5]);
	for (i=1; i<=ngrp; i++)
		fprintf(fp_log,"IniP[%2d] = %12.6g         (phi%d)          %d\n",i+5,IniP[i],i,IniSp[i]);
	fprintf(fp_log,"The Spec and IniSp vectors are just flags to tell if the user gave\n");
	fprintf(fp_log,"a value for the specific parameter.\n");
#endif
	/*obtain observation data into Yp, Yn, Xi, Ls, Xg vectors*/
	fscanf(fp_in,"%s%s%s%s%s", dose_name, posi_name, nega_name, junkname1, junkname2);
	Nmiss = READ_OBSDATA5V(Nobs,Xi,Yp,Yn,Ls,Xg);

	Sort_4_By_Dose(Nobs, Xi, Yn, Yp, Ls);   /* Sort the arrays by dose and move around
											the other arrays appropriately */

	Nobs -= Nmiss;             /* extern variable Nobs has been changed. */

	/* Create ScXi array */
	maxdose = Xi[Nobs];
	for (i = 1; i <= Nobs; i++)
	{
		ScXi[i] = Xi[i]/maxdose;
	}

	junk1 = Xi[1];
	junk = 1;

	for(i = 1; i <= Nobs; i++)
	{
		if(Xi[i] == junk1)			 /* This little loop gets the "group" values */
			Xg[i] = junk;				 /* without explicitly entering them at the  */
		else					 /* spreadsheet level.                       */
		{
			Xg[i] = ++junk;
			junk1 = Xi[i];
		}
	}

	ngrp = junk;
	GXi = DVECTOR(1, ngrp);
	GYp = DVECTOR(1, ngrp);
	GYn = DVECTOR(1, ngrp);

#ifdef LOGGING_ON
	fprintf(fp_log,"\nDose              Pos.Resp.            Neg.Resp.            Lit-Spec. Cov.       Group\n");
	for (i=1; i<=Nobs; i++)
		fprintf(fp_log,"Xi[%3d]=%4g      Yp[%3d]=%4g         Yn[%3d]=%4g         Ls[%3d]=%4g         Xg[%3d]=%3d\n",
		i,Xi[i],i,Yp[i],i,Yn[i],i,Ls[i],i,Xg[i]);
	fprintf(fp_log,"\nNmiss = %3d             (Number of missing values)\n",Nmiss);
	fprintf(fp_log,"Nobs = %4d             (Number of observations)\n",Nobs);
	fflush(fp_log);
#endif

	/*********** end of input data. **************/
	nparm_known = COUNT_SPECVECTOR(nparm, Spec);

	if (Nobs < (nparm-nparm_known))
	{
		FreeUp_mem(Parms, varsum, anasum, vcv, GXi, GYp, GYn, bounded);
		ERRORPRT("Observation # < parameter # for NCTR model.");
	}

	/*output title and summary of intput data  ****************************/
	OUTPUT_TEXT("\n The probability function is: ");
	OUTPUT_TEXT("\n\n Prob. = 1 - exp[-(alpha + th1*Rij) - (beta + th2*Rij)*Dose^rho],");
	OUTPUT_TEXT("\n          where Rij is the centralized litter specific covariate.");
	/*indicates restrict power*/
	gamm0=0;
	if (restrict==Yes)
	{
		gamm0=1.0;
		OUTPUT_TEXT("\n Restrict Power rho >= 1. ");
	}

	fprintf (fp_out, "\n\n\n Total number of observations = %d",Nobs+Nmiss);
	fprintf (fp_out, "\n Total number of records with missing values = %d",Nmiss);
	fprintf (fp_out, "\n Total number of parameters in model = %d", nparm);
	fprintf (fp_out, "\n Total number of specified parameters = %d\n\n", nparm_known);
	fprintf (fp_out, "\n Maximum number of iterations = %d\n", ITMAX);
	fprintf (fp_out, " Relative Function Convergence has been set to: %g\n", Rel_Conv);
	fprintf (fp_out, " Parameter Convergence has been set to: %g\n\n", Parm_Conv);
        fprintf (fp_out, " Number of Bootstrap Iterations per run: %d\n", BSIter);
        fprintf (fp_out, " Bootstrap Seed:  %ld\n\n", BSSeed);
	fflush(fp_out);

	if(nparm_known > 0)
	{
		fprintf (fp_out, " User specifies the following paramters:");
		for (i=1; i<= nparm; i++)
		{
			if(Spec[i] == 1)
				fprintf (fp_out, "\n %15s = %10.5g", Parm_name[i-1], Parms[i]);
		}
		fprintf (fp_out, "\n\n");
		fflush(fp_out);
	}

	/* compute init_lkf for full model and init_lkr for reduced model */
	lkf = 0.0;
	varsum[1].S = 0;
	varsum[2].S = 0;
	for (i=1;i<=Nobs;i++)
	{
		varsum[1].S += Yp[i];
		varsum[2].S += Yn[i];
		W = Yp[i] / (Yp[i]+Yn[i]);
		if (W > 0)   lkf += Yp[i] * log(W);
		if (W < 1)   lkf += Yn[i] * log(1- W);
	}
	W = varsum[1].S / (varsum[1].S + varsum[2].S);
	lkr = varsum[1].S * log(W) + varsum[2].S * log(1- W);

#ifdef LOGGING_ON
	fprintf(fp_log,"\nlkf = %4g          (Likelihood for full model)\n", lkf);
	fprintf(fp_log,"lkr = %4g          (Likelihood for reduced model)\n", lkr);
	fflush(fp_log);

	fprintf(fp_log, "\n********** Call to NCTR_Fit ****************\n");
	fprintf(fp_log, "Variables going in:\n");
	fprintf(fp_log, "nparm = %2d;  ngrp = %2d;  EPS = %6g;  iter = %4d;  xlk = %6g\n",
		nparm, ngrp, EPS, iter, xlk);
	for (i=1; i<=nparm; i++)
		fprintf(fp_log,"Parms[%2d] = %13g\n", i, Parms[i]);
	fflush(fp_log);
#endif
	/* fitting NCTR model and output parameter estimates */

	NCTR_fit(nparm,ngrp, Parms, EPS, &iter, &xlk);

#ifdef LOGGING_ON
	fprintf(fp_log, "Variables coming out:\n");
	fprintf(fp_log, "nparm = %2d;  ngrp = %2d;  EPS = %6g;  iter = %4d;  xlk = %6g\n",
		nparm, ngrp, EPS, iter, xlk);
	for (i=1; i<=nparm; i++)
		fprintf(fp_log,"Parms[%2d] = %13g\n", i, Parms[i]);
	fprintf(fp_log, "*********************************************\n");
	fflush(fp_log);
#endif

	/* Which parameters are on their bounds? */
	Which_Bounded (Spec, Parms, bounded);
#ifdef LOGGING_ON
	fprintf(fp_log, "\nWhich variables are bounded?\n");
	for (i=1; i<=nparm; i++)
		fprintf(fp_log,"Parms[%2d]  ->  bounded[%2d] = %2d\n", i, i, bounded[i]);
	fprintf(fp_log,"\n");
	fflush(fp_log);
#endif

	/* Compute the aprox. covariance matrix */
	INITIALIZE_DMATRIX(vcv, nparm, nparm);
	nvar = NCTR_vcv(nparm,Spec,Parms,vcv);
	nvar = Take_Out_Bounded_Parms(nvar, bounded, vcv);
	if (nvar > 0)
		INVMAT(vcv, nvar);
	OUTPUT_DTMS3PARMS(nparm,Spec,bounded,Parms,Parm_name,vcv,1);
	fflush(fp_out);
#ifndef RBMDS
	fflush(fp_out2);
#endif
	/*compute and output ANOVA table elements*/
	DTMS3ANOVA (nparm,Nobs,Spec,lkf,xlk,lkr,anasum,bounded);

	maxsize = Ls[1];
	for(i = 1; i <= ngrp; i++)
	{
		for( j = 1; j <= Xg[i]; j++)
		{
			if(Ls[i+j] > maxsize)
				maxsize = Ls[i+j];
		}
	}
	maxsize = maxsize - smean;

	TEMP_ANOVA_OUTPUT(anatxt,anasum); /* output ANOVA table */
	fflush(fp_out);

	N_Goodness (ngrp, nparm, Parms, bounded, Nobs, Xi, Yp, Yn,
                    Ls, Xg, SR);
	fflush(fp_out);

#ifndef RBMDS
	/********* output to **.002 ************/
	fprintf (fp_out2, "\n BMD_flag \t %d \n Nosb \t%d \n nparm \t%d",  bmdose, ngrp, 6 );
	fprintf (fp_out2, "\n  Con_lev \t%3.3g ", bmdparm.level);
	fprintf (fp_out2, "\n  RiskType \t%d ", bmdparm.risk);
	fprintf (fp_out2, "\n  Effect \t%3.3g ", bmdparm.effect);
	for (i=1;i<=5; i++) fprintf (fp_out2, "\n %s \t %5.3g", Parm_name[i-1], Parms[i]);
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
		if(Spec[2]==Yes)
		{
#ifndef RBMDS
			fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  No);
#endif
			fprintf (fp_out,"\n\n  %s parameter is fixed. ", Parm_name[1]);
			fprintf (fp_out,"The likelihood function can not be reparameterized in BMD.");
			fflush(fp_out);
#ifndef RBMDS
			fflush(fp_out2);
#endif
			exit(0);
		}

		back=1-exp(-Parms[1]-Parms[3]*sijfixed);
		back1=1-back;
		if (bmdparm.risk==1) back1=1;

#ifdef LOGGING_ON
		{
			fprintf(fp_log,"Before Call to NCTR_BMD\n");
			fprintf(fp_log,"back = %12.6g\n",back);
			fprintf(fp_log,"back1 = %12.6g\n",back1);
			for(i = 0; i <= nparm; i++)
				fprintf(fp_log,"\n\n***In main() Parms[%d]=%f\n", i, Parms[i]);
			fflush(fp_log);
		}
#endif


		NCTR_BMD (nparm, Parms, EPS, &junk, xlk, Rlevel, Bmdl, &BMD);
                SRoI(ngrp, Nobs, GXi, Xi, Xg, SR, Ls, sijfixed + smean, BMD);

#ifdef LOGGING_ON
		{
			fprintf(fp_log,"After Call to NCTR_BMD\n");
			fprintf(fp_log,"Bmdl[1] = %12.6g\n",Bmdl[1]);
			fprintf(fp_log,"BMD = %12.6g\n",BMD);
			fflush(fp_log);
		}
#endif
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
		for (i=1;i<=5;i++)
			fprintf (fp_out2,"\n %f %f", Bmdl[i], Rlevel[i]*back1+back);
#endif
	}
        N_Bootstrap (ngrp, nparm, Parms, bounded, Nobs, Xi, Yp, Yn,
                    Ls, Xg, SR, BSIter, BSSeed);
        if (bmdose==Yes)
        	{
                if(fixedSize==1)
			{
				fprintf(fp_out, "\n\nTo calculate the BMD and BMDL, the litter\n");
				fprintf(fp_out, "specific covariate is fixed at the mean litter\n");
				fprintf(fp_out, "specific covariate of control group: %f", sijfixed+smean);
			}
			else
			{
				fprintf(fp_out, "\n\nTo calculate the BMD and BMDL, the litter\n");
				fprintf(fp_out, "specific covariate is fixed at the overall mean\n");
				fprintf(fp_out, "of the litter specific covariates: %f", sijfixed+smean);
			}
		OUTPUT_BENCHMD(1,BMD);
#ifndef RBMDS
		fprintf(fp_out, "            BMDL =%15.6g\n\n", Bmdl[1]);
#else
		fprintf(fp_out, "            BMDL =%30.22g\n\n", Bmdl[1]);
#endif
		fflush(fp_out);
   		}

#ifndef RBMDS
	fprintf (fp_out2,"\n\n Check_result %d", Yes); /* indicate if all the required data has been output. */
#endif

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
	FREE_DVECTOR (GXi, 1, ngrp);
	FREE_DVECTOR (GYp, 1, ngrp);
	FREE_DVECTOR (GYn, 1, ngrp);

	/*close opened temp files*/
	CLOSE_FILES ();
#ifdef LOGGING_ON
	fflush(fp_log);
	fclose (fp_log);
#endif

	return(0);
} /*end of main*/


/*********************************************************************** */
/* Which_Bounded -- Fills the 1-based vector bounded with 1 if the */
/*                  corresponding parameter is on one of its boundaries or */
/*                  was fixed, 0 otherwise. */
/* */
/*                  Global: */
/*                    nparm, smax, Max_double, SlopeUpperBound*/
/* */
/*                  input: */
/*                    Spec  1-based vector giving which parameters are fixed */
/*                    Parms 1-based vector giving the full vector of */
/*                          parameter values */
/*                  output: */
/*                    bounded 1-based vector whose elements are 1 if */
/*                            the corresponding parameter is either */
/*                            fixed */
/*                            (i.e., Spec[i] == 1) or its value is on */
/*                            a boundary */
/*********************************************************************** */

void Which_Bounded (int Spec[], double Parms[], int bounded[])
{
	int i;
	double maxsize;

	for (i=1; i<=nparm; i++)
	{
		bounded[i] = Spec[i];
	} /* end for */
	maxsize = smax - smean;

	/* There are many parameter restrictions in the Nctr model.*/
	/* The if statements */
	/* below address the parameter space restrictions, and adjust the */
	/* degrees of freedom appropriately */

	/*               _         */
	/*alpha + Theta1*Sij >= 0: */
	if ((Spec[1] != 1 || Spec[3] != 1) && VERYCLOSE(Parms[1], -Parms[3] * maxsize))
	{
		if (Spec[1] == 1) bounded[3] = 1;
		if (Spec[3] == 1) bounded[1] = 1;
	}

	/*              _         */
	/*beta + Theta2*Sij >= 0: */
	if ((Spec[2] != 1 || Spec[4] != 1) && VERYCLOSE(Parms[2], -Parms[4] * maxsize))
	{
		if (Spec[2] == 1) bounded[4] = 1;
		if (Spec[4] == 1) bounded[2] = 1;
	}


	//Rho >= 1 when restricted

	if (restrict == 1 && Spec[5] != 1 && VERYCLOSE(Parms[5], 1))
	{
		bounded[5] = 1;
	}

	//Rho >= 0:

	if (Spec[5] != 1 && VERYCLOSE(Parms[5], 0))
	{
		bounded[5] = 1;
	}

	//beta > = 0:

	if (Spec[2] != 1 && VERYCLOSE(Parms[2], 0))
	{
		bounded[2] = 1;
	}

	//Phi1,..., Phig >= 0:

	for (i = 6; i <= nparm; i++)
	{
		if (Spec[i] != 1 && VERYCLOSE(Parms[i], 0))
		{
			bounded[i] = 1;
		}
	}
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


/*******************************************************************
**NCTR_probs -- computes litter-specific probabilities of an
*               adverse response, and, optionally, the gradient
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
void NCTR_probs (int nobs, double doses[], double sij[], double probs[],
				 int nparm, double p[], int compgrad, double **gradij)
{
	int i, plus5;
	double *pint, ex, ex1, ex2, ex3, ex4;
	double pwx, lamb, ck;

	pint = DVECTOR(1, nparm);

	for (i = 1; i <= nparm; i++)
	{
		pint[i] = p[i];
	}

	for (i = 1; i <= nobs; i++)
	{
		sij[i] = sij[i] - smean;
	}

	if(replace == Yes) /* We're computing a BMDL */
	{
		if (bmdparm.risk==1)
			pint[2] = -log( 1-BMR*exp(pint[1]*(1+pint[3]*sijfixed)))/(1+pint[4]*sijfixed)/pow(tD,pint[5]);
		else
			pint[2]= -log(1-BMR)/(1+pint[4]*sijfixed)/pow(tD, pint[5]);
	}

	for (i = 1; i <= nobs; i++)
	{

		if (doses[i] > 0)

			pwx = pow(doses[i],pint[5]);
		else
			pwx = 0.0;
		ex1 = pint[1]*(1+pint[3]*sij[i]);
		ex2 = pint[2]*(1+pint[4]*sij[i]);
		ex = 1.0 - exp(-(ex1+(ex2*pwx)));
		PROBABILITY_INRANGE (&ex);
		probs[i] = ex;

		if (compgrad == 1)  /* compute partial derivatives */
		{
			lamb = 0.0;
			if (bmdparm.risk == 1)
				ck = -1.0*log((1.0-BMR*lamb));
			else
				ck = -log(1-BMR);
			plus5=5 + Xg[i];

			if (replace == Yes)  /* bmdl calculation */
			{
				ex3 = (1.0-ex);
				ex4 = pow(tD, -pint[5])*( lamb*BMR*(1+pint[2]*sijfixed)/(1-lamb*BMR) );
				gradij[i][2] = (ex3*(1+pint[4]*sij[i])*pwx);
				gradij[i][1] = (ex3*(1+pint[3]*sij[i])+ gradij[i][2]*ex4);
				gradij[i][3] = (ex3*pint[1]*sij[i] + gradij[i][2]*ex4*pint[1]*sijfixed);
				gradij[i][4] = (ex3*pint[2]*sij[i]*pwx - gradij[i][2]*ck*pow(tD,-pint[5])*sijfixed/(1+pint[4]*sijfixed)/(1+pint[4]*sijfixed) );
				if ( doses[i] > 0.0)
					gradij[i][5] = (ex3*ex2*log(doses[i])*pwx - gradij[i][2]*ck/(1+pint[4]*sijfixed) *log(tD)*pow(tD, -pint[5]) );
				else
					gradij[i][5]= 0.0-gradij[i][2]*ck/(1+pint[4]*sijfixed) *log(tD)*pow(tD, -pint[5]);
			}

			else /* replace == No */
			{
				ex3 = (1.0-ex);
				gradij[i][1] = ex3*(1+pint[3]*sij[i]);
				gradij[i][2] = ex3*(1+pint[4]*sij[i])*pwx;
				gradij[i][3] = ex3*sij[i]*pint[1];
				gradij[i][4] = ex3*sij[i]*pint[2]*pwx;
				if (doses[i] > 0)
					gradij[i][5] = ex3*ex2*log(doses[i])*pwx;
				else
					gradij[i][5] = 0.0;

			}  /* end of replace == No */
		}
	}

	/* Change back to Ls array */
	for (i = 1; i <= nobs; i++)
	{
		sij[i] = sij[i] + smean;
	}

	FREE_DVECTOR(pint, 1, nparm);

}


/*******************************************************************
**NCTR_lk -- used to compute the log likelihood for NCTR model.
* 		     Extern var.: smean, smax, Nobs, Xi, Yp, Yn, Ls, Xg.
*
*********************************************************************/
void NCTR_lk(long int *nvar, double *x, long int *nf, double *f,
			 long int *uiparm, double *urparm, void (*ufparm)())
{
	double  xlk;          		  /* log likelihood. */
	int     i, j, k, plus5, jfixed, jvar;
	double  tm1,tm2,tm3,tm;   	          /* temp var. */
	double  *p;				  /* for "untransform" parms. */
	double  *probs;                         /* litter-specific probabilities */
	double  **gradij;                       /* not used right here, but needed for function call */
	int     compgrad;
#ifdef LOGGING_ON
	int     run=-9999;                      /* number of the run that we wish output in the log file */
#endif

#ifdef LOGGING_ON
	if (MxLkCnt == run)
	{
		fprintf(fp_log, "\nFunction: NCTR_lk (Beginning)\n");
		fprintf(fp_log, "The fixed parameters are:\n");
		for (i = 0; i <= (nparm-*nvar-1); i++)
			fprintf(fp_log, "urparm[%2d] = %12.6g   %p\n", i, urparm[i], &urparm[i]);
		fprintf(fp_log, "The varying parameters are:\n");
		for (i = 0; i <= (*nvar-1); i++)
			fprintf(fp_log, "x[%2d] = %12.6g   %p\n", i, x[i], &x[i]);
	}
#endif
	p=DVECTOR(1, nparm);
	probs = DVECTOR(1, Nobs);
	gradij = DMATRIX(1, Nobs, 1, 5);

	jfixed = jvar = 0;
	for(j=1; j<=nparm; j++)  /* reconstruct the parameter vector. */
	{
		if(Spec[j]==Yes)
		{
			p[j]=urparm[jfixed];
			jfixed++;
		}
		else
		{
			p[j]=x[jvar];
			jvar++;
		}
	}

	compgrad = 0;

#ifdef LOGGING_ON
	if (MxLkCnt == run)
	{
		fprintf(fp_log, "\nFunction: NCTR_lk (Before NCTR_probs)\n");
		fprintf(fp_log, "The p array is:\n");
		for (i = 1; i <= nparm; i++)
			fprintf(fp_log, "p[%2d] = %12.6g   %p\n", i, p[i], &p[i]);
		fprintf(fp_log, "The fixed parameters are:\n");
		for (i = 0; i <= (nparm-*nvar-1); i++)
			fprintf(fp_log, "urparm[%2d] = %12.6g   %p\n", i, urparm[i], &urparm[i]);
		fprintf(fp_log, "The varying parameters are:\n");
		for (i = 0; i <= (*nvar-1); i++)
			fprintf(fp_log, "x[%2d] = %12.6g   %p\n", i, x[i], &x[i]);
	}
#endif
	NCTR_probs(Nobs, ScXi, Ls, probs, nparm, p, compgrad, gradij);

#ifdef LOGGING_ON
	if (MxLkCnt == run)
	{
		fprintf(fp_log, "\nFunction: NCTR_lk (After NCTR_probs)\n");
		fprintf(fp_log, "The fixed parameters are:\n");
		for (i = 0; i <= (nparm-*nvar-1); i++)
			fprintf(fp_log, "urparm[%2d] = %12.6g   %p\n", i, urparm[i], &urparm[i]);
		fprintf(fp_log, "The varying parameters are:\n");
		for (i = 0; i <= (*nvar-1); i++)
			fprintf(fp_log, "x[%2d] = %12.6g   %p\n", i, x[i], &x[i]);
	}
#endif

	xlk = 0.0;
	for (i=1; i<=Nobs; i++)
	{
		tm1 = 0.0;
		tm2 = 0.0;
		tm3 = 0.0;
		plus5 = 5 + Xg[i];

		j = (int) Yp[i];
		if ((probs[i] == 0.0) && (j > 0))
		{
			tm1 -= 40.0;
		}
		else
		{
			for (k = 1; k <= j; k++)
			{
				tm = probs[i] + (k - 1)*p[plus5];
				if (tm == 0.0) tm1 += Max_double;
				else tm1 += log(tm);
			}
		}

		j = (int) Yn[i];
		if ((probs[i] >= 1.0) && (j > 0))
		{
			tm2 -= 40.0;
		}
		else
		{
			for (k = 1; k <= j; k++)
			{
				tm = 1.0 - probs[i] + (k - 1) * p[plus5];
				if (tm == 0.0) tm2 += Max_double;
				else tm2 += log (tm);
			}
		}

		j = (int) (Yn[i] + Yp[i]);
		for (k = 1; k <= j; k++)
		{
			tm = 1.0 + (k - 1) * p[plus5];
			if (tm == 0.0) tm3 += Max_double;
			else tm3 += log (tm);
		}

		xlk += (tm1 + tm2 - tm3);
	}

	FREE_DVECTOR(p, 1, nparm);
	FREE_DVECTOR(probs, 1, Nobs);
	FREE_DMATRIX(gradij, 1, Nobs, 1, 5);
	*f = -xlk;

#ifdef LOGGING_ON
	if (MxLkCnt == run)
	{
		fprintf(fp_log, "MxLkCnt = %3d     (Counts Calls to MAX_lk)\n", MxLkCnt);
		fprintf(fp_log, "smean = %12.6g\n", smean);
		fprintf(fp_log, "\nFunction: NCTR_lk (Exiting)\n");
		fprintf(fp_log, "The fixed parameters are:\n");
		for (i = 0; i <= (nparm-*nvar-1); i++)
			fprintf(fp_log, "urparm[%2d] = %12.6g   %p\n", i, urparm[i], &urparm[i]);
		fprintf(fp_log, "The varying parameters are:\n");
		for (i = 0; i <= (*nvar-1); i++)
			fprintf(fp_log, "x[%2d] = %12.6g   %p\n", i, x[i], &(x[i]));
		fprintf(fp_log, "The likelihood value for these is: %12.6g\n", -xlk);
	}
#endif
}


/*******************************************************************
* NCTR_g -- used to compute the gradients for NCTR model.
*		Extern var.: smean, smax, Nobs, Xi, Yp, Yn, Ls, Xg.
*
**********************************************************************/
void NCTR_g(long int *nvar, double *x, long int *nf, double *g,
			long int *uiparm, double *urparm, void (*ufparm)())
{
	double  ex;                                   /* temp var. */
	double  tm1, tm2, tm3, tm1a, tm2a, tm3a, tm, tm12;  /* temp var. */
	double  *dd, *p, *probs, **gradij, *tmp_g;
	int     i, j, k, plus5, jfixed, jvar, compgrad;
#ifdef LOGGING_ON
	int     run=-9999;                          /* number of the run that we wish output in the log file */
#endif

#ifdef LOGGING_ON
	if (MxLkCnt == run)
	{
		fprintf(fp_log, "\nFunction: NCTR_g (Beginning)\n");
		fprintf(fp_log, "The fixed parameters are:\n");
		for (i = 0; i <= (nparm-*nvar-1); i++)
			fprintf(fp_log, "urparm[%2d] = %12.6g\n", i, urparm[i]);
		fprintf(fp_log, "The varying parameters are:\n");
		for (i = 0; i <= (*nvar-1); i++)
			fprintf(fp_log, "x[%2d] = %12.6g\n", i, x[i]);
	}
#endif
	p = DVECTOR(1, nparm);
	probs = DVECTOR (1, Nobs);
	gradij = DMATRIX(1, Nobs, 1, 5);
	dd = DVECTOR(1, nparm);
	tmp_g = DVECTOR(1, nparm);

	jfixed=jvar=0;

	for(j=1; j<=nparm; j++)  /* reconstruct the parmater vector. */
	{
		if(Spec[j]==Yes)
		{
			p[j]=urparm[jfixed];
			jfixed++;
		}
		else
		{
			p[j]=x[jvar];
			jvar++;
		}
	}

	compgrad = 1;
	NCTR_probs(Nobs, ScXi, Ls, probs, nparm, p, compgrad, gradij);


	/**********    initial g[j]'s            **************/
	for (j=1;j<=nparm;j++)
	{
		tmp_g[j]=0.0;
		g[j] = 0.0;
	}
	for (i=1;i<=Nobs;i++)
	{
		ex = probs[i];
		/*Compute first partial derivatives*/
		tm1 = 0.0;
		tm2 = 0.0;
		tm3 = 0.0;
		tm1a = 0.0;
		tm2a = 0.0;
		tm3a = 0.0;
		for(j=6;j<=nparm;j++)
			dd[j] = 0.0;
		plus5=5 + Xg[i];
		j = (int)Yp[i];
		if (ex > 0.0)
		{
			for (k = 1; k <= j; k++)
			{
				tm = ex + (k - 1)*p[plus5];
				tm1 += 1.0/tm;
				tm1a += (1.0/tm)*(k - 1);
			}
		}

		j = (int)Yn[i];
		if (ex < 1.0)
		{
			for (k = 1; k <= j; k++)
			{
				tm = 1.0 - ex + (k - 1) * p[plus5];
				tm2 += 1.0 / tm;
				tm2a += (1.0 / tm) * (k - 1);
			}
		}

		j = (int) (Yn[i] + Yp[i]);
		for (k = 1; k <= j; k++)
		{
			tm = 1.0 + (k - 1) * p[plus5];
			if (tm == 0.0)
				tm = 0.000001;
			tm3 += 1.0 / tm;
			tm3a += (1.0 / tm) * (k - 1);
		}

		tm12 = (tm1 - tm2);
		for (j = 1; j <= 5; j++)
		{
			dd[j] = gradij[i][j] * tm12;
		}

		dd[plus5] = (tm1a + tm2a - tm3a);
		for (j=1;j<=nparm;j++)
		{
			tmp_g[j] -= dd[j];
		}
	}
	/* end of first partial deri. */
	jvar=0;
	for(j=1; j<=nparm; j++)  /* reconstruct the parameter vector. */
		if(Spec[j]==No)
		{
			g[jvar]=tmp_g[j];
			jvar++;
		}


#ifdef LOGGING_ON
		if (MxLkCnt == run)
		{
			fprintf(fp_log, "smean = %12.6g\n", smean);
			fprintf(fp_log, "\nFunction: NCTR_g (Exiting)\n");
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
		FREE_DVECTOR(p, 1,nparm);
		FREE_DVECTOR(probs, 1, Nobs);
		FREE_DMATRIX(gradij, 1, Nobs, 1, 5);
		FREE_DVECTOR(dd, 1, nparm);
		FREE_DVECTOR(tmp_g, 1, nparm);
		return;

}

/*******************************************************************
*  NCTR_grad -- Computes the gradient of the NCTR likelihood function
*  with respect to the user form of the parameters.  This is to be used
*  in the NCTR_vcv, to compute a finite difference approximation to the
*  hessian of the the likelihood function
*
*******************************************************************/
void NCTR_grad(int nparm, int Spec[], double ptf[], double grad[])
{
	double  x, ex;                           /* temp var. */
	double  ex1,ex2, ex3, pwx;          /* temp var. */
	double  tm1,tm2,tm3,tm1a,tm2a,tm3a,tm;   /* temp var. */
	double  *dd, *p;
	int     i,j,k, plus5;

	dd=DVECTOR(1,nparm);
	p=DVECTOR(1,nparm);

	for (j=1; j<=nparm; j++)
		{
		p[j]=ptf[j];  /* copy the ptf[] because it should not be changed. */
		grad[j] = 0.0; /* initialize grad to zero */
		}


	for(j=6; j<=nparm; j++) p[j]=p[j]/(1-p[j]);     /* Phi --> Psi. */
	for(j=3;j<=4;j++)
	{
		if( p[j-2]>0) p[j]= p[j]/p[j-2];
		else p[j]=0;
	}

	for (i=1;i<=Nobs;i++)
	{
		/*Compute first partial derivatives*/
		x = Xi[i];
		pwx = pow(x, p[5]);
		ex1 = p[1]*(1+p[3]*(Ls[i]-smean));
		ex2 = p[2]*(1+p[4]*(Ls[i]-smean));
		ex = 1.0 - exp(-(ex1+(ex2*pwx)));
		PROBABILITY_INRANGE(&ex);
		tm1 = 0.0;
		tm2 = 0.0;
		tm3 = 0.0;
		tm1a = 0.0;
		tm2a = 0.0;
		tm3a = 0.0;
		for (j=6;j<=nparm; j++) dd[j]=0.0;
		j = (int)Yp[i];
		plus5=5+Xg[i];
		for (k=1;k<=j;k++)
		{
			tm = ex + (k-1)*p[plus5];
			if (tm == 0.0)
			{
				Warning ("Warning: Failed in computing information matrix. Can not obtain the Std. Err.");
				FREE_DVECTOR(dd, 1,nparm);
				return;
			}
			tm1 += 1.0/tm;
			tm1a += (1.0/tm)*(k-1);
		}
		j = (int)Yn[i];
		for (k=1;k<=j;k++)
		{
			tm = 1.0 - ex + (k-1)*p[plus5];
			if (tm == 0.0)
			{
				Warning ("Warning: Failed in computing information matrix. Can not obtain the Std. Err.");
				FREE_DVECTOR(dd, 1,nparm);
				return;
			}
			tm2 += 1.0/tm;
			tm2a += (1.0/tm)*(k-1);
		}
		j = (int)(Yp[i]+Yn[i]);
		for (k=1;k<=j;k++)
		{
			tm = 1.0 + (k-1)*p[plus5];
			if (tm == 0.0)
			{
				Warning ("Warning: Failed in computing information matrix. Can not obtain the Std. Err.");
				FREE_DVECTOR(dd, 1,nparm);
				return;
			}
			tm3 += 1.0/tm;
			tm3a += (1.0/tm)*(k-1);
		}

		ex3 = (1.0-ex)*(tm1-tm2);
		dd[1] = ex3;
		dd[2] = ex3*pwx;
		dd[3] = ex3*(Ls[i]-smean);
		dd[4] = ex3*(Ls[i]-smean)*pwx;
		if ( x <=0.0)  dd[5]=0.0;
		else dd[5] = ex3*ex2*log(x)*pwx;

		dd[plus5] = (tm1a + tm2a - tm3a); //*(1+p[plus5])*(1+p[plus5]);
		
		for (j=1;j<=nparm;j++)
		{
			grad[j] += dd[j];
		}

	}

	FREE_DVECTOR(p, 1, nparm);
	FREE_DVECTOR(dd, 1, nparm);
}


/*******************************************************************
**NCTR_vcv -- used to compute the vcv for NCTR model.
*	      Extern var.: smean, smax, Nobs, Xi, Yp, Yn, Ls, Xg.
*
******************************************************************/
int NCTR_vcv(int nparm, int Spec[], double ptf[], double **vcv)
{
	//double  x, ex;                           /* temp var. */
	//double  ex1,ex2, ex3, pwx, nij;          /* temp var. */
	//double  tm1,tm2,tm3,tm1a,tm2a,tm3a,tm;   /* temp var. */
	//double  *dd, *p;
	//int     i,j,k, plus5, indx_j, indx_k;


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
			NCTR_grad(nparm, Spec, saveparms, gradp);
			saveparms[i] = ptemp[i] - h[i];
			NCTR_grad(nparm, Spec, saveparms, gradm);
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

}



/**************************************************************
*MAX_lk -- used to obtain the Maximum log-lilikehood as well as
*           the estimatiors of parameters, given initial p[1..n],
*           object func. , and gradient func. G_func.
*
**************************************************************/
void MAX_lk(int nparm, double p[], double gtol, int *iter, double *fret)
{
	int i, jfixed, jvar;
	long int nvar;
	long int *uiparm;

	double *start;
	double *urparm;
	double *lower, *upper;
	long int dummy;
	void (*ufparm)();

	/* Set up initial parameter array, start.  All parameters go either
	into start (if they are changing to improve the fit) or urparm
	(if they are fixed). */

	nvar =  nparm;
	for (i = 1; i <= nparm; i++)
		nvar = nvar - Spec[i]; /* Count the varying parameters */

	uiparm = LIVECTOR(0,nparm-1);
	urparm = DVECTOR(0, nparm-nvar-1);
	if (nvar > 0)
	{
		start = DVECTOR(0, nvar-1);
		lower = DVECTOR(0, nvar-1);
		upper = DVECTOR(0, nvar-1);
	}
	else
	{
		start = DVECTOR(0,1);
		for(i=1;i<=nparm;i++) urparm[i-1]=p[i];
		NCTR_lk( &nvar, start,  &dummy, fret,
			uiparm, urparm, ufparm);
		*fret= -(*fret);
		ErrorFlag = -1;
		FREE_DVECTOR(start,0,1);
		FREE_LIVECTOR(uiparm,0,nparm-1);
		FREE_DVECTOR(urparm,0,nparm-nvar-1);
		return;
	}

	jfixed = 0;
	jvar = 0;
	for (i = 1; i <= nparm; i++) /* separate the fixed and variable parameters */
	{
		if (Spec[i] == 1)
		{
			urparm[jfixed] = p[i];
			jfixed++;
		}
		else
		{
			start[jvar]= p[i];
			jvar++;
		}
	}

	/* set up the bounds.  Each parameter is unique, here. */
	jvar = 0;
	/* alpha */
	if (Spec[1] != 1)
	{
		lower[jvar] = 0.0;
		upper[jvar] = Max_double ;
		jvar++;
	}
	/* beta */
	if (Spec[2] != 1)
	{
		lower[jvar] = 0.0;
		upper[jvar] = Max_double;
		jvar++;
	}
	/* theta1/alpha */
	if (Spec[3] != 1)
	{
		lower[jvar] = -1/smax;
		upper[jvar] = -1/smin;
		jvar++;
	}
	/* theta2/beta */
	if (Spec[4] != 1)
	{
		lower[jvar] = -1/smax;
		upper[jvar] = -1/smin;
		jvar++;
	}
	/* rho */
	if (Spec[5] != 1)
	{
		if (restrict == 1)
		{
			lower[jvar] = 1.0;
		}
		else
		{
			lower[jvar] = 0.0;
		}
		upper[jvar] = 18.0; /* Maybe this will be a variable someday */
		jvar++;
	}
	/* Phi/(1-Phi) */
	for(i=6; i<=nparm; i++)
	{
		if (Spec[i] != Yes)
		{
			lower[jvar] = 0.0;
			upper[jvar] = Max_double ;
			jvar++;
		}
	}

	Maxloglik = 0.0;

#ifdef LOGGING_ON
	{
		fprintf(fp_log,"\n*************************Inside Max_lk*************************\n");
		fprintf(fp_log,"                        Before run_dmngb\n");
		fprintf(fp_log,"nvar = %3ld             (The number of varying parameters)\n",nvar);
		fprintf(fp_log,"jfixed = %3d           (The number of fixed parameters)\n",jfixed);
		fprintf(fp_log,"jvar = %3d             (The number of restrictions on the parameters)\n",jvar);
		fprintf(fp_log,"Max_double = %6.3g  (The maximum double)\n",Max_double);
		fprintf(fp_log,"smin = %6.3g          (The minimum litter size)\n", smin);
		fprintf(fp_log,"smax = %6.3g          (The maximum litter size)\n", smax);
		fprintf(fp_log,"Maxloglik = %6.3g\n", Maxloglik);
		fprintf(fp_log,"ITMAX = %3d            (The maximum number of iterations)\n",ITMAX);
		fprintf(fp_log,"gtol = %6.3g          (The convergence critria)\n", gtol);
		fprintf(fp_log,"************Bounds Going Into DMNGB*************\n");

		for (i=0; i<jvar; i++)
			fprintf(fp_log,"lower[%2d] = %6.3g         upper[%2d] = %6.3g\n",i,lower[i],i,upper[i]);
		fprintf(fp_log,"************************************************\n");

		for (i=0; i<nparm; i++)
			fprintf(fp_log,"uiparm[%2d] = %12.6ld\n",i,uiparm[i]);
		fprintf(fp_log,"************************************************\n");

		for (i=0; i<nvar; i++)
			fprintf(fp_log,"start[%2d] = %12.6g\n",i,start[i]);
		fprintf(fp_log,"************************************************\n");

		for (i=0; i<jfixed; i++)
			fprintf(fp_log,"urparm[%2d] = %12.6g\n",i,urparm[i]);
		fprintf(fp_log,"**************************************************************\n");
		fflush(fp_log);
	}
#endif
	/* maximize the log-likelihood.  ErrorFlag gives the return code */
	ErrorFlag = run_dmngb((int) nvar, start, lower, upper, Maxloglik,
		gtol, Parm_Conv, ITMAX, 10,
		NCTR_lk,NCTR_g,uiparm,urparm,ufparm,
		DeBuG,fret);

	/* We actually minimized -log-likelihood.  Want to return log-likelihood */
	*fret = -*fret;

#ifdef LOGGING_ON
	{
		fprintf(fp_log,"\n*************************Inside Max_lk*************************\n");
		fprintf(fp_log,"                        After run_dmngb\n");
		fprintf(fp_log,"*fret = %12.6g   (Maximum log-likelihood value)\n",*fret);
		fprintf(fp_log,"ErrorFlag = %2d         (Status of dmngb; Good <= 6)\n",ErrorFlag);
		fprintf(fp_log,"nvar = %3ld             (The number of varying parameters)\n",nvar);
		fprintf(fp_log,"jfixed = %3d           (The number of fixed parameters)\n",jfixed);
		fprintf(fp_log,"jvar = %3d             (The number of restrictions on the parameters)\n",jvar);
		fprintf(fp_log,"Max_double = %6.3g  (The maximum double)\n",Max_double);
		fprintf(fp_log,"smin = %6.3g          (The minimum litter size)\n", smin);
		fprintf(fp_log,"smax = %6.3g          (The maximum litter size)\n", smax);
		fprintf(fp_log,"Maxloglik = %6.3g\n", Maxloglik);
		fprintf(fp_log,"ITMAX = %3d            (The maximum number of iterations)\n",ITMAX);
		fprintf(fp_log,"************Bounds Coming Out DMNGB*************\n");

		for (i=0; i<jvar; i++)
			fprintf(fp_log,"lower[%2d] = %6.3g         upper[%2d] = %6.3g\n",i,lower[i],i,upper[i]);
		fprintf(fp_log,"************************************************\n");

		for (i=0; i<nparm; i++)
			fprintf(fp_log,"uiparm[%2d] = %12.6ld\n",i,uiparm[i]);
		fprintf(fp_log,"************************************************\n");

		for (i=0; i<nvar; i++)
			fprintf(fp_log,"start[%2d] = %12.6g\n",i,start[i]);
		fprintf(fp_log,"************************************************\n");

		for (i=0; i<jfixed; i++)
			fprintf(fp_log,"urparm[%2d] = %12.6g\n",i,urparm[i]);
		fprintf(fp_log,"***************************************************************\n");
		fflush(fp_log);
	}
#endif

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

	FREE_DVECTOR(lower,0,nvar-1);
	FREE_DVECTOR(upper,0,nvar-1);
	FREE_DVECTOR(start,0,nvar-1);
	FREE_LIVECTOR(uiparm,0,nparm-1);
	FREE_DVECTOR(urparm,0,nparm-nvar-1);
}


/**************************************************************
*NCTR_fit -- Used to "prepar" the data for further computation,
*            i.e. compute the extern variables, give the initial
*            parameters, etc. THEN fit the NCTR model.
*            (In fact, these jobs could be done in main().)
*
***************************************************************/
void NCTR_fit(int n, int ngrp, double p[], double gtol,
			  int *iter, double *fret)
{
	int    *SpBak, i, j, junk, count;
	double sum1, nij, ymin, W, xlk, x;
	double *tmYp, *tmYn, *tmXi,*pBak, tmy, tmvcv;
	double *newX, *newY, *newLs, min, max;
	double meanX, meanY, sumXY, sumX2;

	tmYn  = DVECTOR(1, ngrp);
	tmYp  = DVECTOR(1, ngrp);
	tmXi  = DVECTOR(1, ngrp);
	pBak  = DVECTOR(1, n);
	SpBak = IVECTOR(1, n);
	newX  = DVECTOR(1, Nobs);
	newY  = DVECTOR(1, Nobs);
	newLs = DVECTOR(1, Nobs);

	replace = No;
	tmy = 0.0;
	tmvcv = 0.0;

	/* -----------------------------------------------------------------------
	**   Compute statistics for litter-specific covariate:
	**   smean1, smean, smin, smax,   sdif
	** ----------------------------------------------------------------------*/

	sum1 = 0.0;
	nij = 0.0;
	i=1;
	while (Xg[i]==1)
	{
		sum1 += Ls[i];
		nij +=  1.0;
		i++;
	}

	smean1 = sum1 / nij;  /*  the average litter size in group 1. */

	sum1 = 0.0;
	nij = 0.0;
	smax = Ls[1];
	smin = Ls[1];
	xmax=0.0;
	for (i=1;i<=Nobs;i++)
	{
		x = Xi[i];
		sum1 += Ls[i];
		nij  +=  1.0;
		if (Ls[i] > smax)
			smax = Ls[i];
		if (Ls[i] < smin)
			smin = Ls[i];
		if (x > xmax)
			xmax = x;
	}
	smean = sum1/nij;        /*  overall average litter size. */

	/*  for(i=1; i<=Nobs; i++) Xi[i] *=1000; */
	if (smax<=1.0 || sum1 <=0.0)
		ERRORPRT("all litter sizes are zero.");

	smax = smax-smean;
	smin = smin-smean;       /*  note: smin is negative. */

	if(fixedSize==1)
		sijfixed = smean1 - smean;
	else
		sijfixed=0;

	/* -----------------------------------------------------------------------*/
	/*                            Set up parameters                           */

	if(initial==Yes) /* Inserting user-specified initial values */
	{

		for (i=1;i<=n;i++)
			SpBak[i]= 1 - IniSp[i];

		/*  Have to do this because the fun OUTPUT_Init. */
		OUTPUT_TEXT("\n\n                 User Inputs Initial Parameter Values  ");
		OUTPUT_Init(n,SpBak, IniP, Parm_name);
		/*  obtain user input initial values for unspecified parms. */
		for (i=1; i<=n; i++)
		{
			if(IniSp[i]==1)     /*  have been initialized. */
		 {
			 if(Spec[i]==1 )   /*  check if it is for fixed parm. */
				 Warning("The initial value for the fixed parameter is ignored.");
			 /*  p[i]=p[i]. no change. */
			 else   p[i]=IniP[i];
		 }
			else
		 {
			 /*  check if all the unspecified parms were initialized. */
			 if (Spec[i]==0)
				 ERRORPRT("When the initial option is chosen, one has to initial ALL unspecified parameters.");
		 }
		}

		if (p[1] < 0 || p[2]<0 || p[5]<gamm0)
			ERRORPRT("The initial values have to be: Alpha >= 0,  Beta >= 0 and Rho >= 0 (or 1 when there is restriction on Rho). ");
		if (p[1]+p[3]*smin<0 || p[1]+p[3]*smax <0 || p[2]+p[4]*smin <0 || p[2]+p[4]*smax <0)
			ERRORPRT("The initial values have to be: Alpha+Theta1*Rij >= 0 and Beta+Theta2*Rij >= 0.");
		for (j=6; j<=n; j++)
			if (p[j]<0 || p[j]>=1)
				ERRORPRT("The initial values have to be: 0 <= Phi[j] < 1, for the correlation parameters.");
	}
	else  /* Here, we compute default starting values for all unfixed parameters */
	{
		/*  -------- Save the old Spec[] and p[] vectors --------      */
		for (i = 1; i <= n; i++)
		{
			SpBak[i] = Spec[i];
			pBak[i] = p[i];
		}

		/*      ----- Get the starting values in stages --------- */
		/*  First, get starting values for alpha, beta, and rho */
		i = 1;
		for (j = 1; j <= ngrp; j++)   /*  converse nest_data to dicho_data. */
		{
			tmYn[j]=0;
			tmYp[j]=0;
			while (i <= Nobs && Xg[i]==j)
			{
				tmYn[j] += Yn[i];
				tmYp[j] += Yp[i];
				tmXi[j]  = Xi[i];
				i++;
			}
		}

		/*compute initial estimates */
		ymin = 1.0;
		if(Spec[5]==0) p[5]=1.001;
		for (i = 1; i <= ngrp; i++)
		{
			W = tmYp[i]/(tmYp[i]+tmYn[i]);
			PROBABILITY_INRANGE(&W);
			if (ymin > W)
				ymin = W;
			W = -1.0 * log(1.0-W);
			tmvcv += pow(tmXi[i],p[5]);
			tmy += W;
		}
		p[1] = ymin + 0.001;  /* in case ymin=0 */
		p[2] = tmy/tmvcv;

#ifdef LOGGING_ON
		{
			fprintf(fp_log,"p[1]=%f\np[2]=%f\n\n",p[1],p[2]);
			fprintf(fp_log,"SpBak[1]=%d\nSpBak[5]=%d\n\n",SpBak[1],SpBak[5]);
		}
#endif

		if (restrict == Yes)
			min = 1;
		else
			min = 0;
		max = 10;
		for (i=1; i<=Nobs; i++)
		{
			newY[i] = Yp[i]/(Yp[i]+Yn[i]);
			if (newY[i] >= 1)
				newY[i] = 400;
			else
				newY[i] = -log(1.0 - newY[i]);
		}
#ifdef LOGGING_ON
		fprintf(fp_log,"p[1]=%f\np[2]=%f\n\n",p[1],p[2]);
#endif
		gtol = Rel_Conv;

		/*Due to a bad returned p[5] value from fmin subroutine, this has been temporarily
		*disabled by Micheal Ferree 1/7/04
		*/
		/*initialparms(Nobs, min, max, Xi, newY, p, gtol);*/

		for (i=1;i<=2;i++)
			if (SpBak[i]==1)
				p[i]=pBak[i];
#ifdef LOGGING_ON
		{
			fprintf(fp_log,"SpBak[1]=%d\nSpBak[2]=%d\n\n",SpBak[1],SpBak[2]);
			fprintf(fp_log,"p[1]=%f\np[2]=%f\n\n",p[1],p[2]);
		}
#endif
		if (SpBak[5]==1) p[5]=pBak[5];

		for (i = 6; i <= n; i++)  /*  setup NCTR model. */
		{
			Spec[i] = 1;
			p[i] = 0.0;
		}
		Spec[3] = Spec[4] = 1;
		p[3] = p[4] = 0;

#ifdef LOGGING_ON
		{
			fprintf(fp_log,"\n!!!!First Call To Max_lk.  This is to find alpha, beta, and rho.!!!!\n");
			MxLkCnt = 1;
			fprintf(fp_log,"p[2]=%f\nmaxdose=%f\np[5]=%f\n\n",p[2],maxdose,p[5]);
		}
#endif

		/* Scale starting values by max dose level */
		p[2] = p[2]*pow(maxdose, p[5]);


#ifdef LOGGING_ON
		{
			fprintf(fp_log,"\n!!!!First Call To Max_lk.  This is to find alpha, beta, and rho.!!!!\n");
			fprintf(fp_log,"Scaled starting parameters\n");
			for (i=1; i<=nparm; i++)
		 {
			 fprintf(fp_log,"\np[%d]=%f",i,p[i]);
		 }
		}
#endif

		gtol = Rel_Conv;
		MAX_lk(n, p, gtol, &junk, &xlk);

		/* Uncale parameter estimates by max dose level */
		p[2] = p[2]/pow(maxdose, p[5]);

		/* Now p[1], p[2], and p[5] contain starting estimates for */
		/* alpha, beta, and rho  */

		/*** Second, get initial values for Phi's.      ***********/
		count=0;
		for (i=6; i<=n; i++) count += SpBak[i];
		if (count < ngrp)
		{
			for(i=6; i<=n; i++)         /*  Spec[3], [4] remain as 1.  */
			{
				Spec[i] = SpBak[i];      /*  set back to original ones. */
				p[i] = pBak[i];
				if (SpBak[i]==0)  p[i]=0.01;    /*  set init. if unknown. */
			}
			for(j=6;j<=n;j++)
				p[j]=p[j]/(1-p[j]);     /*  Phi --> Psi. */
			for(j=3;j<=4;j++)
			{
				if(p[j-2]>0)
					p[j]= p[j]/p[j-2];
				else p[j]=0;
			}

#ifdef LOGGING_ON
			{
				fprintf(fp_log,"\n!!!!Second Call To Max_lk.  This is to find the Phi's.!!!!\n");
				MxLkCnt = 2;
			}
#endif
			/* Scale starting values by max dose level */
			p[2] = p[2]*pow(maxdose, p[5]);

			gtol = Rel_Conv;
			MAX_lk(n, p, gtol, &junk, &xlk);

			/* Uncale parameter estimates by max dose level */
			p[2] = p[2]/pow(maxdose, p[5]);

			/* Transform parameters to "external" form */
			for(j=6;j<=n;j++)
				p[j]=p[j]/(1+p[j]);     /*  Psi --> Phi. */
			for(j=3;j<=4;j++)
				p[j]= p[j]*p[j-2];      /*  (theta1/alpha) -> theta1; (theta2/beta) -> theta2 */
		}


		/* Find strating values for theta1 and theta2 but also make sure they stay within
		the constraint. */
		meanX = meanY = sumXY = sumX2 = 0.0;
		for (i = 1; i <= Nobs; i++)
		{
			newLs[i] = Ls[i] - smean;
			newX[i] = pow(Xi[i], p[5]);
			if (newLs[i] == 0)
				newLs[i] = 1e-8;
			if (Yn[i] > 0)
				newY[i] = (-log(1-Yp[i]/(Yp[i]+Yn[i]))-p[1]-p[2]*newX[i])/newLs[i];
			else
				newY[i] = (-log(1e-8)-p[1]-p[2]*newX[i])/newLs[i];
			meanX += newX[i];
			meanY += newY[i];
			sumXY += newX[i]*newY[i];
			sumX2 += pow(newX[i], 2);
		}
		meanX = meanX/Nobs;
		meanY = meanY/Nobs;
		p[4] = (sumXY - Nobs*meanX*meanY)/(sumX2 - Nobs*pow(meanX,2));
		p[3] = meanY - p[4]*meanX;
		if (p[3] < (-p[1]/smax))
			p[3] = (-p[1]/smax) + 1e-8;
		if (p[4] < (-p[2]/smax))
			p[4] = (-p[2]/smax) + 1e-8;

		/** Finally, get initial for Theta's.       **************/
		for (i=3; i<=4; i++)
		{
			Spec[i] = SpBak[i];
			/*  	  p[i] = pBak[i]; */
			if (SpBak[i]==1)
				p[i]=pBak[i];
		}

		OUTPUT_TEXT("\n\n                  Default Initial Parameter Values  ");
		OUTPUT_Init(n, Spec, p, Parm_name);
		fflush(fp_out);
	}

	/* ------------------------------------------------------------------ */
	/*                                  Fit the model                     */

	/* ------- Transform the parameters to the "internal" form -------    */
	for(j=6;j<=n;j++)
		p[j]=p[j]/(1-p[j]);      /*  Phi --> Psi. */
	for(j=3;j<=4;j++)
	{
		if( p[j-2]>0) p[j]= p[j]/p[j-2];
		else p[j]=0;
	}


#ifdef LOGGING_ON
	{
		fprintf(fp_log,"\n!!!!Third Call To Max_lk.  This is to optimize over all parms.!!!!\n");
		fflush(fp_log);
		MxLkCnt = 3;
	}
#endif

	/* Scale starting values by max dose level */
	p[2] = p[2]*pow(maxdose, p[5]);

	gtol = Rel_Conv;
	MAX_lk(n, p, gtol, &junk, &xlk);

	/* Uncale parameter estimates by max dose level */
	p[2] = p[2]/pow(maxdose, p[5]);

	/* junk is 'iter' in MAX_lk, which does not get set there.  Just ignore
	this code
	if (junk>=ITMAX)
	Warning("Warning: Maximum iteration may be not large enough. Iterations reach the maxmimum.");
	*/
	/* Print out a warning if ErrorFlag is non-zero */
	do_dmngb_warning(&ErrorFlag);

	*fret = xlk;


	/* ---------- Transform the parameters to the "external" form ------- */
	for(j=6;j<=n;j++)
		p[j]=p[j]/(1+p[j]);      /*  Psi --> Phi. */
	for(j=3;j<=4;j++)
		p[j]= p[j]*p[j-2];       /*  (theta1/alpha) -> theta1; (theta2/beta) -> theta2 */

	FREE_DVECTOR(tmYn, 1, ngrp);
	FREE_DVECTOR(tmYp, 1, ngrp);
	FREE_DVECTOR(tmXi, 1, ngrp);
	FREE_DVECTOR(pBak, 1, n);
	FREE_IVECTOR(SpBak, 1, n);
	FREE_DVECTOR(newX, 1, Nobs);
	FREE_DVECTOR(newY, 1, Nobs);
	FREE_DVECTOR(newLs, 1, Nobs);

}


/************************************************************
* NCTR_BMD -- Used to calculate the BMD and BMDL for NCTR model.
*
*************************************************************/
void NCTR_BMD (int nparm, double p[], double gtol, int *iter, double xlk,
			   double Rlevel[], double Bmdl[], double *BMD)
{
	double   tol;
	double   xa,xb,fa,fb;
	double   D, Bmdjunk,lamb;
	double   *pBak, *porig;
	int      j, k;

	pBak=DVECTOR(1, nparm);
        porig=DVECTOR(1, nparm);

        for(j=1; j<=nparm; j++) porig[j]=p[j]; /* save orig p val */

	/**** compute X^2 value  ************************/
	/* IF ML is the value of the maximized log-likelihood, then ML - LR is the value
	log-likelihood at the BMDL or BMDU */
	if (bmdparm.level<0.5)
		LR = 0.5*QCHISQ(1.0 - 2.0 * bmdparm.level, 1);
	else
		LR = 0.5*QCHISQ(2.0 * bmdparm.level - 1.0, 1);



	for(j=6; j<=nparm; j++) p[j]=p[j]/(1-p[j]);      /*  Phi --> Psi. */
	for(j=3;j<=4;j++)
	{
		if( p[j-2]>0) p[j]= p[j]/p[j-2];
		else p[j]=0;
	}

	for(j=1; j<=nparm; j++) pBak[j]= p[j];         /*  save the p[]. */

	Rlevel[1] = BMR = bmdparm.effect;

	lamb = exp(pBak[1]*(1+pBak[3]*sijfixed));
	/*  Note: have to used pBak[] instead of p[]. */
	if (bmdparm.risk == 1)
		ck = -1.0*log((1.0-BMR*lamb));
	else
	{
		ck = -log(1-BMR);
		/*  Spec[1]=Spec[3]=1; */
	}

	/**** solve the BMD ********************************/
	if (p[5]<= (log(ck/(pBak[2]*(1+pBak[4]*sijfixed)))/log(Max_double)) )
		ERRORPRT("The power parameter is zero (or too close to zero). BMD is undefined.");
	xb = pow(ck/(pBak[2]*(1+pBak[4]*sijfixed)), 1/pBak[5]);
	*BMD = xb;
//	if(fixedSize==1)
//	{
//		fprintf(fp_out, "\n\nTo calculate the BMD and BMDL, the litter\n");
//		fprintf(fp_out, "specific covariate is fixed at the mean litter\n");
//		fprintf(fp_out, "specific covariate of control group: %f", sijfixed+smean);
//	}
//	else
//	{
//		fprintf(fp_out, "\n\nTo calculate the BMD and BMDL, the litter\n");
//		fprintf(fp_out, "specific covariate is fixed at the overall mean\n");
//		fprintf(fp_out, "of the litter specific covariates: %f", sijfixed+smean);
//	}
//
//	OUTPUT_BENCHMD(1,*BMD);

	Spec[2]=1;
	replace=Yes;        /*  replace is extern var., has to changed now. */




	/********* search for BMDL **************************/
	xa = xb*0.5;
	tol = sqrt(DBL_EPSILON);
	BMD_lk = xlk;  /*  get the lk at BMD. */
	fb = -LR;
	fa = BMDL_func(nparm, p, xa, tol);

	while (fa<0.0 && xa > DBL_EPSILON)
	{
		xa *= 0.5; /*  prevent that xa=0.1*BMD is not small enough. */
		fa = BMDL_func(nparm, p, xa, tol);
	}
	if (fa<0.0)
	{
#ifndef RBMDS
		fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  No); /*  computation failed */
#endif
		ERRORPRT("Benchmark dose computation failed.  Lower limit includes zero.");
	}

#ifndef RBMDS
	fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  Yes); /*  computation will succeed. */
#endif
	Bmdl[1] = zeroin(xa, xb, 1.0e-10, BMDL_func, nparm, p, 1.0e-14);
#ifdef MISC_OUT
	printf("           BMDL =%15.6g\n\n", (Bmdl[1]));
#endif
//#ifndef RBMDS
//	fprintf(fp_out, "            BMDL =%15.6g\n\n", Bmdl[1]);
//#else
//	fprintf(fp_out, "            BMDL =%30.22g\n\n", Bmdl[1]);
//#endif
	fflush(fp_out);

	if(bmdlCurve==Yes)
	{
		/****** calculate Bmd[] and Bmdl[] ***********/
		for (k=2; k<=5;k++)
		{

			for(j=1; j<=nparm; j++)
				p[j]= pBak[j];          /*  get the "old" p[]. */

			if (k==2) Rlevel[k]=BMR=0.05;
			else Rlevel[k]= BMR = (k-2)*0.1;

			if (bmdparm.risk == 1)
				ck = -1.0*log((1.0-BMR*lamb));
			else
			{
				ck = -log(1-BMR);
				/*  Spec[1]=Spec[3]=1; */
			}

			/**** solve the BMD ********************************/
			D = 0.0;
			xb = pow(ck/(p[2]*(1+p[4]*sijfixed)), 1/p[5]);
			Bmdjunk = xb;

			/********* search for BMDL **************************/
			xa = xb*0.1;
			tol = FMAX((Bmdjunk)*0.001, 0.0000001);
			BMD_lk = xlk;  /*  get the lk at BMD. */
			fb = -LR;
			fa = BMDL_func(nparm, p, xa, tol);

			if (fa<0.0)
			{
				xa *= 0.01; /*  prevent that xa=0.1*BMD is not small enough. */
				fa = BMDL_func(nparm, p, xa, tol);
				if (fa<0.0)
				{
					Bmdl[k]=-1;
					fprintf(fp_out, "\n BMDL curve computation failed for BMR = %f . \n",BMR);
					fprintf(fp_out, "The BMDL curve appearing in the graph may not be accurate.");
				}
				else
					Bmdl[k] = zeroin(xa, xb, 1.0e-10, BMDL_func, nparm, p, 1.0e-14);
			}
			else
				Bmdl[k] = zeroin(xa, xb, 1.0e-10, BMDL_func, nparm, p, 1.0e-14);
		}
	}
	else for (k=2; k<=5;k++) Bmdl[k] = Rlevel[k]= -1;

        for(j=1; j<=nparm; j++) p[j]=porig[j]; /* reset orig p val */

        FREE_DVECTOR(porig, 1, nparm);
	FREE_DVECTOR(pBak, 1, nparm);
}  /*end BENCHMD*/


/*****************************************************************
* BMDL_func -- used to compute the values of functions BMDL_f (the
*              X^2 value) at the point D, given the parm p[] and
*              number of parm.
*
*              This routine is called by zeroin.
*
*****************************************************************/
double BMDL_func(int nparm, double pBak[], double D, double gtol)
{ 	/*  ck , BMD_lk and LR are calculated in NCTR_BMD() */

	double fD, xlk, *p;
	int j, junk;
	p=DVECTOR(1,nparm);

	for (j=1;j<=nparm;j++) {
		p[j]=pBak[j];  /*  get the "old" p[]. */
#ifdef LOGGING_ON
		fprintf(fp_log,"\n***In BMDL_func(): p[%d]=%f\n",j,p[j]);
#endif
	}

	tD = D;   /*  tD is grobale var. have to change before call MAX_lk(). */
	p[2]= ck/(1+p[4]*sijfixed)/pow(tD, p[5]);

	gtol = Rel_Conv;
	MAX_lk( nparm, p, gtol,  &junk,  &xlk);
	fD = BMD_lk - xlk - LR;

#ifdef LOGGING_ON
	{
		fprintf(fp_log,"\n\n***fD=%f\n", fD);
		fprintf(fp_log,"\n\n***BMD_lk=%f\n", BMD_lk);
		fprintf(fp_log,"\n\n***xlk=%f\n", xlk);
		fprintf(fp_log,"\n\n***LR=%f\n", LR);
	}
#endif
	FREE_DVECTOR(p, 1, nparm);
	return fD;
}


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
	double bkg;
	for (i = 1; i <= nobs; i++)
	{
		bkg = 1-exp(-(Parms[1] + Parms[3]*(Lsc[i]-smean)));
		if (doses[i] <= 0.0)
			P[i] = bkg;
		else
			P[i] = 1-exp(-(Parms[1]+Parms[3]*(Lsc[i]-smean))-
			(Parms[2]+Parms[4]*(Lsc[i]-smean))*pow(doses[i],Parms[5]));
	}
}


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
