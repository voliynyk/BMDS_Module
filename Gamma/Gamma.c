/****************************************************************
*
* IMPORTANT NOTE:  The following variable is the version number for
*                  the current model.  THIS MUST BE CHANGED as
*				   important changes are made to the models.
*
*****************************************************************/
char Version_no[] = "Gamma Model. (Version: 2.16;  Date: 2/28/2013)";
/****************************************************************
char Version_no[] = "Gamma Model. (Version: 2.15;  Date: 10/28/2009)";
char Version_no[] = "Gamma Model. (Version: 2.14;  Date: 11/23/2008)";
char Version_no[] = "Gamma Model. (Version: 2.12;  Date: 03/10/2008)";
char Version_no[] = "Gamma Model. (Version: 2.11;  Date: 10/31/2007)";
char Version_no[] = "Gamma Model. (Version: 2.8;  Date: 02/20/2007)";
*char Version_no[] = "Gamma Model. (Version: 2.6;  Date: 12/08/2006)";
*char Version_no[] = " Revision: 2.2 Date: 2001/07/11 11:55:57";
*****************************************************************/

/****************************************************************
**
* Gamma.c - a ANSI C program for Gamma model fitting with/without
*           a natural background rate in Benchmark Dose.
*
*
* Date: Oct 7, 1997
*
********************************************************************
* Modification Log:
*
* Version Number: 2.1
* Modified By: R. W. Setzer
* Modified Date: 2/25/2000
* Reason:
*
* Version Number: 2.4
* Modified By: Micheal Ferree
* Modified Date: 8/13/2005
* Reason: Took out print out of cancer slope factor
* 
* Version Number: 2.5
* Modified By R. W. Setzer
* Modification Date: 9/29/2005
* Reason: Free all allocated memory
*
* Version Number: 2.6
* Modified By R. W. Setzer
* Modification Date: 12/08/2006
* Reason: Only free vcv_adj if it was allocated
*
* Version Number: 2.7
* Modified By: Geoffrey
* Date: 1/12/2007
* Reason: Incremented version number.
*		 Added last parameter "1" (print SE) in OP_ParmsE().
*
* Version Number: 2.8
* Modified By: Woodrow Setzer
* Date: 2/20/2007
* Reason: Incremented version number to reflect changed compilation options

* Version Number: 2.9
* Modified By: G. Nonato
* Modified Date: 9/23/2007
* Reason:	1. Changed version number to reflect changed.
*		2. Prefix condition with "Xi[1] == 0.0 &&" in Line 441 and beautify code 
*                 {GLN - 09/23/07, PR0823-05} 
*		   It fixes the following:
*		      Different incorrect responses result when running models on data for
*		      which there is no zero dose group and background is specified to be 0.
*		      For the  example data set where doses are 17, 20 and 24, Ns are all 10
*		      and responses are 1, 2 and 6
*
* Version Number: 2.10
* Modified By: Woodrow Setzer
* Date: 10/16/2007
* Reason:  Initial values were incorrect when the lowest dose level was nonzero.
*          In that case, we now set the background to be half the response level
*          at the lowest dose, and use that value when looking for estimates
*          of the other parameters.
*
*
* Version Number: 2.11
* Modified By: G. Nonato
* Modification Date: 10/31/2007
* Reason: (Per Technical Direction #15, 10/31/07 email)
*       Fix the background displaying differently in out file (L1713).
*       Fix the unhandled Win32 Exception (freeing mem variable that was non existing) in L938.
*
* Version Number: Version: 2.12
* Modified By: G. Nonato
* Modification Date: 03/10/2008
* Reason: (Per BMDS 2.0: Problem Report 157 & 147)
*       Fix the Observation # < parameter # for Weibull model problem.
*       Added code to free-up allocated memories before exiting thru ERRORPRT()
*
* Version Number: 2.13
* Modified By: G. Nonato
* Modification Date: 05/16/2008
* Reason: (Per BMDS 2.0: Problem Report 165)
*       Goodness of Fit - Observed column, print values as real numbers.
*
* Version Number: 2.14
* Modified By: R. W. Setzer
* Modification Date: 11/23/2008
* Reason:
*       Changed initial value code to estimate phat with (a + 1)/(a + b + 2)
*       This fixes a rare condition in which the program never returns.
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
* Modified By: Louis Olszyk
* Modification Date: 02/28/2013
* Reason: PR 444 - Fix wording in plot titles
******************************************************************/

#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>

/* Forward-declared functions */

void Gamma_fit(int nparm, double p[], double gtol, int *iter,
			   double *fret);
void Gamma_BMD (int nparm, double p[], double gtol, int *iter,
				double xlk, double Rlevel[], double Bmdl[],double *BMD);
void Gamma_vcv(int nparm, int Spec[], double p[], double **vcv);
void Which_Bounded (int [], double [], int []);
int Model_DF (int []);
double BMDL_func (int, double [], double, double);

#include <benchmark.h>
#include <ERRORPRT.h>
#include <allo_memo.h>
#include <matrix_agb.h>
#include <specialfun.h>
#include <computation.h>
#include <in_outfun.h>



#define EPS 3.0e-8
#define TOLX (10*EPS)
#define STPMX1 10.0
#define ALF 0.000001
#define FREEALL FREE_DVECTOR (xi,1,n); FREE_DVECTOR (pnew, 1,n);	\
	FREE_DVECTOR(hdg,1,n); FREE_DVECTOR(g,1,n);				\
	FREE_DVECTOR(dg,1,n);							\
	FREE_DMATRIX(hessin,1,n,1,n);

#define VERYCLOSE(x,y) (fabs((x) - (y)) < Min_increment)

#define float double
#define GOLD 1.618034
#define GLIMIT 100
#define TINY 1.0e-20
#define SWAP(a,b, junk) (junk)=(a); (a)=(b); (b)=(junk);
#define SHFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d);
#define ZEPS 1.0e-8
#define MOV3(a,b,c, d,e,f)(a)=(d);(b)=(e);(c)=(f);
#define TOL 2.0e-6

//#define DO_LOG		//Uncomment to generate log file code
/*************************************************************
*   giDo_Log = true -> log file is made
*   giDo_Log = false -> no log
*************************************************************/
int    giDo_Log = false;     /*  Uncomment to activate switch for log file   */
char   gacLogFile[FLENGTH], *gcDot2;

/* define input and output files's name */
char  fin[FLENGTH];		/* input temp file */
char  fout[FLENGTH];	/* output temp file */
char  fout2[FLENGTH];
char  plotfilename[FLENGTH];	/* file to pass to GnuPlot */
char  *Parm_name[]={"Background", "Slope", "Power"};
/* define a typedef enum for parameter names */
/* Base1 is because all the parameter vectors, spec, etc, are based at 1 */
typedef enum {Base1, Background, Slope, Power} Model_Parms;
char  *anatxt[]={"Full model", "Fitted model", "Reduced model"};


/* variables will not be changed execpt Spec */
int    *Spec;	/* vector used to identify user input parm.
				will be 1 if user specified that parameter,
				0 otherwise */
int    *IniSp;	/* user initial specified parameter vector */
double *Yp;	/* positive dependent variable data array */
double *Yn;	/* negative dependent variable data array */
double *Xi;	/* independent variable data array */
double *Ypp;	    /* predicted positive dependent values */
double *Ep;	    /* estimated probability  */
double *Rlevel;	    /* risk level (BMR level) */
double *Bmdl;	    /* benchmark dose confidence limit */
double *IniP;	    /* user intialized parameter vector */


int    Nobs;	    /* number of observations */
int    nparm;	    /* number of parameters */
int    restrict;    /* flag for restricting slope>=1 */
int    initial;	    /* flag for user initial parameter estimates */
int    appendix;    /* flag: 0 for append output, 1 for overwrite */
int    smooth;	    /* flag: 0 for unique bmdl curve, 1 for C-spline */
int    bmdlCurve;   /* flag for bmdl curve option */
double xmax, xmin, scale;

double Rel_Conv, Parm_Conv, Maxloglik;
double SlopeUpperBound = 18.0;

double BMDL_Error_Size;
int BMDL_Error;
int DeBuG = 0;


/* changing variable */
int		replace, brat;
double  tD,  BMD_lk, LR, ck, upb=18, BMR;

int ErrorFlag; /* Error States from DMNGB */

/*************************************************************** */
/* fixedParm -- Returns 1 if Spec[(int) parm] is 1, 0 otherwise */
/*              That is returns true if the requested parameter is */
/*              fixed */
/*************************************************************** */

int fixedParm ( Model_Parms i)
{
	return Spec[(int) i];
}

/**************************************************************** */
/* initParm -- Returns 1 if IniSp[(int) parm] is 1, 0 otherwise */
/*             That is, returns true if the requested parameter */
/*             has a user-specified initial value */
/**************************************************************** */

int initParm ( Model_Parms i)
{
	return IniSp[(int) i];
}

/*************************************************************** */
/* allFixed -- returns 1 if all parameters have been specified */
/*             (thus, no fitting is required) */
/*************************************************************** */

int allFixed (void)
{
	int i;
	int res;

	res = 1;
	for ( i = 1; i <= nparm; i++) res = res && Spec[i];
	return res;
}

/****************************************************
**  CLOSE_FILES--used to close input and output files.
*****************************************************/
void CLOSE_FILES (void)
{
	if (fclose (fp_in) != 0 || fclose (fp_out) != 0 
#ifndef RBMDS
		|| fclose (fp_out2) != 0  
#endif
		)
	{
		ERRORPRT ("Error in closing opened files.");

	} /* end if */
}

/* GLN - 03/07/2008
*  Free-up allocated memory before exit
*  upon encountering fatal error.
*/
void FreeUp_mem(double *Parms, VarList *varsum, AnaList *anasum, double  **vcv)
{
	FREE_DVECTOR (Parms,1,nparm);
	FREE_DVECTOR (Ypp,1,Nobs);
	FREE_DVECTOR (Ep,1,Nobs);
	FREE_DVECTOR (Xi,1,Nobs);
	FREE_DVECTOR (Yp,1,Nobs);
	FREE_DVECTOR (Yn,1,Nobs);
	FREE_DVECTOR (IniP,1,nparm);
	FREE_IVECTOR (Spec,1,nparm);
	FREE_IVECTOR (IniSp,1,nparm);
	FREE_VLVECTOR(varsum,1,3);
	FREE_ALVECTOR(anasum,1,3);
	FREE_DVECTOR (Rlevel,1,5);
	FREE_DVECTOR (Bmdl, 1, 5);
	FREE_DMATRIX(vcv,1,nparm,1,nparm);
	if (fp_log != (FILE *) NULL)
		fclose(fp_log);

	return;
}

/****************************************************************
** main--main function used to call Gamma mode fitting program.
Includes: biosubcc.c--common subfunction C program.
*****************************************************************/
int main (int argc, char *argv[])
{

	int	  iter,i, junk;	     /* iteration variable */
	int     bmdose;	     /* flag for computing benchmark dose */
	int     Nmiss;	     /* number of records with missing values */
	int     nparm_known;	     /* number of specified parameters */
	double  lkf,lkr,xlk,W;     /* log likelihoods */
	double  back, BMD, back1;

	double  *Parms;	     /* parameter array */
	VarList *varsum;	     /* info for variables--p. dep.,n. dep., indep. */
	AnaList *anasum;	     /* information for ANONA analysis */
	double  **vcv;	     /* variance and covariance matrix */

	char    model_name[MNLENGTH], user_note[UNLENGTH], data_file[FLENGTH];
	char    dose_name[CNLENGTH], posi_name[CNLENGTH], nega_name[CNLENGTH], junkname[FLENGTH];
	int	  adj_vcv_rows, *bounded;
	double  **vcv_adj;
	char long_path_name[FLENGTH];
	time_t ltime;

	time( &ltime );

	/* These are defined in float.h */
	Min_increment = DBL_EPSILON;
	Max_double = DBL_MAX;

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

	/* Check if input file is open, if not, print error message and exit */
	if (fp_in==NULL)
	{
		fprintf (stderr, "Error in opening input  file.\n");
		fprintf (stderr, "...now exiting to system...\n");

		exit (1);

	} /* end if */

	#ifdef DO_LOG
	strcpy(gacLogFile,argv[1]);
	gcDot2 = strchr(gacLogFile, (int) '.');
	(*gcDot2) = (char) 0;
	strcat(gacLogFile,".log");
	if (giDo_Log)
	{
		fp_log = fopen(gacLogFile, "w");

		if (fp_log == (FILE *) NULL)
			ERRORPRT("Unable to open log for Weibull.");

		fprintf(fp_log,"\n\nargv[1] (before) = %s", argv[1]);
	}
	#endif

	/* begin reading input file from batch file (.(d) ext.) */
//	fscanf(fp_in, "%s",model_name );
//
//	fscanf(fp_in, "%[ ^\n]",user_note );
//	//fscanf(fp_in, "%s",user_note );
//
//	fscanf(fp_in, "%[^\n]", data_file);
//	//fscanf(fp_in, "%s", data_file);
//
//	fscanf(fp_in, "%[^\n]", junkname);
//	//fscanf(fp_in, "%s", junkname);
//
//	fscanf(fp_in, "%d",&Nobs);

	fscanf(fp_in, "%s",model_name );
	fscanf(fp_in, "%[ ^\n]",user_note );
	fscanf(fp_in, "%[^\n]", user_note);
	fscanf(fp_in, "%s", junkname);
	fscanf(fp_in, "%s", junkname);
	fscanf(fp_in, "%d",&Nobs);

	/* assign number of parameters */
	nparm = 3;
	/* allocate memory for arrays */
	Parms = DVECTOR(1, nparm);
	IniP = DVECTOR(1, nparm);
	IniSp= IVECTOR(1, nparm);
	Spec = IVECTOR(1, nparm);
	Ypp = DVECTOR(1, Nobs);
	Ep = DVECTOR(1, Nobs);
	Xi = DVECTOR(1, Nobs);
	Yp = DVECTOR(1, Nobs);
	Yn = DVECTOR(1, Nobs);
	varsum = VLVECTOR(1, 3);
	anasum = ALVECTOR(1, 3);
	Rlevel = DVECTOR(1, 5);
	Bmdl = DVECTOR(1, 5);
	vcv = DMATRIX (1,nparm,1,nparm);


	/* read more values from input file */
	fscanf(fp_in,"%d%lf%lf%d%d%d%d%d", &ITMAX, &Rel_Conv, &Parm_Conv,
		&bmdlCurve, &restrict, &bmdose, &appendix, &smooth);
	fscanf(fp_in,"%lf%d%lf",&bmdparm.effect,&bmdparm.risk,&bmdparm.level);

	/* bmdparm.effect is the benchmark risk level */
	/* bmdparm.risk is 1 is extra risk, 0 if added risk */
	/* bmdparm.level is the confidence level */

	junk = 0;    /* used to see if an extension was added to output file
				 name */

	/* get filenames */
	Get_Names(argv[1], fout, fout2, plotfilename);

	/* open output files */
	if(appendix==Yes)
	{
		fp_out=fopen(fout,"a");		/* append output */
	}
	else
	{
		fp_out=fopen(fout,"w");		/* overwrite output */

	} /* end if */
#ifndef RBMDS
	fp_out2=fopen(fout2,"w");	/* Always overwrite plotting file */
#endif
	/* Check to make sure files are open, if not, print error message and
	exit */
	if (fp_out==NULL 
#ifndef RBMDS
		||fp_out2==NULL 
#endif
		)
	{
		/*        printf("Error in opening  output files.\n"); */
		/*        printf ("...now exiting to system...\n"); */
		FreeUp_mem(Parms, varsum, anasum, vcv);
		fprintf(fp_out,"Error in opening output files.\n");
		fprintf (fp_out,"...Exited to system!\n");
		exit (1);

	} /* end if */


	/* Print model and file information on output page */
	Output_Header(Version_no, argv[1], plotfilename, ctime(&ltime), user_note);

	if (bmdose < 0 || bmdose > 1)
	{
		FreeUp_mem(Parms, varsum, anasum, vcv);
		ERRORPRT("Error in choosing benchmark dose computation.");

	} /* end if */

	/* obtain user input parameters */
	READ_PARAMETERS(nparm,Parms);
	FILL_SPECVECTOR(nparm,Parms,Spec);

	nparm_known = COUNT_SPECVECTOR(nparm,Spec);
	brat=Yes;

	if (Spec[(int) Background ]==1 && Parms[ (int) Background ]<EPS)
	{
		brat=No;

	} /* end if */

	/* obtain user input initial parameters values */
	fscanf(fp_in,"%d", &initial);
	READ_PARAMETERS(nparm,IniP);
	FILL_SPECVECTOR(nparm,IniP,IniSp);

	for(i = 1; i <= nparm; i++)
	{
		if(Spec[i] == 1)
		{
			IniP[i] = 1;

		} /* end if */

	} /* end for */

	/* obtain observation data into Yp, Yn, Xi, vectors
	Yp = number of respondants
	Yn = number on nonrespondants
	Xi = dose level */

	fscanf(fp_in,"%s%s%s", dose_name, posi_name, nega_name);
	Nmiss = READ_OBSDATA3V(Nobs,3,2,3,1, Yp,Yn, Xi);


	Nobs -= Nmiss;             /* extern variable Nobs has been changed */
	#ifdef DO_LOG
	if (giDo_Log)
	{
		fprintf(fp_log,"\nNobs = %d", Nobs);
		fprintf(fp_log,"\nnparm = %d", nparm);
		fprintf(fp_log,"\nnparm_known = %d", nparm_known);
		fprintf(fp_log,"\nNobs < (nparm-nparm_known)");
	}
	#endif

	//if (Nobs < nparm), commented and changed to the code below, GLN 3/7/08
	if (Nobs < (nparm-nparm_known))
	{
		FreeUp_mem(Parms, varsum, anasum, vcv);
		ERRORPRT("Observation # < parameter # for Gamma model.");
	}  /* end if */

	/* end of input data */

	/* output title and summary of input data */
	OUTPUT_TEXT("\n   The form of the probability function is: ");
	OUTPUT_TEXT("\n   P[response]= background+(1-background)*CumGamma[slope*dose,power],");
	OUTPUT_TEXT("   where CumGamma(.) is the cummulative Gamma distribution function");
	fprintf(fp_out,"\n\n   Dependent variable = %s", posi_name);
	fprintf(fp_out,"\n   Independent variable = %s", dose_name);

	if (fixedParm(Background) == Yes)    /* if background parameter is specified */
	{
		if (Parms[(int) Background] <= 0.0000001)
		{
			fprintf(fp_out,"\n   Background parameter is set to zero");
		}
		else
		{
			fprintf(fp_out,"\n   Background parameter is set to %g",
				Parms[1]);

		} /* end if (Parms[1] <= 0.0000001) */
		if(Xi[1] == 0.0 && Parms[(int) Background ] == 0.0 && Yp[1] != 0) //Prefix condition with "Xi[1] == 0.0 &&" {GLN - 09/23/07, PR0823-05}
		{
			FreeUp_mem(Parms, varsum, anasum, vcv);
			ERRORPRT("ERROR:  Background parameter specified as 0, but there were responses\n        at the control dose\n");

		} /* end if */


	}	/* end if (fixedParm(Background) == Yes) */


	if (fixedParm(Slope) == Yes) /* if slope parameter is specified */
	{
		fprintf(fp_out,"\n   Slope parameter is set to %g", Parms[2]);

	} /* end if */


	if (fixedParm(Power) == Yes)		/* if power parameter is specified */
	{
		fprintf(fp_out,"\n   Power parameter is set to %g", Parms[3]);
	}
	else
	{
		if (restrict==Yes)	/* if power parm. is restricted to >= 1 */
		{
			fprintf(fp_out,"\n   Power parameter is restricted as power >=1");
		}
		else
		{
			fprintf(fp_out,"\n   Power parameter is not restricted");

		} /* end if (restrict==Yes) */

	}	/* end if (Spec[3]==Yes) */

	/* end output of specified parameters */

	nparm_known = COUNT_SPECVECTOR(nparm, Spec);
	fprintf (fp_out, "\n\n   Total number of observations = %d",Nobs+Nmiss);
	fprintf (fp_out, "\n   Total number of records with missing values = %d",
		Nmiss);

	fprintf(fp_out, "\n   Maximum number of iterations = %d\n", ITMAX);
	fprintf(fp_out, "   Relative Function Convergence has been set to: %g\n",
		Rel_Conv);
	fprintf(fp_out, "   Parameter Convergence has been set to: %g\n\n",
		Parm_Conv);


	/* Print out user input initialization parameters */
	if(initial==Yes)
	{
		OUTPUT_TEXT("\n\n                 User Inputs Initial Parameter Values  ");
		OUTPUT_Init(nparm, Spec, IniP, Parm_name);

		for (i=1; i<=nparm; i++)
		{
			if(IniSp[i]==1)			/* have been initialized */
			{
				if(Spec[i]==1 )			/* check if it is for fixed parm */
				{
					Warning("The initial value for the fixed parameter is ignored.");

				} /* end if */
			}
			else
			{
				/* check if all the unspecified parms were initialized */
				if (Spec[i]==0)
				{
					ERRORPRT("When the initial option is chosen, one has to initial ALL unspecified parameters.");

				} /* end if */

			} /* end if (IniSp[i]==1) */

		} /* end for */

		if (IniP[1] < 0 || IniP[1]>1 || IniP[2]<0 || IniP[3]<restrict)
		{
			ERRORPRT("The initial values have to be: 1 > background >= 0,  slope >= 0 and 18> power >= 0 (or 1 when there is restriction on power >=1). ");

		} /* end if */

	}	/* end if (initial == yes) */

	/* compute init_lkf for full model and init_lkr for reduced model
	*  The full model assumes that the number of responders at each dose
	*  is independently binomially distributed with individual probability
	*  denoted by p[i] with MLE
	*  W = (number of respondents at dose i)/(total subjects at dose i).
	*  The reduced model assumes that the respondents all come from a
	*  single binomial distribution with probability p and MLE
	*  W = (total number of respondents) / (total subjects)
	*  lkf and lkr are the log-likelihood functions for the binomail dist*/

	lkf = 0.0;
	varsum[1].S = 0;
	varsum[2].S = 0;

	for (i=1;i<=Nobs;i++)
	{
		/* likelihood for full model */
		varsum[1].S += Yp[i];		/* Total number of responders */
		varsum[2].S += Yn[i];		/* Total number of nonresponders */
		W = Yp[i] / (Yp[i]+Yn[i]);	/* Proportion of positive responders */

		if (W > 0)
		{
			lkf += Yp[i] * log(W);

		} /* end if */

		if (W < 1)
		{
			lkf += Yn[i] * log(1- W);

		} /* end if */

	}	/* end for */
	Maxloglik = -lkf;
	/* likelihood for reduced model */
	W = varsum[1].S / (varsum[1].S + varsum[2].S);
	lkr = varsum[1].S * log(W) + varsum[2].S * log(1- W);

	/* fitting Gamma model and output parameter estimators */
	Gamma_fit(nparm, Parms, EPS, &iter, &xlk);

	bounded = IVECTOR(1, nparm);
	Which_Bounded (Spec, Parms, bounded);
	if (ErrorFlag != 0) /* >0 means didn't converge; <0 means all parms fixed */
	{
		bmdose = No;
		for (i = 1; i <= nparm; i++) bounded[i] = 1;
	}
	if (ErrorFlag == 0)
	{
		/* compute the approx. covariance matrix */
		INITIALIZE_DMATRIX(vcv, nparm, nparm);
		Gamma_vcv(nparm,Spec,Parms,vcv);

		/* check for and count the number of bounded and specified
		parameters for use in creating the asymptotic correlation matrix */

		/* count the number of bounded and specified parameters */
		adj_vcv_rows = 0;
		for (i=1; i<=nparm; i++)
		{
			if (bounded[i] == 0)
			{
				adj_vcv_rows++;

			} /* end if */

		} /* end for */

		vcv_adj = DMATRIX(1, adj_vcv_rows, 1, adj_vcv_rows);

		/* create and output asymptotic correlation matrix */
		if (adj_vcv_rows > 0)
		{
			Get_and_OUTPUT_DTMSVCV(nparm, Spec, Parm_name, vcv, vcv_adj, bounded);
		}
	}
	/* Output parameter values and standard errors (if not bounded) */
	if (ErrorFlag != -1)
		OP_ParmsE(nparm,Spec,Parms,Parm_name,vcv_adj, bounded, bmdparm.level, 1);
	/* compute and output ANOVA table elements */
	DTMS3ANOVA (nparm,Nobs,Spec,lkf,xlk,lkr,anasum, bounded);

	/* output ANOVA table */
	OUTPUT_DTMS3ANOVA(anatxt,anasum);
	fflush(fp_out);
	if (ErrorFlag == 0)
	{
		/* print a goodness of fit table */
		Quantal_Goodness(nparm, bounded, Parms, Nobs, Xi, Yp, Yn, scale);
		fflush(fp_out);
		/* compute benchmark dose */

		/* setup plotting file (.002 extension) */
	}
#ifndef RBMDS
	fprintf (fp_out2, "\n BMD_flag \t %d \n Nobs \t%d \n nparm \t%d",
		bmdose, Nobs, nparm );
	fprintf (fp_out2, "\n  Con_lev \t%3.3g ", bmdparm.level);
	fprintf (fp_out2, "\n  RiskType \t%d ", bmdparm.risk);
	fprintf (fp_out2, "\n  Effect \t%3.3g ", bmdparm.effect);

	for (i=1;i<=nparm; i++)
	{
		fprintf (fp_out2, "\n %s \t %5.3g", Parm_name[i-1], Parms[i]);

	} /* end for */
	fflush(fp_out2);

	/* Calculation of 95% confidence intervals at each dose level for
	graphical output */

	{
		double *LL, *UL, *estp;

		LL = DVECTOR(1, Nobs);
		UL = DVECTOR(1, Nobs);
		estp = DVECTOR(1, Nobs);
		Quantal_CI (Nobs, Yp, Yn, 0.95, LL, estp, UL);
		fprintf (fp_out2,"\n\n Data");
		for (i = 1; i <= Nobs; i++)
		{
			fprintf (fp_out2,"\n %f %f %f %f", Xi[i], estp[i], LL[i], UL[i]);
		} /* end for, calculation of confidence intervals */
		FREE_DVECTOR(LL, 1, Nobs);
		FREE_DVECTOR(UL, 1, Nobs);
		FREE_DVECTOR(estp, 1, Nobs);
	}

	fprintf (fp_out2,"\n Max_Min_dose \n  %f %f ", xmax, xmin);
	fflush(fp_out2);
#endif
	/* check if BMD is to be calculated */
	if (bmdose==Yes)
	{
		if(fixedParm(Slope))
		{
#ifndef RBMDS
			fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  No);
#endif
			fprintf (fp_out,"\n\n  %s parameter is fixed. The likelihood function can not be reparameterized in BMD.", Parm_name[1]);
			FreeUp_mem(Parms, varsum, anasum, vcv);
			exit(0);

		} /* end if */

		back=0.0;
		if(brat==1)
		{
			back=Parms[1];

		} /* end if */

		back1=1-back;
		if (bmdparm.risk==1)
		{
			back1=1;

		} /* end if */

		/* Values used to calculate BMDL and possibly BMDL values on the curve
		if specified by user */
		Gamma_BMD (nparm, Parms, EPS, &junk, xlk, Rlevel, Bmdl, &BMD);

		if (BMD >= 0.0)
		{
			OUTPUT_BENCHMD(1, BMD);
		}
		else
		{
			if (BMD == -1.0)
				fprintf(fp_out, "BMD computation failed\n");
			else
				fprintf(fp_out,
				"BMD computation failed.\nModel may be inappropriate for data\n");
		}
		if (Bmdl[1] < 0)
		{
#ifndef RBMDS
			fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  No);
#endif
			fprintf(fp_out,
				"           Benchmark dose computation failed.  Lower limit includes zero.");
		}
		else
		{
			if (fabs (BMDL_Error_Size) > 1e-3)
			{
				fprintf (fp_out, "           Warning: BMDL computation is at best imprecise for these data\n");
			}
#ifndef RBMDS
			fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  Yes);
#endif
			/*  	      printf("           BMDL =%14.6g\n\n", Bmdl[1]);   */
			fprintf(fp_out, 
#ifndef RBMDS
				"            BMDL =%14.6g\n\n"
#else
				"            BMDL =%30.22g\n\n"
#endif
				, Bmdl[1]);
			/* fprintf(fp_out, "Cancer Slope Factor =%14.6g\n\n", bmdparm.effect/Bmdl[1]); */
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
			{
				fprintf (fp_out2,"\n %f %f", Bmdl[i], Rlevel[i]*back1+back);

			} /* end for */
#endif
		}	/* end if */
	} /* end of if (bmdose == yes */

	/* indicate all required data have been output */
#ifndef RBMDS
	fprintf (fp_out2,"\n\n Check_result %d", Yes);
#endif
	/* free allocated memory */
	FREE_DVECTOR (Parms,1,nparm);
	FREE_IVECTOR (bounded, 1, nparm);
	FREE_DVECTOR (Ypp,1,Nobs);
	FREE_DVECTOR (Ep,1,Nobs);
	FREE_DVECTOR (Xi,1,Nobs);
	FREE_DVECTOR (Yp,1,Nobs);
	FREE_DVECTOR (Yn,1,Nobs);
	FREE_DVECTOR (IniP,1,nparm);
	FREE_IVECTOR (Spec,1,nparm);
	FREE_IVECTOR (IniSp,1,nparm);
	FREE_VLVECTOR(varsum,1,3);
	FREE_ALVECTOR(anasum,1,3);
	FREE_DVECTOR (Rlevel,1,5);
	FREE_DVECTOR (Bmdl, 1, 5);
	if (ErrorFlag == 0) {
		FREE_DMATRIX (vcv_adj, 1, adj_vcv_rows, 1, adj_vcv_rows);
	}
	FREE_DMATRIX(vcv,1,nparm,1,nparm);
	CLOSE_FILES ();
	return(0);

}	/* end of main */

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
*                     values
***********************************************************/
void Predict(double doses[], int ndoses, double Parms[], double P[])
{
	int i;
	for (i = 1; i <= ndoses; i++)
		if (doses[i] <= 0.0)
			P[i] = Parms[(int) Background];
		else
			P[i] = Parms[(int) Background] + (1.0 - Parms[(int) Background]) *
			GAMMP(Parms[(int) Power], Parms[(int) Slope] * doses[i]);
}


/***********************************************************
* unpack -- used to construct a full parameter vector out
*           of fixed and varying parts
*           input:
*                 x is the vector of "varying" values
*                 fixed is the vector of "fixed" values
*           output:
*                 p is the vector of parameters
***********************************************************/

void unpack(double x[], double fixed[], double p[])
{
	int j, jfixed, jvar;

	jfixed=jvar=0;
	/* reconstruct the parameter vector */
	for(j=1; j<=nparm; j++)
	{
		if(Spec[j]==Yes)
		{
			p[j]=fixed[jfixed];
			jfixed++;
		}
		else
		{
			p[j]=x[jvar];
			jvar++;
		}
	}
	/* if replace == Yes, then we are computing BMD; replace p[(int) Slope] */
	/* with function of current guess at BMD (i.e., tD) */
	if (replace==Yes)
	{
		if(bmdparm.risk==ADDED)
		{
			/* reparameterize parameter 2 for CI calculation using added risk */
			p[(int) Slope ] =
				XGAMMAI(BMR/(1 - p[(int) Background]),p[(int) Power])/tD;
		}
		else
		{
			/* reparameterize parameter 2 for CI calculation using extra risk */
			/* In Gamma_BMD ck is defined as follows, ck = -log(1-BMR) */
			p[(int) Slope ] = XGAMMAI(BMR,p[(int) Power])/tD;
		}
	}
}

/*********************************************************************** */
/* Which_Bounded -- Fills the 1-based vector bounded with 1 if the */
/*                  corresponding parameter is on one of its boundaries, 0 */
/*                  otherwise. */
/* */
/*                  Global: */
/*                    nparm, logtrans, Max_double, SlopeUpperBound*/
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
	for (i=1; i<=nparm; i++)
	{
		bounded[i] = Spec[i];
	} /* end for */

	if (VERYCLOSE(Parms[1],0.0) || VERYCLOSE(Parms[1],1.0)) bounded[1] = 1;
	if ((restrict == Yes && VERYCLOSE(Parms[3],1.0)) ||
		(restrict != Yes && VERYCLOSE(Parms[3],0.0)) ||
		VERYCLOSE(Parms[3], SlopeUpperBound)) bounded[3] = 1;
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
**Gamma_lk -- used to compute the log likelihood for Gamma model.
*
*		Global var.:
*			Nobs, Xi, Yp, Yn, replace, Spec, ck, bmdparm.risk
*		Input:
*			nvar is the number of parameters
*			nf is not currently used
*			uiparm is a vector of the indices for the parameters
*				of length nparm
*			urparm is a vector of fixed parameters
*				of length jfixed
*			ufparm is not currently used
*		Output:
*			f is the likelihood value
*		Input/Output:
*			x is a vector of non-fixed parameters
*				of length jvar
*
*********************************************************************/

void Gamma_lk(long int *nvar, double *x, long int *nf, double *f,
			  long int *uiparm, double *urparm, void (*ufparm)())
{

	int      i;
	double   xlk;     		    /* log likelihood */
	double   *p;
	double   *Pred;
	/* parameters for calculation */

	p = DVECTOR(1,nparm);
	Pred = DVECTOR(1, Nobs);

	/* construct the parameter vector */
	unpack(x, urparm, p);

	/* Compute the likelihood */
	xlk = 0.0;
	Predict(Xi, Nobs, p, Pred);

	for (i=1;i<=Nobs;i++)
	{
		xlk += Yp[i] * Slog(Pred[i]) +
			Yn[i] * Slog(1.0 - Pred[i]);
	}	/* end for */

	/* free allocated memory */
	FREE_DVECTOR(p, 1, nparm);
	FREE_DVECTOR(Pred, 1, Nobs);
	*f= -xlk;  /* f = -xlk is the maximum likelihood value
			   because DMNGB is a minimization routine */

}	/* end Gamma_lk */

/*******************************************************************
**Gamma_g -- used to compute the gradients for Probit model.
*		Extern var.: smean, smax, Nobs, Xi, Yp, Yn, Ls, Xg.
*
*		input:
*			nvar is the number of parameters
*			nf is something for the fortran maximization--not used
*			uiparm is a vector of the indices for the parameters
*				of length nparm
*			urparm is a vector of fixed parameters
*				of length jparm
*			ufparm is something for the fortran maximization--not used
*		output:
*			f is the likelihood value
*		input/output:
*			x is a vector of non-fixed parameters of length jvar
*
**********************************************************************/

/* This is a finite difference version */

void Gamma_g (long int *nvar, double *x, long int *nf, double *g,
			  long int *uiparm, double *urparm, void (*ufparm) ())
{
	double basefun, *saveparms, *h, tmp, hrat;
	int i, j;
	long int inf;


	saveparms = DVECTOR(0, nparm);
	h = DVECTOR(0, nparm);
	hrat = pow(1.0e-16,0.333333);

	for (i = 0; i < *nvar; i++)
	{
		if (fabs(x[i]) > DBL_EPSILON)
		{
			h[i] = hrat * fabs(x[i]);
			tmp = x[i] + h[i];
			h[i] = tmp - x[i];
		}
		else
			h[i] = hrat;
	}

	for (j = 0; j < *nvar; j++)
		saveparms[j] = x[j];
	for (i = 0; i < *nvar; i++)
	{
		if (i > 0) saveparms[i-1] = x[i-1];
		saveparms[i] = x[i] - h[i];
		Gamma_lk(nvar, saveparms, &inf, &basefun, uiparm, urparm, ufparm);
		saveparms[i] = x[i] + h[i];
		Gamma_lk(nvar, saveparms, &inf, &tmp, uiparm, urparm, ufparm);
		g[i] = (tmp - basefun) / (2.0*h[i]);

	}

	FREE_DVECTOR(saveparms,0,nparm);
	FREE_DVECTOR(h, 0, nparm);
}

/****************************************************************************
**
* Gamma_grad -- Computes the gradient of the Gamma likelihood
* function with respect to the user form of the parameters.  This is to
* be used in Gamma_vcv, to compute a finite difference approximation to
* the hessian of the likelihood function
* Input: nparm -- the number of parameters, the dimension of all the following
arrays.
Spec[] -- if the ith parameter is fixed by the user, Spec[i] == 1,
otherwise Spec[i] == 0.
ptf[] -- vector of parameters, in user form (external form),
based at 1.
* Output: grad[] -- the gradient of the loglikelihood function
(N.B. Gamma_g returns the gradient of -loglikelihood)
Based at 1.

*****************************************************************************/
void Gamma_grad(int nparm, int Spec[], double ptf[],
				double grad[])
{
	long int *uiparm, nvar, i, jfixed, jvar, nf;
	double *urparm, *start, *outgrad, *p;
	void (*ufparm) ();

	/* set up initial parameter array, start.  All parameters go either into */
	/* start (if they are changing to improve the fit) or urparm (if they are */
	/* fixed).*/

	p = DVECTOR(1, nparm);
	for (i = 1; i <= nparm; i++)
	{
		p[i] = ptf[i];
	}
	/* First, count the number of actually varying parameters */
	nvar = nparm;
	for (i = 1; i <= nparm; i++)
	{
		nvar = nvar - Spec[i];
	}
	uiparm = (long int *) malloc ((size_t) (nparm + 2) * sizeof (long int));
	urparm = DVECTOR(0, nparm - nvar);
	start = DVECTOR(0, nparm);
	outgrad = DVECTOR(0, nparm);

	jfixed = 0;
	jvar = 0;
	/* Separate the fixed and variable parameters */
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
	replace = No;
	Gamma_g (&nvar, start, &nf, outgrad, uiparm, urparm, ufparm);
	/* Gamma_g returns the gradient of -loglikelihood */
	for (i = 0; i < nvar; i++)
	{
		grad[i+1] = -outgrad[i];
	}
	FREE_DVECTOR(p, 1, nparm);
	free(uiparm);
	FREE_DVECTOR(urparm, 0, nparm - nvar);
	FREE_DVECTOR(start, 0, nparm);
	FREE_DVECTOR(outgrad, 0, nparm);
}


/****************************************************************************
*	Gamma_vcv -- used to compute the vcv matrix for Gamma model.
*
*		Extern var: Nobs, Xi, Yp, Yn, Ls, Xg.
*
*		Input:	nparm -- the number of parameters
*			Spec[] -- a vector indicating which parameters are
*				  user specified of length nparm
*			ptf[] -- parameter vector of length nparm
*
*		Output:	vcv -- will have the uninverted vcv matrix on exit
*			       of size nvar by nvar, where nvar is the
*                              number of "unfixed" parameters.
*
*
*****************************************************************************/

void Gamma_vcv (int nparm, int Spec[], double ptf[], double **vcv)
{
	double *saveparms, *h, *gradp, *gradm, hrat, tmp;
	int i, j, jvar;

	/* initialize memory for all the arrays */
	saveparms = DVECTOR(1, nparm);
	h = DVECTOR(1, nparm);
	gradp = DVECTOR(1,nparm);
	gradm = DVECTOR(1, nparm);
	/* Get a value of h for each parameter */
	hrat = pow(1.0e-16, 0.333333);

	for (i = 1; i <= nparm; i++)
	{
		if (fabs(ptf[i]) > 1.0e-7)
		{
			h[i] = hrat * fabs(ptf[i]);
			tmp = ptf[i] + h[i];
			h[i] = tmp - ptf[i];
		}
		else
			h[i] = hrat;
	}
	/* initialize saveparms with the parameter values */
	for (i = 1; i <= nparm; i++)
	{
		saveparms[i] = ptf[i];
	}
	/* for each unfixed parameter, compute the second derivative wrt each */
	/* of the others. */
	/* We need to return a full 3x3 matrix, regardless of what is fixed. */
	/* The irrelevant rows and columns will be deleted. */
	for (i = 1; i <= nparm; i++)
	{
		if (i > 1) saveparms[i-1] = ptf[i-1];
		if (Spec[i] != Yes)
		{
			saveparms[i] = ptf[i] + h[i];
			Gamma_grad(nparm, Spec, saveparms, gradp);
			saveparms[i] = ptf[i] - h[i];
			Gamma_grad(nparm, Spec, saveparms, gradm);
			/* Now compute the second derivative */
			jvar = 0;
			for (j = 1; j <= nparm; j++)
			{
				if (Spec[j] != Yes)
				{
					jvar++;
					vcv[i][j] = -(gradp[jvar] - gradm[jvar])/(2.0 * h[i]);
				}
			}
		}
	}

	FREE_DVECTOR(saveparms,1,nparm);
	FREE_DVECTOR(h,1,nparm);
	FREE_DVECTOR(gradp,1,nparm);
	FREE_DVECTOR(gradm,1,nparm);
} /* end of Gamma_vcv */

/**************************************************************
*MAX_lk -- used to obtain the maximum log-likelihood as well as
*           the estimates of parameters, given initial p[1..nparm],
*           object func. , and gradient func. G_func.
*
*		Global var.:
*			Spec, Yp, Min_increment, Max_double, restrict
*		Input:
*			nparm is the number of parameters in the model
*			gtol is the tolerance level
*			iter is the number of iterations
*		Output:
*			fret is the value of the log-likelihood function
*		Input/output:
*			p[] is the parameter vector of length nparm
*				upon output, has the maximum likelihood estimates
*
**************************************************************/
void MAX_lk(int nparm, double p[], double gtol, int *iter, double *fret)

{
	int		i, jfixed, jvar;
	long int		nvar;
	long int	*uiparm;

	double	*start;
	double	*urparm; /* I've hardcoded the maximum number of
					 parameters (3) here */
	double	lower[10], upper[10]; /* hard coded the max number */
	/* of parameters to 3 here */
	void		(*ufparm)();

	/* Set up initial parameter array, start.  All parameters go either
	into start (if they are changing to improve the fit) or urparm
	(if they are fixed).  Set up all elements of Dscale to be 1.0.*/

	nvar = nparm;
	for (i = 1; i <= nparm; i++)
	{
		nvar = nvar - Spec[i]; /* Count the varying parameters */

	} /* end for */

	if (nvar > 0)
	{
		urparm = (double *) malloc((size_t)(nparm-nvar)*sizeof(double));
		start = (double *) malloc((size_t) nvar*sizeof(double));
	}
	else
	{
		long int dummy;
		urparm = (double *) malloc((size_t)(nparm-nvar)*sizeof(double));
		for(i=1;i<=nparm;i++)
		{
			urparm[i-1]=p[i];
		} /* end for */

		/* Nothing to do...calculate likelihood then go back */

		Gamma_lk( &nvar, start, &dummy, fret,
			uiparm, urparm,  ufparm);
		*fret= -(*fret);
		ErrorFlag = -1; /* so I can test in main */
		return;

	} /* end if */


	jfixed = 0;
	jvar = 0;
	for (i = 1; i <= nparm; i++) /* separate the fixed and variable
								 parameters */
	{
		if (Spec[i] == 1)
		{
			urparm[jfixed] = p[i];
			jfixed++;
		}
		else
		{
			start[jvar] = p[i];
			jvar++;

		} /* end if */

	} /* end for */

	/* set up the bounds variable.  Each parameter is unique, here. */
	jvar = 0;
	if (!fixedParm(Background))
	{
		lower[jvar] = 0.0;
		upper[jvar] = 1.0;
		jvar++;

	} /* end if (Spec[1] != 1) */

	if (!fixedParm(Slope))
	{
		lower[jvar] = 0.0;
		upper[jvar] = Max_double;
		jvar++;

	} /* end if */

	if (!fixedParm(Power))
	{
		if (restrict == 1)
		{
			lower[jvar] = 1.0;
		}
		else
		{
			lower[jvar] = 0.0;

		} /* end if */

		upper[jvar] = SlopeUpperBound;

	} /* end if (!fixedParm(Power)) */
	ErrorFlag = run_dmngb((int) nvar, start, lower, upper, Maxloglik,
		Rel_Conv, Parm_Conv, ITMAX, 10,
		Gamma_lk,Gamma_g,uiparm,urparm,ufparm,
		DeBuG,fret);

	*fret = -*fret;

	/* put the parameter values back in p */
	unpack(start, urparm, p);
	/* free allocated memory */
	free(urparm);
	free(start);

}	/* end MAX_lk */

/****************************************************************** */
/* find_A -- function for zeroin to find the shape parameter */
/****************************************************************** */

double find_a(int n, double p[], double a, double gtol)
{
	return XGAMMAI(p[1],a)/p[2] - XGAMMAI(p[3],a)/p[4];
}

/**************************************************************
*Gamma_fit -- Used to "prepare" the data for further computation,
*            i.e. compute the extern variables, give the initial
*            parameters, etc. THEN fit the Gamma model.
*            (In fact, these jobs could be done in main().)
*
*		Global var.:
*			Nobs, scale, initial, IniP, Spec, Yp, Yn, Xi, replace
*		Input:
*			nparm is the number of parameters in the model
*			gtol is the tolerance level
*			iter is the number of iterations
*		Output:
*			fret is the value of the log-likelihood function
*		Input/Output:
*			p[] is the parameter vector of length nparm
*				upon output, has the maximum likelihood estimates
*
***************************************************************/
void Gamma_fit(int nparm, double p[], double gtol,
			   int *iter, double *fret)
{
	int	  i, j, junk;
	double  ymin, W, x, xlk;

	/* Find the minimum and maximum dose values */
	ymin = 1.0;
	xmin = 1000000;
	xmax = 0.0;

	for (i=1;i<=Nobs;i++)
	{
		x=Xi[i];
		if (x < xmin)
		{
			xmin = x;

		} /* end if */

		if (x > xmax)
		{
			xmax = x;

		} /* end if */

	} /* end for */
	if (allFixed())
	{
		scale = 1.0;
	}
	else
	{
		scale=1.0/xmax;
	}
	/* rescale Dose to be: 0 <= Dose <= 1 */
	for (i = 1; i <= Nobs; i++) Xi[i] *= scale;

	/* obtain user-supplied initial estimates for p[] */
	if(initial==Yes)
	{
		for(j=1; j<=nparm; j++)
		{
			if (!fixedParm( (Model_Parms) j))
				p[j]=IniP[j];

		} /* end for */
		OUTPUT_TEXT("\n\n                  Initial (and Specified) Parameter Values  ");
		OUTPUT_Init(nparm, Spec, p, Parm_name);
		/* If these are user-supplied, we have to rescale the Slope */
		p[(int) Slope ] /= scale;
	}
	else
	{
		/* compute initial estimates */
		int contdose, maxdose, bri;
		double ymax, pihat, Dhat;
		/* If the slope is fixed, we have to rescale it, anyway */
		if (fixedParm(Slope))
		{
			p[(int) Slope] /= scale;
		}
		/* Look for the index for the smallest dose, and the largest dose.
		We need to do something a little differently if xmin != 0, but
		first we find the index where Xi[i] == xmin/xmax.
		*/
		contdose = maxdose = 0;
		for (i = 1; i <= Nobs; i++)
		{
			if (fabs(Xi[i] - xmin/xmax) < 1.0e-10)
			{
				contdose = i;
			}
			if (fabs(Xi[i] - 1.0) < 1.0e-10)
			{
				maxdose = i;
			}
		}
		ymax = (Yp[maxdose]+1)/(Yp[maxdose] + Yn[maxdose] + 2);

		if (!fixedParm(Background))
		{
			p[(int) Background ] = (Yp[contdose]+1)/
				(Yp[contdose] + Yn[contdose] + 2);
			/* If xmin > 0, then this is really the response at the lowest non-control
			dose; just halve this estimate.
			*/
			if (xmin > 0) p[(int) Background ] *= 0.5;
		}

		/* There are three possibilities: if one of Power and Slope are fixed */
		/* calculate the other to go through the half-max point.  If neither are */
		/* fixed, compute both so the curve goes through the half-max point and the */
		/* max point. */
		/* Here, compute the half-max point, (Dhat, pihat) */

#define YP(i) ((Yp[i] + 1) / (Yp[i] + Yn[i] + 2))

		pihat = 0.5 * (p[(int) Background ] + YP(maxdose));

		for (i = 1; i < Nobs; i++)
		{
			if (YP(i) <= pihat && YP(i+1) >= pihat)
			{
				bri = i;
				break;
			}
		} /* end for */
		Dhat = Xi[bri] +
			(Xi[bri + 1] - Xi[bri]) *
			((pihat - YP(bri)) / (YP(bri + 1) - YP(bri)));

		if (fixedParm(Power) && !fixedParm(Slope))
		{
			W = (pihat - p[(int) Background])/(1.0 - p[(int) Background]);
			p[(int) Slope] = XGAMMAI(W, p[(int) Power])/Dhat;
		}
		if (!fixedParm(Power) && fixedParm(Slope))
		{
			W = (pihat - p[(int) Background])/(1.0 - p[(int) Background]);
			p[(int) Power] = XGAMMAI_A(W, p[(int) Slope] * Dhat);
		}
		if (!fixedParm(Slope) && !fixedParm(Power))
		{
			double *fittedpoints, ax, bx;
			fittedpoints = DVECTOR(1,4);
			fittedpoints[1] = (pihat - p[(int) Background])/
				(1.0 - p[(int) Background]);
			fittedpoints[2] = Dhat;
			fittedpoints[3] = (YP(maxdose) - p[(int) Background]) /
				(1.0 - p[(int) Background]);
			fittedpoints[4] = 1.0; /* we are working on a scaled dose range */
			if (restrict == Yes)
			{
				ax = 1.0;
			}
			else
			{
				ax = DBL_EPSILON;
			}
			bx = 18;
			if (find_a(4,fittedpoints, ax, 1.0) *
				find_a(4,fittedpoints, bx, 1.0) > 0.0)
			{
			  /* Somethings amiss here; the following line is required to keep
			     the program from hanging.  Memory problem?  */
			  fprintf(stderr,"                                      \n");
			  p[(int) Power] = 1.3;
			  p[(int) Slope] =
			    XGAMMAI(fittedpoints[1],p[(int) Power])/fittedpoints[2];
			}
			else
			{
			  p[(int) Power] =
			    zeroin(ax,bx,1e-10,find_a,4,fittedpoints,1.0);
			  p[(int) Slope] =
			    XGAMMAI(fittedpoints[1],p[(int) Power])/fittedpoints[2];
			}
			FREE_DVECTOR(fittedpoints,1,4);
		}


		/* Now do something silly, so that we can print out the Slope */
		/* on the right scale */
		if (!fixedParm(Slope)) p[(int) Slope] *= scale;
		OUTPUT_TEXT("\n\n                  Default Initial (and Specified) Parameter Values  ");
		OUTPUT_Init(nparm, Spec, p, Parm_name);
		if (!fixedParm(Slope)) p[(int) Slope] /= scale;

	} /* end if (initial == Yes) */

	/* Fit the model. */
	replace = No;

	/* I deleted some code that repeated the fit 5 times, */
	/* perturbing the final point */
	/* a little bit each time (run_dmngb now does this internally)*/
	MAX_lk( nparm, p, gtol,  &junk,  &xlk);

	do_dmngb_warning(&ErrorFlag);

	*fret = xlk;
	/* rescale everything, including the slope */
	for (i = 1; i <= Nobs; i++) Xi[i] /= scale;
	p[(int) Slope] = p[(int) Slope] * scale;
	scale = 1.0;

}	/* end Gamma_fit */

/************************************************************
* Gamma_BMD -- Used to calculate the BMD and BMDL for Gamma model.
*
*		Global var.:
*			bmdparm.level, LR, Spec, replace, bmdparm.risk, bmdparm. effect
*			ck, xmax, Max_double, BMD_lk, bmdlCurve
*		Input:
*			nparm is the number of parameters
*			p[] is the parameter vector of length nparm
*			gtol is the tolerance level
*			iter is the number of iterations
*
*		Output:
*			Rlevel[] is the risk level
*			Bmdl[] is the benchmark dose confidence limit
*			BMD is the benchmark dose
*
*		Input/Output:
*			xlk is the likelihood value of the function
*
*************************************************************/
void Gamma_BMD (int nparm, double p[], double gtol, int *iter, double xlk,
				double Rlevel[], double Bmdl[], double *BMD)
{
	double   tol;
	double   xa,xb,fa,fb, lastfa;
	double   stepsize;
	double   *pBak, *pa, *pb;
	int      i, j, k, bogusBMD;
	double *pred, *doses, effect, BMDtmp;


	pBak=DVECTOR(1, nparm);
	pa = DVECTOR(1, nparm);
	pb = DVECTOR(1, nparm);

	/* This is where stepsize is assigned */

	stepsize = 0.9;

	for(j=1; j<=nparm; j++)
	{
		pBak[j]= p[j];  /* save the p[] */

	} /* end for */

	/* compute X^2 value */
	if (bmdparm.level<0.5)
	{
		LR = QCHISQ(1.0-2*bmdparm.level,1)/2.0;
	}
	else
	{
		LR = QCHISQ(2*bmdparm.level-1.0,1)/2.0;

	} /* end if */

	Rlevel[1] = BMR = bmdparm.effect;

	/***** Solve the BMD ******************************************/

	/* If there is no slope, we want to cap the BMD estimate. */
	/* We still want to compute a lower limit on the BMD, however */
	bogusBMD = 0;

	if (p[3] <= 0.0 || p[(int) Slope] <= 0.0)
	{
		*BMD = xb = 100 * xmax;
		bogusBMD = 1;
		Warning(" Slope or Power parameter essentially zero.  BMD set to 100 * max(Dose).\n");

	}
	else
	{
		if (bmdparm.risk == ADDED)
		{
			xb = XGAMMAI(BMR/(1.0 - p[(int) Background]),
				p[(int) Power])/p[(int) Slope] ;
		}
		else
		{
			xb = XGAMMAI(BMR,p[(int) Power])/p[(int) Slope] ;
		}
		*BMD = xb;
	}
	/* Make sure the response at the BMD is correct */
	pred = DVECTOR(1,4);
	doses = DVECTOR(1,4);
	doses[1] = 0.0;
	doses[2] = *BMD;
	Predict(doses, 2, p, pred);
	if (bmdparm.risk == ADDED)
	{
		effect = pred[2] - pred[1];
	}
	else
	{
		effect = (pred[2] - pred[1])/(1.0 - pred[1]);
	}
	if (!bogusBMD && (fabs(effect - BMR) > 1.0e-3))
	{
#ifndef RBMDS
		fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  No);
#endif
		fprintf(fp_out,
			"\nComputed BMD is %g; response at the computed BMD is %g\n",
			*BMD,pred[2]);
		fprintf (fp_out,
			"Computed effect at this estimate is: %g, requested \
			effect is %g\n",effect, BMR);
		*BMD = -1.0;
		return;
	}
	else
	{
		if (bogusBMD && effect > BMR)
		{
			fprintf(fp_out,
				"\nComputed BMD is %g; response at the computed BMD is %g\n",
				*BMD,
				pred[2]);
			fprintf (fp_out,
				"Computed effect at this estimate is: %g, requested \
				effect is %g\n",effect,BMR);
#ifndef RBMDS
			fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  No);
#endif
			*BMD = -2.0;
			return;
		}
	}
	FREE_DVECTOR(pred,1,4);
	FREE_DVECTOR(doses,1,4);

	/* OK -- We have a legitimate BMD */

	/********* search for BMDL **************************/

	Spec[(int) Slope]=1;
	replace=Yes;

	xa = xb*stepsize;

	tol = FMAX((*BMD)*0.000001, 0.0000001);

	BMD_lk = xlk;				/* get the lk at BMD */
	fb = -LR;
	for (i = 1; i <= nparm; i++) pa[i] = pb[i] = p[i];
	fa = BMDL_func(nparm, pa, xa, tol);  /* if fa > 0, then there is a BMDL
										 somewhere in the range */
	lastfa = DBL_MAX;
	while (fa<0.0 && xa > DBL_MIN && fabs(fa - lastfa) > DBL_EPSILON)
	{
		xb = xa;
		fb = fa;
		for (i = 1; i <= nparm; i++) pb[i] = pa[i];
		xa *= stepsize;
		lastfa = fa;
		fa = BMDL_func(nparm, pa, xa, tol);
	}
	if (fa<0.0)
	{
		/* computation failed */
		Bmdl[1] = -1.0;
		return;
	} /* end if */
	else
	{
		Bmdl[1] = zeroin(xa, xb, 1.0e-10, BMDL_func, nparm, pb, 1.0e-14);
		BMDL_Error_Size = BMDL_func(nparm, pb, Bmdl[1], tol);
	}

	if(bmdlCurve==Yes)
	{
		/* calculate Bmd[] and Bmdl[] */
		for (k=2; k<=5;k++)
		{
			/* 1: reload the original parameter vector */
			for(j=1; j<=nparm; j++)
			{
				p[j]= pBak[j];

			} /* end for */
			/* 2:  Set up the BMR levels to use (may want to change this code) */

			if (k==2)
			{
				Rlevel[k]=BMR=0.05;
			}
			else
			{
				Rlevel[k]= BMR = (k-2)*0.1;

			} /* end if */


			/* 3: solve the BMD[] */
			if (p[3]<= 0 || p[(int) Slope] <= 0)
			{
				BMDtmp = 100 * xmax;
				xb = BMDtmp;
			}
			else
			{
				if (bmdparm.risk == ADDED)
				{
					xb = XGAMMAI(BMR/(1.0 - p[(int) Background]),
						p[(int) Power])/p[(int) Slope] ;
				}
				else
				{
					xb = XGAMMAI(BMR,p[(int) Power])/p[(int) Slope] ;
				}
				BMDtmp = xb;
			}
			/* First get the interval to search over */
			xa = xb * stepsize;
			tol = FMAX(xb*0.0001, 0.0000001);

			BMD_lk = xlk;		/*get the lk at BMD.*/
			fb = -LR;
			for (i = 1; i <= nparm; i++) pa[i] = pb[i] = p[i];
			fa = BMDL_func(nparm, pa, xa, tol);
			lastfa = DBL_MAX;

			while (fa < 0.0 && xa > DBL_MIN && fabs(fa - lastfa) > DBL_EPSILON) /* BMDL does not lie between xa and xb, so reduce xa */
			{
				xb = xa;
				fb = fa;
				for (i = 1; i <= nparm; i++) pb[i] = pa[i];
				xa *= stepsize;   /*prevent that xa=stepsize*BMD is not small eno*/
				lastfa = fa;
				fa = BMDL_func(nparm, pa, xa, tol);
			}
			if (fa<0.0)
			{
				/*computation failed*/
				/*fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  No); */
				Bmdl[k]=-1;
				fprintf(fp_out, "\n BMDL curve computation failed for BMR = %f . \n The BMDL curve appearing in the graph may not be accurate.", BMR);
			}
			else
			{
				/*computation will succeed.*/
				/*fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  Yes); */
				Bmdl[k] = zeroin(xa, xb, 1.0e-10, BMDL_func, nparm, pb, 1.0e-14);
			}
		}
	}
	else
	{
		for (k=2; k<=5;k++)
		{
			Bmdl[k] = Rlevel[k]= -1;
		}
	}

	/* free allocated memory */
	FREE_DVECTOR(pBak,1, nparm);
	FREE_DVECTOR(pa,1, nparm);
	FREE_DVECTOR(pb,1, nparm);

}	/* end Gamma_BMD */


/*****************************************************************
* BMDL_func -- used to compute the values of functions BMDL_f (the
*              X^2 value) at the point D, given the parm p[] and
*              number of parm.
*
*              This routine is called by Binary_root().
*
*		Global var.:
*			tD, ck, BMD_lk, LR
*		Input:
*			nparm is the number of parameters
*			pBak[] is the parameter vector
*			D is the point where we will compute the BMDL
*			gtol is the tolerance level
*
*		Returns:
*			 Likelihood at BMD - Likelihood at BMDL - Chi-Square
*
*****************************************************************/
double BMDL_func(int nparm, double p[], double D, double gtol)
{
	/* ck , BMD_lk and LR are calculated in Gamma_BMD() */

	double fD, xlk;
	int junk;
	tD = D; /* tD is global var. have to change before call MAX_lk().*/

	MAX_lk( nparm, p, gtol,  &junk,  &xlk);

	fD = BMD_lk - xlk - LR;  /* constraint imposed by the likelihood function */

	return fD;

}   /* end BMDL_func */

