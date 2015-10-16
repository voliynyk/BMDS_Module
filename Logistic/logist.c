/****************************************************************
**
* Important NOTE:  The following variable is the version number for
*                  the current model.  THIS MUST BE CHANGED as
*				   important changes are made to the models.
*
*****************************************************************/

char Version_no[] = "Logistic Model. (Version: 2.14; Date: 2/28/2013)";

/****************************************************************
char Version_no[] = "Logistic Model. (Version: 2.13; Date: 10/28/2009)";
char Version_no[] = "Logistic Model. (Version: 2.12; Date: 05/16/2008)";
char Version_no[] = "Logistic Model. (Version: 2.11; Date: 03/10/2008)";
char Version_no[] = "Logistic Model. (Version: 2.10; Date: 09/23/2007)";
char Version_no[] = "Logistic Model. (Version: 2.9; Date: 02/20/2007)";
char Version_no[] = "Logistic Model. (Version: 2.7; Date: 12/08/2006)";
char Version_no[] = "Logistic Model $Revision: 2.2 $ $Date: 2001/07/11 11:53:20 $";
*****************************************************************/

/****************************************************************
**
* Logist.C - a ANSI C program for Logist model fitting with/without
*             a natural background rate in Benchmark Dose.
*
* Date: Feb  23, 1998
*
********************************************************************
* Modification Log:
*
*
* Version Number: 2.1
* Modified By: R. Woodrow Setzer
* Modified Date: 02/26/2000
* Reason:
*
* Version Number: 2.3
* Modified By: Qun He
* Modified Date: 09/03/2003
* Reason:
*
* Version Number: 2.4
* Modified By: Micheal Ferree
* Modified Date: 08/13/2005
* Reason: Took out print out of cancer slope factor
*
* Version Number: 2.5
* Modified By: Woodrow Setzer
* Modified Date: 09/23/2005
* Reason: Free all allocated memory prior to exiting
*
* Version Number: 2.6
* Modified By: Woodrow Setzer
* Modified Date: 03/22/2006
* Reason: Replace calls to MATINVS with INVMAT
*
* Version Number: 2.7
* Modified By: Woodrow Setzer
* Modified Date: 12/08/2006
* Reason: minor correction to internal comment: original was copied from
*         NLogistic code without changing the associated comments
*
*
* Version Number: 2.8
* Modified By: Geoffrey
* Date: 1/12/2007
* Reason: Incremented version number.
*		 Added last parameter "1" (print SE) in OP_ParmsE() if log-logistic.
*
* Version Number: 2.9
* Modified By: Woodrow Setzer
* Date: 2/20/2007
* Reason: Incremented version number to reflect changed compilation options.

* Version Number: 2.10
* Modified By: G. Nonato
* Modified Date: 9/23/2007
* Reason:	1. Changed version number to reflect changed.
*			2. Prefix condition with "Xi[1] == 0.0 &&" in Line 470 and beaufify code {GLN - 09/23/07, PR0823-05} 
			It fixes the following:
				Different incorrect responses result when running models on data for
				which there is no zero dose group and background is specified to be 0.
				For the  example data set where doses are 17, 20 and 24, Ns are all 10
				and responses are 1, 2 and 6

* Version Number: Version: 2.11
* Modified By: G. Nonato
* Modification Date: 03/10/2008
* Reason: (Per BMDS 2.0: Problem Report 157 & 147)
*       Fix the Observation # < parameter # for Weibull model problem.
*       Added code to free-up allocated memories before exiting thru ERRORPRT()
*
* Version Number: 2.12
* Modified By: G. Nonato
* Modification Date: 05/16/2008
* Reason: (Per BMDS 2.0: Problem Report 165)
*       Goodness of Fit - Observed column, print values as real numbers.
*
* Version Number: 2.13
* Modified By: G. Nonato
* Modification Date: 10/28/2009
* Reason:
*      To be able to process files/folders with spaces (PR 257)
*      Fix program freeze due to long variable names (PR 278)
*      Process long files up to 256 characters (PR 303 and 308)
*      Modify code for easy maintenance, all lengths for file name,
*        model name, and column names are placed in benchmark.h
*
* Version Number: 2.14
* Modified By: Louis Olszyk
* Modification Date: 02/28/2013
* Reason: PR 444 - Fix wording in plot titles
*********************************************************************************/

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

void Logist_fit (int nparm, double *p, double gtol, int *iter, double *fret);
void Logist_vcv (int nparm, int *Spec, double *p, double **vcv);
void Which_Bounded (int [], double [], int []);
int Model_DF (int []);
void Logist_BMD (int nparm, double *p, double gtol, int *iter, double xlk,
				 double *Rlevel, double *Bmdl, double *BMD);
double BMDL_func (int nparm, double *pBak, double D, double gtol);



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
#define FREEALL FREE_DVECTOR (xi,1,n); FREE_DVECTOR (pnew, 1,n);\
	FREE_DVECTOR(hdg,1,n); FREE_DVECTOR(g,1,n);\
	FREE_DVECTOR(dg,1,n); \
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

/*** Define input and output files's name  *********************/
char fin[FLENGTH];				/* input temp file */
char fout[FLENGTH];				/* output temp file */
char fout2[FLENGTH];
char plotfilename[FLENGTH];		/* file to send to GnuPlot */
char *Parm_name[]={"background", "intercept", "slope"};
/* define a typedef enum for parameter names */
/* Base1 is because all the parameter vectors, spec, etc, are based at 1 */
typedef enum {Base1, Background, Intercept, Slope} Model_Parms;
char *anatxt[]={"Full model", "Fitted model", "Reduced model"};



/* variables will not be changed execpt Spec */
int    *Spec;    /* vector used to identify user input parm */
int    *IniSp;	 /* user initial specified parameter vector */
double *Yp;      /* positive dependent variable data array */
double *Yn;      /* negative dependent variable data array */
double *Xi;      /* independent variable data array */
double *Ypp;     /* predicted positive depdendent values */
double *Ep;      /* estimated probability */
double *Rlevel;	 /* risk level (BMR level) */
double *Bmdl;	 /* benchmark dose confidence limit */
double *IniP;    /* user initialized parameter vector */
double Rel_Conv, Parm_Conv, Maxloglik;
double SlopeUpperBound = 18.0;

int   Nobs, nparm, restrict, initial, appendix, smooth; /*inputs and options*/
int		bmdlCurve,logtrans; /* flags for bmdl curve and log
							transformation options */
double xmax, xmin, scale, Min_increment, Max_double;

double BMDL_Error_Size;
int BMDL_Error;

/* changing variable */
int	 replace, brat;
double tD,  BMD_lk, LR, ck, upb=18, BMR;

int ErrorFlag; /* Error States from run_dmngb */

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
	}
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
** main--main function used to call Logist mode fitting program.
Includes: biosubcc.c--common subfunction C program.
*****************************************************************/
int main (int argc, char *argv[])
{

	int    iter,i, junk;	/*iteration variable*/

	int    bmdose;		/*flag for computing benchmark dose*/
	int    Nmiss;			/*number of records with missing values*/
	int    nparm_known;	/*number of specified parameters */
	double  lkf,lkr,xlk,W;  /*log likelihoods */
	double  back, BMD, back1, junk2;

	double  *Parms;		/*parameter array*/
	VarList *varsum;		/*info for variables--p. dep.,n. dep., indep.*/
	AnaList *anasum;		/*information for ANOVA analysis*/
	double  **vcv;		/*variance and covariance matrix*/
	double **vcv_adj;		/* adjusted vcv matrix taking into account
							/  any bounded parameters */
	char long_path_name[FLENGTH];
	char   model_name[MNLENGTH], user_note[UNLENGTH];
	char  dose_name[CNLENGTH], posi_name[CNLENGTH], nega_name[CNLENGTH], junkname[FLENGTH];
	int *bounded;
	int adj_vcv_rows;
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
		/*        printf("Error in opening input  file.\n"); */
		/*        printf ("...now exiting to system...\n"); */

		fprintf(stderr,"Error in opening input file.\n");
		fprintf (stderr,"...Exited to system!\n");
		exit (1);
	}


	/* begin reading input from batch file (.(d) ext. */
	fscanf(fp_in, "%s", model_name);
	fscanf(fp_in, "%[ ^\n]", user_note);
	fscanf(fp_in, "%[^\n]", user_note);
	fscanf(fp_in, "%s", junkname);
	fscanf(fp_in, "%s", junkname);
	fscanf(fp_in, "%d",&Nobs);

	/*assign number of parameters*/
	nparm = 3;

	/*allocate memory for arrays*/
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


	/*Read more values from input file*/
	fscanf(fp_in,"%d%lf%lf%d%d%d%d%d%d", &ITMAX, &Rel_Conv, &Parm_Conv,
		&bmdlCurve,&logtrans, &restrict, &bmdose, &appendix, &smooth);
	fscanf(fp_in,"%lf%d%lf",&bmdparm.effect,&bmdparm.risk,&bmdparm.level);
	/* bmdparm.effect is the benchmark risk level */
	/* bmdparm.risk is 0 is extra risk, 1 if added risk */
	/* bmdparm.level is the confidence level */

	junk = 0;    /*Used to see if an extension was added to output file
				 name*/

	Get_Names(argv[1], fout, fout2, plotfilename);

	/* open output files */
	if(appendix==Yes)
	{
		fp_out=fopen(fout,"a");	    /* append output */
	}
	else
	{
		fp_out=fopen(fout,"w");	    /* overwrite output */
	}
#ifndef RBMDS
	fp_out2=fopen(fout2,"w");	    /* always overwrite plotting file */
#endif
	/* check to make sure files are open, if not, print error message and
	exit */
	if (fp_out==NULL 
#ifndef RBMDS
		|| fp_out2==NULL 
#endif
		)
	{

		fprintf(fp_out,"Error in opening output files.\n");
		fprintf (fp_out,"...Exited to system!\n");
		exit (1);
	}

	/* Print model and file information on output page */
	Output_Header(Version_no, argv[1], plotfilename, ctime(&ltime), user_note);

	if (bmdose < 0 || bmdose > 1)
	{
		FreeUp_mem(Parms, varsum, anasum, vcv);
		ERRORPRT("Error in choosing benchmark dose computation.");
	}

	/*obtain user input parameters*/
	READ_PARAMETERS(nparm,Parms);
	FILL_SPECVECTOR(nparm,Parms,Spec);
	if (logtrans==No)
	{
		if (fixedParm(Background))		/* check user specified values */
		{
			if (fixedParm(Intercept))
			{
				Warning("Warning: Both background and intercept were specified. The specified intercept was over written by corresponding specified background. \n");
			}

			if (Parms[(int) Background]>=1)
			{
				ERRORPRT("Background should be in the range of (0, 1).");
			}

			if (Parms[(int) Background ]<=0.0001)
			{
				Parms[(int) Intercept ]= -9.25;
			}
			else
			{
				Parms[(int) Intercept] = -log(1/Parms[1]-1);
			}
			Spec[(int) Intercept] = 1;
		}
		/* end if (Spec[1] ==1) */

		Spec[(int) Background ]=1;
		Parms[(int) Background ]=0;

	}
	/* end if (logtrans == 0) */

	nparm_known = COUNT_SPECVECTOR(nparm,Spec);
	brat=Yes;
	if (fixedParm(Background) && Parms[(int) Background]<EPS)
	{
		brat=No;
	}

	/*obtain user input initial parameters values      *******/
	fscanf(fp_in,"%d", &initial);
	READ_PARAMETERS(nparm,IniP);
	FILL_SPECVECTOR(nparm,IniP,IniSp);

	SWAP(Parms[(int) Intercept], Parms[(int) Slope], junk2); /*The program was reading slope as intercept and     */
	SWAP(Spec[(int) Intercept], Spec[(int) Slope], junk);    /*intercept as slope from the input page.  This      */
	/*will alleviate that without having to drastically
	alter this code or the UI code                     */
	SWAP(IniP[(int) Intercept], IniP[(int) Slope], junk2);

	/*obtain observation data into Yp, Yn, Xi vectors
	Yp = number of respondants
	Yn = number of nonrespondants
	Xi = dose level*/
	fscanf(fp_in,"%s%s%s", dose_name, posi_name, nega_name);
	Nmiss = READ_OBSDATA3V(Nobs,3,2,3,1, Yp,Yn, Xi);


	Nobs -= Nmiss;             /* extern variable Nobs has been changed
							   only use complete observations */

	//if (Nobs < nparm), commented and changed to the code below, GLN 3/7/08
	if (Nobs < (nparm-nparm_known))
	{
		FreeUp_mem(Parms, varsum, anasum, vcv);
		ERRORPRT("Observation # < parameter # for Logist model.");
	}

	/* end of input data */

	/*output title and summary of input data */
	OUTPUT_TEXT("\n   The form of the probability function is: ");
	if (logtrans==Yes)
	{
		OUTPUT_TEXT("\n   P[response] = background+(1-background)/[1+EXP(-intercept-slope*Log(dose))]");
	}
	else
	{
		OUTPUT_TEXT("\n   P[response] = 1/[1+EXP(-intercept-slope*dose)]");
	}

	fprintf(fp_out,"\n\n   Dependent variable = %s", posi_name);
	fprintf(fp_out,"\n   Independent variable = %s", dose_name);
	if(logtrans==Yes)		/* if log(dose) option is selected */
	{
		if (Spec[(int) Background]==Yes)	/* if background parameter is specified */
		{
			if (Parms[(int) Background] <= 0.0000001)
			{
				fprintf(fp_out,"\n   Background parameter is set to zero");
			}
			else
			{
				fprintf(fp_out,"\n   Background parameter is set to %g", Parms[(int) Background]);
			}
		}
	}

	if(Xi[1] == 0.0 && fixedParm(Background) && Parms[(int) Background] == 0.0 && Yp[1] != 0
		&& logtrans != No)		//Prefix condition with "Xi[1] == 0.0 &&" {GLN - 09/23/07, PR0823-05}
	{
		FreeUp_mem(Parms, varsum, anasum, vcv);
		ERRORPRT("ERROR:  Background parameter specified as 0, but there were responses\n        at the control dose\n");

	} /* end if */


	if (fixedParm(Intercept))		/* if intercept parameter is specified */
	{
		fprintf(fp_out,"\n   Intercept parameter is set to %g", Parms[(int) Intercept]);
	}

	if (fixedParm(Slope))		/* if slope parameter is specified */
	{
		fprintf(fp_out,"\n   Slope parameter is set to %g", Parms[(int) Slope]);
	}
	else
	{
		if (restrict==Yes)	/* if slope parameter restricted to >= 1 */
		{
			fprintf(fp_out,"\n   Slope parameter is restricted as slope >= 1");
		}
		else
		{
			fprintf(fp_out,"\n   Slope parameter is not restricted");
		}
	}
	/* end output of specified parameters */

	if(logtrans==No)
	{
		nparm_known -= 1;
	}

	fprintf (fp_out, "\n\n   Total number of observations = %d",Nobs+Nmiss);
	fprintf (fp_out, "\n   Total number of records with missing values = %d",Nmiss);

	fprintf(fp_out, "\n   Maximum number of iterations = %d\n", ITMAX);
	fprintf(fp_out, "   Relative Function Convergence has been set to: %g\n", Rel_Conv);
	fprintf(fp_out, "   Parameter Convergence has been set to: %g\n\n", Parm_Conv);

	if(logtrans == Yes)
	{
		fprintf(fp_out, "\n\n   User has chosen the log transformed model\n");
	}


	if(initial==Yes)
	{
		OUTPUT_TEXT("\n\n                 User Inputs Initial Parameter Values  ");
		OUTPUT_Init(nparm, Spec, IniP, Parm_name);

		for (i=1; i<=nparm; i++)
		{
			if(IniSp[i]==1)		/* have been initialized */
			{
				if(Spec[i]==1 )		/* check if it is for fixed parm */
				{
					if(i > 1)
					{
						Warning("The initial value for the fixed parameter is ignored.");
					}
					else
					{
						if(logtrans==Yes)
						{
							Warning("The initial value for the fixed parameter is ignored.");
						}
					}
				}
			}
			else
			{
				/*check if all the unspecified parms were initialized.*/
				if (Spec[i]==0)
				{
					ERRORPRT("When the initial option is chosen, one has to initial ALL unspecified parameters.");
				}
			}
		}

		for(i = 1; i <= 3; i++)
		{
			if(Spec[i] == 1)
			{
				IniP[i] = 1.0;     /* Make sure the next if statement doesn't
								   bomb when some parameters are specified*/
			}
		}

		if (IniP[1] < 0 || IniP[1]>1 || IniP[3]<restrict )
		{
			FreeUp_mem(Parms, varsum, anasum, vcv);
			ERRORPRT("The initial values have to be: 1 > background >= 0,  slope >0 (or 1 when there is restriction on slope >=1). ");
		}
	}
	/* end if (initial == yes), user initialized values */

	/* compute init_lkf for full model and init_lkr for reduced model
	*  The full model assumes that the number of responders at each dose
	*  is independently binomially distributed with individual probability
	*  denoted by p[i] with MLE
	*  W = (number of respondents at dose i)/(total subjects at dose i).
	*  The reduced model assumes that the respondents all come from a
	*  single binomial distribution with probability p and MLE
	*  W = (total number of respondents) / (total subjects)
	*  lkf and lkr are the log-likelihood functions for the binomial dist*/
	lkf = 0.0;
	varsum[1].S = 0;
	varsum[2].S = 0;
	for (i=1;i<=Nobs;i++)
	{
		/*Likelihood for full model*/
		varsum[1].S += Yp[i];		/* total number of responders */
		varsum[2].S += Yn[i];		/* total number of nonresponders */
		W = Yp[i] / (Yp[i]+Yn[i]);	/* proportion of positive responders */
		if (W > 0)
		{
			lkf += Yp[i] * log(W);
		}

		if (W < 1)
		{
			lkf += Yn[i] * log(1- W);
		}
	}
	Maxloglik = -lkf;
	/* end for, likelihood for full model */

	/* Likelihood for reduced model */
	W = varsum[1].S / (varsum[1].S + varsum[2].S);
	lkr = varsum[1].S * log(W) + varsum[2].S * log(1- W);

	/* fitting Logist model and output parameter estimators */
	Logist_fit(nparm, Parms, EPS, &iter, &xlk);
	bounded = IVECTOR(1, nparm);
	Which_Bounded (Spec, Parms, bounded);
	if (ErrorFlag != 0)
	{
		bmdose = No;
		for (i = 1; i <= nparm; i++) bounded[i] = 1;
	}

	if (ErrorFlag == 0)
	{
		/* Compute the approx. covariance matrix */
		/* Compute the Hessian (in Logist_vcv) for ALL parameters, */
		/* then drop out rows and columns corresponding to fixed/bounded */
		/* parameters, invert what's left and print it out. */
		/* Get_and_OUTPUT_DTMSVCV does all that*/
		INITIALIZE_DMATRIX(vcv, nparm, nparm);
		Logist_vcv(nparm,Spec,Parms,vcv);


		/* Find the parameters on their boundaries */

		adj_vcv_rows = 0;
		for (i=1; i<=nparm; i++)
		{
			if (bounded[i] == 0)
			{
				adj_vcv_rows ++;
			} /* end if */
		} /* end for */

		if (adj_vcv_rows > 0)
			vcv_adj = DMATRIX(1, adj_vcv_rows, 1, adj_vcv_rows);

		/* Output covariance matrix */
		if (adj_vcv_rows > 0)
		{
			Get_and_OUTPUT_DTMSVCV(nparm, Spec, Parm_name, vcv, vcv_adj, bounded);
		}
	}
	/* Output parameter estimates */
	if (ErrorFlag != -1)
		OP_ParmsE(nparm,Spec,Parms,Parm_name,vcv_adj, bounded, bmdparm.level, 1); 
	if (adj_vcv_rows > 0) {
		FREE_DMATRIX(vcv_adj, 1, adj_vcv_rows, 1, adj_vcv_rows);
	}


	/* Compute and output ANOVA table elements */
	DTMS3ANOVA (nparm,Nobs,Spec,lkf,xlk,lkr,anasum, bounded);


	/* Output ANOVA table */
	OUTPUT_DTMS3ANOVA(anatxt,anasum);
	fflush(fp_out);
	if (ErrorFlag == 0)
	{
		/* Print a goodness of fit table */
		Quantal_Goodness(nparm, bounded, Parms, Nobs, Xi, Yp, Yn, scale);

		/* Compute benchmark dose */
	}
	/* Setup plotting file (.002 extension) */

#ifndef RBMDS
	fprintf (fp_out2, "\n BMD_flag \t %d \n Nosb \t%d \n nparm \t%d",  bmdose, Nobs, nparm );
	fprintf (fp_out2, "\n  Con_lev \t%3.3g ", bmdparm.level);
	fprintf (fp_out2, "\n  logtrans \t%d ", logtrans);
	fprintf (fp_out2, "\n  RiskType \t%d ", bmdparm.risk);
	fprintf (fp_out2, "\n  Effect \t%3.3g ", bmdparm.effect);
	for (i=1;i<=nparm; i++)
	{
		fprintf (fp_out2, "\n %s \t %5.3g", Parm_name[i-1], Parms[i]);
	}

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
#endif

	if (bmdose==Yes)
	{
		if (logtrans==Yes)
		{
			back=Parms[(int) Background];
		}
		else
		{
			back=1/(1+exp(-Parms[(int) Intercept])) ;
		}

		if (bmdparm.risk==1)
		{
			back1=1;
		}
		else
		{
			back1=(1-bmdparm.effect);
		}

		/*Values used to calculate BMD and possibly BMDL values on the curve
		if specified by user*/
		Logist_BMD (nparm, Parms, EPS, &junk, xlk, Rlevel, Bmdl, &BMD);

		if (BMD >= 0.0)
		{
			OUTPUT_BENCHMD(1,(BMD)*scale);
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
			/*  	      printf("           BMDL =%14.6g\n\n", (Bmdl[1]));   */
			fprintf(fp_out, 
#ifndef RBMDS
				"            BMDL = %14.6g\n\n"
#else
				"            BMDL = %30.22g\n\n"
#endif

				, Bmdl[1]);
			/*      fprintf(fp_out, "Cancer Slope Factor =%14.6g\n\n", bmdparm.effect/Bmdl[1]); */

#ifndef RBMDS
			fprintf (fp_out2, "\n  RSL \t%f",bmdparm.effect+back1*back);
			fprintf (fp_out2, "\n  BMD \t%f",BMD);
			fprintf (fp_out2, "\n  BMDL \t%f",Bmdl[1]);

			fprintf (fp_out2,"\n\n BMD_line");
			fprintf (fp_out2,"\n %f %f", (xmin-xmax/100), bmdparm.effect+back1*back );
			fprintf (fp_out2,"\n %f %f", BMD, bmdparm.effect+back1*back );
			fprintf (fp_out2,"\n %f %f", BMD, -0.1);

			fprintf (fp_out2,"\n\n BMDL_line");
			fprintf (fp_out2,"\n %f %f", Bmdl[1], -0.1);
			fprintf (fp_out2,"\n %f %f", Bmdl[1], bmdparm.effect*back1+back );

			fprintf (fp_out2, "\n\n BMDL_Curve_flag \t %d  \n smooth_opt  %d", bmdlCurve, smooth);
			fprintf (fp_out2,"\n\n BMDL_curve");
			fprintf (fp_out2,"\n 0.00000 %f", back);

			for (i=1;i<=5;i++)
			{
				if (bmdparm.risk==1)
				{
					back1=1;
				}
				else
				{
					back1=(1-Rlevel[i]);
				}

				fprintf (fp_out2,"\n %f %f", Bmdl[i], Rlevel[i]+back1*back);
			}
#endif
		}
	}


	/*indicate all required data have been output.*/
#ifndef RBMDS
	fprintf (fp_out2,"\n\n Check_result %d", Yes);
#endif
	/* free allocated memory */
	FREE_DVECTOR (Parms,1,nparm);
	FREE_DVECTOR (Ypp,1,Nobs);
	FREE_DVECTOR (Ep,1,Nobs);
	FREE_DVECTOR (Xi,1,Nobs);
	FREE_DVECTOR (Yp,1,Nobs);
	FREE_DVECTOR (Yn,1,Nobs);
	FREE_DVECTOR (IniP,1,nparm);
	FREE_IVECTOR (Spec,1,nparm);
	FREE_IVECTOR (IniSp,1,nparm);
	FREE_IVECTOR (bounded, 1, nparm);
	FREE_VLVECTOR(varsum,1,3);
	FREE_ALVECTOR(anasum,1,3);
	FREE_DVECTOR (Rlevel,1,5);
	FREE_DVECTOR (Bmdl, 1, 5);
	FREE_DMATRIX(vcv,1,nparm,1,nparm);
	CLOSE_FILES ();
	return(0);

}	/*end of main*/

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
	if (logtrans == No)
	{
		for (i = 1; i <= ndoses; i++)
			P[i] = 1.0 / (1.0 + exp(-Parms[(int) Intercept] - Parms[(int) Slope] * doses[i]));
	}
	else
	{
		for (i = 1; i <= ndoses; i++)
			if (doses[i] <= 0.0)
				P[i] = Parms[(int) Background];
			else
				P[i] = Parms[(int) Background] + (1.0 - Parms[(int) Background])/
				(1.0 + exp(-Parms[(int) Intercept] - Parms[(int) Slope] * log(doses[i])));
	}
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
	/* if replace == Yes, then we are computing BMD; if logistic, */
	/* replace p[(int) Slope], and if log-logistic, replace p[(int) Intercept] */
	/* with function of current guess at BMD (i.e., tD) */
	if (replace==Yes)
	{
		if(logtrans==Yes) /* replace p[(int) Intercept]; tD is already on the log scale */
		{
			if(bmdparm.risk==ADDED)
			{
				p[(int) Intercept] = -log((1.0 - p[(int) Background])/BMR - 1.0) - p[(int) Slope] * tD ;
			}
			else
			{
				p[(int) Intercept] = -log((1.0 - BMR)/BMR) - p[(int) Slope] * tD;
			}
		}
		else /* replace p[(int) Slope]; not log dose */
		{
			double Z;
			if (bmdparm.risk == ADDED)
			{
				Z = BMR * (1.0 + exp(-p[(int) Intercept]))/exp(-p[(int) Intercept]);
			}
			else
			{
				Z = BMR;
			}
			p[(int) Slope] = -log((1.0 - Z)/(1.0 + Z * exp(-p[(int) Intercept])))/tD;
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

	if (logtrans == No)
	{
		bounded[1] = 1;
		if (!fixedParm(Slope) && VERYCLOSE(Parms[(int) Slope],(log(Max_double)/xmax)))
			bounded[3] = 1;
	}
	else /* logtrans == Yes */
	{
		if (VERYCLOSE(Parms[(int) Background],0.0) || VERYCLOSE(Parms[(int) Background],1.0)) bounded[1] = 1;
		if ((restrict == Yes && VERYCLOSE(Parms[(int) Slope],1.0)) ||
			(restrict != Yes && VERYCLOSE(Parms[(int) Slope],0.0)) ||
			VERYCLOSE(Parms[(int) Slope], SlopeUpperBound)) bounded[3] = 1;
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
**Logist_lk -- used to compute the (minus) log likelihood for Logist model.
* 		     Extern var.: smean, smax, Nobs, Xi, Yp, Yn, Ls, Xg.
*
*		input:
*			nvar is the number of parameters
*			nf is something for the fortran maximization--not used
*			uiparm is a vector of the indices for the parameters
*				of length nparm
*			urparm is a vector of fixed parameters
*				of length jfixed
*			ufparm is something for the fortran maximization--not used
*		output:
*			f is the likelihood value
*		input/output:
*			x is a vector of non-fixed parameters
*				of length jvar
*
*********************************************************************/

void Logist_lk(long int *nvar, double *x, long int *nf, double *f,
			   long int *uiparm, double *urparm, void (*ufparm)())
{

	int     i;
	double   xlk;          			/*log likelihood;*/
	double   *p;
	double   *Pred;

	/*parms for calculation*/

	p=DVECTOR(1,nparm);
	Pred = DVECTOR(1, Nobs);
	/* construct the parameter vector */
	unpack(x, urparm, p);

	/* Compute the likelihood */
	xlk = 0.0;
	Predict(Xi, Nobs, p, Pred);
	for (i=1;i<=Nobs;i++)
	{
		xlk += Yp[i] * Slog(Pred[i]) + Yn[i] * Slog(1.0 - Pred[i]);
	}
	FREE_DVECTOR(p, 1, nparm);
	FREE_DVECTOR(Pred, 1, Nobs);
	*f= -xlk;

}	/*end of Logist_lk*/

/*******************************************************************
**Logist_g -- used to compute the gradients for Logist model.
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

void Logist_g (long int *nvar, double *x, long int *nf, double *g,
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
		if (fabs(x[i]) > 1.0e-7)
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
		Logist_lk(nvar, saveparms, &inf, &basefun, uiparm, urparm, ufparm);
		saveparms[i] = x[i] + h[i];
		Logist_lk(nvar, saveparms, &inf, &tmp, uiparm, urparm, ufparm);
		g[i] = (tmp - basefun) / (2.0*h[i]);
	}

	FREE_DVECTOR(saveparms,0,nparm);
	FREE_DVECTOR(h, 0, nparm);
}


/****************************************************************************
 **
 * Logist_grad -- Computes the gradient of the logistic likelihood
 * function with respect to the user form of the parameters.  This is to
 * be used in Logist_vcv, to compute a finite difference approximation to
 * the hessian of the likelihood function
 * Input: nparm -- the number of parameters, the dimension of all the following
 arrays.
 Spec[] -- if the ith parameter is fixed by the user, Spec[i] == 1,
 otherwise Spec[i] == 0.
 ptf[] -- vector of parameters, in user form (external form),
 based at 1.
 * Output: grad[] -- the gradient of the loglikelihood function
 (N.B. Logist_g returns the gradient of -loglikelihood)
 Based at 1.

 *****************************************************************************/
void Logist_grad(int nparm, int Spec[], double ptf[],
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
	Logist_g (&nvar, start, &nf, outgrad, uiparm, urparm, ufparm);
	/* Logist_g returns the gradient of -loglikelihood */
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
*	Logist_vcv -- used to compute the vcv matrix for logistic model.
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

void Logist_vcv (int nparm, int Spec[], double ptf[], double **vcv)
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
			Logist_grad(nparm, Spec, saveparms, gradp);
			saveparms[i] = ptf[i] - h[i];
			Logist_grad(nparm, Spec, saveparms, gradm);
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
} /* end of Logist_vcv */






/**************************************************************
 *MAX_lk -- used to obtain the Maximum log-lilikehood as well as
 *           the estimates of parameters, given initial p[1..nparm],
 *           object func. , and gradient func. G_func.
 *
 *		input:
 *			nparm is the number of parameters in the model
 *			gtol is the tolerance level
 *			iter is the number of iterations
 *		output:
 *			fret is the value of the log-likelihood function
 *		input/output:
 *			p[] is the parameter vector of length nparm
 *			    upon output, has the maximum likelihood estimates
 *
 **************************************************************/

void MAX_lk(int nparm, double p[], double gtol, int *iter, double *fret)

{
	int i, jfixed, jvar;
	long int nvar;
	long int *uiparm;

	double *start;
	double *urparm;
	double lower[10],upper[10]; /* hard coded the max number of parameters to 3 here */
	void (*ufparm)();

	/* Set up initial parameter array, start.  All parameters go either
	into start (if they are changing to improve the fit) or urparm
	(if they are fixed).  Set up all elements of Dscale to be 1.0.*/
	nvar = nparm;
	for (i = 1; i <= nparm; i++) nvar = nvar - Spec[i]; /*Count the varying
														parameters */
	if (nvar > 0)
	{
		urparm = (double *) malloc((size_t)(nparm-nvar)*sizeof(double));
		start = (double *) malloc((size_t) nvar*sizeof(double));
	}
	else
	{
		long int dummy;
		urparm = (double *) malloc((size_t)(nparm-nvar)*sizeof(double));
		for(i=1;i<=nparm;i++) urparm[i-1]=p[i];

		Logist_lk( &nvar, start,  &dummy, fret,
			uiparm, urparm,  ufparm);
		ErrorFlag = -1;
		*fret= -(*fret);
		return;
	}

	jfixed = 0;
	jvar = 0;
	for (i = 1; i <= nparm; i++) /*separate the fixed and variable
								 parameters*/
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
		}
	}

	/* set up the bounds variable.  Each parameter is unique, here. */
	jvar = 0;
	if (!fixedParm(Background))
	{
		lower[jvar] = 0.0;
		upper[jvar] = 1.0;
		jvar++;
	}

	if (!fixedParm(Intercept))
	{
		lower[jvar] = -Max_double;
		upper[jvar] = Max_double;
		jvar++;
	}

	if (!fixedParm(Slope))
	{
		if (restrict == 1)
		{
			lower[jvar] = 1.0;
		}
		else
		{
			lower[jvar] = 0.0;
		}
		if (logtrans == Yes)
			upper[jvar] = SlopeUpperBound;
		else
			upper[jvar] = log(Max_double);
	}

	ErrorFlag = run_dmngb((int) nvar, start, lower, upper, Maxloglik,
		Rel_Conv, Parm_Conv, ITMAX, 10,
		Logist_lk,Logist_g,uiparm,urparm,ufparm,
		0,fret);
	*fret = -*fret;

	/* Put the parameter values back in p */
	unpack(start, urparm, p);

	/* free allocated memory */
	free(urparm);
	free(start);

}	/*end of MAX_lk*/


/**************************************************************
*Logist_fit -- Used to "prepare" the data for further computation,
*            i.e. compute the extern variables, give the initial
*            parameters, etc. THEN fit the Logist model.
*            (In fact, these jobs could be done in main().)
*
*		input:
*			nparm is the number of parameters in the model
*			gtol is the tolerance level
*			iter is the number of iterations
*		output:
*			fret is the value of the log-likelihood function
*		input/output:
*			p[] is the parameter vector of length nparm
*				upon output, has the maximum likelihood estimates
*
***************************************************************/
void Logist_fit(int nparm, double p[], double gtol,
				int *iter, double *fret)
{
	int    i, j, junk;
	double  ymin, W, dos, xlk, junk1,junk2, range;
	double *pBak, *tmy, *t, **tmv,**vcv;

	/*Create vectors and matrices*/
	pBak=DVECTOR(1, nparm);
	tmy=DVECTOR(1, nparm);
	t=DVECTOR(1, nparm);
	tmv=DMATRIX(1,nparm,1,nparm);
	vcv=DMATRIX(1,2,1,2);

	ymin = 1.0;
	xmin = 1000000;
	xmax = 0.0;
	for (i=1;i<=Nobs;i++)
	{
		dos=Xi[i];
		if (dos < xmin)
		{
			xmin = dos;
		}

		if (dos > xmax)
		{
			xmax = dos;
		}
	}
	/* For normal logistic model, scale dose. */
	/* Remember to rescale when we're done. */
	if (logtrans != Yes)
	{
		scale=xmax;

		for (i = 1; i <= Nobs; i++)
			Xi[i] = Xi[i] / scale;
	}
	else
		scale = 1.0;

	for (j=1; j<=nparm; j++)
	{
		pBak[j]=p[j];			/*save the input p[].*/
	}

	/****** Obtain initial estimates for p[] ******/
	if(initial==Yes)
	{
		for(j=1; j<=nparm; j++)
		{
			p[j]=IniP[j];
		}
	}
	else
	{
		/*compute initial estimates*/
		if (logtrans == Yes)
		{
			int contdose;
			p[(int) Background] = 0.0;
			contdose = 0;
			for (i = 1; i <= Nobs; i++)
			{
				if (Xi[i] == 0)
				{
					contdose = i;
					p[(int) Background] = Yp[i]/(Yp[i] + Yn[i]);
				}
			}
			INITIALIZE_DMATRIX(vcv,2,2);
			junk1=junk2=0.0;
			for (i=1;i<=Nobs;i++)
			{
				if (i == contdose) continue;
				W = Yp[i] / (Yp[i]+Yn[i]);
				W = (W - p[(int) Background])/(1 - p[(int) Background]);
				if (W <= 0) W = 0.5/(Yp[i]+Yn[i]+1);
				if (W >= 1) W = 1.0 - (0.5/(Yp[i]+Yn[i]+1));
				W = log(W / ( 1.0 - W));

				dos=log(Xi[i]);
				vcv[1][1] += 1.0;
				vcv[1][2] += dos;
				vcv[2][2] += dos*dos;
				junk1 += W;
				junk2 += W*dos;
			}

			vcv[2][1] = vcv[1][2];

			INVMAT(vcv, 2);

			p[(int) Intercept]=vcv[1][1]*junk1+vcv[1][2]*junk2;
			p[(int) Slope]=vcv[2][1]*junk1+vcv[2][2]*junk2;

			if (restrict==Yes && p[(int) Slope]<=1)
			{
				range=(xmax)/2;
				range=log(range);

				p[(int) Intercept]=p[(int) Intercept]+(p[(int) Slope]-1.0)*range;
				p[(int) Slope]=1.0;
			}
		}
		else /* logtrans == No */
		{
			INITIALIZE_DMATRIX(vcv,2,2);
			junk1=junk2=0.0;
			for (i=1;i<=Nobs;i++)
			{
				W = (Yp[i]+0.5) / (Yp[i]+Yn[i]+1);

				W = -1.0*log(1.0/W - 1.0);

				dos=Xi[i];
				vcv[1][1] += 1.0;
				vcv[1][2] += dos;
				vcv[2][2] += dos*dos;
				junk1 += W;
				junk2 += W*dos;
			}

			vcv[2][1] = vcv[1][2];

			INVMAT(vcv, 2);

			p[(int) Background]=0;
			p[(int) Intercept]=vcv[1][1]*junk1+vcv[1][2]*junk2;
			p[(int) Slope]=vcv[2][1]*junk1+vcv[2][2]*junk2;

		} /* end if (logtrans == Yes ){} else{} */
		/* Have to temporarily retransform p[Slope] */
		/**get specified parameters.**/
		for (j=1;j<=nparm;j++)
		{
			if (Spec[j]==Yes)
			{
				p[j]=pBak[j];
			}
		}
		if (!fixedParm(Slope))
			p[(int) Slope] /= scale;
		OUTPUT_TEXT("\n\n                  Default Initial Parameter Values  ");
		OUTPUT_Init(nparm, Spec, p, Parm_name);
		p[(int) Slope] *= scale;
	}



	/***       Fit the model.          ********************/
	replace = No;

	/* I deleted some code that repeated the fit 5 times, */
	/* perturbing the final point */
	/* a little bit each time (run_dmngb now does this internally)*/
	MAX_lk( nparm, p, gtol,  &junk,  &xlk);

	do_dmngb_warning(&ErrorFlag);

	*fret = xlk;
	/* Rescale doses and slope */
	if (logtrans != Yes)
	{
		for (i = 1; i <= Nobs; i++) Xi[i] *= scale;
		p[(int) Slope] /= scale;
		scale = 1.0;
	}

	/* free allocated memory */
	FREE_DVECTOR(pBak, 1, nparm);
	FREE_DVECTOR(tmy, 1, nparm);
	FREE_DVECTOR(t, 1, nparm);
	FREE_DMATRIX(vcv,1,2,1,2);
	FREE_DMATRIX(tmv,1,nparm,1,nparm);

}	/*end Logist_fit*/

/************************************************************
 * Logist_BMD -- Used to calculate the BMD and BMDL for Logist model.
 *
 *		Input:
 *			nparm is the number of parameters
 *			p[] is the parameter vector of length nparm
 *			gtol is the tolerance level
 *			iter is the number of iterations
 *
 *		Output:
 *			xlk is the likelihood value of the function
 *			Rlevel[] is the risk level
 *			Bmdl[] is the benchmark dose confidence limit
 *			BMD is the benchmark dose
 *
 *		Input/Output:
 *
 *************************************************************/
void Logist_BMD (int nparm, double p[], double gtol, int *iter, double xlk,
				 double Rlevel[], double Bmdl[], double *BMD)
{

	double   tol;
	double   xa,xb,fa,fb;
	double   Bmdjunk, BMDtmp;
	double   *pBak, *pa, *pb;
	double   stepsize;
	int      i, j, k, bogusBMD;
	double   *pred, *doses, effect;


	pBak=DVECTOR(1, nparm);
	pa = DVECTOR(1, nparm);
	pb = DVECTOR(1, nparm);
	stepsize = 0.9;
	for(j=1; j<=nparm; j++)
	{
		pBak[j]= p[j];          /*save the p[].*/
	}

	/**** compute X^2 value  ************************/
	/* IF ML is the value of the maximized log-likelihood, then ML - LR is the value
	log-likelihood at the BMDL or BMDU */
	if (bmdparm.level<0.5)
	{
		LR = 0.5*QCHISQ(1.0 - 2.0 * bmdparm.level, 1);
	}
	else
	{
		LR = 0.5*QCHISQ(2.0 * bmdparm.level - 1.0, 1);
	}

	Rlevel[1] = BMR = bmdparm.effect;

	/**** solve the BMD ********************************/

	/* If there is no slope, we want to cap the BMD estimate. */
	/* We still want to compute a lower limit on the BMD, however */
	bogusBMD = 0;
	if (p[(int) Slope] <= DBL_MIN)
	{
		*BMD = 100 * xmax;
		bogusBMD = 1;
		Warning(" Slope parameter essentially zero.  BMD set to 100 * max(Dose).");
		if (logtrans == Yes)
		{
			xb = log(*BMD);
		}
		else
		{
			xb = *BMD;
		}
	}
	else /* there is a substantial (though not necessarily significant) slope */
		if(logtrans==Yes)
		{
			if (bmdparm.risk == ADDED)
			{
				xb = -(p[(int) Intercept] + log((1.0 - p[(int) Background] - BMR)/BMR))/p[(int) Slope];
			}
			else
			{
				xb = -(p[(int) Intercept] + log((1.0 - BMR)/BMR))/p[(int) Slope];
			}
			*BMD = exp(xb);
		}
		else
		{
			double Z;
			if (bmdparm.risk == ADDED)
			{
				Z = BMR * (1.0 + exp(-p[(int) Intercept]))/exp(-p[(int) Intercept]);
			}
			else
			{
				Z = BMR;
			}
			*BMD = -log((1.0 - Z)/(1.0 + Z * exp(-p[(int) Intercept])))/p[(int) Slope];
			xb = *BMD;
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
		FREE_DVECTOR(pred,1,4);
		FREE_DVECTOR(doses,1,4);
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
		/* OK -- We have a legitimate BMD */

		/********* search for BMDL **************************/

		/* Replace either the slope (ordinary logistic) or the intercept
		(log-logistic) with an expression involving the BMD */

		replace=Yes;

		if (logtrans == Yes)
			Spec[(int) Intercept] = 1;
		else
			Spec[(int) Slope] = 1;

		/* First get the interval to search over */
		if (logtrans == Yes)
		{
			xa = xb + log(stepsize);
		}
		else
		{
			xa = xb*stepsize;
		}
		tol = FMAX(xb*0.000001, 0.0000001);

		BMD_lk = xlk;		/*get the lk at BMD.*/
		fb = DBL_MAX;
		for (i = 1; i <= nparm; i++) pa[i] = pb[i] = p[i];
		fa = BMDL_func(nparm, p, xa, tol);

		while (fa<0.0 &&
			((logtrans == No && xa > DBL_MIN) || (logtrans == Yes && exp(xa) > DBL_MIN)) && fabs(fb - fa) > DBL_EPSILON)
			/* BMDL does not lie between xa and xb, so reduce xa */
			/* sneak down from BMD, replacing xb with the previous xa, */
			/* fb with the previous fa */
			/* modify BMDL_func so that it returns in p the current parameter estimates */
			/* so that we are never very far from the (constrained) MLEs */
		{
			xb = xa;
			fb = fa;
			for (i = 1; i <= nparm; i++) pb[i] = pa[i];
			if (logtrans == Yes)
			{
				xa += log(stepsize);
			}
			else
			{
				xa *= stepsize;   /*prevent that xa=0.1*BMD is not small enough.*/
			}
			fa = BMDL_func(nparm, pa, xa, tol);
		}
		if (fa<0.0)
		{
			/*computation failed*/
			Bmdl[1] = -1.0;
		}
		else
		{
			/*computation will succeed.*/
			Bmdl[1] = zeroin(xa, xb, 1.0e-10, BMDL_func, nparm, pb, 1.0e-14);
			BMDL_Error_Size = BMDL_func(nparm, pb, Bmdl[1], tol);
			if(logtrans==Yes)
			{
				Bmdl[1]=exp(Bmdl[1]);
			}
			else
			{
				Bmdl[1] = Bmdl[1] * scale;
			}
		}


		if(bmdlCurve==Yes)
		{
			/****** calculate Bmd[] and Bmdl[] ***********/
			for (k=2; k<=5; k++)
			{
				/* 1: reload the original parameter vector */
				for(j=1; j<=nparm; j++)
				{
					p[j]= pBak[j];
				}
				/* 2:  Set up the BMR levels to use (may want to change this code) */
				if (k==2)
				{
					Rlevel[k]=BMR=0.05;
				}
				else
				{
					Rlevel[k]= BMR = (k-2)*0.1;
				}


				/* 3: solve the BMD[] */
				/* If there is no slope, we want to cap the BMD estimate. */
				/* We still want to compute a lower limit on the BMD, however */
				if (p[(int) Slope] <= DBL_MIN)
				{
					BMDtmp = 100 * xmax;
					if (logtrans == Yes)
					{
						xb = log(BMDtmp);
					}
					else
					{
						xb = BMDtmp;
					}
				}
				else /* there is a substantial (though not necessarily significant) slope */
					if(logtrans==Yes)
					{
						if (bmdparm.risk == ADDED)
						{
							xb = -(p[(int) Intercept] + log((1.0 - p[(int) Background] - BMR)/BMR))/p[(int) Slope];
						}
						else
						{
							xb = -(p[(int) Intercept] + log((1.0 - BMR)/BMR))/p[(int) Slope];
						}
						BMDtmp = exp(xb);
					}
					else
					{
						double Z;
						if (bmdparm.risk == ADDED)
						{
							Z = BMR * (1.0 + exp(-p[(int) Intercept]))/exp(-p[(int) Intercept]);
						}
						else
						{
							Z = BMR;
						}
						BMDtmp = -log((1.0 - Z)/(1.0 + Z * exp(-p[(int) Intercept])))/p[(int) Slope];
						xb = BMDtmp;
					}

					Bmdjunk = xb;


					/********* search for BMDL[] **************************/

					/* First get the interval to search over */
					if (logtrans == Yes)
					{
						xa = xb + log(stepsize);
					}
					else
					{
						xa = xb*stepsize;
					}
					tol = FMAX(xb*0.0001, 0.0000001);

					BMD_lk = xlk;		/*get the lk at BMD.*/
					fb = DBL_MAX;
					for (i = 1; i <= nparm; i++) pa[i] = pb[i] = p[i];
					fa = BMDL_func(nparm, pa, xa, tol);

					while (fa<0.0 &&
						((logtrans == No && xa > DBL_MIN) || (logtrans == Yes && exp(xa) > DBL_MIN)) && fabs(fb - fa) > DBL_EPSILON)
						/* BMDL does not lie between xa and xb, so reduce xa */
					{
						xb = xa;
						fb = fa;
						for (i = 1; i <= nparm; i++) pb[i] = pa[i];
						if (logtrans == Yes)
						{
							xa += log(stepsize);
						}
						else
						{
							xa *= stepsize;   /*prevent that xa=stepsize*BMD is not small enough.*/
						}
						fa = BMDL_func(nparm, pa, xa, tol);
					}
					if (fa<0.0)
					{
						/*computation failed*/
						Bmdl[k]=-1;
						fprintf(fp_out, "\n BMDL curve computation failed for BMR = %f . \n The BMDL curve appearing in the graph may not be accurate.", BMR);
					}
					else
					{
						/*computation will succeed.*/
						Bmdl[k] = zeroin(xa, xb, 1.0e-10, BMDL_func, nparm, pb, 1.0e-14);

						if(logtrans==Yes)
						{
							Bmdl[k]=exp(Bmdl[k]);
						}
						else
							Bmdl[k] = Bmdl[k] * scale;
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

		FREE_DVECTOR(pBak, 1, nparm);
		FREE_DVECTOR(pa, 1, nparm);
		FREE_DVECTOR(pb, 1, nparm);

}  	/*end of Logist_BMD*/

/*****************************************************************
* BMDL_func -- used to compute the values of functions BMDL_f (the
*              X^2 value) at the point D, given the parm p[] and
*              number of parm. It returns in p the (constrained)
*              MLE for p.
*
*              This routine is called by zeroin().
*
*****************************************************************/
double BMDL_func(int nparm, double p[], double D, double gtol)
{
	/* BMD_lk and LR are calculated in Probit_BMD() */
	/* WARNING: p is modified by this function! */

	double fD, xlk;

	int junk;

	tD = D;   /*tD is global var. have to change before call MAX_lk()*/

	MAX_lk( nparm, p, gtol,  &junk,  &xlk);
	fD = BMD_lk - xlk - LR;

	return fD;

}	/*end of BMDL_func*/

