/****************************************************************
**
* IMPORTANT NOTE:  The following variable is the version number for
*                  the current model.  THIS MUST BE CHANGED as
*				   important changes are made to the models.
*
*****************************************************************/
char Version_no[]="Multistage Model. (Version: 3.4;  Date: 05/02/2014)";
/****************************************************************
char Version_no[]="Multistage Model. (Version: 3.3;  Date: 02/28/2013)";
char Version_no[]="Multistage Model. (Version: 3.2;  Date: 05/26/2010)";
char Version_no[]="Multistage Model. (Version: 3.1;  Date: 10/28/2009)";
char Version_no[]="Multistage Model. (Version: 3.0;  Date: 05/16/2008)";
char Version_no[]="Multistage Model. (Version: 2.9;  Date: 03/12/2008)";
char Version_no[]="Multistage Model. (Version: 2.8;  Date: 02/20/2007)";
char Version_no[]="Multistage Model. (Version: 2.6;  Date: 9/13/2006)";
*****************************************************************/

/****************************************************************
**
* Multistage.C - a ANSI C program for Multistage model fitting with/without
*             a natural background rate in Benchmark Dose.
*
* Date: Aug 21, 2000
*
********************************************************************
* Modification Log:
*
* Version Number: 2.3
* Modified By: Qun He
* Modified Date: 3/09/2003
* Reason:
*
* Version Number: 2.4
* Modified By: Micheal Ferree
* Modified Date: 8/13/2005
* Reason: Took out print out of cancer slope factor
*
* Version Number: 2.5
* Modified By: R. Woodrow Setzer
* Modified Date: 10/17/2005
* Reason: Free all allocated memory prior to exiting;
*         Fixed uses of the variable 'bind' from donlp2
*         Changed Chi^2 Res to Scaled Res for consistency with other quantal
*           models
*         Fixed scaled residuals so that
*           divisor is standard error, not variance.
*
* Version Number: 2.6
* Modified By: R. Woodrow Setzer
* Modified Date: 9/07/2006
* Reason: Now based on cancer model 1.3; computes upper bound;
*         printing of cancer slope factor, and plotting of linear
*         extrapolation line, depends on option (currently set to
*         '0': don't do it).
*
* Version Number: 2.7
* Modified By: Geoffrey
* Date: 1/12/2007
* Reason: Incremented version number.
*		 Added last parameter "0" (don't print SE) in OP_ParmsE().
*
* Version Number: 2.8
* Modified By: Woodrow Setzer
* Date: 2/20/2007
* Reason: Incremented version number to reflect changed compilation options.
*
* Version Number: Version: 2.9
* Modified By: G. Nonato
* Modification Date: 03/12/2008
* Reason: (Per BMDS 2.0: Problem Report 157 & 147)
*       Fix the Observation # < parameter # for Weibull model problem.
*       Added code to free-up allocated memories before exiting thru ERRORPRT()
*
* Version Number: 3.0
* Modified By: G. Nonato
* Modification Date: 05/16/2008
* Reason: (Per BMDS 2.0: Problem Report 165)
*       Goodness of Fit - Observed column, print values as real numbers.
*
* Version Number: 3.1
* Modified By: G. Nonato
* Modification Date: 10/28/2009
* Reason:
*      To be able to process files/folders with spaces (PR 257)
*      Fix program freeze due to long variable names (PR 278)
*      Process long files up to 256 characters (PR 303 and 308)
*      Modify code for easy maintenance, all lengths for file name,
*        model name, and column names are placed in benchmark.h
*
* Version Number: 3.2
* Modified By: G. Nonato
* Modification Date: 05/26/2010
* Reason: PR 319
*      Change ERRORPRT("Observation # < parameter # for Multistage Cancer model.");
*      to Warning("Observation # < parameter # for Multistage Cancer model.");
*
* Version Number: 3.3
* Modified By: Louis Olszyk
* Modification Date: 02/28/2013
* Reason: PR 444 - Fix wording in plot titles
*
*
* Version Number: 3.4
* Modified By: Louis Olszyk
* Modification Date: 05/02/2014
* Reason: 1. PRs 439,486 - Allow non-integer values for N and #affected
*         2. Merge in logic for "Cancer" model, which is identical,
*            except for printing out CSL. Code will be conditionally
*            compiled to produce one model vs. the other.
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


extern void loadcommbloc_(long int *ndoses, double *xmax,
			  double nanimals[],double doses[],
			  double affect[], long int *polyord,
			  long int *restrict, double *, double *);

/*  extern void ProfLik(int nparm, double parms[], double Max_LL, double BMD, */
/*  		    double BMDL, char *fname, double alpha, double BMR2, */
/*  		    int *ProfFlag, int *PEFlag, double *PE_bmdl); */


#define EPS 3.0e-8
#define float double
#ifdef BMDS_CANCER_MODEL
#  define BMDS_ENABLE_CSL 1
#else
#  define BMDS_ENABLE_CSL 0
#endif

void Multistage_fit(int nparm, double p[], double gtol,
		int *iter, double *fret);
void Multistage_BMD(int nparm, double p[], double gtol, int *iter, double xlk,
		double Rlevel[], double Bmdl[],double Bmdu[],double *BMD);
void Multistage_vcv(int nparm, int Spec[], double p[], double **vcv);
void GetNewParms(double *p, int size);
void GetMoreParms(double *p, int size);
void GetMLEParms(double *p, int size);
double P_RESPONSE(double q, double p[]);
int Model_DF (int []);

/*** Define input and output files's name  *********************/
char     fin[FLENGTH];  /* input temp file */
char    fout[FLENGTH];  /* output temp file */
char    fout2[FLENGTH]; /* output file for plot */
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


char plotfilename[FLENGTH];	/* plot file name */
char *anatxt[]={"Full model", "Fitted model", "Reduced model"};


/*** variables will not be changed except Spec  *******/
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
int nparm;         /* # of parameters */
long int restrict;      /* flag for restricting the betas */
/* restrict=1 if the betas >= 0 */
int initial;       /* flag for initialized parameters, initialized=1 */
int appendix;      /* flag for append or overwrite */
                   /* append = 1, overwrite = 0 */
int smooth;        /* flag for smooth option, 0=unique, 1=C-spline */
int	bmdlCurve;     /* flag for BMDL curve calculation, 1=yes */
int	bmduCurve;     /* flag for BMDU curve calculation, 1=yes */
double xmax, xmin; /* max and min dose levels */
double scale;      /* used to scale dose */
int *boundary;     /* same as bind array, =1 if a boundary was hit */

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
int PlotLinearExtrapolation, ReportCancerSlopeFactor;

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
	FREE_DVECTOR (Bmdu, 1, 5);
	FREE_DMATRIX(vcv,1,nparm,1,nparm);

	if (fp_log != (FILE *) NULL)
		fclose(fp_log);

	return;
}

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
  double   *Parms;           /* parameter array */
  VarList  *varsum;          /* info for variables--p. dep.,n. dep., indep. */
  AnaList  *anasum;          /* information for ANOVA analysis */
  double   **vcv;            /* variance and covariance matrix */
  char     model_name[MNLENGTH], user_note[UNLENGTH];
  char     dose_name[CNLENGTH], posi_name[CNLENGTH], nega_name[CNLENGTH], junkname[FLENGTH];
  int      *bounded;         /* bounded[i] = 1 if ith parm hits a bound  */
  double   **vcv_adj;        /* adjusted vcv matrix */
  int      adj_vcv_rows;     /* number of rows in adj. vcv matrix */
  double   *parameters;      /* small parameters are set to 0 */
  char long_path_name[FLENGTH];

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
  fscanf(fp_in, "%s", junkname);
  fscanf(fp_in, "%ld%ld",&Nobs,&ndegree);

  /* assign number of parameters */
  nparm = 1+ndegree;

  /* allocate memory for arrays */
  Parms = DVECTOR(1, nparm);
  IniP = DVECTOR(1, nparm);
  IniSp = IVECTOR(1, nparm);
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
  Bmdu = DVECTOR(1, 5);
  vcv = DMATRIX (1,nparm,1,nparm);

  /* read in user input info */
  fscanf(fp_in,"%d%lf%lf%d%ld%d%d%d", &ITMAX, &Rel_Conv, &Parm_Conv, &bmdlCurve, &restrict, &bmdose, &appendix, &smooth);
  fscanf(fp_in,"%lf%d%lf",&bmdparm.effect,&bmdparm.risk,&bmdparm.level);

  /* These flags control whether to report the cancer slope factor,
   * and whether to plot the linear extrapolation. For now, set them
   * here according to conditional compilation. Ideally, these could
   * be set at runtime by program name or some other switch.
   */
  ReportCancerSlopeFactor = BMDS_ENABLE_CSL;
  PlotLinearExtrapolation = BMDS_ENABLE_CSL;

  /* Set bmduCurve equal to bmdlCurve for now. MJF 07AUG05. */
  bmduCurve = bmdlCurve;

  junk = 0;



  /* get filenames */
  Get_Names(argv[1], fout, fout2, plotfilename);


  if(appendix==Yes)   /* if append was selected, then append output */
    fp_out=fopen(fout,"a");
  else                /* otherwise, overwrite the output */
    fp_out=fopen(fout,"w");

#ifndef RBMDS
  /* open the output file used for the plot */
  fp_out2=fopen(fout2,"w");
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
		FreeUp_mem(Parms, varsum, anasum, vcv);
		ERRORPRT("Error in opening output files.\n");
      //exit (1);
    } /* end if */

  /* Print model and file information on output page */
  Output_Header(Version_no, argv[1], plotfilename, ctime(&ltime), user_note);
  fflush(fp_out);

  if (bmdose < 0 || bmdose > 1)
  {
	FreeUp_mem(Parms, varsum, anasum, vcv);
    ERRORPRT("Error in choosing benchmark dose computation.");
  }




  /* obtain user input parameters, -9999 is the default value */
  READ_PARAMETERS(nparm,Parms);
  /* define vector of flags Spec[], 0=unknown, 1=specified */
  FILL_SPECVECTOR(nparm,Parms,Spec);
  /* count the number of specified parameters */
  nparm_known = COUNT_SPECVECTOR(nparm,Spec);
  brat=Yes;
  if (Spec[1]==1 && Parms[1]<EPS) /* if background is spec. and small */
    brat=No;

  /* determine if parameters are to be initialized */
  fscanf(fp_in,"%d", &initial);
  /* obtain the initial parameter values */
  READ_PARAMETERS(nparm,IniP);
  /* define vector of flags IniSp[], 0=initialized */
  FILL_SPECVECTOR(nparm,IniP,IniSp);
  /* define the specified initial parameter values to be 1 */
  for(i = 1; i <= nparm; i++)
    {
      if(Spec[i] == 1)
	IniP[i] = 1.0;
    } /* end for */

  /* read in the Dose name, Response name, and string NEGATIVE_RESPONSE */
  fscanf(fp_in,"%s%s%s", dose_name, posi_name, nega_name);
  /* determine the number of records with missing values */



  Nmiss = READ_OBSDATA3V(Nobs,3,2,3,1, Yp,Yn, Xi);
  /* Xi[]=dose levels, Yp[]= # of positive responses at the ith dose level */
  /* Yn[]= total minus Yp[] */

  Nobs -= Nmiss;       /* extern variable Nobs has been changed */
  /* Nobs is now the # of obs w/o missing values */

  /* print warning if there are too few parameters */
  if (Nobs < (nparm-nparm_known))
  {
	// Commented the following 2 lines of code for PR 319
	//FreeUp_mem(Parms, varsum, anasum, vcv);
	//ERRORPRT("Observation # < parameter # for Multistage model.");
    Warning("Observation # < parameter # for Multistage model.");
  }


  /*********** end of input data *************/

  /************** output title and summary of input data  ***********/

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
  if (Spec[1]==Yes)  /* if the background parameter was specified */
  {
	  fprintf(fp_out,"\n");
	  if (Parms[1] <= 0.0000001) /* if the background is small */
		  fprintf(fp_out,"\n   Background parameter is set to zero");
	  else
		  fprintf(fp_out,"\n   Background parameter is set to %g", Parms[1]);

	  if ((Parms[1] < 1.0e-20) && (Xi[1] < 1.0e-20) && (Yp[1] > 0))
	  {
		  FreeUp_mem(Parms, varsum, anasum, vcv);
		  ERRORPRT("ERROR: Background parameter specified as 0, but there were responses\n     at the control dose");
	  }
  } /* end if */

  if(nparm_known > 0+Spec[1])
    {
      fprintf (fp_out, "\n\n   User specifies the following parameters:");
      for (i=2; i<= nparm; i++)
	{
          if(Spec[i] == 1)
	    fprintf (fp_out, "\n %15s = %10.5g", Parm_name[i-1], Parms[i]);
	}  /* end for */
      fprintf (fp_out, "\n");
    } /* end if */
  /* determine the # of specified parameters */
  nparm_known = COUNT_SPECVECTOR(nparm, Spec);
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
	  OUTPUT_Init(nparm, Spec, IniP, Parm_name);

	  for (i=1; i<=nparm; i++)
	  {
		  if(IniSp[i]==1)    /* if the parameter was initialized */
		  {
			  if(Spec[i]==1 )  /* check to see if the parameter is fixed */
				  Warning("The initial value for the fixed parameter is ignored.");
		  } /* end if */
		  else
		  {
			  /* check if all the unspecified parms were initialized */
			  if (Spec[i]==0)
			  {
				  FreeUp_mem(Parms, varsum, anasum, vcv);
				  ERRORPRT("When the initial option is chosen, one has to initial ALL unspecified parameters.");
			  }
		  }  /* end else */
	  } /* end for */
	  i=0;
	  if(restrict==Yes)  /* if the restrict betas option was chosen */
	  {
		  for(j=2;j<=nparm;j++)
		  {
			  if(IniP[j]<0) i++;  /* if initial value is negative */
		  } /* end for */
	  } /* end if */
	  if (IniP[1] < 0 || IniP[1]>1 || i>0)
	  {
		  FreeUp_mem(Parms, varsum, anasum, vcv);
		  ERRORPRT("The initial values have to be: 1 > Bg >= 0 and Betas >= 0 ). ");
	  }
  }  /* end if (initial==Yes) */

  /* compute init_lkf for full model and init_lkr for reduced model */
  lkf = 0.0;
  varsum[1].S = 0;
  varsum[2].S = 0;
  for (i=1;i<=Nobs;i++)
    {
      varsum[1].S += Yp[i];
      varsum[2].S += Yn[i];
      W = Yp[i] / (Yp[i]+Yn[i]); /* proportion of response */
      if (W > 0)   lkf += Yp[i] * log(W);
      if (W < 1)   lkf += Yn[i] * log(1- W);
    } /* end for */
  W = varsum[1].S / (varsum[1].S + varsum[2].S);
  lkr = varsum[1].S * log(W) + varsum[2].S * log(1- W);

  /* fitting Multistage model and output parameter estimators */
  fflush(fp_out);
  if(nparm_known<nparm)
    Multistage_fit(nparm, Parms, EPS, &iter, &xlk);
  /* Parms[] is now the fitted parameters */
  /* xlk is the log-likelihood            */
  fflush(fp_out);
  /* initialize vcv (varaince-covariance matrix) */
  INITIALIZE_DMATRIX(vcv, nparm, nparm);
  /* compute the approximate variance-covariance matrix */
  /* variance-covariance matrix is the inverse of the info. matrix */
  Multistage_vcv(nparm,Spec,Parms,vcv);

  /* define bounded vector */
  bounded = IVECTOR(1, nparm);
  for (i=1; i<=nparm; i++)
    {
      if(nparm_known<nparm)
	bounded[i] = boundary[i];
      else
	bounded[i]=0;
    }
  if (Parms[1] == 0)
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
    Get_and_OUTPUT_DTMSVCV(nparm,Spec,Parm_name,vcv,vcv_adj,bounded);
  fflush(fp_out);
  parameters = DVECTOR(1, nparm);
  for (i=1; i<=nparm; i++)
    {
      if (bounded[i] == 1)
	parameters[i] = 0;
      else
	parameters[i] = Parms[i];
    }

  /* output model fitting parameters and standard errors */
  OP_ParmsE(nparm,Spec,parameters,Parm_name,vcv_adj, bounded, bmdparm.level, 1);
  fflush(fp_out);
  /* free memory */
  /*  FREE_IVECTOR(bounded, 1, nparm);   need bounded for DTMS3ANOVA */
  FREE_DVECTOR(parameters, 1, nparm);

  /* compute ANOVA table elements */
  DTMS3ANOVA (nparm,Nobs,Spec,lkf,xlk,lkr,anasum, bounded);
  fflush(fp_out);

  anasum[2].TEST = CHISQ(anasum[2].MSE, anasum[2].DF);
  fflush(fp_out);
  /* output ANOVA table */
  OUTPUT_DTMS3ANOVA(anatxt,anasum);
  fflush(fp_out);

  /* print a goodness of fit table */
  Quantal_Goodness(nparm, bounded, Parms, Nobs, Xi, Yp, Yn, scale);
  fflush(fp_out);

  /******************* compute benchmark dose ***********************/

#ifndef RBMDS
  /* print info to file used for plot */
  fprintf (fp_out2, "\n BMD_flag \t %d \n Nobs \t%ld \n nparm \t%d",  bmdose, Nobs, nparm );
  fprintf (fp_out2, "\n  Con_lev \t%3.3g ", bmdparm.level);
  fprintf (fp_out2, "\n  RiskType \t%d ", bmdparm.risk);
  fprintf (fp_out2, "\n  Effect \t%3.3g ", bmdparm.effect);
  for (i=1;i<=nparm; i++) fprintf (fp_out2, "\n %s \t %5.3g", Parm_name[i-1], Parms[i]);

  /* Calculate 95% CIs at each dose level for graphical output */
  {
    double *LL, *UL, *estp;

    LL = DVECTOR(1, Nobs);
    UL = DVECTOR(1, Nobs);
    estp = DVECTOR(1, Nobs);
    Quantal_CI(Nobs, Yp, Yn, 0.95, LL, estp, UL);
    fprintf (fp_out2,"\n\n Data");
    for (i=1;i<=Nobs;i++)
      {
      fprintf (fp_out2,"\n %f %f %f %f", Xi[i], estp[i], LL[i], UL[i]);
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
      /* brat = no if Spec[1] = 1 and Parms[1] < EPS */
      if(brat==1) back=Parms[1];
      back1=1-back;
      if (bmdparm.risk==1) /* if risk type is added */
	back1=1;

      bmdl_bmr_flag = 0;
      bmdu_bmr_flag = 0;
      /* obtain the BMD and BMDL's */
      Multistage_BMD (nparm, Parms, EPS, &junk, xlk, Rlevel, Bmdl, Bmdu, &BMD);

      /* Profile the BMD */

      /*        ProfLik(nparm, Parms, xlk, BMD, Bmdl[1], argv[1],bmdparm.level, bmdparm.effect, */
      /*  	      &ProfFlag,&PEFlag,&PE_bmdl); */
      /* If ProfLik does not plot to the bmdl then try new initial parameters */
      /* ProfFlag    0 = good plot */
      /*             1 = bad plot  */
      /* PEFlag  0 = got an estimate for the bmdl from profile */
      /*         1 = no estimate for bmdl from profile */
      /*        if(ProfFlag == 1) */
      /*  	{ */
      /*  	  for(i=1; i<=nparm; i++) */
      /*  	    { */
      /*  	      Parms[i]= 1; */
      /*  	    } */
      /*  	  ProfLik(nparm, Parms, xlk, BMD, Bmdl[1], argv[1],bmdparm.level, bmdparm.effect,  */
      /*  		  &ProfFlag,&PEFlag,&PE_bmdl);  */
      /*  	} */
      /*        if(ProfFlag == 1)  */
      /*  	{ */
      /*  	  for(i=1; i<=nparm; i++) */
      /*  	    { */
      /*  	      Parms[i]= -1; */
      /*  	    } */
      /*  	  ProfLik(nparm, Parms, xlk, BMD, Bmdl[1], argv[1],bmdparm.level, bmdparm.effect,  */
      /*  		  &ProfFlag,&PEFlag,&PE_bmdl);  */
      /*  	} */
      /*        if (PEFlag == 0) */
      /*  	fprintf(fp_out,"\n         Profile = %14.6g\n",PE_bmdl); */
      /*        if (ProfFlag == 1) */
      /*  	fprintf(fp_out,"\nThe profile plot or the BMDL may be inaccurate.\n\n"); */

    } /* end: if (bmdose==Yes) */

  /* free memory */
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
  FREE_DVECTOR (Bmdu, 1, 5);
  FREE_DMATRIX(vcv,1,nparm,1,nparm);
  FREE_IVECTOR(boundary, 1, nparm);
  FREE_DMATRIX(vcv_adj,1,adj_vcv_rows,1,adj_vcv_rows);

  CLOSE_FILES ();

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
      x = Xi[i];   /* ith dose level */
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
      x = Xi[i];
      ploy=p[nparm];
      for (j=nparm-1; j>=2; j--) ploy=ploy*x+p[j];
      ploy=ploy*x;

      ex = p[1]+(1-p[1])*(1-exp(-ploy)); /* probability of response */
      PROBABILITY_INRANGE(&ex); /* make sure prob is between 0 and 1 */
      ex2=exp(-ploy);

      if (ex != 0)
	 {
         lk1 = Yp[i]/ex - Yn[i]/(1-ex);	 
         lk2 = -Yp[i]/pow(ex, 2) - Yn[i]/pow(1-ex, 2);
         }
      else
         {
         lk1 = 0;
	 lk2 = 0;
         }

      d1[1]= ex2;
      for (j=2;j<=nparm;j++)
	d1[j]= (1-p[1])*ex2*pow(x, j-1); /* first der. of ex wrt beta_j-1   b/c beta_1= parm[2] */


      vcv[1][1] -= lk2*d1[1]*d1[1];
      for (j=2;j<=nparm;j++)
	vcv[1][j] -= lk2*d1[1]*d1[j] - lk1*ex2*pow(x, j-1);

      for (j=2;j<=nparm;j++)
	{
	  if (Spec[j]==0) /* if p[j] is not specified */
	    {
	      for(k=j;k<=nparm; k++)
                {
		vcv[j][k] -= lk2*d1[j]*d1[k] - lk1*(1-p[1])*ex2*pow(x, k+j-2);
                }
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
  int     ii;
  double  ymin, Yt, N, ll, x;
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

  ymin = Xi[1]; /* not used ?????? */
  xmin = Xi[1];
  xmax = 0.0;
  for (i=1;i<=Nobs;i++) /* obtain min and max dose levels */
    {
      x=Xi[i];
      if (x < xmin) xmin = x;
      if (x > xmax) xmax = x;
    } /* end for */
  /** rescale Dose to be: 0 <= Dose <= 1 **/

  scale = xmax;


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
      nanim[i-1] = (Yp[i] + Yn[i]); /* # of animals at ith dose level */
      doses[i-1] = Xi[i]/xmax;                 /* dose levels; rescaled to: 0 <= Dose <= 1  */
      affect[i-1] = (Yp[i]);        /* # of affected animals at ith dose level */
    } /* end for */

  /* calculate the values that will correspond to '0' and 'Infinity'
     in BMD calculations */
  lminbmd = log(DBL_MIN) - log(xmax);
  lmaxbmd = log(DBL_MAX) - log(xmax);
  /* load values for use in getmle, getcl, and getprofile */
  /* into the common block */
  loadcommbloc_(&Nobs,&xmax,nanim,doses,affect,&ndegree,&restrict,
		&lminbmd,&lmaxbmd);


  /* scale parameters */
  for(i = 1; i<=nparm; i++)
    {
      p[i] = p[i]*(pow(xmax,(i-1)));
      if(initial == Yes) IniP[i]= IniP[i]*(pow(xmax,(i-1)));
    }


  /****** Obtain initial estimations for p[] ******/
  if(initial==Yes)
    {
      for(j=1; j<=nparm; j++)
	p[j]=IniP[j];
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
	  if(Yn[i] > 0)
	    Y[i] = log( 1 - Yp[i]/(Yp[i]+Yn[i]) );
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
		    Yt += Yp[i];
		    N += (Yp[i]+Yn[i]);
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
	  if(initial == Yes) IniP[i]= IniP[i]/(pow(xmax,(i-1)));
	}

      /* output initial parameter values */
      OUTPUT_TEXT("\n\n                  Default Initial Parameter Values  ");
      OUTPUT_Init(nparm, Spec, p, Parm_name);

      /* scale back after output */
      for(i = 1; i<=nparm; i++)
	{
	  p[i] = p[i]*(pow(xmax,(i-1)));
	  if(initial == Yes) IniP[i]= IniP[i]*(pow(xmax,(i-1)));
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
    if (Spec[j]==Yes) p[j]=pBak[j];
  /* p[] is now a vector of specified parms and regressed parms */


  for(i = 1; i <= nparm; i++)
    {
      if(initial == Yes && Spec[i] == 0)  /* if the initial option is chosen and */
	{                                   /* if the parameter is not specified   */
	  p[i] = IniP[i];
	  Spec2[i] = 0;
	} /* end if */
      else
	Spec2[i-1] = Spec[i];
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
      doses[i-1]=Xi[i];
    }
  /* unscale parameters */
  for(i = 1; i<=nparm; i++)
    {
      p[i] = p[i]/(pow(xmax,(i-1)));
      if(initial == Yes) IniP[i]= IniP[i]/(pow(xmax,(i-1)));
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
  float MaxLike(int nparm, double p[]);

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
    if (ReportCancerSlopeFactor) {
      fprintf(fp_out, 
#ifndef RBMDS
	      "Cancer Slope Factor =%14.6g\n\n"
#else
	      "Cancer Slope Factor =%30.22g\n\n"
#endif
	      , bmdparm.effect/Bmdl[1]);
    }
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
  double fD, bmdl,  target, *parms, *parms2, xmin, xmax, x;
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

  xmin = Xi[1];
  xmax = 0.0;
  for (i=1;i<=Nobs;i++) /* obtain min and max dose levels */
    {
      x=Xi[i];
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
      Spec2[j-1] = Spec[j];
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
      printf("%f\n",parms2[i]);
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
  double fD, bmdu,  target, *parms, *parms2, xmin, xmax, x;
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

  xmin = Xi[1];
  xmax = 0.0;
  for (i=1;i<=Nobs;i++) /* obtain min and max dose levels */
    {
      x=Xi[i];
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
      Spec2[j-1] = Spec[j];
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
      printf("%f\n",parms2[i]);
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
float MaxLike(int nparm, double p[])
{
  int i, j;
  double like, prob;

  like = 0.0;

  for(i = 1; i <= Nobs; i++)
    {
      prob = p[nparm];
      for(j = nparm-1; j >= 2; j--)
	{
	  prob = Xi[i]*prob + p[j];
	} /* end for */

      prob = Xi[i]*prob;
      prob = p[1] + (1 - p[1])*(1 - exp(-prob));
      if ((prob == 0) || (prob == 1))
	{
	  if (Yp[i] <= 0)
	    like += 0;
	  else
	    like += -1.0e20;
	} /*end if */
      else
	{
	  like += Yp[i]*log(prob) + Yn[i]*log(1 - prob);
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





