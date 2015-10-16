/****************************************************************
 *
 * IMPORTANT NOTE:  The following variable is the version number for
 *                  the current model.  THIS MUST BE CHANGED as
 *				   important changes are made to the models.
 *
 *****************************************************************/

/****************************************************************
char Version_no[] = " $Revision: 2.3 $ $Date: 2001/07/11 22:27:16 $";
char Version_no[]="Weibull Model (Version: 2.7;  Date: 2/20/2007)";
char Version_no[]="Weibull Model (Version: 2.10;  Date: 10/31/2007)";
char Version_no[]="Weibull Model (Version: 2.11;  Date: 03/07/2008)";
char Version_no[]="Weibull Model (Version: 2.12;  Date: 05/16/2008)";
char Version_no[]="Weibull Model (Version: 2.14;  Date: 11/23/2008)";
char Version_no[]="Weibull Model (Version: 2.15;  Date: 10/28/2009)";
*****************************************************************/
char Version_no[]="Weibull Model (Version: 2.16;  Date: 2/28/2013)";

/****************************************************************
 **
 * Weibull.C - a ANSI C program for Weibull model fitting with/without
 *             a natural background rate in Benchmark Dose.
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
 * Version Number: 2.3
 * Modified By: Q. He
 * Modified Date: 09/03/2003
 * Reason:
 *
 * Version Number: 2.5
 * Modified By: Micheal Ferree
 * Modified Date: 8/13/20058
 * Reason: Took out print out of cancer slope factor
 *
 * Version Number: 2.6
 * Modified By: R. W. Setzer
 * Modified Date: 9/29/2005
 * Reason: Free all allocated memory; 
 *
 *
 * Version Number: 2.7
 * Modified By: R. W. Setzer
 * Modified Date: 2/20/2007
 * Reason: Changed version number to reflect changed compilation options.

 * Version Number: 2.8
 * Modified By: G. Nonato
 * Modified Date: 9/23/2007
 * Reason:	1. Changed version number to reflect changed.
 *		2. Prefix condition with "Xi[1] == 0.0 &&" in Line 489 and beaufify code {GLN - 09/23/07, PR0823-05} 
		   It fixes the following:
		       Different incorrect responses result when running models on data for
		       which there is no zero dose group and background is specified to be 0.
		       For the  example data set where doses are 17, 20 and 24, Ns are all 10
		       and responses are 1, 2 and 6
 *
 * Version Number: 2.9
 * Modified By: Woodrow Setzer
 * Modification Date: 10/17/2007
 * Reason:
 *       Changes to allow datasets that have no control dose (that is,
 *       no dose equal to 0.0).
 *
 * Version Number: 2.10
 * Modified By: G. Nonato
 * Modification Date: 10/31/2007
 * Reason: (Per Technical Direction #15, 10/31/07 email)
 *       Fix the background displaying differently in out file (L1713).
 *       Fix the unhandled Win32 Exception (freeing mem variable that was non existing) in L938.
 *
 * Version Number: 2.11
 * Modified By: G. Nonato
 * Modification Date: 03/07/2008
 * Reason: (Per BMDS 2.0: Problem Report 157 & 147)
 *       Fix the Observation # < parameter # for Weibull model problem.
 *       Fix the Invalid model choice.  Defaulting to Weibull model.
 *       Added code to free-up allocated memories before exiting thru ERRORPRT()
 *		Fix the fixed(slope) printing weird # in the Initial/Specified Values
 *
 * Version Number: 2.12
 * Modified By: G. Nonato
 * Modification Date: 05/16/2008
 * Reason: (Per BMDS 2.0: Problem Report 165)
 *       Goodness of Fit - Observed column, print values as real numbers.
 *
 * Version Number: 2.13
 * Modified By: Woodrow Setzer
 * Reason:
 *       Allow lower bound on power to be specified directly (that is, other than
 *       0 or 1).  This allows a direct test of the mabmd software.
 *
 * Version Number: 2.14
 * Modified By: R. Woodrow Setzer
 * Reason:
 *      In calculating initial values, use (a + 1)/(a + b + 2).  This fixes
 *      a case in which the program never returns.
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

/* void ProfLik (int , double *, double, double, double, char *, double); */

/* Forward-declared functions */

void Weibull_fit(int nparm, double p[], double gtol, int *iter,
		 double *fret);
void Weibull_BMD (int nparm, double p[], double gtol, int *iter,
		  double xlk, double Rlevel[], double Bmdl[],double *BMD);
void Weibull_vcv(int nparm, int Spec[], double p[], double **vcv);
void Which_Bounded (int [], double [], int []);
int Model_DF (int []);
double BMDL_func (int, double [], double, double);

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
double    restrict = 0.0;    /* lower bound for shape (power) */
int    initial;	    /* flag for user initial parameter estimates */
int    appendix;    /* flag: 0 for append output, 1 for overwrite */
int    smooth;	    /* flag: 0 for unique bmdl curve, 1 for C-spline */
int    bmdlCurve;   /* flag for bmdl curve option */
double xmax, xmin, scale;

double Rel_Conv, Parm_Conv, Maxloglik;
double SlopeUpperBound = 18.0;

double BMDL_Error_Size;
int BMDL_Error;
int PrintLL = 0; /* if 1, print a line after the anodev table with the */
/* log-likelihood to high precision.  Triggered by a */
/* negative value for ITMAX*/

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
 ** main--main function used to call Weibull mode fitting program.
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

  char    names[3][24] = {"Weibull Model", "Quantal Linear Model", "Quantal Quadratic Model"};
  char    model_name[MNLENGTH], user_note[UNLENGTH];
  char    dose_name[CNLENGTH], posi_name[CNLENGTH], nega_name[CNLENGTH], junkname[FLENGTH];
  int	  select, adj_vcv_rows, *bounded;
  double  **vcv_adj;
  char long_path_name[FLENGTH];
  time_t ltime;

  time( &ltime );

  /* These are defined in float.h */
  //Min_increment = DBL_EPSILON;
  //Max_double = DBL_MAX;

  if(argc < 2)
    {
      fprintf(stderr, "ERROR:  Requires two arguments\nUsage:  %s <file.(d)>\n", argv[0]);
      fprintf (stderr, "   or:  %s -v for version number.\n", argv[0]);
      exit (1);
    } /* end if */

  /********************************************************************
   * {QH 2004/01/14 PR# }
   * Added to show version number if executed from command line with -v
   *********************************************************************/
  if(argc == 2)
    show_version(argv[1], Version_no);

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

    } /* end if */

  /* begin reading input file from batch file (.(d) ext.) */
  fscanf(fp_in, "%s",model_name );
  fscanf(fp_in, "%[ ^\n]",user_note );
  fscanf(fp_in, "%[^\n]", user_note);
  fscanf(fp_in, "%s", junkname);
  fscanf(fp_in, "%s", junkname);
  fscanf(fp_in, "%d",&Nobs);

  /* open the log file if giDo_Log = 1 */
#ifdef DO_LOG
  strcpy(gacLogFile,argv[1]);
  gcDot2 = strchr(gacLogFile, (int) '.');
  (*gcDot2) = (char) 0;
  strcat(gacLogFile,"-Exp.log");
  if (giDo_Log)
    {
      fp_log = fopen(gacLogFile, "w");

      if (fp_log == (FILE *) NULL)
	{
	  ERRORPRT("Unable to open log for Weibull.");
	}

      fprintf(fp_log,"\n\nargv[1] (before) = %s", argv[1]);
    }
#endif

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
  fscanf(fp_in,"%d%lf%lf%d%lf%d%d%d", &ITMAX, &Rel_Conv, &Parm_Conv,
	 &bmdlCurve, &restrict, &bmdose, &appendix, &smooth);
  fscanf(fp_in,"%lf%d%lf",&bmdparm.effect,&bmdparm.risk,&bmdparm.level);

  if (ITMAX < 0)
    {
      ITMAX = -1 * ITMAX;
      PrintLL = 1;
    }
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

      fprintf(fp_out,"Error in opening output files.\n");
      fprintf (fp_out,"...Exited to system!\n");
      FreeUp_mem(Parms, varsum, anasum, vcv);
      exit (1);

    } /* end if */

  select = 1;

#ifdef DO_LOG
  if (giDo_Log)
    {
      fprintf(fp_log,"\n\nmodel_name (before) = %s", model_name);
      fprintf(fp_log,"\nuser_note = %s", user_note);
      fprintf(fp_log,"\nmodel_name[0] = %c", model_name[0]);
      fprintf(fp_log,"\nmodel_name[7] = %c", model_name[7]);
      fprintf(fp_log,"\nmodel_name[8] = %c", model_name[8]);
    }
#endif

  /* Look to see which model the user selected, Weibull, Quantal Linear
     or Quantal Quadratic */
  if (model_name[0] == 'W')
    {
      select = 1;		 /* Case is Weibull */
      strcpy(model_name, names[0]);
#ifndef RBMDS
      fprintf(fp_out2, " %d", 1);
#endif
    }
  else
    {
      if (model_name[0] == 'Q' && model_name[8] == 'L')
	{
	  select = 2;		/* Case is Quantal Linear */
	  strcpy(model_name, names[1]);
#ifndef RBMDS
	  fprintf(fp_out2, " %d", 2);
#endif
	}
      else
	{
	  if(model_name[0] == 'Q' && model_name[8] == 'Q')
	    {
	      select = 3;	/* Case is Quantal Quadratic */
	      strcpy(model_name, names[2]);
#ifndef RBMDS
	      fprintf(fp_out2, " %d", 3);
#endif
	    }
	  else
	    {
	      /* Since the others are just specific cases of the Weibull model,
		 default to Weibull */

	      OUTPUT_TEXT("\n   Invalid model choice.  Defaulting to Weibull model\n");
	      select = 1;
#ifndef RBMDS
	      fprintf(fp_out2, "%d", 1);
#endif

	    } /* end if */

	} /* end if model_name[0] == 'Q' && model_name[7] == 'L'*/

    } /* end if model_name[0] == 'W' */
  strcat(model_name, " using ");
  strcat(model_name, Version_no);

  /* Print model and file information on output page */
  Output_Header(model_name, argv[1], plotfilename, ctime(&ltime), user_note);

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

#ifdef DO_LOG
  if (giDo_Log)
    {
      fflush(fp_log);
      fprintf(fp_log,"\n\nNmiss (before Nobs -= Nmiss) = %d", Nmiss);
      fprintf(fp_log,"\nnparm_known = %d", nparm_known);
      fprintf(fp_log,"\nNobs = %d", Nobs);
    }
#endif

  Nobs -= Nmiss;             /* extern variable Nobs has been changed */

#ifdef DO_LOG
  if (giDo_Log)
    {
      fflush(fp_log);
      fprintf(fp_log,"\n\nNmiss = %d", Nmiss);
      fprintf(fp_log,"\nNobs = %d", Nobs);
      fprintf(fp_log,"\nnparm = %d", nparm);
    }
#endif

  //if (Nobs < nparm), commented and changed to the code below, GLN 2/28/08
  if (Nobs < (nparm-nparm_known))
    {
      FreeUp_mem(Parms, varsum, anasum, vcv);
      switch (select)
	{
	case 1:		/* Weibull model */
	  {
	    ERRORPRT("Observation # < parameter # for Weibull model.");
	  }
	case 2:		/* Quantal Linear model */
	  {
	    ERRORPRT("Observation # < parameter # for Quantal Linear model.");
	  }
	case 3:		/* Quantal Quadratic model */
	  {
	    ERRORPRT("Observation # < parameter # for Quantal Quadratic model.");
	  }

	}  /* end switch */

    }  /* end if */

  /* end of input data */

  /* output title and summary of input data */
  switch (select)
    {
    case 1:		/* Weibull Model */
      {
	OUTPUT_TEXT("\n   The form of the probability function is: ");
	OUTPUT_TEXT("\n   P[response] = background + (1-background)*[1-EXP(-slope*dose^power)]");
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
	    if(Xi[1] == 0.0 && Parms[(int) Background ] == 0.0 && Yp[1] != 0)	//Prefix condition with "Xi[1] == 0.0 &&" {GLN - 09/23/07, PR0823-05}
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
	    if (restrict > 0.0)	/* if lower bound of power parm. is > 0 */
	      {
		fprintf(fp_out,"\n   Power parameter is restricted as power >= %lf", restrict);
	      }
	    else
	      {
		fprintf(fp_out,"\n   Power parameter is not restricted");

	      } /* end if (restrict > 0.0) */

	  }	/* end if (Spec[3]==Yes) */

	break;

      } /* end case 1 */

    case 2:		/* Quantal Linear Model */
      {
	OUTPUT_TEXT("\n   The form of the probability function is: ");
	OUTPUT_TEXT("\n   P[response] = background + (1-background)*[1-EXP(-slope*dose)]");
	fprintf(fp_out,"\n\n   Dependent variable = %s", posi_name);
	fprintf(fp_out,"\n   Independent variable = %s", dose_name);

	/* set power parameter = 1 */
	Spec[(int) Power] = 1;
	Parms[(int) Power] = 1.0;

	if (fixedParm(Background) == Yes)     /* if background parameter is specified */
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
	    if(Parms[(int) Background ] == 0.0 && Yp[1] != 0)
	      {
		FreeUp_mem(Parms, varsum, anasum, vcv);
		ERRORPRT("ERROR:  Background parameter specified as 0, but there were responses\n        at the control dose\n");

	      } /* end if */


	  }	/* end if (Spec[1]==Yes) */

	if (fixedParm(Slope) == Yes)		/* if slope parameter is specified */
	  {
	    fprintf(fp_out,"\n   Slope parameter is set to %g", Parms[2]);

	  } /* end if */

	break;

      } /* end case 2 */

    case 3:		/* Quantal Quadratic Model */
      {
	OUTPUT_TEXT("\n   The form of the probability function is: ");
	OUTPUT_TEXT("\n   P[response] = background + (1-background)*[1-EXP(-slope*dose^2)]");
	fprintf(fp_out,"\n\n   Dependent variable = %s", posi_name);
	fprintf(fp_out,"\n   Independent variable = %s", dose_name);

	/* set power parameter = 2 */
	Spec[(int) Power ] = 1;
	Parms[(int) Power ] = 2.0;

	if (fixedParm(Background) == Yes)    /* if background parameter is specified */
	  {
	    if (Parms[(int) Background ] <= 0.0000001)
	      {
		fprintf(fp_out,"\n   Background parameter is set to zero");
	      }
	    else
	      {
		fprintf(fp_out,"\n   Background parameter is set to %g",
			Parms[1]);

	      } /* end if (Parms[1] <= 0.0000001) */
	    if(Parms[(int) Background ] == 0.0 && Yp[1] != 0)
	      /* assumes control dose comes first */
	      {
		FreeUp_mem(Parms, varsum, anasum, vcv);
		ERRORPRT("ERROR:  Background parameter specified as 0, but there were responses\n        at the control dose\n");

	      } /* end if */

	  }	/* end if (Spec[1]==Yes) */

	if (fixedParm(Slope) ==Yes ) /* if slope parameter is specified */
	  {
	    fprintf(fp_out,"\n   Slope parameter is set to %g", Parms[2]);

	  } /* end if */

	break;

      } /* end case 3 */

    default:
      {
	/* Just in case, default to Weibull */

	OUTPUT_TEXT("\n   Invalid model choice.  Defaulting to Weibull model\n");
	select = 1;
      }
    } /* end of switch */

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
		  FreeUp_mem(Parms, varsum, anasum, vcv);
		  ERRORPRT("When the initial option is chosen, one has to initial ALL unspecified parameters.");

		} /* end if */

	    } /* end if (IniSp[i]==1) */

	} /* end for */

      if (IniP[1] < 0 || IniP[1]>1 || IniP[2]<0 || IniP[3]<restrict)
	{
#define BUFFER_SIZE 128
	  char errbuf[BUFFER_SIZE];
	  snprintf(errbuf, BUFFER_SIZE, 
		   "The initial values have to be: 1 > background >= 0,  slope >= 0 and 18> power >= %lf.", restrict);
	  FreeUp_mem(Parms, varsum, anasum, vcv);
	  ERRORPRT(errbuf);

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

  /* fitting Weibull model and output parameter estimators */
  Weibull_fit(nparm, Parms, EPS, &iter, &xlk);
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
      Weibull_vcv(nparm,Spec,Parms,vcv);

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
    OP_ParmsE(nparm,Spec,Parms,Parm_name,vcv_adj, bounded,bmdparm.level, 1);
  /* compute and output ANOVA table elements */
  DTMS3ANOVA (nparm,Nobs,Spec,lkf,xlk,lkr,anasum, bounded);

  /* output ANOVA table */
  OUTPUT_DTMS3ANOVA(anatxt,anasum);
#ifndef RBMDS
  if (PrintLL == 1) fprintf(fp_out,"\n LogLikelihood: %24.14g\n",xlk);
#else
  if (PrintLL == 1) fprintf(fp_out,"\n LogLikelihood: %30.22g\n",xlk);
#endif
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
#ifdef DO_LOG
	  if (giDo_Log) 
	    {
	      fflush(fp_log);
	      fclose(fp_log);
	    }
#endif

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
      Weibull_BMD (nparm, Parms, EPS, &junk, xlk, Rlevel, Bmdl, &BMD);
      /* Profile the BMD */
      /* reset the BMR first */
      BMR = bmdparm.effect;
      //      ProfLik(nparm, Parms, xlk, BMD, Bmdl[1], argv[1],bmdparm.level);

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
#ifndef RBMDS
	  fprintf(fp_out, "            BMDL =%14.6g\n\n", Bmdl[1]);
#else
	  fprintf(fp_out, "            BMDL =%30.22g\n\n", Bmdl[1]);
#endif
	  /*      fprintf(fp_out, "Cancer Slope Factor =%14.6g\n\n", bmdparm.effect/Bmdl[1]); */
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
  FREE_IVECTOR (bounded, 1, nparm);
  if(ErrorFlag == 0) FREE_DMATRIX (vcv_adj, 1, adj_vcv_rows, 1, adj_vcv_rows);

#ifdef DO_LOG
  if (giDo_Log) 
    {
      fflush(fp_log);
      fclose(fp_log);
    }
#endif

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
	(1 - exp(-(Parms[(int) Slope] * pow(doses[i], Parms[(int) Power]))));
}

/******************************************************************** */
/* D_BMR_Constraint -- compute the derivative of the BMR constraint wrt */
/*                     the other parameters.  That is, if the */
/*                     constraint is Slope = f(Background, Power), */
/*                     compute (df/dBackground, df/dPower). */
/*    input: p[] vector (1-based) of parameter values */
/*    global: tD current guess at BMDL (for BMDL computation) */
/*            BMR BMR value */
/*            bmdparm.risk  either ADDED or EXTRA */
/*    output: dcon vector (length 3, 1 based) such that */
/*            on output, dcon[(int) Background] is df/dBackground , */
/*                       dcon[(int) Power] is df/dPower. */
/*                       dcon[(int) Slope == 0.0 */
/******************************************************************** */
void D_Constraint(double p[], double dcon[])
{
  double powtD, BackComp;

  dcon[(int) Slope] = 0.0;
  powtD = pow(tD,p[(int) Power]);
  BackComp = 1.0 - p[(int) Background];

  if (bmdparm.risk == ADDED)
    {
      dcon[(int) Background ] = BMR / (powtD * BackComp * (BackComp - BMR));
      dcon[(int) Power ] = log(tD) * log(1.0 - BMR/BackComp)/powtD;
    }
  else
    {
      dcon[(int) Background ] = 0.0;
      dcon[(int) Power ] = log(tD) * log(1.0 - BMR)/powtD;
    }
}

/****************************************************************** */
/* D_Predict -- derivatives of the predicted values for each dose in */
/*              doses. */
/*   input: doses - vector (1 based) of doses of length */
/*          ndoses */
/*          Parms - current parameter values at which to compute */
/*                  derivatives; Parms is a vector of 3 elements, */
/*                  regardless of whether a value is specified */
/*   output: Grad[][] - matrix of derivatives whose columns correspond */
/*                      to parameters, and whose rows correspond to */
/*                      doses.  There are only columns for unspecified */
/*                      parameters, in the order Background, Slope, Power. */
/*                      Parameters that have been specified are omitted. */
/****************************************************************** */

void D_Predict(double doses[], int ndoses, double Parms[], double **Grad)
{
  int i, j;
  double DpDx[4], dpow, DSlopeD[4];

  if (replace == Yes) D_Constraint(Parms, DSlopeD);

  for (i = 1; i <= ndoses; i++)
    {
      dpow = pow(doses[i],Parms[(int) Power]);
      if (!fixedParm(Background))
	{
	  if (doses[i] == 0.0)
	    DpDx[(int) Background] = 1.0;
	  else
	    DpDx[(int) Background] =
	      exp(-Parms[(int) Slope] * dpow);
	}

      if (!fixedParm(Slope) || replace == Yes)
	{
	  if (doses[i] == 0.0)
	    DpDx[(int) Slope] = 0.0;
	  else
	    DpDx[(int) Slope] =
	      (1.0 - Parms[(int) Background]) * dpow *
	      exp(-Parms[(int) Slope] * dpow);
	}

      if (!fixedParm(Power))
	{
	  if (doses[i] == 0.0)
	    DpDx[(int) Power] = 0.0;
	  else
	    DpDx[(int) Power] =
	      (1.0 - Parms[(int) Background]) * dpow *
	      Parms[(int) Slope] * log(doses[i]) *
	      exp(-Parms[(int) Slope] * dpow);
	}
      /* Now fill the matrix of Gradients */
      j = 0;
      if (!fixedParm(Background))
	{
	  j++;
	  if (replace == Yes)
	    Grad[i][j] =
	      DpDx[(int) Background] +
	      DpDx[(int) Slope] * DSlopeD[(int) Background];
	  else
	    Grad[i][j] = DpDx[(int) Background];

	}
      if (!fixedParm(Slope))
	{
	  j++;
	  Grad[i][j] = DpDx[(int) Slope];
	}
      if (!fixedParm(Power))
	{
	  j++;
	  if (replace == Yes)
	    Grad[i][j] =
	      DpDx[(int) Power] +
	      DpDx[(int) Slope] * DSlopeD[(int) Power];
	  else
	    Grad[i][j] = DpDx[(int) Power];
	}
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
  /* if replace == Yes, then we are computing BMD; if Probit, */
  /* replace p[3], and if log-Probit, replace p[2] */
  /* with function of current guess at BMD (i.e., tD) */
  if (replace==Yes)
    {
      if(bmdparm.risk==ADDED)
	{
	  /* reparameterize parameter 2 for CI calculation using added risk */
	  p[(int) Slope ] =
	    -log( 1-BMR/(1-p[(int) Background]) )/pow(tD,p[(int) Power]);
	}
      else
	{
	  /* reparameterize parameter 2 for CI calculation using extra risk */
	  /* In Weibull_BMD ck is defined as follows, ck = -log(1-BMR) */
	  p[(int) Slope ] = -log(1 - BMR)/pow(tD,p[(int) Power]);
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
  if (VERYCLOSE(Parms[3],restrict) ||
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
   **Weibull_lk -- used to compute the log likelihood for Weibull model.
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

void Weibull_lk(long int *nvar, double *x, long int *nf, double *f,
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

}	/* end Weibull_lk */

/*******************************************************************
   **Weibull_g -- used to compute the gradients for Weibull model.
   *
   *		Global var.:
   *			Nobs, Xi, Yp, Yn, replace, Spec, ck, bmdparm.risk
   *		Input:
   *			nvar is the number of parameters
   *			nf is not currently used
   *			g is a vector of length nparm to hold the parameters
   *			uiparm is a vector of the indices for the parameters
   *				of length nparm
   *			urparm is a vector of fixed parameters
   *				of length jparm
   *			ufparm is not currently used
   *		Output:
   *			f is the likelihood value
   *		Input/Output:
   *			x is a vector of non-fixed parameters of length jvar
   *
   **********************************************************************/
/* This is the analytic version */

void Weibull_g (long int *nvar, double *x, long int *nf, double *g,
		long int *uiparm, double *urparm, void (*ufparm)())
{
  int      i, j;
  double   xlk;     		    /* log likelihood */
  double   *p;
  double   **Grad, *Pred;

  /* parameters for calculation */

  p = DVECTOR(1,nparm);
  Pred = DVECTOR(1, Nobs);
  Grad = DMATRIX(1, Nobs, 1, (int) *nvar);

  /* construct the parameter vector */
  unpack(x, urparm, p);

  /* Get the current predicted values */
  Predict(Xi, Nobs, p, Pred);

  /* Get the gradient of the model */

  D_Predict(Xi, Nobs, p, Grad);

  for (i = 0; i < *nvar; i++) g[i] = 0.0;

  for (i=1;i<=Nobs;i++)
    {
      xlk = -(Yp[i] * D_Slog(Pred[i]) - Yn[i] * D_Slog(1.0 - Pred[i]));
      j = 0;
      if (!fixedParm(Background))
	{
	  g[j] += xlk * Grad[i][j + 1];
	  j++;
	}
      if (!fixedParm(Slope))
	{
	  g[j] += xlk * Grad[i][j + 1];
	  j++;
	}
      if (!fixedParm(Power))
	{
	  g[j] += xlk * Grad[i][j + 1];
	}
    }	/* end for */

	/* free allocated memory */
  FREE_DVECTOR(p, 1, nparm);
  FREE_DVECTOR(Pred, 1, Nobs);
  FREE_DMATRIX(Grad, 1, Nobs, 1, (int) *nvar);

}	/* end Weibull_g */
/****************************************************************************
   **
   * Weibull_grad -- Computes the gradient of the Weibull likelihood
   * function with respect to the user form of the parameters.  This is to
   * be used in Weibull_vcv, to compute a finite difference approximation to
   * the hessian of the likelihood function
   * Input: nparm -- the number of parameters, the dimension of all the following
 arrays.
 Spec[] -- if the ith parameter is fixed by the user, Spec[i] == 1,
 otherwise Spec[i] == 0.
 ptf[] -- vector of parameters, in user form (external form),
 based at 1.
 * Output: grad[] -- the gradient of the loglikelihood function
 (N.B. Weibull_g returns the gradient of -loglikelihood)
 Based at 1.

  *****************************************************************************/
void Weibull_grad(int nparm, int Spec[], double ptf[],
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
  Weibull_g (&nvar, start, &nf, outgrad, uiparm, urparm, ufparm);
  /* Weibull_g returns the gradient of -loglikelihood */
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
 *	Weibull_vcv -- used to compute the vcv matrix for Weibull model.
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

void Weibull_vcv (int nparm, int Spec[], double ptf[], double **vcv)
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
	  Weibull_grad(nparm, Spec, saveparms, gradp);
	  saveparms[i] = ptf[i] - h[i];
	  Weibull_grad(nparm, Spec, saveparms, gradm);
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
} /* end of Weibull_vcv */

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

      Weibull_lk( &nvar, start, &dummy, fret,
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
      lower[jvar] = restrict;
      upper[jvar] = SlopeUpperBound;

    } /* end if (!fixedParm(Power)) */
  ErrorFlag = run_dmngb((int) nvar, start, lower, upper, Maxloglik,
			Rel_Conv, Parm_Conv, ITMAX, 10,
			Weibull_lk,Weibull_g,uiparm,urparm,ufparm,
			0,fret);

  *fret = -*fret;

  /* put the parameter values back in p */
  unpack(start, urparm, p);
  /* free allocated memory */
  free(urparm);
  free(start);

}	/* end MAX_lk */

/**************************************************************
 *Weibull_fit -- Used to "prepare" the data for further computation,
 *            i.e. compute the extern variables, give the initial
 *            parameters, etc. THEN fit the Weibull model.
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
void Weibull_fit(int nparm, double p[], double gtol,
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
      p[(int) Slope ] = p[(int) Slope ]*pow(scale,-p[(int) Power]);
    }
  else
    {
      /* compute initial estimates */
      int contdose, maxdose;
      double ymax;

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
      ymax = (Yp[maxdose]+1)/(Yp[maxdose] + Yn[maxdose] + 2.0);

      if (!fixedParm(Background))
	{
	  p[(int) Background ] = (Yp[contdose]+1.0)/
	    (Yp[contdose] + Yn[contdose] + 2.0);
	  /* If xmin > 0, then this is really the response at the lowest non-control
	     dose; just halve this estimate.
	  */
	  if (xmin > 0) p[(int) Background ] *= 0.5;
	}
      if (!fixedParm(Slope))
	{
	  W = (ymax - p[(int) Background ]) / (1.0 - p[(int) Background ]);
	  if (W <= 0) W = 0.5/(Yp[maxdose]+Yn[maxdose]+1);
	  if (W >= 1) W = 1.0 - (0.5/(Yp[maxdose]+Yn[maxdose]+1));

	  p[(int) Slope] = -log(1.0 - W);
	}
      if (!fixedParm(Power))
	{
	  double pihat, Dhat;
	  int bri;

#define YP(i) ((Yp[i] + 1) / (Yp[i] + Yn[i] + 2))

	  pihat = 0.5 * (p[(int) Background] + YP(maxdose));

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

	  W = 1.0 - ((pihat - p[(int) Background])/(1.0 - p[(int) Background]));
	    W = -log(W)/p[(int) Slope];
	  if (W > 0.0) {
	    p[(int) Power ] = log(W)/log(Dhat);
	  } else {
	    p[(int) Power ] = 1.0;
	  }

	  if (p[(int) Power ] > 0.8*SlopeUpperBound)
	    p[(int) Power ] = 0.8*SlopeUpperBound;

	  if (restrict == Yes && p[ (int) Power ] <= 1)
	    {
	      p[(int) Power ] = 1.0;
	    }
	}
      /* If the slope is fixed, we have to rescale it, anyway */
      if (fixedParm(Slope))
	{
	  p[(int) Slope] *= pow(scale, -p[(int) Power]);
	}
      /* Now do something silly, so that we can print out the Slope */
      /* on the right scale */
      //GLN commented out the if, as per Woody's suggestion (see next line): 03/07/08
      //if (!fixedParm(Slope)) p[(int) Slope] *= pow(scale,p[(int) Power]);
      p[(int) Slope] *= pow(scale,p[(int) Power]);
      OUTPUT_TEXT("\n\n                  Default Initial (and Specified) Parameter Values  ");
      OUTPUT_Init(nparm, Spec, p, Parm_name);
      //GLN commented out the if, as per Woody's suggestion (see next line): 03/07/08
      //if (!fixedParm(Slope)) p[(int) Slope] *= pow(scale,-p[(int) Power]);
      p[(int) Slope] *= pow(scale,-p[(int) Power]);

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
  p[(int) Slope] = p[(int) Slope] * pow(scale, p[(int) Power]);
  scale = 1.0;

}	/* end Weibull_fit */

/************************************************************
 * Weibull_BMD -- Used to calculate the BMD and BMDL for Weibull model.
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
void Weibull_BMD (int nparm, double p[], double gtol, int *iter, double xlk,
		  double Rlevel[], double Bmdl[], double *BMD)
{
  double   tol;
  double   xa,xb,fa,fb;
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
  /* IF ML is the value of the maximized log-likelihood, then ML - LR is the value
     log-likelihood at the BMDL or BMDU */
  if (bmdparm.level<0.5)
    {
      LR = 0.5*QCHISQ(1.0 - 2.0 * bmdparm.level, 1);
    }
  else
    {
      LR = 0.5*QCHISQ(2.0 * bmdparm.level - 1.0, 1);

    } /* end if */

  Rlevel[1] = BMR = bmdparm.effect;

  /***** Solve the BMD ******************************************/

  /* If there is no slope, we want to cap the BMD estimate. */
  /* We still want to compute a lower limit on the BMD, however */
  bogusBMD = 0;
  if (p[3] <= 0.0 || p[(int) Slope] <= 0.0)
    {
      *BMD = 100 * xmax;
      bogusBMD = 1;
      Warning(" Slope or Power parameter essentially zero.  BMD set to 100 * max(Dose).\n");

    }
  else
    {
      if (bmdparm.risk == ADDED)
	{
	  xb = pow(-log(1.0 - BMR/(1.0 - p[(int) Background]))/p[(int) Slope],1.0/p[3]);
	}
      else
	{
	  xb = pow(-log(1.0 - BMR)/p[(int) Slope],1.0/p[3]);
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
	  /* *BMD = -2.0; */
	  /* return; */
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
  fb = DBL_MAX;
  for (i = 1; i <= nparm; i++) pa[i] = pb[i] = p[i];
  fa = BMDL_func(nparm, pa, xa, tol);  /* if fa > 0, then there is a BMDL
					  somewhere in the range */

  while (fa<0.0 && xa > DBL_MIN && fabs(fb - fa) > DBL_EPSILON)
    {
      xb = xa;
      fb = fa;
      for (i = 1; i <= nparm; i++) pb[i] = pa[i];
      xa *= stepsize;

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
		  xb = pow(-log(1.0 - BMR/(1.0 - p[(int) Background]))/p[(int) Slope],1.0/p[3]);
		}
	      else
		{
		  xb = pow(-log(1.0 - BMR)/p[(int) Slope],1.0/p[3]);
		}
	      BMDtmp = xb;
	    }
	  /* First get the interval to search over */
	  xa = xb * stepsize;
	  tol = FMAX(xb*0.0001, 0.0000001);

	  BMD_lk = xlk;		/*get the lk at BMD.*/
	  fb = DBL_MAX;
	  for (i = 1; i <= nparm; i++) pa[i] = pb[i] = p[i];
	  fa = BMDL_func(nparm, pa, xa, tol);

	  while (fa < 0.0 && xa > DBL_MIN && fabs(fb - fa) > DBL_EPSILON) /* BMDL does not lie between xa and xb, so reduce xa */
	    {
	      xb = xa;
	      fb = fa;
	      for (i = 1; i <= nparm; i++) pb[i] = pa[i];
	      xa *= stepsize;   /*prevent that xa=stepsize*BMD is not small eno*/
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
  for(j=1; j<=nparm; j++)
    {
      p[j]= pBak[j];

    } /* end for */

  /* free allocated memory */
  FREE_DVECTOR(pBak,1, nparm);
  FREE_DVECTOR(pa,1, nparm);
  FREE_DVECTOR(pb,1, nparm);

}	/* end Weibull_BMD */


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
float BMDL_func(int nparm, double p[], double D, double gtol)
{
  /* ck , BMD_lk and LR are calculated in Weibull_BMD() */

  double fD, xlk;
  int junk;
  tD = D; /* tD is global var. have to change before call MAX_lk().*/

  MAX_lk( nparm, p, gtol,  &junk,  &xlk);

  fD = BMD_lk - xlk - LR;  /* constraint imposed by the likelihood function */

  return fD;

}   /* end BMDL_func */

