/*******************************************************************
*
* IMPORTANT NOTE:  The following variable is the version number for
*				   the current model.  THIS MUST BE CHANGED as
*	               important changes are made to the models.
*
********************************************************************/
char Version_no[]="Hill Model. (Version: 2.17;  Date: 01/28/2013)";
/*
char Version_no[]="Hill Model. (Version: 2.16;  Date: 04/06/2011)";
char Version_no[]="Hill Model. (Version: 2.15;  Date: 10/28/2009)";
char Version_no[]="Hill Model. (Version: 2.14;  Date: 06/26/2008)";
char Version_no[]="Hill Model. (Version: 2.12;  Date: 02/20/2007)";
char Version_no[]="Hill Model. (Version: 2.10;  Date: 12/08/2006)";
*/
/*******************************************************************
*
* Hill.C - a ANSI C program for Hill model fitting with/without a
*          natural background rate in Benchmark Dose.
*
* Date: Oct, 2000
*
********************************************************************
* Modification Log:
*
* Version Number: 2.1.1
* Modified By: Qun He
* Modified Date: 3/09/2003
* Reason:
*
* Version Number: 2.4
* Modified By: Micheal Ferree
* Modified Date: 6/24/2005
* Reason: Fixed degrees of freedom problems.  Changed output to
*		  display all Tests of interest with note about Test 2 and 3
*         being the same if rho=0.  Changed confidence level from
*         0.05 to 0.1 for tests 2, 3, and 4.  Will try to clean code
*         for version 2.3.1.
*
* Version Number: 2.5
* Modified By: R. Woodrow Setzer
* Modified Date: 10/21/2005
* Reason: Fixed degrees of freedom problems at the source (bind in getmle.f)
*         Release all allocated memory
*         switched to using cdft to calculate confidence limits
*         allow conditional compilation, so RBMDS version does not produce
*           .002 file
*         allow conditional compilation to print to console (MISC_OUT)
*         cleaned up comments around BMDL_func.
* Outstanding problems: BMDL calculation is still unstable.
*
* Version Number: 2.6
* Modified By: R. Woodrow Setzer
* Modified Date: 03/21/2006
* Reason: move calculation of likelihoods for A1, A2, and R to
*         compute_continuous_liks().
* Outstanding problems: BMDL calculation is still unstable.
*
* Version Number: 2.7
* Modified By: R. Woodrow Setzer
* Modified Date: 03/23/2006
* Reason: Allow A3 likelihood to use fixed variance parameters, if available.
* Outstanding problems: BMDL calculation is still unstable.
*
* Version Number: 2.8
* Modified By: R. Woodrow Setzer
* Modified Date: 07/11/2006
* Reason: Remove redundancy in 'power parameter not restricted' header info
* Outstanding problems: BMDL calculation is still unstable.
* 
* Version Number: 2.9
* Modified By: R. Woodrow Setzer
* Modified Date: 08/03/2006
* Reason: Change variance model parameterization; clean up constraint 
*         specification for donlp2; change from equality to inequality
*         likelihood constraint for BMDL calculation.  Other cosmetic
*         fixups.
*
* Version Number: 2.10
* Modified By: R. Woodrow Setzer
* Modified Date: 12/08/2006
* Reason: fixed bug in initial value for constant variance case:
*         was finding log(alpha), and reporting using it as alpha
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
* Modification Date: 04/07/2008
* Reason: (Per BMDS 2.0: Problem Report 157 & 147)
*       Fix the Observation # < parameter # for Hill model problem.
*       Added code to free-up allocated memories before exiting thru ERRORPRT()
*
* Version Number: 2.14
* Modified By: Woodrow Setzer
* Date: 06/26/2008
* Reason: 1) impossible to fix n to 18 (17.9999 works, though)
*         2) User input specified values not echoing when returning
*            initial parameter values.
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
* Modification Date: 04/06/2011
* Reason: PR384
*      To Fix for the bug in reading parameter file caused by 
*      longer file name length and spaces within the name.
*
* Version Number: 2.17
* Modified by: Louis Olszyk
* Modification Date: 01/18/2013
* Reason: - PR461 - replace "Absolute risk", "Relative risk" and
*           "Point risk" in output with "Absolute deviation",
*           "Relative deviation" and "Point estimate", respectively.
*         - PR444 - Improved plot titles
********************************************************************/

#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
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

extern void getmle_(long int *nparms, double parms[], double parms2[], double *ll,
		    long int *optite, long int *nresm, long int bind[]);

extern void getcl_(long int *which, long int *nparm, double *bmr, double *bmd,
		   double *target, double parms[], long int *risktype, double *bmdl,
		   double parms2[], long int *optite, long int *nresm, long int bind[],
		   long int *flag);

extern void getmlea3_(long int *nparm, double parms[], double parms2[], double *ll,
		      long int *optite, long int *nresm, long int bind[]);

extern void loadcommbloc_(long int *ndoses, double doses[], double means[],
			  long int nanimals[], double svar[], long int *nparm,
			  long int fixed[], double fixedval[], long int *restrict,
			  long int *adverse, long int *model_type, double *xmax,
			  double *xmin);


void GetMoreParms (double *p, int size);
void GetOtherParms (double *p, int size);
double power (double a, double b);
void GetNewParms (double *p, int size);
void Get_BMRS (double *p, double Dose, double bmdl1, double *BMRVals,
	       int sign, int bmr_type);
void HillMeans (int nobs, double p[], double Doses[], double means[]);
void OUTPUT_BENCHMD2 (double BMD);
void GetCLParms (double *p, int size);
void MeanPart (int obs, double *p, double *mg);
void Var2Part (int obs, int const_var, double Vi, double meani, double *p,
	       double *mg, double **mg2, double **vg2);
void VarPart (int obs, int const_var, double Vi, double meani, double *p,
	      double *mg, double *vg);
void Mean2Part (int obs, double *p, double **mg2);
double BMDL_func (int nparm, double xlk, double Dose, double pBak[]);



#define EPS 3.0e-8
#define EPSS 3.0e-7
#define FPMIN 1.0e-30
#define MAXIT 100

#define float double

/* By GLN-04/06/2011, added log file for debugging */

//#define DO_LOG

int    giDo_Log = false;     /*  Uncomment to activate switch for log file   */
char   gacLogFile[FLENGTH], *gcDot2;
//#define MISC_OUT

/* Define input and putput files's name */
char fin[FLENGTH];			/* input temp file */
char fout[FLENGTH];			/* output temp file */
char fout2[FLENGTH];
char plotfilename[FLENGTH];		/* file to pass to GnuPlot */
char *Parm_name[] = { "alpha", "rho", "intercept", "v", "n", "k" };
char *alt_Parm1_name = "lalpha";
char *anatxt[] = { "Full model", "Fitted model", "Reduced model" };


/* variables will not be changed except Spec.
   xxi and yyi will be sorted together */

int *Spec, *SpecBkp;		/* vector used to identify user input parm */
int *IniSp;			/* vector to identify initialized parms */
int *Ni;			/* number of animals in the i-th dose group */
double *Ym;			/* mean response data array */
double *Yd;			/* sample variance (s_i^2) of response data array */
double *Xi;			/* independent variable (dose value) data array */
double *xxi;			/* dose values when not divided to dose groups */
double *yyi;			/* response values when not divided to dose groups */
double *Ysum;			/* sum of responses within a group */
double *Rlevel;
double *Bmdl;
double *IniP;
int in_type, Nobs, ntotal, ndeg, nparm, restrict, initial, appendix, smooth;
int bmdlCurve, sign;
double xmax, xmin, scale;

int cons_var, bmr_type;			/* deviations:  0 if absolute, 1 if std dev,
				   2 if relative, 3 if point, 4 if extra */
double BMD_lk, LR, BMR, react;

/* GLN - 04/07/2008
*  Free-up allocated memory before exit
*  upon encountering fatal error.
*/
//void FreeUp_mem(double *Parms, VarList *varsum, AnaList *anasum, double  **vcv)
void FreeUp_mem(double *Parms, int *bounded, VarList *varsum, AnaList *anasum, int var_type, int tempnparm)
{
	FREE_IVECTOR (bounded, 1, nparm);
	FREE_IVECTOR (Ni, 1, Nobs);
	FREE_DVECTOR (Xi, 1, Nobs);
	FREE_DVECTOR (Ym, 1, Nobs);
	FREE_DVECTOR (Yd, 1, Nobs);

	FREE_DVECTOR (Parms, 1, tempnparm);
	FREE_DVECTOR (IniP, 1, tempnparm);
	FREE_IVECTOR (Spec, 1, tempnparm);
	FREE_IVECTOR (IniSp, 1, tempnparm);
	FREE_VLVECTOR (varsum, 1, 3);

	if (var_type == 0)
	{
		FREE_ALVECTOR (anasum, 1, 5);
	}
	else
	{
		FREE_ALVECTOR (anasum, 1, 5);
	}

	FREE_DVECTOR (Rlevel, 1, 5);
	FREE_DVECTOR (Bmdl, 1, 5);
	
	if (fp_log != (FILE *) NULL)
		fclose(fp_log);

	return;
}

/****************************************************************
 *	main--main function used to call Hill mode fitting program.
 *		  Includes: biosubcc.c--common subfunction C program.
 *****************************************************************/
int main (int argc, char *argv[])
{

  void Hill_vcv (int nparm, int Spec[], double p[], double **vcv);
  int READ_OBSDATA4V (int Nobs, double Xi[], int Ni[], double Ym[],
		      double Yd[]);
  int READ_OBSDATA2V (int ntotal, double xxi[], double yyi[]);
  void AThree_Fit (int nparm, double p[], double gtol, int *iter,
		   double *fret);
  void Hill_BMD (int nparm, double p[], double gtol, int *iter, double xlk,
		 double Rlevel[], double Bmdl[], double *BMD);
  void Get_Second_BMR (double bmdl, double *BMRVals);
  void Hill_fit (int nparm, double p[], int *is_conv,
		 int *iter, double *fret, int *bounded);
  int Model_DF(int []);


  int iter, i, j, jj, junk;	/* iteration variable */
  int group;			/* temp counter eventually equals Nobs */

  int bmdose;			/* flag for computing benchmark dose */
  int Nmiss;			/* number of records with missing values */
  int nparm_known, poly_known;	/* number of specified parameters */
  int adj_vcv_rows;		/* Number of rows/columns in an adjusted vcv matrix */
  double xlk, lkA3, lkA1, lkA2, lkR;	/* log likelihoods */
  double BMD, lep, upep;
  double *stdev;		/* Sample Standard deviation */


  double *Parms, *LKParms;	/* parameter array, parameter array for likelihood test */
  VarList *varsum;		/* info for variables--p. dep.,n. dep., indep. */
  AnaList *anasum;		/* information for ANOVA analysis */
  double **vcv;			/* variance and covariance matrix */
  double **vcv_adj;		/* adjusted vcv matrix taking into account
				   /  any bounded parameters */

  char model_name[MNLENGTH], user_note[UNLENGTH];
  char dose_name[CNLENGTH], no_name[CNLENGTH], mean_name[CNLENGTH], stdev_name[CNLENGTH],
    junkname[FLENGTH], junkname2[FLENGTH];
  char response_name[CNLENGTH];

  int tmpi1, var_type, conv_check;
  int tempnparm, *bounded, temp_sign;
  double Rel_Conv, Parm_Conv;
  double Nd;

  double *doses, *means, *svar, *parms;  /* These are for loading in the common block */
  long int *nanim, *Spec2;
  long int adverse;      /* adverse direction flag to go into the common block */
  long int model_type;   /* model_type = 0 for polynomial */
  long int NobsL, nparmL, restrictL; /* These are just copies to send to FORTRAN */
  char long_path_name[FLENGTH];          /* in long int format                       */
  double *mean, *std;
  time_t ltime;

  /* Set time zone from TZ environment variable. If TZ is not set,
     /  the operating system is queried to obtain the default value
     /  for the variable  */

  /*_tzset();*/
  time (&ltime);

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
    path_name (argv[1]);*/
  if (argc > 2)
  {
	  path_name2(argc, argv, long_path_name);
      argv[1] = long_path_name;
  }


  fp_in = fopen (argv[1], "r");
  if (fp_in == NULL)
    {
      fprintf (stderr, "Error in opening input  file.\n");
      fprintf (stderr, "...now exiting to system...\n");

      exit (1);
    }

	#ifdef DO_LOG
	  strcpy(gacLogFile,argv[1]);
	  gcDot2 = strchr(gacLogFile, (int) '.');
	  (*gcDot2) = (char) 0;
	  strcat(gacLogFile,"-Hil.log");
	  if (giDo_Log)
		{
		  fp_log = fopen(gacLogFile, "w");

		  if (fp_log == (FILE *) NULL)
		{
		  ERRORPRT("Unable to open log for Hill.");
		}

		  fprintf(fp_log,"\n\nargv[1] = %s", argv[1]);
		  fprintf(fp_log,"\nlong_path_name = %s", long_path_name);
		  float flen = FLENGTH;
		  fprintf(fp_log,"\nFLENGTH = %g", flen);
		}
	#endif

  /* Hillmodelfilename.(d) currently contains the wrong number
     /  of parameter values in some instances.  There should be
     /  6 parameters passed, but sometimes only five are passed
     /  on lines 8 and 10 of the .(d) file.  The follwing bit of
     /  code checks to see if line 8 contains 5 or 6 parameters */
/*
  for (i = 1; i <= 7; i++)
    {
      fgets (junkname, 100, fp_in);
    }

  fgets (junkname2, 100, fp_in);

  i = 0;
  nparm = 0;
  while (junkname2[i] != '\0')
    {
      if (junkname2[i] != ' ')
	{
	  nparm++;
	  while ((junkname2[i] != ' ') && (junkname2[i] != '\0'))
	    {
	      i++;
	    }
	}
      else
	{
	  i++;
	}
    }
	rewind (fp_in);
	*/
  /* Go back to the beginning of the file */

  nparm = 6;

  fscanf (fp_in, "%s", model_name);
  fscanf (fp_in, "%[ ^\n]", user_note);
  fscanf (fp_in, "%[^\n]", user_note);
  fscanf (fp_in, "%s", junkname);
  fscanf (fp_in, "%s", junkname);
  fscanf (fp_in, "%d", &in_type);

  /* in_type = 1 if input format is Xi, Ni, Ym, Yd. */
  /* in_type = 0 if input format is Ntotal, Xi, Y_ij. */

  if (in_type == 1)
      fscanf (fp_in, "%d", &Nobs);
  else
      fscanf (fp_in, "%d", &ntotal);

  fscanf (fp_in, "%d", &sign);

  fscanf (fp_in, "%d%lf%lf%d%d%d%d%d", &ITMAX, &Rel_Conv, &Parm_Conv,
	  &bmdlCurve, &restrict, &bmdose, &appendix, &smooth);

  /* restrict=0 if no restriction, +1 if beta_j>0, -1 if beta_j<0 */
  fscanf (fp_in, "%d%lf%d%lf", &bmr_type, &bmdparm.effect, &cons_var,
	  &bmdparm.level);

#ifdef DO_LOG
  if (giDo_Log)
    {
      fprintf(fp_log,"\n\nmodel_name (before) = %s", model_name);
      fprintf(fp_log,"\nuser_note = %s", user_note);
      fprintf(fp_log,"\nin_type = %d", in_type);
      fprintf(fp_log,"\nNobs = %d", Nobs);
      fprintf(fp_log,"\nntotal = %d", ntotal);
      fprintf(fp_log,"\nITMAX = %d", ITMAX);
      fprintf(fp_log,"\nRel_Conv = %g", Rel_Conv);
      fprintf(fp_log,"\nParm_Conv = %g", Parm_Conv);
      fprintf(fp_log,"\nbmdlCurve = %d", bmdlCurve);
      fprintf(fp_log,"\nrestrict = %d", restrict);
      fprintf(fp_log,"\nbmdose = %d", bmdose);
      fprintf(fp_log,"\nappendix = %d", appendix);
      fprintf(fp_log,"\nsmooth = %d", smooth);
      fprintf(fp_log,"\nbmr_type = %d", bmr_type);
      fprintf(fp_log,"\nbmdparm.effect = %g", bmdparm.effect);
      fprintf(fp_log,"\ncons_var = %d", cons_var);
      fprintf(fp_log,"\nbmdparm.level = %g", bmdparm.level);
    }
#endif
  /* allocate memory for arrays, however, for the Parms
     /  array et al, they should be of size 6, but the .(d)
     /  file may only have 5 to read in.  Don't change nparm
     /  yet, but allocate correct space */

  tempnparm = 6;
  Parms = DVECTOR (1, tempnparm);
  Spec = IVECTOR (1, tempnparm);
  IniP = DVECTOR (1, tempnparm);
  IniSp = IVECTOR (1, tempnparm);
  varsum = VLVECTOR (1, 3);
  Rlevel = DVECTOR (1, 5);
  Bmdl = DVECTOR (1, 5);

  bmdparm.risk = 1;

  Get_Names (argv[1], fout, fout2, plotfilename);

  if (appendix == Yes)
    {
      fp_out = fopen (fout, "a");
    }
  else
    {
      fp_out = fopen (fout, "w");
    }
#ifndef RBMDS
  fp_out2 = fopen (fout2, "w");
#endif

  if ((fp_out == NULL) 
#ifndef RBMDS
      || (fp_out2 == NULL)
#endif
      )
    {
      printf ("Error in opening  output files.\n");
      printf ("...now exiting to system...\n");

      fprintf (fp_out, "Error in opening output files.\n");
      fprintf (fp_out, "...Exited to system!\n");
      exit (1);
    }

  /* Print model and file information on output page */
  Output_Header (Version_no, argv[1], plotfilename, ctime (&ltime),
		 user_note);

  if ((bmdose < 0) || (bmdose > 1))
    {
      ERRORPRT ("Error in choosing benchmark dose computation.");
    }

  /* Vectors need to be of size 6, but the input file may only
     /  have 5 coming in */

  /*obtain user input parameters */
  READ_PARAMETERS (nparm, Parms);
#ifdef MISC_OUT
  for (i=0; i < nparm; i++) printf("\n###Specified %s = %f\n", Parm_name[i], Parms[i+1]); /*DEBUG*/
#endif
  FILL_SPECVECTOR (nparm, Parms, Spec);
  poly_known = 0;
#ifdef MISC_OUT
  for (i=0; i < nparm; i++) printf("\n###Is_specified? %s = %d\n", Parm_name[i], Spec[i+1]); /*DEBUG*/
#endif
  for (i = 3; i <= nparm; i++)
    {
      if (Spec[i] == 1)
	{
	  poly_known += 1;
	}
    }

  /* obtain user input initial parameters values */
  fscanf (fp_in, "%d", &initial);
  READ_PARAMETERS (nparm, IniP);
#ifdef MISC_OUT
  for (i=1; i <= nparm; i++) printf("\n###IniP[%d] = %f\n", i, IniP[i]); /*DEBUG*/
#endif
  FILL_SPECVECTOR (nparm, IniP, IniSp);
#ifdef MISC_OUT
  for (i=1; i <= nparm; i++) printf("\n###IniSp[%d] = %d\n", i, IniSp[i]); /*DEBUG*/
#endif

  /* If input file had 5 parameters instead of six, add the sixth one */
  if (nparm == 5)
    {
      Parms[6] = -9999;
      Spec[6] = 0;
      IniP[6] = -9999;
      IniSp[6] = 0;
      nparm++;
    }
  /* Why do we do this? RWS 6/26/2008 */
  for (i = 1; i <= nparm; i++)
    {
      if (Spec[i] == 1)
	{
	  IniP[i] = Parms[i]; /* 1; */
	}
    }

  var_type = 1 - cons_var;	/* vartype is 1 if cons_var is 0, 0 if
				   /  cons_var is 1 */

  /* Set the variance power parameter as specified at 0 if constant
     /  variance is selected */

  if (cons_var == 1)
    {
      Spec[2] = 1;
      Parms[2] = 0;
    } else {
      Parm_name[0] = alt_Parm1_name;
    }

  /* Anasum vector is different in the two variance cases as
     /  there are different tests that one may be interested in */

  if (Spec[2] == 1)
    {
      if (Parms[2] == 0)
	{
	  anasum = ALVECTOR (1, 5);
	}
      else
	{
	  anasum = ALVECTOR (1, 5);
	}
    }
  else
    {
      anasum = ALVECTOR (1, 5);
    }

  /* obtain observation data into Yp, Yn, Xi, Ls, Xg vectors */
  if (in_type == 1)
    {
      fscanf (fp_in, "%s%s%s%s", dose_name, no_name, mean_name, stdev_name);
    }
  else
    {
      fscanf (fp_in, "%s%s", dose_name, response_name);
    }

  if (in_type == 1)
    {
      Ym = DVECTOR (1, Nobs);
      Yd = DVECTOR (1, Nobs);
      Xi = DVECTOR (1, Nobs);
      Ni = IVECTOR (1, Nobs);

      /* The following three lines of code were added/changed because
         /  the user is asked to enter in standard deviation, but the
         /  program treats entry as a variance */

      stdev = DVECTOR (1, Nobs);

      Nmiss = READ_OBSDATA4V (Nobs, Xi, Ni, Ym, stdev);
      for (i = 1; i <= Nobs; i++)
	{
	  Yd[i] = stdev[i] * stdev[i];	/* Yd[i] = ith sample variance */
	}

      Nobs -= Nmiss;		/* extern variable Nobs has been changed */

      FREE_DVECTOR (stdev, 1, Nobs);

    }
  else
    {
      xxi = DVECTOR (1, ntotal);
      yyi = DVECTOR (1, ntotal);
      Ym = DVECTOR (1, ntotal);
      Yd = DVECTOR (1, ntotal);
      Xi = DVECTOR (1, ntotal);
      Ni = IVECTOR (1, ntotal);
      Ysum = DVECTOR (1, ntotal);
      Nmiss = READ_OBSDATA2V (ntotal, xxi, yyi);
      ntotal -= Nmiss;		/* extern variable Nobs has been changed */

      for (i = 1; i <= ntotal; i++)
	{
	  Xi[i] = Ysum[i] = Ym[i] = Yd[i] = 0;
	  Ni[i] = 0;
	}

      Sort_2_By_Dose (ntotal, xxi, yyi);	/* Sort arrays so that the
						   /  following code works.  Located
						   /  in specialfun.h */
      group = 1;

      for (i = 1; i <= ntotal; i++)
	{
	  Xi[group] = xxi[i];
	  Ni[group] += 1;
	  Ysum[group] += yyi[i];
	  if (i < ntotal)
	    {
	      if (xxi[i] != xxi[i + 1])
		{
		  group += 1;
		}
	    }
	}

      Nobs = group;
      jj = 1;

      for (i = 1; i <= Nobs; i++)
	{
	  Ym[i] = Ysum[i] / Ni[i];
	  if (Ni[i] > 1)
	    {
	      for (j = 1; j <= Ni[i]; j++)
		{
		  Yd[i] +=
		    (yyi[jj] - Ym[i]) * (yyi[jj] - Ym[i]) / (Ni[i] - 1);
		  jj += 1;
		}
	    }
	  else
	    {
	      Yd[i] = 0.0;
	    }
	}

      FREE_DVECTOR (xxi, 1, ntotal);
      FREE_DVECTOR (yyi, 1, ntotal);

    }				/* end if (in_type == 1) */


	if (Nobs < (4-poly_known))
	{
	  FreeUp_mem(Parms, bounded, varsum, anasum, var_type, tempnparm);
	  ERRORPRT ("Observation # < parameter # for Hill model.");
	}

  if (in_type == 1)
    {
      for (i = 1; i <= Nobs; i++)
	{
	  if ((Xi[i] < 0) || (Ni[i] < 0) || (Yd[i] < 0))
	    {
	      ERRORPRT
		("Values of dose, group size, and sample variance should be positive...");
	    }
	}
    }
  else
    {
      for (i = 1; i <= Nobs; i++)
	{
	  if (Xi[i] < 0)
	    {
	      ERRORPRT ("Dose value should be positive ...");
	    }
	}
    }

  /* end of input data */


  xmin = Max_double;
  xmax = 0.0;
  for (i = 1; i <= Nobs; i++)
    {
      if (Xi[i] < xmin)
	{
	  xmin = Xi[i];
	}
      if (Xi[i] > xmax)
	{
	  xmax = Xi[i];
	}
    }

  /* Get initial default adverse direction by finding the general
     /  linear trend of the data */

  if ((sign != 1) || (sign != -1))
    {
      sign = Get_Linear_Trend (Nobs, Xi, Ym, Ni);
      temp_sign = 0;
    }
  else
    {
      temp_sign = 9999;
    }				/* end if ((sign != 1) || (sign != -1)) */

  /* This is all for loading the common block for donlp2 */
  doses = DVECTOR(0, Nobs-1);
  means = DVECTOR(0, Nobs-1);
  svar = DVECTOR(0, Nobs-1);
  nanim = LIVECTOR(0, Nobs-1);
  parms = DVECTOR(0, nparm-1);
  Spec2 = LIVECTOR(0, nparm-1);
  for(i = 1; i <= Nobs; i++)
    {
      nanim[i-1] = (long int) Ni[i];
      doses[i-1] = Xi[i]/xmax;
      means[i-1] = Ym[i];
      svar[i-1] = Yd[i];
    }
  for(i=0; i<nparm; i++)
    {
      Spec2[i] = Spec[i+1];
      parms[i] = Parms[i+1];
    }
  /* If k is fixed, we need to rescale it */
  if (Spec2[5] == Yes) parms[5] = parms[5]/xmax;
  NobsL = Nobs;
  nparmL = nparm;
  restrictL = restrict;
  adverse = sign;
  model_type = 2;

  loadcommbloc_(&NobsL, doses, means, nanim, svar, &nparmL, Spec2,
		parms, &restrictL, &adverse, &model_type, &xmax, &xmin);

  FREE_DVECTOR(doses, 0, Nobs-1);
  FREE_DVECTOR(means, 0, Nobs-1);
  FREE_DVECTOR(svar, 0, Nobs-1);
  FREE_LIVECTOR(nanim, 0, Nobs-1);
  FREE_DVECTOR(parms, 0, nparm-1);

  /********************************************************************
   * {QH 2003/12/31 PR# }
   * This redundancy caused a lot of problems.
   *********************************************************************/
  /*FREE_DVECTOR(parms, 0, nparm-1); */

  FREE_LIVECTOR(Spec2, 0, nparm-1);
  /* End of load common block */

  /* output title and summary of input data */
  OUTPUT_TEXT ("\n   The form of the response function is: ");
  OUTPUT_TEXT ("\n   Y[dose] = intercept + v*dose^n/(k^n + dose^n)");

  if (in_type == 1)
    {
      fprintf (fp_out, "\n\n   Dependent variable = %s", mean_name);
    }
  else
    {
      fprintf (fp_out, "\n\n   Dependent variable = %s", response_name);
    }

  fprintf (fp_out, "\n   Independent variable = %s", dose_name);

  for (i = 1; i <= nparm; i++)
    {
      if (Spec[i] == 1)
	{
	  fprintf (fp_out, "\n   %s is set to %g", Parm_name[i - 1],
		   Parms[i]);
	}
    }

  if (restrict == 1)
    {
      fprintf (fp_out,
	       "\n   Power parameter restricted to be greater than 1");
    }
  else if (restrict == 0)
    {
      fprintf (fp_out,
	       "\n   Power parameter is not restricted");
    }

  if (cons_var == 1)
    {
      fprintf (fp_out, "\n   A constant variance model is fit");
    }
  else
    {
      fprintf (fp_out,
       "\n   The variance is to be modeled as Var(i) = exp(lalpha  + rho * ln(mean(i)))");
    }

  nparm_known = COUNT_SPECVECTOR (nparm, Spec);

  fprintf (fp_out, "\n\n   Total number of dose groups = %d", Nobs + Nmiss);
  fprintf (fp_out, "\n   Total number of records with missing values = %d",
	   Nmiss);

  fprintf (fp_out, "\n   Maximum number of iterations = %d\n", ITMAX);
  fprintf (fp_out, "   Relative Function Convergence has been set to: %g\n",
	   Rel_Conv);
  fprintf (fp_out, "   Parameter Convergence has been set to: %g\n\n",
	   Parm_Conv);


  if ((Rel_Conv != 1.0e-8) || (Parm_Conv != 1.0e-8))
    {
      fprintf (fp_out,
	       "****  We are sorry but Relative Function and Parameter Convergence    ****\n");
      fprintf (fp_out,
	       "****  are currently unavailable in this model.  Please keep checking  ****\n");
      fprintf (fp_out,
	       "****  the web sight for model updates which will eventually           ****\n");
      fprintf (fp_out,
	       "****  incorporate these convergence criterion.  Default values used.  ****\n\n");
    }


  if (initial == Yes)
    {
      OUTPUT_TEXT
	("\n\n                 User Inputs Initial Parameter Values  ");
      OUTPUT_Init (nparm, Spec, IniP, Parm_name);
#ifdef MISC_OUT
  for (i=1; i <= nparm; i++) printf("\n###IniP[%d] = %f\n", i, IniP[i]); /*DEBUG*/
#endif
  FILL_SPECVECTOR (nparm, IniP, IniSp);
#ifdef MISC_OUT
  for (i=1; i <= nparm; i++) printf("\n###IniSp[%d] = %d\n", i, IniSp[i]); /*DEBUG*/
#endif

      for (i = 1; i <= nparm; i++)
	{
	  if (IniSp[i] == 1)	/* have been initialized */
	    {
	      if (Spec[i] == 1)	/* check if it is for fixed parm */
		{
		  Warning
		    ("The initial value for the fixed parameter is ignored.");
		}
	    }
	  else
	    {			/*check if all the unspecified parms were initialized */
	      if (Spec[i] == 0)
		{
		  ERRORPRT
		    ("You have to initialize either ALL or NONE of the unspecified parameters.");
		}
	    }
	}			/* end for (i = 1; i <= nparm; i++) */

      if (IniP[1] <= 0)
	{
	  ERRORPRT ("The initial value of variance has to be positive.");
	}
      tmpi1 = 0;
      if (restrict == 1)
	{
	  if (IniP[5] < 1)
	    {
	      tmpi1 = 1;
	    }
	  if (tmpi1 != 0)
	    {
	      ERRORPRT ("You have restricted n > 1, but initialized n < 1.");
	    }
	}
    }				/* end if initial == Yes */


  /* compute likelihoods for A1, A2, and R*/
  lkA1 = lkA2 = lkA3 = lkR = 0.0;
  compute_continuous_liks(Nobs, Ni, Ym, Yd, &lkA1, &lkA2, &lkR);

  /* Compute Likelihood for model A3: Yij = Mu(i) + e(ij)
     Var(e(ij)) = k*(m(xi))^rho
     Only need this if a non-constant variance model is run */
  if (var_type == 1 && (Spec[2] == 0 || (Spec[2] == 1 && Parms[2] != 0.0)))
    {

      LKParms = DVECTOR (1, Nobs + 2);	/* Parameters for fitting the */
                                        /* model above.               */

      nparm = Nobs + 2;

      if (nparm > 25)
	{
	  ERRORPRT ("\n*****\nToo many parameters\nPlease reduce the size of data set!\n*****\n");
	}

      AThree_Fit (nparm, LKParms, EPS, &iter, &lkA3);

      nparm = 6;

      FREE_DVECTOR (LKParms, 1, Nobs + 2);

    }  /* end if var_type == 1 */

  /* If constant variance model then just set equal to model A1. */
  /* MJF, 25MAY2005. */
  else {
    lkA3 = lkA1;
  }

  /* fitting Hill model and output parameter estimates */

  bounded = IVECTOR (1, nparm);

  Hill_fit (nparm, Parms, &conv_check, &iter, &xlk, bounded);

  /* If a bad completion code is returned (conv_check == -1), then
     exit the program */


  if (conv_check == -1)
    {
#ifndef RBMDS
      fprintf (fp_out2, "END  %d", -1);	/* Let plotting program know */
#endif
      /* that it didn't make it    */
      /* need to free allocated memory before exiting: RWS 9/23/2005 */
      /*      exit (0); */
    } else {

    /* Check maximization results to get a better default adverse
       direction */

    if ((Parms[4] < 0) && (temp_sign == 0))
      {
	sign = -1;
      }
    else if ((Parms[4] > 0) && (temp_sign == 0))
      {
	sign = 1;
      }

    /* Compute the aprox. covariance matrix */

    vcv = DMATRIX (1, nparm, 1, nparm);

    INITIALIZE_DMATRIX (vcv, nparm, nparm);

    /* Creates matrix of second partial derivatives of the negative
       of the log likelihood function */

    Hill_vcv (nparm, Spec, Parms, vcv);


    /*  initialize adj_vcv matrix for Get_and_OUTPUTDTMSVCV() */

    adj_vcv_rows = 0;

    for (i = 1; i <= nparm; i++)
      {
	if (bounded[i] == 0)
	  {
	    adj_vcv_rows++;
	  }
      }				/* end for (i = 1; i <= nparm; i++) */

    vcv_adj = DMATRIX (1, adj_vcv_rows, 1, adj_vcv_rows);

    /* Finish creating the asymptotic correlation matrix, and output it */
    Get_and_OUTPUT_DTMSVCV (nparm, Spec, Parm_name, vcv, vcv_adj, bounded);

    /* Output parameter values and standard errors (if not bounded) */
    OP_ParmsE (nparm, Spec, Parms, Parm_name, vcv_adj, bounded, bmdparm.level, 1);


    /* Calculate fitted means and standard deviations, and print them out in a table with */
    /* the data, for comparison. */
    mean = DVECTOR (1, Nobs);
    std = DVECTOR (1, Nobs);

    for (i = 1; i <= Nobs; i++)
      {
	if (Xi[i] == 0)
	  {
	    mean[i] = 0;
	  }
	else
	  {
	    mean[i] =
	      Parms[4] * pow (Xi[i], Parms[5]) / (pow (Parms[6], Parms[5]) + pow (Xi[i], Parms[5]));
	  }

	mean[i] = Parms[3] + mean[i];

	if (Parms[2] == 0)
	  {
	    std[i] = sqrt(Parms[1]);
	  }
	else
	  {
	    std[i] = sqrt(exp(Parms[1] + log(fabs (mean[i])) * Parms[2]));
	  }
      }
    PrintData(mean, std, Xi, Ym, Yd, Ni, Nobs);

    FREE_DVECTOR(mean, 1, Nobs);
    FREE_DVECTOR (std, 1, Nobs);

    /* compute ANOVA table elements */
    /* Always send 1 for var_type to DTMS3ANOVAC, MJF 24JUN05. */
    DTMS3ANOVAC(nparm,Nobs,Spec,lkA3,xlk,lkA2,lkA1,lkR,anasum,1,bounded);


    /* output ANOVA table elements */
    OUTPUT_DTMS3ANOVAC (anatxt,anasum,var_type);

    /* calculate and print a goodness of fit table */
    Goodness (nparm, nparm_known, Parms, var_type, anasum);

    /* compute benchmark dose */
#ifndef RBMDS
    fprintf (fp_out2, "\n BMD_flag \t %d \n Nobs \t%d \n nparm \t%d", bmdose,
	     Nobs, nparm);
    fprintf (fp_out2, "\n  Con_lev \t%3.3g ", bmdparm.level);
  fprintf (fp_out2, "\n  BMRType \t%d ", bmr_type);
  fprintf (fp_out2, "\n  BMRF \t%f ", bmdparm.effect);
#endif

#ifndef RBMDS
    /* print parameters */
    for (i = 1; i <= nparm; i++)
      {
	fprintf (fp_out2, "\n %s \t %5.3g", Parm_name[i - 1], Parms[i]);
      }

    /* print data */
    fprintf (fp_out2, "\n\n Data");
    for (i = 1; i <= Nobs; i++)
      {
	Nd = Ni[i];
	lep = Ym[i] + qstudt(0.025, (Ni[i] - 1)) * sqrt (Yd[i] / Ni[i]);
	upep = Ym[i] + qstudt(0.975, (Ni[i] - 1)) * sqrt (Yd[i] / Ni[i]);
	fprintf (fp_out2, "\n %f %f %f %f", Xi[i], Ym[i], lep, upep);
      }
    
    fprintf (fp_out2, "\n Max_Min_dose \n  %f %f ", xmax, xmin);
#endif
    if (bmdose == Yes)
      {
	Hill_BMD (nparm, Parms, EPS, &junk, xlk, Rlevel, Bmdl, &BMD);

      }				/* end if (bmdose == Yes) */
    /* Memory allocated in this branch */
    FREE_DMATRIX (vcv, 1, nparm, 1, nparm);
    FREE_DMATRIX(vcv_adj, 1, adj_vcv_rows, 1, adj_vcv_rows);
  }  /* end if (conv_check != -1) */

  /* free memory */
  FREE_IVECTOR (bounded, 1, nparm);
  FREE_IVECTOR (Ni, 1, Nobs);
  FREE_DVECTOR (Xi, 1, Nobs);
  FREE_DVECTOR (Ym, 1, Nobs);
  FREE_DVECTOR (Yd, 1, Nobs);

  FREE_DVECTOR (Parms, 1, tempnparm);
  FREE_DVECTOR (IniP, 1, tempnparm);
  FREE_IVECTOR (Spec, 1, tempnparm);
  FREE_IVECTOR (IniSp, 1, tempnparm);
  FREE_VLVECTOR (varsum, 1, 3);

  if (var_type == 0)
    {
      FREE_ALVECTOR (anasum, 1, 5);
    }
  else
    {
      FREE_ALVECTOR (anasum, 1, 5);
    }

  FREE_DVECTOR (Rlevel, 1, 5);
  FREE_DVECTOR (Bmdl, 1, 5);
  

  CLOSE_FILES ();
  return (0);

}				/*end of main */

/******************************************************************
 *	Hill_vcv -- used to compute the vcv for Polynomial model.
 *		Extern var.: smean, smax, Nobs, Xi, Yp, Yn, Ls, Xg.
 *
 ******************************************************************/
void
Hill_vcv (int nparm, int Spec[], double p[], double **vcv)
{
  void F1iDoublePart (int nparm, int const_var, double p[], double **Fn1i,
		      int obs);

  void F2iDoublePart (int nparm, int const_var, double p[], double **Fn2i,
		      int obs);

  void F3iDoublePart (int nparm, int const_var, double p[], double **Fn3i,
		      int obs);

  void MeanPart (int obs, double *p, double *mg);
  void VarPart (int obs, int const_var, double Vi, double meani, double *p,
		double *mg, double *vg);
  void Mean2Part (int obs, double *p, double **mg2);
  void Var2Part (int obs, int const_var, double Vi, double meani, double *p,
		 double *mg, double **mg2, double **vg2);

  double **Fn1i, **Fn2i, **Fn3i, numi;
  int i, j, k, const_var;

  Fn1i = DMATRIX (1, nparm, 1, nparm);
  Fn2i = DMATRIX (1, nparm, 1, nparm);
  Fn3i = DMATRIX (1, nparm, 1, nparm);

  /* Compute partials at parameter estimates */

  if ((Spec[2] == 1) && (p[2] == 0))
    {
      const_var = 1;
    }
  else
    {
      const_var = 0;
    }

  for (i = 1; i <= Nobs; i++)
    {
      numi = Ni[i];
      for (j = 1; j <= nparm; j++)
	{
	  for (k = 1; k <= nparm; k++)
	    {
	      F1iDoublePart (nparm, const_var, p, Fn1i, i);
	      F2iDoublePart (nparm, const_var, p, Fn2i, i);
	      F3iDoublePart (nparm, const_var, p, Fn3i, i);

	      vcv[j][k] = vcv[j][k] + (numi * Fn1i[j][k] / 2);
	      vcv[j][k] += (numi - 1) * Yd[i] * Fn2i[j][k] / 2;
	      vcv[j][k] += numi * Fn3i[j][k] / 2;

	    }
	}
    }

  FREE_DMATRIX (Fn1i, 1, nparm, 1, nparm);
  FREE_DMATRIX (Fn2i, 1, nparm, 1, nparm);
  FREE_DMATRIX (Fn3i, 1, nparm, 1, nparm);

}				/* end Hill_vcv */

/**************************************************************
 *	Hill_fit -- Used to "prepare" the data for further computation,
 *				i.e. give the initial
 *				parameters, etc. THEN fit the Hill model.
 *				(In fact, these jobs could be done in main().)
 *
 ***************************************************************/
void
Hill_fit (int nparm, double p[], int *is_conv,
	  int *iter, double *fret, int *bounded)
{
  int i, j, Ntot, Nisum;
  double *pBak, *parms, *fitparms, *fitparms2, *fitparms3, *doses;
  double ll, ll2, ll3, temp2;
  double *beginp;
  double slope;
  double ymax, ymin, *svar, *means, temp, maxdose, varmean;
  long int nvar, *nanim, nparms, *Spec2, restr, signs;
  long int optite, optite2, optite3, nresm, *bind, *bind2, *bind3, model_type;

  pBak = DVECTOR (1, nparm);

  /* Hill model */

  model_type = 2;

  ymin = Max_double;
  ymax = -1.0e30;
  for (i = 1; i <= Nobs; i++)
    {
      if (Ym[i] < ymin)
	{
	  ymin = Ym[i];
	}
      if (Ym[i] > ymax)
	{
	  ymax = Ym[i];
	}
    }

  /* rescale Dose to be: 0 <= Dose <= 1 */
  scale = 1;


  for (j = 1; j <= nparm; j++)
    {
      pBak[j] = p[j];		/* save the input p[] */
    }

  /* Obtain initial estimations for p[] */
  if (initial == Yes)
    {
      for (j = 1; j <= nparm; j++)
	{
	  p[j] = IniP[j];
	}
    }
  else
    {				/* The following block of code gives initial starting
				   /  values as the ordinary least squares solution */

      if (Spec[3] != 1)
	{
	  p[3] = Ym[1];
	}

      if (Spec[4] != 1)
	{
	  if (sign != -1)
	    {
	      p[4] = ymax - p[3];
	    }
	  else
	    {
	      p[4] = ymin - p[3];
	    }
	}


      if (Spec[6] != 1) {
	i = 1;
	p[6] = 0;
	temp = p[3] + p[4] / 2;
	
	if (sign == 1)
	  {
	    while (Ym[i] < temp)
	      {
		i++;
	      }
	    
	    slope = (Ym[i - 1] - Ym[i]) / (Xi[i] - Xi[i - 1]);
	    p[6] = ((temp - Ym[i]) / slope) + Xi[i];
	  }
	else
	  {
	    while (Ym[i] > temp)
	      {
		i++;
	    }
	    slope = (Ym[i] - Ym[i - 1]) / (Xi[i] - Xi[i - 1]);
	    p[6] = (temp - Ym[i]) / slope + Xi[i];
	  }
	
	if (p[6] < 0)
	  {
	    p[6] = 0;
	  }
      }
      
      if (Spec[5] != 1)
	{
	  p[5] = 0.0;
	  Ntot = 0;
	  for (i = 1; i <= Nobs; i++)
	    {
	      temp = p[4] / (Ym[i] - p[3]) - 1;
	      if ((temp > 0) && (Xi[i] != 0) && (p[6] > 0))
		{
		  temp2 = log (temp) / log (p[6] / Xi[i]);
		  if (temp2 > 0)
		    {
		      p[5] = p[5] + temp2;
		    }
		  Ntot++;
		}
	    }
	  
	  if ((Ntot > 0) && (p[5] > 0))
	    {
	      p[5] = p[5] / Ntot;
	    }
	  else
	    {
	      p[5] = 1.0;
	    }

	  if (p[5] > 18)
	    {
	      p[5] = 18;
	    }
	}

      /* Get initial values for p[1] and p[2]: alpha and rho */
      if (Spec[2] != 1)
	p[2] = 0.0;
      if (Spec[1] != 1)
	{
	  varmean = 0.0;
	  Nisum = 0;
	  for (i = 1; i <= Nobs; i++) 
	    {
	      varmean += Yd[i]*(Ni[i] - 1);
	      Nisum += Ni[i];
	    }
	  if (cons_var == 0)
	    p[1] = log(varmean/(Nisum - Nobs));
	  else
	    p[1] = varmean/(Nisum - Nobs);
	    
	}
      

      OUTPUT_TEXT
	("\n\n                  Default Initial Parameter Values  ");
      OUTPUT_Init (nparm, Spec, p, Parm_name);
    }				/* end if (initial == Yes) */

  maxdose = Xi[1];
  for (i = 2; i <= Nobs; i++)
    {
      if (maxdose < Xi[i])
	{
	  maxdose = Xi[i];
	}
    }

  /*This region scales all the doses down to between 0 and 1 */
  for (i = 1; i <= Nobs; i++)
    {
      Xi[i] = Xi[i] / maxdose;
    }

  /* get specified parameters */
  for (j = 1; j <= nparm; j++)
    {
      if (Spec[j] == Yes)
	{
	  p[j] = pBak[j];
	}
    }
  /* rescale p[6] (k) if it is specified */
  if (Spec[6] == Yes)
    {
      p[6] = p[6] / maxdose;
    }
  nvar = Nobs;
  nparms = nparm;
  restr = restrict;
  signs = sign;
  doses = (double *) malloc ((size_t) (Nobs) * sizeof (double));
  means = (double *) malloc ((size_t) (Nobs) * sizeof (double));
  svar = (double *) malloc ((size_t) (Nobs) * sizeof (double));
  nanim = (long int *) malloc ((size_t) (nvar) * sizeof (long int));
  parms = (double *) malloc ((size_t) (nparm) * sizeof (double));
  fitparms = (double *) malloc ((size_t) (nparm) * sizeof (double));
  fitparms2 = (double *) malloc ((size_t) (nparm) * sizeof (double));
  fitparms3 = (double *) malloc ((size_t) (nparm) * sizeof (double));
  Spec2 = (long int *) malloc ((size_t) (nparm) * sizeof (long int));
  bind = (long int *) malloc ((size_t) (nparm) * sizeof (long int));
  bind2 = (long int *) malloc ((size_t) (nparm) * sizeof (long int));
  bind3 = (long int *) malloc ((size_t) (nparm) * sizeof (long int));
  beginp = (double *) malloc ((size_t) (nparm) * sizeof (double));

  for (i = 1; i <= Nobs; i++)
    {
      nanim[i - 1] = (long int) Ni[i];
      doses[i - 1] = Xi[i];
      means[i - 1] = Ym[i];
      svar[i - 1] = Yd[i];
    }

  for (i = 1; i <= nparm; i++)
    {
      if ((initial == Yes) && (Spec[i] == 0))
	{
	  p[i] = IniP[i];
	  Spec2[i - 1] = 0;
	}
      else
	{
	  Spec2[i - 1] = Spec[i];
	}
      parms[i - 1] = p[i];
    }				/* end for (i = 1; i <= nparm; i++) */

  /* This scales the slope parameter by the maxdose */
  parms[5] = parms[5] / maxdose;

  for (i = 0; i < nparm; i++)
    {
      beginp[i] = parms[i];
    }

  /* This is the first call to getmle and internally in donlp2   */
  /* the parameters are scaled by the abs of their initial value */
  /* So this will give the next call different parameters to start w/ */
  getmle_(&nparms,parms,fitparms2,&ll2,&optite2,&nresm,bind2);

  /* This uses the fitparms from the scaled call as the starting values. */
  getmle_(&nparms,fitparms2,fitparms3,&ll3,&optite3,&nresm,bind3);

  getmle_ (&nparms, parms, fitparms, &ll, &optite, &nresm, bind);
  /*    getmle_ (&nvar, doses, means, nanim, svar, &nparms, parms, Spec2, parms, */
  /*  	   &restr, &signs, fitparms, &ll, &optite, &nresm, bind, &model_type); */

  if ((optite < 0) || (optite > 3)) /* was 2 */
    {
      for (i = 0; i < 30; i++)
	{
#ifdef MISC_OUT
	  printf ("\noptite = %ld", optite);
#endif
	  if (optite != 3)
	    {
	      GetNewParms (beginp, nparm);	/* Get a new starting point */
	    }
	  else
	    {
	      for (j = 0; j < nparm; j++)
		{
		  beginp[j] = fitparms[j];
		}
	    }

	  /* Try again */

	  getmle_ (&nparms, beginp, fitparms, &ll, &optite, &nresm, bind);
	  /*  	  getmle_ (&nvar, doses, means, nanim, svar, &nparms, beginp, Spec2, */
	  /*  		   beginp, &restr, &signs, fitparms, &ll, &optite, &nresm, */
	  /*  		   bind, &model_type); */

	  if ((optite >= 0) && (optite <= 3)) /* was 2 */
	    {
	      break;
	    }
	}
    }				/* end if ((optite < 0) || (optite > 2)) */


  if ((optite < 0) || (optite > 3)) /* was 2 */
    {
      for (i = 0; i < 10; i++)
	{
#ifdef MISC_OUT
	  printf ("\n   optite = %ld", optite);
#endif
	  if (optite != 3)
	    {
	      GetMoreParms (beginp, nparm);	/* Get a new starting point */
	    }
	  else
	    {
	      for (j = 0; j < nparm; j++)
		{
		  beginp[j] = fitparms[j];
		}
	    }

	  /* Try again */

	  getmle_ (&nparms, beginp, fitparms, &ll, &optite, &nresm, bind);
	  /*  	  getmle_ (&nvar, doses, means, nanim, svar, &nparms, beginp, Spec2, */
	  /*  		   beginp, &restr, &signs, fitparms, &ll, &optite, &nresm, */
	  /*  		   bind, &model_type); */

	  if ((optite >= 0) && (optite <= 3)) /* was 2 */
	    {
	      break;
	    }
	}
    }				/* end if ((optite < 0) || (optite > 2)) */

  if ((optite < 0) || (optite > 3)) /* was 2 */
    {
      for (i = 0; i < 10; i++)
	{
	  if (optite != 3)
	    {
	      GetOtherParms (beginp, nparm);	/* Get a new starting point */
	    }
	  else
	    {
	      for (j = 0; j < nparm; j++)
		{
		  beginp[j] = fitparms[j];
		}
	    }

	  /* Try again */

	  getmle_ (&nparms, beginp, fitparms, &ll, &optite, &nresm, bind);
	  /*  	  getmle_ (&nvar, doses, means, nanim, svar, &nparms, beginp, Spec2, */
	  /*  		   beginp, &restr, &signs, fitparms, &ll, &optite, &nresm, */
	  /*  		   bind, &model_type); */

	  if ((optite >= 0) && (optite <= 3 )) /* was 2 */
	    {
	      break;
	    }
	}
    }				/* end if ((optite < 0) || (optite > 2)) */


  /* This decides if the scaling model is better than the unscaled model */
  /* or not. */
  if ((optite2 >= 0) && (optite2 <= 3)) /* was 2 */
    {
      if (ll2 > ll)
	{
	  for (j = 0; j < nparm; j++)
	    {
	      fitparms[j] = fitparms2[j];
	      bind[j] = bind2[j];
	    }
	  optite = optite2;
	  ll = ll2;
	}
    }
  if ((optite3 >= 0) && (optite3 <= 3)) /* was 2 */
    {
      if (ll3 > ll)
	{
	  for (j = 0; j < nparm; j++)
	    {
	      fitparms[j] = fitparms3[j];
	      bind[j] = bind3[j];
	    }
	  optite = optite3;
	  ll = ll3;
	}
    }


  /* Let user know if no optimum was found */
  if ((optite < 0) || (optite > 3)) /* was 2 */
    {
      fprintf (fp_out,
	       "\n\n!!! Warning:  optimum may not have been found.                      !!!");
      fprintf (fp_out,
	       "\n!!! Bad completion code in maximum likelihood optimization routine  !!!");
      fprintf (fp_out,
	       "\n!!! Program halted                                                  !!!\n\n");
      *is_conv = -1;
    }
  else
    {
      *is_conv = 1;
    }

  for (i = 1; i <= nparm; i++)
    {
      p[i] = fitparms[i - 1];
      bounded[i] = bind[i - 1];
    }

  for (i = 1; i <= Nobs; i++)
    {
      Xi[i] = Xi[i] * maxdose;
    }

  p[6] = p[6] * maxdose;

  /* This seems to be code devoted to duplicating what is (or should be) already in bounded.
     Now that info is being transferred from bind correctly, does this get better?
   
     bounded[1] = 0;
     if ((fabs (p[2] - 18) < 1e-20) || (fabs (p[2] + 18) < 1e-20))
     bounded[2] = 1;
     else
     bounded[2] = 0;
     bounded[3] = 0;
     bounded[4] = 0;
     bounded[5] = 0;
     bounded[6] = 0;
     if ((restrict == 1) && (fabs (p[5] - 1) < 1e-20))
     bounded[5] = 1;
     if ((restrict != 1) && (fabs (p[5] - 0) < 1e-20))
     bounded[5] = 1;
     if (fabs (p[6] - .00000001) < 1e-20)
     bounded[6] = 1;
  */
  *fret = ll;

  FREE_DVECTOR (pBak, 1, nparm);
  free(doses);
  free(means);
  free(svar);
  free(nanim);
  free(parms);
  free(fitparms);
  free(fitparms2);
  free(fitparms3);
  free(Spec2);
  free(bind);
  free(bind2);
  free(bind3);
  free(beginp);
  
}				/* end Hill_fit */


/************************************************************
 *	Hill_BMD -- Used to calculate the BMD and BMDL for Hill model.
 *
 *************************************************************/
void
Hill_BMD (int nparm, double p[], double gtol, int *iter, double xlk,
	  double Rlevel[], double Bmdl[], double *BMD)
{
  /*    float Binary_root (float (*O_func) (int, float[], float, float), */
  /*  		     float xa, float xb, float fa, float fb, float tol, */
  /*  		     int nparm, float p[]); */
  float BMD_func (int nparm, float p[], float x, float gtol);


  double tol;
  double xa, xb, fa, fb, bDose;
  double D, bmrtemp, temp;
  double *pBak, *BMRVals;
  double incre, divide, thedose, tempdose[2], bmdlmean[2];
  int j, k, bmrfncsign;


  pBak = DVECTOR (1, nparm);
  BMRVals = DVECTOR (1, 5);

  for (j = 1; j <= nparm; j++)
    {
      pBak[j] = p[j];		/* save the p[] */
    }

  /* compute chi-squared value */
  /* IF ML is the value of the maximized log-likelihood, then ML - LR is the value
     log-likelihood at the BMDL or BMDU */
  if (bmdparm.level < 0.5)
    {
      LR = 0.5*QCHISQ(1.0 - 2.0 * bmdparm.level, 1); 
    }
  else
    {
      LR = 0.5*QCHISQ(2.0 * bmdparm.level - 1.0, 1);
    }
  if ((bmr_type != 3) && (bmr_type != 4))
    {
      Rlevel[1] = BMR = fabs (bmdparm.effect);
    }
  else
    {
      Rlevel[1] = BMR = bmdparm.effect;
    }

  xa = D = 0.0;
  tol = 0.0000001;
  incre = 1.0;
  divide = 50.0;
  bmrfncsign = 1;
  fa = -1;
  thedose = 0;

  /* Get a BMR value that is equal to just Mu(BMD) regardless
     /  of the BMR type to simplify the root finding routine */

  if ((bmr_type == 0) || (bmr_type == 1) || (bmr_type == 2))
    {
      if (sign == -1)
	{
	  bmrtemp = -BMR;
	}
      else
	{
	  bmrtemp = BMR;
	}

      if (bmr_type == 1)
	{
	  if (p[2] != 0)
	    {
	      bmrtemp = bmrtemp * sqrt (exp(p[1] + log(p[3])*p[2]));
	    }
	  else
	    {
	      bmrtemp = bmrtemp * sqrt (p[1]);
	    }
	}
      else if (bmr_type == 2)
	{
	  bmrtemp = bmrtemp * p[3];
	}
    }
  else
    {
      if (bmr_type == 3)
	{
	  bmrtemp = BMR - p[3];
	}
      else
	{
	  bmrtemp = p[4] * BMR;
	}
    }				/*  */

  bmrtemp = bmrtemp + p[3];	/* bmrtype is now a Point BMR type
				   /  regardless of the user input
				   /  type */

  temp = BMR;
  BMR = bmrtemp;
  BMRVals[1] = BMR;

  /* Find the BMD */
  /* first see if the bmr value is appropriate */
  if (p[4] > 0)
    {
      if ((BMR < p[3]) || (BMR > p[3] + p[4]))
	{
	  fprintf (fp_out, "\n\nBMD computation failed");
	  fprintf (fp_out,
		   "\nBMR value is not in the range of the mean function");
	  FREE_DVECTOR(pBak, 1, nparm);
	  FREE_DVECTOR(BMRVals, 1, 5);
	  return ;
	}			/* end if */
    }				/* end if */
  if (p[4] < 0)
    {
      if ((BMR > p[3]) || (BMR < p[3] + p[4]))
	{
	  fprintf (fp_out, "\n\nBMD computation failed");
	  fprintf (fp_out,
		   "\nBMR value is not in the range of the mean function");
	  FREE_DVECTOR(pBak, 1, nparm);
	  FREE_DVECTOR(BMRVals, 1, 5);
	  return ;
	}			/* end if */
    }				/* end if */


  if (BMR - p[3] == 0.0)
    thedose = 0;
  if (BMR - p[3] == p[4]) {
    FREE_DVECTOR(pBak, 1, nparm);
    FREE_DVECTOR(BMRVals, 1, 5);
#ifdef MISC_OUT
    fprintf(stderr, "\n%s\n", "BMD computation failed.");
#endif
    fprintf(fp_out, "\n%s\n", "BMD computation failed.");
    return ;
  } else
    thedose = fabs (p[6] * power (p[4] / (BMR - p[3]) - 1, -1 / p[5]));

  *BMD = thedose;
  OUTPUT_BENCHMD2 (thedose);

  /* Otherwise, BMD is too large */
  /* else
     {
     fprintf(fp_out, "\n Benchmark dose at least 100 times the range of input data. \n Setting BMD to be 100 times the maximum dose.");
     printf("\n Benchmark dose at least 100 times the range of input data. \n Setting BMD to be 100 times the maximum dose.");
     *BMD = xb;
     OUTPUT_BENCHMD2(1, (*BMD) * scale);
     }   *//* end if (fa > 0) */

  if (thedose < 1e-30)
    {
      react = p[3];
#ifndef RBMDS
      fprintf (fp_out2, "\n  RSL \t%f", react);
      fprintf (fp_out2, "\n  BMD \t%f", thedose);
      fprintf (fp_out2, "\n\n BMD_line");
      fprintf (fp_out2, "\n %f %f", -1.0, react);
      fprintf (fp_out2, "\n %f %f", thedose, react);
      fprintf (fp_out2, "\n %f %f", thedose, -0.1);
      fprintf (fp_out2, "\n\n BMDL_comput_ind %d", No);
#endif
      FREE_DVECTOR(pBak, 1, nparm);
      FREE_DVECTOR(BMRVals, 1, 5);
#ifdef MISC_OUT
      fprintf(stderr, "\n%s\n", "BMDL computation failed.");
#endif
      fprintf(fp_out, "\n%s\n", "BMDL computation failed.");
      return ;
    }


  react = p[3] + p[4] / (pow (p[6] / thedose, p[5]) + 1);
#ifndef RBMDS
  fprintf (fp_out2, "\n  RSL \t%f", react);
  fprintf (fp_out2, "\n  BMD \t%f", thedose);
  fprintf (fp_out2, "\n\n BMD_line");
  fprintf (fp_out2, "\n %f %f", -1.0, react);
  fprintf (fp_out2, "\n %f %f", thedose, react);
  fprintf (fp_out2, "\n %f %f", thedose, -0.1);
#endif

  bDose = thedose;
  BMR = temp;

  /* search for BMDL */
  tol = FMAX ((*BMD) * 0.001, 0.0000001);
  // xa = thedose * 0.5;
  xb = thedose;
  BMD_lk = xlk;			/* get the lk at BMD */

  // fb = -LR;

  // BMDL_func returns the BMDL if successful, -1 otherwise
  fa = BMDL_func (nparm, BMD_lk, xb, p);
  if (fa < 0.0)
    {
#ifndef RBMDS
      fprintf (fp_out2, "\n\n BMDL_comput_ind %d", No);	/* computation failed */
#endif
      /* ERRORPRT ("BMDL computation failed."); */
      FREE_DVECTOR(pBak, 1, nparm);
      FREE_DVECTOR(BMRVals, 1, 5);
#ifdef MISC_OUT
      fprintf(stderr, "\n%s\n", "BMDL computation failed.");
#endif
      fprintf(fp_out, "\n%s\n", "BMDL computation failed.");
      return ;
    }
  /***8
      }
  ***/
#ifndef RBMDS
  fprintf (fp_out2, "\n\n BMDL_comput_ind %d", Yes);	/* computation will succeed */
#endif
  Bmdl[1] = fa;
#ifdef MISC_OUT
  printf ("           BMDL =%14.6g\n\n", (Bmdl[1]));
#endif
  fprintf (fp_out, 
#ifndef RBMDS
	   "            BMDL =%14.6g\n\n"
#else
	   "            BMDL =%30.22g\n\n"
#endif
	   , Bmdl[1]);

#ifndef RBMDS
  fprintf (fp_out2, "\n  BMDL \t%f", Bmdl[1]);
  fprintf (fp_out2, "\n\n BMDL_line");
  fprintf (fp_out2, "\n %f %f", Bmdl[1], -0.1);
  fprintf (fp_out2, "\n %f %f", Bmdl[1], react);
#endif


  tempdose[1] = thedose;

  for (j = 1; j <= nparm; j++)
    {
      p[j] = pBak[j];		/* get the "old" p[] */
    }

  HillMeans (1, p, tempdose, bmdlmean);

  Rlevel[1] = bmdlmean[1];

  /* bmdlCurve = No; */

  if (bmdlCurve == Yes)
    {

      Get_BMRS (p, bDose, Bmdl[1], BMRVals, sign, bmr_type);



      /* Now BMRVals[2]...BMRVals[5] contain Point BMRs for BMDL curve
         /  point calculations */

      for (k = 2; k <= 5; k++)
	{
	  /* solve the BMD */

	  for (j = 1; j <= nparm; j++)
	    {
	      p[j] = pBak[j];	/* Get back p[] */
	    }

	  xa = D = 0.0;
	  tol = 0.0000001;
	  incre = 1.0;

	  fa = 1;

	  /* Get a BMR value that is equal to just Mu(BMD) regardless
	     /  of the BMR type to simplify the root finding routine */

	  if ((bmr_type == 0) || (bmr_type == 1) || (bmr_type == 2))
	    {
	      if (sign == -1)
		{
		  bmrtemp = -BMRVals[k];
		}
	      else
		{
		  bmrtemp = BMRVals[k];
		}

	      if (bmr_type == 1)
		{
		  if (p[2] != 0)
		    {
		      bmrtemp = bmrtemp * sqrt (p[1] * pow (p[3], p[2]));
		    }
		  else
		    {
		      bmrtemp = bmrtemp * sqrt (p[1]);
		    }
		}
	      else if (bmr_type == 2)
		{
		  bmrtemp = bmrtemp * p[3];
		}
	    }
	  else
	    {
	      if (bmr_type == 3)
		{
		  bmrtemp = BMRVals[k] - p[3];
		}
	      else
		{
		  bmrtemp = p[4] * BMRVals[k];
		}
	    }			/* end if ((bmr_type == 0) || (bmr_type == 1) || (bmr_type == 2)) */

	  bmrtemp = bmrtemp + p[3];	/* bmrtype is now a Point BMR type
					   /  regardless of the user input
					   /  type */

	  BMR = bmrtemp;


	  if (fa > 0)
	    {
	      thedose =
		fabs (p[6] * power (p[4] / (BMR - p[3]) - 1, -1 / p[5]));
	      /*
		fa = BMD_func(nparm,p,xa,tol);
		fb = BMD_func(nparm,p,xb,tol);
		thedose = Binary_root(BMD_func, xa, xb, fa, fb, tol, nparm, p);
	      */
	    }

	  if (thedose > 100 * (*BMD))
	    {
#ifndef RBMDS
	      fprintf (fp_out2,
		       "\n\n BMDL_Curve_flag \t %d  \n smooth_opt  %d", 0,
		       smooth);
#endif
              FREE_DVECTOR(pBak, 1, nparm);
              FREE_DVECTOR(BMRVals, 1, 5);
#ifdef MISC_OUT
              fprintf(stderr, "\n%s\n", "BMDL computation failed for one or more point(s) on the BMDL curve.  \n          The BMDL curve will not be plotted\n");
#endif
              fprintf(fp_out, "\n%s\n","BMDL computation failed for one or more point(s) on the BMDL curve.  \n          The BMDL curve will not be plotted\n" );
              return ;
	    }

	  /* search for BMDL */
	  tol = FMAX ((*BMD) * 0.001, 0.0000001);
	  xa = thedose / 100.0;
	  xb = thedose;
	  BMD_lk = xlk;		/* get the lk at BMD */
	  fb = -LR;

	  BMR = BMRVals[k];

	  fa = BMDL_func (nparm, BMD_lk, xb, p);

	  Bmdl[k] = fa;
	  Rlevel[k] = bmrtemp;	// bmdlmean[1];

	  if (fa < 0.0)
	    {
	      /* for (j = 1; j <= 5; j++)
	         {
	         fprintf(fp_out2,"\n %f %f", -1.0, -1.0);
	         } */
#ifndef RBMDS
	      fprintf (fp_out2,
		       "\n\n BMDL_Curve_flag \t %d  \n smooth_opt  %d", 0,
		       smooth);
#endif
              FREE_DVECTOR(pBak, 1, nparm);
              FREE_DVECTOR(BMRVals, 1, 5);
#ifdef MISC_OUT
              fprintf(stderr, "\n%s\n", "BMDL computation failed for one or more point(s) on the BMDL curve.  \n          The BMDL curve will not be plotted\n");
#endif
              fprintf(fp_out, "\n%s\n","BMDL computation failed for one or more point(s) on the BMDL curve.  \n          The BMDL curve will not be plotted\n" );
              return ;
	    }


	  tempdose[1] = thedose;

	  for (j = 1; j <= nparm; j++)
	    {
	      p[j] = pBak[j];	/* get the "old" p[] */
	    }

	  HillMeans (1, p, tempdose, bmdlmean);

	}			/* end for (k = 2; k <= 5; k++) */

    }				/* end: if (bmdlCurve == Yes) */
  else
    {

      for (j = 2; j <= 5; j++)
	{
	  Bmdl[j] = Rlevel[j] = -1;
	}

    }				/* end else */
#ifndef RBMDS
  fprintf (fp_out2, "\n\n BMDL_Curve_flag \t %d  \n smooth_opt  %d",
	   bmdlCurve, smooth);
  fprintf (fp_out2, "\n\n BMDL_curve");
  fprintf (fp_out2, "\n 0.00000 %f", p[3]);
  for (k = 1; k <= 5; k++)
    {
      fprintf (fp_out2, "\n %f %f", Bmdl[k], Rlevel[k]);
    }
#endif

  for (j = 1; j <= nparm; j++)
    {
      p[j] = pBak[j];
    }

  FREE_DVECTOR (pBak, 1, nparm);
  FREE_DVECTOR(BMRVals, 1, 5);

}				/* end Hill_BMD */


/*****************************************************************
 *	BMD_func -- assumes the value 0 if Dose equals BMD.
 *
 ******************************************************************/
double
BMD_func (int nparm, float pBak[], float D, float gtol)
{
  double fd, fbmd, tD[2], m[2];
  double *p;
  int j;


  p = DVECTOR (1, nparm);

  for (j = 1; j <= nparm; j++)
    {
      p[j] = pBak[j];		/* get the "old" p[] */
    }

  tD[1] = D;

  HillMeans (1, p, tD, m);

  fd = m[1];

  fbmd = fd - BMR;

  FREE_DVECTOR (1, nparm);

  return fbmd;
}				/* end BMD_func */



/**********************************************************
 *	READ_OBSDATA4V--used to read 4 column data in 4 vectors.
 ***********************************************************/
int
READ_OBSDATA4V (int Nobs, double Xi[], int Ni[], double Ym[], double Yd[])
{
  int Nmiss;			/* number of records with missing values */
  int i, j, n, m;		/* count and iteration control variables */
  double value;			/* temp variable */

  Nmiss = 0;
  for (i = 1; i <= Nobs; i++)
    {
      n = i - Nmiss;
      m = 0;
      for (j = 1; j <= 4; j++)
	{
	  fscanf (fp_in, "%lf", &value);
	  if (value != MISSING)
	    {
	      if (j == 1)
		{
		  Xi[n] = value;
		}
	      if (j == 2)
		{
		  Ni[n] = (int) value;
		}
	      if (j == 3)
		{
		  Ym[n] = value;
		}
	      if (j == 4)
		{
		  Yd[n] = value;
		}
	    }
	  else
	    {
	      m++;
	    }
	}			/* end for (j = 1; j <= 4; j++) */

      if (m != 0)
	{
	  Nmiss++;
	}
      else if (Xi[n] < 0)
	{
	  Nmiss++;
	}
    }				/* end for (i = 1; i <= Nobs; i++) */
  return Nmiss;
}				/* end READ_OBSDATA4V */


/**********************************************************
 *	READ_OBSDATA2V--used to read 2 column data in 2 vectors.
 ***********************************************************/
int
READ_OBSDATA2V (int ntotal, double xxi[], double yyi[])
{
  int Nmiss;			/*number of records with missing values */
  int i, j, n, m;		/*count and iteration control variables */
  double value;			/*temp variable */

  Nmiss = 0;
  for (i = 1; i <= ntotal; i++)
    {
      n = i - Nmiss;
      m = 0;
      for (j = 1; j <= 2; j++)
	{
	  fscanf (fp_in, "%lf", &value);
	  if (value != MISSING)
	    {
	      if (j == 1)
		{
		  xxi[n] = value;
		}
	      if (j == 2)
		{
		  yyi[n] = value;
		}
	    }
	  else
	    {
	      m++;
	    }
	}

      if (m != 0)
	{
	  Nmiss++;
	}
      else if (xxi[n] < 0)
	{
	  Nmiss++;
	}
    }				/* end for (i = 1; i <= ntotal; i++) */

  return Nmiss;
}				/* end READ_OBSDATA2V */

/***************************************************
 *	OUTPUT_BENCHMD2--output specified benchmark dose.
 ****************************************************/
void
OUTPUT_BENCHMD2 (double BMD)
{

  OUTPUT_TEXT (" \n\n        Benchmark Dose Computation");

  /* output to the screen and to bmdswrk.002 temp file */
#ifdef MISC_OUT
  printf ("Specified effect =%14.6g\n\n", bmdparm.effect);
#endif
  fprintf (fp_out, "\nSpecified effect =%14.6g\n\n", bmdparm.effect);

  if (bmr_type == 0)
    {
#ifdef MISC_OUT
      printf ("Risk Type        =     Absolute deviation \n\n");
#endif
      fprintf (fp_out, "Risk Type        =     Absolute deviation \n\n");
    }
  else if (bmr_type == 1)
    {
#ifdef MISC_OUT
      printf
	("Risk Type        =     Estimated standard deviations from the control mean \n\n");
#endif
      fprintf (fp_out,
	       "Risk Type        =     Estimated standard deviations from the control mean \n\n");
    }
  else if (bmr_type == 2)
    {
#ifdef MISC_OUT
      printf ("Risk Type        =     Relative deviation \n\n");
#endif
      fprintf (fp_out, "Risk Type        =     Relative deviation \n\n");
    }
  else if (bmr_type == 3)
    {
#ifdef MISC_OUT
      printf ("Risk Type        =     Point estimate \n\n");
#endif
      fprintf (fp_out, "Risk Type        =     Point estimate \n\n");
    }
  else
    {
#ifdef MISC_OUT
      printf ("Risk Type        =     Extra risk \n\n");
#endif
      fprintf (fp_out, "Risk Type        =     Extra risk \n\n");
    }
#ifdef MISC_OUT
  printf ("Confidence level =%14.6g\n\n", bmdparm.level);
  printf ("             BMD =%14.6g\n\n", BMD);
#endif
  fprintf (fp_out, "Confidence level = %14.6g\n\n", bmdparm.level);
  fprintf (fp_out, 
#ifndef RBMDS
	   "             BMD = %14.6g\n\n"
#else
	   "             BMD = %30.22g\n\n"
#endif
	   , BMD);

}				/* end OUTPUT_BENCHMD2 */


/**********************************************************
 *	READ_OBSDATAV--used to read Ni column data in Ni vectors.
 ***********************************************************/
int
READ_OBSDATAV (int Nobs, double Xi[], int Ni[], double Ym[], double Yd[])
{
  int Nmiss;			/* number of records with missing values */
  int i, j, n, m;		/* count and iteration control variables */
  double value;			/* temp variable */
  double *yy;			/* response of each animal in a group */

  yy = DVECTOR (1, Nobs);
  Nmiss = 0;

  for (i = 1; i <= Nobs; i++)
    {
      n = i - Nmiss;
      m = 0;

      for (j = 1; j <= 2; j++)
	{
	  fscanf (fp_in, "%lf", &value);
	  if (value != MISSING)
	    {
	      if (j == 1)
		{
		  Xi[n] = value;
		}
	      if (j == 2)
		{
		  Ni[n] = (int) value;
		}
	    }
	  else
	    {
	      m++;
	    }
	}			/* end for (j = 1; j <= 2; j++) */

      for (j = 3; j <= Ni[n] + 2; j++)
	{
	  fscanf (fp_in, "%lf", &value);
	  if (value != MISSING)
	    {
	      yy[j - 2] = value;
	      Ym[j - 2] += value / Ni[n];
	    }
	  else
	    {
	      m++;
	    }
	}			/* end for (j = 3; j <= Ni[n] + 2; j++) */

      for (j = 3; j <= Ni[n] + 2; j++)
	{
	  fscanf (fp_in, "%lf", &value);
	  if (value != MISSING)
	    {
	      Yd[j - 2] +=
		(yy[j - 2] - Ym[j - 2]) * (yy[j - 2] - Ym[j - 2]) / (Ni[n] -
								     1.0);
	    }
	  else
	    {
	      m++;
	    }
	}			/* end for (j = 3; j <= Ni[n] + 2; j++) */

      if (m != 0)
	{
	  Nmiss++;
	}
      else if (Xi[n] < 0)
	{
	  Nmiss++;
	}

    }				/* end for (i = 1; i <= Nobs; i++) */

  return Nmiss;
}				/* end READ_OBSDATAV */

/*******************************************************************
 *	AThree_Fit fits a likelihood to a "full" model of the form
 *	Yij = Mu(i) + e(ij)  Var(eij) = k*Mu(i)^b.  The parameters Mu(i)
 *	i = 1,2,...,Nobs, k, and p are estimated using the likelihood
 *	maximization procedure, and the lkA3 = log-likelihood value
 *	upon return.
 *******************************************************************/

void
AThree_Fit (int nparm, double p[], double gtol, int *iter, double *fret)
{

  int i;
  double *parms, *fitparms, *doses, ll, *svar, *means;
  double /**tmy, *t, **tmv,*/ *bsv;
  double **X, **XP, **XPX, *XPY, *Y;
  long int nvar, *nanim, restr, nparms, *bind, nresm, optite;

  nvar = Nobs;
  restr = restrict;
  nparms = nparm;
  X = DMATRIX (1, Nobs, 1, 2);
  XP = DMATRIX (1, 2, 1, Nobs);
  XPX = DMATRIX (1, 2, 1, 2);
  XPY = DVECTOR (1, 2);
  bsv = DVECTOR (1, 2);
  Y = DVECTOR (1, Nobs);

    for (i = 1; i <= Nobs; i++)
    {
      X[i][1] = 1.0;
      X[i][2] = log (Ym[i]);
      Y[i] = log (Yd[i]);
    }

  TRANSPOSE (X, XP, Nobs, 2);

  MATMPYM2 (XP, X, XPX, 2, Nobs, 2);

  INVMAT (XPX, 2);

  MATMPYV2 (2, Nobs, XP, Y, XPY);

  MATMPYV2 (2, 2, XPX, XPY, bsv);

  p[1] = bsv[1];
  p[2] = bsv[2];
 
  for (i = 1; i <= nparm - 2; i++)
    {
      p[i + 2] = Ym[i];
    }

  doses = (double *) malloc ((size_t) (Nobs) * sizeof (double));
  means = (double *) malloc ((size_t) (Nobs) * sizeof (double));
  svar = (double *) malloc ((size_t) (Nobs) * sizeof (double));
  nanim = (long int *) malloc ((size_t) (Nobs) * sizeof (long int));
  parms = (double *) malloc ((size_t) (nparm) * sizeof (double));
  fitparms = (double *) malloc ((size_t) (nparm) * sizeof (double));
  bind = (long int *) malloc ((size_t) (nparm) * sizeof (long int));

  for (i = 1; i <= Nobs; i++)
    {
      nanim[i - 1] = (int) Ni[i];
      doses[i - 1] = Xi[i];
      means[i - 1] = Ym[i];
      svar[i - 1] = Yd[i];
    }

  for (i = 1; i <= nparm; i++)
    {
      parms[i - 1] = p[i];
    }

  getmlea3_ (&nparms, parms, fitparms, &ll, &optite, &nresm, bind);

  for (i = 1; i <= nparm; i++)
    {
      p[i] = fitparms[i - 1];
    }

  *fret = ll;

  FREE_DMATRIX (X, 1, Nobs, 1, 2);
  FREE_DMATRIX (XP, 1, 2, 1, Nobs);
  FREE_DMATRIX (XPX, 1, 2, 1, 2);
  FREE_DVECTOR (XPY, 1, 2);
  FREE_DVECTOR (bsv, 1, 2);
  FREE_DVECTOR (Y, 1, Nobs);
  free(doses);
  free(means);
  free(svar);
  free(nanim);
  free(parms);
  free(fitparms);
  free(bind);
}				/* end AThree_Fit */

/***********************************************************
 *	Given a vector of parameter values, and the number of
 *	parameters in that vector, this function will return three
 *	new parameter values to restart the optimization if a "bad"
 *	completion code is returned from GETCL(), using a uniform
 *	random number centered at p[i]
 ***********************************************************/
void
GetNewParms (double *p, int size)
{
  int i;


  /* Find parameters by randomly selecting new parameters in
     /  a uniform interval of p[i] +/- .005*p[i] */

  for (i = 0; i < size; i++)
    {
      if (Spec[i + 1] != 1)
	{
	  p[i] =
	    ((p[i] * .010) * (double) rand () / (double) RAND_MAX) +
	    p[i] * .995;
	}
    }

  /* If parameters are to be restricted, make sure restrictions
     /  are not violated */

  if (p[0] <= 0)
    {
      p[0] = .00000001;
    }
  if (p[5] <= 0)
    {
      p[5] = .00000001;
    }
  if ((restrict == 1) && (p[4] <= 1))
    {
      p[4] = 1.00000001;
    }
  if ((restrict == 0) && (p[4] < 0))
    {
      p[4] = 0;
    }

}				/* end GetNewParms */


/*****************************************************************
 * BMDL_func -- uses getcl_() to calculate the BMDL using an extra
 * nonlinear equality constraint in donlp2.  Returns the found
 * bmdl if successful, and -1 otherwise.
 *****************************************************************/
double
BMDL_func (int nparm, double xlk, double Dose, double pBak[])
{				/* BMD_lk and LR are calculated in Multistage_BMD() */

  double bmdl, target;

  double *doses, *means, *svar, *parms, *beginp, *fitparms;
  long int *nanim, *Spec2, *bind;
  double  maxdose;
  int i, j, ii;
  long int which,  nresm, optite, model_type;
  long int nvar, signs, nparms,  flag, bmrtype, restr;

  /* Hill model */

  model_type = 2;
  nvar = Nobs;
  signs = sign;
  nparms = nparm;
  bmrtype = bmr_type;
  restr = restrict;



  doses = (double *) malloc ((size_t) (Nobs) * sizeof (double));
  means = (double *) malloc ((size_t) (Nobs) * sizeof (double));
  svar = (double *) malloc ((size_t) (Nobs) * sizeof (double));
  nanim = (long int *) malloc ((size_t) (Nobs) * sizeof (long int));
  parms = (double *) malloc ((size_t) (nparm) * sizeof (double));
  beginp = (double *) malloc ((size_t) (nparm) * sizeof (double));
  fitparms = (double *) malloc ((size_t) (nparm) * sizeof (double));
  Spec2 = (long int *) malloc ((size_t) (nparm) * sizeof (long int));
  bind = (long int *) malloc ((size_t) (nparm) * sizeof (long int));



  for (i = 1; i <= Nobs; i++)
    {
      nanim[i - 1] = (long int) Ni[i];
      doses[i - 1] = Xi[i];
      means[i - 1] = Ym[i];
      svar[i - 1] = Yd[i];
    }

  which = 1;			/* Want an lower confidence limit */

  target = xlk - LR;		/* The value we want the likelihood
				   /  at the BMDL to match */

  for (j = 1; j <= nparm; j++)
    {
      parms[j - 1] = pBak[j];	/* get the "old" p[] */
      Spec2[j - 1] = Spec[j];
      beginp[j - 1] = pBak[j];
    }

  /* Set up and call subroutine that calculates BMDL */

  maxdose = doses[0];

  for (i = 1; i < Nobs; i++)
    {
      if (doses[i] > maxdose)
	{
	  maxdose = doses[i];
	}
    }

  for (i = 0; i < Nobs; i++)
    {
      doses[i] = doses[i] / maxdose;	/* Scale dose to be between 0 and 1 */
    }



  Dose = Dose / maxdose;

  parms[5] = parms[5] / maxdose;
  beginp[5] = beginp[5] / maxdose;	/* K parameter effected by scaling */

  flag = 0;			/* del0 = 2.0, tau0 = 2.0 */

  /* This is the first set of calls to getcl.  The parameters will be scaled */
  /* internally by donlp2. */
  getcl_ (&which, &nparms, &BMR, &Dose, &target, parms, &bmrtype,
	  &bmdl, fitparms, &optite, &nresm, bind, &flag);


  if ((optite < 0) || (optite > 3))
    {
      for (flag = 0; flag <= 1; flag++)
	{

	  for (ii = 0; ii < 15; ii++)
	    {
#ifdef MISC_OUT
	      printf ("\noptite = %ld", optite);
#endif
	      GetNewParms (beginp, nparm);	/* Get a new starting point */

	      /* Try again */

	      getcl_ (&which, &nparms, &BMR, &Dose, &target, beginp, &bmrtype,
		      &bmdl, fitparms, &optite, &nresm, bind, &flag);
	      /*  	      getcl_ (&which, &nvar, doses, means, nanim, svar, &nparms, &BMR, */
	      /*  		      &Dose, &target, beginp, Spec2, beginp, &bmrtype, &restr, */
	      /*  		      &bmdl, fitparms, &optite, &nresm, bind, &signs, */
	      /*  		      &model_type, &flag); */

	      /* if optite >= 0 and < 3, it is successful, and we can stop */

	      if ((optite >= 0) && (optite <= 3))
		{
		  flag = 2;
		  break;
		}
	    }			/* end for (ii = 0; ii < 15; ii++) */

	  for (j = 1; j <= nparm; j++)
	    beginp[j - 1] = pBak[j];
	  beginp[5] = beginp[5] / maxdose;

	}			/* end for (flag = 0; flag<=1; flag++) */
    }

  flag = 0;
  if ((optite < 0) || (optite > 3))
    {
      for (j = 1; j <= nparm; j++)
	beginp[j - 1] = pBak[j];
      beginp[5] = beginp[5] / maxdose;

      for (ii = 0; ii < 5; ii++)
	{
#ifdef MISC_OUT
	  printf ("\n   optite = %ld", optite);
#endif
	  GetCLParms (beginp, nparm);

	  getcl_ (&which, &nparms, &BMR, &Dose, &target, beginp, &bmrtype,
		  &bmdl, fitparms, &optite, &nresm, bind, &flag);
	  /*  	  getcl_ (&which, &nvar, doses, means, nanim, svar, &nparms, &BMR, */
	  /*  		  &Dose, &target, beginp, Spec2, beginp, &bmrtype, &restr, */
	  /*  		  &bmdl, fitparms, &optite, &nresm, bind, &signs, &model_type, */
	  /*  		  &flag); */

	  if ((optite >= 0) && (optite <= 3))
	    {
	      break;
	    }
	}
    }

  /* This is the second set of calls to getcl.  The parameters will not be scaled */
  /* internally by donlp2. */
  if ((optite < 0) || (optite > 3))
    {
      getcl_ (&which, &nparms, &BMR, &Dose, &target, parms, &bmrtype,
	      &bmdl, fitparms, &optite, &nresm, bind, &flag);
      /*    getcl_ (&which, &nvar, doses, means, nanim, svar, &nparms, &BMR, */
      /*  	  &Dose, &target, parms, Spec2, parms, &bmrtype, &restr, */
      /*  	  &bmdl, fitparms, &optite, &nresm, bind, &signs, &model_type, &flag); */


      /* optite is a value that is passed back from GETCL which
	 /  determines whether the optimization was completed successfully
	 /  If optite is less than 0, then it did not, and we want
	 /  to try a different starting point and recompute */

      if ((optite < 0) || (optite > 3))
	{
	  for (flag = 0; flag <= 1; flag++)
	    {

	      for (ii = 0; ii < 15; ii++)
		{
#ifdef MISC_OUT
		  printf ("\noptite = %ld", optite);
#endif

		  GetNewParms (beginp, nparm);	/* Get a new starting point */


		  /* Try again */

		  getcl_ (&which, &nparms, &BMR, &Dose, &target, beginp, &bmrtype,
			  &bmdl, fitparms, &optite, &nresm, bind, &flag);
		  /*  	      getcl_ (&which, &nvar, doses, means, nanim, svar, &nparms, &BMR, */
		  /*  		      &Dose, &target, beginp, Spec2, beginp, &bmrtype, &restr, */
		  /*  		      &bmdl, fitparms, &optite, &nresm, bind, &signs, */
		  /*  		      &model_type, &flag); */

		  /* if optite >= 0 and < 3, it is successful, and we can stop */

		  if ((optite >= 0) && (optite <= 3))
		    {
		      flag = 2;
		      break;
		    }
		}			/* end for (ii = 0; ii < 15; ii++) */

	      for (j = 1; j <= nparm; j++)
		beginp[j - 1] = pBak[j];
	      beginp[5] = beginp[5] / maxdose;

	    }			/* end for (flag = 0; flag<=1; flag++) */
	}

      flag = 0;
      if ((optite < 0) || (optite > 3))
	{
	  for (j = 1; j <= nparm; j++)
	    beginp[j - 1] = pBak[j];
	  beginp[5] = beginp[5] / maxdose;

	  for (ii = 0; ii < 5; ii++)
	    {
#ifdef MISC_OUT
	      printf ("\n   optite = %ld", optite);
#endif
	      GetCLParms (beginp, nparm);


	      getcl_ (&which, &nparms, &BMR, &Dose, &target, beginp, &bmrtype,
		      &bmdl, fitparms, &optite, &nresm, bind, &flag);
	      /*  	  getcl_ (&which, &nvar, doses, means, nanim, svar, &nparms, &BMR, */
	      /*  		  &Dose, &target, beginp, Spec2, beginp, &bmrtype, &restr, */
	      /*  		  &bmdl, fitparms, &optite, &nresm, bind, &signs, &model_type, */
	      /*  		  &flag); */



	      if ((optite >= 0) && (optite <= 3))
		{
		  break;
		}
	    }
	}
    }


  /* Let user know if no optimum was found */

  if ((optite < 0) || (optite > 3))
    {
      //      fprintf (fp_out,
      //       "Warning:  optimum may not have been found.  Bad completion code in Optimization routine.\n");
      // fprintf (fp_out, "\n");
      free(bind);
      free(doses);
      free(svar);
      free(means);
      free(nanim);
      free(parms);
      free(beginp);
      free(fitparms);
      free(Spec2);
      return -1;
    }				/* end if ((optite < 0) || (optite >= 3)) */


  for (j = 1; j <= nparm; j++)
    {
      pBak[j] = fitparms[j - 1];
    }

#ifdef MISC_OUT
  printf ("\n optite = %ld", optite);
#endif
  /* rescale dose back to normal */

  pBak[6] = maxdose * pBak[6];

  Dose = Dose * maxdose;


  bmdl = bmdl * maxdose;

  /*
   * Added to free the memory allocated locally
   * Qun He 12/29/03 {QH 2003/12/31 PR# }
   */

  free(bind);
  free(doses);
  free(svar);
  free(means);
  free(nanim);
  free(parms);
  free(beginp);
  free(fitparms);
  free(Spec2);
  /* end of modification */


  return bmdl;			/* return bmdl to the calling function */
}				/* end BMDL_func */



/****************************************************
 *	Given the number of dose levels, this function
 *	calculates the mean at each dose level and stores
 *	it in the array means[]
 *****************************************************/
void
HillMeans (int nobs, double p[], double Doses[], double means[])
{
  int i;
  double mi;

  for (i = 1; i <= nobs; i++)
    {
      if (Doses[i] == 0)
	{
	  mi = 0;
	}
      else
	{
	  if (p[5] == 0)
	    {
	      mi = 0.5 * p[4];
	    }
	  else
	    {
	      mi = p[4] / (pow (p[6] / Doses[i], p[5]) + 1);
	    }
	}

      means[i] = mi + p[3];
    }				/* end for (i = 1; i <= nobs; i++) */
}				/* end HillMeans */



/***********************************************************
 *	Calculates the second partial derivatives os the function
 *	F = ln(Vi) at dose[obs] where Vi is the estimated variance.  Fn1i[j][k]
 *	contains the second partial with respect to parameters j and
 *	k.
 ***********************************************************/
void
F1iDoublePart (int nparm, int const_var, double p[], double **Fn1i, int obs)
{
  double *mg, *vg, **mg2, **vg2, Vi;
  double meani, temp;
  int j, k;


  for (j = 1; j <= nparm; j++)
    {
      for (k = 1; k <= nparm; k++)
	{
	  Fn1i[j][k] = 0.0;
	}
    }

  mg = DVECTOR (1, nparm);
  vg = DVECTOR (1, nparm);
  mg2 = DMATRIX (1, nparm, 1, nparm);
  vg2 = DMATRIX (1, nparm, 1, nparm);


  /* Compute the estimated mean */

  if (Xi[obs] != 0)
    {
      meani = p[3] + p[4] / (pow (p[6] / Xi[obs], p[5]) + 1);
    }
  else
    {
      meani = p[3];
    }

  if (const_var == 1)
    {
      Vi = p[1];
    }
  else
    {
      Vi = exp(p[1] + log(fabs (meani)) * p[2]);
    }

  /* Get the partial derivatives of the mean function at dose[obs] */
  MeanPart (obs, p, mg);
  /* Get the partial derivatives of the variance function */
  VarPart (obs, const_var, Vi, meani, p, mg, vg);
  /* Get second partials of the mean function */
  Mean2Part (obs, p, mg2);
  /* Get second partials of the variance function */
  Var2Part (obs, const_var, Vi, meani, p, mg, mg2, vg2);


  /* Calculate partial derivative at dose[obs] */

  for (j = 1; j <= nparm; j++)
    {
      for (k = j; k <= nparm; k++)
	{
	  temp = vg2[j][k] - (vg[j] * vg[k] / Vi);
	  Fn1i[j][k] = temp / Vi;
	  Fn1i[k][j] = Fn1i[j][k];
	}
    }


  FREE_DVECTOR (mg, 1, nparm);
  FREE_DVECTOR (vg, 1, nparm);
  FREE_DMATRIX (mg2, 1, nparm, 1, nparm);
  FREE_DMATRIX (vg2, 1, nparm, 1, nparm);

}				/* end F1iDoublePart */



/***********************************************************
 *	Calculates the second partial derivatives of the function
 *	F = 1/Vi at dose[obs] where Vi is the estimated variance.  Fn2i[j][k]
 *	contains the second partial with respect to parameters j and
 *	k.
 ***********************************************************/
void
F2iDoublePart (int nparm, int const_var, double p[], double **Fn2i, int obs)
{
  double *mg, *vg, **mg2, **vg2, Vi;
  double meani, temp;
  int j, k;


  for (j = 1; j <= nparm; j++)
    {
      for (k = 1; k <= nparm; k++)
	{
	  Fn2i[j][k] = 0.0;
	}
    }

  mg = DVECTOR (1, nparm);
  vg = DVECTOR (1, nparm);
  mg2 = DMATRIX (1, nparm, 1, nparm);
  vg2 = DMATRIX (1, nparm, 1, nparm);


  /* Compute the estimated mean */

  if (Xi[obs] != 0)
    {
      meani = p[3] + p[4] / (pow (p[6] / Xi[obs], p[5]) + 1);
    }
  else
    {
      meani = p[3];
    }

  if (const_var == 1)
    {
      Vi = p[1];
    }
  else
    {
      Vi = exp(p[1] + log(fabs (meani)) * p[2]);
    }

  /* Get the partial derivatives of the mean function at dose[obs] */
  MeanPart (obs, p, mg);
  /* Get the partial derivatives of the variance function */
  VarPart (obs, const_var, Vi, meani, p, mg, vg);
  /* Get second partials of the mean function */
  Mean2Part (obs, p, mg2);
  /* Get second partials of the variance function */
  Var2Part (obs, const_var, Vi, meani, p, mg, mg2, vg2);


  /* Calculate partial derivative at dose[obs] */

  for (j = 1; j <= nparm; j++)
    {
      for (k = j; k <= nparm; k++)
	{
	  temp = (2 * vg[j] * vg[k] / Vi) - vg2[j][k];
	  Fn2i[j][k] = temp / (Vi * Vi);
	  Fn2i[k][j] = Fn2i[j][k];
	}
    }


  FREE_DVECTOR (mg, 1, nparm);
  FREE_DVECTOR (vg, 1, nparm);
  FREE_DMATRIX (mg2, 1, nparm, 1, nparm);
  FREE_DMATRIX (vg2, 1, nparm, 1, nparm);

}				/* end F2iDoublePart */



/***********************************************************
 *	Calculates the second partial derivatives of the function
 *	F = (Ybar - Mi)**2/Vi at dose[obs] where Vi is the estimated variance
 *	Ybar is the sample mean, and Mi is the estimated mean.
 *	Fn1i[j][k]
 *	contains the second partial with respect to parameters j and
 *	k.
 ***********************************************************/
void
F3iDoublePart (int nparm, int const_var, double p[], double **Fn3i, int obs)
{
  double *mg, *vg, **mg2, **vg2, Vi;
  double Devi, meani, temp, temp2, temp3;
  int j, k;


  for (j = 1; j <= nparm; j++)
    {
      for (k = 1; k <= nparm; k++)
	{
	  Fn3i[j][k] = 0.0;
	}
    }

  mg = DVECTOR (1, nparm);
  vg = DVECTOR (1, nparm);
  mg2 = DMATRIX (1, nparm, 1, nparm);
  vg2 = DMATRIX (1, nparm, 1, nparm);


  /* Compute the estimated mean */

  if (Xi[obs] != 0)
    {
      meani = p[3] + (p[4] / (pow (p[6] / Xi[obs], p[5]) + 1));
    }
  else
    {
      meani = p[3];
    }

  Devi = Ym[obs] - meani;

  if (const_var == 1)
    {
      Vi = p[1];
    }
  else
    {
      Vi = exp(p[1] + log(fabs (meani)) * p[2]);
    }

  /* Get the partial derivatives of the mean function at dose[obs] */
  MeanPart (obs, p, mg);
  /* Get the partial derivatives of the variance function */
  VarPart (obs, const_var, Vi, meani, p, mg, vg);
  /* Get second partials of the mean function */
  Mean2Part (obs, p, mg2);
  /* Get second partials of the variance function */
  Var2Part (obs, const_var, Vi, meani, p, mg, mg2, vg2);


  /* Calculate partial derivative at dose[obs] */


  for (j = 1; j <= nparm; j++)
    {
      for (k = j; k <= nparm; k++)
	{
	  temp = 2 * Vi * Vi * (mg[j] * mg[k] - Devi * mg2[j][k]);
	  temp2 = 2 * Devi * Vi * (vg[j] * mg[k] + mg[j] * vg[k]);
	  temp3 = Devi * Devi * (2 * vg[j] * vg[k] - Vi * vg2[j][k]);
	  Fn3i[j][k] = (temp + temp2 + temp3) / (Vi * Vi * Vi);
	  Fn3i[k][j] = Fn3i[j][k];
	}
    }


  FREE_DVECTOR (mg, 1, nparm);
  FREE_DVECTOR (vg, 1, nparm);
  FREE_DMATRIX (mg2, 1, nparm, 1, nparm);
  FREE_DMATRIX (vg2, 1, nparm, 1, nparm);
}



/*******************************************************************
 *	First partial derivatives of the mean function at observation obs,
 *	contained in mg[] at exit.
 *******************************************************************/
void
MeanPart (int obs, double *p, double *mg)
{
  double temp, temp2, temp3, p6;

  mg[1] = mg[2] = 0.0;		/* Variance parameters not in mean function */
  mg[3] = 1.0;

  if (p[6] < 1e-8)
    p6 = 1e-8;
  else
    p6 = p[6];
  if (Xi[obs] != 0)
    {
      temp = Xi[obs] / p6;
      temp2 = p6 / Xi[obs];
      temp3 = pow (Xi[obs], p[5]) + pow (p6, p[5]);

      if (temp3 < 1e-20)
	temp3 = 1e-20;
      if (temp < 1)
	{
	  mg[4] = pow (temp, p[5]) / (pow (temp, p[5]) + 1);
	}
      else
	{
	  mg[4] = 1 / (pow (temp2, p[5]) + 1);
	}

      mg[5] = p[4] * pow (Xi[obs] * p6, p[5]) * (Slog (Xi[obs]) - Slog (p6));
      mg[5] = mg[5] / (temp3 * temp3);
      mg[6] =
	-p[4] * p[5] * (pow (Xi[obs] * p6, p[5])) / (p6 * temp3 * temp3);
    }
  else
    {
      mg[4] = mg[5] = mg[6] = 0.0;
    }

}				/* end MeanPart */



/*********************************************************************
 *	First partial derivatives of the variance function.  Vi and meani
 *	are the estimated mean and variance respectively, and should be
 *	passed by the calling function.  const_var = 1 if variance is
 *	constant, = 0 if not.  mg[] are the first partials of the mean
 *	function and should be passed by the calling function.  Partials
 *	stored in vg[] upon exit.
 *********************************************************************/
void
VarPart (int obs, int const_var, double Vi, double meani, double *p,
	 double *mg, double *vg)
{
  int j;

  if (const_var == 1)
    {
      vg[1] = Vi / p[1];
      vg[2] = 0.0;
    }
  else
    {
      vg[1] = Vi;
      vg[2] = Vi * Slog (fabs (meani));
    }

  for (j = 3; j <= nparm; j++)
    {
      if (fabs (meani) > 1e-20)
	vg[j] = p[2] * Vi * mg[j] / fabs (meani);
      else
	vg[j] = 0.0;

    }

}				/* end VarPart */



/********************************************************************
 *	Second partial derivates with respect to the mean function.
 *	mg2[][] contains all second partials upon exit.
 ********************************************************************/
void
Mean2Part (int obs, double *p, double **mg2)
{
  int j;
  double Ri, Si, Wi, gami, temp, p6;

  Ri = 0;
  Si = 0;
  Wi = 0;
  gami = 0;
  temp = 0;

  if (p[6] < 1e-8)
    p6 = 1e-8;
  else
    p6 = p[6];

  if (Xi[obs] != 0)
    {
      Ri = pow (Xi[obs] * p6, p[5]);
      Si = pow (Xi[obs], p[5]) + pow (p6, p[5]);
      Wi = pow (Xi[obs], p[5]) - pow (p6, p[5]);
      gami = Slog (Xi[obs]) - Slog (p6);
      temp = p[4] * Ri / (Si * Si * Si);
    }

  /* Second partials of the mean function involving
     /  the variance and background parameters */

  for (j = 1; j <= nparm; j++)
    {
      mg2[1][j] = mg2[j][1] = mg2[2][j] = mg2[j][2] = mg2[3][j] = mg2[j][3] =
	0.0;
    }

  mg2[4][4] = 0.0;

  if (Xi[obs] != 0)
    {
      mg2[4][5] = Ri * gami / (Si * Si);
      mg2[4][6] = -p[5] * Ri / (p6 * Si * Si);
      mg2[5][5] = -temp * Wi * gami * gami;
      mg2[5][6] = temp * ((p[5] * gami * Wi) - Si) / p6;
      mg2[6][6] = temp * p[5] * (Si - (p[5] * Wi)) / (p6 * p6);
    }
  else
    {
      mg2[4][5] = mg2[4][6] = mg2[5][5] = mg2[5][6] = mg2[6][6] = 0.0;
    }

  mg2[5][4] = mg2[4][5];
  mg2[6][4] = mg2[4][6];
  mg2[6][5] = mg2[5][6];

}				/* end Mean2Part */



/*********************************************************************
 *	Second partial derivatives of the variance function.  Vi and meani
 *	are the estimated mean and variance respectively, and should be
 *	passed by the calling function.  const_var = 1 if variance is
 *	constant, = 0 if not.  mg[] are the first partials of the mean
 *	function, mg2[][] are the second partials of the mean function.
 *	Both should be passed by the calling function.  Partials
 *	stored in vg2[][] upon exit.
 *********************************************************************/
void
Var2Part (int obs, int const_var, double Vi, double meani, double *p,
	  double *mg, double **mg2, double **vg2)
{
  double logam, abmn, temp;
  int j, k, Sign;

  abmn = fabs (meani);
  logam = Slog (abmn);


  if (const_var == 1)
    {
      /* constant variance.  Vi = alpha.  All second partials = 0 */
      for (j = 1; j <= nparm; j++)
	{
	  for (k = 1; k <= nparm; k++)
	    {
	      vg2[j][k] = 0.0;
	    }
	}
    }
  else
    {
      /* non constant variance.  Vi = alpha*|Meani|**rho */

      if (meani < 0)
	{
	  Sign = -1;
	}
      else
	{
	  Sign = 1;
	}

      vg2[1][1] = Vi;
      vg2[1][2] = Vi * logam;
      vg2[2][1] = vg2[1][2];

      for (j = 3; j <= nparm; j++)
	{
	  vg2[1][j] = Sign * p[2] * Vi * mg[j] / abmn;
	  vg2[j][1] = vg2[1][j];
	}
      vg2[2][2] = Vi * logam * logam;

      for (j = 3; j <= nparm; j++)
	{
	  vg2[2][j] = Sign * Vi * mg[j] * ((p[2] * logam) + 1) / abmn;
	  vg2[j][2] = vg2[2][j];
	}

      for (j = 3; j <= nparm; j++)
	{
	  for (k = j; k <= nparm; k++)
	    {
	      temp = ((p[2] - 1) * mg[j] * mg[k] / abmn) + (Sign * mg2[j][k]);
	      vg2[j][k] = p[2] * Vi * temp / abmn;
	      vg2[k][j] = vg2[j][k];
	    }
	}
    }				/* end if (const_var == 1) */

}				/* end Var2Part */



/*********************************************************************
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
void
Get_BMRS (double *p, double Dose, double bmdl1, double *BMRVals, int sign,
	  int bmr_type)
{

  double range, mn;
  double ymax, ymin, bmr_percentage;
  double Percent_Change[4] = { .1, .225, .375, .5 };
  int i;

  ymax = Ym[1];
  ymin = Ym[1];

  /* Find the maximum and minimum observed mean to get BMR's in terms
     /  of a percentage of the range of observed data */

  for (i = 2; i <= Nobs; i++)
    {
      if (Ym[i] > ymax)
	{
	  ymax = Ym[i];
	}
      if (Ym[i] < ymin)
	{
	  ymin = Ym[i];
	}
    }

  range = ymax - ymin;		/* The range of means */

  bmr_percentage = (BMRVals[1] - ymin) / range;	/* bmrtype as a percentage
						   /  of the range */

  /* If adverse direction is "down", then make this as a percentage
     /  of the range, but the pecent decrease over the maximum observed
     /  mean as opposed to the percent increase over the minimum */

  if (sign == -1)
    {
      bmr_percentage = 1 - bmr_percentage;
    }

  if (sign == 1)
    {
      mn = ymin;
    }
  else
    {
      mn = ymax;
    }

  /* Find a set of  Point type BMRs that are spaced throughout
     /  the range, and such that they are not too crowded by the already
     /  computed BMDL */

  if ((.05 < bmr_percentage) && (bmr_percentage < .10))
    {
      Percent_Change[0] = bmr_percentage + .05;
    }
  else if ((.2 < bmr_percentage) && (bmr_percentage < .25))
    {
      Percent_Change[1] = bmr_percentage + .05;
    }
  else if ((.35 < bmr_percentage) && (bmr_percentage < .4))
    {
      Percent_Change[2] = bmr_percentage + .05;
    }
  else if ((.475 < bmr_percentage) && (bmr_percentage < .525))
    {
      Percent_Change[3] = bmr_percentage + .05;
    }


  /* Get BMR values and give them back as the BMR type specified by
     /  the user */

  for (i = 0; i <= 3; i++)
    {
      BMRVals[i + 2] = mn + (sign * Percent_Change[i] * range);

      if (bmr_type == 0)	/* Absolute deviation */
	{
	  BMRVals[i + 2] = fabs (BMRVals[i + 2] - p[3]);
	}
      else if (bmr_type == 1)	/* Std Dev */
	{
	  BMRVals[i + 2] = fabs (BMRVals[i + 2] - p[3]);
	  BMRVals[i + 2] = BMRVals[i + 2] / (sqrt (exp(p[1] + log(p[3])* p[2])));
	}
      else if (bmr_type == 2)	/* Relative deviation */
	{
	  BMRVals[i + 2] = fabs (BMRVals[i + 2] - p[3]) / p[3];
	}
      else if (bmr_type == 4)	/* Extra risk */
	{
	  BMRVals[i + 2] = (BMRVals[i + 2] - p[3]) / p[4];
	}
      else
	/* Point estimate */
	{
	  BMRVals[i + 2] = BMRVals[i + 2];
	}
    }/* end for (i = 0; i <= 3; i++) */

  return;
}				/* end Get_BMRS */

/********************************************************
 * power function to safely compute a^b
 *********************************************************/
double
power (double a, double b)
{
  double log_a, blog_a, a_b;

  if (a >= 1.0e-10)
    log_a = log (a);
  else
    log_a = -24.5258509299404572 + (2.0e10) * a - (5.0e19) * a * a;

  blog_a = b * log_a;
  if (blog_a > 700)
    blog_a = 700.0;

  a_b = exp (blog_a);
  return a_b;
}				/* end power */

/***********************************************************
 *	Given a vector of parameter values, and the number of
 *	parameters in that vector, this function will return three
 *	new parameter values to restart the optimization if a "bad"
 *	completion code is returned from GETCL(), using a uniform
 *	random number centered at p[i]
 ***********************************************************/
void
GetOtherParms (double *p, int size)
{
  int i;


  /* Find parameters by randomly selecting new parameters in
     /  a uniform interval of  */

  for (i = 0; i < size; i++)
    {
      if (Spec[i + 1] != 1)
	{

	  p[i] = 8 * (double) rand () / (double) RAND_MAX;
	}
    }

  /* If parameters are to be restricted, make sure restrictions
     /  are not violated */

  if (p[0] <= 0)
    {
      p[0] = -p[0];
    }
  if (p[5] <= 0)
    {
      p[5] = -p[5];
    }
  if ((restrict == 1) && (p[4] <= 1))
    {
      p[4] = 1 + p[4];
    }
  if ((restrict == 0) && (p[4] < 0))
    {
      p[4] = -p[4];
    }

}				/* end GetOtherParms */

/***********************************************************
 *	Given a vector of parameter values, and the number of
 *	parameters in that vector, this function will return three
 *	new parameter values to restart the optimization if a "bad"
 *	completion code is returned from GETCL(), using a uniform
 *	random number centered at p[i]
 ***********************************************************/
void
GetMoreParms (double *p, int size)
{
  int i;


  /* Find parameters by randomly selecting new parameters in
     /  a uniform interval of p[i] +/-  */

  for (i = 0; i < size; i++)
    {
      if (Spec[i + 1] != 1)
	{
	  p[i] = p[i] + .2 * p[i] * (double) rand () / (double) RAND_MAX;
	}
    }

  /* If parameters are to be restricted, make sure restrictions
     /  are not violated */

  if (p[0] <= 0)
    {
      p[0] = -p[0];
    }
  if (p[5] <= 0)
    {
      p[5] = -p[5];
    }
  if ((restrict == 1) && (p[4] <= 1))
    {
      p[4] = 1 + p[4];
    }
  if ((restrict == 0) && (p[4] < 0))
    {
      p[4] = -p[4];
    }

}				/* end GetMoreParms */

/***********************************************************
 *	Given a vector of parameter values, and the number of
 *	parameters in that vector, this function will return three
 *	new parameter values to restart the optimization if a "bad"
 *	completion code is returned from GETCL(), using a uniform
 *	random number centered at p[i]
 ***********************************************************/
void
GetCLParms (double *p, int size)
{
  int i;


  /* Find parameters by randomly selecting new parameters in
     /  a uniform interval of p[i] +/-  */

  for (i = 0; i < size; i++)
    {
      if (Spec[i + 1] != 1)
	{
	  p[i] = p[i] + .1 * p[i] * (double) rand () / (double) RAND_MAX;
	}
    }

  /* If parameters are to be restricted, make sure restrictions
     /  are not violated */

  if (p[0] <= 0)
    {
      p[0] = -p[0];
    }
  if (p[5] <= 0)
    {
      p[5] = -p[5];
    }
  if ((restrict == 1) && (p[4] <= 1))
    {
      p[4] = 1 + p[4];
    }
  if ((restrict == 0) && (p[4] < 0))
    {
      p[4] = -p[4];
    }

}				/* end GetCLParms */


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
