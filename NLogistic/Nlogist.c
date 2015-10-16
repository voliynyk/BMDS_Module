/****************************************************************
*
* IMPORTANT NOTE:  The following variable is the version number for
*                  the current model.  THIS MUST BE CHANGED as
*		   important changes are made to the models.
*
*****************************************************************/
char Version_no[]="NLogistic Model. (Version: 2.20; Date: 04/27/2015)";

/****************************************************************
*
* Nlogist.C - a ANSI C program for Nlogist model fitting with/without
*             a natural background rate in Benchmark Dose.
*
* Date: March 25, 1997
*
********************************************************************
* Modification Log:
*
* Version Number: 2.9
* Modified By: Micheal Ferree
* Modified Date: 6/21/2005
* Reason: Took out analysis of deviance table
*
* Version Number: 2.10
* Modified By: Woodrow Setzer
* Modified Date: 10/27/2005
* Reason:
*   1) Free allocated memory
*   2) Fix illegal memory reads and writes
*   3) Report only fitted model log-likelihood and AIC (and not "reduced"
*      and "full" models, which were incorrect.
*   4) Added conditional compilation flags for RBMDS and logging.
*
* Version Number: 2.11
* Modified By: Woodrow Setzer
* Modified Date: 03/22/2006
* Reason: Replaced calls to MATINVS with INVMAT
* Version Number: 2.6
*
* Version Number: 2.12
* Modified By: Geoffrey
* Date: 1/12/2007
* Reason: Incremented version number.
*	  Added last parameter "0" (don't print SE) in OUTPUT_DTMS3PARMS().
*
* Version Number: 2.13
* Modified By: Woodrow Setzer
* Date: 2/20/2007
* Reason: Incremented version number to reflect changed compilation options.
*
* Version Number: 2.14
* Modified By: G. Nonato
* Modification Date: 04/07/2008
* Reason: (Per BMDS 2.0: Problem Report 157 & 147)
*       Fix the Observation # < parameter # for NLogistic model problem.
*       Added code to free-up allocated memories before exiting thru ERRORPRT()
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
*
* Version Number: 2.17
* Modified By: Cody Simmons
* Modification Date: 09/08/2014
* Reason: PR232 and PR244
*       Added the bootstrap method for goodness of fit calculation
*
* Version Number: 2.18
* Modified By: Cody Simmons
* Modification Date: 09/15/2014
* Reason: PR300
*       Created scaled residual of interest table for nested models.
*
* Version Number: 2.19
* Modified By: Cody Simmons
* Modification Date: 11/12/2014
* Reason: PR524
*       Fixed calculation of VCV matrix to give correct SE estimates
*
* Version Number: 2.20
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

void Nlogist_fit (int n, int ngrp, double p[], double gtol,
		  int *iter, double *fret);
void Nlogist_BMD (int nparm, double p[], double gtol, int *iter, double xlk,
		  double Rlevel[], double Bmdl[], double *BMD);
int Nlogist_vcv (int nparm, int Spec[], double p[], double **vcv);
void Predict(double [], double [], int, double [], double []);
int Model_DF (int []);
void Which_Bounded (int [], double [], int []);
double BMDL_func (int nparm, double p[], double x, double gtol);
void TEMP_ANOVA_OUTPUT(char *anatxt[], AnaList anasum[]);


#define EPS 3.0e-8
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
/*** Define input and output files's name  *********************/
char fin[FLENGTH];			/*input temp file */
char fout[FLENGTH];			/*output temp file */
char fout2[FLENGTH];

char *Parm_name[] =
{"alpha", "beta", "theta1", "theta2", "rho",
 "phi1", "phi2", "phi3", "phi4", "phi5",
 "phi6", "phi7", "phi8", "phi9", "phi10"};
typedef enum {
  dummy,alpha, beta, theta1, theta2, rho, phi1, phi2, phi3,
  phi4, phi5, phi6, phi7, phi8, phi9, phi10} Model_Parms;
char *anatxt[] =
{"Full model", "Fitted model", "Reduced model"};


/*** variables will not be changed except Spec  *******/
int *Spec;			/*vector used to identify user input parm. */
int *IniSp;
double *Yp;			/*positive dependent variable data array */
double *Yn;			/*negative dependent variable data array */
double *Xi;			/*independent variable data array */
double *Ypp;			/*predicted positive depdendent values */
double *Ep;			/*estimated probability */
double *Ls;
int *Xg;
double *SR; 			/*scaled residual */
double *Rlevel;
double *Bmdl;
double *IniP;
int Nobs, nparm, ngrp, restrict, initial, appendix, smooth, fixedSize,
  bmdlCurve;
double xmax, xmin, scale;
double Min_increment, Max_double;
double Rel_Conv, Parm_Conv, Maxloglik;
double SlopeUpperBound = 18.0;

double BMDL_Error_Size;
int BMDL_Error;
int DeBuG = 0;
FILE *fp_dbg;

int BSIter;      /*number of iterations for bootstrapping method */     
long BSSeed;     /*seed values for bootstrapping method - 0 defaults to time-generated*/

double smean, smax, smin, sijfixed, spfixed, snfixed, sdif;

/*smean: mean of litter-specific covariate (lsc).
   smax: maximum lsc - smean.
   smin: minimun lsc - smean.

   also used in O_func, _g, _vcv. */

double smean1, xmax, gamm0;	// mean litter size in group

/** changing variable **/
int replace, brat;
double tD, BMD_lk, LR, ck, upb = 18, BMR;

int ErrorFlag;			/* Error States from DMNGB */

/* GLN - 04/07/2008
*  Free-up allocated memory before exit
*  upon encountering fatal error.
*/
//void FreeUp_mem(Parms, varsum, anasum, vcv, GXi, GYp, GYn, bounded);
void FreeUp_mem(double *Parms, VarList *varsum, AnaList *anasum, double  **vcv, double *GXi, double *GYp, double *GYn, int *bounded)
{
	FREE_DVECTOR (Parms, 1, nparm);
	FREE_DVECTOR (IniP, 1, nparm);
	FREE_DVECTOR (Xi, 1, Nobs);
	FREE_DVECTOR (Yp, 1, Nobs);
	FREE_DVECTOR (Yn, 1, Nobs);
	FREE_IVECTOR (Ls, 1, Nobs);
	FREE_IVECTOR (Xg, 1, Nobs);
	FREE_IVECTOR (IniSp, 1, nparm);
	FREE_IVECTOR (Spec, 1, nparm);
	FREE_IVECTOR (bounded, 1, nparm);
	FREE_VLVECTOR (varsum, 1, 3);
	FREE_ALVECTOR (anasum, 1, 3);
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
** main--main function used to call Nlogist mode fitting program.
*****************************************************************/
int main (int argc, char *argv[])
{

  int iter, i, j, junk;		/*iteration variable */
  int bmdose;			/*flag for computing benchmark dose */
  int Nmiss;			/*number of records with missing values */
  int nparm_known;		/*number of specified parameters */
  int *bounded;                 /*parameters at boundaries*/
  int nvar;                     /*number of unfixed parameters */
  double BMD, lkf, lkr, xlk, W, junk1;	/*log likelihoods */
  double *Parms;		/*parameter array */
  VarList *varsum;		/*info for variables--p. dep.,n. dep., indep. */
  AnaList *anasum;		/*information for ANONA analysis */
  double **vcv;			/*variance and covariance matrix */
  double back, back1;
  double *GYp, *GYn, *GXi;
  char model_name[MNLENGTH], user_note[UNLENGTH], junkname[FLENGTH];
  char dose_name[CNLENGTH], posi_name[CNLENGTH], nega_name[CNLENGTH], junkname1[FLENGTH], junkname2[FLENGTH];
  char long_path_name[FLENGTH];
  char *extpoint;

  time_t ltime;

  time (&ltime);

  /* Min_increment is a practical smallest positive real for purposes */
  /* of doing arithmetic; Max_double is the largest double */
  //Min_increment = DBL_EPSILON;
  //Max_double = DBL_MAX;

  /*open bmdswrk.001 input and Weibull.out output temp files */

  /* allow spaces in the path name */
  if (argc > 2)
    {
      path_name2(argc, argv, long_path_name);
      argv[1] = long_path_name;
    }

  /* **************************************************************** */
  /* Now comes a bunch of code to input run parameters and data and */
  /* initialize various data structures */
  /* **************************************************************** */

  /* Open the input file, or die */
  fp_in = fopen (argv[1], "rt");

  if (fp_in == NULL)
    {
      printf ("Error in opening input  file.\n");
      printf ("...now exiting to system...\n");

      /*     fprintf (fp_out, "Error in opening input file.\n");*/
      /*     fprintf (fp_out, "...Exited to system!\n");*/
      exit (1);
    }

  /* Read the initial run parameters.  The most important here are */
  /* Nobs and ngrp */

  fscanf (fp_in, "%s", model_name);
  fscanf (fp_in, "%[ ^\n]", user_note);
  fscanf (fp_in, "%[^\n]", user_note);
  fscanf (fp_in, "%s", junkname);
  fscanf (fp_in, "%s", junkname);
  fscanf (fp_in, "%d%d", &Nobs, &ngrp);

  /*assign number of parameters */
  nparm = 5 + ngrp;
  /*allocate memory for arrays */
  Parms = DVECTOR (1, nparm);
  IniP = DVECTOR (1, nparm);
  IniSp = IVECTOR (1, nparm);
  Spec = IVECTOR (1, nparm);
  bounded = IVECTOR (1, nparm);

  Xi = DVECTOR (1, Nobs);
  Yp = DVECTOR (1, Nobs);
  Yn = DVECTOR (1, Nobs);
  Ls = DVECTOR (1, Nobs);
  SR = DVECTOR (1, Nobs);
  Xg = IVECTOR (1, Nobs);

  varsum = VLVECTOR (1, 3);
  anasum = ALVECTOR (1, 3);
  Rlevel = DVECTOR (1, 5);
  Bmdl = DVECTOR (1, 5);
  vcv = DMATRIX (1, nparm, 1, nparm);
  GXi = DVECTOR (1, ngrp);
  GYp = DVECTOR (1, ngrp);
  GYn = DVECTOR (1, ngrp);


  junk = 0; /* Used to see if an extension was added to output file name */

  /* Search the output_name array for any extensions, as we don't
     want to have two extensions attached to the file name */

  strcpy (fout, argv[1]);
  extpoint = strrchr(fout,(int) '.');
  if (extpoint != (char *) NULL)
    {

      *extpoint = '\0';
    }
  strcat(fout,".out");

  strcpy (fout2, fout);

  extpoint = strrchr(fout2,(int) '.');
  if (extpoint != (char *) NULL)
    {
      *extpoint = '\0';
    }
  strcat(fout2,".002");

  /*  input more run parameters from the  .(d) file */

  fscanf (fp_in, "%d%lf%lf%d%d%d%d%d%d", &ITMAX, &Rel_Conv, &Parm_Conv, &bmdlCurve,
	  &restrict, &bmdose, &fixedSize, &appendix, &smooth);

  /* sanity check */
  if (bmdose < 0 || bmdose > 1)
    ERRORPRT ("Error in choosing benchmark dose computation.");

  /* Open the two output files, or die */
  if (appendix == Yes)
    fp_out = fopen (fout, "a");
  else
    fp_out = fopen (fout, "w");
#ifndef RBMDS
  fp_out2 = fopen (fout2, "w");
#endif
  if (fp_out == NULL 
#ifndef RBMDS
      || fp_out2 == NULL
#endif
      )
    {
      printf ("Error in opening  output files.\n");
      printf ("...now exiting to system...\n");

      fprintf (fp_out, "Error in opening output files.\n");
      fprintf (fp_out, "...Exited to system!\n");
      exit (1);
    }


  /* Write out the header for the ".out" file */

  fprintf (fp_out, "\n\n ==================================================================== \n");
  fprintf (fp_out, "   \t  %s ", Version_no);
  fprintf (fp_out, " \n  \t  Input Data File: %s  ", argv[1]);
  fprintf (fp_out, "\n \t\t\t\t\t\t\t%s", ctime (&ltime));
  fprintf (fp_out, " ==================================================================== ");
  fprintf (fp_out, "\n\n %s \n", user_note);
  fprintf (fp_out, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n ");

  /* Read BMD and BMDL computation parameters */

  fscanf (fp_in, "%lf%d%lf%d%ld", &bmdparm.effect, &bmdparm.risk, &bmdparm.level, &BSIter, &BSSeed);

  if (BSSeed == 0)              /* Set seed from time clock if default(BSSeed=0) is specified */
    {
    BSSeed = time (NULL);
    }

  /*read values for fixed (specified) parameters */
  READ_PARAMETERS (nparm, Parms);
  /* Have to do the SWAPs because the UI passes parameters in the wrong */
  /* order, and we have hope of changing it someday */
  SWAP (Parms[4], Parms[5], junk1);
  SWAP (Parms[3], Parms[5], junk1);
  SWAP (Parms[2], Parms[5], junk1);

  /* Here is where Spec is filled; Spec tells us which parameters are fixed */
  FILL_SPECVECTOR (nparm, Parms, Spec);
  nparm_known = COUNT_SPECVECTOR (nparm, Spec);

  /* read user-supplied initial parameter values */
  fscanf (fp_in, "%d", &initial);
  READ_PARAMETERS (nparm, IniP);
  SWAP (IniP[4], IniP[5], junk1);
  SWAP (IniP[3], IniP[5], junk1);
  SWAP (IniP[2], IniP[5], junk1);
  FILL_SPECVECTOR (nparm, IniP, IniSp);

  for (i = 1; i <= nparm; i++)
    {
      if (Spec[i] == 1)
	IniP[i] = 1;
    }

  /* read data into Yp, Yn, Xi, Ls, Xg vectors */
  fscanf (fp_in, "%s%s%s%s%s", dose_name, posi_name, nega_name, junkname1, junkname2);
  Nmiss = READ_OBSDATA5V (Nobs, Xi, Yp, Yn, Ls, Xg);

  /* Sort the arrays by dose and move around */
  /* the other arrays appropriately */
  Sort_4_By_Dose (Nobs, Xi, Yn, Yp, Ls);

  junk1 = Xi[1];
  junk = 1;

  for (i = 1; i <= Nobs; i++)	//Hopefully someday, this "group" variable
    {				//will not have to be entered by the user.
      if (Xi[i] == junk1)	//This little loop gets the "group" values
	Xg[i] = junk;		//without explicitly entering them at the
      else			//spreadsheet level.
	{
	  Xg[i] = ++junk;
	  junk1 = Xi[i];
	}
    }

  ngrp = junk;
  /* Adjust the value of Nobs to account for missing values */
  Nobs -= Nmiss;

  //if (Nobs < nparm), commented and changed to the code below, GLN 04/07/08
  if (Nobs < (nparm-nparm_known))
  {
	FreeUp_mem(Parms, varsum, anasum, vcv, GXi, GYp, GYn, bounded);

    ERRORPRT ("Observation # < parameter # for Nlogist model.");
  }

  /* *********** end of input data. */

  /*output title and summary of intput data  *************************** */
      OUTPUT_TEXT ("\n The probability function is: ");
      OUTPUT_TEXT ("\n\n Prob. = alpha + theta1*Rij + [1 - alpha - theta1*Rij]/ ");
      OUTPUT_TEXT ("\n                       [1+exp(-beta-theta2*Rij-rho*log(Dose))],");
      OUTPUT_TEXT ("\n          where Rij is the litter specific covariate.");

      /*indicates restrict power */
      gamm0 = 0;
      if (restrict == Yes)
	{
	  gamm0 = 1.0;
	  OUTPUT_TEXT ("\n Restrict Power rho >= 1. ");
	}

      nparm_known = COUNT_SPECVECTOR (nparm, Spec);
      fprintf (fp_out, "\n\n\n Total number of observations = %d", Nobs + Nmiss);
      fprintf (fp_out, "\n Total number of records with missing values = %d", Nmiss);
      fprintf (fp_out, "\n Total number of parameters in model = %d", nparm);
      fprintf (fp_out, "\n Total number of specified parameters = %d\n\n", nparm_known);

      fprintf (fp_out, "\n Maximum number of iterations = %d\n", ITMAX);
      fprintf (fp_out, " Relative Function Convergence has been set to: %g\n", Rel_Conv);
      fprintf (fp_out, " Parameter Convergence has been set to: %g\n", Parm_Conv);
      fprintf (fp_out, " Number of Bootstrap Iterations per run: %d\n", BSIter);
      fprintf (fp_out, " Bootstrap Seed:  %ld\n\n", BSSeed);


      if (nparm_known > 0)
	{
	  fprintf (fp_out, " User specifies the following parameters:");
	  for (i = 1; i <= nparm; i++)
	    {
	      if (Spec[i] == 1)
		fprintf (fp_out, "\n %15s = %10.5g", Parm_name[i - 1], Parms[i]);
	    }
	  fprintf (fp_out, "\n\n");
	}


      /*compute init_lkf for full model and init_lkr for reduced model */
      lkf = 0.0;
      varsum[1].S = 0;
      varsum[2].S = 0;
      for (i = 1; i <= Nobs; i++)
	{
	  varsum[1].S += Yp[i];
	  varsum[2].S += Yn[i];
	  W = Yp[i] / (Yp[i] + Yn[i]);
	  if (W > 0)
	    lkf += Yp[i] * log (W);
	  if (W < 1)
	    lkf += Yn[i] * log (1 - W);
	}
      W = varsum[1].S / (varsum[1].S + varsum[2].S);
      lkr = varsum[1].S * log (W) + varsum[2].S * log (1 - W);

      /*fitting Nlogist model and output parameter estimators */
      Nlogist_fit (nparm, ngrp, Parms, EPS, &iter, &xlk);
      /* Which parameters are on their bounds? */
      Which_Bounded (Spec, Parms, bounded);

      /* If We didn't get an acceptable fit, don't compute a bmd and */
      /* don't compute standard errors */
      if (ErrorFlag != 0)
	{
	  bmdose = No;
	  for (i = 1; i <= nparm; i++) bounded[i] = 1;
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
	  nvar = Nlogist_vcv (nparm, Spec, Parms, vcv);
	  /* the first nvar x nvar rows and columns of vcv now contain the */
	  /* hessian of the loglikelihood wrt the unfixed parameters */

	  /* Remove the rows and columns of vcv that correspond to bounded */
	  /* parameters */
	  nvar = Take_Out_Bounded_Parms(nvar, bounded, vcv);

	  /* Do the inverse */
	  INVMAT (vcv, nvar);
	}
      /*output Parameter estimates and standard errors, if not all the */
      /*parameters were fixed */
      if (ErrorFlag != -1)
	OUTPUT_DTMS3PARMS (nparm, Spec, bounded, Parms, Parm_name, vcv, 1);

      /* Now compute and output the analysis of deviance */
      DTMS3ANOVA (nparm, Nobs, Spec, lkf, xlk, lkr, anasum, bounded);
      TEMP_ANOVA_OUTPUT (anatxt, anasum);	/*output ANOVA table */
      /* OUTPUT_DTMS3ANOVA (anatxt, anasum); */	/*output ANOVA table */
      fflush(fp_out);

      /*print a goodness of fit table if we fit something*/
      if (ErrorFlag == 0)
	{
	  N_Goodness (ngrp, nparm, Parms, bounded, Nobs, Xi, Yp, Yn,
		      Ls, Xg, SR);
	  fflush(fp_out);
	}

#ifndef RBMDS
      /* output to **.002 *********** */
      fprintf (fp_out2, "\n BMD_flag \t %d \n Nosb \t%d \n nparm \t%d", 
	       bmdose, ngrp, 6);
      fprintf (fp_out2, "\n  Con_lev \t%3.3g ", bmdparm.level);
      fprintf (fp_out2, "\n  RiskType \t%d ", bmdparm.risk);
      fprintf (fp_out2, "\n  Effect \t%3.3g ", bmdparm.effect);
      for (i = 1; i <= 5; i++)
	fprintf (fp_out2, "\n %s \t %5.3g", Parm_name[i - 1], Parms[i]);
      fprintf (fp_out2, "\n fixedSize \t%f ", sijfixed);
      {
	/* Calculation of 95% confidence intervals at each dose level for */
	/* graphical output */
	double *LL, *UL, *phat;

	LL = DVECTOR(1, ngrp);
	UL = DVECTOR(1, ngrp);
	phat = DVECTOR(1, ngrp);

	Nested_CI(ngrp, Nobs, Yp, Yn, Xg, 0.95, LL, phat, UL);

	/* This seems like a lot of work just to get a list of the unique doses */
	for (j = 1; j <= Nobs; j++)
	  {
	    for (i = 1; i <= ngrp; i++)
	      {
		if (Xg[j] == i)
		  {
		    GXi[i] = Xi[j];
		  }
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
      fprintf (fp_out2, "\n Max_Min_dose \n  %f %f \n", xmax, xmin);
      fflush(fp_out2);
#endif

      /*compute benchmark dose */
      if (bmdose == Yes)
	{
	  if (Spec[2] == Yes)
	    {
#ifndef RBMDS
	      fprintf (fp_out2, "\n\n BMDL_comput_ind %d", No);
#endif
	      fprintf (fp_out, 
		       "\n\n  %s parameter is fixed. The likelihood function can not be reparameterized in BMD.", 
		       Parm_name[1]);
	      exit (0);
	    }

	  back = Parms[1] + Parms[3] * sijfixed;
	  back1 = 1 - back;
	  if (bmdparm.risk == 1)
	    back1 = 1;

	  Nlogist_BMD (nparm, Parms, EPS, &junk, xlk, Rlevel, Bmdl, &BMD);
          SRoI(ngrp, Nobs, GXi, Xi, Xg, SR, Ls, sijfixed, BMD);
	}

        N_Bootstrap (ngrp, nparm, Parms, bounded, Nobs, Xi, Yp, Yn,
		      Ls, Xg, SR, BSIter, BSSeed);

        /*** output BMD & BMDL */
	if (bmdose == Yes)
	{
	  if (BMD >= 0.0)
	    {
	      if (fixedSize == 1)
		fprintf (fp_out,
			 "\n\n\n\nTo calculate the BMD and BMDL, the litter specific covariate is fixed\n at the mean \
litter specific covariate of control group: %f", sijfixed);
	      else
		fprintf (fp_out,
			 "\n\n\n\nTo calculate the BMD and BMDL, the litter specific covariate is fixed\n at the mean \
litter specific covariate of all the data: %f", sijfixed);
	      OUTPUT_BENCHMD (1, BMD);
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
#ifndef RBMDS
	      fprintf (fp_out2, "\n  RSL \t%f", bmdparm.effect * back1 + back);
	      fprintf (fp_out2, "\n  BMD \t%f", BMD);
	      fprintf (fp_out2, "\n  BMDL \t%f", Bmdl[1]);

	      fprintf (fp_out2, "\n\n BMD_line");
	      fprintf (fp_out2, "\n %f %f", (xmin - xmax / 100), 
		       bmdparm.effect * back1 + back);
	      fprintf (fp_out2, "\n %f %f", BMD, bmdparm.effect * back1 + back);
	      fprintf (fp_out2, "\n %f %f", BMD, -0.1);

	      fprintf (fp_out2, "\n\n BMDL_line");
	      fprintf (fp_out2, "\n %f %f", Bmdl[1], -0.1);
	      fprintf (fp_out2, "\n %f %f", Bmdl[1], bmdparm.effect * back1 + back);

	      fprintf (fp_out2, "\n\n BMDL_Curve_flag \t %d  \n smooth_opt  %d", bmdlCurve, smooth);
	      fprintf (fp_out2, "\n\n BMDL_curve");
	      fprintf (fp_out2, "\n 0.00000 %f", back);
	      for (i = 1; i <= 5; i++)
		{
		  fprintf (fp_out2, "\n %f %f", Bmdl[i], Rlevel[i] * back1 + back);
		}
#endif
	    }

	}


      /* indicate all required data have been output. */
#ifndef RBMDS
      fprintf (fp_out2, "\n\n Check_result %d", Yes);
#endif
	FreeUp_mem(Parms, varsum, anasum, vcv, GXi, GYp, GYn, bounded);
      //FREE_DVECTOR (Parms, 1, nparm);
      //FREE_DVECTOR (IniP, 1, nparm);
      //FREE_DVECTOR (Xi, 1, Nobs);
      //FREE_DVECTOR (Yp, 1, Nobs);
      //FREE_DVECTOR (Yn, 1, Nobs);
      //FREE_IVECTOR (Ls, 1, Nobs);
      //FREE_IVECTOR (Xg, 1, Nobs);
      //FREE_IVECTOR (IniSp, 1, nparm);
      //FREE_IVECTOR (Spec, 1, nparm);
      //FREE_IVECTOR (bounded, 1, nparm);
      //FREE_VLVECTOR (varsum, 1, 3);
      //FREE_ALVECTOR (anasum, 1, 3);
      //FREE_DVECTOR (Rlevel, 1, 5);
      //FREE_DVECTOR (Bmdl, 1, 5);
      //FREE_DMATRIX (vcv, 1, nparm, 1, nparm);
      //FREE_DVECTOR (GXi, 1, ngrp);
      //FREE_DVECTOR (GYp, 1, ngrp);
      //FREE_DVECTOR (GYn, 1, ngrp);

      /*close opened temp files */
      if (fclose (fp_in) != 0 || fclose (fp_out) != 0 
#ifndef RBMDS
	  || fclose (fp_out2) != 0
#endif
	  )
	{
	  ERRORPRT ("Error in closing opened files.");
	}

      return 0;

}				/*end of main */

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

  /* There are many parameter restrictions in the Nlogist model.*/
  /* The if statements */
  /* below address the parameter space restrictions, and adjust the */
  /* degrees of freedom appropriately */

  //               _
  //alpha + Theta1*Sij >= 0:

  if ((Spec[1] != 1 || Spec[3] != 1) && VERYCLOSE(Parms[1], -Parms[3] * maxsize))
    {
      if (Spec[1] == 1) bounded[3] = 1;
      if (Spec[3] == 1) bounded[1] = 1;
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
	  bkg = Parms[(int) alpha] + Parms[(int) theta1]*Lsc[i];
	  if (doses[i] <= 0.0)
	    P[i] = bkg;
	  else
	    P[i] = bkg +
	      (1.0 - bkg)/
	      (1.0 + exp(-(Parms[(int) beta] + Parms[(int) theta2] * Lsc[i] +
			   Parms[(int) rho] * log(doses[i])
			   )
			 )
	       );
	}
}
/*******************************************************************
**Nlogist_probs -- computes litter-specific probabilities of an
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
void
Nlogist_probs (double probs[], double p[], int compgrad, double **gradij)
{
  int i;
  double *pint, ex, ex1, ex2, ex3, ex4, dd2, spij, smij;

  pint = DVECTOR(1, nparm);
  for (i = 1; i <= nparm; i++)
    {
      pint[i] = p[i];
    }
  if(replace == Yes) /* We're computing a BMDL */
    {
      if (Spec[1] == Spec[3]) /* alpha and theta1 are transformed */
	{
	  if (bmdparm.risk == 1)
	    {
	      pint[2] = -log ((1 - pint[1] * spfixed - pint[3] * snfixed) /
			     BMR - 1) -
		pint[4] * sijfixed - pint[5] * log (tD);
	    }
	  else
	    {
	      pint[2] = log (BMR/(1.0 - BMR)) - pint[4] * sijfixed -
		pint[5] * log (tD);
	    }
	}
      else /* alpha and theta1 are not transformed */
	{
	  if (bmdparm.risk == 1)
	    {
	      pint[2] = -log ((1 - pint[1] - pint[3] * sijfixed) / BMR - 1) -
		pint[4] * sijfixed - pint[5] * log (tD);
	    }
	  else
	    {
	      pint[2] = log (BMR/(1.0 - BMR)) - pint[4] * sijfixed -
		pint[5] * log (tD);
	    }
	}
    }
  for (i = 1; i <= Nobs; i++)
    {
      spij = (smax - Ls[i])/sdif;
      smij = (Ls[i] - smin)/sdif;
      if (Spec[1] == Spec[3])
	{
	  ex = spij * pint[1] + smij * pint[3];
	}
      else
	{
	  ex = pint[1] +  Ls[i] * pint[3];
	}
      ex2 = 1.0 - ex;
      ex1 = 0.0;
      if (Xi[i] > 0)
	{
	  ex1 = exp ( -(pint[2] + pint[4] * Ls[i] +
			pint[5] * log (Xi[i])));
	  ex = ex + ex2 / (1 + ex1);
	}
      PROBABILITY_INRANGE (&ex);
      probs[i] = ex;
      if (compgrad == 1)
	{
	  ex3 = ex1 / ((1.0 + ex1) * (1.0 + ex1));
	  dd2 = ex2 * ex3;
	  if (replace == Yes)
	    {
	      if (bmdparm.risk == 1)
		{
		  ex4 = 1.0 / (ex2 - BMR);
		}
	      else
		{
		  ex4 = 0.0;
		}
	      if (Spec[1] == Spec[3])
		{
		  if (Xi[i] > 0)
		    {
		      gradij[i][1] = (spij * (1.0 - 1.0/(1.0 + ex1)) +
				      dd2 * ex4 * spfixed);
		      gradij[i][3] = (smij * (1.0 - 1.0 / (1.0 + ex1)) +
				      dd2 * ex4 * snfixed);
		    }
		  else /* Xi[i] == 0 */
		    {
		      gradij[i][1] = spij;
		      gradij[i][3] = smij;
		    }
		}
	      else /* Spec[1] != Spec[3] */
		{
		  if (Xi[i] > 0)
		    {
		      gradij[i][1] = 1.0 - 1.0 / (1.0 + ex1) + dd2 * ex4;
		      gradij[i][3] = Ls[i] * gradij[i][1];
		    }
		  else /* Xi[i] == 0 */
		    {
		      gradij[i][1] = 1.0;
		      gradij[i][3] = Ls[i];
		    }
		} /* end of derivatives that depend on the transformation */
		  /* of alpha and theta1*/
	      gradij[i][2] = 0.0;  /* this is the replaced parameter */
	      if (Xi[i] > 0)
		{
		  gradij[i][4] = dd2 * (Ls[i] - sijfixed);
		  gradij[i][5] = dd2 * (log(Xi[i]) - log(tD));
		}
	      else /* Xi[i] == 0 */
		{
		  gradij[i][4] = 0.0;
		  gradij[i][5] = 0.0;
		}
	    } /* end of replace == Yes */
	  else /* replace == No */
	    {
	      if (Spec[1] == Spec[3])
		{
		  if (Xi[i] > 0)
		    {
		      gradij[i][1] = spij * (1.0 - 1.0 / (1.0 + ex1));
		      gradij[i][3] = smij * (1.0 - 1.0 / (1.0 + ex1));
		    }
		  else /* Xi[i] == 0 */
		    {
		      gradij[i][1] = spij;
		      gradij[i][3] = smij;
		    }
		}
	      else /* Spec[1] != Spec[3] */
		{
		  if (Xi[i] > 0)
		    {
		      gradij[i][1] = 1.0 - 1.0 / (1.0 + ex1);
		      gradij[i][3] = Ls[i] * gradij[i][1];
		    }
		  else /* Xi[i] == 0 */
		    {
		      gradij[i][1] = 1.0;
		      gradij[i][3] = Ls[i];
		    }
		} /* end of derivatives that depend on the transformation */
		  /* of alpha and theta1*/
	      if (Xi[i] > 0)
		{
		  gradij[i][2] = dd2;
		  gradij[i][4] = dd2 * Ls[i];
		  gradij[i][5] = dd2 * log(Xi[i]);
		}
	      else /* Xi[i] == 0 */
		{
		  gradij[i][2] = 0.0;
		  gradij[i][4] = 0.0;
		  gradij[i][5] = 0.0;
		}
	    } /* end of replace == No */
	}
    }

  FREE_DVECTOR(pint, 1, nparm);
} /* end of Nlogist_probs */

/*******************************************************************
**Nlogist_lk -- used to compute the log likelihood for Nlogist model.
*               parameters are in "internal" transformed form.
*
* 		Extern var.: smean, smax, Nobs, Xi, Yp, Yn, Ls, Xg.
*
*********************************************************************/
void
Nlogist_lk (long int *nvar, double *x, long int *nf, double *f,
	    long int *uiparm, double *urparm, void (*ufparm) ())
{
  double xlk;			// log likelihood.
  int    i, j, k, plus5, jfixed, jvar, compgrad;
  double tm1, tm2, tm3, tm;	//temp var.
  double *p;			//for "untransform" parms.
  double *probs;                /* litter-specific probabilities */
  double **gradij;

  p = DVECTOR (1, nparm);
  probs = DVECTOR(1, Nobs);
  gradij = DMATRIX(1, Nobs, 1, 5);

  jfixed = jvar = 0;
  for (j = 1; j <= nparm; j++)	/*reconstruct the parameter vector. */
    {
      if (Spec[j] == Yes)
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
  Nlogist_probs(probs, p, compgrad, gradij);

  xlk = 0.0;
  for (i = 1; i <= Nobs; i++)
    {
      tm1 = 0.0;
      tm2 = 0.0;
      tm3 = 0.0;
      plus5 = 5 + Xg[i];
      j = (int) Yp[i];
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
  FREE_DVECTOR (p, 1, nparm);
  FREE_DVECTOR (probs, 1, Nobs);
  FREE_DMATRIX(gradij, 1, Nobs, 1, 5);
  *f = -xlk;
}


/*******************************************************************
**Nlogist_g -- used to compute the gradients for Nlogist model.
*              Parameters are in "internal" transformed form.
*		Extern var.: smean, smax, Nobs, Xi, Yp, Yn, Ls, Xg.
*
**********************************************************************/
void
Nlogist_g (long int *nvar, double *x, long int *nf, double *g,
	   long int *uiparm, double *urparm, void (*ufparm) ())
{
  double ex;		//temp var.
  double tm1, tm2, tm3, tm1a, tm2a, tm3a, tm, tm12;	//temp var.
  double *dd, *p, *probs, **gradij, *tmp_g;
  int i, j, k, plus5, jfixed, jvar, compgrad;

  jfixed = jvar = 0;

  p = DVECTOR (1, nparm);
  probs = DVECTOR (1, Nobs);
  gradij = DMATRIX(1,Nobs,1,5);
  dd = DVECTOR (1, nparm);
  tmp_g = DVECTOR(1, nparm);

  for (j = 1; j <= nparm; j++)	/* reconstruct the parameter vector. */

    {
      if (Spec[j] == Yes)
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
  Nlogist_probs (probs, p, compgrad, gradij);


  /** initial tmp_g[j]'s            **************/
  for (j = 1; j <= nparm; j++)
    tmp_g[j] = 0.0;
  for (i = 1; i <= Nobs; i++)
    {
      ex = probs[i];

      /*Compute first partial derivatives */
      tm1 = 0.0;
      tm2 = 0.0;
      tm3 = 0.0;
      tm1a = 0.0;
      tm2a = 0.0;
      tm3a = 0.0;
      for (j = 6; j <= nparm; j++)
	dd[j] = 0.0;
      plus5 = 5 + Xg[i];
      j = (int) Yp[i];
      if (ex > 0.0)
	{
	  for (k = 1; k <= j; k++)
	    {
	      tm = ex + (k - 1) * p[plus5];
	      tm1 += 1.0 / tm;
	      tm1a += (1.0 / tm) * (k - 1);
	    }
	}
      j = (int) Yn[i];
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
      for (j = 1; j <= nparm; j++)
	{
	  tmp_g[j] -= dd[j];
	}
    }
  /* end of first partial deri. */
  jvar = 0;
  for (j = 1; j <= nparm; j++)	//reconstruct the parameter vector.

    if (Spec[j] == No)
      {
	g[jvar] = tmp_g[j];
	jvar++;
      }

  FREE_DVECTOR (p, 1, nparm);
  FREE_DVECTOR (probs, 1, Nobs);
  FREE_DMATRIX (gradij, 1, Nobs, 1, 5);
  FREE_DVECTOR (dd, 1, nparm);
  FREE_DVECTOR (tmp_g, 1, nparm);
  return;
}
/****************************************************************************
 **
 * NLogist_grad -- Computes the gradient of the nested logistic likelihood
 * function with respect to the user form of the parameters.  This is to
 * be used in NLogist_vcv, to compute a finite difference approximation to
 * the hessian of the likelihood function
 * Input: nparm -- the number of parameters, the dimension of all the following
                   arrays.
          Spec[] -- if the ith parameter is fixed by the user, Spec[i] == 1,
	            otherwise Spec[i] == 0.
	  ptf[] -- vector of parameters, in user form (external form),
	           based at 1.
 * Output: grad[] -- the gradient of the loglikelihood function
                     (N.B. NLogist_g returns the gradient of -loglikelihood)
		     Based at 1.

 *****************************************************************************/
void NLogist_grad(int nparm, int Spec[], double ptf[],
		  double grad[])
{
  long int *uiparm, nvar, i, j, jfixed, jvar, nf;
  double *urparm, *start, *outgrad, *p, junk1, junk3;
  void (*ufparm) ();

  /* set up initial parameter array, start.  All parameters go either into */
  /* start (if they are changing to improve the fit) or urparm (if they are */
  /* fixed).*/

  p = DVECTOR(1, nparm);
  for (i = 1; i <= nparm; i++)
    {
      p[i] = ptf[i];
    }
  /* ------- Transform the parameters to the "internal" form -------    */
  for (j = 6; j <= nparm; j++)
    p[j] = p[j] / (1 - p[j]);	// Phi --> Psi.

//  if (Spec[1] == Spec[3])
//    {
//      junk1 = p[1];
//      junk3 = p[3];
//      p[1] = junk1 + smin * junk3;
//      p[3] = junk1 + smax * junk3;
//    }

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
  Nlogist_g (&nvar, start, &nf, outgrad, uiparm, urparm, ufparm);
  /* Nlogist_g returns the gradient of -loglikelihood */
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
*	NLogist_vcv -- used to compute the hessian matrix for Nlogist of the
*                      log-likelihood with respect to the parameters.  Uses
*                      a finite-difference (central) derivative of the
*                      gradient.
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

int Nlogist_vcv (int nparm, int Spec[], double ptf[], double **vcv)
{
  double *ptemp, *saveparms, *h, *gradp, *gradm, hrat, tmp, junk1, junk3;
  int i, j, ivar, jvar, nvar;

  /* initialize memory for all the arrays */
  ptemp = DVECTOR(1, nparm);
  saveparms = DVECTOR(1, nparm);
  h = DVECTOR(1, nparm);
  gradp = DVECTOR(1,nparm);
  gradm = DVECTOR(1, nparm);


  ptemp = DVECTOR(1, nparm);
  for (i = 1; i <= nparm; i++)
    {
      ptemp[i] = ptf[i];
    }

  if (Spec[1] == Spec[3])
    {
      junk1 = ptemp[1];
      junk3 = ptemp[3];
      ptemp[1] = junk1 + smin * junk3;
      ptemp[3] = junk1 + smax * junk3;
    }


  /* Get a value of h for each parameter */
  hrat = pow(1.0e-16, 0.333333);

  for (i = 1; i <= nparm; i++)
  {
    if (fabs(ptemp[i]) > 1.0e-7)
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
  /* for each unfixed parameter, compute the second derivative wrt each */
  /* of the others */
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
	NLogist_grad(nparm, Spec, saveparms, gradp);
	saveparms[i] = ptemp[i] - h[i];
	NLogist_grad(nparm, Spec, saveparms, gradm);
	/* Now compute the second derivative */
	jvar = 0;
	for (j = 1; j <= nparm; j++)
	  {
	    if (Spec[j] != Yes)
	      {
		jvar++;
		vcv[ivar][jvar] = -(gradp[jvar] - gradm[jvar])/(2.0 * h[i]);
	      }
	  }
      }
  }

  FREE_DVECTOR(ptemp, 1, nparm);
  FREE_DVECTOR(saveparms,1,nparm);
  FREE_DVECTOR(h,1,nparm);
  FREE_DVECTOR(gradp,1,nparm);
  FREE_DVECTOR(gradm,1,nparm);
  return nvar;
}

/**************************************************************
*MAX_lk -- used to obtain the Maximum log-likelihood as well as
*           the estimates of parameters, given initial p[1..n],
*           object func. , and gradient func. G_func.
*
**************************************************************/
void
MAX_lk (int nparm, double p[], double gtol, int *iter, double *fret)
{
  int i, jfixed, jvar;
  long int nvar;
  long int *uiparm;

  double *start;
  double *urparm, *lower, *upper;
  long int dummy;
  void (*ufparm) ();

  /* Set up initial parameter array, start.  All parameters go either
     into start (if they are changing to improve the fit) or urparm
     (if they are fixed). */

  nvar = nparm;
  for (i = 1; i <= nparm; i++)
    nvar = nvar - Spec[i];	/* Count the varying parameters */


  if (nvar > 0) /* we will do some optimizing */
    {
      urparm = (double *) malloc ((size_t) (nparm - nvar) * sizeof (double));
      start = (double *) malloc ((size_t) nvar * sizeof (double));
      lower = (double *) malloc ((size_t) nvar * sizeof (double));
      upper = (double *) malloc ((size_t) nvar * sizeof (double));
    }
  else /* all parameters are fixed.  Return the likelihood, set */
       /* ErrorFlag to -1, and return */
    {
      urparm = (double *) malloc ((size_t) (nparm - nvar) * sizeof (double));
      for (i = 1; i <= nparm; i++)
	urparm[i - 1] = p[i];

      /* in what follows, Spec was originally placed as the 5th arg. */
      /* Spec is the wrong type, and is never used in Nlogist_lk, anyway, */
      /* so I replaced it with uiparm */
      Nlogist_lk (&nvar, start, &dummy, fret,
		  uiparm, urparm, ufparm);
      *fret = -(*fret);
      ErrorFlag = -1;
      return;
    }

  jfixed = 0;
  jvar = 0;
  for (i = 1; i <= nparm; i++)	//seperate the fixed and variable parameters

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
  if (Spec[1] != 1)
    {
      lower[jvar] = 0.0;
      upper[jvar] = 1;
      jvar++;
    }
  if (Spec[2] != 1)
    {
      lower[jvar] = -Max_double;
      upper[jvar] = Max_double;
      jvar++;
    }
  if (Spec[3] != 1)
    {
      lower[jvar] = 0.0;
      upper[jvar] = 1.0;
      jvar++;
    }

  if (Spec[4] != 1)
    {
      lower[jvar] = -Max_double;
      upper[jvar] = Max_double;
      jvar++;
    }
  if (Spec[5] != 1)
    {
      if (restrict == 1)
	lower[jvar] = 1.0;
      else
	lower[jvar] = 0.0;
      upper[jvar] = SlopeUpperBound;
      jvar++;
    }
  for (i = 6; i <= nparm; i++)
    {
      if (Spec[i] != Yes)
	{
	  lower[jvar] = 0.0;
	  upper[jvar] = Max_double;
	  jvar++;
	}
    }
  Maxloglik = 0.0; /* we don't have a good value here */

  /* maximise the log-likelihood.  ErrorFlag gives the return code */
  ErrorFlag = run_dmngb((int) nvar, start, lower, upper, Maxloglik,
			Rel_Conv, Parm_Conv, ITMAX, 10,
			Nlogist_lk,Nlogist_g,uiparm,urparm,ufparm,
			DeBuG,fret);
  /* We actually minimized -log-likelihood.  Want to return log-likelihood */

  *fret = -*fret;

  /* Put the parameter values back in p */
  jfixed = jvar = 0;
  for (i = 1; i <= nparm; i++)
    {
      if (Spec[i] == Yes)
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

  if (nvar > 0)
    {
      free (start);
      free (urparm);
      free (lower);
      free (upper);
    }
}

/**************************************************************
*Nlogist_fit -- Used to "prepare" the data for further computation,
*            i.e. compute the extern variables, give the initial
*            parameters, etc. THEN fit the Nlogist model.
*            (In fact, these jobs could be done in main().)
*
***************************************************************/
void
Nlogist_fit (int n, int ngrp, double p[], double gtol,
	     int *iter, double *fret)
{
  int *SpBak, *Spec2, i, j, junk, count;
  double sum1, nij, ymin, W, xlk, x, junk1, junk3;
  double *tmYp, *tmYn, *tmXi, *pBak, *tmy, **tmvcv;

  tmvcv = DMATRIX (1, 2, 1, 2);
  tmy = DVECTOR (1, 2);
  Spec2 = IVECTOR (1, 2);
  tmYn = DVECTOR (1, ngrp);
  tmYp = DVECTOR (1, ngrp);
  tmXi = DVECTOR (1, ngrp);
  pBak = DVECTOR (1, n);
  SpBak = IVECTOR (1, n);

  replace = No;

  /* -----------------------------------------------------------------------
  **   Compute statistics for litter-specific covariate:
  **   smean1, smean, smin, smax,   sdif
  ** ----------------------------------------------------------------------*/

  sum1 = 0.0;
  nij = 0.0;
  i = 1;
  while (Xg[i] == 1)
    {
      sum1 += Ls[i];
      nij += 1.0;
      i++;
    }
  smean1 = sum1 / nij;		/*the average lsc in group 1.*/

  sum1 = 0.0;
  nij = 0.0;
  smax = Ls[1];
  smin = Ls[1];
  xmax = 0.0;
  for (i = 1; i <= Nobs; i++)
    {
      x = Xi[i];
      sum1 += Ls[i];
      nij += 1.0;
      if (Ls[i] > smax)
	smax = Ls[i];
      if (Ls[i] < smin)
	smin = Ls[i];
      if (x > xmax)
	xmax = x;
    }
  smean = sum1 / nij;		//overall average lsc.

  sdif = smax - smin;

  if (fixedSize == 1)
    sijfixed = smean1;
  else
    sijfixed = smean;

  /* -----------------------------------------------------------------------*/
  /*                            Set up parameters                           */

  if (initial == Yes) /* Inserting user-specified initial values */
    {
      for (i = 1; i <= n; i++)
	SpBak[i] = 1 - IniSp[i];
      //Have to do this because the fun OUTPUT_Init.
      OUTPUT_TEXT ("\n\n                 User Inputs Initial Parameter Values  ");
      OUTPUT_Init (n, SpBak, IniP, Parm_name);
      //obtain user input initial values for unspecified parms.
      for (i = 1; i <= n; i++)
	{
	  if (IniSp[i] == 1)	// have been initialized.

	    {
	      if (Spec[i] == 1)	// check if it is for fixed parm.

		Warning ("The initial value for the fixed parameter is ignored.");
	      //p[i]=p[i]. no change.
	      else
		p[i] = IniP[i];
	    }
	  else
	    {
	      //check if all the unspecifird parms were initialized.
	      if (Spec[i] == 0)
		ERRORPRT ("When the initial option is chosen, one has to initial ALL unspecified parameters.");
	    }
	}

      if (p[5] < restrict || p[1] + smax * p[3] < 0 || p[1] + smax * p[3] > 1 || p[1] + smin * p[3] < 0 || p[1] + smin * p[3] > 1)
	ERRORPRT ("The initial values have to be:  rho >= 0 (or 1 when there is restriction on rho) and 0 < alpha+theta1*Sij < 1 for ALL Sij. ");

      for (j = 6; j <= n; j++)
	if (p[j] < 0 || p[j] >= 1)
	  ERRORPRT ("The initial values have to be: 0 <= Phi[j] < 1, for the correlation parameters.");
    }
  else /* Here, we compute default starting values for all unfixed parameters */
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
      for (j = 1; j <= ngrp; j++)	/* convert nest_data to dicho_data. */

	{
	  tmYn[j] = 0;
	  tmYp[j] = 0;
	  while (i <= Nobs && Xg[i] == j)
	    {
	      tmYn[j] += Yn[i];
	      tmYp[j] += Yp[i];
	      tmXi[j] = Xi[i];
	      i++;
	    }
	}

      /*compute initial estimates */
      INITIALIZE_DMATRIX (tmvcv, 2, 2);
      INITIALIZE_DVECTOR (tmy, 2, 1);
      ymin = 1.0;
      for (i = 1; i <= ngrp; i++)
	{
	  W = tmYp[i] / (tmYp[i] + tmYn[i]);
	  PROBABILITY_INRANGE (&W);
	  if (W < ymin)
	    ymin = W;
	}

      for (i = 1; i <= ngrp; i++)
	{
	  x = tmXi[i];
	  /*	  if (x > 0) */
	    {
	      W = log ((tmYp[i] - ymin*(tmYp[i] + tmYn[i]) + 0.5)/(tmYn[i] + 0.5));
	      if (x <= CloseToZero)
		x = Log_zero;
	      else
		x = log (x);
	      tmvcv[1][1] += 1.0;
	      tmvcv[1][2] += x;
	      tmvcv[2][2] += x * x;
	      tmy[1] += W;
	      tmy[2] += W * x;
	    }
	}
      tmvcv[2][1] = tmvcv[1][2];
      Spec2[1] = SpBak[2];
      Spec2[2] = SpBak[5];
      p[1] = ymin + 0.001;	/* in case ymin=0 */
      if (Spec2[1] + Spec2[2] == 0) {
	if (0 == INVMAT (tmvcv, 2))
	{
	  p[2] = tmvcv[1][1] * (tmy[1]) + tmvcv[1][2] * (tmy[2]);
	  p[5] = tmvcv[2][1] * (tmy[1]) + tmvcv[2][2] * (tmy[2]);
	}
      else
	{
	  p[5] = 1.001;
	  p[2] = -1.0 ;
	}
      }
      if (restrict == Yes && p[5] < 1.000)
	p[5] = 1.0001;

      for (i = 1; i <= 2; i++)
	if (SpBak[i] == 1)
	  p[i] = pBak[i];
      if (SpBak[5] == 1)
	p[5] = pBak[5];

      for (i = 6; i <= n; i++)	/* setup Logist model */

	{
	  Spec[i] = 1;
	  p[i] = 0.0;
	}
      Spec[3] = Spec[4] = 1;
      p[3] = p[4] = 0.0;

      /* This is a simple, non-nested log-logistic model */

      MAX_lk (n, p, gtol, &junk, &xlk);

      /* Now p[1], p[2], and p[5] contain starting estimates for */
      /* alpha, beta, and rho  */


      /*  Second, get initial values for Phi's. */
      count = 0;
      for (i = 6; i <= n; i++)
	count += SpBak[i];
      if (count < ngrp)
	{
	  for (i = 6; i <= n; i++)	/* Spec[3], [4] remain as 1, so theta1 and theta2 stay 0.0 */

	    {
	      Spec[i] = SpBak[i];	/* set back to original ones. */

	      p[i] = pBak[i];
	      if (SpBak[i] == 0)
		p[i] = 0.01;	/* set init. if not fixed. 0.45 seems large to me. */

	    }
	  for (j = 6; j <= n; j++)
	    p[j] = p[j] / (1 - p[j]);	// Phi --> Psi.



	  MAX_lk (n, p, gtol, &junk, &xlk);

	  /* Transform parameters to "external" form */

	  for (j = 6; j <= n; j++)
	    p[j] = p[j] / (1 + p[j]);	// Psi --> Phi.

	  /* When theta1 == 0, internal and external forms for alpha are the same */
	  /* so we do not need to back transform p[1] */

	}

      /* Finally, get starting values for Theta's. */
      for (i = 3; i <= 4; i++)
	{
	  Spec[i] = SpBak[i];
	  p[i] = pBak[i];
	  if (SpBak[i] == 0)
	    p[i] = 0.0;
	}

      OUTPUT_TEXT ("\n\n                  Default Initial Parameter Values  ");
      OUTPUT_Init (n, Spec, p, Parm_name);
    }
  /* ------------------------------------------------------------------ */
  /*                                  Fit the model                     */

  /* ------- Transform the parameters to the "internal" form -------    */
  for (j = 6; j <= n; j++)
    p[j] = p[j] / (1 - p[j]);	// Phi --> Psi.

  if (Spec[1] == Spec[3])
    {
      junk1 = p[1];
      junk3 = p[3];
      p[1] = junk1 + smin * junk3;
      p[3] = junk1 + smax * junk3;
    }

  /* --------- ML fit and return log-likelihood value --------   */

  MAX_lk (n, p, gtol, &junk, &xlk);

  /* Print out a warning if ErrorFlag is non-zero */
  do_dmngb_warning(&ErrorFlag);

  *fret = xlk;

  /* ---------- Transform the parameters to the "external" form ------- */

  for (j = 6; j <= n; j++)
    p[j] = p[j] / (1 + p[j]);	// Psi --> Phi.

  if (Spec[1] == Spec[3])
    {
      junk1 = p[1];
      junk3 = p[3];
      p[1] = (smax * junk1 - smin * junk3) / (sdif);
      p[3] = (-junk1 + junk3) / (sdif);
    }

  /*             ----------- Free all the memory used --------                                */

  FREE_DVECTOR (tmYn, 1, ngrp);
  FREE_DVECTOR (tmYp, 1, ngrp);
  FREE_DVECTOR (tmXi, 1, ngrp);
  FREE_DVECTOR (pBak, 1, n);
  FREE_IVECTOR (SpBak, 1, n);
  FREE_DMATRIX (tmvcv, 1, 2, 1, 2);
  FREE_DVECTOR (tmy, 1, 2);
  FREE_IVECTOR (Spec2, 1, 2);
}


/************************************************************
* Nlogist_BMD -- Used to calculate the BMD and BMDL for Nlogist model.
*                Input parameters are in "external" untransformed form.
*
*************************************************************/
void
Nlogist_BMD (int nparm, double p[], double gtol, int *iter, double xlk,
	     double Rlevel[], double Bmdl[], double *BMD)
{
  double tol;
  double xa, xb, fa, fb;
  double stepsize;
  double D, Bmdjunk, junk1, junk3, effect;
  double *pBak, *pa, *pb, *pred, *porig, *doses, *lsc;
  int i, j, k, bogusBMD;

  pBak = DVECTOR (1, nparm);
  pa = DVECTOR(1, nparm);
  pb = DVECTOR(1, nparm);
  porig = DVECTOR(1, nparm);

/**** compute X^2 value  ************************/
  /* IF ML is the value of the maximized log-likelihood, then ML - LR is the value
     log-likelihood at the BMDL or BMDU */
  if (bmdparm.level < 0.5)
    LR = QCHISQ(1.0 - 2*bmdparm.level,1)/2.0;
  else
    LR = QCHISQ(2*bmdparm.level - 1.0,1)/2.0;

  for (j = 1; j <= nparm; j++)
    {
    pa[j] = p[j];
    porig[j] = p[j];
    }

  /* Transform parameters into "internal" form */
  for (j = 6; j <= nparm; j++)
    p[j] = p[j] / (1 - p[j]);	// Phi --> Psi.

  if (Spec[1] == Spec[3])
    {
      junk1 = p[1];
      junk3 = p[3];
      p[1] = junk1 + smin * junk3;
      p[3] = junk1 + smax * junk3;
    }
  /* Finished transforming the parameters */

  for (j = 1; j <= nparm; j++)
    pBak[j] = p[j];		//save the p[].

  Rlevel[1] = BMR = bmdparm.effect;

  spfixed = (smax - sijfixed) / (sdif);
  snfixed = (sijfixed - smin) / (sdif);

  if (bmdparm.risk == 1)
    {
      if (Spec[1] == Spec[3])
	{
	  ck = -1.0 * log ((1 - pBak[1] * spfixed - pBak[3] * snfixed) / BMR - 1);
	}
      else
	{
	  ck = -1.0 * log ((1 - pBak[1] - pBak[3] * sijfixed) / BMR - 1);
	}
    }
  else
    ck = -log ((1 - BMR) / BMR);

/**** solve the BMD ********************************/
  bogusBMD = 0;
  if (p[5] <= (ck - pBak[2] - pBak[4] * sijfixed) / 250)
    {
      *BMD = 100 * xmax;
      bogusBMD = 1;
      Warning ("The power parameter is essentially zero. BMD set to 100 * max(Dose).\n");
    }
  else
    xb = exp ((ck - pBak[2] - pBak[4] * sijfixed) / pBak[5]);
  *BMD = xb;

  pred = DVECTOR(1,4);
  doses = DVECTOR(1,4);
  lsc = DVECTOR(1,4);
  doses[1] = 0.0;
  doses[2] = *BMD;
  lsc[1] = lsc[2] = sijfixed;
  Predict(doses, lsc, 2, pa, pred);   /* call this to make Nlogist_probs cleaner */
/*    Nlogist_probs(2, doses, lsc, pred, nparm, p, 0, dummy); */

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
      fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  No);
      fprintf(fp_out,
	      "\nComputed BMD is %g; response at the computed BMD is %g\n",
	      *BMD,pred[2]);
      fprintf (fp_out,
	       "Computed effect at this estimate is: %g, requested \
effect is %g\n",effect, BMR);
      *BMD = -1.0;
      FREE_DVECTOR(pred,1,4);
      FREE_DVECTOR(doses,1,4);
      FREE_DVECTOR(lsc,1,4);
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
	  FREE_DVECTOR(pred,1,4);
	  FREE_DVECTOR(doses,1,4);
	  FREE_DVECTOR(lsc,1,4);
	  return;
	}
    }
  FREE_DVECTOR(pred,1,4);
  FREE_DVECTOR(doses,1,4);
  FREE_DVECTOR(lsc,1,4);

  /* OK -- We have a legitimate BMD */

/********* search for BMDL **************************/
  Spec[2] = 1;
  replace = Yes;		/* replace is extern var., has to be changed now. */
  stepsize = 0.5;  //Start close to the BMD and work downwards.
  xa = xb * stepsize;
  tol = FMAX ((*BMD) * 0.001, 0.0000001);
  // tol = sqrt(DBL_EPSILON);
  BMD_lk = xlk;			//get the lk at BMD.

  fb = DBL_MAX;
  for (i = 1; i <= nparm; i++) pa[i] = pb[i] = p[i];
  fa = BMDL_func (nparm, pa, xa, tol);

  /* Look for a value of xa on the other side of the BMDL.  We know we're there when fa > 0.
   * Stop if xa gets too small, or the profile likelihood gets flat (fabs(fa - fb) too small).
   */ 
  while (fa<0.0 && xa > DBL_MIN && fabs(fa - fb) > DBL_EPSILON)
    {
      xb = xa;
      fb = fa;
      for (i = 1; i <= nparm; i++) pb[i] = pa[i];
      xa *= stepsize;

      fa = BMDL_func (nparm, pa, xa, tol);
    }
  // if we had to stop, and fa is still negative, somethings gone wrong.  Bail out.
  if (fa < 0.0)
    {
      Bmdl[1] = -1.0;
      return;
    }
  else // BMDL is between xa and xb.
    {
      Bmdl[1] = zeroin(xa, xb, 1.0e-10, BMDL_func, nparm, pb, 1.0e-14);
      BMDL_Error_Size = BMDL_func(nparm, pb, Bmdl[1], tol);
    }
 

  if (bmdlCurve == Yes)
    {
/****** calculate Bmd[] and Bmdl[] ***********/
      for (k = 2; k <= 5; k++)
	{         
	  for (j = 1; j <= nparm; j++)
	    p[j] = pBak[j];	//get the "old" p[].

	  if (k == 2)
	    Rlevel[k] = BMR = 0.05;
	  else
	    Rlevel[k] = BMR = (k - 2) * 0.1;

	  if (bmdparm.risk == 1)
	    {
	      if (Spec[1] == Spec[3])
		{
		  ck = -1.0 * log ((1 - pBak[1] * spfixed - pBak[3] * snfixed) / BMR - 1);
		}
	      else
		{
		  ck = -1.0 * log ((1 - pBak[1] - pBak[3] * sijfixed) / BMR - 1);
		}
	    }
	  else
	    ck = -log ((1 - BMR) / BMR);

/**** solve the BMD ********************************/
	  D = 0.0;
	  xb = exp ((ck - p[2] - p[4] * sijfixed) / p[5]);
	  Bmdjunk = xb;

/********* search for BMDL **************************/
	  xa = xb * stepsize;
	  tol = FMAX ((Bmdjunk) * 0.001, 0.0000001);
	  BMD_lk = xlk;		//get the lk at BMD.

	  fb = DBL_MAX;
	  for (i = 1; i <= nparm; i++) pa[i] = pb[i] = p[i];
	  fa = BMDL_func (nparm, pa, xa, tol);

	  while (fa < 0.0 && xa > DBL_MIN && fabs(fa - fb) > DBL_EPSILON)
	    {
	      xb = xa;
	      fb = fa;
	      for (i = 1; i <= nparm; i++) pb[i] = pa[i];
	      xa *= stepsize;
	      fa = BMDL_func (nparm, pa, xa, tol);
	    }
	  if (fa < 0.0)
	    {
	      Bmdl[k] = -1;
	      fprintf (fp_out, "\n BMDL curve computation failed for BMR = %f . \n The BMDL curve appearing in the graph may not be accurate.", BMR);
	    }
	  else
	    {
	      Bmdl[k] = zeroin(xa, xb, 1.0e-10, BMDL_func, nparm, pb, 1.0e-14);
	    }

	}

    }
  else
    for (k = 2; k <= 5; k++)
      Bmdl[k] = Rlevel[k] = -1;

  for (j = 1; j <= nparm; j++)   /* set phi back to original value */
    p[j] = porig[j];       

  FREE_DVECTOR (porig, 1, nparm);
  FREE_DVECTOR (pBak, 1, nparm);
  FREE_DVECTOR (pa, 1, nparm);
  FREE_DVECTOR (pb, 1, nparm);
}				/*end BENCHMD */


/*****************************************************************
* BMDL_func -- used to compute the values of functions BMDL_f (the
*              X^2 value) at the point D, given the parm p[] and
*              number of parm. Input parameters are in "internal" form.
*
*              This routine is called by zeroin().
*
*****************************************************************/
double
BMDL_func (int nparm, double p[], double D, double gtol)
{				// ck , BMD_lk and LR are calculated in Nlogist_BMD()

  double fD, xlk;
  int junk;

  tD = D;			//tD is global var. have to change before call MAX_lk().
  MAX_lk (nparm, p, gtol, &junk, &xlk);
  fD = BMD_lk - xlk - LR;

  return fD;
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
