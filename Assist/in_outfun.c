/****************
*  in_outfun.c
*  Nov. 22 1996
*****************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "benchmark.h"
#include "ERRORPRT.h"
#include "allo_memo.h"
#include "specialfun.h"
#include "matrix_agb.h"
#include "computation.h"
#include "dcdflib/cdflib.h"

/******************************************************
** READ_PARAMETERS--subfunction used to read parameters.
*******************************************************/
void READ_PARAMETERS (int nparm, double Parms[])
{
  int i;

  for (i = 1; i <= nparm; i++)
    fscanf (fp_in, "%lf", &Parms[i]);
}
/**********************************************************
** READ_OBSDATA5V--used to read 5 column data in 5 vectors.
***********************************************************/
int
  READ_OBSDATA5V (int Nobs, double Xi[], double Yp[], double Yn[],
		  double Ls[], int Xg[])

{
  int Nmiss;			/*number of records with missing values */
  int i, j, n, m;		/*count and iteration control variables */
  double value;			/*temp variable */

  Nmiss = 0;
  for (i = 1; i <= Nobs; i++)
    {
      n = i - Nmiss;
      m = 0;
      for (j = 1; j <= 5; j++)
	{
	  fscanf (fp_in, "%lf", &value);
	  if ((j < 5 && value != MISSING) || j == 5 )
	    {
	      if (j == 1)
		Xi[n] = value;
	      if (j == 2)
		Yp[n] = value;
	      if (j == 3)
		Yn[n] = value;
	      if (j == 4)
		Ls[n] = value;
	      if (j == 5)
		Xg[n] = (int) value;
	    }
	  else
	    m++;
	}
      if (m != 0)
	Nmiss++;
      else if ((Yp[n] + Yn[n]) <= 0 || Ls[n] <= 0)
	Nmiss++;
    }

  return Nmiss;
}


/**********************************************************
** READ_OBSDATA3V--used to read 3 column data in 3 vectors.
***********************************************************/
int
  READ_OBSDATA3V (int Nobs, int tcols, int pdcol, int ndcol, int idcol, double Yp[], double Yn[], double X[])
{
  int Nmiss;			/*number of records with missing values */
  int i, j, n, m;		/*count and iteration control variables */
  double value;			/*temp variable */

  Nmiss = 0;
  for (i = 1; i <= Nobs; i++)
    {
      n = i - Nmiss;
      m = 0;
      for (j = 1; j <= tcols; j++)
	{
	  fscanf (fp_in, "%lf", &value);
	  if (value != MISSING)
	    {
	      if (j == idcol)
		X[n] = value;
	      if (j == pdcol)
		Yp[n] = value;
	      if (j == ndcol)
		Yn[n] = value;
	    }
	  else
	    m++;
	}
      if (m != 0)
	Nmiss++;
      else if ((Yp[n] + Yn[n]) <= 0.0)
	Nmiss++;
    }

  return Nmiss;
}

/***********************************************
**  OUTPUT_TEXT--output model fitting text info.
************************************************/
void
  OUTPUT_TEXT (char txt[])
{
#ifdef MISC_OUT
  /*output to the screen */
  printf ("%s\n", txt);
#endif
  /*output to bmdswrk.out */
  fprintf (fp_out, "%s\n", txt);
}


/*********************************************************
** OUTPUT_DTMS3PARMS--output model fitting parameters for a
*                     dichotomous model.
*
*  Modified By: Geoffrey Nonato
*  Date: 12/20/2006
*  Description: Added print_SE plus code to print SE according to print_SEs value
*               (0 means no print, > 0 means print).  An asterisk will be printed instead
*               if print_SE = 0.  And footnote will be printed "Indicates that this value is not calculated."
**********************************************************/
void
  OUTPUT_DTMS3PARMS (int nparm, int Spec[], int bounded[], double Parms[],
		     char *tparms[], double **vcv, int print_SE)
{
	int i, jvar;

	/*output to *.out file */
	OUTPUT_TEXT ("\n\n\n                                Parameter Estimates");
	OUTPUT_TEXT (    "\n       Variable           Estimate             Std. Err. ");
	jvar = 0;
	for (i = 1; i <= nparm; i++)
	{
		if(print_SE > 0)
		{
			if (Spec[i] == 0)
			{
				if (bounded[i] == 0)
				{
					jvar++;
#ifdef MISC_OUT
					printf ("%15s%20.6g%20.6g\n", tparms[i - 1], Parms[i],
						sqrt (vcv[jvar][jvar]));
#endif
					if (vcv[jvar][jvar] > 0)
					{
						fprintf (fp_out, 
#ifndef RBMDS
							"%15s %19.6g %19.6g\n"
#else
							"%15s %30.22g %30.22g\n"
#endif
							, tparms[i - 1], Parms[i],
							sqrt (vcv[jvar][jvar]));
					}
					else
					{
						fprintf (fp_out, 
#ifndef RBMDS
							"%15s %19.6g %19s\n"
#else
							"%15s %30.22g %30s\n"
#endif
							, tparms[i - 1], Parms[i],
							"NA");
					}
				}
				else
				{
					fprintf (fp_out, 
#ifndef RBMDS
						"%15s %19.6g             Bounded\n"
#else
						"%15s %30.22g             Bounded\n"
#endif
						, tparms[i - 1], Parms[i]);
				}
			}
		}
		else //print_SE == 0
		{
			if (Spec[i] == 0)
			{
				if (bounded[i] == 0)
				{
					jvar++;
#ifdef MISC_OUT
					printf ("%15s%20.6g%20s\n", tparms[i - 1], Parms[i],
						"*");
#endif
					fprintf (fp_out, 
#ifndef RBMDS
						"%15s %19.6g %13s\n"
#else
						"%15s %30.22g %24s\n"
#endif
						, tparms[i - 1], Parms[i],
						"*");
				}
				else
				{
					fprintf (fp_out, 
#ifndef RBMDS
						"%15s %19.6g             *\n"
#else
						"%15s %30.22g             *\n"
#endif
						, tparms[i - 1], Parms[i]);
				}
			}
		}
	} //end of for(i = 1; i <= nparm; i++)
	if(print_SE == 0) /* Don't Print Standard Error */
		fprintf (fp_out, "\n* - Indicates that this value is not calculated.\n");
}

/*********************************************************
** OUTPUT_Init--output initial  parameter values.
*
**********************************************************/
void OUTPUT_Init (int nparm, int Spec[], double Parms[], char *tparms[])
{
	int i;
	/*output to screen and bmdswrk.002 tmp file */
	for (i = 1; i <= nparm; i++) {
		if (Spec[i] == 0)
			fprintf (fp_out, 
#ifndef RBMDS
				 "%31s = %12.6g\n"
#else
				 "%31s = %30.22g\n"
#endif
				 , tparms[i - 1], Parms[i]);
		else
			fprintf (fp_out, 
#ifndef RBMDS
				 "%31s = %12.6g   Specified\n"
#else
				 "%31s = %30.22g   Specified\n"
#endif
				 , tparms[i - 1], Parms[i]);
	} /* end for */
}



/*********************************************************************************
*
*  OUTPUT_DTMS3VCV--This function will take the matrix of Likelihood
*		    second derivatives and remove the columns/rows
*		    that correspond to any bounded or fixed parameters,
*		    invert this modified matrix, then output it.
*
*  INPUT:
*      nparm       number of parameters in the model
*      Spec[i]     0 if parameter i is not specified by user
*		   1 if parameter i is specified by user
*      parmtxt[i]  the ith parameters name
*      vcv[][]     contains the second partial derivatives of
*	           the likelihood function.  Element ij is
*		   the second partial with respect to
*		   parameters i and j
*      bounded[i]  0 if the parameter has not been estimated
*		   at a boundary point
*		   1 if the parameter has been estimated at a
*		   boundary point
*  INPUT/OUTPUT:
*	vcv_adj[][]-	Upon entry, is just an uninitialized matrix
*			of size MxM where M = # of non-bounded
*			parameters.  At return, it is the corr.
*			matrix with all bounded and specified
*			parameters removed
**********************************************************************************/
void
  Get_and_OUTPUT_DTMSVCV (int nparm, int Spec[], char *parmtxt[],
			  double **vcv, double **vcv_adj, int *bounded)
{
  int adj_vcv_rows, i, j, checkvar, counter, *notbounded, nboundparm;
  double cor, **temp_vcv;


  /* Create an adjusted variance-covariance matrix that does not
     /  have correlations/variances for parameters that are estimated
     /  at a boundary point */




  temp_vcv = DMATRIX (1, nparm, 1, nparm);


  for (i = 1; i <= nparm; i++)
    {
      for (j = 1; j <= nparm; j++)
	{
	  temp_vcv[i][j] = vcv[i][j];
	}
    }


  nboundparm = 0;
  for (i = 1; i <= nparm; i++)
    {
      if (bounded[i] == 1)
	{
	  nboundparm++;

	}			/* end if */

    }				/* end for */

  notbounded = IVECTOR (1, nparm - nboundparm);
  counter = 0;

  for (i = 1; i <= nparm; i++)
    {
      if (bounded[i] == 0)
	{
	  /* notbounded[j] contains the subscript in bounded that contains the
	     jth non-bounded parameter */

	  counter++;
	  notbounded[counter] = i;

	}			/* end if */

    }				/* end for */

  /* Adjust matrix to only contain non-fixed or bounded parameters */
  adj_vcv_rows = Take_Out_Bounded_Parms (nparm, bounded, vcv);

  for (i = 1; i <= adj_vcv_rows; i++)
    {
      for (j = 1; j <= adj_vcv_rows; j++)
	{
	  vcv_adj[i][j] = vcv[i][j];

	}			/* end for j */

    }				/* end for i */


  /* Invert to get the covariance matrix */

  INVMAT (vcv_adj, adj_vcv_rows);

  /*output to screen and bmdswrk.002 temp file */
  OUTPUT_TEXT ("\n\n           Asymptotic Correlation Matrix of Parameter Estimates\n");

  if (nboundparm > 0)
    {
      /* printf ("           ( *** The model parameter(s) "); */
      fprintf (fp_out, "           ( *** The model parameter(s) ");
      for (i = 1; i <= nparm; i++)
	{
	  if (bounded[i] == 1)
	    {
	      /* printf (" -%s   ", parmtxt[i - 1]); */
	      fprintf (fp_out, " -%s   ", parmtxt[i - 1]);
	    }
	}
      /* printf ("\n                 have been estimated at a boundary point, or have been specified by the user,"); */
      /* printf ("\n                 and do not appear in the correlation matrix )\n\n"); */
      fprintf (fp_out, "\n                 have been estimated at a boundary point, or have been specified by the user,");
      fprintf (fp_out, "\n                 and do not appear in the correlation matrix )\n\n");
    }

  fprintf (fp_out, "          ");

  counter = 0;

  for (i = 1; i <= nparm; i++)
    {
      if (bounded[i] == 0)
	{
	  /* printf ("%15s", parmtxt[i - 1]); */
	  fprintf (fp_out, "%13s", parmtxt[i - 1]);
	}
    }
  /* printf ("\n"); */
  fprintf (fp_out, "\n");


  for (i = 1; i <= adj_vcv_rows; i++)
    {
      /* printf ("\n%10s", parmtxt[notbounded[i] - 1]); */
      fprintf (fp_out, "\n%10s", parmtxt[notbounded[i] - 1]);
      for (j = 1; j <= adj_vcv_rows; j++)
	{
	  if (vcv_adj[i][i] <= 0.0 || vcv_adj[j][j] <= 0.0)
	    {
	      /* printf ("NA"); */
	      fprintf (fp_out, "      NA       ");
	      checkvar = 1;
	    }
	  else
	    {
	      cor = vcv_adj[i][j] /
		(sqrt (vcv_adj[i][i]) * sqrt (vcv_adj[j][j]));
	      /* printf ("%13.2g", cor); */
	      fprintf (fp_out, 
#ifndef RBMDS
		       "%13.2g", 
#else
		       "%30.22g",
#endif
		       cor);
	      checkvar = 0;
	    }
	}
      /* printf ("\n"); */
      fprintf (fp_out, "\n");
    }

  if (checkvar == 1)
    {
      fprintf (fp_out,
      "\n\nNA - This parameter's variance has been estimated as zero or less.\nTHE MODEL HAS PROBABLY NOT CONVERGED!!!\n");
    }

  for (i = 1; i <= nparm; i++)
    {
      for (j = 1; j <= nparm; j++)
	{
	  vcv[i][j] = temp_vcv[i][j];
	}
    }
    
    FREE_IVECTOR(notbounded,1,nparm - nboundparm);
    FREE_DMATRIX(temp_vcv,1,nparm,1,nparm);
}

/******************************************************
**  OUTPUT_DTMS3ANOVA--output ANOVA table values for a
*                      3-parameter dichotomous model.
*******************************************************/
void
  OUTPUT_DTMS3ANOVA (char *anatxt[], AnaList anasum[])
{
  int Nparms;

  /*output to screen */
  OUTPUT_TEXT ("\n\n\n                        Analysis of Deviance Table");
  OUTPUT_TEXT ("\n       Model      Log(likelihood)  # Param's  Deviance  Test d.f.   P-value");

  /*output to *.out file */
  Nparms = anasum[3].DF + 1; /* for full model == # dose groups */
  fprintf (fp_out, 
#ifndef RBMDS
	   "%15s %15.6g %9d\n"
#else
	   "%15s %30.22g %9d\n"
#endif
	   , anatxt[0], anasum[1].SS, Nparms);
  Nparms = Nparms - anasum[2].DF; /* # fitted parameters in fitted model */
  if (anasum[2].DF > 0)
    {
      if (anasum[2].TEST > .0001)
	{
	  fprintf (fp_out, 
#ifndef RBMDS
		   "%15s %15.6g %9d %13.6g %6d %15.4g\n"
#else
		   "%15s %30.22g %9d %30.22g %6d %30.22g\n"
#endif
		   , anatxt[1], anasum[2].SS, Nparms, anasum[2].MSE, anasum[2].DF,
		   anasum[2].TEST);
	}
      else
	{
	  fprintf (fp_out, 
#ifndef RBMDS
		   "%15s %15.6g %9d %13.6g %6d %19.8g\n"
#else
		   "%15s %30.22g %9d %30.22g %6d %30.22g\n"
#endif
		   , anatxt[1], anasum[2].SS, Nparms, anasum[2].MSE, anasum[2].DF,
		   anasum[2].TEST);
	}
    }
  else
    {
      fprintf (fp_out, 
#ifndef RBMDS
	       "%15s %15.6g %9d %13.6g %6d         NA\n"
#else
	       "%15s %30.22g %9d %30.22g %6d         NA\n"
#endif
	       , anatxt[1], anasum[2].SS, Nparms, anasum[2].MSE, anasum[2].DF);
    }
  Nparms = 1; /* reduced model has just 1 parameter */
  if (anasum[1].DF > 1)
    {
      if (anasum[3].TEST > .0001)
	{
	  fprintf (fp_out, 
#ifndef RBMDS
		   "%15s %15.6g %9d %13.6g %6d %15.4g\n"
#else
		   "%15s %30.22g %9d %30.22g %6d %30.22g\n"
#endif
		   , anatxt[2], anasum[3].SS, Nparms, anasum[3].MSE, anasum[3].DF,
		   anasum[3].TEST);
	}
      else
	{
	  fprintf (fp_out, 
#ifndef RBMDS
		   "%15s %15.6g %9d %13.6g %6d         <.0001\n"
#else
		   "%15s %30.22g %9d %30.22g %6d         <.0001\n"
#endif		   
		   , anatxt[2], anasum[3].SS, Nparms, anasum[3].MSE, anasum[3].DF);
	}
    }
  else
    {
      fprintf (fp_out, 
#ifndef RBMDS
	       "%15s %15.6g %9d %13.6g %6d         NA\n"
#else
	       "%15s %30.22g %9d %30.22g %6d         NA\n"
#endif
	       , anatxt[2], anasum[3].SS, Nparms, anasum[3].MSE, anasum[3].DF);
    }

  fprintf (fp_out,
#ifndef RBMDS
	   "\n%15s %15.6g\n"
#else
	   "\n%15s %30.22g\n"
#endif
	   , "AIC:",
	   -2*anasum[2].SS + 2*(1.0 + anasum[3].DF - anasum[2].DF));

}

/***************************************************
**  OUTPUT_BENCHMD--output specified benchmark dose.
****************************************************/
void
  OUTPUT_BENCHMD (int pdcol, double BMD)
{

  OUTPUT_TEXT ("\n\n   Benchmark Dose Computation\n");

  fprintf (fp_out, "Specified effect = %14.6g\n\n", bmdparm.effect);
  if (bmdparm.risk == 0)
    fprintf (fp_out, "Risk Type        =      Extra risk \n\n");
  else
    fprintf (fp_out, "Risk Type        =      Added risk \n\n");
  fprintf (fp_out, "Confidence level = %14.6g\n\n", bmdparm.level);
  fprintf (fp_out, 
#ifndef RBMDS
	   "             BMD = %14.6g\n\n"
#else
	   "             BMD = %30.22g\n\n"
#endif
	   , BMD);

}

/******************************************************************
**  OUTPUT_DTMS3ANOVAC--output likelihood values for continuous
*                       models that may be interesting to the user
*
* Modified By: Micheal Ferree
* Date: 25MAY05
* Purpose: Made to output all tests regardless of type value.
*
******************************************************************/
void
  OUTPUT_DTMS3ANOVAC (char *anatxt[], AnaList anasum[], int type)
{
  char *frmt = 
#ifndef RBMDS 
    "%15s %19.6f %12d %14.6f\n";
#else
    "%15s %30.22f %12d %30.22g\n";
#endif

  fprintf(fp_out, "\n\n\n Model Descriptions for likelihoods calculated\n\n");
  fprintf(fp_out, "\n Model A1:        Yij = Mu(i) + e(ij)\n");
  fprintf(fp_out, "           Var{e(ij)} = Sigma^2\n\n");
  fprintf(fp_out, " Model A2:        Yij = Mu(i) + e(ij)\n");
  fprintf(fp_out, "           Var{e(ij)} = Sigma(i)^2\n\n");
  fprintf(fp_out, " Model A3:        Yij = Mu(i) + e(ij)\n");
  if (type)
    fprintf(fp_out, "           Var{e(ij)} = exp(lalpha + rho*ln(Mu(i)))\n");
  else
    fprintf(fp_out, "           Var{e(ij)} = Sigma^2\n");
    
  fprintf(fp_out, "     Model A3 uses any fixed variance parameters that\n");
  fprintf(fp_out, "     were specified by the user\n\n");
  fprintf(fp_out, " Model  R:         Yi = Mu + e(i)\n");
  fprintf(fp_out, "            Var{e(i)} = Sigma^2\n\n");

  OUTPUT_TEXT ("\n                       Likelihoods of Interest");
  OUTPUT_TEXT ("\n            Model      Log(likelihood)   # Param's      AIC");


  fprintf (fp_out, frmt, "A1", anasum[1].SS, anasum[1].DF, -2*(anasum[1].SS - anasum[1].DF));
  fprintf (fp_out, frmt, "A2", anasum[2].SS, anasum[2].DF, -2*(anasum[2].SS - anasum[2].DF));

  fprintf (fp_out, frmt, "A3", anasum[3].SS, anasum[3].DF, -2*(anasum[3].SS - anasum[3].DF));
  fprintf (fp_out, frmt, "fitted", anasum[5].SS, anasum[5].DF, -2*(anasum[5].SS - anasum[5].DF));
  fprintf (fp_out, frmt, "R", anasum[4].SS, anasum[4].DF, -2*(anasum[4].SS - anasum[4].DF));


}

/*********************************************************
** OP_ParmsE--output model fitting prameters for a
*                    continuous model.
*	INPUT:  nparm is the number of model parameters
*			Spec[] is a vector of zeros and ones, where
*				1 indicates that a parameter is fixed
*			Parms[] is the vector of parameter values
*			tparms[] is the vector of parameter names
*			vcv[][] is the variance-covariance matrix (adjusted
*				for parameters estimated at a boundary point)
*			bounded[] is a vector of zeros and ones.  If
*				bounded[i] = 1, then that parameter has
*				been estimated at some lower or upper bound
*  Modified By: Geoffrey Nonato
*  Date: 12/18/2006
*  Description: Added print_SE  as the last parameter plus code to print SE according to print_SEs value
*               (0 means no print, > 0 means print).  An asterisk will be printed instead
*               if print_SE = 0.  And footnote will be printed "Indicates that this value is not calculated."
**********************************************************/
void
  OP_ParmsE (int nparm, int Spec[], double Parms[], char *tparms[], double **vcv,
	     int *bounded, double conf, int print_SE)
{
	int i, num_bound, print_bound, pos_def;
	double zscore, P, Q, bound, mean = 0, sd = 1;
	int which = 2, status;
	P = 1-(1-conf)/2;
	Q = 1 - P;

	/* Check to make sure variances are all positive; if not, only output estimates, with a stern warning */
	pos_def = 1;
	num_bound = 0;
	for (i = 1; i <= nparm; i++) {
		if (bounded[i] == 0) {
			if (vcv[i - num_bound][i - num_bound] <= 0) {
				pos_def = 0;
				break;
			}
		} else {
			num_bound++;
		}
	}
	OUTPUT_TEXT ("\n\n\n                                 Parameter Estimates");
	fprintf (fp_out,"\n                                                         %4.1f%% Wald Confidence Interval",conf*100);
	OUTPUT_TEXT ("\n       Variable         Estimate        Std. Err.     Lower Conf. Limit   Upper Conf. Limit");
	
	if(print_SE > 0) /* Print Standard Error */
	{
		if (pos_def == 1) 
		{
			cdfnor(&which, &P, &Q, &zscore, &mean, &sd, &status, &bound);

			/*output to screen and Poly.out tmp file */

			num_bound = 0;
			print_bound = 0;

			for (i = 1; i <= nparm; i++)
			{
				if (bounded[i] == 0)
				{
					fprintf (fp_out, 
#ifndef RBMDS
						"%15s%17.6g%17.6g%20.6g%20.6g\n", 
#else
						"%15s %30.22g %30.22g %30.22g %30.22g\n", 
#endif
						tparms[i - 1], Parms[i],
						sqrt (fabs (vcv[i - num_bound][i - num_bound])),
						Parms[i]-zscore*sqrt (fabs (vcv[i - num_bound][i - num_bound])),
						Parms[i]+zscore*sqrt (fabs (vcv[i - num_bound][i - num_bound])));
				}
				else
				{
					if (Spec[i] == 0)
					{
						fprintf (fp_out, 
#ifndef RBMDS
							"%15s%17.6g               NA\n", 
#else
							"%15s%30.22g               NA\n", 

#endif
							tparms[i - 1], Parms[i]);
#ifdef MISC_OUT
						printf ("%15s%20.6g               NA\n", tparms[i - 1], Parms[i]);
#endif
						num_bound++;
						print_bound = 1;
					}
					else
					{
						num_bound++;
					}
				}
				/*fprintf (fp_out, "%15s%20.6g\n",tparms[i-1],Parms[i]);
				printf ("%15s%20.6g\n",tparms[i-1],Parms[i]); */
			}
		} 
		else 
		{
			num_bound = 0;
			print_bound = 0;

			for (i = 1; i <= nparm; i++)
			{
				if (bounded[i] == 0)
				{
					fprintf (fp_out, 
#ifndef RBMDS
						"%15s%17.6g%17s%20s%20s\n", 
#else
						"%15s %30.22g %30s %30s %30s\n", 
#endif
						tparms[i - 1], Parms[i],
						"NA",
						"NA",
						"NA");
				}
				else
				{
					if (Spec[i] == 0)
					{
						fprintf (fp_out, 
#ifndef RBMDS
							"%15s%17.6g               NA\n", 
#else
							"%15s%30.22g               NA\n", 

#endif
							tparms[i - 1], Parms[i]);
#ifdef MISC_OUT
						printf ("%15s%20.6g               NA\n", tparms[i - 1], Parms[i]);
#endif
						num_bound++;
						print_bound = 1;
					}
					else
					{
						num_bound++;
					}
				}
				/*fprintf (fp_out, "%15s%20.6g\n",tparms[i-1],Parms[i]);
				printf ("%15s%20.6g\n",tparms[i-1],Parms[i]); */
			}
		}
	}
	else /* Don't Print Standard Error */
	{
		if (pos_def == 1) 
		{
			cdfnor(&which, &P, &Q, &zscore, &mean, &sd, &status, &bound);

			/*output to screen and Poly.out tmp file */

			num_bound = 0;
			print_bound = 0;

			for (i = 1; i <= nparm; i++)
			{
				if (bounded[i] == 0)
				{
					fprintf (fp_out, 
#ifndef RBMDS
						"%15s%17.6g%13s%17s%19s\n", 
#else
						"%15s %30.22g %18s %22s %24s\n", 
#endif
						tparms[i - 1], Parms[i],
						"*",
						"*",
						"*");
				}
				else
				{
					if (Spec[i] == 0)
					{
						fprintf (fp_out, 
#ifndef RBMDS
							"%15s%17.6g            *                *                  *\n", 
#else
							"%15s%30.22g            *                *                  *\n", 

#endif
							tparms[i - 1], Parms[i]);
#ifdef MISC_OUT
						printf ("%15s%20.6g            *                *                  *\n", tparms[i - 1], Parms[i]);
#endif
						num_bound++;
						print_bound = 1;
					}
					else
					{
						num_bound++;
					}
				}
				/*fprintf (fp_out, "%15s%20.6g\n",tparms[i-1],Parms[i]);
				printf ("%15s%20.6g\n",tparms[i-1],Parms[i]); */
			}
		} 
		else 
		{
			num_bound = 0;
			print_bound = 0;

			for (i = 1; i <= nparm; i++)
			{
				if (bounded[i] == 0)
				{
					fprintf (fp_out, 
#ifndef RBMDS
						"%15s%17.6g%13s%17s%19s\n", 
#else
						"%15s %30.22g %18s %22s %24s\n", 
#endif
						tparms[i - 1], Parms[i],
						"*",
						"*",
						"*");
				}
				else
				{
					if (Spec[i] == 0)
					{
						fprintf (fp_out, 
#ifndef RBMDS
							"%15s%17.6g            *                *                  *\n", 
#else
							"%15s%30.22g            *                *                  *\n", 
#endif
							tparms[i - 1], Parms[i]);
#ifdef MISC_OUT
						printf ("%15s%20.6g            *                *                  *\n", tparms[i - 1], Parms[i]);
#endif
						num_bound++;
						print_bound = 1;
					}
					else
					{
						num_bound++;
					}
				}
				/*fprintf (fp_out, "%15s%20.6g\n",tparms[i-1],Parms[i]);
				printf ("%15s%20.6g\n",tparms[i-1],Parms[i]); */
			}
		}
	}

	if (print_SE > 0 && (pos_def == 1) && (print_bound > 0))
	{
		fprintf (fp_out, "\nNA - Indicates that this parameter has hit a bound\n");
		fprintf (fp_out, "     implied by some inequality constraint and thus\n");
		fprintf (fp_out, "     has no standard error.\n");
	}
	
	if(print_SE == 0) /* Don't Print Standard Error */
		fprintf (fp_out, "\n* - Indicates that this value is not calculated.\n");

	if (pos_def == 0) {
		fprintf(fp_out, "\nAt least some variance estimates are negative.\nTHIS USUALLY MEANS THE MODEL HAS NOT CONVERGED!\nTry again from another starting point.\n");
	}
}				/* end OP_ParmsE */



/******************************************************************
*  Output_Header -- output the header at the top of the output file
******************************************************************/

void
  Output_Header (char *version, char *input_name, char *plot_file_name, char *clocktime, char *note)
{

  fprintf (fp_out, "\n\n ==================================================================== \n");
  fprintf (fp_out, "   \t  %s", version);
  fprintf (fp_out, " \n  \t  Input Data File: %s  ", input_name);
  fprintf (fp_out, "\n  \t  Gnuplot Plotting File:  %s", plot_file_name);
  fprintf (fp_out, "\n \t\t\t\t\t\t\t%s", clocktime);
  fprintf (fp_out, " ==================================================================== ");
  fprintf (fp_out, "\n\n %s \n", note);
  fprintf (fp_out, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n ");

}				/* end Output_Header */


/****************************************************************** */
/* do_dmngb_warning -- print out a warning or error message based on */
/*                     the return code from dmngb, and modify */
/*                     eflag accordingly */
/****************************************************************** */

void do_dmngb_warning( int * eflag)
{
  if (*eflag > 0) Warning("-- NOTE ON MODEL CONVERGENCE:\n");
  if (*eflag == 9)
    {
      Warning("  The allowed number of iterations has been reached with no convergence.\n\
  Try restarting and allowing more iterations.");
      Warning("  Sometimes this is due to too-stringent convergence criteria.");
      Warning("  If the graph indicates a good fit, try increasing the relative");
      Warning("  function convergence criterion, for example, from 2e-16 to");
      Warning("  2e-08, or increasing the parameter convergence parameter.");
      Warning("\n  The tables below indicate the parameter and likelihood");
      Warning("  values that correspond to the graph.");
    } /* end if */
  if (*eflag == -1)
    Warning("-- All parameters have been specified.\n");
  if (*eflag != 0 && *eflag != 9 && *eflag != -1 && *eflag != 7)
    {
      fprintf(fp_out, "\n\n  ErrorFlag is: %d\n",*eflag);
      Warning("  Error: model failed to converge.\n\
Try restarting with new initial values.\n");
      Warning("\n  Tip: Go into Advanced User Mode (on the Options Menu)");
      Warning("  of BMDS).  Specify all the parameter values to be the same");
      Warning("  as the initial values listed in the output from this run.");
      Warning("  Rerun the model, and look at the graph.  Then, repeatedly");
      Warning("  adjust the parameters (mostly slope and power) until the");
      Warning("  curve comes close to the data.  Then, in the next try,");
      Warning("  switch all the \"specified\" labels to \"initialized\", and");
      Warning("  rerun the model.");
    }
  if (*eflag == 7)
    {
      Warning("  The model appears to have converged, in the sense that a");
      Warning("  reasonable change in the parameters is unlikely to increase");
      Warning("  the likelihood any further.  However, the hessian matrix at");
      Warning("  the reported solution appears to be singular, so this may not");
      Warning("  be a true maximum.");
      Warning("     Try restarting from different initial values.  Also,");
      Warning("  check the goodness of fit table and the graph of the function");
      Warning("  before accepting the results of this analysis.  This result can");
      Warning("  happen if the parameters are not uniquely determined, so that,");
      Warning("  for example, a continuum of values for (Slope, Power) in the Gamma");
      Warning("  model yield the same likelihood value.\n");
      *eflag = 0; /* So everything else will work */
    }
}

/*******************************************************************
 **PrintData -- Prints a table of observed values and expected values
 *              for continuous models.
 *********************************************************************/

void PrintData(double *mean, double *std, double *Xi, double *Ym,
	       double *Yd, int *Ni, int Nobs)
{
  int i;
  double obsstd, chi2_res, std0;

  fprintf(fp_out, "\n\n\n     Table of Data and Estimated Values of Interest\n\n");
  fprintf(fp_out, " Dose       N    Obs Mean     Est Mean   Obs Std Dev  Est Std Dev   Scaled Res.\n");
  fprintf(fp_out, "------     ---   --------     --------   -----------  -----------   ----------\n\n");


  for(i = 1; i <= Nobs; i++)
    {
      obsstd = sqrt(Yd[i]);
      std0 = std[i];
      if (std0 < 1e-20)
	std0 = .00000001;
      chi2_res = sqrt(Ni[i])*(Ym[i] - mean[i])/std0;
#ifdef MISC_OUT
      printf("%5.4g%6d%11.3g%13.3g%13.3g%13.3g%15.3g\n",
	     Xi[i], Ni[i], Ym[i], mean[i], obsstd, std[i], chi2_res);
#endif
      fprintf(fp_out,
#ifndef RBMDS
	      "%5.4g%6d%11.3g%13.3g%13.3g%13.3g%15.3g\n",
#else
	      "%5.4g %6d %11.3g %30.22g %13.3g %30.22g %30.22g\n",

#endif
	      Xi[i], Ni[i], Ym[i], mean[i], obsstd, std[i], chi2_res);

    }

}

