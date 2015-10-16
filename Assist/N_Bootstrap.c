#include <stdlib.h>
#include <math.h>
#include "benchmark.h"
#include "specialfun.h"
#include "allo_memo.h"
#include "in_outfun.h"
#include "gsl/gsl_randist.h"

extern int Model_DF (int []);
extern void Predict(double [], double [], int, double [], double []);
int compare_doubles (const void *a, const void *b);
double gsl_cdf_beta_P (double x, double a, double b);


//#define BSDEBUG 1000		/* Bootstrapping debug code (N>0 outputs every N BSDebug Iterations) */

/***********************************************************
*N_Bootstrap -- Used to test the Goodness of Fit
*
************************************************************/
void
N_Bootstrap (int ngrp, int nparm, double Parms[], int bounded[],
	    int Nobs, double Xi[], double Yp[], double Yn[], double Ls[],
	    int Xg[], double SR[], int BSIter, long BSSeed )
{
  double *GYpp;			/*predicted positive depdendent values */
  double *GEp;			/*estimated probability */
  double *P;
  double *Gvar, *GYsum, *Gcount;
  int *GrpSize;			/*numbers of obs. in each group. */
  int i, j, k, l;
  double gfit;
  
  /*** bootstrap variables */
  const gsl_rng_type * type;    /* type for beta distribution */
  gsl_rng * r;                  /* random value for beta distribution */
  double urand;                 /* uniform random value [0,1]  */
  double a, b;                  /* values used in beta distribution function */
  double cutpoint;              /* cut-off probability for bootstrapping */
  double *GYp_new;              /* psuedo-observed values */
  double *GSR_new;              /* SR values based on new variables */
  double **GSR_newsqsum;	/* Sum of SR^2 values for new method */
  int SRsqCount;                /* tracks SR^2 >= SR^2_original */
  int BSLoops = 3;              /* number of times to repeat bootstrap method */
  double *pv;
  double pavg;
  double percentiles[] = {.50, .90, .95, .99};	/* percentile values to calculate */
  int perloc;
  double x, cum_beta, cum_beta_comp, cum_beta_comp_low, cum_beta_comp_hi;  /*cumulative beta variables for NaN method */
  double cum_beta_step = 0.00001;   /* step size for cumulative beta table */
  /*** end of new method variables */
    

  /* Note: have to redefine those var. because the data must be */
  /* changed in order to compute goodness-of-fit statistic. */
  GYpp = DVECTOR (1, Nobs);
  GEp = DVECTOR (1, Nobs);
  GYsum = DVECTOR(1, Nobs);
  Gcount = DVECTOR(1, Nobs);
  Gvar = DVECTOR(1, Nobs);
  GrpSize = IVECTOR (1, ngrp);
  P = DVECTOR (1, nparm);

  /*** for bootstrap method */
  GYp_new = DVECTOR (1, Nobs);  
  GSR_new = DVECTOR (1, Nobs);
  GSR_newsqsum  = DMATRIX (1, BSLoops, 1, BSIter);
  pv = DVECTOR (1,BSLoops);	
#ifdef BSDEBUG
  double *GSR;			/* Add capablility of storing SR values for debug output */   
  FILE *fp_bsdbg;  
  fp_bsdbg = fopen("BSDEBUG.log", "w"); 
  if (fp_bsdbg == NULL) 
    {
    printf ("Error in opening BSDEBUG file\n"); 
    printf ("...now exiting to system...\n"); 
    exit(1);
    }
#endif



/**** compute the GrpSize[] -- # of obs. in each group***/
  for (i = 1; i <= ngrp; i++)
    GrpSize[i] = 0;

  for (j = 1; j <= Nobs; j++)
    GrpSize[Xg[j]] += 1;


  for (j = 6; j <= nparm; j++)
    P[j - 5] = Parms[j]; 

/** compute the estimated prob. and expected obs.  *******/
  Predict(Xi, Ls, Nobs, Parms, GEp);

  

  for (i = 1; i <= Nobs; i++)
    {
      GYsum[i] = Yp[i] + Yn[i];
      GYpp[i] = GEp[i] * GYsum[i];
      Gvar[i] = GYsum[i] * GEp[i] * (1 - GEp[i]) * (1.0 + (GYsum[i] - 1.0) * P[Xg[i]]);
    }


#ifdef BSDEBUG
  /*** Add capability to tally # of observed at each value */
  int max = GYsum[1];
  for (i = 2; i <= Nobs; i++)
    {
      if (max < (int)GYsum[i])
        max = (int)GYsum[i];
    }
  int *Obs_Count;          	
  Obs_Count = IVECTOR(0, max);
#endif




  /*** compute the goodness-of-fit test for litter data. *********/
  gfit = 0.0;

  for (k = 1; k <= Nobs; k++)
    {
      SR[k] = (Yp[k] - GYpp[k])/sqrt(Gvar[k]);
      gfit += SR[k]*SR[k];
    }





  fprintf (fp_out,
#ifndef RBMDS
          "\n\n\n\n Observed Chi-square = %10.4f\n"
#else
          "\n\n\n\n Observed Chi-square = %30.22\n"
#endif
          , gfit);



/*** Bootstrap method for nested models  */

  /*** Set up random number generator for GSL */
  gsl_rng_env_setup();
  type = gsl_rng_default;
  r = gsl_rng_alloc(type);    
  gsl_rng_set (r, BSSeed);

  fprintf (fp_out, "\n \t\t  Bootstrapping Results  \n");
  fprintf (fp_out, "\nNumber of Bootstrap Iterations per run: %d\n", BSIter);
  fprintf (fp_out, "\n\t\t    Bootstrap Chi-square Percentiles");
  fprintf (fp_out, "\n Bootstrap");
  fprintf (fp_out, "\n    Run        P-value    50th     90th     95th     99th");
  fprintf (fp_out, "\n----------------------------------------------------------- ");


#ifdef BSDEBUG
  fprintf (fp_bsdbg, "\n\n ================================================================================================ \n");
  fprintf (fp_bsdbg, "\n \t\t\t\t Bootstrapping Debug Log\n");
  fprintf (fp_bsdbg, " ================================================================================================ ");
#endif

  for (l=1; l<=BSLoops; l++)   /*** bootstrap run loop */
    {
    SRsqCount = 0;
  

#ifdef BSDEBUG
    fprintf (fp_bsdbg, "\n ================================================================================================ ");
    fprintf (fp_bsdbg, "\n \t\t\t\t Bootstrapping Loop # %d\n", l);
    fprintf (fp_bsdbg, "\n GSL random number generator type = %s\n", gsl_rng_name (r));
    fprintf (fp_bsdbg, " Seed = %ld\n\n", BSSeed);
    fprintf (fp_bsdbg, " GSL beta distribution used to calculate Prob. Cutoff\n");
    fprintf (fp_bsdbg, " when dose correlation coefficents are non-zero\n");
    fprintf (fp_bsdbg, " ============================================================================================================================================================================================= ");
    fprintf (fp_bsdbg, "\n\n                                                           Litter       Actual     Pseudo      Actual         Pseudo");
    fprintf (fp_bsdbg, "\n    Dose       Est._Prob.      Cut_Prob.     Expected       Size       Observed   Observed       SR             SR             Phi                         a                        b");
    fprintf (fp_bsdbg, "\n ============================================================================================================================================================================================= ");
#endif

    for (i = 1; i <= BSIter; i++)    /*** bootstrap iteration loop */
      {

#ifdef BSDEBUG 
      if (i % BSDEBUG == 0)
        {
        fprintf(fp_bsdbg, "\n\n==================");
        fprintf(fp_bsdbg, "\n= Iteration %4d =", i);
        fprintf(fp_bsdbg, "\n==================\n");    
        }
      for (k=0; k<=max; k++) Obs_Count[k] = 0;
#endif

        GSR_newsqsum[l][i] = 0;   	/*** Set SR squared sum to zero for new calcs */
        for (j = 1; j <= Nobs; j++)	/*** Loop over each line in litter data */
          {
          GYp_new[j] = 0;   // reset pseudo-observed value to zero 
          
          /*** find cut-off probability (cutpoint) */
          if ((P[Xg[j]] == 0) || (P[Xg[j]] == 1))             // Use Est._Prob. for cut-off probability when phi = 0 or 1
            cutpoint = GEp[j];
          else
            {
            // find cut-off probability from beta distribution
            double tempvalue = GEp[j]*(1.0-GEp[j])/(GEp[j]*(1.0-GEp[j])*P[Xg[j]]) - 1.0;
            a = GEp[j]*tempvalue;
            b = (1.0 - GEp[j])*tempvalue; 
            cutpoint = gsl_ran_beta(r, a, b);
            }

          
          /*** use cumulative method if cutpoint = NaN */
          if (isnan(cutpoint))
              {
              cum_beta_comp_low = gsl_cdf_beta_P(0.00001, a, b);  // low end of cumulative scale
              cum_beta_comp_hi = gsl_cdf_beta_P(0.99999, a, b);   // hi end of cumulative scale
              cum_beta = gsl_rng_uniform(r);                  // random value to compare to cumulative scale
              
              if (cum_beta < cum_beta_comp_low)          // if random value is less than low end 
                cutpoint = 0;
              else if (cum_beta >= cum_beta_comp_hi)        // if random value is greater than high end 
                cutpoint = 1.0 - cum_beta_step;
              else if (fabs(cum_beta_comp_low - cum_beta)< fabs(cum_beta_comp_hi - cum_beta))   // if random value is closer to low end 
                {
                cum_beta_comp = 0;
                x = 0;
                while (x<=1.0 && cum_beta > cum_beta_comp)
                  {
                  cutpoint = x;
                  x += cum_beta_step;
                  cum_beta_comp = gsl_cdf_beta_P(x, a, b);
                  }
                }
              else              // if random value is closer to high end
                {
                cum_beta_comp = 1.0;
                x = 1;
                cutpoint = 1.0;
                while (x>=0.0 && cum_beta < cum_beta_comp)
                  {
                     x -= cum_beta_step;
                     cutpoint = x;
                     cum_beta_comp = gsl_cdf_beta_P(x, a, b);
                  }
                }
              }



          /*** Calculate pseudo-observed value */
          if (P[Xg[j]] != 1)
            {
            for (k = 1; k <= GYsum[j]; k++)
              {
              urand = gsl_rng_uniform (r);	/*generate uniform random number*/
              if (urand <= cutpoint) GYp_new[j]++;    /* Check random number against cutoff-probability.  */
              }
            }
          else
            {
            /*** phi = 1 method */
            urand = gsl_rng_uniform (r);
            if (urand <= cutpoint)
              GYp_new[j] = GYsum[j];
            else
              GYp_new[j] = 0;
            }


          /*** Calculate Scaled Residual based on pseudo-observation (GYp_new[]) */
          GSR_new[j] = Gvar[j] > 0 ? (GYp_new[j] - GYpp[j])/sqrt(Gvar[j]) : 0;
          GSR_newsqsum[l][i] += SQR(GSR_new[j]);         

#ifdef BSDEBUG
  if (i % BSDEBUG == 0)
    {
    Obs_Count[(int)GYp_new[j]]++; 
    //fprintf(fp_bsdbg, "%9.4f\t%6.3f  \t%6.3f   %5.0f      %7.3f\t%5.0f\t     %5.0f\t%9.4f\t%9.4f\n", Xi[j], GEp[j], cutpoint, GYpp[j], GYsum[j], Yp[j], GYp_new[j], SR[j], GSR_new[j]);
    fprintf(fp_bsdbg, "%9.4f       %7.5f        %8.5f      %7.3f      %7.3f      %5.0f      %5.0f      %9.4f      %9.4f      %18.14f      %18.14f      %18.14f\n", 
                      Xi[j], GEp[j], cutpoint, GYpp[j], GYsum[j], Yp[j], GYp_new[j], SR[j], GSR_new[j], P[Xg[j]], a, b);
    }
#endif
          }  

        if (GSR_newsqsum[l][i] >= gfit)
          SRsqCount++;
    
#ifdef BSDEBUG
  if (i % BSDEBUG == 0)
    {
    fprintf(fp_bsdbg, "\n\t\t\t\t\t\t\t\t\t\t\t       ========	      ========\n");
    fprintf(fp_bsdbg, "\t SRsqCount = %5d\t\t\t\t\t\t  Sum of SR^2:       %9.4f	     %9.4f\n\n", SRsqCount, gfit, GSR_newsqsum[l][i]);
    for(k=0; k<=max; k++)
      {
      fprintf(fp_bsdbg, " %d pseudo-obs of value %d\n", Obs_Count[k], k);
      }
    }
#endif

        }

    pv[l] = (double) SRsqCount/BSIter;	/* compute p-value for each BS loop */

  }
  
  /*** Calculate chi-square percentiles for each BS loop */
  pavg = 0;

  for (l=1; l<=BSLoops; l++)
    {
          

      fprintf(fp_out, "\n %5d          %6.4f  ", l, pv[l]);
            
      /*** compute percentile values of sum of SR^2 ***/
      qsort (&GSR_newsqsum[l][1], BSIter, sizeof(GSR_newsqsum[1][1]), compare_doubles);
      for (k=0; k<sizeof(percentiles)/sizeof(percentiles[0]); k++)
        {
    
        double temp = percentiles[k] * BSIter;
        if ((ceilf(temp)==temp) && (floorf(temp)==temp))
          {
          perloc = (int)temp;
          fprintf(fp_out, "%7.4f  ", GSR_newsqsum[l][perloc]);     
          }
        else
          {
          perloc = (int)ceilf(temp);
          fprintf(fp_out, "%7.4f  ", (GSR_newsqsum[l][perloc]+GSR_newsqsum[l][perloc+1])/2.0);
          }
        }


    pavg += pv[l];
    }


  
    /*** combined P-value ***/
    pavg /= 3.0;
    fprintf (fp_out, "\n----------------------------------------------------------- ");
    fprintf (fp_out, "\n Combined       %6.4f  ", pavg);

    /*** combined chi-square percentiles ***/
    qsort (&GSR_newsqsum[1][1], BSLoops*BSIter, sizeof(GSR_newsqsum[1][1]), compare_doubles);
    double *pGSRsum = &GSR_newsqsum[1][1];
    for (k=0; k<sizeof(percentiles)/sizeof(percentiles[0]); k++)
        {      
        double temp = percentiles[k] * BSLoops*BSIter;
        if ((ceilf(temp)==temp) && (floorf(temp)==temp))
          {
          perloc = (int)temp;
          fprintf(fp_out, "%7.4f  ", *(pGSRsum+perloc-1));     
          }
        else
          {
          perloc = (int)ceilf(temp);
          fprintf(fp_out, "%7.4f  ", (*(pGSRsum+perloc-1)+*(pGSRsum+perloc))/2.0);
          }
        }

    fprintf(fp_out, "\n\n\n\nThe results for three separate runs are shown.  If the estimated p-values are sufficiently");
    fprintf(fp_out, "\nstable (do not vary considerably from run to run), then then number of iterations is");
    fprintf(fp_out, "\nconsidered adequate.  The p-value that should be reported is the one that combines");
    fprintf(fp_out, "\nthe results of the three runs.  If sufficient stability is not evident (and especially");
    fprintf(fp_out, "\nif the p-values are close to the critical level for determining adequate fit, e.g., 0.05),");
    fprintf(fp_out, "\nthen the user should consider increasing the number of iterations per run.\n");
    
    


#ifdef BSDEBUG
  if (fclose(fp_bsdbg) != 0)
    printf("Error closing bootstrap file");
#endif

/*** End bootstrap method for nested models - Cody Simmons */


  FREE_DVECTOR (GEp, 1, Nobs);
  FREE_DVECTOR (GYsum, 1, Nobs);
  FREE_DVECTOR (Gcount, 1, Nobs);
  FREE_DVECTOR (Gvar, 1, Nobs);
  FREE_IVECTOR (GrpSize, 1, ngrp);
  FREE_IVECTOR (P, 1, nparm);
  gsl_rng_free(r);    /* free memory from random generator */
}


int compare_doubles (const void *a, const void *b)
{
  if(*(double*)a > *(double*)b) return 1;
  else if (*(double*)a < *(double*)b) return -1;
  else return 0;
}