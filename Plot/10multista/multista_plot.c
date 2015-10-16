/******************************************************************
**
*  IMPORTANT NOTE:  If significant changes are being added to this
*					plotting program, please change the version
*					number in the corresponding model code
*
******************************************************************/

/*********************************************************************
*
*  Multista_plot.c - an ANSI C program that reads the output file 
*                   from Multista.C and then produces a script file for
*					GNU plot.
*  April, 1998
*
**********************************************************************/

#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>

#include <benchmark.h>
#include <specialfun.h>
#include <ERRORPRT.h>

#define  float double
#define  NR_END 1  

char fout[FLENGTH];         /* output file */
FILE  *fp_in, *fp_out; /* file pointers */

void CLOSE_FILES (void)
{
  if (fclose (fp_in) != 0 || fclose (fp_out) != 0 )
    ERRORPRT ("Error in closing opened files.");
}

int main(int argc, char *argv[])
{
  double *DVECTOR(int n1, int n2);
  void FREE_DVECTOR (double *v, int n1, int n2);

  char name1[CNLENGTH];  /* contains the string: BMD_flag */
  char name2[CNLENGTH];  /* contains the string: Nosb */
  char name3[CNLENGTH];  /* contains the string: nparm */
  char name4[CNLENGTH];  /* contains the string: Con_lev */
  char name4a[CNLENGTH];  /* contains the string: RiskType */
  char name4b[CNLENGTH];  /* contains the string: Effect */
  char name8[CNLENGTH];  /* contains the string: Data */
  char name9[CNLENGTH];  /* contains the string: Max_Min_dose */
  char name10[CNLENGTH]; /* contains the string: BMDL_comput_ind */
  char name11[CNLENGTH]; /* contains the string: RSL */
  char name12[CNLENGTH]; /* contains the string: BMD */
  char name13[CNLENGTH]; /* contains the string: BMDL */
  char name14[CNLENGTH]; /* contains the string: BMD_line */ 
  char name15[CNLENGTH]; /* contains the string: BMDL_line */
  char name16[CNLENGTH]; /* contains the string: BMDL_Curve_flag */
  char name17[CNLENGTH]; /* contains the string: smooth_opt */
  char name18[CNLENGTH]; /* contains the string: BMDL_curve */
  char name19[CNLENGTH]; /* contains the string: Check_result */
  char name[CNLENGTH];   /* variable for reading in parameter names */

  int flag1;       /* flag for computing BMD, 1 = yes */
  int flag2;       /* flag for computing BMDL, 1 = yes */
  int flag3;       /* flag for BMDL curve, 1 = yes */
  int smooth;      /* flag for smooth option, 1=unique, 0=C-spline */
  int nparm;       /* number of parameters */
  int nobs;        /* number of observations */
  int result;      /* flag for Check_result, 1 = yes */
  int i, icount, k; /* iteration variables */
  double value;    /* temp variable for storing parameter values */
  double conf;     /* confidence level, .95 is default */
  int rtype;    /* risk type: 0=extra, 1=added */
  double effect = 0;		/* BMR specified effect */
  double xmin;     /* min dose level, can change to smallest "x" value */ 
  double xmax;     /* max dose level, can change to largest "x" value */
  double ymin;     /* min of all responses and lower confidence limits */
  double ymax;     /* max of all responses and upper confidence limits */
  double rsl;      /* rsl = BMR*back1+back, see Multistage.c program */
  double bmd;      /* BMD */
  double bmdl;     /* BMDL */
  double xleft;    /* xmin - xrange/10 */ 
  double xright;   /* xmax + xrange/10 */
  double ybottom;  /* ymin - yrange/10 */
  double ytop;     /* ymax + yrange/10 */
  double xrange;   /* xmax - xmin */ 
  double yrange;   /* ymax - ymin */
  double xlabel;   /* xmin - xrange/12 */
  double ylabel;   /* ymin - yrange/15 */
  double tdose;    /* temp variable for dose */ 
  double tmean;    /* temp variable for response */
  double tmlb;     /* temp variable for lower confidence limit */
  double tmup;     /* temp variable for upper confidence limit */
  double bmd11;    /* xmin - xmax/100 */ 
  double bmd12;    /* same as rsl */
  double bmd21;    /* BMD */
  double bmd22;    /* rsl */
  double bmd31;    /* BMD */
  double bmd32;    /* -0.1 */
  double bmdL11;   /* BMDL */ 
  double bmdL12;   /* -0.1 */
  double bmdL21;   /* BMDL */
  double bmdL22;   /* rsl */
  double *dose;    /* vector for dose levels */
  double *mean;    /* vector of responses */
  double *mlb;     /* vector of lower confidence limits */
  double *mup;     /* vector of upper confidence limits */
  double *ldose;   /* vector of BMR's */
  double *lmean;   /* vector of BMDL's */
  double *beta;    /* vector of parameters */
  char long_path_name[FLENGTH];
  char *extpoint;
  char *pcRiskType = NULL;	/* Risk type for plot title */

  /* allocate memory for arrays */
  ldose=DVECTOR(1, 6);
  lmean=DVECTOR(1, 6);
  /* Allow spaces in path name */
  if (argc > 2)
    {
      path_name2(argc, argv, long_path_name);
      argv[1] = long_path_name;
    } /* end if */
  /* open the input file (has a 002 extension) */
  fp_in=fopen(argv[1], "r");
  /* copy argv[1] to fout */
  strcpy(fout, argv[1]);

  extpoint = strrchr(fout,(int) '.');
  if (extpoint != (char *) NULL)
    {
      *extpoint = '\0';
    }
  strcat(fout, ".plt");
 
  /* open the output file */
  fp_out=fopen(fout, "w");
  if (fp_in==NULL || fp_out==NULL)
    {
      printf("Error in opening input and/or output files. \n");
      printf("...now exiting the system... \n");
      fprintf(fp_out, "Error in opening input and/or output files. \n");
      fprintf(fp_out, "...now exiting the system... \n");
      exit(1);
    } /* end if */

  /*** Read the input data ***/
  fscanf(fp_in, "%s %d", name1, &flag1);
  fscanf(fp_in, "%s %d", name2, &nobs);
  fscanf(fp_in, "%s %d", name3, &nparm);
  fscanf(fp_in, "%s %lg", name4, &conf);
  fscanf(fp_in, "%s %d", name4a, &rtype);
  fscanf(fp_in, "%s %lg", name4b, &effect);
  switch (rtype) {
  case 0:
    pcRiskType = "Extra Risk";
    break;
  case 1:
    pcRiskType = "Added Risk";
    break;

  /* Future risk types would go here before the default case */
  default:
    pcRiskType = "Unknown Risk Type";
    break;
  } /* end switch */

  beta = DVECTOR(0,nparm-1);

  for (k=0; k<=nparm-1; k++)
    { 
      fscanf(fp_in, "%s %lg", name, &value);
      beta[k] = value;
    } /* end for */
  fscanf(fp_in, "%s", name8);
  /* allocate memory */
  dose=DVECTOR(1, nobs);
  mean=DVECTOR(1, nobs);
  mlb=DVECTOR(1, nobs);
  mup=DVECTOR(1, nobs);

  /** Find the max and min of x and y values **/
  ymin=1.0e+66;
  ymax=-1.0e+66;
  for (i=1; i<=nobs; i++)
    {
      fscanf(fp_in, "%lg %lg %lg %lg", &tdose, &tmean, &tmlb, &tmup);
      dose[i]=tdose;
      mean[i]=tmean;
      mlb[i]=tmlb;
      mup[i]=tmup;
      if (tmean < ymin)  ymin=tmean;
      if (tmlb < ymin)   ymin=tmlb;
      if (tmean > ymax)  ymax=tmean;
      if (tmup > ymax)   ymax=tmup;
    } /* end for */

  fscanf(fp_in, "%s", name9);
  fscanf(fp_in, "%lg %lg", &xmax, &xmin);

  if (flag1 == 1) /* if BMD calculation was selected */
    {
      fscanf(fp_in, "%s %lg", name11, &rsl);
      fscanf(fp_in, "%s %lg", name12, &bmd);
      fscanf(fp_in, "%s", name14);
      fscanf(fp_in, "%lg %lg", &bmd11, &bmd12);
      fscanf(fp_in, "%lg %lg", &bmd21, &bmd22);
      fscanf(fp_in, "%lg %lg", &bmd31, &bmd32);
    } /* end if */
  fscanf(fp_in, "%s %d", name10, &flag2);
  if (flag2 == 1) /* if BMDL was computed */
    {
      fscanf(fp_in, "%s %lg", name13, &bmdl);
      if (rsl > ymax) ymax=rsl;
      if (rsl < ymin) ymin=rsl;
      if (bmd > xmax) xmax=bmd;
      if (bmdl < xmin) xmin=bmdl;
      fscanf(fp_in, "%s", name15);
      fscanf(fp_in, "%lg %lg", &bmdL11, &bmdL12);
      fscanf(fp_in, "%lg %lg", &bmdL21, &bmdL22);
      fscanf(fp_in, "%s %d", name16, &flag3);
    } /* end if */

  if (flag3 == 1) /* if all BMDL's were calculated */
    {
      fscanf(fp_in, "%s %d", name17, &smooth);
      fscanf(fp_in, "%s", name18);
      icount=0;
      for (i=1; i<=6; i++)
	{
	  fscanf(fp_in, "%lg %lg", &tdose, &tmean);
	  if (tdose >= 0)
	    {
	      icount++;
	      ldose[icount]=tdose;
	      lmean[icount]=tmean;
	      if (tdose < xmin) xmin=tdose;
	      if (tdose > xmax) xmax=tdose;
	      if (tmean > ymax) ymax=tmean;
	      if (tmean < ymin) ymin=tmean;
	    } /* end if */
	} /* end for */
	
      fscanf(fp_in, "%s %d", name19, &result);
    } /* end  if (flag3 == 1) */
    
  /*** Start to produce the output file ***/

  if (argv[2] == NULL)
    fprintf(fp_out, 
#ifndef TERM_X11
	    "set terminal windows dashed\n"
#else
	    "set terminal x11\n"
#endif
	    );
  else
    {
      fprintf(fp_out, "set terminal postscript monochrome \n");
      fprintf(fp_out, "set output \"Multista.ps\" \n");
    } /* end else */
  fprintf(fp_out, "reset\n");
  fprintf(fp_out, "set time \"\%%H:\%%M \%%m/\%%d \%%Y\" \n");
  fprintf(fp_out, "set bar 3 \n");
  fprintf(fp_out, "set style line 6 linecolor rgb \"black\"\n");
  fprintf(fp_out, "set style line 16 linecolor rgb \"forest-green\" pt 12\n");
  fprintf(fp_out, "set key top left \n");
  fprintf(fp_out, "set xlabel 'dose' \n");
  fprintf(fp_out, "set ylabel 'Fraction Affected' \n");
  fprintf(fp_out, "set mxtics 10 \n");
  fprintf(fp_out, "set mytics 10 \n");
  if (flag1 == 1)
    {
      fprintf(fp_out, "set title 'Multistage Model, with BMR of %lg%% %s for the BMD and %lg Lower Confidence Limit for the BMDL'\n",
	      effect*100., pcRiskType, conf);
      fprintf(fp_out, "rl = %g \n", rsl);
      fprintf(fp_out, "bmd = %g \n", bmd);
    } /* end if */
  if (flag2 == 1)
    fprintf(fp_out, "bmdl = %g \n", bmdl);
  else
    fprintf(fp_out, "set title 'Multistage Model' \n");

  for (k=0; k<=nparm-1; k++)
    fprintf(fp_out, "beta%d = %g \n", k, beta[k]);

  fprintf(fp_out, "f(x)= beta0 + (x>0?( (1- beta0)*(1 - exp(");
  for (k=1; k<=nparm-1; k++)
    fprintf(fp_out, "-beta%d*(x**%d)", k, k);
  fprintf(fp_out, "))) : 0) \n"); 

  xrange=xmax-xmin;
  yrange=ymax-ymin;
  xleft=xmin-xrange/10.0;
  xright=xmax+xrange/10.0;
  ybottom=ymin-yrange/10.0;
  ytop=ymax+yrange/10.0;
  xlabel=xmin-xrange/12.0;
  ylabel=ymin-yrange/15.0;
  fprintf(fp_out, "set offsets %g, %g, 0, 0\n",xrange/20.0, xrange/20.0);
  if (flag1 == 1)
    {
      fprintf(fp_out, "set label 'BMD' at bmd, %g left \n", ylabel);
    }  /* end if */
  if (flag2 ==1)
    {
      fprintf(fp_out, "set label 'BMDL' at bmdl, %g right \n", ylabel);
      fprintf(fp_out, "set label '   ' at %g, rl left \n", xlabel);
    } /* end if */

  if (flag1 == 1 && flag2 == 1 && flag3 == 1)
    {
      fprintf(fp_out, "plot [x=%f:%f] [%f:%f] f(x) title 'Multistage',\\", 
	      xmin>0?0:xmin, xmax, ybottom, ytop);
      fprintf(fp_out, "\n     '-' using 1:2:3:4 notitle with errorbars ls 16,\\");
      fprintf(fp_out, "\n     '-' using 1:2 notitle with lines ls 6,\\");
      fprintf(fp_out, "\n     '-' using 1:2 notitle with lines ls 6,\\");
      if (smooth == 1)
	fprintf(fp_out, "\n     '-' using 1:2 smooth csplines title 'BMD Lower Bound' with lines ls 3,\\");
      else
	fprintf(fp_out, "\n     '-' using 1:2 smooth unique title 'BMD Lower Bound' with lines ls 3,\\");
      fprintf(fp_out, "\n     '-' using 1:2 notitle with points");
   
    } /* end if */

  else 
    if (flag1 == 1 && flag2 == 1)
      {
	fprintf(fp_out, "plot [x=%f:%f] [%f:%f] f(x) title 'Multistage',\\", 
		xmin>0?0:xmin, xmax, ybottom, ytop);
	fprintf(fp_out, "\n     '-' using 1:2:3:4 notitle with errorbars ls 16,\\");
	fprintf(fp_out, "\n     '-' using 1:2 notitle with lines ls 6,\\");
	fprintf(fp_out, "\n     '-' using 1:2 notitle with lines ls 6");
      } /* end if */

    else
      if (flag1 == 1)
	{
	  fprintf(fp_out, "plot [x=%f:%f] [%f:%f] f(x) title 'Multistage',\\", 
		  xmin, xmax, ybottom, ytop);
	  fprintf(fp_out, "\n     '-' using 1:2:3:4 notitle with errorbars ls 16,\\");
	  fprintf(fp_out, "\n     '-' using 1:2 notitle with lines ls 6");
	} /* end if */

      else
	{
	  fprintf(fp_out, "plot [x=%f:%f] [%f:%f] f(x) title 'Multistage',\\", 
		  xmin, xmax, ybottom, ytop);
	  fprintf(fp_out, "\n     '-' using 1:2:3:4 notitle with errorbars ls 16");
	} /* end else */

  for (i=1; i<=nobs; i++)
    fprintf(fp_out, "\n %f %f %f %f", dose[i], mean[i], mlb[i], mup[i]);
  fprintf(fp_out, "\n e");

  bmd11=xleft;   /* to adjust the positions of the BMD and BMDL lines */
  bmd32=bmdL12=ybottom;

  if (flag1 == 1) 
    {
      fprintf(fp_out, "\n %f %f", bmd11, bmd12);
      fprintf(fp_out, "\n %f %f", bmd21, bmd22);
      fprintf(fp_out, "\n %f %f", bmd31, bmd32);
      fprintf(fp_out, "\n e");
    } /* end if */
  if (flag1 == 1 && flag2 == 1)
    {
      fprintf(fp_out, "\n %f %f", bmdL11, bmdL12);
      fprintf(fp_out, "\n %f %f", bmdL21, bmdL22);
      fprintf(fp_out, "\n e");
    } /* end if */

  if (flag1 == 1 && flag2 == 1 && flag3 == 1)
    {
      for (i=1; i<=icount; i++)
	fprintf(fp_out, "\n %f %f", ldose[i], lmean[i]);
      fprintf(fp_out, "\n e");
      for (i=1; i<=icount; i++)
	fprintf(fp_out, "\n %f %f", ldose[i], lmean[i]);
      fprintf(fp_out, "\n e");
    } /* end if */

  /* free memory */
  FREE_DVECTOR(dose, 1, nobs);
  FREE_DVECTOR(mean, 1, nobs);
  FREE_DVECTOR(mlb, 1, nobs);
  FREE_DVECTOR(mup, 1, nobs);
  FREE_DVECTOR(ldose, 1, 6);
  FREE_DVECTOR(lmean, 1, 6);
  FREE_DVECTOR(beta, 0,nparm-1);

  if (fclose (fp_in) != 0 || fclose (fp_out) != 0 )
    ERRORPRT ("Error in closing opened files.");
  return 0;

} /* end of main */

/* Local Variables: */
/* compile-command: "make -f ../rMakefile OBJS=multista_plot.o EXEBASE=10multista" */
/* End: */


