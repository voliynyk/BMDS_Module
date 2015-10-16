/******************************************************************
**
*  IMPORTANT NOTE:  If significant changes are being added to this
*		    plotting program, please change the version
*		    number in the corresponding model code
*
******************************************************************/

/*********************************************************************
*
*  Quanln_plot.c - an ANSI C program that reads the output file 
*                   from Quanln.C and then produces a script file for
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

char  fout[FLENGTH];
FILE  *fp_in, *fp_out;

void CLOSE_FILES (void)
{
  if (fclose (fp_in) != 0 || fclose (fp_out) != 0 )
    ERRORPRT ("Error in closing opened files.");
}




int main(int argc, char *argv[])
{
  double *DVECTOR(int n1, int n2);
  void FREE_DVECTOR (double *v, int n1, int n2);

  char  name1[CNLENGTH], name2[CNLENGTH], name3[CNLENGTH], name4[CNLENGTH], name5[CNLENGTH], name6[CNLENGTH], 
    name8[CNLENGTH], name9[CNLENGTH], name10[CNLENGTH], name11[CNLENGTH], name12[CNLENGTH],
    name13[CNLENGTH], name14[CNLENGTH], name15[CNLENGTH], name16[CNLENGTH], name17[CNLENGTH], name18[CNLENGTH],
    name19[CNLENGTH];

  int   flag1, flag2, flag3, smooth, nparm, nobs, result;
  int   i, icount;

  double  background, slope, conf;
  double  xmin, xmax, ymin, ymax, rsl, bmd, bmdl;
  double  xleft, xright, ybottom, ytop, xrange, yrange, xlabel, ylabel;
  double  tdose, tmean, tmlb, tmup;
  double  bmd11, bmd12, bmd21, bmd22, bmd31, bmd32;
  double  bmdL11, bmdL12, bmdL21, bmdL22;

  double  *dose, *mean, *mlb, *mup, *ldose, *lmean;
  char long_path_name[FLENGTH];
  char *extpoint;

  ldose=DVECTOR(1, 6);
  lmean=DVECTOR(1, 6);
  /* Allow spaces in path name */
  if (argc > 2)
    {
      path_name2(argc, argv, long_path_name);
      argv[1] = long_path_name;
    } /* end if */

  fp_in=fopen(argv[1], "r");

  strcpy(fout, argv[1]);

  extpoint = strrchr(fout,(int) '.');
  if (extpoint != (char *) NULL)
    {
      *extpoint = '\0';
    }
  strcat(fout, ".plt");

  strcat(fout, "Quanln.plt");
 
  fp_out=fopen(fout, "w");

  if (fp_in==NULL || fp_out==NULL)
    {
      printf("Error in opening input and/or output files. \n");
      printf("...now exiting the system... \n");
      fprintf(fp_out, "Error in opening input and/or output files. \n");
      fprintf(fp_out, "...now exiting the system... \n");
      exit(1);
    }

  /*** Read the input data ***/

  fscanf(fp_in, "%s %d", name1, &flag1);
  fscanf(fp_in, "%s %d", name2, &nobs);
  fscanf(fp_in, "%s %d", name3, &nparm);
  fscanf(fp_in, "%s %lg", name4, &conf);
  fscanf(fp_in, "%s %lg", name5, &background);
  fscanf(fp_in, "%s %lg", name6, &slope);

  fscanf(fp_in, "%s", name8);

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
    }

  fscanf(fp_in, "%s", name9);
  fscanf(fp_in, "%lg %lg", &xmax, &xmin);

  if (flag1 == 1)
    {
      fscanf(fp_in, "%s %d", name10, &flag2);
      if (flag2 == 1)
	{
	  fscanf(fp_in, "%s %lg", name11, &rsl);
	  fscanf(fp_in, "%s %lg", name12, &bmd);
	  fscanf(fp_in, "%s %lg", name13, &bmdl);
	  if (rsl > ymax) ymax=rsl;
	  if (rsl < ymin) ymin=rsl;
	  if (bmd > xmax) xmax=bmd;
	  if (bmdl < xmin) xmin=bmdl;

	  fscanf(fp_in, "%s", name14);
	  fscanf(fp_in, "%lg %lg", &bmd11, &bmd12);
	  fscanf(fp_in, "%lg %lg", &bmd21, &bmd22);
	  fscanf(fp_in, "%lg %lg", &bmd31, &bmd32);

	  fscanf(fp_in, "%s", name15);
	  fscanf(fp_in, "%lg %lg", &bmdL11, &bmdL12);
	  fscanf(fp_in, "%lg %lg", &bmdL21, &bmdL22);
	  
	  fscanf(fp_in, "%s %d", name16, &flag3);
	  if (flag3 == 1)
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
		    }
		}
	
	      fscanf(fp_in, "%s %d", name19, &result);
	    }
	}
    }

  /*** Start to produce the output file ***/

  if (argv[2] == NULL)
    fprintf(fp_out, 
#ifndef TERM_X11
	    "set terminal windows \n"
#else
	    "set terminal x11\n"
#endif
	    );
  else
    {
      fprintf(fp_out, "set terminal postscript monochrome \n");
      fprintf(fp_out, "set output \"Quanln.ps\" \n");
    }
  fprintf(fp_out, "reset\n");
  fprintf(fp_out, "set time \"\%%H:\%%M \%%m/\%%d \%%Y\" \n");
  fprintf(fp_out, "set bar 3 \n");
  fprintf(fp_out, "set key top left \n");
  fprintf(fp_out, "set xlabel 'dose' \n");
  fprintf(fp_out, "set ylabel 'Fraction Affected' \n");
  fprintf(fp_out, "set mxtics 10 \n");
  fprintf(fp_out, "set mytics 10 \n");
  if (flag1 == 1 && flag2 == 1)
    {
      fprintf(fp_out, "set title 'Quantal Linear Model, with %lg%% MR and %lg lower confidence limit for BMD (BMDL)' \n", rsl*100., conf);
      fprintf(fp_out, "rl = %g \n", rsl);
      fprintf(fp_out, "bmd = %g \n", bmd);
      fprintf(fp_out, "bmdl = %g \n", bmdl);
    }
  else
    fprintf(fp_out, "set title 'Quantal Linear Model' \n");

  fprintf(fp_out, "back = %g \n", background);
  fprintf(fp_out, "slope = %g \n", slope);

  fprintf(fp_out, "f(x)= back + (x>0?( (1- back)*(1-exp(-slope*x)) ) : 0) \n");

  if (ymin == ymax)
    {
      ymin -= 0.5;
      ymax += 0.5;
    }
  xrange=xmax-xmin;
  yrange=ymax-ymin;
  xleft=xmin-xrange/10.0;
  xright=xmax+xrange/10.0;
  ybottom=ymin-yrange/10.0;
  ytop=ymax+yrange/10.0;
  xlabel=xmin-xrange/12.0;
  ylabel=ymin-yrange/15.0;
  fprintf(fp_out, "set offset %g, %g, 0, 0\n", xrange/20.0, xrange/20.0);
  if (flag1 == 1 && flag2 == 1)
    {
      fprintf(fp_out, "set label 'BMDL' at bmdl, %g right \n", ylabel);
      fprintf(fp_out, "set label 'BMD' at bmd, %g left \n", ylabel);
      fprintf(fp_out, "set label '   ' at %g, rl left \n", xlabel);
    }

  fprintf(fp_out, "plot [x=%f:%f] [%f:%f] f(x) title 'Quantal Linear',\\", 
	  xmin, xmax, ybottom, ytop);
  fprintf(fp_out, "\n     '-' using 1:2:3:4 notitle with errorbars");

  if (flag1 == 1 && flag2 == 1)
    {
      fprintf(fp_out, ",\\");
      if (flag3 == 1)
	{
	  if (smooth == 1)
	    fprintf(fp_out,
		    "\n     '-' using 1:2 smooth csplines title 'BMD Lower Bound',\\");
	  else
	    fprintf(fp_out,
		    "\n     '-' using 1:2 smooth unique title 'BMD Lower Bound',\\");
	  fprintf(fp_out, "\n     '-' using 1:2 notitle with points,\\");
	}
      fprintf(fp_out, "\n     '-' using 1:2 notitle with lines 2,\\");
      fprintf(fp_out, "\n     '-' using 1:2 notitle with lines 2 ");
    }
  else
    {
      fprintf(fp_out, "\n");
    }


  for (i=1; i<=nobs; i++)
    fprintf(fp_out, "\n %f %f %f %f", dose[i], mean[i], mlb[i], mup[i]);
  fprintf(fp_out, "\n e");

  bmd11=xleft;   // to adjust the positions of the BMD and BMDL lines
  bmd32=bmdL12=ybottom;

  if (flag1 == 1 && flag2 == 1 && flag3 == 1)
    {
      for (i=1; i<=icount; i++)
	fprintf(fp_out, "\n %f %f", ldose[i], lmean[i]);
      fprintf(fp_out, "\n e");
      for (i=2; i<=icount; i++)
	fprintf(fp_out, "\n %f %f", ldose[i], lmean[i]);
      fprintf(fp_out, "\n e");
    }

  if (flag1 == 1 && flag2 == 1)
    {
      fprintf(fp_out, "\n %f %f", bmd11, bmd12);
      fprintf(fp_out, "\n %f %f", bmd21, bmd22);
      fprintf(fp_out, "\n %f %f", bmd31, bmd32);
      fprintf(fp_out, "\n e");
      fprintf(fp_out, "\n %f %f", bmdL11, bmdL12);
      fprintf(fp_out, "\n %f %f", bmdL21, bmdL22);
      fprintf(fp_out, "\n e");
    }

  if (argv[2] != NULL)
    fprintf(fp_out, "\n EXIT");

  FREE_DVECTOR(dose, 1, nobs);
  FREE_DVECTOR(mean, 1, nobs);
  FREE_DVECTOR(mlb, 1, nobs);
  FREE_DVECTOR(mup, 1, nobs);
  FREE_DVECTOR(ldose, 1, 6);
  FREE_DVECTOR(lmean, 1, 6);

  // CLOSE_FILES ();
  if (fclose (fp_in) != 0 || fclose (fp_out) != 0 )
    ERRORPRT ("Error in closing opened files.");
  return 0;

}

/* Local Variables: */
/* compile-command: "make -f ../rMakefile OBJS=quanln_plot.o EXEBASE=10quanln" */
/* End: */
