#include<stdio.h>
#include<stdlib.h>
#include<string.h>
/******************************************************************
**
*  IMPORTANT NOTE:  If significant changes are being added to this
*					plotting program, please change the version
*					number in the corresponding model code
*
******************************************************************/

/*********************************************************************
*
*  Hill_plot.c - an ANSI C program that reads the output file 
*                   from Poly.C and then produces a script file for
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

#include "benchmark.h"
#define  float double
#define  NR_END 1  
#define INIVALUE 2012818004

char  fout[FLENGTH];
FILE  *fp_in, *fp_out;

int main(int argc, char *argv[])
{
	double *DVECTOR(int n1, int n2);
	void FREE_DVECTOR (double *v, int n1, int n2);
	void ERRORPRT (char error_text[]);

	char  name1[CNLENGTH], name2[CNLENGTH], name3[CNLENGTH], name4[CNLENGTH], 
		name5[CNLENGTH], name6[CNLENGTH], name7[CNLENGTH], 
		name8[CNLENGTH], name9[CNLENGTH], name10[CNLENGTH], name11[CNLENGTH], name12[CNLENGTH],
		name13[CNLENGTH], name14[CNLENGTH], name15[CNLENGTH], name16[CNLENGTH], name17[CNLENGTH], name18[CNLENGTH];
	char name4a[CNLENGTH];  /* contains the string: BMRType */
	char name4b[CNLENGTH];  /* contains the string: BMRF */
	char long_path_name[FLENGTH];

	int   flag1, flag2, flag3, smooth, nparm, nobs;
	int   i, icount;
	int rtype;    /* bmr type */
	double bmrf;

	double  var1, var2, tmpbeta, conf;
	double  xmin, xmax, ymin, ymax, rsl, bmd, bmdl;
	double  xleft, xright, ybottom, ytop, xrange, yrange, xlabel, ylabel;
	double  tdose, tmean, tmlb, tmup;
	double  bmd11, bmd12, bmd21, bmd22, bmd31, bmd32;
	double  bmdL11, bmdL12, bmdL21, bmdL22;

	double  *parms, *dose, *mean, *mlb, *mup, *ldose, *lmean;
	char    *parm_name[]={"b", "v", "n", "k"};
	int file_flag;
	char *pcRiskType = NULL;	/* Risk type for plot title */

	ldose=DVECTOR(1, 6);
	lmean=DVECTOR(1, 6);
	parms=DVECTOR(1, 4);

	file_flag = 0;
	if (argc > 2)
	{
		path_name2(argc, argv, long_path_name);
		argv[1] = long_path_name;
	} /* end if */

	fp_in=fopen(argv[1], "r");
	strcpy(fout, argv[1]);

	for(i = 0; i < FLENGTH; i++)
	{
		if(fout[i] == '.')
		{
			if(fout[i+1] == '0' && fout[i+2] == '0' && fout[i+3] == '2')
			{
				fout[i+1] = 'p';
				fout[i+2] = 'l';
				fout[i+3] = 't';
				fout[i+4] = '\0';
				break;
			}
		}
		if(fout[i] == '\0')
			break;
	}


	fp_out=fopen(fout, "w");

	if (fp_in==NULL || fp_out==NULL)
	{
		printf("Error in opening input and/or output files. \n");
		printf("...now exiting the system... \n");
		fprintf(fp_out, "Error in opening input and/or output files. \n");
		fprintf(fp_out, "...now exiting the system... \n");
		exit(1);
	}


	for (i=1; i<=4; i++)  parms[i]=0.0;

	/*** Read the input data ***/

	fscanf(fp_in, "%s %d", name1, &flag1);
	fscanf(fp_in, "%s %d", name2, &nobs);
	fscanf(fp_in, "%s %d", name3, &nparm);
	fscanf(fp_in, "%s %lg", name4, &conf);
	fscanf(fp_in, "%s %d", name4a, &rtype);
	switch (rtype) {
	case 0:
	  pcRiskType = "Abs. Dev.";
	  break;
	case 1:
	  pcRiskType = "Std. Dev.";
	  break;
	case 2:
	  pcRiskType = "Rel. Dev.";
	  break;
	case 3:
	  pcRiskType = "Point"; /* point estimate */
	  break;
	case 4:
	  pcRiskType = "Extra"; /* extra risk */
	  break;

	  /* Future risk types would go here before the default case */
	default:
	  pcRiskType = "Unknown Risk Type";
	  break;
	} /* end switch */

	fscanf(fp_in, "%s %lg", name4b, &bmrf);
	fscanf(fp_in, "%s %lg", name5, &var1);
	fscanf(fp_in, "%s %lg", name6, &var2);

	if(flag1 == -1)
	{
		fclose (fp_in);
		fclose (fp_out);
		exit(0);
	}

	for (i=1; i<=4; i++)
	{
		fscanf(fp_in, "%s %lg", name7, &tmpbeta);
		parms[i]=tmpbeta;
	}

	fscanf(fp_in, "%s", name8);

	/*** Add to handle error in input data file 08/30/03 Qun He ***/

	tdose = INIVALUE;
	tmean = INIVALUE;
	tmlb = INIVALUE;
	tmup = INIVALUE;

	if (nobs == INIVALUE)
	{
		ERRORPRT ("\n*****\nNo data from data file.\nProcess aborted.\n*****\n");
		if (fclose (fp_in) != 0 || fclose (fp_out) != 0 )
			ERRORPRT ("\n*****\nError in closing opened files.\nProcess aborted.\n*****\n");
	}
	/*** End of adding codes ***/

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

		/*** Add to handle data errors in input data file. 8/30/03 Qun He ***/
		if (dose[i] == INIVALUE || 
			mean[i] == INIVALUE || 
			mlb[i] == INIVALUE || 
			mup[i] == INIVALUE)
		{
			ERRORPRT ("\n*****\nNo data from data file.\nProcess aborted.\n*****\n");
			if (fclose (fp_in) != 0 || fclose (fp_out) != 0 )
				ERRORPRT ("\n*****\nError in closing opened files.\nProcess aborted.\n*****\n");
		}
		/*** End of adding codes ***/
	}

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
	}

	fscanf(fp_in, "%s %d", name10, &flag2);

	if (flag2 == 1)
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
	}

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
			/*	if ((tdose < 0) && (tmean < 0))
			flag3 = 0; */
		}

		/*   fscanf(fp_in, "%s %d", name19, &result); */

	}

	/*** Start to produce the output file ***/


	//if (argv[2] == NULL)
	fprintf(fp_out, 
#ifndef TERM_X11
		"set terminal windows dashed \n"
#else
		"set terminal x11\n"
#endif
		);
	/*else
	{
	fprintf(fp_out, "set terminal postscript monochrome \n");
	fprintf(fp_out, "set output \"Hill.ps\" \n");
	}
	*/
	fprintf(fp_out, "reset\n");
	fprintf(fp_out, "set time \"\%%H:\%%M \%%m/\%%d \%%Y\" \n");
	fprintf(fp_out, "set bar 3 \n");
	fprintf(fp_out, "set style line 6 linecolor rgb \"black\"\n");
	fprintf(fp_out, "set style line 16 linecolor rgb \"forest-green\" pt 12\n");
	fprintf(fp_out, "set key top left \n");
	fprintf(fp_out, "set xlabel 'dose' \n");
	fprintf(fp_out, "set ylabel 'Mean Response' \n");
	fprintf(fp_out, "set mxtics 10 \n");
	fprintf(fp_out, "set mytics 10 \n");

	if (flag1 == 1)  
	{
	  if (rtype != 3) {
	    fprintf(fp_out, "set title 'Hill Model, with BMR of %lg %s for the BMD and %lg Lower Confidence Limit for the BMDL' \n",
		    bmrf, pcRiskType, conf);
	  } else {
	    /* point estimate */
	    fprintf(fp_out, "set title 'Hill Model, with a Point Estimate BMR of %lg for the BMD and %lg Lower Confidence Limit for the BMDL' \n",
		    bmrf, conf);
	  }
		fprintf(fp_out, "rl = %lg \n", rsl);
		fprintf(fp_out, "bmd = %lg \n", bmd);              
	}

	if (flag2 == 1)
		fprintf(fp_out, "bmdl = %lg \n", bmdl); 
	else
		fprintf(fp_out, "set title 'Hill Model' \n");

	for (i=0; i<=3; i++)
		fprintf(fp_out, "%s = %lg \n", parm_name[i], parms[i+1]);

	fprintf(fp_out, "f(x)= b+(x>0?(v/((k/x)**n + 1)) : 0)\n");

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
	fprintf (fp_out, "set offsets %g, %g, 0, 0\n", xrange/20., xrange/20.);

	if (flag1 == 1) 
	{
		fprintf(fp_out, "set label 'BMD' at bmd, %lg left \n", ylabel);
	}

	if (flag2 == 1)
	{
		fprintf(fp_out, "set label 'BMDL' at bmdl, %lg right \n", ylabel);
		fprintf(fp_out, "set label '   ' at %g, rl left \n", xlabel);
	}

	if (flag1 == 1 && flag2 ==1 && flag3 ==1)
	{
		fprintf(fp_out, "plot [x=%f:%f] [%f:%f] f(x) title 'Hill',\\", 
			xmin, xmax, ybottom, ytop);
		fprintf(fp_out, "\n     '-' using 1:2:3:4 notitle with errorbars ls 16,\\");
		fprintf(fp_out, "\n     '-' using 1:2 notitle with lines ls 6,\\");
		fprintf(fp_out, "\n     '-' using 1:2 notitle with lines ls 6,\\");
		if (smooth == 1)
			fprintf(fp_out, "\n     '-' using 1:2 smooth csplines title 'BMD Lower Bound',\\");
		else
			fprintf(fp_out, "\n     '-' using 1:2 smooth unique title 'BMD Lower Bound',\\");
		fprintf(fp_out, "\n     '-' using 1:2 notitle with points");
	}

	else 
		if (flag1 == 1 && flag2 == 1)
		{
			fprintf(fp_out, "plot [x=%f:%f] [%f:%f] f(x) title 'Hill',\\", 
				xmin, xmax, ybottom, ytop);
			fprintf(fp_out, "\n     '-' using 1:2:3:4 notitle with errorbars ls 16,\\");
			fprintf(fp_out, "\n     '-' using 1:2 notitle with lines ls 6,\\");
			fprintf(fp_out, "\n     '-' using 1:2 notitle with lines ls 6 ");
		}

		else
			if (flag1 == 1)
			{
				fprintf(fp_out, "plot [x=%f:%f] [%f:%f] f(x) title 'Hill',\\", 
					xmin, xmax, ybottom, ytop);
				fprintf(fp_out, "\n     '-' using 1:2:3:4 notitle with errorbars ls 16,\\");
				fprintf(fp_out, "\n     '-' using 1:2 notitle with lines ls 6 ");
			}

			else
			{
				fprintf(fp_out, "plot [x=%f:%f] [%f:%f] f(x) title 'Hill',\\", 
					xmin, xmax, ybottom, ytop);
				fprintf(fp_out, "\n     '-' using 1:2:3:4 notitle with errorbars ls 16");
			}


			for (i=1; i<=nobs; i++)
				fprintf(fp_out, "\n %f %f %f %f", dose[i], mean[i], mlb[i], mup[i]);
			fprintf(fp_out, "\n e");

			bmd11=xleft;   // to adjust the positions of the BMD and BMDL lines
			bmd32=bmdL12=ybottom;

			if (flag1 == 1)
			{
				fprintf(fp_out, "\n %f %f", bmd11, bmd12);
				fprintf(fp_out, "\n %f %f", bmd21, bmd22);
				fprintf(fp_out, "\n %f %f", bmd31, bmd32);
				fprintf(fp_out, "\n e");
			}

			if (flag1 == 1 && flag2 == 1)
			{
				fprintf(fp_out, "\n %f %f", bmdL11, bmdL12);
				fprintf(fp_out, "\n %f %f", bmdL21, bmdL22);
				fprintf(fp_out, "\n e");
			}

			if (flag1 == 1 && flag2 == 1 && flag3 == 1)
			{
				for (i=1; i<=icount; i++)
					fprintf(fp_out, "\n %f %f", ldose[i], lmean[i]);
				fprintf(fp_out, "\n e");
				for (i=2; i<=icount; i++)
					fprintf(fp_out, "\n %f %f", ldose[i], lmean[i]);
				fprintf(fp_out, "\n e");
			}


			FREE_DVECTOR(dose, 1, nobs);
			FREE_DVECTOR(mean, 1, nobs);
			FREE_DVECTOR(mlb, 1, nobs);
			FREE_DVECTOR(mup, 1, nobs);
			FREE_DVECTOR(ldose, 1, 6);
			FREE_DVECTOR(lmean, 1, 6);
			FREE_DVECTOR(parms, 1, 4);

			// CLOSE_FILES ();
			if (fclose (fp_in) != 0 || fclose (fp_out) != 0 )
				ERRORPRT ("Error in closing opened files.");


}

/******************************************
**  ERRORPRT--used to handle error message.
*******************************************/  
void ERRORPRT (char error_text[])
{
	printf("\n%s\n",error_text);
	fprintf(fp_out," \n%s",error_text);

	if (fclose (fp_in) != 0 || fclose (fp_out) != 0 )
	{
		printf("\n Error in closing opened files. \n");
		fprintf(fp_out,"\n Error in closing opened files. \n");
	}
	exit (1);                                
}

/******************************************************************************
**  DVECTOR--allocate a double vector with subscript range [n1,n2].
*******************************************************************************/ 
double *DVECTOR(int n1, int n2)
{
	double *v;

	v=(double *) malloc((size_t) ((n2-n1+1+NR_END)*sizeof(double))); 
	if (!v) ERRORPRT ("Memory allocation failed in DVECTOR()");
	return v-n1+NR_END; 
}

/******************************************************************************
**  FREE_DVECTOR--free a double vector with subscript range [n1,n2].
*******************************************************************************/ 
void FREE_DVECTOR (double *v, int n1, int n2)
{
	free ((FREE_ARG) (v+n1-NR_END));
} 
