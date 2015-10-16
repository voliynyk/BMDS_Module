/*********************************************************************
 *
 *  This routine converts *.002 output from the Exponential model into
 *  a GNU Plot input file.
 *
 *  NOTE: The use of character arrays should be replaced with mallocs to avoid
 *        possible cross contamination of data elements if multiple program
 *        instances call this DLL simultaneously.
 *  Actually, this might not be an issue since these aren't global or static variables.
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
#ifndef WIN32
# include <sys/errno.h>
#endif
#include "benchmark.h"
#include "bmds_plot.h"

#define  float double

char gacProgram[FLENGTH];

//static FILE *gfplog;

void get_input(FILE *fp_in, const eModel_t eModel, char acModel[],
	      bmds_i_parm_t aziParms[], bmds_d_parm_t azdParms[])
{
  int iP_i=0;
  int iP_d=0;
  int iModel; /* Index for model name and function */
  char acTemp[BMDS_NAME_SIZE];
  int i; /* loop variable */

  memset(acTemp, '\0', sizeof(acTemp));
  iModel = eModel;
  //  fprintf(gfplog, "pcModels[%d]= %s", eModel, pcModels[eModel]);
  if (eModel == evExponential) {
    int iSubModel;
    fscanf(fp_in, "%s %d", acTemp, &iSubModel);
    iModel += (iSubModel - iF_EXPO_OFFSET);
  } /* evExponential */
  strcpy(acModel, pcModels[iModel]);
  strcpy(aziParms[iP_i].Name, "Model Index");
  //fprintf(gfplog, "acModel= %s", acModel);
  aziParms[iP_i].iValue = iModel;
  iP_i++;

  /* Read: BMD_flag, Nobs, nparm */
  for (i = 0; i < 3; i++) {
    fscanf(fp_in, "%s %d", aziParms[iP_i].Name, &aziParms[iP_i].iValue);
    iP_i++;
  } /* for */
  if (eModel == evExponential) {
    /* Read: sign, lognorm */
    iP_i++; /* "Save" a spot for BMR Type (don't need to since array) */
    fscanf(fp_in, "%s %d", aziParms[iP_i].Name, &aziParms[iP_i].iValue);
    iP_i++;
    fscanf(fp_in, "%s %d", aziParms[iP_i].Name, &aziParms[iP_i].iValue);
    iP_i++;
  } /* evExponential */

  fscanf(fp_in, "%s %lg", azdParms[iP_d].Name, &azdParms[iP_d].dValue); /* Confidence level */
  iP_d++;
  fscanf(fp_in, "%s %d", aziParms[iBMRTYPE].Name, &aziParms[iBMRTYPE].iValue); /* BMR Type */
  /* Ideally we should validate that BMR Type is valid */

  /* Read: BMTF, lnalpha/alpha, rho,...
   * a/control, b/slope, c/power
   */
  for (i = iBMRF; i < iP4; i++) {
    fscanf(fp_in, "%s %lg", azdParms[i].Name, &azdParms[i].dValue);
    iP_d++;
  } /* for */
  if (eModel == evExponential) {
    fscanf(fp_in, "%s %lg", azdParms[iP_d].Name, &azdParms[iP_d].dValue);
    iP_d++;
  } /* evExponential */
  fscanf(fp_in, "%s", acTemp); /* "Data" */

}

void plot_continuous(const eModel_t eModel, int argc, char *argv[])
{
  double *DVECTOR(int n1, int n2);
  void FREE_DVECTOR (double *v, int n1, int n2);

  FILE *fp_in, *fp_out;
  char acModel[BMDS_NAME_SIZE];
  char acFunc[BMDS_NAME_SIZE]; // Holds modified function for expo 2 & 3
  char fout[FLENGTH];
  char acInFile[FLENGTH];
  char acGraphFile[FLENGTH];
  /* For simplicity, hard code maximum size for array of parmameter
   * names and values. Later, we should allocate the space more
   * intelligently.
   */
  bmds_i_parm_t aziParms[BMDS_MAX_I_PARMS];
  bmds_d_parm_t azdParms[BMDS_MAX_D_PARMS];

  char name9[CNLENGTH], name10[CNLENGTH], name11[CNLENGTH], name12[CNLENGTH],
    name13[CNLENGTH], name14[CNLENGTH], name15[CNLENGTH], name16[CNLENGTH], name17[CNLENGTH], name18[CNLENGTH];

  int i, icount;
  int flag1, flag2, flag3, smooth, nparm, nobs;
  int rtype, iSign, iLognorm;
  /* int iSign; */
  double  conf, bmrf;
  double  xmin, xmax, ymin, ymax, rsl, bmd, bmdl;
  double  xleft, ybottom, ytop, xrange, yrange, xlabel, ylabel;
  double  tdose, tmean, tmlb, tmup;
  double  bmd11, bmd12, bmd21, bmd22, bmd31, bmd32;
  double  bmdL11, bmdL12, bmdL21, bmdL22;

  double  *dose, *mean, *mlb, *mup, *ldose, *lmean;

  /* These are character POINTERS. There is no memory allocated!! */
  const char *pcRiskType;	/* Risk type for plot title */
  const char *pcFunc;		/* Plot function */
  char *pcPtr;
  char cParm;

  memset(acModel, '\0', sizeof(acModel));
  memset(acFunc, '\0', sizeof(acFunc));
  memset(&aziParms, '\0', sizeof(aziParms));
  memset(&azdParms, '\0', sizeof(azdParms));

  strncpy(gacProgram, argv[0], FLENGTH-1);
#ifdef DO_LOG
  /* Generate log file name based on program name. We temporarily */
  /* borrow the output file name buffer to gen the log file name. */
  /* NOTE... This uses the full path name in argv[0]. We should probably
   * just use the base file name component.
   */
  strcpy(fout, gacProgram);
  pcPtr = strrchr(fout, '.');
  strcpy(++pcPtr, "log");
  gfplog = fopen(fout, "w+");
  for (i=1; i < argc; i++) fprintf(gfplog, "argv[%d] = %s\n", i, argv[i]);
#endif

  /* Consider adding a check for too many input arguments. */
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <input file>\n", argv[0]);
    exit(1);
  }

  strncpy(acInFile, argv[1], FLENGTH-1);
  strcpy(fout, acInFile);
  /* Generate output file name based on input name. */
  pcPtr = strrchr(fout, '.');
  strcpy(++pcPtr, "plt");

  /* image file name that gnuplot will create */
  strcpy(acGraphFile, acInFile);
  pcPtr = strrchr(acGraphFile, '.');
  strcpy(++pcPtr, "emf");
#ifdef DO_LOG
  fprintf(gfplog, "Input file name = %s\n", acInFile);
  fprintf(gfplog, "Output file name = %s\n", fout);
  fprintf(gfplog, "Image file name = %s\n", acGraphFile);
#endif

  fp_in=fopen(acInFile, "r");
  if (fp_in == NULL) {
    fprintf(stderr, "%s could not open input file %s. Errno=%d\n",
	    argv[0], acInFile, errno);
    exit(2);
  }
  fp_out=fopen(fout, "w");
  if (fp_out == NULL) {
    fprintf(stderr, "%s could not open output file %s. Errno=%d\n",
	    argv[0], fout, errno);
    exit(2);
  }

  /*** Get first batch of input parameters ***/
  (void) get_input(fp_in, eModel, acModel, aziParms, azdParms);

  /* Store values in convenience variables */
  flag1 = aziParms[iBMDFLAG].iValue;
  nobs = aziParms[iNOBS].iValue;
  nparm = aziParms[iNPARMS].iValue;
  rtype = aziParms[iBMRTYPE].iValue;
  iSign = aziParms[iSIGN].iValue;
  iLognorm = aziParms[iLOGNORM].iValue;
  conf = azdParms[iCONF].dValue;
  bmrf = azdParms[iBMRF].dValue;
  pcRiskType = pcRiskTypes[rtype];
  pcFunc = pcModelFuncs[aziParms[iFUNC].iValue];
  if (aziParms[iFUNC].iValue == iF_EXPO2
      || aziParms[iFUNC].iValue == iF_EXPO3) {
    sprintf(acFunc, pcFunc, iSign);
    pcFunc = acFunc;
  } /* iF_EXPO2 || iF_EXPO3 */
  /*** Read the input data ***/

  dose=DVECTOR(1, nobs);
  mean=DVECTOR(1, nobs);
  mlb=DVECTOR(1, nobs);
  mup=DVECTOR(1, nobs);
  ldose=DVECTOR(1, nobs);
  lmean=DVECTOR(1, nobs);

  /** Find the max and min of x and y values **/
  ymin=1.0e+66;
  ymax=-1.0e+66;
  for (i=1; i<=nobs; i++) {
    fscanf(fp_in, "%lg %lg %lg %lg", &tdose, &tmean, &tmlb, &tmup);
    if (iLognorm) {
      /* Convert from log-scale to regular. */
      /* Having condition inside loop is OK since only a few iterations. */
      tmean=exp(tmean);
      tmlb=exp(tmlb);
      tmup=exp(tmup);
    }	/* end if iLognorm */
    dose[i]=tdose;
    mean[i]=tmean;
    mlb[i]=tmlb;
    mup[i]=tmup;
    if (tmean < ymin)  ymin=tmean;
    if (tmlb < ymin)   ymin=tmlb;
    if (tmean > ymax)  ymax=tmean;
    if (tmup > ymax)   ymax=tmup;
  } /* end for i = 1 to nobs */

  fscanf(fp_in, "%s", name9); /* "Max_Min_dose" */
  fscanf(fp_in, "%lg %lg", &xmax, &xmin);

  if (flag1 == 1) { /* if BMD calculation was selected */
    fscanf(fp_in, "%s %lg", name11, &rsl); /* RSL */
    fscanf(fp_in, "%s %lg", name12, &bmd); /* BMD */
    fscanf(fp_in, "%s", name14);	       /* "BMD_line" */
    fscanf(fp_in, "%lg %lg", &bmd11, &bmd12);
    fscanf(fp_in, "%lg %lg", &bmd21, &bmd22);
    fscanf(fp_in, "%lg %lg", &bmd31, &bmd32);
  } /* end if */
  fscanf(fp_in, "%s %d", name10, &flag2); /* BMDL_comput_ind */
  if (flag2 == 1) { /* if bmdl was computed */
    fscanf(fp_in, "%s %lg", name13, &bmdl); /* BMDL */
    if (rsl > ymax) ymax=rsl;
    if (rsl < ymin) ymin=rsl;
    if (bmd > xmax) xmax=bmd;
    if (bmdl < xmin) xmin=bmdl;
    fscanf(fp_in, "%s", name15); /* "BMDL_line" */
    fscanf(fp_in, "%lg %lg", &bmdL11, &bmdL12);
    fscanf(fp_in, "%lg %lg", &bmdL21, &bmdL22);	  
  } /* if (flag2 == 1) */

  fscanf(fp_in, "%s %d", name16, &flag3); /* BMDL_Curve_flag */
  fscanf(fp_in, "%s %d", name17, &smooth); /* smooth_opt */

  if (flag3 == 1) { /* if all bmdls were calculated */
    fscanf(fp_in, "%s", name18); /* "BMDL_curve" */
    icount=0;
    for (i=1; i<=6; i++) {
      fscanf(fp_in, "%lg %lg", &tdose, &tmean);
      if (tdose >= 0) {
	icount++;
	ldose[icount]=tdose;
	lmean[icount]=tmean;
	if (tdose < xmin) xmin=tdose;
	if (tdose > xmax) xmax=tdose;
	if (tmean > ymax) ymax=tmean;
	if (tmean < ymin) ymin=tmean;
      } /* end if tdose >= 0 */
    } /* end for */
  } /* end if (flag3 == 1) */



  /*** Start to produce the output file ***/

  fprintf(fp_out, "set terminal windows dashed\n");
  fprintf(fp_out, "reset\n");
  //fprintf(fp_out, MYTERM);
  //fprintf(fp_out, "set output \"%s\"\n", acGraphFile);
  fprintf(fp_out, "set time \"\%%H:\%%M \%%m/\%%d \%%Y\"\n");
  fprintf(fp_out, "set bar 3\n");
  fprintf(fp_out, acSetBMDLineStyle);
  fprintf(fp_out, acSetErrorBarStyle);
  fprintf(fp_out, "set key top left\n");
  fprintf(fp_out, "set xlabel 'dose'\n");
  fprintf(fp_out, "set ylabel 'Mean Response'\n");
  fprintf(fp_out, "set mxtics 10\n");
  fprintf(fp_out, "set mytics 10\n");

  if (flag1 == 1) {
    if (rtype != 3) {
      fprintf(fp_out, "set title '%s Model, with BMR of %lg %s for the BMD and %lg Lower Confidence Limit for the BMDL'\n",
	      acModel, bmrf, pcRiskType, conf);
    } else {
      /* Point Estimate */
      fprintf(fp_out, "set title '%s Model, with Point Estimate BMR of %lg for the BMD and %lg Lower Confidence Limit for the BMDL'\n",
	      acModel, bmrf, conf);
    } /* if (rtype != 3) - else */
    fprintf(fp_out, "rl = %lg \n", rsl);
    fprintf(fp_out, "bmd = %lg \n", bmd);
  } /* end if flag1 == 1 */
  if (flag2 == 1) fprintf(fp_out, "bmdl = %lg\n", bmdl);
  else fprintf(fp_out, "set title '%s'\n", acModel);

  /* Output lines for variance and model parameters */
  for (cParm='a',i=0; i < nparm; i++) {
    fprintf(fp_out, "%c = %lg\n", cParm, azdParms[iALPHA + i].dValue);
    cParm++;
  }
  fprintf(fp_out, "f(x)= %s\n", pcFunc);

  if (ymin == ymax) {
    ymin -= 0.5;
    ymax += 0.5;
  }
  xrange=xmax-xmin;
  yrange=ymax-ymin;
  xleft=xmin-xrange/10.0;
  /*xright=xmax+xrange/10.0;*/
  ybottom=ymin-yrange/10.0;
  ytop=ymax+yrange/10.0;
  xlabel=xmin-xrange/12.0;
  ylabel=ymin-yrange/15.0;
  fprintf(fp_out, "set offset %g, %g, 0, 0\n", xrange/20.0, xrange/20.0);
  if (flag1 == 1) {
    fprintf(fp_out, "set label 'BMD' at bmd, %lg left\n", ylabel);
  } /* end if */

  if (flag2 == 1) {
    fprintf(fp_out, "set label 'BMDL' at bmdl, %lg right\n", ylabel);
    fprintf(fp_out, "set label '   ' at %g, rl left\n", xlabel);
  } /* end if */

  if (flag1 == 1 && flag2 == 1 && flag3 == 1) {
    fprintf(fp_out, "plot [x=%f:%f] [%f:%f] f(x) title '%s',\\", 
	    xmin, xmax, ybottom, ytop, acModel);
    fprintf(fp_out, "\n\t'-' using 1:2:3:4 notitle %s,\\", acErrorBars);
    fprintf(fp_out, "\n\t'-' using 1:2 notitle %s,\\", acLines02);
    fprintf(fp_out, "\n\t'-' using 1:2 notitle %s, \\", acLines02);
    if (smooth == 1) fprintf(fp_out, "\n\t'-' using 1:2 smooth csplines title 'BMD Lower Bound',\\");
    else fprintf(fp_out, "\n\t'-' using 1:2 smooth unique title 'BMD Lower Bound',\\");
    fprintf(fp_out, "\n\t'-' using 1:2 notitle with points");
  } /* end if */
  else if (flag1 == 1 && flag2 ==1) {
    fprintf(fp_out, "plot [x=%f:%f] [%f:%f] f(x) title '%s',\\", 
	    xmin, xmax, ybottom, ytop, acModel);
    fprintf(fp_out, "\n\t'-' using 1:2:3:4 notitle %s,\\", acErrorBars);
    fprintf(fp_out, "\n\t'-' using 1:2 notitle %s,\\", acLines02);
    fprintf(fp_out, "\n\t'-' using 1:2 notitle %s", acLines02);
  } /* end if */
  else if (flag1 == 1) { 
    fprintf(fp_out, "plot [x=%f:%f] [%f:%f] f(x) title '%s',\\", 
	    xmin, xmax, ybottom, ytop, acModel);
    fprintf(fp_out, "\n\t'-' using 1:2:3:4 notitle %s,\\", acErrorBars);
    fprintf(fp_out, "\n\t'-' using 1:2 notitle %s", acLines02);
  } /* end if */
  else {
    fprintf(fp_out, "plot [x=%f:%f] [%f:%f] f(x) title '%s',\\", 
	    xmin, xmax, ybottom, ytop, acModel);
    fprintf(fp_out, "\n\t'-' using 1:2:3:4 notitle %s", acErrorBars);
  } /* end else */

  for (i=1; i<=nobs; i++) fprintf(fp_out, "\n %f %f %f %f", dose[i], mean[i], mlb[i], mup[i]);
  fprintf(fp_out, "\ne");

  bmd11=xleft;   // to adjust the positions of the BMD and BMDL lines
  bmd32=bmdL12=ybottom;

  if (flag1 == 1) {
    fprintf(fp_out, "\n%f %f", bmd11, bmd12);
    fprintf(fp_out, "\n%f %f", bmd21, bmd22);
    fprintf(fp_out, "\n%f %f", bmd31, bmd32);
    fprintf(fp_out, "\ne");
  } /* end if */

  if (flag1 == 1 && flag2 == 1) {	
    fprintf(fp_out, "\n%f %f", bmdL11, bmdL12);
    fprintf(fp_out, "\n%f %f", bmdL21, bmdL22);
    fprintf(fp_out, "\ne");
  } /* end if */

  if (flag1 == 1 && flag2 == 1 && flag3 == 1) {
    for (i=1; i<=icount; i++) fprintf(fp_out, "\n%f %f", ldose[i], lmean[i]);
    fprintf(fp_out, "\ne");
    for (i=2; i<=icount; i++) fprintf(fp_out, "\n%f %f", ldose[i], lmean[i]);
    fprintf(fp_out, "\ne");
  } /* end if */
  fprintf(fp_out, "\n");

  FREE_DVECTOR(dose, 1, nobs);
  FREE_DVECTOR(mean, 1, nobs);
  FREE_DVECTOR(mlb, 1, nobs);
  FREE_DVECTOR(mup, 1, nobs);
  FREE_DVECTOR(ldose, 1, 6);
  FREE_DVECTOR(lmean, 1, 6);

  if (fclose (fp_in)) {
    fprintf(stderr, "%s: Error closing input file %s. Errno=%d\n",
	    argv[0], acInFile, errno);
    exit(2);
  }
  if (fclose(fp_out)) {
    fprintf(stderr, "%s: Error closing output file %s. Errno=%d\n",
	    argv[0], fout, errno);
    exit(2);
  }
#ifdef DO_LOG
  fclose(gfplog);
#endif
}

/******************************************************************************
 **  DVECTOR--allocate a double vector with subscript range [n1,n2].
 *******************************************************************************/ 
#define  NR_END 1
double *DVECTOR(int n1, int n2)
{
  size_t iBytes = (size_t) ((n2-n1+1+NR_END)*sizeof(double));
  double *v;

  v=(double *) malloc(iBytes); 
  if (!v) {
    fprintf(stderr, "%s: Failed to allocate %d bytes of memory. Errno=%d\n",
	    gacProgram, iBytes, errno);
    exit(2);
  }
  return v-n1+NR_END;
}

/******************************************************************************
 **  FREE_DVECTOR--free a double vector with subscript range [n1,n2].
 *******************************************************************************/ 
void FREE_DVECTOR (double *v, int n1, int n2)
{
  free ((FREE_ARG) (v+n1-NR_END));
} 
