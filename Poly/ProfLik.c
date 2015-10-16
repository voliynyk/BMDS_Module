/***************************************************************** */
/* ProfLik.c -- a prototype function to compute points for profile */
/*              likelihoods for the */
/*              BMD.  This will work with all our models that are  */
/*              fit with dmngb, and thus have a function BMDL_func */
/*                                                                 */
/***************************************************************** */

#include <benchmark.h>
#include <allo_memo.h>
#include <string.h>
#include <ERRORPRT.h>
#include <specialfun.h>

typedef enum {Base1, Parm1, Parm2, Parm3, Parm4, Parm5, Parm6,
              Parm7, Parm8, Parm9, Parm10, Parm11, Parm12, Parm13,
              Parm14,Parm15, Parm16, Parm17, Parm18} Model_Parms;
extern char *Parm_name[];
/*  extern int ErrorFlag; */
/*  extern double tD; */

extern void getprofile_(double parms[],double *BMD, double *Max_LL,
			long int *optite,long int *tmpcnt,double *BMR2); 

/* void MAX_lk(int nparm, double parm[], double gtol, int *iter, double *fret); */
int fixedParm(int);

/* will need to change value of replace, so be sure to save the old value: */
/* this code will be used at the end of the run, but the state of the */
/* system should not be changed by running it, just to be safe. */
/* extern int replace; */
extern int *Spec;
/* int save_replace; */
#define NSTEPSD 20
#define NSTEPSU 20
#define BUFFSIZE 1000
#define MaxIter 40


void ProfLik(int nparm, double Parms[], double Max_LL, double BMD, double BMDL,
	     char *fname, double alpha, double BMR2, int *ProfFlag, int *PEFlag,
	     double *PE_bmdl, int *EmptyFlag)
{
  FILE *fp_out3;
  FILE *fp_log2;
  double *pvals, *psav, *LLvals, *BMDvals, delta, gtol, xlk, LLdeltalim,
    LLdelta1s, LLdelta2s, deltaLL, *BMDdown, *BMDup, *LLdown, *LLup;
  int i, nsd, nsu, mlpos, *LLok, NSTEPS,itercount;
  char outfile[186], *dot, fname2[186], logfile[186],*dot2;
  double tD;
  long int tmpcnt,tmpcnt2,ErrorFlag = -9999;
  double ab_bmdl,be_bmdl,ab_LL,be_LL,target_tmp;

  /* silent produces a log out file of what the process values     */
  /* are through the process, if silent = false log is created     */
  /*                             silent = true  log is not created */
  int silent = true;
  /*silent = false;*/

  /* open the log file if silent = 0 */
  if (silent == false){
    strcpy(logfile,fname);
    dot2 = strchr(logfile, (int) '.');
    (*dot2) = (char) 0;
    strcpy(fname2,logfile);
    strcat(logfile,"-BMD-profile.log");
    fp_log2 = fopen(logfile, "w");
    if (fp_log2 == (FILE *) NULL) {
      ERRORPRT("Unable to open log for Profile.");
      return;
    }
  }

  /* save_replace = replace; */
  /* replace = Yes; */
  NSTEPS = NSTEPSD + NSTEPSU + 1;
  mlpos = NSTEPSD + 1;
  if (silent == false){
    fprintf(fp_log2,"NSTEPS = %d\n",NSTEPS);
    fprintf(fp_log2,"BMD = %f\n",BMD);
    fprintf(fp_log2,"BMDL = %f\n",BMDL);
    fprintf(fp_log2,"Max_LL = %f\n", Max_LL);
    fprintf(fp_log2,"BMR = %f\n",BMR2);
    fflush(fp_log2);
  }

  /* allocate memory */
  BMDvals = DVECTOR(1,NSTEPS);
  LLvals = DVECTOR(1, NSTEPS);
  pvals = DVECTOR(1, nparm);
  psav = DVECTOR(1, nparm);
  LLok = IVECTOR(1, NSTEPS);
  BMDdown = DVECTOR(0,BUFFSIZE);
  BMDup = DVECTOR(0,BUFFSIZE);
  LLdown = DVECTOR(0,BUFFSIZE);
  LLup = DVECTOR(0,BUFFSIZE);

  gtol = 1e-7; /* I think this doesn't get used */
  /* Set up the likelihood decrement to target */

  if (alpha > 0.5) alpha = 1.0 - alpha;
  LLdelta2s = QCHISQ(1.0 - alpha, 1)/2.0;      /* 2 sided test */
  LLdelta1s = QCHISQ(1.0 - 2.0*alpha, 1)/2.0;  /* 1 sided test */
  LLdeltalim = QCHISQ(0.99,1)/2.0;             /* limit on the LL */

  if (silent == false){
    fprintf(fp_log2,"\nLLdelta2s = %f\nLLdelta1s = %f\nLLdeltalim = %f\n\n",
	    LLdelta2s,LLdelta1s,LLdeltalim);
    fflush(fp_log2);
  }
  /* Open the output file.  It has two columns, BMD and LL */
  strcpy(outfile,fname);
  dot = strchr(outfile, (int) '.');
  (*dot) = (char) 0;
  strcpy(fname2,outfile);
  strcat(outfile,"-BMD-profile.plt");
  fp_out3 = fopen(outfile, "w");
  if (fp_out3 == (FILE *) NULL) {
    ERRORPRT("Unable to open output for Profile.");
    return;
  }
  /* Write The preamble stuff */
  fprintf(fp_out3,"set terminal windows\nreset\nset xlabel 'BMD'\n");
  fprintf(fp_out3,"set size square\n");
  fprintf(fp_out3,"set ylabel 'Log Likelihood'\n");
  fprintf(fp_out3,"set mxtics 10\n");
  fprintf(fp_out3,"set mytics 10\n");
  fprintf(fp_out3,"set x2tics (\"BMDL\" %g, \"BMD\" %g)\n",BMDL, BMD);
  fprintf(fp_out3,"set y2tics (%g, %g)\n",
	  Max_LL - LLdelta1s, Max_LL - LLdelta2s);
  fprintf(fp_out3,"set grid x2tics y2tics\n");
  fprintf(fp_out3,"set title '%s:Profile Likelihood for the BMD'\n",fname2);

  /* set starting values for BMD, BMDL, and Log Likelihood  */
  if (BMDL < 0) BMDL = 0.0;
  BMDdown[0] = BMDup[0] = BMD;

  /* Compute LL values for decreasing BMD */
  deltaLL = LLdelta2s / NSTEPSD;
  for (i = 1; i <= nparm; i++) pvals[i] = Parms[i];
  nsd = 0;

  if (silent == false){
    fprintf(fp_log2,"\n**********************************************\n\n");
    fprintf(fp_log2,"Find LL using getprofile at the BMD\n\n");
    fprintf(fp_log2,"Before getprofile run\n");
    for (i=1; i <= nparm; i++)
      fprintf(fp_log2,"parameter %d = %g\n",i,pvals[i]);
    fflush(fp_log2);
  }
  tmpcnt = 0;
  tD = BMD;
  getprofile_(pvals,&tD,&xlk,&ErrorFlag,&tmpcnt,&BMR2);
  if (silent == false) {
    fprintf(fp_log2,"After getprofile run\n");
    for (i=1; i <= nparm; i++)
      fprintf(fp_log2,"parameter %d = %g\n",i,pvals[i]);
    fprintf(fp_log2,"xlk = %13.7g     (This is the max LL found by getprofile at the BMD)\n",
	    xlk);
    fprintf(fp_log2,"ErrorFlag = %2ld\n",ErrorFlag);
    fprintf(fp_log2,"tD = %13.4g      (The BMD for the first run of getprofile)\n",tD);
  }
  Max_LL = xlk;
  LLdown[0] = LLup[0] = Max_LL;
/*    tmpbmd = BMD; */
  tmpcnt = 0;
  tmpcnt2 = 1;
/*    pvals[1]= 1.5;  */
/*    getprofile_(pvals,&tmpbmd,&xlk,&ErrorFlag,&tmpcnt,&BMR2); */
/*    fprintf(fp_out3,"%f  %f  %ld  %f\n",tmpbmd,xlk,ErrorFlag,deltaLL); */
/*    for (i=1; i <= nparm; i++) */
/*      fprintf(fp_out3,"parameter %d = %f\n",i,pvals[i]); */
/*    fprintf(fp_out3,"BMD = %f, BMDL = %f, xlk = %f, BMR = %f",BMD,BMDL,xlk,BMR2); */
  if (silent == false){
    fprintf(fp_log2,"*****************************************************************\n");
    fprintf(fp_log2,"****                    Downward Loops                       ****\n");
    fprintf(fp_log2,"*****************************************************************\n");
  }

  do
    {
      delta = (BMD - BMDL)/NSTEPSD;
      for (i = 1; i <= nparm; i++) psav[i] = pvals[i];
      itercount = 0;
      if (silent == false){
	fprintf(fp_log2,"\n**********************************************\n\n");
	fprintf(fp_log2,"Main Loop %ld\n\n",tmpcnt2);
	for (i=1; i <= nparm; i++)
	  fprintf(fp_log2,"parameter %d = %g\n",i,pvals[i]);
	fflush(fp_log2);
      }
      do
	{
	  tD = BMDdown[nsd] - delta;
	  tmpcnt++;
	  if (silent == false){
	    fprintf(fp_log2,"\n\nInner Loop: %d\n",itercount+1);
	    fprintf(fp_log2,"BMD value going into getprofile: %f\n",tD);
	    fflush(fp_log2);
	  }
	  getprofile_(pvals,&tD,&xlk,&ErrorFlag,&tmpcnt,&BMR2);
	  /*	  MAX_lk(nparm, pvals, gtol, &junk, &xlk);   */
	  if (silent == false){
	    fprintf(fp_log2,"The ErrorFlag value from getprofile is: %ld\n",ErrorFlag);
	    fprintf(fp_log2,"The Likehood value produced is: %f\n",xlk);
	    fflush(fp_log2);
	  }
	  if (LLdown[nsd] - xlk > deltaLL ||
	      (ErrorFlag !=0 && ErrorFlag != -1))
	    {
	      delta *= 0.9;
	      for (i = 1; i <= nparm; i++) pvals[i] = psav[i];
	      itercount++;
	    }
	      
	} while ((LLdown[nsd] - xlk > deltaLL ||
		  (ErrorFlag !=0 && ErrorFlag != -1))
		 && itercount < MaxIter);
      if (itercount < MaxIter && tD > 0)
	{
	  nsd++;
	  BMDdown[nsd] = tD;
	  LLdown[nsd] = xlk;
	}
      tmpcnt2++;
      
    } while (nsd < BUFFSIZE-1 && LLdown[nsd] > (Max_LL - LLdeltalim) &&
	     itercount < MaxIter && tD > 0);
  if (silent == false){
    fprintf(fp_log2,"\n\n\n*****************************************************************\n");
    fprintf(fp_log2,"****                    Upward Loops                         ****\n");
    fprintf(fp_log2,"*****************************************************************\n");
  }

  tmpcnt2 = 1;
  /* Compute LL values for increasing BMD */
  deltaLL = LLdelta2s / NSTEPSU;
  for (i = 1; i <= nparm; i++) pvals[i] = Parms[i];
  nsu = 0;
  do
    {
      delta = (BMD - BMDL)/NSTEPSU;
      for (i = 1; i <= nparm; i++) psav[i] = pvals[i];
      itercount = 0;
      if (silent == false){
	fprintf(fp_log2,"\n**********************************************\n\n");
	fprintf(fp_log2,"Main Loop %ld\n\n",tmpcnt2);
	for (i=1; i <= nparm; i++)
	  fprintf(fp_log2,"parameter %d = %g\n",i,pvals[i]);
	fflush(fp_log2);
      }
      do
	{
	  tD = BMDup[nsu] + delta;
	  tmpcnt++;
	  if (silent == false){
	    fprintf(fp_log2,"\n\nInner Loop: %d\n",itercount+1);
	    fprintf(fp_log2,"BMD value going into getprofile: %f\n",tD);
	    fflush(fp_log2);
	  }
	  getprofile_(pvals,&tD,&xlk,&ErrorFlag,&tmpcnt,&BMR2); 
	  /*  MAX_lk(nparm, pvals, gtol, &junk, &xlk); */
	  if (silent == false){
	    fprintf(fp_log2,"The ErrorFlag value from getprofile is: %ld\n",ErrorFlag);
	    fprintf(fp_log2,"The Likehood value produced is: %f\n",xlk);
	    fflush(fp_log2);
	  }
	  if (LLup[nsu] - xlk > deltaLL ||
	      (ErrorFlag !=0 && ErrorFlag != -1))
	    {
	      delta *= 0.9;
	      for (i = 1; i <= nparm; i++) pvals[i] = psav[i];
	      itercount++;
	    }
	      
	} while ((LLup[nsu] - xlk > deltaLL ||
		  (ErrorFlag !=0 && ErrorFlag != -1))
		 && itercount < MaxIter);
      if (itercount < MaxIter)
	{
	  nsu++;
	  BMDup[nsu] = tD;
	  LLup[nsu] = xlk;
	}
      
      tmpcnt2++;
    } while (nsu < BUFFSIZE-500 && LLup[nsu] > (Max_LL - LLdeltalim) &&
	     itercount < MaxIter);
  fprintf(fp_out3,"set xrange [%g:%g]\n",BMDdown[nsd],BMDup[nsu]);
  fprintf(fp_out3,"set x2range [%g:%g]\n",BMDdown[nsd],BMDup[nsu]);
  fprintf(fp_out3,"plot '-' using 1:2 smooth unique notitle\n");
  for (i = nsd; i >= 0; i--)
    {
      fprintf(fp_out3,"%g %g\n",BMDdown[i], LLdown[i]);
    }
  for (i = 1; i <= nsu; i++)
    {
      fprintf(fp_out3,"%g %g\n",BMDup[i], LLup[i]);
    }
  fprintf(fp_out3,"e\n");
  fprintf(fp_out3,"reset\n");
  fclose(fp_out3);

  /*This part gives an interpolation estimate to the bmdl */
  /*to check the estimate from the main program */
  target_tmp = Max_LL-LLdelta1s;
  ab_bmdl = -1;
  be_bmdl = -1;
  ab_LL = -1;
  be_LL = -1;

  for (i = 0; i<nsd; i++)
    {
      if (LLdown[i] > target_tmp)
	{
	  ab_bmdl = BMDdown[i];
	  be_bmdl = BMDdown[i+1];
	  ab_LL = LLdown[i];
	  be_LL = LLdown[i+1];
	}
    }
  if (ab_bmdl>0 && be_bmdl>0 && (ab_LL-be_LL)>0 && (ab_bmdl-be_bmdl)>0 && nsd>10)
    {
      *PE_bmdl = ab_bmdl - ((ab_LL-target_tmp)/(ab_LL-be_LL))*(ab_bmdl-be_bmdl);
      *PEFlag = 0;
    }
  else
    *PEFlag = 1;
  if(BMDdown[nsd]>BMDL) 
    *ProfFlag = 1;
  else
    *ProfFlag =0;
  if(nsd < 3)
    *EmptyFlag = 1;
  else
    *EmptyFlag = 0;
  if (silent == false){
    fprintf(fp_log2,"\n\n*********************************************************\n");
    fprintf(fp_log2,"Profile Estimate Flag = %d\n",*PEFlag);
    fprintf(fp_log2,"Profile Estimate of BMDL = %g\n",*PE_bmdl);
    fprintf(fp_log2,"Profile Plot Flag = %d\n",*ProfFlag);
    fprintf(fp_log2,"Empty Plot Flag = %d\n",*EmptyFlag);
    fflush(fp_log2);
  }
  fclose(fp_log2);

  FREE_DVECTOR(BMDvals, 1, NSTEPS);
  FREE_DVECTOR(LLvals, 1, NSTEPS);
  FREE_DVECTOR(pvals, 1, nparm);
  FREE_DVECTOR(psav, 1, nparm);
  FREE_IVECTOR(LLok, 1, NSTEPS);
  FREE_DVECTOR(BMDdown,0,BUFFSIZE);
  FREE_DVECTOR(LLdown,0,BUFFSIZE);
  FREE_DVECTOR(BMDup,0,BUFFSIZE);
  FREE_DVECTOR(LLup,0,BUFFSIZE);
  /*  replace = save_replace; */
}
