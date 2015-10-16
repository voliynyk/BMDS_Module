/***************************************************************** */
/* ProfLik.c -- a prototype function to compute points for profile */
/*              likelihoods for the */
/*              BMD.  This will work with all our models that are  */
/*              fit with dmngb, and thus have a function BMDL_func */
/*                                                                 */
/***************************************************************** */

#define ASSIST
#include <benchmark.h>
#include <allo_memo.h>
#include <string.h>
#include <ERRORPRT.h>
#include <specialfun.h>

typedef enum {Base1, Parm1, Parm2, Parm3, Parm4, Parm5, Parm6,
              Parm7, Parm8, Parm9, Parm10, Parm11, Parm12, Parm13,
              Parm14,Parm15, Parm16, Parm17, Parm18} Model_Parms;
extern char *Parm_name[];
extern int ErrorFlag;
extern double tD;

void MAX_lk(int nparm, double parm[], double gtol, int *iter, double *fret);
int fixedParm(int);

/* will need to change value of replace, so be sure to save the old value: */
/* this code will be used at the end of the run, but the state of the */
/* system should not be changed by running it, just to be safe. */
extern int replace;
extern int *Spec;
int save_replace;
#define NSTEPSD 20
#define NSTEPSU 20
#define BUFFSIZE 1000
#define MaxIter 40

void ProfLik(int nparm, double Parms[], double Max_LL,
	     double BMD, double BMDL, char *fname, double alpha)
{
  FILE *fp_out3;
  double *pvals, *psav, *LLvals, *BMDvals, delta, gtol, xlk, LLdeltalim,
    LLdelta1s, LLdelta2s, deltaLL, *BMDdown, *BMDup, *LLdown, *LLup;
  int i, nsd, nsu, mlpos, junk, *LLok, NSTEPS, itercount;
  char outfile[186], *dot, fname2[186];

  save_replace = replace;
  replace = Yes;
  NSTEPS = NSTEPSD + NSTEPSU + 1;
  mlpos = NSTEPSD + 1;
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
  LLdelta2s = QCHISQ(1.0 - alpha, 1)/2.0;
  LLdelta1s = QCHISQ(1.0 - 2.0*alpha, 1)/2.0;
  LLdeltalim = QCHISQ(0.99,1)/2.0;
  /* Open the output file.  It has two columns, BMD and LL */
  strcpy(outfile,fname);
  dot = strchr(outfile, (int) '.');
  (*dot) = (char) 0;
  strcpy(fname2,outfile);
  strcat(outfile,"-BMD-profile.plt");
  fp_out3 = fopen(outfile, "w");
  if (fp_out3 == (FILE *) NULL) {
    Warning("Unable to open output for Profile.");
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

  if (BMDL < 0) BMDL = 0.0;
  BMDdown[0] = BMDup[0] = BMD;
  LLdown[0] = LLup[0] = Max_LL;
  /* Compute LL values for decreasing BMD */
  deltaLL = LLdelta2s / NSTEPSD;
  for (i = 1; i <= nparm; i++) pvals[i] = Parms[i];
  nsd = 0;
  do
    {
      delta = (BMD - BMDL)/NSTEPSD;
      for (i = 1; i <= nparm; i++) psav[i] = pvals[i];
      itercount = 0;
      do
	{
	  tD = BMDdown[nsd] - delta;
	  MAX_lk(nparm, pvals, gtol, &junk, &xlk);
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
      
      
    } while (nsd < BUFFSIZE-1 && LLdown[nsd] > (Max_LL - LLdeltalim) &&
	     itercount < MaxIter && tD > 0);
  /* Compute LL values for increasing BMD */
  deltaLL = LLdelta2s / NSTEPSU;
  for (i = 1; i <= nparm; i++) pvals[i] = Parms[i];
  nsu = 0;
  do
    {
      delta = (BMD - BMDL)/NSTEPSU;
      for (i = 1; i <= nparm; i++) psav[i] = pvals[i];
      itercount = 0;
      do
	{
	  tD = BMDup[nsu] + delta;
	  MAX_lk(nparm, pvals, gtol, &junk, &xlk);
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
      
      
    } while (nsu < BUFFSIZE-1 && LLup[nsu] > (Max_LL - LLdeltalim) &&
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
  
  
  
  FREE_DVECTOR(BMDvals, 1, NSTEPS);
  FREE_DVECTOR(LLvals, 1, NSTEPS);
  FREE_DVECTOR(pvals, 1, nparm);
  FREE_DVECTOR(psav, 1, nparm);
  FREE_IVECTOR(LLok, 1, NSTEPS);
  FREE_DVECTOR(BMDdown,0,BUFFSIZE);
  FREE_DVECTOR(LLdown,0,BUFFSIZE);
  FREE_DVECTOR(BMDup,0,BUFFSIZE);
  FREE_DVECTOR(LLup,0,BUFFSIZE);
  replace = save_replace;
}
