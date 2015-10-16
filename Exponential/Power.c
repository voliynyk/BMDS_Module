/****************************************************************
**
* IMPORTANT NOTE:  The following variable is the version number for
*                  the current model.  THIS MUST BE CHANGED as
*	    	   important changes are made to the models.
*
*****************************************************************/

char Version_no[]="Power Model. (Version: 2.14;  Date: 02/20/2007)";
/*
char Version_no[]="Power Model. (Version: 2.12;  Date: 12/08/2006)";
char Version_no[] = "Power Model. Revision: 2.3 Date: 2005/07/24 20:57:36";
*/

/****************************************************************
*
*
* Power.C - a ANSI C program for the Continuous Power model fitting
*
* Date: Oct, 2000
*
********************************************************************
* Modification Log:
*
* Version Number: 2.3
* Modified By: Q. He
* Modified Date: 09/08/2003
* Reason:
*
* Version Number: 2.4
* Modified By: Micheal Ferree
* Modified Date: 08/13/2005
* Reason:
*
* Version Number: 2.5
* Modified By: R. Woodrow Setzer
* Modified Date: 10/17/2005
* Reason: Free all allocated memory, fix bugs due to using uninitialized data.
* Remaining Problems: remove Stud_t_func, use the t routines in Assist
*
* Version Number: 2.6
* Modified By: R. Woodrow Setzer
* Modified Date: 12/06/2005
* Reason: Always allocate anasum with 5 elements.
*
* Version Number: 2.7
* Modified By: R. Woodrow Setzer
* Modified Date: 03/21/2006
* Reason: moved calculation of likelihoods for models A1, A2, and R to
*         calculate_continuous_liks() (fixes errors in likelihoods for
*         A1 and R)
*
*
* Version Number: 2.8
* Modified By: R. Woodrow Setzer
* Modified Date: 03/24/2006
* Reason: Calculate likelihood for A3 using fixed values for alpha and rho,
*         if supplied
*
* Version Number: 2.9
* Modified By: R. Woodrow Setzer
* Modified Date: 05/04/2006
* Reason: Confidence intervals were calculated incorrectly when using
*         a constant
*         variance model (Ntot was not set)
*
* Version Number: 2.10
* Modified By: R. Woodrow Setzer
* Modified Date: 07/11/2006
* Reason: This time, really use user supplied fixed values for alpha and
*         rho in calculating likelihood for A3.
*
* Version Number: 2.11
* Modified By: R. Woodrow Setzer
* Modified Date: 08/04/2006
* Reason: i) Change variance model parameterization
*         ii) Clean up constraint specification in donlp2
*         iii) In BMDL calculation, change equality likelihood constraint
*              to inequality constraint.
*
* Version Number: 2.12
* Modified By: R. Woodrow Setzer
* Modified Date: 12/08/2006
* Reason: fixed bug in initial value for constant variance case:
*         was finding log(alpha), and reporting using it as alpha
*
* Version Number: 2.13
* Modified By: Geoffrey
* Date: 1/12/2007
* Reason: Incremented version number.
*		 Added last parameter "1" (print SE) in OP_ParmsE().
*
* Version Number: 2.14
* Modified By: Woodrow Setzer
* Date: 2/20/2007
* Reason: Incremented version number to reflect changed compilation options.
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


extern void getmle_(long int *ndoses, double doses[], double means[],
		    long int nanimals[], double svar[], long int *nparm,
		    double parms[], long int fixed[], double fixedval[],
		    long int *restrict,long int *adverse, double parms2[],double *ll,
		    long int *optite, long int *nresm, long int bind[],
		    long int *model_type, long int *flag);

extern void getcl_(long int *which, long int *ndoses, double doses[], double means[],
		   long int nanimals[], double svar[], long int *nparm, double *bmr,
		   double *bmd, double *target, double parms[],
		   long int fixed[], double fixedval[], long int *risktype,
		   long int *restrict, double *bmdl, double parms2[],
		   long int *optite, long int *nresm, long int bind[], long int *adv,
		   long int *model_type, long int *flag);

extern void getmlea3_(long int *ndoses, double doses[], double means[],
		      long int nanimals[], double svar[], long int *nparm,
		      double parms[], long int fixed[], double fixedval[],
		      long int *restrict, double parms2[],double *ll,
		      long int *optite, long int *nresm, long int bind[]);

extern void loadcommbloc_(double *maxYm, double *maxYd);

void PowMeans(int nobs, double p[], double Doses[], double means[]);
void Get_BMRS(double *p, double Dose, double bmdl1, double *BMRVals, int sign, int bmr_type);
void OUTPUT_BENCHMD2(int pdcol, double BMD);
void GetMoreParms(double *p, int size);
double power(double a, double b);
void GetOtherParms(double *p, int size);
void GetMLEParms(double *p, int size);
void GetNewParms(double *p, int size);
int Model_DF(int []);


#define EPS 3.0e-8
#define EPSS 3.0e-7
#define FPMIN 1.0e-30
#define MAXIT 100
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
char    fin[128];   /*input temp file*/
char    fout[128];  /*output temp file*/
char    fout2[128];
char	plotfilename[128];  /* file to pass to GnuPlot */
char    *Parm_name[]={"alpha", "rho", "control", "slope", "power"};
char    *alt_Parm1_name = "lalpha";
char    *anatxt[]={"Full model", "Fitted model", "Reduced model"};
char    fname2[186], logfile[186], *dot2;
/* moved to benchmark.h by GLN - 03/23/2005: 
   FILE    *fp_log;             log file if requested */
/*************************************************************
 *   bNo_Log = true -> no log file
 *   bNo_Log = false -> log file is made
 *************************************************************/
int     bNo_Log = true;     /*  switch for log file   */


/*** variables will not be changed except Spec, and xxi and yyi
     will be sorted together *******/

int    *Spec;     /*vector used to identify user input parm.*/
int    *IniSp;
int    *Ni;       /*number of animals in the i-th dose group*/
double *Ym;       /*mean response data array*/
double *ScYm;     /*scaled mean response data */
double maxYm;     /*max mean response */
double *Yd;       /*sample standard variance (s_i^2) of response data array*/
double *ScYd;     /*scaled sample variances of responses */
double maxYd;     /*max sample variance */
double *Xi;       /*independent variable (dose value) data array*/
double *xxi;      /*dose values when not divided to dose groups*/
double *yyi;      /*response values when not divided to dose groups*/
double *Ysum;     /*sum of responses within a group*/
double *Rlevel;
double *Bmdl;
double *IniP;
int    in_type, Nobs, ntotal, nparm, restrict, initial, appendix, smooth;
int    bmdlCurve, sign, BMDL_RUN;
double xmax, xmin, yymin, yymax, scale;
int    bmr_type;

/** changing variable **/
int    replace, brat;
double tD, BMD_lk, LR, ck, upb=18, BMR, react;
int    const_var;



/****************************************************************
 *  main--main function used to call Power mode fitting program.
  *****************************************************************/
int  main (int argc, char *argv[])
{
  void Pow_fit(int nparm, double p[], double gtol,int *iter, double *fret, int *bounded, int *is_conv);
  void Pow_vcv(int nparm, int Spec[], double p[], double **vcv);
  void Goodness(int nparm, int nparm_known, double Parms[], int type, AnaList anasum[]);
  int READ_OBSDATA4V(int Nobs,double Xi[],int Ni[],double Ym[],double Yd[]);
  int READ_OBSDATA2V(int ntotal, double xxi[], double yyi[]);
  void AThree_Fit(int nparm, double p[], double gtol, int *iter, double *fret);
  double BMDL_func(int nparm, double xlk, double Dose, double pBak[], double gtol);
  void Pow_BMD (int nparm, double p[], double gtol, int *iter, double xlk, 
		double Rlevel[], double Bmdl[],double *BMD);

  int     iter,i, j, jj, junk;         /*iteration variable*/
  int     group;                       /*temp counter eventually equals Nobs*/
  int     bmdose;                      /*flag for computing benchmark dose*/
  int     Nmiss;                       /*number of records with missing values*/
  int     nparm_known, power_known;    /*number of specified parameters */
  double  xlk, lkA3, lkA1, lkA2, lkR;  /*log likelihoods */
  double  BMD, lep, upep, Rel_Conv, Parm_Conv;
  double  *stdev;                      /*user input sample standard deviation*/
  double  *Parms, *LKParms;            /*parameter array*/
  VarList *varsum;                     /*info for variables--p. dep.,n. dep., indep.*/
  AnaList *anasum;                     /*information for ANONA analysis*/
  double  **vcv;                       /*variance and covariance matrix*/
  char    model_name[82], user_note[82];
  char    dose_name[20], no_name[20], mean_name[20], stdev_name[20], junkname[82];
  char    response_name[20];
  int     Ntot, var_type;    // const_var;
  double  bmr_root, sp;
  double  Nd;
  int     *bounded;
  int     conv_check;
  int     adj_vcv_rows;
  double  **vcv_adj;
  int     temp_sign;
  double  maxYm, maxYd;
  double  *mean, *std;                 /* fitted means and standard deviations */
  time_t  ltime;
  char long_path_name[300];

  /* Set time zone from TZ environment variable. If TZ is not set,
   * the operating system is queried to obtain the default value
   * for the variable.    */
  time( &ltime );


  /********************************************************************
   * {QH 2004/01/14 PR# }
   * Added to show version number if executed from command line with -v
   *********************************************************************/
  if(argc == 2)
    show_version(argv[1], Version_no);

  if(argc < 2)
    {
      fprintf(stderr, "ERROR:  Requires two arguments\nUsage:  %s <file.(d)>\n", argv[0]);
      fprintf (stderr, "   or:  %s -v for version number.\n", argv[0]);
      exit (1);
    } /* end if */


  if (argc > 2)
    {
      path_name2(argc, argv, long_path_name);
      argv[1] = long_path_name;
    }

  fp_in=fopen(argv[1], "r");
  if (fp_in==NULL)
    {
      printf("Error in opening input  file.\n");
      printf ("...now exiting to system...\n");
      exit (1);
    }

  /* open the log file if bNo_Log = 0 */
  if (bNo_Log == false){
    strcpy(logfile,argv[1]);
    dot2 = strchr(logfile, (int) '.');
    (*dot2) = (char) 0;
    strcpy(fname2,logfile);
    strcat(logfile,"-pow.log");
    fp_log = fopen(logfile, "w");

    if (fp_log == (FILE *) NULL)
      ERRORPRT("Unable to open log for Power.C.");
  }

  fscanf(fp_in, "%s", model_name);     /* scans the model name */
  fscanf(fp_in, "%[ ^\n]", user_note);
  fscanf(fp_in, "%[^\n]", user_note);
  fscanf(fp_in, "%s", junkname);
  fscanf(fp_in, "%s", junkname);
  fscanf(fp_in, "%d",&in_type);
  /* in_type=1 if input format is Xi, Ni, Ym, Yd. */
  /* in_type=0 if input format is Ntotal, Xi, Y_ij. */

  if (in_type==1)
    fscanf(fp_in, "%d", &Nobs);
  else
    fscanf(fp_in, "%d", &ntotal);

  fscanf(fp_in, "%d", &sign);
  if ((sign != 1) || (sign != -1))
    sign = 0;

  /*assign number of parameters*/
  nparm = 5;

  if (bNo_Log == false){
    fprintf(fp_log,"model_name: %s\n\n",model_name);
    fprintf(fp_log,"INPUT VALUES FOR DATA SET\n\n");
    fprintf(fp_log,"in_type = %d        (Input Format  in_type=1 if input format is Xi, Ni, Ym, Yd\n",in_type);
    fprintf(fp_log,"                                  in_type=0 if input format is Ntotal, Xi, Y_ij)\n");

    if (in_type==1)
      fprintf(fp_log,"Nobs = %d           (Number of Observations)\n",Nobs);
    else
      fprintf(fp_log,"ntotal = %d         (Total Number of Subjects)\n",ntotal);

    fprintf(fp_log,"sign = %2d          (Adverse Direction; 0=Automatic, 1=Up, -1=Down)\n",sign);
    fprintf(fp_log,"nparm = %d          (The Number of Parameters in the Model)\n",nparm);
    fflush(fp_log);
  }

  /*allocate memory for arrays*/
  Parms = DVECTOR(1, nparm);
  Spec = IVECTOR(1, nparm);
  IniP = DVECTOR(1, nparm);
  IniSp = IVECTOR(1, nparm);
  varsum = VLVECTOR(1, 3);
  Rlevel = DVECTOR(1, 5);
  Bmdl = DVECTOR(1, 5);

  fscanf(fp_in,"%d%lf%lf%d%d%d%d%d", &ITMAX, &Rel_Conv, &Parm_Conv, \
	 &bmdlCurve, &restrict, &bmdose, &appendix, &smooth);

  /* restrict=0 if no restriction, 1 if power restricted >= 1 */
  fscanf(fp_in,"%d%lf%d%lf",&bmr_type,&bmdparm.effect,&const_var,&bmdparm.level);

  bmdparm.risk = 1;
  junk = 0;    /* Used to see if an extension was added to output file name */

  if (bNo_Log == false){
    fprintf(fp_log,"ITMAX = %d                 (Maximum number of iterations)\n",ITMAX);
    fprintf(fp_log,"Rel_Conv = %g     (Rel Function Convergence, default=2.22045e-16)\n",Rel_Conv);
    fprintf(fp_log,"Parm_Conv = %g    (Parameter Convergence)\n",Parm_Conv);
    fprintf(fp_log,"bmdlCurve = %d               (BMDL Curve Calculation; 1=yes, 0=no)\n",bmdlCurve);
    fprintf(fp_log,"restrict = %2d               (Restriction on power coefficients)\n",restrict);
    fprintf(fp_log,"bmdose = %d                  (BMD Calculation; 1=yes, 0=no)\n",bmdose);
    fprintf(fp_log,"appendix = %d                (Append or Overwrite output file)\n",appendix);
    fprintf(fp_log,"smooth = %2d                 (Smooth Option)\n",smooth);
    fprintf(fp_log,"bmr_type = %d                (BMR Type)\n",bmr_type);
    fprintf(fp_log,"bmdparm.effect = %4g       (BMR Factor)\n",bmdparm.effect);
    fprintf(fp_log,"const_var = %d                (Constant Variance)\n",const_var);
    fprintf(fp_log,"bmdparm.level = %4g        (Confidence level)\n",bmdparm.level);
  }
  /* get filenames */
  Get_Names(argv[1], fout, fout2, plotfilename);

  if(appendix==Yes)
    fp_out=fopen(fout,"a");
  else
    fp_out=fopen(fout,"w");
#ifndef RBMDS
  fp_out2=fopen(fout2,"w");
#endif
  if (fp_out==NULL 
#ifndef RBMDS
      ||fp_out2==NULL
#endif
      )
    {
      printf("Error in opening  output files.\n");
      printf ("...now exiting to system...\n");

      fprintf(fp_out,"Error in opening output files.\n");
      fprintf (fp_out,"...Exited to system!\n");
      if (bNo_Log == false) {
	fflush(fp_log);
	fclose(fp_log);
      }
      CLOSE_FILES();
      exit (1);
    }

  /* Print model and file information on output page */
  Output_Header(Version_no, argv[1], plotfilename, ctime(&ltime), user_note);

  if (bmdose < 0 || bmdose > 1)
    ERRORPRT("Error in choosing benchmark dose computation.");

  /*obtain user input parameters*/
  READ_PARAMETERS(nparm,Parms);
  FILL_SPECVECTOR(nparm,Parms,Spec);

  nparm_known = COUNT_SPECVECTOR(nparm,Spec);
  power_known=0;

  for (i=3; i<=5; i++){
    if (Spec[i] == 1)
      power_known += 1;
  }

  brat=Yes;

  if (Spec[1]==1 && Parms[1]<EPS)
    brat=No;

  /*obtain user input initial parameters values      *******/
  fscanf(fp_in,"%d", &initial);
  READ_PARAMETERS(nparm,IniP);
  FILL_SPECVECTOR(nparm,IniP,IniSp);

  for(i = 1; i <= nparm; i++){
    if(Spec[i] == 1)
      IniP[i] = 1;
  }

  if (bNo_Log == false) {
    fprintf(fp_log,"\nThe specified or default parameters are:\n");
    fprintf(fp_log,"                                           Spec Value\n");
    fprintf(fp_log,"Parms[1] = %12.6g        (alpha)        %d\n",Parms[1],Spec[1]);
    fprintf(fp_log,"Parms[2] = %12.6g        (rho)          %d\n",Parms[2],Spec[2]);
    fprintf(fp_log,"Parms[3] = %12.6g        (control)      %d\n",Parms[3],Spec[3]);
    fprintf(fp_log,"Parms[4] = %12.6g        (slope)        %d\n",Parms[4],Spec[4]);
    fprintf(fp_log,"Parms[5] = %12.6g        (power)        %d\n",Parms[5],Spec[5]);
    fprintf(fp_log,"\ninitial = %d                (Number of initialized parameters)\n",initial);
    fprintf(fp_log,"The initial parameter values are:\n");
    fprintf(fp_log,"                                           IniSp Value\n");
    fprintf(fp_log,"IniP[1] = %12.6g         (alpha)         %d\n",IniP[1],IniSp[1]);
    fprintf(fp_log,"IniP[2] = %12.6g         (rho)           %d\n",IniP[2],IniSp[2]);
    fprintf(fp_log,"IniP[3] = %12.6g         (control)       %d\n",IniP[3],IniSp[3]);
    fprintf(fp_log,"IniP[4] = %12.6g         (slope)         %d\n",IniP[4],IniSp[4]);
    fprintf(fp_log,"IniP[5] = %12.6g         (power)         %d\n",IniP[5],IniSp[5]);
    fprintf(fp_log,"The Spec and IniSp vectors are just flags to tell if the user gave\n");
    fprintf(fp_log,"a value for the specific parameter.\n");
  }

  /*obtain observation data into Yp, Yn, Xi, Ls, Xg vectors*/
  if (in_type == 1)
    fscanf(fp_in,"%s%s%s%s", dose_name, no_name, mean_name, stdev_name);
  else
    fscanf(fp_in,"%s%s", dose_name, response_name);

  if(Spec[2] == 1){ /* Determine whether it a heterogeneous or homogeneous
		       model for likelihood tests */
    if(Parms[2] == 0) {
      anasum = ALVECTOR(1, 5);
      const_var = 1;
    }
    else{
      anasum = ALVECTOR(1, 5);
      const_var = 0;
    }
  }
  else {
    anasum = ALVECTOR(1, 5);
    const_var = 0;
  }
  if (const_var == 0) {
    Parm_name[0] = alt_Parm1_name;
  }

  if(bNo_Log == false) {
    fprintf(fp_log,"const_var = %d\n",const_var);
    fflush(fp_log);
  }

  if (in_type==1) {
    Ym = DVECTOR(1, Nobs);
    Yd = DVECTOR(1, Nobs);
    Xi = DVECTOR(1, Nobs);
    Ni = IVECTOR(1, Nobs);
    stdev = DVECTOR(1, Nobs);
    Nmiss = READ_OBSDATA4V(Nobs, Xi, Ni, Ym, stdev);
    for(i = 1; i <= Nobs; i++)
      Yd[i] = stdev[i]*stdev[i];

    /* extern variable Nobs has been changed. */
    Nobs -= Nmiss;
  }
  else {
    xxi = DVECTOR(1, ntotal);
    yyi = DVECTOR(1, ntotal);
    Ym = DVECTOR(1, ntotal);
    Yd = DVECTOR(1, ntotal);
    Xi = DVECTOR(1, ntotal);
    Ni = IVECTOR(1, ntotal);
    Ysum = DVECTOR(1, ntotal);
    Nmiss = READ_OBSDATA2V(ntotal, xxi, yyi);
    ntotal -= Nmiss;  /* extern variable Nobs has been changed. */

    for (i=1; i<=ntotal; i++) {
      Xi[i] = Ysum[i] = Ym[i] = Yd[i] = 0;
      Ni[i] = 0;
    }
    Sort_2_By_Dose(ntotal, xxi, yyi);  /*Sort xxi and make appropriate changes to yyi */
    group = 1;

    for (i=1; i<=ntotal; i++) {
      Xi[group]=xxi[i];
      Ni[group] += 1;
      Ysum[group] += yyi[i];
      if(i < ntotal) {
	if(xxi[i] != xxi[i + 1])
	  group +=1;
      }
    }
    Nobs=group;
    jj=1;

    for (i=1; i<=Nobs; i++) {
      Ym[i] = Ysum[i]/Ni[i];
      for (j=1; j<=Ni[i]; j++){
	Yd[i] += (yyi[jj]-Ym[i])*(yyi[jj]-Ym[i])/(Ni[i]-1);
	jj += 1;
      }
    }
  } // end else (in_type==1)

  ScYd = DVECTOR(1, Nobs);
  ScYm = DVECTOR(1, Nobs);
  maxYd = Yd[1];
  maxYm = Ym[1];

  for (i=2; i<=Nobs; i++){
    if(Yd[i] > maxYd)
      maxYd = Yd[i];
    if(Ym[i] > maxYm)
      maxYm = Ym[i];
  }

  for (i=1; i<=Nobs; i++){
    ScYd[i] = Yd[i]/maxYd;
    ScYm[i] = Ym[i]/maxYm;
  }

  if (bNo_Log == false) {
    fprintf(fp_log,"\nDOSE's          N's                MEAN's             VAR's\n");

    for (i=1; i<=Nobs; i++)
      fprintf(fp_log,"Xi[%d]=%4g      Ni[%d]=%4d         Ym[%d]=%4g         Yd[%d]=%4g\n", \
	      i,Xi[i],i,Ni[i],i,Ym[i],i,Yd[i]);

    fprintf(fp_log,"\nNmiss = %2d             (Number of missing values)\n",Nmiss);
    fprintf(fp_log,"Nobs = %2d              (Number of observations)\n",Nobs);
    fflush(fp_log);
  }

  if (Nobs < 3-power_known)
    ERRORPRT("Observation # < parameter # for Power model.");

  if (in_type == 1) {
    for (i=1; i<=Nobs; i++) {
      if ( (Xi[i] < 0) || (Ni[i] < 0) || (Yd[i] < 0) )
	ERRORPRT("Values of dose, group size, and sample variance should be positive...");
    }
  }
  else {
    for (i=1; i<=Nobs; i++) {
      if (Xi[i] < 0)
	ERRORPRT("Dose value should be positive ...");
    }
  }
  /*********** end of input data *****/

  maxYm = Ym[1];                   /***Previously done !!!  Q.H. 1/28/04 ***/
  maxYd = Yd[1];
  for (i = 2; i <= Nobs; i++){
    if (maxYm < Ym[i])
      maxYm = Ym[i];
    if (maxYd < Yd[i])
      maxYd = Yd[i];
  }

  loadcommbloc_(&maxYm,&maxYd);        /*This should be implemented at a later time!!!!*/
  /*MJF 10/11/00 */

  xmin = yymin = Max_double;
  xmax = yymax = 0.0;
  for (i=1;i<=Nobs;i++){
    if (Xi[i] < xmin) {
      xmin = Xi[i];
      bmr_root = sqrt(Yd[i]);
    }
    if (Xi[i] > xmax)
      xmax = Xi[i];
    if (Ym[i] > yymax)
      yymax = Ym[i];
    if (Ym[i] < yymin)
      yymin = Ym[i];
  }
  /*output title and summary of intput data  ****************************/
  OUTPUT_TEXT("\n   The form of the response function is: ");

  OUTPUT_TEXT("\n   Y[dose] = control + slope * dose^power");
  if (in_type == 1)
    fprintf(fp_out,"\n\n   Dependent variable = %s", mean_name);
  else
    fprintf(fp_out,"\n\n   Dependent variable = %s", response_name);

  fprintf(fp_out,"\n   Independent variable = %s", dose_name);

  for (i=1; i<=nparm; i++){
    if (Spec[i] == 1)
      fprintf(fp_out,"\n   %s is set to %g", Parm_name[i-1],Parms[i]);
  }

  if ((restrict==1) && (Spec[5] != 1))
    fprintf(fp_out,"\n   The power is restricted to be greater than or equal to 1");
  if ((restrict!=1) && (Spec[5] != 1))
    fprintf(fp_out,"\n   The power is not restricted");

  if(const_var == 1)
    fprintf(fp_out, "\n   A constant variance model is fit");
  else
    fprintf(fp_out, 
     "\n   The variance is to be modeled as Var(i) = exp(lalpha + log(mean(i)) * rho)");

  nparm_known = COUNT_SPECVECTOR(nparm, Spec);
  fprintf (fp_out, "\n\n   Total number of dose groups = %d",Nobs+Nmiss);
  fprintf (fp_out, "\n   Total number of records with missing values = %d",Nmiss);

  fprintf(fp_out, "\n   Maximum number of iterations = %d\n", ITMAX);
  fprintf(fp_out, "   Relative Function Convergence has been set to: %g\n", Rel_Conv);
  fprintf(fp_out, "   Parameter Convergence has been set to: %g\n\n", Parm_Conv);

  if ( (fabs(Rel_Conv - 1e-8) > 1e-10) || (fabs(Parm_Conv - 1e-8) > 1e-10) ) {
    fprintf(fp_out, "****  We are sorry but Relative Function and Parameter Convergence    ****\n");
    fprintf(fp_out, "****  are currently unavailable in this model.  Please keep checking  ****\n");
    fprintf(fp_out, "****  the web sight for model updates which will eventually           ****\n");
    fprintf(fp_out, "****  incorporate these convergence criterion.  Default values used.  ****\n\n");
  }

  if(initial==Yes) {
    OUTPUT_TEXT("\n\n                 User Inputs Initial Parameter Values  ");
    OUTPUT_Init(nparm, Spec, IniP, Parm_name);

    for (i=1; i<=nparm; i++){
      if(IniSp[i]==1) {   // have been initialized.
	if(Spec[i]==1 )  /* check if it is for fixed parm. */
	  Warning("The initial value for the fixed parameter is ignored.");
      }
      else {
	/* check if all the unspecified parms were initialized. */
	if (Spec[i]==0)
	  ERRORPRT("You have to initialize either ALL or NONE of the unspecified parameters.");
      }
    }

    if (IniP[1] <= 0)
      ERRORPRT("The initial value of variance has to be positive.");

    if (restrict == 1 && IniP[5] < 1)
      ERRORPRT("Initial value of power conflicts with sign restriction.");
  }

  /*compute likelihoods for A1, A2, and R*/
  lkA1 = lkA2 = lkA3 = lkR = 0.0;
  compute_continuous_liks(Nobs, Ni, Ym, Yd, &lkA1, &lkA2, &lkR);
  /* Compute Likelihood for model A3: Yij = Mu(i) + e(ij)
                                      Var(e(ij)) = k*(m(xi))^rho*/

  /* Calculate A3 likelihood regardless of constant variance choice. */
  /* MJF, 25MAY2005. */
  if (const_var == 0 && (Spec[2] == 0 || (Spec[2] == 1 && Parms[2] != 0.0))) { 
    // Parameters for fitting the model above
    LKParms = DVECTOR(1, Nobs+2);
    LKParms[1] = Parms[1];
    LKParms[2] = Parms[2];

    /******* Fit the A3 model above *********/
    nparm = Nobs+2;

    if (bNo_Log == false) {
      fprintf(fp_log, "\n********** Call to AThree_Fit ****************\n");
      fprintf(fp_log, "Variables going in:\n");
      fprintf(fp_log, "nparm = %2d;  EPS = %4g;  iter = %4d;  lkA3 = %6g\n", nparm, EPS, iter, lkA3);
      for (i=1; i<=Nobs+2; i++) {
	fprintf(fp_log,"LKParms[%d] = %13g\n",	i, LKParms[i]);
      }
    }
    AThree_Fit(nparm, LKParms, EPS, &iter, &lkA3);

    if (bNo_Log == false) {
      fprintf(fp_log, "\n********** Values out from AThree_Fit ****************\n");
      fprintf(fp_log, "Variables going in:\n");
      fprintf(fp_log, "nparm = %2d;  EPS = %4g;  iter = %4d;  lkA3 = %6g\n", nparm, EPS, iter, lkA3);
      for (i=1; i<=Nobs+2; i++) {
	fprintf(fp_log,"LKParms[%d] = %13g\n",	i, LKParms[i]);
      }
    }

    nparm = 5;

    FREE_DVECTOR(LKParms, 1, Nobs+2);
  }
  /* If constant variance model then just set equal to model A1. */
  /* MJF, 25MAY2005. */
  else {
    lkA3 = lkA1;
  }

  /* Get initial default adverse direction by finding the general
     linear trend of the data */
  if ((sign != 1) || (sign != -1)){
    sign = Get_Linear_Trend(Nobs, Xi, Ym, Ni);
    temp_sign = 0;
  }
  else {
    temp_sign = 9999;
  }/* end if ((sign != 1) || (sign != -1)) */

  if (bNo_Log == false) {
    fprintf(fp_log,"\nlkA1 = %4g          (Likelihood for model A1)\n", lkA1);
    fprintf(fp_log,"lkA2 = %4g          (Likelihood for model A2)\n", lkA2);
    fprintf(fp_log,"lkA3 = %4g          (Likelihood for model A3)\n", lkA3);
    fprintf(fp_log,"lkR = %4g           (Likelihood for model R)\n",lkR);
    fprintf(fp_log,"temp_sign = %4d     (This is 9999 if automatic is not selected and 0 if auto is selected)\n",temp_sign);
    fprintf(fp_log,"sign = %4d          (Sign is assigned here if auto is chosen)\n",sign);
  }

  bounded = IVECTOR(1, nparm);

  /*fitting Power model and output parameter estimates */
  if (bNo_Log == false) {
    fprintf(fp_log, "\n********** Call to Pow_Fit ****************\n");
    fprintf(fp_log, "Variables going in:\n");
    fprintf(fp_log, "nparm = %2d;  conv_check = %2d;  iter = %4d;  xlk = %6g\n", nparm, conv_check, iter, xlk);
    for (i=1; i<=5; i++) {
      fprintf(fp_log,"Parms[%d] = %13g          bounded[%d] = %2d\n",	i, Parms[i],i,bounded[i]);
    }
  }

  Pow_fit(nparm, Parms, EPS, &iter, &xlk, bounded, &conv_check);

  if (bNo_Log == false) {
    fprintf(fp_log, "Variables coming out of Pow_Fit:\n");
    fprintf(fp_log, "nparm = %2d;  conv_check = %2d;  iter = %4d;  xlk = %6g\n", nparm, conv_check, iter, xlk);
    fprintf(fp_log, "restrict = %2d\n", restrict);
    for (i=1; i<=5; i++) {
      fprintf(fp_log,"Parms[%d] = %13g          bounded[%d] = %2d\n",
	      i,Parms[i],i,bounded[i]);
    }
  }

  if (bNo_Log == false) {
    fprintf(fp_log, "\n******** bounded array after being corrected ************\n");
    for (i=1; i<=5; i++) {
      fprintf(fp_log,"bounded[%d] = %2d\n",i,bounded[i]);
    }
  }
  

  /* Check maximization results to get a better default adverse direction */
  if ((Parms[4] < 0) && (sign == 0))
    sign = -1;
  else
    if ((Parms[4] > 0) && (sign == 0))
      sign = 1;

  /* Compute the approx. covariance matrix */
  vcv = DMATRIX (1,nparm,1,nparm);

  INITIALIZE_DMATRIX(vcv, nparm, nparm);
  /* Creates matrix of second partial derivatives of the negative
     of the log likelihood function */

  Pow_vcv(nparm,Spec,Parms,vcv);

  /*  initialize adj_vcv matrix for Get_and_OUTPUTDTMSVCV() */

  adj_vcv_rows = 0;

  for (i = 1; i <= nparm; i++) {
    if (bounded[i] == 0) {
      adj_vcv_rows++;
    } /* end if */
  } /* end for */

  vcv_adj = DMATRIX(1, adj_vcv_rows, 1, adj_vcv_rows);

  /*output covariance matrix*/

  Get_and_OUTPUT_DTMSVCV(nparm, Spec, Parm_name, vcv, vcv_adj, bounded);

  OP_ParmsE(nparm,Spec,Parms,Parm_name,vcv_adj, bounded, bmdparm.level, 1);


  if (bNo_Log == false) {
    fprintf(fp_log, "\n******** bounded array after Get_and_OUTPUT_DTMSVCV ************\n");
    fprintf(fp_log, "\n*************** and going into OP_ParmsE ***********************\n");
    for (i=1; i<=5; i++) {
      fprintf(fp_log,"bounded[%d] = %2d\n",i,bounded[i]);
    }
  }

  /* Calculate fitted means and standard deviations, and print them out in a table with */
  /* the data, for comparison. */
  mean = DVECTOR (1, Nobs);
  std = DVECTOR (1, Nobs);
  for (i =1 ; i <= Nobs; i++)
    {
      if(Xi[i] == 0 && Parms[5] < 0)
	mean[i] = Parms[3] + Parms[4]*pow(EPS, Parms[5]);
      else
	mean[i] = Parms[3] + Parms[4]*pow(Xi[i], Parms[5]);
      if(Parms[2] == 0)
	std[i] = sqrt(Parms[1]);
      else
	std[i] = sqrt(exp(Parms[1] + log(fabs(mean[i])) * Parms[2]));
    } 

  PrintData(mean, std, Xi, Ym, Yd, Ni, Nobs);

  FREE_DVECTOR(mean, 1, Nobs);
  FREE_DVECTOR (std, 1, Nobs);

  /*compute and output ANOVA table elements*/
  var_type = -1*const_var + 1;

  if (bNo_Log == false) {
    fprintf(fp_log, "\n********** Call to DTMS3ANOVAC ****************\n");
    fprintf(fp_log, "Variables going in:\n");
    fprintf(fp_log, "nparm=%2d;  Nobs=%2d; lkA3=%6g; xlk=%6g; lkA2=%6g; lkA1=%6g; lkR=%6g; var_type=%2d;\n",
	    nparm, Nobs, lkA3, xlk, lkA2, lkA1, lkR, var_type);
    for (i=1; i<=5; i++) {
      fprintf(fp_log,"Spec[%2d] = %2d          bounded[%2d] = %2d\n",i,Spec[i],i,bounded[i]);
    }
  }

  /* Send 1 for var_type so we always assume modeled variance. */
  /* MJF, 20MAY2005. */
  DTMS3ANOVAC (nparm,Nobs,Spec,lkA3,xlk,lkA2, lkA1, lkR,anasum, 1, bounded);

  /* free memory */
  FREE_IVECTOR (bounded,1,nparm);

  if (bNo_Log == false) {
    fprintf(fp_log, "\n********** After Call to DTMS3ANOVAC ****************\n");
    fprintf(fp_log, "Variables coming out:\n");
    for (i=1; i<=5; i++) {
      fprintf(fp_log, "anasum[%2d].MSE=%12.5g    anasum[%2d].TEST=%f    anasum[%2d].DF=%d\n",
	      i,anasum[i].MSE,i,anasum[i].TEST,i,anasum[i].DF);
    }
  }

  /* MJF, 20MAY2005. */
  OUTPUT_DTMS3ANOVAC(anatxt,anasum, var_type); /*output ANOVA table*/

  /* MJF, 20MAY2005. */
  /*print a goodness of fit table*/
  Goodness(nparm, nparm_known, Parms, var_type, anasum);

  /*compute benchmark dose*/
#ifndef RBMDS
  fprintf (fp_out2, "\n BMD_flag \t %d \n Nobs \t%d \n nparm \t%d",  bmdose, Nobs, nparm );
  fprintf (fp_out2, "\n  Con_lev \t%3.3g ", bmdparm.level);
#endif
  sp = 0.0;
  Ntot = 0;
  if(Spec[2] == 1 && Parms[2] == 0) {
    for(i = 1; i <= Nobs; i++) {
      sp += (Ni[i] - 1)*Yd[i];
      Ntot += Ni[i];
    }
    sp = sp/(Ntot - Nobs);
  }

#ifndef RBMDS
  for (i=1;i<=nparm; i++)
    fprintf (fp_out2, "\n %s \t %5.3g", Parm_name[i-1], Parms[i]);

  fprintf (fp_out2,"\n\n Data");

  for (i=1;i<=Nobs;i++){
    Nd = Ni[i];
    if(Spec[2] == 1 && Parms[2] == 0) {
      lep=Ym[i] + qstudt(0.025, Nd - 1.) * sqrt(sp/Nd);
      upep=Ym[i] + qstudt(0.975, Nd - 1.) * sqrt(sp/Nd);
    }
    else {
      lep=Ym[i] + qstudt(0.025, Nd - 1.) * sqrt(Yd[i]/Nd);
      upep=Ym[i] + qstudt(0.975, Nd - 1.) * sqrt(Yd[i]/Nd);
    }

    fprintf (fp_out2,"\n %f %f %f %f", Xi[i], Ym[i], lep, upep);
  }

  fprintf (fp_out2,"\n Max_Min_dose \n  %f %f ", xmax, xmin);
#endif
  if (bmdose==Yes) {
    if(Spec[4]==Yes) {
#ifndef RBMDS
      fprintf (fp_out2, "\n\n BMDL_comput_ind %d",  No);
#endif
      fprintf (fp_out,"\n\n %s parameter is fixed.  ", Parm_name[3]);
      fprintf (fp_out, "The likelihood function can not be reparameterized in BMD.");
    } else /* end if */
      Pow_BMD (nparm, Parms, EPS, &junk, xlk, Rlevel, Bmdl, &BMD);

	//void Pow_BMD (int nparm, double p[], double gtol, int *iter, double xlk,
	//      double Rlevel[], double Bmdl[], double *BMD)	
  } /* end if (bmdose == Yes) */

  if (in_type == 1) {
    FREE_DVECTOR (Yd,1,Nobs);
    FREE_DVECTOR (stdev, 1, Nobs);
    FREE_DVECTOR (Ym,1,Nobs);
    FREE_DVECTOR (Xi,1,Nobs);
    FREE_IVECTOR (Ni,1,Nobs);
  } else {
    FREE_DVECTOR (Yd, 1, ntotal);
    FREE_DVECTOR (xxi, 1, ntotal);
    FREE_DVECTOR (yyi, 1, ntotal);
    FREE_DVECTOR (Ym, 1, ntotal);
    FREE_DVECTOR (Xi, 1, ntotal);
    FREE_IVECTOR (Ni, 1, ntotal);
    FREE_DVECTOR (Ysum, 1, ntotal);
  }

  FREE_DVECTOR (Parms,1,nparm);
  FREE_DVECTOR (ScYm,1,Nobs);

  FREE_DMATRIX (vcv_adj, 1, adj_vcv_rows, 1, adj_vcv_rows);

  FREE_DVECTOR (ScYd,1,Nobs);
  FREE_DVECTOR (IniP,1,nparm);
  FREE_IVECTOR (Spec,1,nparm);
  FREE_IVECTOR (IniSp,1,nparm);
  FREE_VLVECTOR(varsum,1,3);

  if(var_type == 0)
    FREE_ALVECTOR(anasum,1,4);
  else
    FREE_ALVECTOR(anasum,1,5);

  FREE_DVECTOR (Rlevel,1,5);
  FREE_DVECTOR (Bmdl, 1, 5);
  FREE_DMATRIX(vcv,1,nparm,1,nparm);
  if (bNo_Log == false) {
    fflush(fp_log);
    fclose(fp_log);
  }

  CLOSE_FILES ();
  return(0);

} /*end of main*/

/*******************************************************************
 **Pow_lk -- used to compute the log likelihood for Power model.
 * 		     Extern var.: smean, smax, Nobs, Xi, Yp, Yn, Ls, Xg.
 *
 *********************************************************************/

void Pow_lk(long int *nvar, double *x, long int *nf, double *f,
	    long int *uiparm, double *urparm, int (*ufparm)())
{

  int     i,j, jfixed,jvar;
  double   xlk;          		 /*  log likelihood */
  double   ex1, ex2, ex3, em, sig;        /* temp var */
  double   *p;

  /* parms for calculation */
  p=DVECTOR(1,nparm);
  jfixed=jvar=0;
  for(j=1; j<=nparm; j++)  /*reconstruct the parameter vector. */
    {
      if(Spec[j]==Yes)
	{
	  p[j]=urparm[jfixed];
	  jfixed++;
	}
      else
	{
	  p[j]=x[jvar];
	  jvar++;
	}
    }

  /* Restriction equation in BMDL calculation */
  if (replace==Yes)
    {
      if (bmdparm.risk==1)
	p[4]=BMR*pow(fabs(tD), -p[5]);
      else
	p[4]=p[3]*(BMR-1)*pow(fabs(tD), -p[5]);
    }

  /* To compute the log likelihood as a function of p[i] */

  xlk = 0.0;

  for (i=1;i<=Nobs;i++)
    {
      if (Xi[i] != 0)
	em=p[3]+p[4]*pow(fabs(Xi[i]), p[5]);
      else em=p[3];
      if (em != 0)
	sig=p[1]*pow(fabs(em), p[2]);
      else
	sig=0;
      ex1=0.5*Ni[i]*log(fabs(sig));
      ex2=(Ni[i]-1)*Yd[i]/(2.0*sig);
      ex3=Ni[i]*(Ym[i]-em)*(Ym[i]-em)/(2.0*sig);
      xlk -= ex1 + ex2 + ex3;
    }


  FREE_DVECTOR(p, 1, nparm);
  *f=-xlk;
}


/*******************************************************************
 **Pow_g -- used to compute the gradients for Power model.
 *		Extern var.: smean, smax, Nobs, Xi, Yp, Yn, Ls, Xg.
 *
 **********************************************************************/
void Pow_g(long int *nvar, double *x, long int *nf, double *g,
	   long int *uiparm, double *urparm, int (*ufparm)())
{
  int     i,j,k, jfixed, jvar;
  double  em;      /* temp var */
  double  sig, Ai, Bi, Hi, HHi;
  double  *dd, *p;		/* parms and gradients for calculation */

  dd=DVECTOR(1,nparm);
  p=DVECTOR(1,nparm);

  jfixed=jvar=0;

  if(BMDL_RUN == 1)
    {
      if(replace == Yes)
	{
	  Spec[2] = 1;
	  p[2] = 0.0;
	}
    }
  if(BMDL_RUN == 2)
    {
      if(replace == Yes)
	{
	  Spec[2] = 0;
	  Spec[3] = Spec[5] = Spec[4] = 1;
	}
    }

  if(BMDL_RUN == 3)
    {
      if(replace == Yes)
	Spec[3] = Spec[5] = Spec[4] = 0;
    }

  for(j=1; j<=nparm; j++)  /* reconstruct the parameter vector. */
    {
      if(Spec[j]==Yes)
	{
	  p[j]=urparm[jfixed];
	  jfixed++;
	}
      else
	{
	  p[j]=x[jvar];
	  jvar++;
	}
    }


  /* Restriction equation in BMDL calculation */
  if (replace==Yes)
    {
      if (bmdparm.risk==1)
	p[4]=BMR*pow(fabs(tD), -p[5]);
      else
	p[4]=p[3]*(BMR-1)*pow(fabs(tD), -p[5]);
    }

  /* initial g[j]s */
  for (j=1;j<=nparm;j++) g[j]=dd[j]=0.0;

  /* Computation of first partial derivatives */
  for (i=1; i<=Nobs; i++)
    {
      if (Xi[i] != 0)
	em = p[3]+p[4]*pow(fabs(Xi[i]), p[5]);
      else em=p[3];

      if (em != 0)
	sig = p[1]*pow(fabs(em), p[2]);
      else sig=0;

      Ai=(Ni[i]-1)*Yd[i]+Ni[i]*Ym[i]*Ym[i];
      Bi=Ni[i]*Ym[i];

      Hi=Ai/(2*sig) - Bi*em/(sig) + Ni[i]*em*em/(2*sig);
      HHi=-p[2]*Ai/(2*sig*em) + (p[2]-1)*Bi/sig - (p[2]-2)*Ni[i]*em/(2*sig);


      g[1] += Ni[i]/(2*p[1]) - Hi/p[1];
      g[2] += (0.5*Ni[i]-Hi)*log(fabs(em));

      dd[3]=1.0;
      if (Xi[i] != 0)
	{
	  dd[4]=pow(fabs(Xi[i]), p[5]);
	  dd[5]=p[4]*pow(fabs(Xi[i]), p[5])*log(fabs(Xi[i]));
	}
      else
	dd[4]=dd[5]=0.0;
      if (replace == Yes)
	{
	  if (bmdparm.risk==1 && tD!=0)
	    dd[5] -= dd[4]*BMR*pow(fabs(tD), -p[5])*log(fabs(tD));
	  else if (bmdparm.risk==0 && tD !=0)
	    {
	      dd[3] += dd[4]*(BMR-1)/pow(fabs(tD), p[5]);
	      dd[5] -= dd[4]*p[3]*(BMR-1)*log(fabs(tD))/pow(fabs(tD), p[5]);
	    }
	}

      for (k=3; k<=nparm; k++)
	g[k] += dd[k]*(p[2]*Ni[i]/(2*em) + HHi);


    }


  /* end of first partial deri. */

  jvar=0;
  for(j=1; j<=nparm; j++)  /* reconstruct the parameter vector. */
    if(Spec[j]==No)
      {
	g[jvar]=g[j];
	jvar++;
      }

  FREE_DVECTOR(dd, 1,nparm);
  FREE_DVECTOR(p, 1,nparm);
  return;
}


/******************************************************************
 *	Pow_vcv -- used to compute the vcv for Power model.
 *		Extern var.: Nobs, Xi, Yp, Yn, Ls, Xg.
 *
 ******************************************************************/
void
Pow_vcv (int nparm, int Spec[], double p[], double **vcv)
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

/*******************************************************************
 *	First partial derivatives of the mean function at observation obs,
 *	contained in mg[] at exit.  This is for Power model.
 *******************************************************************/
void
MeanPart (int obs, double *p, double *mg)
{
  mg[1] = mg[2] = 0.0;		/* Variance parameters not in mean function */
  mg[3] = 1.0;

  if (Xi[obs] != 0)
    {
      mg[4] = pow (Xi[obs], p[5]);
      mg[5] = Slog(Xi[obs])*p[4]*mg[4];
    }
  else
    {
      mg[4] = mg[5] = 0.0;
    }

}				/* end MeanPart */
/********************************************************************
 *	Second partial derivates with respect to the mean function.
 *	mg2[][] contains all second partials upon exit.
 *      This is for Power model.
 ********************************************************************/
void
Mean2Part (int obs, double *p, double **mg2)
{
  int j;
  double Drho, lD;
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
      Drho = pow(Xi[obs], p[5]);
      lD = log(Xi[obs]);
      mg2[4][5] = lD * Drho;
      mg2[5][5] = p[4]*lD*lD*Drho;
    }
  else
    {
      mg2[4][5] = mg2[5][5] = 0.0;
    }

  mg2[5][4] = mg2[4][5];

}				/* end Mean2Part */
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
      meani = p[3] + p[4] * pow (Xi[obs], p[5]);
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
      meani = p[3] + p[4] * pow (Xi[obs], p[5]);
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
       meani = p[3] + p[4] * pow (Xi[obs], p[5]);
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

/**************************************************************
 *Pow_fit -- Used to "prepare" the data for further computation,
 *            i.e. compute the extern variables, give the initial
 *            parameters, etc. THEN fit the Power model.
 *            (In fact, these jobs could be done in main().)
 *
 ***************************************************************/
void Pow_fit(int nparm, double p[], double gtol,
	     int *iter, double *fret, int *bounded, int *is_conv)
{
	int   i, j, junk, count_neg, count_max, ii, key, size, count_min;

	double *pBak, /**tmy, *t, **tmv,*/ *Y_Less_Bg, *LogDose;
	double est1, est2, pe, dfpe, **x, **xp, **xpx, *betav, *xpy;
	double maxdose, *doses, *means, *svar, *parms, *fitparms, *fitparms2;
	double *fitparms3, *beginp, ll, ll2, ll3;
	int inc_flag, dec_flag, gmcnt;
	long int *nanim, *Spec2, *bind, *bind2, *bind3, optite, optite2, optite3;
	long int nresm, model_type, flag;
	long int nvar, signs, nparms, restr;
	double maxYm, maxYd, *mean, *var, LL;


	if(bNo_Log == false) {
		fprintf(fp_log,"The variables inside Pow_fit\n");
		fprintf(fp_log,"nparm = %d\n",nparm);
		fprintf(fp_log,"sign = %d\n",sign);
		fflush(fp_log);
	}

	pBak=DVECTOR(1, nparm);

	xpx = DMATRIX(1,2,1,2);
	betav = DVECTOR(1,2);
	xpy = DVECTOR(1,2);

	model_type = 1;
	pe = 0; /* this was added to get correct initial values */
	dfpe = 0;

	/** rescale Dose to be: 0 <= Dose <= 1 **/
	scale=1;

	for (j=1; j<=nparm; j++)
		pBak[j]=p[j];    /* save the input p[]. */

	/****** Obtain initial estimations for p[] ******/
	if(initial==Yes)
	{
		for(j=1; j<=nparm; j++)
			p[j]=IniP[j];
		size = 0;
	}
	else
	{
		/*compute initial estimates*/  /* New */
		est1=est2=0.0;
		if (Spec[1] == 0) {
			for(i = 1; i <= Nobs; i++)
			{
				pe += (Ni[i]-1)*Yd[i];
				dfpe += Ni[i]-1;

			}

			if (const_var == 0)
				p[1] = log(pe/dfpe);
			else 
				p[1] = pe/dfpe;

		} 
		if (Spec[2] == 0) {
			p[2]=0;
		}

		count_neg = 0;
		junk = 0;
		ii = 2;
		count_min = 0;
		count_max = 0;
		key = 1;


		p[3] = Ym[1];

		if(restrict == -1)
			key = -1;

		for(i = 2; i <= Nobs; i++)
		{
			if(Ym[i] - p[3] < 0)
				count_neg++;
		}

		for(i = 2; i <= Nobs; i++)
		{
			if(Ym[i] == yymax)
				count_max++;
			if(Ym[i] == yymin)
				count_min++;
		}

		if(restrict == 0)
		{
			if(count_neg > Nobs - count_neg - count_max)
				key = -1;
			if(count_neg == (Nobs - count_neg- count_max))
				if(Ym[Nobs] < Ym[1])
					key = -1;
		}

		if(key == 1)
		{
			size = Nobs - 1 - count_min;
			p[3] = yymin;
		}
		else
		{
			size = Nobs - 1 - count_max;
			p[3] = yymax;
		}


		if(Nobs <= 1 || size < 0)
			ERRORPRT("Error in selecting default parameters\n");

		if(size == 0)
		{
			p[4] = 0.0;
			p[5] = 1.0;
		}

		if(size == 1)
		{
			p[5] = 1.0;
			for(i = 2; i <= Nobs; i++)
			{
				if(Ym[i] != p[3])
				{
					p[4] = (Ym[i] - p[3])/Xi[i];
					break;
				}
			}
		}
		if(size >= 2)
		{
			Y_Less_Bg=DVECTOR(1, size);
			LogDose=DVECTOR(1, size);
			x = DMATRIX(1, size, 1, 2);
			xp = DMATRIX(1, 2, 1, size);

			i = 1;

			while (i<= size)
			{
				if((Ym[ii]-p[3]) == 0)
					ii++;
				else
				{
					Y_Less_Bg[i] = key*(Ym[ii] - p[3]);
					LogDose[i] = log(Xi[ii]);
					Y_Less_Bg[i] = log(Y_Less_Bg[i]);
					ii++;
					i++;
				}
			}

			for(i = 1; i <= size; i++)
			{
				x[i][1] = 1.0;
				x[i][2] = LogDose[i];
			}

			TRANSPOSE(x, xp, size, 2);

			MATMPYM2(xp, x, xpx, 2, size, 2);

			MATMPYV2(2, size, xp, Y_Less_Bg, xpy);

			INVMAT(xpx, 2);

			MATMPYV2(2,2,xpx, xpy, betav);

			p[5] = betav[2];

			p[4] = key*exp(betav[1]);

			FREE_DMATRIX(x, 1, size, 1, 2);
			FREE_DMATRIX(xp, 1, 2, 1, size);

			FREE_DVECTOR(Y_Less_Bg, 1, size);
			FREE_DVECTOR(LogDose, 1, size);

		} /* end  if(size >= 2) */

		inc_flag = 1;
		dec_flag = 1;
		i = 1;
		while ((i<Nobs) && (inc_flag == 1))
		{
			if (Ym[i+1] < Ym[i])
				inc_flag = 0;
			i++;
		}
		i = 1;
		while ((i<Nobs) && (dec_flag == 1))
		{
			if (Ym[i+1] > Ym[i])
				dec_flag = 0;
			i++;
		}

		if ((inc_flag == 1) && (dec_flag == 0))
		{
			if (Spec[3] != 1)
				p[3] = Ym[1];
			if (Spec[4] != 1)
				p[4] = (Ym[Nobs]-Ym[1])/pow(Xi[Nobs],p[5]);
		}
		if ((inc_flag == 0) && (dec_flag == 1))
		{
			if (Spec[3] != 1)
				p[3] = Ym[1];
			if (Spec[4] != 1)
				p[4] = (Ym[Nobs]-Ym[1])/pow(Xi[Nobs],p[5]);
		}


		OUTPUT_TEXT("\n\n                  Default Initial Parameter Values  ");
		OUTPUT_Init(nparm, Spec, p, Parm_name);
	} /* end if(initial==Yes) */

	FREE_DMATRIX(xpx, 1, 2, 1, 2);
	FREE_DVECTOR(betav, 1, 2);
	FREE_DVECTOR(xpy, 1, 2);

	/* the code below this point was taken from  Hill_fit */
	replace = No;

	maxdose = Xi[1];
	maxYm = Ym[1];
	maxYd = Yd[1];
	for (i = 2; i <= Nobs; i++)
	{
		if (maxdose < Xi[i])
			maxdose = Xi[i];
		if (maxYm < Ym[i])
			maxYm = Ym[i];
		if (maxYd < Yd[i])
			maxYd = Yd[i];
	}

	for (i = 1; i <= Nobs; i++)
	{
		Xi[i] = Xi[i] / maxdose;
	}

	/* get specified parameters */
	for (j = 1; j <= nparm; j++)
	{
		if (Spec[j] == Yes)
		{
			p[j]=pBak[j];
		}
	}
	nvar = Nobs;
	nparms = nparm;
	restr = restrict;
	signs = sign;

	doses = DVECTOR(0, Nobs-1);
	means = DVECTOR(0, Nobs-1);
	svar =  DVECTOR(0, Nobs-1);
	nanim = LIVECTOR(0, Nobs-1);
	parms = DVECTOR(0, nparm-1);
	fitparms = DVECTOR(0, nparm-1);
	fitparms2 = DVECTOR(0, nparm-1);
	fitparms3 = DVECTOR(0, nparm-1);
	Spec2 = LIVECTOR(0, nparm-1);
	bind = LIVECTOR(0, nparm-1);
	bind2 = LIVECTOR(0, nparm-1);
	bind3 = LIVECTOR(0, nparm-1);
	beginp = DVECTOR(0, nparm-1);

	for (i = 1; i <= Nobs; i++)
	{
		nanim[i-1] = (long int) Ni[i];
		doses[i-1] = Xi[i];
		means[i-1] = Ym[i]/maxYm;
		svar[i-1] = Yd[i]/pow(maxYm,2);
	}

	for (i = 1; i <= nparm; i++)
	{
		if ((initial == Yes) && (Spec[i] == 0))
		{
			p[i] = IniP[i];
			Spec2[i-1] = 0;
		}
		else
		{
			Spec2[i-1] = Spec[i];
		}
		parms[i-1] = p[i];
	}	/* end for (i = 1; i <= nparm; i++) */

	/* rescale the slope */
	parms[3] = parms[3]*pow(maxdose, parms[4]);

	/* Scale parms by max mean and variance */
	if (parms[1] == 0)
		parms[0] = parms[0]/pow(maxYm,2);
	else
		parms[0] = parms[0] + (parms[1] - 2.0)*log(maxYm);
	parms[2] = parms[2]/maxYm;
	parms[3] = parms[3]/maxYm;

	for (i = 0; i < nparm; i++)
	{
		beginp[i] = parms[i];
	}
	flag = -1;

	if (bNo_Log == false)
	{
		fprintf(fp_log,"\n***********First call to getmle.****************\n");
		fprintf(fp_log,"(this call scales the parameters by their starting values in donlp2)\n");
		fprintf(fp_log,"nparms = %ld\n",nparms);
		fprintf(fp_log,"flag = %ld\n",flag);
		fprintf(fp_log,"These are the parameters going into getmle.  The slope\n");
		fprintf(fp_log,"has been scaled by maxdose^power.\n");
		fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
		fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
		fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
		fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
		fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
		fprintf(fp_log,"************************************************\n");
		fflush(fp_log);
	}

	/* This is the first call to getmle.  The parameters will be scaled */
	/* internally by donlp2 by their starting values.                   */
	gmcnt = 1;
	getmle_(&nvar, doses, means, nanim, svar, &nparms, beginp, Spec2, beginp,
		&restr, &signs, fitparms2, &ll2, &optite2, &nresm, bind2,
		&model_type, &flag);

	if (bNo_Log == false)
	{
		fprintf(fp_log,"\n*******After the first call to getmle.**********\n");
		fprintf(fp_log,"nparms = %ld\n",nparms);
		fprintf(fp_log,"flag = %ld\n",flag);
		fprintf(fp_log,"ll2 = %10.5g   (likelihood)\n",ll2);
		fprintf(fp_log,"optite2 = %ld    (good optimum if 0<=optite2<=2)\n",optite2);
		fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
		fprintf(fp_log,"These are the parameters coming out of the first run\n");
		fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");
		fprintf(fp_log,"fitparms2[0] = %12.5g  (alpha)   bind[0] = %ld\n",fitparms2[0], bind2[0]);
		fprintf(fp_log,"fitparms2[1] = %12.5g  (rho)     bind[1] = %ld\n",fitparms2[1], bind2[1]);
		fprintf(fp_log,"fitparms2[2] = %12.5g  (control) bind[2] = %ld\n",fitparms2[2], bind2[2]);
		fprintf(fp_log,"fitparms2[3] = %12.5g  (slope)   bind[3] = %ld\n",fitparms2[3], bind2[3]);
		fprintf(fp_log,"fitparms2[4] = %12.5g  (power)   bind[4] = %ld\n",fitparms2[4], bind2[4]);
		fprintf(fp_log,"*************************************************\n");
		fflush(fp_log);
	}

	if ((optite2 < 0) || (optite2 > 2))
	{
		for (i = 0; i < 4; i++)
		{
			GetNewParms(beginp, nparm);

			gmcnt = gmcnt + 1;
			if (bNo_Log == false)
			{
				fprintf(fp_log,"\n***********Call #%d to getmle.****************\n",gmcnt);
				fprintf(fp_log,"(this call scales the parameters by their starting values in donlp2)\n");
				fprintf(fp_log,"nparms = %ld\n",nparms);
				fprintf(fp_log,"flag = %ld\n",flag);
				fprintf(fp_log,"These are the parameters going into getmle.  The slope\n");
				fprintf(fp_log,"has been scaled by maxdose^power.\n");
				fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
				fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
				fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
				fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
				fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
				fprintf(fp_log,"************************************************\n");
				fflush(fp_log);
			}

			getmle_(&nvar, doses, means, nanim, svar, &nparms, beginp, Spec2, beginp,
				&restr, &signs, fitparms2, &ll2, &optite2, &nresm, bind2,
				&model_type, &flag);

			if (bNo_Log == false)
			{
				fprintf(fp_log,"\n*******After Call #%d to getmle.**********\n",gmcnt);
				fprintf(fp_log,"nparms = %ld\n",nparms);
				fprintf(fp_log,"flag = %ld\n",flag);
				fprintf(fp_log,"ll2 = %10.5g   (likelihood)\n",ll2);
				fprintf(fp_log,"optite2 = %ld    (good optimum if 0<=optite2<=2)\n",optite2);
				fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
				fprintf(fp_log,"These are the parameters coming out of the first run\n");
				fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");
				fprintf(fp_log,"fitparms2[0] = %12.5g  (alpha)\n",fitparms2[0]);
				fprintf(fp_log,"fitparms2[1] = %12.5g  (rho)\n",fitparms2[1]);
				fprintf(fp_log,"fitparms2[2] = %12.5g  (control)\n",fitparms2[2]);
				fprintf(fp_log,"fitparms2[3] = %12.5g  (slope)\n",fitparms2[3]);
				fprintf(fp_log,"fitparms2[4] = %12.5g  (power)\n",fitparms2[4]);
				fprintf(fp_log,"*************************************************\n");
				fflush(fp_log);
			}

			if ((optite2 >= 0) && (optite2 <= 2))
			{
				break;
			}
		}
	}

	if ((optite2 < 0) || (optite2 > 2))
	{
		for (i = 0; i < 4; i++)
		{
			GetMoreParms(beginp, nparm);

			gmcnt = gmcnt + 1;
			if (bNo_Log == false)
			{
				fprintf(fp_log,"\n***********Call #%d to getmle.****************\n",gmcnt);
				fprintf(fp_log,"(this call scales the parameters by their starting values in donlp2)\n");
				fprintf(fp_log,"nparms = %ld\n",nparms);
				fprintf(fp_log,"flag = %ld\n",flag);
				fprintf(fp_log,"These are the parameters going into getmle.  The slope\n");
				fprintf(fp_log,"has been scaled by maxdose^power.\n");
				fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
				fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
				fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
				fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
				fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
				fprintf(fp_log,"************************************************\n");
				fflush(fp_log);
			}

			getmle_(&nvar, doses, means, nanim, svar, &nparms, beginp, Spec2, beginp,
				&restr, &signs, fitparms2, &ll2, &optite2, &nresm, bind2,
				&model_type, &flag);

			if (bNo_Log == false)
			{
				fprintf(fp_log,"\n*******After Call #%d to getmle.**********\n",gmcnt);
				fprintf(fp_log,"nparms = %ld\n",nparms);
				fprintf(fp_log,"flag = %ld\n",flag);
				fprintf(fp_log,"ll2 = %10.5g   (likelihood)\n",ll2);
				fprintf(fp_log,"optite2 = %ld    (good optimum if 0<=optite2<=2)\n",optite2);
				fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
				fprintf(fp_log,"These are the parameters coming out of the first run\n");
				fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");
				fprintf(fp_log,"fitparms2[0] = %12.5g  (alpha)\n",fitparms2[0]);
				fprintf(fp_log,"fitparms2[1] = %12.5g  (rho)\n",fitparms2[1]);
				fprintf(fp_log,"fitparms2[2] = %12.5g  (control)\n",fitparms2[2]);
				fprintf(fp_log,"fitparms2[3] = %12.5g  (slope)\n",fitparms2[3]);
				fprintf(fp_log,"fitparms2[4] = %12.5g  (power)\n",fitparms2[4]);
				fprintf(fp_log,"*************************************************\n");
				fflush(fp_log);
			}

			if ((optite2 >= 0) && (optite2 <= 2))
			{
				break;
			}
		}
	}

	if ((optite2 < 0) || (optite2 > 2))
	{
		for (i = 0; i < 4; i++)
		{
			GetMLEParms(beginp, nparm);

			gmcnt = gmcnt + 1;
			if (bNo_Log == false)
			{
				fprintf(fp_log,"\n***********Call #%d to getmle.****************\n",gmcnt);
				fprintf(fp_log,"(this call scales the parameters by their starting values in donlp2)\n");
				fprintf(fp_log,"nparms = %ld\n",nparms);
				fprintf(fp_log,"flag = %ld\n",flag);
				fprintf(fp_log,"These are the parameters going into getmle.  The slope\n");
				fprintf(fp_log,"has been scaled by maxdose^power.\n");
				fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
				fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
				fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
				fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
				fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
				fprintf(fp_log,"************************************************\n");
				fflush(fp_log);
			}

			getmle_(&nvar, doses, means, nanim, svar, &nparms, beginp, Spec2, beginp,
				&restr, &signs, fitparms2, &ll2, &optite2, &nresm, bind2,
				&model_type, &flag);

			if (bNo_Log == false)
			{
				fprintf(fp_log,"\n*******After Call #%d to getmle.**********\n",gmcnt);
				fprintf(fp_log,"nparms = %ld\n",nparms);
				fprintf(fp_log,"flag = %ld\n",flag);
				fprintf(fp_log,"ll2 = %10.5g   (likelihood)\n",ll2);
				fprintf(fp_log,"optite2 = %ld    (good optimum if 0<=optite2<=2)\n",optite2);
				fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
				fprintf(fp_log,"These are the parameters coming out of the first run\n");
				fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");
				fprintf(fp_log,"fitparms2[0] = %12.5g  (alpha)\n",fitparms2[0]);
				fprintf(fp_log,"fitparms2[1] = %12.5g  (rho)\n",fitparms2[1]);
				fprintf(fp_log,"fitparms2[2] = %12.5g  (control)\n",fitparms2[2]);
				fprintf(fp_log,"fitparms2[3] = %12.5g  (slope)\n",fitparms2[3]);
				fprintf(fp_log,"fitparms2[4] = %12.5g  (power)\n",fitparms2[4]);
				fprintf(fp_log,"*************************************************\n");
				fflush(fp_log);
			}

			if ((optite2 >= 0) && (optite2 <= 2))
			{
				break;
			}
		}
	}

	flag = 0;

	/* This is at most the 14th call to getmle.  The initial parameters will be */
	/* the fit parameters from the scaled (previous) call                       */
	gmcnt = gmcnt + 1;
	if (bNo_Log == false)
	{
		fprintf(fp_log,"\n***********Call #%d to getmle.****************\n",gmcnt);
		fprintf(fp_log,"(This call doesn't scale the parameters by their starting values in\n");
		fprintf(fp_log," donlp2.  This call is using the parms coming out of the previous calls.)\n");
		fprintf(fp_log,"nparms = %ld\n",nparms);
		fprintf(fp_log,"flag = %ld\n",flag);
		fprintf(fp_log,"These are the parameters going into getmle.  The slope\n");
		fprintf(fp_log,"has been scaled by maxdose^power.\n");
		fprintf(fp_log,"fitparms2[0] = %12.5g  (alpha)\n",fitparms2[0]);
		fprintf(fp_log,"fitparms2[1] = %12.5g  (rho)\n",fitparms2[1]);
		fprintf(fp_log,"fitparms2[2] = %12.5g  (control)\n",fitparms2[2]);
		fprintf(fp_log,"fitparms2[3] = %12.5g  (slope)\n",fitparms2[3]);
		fprintf(fp_log,"fitparms2[4] = %12.5g  (power)\n",fitparms2[4]);
		fprintf(fp_log,"************************************************\n");
		fflush(fp_log);
	}

	getmle_(&nvar, doses, means, nanim, svar, &nparms, fitparms2, Spec2, beginp,
		&restr, &signs, fitparms3, &ll3, &optite3, &nresm, bind3,
		&model_type, &flag);

	if (bNo_Log == false)
	{
		fprintf(fp_log,"\n*******After Call #%d to getmle.**********\n",gmcnt);
		fprintf(fp_log,"nparms = %ld\n",nparms);
		fprintf(fp_log,"flag = %ld\n",flag);
		fprintf(fp_log,"ll3 = %10.5g   (likelihood)\n",ll3);
		fprintf(fp_log,"optite3 = %ld    (good optimum if 0<=optite3<=2)\n",optite3);
		fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
		fprintf(fp_log,"These are the parameters coming out of the first run\n");
		fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");
		fprintf(fp_log,"fitparms3[0] = %12.5g  (alpha)\n",fitparms3[0]);
		fprintf(fp_log,"fitparms3[1] = %12.5g  (rho)\n",fitparms3[1]);
		fprintf(fp_log,"fitparms3[2] = %12.5g  (control)\n",fitparms3[2]);
		fprintf(fp_log,"fitparms3[3] = %12.5g  (slope)\n",fitparms3[3]);
		fprintf(fp_log,"fitparms3[4] = %12.5g  (power)\n",fitparms3[4]);
		fprintf(fp_log,"*************************************************\n");
		fflush(fp_log);
	}


	/* This starts the loops without scaling in donlp2 */
	for (flag=0; flag<=1; flag++)
	{ /* flag changes the value of delta in powmeans */

		gmcnt = gmcnt + 1;
		if (bNo_Log == false)
		{
			fprintf(fp_log,"\n***********Call #%d to getmle.****************\n",gmcnt);
			fprintf(fp_log,"(this call doesn't scale the parameters by their starting values in donlp2)\n");
			fprintf(fp_log,"nparms = %ld\n",nparms);
			fprintf(fp_log,"flag = %ld\n",flag);
			fprintf(fp_log,"These are the parameters going into getmle.  The slope\n");
			fprintf(fp_log,"has been scaled by maxdose^power.\n");
			fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
			fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
			fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
			fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
			fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
			fprintf(fp_log,"************************************************\n");
			fflush(fp_log);
		}

		getmle_(&nvar, doses, means, nanim, svar, &nparms, beginp, Spec2, beginp,
			&restr, &signs, fitparms, &ll, &optite, &nresm, bind,
			&model_type, &flag);

		if (bNo_Log == false)
		{
			fprintf(fp_log,"\n*******After Call #%d to getmle.**********\n",gmcnt);
			fprintf(fp_log,"nparms = %ld\n",nparms);
			fprintf(fp_log,"flag = %ld\n",flag);
			fprintf(fp_log,"ll = %10.5g   (likelihood)\n",ll);
			fprintf(fp_log,"optite = %ld    (good optimum if 0<=optite<=2)\n",optite);
			fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
			fprintf(fp_log,"These are the parameters coming out of the first run\n");
			fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");
			fprintf(fp_log,"fitparms[0] = %12.5g  (alpha)\n",fitparms[0]);
			fprintf(fp_log,"fitparms[1] = %12.5g  (rho)\n",fitparms[1]);
			fprintf(fp_log,"fitparms[2] = %12.5g  (control)\n",fitparms[2]);
			fprintf(fp_log,"fitparms[3] = %12.5g  (slope)\n",fitparms[3]);
			fprintf(fp_log,"fitparms[4] = %12.5g  (power)\n",fitparms[4]);
			fprintf(fp_log,"*************************************************\n");
			fflush(fp_log);
		}

		if ((optite >= 0) && (optite <=2))
		{
			flag = 2;
			break;
		}

		if ((optite < 0) || (optite > 2))
		{
			for (i = 0; i < 30; i++)
			{
				if (optite != 3)
				{
					GetNewParms(beginp, nparm);	/* Get a new starting point */
				}
				else
				{
					for (j = 0; j < nparm; j++)
					{
						beginp[j] = fitparms[j];
					}
				}

				/* Try again */
				gmcnt = gmcnt + 1;
				if (bNo_Log == false)
				{
					fprintf(fp_log,"\n***********Call #%d to getmle.****************\n",gmcnt);
					fprintf(fp_log,"(this call doesn't scale the parameters by their starting values in donlp2)\n");
					fprintf(fp_log,"nparms = %ld\n",nparms);
					fprintf(fp_log,"flag = %ld\n",flag);
					fprintf(fp_log,"These are the parameters going into getmle.  The slope\n");
					fprintf(fp_log,"has been scaled by maxdose^power.\n");
					fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
					fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
					fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
					fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
					fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
					fprintf(fp_log,"************************************************\n");
					fflush(fp_log);
				}

				getmle_(&nvar, doses, means, nanim, svar, &nparms, beginp, Spec2, beginp,
					&restr, &signs, fitparms, &ll, &optite, &nresm, bind,
					&model_type, &flag);

				if (bNo_Log == false)
				{
					fprintf(fp_log,"\n*******After Call #%d to getmle.**********\n",gmcnt);
					fprintf(fp_log,"nparms = %ld\n",nparms);
					fprintf(fp_log,"flag = %ld\n",flag);
					fprintf(fp_log,"ll = %10.5g   (likelihood)\n",ll);
					fprintf(fp_log,"optite = %ld    (good optimum if 0<=optite<=2)\n",optite);
					fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
					fprintf(fp_log,"These are the parameters coming out of the first run\n");
					fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");
					fprintf(fp_log,"fitparms[0] = %12.5g  (alpha)\n",fitparms[0]);
					fprintf(fp_log,"fitparms[1] = %12.5g  (rho)\n",fitparms[1]);
					fprintf(fp_log,"fitparms[2] = %12.5g  (control)\n",fitparms[2]);
					fprintf(fp_log,"fitparms[3] = %12.5g  (slope)\n",fitparms[3]);
					fprintf(fp_log,"fitparms[4] = %12.5g  (power)\n",fitparms[4]);
					fprintf(fp_log,"*************************************************\n");
					fflush(fp_log);
				}

				if ((optite >= 0) && (optite <= 2))
				{
					flag = 2;
					break;
				}
			}
		}	/* end if ((optite < 0) || (optite > 2)) */


		if ((optite < 0) || (optite > 2))
		{
			for (i = 0; i < 30; i++)
			{
				if (optite != 3)
				{
					GetMoreParms(beginp, nparm);		/* Get a new starting point */
				}
				else
				{
					for (j = 0; j < nparm; j++)
					{
						beginp[j] = fitparms[j];
					}
				}

				/* Try again */
				gmcnt = gmcnt + 1;
				if (bNo_Log == false)
				{
					fprintf(fp_log,"\n***********Call #%d to getmle.****************\n",gmcnt);
					fprintf(fp_log,"(this call doesn't scale the parameters by their starting values in donlp2)\n");
					fprintf(fp_log,"nparms = %ld\n",nparms);
					fprintf(fp_log,"flag = %ld\n",flag);
					fprintf(fp_log,"These are the parameters going into getmle.  The slope\n");
					fprintf(fp_log,"has been scaled by maxdose^power.\n");
					fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
					fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
					fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
					fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
					fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
					fprintf(fp_log,"************************************************\n");
					fflush(fp_log);
				}

				getmle_(&nvar, doses, means, nanim, svar, &nparms, beginp, Spec2, beginp,
					&restr, &signs, fitparms, &ll, &optite, &nresm, bind,
					&model_type, &flag);

				if (bNo_Log == false)
				{
					fprintf(fp_log,"\n*******After Call #%d to getmle.**********\n",gmcnt);
					fprintf(fp_log,"nparms = %ld\n",nparms);
					fprintf(fp_log,"flag = %ld\n",flag);
					fprintf(fp_log,"ll = %10.5g   (likelihood)\n",ll);
					fprintf(fp_log,"optite = %ld    (good optimum if 0<=optite<=2)\n",optite);
					fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
					fprintf(fp_log,"These are the parameters coming out of the first run\n");
					fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");
					fprintf(fp_log,"fitparms[0] = %12.5g  (alpha)\n",fitparms[0]);
					fprintf(fp_log,"fitparms[1] = %12.5g  (rho)\n",fitparms[1]);
					fprintf(fp_log,"fitparms[2] = %12.5g  (control)\n",fitparms[2]);
					fprintf(fp_log,"fitparms[3] = %12.5g  (slope)\n",fitparms[3]);
					fprintf(fp_log,"fitparms[4] = %12.5g  (power)\n",fitparms[4]);
					fprintf(fp_log,"*************************************************\n");
					fflush(fp_log);
				}

				if ((optite >= 0) && (optite <= 2))
				{
					flag = 2;
					break;
				}
			}
		}	/* end if ((optite < 0) || (optite > 2)) */

		if ((optite < 0) || (optite > 2))
		{
			for (i = 0; i < nparm; i++)
			{
				parms[i] = p[i+1];
			}
			parms[3] = parms[3]*pow(maxdose, parms[4]);
			for (i = 0; i < nparm; i++)
			{
				beginp[i] = parms[i];
			}

			gmcnt = gmcnt + 1;
			if (bNo_Log == false)
			{
				fprintf(fp_log,"\n***********Call #%d to getmle.****************\n",gmcnt);
				fprintf(fp_log,"(this call doesn't scale the parameters by their starting values in donlp2)\n");
				fprintf(fp_log,"nparms = %ld\n",nparms);
				fprintf(fp_log,"flag = %ld\n",flag);
				fprintf(fp_log,"These are the parameters going into getmle.  The slope\n");
				fprintf(fp_log,"has been scaled by maxdose^power.\n");
				fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
				fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
				fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
				fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
				fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
				fprintf(fp_log,"************************************************\n");
				fflush(fp_log);
			}

			getmle_(&nvar, doses, means, nanim, svar, &nparms, parms, Spec2, parms,
				&restr, &signs, fitparms, &ll, &optite, &nresm, bind,
				&model_type, &flag);

			if (bNo_Log == false)
			{
				fprintf(fp_log,"\n*******After Call #%d to getmle.**********\n",gmcnt);
				fprintf(fp_log,"nparms = %ld\n",nparms);
				fprintf(fp_log,"flag = %ld\n",flag);
				fprintf(fp_log,"ll = %10.5g   (likelihood)\n",ll);
				fprintf(fp_log,"optite = %ld    (good optimum if 0<=optite<=2)\n",optite);
				fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
				fprintf(fp_log,"These are the parameters coming out of the first run\n");
				fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");
				fprintf(fp_log,"fitparms[0] = %12.5g  (alpha)\n",fitparms[0]);
				fprintf(fp_log,"fitparms[1] = %12.5g  (rho)\n",fitparms[1]);
				fprintf(fp_log,"fitparms[2] = %12.5g  (control)\n",fitparms[2]);
				fprintf(fp_log,"fitparms[3] = %12.5g  (slope)\n",fitparms[3]);
				fprintf(fp_log,"fitparms[4] = %12.5g  (power)\n",fitparms[4]);
				fprintf(fp_log,"*************************************************\n");
				fflush(fp_log);
			}

		}

		if ((optite < 0) || (optite > 2))
		{
			for (i = 0; i < 30; i++)
			{
				if (optite != 3)
				{
					GetMLEParms(beginp, nparm);	/* Get a new starting point */
				}
				else
				{
					for (j = 0; j < nparm; j++)
					{
						beginp[j] = fitparms[j];
					}
				}

				/* Try again */
				gmcnt = gmcnt + 1;
				if (bNo_Log == false)
				{
					fprintf(fp_log,"\n***********Call #%d to getmle.****************\n",gmcnt);
					fprintf(fp_log,"(this call doesn't scale the parameters by their starting values in donlp2)\n");
					fprintf(fp_log,"nparms = %ld\n",nparms);
					fprintf(fp_log,"flag = %ld\n",flag);
					fprintf(fp_log,"These are the parameters going into getmle.  The slope\n");
					fprintf(fp_log,"has been scaled by maxdose^power.\n");
					fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
					fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
					fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
					fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
					fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
					fprintf(fp_log,"************************************************\n");
					fflush(fp_log);
				}

				getmle_(&nvar, doses, means, nanim, svar, &nparms, beginp, Spec2, beginp,
					&restr, &signs, fitparms, &ll, &optite, &nresm, bind,
					&model_type, &flag);

				if (bNo_Log == false)
				{
					fprintf(fp_log,"\n*******After Call #%d to getmle.**********\n",gmcnt);
					fprintf(fp_log,"nparms = %ld\n",nparms);
					fprintf(fp_log,"flag = %ld\n",flag);
					fprintf(fp_log,"ll = %10.5g   (likelihood)\n",ll);
					fprintf(fp_log,"optite = %ld    (good optimum if 0<=optite<=2)\n",optite);
					fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
					fprintf(fp_log,"These are the parameters coming out of the first run\n");
					fprintf(fp_log,"of getmle.  The slope is still in scaled form.\n");
					fprintf(fp_log,"fitparms[0] = %12.5g  (alpha)\n",fitparms[0]);
					fprintf(fp_log,"fitparms[1] = %12.5g  (rho)\n",fitparms[1]);
					fprintf(fp_log,"fitparms[2] = %12.5g  (control)\n",fitparms[2]);
					fprintf(fp_log,"fitparms[3] = %12.5g  (slope)\n",fitparms[3]);
					fprintf(fp_log,"fitparms[4] = %12.5g  (power)\n",fitparms[4]);
					fprintf(fp_log,"*************************************************\n");
					fflush(fp_log);
				}

				if ((optite >= 0) && (optite <= 2))
				{
					flag = 2;
					break;
				}
			}
		}	/* end if ((optite < 0) || (optite > 2)) */

	} /* end  for (flag=0; flag<=1; flag++) */


	/* This decides if the scaling model is better than the unscaled model */
	/* or not. */
	if ((optite2 >= 0) && (optite2 <= 2))
	{
		if ((ll2 > ll) || ((optite < 0) || (optite > 2)))
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
	if ((optite3 >= 0) && (optite3 <= 2))
	{
		if ((ll3 > ll) || ((optite < 0) || (optite > 2)))
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

	/* Unscale by max mean and max variance */
	if (fitparms[1] == 0)
		fitparms[0] = fitparms[0]*pow(maxYm,2);
	else
		fitparms[0] = fitparms[0] + (2 - fitparms[1]) * log(maxYm);
	fitparms[2] = fitparms[2]*maxYm;
	fitparms[3] = fitparms[3]*maxYm;

	mean = DVECTOR(1,Nobs);
	var = DVECTOR(1,Nobs);
	LL = 0;
	/* Find Likelihood with the corrected parmeters */
	for (i = 1; i <= Nobs; i++)
	{
		if (Xi[i] == 0)
			mean[i] = fitparms[2];
		else
			mean[i] = fitparms[2] + fitparms[3]*pow(Xi[i],fitparms[4]);
		if (fitparms[1] == 0)
			var[i] = fitparms[0];
		else
			var[i] = exp(fitparms[0] + log(fabs(mean[i])) * fitparms[1]);
		if (var[i] <= 0)
			var[i] = 0.000000001;
		LL = LL + Ni[i]*Slog(var[i])/2 + (Ni[i]-1)*Yd[i]/(2*var[i]);
		LL = LL + Ni[i]*(pow((Ym[i]-mean[i]),2))/(2*var[i]);
	}
	ll = -LL;

	FREE_DVECTOR(mean, 1, Nobs);
	FREE_DVECTOR(var, 1, Nobs);

	/* Let user know if no optimum was found */
	if ((optite < 0) || (optite > 2))
	{
		fprintf(fp_out, "\n\n!!! Warning:  optimum may not have been found.                      !!!");
		fprintf(fp_out, "\n!!! Bad completion code in maximum likelihood optimization routine  !!!");
		fprintf(fp_out, "\n!!! Program halted                                                  !!!\n");
		fprintf(fp_out, "\n!!! Try choosing different initial values                           !!!\n\n");
		*is_conv = -1;
		exit(0);
	} /* end if */
	else
	{
		*is_conv = 1;
	}

	for (i = 1; i <= nparm; i++)
	{
		p[i] = fitparms[i-1];
		bounded[i] = bind[i-1];
	}

	for (i = 1; i <= Nobs; i++)
	{
		Xi[i] = Xi[i] * maxdose;
	}
	/* rescale the slope */
	p[4] = p[4]/pow(maxdose, p[5]);

	/* Let's try donlp2's idea about these
	if ((restrict == 1) && (fabs(p[5]-1) < 1e-20))
	bounded[5] = 1;
	if ((restrict == 0) && (fabs(p[5]) < 1e-20))
	bounded[5] = 1;
	if ((fabs(p[2]-18) < 1e-20) || (fabs(p[2]+18) < 1e-20))
	bounded[2] = 1;
	*/


	*fret = ll;

	FREE_DVECTOR(pBak, 1, nparm);
	FREE_DVECTOR(doses, 0, Nobs-1);
	FREE_DVECTOR(means, 0, Nobs-1);
	FREE_DVECTOR(svar, 0, Nobs-1);
	FREE_LIVECTOR(nanim, 0, Nobs-1);
	FREE_DVECTOR(parms, 0, nparm-1);
	FREE_DVECTOR(fitparms, 0, nparm-1);
	FREE_DVECTOR(fitparms2, 0, nparm-1);
	FREE_DVECTOR(fitparms3, 0, nparm-1);
	FREE_LIVECTOR(Spec2, 0, nparm-1);
	FREE_LIVECTOR(bind, 0, nparm-1);
	FREE_LIVECTOR(bind2, 0, nparm-1);
	FREE_LIVECTOR(bind3, 0, nparm-1);
	FREE_DVECTOR(beginp, 0, nparm-1);


}	/* end Pow_fit */

/******************************************************************
 * Pow_func -- used to compute the value of the power function
 *              at the point xx, given parm p[] and nparm
 *******************************************************************/
float Pow_func (int nparm, float pBak[], float xx, float gtol)
{
  double fP, *p;
  int j;
  p=DVECTOR(1, nparm);

  fP = p[nparm];
  for (j=nparm-1; j>=1; j--) fP=fP*xx+p[j];

  return fP;

  FREE_DVECTOR(p, 1, nparm);
}

/*****************************************************************
 * BMD_func -- assumes the value 0 if Dose equals BMD.
 *			  This routine is called by Binary_root()
 ******************************************************************/
float BMD_func(int nparm, float pBak[], float D, float gtol)
{
  double fd, fbmd, tD[2], m[2];
  double *p;
  int    j;

  p = DVECTOR(1, nparm);
  for (j=1;j<=nparm;j++)
    p[j]=pBak[j];      /* get the "old" p[]. */

  tD[1] = D;

  PowMeans(1, p, tD, m);

  fd = m[1];

  fbmd = fd - BMR;

  return fbmd;

  FREE_DVECTOR(p, 1, nparm);

}


/**********************************************************
 ** READ_OBSDATA4V--used to read 4 column data in 4 vectors.
 ***********************************************************/
int READ_OBSDATA4V(int Nobs,double Xi[],int Ni[],double Ym[],double Yd[])

{
  int     Nmiss;          /*number of records with missing values*/
  int     i,j,n,m;        /*count and iteration control variables*/
  double  value;          /*temp variable*/

  Nmiss = 0;
  for(i=1;i<=Nobs;i++)
    {
      n = i-Nmiss;
      m = 0;
      for (j=1;j<=4;j++)
	{
	  fscanf(fp_in,"%lf",&value);
	  if (value != MISSING)
	    {
	      if (j==1)     Xi[n]=value;
 	      if (j==2)     Ni[n]=(int)value;
 	      if (j==3)     Ym[n]=value;
 	      if (j==4)     Yd[n]=value;
	    }
	  else
	    m++;
	}
      if (m != 0)	   Nmiss++;
      else if (Xi[n] < 0)      Nmiss++;
    }

  return Nmiss;
}


/**********************************************************
 ** READ_OBSDATA2V--used to read 2 column data in 2 vectors.
 ***********************************************************/
int READ_OBSDATA2V(int ntotal, double xxi[], double yyi[])

{
  int     Nmiss;          /*number of records with missing values*/
  int     i,j,n,m;        /*count and iteration control variables*/
  double  value;          /*temp variable*/

  Nmiss = 0;
  for(i=1;i<=ntotal;i++)
    {
      n = i-Nmiss;
      m = 0;
      for (j=1;j<=2;j++)
	{
	  fscanf(fp_in,"%lf",&value);
	  if (value != MISSING)
	    {
	      if (j==1)     xxi[n]=value;
 	      if (j==2)     yyi[n]=value;
	    }
	  else
	    m++;
	}
      if (m != 0)	   Nmiss++;
      else if (xxi[n] < 0)      Nmiss++;
    }

  return Nmiss;
}

/**********************************************************
 ** READ_OBSDATAV--used to read Ni column data in Ni vectors.
 ***********************************************************/
int READ_OBSDATAV(int Nobs,double Xi[],int Ni[],double Ym[],double Yd[])

{
  int     Nmiss;          /*number of records with missing values*/
  int     i,j,n,m;        /*count and iteration control variables*/
  double  value;          /*temp variable*/
  double  *yy;            /*response of each animal in a group*/

  yy = DVECTOR(1,Nobs+2);/*this line was added to define yy array*/
  Nmiss = 0;
  for(i=1;i<=Nobs;i++)
    {
      n = i-Nmiss;
      m = 0;
      for (j=1;j<=2;j++)
	{
	  fscanf(fp_in,"%lf",&value);
	  if (value != MISSING)
	    {
	      if (j==1)     Xi[n]=value;
	      if (j==2)     Ni[n]=(int)value;
	    }
	  else
	    m++;
	}
      for (j=3; j<=Ni[n]+2; j++)
	{
	  fscanf(fp_in,"%lf",&value);
	  if (value != MISSING)
	    {
	      yy[j-2]=value;
	      Ym[j-2] += value/Ni[n];
	    }
	  else  m++;
	}
      for (j=3; j<=Ni[n]+2; j++)
	{
	  fscanf(fp_in,"%lf",&value);
	  if (value != MISSING)
	    Yd[j-2] += (yy[j-2]-Ym[j-2])*(yy[j-2]-Ym[j-2])/(Ni[n]-1.0);
	  else m++;
	}

      if (m != 0)	  	                  Nmiss++;
      else if (Xi[n] < 0)      Nmiss++;
    }

  return Nmiss;
}

/*******************************************************************
 *	AThree_Fit fits a likelihood to a "full" model of the form
 *	Yij = Mu(i) + e(ij)  Var(eij) = k*Mu(i)^b.  The parameters Mu(i)
 *	i = 1,2,...,Nobs, k, and p are estimated using the likelihood
 *	maximization procedure, and the lkA3 = log-likelihood value
 *	upon return.
 *******************************************************************/

void AThree_Fit(int nparm, double p[], double gtol, int *iter, double *fret)
{

  int i;
  double *parms, *fitparms, *doses, ll, *svar, *means;
  double /**tmy, *t, **tmv,*/ *bsv;
  double **X, **XP, **XPX, *XPY, *Y;
  long int *Spec2, *bind, *nanim, nresm, optite;
  long int restr, nparms, nvar;

  switch (Spec[1] * 10 + Spec[2]) {
  case 11: /* Both alpha and rho fixed */
    break;
  case 1: /* Estimating alpha; rho fixed */
    /* estimate alpha as mean of Y/pow(mean, rho) */
    p[1] = 0.0;
    for (i = 1; i <= Nobs; i++)
      p[1] += Yd[i] / pow(Ym[i], p[2]);
    p[1] = log(p[1] / (double) Nobs);
    break;
  case 10: /* Estimate rho, alpha fixed */
    /* estimate rho as mean of (log(Var) - log(alpha)) / log(mu) */
    p[2] = 0.0;
    for (i = 1; i <= Nobs; i++)
      p[2] += (Slog(Yd[i]) - Slog(p[1])) / Slog(Ym[i]);
    p[2] = p[2] / (double) Nobs;
    break;
  default:  /* estimating both variance parameters */
    /* Allocate memory for linear regression */
    X = DMATRIX(1, Nobs, 1, 2);
    XP = DMATRIX(1, 2, 1, Nobs);
    XPX = DMATRIX(1, 2, 1, 2);
    XPY = DVECTOR(1, 2);
    bsv = DVECTOR(1, 2);
    Y = DVECTOR(1, Nobs);
    
    /* set up regression problem; regress log(variance)
       on log means */
    for (i = 1; i <= Nobs; i++)
      {
	X[i][1] = 1.0;
	X[i][2] = Slog(Ym[i]);
	Y[i] = Slog(Yd[i]);
      }
    /* linear regression */
    TRANSPOSE(X, XP, Nobs, 2);
    MATMPYM2(XP, X, XPX, 2, Nobs, 2);
    INVMAT(XPX, 2);
    MATMPYV2(2, Nobs, XP, Y, XPY);
    MATMPYV2(2, 2, XPX, XPY, bsv);
    /* transfer the results to p[1:2] */
    p[1] = bsv[1];
    p[2] = bsv[2];
    /* Free memory */
    FREE_DMATRIX(X, 1, Nobs, 1, 2);
    FREE_DMATRIX(XP, 1, 2, 1, Nobs);
    FREE_DMATRIX(XPX, 1, 2, 1, 2);
    FREE_DVECTOR(XPY, 1, 2);
    FREE_DVECTOR(bsv, 1, 2);
    FREE_DVECTOR(Y, 1, Nobs);
    break;
    
    
  }

  for (i = 1; i <= nparm-2; i++)
    {
      p[i + 2] = Ym[i];
    }

  nvar = Nobs;
  restr = restrict;
  nparms = nparm;

  doses = DVECTOR(0, Nobs-1);
  means = DVECTOR(0, Nobs-1);
  svar =  DVECTOR(0, Nobs-1);
  nanim = LIVECTOR(0, Nobs-1);
  parms = DVECTOR(0, nparm-1);
  fitparms = DVECTOR(0, nparm-1);
  Spec2 = LIVECTOR(0, nparm-1);
  bind = LIVECTOR(0, nparm-1);

  for (i = 1; i <= Nobs; i++)
    {
      nanim[i - 1] = (long int) Ni[i];
      doses[i - 1] = Xi[i];
      means[i - 1] = Ym[i];
      svar[i - 1] = Yd[i];
    }

  for (i = 1; i <= nparm; i++)
    {
      Spec2[i - 1] = Spec[i];
      parms[i - 1] = p[i];
    }

  getmlea3_(&nvar, doses, means, nanim, svar, &nparms, parms,
	    Spec2, parms, &restr, fitparms, &ll, &optite,
	    &nresm, bind);

  for (i = 1; i <= nparm; i++)
    {
      p[i] = fitparms[i - 1];
    }

  *fret = ll;


  FREE_DVECTOR(doses, 0, Nobs-1);
  FREE_DVECTOR(means, 0, Nobs-1);
  FREE_DVECTOR(svar, 0, Nobs-1);
  FREE_LIVECTOR(nanim, 0, Nobs-1);
  FREE_DVECTOR(parms, 0, nparm-1);
  FREE_DVECTOR(fitparms, 0, nparm-1);
  FREE_LIVECTOR(Spec2, 0, nparm-1);
  FREE_LIVECTOR(bind, 0, nparm-1);


}	/* end AThree_Fit */

/***********************************************************
 *	Given a vector of parameter values, and the number of
 *	parameters in that vector, this function will return three
 *	new parameter values to restart the optimization if a "bad"
 *	completion code is returned from GETCL(), using a uniform
 *	random number centered at p[i]
 ***********************************************************/
void GetNewParms(double *p, int size)
{
  int i;


  /* Find parameters by randomly selecting new parameters in
     /  a uniform interval of p[i] +/- .005*p[i] */

  for (i = 0; i < size; i++)
    {
      if (Spec[i + 1] != 1)
	{
	  p[i] = ((p[i] * .010) * (double)rand() / 32767.0) + p[i] * .995;
	}
    }

  /* If parameters are to be restricted, make sure restrictions
     /  are not violated */

  if (p[0] <= 0)
    {
      p[0] = -p[0];
    }

  if ((restrict == 1) && (p[4] <= 1))
    {
      p[4] = 1.00000001;
    }
  if ((restrict == 0) && (p[4] < 0))
    {
      p[4] = -p[4];
    }

}	/* end GetNewParms */

/*****************************************************************
 *	BMDL_func -- used to compute the values of functions BMDL_f (the
 *		     X^2 value) at the point D, given the parm p[] and
 *		     number of parm. If GETCL finishes with bad
 *		     convergence, then the return value is -1
 *****************************************************************/
double BMDL_func(int nparm, double xlk, double Dose, double pBak[],
		 double gtol)
{	/* BMD_lk and LR are calculated in Pow_BMD() */

  double fD, bmdl, *doses, target, *parms, *means, *beginp;
  double *svar, *fitparms, maxdose;
  int i, j, ii, gccnt;
  long int *Spec2, flag, which, *nanim, *bind, nresm, optite;
  long int model_type, nvar, nparms, signs, restr, bmrtype;

  /* Power model */

  model_type = 1;
  nvar = Nobs;
  nparms = nparm;
  signs = sign;
  restr = restrict;
  bmrtype = bmr_type;

  doses = DVECTOR(0, Nobs-1);
  means = DVECTOR(0, Nobs-1);
  svar =  DVECTOR(0, Nobs-1);
  nanim = LIVECTOR(0, Nobs-1);
  parms = DVECTOR(0, nparm-1);
  fitparms = DVECTOR(0, nparm-1);
  Spec2 = LIVECTOR(0, nparm-1);
  bind = LIVECTOR(0, nparm-1);
  beginp = DVECTOR(0, nparm-1);

  for (i = 1; i <= Nobs; i++)
    {
      nanim[i - 1] = (long int) Ni[i];
      doses[i - 1] = Xi[i];
      means[i - 1] = Ym[i];
      svar[i - 1] = Yd[i];
    }
  which = 1;			/* Want an lower confidence limit */
  target = xlk - LR;	/* The value we want the likelihood
			   /  at the BMDL to match */
  flag = 0;
  for (j = 1; j <= nparm; j++)
    {
      parms[j - 1] = pBak[j];		/* get the "old" p[] */
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
      doses[i] = doses[i] / maxdose;
    }
  Dose = Dose / maxdose;

  /* slope parameter effected by scaling */
  parms[3] = parms[3]*pow(maxdose, parms[4]);
  beginp[3] = beginp[3]*pow(maxdose, parms[4]);

  /* This is the first set of calls to getcl.  These will all have    */
  /* donlp2 scale the parameters internally.  The max number of calls */
  /* to getcl for this set is 17.                                     */

  gccnt = 1;
  if (bNo_Log == false)
    {
      fprintf(fp_log,"\n***********Call #%d to getcl.****************\n",gccnt);
      fprintf(fp_log,"(this call scales the parameters by their starting values in donlp2)\n");
      fprintf(fp_log,"nparms = %ld\n",nparms);
      fprintf(fp_log,"flag = %ld\n",flag);
      fprintf(fp_log,"BMR = %10.5g\n",BMR);
      fprintf(fp_log,"bmrtype = %ld\n",bmrtype);
      fprintf(fp_log,"target = %10.5g\n",target);
      fprintf(fp_log,"bmdl = %10.5g  (This should be the BMD)\n",bmdl);
      fprintf(fp_log,"These are the parameters going into getcl.  The slope\n");
      fprintf(fp_log,"has been scaled by maxdose^power.\n");
      fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
      fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
      fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
      fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
      fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
      fprintf(fp_log,"************************************************\n");
      fflush(fp_log);
    }

  getcl_(&which, &nvar, doses, means, nanim, svar, &nparms, &BMR,
	 &Dose, &target, parms, Spec2, parms, &bmrtype, &restr,
	 &bmdl, fitparms, &optite, &nresm, bind, &signs, &model_type, &flag);

  if (bNo_Log == false)
    {
      fprintf(fp_log,"\n*******After Call #%d to getcl.**********\n",gccnt);
      fprintf(fp_log,"nparms = %ld\n",nparms);
      fprintf(fp_log,"flag = %ld\n",flag);
      fprintf(fp_log,"bmdl = %10.5g   (BMDL)\n",bmdl);
      fprintf(fp_log,"optite = %ld    (good optimum if 0<=optite<=2)\n",optite);
      fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
      fprintf(fp_log,"These are the parameters coming out of the first run\n");
      fprintf(fp_log,"of getcl.  The slope is still in scaled form.\n");
      fprintf(fp_log,"fitparms[0] = %12.5g  (alpha)\n",fitparms[0]);
      fprintf(fp_log,"fitparms[1] = %12.5g  (rho)\n",fitparms[1]);
      fprintf(fp_log,"fitparms[2] = %12.5g  (control)\n",fitparms[2]);
      fprintf(fp_log,"fitparms[3] = %12.5g  (slope)\n",fitparms[3]);
      fprintf(fp_log,"fitparms[4] = %12.5g  (power)\n",fitparms[4]);
      fprintf(fp_log,"*************************************************\n");
      fflush(fp_log);
    }

  /* optite is a value that is passed back from GETCL which
     /  determines whether the optimization was completed successfully
     /  If optite is less than 0, then it did not, and we want
     /  to try a different starting point and recompute */

  if ((optite < 0) || (optite > 3))
    {

      for (ii = 0; ii < 10; ii++)
	{
#ifdef MISC_OUT
	  printf("\n optite = %ld", optite);
#endif
	  if (optite != 4)
	    {
	      GetNewParms(beginp, nparm);		/* Get a new starting point */
	    } /* end if */
	  else
	    {
	      for (i = 0; i < nparm; i++)
		{
		  beginp[i] = fitparms[i];
		}
	    } /* end else */

	  /* Try again */
	  gccnt = gccnt + 1;
	  if (bNo_Log == false)
	    {
	      fprintf(fp_log,"\n***********Call #%d to getcl.****************\n",gccnt);
	      fprintf(fp_log,"(this call scales the parameters by their starting values in donlp2)\n");
	      fprintf(fp_log,"nparms = %ld\n",nparms);
	      fprintf(fp_log,"flag = %ld\n",flag);
	      fprintf(fp_log,"BMR = %10.5g\n",BMR);
	      fprintf(fp_log,"bmrtype = %ld\n",bmrtype);
	      fprintf(fp_log,"target = %10.5g\n",target);
	      fprintf(fp_log,"bmdl = %10.5g  (This should be the BMD)\n",bmdl);
	      fprintf(fp_log,"These are the parameters going into getcl.  The slope\n");
	      fprintf(fp_log,"has been scaled by maxdose^power.\n");
	      fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
	      fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
	      fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
	      fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
	      fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
	      fprintf(fp_log,"************************************************\n");
	      fflush(fp_log);
	    }

	  getcl_(&which, &nvar, doses, means, nanim, svar, &nparms, &BMR,
		 &Dose, &target, beginp, Spec2, beginp, &bmrtype, &restr,
		 &bmdl, fitparms, &optite, &nresm, bind, &signs, &model_type, &flag);

	  if (bNo_Log == false)
	    {
	      fprintf(fp_log,"\n*******After Call #%d to getcl.**********\n",gccnt);
	      fprintf(fp_log,"nparms = %ld\n",nparms);
	      fprintf(fp_log,"flag = %ld\n",flag);
	      fprintf(fp_log,"bmdl = %10.5g   (BMDL)\n",bmdl);
	      fprintf(fp_log,"optite = %ld    (good optimum if 0<=optite<=2)\n",optite);
	      fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
	      fprintf(fp_log,"These are the parameters coming out of the first run\n");
	      fprintf(fp_log,"of getcl.  The slope is still in scaled form.\n");
	      fprintf(fp_log,"fitparms[0] = %12.5g  (alpha)\n",fitparms[0]);
	      fprintf(fp_log,"fitparms[1] = %12.5g  (rho)\n",fitparms[1]);
	      fprintf(fp_log,"fitparms[2] = %12.5g  (control)\n",fitparms[2]);
	      fprintf(fp_log,"fitparms[3] = %12.5g  (slope)\n",fitparms[3]);
	      fprintf(fp_log,"fitparms[4] = %12.5g  (power)\n",fitparms[4]);
	      fprintf(fp_log,"*************************************************\n");
	      fflush(fp_log);
	    }

	  /* if optite >= 0 and < 3, it is successful, and we can stop */

	  if ((optite >= 0) && (optite <= 3))
	    {
	      break;
	    }
	}  /* end for (ii = 0; ii < 15; ii++) */
    }  /* end if ((optite < 0) || (optite >= 3)) */


  if ((optite < 0) || (optite > 3))
    {
      for (j = 1; j <= nparm; j++)
	beginp[j - 1] = pBak[j];
      beginp[3] = beginp[3]*pow(maxdose, parms[4]);

      gccnt = gccnt + 1;
      if (bNo_Log == false)
	{
	  fprintf(fp_log,"\n***********Call #%d to getcl.****************\n",gccnt);
	  fprintf(fp_log,"(this call scales the parameters by their starting values in donlp2)\n");
	  fprintf(fp_log,"nparms = %ld\n",nparms);
	  fprintf(fp_log,"flag = %ld\n",flag);
	  fprintf(fp_log,"BMR = %10.5g\n",BMR);
	  fprintf(fp_log,"bmrtype = %ld\n",bmrtype);
	  fprintf(fp_log,"target = %10.5g\n",target);
	  fprintf(fp_log,"bmdl = %10.5g  (This should be the BMD)\n",bmdl);
	  fprintf(fp_log,"These are the parameters going into getcl.  The slope\n");
	  fprintf(fp_log,"has been scaled by maxdose^power.\n");
	  fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
	  fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
	  fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
	  fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
	  fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
	  fprintf(fp_log,"************************************************\n");
	  fflush(fp_log);
	}

      getcl_(&which, &nvar, doses, means, nanim, svar, &nparms, &BMR,
	     &Dose, &target, beginp, Spec2, beginp, &bmrtype, &restr,
	     &bmdl, fitparms, &optite, &nresm, bind, &signs, &model_type, &flag);

      if (bNo_Log == false)
	{
	  fprintf(fp_log,"\n*******After Call #%d to getcl.**********\n",gccnt);
	  fprintf(fp_log,"nparms = %ld\n",nparms);
	  fprintf(fp_log,"flag = %ld\n",flag);
	  fprintf(fp_log,"bmdl = %10.5g   (BMDL)\n",bmdl);
	  fprintf(fp_log,"optite = %ld    (good optimum if 0<=optite<=2)\n",optite);
	  fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
	  fprintf(fp_log,"These are the parameters coming out of the first run\n");
	  fprintf(fp_log,"of getcl.  The slope is still in scaled form.\n");
	  fprintf(fp_log,"fitparms[0] = %12.5g  (alpha)\n",fitparms[0]);
	  fprintf(fp_log,"fitparms[1] = %12.5g  (rho)\n",fitparms[1]);
	  fprintf(fp_log,"fitparms[2] = %12.5g  (control)\n",fitparms[2]);
	  fprintf(fp_log,"fitparms[3] = %12.5g  (slope)\n",fitparms[3]);
	  fprintf(fp_log,"fitparms[4] = %12.5g  (power)\n",fitparms[4]);
	  fprintf(fp_log,"*************************************************\n");
	  fflush(fp_log);
	}
    }

  if ((optite < 0) || (optite > 3))
    {
      for (ii = 0; ii < 5; ii++)
	{
#ifdef MISC_OUT
	  printf("\n     optite = %ld", optite);
#endif
	  if (optite != 4)
	    {
	      GetOtherParms(beginp, nparm);		/* Get a new starting point */
	    } /* end if */
	  else
	    {
	      for (i = 0; i < nparm; i++)
		{
		  beginp[i] = fitparms[i];
		}
	    } /* end else */

	  /* Try again */
	  gccnt = gccnt + 1;
	  if (bNo_Log == false)
	    {
	      fprintf(fp_log,"\n***********Call #%d to getcl.****************\n",gccnt);
	      fprintf(fp_log,"(this call scales the parameters by their starting values in donlp2)\n");
	      fprintf(fp_log,"nparms = %ld\n",nparms);
	      fprintf(fp_log,"flag = %ld\n",flag);
	      fprintf(fp_log,"BMR = %10.5g\n",BMR);
	      fprintf(fp_log,"bmrtype = %ld\n",bmrtype);
	      fprintf(fp_log,"target = %10.5g\n",target);
	      fprintf(fp_log,"bmdl = %10.5g  (This should be the BMD)\n",bmdl);
	      fprintf(fp_log,"These are the parameters going into getcl.  The slope\n");
	      fprintf(fp_log,"has been scaled by maxdose^power.\n");
	      fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
	      fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
	      fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
	      fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
	      fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
	      fprintf(fp_log,"************************************************\n");
	      fflush(fp_log);
	    }

	  getcl_(&which, &nvar, doses, means, nanim, svar, &nparms, &BMR,
		 &Dose, &target, beginp, Spec2, beginp, &bmrtype, &restr,
		 &bmdl, fitparms, &optite, &nresm, bind, &signs, &model_type, &flag);

	  if (bNo_Log == false)
	    {
	      fprintf(fp_log,"\n*******After Call #%d to getcl.**********\n",gccnt);
	      fprintf(fp_log,"nparms = %ld\n",nparms);
	      fprintf(fp_log,"flag = %ld\n",flag);
	      fprintf(fp_log,"bmdl = %10.5g   (BMDL)\n",bmdl);
	      fprintf(fp_log,"optite = %ld    (good optimum if 0<=optite<=2)\n",optite);
	      fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
	      fprintf(fp_log,"These are the parameters coming out of the first run\n");
	      fprintf(fp_log,"of getcl.  The slope is still in scaled form.\n");
	      fprintf(fp_log,"fitparms[0] = %12.5g  (alpha)\n",fitparms[0]);
	      fprintf(fp_log,"fitparms[1] = %12.5g  (rho)\n",fitparms[1]);
	      fprintf(fp_log,"fitparms[2] = %12.5g  (control)\n",fitparms[2]);
	      fprintf(fp_log,"fitparms[3] = %12.5g  (slope)\n",fitparms[3]);
	      fprintf(fp_log,"fitparms[4] = %12.5g  (power)\n",fitparms[4]);
	      fprintf(fp_log,"*************************************************\n");
	      fflush(fp_log);
	    }
	  /* if optite >= 0 and < 3, it is successful, and we can stop */

	  if ((optite >= 0) && (optite <= 3))
	    {
	      break;
	    }
	}	/* end for (ii = 0; ii < 15; ii++) */
    }	/* end if ((optite < 0) || (optite >= 3)) */
  /*END OF FIRST SET OF CALLS TO getcl. */


  /* This is the second set of calls to getcl.  These will not have    */
  /* donlp2 scale the parameters internally.  The max number of calls */
  /* to getcl for this set is 17.                                     */

  if ((optite < 0) || (optite > 3))
    {
      gccnt = gccnt + 1;
      if (bNo_Log == false)
	{
	  fprintf(fp_log,"\n***********Call #%d to getcl.****************\n",gccnt);
	  fprintf(fp_log,"(this call doesn't scale the parameters by their starting values in donlp2)\n");
	  fprintf(fp_log,"nparms = %ld\n",nparms);
	  fprintf(fp_log,"flag = %ld\n",flag);
	  fprintf(fp_log,"BMR = %10.5g\n",BMR);
	  fprintf(fp_log,"bmrtype = %ld\n",bmrtype);
	  fprintf(fp_log,"target = %10.5g\n",target);
	  fprintf(fp_log,"bmdl = %10.5g  (This should be the BMD)\n",bmdl);
	  fprintf(fp_log,"These are the parameters going into getcl.  The slope\n");
	  fprintf(fp_log,"has been scaled by maxdose^power.\n");
	  fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
	  fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
	  fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
	  fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
	  fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
	  fprintf(fp_log,"************************************************\n");
	  fflush(fp_log);
	}

      getcl_(&which, &nvar, doses, means, nanim, svar, &nparms, &BMR,
	     &Dose, &target, parms, Spec2, parms, &bmrtype, &restr,
	     &bmdl, fitparms, &optite, &nresm, bind, &signs, &model_type, &flag);

      if (bNo_Log == false)
	{
	  fprintf(fp_log,"\n*******After Call #%d to getcl.**********\n",gccnt);
	  fprintf(fp_log,"nparms = %ld\n",nparms);
	  fprintf(fp_log,"flag = %ld\n",flag);
	  fprintf(fp_log,"bmdl = %10.5g   (BMDL)\n",bmdl);
	  fprintf(fp_log,"optite = %ld    (good optimum if 0<=optite<=2)\n",optite);
	  fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
	  fprintf(fp_log,"These are the parameters coming out of the first run\n");
	  fprintf(fp_log,"of getcl.  The slope is still in scaled form.\n");
	  fprintf(fp_log,"fitparms[0] = %12.5g  (alpha)\n",fitparms[0]);
	  fprintf(fp_log,"fitparms[1] = %12.5g  (rho)\n",fitparms[1]);
	  fprintf(fp_log,"fitparms[2] = %12.5g  (control)\n",fitparms[2]);
	  fprintf(fp_log,"fitparms[3] = %12.5g  (slope)\n",fitparms[3]);
	  fprintf(fp_log,"fitparms[4] = %12.5g  (power)\n",fitparms[4]);
	  fprintf(fp_log,"*************************************************\n");
	  fflush(fp_log);
	}

      /* optite is a value that is passed back from GETCL which
	 /  determines whether the optimization was completed successfully
	 /  If optite is less than 0, then it did not, and we want
	 /  to try a different starting point and recompute */

      if ((optite < 0) || (optite > 3))
	{
	  for (ii = 0; ii < 10; ii++)
	    {
#ifdef MISC_OUT
	      printf("\n optite = %ld", optite);
#endif
	      if (optite != 4)
		{
		  GetNewParms(beginp, nparm);		/* Get a new starting point */
		} /* end if */
	      else
		{
		  for (i = 0; i < nparm; i++)
		    {
		      beginp[i] = fitparms[i];
		    }
		} /* end else */

	      /* Try again */
	      gccnt = gccnt + 1;
	      if (bNo_Log == false)
		{
		  fprintf(fp_log,"\n***********Call #%d to getcl.****************\n",gccnt);
		  fprintf(fp_log,"(this call doesn't scale the parameters by their starting values in donlp2)\n");
		  fprintf(fp_log,"nparms = %ld\n",nparms);
		  fprintf(fp_log,"flag = %ld\n",flag);
		  fprintf(fp_log,"BMR = %10.5g\n",BMR);
		  fprintf(fp_log,"bmrtype = %ld\n",bmrtype);
		  fprintf(fp_log,"target = %10.5g\n",target);
		  fprintf(fp_log,"bmdl = %10.5g  (This should be the BMD)\n",bmdl);
		  fprintf(fp_log,"These are the parameters going into getcl.  The slope\n");
		  fprintf(fp_log,"has been scaled by maxdose^power.\n");
		  fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
		  fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
		  fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
		  fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
		  fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
		  fprintf(fp_log,"************************************************\n");
		  fflush(fp_log);
		}

	      getcl_(&which, &nvar, doses, means, nanim, svar, &nparms, &BMR,
		     &Dose, &target, beginp, Spec2, beginp, &bmrtype, &restr,
		     &bmdl, fitparms, &optite, &nresm, bind, &signs, &model_type, &flag);

	      if (bNo_Log == false)
		{
		  fprintf(fp_log,"\n*******After Call #%d to getcl.**********\n",gccnt);
		  fprintf(fp_log,"nparms = %ld\n",nparms);
		  fprintf(fp_log,"flag = %ld\n",flag);
		  fprintf(fp_log,"bmdl = %10.5g   (BMDL)\n",bmdl);
		  fprintf(fp_log,"optite = %ld    (good optimum if 0<=optite<=2)\n",optite);
		  fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
		  fprintf(fp_log,"These are the parameters coming out of the first run\n");
		  fprintf(fp_log,"of getcl.  The slope is still in scaled form.\n");
		  fprintf(fp_log,"fitparms[0] = %12.5g  (alpha)\n",fitparms[0]);
		  fprintf(fp_log,"fitparms[1] = %12.5g  (rho)\n",fitparms[1]);
		  fprintf(fp_log,"fitparms[2] = %12.5g  (control)\n",fitparms[2]);
		  fprintf(fp_log,"fitparms[3] = %12.5g  (slope)\n",fitparms[3]);
		  fprintf(fp_log,"fitparms[4] = %12.5g  (power)\n",fitparms[4]);
		  fprintf(fp_log,"*************************************************\n");
		  fflush(fp_log);
		}
	      /* if optite >= 0 and < 3, it is successful, and we can stop */

	      if ((optite >= 0) && (optite <= 3))
		{
		  break;
		}
	    }	/* end for (ii = 0; ii < 15; ii++) */
	}	/* end if ((optite < 0) || (optite >= 3)) */

      if ((optite < 0) || (optite > 3))
	{
	  for (j = 1; j <= nparm; j++)
	    beginp[j - 1] = pBak[j];
	  beginp[3] = beginp[3]*pow(maxdose, parms[4]);

	  gccnt = gccnt + 1;
	  if (bNo_Log == false)
	    {
	      fprintf(fp_log,"\n***********Call #%d to getcl.****************\n",gccnt);
	      fprintf(fp_log,"(this call doesn't scale the parameters by their starting values in donlp2)\n");
	      fprintf(fp_log,"nparms = %ld\n",nparms);
	      fprintf(fp_log,"flag = %ld\n",flag);
	      fprintf(fp_log,"BMR = %10.5g\n",BMR);
	      fprintf(fp_log,"bmrtype = %ld\n",bmrtype);
	      fprintf(fp_log,"target = %10.5g\n",target);
	      fprintf(fp_log,"bmdl = %10.5g  (This should be the BMD)\n",bmdl);
	      fprintf(fp_log,"These are the parameters going into getcl.  The slope\n");
	      fprintf(fp_log,"has been scaled by maxdose^power.\n");
	      fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
	      fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
	      fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
	      fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
	      fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
	      fprintf(fp_log,"************************************************\n");
	      fflush(fp_log);
	    }

	  getcl_(&which, &nvar, doses, means, nanim, svar, &nparms, &BMR,
		 &Dose, &target, beginp, Spec2, beginp, &bmrtype, &restr,
		 &bmdl, fitparms, &optite, &nresm, bind, &signs, &model_type, &flag);

	  if (bNo_Log == false)
	    {
	      fprintf(fp_log,"\n*******After Call #%d to getcl.**********\n",gccnt);
	      fprintf(fp_log,"nparms = %ld\n",nparms);
	      fprintf(fp_log,"flag = %ld\n",flag);
	      fprintf(fp_log,"bmdl = %10.5g   (BMDL)\n",bmdl);
	      fprintf(fp_log,"optite = %ld    (good optimum if 0<=optite<=2)\n",optite);
	      fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
	      fprintf(fp_log,"These are the parameters coming out of the first run\n");
	      fprintf(fp_log,"of getcl.  The slope is still in scaled form.\n");
	      fprintf(fp_log,"fitparms[0] = %12.5g  (alpha)\n",fitparms[0]);
	      fprintf(fp_log,"fitparms[1] = %12.5g  (rho)\n",fitparms[1]);
	      fprintf(fp_log,"fitparms[2] = %12.5g  (control)\n",fitparms[2]);
	      fprintf(fp_log,"fitparms[3] = %12.5g  (slope)\n",fitparms[3]);
	      fprintf(fp_log,"fitparms[4] = %12.5g  (power)\n",fitparms[4]);
	      fprintf(fp_log,"*************************************************\n");
	      fflush(fp_log);
	    }
	}

      if ((optite < 0) || (optite > 3))
	{
	  for (ii = 0; ii < 5; ii++)
	    {
#ifdef MISC_OUT
	      printf("\n     optite = %ld", optite);
#endif
	      if (optite != 4)
		{
		  GetOtherParms(beginp, nparm);		/* Get a new starting point */
		} /* end if */
	      else
		{
		  for (i = 0; i < nparm; i++)
		    {
		      beginp[i] = fitparms[i];
		    }
		} /* end else */

	      /* Try again */
	      gccnt = gccnt + 1;
	      if (bNo_Log == false)
		{
		  fprintf(fp_log,"\n***********Call #%d to getcl.****************\n",gccnt);
		  fprintf(fp_log,"(this call doesn't scale the parameters by their starting values in donlp2)\n");
		  fprintf(fp_log,"nparms = %ld\n",nparms);
		  fprintf(fp_log,"flag = %ld\n",flag);
		  fprintf(fp_log,"BMR = %10.5g\n",BMR);
		  fprintf(fp_log,"bmrtype = %ld\n",bmrtype);
		  fprintf(fp_log,"target = %10.5g\n",target);
		  fprintf(fp_log,"bmdl = %10.5g  (This should be the BMD)\n",bmdl);
		  fprintf(fp_log,"These are the parameters going into getcl.  The slope\n");
		  fprintf(fp_log,"has been scaled by maxdose^power.\n");
		  fprintf(fp_log,"beginp[0] = %12.5g  (alpha)\n",beginp[0]);
		  fprintf(fp_log,"beginp[1] = %12.5g  (rho)\n",beginp[1]);
		  fprintf(fp_log,"beginp[2] = %12.5g  (control)\n",beginp[2]);
		  fprintf(fp_log,"beginp[3] = %12.5g  (slope)\n",beginp[3]);
		  fprintf(fp_log,"beginp[4] = %12.5g  (power)\n",beginp[4]);
		  fprintf(fp_log,"************************************************\n");
		  fflush(fp_log);
		}

	      getcl_(&which, &nvar, doses, means, nanim, svar, &nparms, &BMR,
		     &Dose, &target, beginp, Spec2, beginp, &bmrtype, &restr,
		     &bmdl, fitparms, &optite, &nresm, bind, &signs, &model_type, &flag);

	      if (bNo_Log == false)
		{
		  fprintf(fp_log,"\n*******After Call #%d to getcl.**********\n",gccnt);
		  fprintf(fp_log,"nparms = %ld\n",nparms);
		  fprintf(fp_log,"flag = %ld\n",flag);
		  fprintf(fp_log,"bmdl = %10.5g   (BMDL)\n",bmdl);
		  fprintf(fp_log,"optite = %ld    (good optimum if 0<=optite<=2)\n",optite);
		  fprintf(fp_log,"nresm = %ld    (no idea what this does)\n",nresm);
		  fprintf(fp_log,"These are the parameters coming out of the first run\n");
		  fprintf(fp_log,"of getcl.  The slope is still in scaled form.\n");
		  fprintf(fp_log,"fitparms[0] = %12.5g  (alpha)\n",fitparms[0]);
		  fprintf(fp_log,"fitparms[1] = %12.5g  (rho)\n",fitparms[1]);
		  fprintf(fp_log,"fitparms[2] = %12.5g  (control)\n",fitparms[2]);
		  fprintf(fp_log,"fitparms[3] = %12.5g  (slope)\n",fitparms[3]);
		  fprintf(fp_log,"fitparms[4] = %12.5g  (power)\n",fitparms[4]);
		  fprintf(fp_log,"*************************************************\n");
		  fflush(fp_log);
		}
	      /* if optite >= 0 and < 3, it is successful, and we can stop */

	      if ((optite >= 0) && (optite <= 3))
		{
		  break;
		}
	    }	/* end for (ii = 0; ii < 15; ii++) */
	}   /* end if ((optite < 0) || (optite >= 3)) */
    } /* END OF SECOND SET OF CALLS TO getcl. */


  /* Let user know if no optimum was found */
  if ((optite < 0) || (optite > 3))
    {
      fprintf(fp_out, "Warning:  optimum may not have been found.  Bad completion code in Optimization routine.\n");
      fD = -1;
    } else {   /* end if ((optite < 0) || (optite >= 3)) */
    for (j = 1; j <= nparm; j++)
      {
	pBak[j] = fitparms[j-1];
      }

    /* rescale dose and slope back to normal */
    pBak[4] = pBak[4]/pow(maxdose,pBak[5]);
    Dose = Dose * maxdose;
    
    bmdl = bmdl * maxdose;
    
#ifdef MISC_OUT
    printf("\n optite = %ld", optite);
#endif
    fD = bmdl;
  }

  FREE_DVECTOR(doses, 0, Nobs-1);
  FREE_DVECTOR(means, 0, Nobs-1);
  FREE_DVECTOR(svar, 0, Nobs-1);
  FREE_LIVECTOR(nanim, 0, Nobs-1);
  FREE_DVECTOR(parms, 0, nparm-1);
  FREE_DVECTOR(fitparms, 0, nparm-1);
  FREE_LIVECTOR(Spec2, 0, nparm-1);
  FREE_LIVECTOR(bind, 0, nparm-1);
  FREE_DVECTOR(beginp, 0, nparm-1);

  return fD;			/* return bmdl to the calling function */

}	/* end BMDL_func */

/************************************************************
 * Pow_BMD -- Used to calculate the BMD and BMDL for Power model.
 * nparm, BMD_lk, xb, p, tol
 * If the calculation fails, Pow_BMD prints a message and returns.
 *************************************************************/
void Pow_BMD (int nparm, double p[], double gtol, int *iter, double xlk,
	      double Rlevel[], double Bmdl[], double *BMD)
{
  float Binary_root(float (*O_func)(int, float [], float, float),
		    float xa, float xb, float fa, float fb, float tol,
		    int nparm, float p[]);
  double BMD_func(int nparm, float p[], float x, float gtol);

  double   tol;
  double   xa,xb,fa,fb, bDose;
  double   D, bmrtemp, temp;
  double   *pBak, *BMRVals;
  double   incre, divide, thedose, tempdose[2], bmdlmean[2];
  int      j, k, bmrfncsign;


  pBak=DVECTOR(1, nparm);
  BMRVals=DVECTOR(1,5);

  for(j=1; j<=nparm; j++)
    pBak[j]= p[j];          /* save the p[]. */

  /**** compute X^2 value  ************************/
  /* IF ML is the value of the maximized log-likelihood, then ML - LR is the value
     log-likelihood at the BMDL or BMDU */
  if (bmdparm.level<0.5)
    LR = 0.5*QCHISQ(1.0 - 2.0 * bmdparm.level, 1);
  else
    LR = 0.5*QCHISQ(2.0 * bmdparm.level - 1.0, 1);

  if(bmr_type != 3 && bmr_type != 4)
    Rlevel[1] = BMR = fabs(bmdparm.effect);
  else
    Rlevel[1] = BMR = bmdparm.effect;

  xa = D = 0.0;
  tol=0.0000001;
  incre=1.0;
  divide=50.0;
  bmrfncsign = 1;
  fa = -1;
  temp = 0;

  if (p[5] == 0)
    {
      fprintf(fp_out, "\n Since the power was estimated to be 0,");
      fprintf(fp_out, " the BMD is infinite.");
      fprintf(fp_out, "\n Setting BMD = 100*(maximum dose). ");
      thedose = 100*xmax;
    } /* end if */
  else if (p[4] == 0)
    {
      fprintf(fp_out, "\n Since the slope was estimated to be 0,");
      fprintf(fp_out, " the BMD is infinite.");
      fprintf(fp_out, "\n Setting BMD = 100*(maximum dose).");
      thedose = 100*xmax;
    } /* end if */
  else
    {

      /* Get a BMR value that is equal to just Mu(BMD) regardless
	 of the BMR type to simplify the root finding routine */

      if(bmr_type == 0 || bmr_type == 1 || bmr_type == 2)
	{
	  if(sign == -1)
	    bmrtemp = -BMR;
	  else
	    bmrtemp = BMR;

	  if(bmr_type == 1)
	    {
	      if(p[2] != 0)
		bmrtemp = bmrtemp*sqrt(exp(p[1] + log(p[3]) * p[2]));
	      else
		bmrtemp = bmrtemp*sqrt(p[1]);
	    }
	  else if(bmr_type == 2)
	    {
	      bmrtemp = bmrtemp*p[3];

	    }
	}
      else
	{
	  if(bmr_type == 3)
	    {
	      bmrtemp = BMR - p[3];
	    }
	  else
	    bmrtemp = p[4]*BMR;
	}

      bmrtemp = bmrtemp + p[3];
      /* bmrtype is now a Point BMR type
	 regardless of the user input type */

      temp = BMR;
      BMR = bmrtemp;
      BMRVals[1] = BMR;

      if ((p[4] > 0) && (BMR < p[3]))
	{
	  fprintf (fp_out, "\n\nERROR: BMR value is outside the range of the mean function.\n");
	  FREE_DVECTOR(pBak, 1, nparm);
	  FREE_DVECTOR(BMRVals, 1, 5);
	  return;
	}
      if ((p[4] < 0) && (BMR > p[3]))
	{
	  fprintf (fp_out, "\n\nERROR: BMR value is outside the range of the mean function.\n");
	  FREE_DVECTOR(pBak, 1, nparm);
	  FREE_DVECTOR(BMRVals, 1, 5);
	  return;
	}


      /* Find the BMD */
      if (fabs(BMR-p[3]) < 1e-30)
	thedose = 0.0;
      else
	thedose=power(fabs((BMR-p[3])/p[4]), 1/p[5]);

    } /* end else */

  *BMD = thedose;
  OUTPUT_BENCHMD2(1,(*BMD));

  bDose = thedose;
  BMR = temp;
  /* BMRVals[1] = BMR; */

  /* output bmd to .002 file */
  /* if (bmdparm.risk==1) react=fabs(bmdparm.effect)*sign+p[3];
     else react=bmdparm.effect*p[3]; */

  react = p[3] + p[4]*pow(bDose,p[5]);
#ifndef RBMDS
  fprintf (fp_out2, "\n  RSL \t%f",react);
  fprintf (fp_out2, "\n  BMD \t%f",bDose);

  fprintf (fp_out2,"\n\n BMD_line");
  fprintf (fp_out2,"\n %f %f", -1.0, react);
  fprintf (fp_out2,"\n %f %f", bDose, react);
  fprintf (fp_out2,"\n %f %f", bDose, -0.1);
#endif
  if (bDose <= 0.0)
    {
      if ((bmr_type == 2) && (p[3] < 1e-30))
	fprintf (fp_out, "\n\nERROR: Control = 0 and Relative Risk was selected\n");
#ifndef RBMDS
      fprintf (fp_out2, "\n\n BMDL_comput_ind %d", No);  //computation failed.
#endif
      Warning("BMDL computation failed.");
      FREE_DVECTOR(pBak, 1, nparm);
      FREE_DVECTOR(BMRVals, 1, 5);
      return;
    } /* end if */
  if (bDose > 1000*xmax)
    {
#ifndef RBMDS
      fprintf (fp_out2, "\n\n BMDL_comput_ind %d", No);  //computation failed.
#endif
      Warning("BMDL computation failed.");
      return;
    } /* end if */



  /********* search for BMDL **************************/
  tol = FMAX((*BMD)*0.001, 0.0000001);

  xa = thedose/100.0;
  xb = thedose;
  BMD_lk = xlk;       /* get the lk at BMD. */

  fb = -LR;

  fa = BMDL_func(nparm, BMD_lk, xb, p, tol);

  if (fa<0.0)
    {
#ifndef RBMDS
      fprintf (fp_out2, "\n\n BMDL_comput_ind %d", No);  /* computation failed. */
#endif
      Warning ("BMDL computation failed.");
      FREE_DVECTOR(pBak, 1, nparm);
      FREE_DVECTOR(BMRVals, 1, 5);
      return;
    }

  /* computation succeeded. */
#ifndef RBMDS
  fprintf (fp_out2, "\n\n BMDL_comput_ind %d", Yes);
#endif

  Bmdl[1] = fa;

  tempdose[1] = thedose;

  for(j=1; j<=nparm; j++)
    p[j] = pBak[j];   /* get the "old" p[]. */



  PowMeans(1, p, tempdose, bmdlmean);

  Rlevel[1]= bmdlmean[1];

#ifdef MISC_OUT
  printf("           BMDL =%14.6g\n\n", (Bmdl[1]));
#endif
  fprintf(fp_out, 
#ifndef RBMDS
	  "\n            BMDL = %-14.6g\n\n"
#else
	  "\n            BMDL = %-30.22g\n\n"
#endif
	  , Bmdl[1]);

#ifndef RBMDS
  /* output bmdl data to .002 file */
  fprintf (fp_out2, "\n\n  BMDL \t%f",Bmdl[1]);
  fprintf (fp_out2,"\n\n BMDL_line");
  fprintf (fp_out2,"\n %f %f", Bmdl[1], -0.1);
  fprintf (fp_out2,"\n %f %f", Bmdl[1], react);
#endif
  bmdlCurve = No;
  if (bmdlCurve==Yes)
    {

      Get_BMRS(p, bDose, Bmdl[1], BMRVals, sign, bmr_type);
      /* Now BMRVals[2]...BMRVals[5] contain Point BMRs for BMDL curve
	 point calculations.  */

      for (k=2; k<=5;k++)
	{
	  /**** solve the BMD ********************************/

	  for(j=1; j<=nparm; j++)
	    {
	      p[j]= pBak[j];	 /* Get back p[] */
	    }

	  xa = D = 0.0;
	  tol=0.0000001;
	  incre=1.0;

	  fa = -1;


	  /* Get a BMR value that is equal to just Mu(BMD) regardless
	     of the BMR type to simplify the root finding routine */

	  if(bmr_type == 0 || bmr_type == 1 || bmr_type == 2)
	    {
	      if(sign == -1)
		bmrtemp = -BMRVals[k];
	      else
		bmrtemp = BMRVals[k];

	      if(bmr_type == 1)
		{
		  if(p[2] != 0)
		    bmrtemp = bmrtemp*sqrt(p[1]*pow(p[3],p[2]));
		  else
		    bmrtemp = bmrtemp*sqrt(p[1]);
		}
	      else if(bmr_type == 2)
		{
		  bmrtemp = bmrtemp*p[3];
		}
	    }
	  else
	    {
	      if(bmr_type == 3)
		{
		  bmrtemp = BMRVals[k] - p[3];
		}
	      else
		bmrtemp = p[4]*BMRVals[k];
	    }

	  bmrtemp = bmrtemp + p[3];
	  /* bmrtype is now a Point BMR type
	     regardless of the user input type */

	  BMR = bmrtemp;

	  /*  BMR = BMRVals[k]; */
	  thedose=pow(fabs((BMR-p[3])/p[4]), 1/p[5]);


	  /********* search for BMDL **************************/
	  tol = FMAX((*BMD)*0.001, 0.0000001);

	  xa = thedose/100.0;
	  xb = thedose;
	  BMD_lk = xlk;       /* get the lk at BMD. */

	  fb = -LR;

	  BMR = BMRVals[k];

	  fa = BMDL_func(nparm, BMD_lk, xb, p, tol);

	  if (fa<0.0)
	    {
	      /*	fprintf (fp_out2, "\n\n BMDL_comput_ind %d", No); */
	      Warning ("BMDL computation failed for one or more point on the BMDL curve.  \n          The BMDL curve will not be plotted\n");
	      FREE_DVECTOR(pBak, 1, nparm);
	      FREE_DVECTOR(BMRVals, 1, 5);
	      return;
	    }


	  Bmdl[k] = fa;
	  tempdose[1] = thedose;

	  for(j=1; j<=nparm; j++)
	    p[j] = pBak[j];          /* get the "old" p[]. */

	  PowMeans(1, p, tempdose, bmdlmean);
	  Rlevel[k]= bmdlmean[1];
	}
    }
  else for (k=2; k<=5;k++) Bmdl[k] = Rlevel[k]= -1;

  for(j=1; j<=nparm; j++)
    p[j]= pBak[j];

#ifndef RBMDS
  fprintf (fp_out2, "\n\n BMDL_Curve_flag \t %d  \n smooth_opt  %d", bmdlCurve, smooth);
  fprintf (fp_out2,"\n\n BMDL_curve");
  fprintf (fp_out2,"\n 0.00000 %f", p[3]);
  for (k=1;k<=5;k++)
    if (Bmdl[k] < xmax)
      {
	fprintf (fp_out2,"\n %f %f", Bmdl[k], Rlevel[k]/* +Parms[3] */);
      } /* end if*/
#endif

  FREE_DVECTOR(pBak, 1, nparm);
  FREE_DVECTOR(BMRVals, 1, 5);
}  /*end Pow_BMD*/

void Get_BMRS(double *p, double Dose, double bmdl1, double *BMRVals, int sign, int bmr_type)
/*********************************************************************
 *
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
{
  double range, mn;
  double ymax, ymin, bmr_percentage;
  double Percent_Change[4] = {.1, .225, .375, .5};
  int i;

  ymax = Ym[1];
  ymin = Ym[1];

  /* Find the maximum and minimum observed mean to get BMR's in terms
     of a percentage of the range of observed data */

  for(i = 2; i <= Nobs; i++)
    {
      if(Ym[i] > ymax)
	ymax = Ym[i];
      if(Ym[i] < ymin)
	ymin = Ym[i];
    }

  range = ymax - ymin;	/* The range of means */
  /* bmrtype as a percentage of the range */
  bmr_percentage = (BMRVals[1]-ymin)/range;


  /* If adverse direction is "down", then make this as a percentage
     of the range, but the pecent decrease over the maximum observed
     mean as opposed to the percent increase over the minimum */

  if(sign == -1)
    bmr_percentage = 1 - bmr_percentage;

  if(sign == 1)
    mn = ymin;
  else
    mn = ymax;

  /* Find a set of  Point type BMRs that are spaced throughout
     the range, and such that they are not too crowded by the already
     computed BMDL */

  if(.05 < bmr_percentage && bmr_percentage < .10)
    {
      Percent_Change[0] = bmr_percentage + .05;
    }
  else
    {
      if(.2 < bmr_percentage && bmr_percentage < .25)
	{
	  Percent_Change[1] = bmr_percentage + .05;
	}
      else
	{
	  if(.35 < bmr_percentage && bmr_percentage < .4)
	    {
	      Percent_Change[2] = bmr_percentage + .05;
	    }
	  else
	    {
	      if(.475 < bmr_percentage && bmr_percentage < .525)
		{
		  Percent_Change[3] = bmr_percentage + .05;
		}
	    }
	}
    }

  /* Get BMR values and give them back as the BMR type specified by
     the user */

  for(i = 0; i <= 3; i++)
    {
      BMRVals[i+2] = mn + sign*Percent_Change[i]*range;

      if(bmr_type == 0)	/* Absolute risk */
	{
	  BMRVals[i+2] = fabs(BMRVals[i+2] - p[3]);
	}
      else
	{
	  if(bmr_type == 1)		/* Std Dev risk */
	    {
	      BMRVals[i+2] = fabs(BMRVals[i+2] - p[3]);
	      BMRVals[i+2] = BMRVals[i+2]/sqrt(exp(p[1] + log(p[3]) * p[2]));
	    }
	  else
	    {
	      if(bmr_type == 2)		/* Relative risk */
		{
		  BMRVals[i+2] = fabs(BMRVals[i+2] - p[3])/p[3];
		}
	      else
		{
		  if(bmr_type == 4)	/* Extra risk */
		    {
		      BMRVals[i+2] = (BMRVals[i+2] - p[3])/p[4];
		    }
		  else		/* Point risk */
		    BMRVals[i+2] = BMRVals[i+2];
		}
	    }
	}
    }

} /* end GetBMRS */

void PowMeans(int nobs, double p[], double Doses[], double means[])
/****************************************************
/   Given the number of dose levels, this function
/	calculates the mean at each dose level and stores
/	it in the array means[]
*****************************************************/
{
  int i;

  for(i = 1; i <= nobs; i++)
    means[i] = p[3] + p[4]*pow(Doses[i],p[5]);

} /* powmeans */

/***************************************************
 *	OUTPUT_BENCHMD2--output specified benchmark dose.
 ****************************************************/
void OUTPUT_BENCHMD2(int pdcol, double BMD)
{

  OUTPUT_TEXT(" \n\n               Benchmark Dose Computation");

  /* output to the screen and to bmdswrk.002 temp file */
#ifdef MISC_OUT
  printf("Specified effect =%14.6g\n\n", bmdparm.effect);
#endif
  fprintf(fp_out, "\nSpecified effect =%14.6g\n\n", bmdparm.effect);

  if (bmr_type == 0)
    {
#ifdef MISC_OUT
      printf("Risk Type        =     Absolute risk \n\n");
#endif
      fprintf(fp_out,"Risk Type        =     Absolute risk \n\n");
    }
  else if (bmr_type == 1)
    {
#ifdef MISC_OUT
      printf("Risk Type        =     Estimated standard deviations from the control mean \n\n");
#endif
      fprintf(fp_out,"Risk Type        =     Estimated standard deviations from the control mean \n\n");
    }
  else if (bmr_type == 2)
    {
#ifdef MISC_OUT
      printf("Risk Type        =     Relative risk \n\n");
#endif
      fprintf(fp_out,"Risk Type        =     Relative risk \n\n");
    }
  else if (bmr_type == 3)
    {
#ifdef MISC_OUT
      printf("Risk Type        =     Point risk \n\n");
#endif
      fprintf(fp_out,"Risk Type        =     Point risk \n\n");
    }
  else
    {
#ifdef MISC_OUT
      printf("Risk Type        =     Extra risk \n\n");
#endif
      fprintf(fp_out,"Risk Type        =     Extra risk \n\n");
    }

#ifdef MISC_OUT
  printf("Confidence level =%14.6g\n\n",bmdparm.level);
  printf("             BMD =%14.6g\n\n",BMD);
#endif
  fprintf(fp_out, "Confidence level =%14.6g\n\n",bmdparm.level);
  fprintf(fp_out, 
#ifndef RBMDS
	  "             BMD = %-14.6g\n\n"
#else
	  "             BMD = %-30.22g\n\n"
#endif
	  ,BMD);

}	/* end OUTPUT_BENCHMD2 */

/***********************************************************
 *	Given a vector of parameter values, and the number of
 *	parameters in that vector, this function will return three
 *	new parameter values to restart the optimization if a "bad"
 *	completion code is returned from GETCL(), using a uniform
 *	random number centered at p[i]
 ***********************************************************/
void GetMoreParms(double *p, int size)
{


  if (Spec[1] != 1)
    p[0] = 1*((double)rand() / (double) RAND_MAX);
  if (Spec[2] != 1)
    p[1] = 1*((double)rand() / (double) RAND_MAX);
  if (Spec[3] != 1)
    p[2] = 1*((double)rand() / (double) RAND_MAX);
  if (Spec[4] != 1)
    {
      if (p[3] >= 0)
	p[3] = 1*((double)rand() / (double) RAND_MAX);
      else
	p[3] = -1*((double)rand() / (double) RAND_MAX);
    } /* end if */
  if (Spec[5] != 1)
    {
      if (restrict == 1)
	p[4] = 1 + ((double)rand() / (double) RAND_MAX);
      else
	p[4] = ((double)rand() / (double) RAND_MAX);
    } /* end if */

}	/* end GetMoreParms */

/***********************************************************
 *	Given a vector of parameter values, and the number of
 *	parameters in that vector, this function will return three
 *	new parameter values to restart the optimization if a "bad"
 *	completion code is returned from GETCL(), using a uniform
 *	random number centered at p[i]
 ***********************************************************/
void GetMLEParms(double *p, int size)
{


  if (Spec[1] != 1)
    p[0] = p[0] + 2*((double)rand() / (double) RAND_MAX);
  if (Spec[2] != 1)
    p[1] = p[1] + 1*((double)rand() / (double) RAND_MAX);
  if (Spec[3] != 1)
    p[2] = p[2] + 2*((double)rand() / (double) RAND_MAX);
  if (Spec[4] != 1)
    {
      if (p[3] >= 0)
	p[3] = p[3] + 2*((double)rand() / (double) RAND_MAX);
      else
	p[3] = p[3] -2*((double)rand() / (double) RAND_MAX);
    } /* end if */
  if (Spec[5] != 1)
    {
      if (restrict == 1)
	p[4] = p[4] + 1 + ((double)rand() / (double) RAND_MAX);
      else
	p[4] = p[4] + ((double)rand() / (double) RAND_MAX);
    } /* end if */

}	/* end GetOtherParms */

/***********************************************************
 *	Given a vector of parameter values, and the number of
 *	parameters in that vector, this function will return three
 *	new parameter values to restart the optimization if a "bad"
 *	completion code is returned from GETCL(), using a uniform
 *	random number centered at p[i]
 ***********************************************************/
void GetOtherParms(double *p, int size)
{


  if (Spec[1] != 1)
    p[0] = p[0] + .2*((double)rand() / (double) RAND_MAX);
  if (Spec[2] != 1)
    p[1] = p[1] + .1*((double)rand() / (double) RAND_MAX);
  if (Spec[3] != 1)
    p[2] = p[2] + .2*((double)rand() / (double) RAND_MAX);
  if (Spec[4] != 1)
    {
      if (p[3] >= 0)
	p[3] = p[3] + .2*((double)rand() / (double) RAND_MAX);
      else
	p[3] = p[3] -.2*((double)rand() / (double) RAND_MAX);
    } /* end if */
  if (Spec[5] != 1)
    {
      if (restrict == 1)
	p[4] = p[4] + .1 + ((double)rand() / (double) RAND_MAX);
      else
	p[4] = p[4] + ((double)rand() / (double) RAND_MAX);
    } /* end if */

}	/* end GetOtherParms */


/********************************************************
 * power function to safely compute a^b
 *********************************************************/
double power(double a, double b)
{
  double log_a, blog_a, a_b;

  if (a >= 1.0e-10)
    log_a = log(a);
  else
    log_a = -24.5258509299404572 + (2.0e10)*a - (5.0e19)*a*a;

  blog_a = b*log_a;
  if (blog_a > 700)
    blog_a = 700.0;

  a_b = exp(blog_a);
  return a_b;
} /* end power */


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
