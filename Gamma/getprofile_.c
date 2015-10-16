/********************************************
 *  getprofile_       1/10/01               *
 *  This function is used in the discrete   *
 *  models to call Max_lk instead of donlp2.*
 ********************************************/

extern int nparm;         // This is the # of parameters

void MAX_lk(int , double [], double , int *, double *,
	    int *, double);


void getprofile_(double parms[], double *BMD, double *Max_LL,
		 long int *optite, long int *tmpcnt, long int *BMRTYPE,
		 double *BMR2) 
{
  int junk;
  int eflag = -9999;
  double tDose, gtol = 1e-7;
 

  junk = 0;
  tDose = *BMD;
  MAX_lk(nparm, parms, gtol, &junk, Max_LL, &eflag, tDose); 
  *optite = eflag;

}
