#include <stdlib.h>
#include <float.h>
#include "benchmark.h"
#include "allo_memo.h"




void SortByLs (int, int, int [], double [], double [], double [], double [], double []);

/*************************************************
*SRoI - Used to find the scaled residual of interest
*
**************************************************/
void SRoI (int ngrp, int Nobs, double GXi[], double Xi[], int Xg[], double SR[], double LSC[], double meanLSC, double BMD)
{
  int i, j;
  int locDose = 0;	/*location of closest dosegroup */
  int locLSC = 0;       /*location of closest LSC(first found) */
  int litSR = 1;    	/*number of litters with closest LSC*/
  double maxSR, minSR, avgSR, maxabsSR, minabsSR, avgabsSR;
  double closeLSC = 0;
  double closeDose = 0;
  double idiff, diff;


  /*** choose dose group closest to BMD */
  diff = DBL_MAX;
  for (i=1; i<=ngrp; i++)
    {
    idiff = fabs(BMD - GXi[i]);
    if (idiff < diff)
      {
      diff = idiff;
      locDose = i;
      closeDose = GXi[i];
      }
    }


  /*** choose LSC closest to mean LSC */
  diff = DBL_MAX;
  for (i=1; i<=Nobs; i++)
    {
    if (Xi[i] == closeDose)  // picks out dose group that is closest to BMD
      {
      idiff = fabs(LSC[i] - meanLSC);
      if (idiff == diff)
        litSR++;
      if (idiff < diff)
        {
        litSR = 1;             //reset litter count
        diff = idiff;
        closeLSC = LSC[i];     //value of closest LSC
        locLSC = i;        
        }
      }
    }
  

  /*** calculate max, min, average SRoI and |SRoi|*/
  maxSR = SR[locLSC];
  minSR = SR[locLSC];
  avgSR = SR[locLSC];
  maxabsSR = fabs(SR[locLSC]);
  minabsSR = fabs(SR[locLSC]);
  avgabsSR = fabs(SR[locLSC]);

  if (litSR != 1)
    {
    avgSR = 0;
    avgabsSR = 0;
    for (i=locLSC; i<=locLSC+litSR-1; i++)                //relies on dose groups being sorted by LSC
      {
      if (Xi[i] == GXi[locDose] && LSC[i] == closeLSC)
        {
        if (SR[i] > maxSR) maxSR = SR[i];
        if (SR[i] < minSR) minSR = SR[i];
        if (fabs(SR[i]) > maxabsSR) maxabsSR = fabs(SR[i]);
        if (fabs(SR[i]) < minabsSR) minabsSR = fabs(SR[i]);
        avgSR += SR[i];
        avgabsSR += fabs(SR[i]);
        printf("SR[i] = %9.4f, avgabsSR = %9.4f\n", SR[i], avgabsSR);
        }
      }
    avgSR /= litSR;
    avgabsSR /= litSR;
    }




  fprintf(fp_out, "\n\n\nScaled Residual(s) for Dose Group Nearest the BMD");
  fprintf(fp_out, "\n------------------------------");


  fprintf(fp_out, "\nMinimum scaled residual for dose group nearest the BMD =      %9.4f", minSR);
  fprintf(fp_out, "\nMinimum ABS(scaled residual) for dose group nearest the BMD = %9.4f", minabsSR);
  fprintf(fp_out, "\nAverage scaled residual for dose group nearest the BMD =      %9.4f", avgSR);
  fprintf(fp_out, "\nAverage ABS(scaled residual) for dose group nearest the BMD = %9.4f", avgabsSR);
  fprintf(fp_out, "\nMaximum scaled residual for dose group nearest the BMD =      %9.4f", maxSR);
  fprintf(fp_out, "\nMaximum ABS(scaled residual) for dose group nearest the BMD = %9.4f", maxabsSR);
  fprintf(fp_out, "\nNumber of litters used for scaled residual for dose group nearest the BMD = %d", litSR);




  
}