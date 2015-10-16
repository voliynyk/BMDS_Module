#include "benchmark.h"
#include "in_outfun.h"

/***********************************************************
 *Goodness -- Used to test the Goodness of Fit
 *
 * Modified By: Micheal Ferree
 * Date: 02JUN05
 * Reason: Fixed DF output.
 *
 ************************************************************/

void Goodness(int nparm, int nparm_known, double Parms[], int type, AnaList anasum[])
{
  int df1, df2, df3, df4;
  char small_pv[2][9] = {"<.0001", "NA"};

  df4 = anasum[2].DF - anasum[4].DF;
  df1 = anasum[2].DF - anasum[1].DF;
  df2 = anasum[2].DF - anasum[3].DF;
  df3 = anasum[3].DF - anasum[5].DF;

  OUTPUT_TEXT("\n\n                   Explanation of Tests  ");
  fprintf(fp_out, "\n Test 1:  Do responses and/or variances differ among Dose levels? \n          (A2 vs. R)");
  fprintf(fp_out, "\n Test 2:  Are Variances Homogeneous? (A1 vs A2)");
  fprintf(fp_out, "\n Test 3:  Are variances adequately modeled? (A2 vs. A3)");
  fprintf(fp_out, "\n Test 4:  Does the Model for the Mean Fit? (A3 vs. fitted)");
  fprintf(fp_out, "\n (Note:  When rho=0 the results of Test 3 and Test 2 will be the same.)");

  OUTPUT_TEXT("\n\n                     Tests of Interest    ");
  OUTPUT_TEXT("\n   Test    -2*log(Likelihood Ratio)  Test df        p-value    ");

  if(anasum[4].TEST >= .0001 && df4 > 0)
    {
      fprintf(fp_out, "\n   Test 1 %20.6g %10d %15.4g", anasum[4].MSE, df4,anasum[4].TEST);
    }
  else
    {
      /* Check degrees of freedom.  Since Chi-Square function returns 0 if
	 degrees of freedom is less than or equal to zero, checking the degrees
	 of freedom here is valid.  If there is an invalid number for the d.f.
	 then output *** as the p-value */
      if(df4 > 0)
	{
	  fprintf(fp_out, "\n   Test 1 %20.6g %10d %15s", anasum[4].MSE,df4, small_pv[0]);
	}
      else
	{
	  fprintf(fp_out, "\n   Test 1 %20.6g %10d %15s", anasum[4].MSE,df4, small_pv[1]);
	}
    }
  if(anasum[1].TEST >= .0001 && df1 > 0)
    {
      fprintf(fp_out, "\n   Test 2 %20.6g %10d %15.4g", anasum[1].MSE, df1,anasum[1].TEST);
    }
  else
    {
      if(df1 > 0)
	{
	  fprintf(fp_out, "\n   Test 2 %20.6g %10d %15s", anasum[1].MSE,df1, small_pv[0]);
	}
      else
	{
	  fprintf(fp_out, "\n   Test 2 %20.6g %10d %15s", anasum[1].MSE,df1, small_pv[1]);
	}
    }
  if(anasum[2].TEST >= .0001 && df2 > 0)
    {
      fprintf(fp_out, "\n   Test 3 %20.6g %10d %15.4g", anasum[2].MSE, df2,anasum[2].TEST);
    }
  else
    {
      if(df2 > 0)
	{
	  fprintf(fp_out, "\n   Test 3 %20.6g %10d %15s", anasum[2].MSE,df2, small_pv[0]);
	}
      else
	{
	  fprintf(fp_out, "\n   Test 3 %20.6g %10d %15s", anasum[2].MSE,df2, small_pv[1]);
	}
    }
  if(anasum[3].TEST >= .0001 && df3 > 0)
    {
      fprintf(fp_out, "\n   Test 4 %20.6g %10d %15.4g", anasum[3].MSE, df3,anasum[3].TEST);
    }
  else
    {
      if(df3 > 0)
	{
	  fprintf(fp_out, "\n   Test 4 %20.6g %10d %15s", anasum[3].MSE,df3, small_pv[0]);
	}
      else
	{
	  fprintf(fp_out, "\n   Test 4 %20.6g %10d %15s", anasum[3].MSE,df3, small_pv[1]);
	}
    }

  /* Output conclusion for Test 1. */
  if(df4 <= 0)
    {
      fprintf(fp_out, "\n\nNA - Degrees of freedom for Test 1 are less than or equal to 0.  The Chi-Square\n");
      fprintf(fp_out, "     test for fit is not valid");
    }
  else if(anasum[4].TEST < .05)
    {
      fprintf(fp_out, "\n\nThe p-value for Test 1 is less than .05.  There appears to be a\n");
      fprintf(fp_out, "difference between response and/or variances among the dose levels");
      fprintf(fp_out, "\nIt seems appropriate to model the data");
    }
  else
    {
      fprintf(fp_out, "\n\nThe p-value for Test 1 is greater than .05.  There may not be a\n");
      fprintf(fp_out, "diffence between responses and/or variances among the dose levels");
      fprintf(fp_out, "\nModelling the data with a dose/response curve may not be appropriate");
    }

  /* Output conclusion for Test 2. */
  if(type == 1)
    {
      if(df1 <= 0)
	{
	  fprintf(fp_out, "\n\nNA - Degrees of freedom for Test 2 are less than or equal to 0.  The Chi-Square\n");
	  fprintf(fp_out, "     test for fit is not valid");
	}
      else if(anasum[1].TEST < .1)
	{
	  fprintf(fp_out, "\n\nThe p-value for Test 2 is less than .1.  A non-homogeneous variance \nmodel appears to be appropriate");
	}
      else
	{
	  fprintf(fp_out, "\n\nThe p-value for Test 2 is greater than .1.  Consider running a \nhomogeneous model");
	}
    }
  else
    {
      if(df1 <= 0)
	{
	  fprintf(fp_out, "\n\nNA - Degrees of freedom for Test 2 are less than or equal to 0.  The Chi-Square\n");
	  fprintf(fp_out, "     test for fit is not valid");
	}
      else if(anasum[1].TEST < .1)
	{
	  fprintf(fp_out, "\n\nThe p-value for Test 2 is less than .1.  Consider running a \nnon-homogeneous variance model");
	}
      else
	{
	  fprintf(fp_out, "\n\nThe p-value for Test 2 is greater than .1.  A homogeneous variance \nmodel appears to be appropriate here\n");
	}
    }

  /* Output conclusion for Test 3. */
  if(df2 <= 0)
    {
      fprintf(fp_out, "\n\nNA - Degrees of freedom for Test 3 are less than or equal to 0.  The Chi-Square\n");
      fprintf(fp_out, "     test for fit is not valid");
    }
  else if(anasum[2].TEST < .1)
    {
      fprintf(fp_out, "\n\nThe p-value for Test 3 is less than .1.  You may want to consider a \ndifferent variance model");
    }
  else
    {
      fprintf(fp_out, "\n\nThe p-value for Test 3 is greater than .1.  The modeled variance appears \n to be appropriate here");
    }

  /* Output conclusion for Test 4. */
  if(df3 <= 0)
    {
      fprintf(fp_out, "\n\nNA - Degrees of freedom for Test 4 are less than or equal to 0.  The Chi-Square\n");
      fprintf(fp_out, "     test for fit is not valid");
    }
  else if(anasum[3].TEST < .1)
    {
      fprintf(fp_out, "\n\nThe p-value for Test 4 is less than .1.  You may want to try a different \nmodel");
    }
  else
    {
      fprintf(fp_out, "\n\nThe p-value for Test 4 is greater than .1.  The model chosen seems \nto adequately describe the data");
    }

  fprintf(fp_out,"\n");

}
