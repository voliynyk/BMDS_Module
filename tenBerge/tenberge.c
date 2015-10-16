// *******************************************************************************
//    Converted to C from VisualBasic program of Dosresp (Dose-Response Analysis)
//    software (W.F. ten Berge, Nov. 2001) without the visual components
//    
//    by Qun He (12/26/2006)
// ****************************************************************

// ***********************************************************************
//               General information of the model
// ***********************************************************************

char Version_no[] = "Ten Berge Model. (Version: 1.0; Date: 12/26/2006)";
char Model_info[] = "Dose-Response Analysis\n\n Method of Maximum Likelihood according to: \
					\n D.J. Finney, 1977. Probit Analysis. Cambridge University Press.";
char Formula[] = "Model: P(v1, v2, ...) = Link(B0 + B1*v1 + B2*v2 + ...)\n \
                  \n   Link is either Logit or Probit \
				  \n   v1, v2, ... are the variables (transformations of the input parameters)\n";


// ********************************************
//      Header files
// ********************************************
#include <windows.h>    // Win32 Header File 
#include <windowsx.h>   // Win32 Header File 
#include <commctrl.h>   // Win32 Header File 
#include <mmsystem.h>   // Win32 Header File 
#include <shellapi.h>   // Win32 Header File 
#include <shlobj.h>     // Win32 Header File 
#include <richedit.h>   // Win32 Header File 
#include <wchar.h>      // Win32 Header File 
#include <objbase.h>    // Win32 Header File 
#include <ocidl.h>      // Win32 Header File 
#include <winuser.h>    // Win32 Header File 
#include <olectl.h>     // Win32 Header File 
#include <conio.h>
#include <direct.h>
#include <ctype.h>
#include <io.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#include <setjmp.h>
#include <time.h>
#include <stdarg.h>
#include <process.h>

#include <benchmark.h>
#include <ERRORPRT.h>
#include <allo_memo.h>
#include <matrix_agb.h>
#include <specialfun.h>
#include <computation.h>
#include <in_outfun.h>

// *************************************************
//                System Variables
// *************************************************

char    BCX_STR [1024*1024];
jmp_buf GosubStack[32];
int     GosubNdx;
#define DEBUG 1;

// *************************************************
//            User Global Variables
// *************************************************

static int     I;
static int     J;
static int     SW;
static FILE*  FP2;
static FILE* logfile;
static FILE* response;

#define NUMOFVAR 14;           /* See variable Vrx[][] */ 
#define LENOFNAME 20;          /* See variable Vrx[][] */

int NZ;                        /* No. of variables */
int NW;                        /* No. of observations */
int NY;                        /* No. of selected single transformed variable for analysis */
int NC;                        /* No. of pairs of product of transformed variables */
int NX;                        /* Total No. (selected single variable + selected pair of variables) */
int KBZ;                       /* (1=without; 2=with) background response correction */
int KPZ;                       /* (1=Probit; 2=Logit) model */
int KZ, KZY, HQ, MQ, NXT;
int NB;                        /* Flag of error in transform variables */
int MV;                        /* Index of requested variable for evaluation */
int NS, MX, DF;
int N1;                        /* First sequence number of trials for analysis */
int N2;                        /* Last sequence number of trials for analysis */
int TW;
int NXC;                       /* No. of variables selected for product analysis */
double CH;                     /* Chi-Squares value */
double CX;                     /* Correction for variances (Chi-Squares/Degrees of Freedom) */
double ZM, P;
double PW;                     /* Requested response in percentage (set to 50 when user specify <= 0) */
double Q, SQ, Y, Z;
double XZ; 
double TX;                     /* Value of Student t set by user */
double CZ, DC, DW;
double AW, BW, CW, Cy, RA, RB;
int K0[10];                    /* List of selected single transformed variable for analysis */
int K1[5];                     /* List of the left variables for the product analysis */
int K2[5];                     /* List of the right variables for the product analysis */
int WF[12];
int NT[10];                    /* List of variables' type of transformation (logrithmic, reciprocal, none) */
double XU[10];                 /* List of transformed values for each variable */
double XW[10], PA[10];
double AX[10], BX[10], FX[10];
double ZA[10][10], BV[10][10], ZC[10][10], ZB[10];
int TR;                        /* Number of variables chosen for ratio regression coefficient evaluation */
int NRC1;                      /* First variable chosen for ration regression coefficient evaluation */
int NRC2;                      /* Second variable chosen for ration regression coefficient evaluation */
int FlEr;                      /* Error flag */
char Fina[];                /* Input file name */
int *N;                        /* No. of subjects array */
int *LR;                       /* No. of responders array */
double **X;                    /* Input data matrix */
char Vrx[14][40];              /* Input data column names */
char conc_col[CNLENGTH];
char expos_time_col[CNLENGTH];
char subj_prop_col[CNLENGTH];
char num_expos_col[CNLENGTH];
char num_respd_col[CNLENGTH];
char comment_line[128];
int num_missing_value;         /* No. of records with missing values */
int stdchisqr;                 /* Flag indicates Chisquare() calculated */
int chisqr_skip = 1;           /* Flag the print out of Chisquare() calculattion */
time_t  ltime;

// *************************************************
//       Define input and output files's name  
// *************************************************

char fout[FLENGTH];				/* output temp file */
char fout2[FLENGTH];
char plotfilename[FLENGTH];		/* file to send to GnuPlot */
char infilestem[FLENGTH];       /* input file name stem */
char respgraphfile[FLENGTH];     /* graphics file name */
  
// *************************************************
//               Standard Macros
// *************************************************

#define VAL(a)(double)atof(a)


// *************************************************
//               Standard Prototypes
// *************************************************

char*   BCX_TmpStr(size_t);
char*   str (double);
char*   join (int, ... );
double  Abs (double);
double  Exp (double);

// *************************************************
//               User's Prototypes
// *************************************************

void    RetrDat (void);
void    SuPPDat (void);
void    Filcomb3 (void);
void    nrm21 (int);
void    nrm22 (void);
void    nrm23 (void);
void    nrm24 (void);
void    nrm25 (int);
void    nrm26 (void);
void    nrm28 (void);
void    nrm29 (void);
void    pqa (void);
void    pqb (void);
void    matrinv (double A[10][10], double B[10][10]);
void    matrix (double A[10][10], double B[10], double C[10]);
void    CalcML (void);
void    CalcDoseResponse(void);
void    CalcResponseDose(void);
void    CalcResponseGraph(void);
void    Graphic2(double Dx,  double *Px1, double *Px2, double *Px3, double SV[]);
void    CalcRatio(void);
void    nrm41 (int);
void    nrm42 (void);
void    nrm43 (void);
void    nrm44 (void);
void    Chisquare (void);
void    Warn1 (double QZ);
void    Warn2 (double QZ);
void    WriDat (void);
void    CalcErr (void);
void    check_args(int argc, char *argv[]);
void    OPEN_FILES(int argc, char *argv[]);
void    read_modeling_data(void);
int     READ_OBSDATA (int r, int c, double **m);
void    output_read_data(void);
void    CLOSE_FILES(void);
void    calML_H1(void);
void    NRM2BGR(void);
void    Get_File_Stem(char *argv, char *infilestem);
void    Derive_File_Name(char *stem, char *newfile, const char *ext);

int     check_data_sectors(const char *sector);


// *************************************************
//                  Main Program
// *************************************************

int main(int argc, char *argv[])
{
  char sector[20];
  int compare;
  time(&ltime);

  check_args(argc, argv);
  Get_Names(argv[1], fout, fout2, plotfilename);    /* get input and output file names */
  Get_File_Stem(argv[1], infilestem);               /* get input file name stemp */
  Derive_File_Name(infilestem, respgraphfile, "Asc");/* derive graphics file name from infile stem */              
  OPEN_FILES(argc, argv);                           /* open output files */
  read_modeling_data();                                      /* read input file */

  Output_Header(Version_no, argv[1], plotfilename, ctime(&ltime), Model_info);
  output_read_data();                               /* output some input data */
  
  CalcML();
  WriDat();

  fscanf(fp_in, "%s", sector);
  compare = check_data_sectors(sector);

  while (compare != 0) {
	  if (compare == 1)
		  CalcDoseResponse();
	  else if (compare == 2)
		  CalcResponseDose();
	  else if (compare == 3)
		  CalcRatio();
	  else if (compare == 4)
		  CalcResponseGraph();

	  fscanf(fp_in, "%s", sector);
	  compare = check_data_sectors(sector);
  }

  FREE_IVECTOR (N,1,NW);
  FREE_IVECTOR (LR,1,NW);
  FREE_DMATRIX(X, 1, NW, 1, NZ+2);
  CLOSE_FILES();                                    /* close all opened files */

  return 0;
}                                                   /* End of main program */


// *************************************************
//                 Run Time Functions
// *************************************************

char *BCX_TmpStr (size_t Bites)
{
  static int   StrCnt;
  static char *StrFunc[2048];
  StrCnt=(StrCnt + 1) & 2047;
  if(StrFunc[StrCnt]) free (StrFunc[StrCnt]);
  return StrFunc[StrCnt]=(char*)calloc(Bites+128,1); 
}


char *str (double d)
{
  register char *strtmp = BCX_TmpStr(16);
  sprintf(strtmp,"% .15G",d);
  return strtmp;
}

char * join(int n, ...)
{
  register int i = n, tmplen = 0;
  register char *s_;
  register char *strtmp = 0;
  va_list marker;
  va_start(marker, n); // Initialize variable arguments
  while(i-- > 0)
  {
    s_ = va_arg(marker, char *);
    tmplen += strlen(s_);
  }
  strtmp = BCX_TmpStr(tmplen);
  va_end(marker); // Reset variable arguments
  i = n;
  va_start(marker, n); // Initialize variable arguments
  while(i-- > 0)
  {
    s_ = va_arg(marker, char *);
    strtmp = strcat(strtmp, s_);
  }
  va_end(marker); // Reset variable arguments
  return strtmp;
}


double Exp (double arg)
{
  return pow(2.718281828459045,arg);
}


double Abs (double a)
{
  if(a<0) return -a;
  return  a;
}

// ************************************
//       User Functions
// ************************************

void check_args(int argc, char *argv[]) {
  if(argc == 2)
	  show_version(argv[1], Version_no);

  if(argc < 2) {
      fprintf(stderr, "ERROR:  Requires two arguments\nUsage:  %s <file.(d)>\n", argv[0]);
      fprintf (stderr, "   or:  %s -v for version number.\n", argv[0]);
      exit (1);
  }
}

int check_data_sectors(const char *sector) 
{
	int compare = strcmp("end", sector);
		
	if (compare == 0)
		return 0;

	compare = strcmp("dose", sector);
		
	if (compare == 0)
		return 1;

	compare = strcmp("response", sector);

	if (compare == 0)
		return 2;

	compare = strcmp("ratio", sector);

	if (compare == 0)
		return 3;

	compare = strcmp("graph/response", sector);

	if (compare == 0)
		return 4;
}

// *****************************************************
//  OPEN_FILES--used to open input and output files.
// *****************************************************

void OPEN_FILES (int argc, char *argv[]) {
  logfile = fopen("log.txt", "a");  /* open log file */
  fp_in = fopen(argv[1], "r");      /* open input file */

  // open output files 
  fp_out = fopen(fout, "w");	    /* overwrite output */
  fp_out2 = fopen(fout2, "w");	    /* always overwrite plotting file */

  // check to make sure files are open, if not, print error message and exit 
  if (fp_in == NULL || fp_out == NULL || fp_out2 == NULL) {
      fprintf(logfile, "Error in opening input and output files.");
      ERRORPRT ("Error in opening input and output files.");
  }

  response = fopen(respgraphfile, "w");  /* open response graph file */
}

// *****************************************************
//  CLOSE_FILES--used to close input and output files.
// *****************************************************

void CLOSE_FILES (void) {
  if (fclose(fp_in) != 0 || fclose(fp_out) != 0 || fclose(fp_out2) != 0) {
	  fprintf(logfile, "Error in closing opened files.");
	  fclose(response);
	  fclose(logfile);
	  ERRORPRT ("Error in closing opened files.");
  }

  fclose(response);
  fclose(logfile);
}


void read_modeling_data () {
  char junk[255]; 
  fscanf(fp_in, "%d", &NZ);                       /* begin reading input from input .(d) file */
  fscanf(fp_in, "%d", &NW);
  fgets(junk, 128, fp_in);                        /* to get the junk out of the previous line */   

  int i;
  Vrx[0][0] = 0;
  for (i = 1; i < 14; i++) {
	  if (i <= NZ + 2) {
		  fgets(Vrx[i], 40, fp_in);
		  Vrx[i][strlen(Vrx[i]) - 1] = 0;         /* strip carriage return */
	  } else {
		  Vrx[i][0] = 0;
	  }
  }

  #ifdef DEBUG
  fprintf(logfile, "\n########################### %s\n", ctime(&ltime));
  fprintf(logfile, "Reading data from input .(d) file...\n\n");
  fprintf(logfile, "number of variables: %d \n", NZ);
  fprintf(logfile, "number of observations: %d \n", NW);
  fprintf(logfile, "concentration column: %s \n", Vrx[1]);
  fprintf(logfile, "exposure time: %s \n", Vrx[2]);
  fprintf(logfile, "subject property: %s \n", Vrx[3]);
  fprintf(logfile, "number of exposed: %s \n", Vrx[4]);
  fprintf(logfile, "number of reponded: %s \n", Vrx[5]);
  #endif

  N = IVECTOR(1, NW);
  LR = IVECTOR(1, NW);
  X = DMATRIX (1, NW, 1, NZ);          /* init data variables */
  num_missing_value = READ_OBSDATA(NW, NZ, X);

  #ifdef DEBUG
  fprintf(logfile, "\n\nNumber of missing value = %d\n", num_missing_value);
  fprintf(logfile, "Reading input data file finished.\n");
  #endif

  fscanf(fp_in, "%s", comment_line);
  fscanf(fp_in, "%d%d", &KPZ, &KBZ);

  for (i = 1; i <= NZ; i++) {
	  fscanf(fp_in, "%d", &NT[i]);
  }

  fscanf(fp_in, "%d", &NY);

  for(i = 1; i <= NY; i++) {
	  fscanf(fp_in, "%d", &K0[i]);
  }

  fscanf(fp_in, "%d", &NC);

  for(i = 1; i <= NC; i++) {
	  fscanf(fp_in, "%d%d", &K1[i], &K2[i]);
  }

  fscanf(fp_in, "%d%d", &N1, &N2);                /* first & last sequence number of trials for analysis */

  #ifdef DEBUG
  fprintf(logfile, "\n\nUser selections\n");
  fprintf(logfile, "%s\n", comment_line);
  fprintf(logfile, "KPZ=%d\n", KPZ);
  fprintf(logfile, "KBZ=%d\n", KBZ);
  fprintf(logfile, "NT[1]=%d\n", NT[1]);
  fprintf(logfile, "NT[2]=%d\n", NT[2]);
  fprintf(logfile, "NT[3]=%d\n", NT[3]);
  fprintf(logfile, "NY=%d\n", NY);
  fprintf(logfile, "K0[1]=%d\n", K0[1]);
  fprintf(logfile, "K0[2]=%d\n", K0[2]);
  fprintf(logfile, "K0[3]=%d\n", K0[3]);
  fprintf(logfile, "NC=%d\n", NC);
  fprintf(logfile, "K1[1]=%d\n", K1[1]);
  fprintf(logfile, "K2[1]=%d\n", K2[1]);
  fprintf(logfile, "N1=%d\n", N1);
  fprintf(logfile, "N2=%d\n", N2);
  #endif
}

void output_read_data(void) {
  fprintf(fp_out, "\n   %s", Formula);
  fprintf(fp_out, "\n\n   Number of input parameters = %d", NZ);
  fprintf(fp_out, "\n   Total number of observations = %d", NW);
  fprintf(fp_out, "\n   Total number of records with missing values = %d\n\n", num_missing_value);

  fprintf(fp_out2,"\n\n Model calculations done.");
}

void nrm21 (int JM)
{
  NB = 0;
  int IV;
  for(IV = 1; IV <= NZ; IV++) {
	  if (NT[IV] == 1) {
		  if (X[JM][IV] == 0) {
			  NB = 1;
			  return;
		  } else {
			  XU[IV] = log(X[JM][IV]);
		  }
	  }
		  
	  if (NT[IV] == 2) {
		  if (X[JM][IV] == 0) {
			  NB = 1;
			  return;
		  } else {
			  XU[IV] = 1 / X[JM][IV];
		  }
	  }
		  
	  if (NT[IV] == 3) {
		  XU[IV] = X[JM][IV];
	  }

#ifdef DEBUG
fprintf(logfile, "X[%d][%d]=%f XU[%d]=%f\n", JM, IV, X[JM][IV], IV, XU[IV]);
#endif
  }//end for(IV = 1; IV <= NZ; ++IV)
}


void nrm22 (void)
{
  int IW, NQ;

  if (NY == 0) {
      return;
  }

  for (IW = 1; IW <= NY; IW++) {
      NQ = K0[IW];
      XW[IW] = XU[NQ];
  }
}


void nrm23 (void)
{
  int IX, H1, H2;

  if (NC == 0) {
	  return;
  }

  for (IX = 1; IX <= NC; IX++) {
	  H1 = K1[IX];
      H2 = K2[IX];
      XW[NY+IX] = XU[H1] * XU[H2];
  }
}


void nrm24 (void)
{
  int I, J;
  for(I = MX; I <= NXT; I++) {
	  for(J = MX; J <= NXT; J++) {
		  ZA[I][J] = 0; 
	  }
      
	  ZB[I] = 0;  
  }

  CH = 0;
  NS = 0;
  CX = 0;
}

void nrm25 (int JM)
{
  int H, K;
	  
  if (KBZ == 1) {
	  RA = N[JM] / (P * Q);
	  RB = (LR[JM] + 0.0) / N[JM] - P;    /* we don't realy want a int / int result */
	    
#ifdef DEBUG
	  fprintf(logfile, "RA=%f RB=%f N[%d]=%d LR[%d]=%d P=%f Q=%f\n", RA, RB, JM, N[JM], JM, LR[JM], P, Q);
#endif
  } else if(KBZ == 2) {
	  RA = N[JM] / (PW * Q * SQ);
	  RB = (LR[JM] + 0.0) / N[JM] - PW;
  }

  CH = CH + RA * pow(RB,2);
  #ifdef DEBUG
  fprintf(logfile, "CH=%f\n", CH);
  #endif
  NS += 1;

  for(H = MX; H <= NXT; H++) {
	  for(K = MX; K <= NXT; K++) {
		  //ZA[H][K] = ZA[H][K] + RA * PA[H] * PA[K];
		  ZA[K][H] = ZA[K][H] + RA * PA[H] * PA[K];
	  }
    
	  ZB[H] = ZB[H] + RA * RB * PA[H];
  }
}


void nrm26 (void)
{
  int I, J;

  for (I = MX; I <= NXT; I++) {
	  for (J = MX; J <= NXT; J++) {
		  ZC[J][I] = ZA[J][I];
	  }
    }
}


void nrm28 (void)
{
  int I;
  ZM = 0;
	 
  for(I = MX; I <= NXT; I++) {
      BX[I] = BX[I] + FX[I];
      ZM = ZM + Abs(FX[I] / BX[I]);
	  //fprintf(fp_out,"\n   BX[%d]=%f",(int)I, BX[I]);
  }
  
  CX = 1;
  DF = NS - NXT + MX - 1;
}

void nrm29 (void)
{
//print out
}


void pqa (void)
{
  double S;

  if (KPZ == 1) {
	  S = 1 / (1 + 0.2316419 * Abs(Y));
      Z = Exp(-(pow(Y,2)) / 2) / sqrt(6.283185307);
      P = Z * (S * 0.31938153 - pow(S,2) * 0.356563782 + pow(S,3) * 1.781477937 - pow(S,4) * 1.821255978 + pow(S,5) * 1.330274429);
	  
	  if (Y > 0) {
		  P = 1 - P;
	  }
  } else if (KPZ == 2) {
      P = Exp(Y) / (Exp(Y) + 1);
      Z = P * (1 - P);
  }
}


void pqb (void)
{
  double QT, TQ, AQ, BQ;
  QT = PW / 100;
	  
  if (KPZ == 1) {
	  if (QT > 0.5) {
		  QT = 1 - QT;
	  }

      TQ = sqrt(log(1 / pow(QT,2)));
      AQ = 2.515517 + 0.802853 * TQ + 0.010328 * pow(TQ,2);
      BQ = 1 + 1.432788 * TQ + 0.189269 * pow(TQ,2) + 0.001308 * pow(TQ,3);
      Y = TQ - AQ / BQ;

	  if (PW < 50) {
          Y = -Y;
      }
    
  } else if (KPZ == 2) {
	  Y = log(QT / (1 - QT));
  }
}


void matrinv (double A[10][10], double BZ[10][10])
{
  int HA, HB, HS, HX, HY, I, J;
  double G;

  #ifdef DEBUG
  fprintf(logfile, "\n\nmatrinv function:\n");
  #endif

  for(I = MX; I <= NXT; I++) {
	  for(J = MX; J <= NXT; J++) {
          BZ[J][I] = 0;
		  if(I == J) {
              BZ[J][I] = 1;
		  }
	  }
  }

  HS = 1;
  HA = MX;
  HB = NXT - 1;

  for(HY = 1; HY <= 2; HY++) {
	  #ifdef DEBUG
	  fprintf(logfile, "HS=%d\n", HS);
	  fprintf(logfile, "HY=%d\n", HY);
	  #endif

	  for(HX = HA; HS >= 0 ? HX <= HB : HX >= HB; HX += HS) {
		  #ifdef DEBUG
		  fprintf(logfile, "HX=%d\n", HX);
		  #endif

		  for(I = HX+HS; HS >= 0 ? I <= HB+HS : I >= HB+HS; I += HS) {
			  #ifdef DEBUG
			  fprintf(logfile, "I=%d\n", I);
			  fprintf(logfile, "A[%d][%d]=%f\n", HX, HX, A[HX][HX]);
			  #endif
  
			  if(A[HX][HX] != 0) {
                  //G = -A[I][HX] / A[HX][HX];
				  G = -A[HX][I] / A[HX][HX];
				  for(J = MX; J <= NXT; J++) {
                      //A[I][J] = A[I][J] + G * A[HX][J];
                      //BZ[I][J] = BZ[I][J] + G * BZ[HX][J];
					  A[J][I] = A[J][I] + G * A[J][HX];
                      BZ[J][I] = BZ[J][I] + G * BZ[J][HX];
				  }
			  }
		  }
	  }
      
	  HS = -1;
      HA = NXT;
      HB = MX + 1;
  }

  for(I = MX; I <= NXT; I++) {
      G = A[I][I];

	  for(J = MX; J <= NXT; J++) {
          //A[I][J] = A[I][J] / G;
          //BZ[I][J] = BZ[I][J] / G;
		  A[J][I] = A[J][I] / G;
          BZ[J][I] = BZ[J][I] / G;
	  }
  }
}

void matrix (double A[10][10], double B[10], double C[10])
{
  int K, L, HM, I, J;
  double D, E, U, UH;

  for (I = MX; I <= NXT; I++) {
	  C[I] = 0;
  }

  for (I = MX; I <= NXT - 1; I++) {
	  UH = A[I][I];
      HM = I;

	  for (K = I+1; K <= NXT; K++) {
		  //if (Abs(UH) < Abs(A[K][I])) {
		  if (Abs(UH) < Abs(A[I][K])) {
              HM = K;
		  }
	  }

	  if (HM != I) {
		  for (K = MX; K <= NXT; K++) {
              //U = A[HM][K];
              //A[HM][K] = A[I][K];
              //A[I][K] = U;
			  U = A[K][HM];
              A[K][HM] = A[K][I];
              A[K][I] = U;
		  }

          U = B[HM];
          B[HM] = B[I];
          B[I] = U;
	  }

	  for (K = I+1; K <= NXT; K++) {
		  if (Abs(A[I][I]) > 0) {
              //D = A[K][I] / A[I][I];
			  D = A[I][K] / A[I][I];
		  }

		  for (L = MX; L <= NXT; L++) {
              //A[K][L] = A[K][L] - D * A[I][L];
			  A[L][K] = A[L][K] - D * A[L][I];
		  }

		  B[K] = B[K] - D * B[I];
	  }
  }

  for (I = NXT; I >= MX; I--) {
      E = 0;

	  for (K = MX; K <= NXT; K++) {
          //E = E + C[K] * A[I][K];
		  E = E + C[K] * A[K][I];
	  }

	  if (Abs(A[I][I]) > 0) {
          C[I] = (B[I] - E) / A[I][I];
	  }
  }
}

void CalcML (void)
{
  int I, J, TRX;
  TRX = 0;

  if (N1 < 1 || N2 > NW || N2 < N1){
    return;
  }
  
  NX = NY + NC;

  if(NX == 0){
    return;
  }

  MX = 0;
  XW[0] = 1;
	 
  if(KBZ == 1){   		  
	  NXT = NX;  
  }

  if(KBZ == 2){  
	  NXT = NX + 1; 
  }  
 
  #ifdef DEBUG
  fprintf(logfile, "\n\nKBZ=%d MX=%d NXT=%d\n", KBZ, MX, NXT);
  #endif

  for(I = 0; I <= NX; I++){
	  BX[I] = 0;
  }

  if(KPZ == 1){
	  BX[0] = 5;
  }

  BX[NX + 1] = 0.01;

  for (I = 0; I < 10; I++) {
      #ifdef DEBUG
	  fprintf(logfile, "Initital BX[%d]=%f\n", I, BX[I]);
      #endif
  }

  do {
	  //fprintf(fp_out, "\nCalML called %d times. ZM=%f", TRX, ZM);
	  nrm24();
	  SQ = 1 - BX[NX + 1];
    
	  for(J = N1; J <= N2; J++) {
		  if(KBZ == 1) {
			  nrm21(J);
			  calML_H1();
			  Q = 1 - P;
    
			  for(I = 0; I <= NX; I++) {
				  PA[I] = Z * XW[I];
			  }
		  } else if(KBZ == 2) {//end if(KBZ == 1)
			  nrm21(J);
			  if(NB == 1) {
				  P = 0;
				  Z = 0;
			  } else {
				  calML_H1();
			  }
			  NRM2BGR();
		  }//end if(KBZ == 2)
		
		  nrm25(J);
	  }//end for(J = N1; J <= N2; ++J) 

	  nrm26();
	  matrix(ZA, ZB, FX);
	  nrm28();
	  TRX += 1;
	  //fprintf(fp_out, "\nAt loops end - CalML called %d times. ZM=%f", TRX, ZM);
  } while (TRX <= 100 && ZM >= 0.0001); //end while(TRX <= 20 && ZM >= .0001)

  matrinv(ZC, BV);
  Chisquare();
  nrm29();

  #ifdef DEBUG
  fprintf(logfile, "\n\nCalcML() finished.\n");
  #endif

  return;
}

void calML_H1(void) {
  nrm22();
  nrm23();
  Y = -5;
  
  if(KPZ == 2) {
	  Y = 0;
  }
  
  for(I = 0; I <= NX; I++) {
	  Y = Y + BX[I] * XW[I];
  }

  pqa();
}

void NRM2BGR(void) {
  Q = 1 - P;
  PW = BX[NX + 1] + P * SQ;

  for(I = 0; I <= NX; I++) {
	  PA[I] = Z * XW[I] * SQ;
  }

  PA[NX + 1] = Q;
}


void CalcDoseResponse(void)
{
	int i, j, IQ, JQ;
	int Plus = 0;
	int student;                       /*Student t indicator*/

	fscanf(fp_in, "%d%lf%lf%d", &student, &TX, &PW, &MV);

	for (i = 1; i <= NY; i++) {
		if (K0[i] == MV)
			Plus = 1;
	}

	for (i = 1; i <= NC; i++) {
		if (K1[i] == MV || K2[i] == MV)
			Plus = 1;
	}

	if (Plus == 0)
		return;

	for (i = 1; i <= NZ; i++)
		XU[i] = 0;                     /* initialize XU[3] */

		
	for (i = 1; i <= NY; i++) {
		for (j = 1; j <= NZ; j++) {
			if (K0[i] == j) {
				XU[j] = j;
			}
		}
	}

	for (i = 1; i <= NC; i++) {
		for (j = 1; j <= NZ; j++) {
			if (K1[i] == j || K2[i] == j) {
				XU[j] = j;
			}
		}
	}

	for (i = 1; i <= NZ; i++) {
		if (i != MV && XU[i] != 0)
			fscanf(fp_in, "%lf", &XU[i]);
	}

	if (PW <= 0)
		PW = 50;

	if (student == 1) {
		Chisquare();
	} else {
		TX = 0;
	}

	fprintf(fp_out, "\n   Estimation of %s", Vrx[MV]);
	fprintf(fp_out, "\n   Response\t = %lf  percent", PW);

	pqb();

    if (KPZ == 1) 
		Y = Y + 5;
		
	AX[0] = BX[0] - Y;

	for (i = 1; i <= NX; i++)
		AX[i] = BX[i];

	XU[MV] = 1;
	
	for (i = 1; i <= NZ; i++) {
		if (i != MV && XU[i] != 0) {
			fprintf(fp_out, "\n   %s\t = %lf", Vrx[i], XU[i]);
			nrm41(i);
		}
	}

	nrm22();
	nrm23();

	XW[0] = 1;
	AW = 0;
	BW = 0;
	CW = 0;
	XZ = 0;
	Cy = 0;

	for (i = 0; i <= NX; i++) {
		IQ = 0;

		if (i > 0 && i < NY + 1)
			IQ = K0[i];

		MQ = i;
		nrm42();

		if (IQ == MV && HQ == 1) {
			Cy += AX[i];
			goto nrm4c;
		}

		if (HQ == 2) {
			Cy += AX[i] * XW[i];
			goto nrm4c;
		}

		XZ = XZ - AX[i] * XW[i];

nrm4c:
		for (j = 0; j <= NX; j++) {
			JQ = 0;

			if (j > 0 && j < NY + 1)
				JQ = K0[j];

			if (HQ == 2) 
				goto nrm4d;

			MQ = j;
			nrm42();

nrm4d:
			CZ = (AX[i] * AX[j] - pow(TX, 2) * BV[j][i] * CX) * XW[i] * XW[j];

			if (HQ == 1) {
				if (IQ != MV && JQ != MV)
					CW += CZ;

				if (IQ == MV && JQ != MV)
					BW += CZ;

				if (JQ == MV && IQ != MV)
					BW += CZ;

				if (IQ == MV && JQ == MV)
					AW += CZ;
			}

			if (HQ == 2) {
				if (IQ != MV && JQ != MV && i != j)
					BW += CZ;
				else
					AW += CZ;
			}
		}
	}

	XZ = XZ / Cy;
	fprintf(fp_out, "\n\n   Estimated %s  %lf percent = ", Vrx[MV], PW);
	nrm43();

	if (student == 1)
		fprintf(fp_out, "   Deviate Corresponding to Confidence Level of Interest = %f\n", TX);

	if (DF == 0) {
		char Msg1[] = "number of degrees of freedom = ";
		char Msg2[] = "Confidence limits cannot be calculated!";
		fprintf(stderr,"\n%s%d\n%s\n", Msg1, DF, Msg2);
		fprintf(fp_out,"\n   %s%d\n%s\n", Msg1, DF, Msg2);
	
		return;
	}

	DC = pow(BW, 2) - 4 * AW * CW;

	if (DC <= 0) {
		char Msg1[] = "95% confidence limits cannot be calculated!";
		char Msg2[] = "Try a smaller Student t or standard normal";
        char Msg3[] = "deviate with smaller confidence probability!";
		fprintf(stderr,"\n%s\n%s\n%s\n", Msg1, Msg2, Msg3);
		fprintf(fp_out,"\n   %s\n%s\n%s\n", Msg1, Msg2, Msg3);
        
		return;
	}

	DW = sqrt(DC);
	XZ = (-BW - DW) / (2 * AW);
	fprintf(fp_out,"   Lower limit %s  %lf percent  = ", Vrx[MV], PW);
	nrm43();

	XZ = (-BW + DW) / (2 * AW);
	fprintf(fp_out,"   Upper limit %s  %lf percent  = ", Vrx[MV], PW);
	nrm43();
}

void CalcResponseDose(void)
{
	int I, J;
	double VY, SY, Y1;
	int student;                       /*Student t indicator*/

	if (NX == 0)
		return;

	fscanf(fp_in, "%d%lf", &student, &TX);

	for (I = 1; I <= NZ; I++)
		XU[I] = 0;                     /* initialize XU[3] */

	for (I = 1; I <= NY; I++) {
		for (J = 1; J <= NZ; J++) {
			if (K0[I] == J) {
				XU[J] = J;
			}
		}
	}

	for (I = 1; I <= NC; I++) {
		for (J = 1; J <= NZ; J++) {
			if (K1[I] == J || K2[I] == J) {
				XU[J] = J;
			}
		}
	}

	for (I = 1; I <= NZ; I++) {
		if (XU[I] != 0)
			fscanf(fp_in, "%lf", &XU[I]);
	}

	if (student == 1) {
		Chisquare();
	} else {
		TX = 0;
	}

	fprintf(fp_out, "\n   Estimation of response");

	for (I = 1; I <= NZ; I++) {
		if (XU[I] != 0)
			fprintf(fp_out, "\n   %s\t = %lf", Vrx[I], XU[I]);

		nrm41(I);
	}

	nrm22();
	nrm23();

	XW[0] = 1;
	VY = 0;
	Y = 0;

	if (KPZ == 1)
		Y = -5;

	for (I = 0; I <= NX; I++) {
		for (J = 0; J <= NX; J++)
			VY += XW[I] * XW[J] * BV[J][I] * CX;
		Y += BX[I] * XW[I];
	}

	SY = sqrt(VY);
	Y1 = Y;

	pqa();
	fprintf(fp_out, "\n\n   Response   = ");
	nrm44();
		
	if (student == 1)
		fprintf(fp_out, "   Deviate Corresponding to Confidence Level of Interest = %f\n", TX);

	if (DF == 0) {
		fprintf(fp_out, "   Degrees of freedom = 0");
		return;
	}

	Y = Y1 - TX * SY;

	pqa();
	fprintf(fp_out, "   LL-response   = ");
	nrm44();

	Y = Y1 + TX * SY;

	pqa();
	fprintf(fp_out, "   UL-response   = ");
	nrm44();

	return;
}

void CalcResponseGraph(void)
{
	double Wdx, Dx, DX1, DX2;
	double Px1, Px2, Px3;
	double SV[NZ];
	int student;                       /*Student t indicator*/
	int index;

	fscanf(fp_in, "%d%lf%d%lf%lf", &student, &TX, &MV, &DX1, &DX2);

	if (student == 0)
		TX = 0;
	
	for (I = 1; I <= NZ; I++)
		SV[I] = 0;                     /* initialize XU[3] */

	for (I = 1; I <= NY; I++) {
		for (J = 1; J <= NZ; J++) {
			if (K0[I] == J) {
				SV[J] = J;
			}
		}
	}

	for (I = 1; I <= NC; I++) {
		for (J = 1; J <= NZ; J++) {
			if (K1[I] == J || K2[I] == J) {
				SV[J] = J;
			}
		}
	}

	for (I = 1; I <= NZ; I++) {
		if (SV[I] != 0 && I != MV) {
			fscanf(fp_in, "%d%lf", &index, &SV[I]);
		}
	}

	Wdx = DX2 - DX1;

	if (Wdx <= 0)
		return;

	if (DX1 = 0)
		DX1 = Wdx / 100;

	Graphic2(DX1, &Px1, &Px2, &Px3, SV);

	for (Dx = DX1; Dx <= DX2; Dx += Wdx/100) {
		Graphic2(Dx, &Px1, &Px2, &Px3, SV);
    
		fprintf(response, "   %f", Dx);
		fprintf(response, "   %f", Px1);
		fprintf(response, "   %f", Px2);
		fprintf(response, "   %f\n", Px3);
	}
}

void Graphic2(double Dx,  double *Px1, double *Px2, double *Px3, double SV[])
{
	int I, J;
	double SY, VY, Y1;

	for (I = 1; I <= NZ; I++) {
		if (MV == I)
			XU[I] = Dx;
		else
			XU[I] = SV[I];

		nrm41(I);
	}

	nrm22();
	nrm23();

	XW[0] = 1;
	VY = 0;
	Y = 0;

	if (KPZ == 1)
		Y = -5;

	for (I = 0; I <= NX; I++) {
		for (J = 0; J <= NX; J++ )
			VY += XW[I] * XW[J] * BV[J][I] * CX;

		Y += BX[I] * XW[I];
	}

	SY = sqrt(VY);
	Y1 = Y;


	pqa();

	*Px1 = P;

	if (DF == 0)
		return;

	Y = Y1 - TX * SY;

	pqa();

	*Px2 = P;

	Y = Y1 + TX * SY;

	pqa();

	*Px3 = P;
}

void CalcRatio(void) 
{
	int H, I, J;
	double RL, RV;
	int student;                       /* Student t indicator */
	char Temp[14][40];                 /* Temporary column names array */

	NX = NY + NC;
	fscanf(fp_in, "%d%lf%d", &student, &TX, &TR);

	if (student == 0)
		TX = 0;

	if (NX < 2 || TR != 2)
		return;

	fscanf(fp_in, "%d%d", &NRC1, &NRC2);

	for (I = 1; I <= NY; I++) {
		strcpy(Temp[I], Vrx[K0[I]]);
	}

	for (I = NY + 1; I <= NX; I++) {
		strncat(Vrx[K1[I-NY]], " ; ", 3);
		int leng = strlen(Vrx[K2[I-NY]]);
		strncat(Vrx[K1[I-NY]], Vrx[K2[I-NY]], leng);
		strcpy(Temp[I], Vrx[K1[I-NY]]);
	}
		
	if (student == 1) {
		Chisquare();
	} else {
		TX = 0;
	}
	
	fprintf(fp_out, "\n   Estimation of ratio between regression coefficients");
	fprintf(fp_out, "\n   Ratio between regression coefficients\n   %s and %s", Temp[NRC1], Temp[NRC2]);

	I = NRC1;
	J = NRC2;

	if (I > NZ || J > NZ)
		return;

	RL = BX[I] / BX[J];
	RV = BV[I][I] / pow(BX[I], 2) + BV[J][J] / pow(BX[J], 2) - 2 * BV[J][I] / BX[I] / BX[J];
	RV = RV * CX * pow(RL, 2);

	if (student == 1)
		fprintf(fp_out, "\n\n   Deviate Corresponding to Confidence Level of Interest = %f", TX);

	fprintf(fp_out, "\n\n   Ratio      =    %f", RL);
	fprintf(fp_out, "\n\n   Confidence limits");
	fprintf(fp_out, "\n     %f     %f", (RL - TX * sqrt(RV)), (RL + TX * sqrt(RV)));

	return;
}

void nrm41 (int IM)
{
	if (NT[IM] == 1) {
		if (XU[IM] > 0)
			XU[IM] = log(XU[IM]);

		return;
	}

	if (NT[IM] == 2)
		if(XU[IM]>0)
			XU[IM]=1/XU[IM];

	return;
}


void nrm42 (void)
{
  HQ = 1;
  if (MQ > NY) {
	  if (K1[MQ-NY] == MV || K2[MQ-NY] == MV)
		  HQ = 2;
  }
}


void nrm43 (void)
{
	if (NT[MV] == 1) {
		fprintf(fp_out,"%9.3e\n", Exp(XZ));
		return;
	}
	  
	if (NT[MV] == 2) {
		fprintf(fp_out,"%9.3e\n", 1/XZ);
		return;
	}
	  
	if (NT[MV] == 3)
		fprintf(fp_out,"%9.3e\n", XZ);
  
	return;
}


void nrm44 (void)
{
  fprintf(fp_out, "%8.2e percent\n", P*100);
}


void Chisquare (void)
{
  double XA, CC, CP, CQ, QP, QZ, S, ZL, ZI;
  int I;

  XA = sqrt(CH);
  
  if(DF % 2 == 0)
      goto CHI1;

  Y = 0;
  if (setjmp(GosubStack[GosubNdx++])==0) goto CHI9;

  if(DF == 1)
      goto CHI2;

  CQ = -log(XA);

  for(I = 1; I <= DF - 2; I += 2) {
      CQ = CQ + 2 * log(XA) - log(I);
      CP = CQ + ZL;
      CC = Exp(CP);
      Y += CC;
  }

CHI2:
  QZ=2*(Z*QP+Y);
  if (setjmp(GosubStack[GosubNdx++])==0) goto CHI8;
#ifdef DEBUG
  fprintf(logfile, "\njumped back to CHI2, QZ=%f QP=%f Y=%f\n\n", QZ, QP, Y);
#endif
  return;

CHI1:
  if (setjmp(GosubStack[GosubNdx++])==0) goto CHI9;
  Y=Z;
  if(DF==2)
    {
      goto CHI3;
    }
  CQ=0;
  for(I=2; I<=DF-2; I+=2)
    {
      CQ=CQ+2*log(XA)-log(I);
      CP=CQ+ZL;
      CC=Exp(CP);
      Y+=CC;
    }

CHI3:
  QZ = ZI * Y;
  if (setjmp(GosubStack[GosubNdx++])==0) goto CHI8;
#ifdef DEBUG
  fprintf(logfile, "\njumped back to CHI3, QZ=%f ZI=%f Y=%f\n\n", QZ, ZI, Y);
#endif
  return;

CHI8:
#ifdef DEBUG
  fprintf(logfile, "\njumped to CHI8, QZ=%f\n\n", QZ);
#endif
  if(QZ>.05 && chisqr_skip == 0)
      Warn1(QZ);
  
  if(QZ<.05 && chisqr_skip == 0)
      Warn2(QZ);

  chisqr_skip = 0;
  longjmp(GosubStack[--GosubNdx],1);


CHI9:
  S = 1 / (1 + .2316419 * Abs(XA));
  ZI = sqrt(8 * atan(1));
  ZL = -(pow(XA,2)) / 2 - log(ZI);
  Z = Exp(ZL);
  QP=.31938153*S-.356563782*pow(S,2)+1.781477937*pow(S,3)-1.821255978*pow(S,4)+1.330274429*pow(S,5);
#ifdef DEBUG
  fprintf(logfile, "\njumped to CHI9, S=%f ZI=%f ZL=%f Z=%f QP=%f", S, ZI, ZL, Z, QP);
#endif
  longjmp(GosubStack[--GosubNdx],1);
}



void Warn1 (double QZ)
{
  CX = 1;

  fprintf(fp_out, "\n\n   Probability of correct model (p-value) is %f\n", QZ);
  fprintf(fp_out,"%s\n","   The prediction of the model is sufficient. Use for estimation of the");
  fprintf(fp_out,"%s\n\n","   95% confidence limits the Standard Normal Deviate");
  fprintf(fp_out,"%s\n\n","   No correction for variances required!");
}


void Warn2 (double QZ)
{
  CX=CH/DF;
  
  fprintf(fp_out, "\n\n   Probability of correct model(p-value) is %f\n", QZ);
  fprintf(fp_out,"%s\n","   The prediction of the model is not sufficient. Use for estimation of the");
  fprintf(fp_out,"%s% d%s\n\n","   95% confidence limits Student t with ",(int)DF," degrees of freedom");
  fprintf(fp_out,"%s%3.3f\n\n","   Correction for variances Chi-Squares/Degrees of Freedom = ",CX);
}


void WriDat (void) {
  int I, J, TM, TX;
  
  NX = NY + NC;

  if( NX == 0) {
      return;
  }

  for(I = 1; I <= NY; I++) {
    WF[I] = K0[I];
  }

  for(I = 1; I <= NC; I++) {
    WF[NY + I] = K1[I];
    WF[NY + NC + I] = K2[I];
  }

  TW = NY + 2 * NC;
  //fprintf(fp_out, "\nInitial TW=%d\n", TW);

  for(I = TW; I >= 2; I--) {
	  for(J = 1; J <= I; J++) {
		  if(WF[I] < WF[J]) {
            TM = WF[I];
            WF[I] = WF[J];
            WF[J] = TM;
          }
      }
  }

  for (I = 1; I <= TW; I++) {
   #ifdef DEBUG
   fprintf(logfile, "WF[%d]=%d ", I, WF[I]);
   #endif
  }

  #ifdef DEBUG
  fprintf(logfile, "\n\n");
  #endif

  I = 1;

  while(I < TW) {
    I += 1;

	if(WF[I] == WF[I - 1]) {
        TW -= 1;
		for(J = I; J <= TW; J++) {
            WF[J] = WF[J + 1];
			#ifdef DEBUG
            fprintf(logfile, "WF[%d]=%d ", J, WF[J]);
            #endif
        }

		#ifdef DEBUG
        fprintf(logfile, "\nTW=%d\n", TW);
        #endif
     }
  }

  if(FlEr == TRUE) {
    return;
  }

time( &ltime );


fprintf(fp_out, "\n");

for(I = 1; I <= TW; I++){
    fprintf(fp_out, "%15s", Vrx[WF[I]]);
}

fprintf(fp_out,"%15s%15s\n\n", Vrx[NZ+1], Vrx[NZ+2]);

// Print out the raw data
TX = 0;

for(I = N1; I <= N2; I++) {

    for(J = 1; J <= TW; J++) {
		fprintf(fp_out,"%15.2f",X[I][WF[J]]);
    }

    fprintf(fp_out,"%14d.",N[I]);
    fprintf(fp_out,"%14d.\n",LR[I]);
    TX += 1;

    if(TX % 5 == 0)
        fprintf(fp_out,"\n");
} 
// End of printing the raw data

fprintf(fp_out,"\n   %s% d%s% d\n","Selection of observations from number ",(int)N1," through ",(int)N2);
fprintf(fp_out,"\n   %s\n","Transformation of input parameters");

for(I=1; I<=NZ; I+=1)
  {
    fprintf(fp_out,"   %-15s ", Vrx[I]);
    while(1)
    {
    if(NT[I]==1)
      {
        fprintf(fp_out,"   %s\n"," is transformed logaritmically!");
        break;
      }
    if(NT[I]==2)
      {
        fprintf(fp_out,"   %s\n"," is transformed reciprocally!");
        break;
      }
    if(NT[I]==3)
      {
        fprintf(fp_out,"   %s\n"," is not transformed at all!");
      }
    break;
    }
}
fprintf(fp_out,"\n");
while(1)
{
if(KPZ==1)
  {
    fprintf(fp_out,"   %s","Probit link used ");
    break;
  }
if(KPZ==2)
  {
    fprintf(fp_out,"   %s","Logit link used ");
  }
break;
}
while(1)
{
if(KBZ==1)
  {
    fprintf(fp_out,"%s\n","without background response correction!");
    break;
  }
if(KBZ==2)
  {
    fprintf(fp_out,"%s\n","with background response correction!");
  }
break;
}
fprintf(fp_out,"\n");
for(I=1; I<=NY; I+=1)
  {
	  fprintf(fp_out,"   %s% d%s%s%s\n","Variable ",(int)I,"  =  ",((NT[I]==3) ? "" : "Transformed "),Vrx[K0[I]]);
  }
if(NC>0)
  {
    for(I=1; I<=NC; I+=1)
      {
        fprintf(fp_out,"   %s% d%s%s%s%s%s%s\n","Variable ",(int)NY+I,"  =  Product of ",((NT[I]==3) ? "" : "transformed "), Vrx[K1[I]]," and ",((NT[I]==3) ? "" : "transformed "),Vrx[K2[I]]);
      }
  }
	
fprintf(fp_out,"\n\n   Chi-Square          =  %-6.2f", CH);
fprintf(fp_out,"\n   Degrees of Freedom  =  %-5d\n\n", DF);

for(I=MX; I<=NXT; I+=1)
  {
    fprintf(fp_out,"   %s%d%s","B",(int)I," = ");
    fprintf(fp_out,"%9.3e\t",BX[I]);
    if(DF>0)
      {
        fprintf(fp_out,"   %s%d%s","Student t  for B",(int)I," = ");
        fprintf(fp_out,"%2.2f\n",BX[I]/sqrt(BV[I][I]*CX));
      }
    else
      {
        fprintf(fp_out,"\n");
      }
  }
if(DF==0)
  {
    return;
  }
fprintf(fp_out,"\n");
for(I=MX; I<=NXT; I+=1)
  {
    for(J=I; J<=NXT; J+=1)
      {
        while(1)
        {
        if(I==J)
          {
            fprintf(fp_out,"   %s%d%d%s","  variance  B",(int)I,(int)J," = ");
            break;
          }
        if(I!=J)
          {
            fprintf(fp_out,"   %s%d%d%s","covariance  B",(int)I,(int)J," = ");
          }
        break;
        }
      fprintf(fp_out,"%9.3e\n",BV[I][J]*CX);
    }
}

return;

FAULT1:;
CalcErr();
}


void CalcErr (void)
{
  time( &ltime );
  static char Msg[2048];
  memset(&Msg,0,sizeof(Msg));
  NY=0;
  NC=0;
  NX=0;
  if(FP2)
   {
     fflush(FP2);
     fclose(FP2);
   }
  FlEr=TRUE;
  sprintf(BCX_STR,"%s","The combination of data\n");
  strcpy(Msg,BCX_STR);
  sprintf(Msg,"%s%s",Msg,"and design for mathematical\n");
  sprintf(Msg,"%s%s",Msg,"analysis produced an error\n");
  sprintf(Msg,"%s%s",Msg,"in the calculating procedures!\n");
  sprintf(Msg,"%s%s",Msg,"Please, revise your design\n");
  sprintf(Msg,"%s%s",Msg,"for mathematical analysis!");
  MessageBox (GetActiveWindow(),Msg,"",0);

  if((FP2=fopen("DoseResp.Log","a"))==0){
	fprintf(stderr,"Can't open file %s\n","DoseResp.Log");
	exit(1);
  }

  fprintf(FP2,"\n");
  fprintf(FP2,"%s%s\n","Filename = ",Fina);
  fprintf(FP2,"%s\n",ctime(&ltime));
  fprintf(FP2,"%s\n","An error occurred in the calculation!");
  fprintf(FP2,"\n");
  if(FP2)
   {
     fflush(FP2);
     fclose(FP2);
   }
}

// ************************************************************
//    READ_OBSDATA--used to read data in a matrix(row, col).
// ************************************************************

int READ_OBSDATA (int r, int c, double **matrix)
{
  int Nmiss;			/*number of records with missing values */
  int i, j, n, m;		/*count and iteration control variables */
  double dvalue;			/*temp variable */
  int ivalue;

  Nmiss = 0;
  for (i = 1; i <= r; i++) {
      n = i - Nmiss;
      m = 0;

      for (j = 1; j <= c + 2; j++) {
		  if (j <= c)
			  fscanf(fp_in, "%lf", &dvalue);
		  else
			  fscanf(fp_in, "%d", &ivalue);

		if ((j <=c && dvalue != MISSING) || (j > c && ivalue != MISSING)) {
			if(j == NZ + 1) {
				N[n] = ivalue;
			#ifdef DEBUG
				fprintf(logfile, "\nN[%d] = %d", n, N[n]);
            #endif
			}

			if(j == NZ + 2) {
				LR[n] = ivalue;
			#ifdef DEBUG
				fprintf(logfile, "\nLR[%d] = %d", n, LR[n]);
            #endif
			}

			if (j <= c) {
				matrix[n][j] = dvalue;
			#ifdef DEBUG
				fprintf(logfile, "\nData matrix[%d,%d] = %lf", n, j, matrix[n][j]);
            #endif
			}
		}
	    else
			m++;
	  }

      if (m != 0)
		Nmiss++;
  }

  return Nmiss;
}


void Get_File_Stem(char *argv, char *filestem)
{
	int i = 0;

	for (; i <= 127; i++) { /* TODO: change to dynamically calculate the length of argv */
		filestem[i] = argv[i];

		if (argv[i] == '.')
			break;
	}

	infilestem[i] = 0; /* ends the file stem properly*/
}

void Derive_File_Name(char *stem, char *newfile, const char *ext)
{
	int len = strlen(stem);
	strcpy(newfile, stem);
	newfile[len] = '.';
	strncat(newfile, ext, strlen(ext));
	newfile[len + strlen(ext) + 1] = 0; /* to end the new file name properly */
}
