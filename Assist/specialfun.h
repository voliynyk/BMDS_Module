/****************
*  specialfun.h
*  Nov. 22 1996
*****************/

double GAMMP(double, double);
double XGAMMAI_A(double, double);
double XGAMMAI(double, double);
double CNORM(double);
double NORM(double);
double RNORM(double);
double CHISQ(double, int);
double QCHISQ(double, int);
double qstudt(double, int);
void PROBABILITY_W();
void PROBABILITY_INRANGE();
void Sort_2_By_Dose(int size, double xxi[], double yyi[]);
void SortByLs (int nobs, int ngrp, int GrpSize[], double GLs[],
	       double GYp[], double GYn[], double GYpp[],
	       double GEp[]);
void Sort_4_By_Dose(int size, double dos[], double tot[], double res[], double litter[]);
void Get_Names(char *input_name, char *fileout_name, char *file002_name, char *plot_file_name);
void path_name(char *input_name);
void path_name2(int argc, char **argv, char *out);
double Slog(double);
double D_Slog(double);
void show_version(char argv[], char version[]);

#define SQR(a)((a) == 0.0 ?  0.0 : (a)*(a))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define FMAX(a,b)((a)>(b) ? (a):(b))
#define FMIN(a,b)((a)<(b) ? (a):(b))
