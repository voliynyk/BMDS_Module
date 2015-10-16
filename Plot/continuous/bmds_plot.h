/* Copyright 2014 - Louis Olszyk on behalf of Lockheed Martin for the U.S. EPA
 * All rights reserved.
 */

#ifndef __BMDS_PLOT_H__
#define __BMDS_PLOT_H__

#ifndef TERM_X11
/*#define MYTERM "set terminal windows\n"*/
#define MYTERM "set terminal emf color dashed\n"
#else
#define MYTERM "set terminal x11\n"
#endif

#define BMDS_NAME_SIZE 64

/*****************************
typedef struct {
  eModel_t Number;
  char *Name;
} model_t;
const char zModels = {{evHill, "Hill"},
			 {evExponential, "Exponential"},
			 {evLinear, "Linear"},
			 {evPolynomial, "Polynomial"},
			 {evPower, "Power"}
};
*****************************/
#define BMDS_MODEL_NAME(mNum) pcModels[(mNum)]

typedef enum {
    evAbs=0, evStd=1, evRel=2, evPoint=3, evExtra=4
} eRisk_t;

/* Model identifiers */
typedef enum {
    evHill=0,
    evExponential=1,
    evLinear=5,
    evPolynomial=6,
    evPower=7
} eModel_t;
/* Macros for indexing into list of func strings */
#define iF_HILL  0
#define iF_EXPO2 1
#define iF_EXPO3 2
#define iF_EXPO4 3
#define iF_EXPO5 4
#define iF_LIN   5
#define iF_POLY  6
#define iF_POW   7
/* Offset for exponential sub-model, which starts at "2" */
#define iF_EXPO_OFFSET 2
/* a and b are used for the variance parameters */

/* Store parameter name/value pairs */
typedef struct {
  char Name[BMDS_NAME_SIZE];
  int iValue;
} bmds_i_parm_t;

typedef struct {
  char Name[BMDS_NAME_SIZE];
  double dValue;
} bmds_d_parm_t;

#define BMDS_MAX_I_PARMS 10 /* Arbitrary max size for array of integer parms */
#define BMDS_MAX_D_PARMS 10 /* Arbitrary max size for array of double parms */

/* Indices for integral parameters */
#define iFUNC     0 /* Model function index */
#define iBMDFLAG  1
#define iNOBS     2 /* # observations */
#define iNPARMS   3 /* # model parameters */
#define iBMRTYPE  4
#define iSIGN     5
#define iLOGNORM  6

/* Indices for double parameters */
#define iCONF   0
#define iBMRF   1
#define iALPHA  2 /* or lnalpha */
#define iRHO    3
#define iP1     4
#define iP2     5
#define iP3     6
#define iP4     7

/* Common strings used to create the plots */
static const char acSetErrorBarStyle[] = "set style line 16 linecolor rgb \"forest-green\" pt 12\n";
static const char acSetBMDLineStyle[] = "set style line 6 linecolor rgb \"black\"\n";
static const char acErrorBars[] = "with errorbars ls 16";
static const char acLines02[] = "with lines ls 6";
static const char *pcRiskTypes[] = {"Abs. Dev.", "Std. Dev.", "Rel. Dev.",
			    "Point", "Extra"
};

/* Model names - For simplicity treat the exponential sub-functions
 * as separate models.
 */
static const char *pcModels[] = {
  "Hill",
  "Exponential 2", "Exponential 3", "Exponential 4", "Exponential 5",
  "Linear",
  "Polynomial",
  "Power"
};

/* Model function strings for plotting */
/* - Note that Expo 2&3 have a substitution parameter.  */
static const char *pcModelFuncs[] = {
  "HILL",
  "c * exp(%d * d * x)",
  "c * exp(%d * (d * x)**f)",
  "c * (e - ((e - 1) * exp(-d * x)))",
  "c * (e - ((e - 1) * exp(-(d * x)**f)))",
  "LINEAR",
  "POLYNOMIAL",
  "c+(x>0?( d*x**e ) : 0)"
};

/******** Function prototypes *********/
void plot_continuous(const eModel_t eModel, int argc, char *argv[]);

#endif
/* end __BMDS_PLOT_H__ */
