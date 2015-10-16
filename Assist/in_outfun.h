/****************
*  in_outfun.h
*  Nov. 22 1996
* Modified by: Geoffrey Nonato
* Date Modified: 12/15/2006
* Modification: Added "int print_SE" as the last parameter for function OP_ParmsE, to facilitate print/no print SE for each model.
*****************/
void READ_PARAMETERS(int nparm, double Parms[]);
int READ_OBSDATA5V(int Nobs,double X[],double Yp[],double Yn[],double Ls[],    
                   int Xg[]);
int READ_OBSDATA3V(int Nobs,int tcols,int pdcol,int ndcol,int idcol,
                   double Yp[],double Yn[],double X[]);
void OUTPUT_TEXT (char txt[]);

void OUTPUT_NDTMSPARMS(int nparm,int dgrp,double Parms[],
                       char *tparms[],double **vcv);
void OUTPUT_Init(int nparm, int Spec[],double Parms[],char *tparms[]);

void OUTPUT_DTMS3PARMS(int nparm, int Spec[], int bounded[], double Parms[],
                       char *tparms[],double **vcv, int print_SE);

void OUTPUT_DTMS3ANOVA(char *anatxt[], AnaList anasum[]);

void OUTPUT_DTMSGNFIT(int Nobs,int Spec[],
                      double Yp[],double Yn[],double X[],double ypp[]);
void OUTPUT_BENCHMD(int pdcol, double BMD);

void OUTPUT_DTMS3ANOVAC(char *anatxt[], AnaList anasum[], int type);
void Get_and_OUTPUT_DTMSVCV(int nparm, int Spec[], char *parmtxt[],double **vcv, double **vcv_adj, int *bounded);
void OP_ParmsE(int nparm, int Spec[],double Parms[],char *tparms[],double **vcv,int *bounded, double conf, int print_SE);
void Output_Header(char *version, char *input_name, char *plot_file_name, char *clocktime, char *note);
void do_dmngb_warning( int * );

