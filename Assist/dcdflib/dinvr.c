#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cdflib.h"

/* DEFINE DINVR */
static void E0000(int IENTRY,int *status,double *x,double *fx,
		  unsigned long *qleft,unsigned long *qhi,double *zabsst,
		  double *zabsto,double *zbig,double *zrelst,
		  double *zrelto,double *zsmall,double *zstpmu)
{
#define qxmon(zx,zy,zz) (int)((zx) <= (zy) && (zy) <= (zz))
static double absstp,abstol,big,fbig,fsmall,relstp,reltol,small,step,stpmul,xhi,
    xlb,xlo,xsave,xub,yy;
static int i99999;
static unsigned long qbdd,qcond,qdum1,qdum2,qincr,qlim,qok,qup;
    switch(IENTRY){case 0: goto DINVR; case 1: goto DSTINV;}
DINVR:
    if(*status > 0) goto S310;
    qcond = !qxmon(small,*x,big);
    if(qcond) ftnstop(" SMALL, X, BIG not monotone in INVR");
    xsave = *x;
/*
     See that SMALL and BIG bound the zero and set QINCR
*/
    *x = small;
/*
     GET-FUNCTION-VALUE
*/
    i99999 = 1;
    goto S300;
S10:
    fsmall = *fx;
    *x = big;
/*
     GET-FUNCTION-VALUE
*/
    i99999 = 2;
    goto S300;
S20:
    fbig = *fx;
    qincr = fbig > fsmall;
    if(!qincr) goto S50;
    if(fsmall <= 0.0e0) goto S30;
    *status = -1;
    *qleft = *qhi = 1;
    return;
S30:
    if(fbig >= 0.0e0) goto S40;
    *status = -1;
    *qleft = *qhi = 0;
    return;
S40:
    goto S80;
S50:
    if(fsmall >= 0.0e0) goto S60;
    *status = -1;
    *qleft = 1;
    *qhi = 0;
    return;
S60:
    if(fbig <= 0.0e0) goto S70;
    *status = -1;
    *qleft = 0;
    *qhi = 1;
    return;
S80:
S70:
    *x = xsave;
    step = fifdmax1(absstp,relstp*fabs(*x));
/*
      YY = F(X) - Y
     GET-FUNCTION-VALUE
*/
    i99999 = 3;
    goto S300;
S90:
    yy = *fx;
    if(!(yy == 0.0e0)) goto S100;
    *status = 0;
    qok = 1;
    return;
S100:
    qup = ((qincr && (yy < 0.0e0)) || (!qincr && (yy > 0.0e0)));
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     HANDLE CASE IN WHICH WE MUST STEP HIGHER
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
    if(!qup) goto S170;
    xlb = xsave;
    xub = fifdmin1(xlb+step,big);
    goto S120;
S110:
    if(qcond) goto S150;
S120:
/*
      YY = F(XUB) - Y
*/
    *x = xub;
/*
     GET-FUNCTION-VALUE
*/
    i99999 = 4;
    goto S300;
S130:
    yy = *fx;
    qbdd = ((qincr && (yy >= 0.0e0)) || (!qincr && (yy <= 0.0e0)));
    qlim = xub >= big;
    qcond = qbdd || qlim;
    if(qcond) goto S140;
    step = stpmul*step;
    xlb = xub;
    xub = fifdmin1(xlb+step,big);
S140:
    goto S110;
S150:
    if(!(qlim && !qbdd)) goto S160;
    *status = -1;
    *qleft = 0;
    *qhi = !qincr;
    *x = big;
    return;
S160:
    goto S240;
S170:
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     HANDLE CASE IN WHICH WE MUST STEP LOWER
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
    xub = xsave;
    xlb = fifdmax1(xub-step,small);
    goto S190;
S180:
    if(qcond) goto S220;
S190:
/*
      YY = F(XLB) - Y
*/
    *x = xlb;
/*
     GET-FUNCTION-VALUE
*/
    i99999 = 5;
    goto S300;
S200:
    yy = *fx;
    qbdd =((qincr && (yy <= 0.0e0)) || (!qincr && (yy >= 0.0e0)));
    qlim = xlb <= small;
    qcond = qbdd || qlim;
    if(qcond) goto S210;
    step = stpmul*step;
    xub = xlb;
    xlb = fifdmax1(xub-step,small);
S210:
    goto S180;
S220:
    if(!(qlim && !qbdd)) goto S230;
    *status = -1;
    *qleft = 1;
    *qhi = qincr;
    *x = small;
    return;
S240:
S230:
    dstzr(&xlb,&xub,&abstol,&reltol);
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
    *status = 0;
    goto S260;
S250:
    if(!(*status == 1)) goto S290;
S260:
    dzror(status,x,fx,&xlo,&xhi,&qdum1,&qdum2);
    if(!(*status == 1)) goto S280;
/*
     GET-FUNCTION-VALUE
*/
    i99999 = 6;
    goto S300;
S280:
S270:
    goto S250;
S290:
    *x = xlo;
    *status = 0;
    return;
DSTINV:
    small = *zsmall;
    big = *zbig;
    absstp = *zabsst;
    relstp = *zrelst;
    stpmul = *zstpmu;
    abstol = *zabsto;
    reltol = *zrelto;
    return;
S300:
/*
     TO GET-FUNCTION-VALUE
*/
    *status = 1;
    return;
S310:
    switch((int)i99999){case 1: goto S10;case 2: goto S20;case 3: goto S90;case 
      4: goto S130;case 5: goto S200;case 6: goto S270;default: break;}
#undef qxmon
}


void dinvr(int *status,double *x,double *fx,
	   unsigned long *qleft,unsigned long *qhi)
/*
**********************************************************************
 
     void dinvr(int *status,double *x,double *fx,
           unsigned long *qleft,unsigned long *qhi)

          Double precision
          bounds the zero of the function and invokes zror
                    Reverse Communication
 
 
                              Function
 
 
     Bounds the    function  and  invokes  ZROR   to perform the   zero
     finding.  STINVR  must  have   been  called  before this   routine
     in order to set its parameters.
 
 
                              Arguments
 
 
     STATUS <--> At the beginning of a zero finding problem, STATUS
                 should be set to 0 and INVR invoked.  (The value
                 of parameters other than X will be ignored on this cal
 
                 When INVR needs the function evaluated, it will set
                 STATUS to 1 and return.  The value of the function
                 should be set in FX and INVR again called without
                 changing any of its other parameters.
 
                 When INVR has finished without error, it will return
                 with STATUS 0.  In that case X is approximately a root
                 of F(X).
 
                 If INVR cannot bound the function, it returns status
                 -1 and sets QLEFT and QHI.
                         INTEGER STATUS
 
     X <-- The value of X at which F(X) is to be evaluated.
                         DOUBLE PRECISION X
 
     FX --> The value of F(X) calculated when INVR returns with
            STATUS = 1.
                         DOUBLE PRECISION FX
 
     QLEFT <-- Defined only if QMFINV returns .FALSE.  In that
          case it is .TRUE. If the stepping search terminated
          unsucessfully at SMALL.  If it is .FALSE. the search
          terminated unsucessfully at BIG.
                    QLEFT is LOGICAL
 
     QHI <-- Defined only if QMFINV returns .FALSE.  In that
          case it is .TRUE. if F(X) .GT. Y at the termination
          of the search and .FALSE. if F(X) .LT. Y at the
          termination of the search.
                    QHI is LOGICAL
 
**********************************************************************
*/
{

  E0000(0,status,x,fx,qleft,qhi,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
}



void dstinv(double *zsmall,double *zbig,double *zabsst,
	    double *zrelst,double *zstpmu,double *zabsto,
	    double *zrelto)
/*
**********************************************************************
      void dstinv(double *zsmall,double *zbig,double *zabsst,
            double *zrelst,double *zstpmu,double *zabsto,
            double *zrelto)

      Double Precision - SeT INverse finder - Reverse Communication
                              Function
     Concise Description - Given a monotone function F finds X
     such that F(X) = Y.  Uses Reverse communication -- see invr.
     This routine sets quantities needed by INVR.
          More Precise Description of INVR -
     F must be a monotone function, the results of QMFINV are
     otherwise undefined.  QINCR must be .TRUE. if F is non-
     decreasing and .FALSE. if F is non-increasing.
     QMFINV will return .TRUE. if and only if F(SMALL) and
     F(BIG) bracket Y, i. e.,
          QINCR is .TRUE. and F(SMALL).LE.Y.LE.F(BIG) or
          QINCR is .FALSE. and F(BIG).LE.Y.LE.F(SMALL)
     if QMFINV returns .TRUE., then the X returned satisfies
     the following condition.  let
               TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
     then if QINCR is .TRUE.,
          F(X-TOL(X)) .LE. Y .LE. F(X+TOL(X))
     and if QINCR is .FALSE.
          F(X-TOL(X)) .GE. Y .GE. F(X+TOL(X))
                              Arguments
     SMALL --> The left endpoint of the interval to be
          searched for a solution.
                    SMALL is DOUBLE PRECISION
     BIG --> The right endpoint of the interval to be
          searched for a solution.
                    BIG is DOUBLE PRECISION
     ABSSTP, RELSTP --> The initial step size in the search
          is MAX(ABSSTP,RELSTP*ABS(X)). See algorithm.
                    ABSSTP is DOUBLE PRECISION
                    RELSTP is DOUBLE PRECISION
     STPMUL --> When a step doesn't bound the zero, the step
                size is multiplied by STPMUL and another step
                taken.  A popular value is 2.0
                    DOUBLE PRECISION STPMUL
     ABSTOL, RELTOL --> Two numbers that determine the accuracy
          of the solution.  See function for a precise definition.
                    ABSTOL is DOUBLE PRECISION
                    RELTOL is DOUBLE PRECISION
                              Method
     Compares F(X) with Y for the input value of X then uses QINCR
     to determine whether to step left or right to bound the
     desired x.  the initial step size is
          MAX(ABSSTP,RELSTP*ABS(S)) for the input value of X.
     Iteratively steps right or left until it bounds X.
     At each step which doesn't bound X, the step size is doubled.
     The routine is careful never to step beyond SMALL or BIG.  If
     it hasn't bounded X at SMALL or BIG, QMFINV returns .FALSE.
     after setting QLEFT and QHI.
     If X is successfully bounded then Algorithm R of the paper
     'Two Efficient Algorithms with Guaranteed Convergence for
     Finding a Zero of a Function' by J. C. P. Bus and
     T. J. Dekker in ACM Transactions on Mathematical
     Software, Volume 1, No. 4 page 330 (DEC. '75) is employed
     to find the zero of the function F(X)-Y. This is routine
     QRZERO.
**********************************************************************
*/
{
    E0000(1,NULL,NULL,NULL,NULL,NULL,zabsst,zabsto,zbig,zrelst,zrelto,zsmall,
    zstpmu);
}
