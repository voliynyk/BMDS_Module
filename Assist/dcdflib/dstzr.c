#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cdflib.h"


/* DEFINE DZROR */
static void E0001(int IENTRY,int *status,double *x,double *fx,
		  double *xlo,double *xhi,unsigned long *qleft,
		  unsigned long *qhi,double *zabstl,double *zreltl,
		  double *zxhi,double *zxlo)
{
#define ftol(zx) (0.5e0*fifdmax1(abstol,reltol*fabs((zx))))
static double a,abstol,b,c,d,fa,fb,fc,fd,fda,fdb,m,mb,p,q,reltol,tol,w,xxhi,xxlo;
static int ext,i99999;
static unsigned long first,qrzero;
    switch(IENTRY){case 0: goto DZROR; case 1: goto DSTZR;}
DZROR:
    if(*status > 0) goto S280;
    *xlo = xxlo;
    *xhi = xxhi;
    b = *x = *xlo;
/*
     GET-FUNCTION-VALUE
*/
    i99999 = 1;
    goto S270;
S10:
    fb = *fx;
    *xlo = *xhi;
    a = *x = *xlo;
/*
     GET-FUNCTION-VALUE
*/
    i99999 = 2;
    goto S270;
S20:
/*
     Check that F(ZXLO) < 0 < F(ZXHI)  or
                F(ZXLO) > 0 > F(ZXHI)
*/
    if(!(fb < 0.0e0)) goto S40;
    if(!(*fx < 0.0e0)) goto S30;
    *status = -1;
    *qleft = *fx < fb;
    *qhi = 0;
    return;
S40:
S30:
    if(!(fb > 0.0e0)) goto S60;
    if(!(*fx > 0.0e0)) goto S50;
    *status = -1;
    *qleft = *fx > fb;
    *qhi = 1;
    return;
S60:
S50:
    fa = *fx;
    first = 1;
S70:
    c = a;
    fc = fa;
    ext = 0;
S80:
    if(!(fabs(fc) < fabs(fb))) goto S100;
    if(!(c != a)) goto S90;
    d = a;
    fd = fa;
S90:
    a = b;
    fa = fb;
    *xlo = c;
    b = *xlo;
    fb = fc;
    c = a;
    fc = fa;
S100:
    tol = ftol(*xlo);
    m = (c+b)*.5e0;
    mb = m-b;
    if(!(fabs(mb) > tol)) goto S240;
    if(!(ext > 3)) goto S110;
    w = mb;
    goto S190;
S110:
    tol = fifdsign(tol,mb);
    p = (b-a)*fb;
    if(!first) goto S120;
    q = fa-fb;
    first = 0;
    goto S130;
S120:
    fdb = (fd-fb)/(d-b);
    fda = (fd-fa)/(d-a);
    p = fda*p;
    q = fdb*fa-fda*fb;
S130:
    if(!(p < 0.0e0)) goto S140;
    p = -p;
    q = -q;
S140:
    if(ext == 3) p *= 2.0e0;
    if(!(p*1.0e0 == 0.0e0 || p <= q*tol)) goto S150;
    w = tol;
    goto S180;
S150:
    if(!(p < mb*q)) goto S160;
    w = p/q;
    goto S170;
S160:
    w = mb;
S190:
S180:
S170:
    d = a;
    fd = fa;
    a = b;
    fa = fb;
    b += w;
    *xlo = b;
    *x = *xlo;
/*
     GET-FUNCTION-VALUE
*/
    i99999 = 3;
    goto S270;
S200:
    fb = *fx;
    if(!(fc*fb >= 0.0e0)) goto S210;
    goto S70;
S210:
    if(!(w == mb)) goto S220;
    ext = 0;
    goto S230;
S220:
    ext += 1;
S230:
    goto S80;
S240:
    *xhi = c;
    qrzero = (((fc >= 0.0e0) && (fb <= 0.0e0)) || ((fc < 0.0e0) && (fb >= 0.0e0)));
    if(!qrzero) goto S250;
    *status = 0;
    goto S260;
S250:
    *status = -1;
S260:
    return;
DSTZR:
    xxlo = *zxlo;
    xxhi = *zxhi;
    abstol = *zabstl;
    reltol = *zreltl;
    return;
S270:
/*
     TO GET-FUNCTION-VALUE
*/
    *status = 1;
    return;
S280:
    switch((int)i99999){case 1: goto S10;case 2: goto S20;case 3: goto S200;
      default: break;}
#undef ftol
}


void dstzr(double *zxlo,double *zxhi,double *zabstl,double *zreltl)
/*
**********************************************************************
     void dstzr(double *zxlo,double *zxhi,double *zabstl,double *zreltl)
     Double precision SeT ZeRo finder - Reverse communication version
                              Function
     Sets quantities needed by ZROR.  The function of ZROR
     and the quantities set is given here.
     Concise Description - Given a function F
     find XLO such that F(XLO) = 0.
          More Precise Description -
     Input condition. F is a double precision function of a single
     double precision argument and XLO and XHI are such that
          F(XLO)*F(XHI)  .LE.  0.0
     If the input condition is met, QRZERO returns .TRUE.
     and output values of XLO and XHI satisfy the following
          F(XLO)*F(XHI)  .LE. 0.
          ABS(F(XLO)  .LE. ABS(F(XHI)
          ABS(XLO-XHI)  .LE. TOL(X)
     where
          TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
     If this algorithm does not find XLO and XHI satisfying
     these conditions then QRZERO returns .FALSE.  This
     implies that the input condition was not met.
                              Arguments
     XLO --> The left endpoint of the interval to be
           searched for a solution.
                    XLO is DOUBLE PRECISION
     XHI --> The right endpoint of the interval to be
           for a solution.
                    XHI is DOUBLE PRECISION
     ABSTOL, RELTOL --> Two numbers that determine the accuracy
                      of the solution.  See function for a
                      precise definition.
                    ABSTOL is DOUBLE PRECISION
                    RELTOL is DOUBLE PRECISION
                              Method
     Algorithm R of the paper 'Two Efficient Algorithms with
     Guaranteed Convergence for Finding a Zero of a Function'
     by J. C. P. Bus and T. J. Dekker in ACM Transactions on
     Mathematical Software, Volume 1, no. 4 page 330
     (Dec. '75) is employed to find the zero of F(X)-Y.
**********************************************************************
*/
{
    E0001(1,NULL,NULL,NULL,NULL,NULL,NULL,NULL,zabstl,zreltl,zxhi,zxlo);
}



void dzror(int *status,double *x,double *fx,double *xlo,
	   double *xhi,unsigned long *qleft,unsigned long *qhi)
/*
**********************************************************************
 
     void dzror(int *status,double *x,double *fx,double *xlo,
           double *xhi,unsigned long *qleft,unsigned long *qhi)

     Double precision ZeRo of a function -- Reverse Communication
 
 
                              Function
 
 
     Performs the zero finding.  STZROR must have been called before
     this routine in order to set its parameters.
 
 
                              Arguments
 
 
     STATUS <--> At the beginning of a zero finding problem, STATUS
                 should be set to 0 and ZROR invoked.  (The value
                 of other parameters will be ignored on this call.)
 
                 When ZROR needs the function evaluated, it will set
                 STATUS to 1 and return.  The value of the function
                 should be set in FX and ZROR again called without
                 changing any of its other parameters.
 
                 When ZROR has finished without error, it will return
                 with STATUS 0.  In that case (XLO,XHI) bound the answe
 
                 If ZROR finds an error (which implies that F(XLO)-Y an
                 F(XHI)-Y have the same sign, it returns STATUS -1.  In
                 this case, XLO and XHI are undefined.
                         INTEGER STATUS
 
     X <-- The value of X at which F(X) is to be evaluated.
                         DOUBLE PRECISION X
 
     FX --> The value of F(X) calculated when ZROR returns with
            STATUS = 1.
                         DOUBLE PRECISION FX
 
     XLO <-- When ZROR returns with STATUS = 0, XLO bounds the
             inverval in X containing the solution below.
                         DOUBLE PRECISION XLO
 
     XHI <-- When ZROR returns with STATUS = 0, XHI bounds the
             inverval in X containing the solution above.
                         DOUBLE PRECISION XHI
 
     QLEFT <-- .TRUE. if the stepping search terminated unsucessfully
                at XLO.  If it is .FALSE. the search terminated
                unsucessfully at XHI.
                    QLEFT is LOGICAL
 
     QHI <-- .TRUE. if F(X) .GT. Y at the termination of the
              search and .FALSE. if F(X) .LT. Y at the
              termination of the search.
                    QHI is LOGICAL
 
**********************************************************************
*/
{
    E0001(0,status,x,fx,xlo,xhi,qleft,qhi,NULL,NULL,NULL,NULL);
}


