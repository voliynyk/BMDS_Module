#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cdflib.h"

void cdftnc(int *which,double *p,double *q,double *t,double *df,
            double *pnonc,int *status,double *bound)
/**********************************************************************
 
   void cdftnc(int *which,double *p,double *q,double *t,double *df,
               double *pnonc,int *status,double *bound)

                Cumulative Distribution Function
                   Non-Central T distribution
 
                                Function
 
      Calculates any one parameter of the noncentral t distribution give
      values for the others.
 
                                Arguments
 
      WHICH --> Integer indicating which  argument
                values is to be calculated from the others.
                Legal range: 1..3
                iwhich = 1 : Calculate P and Q from T,DF,PNONC
                iwhich = 2 : Calculate T from P,Q,DF,PNONC
                iwhich = 3 : Calculate DF from P,Q,T
                iwhich = 4 : Calculate PNONC from P,Q,DF,T
 
         P <--> The integral from -infinity to t of the noncentral t-den
               Input range: (0,1].
 
         Q <--> 1-P.
               Input range: (0, 1].
                P + Q = 1.0.
 
         T <--> Upper limit of integration of the noncentral t-density.
                Input range: ( -infinity, +infinity).
                Search range: [ -1E100, 1E100 ]
 
         DF <--> Degrees of freedom of the noncentral t-distribution.
                 Input range: (0 , +infinity).
                 Search range: [1e-100, 1E10]
 
      PNONC <--> Noncentrality parameter of the noncentral t-distributio
                 Input range: [-infinity , +infinity).
                 Search range: [-1e4, 1E4]
 
      STATUS <-- 0 if calculation completed correctly
                -I if input parameter number I is out of range
                 1 if answer appears to be lower than lowest
                   search bound
                 2 if answer appears to be higher than greatest
                   search bound
                 3 if P + Q .ne. 1
 
      BOUND <-- Undefined if STATUS is 0
 
                Bound exceeded by parameter number I if STATUS
                is negative.
 
                Lower search bound if STATUS is 1.
 
                Upper search bound if STATUS is 2.
 
                                 Method
 
      Upper tail    of  the  cumulative  noncentral t is calculated usin
      formulae  from page 532  of Johnson, Kotz,  Balakrishnan, Coninuou
      Univariate Distributions, Vol 2, 2nd Edition.  Wiley (1995)
 
      Computation of other parameters involve a seach for a value that
      produces  the desired  value  of P.   The search relies  on  the
      monotinicity of P with the other parameter.
 
**********************************************************************/
{
#define tent4 1.0e4
#define tol 1.0e-8
#define atol 1.0e-50
#define zero 1.0e-100
#define one ( 1.0e0 - 1.0e-16 )
#define inf 1.0e100
static double K3 = 0.5e0;
static double K4 = 5.0e0;
static double ccum,cum,fx;
static unsigned long qhi,qleft;
static double T1,T2,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14;
/*
     ..
     .. Executable Statements ..
*/
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 5.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
    if(!(*p < 0.0e0 || *p > one)) goto S60;
    if(!(*p < 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = one;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 3) goto S90;
    if(!(*df <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    *status = -5;
    return;
S90:
S80:
    if(*which == 4) goto S100;
S100:
    if(1 == *which) {
        cumtnc(t,df,pnonc,p,q);
        *status = 0;
    }
    else if(2 == *which) {
        *t = 5.0e0;
        T1 = -inf;
        T2 = inf;
        T5 = atol;
        T6 = tol;
        dstinv(&T1,&T2,&K3,&K3,&K4,&T5,&T6);
        *status = 0;
        dinvr(status,t,&fx,&qleft,&qhi);
S110:
        if(!(*status == 1)) goto S120;
        cumtnc(t,df,pnonc,&cum,&ccum);
        fx = cum - *p;
        dinvr(status,t,&fx,&qleft,&qhi);
        goto S110;
S120:
        if(!(*status == -1)) goto S150;
        if(!qleft) goto S130;
        *status = 1;
        *bound = -inf;
        goto S140;
S130:
        *status = 2;
        *bound = inf;
S150:
S140:
        ;
    }
    else if(3 == *which) {
        *df = 5.0e0;
        T7 = zero;
        T8 = tent4;
        T9 = atol;
        T10 = tol;
        dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
        *status = 0;
        dinvr(status,df,&fx,&qleft,&qhi);
S160:
        if(!(*status == 1)) goto S170;
        cumtnc(t,df,pnonc,&cum,&ccum);
        fx = cum - *p;
        dinvr(status,df,&fx,&qleft,&qhi);
        goto S160;
S170:
        if(!(*status == -1)) goto S200;
        if(!qleft) goto S180;
        *status = 1;
        *bound = zero;
        goto S190;
S180:
        *status = 2;
        *bound = inf;
S200:
S190:
        ;
    }
    else if(4 == *which) {
        *pnonc = 5.0e0;
        T11 = -tent4;
        T12 = tent4;
        T13 = atol;
        T14 = tol;
        dstinv(&T11,&T12,&K3,&K3,&K4,&T13,&T14);
        *status = 0;
        dinvr(status,pnonc,&fx,&qleft,&qhi);
S210:
        if(!(*status == 1)) goto S220;
        cumtnc(t,df,pnonc,&cum,&ccum);
        fx = cum - *p;
        dinvr(status,pnonc,&fx,&qleft,&qhi);
        goto S210;
S220:
        if(!(*status == -1)) goto S250;
        if(!qleft) goto S230;
        *status = 1;
        *bound = 0.0e0;
        goto S240;
S230:
        *status = 2;
        *bound = tent4;
S240:
        ;
    }
S250:
    return;
#undef tent4
#undef tol
#undef atol
#undef zero
#undef one
#undef inf
}
