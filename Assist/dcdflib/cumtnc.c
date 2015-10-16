#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cdflib.h"

void cumtnc(double *t,double *df,double *pnonc,double *cum,
            double *ccum)
/**********************************************************************
 
     void cumtnc(double *t,double *df,double *pnonc,double *cum,
                 double *ccum)
 
                  CUMulative Non-Central T-distribution
 
 
                               Function
 
 
      Computes the integral from -infinity to T of the non-central
      t-density.
 
 
                               Arguments
 
 
      T --> Upper limit of integration of the non-central t-density.
 
      DF --> Degrees of freedom of the non-central t-distribution.
 
      PNONC --> Non-centrality parameter of the non-central t distibutio
 
      CUM <-- Cumulative t-distribution.
 
      CCUM <-- Compliment of Cumulative t-distribution.
 
 
                               Method
 
      Upper tail    of  the  cumulative  noncentral t   using
      formulae from page 532  of Johnson, Kotz,  Balakrishnan, Coninuous
      Univariate Distributions, Vol 2, 2nd Edition.  Wiley (1995)
 
      This implementation starts the calculation at i = lambda,
      which is near the largest Di.  It then sums forward and backward.
**********************************************************************/
{
#define one 1.0e0
#define zero 0.0e0
#define half 0.5e0
#define two 2.0e0
#define onep5 1.5e0
#define conv 1.0e-7
#define tiny 1.0e-10
static double alghdf,b,bb,bbcent,bcent,cent,d,dcent,dpnonc,dum1,dum2,e,ecent,
    halfdf,lambda,lnomx,lnx,omx,pnonc2,s,scent,ss,sscent,t2,term,tt,twoi,x,xi,
    xlnd,xlne;
static int ierr;
static unsigned long qrevs;
static double T1,T2,T3,T4,T5,T6,T7,T8,T9,T10;
/*
     ..
     .. Executable Statements ..
*/
/*
     Case pnonc essentially zero
*/
    if(fabs(*pnonc) <= tiny) {
        cumt(t,df,cum,ccum);
        return;
    }
    qrevs = *t < zero;
    if(qrevs) {
        tt = -*t;
        dpnonc = -*pnonc;
    }
    else  {
        tt = *t;
        dpnonc = *pnonc;
    }
    pnonc2 = dpnonc * dpnonc;
    t2 = tt * tt;
    if(fabs(tt) <= tiny) {
        T1 = -*pnonc;
        cumnor(&T1,cum,ccum);
        return;
    }
    lambda = half * pnonc2;
    x = *df / (*df + t2);
    omx = one - x;
    lnx = log(x);
    lnomx = log(omx);
    halfdf = half * *df;
    alghdf = gamln(&halfdf);
/*
     ******************** Case i = lambda
*/
    cent = fifidint(lambda);
    if(cent < one) cent = one;
/*
     Compute d=T(2i) in log space and offset by exp(-lambda)
*/
    T2 = cent + one;
    xlnd = cent * log(lambda) - gamln(&T2) - lambda;
    dcent = exp(xlnd);
/*
     Compute e=t(2i+1) in log space offset by exp(-lambda)
*/
    T3 = cent + onep5;
    xlne = (cent + half) * log(lambda) - gamln(&T3) - lambda;
    ecent = exp(xlne);
    if(dpnonc < zero) ecent = -ecent;
/*
     Compute bcent=B(2*cent)
*/
    T4 = cent + half;
    bratio(&halfdf,&T4,&x,&omx,&bcent,&dum1,&ierr);
/*
     compute bbcent=B(2*cent+1)
*/
    T5 = cent + one;
    bratio(&halfdf,&T5,&x,&omx,&bbcent,&dum2,&ierr);
/*
     Case bcent and bbcent are essentially zero
     Thus t is effectively infinite
*/
    if(bcent + bbcent < tiny) {
        if(qrevs) {
            *cum = zero;
            *ccum = one;
        }
        else  {
            *cum = one;
            *ccum = zero;
        }
        return;
    }
/*
     Case bcent and bbcent are essentially one
     Thus t is effectively zero
*/
    if(dum1 + dum2 < tiny) {
        T6 = -*pnonc;
        cumnor(&T6,cum,ccum);
        return;
    }
/*
     First term in ccum is D*B + E*BB
*/
    *ccum = dcent * bcent + ecent * bbcent;
/*
     compute s(cent) = B(2*(cent+1)) - B(2*cent))
*/
    T7 = halfdf + cent + half;
    T8 = cent + onep5;
    scent = gamln(&T7) - gamln(&T8) - alghdf + halfdf * lnx + (cent + half) * 
      lnomx;
    scent = exp(scent);
/*
     compute ss(cent) = B(2*cent+3) - B(2*cent+1)
*/
    T9 = halfdf + cent + one;
    T10 = cent + two;
    sscent = gamln(&T9) - gamln(&T10) - alghdf + halfdf * lnx + (cent + one) * 
      lnomx;
    sscent = exp(sscent);
/*
     ******************** Sum Forward
*/
    xi = cent + one;
    twoi = two * xi;
    d = dcent;
    e = ecent;
    b = bcent;
    bb = bbcent;
    s = scent;
    ss = sscent;
S10:
    b += s;
    bb += ss;
    d = lambda / xi * d;
    e = lambda / (xi + half) * e;
    term = d * b + e * bb;
    *ccum += term;
    s = s * omx * (*df + twoi - one) / (twoi + one);
    ss = ss * omx * (*df + twoi) / (twoi + two);
    xi += one;
    twoi = two * xi;
    if(fabs(term) > conv * *ccum) goto S10;
/*
     ******************** Sum Backward
*/
    xi = cent;
    twoi = two * xi;
    d = dcent;
    e = ecent;
    b = bcent;
    bb = bbcent;
    s = scent * (one + twoi) / ((*df + twoi - one) * omx);
    ss = sscent * (two + twoi) / ((*df + twoi) * omx);
S20:
    b -= s;
    bb -= ss;
    d *= (xi / lambda);
    e *= ((xi + half) / lambda);
    term = d * b + e * bb;
    *ccum += term;
    xi -= one;
    if(xi < half) goto S30;
    twoi = two * xi;
    s = s * (one + twoi) / ((*df + twoi - one) * omx);
    ss = ss * (two + twoi) / ((*df + twoi) * omx);
    if(fabs(term) > conv * *ccum) goto S20;
S30:
    if(qrevs) {
        *cum = half * *ccum;
        *ccum = one - *cum;
    }
    else  {
        *ccum = half * *ccum;
        *cum = one - *ccum;
    }
/*
     Due to roundoff error the answer may not lie between zero and one
     Force it to do so
*/
    *cum = fifdmax1(fifdmin1(*cum,one),zero);
    *ccum = fifdmax1(fifdmin1(*ccum,one),zero);
    return;
#undef one
#undef zero
#undef half
#undef two
#undef onep5
#undef conv
#undef tiny
}
