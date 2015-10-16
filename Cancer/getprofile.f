C     GETPROFILE : Computes the maximum likelihood for a given BMD for the
c     multistage model using donlp2.
c     INPUT:
c         XNDOSES     number of dose levels (INTEGER)
c         XDOSES(*)   dose levels used      (DOUBLE)
c         XAFFECT(*)  number of affected animals at Ith dose (INTEGER)
c         XNANIM(*)   number of animals at Ith dose (INTEGER)
c         XPOLYORD    polynomial order requested (INTEGER)
c         XPARMS(*)   ML estimates should not be changed, and may
c                     not be used (0:POLYORD) (DOUBLE)
c         XBMR        Benchmark Response level requested, 
c                     0 < XBMR < 1 (DOUBLE)
c         XBMD        ML estimate of BMD (DOUBLE)
c         XTARG       Target likelihood value (MaxLL - X2/2) (DOUBLE)
c         XFIXD(*)    1 if I'th (0:polyord) parameter is specified
c                     0 otherwise (INTEGER)
c         XVAL(*)     Value of fixed parameter (0:polyord) (DOUBLE)
c         XRISK       1 for extra risk, 2 for additional risk (INTEGER)
c         XRESTR      1 if parameters are restricted to non-negative values
c                     0 otherwise 
c     OUTPUT:
c         XOPTITE     termination code (INTEGER)
c         PLL         profile log likelihood 
c         XPARMS3     (0:XPOLYORD) parameter estimates at max likelihood given BMD
C
C     Get  the maximum likelihood profile
      SUBROUTINE GETPROFILE(XPARMS,XBMD,PLL,XOPTITE,TMPCNT,XBMR)
      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
c      INTEGER XNDOSES,XAFFECT(*),XNANIM(*),XPOLYORD
c     $     ,XFIXD(0:XPOLYORD),XRISK,XRESTR,XOPTITE,XNRESM,
c     $     XBIND(0:XPOLYORD) 
c      DOUBLE PRECISION XPARMS(0:XPOLYORD),XDOSES(*),XBMR,XBMD,XTARG,XVAL(0:XPOLYORD),XPARMS3(0:XPOLYORD),temp1
      INTEGER XOPTITE, TMPCNT, XNRESM, XBIND(0:polyord)
      DOUBLE PRECISION XPARMS(0:polyord+1),XBMD,PLL, temp1,XBMR
      INTEGER I, temp2, temp3
c      REAL maxdose
C -------------------------------------
C Define the problem type
      probtype = 4

C Load values into common blocks
      counter = TMPCNT
c      ndoses = XNDOSES
c      PRINT *," "
c      DO I = 1, ndoses
c         PRINT *,"DOSE",I,"= ",doses(I)
c         PRINT *,"AFFECTED",I,"= ",affected(I)
c         PRINT *,"NANIM",I,"= ",nanimals(I)
c      ENDDO
c      polyord = XPOLYORD
      N = polyord + 1
      bmr=XBMR
      bmd=XBMD
c      PRINT *,"N =",N
c      PRINT *,"BMR =",bmr
c      PRINT *,"BMD =",bmd
      DO I = 1, N
         X(I) = XPARMS(I)
c         parmfixd(I) = XFIXD(I)
c         parmval(I) = XVAL(I)
      ENDDO
c      risktype = XRISK
c      restrict = XRESTR

C     Make sure control dose level is first in unrestricted case
c     This is needed for donlp2usrfc.f to work correctly 
c      IF (restrict .NE. 1) THEN
c	   DO I = 2, ndoses 
c	      IF (doses(I) .LT. doses(1)) THEN
c	         temp1 = doses(1)
c	         temp2 = affected(1)
c	         temp3 = nanimals(1)
c	         doses(1) = doses(I)
c	         affected(1) = affected(I)
c	         nanimals(1) = nanimals(I)
c                 doses(I) = temp1
c	         affected(I) = temp2
c	         nanimals(I) = temp3
c	       ENDIF
c	    ENDDO
c	ENDIF       
C ------------------------------------------
C
C Rescale the doses by the largest dose
c      maxdose = doses(1)
c      DO I = 2, ndoses
c         IF (doses(I) .GT. maxdose) maxdose = doses(I)
c      ENDDO
c      DO I = 1, ndoses
c         doses(I) = doses(I)/maxdose
c      ENDDO
      bmd = bmd/maxdose
      DO I = 2, N
         X(I) = X(I)*(maxdose ** (I-1))
c         IF (parmfixd(I-1) .GT. 0) parmval(I-1) = parmval(I-1)*(maxdose
c     $    ** (I - 1)) 
      ENDDO
c      DO I=1,polyord+1
c         PRINT *,X(I),parmfixd(I-1)
c      ENDDO
C Scale by the starting values
c      DO I = 1, N
c         IF (X(I) .GT. 0.) THEN
c            XSC(I) = X(I)
c         ELSE
c            XSC(I) = 1.0
c         ENDIF
c      ENDDO
C Finally, do the work
      call donlp2

C Unscale by starting values
c      DO I=1,N
c         X(I) = X(I)*XSC(I)
c      ENDDO

c Unscale by maxdose
      DO I=2,polyord+1
         X(I) = X(I)/(maxdose**(I-1))
      ENDDO
      DO I=1, polyord+1
         XPARMS(I) = X(I)
      ENDDO
      IF ((optite .EQ. 0) .OR. (optite .EQ. 1) .OR. 
     $     (optite .EQ. 2)) THEN
         XOPTITE = 0;
      ELSE
         XOPTITE = 1;
       ENDIF
c       PRINT *,"OPTITE =",OPTITE
C      XNRESM = NG+NH
c      DO I=0, XNRESM-1
c         PRINT *,I,bind(1),XNRESM,polyord
c         XBIND(I) = bind(I+1)
c      ENDDO
      PLL = -FX
c      PRINT *,"The LL =",PLL
C     Is the Ith parameter fixed or stuck at a lower bound?
c      DO I=1, NG+NH
c         IF (GUNIT(1,I) .EQ. 1 .AND. GUNIT(2,I) .GT. 1) XBIND(GUNIT(2,I)
c     $        -1) = bind(I)  
c      ENDDO
      RETURN
      END
