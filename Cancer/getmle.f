C     GETMLE : Estimates model parameters
c     INPUT:
c         XNDOSES     number of dose levels (INTEGER)
c         XDOSES(*)   dose levels used      (DOUBLE)
c         XAFFECT(*)  number of affected animals at Ith dose (INTEGER)
c         XNANIM(*)   number of animals at Ith dose (INTEGER)
c         XPOLYORD    polynomial order requested (INTEGER)
c         XPARMS(*)   Start values for model parameters (0:polyord) (DOUBLE)
c         XFIXD(*)    1 if I'th (0:polyord) parameter is specified
c                     0 otherwise (INTEGER)
c         XVAL(*)     Value of fixed parameter (0:polyord) (DOUBLE)
c         XRESTR      1 if parameters are restricted to non-negative values
c                     0 otherwise 
c     OUTPUT:
c         XOPTITE     termination code (INTEGER)
c         XPARMS2     (0:XPOLYORD) parameter estimates
C         XLL         value of the log likelihood function at the maximum
C
      SUBROUTINE GETMLE(XPARMS,XFIXD,XVAL,XPARMS2,XLL
     $     ,XOPTITE,XNRESM,XBIND)
      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER XFIXD(0:polyord),XOPTITE,XNRESM, XBIND(0:polyord) 
      DOUBLE PRECISION XPARMS(0:polyord),XVAL(0:polyord),
     $      XPARMS2(0:polyord),XLL, temp1
      INTEGER I, temp2, temp3
c      REAL maxdose
C -------------------------------------
C Define the problem type
      probtype = 1
C Load values into common blocks
c      ndoses = XNDOSES
c      DO I = 1, ndoses
c         doses(I) = XDOSES(I)
c         affected(I) = XAFFECT(I)
c         nanimals(I) = XNANIM(I)
c      ENDDO
c      polyord = XPOLYORD
      N = polyord + 1
      DO I = 0, polyord
         X(I+1) = XPARMS(I)
         parmfixd(I) = XFIXD(I)
         parmval(I) = XVAL(I)
      ENDDO
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
c      
c      DO I = 1, ndoses
c         doses(I) = doses(I)/maxdose
c      ENDDO
c      bmd = bmd/maxdose
c      DO I = 2, N
c         X(I) = X(I)*(maxdose ** (I-1))
c         IF (parmfixd(I) .GT. 0) parmval(I) = parmval(I)*(maxdose ** (I
c     $        - 1)) 
c      ENDDO
C Scale by the starting values
c      DO I = 1, N
c         IF (X(I) .GT. 0.) THEN
c            XSC(I) = X(I)
c         ELSE
c            XSC(I) = 1
c         ENDIF
c      ENDDO
	
C Finally, do the work
      call donlp2

C Unscale by starting values	
c      DO I=1,N
c         X(I) = X(I)*XSC(I)
c      ENDDO
c Unscale by maxdose
c      DO I=2,N
c         X(I) = X(I)/(maxdose**(I-1))
c      ENDDO
C     Transfer the parameter estimates, and initialize XBIND(I)
      DO I=0, polyord
         XPARMS2(I) = X(I+1)
         XBIND(I) = 0
c         PRINT *,XPARMS2(I)
      ENDDO

C     Set the appropriate values in XBIND (donlp2 indexes parameters from 1
C     but XBIND indexes parameters from 0
      XNRESM = NG+NH
      DO I = 1, XNRESM
         IF (BIND(I) .NE. 0) XBIND(GUNIT(2 , I) - 1) = 1
      ENDDO
      XOPTITE = OPTITE
      XLL = -FX
      RETURN
      END
