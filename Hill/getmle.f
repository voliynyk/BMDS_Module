C     GETMLE : Estimates model parameters
c     INPUT:
c         XNDOSES     number of dose levels (INTEGER)
c         XDOSES(*)   dose levels used      (DOUBLE)
c         XMEAN(*)	the sample mean at Ith dose (DOUBLE)
c         XNANIM(*)   number of animals at Ith dose (INTEGER)
c		XVAR(*)		sample variance at Ith dose(DOUBLE)
c         XNPARMS		the number of parameters in the model (INTEGER)
c         XPARMS(*)   Start values for model parameters (0:polyord) (DOUBLE)
c         XFIXD(*)    1 if I'th (0:polyord) parameter is specified
c                     0 otherwise (INTEGER)
c         XVAL(*)     Value of fixed parameter (0:polyord) (DOUBLE)
c         XRESTR      1 if parameters are restricted to non-negative values
c                     0 if there is no restriction
c					-1 if restricted to negative values
c					-1 if parameters are restricted to non-positive values
c		XADVERSE	1 if adverse direction is larger
c					-1 if adverse direction is smaller
c		XMODEL		0 if Polynomial
c					1 if Power
c					2 if Hill
c					
c     OUTPUT:
c         XOPTITE     termination code (INTEGER)
c         XPARMS2     (0:XPOLYORD) parameter estimates
C         XLL         value of the log likelihood function at the maximum
C
      SUBROUTINE GETMLE(XNPARMS,XPARMS,XPARMS2,XLL,XOPTITE,XNRESM,XBIND)
      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER XWHICH,XNPARMS,XOPTITE,XNRESM, XBIND(0:XNPARMS-1)
      DOUBLE PRECISION XPARMS(0:XNPARMS-1),XLL,XPARMS2(0:XNPARMS-1)
      INTEGER I,J
      DOUBLE PRECISION maxdose
C -------------------------------------
C Define the problem type
      probtype = 1
c	modtype = XMODEL
C Load values into common blocks
c      ndoses = XNDOSES
c      DO I = 1, ndoses
c	   var(I) = XVAR(I)
c         doses(I) = XDOSES(I)
c         mean(I) = XMEAN(I)
c         nanimals(I) = XNANIM(I)
c      ENDDO
c	adverse = XADVERSE
c
c     This counts the number of calls to getmle
      mlecnt = mlecnt + 1
      nparm = XNPARMS
      N = nparm
      DO I = 0, nparm-1
         X(I+1) = XPARMS(I)
c         parmfixd(I) = XFIXD(I)
c         parmval(I) = XVAL(I)
      ENDDO
c	IF(parmfixd(1).EQ.1.AND.parmval(1).EQ.0) THEN
c		constvar = 1
c	ELSE
c		constvar = 0
c	ENDIF
c      restrict = XRESTR
C
c	Scale by the starting values
C      DO I = 1, N
C         IF (X(I) .NE. 0.) THEN
C            XSC(I) = X(I)
C         ELSE
C            XSC(I) = 1
C         ENDIF
C      ENDDO
C      WRITE(*, 1000)
C 1000 FORMAT(1x)
C      WRITE(*, 1001) N
C 1001 FORMAT(1x,'N =',I4)
C      DO I = 1, N
C         WRITE(*, 1002) I, X(I), PARMFIXD(I-1),PARMVAL(I-1)
C 1002    FORMAT(1X,'I = ',I2,' X(I)= ',G10.4,' FIXED(I) = ',I2, 
C     1   ' VAL(I) =', G10.4)
C      ENDDO
C Finally, do the work
      call donlp2


c     unscale parameters
      DO I = 1, nparm
         X(I) = X(I)*XSC(I)
      ENDDO

C     Transfer the parameter estimates, and initialize XBIND(I)
      DO I=0, XNPARMS-1
         XPARMS2(I) = X(I+1)
         XBIND(I) = 0
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
