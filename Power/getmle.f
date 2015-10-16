C     GETMLE : Estimates model parameters
c     INPUT:
c         XNDOSES     number of dose levels (INTEGER)
c         XDOSES(*)   dose levels used      (DOUBLE)
c         XMEAN(*)    the sample mean at Ith dose (DOUBLE)
c         XNANIM(*)   number of animals at Ith dose (INTEGER)
c	  XVAR(*)     sample variance at Ith dose(DOUBLE)
c         XNPARMS     the number of parameters in the model (INTEGER)
c         XPARMS(*)   Start values for model parameters (0:polyord) (DOUBLE)
c         XFIXD(*)    1 if I'th (0:polyord) parameter is specified
c                     0 otherwise (INTEGER)
c         XVAL(*)     Value of fixed parameter (0:polyord) (DOUBLE)
c         XRESTR      1 if parameters are restricted to non-negative values
c                     0 if there is no restriction
c                     -1 if parameters are restricted to non-positive values
c	  XADVERSE    1 if adverse direction is larger
c		      -1 if adverse direction is smaller
c	  XMODEL      0 if Polynomial
c		      1 if Power
c		      2 if Hill
c					
c     OUTPUT:
c         XOPTITE     termination code (INTEGER)
c         XPARMS2     (0:XPOLYORD) parameter estimates
C         XLL         value of the log likelihood function at the maximum
C
      SUBROUTINE GETMLE(XNDOSES,XDOSES,XMEAN,XNANIM,XVAR,
     $     XNPARMS,XPARMS,XFIXD,XVAL,XRESTR,XADVERSE,
     $     XPARMS2,XLL,XOPTITE,XNRESM,XBIND,XMODEL,XFLAG)
      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER XNDOSES,XNANIM(*),XNPARMS,XADVERSE,XMODEL,XFLAG,
     $     XFIXD(0:XNPARMS-1),XRESTR,XOPTITE,XNRESM, XBIND(0:XNPARMS-1)
      DOUBLE PRECISION XDOSES(*),XPARMS(0:XNPARMS-1),XLL,
     $     XVAL(0:XNPARMS-1),XPARMS2(0:XNPARMS-1),XMEAN(*),XVAR(*)
      INTEGER I
      
C     -------------------------------------
      mlecnt = mlecnt + 1

C     Define the problem type
      probtype = 1
      modtype = XMODEL

C     Load values into common blocks
      ndoses = XNDOSES
      DO I = 1, ndoses
         var(I) = XVAR(I)
         doses(I) = XDOSES(I)
         mean(I) = XMEAN(I)
         nanimals(I) = XNANIM(I)
      ENDDO
      adverse = XADVERSE
      flag = XFLAG
      nparm = XNPARMS
      N = nparm
      DO I = 0, nparm-1
         X(I+1) = XPARMS(I)
         parmfixd(I) = XFIXD(I)
         parmval(I) = XVAL(I)
      ENDDO
      IF(parmfixd(1).EQ.1.AND.parmval(1).EQ.0) THEN
         constvar = 1
      ELSE
         constvar = 0
      ENDIF
      restrict = XRESTR
C     
c     Scale by the starting values
      DO I = 1, N
         IF (X(I) .NE. 0.) THEN
            XSC(I) = X(I)
         ELSE
            XSC(I) = 1
         ENDIF
      ENDDO
      
C     Finally, do the work
      call donlp2
c     
c     Unscale the parameters by starting values
      DO I=1, nparm
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
