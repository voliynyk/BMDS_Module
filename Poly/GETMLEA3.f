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
c					-1 if parameters are restricted to non-positive values
c		XADVERSE	1 if adverse direction is larger
c					-1 if adverse direction is smaller
c					1 if the Power model is desired
c					2 if the Hill model is desired
c     OUTPUT:
c         XOPTITE     termination code (INTEGER)
c         XPARMS2     (0:XPOLYORD) parameter estimates
C         XLL         value of the log likelihood function at the maximum
C
      SUBROUTINE GETMLEA3(XNDOSES,XDOSES,XMEAN,XNANIM,XVAR
     $     ,XNPARMS,XPARMS,XFIXD,XVAL,XRESTR
     $     ,XPARMS2,XLL,XOPTITE,XNRESM,XBIND)
      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER XNDOSES,XNANIM(*),XNPARMS
     $     ,XFIXD(0:XNPARMS),XRESTR,XOPTITE,XNRESM, XBIND(0:XNPARMS)
      DOUBLE PRECISION XDOSES(*),XPARMS(0:XNPARMS),XLL
     $     ,XVAL(0:XNPARMS),XPARMS2(0:XNPARMS),XMEAN(*),XVAR(*)
      INTEGER I
      
C -------------------------------------
      a3cnt = a3cnt + 1
      
C     Define the problem type
      probtype = 3
C     Load values into common blocks
      ndoses = XNDOSES
      DO I = 1, ndoses
         var(I) = XVAR(I)
         doses(I) = XDOSES(I)
         mean(I) = XMEAN(I)
         nanimals(I) = XNANIM(I)
      ENDDO
      nparm = XNPARMS
      N = nparm
      DO I = 0, nparm-1
         X(I+1) = XPARMS(I)
         parmfixd(I) = XFIXD(I)
         parmval(I) = XVAL(I)
      ENDDO
      constvar = 0
      restrict = XRESTR
c
C     Finally, do the work
      call donlp2
C
      DO I=0, XNPARMS-1
         XPARMS2(I) = X(I+1)
      ENDDO
      XOPTITE = OPTITE
      XNRESM = NG+NH
      XLL = -FX
C     Is the Ith parameter fixed or stuck at a lower bound?
      RETURN
      END
