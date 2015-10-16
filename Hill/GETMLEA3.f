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
      SUBROUTINE GETMLEA3(XNPARMS,XPARMS,XPARMS2,XLL,XOPTITE,XNRESM,
     $     XBIND)
      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER XNPARMS,XOPTITE,XNRESM, XBIND(0:XNPARMS)
      DOUBLE PRECISION XPARMS(0:XNPARMS),XLL,XPARMS2(0:XNPARMS)
      INTEGER I, temp_constvar
      
C -------------------------------------
C Define the problem type
      probtype = 3
c     This counts the calls the getmlea3
      a3cnt = a3cnt + 1
      nparm = XNPARMS
      N = nparm
      DO I = 0, nparm-1
         X(I+1) = XPARMS(I)
      ENDDO
      temp_constvar = constvar
C Finally, do the work
      call donlp2
	DO I=0, XNPARMS-1
         XPARMS2(I) = X(I+1)
      ENDDO
      XOPTITE = OPTITE
      XNRESM = NG+NH
      DO I=0, XNRESM-1
         XBIND(I) = bind(I+1)
      ENDDO
      XLL = -FX
C     Is the Ith parameter fixed or stuck at a lower bound?
      DO I=1, NG+NH
         IF (GUNIT(1,I) .EQ. 1 .AND. GUNIT(2,I) .GT. 1) XBIND(GUNIT(2,I)
     $        -1) = bind(I)  
      ENDDO
      constvar = temp_constvar
      RETURN
      END
