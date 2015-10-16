C     GETMLE : Estimates model parameters
c     INPUT:
c         XNPARMS		the number of parameters in the model (INTEGER)
c         XPARMS(*)   Start values for model parameters (0:polyord) (DOUBLE)
c					
c     OUTPUT:
c         XOPTITE     termination code (INTEGER)
c         XPARMS2     (0:XPOLYORD) parameter estimates
C         XLL         value of the log likelihood function at the maximum
C
c
      SUBROUTINE GETMLE(XNPARMS,XPARMS,XPARMS2,XLL,XOPTITE,
     $     XNRESM,XBIND,XFLAG)
      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER XNPARMS,XFLAG,XOPTITE,XNRESM, XBIND(0:XNPARMS)
      DOUBLE PRECISION XPARMS(0:XNPARMS-1),XLL,XPARMS2(0:XNPARMS-1)
      INTEGER I, iRunCount
      DATA iDebug/0/
      DATA iRunCount/0/
      SAVE iRunCount
C     -------------------------------------
C     Open a file to print debugging information
      if (iDebug .gt. 0) then
         iRunCount = iRunCount + 1
         open(31, file='poly_debug.out', access='APPEND')
         write (31, *) 'Entered GETMLE. RunCount= ', iRunCount
      endif
C     Define the problem type
      mlecnt = mlecnt + 1
      probtype = 1
      flag = XFLAG
      nparm = XNPARMS
      N = nparm

      DO I = 0, nparm-1
         X(I+1) = XPARMS(I)
      ENDDO
C
C     Finally, do the work
      call donlp2

C     This is to unscale the parameters by their starting values
C     if they get scaled in donlp2 SETUP0
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

      if (iDebug .gt. 0) then
         write (31, *) 'Exiting GETMLE'
         close(31)
         iDebug = 0
      endif
      RETURN
      END
