C     GETCL : Computes the confidence limits on the BMD for the
c     multistage model using donlp2.
c     INPUT:
C         XWHICH      1 = lower confidence limit
C                     2 = upper confidence limit
c         XNDOSES     number of dose levels (INTEGER)
c         XDOSES(*)   dose levels used      (DOUBLE)
c         XAFFECT(*)  number of affected animals at Ith dose (INTEGER)
c         XNANIM(*)   number of animals at Ith dose (INTEGER)
c         XPOLYORD    polynomial order requested (INTEGER)
c         XBMR        Benchmark Response level requested, 
c                     0 < XBMR < 1 (DOUBLE)
c         XBMD        ML estimate of BMD (DOUBLE)
c         XTARG       Target likelihood value (MaxLL - X2/2) (DOUBLE)
c         XPARMS(*)   ML estimates of model parameters (0:polyord) (DOUBLE)
c         XFIXD(*)    1 if I'th (0:polyord) parameter is specified
c                     0 otherwise (INTEGER)
c         XVAL(*)     Value of fixed parameter (0:polyord) (DOUBLE)
c         XRISK       1 for extra risk, 2 for additional risk (INTEGER)
c         XRESTR      1 if parameters are restricted to non-negative values
c                     0 otherwise 
c     OUTPUT:
c         XOPTITE     termination code (INTEGER)
c         BMDL        estimate of BMDL
c         XPARMS2     (0:XPOLYORD) parameter estimates at BMDL
C
C     Get  the lower confidence limit
c      SUBROUTINE GETCL(XWHICH,XNDOSES,XDOSES,XAFFECT,XNANIM,XPOLYORD
c     $     ,XBMR,XBMD,XTARG,XPARMS,XFIXD,XVAL,XRISK,XRESTR,BMDL,XPARMS2
c     $     ,XOPTITE,XNRESM, XBIND)
c      INCLUDE 'O8COMM.INC'
c      INCLUDE 'O8FINT.INC'
c      INCLUDE 'PROBLEM.INC'
c      INTEGER XWHICH,XNDOSES,XAFFECT(*),XNANIM(*),XPOLYORD
c     $     ,XFIXD(0:XPOLYORD),XRISK,XRESTR,XOPTITE,XNRESM,
c     $     XBIND(0:XPOLYORD) 
c      DOUBLE PRECISION XDOSES(*),XBMR,XBMD,XTARG,XPARMS(0:XPOLYORD)
c     $     ,XVAL(0:XPOLYORD), BMDL, XPARMS2(0:XPOLYORD), temp1
c      INTEGER I, temp2, temp3

      SUBROUTINE GETCL(XWHICH,XBMR,XBMD,XTARG,XPARMS,XFIXD,XVAL,XRISK,
     $     BMDL,XPARMS2,XOPTITE,XNRESM, XBIND, XBMDbnd)
      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER XWHICH,XFIXD(0:polyord),XRISK,XOPTITE,XNRESM,
     $     XBIND(0:polyord), XBMDbnd
      DOUBLE PRECISION XBMR,XBMD,XTARG,XPARMS(0:polyord),
     $     XVAL(0:polyord),BMDL,XPARMS2(0:polyord)
      INTEGER I
      DOUBLE PRECISION temp1, temp2, temp3
      DATA iDebug/0/
c      REAL maxdose
C -------------------------------------
C Define the problem type
      probtype = XWHICH + 1
C Load values into common blocks
      N = polyord + 2
      bmr=XBMR
C Work with log(bmd)
      bmd=LOG(XBMD)
      target = XTARG
      X(1) = bmd
      DO I = 0, polyord
         X(I+2) = XPARMS(I)
         parmfixd(I) = XFIXD(I)
         parmval(I) = XVAL(I)
      ENDDO
      risktype = XRISK

	
C Finally, do the work
      call donlp2

      BMDL = EXP(X(1))
      DO I=0, polyord
         XPARMS2(I) = X(I+2)
      ENDDO
      XOPTITE = OPTITE
C     Is the Ith parameter fixed or stuck at a lower bound?

      DO I=3, NG+NH
         IF (GUNIT(1,I) .EQ. 1 .AND. GUNIT(2,I) .GT. 1) XBIND(GUNIT(2,I)
     $        -2) = bind(I)  
      ENDDO
C     Foolishly, the upper and lower bounds for bmdu, bmdl correspond
C     to different constraints.
C     Is the BMD at its lower (for BMDL: XWHICH = 1) or upper
C     (for BMDU: XWHICH = 2) bound?
      IF (XWHICH .EQ. 1) THEN
         XBMDbnd = bind(NH + 2)
      ELSE
         XBMDbnd = bind(NH + 3)
      ENDIF
C     Reset the debug flag
      if (iDebug .gt. 0) then
         close(31)
         iDebug = 0
      endif
      RETURN
      END
