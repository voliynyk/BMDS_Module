C     GETCL : Computes the confidence limits on the BMD for the
c     multistage model using donlp2.
c     INPUT:
C         XWHICH      1 = lower confidence limit
C                     2 = upper confidence limit
c         XNDOSES     number of dose levels (INTEGER)
c         XDOSES(*)   dose levels used      (DOUBLE)
c         XAFFECT(*)  number of affected animals at Ith dose (INTEGER)
c         XNANIM(*)   number of animals at Ith dose (INTEGER)
c         XNPARMS     number of parameters in the model (INTEGER)
c         XBMR        Benchmark Response level requested, 
c                     0 < XBMR < 1 (DOUBLE)
c         XBMD        ML estimate of BMD (DOUBLE)
c         XTARG       Target likelihood value (MaxLL - X2/2) (DOUBLE)
c         XPARMS(*)   ML estimates of model parameters (0:polyord) (DOUBLE)
c         XFIXD(*)    1 if I'th (0:polyord) parameter is specified
c                     0 otherwise (INTEGER)
c         XVAL(*)     Value of fixed parameter (0:polyord) (DOUBLE)
c         XRISK       0 for absolute risk
C		      1 for standard deviations from control mean risk
c		      2 for relative risk
c		      3 for point risk
c		      4 for extra risk (Hill model ONLY)
c         XRESTR      1 if parameters are restricted to non-negative values
c                     0 otherwise 
c     OUTPUT:
c         XOPTITE     termination code (INTEGER)
c         BMDL        estimate of BMDL
c         XPARMS2     (0:XNPARMS) parameter estimates at BMDL
c		XADVERSE    1 if larger
c					-1 if smaller
c		XMODEL		0 if Polynomial
c					1 if Power
c					2 if Hill
C
C     Get  the lower confidence limit
      SUBROUTINE GETCL(XWHICH,XNDOSES,XDOSES,XMEAN,XNANIM,XVAR,XNPARMS
     $     ,XBMR,XBMD,XTARG,XPARMS,XFIXD,XVAL,XRISK,XRESTR,BMDL,XPARMS2
     $     ,XOPTITE,XNRESM, XBIND, XADVERSE, XMODEL, XFLAG, XMLEPARM)
      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER XWHICH,XNDOSES,XNANIM(*),XNPARMS
     $     ,XFIXD(0:XNPARMS-1),XRISK,XRESTR,XOPTITE,XNRESM,
     $     XBIND(0:XNPARMS-1),XADVERSE, XMODEL, XFLAG
      DOUBLE PRECISION XDOSES(*),XBMR,XBMD,XTARG,XPARMS(0:XNPARMS-1)
     $     ,XVAL(0:XNPARMS-1), BMDL, XPARMS2(0:XNPARMS-1),XMEAN(*),
     $	 XVAR(*), XMLEPARM(0:XNPARMS-1)
      INTEGER I
      REAL EPS
      EPS = .00000001
C -------------------------------------
C     Open a file to print debugging information
      if (iDebug .gt. 0) then
         open(31, file='poly_debug.out', access='APPEND')
         write (31, *) 'Entered GETCL.'
      endif
C Define the problem type
      probtype = 2
	modtype = XMODEL
	flag = XFLAG

	clcnt = clcnt + 1

C Load values into common blocks
      ndoses = XNDOSES
      DO I = 1, ndoses
	   var(I) = XVAR(I)
         doses(I) = XDOSES(I)
         mean(I) = XMEAN(I)
         nanimals(I) = XNANIM(I)
      ENDDO
      nparm = XNPARMS + 1
      N = nparm
      bmr=XBMR
      bmd=XBMD
      target = XTARG
      X(1) = bmd 

      DO I = 0, nparm - 2
        X(I+2) = XPARMS(I)
         parmfixd(I) = XFIXD(I)
         parmval(I) = XVAL(I)
      ENDDO

      risktype = XRISK
      restrict = XRESTR
      adverse = XADVERSE
     

c
C Finally, do the work
      call donlp2


      DO I = 1, nparm
         X(I) = X(I) *XSC(I)
      ENDDO

      BMDL = X(1)
      DO I=0, XNPARMS-1
         XPARMS2(I) = X(I+2)
      ENDDO
      XOPTITE = OPTITE
c Adjust nparm back to nparm-1 
	  nparm = nparm -1
C     Is the Ith parameter fixed or stuck at a lower bound?
c      DO I=3, NG+NH
c         IF (GUNIT(1,I) .EQ. 1 .AND. GUNIT(2,I) .GT. 1) XBIND(GUNIT(2,I)
c     $        -2) = bind(I)  
c      ENDDO

      if (iDebug .gt. 0) then
         write (31, *) 'Exiting GETCL'
         close(31)
         iDebug = 0
      endif
      RETURN
      END
