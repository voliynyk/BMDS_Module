C     GETPROFILE: Compute the maximum likelihood for a give BMD for 
C				   the polynominal model . 
c     INPUT:
c         XNDOSES     number of dose levels (INTEGER)
c         XDOSES(*)   dose levels used      (DOUBLE)
c         XMEAN(*)	the sample mean at Ith dose (DOUBLE)
c         XNANIM(*)   number of animals at Ith dose (INTEGER)
c		  XVAR(*)		sample variance at Ith dose(DOUBLE)
c         XNPARMS		the number of parameters in the model (INTEGER)
c         XPARMS(*)   Start values for model parameters (0:polyord) (DOUBLE)
c         XFIXD(*)    1 if I'th (0:polyord) parameter is specified
c                     0 otherwise (INTEGER)
c         XVAL(*)     Value of fixed parameter (0:polyord) (DOUBLE)
c         XBMD		   A GIVEN BMD
c		  XRESTR      1 if parameters are restricted to non-negative values
c                     0 if there is no restriction
c					-1 if restricted to negative values
c					-1 if parameters are restricted to non-positive values
c		XADVERSE	1 if adverse direction is larger
c					-1 if adverse direction is smaller
c         XRISK       0 for absolute risk 1
C					1 for standard deviations from control mean risk
c					2 for relative risk
c					3 for point risk
c					4 for extra risk (Hill model ONLY)
c		XMODEL		0 if Polynomial
c					1 if Power
c					2 if Hill
c		XFLAG		flag used for getprofile is the same as that for getmle
c		XITER		iteration number used to control flag			
c     OUTPUT:
c         XOPTITE     termination code (INTEGER)
c         XPARMS2     (0:XPOLYORD) parameter estimates
C         XLL         value of the reduced log likelihood function at the maximum
C
      SUBROUTINE GETPROFILE (XPARMS, XBMD, XLL,XOPTITE, XITER, XBMR)

      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER XOPTITE, XITER, I
      double precision XPARMS(0:nparm), XBMD, XLL, XBMR
c      SUBROUTINE ( XNDOSES,XDOSES,XMEAN,XNANIM,XVAR
C    $     ,XNPARMS,XPARMS,XFIXD,XVAL,XRESTR,XRISK, XADVERSE
C     $     ,XPARMS2,XLL,XOPTITE,XNRESM,XBIND,XMODEL,XFLAG)
c      INTEGER XNDOSES,XNANIM(*),XNPARMS,XADVERSE,XRISK, XMODEL,XFLAG
c     $     ,XFIXD(0:XNPARMS-1),XRESTR,XOPTITE,XNRESM, XBIND(0:XNPARMS-1)
c      DOUBLE PRECISION XBMD, XDOSES(*),XPARMS(0:XNPARMS-1),PLL
c    $     ,XVAL(0:XNPARMS-1),XPARMS2(0:XNPARMS-1),XMEAN(*),XVAR(*)
      
      
C -------------------------------------
     
C Define the problem type
      probtype = 4
      counter = XITER
      modtype = 0
c      if (xiter .le. 14) then 
         flag =0
c      elseif (xiter .ge. 15 .and. xiter .le. 24) then
c         flag =1
c      elseif (xiter .ge. 25 .and. xiter .le. 34) then
c         flag =2
c      else 
c         flag = 3
c      endif 
      
c     flag = XFLAG
C     Load values into common blocks
c     ndoses = XNDOSES
c     DO I = 1, ndoses
c     var(I) = XVAR(I)
c     doses(I) = XDOSES(I)
c     mean(I) = XMEAN(I)
c     nanimals(I) = XNANIM(I)
c     ENDDO
c     adverse = XADVERSE
c     
c     nparm = XNPARMS

      bmr = XBMR
      N = nparm
      
      IF(parmfixd(1).EQ.1.AND.parmval(1).EQ.0) THEN
         constvar = 1
      ELSE
         constvar = 0
      ENDIF
      
c     restrict = XRESTR
c     risktype = XRISK
      bmd = XBMD/maxdose 
c      PRINT *, maxdose
c     Finally, do the work
      
      
c     print*, 'print modty,flag bmr restrict risktype bmd probtype modtype adverse'
c     print*, modtype,flag, bmr, restrict, risktype, bmd, probtype, modtype, adverse 
      
      
      
      DO I = 1, 3
         X(I) = XPARMS(I)
c     parmfixd(I) = XFIXD(I)
c     parmval(I) = XVAL(I)
      ENDDO
      DO I = 4, nparm
         X(I) = XPARMS(I)*maxdose**(I-3)
      ENDDO
      
      call donlp2
      
      DO I = 1, nparm
         X(I) = X(I)*XSC(I)
      ENDDO
      
c     Passing back newly estimated parameters
      DO I=1, 3
         XPARMS(I) = X(I)
      ENDDO
      DO I = 4, nparm
         XPARMS(I) = X(I)/maxdose**(I-3)
      ENDDO

      IF ((optite .EQ. 0) .OR. (optite .EQ. 1) .OR. 
     $     (optite .EQ. 2)) THEN
         XOPTITE = 0;
      ELSE
         XOPTITE = 1;
      ENDIF

      XLL = -FX
      
c     XNRESM = NG+NH
c     
C     Is the Ith parameter fixed or stuck at a lower bound?
c     
c     DO I=0, XNRESM-1
c     XBIND(I) = bind(I+1)
c     ENDDO
c      DO I=1, NG+NH
c     IF (GUNIT(1,I) .EQ. 1 .AND. GUNIT(2,I) .GT. 1) XBIND(GUNIT(2,I)
c     $        -1) = bind(I)
c     ENDDO
      
      RETURN
      END
