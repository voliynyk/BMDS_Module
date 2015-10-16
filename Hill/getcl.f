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
c         XRISK       0 for absolute risk 1
C					1 for standard deviations from control mean risk
c					2 for relative risk
c					3 for point risk
c					4 for extra risk (Hill model ONLY)
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
      SUBROUTINE GETCL(XWHICH,XNPARMS,XBMR,XBMD,XTARG,XPARMS,XRISK,
     $     BMDL,XPARMS2,XOPTITE,XNRESM,XBIND,XFLAG)
      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER XWHICH,XNPARMS,XRISK,XOPTITE,XNRESM,
     $     XBIND(0:XNPARMS), XFLAG
      DOUBLE PRECISION XBMR,XBMD,XTARG,XPARMS(0:XNPARMS),
     $     BMDL, XPARMS2(0:XNPARMS)
      INTEGER I,J,K
      REAL maxdose,SIGN,Vi,ABMN,Ai,Ni,Bi,Gi,Hi,TEMP,Devi,GRAD(23),EPS
      EPS = 2.220446D-16
C -------------------------------------
C Define the problem type
      probtype = 2
c	modtype = XMODEL
C Load values into common blocks
c      ndoses = XNDOSES
c      DO I = 1, ndoses
c	   var(I) = XVAR(I)
c         doses(I) = XDOSES(I)
c         mean(I) = XMEAN(I)
c         nanimals(I) = XNANIM(I)
c      ENDDO
c     This counts the calls to getcl
      clcnt = clcnt + 1
      nparm = XNPARMS + 1
      N = nparm
      bmr=XBMR
      bmd=XBMD
	flag = XFLAG
      target = XTARG
      X(1) = bmd
      DO I = 0, nparm - 2
         X(I+2) = XPARMS(I)
c         parmfixd(I) = XFIXD(I)
c         parmval(I) = XVAL(I)
      ENDDO
      risktype = XRISK
c      restrict = XRESTR
c	adverse = XADVERSE
c
	CALL hillpart(X)
c
c    Get gradients for scaling factors
c
	DO J = 1,nparm
		GRAD(J) = 0.D0
	ENDDO
		DO K = 1, ndoses
			IF(means(K).LT.0) THEN
			  SIGN = -1
			ELSE
			  SIGN = 1
			ENDIF
			ABMN = ABS(means(K))
			IF(constvar.EQ.1) THEN
				Vi = X(2)
			ELSE
				IF(ABMN.EQ.0) THEN
					Vi = EPS
				ELSE
					Vi = X(2)*(ABMN**X(3))
				ENDIF
			ENDIF
			IF(Vi.EQ.0) THEN
				Vi = EPS
			ENDIF
			Ni = nanimals(K)
			Devi = mean(K) - means(K)
			IF(constvar.EQ.0) THEN
				Ai = (Ni-1)*var(K) + Ni*mean(K)*mean(K)
				Bi = SIGN*Ni*mean(K)
				Hi = Ai/(2*Vi) - Bi*ABMN/Vi + Ni*ABMN*ABMN/(2*Vi)
				Gi = -Ai*X(3)/(2*Vi*ABMN) + Bi*(X(3)-1)/Vi
				Gi = Gi - Ni*(X(3)-2)*ABMN/(2*Vi)
				GRAD(2) = GRAD(2) + Ni/(2*X(2)) - Hi/X(2)
				GRAD(3) = GRAD(3) + (Ni/2 - Hi)*LOG(ABMN)
				DO J = 4, nparm
					TEMP = SIGN*grads(K,J)*(Ni*X(3)/(2*ABMN) + Gi)
					GRAD(J) = GRAD(J) + TEMP
				ENDDO
			ELSE
				TEMP = Ni*Devi*Devi + (Ni-1)*var(K)
				GRAD(2) = GRAD(2) + (Ni - TEMP/Vi)/(2*Vi)
				GRAD(3) = 0.0
				DO J = 4, nparm
					GRAD(J) = GRAD(J) - grads(K,J)*Devi*Ni/Vi
				ENDDO
			ENDIF
		ENDDO

C Finally, do the work
      call donlp2

c     unscale parameters
      DO I = 1, nparm
         X(I) = X(I)*XSC(I)
      ENDDO

      BMDL = X(1)
      DO I=0, XNPARMS-1
         XPARMS2(I) = X(I+2)
      ENDDO
      XOPTITE = OPTITE
C     Is the Ith parameter fixed or stuck at a lower bound?
      DO I=3, NG+NH
         IF (GUNIT(1,I) .EQ. 1 .AND. GUNIT(2,I) .GT. 1) XBIND(GUNIT(2,I)
     $        -2) = bind(I)  
      ENDDO
      RETURN
      END
