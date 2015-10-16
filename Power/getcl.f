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
      SUBROUTINE GETCL(XWHICH,XNDOSES,XDOSES,XMEAN,XNANIM,XVAR,XNPARMS
     $     ,XBMR,XBMD,XTARG,XPARMS,XFIXD,XVAL,XRISK,XRESTR,BMDL,XPARMS2
     $     ,XOPTITE,XNRESM, XBIND, XADVERSE, XMODEL, XFLAG)
      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER XWHICH,XNDOSES,XNANIM(*),XNPARMS
     $     ,XFIXD(0:XNPARMS),XRISK,XRESTR,XOPTITE,XNRESM,
     $     XBIND(0:XNPARMS),XADVERSE, XMODEL, XFLAG
      DOUBLE PRECISION XDOSES(*),XBMR,XBMD,XTARG,XPARMS(0:XNPARMS)
     $     ,XVAL(0:XNPARMS), BMDL, XPARMS2(0:XNPARMS),XMEAN(*),
     $	 XVAR(*)
      INTEGER I
cc      REAL SIGN,Vi,ABMN,Ai,Ni,Bi,Gi,Hi,TEMP,Devi,GRAD(23),EPS
	REAL EPS
	EPS = .00000001
C -------------------------------------
        mlecnt = -9999
        clcnt = clcnt + 1

C Define the problem type
      probtype = 2
	modtype = XMODEL
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
	flag = XFLAG
c
cccccc     rescale
c      maxdose = doses(1)
c	DO I = 2, ndoses
c	   IF (doses(I) .GT. maxdose) maxdose = doses(I)
c	ENDDO
c	DO I = 1, ndoses
c	   doses(I) = doses(I)/maxdose
c	ENDDO
c	bmd = bmd/maxdose
c	X(1) = X(1)/maxdose
c	X(5) = X(5)*(maxdose**X(6))
c	IF (parmfixd(4) .EQ. 1) parmval(4) = parmval(4)*(maxdose**X(6))
CCCCC	      
cc	CALL POWMEAN(X)
	
c
c    Get gradients for scaling factors
c
      
cc	DO J = 1,nparm
cc		GRAD(J) = 0.D0
cc	ENDDO
cc		DO K = 1, ndoses
cc			IF(means(K).LT.0) THEN
cc			  SIGN = -1
cc			ELSE
cc			  SIGN = 1
cc			ENDIF
cc			ABMN = ABS(means(K))
cc			IF(constvar.EQ.1) THEN
cc				Vi = X(2)
cc			ELSE
cc				IF(ABMN.EQ.0) THEN
cc					Vi = EPS
cc				ELSE
cc					Vi = X(2)*(ABMN**X(3))
cc				ENDIF
cc			ENDIF
cc			IF(Vi.EQ.0) THEN
cc				Vi = EPS
cc			ENDIF
cc			Ni = nanimals(K)
cc			Devi = mean(K) - means(K)
cc			IF(constvar.EQ.0) THEN
cc				Ai = (Ni-1)*var(K) + Ni*mean(K)*mean(K)
cc				Bi = SIGN*Ni*mean(K)
cc				Hi = Ai/(2*Vi) - Bi*ABMN/Vi + Ni*ABMN*ABMN/(2*Vi)
cc				Gi = -Ai*X(3)/(2*Vi*ABMN) + Bi*(X(3)-1)/Vi
cc				Gi = Gi - Ni*(X(3)-2)*ABMN/(2*Vi)
cc				GRAD(2) = GRAD(2) + Ni/(2*X(2)) - Hi/X(2)
c	print*, 1
cc				GRAD(3) = GRAD(3) + (Ni/2 - Hi)*LOG(ABMN+EPS)
c	print*, 2
cc				DO J = 4, nparm
cc					TEMP = SIGN*grads(K,J)*(Ni*X(3)/(2*ABMN) + Gi)
cc					GRAD(J) = GRAD(J) + TEMP
cc				ENDDO
cc			ELSE
cc				TEMP = Ni*Devi*Devi + (Ni-1)*var(K)
cc				GRAD(2) = GRAD(2) + (Ni - TEMP/Vi)/(2*Vi)
cc				GRAD(3) = 0.0
cc				DO J = 4, nparm
cc					GRAD(J) = GRAD(J) - grads(K,J)*Devi*Ni/Vi
cc				ENDDO
cc			ENDIF
cc		ENDDO
c
      
c
      DO I = 1, N
cc         IF (GRAD(I) .EQ. 0.) THEN
         IF (X(I) .EQ. 0.) THEN
            XSC(I) = 1.0
         ELSE
            XSC(I) = X(I)
         ENDIF
      ENDDO
	
c
C Finally, do the work
      call donlp2
c Unscale by maxdose
cc      X(1) = X(1)*maxdose
cc
cc      X(5) = X(5)/(maxdose**X(6))
      
c     Unscale parameters by starting values
      DO I=1, nparm
         X(I) = X(I)*XSC(I)
      ENDDO

      BMDL = X(1)
      DO I = 0, XNPARMS-1
         XPARMS2(I) = X(I+2)
      ENDDO

	
      XOPTITE = OPTITE
C     Is the Ith parameter fixed or stuck at a lower bound?
c      DO I=3, NG+NH
c         IF (GUNIT(1,I) .EQ. 1 .AND. GUNIT(2,I) .GT. 1) XBIND(GUNIT(2,I)
c     $        -2) = bind(I)  
c      ENDDO
      RETURN
      END
