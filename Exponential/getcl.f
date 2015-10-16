C     GETCL : Computes the confidence limits on the BMD for the
c     multistage model using donlp2.
c     INPUT:
C     XWHICH      1 = lower confidence limit
C     2 = upper confidence limit
c     XNDOSES     number of dose levels (INTEGER)
c     XDOSES(*)   dose levels used      (DOUBLE)
c     XAFFECT(*)  number of affected animals at Ith dose (INTEGER)
c     XNANIM(*)   number of animals at Ith dose (INTEGER)
c     XNPARMS     number of parameters in the model (INTEGER)
c     XBMR        Benchmark Response level requested, 
c     0 < XBMR < 1 (DOUBLE)
c     XBMD        ML estimate of BMD (DOUBLE)
c     XTARG       Target likelihood value (MaxLL - X2/2) (DOUBLE)
c     XPARMS(*)   ML estimates of model parameters (0:polyord) (DOUBLE)
c     XFIXD(*)    1 if I'th (0:polyord) parameter is specified
c     0 otherwise (INTEGER)
c     XVAL(*)     Value of fixed parameter (0:polyord) (DOUBLE)
c     XRISK       0 for absolute risk 1
C     1 for standard deviations from control mean risk
c     2 for relative risk
c     3 for point risk
c     4 for extra risk (Hill model ONLY)
c     XRESTR      1 if parameters are restricted to non-negative values
c     0 otherwise 
c     XLOGNORM    1 if parameters are for lognormal computation
c     0 otherwise 
c     OUTPUT:
c     XOPTITE     termination code (INTEGER)
c     BMDL        estimate of BMDL
c     XPARMS2     (0:XNPARMS) parameter estimates at BMDL
c     XADVERSE    1 if larger
c     -1 if smaller
c     XMODEL      0 if Polynomial
c     1 if Power
c     2 if Hill
C     
C     Get  the lower confidence limit
      SUBROUTINE GETCL(XWHICH,XNDOSES,XDOSES,XMEAN,XNANIM,XVAR,XNPARMS
     $     ,XBMR,XBMD,XTARG,XPARMS,XFIXD,XVAL,XRISK,XRESTR,BMDL,XPARMS2
     $     ,XOPTITE,XNRESM, XBIND, XADVERSE, XMODEL, XFLAG, XLOGNORM)
      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER XWHICH,XNDOSES,XNANIM(*),XNPARMS
     $     ,XFIXD(0:XNPARMS),XRISK,XRESTR,XOPTITE,XNRESM,
     $     XBIND(0:XNPARMS),XADVERSE, XMODEL, XFLAG, XLOGNORM
      DOUBLE PRECISION XDOSES(*),XBMR,XBMD,XTARG,XPARMS(0:XNPARMS)
     $     ,XVAL(0:XNPARMS), BMDL, XPARMS2(0:XNPARMS),XMEAN(*),
     $     XVAR(*)
      INTEGER I
c     c      REAL SIGN,Vi,ABMN,Ai,Ni,Bi,Gi,Hi,TEMP,Devi,GRAD(23),EPS
      REAL EPS
      EPS = .00000001
C     -------------------------------------
      mlecnt = -9999
      clcnt = clcnt + 1

C     Define the problem type
      probtype = 2
      modtype = XMODEL
C     Load values into common blocks
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
      lognorm = XLOGNORM
c     
ccccccrescale
c     maxdose = doses(1)
c     DO I = 2, ndoses
c     IF (doses(I) .GT. maxdose) maxdose = doses(I)
c     ENDDO
c     DO I = 1, ndoses
c     doses(I) = doses(I)/maxdose
c     ENDDO
c     bmd = bmd/maxdose
c     X(1) = X(1)/maxdose
c     X(5) = X(5)*(maxdose**X(6))
c     IF (parmfixd(4) .EQ. 1) parmval(4) = parmval(4)*(maxdose**X(6))
CCCCC 
c     c CALL POWMEAN(X)
      
c     
c     Get gradients for scaling factors
c     
      
c     c DO J = 1,nparm
c     c         GRAD(J) = 0.D0
c     c ENDDO
c     c         DO K = 1, ndoses
c     c                 IF(means(K).LT.0) THEN
c     c                   SIGN = -1
c     c                 ELSE
c     c                   SIGN = 1
c     c                 ENDIF
c     c                 ABMN = ABS(means(K))
c     c                 IF(constvar.EQ.1) THEN
c     c                         Vi = X(2)
c     c                 ELSE
c     c                         IF(ABMN.EQ.0) THEN
c     c                                 Vi = EPS
c     c                         ELSE
c     c                                 Vi = X(2)*(ABMN**X(3))
c     c                         ENDIF
c     c                 ENDIF
c     c                 IF(Vi.EQ.0) THEN
c     c                         Vi = EPS
c     c                 ENDIF
c     c                 Ni = nanimals(K)
c     c                 Devi = mean(K) - means(K)
c     c                 IF(constvar.EQ.0) THEN
c     c                         Ai = (Ni-1)*var(K) + Ni*mean(K)*mean(K)
c     c                         Bi = SIGN*Ni*mean(K)
c     c                         Hi = Ai/(2*Vi) - Bi*ABMN/Vi + Ni*ABMN*ABMN/(2*Vi)
c     c                         Gi = -Ai*X(3)/(2*Vi*ABMN) + Bi*(X(3)-1)/Vi
c     c                         Gi = Gi - Ni*(X(3)-2)*ABMN/(2*Vi)
c     c                         GRAD(2) = GRAD(2) + Ni/(2*X(2)) - Hi/X(2)
c     print*, 1
c     c                         GRAD(3) = GRAD(3) + (Ni/2 - Hi)*LOG(ABMN+EPS)
c     print*, 2
c     c                         DO J = 4, nparm
c     c                                 TEMP = SIGN*grads(K,J)*(Ni*X(3)/(2*ABMN) + Gi)
c     c                                 GRAD(J) = GRAD(J) + TEMP
c     c                         ENDDO
c     c                 ELSE
c     c                         TEMP = Ni*Devi*Devi + (Ni-1)*var(K)
c     c                         GRAD(2) = GRAD(2) + (Ni - TEMP/Vi)/(2*Vi)
c     c                         GRAD(3) = 0.0
c     c                         DO J = 4, nparm
c     c                                 GRAD(J) = GRAD(J) - grads(K,J)*Devi*Ni/Vi
c     c                         ENDDO
c     c                 ENDIF
c     c         ENDDO
c     
      
c     
      DO I = 1, N
c     c         IF (GRAD(I) .EQ. 0.) THEN
         IF (X(I) .EQ. 0.) THEN
            XSC(I) = 1.0
         ELSE
            XSC(I) = X(I)
         ENDIF
      ENDDO
      
c     
C     Finally, do the work
      call donlp2
c     Unscale by maxdose
c     c      X(1) = X(1)*maxdose
c     c
c     c      X(5) = X(5)/(maxdose**X(6))
      
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
c     DO I=3, NG+NH
c     IF (GUNIT(1,I) .EQ. 1 .AND. GUNIT(2,I) .GT. 1) XBIND(GUNIT(2,I)
c     $        -2) = bind(I)  
c     ENDDO
      RETURN
      END
