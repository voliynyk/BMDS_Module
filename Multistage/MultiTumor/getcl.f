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
      SUBROUTINE GETCL(XWHICH,XBMR,XBMD,XTARG,XPARMS,XFIXD,XVAL,XRISK,
     $     BMDL,XPARMS2,XOPTITE,XNRESM, XBIND, XBMDbnd)
      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER XWHICH,XFIXD(0:polyord),XRISK,XOPTITE,XNRESM,
     $     XBIND(0:polyord), XBMDbnd
      DOUBLE PRECISION XBMR,XBMD,XTARG,XPARMS(0:polyord),
     $     XVAL(0:polyord),BMDL,XPARMS2(0:polyord),temp1
      INTEGER I, temp2, temp3
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
      RETURN
      END
c       *** Modified by gln 01/28/2009, to impliment multiple tumor      
C
C       *** Added by cvl 8/2007 - start block
C     GETCL2 : Computes the confidence limits on the BMD for two combined 
c     endpoints fit to a multistage model using donlp2.
c     INPUT:
C         XWHICH      4 = combined lower confidence limit
c         XBMR        Benchmark Response level requested, 
c                     0 < XBMR < 1 (DOUBLE)
c         XBMD        ML estimate of BMD (DOUBLE)
c         XTARG       Target likelihood value (MaxLL - X2/2) (DOUBLE)
c         XPARMS(*)   ML estimates of model parameters (0:polyord) (DOUBLE)
c         AXPARMS(*)  ML estimates of model parameters (0:Apolyord) (DOUBLE)
c         XRISK       1 for extra risk (only exrisk allowed)
c         XRESTR      1 if parameters are restricted to non-negative values
c                     0 otherwise 
c     OUTPUT:
c         XOPTITE     termination code (INTEGER)
c         BMDL        estimate of BMDL
c         XPARMS2     (0:XPOLYORD+Apolyord+1) parameter estimates at BMDL
C
c      SUBROUTINE GETCL2(XWHICH,XBMR,XBMD,XTARG,XPARMS,AXPARMS,
c     $   XFIXD,XVAL,AXFIXD, AXVAL, XRISK,BMDL,XPARMS2,XOPTITE,XNRESM,
c     $    XBIND, XBMDbnd,CB,xlk,axmax,xmax,scale)
c      INCLUDE 'O8COMM.INC'
c      INCLUDE 'O8FINT.INC'
c      INCLUDE 'PROBLEM.INC'
c      INTEGER XWHICH,XFIXD(0:polyord),XRISK,XOPTITE,XNRESM,
c     $     XBIND(0:polyord), XBMDbnd,CB,
c     $     AXBIND(0:Apolyord), AXFIXD(0:Apolyord)
c      DOUBLE PRECISION XBMR,XBMD,XTARG,XPARMS(0:polyord),
c     $     XVAL(0:polyord),BMDL,XPARMS2(0:CB),temp1,
c     $     AXPARMS(0:Apolyord),AXVAL(0:polyord),large
c      Double Precision probs(25),SLOGF,sum,sum2,xlk,axmax,xmax,scale
c      INTEGER I, temp2, temp3, K, J
c      REAL maxdose
c
c       *** Start of modification gln 01/28/2009, to impliment multiple tumor

      SUBROUTINE GETCLMT(XWHICH,lnMaxNParms,XBMR,XBMD,XTARG,XPARMS,
     $     XFIXD,XVAL,XRISK,BMDL,XPARMS2,
     $     XOPTITE,XNRESM,XBIND,XBMDbnd)
      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
      INCLUDE 'MULTITUMOR.INC'
      INTEGER XWHICH,lnMaxNParms,XFIXD(0:nT-1,0:lnMaxNParms-1),XRISK,
     $     XOPTITE,XNRESM,XBIND(0:239),XBMDbnd
      DOUBLE PRECISION XBMR,XBMD,XTARG,XPARMS(0:nT-1,0:lnMaxNParms-1),
     $     XVAL(0:nT-1,0:lnMaxNParms-1),
     $     BMDL,XPARMS2(0:nT-1,0:lnMaxNParms-1)
      INTEGER I, temp2, temp3, L, K, J, ioffset
      DOUBLE PRECISION SLOGF, temp1, xlk, SUM2, SUM, 
     $     PROBS(1:nMaxObs,0:nT-1),correction

c       *** End of modification gln 01/28/2009, to impliment multiple tumor
C -------------------------------------
C Define the problem type
c      print *,'Inside getcl2'
      probtype = XWHICH + 1
      print *,'probtype = ', probtype
C  BCA added 1/09 start *******
C     T must be initialized or passed to be the number of tumors, replaced by nT (gln-01/28/2009)
      N = 0
      DO J=0,nT-1
	     N = N + arnPoly(J)
      ENDDO
      N= N + nT + 1
C  BCA added 1/09 end *******
C      N = total length of parameter vector X that starts with X(1) being the log(BMDL)
      print *,'N= ', N
      bmr=XBMR
C Work with log(bmd)
      print *,'Xbmd= ', XBMD
      bmd=LOG(XBMD)
      print *,'XTARG= ', XTARG
      target = XTARG
C   added by gln 02/04/09 start Check values of XPARMS ******
      DO I=0, nT-1
         Print *,"Tumor ",I
         DO J=0, arnPoly(I)
            print *, 'XPARMS(',I,',',J,')=',XPARMS(I,J)
         ENDDO
         print *
      ENDDO
C   added by gln 02/04/09 end Check values of XPARMS ******

C       *** Added by cvl 3/2008 - end block
C   Changed by bca 1/09 start ***********************
      X(1) = bmd
      ioffset = 0
      DO J = 0, nT-1
      if(J.gt.0) ioffset = ioffset + arnPoly(J-1)
      	DO I = 0, arnPoly(J)
           X(I+J+ioffset+2) = XPARMS(J,I)
c           X(I+J+ioffset+2) = XPARMS(J+1,I+1)  changed to above -- bca 3/17
c           X(I+J+ioffset+2) = PARMSX(J+1,I+1)
        ENDDO
      ENDDO
      DO J = 1, nT-1
         X(3) = X(3) + XPARMS(J,1)
      ENDDO
C   Changed by bca 1/09 end ***********************

	
C Finally, do the work
      print *, 'Before call to donlp2 in GETCLMT'
      DO J = 1, N
         print *, 'X(',J,') = ', X(J)
      ENDDO
      call donlp2

      BMDL = EXP(X(1))
      print *,' scaled BMDL est = ' , BMDL
C   Changed by bca 1/09 start ***********************
      ioffset = 0 
      DO J = 0, nT-1
        if (J.gt.0) ioffset = ioffset + arnPoly(J-1)
        DO I = 0, arnPoly(J)
           XPARMS2(J,I) = X(I+J+ioffset+2)
           print *, 'XPARMS2(',J,',',I,')=',XPARMS2(J,I)
        ENDDO
      ENDDO
C   Changed by bca 1/09 end ***********************      
      XOPTITE = OPTITE

C     Is the Ith parameter fixed or stuck at a lower bound?

      DO I=3, NG+NH
         IF (GUNIT(1,I) .EQ. 1 .AND. GUNIT(2,I) .GT. 1) XBIND(GUNIT(2,I)
     $        -2) = bind(I)  
      ENDDO
C     Is the BMD at its lower (for BMDL: XWHICH = 4)
      XBMDbnd = bind(NH + 2)


C       *** Added by cvl 3/2008 - start block 
C     Changed by BCA 1/09 start *************************************
C     Calculation of the combined log-lieklihood
      xlk = 0.0
      SUM2 = 0.0
	  ioffset = 0
      correction = 0.0
      DO L=1, nT-1
         ioffset = ioffset + arnPoly(L-1)
         correction = correction + X(ioffset+L+3)
      ENDDO
      DO L=nT-1,0,-1
         IF(L.lt.nT-1) ioffset = ioffset-arnPoly(L)
            DO K=1, arnObs(L)
               SUM=X(ioffset+L+arnPoly(L)+2)
	       IF(ioffset+L+arnPoly(L)+2.eq.3)
     $              SUM=X(ioffset+L+arnPoly(L)+2)-correction
               DO J = arnPoly(L)-1, 0,-1
	           IF(L.eq.0 .AND. J.eq.1) THEN
		          SUM = SUM * mdDoses(L,K) + X(J)-correction
	           ELSE
                  SUM = SUM * mdDoses(L,K) + X(J)
		       ENDIF
            ENDDO
            IF (SUM .LT. 0) SUM = 0.0D0
            PROBS(K,L) = 1 - DEXP(-SUM)
      ENDDO
            
      DO K=1,arnObs(L)
         SUM2 = SUM2+DBLE(mnAffected(L,K))*SLOGF(PROBS(K,L))+
     $     DBLE(mnAnim(L,K)-mnAffected(L,K)) * SLOGF(1.0D0- 
     $     PROBS(K,L))
      ENDDO
      xlk = xlk + SUM2
	  SUM2 = 0.0
	  ENDDO
c	end of do over L
C     Changed by BCA 1/09 end *************************************
      
C       *** Added by cvl 3/2008 - end block
      RETURN
      END
C       *** Added by cvl 8/2007 - end block




 
