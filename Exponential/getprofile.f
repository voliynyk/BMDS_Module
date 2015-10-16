C     GETPROFILE: Compute the maximum likelihood for a given BMD for 
C     the power model using donlp2. 
c     INPUT:
c         XPARMS(*)   Start values for model parameters (0:polyord) (DOUBLE)
c         XBMD             A GIVEN BMD
c         XMODEL                0 if Polynomial
c                               1 if Power
c                               2 if Hill
c         XFLAG         flag used for getprofile is the same as that for getmle
c         XITER         iteration number used to control flag                   
c     OUTPUT:
c         XOPTITE     termination code (INTEGER)
c         XPARMS2     (0:XPOLYORD) parameter estimates
C         XLL         value of the reduced log likelihood function at the maximum
C
      SUBROUTINE GETPROFILE(XPARMS,XBMD,XLL,XOPTITE,XITER,XBMRTYPE,XBMR)

      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER XOPTITE, XITER, XBMRTYPE, I
      DOUBLE PRECISION XPARMS(0:5), XBMD, XLL, XBMR
c      OPEN(1,FILE='profile_info.log',STATUS='UNKNOWN')

C -------------------------------------
      mlecnt = 9999
      clcnt = 9999
      procnt = procnt + 1
     
C     Define the problem type
      probtype = 4
      modtype = 1

      nparm = 5
      flag = 0
      bmr = XBMR
      risktype = XBMRTYPE
      N = nparm
c      bmd = XBMD
      bmd = XBMD/maxdose 
c     For right now try to optimize without dose being scaled
c      DO I = 1, ndoses
c         doses(I) = doses(I)*maxdose
c      ENDDO

      DO I = 1, nparm
         X(I) = XPARMS(I)
      ENDDO

      X(4) = XPARMS(4)*(maxdose**XPARMS(5))

c      WRITE(1,*) "procnt = ",procnt
c      WRITE(1,*) "probtype = ",probtype
c      WRITE(1,*) "modtype = ",modtype
c      WRITE(1,*) "nparm = ",nparm
c      WRITE(1,*) "flag = ",flag
c      WRITE(1,*) "bmr = ",bmr
c      WRITE(1,*) "risktype = ",risktype
c      WRITE(1,*) "N = ",N
c      WRITE(1,*) "bmd = ",bmd
c      WRITE(1,*) "adverse = ",adverse
c      WRITE(1,*) " "

c     Finally, do the work      
      call donlp2
      
      DO I = 1, nparm
         X(I) = X(I)*XSC(I)
      ENDDO
      
c     Passing back newly estimated parameters
      DO I=1, nparm
         XPARMS(I) = X(I)
      ENDDO
      XPARMS(4) = X(4)/(maxdose**X(5))

C      PRINT *,"optite=",optite

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
