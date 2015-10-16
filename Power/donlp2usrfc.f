
C     The problem is to minimize f(D) = D, subject to the equality
c     constraints LIK(D, beta) - TARGET = 0, a constraint on D and beta
c     based on the benchmark response level and the type of risk used,
c     and possible equality constraints on some of the betas, and
c     possible inequality constraints of the form beta_i >= 0 for the
c     remaining betas.  The vector X contains all the parameters.  X(1)
c     contains D (the benchmark dose), X(2) ... X(polyord + 1) contain
c     the parameters for the dose-response model: X(2) contains b_0, etc.
C     Setup of X:
C         ML Estimation:
C             X(i) = beta[i-1]
C         Confidence Limit estimation
C             X(1) = D (bmdl)
C             X(i) = beta[i-2]
C
      SUBROUTINE SETUP0
      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER J,K
C
C   NAME IS IDENT OF THE EXAMPLE/USER AND CAN BE SET AT USERS WILL
C   THE FIRST CHARACTER MUST BE ALPHABETIC.  40 CHARACTERS MAXIMUM
c     WRITE ONTO NAME
      IF (probtype .EQ. 1) THEN
c     MLE
         IF (mlecnt .LE. 9) THEN
            WRITE (NAME, '(A,I1)') 'MLE00',mlecnt
         ELSEIF (mlecnt .LE. 99) THEN
            WRITE (NAME, '(A,I2)') 'MLE0',mlecnt
         ELSEIF (mlecnt .LE. 999) THEN
            WRITE (NAME, '(A,I3)') 'MLE',mlecnt
         ELSE
            WRITE (NAME, '(A)') 'ERRORM'
         ENDIF
      ELSEIF (probtype .EQ. 2) THEN
c     BMDL
         IF (clcnt .LE. 9) THEN
            WRITE (NAME, '(A,I1)') 'BMDL00',clcnt
         ELSEIF (clcnt .LE. 99) THEN
            WRITE (NAME, '(A,I2)') 'BMDL0',clcnt
         ELSEIF (clcnt .LE. 999) THEN
            WRITE (NAME, '(A,I3)') 'BMDL',clcnt
         ELSE
            WRITE (NAME, '(A)') 'ERRORB'
         ENDIF
      ELSEIF (probtype .EQ. 3) THEN
c     MLEA3
      ELSEIF (probtype .EQ. 4) THEN
c     Profile
         IF (procnt .LE. 9) THEN
            WRITE (NAME, '(A,I1)') 'PROF00',procnt
         ELSEIF (procnt .LE. 99) THEN
            WRITE (NAME, '(A,I2)') 'PROF0',procnt
         ELSEIF (procnt .LE. 999) THEN
            WRITE (NAME, '(A,I3)') 'PROF',procnt
         ELSE
            WRITE (NAME, '(A)') 'ERRORB'
         ENDIF           
      ELSE
         WRITE (NAME, '(A)') 'OTHER'
      ENDIF
C     
C     Scale by the starting values
      IF (((probtype.EQ.1).AND.(flag.EQ.-1)) .OR. 
     $     ((clcnt.GE.1).AND.(clcnt.LE.17)))THEN
         DO K = 1, nparm
            IF (X(K) .NE. 0) THEN
               XSC(K) = ABS(X(K))
            ELSE
               XSC(K) = 1.0D0
            ENDIF
         ENDDO
      ENDIF
C
C
      IF(probtype .EQ. 3) THEN
         CALL a3fillgunit
      ELSE
         CALL POWfillgunit
      ENDIF
c     
      ANALYT=.TRUE.
      COLD=.TRUE.
      IF (probtype .EQ. 4)THEN
         SILENT=.FALSE. 
      ELSE
         SILENT=.TRUE.
      ENDIF
      EPSDIF=0.00001
C
      PROU=10
      MEU=20
C     
C     DEL0 AND TAU0: SEE BELOW
      IF (probtype .EQ. 1) THEN
C     MLE
c         IF (flag .LE. 0) THEN
            DEL0 = 1.0D-4
            TAU0 = 1.0D0
c         ELSE
c            DEL0 = 1.5D-4
c            TAU0 = 1.0D0
c         ENDIF
      ELSEIF (probtype .EQ. 2) THEN
C     BMDL
         DEL0 = 2.0
         TAU0 = 2.0
      ELSEIF (probtype .EQ. 3) THEN
C     MLEA3
         DEL0 = 1.0D0
         TAU0 = 1.0D0	 
      ELSEIF (probtype .EQ. 4) THEN
C     PROFILE LIKELIHOOD
         DEL0 = 2.0D-4
         TAU0 = 1.0D0
      ELSE
         DEL0 = 2.0
         TAU0 = 2.0
      ENDIF
      
C     GCONST-ARRAY:
      DO J=0,NG+NH
         GCONST(J)=.FALSE.
C     IF THE J-TH FUNCTION IS AFFINE LINEAR
      ENDDO
      RETURN
      END
      
C     OBJECTIVE FUNCTION
      SUBROUTINE EF(X,FX)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      DOUBLE PRECISION LK,X(*),FX
      INTEGER J
      IF((probtype.EQ.4) .AND. (procnt.EQ.1))THEN
         OPEN(11,FILE='objfunc.txt',STATUS='UNKNOWN')
      ENDIF
c	
      ICF=ICF+1
      IF((probtype.EQ.4) .AND. (procnt.EQ.1))THEN
         WRITE(11,*) NAME
         DO J = 1, nparm
            WRITE(11,*) "X(",J,")= ",X(J)
         ENDDO
      ENDIF

      IF ((probtype.EQ.1) .OR. (probtype.EQ.3) .OR. (probtype.EQ.4))
     $     THEN
C     MLE, MLEA3, or PROFILE
         IF ((probtype .EQ. 1) .OR. (probtype .EQ. 4)) THEN
            CALL POWmean(X)
         ELSE
            CALL a3mean(X)
         ENDIF
C     
         CALL NegLogLike(X,LK)
         FX = LK
c
      IF((probtype.EQ.4) .AND. (procnt.EQ.1))THEN
         DO J = 1, nparm
            WRITE(11,*) "X(",J,")= ",X(J)
         ENDDO
         WRITE(11,*) "FX=    ",FX
      ENDIF
C
      ELSEIF (probtype .EQ. 2) THEN
C     Lower Confidence Limit
C     Find the smallest D less than MLE s.t. likelihood and BMR constraints hold
	   
	   FX = X(1)

      ELSE
C Upper Confidence Limit
         FX = -X(1)
      ENDIF
      RETURN
      END

C     GRADIENT OF OBJECTIVE FUNCTION
      SUBROUTINE EGRADF(X,GRADF)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      DOUBLE PRECISION X(*),GRADF(*)
      INTEGER J
C     USER DECLARATIONS, IF ANY, FOLLOW
      ICGF=ICGF+1
      IF((probtype.EQ.4) .AND. (procnt.EQ.1))THEN
         OPEN(12,FILE='gradobjfunc.txt',STATUS='UNKNOWN')
      ENDIF

      DO J = 1, nparm
         GRADF(J) = 0.0D0
      ENDDO

      IF((probtype.EQ.4) .AND. (procnt.EQ.1))THEN
         WRITE(12,*) NAME
         DO J = 1, nparm
            WRITE(12,*) "X(",J,")= ",X(J)
         ENDDO
      ENDIF

      IF ((probtype.EQ.1) .OR. (probtype.EQ.3) .OR. (probtype.EQ.4))
     $     THEN
c     MLE, MLEA3, or PROFILE
         IF ((probtype .EQ. 1) .OR. (probtype .EQ. 4))	THEN
            CALL POWmean(X)
            CALL POWpart(X)
         ELSE
            CALL a3mean(X)
            CALL a3part(X)
         ENDIF
         CALL DNegLogLike(X,GRADF)

c****************************************************
c         GRADF(5) = -GRADF(5)
c****************************************************
C     
      ELSEIF (probtype .EQ. 2) THEN
C     Lower Confidence Limit
         GRADF(1) = 1.D0
      ELSE
C     Upper Confidence Limit
C        GRADF(1) = -1.D0
      ENDIF

      IF((probtype.EQ.4) .AND. (procnt.EQ.1))THEN
         DO J = 1, nparm
            WRITE(12,*) "X(",J,")= ",X(J)
         ENDDO
         DO J = 1, nparm
            WRITE(12,*) "GRADF",J,"= ",GRADF(J)
         ENDDO
      ENDIF

      RETURN
      END

C  COMPUTE THE I-TH EQUALITY CONSTAINT, VALUE IS HXI
      SUBROUTINE EH(I,X,HXI)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      DOUBLE PRECISION X(*),HXI,LK,bmrval
      INTEGER I,J
      IF((probtype.EQ.4) .AND. (procnt.EQ.1))THEN
         OPEN(13,FILE='eqconst.txt',STATUS='UNKNOWN')
      ENDIF

      IF((probtype.EQ.4) .AND. (procnt.EQ.1))THEN
         WRITE(13,*) NAME
         DO J = 1, nparm
            WRITE(13,*) "X(",J,")= ",X(J)
         ENDDO
         WRITE(13,*) "X(4)*bmd**X(5) = ",X(4)*bmd**X(5)
      ENDIF

      CRES(I)=CRES(I)+1
      IF (probtype .EQ. 1 .OR. probtype .EQ. 3) THEN
C     ML Estimation
         HXI = X(GUNIT(2,I)) - parmval(GUNIT(2,I) - 1)
      ELSEIF (probtype .EQ. 4) THEN
c     Profile
         IF (I .LT. NH) THEN
            HXI = X(GUNIT(2,I)) - parmval(GUNIT(2,I) - 1)
         ELSE
c           Negative linear trend or slope
            IF (adverse .LT. 0) THEN
c              Absolute Deviation
               IF (risktype .EQ. 0) THEN
                  HXI = X(4)*(bmd**X(5))+bmr
c              Standard Deviation
               ELSEIF (risktype .EQ. 1) THEN
                  HXI = X(4)*(bmd**X(5))+bmr*(var(1)**(0.5))
c              Relative Deviation
               ELSEIF (risktype .EQ. 2) THEN
                  HXI = X(4)*(bmd**X(5))+bmr*X(3)
c              Point Deviation
               ELSEIF (risktype .EQ. 3) THEN
                  HXI = X(3)+X(4)*(bmd**X(5))-bmr
               ELSE
                  PRINT *,"ERROR IN RISKTYPE VALUE"
               ENDIF
c           Positive linear trend or slope
            ELSEIF (adverse .GT. 0) THEN
c              Absolute Deviation
               IF (risktype .EQ. 0) THEN
                  HXI = X(4)*(bmd**X(5))-bmr
c              Standard Deviation
               ELSEIF (risktype .EQ. 1) THEN
                  HXI = X(4)*(bmd**X(5))-bmr*(var(1)**(0.5))
c              Relative Deviation
               ELSEIF (risktype .EQ. 2) THEN
                  HXI = X(4)*(bmd**X(5))-bmr*X(3)
c              Point Deviation
               ELSEIF (risktype .EQ. 3) THEN
                  HXI = X(3)+X(4)*(bmd**X(5))-bmr
               ELSE
                  PRINT *,"ERROR IN RISKTYPE VALUE"
               ENDIF 
            ELSE
               PRINT *,"ERROR IN ADVERSE VALUE"
            ENDIF
         ENDIF
      ELSE
C     Confidence Limits
      
         CALL POWmean(X)
         IF (I .EQ. 1) THEN
c     BMR-derived constraint
            CALL BMRComp(X,risktype,bmrval)
            HXI = bmrval - bmr
         ELSE
c     Parameter fixed to a given value
            HXI=X(GUNIT(2,I)) - parmval(GUNIT(2,I)-2)
         ENDIF
      ENDIF

      IF((probtype.EQ.4) .AND. (procnt.EQ.1))THEN
         DO J = 1, nparm
            WRITE(13,*) "X(",J,")= ",X(J)
         ENDDO
         WRITE(13,*) "HXI = ",HXI
      ENDIF
      
      RETURN
      END
C
C  COMPUTE THE GRADIENT OF THE I-TH EQUALITY CONSTRAINT
      SUBROUTINE EGRADH(I,X,GRADHI)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      DOUBLE PRECISION X(*),GRADHI(*)
      INTEGER I,J
      IF((probtype.EQ.4) .AND. (procnt.EQ.1))THEN
         OPEN(14,FILE='gradeqconst.txt',STATUS='UNKNOWN')
      ENDIF
      IF((probtype.EQ.4) .AND. (procnt.EQ.1))THEN
         WRITE(14,*) NAME
         DO J = 1, nparm
            WRITE(14,*) "X(",J,")= ",X(J)
         ENDDO
      ENDIF

      IF ( GUNIT(1,I) .NE. 1 ) CGRES(I)=CGRES(I)+1
      DO  J=1,NX
	GRADHI(J)=0.D0
      ENDDO
C
      IF ((probtype .EQ. 1) .OR. (probtype .EQ. 3)) THEN
C     ML Estimation
C     Do nothing, since these are taken care of in GUNIT

      ELSEIF (probtype .EQ. 4) THEN
c     PROFILE

         IF (I .LT. NH) THEN
c     Do nothing for the first equality contraints taken
c     care of in GUNIT
         ELSEIF (I .EQ. NH) THEN
            GRADHI(1) = 0
            GRADHI(2) = 0
            GRADHI(4) = (bmd**X(5))
            GRADHI(5) = X(4)*(bmd**X(5))*LOG(bmd)
c           Negative linear trend or slope
            IF (adverse .LT. 0) THEN
c              Absolute Deviation
               IF (risktype .EQ. 0) THEN
                  GRADHI(3) = 0
c              Standard Deviation
               ELSEIF (risktype .EQ. 1) THEN
                  GRADHI(3) = 0
c              Relative Deviation
               ELSEIF (risktype .EQ. 2) THEN
                  GRADHI(3) = bmr
c              Point Deviation
               ELSEIF (risktype .EQ. 3) THEN
                  GRADHI(3) = 1.0
               ELSE
                  PRINT *,"ERROR IN RISKTYPE VALUE"
               ENDIF
c           Positive linear trend or slope
            ELSEIF (adverse .GT. 0) THEN
c              Absolute Deviation
               IF (risktype .EQ. 0) THEN
                  GRADHI(3) = 0
c              Standard Deviation
               ELSEIF (risktype .EQ. 1) THEN
                  GRADHI(3) = 0
c              Relative Deviation
               ELSEIF (risktype .EQ. 2) THEN
                  GRADHI(3) = -bmr
c              Point Deviation
               ELSEIF (risktype .EQ. 3) THEN
                  GRADHI(3) = 1.0
               ELSE
                  PRINT *,"ERROR IN RISKTYPE VALUE"
               ENDIF 
            ELSE
               PRINT *,"ERROR IN ADVERSE VALUE"
            ENDIF            
         ELSE
            PRINT *,"ERROR IN I VALUE"
         ENDIF
C
      ELSE
         CALL POWmean(X)
         CALL POWpart(X)
         IF (I .EQ. 1) THEN
c     BMD/BMR constraint
            CALL DBMRComp(X,risktype,GRADHI)
         ELSE
c     Individual parameter constraint
            GRADHI(GUNIT(2,I))=1.0D0
         ENDIF
      ENDIF

      IF((probtype.EQ.4) .AND. (procnt.EQ.1))THEN
         DO J = 1, nparm
            WRITE(14,*) "X(",J,")= ",X(J)
         ENDDO
         DO J = 1, nparm
            WRITE(14,*) "GRADHI",J,"= ",GRADHI(J)
         ENDDO
      ENDIF

      RETURN
      END

C COMPUTE THE I-TH INEQUALITY CONSTRAINT, BOUNDS INCLUDED

      SUBROUTINE EG(I,X,GXI)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      DOUBLE PRECISION X(*)
	DOUBLE PRECISION GXI	    
      INTEGER I

	CALL POWCompIneq(X,GXI,I)
      
      RETURN
      END
C COMPUTE THE GRADIENT OF THE I-TH INEQUALITY CONSTRAINT
C NOT NECESSARY FOR BOUNDS, BUT CONSTANT GRADIENTS MUST BE SET
C HERE E.G. USING DCOPY FROM A DATA-FIELD
      SUBROUTINE EGRADG(I,X,GRADGI)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      DOUBLE PRECISION X(*) ,GRADGI(*)
      INTEGER I,J
      IF ( GUNIT(1,I+NH) .NE. 1 ) CGRES(I+NH)=CGRES(I+NH)+1

      DO  J=1,NX
	GRADGI(J)=0.D0
      ENDDO

	CALL POWCompIneqGrad(X,GRADGI,I)
     
      RETURN
      END
c
      SUBROUTINE SETUP
      INCLUDE 'O8COMM.INC'
      INCLUDE 'PROBLEM.INC'
      
C     CHANGE TERMINATION CRITERION
c      EPSX = TM6
C     CHANGE I/O-CONTROL
C	INTAKT = .TRUE.
c      TE0=.TRUE.
c	 TE1=.TRUE.
c      TE2=.TRUE.
C      TE3=.TRUE.
C*** NOW YOU GET FOR EVERY STEP A ONE-LINE-OUTPUT ON STDOUT
C      ......
      RETURN
      END
