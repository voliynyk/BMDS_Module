
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
      INTEGER I,J,counter
C
C   NAME IS IDENT OF THE EXAMPLE/USER AND CAN BE SET AT USERS WILL
C   THE FIRST CHARACTER MUST BE ALPHABETIC.  40 CHARACTERS MAXIMUM
      IF (probtype .EQ. 1) THEN 
         IF (flag .LE. 9) THEN
            WRITE (NAME, '(A,I1)') 'MLE0',flag
         ELSE
            WRITE (NAME, '(A,I2)') 'MLE',flag
         ENDIF
      ELSEIF (probtype .EQ. 2) THEN
         NAME = 'BMDL'
      ELSEIF (probtype .EQ. 3) THEN
         NAME = 'MLEA3'
      ELSEIF (probtype .EQ. 4) THEN
         IF (counter .LE. 9) THEN
            WRITE (NAME, '(A,I1)') 'LLPro0',counter
         ELSE
            WRITE (NAME,'(A,I2)') 'LLPro',counter
         ENDIF
      ENDIF	
c
C     Scale by the starting values
      IF ((mlecnt.EQ.1) .OR. ((clcnt.GE.1).AND.(clcnt.LE.36)))THEN
         DO I = 1, nparm
            IF (X(I) .GT. 0.) THEN
               XSC(I) = ABS(X(I))
            ELSE
               XSC(I) = 1.0
            ENDIF
         ENDDO
      ENDIF
C
C     FILL GUNIT ARRAYS
      IF(probtype.EQ.3) THEN
         CALL a3fillgunit
      ELSE		
         CALL hillfillgunit
      ENDIF
c
      ANALYT=.TRUE.
      COLD=.TRUE. 

C     THIS CONTROLS THE OUTPUT OF A LOG FILE FOR DONLP2
      SILENT=.TRUE.
      EPSDIF=0.00001
      PROU=10
      MEU=20
C
C     DEL0 AND TAU0: SEE BELOW
      IF (probtype .EQ. 1. OR. probtype. EQ. 3) THEN
          DEL0 = 1.0D-4
          TAU0 = 1.0D0
	 
      ELSEIF (flag .EQ. 0) THEN
         DEL0 = 2.0D0
         TAU0 = 2.0D0
      ELSE
         DEL0 = 1.0D1
         TAU0 = 1.0D0
      ENDIF
C
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
      INTEGER I, J, K
c     
      ICF=ICF+1
      IF (probtype .EQ. 1 .OR. probtype.EQ.3) THEN
C     Maximum Likelihood Estimation
         IF(probtype.EQ.1) THEN
            CALL hillmean(X)
         ELSE
	    CALL a3mean(X)
         ENDIF
C     
         CALL NegLogLike(X,LK)
         FX = LK
C     
      ELSEIF (probtype .EQ. 2) THEN
C     Lower Confidence Limit
C     Find the smallest D less than MLE s.t. likelihood and BMR constraints hold
         FX = X(1)
C         PRINT *, 'BMDL: ',FX
      ELSE
C     Upper Confidence Limit
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
      DO J = 1, nparm
         GRADF(J) = 0.0D0
      ENDDO
      IF (probtype .EQ. 1 . OR. probtype. EQ. 3) THEN
C     Maximum Likelihood Estimation
         IF(probtype.EQ.1) THEN
            CALL hillmean(X)
            CALL hillpart(X)
         ELSE
            CALL a3mean(X)
            CALL a3part(X)
         ENDIF     
         CALL DNegLogLike(X,GRADF)
C
      ELSEIF (probtype .EQ. 2) THEN
C     Lower Confidence Limit
         GRADF(1) = 1.D0
C
      ELSE
C     Upper Confidence Limit
         GRADF(1) = -1.D0
      ENDIF
      RETURN
      END


C     COMPUTE THE I-TH EQUALITY CONSTRAINT, VALUE IS HXI
      SUBROUTINE EH(I,X,HXI)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      DOUBLE PRECISION X(*),HXI,LK,bmrval
      INTEGER I
      CRES(I)=CRES(I)+1
      IF (probtype .EQ. 1 .OR. probtype .EQ. 3) THEN
C     ML Estimation
         HXI = X(GUNIT(2,I)) - parmval(GUNIT(2,I) - 1)
      ELSE
C     Confidence Limits
         CALL hillmean(X)
         IF (I .EQ. 1) THEN
c     BMR-derived constraint
            CALL BMRComp(X,risktype,bmrval)
            HXI = bmrval - bmr
C           PRINT *, 'delta BMR', HXI
         ELSE
c     Parameter fixed to a given value
            HXI=X(GUNIT(2,I)) - parmval(GUNIT(2,I)-2)
         ENDIF
      ENDIF
      RETURN
      END


C     COMPUTE THE GRADIENT OF THE I-TH EQUALITY CONSTRAINT
      SUBROUTINE EGRADH(I,X,GRADHI)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      DOUBLE PRECISION X(*),GRADHI(*)
      INTEGER I,J
      IF ( GUNIT(1,I) .NE. 1 ) CGRES(I)=CGRES(I)+1
      DO  J=1,NX
         GRADHI(J)=0.D0
      ENDDO
C     
      IF (probtype .EQ. 1. OR. probtype.EQ.3) THEN
C     ML Estimation
C     Do nothing, since these are taken care of in GUNIT
      ELSE
         CALL hillmean(X)
         CALL hillpart(X)
         IF (I .EQ. 1) THEN
c     BMD/BMR constraint
            CALL DBMRComp(X,risktype,GRADHI)
         ELSE
c     Individual parameter constraint
            GRADHI(GUNIT(2,I))=1.0D0
         ENDIF
      ENDIF
      RETURN
      END


C     COMPUTE THE I-TH INEQUALITY CONSTRAINT, BOUNDS INCLUDED
      SUBROUTINE EG(I,X,GXI)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      DOUBLE PRECISION X(*)
      DOUBLE PRECISION GXI	    
      INTEGER I
      
      CALL HillCompIneq(X,GXI,I)
	
      RETURN
      END


C     COMPUTE THE GRADIENT OF THE I-TH INEQUALITY CONSTRAINT
C     NOT NECESSARY FOR BOUNDS, BUT CONSTANT GRADIENTS MUST BE SET
C     HERE E.G. USING DCOPY FROM A DATA-FIELD
      SUBROUTINE EGRADG(I,X,GRADGI)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      DOUBLE PRECISION X(*) ,GRADGI(*)
      INTEGER I,J
      IF ( GUNIT(1,I+NH) .NE. 1 ) CGRES(I+NH)=CGRES(I+NH)+1
      DO  J=1,NX
         GRADGI(J)=0.D0
      ENDDO
      CALL HillCompIneqGrad(X,GRADGI,I)
      RETURN
      END


      SUBROUTINE SETUP
      INCLUDE 'O8COMM.INC'
      INCLUDE 'O8CONS.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER I
C     CHANGE TERMINATION CRITERION
C      EPSX=TM6
C       EPSX = 1.0D-300
C       DELMIN = 1.0D-300
C       SMALLW = 1.0D-300
C      TAU0 = 0.1
C     CHANGE I/O-CONTROL
C     INTAKT = .TRUE.
C     TE0=.TRUE.
C     TE1=.TRUE.
C     TE2=.TRUE.
C     TE3=.TRUE.
C***  NOW YOU GET FOR EVERY STEP A ONE-LINE-OUTPUT ON STDOUT
C     ......
      RETURN
      END
