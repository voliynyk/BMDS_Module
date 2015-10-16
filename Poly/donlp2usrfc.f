c       For probtype=2 (for finding the bmdl)
C       The problem is to minimize f(D) = D, subject to the equality
c       constraints LIK(D, beta) - TARGET = 0, a constraint on D and beta
c       based on the benchmark response level and the type of risk used,
c       and possible equality constraints on some of the betas, and
c       possible inequality constraints of the form beta_i >= 0 for the
c       remaining betas.  The vector X contains all the parameters.  X(1)
c       contains D (the benchmark dose), X(2) contains alpha, X(3) contains
c       rho, X(4)... X(polyord + 1) contain the parameters for the 
c       dose-response model: X(4) contains b_0, etc.
c
c
C       Setup of X:
C       ML Estimation:  probtype = 1
c       X(1) = lalpha
c       X(2) = rho
C       X(i) = beta[i-3]
C       Confidence Limit estimation: probtype = 2
C       X(1) = D (bmdl)
c       X(2) = lalpha
c       X(3) = rho
C       X(i) = beta[i-4]
c       ML Estimation for model A3: probtype = 3
c       X(1) = lalpha
c       X(2) = rho
c       X(i) = beta[i-3]
c       Profile Likelihood: probtype = 4
c       X(1) = lalpha
c       X(2) = rho
c       X(i) = beta[i-3]
C       
C
        SUBROUTINE SETUP0
        INCLUDE 'O8COMM.INC'
        INCLUDE 'O8FINT.INC'
        INCLUDE 'PROBLEM.INC'
        INTEGER J
        INTEGER K, I
        if (iDebug .gt. 0) then
           write(31, *) 'Problem_type= ', probtype
        endif
C       Counter for scaling if the X's are scaled by 
C       their initial values.
C    
C       NAME IS IDENT OF THE EXAMPLE/USER AND CAN BE SET AT USERS WILL
C       THE FIRST CHARACTER MUST BE ALPHABETIC.  40 CHARACTERS MAXIMUM
        IF (probtype .EQ. 1) THEN 
           WRITE (NAME, '(A,I3)') 'MLE',mlecnt
        ELSEIF (probtype .EQ. 2) THEN
           WRITE (NAME, '(A,I3)') 'BMDL',clcnt
        ELSEIF (probtype .EQ. 4) THEN
c       WRITE ONTO NAME
           IF (counter .LE. 9) THEN
              WRITE (NAME, '(A,I1)') 'LLPro0',counter
           ELSE
              WRITE (NAME,'(A,I2)') 'LLPro',counter
           ENDIF
        ELSE
           NAME = 'OTHER'
        ENDIF
C       
        IF(probtype .EQ. 3) THEN
           CALL a3fillgunit
        ELSE
           CALL polyfillgunit
        ENDIF
c       
        ANALYT=.TRUE.
        COLD=.TRUE. 
        SILENT=.TRUE.

        EPSDIF=0.00001
        PROU=10
        MEU=20
c
C       DEL0 AND TAU0: SEE BELOW
        IF ((probtype .EQ. 1) .OR. (probtype .EQ.3) .OR. 
     1    (probtype .EQ. 4)) THEN
           IF(modtype.EQ.0) THEN
              IF (flag .EQ. 0) THEN             
                 DEL0 = 1.0D-4
                 TAU0 = 1.0D0
              ELSE
                 DEL0 = 1.0D-2
                 TAU0 = 1.0D-6
              ENDIF
           ELSE
              DEL0 = 1.0D-2
              TAU0 = 1.0D-6
           ENDIF
        ELSE
           IF (modtype .EQ. 0) THEN
              IF (flag .EQ. 0) THEN
                 DEL0 = 2.0D0
                 TAU0 = 2.0D0
              ELSE
                 DEL0 = 1.0D1
                 TAU0 = 1.0D0
              ENDIF
           ELSE
              DEL0 = 1.0D1
              TAU0 = 1.0D0
           ENDIF
        ENDIF
        
C       GCONST-ARRAY:
        DO J=0,NG+NH
           GCONST(J)=.FALSE.
C       IF THE J-TH FUNCTION IS AFFINE LINEAR
        ENDDO

c       This scales the parameters internally 
c       IF ((mlecnt.EQ.1).OR.((clcnt.GE.1).AND.(clcnt.LE.48))) THEN
c       IF ((mlecnt.EQ.1).OR.(clcnt.GT.48)) THEN
c       IF ((clcnt.GE.1).AND.(clcnt.LE.48)) THEN
c       IF (clcnt.GE.1) THEN
        IF ((mlecnt .EQ. 1).OR.(clcnt .GE. 1))THEN
          DO K = 1, nparm
             IF (X(K) .NE. 0.) THEN
                XSC(K) = ABS(X(K))
             ELSE
                XSC(K) = 1.0D0
             ENDIF
          ENDDO
        ENDIF
      if (iDebug .gt. 0) then
         write(31, *) 'I ','GUNIT (1,I) ','GUNIT(2,I) ','GUNIT(3,I)'
         write(31, *) 'Equality constraints:'
         DO I=1,NH
            write(31, *) I,GUNIT(1,I),GUNIT(2,I),GUNIT(3,I)
         ENDDO
         write(31, *) 'Inequality constraints:'
         DO I=1,NG
            write(31, *) I,GUNIT(1,I+NH),GUNIT(2,I+NH),GUNIT(3,I+NH)
         ENDDO
      endif
C     End of debug
        RETURN
        END

C       OBJECTIVE FUNCTION
        SUBROUTINE EF(X,FX)
        INCLUDE 'O8FUCO.INC'
        INCLUDE 'PROBLEM.INC'
        DOUBLE PRECISION LK,X(*),FX
        INTEGER K, J

        ICF=ICF+1
        
        IF ((probtype .EQ. 1) .OR. (probtype.EQ.3) .OR. 
     1    (probtype .EQ. 4)) THEN
C       
           IF(probtype.EQ.1 .OR. probtype.EQ.4) THEN
              CALL polymean(X)
           ELSE
              CALL a3mean(X)
           ENDIF
C       
           CALL NegLogLike(X,LK)
           FX = LK
C
        ELSEIF (probtype .EQ. 2) THEN
C       Lower Confidence Limit
C       Find the smallest D less than MLE s.t. likelihood and BMR constraints hold
           FX = X(1)
C
        ELSE
C       Upper Confidence Limit
           FX = -X(1)
        ENDIF
        
      if (iDebug .gt. 0) then
C         write(31, *) "EF() Results"
         write(31, 900) (X(J),J=1,N)
         write(31, 901) FX
 900     format ('X()= ',(10(D24.18,1X)))
 901     format ('FX= ', D24.18)
      endif
        RETURN
        END

C       GRADIENT OF OBJECTIVE FUNCTION
        SUBROUTINE EGRADF(X,GRADF)
        INCLUDE 'O8FUCO.INC'
        INCLUDE 'PROBLEM.INC'
        DOUBLE PRECISION X(*),GRADF(*)
        INTEGER J
C       USER DECLARATIONS, IF ANY ,FOLLOW
        
        ICGF=ICGF+1
c       
        DO J = 1, nparm
           GRADF(J) = 0.0D0
        ENDDO
        IF(probtype .EQ. 1 .OR. probtype .EQ. 3 .OR. 
     1    probtype .EQ. 4) THEN
           IF (probtype.EQ.1 .OR. probtype.EQ.4) THEN
c       ML Estimation and Profile Calculation
              CALL polymean(X)
              CALL polypart(X)

           ELSE
c       ML Estimation A3
              CALL a3mean(X)
              CALL a3part(X)
           ENDIF
           CALL DNegLogLike(X,GRADF) 

        ELSEIF (probtype .EQ. 2) THEN
C       Lower Confidence Limit
           GRADF(1) = 1.D0

        ELSE
C       Upper Confidence Limit
        ENDIF
        RETURN
        END
C       
C       COMPUTE THE I-TH EQUALITY CONSTRAINT, VALUE IS HXI
        SUBROUTINE EH(I,X,HXI)
        INCLUDE 'O8FUCO.INC'
        INCLUDE 'PROBLEM.INC'
        DOUBLE PRECISION X(*),HXI,LK,bmrval
        INTEGER I, J
        
        CRES(I)=CRES(I)+1
        
        IF (probtype .EQ. 1) THEN
C       ML Estimation and Profile
              HXI = X(GUNIT(2,I)) - parmval(GUNIT(2,I) - 1)
        ELSEIF (probtype .EQ. 2) THEN
C       Confidence Limits
           IF (I .EQ. 1) THEN
c       BMR-derived constraint
              CALL polymean(X)
              CALL BMRComp(X,risktype,bmrval)
              HXI = bmrval - bmr
           ELSE
c       Parameter fixed to a given value
              HXI=X(GUNIT(2,I)) - parmval(GUNIT(2,I)-2)
           ENDIF
        ELSEIF (probtype .EQ. 3) THEN
C       ML Estimation and Profile
           HXI = X(GUNIT(2,I)) - parmval(GUNIT(2,I) - 1)
        ELSEIF (probtype .EQ. 4) THEN
C       ML Estimation and Profile
           IF (I .EQ. NH) THEN
c       The last equality constriant is for profile calculation 
                 CALL polymean(X)
                 CALL BMRComp(X,risktype,bmrval)
                 HXI = bmrval - bmr
           ELSE 
              HXI = X(GUNIT(2,I)) - parmval(GUNIT(2,I) - 1)
           ENDIF        
        ENDIF

      if (iDebug .gt. 0) then
C         write(31, *) "EH() Results"
         write(31, 910) i, Gunit(2,I), HXI
C         write(31, 920) (X(J),J=1,N)
 910     format ('I=', I3, ' G(2,I)= ', I3, ' HXI= ', D24.18)
 920     format ('X()='(10(D24.18,1X)))
      endif
        RETURN
        END
C       
C       COMPUTE THE GRADIENT OF THE I-TH EQUALITY CONSTRAINT
        SUBROUTINE EGRADH(I,X,GRADHI)
        INCLUDE 'O8FUCO.INC'
        INCLUDE 'PROBLEM.INC'
        DOUBLE PRECISION X(*),GRADHI(*)
        INTEGER I,J,K
        
        IF (GUNIT(1,I) .NE. 1) CGRES(I)=CGRES(I)+1
        DO  J=1,NX
           GRADHI(J)=0.D0
        ENDDO
C  Only two gradients need to be explicitly calculated for the equality
C  constraints: when probtype .EQ. 2 and I .EQ. 1, the BMR constraint
C  and when probtype .EQ. 4 and I .EQ. NH, the BMR constraint.  For all the
C  rest, the constraints are a linear function of a single parameter, and
C  are handled in GUNIT()
        IF (probtype .EQ. 2 ) THEN
           IF (I .EQ. 1) THEN
              CALL polymean(X)
              CALL polypart(X)
c       BMD/BMR constraint
              CALL DBMRComp(X,risktype,GRADHI)
           ENDIF
        ELSEIF (probtype .EQ. 4) THEN
           IF (I .EQ. NH) THEN
              CALL polymean(X)
              call polypart(X)
c       BMD/BMR constraint
              CALL DBMRComp(X,risktype,GRADHI)
           ENDIF                    
        ENDIF
      if (iDebug .gt. 0) then
C         write(31, *) "EGRADH() Results"
C         write(31, 920) I, (X(J),J=1,N)
         write(31, 921) (GRADHI(J),J=1,N)
 920     format ('I=', I3, ' X()= '(10(D24.18,1X)))
 921     format ('GRADHI()= '(10(D24.18,1X)))
      endif
        RETURN
        END
C
C COMPUTE THE I-TH INEQUALITY CONSTRAINT, BOUNDS INCLUDED
        SUBROUTINE EG(I,X,GXI)
        INCLUDE 'O8FUCO.INC'
        INCLUDE 'PROBLEM.INC'
        DOUBLE PRECISION X(*)
        DOUBLE PRECISION GXI        
        INTEGER I
        
        CALL PolyCompIneq(X,GXI,I)

      if (iDebug .gt. 0) then
C         write(31, *) "EG() Results"
         write(31, 930) i, Gunit(2,I+NH), GXI
 930     format ('I=', I3, ' G(2,I+NH)= ', I3, ' GXI = ', D24.18)
      endif
        RETURN
        END
C
C       COMPUTE THE GRADIENT OF THE I-TH INEQUALITY CONSTRAINT
C       NOT NECESSARY FOR BOUNDS, BUT CONSTANT GRADIENTS MUST BE SET
C       HERE E.G. USING DCOPY FROM A DATA-FIELD
        SUBROUTINE EGRADG(I,X,GRADGI)
        INCLUDE 'O8FUCO.INC'
        INCLUDE 'PROBLEM.INC'
        DOUBLE PRECISION X(*) ,GRADGI(*)
        INTEGER I,J,K
        
        IF ( GUNIT(1,I+NH) .NE. 1 ) CGRES(I+NH)=CGRES(I+NH)+1
        
        DO  J=1,NX
           GRADGI(J)=0.D0
        ENDDO
        
        IF (probtype.EQ.3) THEN
C       ML Estimation
C       Do nothing; taken care of by GUNIT
        ELSEIF (probtype .EQ. 1 .OR. probtype .EQ. 4) THEN
           CALL PolyCompIneqGrad(X,GRADGI,I)
        ELSEIF (probtype .EQ. 2) THEN
           CALL PolyCompIneqGrad(X,GRADGI,I)
        ENDIF
      if (iDebug .gt. 0) then
C         write(31, *) "EGRADG() Results"
C         write(31, 940) I, probtype
         write(31, 941) (GRADGI(J),J=1,N)
 940     format ('I=',I3,' probtype=',I2)
 941     format ('GRADGI()= '(10(D24.18,1X)))
      endif
        RETURN
        END
C
C
      SUBROUTINE SETUP
      INCLUDE 'O8COMM.INC'
      INCLUDE 'PROBLEM.INC'
C     CHANGE TERMINATION CRITERION
C       EPSX = 1.0D-6
C     CHANGE I/O-CONTROL
C       INTAKT = .TRUE.
        TE0=.FALSE.
        TE1=.FALSE.
        TE2=.FALSE.
        TE3=.FALSE.
C*** NOW YOU GET FOR EVERY STEP A ONE-LINE-OUTPUT ON STDOUT
C*** IF ANY OF THE TE'S ARE SET TO TRUE
C      ......
      RETURN
      END
