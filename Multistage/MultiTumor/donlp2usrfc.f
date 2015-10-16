C     This function is used to setup donlp2.  There are 4 types of problems
c     can be setup with this donlp2usrfc file
c     probtype = 1  for ML estimation
c     2  for Lower Confidence
c     3  for Upper Confidence
c     4  for Log Likelihood Profile
c     5  for combined lower confidence limit
c     
C     The problem is to minimize f(D) = D, subject to the equality
c     constraints LIK(D, beta) - TARGET = 0, a constraint on D and beta
c     based on the benchmark response level and the type of risk used,
c     and possible equality constraints on some of the betas, and
c     possible inequality constraints of the form beta_i >= 0 for the
c     remaining betas.  The vector X contains all the parameters.  X(1)
c     contains D (the benchmark dose), X(2) ... X(polyord + 1) contain
c     the parameters for the dose-response model: X(2) contains b_0, etc.
C     Setup of X:
C     ML Estimation:
C     X(i) = beta[i-1]
C     Confidence Limit estimation
C     X(1) = log(D) (log(bmdl))
C     X(i) = beta[i-2]
C     LL Profile
c     X(i) = beta[i-1]
c     
      SUBROUTINE SETUP0
      INCLUDE 'O8COMM.INC'
      INCLUDE 'PROBLEM.INC'
      INCLUDE 'O8FINT.INC'
      INCLUDE 'MULTITUMOR.INC'
      INTEGER I,J,ioffset, K
C     NAME IS IDENT OF THE EXAMPLE/USER AND CAN BE SET AT USERS WILL
C     THE FIRST CHARACTER MUST BE ALPHABETIC.  40 CHARACTERS MAXIMUM
c     NAME = 'Multistage BMDL'
      IF (probtype .EQ. 1) THEN 
         NAME = 'ML'
      ELSEIF ((probtype .EQ. 2) .OR. (probtype .EQ. 3)) THEN
         NAME = 'BMDL'
      ELSEIF (probtype .EQ. 4) THEN
c     WRITE ONTO NAME
         IF (counter .LE. 9) THEN
            WRITE (NAME, '(A,I1)') 'LLPro0',counter
         ELSE
            WRITE (NAME,'(A,I2)') 'LLPro',counter
         ENDIF
C     *** added by CVL 7/2007 -start block
      ELSEIF (probtype .EQ. 5) THEN
         NAME = 'Combo BMDL'
C        Print *,'Name =',Name
C     *** added by CVL 7/2007 end block
      ENDIF

C     Scale by the starting values
      IF (probtype .EQ. 4) THEN
         DO I = 1, N
            IF (X(I) .GT. 0.) THEN
               XSC(I) = X(I)
            ELSE
               XSC(I) = 1.0
            ENDIF
         ENDDO
      ENDIF
C     X IS INITIAL GUESS AND ALSO HOLDS THE CURRENT SOLUTION
C     PROBLEM DIMENSION N=DIM(X), NH=DIM(H), NG=DIM(G)
C     
      IF (probtype .EQ. 1) THEN
C     ML Estimation
C     
C     Objective function is a function of all the parameters
         GUNIT(1,0) = -1
         GUNIT(2,0) = 0
         GUNIT(3,0) = 0
C     The only equality constraints are those imposed by presetting some
C     parameters, as determined by parmfixed(i)
         NH = 0
         J = 0
         DO I = 0, polyord
            IF (parmfixd(I) .EQ. 1) THEN
               J = J + 1
               NH = NH + 1
               GUNIT(1, J) = 1
               GUNIT(2, J) = I + 1
               GUNIT(3, J) = 1
            ENDIF
         ENDDO
c     Below, I is the argument to EG()
C     If beta[0] is not fixed, then it must be non-negative: I = 1
         NG = 0
         IF (parmfixd(0) .NE. 1) THEN
            J = J + 1
            NG = NG + 1
            GUNIT(1, J) = 1
            GUNIT(2, J) = 1
            GUNIT(3, J) = 1
         ENDIF
C     If restrict .EQ. 1, then all other non-fixed betas must be non
c     -negative: I = 2:(polyord+1)
         IF (restrict .EQ. 1) THEN
            DO I = 1, polyord
               IF (parmfixd(I) .NE. 1) THEN
                  J = J + 1
                  NG = NG + 1
                  GUNIT(1, J) = 1
                  GUNIT(2, J) = I + 1
                  GUNIT(3, J) = 1
               ENDIF
            ENDDO
         ENDIF
      ELSEIF (probtype .EQ. 2) THEN
C     Lower Confidence Limit
C     
C     Objective function is linear in X(1)
         GUNIT(1,0)=1
         GUNIT(2,0)=1
         GUNIT(3,0)=1

c     There is always at least one equality constraint,
c     imposed by the benchmark response level.  Additional equality
c     constraints are imposed by fixing specific betas to particular
c     values.
         NH=0
C     GUNIT-ARRAY, SEE DONLP2DOC.TXT
C     First equality constraint is computed
         J = 1
         NH = NH + 1
         GUNIT(1,J)=-1
         GUNIT(2,J)=0
         GUNIT(3,J)=0
C     Fixed parameters
         DO I = 0, polyord
            IF (parmfixd(I) .EQ. 1) THEN
               NH = NH + 1
               GUNIT(1,J)=1
               GUNIT(2,J)=I + 2
               GUNIT(3,J)=1
               J=J+1
            ENDIF
         ENDDO
c     There are always at least three inequality constraints:
c       ~ log-likelihood > target 
c       ~ X(1) - double_eps .GE. 0 
c       ~ X(1) .LT. bmd.  
c     In addition, if the background is
c     not fixed, then X(2) must be .GE. 0.  Finally, if restrict .EQ. 1,
c     then there are an additional polyord inequality constraints.     
         NG = 0
         J = NH
c     log-likelihood - target > 0: I = 1
         J = J + 1
         NG = NG + 1
         GUNIT(1, J) = -1
         GUNIT(2, J) = 0
         GUNIT(3, J) = 0
c     X(1) - (log(double_min) - log(scale)) >=0: I = 2
         J = J + 1
         NG = NG + 1
         GUNIT(1,J)=1
         GUNIT(2,J)=1
         GUNIT(3,J)=1
c     BMD - X(1) >= 0: I = 3
         J = J + 1
         NG = NG + 1
         GUNIT(1,J)=1
         GUNIT(2,J)=1
         GUNIT(3,J)=-1
c     X(2) >= 0: I = 3 + parmfixd
         IF (parmfixd(0) .NE. 1) THEN
            NG = NG + 1
            J=J+1
            GUNIT(1,J)=1
            GUNIT(2,J)=2
            GUNIT(3,J)=1
         ENDIF
c     X(J) >= 0: I > 3, if restrict .EQ. 1
         IF (restrict .EQ. 1) THEN
            DO I = 1, polyord
               IF (parmfixd(I) .NE. 1) THEN
                  NG = NG + 1
                  J=J+1
                  GUNIT(1,J)=1
                  GUNIT(2,J)=I+2
                  GUNIT(3,J)=1
               ENDIF
            ENDDO
         ENDIF
      ELSEIF (probtype .EQ. 3) THEN
C     
C     Upper confidence limit
C     
C     Objective function is linear in X(1)
         GUNIT(1,0)=1
         GUNIT(2,0)=1
         GUNIT(3,0)=-1
c     There is always at least one equality constraint,
c     imposed by the benchmark response level.  Additional equality
c     constraints are imposed by fixing specific betas to particular
c     values.
         NH=0
         J = 0
C     GUNIT-ARRAY, SEE DONLP2DOC.TXT
C     First equality constraint is computed
         NH = NH + 1
         J = J + 1
         GUNIT(1,J)=-1
         GUNIT(2,J)=0
         GUNIT(3,J)=0
C     Fixed Parameters
         DO I = 0, polyord
            IF (parmfixd(I) .EQ. 1) THEN
               NH = NH + 1
               J = J + 1
               GUNIT(1,J)=1
               GUNIT(2,J)=I + 2
               GUNIT(3,J)=1
            ENDIF
         ENDDO
c     There are always at least three inequality constraints:
C       ~ log-likelihood - target > 0
C       ~ X(1) .GE. log(bmd)
C       ~ log(double_max) - log(scale) - X(1) > 0
C     In addition, if the background is
c     not fixed, then X(2) must be .GE. 0.  Finally, if restrict .EQ. 1,
c     then there are an additional polyord inequality constraints.     
         NG = 0
c     log-likelihood - target > 0: I = 1
         NG = NG + 1
         J = J + 1
         GUNIT(1, J) = -1
         GUNIT(2, J) = 0
         GUNIT(3, J) = 0
c     X(1) - log(BMD) >= 0: I = 2
         NG = NG + 1
         J = J + 1
         GUNIT(1,J)=1
         GUNIT(2,J)=1
         GUNIT(3,J)=1
C     log(max_double) - log(scale) - X(1) > 0: I = 3
         NG = NG + 1
         J = J + 1
         GUNIT(1,J)=1
         GUNIT(2,J)=1
         GUNIT(3,J)=-1
c     X(2) >= 0: I = 3 + parmfixd(0)
         IF (parmfixd(0) .NE. 1) THEN
            NG = NG + 1
            J = J + 1
            GUNIT(1,J)=1
            GUNIT(2,J)=2
            GUNIT(3,J)=1
         ENDIF
c     X(J) >= 0: I > 3 + parmfixd(0)
         IF (restrict .EQ. 1) THEN
            DO I = 1, polyord
               IF (parmfixd(I) .NE. 1) THEN
                  NG = NG + 1
                  J=J+1
                  GUNIT(1,J)=1
                  GUNIT(2,J)=I+2
                  GUNIT(3,J)=1
               ENDIF
            ENDDO
         ENDIF
         
      ELSEIF (probtype .EQ. 4) THEN
c     
c     Find Profile
c     
C     Objective function is a function of all the parameters
         GUNIT(1,0) = -1
         GUNIT(2,0) = 0
         GUNIT(3,0) = 0

C     The equality constraints are those imposed by presetting some
C     parameters, as determined by parmfixed(i), and the one imposed
C     by the given BMD
         NH = 0
         J = 0
         DO I = 0, polyord
            IF (parmfixd(I) .EQ. 1) THEN
               J = J + 1
               NH = NH + 1
               GUNIT(1, J) = 1
               GUNIT(2, J) = I + 1
               GUNIT(3, J) = 1
            ENDIF
         ENDDO

C     The constraint imposed by the BMD is a function of all the
C     parameters
c     beta[1]*BMD + beta[2]*BMD^2 +...+ beta[n]*BMD^n - ln(1-A) = 0
c     
         NH = NH + 1
         J = J + 1
         GUNIT(1, J) = -1
         GUNIT(2, J) = 0
         GUNIT(3, J) = 0
C     If restrict .EQ. 1, then all non-fixed betas must be non
c     -negative
         NG = 0
C     If beta[0] is not fixed, then it must be non-negative
         IF (parmfixd(0) .NE. 1) THEN
            J = J + 1
            NG = NG + 1
            GUNIT(1, J) = 1
            GUNIT(2, J) = 1
            GUNIT(3, J) = 1
         ENDIF
         IF (restrict .EQ. 1) THEN
            DO I = 1, polyord
               IF (parmfixd(I) .NE. 1) THEN
                  J = J + 1
                  NG = NG + 1
                  GUNIT(1, J) = 1
                  GUNIT(2, J) = I + 1
                  GUNIT(3, J) = 1
               ENDIF
            ENDDO
         ENDIF
C     *** added by CVL 8/2007 -start block
      ELSEIF (probtype .EQ. 5) THEN
C
C     Combination Lower Confidence Limit 
C     
C     Objective function is linear in X(1)
         GUNIT(1,0)=1
         GUNIT(2,0)=1
         GUNIT(3,0)=1

c     There is always at least one equality constraint,
c     imposed by the benchmark response level.  Additional equality
c     constraints are imposed by fixing specific betas to particular
c     values.
         NH=0
C     GUNIT-ARRAY, SEE DONLP2DOC.TXT
C     First equality constraint is computed
         J = 1
         NH = NH + 1
         GUNIT(1,J)=-1
         GUNIT(2,J)=0
         GUNIT(3,J)=0
C     Fixed parameters but not allowing linear terms to be fixed
C  Changed by BCA 01/09 start *****************************
      ioffset = 0
      DO K = 0, nT-1
         if (K.gt.0) ioffset = ioffset + arnPoly(K-1)
         DO I = 0, arnPoly(K)
            IF (mnParmFixd(K,I) .EQ. 1 .AND. I .ne. 1) THEN
               NH = NH + 1
               J = J + 1
               GUNIT(1,J)=1
               GUNIT(2,J)=I + ioffset + K + 2
               GUNIT(3,J)=1
               parmval(J) = mdParmVal(K,I)
            ENDIF
         ENDDO
      ENDDO
C  Changed by BCA 01/09 end *****************************

c     There are always at least three inequality constraints:
c       ~ log-likelihood > target 
c       ~ X(1) - double_eps .GE. 0 
c       ~ X(1) .LT. bmd.  
c     In addition, if the background is
c     not fixed, then X(2) and X(polyorb + 3) must be .GE. 0.  
c     Finally, if restrict .EQ. 1,
c       ~ X(3)  .GE. X(polyorb + 4)
c     and there are an additional polyord and apolyord inequality constraints.     
         NG = 0
         J = NH
c     log-likelihood - target > 0: I = NH + 1
         J = J + 1
         NG = NG + 1
         GUNIT(1, J) = -1
         GUNIT(2, J) = 0
         GUNIT(3, J) = 0
c     X(1) - (log(double_min) - log(scale)) >=0: I = NH + 2
         J = J + 1
         NG = NG + 1
         GUNIT(1,J)=1
         GUNIT(2,J)=1
         GUNIT(3,J)=1
c     BMD - X(1) >= 0: I = 3 (on log scale)
         J = J + 1
         NG = NG + 1
         GUNIT(1,J)=1
         GUNIT(2,J)=1
         GUNIT(3,J)=-1
         IF (restrict .EQ. 1) THEN
c     X(3) > X(polyorb + 4): I=4
            J = J + 1
            NG = NG + 1
            GUNIT(1, J) = -1
            GUNIT(2, J) = 0
            GUNIT(3, J) = 0
         endif
C   Changed by BCA 01/09 start *********************
C   Background terms greater than or equal to zero
      ioffset = 0
      DO K = 0, nT-1
         if (K.gt.0) ioffset = ioffset + arnPoly(K-1)
         IF (mnParmFixd(K,0) .ne. 1 ) THEN
            NG = NG + 1
            J = J + 1
            GUNIT(1,J)=1
            GUNIT(2,J)=ioffset + K + 2
            GUNIT(3,J)=1
         ENDIF
      ENDDO
C   changed by BCA 01/09 end **************************
c     X(J) >= 0: I in (3,..,poly+2,polyord+4 .. N) if restrict .EQ. 1
c    This is a computed constraint for the reparmeratized model S variable
C   changed by BCA 01/09 start **************************
      IF (restrict .EQ. 1) THEN
         DO I = 2, arnPoly(0)
            IF (mnParmFixd(0,I) .NE. 1) THEN
               NG = NG + 1
               J=J+1
               GUNIT(1,J)=1
               GUNIT(2,J)=I+2
               GUNIT(3,J)=1
            ENDIF
         ENDDO
         ioffset = 0
         DO K = 1, nT-1
            ioffset = ioffset + arnPoly(K-1)
            DO I = 1, arnPoly(K)
               IF (mnParmFixd(K,I) .ne. 1 ) THEN
                  NG = NG + 1
                  J = J + 1
                  GUNIT(1,J)=1
                  GUNIT(2,J)=I+ioffset + K + 2
                  GUNIT(3,J)=1
               ENDIF
            ENDDO
         ENDDO
      ENDIF
C   changed by BCA 01/09 end **************************
       PRINT *,'I','GUNIT (1,I)','GUNIT(2,I)','GUNIT(3,I)'
       DO I=0,NH+NG
       PRINT *,I,GUNIT(1,I),GUNIT(2,I),GUNIT(3,I)
       ENDDO
C     *** added by CVL 8/2007 - end block
      ENDIF
      ANALYT=.TRUE.
      COLD=.TRUE.  
C **** CVL 8/2007 modified IF block below for silent running
C      of donlp2 for combo BMDL
      IF (probtype .EQ. 5) THEN
          SILENT = .TRUE.
          TE0=.TRUE.
          TE1=.TRUE.
          TE2=.TRUE.
          TE3=.TRUE.
      ELSE
          SILENT=.TRUE.
      ENDIF
      EPSDIF=0.00001
      PROU=10
      MEU=20
      
C     DEL0 AND TAU0: SEE BELOW
      IF ((probtype .EQ. 1) .OR. (probtype .EQ. 4)) THEN
         DEL0 = 1.0D0
         TAU0 = 1.0D0
      ELSE
         DEL0=1.0D-1
         TAU0=2.0D0
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
      INCLUDE 'MULTITUMOR.INC'
      DOUBLE PRECISION X(*),FX, SUM, PROBS(MAXORDER + 1), SLOGF
      INTEGER K, J
      ICF=ICF+1
      IF (probtype .EQ. 1) THEN
C     ML Estimation
C     PRINT *,N,'ML ESTIMATION N'
         DO K=1, ndoses
            SUM=X(N)
            DO J=N-1,1,-1
               SUM = SUM * doses(K) + X(J)
            ENDDO
            PROBS(K) = 1 - DEXP(-SUM)
         ENDDO
         SUM = 0.D0
         DO K=1,ndoses
            SUM = SUM + DBLE(affected(K)) * SLOGF(PROBS(K)) +
     $           DBLE(nanimals(K)-affected(K)) * SLOGF(1.0D0 - PROBS(K
     $           ))
         ENDDO
         FX = -SUM
      ELSEIF (probtype .EQ. 2) THEN
C     Lower Confidence Limit
C     Find the largest D less than MLE s.t. likelihood and BMR constraints hold
         FX=X(1)
      ELSEIF (probtype .EQ. 3) THEN
C     Upper Confidence Limit
         FX = -X(1)
      ELSEIF (probtype .EQ. 4) THEN
C     BMD Profile
c     PRINT *,polyord,ndoses
         DO K=1, ndoses
            SUM=X(N)
            DO J=N-1,1,-1
               SUM = SUM * doses(K) + X(J)
            ENDDO
            PROBS(K) = 1 - DEXP(-SUM)
         ENDDO
         SUM = 0.D0
         DO K=1,ndoses
            SUM = SUM + DBLE(affected(K)) * SLOGF(PROBS(K)) +
     $           DBLE(nanimals(K)-affected(K)) * SLOGF(1.0D0 - PROBS(K
     $           ))
         ENDDO
         FX = -SUM         
C     *** Added by CVL 8/2007  ** start block
      ELSEIF (probtype .EQ. 5) THEN
C     Combined Lower Confidence Limit
C     Find the largest smallest D less than MLE s.t. likelihood and BMR constraints hold
         FX=X(1)
C     *** Added by CVL 8/2007  ** end block
      ENDIF
      RETURN
      END

C     GRADIENT OF OBJECTIVE FUNCTION
      SUBROUTINE EGRADF(X,GRADF)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      INCLUDE 'MULTITUMOR.INC'
      DOUBLE PRECISION X(*),GRADF(*), DSLOG, P, RESID
C     USER DECLARATIONS, IF ANY ,FOLLOW
      INTEGER J, K
      ICGF=ICGF+1
      DO J = 1, N
         GRADF(J) = 0.0D0
      ENDDO
      IF (probtype .EQ. 1) THEN
C     ML Estimation
         DO K = 1, ndoses
c     Compute the probability of a response
            P = X(N)
            DO J = N-1,1,-1
               P = P * doses(K) + X(J)
            ENDDO
            P = 1.0D0 - DEXP(-P)
c     Compute the 'residual': affected(K)/P - nanimals(K), but do it so
c     that it won't blow up at P == 0.
            RESID = (DBLE(affected(K)) * DSLOG(P) - DBLE(nanimals(K)
     $           -affected(K)) * DSLOG(1.0D0 - P)) * (1.0D0 - P) 
c     Components of the gradient
            GRADF(1) = GRADF(1) - RESID
            DO J = 2, N
               GRADF(J) = GRADF(J) - RESID * (doses(K)**(J-1))
            ENDDO   
         ENDDO
C     
      ELSEIF (probtype .EQ. 2) THEN
C     Lower Confidence Limit
         GRADF(1) = 1.D0
      ELSEIF (probtype .EQ. 3) THEN
C     Upper Confidence Limit
         GRADF(1) = -1.D0
c     
C     *** Modified by CVL 8/2007 to make Else if
C     ELSE
      ELSEIF (probtype .EQ. 4) THEN
C     LL Profile
         DO K = 1, ndoses
c     Compute the probability of a response
            P = X(N)
            DO J = N-1,1,-1
               P = P * doses(K) + X(J)
            ENDDO
            P = 1.0D0 - DEXP(-P)
c     Compute the 'residual': affected(K)/P - nanimals(K), but do it so
c     that it won't blow up at P == 0.
            RESID = (DBLE(affected(K)) * DSLOG(P) - DBLE(nanimals(K)
     $           -affected(K)) * DSLOG(1.0D0 - P)) * (1.0D0 - P) 
c     Components of the gradient
            GRADF(1) = GRADF(1) - RESID
            DO J = 2, N
               GRADF(J) = GRADF(J) - RESID * (doses(K)**(J-1))
            ENDDO   
         ENDDO
C     *** Added by CVL 8/2007  ** start block
      ELSEIF (probtype .EQ. 5) THEN
C     Lower Confidence Limit
         GRADF(1) = 1.D0
C     *** Added by CVL 8/2007  ** end block
      ENDIF
      RETURN
      END
      
c     Derivative of slog wrt p
      DOUBLE PRECISION FUNCTION DSLOG(P)
      DOUBLE PRECISION P
      IF (P .GE. 1.D-10) THEN
         DSLOG = 1.0D0/P
      ELSE
         DSLOG = 2.0D10 - 1.0D20 * P
      ENDIF
      RETURN
      END

C     COMPUTE THE I-TH EQUALITY CONSTAINT, VALUE IS HXI
      SUBROUTINE EH(I,X,HXI)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      INCLUDE 'MULTITUMOR.INC'
      INTEGER I
      DOUBLE PRECISION X(*),HXI
      DOUBLE PRECISION PROBS(MAXORDER + 1), SUM, SLOGF 
      DOUBLE PRECISION S, S2, SUM3, BD, SUM2, D
      INTEGER J, K, ioffset, L
      CRES(I)=CRES(I)+1
      IF (probtype .EQ. 1) THEN
C     ML Estimation
         HXI = X(GUNIT(2,I)) - parmval(GUNIT(2,I) - 1)
      ELSEIF (probtype .EQ. 2) THEN
C     Confidence Limits
         IF (I .EQ. 1) THEN
c     BMR-derived constraint
            D = DEXP(X(1))
            IF (risktype .EQ. 1) THEN
               SUM =  LOG(1.0 - bmr)
            ELSE
               SUM =  SLOGF(1.0 - bmr*DEXP(X(2)))
            ENDIF
            SUM2 = X(N)*D
            DO K=N-1, 3, -1
               SUM2 = SUM2*D+X(K)*D
            ENDDO
            HXI = SUM + SUM2
         ELSE
c     Parameter fixed to a given value
            HXI=X(GUNIT(2,I)) - parmval(GUNIT(2,I)-2)
         ENDIF
      ELSEIF (probtype .EQ. 3) THEN
C     Confidence Limits
         IF (I .EQ. 1) THEN
c     BMR-derived constraint
            D = DEXP(X(1))
            IF (risktype .EQ. 1) THEN
               SUM =  LOG(1.0 - bmr)
            ELSE
               SUM =  SLOGF(1.0 - bmr*DEXP(X(2)))
            ENDIF
            SUM2 = X(N)*D
            DO K=N-1, 3, -1
               SUM2 = SUM2*D+X(K)*D
            ENDDO
            HXI = SUM + SUM2
         ELSE
c     Parameter fixed to a given value
            HXI=X(GUNIT(2,I)) - parmval(GUNIT(2,I)-2)
         ENDIF
      ELSEIF (probtype .EQ. 4) THEN
C     LL Profile
         IF (I .NE. NH) THEN
c     Parameter fixed to a given value
            HXI = X(GUNIT(2,I)) - parmval(GUNIT(2,I) - 1)
         ELSE          
c     BMD/BMR-derived constraint
            IF (risktype .EQ. 1) THEN
               SUM = LOG(1.0 - bmr)
            ELSE
               IF (X(1) .GT. 700) X(1) = 700.0D0
               SUM = SLOGF(1.0 - bmr*(DEXP(X(1))))
            ENDIF
            DO K=2, N
               SUM = SUM + X(K) * bmd**(K-1)
            ENDDO
            HXI = SUM         
         ENDIF
C     *** Added by CVL 8/2007  ** start block
C    Note that the S is the reparamatized linear term 
C
      ELSEIF (probtype .EQ. 5) THEN
C     Confidence Limits
      IF (I .EQ. 1) THEN
c     BMR-derived constraint
C     Changed by BCA 01/09 start *************************
      D = DEXP(X(1))
      HXI =  LOG(1.0 - bmr)
      ioffset = 0
      DO L=1,nT-1
         ioffset = ioffset + arnPoly(L-1)
      ENDDO
      DO L=nT-1,0,-1
         IF(L.lt.nT-1) ioffset = ioffset-arnPoly(L)
         SUM = X(ioffset+L+arnPoly(L)+2)*D
         DO K=arnPoly(L)-1,2,-1
            SUM = SUM*D+X(ioffset+L+K+2)*D
         ENDDO
         If (arnPoly(L) .eq. 1)  SUM = 0
         HXI = HXI + SUM*D
      ENDDO
      HXI = HXI + X(3)*D
      ELSE
c     Parameter fixed to a given value
            HXI=X(GUNIT(2,I)) - parmval(I)
      ENDIF
        print *, 'I=',i,' HXI =', HXI,' G(2,I)=', Gunit(2,I)
C     *** Added by CVL 8/2007  ** end block
c	end of probtype = 5 in EH subroutine
C     Changed by BCA 01/09 end ***********************************

c        print *, 'I=',i,' HXI =', HXI,' G(2,I)=', Gunit(2,I)
C     *** Added by CVL 8/2007  ** end block
      ENDIF
      RETURN
      END
C     
C     COMPUTE THE GRADIENT OF THE I-TH EQUALITY CONSTRAINT
      SUBROUTINE EGRADH(I,X,GRADHI)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      INCLUDE 'MULTITUMOR.INC'
      DOUBLE PRECISION X(*),GRADHI(*), SUM, RESID, P, DSLOG, D
      DOUBLE PRECISION S, S2, SUM2, SUM3, BD, tempsum
      INTEGER I,J,K,ioffset,L
      IF ( GUNIT(1,I) .NE. 1 ) CGRES(I)=CGRES(I)+1
      DO  J=1,NX
         GRADHI(J)=0.D0
      ENDDO
C     
      IF (probtype .EQ. 1) THEN
C     ML Estimation
C     Do nothing, since these are taken care of in GUNIT
      ELSEIF(probtype .EQ. 2) THEN
C     Lower Confidence Limit
         D = DEXP(X(1))
         IF (I .EQ. 1) THEN
c     BMD/BMR constraint
c     J = 1                 
            SUM = DBLE(N-2)*X(N)
            DO K = N-1,3,-1
               SUM = SUM * D + DBLE(K-2)*X(K)
            ENDDO 
            GRADHI(1) = SUM * D
c     J = 2
            IF (risktype .EQ. 1) THEN
               GRADHI(2) = 0.D0
            ELSE
               GRADHI(2) = - DBLE(bmr) * DEXP(X(2)) / (1.D0 - DBLE(bmr)
     $              * DEXP(X(2)))
            ENDIF
c     J >= 3
            DO J = 3, N
               GRADHI(J) = D**(J - 2)
            ENDDO
         ELSE
c     Individual parameter constraint
            GRADHI(GUNIT(2,I))=1.0D0
         ENDIF
c     
      ELSEIF(probtype .EQ. 3) THEN
         D = EXP(X(1))
C     Upper Confidence Limit
         IF (I .EQ. 1) THEN
c     BMD/BMR constraint
c     J = 1
            SUM = DBLE(N-2)*X(N)
            DO K = N-1,3,-1
               SUM = SUM * D + DBLE(K-2)*X(K)
            ENDDO 
            GRADHI(1) = SUM * D
c     J = 2
            IF (risktype .EQ. 1) THEN
               GRADHI(2) = 0.D0
            ELSE
               GRADHI(2) = - DBLE(bmr) * DEXP(X(2)) / (1.D0 - DBLE(bmr)
     $              * DEXP(X(2)))
            ENDIF
c     J >= 3
            DO J = 3, N
               GRADHI(J) = D**(J - 2)
            ENDDO
         ELSE
c     Individual parameter constraint
            GRADHI(GUNIT(2,I))=1.0D0
         ENDIF
c     
c     Profile Likelihood
      ELSEIF (probtype .EQ. 4) THEN
         IF (I .EQ. NH) THEN
            IF (risktype .EQ. 1) THEN
               GRADHI(1) = 0.D0
            ELSE
               GRADHI(1) = - DBLE(bmr) * DEXP(X(1)) / (1.D0 - DBLE(bmr)
     $              * DEXP(X(1)))
            ENDIF
            DO J = 2, N
               GRADHI(J) = bmd**(J-1)
            ENDDO
         ENDIF
C     *** Added by CVL 8/2007  ** start block
C           - the X(polyord + 4) = Model A q1 
C           - the X(3) = Model B q1 or rather should be S
C
C     Changed by BCA 01/09 start *************************
      ELSEIF(probtype .EQ. 5) THEN
C     Lower Confidence Limit
         D = DEXP(X(1))
        IF (I .EQ. 1) THEN
c     BMD/BMR constraint
c     J = 1
      ioffset = 0
      DO L=1,nT-1
         ioffset = ioffset + arnPoly(L-1)
      ENDDO
      tempsum = 0.0
      DO L=nT-1,0,-1
         IF(L.lt.nT-1) ioffset = ioffset-arnPoly(L)
         SUM = X(ioffset+L+arnPoly(L)+2)*DBLE(arnPoly(L))
         DO K=arnPoly(L)-1,2,-1
            SUM = SUM*D+X(ioffset+L+K+2)*DBLE(K)
         ENDDO
         If (arnPoly(L) .eq. 1)  SUM = 0
         tempsum = tempsum + SUM*D
      ENDDO
      tempsum = tempsum + X(3)
      GRADHI(1) = tempsum*D
       
c     J = 2 and higher
C            Only for risktype .EQ. 1 (Extra Risk) Background gradients 
C            and the gradient for linear terms on model B on up 
C            are zero
      ioffset = 0
      GRADHI(2) = 0.D0
      GRADHI(3) = D
      DO J = 4, arnPoly(0) + 2
         GRADHI(J) =  D**(J - 2)
      ENDDO
      DO K=1, nT-1
         ioffset = ioffset+arnPoly(K-1)
         GRADHI(ioffset+K+2) = 0.0
         GRADHI(ioffset+K+3) = 0.0
         DO J = 2, arnPoly(K)
            GRADHI(ioffset+K+J+2) = D**J
         ENDDO
      ENDDO
              
      ELSE
c     Individual parameter constraint
            GRADHI(GUNIT(2,I))=1.0D0
C    Changed by BCA 01/09 end ******************************
      ENDIF
C     *** Added by CVL 8/2007  ** end block
      ENDIF
      RETURN
      END

C     COMPUTE THE I-TH INEQUALITY CONSTAINT, BOUNDS INCLUDED

      SUBROUTINE EG(I,X,GXI)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      INCLUDE 'MULTITUMOR.INC'
      DOUBLE PRECISION PROBS(MAXORDER + 1),X(*),GXI,SUM, S, S2
      DOUBLE PRECISION SLOGF, SUM2, SUM1(MAXORDER + 1), BD,
     $     correction, mPROBS(1:100,0:10)
      INTEGER I, J, K, ioffset, L
      IF ( GUNIT(1,I+NH) .EQ. -1 ) CRES(I+NH)=CRES(I+NH)+1
      
      IF (probtype .EQ. 1) THEN
C     ML Estimation
         GXI = X(GUNIT(2, I+NH))
      ELSEIF (probtype .EQ. 2) THEN
C     Lower Confidence Limit
         IF (I .EQ. 1) THEN
c     Likelihood - target
            DO K=1, ndoses
               SUM=X(N)
               DO J=N-1,2,-1
                  SUM = SUM * doses(K) + X(J)
               ENDDO
               IF (SUM .LT. 0) THEN
                  SUM = 0.0D0
               ELSE
               ENDIF
               PROBS(K) = 1 - DEXP(-SUM)
            ENDDO
            SUM = -target
            DO K=1,ndoses
               SUM = SUM + DBLE(affected(K)) * SLOGF(PROBS(K)) +
     $              DBLE(nanimals(K)-affected(K)) * SLOGF(1.0D0 - 
     $              PROBS(K))
            ENDDO
            GXI=SUM
         ELSEIF (I .EQ. 2) THEN
            GXI = X(1) - lminbmd
         ELSEIF (I .EQ. 3) THEN
            GXI = bmd - X(1)
         ELSE
            GXI=X(GUNIT(2,I+NH))
         ENDIF
      ELSEIF (probtype .EQ. 3) THEN
C     Upper Confidence Limit
         IF (I .EQ. 1) THEN
c     Likelihood - target
            DO K=1, ndoses
               SUM=X(N)
               DO J=N-1,2,-1
                  SUM = SUM * doses(K) + X(J)
               ENDDO
               IF (SUM .LT. 0) THEN
                  SUM = 0.0D0
               ELSE
               ENDIF
               PROBS(K) = 1 - DEXP(-SUM)
            ENDDO
            SUM = -target
            DO K=1,ndoses
               SUM = SUM + DBLE(affected(K)) * SLOGF(PROBS(K)) +
     $              DBLE(nanimals(K)-affected(K)) * SLOGF(1.0D0 - 
     $              PROBS(K))
            ENDDO
            GXI=SUM
         ELSEIF (I .EQ. 2) THEN
            GXI = X(1) - bmd
         ELSEIF (I .EQ. 3) THEN
            GXI = lmaxbmd - X(1)
         ELSE
            GXI=X(GUNIT(2,I+NH))
         ENDIF
      ELSEIF ((probtype .EQ. 4) .AND. (restrict .EQ. 1)) THEN
C     LL Profile
         GXI = X(GUNIT(2, I+NH))
C     *** Added by CVL 8/2007  ** start block
      ELSEIF (probtype .EQ. 5) THEN
C
C     Combined Lower Confidence Limit
C
        IF (I .EQ. 1) THEN
C
C     Likelihood - target
C
C    Changed by BCA 01/09 start*******************************
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
   			IF(ioffset+L+arnPoly(L)+2.eq.3) THEN
   			  SUM = X(ioffset+L+arnPoly(L)+2) - correction
   			ENDIF
   			DO J = arnPoly(L)-1, 0,-1
   			   IF(L.eq.0 .AND. J.eq.1) THEN
   			     SUM = SUM * mdDoses(L,K)+X(ioffset+L+J+2) - correction
               ELSE
                 SUM = SUM * mdDoses(L,K) + X(ioffset+L+J+2)
               ENDIF
            ENDDO
   			IF (SUM .LT. 0) SUM = 0.0D0
   			mPROBS(K,L) = 1 - DEXP(-SUM)
		ENDDO
        
		DO K=1,arnObs(L)
		   SUM2 = SUM2 + DBLE(mnAffected(L,K)) * SLOGF(mPROBS(K,L)) +
     $              DBLE(mnAnim(L,K)-mnAffected(L,K)) * SLOGF(1.0D0 - 
     $              mPROBS(K,L))
		ENDDO
      ENDDO
c	end of do over L
        GXI=-target+SUM2
        ELSEIF (I .EQ. 2) THEN
            GXI = X(1) - lminbmd
        ELSEIF (I .EQ. 3) THEN
            GXI = bmd - X(1)
        ELSEIF (I .EQ. 4 .and. restrict .eq. 1) THEN
		ioffset = 0
		correction = 0.0
c			added the preceding 2 lines to reset ioffset and correction -- bca 3/18/09
		DO L=1, nT-1
			ioffset = ioffset + arnPoly(L-1)
			correction = correction + X(ioffset+L+3)
		ENDDO            
		GXI = X(3) - correction
        ELSE
            GXI=X(GUNIT(2,I+NH))
        ENDIF
      print *, 'In EG: I=',i,' GXI =', GXI,' G(2,I+NH)=', Gunit(2,I+NH)
C     *** Added by CVL 8/2007  ** end block
C    Changed by BCA 01/09 end ****************************************
      ENDIF
      RETURN
      END

C     COMPUTE THE GRADIENT OF THE I-TH INEQUALITY CONSTRAINT
C     NOT NECESSARY FOR BOUNDS, BUT CONSTANT GRADIENTS MUST BE SET
C     HERE E.G. USING DCOPY FROM A DATA-FIELD
C     This code looks wrong, but all the inequality constraints are
C     linear functions of single parameters, and are thus taken care of
C     by GUNIT, so this should not even be called.
      SUBROUTINE EGRADG(I,X,GRADGI)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      INCLUDE 'MULTITUMOR.INC'
      DOUBLE PRECISION X(*) ,GRADGI(*), P, RESID, DSLOG
      DOUBLE PRECISION S, S2, SUM, SUM2, BD, correction
      INTEGER I, J, K, ioffset, L
      IF ( GUNIT(1,I+NH) .NE. 1 ) CGRES(I+NH)=CGRES(I+NH)+1
      DO  J=1,NX
         GRADGI(J)=0.D0
      ENDDO
      IF (probtype .EQ. 1) THEN
C     ML Estimation
C     Do nothing; taken care of by GUNIT
      ELSEIF (probtype .EQ. 2) THEN
C     Lower Confidence Limit
         IF (I .EQ. 1) THEN
c     Likelihood constraint: D is not a parameter
            GRADGI(1) = 0.D0
            DO K = 1, ndoses
c     Compute the probability of a response
               P = X(N)
               DO J = N-1,2,-1
                  P = P * doses(K) + X(J)
               ENDDO
               P = 1.0D0 - DEXP(-P)
c     Compute the 'residual': affected(K)/P - nanimals(K), but do it so
c     that it won't blow up at P == 0.
               RESID = (DBLE(affected(K)) * DSLOG(P) - DBLE(nanimals(K)
     $              -affected(K)) * DSLOG(1.0D0 - P)) * (1.0D0 - P) 
c     Components of the gradient
               GRADGI(2) = GRADGI(2) + RESID
               DO J = 3, N
                  GRADGI(J) = GRADGI(J) + RESID * (doses(K)**(J-2))
               ENDDO   
            ENDDO
         ENDIF
      ELSEIF (probtype .EQ. 3) THEN
C     Upper Confidence Limit
         IF (I .EQ. 1) THEN
c     Likelihood constraint: D is not a parameter
            GRADGI(1) = 0.D0
            DO K = 1, ndoses
c     Compute the probability of a response
               P = X(N)
               DO J = N-1,2,-1
                  P = P * doses(K) + X(J)
               ENDDO
               P = 1.0D0 - EXP(-P)
c     Compute the 'residual': affected(K)/P - nanimals(K), but do it so
c     that it won't blow up at P == 0.
               RESID = (DBLE(affected(K)) * DSLOG(P) - DBLE(nanimals(K)
     $              -affected(K)) * DSLOG(1.0D0 - P)) * (1.0D0 - P) 
c     Components of the gradient
               GRADGI(2) = GRADGI(2) + RESID
               DO J = 3, N
                  GRADGI(J) = GRADGI(J) + RESID * (doses(K)**(J-2))
               ENDDO   
            ENDDO
         ENDIF
C     *** Added by CVL 8/2007  ** start block
      ELSEIF (probtype .EQ. 5) THEN
C     Combined Lower Confidence Limit
C
         IF (I .EQ. 1) THEN
c     Likelihood constraint: D is not a parameter
            GRADGI(1) = 0.D0
C  Changed by BCA  01/09 start ********************
	    ioffset = 0
	    correction = 0.0
c			added the preceding 2 lines to reset ioffset and correction -- bca 3/18/09	    
            DO L=1, nT-1
		ioffset = ioffset + arnPoly(L-1)
		correction = correction + X(ioffset+L+3)
	    ENDDO            
            DO L=nT-1, 0, -1
	    	IF(L.lt.nT-1) ioffset = ioffset - arnPoly(L)
            	DO K=1, arnObs(L)
               		P=X(ioffset+L+arnPoly(L)+2)
	   		IF(ioffset+L+arnPoly(L)+2.EQ.3)
     $                       P= X(ioffset+L+arnPoly(L)+2)-correction
               DO J = arnPoly(L)-1, 0,-1
                  P = P * mdDoses(L,K) + X(ioffset+L+J+2)
                  IF(ioffset+L+arnPoly(L)+2.eq.3) THEN
                     P= P* mdDoses(L,K) +X(ioffset+L+arnPoly(L)+2)-
     $                  correction
                  ELSE
                     P = P * mdDoses(L,K) + X(ioffset+L+J+2)
                  ENDIF
                ENDDO
                P = 1.0D0 - DEXP(-P)
                RESID = (DBLE(mnAffected(L,K)) * DSLOG(P)-
     $                DBLE(mnAnim(L,K)-mnAffected(L,K)) * 
     $                DSLOG(1.0D0 - P)) * (1.0D0 - P) 
c     Components of the gradient
                	DO J = arnPoly(L), 1, -1
                		GRADGI(ioffset+L+J+2)=GRADGI(ioffset+L+J+2)+
     $                    RESID*(mdDoses(L,K)**J)
                	ENDDO   
                	GRADGI(ioffset+L+2) = GRADGI(ioffset+L+2) + RESID
            	ENDDO
            ENDDO
         ELSE 
            IF (restrict .EQ. 1 .and. I .eq. 4 ) THEN
		ioffset = 0
		DO L=1,nT-1
	    		ioffset = ioffset + arnPoly(L-1)
                	GRADGI(ioffset+L+3) = -1
		ENDDO
                GRADGI(3) = 1
            ENDIF
         ENDIF
c         print *, 'Inequality Constraint ',I
c         DO j = 1, N
c           print *, ' X(',j,') =', X(J),' Gradgi =', GRADGI(J)
c         ENDDO
C    Changed by BCA 01/09 end **********************************

C     *** Added by CVL 8/2007  ** end block
      ELSEIF ((probtype .EQ. 4) .AND. (restrict .EQ. 1)) THEN
C     ML Estimation
C     Do nothing; taken care of by GUNIT
c     ENDIF
      ENDIF                       
      RETURN
      END
c     
      SUBROUTINE SETUP
      INCLUDE 'O8COMM.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER I
C     CHANGE TERMINATION CRITERION
C     EPSX=..
C     CHANGE I/O-CONTROL
      INTAKT=.FALSE.
      TE0=.FALSE.
      TE2=.FALSE.
      TE3=.FALSE.
C***  NOW YOU GET FOR EVERY STEP A ONE-LINE-OUTPUT ON STDOUT
C     ......
      RETURN
      END
