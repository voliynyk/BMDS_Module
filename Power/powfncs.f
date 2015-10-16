c               This file contains three functions that are specific to the
c               to the power model when doing maximum likelihood and confidence
c               limit calculations using donlp2 and the associated user 
c               functions which are non-specific to the continuous model
c               that is being run.
c               POWGUNIT fills the GUNIT array in donlp2 with the model 
c                       specific equality, non-equality constraint as well as
c                       the objective function.
c               POWMEAN just calculates the mean specific to that model as a
c                       vector of means, one at each dose level
c               POWPART calculates the vector of partial derivates of the
c                       mean with respect to all model parameters.
c               

        SUBROUTINE POWFILLGUNIT
        INCLUDE 'O8COMM.INC'
        INCLUDE 'PROBLEM.INC'
        INTEGER I, J
C       
C       X IS INITIAL GUESS AND ALSO HOLDS THE CURRENT SOLUTION
C       PROBLEM DIMENSION N=DIM(X), NH=DIM(H), NG=DIM(G)
C       
C   Parameters when probtype .eq. 1:
C   X(1)  alpha (constvar .eq.1) or lalpha (constvar .eq. 0)
C   X(2)  rho
C   X(3)  control
C   X(4)  slope
C   X(5)  n (power)
C
C   Parameters when probtype .eq, 2:
C   X(1)  BMDL
C   X(2)  alpha (constvar .eq.1) or lalpha (constvar .eq. 0)
C   X(3)  rho
C   X(4)  control
C   X(5)  slope
C   X(6)  n (power parameter)

C   
        IF (probtype .EQ. 1) THEN
C       ML Estimation
C       
C       Objective function is a function of all the parameters
           GUNIT(1,0) = -1
           GUNIT(2,0) = 0
           GUNIT(3,0) = 0
C       The only equality constraints are those imposed by presetting some
C       parameters, as determined by parmfixed(i)
           NH = 0
           J = 0
           DO I = 0, 4
              IF (parmfixd(I) .EQ. 1) THEN
                 J = J + 1
                 NH = NH + 1
                 GUNIT(1, J) = 1
                 GUNIT(2, J) = I + 1
                 GUNIT(3, J) = 1
              ENDIF
           ENDDO
C       
           NG = 0
C       
C Inequality constraints: Numbers depend on constvar
C Constraint #            Constvar .EQ. 1          Constvar .EQ. 0
C     1                      alpha > 0               n > 0 
C                                                   (or >= 1, 
C                                                    depending on restrict)
C     2                    n > 0 (or >=1, depending) n <= 18
C     3                    n <= 18                   slope <> 0 (adverse)      
C     4                    slope <> 0                rho <= 18
C     5                    rho <= 18                
C       Variance coefficient parameter, alpha, must be positive
C       I: 1 (If constvar .eq .1)
           IF(constvar.EQ.1) THEN
              NG = NG + 1
              J = J + 1
              GUNIT(1, J) = 1
              GUNIT(2, J) = 1
              GUNIT(3, J) = 1 
           ENDIF
C       If restrict .EQ. 1, then the power parameter, n, must be >= 1
c       Otherwise, n >= 0
c       I: 1 + constvar
           NG = NG + 1
           J = J + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 5
           GUNIT(3, J) = 1         
C       
c       n parameter <= 18
c       I: 2 + constvar
           J = J + 1
           NG = NG + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 5
           GUNIT(3, J) = -1 
c       adverse direction
c       I: 3 + constvar
           J = J + 1
           NG = NG + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 4
           IF (adverse .EQ. 1) THEN
              GUNIT(3, J) = 1 
           ELSE
              GUNIT(3, J) = -1 
           ENDIF

c       rho variance parameter <= 18
c       I: 4 + constvar
           J = J + 1
           NG = NG + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 2
           GUNIT(3, J) = -1 

C       
C       
        ELSE
C       Lower Confidence Limit
C       
c       There are always at least one equality constraint: imposed by 
c       the benchmark response level.  Additional equality
c       constraints are imposed by fixing specific betas to particular
c       values.
           NH=1
C       GUNIT-ARRAY, SEE DONLP2DOC.TXT
C       Objective function is linear in X(1)
           GUNIT(1,0)=1
           GUNIT(2,0)=1
           GUNIT(3,0)=1
C       First equality constraint is computed
           J=1
           GUNIT(1,J)=-1
           GUNIT(2,J)=0
           GUNIT(3,J)=0
           DO I = 0, 4
              IF (parmfixd(I) .EQ. 1) THEN
                 J = J+1
                 NH = NH + 1
                 GUNIT(1,J)=1
                 GUNIT(2,J)=I + 2
                 GUNIT(3,J)=1
              ENDIF
           ENDDO
c       There are always at least three inequality constraints, since 
c       1)  the loglikelihood must be >= target loglikelihood
c       2) X(1) must be .GE. 0 and 
c       3) X(1) must be .LT. bmd.  

c       In addition, if the variance coefficent
c       is not fixed, then X(2) must be .GE. 0.  Finally, if restrict .EQ. 1,
c       then there are an additional inequality constraints.     
           NG=0
C Inequality constraints
C Constraint #            Constvar .EQ. 1          Constvar .EQ. 0
C     1                     loglik >= target         loglik >= target
C     2                     BMDL > 0                 BMDL > 0 
C     3                     BMDL < BMD               BMDL < BMD
C     4                      alpha > 0               n > 0 
C                                                   (or >= 1, 
C                                                    depending on restrict)
C     5                    n > 0 (or >=1, depending) n <= 18
C     6                    n <= 18                   slope <> 0 (adverse)      
C     7                    slope <> 0                rho <= 18
C     8                    rho <= 18                
c       
C In the comments below, I refers to the argument passed to EG()
c       likelihood >= target I = 1
           NG = NG + 1
           J = J + 1
           GUNIT(1, J) = -1
           GUNIT(2, J) = 0
           GUNIT(3, J) = 0
c       
c       X(1) >= 0  I = 2
           NG = NG + 1
           J = J+1
           GUNIT(1,J) = 1
           GUNIT(2,J) = 1
           GUNIT(3,J) = 1
c       BMD - X(1) >= 0 I = 3
           NG = NG + 1
           J = J+1
           GUNIT(1,J)=1
           GUNIT(2,J)=1
           GUNIT(3,J)=-1

C       Variance coefficient parameter, alpha, must be positive
C       I = 4 (if constvar .eq. 1)
           IF(constvar .EQ. 1) THEN
              J = J + 1
              NG = NG + 1
              GUNIT(1, J) = 1
              GUNIT(2, J) = 2
              GUNIT(3, J) = 1 
           ENDIF
C       If restrict .EQ. 1, then the power parameter, n, must be >= 1
c       otherwise, n >= 0 I = 4 + constvar
           J = J + 1
           NG = NG + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 6
           GUNIT(3, J) = 1 
c       
           
c       n parameter <= 18
c       I = 5 + constvar
           J = J + 1
           NG = NG + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 6
           GUNIT(3, J) = -1 
c       adverse direction
c       I: 6 + constvar
           J = J + 1
           NG = NG + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 5
           IF (adverse .EQ. 1) THEN
              GUNIT(3, J) = 1 
           ELSE
              GUNIT(3, J) = -1 
           ENDIF
c       rho variance parameter <= 18
c       I: 7 + constvar
           J = J + 1
           NG = NG + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 3
           GUNIT(3, J) = -1        
C       
        ENDIF
        RETURN
        END
c       
        SUBROUTINE POWMEAN(X)
C       INCLUDE 'O8COMM.INC'
        INCLUDE 'PROBLEM.INC'
        INTEGER I
        DOUBLE PRECISION X(*), varmean, delta
c       
        IF(probtype.EQ.1) THEN
           DO I = 1, ndoses
              IF (doses(I).NE.0) THEN
                 means(I) = X(3) + X(4)*((doses(I))**X(5))
              ELSE
                 means(I) = X(3)
              ENDIF
           ENDDO
           varmean = 0.0
           DO I = 1, ndoses
              varmean = varmean + var(I)
           ENDDO
           varmean = varmean/ndoses
           
        ELSE
           
           IF (X(6) .LT. 0) X(6) = 0.0
           IF (X(6) .GT. 18) X(6) = 18.0
           IF (X(1) .GT. bmd) X(1) = bmd
           IF (risktype .EQ. 2) THEN
              IF (X(4) .LT. 1.0D-20) THEN
                 X(4) = 1.0D-20
              ELSE
              ENDIF
           ENDIF
           IF ((risktype .EQ. 4) .AND. (ABS(X(5)) .LT. 1.0D-20)) THEN
              IF (X(5) .GE. 0) THEN
                 X(5) = 1.0D-20
              ELSE
                 X(5) = -1.0D-20
              ENDIF
           ENDIF
           
           DO I = 1, ndoses
              IF (doses(I).NE.0) THEN
                 means(I) = X(4) + X(5)*((doses(I))**X(6))
              ELSE
                 means(I) = X(4)
              ENDIF
           ENDDO

           DoseMean = X(4)      

           IF(X(1).GT.0) THEN
              DoseMean = DoseMean + X(5)*(X(1)**X(6))
           ELSE
           ENDIF  
           
        ENDIF
c       
        RETURN
        END
c       
c       
c       POWpart gives the partial derivatives of the Power model mean function
c       and also, in grads(ndoses+1, j) gives the partial derivatives
c       of the mean function evaluated at the estimated dose D.
        SUBROUTINE POWPART(X)
C       INCLUDE 'O8COMM.INC'
        INCLUDE 'PROBLEM.INC'
        INTEGER I
        DOUBLE PRECISION EPS,X(*)
        EPS = 1.0d-8
c       
c       Maximum likelihood estimation Hill model mean partial derivatives
        IF(probtype.EQ.1) THEN
           DO I = 1, ndoses
              grads(I,1) = 0.0
              grads(I,2) = 0.0
              grads(I,3) = 1.0
              IF(doses(I).NE.0) THEN
                 grads(I,4) = (doses(I))**X(5) 
                 
                 grads(I,5) = X(4)*LOG(doses(I))*(doses(I))**X(5)
                 
              ELSE
                 grads(I,4) = 0.0
                 grads(I,5) = 0.0
              ENDIF
           ENDDO

        ELSE

           DO I = 1, ndoses
              grads(I,1) = 0.0
              grads(I,2) = 0.0
              grads(I,3) = 0.0
              grads(I,4) = 1.0
              IF(doses(I).NE.0) THEN
                 grads(I,5) = (doses(I))**X(6)  
                 grads(I,6) = X(5)*LOG(doses(I))*(doses(I))**X(6)
                 

              ELSE
                 grads(I,5) = 0.0
                 grads(I,6) = 0.0

              ENDIF

           ENDDO

           DoseMeanGrad(2) = 0.0
           DoseMeanGrad(3) = 0.0
           DoseMeanGrad(4) = 1.0
           IF(X(1).GT.0) THEN
              DoseMeanGrad(1)=X(5)*X(6)*( (X(1))**(X(6)-1) )
              
              
              DoseMeanGrad(5) = (X(1))**X(6)


              DoseMeanGrad(6) = X(5)*((X(1))**X(6))*LOG(X(1))
              
           ELSE
              DoseMeanGrad(1) = 0.0
              DoseMeanGrad(5) = 0.0
              DoseMeanGrad(6) = 0.0

           ENDIF

        ENDIF
        RETURN
        END

        SUBROUTINE POWCompIneq(X,GXI,I)
        INCLUDE 'O8FUCO.INC'
        INCLUDE 'PROBLEM.INC'
        INTEGER I
        DOUBLE PRECISION EPS, X(*), GXI, LK
C The following constant controls how close we can get to our bound, to
C simulate a strict inequality constraint (donlp2 enforces a '<=' 
C constraint)
        EPS = 1.0d-15

        IF (probtype .EQ. 1) THEN
C       ML Estimation
C       alpha >= EPS > 0
           IF(I .EQ. 1 .AND. constvar .EQ. 1) THEN
              GXI = X(GUNIT(2, I+NH)) - EPS
C       n >= 1, or n >= 0
           ELSEIF(I.EQ.1 + constvar) THEN
              IF(restrict.EQ.1) THEN
                 GXI = X(GUNIT(2, I+NH)) - 1
              ELSE
                 GXI = X(GUNIT(2, I+NH)) 
              ENDIF
c       n parameter <= 18
           ELSEIF (I .EQ. 2 + constvar) THEN
              GXI = 18 - X(GUNIT(2, I+NH))
c       adverse direction
           ELSEIF (I .EQ. 3 + constvar) THEN 
              IF (adverse .EQ. 1) THEN            
                 GXI = X(GUNIT(2, I+NH)) 
              ELSE 
                 GXI = -X(GUNIT(2, I+NH)) 
              ENDIF
c       rho variance parameter <= 18
           ELSEIF (I .EQ. 4 + constvar) THEN 
              GXI = 18 - X(GUNIT(2, I+NH)) 
           ENDIF
           
        ELSEIF (probtype .EQ. 2) THEN
C       Lower Confidence Limit
           IF (I .EQ. 1) THEN
c       LogLikelihood >= target
              CALL POWMEAN(X)
              CALL NegLogLike(X, LK)
C NegLogLike returns (in LK) the negative log likelihood.
              GXI = - LK - target
           ELSEIF (I .EQ. 2) THEN
c       X(1) >= 0
              GXI = X(GUNIT(2, I+NH)) - EPS
           ELSEIF (I .EQ. 3) THEN
c       BMD - X(1) >= 0
              GXI=bmd-X(1)
           ELSEIF(I .EQ. 4 .AND. constvar .eq. 1) THEN
C       alpha > 0
              GXI = X(GUNIT(2, I+NH)) - EPS
           ELSEIF(I .EQ. 4 + constvar) THEN
C       If restrict .EQ. 1, then the power parameter, n, must be >= 1
c       otherwise, n >= 0
              IF (restrict .EQ. 1) THEN
                 GXI = X(GUNIT(2, I+NH)) - 1
              ELSE
                 GXI = X(GUNIT(2, I+NH)) 
              ENDIF
           ELSEIF(I .EQ. 5 + constvar) THEN
c       n parameter <= 18
              GXI = 18 - X(GUNIT(2,I+NH))
           ELSEIF(I .EQ. 6 + constvar) THEN
c       adverse direction
              IF (adverse .EQ. 1) THEN            
                 GXI = X(GUNIT(2, I+NH)) 
              ELSE 
                 GXI = -X(GUNIT(2, I+NH)) 
              ENDIF
           ELSEIF(I .EQ. 7 + constvar) THEN
c       rho variance parameter <= 18
              GXI = 18 - X(GUNIT(2,I+NH))
           ENDIF
        ELSEIF (probtype .EQ. 3) THEN
          IF (constvar .EQ. 1 .AND. I .EQ. 1) THEN
                GXI = X(GUNIT(2, I+NH)) - EPS
          ELSEIF (I .EQ. 1 + constvar) THEN
                GXI = 18 - X(GUNIT(2, I + NH))
          ENDIF
        ELSE
C       Upper Confidence Limit
           GXI=X(GUNIT(2,I+NH))
        ENDIF
        
        RETURN
        END

        SUBROUTINE POWCompIneqGrad(X,GRAD,I)
        INCLUDE 'O8FUCO.INC'
        INCLUDE 'PROBLEM.INC'
        INTEGER I, J
        DOUBLE PRECISION X(*), GRAD(*)
        IF (probtype .EQ. 2) THEN
           IF (I .EQ. 1) THEN
              CALL POWMEAN(X)
              CALL POWPART(X)
              CALL DNegLogLike(X, GRAD)
C We are working with -NegLogLike, so need to change signs of gradient
              DO J=1, NX
                 GRAD(J) = -1.0D0*GRAD(J)
              ENDDO
           ENDIF
        ENDIF
        
        RETURN
        END

        DOUBLE PRECISION FUNCTION powerf(a,b)
c       INCLUDE 'PROBLEM.INC'
        DOUBLE PRECISION loga, bloga, a, b
        if (a .GE. .0000000001) then
           loga = log(a)
        else
           loga = -24.5258509299404572 + (2.0D10)*a
     1       -(5.0D19)*a*a
        endif
        bloga = b*loga
        if (bloga .GT. 700) then
           bloga = 700
        else
        endif

        powerf = dexp(bloga)

        RETURN
        END
