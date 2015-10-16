c     This file contains three functions that are specific to the
c     to the exponential model when doing maximum likelihood and confidence
c     limit calculations using donlp2 and the associated user 
c     functions which are non-specific to the continuous model
c     that is being run.
c     EXPFILLGUNIT fills the GUNIT array in donlp2 with the model 
c     specific equality, non-equality constraint as well as
c     the objective function.
c     EXPMEAN just calculates the mean specific to that model as a
c     vector of means, one at each dose level
c     EXPPART calculates the vector of partial derivates of the
c     mean with respect to all model parameters.
c     

      SUBROUTINE EXPFILLGUNIT
      INCLUDE 'O8COMM.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER I, J
C     
C     X IS INITIAL GUESS AND ALSO HOLDS THE CURRENT SOLUTION
C     PROBLEM DIMENSION N=DIM(X), NH=DIM(H), NG=DIM(G)
C     
C     Parameters when probtype .eq. 1:
C     X(1)  lalpha 
C     X(2)  rho
C     X(3)  a
C     X(4)  b
C     X(5)  c
C     X(6)  d
C     
C     Parameters when probtype .eq, 2:
C     X(1)  BMDL
C     X(2)  lalpha
C     X(3)  rho
C     X(4)  a
C     X(5)  b
C     X(6)  c
C     X(7)  d

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
         DO I = 0, 5
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
C     Inequality constraints: Numbers DO NOT depend on constvar because always use lalpha (unconstrained)
C     but they do depend on modtype 
C     modtype=3: model 2, no constraint on c or d (handled by equality) so NG=3;
C     modtype=4: model 3, no constraint on c so NG=4; 
C     modtype=5: model 4, no constraint on d, but on c add 1 (if increasing) or 2 (if decreasing) constraints, 
C     so NG=4 or 5
C     modtype=6: model 5, contraints on c (1 or 2) and on d (2), so NG = 6 or 7
C     Constraint # 
C     1                    rho (x(2)) <= 18           
C     2                    a (x(3)) > 0  
C     3                    b (x(4)) > 0
C     4                    d (x(6)) >= 1
C     5                    d (x(6)) <= 18
C     6                    c (x(5)) > 0 if adverse < 0 (down bad) or c > 1 if adverse > 0   
C     7                    c < 1        if adverse < 0 (down bad)
C     
C     
c     rho variance parameter <= 18
c     I: 1
         J = J + 1
         NG = NG + 1
         GUNIT(1, J) = 1
         GUNIT(2, J) = 2
         GUNIT(3, J) = -1 
c     
c     I: 2
         NG = NG + 1
         J = J + 1
         GUNIT(1, J) = 1
         GUNIT(2, J) = 3
         GUNIT(3, J) = 1         
C     
c     
c     I: 3
         NG = NG + 1
         J = J + 1
         GUNIT(1, J) = 1
         GUNIT(2, J) = 4
         GUNIT(3, J) = 1

c     d parameter >= 1 (modtype = 4 or 6)
c     I: 4
         IF (modtype .EQ. 4 .OR. modtype .EQ. 6) THEN
            J = J + 1
            NG = NG + 1
            GUNIT(1, J) = 1
            GUNIT(2, J) = 6
            GUNIT(3, J) = 1 
c     d parameter <= 18
c     I: 5
            J = J + 1
            NG = NG + 1
            GUNIT(1, J) = 1
            GUNIT(2, J) = 6
            GUNIT(3, J) = -1
         ENDIF

C     Now constraints only if models 4 or 5 (modtype=5 or 6)
         IF ((modtype .EQ. 5) .OR. (modtype .EQ. 6)) THEN       
c     c parameter lower bound (always same gunit specification regardless of lower bound 0 or 1)
c     I: 4 or 6
            J = J + 1
            NG = NG + 1
            GUNIT(1, J) = 1
            GUNIT(2, J) = 5
            GUNIT(3, J) = 1
C     c parameter upper bound if adverse direction down; I: 5 or 7
            IF (adverse .EQ. -1) THEN
               J = J + 1
               NG = NG + 1
               GUNIT(1, J) = 1
               GUNIT(2, J) = 5
               GUNIT(3, J) = -1
            ENDIF
         ENDIF

C     ENDIF
C     
C     
C     
      ELSEIF (probtype .EQ. 2) THEN
C     C       Lower Confidence Limit  
C     
C     c       There are always at least one equality constraint: imposed by 
C     c       the benchmark response level.  Additional equality
C     c       constraints are imposed by fixing specific betas to particular
C     c       values.
         NH=1
C     C       GUNIT-ARRAY, SEE DONLP2DOC.TXT
C     C       Objective function is linear in X(1)
         GUNIT(1,0)=1
         GUNIT(2,0)=1
         GUNIT(3,0)=1
C     C       First equality constraint is computed
         J=1
         GUNIT(1,J)=-1
         GUNIT(2,J)=0
         GUNIT(3,J)=0
         DO I = 0, 5
            IF (parmfixd(I) .EQ. 1) THEN
               J = J+1
               NH = NH + 1
               GUNIT(1,J)=1
               GUNIT(2,J)=I + 2
               GUNIT(3,J)=1
            ENDIF
         ENDDO
C     c *******      There are always at least three inequality constraints, since 
C     c       1)  the loglikelihood must be >= target loglikelihood
C     c       2) X(1) must be .GE. 0 and 
C     c       3) X(1) must be .LT. bmd.  
C     
C     c       NO added contraint for constant variance, since always use lnalpha now.
C     c       Other constraints are listed below, like the case for the MLE
C     c            
         NG=0
C     C Inequality constraints
c     C Constraint #            
C     C     1                     loglik >= target         
c     C     2                     BMDL > 0                 
C     C     3                     BMDL < BMD               
C     c     4                    rho (x(3)) <= 18           
C     c     5                    a (x(4)) > 0  
C     c     6                    b (x(5)) > 0
C     c     7                    d (x(7)) >= 1
C     c     8                    d (x(7)) <= 18
C     c     9                    c (x(6)) > 0 if adverse < 0 (down bad) or c > 1 if adverse > 0   
C     c     10                    c < 1        if adverse < 0 (down bad)
C     C In the comments below, I refers to the argument passed to EG()
C     c       likelihood >= target I = 1
         NG = NG + 1
         J = J + 1
         GUNIT(1, J) = -1
         GUNIT(2, J) = 0
         GUNIT(3, J) = 0
C     c       
C     c       X(1) >= 0  I = 2
         NG = NG + 1
         J = J+1
         GUNIT(1,J) = 1
         GUNIT(2,J) = 1
         GUNIT(3,J) = 1
c     c       BMD - X(1) >= 0 I = 3
         NG = NG + 1
         J = J+1
         GUNIT(1,J)=1
         GUNIT(2,J)=1
         GUNIT(3,J)=-1
C     
C     C       rho less than 18
C     C       I = 4 
         J = J + 1
         NG = NG + 1
         GUNIT(1, J) = 1
         GUNIT(2, J) = 3
         GUNIT(3, J) = -1 

C     C       a greater than 0
c     c       I = 5
         J = J + 1
         NG = NG + 1
         GUNIT(1, J) = 1
         GUNIT(2, J) = 4
         GUNIT(3, J) = 1 
C     c       
C     
C     c       b parameter > 0
C     c       I = 6
         J = J + 1
         NG = NG + 1
         GUNIT(1, J) = 1
         GUNIT(2, J) = 5
         GUNIT(3, J) = 1
         
c     c       d >= 1
C     c       I= 7
         IF (modtype .EQ. 4 .OR. modtype .EQ. 6) THEN
            J = J + 1
            NG = NG + 1
            GUNIT(1, J) = 1
            GUNIT(2, J) = 7
            GUNIT(3, J) = 1 

C     c       d  <= 18
C     c       I: 8 
            J = J + 1
            NG = NG + 1
            GUNIT(1, J) = 1
            GUNIT(2, J) = 7
            GUNIT(3, J) = -1        
         ENDIF

C     Now constraints only if models 4 or 5 (modtype=5 or 6)
         IF ((modtype .EQ. 5) .OR. (modtype .EQ. 6)) THEN       
c     c parameter lower bound (always same gunit specification regardless of lower bound 0 or 1)
C     I: 7 or 9
            J = J + 1
            NG = NG + 1
            GUNIT(1, J) = 1
            GUNIT(2, J) = 6
            GUNIT(3, J) = 1
C     c parameter upper bound if adverse direction down; I: 8 or 10
            IF (adverse .EQ. -1) THEN
               J = J + 1
               NG = NG + 1
               GUNIT(1, J) = 1
               GUNIT(2, J) = 6
               GUNIT(3, J) = -1
            ENDIF
         ENDIF


      ENDIF
      RETURN
      END
c     

      SUBROUTINE EXPMEAN(X)
C     INCLUDE 'O8COMM.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER I
      DOUBLE PRECISION X(*), varmean, delta
c     
      IF(probtype.EQ.1) THEN
         varmean = 0.0
         IF((modtype.EQ.3) .OR. (modtype.EQ.4)) THEN
            DO I = 1, ndoses
               means(I) = X(3)*EXP(adverse*(X(4)*doses(I))**X(6))
               IF(lognorm.EQ.1) means(I) = LOG(means(I))
            ENDDO
         ELSE
            DO I = 1, ndoses
               varmean = EXP(-1.0*(X(4)*doses(I))**X(6))
c     means(I) = X(3)*(X(5)-(X(5)-1)*EXP(-1.0*(X(4)*doses(I))**X(6)))
               means(I) = X(3)* (X(5)-(X(5)-1) * varmean)
               IF(lognorm.EQ.1) means(I) = LOG(means(I))
            ENDDO
         ENDIF
         varmean = 0.0
         DO I = 1, ndoses
            varmean = varmean + var(I)
         ENDDO
         varmean = varmean/ndoses
      ELSEIF(probtype.EQ.2) THEN
         IF((modtype.EQ.3) .OR. (modtype.EQ.4)) THEN
            DO I = 1, ndoses
               means(I) = X(4)*EXP(adverse*(X(5)*doses(I))**X(7))
               IF(lognorm.EQ.1) means(I) = LOG(means(I))
            ENDDO
         ELSE
            DO I = 1, ndoses
               varmean = EXP(-1.0*(X(5)*doses(I))**X(7))
c     means(I) = X(4)*(X(6)-(X(6)-1)*EXP(-1.0*(X(5)*doses(I))**X(7)))
               means(I) = X(4) * (X(6)-(X(6)-1) * varmean)
               IF(lognorm.EQ.1) means(I) = LOG(means(I))
            ENDDO
         ENDIF
         varmean = 0.0
         DO I = 1, ndoses
            varmean = varmean + var(I)
         ENDDO
         varmean = varmean/ndoses

         DoseMean = X(4)      
         IF(X(1).GT.0) THEN
            IF((modtype.EQ.3) .OR. (modtype.EQ.4)) THEN
               DoseMean = X(4)*EXP(adverse*(X(5)*X(1))**X(7))
            ELSE
               DoseMean = X(4)*(X(6)-(X(6)-1)*
     $              EXP(-1.0*(X(5)*X(1))**X(7)))
            ENDIF
         ENDIF
         IF(lognorm.EQ.1) DoseMean = LOG(DoseMean)
      ENDIF

c     
      RETURN
      END


c     SUBROUTINE EXPMEAN(X)
c     INCLUDE 'O8COMM.INC'
c     INCLUDE 'PROBLEM.INC'
c     INTEGER I
c     DOUBLE PRECISION X(*), varmean, delta
c     
c     IF(probtype.EQ.1) THEN
c     varmean = 0.0
c     IF((modtype.EQ.3) .OR. (modtype.EQ.4)) THEN
c     DO I = 1, ndoses
c     means(I) = X(3)*EXP(adverse*(X(4)*doses(I))**X(6))
c     IF(lognorm.EQ.1) means(I) = LOG(means(I))
c     ENDDO
c     ELSE
c     DO I = 1, ndoses
c     varmean = EXP(-1.0*(X(4)*doses(I))**X(6))
c     c            means(I) = X(3)*(X(5)-(X(5)-1)*EXP(-1.0*(X(4)*doses(I))**X(6)))
c     means(I) = X(3)* (X(5)-(X(5)-1) * varmean)
c     IF(lognorm.EQ.1) means(I) = LOG(means(I))
c     ENDDO
c     ENDIF
c     varmean = 0.0
c     DO I = 1, ndoses
c     varmean = varmean + var(I)
c     ENDDO
c     varmean = varmean/ndoses
c     ELSEIF(probtype.EQ.2) THEN
c     IF((modtype.EQ.3) .OR. (modtype.EQ.4)) THEN
c     DO I = 1, ndoses
c     means(I) = X(4)*EXP(adverse*(X(5)*doses(I))**X(7))
c     IF(lognorm.EQ.1) means(I) = LOG(means(I))
c     ENDDO
c     ELSE
c     DO I = 1, ndoses
c     varmean = EXP(-1.0*(X(5)*doses(I))**X(7))
c     c           means(I) = X(4)*(X(6)-(X(6)-1)*EXP(-1.0*(X(5)*doses(I))**X(7)))
c     means(I) = X(4) * (X(6)-(X(6)-1) * varmean)
c     IF(lognorm.EQ.1) means(I) = LOG(means(I))
c     ENDDO
c     ENDIF
c     varmean = 0.0
c     DO I = 1, ndoses
c     varmean = varmean + var(I)
c     ENDDO
c     varmean = varmean/ndoses
c     
c     DoseMean = X(4)      
c     IF(X(1).GT.0) THEN
c     IF((modtype.EQ.3) .OR. (modtype.EQ.4)) THEN
c     DoseMean = X(4)*EXP(adverse*(X(5)*X(1))**X(7))
c     ELSE
c     DoseMean = X(4)*(X(6)-(X(6)-1)*EXP(-1.0*(X(5)*X(1))**X(7)))
c     ENDIF
c     IF(lognorm.EQ.1) DoseMean = LOG(DoseMean)
c     ENDIF
c     ENDIF
c     
c     c           IF (X(6) .LT. 0) X(6) = 0.0
c     c           IF (X(6) .GT. 18) X(6) = 18.0
c     c           IF (X(1) .GT. bmd) X(1) = bmd
c     c           IF (risktype .EQ. 2) THEN
c     c              IF (X(4) .LT. 1.0D-20) THEN
c     c                 X(4) = 1.0D-20
c     c              ELSE
c     c              ENDIF
c     c           ENDIF
c     c           IF ((risktype .EQ. 4) .AND. (ABS(X(5)) .LT. 1.0D-20)) THEN
c     c              IF (X(5) .GE. 0) THEN
c     c                 X(5) = 1.0D-20
c     c              ELSE
c     c                 X(5) = -1.0D-20
c     c              ENDIF
c     c           ENDIF
c     c        ENDIF
c     c       
c     RETURN
c     END
c     
c     
c     EXPpart gives the partial derivatives of the exponential model
c     mean function
c     and also, in grads(ndoses+1, j) gives the partial derivatives
c     of the mean function evaluated at the estimated dose D.

      SUBROUTINE EXPPART(X)
C     INCLUDE 'O8COMM.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER I
      DOUBLE PRECISION EPS,X(*)
      EPS = 1.0d-8
c     
C     ExpMean called before this subroutine in DONLP2usrfcs – I assume that the 
c     current values of means(I) and DoseMean are available to this subroutine.
c     Maximum likelihood estimation exponential models mean partial derivatives
      IF(probtype.EQ.1) THEN
         DO I = 1, ndoses
            grads(I,1) = 0.0
            grads(I,2) = 0.0
            IF((modtype.EQ.3) .OR. (modtype.EQ.4)) THEN
               grads(I,3) = EXP(adverse*(doses(I)*X(4))**X(6))
               grads(I,4) = adverse*doses(I)*
     $              EXP(adverse*(doses(I)*X(4))**X(6))*X(3)*
     $              (doses(I)*X(4))**(X(6)-1)*X(6)
               grads(I,5) = 0.0
               grads(I,6) = 0.0
               IF((doses(I).NE.0.0) .AND. (modtype .EQ. 4)) THEN
                  grads(I,6) = adverse*
     $                 EXP(adverse*(doses(I)*X(4))**X(6))*X(3)*
     $                 (doses(I)*X(4))**X(6)*LOG(doses(I)*X(4))
               ENDIF
            ELSE
               grads(I,3) = (X(5)-(X(5)-1)*
     $              EXP(-1.0*(X(4)*doses(I))**X(6)))
c     BCA -- fixed the following 2/27/08 by adding a factor of doses(I)
c     grads(I,4) = X(3)*(X(5)-1)*EXP(-1.0*(X(4)*doses(I))**X(6))*
c     (doses(I)*X(4))**(X(6)-1)*X(6)
               grads(I,4) = X(3)*(X(5)-1)*doses(I)*
     $              EXP(-1.0*(X(4)*doses(I))**X(6))*(doses(I)*X(4))**
     $              (X(6)-1)*X(6)
               grads(I,5) = X(3)*(1.0-EXP(-1.0*(X(4)*doses(I))**X(6)))
               grads(I,6) = 0.0
               IF((doses(I).NE.0.0) .AND. (modtype .EQ. 6)) THEN
                  grads(I,6) = X(3)*(X(5)-1)*EXP(-1.0*
     $                 (X(4)*doses(I))**X(6))*(doses(I)*X(4))**(X(6))*
     $                 LOG(doses(I)*X(4))
               ENDIF
            ENDIF
            DO J = 3, 6
               IF (lognorm.EQ.1) grads(I,J) = grads(I,J)/EXP(means(I))
            ENDDO
         ENDDO

      ELSEIF(probtype.EQ.2) THEN

         DO I = 1, ndoses
            grads(I,1) = 0.0
            grads(I,2) = 0.0
            grads(I,3) = 0.0
            IF((modtype.EQ.3) .OR. (modtype.EQ.4)) THEN
               grads(I,4) = EXP(adverse*(doses(I)*X(5))**X(7))
               grads(I,5) = adverse*doses(I)*EXP(adverse*
     $              (doses(I)*X(5))**X(7))*X(4)*(doses(I)*X(5))**
     $              (X(7)-1)*X(7)
               grads(I,6) = 0.0
               grads(I,7) = 0.0
               IF((doses(I).NE.0.0) .AND. (modtype .EQ. 4)) THEN
                  grads(I,7) = adverse*EXP(adverse*(doses(I)*X(5))**
     $                 X(7))*X(4)*(doses(I)*X(5))**X(7)*LOG(doses(I)*
     $                 X(5))
               ENDIF
            ELSE
               grads(I,4) = (X(6)-(X(6)-1)*
     $              EXP(-1.0*(X(5)*doses(I))**X(7)))

               grads(I,5) = X(4)*(X(6)-1)*doses(I)*
     $              EXP(-1.0*(X(5)*doses(I))**X(7))*(doses(I)*X(5))**
     $              (X(7)-1)*X(7)
               grads(I,6) = X(4)*(1.0-EXP(-1.0*(X(5)*doses(I))**X(7)))
               grads(I,7) = 0.0
               IF((doses(I).NE.0.0) .AND. (modtype .EQ. 6)) THEN
                  grads(I,7) = X(4)*(X(6)-1)*EXP(-1.0*(X(5)*
     $                 doses(I))**X(7))*(doses(I)*X(5))**(X(7))*
     $                 LOG(doses(I)*X(5))
               ENDIF
            ENDIF
            DO J = 4, 7
               IF (lognorm.EQ.1) grads(I,J) = grads(I,J)/EXP(means(I))
            ENDDO
         ENDDO

         DoseMeanGrad(2) = 0.0
         DoseMeanGrad(3) = 0.0

         IF((modtype.EQ.3) .OR. (modtype.EQ.4)) THEN
            DoseMeanGrad(4) = EXP(adverse*(X(1)*X(5))**X(7))
         ELSE
            DoseMeanGrad(4) = (X(6)-(X(6)-1)*
     $           EXP(-1.0*(X(5)*X(1))**X(7)))
         ENDIF

         IF(X(1).GT.0) THEN
            IF((modtype.EQ.3) .OR. (modtype.EQ.4)) THEN
               DoseMeanGrad(1) = adverse*X(5)*EXP(adverse*(X(1)*
     $              X(5))**X(7))*X(4)*(X(1)*X(5))**(X(7)-1)*X(7)
            ELSE
               DoseMeanGrad(1) = X(4)*(X(6)-1)*X(5)*
     $              EXP(-1.0*(X(5)*X(1))**X(7))*(X(1)*X(5))**(X(7)-1)*
     $              X(7)
            ENDIF
            
            
            IF((modtype.EQ.3) .OR. (modtype.EQ.4)) THEN
               DoseMeanGrad(5) = adverse*X(1)*
     $              EXP(adverse*(X(1)*X(5))**X(7))*X(4)*(X(1)*X(5))**
     $              (X(7)-1)*X(7)
            ELSE
               DoseMeanGrad(5) = X(4)*(X(6)-1)*X(1)*EXP(-1.0*(X(5)*
     $              X(1))**X(7))*(X(1)*X(5))**(X(7)-1)*X(7)
            ENDIF
            
            IF((modtype.EQ.3) .OR. (modtype.EQ.4)) THEN
               DoseMeanGrad(6) = 0.0
            ELSE
               DoseMeanGrad(6) = X(4)*(1.0 -
     $              EXP(-1.0*(X(5)*X(1))**X(7)))
            ENDIF

            IF((modtype.EQ.3) .OR. (modtype.EQ.5)) THEN
               DoseMeanGrad(7) = 0.0
            ELSEIF(modtype.EQ.4) THEN
               DoseMeanGrad(7) = adverse*EXP(adverse*(X(1)*
     $              X(5))**X(7))*X(4)*(X(1)*X(5))**X(7)*LOG(X(1)*X(5))
            ELSEIF(modtype.EQ.6) THEN
               DoseMeanGrad(7) = X(4)*(X(6)-1)*EXP(-1.0*(X(5)*
     $              X(1))**X(7))*(X(1)*X(5))**(X(7))*LOG(X(1)*X(5))
            ENDIF
            
         ELSE
            DoseMeanGrad(1) = 0.0
            DoseMeanGrad(5) = 0.0
            DoseMeanGrad(6) = 0.0
            DoseMeanGrad(7) = 0.0
         ENDIF
         IF(lognorm.EQ.1) DoseMeanGrad(1) = DoseMeanGrad(1)/
     $        EXP(DoseMean)
         DO J = 4, 7
            IF(lognorm .EQ. 1) DoseMeanGrad(J) = DoseMeanGrad(J)/
     $           EXP(DoseMean)
         ENDDO
      ENDIF
      RETURN
      END

      SUBROUTINE EXPCompIneq(X,GXI,I)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER I
      DOUBLE PRECISION EPS, X(*), GXI, LK
C     The following constant controls how close we can get to our bound, to
C     simulate a strict inequality constraint (donlp2 enforces a '<=' 
C     constraint)
      EPS = 1.0d-15

      IF (probtype .EQ. 1) THEN
C     ML Estimation
C     rho <= 18
         IF(I .EQ. 1) THEN
            GXI = 18.0 - X(GUNIT(2, I+NH))
C     a > 0
         ELSEIF(I.EQ.2) THEN
            GXI = X(GUNIT(2, I+NH)) - EPS
c     b >= 0
         ELSEIF (I .EQ. 3) THEN
            GXI = X(GUNIT(2, I+NH))
c     d => 1
         ELSEIF ((I .EQ. 4) .AND. ((modtype .EQ. 4) .OR.
     $           (modtype .EQ. 6))) THEN
            GXI = X(GUNIT(2, I+NH)) - 1.0
c     d <= 18
         ELSEIF ((I .EQ. 5) .AND. ((modtype .EQ. 4) .OR.
     $           (modtype .EQ. 6))) THEN
            GXI = 18.0 - X(GUNIT(2, I+NH))
C     The rest only for modtype =5 or 6 (models 4 and 5), should only get to these I values when that is true.
            
c     c > 1 or >= 0 (adverse up or down, respectively)
         ELSEIF (((I .EQ. 4) .AND. (modtype .EQ. 5)) .OR.
     $           ((I .EQ. 6) .AND. (modtype .EQ. 6))) THEN 
            IF (adverse .EQ. 1) THEN            
               GXI = X(GUNIT(2, I+NH)) - (1.0 + EPS)
            ELSE 
               GXI = X(GUNIT(2, I+NH)) 
            ENDIF
c     c < 1 if adverse = -1 
         ELSEIF (adverse .EQ. -1 .AND. ((I .EQ. 5 .AND. modtype .EQ.
     $           5) .OR. (I .EQ. 7 .AND. modtype .EQ. 6))) THEN 
            GXI = (1.0 - EPS) - X(GUNIT(2, I+NH)) 
         ENDIF
         
      ELSEIF (probtype .EQ. 2) THEN
C     Lower Confidence Limit                    
         IF (I .EQ. 1) THEN
c     LogLikelihood >= target
            CALL EXPMEAN(X)
            CALL NegLogLike(X, LK)
C     NegLogLike returns (in LK) the negative log likelihood.
            GXI = - LK - target
         ELSEIF (I .EQ. 2) THEN
c     X(1) > 0
            GXI = X(GUNIT(2, I+NH)) - EPS
         ELSEIF (I .EQ. 3) THEN
c     BMD - X(1) >= 0
            GXI=bmd-X(1)
         ELSEIF(I .EQ. 4) THEN
C     rho <= 18
            GXI = 18.0 - X(GUNIT(2, I+NH))
         ELSEIF(I .EQ. 5) THEN
C     a > 0
            GXI = X(GUNIT(2, I+NH)) - EPS
         ELSE IF(I.EQ.6) THEN
c     b parameter >= 0
            GXI = X(GUNIT(2,I+NH))
         ELSE IF((I .EQ. 7).AND.
     $           ((modtype.EQ.4).OR.(modtype.EQ.6))) THEN
c     d constraints >= 1:
            GXI = X(GUNIT(2,I+NH)) - 1
         ELSE IF ((I .EQ. 8) .AND. ((modtype .EQ. 4) .OR.
     $           (modtype .EQ. 6))) THEN
c     d constraints <= 18:
            GXI = 18.0 - X(GUNIT(2, I+NH))
C     The rest only for modtype =5 or 6 (models 4 and 5), should only get to these I values when that is true.
            
c     c > 1 or >= 0 (adverse up or down, respectively)
         ELSE IF (((I .EQ. 7) .AND. (modtype .EQ. 5)) .OR.
     $           ((I .EQ. 9) .AND. (modtype .EQ. 6))) THEN 
            IF (adverse .EQ. 1) THEN            
               GXI = X(GUNIT(2, I+NH)) - (1.0 + EPS)
            ELSE 
               GXI = X(GUNIT(2, I+NH)) 
            ENDIF
c     c < 1 if adverse = -1 
         ELSE IF (adverse .EQ. -1 .AND. ((I .EQ. 8 .AND.
     $           modtype .EQ. 5) .OR. (I .EQ. 10 .AND.
     $           modtype .EQ. 6))) THEN 
            GXI = (1.0 - EPS) - X(GUNIT(2, I+NH))       

         ENDIF

      ELSE IF (probtype .EQ. 3) THEN
         IF (constvar .EQ. 1 .AND. I .EQ. 1) THEN
            GXI = X(GUNIT(2, I+NH)) - EPS
         ELSE IF (I .EQ. 1 + constvar) THEN
            GXI = 18 - X(GUNIT(2, I + NH))
         ENDIF
      ELSE
C     Upper Confidence Limit
         GXI=X(GUNIT(2,I+NH))
      ENDIF
      
      RETURN
      END

      SUBROUTINE EXPCompIneqGrad(X,GRAD,I)
      INCLUDE 'O8FUCO.INC'
      INCLUDE 'PROBLEM.INC'
      INTEGER I, J
      DOUBLE PRECISION X(*), GRAD(*)
      IF (probtype .EQ. 2) THEN
         IF (I .EQ. 1) THEN
            CALL EXPMEAN(X)
            CALL EXPPART(X)
            CALL DNegLogLike(X, GRAD)
C     We are working with -NegLogLike, so need to change signs of gradient
            DO J=1, NX
               GRAD(J) = -1.0D0*GRAD(J)
            ENDDO
         ENDIF
      ENDIF
      
      RETURN
      END

      DOUBLE PRECISION FUNCTION powerf(a,b)
c     INCLUDE 'PROBLEM.INC'
      DOUBLE PRECISION loga, bloga, a, b
      if (a .GE. .0000000001) then
         loga = log(a)
      else
         loga = -24.5258509299404572 + (2.0D10)*a
     1        -(5.0D19)*a*a
      endif
      bloga = b*loga
      if (bloga .GT. 700) then
         bloga = 700
      else
      endif

      powerf = dexp(bloga)

      RETURN
      END
