c		This file contains three functions that are specific to the
c		to the poly model when doing maximum likelihood and confidence
c		limit calculations using donlp2 and the associated user 
c		functions which are non-specific to the continuous model
c		that is being run.
c		FILLGUNIT fills the GUNIT array in donlp2 with the model 
c			specific equality, non-equality constraint as well as
c			the objective function.
c		HILLMEAN just calculates the mean specific to that model as a
c			vector of means, one at each dose level
c		HILLPART calculates the vector of partial derivates of the
c			mean with respect to all model parameters.

	SUBROUTINE POLYFILLGUNIT
	INCLUDE 'O8COMM.INC'
	INCLUDE 'PROBLEM.INC'
	INTEGER I, J
C       
C       X IS INITIAL GUESS AND ALSO HOLDS THE CURRENT SOLUTION
C       PROBLEM DIMENSION N=DIM(X), NH=DIM(H), NG=DIM(G)
C       probtype    Problem
C          1        Maximum likelihood estimation
C          2        Lower confidence limit on BMD
C          3        (not used here) Fit Model A3
C          4        Profile Likelihood
C          5        (not yet implemented) Upper Confidence
C                   limit for BMD
C       
C   Parameters when probtype .eq. 1 or 4:
C   X(1) alpha (constvar .eq. 1) or lalpha (constvar .eq. 0)
C   X(2) rho
C   X(3) beta(0)
C   X(4) beta(1)
C        .
C        .
C        .
C   X(K+3) beta(K)
C   Probtype .eq. 2
C   X(1) alpha (constvar .eq. 1) or lalpha (constvar .eq. 0)
C   X(2) rho
C   X(3) BMDL
C   X(4) beta(0)
C   X(5) beta(1)
C        .
C        .
C        .
C   X(K+4) beta(K)
	IF(probtype .EQ. 1) THEN
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
	   DO I = 0, nparm-1
	      IF (parmfixd(I) .EQ. 1) THEN
		 J = J + 1
		 NH = NH + 1
		 GUNIT(1, J) = 1
		 GUNIT(2, J) = I + 1
		 GUNIT(3, J) = 1
	      ENDIF
	   ENDDO

	   NG = 0
C Inequality constraints: Numbers depend on constvar
C Constraint #        Constvar .eq. 1          Constvar .eq. 0
C     1                 alpha > 0                rho < 18
C     2                 rho < 18                 
C       Remaining constraints happen if polynomial coefficients
C       are constrained to be positive or negative (depends on the
C       sign of restrict.
C
C       If constvar .eq. 1, then 
C	Variance coefficient parameter, alpha, must be positive
C	else 
C       rho <= 18.  This is
C       I = 1
	   NG = NG + 1
	   J = J + 1
	   GUNIT(1, J) = 1
           IF (constvar .EQ. 1) THEN
	      GUNIT(2, J) = 1
	      GUNIT(3, J) = 1 
	   ELSE
	      GUNIT(2, J) = 2
	      GUNIT(3, J) = -1 
	   ENDIF
	   
c	If there is a restriction, then bound the parameters
c       
	   IF (restrict .NE. 0) THEN
c	Restricted to be positive or negative
	      DO I = 4, nparm
                 J = J + 1
                 NG = NG + 1
                 GUNIT(1, J) = 1
                 GUNIT(2, J) = I
                 IF (restrict .EQ. -1) THEN
                    GUNIT(3, J) = -1 
                 ELSE
                    GUNIT(3, J) = 1
                 ENDIF
	      ENDDO
	   ENDIF       
	ELSEIF (probtype .EQ. 2) THEN
C       Lower Confidence Limit
C       
c       There is always at least one equality constraint, imposed by 
C       the benchmark response level.  Additional equality
c       constraints are imposed by fixing specific betas to particular
c       values.
C       GUNIT-ARRAY, SEE DONLP2DOC.TXT
C       Objective function is linear in X(1)
	   GUNIT(1,0)=1
	   GUNIT(2,0)=1
	   GUNIT(3,0)=1
C       First equality constraint is computed
	   NH = 0
           J = 0
           
           NH = NH + 1
           J = J + 1
           GUNIT(1,J)=-1
           GUNIT(2,J)=0
           GUNIT(3,J)=0

	   DO I = 0, nparm - 1
	      IF (parmfixd(I) .EQ. 1) THEN
		 J = J + 1
		 NH = NH + 1
		 GUNIT(1,J)=1
		 GUNIT(2,J)=I + 2
		 GUNIT(3,J)=1
	      ENDIF
	   ENDDO
c       There are always at least three inequality constraints:
C        1) loglikelihood must be >= target
C        2) X(1) must be .GE. 0, and
C        3) X(1) must be .LT. BMD

c       In addition, if the variance is constant,
c       then X(2) must be .GE. 0.  Also, X(3) must be .LE. 18
C       Finally, if restrict .EQ. 1,
c       then there are an additional inequality constraints.     
	   NG=0
C       loglikelihood >= target : I = 1
           NG = NG + 1
           J = J + 1
           GUNIT(1, J) = -1
           GUNIT(2, J) = 0
           GUNIT(3, J) = 0
c	X(1) >= 0: I = 2
           NG = NG + 1
	   J = J + 1
	   GUNIT(1,J) = 1
	   GUNIT(2,J) = 1
	   GUNIT(3,J) = 1
c       BMD - X(1) >= 0: I = 3
           NG = NG + 1
	   J = J + 1
	   GUNIT(1,J)=1
	   GUNIT(2,J)=1
	   GUNIT(3,J)=-1

c       restrict rho <= 18: I = 4

           NG = NG + 1
           J = J + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 3
           GUNIT(3, J) = -1 

C	Variance coefficient parameter, alpha, must be positive
C       I = 5 (if)
	   IF(constvar .EQ. 1) THEN
	      NG = NG + 1
	      J = J + 1
	      GUNIT(1, J) = 1
	      GUNIT(2, J) = 2
	      GUNIT(3, J) = 1 
	   ENDIF

c	If there is a restriction, then bound the parameters
c       
	   IF(restrict.NE.0) THEN
c       I = Index + 4 + constvar for 
c	Restricted to be positive or negative
	      DO I = 5, nparm
		 IF(parmfixd(I-2).NE.1) THEN
		    J = J + 1
		    NG = NG + 1
		    GUNIT(1, J) = 1
		    GUNIT(2, J) = I
		    IF (restrict .EQ. -1) THEN
		       GUNIT(3, J) = -1 
		    ELSE
		       GUNIT(3, J) = 1
		    ENDIF 
		 ENDIF
	      ENDDO

	   ENDIF
C       
	ELSEIF (probtype .EQ. 4) THEN
C       Profile Calc.
C       Like probtype .eq. 1, with one more equality constraint.
C       
C       Objective function is a function of all the parameters
	   
	   GUNIT(1,0) = -1
	   GUNIT(2,0) = 0
	   GUNIT(3,0) = 0
C       constraints imposed by presetting some
C       parameters, as determined by parmfixed(i)
	   NH = 0
	   J = 0
	   DO I = 0, nparm-1
	      IF (parmfixd(I) .EQ. 1) THEN
		 J = J + 1
		 NH = NH + 1
		 GUNIT(1, J) = 1
		 GUNIT(2, J) = I + 1
		 GUNIT(3, J) = 1
	      ENDIF
	   ENDDO
C       For the Profile Calcs. there is a bmr equality constraint also
           NH = NH+1
           J = J+1
           GUNIT(1, J) = -1
           GUNIT(2, J) = 0
           GUNIT(3, J) = 0

C       
	   NG = 0
C Inequality constraints: Numbers depend on constvar
C Constraint #        Constvar .eq. 1          Constvar .eq. 0
C     1                 alpha > 0                rho < 18
C     2                 rho < 18                 
C       Remaining constraints happen if polynomial coefficients
C       are constrained to be positive or negative (depends on the
C       sign of restrict.
C
C	Variance coefficient parameter, alpha, must be positive
C	
           IF (constvar .EQ. 1) THEN
              IF(parmfixd(0).NE.1) THEN
                 NG = NG + 1
                 J = J + 1
                 GUNIT(1, J) = 1
                 GUNIT(2, J) = 1
                 GUNIT(3, J) = 1 
              ENDIF
           ENDIF

c       restrict rho <= 18: I = 1 + constvar
           NG = NG + 1
           J = J + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 2
           GUNIT(3, J) = -1 
	   
c	If there is a restriction, then bound the parameters
c       
	   IF (restrict .NE. 0) THEN
c	Restricted to be positive or negative
	      DO I = 4, nparm
                 J = J + 1
                 NG = NG + 1
                 GUNIT(1, J) = 1
                 GUNIT(2, J) = I
                 IF (restrict .EQ. -1) THEN
                    GUNIT(3, J) = -1 
                 ELSE
                    GUNIT(3, J) = 1
                 ENDIF
	      ENDDO
	   ENDIF       

	ENDIF


	
	RETURN
	END
c       
	SUBROUTINE POLYMEAN(X)
C	INCLUDE 'O8COMM.INC'
	INCLUDE 'PROBLEM.INC'
	INTEGER I, J , K
	DOUBLE PRECISION X(*)
c	OPEN(11,FILE='CheckMeans.txt',STATUS='UNKNOWN')
c       
	IF ((probtype .EQ. 1).OR.(probtype .EQ. 4)) THEN
	   DO I = 1, ndoses
	      means(I) = X(nparm)
	      DO J = nparm-1, 3, -1
		 means(I) = means(I)*doses(I) + X(J)
	      ENDDO
	   ENDDO
	   IF (probtype .EQ. 4) THEN
	      bmdmean = X(nparm)
	      DO J = nparm-1, 3, -1
		 bmdmean = bmdmean*bmd + X(J)
	      ENDDO
	   ENDIF
	ELSE
	   DO I = 1, ndoses
	      means(I) = X(nparm)
	      DO J = nparm-1, 4, -1
		 means(I) = means(I)*doses(I) + X(J)
	      ENDDO
	   ENDDO
C       
	   DoseMean = X(nparm)
	   DO J = nparm-1, 4, -1
	      DoseMean = DoseMean*X(1) + X(J)
	   ENDDO
	ENDIF

	RETURN
	END
c	
c       
c	polypart gives the partial derivatives of the poly model mean function
c	and also, in grads(ndoses+1, j) gives the partial derivatives
c	of the mean function evaluated at the estimated dose D.
	SUBROUTINE POLYPART(X)
C	INCLUDE 'O8COMM.INC'
	INCLUDE 'PROBLEM.INC'
	INTEGER I,J
	DOUBLE PRECISION EPS,X(*)
	EPS = 1.0d-8
c	OPEN(21,FILE='CheckMeans.txt',STATUS='UNKNOWN')
	
c	Maximum likelihood estimation and Profile calculation for polynomial
c       model mean partial derivatives
	IF(probtype.EQ.1 .OR. probtype .EQ. 4) THEN
	   DO I = 1, ndoses
	      grads(I,1) = 0.0
	      grads(I,2) = 0.0
	      grads(I,3) = 1.0
	      DO J = 4, nparm
		 grads(I, J) = doses(I)**(J-3)
	      ENDDO
	   ENDDO
	   IF (probtype .EQ. 4) THEN
	      bmdmeangrad(1) = 0.0
	      bmdmeangrad(2) = 0.0
	      bmdmeangrad(3) = 1.0
	      DO I = 4,nparm
		 bmdmeangrad(I) = bmd**(I-3)
	      ENDDO 		
	   ENDIF
	ELSE
	   DO I = 1, ndoses
	      grads(I,1) = 0.0
	      grads(I,2) = 0.0
	      grads(I,3) = 0.0
	      grads(I,4) = 1.0
	      DO J = 5, nparm
		 grads(I, J) = doses(I)**(J-4)
	      ENDDO
	   ENDDO
c	Gradient of the mean function at X(1) with respect to all parameters
	   DoseMeanGrad(1) = X(5)		
	   DO I = 6, nparm
	      DoseMeanGrad(1) = 
     1          DoseMeanGrad(1) + (I-4)*X(I)*(X(1)**(I-5))
	   ENDDO
	   DoseMeanGrad(2) = 0.0
	   DoseMeanGrad(3) = 0.0
	   DoseMeanGrad(4) = 1.0
	   DO I = 5,nparm
	      DoseMeanGrad(I) = X(1)**(I-4)
	   ENDDO 	
	   
	ENDIF
	RETURN
	END

	SUBROUTINE PolyCompIneq(X,GXI,I)
	INCLUDE 'O8FUCO.INC'
	INCLUDE 'PROBLEM.INC'
	INTEGER I
	INTEGER I_EG
	DOUBLE PRECISION EPS, X(*), GXI, LK
	I_EG = I + NH
	EPS = 1.0D-15
	IF (probtype .EQ. 1) THEN
	   IF (I .EQ. 1 ) THEN
	      IF (constvar .EQ. 1) THEN
		 GXI = X(GUNIT(2, I_EG)) - EPS
	      ELSE
		 GXI = 18 - X(GUNIT(2, I_EG))
	      ENDIF
	   ELSE 
	      GXI = restrict * X(GUNIT(2, I_EG))
	   ENDIF
	ELSEIF (probtype .EQ. 2) THEN
C       Lower Confidence Limit
	   IF (I .EQ. 1) THEN
C       loglikelihood >= target
	      CALL polymean(X)
	      CALL NegLogLike(X, LK)
	      GXI = - LK - target
	   ELSEIF (I .EQ. 2) THEN
	      GXI = X(1) - EPS
	   ELSEIF(I .EQ. 3) THEN
	      GXI = bmd - X(1)
	   ELSEIF (I .EQ. 4) THEN
	      GXI = 18 - X(3)
	   ELSEIF (I .EQ. 5 .AND. constvar .EQ. 1) THEN
	      GXI = X(2) - EPS
	   ELSE
	      GXI = X(GUNIT(2, I_EG)) * restrict
	   ENDIF
	ELSEIF (probtype .EQ. 3) THEN
	   IF (constvar .EQ. 1 .AND. I .EQ. 1) THEN
	      GXI = X(1) - EPS
	   ELSEIF (I .EQ. 1 + constvar) THEN
		 GXI = 18 - X(2)
	   ENDIF
	ELSEIF(probtype .EQ. 4) THEN
	   IF (constvar .EQ. 1) THEN
	      IF (I .EQ. 1) THEN
		 GXI = X(1) - EPS
	      ENDIF
	   ELSEIF (I .EQ. 1 + constvar) THEN
	      GXI = 18 - X(2)
	   ELSE 
	      GXI = restrict * X(GUNIT(2, I + constvar + NH))
	   ENDIF

	ENDIF

	RETURN
	END

	SUBROUTINE PolyCompIneqGrad(X,GRAD,I)
	INCLUDE 'O8FUCO.INC'
	INCLUDE 'PROBLEM.INC'
	INTEGER I, J
	DOUBLE PRECISION X(*), GRAD(*)

C  We need to compute gradients in the following case:
C  probtype .EQ. 2: I .EQ. 1

	IF (probtype .EQ. 2 .AND. I .EQ. 1) THEN
	   CALL polymean(X)
	   CALL polypart(X)
	   CALL DNegLogLike(X, GRAD)
	   DO J = 1, NX
	      GRAD(J) = -1.0D0*GRAD(J)
	   ENDDO
	ENDIF
	RETURN
	END

