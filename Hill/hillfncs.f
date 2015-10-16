c       This file contains three functions that are specific to the
c       to the hill model when doing maximum likelihood and confidence
c       limit calculations using donlp2 and the associated user 
c       functions which are non-specific to the continuous model
c       that is being run.
c       FILLGUNIT fills the GUNIT array in donlp2 with the model 
c       specific equality, non-equality constraint as well as
c       the objective function.
c       HILLMEAN just calculates the mean specific to that model as a
c       vector of means, one at each dose level
c       HILLPART calculates the vector of partial derivates of the
c       mean with respect to all model parameters.
c       

	SUBROUTINE HILLFILLGUNIT
	INCLUDE 'O8COMM.INC'
	INCLUDE 'PROBLEM.INC'
	INTEGER I, J
C       
C       X IS INITIAL GUESS AND ALSO HOLDS THE CURRENT SOLUTION
C       PROBLEM DIMENSION N=DIM(X), NH=DIM(H), NG=DIM(G)
C       
C       Parameters when probtype .eq. 1:
C       X(1)  alpha (constvar .eq.1) or lalpha (constvar .eq. 0)
C       X(2)  rho
C       X(3)  intercept
C       X(4)  v
C       X(5)  n (power parameter)
C       X(6)  k
C       
C       Parameters when probtype .eq, 2:
C       X(1)  BMDL
C       X(2)  alpha (constvar .eq.1) or lalpha (constvar .eq. 0)
C       X(3)  rho
C       X(4)  intercept
C       X(5)  v
C       X(6)  n (power parameter)
C       X(7)  k
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
	   DO I = 0, nparm-1
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
C       Inequality constraints: Numbers depend on constvar
C       Constraint #            Constvar .EQ. 1          Constvar .EQ. 0
C       1                      alpha > 0               n > 0 
C       (or >= 1, 
C       depending on restrict)
C       2                    n > 0 (or >=1, depending) k > 0
C       3                    k > 0                     n <= 18
C       4                    n <= 18                   rho <= 18
C       5                    rho <= 18                

C	Variance coefficient parameter, alpha, must be positive
C
           IF (constvar .EQ. 1) THEN
              NG = NG + 1
              J = J + 1
              GUNIT(1, J) = 1
              GUNIT(2, J) = 1
              GUNIT(3, J) = 1 
           ENDIF
C       If restrict .EQ. 1, then the power parameter, n, must be >= 1
c	Otherwise, n >= 0
           NG = NG + 1
           J = J + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 5
           GUNIT(3, J) = 1         
C	
C	K parameter is greater than 0
           NG = NG + 1
           J = J + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 6
           GUNIT(3, J) = 1

c	n parameter <= 18
           J = J + 1
           NG = NG + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 5
           GUNIT(3, J) = -1 

c	rho variance parameter <= 18
           J = J + 1
           NG = NG + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 2
           GUNIT(3, J) = -1 
C       
C       
	ELSEIF (probtype .EQ. 2) THEN
C       Lower Confidence Limit
C       
c       There is always at least one equality constraint: imposed by
c       the benchmark response level.  Additional equality 
c       constraints are imposed by fixing specific betas to particular
c       values.

C       GUNIT-ARRAY, SEE DONLP2DOC.TXT
C       Objective function is linear in X(1)
	   GUNIT(1,0)=1
	   GUNIT(2,0)=1
	   GUNIT(3,0)=1
C       First equality constraint is computed: BMR constraint
	   NH=1
	   J = 1
	   GUNIT(1,J)=-1
	   GUNIT(2,J)=0
	   GUNIT(3,J)=0
C       Equality constraints imposed by fixed parameters
	   DO I = 0, nparm - 1
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
C       Inequality constraints
C       Constraint #            Constvar .EQ. 1          Constvar .EQ. 0
C       1                     loglik >= target         loglik >= target
C       2                     BMDL > 0                 BMDL > 0 
C       3                     BMDL < BMD               BMDL < BMD
C       4                      alpha > 0               n > 0 (or >= 1, 
C                                                      depending on restrict)
C       5                    n > 0 (or >=1, depending) k > 0
C       6                    k > 0                     n <= 18
C       7                    n <= 18                   rho <= 18
C       8                    rho <= 18                
c       
C       In the comments below, I refers to the argument passed to EG()
c       likelihood >= target I = 1
           NG = NG + 1
           J = J + 1
           GUNIT(1, J) = -1
           GUNIT(2, J) = 0
           GUNIT(3, J) = 0
c       
c	X(1) >= 0  I = 2
           NG = NG + 1
	   J = J + 1
	   GUNIT(1,J) = 1
	   GUNIT(2,J) = 1
	   GUNIT(3,J) = 1
c       
c       BMD - X(1) >= 0  I = 3
           NG = NG + 1
	   J = J + 1
	   GUNIT(1,J)=1
	   GUNIT(2,J)=1
	   GUNIT(3,J)=-1
c       
C	Variance coefficient parameter, alpha, must be positive
C       I = 4 (if constvar .eq. 1)
	   IF(constvar.EQ.1) THEN
              NG = NG + 1
	      J = J + 1
	      GUNIT(1, J) = 1
	      GUNIT(2, J) = 2
	      GUNIT(3, J) = 1 
	   ENDIF
C       If restrict .EQ. 1, then the power parameter, n, must be >= 1
c	otherwise, n >= 0       I = 4 + constvar 
           J = J + 1
           NG = NG + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 6
           GUNIT(3, J) = 1 
c       K > 0     I = 5 + constvar
           J = J + 1
           NG = NG + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 7
           GUNIT(3, J) = 1        
c	n parameter <= 18         I = 6 + constvar
           J = J + 1
           NG = NG + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 6
           GUNIT(3, J) = -1 
c	rho variance parameter <= 18       I = 7 + constvar
           J = J + 1
           NG = NG + 1
           GUNIT(1, J) = 1
           GUNIT(2, J) = 3
           GUNIT(3, J) = -1 
	   ENDIF
	   RETURN
	   END
c       
	SUBROUTINE HILLMEAN(X)
	INCLUDE 'PROBLEM.INC'
	INTEGER I
	REAL TEMP1, TEMP2, Kparm
	DOUBLE PRECISION X(*)
c       
	IF(probtype.EQ.1) THEN
	   IF(X(6).GT.0) THEN
	      Kparm = X(6)
	   ELSE
	      Kparm = .00000001
	   ENDIF

	   DO I = 1, ndoses
	      means(I) = X(3)
	      IF(doses(I).NE.0) THEN
		 IF((doses(I)/Kparm).LT.1) THEN
		    TEMP1 = (doses(I)/Kparm)**X(5)
		    TEMP2 = (doses(I)/Kparm)**X(5) + 1
		 ELSE
		    TEMP1 = 1.0
		    TEMP2 = (Kparm/doses(I))**X(5) + 1
		 ENDIF
		 means(I) = means(I) + X(4)*TEMP1/TEMP2

	      ENDIF
	   ENDDO
	ELSE
	   IF (X(1) .LE. 1.0D-20) X(1) = 1.0D-20
	   IF (X(7) .LE. 0) X(7) = .00000001
	   IF (X(6) .LT. -18) X(6) = -18.0
	   IF ((risktype .EQ. 2).AND.(X(4) .LT. 1.0D-20)) X(4) = 1.0D-20
	   IF ((risktype .EQ. 4).AND.(ABS(X(5)) .LT. 1.0D-20)) THEN
	      IF (X(5) .GE. 0) THEN
		 X(5) = 1.0D-20
	      ELSE
		 X(5) = -1.0D-20
	      ENDIF
	   ENDIF

	   IF(X(7).LE.0) THEN
	      Kparm = .00000001
	   ELSE
	      Kparm = X(7)
	   ENDIF
	   DO I = 1, ndoses
	      means(I) = X(4)
	      IF(doses(I).NE.0) THEN
		 IF ((((doses(I)/Kparm).LT.1) .AND. (X(6) .GT.0))
     1      .OR. (((doses(I)/Kparm).GT.1) .AND. (X(6) .LT.0)) ) THEN
		    TEMP1 = (doses(I)/Kparm)**X(6)
		    TEMP2 = (doses(I)/Kparm)**X(6) + 1
		 ELSE
		    TEMP1 = 1.0
		    TEMP2 = (Kparm/doses(I))**X(6) + 1
		 ENDIF
		 means(I) = means(I) + X(5)*TEMP1/TEMP2
	      ENDIF
	      
	   ENDDO
	   DoseMean = X(4)	
	   IF(X(1).GT.0) THEN
	      IF ( (((X(1)/Kparm).LT.1) .AND. (X(6) .GT.0)) 
     1	   .OR. (((X(1)/Kparm).GT.1) .AND. (X(6) .LT.0)) )THEN
		 TEMP1 = (X(1)/Kparm)**X(6)
		 TEMP2 = (X(1)/Kparm)**X(6) + 1
	      ELSE
		 TEMP1 = 1.0
		 TEMP2 = (Kparm/X(1))**X(6) + 1
	      ENDIF
	      DoseMean = DoseMean + X(5)*TEMP1/TEMP2
	   ENDIF
	ENDIF
c       
	RETURN
	END
c	
c       
c	Hillpart gives the partial derivatives of the Hill model mean function
c	and also, in grads(ndoses+1, j) gives the partial derivatives
c	of the mean function evaluated at the estimated dose D.
	SUBROUTINE HILLPART(X)
C	INCLUDE 'O8COMM.INC'
	INCLUDE 'PROBLEM.INC'
	INTEGER I, ChangeK
	DOUBLE PRECISION TEMP,TEMP2,EPS,X(*),Kparm
	EPS = 1.0d-8
c       
c	Maximum likelihood estimation Hill model mean partial derivatives
	IF(probtype.EQ.1) THEN
	   IF(X(6).LE.0) THEN
	      Kparm = .00000001
	      ChangeK = 1
	   ELSE
	      Kparm = X(6)
	      ChangeK = 0
	   ENDIF
	   DO I = 1, ndoses
	      grads(I,1) = 0.0
	      grads(I,2) = 0.0
	      grads(I,3) = 1.0
	      IF(doses(I).NE.0) THEN
		 TEMP = doses(I)/Kparm
		 TEMP2 = Kparm/doses(I)
		 IF(TEMP.LT.1) THEN
		    grads(I,4) = TEMP**X(5)/(1+TEMP**X(5))
		 ELSE
		    grads(I,4) = 1/(TEMP2**X(5) + 1)
		 ENDIF
		 grads(I,5) = X(4)*(doses(I)*Kparm)**X(5)
		 grads(I,5) = grads(I,5)*(LOG(doses(I)) - LOG(Kparm))
		 grads(I,5) = grads(I,5)/(doses(I)**X(5) + Kparm**X(5))**2
		 grads(I,6) = -X(4)*X(5)*(doses(I)*Kparm)**X(5)/Kparm
		 grads(I,6) = grads(I,6)/(doses(I)**X(5) + Kparm**X(5))**2
	      ELSE
		 grads(I,4) = 0.0
		 grads(I,5) = 0.0
		 grads(I,6) = 0.0
	      ENDIF
	      IF(ChangeK.EQ.1) THEN
		 grads(I,6) = 0.0
	      ENDIF
	   ENDDO
	ELSE
	   
	   IF(X(7).LE.0) THEN
	      Kparm = .00000001
	      ChangeK = 1
	   ELSE
	      Kparm = X(7)
	      ChangeK = 0
	   ENDIF
	   
	   DO I = 1, ndoses
	      grads(I,1) = 0.0
	      grads(I,2) = 0.0
	      grads(I,3) = 0.0
	      grads(I,4) = 1.0
	      IF(doses(I).NE.0) THEN
		 TEMP = doses(I)/Kparm
		 TEMP2 = Kparm/doses(I)
		 IF(TEMP.LT.1) THEN
		    grads(I,5) = TEMP**X(6)/(1+TEMP**X(6))
		 ELSE
		    grads(I,5) = 1/(TEMP2**X(6) + 1)
		 ENDIF
		 grads(I,6) = X(5)*(doses(I)*Kparm)**X(6)
		 grads(I,6) = grads(I,6)*(LOG(doses(I)) - LOG(Kparm))
		 grads(I,6) = grads(I,6)/(doses(I)**X(6) + Kparm**X(6))**2
		 grads(I,7) = -X(5)*X(6)*(doses(I)*Kparm)**X(6)/Kparm
		 grads(I,7) = grads(I,7)/(doses(I)**X(6)+Kparm**X(6))**2

	      ELSE
		 grads(I,5) = 0.0
		 grads(I,6) = 0.0
		 grads(I,7) = 0.0
	      ENDIF
	      IF(ChangeK.EQ.1) THEN
		 grads(I,7) = 0.0
	      ENDIF
	   ENDDO
	   TEMP = X(1)/Kparm
	   TEMP2 = Kparm/X(1)

	   DoseMeanGrad(2) = 0.0
	   DoseMeanGrad(3) = 0.0
	   DoseMeanGrad(4) = 1.0

	   IF( (X(1)**X(6)+Kparm**X(6)) .GT. 1.0D-15) THEN
	      DoseMeanGrad(1)=((Kparm)**X(6))/(X(1)**X(6)+Kparm**X(6))**2
	      DoseMeanGrad(1)=DoseMeanGrad(1)*X(5)*X(6)*(X(1)**(X(6)-1))
	      IF(TEMP.LT.1) THEN

		 DoseMeanGrad(5) = TEMP**X(6)/(1+TEMP**X(6))

	      ELSE
		 DoseMeanGrad(5) = 1/(TEMP2 + 1)
	      ENDIF

	      DoseMeanGrad(6)=(Kparm**X(6))/(X(1)**X(6)+Kparm**X(6))**2
	      DoseMeanGrad(6)=DoseMeanGrad(6)*X(5)*(X(1)**X(6))
	      DoseMeanGrad(6)=DoseMeanGrad(6)*(LOG(X(1))-LOG(Kparm))
	      DoseMeanGrad(7)=
     1          -(Kparm**(X(6)-1))/(X(1)**X(6)+Kparm**X(6))**2
	      DoseMeanGrad(7)=DoseMeanGrad(7)*X(5)*X(6)*(X(1)**X(6))

	   ELSE
	      DoseMeanGrad(1) = ((Kparm)**X(6))/((1.0D-15)**2)
	      DoseMeanGrad(1) = 
     1          DoseMeanGrad(1)*X(5)*X(6)*(X(1)**(X(6)-1))
	      DoseMeanGrad(6)=(Kparm**X(6))/(1.0D-15)**2
	      DoseMeanGrad(6)=DoseMeanGrad(6)*X(5)*(X(1)**X(6))
	      DoseMeanGrad(6)=DoseMeanGrad(6)*(LOG(X(1))-LOG(Kparm))
	      DoseMeanGrad(7)=-(Kparm**(X(6)-1))/(1.0D-15)**2
	      DoseMeanGrad(7)=DoseMeanGrad(7)*X(5)*X(6)*(X(1)**X(6))
	      

	   ENDIF
	   
	   IF(ChangeK.EQ.1) THEN
	      DoseMeanGrad(7) = 0.0
	   ENDIF
	ENDIF
	RETURN
	END

	SUBROUTINE HillCompIneq(X,GXI,I)
	INCLUDE 'O8FUCO.INC'
	INCLUDE 'PROBLEM.INC'
	INTEGER I
	DOUBLE PRECISION EPS, X(*), GXI, LK
C The following constant controls how close we can get to our bound, to
C simulate a strict inequality constraint (donlp2 enforces a '<=' 
C constraint)
	EPS = 1.0D-15

	IF (probtype .EQ. 1) THEN
C       ML Estimation
C       alpha >= EPS > 0
	   IF(constvar .EQ. 1 .AND. I .EQ. 1) THEN
	      GXI = X(GUNIT(2, I+NH)) - EPS
C       n >= 1, or n >= 0
	   ELSEIF(I .EQ. 1 + constvar) THEN
	      IF(restrict.EQ.1) THEN
		 GXI = X(GUNIT(2, I+NH)) - 1
	      ELSE
		 GXI = X(GUNIT(2, I+NH))
	      ENDIF
           ELSEIF (I .EQ. 2 + constvar) THEN
C       K >= EPS > 0
	      GXI = X(GUNIT(2, I+NH)) - EPS
C       n <= 18
	   ELSEIF (I .EQ. 3 + constvar) THEN
	      GXI = 18 - X(GUNIT(2, I+NH))
	      
C       rho <= 18
	   ELSEIF (I .EQ. 4 + constvar) THEN
	      GXI = 18 - X(GUNIT(2, I+NH))
	   ENDIF

	ELSEIF (probtype .EQ. 2) THEN
C       Lower Confidence Limit
c       LogLikelihood >= target
	   IF (I .EQ. 1) THEN
              CALL hillmean(X)
              CALL NegLogLike(X, LK)
C NegLogLike returns (in LK) the negative log likelihood.
              GXI = - LK - target
C	      PRINT *, 'BMDL, LL, target, and delta LL: ', X(1), -LK, target, GXI
	   ELSEIF (I .EQ. 2) THEN
C       BMDL > 0
	      GXI = X(GUNIT(2, I+NH)) - EPS
	   ELSEIF (I .EQ. 3) THEN
C       BMDL < BMD
	      GXI = bmd-X(1)
	      
	   ELSEIF (I .EQ. 4 .AND. constvar .EQ. 1) THEN
C       alpha > 0
	      GXI = X(GUNIT(2, I+NH)) - EPS

	   ELSEIF(I .EQ. 4 + constvar) THEN
C       n >= 0 or 1
              IF (restrict .EQ. 1) THEN
                 GXI = X(GUNIT(2, I+NH)) - 1
	      ELSE
                 GXI = X(GUNIT(2, I+NH)) 
              ENDIF
	   ELSEIF (I .EQ. 5 + constvar) THEN
C       K > 0
	      GXI = X(GUNIT(2, I+NH)) - EPS

	   ELSEIF(I .EQ. 6 + constvar) THEN
C       n <= 18
	      GXI = 18 - X(GUNIT(2,I+NH))

           ELSEIF (I .EQ. 7 + constvar) THEN
C       rho <= 18
              GXI = 18 - X(GUNIT(2, I + NH))
	   ENDIF
	   
	ELSEIF (probtype .EQ. 3) THEN
	  IF (constvar .EQ. 1 .AND. I .EQ. 1) THEN
	        GXI = X(GUNIT(2, I+NH)) - EPS
	  ELSEIF (I .EQ. 1 + constvar) THEN
                GXI = 18 - X(GUNIT(2, I + NH))
          ENDIF
        ENDIF
	
	RETURN
	END

	SUBROUTINE HillCompIneqGrad(X,GRAD,I)
	INCLUDE 'O8FUCO.INC'
	INCLUDE 'PROBLEM.INC'
	INTEGER I, J
	DOUBLE PRECISION X(*), GRAD(*)
        IF (probtype .EQ. 2) THEN
           IF (I .EQ. 1) THEN
	      CALL hillmean(X)
	      CALL hillpart(X)
              CALL DNegLogLike(X, GRAD)
C We are working with -NegLogLike, so need to change signs of gradient
              DO J=1, NX
                 GRAD(J) = -1.0D0*GRAD(J)
              ENDDO
           ENDIF
        ENDIF
	
	RETURN
	END
