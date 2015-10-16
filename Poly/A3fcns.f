c       This file contains three functions that are specific to the
c       to the hill model when doing maximum likelihood and confidence
c       limit calculations using donlp2 and the associated user 
c       functions which are non-specific to the continuous model
c       that is being run.
c       FILLGUNIT fills the GUNIT array in donlp2 with the model 
c       specific equality, non-equality constraint as well as
c       the objective function.
c       A3MEAN just calculates the mean specific to that model as a
c       vector of means, one at each dose level
c       A3PART calculates the vector of partial derivates of the
c       mean with respect to all model parameters.
c       

	SUBROUTINE A3FILLGUNIT
	INCLUDE 'O8COMM.INC'
	INCLUDE 'PROBLEM.INC'
	INTEGER I,J
C       
C       NAME IS IDENT OF THE EXAMPLE/USER AND CAN BE SET AT USERS WILL
C       THE FIRST CHARACTER MUST BE ALPHABETIC.  40 CHARACTERS MAXIMUM
c       NAME='A3 LIKELIHOOD'
C       X IS INITIAL GUESS AND ALSO HOLDS THE CURRENT SOLUTION
C       PROBLEM DIMENSION N=DIM(X), NH=DIM(H), NG=DIM(G)
C       
C       ML Estimation
C       
C       Objective function is a function of all the parameters
	
	GUNIT(1,0) = -1
	GUNIT(2,0) = 0
	GUNIT(3,0) = 0
C       The only equality constraints in the A3 model are those
C       from presetting alpha and rho (X[1] and X[2])
	
	NH = 0
	J = 0
	DO I = 0, 1
	   IF (parmfixd(I) .EQ. 1) THEN
	      J = J + 1
	      NH = NH + 1
	      GUNIT(1, J) = 1
	      GUNIT(2, J) = I + 1
	      GUNIT(3, J) = 1
	   ENDIF
	ENDDO
	NG = 0
C       
C	Variance coefficient parameter, alpha, must be positive, when
C	using the constant variance model
C       
	IF (constvar .EQ. 1) THEN
	   J = J + 1
	   NG = NG + 1
	   GUNIT(1, J) = 1
	   GUNIT(2, J) = 1
	   GUNIT(3, J) = 1 
	ENDIF
C       Constrain rho <= 18
	J = J + 1
	NG = NG + 1
	GUNIT(1, J) = 1
	GUNIT(2, J) = 2
	GUNIT(3, J) = -1
	RETURN
	END
c       
	SUBROUTINE A3MEAN(X)
c	INCLUDE 'O8COMM.INC'
	INCLUDE 'PROBLEM.INC'
	INTEGER I
	DOUBLE PRECISION X(*)
c       
	DO I = 1, nparm
	   means(I) = X(I+2)
	ENDDO
c       
	RETURN
	END
c	
c       
c	A3part gives the partial derivatives of the A3 model mean function
c	and also, in grads(ndoses+1, j) gives the partial derivatives
c	of the mean function evaluated at the estimated dose D.
	SUBROUTINE A3PART(X)
c	INCLUDE 'O8COMM.INC'
	INCLUDE 'PROBLEM.INC'
	INTEGER I,J
	DOUBLE PRECISION X(*)
c       
	DO I = 1, ndoses
c	Maximum likelihood estimation A3 model mean partial derivatives
	   grads(I,1) = 0.0
	   grads(I,2) = 0.0
	   DO J = 3, nparm
	      IF(J.EQ.(I+2)) THEN
		 grads(I, J) = 1.0
	      ELSE
		 grads(I, J) = 0.0
	      ENDIF
	   ENDDO
	ENDDO
	RETURN
	END
