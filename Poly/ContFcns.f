c	This file contains various functions that are common
c	to all of the Continuous models.  
c       
c	Calculates the log likelihood at the current
c	estimated vector of parameters	
c	The model mean() function (i.e., Hillmean, PolyMean, etc) should
c	have been called prior to invoking this routine
	
	SUBROUTINE NegLogLike(X,LK)
	INCLUDE 'PROBLEM.INC'
	INTEGER I
	DOUBLE PRECISION X(*),Vi(MAXDOSES),Ni,LK,SUM,DLOG
c       
	SUM= 0
	CALL EstVar(X,Vi)
	DO I = 1,ndoses
	   IF(Vi(I).LE.0) Vi(I) = 0.00000001
	   Ni = nanimals(I)
	   SUM = SUM + Ni*DLOG(Vi(I))/2 + (Ni-1)*var(I)/(2*Vi(I))
	   SUM = SUM + Ni*((mean(I)-means(I))**2)/(2*Vi(I))
	ENDDO
	LK = SUM
c       
	RETURN
	END
c	
c       
c	Returns the estimated variances at each dose.  The model
c	mean() function (i.e., Hillmean, PolyMean, etc)
c	should have been called prior to this function call
c       
	SUBROUTINE EstVar(X,Vi)
	INCLUDE 'PROBLEM.INC'
	INTEGER I,ptype
	DOUBLE PRECISION Vi(*),X(*)
c       
	ptype = probtype - 1
	IF(probtype.EQ.3) ptype = 0	
	IF(probtype.EQ.4) ptype = 0	
	DO I = 1,ndoses
	   IF(constvar.EQ.1) THEN
	      Vi(I) = X(1+ptype)
	   ELSEIF (ABS(means(I)) .GT. 1.0D-10) THEN
	      Vi(I) = EXP(X(1+ptype) + LOG((ABS(means(I)))) * X(2+ptype))
	   ELSE 
	      Vi(I) = 0.0
	   ENDIF
	ENDDO
	RETURN
	END
c       
c	calculates the partials of the -LogLikelihood function
c	with respect to all parameters.  Upon return to the 
c	calling function, G(I) = d(-L)/dX(I)
c	The model mean() function (i.e., Hillmean, PolyMean, etc) 
c	and the appropriate partial derivative function should
c	have been called prior to invoking this routine
	SUBROUTINE DNegLogLike(X,G)
	INCLUDE 'PROBLEM.INC'
	INTEGER I,ptype
	DOUBLE PRECISION Vi(MAXDOSES),dFi1,dFi2,dFi3
	DOUBLE PRECISION DVi(MAXDOSES,1:MAXORDER),DEVi,G(*),Ni,TEMP,X(*)
c
c	Get estimated variance
	DO I = 1,nparm
	   G(I) = 0.0
	ENDDO
	ptype = probtype - 1
	IF(probtype.EQ.3) ptype = 0
	IF(probtype.EQ.4) ptype = 0
	ChangeV = 0
	CALL EstVar(X,Vi)
	CALL VarPart(X,Vi,DVi)
	DO I = 1,ndoses
	   IF(Vi(I).LE.0) THEN
	      Vi(I) = .00000001
	   ENDIF
	   DEVi = mean(I) - means(I)
	   Ni = nanimals(I)
	   DO J = 1+ptype,nparm
	      dFi1 = DVi(I,J)/Vi(I)
	      dFi2 = -DVi(I,J)/(Vi(I)**2)
	      dFi3 = -(2*DEVi*Vi(I)*grads(I,J)+DEVi*DEVi*DVi(I,J))
	      dFi3 = dFi3/(Vi(I)*Vi(I))
	      TEMP = Ni*dFi1/2 + (Ni-1)*var(I)*dFi2/2 + Ni*dFi3/2
	      G(J) = G(J) + TEMP
	   ENDDO
	   IF(ptype.EQ.1) THEN
	      G(1) = 0.0
	   ENDIF		
	ENDDO
	RETURN
	END
c
c	Returns the partials of the variance function at each dose.  
c	The model mean() function (i.e., Hillmean, PolyMean, etc) should
c	have been called prior to invoking this routine
c	Upon return to calling function, DVi(I,J) = dV/dX(J) at dose I
	SUBROUTINE VarPart(X,Vi,DVi)
	INCLUDE 'PROBLEM.INC'
	INTEGER I,ptype,sign
	DOUBLE PRECISION Vi(*),DVi(MAXDOSES,*),X(*)
c       
	ptype = probtype - 1
	IF(probtype.EQ.3) ptype = 0
	IF(probtype.EQ.4) ptype = 0	
	sign = 1
	DO I = 1,ndoses
c	   IF(means(I).LT.0) sign = -1
c       Constant Variance
	   IF(constvar.EQ.1) THEN
	      DVi(I,1+ptype) = 1.0
	      DO J = 2+ptype, nparm
		 DVi(I,J+ptype) = 0.0
	      ENDDO
c       Modeled Variance
	   ELSE
	      DVi(I,1+ptype) = Vi(I)
	      IF (ABS(means(I)) .GT. 1.0D-10) THEN
		 DVi(I,2+ptype) = Vi(I)*DLOG(ABS(means(I)))
	      ELSE
		 DVi(I,2+ptype) = 0.0
	      ENDIF
	      IF (parmfixd(1) .EQ. 1) DVi(I,2+ptype) = 0.0
c       
	      DO J = 3+ptype,nparm
		 DVi(I,J) = sign*X(2+ptype)*Vi(I)*grads(I,J)
		 IF (means(I) .EQ. 0) THEN
		    DVi(I,J) = 0.0
		 ELSE 
		    DVi(I,J) = DVi(I,J)/ABS(means(I))
		 ENDIF  
	      ENDDO
c       
	      
	      IF(ptype.EQ.1)  DVi(I,1) = 0.0
	   ENDIF
c       
	ENDDO	
	RETURN
	END

c	Computes the value of the bmr function.
c	 0 - absolute deviation
C	 1 - Std Dev
c	 2 - Relative
c	 3 - Point
c	 4 - extra (Hill only)
c
c	BMRVal is returned with the current iterations BMR value
c	i.e.,
c	if 0, then BMRVal = ABS(Mean(Dose) - Mean(0))
c	if 1, then BMRVal = ABS(Mean(Dose) - Mean(0))/Sqrt(Var(0))
c	if 2, then BMRVal = ABS(Mean(Dose) - Mean(0))/Mean(0)
c	if 3, then BMRVal = Mean(Dose)
c	if 4, then BMRVal = (Mean(Dose) - Mean(0))/(Mean(max) - Mean(0))
c
c	The model mean routine should be called prior to invoking this
c	routine
c
	SUBROUTINE BMRComp(X,bmrtype,bmrval)
	INCLUDE 'PROBLEM.INC'
	INTEGER bmrtype 
	DOUBLE PRECISION V,STD,X(*),bmrval, DLOG

c       This is the bmr constraint for Profile Calc.
	IF(probtype .EQ. 4) THEN
	   IF(bmrtype .EQ. 0) THEN
	      bmrval = DABS(bmdmean - X(3))
	   ELSEIF(bmrtype.EQ.1) THEN
	      IF(constvar.EQ.1) THEN
		 V = X(1)
	      ELSE
		 V = EXP(X(1) + DLOG(DABS(X(3)))*X(2))
	      ENDIF
	      IF(V.LE.0) THEN
		 V = .00000001
	      ENDIF
	      STD = DSQRT(V)
	      bmrval = DABS(bmdmean - X(3))/STD
	   ELSEIF(bmrtype.EQ.2) THEN
	      bmrval = DABS(bmdmean - X(3))/X(3)
	   ELSEIF(bmrtype.EQ.3) THEN
	      bmrval = bmdmean
	   ELSE
	      bmrval = (bmdmean - X(3))/X(4)
	   ENDIF
	ELSEIF(probtype .EQ. 2) THEN
c       This is the bmr constraint for lower confidence limit
	   IF(bmrtype.EQ.0) THEN
	      bmrval = DABS(DoseMean - X(4))
	   ELSEIF(bmrtype.EQ.1) THEN
	      IF(constvar.EQ.1) THEN
		 V = X(2)
	      ELSE
		 V = EXP(X(2) + DLOG(DABS(X(4)))*X(3))
	      ENDIF
	      IF(V.LE.0) THEN
		 V = .00000001
	      ENDIF
	      STD = DSQRT(V)
	      bmrval = DABS(DoseMean - X(4))/STD
	   ELSEIF(bmrtype.EQ.2) THEN
	      bmrval = DABS(DoseMean - X(4))/X(4)
	   ELSEIF(bmrtype.EQ.3) THEN
	      bmrval = DoseMean
	   ELSE
	      bmrval = (DoseMean - X(4))/X(5)
	   ENDIF
	ENDIF
	RETURN
	END
c
c	DBMRComp computes the partial derivatives of the BMR
c	equality constraint in the Confidence Limit problem.
c	Upon return, G(J) = d(BMDfunc)/dX(J)
c	The model mean routine and partial derivative routine
c	should be called prior to invoking this
c	routine.  NOTE:  This needs to change by model!!!!
c
	SUBROUTINE DBMRComp(X,bmrtype,G)
	INCLUDE 'PROBLEM.INC'
	INTEGER sign, J, ChangeV,bmrtype,sign2
	DOUBLE PRECISION V,STD,DABS,DLOG,MeanDev,G(*),X(*)
c	OPEN(11,FILE='BMRConst.txt',STATUS='UNKNOWN')
c
c       Profile Calc.
	IF (probtype .EQ. 4) THEN
	   MeanDev = bmdmean - X(3)
	   IF(MeanDev.LT.0) THEN
	      sign = -1
	   ELSE
	      sign = 1
	   ENDIF
	   IF(bmrtype.EQ.0) THEN
	      G(1) = 0.0
	      G(2) = 0.0
	      G(3) = 0.0
	      DO J = 4, nparm
		 G(J) = sign*bmdmeangrad(J)
	      ENDDO
	   ELSEIF(bmrtype.EQ.1) THEN
	      IF(X(3).LT.0) THEN
		 sign2 = -1
	      ELSE
		 sign2 = 1
	      ENDIF
	      IF(constvar.EQ.1) THEN
		 V = X(1)
	      ELSE
		 V = EXP(X(1) + LOG(DABS(X(3))) * X(2))
	      ENDIF
	      IF(V.LE.0) THEN
		 V = .00000001
		 ChangeV = 1
	      ENDIF
	      STD = DSQRT(V)
	      IF(constvar.EQ.0) THEN
		 IF(ChangeV.EQ.1) THEN
		    G(1) = 0.0
		    G(2) = 0.0
		 ELSE
		    G(1) = -sign*MeanDev/(2*STD)
		    G(2) = -sign*MeanDev*DLOG(DABS(X(3)))/(2*STD)
		 ENDIF
		 G(3) = -sign*sign2*X(2)*MeanDev
		 G(3) =	G(3)/(2*DABS(X(3))*STD)
	      ELSE
		 G(2) = 0.0
		 IF(ChangeV.EQ.1) THEN
		    G(1) = 0.0
		 ELSE
		    G(1) = -sign*MeanDev/(2*V*STD)
		 ENDIF
		 G(3) = 0
	      ENDIF
	      DO J = 4, nparm
		 G(J) = sign*bmdmeangrad(J)/STD
	      ENDDO
	      
	   ELSEIF(bmrtype .EQ. 2) THEN
	      G(1) = 0.0
	      G(2) = 0.0
	      G(3) = -sign*MeanDev/(X(3)*X(3))
	      DO J = 4, nparm
		 G(J) = sign*bmdmeangrad(J)/X(3)
	      ENDDO
	   ELSEIF(bmrtype .EQ. 3) THEN
	      G(1) = 0.0
	      G(2) = 0.0
	      G(3) = bmdmeangrad(4)
	      DO J = 4, nparm
		 G(J) = bmdmeangrad(J)
	      ENDDO
	   ELSE
	   ENDIF

c       Lower Confidence Limit
	ELSEIF (probtype .EQ. 2) THEN
	   MeanDev = DoseMean - X(4)
	   IF(MeanDev.LT.0) THEN
	      sign = -1
	   ELSE
	      sign = 1
	   ENDIF
	   IF(bmrtype.EQ.0) THEN
	      G(1) = sign*DoseMeanGrad(1)
	      G(2) = 0.0
	      G(3) = 0.0
	      G(4) = 0.0
	      DO J = 5, nparm
		 G(J) = sign*DoseMeanGrad(J)
	      ENDDO
	   ELSEIF(bmrtype.EQ.1) THEN
	      IF(X(4).LT.0) THEN
		 sign2 = -1
	      ELSE
		 sign2 = 1
	      ENDIF
	      IF(constvar.EQ.1) THEN
		 V = X(2)
	      ELSE
		 V = EXP(X(2) + LOG(DABS(X(4)))*X(3))
	      ENDIF
	      IF(V.LE.0) THEN
		 V = .00000001
		 ChangeV = 1
	      ENDIF
	      STD = DSQRT(V)
	      IF(constvar.EQ.0) THEN
		 IF(ChangeV.EQ.1) THEN
		    G(2) = 0.0
		    G(3) = 0.0
		 ELSE
		    G(2) = -sign*MeanDev/(2*STD)
		    G(3) = -sign*MeanDev*DLOG(DABS(X(4)))/(2*STD)
		 ENDIF
		 G(4) = -sign*sign2*X(3)*MeanDev
		 G(4) =	G(4)/(2*DABS(X(4))*STD)
	      ELSE
		 G(3) = 0.0
		 IF(ChangeV.EQ.1) THEN
		    G(2) = 0.0
		 ELSE
		    G(2) = -sign*MeanDev/(2*V*STD)
		 ENDIF
		 G(4) = 0
	      ENDIF
	      G(1) = sign*DoseMeanGrad(1)/STD
	      DO J = 5, nparm
		 G(J) = sign*DoseMeanGrad(J)/STD
	      ENDDO
	   ELSEIF(bmrtype .EQ. 2) THEN
	      G(1) = sign*DoseMeanGrad(1)/X(4)
	      G(2) = 0.0
	      G(3) = 0.0
	      G(4) = -sign*MeanDev/(X(4)*X(4))
	      DO J = 5, nparm
		 G(J) = sign*DoseMeanGrad(J)/X(4)
	      ENDDO
	   ELSEIF(bmrtype .EQ. 3) THEN
	      G(1) = DoseMeanGrad(1)
	      G(2) = 0.0
	      G(3) = 0.0
	      G(4) = DoseMeanGrad(4)
	      DO J = 5, nparm
		 G(J) = DoseMeanGrad(J)
	      ENDDO
	   ELSE
	   ENDIF
	ENDIF
	RETURN
	END
