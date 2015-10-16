c       This file contains various functions that are common
c       to all of the Continuous models.  
c       
c       Calculates the log likelihood at the current
c       estimated vector of parameters  
c       The model mean() function (i.e., Hillmean, PolyMean, etc) should
c       have been called prior to invoking this routine
        
c       DOUBLE PRECISION FUNCTION powerf(a,b)
c       INCLUDE 'PROBLEM.INC'
c       DOUBLE PRECISION loga, bloga, a, b
        
c       if (a .GE. .0000000001) then
c          loga = log(a)
c       else
c          loga = -24.5258509299404572 + (2.0D10)*a - (5.0D19)*a*a
c       endif
c       bloga = b*loga
c       if (bloga .GT. 700) then
c           bloga = 700
c       endif
c       powerf = dexp(bloga)
c       RETURN
c       END

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
c       Returns the estimated variances at each dose.  The model
c       mean() function (i.e., Hillmean, PolyMean, etc)
c       should have been called prior to this function call
c       
        SUBROUTINE EstVar(X,Vi)
        INCLUDE 'PROBLEM.INC'
        INTEGER I,ptype
        DOUBLE PRECISION Vi(*),X(*)
c       
        IF (probtype .EQ. 1 .OR. probtype .EQ. 3 .OR.
     1    probtyp .EQ. 4) THEN
           ptype = 0
        ELSE
           ptype = 1
        ENDIF
        DO I = 1,ndoses
           IF (ABS(means(I)) .GT. 1.0D-10) THEN
              Vi(I) = EXP(X(1+ptype) + LOG(ABS(means(I))) * X(2+ptype))
           ELSE 
              Vi(I) = 0.0
           ENDIF
        ENDDO

        RETURN
        END
c       
c       calculates the partials of the -LogLikelihood function
c       with respect to all parameters.  Upon return to the 
c       calling function, G(I) = d(-L)/dX(I)
c       The model mean() function (i.e., Hillmean, PolyMean, etc) 
c       and the appropriate partial derivative function should
c       have been called prior to invoking this routine
        SUBROUTINE DNegLogLike(X,G)
        INCLUDE 'PROBLEM.INC'
        INTEGER I,ptype
        DOUBLE PRECISION Vi(MAXDOSES),dFi1,dFi2,dFi3
        DOUBLE PRECISION DVi(MAXDOSES,1:MAXORDER),DEVi,G(*),Ni,TEMP,X(*)
c
c       Get estimated variance
        DO I = 1,nparm
           G(I) = 0.0
        ENDDO
        ptype = probtype - 1
        IF(probtype.EQ.3) ptype = 0
        IF(probtype.EQ.4) ptype = 0     
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
        DO J = 1+ptype,nparm
              IF(parmfixd(J-1-ptype).EQ.1) THEN
                 G(J) = 0.0
          ENDIF
                ENDDO
        RETURN
        END
c
c       Returns the partials of the variance function at each dose.  
c       The model mean() function (i.e., Hillmean, PolyMean, etc) should
c       have been called prior to invoking this routine
c       Upon return to calling function, DVi(I,J) = dV/dX(J) at dose I
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
                IF(means(I).LT.0)       sign = -1
c       
c       Modeled Variance
      DVi(I,1+ptype) = Vi(I)
      IF ((parmfixd(1) .EQ. 1).AND.(X(2+ptype).EQ.0.0)) THEN
                DO J = 2+ptype, nparm
          DVi(I,J) = 0.0
        ENDDO
      ELSE
        IF (parmfixd(1) .EQ. 1) THEN
               DVi(I,2+ptype) = 0.0
        ELSE
               IF (ABS(means(I)) .GT. 1.0D-10) THEN
                      DVi(I,2+ptype) = Vi(I)*DLOG(ABS(means(I)))
               ELSE
                      DVi(I,2+ptype) = 0.0
               ENDIF
        ENDIF
      DO J = 3+ptype,nparm
          DVi(I,J) = sign*X(2+ptype)*Vi(I)*grads(I,J)
                IF (means(I) .EQ. 0) THEN
                        DVi(I,J) = 0.0
                ELSE 
                        DVi(I,J) = DVi(I,J)/ABS(means(I))
                ENDIF
      ENDDO
                ENDIF
                IF(ptype.EQ.1)  DVi(I,1) = 0.0
          ENDDO
        RETURN
        END

        
c               DO I = 1,ndoses
c                       IF(means(I).LT.0)       sign = -1
cc       Constant Variance
cc           IF(constvar.EQ.1) THEN
cc              DVi(I,1+ptype) = 1.0
cc              DO J = 2+ptype, nparm
cc                 DVi(I,J+ptype) = 0.0
cc              ENDDO
cc       Modeled Variance
c                       IF(constvar.NE.1) THEN
c                               DVi(I,1+ptype) = Vi(I)
c                               IF (ABS(means(I)) .GT. 1.0D-10) THEN
c                                       DVi(I,2+ptype) = Vi(I)*DLOG(ABS(means(I)))
c                               ELSE
c                                       DVi(I,2+ptype) = 0.0
c                               ENDIF
c                       ENDIF
c                       IF (parmfixd(1) .EQ. 1) THEN
c                               DVi(I,2+ptype) = 0.0
c                               DO J = 3+ptype,nparm
c                                       DVi(I,J) = sign*X(2+ptype)*Vi(I)*grads(I,J)
c                                       IF (means(I) .EQ. 0) THEN
c                                               DVi(I,J) = 0.0
c                                       ELSE 
c                                               DVi(I,J) = DVi(I,J)/ABS(means(I))
c                                       ENDIF  
c                               ENDDO
c                               IF(ptype.EQ.1)  DVi(I,1) = 0.0
c                       ENDIF
c               ENDDO   
c        RETURN
c        END

c       Computes the value of the bmr function.
c        0 - absolute deviation
C        1 - Std Dev
c        2 - Relative
c        3 - Point
c        4 - extra (Hill only)
c
c       BMRVal is returned with the current iterations BMR value
c       i.e.,
c       if 0, then BMRVal = ABS(Mean(Dose) - Mean(0))
c       if 1, then BMRVal = ABS(Mean(Dose) - Mean(0))/Sqrt(Var(0))
c       if 2, then BMRVal = ABS(Mean(Dose) - Mean(0))/Mean(0)
c       if 3, then BMRVal = Mean(Dose)
c       if 4, then BMRVal = (Mean(Dose) - Mean(0))/(Mean(max) - Mean(0))
c
c       The model mean routine should be called prior to invoking this
c       routine
c
        SUBROUTINE BMRComp(X,bmrtype,bmrval)
        INCLUDE 'PROBLEM.INC'
        INTEGER bmrtype 
        DOUBLE PRECISION V,STD,X(*),bmrval

c       This is the bmr constraint for Profile Calc.
        IF(probtype .EQ. 4) THEN
           IF(bmrtype .EQ. 0) THEN
              bmrval = DABS(bmdmean - X(3))
           ELSEIF(bmrtype.EQ.1) THEN
              V = EXP(X(1) + LOG(DABS(X(3))) * X(2))
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

c       This is the bmr constraint for lower confidence limit
        ELSEIF(bmrtype.EQ.0) THEN
           bmrval = DABS(DoseMean - X(4))
           IF(lognorm.EQ.1) bmrval = DABS(EXP(DoseMean) - X(4))
        ELSEIF(bmrtype.EQ.1) THEN
           
              V = EXP(X(2) + LOG(DABS(X(4))) * X(3))
           
           IF(V.LE.0) THEN
              V = .00000001
           ENDIF
           STD = DSQRT(V)
           bmrval = DABS(DoseMean - X(4))/STD
           IF(lognorm.EQ.1) bmrval = DABS((DoseMean - LOG(X(4))))/STD
        ELSEIF(bmrtype.EQ.2) THEN
           bmrval = DABS(DoseMean - X(4))/X(4)
           IF(lognorm.EQ.1) bmrval = DABS(EXP(DoseMean) - X(4))/X(4)
        ELSEIF(bmrtype.EQ.3) THEN
           bmrval = DoseMean
           IF(lognorm.EQ.1) bmrval = EXP(DoseMean)
        ELSE
           bmrval = DABS(DoseMean - X(4))/DABS(X(4) - X(4)*X(6))
           IF(lognorm.EQ.1) THEN
              bmrval=DABS((EXP(DoseMean))-X(4))/DABS(X(4)-X(4)*X(6))
           ENDIF
        ENDIF
        
        RETURN
        END
c
c       DBMRComp computes the partial derivatives of the BMR
c       equality constraint in the Confidence Limit problem.
c       Upon return, G(J) = d(BMDfunc)/dX(J)
c       The model mean routine and partial derivative routine
c       should be called prior to invoking this
c       routine.  NOTE:  This needs to change by model!!!!
c
        SUBROUTINE DBMRComp(X,bmrtype,G)
        INCLUDE 'PROBLEM.INC'
        INTEGER sign, J, ChangeV,bmrtype,sign2
        DOUBLE PRECISION V,STD,DABS,DLOG,MeanDev,G(*),X(*)
c       OPEN(11,FILE='BMRConst.txt',STATUS='UNKNOWN')
c
        ChangeV = 0
c       Profile Calc.
C        IF (probtype .EQ. 4) THEN  NOT DONE


c       Lower Confidence Limit
        IF (probtype.EQ.2) THEN
           MeanDev = DoseMean - X(4)
           IF(lognorm.EQ.1) MeanDev = EXP(DoseMean)-X(4)
           IF(MeanDev.LT.0) THEN
              sign = -1
           ELSE
              sign = 1
           ENDIF
           IF(bmrtype.EQ.0) THEN
              G(1) = sign*DoseMeanGrad(1)
              IF(lognorm.EQ.1) G(1) = sign*DoseMeanGrad(1)*EXP(DoseMean)
              G(2) = 0.0
              G(3) = 0.0
              G(4) = sign*(DoseMeanGrad(4) - 1)
              IF(lognorm.EQ.1) THEN
                 G(4) = sign*(((EXP(DoseMean))*DoseMeanGrad(4))-1)
              ENDIF
              DO J = 5, nparm
                 G(J) = sign*DoseMeanGrad(J)
                 IF(lognorm.EQ.1) THEN
                    G(J) = sign*DoseMeanGrad(J)*EXP(DoseMean)
                 ENDIF
              ENDDO
           ELSEIF(bmrtype.EQ.1) THEN
                 V = EXP(X(2) + LOG(DABS(X(4))) * X(3))
              
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
c                 G(4) = -sign*X(3)*MeanDev
c                 G(4) = G(4)/(2*DABS(X(4))*STD)
                  G(4) = DoseMeanGrad(4) - 1.0
                  G(4) = G(4) - MeanDev*X(3)/(2.0 * X(4))
                  G(4) = (sign/STD)* G(4)
              ELSE
                 G(3) = 0.0
                 IF(ChangeV.EQ.1) THEN
                    G(2) = 0.0
                 ELSE
                    G(2) = -sign*MeanDev/(2*STD)
                    IF(lognorm.EQ.1) THEN
                       G(2)=-sign*(DoseMean-LOG(X(4)))/(2*STD)
                    ENDIF
                 ENDIF
                 G(4) = sign*(DoseMeanGrad(4) - 1)/STD
                 IF(lognorm.EQ.1)
     $                G(4) = sign*(DoseMeanGrad(4)-1.0/X(4))/STD
              ENDIF
              G(1) = sign*DoseMeanGrad(1)/STD
              DO J = 5, nparm
                 G(J) = sign*DoseMeanGrad(J)/STD
              ENDDO
           ELSEIF(bmrtype .EQ. 2) THEN
              G(1) = sign*DoseMeanGrad(1)/X(4)
              IF(lognorm.EQ.1)
     $             G(1) = sign*EXP(DoseMean)*DoseMeanGrad(1)/X(4)
              G(2) = 0.0
              G(3) = 0.0
      G(4) = sign*(((DoseMeanGrad(4)-1)/X(4))-(MeanDev/(X(4)*X(4))))
      IF(lognorm.EQ.1) THEN
c     Following code was commented and divided into segment, coz of compiler error 
c     G(4) = sign*(((DoseMeanGrad(4)*EXP(DoseMean)-1)/X(4))-(MeanDev/(X(4)*X(4))))      
      G(4) = DoseMeanGrad(4)/X(4) - 1.0/(X(4)*X(4))
      G(4) = sign * EXP(DoseMean) * G(4)
      ENDIF
              DO J = 5, nparm
                 G(J) = sign*DoseMeanGrad(J)/X(4)
                 IF(lognorm.EQ.1) G(J) = sign*(EXP(DoseMean))*
     $                DoseMeanGrad(J)/X(4)
              ENDDO
           ELSEIF(bmrtype .EQ. 3) THEN
              G(1) = DoseMeanGrad(1)
              IF(lognorm.EQ.1) G(1) = DoseMeanGrad(1)*EXP(DoseMean)
              G(2) = 0.0
              G(3) = 0.0
              G(4) = DoseMeanGrad(4)
              IF(lognorm.EQ.1) G(4) = DoseMeanGrad(4)*EXP(DoseMean)
              DO J = 5, nparm
                 G(J) = DoseMeanGrad(J)
                 IF(lognorm.EQ.1) G(J) = DoseMeanGrad(J)*EXP(DoseMean)
              ENDDO
           ELSE
              G(1) = DoseMeanGrad(1)/(X(6)*X(4) - X(4))
              IF(lognorm.EQ.1)
     $             G(1) = EXP(DoseMean)*DoseMeanGrad(1)/
     $             (X(6)*X(4) - X(4))
              G(4) = (DoseMeanGrad(4) - 1)/(X(6)*X(4) - X(4)) -
     $             MeanDev * ((X(6) - 1)/((X(6)*X(4) -
     $             X(4))*(X(6)*X(4) - X(4))))
              IF(lognorm.EQ.1) THEN
c     Following code was commented and divided into segment, coz of compiler error 
c      G(4) = (DoseMeanGrad(4)*EXP(DoseMean) - 1)/(X(6)*X(4) - X(4)) - MeanDev * ((X(6) - 1)/((X(6)*X(4) - X(4))*(X(6)*X(4) - X(4))))
                 G(4) = (DoseMeanGrad(4)*(EXP(DoseMean)) - 1)
                 G(4) = G(4)/(X(6)*X(4) - X(4))
                 G(2) = (X(6)*X(4) - X(4))
                 G(3) = G(2) * (X(6)*X(4) - X(4))
                 G(2) = MeanDev * (X(6) - 1)/G(3)
                 G(4) = G(4) - G(2)
              ENDIF
              G(2) = 0.0
              G(3) = 0.0
              G(5) = DoseMeanGrad(5)/(X(6)*X(4) - X(4))
              IF(lognorm.EQ.1) G(5) = DoseMeanGrad(5)*EXP(DoseMean)/
     $             (X(6)*X(4) - X(4))
              G(6) = DoseMeanGrad(6)/(X(6)*X(4) - X(4)) -
     $             MeanDev*X(4)/((X(6)*X(4) - X(4))*(X(6)*X(4) - X(4)))
              IF(lognorm.EQ.1) G(6) = DoseMeanGrad(6)*EXP(DoseMean)/
     $             (X(6)*X(4) - X(4)) - MeanDev*X(4)/
     $             ((X(6)*X(4) - X(4))*(X(6)*X(4) - X(4)))
              G(7) = DoseMeanGrad(7)/(X(6)*X(4) - X(4))
              IF(lognorm.EQ.1) G(7) = DoseMeanGrad(7)*EXP(DoseMean)/
     $             (X(6)*X(4) - X(4))
           ENDIF
        ENDIF
        RETURN
        END
