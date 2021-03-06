C**** IF BLOC=.TRUE. THEN IT IS ASSUMED THAT FUNCTIONEVALUATION TAKES PLACE
C**** AT ONCE IN AN EXTERNAL ROUTINE AND 
C**** THAT USER_EVAL HAS BEEN CALLED BEFORE CALLING FOR EVALUATION OF FUNCTIONS
C**** THE LATTER THEN SIMPLY CONSISTS IN COPYING DATA FROM COMMON O8FUEXT
C**** TO DONLP2'S OWN DATA
C**** USER_EVAL MUST SET VALID=TRUE, IF FUNCTIONVALUES ARE VALID FOR THE 
C**** CURRENT XTR
C**** CORR IS SET TO TRUE BY DONLP2, IF THE INITIAL X DOES NOT SATISFY
C**** THE BOUND CONSTRAINTS. X IS MODIFIED IN THIS CASE
C**** DIFFTYPE=1,2,3 NUMERICAL DIFFERENTIATION BY THE ORDINARY FORWARD
C**** DIFFERENCES, BY CENTRAL DIFFERENCES OR BY RICHARDSON-EXTRAPOLATION
C**** OF ORDER 6, REQUIRING N, 2N , 6N ADDITIONAL FUNCTION EVALUATIONS
C**** RESPECTIVELY
C**** EPSFCN IS THE ASSUMED PRECISION OF THE FUNCTION EVALUATION, TO BE
C**** SET BY THE USER
C**** TAUBND: AMOUNT BY WHICH BOUND CONSTRAINTS MAY BE VIOLATED DURING
C**** FINITE DIFFERENCING, SET BY THE USER
      LOGICAL BLOC,VALID,CORR
      INTEGER DIFFTYPE
      DOUBLE PRECISION XTR(NX),XSC(NX),FU(0:NRESM),FUGRAD(NX,0:NRESM),
     F                 FUD(0:NRESM,1:6),EPSFCN,TAUBND
      COMMON/O8FUEXT/XTR,XSC,FU,FUGRAD,FUD,EPSFCN,TAUBND
      COMMON/O8FUPAR/BLOC,VALID,CORR,DIFFTYPE
