c     Current limits: 100 dose groups, 25 total parameters.  
c     parmfixd(i) is 1 if the ith coefficient is fixed, 0 otherwise 
c     parmval(i) is the value to which the ith coefficient is fixed
c     risktype is	1 for absolute risk 
c			2 for std. dev. risk
c			3 for relative risk
c			4 for point risk 
c     restrict is 1 if coefficients are non-negative, -1 non-positive
c				  0 otherwise
c     probtype is 1 for ML estimation
c                 2 for Lower confidence limit
c                 3 for upper confidence limit
c		  4 for bmd ML profile

      INTEGER MAXDOSES, MAXORDER, MAXTUMOR
      PARAMETER (MAXDOSES = 100, MAXORDER = 23, MAXTUMOR = 10)
      DOUBLE PRECISION doses(MAXDOSES), bmr, bmd, target,
     $      parmval(0:MAXORDER), maxdose, lminbmd, lmaxbmd
      DOUBLE PRECISION mdDoses(0:MAXTUMOR,1:MAXDOSES), mdParmVal(0:MAXTUMOR, 0:MAXORDER)
      INTEGER ndoses, affected(MAXDOSES), nanimals(MAXDOSES), polyord,
     $      risktype,nparms, parmfixd(0:MAXORDER), restrict, probtype
     $      ,counter
      INTEGER nT, nObs, mnAnim(0:MAXTUMOR,1:MAXDOSES), arnRestrict(0:MAXTUMOR),
     $      mnParmFixd(0:MAXTUMOR, 0:MAXORDER), mnAffected(0:MAXTUMOR,1:MAXDOSES),
     $      arnParms(0:MAXTUMOR), arnObs(0:MAXTUMOR), arnPoly(0:MAXTUMOR)
      COMMON /PROBDATA/ doses, bmr, bmd, target, parmval, 
     $     maxdose, lminbmd, lmaxbmd, ndoses,
     $     affected, nanimals, polyord, risktype, nparms, parmfixd, 
     $     restrict, probtype, counter 
      COMMON /PROBDATA3/mdDoses, mdParmVal,
     $     nT, nObs, mnAnim, arnRestrict, mnParmFixd,
     $     mnAffected, arnParms, arnObs, arnPoly