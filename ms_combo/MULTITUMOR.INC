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
c     Added (01/23/2009) by Geoffrey to facilitate multiple (more than 2) tumors processing
c         MNMAXTUMOR - maximum number of tumors, default to 10
c         MNDOSES - maximum number of observations/doses, default to 100
c         MNMAXORDER - maximum of parameters, default to 23
c         nT - actual number of tumors
c         nMaxObs - actual maximum number of observations
c         nMaxNParms - actual maximum number of parameters
c         mdDoses(0:MNMAXTUMOR,1:MNDOSES) - double matrix for doses/tumor, index 0,1
c         mdParmVal(0:MNMAXTUMOR,0:MNMAXORDER) - double matrix for parmval/tumor, index 0
c         mnAnim(0:MNMAXTUMOR,1:MNDOSES) - int matrix of MNMAXTUMOR X MNDOSES for the animals, index 0,1
c         mnParmFixd(0:MNMAXTUMOR,0:MNMAXORDER) - int matrix of fixed parameter, index 0
c         mnAffected(0:MNMAXTUMOR,1:MNDOSES) - int matrix of affected values, index 0,1
c         arnParms(0:MNMAXTUMOR) - an array of integer for actual nParms/tumor, index 0
c         arnObs(0:MNMAXTUMOR) - an array of integer for actual nObs/tumor, index 0
c         arnPoly(0:MNMAXTUMOR) - an array of integer for number of poly/tumor, index 0
c         arnRestrict(0:MNMAXTUMOR) - integer restriction array, index 0

c      INTEGER MNDOSES, MNMAXORDER, MNMAXTUMOR
c      PARAMETER (MNDOSES = 100, MNMAXORDER = 23, MNMAXTUMOR = 10)
      DOUBLE PRECISION mdDoses(0:9,1:100),mdParmVal(0:9,0:22)
      INTEGER mnAnim(0:9,1:100),mnParmFixd(0:9,0:22),
     $      mnAffected(0:9,1:100),arnRestrict(0:9),
     $      arnParms(0:9),arnObs(0:9),arnPoly(0:9),nT, nMaxObs
      COMMON /MULTIDATA/ mdDoses, mdParmVal, mnAnim, 
     $     mnParmFixd, mnAffected, arnRestrict, arnParms, arnObs, 
     $     arnPoly, nT, nMaxObs

C     Louis Olszyk 07/2010 - Import debug log level from C code
C     geLogLevel specifies the amount of logging
      INTEGER geLogLevel, LOG_NONE, LOG_DEVELOPER
      COMMON /LOGINFO/ geLogLevel
C     Log level values
      PARAMETER (LOG_NONE = 0, LOG_DEVELOPER = 32768)
