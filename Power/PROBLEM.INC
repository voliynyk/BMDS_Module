c     Current limits: 100 dose groups, 24 parameters.
c     Since this is for the power model there are only 5 parameters
c     at most.
c  
c     parmfixd(i) = 1 if the ith coefficient is fixed
c                   0 otherwise 
c     parmval(i): is the value to which the ith coefficient is fixed
c     risktype = 0 for absolute risk 
c		 1 for std. dev. risk
c		 2 for relative risk
c	         3 for point risk
c                4 for extra risk(Hill model only)   
c     restrict = 1 if coefficients are non-negative
c                -1 non-positive
c		 0 otherwise
c     probtype = 1 for ML estimation
c                2 for Lower confidence limit
c                3 for model A3 (see C code for details) estimation
c                4 for Profile plot
c     DoseMean: is the mean at the estimated dose level in confidence 
c       limit calculation
c     DoseMeanGrad: is the vector of partial derivatives of DoseMean
c	with respect to all the parameters
c     ndoses: is the number of dose levels
c     nanimals: is the number of animals in each dose group
c     nparm: is the number of parameters
c            5 for getmle and getprofile
c            6 for getcl
c     adverse = 1 if the linear trend of the data is positive
c               -1 if the linear trend of the data is negative
c     modtype = 0 for polynomial model
c               1 for power model
c               2 for Hill model
c     constvar = 1 for constant variance model: alpha=sigma(i)^2
c                0 for modeled variance: exp(lalpha)*mu(dose(i))^rho=sigma(i)^2
c     flag: determines whether the parameters are scaled by their starting
c       values in donlp2
c     mlecnt: counts the number of calls to getmle
c     clcnt: counts the number of calls to getcl
c     a3cnt: counts the number of calls to getmlea3
c     procnt: counts the number of calls to getprofile
c     bmr: benchmark risk factor
c     bmd: benchmark dose
c     target: likelihood target value for the bmdl
c     mean: the user supplied mean at each dose level
c     means: the estimated mean at each dose level
c     grads: the gradients of the mean function at each dose level
c       with respect to each parameter
c     doses: the dose levels scaled by the max dose level (0<=doses<=1)
c     maxdose: the maximum dose level (before scaling)
c     mindose: the minimum dose level (before scaling)
c     bmdmean: the estimated mean at the bmd
c     bmdmeangrad: the gradients of the mean function at the bmd 
c	with respect to each parameter

      INTEGER MAXDOSES, MAXORDER	
      PARAMETER (MAXDOSES = 100, MAXORDER = 23)
      INTEGER ndoses,nanimals(MAXDOSES),nparm,
     $      risktype,parmfixd(0:MAXORDER),
     $      restrict,probtype,adverse,modtype,constvar,flag,
     $      mlecnt,clcnt,a3cnt,procnt		
      DOUBLE PRECISION bmr,bmd,target,parmval(0:MAXORDER),
     $      mean(MAXDOSES),var(MAXDOSES),means(MAXDOSES),
     $	    grads(MAXDOSES,1:MAXORDER),DoseMean,
     $	    DoseMeanGrad(1:MAXORDER), doses(MAXDOSES),
     $	    maxdose, mindose,bmdmean,bmdmeangrad(1:MAXORDER),
     $      maxmean, maxvar

      COMMON /PROBDATA/bmr,bmd,target,parmval,mean,var,
     $      means,grads,DoseMean,DoseMeanGrad,doses,maxdose,
     $      mindose,bmdmean,bmdmeangrad,maxmean,maxvar,
     $      probtype,ndoses,nanimals,nparm,risktype,parmfixd,
     $      restrict,modtype,constvar,adverse,flag,mlecnt,
     $      clcnt,a3cnt,procnt      
