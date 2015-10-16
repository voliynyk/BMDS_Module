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
c                 3 for model A3 (see C code for details) estimation
c	  DoseMean is the mean at the estimated dose level in confidence 
c				limit calculation
c	  DoseMeanGrad is the vector of partial derivatives of DoseMean
c				with respect to all the parameters
c     mlecnt counts the number of calls to getmle_
c     clcnt counts the number of calls to getcl_
c     a3cnt counts the number of calls to getmlea3_
c
      INTEGER MAXDOSES, MAXORDER	
      PARAMETER (MAXDOSES = 100, MAXORDER = 23)
      INTEGER ndoses,nanimals(MAXDOSES),nparm
     $      ,risktype,nparms,parmfixd(0:MAXORDER)
     $      ,restrict,probtype,adverse,modtype,constvar
     $      ,flag,mlecnt,clcnt,a3cnt		
      DOUBLE PRECISION bmr,bmd,target,parmval(0:MAXORDER)
     $      ,mean(MAXDOSES),var(MAXDOSES),means(MAXDOSES)
     $	    ,grads(MAXDOSES,1:MAXORDER),DoseMean
     $	    ,DoseMeanGrad(1:MAXORDER), doses(MAXDOSES) 
      COMMON /PROBDATA/
     $      doses,bmr,bmd,target,parmval,
     $      means,grads,DoseMean,DoseMeanGrad,
     $      mean,var, 
     $      probtype,ndoses,nanimals,nparm,
     $      risktype,nparms,parmfixd,constvar,
     $      adverse,flag,mlecnt,clcnt,a3cnt,
     $      restrict,modtype
      
