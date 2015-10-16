c     LOADCOMMBLOC : Loads values into the common block for
c     donlp2.
c     INPUT:
c         XNDOSES     number of dose levels (INTEGER)
c         XDOSES(*)   dose levels used      (DOUBLE)
c         XMEAN(*)    the sample mean at Ith dose (DOUBLE)
c         XNANIM(*)   number of animals at Ith dose (INTEGER)
c         XVAR(*)     sample variance at Ith dose(DOUBLE)
c         XNPARMS     the number of parameters in the model (INTEGER)
c         XFIXD(*)    1 if I'th (0:polyord) parameter is specified
c                     0 otherwise (INTEGER)
c         XVAL(*)     Value of fixed parameter (0:polyord) (DOUBLE)
c         XRESTR      1 if parameters are restricted to non-negative values
c                     0 if there is no restriction
c     	              -1 if parameters are restricted to non-positive values
c         XADVERSE    1 if adverse direction is larger
c		      -1 if adverse direction is smaller
c	  XMODEL      0 if Polynomial
c		      1 if Power
c		      2 if Hill
c         XMAX        the max dose level (DOUBLE)
c         XMIN        the min dose level (DOUBLE)
c
c
c      SUBROUTINE LOADCOMMBLOC(XNDOSES,XDOSES,XMEAN,XNANIM,XVAR,
c     $     XNPARMS,XFIXD,XVAL,XRESTR,XADVERSE,XMODEL,XMAX,XMIN)
      SUBROUTINE LOADCOMMBLOC(YMMAX,YDMAX)

      INCLUDE 'PROBLEM.INC'     
c      INTEGER XNDOSES,XNANIM(*),XNPARMS,XFIXD(0:XNPARMS-1),
c     $     XRESTR,XADVERSE,XMODEL,
c     $     I,temp3
c      DOUBLE PRECISION XDOSES(*),XMEAN(*),XVAR(*),XVAL(0:XNPARMS-1),
c     $     XMAX,XMIN,
c     $     temp1,temp2,temp4
      DOUBLE PRECISION YMMAX,YDMAX

      mlecnt = 0
      clcnt =0
      a3cnt =0

      maxmean = YMMAX
      maxvar = YDMAX
c      ndoses = XNDOSES
c      nparm = XNPARMS
c      restrict = XRESTR
c      adverse = XADVERSE
c      modtype = XMODEL
c      maxdose = XMAX
c      mindose = XMIN
c      DO I = 1, ndoses
c         doses(I) = XDOSES(I)
c         mean(I) = XMEAN(I)
c         nanimals(I) = XNANIM(I)
c         var(I) = XVAR(I)
c      ENDDO
c      DO I = 0, nparm-1
c         parmfixd(I) = XFIXD(I)
c         parmval(I) = XVAL(I)
c      ENDDO
c      IF(parmfixd(1).EQ.1.AND.parmval(1).EQ.0) THEN
c         constvar = 1
c      ELSE
c         constvar = 0
c      ENDIF

C     Make sure control dose level is first in unrestricted case
c     This is needed for donlp2usrfc.f to work correctly 
c      IF (restrict .EQ. 0) THEN
c         DO I = 2, ndoses 
c            IF (doses(I) .LT. doses(1)) THEN
c               temp1 = doses(1)
c               temp2 = mean(1)
c               temp3 = nanimals(1)
c               temp4 = var(1)
c               doses(1) = doses(I)
c               mean(1) = mean(I)
c               nanimals(1) = nanimals(I)
c               var(1) = var(I)
c               doses(I) = temp1
c               mean(I) = temp2
c               nanimals(I) = temp3
c               var(I) = temp4
c            ENDIF
c         ENDDO
c      ENDIF       


      RETURN
      END
