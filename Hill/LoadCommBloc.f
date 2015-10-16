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
      SUBROUTINE LOADCOMMBLOC(XNDOSES,XDOSES,XMEAN,XNANIM,XVAR,
     $     XNPARMS,XFIXD,XVAL,XRESTR,XADVERSE,XMODEL,XMAX,XMIN)

      INCLUDE 'PROBLEM.INC'     
      INTEGER XNDOSES,XNANIM(*),XNPARMS,XFIXD(0:XNPARMS-1),
     $     XRESTR,XADVERSE,XMODEL,
     $     I,temp3
      DOUBLE PRECISION XDOSES(*),XMEAN(*),XVAR(*),XVAL(0:XNPARMS-1),
     $     XMAX,XMIN,
     $     temp1,temp2,temp4

      mlecnt = 0
      clcnt = 0
      a3cnt = 0
      ndoses = XNDOSES
      nparm = XNPARMS
      restrict = XRESTR
      adverse = XADVERSE
      modtype = XMODEL
      maxdose = XMAX
      mindose = XMIN
      DO I = 1, ndoses
         doses(I) = XDOSES(I)
         mean(I) = XMEAN(I)
         nanimals(I) = XNANIM(I)
         var(I) = XVAR(I)
      ENDDO
      DO I = 0, nparm-1
         parmfixd(I) = XFIXD(I)
         parmval(I) = XVAL(I)
      ENDDO
      IF(parmfixd(1).EQ.1.AND.parmval(1).EQ.0) THEN
         constvar = 1
      ELSE
         constvar = 0
      ENDIF

C     Make sure control dose level is first in unrestricted case
c     This is needed for donlp2usrfc.f to work correctly 
      IF (restrict .EQ. 0) THEN
         DO I = 2, ndoses 
            IF (doses(I) .LT. doses(1)) THEN
               temp1 = doses(1)
               temp2 = mean(1)
               temp3 = nanimals(1)
               temp4 = var(1)
               doses(1) = doses(I)
               mean(1) = mean(I)
               nanimals(1) = nanimals(I)
               var(1) = var(I)
               doses(I) = temp1
               mean(I) = temp2
               nanimals(I) = temp3
               var(I) = temp4
            ENDIF
         ENDDO
      ENDIF       


      RETURN
      END
