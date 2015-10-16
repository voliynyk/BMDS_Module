c     LOADCOMMBLOC : Loads values into the common block for
c     donlp2.
c     INPUT:
c     XMAX       the max dose level (DOUBLE)
c
c
      SUBROUTINE LOADCOMMBLOC(XNDOSES,XMAX,XNANIM,XDOSES,XAFFECT,
     $     XPOLYORD,XRESTR, XLMIND, XLMAXD)

      INCLUDE 'PROBLEM.INC'     
      INTEGER XNDOSES,XPOLYORD,XRESTR
      DOUBLE PRECISION XAFFECT(*),XNANIM(*)
      DOUBLE PRECISION XMAX,XDOSES(*), XLMIND, XLMAXD
      INTEGER I
      DOUBLE PRECISION temp1,temp2,temp3
      maxdose = XMAX
      polyord = XPOLYORD
      ndoses = XNDOSES
      restrict = XRESTR
      lminbmd = XLMIND
      lmaxbmd = XLMAXD
      DO I = 1, ndoses
         doses(I) = XDOSES(I)
         affected(I) = XAFFECT(I)
         nanimals(I) = XNANIM(I)
      ENDDO

C     Make sure control dose level is first in unrestricted case
c     This is needed for donlp2usrfc.f to work correctly 
      IF (restrict .NE. 1) THEN
         DO I = 2, ndoses 
            IF (doses(I) .LT. doses(1)) THEN
               temp1 = doses(1)
               temp2 = affected(1)
               temp3 = nanimals(1)
               doses(1) = doses(I)
               affected(1) = affected(I)
               nanimals(1) = nanimals(I)
               doses(I) = temp1
               affected(I) = temp2
               nanimals(I) = temp3
            ENDIF
         ENDDO
      ENDIF       

      RETURN
      END
