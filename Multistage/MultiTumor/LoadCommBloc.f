c     LOADCOMMBLOC : Loads values into the common block for
c     donlp2.
c     INPUT:
c     XMAX       the max dose level (DOUBLE)
c
c
      SUBROUTINE LOADCOMMBLOC(XNDOSES,XMAX,XNANIM,XDOSES,XAFFECT,
     $     XPOLYORD,XRESTR, XLMIND, XLMAXD, XNT)
     
      INCLUDE 'PROBLEM.INC'
      INCLUDE 'MULTITUMOR.INC'
      INTEGER XNDOSES,XAFFECT(*),XNANIM(*),XPOLYORD,XRESTR, XNT
c     $     ,I,temp2,temp3
      DOUBLE PRECISION XMAX,XDOSES(*),temp1, XLMIND, XLMAXD
      INTEGER I,temp2,temp3
      maxdose = XMAX
      polyord = XPOLYORD
      ndoses = XNDOSES
      restrict = XRESTR
      lminbmd = XLMIND
      lmaxbmd = XLMAXD
      Print *,"Inside Load Common BLock"
      Print *,"I, Dose, Affected, Number Animals"
      DO I = 1, ndoses
         doses(I) = XDOSES(I)
         affected(I) = XAFFECT(I)
         nanimals(I) = XNANIM(I)
         mdDoses(XNT,I) = XDOSES(I)
         mnAffected(XNT,I) = XAFFECT(I)
         mnAnim(XNT,I) = XNANIM(I)
      Print *,I,doses(I),affected(i),nanimals(I)
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


c     LOADCOMMBL2 : Loads values into the common block for
c     donlp2.
c     INPUT:
c     XMAX       the max dose level (DOUBLE)
c
c
      SUBROUTINE MULTILOADCOMMBLOC(lnT, lnObs, lnMaxNParms, 
     $     anXDOSES, anXPOLYORD, anXRESTR, anParms)

      INCLUDE 'PROBLEM.INC'  
      INCLUDE 'MULTITUMOR.INC'
      INTEGER lnT, lnObs, lnMaxNParms, anXDOSES(0:9), 
     $     anXPOLYORD(0:9), anXRESTR(0:9),
     $     anParms(0:9)
      INTEGER I, O, temp2, temp3, ndose, temp1
      nT = lnT
      nMaxObs = lnObs
      nMaxNParms = lnMaxNParms
      Print *,"Inside MULTILOADCOMMBLOC"
      Print *,"I, Obs, PolyDeg, Number Parms"      
      DO I = 0, lnT-1
         arnObs(I) = anXDOSES(I)
         arnPoly(I) = anXPOLYORD(I)
         arnParms(I) = anParms(I)
         arnRestrict(I) = anXRESTR(I)
         Print *,I,arnObs(I),arnPoly(I),arnParms(I)
      ENDDO
      Print *,"I, Dose, Affected, Number Animals"
      DO I = 0, nT-1
         DO O = 1, nMaxObs
            Print *,I,mdDoses(I,O),mnAffected(I,O),mnAnim(I,O)
         ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE PARMSLOADCOMMBLOC(XnT, XPARMS)
      INCLUDE 'MULTITUMOR.INC'
      INTEGER Xnt, I
      DOUBLE PRECISION XPARMS(0:*)
      Print *,"XnT=",Xnt,"  arnPoly(",XnT,")=",arnPoly(XnT)
      DO I=0, arnPoly(XnT)
         mdParmVal(XnT, I) = XPARMS(I)
         print *, 'XPARMS(',I,')=',XPARMS(I)
      ENDDO
C   added by gln 02/04/09 start Check values of XPARMS ******
      Print *,"PARMSLOADCOMMBLOC, Tumor ",XnT
      DO I=0, arnPoly(XnT)
         print *, 'XPARMS(',I,')=',mdParmVal(XnT, I)
      ENDDO
      print *
C   added by gln 02/04/09 end Check values of XPARMS ******
      RETURN
      END

      SUBROUTINE GETPARMSCOMMBLOC(XnT, XPARMS)
      INCLUDE 'MULTITUMOR.INC'
      INTEGER Xnt, I
      DOUBLE PRECISION XPARMS(*)
      DO I=0, arnPoly(XnT)
         XPARMS(I) = mdParmVal(XnT, I)
      ENDDO
      RETURN
      END      