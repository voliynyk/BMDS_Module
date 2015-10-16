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
      INTEGER XNDOSES,XPOLYORD,XRESTR, XNT
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
      if (geLogLevel .ge. LOG_DEVELOPER) then
         Print *,"Inside Load Common BLock"
         Print *,"maxdose =", XMAX
         Print *,"polyord =", XPOLYORD
         Print *,"ndoses =", XNDOSES
         Print *,"restrict =", XRESTR
         Print *,"lminbmd =", XLMIND
         Print *,"lmaxbmd =", XLMAXD
      endif

      DO I = 1, ndoses
         doses(I) = XDOSES(I)
         affected(I) = XAFFECT(I)
         nanimals(I) = XNANIM(I)
         mdDoses(XNT,I) = XDOSES(I)
         mnAffected(XNT,I) = XAFFECT(I)
         mnAnim(XNT,I) = XNANIM(I)
      ENDDO

      if (geLogLevel .ge. LOG_DEVELOPER) then
         Print *,"I, Dose, Affected, Number Animals"
         DO I = 1, ndoses
            Print *,I,doses(I),affected(i),nanimals(I)
         ENDDO
      endif

C     Make sure control dose level is first in unrestricted case
c     This is needed for donlp2usrfc.f to work correctly 
      IF (restrict .NE. 1) THEN
         DO I = 2, ndoses 
            IF (doses(I) .LT. doses(1)) THEN
c         		print *,"swapping doses(",I,") < doses(1)"
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


c     MULTILOADCOMMBLOC : Loads values into the common block for
c     donlp2.
C     Copy tumor dose group info into fortran common data structures
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
      INTEGER I, O, J, ndose, iSwaps,
     $     aiAffected(1:100), aiAnimals(1:100)
      DOUBLE PRECISION adDoses(1:100)

      nT = lnT
      nMaxObs = lnObs
      nMaxNParms = lnMaxNParms
      if (geLogLevel .ge. LOG_DEVELOPER) then
         Print *,"Inside MULTILOADCOMMBLOC"
         Print *,"nT =", lnT
         Print *,"nMaxObs =", lnObs
         Print *,"nMaxNParms =", lnMaxNParms
      endif
C      Print *,"I, Obs, PolyDeg, Number Parms, Restrict"
C      DO I = 0, lnT-1
C         arnObs(I) = anXDOSES(I)
C         arnPoly(I) = anXPOLYORD(I)
C         arnParms(I) = anParms(I)
C         arnRestrict(I) = anXRESTR(I)
C         Print *,I,arnObs(I),arnPoly(I),arnParms(I),arnRestrict(I)
C      ENDDO

C     LCO 3/2010 - Reverse order to match the parms passed later into getclmt
      o = lnT-1
      DO I = 0, lnT-1
         arnObs(I) = anXDOSES(o)
         arnPoly(I) = anXPOLYORD(o)
         arnParms(I) = anParms(o)
         arnRestrict(I) = anXRESTR(o)
C         Print *,I,arnObs(I),arnPoly(I),arnParms(I),arnRestrict(I)
         o = o - 1
      ENDDO

      if (geLogLevel .ge. LOG_DEVELOPER) then
         Print *,"I, O, Dose(I,O), Affected(I,O), #Animals(I,O)"
         DO I = 0, nT-1
            DO O = 1, arnObs(I)
               Print *,I,O,mdDoses(I,O),mnAffected(I,O),mnAnim(I,O)
            ENDDO
         ENDDO
      endif
      RETURN
      END

      SUBROUTINE PARMSLOADCOMMBLOC(XnT, XPARMS)
      INCLUDE 'MULTITUMOR.INC'
      INTEGER Xnt, I
      DOUBLE PRECISION XPARMS(0:*)
C      Print *,"Entered parmsloadcommbloc"
C      Print *,"XnT=",Xnt,"  arnPoly(",XnT,")=",arnPoly(XnT)
      DO I=0, arnPoly(XnT)
         mdParmVal(XnT, I) = XPARMS(I)
C         print *, "mdParmVal(",XnT,",",I,") =XPARMS(",I,")=",XPARMS(I)
      ENDDO

      if (geLogLevel .ge. LOG_DEVELOPER) then
         Print *,"PARMSLOADCOMMBLOC, Tumor ",XnT
         DO I=0, arnPoly(XnT)
            print *, 'XPARMS(',I,')=',mdParmVal(XnT, I)
         ENDDO
         print *
      endif

      RETURN
      END

      SUBROUTINE GETPARMSCOMMBLOC(XnT, XPARMS)
      INCLUDE 'MULTITUMOR.INC'
      INTEGER Xnt, I
      DOUBLE PRECISION XPARMS(*)
c      print *,"Entered getparmscommbloc: XnT=", XnT
      DO I=0, arnPoly(XnT)
         XPARMS(I) = mdParmVal(XnT, I)
c         print *,"XPARMS(",I,")=mdParmVal(",XnT,",",I,") = ",
c     1         mdParmVal(XnT, I)
      ENDDO
      RETURN
      END      
