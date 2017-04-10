
*$ CREATE MGDRAW.FOR
*COPY MGDRAW
*                                                                      *
*=== mgdraw ===========================================================*
*                                                                      *
      SUBROUTINE MGDRAW ( ICODE, MREG )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1990-2005      by        Alfredo Ferrari           *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     MaGnetic field trajectory DRAWing: actually this entry manages   *
*                                        all trajectory dumping for    *
*                                        drawing                       *
*                                                                      *
*     Created by                     Alfredo Ferrari                   *
*                                    INFN - Milan                      *
*     last change  09-jul-05  by     Alfredo Ferrari                   *
*                                    INFN - Milan                      *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(CASLIM)'
      INCLUDE '(COMPUT)'
      INCLUDE '(SOURCM)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(GENSTK)'
      INCLUDE '(MGDDCM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(SUMCOU)'
      INCLUDE '(TRACKR)'
      INCLUDE '(SOUEVT)'
      COMMON /LOCALDEPOSIT/xsumde1

*
      CHARACTER*20 FILNAM
      LOGICAL LFCOPE
      SAVE LFCOPE
      DATA LFCOPE / .FALSE. /
      DATA RdErth/637814000./
*
*----------------------------------------------------------------------*
*                                                                      *
*     Icode = 1: call from Kaskad                                      *
*     Icode = 2: call from Emfsco                                      *
*     Icode = 3: call from Kasneu                                      *
*     Icode = 4: call from Kashea                                      *
*     Icode = 5: call from Kasoph                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
c      IF ( .NOT. LFCOPE ) THEN
c         LFCOPE = .TRUE.
c         IF ( KOMPUT .EQ. 2 ) THEN
c            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
c         ELSE
c            FILNAM = CFDRAW
c         END IF
c         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
c     &          'UNFORMATTED' )
c      END IF
c      WRITE (IODRAW) NTRACK, MTRACK, JTRACK, SNGL (ETRACK),
c     &               SNGL (WTRACK)
c      WRITE (IODRAW) ( SNGL (XTRACK (I)), SNGL (YTRACK (I)),
c     &                 SNGL (ZTRACK (I)), I = 0, NTRACK ),
c     &               ( SNGL (DTRACK (I)), I = 1,MTRACK ),
c     &                 SNGL (CTRACK)
d
d      do i=0,NTRACK-1
d        RADXY0 = SQRT(XTRACK(i)*XTRACK(i)+ YTRACK(i)*YTRACK(i))
d        RADXY1 = SQRT(XTRACK(i+1)*XTRACK(i+1)
d     &              + YTRACK(i+1)*YTRACK(i+1))
d        IF(RADXY0.LE.40.0 
d     & .AND. ZTRACK(I).GE.30.0.AND. ZTRACK(I).LE.31.0
d     & .AND. RADXY1.LE.40.0 
d     & .AND. ZTRACK(I+1).GE.30.0 .AND. ZTRACK(I+1).LE.31.0)  THEN
c         WRITE(99,*) JTRACK, ichrge(jtrack), 
c     & XTRACK(i), YTRACK(i),ZTRACK(I),DTRACK(I+1)
d       if(jtrack.gt.0) then 
d       if(ichrge(jtrack).ne.0)xsumde1=xsumde1+DTRACK(I+1)
                  
d        ENDIF
d        ENDIF
d      enddo
      
      IF(MREG.EQ.1.or.MREG.EQ.6.or.MREG.EQ.7.or.MREG.EQ.18.or.
     &  MREG.EQ.11.or.MREG.EQ.12.or.MREG.EQ.13.or.MREG.EQ.14 .or.
     &  MREG.EQ.15.or.MREG.EQ.16.or.MREG.EQ.17) THEN 
      DO I=0,NTRACK-1
        DO J = 1,MTRACK 
         WRITE (99,100) NCASE,MREG,MTRACK,
     &   JTRACK,ATRACK, etrack,
     &   XTRACK(0),YTRACK(0),ZTRACK(0),
     &   CXTRCK,CYTRCK,CZTRCK,
     &   DTRACK(J), 1
        ENDDO
      ENDDO
      ENDIF 
      RETURN
 100   FORMAT (1X,I7,2(1X,I3),
     & 1x,I4,2(1x,e12.6),
     & 3(1x,e12.6),
     & 3(1x,e12.6),  
     & 1x,e12.6,1x,I2)        
*
*======================================================================*
*                                                                      *
*     Boundary-(X)crossing DRAWing:                                    *
*                                                                      *
*     Icode = 1x: call from Kaskad                                     *
*             19: boundary crossing                                    *
*     Icode = 2x: call from Emfsco                                     *
*             29: boundary crossing                                    *
*     Icode = 3x: call from Kasneu                                     *
*             39: boundary crossing                                    *
*     Icode = 4x: call from Kashea                                     *
*             49: boundary crossing                                    *
*     Icode = 5x: call from Kasoph                                     *
*             59: boundary crossing                                    *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY BXDRAW ( ICODE, MREG, NEWREG, XSCO, YSCO, ZSCO )

          IF(
     &      (MREG.eq.19.and.NEWREG.eq.1).or.
     &      (MREG.eq.1.and.NEWREG.eq.19).or.
     &      (MREG.eq.19.and.NEWREG.eq.6).or.
     &      (MREG.eq.6.and.NEWREG.eq.19).or.
     &      (MREG.eq.19.and.NEWREG.eq.7).or.
     &      (MREG.eq.7.and.NEWREG.eq.19).or.
     &      (MREG.eq.19.and.NEWREG.eq.11).or.
     &      (MREG.eq.11.and.NEWREG.eq.19).or.
     &      (MREG.eq.19.and.NEWREG.eq.12).or.
     &      (MREG.eq.12.and.NEWREG.eq.19).or.
     &      (MREG.eq.10.and.NEWREG.eq.13).or.
     &      (MREG.eq.13.and.NEWREG.eq.10).or.
     &      (MREG.eq.19.and.NEWREG.eq.14).or.
     &      (MREG.eq.14.and.NEWREG.eq.19).or.
     &      (MREG.eq.19.and.NEWREG.eq.15).or.
     &      (MREG.eq.15.and.NEWREG.eq.19).or.
     &      (MREG.eq.19.and.NEWREG.eq.16).or.
     &      (MREG.eq.16.and.NEWREG.eq.19).or.
     &      (MREG.eq.19.and.NEWREG.eq.17).or.
     &      (MREG.eq.17.and.NEWREG.eq.19).or.
     &      (MREG.eq.19.and.NEWREG.eq.18).or.
     &      (MREG.eq.18.and.NEWREG.eq.19)
     &        ) THEN
c           WRITE(99,100) ncase,mreg,NEWREG,jtrack,
c     &     etrack-am(jtrack),
c     &     XSCO,YSCO,ZSCO,atrack,
c     &     CXTRCK,CYTRCK,CZTRCK
       endif
c      IF(MREG.EQ.1.and.NEWREG.eq.19) THEN 
c       WRITE (99,*) XTRACK(0),' ',YTRACK(0),' ',ZTRACK(0),' ',
c     & DTRACK(1),' ' ,ETRACK 
c       WRITE (99,*) XTRACK(1),' ',YTRACK(1),' ',ZTRACK(1),' ',
c     & DTRACK(1),' ' ,ETRACK
  
c      ENDIF
      RETURN


*
*======================================================================*
*                                                                      *
*     Event End DRAWing:                                               *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY EEDRAW ( ICODE )
      RETURN
*
*======================================================================*
*                                                                      *
*     ENergy deposition DRAWing:                                       *
*                                                                      *
*     Icode = 1x: call from Kaskad                                     *
*             10: elastic interaction recoil                           *
*             11: inelastic interaction recoil                         *
*             12: stopping particle                                    *
*             13: pseudo-neutron deposition                            *
*             14: escape                                               *
*             15: time kill                                            *
*     Icode = 2x: call from Emfsco                                     *
*             20: local energy deposition (i.e. photoelectric)         *
*             21: below threshold, iarg=1                              *
*             22: below threshold, iarg=2                              *
*             23: escape                                               *
*             24: time kill                                            *
*     Icode = 3x: call from Kasneu                                     *
*             30: target recoil                                        *
*             31: below threshold                                      *
*             32: escape                                               *
*             33: time kill                                            *
*     Icode = 4x: call from Kashea                                     *
*             40: escape                                               *
*             41: time kill                                            *
*             42: delta ray stack overflow                             *
*     Icode = 5x: call from Kasoph                                     *
*             50: optical photon absorption                            *
*             51: escape                                               *
*             52: time kill                                            *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY ENDRAW ( ICODE, MREG, RULL, XSCO, YSCO, ZSCO )
      
      IF(MREG.EQ.1.or.MREG.EQ.6.or.MREG.EQ.7.or.MREG.EQ.18.or.
     &  MREG.EQ.11.or.MREG.EQ.12.or.MREG.EQ.13.or.MREG.EQ.14 .or.
     &  MREG.EQ.15.or.MREG.EQ.16.or.MREG.EQ.17) THEN 
        WRITE (99,100) NCASE,MREG,MTRACK,
     &   JTRACK,ATRACK,etrack,
     &   XSCO,YSCO,ZSCO,
     &   CXTRCK,CYTRCK,CZTRCK,
     &   RULL, 0
      ENDIF 
      
c      IF ( .NOT. LFCOPE ) THEN
c         LFCOPE = .TRUE.
c         IF ( KOMPUT .EQ. 2 ) THEN
c            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
c         ELSE
c            FILNAM = CFDRAW
c         END IF
c         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
c     &          'UNFORMATTED' )
c      END IF
c      WRITE (IODRAW)  0, ICODE, JTRACK, SNGL (ETRACK), SNGL (WTRACK)
c      WRITE (IODRAW)  SNGL (XSCO), SNGL (YSCO), SNGL (ZSCO), SNGL (RULL)
      RETURN
*
*======================================================================*
*                                                                      *
*     SOurce particle DRAWing:                                         *
*                                                                      *
*======================================================================*
*
      ENTRY SODRAW
c        write (99,*) ICODE, MREG, RULL, XSCO, YSCO, ZSCO
c      IF ( .NOT. LFCOPE ) THEN
c         LFCOPE = .TRUE.
c         IF ( KOMPUT .EQ. 2 ) THEN
c            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
c         ELSE
c            FILNAM = CFDRAW
c         END IF
c         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
c     &          'UNFORMATTED' )
c      END IF
c      WRITE (IODRAW) -NCASE, NPFLKA, NSTMAX, SNGL (TKESUM),
c     &                SNGL (WEIPRI)
*  +-------------------------------------------------------------------*
*  |  Patch for heavy ions: it works only for 1 source particle on
*  |  the stack for the time being
c      IF ( ABS (ILOFLK (NPFLKA)) .GE. 10000 ) THEN
c         IONID = ILOFLK (NPFLKA)
c         CALL DCDION ( IONID )
c         WRITE (IODRAW) ( IONID,SNGL(TKEFLK(I)+AMNHEA(-IONID)),
c     &                    SNGL (WTFLK(I)), SNGL (XFLK (I)),
c     &                    SNGL (YFLK (I)), SNGL (ZFLK (I)),
c     &                    SNGL (TXFLK(I)), SNGL (TYFLK(I)),
c     &                    SNGL (TZFLK(I)), I = 1, NPFLKA )
*  |
*  +-------------------------------------------------------------------*
*  |
c      ELSE
c         WRITE (IODRAW) ( ILOFLK(I), SNGL (TKEFLK(I)+AM(ILOFLK(I))),
c     &                    SNGL (WTFLK(I)), SNGL (XFLK (I)),
c     &                    SNGL (YFLK (I)), SNGL (ZFLK (I)),
c     &                    SNGL (TXFLK(I)), SNGL (TYFLK(I)),
c     &                    SNGL (TZFLK(I)), I = 1, NPFLKA )
c      END IF
*  |
*  +-------------------------------------------------------------------*
      RETURN
*
*======================================================================*
*                                                                      *
*     USer dependent DRAWing:                                          *
*                                                                      *
*     Icode = 10x: call from Kaskad                                    *
*             100: elastic   interaction secondaries                   *
*             101: inelastic interaction secondaries                   *
*             102: particle decay  secondaries                         *
*             103: delta ray  generation secondaries                   *
*             104: pair production secondaries                         *
*             105: bremsstrahlung  secondaries                         *
*     Icode = 20x: call from Emfsco                                    *
*             208: bremsstrahlung secondaries                          *
*             210: Moller secondaries                                  *
*             212: Bhabha secondaries                                  *
*             214: in-flight annihilation secondaries                  *
*             215: annihilation at rest   secondaries                  *
*             217: pair production        secondaries                  *
*             219: Compton scattering     secondaries                  *
*             221: photoelectric          secondaries                  *
*             225: Rayleigh scattering    secondaries                  *
*     Icode = 30x: call from Kasneu                                    *
*             300: interaction secondaries                             *
*     Icode = 40x: call from Kashea                                    *
*             400: delta ray  generation secondaries                   *
*  For all interactions secondaries are put on GENSTK common (kp=1,np) *
*  but for KASHEA delta ray generation where only the secondary elec-  *
*  tron is present and stacked on FLKSTK common for kp=lstack          *
*                                                                      *
*======================================================================*
*
      ENTRY USDRAW ( ICODE, MREG, XSCO, YSCO, ZSCO )
c      IF ( .NOT. LFCOPE ) THEN
c         LFCOPE = .TRUE.
c         IF ( KOMPUT .EQ. 2 ) THEN
c            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
c         ELSE
c            FILNAM = CFDRAW
c         END IF
c         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
c     &          'UNFORMATTED' )
c      END IF
* No output by default:
      RETURN
*=== End of subrutine Mgdraw ==========================================*
      END

