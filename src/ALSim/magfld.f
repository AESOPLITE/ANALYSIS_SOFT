*
*===magfld=============================================================*
*
      SUBROUTINE MAGFLD ( X, Y, Z, BTX, BTY, BTZ, B, NREG, IDISC )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1988-2010      by Alberto Fasso` & Alfredo Ferrari *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     Created  in     1988         by    Alberto Fasso`                *
*                                                                      *
*                                                                      *
*     Last change on 06-Nov-10     by    Alfredo Ferrari               *
*                                                                      *
*     Input variables:                                                 *
*            x,y,z = current position                                  *
*            nreg  = current region                                    *
*     Output variables:                                                *
*            btx,bty,btz = cosines of the magn. field vector           *
*            B = magnetic field intensity (Tesla)                      *
*            idisc = set to 1 if the particle has to be discarded      *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(CMEMFL)'
      INCLUDE '(CSMCRY)'
      LOGICAL FIRSTCALL
      DATA FIRSTCALL/.TRUE./ 
      DIMENSION XX(81,81,81),YY(81,81,81),ZZ(81,81,81)
      DIMENSION BXX(81,81,81),BYY(81,81,81),BZZ(81,81,81)
      DIMENSION BTT(81,81,81)
        

*
*  +-------------------------------------------------------------------*
*  |  Earth geomagnetic field:
      IF ( LGMFLD ) THEN
         CALL GEOFLD ( X, Y, Z, BTX, BTY, BTZ, B, NREG, IDISC )
         RETURN
      END IF
      IF ( firstcall ) THEN
      firstcall=.FALSE.
      OPEN(56,
     & FILE='/home/psm/Documents/UDEL/AESOPLITE/fieldmap.txt',
     & STATUS='OLD')
      DO II=1,81
       DO JJ=1,81
        DO KK=1,81
           READ(56,*) XX(II,JJ,KK),YY(II,JJ,KK),ZZ(II,JJ,KK),
     &   BXX(II,JJ,KK),BYY(II,JJ,KK),BZZ(II,JJ,KK)
         XX(II,JJ,KK)=XX(II,JJ,KK)/10.d0
         YY(II,JJ,KK)=YY(II,JJ,KK)/10.d0
         ZZ(II,JJ,KK)=ZZ(II,JJ,KK)/10.d0

        BTT(II,JJ,KK)
     &  =SQRT(BXX(II,JJ,KK)**2 + BYY(II,JJ,KK)**2 + BZZ(II,JJ,KK)**2)
        BXX(II,JJ,KK) = BXX(II,JJ,KK)/BTT(II,JJ,KK)
        BYY(II,JJ,KK) = BYY(II,JJ,KK)/BTT(II,JJ,KK)
        BZZ(II,JJ,KK) = BZZ(II,JJ,KK)/BTT(II,JJ,KK)
        BTT(II,JJ,KK) = BTT(II,JJ,KK)/1.0D4

        ENDDO
       ENDDO
      ENDDO  
      close(56)
      ENDIF


*  |
*  +-------------------------------------------------------------------*
      IDISC = 0
      IXX = INT(2.D0 * (X + 20.0D0)) + 1
      IYY = INT(2.D0 * (Y + 20.0D0)) + 1
      IZZ = INT(2.D0 * (Z + 29.5D0)) + 1
      IF(IXX.LT.1.OR.IYY.LT.1.OR.IZZ.LT.1.OR.
     &   IXX.GT.81.OR.IYY.GT.81.OR.IZZ.GT.81) THEN
       BTX   = 0.0D0
       BTY   = 0.0D0
       BTZ   = 0.0D0
       B     = 0.0D0
      ELSE
       BTX   = BXX(IXX,IYY,IZZ)
       BTY   = BYY(IXX,IYY,IZZ)
       BTZ   = BZZ(IXX,IYY,IZZ)
       B     = BTT(IXX,IYY,IZZ)
      ENDIF
      RETURN
*=== End of subroutine Magfld =========================================*
      END