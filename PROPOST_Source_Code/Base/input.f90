!-----------------------------------------------------------------------------------------------!
!    Read files                                                                                 !
!    Copyright (C) 2021  J. Baltazar                                                            !
!                                                                                               !
!    This program is free software: you can redistribute it and/or modify it under the terms of !
!    the GNU Affero General Public License as published by the Free Software Foundation, either !
!    version 3 of the License, or (at your option) any later version.                           !
!                                                                                               !
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;  !
!    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  !
!    See the GNU Affero General Public License for more details.                                !
!                                                                                               !
!    You should have received a copy of the GNU Affero General Public License                   !
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.                     !
!-----------------------------------------------------------------------------------------------!
SUBROUTINE INPUT(IDFILE,JJ)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 21052014, J. Baltazar, 2014 version 1.0                                        !
!    Modified  : 05062014, J. Baltazar, read PROPANEL_//IDENTU(JJ)//'_T'//ICHAR//'.DAT'         !
!    Modified  : 20012015, J. Baltazar, 2015 version 1.0                                        !
!    Modified  : 09032015, J. Baltazar, 2015 version 1.1 Unsteady Flow                          !
!    Modified  : 07042015, J. Baltazar, 2015 version 1.2                                        !
!    Modified  : 06072016, J. Baltazar, 2016 version 1.0                                        !
!    Modified  : 09022017, J. Baltazar, 2017 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPOST_MOD
IMPLICIT NONE
CHARACTER INFILE*50
CHARACTER TFILE*9,TFILEP*9,TFILEN*10,TFILEH*7,FPOINT*9,ZONE*4,WAKE*4,TZONE*5
CHARACTER*3  NUM,IFILE,JFILE,ICHAR
INTEGER :: I,J,K,JJ,TT,IDFILE
!-----------------------------------------------------------------------------------------------!
WRITE(ICHAR,'(I3)') NT
SELECT CASE (IDFILE)
!-----------------------------------------------------------------------------------------------!
CASE(11)
!-----------------------------------------------------------------------------------------------!
!    Read Panel Coordinates                                                                     !
!-----------------------------------------------------------------------------------------------!
   IF (NWA == 0) INFILE='PROPANEL.DAT'
   IF (NWA /= 0) INFILE='PROPANEL_'//IDENTU(JJ)//'_T'//TRIM(ADJUSTL(ICHAR))//'.DAT'
   OPEN(UNIT=11,FILE=INFILE,STATUS='UNKNOWN')
   READ(11,*)
   READ(11,*)
!-----------------------------------------------------------------------------------------------!
!    Read Panel Coordinates on the Blade                                                        !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      READ(11,*) ZONE,TFILE,NUM,FPOINT,IFILE,NCP1,JFILE,NRP1
      ALLOCATE(XP(NCP1,NRP1),YP(NCP1,NRP1),ZP(NCP1,NRP1),RP(NCP1,NRP1),TP(NCP1,NRP1))
      XP=0.D0
      YP=0.D0
      ZP=0.D0
      RP=0.D0
      TP=0.D0
      DO J=1,NRP1
         DO I=1,NCP1
            READ(11,*) XP(I,J),YP(I,J),ZP(I,J)
         END DO !I=1,NCP1
      END DO !J=1,NRP1
      DO K=2,NB
         READ(11,*)
         DO J=1,NRP1
            DO I=1,NCP1
               READ(11,*) WORK1,WORK1,WORK1
            END DO !I=1,NCP1
         END DO !J=1,NRP1
      END DO !K=2,NB
!-----------------------------------------------------------------------------------------------!
!    Compute Cylindrical Coordinates                                                            !
!-----------------------------------------------------------------------------------------------!
      RP=DSQRT(YP*YP+ZP*ZP)
      TP=DATAN2(ZP,YP)
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Read Panel Coordinates on the Blade Wake                                                   !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      READ(11,*) ZONE,TFILE,WAKE,NUM,FPOINT,IFILE,NPW1,JFILE,NRW1
      ALLOCATE(XPW(NPW1,NRW1),YPW(NPW1,NRW1),ZPW(NPW1,NRW1),RPW(NPW1,NRW1),TPW(NPW1,NRW1))
      XPW=0.D0
      YPW=0.D0
      ZPW=0.D0
      RPW=0.D0
      TPW=0.D0
      DO J=1,NRW1
         DO I=1,NPW1
            READ(11,*) XPW(I,J),YPW(I,J),ZPW(I,J)
         END DO !I=1,NPW1
      END DO !J=1,NRW1
      DO K=2,NB
         READ(11,*)
         DO J=1,NRW1
            DO I=1,NPW1
               READ(11,*) WORK1,WORK1,WORK1
            END DO !I=1,NPW1
         END DO !J=1,NRW1
      END DO !K=2,NB
!-----------------------------------------------------------------------------------------------!
!    Compute Cylindrical Coordinates                                                            !
!-----------------------------------------------------------------------------------------------!
      RPW=DSQRT(YPW*YPW+ZPW*ZPW)
      TPW=DATAN2(ZPW,YPW)
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Read Panel Coordinates on the Nozzle                                                       !
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      READ(11,*) ZONE,TFILE,NUM,FPOINT,IFILE,NNXT,JFILE,NNTT
      ALLOCATE(XN(NNXT,NNTT),YN(NNXT,NNTT),ZN(NNXT,NNTT),RN(NNXT,NNTT),TN(NNXT,NNTT))
      XN=0.D0
      YN=0.D0
      ZN=0.D0
      RN=0.D0
      TN=0.D0
      DO J=1,NNTT
         DO I=1,NNXT
            READ(11,*) XN(I,J),YN(I,J),ZN(I,J)
         END DO !I=1,NNXT
      END DO !J=1,NNTT
      DO K=2,NB
         READ(11,*)
         DO J=1,NNTT
            DO I=1,NNXT
               READ(11,*) WORK1,WORK1,WORK1
            END DO !I=1,NNXT
         END DO !J=1,NNTT
      END DO !K=2,NB
!-----------------------------------------------------------------------------------------------!
!    Compute Cylindrical Coordinates                                                            !
!-----------------------------------------------------------------------------------------------!
      RN=DSQRT(YN*YN+ZN*ZN)
      TN=DATAN2(ZN,YN)
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Read Panel Coordinates on the Nozzle Wake                                                  !
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      READ(11,*) ZONE,TFILE,WAKE,NUM,FPOINT,IFILE,NNW1,JFILE,NNTT
      ALLOCATE(XNW(NNW1,NNTT),YNW(NNW1,NNTT),ZNW(NNW1,NNTT),RNW(NNW1,NNTT),TNW(NNW1,NNTT))
      XNW=0.D0
      YNW=0.D0
      ZNW=0.D0
      RNW=0.D0
      TNW=0.D0
      DO J=1,NNTT
         DO I=1,NNW1
            READ(11,*) XNW(I,J),YNW(I,J),ZNW(I,J)
         END DO !I=1,NNW1
      END DO !J=1,NNTT
      DO K=2,NB
         READ(11,*)
         DO J=1,NNTT
            DO I=1,NNW1
               READ(11,*) WORK1,WORK1,WORK1
            END DO !I=1,NNW1
         END DO !J=1,NNTT
      END DO !K=2,NB
!-----------------------------------------------------------------------------------------------!
!    Compute Cylindrical Coordinates                                                            !
!-----------------------------------------------------------------------------------------------!
      RNW=DSQRT(YNW*YNW+ZNW*ZNW)
      TNW=DATAN2(ZNW,YNW)
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Read Panel Coordinates on the Hub                                                          !
!-----------------------------------------------------------------------------------------------!
   IF (IABS(IH) == 1) THEN
      READ(11,*) ZONE,TFILE,NUM,FPOINT,IFILE,NHX1,JFILE,NHTT
      ALLOCATE(XH(NHX1,NHTT),YH(NHX1,NHTT),ZH(NHX1,NHTT),RH(NHX1,NHTT),TH(NHX1,NHTT))
      XH=0.D0
      YH=0.D0
      ZH=0.D0
      RH=0.D0
      TH=0.D0
      DO J=1,NHTT
         DO I=1,NHX1
            READ(11,*) XH(I,J),YH(I,J),ZH(I,J)
         END DO !I=1,NHX1
      END DO !J=1,NHTT
      DO K=2,NB
         READ(11,*)
         DO J=1,NHTT
            DO I=1,NHX1
               READ(11,*) WORK1,WORK1,WORK1
            END DO !I=1,NHX1
         END DO !J=1,NHTT
      END DO !K=2,NB
!-----------------------------------------------------------------------------------------------!
!    Compute Cylindrical Coordinates                                                            !
!-----------------------------------------------------------------------------------------------!
      RH=DSQRT(YH*YH+ZH*ZH)
      TH=DATAN2(ZH,YH)
   END IF !(IABS(IH) == 1)
   CLOSE(UNIT=11)
!-----------------------------------------------------------------------------------------------!
CASE(21)
!-----------------------------------------------------------------------------------------------!
!    Read Geometry                                                                              !
!-----------------------------------------------------------------------------------------------!
   IF (NWA == 0) INFILE='GEOMETRY_'//IDENTU(JJ)//'.DAT'
   IF (NWA /= 0) INFILE='GEOMETRY_'//IDENTU(JJ)//'_T'//TRIM(ADJUSTL(ICHAR))//'.DAT'
   OPEN(UNIT=21,FILE=INFILE,STATUS='UNKNOWN')
   READ(21,*)
   READ(21,*)
!-----------------------------------------------------------------------------------------------!
!    Blade Geometry                                                                             !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      READ(21,*) ZONE,TFILEP,FPOINT,IFILE,NCP,JFILE,NRP
      ALLOCATE(XP0  (NCP,NRP),YP0  (NCP,NRP),ZP0  (NCP,NRP))
      ALLOCATE(UNXP0(NCP,NRP),UNYP0(NCP,NRP),UNZP0(NCP,NRP))
      ALLOCATE(ET1XP(NCP,NRP),ET1YP(NCP,NRP),ET1ZP(NCP,NRP))
      ALLOCATE(ET2XP(NCP,NRP),ET2YP(NCP,NRP),ET2ZP(NCP,NRP))
      ALLOCATE(AT1XP(NCP,NRP),AT1YP(NCP,NRP),AT1ZP(NCP,NRP))
      ALLOCATE(AT2XP(NCP,NRP),AT2YP(NCP,NRP),AT2ZP(NCP,NRP))
      ALLOCATE(AP0  (NCP,NRP))
      DO J=1,NRP
         DO I=1,NCP
            READ(21,*) XP0(I,J),YP0(I,J),ZP0(I,J),UNXP0(I,J),UNYP0(I,J),UNZP0(I,J), &
                 ET1XP(I,J),ET1YP(I,J),ET1ZP(I,J),ET2XP(I,J),ET2YP(I,J),ET2ZP(I,J), &
                 AT1XP(I,J),AT1YP(I,J),AT1ZP(I,J),AT2XP(I,J),AT2YP(I,J),AT2ZP(I,J),AP0(I,J)
         END DO !I=1,NCP
      END DO !J=1,NRP
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Geometry                                                                            !
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      READ(21,*) ZONE,TFILEN,FPOINT,IFILE,NNXT1,JFILE,NNTP
      ALLOCATE(XN0  (NNXT1,NNTP),YN0  (NNXT1,NNTP),ZN0  (NNXT1,NNTP))
      ALLOCATE(UNXN0(NNXT1,NNTP),UNYN0(NNXT1,NNTP),UNZN0(NNXT1,NNTP))
      ALLOCATE(ET1XN(NNXT1,NNTP),ET1YN(NNXT1,NNTP),ET1ZN(NNXT1,NNTP))
      ALLOCATE(ET2XN(NNXT1,NNTP),ET2YN(NNXT1,NNTP),ET2ZN(NNXT1,NNTP))
      ALLOCATE(AT1XN(NNXT1,NNTP),AT1YN(NNXT1,NNTP),AT1ZN(NNXT1,NNTP))
      ALLOCATE(AT2XN(NNXT1,NNTP),AT2YN(NNXT1,NNTP),AT2ZN(NNXT1,NNTP))
      ALLOCATE(AN0  (NNXT1,NNTP))
      DO J=1,NNTP
         DO I=1,NNXT1
            READ(21,*) XN0(I,J),YN0(I,J),ZN0(I,J),UNXN0(I,J),UNYN0(I,J),UNZN0(I,J), &
                 ET1XN(I,J),ET1YN(I,J),ET1ZN(I,J),ET2XN(I,J),ET2YN(I,J),ET2ZN(I,J), &
                 AT1XN(I,J),AT1YN(I,J),AT1ZN(I,J),AT2XN(I,J),AT2YN(I,J),AT2ZN(I,J),AN0(I,J)
         END DO !I=1,NNXT1
      END DO !J=1,NNTP
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub Geometry                                                                               !
!-----------------------------------------------------------------------------------------------!
   IF (IABS(IH) == 1) THEN
      READ(21,*) ZONE,TFILEH,FPOINT,IFILE,NHX,JFILE,NHTP
      ALLOCATE(XH0  (NHX,NHTP),YH0  (NHX,NHTP),ZH0  (NHX,NHTP))
      ALLOCATE(UNXH0(NHX,NHTP),UNYH0(NHX,NHTP),UNZH0(NHX,NHTP))
      ALLOCATE(ET1XH(NHX,NHTP),ET1YH(NHX,NHTP),ET1ZH(NHX,NHTP))
      ALLOCATE(ET2XH(NHX,NHTP),ET2YH(NHX,NHTP),ET2ZH(NHX,NHTP))
      ALLOCATE(AT1XH(NHX,NHTP),AT1YH(NHX,NHTP),AT1ZH(NHX,NHTP))
      ALLOCATE(AT2XH(NHX,NHTP),AT2YH(NHX,NHTP),AT2ZH(NHX,NHTP))
      ALLOCATE(AH0  (NHX,NHTP))
      DO J=1,NHTP
         DO I=1,NHX
            READ(21,*) XH0(I,J),YH0(I,J),ZH0(I,J),UNXH0(I,J),UNYH0(I,J),UNZH0(I,J), &
                 ET1XH(I,J),ET1YH(I,J),ET1ZH(I,J),ET2XH(I,J),ET2YH(I,J),ET2ZH(I,J), &
                 AT1XH(I,J),AT1YH(I,J),AT1ZH(I,J),AT2XH(I,J),AT2YH(I,J),AT2ZH(I,J),AH0(I,J)
         END DO !I=1,NHX
      END DO !J=1,NHTP
   END IF !(IABS(IH) == 1)
   CLOSE(UNIT=21)
!-----------------------------------------------------------------------------------------------!
CASE(22)
!-----------------------------------------------------------------------------------------------!
!    Read Wake Solution                                                                         !
!-----------------------------------------------------------------------------------------------!
   INFILE='CIRC_'//IDENTU(JJ)//'.DAT'
   OPEN(UNIT=22,FILE=INFILE,STATUS='UNKNOWN')
   READ(22,*)
   READ(22,*)
   IF (IP == 1) THEN
      ALLOCATE(RR(NRW),POTPW(NPW,NRW,0:NT))
      POTPW=0.D0
   END IF !(IP == 1)
   DO TT=0,NT
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
      IF (IP == 1) THEN
         READ(22,*) ZONE,TFILE,WAKE,TZONE,IFILE,NRW,FPOINT
         DO J=1,NRW
            READ(22,*) RR(J),POTPW(1,J,TT),WORK1,WORK1,WORK1,WORK1,WORK1
         END DO !J=1,NRW
      END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
      IF (IN == 1) THEN
         READ(22,*) ZONE,TFILEN,WAKE,TZONE,IFILE,NNTP,FPOINT
         DO J=1,NNTP
            READ(22,*) WORK1,WORK1,WORK1,WORK1,WORK1,WORK1,WORK1
         END DO !J=1,NNTP
      END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
   END DO !TT=0,NT
   CLOSE(UNIT=22)
!-----------------------------------------------------------------------------------------------!
CASE(23)
!-----------------------------------------------------------------------------------------------!
!    Read Solution in Tecplot Format                                                            !
!-----------------------------------------------------------------------------------------------!
   INFILE='SOL_'//IDENTU(JJ)//'.DAT'
   OPEN(UNIT=23,FILE=INFILE,STATUS='UNKNOWN')
   READ(23,*)
   READ(23,*)
   IF (IP == 1) THEN
      ALLOCATE(POTP(NCP,NRP,0:NT),SOURCEP(NCP,NRP,0:NT))
      ALLOCATE(CPP(NCP,NRP,0:NT),CPNP(NCP,NRP,0:NT))
      ALLOCATE(VT1P(NCP,NRP,0:NT),VT2P(NCP,NRP,0:NT))
      ALLOCATE(VTT1P(NCP,NRP,0:NT),VTT2P(NCP,NRP,0:NT))
   END IF !(IP == 1)
   IF ((IP == 1).AND.(NCPW /= 0)) THEN
      ALLOCATE(XPW0(IABS(NCPW),NRW),YPW0(IABS(NCPW),NRW),ZPW0(IABS(NCPW),NRW))
      ALLOCATE(SOURCEPWCAV(IABS(NCPW),NRW,0:NT),VT1PW(IABS(NCPW),NRW))
      IF (NCPW < 0) ALLOCATE(POTPWP(IABS(NCPW),NRW,0:NT))
      IF (NCPW > 0) ALLOCATE(POTPWS(IABS(NCPW),NRW,0:NT))
   END IF !((IP == 1).AND.(NCPW /= 0))
   IF (IN == 1) THEN
      ALLOCATE(POTN(NNXT1,NNTP,0:NT),SOURCEN(NNXT1,NNTP,0:NT))
      ALLOCATE(CPN(NNXT1,NNTP,0:NT),CPNN(NNXT1,NNTP,0:NT))
      ALLOCATE(VT1N(NNXT1,NNTP,0:NT),VT2N(NNXT1,NNTP,0:NT))
      ALLOCATE(VTT1N(NNXT1,NNTP,0:NT),VTT2N(NNXT1,NNTP,0:NT))
   END IF !(IN == 1)
   IF (IABS(IH) == 1) THEN
      ALLOCATE(POTH(NHX,NHTP,0:NT),SOURCEH(NHX,NHTP,0:NT))
      ALLOCATE(CPH(NHX,NHTP,0:NT),CPNH(NHX,NHTP,0:NT))
      ALLOCATE(VT1H(NHX,NHTP,0:NT),VT2H(NHX,NHTP,0:NT))
      ALLOCATE(VTT1H(NHX,NHTP,0:NT),VTT2H(NHX,NHTP,0:NT))
   END IF !(IABS(IH) == 1)
   DO TT=0,NT
!-----------------------------------------------------------------------------------------------!
!    Blade Solution                                                                             !
!-----------------------------------------------------------------------------------------------!
      IF (IP == 1) THEN
         READ(23,*) ZONE,TFILEP,TZONE,FPOINT,IFILE,NCP,JFILE,NRP
         DO J=1,NRP
            DO I=1,NCP
               READ(23,*) XP0(I,J),YP0(I,J),ZP0(I,J),POTP(I,J,TT),SOURCEP(I,J,TT), &
                          VT1P(I,J,TT),VT2P(I,J,TT),VTT1P(I,J,TT),VTT2P(I,J,TT),CPP(I,J,TT), &
                          CPNP(I,J,TT)
            END DO !I=1,NCP
         END DO !J=1,NRP
      END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake Solution                                                                        !
!-----------------------------------------------------------------------------------------------!
      IF ((IP == 1).AND.(NCPW < 0)) THEN
         READ(23,*)
         DO J=1,NRW
            DO I=1,IABS(NCPW)
               READ(23,*) XPW0(I,J),YPW0(I,J),ZPW0(I,J),POTPWP(I,J,TT), &
                          SOURCEPWCAV(I,J,TT),VT1PW(I,J),WORK1,WORK1,WORK1,WORK1,WORK1
            END DO !I=1,IABS(NCPW)
         END DO !J=1,NRW
      END IF !((IP == 1).AND.(NCPW < 0))
      IF ((IP == 1).AND.(NCPW > 0)) THEN
         READ(23,*)
         DO J=1,NRW
            DO I=1,IABS(NCPW)
               READ(23,*) XPW0(I,J),YPW0(I,J),ZPW0(I,J),POTPWS(I,J,TT), &
                          SOURCEPWCAV(I,J,TT),VT1PW(I,J),WORK1,WORK1,WORK1,WORK1,WORK1
            END DO !I=1,IABS(NCPW)
         END DO !J=1,NRW
      END IF !((IP == 1).AND.(NCPW > 0))
!-----------------------------------------------------------------------------------------------!
!    Nozzle Solution                                                                            !
!-----------------------------------------------------------------------------------------------!
      IF (IN == 1) THEN
         READ(23,*) ZONE,TFILEN,TZONE,FPOINT,IFILE,NNXT1,JFILE,NNTP
         DO J=1,NNTP
            DO I=1,NNXT1
               READ(23,*) XN0(I,J),YN0(I,J),ZN0(I,J),POTN(I,J,TT),SOURCEN(I,J,TT), &
                          VT1N(I,J,TT),VT2N(I,J,TT),VTT1N(I,J,TT),VTT2N(I,J,TT),CPN(I,J,TT), &
                          CPNN(I,J,TT)
            END DO !I=1,NNXT1
         END DO !J=1,NNTP
      END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub Solution                                                                               !
!-----------------------------------------------------------------------------------------------!
      IF (IABS(IH) == 1) THEN
         READ(23,*) ZONE,TFILEH,TZONE,FPOINT,IFILE,NHX,JFILE,NHTP
         DO J=1,NHTP
            DO I=1,NHX
               READ(23,*) XH0(I,J),YH0(I,J),ZH0(I,J),POTH(I,J,TT),SOURCEH(I,J,TT), &
                          VT1H(I,J,TT),VT2H(I,J,TT),VTT1H(I,J,TT),VTT2H(I,J,TT),CPH(I,J,TT), &
                          CPNH(I,J,TT)
            END DO !I=1,NHX
         END DO !J=1,NHTP
      END IF !(IABS(IH) == 1)
!-----------------------------------------------------------------------------------------------!
   END DO !TT=0,NT
   CLOSE(UNIT=23)
!-----------------------------------------------------------------------------------------------!
END SELECT !(IDFILE)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE INPUT
!-----------------------------------------------------------------------------------------------!
