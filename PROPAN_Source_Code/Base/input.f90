!-----------------------------------------------------------------------------------------------!
!    Input File                                                                                 !
!    Copyright (C) 2021  J. Baltazar and J.A.C. Falcão de Campos                                !
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
SUBROUTINE INPUT(INPFILE)
!-----------------------------------------------------------------------------------------------!
!    Modified  : 03122013, J. Baltazar, version 1.0                                             !
!    Modified  : 29052014, J. Baltazar, Wake Alignment Model 1                                  !
!    Modified  : 09102014, J. Baltazar, Wake Alignment Model 2                                  !
!    Modified  : 17102014, J. Baltazar, Status of allocated arrays                              !
!    Modified  : 28102014, J. Baltazar, Interpolation scheme at nozzle t.e.                     !
!    Modified  : 30102014, J. Baltazar, version 2.2, Interpolation scheme at nozzle t.e.        !
!    Modified  : 05112014, J. Baltazar, version 3.0, Steady Cavitation Model                    !
!    Modified  : 18112014, J. Baltazar, version 3.2, Mid-Chord Cavitation Model                 !
!    Modified  : 19112014, J. Baltazar, version 3.3, Steady Super-Cavitation Model              !
!    Modified  : 12122014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!    Modified  : 27052016, J. Baltazar, 2016 version 1.2                                        !
!    Modified  : 22062016, J. Baltazar, 2016 version 1.3                                        !
!    Modified  : 20102016, J. Baltazar, 2016 version 1.4, Broyden's Method for IPKC             !
!    Modified  : 03012018, J. Baltazar, 2018 version 1.0, Gap strip alignment                   !
!    Modified  : 02102018, J. Baltazar, 2018 version 1.1, Duct Kutta points at inner and outer  !
!    Modified  : 04112020, J. Baltazar, 2020 version 1.2, Dynamic inflow                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
CHARACTER*50 INFILE
CHARACTER*10 TFILEN
CHARACTER*9  TFILEP,TFILE,FPOINT
CHARACTER*7  TFILEH
CHARACTER*5  TZONE
CHARACTER*4  ZONE,WAKE
CHARACTER*3  NUM,IFILE,JFILE
INTEGER :: I,J,K,INPFILE
NAMELIST /PROJ_INPUT/ COMMENT,IP,IN,IH,NB,JI,JF,IROTOR,NU,IDENTU,UU,NKIT,MKIT,IK,TOLK,BETA,     &
                      XITE,XOTE,NWA,IWA,NWALIGN,ITE,NRV,JV,NWV,XV,NBASIS,CCI,CCF,EPS,IFARV,     &
                      INTERPW,IBD,DELTA,NN,ALFV,IHUB,ITIP,REXTRA,IDENTN,LD,CR,IGRIDI,IGRIDO,    &
                      INTERN,NRNI,XIL,YIL,NRNO,XOL,YOL,INTERH,NHI,XHI,RHI,NHP,XHP,PTH,PGAP,     &
                      PREX,NGAP,CQ,FREX,TOLG,NCAV,SIGMA,FN,IRED,IMCP,IMCS,NCPW,CREX,TOLC,NREV,  &
                      NTETA,IFREQ,IFILESOL,IPAN,INTE,TOLS,MMAX,IFARP,ISTRIP,ISOLVER,IFIELD,ICP
!-----------------------------------------------------------------------------------------------!
SELECT CASE (INPFILE)
!-----------------------------------------------------------------------------------------------!
CASE(10)
!-----------------------------------------------------------------------------------------------!
!    Operational Conditions and Control Data                                                    !
!-----------------------------------------------------------------------------------------------!
   INFILE='PROPAN.INP'
   OPEN(UNIT=10,FILE=INFILE,STATUS='UNKNOWN')
   READ(10,PROJ_INPUT)
   CLOSE(UNIT=10)
!-----------------------------------------------------------------------------------------------!
   IF (IROTOR == 1) THEN !Adapt to Code Nomenclature
      WORK1=IMCP
      IMCP =IMCS
      IMCS =WORK1
   END IF !(IROTOR == 1)
!-----------------------------------------------------------------------------------------------!
CASE(11)
!-----------------------------------------------------------------------------------------------!
!    Read Panel Coordinates                                                                     !
!-----------------------------------------------------------------------------------------------!
   INFILE='PROPANEL.DAT'
   OPEN(UNIT=11,FILE=INFILE,STATUS='UNKNOWN')
   READ(11,*)
   READ(11,*)
!-----------------------------------------------------------------------------------------------!
!    Read Panel Coordinates on the Blade                                                        !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      READ(11,*) ZONE,TFILE,NUM,FPOINT,IFILE,NCP1,JFILE,NRP1
      IF (.NOT. ALLOCATED(XP)) &
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
      IF (.NOT. ALLOCATED(XPW)) &
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
      IF (.NOT. ALLOCATED(XN)) &
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
      IF (.NOT. ALLOCATED(XNW)) &
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
      IF (.NOT. ALLOCATED(XH)) &
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
CASE(12)
!-----------------------------------------------------------------------------------------------!
!    Read Input Wake                                                                            !
!-----------------------------------------------------------------------------------------------!
   INFILE='WAKE.INP'
   OPEN(UNIT=12,FILE=INFILE,STATUS='UNKNOWN')
   READ(12,*) NWR,NWT
   ALLOCATE(WRR(NWR),WTT(NWT),VX(NWT,NWR),VT(NWT,NWR),VR(NWT,NWR))
   WRR=0.D0
   WTT=0.D0
   VX =0.D0
   VT =0.D0
   VR =0.D0
   READ(12,*) (WRR(I),I=1,NWR)
   DO J=1,NWT
      READ(12,*) WTT(J),(VX(J,I),I=1,NWR)
   END DO !J=1,NWT
   READ(12,*)
   DO J=1,NWT
      READ(12,*) WTT(J),(VT(J,I),I=1,NWR)
   END DO !J=1,NWT
   READ(12,*)
   DO J=1,NWT
      READ(12,*) WTT(J),(VR(J,I),I=1,NWR)
   END DO !J=1,NWT
   WTT=WTT*PI/180.D0 !Convert to Radians
   CLOSE(UNIT=12)
!-----------------------------------------------------------------------------------------------!
!    Compute Fourier Analysis                                                                   !
!-----------------------------------------------------------------------------------------------!
   ALLOCATE(ANVX(0:IFREQ,NWR),BNVX(1:IFREQ,NWR))
   ALLOCATE(ANVT(0:IFREQ,NWR),BNVT(1:IFREQ,NWR))
   ALLOCATE(ANVR(0:IFREQ,NWR),BNVR(1:IFREQ,NWR))
   ANVX=0.D0
   BNVX=0.D0
   ANVT=0.D0
   BNVT=0.D0
   ANVR=0.D0
   BNVR=0.D0
!-----------------------------------------------------------------------------------------------!
   DO J=1,NWR
      CALL FOURIER_COEF(IFREQ,NWT,WTT,VX(:,J),ANVX(:,J),BNVX(:,J))
      CALL FOURIER_COEF(IFREQ,NWT,WTT,VT(:,J),ANVT(:,J),BNVT(:,J))
      CALL FOURIER_COEF(IFREQ,NWT,WTT,VR(:,J),ANVR(:,J),BNVR(:,J))
   END DO !J=1,NWR
!-----------------------------------------------------------------------------------------------!
CASE(13)
!-----------------------------------------------------------------------------------------------!
!    Read Field Points                                                                          !
!-----------------------------------------------------------------------------------------------!
   IF (IFIELD == 1) THEN
      INFILE='VFIELD.INP'
      OPEN(UNIT=13,FILE=INFILE,STATUS='UNKNOWN')
      READ(13,*)
      READ(13,*)
      READ(13,*) ZONE,TFILE,FPOINT,IFILE,NFX,JFILE,NFY
      ALLOCATE(XVF(NFX,NFY),YVF(NFX,NFY),ZVF(NFX,NFY))
      DO J=1,NFY
         DO I=1,NFX
            READ(13,*) XVF(I,J),YVF(I,J),ZVF(I,J)
         END DO !I=1,NFX
      END DO !J=1,NFY
      CLOSE(UNIT=13)
   END IF !(IFIELD == 1)
!-----------------------------------------------------------------------------------------------!
CASE(14)
!-----------------------------------------------------------------------------------------------!
!    Read Solution File                                                                         !
!-----------------------------------------------------------------------------------------------!
   INFILE='SOL.INP'
   OPEN(UNIT=14,FILE=INFILE,STATUS='UNKNOWN')
   READ(14,*)
   READ(14,*)
!-----------------------------------------------------------------------------------------------!
!    Blade Solution                                                                             !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      READ(14,*) ZONE,TFILEP,TZONE,FPOINT,IFILE,NCP,JFILE,NRP
      DO J=1,NRP
         DO I=1,NCP
            READ(14,*) WORK1,WORK1,WORK1,POTP(I,J,0),SOURCEP(I,J,1), & !XP0(I,J),YP0(I,J),ZP0(I,J)
                       VT1P(I,J),VT2P(I,J),VTT1P(I,J),VTT2P(I,J),CPP(I,J,0), &
                       CPNP(I,J,0)
         END DO !I=1,NCP
      END DO !J=1,NRP
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake Solution                                                                        !
!-----------------------------------------------------------------------------------------------!
   IF ((IP == 1).AND.(NCPW < 0)) THEN
      READ(14,*)
      DO J=1,NRW
         DO I=1,IABS(NCPW)
            READ(14,*) WORK1,WORK1,WORK1,POTPWP(I,J,0), & !XPW0(I,J),YPW0(I,J),ZPW0(I,J)
                       SOURCEPWCAV(I,J,1),VT1PW(I,J),WORK1,WORK1,WORK1,WORK1,WORK1
         END DO !I=1,IABS(NCPW)
      END DO !J=1,NRW
   END IF !((IP == 1).AND.(NCPW < 0))
   IF ((IP == 1).AND.(NCPW > 0)) THEN
      READ(14,*)
      DO J=1,NRW
         DO I=1,IABS(NCPW)
            READ(14,*) WORK1,WORK1,WORK1,POTPWS(I,J,0), & !XPW0(I,J),YPW0(I,J),ZPW0(I,J)
                       SOURCEPWCAV(I,J,1),VT1PW(I,J),WORK1,WORK1,WORK1,WORK1,WORK1
         END DO !I=1,IABS(NCPW)
      END DO !J=1,NRW
   END IF !((IP == 1).AND.(NCPW > 0))
!-----------------------------------------------------------------------------------------------!
!    Nozzle Solution                                                                            !
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      READ(14,*) ZONE,TFILEN,TZONE,FPOINT,IFILE,NNXT1,JFILE,NNTP
      DO J=1,NNTP
         DO I=1,NNXT1
            READ(14,*) WORK1,WORK1,WORK1,POTN(I,J,0),SOURCEN(I,J,1), & !XN0(I,J),YN0(I,J),ZN0(I,J)
                       VT1N(I,J),VT2N(I,J),VTT1N(I,J),VTT2N(I,J),CPN(I,J,0), &
                       CPNN(I,J,0)
         END DO !I=1,NNXT1
      END DO !J=1,NNTP
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub Solution                                                                               !
!-----------------------------------------------------------------------------------------------!
   IF (IABS(IH) == 1) THEN
      READ(14,*) ZONE,TFILEH,TZONE,FPOINT,IFILE,NHX,JFILE,NHTP
      DO J=1,NHTP
         DO I=1,NHX
            READ(14,*) WORK1,WORK1,WORK1,POTH(I,J,0),SOURCEH(I,J,1), & !XH0(I,J),YH0(I,J),ZH0(I,J)
                       VT1H(I,J),VT2H(I,J),VTT1H(I,J),VTT2H(I,J),CPH(I,J,0), &
                       CPNH(I,J,0)
         END DO !I=1,NHX
      END DO !J=1,NHTP
   END IF !(IABS(IH) == 1)
!-----------------------------------------------------------------------------------------------!
   CLOSE(UNIT=14)
!-----------------------------------------------------------------------------------------------!
!    Read Solution File                                                                         !
!-----------------------------------------------------------------------------------------------!
   INFILE='SOLW.INP'
   OPEN(UNIT=15,FILE=INFILE,STATUS='UNKNOWN')
   READ(15,*)
   READ(15,*)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      READ(15,*) ZONE,TFILE,WAKE,TZONE,IFILE,NRW,FPOINT
      DO J=1,NRW
         READ(15,*) WORK1,POTPW(1,J,0),WORK1,WORK1,WORK1,WORK1,DCP(J)
      END DO !J=1,NRW
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      READ(15,*) ZONE,TFILEN,WAKE,TZONE,IFILE,NNTP,FPOINT
      DO J=1,NNTP
         READ(15,*) WORK1,POTNW(1,J,0),WORK1,WORK1,WORK1,WORK1,DCN(J)
      END DO !J=1,NNTP
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
   CLOSE(UNIT=15)
!-----------------------------------------------------------------------------------------------!
END SELECT !(INPFILE)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE INPUT
!-----------------------------------------------------------------------------------------------!
