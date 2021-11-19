!-----------------------------------------------------------------------------------------------!
!    Generate Blade Grid                                                                        !
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
SUBROUTINE BLADEGRID
!-----------------------------------------------------------------------------------------------!
!    Modified  : 14102013, J. Baltazar, 2013 version 1.0                                        !
!    Modified  : 09052014, J. Baltazar, 2014 version 1.0                                        !
!    Modified  : 02122014, J. Baltazar, 2014 version 1.2                                        !
!    Modified  : 06122014, J. Baltazar, 2014 version 1.3                                        !
!    Modified  : 29062015, J. Baltazar, 2015 version 1.3                                        !
!    Modified  : 05072017, J. Baltazar, 2017 version 1.0                                        !
!    Modified  : 05012018, J. Baltazar, 2018 version 1.0                                        !
!    Modified  : 06042020, J. Baltazar, 2020 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPANEL_MOD
IMPLICIT NONE
INTEGER*4 :: I,IM1,J,K,IOERR
INTEGER*4 :: LIDENT,NCI1
DOUBLE PRECISION :: PHI,SPHI,CPHI,TK,X0,TTE,TLE,STE,FX
DOUBLE PRECISION :: TTETA,TETA,STETA,CTETA,SLE,RLE,TC,XC,RR,LIDENTN,PITCH
DOUBLE PRECISION :: CGAP,PTGAP,XRGAP,SKEWGAP
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: RI,CI,PTI,XRI,SKEWI,T0I,F0I
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: R,C,PT,XR,SKEW,T0,F0
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: S,SCN,SCII,TN0,TNF,TNN,FN0,FNN
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: SCI,YBI,YFI,SC,SB,YB,SF,YF,YBN,YFN,SC1,YB1,YF1
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: SBI,SFI,SBN,SFN,RROT,TROT
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: JUNK
!-----------------------------------------------------------------------------------------------!
TK=2.5980762D0
X0=1.D0
!-----------------------------------------------------------------------------------------------!
!    Read Blade Data for Conventional Grid                                                      !
!-----------------------------------------------------------------------------------------------!
I=1
DO WHILE (IDENTP(I:I).NE.' ')
   I=I+1
END DO !(IDENTP(I:I).NE.' ')
LIDENT=I-1
!-----------------------------------------------------------------------------------------------!
OPEN(UNIT=30,FILE=IDENTP(1:LIDENT)//'.DAT',STATUS='UNKNOWN')
READ(30,*) IGEOM
READ(30,*)
READ(30,*)
READ(30,*)
READ(30,*) NRI,NCI
!-----------------------------------------------------------------------------------------------!
ALLOCATE(RI(NRI),CI(NRI),PTI(NRI),XRI(NRI),SKEWI(NRI),T0I(NRI),F0I(NRI))
RI   =0.D0
CI   =0.D0
PTI  =0.D0
XRI  =0.D0
SKEWI=0.D0
T0I  =0.D0
F0I  =0.D0
!-----------------------------------------------------------------------------------------------!
DO I=1,NRI
   READ(30,*) RI(I),CI(I),PTI(I),XRI(I),SKEWI(I),T0I(I),F0I(I)
END DO !I=1,NRI
!-----------------------------------------------------------------------------------------------!
!    Test File Length                                                                           !
!-----------------------------------------------------------------------------------------------!
ALLOCATE(JUNK(4*NRI*NCI))
READ(30,*,IOSTAT=IOERR) (JUNK(K),K=1,(4*NRI*NCI))
!-----------------------------------------------------------------------------------------------!
IF (IOERR == 0) THEN
   ALLOCATE(SBI(NCI,NRI),YBI(NCI,NRI),SFI(NCI,NRI),YFI(NCI,NRI))
   SBI=0.D0
   YBI=0.D0
   SFI=0.D0
   YFI=0.D0
   K=1
   DO J=1,NRI
      DO I=1,NCI
         SBI(I,J)=JUNK(K  )
         YBI(I,J)=JUNK(K+1)
         SFI(I,J)=JUNK(K+2)
         YFI(I,J)=JUNK(K+3)
         K=K+4
      END DO !I=1,NCI
   END DO !J=1,NRI
END IF !(IOERR == 0)
!-----------------------------------------------------------------------------------------------!
IF (IOERR < 0) THEN
   ALLOCATE(SCI(NCI,NRI),YBI(NCI,NRI),YFI(NCI,NRI))
   SCI=0.D0
   YBI=0.D0
   YFI=0.D0
   K=1
   DO J=1,NRI
      DO I=1,NCI
         SCI(I,J)=JUNK(K  )
         YBI(I,J)=JUNK(K+1)
         YFI(I,J)=JUNK(K+2)
         K=K+3
      END DO !I=1,NCI
   END DO !J=1,NRI
END IF !(IOERR < 0)
!-----------------------------------------------------------------------------------------------!
DEALLOCATE(JUNK)
CLOSE(UNIT=30)
!-----------------------------------------------------------------------------------------------!
!    Counters                                                                                   !
!-----------------------------------------------------------------------------------------------!
NC1 =NC+1
NCP =2*NC
NCP1=NCP+1
NRP1=NRP+1
!-----------------------------------------------------------------------------------------------!
IF (RMAX == 0.D0) RMAX=1.D0
ALPHAH =ALPHAH *PI/180.D0
ALPHAT =ALPHAT *PI/180.D0
ALPHALE=ALPHALE*PI/180.D0
ALPHATE=ALPHATE*PI/180.D0
!-----------------------------------------------------------------------------------------------!
ALLOCATE(R(NRP1),C(NRP1),PT(NRP1),XR(NRP1),SKEW(NRP1),T0(NRP1),F0(NRP1))
R   =0.D0
C   =0.D0
PT  =0.D0
XR  =0.D0
SKEW=0.D0
T0  =0.D0
F0  =0.D0
!-----------------------------------------------------------------------------------------------!
DO I=1,NRP1
   PHI=DFLOAT(I-1)*(PI-ALPHAT-ALPHAH)/DFLOAT(NRP)+ALPHAH
   R(I)=((RPH*DCOS(ALPHAT)+RMAX*DCOS(ALPHAH))-(RMAX-RPH)*DCOS(PHI))/(DCOS(ALPHAH)+DCOS(ALPHAT))
END DO !I=1,NRP1
!-----------------------------------------------------------------------------------------------!
IF (INTERP == 0) THEN
   CALL LININT(NRI,RI,CI   ,NRP1,R,C   )
   CALL LININT(NRI,RI,PTI  ,NRP1,R,PT  )
   CALL LININT(NRI,RI,XRI  ,NRP1,R,XR  )
   CALL LININT(NRI,RI,SKEWI,NRP1,R,SKEW)
   CALL LININT(NRI,RI,T0I  ,NRP1,R,T0  )
   CALL LININT(NRI,RI,F0I  ,NRP1,R,F0  )
ELSEIF (INTERP == 1) THEN
   CALL INTK1 (NRI,RI,CI   ,NRP1,R,C   )
   CALL INTK1 (NRI,RI,PTI  ,NRP1,R,PT  )
   CALL INTK1 (NRI,RI,XRI  ,NRP1,R,XR  )
   CALL INTK1 (NRI,RI,SKEWI,NRP1,R,SKEW)
   CALL INTK1 (NRI,RI,T0I  ,NRP1,R,T0  )
   CALL INTK1 (NRI,RI,F0I  ,NRP1,R,F0  )
ELSEIF (INTERP == 2) THEN
   CALL SPLINT(NRI,RI,CI   ,NRP1,R,C   )
   CALL SPLINT(NRI,RI,PTI  ,NRP1,R,PT  )
   CALL SPLINT(NRI,RI,XRI  ,NRP1,R,XR  )
   CALL SPLINT(NRI,RI,SKEWI,NRP1,R,SKEW)
   CALL SPLINT(NRI,RI,T0I  ,NRP1,R,T0  )
   CALL SPLINT(NRI,RI,F0I  ,NRP1,R,F0  )
END IF !(INTERP)
!-----------------------------------------------------------------------------------------------!
IF (IOERR == 0) THEN
   ALLOCATE(SB(NCI,NRP1),YB(NCI,NRP1),SF(NCI,NRP1),YF(NCI,NRP1))
   SB=0.D0
   YB=0.D0
   SF=0.D0
   YF=0.D0
   DO I=1,NCI
      IF (INTERP == 0) THEN
         CALL LININT(NRI,RI,SBI(I,:),NRP1,R,SB(I,:))
         CALL LININT(NRI,RI,YBI(I,:),NRP1,R,YB(I,:))
         CALL LININT(NRI,RI,SFI(I,:),NRP1,R,SF(I,:))
         CALL LININT(NRI,RI,YFI(I,:),NRP1,R,YF(I,:))
      ELSEIF (INTERP == 1) THEN
         CALL INTK1 (NRI,RI,SBI(I,:),NRP1,R,SB(I,:))
         CALL INTK1 (NRI,RI,YBI(I,:),NRP1,R,YB(I,:))
         CALL INTK1 (NRI,RI,SFI(I,:),NRP1,R,SF(I,:))
         CALL INTK1 (NRI,RI,YFI(I,:),NRP1,R,YF(I,:))
      ELSEIF (INTERP == 2) THEN
         CALL SPLINT(NRI,RI,SBI(I,:),NRP1,R,SB(I,:))
         CALL SPLINT(NRI,RI,YBI(I,:),NRP1,R,YB(I,:))
         CALL SPLINT(NRI,RI,SFI(I,:),NRP1,R,SF(I,:))
         CALL SPLINT(NRI,RI,YFI(I,:),NRP1,R,YF(I,:))
      END IF !(INTERP)
   END DO !I=1,NCI
   ALLOCATE(SC1(NCI,NRP1),YB1(NCI,NRP1),YF1(NCI,NRP1))
   SC1=0.D0
   YB1=0.D0
   YF1=0.D0
   SC1(:,:)=0.5D0*(SB(:,:)+SF(:,:))
   DO J=1,NRP1
      DO I=1,NCI
         IF (INTERP == 0) THEN
            CALL LININT(NCI,SC1(:,J),YB(:,J),1,SC1(I,J),YB1(I,J))
            CALL LININT(NCI,SC1(:,J),YF(:,J),1,SC1(I,J),YF1(I,J))
         ELSEIF (INTERP == 1) THEN
            CALL INTK1 (NCI,SC1(:,J),YB(:,J),1,SC1(I,J),YB1(I,J))
            CALL INTK1 (NCI,SC1(:,J),YF(:,J),1,SC1(I,J),YF1(I,J))
         ELSEIF (INTERP == 2) THEN
            CALL SPLINT(NCI,SC1(:,J),YB(:,J),1,SC1(I,J),YB1(I,J))
            CALL SPLINT(NCI,SC1(:,J),YF(:,J),1,SC1(I,J),YF1(I,J))
         END IF !(INTERP)
      END DO !I=1,NCI
   END DO !J=1,NRP1
END IF !(IOERR == 0)
!-----------------------------------------------------------------------------------------------!
IF (IOERR < 0) THEN
   ALLOCATE(SC(NCI,NRP1),YB(NCI,NRP1),YF(NCI,NRP1))
   SC=0.D0
   YB=0.D0
   YF=0.D0
   DO I=1,NCI
      IF (INTERP == 0) THEN
         CALL LININT(NRI,RI,SCI(I,:),NRP1,R,SC(I,:))
         CALL LININT(NRI,RI,YBI(I,:),NRP1,R,YB(I,:))
         CALL LININT(NRI,RI,YFI(I,:),NRP1,R,YF(I,:))
      ELSEIF (INTERP == 1) THEN
         CALL INTK1 (NRI,RI,SCI(I,:),NRP1,R,SC(I,:))
         CALL INTK1 (NRI,RI,YBI(I,:),NRP1,R,YB(I,:))
         CALL INTK1 (NRI,RI,YFI(I,:),NRP1,R,YF(I,:))
      ELSEIF (INTERP == 2) THEN
         CALL SPLINT(NRI,RI,SCI(I,:),NRP1,R,SC(I,:))
         CALL SPLINT(NRI,RI,YBI(I,:),NRP1,R,YB(I,:))
         CALL SPLINT(NRI,RI,YFI(I,:),NRP1,R,YF(I,:))
      END IF !(INTERP)
   END DO !I=1,NCI
END IF !(IOERR < 0)
!-----------------------------------------------------------------------------------------------!
!    Define Chordwise Positions and Interpolate Section                                         !
!-----------------------------------------------------------------------------------------------!
ALLOCATE(SCN(NC1),YBN(NC1,NRP1),YFN(NC1,NRP1))
SCN=0.D0
YBN=0.D0
YFN=0.D0
IF (ISC == 1) THEN
   ALLOCATE(S(NCI),SBN(NC1,NRP1),SFN(NC1,NRP1))
   S  =0.D0
   SBN=0.D0
   SFN=0.D0
END IF !(ISC == 1)
!-----------------------------------------------------------------------------------------------!
DO I=1,NC1
   PHI=DFLOAT(I-1)*(PI-ALPHATE-ALPHALE)/DFLOAT(NC)+ALPHALE
   SCN(I)=(DCOS(ALPHALE)-DCOS(PHI))/(DCOS(ALPHALE)+DCOS(ALPHATE))
END DO !I=1,NC1
IF (DABS(SCN(1)) < TOL)  SCN(1)=0.D0
!-----------------------------------------------------------------------------------------------!
DO J=1,NRP1
   TTE=DABS(YB(NCI,J)-YF(NCI,J))
   IF (TTE < TOL) TTE=0.D0
!-----------------------------------------------------------------------------------------------!
   TLE=DABS(YB(1,J)-YF(1,J))
   IF (TLE < TOL) TLE=0.D0
!-----------------------------------------------------------------------------------------------!
!    Interpolate Section                                                                        !
!-----------------------------------------------------------------------------------------------!
   IF ((TLE /= 0.D0).AND.(TTE /= 0.D0)) THEN
      IF (IOERR == 0) STOP 'LE AND TE CORRECTION NOT AVAILABLE FOR CURRENT GEOMETRY DEFINITION!'
      NCI1=NCI-1
!-----------------------------------------------------------------------------------------------!
!    Trailing Edge                                                                              !
!-----------------------------------------------------------------------------------------------!
      ALLOCATE(SCII(NCI1),TN0(NCI),TNF(NCI1),TNN(NC1),FN0(NCI),FNN(NC1))
      SCII=0.D0
      TN0 =0.D0
      TNF =0.D0
      TNN =0.D0
      FN0 =0.D0
      FNN =0.D0
!-----------------------------------------------------------------------------------------------!
      TN0(:)=YB(:,J)-YF(:,J)
      FN0(:)=0.5D0*(YB(:,J)+YF(:,J))
!-----------------------------------------------------------------------------------------------!
      STE=SC(NCI,J)+(SC(NCI,J)-SC(NCI1,J))*TN0(NCI)/(TN0(NCI1)-TN0(NCI))
      DO I=2,NCI1
         SCII(I-1)=SC(I,J)*SC(NCI,J)/STE
         FX=TK*DSQRT(SCII(I-1))*(X0-SCII(I-1))
         TNF(I-1)=TN0(I)/FX
      END DO !I=2,NCI1
      SCII(NCI1)=SC(NCI,J)*SC(NCI,J)/STE
      FX=TK*DSQRT(SCII(NCI1))*(X0-SCII(NCI1))
      TNF(NCI1)=TN0(NCI)/FX
!-----------------------------------------------------------------------------------------------!
      IF (INTERP == 0) THEN
         CALL LININT(NCI1,SCII   ,TNF,NC1,SCN,TNN)
         CALL LININT(NCI ,SC(:,J),FN0,NC1,SCN,FNN)
      ELSEIF (INTERP == 1) THEN
         CALL INTK1 (NCI1,SCII   ,TNF,NC1,SCN,TNN)
         CALL INTK1 (NCI ,SC(:,J),FN0,NC1,SCN,FNN)
      ELSEIF (INTERP == 2) THEN
         CALL SPLINT(NCI1,SCII   ,TNF,NC1,SCN,TNN)
         CALL SPLINT(NCI ,SC(:,J),FN0,NC1,SCN,FNN)
      END IF !(INTERP)
!-----------------------------------------------------------------------------------------------!
      TNN(:)=TNN(:)*TK*DSQRT(SCN(:))*(X0-SCN(:))
!-----------------------------------------------------------------------------------------------!
!    Leading Edge                                                                               !
!-----------------------------------------------------------------------------------------------!
      TTETA=(TN0(2)-TN0(1))/(SC(2,J)-SC(1,J))
      TETA =DATAN(TTETA)
      STETA=DSIN(TETA)
      CTETA=DCOS(TETA)
      SLE  =TN0(1)*(1.D0-STETA)/(CTETA+TTETA**2*CTETA-TTETA)
      RLE  =(TN0(1)+SLE*TTETA)/CTETA
!-----------------------------------------------------------------------------------------------!
      DO I=1,NC1
         IF (SCN(I) < SLE) THEN
            TNN(I)=DSQRT(RLE**2-(SCN(I)-RLE)**2)
         END IF !(SCN(I) < SLE)
      END DO !I=1,NC1
!-----------------------------------------------------------------------------------------------!
      YBN(:,J)=FNN(:)+0.5D0*TNN(:)
      YFN(:,J)=FNN(:)-0.5D0*TNN(:)
!-----------------------------------------------------------------------------------------------!
      DEALLOCATE(SCII,TN0,TNF,TNN,FN0,FNN)
!-----------------------------------------------------------------------------------------------!
!    Leading Edge Roundoff by a Circle                                                          !
!-----------------------------------------------------------------------------------------------!
   ELSEIF (TLE /= 0.D0) THEN
      IF (IOERR == 0) STOP 'LE CORRECTION NOT AVAILABLE FOR CURRENT GEOMETRY DEFINITION!'
      ALLOCATE(TN0(NCI),TNN(NC1),FN0(NCI),FNN(NC1))
      TN0=0.D0
      TNN=0.D0
      FN0=0.D0
      FNN=0.D0
!-----------------------------------------------------------------------------------------------!
      TN0(:)=YB(:,J)-YF(:,J)
      FN0(:)=0.5D0*(YB(:,J)+YF(:,J))
!-----------------------------------------------------------------------------------------------!
      TTETA=(TN0(2)-TN0(1))/(SC(2,J)-SC(1,J))
      TETA =DATAN(TTETA)
      STETA=DSIN(TETA)
      CTETA=DCOS(TETA)
      SLE  =TN0(1)*(1.D0-STETA)/(CTETA+TTETA**2*CTETA-TTETA)
      RLE  =(TN0(1)+SLE*TTETA)/CTETA
!-----------------------------------------------------------------------------------------------!
      DO I=1,NC1
         IF (SCN(I) < SLE) THEN
            TNN(I)=DSQRT(RLE**2-(SCN(I)-RLE)**2)
         ELSE !(SCN(I) < SLE)
            IF (INTERP == 0) THEN
               CALL LININT(NCI,SC(:,J),TN0,1,SCN(I),TNN(I))
            ELSEIF (INTERP == 1) THEN
               CALL INTK1 (NCI,SC(:,J),TN0,1,SCN(I),TNN(I))
            ELSEIF (INTERP == 2) THEN
               CALL SPLINT(NCI,SC(:,J),TN0,1,SCN(I),TNN(I))
            END IF !(INTERP)
         END IF !(SCN(I) < SLE)
      END DO !I=1,NC1
!-----------------------------------------------------------------------------------------------!
      IF (INTERP == 0) THEN
         CALL LININT(NCI,SC(:,J),FN0,NC1,SCN,FNN)
      ELSEIF (INTERP == 1) THEN
         CALL INTK1 (NCI,SC(:,J),FN0,NC1,SCN,FNN)
      ELSEIF (INTERP == 2) THEN
         CALL SPLINT(NCI,SC(:,J),FN0,NC1,SCN,FNN)
      END IF !(INTERP)
!-----------------------------------------------------------------------------------------------!
      YBN(:,J)=FNN(:)+0.5D0*TNN(:)
      YFN(:,J)=FNN(:)-0.5D0*TNN(:)
!-----------------------------------------------------------------------------------------------!
      DEALLOCATE(TN0,TNN,FN0,FNN)
!-----------------------------------------------------------------------------------------------!
!    Linear Extrapolation to the Trailing Edge                                                  !
!-----------------------------------------------------------------------------------------------!
   ELSEIF (TTE /= 0.D0) THEN
      !IF (IOERR == 0) STOP 'TE CORRECTION NOT AVAILABLE FOR CURRENT GEOMETRY DEFINITION!'
      IF (IOERR == 0) THEN
         NCI1=NCI-1
!-----------------------------------------------------------------------------------------------!
         ALLOCATE(SCII(NCI1),TN0(NCI),TNF(NCI1),TNN(NC1),FN0(NCI),FNN(NC1))
         SCII=0.D0
         TN0 =0.D0
         TNF =0.D0
         TNN =0.D0
         FN0 =0.D0
         FNN =0.D0
!-----------------------------------------------------------------------------------------------!
         TN0(:)=YB1(:,J)-YF1(:,J)
         FN0(:)=0.5D0*(YB1(:,J)+YF1(:,J))
!-----------------------------------------------------------------------------------------------!
         STE=SC1(NCI,J)+(SC1(NCI,J)-SC1(NCI1,J))*TN0(NCI)/(TN0(NCI1)-TN0(NCI))
         DO I=2,NCI1
            SCII(I-1)=SC1(I,J)*SC1(NCI,J)/STE
            FX=TK*DSQRT(SCII(I-1))*(X0-SCII(I-1))
            TNF(I-1)=TN0(I)/FX
         END DO !I=2,NCI1
         SCII(NCI1)=SC1(NCI,J)*SC1(NCI,J)/STE
         FX=TK*DSQRT(SCII(NCI1))*(X0-SCII(NCI1))
         TNF(NCI1)=TN0(NCI)/FX
!-----------------------------------------------------------------------------------------------!
         IF (INTERP == 0) THEN
            CALL LININT(NCI1,SCII    ,TNF,NC1,SCN,TNN)
            CALL LININT(NCI ,SC1(:,J),FN0,NC1,SCN,FNN)
         ELSEIF (INTERP == 1) THEN
            CALL INTK1 (NCI1,SCII    ,TNF,NC1,SCN,TNN)
            CALL INTK1 (NCI ,SC1(:,J),FN0,NC1,SCN,FNN)
         ELSEIF (INTERP == 2) THEN
            CALL SPLINT(NCI1,SCII    ,TNF,NC1,SCN,TNN)
            CALL SPLINT(NCI ,SC1(:,J),FN0,NC1,SCN,FNN)
         END IF !(INTERP)
!-----------------------------------------------------------------------------------------------!
         TNN(:)  =TNN(:)*TK*DSQRT(SCN(:))*(X0-SCN(:))
         YBN(:,J)=FNN(:)+0.5D0*TNN(:)
         YFN(:,J)=FNN(:)-0.5D0*TNN(:)
!-----------------------------------------------------------------------------------------------!
         DEALLOCATE(SCII,TN0,TNF,TNN,FN0,FNN)
      ELSE
         NCI1=NCI-1
!-----------------------------------------------------------------------------------------------!
         ALLOCATE(SCII(NCI1),TN0(NCI),TNF(NCI1),TNN(NC1),FN0(NCI),FNN(NC1))
         SCII=0.D0
         TN0 =0.D0
         TNF =0.D0
         TNN =0.D0
         FN0 =0.D0
         FNN =0.D0
!-----------------------------------------------------------------------------------------------!
         TN0(:)=YB(:,J)-YF(:,J)
         FN0(:)=0.5D0*(YB(:,J)+YF(:,J))
!-----------------------------------------------------------------------------------------------!
         STE=SC(NCI,J)+(SC(NCI,J)-SC(NCI1,J))*TN0(NCI)/(TN0(NCI1)-TN0(NCI))
         DO I=2,NCI1
            SCII(I-1)=SC(I,J)*SC(NCI,J)/STE
            FX=TK*DSQRT(SCII(I-1))*(X0-SCII(I-1))
            TNF(I-1)=TN0(I)/FX
         END DO !I=2,NCI1
         SCII(NCI1)=SC(NCI,J)*SC(NCI,J)/STE
         FX=TK*DSQRT(SCII(NCI1))*(X0-SCII(NCI1))
         TNF(NCI1)=TN0(NCI)/FX
!-----------------------------------------------------------------------------------------------!
         IF (INTERP == 0) THEN
            CALL LININT(NCI1,SCII   ,TNF,NC1,SCN,TNN)
            CALL LININT(NCI ,SC(:,J),FN0,NC1,SCN,FNN)
         ELSEIF (INTERP == 1) THEN
            CALL INTK1 (NCI1,SCII   ,TNF,NC1,SCN,TNN)
            CALL INTK1 (NCI ,SC(:,J),FN0,NC1,SCN,FNN)
         ELSEIF (INTERP == 2) THEN
            CALL SPLINT(NCI1,SCII   ,TNF,NC1,SCN,TNN)
            CALL SPLINT(NCI ,SC(:,J),FN0,NC1,SCN,FNN)
         END IF !(INTERP)
!-----------------------------------------------------------------------------------------------!
         TNN(:)  =TNN(:)*TK*DSQRT(SCN(:))*(X0-SCN(:))
         YBN(:,J)=FNN(:)+0.5D0*TNN(:)
         YFN(:,J)=FNN(:)-0.5D0*TNN(:)
!-----------------------------------------------------------------------------------------------!
         DEALLOCATE(SCII,TN0,TNF,TNN,FN0,FNN)
      END IF !(IOERR == 0)
!-----------------------------------------------------------------------------------------------!
!    Section Interpolation along the chord                                                      !
!-----------------------------------------------------------------------------------------------!
   ELSE
      IF (IOERR == 0) THEN
         IF (ISC == 0) THEN
            IF (INTERP == 0) THEN
               CALL LININT(NCI,SB(:,J),YB(:,J),NC1,SCN,YBN(:,J))
               CALL LININT(NCI,SF(:,J),YF(:,J),NC1,SCN,YFN(:,J))
            ELSEIF (INTERP == 1) THEN
               CALL INTK1 (NCI,SB(:,J),YB(:,J),NC1,SCN,YBN(:,J))
               CALL INTK1 (NCI,SF(:,J),YF(:,J),NC1,SCN,YFN(:,J))
            ELSEIF (INTERP == 2) THEN
               CALL SPLINT(NCI,SB(:,J),YB(:,J),NC1,SCN,YBN(:,J))
               CALL SPLINT(NCI,SF(:,J),YF(:,J),NC1,SCN,YFN(:,J))
            END IF !(INTERP)
         ELSEIF (ISC == 1) THEN
            S=0.D0
            DO I=2,NCI
               IM1=I-1
               S(I)=S(IM1)+DSQRT((SB(I,J)-SB(IM1,J))**2+(YB(I,J)-YB(IM1,J))**2)
            END DO !I=2,NCI
            S(:)=S(:)/S(NCI)
            IF (INTERP == 0) THEN
               CALL LININT(NCI,S,SB(:,J),NC1,SCN,SBN(:,J))
               CALL LININT(NCI,S,YB(:,J),NC1,SCN,YBN(:,J))
            ELSEIF (INTERP == 1) THEN
               CALL INTK1 (NCI,S,SB(:,J),NC1,SCN,SBN(:,J))
               CALL INTK1 (NCI,S,YB(:,J),NC1,SCN,YBN(:,J))
            ELSEIF (INTERP == 2) THEN
               CALL SPLINT(NCI,S,SB(:,J),NC1,SCN,SBN(:,J))
               CALL SPLINT(NCI,S,YB(:,J),NC1,SCN,YBN(:,J))
            END IF !(INTERP)
!-----------------------------------------------------------------------------------------------!
            DO I=2,NCI
               IM1=I-1
               S(I)=S(IM1)+DSQRT((SF(I,J)-SF(IM1,J))**2+(YF(I,J)-YF(IM1,J))**2)
            END DO !I=2,NCI
            S(:)=S(:)/S(NCI)
            IF (INTERP == 0) THEN
               CALL LININT(NCI,S,SF(:,J),NC1,SCN,SFN(:,J))
               CALL LININT(NCI,S,YF(:,J),NC1,SCN,YFN(:,J))
            ELSEIF (INTERP == 1) THEN
               CALL INTK1 (NCI,S,SF(:,J),NC1,SCN,SFN(:,J))
               CALL INTK1 (NCI,S,YF(:,J),NC1,SCN,YFN(:,J))
            ELSEIF (INTERP == 2) THEN
               CALL SPLINT(NCI,S,SF(:,J),NC1,SCN,SFN(:,J))
               CALL SPLINT(NCI,S,YF(:,J),NC1,SCN,YFN(:,J))
            END IF !(INTERP)
         END IF !(ISC)
      END IF !(IOERR == 0)
!-----------------------------------------------------------------------------------------------!
      IF (IOERR < 0) THEN
         IF (ISC == 0) THEN
            IF (INTERP == 0) THEN
               CALL LININT(NCI,SC(:,J),YB(:,J),NC1,SCN,YBN(:,J))
               CALL LININT(NCI,SC(:,J),YF(:,J),NC1,SCN,YFN(:,J))
            ELSEIF (INTERP == 1) THEN
               CALL INTK1 (NCI,SC(:,J),YB(:,J),NC1,SCN,YBN(:,J))
               CALL INTK1 (NCI,SC(:,J),YF(:,J),NC1,SCN,YFN(:,J))
            ELSEIF (INTERP == 2) THEN
               CALL SPLINT(NCI,SC(:,J),YB(:,J),NC1,SCN,YBN(:,J))
               CALL SPLINT(NCI,SC(:,J),YF(:,J),NC1,SCN,YFN(:,J))
            END IF !(INTERP)
         ELSEIF (ISC == 1) THEN
            STOP 'NOT DEFINED FOR ISC=1!'
         END IF !(ISC)
      END IF !(IOERR < 0)
   END IF
END DO !J=1,NRP1
!-----------------------------------------------------------------------------------------------!
!    Non-Dimensional Quantities by Propeller Radius                                             !
!-----------------------------------------------------------------------------------------------!
C   =2.D0*C
PT  =2.D0*PT
XR  =2.D0*XR
SKEW=SKEW*PI/180.D0
T0  =T0*C
F0  =F0*C
IF (IGEOM == 'WINGGEOM') PT=0.5D0*PT !Wing Case
!-----------------------------------------------------------------------------------------------!
!    Calculate Panel Corner Points                                                              !
!-----------------------------------------------------------------------------------------------!
IF (ISTRIP == 0) THEN
   ALLOCATE(XP(NCP1,NRP1  ),YP(NCP1,NRP1  ),ZP(NCP1,NRP1  ),RP(NCP1,NRP1  ),TP(NCP1,NRP1  ))
ELSEIF ((ISTRIP == 1).OR.(ISTRIP == 2)) THEN
   ALLOCATE(XP(NCP1,NRP1+1),YP(NCP1,NRP1+1),ZP(NCP1,NRP1+1),RP(NCP1,NRP1+1),TP(NCP1,NRP1+1))
END IF !(ISTRIP)
XP=0.D0
YP=0.D0
ZP=0.D0
RP=0.D0
TP=0.D0
!-----------------------------------------------------------------------------------------------!
DO J=1,NRP1
!-----------------------------------------------------------------------------------------------!
!    Pitch Angle of Blade Section                                                               !
!-----------------------------------------------------------------------------------------------!
   PHI =DATAN2(PT(J),2.D0*PI*R(J)) !DATAN(PT(J)/(2.D0*PI*R(J)))
   SPHI=DSIN(PHI)
   CPHI=DCOS(PHI)
!-----------------------------------------------------------------------------------------------!
!    Loop on Chordwise Stations                                                                 !
!-----------------------------------------------------------------------------------------------!
   DO I=1,NC1
      IF (ISC == 0) THEN
!-----------------------------------------------------------------------------------------------!
!    Wing Case                                                                                  !
!-----------------------------------------------------------------------------------------------!
         IF (IGEOM == 'WINGGEOM') THEN
            XC=XR(J)+C(J)*(SCN(I)-0.5D0)*DCOSD(PT(J))
            TC=     -C(J)*(SCN(I)-0.5D0)*DSIND(PT(J))
            XP(NC1-I+1,J)=XC+YFN(I,J)*C(J)*DSIND(PT(J))
            XP(NC1+I-1,J)=XC+YBN(I,J)*C(J)*DSIND(PT(J))
            YP(NC1-I+1,J)=R(J)
            YP(NC1+I-1,J)=R(J)
            ZP(NC1-I+1,J)=TC+YFN(I,J)*C(J)*DCOSD(PT(J))
            ZP(NC1+I-1,J)=TC+YBN(I,J)*C(J)*DCOSD(PT(J))
!-----------------------------------------------------------------------------------------------!
!    Rotor Case                                                                                 !
!-----------------------------------------------------------------------------------------------!
         ELSE !(IGEOM)
            XC=  XR(J)+C(J)*(SCN(I)-0.5D0)*SPHI
            TC=SKEW(J)+C(J)*(SCN(I)-0.5D0)*CPHI/R(J)
            XP(NC1-I+1,J)=XC-YFN(I,J)*C(J)*CPHI
            XP(NC1+I-1,J)=XC-YBN(I,J)*C(J)*CPHI
            TP(NC1-I+1,J)=TC+YFN(I,J)*C(J)*SPHI/R(J)
            TP(NC1+I-1,J)=TC+YBN(I,J)*C(J)*SPHI/R(J)
            YP(NC1-I+1,J)=R(J)*DCOS(TP(NC1-I+1,J))
            YP(NC1+I-1,J)=R(J)*DCOS(TP(NC1+I-1,J))
            ZP(NC1-I+1,J)=R(J)*DSIN(TP(NC1-I+1,J))
            ZP(NC1+I-1,J)=R(J)*DSIN(TP(NC1+I-1,J))
         END IF !(IGEOM)
      ELSEIF (ISC == 1) THEN
!-----------------------------------------------------------------------------------------------!
!    Wing Case                                                                                  !
!-----------------------------------------------------------------------------------------------!
         IF (IGEOM == 'WINGGEOM') THEN
            XP(NC1-I+1,J)=XR(J)+C(J)*(SFN(I,J)-0.5D0)*DCOSD(PT(J))+YFN(I,J)*C(J)*DSIND(PT(J))
            XP(NC1+I-1,J)=XR(J)+C(J)*(SBN(I,J)-0.5D0)*DCOSD(PT(J))+YBN(I,J)*C(J)*DSIND(PT(J))
            YP(NC1-I+1,J)=R(J)
            YP(NC1+I-1,J)=R(J)
            ZP(NC1-I+1,J)=-C(J)*(SFN(I,J)-0.5D0)*DSIND(PT(J))+YFN(I,J)*C(J)*DCOSD(PT(J))
            ZP(NC1+I-1,J)=-C(J)*(SBN(I,J)-0.5D0)*DSIND(PT(J))+YBN(I,J)*C(J)*DCOSD(PT(J))
!-----------------------------------------------------------------------------------------------!
!    Rotor Case                                                                                 !
!-----------------------------------------------------------------------------------------------!
         ELSE !(IGEOM)
            XP(NC1-I+1,J)=  XR(J)+C(J)*(SFN(I,J)-0.5D0)*SPHI-YFN(I,J)*C(J)*CPHI
            XP(NC1+I-1,J)=  XR(J)+C(J)*(SBN(I,J)-0.5D0)*SPHI-YBN(I,J)*C(J)*CPHI
            TP(NC1-I+1,J)=SKEW(J)+C(J)*(SFN(I,J)-0.5D0)*CPHI/R(J)+YFN(I,J)*C(J)*SPHI/R(J)
            TP(NC1+I-1,J)=SKEW(J)+C(J)*(SBN(I,J)-0.5D0)*CPHI/R(J)+YBN(I,J)*C(J)*SPHI/R(J)
            YP(NC1-I+1,J)=   R(J)*DCOS(TP(NC1-I+1,J))
            YP(NC1+I-1,J)=   R(J)*DCOS(TP(NC1+I-1,J))
            ZP(NC1-I+1,J)=   R(J)*DSIN(TP(NC1-I+1,J))
            ZP(NC1+I-1,J)=   R(J)*DSIN(TP(NC1+I-1,J))
         END IF !(IGEOM)
      END IF !(ISC)
   END DO !I=1,NC1
!-----------------------------------------------------------------------------------------------!
!    Hub Matching With Blade Root                                                               !
!-----------------------------------------------------------------------------------------------!
   IF ((IHCORR == 1).AND.(J == 1)) THEN
      DO I=1,NCP1
         IF (INTERH == 0) THEN
            CALL LININT(NHI,XHI,RHI,1,XP(I,1),RP(I,1))
         ELSEIF (INTERH == 1) THEN
            CALL INTK1 (NHI,XHI,RHI,1,XP(I,1),RP(I,1))
         ELSEIF (INTERH == 2) THEN
            CALL SPLINT(NHI,XHI,RHI,1,XP(I,1),RP(I,1))
         END IF !(INTERH)
         YP(I,1)=RP(I,1)*DCOS(TP(I,1))
         ZP(I,1)=RP(I,1)*DSIN(TP(I,1))
      END DO !I=1,NCP1
   END IF !((IHCORR == 1).AND.(J == 1))
END DO !J=1,NRP1
!-----------------------------------------------------------------------------------------------!
!    Controlable Blade Pitch                                                                    !
!-----------------------------------------------------------------------------------------------!
IF (ANGPITCH /= 0.D0) THEN
   IF (ISTRIP == 0) THEN
      ALLOCATE(RROT(NCP1,NRP1  ),TROT(NCP1,NRP1  ))
   ELSEIF ((ISTRIP == 1).OR.(ISTRIP == 2)) THEN
      ALLOCATE(RROT(NCP1,NRP1+1),TROT(NCP1,NRP1+1))
   END IF !(ISTRIP)
!-----------------------------------------------------------------------------------------------!
!    Rotation Along Y-Axis                                                                      !
!-----------------------------------------------------------------------------------------------!
   RROT=DSQRT(XP*XP+ZP*ZP)
   TROT=DATAN2(XP,ZP)
   TROT=TROT+ANGPITCH*PI/180.D0
   XP=RROT*DSIN(TROT)
   ZP=RROT*DCOS(TROT)
!-----------------------------------------------------------------------------------------------!
!    Original Cylindrical Coordinates Definition                                                !
!-----------------------------------------------------------------------------------------------!
   RP=DSQRT(YP*YP*ZP*ZP)
   TP=DATAN2(ZP,YP)
!-----------------------------------------------------------------------------------------------!
!    Root Section                                                                               !
!-----------------------------------------------------------------------------------------------!
   RP(:,1)=RPH
   YP(:,1)=RP(:,1)*DCOS(TP(:,1))
   ZP(:,1)=RP(:,1)*DSIN(TP(:,1))
   DEALLOCATE(RROT,TROT)
END IF !(ANGPITCH /= 0.D0)
!-----------------------------------------------------------------------------------------------!
!    Gap Model                                                                                  !
!-----------------------------------------------------------------------------------------------!
IF ((IGRIDI == 0).AND.(ISTRIP == 1)) STOP 'IGRIDI = 0 AND ISTRIP = 1'
IF ((IGRIDI == 1).AND.(ISTRIP == 1)) THEN
!-----------------------------------------------------------------------------------------------!
!    Gap Pitch                                                                                  !
!-----------------------------------------------------------------------------------------------!
   RR=RMAX+CR
   IF (INTERP == 0) THEN
      CALL LININT(NRI,RI,CI   ,1,RR,CGAP   )
      CALL LININT(NRI,RI,SKEWI,1,RR,SKEWGAP)
      CALL LININT(NRI,RI,PTI  ,1,RR,PTGAP  )
      CALL LININT(NRI,RI,XRI  ,1,RR,XRGAP  )
   ELSEIF (INTERP == 1) THEN
      CALL INTK1 (NRI,RI,CI   ,1,RR,CGAP   )
      CALL INTK1 (NRI,RI,SKEWI,1,RR,SKEWGAP)
      CALL INTK1 (NRI,RI,PTI  ,1,RR,PTGAP  )
      CALL INTK1 (NRI,RI,XRI  ,1,RR,XRGAP  )
   ELSEIF (INTERP == 2) THEN
      CALL SPLINT(NRI,RI,CI   ,1,RR,CGAP   )
      CALL SPLINT(NRI,RI,SKEWI,1,RR,SKEWGAP)
      CALL SPLINT(NRI,RI,PTI  ,1,RR,PTGAP  )
      CALL SPLINT(NRI,RI,XRI  ,1,RR,XRGAP  )
   END IF !(INTERP)
   CGAP   =2.D0*CGAP !Non-Dimensional Quantities by Propeller Radius
   PTGAP  =2.D0*PTGAP
   PGAP   =2.D0*PGAP
   XRGAP  =2.D0*XRGAP
   SKEWGAP=SKEWGAP*PI/180.D0
!-----------------------------------------------------------------------------------------------!
   IF (PGAP == 0.D0) THEN
      PHI  =DATAN2(PTGAP,2.D0*PI*RR)
   ELSE !(PGAP)
      PHI  =DATAN2(PGAP ,2.D0*PI*RR)
      XRGAP  =XRGAP  -0.5D0*CGAP*DSIN(DATAN2(PTGAP,2.D0*PI*RR))   +0.5D0*CGAP*DSIN(PHI)
      SKEWGAP=SKEWGAP-0.5D0*CGAP*DCOS(DATAN2(PTGAP,2.D0*PI*RR))/RR+0.5D0*CGAP*DCOS(PHI)/RR
   END IF !(PGAP)
!-----------------------------------------------------------------------------------------------!
!    Loop on Chordwise Stations                                                                 !
!-----------------------------------------------------------------------------------------------!
   DO I=1,NC1
      IF (ISC == 0) THEN
         XC=  XRGAP+CGAP*(SCN(I)-0.5D0)*DSIN(PHI)
         TC=SKEWGAP+CGAP*(SCN(I)-0.5D0)*DCOS(PHI)/RR
         XP(NC1-I+1,NRP1+1)=XC-YFN(I,NRP1)*CGAP*DCOS(PHI)
         XP(NC1+I-1,NRP1+1)=XC-YBN(I,NRP1)*CGAP*DCOS(PHI)
         CALL NOZZLEDEF('INNER',XP(NC1-I+1,NRP1+1),RP(NC1-I+1,NRP1+1))
         CALL NOZZLEDEF('INNER',XP(NC1+I-1,NRP1+1),RP(NC1+I-1,NRP1+1))
         TP(NC1-I+1,NRP1+1)=TC+YFN(I,NRP1)*CGAP*DSIN(PHI)/RP(NC1-I+1,NRP1+1)
         TP(NC1+I-1,NRP1+1)=TC+YBN(I,NRP1)*CGAP*DSIN(PHI)/RP(NC1+I-1,NRP1+1)
         YP(NC1-I+1,NRP1+1)=RP(NC1-I+1,NRP1+1)*DCOS(TP(NC1-I+1,NRP1+1))
         YP(NC1+I-1,NRP1+1)=RP(NC1+I-1,NRP1+1)*DCOS(TP(NC1+I-1,NRP1+1))
         ZP(NC1-I+1,NRP1+1)=RP(NC1-I+1,NRP1+1)*DSIN(TP(NC1-I+1,NRP1+1))
         ZP(NC1+I-1,NRP1+1)=RP(NC1+I-1,NRP1+1)*DSIN(TP(NC1+I-1,NRP1+1))
      ELSEIF (ISC == 1) THEN
         XP(NC1-I+1,NRP1+1)=XRGAP+CGAP*(SFN(I,NRP1)-0.5D0)*DSIN(PHI)-YFN(I,NRP1)*CGAP*DCOS(PHI)
         XP(NC1+I-1,NRP1+1)=XRGAP+CGAP*(SBN(I,NRP1)-0.5D0)*DSIN(PHI)-YBN(I,NRP1)*CGAP*DCOS(PHI)
         CALL NOZZLEDEF('INNER',XP(NC1-I+1,NRP1+1),RP(NC1-I+1,NRP1+1))
         CALL NOZZLEDEF('INNER',XP(NC1+I-1,NRP1+1),RP(NC1+I-1,NRP1+1))
         TP(NC1-I+1,NRP1+1)=SKEWGAP+CGAP*(SFN(I,NRP1)-0.5D0)*DCOS(PHI)/RP(NC1-I+1,NRP1+1)+ &
                                            YFN(I,NRP1)*CGAP*DSIN(PHI)/RP(NC1-I+1,NRP1+1)
         TP(NC1+I-1,NRP1+1)=SKEWGAP+CGAP*(SBN(I,NRP1)-0.5D0)*DCOS(PHI)/RP(NC1+I-1,NRP1+1)+ &
                                            YBN(I,NRP1)*CGAP*DSIN(PHI)/RP(NC1+I-1,NRP1+1)
         YP(NC1-I+1,NRP1+1)=RP(NC1-I+1,NRP1+1)*DCOS(TP(NC1-I+1,NRP1+1))
         YP(NC1+I-1,NRP1+1)=RP(NC1+I-1,NRP1+1)*DCOS(TP(NC1+I-1,NRP1+1))
         ZP(NC1-I+1,NRP1+1)=RP(NC1-I+1,NRP1+1)*DSIN(TP(NC1-I+1,NRP1+1))
         ZP(NC1+I-1,NRP1+1)=RP(NC1+I-1,NRP1+1)*DSIN(TP(NC1+I-1,NRP1+1))
      END IF !(ISC)
   END DO !I=1,NC1
END IF !((IGRIDI == 1).AND.(ISTRIP == 1))
!-----------------------------------------------------------------------------------------------!
!    Closed Tip                                                                                 !
!-----------------------------------------------------------------------------------------------!
IF (ISTRIP == 2) THEN
   DO I=1,NC1
      XP(NC1-I+1,NRP1+1)=0.5D0*(XP(NC1-I+1,NRP1)+XP(NC1+I-1,NRP1))
      XP(NC1+I-1,NRP1+1)=XP(NC1-I+1,NRP1+1)
      YP(NC1-I+1,NRP1+1)=0.5D0*(YP(NC1-I+1,NRP1)+YP(NC1+I-1,NRP1))
      YP(NC1+I-1,NRP1+1)=YP(NC1-I+1,NRP1+1)
      ZP(NC1-I+1,NRP1+1)=0.5D0*(ZP(NC1-I+1,NRP1)+ZP(NC1+I-1,NRP1))
      ZP(NC1+I-1,NRP1+1)=ZP(NC1-I+1,NRP1+1)
      TP(NC1-I+1,NRP1+1)=0.5D0*(TP(NC1-I+1,NRP1)+TP(NC1+I-1,NRP1))
      TP(NC1+I-1,NRP1+1)=TP(NC1-I+1,NRP1+1)
   END DO !I=1,NC1
END IF !(ISTRIP == 2)
!-----------------------------------------------------------------------------------------------!
!    Write Blade Grid in Tecplot Format                                                         !
!-----------------------------------------------------------------------------------------------!
IF (ISTRIP == 0) THEN
   K=1
   WRITE(20,100) ' ZONE T="BLADE',K,'" F=POINT, I=',NCP1,' J=',NRP1
   DO J=1,NRP1
      DO I=1,NCP1
         RP(I,J)=DSQRT(YP(I,J)*YP(I,J)+ZP(I,J)*ZP(I,J))
         WRITE(20,110) XP(I,J),YP(I,J),ZP(I,J)
      END DO !I=1,NCP1
   END DO !J=1,NRP1
   DO K=2,NB
      PHI=DFLOAT(K-1)/DFLOAT(NB)*2.D0*PI
      WRITE(20,100) ' ZONE T="BLADE',K,'" F=POINT, I=',NCP1,' J=',NRP1
      DO J=1,NRP1
         DO I=1,NCP1
            IF (IGEOM /= 'WINGGEOM') WRITE(20,110) XP(I,J),RP(I,J)*DCOS(TP(I,J)+PHI), &
                                                                       RP(I,J)*DSIN(TP(I,J)+PHI)
            IF (IGEOM == 'WINGGEOM') WRITE(20,110) XP(I,J),-YP(I,J),ZP(I,J)
         END DO !I=1,NCP1
      END DO !J=1,NRP1
   END DO !K=2,NB
ELSE !(ISTRIP)
   K=1
   WRITE(20,100) ' ZONE T="BLADE',K,'" F=POINT, I=',NCP1,' J=',(NRP1+1)
   DO J=1,(NRP1+1)
      DO I=1,NCP1
         RP(I,J)=DSQRT(YP(I,J)*YP(I,J)+ZP(I,J)*ZP(I,J))
         WRITE(20,110) XP(I,J),YP(I,J),ZP(I,J)
      END DO !I=1,NCP1
   END DO !J=1,(NRP1+1)
   DO K=2,NB
      PHI=DFLOAT(K-1)/DFLOAT(NB)*2.D0*PI
      WRITE(20,100) ' ZONE T="BLADE',K,'" F=POINT, I=',NCP1,' J=',(NRP1+1)
      DO J=1,(NRP1+1)
         DO I=1,NCP1
            IF (IGEOM /= 'WINGGEOM') WRITE(20,110) XP(I,J),RP(I,J)*DCOS(TP(I,J)+PHI), &
                                                                       RP(I,J)*DSIN(TP(I,J)+PHI)
            IF (IGEOM == 'WINGGEOM') WRITE(20,110) XP(I,J),-YP(I,J),ZP(I,J)
         END DO !I=1,NCP1
      END DO !J=1,(NRP1+1)
   END DO !K=2,NB
END IF !(ISTRIP)
!-----------------------------------------------------------------------------------------------!
!    Deallocate Variables                                                                       !
!-----------------------------------------------------------------------------------------------!
DEALLOCATE(RI,CI,PTI,XRI,SKEWI,T0I,F0I)
DEALLOCATE(R,C,PT,XR,SKEW,T0,F0)
DEALLOCATE(SCN,YBN,YFN,YBI,YFI,YB,YF)
IF (IOERR == 0) DEALLOCATE(SB,SF)
IF (IOERR <  0) DEALLOCATE(SC,SCI)
IF (ISC   == 1) DEALLOCATE(S,SBN,SFN)
!-----------------------------------------------------------------------------------------------!
!    Formats                                                                                    !
!-----------------------------------------------------------------------------------------------!
100 FORMAT(A,I4,A,I4,A,I4)
110 FORMAT(3(2X,E23.16))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE BLADEGRID
!-----------------------------------------------------------------------------------------------!
