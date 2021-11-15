!-----------------------------------------------------------------------------------------------!
!    Wetted right-hand-side                                                                     !
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
SUBROUTINE SOLVERHSWET(JJ,II,TT)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 25112013, J. Baltazar, version 1.0                                             !
!    Modified  : 10112013, J. Baltazar, bounds checked                                          !
!    Modified  : 04122014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!    Modified  : 02022015, J. Baltazar, Gap sources                                             !
!    Modified  : 03062016, J. Baltazar, 2016 version 1.3                                        !
!    Modified  : 25102016, J. Baltazar, 2016 version 1.4                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,J2,K,L,KB,JJ,II,TT
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4),TE(4),TL
DOUBLE PRECISION :: VWX,VWY,VWZ,NCNP,DCPT,UINFSQ,RTEX,RTEY,RTEZ
DOUBLE PRECISION :: X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: POT
!-----------------------------------------------------------------------------------------------!
!    Rewind Matrices                                                                            !
!-----------------------------------------------------------------------------------------------!
REWIND(21)
READ  (21) ((DIJ(I,J,1),I=1,NPAN),J=1,NPAN)
!-----------------------------------------------------------------------------------------------!
SI=0.D0
!-----------------------------------------------------------------------------------------------!
!    Hub Panels                                                                                 !
!-----------------------------------------------------------------------------------------------!
!    First Half Sector                                                                          !
!-----------------------------------------------------------------------------------------------!
DO KB=1,NB
   DO J=1,NHT
      DO I=1,NHX
         K=(J-1)*NHX+I
!-----------------------------------------------------------------------------------------------!
!    Define Panel                                                                               !
!-----------------------------------------------------------------------------------------------!
         TL=DFLOAT(KB-1)*2.D0*PI/DFLOAT(NB)
         TE(1)=TH(I+1,J  )+TL
         TE(2)=TH(I  ,J  )+TL
         TE(3)=TH(I  ,J+1)+TL
         TE(4)=TH(I+1,J+1)+TL
         XX(1)=XH(I+1,J  )
         XX(2)=XH(I  ,J  )
         XX(3)=XH(I  ,J+1)
         XX(4)=XH(I+1,J+1)
         YY(1)=RH(I+1,J  )*DCOS(TE(1))
         YY(2)=RH(I  ,J  )*DCOS(TE(2))
         YY(3)=RH(I  ,J+1)*DCOS(TE(3))
         YY(4)=RH(I+1,J+1)*DCOS(TE(4))
         ZZ(1)=RH(I+1,J  )*DSIN(TE(1))
         ZZ(2)=RH(I  ,J  )*DSIN(TE(2))
         ZZ(3)=RH(I  ,J+1)*DSIN(TE(3))
         ZZ(4)=RH(I+1,J+1)*DSIN(TE(4))
!-----------------------------------------------------------------------------------------------!
!    Compute Panel Centroid Data                                                                !
!-----------------------------------------------------------------------------------------------!
         CALL PANEL(XX,YY,ZZ,X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0)
!-----------------------------------------------------------------------------------------------!
!    For Flat Panel Redefine Corner Points                                                      !
!-----------------------------------------------------------------------------------------------!
         IF (IPAN == 0) CALL PANELFLAT(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0)
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
         IF (IROTOR == 0) THEN
            CALL VWAKE(TT,X0,Y0,Z0,KB,IFREQ,VWX,VWY,VWZ)
            IF (TT == 0) CALL VWAKE(0,X0,Y0,Z0,KB,0,VWX,VWY,VWZ)
            SOURCEH(I,J,KB)=VWX*UU(JJ)/PI*UNX0+ & !x component
                            VWY*UU(JJ)/PI*UNY0+ & !y component
                            VWZ*UU(JJ)/PI*UNZ0    !z component
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
         ELSEIF (IROTOR == 1) THEN
            CALL VWAKE(TT,X0,Y0,Z0,KB,IFREQ,VWX,VWY,VWZ)
            IF (TT == 0) CALL VWAKE(0,X0,Y0,Z0,KB,0,VWX,VWY,VWZ)
            SOURCEH(I,J,KB)=VWX/UU(JJ)*UNX0+ & !x component
                            VWY/UU(JJ)*UNY0+ & !y component
                            VWZ/UU(JJ)*UNZ0    !z component
         END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
         DO L=1,NPAN
            SI(L)=SI(L)+SIJ(L,K,KB)*SOURCEH(I,J,KB)
         END DO !L=1,NPAN
!-----------------------------------------------------------------------------------------------!
      END DO !I=1,NHX
   END DO !J=1,NHT
END DO !KB=1,NB
!-----------------------------------------------------------------------------------------------!
!    Second Half Sector                                                                         !
!-----------------------------------------------------------------------------------------------!
DO KB=1,NB
   DO J=NHT2,NHTT1
      DO I=1,NHX
         K=(J-2)*NHX+I
!-----------------------------------------------------------------------------------------------!
!    Define Panel                                                                               !
!-----------------------------------------------------------------------------------------------!
         TL=DFLOAT(KB-1)*2.D0*PI/DFLOAT(NB)
         TE(1)=TH(I+1,J  )+TL
         TE(2)=TH(I  ,J  )+TL
         TE(3)=TH(I  ,J+1)+TL
         TE(4)=TH(I+1,J+1)+TL
         XX(1)=XH(I+1,J  )
         XX(2)=XH(I  ,J  )
         XX(3)=XH(I  ,J+1)
         XX(4)=XH(I+1,J+1)
         YY(1)=RH(I+1,J  )*DCOS(TE(1))
         YY(2)=RH(I  ,J  )*DCOS(TE(2))
         YY(3)=RH(I  ,J+1)*DCOS(TE(3))
         YY(4)=RH(I+1,J+1)*DCOS(TE(4))
         ZZ(1)=RH(I+1,J  )*DSIN(TE(1))
         ZZ(2)=RH(I  ,J  )*DSIN(TE(2))
         ZZ(3)=RH(I  ,J+1)*DSIN(TE(3))
         ZZ(4)=RH(I+1,J+1)*DSIN(TE(4))
!-----------------------------------------------------------------------------------------------!
!    Compute Panel Centroid Data                                                                !
!-----------------------------------------------------------------------------------------------!
         CALL PANEL(XX,YY,ZZ,X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0)
!-----------------------------------------------------------------------------------------------!
!    For Flat Panel Redefine Corner Points                                                      !
!-----------------------------------------------------------------------------------------------!
         IF (IPAN == 0) CALL PANELFLAT(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0)
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
         IF (IROTOR == 0) THEN
            CALL VWAKE(TT,X0,Y0,Z0,KB,IFREQ,VWX,VWY,VWZ)
            IF (TT == 0) CALL VWAKE(0,X0,Y0,Z0,KB,0,VWX,VWY,VWZ)
            SOURCEH(I,J-1,KB)=VWX*UU(JJ)/PI*UNX0+ & !x component
                              VWY*UU(JJ)/PI*UNY0+ & !y component
                              VWZ*UU(JJ)/PI*UNZ0    !z component
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
         ELSEIF (IROTOR == 1) THEN
            CALL VWAKE(TT,X0,Y0,Z0,KB,IFREQ,VWX,VWY,VWZ)
            IF (TT == 0) CALL VWAKE(0,X0,Y0,Z0,KB,0,VWX,VWY,VWZ)
            SOURCEH(I,J-1,KB)=VWX/UU(JJ)*UNX0+ & !x component
                              VWY/UU(JJ)*UNY0+ & !y component
                              VWZ/UU(JJ)*UNZ0    !z component
         END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
         DO L=1,NPAN
            SI(L)=SI(L)+SIJ(L,K,KB)*SOURCEH(I,J-1,KB)
         END DO !L=1,NPAN
!-----------------------------------------------------------------------------------------------!
      END DO !I=1,NHX
   END DO !J=NHT2,NHTT1
END DO !KB=1,NB
!-----------------------------------------------------------------------------------------------!
IF (TT > 0) THEN
   ALLOCATE(POT(NHX,NHTP))
   DO KB=2,NB
      IF (IRED == 0) CALL PERIODICFLOW(TT,KB,NHX,NHTP,NT,POTH   ,POT)
      IF (IRED == 1) CALL PERIODICFLOW(TT,KB,NHX,NHTP,NT,POTHWET,POT)
      DO J=1,NHTP
         DO I=1,NHX
            K=(J-1)*NHX+I
!-----------------------------------------------------------------------------------------------!
            DO L=1,NPAN
               SI(L)=SI(L)-DIJ(L,K,KB)*POT(I,J)
            END DO !DO L=1,NPAN
!-----------------------------------------------------------------------------------------------!
         END DO !I=1,NHX
      END DO !J=1,NHTP
   END DO !KB=2,NB
   DEALLOCATE(POT)
END IF !(TT > 0)
!-----------------------------------------------------------------------------------------------!
!    Blade Panels                                                                               !
!-----------------------------------------------------------------------------------------------!
DO KB=1,NB
   DO J=1,(NRP-1)
      DO I=1,NCP
         K=NHPAN+(J-1)*NCP+I
!-----------------------------------------------------------------------------------------------!
!    Define Panel                                                                               !
!-----------------------------------------------------------------------------------------------!
         TL=DFLOAT(KB-1)*2.D0*PI/DFLOAT(NB)
         TE(1)=TP(I  ,J  )+TL
         TE(2)=TP(I+1,J  )+TL
         TE(3)=TP(I+1,J+1)+TL
         TE(4)=TP(I  ,J+1)+TL
         XX(1)=XP(I  ,J  )
         XX(2)=XP(I+1,J  )
         XX(3)=XP(I+1,J+1)
         XX(4)=XP(I  ,J+1)
         YY(1)=RP(I  ,J  )*DCOS(TE(1))
         YY(2)=RP(I+1,J  )*DCOS(TE(2))
         YY(3)=RP(I+1,J+1)*DCOS(TE(3))
         YY(4)=RP(I  ,J+1)*DCOS(TE(4))
         ZZ(1)=RP(I  ,J  )*DSIN(TE(1))
         ZZ(2)=RP(I+1,J  )*DSIN(TE(2))
         ZZ(3)=RP(I+1,J+1)*DSIN(TE(3))
         ZZ(4)=RP(I  ,J+1)*DSIN(TE(4))
	 IF ((IROTOR == -1).AND.(KB == 2)) THEN
            XX(1)= XP(I+1,J  )
            XX(2)= XP(I  ,J  )
            XX(3)= XP(I  ,J+1)
            XX(4)= XP(I+1,J+1)
            YY(1)=-YP(I+1,J  )
            YY(2)=-YP(I  ,J  )
            YY(3)=-YP(I  ,J+1)
            YY(4)=-YP(I+1,J+1)
            ZZ(1)= ZP(I+1,J  )
            ZZ(2)= ZP(I  ,J  )
            ZZ(3)= ZP(I  ,J+1)
            ZZ(4)= ZP(I+1,J+1)
         END IF !((IROTOR == -1).AND.(KB == 2))
!-----------------------------------------------------------------------------------------------!
!    Compute Panel Centroid Data                                                                !
!-----------------------------------------------------------------------------------------------!
         CALL PANEL(XX,YY,ZZ,X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0)
!-----------------------------------------------------------------------------------------------!
!    For Flat Panel Redefine Corner Points                                                      !
!-----------------------------------------------------------------------------------------------!
         IF (IPAN == 0) CALL PANELFLAT(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0)
!-----------------------------------------------------------------------------------------------!
!    Wing Case                                                                                  !
!-----------------------------------------------------------------------------------------------!
         IF (IROTOR == -1) THEN
!*          CALL VWAKE(0,XP0(I,J),YP0(I,J),ZP0(I,J),1,0,VWX,VWY,VWZ)
            SOURCEP(I,J,KB)=DCOSD(UU(JJ))*UNX0+DSIND(UU(JJ))*UNZ0
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
         ELSEIF (IROTOR == 0) THEN
            CALL VWAKE(TT,X0,Y0,Z0,KB,IFREQ,VWX,VWY,VWZ)
            IF (TT == 0) CALL VWAKE(0,X0,Y0,Z0,KB,0,VWX,VWY,VWZ)
            SOURCEP(I,J,KB)=VWX*UU(JJ)/PI*UNX0+ &     !x component
                           (VWY*UU(JJ)/PI-Z0)*UNY0+ & !y component
                           (VWZ*UU(JJ)/PI+Y0)*UNZ0    !z component
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
         ELSEIF (IROTOR == 1) THEN
            CALL VWAKE(TT,X0,Y0,Z0,KB,IFREQ,VWX,VWY,VWZ)
            IF (TT == 0) CALL VWAKE(0,X0,Y0,Z0,KB,0,VWX,VWY,VWZ)
            SOURCEP(I,J,KB)=VWX/UU(JJ)*UNX0+ &     !x component
                           (VWY/UU(JJ)-Z0)*UNY0+ & !y component
                           (VWZ/UU(JJ)+Y0)*UNZ0    !z component
         END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
         DO L=1,NPAN
            SI(L)=SI(L)+SIJ(L,K,KB)*SOURCEP(I,J,KB)
         END DO !L=1,NPAN
!-----------------------------------------------------------------------------------------------!
      END DO !I=1,NCP
   END DO !J=1,(NRP-1)
END DO !KB=1,NB
!-----------------------------------------------------------------------------------------------!
!    Blade Gap Panels                                                                           !
!-----------------------------------------------------------------------------------------------!
DO KB=1,NB
   J=NRP
   DO I=1,NCP
      K=NHPAN+(J-1)*NCP+I
!-----------------------------------------------------------------------------------------------!
!    Define Panel                                                                               !
!-----------------------------------------------------------------------------------------------!
      TL=DFLOAT(KB-1)*2.D0*PI/DFLOAT(NB)
      TE(1)=TP(I  ,J  )+TL
      TE(2)=TP(I+1,J  )+TL
      TE(3)=TP(I+1,J+1)+TL
      TE(4)=TP(I  ,J+1)+TL
      XX(1)=XP(I  ,J  )
      XX(2)=XP(I+1,J  )
      XX(3)=XP(I+1,J+1)
      XX(4)=XP(I  ,J+1)
      YY(1)=RP(I  ,J  )*DCOS(TE(1))
      YY(2)=RP(I+1,J  )*DCOS(TE(2))
      YY(3)=RP(I+1,J+1)*DCOS(TE(3))
      YY(4)=RP(I  ,J+1)*DCOS(TE(4))
      ZZ(1)=RP(I  ,J  )*DSIN(TE(1))
      ZZ(2)=RP(I+1,J  )*DSIN(TE(2))
      ZZ(3)=RP(I+1,J+1)*DSIN(TE(3))
      ZZ(4)=RP(I  ,J+1)*DSIN(TE(4))
      IF ((IROTOR == -1).AND.(KB == 2)) THEN
         XX(1)= XP(I+1,J  )
         XX(2)= XP(I  ,J  )
         XX(3)= XP(I  ,J+1)
         XX(4)= XP(I+1,J+1)
         YY(1)=-YP(I+1,J  )
         YY(2)=-YP(I  ,J  )
         YY(3)=-YP(I  ,J+1)
         YY(4)=-YP(I+1,J+1)
         ZZ(1)= ZP(I+1,J  )
         ZZ(2)= ZP(I  ,J  )
         ZZ(3)= ZP(I  ,J+1)
         ZZ(4)= ZP(I+1,J+1)
      END IF !((IROTOR == -1).AND.(KB == 2))
!-----------------------------------------------------------------------------------------------!
!    Compute Panel Centroid Data                                                                !
!-----------------------------------------------------------------------------------------------!
      CALL PANEL(XX,YY,ZZ,X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0)
!-----------------------------------------------------------------------------------------------!
!    For Flat Panel Redefine Corner Points                                                      !
!-----------------------------------------------------------------------------------------------!
      IF (IPAN == 0) CALL PANELFLAT(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0)
!-----------------------------------------------------------------------------------------------!
      NCNP =UNXC0(I)*UNXP0(I,J)+UNYC0(I)*UNYP0(I,J)+UNZC0(I)*UNZP0(I,J) !dot product for KB=1
      IF (TT < IDINT(P)*(KB-1)) THEN
         DCPT=CPP(I,J,0)-CPP(NCP-I+1,J,0)
      ELSE
         DCPT=CPP(I,J,(TT-IDINT(P)*(KB-1)))-CPP(NCP-I+1,J,(TT-IDINT(P)*(KB-1)))
      END IF !(TT < IDINT(P)*(KB-1))
!-----------------------------------------------------------------------------------------------!
!    Wing Case                                                                                  !
!-----------------------------------------------------------------------------------------------!
      IF (IROTOR == -1) THEN
!*       CALL VWAKE(0,XP0(I,J),YP0(I,J),ZP0(I,J),1,0,VWX,VWY,VWZ)
         SOURCEP(I,J,KB)=DCOSD(UU(JJ))*UNX0+DSIND(UU(JJ))*UNZ0
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
      ELSEIF (IROTOR == 0) THEN
         CALL VWAKE(TT,X0,Y0,Z0,KB,IFREQ,VWX,VWY,VWZ)
         IF (TT == 0) CALL VWAKE(0,X0,Y0,Z0,KB,0,VWX,VWY,VWZ)
         UINFSQ=(VWX*UU(JJ)/PI)**2+(VWY*UU(JJ)/PI-Z0)**2+(VWZ*UU(JJ)/PI+Y0)**2
         IF (II == 0) THEN
            SOURCEP(I,J,KB)=VWX*UU(JJ)/PI*UNX0+ &     !x component
                           (VWY*UU(JJ)/PI-Z0)*UNY0+ & !y component
                           (VWZ*UU(JJ)/PI+Y0)*UNZ0    !z component
!*       ELSEIF (II == 1) THEN
!*          SOURCEP(I,J,KB)=VWX*UU(JJ)/PI*UNX0+ &     !x component
!*                         (VWY*UU(JJ)/PI-Z0)*UNY0+ & !y component
!*                         (VWZ*UU(JJ)/PI+Y0)*UNZ0+ & !z component
!*                          DSIGN(DSQRT(DABS(DCPT)*UINFSQ),DCPT)*CQ*NCNP
         ELSE !(II)
            SOURCEP(I,J,KB)=(VWX*UU(JJ)/PI*UNX0+ &     !x component
                            (VWY*UU(JJ)/PI-Z0)*UNY0+ & !y component
                            (VWZ*UU(JJ)/PI+Y0)*UNZ0+ & !z component
                             DSIGN(DSQRT(DABS(DCPT)*UINFSQ),DCPT)*CQ*NCNP)*FREX+ &
                             SOURCEP(I,J,KB)*(1.D0-FREX)
         END IF !(II)
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
      ELSEIF (IROTOR == 1) THEN
         CALL VWAKE(TT,X0,Y0,Z0,KB,IFREQ,VWX,VWY,VWZ)
         IF (TT == 0) CALL VWAKE(0,X0,Y0,Z0,KB,0,VWX,VWY,VWZ)
         UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ)-Z0)**2+(VWZ/UU(JJ)+Y0)**2
         IF (II == 0) THEN
            SOURCEP(I,J,KB)=VWX/UU(JJ)*UNX0+ &     !x component
                           (VWY/UU(JJ)-Z0)*UNY0+ & !y component
                           (VWZ/UU(JJ)+Y0)*UNZ0    !z component
!*       ELSEIF (II == 1) THEN
!*          SOURCEP(I,J,KB)=VWX/UU(JJ)*UNX0+ &     !x component
!*                         (VWY/UU(JJ)-Z0)*UNY0+ & !y component
!*                         (VWZ/UU(JJ)+Y0)*UNZ0+ & !z component
!*                          DSIGN(DSQRT(DABS(DCPT)*UINFSQ),DCPT)*CQ*NCNP
         ELSE !(II)
            SOURCEP(I,J,KB)=(VWX/UU(JJ)*UNX0+ &     !x component
                            (VWY/UU(JJ)-Z0)*UNY0+ & !y component
                            (VWZ/UU(JJ)+Y0)*UNZ0+ & !z component
                             DSIGN(DSQRT(DABS(DCPT)*UINFSQ),DCPT)*CQ*NCNP)*FREX+ &
                             SOURCEP(I,J,KB)*(1.D0-FREX)
         END IF !(II)
      END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
      DO L=1,NPAN
         SI(L)=SI(L)+SIJ(L,K,KB)*SOURCEP(I,J,KB)
      END DO !L=1,NPAN
!-----------------------------------------------------------------------------------------------!
   END DO !I=1,NCP
END DO !KB=1,NB
!-----------------------------------------------------------------------------------------------!
IF (TT > 0) THEN
   ALLOCATE(POT(NCP,NRP))
   DO KB=2,NB
      IF (IRED == 0) CALL PERIODICFLOW(TT,KB,NCP,NRP,NT,POTP   ,POT)
      IF (IRED == 1) CALL PERIODICFLOW(TT,KB,NCP,NRP,NT,POTPWET,POT)
      DO J=1,NRP
         DO I=1,NCP
            K=NHPAN+(J-1)*NCP+I
!-----------------------------------------------------------------------------------------------!
            DO L=1,NPAN
               SI(L)=SI(L)-DIJ(L,K,KB)*POT(I,J)
            END DO !DO L=1,NPAN
!-----------------------------------------------------------------------------------------------!
         END DO !I=1,NCP
      END DO !J=1,NRP
   END DO !KB=2,NB
   DEALLOCATE(POT)
END IF !(TT > 0)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake Panels                                                                          !
!-----------------------------------------------------------------------------------------------!
IF (TT > 0) THEN
   DO J=1,NRW
      DO I=2,NPW
         K=(J-1)*NPW+I
!-----------------------------------------------------------------------------------------------!
         DO L=1,NPAN
            IF (I == 2) THEN
               IF (IRED == 0) SI(L)=SI(L)+(KIJ(L,J,2)+WIJ(L,K,1)*0.5D0)*POTPW   (I,J,TT)
               IF (IRED == 1) SI(L)=SI(L)+(KIJ(L,J,2)+WIJ(L,K,1)*0.5D0)*POTPWWET(I,J,TT)
            ELSE !(I == 2)
               IF (IRED == 0) SI(L)=SI(L)+(WIJ(L,K-1,1)+WIJ(L,K,1))*0.5D0*POTPW   (I,J,TT)
               IF (IRED == 1) SI(L)=SI(L)+(WIJ(L,K-1,1)+WIJ(L,K,1))*0.5D0*POTPWWET(I,J,TT)
            END IF !(I == 2)
         END DO !L=1,NPAN
!-----------------------------------------------------------------------------------------------!
      END DO !I=2,NPW
   END DO !J=1,NRW
!-----------------------------------------------------------------------------------------------!
   ALLOCATE(POT(NPW,NRW))
   DO KB=2,NB
      IF (IRED == 0) CALL PERIODICFLOW(TT,KB,NPW,NRW,NT,POTPW  ,POT)
      IF (IRED == 1) CALL PERIODICFLOW(TT,KB,NPW,NRW,NT,POTPWWET,POT)
      DO J=1,NRW
         DO I=1,NPW
            K=(J-1)*NPW+I
!-----------------------------------------------------------------------------------------------!
            DO L=1,NPAN
               SI(L)=SI(L)+WIJ(L,K,KB)*POT(I,J)
            END DO !L=1,NPAN
!-----------------------------------------------------------------------------------------------!
         END DO !I=1,NPW
      END DO !J=1,NRW
   END DO !KB=2,NB
   DEALLOCATE(POT)
END IF !(TT > 0)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Panels                                                                              !
!-----------------------------------------------------------------------------------------------!
!    First Half Sector                                                                          !
!-----------------------------------------------------------------------------------------------!
DO KB=1,NB
   DO J=1,NNT
      DO I=1,NNXT1
         K=NHPAN+NPPAN+(J-1)*NNXT1+I
!-----------------------------------------------------------------------------------------------!
!    Define Panel                                                                               !
!-----------------------------------------------------------------------------------------------!
         TL=DFLOAT(KB-1)*2.D0*PI/DFLOAT(NB)
         TE(1)=TN(I  ,J  )+TL
         TE(2)=TN(I+1,J  )+TL
         TE(3)=TN(I+1,J+1)+TL
         TE(4)=TN(I  ,J+1)+TL
         XX(1)=XN(I  ,J  )
         XX(2)=XN(I+1,J  )
         XX(3)=XN(I+1,J+1)
         XX(4)=XN(I  ,J+1)
         YY(1)=RN(I  ,J  )*DCOS(TE(1))
         YY(2)=RN(I+1,J  )*DCOS(TE(2))
         YY(3)=RN(I+1,J+1)*DCOS(TE(3))
         YY(4)=RN(I  ,J+1)*DCOS(TE(4))
         ZZ(1)=RN(I  ,J  )*DSIN(TE(1))
         ZZ(2)=RN(I+1,J  )*DSIN(TE(2))
         ZZ(3)=RN(I+1,J+1)*DSIN(TE(3))
         ZZ(4)=RN(I  ,J+1)*DSIN(TE(4))
!-----------------------------------------------------------------------------------------------!
!    Compute Panel Centroid Data                                                                !
!-----------------------------------------------------------------------------------------------!
         CALL PANEL(XX,YY,ZZ,X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0)
!-----------------------------------------------------------------------------------------------!
!    For Flat Panel Redefine Corner Points                                                      !
!-----------------------------------------------------------------------------------------------!
         IF (IPAN == 0) CALL PANELFLAT(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0)
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
         IF (IROTOR == 0) THEN
            CALL VWAKE(TT,X0,Y0,Z0,KB,IFREQ,VWX,VWY,VWZ)
            IF (TT == 0) CALL VWAKE(0,X0,Y0,Z0,KB,0,VWX,VWY,VWZ)
            SOURCEN(I,J,KB)=VWX*UU(JJ)/PI*UNX0+ & !x component
                            VWY*UU(JJ)/PI*UNY0+ & !y component
                            VWZ*UU(JJ)/PI*UNZ0    !z component
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
         ELSEIF (IROTOR == 1) THEN
            CALL VWAKE(TT,X0,Y0,Z0,KB,IFREQ,VWX,VWY,VWZ)
            IF (TT == 0) CALL VWAKE(0,X0,Y0,Z0,KB,0,VWX,VWY,VWZ)
            SOURCEN(I,J,KB)=VWX/UU(JJ)*UNX0+ & !x component
                            VWY/UU(JJ)*UNY0+ & !y component
                            VWZ/UU(JJ)*UNZ0    !z component
         END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
         DO L=1,NPAN
            SI(L)=SI(L)+SIJ(L,K,KB)*SOURCEN(I,J,KB)
         END DO !L=1,NPAN
!-----------------------------------------------------------------------------------------------!
      END DO !I=1,NNXT1
   END DO !J=1,NNT
END DO !KB=1,NB
!-----------------------------------------------------------------------------------------------!
!    Second Half Sector                                                                         !
!-----------------------------------------------------------------------------------------------!
DO KB=1,NB
   DO J=NNT2,NNTT1
      DO I=1,NNXT1
         K=NHPAN+NPPAN+(J-2)*NNXT1+I
!-----------------------------------------------------------------------------------------------!
!    Define Panel                                                                               !
!-----------------------------------------------------------------------------------------------!
         TL=DFLOAT(KB-1)*2.D0*PI/DFLOAT(NB)
         TE(1)=TN(I  ,J  )+TL
         TE(2)=TN(I+1,J  )+TL
         TE(3)=TN(I+1,J+1)+TL
         TE(4)=TN(I  ,J+1)+TL
         XX(1)=XN(I  ,J  )
         XX(2)=XN(I+1,J  )
         XX(3)=XN(I+1,J+1)
         XX(4)=XN(I  ,J+1)
         YY(1)=RN(I  ,J  )*DCOS(TE(1))
         YY(2)=RN(I+1,J  )*DCOS(TE(2))
         YY(3)=RN(I+1,J+1)*DCOS(TE(3))
         YY(4)=RN(I  ,J+1)*DCOS(TE(4))
         ZZ(1)=RN(I  ,J  )*DSIN(TE(1))
         ZZ(2)=RN(I+1,J  )*DSIN(TE(2))
         ZZ(3)=RN(I+1,J+1)*DSIN(TE(3))
         ZZ(4)=RN(I  ,J+1)*DSIN(TE(4))
!-----------------------------------------------------------------------------------------------!
!    Compute Panel Centroid Data                                                                !
!-----------------------------------------------------------------------------------------------!
         CALL PANEL(XX,YY,ZZ,X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0)
!-----------------------------------------------------------------------------------------------!
!    For Flat Panel Redefine Corner Points                                                      !
!-----------------------------------------------------------------------------------------------!
         IF (IPAN == 0) CALL PANELFLAT(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0)
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
         IF (IROTOR == 0) THEN
            CALL VWAKE(TT,X0,Y0,Z0,KB,IFREQ,VWX,VWY,VWZ)
            IF (TT == 0) CALL VWAKE(0,X0,Y0,Z0,KB,0,VWX,VWY,VWZ)
            SOURCEN(I,J-1,KB)=VWX*UU(JJ)/PI*UNX0+ & !x component
                              VWY*UU(JJ)/PI*UNY0+ & !y component
                              VWZ*UU(JJ)/PI*UNZ0    !z component
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
         ELSEIF (IROTOR == 1) THEN
            CALL VWAKE(TT,X0,Y0,Z0,KB,IFREQ,VWX,VWY,VWZ)
            IF (TT == 0) CALL VWAKE(0,X0,Y0,Z0,KB,0,VWX,VWY,VWZ)
            SOURCEN(I,J-1,KB)=VWX/UU(JJ)*UNX0+ & !x component
                              VWY/UU(JJ)*UNY0+ & !y component
                              VWZ/UU(JJ)*UNZ0    !z component
         END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
         DO L=1,NPAN
            SI(L)=SI(L)+SIJ(L,K,KB)*SOURCEN(I,J-1,KB)
         END DO !L=1,NPAN
!-----------------------------------------------------------------------------------------------!
      END DO !I=1,NNXT1
   END DO !J=NNT2,NNTT1
END DO !KB=1,NB
!-----------------------------------------------------------------------------------------------!
IF (TT > 0) THEN
   ALLOCATE(POT(NNXT1,NNTP))
   DO KB=2,NB
      IF (IRED == 0) CALL PERIODICFLOW(TT,KB,NNXT1,NNTP,NT,POTN   ,POT)
      IF (IRED == 1) CALL PERIODICFLOW(TT,KB,NNXT1,NNTP,NT,POTNWET,POT)
      DO J=1,NNTP
         DO I=1,NNXT1
            K=NHPAN+NPPAN+(J-1)*NNXT1+I
!-----------------------------------------------------------------------------------------------!
            DO L=1,NPAN
               SI(L)=SI(L)-DIJ(L,K,KB)*POT(I,J)
            END DO !DO L=1,NPAN
!-----------------------------------------------------------------------------------------------!
         END DO !I=1,NNXT1
      END DO !J=1,NNTP
   END DO !KB=2,NB
   DEALLOCATE(POT)
END IF !(TT > 0)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake Panels                                                                         !
!-----------------------------------------------------------------------------------------------!
IF (TT > 0) THEN
   DO J=1,NNTP
      DO I=2,NNW
         K=NPWPAN+(J-1)*NNW+I
!-----------------------------------------------------------------------------------------------!
         DO L=1,NPAN
            IF (I == 2) THEN
               IF (IRED == 0) SI(L)=SI(L)+(KIJ(L,NRW+J,2)+WIJ(L,K,1)*0.5D0)*POTNW   (I,J,TT)
               IF (IRED == 1) SI(L)=SI(L)+(KIJ(L,NRW+J,2)+WIJ(L,K,1)*0.5D0)*POTNWWET(I,J,TT)
            ELSE !(I == 2)
               IF (IRED == 0) SI(L)=SI(L)+(WIJ(L,K-1,1)+WIJ(L,K,1))*0.5D0*POTNW   (I,J,TT)
               IF (IRED == 1) SI(L)=SI(L)+(WIJ(L,K-1,1)+WIJ(L,K,1))*0.5D0*POTNWWET(I,J,TT)
            END IF !(I == 2)
         END DO !L=1,NPAN
!-----------------------------------------------------------------------------------------------!
      END DO !I=2,NNW
   END DO !J=1,NNTP
!-----------------------------------------------------------------------------------------------!
   ALLOCATE(POT(NNW,NNTP))
   DO KB=2,NB
      IF (IRED == 0) CALL PERIODICFLOW(TT,KB,NNW,NNTP,NT,POTNW   ,POT)
      IF (IRED == 1) CALL PERIODICFLOW(TT,KB,NNW,NNTP,NT,POTNWWET,POT)
      DO J=1,NNTP
         DO I=1,NNW
            K=(J-1)*NNW+I
!-----------------------------------------------------------------------------------------------!
            DO L=1,NPAN
               SI(L)=SI(L)+WIJ(L,K,KB)*POT(I,J)
            END DO !L=1,NPAN
!-----------------------------------------------------------------------------------------------!
         END DO !I=1,NNW
      END DO !J=1,NNTP
   END DO !KB=2,NB
   DEALLOCATE(POT)
END IF !(TT > 0)
!-----------------------------------------------------------------------------------------------!
!    Correction For Trailing Edge Thickness in the Kutta Condition                              !
!-----------------------------------------------------------------------------------------------!
DO L=1,NPAN
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
   DO J=1,NRW
      J2=JI-1+J
      RTEX=XP0(NCP,J2)-XP0(1,J2)
      RTEY=YP0(NCP,J2)-YP0(1,J2)
      RTEZ=ZP0(NCP,J2)-ZP0(1,J2)
!-----------------------------------------------------------------------------------------------!
!    Wing Case                                                                                  !
!-----------------------------------------------------------------------------------------------!
      IF (IROTOR == -1) THEN
!*       CALL VWAKE(0,XP0(NCP,J2),YP0(NCP,J2),ZP0(NCP,J2),1,0,VWX,VWY,VWZ)
         DO KB=1,NB
            DO K=((J-1)*NPW+1),(J*NPW)
               SI(L)=SI(L)+WIJ(L,K,1)*(DCOSD(UU(JJ))*RTEX+DSIND(UU(JJ))*RTEZ)
            END DO !K=((J-1)*NPW+1),(J*NPW)
         END DO !KB=1,NB
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
      ELSEIF (IROTOR == 0) THEN
         IF (NT == 0) THEN !Steady Case
            DO KB=1,NB
               CALL VWAKE(0,XP0(NCP,J2),YP0(NCP,J2),ZP0(NCP,J2),KB,0,VWX,VWY,VWZ)
               SI(L)=SI(L)+WIJ(L,J,KB)*(VWX*UU(JJ)/PI*RTEX+ &              !x component
                                       (VWY*UU(JJ)/PI-ZP0(NCP,J2))*RTEY+ & !y component
                                       (VWZ*UU(JJ)/PI+YP0(NCP,J2))*RTEZ)   !z component
            END DO !KB=1,NB
         ELSEIF (TT == 0) THEN
            DO KB=1,NB
               CALL VWAKE(0,XP0(NCP,J2),YP0(NCP,J2),ZP0(NCP,J2),KB,0,VWX,VWY,VWZ)
               DO K=((J-1)*NPW+1),(J*NPW)
                  SI(L)=SI(L)+WIJ(L,K,KB)*(VWX*UU(JJ)/PI*RTEX+ &              !x component
                                          (VWY*UU(JJ)/PI-ZP0(NCP,J2))*RTEY+ & !y component
                                          (VWZ*UU(JJ)/PI+YP0(NCP,J2))*RTEZ)   !z component
               END DO !K=((J-1)*NPW+1),(J*NPW)
            END DO !KB=1,NB
         ELSE !(TT)
            CALL VWAKE(TT,XP0(NCP,J2),YP0(NCP,J2),ZP0(NCP,J2),1,IFREQ,VWX,VWY,VWZ)
            SI(L)=SI(L)+KIJ(L,J,1)*(VWX*UU(JJ)/PI*RTEX+ &              !x component
                                   (VWY*UU(JJ)/PI-ZP0(NCP,J2))*RTEY+ & !y component
                                   (VWZ*UU(JJ)/PI+YP0(NCP,J2))*RTEZ)   !z component
         END IF !(TT)
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
      ELSEIF (IROTOR == 1) THEN
         IF (NT == 0) THEN
            DO KB=1,NB
               CALL VWAKE(0,XP0(NCP,J2),YP0(NCP,J2),ZP0(NCP,J2),KB,0,VWX,VWY,VWZ)
               SI(L)=SI(L)+WIJ(L,J,KB)*(VWX/UU(JJ)*RTEX+ &              !x component
                                       (VWY/UU(JJ)-ZP0(NCP,J2))*RTEY+ & !y component
                                       (VWZ/UU(JJ)+YP0(NCP,J2))*RTEZ)   !z component
            END DO !(KB=1,NB)
         ELSEIF (TT == 0) THEN
            DO KB=1,NB
               CALL VWAKE(0,XP0(NCP,J2),YP0(NCP,J2),ZP0(NCP,J2),KB,0,VWX,VWY,VWZ)
               DO K=((J-1)*NPW+1),(J*NPW)
                  SI(L)=SI(L)+WIJ(L,K,KB)*(VWX/UU(JJ)*RTEX+ &              !x component
                                          (VWY/UU(JJ)-ZP0(NCP,J2))*RTEY+ & !y component
                                          (VWZ/UU(JJ)+YP0(NCP,J2))*RTEZ)   !z component
               END DO !K=((J-1)*NPW+1),(J*NPW)
            END DO !(KB=1,NB)
         ELSE !(TT)
            CALL VWAKE(TT,XP0(NCP,J2),YP0(NCP,J2),ZP0(NCP,J2),1,IFREQ,VWX,VWY,VWZ)
            SI(L)=SI(L)+KIJ(L,J,1)*(VWX/UU(JJ)*RTEX+ &              !x component
                                   (VWY/UU(JJ)-ZP0(NCP,J2))*RTEY+ & !y component
                                   (VWZ/UU(JJ)+YP0(NCP,J2))*RTEZ)   !z component
         END IF !(TT)
      END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
   END DO !J=1,NRW
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
   DO J=1,NNTP
      RTEX=XN0(NNX,J)-XN0(NNX1,J)
      RTEY=YN0(NNX,J)-YN0(NNX1,J)
      RTEZ=ZN0(NNX,J)-ZN0(NNX1,J)
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
      IF (IROTOR == 0) THEN
         IF (NT == 0) THEN
            DO KB=1,NB
               CALL VWAKE(0,XN0(NNX,J),YN0(NNX,J),ZN0(NNX,J),KB,0,VWX,VWY,VWZ)
               SI(L)=SI(L)+WIJ(L,NRW+J,KB)*(VWX*UU(JJ)/PI*RTEX+ &             !x component
                                           (VWY*UU(JJ)/PI-ZN0(NNX,J))*RTEY+ & !y component
                                           (VWZ*UU(JJ)/PI+YN0(NNX,J))*RTEZ)   !z component
            END DO !KB=1,NB
         ELSEIF (TT == 0) THEN
            DO KB=1,NB
               CALL VWAKE(0,XN0(NNX,J),YN0(NNX,J),ZN0(NNX,J),KB,0,VWX,VWY,VWZ)
               DO K=((J-1)*NNW+1),(J*NNW)
                  SI(L)=SI(L)+WIJ(L,NPWPAN+K,KB)*(VWX*UU(JJ)/PI*RTEX+ &             !x component
                                                 (VWY*UU(JJ)/PI-ZN0(NNX,J))*RTEY+ & !y component
                                                 (VWZ*UU(JJ)/PI+YN0(NNX,J))*RTEZ)   !z component
               END DO !K=((J-1)*NNW+1),(J*NNW)
            END DO !KB=1,NB
         ELSE !(TT)
            CALL VWAKE(TT,XN0(NNX,J),YN0(NNX,J),ZN0(NNX,J),1,IFREQ,VWX,VWY,VWZ)
            SI(L)=SI(L)+KIJ(L,NRW+J,1)*(VWX*UU(JJ)/PI*RTEX+ &             !x component
                                       (VWY*UU(JJ)/PI-ZN0(NNX,J))*RTEY+ & !y component
                                       (VWZ*UU(JJ)/PI+YN0(NNX,J))*RTEZ)   !z component
         END IF !(TT)
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
      ELSEIF (IROTOR == 1) THEN
         IF (NT == 0) THEN
            DO KB=1,NB
               CALL VWAKE(0,XN0(NNX,J),YN0(NNX,J),ZN0(NNX,J),KB,0,VWX,VWY,VWZ)
               SI(L)=SI(L)+WIJ(L,NRW+J,KB)*(VWX/UU(JJ)*RTEX+ &             !x component
                                           (VWY/UU(JJ)-ZN0(NNX,J))*RTEY+ & !y component
                                           (VWZ/UU(JJ)+YN0(NNX,J))*RTEZ)   !z component
            END DO !KB=1,NB
         ELSEIF (TT == 0) THEN
            DO KB=1,NB
               CALL VWAKE(0,XN0(NNX,J),YN0(NNX,J),ZN0(NNX,J),KB,0,VWX,VWY,VWZ)
               DO K=((J-1)*NNW+1),(J*NNW)
                  SI(L)=SI(L)+WIJ(L,NPWPAN+K,KB)*(VWX/UU(JJ)*RTEX+ &             !x component
                                                 (VWY/UU(JJ)-ZN0(NNX,J))*RTEY+ & !y component
                                                 (VWZ/UU(JJ)+YN0(NNX,J))*RTEZ)   !z component
               END DO !K=((J-1)*NNW+1),(J*NNW)
            END DO !KB=1,NB
         ELSE !(TT)
            CALL VWAKE(TT,XN0(NNX,J),YN0(NNX,J),ZN0(NNX,J),1,IFREQ,VWX,VWY,VWZ)
            SI(L)=SI(L)+KIJ(L,NRW+J,1)*(VWX/UU(JJ)*RTEX+ &             !x component
                                       (VWY/UU(JJ)-ZN0(NNX,J))*RTEY+ & !y component
                                       (VWZ/UU(JJ)+YN0(NNX,J))*RTEZ)   !z component
         END IF !(TT)
      END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
   END DO !J=1,NNTP
!-----------------------------------------------------------------------------------------------!
!    Save the Right Hand Side                                                                   !
!-----------------------------------------------------------------------------------------------!
   RHS(L)=SI(L)
!-----------------------------------------------------------------------------------------------!
END DO !L=1,NPAN
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE SOLVERHSWET
!-----------------------------------------------------------------------------------------------!
