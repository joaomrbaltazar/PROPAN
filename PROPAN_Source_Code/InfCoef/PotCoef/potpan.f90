!-----------------------------------------------------------------------------------------------!
!    Subroutine POTPANH                                                                         !
!                                                                                               !
!    Computation of the potentials at a field point due to a constant source/dipole             !
!    distribution of strength -4*pi on a quadrilateral (hyperboloidal) panel.                   !
!                                                                                               !
!    Near field formulas:                                                                       !
!                                                                                               !
!    Exact analytical integration for the dipole potential. Approximate analytical integration  !
!    for the source potential. (exact for a plane panel)                                        !
!                                                                                               !
!    Author  : J.A.C. Falcão de Campos                                                          !
!    Revised : Joao Baltazar, IST, May 2005                                                     !
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
SUBROUTINE POTPANH(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0,X,Y,Z,PS,PD)
!-----------------------------------------------------------------------------------------------!
!    INPUT:                                                                                     !
!                                                                                               !
!    XX(4),YY(4),ZZ(4):      CARTESIAN COORDINATES OF PANEL CORNER POINTS                       !
!                                                                                               !
!    X0,Y0,Z0         :      CARTESIAN COORDINATES OF PANEL CENTROID                            !
!                                                                                               !
!    UNX0,UNY0,UNZ0   :      CARTESIAN COMPONENTS OF UNIT NORMAL AT THE CENTROID                !
!                                                                                               !
!    X,Y,Z            :      CARTESIAN COORDINATES OF FIELD POINT                               !
!                                                                                               !
!    OUTPUT:                                                                                    !
!                                                                                               !
!    PS               :      POTENTIAL DUE TO THE SOURCE DISTRIBUTION                           !
!                                                                                               !
!    PD               :      POTENTIAL DUE TO THE DIPOLE DISTRIBUTION                           !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4),X0,Y0,Z0,UNX0,UNY0,UNZ0,X,Y,Z,PS,PD
DOUBLE PRECISION :: RX(4),RY(4),RZ(4),RM(4)
DOUBLE PRECISION :: A1X(4),A1Y(4),A1Z(4),A1M(4),A2X(4),A2Y(4),A2Z(4),A2M(4)
DOUBLE PRECISION :: VNX(4),VNY(4),VNZ(4),VNM(4)
DOUBLE PRECISION :: RVA1X(4),RVA1Y(4),RVA1Z(4),RVA1M(4)
DOUBLE PRECISION :: RVA2X(4),RVA2Y(4),RVA2Z(4),RVA2M(4)
DOUBLE PRECISION :: RDVN(4),RDA1(4),RDA2(4),RA12(4),PD1,PD2,PD3,PD4
DOUBLE PRECISION :: RVA1DN(4),RVA2DN(4),ARG1,ARG2,AS1,AS2,PS1,PS2,PS3,PS4
DOUBLE PRECISION :: H,RN0,ZETA(4),TOL,PI,TWOPI
!-----------------------------------------------------------------------------------------------!
!    Computations                                                                               !
!-----------------------------------------------------------------------------------------------!
PI=4.D0*DATAN(1.D0)
TWOPI=2.D0*PI
TOL=1.D-16
!-----------------------------------------------------------------------------------------------!
!    Corner Vectors                                                                             !
!-----------------------------------------------------------------------------------------------!
RX(:)=XX(:)-X
RY(:)=YY(:)-Y
RZ(:)=ZZ(:)-Z
RM(:)=DSQRT(RX(:)*RX(:)+RY(:)*RY(:)+RZ(:)*RZ(:))
!-----------------------------------------------------------------------------------------------!
!    Tangent Vectors (A1,A2)                                                                    !
!-----------------------------------------------------------------------------------------------!
!    CSI=-1 ETA=-1                                                                              !
!-----------------------------------------------------------------------------------------------!
A1X(1)=0.5D0*(XX(2)-XX(1))
A1Y(1)=0.5D0*(YY(2)-YY(1))
A1Z(1)=0.5D0*(ZZ(2)-ZZ(1))
A2X(1)=0.5D0*(XX(4)-XX(1))
A2Y(1)=0.5D0*(YY(4)-YY(1))
A2Z(1)=0.5D0*(ZZ(4)-ZZ(1))
!-----------------------------------------------------------------------------------------------!
!    CSI=1 ETA=-1                                                                               !
!-----------------------------------------------------------------------------------------------!
A1X(2)=0.5D0*(XX(2)-XX(1))
A1Y(2)=0.5D0*(YY(2)-YY(1))
A1Z(2)=0.5D0*(ZZ(2)-ZZ(1))
A2X(2)=0.5D0*(XX(3)-XX(2))
A2Y(2)=0.5D0*(YY(3)-YY(2))
A2Z(2)=0.5D0*(ZZ(3)-ZZ(2))
!-----------------------------------------------------------------------------------------------!
!    CSI=1 ETA=1                                                                                !
!-----------------------------------------------------------------------------------------------!
A1X(3)=0.5D0*(XX(3)-XX(4))
A1Y(3)=0.5D0*(YY(3)-YY(4))
A1Z(3)=0.5D0*(ZZ(3)-ZZ(4))
A2X(3)=0.5D0*(XX(3)-XX(2))
A2Y(3)=0.5D0*(YY(3)-YY(2))
A2Z(3)=0.5D0*(ZZ(3)-ZZ(2))
!-----------------------------------------------------------------------------------------------!
!    CSI=-1 ETA=1                                                                               !
!-----------------------------------------------------------------------------------------------!
A1X(4)=0.5D0*(XX(3)-XX(4))
A1Y(4)=0.5D0*(YY(3)-YY(4))
A1Z(4)=0.5D0*(ZZ(3)-ZZ(4))
A2X(4)=0.5D0*(XX(4)-XX(1))
A2Y(4)=0.5D0*(YY(4)-YY(1))
A2Z(4)=0.5D0*(ZZ(4)-ZZ(1))
!-----------------------------------------------------------------------------------------------!
A1M=DSQRT(A1X*A1X+A1Y*A1Y+A1Z*A1Z)
A2M=DSQRT(A2X*A2X+A2Y*A2Y+A2Z*A2Z)
!-----------------------------------------------------------------------------------------------!
!    Normal Vectors                                                                             !
!-----------------------------------------------------------------------------------------------!
VNX=A1Y*A2Z-A1Z*A2Y
VNY=A1Z*A2X-A1X*A2Z
VNZ=A1X*A2Y-A1Y*A2X
VNM=DSQRT(VNX*VNX+VNY*VNY+VNZ*VNZ)
!-----------------------------------------------------------------------------------------------!
!    R X A1 and R X A2                                                                          !
!-----------------------------------------------------------------------------------------------!
RVA1X=RY*A1Z-RZ*A1Y
RVA1Y=RZ*A1X-RX*A1Z
RVA1Z=RX*A1Y-RY*A1X
RVA1M=DSQRT(RVA1X*RVA1X+RVA1Y*RVA1Y+RVA1Z*RVA1Z)
!-----------------------------------------------------------------------------------------------!
RVA2X=RY*A2Z-RZ*A2Y
RVA2Y=RZ*A2X-RX*A2Z
RVA2Z=RX*A2Y-RY*A2X
RVA2M=DSQRT(RVA2X*RVA2X+RVA2Y*RVA2Y+RVA2Z*RVA2Z)
!-----------------------------------------------------------------------------------------------!
!    R . (A1 X A2)                                                                              !
!-----------------------------------------------------------------------------------------------!
RDVN=RX*VNX+RY*VNY+RZ*VNZ
!-----------------------------------------------------------------------------------------------!
!    R . A1 and R . A2                                                                          !
!-----------------------------------------------------------------------------------------------!
RDA1=RX*A1X+RY*A1Y+RZ*A1Z
RDA2=RX*A2X+RY*A2Y+RZ*A2Z
!-----------------------------------------------------------------------------------------------!
!    (R X A1) . (R X A2)                                                                        !
!-----------------------------------------------------------------------------------------------!
RA12=RVA1X*RVA2X+RVA1Y*RVA2Y+RVA1Z*RVA2Z
!-----------------------------------------------------------------------------------------------!
!    (R X A1).N and (R X A2).N                                                                  !
!-----------------------------------------------------------------------------------------------!
RVA1DN=RVA1X*UNX0+RVA1Y*UNY0+RVA1Z*UNZ0
RVA2DN=RVA2X*UNX0+RVA2Y*UNY0+RVA2Z*UNZ0
!-----------------------------------------------------------------------------------------------!
!    Vertex 1                                                                                   !
!-----------------------------------------------------------------------------------------------!
!    Dipole Potential due to the Vertex                                                         !
!-----------------------------------------------------------------------------------------------!
H=RM(1)*RDVN(1)
PD1=DATAN2(H,RA12(1))
IF ((H == 0.D0).AND.(RA12(1) == 0.D0)) PD1=0.D0
!-----------------------------------------------------------------------------------------------!
!    Approximate Source Potential due to the Vertex                                             !
!-----------------------------------------------------------------------------------------------!
IF (RVA1M(1) > TOL) THEN
   ARG1=RM(1)+RM(2)+2.D0*A1M(1)
   ARG2=RM(1)+RM(2)-2.D0*A1M(1)
   IF (DABS(ARG2) > TOL) THEN
      AS1=DLOG(ARG1/ARG2)
      PS1=RVA1DN(1)*AS1/A1M(1)
   ELSE !(ARG2 > TOL)
      PS1=0.D0
   END IF !(ARG2 > TOL)
ELSE !(RVA1M(1) > TOL)
   PS1=0.D0
END IF !(RVA1M(1) > TOL)
!-----------------------------------------------------------------------------------------------!
!    Vertex 2                                                                                   !
!-----------------------------------------------------------------------------------------------!
!    Dipole Potential due to the Vertex                                                         !
!-----------------------------------------------------------------------------------------------!
H=RM(2)*RDVN(2)
PD2=DATAN2(H,RA12(2))
IF ((H == 0.D0).AND.(RA12(2) == 0.D0)) PD2=0.D0
!-----------------------------------------------------------------------------------------------!
!    Approximate Source Potential due to the Vertex                                             !
!-----------------------------------------------------------------------------------------------!
IF (RVA2M(2) > TOL) THEN
   ARG1=RM(2)+RM(3)+2.D0*A2M(2)
   ARG2=RM(2)+RM(3)-2.D0*A2M(2)
   IF (DABS(ARG2) > TOL) THEN
      AS2=DLOG(ARG1/ARG2)
      PS2=RVA2DN(2)*AS2/A2M(2)
   ELSE !(ARG2 > TOL)
      PS2=0.D0
   END IF !(ARG2 > TOL)
ELSE !(RVA2M(2) > TOL)
   PS2=0.D0
END IF !(RVA2M(2) > TOL)
!-----------------------------------------------------------------------------------------------!
!    Vertex 3                                                                                   !
!-----------------------------------------------------------------------------------------------!
!    Dipole Potential due to the Vertex                                                         !
!-----------------------------------------------------------------------------------------------!
H=RM(3)*RDVN(3)
PD3=DATAN2(H,RA12(3))
IF ((H == 0.D0).AND.(RA12(3) == 0.D0)) PD3=0.D0
!-----------------------------------------------------------------------------------------------!
!    Approximate Source Potential due to the Vertex                                             !
!-----------------------------------------------------------------------------------------------!
IF (RVA1M(3) > TOL) THEN
   ARG1=RM(3)+RM(4)+2.D0*A1M(3)
   ARG2=RM(3)+RM(4)-2.D0*A1M(3)
   IF (DABS(ARG2) > TOL) THEN
      AS1=DLOG(ARG1/ARG2)
      PS3=RVA1DN(3)*AS1/A1M(3)
   ELSE !(ARG2 > TOL)
      PS3=0.D0
   END IF !(ARG2 > TOL)
ELSE !(RVA1M(3) > TOL)
   PS3=0.D0
END IF !(RVA1M(3) > TOL)
!-----------------------------------------------------------------------------------------------!
!    Vertex 4                                                                                   !
!-----------------------------------------------------------------------------------------------!
!    Dipole Potential due to the Vertex                                                         !
!-----------------------------------------------------------------------------------------------!
H=RM(4)*RDVN(4)
PD4=DATAN2(H,RA12(4))
IF ((H == 0.D0).AND.(RA12(4) == 0.D0)) PD4=0.D0
!-----------------------------------------------------------------------------------------------!
!    Approximate Source Potential due to the Vertex                                             !
!-----------------------------------------------------------------------------------------------!
IF (RVA2M(4) > TOL) THEN
   ARG1=RM(4)+RM(1)+2.D0*A2M(4)
   ARG2=RM(4)+RM(1)-2.D0*A2M(4)
   IF (DABS(ARG2) > TOL) THEN
      AS2=DLOG(ARG1/ARG2)
      PS4=RVA2DN(4)*AS2/A2M(4)
   ELSE !(ARG2 > TOL)
      PS4=0.D0
   END IF !(ARG2 > TOL)
ELSE !(RVA2M(4) > TOL)
   PS4=0.D0
END IF !(RVA2M(4) > TOL)
!-----------------------------------------------------------------------------------------------!
!    Total Potentials                                                                           !
!-----------------------------------------------------------------------------------------------!
!    Dipole Corrections                                                                         !
!-----------------------------------------------------------------------------------------------!
RN0=(X0-X)*UNX0+(Y0-Y)*UNY0+(Z0-Z)*UNZ0
IF (DABS(RN0) > 1.D-10) THEN
   IF ((PD2 > 0.D0).AND.(PD4 > 0.D0)) THEN
      IF ((PD1 < 0.D0).AND.(PD3 < 0.D0)) THEN
         PD1=PD1+TWOPI
         PD3=PD3+TWOPI
         PD=PD1-PD2+PD3-PD4
         IF (DABS(PD) > TWOPI) THEN
            PD1=PD1-TWOPI
            PD3=PD3-TWOPI
         END IF
      ELSEIF (PD1 < 0.D0) THEN
         PD1=PD1+TWOPI
         PD=PD1-PD2+PD3-PD4
         IF (DABS(PD) > TWOPI) PD1=PD1-TWOPI
      ELSEIF (PD3 < 0.D0) THEN
         PD3=PD3+TWOPI
         PD=PD1-PD2+PD3-PD4
         IF (DABS(PD) > TWOPI) PD3=PD3-TWOPI
      END IF
   END IF
   IF ((PD2 < 0.D0).AND.(PD4 < 0.D0)) THEN
      IF ((PD1 > 0.D0).AND.(PD3 > 0.D0)) THEN
         PD1=PD1-TWOPI
         PD3=PD3-TWOPI
         PD=PD1-PD2+PD3-PD4
         IF (DABS(PD) > TWOPI) THEN
            PD1=PD1+TWOPI
            PD3=PD3+TWOPI
         END IF
      ELSEIF (PD1 > 0.D0) THEN
         PD1=PD1-TWOPI
         PD=PD1-PD2+PD3-PD4
         IF (DABS(PD) > TWOPI) PD1=PD1+TWOPI
      ELSEIF (PD3 > 0.D0) THEN
         PD3=PD3-TWOPI
         PD=PD1-PD2+PD3-PD4
         IF (DABS(PD) > TWOPI) PD3=PD3+TWOPI
      END IF
   END IF
ELSE !(DABS(RN0) > 1.D-10)
   IF (RDVN(1) > 0.D0) THEN
      IF (PD1 > 0.D0) PD1=PD1-PI
   END IF !(RDVN(1) > 0.D0)
   IF (RDVN(2) > 0.D0) THEN
      IF (PD2 > 0.D0) PD2=PD2-PI
   END IF !(RDVN(2) > 0.D0)
   IF (RDVN(3) > 0.D0) THEN
      IF (PD3 > 0.D0) PD3=PD3-PI
   END IF !(RDVN(3) > 0.D0)
   IF (RDVN(4) > 0.D0) THEN
      IF (PD4 > 0.D0) PD4=PD4-PI
   END IF !(RDVN(4) > 0.D0)
!-----------------------------------------------------------------------------------------------!
!    Flat Panel                                                                                 !
!-----------------------------------------------------------------------------------------------!
   ZETA(:)=DABS((XX(:)-X0)*UNX0+(YY(:)-Y0)*UNY0+(ZZ(:)-Z0)*UNZ0)
   IF (SUM(ZETA) < 1.D-10) THEN
      PD1=0.D0
      PD2=0.D0
      PD3=0.D0
      PD4=0.D0
   END IF !(SUM(ZETA) < 1.D-10)
END IF !(DABS(RN0) > 1.D-10)
!-----------------------------------------------------------------------------------------------!
!    Dipole                                                                                     !
!-----------------------------------------------------------------------------------------------!
PD=PD1-PD2+PD3-PD4
PD=-PD
!-----------------------------------------------------------------------------------------------!
!    Source                                                                                     !
!-----------------------------------------------------------------------------------------------!
PS=PS1+PS2-PS3-PS4+RN0*PD
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE POTPANH
!===============================================================================================!
!    Subroutine POTPAN3                                                                         !
!                                                                                               !
!    Computation of the potentials at a field point due to a constant source/dipole             !
!    distribution of strength -4*pi on a triangular panel.                                      !
!                                                                                               !
!    Near field formulas:                                                                       !
!                                                                                               !
!    Exact analytical integration                                                               !
!                                                                                               !
!    Created by: Joao Baltazar, IST, May 2005                                                   !
!    Modified  : 03072017, J. Baltazar, 2017 version 1.0                                        !
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
SUBROUTINE POTPAN3(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0,X,Y,Z,PS,PD)
!-----------------------------------------------------------------------------------------------!
!    INPUT:                                                                                     !
!                                                                                               !
!    XX(3),YY(3),ZZ(3): CARTESIAN COORDINATES OF PANEL CORNER POINTS                            !
!                                                                                               !
!    X0,Y0,Z0         : CARTESIAN COORDINATES OF PANEL CENTROID                                 !
!                                                                                               !
!    UNX0,UNY0,UNZ0   : CARTESIAN COMPONENTS OF UNIT NORMAL AT THE CENTROID                     !
!                                                                                               !
!    X,Y,Z            : CARTESIAN COORDINATES OF FIELD POINT                                    !
!                                                                                               !
!    OUTPUT:                                                                                    !
!                                                                                               !
!    PS               : POTENTIAL DUE TO THE SOURCE DISTRIBUTION                                !
!                                                                                               !
!    PD               : POTENTIAL DUE TO THE DIPOLE DISTRIBUTION                                !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
DOUBLE PRECISION :: XX(3),YY(3),ZZ(3),X,Y,Z,PS,PD,UNX0,UNY0,UNZ0,X0,Y0,Z0
DOUBLE PRECISION :: RX(3),RY(3),RZ(3),RM(3)
DOUBLE PRECISION :: A1X(3),A1Y(3),A1Z(3),A1M(3),A2X(3),A2Y(3),A2Z(3),A2M(3)
DOUBLE PRECISION :: VNX(3),VNY(3),VNZ(3),VNM(3)
DOUBLE PRECISION :: RVA1X(3),RVA1Y(3),RVA1Z(3),RVA1M(3)
DOUBLE PRECISION :: RVA2X(3),RVA2Y(3),RVA2Z(3),RVA2M(3)
DOUBLE PRECISION :: RDVN(3),RDA1(3),RDA2(3),RA12(3),PD1,PD2,PD3
DOUBLE PRECISION :: RVA1DN(3),RVA2DN(3),ARG1,ARG2,AS1,AS2,PS1,PS2,PS3
DOUBLE PRECISION :: H,RN0,TOL,PI
!-----------------------------------------------------------------------------------------------!
!    Computations                                                                               !
!-----------------------------------------------------------------------------------------------!
PI=4.D0*DATAN(1.D0)
TOL=1.D-16 !30
!-----------------------------------------------------------------------------------------------!
!    Corner Vectors                                                                             !
!-----------------------------------------------------------------------------------------------!
RX(:)=XX(:)-X
RY(:)=YY(:)-Y
RZ(:)=ZZ(:)-Z
RM(:)=DSQRT(RX(:)*RX(:)+RY(:)*RY(:)+RZ(:)*RZ(:))
!-----------------------------------------------------------------------------------------------!
!    Tangent Vectors (A1,A2)                                                                    !
!-----------------------------------------------------------------------------------------------!
!    CSI=0 ETA=0                                                                                !
!-----------------------------------------------------------------------------------------------!
A1X(1)=XX(2)-XX(1)
A1Y(1)=YY(2)-YY(1)
A1Z(1)=ZZ(2)-ZZ(1)
A2X(1)=XX(3)-XX(1)
A2Y(1)=YY(3)-YY(1)
A2Z(1)=ZZ(3)-ZZ(1)
!-----------------------------------------------------------------------------------------------!
!    CSI=1 ETA=0                                                                                !
!-----------------------------------------------------------------------------------------------!
A1X(2)=XX(2)-XX(1)
A1Y(2)=YY(2)-YY(1)
A1Z(2)=ZZ(2)-ZZ(1)
A2X(2)=XX(3)-XX(2)
A2Y(2)=YY(3)-YY(2)
A2Z(2)=ZZ(3)-ZZ(2)
!-----------------------------------------------------------------------------------------------!
!    CSI=0 ETA=1                                                                                !
!-----------------------------------------------------------------------------------------------!
A1X(3)=XX(3)-XX(1)
A1Y(3)=YY(3)-YY(1)
A1Z(3)=ZZ(3)-ZZ(1)
A2X(3)=XX(3)-XX(2)
A2Y(3)=YY(3)-YY(2)
A2Z(3)=ZZ(3)-ZZ(2)
!-----------------------------------------------------------------------------------------------!
A1M=DSQRT(A1X*A1X+A1Y*A1Y+A1Z*A1Z)
A2M=DSQRT(A2X*A2X+A2Y*A2Y+A2Z*A2Z)
!-----------------------------------------------------------------------------------------------!
!    Normal Vector                                                                              !
!-----------------------------------------------------------------------------------------------!
VNX=A1Y*A2Z-A1Z*A2Y
VNY=A1Z*A2X-A1X*A2Z
VNZ=A1X*A2Y-A1Y*A2X
VNM=DSQRT(VNX*VNX+VNY*VNY+VNZ*VNZ)
!-----------------------------------------------------------------------------------------------!
!    R X A1 and R X A2                                                                          !
!-----------------------------------------------------------------------------------------------!
RVA1X=RY*A1Z-RZ*A1Y
RVA1Y=RZ*A1X-RX*A1Z
RVA1Z=RX*A1Y-RY*A1X
RVA1M=DSQRT(RVA1X*RVA1X+RVA1Y*RVA1Y+RVA1Z*RVA1Z)
!-----------------------------------------------------------------------------------------------!
RVA2X=RY*A2Z-RZ*A2Y
RVA2Y=RZ*A2X-RX*A2Z
RVA2Z=RX*A2Y-RY*A2X
RVA2M=DSQRT(RVA2X*RVA2X+RVA2Y*RVA2Y+RVA2Z*RVA2Z)
!-----------------------------------------------------------------------------------------------!
!    R . (A1 X A2)                                                                              !
!-----------------------------------------------------------------------------------------------!
RDVN=RX*VNX+RY*VNY+RZ*VNZ
!-----------------------------------------------------------------------------------------------!
!    R . A1 and R . A2                                                                          !
!-----------------------------------------------------------------------------------------------!
RDA1=RX*A1X+RY*A1Y+RZ*A1Z
RDA2=RX*A2X+RY*A2Y+RZ*A2Z
!-----------------------------------------------------------------------------------------------!
!    (R X A1) . (R X A2)                                                                        !
!-----------------------------------------------------------------------------------------------!
RA12=RVA1X*RVA2X+RVA1Y*RVA2Y+RVA1Z*RVA2Z
!-----------------------------------------------------------------------------------------------!
!    (R X A1).N and (R X A2).N                                                                  !
!-----------------------------------------------------------------------------------------------!
RVA1DN=RVA1X*UNX0+RVA1Y*UNY0+RVA1Z*UNZ0
RVA2DN=RVA2X*UNX0+RVA2Y*UNY0+RVA2Z*UNZ0
!-----------------------------------------------------------------------------------------------!
!    Vertex 1                                                                                   !
!-----------------------------------------------------------------------------------------------!
!    Dipole Potential due to the Vertex                                                         !
!-----------------------------------------------------------------------------------------------!
H=RM(1)*RDVN(1)
PD1=DATAN2(H,RA12(1))
IF (DSQRT(H**2+RA12(1)**2) < TOL) PD1=0.d0
!-----------------------------------------------------------------------------------------------!
!    Source Potetinal due to the Vertex                                                         !
!-----------------------------------------------------------------------------------------------!
IF (RVA1M(1) > TOL) THEN
   ARG1=RM(1)+RM(2)+A1M(1)
   ARG2=RM(1)+RM(2)-A1M(1)
   IF (DABS(ARG2) > TOL) THEN
      AS1=DLOG(ARG1/ARG2)
      PS1=RVA1DN(1)*AS1/A1M(1)
   ELSE !(ARG2 > TOL)
      PS1=0.D0
   END IF !(ARG2 > TOL)
ELSE !(RVA1M(1) > TOL)
   PS1=0.D0
END IF !(RVA1M(1) > TOL)
!-----------------------------------------------------------------------------------------------!
!    Vertex 2                                                                                   !
!-----------------------------------------------------------------------------------------------!
!    Dipole Potential due to the Vertex                                                         !
!-----------------------------------------------------------------------------------------------!
H=RM(2)*RDVN(2)
PD2=DATAN2(H,RA12(2))
IF (DSQRT(H**2+RA12(2)**2) < TOL) PD2=0.d0
!-----------------------------------------------------------------------------------------------!
!    Approximate Source Potential due to the Vertex                                             !
!-----------------------------------------------------------------------------------------------!
IF (RVA2M(2) > TOL) THEN
   ARG1=RM(2)+RM(3)+A2M(2)
   ARG2=RM(2)+RM(3)-A2M(2)
   IF (DABS(ARG2) > TOL) THEN
      AS2=DLOG(ARG1/ARG2)
      PS2=RVA2DN(2)*AS2/A2M(2)
   ELSE !(ARG2 > TOL)
      PS2=0.D0
   END IF !(ARG2 > TOL)
ELSE !(RVA2M(2) > TOL)
   PS2=0.D0
END IF !(RVA2M(2) > TOL)
!-----------------------------------------------------------------------------------------------!
!    Vertex 3                                                                                   !
!-----------------------------------------------------------------------------------------------!
!    Dipole Potential due to the Vertex                                                         !
!-----------------------------------------------------------------------------------------------!
H=RM(3)*RDVN(3)
PD3=DATAN2(H,RA12(3))
IF (DSQRT(H**2+RA12(3)**2) < TOL) PD3=0.d0
!-----------------------------------------------------------------------------------------------!
!    Approximate Source Potential due to the Vertex                                             !
!-----------------------------------------------------------------------------------------------!
IF (RVA1M(3) > TOL) THEN
   ARG1=RM(3)+RM(1)+A1M(3)
   ARG2=RM(3)+RM(1)-A1M(3)
   IF (DABS(ARG2) > TOL) THEN
      AS1=DLOG(ARG1/ARG2)
      PS3=RVA1DN(3)*AS1/A1M(3)
   ELSE !(ARG2 > TOL)
      PS3=0.D0
   END IF !(ARG2 > TOL)
ELSE !(RVA1M(3) > TOL)
   PS3=0.D0
END IF !(RVA1M(3) > TOL)
!-----------------------------------------------------------------------------------------------!
!    Total Potentials                                                                           !
!-----------------------------------------------------------------------------------------------!
!    Dipole                                                                                     !
!-----------------------------------------------------------------------------------------------!
PD=PD1-PD2+PD3
PD=-PD
!-----------------------------------------------------------------------------------------------!
!    Dipole Correction                                                                          !
!-----------------------------------------------------------------------------------------------!
RN0=(X0-X)*UNX0+(Y0-Y)*UNY0+(Z0-Z)*UNZ0
IF (DABS(RN0) < TOL) PD=0.D0
!-----------------------------------------------------------------------------------------------!
!    Source                                                                                     !
!-----------------------------------------------------------------------------------------------!
PS=PS1+PS2-PS3+RN0*PD
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE POTPAN3
!-----------------------------------------------------------------------------------------------!
