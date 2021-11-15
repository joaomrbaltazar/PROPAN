!-----------------------------------------------------------------------------------------------!
!    Computation of the velocities at a field point due to a constant source/dipole             !
!    distribution of strenth -4*pi on a quadrilateral (hyperboloidal) panel.                    !
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
!                                                                                               !                                                                                               !
!    Created by:   022007, J. Baltazar,                                                         !
!    Modified  :   122010, J. Baltazar, numerical tolerance                                     !
!    Modified  : 12052016, J. Baltazar, revised                                                 !
!-----------------------------------------------------------------------------------------------!
SUBROUTINE VELPAN(XX,YY,ZZ,X0,Y0,Z0,T1X,T1Y,T1Z,T2X,T2Y,T2Z,X,Y,Z,EPS,VS,VD)
!-----------------------------------------------------------------------------------------------!
!    INPUT:                                                                                     !
!                                                                                               !
!    XX(4),YY(4),ZZ(4):      CARTESIAN COORDINATES OF PANEL CORNER POINTS                       !
!                                                                                               !
!    X0,Y0,Z0         :      CARTESIAN COORDINATES OF PANEL CENTROID                            !
!                                                                                               !
!    T1X,T1Y,T1Z      :      COMPONENTS OF THE TANGENT VECTOR T1 AT THE CENTER POINT            !
!                                                                                               !
!    T2X,T2Y,T2Z      :      COMPONENTS OF THE TANGENT VECTOR T2 AT THE CENTER POINT            !
!                                                                                               !
!    X,Y,Z            :      CARTESIAN COORDINATES OF FIELD POINT                               !
!                                                                                               !
!    EPS              :      NUMERICAL TOLERANCE                                                !
!                                                                                               !
!    OUTPUT:                                                                                    !
!                                                                                               !
!    VS(3)            :      VELOCITY INDUCED BY A SOURCE DISTRIBUTION                          !
!                                                                                               !
!    VD(3)            :      VELOCITY INDUCED BY A DIPOLE DISTRIBUTION                          !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
INTEGER          :: K,K1
logical          :: flag
DOUBLE PRECISION :: PI,EPS,LIM,SDIV
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4),X0,Y0,Z0,X,Y,Z,XP(5),YP(5),ZP(5)
DOUBLE PRECISION :: T1X,T1Y,T1Z,T2X,T2Y,T2Z,XS,YS,ZS,XI(5),YI(5)
DOUBLE PRECISION :: A11,A12,A13,A21,A22,A23,A31,A32,A33
DOUBLE PRECISION :: D12,XI0,YI0,X0H,Y0H,Z0H
DOUBLE PRECISION :: XC1,YC1,ZC1,XC2,YC2,ZC2,XC3,YC3,ZC3,XC4,YC4,ZC4
DOUBLE PRECISION :: ARG1,ARG2,SOLUT,DD,MM,R1,R2,E1,E2,H1,H2,GAMMA
DOUBLE PRECISION :: DX,DY,DZ,R0,XH,YH,ZH,VX,VY,VZ
DOUBLE PRECISION :: DX1,DY1,DZ1,R12X,R12Y,R12Z,R12S2,R0R1,R0R2,KA
DOUBLE PRECISION :: VS(3),VD(3)
!-----------------------------------------------------------------------------------------------!
PI=4.D0*DATAN(1.D0)
LIM=1.D19
flag=.true.
!-----------------------------------------------------------------------------------------------!
!    Tangent Vector T1                                                                          !
!-----------------------------------------------------------------------------------------------!
D12=DSQRT(T1X*T1X+T1Y*T1Y+T1Z*T1Z)
A11=T1X/D12
A12=T1Y/D12
A13=T1Z/D12
!-----------------------------------------------------------------------------------------------!
!    Tangent Vector T2                                                                          !
!-----------------------------------------------------------------------------------------------!
D12=DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
A21=T2X/D12
A22=T2Y/D12
A23=T2Z/D12
!-----------------------------------------------------------------------------------------------!
!    Unit Normal Vector                                                                         !
!-----------------------------------------------------------------------------------------------!
A31=T1Y*T2Z-T2Y*T1Z
A32=T2X*T1Z-T1X*T2Z
A33=T1X*T2Y-T2X*T1Y
D12=DSQRT(A31*A31+A32*A32+A33*A33)
A31=A31/D12
A32=A32/D12
A33=A33/D12
!-----------------------------------------------------------------------------------------------!
!    Compute Cornerpoint (XC,YC,ZC) for the Average Surface Panel                               !
!    with Respect to the Reference Coordinate System                                            !
!-----------------------------------------------------------------------------------------------!
XC1=XX(1)-X0
YC1=YY(1)-Y0
ZC1=ZZ(1)-Z0
XC2=XX(2)-X0
YC2=YY(2)-Y0
ZC2=ZZ(2)-Z0
XC3=XX(3)-X0
YC3=YY(3)-Y0
ZC3=ZZ(3)-Z0
XC4=XX(4)-X0
YC4=YY(4)-Y0
ZC4=ZZ(4)-Z0
!-----------------------------------------------------------------------------------------------!
!    Transform Cornerpoints of the Average Panel to the Local Coordinate System                 !
!-----------------------------------------------------------------------------------------------!
XI(1)=A11*XC1+A12*YC1+A13*ZC1
YI(1)=A21*XC1+A22*YC1+A23*ZC1
XI(2)=A11*XC2+A12*YC2+A13*ZC2
YI(2)=A21*XC2+A22*YC2+A23*ZC2
XI(3)=A11*XC3+A12*YC3+A13*ZC3
YI(3)=A21*XC3+A22*YC3+A23*ZC3
XI(4)=A11*XC4+A12*YC4+A13*ZC4
YI(4)=A21*XC4+A22*YC4+A23*ZC4
XI(5)=XI(1)
YI(5)=YI(1)
!-----------------------------------------------------------------------------------------------!
!    Compute Distance R0 from Offbody Point to Panel Centroid                                   !
!-----------------------------------------------------------------------------------------------!
DX=X-X0
DY=Y-Y0
DZ=Z-Z0
R0=DSQRT(DX*DX+DY*DY+DZ*DZ)
!-----------------------------------------------------------------------------------------------!
XH=A11*DX+A12*DY+A13*DZ
YH=A21*DX+A22*DY+A23*DZ
ZH=A31*DX+A32*DY+A33*DZ
!-----------------------------------------------------------------------------------------------!
XP(1:4)=XX(1:4)
YP(1:4)=YY(1:4)
ZP(1:4)=ZZ(1:4)
XP(5)=XP(1)
YP(5)=YP(1)
ZP(5)=ZP(1)
!-----------------------------------------------------------------------------------------------!
!    Source Distribution                                                                        !
!-----------------------------------------------------------------------------------------------!
SOLUT=1.D0 !Different numeration of the corner points in comparison to Hess and Smith (1962)
VS=0.D0
!-----------------------------------------------------------------------------------------------!
VX=0.D0
VY=0.D0
VZ=0.D0
DO K=1,4
   K1=K+1
   dx =xh-xi(k)
   dy =yh-yi(k)
   dz =zh
   dx1=xh-xi(k1)
   dy1=yh-yi(k1)
   dz1=zh
!-----------------------------------------------------------------------------------------------!
   R12X= DY*DZ1-DZ*DY1
   R12Y=-DX*DZ1+DZ*DX1
   R12Z= DX*DY1-DY*DX1
!-----------------------------------------------------------------------------------------------!
   R12S2=R12X*R12X+R12Y*R12Y+R12Z*R12Z
   R1=DSQRT(DX *DX +DY *DY +DZ *DZ )
   R2=DSQRT(DX1*DX1+DY1*DY1+DZ1*DZ1)
   IF ((R1 < EPS).AND.(R2 < EPS).AND.(R12S2 < EPS)) flag=.false.
!-----------------------------------------------------------------------------------------------!
   E1=dx*dx+ZH*ZH
   E2=dx1*dx1+ZH*ZH
   H1=dx*dy
   H2=dx1*dy1
!-----------------------------------------------------------------------------------------------!
   DD=DSQRT((XI(K1)-XI(K))*(XI(K1)-XI(K))+(YI(K1)-YI(K))*(YI(K1)-YI(K)))
!  MM=(YI(K1)-YI(K))/(XI(K1)-XI(K))
   MM=SDIV((YI(K1)-YI(K)),(XI(K1)-XI(K)),0.D0)
!-----------------------------------------------------------------------------------------------!
   IF (DD > EPS) THEN
      ARG1=R1+R2-DD
      ARG2=R1+R2+DD
      IF (ARG1 > EPS) THEN
         VX=VX+(YI(K1)-YI(K))/DD*DLOG(ARG1/ARG2)
         VY=VY+(XI(K)-XI(K1))/DD*DLOG(ARG1/ARG2)
      END IF !(ARG1 > EPS)
!-----------------------------------------------------------------------------------------------!
      IF (DABS(MM) >= LIM) THEN
         ARG1=0.D0
         ARG2=0.D0
         VZ=VZ+ARG1-ARG2
      ELSEIF (DABS(ZH) <= EPS) THEN
         ARG1=DSIGN(0.5D0*PI,MM*E1-H1)
         ARG2=DSIGN(0.5D0*PI,MM*E2-H2)
         VZ=VZ+ARG1-ARG2
      ELSE
         ARG1=(MM*E1-H1)/(ZH*R1)
         ARG2=(MM*E2-H2)/(ZH*R2)
         VZ=VZ+DATAN(ARG1)-DATAN(ARG2)
      END IF
   END IF !(DD > EPS)
!-----------------------------------------------------------------------------------------------!
!  DD=DSQRT((XI(K1)-XI(K))**2+(YI(K1)-YI(K))**2)
!  MM=(YI(K1)-YI(K))/(XI(K1)-XI(K))
!  R1=DSQRT((XH-XI(K ))**2+(YH-YI(K ))**2+ZH**2)
!  R2=DSQRT((XH-XI(K1))**2+(YH-YI(K1))**2+ZH**2)
!  E1=(XH-XI(K ))**2+ZH**2
!  E2=(XH-XI(K1))**2+ZH**2
!  H1=(XH-XI(K ))*(YH-YI(K ))
!  H2=(XH-XI(K1))*(YH-YI(K1))
!-----------------------------------------------------------------------------------------------!
END DO !K=1,4
!-----------------------------------------------------------------------------------------------!
if (flag == .false.) then
   vx=0.d0
   vy=0.d0
   vz=2.d0*pi
end if
!-----------------------------------------------------------------------------------------------!
VS(1)=SOLUT*(A11*VX+A21*VY+A31*VZ)
VS(2)=SOLUT*(A12*VX+A22*VY+A32*VZ)
VS(3)=SOLUT*(A13*VX+A23*VY+A33*VZ)
!-----------------------------------------------------------------------------------------------!
!    Vortex Ring - Doublet Distribution                                                         !
!-----------------------------------------------------------------------------------------------!
GAMMA=-1.D0
VD=0.D0
!-----------------------------------------------------------------------------------------------!
DO K=1,4
   K1 =K+1
   DX =X-XP(K )
   DY =Y-YP(K )
   DZ =Z-ZP(K )
   DX1=X-XP(K1)
   DY1=Y-YP(K1)
   DZ1=Z-ZP(K1)
!-----------------------------------------------------------------------------------------------!
   R12X= DY*DZ1-DZ*DY1
   R12Y=-DX*DZ1+DZ*DX1
   R12Z= DX*DY1-DY*DX1
!-----------------------------------------------------------------------------------------------!
   R12S2=R12X*R12X+R12Y*R12Y+R12Z*R12Z
   R1=DSQRT(DX *DX +DY *DY +DZ *DZ )
   R2=DSQRT(DX1*DX1+DY1*DY1+DZ1*DZ1)
!-----------------------------------------------------------------------------------------------!
   IF ((R1 > EPS).AND.(R2 > EPS).AND.(R12S2 > EPS)) THEN
      R0R1=(XP(K1)-XP(K))*(X-XP(K ))+(YP(K1)-YP(K))*(Y-YP(K ))+(ZP(K1)-ZP(K))*(Z-ZP(K ))
      R0R2=(XP(K1)-XP(K))*(X-XP(K1))+(YP(K1)-YP(K))*(Y-YP(K1))+(ZP(K1)-ZP(K))*(Z-ZP(K1))
!-----------------------------------------------------------------------------------------------!
      KA=GAMMA/(R12S2)*(R0R1/R1-R0R2/R2)
      VD(1)=VD(1)+KA*R12X
      VD(2)=VD(2)+KA*R12Y
      VD(3)=VD(3)+KA*R12Z
   END IF !((R1 > EPS).AND.(R2 > EPS).AND.(R12S2 > EPS))
!-----------------------------------------------------------------------------------------------!
END DO !K=1,4
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE VELPAN
!===============================================================================================!
!    Image induced velocities                                                                   !
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
SUBROUTINE IMAGEVELO(N,XX,YY,ZZ,X,Y,Z,EPS,VD)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 03122013, J. Baltazar, version 1.0                                             !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
INTEGER :: K,K1,N
DOUBLE PRECISION :: X,Y,Z
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4),XP(5),YP(5),ZP(5),VD(3)
DOUBLE PRECISION :: DX,DY,DZ,DX1,DY1,DZ1,R12X,R12Y,R12Z,R12S2,R1,R2,R0R1,R0R2,KA,GAMMA,EPS
!-----------------------------------------------------------------------------------------------!
XP(1:4)=XX(1:4)
YP(1:4)=YY(1:4)
ZP(1:4)=ZZ(1:4)
XP(5)=XP(1)
YP(5)=YP(1)
ZP(5)=ZP(1)
!-----------------------------------------------------------------------------------------------!
GAMMA=1.D0
VD=0.D0
!-----------------------------------------------------------------------------------------------!
DO K=1,4,N
   K1 =K+1
   DX =X-XP(K )
   DY =Y-YP(K )
   DZ =Z-ZP(K )
   DX1=X-XP(K1)
   DY1=Y-YP(K1)
   DZ1=Z-ZP(K1)
!-----------------------------------------------------------------------------------------------!
   R12X= DY*DZ1-DZ*DY1
   R12Y=-DX*DZ1+DZ*DX1
   R12Z= DX*DY1-DY*DX1
!-----------------------------------------------------------------------------------------------!
   R12S2=R12X*R12X+R12Y*R12Y+R12Z*R12Z
   R1=DSQRT(DX *DX +DY *DY +DZ *DZ )
   R2=DSQRT(DX1*DX1+DY1*DY1+DZ1*DZ1)
!-----------------------------------------------------------------------------------------------!
   IF ((R1 > EPS).AND.(R2 > EPS).AND.(R12S2 > EPS)) THEN
      R0R1=(XP(K1)-XP(K))*(X-XP(K ))+(YP(K1)-YP(K))*(Y-YP(K ))+(ZP(K1)-ZP(K))*(Z-ZP(K ))
      R0R2=(XP(K1)-XP(K))*(X-XP(K1))+(YP(K1)-YP(K))*(Y-YP(K1))+(ZP(K1)-ZP(K))*(Z-ZP(K1))
!-----------------------------------------------------------------------------------------------!
      KA=GAMMA/(R12S2)*(R0R1/R1-R0R2/R2)
      VD(1)=VD(1)+KA*R12X
      VD(2)=VD(2)+KA*R12Y
      VD(3)=VD(3)+KA*R12Z
   END IF !((R1 > EPS).AND.(R2 > EPS).AND.(R12S2 > EPS))
!-----------------------------------------------------------------------------------------------!
END DO !K=1,4,N
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE IMAGEVELO
!-----------------------------------------------------------------------------------------------!
