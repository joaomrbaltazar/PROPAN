!-----------------------------------------------------------------------------------------------!
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
SUBROUTINE GAUSSVEL(N,XX,YY,ZZ,X,Y,Z,VS,VD)
!-----------------------------------------------------------------------------------------------!
!    INPUT:                                                                                     !
!                                                                                               !
!    N                    :      ORDER OF GAUSSIAN RULE FOR PANEL INTEGRATION                   !
!                                                                                               !
!    XX(4),YY(4),ZZ(4)    :      CARTESIAN COORDINATES OF PANEL CORNER POINTS                   !
!                                                                                               !
!    X,Y,Z                :      CARTESIAN COORDINATES OF FIELD POINT                           !
!                                                                                               !
!    OUTPUT:                                                                                    !
!                                                                                               !
!    VS(3)                :      CARTESIAN VELOCITY COMPONENTS INDUCED BY                       !
!                                SOURCE DISTRIBUTION ON THE ELEMENT OF STRENGTH -4*PI           !
!                                                                                               !
!    VD(3)                :      CARTESIAN VELOCITY COMPONENTS INDUCED BY                       !
!                                UNIT DIPOLE DISTRIBUTION ON THE ELEMENT OF STRENGHT -4*PI      !
!-----------------------------------------------------------------------------------------------!
!    Declarations                                                                               !
!-----------------------------------------------------------------------------------------------!
INTEGER :: I,J,N
DOUBLE PRECISION XG(N),W(N)
DOUBLE PRECISION CSI,ETA
DOUBLE PRECISION XX(4),YY(4),ZZ(4),X,Y,Z
DOUBLE PRECISION VS(3),VD(3)
DOUBLE PRECISION RM,A1X,A1Y,A1Z,A2X,A2Y,A2Z,VNX,VNY,VNZ,VNM,RDVN
DOUBLE PRECISION FS(3),FD(3),SUMS(3),SUMD(3),XQ,YQ,ZQ
!-----------------------------------------------------------------------------------------------!
!    Computations                                                                               !
!-----------------------------------------------------------------------------------------------!
CALL GAUSS2M(N,XG,W)
!-----------------------------------------------------------------------------------------------!
!    Loop on the Element                                                                        !
!-----------------------------------------------------------------------------------------------!
SUMS=0.D0
SUMD=0.D0
DO I=1,N
   CSI=XG(I)
   DO J=1,N
      ETA=XG(J)
      CALL PANCOOR(XX,YY,ZZ,CSI,ETA,XQ,YQ,ZQ)
!-----------------------------------------------------------------------------------------------!
!    Tangent vectors                                                                            !
!-----------------------------------------------------------------------------------------------!
      A1X=0.25D0*(XX(2)-XX(1)+XX(3)-XX(4))+0.25D0*ETA*(XX(3)-XX(2)+XX(1)-XX(4))
      A1Y=0.25D0*(YY(2)-YY(1)+YY(3)-YY(4))+0.25D0*ETA*(YY(3)-YY(2)+YY(1)-YY(4))
      A1Z=0.25D0*(ZZ(2)-ZZ(1)+ZZ(3)-ZZ(4))+0.25D0*ETA*(ZZ(3)-ZZ(2)+ZZ(1)-ZZ(4))
      A2X=0.25D0*(XX(4)-XX(1)+XX(3)-XX(2))+0.25D0*CSI*(XX(3)-XX(2)+XX(1)-XX(4))
      A2Y=0.25D0*(YY(4)-YY(1)+YY(3)-YY(2))+0.25D0*CSI*(YY(3)-YY(2)+YY(1)-YY(4))
      A2Z=0.25D0*(ZZ(4)-ZZ(1)+ZZ(3)-ZZ(2))+0.25D0*CSI*(ZZ(3)-ZZ(2)+ZZ(1)-ZZ(4))
!-----------------------------------------------------------------------------------------------!
!    Normal vector                                                                              !
!-----------------------------------------------------------------------------------------------!
      VNX=A1Y*A2Z-A1Z*A2Y
      VNY=A1Z*A2X-A1X*A2Z
      VNZ=A1X*A2Y-A1Y*A2X
!-----------------------------------------------------------------------------------------------!
!    Jacobian                                                                                   !
!-----------------------------------------------------------------------------------------------!
      VNM=DSQRT(VNX*VNX+VNY*VNY+VNZ*VNZ)
!-----------------------------------------------------------------------------------------------!
!    Dot product                                                                                !
!-----------------------------------------------------------------------------------------------!
      RDVN=VNX*(X-XQ)+VNY*(Y-YQ)+VNZ*(Z-ZQ)
!-----------------------------------------------------------------------------------------------!
!    Cartesian distance                                                                         !
!-----------------------------------------------------------------------------------------------!
      RM=DSQRT((X-XQ)**2+(Y-YQ)**2+(Z-ZQ)**2)
!-----------------------------------------------------------------------------------------------!
!    Source                                                                                     !
!-----------------------------------------------------------------------------------------------!
      FS(1)=(X-XQ)*VNM/RM**3
      FS(2)=(Y-YQ)*VNM/RM**3
      FS(3)=(Z-ZQ)*VNM/RM**3
!-----------------------------------------------------------------------------------------------!
!    Dipole                                                                                     !
!-----------------------------------------------------------------------------------------------!
      FD(1)=(3.D0*RDVN*(X-XQ)/RM**2-VNX)/RM**3
      FD(2)=(3.D0*RDVN*(Y-YQ)/RM**2-VNY)/RM**3
      FD(3)=(3.D0*RDVN*(Z-ZQ)/RM**2-VNZ)/RM**3
!-----------------------------------------------------------------------------------------------!
      SUMS(1:3)=SUMS(1:3)+W(I)*W(J)*FS(1:3)
      SUMD(1:3)=SUMD(1:3)+W(I)*W(J)*FD(1:3)
   END DO !J=1,N
END DO !I=1,N
!-----------------------------------------------------------------------------------------------!
!    Source                                                                                     !
!-----------------------------------------------------------------------------------------------!
VS=-SUMS
!-----------------------------------------------------------------------------------------------!
!    Dipole                                                                                     !
!-----------------------------------------------------------------------------------------------!
VD=-SUMD
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE GAUSSVEL
!-----------------------------------------------------------------------------------------------!
