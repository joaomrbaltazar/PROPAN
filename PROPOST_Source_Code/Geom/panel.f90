!-----------------------------------------------------------------------------------------------!
!    Subroutine PANEL                                                                           !
!                                                                                               !
!    Computation of the coordinates of element center, components of the unit normal and an     !
!    approximate area of a quadrilateral (hyperboloidal) panel. The area is exact for a         !
!    quadrilateral plane panel.                                                                 !
!                                                                                               !
!    Author: J.A.C. Falcao de Campos, IST, Jan 1995                                             !
!                                                                                               !
!    Copyright (C) 2021  J.A.C. Falcão de Campos                                                !
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
SUBROUTINE PANEL(XX,YY,ZZ,X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UN0X,UN0Y,UN0Z,A0)
!-----------------------------------------------------------------------------------------------!
!    INPUT:                                                                                     !
!                                                                                               !
!    XX(4),YY(4),ZZ(4):   CARTESIAN COORDINATES OF PANEL CORNER POINTS                          !
!                                                                                               !
!    OUTPUT:                                                                                    !
!                                                                                               !
!    X0,Y0,Z0         :   CARTESIAN COORDINATES OF PANEL CENTER POINT                           !
!                                                                                               !
!    A1X,A1Y,A1Z      :   COMPONENTS OF THE TANGENT VECTOR A1 AT THE CENTER POINT               !
!                                                                                               !
!    A2X,A2Y,A2Z      :   COMPONENTS OF THE TANGENT VECTOR A2 AT THE CENTER POINT               !
!                                                                                               !
!    UN0X,UN0Y,UN0Z   :   COMPONENTS OF THE UNIT NORMAL AT THE CENTER POINT                     !
!                                                                                               !
!    A0               :   PANEL AREA                                                            !
!-----------------------------------------------------------------------------------------------!
!    Declarations                                                                               !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4)
DOUBLE PRECISION :: X0,Y0,Z0,UN0X,UN0Y,UN0Z,A0
DOUBLE PRECISION :: A1X,A1Y,A1Z,A2X,A2Y,A2Z,VNX,VNY,VNZ,VNM
!-----------------------------------------------------------------------------------------------!
!    Center Point                                                                               !
!-----------------------------------------------------------------------------------------------!
X0=(XX(1)+XX(2)+XX(3)+XX(4))/4.D0
Y0=(YY(1)+YY(2)+YY(3)+YY(4))/4.D0
Z0=(ZZ(1)+ZZ(2)+ZZ(3)+ZZ(4))/4.D0
!-----------------------------------------------------------------------------------------------!
!    Tangent Vectors at the Center                                                              !
!-----------------------------------------------------------------------------------------------!
A1X=(XX(2)-XX(1)+XX(3)-XX(4))/4.D0
A1Y=(YY(2)-YY(1)+YY(3)-YY(4))/4.D0
A1Z=(ZZ(2)-ZZ(1)+ZZ(3)-ZZ(4))/4.D0
A2X=(XX(4)-XX(1)+XX(3)-XX(2))/4.D0
A2Y=(YY(4)-YY(1)+YY(3)-YY(2))/4.D0
A2Z=(ZZ(4)-ZZ(1)+ZZ(3)-ZZ(2))/4.D0
!-----------------------------------------------------------------------------------------------!
!    Unit Normal at the Center Point                                                            !
!-----------------------------------------------------------------------------------------------!
VNX=A1Y*A2Z-A1Z*A2Y
VNY=A1Z*A2X-A1X*A2Z
VNZ=A1X*A2Y-A1Y*A2X
VNM=DSQRT(VNX*VNX+VNY*VNY+VNZ*VNZ)
UN0X=VNX/VNM
UN0Y=VNY/VNM
UN0Z=VNZ/VNM
!-----------------------------------------------------------------------------------------------!
!    Panel Area                                                                                 !
!-----------------------------------------------------------------------------------------------!
A0=4.D0*VNM
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PANEL
!===============================================================================================!
!    Subroutine PANEL3                                                                          !
!                                                                                               !
!    Computation of the coordinates of element center, components of the unit normal of a       !
!    triangular panel.                                                                          !
!                                                                                               !
!    Author: Joao Baltazar, IST, Dez 2000                                                       !
!                                                                                               !
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
SUBROUTINE PANEL3(XX,YY,ZZ,X0,Y0,Z0,UN0X,UN0Y,UN0Z,A0)
!-----------------------------------------------------------------------------------------------!
!    INPUT:                                                                                     !
!                                                                                               !
!    XX(3),YY(3),ZZ(3):     CARTESIAN COORDINATES OF PANEL CORNER POINTS                        !
!                                                                                               !
!    OUTPUT:                                                                                    !
!                                                                                               !
!    X0,Y0,Z0         :     CARTESIAN COORDINATES OF PANEL CENTER POINT                         !
!                                                                                               !
!    UN0X,UN0Y,UN0Z   :     COMPONENTS OF THE UNIT NORMAL AT THE CENTER POINT                   !
!                                                                                               !
!    A0               :     PANEL AREA                                                          !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
DOUBLE PRECISION :: XX(3),YY(3),ZZ(3)
DOUBLE PRECISION :: X0,Y0,Z0,UN0X,UN0Y,UN0Z,A0
DOUBLE PRECISION :: A1X,A1Y,A1Z,A2X,A2Y,A2Z,VNX,VNY,VNZ,VNM
!-----------------------------------------------------------------------------------------------!
!    Center Point                                                                               !
!-----------------------------------------------------------------------------------------------!
X0=(XX(1)+XX(2)+XX(3))/3.D0
Y0=(YY(1)+YY(2)+YY(3))/3.D0
Z0=(ZZ(1)+ZZ(2)+ZZ(3))/3.D0
!-----------------------------------------------------------------------------------------------!
!    Tangent Vectors at the Center                                                              !
!-----------------------------------------------------------------------------------------------!
A1X=XX(2)-XX(1)
A1Y=YY(2)-YY(1)
A1Z=ZZ(2)-ZZ(1)
A2X=XX(3)-XX(1)
A2Y=YY(3)-YY(1)
A2Z=ZZ(3)-ZZ(1)
!-----------------------------------------------------------------------------------------------!
!    Unit Normal at the Center Point                                                            !
!-----------------------------------------------------------------------------------------------!
VNX=A1Y*A2Z-A1Z*A2Y
VNY=A1Z*A2X-A1X*A2Z
VNZ=A1X*A2Y-A1Y*A2X
VNM=DSQRT(VNX*VNX+VNY*VNY+VNZ*VNZ)
UN0X=VNX/VNM
UN0Y=VNY/VNM
UN0Z=VNZ/VNM
!-----------------------------------------------------------------------------------------------!
!    Panel Area                                                                                 !
!-----------------------------------------------------------------------------------------------!
A0=2.D0*VNM
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PANEL3
!===============================================================================================!
!    Computation of the coordinates of a flat panel                                             !
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
SUBROUTINE PANELFLAT(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0)
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4)
DOUBLE PRECISION :: X0,Y0,Z0,UNX0,UNY0,UNZ0,ZETA
!-----------------------------------------------------------------------------------------------!
!    For Flat Panel Redefine Corner Points                                                      !
!-----------------------------------------------------------------------------------------------!
ZETA=(XX(1)-X0)*UNX0+(YY(1)-Y0)*UNY0+(ZZ(1)-Z0)*UNZ0
XX(1)=XX(1)-ZETA*UNX0
YY(1)=YY(1)-ZETA*UNY0
ZZ(1)=ZZ(1)-ZETA*UNZ0
!-----------------------------------------------------------------------------------------------!
ZETA=(XX(2)-X0)*UNX0+(YY(2)-Y0)*UNY0+(ZZ(2)-Z0)*UNZ0
XX(2)=XX(2)-ZETA*UNX0
YY(2)=YY(2)-ZETA*UNY0
ZZ(2)=ZZ(2)-ZETA*UNZ0
!-----------------------------------------------------------------------------------------------!
ZETA=(XX(3)-X0)*UNX0+(YY(3)-Y0)*UNY0+(ZZ(3)-Z0)*UNZ0
XX(3)=XX(3)-ZETA*UNX0
YY(3)=YY(3)-ZETA*UNY0
ZZ(3)=ZZ(3)-ZETA*UNZ0
!-----------------------------------------------------------------------------------------------!
ZETA=(XX(4)-X0)*UNX0+(YY(4)-Y0)*UNY0+(ZZ(4)-Z0)*UNZ0
XX(4)=XX(4)-ZETA*UNX0
YY(4)=YY(4)-ZETA*UNY0
ZZ(4)=ZZ(4)-ZETA*UNZ0
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PANELFLAT
!-----------------------------------------------------------------------------------------------!
