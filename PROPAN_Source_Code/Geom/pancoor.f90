!-----------------------------------------------------------------------------------------------!
!    Computes the cartesian coordinates of a point on a hyperboloidal panel                     !
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
SUBROUTINE PANCOOR(XX,YY,ZZ,CSI,ETA,XQ,YQ,ZQ)
!-----------------------------------------------------------------------------------------------!
!    INPUT:                                                                                     !
!                                                                                               !
!    XX(4),YY(4),ZZ(4):   CARTESIAN COORDINATES OF PANEL CORNER POINTS                          !
!                                                                                               !
!    CSI, ETA         :   PARAMETRIC COORDINATES ON THE PANEL                                   !
!                                                                                               !
!    OUTPUT:                                                                                    !
!                                                                                               !
!    XQ,YQ,ZQ         :   CARTESIAN COORDINATE OF THE POINT                                     !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4),CSI,ETA,XQ,YQ,ZQ,PSI1,PSI2,PSI3,PSI4
!-----------------------------------------------------------------------------------------------!
PSI1=0.25D0*(1.D0-CSI)*(1.D0-ETA)
PSI2=0.25D0*(1.D0+CSI)*(1.D0-ETA)
PSI3=0.25D0*(1.D0+CSI)*(1.D0+ETA)
PSI4=0.25D0*(1.D0-CSI)*(1.D0+ETA)
!-----------------------------------------------------------------------------------------------!
XQ=PSI1*XX(1)+PSI2*XX(2)+PSI3*XX(3)+PSI4*XX(4)
YQ=PSI1*YY(1)+PSI2*YY(2)+PSI3*YY(3)+PSI4*YY(4)
ZQ=PSI1*ZZ(1)+PSI2*ZZ(2)+PSI3*ZZ(3)+PSI4*ZZ(4)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PANCOOR
!===============================================================================================!
!    Computes the cartesian coordinates of a point on a triangular panel                        !
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
SUBROUTINE PANCOOR3(XX,YY,ZZ,CSI,ETA,XQ,YQ,ZQ)
!-----------------------------------------------------------------------------------------------!
!    INPUT:                                                                                     !
!                                                                                               !
!    XX(3),YY(3),ZZ(3):      CARTESIAN COORDINATES OF PANEL CORNER POINTS                       !
!                                                                                               !
!    CSI, ETA         :      PARAMETRIC COORDINATES ON THE PANEL                                !
!                                                                                               !
!    OUTPUT:                                                                                    !
!                                                                                               !
!    XQ,YQ,ZQ         :      CARTESIAN COORDINATE OF THE POINT                                  !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
DOUBLE PRECISION :: XX(3),YY(3),ZZ(3),CSI,ETA,XQ,YQ,ZQ,PSI1,PSI2,PSI3
!-----------------------------------------------------------------------------------------------!
PSI1=1.D0-CSI-ETA
PSI2=CSI
PSI3=ETA
!-----------------------------------------------------------------------------------------------!
XQ=PSI1*XX(1)+PSI2*XX(2)+PSI3*XX(3)
YQ=PSI1*YY(1)+PSI2*YY(2)+PSI3*YY(3)
ZQ=PSI1*ZZ(1)+PSI2*ZZ(2)+PSI3*ZZ(3)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PANCOOR3
!===============================================================================================!
!    Compute tangent vectors on hyperboloidal panel                                             !
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
SUBROUTINE PANVECTORS(XX,YY,ZZ,CSI,ETA,A1X,A1Y,A1Z,A2X,A2Y,A2Z)
!-----------------------------------------------------------------------------------------------!
!    Created by: Joao Baltazar, IST, October 2004                                               !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4),CSI,ETA
DOUBLE PRECISION :: A1X,A1Y,A1Z,A2X,A2Y,A2Z
!-----------------------------------------------------------------------------------------------!
A1X=0.25D0*(XX(2)-XX(1)+XX(3)-XX(4))+0.25D0*ETA*(XX(3)-XX(2)+XX(1)-XX(4))
A1Y=0.25D0*(YY(2)-YY(1)+YY(3)-YY(4))+0.25D0*ETA*(YY(3)-YY(2)+YY(1)-YY(4))
A1Z=0.25D0*(ZZ(2)-ZZ(1)+ZZ(3)-ZZ(4))+0.25D0*ETA*(ZZ(3)-ZZ(2)+ZZ(1)-ZZ(4))
A2X=0.25D0*(XX(4)-XX(1)+XX(3)-XX(2))+0.25D0*CSI*(XX(3)-XX(2)+XX(1)-XX(4))
A2Y=0.25D0*(YY(4)-YY(1)+YY(3)-YY(2))+0.25D0*CSI*(YY(3)-YY(2)+YY(1)-YY(4))
A2Z=0.25D0*(ZZ(4)-ZZ(1)+ZZ(3)-ZZ(2))+0.25D0*CSI*(ZZ(3)-ZZ(2)+ZZ(1)-ZZ(4))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PANVECTORS
!-----------------------------------------------------------------------------------------------!
