!-----------------------------------------------------------------------------------------------!
!    Compute Centroid Unit Normal and Area                                                      !
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
SUBROUTINE BLADEWAKEGEOM
!-----------------------------------------------------------------------------------------------!
!    Created by: 04112013, J. Baltazar, IST                                                     !
!    Modified  : 04112013, J. Baltazar, version 1.0                                             !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4)
DOUBLE PRECISION :: X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0
!-----------------------------------------------------------------------------------------------!
!    Loop on the Blade Wake Panels                                                              !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRW
   DO I=1,NPW
      XX(1)=XPW(I  ,J  )
      XX(2)=XPW(I+1,J  )
      XX(3)=XPW(I+1,J+1)
      XX(4)=XPW(I  ,J+1)
      YY(1)=YPW(I  ,J  )
      YY(2)=YPW(I+1,J  )
      YY(3)=YPW(I+1,J+1)
      YY(4)=YPW(I  ,J+1)
      ZZ(1)=ZPW(I  ,J  )
      ZZ(2)=ZPW(I+1,J  )
      ZZ(3)=ZPW(I+1,J+1)
      ZZ(4)=ZPW(I  ,J+1)
!-----------------------------------------------------------------------------------------------!
!    For Flat Panel Redefine Corner Points                                                      !
!-----------------------------------------------------------------------------------------------!
      IF (IPAN == 0) THEN
         CALL PANEL(XX,YY,ZZ,X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0)
         CALL PANELFLAT(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0)
      END IF !(IPAN == 0)
!-----------------------------------------------------------------------------------------------!
!    Compute Panel Centroid Data                                                                !
!-----------------------------------------------------------------------------------------------!
      CALL PANEL(XX,YY,ZZ,XPW0(I,J),YPW0(I,J),ZPW0(I,J),AT1XPW(I,J),AT1YPW(I,J),AT1ZPW(I,J), &
              AT2XPW(I,J),AT2YPW(I,J),AT2ZPW(I,J),UNXPW0(I,J),UNYPW0(I,J),UNZPW0(I,J),APW0(I,J))
   END DO !I=1,NPW
END DO !J=1,NRW
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE BLADEWAKEGEOM
!-----------------------------------------------------------------------------------------------!
