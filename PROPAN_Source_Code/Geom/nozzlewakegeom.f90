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
SUBROUTINE NOZZLEWAKEGEOM
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
!    Loop on the Nozzle Wake Panels                                                             !
!-----------------------------------------------------------------------------------------------!
!    First Half Sector                                                                          !
!-----------------------------------------------------------------------------------------------!
DO J=1,NNT
   DO I=1,NNW
      XX(1)=XNW(I  ,J  )
      XX(2)=XNW(I+1,J  )
      XX(3)=XNW(I+1,J+1)
      XX(4)=XNW(I  ,J+1)
      YY(1)=YNW(I  ,J  )
      YY(2)=YNW(I+1,J  )
      YY(3)=YNW(I+1,J+1)
      YY(4)=YNW(I  ,J+1)
      ZZ(1)=ZNW(I  ,J  )
      ZZ(2)=ZNW(I+1,J  )
      ZZ(3)=ZNW(I+1,J+1)
      ZZ(4)=ZNW(I  ,J+1)
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
      CALL PANEL(XX,YY,ZZ,XNW0(I,J),YNW0(I,J),ZNW0(I,J), &
                    AT1XNW(I,J),AT1YNW(I,J),AT1ZNW(I,J), &
                    AT2XNW(I,J),AT2YNW(I,J),AT2ZNW(I,J), &
                    UNXNW0(I,J),UNYNW0(I,J),UNZNW0(I,J),ANW0(I,J))
   END DO !I=1,NNW
END DO !J=1,NNT
!-----------------------------------------------------------------------------------------------!
!    Second Half Sector                                                                         !
!-----------------------------------------------------------------------------------------------!
DO J=NNT2,NNTT1
   DO I=1,NNW
      XX(1)=XNW(I  ,J  )
      XX(2)=XNW(I+1,J  )
      XX(3)=XNW(I+1,J+1)
      XX(4)=XNW(I  ,J+1)
      YY(1)=YNW(I  ,J  )
      YY(2)=YNW(I+1,J  )
      YY(3)=YNW(I+1,J+1)
      YY(4)=YNW(I  ,J+1)
      ZZ(1)=ZNW(I  ,J  )
      ZZ(2)=ZNW(I+1,J  )
      ZZ(3)=ZNW(I+1,J+1)
      ZZ(4)=ZNW(I  ,J+1)
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
      CALL PANEL(XX,YY,ZZ,XNW0(I,J-1),YNW0(I,J-1),ZNW0(I,J-1), &
                    AT1XNW(I,J-1),AT1YNW(I,J-1),AT1ZNW(I,J-1), &
                    AT2XNW(I,J-1),AT2YNW(I,J-1),AT2ZNW(I,J-1), &
                    UNXNW0(I,J-1),UNYNW0(I,J-1),UNZNW0(I,J-1),ANW0(I,J-1))
   END DO !I=1,NNW
END DO !J=NNT2,NNTT1
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE NOZZLEWAKEGEOM
!-----------------------------------------------------------------------------------------------!
