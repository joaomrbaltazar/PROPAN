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
SUBROUTINE NOZZLEGEOM
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
!    Loop on the Nozzle Panels                                                                  !
!-----------------------------------------------------------------------------------------------!
!    First Half Sector                                                                          !
!-----------------------------------------------------------------------------------------------!
DO J=1,NNT
   DO I=1,NNXT1
      XX(1)=XN(I  ,J  )
      XX(2)=XN(I+1,J  )
      XX(3)=XN(I+1,J+1)
      XX(4)=XN(I  ,J+1)
      YY(1)=YN(I  ,J  )
      YY(2)=YN(I+1,J  )
      YY(3)=YN(I+1,J+1)
      YY(4)=YN(I  ,J+1)
      ZZ(1)=ZN(I  ,J  )
      ZZ(2)=ZN(I+1,J  )
      ZZ(3)=ZN(I+1,J+1)
      ZZ(4)=ZN(I  ,J+1)
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
      CALL PANEL(XX,YY,ZZ,XN0(I,J),YN0(I,J),ZN0(I,J),AT1XN(I,J),AT1YN(I,J),AT1ZN(I,J), &
                 AT2XN(I,J),AT2YN(I,J),AT2ZN(I,J),UNXN0(I,J),UNYN0(I,J),UNZN0(I,J),AN0(I,J))
   END DO !I=1,NNXT1
END DO !J=1,NNT
!-----------------------------------------------------------------------------------------------!
!    Second Half Sector                                                                         !
!-----------------------------------------------------------------------------------------------!
DO J=NNT2,NNTT1
   DO I=1,NNXT1
      XX(1)=XN(I  ,J  )
      XX(2)=XN(I+1,J  )
      XX(3)=XN(I+1,J+1)
      XX(4)=XN(I  ,J+1)
      YY(1)=YN(I  ,J  )
      YY(2)=YN(I+1,J  )
      YY(3)=YN(I+1,J+1)
      YY(4)=YN(I  ,J+1)
      ZZ(1)=ZN(I  ,J  )
      ZZ(2)=ZN(I+1,J  )
      ZZ(3)=ZN(I+1,J+1)
      ZZ(4)=ZN(I  ,J+1)
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
      CALL PANEL(XX,YY,ZZ,XN0(I,J-1),YN0(I,J-1),ZN0(I,J-1), &
                 AT1XN(I,J-1),AT1YN(I,J-1),AT1ZN(I,J-1), &
                 AT2XN(I,J-1),AT2YN(I,J-1),AT2ZN(I,J-1), &
                 UNXN0(I,J-1),UNYN0(I,J-1),UNZN0(I,J-1),AN0(I,J-1))
   END DO !I=1,NNXT1
END DO !J=NNT2,NNTT1
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE NOZZLEGEOM
!-----------------------------------------------------------------------------------------------!
