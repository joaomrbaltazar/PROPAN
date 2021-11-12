!-----------------------------------------------------------------------------------------------!
!    Compute Centroid Unit Normal and Area                                                      !
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
SUBROUTINE BLADEGEOM
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 14112013, J. Baltazar, version 1.0                                             !
!    Modified  : 30062017, J. Baltazar, 2017 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4),XA(NCP1,NRP:NRP1),YA(NCP1,NRP:NRP1),ZA(NCP1,NRP:NRP1)
DOUBLE PRECISION :: X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0
!-----------------------------------------------------------------------------------------------!
!    Loop on the Blade Panels                                                                   !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRP
   DO I=1,NCP
      XX(1)=XP(I  ,J  )
      XX(2)=XP(I+1,J  )
      XX(3)=XP(I+1,J+1)
      XX(4)=XP(I  ,J+1)
      YY(1)=YP(I  ,J  )
      YY(2)=YP(I+1,J  )
      YY(3)=YP(I+1,J+1)
      YY(4)=YP(I  ,J+1)
      ZZ(1)=ZP(I  ,J  )
      ZZ(2)=ZP(I+1,J  )
      ZZ(3)=ZP(I+1,J+1)
      ZZ(4)=ZP(I  ,J+1)
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
      CALL PANEL(XX,YY,ZZ,XP0(I,J),YP0(I,J),ZP0(I,J),AT1XP(I,J),AT1YP(I,J),AT1ZP(I,J), &
                 AT2XP(I,J),AT2YP(I,J),AT2ZP(I,J),UNXP0(I,J),UNYP0(I,J),UNZP0(I,J),AP0(I,J))
   END DO !I=1,NCP
END DO !I=1,NRP
!-----------------------------------------------------------------------------------------------!
!    Loop on the Blade Gap Panels                                                               !
!-----------------------------------------------------------------------------------------------!
IF (ISTRIP == 1) THEN
   DO J=NRP,NRP1
      DO I=1,NC1
         XA(NC1-I+1,J)=0.5D0*(XP(NC1-I+1,J)+XP(NC1+I-1,J))
         YA(NC1-I+1,J)=0.5D0*(YP(NC1-I+1,J)+YP(NC1+I-1,J))
         ZA(NC1-I+1,J)=0.5D0*(ZP(NC1-I+1,J)+ZP(NC1+I-1,J))
         XA(NC1+I-1,J)=XA(NC1-I+1,J)
         YA(NC1+I-1,J)=YA(NC1-I+1,J)
         ZA(NC1+I-1,J)=ZA(NC1-I+1,J)
      END DO !I=1,NC1
   END DO !J=NRP,NRP1
   DO I=1,NCP
      XX(1)=XA(I  ,NRP )
      XX(2)=XA(I+1,NRP )
      XX(3)=XA(I+1,NRP1)
      XX(4)=XA(I  ,NRP1)
      YY(1)=YA(I  ,NRP )
      YY(2)=YA(I+1,NRP )
      YY(3)=YA(I+1,NRP1)
      YY(4)=YA(I  ,NRP1)
      ZZ(1)=ZA(I  ,NRP )
      ZZ(2)=ZA(I+1,NRP )
      ZZ(3)=ZA(I+1,NRP1)
      ZZ(4)=ZA(I  ,NRP1)
!-----------------------------------------------------------------------------------------------!
!    Compute Panel Centroid Data                                                                !
!-----------------------------------------------------------------------------------------------!
      CALL PANEL(XX,YY,ZZ,X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNXC0(I),UNYC0(I),UNZC0(I),A0)
   END DO !I=1,NCP
END IF !(ISTRIP == 1)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE BLADEGEOM
!-----------------------------------------------------------------------------------------------!
