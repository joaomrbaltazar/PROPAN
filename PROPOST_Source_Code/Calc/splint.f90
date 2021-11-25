!-----------------------------------------------------------------------------------------------!
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
SUBROUTINE SPLINT(NDATA,X,Z,N,XX,ZZ)
!-----------------------------------------------------------------------------------------------!
!    Author: 07060213, J. Baltazar                                                              !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
INTEGER :: I,NDATA,N
DOUBLE PRECISION :: X(NDATA),Z(NDATA),XX(N),ZZ(N)
DOUBLE PRECISION :: B(NDATA),C(NDATA),D(NDATA),ISPLINE
!-----------------------------------------------------------------------------------------------!
!    Call Spline to Calculate Spline Coeficients                                                !
!-----------------------------------------------------------------------------------------------!
CALL SPLINE(X,Z,B,C,D,NDATA)
!-----------------------------------------------------------------------------------------------!
DO I=1,N
   ZZ(I)=ISPLINE(XX(I),X,Z,B,C,D,NDATA)
END DO !I=1,N
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE SPLINT
!-----------------------------------------------------------------------------------------------!
