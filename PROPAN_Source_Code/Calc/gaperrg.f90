!-----------------------------------------------------------------------------------------------!
!    Write Error in Gap Model                                                                   !
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
SUBROUTINE GAPERRG(TT,ERRG)
!-----------------------------------------------------------------------------------------------!
!    Created by: 19112013, J. Baltazar, version 1.0                                             !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,TT
DOUBLE PRECISION :: ERRG
!-----------------------------------------------------------------------------------------------!
ERRG=0.D0
DO I=1,NCP
   WORK1=CPP(I,NRP,TT)-CPP(NCP-I+1,NRP,TT)
   DCPG(I)=DABS(DCPG(I)-WORK1)
   IF (DABS(DCPG(I)) > ERRG) ERRG=DABS(DCPG(I))
END DO !I=1,NCP
DO I=1,NCP
   DCPG(I)=CPP(I,NRP,TT)-CPP(NCP-I+1,NRP,TT)
END DO !I=1,NCP
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE GAPERRG
!-----------------------------------------------------------------------------------------------!
