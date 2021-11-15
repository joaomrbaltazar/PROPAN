!-----------------------------------------------------------------------------------------------!
!    Unsteady wake model                                                                        !
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
SUBROUTINE SOLVEWAKE(TT,CC)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 25112013, J. Baltazar, version 1.0                                             !
!    Modified  : 04122014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,TT,CC
!-----------------------------------------------------------------------------------------------!
!    Dipole Strength in the Blade Wake                                                          !
!-----------------------------------------------------------------------------------------------!
DO I=NPW,2,-1
   POTPW(I,:,TT)=POTPW(I-1,:,TT-1)
END DO !I=NPW,2,-1
POTPW(1,:,TT)=0.D0 !Unknown
!-----------------------------------------------------------------------------------------------!
IF ((IRED == 1).AND.(CC == 0)) THEN
   DO I=NPW,2,-1
      POTPWWET(I,:,TT)=POTPWWET(I-1,:,TT-1)
   END DO !I=NPW,2,-1
   POTPWWET(1,:,TT)=0.D0 !Unknown
END IF !((IRED == 1).AND.(CC == 0))
!-----------------------------------------------------------------------------------------------!
!    Dipole Strength in the Nozzle Wake                                                         !
!-----------------------------------------------------------------------------------------------!
DO I=NNW,2,-1
   POTNW(I,:,TT)=POTNW(I-1,:,TT-1)
END DO !I=NNW,2,-1
POTNW(1,:,TT)=0.D0 !Unknown
!-----------------------------------------------------------------------------------------------!
IF ((IRED == 1).AND.(CC == 0)) THEN
   DO I=NNW,2,-1
      POTNWWET(I,:,TT)=POTNWWET(I-1,:,TT-1)
   END DO !I=NNW,2,-1
   POTNWWET(1,:,TT)=0.D0 !Unknown
END IF !((IRED == 1).AND.(CC == 0))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE SOLVEWAKE
!-----------------------------------------------------------------------------------------------!
