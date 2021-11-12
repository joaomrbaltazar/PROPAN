!-----------------------------------------------------------------------------------------------!
!    Fourier Function                                                                           !
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
SUBROUTINE FOURIER_FUNCTION(NI,AN,BN,X,F)
!-----------------------------------------------------------------------------------------------!
!    Created by: Joao Baltazar, IST                                                             !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
INTEGER :: I,NI
DOUBLE PRECISION :: AN(0:NI),BN(1:NI),X,F
!-----------------------------------------------------------------------------------------------!
F=AN(0)*0.5D0
DO I=1,NI
   F=F+AN(I)*DCOS(DFLOAT(I)*X)+BN(I)*DSIN(DFLOAT(I)*X)
END DO !I=1,NI
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE FOURIER_FUNCTION
!-----------------------------------------------------------------------------------------------!
