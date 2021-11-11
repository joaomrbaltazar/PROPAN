!-----------------------------------------------------------------------------------------------!
!    Periodic Flow                                                                              !
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
SUBROUTINE PERIODICFLOW(TT,KB,N1,N2,N3,VART,VAR)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 25112013, J. Baltazar, version 1.0                                             !
!    Modified  : 03062016, J. Baltazar, 2016 version 1.3                                        !
!    Modified  : 23122020, J. Baltazar, 2020 version 1.2                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
INTEGER :: TT,KB,N1,N2,N3
DOUBLE PRECISION :: VART(N1,N2,0:N3),VAR(N1,N2)
!-----------------------------------------------------------------------------------------------!
IF (KB == 1) THEN
   VAR(:,:)=VART(:,:,TT)
ELSEIF (IFILESOL == 1) THEN
   VAR(:,:)=VART(:,:,TT)
ELSEIF (TT < IDINT(P)*(KB-1)) THEN
   VAR(:,:)=VART(:,:,0)
ELSE !(TT)
   VAR(:,:)=VART(:,:,(TT-IDINT(P)*(KB-1)))
END IF !(TT)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PERIODICFLOW
!-----------------------------------------------------------------------------------------------!
