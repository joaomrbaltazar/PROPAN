!-----------------------------------------------------------------------------------------------!
!    Periodic flow                                                                              !
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
SUBROUTINE PERIODICFLOW(TT,KB,N1,N2,N3,VART,VAR)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 23052016, J. Baltazar, 2016 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPOST_MOD
INTEGER :: TT,KB,N1,N2,N3
DOUBLE PRECISION :: VART(N1,N2,0:N3),VAR(N1,N2)
!-----------------------------------------------------------------------------------------------!
IF (TT == 0) THEN
   VAR(:,:)=VART(:,:,0)
ELSEIF (TT < IDINT(P)*(KB-1)) THEN
   VAR(:,:)=VART(:,:,0)
ELSE
   VAR(:,:)=VART(:,:,(TT-IDINT(P)*(KB-1)))
END IF !(TT < IDINT(P)*(KB-1))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PERIODICFLOW
!-----------------------------------------------------------------------------------------------!
