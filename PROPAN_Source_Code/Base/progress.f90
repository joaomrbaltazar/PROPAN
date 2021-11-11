!-----------------------------------------------------------------------------------------------!
!    Progress Bar                                                                               !
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
SUBROUTINE PROGRESS(J)
!-----------------------------------------------------------------------------------------------!
!    Modified  : 11102013, J. Baltazar, version 1.0                                             !
!                06112013, J. Baltazar, include carriagecontrol                                 !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
INTEGER      :: J,K
CHARACTER*68 :: BAR=' ???% [                                                            ]'
!-----------------------------------------------------------------------------------------------!
WRITE(UNIT=BAR(2:4),FMT='(I3)') J
IF (J == 100) THEN
   DO K=1,NINT(J*0.6)
      BAR(7+K:7+K)='='
   END DO !K=1,NINT(J*0.6)
ELSEIF (J > 0) THEN
   DO K=1,NINT(J*0.6)-1
      BAR(7+K:7+K)='='
   END DO !K=1,NINT(J*0.6)-1
   BAR(7+K:7+K)='>'
END IF !(J)
!-----------------------------------------------------------------------------------------------!
!    Print the Progress Bar                                                                     !
!-----------------------------------------------------------------------------------------------!
OPEN(UNIT=6,CARRIAGECONTROL='FORTRAN')
IF (J == 100) THEN
   WRITE(UNIT=6,FMT='(A1,A1,A68)',ADVANCE='YES') ' ',CHAR(13),BAR
ELSE
   WRITE(UNIT=6,FMT='(A1,A1,A68)',ADVANCE='NO' ) ' ',CHAR(13),BAR
END IF !(J)
CLOSE(UNIT=6)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PROGRESS
!-----------------------------------------------------------------------------------------------!
