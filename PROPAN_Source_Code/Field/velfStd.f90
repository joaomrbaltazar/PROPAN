!-----------------------------------------------------------------------------------------------!
!    Calculation of the Steady Velocity Field                                                   !
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
SUBROUTINE VELFSTD(CREF,JJ,TT,N1,N2,XTMP,YTMP,ZTMP,UTMP,VTMP,WTMP)
!-----------------------------------------------------------------------------------------------!
!    Created by: 09052016, J. Baltazar, version 1.2                                             !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
CHARACTER*8 CREF
INTEGER :: I,J,JJ,TT,T1,N1,N2
DOUBLE PRECISION :: TETA
DOUBLE PRECISION :: XTMP(N1,N2),YTMP(N1,N2),ZTMP(N1,N2)
DOUBLE PRECISION :: UTMP(N1,N2),VTMP(N1,N2),WTMP(N1,N2)
DOUBLE PRECISION :: UNEW(N1,N2),VNEW(N1,N2),WNEW(N1,N2)
!-----------------------------------------------------------------------------------------------!
!    Undisturbed Velocity                                                                       !
!-----------------------------------------------------------------------------------------------!
CALL VELINF(CREF,JJ,TT,N1,N2,XTMP,YTMP,ZTMP,UNEW,VNEW,WNEW)
UTMP=UNEW
VTMP=VNEW
WTMP=WNEW
!-----------------------------------------------------------------------------------------------!
!    Blade Velocities                                                                           !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   CALL BLADEVELO(TT,N1,N2,XTMP,YTMP,ZTMP,UNEW,VNEW,WNEW)
   UTMP=UTMP+UNEW
   VTMP=VTMP+VNEW
   WTMP=WTMP+WNEW
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake Velocities                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   CALL BLADEWAKEVELO(TT,N1,N2,XTMP,YTMP,ZTMP,UNEW,VNEW,WNEW)
   UTMP=UTMP+UNEW
   VTMP=VTMP+VNEW
   WTMP=WTMP+WNEW
END IF !(IP ==1 )
!-----------------------------------------------------------------------------------------------!
!    Nozzle Velocities                                                                          !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   CALL NOZZLEVELO(TT,N1,N2,XTMP,YTMP,ZTMP,UNEW,VNEW,WNEW)
   UTMP=UTMP+UNEW
   VTMP=VTMP+VNEW
   WTMP=WTMP+WNEW
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake Velocities                                                                     !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   CALL NOZZLEWAKEVELO(TT,N1,N2,XTMP,YTMP,ZTMP,UNEW,VNEW,WNEW)
   UTMP=UTMP+UNEW
   VTMP=VTMP+VNEW
   WTMP=WTMP+WNEW
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub Velocities                                                                             !
!-----------------------------------------------------------------------------------------------!
IF (IABS(IH) == 1) THEN
   IF (NWA == 0) CALL      HUBVELO(TT,N1,N2,XTMP,YTMP,ZTMP,UNEW,VNEW,WNEW)
   IF (NWA /= 0) CALL IMAGEHUBVELO(TT,N1,N2,XTMP,YTMP,ZTMP,UNEW,VNEW,WNEW)
   UTMP=UTMP+UNEW
   VTMP=VTMP+VNEW
   WTMP=WTMP+WNEW
END IF !(IABS(IH) == 1)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE VELFSTD
!-----------------------------------------------------------------------------------------------!
