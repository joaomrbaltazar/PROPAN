!-----------------------------------------------------------------------------------------------!
!    Calculation of the Steady Pressure Field                                                   !
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
SUBROUTINE PRESFSTD(CREF,JJ,TT,N1,N2,XTMP,YTMP,ZTMP,DPOTDT)
!-----------------------------------------------------------------------------------------------!
!    Created by: 09052016, J. Baltazar, version 1.2                                             !
!    Modified  : 26102016, J. Baltazar, version 1.4                                             !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
CHARACTER*8 CREF
INTEGER :: I,J,JJ,TT,N1,N2
DOUBLE PRECISION :: XTMP(N1,N2),YTMP(N1,N2),ZTMP(N1,N2)
DOUBLE PRECISION :: PTMP(N1,N2),DPOTDT(N1,N2)
!-----------------------------------------------------------------------------------------------!
!    Blade Potential                                                                            !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   CALL BLADECOEFF(TT,N1,N2,XTMP,YTMP,ZTMP,PTMP)
   POTF(:,:,TT)=PTMP(:,:)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake Potential                                                                       !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   CALL BLADEWAKECOEFF(TT,N1,N2,XTMP,YTMP,ZTMP,PTMP)
   POTF(:,:,TT)=POTF(:,:,TT)+PTMP(:,:)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Potential                                                                           !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   CALL NOZZLECOEFF(TT,N1,N2,XTMP,YTMP,ZTMP,PTMP)
   POTF(:,:,TT)=POTF(:,:,TT)+PTMP(:,:)
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake Potential                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   CALL NOZZLEWAKECOEFF(TT,N1,N2,XTMP,YTMP,ZTMP,PTMP)
   POTF(:,:,TT)=POTF(:,:,TT)+PTMP(:,:)
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub Potential                                                                              !
!-----------------------------------------------------------------------------------------------!
IF (IABS(IH) == 1) THEN
   CALL HUBCOEFF(TT,N1,N2,XTMP,YTMP,ZTMP,PTMP)
   POTF(:,:,TT)=POTF(:,:,TT)+PTMP(:,:)
END IF !(IABS(IH) == 1)
!-----------------------------------------------------------------------------------------------!
DPOTDT(:,:)=0.D0
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PRESFSTD
!-----------------------------------------------------------------------------------------------!
