!-----------------------------------------------------------------------------------------------!
!    Delete variables                                                                           !
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
SUBROUTINE DELVARS
!-----------------------------------------------------------------------------------------------!
!    Created by: 14102013, J. Baltazar, version 1.0                                             !
!    Modified  : 28102013, J. Baltazar, version 1.0                                             !
!-----------------------------------------------------------------------------------------------!
USE PROPANEL_MOD
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------!
!    Blade Variables                                                                            !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) DEALLOCATE(XP,YP,ZP,RP,TP)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake Variables                                                                       !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) DEALLOCATE(XPW,YPW,ZPW,RPW,TPW)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Variables                                                                           !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) DEALLOCATE(XN,YN,ZN,RN,TN)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake Variables                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) DEALLOCATE(XNW,YNW,ZNW,RNW,TNW)
!-----------------------------------------------------------------------------------------------!
!    Hub Variables                                                                              !
!-----------------------------------------------------------------------------------------------!
IF (ABS(IH) == 1) DEALLOCATE(XH,YH,ZH,RH,TH)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE DELVARS
!-----------------------------------------------------------------------------------------------!
