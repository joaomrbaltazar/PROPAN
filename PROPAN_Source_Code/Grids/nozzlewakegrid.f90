!-----------------------------------------------------------------------------------------------!
!    Generate nozzle wake grid                                                                  !
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
SUBROUTINE NOZZLEWAKEGRID
!-----------------------------------------------------------------------------------------------!
!    Created by: 21102013, J. Baltazar, version 1.0                                             !
!    Modified  : 08102014, J. Baltazar, version v2.1                                            !
!-----------------------------------------------------------------------------------------------!
!    Declarations                                                                               !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake Grid                                                                           !
!-----------------------------------------------------------------------------------------------!
DO J=1,NNTT
   XNW(1,J)=XN(NNX1,J)
   RNW(1,J)=RN(NNX1,J)
   TNW(1,J)=TN(NNX1,J)
   YNW(1,J)=RNW(1,J)*DCOS(TNW(1,J))
   ZNW(1,J)=RNW(1,J)*DSIN(TNW(1,J))
!-----------------------------------------------------------------------------------------------!
!    Loop on Axial Stations                                                                     !
!-----------------------------------------------------------------------------------------------!
   DO I=2,NNW1
      XNW(I,J)=XPW(NND+I,NRW1)
      TNW(I,J)=TNW(I-1,J)+(TPW(NND+I,NRW1)-TPW(NND+I-1,NRW1))
      RNW(I,J)=RNW(1,J)
      YNW(I,J)=RNW(I,J)*DCOS(TNW(I,J))
      ZNW(I,J)=RNW(I,J)*DSIN(TNW(I,J))
   END DO !I=2,NNW1
END DO !J=1,NNTT
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE NOZZLEWAKEGRID
!-----------------------------------------------------------------------------------------------!
