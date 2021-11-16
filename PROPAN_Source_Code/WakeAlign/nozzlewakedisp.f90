!-----------------------------------------------------------------------------------------------!
!    Nozzle wake displacement                                                                   !
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
SUBROUTINE NOZZLEWAKEDISP(I,XXNEW,RRNEW,TANEW,XPW_TMP,RPW_TMP,TPW_TMP,XNW_TMP,RNW_TMP,TNW_TMP)
!-----------------------------------------------------------------------------------------------!
!    Created by: Joao Baltazar, IST, August 2012                                                !
!    Modified  : 26052014, J. Baltazar, Wake Alignment Module                                   !
!    Modified  : 15102014, J. Baltazar, version 2.1, Wake Alignment Model 2                     !
!    Modified  : 18042019, J. Baltazar, 2019 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,II,IP1,IIP1
DOUBLE PRECISION :: XXNEW(NRW1),RRNEW(NRW1),TANEW(NRW1),XDIS(NRW1),RDIS(NRW1),TDIS(NRW1)
DOUBLE PRECISION :: XPW_TMP(NPW1,NRW1),RPW_TMP(NPW1,NRW1),TPW_TMP(NPW1,NRW1)
DOUBLE PRECISION :: XNW_TMP(NNW1,NNTT),RNW_TMP(NNW1,NNTT),TNW_TMP(NNW1,NNTT)
!-----------------------------------------------------------------------------------------------!
!    Displacement                                                                               !
!-----------------------------------------------------------------------------------------------!
IP1=I+1
XDIS(:)=XXNEW(:)-XPW_TMP(IP1,:)
RDIS(:)=RRNEW(:)-RPW_TMP(IP1,:)
TDIS(:)=TANEW(:)-TPW_TMP(IP1,:)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Geometry Update                                                                     !
!-----------------------------------------------------------------------------------------------!
IF (IP1 > NND) THEN
   IP1=I+1-NND
ELSE
   IP1=1
END IF !(IP1 > NND)
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   XNW(IP1,:)=XNW_TMP(IP1,:)
   RNW(IP1,:)=RNW_TMP(IP1,:)
   TNW(IP1,:)=TNW_TMP(IP1,:)+TDIS(NRW1)
   YNW(IP1,:)=RNW(IP1,:)*DCOS(TNW(IP1,:))
   ZNW(IP1,:)=RNW(IP1,:)*DSIN(TNW(IP1,:))
!-----------------------------------------------------------------------------------------------!
!    Passive Update                                                                             !
!-----------------------------------------------------------------------------------------------!
   IIP1=IP1+1
   DO II=IIP1,NNW1
      XNW(II,:)=XNW_TMP(II,:)
      RNW(II,:)=RNW_TMP(II,:)
      TNW(II,:)=TNW_TMP(II,:)+TDIS(NRW1)
      YNW(II,:)=RNW(II,:)*DCOS(TNW(II,:))
      ZNW(II,:)=RNW(II,:)*DSIN(TNW(II,:))
   END DO !II=IIP1,NNW1
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE NOZZLEWAKEDISP
!-----------------------------------------------------------------------------------------------!
