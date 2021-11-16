!-----------------------------------------------------------------------------------------------!
!    Blade wake displacement                                                                    !
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
SUBROUTINE BLADEWAKEDISP(I,XXNEW,RRNEW,TANEW,XPW_TMP,RPW_TMP,TPW_TMP)
!-----------------------------------------------------------------------------------------------!
!    Created by: Joao Baltazar, IST, July 2009                                                  !
!    Modified  : 26052014, J. Baltazar, Wake Alignment Module                                   !
!    Modified  : 15102014, J. Baltazar, version 2.1, Wake Alignment Model 2                     !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,II,IP1,IIP1
DOUBLE PRECISION :: XXNEW(NRW1),RRNEW(NRW1),TANEW(NRW1),XDIS(NRW1),RDIS(NRW1),TDIS(NRW1)
DOUBLE PRECISION :: XPW_TMP(NPW1,NRW1),RPW_TMP(NPW1,NRW1),TPW_TMP(NPW1,NRW1)
!-----------------------------------------------------------------------------------------------!
!    Displacement                                                                               !
!-----------------------------------------------------------------------------------------------!
IP1=I+1
XDIS(:)=XXNEW(:)-XPW_TMP(IP1,:)
RDIS(:)=RRNEW(:)-RPW_TMP(IP1,:)
TDIS(:)=TANEW(:)-TPW_TMP(IP1,:)
!-----------------------------------------------------------------------------------------------!
!    Wake Geometry Update                                                                       !
!-----------------------------------------------------------------------------------------------!
IP1=I+1
XPW(IP1,:)=XXNEW(:)
RPW(IP1,:)=RRNEW(:)
TPW(IP1,:)=TANEW(:)
YPW(IP1,:)=RRNEW(:)*DCOS(TANEW(:))
ZPW(IP1,:)=RRNEW(:)*DSIN(TANEW(:))
!-----------------------------------------------------------------------------------------------!
!    Passive Update                                                                             !
!-----------------------------------------------------------------------------------------------!
DO II=I+1,NPW
   IIP1=II+1
   XPW(IIP1,:)=XPW_TMP(IIP1,:) !+XDIS(:)
   RPW(IIP1,:)=RPW_TMP(IIP1,:) !+RDIS(:)
   TPW(IIP1,:)=TPW_TMP(IIP1,:)+TDIS(:)
   YPW(IIP1,:)=RPW(IIP1,:)*DCOS(TPW(IIP1,:))
   ZPW(IIP1,:)=RPW(IIP1,:)*DSIN(TPW(IIP1,:))
END DO !II=I+1,NPW
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE BLADEWAKEDISP
!-----------------------------------------------------------------------------------------------!
