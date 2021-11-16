!-----------------------------------------------------------------------------------------------!
!    Nozzle displacement                                                                        !
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
SUBROUTINE NOZZLEDISP(I,XXNEW,RRNEW,TANEW,XPW_TMP,RPW_TMP,TPW_TMP,XN_TMP,RN_TMP,TN_TMP)
!-----------------------------------------------------------------------------------------------!
!    Created by: Joao Baltazar, IST, August 2012                                                !
!    Modified  : 26052014, J. Baltazar, Wake Alignment Module                                   !
!    Modified  : 15102014, J. Baltazar, version 2.1, Wake Alignment Model 2                     !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,II,IP1,IIP,IIP1
DOUBLE PRECISION :: XXNEW(NRW1),RRNEW(NRW1),TANEW(NRW1),XDIS(NRW1),RDIS(NRW1),TDIS(NRW1)
DOUBLE PRECISION :: XPW_TMP(NPW1,NRW1),RPW_TMP(NPW1,NRW1),TPW_TMP(NPW1,NRW1)
DOUBLE PRECISION ::  XN_TMP(NNXT,NNTT), RN_TMP(NNXT,NNTT), TN_TMP(NNXT,NNTT)
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
IP1=NNU+NC+I+1
IF (IP1 <= NNX1) THEN
   XN(IP1,:)=XN_TMP(IP1,:)
   RN(IP1,:)=RN_TMP(IP1,:)
   TN(IP1,:)=TN_TMP(IP1,:)+TDIS(NRW1)
   YN(IP1,:)=RN(IP1,:)*DCOS(TN(IP1,:))
   ZN(IP1,:)=RN(IP1,:)*DSIN(TN(IP1,:))
!-----------------------------------------------------------------------------------------------!
!    Passive Update                                                                             !
!-----------------------------------------------------------------------------------------------!
   IIP=1
   IIP1=IP1+1
   DO II=IIP1,NNX1
      XN(II,:)=XN_TMP(II,:)
      RN(II,:)=RN_TMP(II,:)
      TN(II,:)=TN_TMP(II,:)+TDIS(NRW1)
      YN(II,:)=RN(II,:)*DCOS(TN(II,:))
      ZN(II,:)=RN(II,:)*DSIN(TN(II,:))
      IIP=IIP+1
   END DO !II=IIP1,NNX1
!-----------------------------------------------------------------------------------------------!
   IIP=IIP+NNX1
   DO II=NNX1+1,IIP
      XN(II,:)=XN_TMP(II,:)
      RN(II,:)=RN_TMP(II,:)
      TN(II,:)=TN_TMP(II,:)+TDIS(NRW1)
      YN(II,:)=RN(II,:)*DCOS(TN(II,:))
      ZN(II,:)=RN(II,:)*DSIN(TN(II,:))
   END DO !II=NNX1+1,IIP
END IF !(IP1 <= NNX1)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE NOZZLEDISP
!-----------------------------------------------------------------------------------------------!
