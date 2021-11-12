!-----------------------------------------------------------------------------------------------!
!    Ship Wake Field                                                                            !
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
SUBROUTINE VWAKE(TT,X,Y,Z,KB,NFT,VWX,VWY,VWZ)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 28102013, J. Baltazar, version 1.0                                             !
!                16112013: J. Baltazar, case NTETA=0                                            !
!    Modified  : 05022016, J. Baltazar, 2016 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: J,KB,TT,NFT
DOUBLE PRECISION :: X,Y,Z,XNEW,YNEW,ZNEW,RNEW,TNEW
DOUBLE PRECISION :: AN(0:NFT),BN(1:NFT),VXX(NWR),VTT(NWR),VRR(NWR),VWX,VWY,VWZ,VWT,VWR
!-----------------------------------------------------------------------------------------------!
!    Rotating Reference Frame to Inertial Reference Frame                                       !
!-----------------------------------------------------------------------------------------------!
CALL FRAME('INERTIAL',TT,KB,X,Y,Z,XNEW,YNEW,ZNEW,RNEW,TNEW)
!-----------------------------------------------------------------------------------------------!
!    Axial Velocity                                                                             !
!-----------------------------------------------------------------------------------------------!
DO J=1,NWR
   CALL FOURIER_FUNCTION(NFT,ANVX(:,J),BNVX(:,J),TNEW,VXX(J))
END DO !J=1,NWR
CALL LININT(NWR,WRR,VXX,1,RNEW,VWX)
!-----------------------------------------------------------------------------------------------!
!    Tangential Velocity                                                                        !
!-----------------------------------------------------------------------------------------------!
DO J=1,NWR
   CALL FOURIER_FUNCTION(NFT,ANVT(:,J),BNVT(:,J),TNEW,VTT(J))
END DO !J=1,NWR
CALL LININT(NWR,WRR,VTT,1,RNEW,VWT)
!-----------------------------------------------------------------------------------------------!
!    Radial Velocity                                                                            !
!-----------------------------------------------------------------------------------------------!
DO J=1,NWR
   CALL FOURIER_FUNCTION(NFT,ANVR(:,J),BNVR(:,J),TNEW,VRR(J))
END DO !J=1,NWR
CALL LININT(NWR,WRR,VRR,1,RNEW,VWR)
!-----------------------------------------------------------------------------------------------!
!    Cartesian Components in the Inertial Reference Frame                                       !
!-----------------------------------------------------------------------------------------------!
VWY=VWR*DCOS(DATAN2(Z,Y))-VWT*DSIN(DATAN2(Z,Y))
VWZ=VWR*DSIN(DATAN2(Z,Y))+VWT*DCOS(DATAN2(Z,Y))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE VWAKE
!-----------------------------------------------------------------------------------------------!
