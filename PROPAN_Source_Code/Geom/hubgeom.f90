!-----------------------------------------------------------------------------------------------!
!    Compute Centroid Unit Normal and Area                                                      !
!    Copyright (C) 2021  J. Baltazar and J.A.C. Falcão de Campos                                !
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
SUBROUTINE HUBGEOM
!-----------------------------------------------------------------------------------------------!
!    Created by: 04112013, J. Baltazar, IST                                                     !
!    Modified  : 04112013, J. Baltazar, version 1.0                                             !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4)
DOUBLE PRECISION :: X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0
!-----------------------------------------------------------------------------------------------!
!    Loop on the Hub Panels                                                                     !
!-----------------------------------------------------------------------------------------------!
!    First Half Sector                                                                          !
!-----------------------------------------------------------------------------------------------!
DO J=1,NHT
   DO I=1,NHX
      XX(1)=XH(I+1,J  )
      XX(2)=XH(I  ,J  )
      XX(3)=XH(I  ,J+1)
      XX(4)=XH(I+1,J+1)
      YY(1)=YH(I+1,J  )
      YY(2)=YH(I  ,J  )
      YY(3)=YH(I  ,J+1)
      YY(4)=YH(I+1,J+1)
      ZZ(1)=ZH(I+1,J  )
      ZZ(2)=ZH(I  ,J  )
      ZZ(3)=ZH(I  ,J+1)
      ZZ(4)=ZH(I+1,J+1)
!-----------------------------------------------------------------------------------------------!
!    For Flat Panel Redefine Corner Points                                                      !
!-----------------------------------------------------------------------------------------------!
      IF (IPAN == 0) THEN
         CALL PANEL(XX,YY,ZZ,X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0)
         CALL PANELFLAT(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0)
      END IF !(IPAN == 0)
!-----------------------------------------------------------------------------------------------!
!    Compute Panel Centroid Data                                                                !
!-----------------------------------------------------------------------------------------------!
      CALL PANEL(XX,YY,ZZ,XH0(I,J),YH0(I,J),ZH0(I,J),AT1XH(I,J),AT1YH(I,J),AT1ZH(I,J), &
                 AT2XH(I,J),AT2YH(I,J),AT2ZH(I,J),UNXH0(I,J),UNYH0(I,J),UNZH0(I,J),AH0(I,J))
   END DO !I=1,NHX
END DO !J=1,NHT
!-----------------------------------------------------------------------------------------------!
!    Second Half Sector                                                                         !
!-----------------------------------------------------------------------------------------------!
DO J=NHT2,NHTT1
   DO I=1,NHX
      XX(1)=XH(I+1,J  )
      XX(2)=XH(I  ,J  )
      XX(3)=XH(I  ,J+1)
      XX(4)=XH(I+1,J+1)
      YY(1)=YH(I+1,J  )
      YY(2)=YH(I  ,J  )
      YY(3)=YH(I  ,J+1)
      YY(4)=YH(I+1,J+1)
      ZZ(1)=ZH(I+1,J  )
      ZZ(2)=ZH(I  ,J  )
      ZZ(3)=ZH(I  ,J+1)
      ZZ(4)=ZH(I+1,J+1)
!-----------------------------------------------------------------------------------------------!
!    For Flat Panel Redefine Corner Points                                                      !
!-----------------------------------------------------------------------------------------------!
      IF (IPAN == 0) THEN
         CALL PANEL(XX,YY,ZZ,X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0)
         CALL PANELFLAT(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0)
      END IF !(IPAN == 0)
!-----------------------------------------------------------------------------------------------!
!    Compute Panel Centroid Data                                                                !
!-----------------------------------------------------------------------------------------------!
      CALL PANEL(XX,YY,ZZ,XH0(I,J-1),YH0(I,J-1),ZH0(I,J-1), &
                 AT1XH(I,J-1),AT1YH(I,J-1),AT1ZH(I,J-1), &
                 AT2XH(I,J-1),AT2YH(I,J-1),AT2ZH(I,J-1), &
                 UNXH0(I,J-1),UNYH0(I,J-1),UNZH0(I,J-1),AH0(I,J-1))
   END DO !I=1,NHX
END DO !J=NHT2,NHTT1
!-----------------------------------------------------------------------------------------------!
!    Fill in Extra J Lines for Periodic Boundary Condition on the Hub                           !
!-----------------------------------------------------------------------------------------------!
DO I=1,NHX
   WORK1=DSQRT(YH0(I,NHTP)*YH0(I,NHTP)+ZH0(I,NHTP)*ZH0(I,NHTP))
   WORK2=DATAN2(ZH0(I,NHTP),YH0(I,NHTP))
   WORK2=WORK2-2.D0*PI/DFLOAT(NB)
   XH0(I,0)=XH0(I,NHTP)
   YH0(I,0)=WORK1*DCOS(WORK2)
   ZH0(I,0)=WORK1*DSIN(WORK2)
!-----------------------------------------------------------------------------------------------!
   AT1XH(I,0)=AT1XH(I,NHTP)
   AT1YH(I,0)=AT1YH(I,NHTP)
   AT1ZH(I,0)=AT1ZH(I,NHTP)
   AT2XH(I,0)=AT2XH(I,NHTP)
   AT2YH(I,0)=AT2YH(I,NHTP)
   AT2ZH(I,0)=AT2ZH(I,NHTP)
!-----------------------------------------------------------------------------------------------!
   WORK1=DSQRT(YH0(I,1)*YH0(I,1)+ZH0(I,1)*ZH0(I,1))
   WORK2=DATAN2(ZH0(I,1),YH0(I,1))
   WORK2=WORK2+2.D0*PI/DFLOAT(NB)
   XH0(I,NHTP1)=XH0(I,1)
   YH0(I,NHTP1)=WORK1*DCOS(WORK2)
   ZH0(I,NHTP1)=WORK1*DSIN(WORK2)
!-----------------------------------------------------------------------------------------------!
   AT1XH(I,NHTP1)=AT1XH(I,1)
   AT1YH(I,NHTP1)=AT1YH(I,1)
   AT1ZH(I,NHTP1)=AT1ZH(I,1)
   AT2XH(I,NHTP1)=AT2XH(I,1)
   AT2YH(I,NHTP1)=AT2YH(I,1)
   AT2ZH(I,NHTP1)=AT2ZH(I,1)
END DO !I=1,NHX
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE HUBGEOM
!-----------------------------------------------------------------------------------------------!
