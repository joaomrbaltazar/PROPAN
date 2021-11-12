!-----------------------------------------------------------------------------------------------!
!    Subroutine NUMDIFBLADEWAKE                                                                 !
!                                                                                               !
!    The subroutine NUMDIFBLADEWAKE computes the coefficients for finite differencing on the    !
!    blade wake mesh.                                                                           !
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
SUBROUTINE NUMDIFBLADEWAKE
!-----------------------------------------------------------------------------------------------!
!    Created by: 27112014, J. Baltazar, version 3.3, Super-Cavitation Model                     !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J
DOUBLE PRECISION :: SAX,SAY,SAZ,SBX,SBY,SBZ,SCX,SCY,SCZ,SA1,SB1,SC1,SA,SB,SS,SP
DOUBLE PRECISION :: T1,T1X,T1Y,T1Z,AA1,BB1,CC1
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRW
   DO I=1,IABS(NCPW)
!-----------------------------------------------------------------------------------------------!
!    u1 Derivative                                                                              !
!-----------------------------------------------------------------------------------------------!
      IF (I == 1) THEN
         SAX=AT1XPW(1,J)
         SAY=AT1YPW(1,J)
         SAZ=AT1ZPW(1,J)
         SBX=AT1XPW(2,J)
         SBY=AT1YPW(2,J)
         SBZ=AT1ZPW(2,J)
         SCX=AT1XPW(3,J)
         SCY=AT1YPW(3,J)
         SCZ=AT1ZPW(3,J)
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA =SA1+SB1
         SB =SB1+SC1
         T1X=SAX
         T1Y=SAY
         T1Z=SAZ
!-----------------------------------------------------------------------------------------------!
         SP =SA*SB
         SS =SA+SB
         AA1=-(2.D0*SA+SB)/(SA*SS)
         BB1=SS/SP
         CC1=-SA/(SB*SS)
         A1PW(I,J)=AA1
         B1PW(I,J)=BB1
         C1PW(I,J)=CC1
      END IF !(I == 1)
!-----------------------------------------------------------------------------------------------!
      IF (I == IABS(NCPW)) THEN
         SAX=AT1XPW(IABS(NCPW)-2,J)
         SAY=AT1YPW(IABS(NCPW)-2,J)
         SAZ=AT1ZPW(IABS(NCPW)-2,J)
         SBX=AT1XPW(IABS(NCPW)-1,J)
         SBY=AT1YPW(IABS(NCPW)-1,J)
         SBZ=AT1ZPW(IABS(NCPW)-1,J)
         SCX=AT1XPW(IABS(NCPW)  ,J)
         SCY=AT1YPW(IABS(NCPW)  ,J)
         SCZ=AT1ZPW(IABS(NCPW)  ,J)
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA =SA1+SB1
         SB =SB1+SC1
         T1X=SCX
         T1Y=SCY
         T1Z=SCZ
!-----------------------------------------------------------------------------------------------!
         SP =SA*SB
         SS =SA+SB
         AA1=SB/(SA*SS)
         BB1=-SS/SP
         CC1=(SA+2.D0*SB)/(SB*SS)
         A1PW(I,J)=AA1
         B1PW(I,J)=BB1
         C1PW(I,J)=CC1
      END IF !(I == IABS(NCPW))
!-----------------------------------------------------------------------------------------------!
      IF ((I /= 1).AND.(I /= IABS(NCPW))) THEN
         SAX=AT1XPW(I-1,J)
         SAY=AT1YPW(I-1,J)
         SAZ=AT1ZPW(I-1,J)
         SBX=AT1XPW(I  ,J)
         SBY=AT1YPW(I  ,J)
         SBZ=AT1ZPW(I  ,J)
         SCX=AT1XPW(I+1,J)
         SCY=AT1YPW(I+1,J)
         SCZ=AT1ZPW(I+1,J)
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA =SA1+SB1
         SB =SB1+SC1
         T1X=SBX
         T1Y=SBY
         T1Z=SBZ
!-----------------------------------------------------------------------------------------------!
         SP =SA*SB
         SS =SA+SB
         AA1=-SB/(SA*SS)
         BB1=(SB-SA)/SP
         CC1= SA/(SB*SS)
         A1PW(I,J)=AA1
         B1PW(I,J)=BB1
         C1PW(I,J)=CC1
      END IF !((I /= 1).AND.(I /= IABS(NCPW)))
!-----------------------------------------------------------------------------------------------!
!    Compute t1                                                                                 !
!-----------------------------------------------------------------------------------------------!
      T1=DSQRT(T1X*T1X+T1Y*T1Y+T1Z*T1Z)
      T1X=T1X/T1
      T1Y=T1Y/T1
      T1Z=T1Z/T1
      ET1XPW(I,J)=T1X
      ET1YPW(I,J)=T1Y
      ET1ZPW(I,J)=T1Z
!-----------------------------------------------------------------------------------------------!
   END DO !I=1,IABS(NCPW)
END DO !J=1,NRW
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE NUMDIFBLADEWAKE
!-----------------------------------------------------------------------------------------------!
