!-----------------------------------------------------------------------------------------------!
!    The subroutine NUMDIFBLADE computes the coefficients for finite differencing on the        !
!    blade mesh.                                                                                !
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
SUBROUTINE NUMDIFBLADE
!-----------------------------------------------------------------------------------------------!
!    Created by: J.A.C. Falcao de Campos, IST                                                   !
!    Modified  : 05112013, J. Baltazar, version 1.0                                             !
!    Modified  : 26052014, J. Baltazar, Wake Alignment Module                                   !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J
DOUBLE PRECISION :: SAX,SAY,SAZ,SBX,SBY,SBZ,SCX,SCY,SCZ,SA1,SB1,SC1,SA,SB,SS,SP
DOUBLE PRECISION :: T1,T2,ET2X,ET2Y,ET2Z,ET2,ET12,ET22
DOUBLE PRECISION :: T1X,T1Y,T1Z,T2X,T2Y,T2Z,AA1,BB1,CC1,AA2,BB2,CC2
!-----------------------------------------------------------------------------------------------!
!    Propeller Blade                                                                            !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRP
   DO I=1,NCP
!-----------------------------------------------------------------------------------------------!
!    u1 Derivative                                                                              !
!-----------------------------------------------------------------------------------------------!
      IF (I == 1) THEN
         SAX=AT1XP(1,J)
         SAY=AT1YP(1,J)
         SAZ=AT1ZP(1,J)
         SBX=AT1XP(2,J)
         SBY=AT1YP(2,J)
         SBZ=AT1ZP(2,J)
         SCX=AT1XP(3,J)
         SCY=AT1YP(3,J)
         SCZ=AT1ZP(3,J)
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
         A1P(I,J)=AA1
         B1P(I,J)=BB1
         C1P(I,J)=CC1
      END IF !(I == 1)
!-----------------------------------------------------------------------------------------------!
      IF (I == NCP) THEN
         SAX=AT1XP(NCP-2,J)
         SAY=AT1YP(NCP-2,J)
         SAZ=AT1ZP(NCP-2,J)
         SBX=AT1XP(NCP-1,J)
         SBY=AT1YP(NCP-1,J)
         SBZ=AT1ZP(NCP-1,J)
         SCX=AT1XP(NCP  ,J)
         SCY=AT1YP(NCP  ,J)
         SCZ=AT1ZP(NCP  ,J)
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
         A1P(I,J)=AA1
         B1P(I,J)=BB1
         C1P(I,J)=CC1
      END IF !(I == NCP)
!-----------------------------------------------------------------------------------------------!
      IF ((I /= 1).AND.(I /= NCP)) THEN
         SAX=AT1XP(I-1,J)
         SAY=AT1YP(I-1,J)
         SAZ=AT1ZP(I-1,J)
         SBX=AT1XP(I  ,J)
         SBY=AT1YP(I  ,J)
         SBZ=AT1ZP(I  ,J)
         SCX=AT1XP(I+1,J)
         SCY=AT1YP(I+1,J)
         SCZ=AT1ZP(I+1,J)
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
         A1P(I,J)=AA1
         B1P(I,J)=BB1
         C1P(I,J)=CC1
      END IF !((I /= 1).AND.(I /= NCP))
!-----------------------------------------------------------------------------------------------!
      IF (J == NRP.AND.ISTRIP == 2) THEN
!-----------------------------------------------------------------------------------------------!
         IF (I == 1.OR.I == NC+1) THEN
            SAX=AT1XP(I  ,J)
            SAY=AT1YP(I  ,J)
            SAZ=AT1ZP(I  ,J)
            SBX=AT1XP(I+1,J)
            SBY=AT1YP(I+1,J)
            SBZ=AT1ZP(I+1,J)
            SCX=AT1XP(I+2,J)
            SCY=AT1YP(I+2,J)
            SCZ=AT1ZP(I+2,J)
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
            A1P(I,J)=AA1
            B1P(I,J)=BB1
            C1P(I,J)=CC1
!-----------------------------------------------------------------------------------------------!
         ELSEIF (I == NC.OR.I == NCP) THEN
            SAX=AT1XP(I-2,J)
            SAY=AT1YP(I-2,J)
            SAZ=AT1ZP(I-2,J)
            SBX=AT1XP(I-1,J)
            SBY=AT1YP(I-1,J)
            SBZ=AT1ZP(I-1,J)
            SCX=AT1XP(I  ,J)
            SCY=AT1YP(I  ,J)
            SCZ=AT1ZP(I  ,J)
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
            A1P(I,J)=AA1
            B1P(I,J)=BB1
            C1P(I,J)=CC1
!-----------------------------------------------------------------------------------------------!
         ELSE
            SAX=AT1XP(I-1,J)
            SAY=AT1YP(I-1,J)
            SAZ=AT1ZP(I-1,J)
            SBX=AT1XP(I  ,J)
            SBY=AT1YP(I  ,J)
            SBZ=AT1ZP(I  ,J)
            SCX=AT1XP(I+1,J)
            SCY=AT1YP(I+1,J)
            SCZ=AT1ZP(I+1,J)
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
            A1P(I,J)=AA1
            B1P(I,J)=BB1
            C1P(I,J)=CC1
         END IF
!-----------------------------------------------------------------------------------------------!
      END IF !(J == NRP.AND.ISTRIP == 2)
!-----------------------------------------------------------------------------------------------!
!    Compute t1                                                                                 !
!-----------------------------------------------------------------------------------------------!
      T1=DSQRT(T1X*T1X+T1Y*T1Y+T1Z*T1Z)
      T1X=T1X/T1
      T1Y=T1Y/T1
      T1Z=T1Z/T1
      ET1XP(I,J)=T1X
      ET1YP(I,J)=T1Y
      ET1ZP(I,J)=T1Z
!-----------------------------------------------------------------------------------------------!
!    u2 Derivative                                                                              !
!-----------------------------------------------------------------------------------------------!
      IF (J == 1) THEN
         SAX=AT2XP(I,1)
         SAY=AT2YP(I,1)
         SAZ=AT2ZP(I,1)
         SBX=AT2XP(I,2)
         SBY=AT2YP(I,2)
         SBZ=AT2ZP(I,2)
         SCX=AT2XP(I,3)
         SCY=AT2YP(I,3)
         SCZ=AT2ZP(I,3)
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA=SA1+SB1
         SB=SB1+SC1
!-----------------------------------------------------------------------------------------------!
         T2X=SAX
         T2Y=SAY
         T2Z=SAZ
         SP =SA*SB
         SS =SA+SB
         AA2=-(2.D0*SA+SB)/(SA*SS)
         BB2=SS/SP
         CC2=-SA/(SB*SS)
         AA2P(I,J)=AA2
         BB2P(I,J)=BB2
         CC2P(I,J)=CC2
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
         T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
         T2X=T2X/T2
         T2Y=T2Y/T2
         T2Z=T2Z/T2
         ET2X=T1Z*UNYP0(I,J)-T1Y*UNZP0(I,J)
         ET2Y=T1X*UNZP0(I,J)-T1Z*UNXP0(I,J)
         ET2Z=T1Y*UNXP0(I,J)-T1X*UNYP0(I,J)
         ET2=DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
         ET2X=ET2X/ET2
         ET2Y=ET2Y/ET2
         ET2Z=ET2Z/ET2
         ET2XP(I,J)=ET2X
         ET2YP(I,J)=ET2Y
         ET2ZP(I,J)=ET2Z
         ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
         ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
         IF (I == 1) THEN
            A2P(I,J)=(AA2-ET12*AA1)/ET22
            B2P(I,J)=BB2/ET22
            C2P(I,J)=CC2/ET22
            D2P(I,J)=-BB1*ET12/ET22
            E2P(I,J)=-CC1*ET12/ET22
         END IF !(I == 1)
         IF (I == NCP) THEN
            A2P(I,J)=(AA2-ET12*CC1)/ET22
            B2P(I,J)=BB2/ET22
            C2P(I,J)=CC2/ET22
            D2P(I,J)=-BB1*ET12/ET22
            E2P(I,J)=-AA1*ET12/ET22
         END IF !(I == NCP)
         IF ((I /= 1).AND.(I /= NCP)) THEN
            A2P(I,J)=(AA2-ET12*BB1)/ET22
            B2P(I,J)=BB2/ET22
            C2P(I,J)=CC2/ET22
            D2P(I,J)=-AA1*ET12/ET22
            E2P(I,J)=-CC1*ET12/ET22
         END IF !((I /= 1).AND.(I /= NCP))
      END IF !(J == 1)
!-----------------------------------------------------------------------------------------------!
      IF (J == NRP) THEN
         SAX=AT2XP(I,NRP-2)
         SAY=AT2YP(I,NRP-2)
         SAZ=AT2ZP(I,NRP-2)
         SBX=AT2XP(I,NRP-1)
         SBY=AT2YP(I,NRP-1)
         SBZ=AT2ZP(I,NRP-1)
         SCX=AT2XP(I,NRP  )
         SCY=AT2YP(I,NRP  )
         SCZ=AT2ZP(I,NRP  )
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA=SA1+SB1
         SB=SB1+SC1
!-----------------------------------------------------------------------------------------------!
         T2X=SCX
         T2Y=SCY
         T2Z=SCZ
         SP =SA*SB
         SS =SA+SB
         AA2=SB/(SA*SS)
         BB2=-SS/SP
         CC2=(SA+2.D0*SB)/(SB*SS)
         AA2P(I,J)=AA2
         BB2P(I,J)=BB2
         CC2P(I,J)=CC2
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
         T2=DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
         T2X=T2X/T2
         T2Y=T2Y/T2
         T2Z=T2Z/T2
         ET2X=T1Z*UNYP0(I,J)-T1Y*UNZP0(I,J)
         ET2Y=T1X*UNZP0(I,J)-T1Z*UNXP0(I,J)
         ET2Z=T1Y*UNXP0(I,J)-T1X*UNYP0(I,J)
         ET2=DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
         ET2X=ET2X/ET2
         ET2Y=ET2Y/ET2
         ET2Z=ET2Z/ET2
         ET2XP(I,J)=ET2X
         ET2YP(I,J)=ET2Y
         ET2ZP(I,J)=ET2Z
         ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
         ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
         IF (I == 1) THEN
            A2P(I,J)=AA2/ET22
            B2P(I,J)=BB2/ET22
            C2P(I,J)=(CC2-ET12*AA1)/ET22
            D2P(I,J)=-BB1*ET12/ET22
            E2P(I,J)=-CC1*ET12/ET22
         END IF !(I == 1)
         IF (I == NCP) THEN
            A2P(I,J)=AA2/ET22
            B2P(I,J)=BB2/ET22
            C2P(I,J)=(CC2-ET12*CC1)/ET22
            D2P(I,J)=-BB1*ET12/ET22
            E2P(I,J)=-AA1*ET12/ET22
         END IF !(I == NCP)
         IF ((I /= 1).AND.(I /= NCP)) THEN
            A2P(I,J)=AA2/ET22
            B2P(I,J)=BB2/ET22
            C2P(I,J)=(CC2-ET12*BB1)/ET22
            D2P(I,J)=-AA1*ET12/ET22
            E2P(I,J)=-CC1*ET12/ET22
         END IF !((I /= 1).AND.(I /= NCP))
      END IF !(J == NRP)
!-----------------------------------------------------------------------------------------------!
      IF ((J > 1).AND.(J < NRP)) THEN
         SAX=AT2XP(I,J-1)
         SAY=AT2YP(I,J-1)
         SAZ=AT2ZP(I,J-1)
         SBX=AT2XP(I,J  )
         SBY=AT2YP(I,J  )
         SBZ=AT2ZP(I,J  )
         SCX=AT2XP(I,J+1)
         SCY=AT2YP(I,J+1)
         SCZ=AT2ZP(I,J+1)
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA =SA1+SB1
         SB =SB1+SC1
!-----------------------------------------------------------------------------------------------!
         SP =SA*SB
         SS =SA+SB
         T2X=SBX
         T2Y=SBY
         T2Z=SBZ
         AA2=-SB/(SA*SS)
         BB2=(SB-SA)/SP
         CC2= SA/(SB*SS)
         AA2P(I,J)=AA2
         BB2P(I,J)=BB2
         CC2P(I,J)=CC2
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
         T2=DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
         T2X=T2X/T2
         T2Y=T2Y/T2
         T2Z=T2Z/T2
         ET2X=T1Z*UNYP0(I,J)-T1Y*UNZP0(I,J)
         ET2Y=T1X*UNZP0(I,J)-T1Z*UNXP0(I,J)
         ET2Z=T1Y*UNXP0(I,J)-T1X*UNYP0(I,J)
         ET2=DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
         ET2X=ET2X/ET2
         ET2Y=ET2Y/ET2
         ET2Z=ET2Z/ET2
         ET2XP(I,J)=ET2X
         ET2YP(I,J)=ET2Y
         ET2ZP(I,J)=ET2Z
         ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
         ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
         IF (I == 1) THEN
            A2P(I,J)=AA2/ET22
            B2P(I,J)=(BB2-ET12*AA1)/ET22
            C2P(I,J)=CC2/ET22
            D2P(I,J)=-BB1*ET12/ET22
            E2P(I,J)=-CC1*ET12/ET22
         END IF !(I == 1)
         IF (I == NCP) THEN
            A2P(I,J)=AA2/ET22
            B2P(I,J)=(BB2-ET12*CC1)/ET22
            C2P(I,J)=CC2/ET22
            D2P(I,J)=-BB1*ET12/ET22
            E2P(I,J)=-AA1*ET12/ET22
         END IF !(I == NCP)
         IF ((I /= 1).AND.(I /= NCP)) THEN
            A2P(I,J)=AA2/ET22
            B2P(I,J)=(BB2-ET12*BB1)/ET22
            C2P(I,J)=CC2/ET22
            D2P(I,J)=-AA1*ET12/ET22
            E2P(I,J)=-CC1*ET12/ET22
         END IF !((I /= 1).AND.(I /= NCP))
      END IF !((J > 1).AND.(J < NRP))
!-----------------------------------------------------------------------------------------------!
      IF (J == NRP.AND.ISTRIP == 2) THEN
         SBX=AT2XP(I      ,J)
         SBY=AT2YP(I      ,J)
         SBZ=AT2ZP(I      ,J)
         SCX=AT2XP(NCP-I+1,J)
         SCY=AT2YP(NCP-I+1,J)
         SCZ=AT2ZP(NCP-I+1,J)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA =SB1
         SB =SC1
!-----------------------------------------------------------------------------------------------!
         SS =SA+SB
         T2X=SBX
         T2Y=SBY
         T2Z=SBZ
         BB2=-1.D0/SS
         CC2= 1.D0/SS
         AA2P(I,J)=0.D0
         BB2P(I,J)=BB2
         CC2P(I,J)=CC2
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
         T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
         T2X=T2X/T2
         T2Y=T2Y/T2
         T2Z=T2Z/T2
         ET2X=T1Z*UNYP0(I,J)-T1Y*UNZP0(I,J)
         ET2Y=T1X*UNZP0(I,J)-T1Z*UNXP0(I,J)
         ET2Z=T1Y*UNXP0(I,J)-T1X*UNYP0(I,J)
         ET2=DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
         ET2X=ET2X/ET2
         ET2Y=ET2Y/ET2
         ET2Z=ET2Z/ET2
         ET2XP(I,J)=ET2X
         ET2YP(I,J)=ET2Y
         ET2ZP(I,J)=ET2Z
         ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
         ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
         IF (I == 1.OR.I == NC+1) THEN
            A2P(I,J)=0.D0
            B2P(I,J)=(BB2-ET12*AA1)/ET22
            C2P(I,J)=CC2/ET22
            D2P(I,J)=-BB1*ET12/ET22
            E2P(I,J)=-CC1*ET12/ET22
         ELSEIF (I == NC.OR.I == NCP) THEN
            A2P(I,J)=0.D0
            B2P(I,J)=(BB2-ET12*CC1)/ET22
            C2P(I,J)=CC2/ET22
            D2P(I,J)=-BB1*ET12/ET22
            E2P(I,J)=-AA1*ET12/ET22
         ELSE
            A2P(I,J)=0.D0
            B2P(I,J)=(BB2-ET12*BB1)/ET22
            C2P(I,J)=CC2/ET22
            D2P(I,J)=-AA1*ET12/ET22
            E2P(I,J)=-CC1*ET12/ET22
         END IF
      END IF !(J == NRP.AND.ISTRIP == 2)
!-----------------------------------------------------------------------------------------------!
   END DO !I=1,NCP
END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE NUMDIFBLADE
!-----------------------------------------------------------------------------------------------!
