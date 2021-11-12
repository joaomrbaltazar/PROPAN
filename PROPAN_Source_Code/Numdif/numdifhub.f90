!-----------------------------------------------------------------------------------------------!
!    The subroutine NUMDIFH computes the coefficients for finite differencing on the hub mesh   !
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
SUBROUTINE NUMDIFHUB
!-----------------------------------------------------------------------------------------------!
!    Created by: J.A.C. Falcao de Campos, IST                                                   !
!    Modified  : 05112013, J. Baltazar, version 1.0                                             !
!    Modified  : 26052014, J. Baltazar, Wake Alignment Module                                   !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J
DOUBLE PRECISION :: ET2X,ET2Y,ET2Z,ET2,ET12,ET22
DOUBLE PRECISION :: SAX,SAY,SAZ,SBX,SBY,SBZ,SCX,SCY,SCZ,SA1,SB1,SC1,SA,SB,SS,SP
DOUBLE PRECISION :: T1X,T1Y,T1Z,T2X,T2Y,T2Z,AA1,BB1,CC1,AA2,BB2,CC2,T1,T2
!-----------------------------------------------------------------------------------------------!
DO J=1,NHTP
   DO I=1,NHX
!-----------------------------------------------------------------------------------------------!
!    u1 Derivative                                                                              !
!-----------------------------------------------------------------------------------------------!
      IF (I == 1) THEN
         SAX=AT1XH(1,J)
         SAY=AT1YH(1,J)
         SAZ=AT1ZH(1,J)
         SBX=AT1XH(2,J)
         SBY=AT1YH(2,J)
         SBZ=AT1ZH(2,J)
         SCX=AT1XH(3,J)
         SCY=AT1YH(3,J)
         SCZ=AT1ZH(3,J)
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA =SA1+SB1
         SB =SB1+SC1
         T1X=-SAX
         T1Y=-SAY
         T1Z=-SAZ
         SP =SA*SB
         SS =SA+SB
         AA1=-(2.D0*SA+SB)/(SA*SS)
         BB1=SS/SP
         CC1=-SA/(SB*SS)
         A1H(I,J)=AA1
         B1H(I,J)=BB1
         C1H(I,J)=CC1
      END IF !(I == 1)
!-----------------------------------------------------------------------------------------------!
      IF (I == NHX) THEN
         SAX=AT1XH(NHX-2,J)
         SAY=AT1YH(NHX-2,J)
         SAZ=AT1ZH(NHX-2,J)
         SBX=AT1XH(NHX-1,J)
         SBY=AT1YH(NHX-1,J)
         SBZ=AT1ZH(NHX-1,J)
         SCX=AT1XH(NHX  ,J)
         SCY=AT1YH(NHX  ,J)
         SCZ=AT1ZH(NHX  ,J)
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA =SA1+SB1
         SB =SB1+SC1
         T1X=-SCX
         T1Y=-SCY
         T1Z=-SCZ
         SP =SA*SB
         SS =SA+SB
         AA1=SB/(SA*SS)
         BB1=-SS/SP
         CC1=(SA+2.D0*SB)/(SB*SS)
         A1H(I,J)=AA1
         B1H(I,J)=BB1
         C1H(I,J)=CC1
      END IF !(I == NHX)
!-----------------------------------------------------------------------------------------------!
      IF ((I /= 1).AND.(I /= NHX)) THEN
         SAX=AT1XH(I-1,J)
         SAY=AT1YH(I-1,J)
         SAZ=AT1ZH(I-1,J)
         SBX=AT1XH(I  ,J)
         SBY=AT1YH(I  ,J)
         SBZ=AT1ZH(I  ,J)
         SCX=AT1XH(I+1,J)
         SCY=AT1YH(I+1,J)
         SCZ=AT1ZH(I+1,J)
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA =SA1+SB1
         SB =SB1+SC1
         T1X=-SBX
         T1Y=-SBY
         T1Z=-SBZ
         SP =SA*SB
         SS =SA+SB
         AA1=-SB/(SA*SS)
         BB1=(SB-SA)/SP
         CC1=SA/(SB*SS)
         A1H(I,J)=AA1
         B1H(I,J)=BB1
         C1H(I,J)=CC1
      END IF !((I /= 1).AND.(I /= NHX))
!-----------------------------------------------------------------------------------------------!
!    Compute t1                                                                                 !
!-----------------------------------------------------------------------------------------------!
      T1=DSQRT(T1X*T1X+T1Y*T1Y+T1Z*T1Z)
      T1X=T1X/T1
      T1Y=T1Y/T1
      T1Z=T1Z/T1
      ET1XH(I,J)=T1X
      ET1YH(I,J)=T1Y
      ET1ZH(I,J)=T1Z
!-----------------------------------------------------------------------------------------------!
!    u2 Derivative                                                                              !
!-----------------------------------------------------------------------------------------------!
      IF (J == 1) THEN
         SAX=AT2XH(I,NHTP)
         SAY=AT2YH(I,NHTP)
         SAZ=AT2ZH(I,NHTP)
         SBX=AT2XH(I,1   )
         SBY=AT2YH(I,1   )
         SBZ=AT2ZH(I,1   )
         SCX=AT2XH(I,2   )
         SCY=AT2YH(I,2   )
         SCZ=AT2ZH(I,2   )
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA =SA1+SB1
         SB =SB1+SC1
         T2X=SBX
         T2Y=SBY
         T2Z=SBZ
         SP =SA*SB
         SS =SA+SB
         AA2=-SB/(SA*SS)
         BB2=(SB-SA)/SP
         CC2=SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
         T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
         T2X=T2X/T2
         T2Y=T2Y/T2
         T2Z=T2Z/T2
         ET2X=T1Z*UNYH0(I,J)-T1Y*UNZH0(I,J)
         ET2Y=T1X*UNZH0(I,J)-T1Z*UNXH0(I,J)
         ET2Z=T1Y*UNXH0(I,J)-T1X*UNYH0(I,J)
         ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
         ET2X=ET2X/ET2
         ET2Y=ET2Y/ET2
         ET2Z=ET2Z/ET2
         ET2XH(I,J)=ET2X
         ET2YH(I,J)=ET2Y
         ET2ZH(I,J)=ET2Z
         ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
         ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
         IF (I == 1) THEN
            A2H(I,J)=AA2/ET22
            B2H(I,J)=(BB2-ET12*AA1)/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-BB1*ET12/ET22
            E2H(I,J)=-CC1*ET12/ET22
         END IF !(I == 1)
         IF (I == NHX) THEN
            A2H(I,J)=AA2/ET22
            B2H(I,J)=(BB2-ET12*CC1)/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-BB1*ET12/ET22
            E2H(I,J)=-AA1*ET12/ET22
         END IF !(I == NHX)
         IF ((I /= 1).AND.(I /= NHX)) THEN
            A2H(I,J)=AA2/ET22
            B2H(I,J)=(BB2-ET12*BB1)/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-AA1*ET12/ET22
            E2H(I,J)=-CC1*ET12/ET22
         END IF !((I /= 1).AND.(I /= NHX))
      END IF !(J == 1)
!-----------------------------------------------------------------------------------------------!
      IF (J == NHTP) THEN
         SAX=AT2XH(I,NHTP-1)
         SAY=AT2YH(I,NHTP-1)
         SAZ=AT2ZH(I,NHTP-1)
         SBX=AT2XH(I,NHTP  )
         SBY=AT2YH(I,NHTP  )
         SBZ=AT2ZH(I,NHTP  )
         SCX=AT2XH(I,1     )
         SCY=AT2YH(I,1     )
         SCZ=AT2ZH(I,1     )
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA =SA1+SB1
         SB =SB1+SC1
         T2X=SBX
         T2Y=SBY
         T2Z=SBZ
         SP =SA*SB
         SS =SA+SB
         AA2=-SB/(SA*SS)
         BB2=(SB-SA)/SP
         CC2=SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
         T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
         T2X=T2X/T2
         T2Y=T2Y/T2
         T2Z=T2Z/T2
         ET2X=T1Z*UNYH0(I,J)-T1Y*UNZH0(I,J)
         ET2Y=T1X*UNZH0(I,J)-T1Z*UNXH0(I,J)
         ET2Z=T1Y*UNXH0(I,J)-T1X*UNYH0(I,J)
         ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
         ET2X=ET2X/ET2
         ET2Y=ET2Y/ET2
         ET2Z=ET2Z/ET2
         ET2XH(I,J)=ET2X
         ET2YH(I,J)=ET2Y
         ET2ZH(I,J)=ET2Z
         ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
         ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
         IF (I == 1) THEN
            A2H(I,J)=AA2/ET22
            B2H(I,J)=(BB2-ET12*AA1)/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-BB1*ET12/ET22
            E2H(I,J)=-CC1*ET12/ET22
         END IF !(I == 1)
         IF (I == NHX) THEN
            A2H(I,J)=AA2/ET22
            B2H(I,J)=(BB2-ET12*CC1)/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-BB1*ET12/ET22
            E2H(I,J)=-AA1*ET12/ET22
         END IF !(I == NHX)
         IF ((I /= 1).AND.(I /= NHX)) THEN
            A2H(I,J)=AA2/ET22
            B2H(I,J)=(BB2-ET12*BB1)/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-AA1*ET12/ET22
            E2H(I,J)=-CC1*ET12/ET22
         END IF !((I /= 1).AND.(I /= NHX))
      END IF !(J == NHTP)
!-----------------------------------------------------------------------------------------------!
      IF (J == NHT) THEN
         IF (I == 1) THEN
            SAX=AT2XH(1,NHT-1)
            SAY=AT2YH(1,NHT-1)
            SAZ=AT2ZH(1,NHT-1)
            SBX=AT2XH(1,NHT  )
            SBY=AT2YH(1,NHT  )
            SBZ=AT2ZH(1,NHT  )
            SCX=AT2XH(1,NHT+1)
            SCY=AT2YH(1,NHT+1)
            SCZ=AT2ZH(1,NHT+1)
            SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
            SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
            SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
            SA =SA1+SB1
            SB =SB1+SC1
            T2X=SBX
            T2Y=SBY
            T2Z=SBZ
            SP =SA*SB
            SS =SA+SB
            AA2=-SB/(SA*SS)
            BB2=(SB-SA)/SP
            CC2=SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
            T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
            T2X=T2X/T2
            T2Y=T2Y/T2
            T2Z=T2Z/T2
            ET2X=T1Z*UNYH0(I,J)-T1Y*UNZH0(I,J)
            ET2Y=T1X*UNZH0(I,J)-T1Z*UNXH0(I,J)
            ET2Z=T1Y*UNXH0(I,J)-T1X*UNYH0(I,J)
            ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
            ET2X=ET2X/ET2
            ET2Y=ET2Y/ET2
            ET2Z=ET2Z/ET2
            ET2XH(I,J)=ET2X
            ET2YH(I,J)=ET2Y
            ET2ZH(I,J)=ET2Z
            ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
            ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
            A2H(I,J)=AA2/ET22
            B2H(I,J)=(BB2-ET12*AA1)/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-BB1*ET12/ET22
            E2H(I,J)=-CC1*ET12/ET22
         END IF !(I == 1)
         IF ((I > 1).AND.(I <= NHU)) THEN
            SAX=AT2XH(I,NHT-1)
            SAY=AT2YH(I,NHT-1)
            SAZ=AT2ZH(I,NHT-1)
            SBX=AT2XH(I,NHT  )
            SBY=AT2YH(I,NHT  )
            SBZ=AT2ZH(I,NHT  )
            SCX=AT2XH(I,NHT+1)
            SCY=AT2YH(I,NHT+1)
            SCZ=AT2ZH(I,NHT+1)
            SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
            SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
            SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
            SA =SA1+SB1
            SB =SB1+SC1
            T2X=SBX
            T2Y=SBY
            T2Z=SBZ
            SP =SA*SB
            SS =SA+SB
            AA2=-SB/(SA*SS)
            BB2=(SB-SA)/SP
            CC2=SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
            T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
            T2X=T2X/T2
            T2Y=T2Y/T2
            T2Z=T2Z/T2
            ET2X=T1Z*UNYH0(I,J)-T1Y*UNZH0(I,J)
            ET2Y=T1X*UNZH0(I,J)-T1Z*UNXH0(I,J)
            ET2Z=T1Y*UNXH0(I,J)-T1X*UNYH0(I,J)
            ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
            ET2X=ET2X/ET2
            ET2Y=ET2Y/ET2
            ET2Z=ET2Z/ET2
            ET2XH(I,J)=ET2X
            ET2YH(I,J)=ET2Y
            ET2ZH(I,J)=ET2Z
            ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
            ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
            A2H(I,J)=AA2/ET22
            B2H(I,J)=(BB2-ET12*BB1)/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-AA1*ET12/ET22
            E2H(I,J)=-CC1*ET12/ET22
         END IF !((I > 1).AND.(I <= NHU))
         IF ((I > NHU).AND.(I < NH2)) THEN
            SAX=AT2XH(I,NHT-2)
            SAY=AT2YH(I,NHT-2)
            SAZ=AT2ZH(I,NHT-2)
            SBX=AT2XH(I,NHT-1)
            SBY=AT2YH(I,NHT-1)
            SBZ=AT2ZH(I,NHT-1)
            SCX=AT2XH(I,NHT  )
            SCY=AT2YH(I,NHT  )
            SCZ=AT2ZH(I,NHT  )
            SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
            SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
            SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
            SA =SA1+SB1
            SB =SB1+SC1
            T2X=SBX
            T2Y=SBY
            T2Z=SBZ
            SP =SA*SB
            SS =SA+SB
            AA2=SB/(SA*SS)
            BB2=-SS/SP
            CC2=(SA+2.D0*SB)/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
            T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
            T2X=T2X/T2
            T2Y=T2Y/T2
            T2Z=T2Z/T2
            ET2X=T1Z*UNYH0(I,J)-T1Y*UNZH0(I,J)
            ET2Y=T1X*UNZH0(I,J)-T1Z*UNXH0(I,J)
            ET2Z=T1Y*UNXH0(I,J)-T1X*UNYH0(I,J)
            ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
            ET2X=ET2X/ET2
            ET2Y=ET2Y/ET2
            ET2Z=ET2Z/ET2
            ET2XH(I,J)=ET2X
            ET2YH(I,J)=ET2Y
            ET2ZH(I,J)=ET2Z
            ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
            ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
            A2H(I,J)=AA2/ET22
            B2H(I,J)=BB2/ET22
            C2H(I,J)=(CC2-ET12*BB1)/ET22
            D2H(I,J)=-AA1*ET12/ET22
            E2H(I,J)=-CC1*ET12/ET22
         END IF !((I > NHU).AND.(I < NH2))
         IF ((I >= NH2).AND.(I < NHX).AND.(NHU /= 0)) THEN
            SAX=AT2XH(I,NHT-2)
            SAY=AT2YH(I,NHT-2)
            SAZ=AT2ZH(I,NHT-2)
            SBX=AT2XH(I,NHT-1)
            SBY=AT2YH(I,NHT-1)
            SBZ=AT2ZH(I,NHT-1)
            SCX=AT2XH(I,NHT  )
            SCY=AT2YH(I,NHT  )
            SCZ=AT2ZH(I,NHT  )
            SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
            SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
            SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
            SA =SA1+SB1
            SB =SB1+SC1
            T2X=SBX
            T2Y=SBY
            T2Z=SBZ
            SP =SA*SB
            SS =SA+SB
            AA2=SB/(SA*SS)
            BB2=-SS/SP
            CC2=(SA+2.D0*SB)/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
            T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
            T2X=T2X/T2
            T2Y=T2Y/T2
            T2Z=T2Z/T2
            ET2X=T1Z*UNYH0(I,J)-T1Y*UNZH0(I,J)
            ET2Y=T1X*UNZH0(I,J)-T1Z*UNXH0(I,J)
            ET2Z=T1Y*UNXH0(I,J)-T1X*UNYH0(I,J)
            ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
            ET2X=ET2X/ET2
            ET2Y=ET2Y/ET2
            ET2Z=ET2Z/ET2
            ET2XH(I,J)=ET2X
            ET2YH(I,J)=ET2Y
            ET2ZH(I,J)=ET2Z
            ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
            ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
            A2H(I,J)=AA2/ET22
            B2H(I,J)=BB2/ET22
            C2H(I,J)=(CC2-ET12*BB1)/ET22
            D2H(I,J)=-AA1*ET12/ET22
            E2H(I,J)=-CC1*ET12/ET22
         END IF !((I >= NH2).AND.(I < NHX).AND.(NHU /= 0))
         IF ((I == NHX).AND.(NHU == 0)) THEN
            SAX=AT2XH(NHX,NHT-1)
            SAY=AT2YH(NHX,NHT-1)
            SAZ=AT2ZH(NHX,NHT-1)
            SBX=AT2XH(NHX,NHT  )
            SBY=AT2YH(NHX,NHT  )
            SBZ=AT2ZH(NHX,NHT  )
            SCX=AT2XH(NHX,NHT+1)
            SCY=AT2YH(NHX,NHT+1)
            SCZ=AT2ZH(NHX,NHT+1)
            SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
            SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
            SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
            SA =SA1+SB1
            SB =SB1+SC1
            T2X=SBX
            T2Y=SBY
            T2Z=SBZ
            SP =SA*SB
            SS =SA+SB
            AA2=-SB/(SA*SS)
            BB2=(SB-SA)/SP
            CC2=SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
            T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
            T2X=T2X/T2
            T2Y=T2Y/T2
            T2Z=T2Z/T2
            ET2X=T1Z*UNYH0(I,J)-T1Y*UNZH0(I,J)
            ET2Y=T1X*UNZH0(I,J)-T1Z*UNXH0(I,J)
            ET2Z=T1Y*UNXH0(I,J)-T1X*UNYH0(I,J)
            ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
            ET2X=ET2X/ET2
            ET2Y=ET2Y/ET2
            ET2Z=ET2Z/ET2
            ET2XH(I,J)=ET2X
            ET2YH(I,J)=ET2Y
            ET2ZH(I,J)=ET2Z
            ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
            ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
            A2H(I,J)=AA2/ET22
            B2H(I,J)=(BB2-ET12*CC1)/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-BB1*ET12/ET22
            E2H(I,J)=-AA1*ET12/ET22
         END IF !((I == NHX).AND.(NHU == 0))
         IF ((I == NHX).AND.(NHU /= 0)) THEN
            SAX=AT2XH(NHX,NHT-2)
            SAY=AT2YH(NHX,NHT-2)
            SAZ=AT2ZH(NHX,NHT-2)
            SBX=AT2XH(NHX,NHT-1)
            SBY=AT2YH(NHX,NHT-1)
            SBZ=AT2ZH(NHX,NHT-1)
            SCX=AT2XH(NHX,NHT  )
            SCY=AT2YH(NHX,NHT  )
            SCZ=AT2ZH(NHX,NHT  )
            SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
            SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
            SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
            SA =SA1+SB1
            SB =SB1+SC1
            T2X=SBX
            T2Y=SBY
            T2Z=SBZ
            SP =SA*SB
            SS =SA+SB
            AA2=SB/(SA*SS)
            BB2=-SS/SP
            CC2=(SA+2.D0*SB)/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
            T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
            T2X=T2X/T2
            T2Y=T2Y/T2
            T2Z=T2Z/T2
            ET2X=T1Z*UNYH0(I,J)-T1Y*UNZH0(I,J)
            ET2Y=T1X*UNZH0(I,J)-T1Z*UNXH0(I,J)
            ET2Z=T1Y*UNXH0(I,J)-T1X*UNYH0(I,J)
            ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
            ET2X=ET2X/ET2
            ET2Y=ET2Y/ET2
            ET2Z=ET2Z/ET2
            ET2XH(I,J)=ET2X
            ET2YH(I,J)=ET2Y
            ET2ZH(I,J)=ET2Z
            ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
            ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
            A2H(I,J)=AA2/ET22
            B2H(I,J)=BB2/ET22
            C2H(I,J)=(CC2-ET12*CC1)/ET22
            D2H(I,J)=-BB1*ET12/ET22
            E2H(I,J)=-AA1*ET12/ET22
         END IF !((I == NHX).AND.(NHU /= 0))
      END IF !(J == NHT)
!-----------------------------------------------------------------------------------------------!
      IF (J == NHT1) THEN
         IF (I == 1) THEN
            SAX=AT2XH(1,NHT1-1)
            SAY=AT2YH(1,NHT1-1)
            SAZ=AT2ZH(1,NHT1-1)
            SBX=AT2XH(1,NHT1  )
            SBY=AT2YH(1,NHT1  )
            SBZ=AT2ZH(1,NHT1  )
            SCX=AT2XH(1,NHT1+1)
            SCY=AT2YH(1,NHT1+1)
            SCZ=AT2ZH(1,NHT1+1)
            SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
            SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
            SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
            SA =SA1+SB1
            SB =SB1+SC1
            SP =SA*SB
            SS =SA+SB
            T2X=SBX
            T2Y=SBY
            T2Z=SBZ
            AA2=-SB/(SA*SS)
            BB2=(SB-SA)/SP
            CC2= SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
            T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
            T2X=T2X/T2
            T2Y=T2Y/T2
            T2Z=T2Z/T2
            ET2X=T1Z*UNYH0(I,J)-T1Y*UNZH0(I,J)
            ET2Y=T1X*UNZH0(I,J)-T1Z*UNXH0(I,J)
            ET2Z=T1Y*UNXH0(I,J)-T1X*UNYH0(I,J)
            ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
            ET2X=ET2X/ET2
            ET2Y=ET2Y/ET2
            ET2Z=ET2Z/ET2
            ET2XH(I,J)=ET2X
            ET2YH(I,J)=ET2Y
            ET2ZH(I,J)=ET2Z
            ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
            ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
            A2H(I,J)=AA2/ET22
            B2H(I,J)=(BB2-ET12*AA1)/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-BB1*ET12/ET22
            E2H(I,J)=-CC1*ET12/ET22
         END IF !(I == 1)
         IF ((I > 1).AND.(I <= NHU)) THEN
            SAX=AT2XH(I,NHT1-1)
            SAY=AT2YH(I,NHT1-1)
            SAZ=AT2ZH(I,NHT1-1)
            SBX=AT2XH(I,NHT1  )
            SBY=AT2YH(I,NHT1  )
            SBZ=AT2ZH(I,NHT1  )
            SCX=AT2XH(I,NHT1+1)
            SCY=AT2YH(I,NHT1+1)
            SCZ=AT2ZH(I,NHT1+1)
            SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
            SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
            SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
            SA =SA1+SB1
            SB =SB1+SC1
            SP =SA*SB
            SS =SA+SB
            T2X=SBX
            T2Y=SBY
            T2Z=SBZ
            AA2=-SB/(SA*SS)
            BB2=(SB-SA)/SP
            CC2= SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
            T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
            T2X=T2X/T2
            T2Y=T2Y/T2
            T2Z=T2Z/T2
            ET2X=T1Z*UNYH0(I,J)-T1Y*UNZH0(I,J)
            ET2Y=T1X*UNZH0(I,J)-T1Z*UNXH0(I,J)
            ET2Z=T1Y*UNXH0(I,J)-T1X*UNYH0(I,J)
            ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
            ET2X=ET2X/ET2
            ET2Y=ET2Y/ET2
            ET2Z=ET2Z/ET2
            ET2XH(I,J)=ET2X
            ET2YH(I,J)=ET2Y
            ET2ZH(I,J)=ET2Z
            ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
            ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
            A2H(I,J)=AA2/ET22
            B2H(I,J)=(BB2-ET12*BB1)/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-AA1*ET12/ET22
            E2H(I,J)=-CC1*ET12/ET22
         END IF !((I > 1).AND.(I <= NHU))
         IF ((I > NHU).AND.(I < NH2)) THEN
            SAX=AT2XH(I,NHT1  )
            SAY=AT2YH(I,NHT1  )
            SAZ=AT2ZH(I,NHT1  )
            SBX=AT2XH(I,NHT1+1)
            SBY=AT2YH(I,NHT1+1)
            SBZ=AT2ZH(I,NHT1+1)
            SCX=AT2XH(I,NHT1+2)
            SCY=AT2YH(I,NHT1+2)
            SCZ=AT2ZH(I,NHT1+2)
            SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
            SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
            SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
            SA =SA1+SB1
            SB =SB1+SC1
            SP =SA*SB
            SS =SA+SB
            T2X=SBX
            T2Y=SBY
            T2Z=SBZ
            AA2=-(2.D0*SA+SB)/(SA*SS)
            BB2=SS/SP
            CC2=-SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
            T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
            T2X=T2X/T2
            T2Y=T2Y/T2
            T2Z=T2Z/T2
            ET2X=T1Z*UNYH0(I,J)-T1Y*UNZH0(I,J)
            ET2Y=T1X*UNZH0(I,J)-T1Z*UNXH0(I,J)
            ET2Z=T1Y*UNXH0(I,J)-T1X*UNYH0(I,J)
            ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
            ET2X=ET2X/ET2
            ET2Y=ET2Y/ET2
            ET2Z=ET2Z/ET2
            ET2XH(I,J)=ET2X
            ET2YH(I,J)=ET2Y
            ET2ZH(I,J)=ET2Z
            ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
            ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
            A2H(I,J)=(AA2-ET12*BB1)/ET22
            B2H(I,J)=BB2/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-AA1*ET12/ET22
            E2H(I,J)=-CC1*ET12/ET22
         END IF !((I > NHU).AND.(I < NH2))
         IF ((I >= NH2).AND.(I < NHX).AND.(NHU /= 0)) THEN
            SAX=AT2XH(I,NHT1  )
            SAY=AT2YH(I,NHT1  )
            SAZ=AT2ZH(I,NHT1  )
            SBX=AT2XH(I,NHT1+1)
            SBY=AT2YH(I,NHT1+1)
            SBZ=AT2ZH(I,NHT1+1)
            SCX=AT2XH(I,NHT1+2)
            SCY=AT2YH(I,NHT1+2)
            SCZ=AT2ZH(I,NHT1+2)
            SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
            SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
            SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
            SA =SA1+SB1
            SB =SB1+SC1
            SP =SA*SB
            SS =SA+SB
            T2X=SBX
            T2Y=SBY
            T2Z=SBZ
            AA2=-(2.D0*SA+SB)/(SA*SS)
            BB2=SS/SP
            CC2=-SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
            T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
            T2X=T2X/T2
            T2Y=T2Y/T2
            T2Z=T2Z/T2
            ET2X=T1Z*UNYH0(I,J)-T1Y*UNZH0(I,J)
            ET2Y=T1X*UNZH0(I,J)-T1Z*UNXH0(I,J)
            ET2Z=T1Y*UNXH0(I,J)-T1X*UNYH0(I,J)
            ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
            ET2X=ET2X/ET2
            ET2Y=ET2Y/ET2
            ET2Z=ET2Z/ET2
            ET2XH(I,J)=ET2X
            ET2YH(I,J)=ET2Y
            ET2ZH(I,J)=ET2Z
            ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
            ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
            A2H(I,J)=(AA2-ET12*BB1)/ET22
            B2H(I,J)=BB2/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-AA1*ET12/ET22
            E2H(I,J)=-CC1*ET12/ET22
         END IF !((I >= NH2).AND.(I < NHX).AND.(NHU /= 0))
         IF ((I == NHX).AND.(NHU == 0)) THEN
            SAX=AT2XH(NHX,NHT1-1)
            SAY=AT2YH(NHX,NHT1-1)
            SAZ=AT2ZH(NHX,NHT1-1)
            SBX=AT2XH(NHX,NHT1  )
            SBY=AT2YH(NHX,NHT1  )
            SBZ=AT2ZH(NHX,NHT1  )
            SCX=AT2XH(NHX,NHT1+1)
            SCY=AT2YH(NHX,NHT1+1)
            SCZ=AT2ZH(NHX,NHT1+1)
            SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
            SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
            SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
            SA =SA1+SB1
            SB =SB1+SC1
            SP =SA*SB
            SS =SA+SB
            T2X=SBX
            T2Y=SBY
            T2Z=SBZ
            AA2=-SB/(SA*SS)
            BB2=(SB-SA)/SP
            CC2= SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
            T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
            T2X=T2X/T2
            T2Y=T2Y/T2
            T2Z=T2Z/T2
            ET2X=T1Z*UNYH0(I,J)-T1Y*UNZH0(I,J)
            ET2Y=T1X*UNZH0(I,J)-T1Z*UNXH0(I,J)
            ET2Z=T1Y*UNXH0(I,J)-T1X*UNYH0(I,J)
            ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
            ET2X=ET2X/ET2
            ET2Y=ET2Y/ET2
            ET2Z=ET2Z/ET2
            ET2XH(I,J)=ET2X
            ET2YH(I,J)=ET2Y
            ET2ZH(I,J)=ET2Z
            ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
            ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
            A2H(I,J)=AA2/ET22
            B2H(I,J)=(BB2-ET12*CC1)/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-BB1*ET12/ET22
            E2H(I,J)=-AA1*ET12/ET22
         END IF !((I == NHX).AND.(NHU == 0))
         IF ((I == NHX).AND.(NHU /= 0)) THEN
            SAX=AT2XH(NHX,NHT1  )
            SAY=AT2YH(NHX,NHT1  )
            SAZ=AT2ZH(NHX,NHT1  )
            SBX=AT2XH(NHX,NHT1+1)
            SBY=AT2YH(NHX,NHT1+1)
            SBZ=AT2ZH(NHX,NHT1+1)
            SCX=AT2XH(NHX,NHT1+2)
            SCY=AT2YH(NHX,NHT1+2)
            SCZ=AT2ZH(NHX,NHT1+2)
            SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
            SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
            SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
            SA =SA1+SB1
            SB =SB1+SC1
            SP =SA*SB
            SS =SA+SB
            T2X=SBX
            T2Y=SBY
            T2Z=SBZ
            AA2=-(2.D0*SA+SB)/(SA*SS)
            BB2=SS/SP
            CC2=-SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
            T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
            T2X=T2X/T2
            T2Y=T2Y/T2
            T2Z=T2Z/T2
            ET2X=T1Z*UNYH0(I,J)-T1Y*UNZH0(I,J)
            ET2Y=T1X*UNZH0(I,J)-T1Z*UNXH0(I,J)
            ET2Z=T1Y*UNXH0(I,J)-T1X*UNYH0(I,J)
            ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
            ET2X=ET2X/ET2
            ET2Y=ET2Y/ET2
            ET2Z=ET2Z/ET2
            ET2XH(I,J)=ET2X
            ET2YH(I,J)=ET2Y
            ET2ZH(I,J)=ET2Z
            ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
            ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
            A2H(I,J)=(AA2-ET12*CC1)/ET22
            B2H(I,J)=BB2/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-BB1*ET12/ET22
            E2H(I,J)=-AA1*ET12/ET22
         END IF !((I == NHX).AND.(NHU /= 0))
      END IF !(J == NHT1)
!-----------------------------------------------------------------------------------------------!
      IF ((J /= 1).AND.(J /= NHTP).AND.(J /= NHT).AND.(J /= NHT1)) THEN
         SAX=AT2XH(I,J-1)
         SAY=AT2YH(I,J-1)
         SAZ=AT2ZH(I,J-1)
         SBX=AT2XH(I,J  )
         SBY=AT2YH(I,J  )
         SBZ=AT2ZH(I,J  )
         SCX=AT2XH(I,J+1)
         SCY=AT2YH(I,J+1)
         SCZ=AT2ZH(I,J+1)
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA =SA1+SB1
         SB =SB1+SC1
         SP =SA*SB
         SS =SA+SB
         T2X=SBX
         T2Y=SBY
         T2Z=SBZ
         AA2=-SB/(SA*SS)
         BB2=(SB-SA)/SP
         CC2= SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!

!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
         T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
         T2X=T2X/T2
         T2Y=T2Y/T2
         T2Z=T2Z/T2
         ET2X=T1Z*UNYH0(I,J)-T1Y*UNZH0(I,J)
         ET2Y=T1X*UNZH0(I,J)-T1Z*UNXH0(I,J)
         ET2Z=T1Y*UNXH0(I,J)-T1X*UNYH0(I,J)
         ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
         ET2X=ET2X/ET2
         ET2Y=ET2Y/ET2
         ET2Z=ET2Z/ET2
         ET2XH(I,J)=ET2X
         ET2YH(I,J)=ET2Y
         ET2ZH(I,J)=ET2Z
         ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
         ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
         IF (I == 1) THEN
            A2H(I,J)=AA2/ET22
            B2H(I,J)=(BB2-ET12*AA1)/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-BB1*ET12/ET22
            E2H(I,J)=-CC1*ET12/ET22
         END IF !(I == 1)
         IF (I == NHX) THEN
            A2H(I,J)=AA2/ET22
            B2H(I,J)=(BB2-ET12*CC1)/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-BB1*ET12/ET22
            E2H(I,J)=-AA1*ET12/ET22
         END IF !(I == NHX)
         IF ((I /= 1).AND.(I /= NHX)) THEN
            A2H(I,J)=AA2/ET22
            B2H(I,J)=(BB2-ET12*BB1)/ET22
            C2H(I,J)=CC2/ET22
            D2H(I,J)=-AA1*ET12/ET22
            E2H(I,J)=-CC1*ET12/ET22
         END IF !((I /= 1).AND.(I /= NHX))
      END IF !((J /= 1).AND.(J /= NHTP).AND.(J /= NHT).AND.(J /= NHT1))
!-----------------------------------------------------------------------------------------------!
   END DO !I=1,NHX
END DO !J=1,NHTP
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE NUMDIFHUB
!-----------------------------------------------------------------------------------------------!
