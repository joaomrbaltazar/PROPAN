!-----------------------------------------------------------------------------------------------!
!    The subroutine NUMDIFNOZZLE computes the coefficients for finite differencing on the       !
!    nozzle mesh                                                                                !
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
SUBROUTINE NUMDIFNOZZLE
!-----------------------------------------------------------------------------------------------!
!    Created by: Joao Baltazar, IST                                                             !
!    Modified  : 05112013, J. Baltazar, version 1.0                                             !
!    Modified  : 26052014, J. Baltazar, Wake Alignment Module                                   !
!    Modified  : 06102015, J. Baltazar, revision                                                !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J
DOUBLE PRECISION :: SAX,SAY,SAZ,SBX,SBY,SBZ,SCX,SCY,SCZ,SA1,SB1,SC1,SA,SB,SS,SP
DOUBLE PRECISION :: T1,T2,ET2X,ET2Y,ET2Z,ET2,ET12,ET22
DOUBLE PRECISION :: T1X,T1Y,T1Z,T2X,T2Y,T2Z,AA1,BB1,CC1,AA2,BB2,CC2
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
DO J=1,NNTP
   DO I=1,NNXT1
!-----------------------------------------------------------------------------------------------!
!    u1 Derivative                                                                              !
!-----------------------------------------------------------------------------------------------!
      IF (I == 1) THEN
         SAX=AT1XN(NNXT1,J)
         SAY=AT1YN(NNXT1,J)
         SAZ=AT1ZN(NNXT1,J)
         SBX=AT1XN(1    ,J)
         SBY=AT1YN(1    ,J)
         SBZ=AT1ZN(1    ,J)
         SCX=AT1XN(2    ,J)
         SCY=AT1YN(2    ,J)
         SCZ=AT1ZN(2    ,J)
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
         A1N(I,J)=AA1
         B1N(I,J)=BB1
         C1N(I,J)=CC1
!-----------------------------------------------------------------------------------------------!
!!    ELSEIF ((I >= NNU).AND.(I <= NNU+2)) THEN !NNU+1
!!       SAX=AT1XN(I-1,J)
!!       SAY=AT1YN(I-1,J)
!!       SAZ=AT1ZN(I-1,J)
!!       SBX=AT1XN(I  ,J)
!!       SBY=AT1YN(I  ,J)
!!       SBZ=AT1ZN(I  ,J)
!!       SCX=AT1XN(I+1,J)
!!       SCY=AT1YN(I+1,J)
!!       SCZ=AT1ZN(I+1,J)
!!       SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
!!       SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
!!       SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
!!       SA =SA1+SB1
!!       SB =SB1+SC1
!!       T1X=SBX
!!       T1Y=SBY
!!       T1Z=SBZ
!-----------------------------------------------------------------------------------------------!
!!       SP =SA*SB
!!       SS =SA+SB
!!       IF (I == NNU) THEN
!!          AA1=-1.D0/SA
!!          BB1= 1.D0/SA
!!          CC1= 0.D0
!!       ELSEIF (I == (NNU+1)) THEN
!!          AA1=-SB/(SA*SS)
!!          BB1=(SB-SA)/SP
!!          CC1= SA/(SB*SS)
!!       ELSEIF (I == (NNU+2)) THEN
!!          AA1= 0.D0
!!          BB1=-1.D0/SB
!!          CC1= 1.D0/SB
!!       END IF !(I == NNU)
!!       A1N(I,J)=AA1
!!       B1N(I,J)=BB1
!!       C1N(I,J)=CC1
!-----------------------------------------------------------------------------------------------!
!!    ELSEIF ((I >= NN2-1).AND.(I <= NN2+1)) THEN !(NN2)
!!       SAX=AT1XN(I-1,J)
!!       SAY=AT1YN(I-1,J)
!!       SAZ=AT1ZN(I-1,J)
!!       SBX=AT1XN(I  ,J)
!!       SBY=AT1YN(I  ,J)
!!       SBZ=AT1ZN(I  ,J)
!!       SCX=AT1XN(I+1,J)
!!       SCY=AT1YN(I+1,J)
!!       SCZ=AT1ZN(I+1,J)
!!       SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
!!       SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
!!       SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
!!       SA =SA1+SB1
!!       SB =SB1+SC1
!!       T1X=SBX
!!       T1Y=SBY
!!       T1Z=SBZ
!-----------------------------------------------------------------------------------------------!
!!       SP =SA*SB
!!       SS =SA+SB
!!       IF (I == (NN2-1)) THEN
!!          AA1=-1.D0/SA
!!          BB1= 1.D0/SA
!!          CC1= 0.D0
!!       ELSEIF (I == NN2) THEN
!!          AA1=-SB/(SA*SS)
!!          BB1=(SB-SA)/SP
!!          CC1= SA/(SB*SS)
!!       ELSEIF (I == (NN2+1)) THEN
!!          AA1= 0.D0
!!          BB1=-1.D0/SB
!!          CC1= 1.D0/SB
!!       END IF !(I == NNU)
!!       A1N(I,J)=AA1
!!       B1N(I,J)=BB1
!!       C1N(I,J)=CC1
!-----------------------------------------------------------------------------------------------!
      ELSEIF ((I == NNX).OR.(I == NNU).OR.(I == (NN2-1))) THEN
         SAX=AT1XN(I-2,J)
         SAY=AT1YN(I-2,J)
         SAZ=AT1ZN(I-2,J)
         SBX=AT1XN(I-1,J)
         SBY=AT1YN(I-1,J)
         SBZ=AT1ZN(I-1,J)
         SCX=AT1XN(I  ,J)
         SCY=AT1YN(I  ,J)
         SCZ=AT1ZN(I  ,J)
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
         AA1= SB/(SA*SS)
         BB1=-SS/SP
         CC1=(SA+2.D0*SB)/(SB*SS)
         A1N(I,J)=AA1
         B1N(I,J)=BB1
         C1N(I,J)=CC1
!-----------------------------------------------------------------------------------------------!
      ELSEIF ((I == NNX1).OR.(I == (NNU+2)).OR.(I == (NN2+1))) THEN
         SAX=AT1XN(I  ,J)
         SAY=AT1YN(I  ,J)
         SAZ=AT1ZN(I  ,J)
         SBX=AT1XN(I+1,J)
         SBY=AT1YN(I+1,J)
         SBZ=AT1ZN(I+1,J)
         SCX=AT1XN(I+2,J)
         SCY=AT1YN(I+2,J)
         SCZ=AT1ZN(I+2,J)
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
         A1N(I,J)=AA1
         B1N(I,J)=BB1
         C1N(I,J)=CC1
!-----------------------------------------------------------------------------------------------!
      ELSEIF (I == NNXT1) THEN
         SAX=AT1XN(NNXT1-1,J)
         SAY=AT1YN(NNXT1-1,J)
         SAZ=AT1ZN(NNXT1-1,J)
         SBX=AT1XN(NNXT1  ,J)
         SBY=AT1YN(NNXT1  ,J)
         SBZ=AT1ZN(NNXT1  ,J)
         SCX=AT1XN(1      ,J)
         SCY=AT1YN(1      ,J)
         SCZ=AT1ZN(1      ,J)
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
         A1N(I,J)=AA1
         B1N(I,J)=BB1
         C1N(I,J)=CC1
!-----------------------------------------------------------------------------------------------!
      ELSE !(I)
         SAX=AT1XN(I-1,J)
         SAY=AT1YN(I-1,J)
         SAZ=AT1ZN(I-1,J)
         SBX=AT1XN(I  ,J)
         SBY=AT1YN(I  ,J)
         SBZ=AT1ZN(I  ,J)
         SCX=AT1XN(I+1,J)
         SCY=AT1YN(I+1,J)
         SCZ=AT1ZN(I+1,J)
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
         A1N(I,J)=AA1
         B1N(I,J)=BB1
         C1N(I,J)=CC1
      END IF !(I)
!-----------------------------------------------------------------------------------------------!
!    Compute t1                                                                                 !
!-----------------------------------------------------------------------------------------------!
      T1 =DSQRT(T1X*T1X+T1Y*T1Y+T1Z*T1Z)
      T1X=T1X/T1
      T1Y=T1Y/T1
      T1Z=T1Z/T1
      ET1XN(I,J)=T1X
      ET1YN(I,J)=T1Y
      ET1ZN(I,J)=T1Z
!-----------------------------------------------------------------------------------------------!
!    u2 Derivative                                                                              !
!-----------------------------------------------------------------------------------------------!
      IF (J == 1) THEN
         SAX=AT2XN(I,NNTP)
         SAY=AT2YN(I,NNTP)
         SAZ=AT2ZN(I,NNTP)
         SBX=AT2XN(I,1   )
         SBY=AT2YN(I,1   )
         SBZ=AT2ZN(I,1   )
         SCX=AT2XN(I,2   )
         SCY=AT2YN(I,2   )
         SCZ=AT2ZN(I,2   )
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA =SA1+SB1
         SB =SB1+SC1
         T2X=SBX
         T2Y=SBY
         T2Z=SBZ
!-----------------------------------------------------------------------------------------------!
         SP =SA*SB
         SS =SA+SB
         AA2=-SB/(SA*SS)
         BB2=(SB-SA)/SP
         CC2= SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
         T2  =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
         T2X =T2X/T2
         T2Y =T2Y/T2
         T2Z =T2Z/T2
         ET2X=T1Z*UNYN0(I,J)-T1Y*UNZN0(I,J)
         ET2Y=T1X*UNZN0(I,J)-T1Z*UNXN0(I,J)
         ET2Z=T1Y*UNXN0(I,J)-T1X*UNYN0(I,J)
         ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
         ET2X=ET2X/ET2
         ET2Y=ET2Y/ET2
         ET2Z=ET2Z/ET2
         ET2XN(I,J)=ET2X
         ET2YN(I,J)=ET2Y
         ET2ZN(I,J)=ET2Z
         ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
         ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
         IF ((I == NNX).OR.(I == NNU).OR.(I == (NN2-1))) THEN
            A2N(I,J)=AA2/ET22
            B2N(I,J)=(BB2-ET12*CC1)/ET22
            C2N(I,J)=CC2/ET22
            D2N(I,J)=-BB1*ET12/ET22
            E2N(I,J)=-AA1*ET12/ET22
         ELSEIF ((I == NNX1).OR.(I == (NNU+2)).OR.(I == (NN2+1))) THEN
            A2N(I,J)=AA2/ET22
            B2N(I,J)=(BB2-ET12*AA1)/ET22
            C2N(I,J)=CC2/ET22
            D2N(I,J)=-BB1*ET12/ET22
            E2N(I,J)=-CC1*ET12/ET22
         ELSE !(I)
            A2N(I,J)=AA2/ET22
            B2N(I,J)=(BB2-ET12*BB1)/ET22
            C2N(I,J)=CC2/ET22
            D2N(I,J)=-AA1*ET12/ET22
            E2N(I,J)=-CC1*ET12/ET22
         END IF !(I)
!-----------------------------------------------------------------------------------------------!
      ELSEIF (J == NNT) THEN
         SAX=AT2XN(I,NNT-1)
         SAY=AT2YN(I,NNT-1)
         SAZ=AT2ZN(I,NNT-1)
         SBX=AT2XN(I,NNT  )
         SBY=AT2YN(I,NNT  )
         SBZ=AT2ZN(I,NNT  )
         SCX=AT2XN(I,NNT+1)
         SCY=AT2YN(I,NNT+1)
         SCZ=AT2ZN(I,NNT+1)
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA =SA1+SB1
         SB =SB1+SC1
         T2X=SBX
         T2Y=SBY
         T2Z=SBZ
!-----------------------------------------------------------------------------------------------!
         SP =SA*SB
         SS =SA+SB
         AA2=-SB/(SA*SS)
         BB2=(SB-SA)/SP
         CC2= SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!
         IF (ISTRIP == 1) THEN
            IF (((I >= NNU).AND.(I < NNX1)).or.((i > nnx).and.(i <= nntp-nnu+1))) THEN
               SAX=AT2XN(I,NNT-2)
               SAY=AT2YN(I,NNT-2)
               SAZ=AT2ZN(I,NNT-2)
               SBX=AT2XN(I,NNT-1)
               SBY=AT2YN(I,NNT-1)
               SBZ=AT2ZN(I,NNT-1)
               SCX=AT2XN(I,NNT  )
               SCY=AT2YN(I,NNT  )
               SCZ=AT2ZN(I,NNT  )
               SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
               SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
               SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
               SA =SA1+SB1
               SB =SB1+SC1
               T2X=SCX
               T2Y=SCY
               T2Z=SCZ
!-----------------------------------------------------------------------------------------------!
               SP =SA*SB
               SS =SA+SB
               AA2=SB/(SA*SS)
               BB2=-SS/SP
               CC2=(SA+2.D0*SB)/(SB*SS)
            END IF !((I >= NNU).AND.(I < NNX1))
         END IF !(ISTRIP == 1)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
         T2  =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
         T2X =T2X/T2
         T2Y =T2Y/T2
         T2Z =T2Z/T2
         ET2X=T1Z*UNYN0(I,J)-T1Y*UNZN0(I,J)
         ET2Y=T1X*UNZN0(I,J)-T1Z*UNXN0(I,J)
         ET2Z=T1Y*UNXN0(I,J)-T1X*UNYN0(I,J)
         ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
         ET2X=ET2X/ET2
         ET2Y=ET2Y/ET2
         ET2Z=ET2Z/ET2
         ET2XN(I,J)=ET2X
         ET2YN(I,J)=ET2Y
         ET2ZN(I,J)=ET2Z
         ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
         ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
         IF ((I == NNX).OR.(I == NNU).OR.(I == (NN2-1))) THEN
            A2N(I,J)=AA2/ET22
            B2N(I,J)=(BB2-ET12*CC1)/ET22
            C2N(I,J)=CC2/ET22
            D2N(I,J)=-BB1*ET12/ET22
            E2N(I,J)=-AA1*ET12/ET22
         ELSEIF ((I == NNX1).OR.(I == (NNU+2)).OR.(I == (NN2+1))) THEN
            A2N(I,J)=AA2/ET22
            B2N(I,J)=(BB2-ET12*AA1)/ET22
            C2N(I,J)=CC2/ET22
            D2N(I,J)=-BB1*ET12/ET22
            E2N(I,J)=-CC1*ET12/ET22
         ELSE !(I)
            A2N(I,J)=AA2/ET22
            B2N(I,J)=(BB2-ET12*BB1)/ET22
            C2N(I,J)=CC2/ET22
            D2N(I,J)=-AA1*ET12/ET22
            E2N(I,J)=-CC1*ET12/ET22
         END IF !(I)
         IF (ISTRIP == 1) THEN
            IF (((I > NNU).AND.(I < NNX)).or.((i > nnx1).and.(i <= nntp-nnu+1))) THEN
               A2N(I,J)=AA2/ET22
               B2N(I,J)=BB2/ET22
               C2N(I,J)=(CC2-ET12*BB1)/ET22
               D2N(I,J)=-AA1*ET12/ET22
               E2N(I,J)=-CC1*ET12/ET22
            END IF !((I > NNU).AND.(I < NNX))
            IF ((I == NNX).OR.(I == NNU).OR.(I == (NN2-1))) THEN
               A2N(I,J)=AA2/ET22
               B2N(I,J)=BB2/ET22
               C2N(I,J)=(CC2-ET12*CC1)/ET22
               D2N(I,J)=-BB1*ET12/ET22
               E2N(I,J)=-AA1*ET12/ET22
            END IF !(I == NNX)
            IF ((I == (NNU+2)).OR.(I == (NN2+1)).or.(i == nnx1)) THEN
               A2N(I,J)=AA2/ET22
               B2N(I,J)=BB2/ET22
               C2N(I,J)=(CC2-ET12*AA1)/ET22
               D2N(I,J)=-BB1*ET12/ET22
               E2N(I,J)=-CC1*ET12/ET22
            END IF !((I == (NNU+2)).OR.(I == (NN2+1)))
         END IF !(ISTRIP == 1)
!-----------------------------------------------------------------------------------------------!
      ELSEIF (J == NNT+1) THEN
         SAX=AT2XN(I,NNT  )
         SAY=AT2YN(I,NNT  )
         SAZ=AT2ZN(I,NNT  )
         SBX=AT2XN(I,NNT+1)
         SBY=AT2YN(I,NNT+1)
         SBZ=AT2ZN(I,NNT+1)
         SCX=AT2XN(I,NNT+2)
         SCY=AT2YN(I,NNT+2)
         SCZ=AT2ZN(I,NNT+2)
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA =SA1+SB1
         SB =SB1+SC1
         T2X=SBX
         T2Y=SBY
         T2Z=SBZ
!-----------------------------------------------------------------------------------------------!
         SP =SA*SB
         SS =SA+SB
         AA2=-SB/(SA*SS)
         BB2=(SB-SA)/SP
         CC2= SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!
         IF (ISTRIP == 1) THEN
            IF (((I >= NNU).AND.(I < NNX1)).or.((i > nnx).and.(i <= nntp-nnu+1))) THEN
               SAX=AT2XN(I,NNT+1)
               SAY=AT2YN(I,NNT+1)
               SAZ=AT2ZN(I,NNT+1)
               SBX=AT2XN(I,NNT+2)
               SBY=AT2YN(I,NNT+2)
               SBZ=AT2ZN(I,NNT+2)
               SCX=AT2XN(I,NNT+3)
               SCY=AT2YN(I,NNT+3)
               SCZ=AT2ZN(I,NNT+3)
               SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
               SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
               SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
               SA =SA1+SB1
               SB =SB1+SC1
               T2X=SAX
               T2Y=SAY
               T2Z=SAZ
!-----------------------------------------------------------------------------------------------!
               SP =SA*SB
               SS =SA+SB
               AA2=-(2.D0*SA+SB)/(SA*SS)
               BB2=SS/SP
               CC2=-SA/(SB*SS)
            END IF !((I >= NNU).AND.(I < NNX1))
         END IF !(ISTRIP == 1)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
         T2  =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
         T2X =T2X/T2
         T2Y =T2Y/T2
         T2Z =T2Z/T2
         ET2X=T1Z*UNYN0(I,J)-T1Y*UNZN0(I,J)
         ET2Y=T1X*UNZN0(I,J)-T1Z*UNXN0(I,J)
         ET2Z=T1Y*UNXN0(I,J)-T1X*UNYN0(I,J)
         ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
         ET2X=ET2X/ET2
         ET2Y=ET2Y/ET2
         ET2Z=ET2Z/ET2
         ET2XN(I,J)=ET2X
         ET2YN(I,J)=ET2Y
         ET2ZN(I,J)=ET2Z
         ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
         ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
         IF ((I == NNX).OR.(I == NNU).OR.(I == (NN2-1))) THEN
            A2N(I,J)=AA2/ET22
            B2N(I,J)=(BB2-ET12*CC1)/ET22
            C2N(I,J)=CC2/ET22
            D2N(I,J)=-BB1*ET12/ET22
            E2N(I,J)=-AA1*ET12/ET22
         ELSEIF ((I == NNX1).OR.(I == (NNU+2)).OR.(I == (NN2+1))) THEN
            A2N(I,J)=AA2/ET22
            B2N(I,J)=(BB2-ET12*AA1)/ET22
            C2N(I,J)=CC2/ET22
            D2N(I,J)=-BB1*ET12/ET22
            E2N(I,J)=-CC1*ET12/ET22
         ELSE !(I)
            A2N(I,J)=AA2/ET22
            B2N(I,J)=(BB2-ET12*BB1)/ET22
            C2N(I,J)=CC2/ET22
            D2N(I,J)=-AA1*ET12/ET22
            E2N(I,J)=-CC1*ET12/ET22
         END IF !(I)
         IF (ISTRIP == 1) THEN
            IF (((I > NNU).AND.(I < NNX)).or.((i > nnx1).and.(i <= nntp-nnu+1))) THEN
               A2N(I,J)=(AA2-ET12*BB1)/ET22
               B2N(I,J)=BB2/ET22
               C2N(I,J)=CC2/ET22
               D2N(I,J)=-AA1*ET12/ET22
               E2N(I,J)=-CC1*ET12/ET22
            END IF !((I > NNU).AND.(I < NNX))
            IF ((I == NNX).OR.(I == NNU).OR.(I == (NN2-1))) THEN
               A2N(I,J)=(AA2-ET12*CC1)/ET22
               B2N(I,J)=BB2/ET22
               C2N(I,J)=CC2/ET22
               D2N(I,J)=-BB1*ET12/ET22
               E2N(I,J)=-AA1*ET12/ET22
            END IF !(I == NNX)
            IF ((I == (NNU+2)).OR.(I == (NN2+1)).or.(i == nnx1)) THEN
               A2N(I,J)=(AA2-ET12*AA1)/ET22
               B2N(I,J)=BB2/ET22
               C2N(I,J)=CC2/ET22
               D2N(I,J)=-BB1*ET12/ET22
               E2N(I,J)=-CC1*ET12/ET22
            END IF !((I == (NNU+2)).OR.(I == (NN2+1)))
         END IF !(ISTRIP == 1)
!-----------------------------------------------------------------------------------------------!
      ELSEIF (J == NNTP) THEN
         SAX=AT2XN(I,NNTP-1)
         SAY=AT2YN(I,NNTP-1)
         SAZ=AT2ZN(I,NNTP-1)
         SBX=AT2XN(I,NNTP  )
         SBY=AT2YN(I,NNTP  )
         SBZ=AT2ZN(I,NNTP  )
         SCX=AT2XN(I,1     )
         SCY=AT2YN(I,1     )
         SCZ=AT2ZN(I,1     )
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA =SA1+SB1
         SB =SB1+SC1
         T2X=SBX
         T2Y=SBY
         T2Z=SBZ
!-----------------------------------------------------------------------------------------------!
         SP =SA*SB
         SS =SA+SB
         AA2=-SB/(SA*SS)
         BB2=(SB-SA)/SP
         CC2= SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
         T2  =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
         T2X =T2X/T2
         T2Y =T2Y/T2
         T2Z =T2Z/T2
         ET2X=T1Z*UNYN0(I,J)-T1Y*UNZN0(I,J)
         ET2Y=T1X*UNZN0(I,J)-T1Z*UNXN0(I,J)
         ET2Z=T1Y*UNXN0(I,J)-T1X*UNYN0(I,J)
         ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
         ET2X=ET2X/ET2
         ET2Y=ET2Y/ET2
         ET2Z=ET2Z/ET2
         ET2XN(I,J)=ET2X
         ET2YN(I,J)=ET2Y
         ET2ZN(I,J)=ET2Z
         ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
         ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
         IF ((I == NNX).OR.(I == NNU).OR.(I == (NN2-1))) THEN
            A2N(I,J)=AA2/ET22
            B2N(I,J)=(BB2-ET12*CC1)/ET22
            C2N(I,J)=CC2/ET22
            D2N(I,J)=-BB1*ET12/ET22
            E2N(I,J)=-AA1*ET12/ET22
         ELSEIF ((I == NNX1).OR.(I == (NNU+2)).OR.(I == (NN2+1))) THEN
            A2N(I,J)=AA2/ET22
            B2N(I,J)=(BB2-ET12*AA1)/ET22
            C2N(I,J)=CC2/ET22
            D2N(I,J)=-BB1*ET12/ET22
            E2N(I,J)=-CC1*ET12/ET22
         ELSE !(I)
            A2N(I,J)=AA2/ET22
            B2N(I,J)=(BB2-ET12*BB1)/ET22
            C2N(I,J)=CC2/ET22
            D2N(I,J)=-AA1*ET12/ET22
            E2N(I,J)=-CC1*ET12/ET22
         END IF !(I)
!-----------------------------------------------------------------------------------------------!
      ELSE !(J)
         SAX=AT2XN(I,J-1)
         SAY=AT2YN(I,J-1)
         SAZ=AT2ZN(I,J-1)
         SBX=AT2XN(I,J  )
         SBY=AT2YN(I,J  )
         SBZ=AT2ZN(I,J  )
         SCX=AT2XN(I,J+1)
         SCY=AT2YN(I,J+1)
         SCZ=AT2ZN(I,J+1)
         SA1=DSQRT(SAX*SAX+SAY*SAY+SAZ*SAZ)
         SB1=DSQRT(SBX*SBX+SBY*SBY+SBZ*SBZ)
         SC1=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)
         SA =SA1+SB1
         SB =SB1+SC1
         T2X=SBX
         T2Y=SBY
         T2Z=SBZ
!-----------------------------------------------------------------------------------------------!
         SP =SA*SB
         SS =SA+SB
         AA2=-SB/(SA*SS)
         BB2=(SB-SA)/SP
         CC2= SA/(SB*SS)
!-----------------------------------------------------------------------------------------------!
!    Compute Coefficients                                                                       !
!-----------------------------------------------------------------------------------------------!
         T2  =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
         T2X =T2X/T2
         T2Y =T2Y/T2
         T2Z =T2Z/T2
         ET2X=T1Z*UNYN0(I,J)-T1Y*UNZN0(I,J)
         ET2Y=T1X*UNZN0(I,J)-T1Z*UNXN0(I,J)
         ET2Z=T1Y*UNXN0(I,J)-T1X*UNYN0(I,J)
         ET2 =DSQRT(ET2X*ET2X+ET2Y*ET2Y+ET2Z*ET2Z)
         ET2X=ET2X/ET2
         ET2Y=ET2Y/ET2
         ET2Z=ET2Z/ET2
         ET2XN(I,J)=ET2X
         ET2YN(I,J)=ET2Y
         ET2ZN(I,J)=ET2Z
         ET12=T1X*T2X+T1Y*T2Y+T1Z*T2Z
         ET22=T2X*ET2X+T2Y*ET2Y+T2Z*ET2Z
         IF ((I == NNX).OR.(I == NNU).OR.(I == (NN2-1))) THEN
            A2N(I,J)=AA2/ET22
            B2N(I,J)=(BB2-ET12*CC1)/ET22
            C2N(I,J)=CC2/ET22
            D2N(I,J)=-BB1*ET12/ET22
            E2N(I,J)=-AA1*ET12/ET22
         ELSEIF ((I == NNX1).OR.(I == (NNU+2)).OR.(I == (NN2+1))) THEN
            A2N(I,J)=AA2/ET22
            B2N(I,J)=(BB2-ET12*AA1)/ET22
            C2N(I,J)=CC2/ET22
            D2N(I,J)=-BB1*ET12/ET22
            E2N(I,J)=-CC1*ET12/ET22
         ELSE !(I)
            A2N(I,J)=AA2/ET22
            B2N(I,J)=(BB2-ET12*BB1)/ET22
            C2N(I,J)=CC2/ET22
            D2N(I,J)=-AA1*ET12/ET22
            E2N(I,J)=-CC1*ET12/ET22
         END IF !(I)
      END IF !(J)
!-----------------------------------------------------------------------------------------------!
   END DO !I=1,NNXT1
END DO !J=1,NNTP
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE NUMDIFNOZZLE
!-----------------------------------------------------------------------------------------------!
