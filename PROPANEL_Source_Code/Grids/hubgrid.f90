!-----------------------------------------------------------------------------------------------!
!    Generate Hub Grid                                                                          !
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
SUBROUTINE HUBGRID
!-----------------------------------------------------------------------------------------------!
!    Created by: 28102013, J. Baltazar, version 1.0                                             !
!    Modified  : 09052014, J. Baltazar, new TH at downstream and Cosine on Theta                !
!                08102014, J. Baltazar, version 1.1                                             !
!-----------------------------------------------------------------------------------------------!
!    Declarations                                                                               !
!-----------------------------------------------------------------------------------------------!
USE PROPANEL_MOD
IMPLICIT NONE
INTEGER :: I,J,K
INTEGER :: NH2,NH3
DOUBLE PRECISION :: XH1,XH2,TH0,TH1,TH2,TH3,XBM,TBM,PHI
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: PH
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: XXG,YYG
!-----------------------------------------------------------------------------------------------!
IF (IH == 1) THEN
!-----------------------------------------------------------------------------------------------!
!    Counters                                                                                   !
!-----------------------------------------------------------------------------------------------!
   NHX  =NHU+NC+NHD
   NHX1 =NHX+1
   NHU1 =NHU+1
   NH2  =NHU1+NC
   NH3  =NH2+NHD
   NHT1 =NHT+1
   NHTT =2*NHT1
   NHTT1=NHTT-1
   NCH  =NC
!-----------------------------------------------------------------------------------------------!
   ALLOCATE(XH(NH3,NHTT),YH(NH3,NHTT),ZH(NH3,NHTT),TH(NH3,NHTT),RH(NH3,NHTT),PH(NH3))
   XH=0.D0
   YH=0.D0
   ZH=0.D0
   TH=0.D0
   RH=0.D0
   PH=0.D0
!-----------------------------------------------------------------------------------------------!
!    Define Coordinates of Leading and Trailing Edge of Blade                                   !
!-----------------------------------------------------------------------------------------------!
   XH1=XP(NC1,1)
   XH2=0.5D0*(XP(1,1)+XP(NCP1,1))
   TH1=TP(NC1,1)
   TH2=0.5D0*(TP(1,1)+TP(NCP1,1))
   PTH=2.D0*PTH !non-dimensional by propeller radius
!-----------------------------------------------------------------------------------------------!
!    Panel Upstream Section of Hub                                                              !
!-----------------------------------------------------------------------------------------------!
   DO I=1,NHU1
!-----------------------------------------------------------------------------------------------!
!    Cosine Distribution on X                                                                   !
!-----------------------------------------------------------------------------------------------!
      PHI=DFLOAT(I-1)*PI/DFLOAT(NHU)
      XH(I,1    )=0.5D0*(XH0+XH1)-0.5D0*(XH1-XH0)*DCOS(PHI)
      XH(I,NHTT1)=XH(I,1)
   END DO !I=1,NHU1
!-----------------------------------------------------------------------------------------------!
!    Interpolate the Pitch                                                                      !
!-----------------------------------------------------------------------------------------------!
   IF (INTERH == 0) THEN
      CALL LININT(NHP,XHP,PTH,NHU1,XH(1:NHU1,1),PH(1:NHU1))
   ELSEIF (INTERH == 1) THEN
      CALL INTK1 (NHP,XHP,PTH,NHU1,XH(1:NHU1,1),PH(1:NHU1))
   ELSEIF (INTERH == 2) THEN
      CALL SPLINT(NHP,XHP,PTH,NHU1,XH(1:NHU1,1),PH(1:NHU1))
   END IF !(INTERH)
!-----------------------------------------------------------------------------------------------!
!    Define Theta Coordinate                                                                    !
!-----------------------------------------------------------------------------------------------!
   TH(NHU1,1    )=TH1
   TH(NHU1,NHTT1)=TH(NHU1,1)+2.D0*PI/DFLOAT(NB)
   DO I=NHU,1,-1
!!    TH(I,1    )=TH0+(XH(I,1)-XH0)*2.D0*PI/PTH
      TH(I,1    )=TH(I+1,1)+(XH(I,1)-XH(I+1,1))*2.D0*PI/PH(I+1)
      TH(I,NHTT1)=TH(I,1)+2.D0*PI/DFLOAT(NB)
   END DO !I=NHU,1,-1
   TH0=TH(1,1)
!-----------------------------------------------------------------------------------------------!
!    Panel Distribution Between Blades                                                          !
!-----------------------------------------------------------------------------------------------!
   DO I=1,NC
      XH(NHU1+I,1    )=XP(NC1+I,1)
      TH(NHU1+I,1    )=TP(NC1+I,1)
      XH(NHU1+I,NHTT1)=XP(NC1-I,1)
      TH(NHU1+I,NHTT1)=TP(NC1-I,1)+2.D0*PI/DFLOAT(NB) 
   END DO !I=1,NC
!-----------------------------------------------------------------------------------------------!
!    Panel Downstream Section of Hub                                                            !
!-----------------------------------------------------------------------------------------------!
   DO I=1,NHD
!-----------------------------------------------------------------------------------------------!
!    Cosine Distribution on X                                                                   !
!-----------------------------------------------------------------------------------------------!
      IF (IHR == 0) THEN
         PHI=DFLOAT(I)/DFLOAT(NHD)*PI/2.D0
         XH(NH2+I,1    )=XH3+(XH2-XH3)*DCOS(PHI)
         XH(NH2+I,NHTT1)=XH(NH2+I,1)
!-----------------------------------------------------------------------------------------------!
!    Interpolate the Pitch                                                                      !
!-----------------------------------------------------------------------------------------------!
         IF (INTERH == 0) THEN
            CALL LININT(NHP,XHP,PTH,1,XH(NH2+I,1),PH(NH2+I))
         ELSEIF (INTERH == 1) THEN
            CALL INTK1 (NHP,XHP,PTH,1,XH(NH2+I,1),PH(NH2+I))
         ELSEIF (INTERH == 2) THEN
            CALL SPLINT(NHP,XHP,PTH,1,XH(NH2+I,1),PH(NH2+I))
         END IF !(INTERH)
         TH(NH2+I,1)=TH(NH2+I-1,1)+(XH(NH2+I,1)-XH(NH2+I-1,1))*2.D0*PI/PH(NH2+I)
!-----------------------------------------------------------------------------------------------!
!    Wake Distribution on X                                                                     !
!-----------------------------------------------------------------------------------------------!
      ELSEIF (IHR == 1) THEN
         XH(NH2+I,1    )=XPW(ISTEP*I+1,1)
         XH(NH2+I,NHTT1)=XPW(ISTEP*I+1,1)
         CALL LININT(NPW1,XPW(:,1),TPW(:,1),1,XH(NH2+I,1),TH(NH2+I,1))
      END IF !(IHR)
   END DO !I=1,NHD
!-----------------------------------------------------------------------------------------------!
!    Define Theta Coordinate                                                                    !
!-----------------------------------------------------------------------------------------------!
   DO I=1,NHD
      TH(NH2+I,NHTT1)=TH(NH2+I,1)+2.D0*PI/DFLOAT(NB)
   END DO !I=1,NHD
   TH3=TH(NH3,1)
!-----------------------------------------------------------------------------------------------!
!    Define Point Distribution on Boundaries 2 and 4 of Domain                                  !
!-----------------------------------------------------------------------------------------------!
   ALPHAHT=ALPHAHT*PI/180.D0
   DO J=1,NHTT1
      XH(1   ,J)=XH0
      XH(NHX1,J)=XH3
!-----------------------------------------------------------------------------------------------!
!    Cosine Distribution on Theta                                                               !
!-----------------------------------------------------------------------------------------------!
      IF (ITHETA == 0) THEN
         PHI    =ALPHAHT+DFLOAT(J-1)/DFLOAT(NHTT1-1)*(PI-2.D0*ALPHAHT)
         TH(1   ,J)=(TH0+PI/DFLOAT(NB))-(PI/DFLOAT(NB))*DCOS(PHI)/DCOS(ALPHAHT)
         TH(NHX1,J)=(TH3+PI/DFLOAT(NB))-(PI/DFLOAT(NB))*DCOS(PHI)/DCOS(ALPHAHT)
!-----------------------------------------------------------------------------------------------!
!    Equidistant Distribution on Theta                                                          !
!-----------------------------------------------------------------------------------------------!
      ELSEIF (ITHETA == 1) THEN
         PHI       =1.D0-2.D0*DFLOAT(J-1)/DFLOAT(NHTT1-1)
         TH(1   ,J)=(TH0+PI/DFLOAT(NB))-(PI/DFLOAT(NB))*PHI
         TH(NHX1,J)=(TH3+PI/DFLOAT(NB))-(PI/DFLOAT(NB))*PHI
      END IF !(ITHETA)
   END DO !J=1,NHTT1
!-----------------------------------------------------------------------------------------------!
!    Generate Grid on the Hub with the Program GRAPE                                            !
!-----------------------------------------------------------------------------------------------!
!    Define Arrays with Boundary Points                                                         !
!-----------------------------------------------------------------------------------------------!
   ALLOCATE(XXG(NHX1,NHTT1),YYG(NHX1,NHTT1))
   XXG=0.D0
   YYG=0.D0
!-----------------------------------------------------------------------------------------------!
   DO I=1,NHX1
      XXG(I,1    )=XH(I,1    )
      YYG(I,1    )=TH(I,1    )
      XXG(I,NHTT1)=XH(I,NHTT1)
      YYG(I,NHTT1)=TH(I,NHTT1)
   END DO !I=1,NHX1
!-----------------------------------------------------------------------------------------------!
   DO J=1,NHTT1
      XXG(1   ,J)=XH(1   ,J)
      YYG(1   ,J)=TH(1   ,J)
      XXG(NHX1,J)=XH(NHX1,J)
      YYG(NHX1,J)=TH(NHX1,J)
   END DO !J=1,NHTT1
!-----------------------------------------------------------------------------------------------!
   CALL GRAPE(NHX1,NHTT1,XXG,YYG,ITERH)
!-----------------------------------------------------------------------------------------------!
   DO I=1,NHX1
!-----------------------------------------------------------------------------------------------!
!    Fill in First Half-Sector                                                                  !
!-----------------------------------------------------------------------------------------------!
      DO J=1,NHT1
         XH(I,J)=XXG(I,J+NHT1-1)
         TH(I,J)=YYG(I,J+NHT1-1)-2.D0*PI/DFLOAT(NB)
!-----------------------------------------------------------------------------------------------!
!    Interpolate the Radius                                                                     !
!-----------------------------------------------------------------------------------------------!
         IF (INTERH == 0) THEN
            CALL LININT(NHI,XHI,RHI,1,XH(I,J),RH(I,J))
         ELSEIF (INTERH == 1) THEN
            CALL INTK1 (NHI,XHI,RHI,1,XH(I,J),RH(I,J))
         ELSEIF (INTERH == 2) THEN
            CALL SPLINT(NHI,XHI,RHI,1,XH(I,J),RH(I,J))
         END IF !(INTERH)
      END DO !J=1,NHT1
!-----------------------------------------------------------------------------------------------!
!    Fill in Second Half-Sector                                                                 !
!-----------------------------------------------------------------------------------------------!
      DO J=1,NHT1
         XH(I,NHT1+J)=XXG(I,J)
         TH(I,NHT1+J)=YYG(I,J)
!-----------------------------------------------------------------------------------------------!
!    Interpolate the Radius                                                                     !
!-----------------------------------------------------------------------------------------------!
         IF (INTERH == 0) THEN
            CALL LININT(NHI,XHI,RHI,1,XH(I,NHT1+J),RH(I,NHT1+J))
         ELSEIF (INTERH == 1) THEN
            CALL INTK1 (NHI,XHI,RHI,1,XH(I,NHT1+J),RH(I,NHT1+J))
         ELSEIF (INTERH == 2) THEN
            CALL SPLINT(NHI,XHI,RHI,1,XH(I,NHT1+J),RH(I,NHT1+J))
         END IF !(INTERH)
      END DO !J=1,NHT1
!-----------------------------------------------------------------------------------------------!
!    Compute Cartesian Coordinates Complete Hub                                                 !
!-----------------------------------------------------------------------------------------------!
      DO J=1,NHTT
         YH(I,J)=RH(I,J)*DCOS(TH(I,J))
         ZH(I,J)=RH(I,J)*DSIN(TH(I,J))
      END DO !J=1,NHTT
   END DO !I=1,NHX1
!-----------------------------------------------------------------------------------------------!
!    Deallocate Variables                                                                       !
!-----------------------------------------------------------------------------------------------!
   DEALLOCATE(PH,XXG,YYG)
END IF !(IH == 1)
!-----------------------------------------------------------------------------------------------!
!    Fictious hub to close root of blade sections                                               !
!-----------------------------------------------------------------------------------------------!
IF (IH == -1) THEN
!-----------------------------------------------------------------------------------------------!
!    Counters                                                                                   !
!-----------------------------------------------------------------------------------------------!
   NHU =0
   NHD =0
   NCH =NC
   NHX =NCH
   NHX1=NHX+1
   NHU1=NHU+1
   NH2 =NHU1+NCH
   NHT =2
   NHT1=NHT+1
   NHT2=NHT1+1
   NHTT=2*NHT1
!-----------------------------------------------------------------------------------------------!
   ALLOCATE(XH(NHX1,NHTT),YH(NHX1,NHTT),ZH(NHX1,NHTT),RH(NHX1,NHTT),TH(NHX1,NHTT))
   XH=0.D0
   YH=0.D0
   ZH=0.D0
   RH=0.D0
   TH=0.D0
!-----------------------------------------------------------------------------------------------!
   DO I=1,NHX1
      XBM=0.5D0*(XP(NC1-I+1,1)+XP(NC1+I-1,1))
      TBM=0.5D0*(TP(NC1-I+1,1)+TP(NC1+I-1,1))
      XH(I,1   )=XBM
      XH(I,NHT1)=XP(NC1-I+1,1)
      DO J=2,NHT
         XH(I,J)=XH(I,1)+DFLOAT(J-1)*(XH(I,NHT1)-XH(I,1))/DFLOAT(NHT)
      END DO !J=2,NHT
      XH(I,NHT2)=XP(NC1+I-1,1)
      XH(I,NHTT)=XBM
      DO J=2,NHT
         XH(I,NHT1+J)=XH(I,NHT2)+DFLOAT(J-1)*(XH(I,NHTT)-XH(I,NHT2))/DFLOAT(NHT)
      END DO !J=2,NHT
      DO J=1,NHTT
         RH(I,J)=RPH
      END DO !J=1,NHTT
!-----------------------------------------------------------------------------------------------!
      TH(I,1   )=TBM
      TH(I,NHT1)=TP(NC1-I+1,1)
      DO J=2,NHT
         TH(I,J)=TH(I,1)+DFLOAT(J-1)*(TH(I,NHT1)-TH(I,1))/DFLOAT(NHT)
      END DO !J=2,NHT
      TH(I,NHT2)=TP(NC1+I-1,1)
      TH(I,NHTT)=TBM
      DO J=2,NHT
         TH(I,NHT1+J)=TH(I,NHT2)+DFLOAT(J-1)*(TH(I,NHTT)-TH(I,NHT2))/DFLOAT(NHT)
      END DO !J=2,NHT
      DO J=1,NHTT
         YH(I,J)=RH(I,J)*DCOS(TH(I,J))
         ZH(I,J)=RH(I,J)*DSIN(TH(I,J))
      END DO !J=1,NHTT
   END DO !I=1,NHX1
!-----------------------------------------------------------------------------------------------!
END IF !(IH == -1)
!-----------------------------------------------------------------------------------------------!
!    Write Hub Grid in Tecplot Format                                                           !
!-----------------------------------------------------------------------------------------------!
K=1
WRITE(20,100) ' ZONE T="HUB',K,'" F=POINT, I=',NHX1,' J=',NHTT
DO J=1,NHT1
   DO I=1,NHX1
      WRITE(20,110) XH(I,J),YH(I,J),ZH(I,J)
   END DO !I=1,NHX1
END DO !J=1,NHT1
DO J=(NHT1+1),NHTT
   DO I=1,NHX1
      WRITE(20,110) XH(I,J),YH(I,J),ZH(I,J)
   END DO !I=1,NHX1
END DO !J=(NHT1+1),NHTT
DO K=2,NB
   PHI=DFLOAT(K-1)/DFLOAT(NB)*2.D0*PI
   WRITE(20,100) ' ZONE T="HUB',K,'" F=POINT, I=',NHX1,' J=',NHTT
   DO J=1,NHT1
      DO I=1,NHX1
         WRITE(20,110) XH(I,J),RH(I,J)*DCOS(TH(I,J)+PHI),RH(I,J)*DSIN(TH(I,J)+PHI)
      END DO !I=1,NHX1
   END DO !J=1,NHT1
   DO J=(NHT1+1),NHTT
      DO I=1,NHX1
         WRITE(20,110) XH(I,J),RH(I,J)*DCOS(TH(I,J)+PHI),RH(I,J)*DSIN(TH(I,J)+PHI)
      END DO !I=1,NHX1
   END DO !J=(NHT1+1),NHTT
END DO !K=2,NB
!-----------------------------------------------------------------------------------------------!
!    Formats                                                                                    !
!-----------------------------------------------------------------------------------------------!
100 FORMAT(A,I4,A,I4,A,I4)
110 FORMAT(3(2X,E23.16))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE HUBGRID
!-----------------------------------------------------------------------------------------------
