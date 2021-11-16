!-----------------------------------------------------------------------------------------------!
!    Generate nozzle grid                                                                       !
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
SUBROUTINE NOZZLEGRID
!-----------------------------------------------------------------------------------------------!
!    Created by: 01042011, J. Baltazar, version 1.0                                             !
!    Modified  : 21102013, J. Baltazar, version 1.0                                             !
!    Modified  : 23102014, J. Baltazar, version v2.1                                            !
!-----------------------------------------------------------------------------------------------!
!    Declarations                                                                               !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J
INTEGER :: LIDENTN
DOUBLE PRECISION :: XNO,XN1,XN2,XN3,XNP,TNO,TN1,TN2,TN3,TNP,PHI,S0,S1,PTN2
DOUBLE PRECISION :: ETA,XKS,ETAM,XKSM
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: S
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: XXG,YYG
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: XXW,YYW
!-----------------------------------------------------------------------------------------------!
I=1
DO WHILE (IDENTN(I:I).NE.' ')
   I=I+1
END DO !(IDENTN(I:I).NE.' ')
LIDENTN=I-1
!-----------------------------------------------------------------------------------------------!
!    Define Coordinates of Leading and Trailing Edge of Tip Blade                               !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
!! IF (ISTRIP > 0) NRP1=NRP1+1
   XN1=XP(NC1,NRP1)
   XN2=0.5D0*(XP(1     ,NRP1)+XP(NCP1    ,NRP1))
   XNP=0.5D0*(XP(NC/2+1,NRP1)+XP(NC/2+NC1,NRP1))
   TNP=0.5D0*(TP(NC/2+1,NRP1)+TP(NC/2+NC1,NRP1))
   PNN=2.D0*PNN !non-dimensional by propeller radius
!! IF (IGRIDI == 0) THEN
!!    TN1=TNP+(XN1-XNP)*2.D0*PI/PTN
!!    TN2=TNP+(XN2-XNP)*2.D0*PI/PTN
!! ELSEIF (IGRIDI == 1) THEN
      TN1=TP(NC1,NRP1)
      TN2=0.5D0*(TP(1,NRP1)+TP(NCP1,NRP1))
!! END IF !(IGRIDI)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Define Coordinates of End-Points of Boundary 1 of Domain                                   !
!-----------------------------------------------------------------------------------------------!
XNO=-LD
XN3= LD
IF (IP == 1) THEN
   IF (INTERN == 0) THEN
      CALL LININT(NNP,XNN,PNN,1,XNO,PTN2)
   ELSEIF (INTERN == 1) THEN
      CALL INTK1 (NNP,XNN,PNN,1,XNO,PTN2)
   ELSEIF (INTERN == 2) THEN
      CALL SPLINT(NNP,XNN,PNN,1,XNO,PTN2)
   END IF !(INTERN)
   TNO=TN1+(XNO-XN1)*2.D0*PI/PTN2
   IF (INTERN == 0) THEN
      CALL LININT(NNP,XNN,PNN,1,XN3,PTN2)
   ELSEIF (INTERN == 1) THEN
      CALL INTK1 (NNP,XNN,PNN,1,XN3,PTN2)
   ELSEIF (INTERN == 2) THEN
      CALL SPLINT(NNP,XNN,PNN,1,XN3,PTN2)
   END IF !(INTERN)
   TN3=TN2+(XN3-XN2)*2.D0*PI/PTN2
   IF (ISTRIP == 1) TN3=TPW(NND+1,NRW1)
ELSEIF (IP == 0) THEN
   IF (PTN < TOL) THEN
      TNO=0.D0
      TN3=0.D0
   ELSE !(PTN < TOL)
      TNO=0.D0
      IF (INTERN == 0) THEN
         CALL LININT(NNP,XNN,PNN,1,XN3,PTN2)
      ELSEIF (INTERN == 1) THEN
         CALL INTK1 (NNP,XNN,PNN,1,XN3,PTN2)
      ELSEIF (INTERN == 2) THEN
         CALL SPLINT(NNP,XNN,PNN,1,XN3,PTN2)
      END IF !(INTERN)
      TN3=TNO+(XN3-XNO)*2.D0*PI/PTN2
   END IF !(PTN < TOL)
END IF !(IP)
!-----------------------------------------------------------------------------------------------!
!    Define Point Distribution on Boundaries 1 and 3 of Domain - With Propeller                 !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
!-----------------------------------------------------------------------------------------------!
!    Panel Upstream Section of Nozzle                                                           !
!-----------------------------------------------------------------------------------------------!
   ALLOCATE(S(NNU1))
!! S0=DABS(XN1-XNO)/DABS(XN1-0.5D0*(XP(NC1+1,NRP1)+XP(NC1-1,NRP1)))/DFLOAT(NNU)
   S0=DABS(XN1-XNO)/DABS(0.5D0*(XN1-XNO)*(DCOS(PI/DFLOAT(NNU))-1.D0))/DFLOAT(NNU)
   S1=DABS(XN1-XNO)/DABS(0.5D0*(XP(NC1+2,NRP1)+XP(NC1-2,NRP1))- &
                         0.5D0*(XP(NC1+1,NRP1)+XP(NC1-1,NRP1)))/DFLOAT(NNU)
   CALL STRET2(S,S0,S1,NNU1)
   DO I=1,NNU1
!-----------------------------------------------------------------------------------------------!
!    Cosine Distribution on X                                                                   !
!-----------------------------------------------------------------------------------------------!
!!    PHI        =DFLOAT(I-1)*PI/DFLOAT(NNU)
!!    XN(I,1    )=0.5D0*(XNO+XN1)-0.5D0*(XN1-XNO)*DCOS(PHI)
!-----------------------------------------------------------------------------------------------!
!    Vinokur Distribution on X                                                                  !
!-----------------------------------------------------------------------------------------------!
      XN(I,1    )=XNO+(XN1-XNO)*S(I)
      XN(I,NNTT1)=XN(I,1)
      IF (INTERN == 0) THEN
         CALL LININT(NNP,XNN,PNN,1,XN(I,1),PTN2)
      ELSEIF (INTERN == 1) THEN
         CALL INTK1 (NNP,XNN,PNN,1,XN(I,1),PTN2)
      ELSEIF (INTERN == 2) THEN
         CALL SPLINT(NNP,XNN,PNN,1,XN(I,1),PTN2)
      END IF !(INTERN)
      TN(I,1    )=TNO+(XN(I,1)-XNO)*2.D0*PI/PTN2
      TN(I,NNTT1)=TN(I,1)+2.D0*PI/DFLOAT(NB)
   END DO !I=1,NNU1
   DEALLOCATE(S)
!-----------------------------------------------------------------------------------------------!
   IF (IGRIDI == 0) THEN
      DO I=1,NC
         XN(NNU1+I,1    )=0.5D0*(XP(NC1+I,NRP1)+XP(NC1-I,NRP1))
         IF (INTERN == 0) THEN
            CALL LININT(NNP,XNN,PNN,1,XN(NNU1+I,1),PTN2)
         ELSEIF (INTERN == 1) THEN
            CALL INTK1 (NNP,XNN,PNN,1,XN(NNU1+I,1),PTN2)
         ELSEIF (INTERN == 2) THEN
            CALL SPLINT(NNP,XNN,PNN,1,XN(NNU1+I,1),PTN2)
         END IF !(INTERN)
         TN(NNU1+I,1    )=TN1+(XN(NNU1+I,1)-XN1)*2.D0*PI/PTN2
         XN(NNU1+I,NNTT1)=0.5D0*(XP(NC1+I,NRP1)+XP(NC1-I,NRP1))
         TN(NNU1+I,NNTT1)=TN(NNU1+I,1)+2.D0*PI/DFLOAT(NB) 
      END DO !I=1,NC
   ELSEIF (IGRIDI == 1) THEN
      DO I=1,NC
         XN(NNU1+I,1    )=XP(NC1+I,NRP1)
         TN(NNU1+I,1    )=TP(NC1+I,NRP1)
         XN(NNU1+I,NNTT1)=XP(NC1-I,NRP1)
         TN(NNU1+I,NNTT1)=TP(NC1-I,NRP1)+2.D0*PI/DFLOAT(NB) 
      END DO !I=1,NC
   END IF !(IGRIDI)
!-----------------------------------------------------------------------------------------------!
!    Panel Downstream Section of Nozzle                                                         !
!-----------------------------------------------------------------------------------------------!
   ALLOCATE(S(NND+1))
!! S0=DABS(XN3-XN2)/DABS(XN2-0.5D0*(XP(2,NRP1)+XP(NCP1-1,NRP1)))/DFLOAT(NND)
   S0=DABS(XN3-XN2)/DABS(0.5D0*(XN3-XN2)*(DCOS(DFLOAT(NND-1)/DFLOAT(NND)*PI)+1.D0))/DFLOAT(NND)
   S1=DABS(XN3-XN2)/DABS(0.5D0*(XP(3,NRP1)+XP(NCP1-2,NRP1))- &
                         0.5D0*(XP(2,NRP1)+XP(NCP1-1,NRP1)))/DFLOAT(NND)
   CALL STRET2(S,S1,S0,NND+1)
   DO I=1,NND
!-----------------------------------------------------------------------------------------------!
!    Gap Model                                                                                  !
!-----------------------------------------------------------------------------------------------!
      IF (ISTRIP == 1) THEN
         XN(NN2+I,1    )=XPW(I+1,NRW1)
         XN(NN2+I,NNTT1)=XN (NN2+I,1)
         TN(NN2+I,1    )=TPW(I+1,NRW1)
         TN(NN2+I,NNTT1)=TN (NN2+I,1)+2.D0*PI/DFLOAT(NB)
!-----------------------------------------------------------------------------------------------!
!    Vinokur Distribution on X                                                                  !
!-----------------------------------------------------------------------------------------------!
      ELSE !(ISTRIP == 1)
         XN(NN2+I,1    )=XN2+(XN3-XN2)*S(I+1)
         XN(NN2+I,NNTT1)=XN(NN2+I,1)
         TN(NN2+I,1    )=TN2+(XN(NN2+I,1)-XN2)*2.D0*PI/PTN
         TN(NN2+I,NNTT1)=TN(NN2+I,1)+2.D0*PI/DFLOAT(NB)
!-----------------------------------------------------------------------------------------------!
!    Cosine Distribution on X                                                                   !
!-----------------------------------------------------------------------------------------------!
!!       PHI            =DFLOAT(I)*PI/DFLOAT(NND)
!!       XN(NN2+I,1    )=0.5D0*(XN2+XN3)-0.5D0*(XN3-XN2)*DCOS(PHI)
!-----------------------------------------------------------------------------------------------!
      END IF !(ITSRIP == 1)
   END DO !I=1,NND
   DEALLOCATE(S)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Define Point Distribution on Boundaries 1 and 3 of Domain - Without Propeller              !
!-----------------------------------------------------------------------------------------------!
IF (IP == 0) THEN
   DO I=1,NNX1
!-----------------------------------------------------------------------------------------------!
!    Cosine Distribution on X                                                                   !
!-----------------------------------------------------------------------------------------------!
      PHI=DFLOAT(I-1)*PI/DFLOAT(NNX)
      XN(I,1    )=0.5D0*(XNO+XN3)-0.5D0*(XN3-XNO)*DCOS(PHI)
      XN(I,NNTT1)=XN(I,1)
      IF (PTN < TOL) THEN
         TN(I,1)=TNO
      ELSE !(PTN < TOL)
         IF (INTERN == 0) THEN
            CALL LININT(NNP,XNN,PNN,1,XN(I,1),PTN2)
         ELSEIF (INTERN == 1) THEN
            CALL INTK1 (NNP,XNN,PNN,1,XN(I,1),PTN2)
         ELSEIF (INTERN == 2) THEN
            CALL SPLINT(NNP,XNN,PNN,1,XN(I,1),PTN2)
         END IF !(INTERN)
         TN(I,1)=TNO+(XN(I,1)-XNO)*2.D0*PI/PTN2
      END IF !(PTN < TOL)
      TN(I,NNTT1)=TN(I,1)+2.D0*PI/DFLOAT(NB)
   END DO !I=1,NNX1
END IF !(IP == 0)
!-----------------------------------------------------------------------------------------------!
!    Define Point Distribution on Boundaries 2 and 4 of Domain                                  !
!-----------------------------------------------------------------------------------------------!
DO J=1,NNTT1
   XN(1,J)=XNO
!-----------------------------------------------------------------------------------------------!
!    Cosine Distribution on Teta                                                                !
!-----------------------------------------------------------------------------------------------!
!! ALPHAN=ALPHAN*PI/180.D0
!! PHI       =ALPHAN+DFLOAT(J-1)/DFLOAT(NNTT1-1)*(PI-2.D0*ALPHAN)
!! TN(1   ,J)=(TNO+PI/DFLOAT(NB))-(PI/DFLOAT(NB))*DCOS(PHI)/DCOS(ALPHAN)
!! XN(NNX1,J)=XN3
!! TN(NNX1,J)=(TN3+PI/DFLOAT(NB))-(PI/DFLOAT(NB))*DCOS(PHI)/DCOS(ALPHAN)
!-----------------------------------------------------------------------------------------------!
!    Equidistant Distribution on Teta                                                           !
!-----------------------------------------------------------------------------------------------!
   PHI       =1.D0-2.D0*DFLOAT(J-1)/DFLOAT(NNTT1-1)
   TN(1   ,J)=(TNO+PI/DFLOAT(NB))-(PI/DFLOAT(NB))*PHI
   XN(NNX1,J)=XN3
   TN(NNX1,J)=(TN3+PI/DFLOAT(NB))-(PI/DFLOAT(NB))*PHI
END DO !J=1,NNTT1
!-----------------------------------------------------------------------------------------------!
!    Generate Grid on the Nozzle with the Program GRAPE                                         !
!-----------------------------------------------------------------------------------------------!
!    Define Arrays with Boundary Points                                                         !
!-----------------------------------------------------------------------------------------------!
ALLOCATE(XXG(NNX1,NNTT1),YYG(NNX1,NNTT1))
ALLOCATE(XXW(NNX1,NNTT1),YYW(NNX1,NNTT1))
XXG=0.D0
YYG=0.D0
XXW=0.D0
YYW=0.D0
!-----------------------------------------------------------------------------------------------!
DO I=1,NNX1
   XXG(I,1    )=XN(I,1    )
   YYG(I,1    )=TN(I,1    )
   XXG(I,NNTT1)=XN(I,NNTT1)
   YYG(I,NNTT1)=TN(I,NNTT1)
END DO !I=1,NNX1
!-----------------------------------------------------------------------------------------------!
DO J=1,NNTT1
   XXG(1   ,J)=XN(1   ,J)
   YYG(1   ,J)=TN(1   ,J)
   XXG(NNX1,J)=XN(NNX1,J)
   YYG(NNX1,J)=TN(NNX1,J)
END DO !J=1,NNTT1
!-----------------------------------------------------------------------------------------------!
!    Transfinite Interpolation                                                                  !
!-----------------------------------------------------------------------------------------------!
DO J=2,(NNTT1-1)
   ETAM=DFLOAT(J-1)/DFLOAT(NNTT1-1)
   ETA =1.D0-ETAM
   DO I=2,(NNX1-1)
      XKSM=DFLOAT(I-1)/DFLOAT(NNX1-1)
      XKS =1.D0-XKSM
      XXG(I,J)=XKS*XXG(1,J)+XKSM*XXG(NNX1,J)+ETA*XXG(I,1)+ETAM*XXG(I,NNTT1) &
              -XKS*ETA *XXG(1,1    )-XKSM*ETA *XXG(NNX1,1    ) &
              -XKS*ETAM*XXG(1,NNTT1)-XKSM*ETAM*XXG(NNX1,NNTT1)
      YYG(I,J)=XKS*YYG(1,J)+XKSM*YYG(NNX1,J)+ETA*YYG(I,1)+ETAM*YYG(I,NNTT1) &
              -XKS*ETA *YYG(1,1    )-XKSM*ETA *YYG(NNX1,1    ) &
              -XKS*ETAM*YYG(1,NNTT1)-XKSM*ETAM*YYG(NNX1,NNTT1)
   END DO !I=2,(NNX1-1)
END DO !J=2,(NNTT1-1)
!-----------------------------------------------------------------------------------------------!
!    Local Grid Correction                                                                      !
!-----------------------------------------------------------------------------------------------!
XXW=XXG
YYW=YYG
DO I=NNU1-5,NNU1
   DO J=2,NNT1
      XXG(I,J)=XXG(I,J)-0.16D0*DFLOAT(J-NNT1)/DFLOAT(2-NNT1)*(XXW(I,J)-XXW(NNU1-5,J))
      IF (INTERN == 0) THEN
         CALL LININT(NNP,XNN,PNN,1,XXG(I,J),PTN2)
      ELSEIF (INTERN == 1) THEN
         CALL INTK1 (NNP,XNN,PNN,1,XXG(I,J),PTN2)
      ELSEIF (INTERN == 2) THEN
         CALL SPLINT(NNP,XNN,PNN,1,XXG(I,J),PTN2)
      END IF !(INTERN)
      YYG(I,J)=YYG(NNU1-5,J)+(XXG(I,J)-XXG(NNU1-5,J))*2.D0*PI/PTN2
   END DO !J=2,NNT1
END DO !I=NNU1-5,NNU1
DO I=NNU1+1,NNU1+1
   DO J=2,NNT1
      XXG(I,J)=XXG(I,J)-0.08D0*DFLOAT(J-NNT1)/DFLOAT(2-NNT1)*(XXW(I,J)-XXW(NNU1-5,J))
      IF (INTERN == 0) THEN
         CALL LININT(NNP,XNN,PNN,1,XXG(I,J),PTN2)
      ELSEIF (INTERN == 1) THEN
         CALL INTK1 (NNP,XNN,PNN,1,XXG(I,J),PTN2)
      ELSEIF (INTERN == 2) THEN
         CALL SPLINT(NNP,XNN,PNN,1,XXG(I,J),PTN2)
      END IF !(INTERN)
      YYG(I,J)=YYG(NNU1,J)+(XXG(I,J)-XXG(NNU1,J))*2.D0*PI &
              /(PTN2*0.4D0+PTN2*0.6D0*DFLOAT(J-2)/DFLOAT(NNT1-2))
   END DO !J=2,NNT1
END DO !I=NNU1+1,NNU1+1
!-----------------------------------------------------------------------------------------------!
!    Reorder the nozzle grid coordinates: j=1       :   line betwwen blades                     !
!                                         j= NNT+1  :   pressure side                           !
!                                         j= NNT+2  :   suction side                            !
!                                         j=2(NNT+1):   line between blades (j=1)               !
!-----------------------------------------------------------------------------------------------!
!    Inner Side of the Nozzle                                                                   !
!-----------------------------------------------------------------------------------------------!
DO I=1,NNX1
!-----------------------------------------------------------------------------------------------!
!    Fill in First Half-Sector                                                                  !
!-----------------------------------------------------------------------------------------------!
   DO J=1,NNT1
      XN(I,J)=XXG(I,J+NNT1-1)
      TN(I,J)=YYG(I,J+NNT1-1)-2.D0*PI/DFLOAT(NB)
      CALL NOZZLEDEF('INNER',XN(I,J),RN(I,J))
   END DO !J=1,NNT1
!-----------------------------------------------------------------------------------------------!
!    Fill in Second Half-Sector                                                                 !
!-----------------------------------------------------------------------------------------------!
   DO J=1,NNT1
      XN(I,NNT1+J)=XXG(I,J)
      TN(I,NNT1+J)=YYG(I,J)
      CALL NOZZLEDEF('INNER',XN(I,NNT1+J),RN(I,NNT1+J))
   END DO !J=1,NNT1
!-----------------------------------------------------------------------------------------------!
!    Compute Cartesian Coordinates Complete Nozzle                                              !
!-----------------------------------------------------------------------------------------------!
   DO J=1,NNTT
      YN(I,J)=RN(I,J)*DCOS(TN(I,J)) 
      ZN(I,J)=RN(I,J)*DSIN(TN(I,J))
   END DO !J=1,NNTT
END DO !I=1,NNX1
!-----------------------------------------------------------------------------------------------!
!    Outer Side of the Nozzle                                                                   !
!-----------------------------------------------------------------------------------------------!
IF (IGRIDO == 0) THEN
!-----------------------------------------------------------------------------------------------!
   DO I=1,NNX
      XN(NNX1+I,1   )=XN(NNX1-I,1   )
      TN(NNX1+I,1   )=TN(NNX1-I,1   )
      XN(NNX1+I,NNTT)=XN(NNX1-I,NNTT)
      TN(NNX1+I,NNTT)=TN(NNX1-I,NNTT)
   END DO !I=1,NNX
!-----------------------------------------------------------------------------------------------!
   DO J=1,NNTT
      XN(NNXT,J)=XN(1,J)
      TN(NNXT,J)=TN(1,J)
   END DO !J=1,NNTT
!-----------------------------------------------------------------------------------------------!
!    Transfinite Interpolation                                                                  !
!-----------------------------------------------------------------------------------------------!
   DO J=2,(NNTT-1)
      ETAM=DFLOAT(J-1)/DFLOAT(NNTT-1)
      ETA =1.D0-ETAM
      DO I=(NNX1+1),(NNXT-1)
         XKSM=DFLOAT(I-NNX1)/DFLOAT(NNXT-NNX1)
         XKS =1.D0-XKSM
         XN(I,J)=XKS*XN(NNX1,J)+XKSM*XN(NNXT,J)+ETA*XN(I,1)+ETAM*XN(I,NNTT) &
                -XKS*ETA *XN(NNX1,1   )-XKSM*ETA *XN(NNXT,1   ) &
                -XKS*ETAM*XN(NNX1,NNTT)-XKSM*ETAM*XN(NNXT,NNTT)
         TN(I,J)=XKS*TN(NNX1,J)+XKSM*TN(NNXT,J)+ETA*TN(I,1)+ETAM*TN(I,NNTT) &
                -XKS*ETA *TN(NNX1,1   )-XKSM*ETA *TN(NNXT,1   ) &
                -XKS*ETAM*TN(NNX1,NNTT)-XKSM*ETAM*TN(NNXT,NNTT)
      END DO !I=(NNX1+1),(NNXT-1)
   END DO !J=2,(NNTT-1)
!-----------------------------------------------------------------------------------------------!
   DO I=1,NNX
!-----------------------------------------------------------------------------------------------!
!    Fill in First Half-Sector                                                                  !
!-----------------------------------------------------------------------------------------------!
      DO J=1,NNT1
         CALL NOZZLEDEF('OUTER',XN(NNX1+I,J),RN(NNX1+I,J))
      END DO !J=1,NNT1
!-----------------------------------------------------------------------------------------------!
!    Fill in Second Half-Sector                                                                 !
!-----------------------------------------------------------------------------------------------!
      DO J=1,NNT1
         CALL NOZZLEDEF('OUTER',XN(NNX1+I,NNT1+J),RN(NNX1+I,NNT1+J))
      END DO !J=1,NNT1
!-----------------------------------------------------------------------------------------------!
!    Compute Cartesian Coordinates Complete Nozzle                                              !
!-----------------------------------------------------------------------------------------------!
      DO J=1,NNTT
         YN(NNX1+I,J)=RN(NNX1+I,J)*DCOS(TN(NNX1+I,J)) 
         ZN(NNX1+I,J)=RN(NNX1+I,J)*DSIN(TN(NNX1+I,J))
      END DO !J=1,NNTT
   END DO !I=1,NNX
!-----------------------------------------------------------------------------------------------!
ELSEIF (IGRIDO == 1) THEN
!-----------------------------------------------------------------------------------------------!
   DO I=1,NNX
!-----------------------------------------------------------------------------------------------!
!    Fill in First Half-Sector                                                                  !
!-----------------------------------------------------------------------------------------------!
      DO J=1,NNT1
         XN(NNX1+I,J)=XN(NNX1-I,J)
         TN(NNX1+I,J)=TN(NNX1-I,J)
         CALL NOZZLEDEF('OUTER',XN(NNX1+I,J),RN(NNX1+I,J))
      END DO !J=1,NNT1
!-----------------------------------------------------------------------------------------------!
!    Fill in Second Half-Sector                                                                 !
!-----------------------------------------------------------------------------------------------!
      DO J=1,NNT1
         XN(NNX1+I,NNT1+J)=XN(NNX1-I,NNT1+J)
         TN(NNX1+I,NNT1+J)=TN(NNX1-I,NNT1+J)
         CALL NOZZLEDEF('OUTER',XN(NNX1+I,NNT1+J),RN(NNX1+I,NNT1+J))
      END DO !J=1,NNT1
!-----------------------------------------------------------------------------------------------!
!    Compute Cartesian Coordinates Complete Nozzle                                              !
!-----------------------------------------------------------------------------------------------!
      DO J=1,NNTT
         YN(NNX1+I,J)=RN(NNX1+I,J)*DCOS(TN(NNX1+I,J)) 
         ZN(NNX1+I,J)=RN(NNX1+I,J)*DSIN(TN(NNX1+I,J))
      END DO !J=1,NNTT
   END DO !I=1,NNX
END IF !(IGRIDO)
!-----------------------------------------------------------------------------------------------!
!    Deallocate Variables                                                                       !
!-----------------------------------------------------------------------------------------------!
DEALLOCATE(XXG,YYG)
DEALLOCATE(XXW,YYW)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE NOZZLEGRID
!-----------------------------------------------------------------------------------------------!
