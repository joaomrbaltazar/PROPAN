!-----------------------------------------------------------------------------------------------!
!    Generate Nozzle Grid                                                                       !
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
SUBROUTINE NOZZLEGRID
!-----------------------------------------------------------------------------------------------!
!    Created by: 01042011, J. Baltazar, version 1.0                                             !
!    Modified  : 21102013, J. Baltazar, version 1.0                                             !
!    Modified  : 02122014, J. Baltazar, 2014 version 1.2                                        !
!    Modified  : 19012015, J. Baltazar, 2015 version 1.1, (ALPHAN) cosine on theta              !
!-----------------------------------------------------------------------------------------------!
!    Declarations                                                                               !
!-----------------------------------------------------------------------------------------------!
USE PROPANEL_MOD
IMPLICIT NONE
INTEGER :: I,J,K
INTEGER :: LIDENTN,NNU1,NN2
DOUBLE PRECISION :: XKS,XKSM,ETA,ETAM
DOUBLE PRECISION :: XN0,XN1,XN2,XN3,XNP,TN0,TN1,TN2,TN3,TNP,PHI,S0,S1
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: S
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: XXG,YYG
!-----------------------------------------------------------------------------------------------!
I=1
DO WHILE (IDENTN(I:I).NE.' ')
   I=I+1
END DO !(IDENTN(I:I).NE.' ')
LIDENTN=I-1
!-----------------------------------------------------------------------------------------------!
!    Counters                                                                                   !
!-----------------------------------------------------------------------------------------------!
NNX  =NNU+NC+NND
NNX1 =NNX+1
NNU1 =NNU+1
NN2  =NNU1+NC
NNT1 =NNT+1
NNTT =2*NNT1
NNTT1=NNTT-1
NNXT =NNX1+NNX
NCN  =NC
!-----------------------------------------------------------------------------------------------!
ALLOCATE(XN(NNXT,NNTT),YN(NNXT,NNTT),ZN(NNXT,NNTT),RN(NNXT,NNTT),TN(NNXT,NNTT))
XN=0.D0
YN=0.D0
ZN=0.D0
RN=0.D0
TN=0.D0
!-----------------------------------------------------------------------------------------------!
!    Define Coordinates of Leading and Trailing Edge of Tip Blade                               !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   IF (ISTRIP > 0) NRP1=NRP1+1
   XN1=XP(NC1,NRP1)
   XN2=0.5D0*(XP(1     ,NRP1)+XP(NCP1    ,NRP1))
   XNP=0.5D0*(XP(NC/2+1,NRP1)+XP(NC/2+NC1,NRP1))
   TNP=0.5D0*(TP(NC/2+1,NRP1)+TP(NC/2+NC1,NRP1))
   PTN=2.D0*PTN !non-dimensional by propeller radius
   IF (IGRIDI == 0) THEN
      TN1=TNP+(XN1-XNP)*2.D0*PI/PTN
      TN2=TNP+(XN2-XNP)*2.D0*PI/PTN
   ELSEIF (IGRIDI == 1) THEN
      TN1=TP(NC1,NRP1)
      TN2=0.5D0*(TP(1,NRP1)+TP(NCP1,NRP1))
   END IF !(IGRIDI)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Define Coordinates of End-Points of Boundary 1 of Domain                                   !
!-----------------------------------------------------------------------------------------------!
XN0=-LD
XN3= LD
IF (IP == 1) THEN
   TN0=TN1+(XN0-XN1)*2.D0*PI/PTN
   TN3=TN2+(XN3-XN2)*2.D0*PI/PTN
   IF (ISTRIP == 1) TN3=TPW(NND+1,NRW1)
ELSEIF (IP == 0) THEN
   IF (PTN < TOL) THEN
      TN0=0.D0
      TN3=0.D0
   ELSE !(PTN < TOL)
      TN0=0.D0
      TN3=TN0+(XN3-XN0)*2.D0*PI/PTN
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
!! S0=DABS(XN1-XN0)/DABS(XN1-0.5D0*(XP(NC1+1,NRP1)+XP(NC1-1,NRP1)))/DFLOAT(NNU)
   S0=DABS(XN1-XN0)/DABS(0.5D0*(XN1-XN0)*(DCOS(PI/DFLOAT(NNU))-1.D0))/DFLOAT(NNU)
   S1=DABS(XN1-XN0)/DABS(0.5D0*(XP(NC1+2,NRP1)+XP(NC1-2,NRP1))- &
                         0.5D0*(XP(NC1+1,NRP1)+XP(NC1-1,NRP1)))/DFLOAT(NNU)
   CALL STRET2(S,S0,S1,NNU1)
   DO I=1,NNU1
!-----------------------------------------------------------------------------------------------!
!    Cosine Distribution on X                                                                   !
!-----------------------------------------------------------------------------------------------!
!!    PHI        =DFLOAT(I-1)*PI/DFLOAT(NNU)
!!    XN(I,1    )=0.5D0*(XN0+XN1)-0.5D0*(XN1-XN0)*DCOS(PHI)
!-----------------------------------------------------------------------------------------------!
!    Vinokur Distribution on X                                                                  !
!-----------------------------------------------------------------------------------------------!
      XN(I,1    )=XN0+(XN1-XN0)*S(I)
      XN(I,NNTT1)=XN(I,1)
      TN(I,1    )=TN0+(XN(I,1)-XN0)*2.D0*PI/PTN
      TN(I,NNTT1)=TN(I,1)+2.D0*PI/DFLOAT(NB)
   END DO !I=1,NNU1
   DEALLOCATE(S)
!-----------------------------------------------------------------------------------------------!
   IF (IGRIDI == 0) THEN
      DO I=1,NC
         XN(NNU1+I,1    )=0.5D0*(XP(NC1+I,NRP1)+XP(NC1-I,NRP1))
         TN(NNU1+I,1    )=TN1+(XN(NNU1+I,1)-XN1)*2.D0*PI/PTN
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
      XN(I,1    )=0.5D0*(XN0+XN3)-0.5D0*(XN3-XN0)*DCOS(PHI)
      XN(I,NNTT1)=XN(I,1)
      IF (PTN < TOL) THEN
         TN(I,1)=TN0
      ELSE !(PTN < TOL)
         TN(I,1)=TN0+(XN(I,1)-XN0)*2.D0*PI/PTN
      END IF !(PTN < TOL)
      TN(I,NNTT1)=TN(I,1)+2.D0*PI/DFLOAT(NB)
   END DO !I=1,NNX1
END IF !(IP == 0)
!-----------------------------------------------------------------------------------------------!
!    Define Point Distribution on Boundaries 2 and 4 of Domain                                  !
!-----------------------------------------------------------------------------------------------!
ALPHAN=ALPHAN*PI/180.D0
DO J=1,NNTT1
   XN(1,J)=XN0
!-----------------------------------------------------------------------------------------------!
!    Cosine Distribution on Theta                                                               !
!-----------------------------------------------------------------------------------------------!
   IF (ALPHAN > TOL) THEN
      PHI       =ALPHAN+DFLOAT(J-1)/DFLOAT(NNTT1-1)*(PI-2.D0*ALPHAN)
      TN(1   ,J)=(TN0+PI/DFLOAT(NB))-(PI/DFLOAT(NB))*DCOS(PHI)/DCOS(ALPHAN)
      XN(NNX1,J)=XN3
      TN(NNX1,J)=(TN3+PI/DFLOAT(NB))-(PI/DFLOAT(NB))*DCOS(PHI)/DCOS(ALPHAN)
!-----------------------------------------------------------------------------------------------!
!    Equidistant Distribution on Theta                                                          !
!-----------------------------------------------------------------------------------------------!
   ELSE !(ALPHAN > TOL)
      PHI       =1.D0-2.D0*DFLOAT(J-1)/DFLOAT(NNTT1-1)
      TN(1   ,J)=(TN0+PI/DFLOAT(NB))-(PI/DFLOAT(NB))*PHI
      XN(NNX1,J)=XN3
      TN(NNX1,J)=(TN3+PI/DFLOAT(NB))-(PI/DFLOAT(NB))*PHI
   END IF !(ALPHAN > TOL)
END DO !J=1,NNTT1
!-----------------------------------------------------------------------------------------------!
!    Generate Grid on the Nozzle with the Program GRAPE                                         !
!-----------------------------------------------------------------------------------------------!
!    Define Arrays with Boundary Points                                                         !
!-----------------------------------------------------------------------------------------------!
ALLOCATE(XXG(NNX1,NNTT1),YYG(NNX1,NNTT1))
XXG=0.D0
YYG=0.D0
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
IF (IGRIDI == 0) THEN
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
ELSEIF (IGRIDI == 1) THEN
   CALL GRAPE(NNX1,NNTT1,XXG,YYG,5000)
END IF !(IGRIDI)
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
                -XKS*ETAM*XN(NNX1,NNTT)-ETAM*XKSM*XN(NNXT,NNTT)
         TN(I,J)=XKS*TN(NNX1,J)+XKSM*TN(NNXT,J)+ETA*TN(I,1)+ETAM*TN(I,NNTT) &
                -XKS*ETA *TN(NNX1,1   )-XKSM*ETA *TN(NNXT,1   ) &
                -XKS*ETAM*TN(NNX1,NNTT)-ETAM*XKSM*TN(NNXT,NNTT)
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
!    Write Nozzle Grid in Tecplot Format                                                        !
!-----------------------------------------------------------------------------------------------!
K=1
WRITE(20,100) ' ZONE T="NOZZLE',K,'" F=POINT, I=',NNXT,' J=',NNTT
DO J=1,NNTT
   DO I=1,NNXT
      WRITE(20,110) XN(I,J),YN(I,J),ZN(I,J)
   END DO !I=1,NNXT
END DO !J=1,NNTT
DO K=2,NB
   PHI=DFLOAT(K-1)/DFLOAT(NB)*2.D0*PI
   WRITE(20,100) ' ZONE T="NOZZLE',K,'" F=POINT, I=',NNXT,' J=',NNTT
   DO J=1,NNTT
      DO I=1,NNXT
         WRITE(20,110) XN(I,J),RN(I,J)*DCOS(TN(I,J)+PHI),RN(I,J)*DSIN(TN(I,J)+PHI)
      END DO !I=1,NNXT
   END DO !J=1,NNTT
END DO !K=2,NB
!-----------------------------------------------------------------------------------------------!
!    Deallocate Variables                                                                       !
!-----------------------------------------------------------------------------------------------!
DEALLOCATE(XXG,YYG)
!-----------------------------------------------------------------------------------------------!
!    Formats                                                                                    !
!-----------------------------------------------------------------------------------------------!
100 FORMAT(A,I4,A,I4,A,I4)
110 FORMAT(3(2X,E23.16))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE NOZZLEGRID
!-----------------------------------------------------------------------------------------------!
