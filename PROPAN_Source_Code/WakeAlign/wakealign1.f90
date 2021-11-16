!-----------------------------------------------------------------------------------------------!
!    Wake alignment model 1                                                                     !
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
SUBROUTINE WAKEALIGN1(JJ,TT)
!-----------------------------------------------------------------------------------------------!
!    Created by: 29052014, J. Baltazar, version 2.0, Wake Alignment Model 1                     !
!    Modified  : 05022016, J. Baltazar, 2016 version 1.0                                        !
!    Modified  : 09052016, J. Baltazar, 2016 version 1.2                                        !
!    Modified  : 17042019, J. Baltazar, 2019 version 1.0, Wake alignment robustness             !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,IP1,J,JM1,JP1,JJ,TT
DOUBLE PRECISION :: XPW_TMP(NPW1,NRW1),YPW_TMP(NPW1,NRW1),ZPW_TMP(NPW1,NRW1)
DOUBLE PRECISION :: RPW_TMP(NPW1,NRW1),TPW_TMP(NPW1,NRW1)
DOUBLE PRECISION ::  XN_TMP(NNXT,NNTT), YN_TMP(NNXT,NNTT), ZN_TMP(NNXT,NNTT)
DOUBLE PRECISION ::  RN_TMP(NNXT,NNTT), TN_TMP(NNXT,NNTT)
DOUBLE PRECISION :: XNW_TMP(NNW1,NNTT),YNW_TMP(NNW1,NNTT),ZNW_TMP(NNW1,NNTT)
DOUBLE PRECISION :: RNW_TMP(NNW1,NNTT),TNW_TMP(NNW1,NNTT)
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: XXNEW,YYNEW,ZZNEW,RRNEW,TANEW,S,STMP
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: XX,YY,ZZ,RR,TA,DX,DR,DT,BETAI
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: UX,UY,UZ,UR,UT,UXNEW,UYNEW,UZNEW
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: XTMP,YTMP,ZTMP,UTMP,VTMP,WTMP
!-----------------------------------------------------------------------------------------------!
!    Allocate Variables                                                                         !
!-----------------------------------------------------------------------------------------------!
ALLOCATE(XXNEW(NRW1),YYNEW(NRW1),ZZNEW(NRW1),RRNEW(NRW1),TANEW(NRW1),S(NRW1),STMP(NRV))
ALLOCATE(XX(NRW1),YY(NRW1),ZZ(NRW1),RR(NRW1),TA(NRW1),DX(NRW1),DR(NRW1),DT(NRW1),BETAI(NRW1))
ALLOCATE(XTMP(NRV,1),YTMP(NRV,1),ZTMP(NRV,1),UTMP(NRV,1),VTMP(NRV,1),WTMP(NRV,1))
ALLOCATE(UX(NRW1),UY(NRW1),UZ(NRW1),UR(NRW1),UT(NRW1),UXNEW(NRW1),UYNEW(NRW1),UZNEW(NRW1))
!-----------------------------------------------------------------------------------------------!
!    Wake Alignment                                                                             !
!-----------------------------------------------------------------------------------------------!
DO I=1,NWALIGN
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      XPW_TMP=XPW
      YPW_TMP=YPW
      ZPW_TMP=ZPW
      RPW_TMP=RPW
      TPW_TMP=TPW
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      XN_TMP=XN
      YN_TMP=YN
      ZN_TMP=ZN
      RN_TMP=RN
      TN_TMP=TN
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      XNW_TMP=XNW
      YNW_TMP=YNW
      ZNW_TMP=ZNW
      RNW_TMP=RNW
      TNW_TMP=TNW
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Prediction                                                                                 !
!-----------------------------------------------------------------------------------------------!
   XX=0.D0
   YY=0.D0
   ZZ=0.D0
   RR=0.D0
   TA=0.D0
!-----------------------------------------------------------------------------------------------!
   XX(:)=XPW(I,:)
   YY(:)=YPW(I,:)
   ZZ(:)=ZPW(I,:)
   RR(:)=RPW(I,:)
   TA(:)=TPW(I,:)
!-----------------------------------------------------------------------------------------------!
   DX=0.D0
   DR=0.D0
   DT=0.D0
!-----------------------------------------------------------------------------------------------!
   IP1=I+1
   DX(:)=XPW(IP1,:)-XPW(I,:)
   DR(:)=RPW(IP1,:)-RPW(I,:)
   DT(:)=TPW(IP1,:)-TPW(I,:)
!-----------------------------------------------------------------------------------------------!
!    Velocity Calculation Points                                                                !
!-----------------------------------------------------------------------------------------------!
   DO J=1,NRV
      JP1=JV(J)+1
      XTMP(J,1)=0.5D0*(XX(JV(J))+XX(JP1))
      YTMP(J,1)=0.5D0*(YY(JV(J))+YY(JP1))
      ZTMP(J,1)=0.5D0*(ZZ(JV(J))+ZZ(JP1))
   END DO !J=1,NRV
!-----------------------------------------------------------------------------------------------!
!    Velocity Calculation on Each Vertex                                                        !
!-----------------------------------------------------------------------------------------------!
   CALL VELFSTD('ROTATING',JJ,TT,NRV,1,XTMP,YTMP,ZTMP,UTMP,VTMP,WTMP)
!-----------------------------------------------------------------------------------------------!
!    Interpolation Variables                                                                    !
!-----------------------------------------------------------------------------------------------!
   S(1)=0.D0
   DO J=2,NRW1
      JM1=J-1
      S(J)=S(JM1)+DSQRT((XX(JM1)-XX(J))**2+(YY(JM1)-YY(J))**2+(ZZ(JM1)-ZZ(J))**2)
   END DO !J=2,NRW1
   S(:)=S(:)/S(NRW1)
!-----------------------------------------------------------------------------------------------!
   DO J=1,NRV
      JP1=JV(J)+1
      STMP(J)=0.5D0*S(JV(J))+0.5D0*S(JP1)
   END DO !J=1,NRV
!-----------------------------------------------------------------------------------------------!
!    Velocity Interpolation                                                                     !
!-----------------------------------------------------------------------------------------------!
   IF (INTERPW == 0) THEN
      CALL LININT(NRV,STMP,UTMP(:,1),NRW1,S,UX)
      CALL LININT(NRV,STMP,VTMP(:,1),NRW1,S,UY)
      CALL LININT(NRV,STMP,WTMP(:,1),NRW1,S,UZ)
   ELSEIF (INTERPW == 1) THEN
      CALL INTK1 (NRV,STMP,UTMP(:,1),NRW1,S,UX)
      CALL INTK1 (NRV,STMP,VTMP(:,1),NRW1,S,UY)
      CALL INTK1 (NRV,STMP,WTMP(:,1),NRW1,S,UZ)
   ELSEIF (INTERPW == 2) THEN
      CALL SPLINT(NRV,STMP,UTMP(:,1),NRW1,S,UX)
      CALL SPLINT(NRV,STMP,VTMP(:,1),NRW1,S,UY)
      CALL SPLINT(NRV,STMP,WTMP(:,1),NRW1,S,UZ)
   END IF !(INTERPW)
!-----------------------------------------------------------------------------------------------!
!    Cylindrical Coordinates                                                                    !
!-----------------------------------------------------------------------------------------------!
   UR= UY*DCOS(TA)+UZ*DSIN(TA)
   UT=-UY*DSIN(TA)+UZ*DCOS(TA)
!-----------------------------------------------------------------------------------------------!
!    Induced Hydrodynamic Pitch Angle                                                           !
!-----------------------------------------------------------------------------------------------!
   BETAI=DATAN2(UX,UT)
!-----------------------------------------------------------------------------------------------!
!    Mangler Condition at T.E.                                                                  !
!-----------------------------------------------------------------------------------------------!
   IF (I <= ITE) THEN
      DT(:)=DABS(TP(NC1,:)-TP(1,:))
      BETAI(:)=DATAN2(DABS(XP(NC1,:)-XP(1,:)),RR(:)*DT(:))
   END IF !(I <= ITE)
!-----------------------------------------------------------------------------------------------!
!    Correction                                                                                 !
!-----------------------------------------------------------------------------------------------!
!! DX=DT*RR*DTAN(BETAI)
   DT=DX/RR/DTAN(BETAI)
!-----------------------------------------------------------------------------------------------!
!    Predicted XXNEW                                                                            !
!-----------------------------------------------------------------------------------------------!
   XXNEW=XX+DX
!-----------------------------------------------------------------------------------------------!
!    Predicted RRNEW                                                                            !
!-----------------------------------------------------------------------------------------------!
!! RRNEW=RR+DR
   DO J=1,NRW1
      IF (INTERPW == 0) THEN
         CALL LININT(NPW1,XPW_TMP(:,J),RPW_TMP(:,J),1,XXNEW(J),RRNEW(J))
      ELSEIF (INTERPW == 1) THEN
         CALL INTK1 (NPW1,XPW_TMP(:,J),RPW_TMP(:,J),1,XXNEW(J),RRNEW(J))
      ELSEIF (INTERPW == 2) THEN
         CALL SPLINT(NPW1,XPW_TMP(:,J),RPW_TMP(:,J),1,XXNEW(J),RRNEW(J))
      END IF !(INTERPW)
   END DO !J=1,NRW1
   IF (ISTRIP == 1) THEN
      IF (XXNEW(NRW1) <= LD) THEN
         CALL NOZZLEDEF('INNER',XXNEW(NRW1),RRNEW(NRW1))
      ELSE
         CALL NOZZLEDEF('INNER',LD,RRNEW(NRW1))
      END IF
   END IF !(ISTRIP == 1)
!-----------------------------------------------------------------------------------------------!
!    Predicted TANEW                                                                            !
!-----------------------------------------------------------------------------------------------!
   TANEW=TA+DT
!-----------------------------------------------------------------------------------------------!
!    Predicted YYNEW and ZZNEW                                                                  !
!-----------------------------------------------------------------------------------------------!
   YYNEW=RRNEW*DCOS(TANEW)
   ZZNEW=RRNEW*DSIN(TANEW)
!-----------------------------------------------------------------------------------------------!
!    Displacement                                                                               !
!-----------------------------------------------------------------------------------------------!
   CALL BLADEWAKEDISP (I,XXNEW,RRNEW,TANEW,XPW_TMP,RPW_TMP,TPW_TMP)
   CALL NOZZLEDISP    (I,XXNEW,RRNEW,TANEW,XPW_TMP,RPW_TMP,TPW_TMP, XN_TMP, RN_TMP, TN_TMP)
   CALL NOZZLEWAKEDISP(I,XXNEW,RRNEW,TANEW,XPW_TMP,RPW_TMP,TPW_TMP,XNW_TMP,RNW_TMP,TNW_TMP)
!-----------------------------------------------------------------------------------------------!
!    Velocity Calculation Points                                                                !
!-----------------------------------------------------------------------------------------------!
   DO J=1,NRV
      JP1=JV(J)+1
      XTMP(J,1)=0.5D0*(XXNEW(JV(J))+XXNEW(JP1))
      YTMP(J,1)=0.5D0*(YYNEW(JV(J))+YYNEW(JP1))
      ZTMP(J,1)=0.5D0*(ZZNEW(JV(J))+ZZNEW(JP1))
   END DO !J=1,NRV
!-----------------------------------------------------------------------------------------------!
!    Velocity Calculation on Each Vertex                                                        !
!-----------------------------------------------------------------------------------------------!
   CALL VELFSTD('ROTATING',JJ,TT,NRV,1,XTMP,YTMP,ZTMP,UTMP,VTMP,WTMP)
!-----------------------------------------------------------------------------------------------!
!    Interpolation Variable                                                                     !
!-----------------------------------------------------------------------------------------------!
   S(1)=0.D0
   DO J=2,NRW1
      JM1=J-1
      S(J)=S(JM1)+DSQRT((XXNEW(JM1)-XXNEW(J))**2+ &
                        (YYNEW(JM1)-YYNEW(J))**2+ &
                        (ZZNEW(JM1)-ZZNEW(J))**2)
   END DO !J=2,NRW1
   S(:)=S(:)/S(NRW1)
!-----------------------------------------------------------------------------------------------!
   DO J=1,NRV
      JP1=JV(J)+1
      STMP(J)=0.5D0*S(JV(J))+0.5D0*S(JP1)
   END DO !J=1,NRV
!-----------------------------------------------------------------------------------------------!
!    Velocity Interpolation                                                                     !
!-----------------------------------------------------------------------------------------------!
   IF (INTERPW == 0) THEN
      CALL LININT(NRV,STMP,UTMP(:,1),NRW1,S,UXNEW)
      CALL LININT(NRV,STMP,VTMP(:,1),NRW1,S,UYNEW)
      CALL LININT(NRV,STMP,WTMP(:,1),NRW1,S,UZNEW)
   ELSEIF (INTERPW == 1) THEN
      CALL INTK1 (NRV,STMP,UTMP(:,1),NRW1,S,UXNEW)
      CALL INTK1 (NRV,STMP,VTMP(:,1),NRW1,S,UYNEW)
      CALL INTK1 (NRV,STMP,WTMP(:,1),NRW1,S,UZNEW)
   ELSEIF (INTERPW == 2) THEN
      CALL SPLINT(NRV,STMP,UTMP(:,1),NRW1,S,UXNEW)
      CALL SPLINT(NRV,STMP,VTMP(:,1),NRW1,S,UYNEW)
      CALL SPLINT(NRV,STMP,WTMP(:,1),NRW1,S,UZNEW)
   END IF !(INTERPW)
!-----------------------------------------------------------------------------------------------!
!    Average of the Two Velocities                                                              !
!-----------------------------------------------------------------------------------------------!
   UX=0.5D0*(UX+UXNEW)
   UY=0.5D0*(UY+UYNEW)
   UZ=0.5D0*(UZ+UZNEW)
!-----------------------------------------------------------------------------------------------!
!    Cylindrical Coordinates                                                                    !
!-----------------------------------------------------------------------------------------------!
   UR= UY*DCOS(TA)+UZ*DSIN(TA)
   UT=-UY*DSIN(TA)+UZ*DCOS(TA)
!-----------------------------------------------------------------------------------------------!
!    Induced Hydrodynamic Pitch Angle                                                           !
!-----------------------------------------------------------------------------------------------!
   BETAI=DATAN2(UX,UT)
!-----------------------------------------------------------------------------------------------!
!    Mangler Condition at T.E.                                                                  !
!-----------------------------------------------------------------------------------------------!
   IF (I <= ITE) THEN
      DT(:)=DABS(TP(NC1,:)-TP(1,:))
      BETAI(:)=DATAN2(DABS(XP(NC1,:)-XP(1,:)),RR(:)*DT(:))
   END IF !(I <= ITE)
!-----------------------------------------------------------------------------------------------!
!    New Position                                                                               !
!-----------------------------------------------------------------------------------------------!
!! DX=DT*RR*DTAN(BETAI)
   DT=DX/RR/DTAN(BETAI)
!-----------------------------------------------------------------------------------------------!
!    Corrected XXNEW                                                                            !
!-----------------------------------------------------------------------------------------------!
   XXNEW=XX+DX
!-----------------------------------------------------------------------------------------------!
!    Corrected RRNEW                                                                            !
!-----------------------------------------------------------------------------------------------!
!! RRNEW=RR+DR
   DO J=1,NRW1
      IF (INTERPW == 0) THEN
         CALL LININT(NPW1,XPW_TMP(:,J),RPW_TMP(:,J),1,XXNEW(J),RRNEW(J))
      ELSEIF (INTERPW == 1) THEN
         CALL INTK1 (NPW1,XPW_TMP(:,J),RPW_TMP(:,J),1,XXNEW(J),RRNEW(J))
      ELSEIF (INTERPW == 2) THEN
         CALL SPLINT(NPW1,XPW_TMP(:,J),RPW_TMP(:,J),1,XXNEW(J),RRNEW(J))
      END IF !(INTERPW)
   END DO !J=1,NRW1
   IF (ISTRIP == 1) THEN
      IF (XXNEW(NRW1) <= LD) THEN
         CALL NOZZLEDEF('INNER',XXNEW(NRW1),RRNEW(NRW1))
      ELSE
         CALL NOZZLEDEF('INNER',LD,RRNEW(NRW1))
      END IF
   END IF !(ISTRIP == 1)
!-----------------------------------------------------------------------------------------------!
!    Corrected TANEW                                                                            !
!-----------------------------------------------------------------------------------------------!
   TANEW=TA+DT
!-----------------------------------------------------------------------------------------------!
!    Corrected YYNEW and ZZNEW                                                                  !
!-----------------------------------------------------------------------------------------------!
   YYNEW=RRNEW*DCOS(TANEW)
   ZZNEW=RRNEW*DSIN(TANEW)
!-----------------------------------------------------------------------------------------------!
!    Displacement                                                                               !
!-----------------------------------------------------------------------------------------------!
   CALL BLADEWAKEDISP (I,XXNEW,RRNEW,TANEW,XPW_TMP,RPW_TMP,TPW_TMP)
   CALL NOZZLEDISP    (I,XXNEW,RRNEW,TANEW,XPW_TMP,RPW_TMP,TPW_TMP, XN_TMP, RN_TMP, TN_TMP)
   CALL NOZZLEWAKEDISP(I,XXNEW,RRNEW,TANEW,XPW_TMP,RPW_TMP,TPW_TMP,XNW_TMP,RNW_TMP,TNW_TMP)
!-----------------------------------------------------------------------------------------------!
END DO !I=1,NWALIGN
!-----------------------------------------------------------------------------------------------!
!    Frozen Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
DO I=NWALIGN+1,NPW
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      XPW_TMP=XPW
      YPW_TMP=YPW
      ZPW_TMP=ZPW
      RPW_TMP=RPW
      TPW_TMP=TPW
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      XN_TMP=XN
      YN_TMP=YN
      ZN_TMP=ZN
      RN_TMP=RN
      TN_TMP=TN
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      XNW_TMP=XNW
      YNW_TMP=YNW
      ZNW_TMP=ZNW
      RNW_TMP=RNW
      TNW_TMP=TNW
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Prediction                                                                                 !
!-----------------------------------------------------------------------------------------------!
   XX=0.D0
   YY=0.D0
   ZZ=0.D0
   RR=0.D0
   TA=0.D0
!-----------------------------------------------------------------------------------------------!
   XX(:)=XPW(I,:)
   YY(:)=YPW(I,:)
   ZZ(:)=ZPW(I,:)
   RR(:)=RPW(I,:)
   TA(:)=TPW(I,:)
!-----------------------------------------------------------------------------------------------!
   DX=0.D0
   DR=0.D0
   DT=0.D0
!-----------------------------------------------------------------------------------------------!
   IP1=I+1
   DX(:)=XPW(IP1,:)-XPW(I,:)
   DR(:)=RPW(IP1,:)-RPW(I,:)
   DT(:)=TPW(IP1,:)-TPW(I,:)
!-----------------------------------------------------------------------------------------------!
!    New Position                                                                               !
!-----------------------------------------------------------------------------------------------!
!! DX=DT*RR*DTAN(BETAI)
   DT=DX/RR/DTAN(BETAI)
!-----------------------------------------------------------------------------------------------!
!    Frozen XXNEW                                                                               !
!-----------------------------------------------------------------------------------------------!
   XXNEW=XX+DX
!-----------------------------------------------------------------------------------------------!
!    Frozen RRNEW                                                                               !
!-----------------------------------------------------------------------------------------------!
!! RRNEW=RR+DR
   DO J=1,NRW1
      IF (INTERPW == 0) THEN
         CALL LININT(NPW1,XPW_TMP(:,J),RPW_TMP(:,J),1,XXNEW(J),RRNEW(J))
      ELSEIF (INTERPW == 1) THEN
         CALL INTK1 (NPW1,XPW_TMP(:,J),RPW_TMP(:,J),1,XXNEW(J),RRNEW(J))
      ELSEIF (INTERPW == 2) THEN
         CALL SPLINT(NPW1,XPW_TMP(:,J),RPW_TMP(:,J),1,XXNEW(J),RRNEW(J))
      END IF !(INTERPW)
   END DO !J=1,NRW1
   IF (ISTRIP == 1) THEN
      IF (XXNEW(NRW1) <= LD) THEN
         CALL NOZZLEDEF('INNER',XXNEW(NRW1),RRNEW(NRW1))
      ELSE
         CALL NOZZLEDEF('INNER',LD,RRNEW(NRW1))
      END IF
   END IF !(ISTRIP == 1)
!-----------------------------------------------------------------------------------------------!
!    Frozen TANEW                                                                               !
!-----------------------------------------------------------------------------------------------!
   TANEW=TA+DT
!-----------------------------------------------------------------------------------------------!
!    Frozen YYNEW and ZZNEW                                                                     !
!-----------------------------------------------------------------------------------------------!
   YYNEW=RRNEW*DCOS(TANEW)
   ZZNEW=RRNEW*DSIN(TANEW)
!-----------------------------------------------------------------------------------------------!
!    Displacement                                                                               !
!-----------------------------------------------------------------------------------------------!
   CALL BLADEWAKEDISP (I,XXNEW,RRNEW,TANEW,XPW_TMP,RPW_TMP,TPW_TMP)
   CALL NOZZLEDISP    (I,XXNEW,RRNEW,TANEW,XPW_TMP,RPW_TMP,TPW_TMP, XN_TMP, RN_TMP, TN_TMP)
   CALL NOZZLEWAKEDISP(I,XXNEW,RRNEW,TANEW,XPW_TMP,RPW_TMP,TPW_TMP,XNW_TMP,RNW_TMP,TNW_TMP)
END DO !I=NWALIGN+1,NPW
!-----------------------------------------------------------------------------------------------!
!    Deallocate Variables                                                                       !
!-----------------------------------------------------------------------------------------------!
DEALLOCATE(XXNEW,YYNEW,ZZNEW,RRNEW,TANEW,S,STMP)
DEALLOCATE(XX,YY,ZZ,RR,TA,DX,DR,DT,BETAI)
DEALLOCATE(XTMP,YTMP,ZTMP,UTMP,VTMP,WTMP)
DEALLOCATE(UX,UY,UZ,UR,UT,UXNEW,UYNEW,UZNEW)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE WAKEALIGN1
!-----------------------------------------------------------------------------------------------!
