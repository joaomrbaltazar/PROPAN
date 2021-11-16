!-----------------------------------------------------------------------------------------------!
!    Kinematic boundary condition for cavity flow                                               !
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
SUBROUTINE CAVTHICKS(JJ,TT,CC)
!-----------------------------------------------------------------------------------------------!
!    Created by: 07112014, J. Baltazar, Cavitation Model                                        !
!    Modified  : 17112014, J. Baltazar, version 3.1, Unsteady Cavitation Model                  !
!    Modified  : 25112014, J. Baltazar, version 3.2, Mid-Chord Cavitation                       !
!    Modified  : 27112014, J. Baltazar, version 3.3, Super-Cavitation Model                     !
!    Modified  : 09122014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!    Modified  : 06032015, J. Baltazar, Robustness for Kutta Condition                          !
!    Modified  : 29022016, J. Baltazar, Revision                                                !
!    Modified  : 22042016, J. Baltazar, 2016 version 1.1                                        !
!    Modified  : 31052016, J. Baltazar, 2016 version 1.2                                        !
!    Modified  : 08072016, J. Baltazar, 2016 version 1.3, new relaxation for unsteady terms     !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,J2,JJ,TT,CC,NCC
DOUBLE PRECISION :: VWX,VWY,VWZ
DOUBLE PRECISION :: T2X,T2Y,T2Z,T2,ET12
DOUBLE PRECISION :: VTT1S,VTT2S,VTTNP
DOUBLE PRECISION :: A,B,C,DETADT,DCAVDT,FRELAX
DOUBLE PRECISION :: AA1,BB1,AA2,BB2,CC2,SA1,SB1,SC1,SA,SB,SP,SS
DOUBLE PRECISION :: CAV1,CAV2,CAM1,CAM2
!-----------------------------------------------------------------------------------------------!
DO J=1,NRP
   IF (IDS(J,TT) /= 0) THEN
!-----------------------------------------------------------------------------------------------!
      CAV1=0.D0
      DO I=IDS(J,TT),IRS(J,TT)
         T2X=AT2XP(I,J)
         T2Y=AT2YP(I,J)
         T2Z=AT2ZP(I,J)
         T2 =DSQRT(T2X*T2X+T2Y*T2Y+T2Z*T2Z)
         T2X=T2X/T2
         T2Y=T2Y/T2
         T2Z=T2Z/T2
         ET12=ET1XP(I,J)*T2X+ET1YP(I,J)*T2Y+ET1ZP(I,J)*T2Z
!-----------------------------------------------------------------------------------------------!
         IF (TT == 0) THEN
            CALL VWAKE(TT,XP0(I,J),YP0(I,J),ZP0(I,J),1,0,VWX,VWY,VWZ)
         ELSE !(TT)
            CALL VWAKE(TT,XP0(I,J),YP0(I,J),ZP0(I,J),1,IFREQ,VWX,VWY,VWZ)
         END IF !(TT)
         DETADT=0.D0
         FRELAX=0.D0
         IF (TT > 0) THEN
            IF (THICKP(I,J,TT) == 0.D0) THEN
               DETADT=0.D0
            ELSEIF (TT == 1) THEN
               DETADT=THICKP(I,J,1)/DTETA-THICKP(I,J,0)/DTETA
            ELSE !(TT == 1)
               DETADT=THICKP(I,J,TT-2)*0.5D0/DTETA- &
                      THICKP(I,J,TT-1)*2.0D0/DTETA+ &
                      THICKP(I,J,TT  )*1.5D0/DTETA
            END IF !(TT == 1)
            IF (TT <= (2*NTETA)) THEN
               FRELAX=DFLOAT(TT-1)/DFLOAT(2*NTETA-1)
            ELSE !(TT <= (2*NTETA))
               FRELAX=1.D0
            END IF !(TT <= (2*NTETA))
         END IF !(TT > 0)
!-----------------------------------------------------------------------------------------------!
         IF (IROTOR == 0) THEN
            VTT2S=VWX*UU(JJ)/PI*T2X+ &                  !x1 component of Total Velocity
                 (VWY*UU(JJ)/PI-ZP0(I,J))*T2Y+ &        !y1 component of Total Velocity
                 (VWZ*UU(JJ)/PI+YP0(I,J))*T2Z+ &        !z1 component of Total Velocity
                 +VT2S(I,J)                             !s1 perturbation component
            VTTNP=VWX*UU(JJ)/PI*UNXP0(I,J)+ &           !x2 component of Total Velocity
                 (VWY*UU(JJ)/PI-ZP0(I,J))*UNYP0(I,J)+ & !y2 component of Total Velocity
                 (VWZ*UU(JJ)/PI+YP0(I,J))*UNZP0(I,J)+ & !z2 component of Total Velocity
                 -SOURCEP(I,J,1)                        !s2 perturbation component
         ELSEIF (IROTOR == 1) THEN
            VTT2S=VWX/UU(JJ)*T2X+ &                     !x1 component of Total Velocity
                 (VWY/UU(JJ)-ZP0(I,J))*T2Y+ &           !y1 component of Total Velocity
                 (VWZ/UU(JJ)+YP0(I,J))*T2Z+ &           !z1 component of Total Velocity
                 +VT2S(I,J)                             !s1 perturbation component
            VTTNP=VWX/UU(JJ)*UNXP0(I,J)+ &              !x2 component of Total Velocity
                 (VWY/UU(JJ)-ZP0(I,J))*UNYP0(I,J)+ &    !y2 component of Total Velocity
                 (VWZ/UU(JJ)+YP0(I,J))*UNZP0(I,J)+ &    !z2 component of Total Velocity
                 -SOURCEP(I,J,1)                        !s2 perturbation component
         END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
         A=VTT1P(I,J)-VTT2S*ET12
         B=VTT2S-VTT1P(I,J)*ET12
         C=1.D0-ET12*ET12
!-----------------------------------------------------------------------------------------------!
         SB =2.D0*DSQRT(AT1XP(I,J)**2+AT1YP(I,J)**2+AT1ZP(I,J)**2)
         AA1= 1.D0/SB
         BB1=-1.D0/SB
!-----------------------------------------------------------------------------------------------!
         IF (J == 1) THEN
            WORK1=AA1*A
            WORK2=C*(VTTNP-DETADT*FRELAX)-BB1*A*CAV1
            WORK3=0.D0
         ELSEIF (J == 2) THEN
            IF ((I >= IDS(J-1,TT)).AND.(I <= IRS(J-1,TT))) THEN
               SB1=DSQRT(AT2XP(I,J-1)**2+AT2YP(I,J-1)**2+AT2ZP(I,J-1)**2)
               SC1=DSQRT(AT2XP(I,J  )**2+AT2YP(I,J  )**2+AT2ZP(I,J  )**2)
               AA2= 1.D0/(SB1+SC1)
               BB2=-1.D0/(SB1+SC1)
               WORK1=AA1*A+AA2*B*0.5D0
               WORK2=C*(VTTNP-DETADT*FRELAX)-(BB1*A+AA2*B*0.5D0)*CAV1
               WORK3=BB2*B*THICKP(I,J-1,TT)
            ELSE
               WORK1=AA1*A
               WORK2=C*(VTTNP-DETADT*FRELAX)-BB1*A*CAV1
               WORK3=0.D0
            END IF
         ELSEIF ((I >= IDS(J-1,TT)).AND.(I <= IRS(J-1,TT)).AND. &
                                                 (I >= IDS(J-2,TT)).AND.(I <= IRS(J-2,TT))) THEN
            SA1=DSQRT(AT2XP(I,J-2)**2+AT2YP(I,J-2)**2+AT2ZP(I,J-2)**2)
            SB1=DSQRT(AT2XP(I,J-1)**2+AT2YP(I,J-1)**2+AT2ZP(I,J-1)**2)
            SC1=DSQRT(AT2XP(I,J  )**2+AT2YP(I,J  )**2+AT2ZP(I,J  )**2)
            SA =SA1+SB1
            SB =SB1+SC1
            SP =SA*SB
            SS =SA+SB
            AA2=(SA+2.D0*SB)/(SB*SS)
            BB2=-SS/SP
            CC2=SB/(SA*SS)
            WORK1=AA1*A+AA2*B*0.5D0
            WORK2=C*(VTTNP-DETADT*FRELAX)-(BB1*A+AA2*B*0.5D0)*CAV1
            WORK3=BB2*B*THICKP(I,J-1,TT)+CC2*B*THICKP(I,J-2,TT)
         ELSEIF ((I >= IDS(J-1,TT)).AND.(I <= IRS(J-1,TT))) THEN
            SB1=DSQRT(AT2XP(I,J-1)**2+AT2YP(I,J-1)**2+AT2ZP(I,J-1)**2)
            SC1=DSQRT(AT2XP(I,J  )**2+AT2YP(I,J  )**2+AT2ZP(I,J  )**2)
            AA2= 1.D0/(SB1+SC1)
            BB2=-1.D0/(SB1+SC1)
            WORK1=AA1*A+AA2*B*0.5D0
            WORK2=C*(VTTNP-DETADT*FRELAX)-(BB1*A+AA2*B*0.5D0)*CAV1
            WORK3=BB2*B*THICKP(I,J-1,TT)
         ELSE
            WORK1=AA1*A
            WORK2=C*(VTTNP-DETADT*FRELAX)-BB1*A*CAV1
            WORK3=0.D0
         END IF
!-----------------------------------------------------------------------------------------------!
!    Cavity Thickness                                                                           !
!-----------------------------------------------------------------------------------------------!
         CAV2=(WORK2-WORK3)/WORK1
         THICKP(I,J,TT)=0.5D0*(CAV1+CAV2)
         if ((i >= ids(j,tt)).and.(i <= ids(j,tt)+imcs).and.(thickp(i,j,tt) < 0.d0)) then !***Correction***
            write(30,'(A,I4,A,I4)') ' CC=',cc,', TT=',tt
            write(30,'(A,I4,I4)  ') ' Eta > 0 at the K.B.C., (I,J)=',i,j
            cav2=dabs(cav2)
            thickp(i,j,tt)=dabs(thickp(i,j,tt))
         end if
         CAV1=CAV2
      END DO !I=IDS(J,TT),IRS(J,TT)
!-----------------------------------------------------------------------------------------------!
!    Cavity Extent                                                                              !
!-----------------------------------------------------------------------------------------------!
      IF ((THICKP(IRS(J,TT),J,TT) > 0.D0).AND. &
                                  (THICKP(IRS(J,TT),J,TT)-THICKP(IRS(J,TT)-1,J,TT) < 0.D0)) THEN
         NCC=MIN(IRS(J,TT)+3+2*NCP/100,NCP)
      ELSE !((THICKP(IRS(J,TT),J,TT) > 0.D0).AND. &
!                                      (THICKP(IRS(J,TT),J,TT)-THICKP(IRS(J,TT)-1,J,TT) < 0.D0))
         NCC=MIN(IRS(J,TT)+2+2*NCP/100,NCP)
      END IF !((THICKP(IRS(J,TT),J,TT) > 0.D0).AND. &
!                                      (THICKP(IRS(J,TT),J,TT)-THICKP(IRS(J,TT)-1,J,TT) < 0.D0))
!-----------------------------------------------------------------------------------------------!
      IF (IRS(J,TT) < NCP) THEN
         DO I=IRS(J,TT)+1,NCP
            SA1=DSQRT(AT2XP(I-2,J)**2+AT2YP(I-2,J)**2+AT2ZP(I-2,J)**2)
            SB1=DSQRT(AT2XP(I-1,J)**2+AT2YP(I-1,J)**2+AT2ZP(I-1,J)**2)
            SC1=DSQRT(AT2XP(I  ,J)**2+AT2YP(I  ,J)**2+AT2ZP(I  ,J)**2)
            SA =SA1+SB1
            SB =SB1+SC1
            THICKP(I,J,TT)=THICKP(I-1,J,TT)+(THICKP(I-1,J,TT)-THICKP(I-2,J,TT))/SA*SB
            IF (THICKP(I,J,TT) < 0.D0) THICKP(I,J,TT)=0.D0
         END DO !I=IRS(J,TT)+1,NCP
      END IF !(IRS(J,TT) < NCP)
!-----------------------------------------------------------------------------------------------!
      IF (IDS(J,TT) > NC1) THEN
         DO I=IDS(J,TT)-1,NC1,-1
            SA1=DSQRT(AT1XP(I+2,J)**2+AT1YP(I+2,J)**2+AT1ZP(I+2,J)**2) !***
            SB1=DSQRT(AT1XP(I+1,J)**2+AT1YP(I+1,J)**2+AT1ZP(I+1,J)**2)
            SC1=DSQRT(AT1XP(I  ,J)**2+AT1YP(I  ,J)**2+AT1ZP(I  ,J)**2)
            SA =SA1+SB1
            SB =SB1+SC1
            THICKP(I,J,TT)=THICKP(I+1,J,TT)+(THICKP(I+1,J,TT)-THICKP(I+2,J,TT))/SA*SB
            IF (THICKP(I,J,TT) < 0.D0) THICKP(I,J,TT)=0.D0
         END DO !I=IDS(J,TT)-1,NC1,-1
      END IF !(IDS(J,TT) > NC1)
!-----------------------------------------------------------------------------------------------!
!    Mid-Chord Cavitation                                                                       !
!-----------------------------------------------------------------------------------------------!
      IF (IMCS > 0) THEN
         IF (CPNP(IDP(J,TT)-1,J,TT) <= -SIGMA) THEN
            i=idp(j,TT)-1
!*          irp(j,TT)=idp(j,TT)+1
            do while ((cpnp(i,j,TT) <= -sigma).and.(i >= nc1))
               i=i-1
            end do
            idp(j,TT)=max(i+1,nc1)
            ncc=min(irp(j,TT)+2,ncp)
            i=irp(j,TT)
            do while ((thickp(i,j,TT) > 0.d0).and.(i < ncc))
               i=i+1
            end do
            irp(j,TT)=i+1
!-----------------------------------------------------------------------------------------------!
         ELSEIF (THICKP(IDP(J,TT),J,TT) <= 0.D0) THEN
            i=idp(j,TT)
            do while ((thickp(i,j,TT) <= 0.d0).and.(i < ncp))
               i=i+1
            end do
            idp(j,TT)=i
            ncc=min(irp(j,TT)+2,ncp)
!*          i=irp(j)
            do while ((thickp(i,j,TT) > 0.d0).and.(i < ncc))
               i=i+1
            end do
            irp(j,TT)=i-1
!-----------------------------------------------------------------------------------------------!
         ELSE
            i=idp(j,TT)
            ncc=min(irp(j,TT)+2,ncp)
            do while ((thickp(i,j,TT) > 0.d0).and.(i < ncc))
               i=i+1
            end do
            irp(j,TT)=i !+1
         END IF
!-----------------------------------------------------------------------------------------------!
!    Sheet Cavitation                                                                           !
!-----------------------------------------------------------------------------------------------!
      ELSE !(IMCS > 0)
         I=NC1
         DO WHILE ((THICKP(I,J,TT) <= 0.D0).AND.(I < NCC))
            I=I+1
         END DO !((THICKP(I,J,TT) <= 0.D0).AND.(I < NCC))
         IDS(J,TT)=I
         IF (IDS(J,TT) == NCP) IDS(J,TT)=0
!-----------------------------------------------------------------------------------------------!
         IRS(J,TT)=0
         IF (IDS(J,TT) /= 0) THEN
            I=IDS(J,TT)
            DO WHILE ((THICKP(I,J,TT) > 0.D0).AND.(I < NCC))
               I=I+1
            END DO !((THICKP(I,J,TT) > 0.D0).AND.(I < NCC))
            IF ((THICKP(I,J,TT) > 0.D0).AND.(I == NCP)) THEN
               IRS(J,TT)=I
            ELSE !((THICKP(I,J,TT) > 0.D0).AND.(I == NCP))
               IRS(J,TT)=I-1
               if (irs(j,tt) < ids(j,tt)) then
                  ids(j,tt)=0
                  irs(j,tt)=0
               end if
            END IF !((THICKP(I,J,TT) > 0.D0).AND.(I == NCP))
         END IF !(IDS(J,TT) /= 0)
      END IF !(IMCS > 0)
!-----------------------------------------------------------------------------------------------!
!    Robustness for Kutta Condition                                                             !
!-----------------------------------------------------------------------------------------------!
!*    if (irs(j,tt) == ncp-1) then
!*       irs(j,tt)=ncp
!*       thickp(ncp,j,tt)=1.d-3
!*       cav1=1.d-3
!*    end if !(irp(j,tt) == 2)
      if ((irs(j,tt) > ncp-3).and.(irs(j,tt) /= ncp)) then
         do i=irs(j,tt)+1,ncp
            thickp(i,j,tt)=max(thickp(i,j,tt),1.d-3)
         end do !i=irs(j,tt)+1,ncp
         cav1=1.d-3
         irs(j,tt)=ncp
      end if !((irs(j,tt) > ncp-3).and.(irs(j,tt) /= ncp))
!-----------------------------------------------------------------------------------------------!
!    Cavity Thickness                                                                           !
!-----------------------------------------------------------------------------------------------!
      IF ((IDS(J,TT) == 0).AND.(IRS(J,TT) == 0)) THICKP(NC1:NCP,J,TT)=0.D0
      IF (IDS(J,TT) > NC1) THICKP(NC1:IDS(J,TT)-1,J,TT)=0.D0
      IF (IRS(J,TT) < NCP) THICKP(IRS(J,TT)+1:NCP,J,TT)=0.D0
!-----------------------------------------------------------------------------------------------!
!    Super Cavitation                                                                           !
!-----------------------------------------------------------------------------------------------!
      IF ((IRS(J,TT) == NCP).AND.(THICKP(NCP,J,TT) > 0.D0)) THEN
         J2=J-JI+1
         CAM1=CAV1/2.D0
         IF (IDPWS(J2,TT) == 0) THEN
            IDPWS(J2,TT)=1
            IRPWS(J2,TT)=2
         END IF !(IDPWS(J2,TT) == 0)
         DO I=IDPWS(J2,TT),IRPWS(J2,TT)
            IF (TT == 0) THEN
               CALL VWAKE(TT,XPW0(I,J2),YPW0(I,J2),ZPW0(I,J2),1,0,VWX,VWY,VWZ)
            ELSE !(TT)
               CALL VWAKE(TT,XPW0(I,J2),YPW0(I,J2),ZPW0(I,J2),1,IFREQ,VWX,VWY,VWZ)
            END IF !(TT)
            DETADT=0.D0
            DCAVDT=0.D0
            FRELAX=0.D0
            IF (TT > 0) THEN
               IF (THICKPWS(I,J,TT) == 0.D0) THEN
                  DETADT=0.D0
                  DCAVDT=0.D0
               ELSEIF (TT == 1) THEN
                  DETADT=THICKPWS (I,J2,1)/DTETA-THICKPWS (I,J2,0)/DTETA
                  DCAVDT=CAMBERPWS(I,J2,1)/DTETA-CAMBERPWS(I,J2,0)/DTETA
               ELSE !(TT == 1)
                  DETADT=THICKPWS (I,J2,TT-2)*0.5D0/DTETA- &
                         THICKPWS (I,J2,TT-1)*2.0D0/DTETA+ &
                         THICKPWS (I,J2,TT  )*1.5D0/DTETA
                  DCAVDT=CAMBERPWS(I,J2,TT-2)*0.5D0/DTETA- &
                         CAMBERPWS(I,J2,TT-1)*2.0D0/DTETA+ &
                         CAMBERPWS(I,J2,TT  )*1.5D0/DTETA
               END IF !(TT == 1)
               IF (TT <= (2*NTETA)) THEN
                  FRELAX=DFLOAT(TT-1)/DFLOAT(2*NTETA-1)
               ELSE !(TT <= (2*NTETA))
                  FRELAX=1.D0
               END IF !(TT <= (2*NTETA))
            END IF !(TT > 0)
!-----------------------------------------------------------------------------------------------!
            IF (IROTOR == 0) THEN
               VTT1S=VWX*UU(JJ)/PI*ET1XPW(I,J2)+ &             !x1 component of Total Velocity
                    (VWY*UU(JJ)/PI-ZPW0(I,J2))*ET1YPW(I,J2)+ & !y1 component of Total Velocity
                    (VWZ*UU(JJ)/PI+YPW0(I,J2))*ET1ZPW(I,J2)+ & !z1 component of Total Velocity
                    +VT1PW(I,J2)                               !s1 perturbation component
            ELSEIF (IROTOR == 1) THEN
               VTT1S=VWX/UU(JJ)*ET1XPW(I,J2)+ &             !x1 component of Total Velocity
                    (VWY/UU(JJ)-ZPW0(I,J2))*ET1YPW(I,J2)+ & !y1 component of Total Velocity
                    (VWZ/UU(JJ)+YPW0(I,J2))*ET1ZPW(I,J2)+ & !z1 component of Total Velocity
                    +VT1PW(I,J2)                            !s1 perturbation component
            END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
            SB=2.D0*DSQRT(AT1XPW(I,J2)**2+AT1YPW(I,J2)**2+AT1ZPW(I,J2)**2)
            CAV2=CAV1-SB/VTT1S*(SOURCEPWCAV(I,J2,TT)+DETADT*FRELAX) !Sign for suction side
            THICKPWS(I,J2,TT)=0.5D0*(CAV1+CAV2)
            CAV1=CAV2
!-----------------------------------------------------------------------------------------------!
            CAM2=CAM1-SB/VTT1S*DCAVDT*FRELAX
            CAMBERPWS(I,J2,TT)=0.5D0*(CAM1+CAM2)
            CAM1=CAM2
         END DO !I=IDPWS(J2,TT),IRPWS(J2,TT)
!-----------------------------------------------------------------------------------------------!
!    Cavity Extend                                                                              !
!-----------------------------------------------------------------------------------------------!
         I=MIN(IRPWS(J2,TT)+1,IABS(NCPW))
         IF (I == 2) THEN
            SA1=DSQRT(AT1XP (NCP,J )**2+AT1YP (NCP,J )**2+AT1ZP (NCP,J )**2)
            SB1=DSQRT(AT1XPW(I-1,J2)**2+AT1YPW(I-1,J2)**2+AT1ZPW(I-1,J2)**2)
            SC1=DSQRT(AT1XPW(I  ,J2)**2+AT1YPW(I  ,J2)**2+AT1ZPW(I  ,J2)**2)
            SA =SA1+SB1
            SB =SB1+SC1
            CAMBERPWS(I,J2,TT)=CAMBERPWS(I-1,J2,TT)
            THICKPWS (I,J2,TT)=0.5D0*THICKPWS(I-1,J2,TT)+CAMBERPWS(I-1,J2,TT)+ &
                              (0.5D0*THICKPWS(I-1,J2,TT)+CAMBERPWS(I-1,J2,TT)- &
                               THICKP(1,J,TT))/SA*SB
            THICKPWS (I,J2,TT)=2.D0*(THICKPWS(I  ,J2,TT)-CAMBERPWS(I  ,J2,TT))
         ELSE !(I == 2)
            SA1=DSQRT(AT1XPW(I-2,J2)**2+AT1YPW(I-2,J2)**2+AT1ZPW(I-2,J2)**2)
            SB1=DSQRT(AT1XPW(I-1,J2)**2+AT1YPW(I-1,J2)**2+AT1ZPW(I-1,J2)**2)
            SC1=DSQRT(AT1XPW(I  ,J2)**2+AT1YPW(I  ,J2)**2+AT1ZPW(I  ,J2)**2)
            SA =SA1+SB1
            SB =SB1+SC1
            CAMBERPWS(I,J2,TT)=CAMBERPWS(I-1,J2,TT)
            THICKPWS (I,J2,TT)=0.5D0*THICKPWS(I-1,J2,TT)+CAMBERPWS(I-1,J2,TT)+ &
                              (0.5D0*THICKPWS(I-1,J2,TT)+CAMBERPWS(I-1,J2,TT)- &
                               0.5D0*THICKPWS(I-2,J2,TT)-CAMBERPWS(I-2,J2,TT))/SA*SB
            THICKPWS (I,J2,TT)=2.D0*(THICKPWS(I  ,J2,TT)-CAMBERPWS(I  ,J2,TT))
         END IF !(I == 2)
!-----------------------------------------------------------------------------------------------!
!    Cavity Re-Attachment                                                                       !
!-----------------------------------------------------------------------------------------------!
         IDPWS(J2,TT)=1
         IF (THICKPWS(IRPWS(J2,TT),J2,TT) <= 0.D0) THEN
            THICKPWS(IRPWS(J2,TT),J2,TT)=0.D0
            IRPWS(J2,TT)=MAX(IRPWS(J2,TT)-1,1)
         ELSE !(THICKPWS(IRPWS(J,TT),J2,TT) <= 0.D0)
            IRPWS(J2,TT)=MIN(IRPWS(J2,TT)+1,IABS(NCPW))
         END IF !(THICKPWS(IRPWP(J2,TT),J2,TT) <= 0.D0)
         IF (IRPWS(J2,TT) < IABS(NCPW)) THEN
            THICKPWS (IRPWS(J2,TT)+1:IABS(NCPW),J2,TT)=0.D0
            CAMBERPWS(IRPWS(J2,TT)+1:IABS(NCPW),J2,TT)=0.D0
         END IF
      ELSE !((IRS(J,TT) == NCP).AND.(THICKP(NCP,J,TT) > 0.D0))
         J2=J-JI+1
         IDPWS(J2,TT)=0
         IRPWS(J2,TT)=0
         IF (NCPW > 0) THICKPWS (1:IABS(NCPW),J2,TT)=0.D0
         IF (NCPW > 0) CAMBERPWS(1:IABS(NCPW),J2,TT)=0.D0
      END IF !((IRS(J,TT) == NCP).AND.(THICKS(NCP,J,TT) > 0.D0))
!-----------------------------------------------------------------------------------------------!
   END IF !(IDS(J,TT) /= 0)
END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
DO J=1,NRP
   IF (IDS(J,TT) == 0) IRS(J,TT)=0
END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
!    Number of Cavitating Panels on Suction Side                                                !
!-----------------------------------------------------------------------------------------------!
NSCAV=0
DO J=1,NRP
   DO I=IDS(J,TT),IRS(J,TT)
      IF (I /= 0) NSCAV=NSCAV+1
   END DO !I=IDS(J,TT),IRS(J,TT)
END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
NWCAV=0
IF (NCPW > 0) THEN
   DO J=1,NRP
      J2=J-JI+1
      DO I=IDPWS(J2,TT),IRPWS(J2,TT)
         IF (I /= 0) NWCAV=NWCAV+1
      END DO !I=IDPWS(J2,TT),IRPWS(J2,TT)
   END DO !J=1,NRP
END IF !(NCPW > 0)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE CAVTHICKS
!-----------------------------------------------------------------------------------------------!
