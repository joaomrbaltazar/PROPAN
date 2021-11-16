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
SUBROUTINE CAVTHICKP(JJ,TT,CC)
!-----------------------------------------------------------------------------------------------!
!    Created by: 07112014, J. Baltazar, Cavitation Model                                        !
!    Modified  : 11112014, J. Baltazar, version 3.1, Unsteady Cavitation Model                  !
!    Modified  : 19112014, J. Baltazar, version 3.2, Mid-chord cavitation                       !
!    Modified  : 27112014, J. Baltazar, version 3.3, Super-Cavitation Model                     !
!    Modified  : 09122014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!    Modified  : 06032015, J. Baltazar, Robustness for Kutta Condition                          !
!    Modified  : 29022016, J. Baltazar, Revision                                                !
!    Modified  : 12042016, J. Baltazar, 2016 version 1.1                                        !
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
   IF (IDP(J,TT) /= 0) THEN
!-----------------------------------------------------------------------------------------------!
      CAV1=0.D0
      DO I=IDP(J,TT),IRP(J,TT),-1
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
            ELSE
               DETADT=THICKP(I,J,TT-2)*0.5D0/DTETA- &
                      THICKP(I,J,TT-1)*2.0D0/DTETA+ &
                      THICKP(I,J,TT  )*1.5D0/DTETA
            END IF
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
            VTT2S=VWX/UU(JJ)*T2X+ &                  !x1 component of Total Velocity
                 (VWY/UU(JJ)-ZP0(I,J))*T2Y+ &        !y1 component of Total Velocity
                 (VWZ/UU(JJ)+YP0(I,J))*T2Z+ &        !z1 component of Total Velocity
                 +VT2S(I,J)                          !s1 perturbation component
            VTTNP=VWX/UU(JJ)*UNXP0(I,J)+ &           !x2 component of Total Velocity
                 (VWY/UU(JJ)-ZP0(I,J))*UNYP0(I,J)+ & !y2 component of Total Velocity
                 (VWZ/UU(JJ)+YP0(I,J))*UNZP0(I,J)+ & !z2 component of Total Velocity
                 -SOURCEP(I,J,1)                     !s2 perturbation component
         END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
         A=VTT2S*ET12-VTT1P(I,J) !Pressure Side
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
            IF ((I <= IDP(J-1,TT)).AND.(I >= IRP(J-1,TT))) THEN
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
         ELSEIF ((I <= IDP(J-1,TT)).AND.(I >= IRP(J-1,TT)).AND. &
                                                 (I <= IDP(J-2,TT)).AND.(I >= IRP(J-2,TT))) THEN
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
         ELSEIF ((I <= IDP(J-1,TT)).AND.(I >= IRP(J-1,TT))) THEN
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
         if ((i <= idp(j,tt)).and.(i >= idp(j,tt)-imcp).and.(thickp(i,j,tt) < 0.d0)) then !***Correction***
            write(30,'(A,I4,A,I4)') ' CC=',cc,', TT=',tt
            write(30,'(A,I4,I4)  ') ' Eta > 0 at the K.B.C., (I,J)=',i,j
            cav2=dabs(cav2)
            thickp(i,j,tt)=dabs(thickp(i,j,tt))
         end if
         CAV1=CAV2
      END DO !I=IDP(J,TT),IRP(J,TT),-1
!-----------------------------------------------------------------------------------------------!
!    Cavity Extend                                                                              !
!-----------------------------------------------------------------------------------------------!
      IF ((THICKP(IRP(J,TT),J,TT) > 0.D0).AND. &
                                  (THICKP(IRP(J,TT),J,TT)-THICKP(IRP(J,TT)+1,J,TT) < 0.D0)) THEN
         NCC=MAX(IRP(J,TT)-3-2*NCP/100,1)
      ELSE !((THICKP(IRP(J,TT),J,TT) > 0.D0).AND. &
!                                      (THICKP(IRP(J,TT),J,TT)-THICKP(IRP(J,TT)+1,J,TT) < 0.D0))
         NCC=MAX(IRP(J,TT)-2-2*NCP/100,1)
      END IF !((THICKP(IRP(J,TT),J,TT) > 0.D0).AND. &
!                                      (THICKP(IRP(J,TT),J,TT)-THICKP(IRP(J,TT)+1,J,TT) < 0.D0))
!-----------------------------------------------------------------------------------------------!
      IF (IRP(J,TT) > 1) THEN
         DO I=(IRP(J,TT)-1),1,-1
            SA1=DSQRT(AT1XP(I+2,J)**2+AT1YP(I+2,J)**2+AT1ZP(I+2,J)**2) !***
            SB1=DSQRT(AT1XP(I+1,J)**2+AT1YP(I+1,J)**2+AT1ZP(I+1,J)**2)
            SC1=DSQRT(AT1XP(I  ,J)**2+AT1YP(I  ,J)**2+AT1ZP(I  ,J)**2)
            SA =SA1+SB1
            SB =SB1+SC1
            THICKP(I,J,TT)=THICKP(I+1,J,TT)+(THICKP(I+1,J,TT)-THICKP(I+2,J,TT))/SA*SB
            IF (THICKP(I,J,TT) < 0.D0) THICKP(I,J,TT)=0.D0
         END DO !I=(IRP(J,TT)-1),1,-1
      END IF !(IRP(J,TT) > 1)
!-----------------------------------------------------------------------------------------------!
      IF (IDP(J,TT) < NC) THEN
         DO I=IDP(J,TT)+1,NC
            SA1=DSQRT(AT1XP(I-2,J)**2+AT1YP(I-2,J)**2+AT1ZP(I-2,J)**2)
            SB1=DSQRT(AT1XP(I-1,J)**2+AT1YP(I-1,J)**2+AT1ZP(I-1,J)**2)
            SC1=DSQRT(AT1XP(I  ,J)**2+AT1YP(I  ,J)**2+AT1ZP(I  ,J)**2)
            SA =SA1+SB1
            SB =SB1+SC1
            THICKP(I,J,TT)=THICKP(I-1,J,TT)+(THICKP(I-1,J,TT)-THICKP(I-2,J,TT))/SA*SB
            IF (THICKP(I,J,TT) < 0.D0) THICKP(I,J,TT)=0.D0
         END DO !I=IDP(J,TT)+1,NC
      END IF !(IDP(J,TT) < NC)
!-----------------------------------------------------------------------------------------------!
!    Mid-Chord Cavitation                                                                       !
!-----------------------------------------------------------------------------------------------!
      IF (IMCP > 0) THEN
         IF (CPNP(IDP(J,TT)+1,J,TT) <= -SIGMA) THEN
            i=idp(j,TT)+1
!*          irp(j,TT)=idp(j,TT)-1
            do while ((cpnp(i,j,TT) <= -sigma).and.(i <= nc))
               i=i+1
            end do
            idp(j,TT)=min(i-1,nc)
            ncc=max(irp(j,TT)-2,1)
            i=irp(j,TT)
            do while ((thickp(i,j,TT) > 0.d0).and.(i > ncc))
               i=i-1
            end do
            irp(j,TT)=i+1
!-----------------------------------------------------------------------------------------------!
         ELSEIF (THICKP(IDP(J,TT),J,TT) <= 0.D0) THEN
            i=idp(j,TT)
            do while ((thickp(i,j,TT) <= 0.d0).and.(i > 1))
               i=i-1
            end do
            idp(j,TT)=i
            ncc=max(irp(j,TT)-2,1)
!*          i=irp(j)
            do while ((thickp(i,j,TT) > 0.d0).and.(i > ncc))
               i=i-1
            end do
            irp(j,TT)=i+1
!-----------------------------------------------------------------------------------------------!
         ELSE
            i=idp(j,TT)
            ncc=max(irp(j,TT)-2,1)
            do while ((thickp(i,j,TT) > 0.d0).and.(i > ncc))
               i=i-1
            end do
            irp(j,TT)=i !+1
         END IF
!-----------------------------------------------------------------------------------------------!
!    Sheet Cavitation                                                                           !
!-----------------------------------------------------------------------------------------------!
      ELSE !(IMCP > 0)
         I=NC
         DO WHILE ((THICKP(I,J,TT) <= 0.D0).AND.(I > NCC))
            I=I-1
         END DO !((THICKP(I,J,TT) <= 0.D0).AND.(I > NCC))
         IDP(J,TT)=I
         IF (IDP(J,TT) == 1) IDP(J,TT)=0
!-----------------------------------------------------------------------------------------------!
         IRP(J,TT)=0
         IF (IDP(J,TT) /= 0) THEN
            I=IDP(J,TT)
            DO WHILE ((THICKP(I,J,TT) > 0.D0).AND.(I > NCC))
               I=I-1
            END DO !((THICKP(I,J,TT) > 0.D0).AND.(I > NCC))
            IF ((THICKP(I,J,TT) > 0.D0).AND.(I == 1)) THEN
               IRP(J,TT)=I
            ELSE !((THICKP(I,J,TT) > 0.D0).AND.(I == 1))
               IRP(J,TT)=I+1
               if (irp(j,tt) > idp(j,tt)) then
                  idp(j,tt)=0
                  irp(j,tt)=0
               end if
            END IF !((THICKP(I,J,TT) > 0.D0).AND.(I == 1))
         END IF !(IDS(J,TT) /= 0)
      END IF !(IMCP > 0)
!-----------------------------------------------------------------------------------------------!
!    Robustness for Kutta Condition                                                             !
!-----------------------------------------------------------------------------------------------!
!*    if (irp(j,tt) == 2) then
!*       irp(j,tt)=1
!*       thickp(1,j,tt)=1.d-3
!*       cav1=1.d-3
!*    end if !(irp(j,tt) == 2)
      if ((irp(j,tt) < 4).and.(irp(j,tt) /= 1)) then
         do i=irp(j,tt)-1,1,-1
            thickp(i,j,tt)=max(thickp(i,j,tt),1.d-3)
         end do !i=irp(j,tt)-1,1,-1
         cav1=1.d-3
         irp(j,tt)=1
      end if !((irp(j,tt) < 4).and.(irp(j,tt) /= 1))
!-----------------------------------------------------------------------------------------------!
!    Cavity Thickness                                                                           !
!-----------------------------------------------------------------------------------------------!
      IF ((IDP(J,TT) == 0).AND.(IRP(J,TT) == 0)) THICKP(NC:1:-1,J,TT)=0.D0
      IF (IDP(J,TT) < NC) THICKP(NC:IDP(J,TT)+1:-1,J,TT)=0.D0
      IF (IRP(J,TT) > 1 ) THICKP(IRP(J,TT)-1 :1:-1,J,TT)=0.D0
!-----------------------------------------------------------------------------------------------!
!    Super Cavitation                                                                           !
!-----------------------------------------------------------------------------------------------!
      IF ((IRP(J,TT) == 1).AND.(THICKP(1,J,TT) > 0.D0)) THEN
         J2=J-JI+1
         CAM1=CAV1/2.D0
         IF (IDPWP(J2,TT) == 0) THEN
            IDPWP(J2,TT)=1
            IRPWP(J2,TT)=2
         END IF !(IDPWP(J2,TT) == 0)
         DO I=IDPWP(J2,TT),IRPWP(J2,TT)
            IF (TT == 0) THEN
               CALL VWAKE(TT,XPW0(I,J2),YPW0(I,J2),ZPW0(I,J2),1,0,VWX,VWY,VWZ)
            ELSE !(TT)
               CALL VWAKE(TT,XPW0(I,J2),YPW0(I,J2),ZPW0(I,J2),1,IFREQ,VWX,VWY,VWZ)
            END IF !(TT)
            DETADT=0.D0
            DCAVDT=0.D0
            FRELAX=0.D0
            IF (TT > 0) THEN
               IF (THICKPWP(I,J2,TT) == 0.D0) THEN
                  DETADT=0.D0
                  DCAVDT=0.D0
               ELSEIF (TT == 1) THEN
                  DETADT=THICKPWP (I,J2,1)/DTETA-THICKPWP (I,J2,0)/DTETA
                  DCAVDT=CAMBERPWP(I,J2,1)/DTETA-CAMBERPWP(I,J2,0)/DTETA
               ELSE !(TT == 1)
                  DETADT=THICKPWP (I,J2,TT-2)*0.5D0/DTETA- &
                         THICKPWP (I,J2,TT-1)*2.0D0/DTETA+ &
                         THICKPWP (I,J2,TT  )*1.5D0/DTETA
                  DCAVDT=CAMBERPWP(I,J2,TT-2)*0.5D0/DTETA- &
                         CAMBERPWP(I,J2,TT-1)*2.0D0/DTETA+ &
                         CAMBERPWP(I,J2,TT  )*1.5D0/DTETA
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
            CAV2=CAV1+SB/VTT1S*(SOURCEPWCAV(I,J2,TT)-DETADT*FRELAX) !Sign for pressure side
            THICKPWP(I,J2,TT)=0.5D0*(CAV1+CAV2)
            CAV1=CAV2
!-----------------------------------------------------------------------------------------------!
            CAM2=CAM1-SB/VTT1S*DCAVDT*FRELAX
            CAMBERPWP(I,J2,TT)=0.5D0*(CAM1+CAM2)
            CAM1=CAM2
         END DO !I=IDPWP(J2,TT),IRPWP(J2,TT)
!-----------------------------------------------------------------------------------------------!
!    Cavity Extend                                                                              !
!-----------------------------------------------------------------------------------------------!
         I=MIN(IRPWP(J2,TT)+1,IABS(NCPW))
         IF (I == 2) THEN
            SA1=DSQRT(AT1XP (1  ,J )**2+AT1YP (1  ,J )**2+AT1ZP (1  ,J )**2)
            SB1=DSQRT(AT1XPW(I-1,J2)**2+AT1YPW(I-1,J2)**2+AT1ZPW(I-1,J2)**2)
            SC1=DSQRT(AT1XPW(I  ,J2)**2+AT1YPW(I  ,J2)**2+AT1ZPW(I  ,J2)**2)
            SA =SA1+SB1
            SB =SB1+SC1
            CAMBERPWP(I,J2,TT)=CAMBERPWP(I-1,J2,TT)
            THICKPWP (I,J2,TT)=0.5D0*THICKPWP(I-1,J2,TT)+CAMBERPWP(I-1,J2,TT)+ &
                              (0.5D0*THICKPWP(I-1,J2,TT)+CAMBERPWP(I-1,J2,TT)- &
                               THICKP(1,J,TT))/SA*SB
            THICKPWP (I,J2,TT)=2.D0*(THICKPWP(I  ,J2,TT)-CAMBERPWP(I  ,J2,TT))
         ELSE !(I == 2)
            SA1=DSQRT(AT1XPW(I-2,J2)**2+AT1YPW(I-2,J2)**2+AT1ZPW(I-2,J2)**2)
            SB1=DSQRT(AT1XPW(I-1,J2)**2+AT1YPW(I-1,J2)**2+AT1ZPW(I-1,J2)**2)
            SC1=DSQRT(AT1XPW(I  ,J2)**2+AT1YPW(I  ,J2)**2+AT1ZPW(I  ,J2)**2)
            SA =SA1+SB1
            SB =SB1+SC1
            CAMBERPWP(I,J2,TT)=CAMBERPWP(I-1,J2,TT)
            THICKPWP (I,J2,TT)=0.5D0*THICKPWP(I-1,J2,TT)+CAMBERPWP(I-1,J2,TT)+ &
                              (0.5D0*THICKPWP(I-1,J2,TT)+CAMBERPWP(I-1,J2,TT)- &
                               0.5D0*THICKPWP(I-2,J2,TT)-CAMBERPWP(I-2,J2,TT))/SA*SB
            THICKPWP (I,J2,TT)=2.D0*(THICKPWP(I  ,J2,TT)-CAMBERPWP(I  ,J2,TT))
         END IF !(I == 2)
!-----------------------------------------------------------------------------------------------!
!    Cavity Re-Attachment                                                                       !
!-----------------------------------------------------------------------------------------------!
         IDPWP(J2,TT)=1
         IF (THICKPWP(IRPWP(J2,TT),J2,TT) <= 0.D0) THEN
            THICKPWP(IRPWP(J2,TT),J2,TT)=0.D0
            IRPWP(J2,TT)=MAX(IRPWP(J2,TT)-1,1)
         ELSE !(THICKPWP(IRPWP(J2,TT),J2,TT) <= 0.D0)
            IRPWP(J2,TT)=MIN(IRPWP(J2,TT)+1,IABS(NCPW))
         END IF !(THICKPWP(IRPWP(J2,TT),J2,TT) <= 0.D0)
         IF (IRPWP(J2,TT) < IABS(NCPW)) THEN
            THICKPWP (IRPWP(J2,TT)+1:IABS(NCPW),J2,TT)=0.D0
            CAMBERPWP(IRPWP(J2,TT)+1:IABS(NCPW),J2,TT)=0.D0
         END IF
      ELSE !((IRP(J,TT) == 1).AND.(THICKP(1,J,TT) > 0.D0))
         J2=J-JI+1
         IDPWP(J2,TT)=0
         IRPWP(J2,TT)=0
         IF (NCPW < 0) THICKPWP (1:IABS(NCPW),J2,TT)=0.D0
         IF (NCPW < 0) CAMBERPWP(1:IABS(NCPW),J2,TT)=0.D0
      END IF !((IRP(J,TT) == 1).AND.(THICKP(1,J,TT) > 0.D0))
!-----------------------------------------------------------------------------------------------!
   END IF !(IDP(J,TT) /= 0)
END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
DO J=1,NRP
   IF (IDP(J,TT) == 0) IRP(J,TT)=0
END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
!    Number of Cavitating Panels on Pressure Side                                               !
!-----------------------------------------------------------------------------------------------!
NPCAV=0
DO J=1,NRP
   DO I=IDP(J,TT),IRP(J,TT),-1
      IF (I /= 0) NPCAV=NPCAV+1
   END DO !I=IDP(J,TT),IRP(J,TT),-1
END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
NWCAV=0
IF (NCPW < 0) THEN
   DO J=1,NRP
      J2=J-JI+1
      DO I=IDPWP(J2,TT),IRPWP(J2,TT)
         IF (I /= 0) NWCAV=NWCAV+1
      END DO !I=IDPWP(J2,TT),IRPWP(J2,TT)
   END DO !J=1,NRP
END IF !(NCPW < 0)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE CAVTHICKP
!-----------------------------------------------------------------------------------------------!
