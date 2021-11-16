!-----------------------------------------------------------------------------------------------!
!    Computation of total velocities and pressure on the propeller blade                        !
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
SUBROUTINE PRESP(JJ,TT)
!-----------------------------------------------------------------------------------------------!
!    Created by: J.A.C. Falcao de Campos, IST                                                   !
!    Modified  : 15112013, J. Baltazar, version 1.0                                             !
!    Modified  : 28102014, J. Baltazar, Cavitation Model                                        !
!    Modified  : 10122014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!    Modified  : 07012015, J. Baltazar, 2015 version 1.0, Wing Case                             !
!    Modified  : 09092015, J. Baltazar, correction in the gravity term                          !
!    Modified  : 22042016, J. Baltazar, cavity recover                                          !
!    Modified  : 06062016, J. Baltazar, 2016 version 1.3, new relation for unsteady terms       !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,JJ,TT
DOUBLE PRECISION :: VWX,VWY,VWZ,RP0,TP0,DPHIDT,VTTSQ,UINFSQ,FRELAX,SDIV,SIGMAREC
!-----------------------------------------------------------------------------------------------!
DO J=1,NRP
   DO I=1,NCP
      IF (TT == 0) THEN
         CALL VWAKE(TT,XP0(I,J),YP0(I,J),ZP0(I,J),1,0,VWX,VWY,VWZ)
!-----------------------------------------------------------------------------------------------!
!    Wing Case                                                                                  !
!-----------------------------------------------------------------------------------------------!
         IF (IROTOR == -1) THEN
            VTT1P(I,J)=DCOSD(UU(JJ))*ET1XP(I,J)+DSIND(UU(JJ))*ET1ZP(I,J)+VT1P(I,J)
            VTT2P(I,J)=DCOSD(UU(JJ))*ET2XP(I,J)+DSIND(UU(JJ))*ET2ZP(I,J)+VT2P(I,J)
            VTTSQ=VTT1P(I,J)*VTT1P(I,J)+VTT2P(I,J)*VTT2P(I,J)
            UINFSQ=1.D0
            IF (ICP == 1) CPP(I,J,TT)=1.D0-VTTSQ/UINFSQ
            IF (ICP == 2) CPP(I,J,TT)=0.D0
            CPNP(I,J,TT)=CPP(I,J,TT) !Forces Computation
!*          IF ((I <= IDP(J,TT)).AND.(I >= IRP(J,TT))) THEN !Pressure on the Cavity
!*             IF (ICP == 1) CPP(I,J,TT)=-SIGMA/(PI*PI*UINFSQ)
!*             IF (ICP == 2) CPP(I,J,TT)=-SIGMA/(PI*PI) !SIGMAREC(I,J)
!*             CPNP(I,J,TT)=-SIGMA
!*          ELSEIF ((I >= IDS(J,TT)).AND.(I <= IRS(J,TT))) THEN
!*             IF (ICP == 1) CPP(I,J,TT)=-SIGMA/(PI*PI*UINFSQ)
!*             IF (ICP == 2) CPP(I,J,TT)=-SIGMA/(PI*PI) !SIGMAREC(I,J)
!*             CPNP(I,J,TT)=-SIGMA
!*          END IF
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
         ELSEIF (IROTOR == 0) THEN
            VTT1P(I,J)=VWX*UU(JJ)/PI*ET1XP(I,J)+ &           !x1 component of Total Velocity
                      (VWY*UU(JJ)/PI-ZP0(I,J))*ET1YP(I,J)+ & !y1 component of Total Velocity
                      (VWZ*UU(JJ)/PI+YP0(I,J))*ET1ZP(I,J)+ & !z1 component of Total Velocity
                      +VT1P(I,J)                             !s1 perturbation component
            VTT2P(I,J)=VWX*UU(JJ)/PI*ET2XP(I,J)+ &           !x2 component of Total Velocity
                      (VWY*UU(JJ)/PI-ZP0(I,J))*ET2YP(I,J)+ & !y2 component of Total Velocity
                      (VWZ*UU(JJ)/PI+YP0(I,J))*ET2ZP(I,J)+ & !z2 component of Total Velocity
                      +VT2P(I,J)                             !s2 perturbation component
            VTTSQ=VTT1P(I,J)*VTT1P(I,J)+VTT2P(I,J)*VTT2P(I,J)
            UINFSQ=(VWX*UU(JJ)/PI)**2+(VWY*UU(JJ)/PI-ZP0(I,J))**2+(VWZ*UU(JJ)/PI+YP0(I,J))**2
            IF (ICP == 1) CPP(I,J,TT)=1.D0-VTTSQ/UINFSQ-SDIV(2.D0,FN*FN,0.D0)*YP0(I,J)/UINFSQ
            IF (ICP == 2) CPP(I,J,TT)=UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YP0(I,J)
            CPNP(I,J,TT)=PI*PI*(UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YP0(I,J))
            IF ((I <= IDP(J,TT)).AND.(I >= IRP(J,TT))) THEN !Pressure on the Cavity
               IF (ICP == 1) CPP(I,J,TT)=-SIGMAREC(I,J,TT)/(PI*PI*UINFSQ) !SIGMA
               IF (ICP == 2) CPP(I,J,TT)=-SIGMAREC(I,J,TT)/(PI*PI) !SIGMA
               CPNP(I,J,TT)=-SIGMAREC(I,J,TT) !SIGMA
            ELSEIF ((I >= IDS(J,TT)).AND.(I <= IRS(J,TT))) THEN
               IF (ICP == 1) CPP(I,J,TT)=-SIGMAREC(I,J,TT)/(PI*PI*UINFSQ) !SIGMA
               IF (ICP == 2) CPP(I,J,TT)=-SIGMAREC(I,J,TT)/(PI*PI) !SIGMA
               CPNP(I,J,TT)=-SIGMAREC(I,J,TT) !SIGMA
            END IF
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
         ELSEIF (IROTOR == 1) THEN
            VTT1P(I,J)=VWX/UU(JJ)*ET1XP(I,J)+ &           !x1 component of Total Velocity
                      (VWY/UU(JJ)-ZP0(I,J))*ET1YP(I,J)+ & !y1 component of Total Velocity
                      (VWZ/UU(JJ)+YP0(I,J))*ET1ZP(I,J)+ & !z1 component of Total Velocity
                      +VT1P(I,J)                          !s1 perturbation component
            VTT2P(I,J)=VWX/UU(JJ)*ET2XP(I,J)+ &           !x2 component of Total Velocity
                      (VWY/UU(JJ)-ZP0(I,J))*ET2YP(I,J)+ & !y2 component of Total Velocity
                      (VWZ/UU(JJ)+YP0(I,J))*ET2ZP(I,J)+ & !z2 component of Total Velocity
                      +VT2P(I,J)                          !s2 perturbation component
            VTTSQ=VTT1P(I,J)*VTT1P(I,J)+VTT2P(I,J)*VTT2P(I,J)
            UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ)-ZP0(I,J))**2+(VWZ/UU(JJ)+YP0(I,J))**2
            IF (ICP == 1) CPP(I,J,TT)=1.D0-VTTSQ/UINFSQ-SDIV(2.D0,FN*FN,0.D0)*YP0(I,J)/UINFSQ
            IF (ICP == 2) CPP(I,J,TT)=UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YP0(I,J)
            CPNP(I,J,TT)=PI*PI*(UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YP0(I,J))
            IF ((I <= IDP(J,TT)).AND.(I >= IRP(J,TT))) THEN !Pressure on the Cavity
               IF (ICP == 1) CPP(I,J,TT)=-SIGMAREC(I,J,TT)/(PI*PI*UINFSQ) !SIGMA
               IF (ICP == 2) CPP(I,J,TT)=-SIGMAREC(I,J,TT)/(PI*PI) !SIGMA
               CPNP(I,J,TT)=-SIGMAREC(I,J,TT) !SIGMA
            ELSEIF ((I >= IDS(J,TT)).AND.(I <= IRS(J,TT))) THEN
               IF (ICP == 1) CPP(I,J,TT)=-SIGMAREC(I,J,TT)/(PI*PI*UINFSQ) !SIGMA
               IF (ICP == 2) CPP(I,J,TT)=-SIGMAREC(I,J,TT)/(PI*PI) !SIGMA
               CPNP(I,J,TT)=-SIGMAREC(I,J,TT) !SIGMA
            END IF
         END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
      ELSE !(TT == 0)
         CALL VWAKE(TT,XP0(I,J),YP0(I,J),ZP0(I,J),1,IFREQ,VWX,VWY,VWZ)
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
         IF (IROTOR == 0) THEN
            VTT1P(I,J)=VWX*UU(JJ)/PI*ET1XP(I,J)+ &           !x1 component of Total Velocity
                      (VWY*UU(JJ)/PI-ZP0(I,J))*ET1YP(I,J)+ & !y1 component of Total Velocity
                      (VWZ*UU(JJ)/PI+YP0(I,J))*ET1ZP(I,J)+ & !z1 component of Total Velocity
                      +VT1P(I,J)                             !s1 perturbation component
            VTT2P(I,J)=VWX*UU(JJ)/PI*ET2XP(I,J)+ &           !x2 component of Total Velocity
                      (VWY*UU(JJ)/PI-ZP0(I,J))*ET2YP(I,J)+ & !y2 component of Total Velocity
                      (VWZ*UU(JJ)/PI+YP0(I,J))*ET2ZP(I,J)+ & !z2 component of Total Velocity
                      +VT2P(I,J)                             !s2 perturbation component
            RP0=DSQRT(YP0(I,J)*YP0(I,J)+ZP0(I,J)*ZP0(I,J))
            TP0=DATAN2(ZP0(I,J),YP0(I,J))
            VTTSQ=VTT1P(I,J)*VTT1P(I,J)+VTT2P(I,J)*VTT2P(I,J)
            UINFSQ=(VWX*UU(JJ)/PI)**2+(VWY*UU(JJ)/PI-ZP0(I,J))**2+(VWZ*UU(JJ)/PI+YP0(I,J))**2
            IF (TT == 1) THEN
               DPHIDT=POTP(I,J,1)/DTETA-POTP(I,J,0)/DTETA
            ELSE !(TT == 1)
               DPHIDT=POTP(I,J,TT-2)*0.5D0/DTETA- &
                      POTP(I,J,TT-1)*2.D0 /DTETA+ &
                      POTP(I,J,TT  )*1.5D0/DTETA
            END IF !(TT == 1)
            IF (TT <= (2*NTETA)) THEN
               FRELAX=DFLOAT(TT-1)/DFLOAT(2*NTETA-1)
            ELSE !(TT <= (2*NTETA))
               FRELAX=1.D0
            END IF !(TT <= (2*NTETA))
            IF (ICP == 1) CPP(I,J,TT)=1.D0-VTTSQ/UINFSQ-2.D0*DPHIDT*FRELAX/UINFSQ &
                             -SDIV(2.D0,FN*FN,0.D0)*RP0*DCOS(TP0-DTETA*DFLOAT(TT-1))/UINFSQ
            IF (ICP == 2) CPP(I,J,TT)=UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RP0*DCOS(TP0-DTETA*DFLOAT(TT-1))
            CPNP(I,J,TT)=PI*PI*(UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                   -SDIV(2.D0,FN*FN,0.D0)*RP0*DCOS(TP0-DTETA*DFLOAT(TT-1)))
            IF ((I <= IDP(J,TT)).AND.(I >= IRP(J,TT))) THEN !Pressure on the Cavity
               IF (ICP == 1) CPP(I,J,TT)=-SIGMA/(PI*PI*UINFSQ) !SIGMAREC(I,J,TT)
               IF (ICP == 2) CPP(I,J,TT)=-SIGMA/(PI*PI) !SIGMAREC(I,J,TT)
               CPNP(I,J,TT)=-SIGMA !SIGMAREC(I,J,TT)
            ELSEIF ((I >= IDS(J,TT)).AND.(I <= IRS(J,TT))) THEN
               IF (ICP == 1) CPP(I,J,TT)=-SIGMA/(PI*PI*UINFSQ) !SIGMAREC(I,J,TT)
               IF (ICP == 2) CPP(I,J,TT)=-SIGMA/(PI*PI) !SIGMAREC(I,J,TT)
               CPNP(I,J,TT)=-SIGMA !SIGMAREC(I,J,TT)
            END IF
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
         ELSEIF (IROTOR == 1) THEN
            VTT1P(I,J)=VWX/UU(JJ)*ET1XP(I,J)+ &           !x1 component of Total Velocity
                      (VWY/UU(JJ)-ZP0(I,J))*ET1YP(I,J)+ & !y1 component of Total Velocity
                      (VWZ/UU(JJ)+YP0(I,J))*ET1ZP(I,J)+ & !z1 component of Total Velocity
                      +VT1P(I,J)                          !s1 perturbation component
            VTT2P(I,J)=VWX/UU(JJ)*ET2XP(I,J)+ &           !x2 component of Total Velocity
                      (VWY/UU(JJ)-ZP0(I,J))*ET2YP(I,J)+ & !y2 component of Total Velocity
                      (VWZ/UU(JJ)+YP0(I,J))*ET2ZP(I,J)+ & !z2 component of Total Velocity
                      +VT2P(I,J)                          !s2 perturbation component
            RP0=DSQRT(YP0(I,J)*YP0(I,J)+ZP0(I,J)*ZP0(I,J))
            TP0=DATAN2(ZP0(I,J),YP0(I,J))
            VTTSQ=VTT1P(I,J)*VTT1P(I,J)+VTT2P(I,J)*VTT2P(I,J)
            UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ)-ZP0(I,J))**2+(VWZ/UU(JJ)+YP0(I,J))**2
            IF (TT == 1) THEN
               DPHIDT=POTP(I,J,1)/DTETA-POTP(I,J,0)/DTETA
            ELSE !(TT == 1)
               DPHIDT=POTP(I,J,TT-2)*0.5D0/DTETA- &
                      POTP(I,J,TT-1)*2.D0 /DTETA+ &
                      POTP(I,J,TT  )*1.5D0/DTETA
            END IF !(TT == 1)
            IF (TT <= (2*NTETA)) THEN
               FRELAX=DFLOAT(TT-1)/DFLOAT(2*NTETA-1)
            ELSE !(TT <= (2*NTETA))
               FRELAX=1.D0
            END IF !(TT <= (2*NTETA))
            IF (ICP == 1) CPP(I,J,TT)=1.D0-VTTSQ/UINFSQ-2.D0*DPHIDT*FRELAX/UINFSQ &
                             -SDIV(2.D0,FN*FN,0.D0)*RP0*DCOS(TP0-DTETA*DFLOAT(TT-1))/UINFSQ
            IF (ICP == 2) CPP(I,J,TT)=UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RP0*DCOS(TP0-DTETA*DFLOAT(TT-1))
            CPNP(I,J,TT)=PI*PI*(UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                   -SDIV(2.D0,FN*FN,0.D0)*RP0*DCOS(TP0-DTETA*DFLOAT(TT-1)))
            IF ((I <= IDP(J,TT)).AND.(I >= IRP(J,TT))) THEN !Pressure on the Cavity
               IF (ICP == 1) CPP(I,J,TT)=-SIGMA/(PI*PI*UINFSQ) !SIGMAREC(I,J,TT)
               IF (ICP == 2) CPP(I,J,TT)=-SIGMA/(PI*PI) !SIGMAREC(I,J,TT)
               CPNP(I,J,TT)=-SIGMA !SIGMAREC(I,J,TT)
            ELSEIF ((I >= IDS(J,TT)).AND.(I <= IRS(J,TT))) THEN
               IF (ICP == 1) CPP(I,J,TT)=-SIGMA/(PI*PI*UINFSQ) !SIGMAREC(I,J,TT)
               IF (ICP == 2) CPP(I,J,TT)=-SIGMA/(PI*PI) !SIGMAREC(I,J,TT)
               CPNP(I,J,TT)=-SIGMA !SIGMAREC(I,J,TT)
            END IF
         END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
      END IF !(TT == 0)
   END DO !I=1,NCP
END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PRESP
!-----------------------------------------------------------------------------------------------!
