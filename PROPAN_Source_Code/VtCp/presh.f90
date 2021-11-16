!-----------------------------------------------------------------------------------------------!
!    Computation of total velocities and pressure on the hub                                    !
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
SUBROUTINE PRESH(JJ,TT)
!-----------------------------------------------------------------------------------------------!
!    Created by: J.A.C. Falcao de Campos, IST                                                   !
!    Modified  : 22112013, J. Baltazar, version 1.0                                             !
!    Modified  : 28102014, J. Baltazar, correction for unsteady inflow                          !
!    Modified  : 10122014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!    Modified  : 09092015, J. Baltazar, correction in the gravity term                          !
!    Modified  : 06062016, J. Baltazar, 2016 version 1.3, new relation for unsteady terms       !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,JJ,TT
DOUBLE PRECISION :: VWX,VWY,VWZ,RH0,TH0,DPHIDT,VTTSQ,UINFSQ,FRELAX,SDIV
!-----------------------------------------------------------------------------------------------!
DO J=1,NHTP
   DO I=1,NHX
      IF (TT == 0) THEN
         CALL VWAKE(TT,XH0(I,J),YH0(I,J),ZH0(I,J),1,0,VWX,VWY,VWZ)
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
         IF (IROTOR == 0) THEN
            VTT1H(I,J)=VWX*UU(JJ)/PI*ET1XH(I,J)+ &           !x1 component of Total Velocity
                      (VWY*UU(JJ)/PI-ZH0(I,J))*ET1YH(I,J)+ & !y1 component of Total Velocity
                      (VWZ*UU(JJ)/PI+YH0(I,J))*ET1ZH(I,J)+ & !z1 component of Total Velocity
                      +VT1H(I,J)                             !s1 perturbation component
            VTT2H(I,J)=VWX*UU(JJ)/PI*ET2XH(I,J)+ &           !x2 component of Total Velocity
                      (VWY*UU(JJ)/PI-ZH0(I,J))*ET2YH(I,J)+ & !y2 component of Total Velocity
                      (VWZ*UU(JJ)/PI+YH0(I,J))*ET2ZH(I,J)+ & !z2 component of Total Velocity
                      +VT2H(I,J)                             !s2 perturbation component
            VTTSQ=VTT1H(I,J)*VTT1H(I,J)+VTT2H(I,J)*VTT2H(I,J)
            UINFSQ=(VWX*UU(JJ)/PI)**2+(VWY*UU(JJ)/PI-ZH0(I,J))**2+(VWZ*UU(JJ)/PI+YH0(I,J))**2
            IF (ICP == 1) CPH(I,J,TT)=1.D0-VTTSQ/UINFSQ-SDIV(2.D0,FN*FN,0.D0)*YH0(I,J)/UINFSQ
            IF (ICP == 2) CPH(I,J,TT)=UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YH0(I,J)
            CPNH(I,J,TT)=PI*PI*(UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YH0(I,J))
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
         ELSEIF (IROTOR == 1) THEN
            VTT1H(I,J)=VWX/UU(JJ)*ET1XH(I,J)+ &           !x1 component of Total Velocity
                      (VWY/UU(JJ)-ZH0(I,J))*ET1YH(I,J)+ & !y1 component of Total Velocity
                      (VWZ/UU(JJ)+YH0(I,J))*ET1ZH(I,J)+ & !z1 component of Total Velocity
                      +VT1H(I,J)                          !s1 perturbation component
            VTT2H(I,J)=VWX/UU(JJ)*ET2XH(I,J)+ &           !x2 component of Total Velocity
                      (VWY/UU(JJ)-ZH0(I,J))*ET2YH(I,J)+ & !y2 component of Total Velocity
                      (VWZ/UU(JJ)+YH0(I,J))*ET2ZH(I,J)+ & !z2 component of Total Velocity
                      +VT2H(I,J)                          !s2 perturbation component
            VTTSQ=VTT1H(I,J)*VTT1H(I,J)+VTT2H(I,J)*VTT2H(I,J)
            UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ)-ZH0(I,J))**2+(VWZ/UU(JJ)+YH0(I,J))**2
            IF (ICP == 1) CPH(I,J,TT)=1.D0-VTTSQ/UINFSQ-SDIV(2.D0,FN*FN,0.D0)*YH0(I,J)/UINFSQ
            IF (ICP == 2) CPH(I,J,TT)=UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YH0(I,J)
            CPNH(I,J,TT)=PI*PI*(UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YH0(I,J))
         END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
      ELSE !(TT == 0)
         CALL VWAKE(TT,XH0(I,J),YH0(I,J),ZH0(I,J),1,IFREQ,VWX,VWY,VWZ)
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
         IF (IROTOR == 0) THEN
            VTT1H(I,J)=VWX*UU(JJ)/PI*ET1XH(I,J)+ &           !x1 component of Total Velocity
                      (VWY*UU(JJ)/PI-ZH0(I,J))*ET1YH(I,J)+ & !y1 component of Total Velocity
                      (VWZ*UU(JJ)/PI+YH0(I,J))*ET1ZH(I,J)+ & !z1 component of Total Velocity
                      +VT1H(I,J)                             !s1 perturbation component
            VTT2H(I,J)=VWX*UU(JJ)/PI*ET2XH(I,J)+ &           !x2 component of Total Velocity
                      (VWY*UU(JJ)/PI-ZH0(I,J))*ET2YH(I,J)+ & !y2 component of Total Velocity
                      (VWZ*UU(JJ)/PI+YH0(I,J))*ET2ZH(I,J)+ & !z2 component of Total Velocity
                      +VT2H(I,J)                             !s2 perturbation component
            RH0=DSQRT(YH0(I,J)*YH0(I,J)+ZH0(I,J)*ZH0(I,J))
            TH0=DATAN2(ZH0(I,J),YH0(I,J))
            VTTSQ=VTT1H(I,J)*VTT1H(I,J)+VTT2H(I,J)*VTT2H(I,J)
            UINFSQ=(VWX*UU(JJ)/PI)**2+(VWY*UU(JJ)/PI-ZH0(I,J))**2+(VWZ*UU(JJ)/PI+YH0(I,J))**2
            IF (TT == 1) THEN
               DPHIDT=POTH(I,J,1)/DTETA-POTH(I,J,0)/DTETA
            ELSE !(TT == 1)
               DPHIDT=POTH(I,J,TT-2)*0.5D0/DTETA- &
                      POTH(I,J,TT-1)*2.D0 /DTETA+ &
                      POTH(I,J,TT  )*1.5D0/DTETA
            END IF !(TT == 1)
            IF (TT <= (2*NTETA)) THEN
               FRELAX=(DFLOAT(TT)-1.D0)/DFLOAT(2*NTETA-1)
            ELSE !(TT <= (2*NTETA))
               FRELAX=1.D0
            END IF !(TT <= (2*NTETA))
            IF (ICP == 1) CPH(I,J,TT)=1.D0-VTTSQ/UINFSQ-2.D0*DPHIDT*FRELAX/UINFSQ &
                             -SDIV(2.D0,FN*FN,0.D0)*RH0*DCOS(TH0-DTETA*(DFLOAT(TT)-1.D0))/UINFSQ
            IF (ICP == 2) CPH(I,J,TT)=UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RH0*DCOS(TH0-DTETA*(DFLOAT(TT)-1.D0))
            CPNH(I,J,TT)=PI*PI*(UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                   -SDIV(2.D0,FN*FN,0.D0)*RH0*DCOS(TH0-DTETA*(DFLOAT(TT)-1.D0)))
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
         ELSEIF (IROTOR == 1) THEN
            VTT1H(I,J)=VWX/UU(JJ)*ET1XH(I,J)+ &           !x1 component of Total Velocity
                      (VWY/UU(JJ)-ZH0(I,J))*ET1YH(I,J)+ & !y1 component of Total Velocity
                      (VWZ/UU(JJ)+YH0(I,J))*ET1ZH(I,J)+ & !z1 component of Total Velocity
                      +VT1H(I,J)                          !s1 perturbation component
            VTT2H(I,J)=VWX/UU(JJ)*ET2XH(I,J)+ &           !x2 component of Total Velocity
                      (VWY/UU(JJ)-ZH0(I,J))*ET2YH(I,J)+ & !y2 component of Total Velocity
                      (VWZ/UU(JJ)+YH0(I,J))*ET2ZH(I,J)+ & !z2 component of Total Velocity
                      +VT2H(I,J)                          !s2 perturbation component
            RH0=DSQRT(YH0(I,J)*YH0(I,J)+ZH0(I,J)*ZH0(I,J))
            TH0=DATAN2(ZH0(I,J),YH0(I,J))
            VTTSQ=VTT1H(I,J)*VTT1H(I,J)+VTT2H(I,J)*VTT2H(I,J)
            UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ)-ZH0(I,J))**2+(VWZ/UU(JJ)+YH0(I,J))**2
            IF (TT == 1) THEN
               DPHIDT=POTH(I,J,1)/DTETA-POTH(I,J,0)/DTETA
            ELSE !(TT == 1)
               DPHIDT=POTH(I,J,TT-2)*0.5D0/DTETA- &
                      POTH(I,J,TT-1)*2.D0 /DTETA+ &
                      POTH(I,J,TT  )*1.5D0/DTETA
            END IF !(TT == 1)
            IF (TT <= (2*NTETA)) THEN
               FRELAX=(DFLOAT(TT)-1.D0)/DFLOAT(2*NTETA-1)
            ELSE !(TT <= (2*NTETA))
               FRELAX=1.D0
            END IF !(TT <= (2*NTETA))
            IF (ICP == 1) CPH(I,J,TT)=1.D0-VTTSQ/UINFSQ-2.D0*DPHIDT*FRELAX/UINFSQ &
                             -SDIV(2.D0,FN*FN,0.D0)*RH0*DCOS(TH0-DTETA*(DFLOAT(TT)-1.D0))/UINFSQ
            IF (ICP == 2) CPH(I,J,TT)=UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RH0*DCOS(TH0-DTETA*(DFLOAT(TT)-1.D0))
            CPNH(I,J,TT)=PI*PI*(UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                   -SDIV(2.D0,FN*FN,0.D0)*RH0*DCOS(TH0-DTETA*(DFLOAT(TT)-1.D0)))
         END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
      END IF !(TT == 0)
   END DO !I=1,NHX
END DO !J=1,NHTP
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PRESH
!-----------------------------------------------------------------------------------------------!
