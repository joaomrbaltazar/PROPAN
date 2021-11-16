!-----------------------------------------------------------------------------------------------!
!    Computation of total velocities and pressure on the nozzle                                 !
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
SUBROUTINE PRESN(JJ,TT)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 20112013, J. Baltazar, 2014 version 1.0                                        !
!    Modified  : 28102014, J. Baltazar, Correction for unsteady inflow                          !
!    Modified  : 10122014, J. Baltazar, 2014 version 3.4, Unsteady Super-Cavitation Model       !
!    Modified  : 06012015, J. Baltazar, 2015 version 1.0, Revision                              !
!    Modified  : 09092015, J. Baltazar, correction in the gravity term                          !
!    Modified  : 06062016, J. Baltazar, 2016 version 1.3, new relation for unsteady terms       !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,JJ,TT
DOUBLE PRECISION :: VWX,VWY,VWZ,RN0,TN0,DPHIDT,VTTSQ,UINFSQ,FRELAX,SDIV
!-----------------------------------------------------------------------------------------------!
DO J=1,NNTP
   DO I=1,NNXT1
      IF (TT == 0) THEN
         CALL VWAKE(TT,XN0(I,J),YN0(I,J),ZN0(I,J),1,0,VWX,VWY,VWZ)
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
         IF (IROTOR == 0) THEN
            VTT1N(I,J)=VWX*UU(JJ)/PI*ET1XN(I,J)+ &           !x1 component of Total Velocity
                      (VWY*UU(JJ)/PI-ZN0(I,J))*ET1YN(I,J)+ & !y1 component of Total Velocity
                      (VWZ*UU(JJ)/PI+YN0(I,J))*ET1ZN(I,J)+ & !z1 component of Total Velocity
                      +VT1N(I,J)                             !s1 perturbation component
            VTT2N(I,J)=VWX*UU(JJ)/PI*ET2XN(I,J)+ &           !x2 component of Total Velocity
                      (VWY*UU(JJ)/PI-ZN0(I,J))*ET2YN(I,J)+ & !y2 component of Total Velocity
                      (VWZ*UU(JJ)/PI+YN0(I,J))*ET2ZN(I,J)+ & !z2 component of Total Velocity
                      +VT2N(I,J)                             !s2 perturbation component
            VTTSQ=VTT1N(I,J)*VTT1N(I,J)+VTT2N(I,J)*VTT2N(I,J)
            UINFSQ=(VWX*UU(JJ)/PI)**2+(VWY*UU(JJ)/PI-ZN0(I,J))**2+(VWZ*UU(JJ)/PI+YN0(I,J))**2
            IF (ICP == 1) CPN(I,J,TT)=1.D0-VTTSQ/UINFSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(I,J)/UINFSQ
            IF (ICP == 2) CPN(I,J,TT)=UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(I,J)
            CPNN(I,J,TT)=PI*PI*(UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(I,J))
!-----------------------------------------------------------------------------------------------!
!    Without Rotational Velocity                                                                !
!-----------------------------------------------------------------------------------------------!
            IF (IP == 0) THEN
               UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ))**2+(VWZ/UU(JJ))**2
               IF (ICP == 1) CPN(I,J,TT)=1.D0-VTTSQ/UINFSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(I,J)/UINFSQ
               IF (ICP == 2) CPN(I,J,TT)=UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(I,J)
               CPNN(I,J,TT)=PI*PI*(UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(I,J))
            END IF !(IP == 0)
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
         ELSEIF (IROTOR == 1) THEN
            VTT1N(I,J)=VWX/UU(JJ)*ET1XN(I,J)+ &           !x1 component of Total Velocity
                      (VWY/UU(JJ)-ZN0(I,J))*ET1YN(I,J)+ & !y1 component of Total Velocity
                      (VWZ/UU(JJ)+YN0(I,J))*ET1ZN(I,J)+ & !z1 component of Total Velocity
                      +VT1N(I,J)                          !s1 perturbation component
            VTT2N(I,J)=VWX/UU(JJ)*ET2XN(I,J)+ &           !x2 component of Total Velocity
                      (VWY/UU(JJ)-ZN0(I,J))*ET2YN(I,J)+ & !y2 component of Total Velocity
                      (VWZ/UU(JJ)+YN0(I,J))*ET2ZN(I,J)+ & !z2 component of Total Velocity
                      +VT2N(I,J)                          !s2 perturbation component
            VTTSQ=VTT1N(I,J)*VTT1N(I,J)+VTT2N(I,J)*VTT2N(I,J)
            UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ)-ZN0(I,J))**2+(VWZ/UU(JJ)+YN0(I,J))**2
            IF (ICP == 1) CPN(I,J,TT)=1.D0-VTTSQ/UINFSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(I,J)/UINFSQ
            IF (ICP == 2) CPN(I,J,TT)=UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(I,J)
            CPNN(I,J,TT)=PI*PI*(UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(I,J))
!-----------------------------------------------------------------------------------------------!
!    Without Rotational Velocity                                                                !
!-----------------------------------------------------------------------------------------------!
            IF (IP == 0) THEN
               UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ))**2+(VWZ/UU(JJ))**2
               IF (ICP == 1) CPN(I,J,TT)=1.D0-VTTSQ/UINFSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(I,J)/UINFSQ
               IF (ICP == 2) CPN(I,J,TT)=UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(I,J)
               CPNN(I,J,TT)=PI*PI*(UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(I,J))
            END IF !(IP == 0)
         END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
      ELSE !(TT == 0)
         CALL VWAKE(TT,XN0(I,J),YN0(I,J),ZN0(I,J),1,IFREQ,VWX,VWY,VWZ)
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
         IF (IROTOR == 0) THEN
            VTT1N(I,J)=VWX*UU(JJ)/PI*ET1XN(I,J)+ &           !x1 component of Total Velocity
                      (VWY*UU(JJ)/PI-ZN0(I,J))*ET1YN(I,J)+ & !y1 component of Total Velocity
                      (VWZ*UU(JJ)/PI+YN0(I,J))*ET1ZN(I,J)+ & !z1 component of Total Velocity
                      +VT1N(I,J)                             !s1 perturbation component
            VTT2N(I,J)=VWX*UU(JJ)/PI*ET2XN(I,J)+ &           !x2 component of Total Velocity
                      (VWY*UU(JJ)/PI-ZN0(I,J))*ET2YN(I,J)+ & !y2 component of Total Velocity
                      (VWZ*UU(JJ)/PI+YN0(I,J))*ET2ZN(I,J)+ & !z2 component of Total Velocity
                      +VT2N(I,J)                             !s2 perturbation component
            RN0=DSQRT(YN0(I,J)*YN0(I,J)+ZN0(I,J)*ZN0(I,J))
            TN0=DATAN2(ZN0(I,J),YN0(I,J))
            VTTSQ=VTT1N(I,J)*VTT1N(I,J)+VTT2N(I,J)*VTT2N(I,J)
            UINFSQ=(VWX*UU(JJ)/PI)**2+(VWY*UU(JJ)/PI-ZN0(I,J))**2+(VWZ*UU(JJ)/PI+YN0(I,J))**2
            IF (TT == 1) THEN
               DPHIDT=POTN(I,J,1)/DTETA-POTN(I,J,0)/DTETA
            ELSE !(TT == 1)
               DPHIDT=POTN(I,J,TT-2)*0.5D0/DTETA- &
                      POTN(I,J,TT-1)*2.D0 /DTETA+ &
                      POTN(I,J,TT  )*1.5D0/DTETA
            END IF !(TT == 1)
            IF (TT <= (2*NTETA)) THEN
               FRELAX=(DFLOAT(TT)-1.D0)/DFLOAT(2*NTETA-1)
            ELSE !(TT <= (2*NTETA))
               FRELAX=1.D0
            END IF !(TT <= (2*NTETA))
            IF (ICP == 1) CPN(I,J,TT)=1.D0-VTTSQ/UINFSQ-2.D0*DPHIDT*FRELAX/UINFSQ &
                             -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))/UINFSQ
            IF (ICP == 2) CPN(I,J,TT)=UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))
            CPNN(I,J,TT)=PI*PI*(UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                   -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0)))
!-----------------------------------------------------------------------------------------------!
!    Without Rotational Velocity                                                                !
!-----------------------------------------------------------------------------------------------!
            IF (IP == 0) THEN
               UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ))**2+(VWZ/UU(JJ))**2
               IF (TT == 1) THEN
                  DPHIDT=POTN(I,J,1)/DTETA-POTN(I,J,0)/DTETA
               ELSE !(TT == 1)
                  DPHIDT=POTN(I,J,TT-2)*0.5D0/DTETA- &
                         POTN(I,J,TT-1)*2.D0 /DTETA+ &
                         POTN(I,J,TT  )*1.5D0/DTETA
               END IF !(TT == 1)
               IF (TT <= (2*NTETA)) THEN
                  FRELAX=(DFLOAT(TT)-1.D0)/DFLOAT(2*NTETA-1)
               ELSE !(TT <= (2*NTETA))
                  FRELAX=1.D0
               END IF !(TT <= (2*NTETA))
               IF (ICP == 1) CPN(I,J,TT)=1.D0-VTTSQ/UINFSQ-2.D0*DPHIDT*FRELAX/UINFSQ &
                             -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))/UINFSQ
               IF (ICP == 2) CPN(I,J,TT)=UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))
               CPNN(I,J,TT)=PI*PI*(UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                   -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0)))
            END IF !(IP == 0)
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
         ELSEIF (IROTOR == 1) THEN
            VTT1N(I,J)=VWX/UU(JJ)*ET1XN(I,J)+ &           !x1 component of Total Velocity
                      (VWY/UU(JJ)-ZN0(I,J))*ET1YN(I,J)+ & !y1 component of Total Velocity
                      (VWZ/UU(JJ)+YN0(I,J))*ET1ZN(I,J)+ & !z1 component of Total Velocity
                      +VT1N(I,J)                          !s1 perturbation component
            VTT2N(I,J)=VWX/UU(JJ)*ET2XN(I,J)+ &           !x2 component of Total Velocity
                      (VWY/UU(JJ)-ZN0(I,J))*ET2YN(I,J)+ & !y2 component of Total Velocity
                      (VWZ/UU(JJ)+YN0(I,J))*ET2ZN(I,J)+ & !z2 component of Total Velocity
                      +VT2N(I,J)                          !s2 perturbation component
            RN0=DSQRT(YN0(I,J)*YN0(I,J)+ZN0(I,J)*ZN0(I,J))
            TN0=DATAN2(ZN0(I,J),YN0(I,J))
            VTTSQ=VTT1N(I,J)*VTT1N(I,J)+VTT2N(I,J)*VTT2N(I,J)
            UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ)-ZN0(I,J))**2+(VWZ/UU(JJ)+YN0(I,J))**2
            IF (TT == 1) THEN
               DPHIDT=POTN(I,J,1)/DTETA-POTN(I,J,0)/DTETA
            ELSE !(TT == 1)
               DPHIDT=POTN(I,J,TT-2)*0.5D0/DTETA- &
                      POTN(I,J,TT-1)*2.D0 /DTETA+ &
                      POTN(I,J,TT  )*1.5D0/DTETA
            END IF !(TT == 1)
            IF (TT <= (2*NTETA)) THEN
               FRELAX=(DFLOAT(TT)-1.D0)/DFLOAT(2*NTETA-1)
            ELSE !(TT <= (2*NTETA))
               FRELAX=1.D0
            END IF !(TT <= (2*NTETA))
            IF (ICP == 1) CPN(I,J,TT)=1.D0-VTTSQ/UINFSQ-2.D0*DPHIDT*FRELAX/UINFSQ &
                             -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))/UINFSQ
            IF (ICP == 2) CPN(I,J,TT)=UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))
            CPNN(I,J,TT)=PI*PI*(UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                   -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0)))
!-----------------------------------------------------------------------------------------------!
!    Without Rotational Velocity                                                                !
!-----------------------------------------------------------------------------------------------!
            IF (IP == 0) THEN
               UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ))**2+(VWZ/UU(JJ))**2
               IF (TT == 1) THEN
                  DPHIDT=POTN(I,J,1)/DTETA-POTN(I,J,0)/DTETA
               ELSE !(TT == 1)
                  DPHIDT=POTN(I,J,TT-2)*0.5D0/DTETA- &
                         POTN(I,J,TT-1)*2.D0 /DTETA+ &
                         POTN(I,J,TT  )*1.5D0/DTETA
               END IF !(TT == 1)
               IF (TT <= (2*NTETA)) THEN
                  FRELAX=(DFLOAT(TT)-1.D0)/DFLOAT(2*NTETA-1)
               ELSE !(TT <= (2*NTETA))
                  FRELAX=1.D0
               END IF !(TT <= (2*NTETA))
               IF (ICP == 1) CPN(I,J,TT)=1.D0-VTTSQ/UINFSQ-2.D0*DPHIDT*FRELAX/UINFSQ &
                             -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))/UINFSQ
               IF (ICP == 2) CPN(I,J,TT)=UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))
               CPNN(I,J,TT)=PI*PI*(UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                   -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0)))
            END IF !(IP == 0)
         END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
      END IF !(TT == 0)
   END DO !I=1,NNXT1
END DO !J=1,NNTP
!-----------------------------------------------------------------------------------------------!
!    Velocity and Pressure Correction                                                           !
!-----------------------------------------------------------------------------------------------!
IF (TT == 0) THEN
   DO J=1,NNTP
      CALL VWAKE(TT,XN0(NNU+1,J),YN0(NNU+1,J),ZN0(NNU+1,J),1,0,VWX,VWY,VWZ)
      VTT1N(NNU+1,J)=0.5D0*(VTT1N(NNU,J)+VTT1N(NNU+2,J))
      VTT2N(NNU+1,J)=0.5D0*(VTT2N(NNU,J)+VTT2N(NNU+2,J))
      VTTSQ=VTT1N(NNU+1,J)*VTT1N(NNU+1,J)+VTT2N(NNU+1,J)*VTT2N(NNU+1,J)
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
      IF (IROTOR == 0) THEN
         UINFSQ=(VWX*UU(JJ)/PI             )**2+ &
                (VWY*UU(JJ)/PI-ZN0(NNU+1,J))**2+ &
                (VWZ*UU(JJ)/PI+YN0(NNU+1,J))**2
         IF (ICP == 1) CPN(NNU+1,J,TT)=1.D0-VTTSQ/UINFSQ &
                                                      -SDIV(2.D0,FN*FN,0.D0)*YN0(NNU+1,J)/UINFSQ
         IF (ICP == 2) CPN(NNU+1,J,TT)=UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NNU+1,J)
         CPNN(NNU+1,J,TT)=PI*PI*(UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NNU+1,J))
!-----------------------------------------------------------------------------------------------!
!    Without Rotational Velocity                                                                !
!-----------------------------------------------------------------------------------------------!
         IF (IP == 0) THEN
            UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ))**2+(VWZ/UU(JJ))**2
            IF (ICP == 1) CPN(NNU+1,J,TT)=1.D0-VTTSQ/UINFSQ &
                                                      -SDIV(2.D0,FN*FN,0.D0)*YN0(NNU+1,J)/UINFSQ
            IF (ICP == 2) CPN(NNU+1,J,TT)=UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NNU+1,J)
            CPNN(NNU+1,J,TT)=PI*PI*(UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NNU+1,J))
         END IF !(IP == 0)
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
      ELSEIF (IROTOR == 1) THEN
         UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ)-ZN0(NNU+1,J))**2+(VWZ/UU(JJ)+YN0(NNU+1,J))**2
         IF (ICP == 1) CPN(NNU+1,J,TT)=1.D0-VTTSQ/UINFSQ &
                                                      -SDIV(2.D0,FN*FN,0.D0)*YN0(NNU+1,J)/UINFSQ
         IF (ICP == 2) CPN(NNU+1,J,TT)=UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NNU+1,J)
         CPNN(NNU+1,J,TT)=PI*PI*(UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NNU+1,J))
!-----------------------------------------------------------------------------------------------!
!    Without Rotational Velocity                                                                !
!-----------------------------------------------------------------------------------------------!
         IF (IP == 0) THEN
            UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ))**2+(VWZ/UU(JJ))**2
            IF (ICP == 1) CPN(NNU+1,J,TT)=1.D0-VTTSQ/UINFSQ &
                                                      -SDIV(2.D0,FN*FN,0.D0)*YN0(NNU+1,J)/UINFSQ
            IF (ICP == 2) CPN(NNU+1,J,TT)=UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NNU+1,J)
            CPNN(NNU+1,J,TT)=PI*PI*(UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NNU+1,J))
         END IF !(IP == 0)
      END IF !(IROTOR)
   END DO !J=1,NNTP
ELSE !(TT == 0)
   DO J=1,NNTP
      CALL VWAKE(TT,XN0(NNU+1,J),YN0(NNU+1,J),ZN0(NNU+1,J),1,IFREQ,VWX,VWY,VWZ)
      RN0=DSQRT(YN0(NNU+1,J)*YN0(NNU+1,J)+ZN0(NNU+1,J)*ZN0(NNU+1,J))
      TN0=DATAN2(ZN0(NNU+1,J),YN0(NNU+1,J))
      VTT1N(NNU+1,J)=0.5D0*(VTT1N(NNU,J)+VTT1N(NNU+2,J))
      VTT2N(NNU+1,J)=0.5D0*(VTT2N(NNU,J)+VTT2N(NNU+2,J))
      VTTSQ=VTT1N(NNU+1,J)*VTT1N(NNU+1,J)+VTT2N(NNU+1,J)*VTT2N(NNU+1,J)
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
      IF (IROTOR == 0) THEN
         UINFSQ=(VWX*UU(JJ)/PI             )**2+ &
                (VWY*UU(JJ)/PI-ZN0(NNU+1,J))**2+ &
                (VWZ*UU(JJ)/PI+YN0(NNU+1,J))**2
         IF (TT == 1) THEN
            DPHIDT=POTN(NNU+1,J,1)/DTETA-POTN(NNU+1,J,0)/DTETA
         ELSE !(TT == 1)
            DPHIDT=POTN(NNU+1,J,TT-2)*0.5D0/DTETA- &
                   POTN(NNU+1,J,TT-1)*2.D0 /DTETA+ &
                   POTN(NNU+1,J,TT  )*1.5D0/DTETA
         END IF !(TT == 1)
         IF (TT <= (2*NTETA)) THEN
            FRELAX=(DFLOAT(TT)-1.D0)/DFLOAT(2*NTETA-1)
         ELSE !(TT <= (2*NTETA))
            FRELAX=1.D0
         END IF !(TT <= (2*NTETA))
         IF (ICP == 1) CPN(NNU+1,J,TT)=1.D0-VTTSQ/UINFSQ-2.D0*DPHIDT*FRELAX/UINFSQ &
                             -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))/UINFSQ
         IF (ICP == 2) CPN(NNU+1,J,TT)=UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))
         CPNN(NNU+1,J,TT)=PI*PI*(UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                   -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0)))
!-----------------------------------------------------------------------------------------------!
!    Without Rotational Velocity                                                                !
!-----------------------------------------------------------------------------------------------!
         IF (IP == 0) THEN
            UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ))**2+(VWZ/UU(JJ))**2
            IF (ICP == 1) CPN(NNU+1,J,TT)=1.D0-VTTSQ/UINFSQ-2.D0*DPHIDT*FRELAX/UINFSQ &
                             -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))/UINFSQ
            IF (ICP == 2) CPN(NNU+1,J,TT)=UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))
            CPNN(NNU+1,J,TT)=PI*PI*(UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                   -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0)))
         END IF !(IP == 0)
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
      ELSEIF (IROTOR == 1) THEN
         UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ)-ZN0(NNU+1,J))**2+(VWZ/UU(JJ)+YN0(NNU+1,J))**2
         IF (TT == 1) THEN
            DPHIDT=POTN(NNU+1,J,1)/DTETA-POTN(NNU+1,J,0)/DTETA
         ELSE !(TT == 1)
            DPHIDT=POTN(NNU+1,J,TT-2)*0.5D0/DTETA- &
                   POTN(NNU+1,J,TT-1)*2.D0 /DTETA+ &
                   POTN(NNU+1,J,TT  )*1.5D0/DTETA
         END IF !(TT == 1)
         IF (TT <= (2*NTETA)) THEN
            FRELAX=(DFLOAT(TT)-1.D0)/DFLOAT(2*NTETA-1)
         ELSE !(TT <= (2*NTETA))
            FRELAX=1.D0
         END IF !(TT <= (2*NTETA))
         IF (ICP == 1) CPN(NNU+1,J,TT)=1.D0-VTTSQ/UINFSQ-2.D0*DPHIDT*FRELAX/UINFSQ &
                             -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))/UINFSQ
         IF (ICP == 2) CPN(NNU+1,J,TT)=UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))
         CPNN(NNU+1,J,TT)=PI*PI*(UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                   -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0)))
!-----------------------------------------------------------------------------------------------!
!    Without Rotational Velocity                                                                !
!-----------------------------------------------------------------------------------------------!
         IF (IP == 0) THEN
            UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ))**2+(VWZ/UU(JJ))**2
            IF (ICP == 1) CPN(NNU+1,J,TT)=1.D0-VTTSQ/UINFSQ-2.D0*DPHIDT*FRELAX/UINFSQ &
                             -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))/UINFSQ
            IF (ICP == 2) CPN(NNU+1,J,TT)=UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))
            CPNN(NNU+1,J,TT)=PI*PI*(UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                   -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0)))
         END IF !(IP == 0)
      END IF !(IROTOR)
   END DO !J=1,NNTP
END IF !(TT == 0)
!-----------------------------------------------------------------------------------------------!
IF (TT == 0) THEN
   DO J=1,NNTP
      CALL VWAKE(TT,XN0(NN2,J),YN0(NN2,J),ZN0(NN2,J),1,0,VWX,VWY,VWZ)
      RN0=DSQRT(YN0(NN2,J)*YN0(NN2,J)+ZN0(NN2,J)*ZN0(NN2,J))
      TN0=DATAN2(ZN0(NN2,J),YN0(NN2,J))
      VTT1N(NN2,J)=0.5D0*(VTT1N(NN2-1,J)+VTT1N(NN2+1,J))
      VTT2N(NN2,J)=0.5D0*(VTT2N(NN2-1,J)+VTT2N(NN2+1,J))
      VTTSQ=VTT1N(NN2,J)*VTT1N(NN2,J)+VTT2N(NN2,J)*VTT2N(NN2,J)
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
      IF (IROTOR == 0) THEN
         UINFSQ=(VWX*UU(JJ)/PI)**2+(VWY*UU(JJ)/PI-ZN0(NN2,J))**2+(VWZ*UU(JJ)/PI+YN0(NN2,J))**2
         IF (ICP == 1) CPN(NN2,J,TT)=1.D0-VTTSQ/UINFSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NN2,J)/UINFSQ
         IF (ICP == 2) CPN(NN2,J,TT)=UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NN2,J)
         CPNN(NN2,J,TT)=PI*PI*(UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NN2,J))
!-----------------------------------------------------------------------------------------------!
!    Without Rotational Velocity                                                                !
!-----------------------------------------------------------------------------------------------!
         IF (IP == 0) THEN
            UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ))**2+(VWZ/UU(JJ))**2
            IF (ICP == 1) CPN(NN2,J,TT)=1.D0-VTTSQ/UINFSQ &
                                                        -SDIV(2.D0,FN*FN,0.D0)*YN0(NN2,J)/UINFSQ
            IF (ICP == 2) CPN(NN2,J,TT)=UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NN2,J)
            CPNN(NN2,J,TT)=PI*PI*(UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NN2,J))
         END IF !(IP == 0)
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
      ELSEIF (IROTOR == 1) THEN
         UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ)-ZN0(NN2,J))**2+(VWZ/UU(JJ)+YN0(NN2,J))**2
         IF (ICP == 1) CPN(NN2,J,TT)=1.D0-VTTSQ/UINFSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NN2,J)/UINFSQ
         IF (ICP == 2) CPN(NN2,J,TT)=UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NN2,J)
         CPNN(NN2,J,TT)=PI*PI*(UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NN2,J))
!-----------------------------------------------------------------------------------------------!
!    Without Rotational Velocity                                                                !
!-----------------------------------------------------------------------------------------------!
         IF (IP == 0) THEN
            UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ))**2+(VWZ/UU(JJ))**2
            IF (ICP == 1) CPN(NN2,J,TT)=1.D0-VTTSQ/UINFSQ &
                                                        -SDIV(2.D0,FN*FN,0.D0)*YN0(NN2,J)/UINFSQ
            IF (ICP == 2) CPN(NN2,J,TT)=UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NN2,J)
            CPNN(NN2,J,TT)=PI*PI*(UINFSQ-VTTSQ-SDIV(2.D0,FN*FN,0.D0)*YN0(NN2,J))
         END IF !(IP == 0)
      END IF !(IROTOR)
   END DO !J=1,NNTP
ELSE !(TT == 0)
   DO J=1,NNTP
      CALL VWAKE(TT,XN0(NN2,J),YN0(NN2,J),ZN0(NN2,J),1,IFREQ,VWX,VWY,VWZ)
      RN0=DSQRT(YN0(NN2,J)*YN0(NN2,J)+ZN0(NN2,J)*ZN0(NN2,J))
      TN0=DATAN2(ZN0(NN2,J),YN0(NN2,J))
      VTT1N(NN2,J)=0.5D0*(VTT1N(NN2-1,J)+VTT1N(NN2+1,J))
      VTT2N(NN2,J)=0.5D0*(VTT2N(NN2-1,J)+VTT2N(NN2+1,J))
      VTTSQ=VTT1N(NN2,J)*VTT1N(NN2,J)+VTT2N(NN2,J)*VTT2N(NN2,J)
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
      IF (IROTOR == 0) THEN
         UINFSQ=(VWX*UU(JJ)/PI)**2+(VWY*UU(JJ)/PI-ZN0(NN2,J))**2+(VWZ*UU(JJ)/PI+YN0(NN2,J))**2
         IF (TT == 1) THEN
            DPHIDT=POTN(NN2,J,1)/DTETA-POTN(NN2,J,0)/DTETA
         ELSE !(TT == 1)
            DPHIDT=POTN(NN2,J,TT-2)*0.5D0/DTETA- &
                   POTN(NN2,J,TT-1)*2.D0 /DTETA+ &
                   POTN(NN2,J,TT  )*1.5D0/DTETA
         END IF !(TT == 1)
         IF (TT <= (2*NTETA)) THEN
            FRELAX=(DFLOAT(TT)-1.D0)/DFLOAT(2*NTETA-1)
         ELSE !(TT <= (2*NTETA))
            FRELAX=1.D0
         END IF !(TT <= (2*NTETA))
         IF (ICP == 1) CPN(NN2,J,TT)=1.D0-VTTSQ/UINFSQ-2.D0*DPHIDT*FRELAX/UINFSQ &
                             -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))/UINFSQ
         IF (ICP == 2) CPN(NN2,J,TT)=UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))
         CPNN(NN2,J,TT)=PI*PI*(UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX&
                                   -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0)))
!-----------------------------------------------------------------------------------------------!
!    Without Rotational Velocity                                                                !
!-----------------------------------------------------------------------------------------------!
         IF (IP == 0) THEN
            UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ))**2+(VWZ/UU(JJ))**2
            IF (ICP == 1) CPN(NN2,J,TT)=1.D0-VTTSQ/UINFSQ-2.D0*DPHIDT*FRELAX/UINFSQ &
                             -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))/UINFSQ
            IF (ICP == 2) CPN(NN2,J,TT)=UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))
            CPNN(NN2,J,TT)=PI*PI*(UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                   -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0)))
         END IF !(IP == 0)
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
      ELSEIF (IROTOR == 1) THEN
         UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ)-ZN0(NN2,J))**2+(VWZ/UU(JJ)+YN0(NN2,J))**2
         IF (TT == 1) THEN
            DPHIDT=POTN(NN2,J,1)/DTETA-POTN(NN2,J,0)/DTETA
         ELSE !(TT == 1)
            DPHIDT=POTN(NN2,J,TT-2)*0.5D0/DTETA- &
                   POTN(NN2,J,TT-1)*2.D0 /DTETA+ &
                   POTN(NN2,J,TT  )*1.5D0/DTETA
         END IF !(TT == 1)
         IF (TT <= (2*NTETA)) THEN
            FRELAX=(DFLOAT(TT)-1.D0)/DFLOAT(2*NTETA-1)
         ELSE !(TT <= (2*NTETA))
            FRELAX=1.D0
         END IF !(TT <= (2*NTETA))
         IF (ICP == 1) CPN(NN2,J,TT)=1.D0-VTTSQ/UINFSQ-2.D0*DPHIDT*FRELAX/UINFSQ &
                             -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))/UINFSQ
         IF (ICP == 2) CPN(NN2,J,TT)=UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))
         CPNN(NN2,J,TT)=PI*PI*(UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                   -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0)))
!-----------------------------------------------------------------------------------------------!
!    Without Rotational Velocity                                                                !
!-----------------------------------------------------------------------------------------------!
         IF (IP == 0) THEN
            UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ))**2+(VWZ/UU(JJ))**2
            IF (ICP == 1) CPN(NN2,J,TT)=1.D0-VTTSQ/UINFSQ-2.D0*DPHIDT*FRELAX/UINFSQ &
                             -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))/UINFSQ
            IF (ICP == 2) CPN(NN2,J,TT)=UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0))
            CPNN(NN2,J,TT)=PI*PI*(UINFSQ-VTTSQ-2.D0*DPHIDT*FRELAX &
                                   -SDIV(2.D0,FN*FN,0.D0)*RN0*DCOS(TN0-DTETA*(DFLOAT(TT)-1.D0)))
         END IF !(IP == 0)
      END IF !(IROTOR)
   END DO !J=1,NNTP
END IF !(TT == 0)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PRESN
!-----------------------------------------------------------------------------------------------!
