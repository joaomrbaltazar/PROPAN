!-----------------------------------------------------------------------------------------------!
!    Dynamic boundary condition for cavity flow                                                 !
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
SUBROUTINE CAVPOTS(JJ,TT,CC)
!-----------------------------------------------------------------------------------------------!
!    Created by: 27102014, J. Baltazar, Cavitation Model                                        !
!    Modified  : 17112014, J. Baltazar, version 3.1, Unsteady Cavitation Model                  !
!    Modified  : 20112014, J. Baltazar, version 3.3, Super-Cavitation                           !
!    Modified  : 10122014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!    Modified  : 04032015, J. Baltazar, correction for super-cavitation                         !
!    Modified  : 09092015, J. Baltazar, correction in the gravity term                          !
!    Modified  : 11042016, J. Baltazar, 2016 version 1.1                                        !
!    Modified  : 16052016, J. Baltazar, 2016 version 1.2                                        !
!    Modified  : 21062016, J. Baltazar, 2016 version 1.3, new relaxation for unsteady terms     !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,J2,JJ,TT,CC
DOUBLE PRECISION :: RP0,TP0,RPW0,TPW0
DOUBLE PRECISION :: A0,A1,A2,S1,S2,S3,POT1,POT2
DOUBLE PRECISION :: VWX,VWY,VWZ,UINFSQ,DPHIDT,FRELAX,SDIV
!-----------------------------------------------------------------------------------------------!
DO J=1,NRP
   IF (IDS(J,TT) /= 0) THEN
!-----------------------------------------------------------------------------------------------!
!    Potential at the detachment                                                                !
!-----------------------------------------------------------------------------------------------!
      S1=DSQRT(AT1XP(IDS(J,TT)-1,J)**2+AT1YP(IDS(J,TT)-1,J)**2+AT1ZP(IDS(J,TT)-1,J)**2)
      S2=DSQRT(AT1XP(IDS(J,TT)-2,J)**2+AT1YP(IDS(J,TT)-2,J)**2+AT1ZP(IDS(J,TT)-2,J)**2)+S1*2.D0
      S3=DSQRT(AT1XP(IDS(J,TT)-2,J)**2+AT1YP(IDS(J,TT)-2,J)**2+AT1ZP(IDS(J,TT)-2,J)**2)+S2
      S3=DSQRT(AT1XP(IDS(J,TT)-3,J)**2+AT1YP(IDS(J,TT)-3,J)**2+AT1ZP(IDS(J,TT)-3,J)**2)+S3
!-----------------------------------------------------------------------------------------------!
      A0=S2*S3/(S1-S2)/(S1-S3)
      A1=S1*S3/(S2-S1)/(S2-S3)
      A2=S1*S2/(S3-S1)/(S3-S2)
      POT1=A0*POTP(IDS(J,TT)-1,J,TT)+A1*POTP(IDS(J,TT)-2,J,TT)+A2*POTP(IDS(J,TT)-3,J,TT)
!-----------------------------------------------------------------------------------------------!
!    Potential at Cavity from the Dynamic Boundary Condition                                    !
!-----------------------------------------------------------------------------------------------!
      DO I=IDS(J,TT),IRS(J,TT)
         IF (TT == 0) THEN
            CALL VWAKE(TT,XP0(I,J),YP0(I,J),ZP0(I,J),1,0,VWX,VWY,VWZ)
         ELSE !(TT == 0)
            CALL VWAKE(TT,XP0(I,J),YP0(I,J),ZP0(I,J),1,IFREQ,VWX,VWY,VWZ)
         END IF !(TT == 0)
         IF (IROTOR == 0) THEN
            UINFSQ=(VWX*UU(JJ)/PI)**2+(VWY*UU(JJ)/PI-ZP0(I,J))**2+(VWZ*UU(JJ)/PI+YP0(I,J))**2
         ELSEIF (IROTOR == 1) THEN
            UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ)-ZP0(I,J))**2+(VWZ/UU(JJ)+YP0(I,J))**2
         END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
         DPHIDT=0.D0
         FRELAX=0.D0
         IF (TT > 0) THEN
            IF (TT == 1) THEN
               DPHIDT=POTP(I,J,1)/DTETA-POTP(I,J,0)/DTETA
            ELSE !(TT == 1)
               DPHIDT=POTP(I,J,TT-2)*0.5D0/DTETA- &
                      POTP(I,J,TT-1)*2.0D0/DTETA+ &
                      POTP(I,J,TT  )*1.5D0/DTETA
            END IF !(TT == 1)
            IF (TT <= (2*NTETA)) THEN
               FRELAX=CREX*DFLOAT(TT-1)/DFLOAT(2*NTETA-1)
            ELSE !(TT <= (2*NTETA))
               FRELAX=CREX
            END IF !(TT <= (2*NTETA))
         END IF !(TT > 0)
!-----------------------------------------------------------------------------------------------!
         RP0=DSQRT(YP0(I,J)*YP0(I,J)+ZP0(I,J)*ZP0(I,J))
         TP0=DATAN2(ZP0(I,J),YP0(I,J))
         S1 =2.D0*DSQRT(AT1XP(I,J)**2+AT1YP(I,J)**2+AT1ZP(I,J)**2)
         IF (TT == 0) THEN
            WORK1=UINFSQ+SIGMA/PI/PI-SDIV(2.D0,FN*FN,0.D0)*YP0(I,J)-VTT2P(I,J)**2
         ELSE !(TT == 0)
            WORK1=UINFSQ+SIGMA/PI/PI-2.D0*DPHIDT*FRELAX-VTT2P(I,J)**2 &
                                    -SDIV(2.D0,FN*FN,0.D0)*RP0*DCOS(TP0-DTETA*(DFLOAT(TT)-1.D0))
         END IF !(TT == 0)
         IF (WORK1 < 0.D0) THEN
            IF (TT == 0) THEN
               WORK1=UINFSQ+SIGMA/PI/PI-SDIV(2.D0,FN*FN,0.D0)*YP0(I,J)
            ELSE !(TT == 0)
               WORK1=UINFSQ+SIGMA/PI/PI-2.D0*DPHIDT*FRELAX &
                                    -SDIV(2.D0,FN*FN,0.D0)*RP0*DCOS(TP0-DTETA*(DFLOAT(TT)-1.D0))
            END IF !(TT == 0)
            WRITE(30,'(A,I4,A,I4)') ' CC=',CC,', TT=',TT
            WRITE(30,'(A,I4,I4)') ' VU2=0 at the D.B.C., (I,J)=',I,J
            if (work1 < 0.d0) then
               write(30,'(a,e23.16)') ' WORK1=',work1
               work1=dabs(work1)
            end if !(work1 < 0.d0)
         END IF !(WORK1 < 0.D0)
         WORK1=DSQRT(WORK1)*S1
         IF (IROTOR == 0) THEN
            WORK2=(VWX*UU(JJ)/PI         )*ET1XP(I,J)+ &
                  (VWY*UU(JJ)/PI-ZP0(I,J))*ET1YP(I,J)+ &
                  (VWZ*UU(JJ)/PI+YP0(I,J))*ET1ZP(I,J)
         ELSEIF (IROTOR == 1) THEN
            WORK2=(VWX/UU(JJ)         )*ET1XP(I,J)+ &
                  (VWY/UU(JJ)-ZP0(I,J))*ET1YP(I,J)+ &
                  (VWZ/UU(JJ)+YP0(I,J))*ET1ZP(I,J)
         END IF !(IROTOR)
         WORK2=WORK2*S1
!-----------------------------------------------------------------------------------------------!
         POT2        =POT1+WORK1-WORK2
         POTP(I,J,TT)=(POT1+POT2)*0.5D0
         POT1        =POT2
      END DO !I=IDS(J,TT),IRS(J,TT)
   END IF !(IDS(J,TT) /= 0)
!-----------------------------------------------------------------------------------------------!
!    Super-Cavitation                                                                           !
!-----------------------------------------------------------------------------------------------!
   J2=J-JI+1
   IF (J2 > 0) THEN
      IF (IDPWS(J2,TT) /= 0) THEN
         DO I=IDPWS(J2,TT),IRPWS(J2,TT)
            IF (TT == 0) THEN
               CALL VWAKE(TT,XPW0(I,J2),YPW0(I,J2),ZPW0(I,J2),1,0,VWX,VWY,VWZ)
            ELSE !(TT == 0)
               CALL VWAKE(TT,XPW0(I,J2),YPW0(I,J2),ZPW0(I,J2),1,IFREQ,VWX,VWY,VWZ)
            END IF !(TT == 0)
            IF (IROTOR == 0) THEN
               UINFSQ=(VWX*UU(JJ)/PI)**2+(VWY*UU(JJ)/PI-ZPW0(I,J2))**2+ &
                                                                   (VWZ*UU(JJ)/PI+YPW0(I,J2))**2
            ELSEIF (IROTOR == 1) THEN
               UINFSQ=(VWX/UU(JJ))**2+(VWY/UU(JJ)-ZPW0(I,J2))**2+(VWZ/UU(JJ)+YPW0(I,J2))**2
            END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
            DPHIDT=0.D0
            FRELAX=0.D0
            IF (TT > 0) THEN
               IF (TT == 1) THEN
                  DPHIDT=POTPWS(I,J2,1)/DTETA-POTPWS(I,J2,0)/DTETA
               ELSE !(TT == 1)
                  DPHIDT=POTPWS(I,J2,TT-2)*0.5D0/DTETA- &
                         POTPWS(I,J2,TT-1)*2.0D0/DTETA+ &
                         POTPWS(I,J2,TT  )*1.5D0/DTETA
               END IF !(TT == 1)
               IF (TT <= (2*NTETA)) THEN
                  FRELAX=CREX*DFLOAT(TT-1)/DFLOAT(2*NTETA-1)
               ELSE !(TT <= (2*NTETA))
                  FRELAX=CREX
               END IF !(TT <= (2*NTETA))
            END IF !(TT > 0)
!-----------------------------------------------------------------------------------------------!
            RPW0=DSQRT(YPW0(I,J2)*YPW0(I,J2)+ZPW0(I,J2)*ZPW0(I,J2))
            TPW0=DATAN2(ZPW0(I,J2),YPW0(I,J2))
            S1  =2.D0*DSQRT(AT1XPW(I,J2)**2+AT1YPW(I,J2)**2+AT1ZPW(I,J2)**2)
            IF (TT == 0) THEN
               WORK1=UINFSQ+SIGMA/PI/PI-SDIV(2.D0,FN*FN,0.D0)*YPW0(I,J2)
            ELSE !(TT == 0)
               WORK1=UINFSQ+SIGMA/PI/PI-2.D0*DPHIDT*FRELAX &
                                  -SDIV(2.D0,FN*FN,0.D0)*RPW0*DCOS(TPW0-DTETA*(DFLOAT(TT)-1.D0))
            END IF !(TT == 0)
            if (work1 < 0.d0) then
               write(30,'(a,e23.16)') ' WORK1=',work1
               work1=dabs(work1)
            end if !(work1 < 0.d0)
            WORK1=DSQRT(WORK1)*S1
            IF (IROTOR == 0) THEN
               WORK2=(VWX*UU(JJ)/PI           )*ET1XPW(I,J2)+ &
                     (VWY*UU(JJ)/PI-ZPW0(I,J2))*ET1YPW(I,J2)+ &
                     (VWZ*UU(JJ)/PI+YPW0(I,J2))*ET1ZPW(I,J2)
            ELSEIF (IROTOR == 1) THEN
               WORK2=(VWX/UU(JJ)           )*ET1XPW(I,J2)+ &
                     (VWY/UU(JJ)-ZPW0(I,J2))*ET1YPW(I,J2)+ &
                     (VWZ/UU(JJ)+YPW0(I,J2))*ET1ZPW(I,J2)
            END IF !(IROTOR)
            WORK2=WORK2*S1
!-----------------------------------------------------------------------------------------------!
            POT2           =POT1+WORK1-WORK2
            POTPWS(I,J2,TT)=(POT1+POT2)*0.5D0
            POT1           =POT2
         END DO !I=IDPWS(J2,TT),IRPWS(J2,TT)
      END IF !(J2 > 0)
   END IF !(IDPWS(J2,TT) /= 0)
END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE CAVPOTS
!-----------------------------------------------------------------------------------------------!
