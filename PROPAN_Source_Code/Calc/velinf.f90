!-----------------------------------------------------------------------------------------------!
!    Undisturbed Inflow Velocity in Inertial or Rotating Reference Frames                       !
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
SUBROUTINE VELINF(CREF,JJ,TT,N1,N2,XTMP,YTMP,ZTMP,UTMP,VTMP,WTMP)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar                                                                    !
!    Modified  : 02122013, J. Baltazar, version 1.0                                             !
!    Modified  : 27052014, J. Baltazar, Wake Alignment Module                                   !
!    Modified  : 05022016, J. Baltazar, rotating and inertial reference frames                  !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
CHARACTER*8 CREF
INTEGER :: I,J,J2,JJ,TT,N1,N2
DOUBLE PRECISION :: VWX,VWY,VWZ,VWT,VWR,VXX(NWR),VTT(NWR),VRR(NWR)
DOUBLE PRECISION :: XTMP(N1,N2),YTMP(N1,N2),ZTMP(N1,N2),RTMP(N1,N2),TTMP(N1,N2)
DOUBLE PRECISION :: UTMP(N1,N2),VTMP(N1,N2),WTMP(N1,N2)
!-----------------------------------------------------------------------------------------------!
RTMP=DSQRT(YTMP*YTMP+ZTMP*ZTMP)
TTMP=DATAN2(ZTMP,YTMP)
!-----------------------------------------------------------------------------------------------!
UTMP=0.D0
VTMP=0.D0
WTMP=0.D0
!-----------------------------------------------------------------------------------------------!
!    Undisturbed Inflow Velocity in Rotating Reference Frame                                    !
!-----------------------------------------------------------------------------------------------!
IF (CREF == 'ROTATING') THEN
!-----------------------------------------------------------------------------------------------!
!    Loop on Field Points                                                                       !
!-----------------------------------------------------------------------------------------------!
   DO J=1,N2
      DO I=1,N1
         IF (TT == 0) THEN
            CALL VWAKE(TT,XTMP(I,J),YTMP(I,J),ZTMP(I,J),1,    0,VWX,VWY,VWZ)
         ELSE !(TT == 0)
            CALL VWAKE(TT,XTMP(I,J),YTMP(I,J),ZTMP(I,J),1,IFREQ,VWX,VWY,VWZ)
         END IF !(TT == 0)
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
         IF (IROTOR == 0) THEN
            UTMP(I,J)=VWX*UU(JJ)/PI
            VTMP(I,J)=VWY*UU(JJ)/PI-ZTMP(I,J)
            WTMP(I,J)=VWZ*UU(JJ)/PI+YTMP(I,J)
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
         ELSEIF (IROTOR == 1) THEN
            UTMP(I,J)=VWX/UU(JJ)
            VTMP(I,J)=VWY/UU(JJ)-ZTMP(I,J)
            WTMP(I,J)=VWZ/UU(JJ)+YTMP(I,J)
         END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
      END DO !I=1,N1
   END DO !J=1,N2
END IF !(CREF == 'ROTATING')
!-----------------------------------------------------------------------------------------------!
!    Undisturbed Inflow Velocity in Inertial Reference Frame                                    !
!-----------------------------------------------------------------------------------------------!
IF (CREF == 'INERTIAL') THEN
!-----------------------------------------------------------------------------------------------!
!    Loop on Field Points                                                                       !
!-----------------------------------------------------------------------------------------------!
   DO J=1,N2
      DO I=1,N1
         IF (TT == 0) THEN
!-----------------------------------------------------------------------------------------------!
!    Axial Velocity                                                                             !
!-----------------------------------------------------------------------------------------------!
            DO J2=1,NWR
               CALL FOURIER_FUNCTION(0,ANVX(:,J2),BNVX(:,J2),TTMP(I,J),VXX(J2))
            END DO !J2=1,NWR
            CALL LININT(NWR,WRR,VXX,1,RTMP(I,J),VWX)
!-----------------------------------------------------------------------------------------------!
!    Tangential Velocity                                                                        !
!-----------------------------------------------------------------------------------------------!
            DO J2=1,NWR
               CALL FOURIER_FUNCTION(0,ANVT(:,J2),BNVT(:,J2),TTMP(I,J),VTT(J2))
            END DO !J2=1,NWR
            CALL LININT(NWR,WRR,VTT,1,RTMP(I,J),VWT)
!-----------------------------------------------------------------------------------------------!
!    Radial Velocity                                                                            !
!-----------------------------------------------------------------------------------------------!
            DO J2=1,NWR
               CALL FOURIER_FUNCTION(0,ANVR(:,J2),BNVR(:,J2),TTMP(I,J),VRR(J2))
            END DO !J2=1,NWR
            CALL LININT(NWR,WRR,VRR,1,RTMP(I,J),VWR)
         ELSE !(TT == 0)
!-----------------------------------------------------------------------------------------------!
!    Axial Velocity                                                                             !
!-----------------------------------------------------------------------------------------------!
            DO J2=1,NWR
               CALL FOURIER_FUNCTION(IFREQ,ANVX(:,J2),BNVX(:,J2),TTMP(I,J),VXX(J2))
            END DO !J2=1,NWR
            CALL LININT(NWR,WRR,VXX,1,RTMP(I,J),VWX)
!-----------------------------------------------------------------------------------------------!
!    Tangential Velocity                                                                        !
!-----------------------------------------------------------------------------------------------!
            DO J2=1,NWR
               CALL FOURIER_FUNCTION(IFREQ,ANVT(:,J2),BNVT(:,J2),TTMP(I,J),VTT(J2))
            END DO !J2=1,NWR
            CALL LININT(NWR,WRR,VTT,1,RTMP(I,J),VWT)
!-----------------------------------------------------------------------------------------------!
!    Radial Velocity                                                                            !
!-----------------------------------------------------------------------------------------------!
            DO J2=1,NWR
               CALL FOURIER_FUNCTION(IFREQ,ANVR(:,J2),BNVR(:,J2),TTMP(I,J),VRR(J2))
            END DO !J2=1,NWR
            CALL LININT(NWR,WRR,VRR,1,RTMP(I,J),VWR)
         END IF !(TT == 0)
!-----------------------------------------------------------------------------------------------!
!    Cartesian Components                                                                       !
!-----------------------------------------------------------------------------------------------!
         VWY=VWR*DCOS(TTMP(I,J))-VWT*DSIN(TTMP(I,J))
         VWZ=VWR*DSIN(TTMP(I,J))+VWT*DCOS(TTMP(I,J))
!-----------------------------------------------------------------------------------------------!
!    Propeller Case                                                                             !
!-----------------------------------------------------------------------------------------------!
         IF (IROTOR == 0) THEN
            UTMP(I,J)=VWX*UU(JJ)/PI
            VTMP(I,J)=VWY*UU(JJ)/PI
            WTMP(I,J)=VWZ*UU(JJ)/PI
!-----------------------------------------------------------------------------------------------!
!    Turbine Case                                                                               !
!-----------------------------------------------------------------------------------------------!
         ELSEIF (IROTOR == 1) THEN
            UTMP(I,J)=VWX/UU(JJ)
            VTMP(I,J)=VWY/UU(JJ)
            WTMP(I,J)=VWZ/UU(JJ)
         END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
      END DO !I=1,N1
   END DO !J=1,N2
END IF !(CREF == 'INERTIAL')
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE VELINF
!-----------------------------------------------------------------------------------------------!
