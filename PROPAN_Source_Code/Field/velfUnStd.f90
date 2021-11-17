!-----------------------------------------------------------------------------------------------!
!    Calculation of the Unsteady Velocity Field                                                 !
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
SUBROUTINE VELFUNSTD(CREF,JJ,TT,N1,N2,XTMP,YTMP,ZTMP,UTMP,VTMP,WTMP)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar                                                                    !
!    Modified  : 03122013, J. Baltazar, version 1.0                                             !
!    Modified  : 28052014, J. Baltazar, revised                                                 !
!    Modified  : 10022016, J. Baltazar, rotating and inertial reference frames                  !
!    Modified  : 28042016, J. Baltazar, version 1.1                                             !
!    Modified  : 09052016, J. Baltazar, version 1.2                                             !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
CHARACTER*8 CREF
INTEGER :: I,J,JJ,TT,T1,N1,N2
DOUBLE PRECISION :: TETA
DOUBLE PRECISION :: XTMP(N1,N2),YTMP(N1,N2),ZTMP(N1,N2)
DOUBLE PRECISION :: UTMP(N1,N2),VTMP(N1,N2),WTMP(N1,N2)
DOUBLE PRECISION :: XNEW(N1,N2),YNEW(N1,N2),ZNEW(N1,N2),RNEW(N1,N2),TNEW(N1,N2)
DOUBLE PRECISION :: UNEW(N1,N2),VNEW(N1,N2),WNEW(N1,N2)
!-----------------------------------------------------------------------------------------------!
XNEW=0.D0
YNEW=0.D0
ZNEW=0.D0
RNEW=0.D0
TNEW=0.D0
IF (CREF == 'INERTIAL') THEN
   DO J=1,N2
      DO I=1,N1
         CALL FRAME('ROTATING',TT,1,XTMP(I,J),YTMP(I,J),ZTMP(I,J), &
                                              XNEW(I,J),YNEW(I,J),ZNEW(I,J),RNEW(I,J),TNEW(I,J))
      END DO !I=1,N1
   END DO !J=1,N2
END IF !(CREF == 'INERTIAL')
!-----------------------------------------------------------------------------------------------!
!    Undisturbed Velocity                                                                       !
!-----------------------------------------------------------------------------------------------!
CALL VELINF(CREF,JJ,TT,N1,N2,XTMP,YTMP,ZTMP,UNEW,VNEW,WNEW)
UTMP=UNEW
VTMP=VNEW
WTMP=WNEW
!-----------------------------------------------------------------------------------------------!
!    Blade Velocities                                                                           !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   IF (CREF == 'ROTATING') THEN
      CALL BLADEVELO(TT,N1,N2,XTMP,YTMP,ZTMP,UNEW,VNEW,WNEW)
      UTMP=UTMP+UNEW
      VTMP=VTMP+VNEW
      WTMP=WTMP+WNEW
   ELSEIF (CREF == 'INERTIAL') THEN
      CALL BLADEVELO(TT,N1,N2,XNEW,YNEW,ZNEW,UNEW,VNEW,WNEW)
      T1=TT
      teta=0.d0
      IF (NTETA /= 0) THEN
         T1=TT-(TT-1)/NTETA*DFLOAT(NTETA)
         TETA=2.D0*PI-DTETA*(DFLOAT(T1)-1.D0)
      END IF !(NTETA /= 0)
      UTMP=UTMP+UNEW
      VTMP=VTMP+VNEW*DCOS(TETA)-WNEW*DSIN(TETA)
      WTMP=WTMP+VNEW*DSIN(TETA)+WNEW*DCOS(TETA)
   END IF !(CREF)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake Velocities                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   IF (CREF == 'ROTATING') THEN
      CALL BLADEWAKEVELO(TT,N1,N2,XTMP,YTMP,ZTMP,UNEW,VNEW,WNEW)
      UTMP=UTMP+UNEW
      VTMP=VTMP+VNEW
      WTMP=WTMP+WNEW
   ELSEIF (CREF == 'INERTIAL') THEN
      CALL BLADEWAKEVELO(TT,N1,N2,XNEW,YNEW,ZNEW,UNEW,VNEW,WNEW)
      T1=TT
      teta=0.d0
      IF (NTETA /= 0) THEN
         T1=TT-(TT-1)/NTETA*DFLOAT(NTETA)
         TETA=2.D0*PI-DTETA*(DFLOAT(T1)-1.D0)
      END IF !(NTETA /= 0)
      UTMP=UTMP+UNEW
      VTMP=VTMP+VNEW*DCOS(TETA)-WNEW*DSIN(TETA)
      WTMP=WTMP+VNEW*DSIN(TETA)+WNEW*DCOS(TETA)
   END IF !(CREF)
END IF !(IP ==1 )
!-----------------------------------------------------------------------------------------------!
!    Nozzle Velocities                                                                          !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   IF (CREF == 'ROTATING') THEN
      CALL NOZZLEVELO(TT,N1,N2,XTMP,YTMP,ZTMP,UNEW,VNEW,WNEW)
      UTMP=UTMP+UNEW
      VTMP=VTMP+VNEW
      WTMP=WTMP+WNEW
   ELSEIF (CREF == 'INERTIAL') THEN
      CALL NOZZLEVELO(TT,N1,N2,XNEW,YNEW,ZNEW,UNEW,VNEW,WNEW)
      T1=TT
      teta=0.d0
      IF (NTETA /= 0) THEN
         T1=TT-(TT-1)/NTETA*DFLOAT(NTETA)
         TETA=2.D0*PI-DTETA*(DFLOAT(T1)-1.D0)
      END IF !(NTETA /= 0)
      UTMP=UTMP+UNEW
      VTMP=VTMP+VNEW*DCOS(TETA)-WNEW*DSIN(TETA)
      WTMP=WTMP+VNEW*DSIN(TETA)+WNEW*DCOS(TETA)
   END IF !(CREF)
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake Velocities                                                                     !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   IF (CREF == 'ROTATING') THEN
      CALL NOZZLEWAKEVELO(TT,N1,N2,XTMP,YTMP,ZTMP,UNEW,VNEW,WNEW)
      UTMP=UTMP+UNEW
      VTMP=VTMP+VNEW
      WTMP=WTMP+WNEW
   ELSEIF (CREF == 'INERTIAL') THEN
      CALL NOZZLEWAKEVELO(TT,N1,N2,XNEW,YNEW,ZNEW,UNEW,VNEW,WNEW)
      T1=TT
      IF (NTETA /= 0) THEN
         T1=TT-(TT-1)/NTETA*DFLOAT(NTETA)
         TETA=2.D0*PI-DTETA*(DFLOAT(T1)-1.D0)
      END IF !(NTETA /= 0)
      UTMP=UTMP+UNEW
      VTMP=VTMP+VNEW*DCOS(TETA)-WNEW*DSIN(TETA)
      WTMP=WTMP+VNEW*DSIN(TETA)+WNEW*DCOS(TETA)
   END IF !(CREF)
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub Velocities                                                                             !
!-----------------------------------------------------------------------------------------------!
IF (IABS(IH) == 1) THEN
   IF (CREF == 'ROTATING') THEN
      IF (NWA == 0) CALL      HUBVELO(TT,N1,N2,XTMP,YTMP,ZTMP,UNEW,VNEW,WNEW)
      IF (NWA /= 0) CALL IMAGEHUBVELO(TT,N1,N2,XTMP,YTMP,ZTMP,UNEW,VNEW,WNEW)
      UTMP=UTMP+UNEW
      VTMP=VTMP+VNEW
      WTMP=WTMP+WNEW
   ELSEIF (CREF == 'INERTIAL') THEN
      IF (NWA == 0) CALL      HUBVELO(TT,N1,N2,XNEW,YNEW,ZNEW,UNEW,VNEW,WNEW)
      IF (NWA /= 0) CALL IMAGEHUBVELO(TT,N1,N2,XNEW,YNEW,ZNEW,UNEW,VNEW,WNEW)
      T1=TT
      teta=0.d0
      IF (NTETA /= 0) THEN
         T1=TT-(TT-1)/NTETA*DFLOAT(NTETA)
         TETA=2.D0*PI-DTETA*(DFLOAT(T1)-1.D0)
      END IF !(NTETA /= 0)
      UTMP=UTMP+UNEW
      VTMP=VTMP+VNEW*DCOS(TETA)-WNEW*DSIN(TETA)
      WTMP=WTMP+VNEW*DSIN(TETA)+WNEW*DCOS(TETA)
   END IF !(CREF)
END IF !(IABS(IH) == 1)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE VELFUNSTD
!-----------------------------------------------------------------------------------------------!
