!-----------------------------------------------------------------------------------------------!
!    Calculation of the cavity properties                                                       !
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
SUBROUTINE CAVPROP(TT)
!-----------------------------------------------------------------------------------------------!
!    Created by: 07112014, J. Baltazar, Cavitation Model                                        !
!    Modified  : 12112014, J. Baltazar, Cavitation on suction and pressure sides                !
!    Modified  : 27112014, J. Baltazar, version 3.3, Super-Cavitation Model                     !
!    Modified  : 11042016, J. Baltazar, 2016 version 1.1                                        !
!    Modified  : 24052016, J. Baltazar, 2016 version 1.2                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,J2,TT
DOUBLE PRECISION :: LL,LC,AES,ACS,VCS,AEP,ACP,VCP
!-----------------------------------------------------------------------------------------------!
!    Cavitation on the Suction Side of the Blade                                                !
!-----------------------------------------------------------------------------------------------!
LCMAXS=0.D0
AES   =0.D0
ACS   =0.D0
VCS   =0.D0
RACAES=0.D0
ETAS  =0.D0
IF (NSCAV /= 0) THEN
   DO J=1,NRP
      LL=0.D0
      LC=0.D0
      DO I=NC1,NCP
         LL=LL+2.D0*DSQRT(AT1XP(I,J)**2+AT1YP(I,J)**2+AT1ZP(I,J)**2)
         AES=AES+AP0(I,J)
         IF ((I >= IDS(J,TT)).AND.(I <= IRS(J,TT))) THEN
            LC=LC+2.D0*DSQRT(AT1XP(I,J)**2+AT1YP(I,J)**2+AT1ZP(I,J)**2)
            ACS=ACS+AP0(I,J)
            VCS=VCS+AP0(I,J)*THICKP(I,J,TT)
         END IF !((I >= IDS(J,TT)).AND.(I <= IRS(J,TT)))
      END DO !I=NC1,NCP
!-----------------------------------------------------------------------------------------------!
!    Cavitation on the Blade Wake                                                               !
!-----------------------------------------------------------------------------------------------!
      IF (IDPWS(J,TT) /= 0) THEN
         J2=J-JI+1
         DO I=1,IABS(NCPW)
            LL=LL+2.D0*DSQRT(AT1XPW(I,J2)**2+AT1YPW(I,J2)**2+AT1ZPW(I,J2)**2)
            AEP=AEP+APW0(I,J2)
            IF ((I >= IDPWS(J2,TT)).AND.(I <= IRPWS(J2,TT))) THEN
               LC=LC+2.D0*DSQRT(AT1XPW(I,J2)**2+AT1YPW(I,J2)**2+AT1ZPW(I,J2)**2)
               ACP=ACP+APW0(I,J2)
               VCP=VCP+APW0(I,J2)*THICKPWS(I,J2,TT)
            END IF !((I >= IDPWS(J2,TT)).AND.(I <= IRPWS(J2,TT)))
         END DO !I=1,IABS(NCPW)
      END IF !(IDPWS(J,TT) /= 0)
!-----------------------------------------------------------------------------------------------!
      LC=LC/LL
      IF (LC > LCMAXS) LCMAXS=LC
   END DO !J=1,NRP
   RACAES=ACS/AES
   ETAS  =VCS/ACS/LCMAXS*10.D0
END IF !(NSCAV /= 0)
!-----------------------------------------------------------------------------------------------!
!    Cavitation on the Pressure Side of the Blade                                               !
!-----------------------------------------------------------------------------------------------!
LCMAXP=0.D0
AEP   =0.D0
ACP   =0.D0
VCP   =0.D0
RACAEP=0.D0
ETAP  =0.D0
IF (NPCAV /= 0) THEN
   DO J=1,NRP
      LL=0.D0
      LC=0.D0
      DO I=NC,1,-1
         LL=LL+2.D0*DSQRT(AT1XP(I,J)**2+AT1YP(I,J)**2+AT1ZP(I,J)**2)
         AEP=AEP+AP0(I,J)
         IF ((I <= IDP(J,TT)).AND.(I >= IRP(J,TT))) THEN
            LC=LC+2.D0*DSQRT(AT1XP(I,J)**2+AT1YP(I,J)**2+AT1ZP(I,J)**2)
            ACP=ACP+AP0(I,J)
            VCP=VCP+AP0(I,J)*THICKP(I,J,TT)
         END IF !((I <= IDP(J,TT)).AND.(I >= IRP(J,TT)))
      END DO !I=NC,1,-1
!-----------------------------------------------------------------------------------------------!
!    Cavitation on the Blade Wake                                                               !
!-----------------------------------------------------------------------------------------------!
      IF (IDPWP(J,TT) /= 0) THEN
         J2=J-JI+1
         DO I=1,IABS(NCPW)
            LL=LL+2.D0*DSQRT(AT1XPW(I,J2)**2+AT1YPW(I,J2)**2+AT1ZPW(I,J2)**2)
            AEP=AEP+APW0(I,J2)
            IF ((I >= IDPWP(J2,TT)).AND.(I <= IRPWP(J2,TT))) THEN
               LC=LC+2.D0*DSQRT(AT1XPW(I,J2)**2+AT1YPW(I,J2)**2+AT1ZPW(I,J2)**2)
               ACP=ACP+APW0(I,J2)
               VCP=VCP+APW0(I,J2)*THICKPWP(I,J2,TT)
            END IF !((I >= IDPWP(J2,TT)).AND.(I <= IRPWP(J2,TT)))
         END DO !I=1,IABS(NCPW)
      END IF !(IDPWP(J,TT) /= 0)
!-----------------------------------------------------------------------------------------------!
      LC=LC/LL
      IF (LC > LCMAXP) LCMAXP=LC
   END DO !J=1,NRP
   RACAEP=ACP/AEP
   ETAP  =VCP/ACP/LCMAXP*10.D0
END IF !(NPCAV /=0)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE CAVPROP
!-----------------------------------------------------------------------------------------------!
