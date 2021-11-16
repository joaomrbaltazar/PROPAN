!-----------------------------------------------------------------------------------------------!
!    Cavitation Check                                                                           !
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
SUBROUTINE CAVCHECK(TT)
!-----------------------------------------------------------------------------------------------!
!    Created by: 07112014, J. Baltazar, version 3.0, Steady Cavitation Model                    !
!    Modified  : 13112014, J. Baltazar, version 3.1, Unsteady Cavitation Model                  !
!    Modified  : 27112014, J. Baltazar, version 3.3, Super-Cavitation Model                     !
!    Modified  : 11122014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!    Modified  : 22042016, J. Baltazar, 2016 version 1.1                                        !
!    Modified  : 31052016, J. Baltazar, 2016 version 1.2                                        !
!    Modified  : 22062016, J. Baltazar, 2016 version 1.3                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,TT
!-----------------------------------------------------------------------------------------------!
!    Cavitation Check                                                                           !
!-----------------------------------------------------------------------------------------------!
NPCAV=0
NSCAV=0
NWCAV=0
!-----------------------------------------------------------------------------------------------!
!    Cavity on the Pressure Side                                                                !
!-----------------------------------------------------------------------------------------------!
IF (NCPW <= 0) THEN
   DO J=1,(NRP-1) !Last Strip Excluded
      DO I=NC,1,-1
         IF ((CPNP(I,J,TT)+SIGMA) <= TOL) NPCAV=NPCAV+1
      END DO !I=NC,1,-1
   END DO !J=1,(NRP-1)
END IF !(NCPW <= 0)
!-----------------------------------------------------------------------------------------------!
IF (NPCAV /= 0) THEN
   DO J=1,(NRP-1) !Last Strip Excluded
      I=NC
      DO WHILE (((CPNP(I,J,TT)+SIGMA) > TOL).AND.(I > 1))
         I=I-1
      END DO !(((CPNP(I,J,TT)+SIGMA) > TOL).AND.(I > 1))
      IDP(J,TT)=I
      IF ((IDP(J,TT) == 1).AND.((CPNP(1,J,TT)+SIGMA) > TOL)) IDP(J,TT)=0
!-----------------------------------------------------------------------------------------------!
      IRP(J,TT)=0
      IF (IDP(J,TT) /= 0) THEN
         I=IDP(J,TT)
         DO WHILE (((CPNP(I,J,TT)+SIGMA) <= TOL).AND.(I > 1))
            I=I-1
         END DO !(((CPNP(I,J,TT)+SIGMA) <= TOL).AND.(I > 1))
         IRP(J,TT)=I+1
      END IF !(IDP(J,TT) /= 0)
   END DO !J=1,(NRP-1)
!-----------------------------------------------------------------------------------------------!
!    Number of Cavitating Panels                                                                !
!-----------------------------------------------------------------------------------------------!
   NPCAV=0
   DO J=1,(NRP-1) !Last Strip Excluded
      DO I=IDP(J,TT),IRP(J,TT),-1
         IF (I /= 0) NPCAV=NPCAV+1
      END DO !I=IDP(J,TT),IRP(J,TT),-1
   END DO !J=1,(NRP-1) !Last Strip Excluded
END IF !(NPCAV /= 0)
!-----------------------------------------------------------------------------------------------!
!    Cavity on the Suction Side                                                                 !
!-----------------------------------------------------------------------------------------------!
IF (NCPW >= 0) THEN
   DO J=1,(NRP-1) !Last Strip Excluded
      DO I=NC1,NCP
         IF ((CPNP(I,J,TT)+SIGMA) <= TOL) NSCAV=NSCAV+1
      END DO !I=NC1,NCP
   END DO !J=1,(NRP-1)
END IF !(NCPW >= 0)
!-----------------------------------------------------------------------------------------------!
IF (NSCAV /= 0) THEN
   DO J=1,(NRP-1) !Last Strip Excluded
      I=NC1
      DO WHILE (((CPNP(I,J,TT)+SIGMA) > TOL).AND.(I < NCP))
         I=I+1
      END DO !(((CPNP(I,J,TT)+SIGMA) > TOL).AND.(I < NCP))
      IDS(J,TT)=I
      IF ((IDS(J,TT) == NCP).AND.((CPNP(NCP,J,TT)+SIGMA) > TOL)) IDS(J,TT)=0
!-----------------------------------------------------------------------------------------------!
      IRS(J,TT)=0
      IF (IDS(J,TT) /= 0) THEN
         I=IDS(J,TT)
         DO WHILE (((CPNP(I,J,TT)+SIGMA) <= TOL).AND.(I < NCP))
            I=I+1
         END DO !(((CPNP(I,J,TT)+SIGMA) <= TOL).AND.(I < NCP))
         IRS(J,TT)=I-1
      END IF !(IDS(J,TT) /= 0)
   END DO !J=1,(NRP-1)
!-----------------------------------------------------------------------------------------------!
!    Number of Cavitating Panels                                                                !
!-----------------------------------------------------------------------------------------------!
   NSCAV=0
   DO J=1,(NRP-1) !Last Strip Excluded
      DO I=IDS(J,TT),IRS(J,TT)
         IF (I /= 0) NSCAV=NSCAV+1
      END DO !I=IDS(J,TT),IRS(J,TT)
   END DO !J=1,(NRP-1) !Last Strip Excluded
END IF !(NSCAV /= 0)
!-----------------------------------------------------------------------------------------------!
!    Unsteady Cavitation                                                                        !
!-----------------------------------------------------------------------------------------------!
IF ((TT > 0).AND.((NPCAV /= 0).OR.(NSCAV /= 0))) THEN
!* IDP   (  :,TT)=IDP   (  :,TT-1)
!* IRP   (  :,TT)=IRP   (  :,TT-1)
!* IDS   (  :,TT)=IDS   (  :,TT-1)
!* IRS   (  :,TT)=IRS   (  :,TT-1)
   IF (IRED == 1) POTP(:,:,TT)=POTP(:,:,TT-1)
!! THICKP(:,:,TT)=THICKP(:,:,TT-1)
!-----------------------------------------------------------------------------------------------!
!* IDPWP (  :,TT)=IDPWP (  :,TT-1)
!* IRPWP (  :,TT)=IRPWP (  :,TT-1)
!* IDPWS (  :,TT)=IDPWS (  :,TT-1)
!* IRPWS (  :,TT)=IRPWS (  :,TT-1)
   IF (NCPW < 0) THEN
      IF (IRED == 1) POTPWP     (:,:,TT)=POTPWP     (:,:,TT-1)
      IF (IRED == 1) SOURCEPWCAV(:,:,TT)=SOURCEPWCAV(:,:,TT-1)
!!    THICKPWP (:,:,TT)=THICKPWP (:,:,TT-1)
!!    CAMBERPWP(:,:,TT)=CAMBERPWP(:,:,TT-1)
   END IF !(NCPW < 0)
   IF (NCPW > 0) THEN
      IF (IRED == 1) POTPWS     (:,:,TT)=POTPWS     (:,:,TT-1)
      IF (IRED == 1) SOURCEPWCAV(:,:,TT)=SOURCEPWCAV(:,:,TT-1)
!!    THICKPWS (:,:,TT)=THICKPWS (:,:,TT-1)
!!    CAMBERPWS(:,:,TT)=CAMBERPWS(:,:,TT-1)
   END IF !(NCPW > 0)
END IF !((TT > 0).AND.((NPCAV /= 0).OR.(NSCAV /= 0)))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE CAVCHECK
!-----------------------------------------------------------------------------------------------!
