!-----------------------------------------------------------------------------------------------!
!    Calculation of the Unsteady Pressure Field                                                 !
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
SUBROUTINE PRESFUNSTD(CREF,JJ,TT,N1,N2,XTMP,YTMP,ZTMP,DPOTDT)
!-----------------------------------------------------------------------------------------------!
!    Created by: 11022016, J. Baltazar, version 1.0                                             !
!    Modified  : 22042016, J. Baltazar, version 1.1                                             !
!    Modified  : 25102016, J. Baltazar, version 1.4                                             !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
CHARACTER*8 CREF
INTEGER :: I,J,JJ,TT,N1,N2
DOUBLE PRECISION :: XTMP(N1,N2),YTMP(N1,N2),ZTMP(N1,N2)
DOUBLE PRECISION :: XNEW(N1,N2),YNEW(N1,N2),ZNEW(N1,N2),RNEW(N1,N2),TNEW(N1,N2)
DOUBLE PRECISION :: PTMP(N1,N2),DPOTDT(N1,N2)
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
!    Blade Potential                                                                            !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   IF (CREF == 'ROTATING') THEN
   ELSEIF (CREF == 'INERTIAL') THEN
      CALL BLADECOEFF(TT,N1,N2,XNEW,YNEW,ZNEW,PTMP)
      POTF(:,:,TT)=PTMP(:,:)
   END IF !(CREF)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake Potential                                                                       !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   IF (CREF == 'ROTATING') THEN
   ELSEIF (CREF == 'INERTIAL') THEN
      CALL BLADEWAKECOEFF(TT,N1,N2,XNEW,YNEW,ZNEW,PTMP)
      POTF(:,:,TT)=POTF(:,:,TT)+PTMP(:,:)
   END IF !(CREF)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Potential                                                                           !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   IF (CREF == 'ROTATING') THEN
   ELSEIF (CREF == 'INERTIAL') THEN
      CALL NOZZLECOEFF(TT,N1,N2,XNEW,YNEW,ZNEW,PTMP)
      POTF(:,:,TT)=POTF(:,:,TT)+PTMP(:,:)
   END IF !(CREF)
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake Potential                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   IF (CREF == 'ROTATING') THEN
   ELSEIF (CREF == 'INERTIAL') THEN
      CALL NOZZLEWAKECOEFF(TT,N1,N2,XNEW,YNEW,ZNEW,PTMP)
      POTF(:,:,TT)=POTF(:,:,TT)+PTMP(:,:)
   END IF !(CREF)
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub Potential                                                                              !
!-----------------------------------------------------------------------------------------------!
IF (IABS(IH) == 1) THEN
   IF (CREF == 'ROTATING') THEN
   ELSEIF (CREF == 'INERTIAL') THEN
      CALL HUBCOEFF(TT,N1,N2,XNEW,YNEW,ZNEW,PTMP)
      POTF(:,:,TT)=POTF(:,:,TT)+PTMP(:,:)
   END IF !(CREF)
END IF !(IABS(IH) == 1)
!-----------------------------------------------------------------------------------------------!
IF (TT == 0) THEN
   DPOTDT(:,:)=0.D0
ELSEIF (TT == 1) THEN
   DPOTDT(:,:)=POTF(:,:,1)/DTETA-POTF(:,:,0)/DTETA
ELSE !(TT == 1)
   DPOTDT(:,:)=POTF(:,:,TT-2)*0.5D0/DTETA- &
               POTF(:,:,TT-1)*2.0D0/DTETA+ &
               POTF(:,:,TT  )*1.5D0/DTETA
END IF !(TT == 1)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PRESFUNSTD
!-----------------------------------------------------------------------------------------------!
