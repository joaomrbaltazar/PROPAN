!-----------------------------------------------------------------------------------------------!
!    Write Error in Cavitation Model                                                            !
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
SUBROUTINE CAVERRC(CC,ERRC)
!-----------------------------------------------------------------------------------------------!
!    Created by: 19112014, J. Baltazar, version 3.2                                             !
!    Modified  : 04122014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: CC
DOUBLE PRECISION :: ERRC
!-----------------------------------------------------------------------------------------------!
ERRC=0.D0
IF (CC > 0) THEN
   ERRLCS=DABS(ERRLCS-LCMAXS)/(ERRLCS+1.D-1)
   ERRACS=DABS(ERRACS-RACAES)/(ERRACS+1.D-1)
   ERRTCS=DABS(ERRTCS-ETAS  )/(ERRTCS+1.D-1)
   ERRLCP=DABS(ERRLCP-LCMAXP)/(ERRLCP+1.D-1)
   ERRACP=DABS(ERRACP-RACAEP)/(ERRACP+1.D-1)
   ERRTCP=DABS(ERRTCP-ETAP  )/(ERRTCP+1.D-1)
   IF (IROTOR == 0) THEN
      IF (NSCAV /= 0) THEN
         WRITE(20,*) 'ERRLCS=',ERRLCS
         WRITE(20,*) 'ERRACS=',ERRACS
         WRITE(20,*) 'ERRTCS=',ERRTCS
      END IF !(NSCAV /= 0)
      IF (NPCAV /= 0) THEN
         WRITE(20,*) 'ERRLCP=',ERRLCP
         WRITE(20,*) 'ERRACP=',ERRACP
         WRITE(20,*) 'ERRTCP=',ERRTCP
      END IF !(NPCAV /= 0)
   ELSEIF (IROTOR == 1) THEN
      IF (NSCAV /= 0) THEN
         WRITE(20,*) 'ERRLCP=',ERRLCS
         WRITE(20,*) 'ERRACP=',ERRACS
         WRITE(20,*) 'ERRTCP=',ERRTCS
      END IF !(NSCAV /= 0)
      IF (NPCAV /= 0) THEN
         WRITE(20,*) 'ERRLCS=',ERRLCP
         WRITE(20,*) 'ERRACS=',ERRACP
         WRITE(20,*) 'ERRTCS=',ERRTCP
      END IF !(NPCAV /= 0)
   END IF !(IROTOR)
   ERRC=DMAX1(ERRLCS,ERRACS,ERRTCS,ERRLCP,ERRACP,ERRTCP)
END IF !(CC > 0)
!-----------------------------------------------------------------------------------------------!
ERRLCS=LCMAXS
ERRACS=RACAES
ERRTCS=ETAS
ERRLCP=LCMAXP
ERRACP=RACAEP
ERRTCP=ETAP
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE CAVERRC
!-----------------------------------------------------------------------------------------------!
