!-----------------------------------------------------------------------------------------------!
!    Relation between rotating and inertial reference frames                                    !
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
SUBROUTINE FRAME(CREF,TT,KB,X,Y,Z,XNEW,YNEW,ZNEW,RNEW,TNEW)
!-----------------------------------------------------------------------------------------------!
!    Created by: 05022016, J. Baltazar, version 1.0                                             !
!    Modified  : 27042016, J. Baltazar, version 1.1                                             !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
CHARACTER*8 CREF
INTEGER :: TT,T1,KB
DOUBLE PRECISION :: X,Y,Z,R,T,XNEW,YNEW,ZNEW,RNEW,TNEW
!-----------------------------------------------------------------------------------------------!
!    Rotating Reference Frame to Inertial Reference Frame                                       !
!-----------------------------------------------------------------------------------------------!
IF (CREF == 'INERTIAL') THEN
   T1=TT
   IF (NTETA /= 0) T1=TT-(TT-1)/NTETA*DFLOAT(NTETA)
   R=DSQRT(Y*Y+Z*Z)
   T=DATAN2(Z,Y)-DTETA*(DFLOAT(T1)-1.D0-P*(DFLOAT(KB)-1.D0)) !2.D0*PI+
   do while (t < 0.d0)
      t=t+2.d0*pi
   end do
!-----------------------------------------------------------------------------------------------!
   XNEW=X
   RNEW=R
   TNEW=T
   YNEW=RNEW*DCOS(TNEW)
   ZNEW=RNEW*DSIN(TNEW)
END IF !(CREF == 'INERTIAL')
!-----------------------------------------------------------------------------------------------!
!    Inertial Reference Frame to Rotating Reference Frame                                       !
!-----------------------------------------------------------------------------------------------!
IF (CREF == 'ROTATING') THEN
   T1=TT
   IF (NTETA /= 0) T1=TT-(TT-1)/NTETA*DFLOAT(NTETA)
   R=DSQRT(Y*Y+Z*Z)
   T=DATAN2(Z,Y)+DTETA*(DFLOAT(T1)-1.D0-P*(DFLOAT(KB)-1.D0)) !-2.D0*PI
   do while (t < 0.d0)
      t=t+2.d0*pi
   end do
!-----------------------------------------------------------------------------------------------!
   XNEW=X
   RNEW=R
   TNEW=T
   YNEW=RNEW*DCOS(TNEW)
   ZNEW=RNEW*DSIN(TNEW)
END IF !(CREF == 'ROTATING')
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE FRAME
!-----------------------------------------------------------------------------------------------!
