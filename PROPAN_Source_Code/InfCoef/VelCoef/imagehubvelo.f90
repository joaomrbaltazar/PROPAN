!-----------------------------------------------------------------------------------------------!
!    Image hub influence coefficients                                                           !
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
SUBROUTINE IMAGEHUBVELO(TT,N1,N2,XTMP,YTMP,ZTMP,UTMP,VTMP,WTMP)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 03122013, J. Baltazar, version 1.0                                             !
!    Modified  : 01122014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,II,JJ,KB,N1,N2,TT
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4),TE(4),TL,RR
DOUBLE PRECISION :: X0,Y0,Z0,UNX0,UNY0,UNZ0,A0
DOUBLE PRECISION :: A1X,A1Y,A1Z,A2X,A2Y,A2Z
DOUBLE PRECISION :: VS(3),VD(3)
DOUBLE PRECISION :: POT(NPW,NRW)
DOUBLE PRECISION :: XTMP(N1,N2),YTMP(N1,N2),ZTMP(N1,N2),UTMP(N1,N2),VTMP(N1,N2),WTMP(N1,N2)
!-----------------------------------------------------------------------------------------------!
UTMP=0.D0
VTMP=0.D0
WTMP=0.D0
!-----------------------------------------------------------------------------------------------!
!    Loop on Blade Wake Panels                                                                  !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRW
   DO I=1,NPW
!-----------------------------------------------------------------------------------------------!
!    Loop on the Number of Blades                                                               !
!-----------------------------------------------------------------------------------------------!
      DO KB=1,NB
         IF (KB > 1) CALL PERIODICFLOW(TT,KB,NPW,NRW,NT,POTPW,POT)
!-----------------------------------------------------------------------------------------------!
!    Interpolate the Radius                                                                     !
!-----------------------------------------------------------------------------------------------!
         IF (INTERH == 0) THEN
            CALL LININT(NHI,XHI,RHI,1,XPW(I,J),RR)
         ELSEIF (INTERH == 1) THEN
            CALL INTK1 (NHI,XHI,RHI,1,XPW(I,J),RR)
         ELSEIF (INTERH == 2) THEN
            CALL SPLINT(NHI,XHI,RHI,1,XPW(I,J),RR)
         END IF !(INTERH)
!-----------------------------------------------------------------------------------------------!
!    Define Segment                                                                             !
!-----------------------------------------------------------------------------------------------!
         TL=DFLOAT(KB-1)*2.D0*PI/DFLOAT(NB)
         TE(1)=TPW(I  ,J  )+TL
         TE(2)=TPW(I+1,J  )+TL
         TE(3)=TPW(I+1,J+1)+TL
         TE(4)=TPW(I  ,J+1)+TL
         XX(1)=XPW(I  ,J  )
         XX(2)=XPW(I+1,J  )
         XX(3)=XPW(I+1,J+1)
         XX(4)=XPW(I  ,J+1)
         YY(1)=RR*RR/RPW(I  ,J  )*DCOS(TE(1))
         YY(2)=RR*RR/RPW(I+1,J  )*DCOS(TE(2))
         YY(3)=RR*RR/RPW(I+1,J+1)*DCOS(TE(3))
         YY(4)=RR*RR/RPW(I  ,J+1)*DCOS(TE(4))
         ZZ(1)=RR*RR/RPW(I  ,J  )*DSIN(TE(1))
         ZZ(2)=RR*RR/RPW(I+1,J  )*DSIN(TE(2))
         ZZ(3)=RR*RR/RPW(I+1,J+1)*DSIN(TE(3))
         ZZ(4)=RR*RR/RPW(I  ,J+1)*DSIN(TE(4))
!-----------------------------------------------------------------------------------------------!
!    Compute Panel Centroid Data                                                                !
!-----------------------------------------------------------------------------------------------!
         CALL PANEL(XX,YY,ZZ,X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0)
!-----------------------------------------------------------------------------------------------!
!    For Flat Panel Redefine Corner Points                                                      !
!-----------------------------------------------------------------------------------------------!
         IF (IPAN == 0) CALL PANELFLAT(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0)
!-----------------------------------------------------------------------------------------------!
!    Loop on Wake Control Points                                                                !
!-----------------------------------------------------------------------------------------------!
         DO II=1,N1
            DO JJ=1,N2
!-----------------------------------------------------------------------------------------------!
!    Velocity Influence Coefficient                                                             !
!-----------------------------------------------------------------------------------------------!
               CALL IMAGEVELO(2,XX,YY,ZZ,XTMP(II,JJ),YTMP(II,JJ),ZTMP(II,JJ),EPS,VD)
!-----------------------------------------------------------------------------------------------!
!    Wake Velocity Calculation                                                                  !
!-----------------------------------------------------------------------------------------------!
               IF (KB == 1) THEN
                  UTMP(II,JJ)=UTMP(II,JJ)+VD(1)*POTPW(I,J,TT)/(4.D0*PI)
                  VTMP(II,JJ)=VTMP(II,JJ)+VD(2)*POTPW(I,J,TT)/(4.D0*PI)
                  WTMP(II,JJ)=WTMP(II,JJ)+VD(3)*POTPW(I,J,TT)/(4.D0*PI)
               ELSE !(KB == 1)
                  UTMP(II,JJ)=UTMP(II,JJ)+VD(1)*POT(I,J)/(4.D0*PI)
                  VTMP(II,JJ)=VTMP(II,JJ)+VD(2)*POT(I,J)/(4.D0*PI)
                  WTMP(II,JJ)=WTMP(II,JJ)+VD(3)*POT(I,J)/(4.D0*PI)
               END IF !(KB == 1)
!-----------------------------------------------------------------------------------------------!
            END DO !JJ=1,N2
         END DO !II=1,N1
!-----------------------------------------------------------------------------------------------!
      END DO !KB=1,NB
!-----------------------------------------------------------------------------------------------!
   END DO !I=1,NPW
END DO !J=1,NRW
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE IMAGEHUBVELO
!-----------------------------------------------------------------------------------------------!
