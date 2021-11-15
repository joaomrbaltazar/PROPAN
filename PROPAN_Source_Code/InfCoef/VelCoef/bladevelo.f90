!-----------------------------------------------------------------------------------------------!
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
SUBROUTINE BLADEVELO(TT,N1,N2,XTMP,YTMP,ZTMP,UTMP,VTMP,WTMP)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 02122013, J. Baltazar, version 1.0                                             !
!    Modified  : 30052014, J. Baltazar, revised                                                 !
!    Modified  : 30112014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!    Modified  : 07012015, J. Baltazar, 2015 version 1.0, Wing Case                             !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,II,JJ,KB,N1,N2,TT
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4),TE(4),TL,D1,D2,D,R0
DOUBLE PRECISION :: X0,Y0,Z0,UNX0,UNY0,UNZ0,A0
DOUBLE PRECISION :: A1X,A1Y,A1Z,A2X,A2Y,A2Z
DOUBLE PRECISION :: VS(3),VD(3)
DOUBLE PRECISION :: POT(NCP,NRP),SOURCE(NCP,NRP)
DOUBLE PRECISION :: XTMP(N1,N2),YTMP(N1,N2),ZTMP(N1,N2),UTMP(N1,N2),VTMP(N1,N2),WTMP(N1,N2)
!-----------------------------------------------------------------------------------------------!
UTMP=0.D0
VTMP=0.D0
WTMP=0.D0
!-----------------------------------------------------------------------------------------------!
!    Loop on Blade Panels                                                                       !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRP
   DO I=1,NCP
!-----------------------------------------------------------------------------------------------!
!    Loop on the Number of Blades                                                               !
!-----------------------------------------------------------------------------------------------!
      DO KB=1,NB
         IF (KB > 1) CALL PERIODICFLOW(TT,KB,NCP,NRP,NT,POTP,POT)
         IF ((KB > 1).AND.(NCAV > 0)) CALL PERIODICFLOW(TT,KB,NCP,NRP,NT,SOURCEPCAV,SOURCE)
!-----------------------------------------------------------------------------------------------!
!    Define Panel                                                                               !
!-----------------------------------------------------------------------------------------------!
         TL=DFLOAT(KB-1)*2.D0*PI/DFLOAT(NB)
         TE(1)=TP(I  ,J  )+TL
         TE(2)=TP(I+1,J  )+TL
         TE(3)=TP(I+1,J+1)+TL
         TE(4)=TP(I  ,J+1)+TL
         XX(1)=XP(I  ,J  )
         XX(2)=XP(I+1,J  )
         XX(3)=XP(I+1,J+1)
         XX(4)=XP(I  ,J+1)
         YY(1)=RP(I  ,J  )*DCOS(TE(1))
         YY(2)=RP(I+1,J  )*DCOS(TE(2))
         YY(3)=RP(I+1,J+1)*DCOS(TE(3))
         YY(4)=RP(I  ,J+1)*DCOS(TE(4))
         ZZ(1)=RP(I  ,J  )*DSIN(TE(1))
         ZZ(2)=RP(I+1,J  )*DSIN(TE(2))
         ZZ(3)=RP(I+1,J+1)*DSIN(TE(3))
         ZZ(4)=RP(I  ,J+1)*DSIN(TE(4))
	       IF ((IROTOR == -1).AND.(KB == 2)) THEN
            XX(1)= XP(I+1,J  )
            XX(2)= XP(I  ,J  )
            XX(3)= XP(I  ,J+1)
            XX(4)= XP(I+1,J+1)
            YY(1)=-YP(I+1,J  )
            YY(2)=-YP(I  ,J  )
            YY(3)=-YP(I  ,J+1)
            YY(4)=-YP(I+1,J+1)
            ZZ(1)= ZP(I+1,J  )
            ZZ(2)= ZP(I  ,J  )
            ZZ(3)= ZP(I  ,J+1)
            ZZ(4)= ZP(I+1,J+1)
         END IF !((IROTOR == -1).AND.(KB == 2))
!-----------------------------------------------------------------------------------------------!
!    Compute Panel Centroid Data                                                                !
!-----------------------------------------------------------------------------------------------!
         CALL PANEL(XX,YY,ZZ,X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0)
!-----------------------------------------------------------------------------------------------!
!    For Flat Panel Redefine Corner Points                                                      !
!-----------------------------------------------------------------------------------------------!
         IF (IPAN == 0) CALL PANELFLAT(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0)
!-----------------------------------------------------------------------------------------------!
!    Compute Panel Diagonals                                                                    !
!-----------------------------------------------------------------------------------------------!
         D1=DSQRT((XX(3)-XX(1))**2+(YY(3)-YY(1))**2+(ZZ(3)-ZZ(1))**2)
         D=D1
         D2=DSQRT((XX(4)-XX(2))**2+(YY(4)-YY(2))**2+(ZZ(4)-ZZ(2))**2)
         IF (D2 > D1) D=D2
!-----------------------------------------------------------------------------------------------!
!    Loop on Wake Control Points                                                                !
!-----------------------------------------------------------------------------------------------!
         DO II=1,N1
            DO JJ=1,N2
               R0=DSQRT((XTMP(II,JJ)-X0)**2+(YTMP(II,JJ)-Y0)**2+(ZTMP(II,JJ)-Z0)**2)
!-----------------------------------------------------------------------------------------------!
!    Velocity Influence Coefficient                                                             !
!-----------------------------------------------------------------------------------------------!
               IF (IFARV == 1) THEN
                  IF (R0/D >= FARTOL2) CALL GAUSSVEL(4,XX,YY,ZZ, &
                                                     XTMP(II,JJ),YTMP(II,JJ),ZTMP(II,JJ),VS,VD)
               END IF !(IFARV == 1)
!-----------------------------------------------------------------------------------------------!
               IF (R0/D < FARTOL2.OR.IFARV == 0) THEN
                  CALL VELPAN(XX,YY,ZZ,X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z, &
                              XTMP(II,JJ),YTMP(II,JJ),ZTMP(II,JJ),EPS,VS,VD)
               END IF !(R0/D < FARTOL2.OR.IFARV == 0)
!-----------------------------------------------------------------------------------------------!
!    Blade Velocity Calculation                                                                 !
!-----------------------------------------------------------------------------------------------!
               IF (KB == 1) THEN
                  UTMP(II,JJ)=UTMP(II,JJ)+(VD(1)*POTP(I,J,TT)+VS(1)*SOURCEP(I,J,1))/(4.D0*PI)
                  VTMP(II,JJ)=VTMP(II,JJ)+(VD(2)*POTP(I,J,TT)+VS(2)*SOURCEP(I,J,1))/(4.D0*PI)
                  WTMP(II,JJ)=WTMP(II,JJ)+(VD(3)*POTP(I,J,TT)+VS(3)*SOURCEP(I,J,1))/(4.D0*PI)
               ELSEIF (NCAV > 0) THEN
                  UTMP(II,JJ)=UTMP(II,JJ)+(VD(1)*POT(I,J)+VS(1)*SOURCE(I,J))/(4.D0*PI)
                  VTMP(II,JJ)=VTMP(II,JJ)+(VD(2)*POT(I,J)+VS(2)*SOURCE(I,J))/(4.D0*PI)
                  WTMP(II,JJ)=WTMP(II,JJ)+(VD(3)*POT(I,J)+VS(3)*SOURCE(I,J))/(4.D0*PI)
               ELSE
                  UTMP(II,JJ)=UTMP(II,JJ)+(VD(1)*POT(I,J)+VS(1)*SOURCEP(I,J,KB))/(4.D0*PI)
                  VTMP(II,JJ)=VTMP(II,JJ)+(VD(2)*POT(I,J)+VS(2)*SOURCEP(I,J,KB))/(4.D0*PI)
                  WTMP(II,JJ)=WTMP(II,JJ)+(VD(3)*POT(I,J)+VS(3)*SOURCEP(I,J,KB))/(4.D0*PI)
               END IF
!-----------------------------------------------------------------------------------------------!
            END DO !JJ=1,N2
         END DO !II=1,N1
!-----------------------------------------------------------------------------------------------!
      END DO !KB=1,NB
!-----------------------------------------------------------------------------------------------!
   END DO !I=1,NCP
END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE BLADEVELO
!-----------------------------------------------------------------------------------------------!
