!-----------------------------------------------------------------------------------------------!
!    Blade Wake Potential Calculation                                                           !
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
SUBROUTINE BLADEWAKECOEFF(TT,N1,N2,XTMP,YTMP,ZTMP,PTMP)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 25102016, J. Baltazar, 2016 version 1.4                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,II,JJ,K,KK,L,KB,TT,N1,N2
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4),TE(4),TL,D1,D2,D,R0
DOUBLE PRECISION :: X0,Y0,Z0,UNX0,UNY0,UNZ0,A0
DOUBLE PRECISION :: A1X,A1Y,A1Z,A2X,A2Y,A2Z
DOUBLE PRECISION :: PHIS,PHID,PDL,PDR,POT(NPW,NRW)
DOUBLE PRECISION :: XTMP(N1,N2),YTMP(N1,N2),ZTMP(N1,N2),PTMP(N1,N2)
!-----------------------------------------------------------------------------------------------!
PTMP=0.D0
!-----------------------------------------------------------------------------------------------!
!    Loop on Blade Wake Panels                                                                  !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRW
   DO I=1,NPW
      K=(J-1)*NPW+I
!-----------------------------------------------------------------------------------------------!
!    Loop on the Number of Blades                                                               !
!-----------------------------------------------------------------------------------------------!
      DO KB=1,NB
         IF (KB > 1) CALL PERIODICFLOW(TT,KB,NPW,NRW,NT,POTPW,POT)
!-----------------------------------------------------------------------------------------------!
!    Define Panel                                                                               !
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
         YY(1)=RPW(I  ,J  )*DCOS(TE(1))
         YY(2)=RPW(I+1,J  )*DCOS(TE(2))
         YY(3)=RPW(I+1,J+1)*DCOS(TE(3))
         YY(4)=RPW(I  ,J+1)*DCOS(TE(4))
         ZZ(1)=RPW(I  ,J  )*DSIN(TE(1))
         ZZ(2)=RPW(I+1,J  )*DSIN(TE(2))
         ZZ(3)=RPW(I+1,J+1)*DSIN(TE(3))
         ZZ(4)=RPW(I  ,J+1)*DSIN(TE(4))
	 IF ((IROTOR == -1).AND.(KB == 2)) THEN
            XX(1)= XPW(I+1,J  )
            XX(2)= XPW(I  ,J  )
            XX(3)= XPW(I  ,J+1)
            XX(4)= XPW(I+1,J+1)
            YY(1)=-YPW(I+1,J  )
            YY(2)=-YPW(I  ,J  )
            YY(3)=-YPW(I  ,J+1)
            YY(4)=-YPW(I+1,J+1)
            ZZ(1)= ZPW(I+1,J  )
            ZZ(2)= ZPW(I  ,J  )
            ZZ(3)= ZPW(I  ,J+1)
            ZZ(4)= ZPW(I+1,J+1)
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
!    Loop on Field Points                                                                       !
!-----------------------------------------------------------------------------------------------!
         DO II=1,N1
            DO JJ=1,N2
               R0=DSQRT((XTMP(II,JJ)-X0)**2+(YTMP(II,JJ)-Y0)**2+(ZTMP(II,JJ)-Z0)**2)
!-----------------------------------------------------------------------------------------------!
!    Potential Influence Coefficients                                                           !
!-----------------------------------------------------------------------------------------------!
               IF (IFARP == 1) THEN
                  IF (R0/D >= FARTOL2) CALL FARFIELD(X0,Y0,Z0,UNX0,UNY0,UNZ0,A0, &
                                                  XTMP(II,JJ),YTMP(II,JJ),ZTMP(II,JJ),PHIS,PHID)
                  IF (R0/D >= FARTOL1.AND.R0/D < FARTOL2) CALL GAUSSPAN(2,XX,YY,ZZ, &
                                                  XTMP(II,JJ),YTMP(II,JJ),ZTMP(II,JJ),PHIS,PHID)
               END IF !(IFARP == 1)
!-----------------------------------------------------------------------------------------------!
               IF (R0/D < FARTOL1.OR.IFARP == 0) THEN
                  IF (INTE == 0) CALL POTPANH(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                  XTMP(II,JJ),YTMP(II,JJ),ZTMP(II,JJ),PHIS,PHID)
                  IF (INTE /= 0) CALL MSTRPAN_ADP(TOLS,MMAX,XX,YY,ZZ, &
                                                  XTMP(II,JJ),YTMP(II,JJ),ZTMP(II,JJ),PHIS,PHID)
               END IF !(R0/D < FARTOL1.OR.IFARP == 0)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake Potential Calculation                                                           !
!-----------------------------------------------------------------------------------------------!
               IF (KB == 1) THEN
                  PTMP(II,JJ)=PTMP(II,JJ)+PHID*POTPW(I,J,TT)/(4.D0*PI)
               ELSE !(KB == 1)
                  PTMP(II,JJ)=PTMP(II,JJ)+PHID*POT(I,J)/(4.D0*PI)
               END IF !(KB ==1 )
!-----------------------------------------------------------------------------------------------!
            END DO !JJ=1,N2
         END DO !II=1,N1
!-----------------------------------------------------------------------------------------------!
      END DO !KB=1,NB
!-----------------------------------------------------------------------------------------------!
   END DO !I=1,NPW
END DO !J=1,NRW
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE BLADEWAKECOEFF
!-----------------------------------------------------------------------------------------------!
