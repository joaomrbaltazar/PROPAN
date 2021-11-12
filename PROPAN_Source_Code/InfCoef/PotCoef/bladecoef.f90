!-----------------------------------------------------------------------------------------------!
!    Blade Influence Coefficients                                                               !
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
SUBROUTINE BLADECOEF
!-----------------------------------------------------------------------------------------------!
!    Created by: J.A.C. Falcao de Campos, IST                                                   !
!    Modified  : 02122013, J. Baltazar, 2014 version 1.0                                        !
!    Modified  : 29052014, J. Baltazar, Revised                                                 !
!    Modified  : 27112014, J. Baltazar, 2014 version 3.3, Super-Cavitation Model                !
!    Modified  : 07012015, J. Baltazar, 2015 version 1.0, Wing Case                             !
!    Modified  : 11022016, J. Baltazar, Potential at field points included                      !
!    Modified  : 25102016, J. Baltazar, version 1.4                                             !
!    Modified  : 03072017, J. Baltazar, 2017 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,II,JJ,K,L,KB
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4),TE(4),XT(3),YT(3),ZT(3),TL,D1,D2,D,R0
DOUBLE PRECISION :: X0,Y0,Z0,UNX0,UNY0,UNZ0,A0
DOUBLE PRECISION :: A1X,A1Y,A1Z,A2X,A2Y,A2Z
DOUBLE PRECISION :: PHIS,PHID
!-----------------------------------------------------------------------------------------------!
!    Loop on Blade Panels                                                                       !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRP
   DO I=1,NCP
      K=NHPAN+(J-1)*NCP+I
!-----------------------------------------------------------------------------------------------!
!    Loop on the Number of Blades                                                               !
!-----------------------------------------------------------------------------------------------!
      DO KB=1,NB
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
!    Loop on Hub Control Points                                                                 !
!-----------------------------------------------------------------------------------------------!
         DO JJ=1,NHTP
            DO II=1,NHX
               L=(JJ-1)*NHX+II
!-----------------------------------------------------------------------------------------------!
!    Influence Coefficients                                                                     !
!-----------------------------------------------------------------------------------------------!
               R0=DSQRT((XH0(II,JJ)-X0)**2+(YH0(II,JJ)-Y0)**2+(ZH0(II,JJ)-Z0)**2)
               IF (IFARP == 1) THEN
!-----------------------------------------------------------------------------------------------!
                  IF (R0/D >= FARTOL2) CALL FARFIELD(X0,Y0,Z0,UNX0,UNY0,UNZ0,A0, &
                                                     XH0(II,JJ),YH0(II,JJ),ZH0(II,JJ),PHIS,PHID)
!-----------------------------------------------------------------------------------------------!
                  IF (R0/D >= FARTOL1.AND.R0/D < FARTOL2) THEN
!-----------------------------------------------------------------------------------------------!
                     IF (DABS(XX(3)-XX(4)) <= TOL.AND.J == NRP) THEN
                        XT(1:3)=XX(1:3)
                        YT(1:3)=YY(1:3)
                        ZT(1:3)=ZZ(1:3)
                        CALL GAUSSPAN3(4,XT,YT,ZT,XH0(II,JJ),YH0(II,JJ),ZH0(II,JJ),PHIS,PHID)
                     ELSEIF (ISTRIP == 2.AND.J == NRP) THEN
                        IF (I == 1.OR.I == NC+1) THEN
                           XT(1:3)=XX(1:3)
                           YT(1:3)=YY(1:3)
                           ZT(1:3)=ZZ(1:3)
                           CALL GAUSSPAN3(4,XT,YT,ZT,XH0(II,JJ),YH0(II,JJ),ZH0(II,JJ),PHIS,PHID)
                        ELSEIF (I == NC.OR.I == NCP) THEN
                           XT(1)=XX(1)
                           YT(1)=YY(1)
                           ZT(1)=ZZ(1)
                           XT(2)=XX(2)
                           YT(2)=YY(2)
                           ZT(2)=ZZ(2)
                           XT(3)=XX(4)
                           YT(3)=YY(4)
                           ZT(3)=ZZ(4)
                           CALL GAUSSPAN3(4,XT,YT,ZT,XH0(II,JJ),YH0(II,JJ),ZH0(II,JJ),PHIS,PHID)
                        ELSE
                           CALL GAUSSPAN(2,XX,YY,ZZ,XH0(II,JJ),YH0(II,JJ),ZH0(II,JJ),PHIS,PHID)
                        END IF
                     ELSE
                        CALL GAUSSPAN(2,XX,YY,ZZ,XH0(II,JJ),YH0(II,JJ),ZH0(II,JJ),PHIS,PHID)
                     END IF
!-----------------------------------------------------------------------------------------------!
                  END IF !(R0/D >= FARTOL1.AND.R0/D < FARTOL2)
!-----------------------------------------------------------------------------------------------!
               END IF !(IFARP == 1)
!-----------------------------------------------------------------------------------------------!
               IF (R0/D < FARTOL1.OR.IFARP == 0) THEN
!-----------------------------------------------------------------------------------------------!
                  IF (DABS(XX(3)-XX(4)) <= TOL.AND.J == NRP) THEN
                     XT(1:3)=XX(1:3)
                     YT(1:3)=YY(1:3)
                     ZT(1:3)=ZZ(1:3)
                     CALL POTPAN3(XT,YT,ZT,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                     XH0(II,JJ),YH0(II,JJ),ZH0(II,JJ),PHIS,PHID)
                  ELSEIF (ISTRIP == 2.AND.J == NRP) THEN
                     IF (I == 1.OR.I == NC+1) THEN
                        XT(1:3)=XX(1:3)
                        YT(1:3)=YY(1:3)
                        ZT(1:3)=ZZ(1:3)
                        CALL POTPAN3(XT,YT,ZT,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                     XH0(II,JJ),YH0(II,JJ),ZH0(II,JJ),PHIS,PHID)
                     ELSEIF (I == NC.OR.I == NCP) THEN
                        XT(1)=XX(1)
                        YT(1)=YY(1)
                        ZT(1)=ZZ(1)
                        XT(2)=XX(2)
                        YT(2)=YY(2)
                        ZT(2)=ZZ(2)
                        XT(3)=XX(4)
                        YT(3)=YY(4)
                        ZT(3)=ZZ(4)
                        CALL POTPAN3(XT,YT,ZT,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                     XH0(II,JJ),YH0(II,JJ),ZH0(II,JJ),PHIS,PHID)
                     ELSE
                        IF (INTE /= 1) CALL POTPANH(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                     XH0(II,JJ),YH0(II,JJ),ZH0(II,JJ),PHIS,PHID)
                        IF (INTE == 1) CALL MSTRPAN_ADP(TOLS,MMAX,XX,YY,ZZ, &
                                                     XH0(II,JJ),YH0(II,JJ),ZH0(II,JJ),PHIS,PHID)
                     END IF
                  ELSE
                     IF (INTE /= 1) CALL POTPANH(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                     XH0(II,JJ),YH0(II,JJ),ZH0(II,JJ),PHIS,PHID)
                     IF (INTE == 1) CALL MSTRPAN_ADP(TOLS,MMAX,XX,YY,ZZ, &
                                                     XH0(II,JJ),YH0(II,JJ),ZH0(II,JJ),PHIS,PHID)
                  END IF
!-----------------------------------------------------------------------------------------------!
               END IF !(R0/D < FARTOL1.OR.IFARP == 0)
!-----------------------------------------------------------------------------------------------!
!    Define Influence Coefficient                                                               !
!-----------------------------------------------------------------------------------------------!
               DIJ(L,K,KB)=PHID
               SIJ(L,K,KB)=PHIS
            END DO !II=1,NHX
         END DO !JJ=1,NHTP
!-----------------------------------------------------------------------------------------------!
!    Loop on Blade Control Points                                                               !
!-----------------------------------------------------------------------------------------------!
         DO JJ=1,NRP
            DO II=1,NCP
               L=NHPAN+(JJ-1)*NCP+II
!-----------------------------------------------------------------------------------------------!
!    Influence Coefficients                                                                     !
!-----------------------------------------------------------------------------------------------!
               R0=DSQRT((XP0(II,JJ)-X0)**2+(YP0(II,JJ)-Y0)**2+(ZP0(II,JJ)-Z0)**2)
               IF (IFARP == 1) THEN
!-----------------------------------------------------------------------------------------------!
                  IF (R0/D >= FARTOL2) CALL FARFIELD(X0,Y0,Z0,UNX0,UNY0,UNZ0,A0, &
                                                     XP0(II,JJ),YP0(II,JJ),ZP0(II,JJ),PHIS,PHID)
!-----------------------------------------------------------------------------------------------!
                  IF (R0/D >= FARTOL1.AND.R0/D < FARTOL2) THEN
!-----------------------------------------------------------------------------------------------!
                     IF (DABS(XX(3)-XX(4)) <= TOL.AND.J == NRP) THEN
                        XT(1:3)=XX(1:3)
                        YT(1:3)=YY(1:3)
                        ZT(1:3)=ZZ(1:3)
                        CALL GAUSSPAN3(4,XT,YT,ZT,XP0(II,JJ),YP0(II,JJ),ZP0(II,JJ),PHIS,PHID)
                     ELSEIF (ISTRIP == 2.AND.J == NRP) THEN
                        IF (I == 1.OR.I == NC+1) THEN
                           XT(1:3)=XX(1:3)
                           YT(1:3)=YY(1:3)
                           ZT(1:3)=ZZ(1:3)
                           CALL GAUSSPAN3(4,XT,YT,ZT,XP0(II,JJ),YP0(II,JJ),ZP0(II,JJ),PHIS,PHID)
                        ELSEIF (I == NC.OR.I == NCP) THEN
                           XT(1)=XX(1)
                           YT(1)=YY(1)
                           ZT(1)=ZZ(1)
                           XT(2)=XX(2)
                           YT(2)=YY(2)
                           ZT(2)=ZZ(2)
                           XT(3)=XX(4)
                           YT(3)=YY(4)
                           ZT(3)=ZZ(4)
                           CALL GAUSSPAN3(4,XT,YT,ZT,XP0(II,JJ),YP0(II,JJ),ZP0(II,JJ),PHIS,PHID)
                        ELSE
                           CALL GAUSSPAN(2,XX,YY,ZZ,XP0(II,JJ),YP0(II,JJ),ZP0(II,JJ),PHIS,PHID)
                        END IF
                     ELSE
                        CALL GAUSSPAN(2,XX,YY,ZZ,XP0(II,JJ),YP0(II,JJ),ZP0(II,JJ),PHIS,PHID)
                     END IF
!-----------------------------------------------------------------------------------------------!
                  END IF !(R0/D >= FARTOL1.AND.R0/D < FARTOL2)
!-----------------------------------------------------------------------------------------------!
               END IF !(IFARP == 1)
!-----------------------------------------------------------------------------------------------!
               IF (R0/D < FARTOL1.OR.IFARP == 0) THEN
!-----------------------------------------------------------------------------------------------!
                  IF (DABS(XX(3)-XX(4)) <= TOL.AND.J == NRP) THEN
                     XT(1:3)=XX(1:3)
                     YT(1:3)=YY(1:3)
                     ZT(1:3)=ZZ(1:3)
                     CALL POTPAN3(XT,YT,ZT,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                     XP0(II,JJ),YP0(II,JJ),ZP0(II,JJ),PHIS,PHID)
                  ELSEIF (ISTRIP == 2.AND.J == NRP) THEN
                     IF (I == 1.OR.I == NC+1) THEN
                        XT(1:3)=XX(1:3)
                        YT(1:3)=YY(1:3)
                        ZT(1:3)=ZZ(1:3)
                        CALL POTPAN3(XT,YT,ZT,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                     XP0(II,JJ),YP0(II,JJ),ZP0(II,JJ),PHIS,PHID)
                     ELSEIF (I == NC.OR.I == NCP) THEN
                        XT(1)=XX(1)
                        YT(1)=YY(1)
                        ZT(1)=ZZ(1)
                        XT(2)=XX(2)
                        YT(2)=YY(2)
                        ZT(2)=ZZ(2)
                        XT(3)=XX(4)
                        YT(3)=YY(4)
                        ZT(3)=ZZ(4)
                        CALL POTPAN3(XT,YT,ZT,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                     XP0(II,JJ),YP0(II,JJ),ZP0(II,JJ),PHIS,PHID)
                     ELSE
                        IF (INTE /= 1) CALL POTPANH(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                     XP0(II,JJ),YP0(II,JJ),ZP0(II,JJ),PHIS,PHID)
                        IF (INTE == 1) CALL MSTRPAN_ADP(TOLS,MMAX,XX,YY,ZZ, &
                                                     XP0(II,JJ),YP0(II,JJ),ZP0(II,JJ),PHIS,PHID)
                     END IF
                  ELSE
                     IF (INTE /= 1) CALL POTPANH(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                     XP0(II,JJ),YP0(II,JJ),ZP0(II,JJ),PHIS,PHID)
                     IF (INTE == 1) CALL MSTRPAN_ADP(TOLS,MMAX,XX,YY,ZZ, &
                                                     XP0(II,JJ),YP0(II,JJ),ZP0(II,JJ),PHIS,PHID)
                  END IF
!-----------------------------------------------------------------------------------------------!
               END IF !(R0/D < FARTOL1.OR.IFARP == 0)
!-----------------------------------------------------------------------------------------------!
!    Define Influence Coefficient                                                               !
!-----------------------------------------------------------------------------------------------!
               DIJ(L,K,KB)=PHID
               SIJ(L,K,KB)=PHIS
            END DO !II=1,NCP
         END DO !JJ=1,NRP
!-----------------------------------------------------------------------------------------------!
!    Loop on Nozzle Control Points                                                              !
!-----------------------------------------------------------------------------------------------!
         DO JJ=1,NNTP
            DO II=1,NNXT1
               L=NHPAN+NPPAN+(JJ-1)*NNXT1+II
!-----------------------------------------------------------------------------------------------!
!    Influence Coefficients                                                                     !
!-----------------------------------------------------------------------------------------------!
               R0=DSQRT((XN0(II,JJ)-X0)**2+(YN0(II,JJ)-Y0)**2+(ZN0(II,JJ)-Z0)**2)
               IF (IFARP == 1) THEN
!-----------------------------------------------------------------------------------------------!
                  IF (R0/D >= FARTOL2) CALL FARFIELD(X0,Y0,Z0,UNX0,UNY0,UNZ0,A0, &
                                                     XN0(II,JJ),YN0(II,JJ),ZN0(II,JJ),PHIS,PHID)
!-----------------------------------------------------------------------------------------------!
                  IF (R0/D >= FARTOL1.AND.R0/D < FARTOL2) THEN
!-----------------------------------------------------------------------------------------------!
                     IF (DABS(XX(3)-XX(4)) <= TOL.AND.J == NRP) THEN
                        XT(1:3)=XX(1:3)
                        YT(1:3)=YY(1:3)
                        ZT(1:3)=ZZ(1:3)
                        CALL GAUSSPAN3(4,XT,YT,ZT,XN0(II,JJ),YN0(II,JJ),ZN0(II,JJ),PHIS,PHID)
                     ELSEIF (ISTRIP == 2.AND.J == NRP) THEN
                        IF (I == 1.OR.I == NC+1) THEN
                           XT(1:3)=XX(1:3)
                           YT(1:3)=YY(1:3)
                           ZT(1:3)=ZZ(1:3)
                           CALL GAUSSPAN3(4,XT,YT,ZT,XN0(II,JJ),YN0(II,JJ),ZN0(II,JJ),PHIS,PHID)
                        ELSEIF (I == NC.OR.I == NCP) THEN
                           XT(1)=XX(1)
                           YT(1)=YY(1)
                           ZT(1)=ZZ(1)
                           XT(2)=XX(2)
                           YT(2)=YY(2)
                           ZT(2)=ZZ(2)
                           XT(3)=XX(4)
                           YT(3)=YY(4)
                           ZT(3)=ZZ(4)
                           CALL GAUSSPAN3(4,XT,YT,ZT,XN0(II,JJ),YN0(II,JJ),ZN0(II,JJ),PHIS,PHID)
                        ELSE
                           CALL GAUSSPAN(2,XX,YY,ZZ,XN0(II,JJ),YN0(II,JJ),ZN0(II,JJ),PHIS,PHID)
                        END IF
                     ELSE
                        CALL GAUSSPAN(2,XX,YY,ZZ,XN0(II,JJ),YN0(II,JJ),ZN0(II,JJ),PHIS,PHID)
                     END IF
!-----------------------------------------------------------------------------------------------!
                  END IF !(R0/D >= FARTOL1.AND.R0/D < FARTOL2)
!-----------------------------------------------------------------------------------------------!
               END IF !(IFARP == 1)
!-----------------------------------------------------------------------------------------------!
               IF (R0/D < FARTOL1.OR.IFARP == 0) THEN
!-----------------------------------------------------------------------------------------------!
                  IF (DABS(XX(3)-XX(4)) <= TOL.AND.J == NRP) THEN
                     XT(1:3)=XX(1:3)
                     YT(1:3)=YY(1:3)
                     ZT(1:3)=ZZ(1:3)
                     CALL POTPAN3(XT,YT,ZT,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                     XN0(II,JJ),YN0(II,JJ),ZN0(II,JJ),PHIS,PHID)
                  ELSEIF (ISTRIP == 2.AND.J == NRP) THEN
                     IF (I == 1.OR.I == NC+1) THEN
                        XT(1:3)=XX(1:3)
                        YT(1:3)=YY(1:3)
                        ZT(1:3)=ZZ(1:3)
                        CALL POTPAN3(XT,YT,ZT,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                     XN0(II,JJ),YN0(II,JJ),ZN0(II,JJ),PHIS,PHID)
                     ELSEIF (I == NC.OR.I == NCP) THEN
                        XT(1)=XX(1)
                        YT(1)=YY(1)
                        ZT(1)=ZZ(1)
                        XT(2)=XX(2)
                        YT(2)=YY(2)
                        ZT(2)=ZZ(2)
                        XT(3)=XX(4)
                        YT(3)=YY(4)
                        ZT(3)=ZZ(4)
                        CALL POTPAN3(XT,YT,ZT,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                     XN0(II,JJ),YN0(II,JJ),ZN0(II,JJ),PHIS,PHID)
                     ELSE
                        IF (INTE /= 1) CALL POTPANH(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                     XN0(II,JJ),YN0(II,JJ),ZN0(II,JJ),PHIS,PHID)
                        IF (INTE == 1) CALL MSTRPAN_ADP(TOLS,MMAX,XX,YY,ZZ, &
                                                     XN0(II,JJ),YN0(II,JJ),ZN0(II,JJ),PHIS,PHID)
                     END IF
                  ELSE
                     IF (INTE /= 1) CALL POTPANH(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                     XN0(II,JJ),YN0(II,JJ),ZN0(II,JJ),PHIS,PHID)
                     IF (INTE == 1) CALL MSTRPAN_ADP(TOLS,MMAX,XX,YY,ZZ, &
                                                     XN0(II,JJ),YN0(II,JJ),ZN0(II,JJ),PHIS,PHID)
                  END IF
!-----------------------------------------------------------------------------------------------!
               END IF !(R0/D < FARTOL1.OR.IFARP == 0)
!-----------------------------------------------------------------------------------------------!
!    Define Influence Coefficient                                                               !
!-----------------------------------------------------------------------------------------------!
               DIJ(L,K,KB)=PHID
               SIJ(L,K,KB)=PHIS
            END DO !II=1,NNXT1
         END DO !JJ=1,NNTP
!-----------------------------------------------------------------------------------------------!
!    Loop on Blade Wake Control Points                                                          !
!-----------------------------------------------------------------------------------------------!
         DO JJ=1,NRW
            DO II=1,IABS(NCPW)
               L=NHPAN+NPPAN+NNPAN+(JJ-1)*IABS(NCPW)+II
!-----------------------------------------------------------------------------------------------!
!    Influence Coefficients                                                                     !
!-----------------------------------------------------------------------------------------------!
               R0=DSQRT((XPW0(II,JJ)-X0)**2+(YPW0(II,JJ)-Y0)**2+(ZPW0(II,JJ)-Z0)**2)
               IF (IFARP == 1) THEN
!-----------------------------------------------------------------------------------------------!
                  IF (R0/D >= FARTOL2) CALL FARFIELD(X0,Y0,Z0,UNX0,UNY0,UNZ0,A0, &
                                                  XPW0(II,JJ),YPW0(II,JJ),ZPW0(II,JJ),PHIS,PHID)
!-----------------------------------------------------------------------------------------------!
                  IF (R0/D >= FARTOL1.AND.R0/D < FARTOL2) THEN
!-----------------------------------------------------------------------------------------------!
                     IF (DABS(XX(3)-XX(4)) <= TOL.AND.J == NRP) THEN
                        XT(1:3)=XX(1:3)
                        YT(1:3)=YY(1:3)
                        ZT(1:3)=ZZ(1:3)
                        CALL GAUSSPAN3(4,XT,YT,ZT,XPW0(II,JJ),YPW0(II,JJ),ZPW0(II,JJ),PHIS,PHID)
                     ELSEIF (ISTRIP == 2.AND.J == NRP) THEN
                        IF (I == 1.OR.I == NC+1) THEN
                           XT(1:3)=XX(1:3)
                           YT(1:3)=YY(1:3)
                           ZT(1:3)=ZZ(1:3)
                           CALL GAUSSPAN3(4,XT,YT,ZT, &
                                                  XPW0(II,JJ),YPW0(II,JJ),ZPW0(II,JJ),PHIS,PHID)
                        ELSEIF (I == NC.OR.I == NCP) THEN
                           XT(1)=XX(1)
                           YT(1)=YY(1)
                           ZT(1)=ZZ(1)
                           XT(2)=XX(2)
                           YT(2)=YY(2)
                           ZT(2)=ZZ(2)
                           XT(3)=XX(4)
                           YT(3)=YY(4)
                           ZT(3)=ZZ(4)
                           CALL GAUSSPAN3(4,XT,YT,ZT, &
                                                   XPW0(II,JJ),YPW0(II,JJ),ZPW0(II,JJ),PHIS,PHID)
                        ELSE
                           CALL GAUSSPAN(2,XX,YY,ZZ, &
                                                   XPW0(II,JJ),YPW0(II,JJ),ZPW0(II,JJ),PHIS,PHID)
                        END IF
                     ELSE
                        CALL GAUSSPAN(2,XX,YY,ZZ,XPW0(II,JJ),YPW0(II,JJ),ZPW0(II,JJ),PHIS,PHID)
                     END IF
!-----------------------------------------------------------------------------------------------!
                  END IF !(R0/D >= FARTOL1.AND.R0/D < FARTOL2)
!-----------------------------------------------------------------------------------------------!
               END IF !(IFARP == 1)
!-----------------------------------------------------------------------------------------------!
               IF (R0/D < FARTOL1.OR.IFARP == 0) THEN
!-----------------------------------------------------------------------------------------------!
                  IF (DABS(XX(3)-XX(4)) <= TOL.AND.J == NRP) THEN
                     XT(1:3)=XX(1:3)
                     YT(1:3)=YY(1:3)
                     ZT(1:3)=ZZ(1:3)
                     CALL POTPAN3(XT,YT,ZT,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                  XPW0(II,JJ),YPW0(II,JJ),ZPW0(II,JJ),PHIS,PHID)
                  ELSEIF (ISTRIP == 2.AND.J == NRP) THEN
                     IF (I == 1.OR.I == NC+1) THEN
                        XT(1:3)=XX(1:3)
                        YT(1:3)=YY(1:3)
                        ZT(1:3)=ZZ(1:3)
                        CALL POTPAN3(XT,YT,ZT,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                  XPW0(II,JJ),YPW0(II,JJ),ZPW0(II,JJ),PHIS,PHID)
                     ELSEIF (I == NC.OR.I == NCP) THEN
                        XT(1)=XX(1)
                        YT(1)=YY(1)
                        ZT(1)=ZZ(1)
                        XT(2)=XX(2)
                        YT(2)=YY(2)
                        ZT(2)=ZZ(2)
                        XT(3)=XX(4)
                        YT(3)=YY(4)
                        ZT(3)=ZZ(4)
                        CALL POTPAN3(XT,YT,ZT,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                  XPW0(II,JJ),YPW0(II,JJ),ZPW0(II,JJ),PHIS,PHID)
                     ELSE
                        IF (INTE /= 1) CALL POTPANH(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                  XPW0(II,JJ),YPW0(II,JJ),ZPW0(II,JJ),PHIS,PHID)
                        IF (INTE == 1) CALL MSTRPAN_ADP(TOLS,MMAX,XX,YY,ZZ, &
                                                  XPW0(II,JJ),YPW0(II,JJ),ZPW0(II,JJ),PHIS,PHID)
                     END IF
                  ELSE
                     IF (INTE /= 1) CALL POTPANH(XX,YY,ZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0, &
                                                  XPW0(II,JJ),YPW0(II,JJ),ZPW0(II,JJ),PHIS,PHID)
                     IF (INTE == 1) CALL MSTRPAN_ADP(TOLS,MMAX,XX,YY,ZZ, &
                                                  XPW0(II,JJ),YPW0(II,JJ),ZPW0(II,JJ),PHIS,PHID)
                  END IF
!-----------------------------------------------------------------------------------------------!
               END IF !(R0/D < FARTOL1.OR.IFARP == 0)
!-----------------------------------------------------------------------------------------------!
!    Define Influence Coefficient                                                               !
!-----------------------------------------------------------------------------------------------!
               DIJ(L,K,KB)=PHID
               SIJ(L,K,KB)=PHIS
            END DO !II=1,IABS(NCPW)
         END DO !JJ=1,NRW
!-----------------------------------------------------------------------------------------------!
      END DO !KB=1,NB
!-----------------------------------------------------------------------------------------------!
   END DO !I=1,NCP
END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE BLADECOEF
!-----------------------------------------------------------------------------------------------!
