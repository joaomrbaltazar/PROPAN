!-----------------------------------------------------------------------------------------------!
!    Induced potential by numerical integration                                                 !
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
SUBROUTINE MSTRPAN_ADP(TOL,MMAX,XX,YY,ZZ,X,Y,Z,PS,PD)
!-----------------------------------------------------------------------------------------------!
!    Original : MSTRPAN                                                                         !
!    Modified by : Joao Baltazar, IST, April 2001                                               !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
INTEGER :: I,NP,N,INT,MMAX
DOUBLE PRECISION :: TOL,XX(4),YY(4),ZZ(4),X,Y,Z,PS,PD
DOUBLE PRECISION :: RL1,RL2,ASP,STEP,THETA,THEP,PS1,PD1
DOUBLE PRECISION :: XN(4),YN(4),ZN(4)
!-----------------------------------------------------------------------------------------------!
!    Define Sides                                                                               !
!-----------------------------------------------------------------------------------------------!
RL1=(XX(2)-XX(1))**2+(YY(2)-YY(1))**2+(ZZ(2)-ZZ(1))**2
RL1=RL1+(XX(4)-XX(3))**2+(YY(4)-YY(3))**2+(ZZ(4)-ZZ(3))**2      
RL2=(XX(2)-XX(3))**2+(YY(2)-YY(3))**2+(ZZ(2)-ZZ(3))**2
RL2=RL2+(XX(4)-XX(1))**2+(YY(4)-YY(1))**2+(ZZ(4)-ZZ(1))**2
!-----------------------------------------------------------------------------------------------!
!    IF L1 > L2 Divide Along L1                                                                 !
!-----------------------------------------------------------------------------------------------!
PS=0.D0
PD=0.D0
!-----------------------------------------------------------------------------------------------!
IF (RL1 > RL2) THEN
   ASP=DSQRT(RL1/RL2)
   NP=ASP+0.01
!-----------------------------------------------------------------------------------------------!
!    Intermediate Points                                                                        !
!-----------------------------------------------------------------------------------------------!
   STEP=1.D0/NP
   THETA=0.D0
   DO I=1,NP
      XN(1)=XX(1)+THETA*(XX(2)-XX(1))
      THEP=THETA+STEP       
      XN(2)=XX(1)+THEP *(XX(2)-XX(1))
      XN(4)=XX(4)+THETA*(XX(3)-XX(4))
      XN(3)=XX(4)+THEP *(XX(3)-XX(4))
      YN(1)=YY(1)+THETA*(YY(2)-YY(1))
      YN(2)=YY(1)+THEP *(YY(2)-YY(1))
      YN(4)=YY(4)+THETA*(YY(3)-YY(4))
      YN(3)=YY(4)+THEP *(YY(3)-YY(4))
      ZN(1)=ZZ(1)+THETA*(ZZ(2)-ZZ(1))
      ZN(2)=ZZ(1)+THEP *(ZZ(2)-ZZ(1))
      ZN(4)=ZZ(4)+THETA*(ZZ(3)-ZZ(4))
      ZN(3)=ZZ(4)+THEP *(ZZ(3)-ZZ(4))
      THETA=THEP
!-----------------------------------------------------------------------------------------------!
!    Sub-Panels                                                                                 !
!-----------------------------------------------------------------------------------------------!
!     CALL POTANA_ADP(TOL,MMAX,XN,YN,ZN,X,Y,Z,PS1,PD1) !Analytical Integration
      CALL POTNUM_ADP(TOL,MMAX,XN,YN,ZN,X,Y,Z,PS1,PD1) !Numerical Integration
!-----------------------------------------------------------------------------------------------!
!    Source Potential (SUMS)                                                                    !
!-----------------------------------------------------------------------------------------------!
      PS=PS+PS1
!-----------------------------------------------------------------------------------------------!
!    Dipole Potential (SUMD)                                                                    !
!-----------------------------------------------------------------------------------------------!
      PD=PD+PD1
   END DO !I=1,NP
ELSE !(RL1 > RL2)
!-----------------------------------------------------------------------------------------------!
!    IF L2 >= L1 Divide Along L2                                                                !
!-----------------------------------------------------------------------------------------------!
   ASP=DSQRT(RL2/RL1)
   NP=ASP+0.01
!-----------------------------------------------------------------------------------------------!
!    Intermediate Points                                                                        !
!-----------------------------------------------------------------------------------------------!
   STEP=1.D0/NP
   THETA=0.D0
   DO I=1,NP
      XN(1)=XX(1)+THETA*(XX(4)-XX(1))
      THEP=THETA+STEP      
      XN(2)=XX(2)+THETA*(XX(3)-XX(2))
      XN(3)=XX(2)+THEP *(XX(3)-XX(2))
      XN(4)=XX(1)+THEP *(XX(4)-XX(1))
      YN(1)=YY(1)+THETA*(YY(4)-YY(1))
      YN(2)=YY(2)+THETA*(YY(3)-YY(2))
      YN(3)=YY(2)+THEP *(YY(3)-YY(2))
      YN(4)=YY(1)+THEP *(YY(4)-YY(1))
      ZN(1)=ZZ(1)+THETA*(ZZ(4)-ZZ(1))
      ZN(2)=ZZ(2)+THETA*(ZZ(3)-ZZ(2))
      ZN(3)=ZZ(2)+THEP *(ZZ(3)-ZZ(2))
      ZN(4)=ZZ(1)+THEP *(ZZ(4)-ZZ(1))
      THETA=THEP
!-----------------------------------------------------------------------------------------------!
!    Sub-Panels                                                                                 !
!-----------------------------------------------------------------------------------------------!
!     CALL POTANA_ADP(TOL,MMAX,XN,YN,ZN,X,Y,Z,PS1,PD1) !Analytical Integration
      CALL POTNUM_ADP(TOL,MMAX,XN,YN,ZN,X,Y,Z,PS1,PD1) !Numerical Integration
!-----------------------------------------------------------------------------------------------!
!    Source Potential (SUMS)                                                                    !
!-----------------------------------------------------------------------------------------------!
      PS=PS+PS1
!-----------------------------------------------------------------------------------------------!
!    Dipole Potential (SUMD)                                                                    !
!-----------------------------------------------------------------------------------------------!
      PD=PD+PD1
   END DO !I=1,NP
END IF !(RL1 > RL2)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE MSTRPAN_ADP
!===============================================================================================!
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
SUBROUTINE POTANA_ADP(TOL,MMAX,XX,YY,ZZ,X,Y,Z,PS,PD)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST, August 2012                                                  !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
INTEGER :: I,J,M,MMAX
LOGICAL :: OK
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4),X,Y,Z,PS,PD
DOUBLE PRECISION :: PRVPS,PRVPD,SUMS,SUMD,DELTAPS,DELTAPD,TOL
DOUBLE PRECISION :: CSI,CSIP,DCSI,ETA,ETAP,DETA
DOUBLE PRECISION :: XXX(4),YYY(4),ZZZ(4),X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0
!-----------------------------------------------------------------------------------------------!
OK=.TRUE.
PRVPS=1.D7
PRVPD=1.D7
M=1
DO WHILE (OK.AND.((M*M) <= MMAX))
!-----------------------------------------------------------------------------------------------!
!    Subdivision With M*M Panels                                                                !
!-----------------------------------------------------------------------------------------------!
   PS=0.D0
   PD=0.D0
!-----------------------------------------------------------------------------------------------!
   DCSI = 2.D0/M
   DETA = 2.D0/M
   ETA = -1.D0
!-----------------------------------------------------------------------------------------------!
   DO I=1,M
      CSI=-1.D0
      DO J=1,M
!-----------------------------------------------------------------------------------------------!
         CSIP=CSI+DCSI
         ETAP=ETA+DETA
         CALL PANCOOR(XX,YY,ZZ,CSI ,ETA ,XXX(1),YYY(1),ZZZ(1))
         CALL PANCOOR(XX,YY,ZZ,CSIP,ETA ,XXX(2),YYY(2),ZZZ(2))
         CALL PANCOOR(XX,YY,ZZ,CSIP,ETAP,XXX(3),YYY(3),ZZZ(3))
         CALL PANCOOR(XX,YY,ZZ,CSI ,ETAP,XXX(4),YYY(4),ZZZ(4))
!-----------------------------------------------------------------------------------------------!
!    Analytical Integration                                                                     !
!-----------------------------------------------------------------------------------------------!
         CALL PANEL(XXX,YYY,ZZZ,X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0)
         CALL POTPANH(XXX,YYY,ZZZ,X0,Y0,Z0,UNX0,UNY0,UNZ0,X,Y,Z,SUMS,SUMD)
         PS=PS+SUMS
         PD=PD+SUMD
!-----------------------------------------------------------------------------------------------!
         CSI=CSIP
      END DO !J=1,M
      ETA=ETAP
   END DO !I=1,M
!-----------------------------------------------------------------------------------------------!
!    Delta Source/Dipole Potential                                                              !
!-----------------------------------------------------------------------------------------------!
   DELTAPS=DABS(PS-PRVPS)/(1.D0+DABS(PS))
   DELTAPD=DABS(PD-PRVPD)/(1.D0+DABS(PD))
!-----------------------------------------------------------------------------------------------!
!  IF ((DELTAPS < TOL).AND.(DELTAPD < TOL)) OK=.FALSE.
   IF (DELTAPS < TOL) OK=.FALSE.
   PRVPS=PS
   PRVPD=PD
   M=3*M
!-----------------------------------------------------------------------------------------------!
END DO !(OK.AND.((M*M) <= MMAX))
!-----------------------------------------------------------------------------------------------!
!IF ((DELTAPS > TOL).OR.(DELTAPD > TOL)) THEN
IF (DELTAPS > TOL) THEN
   WRITE(10,*) 'POTANA_ADP=',X,Y,Z
   WRITE(10,*) 'TOL=',DELTAPS,DELTAPD
END IF !((TOLS > TOL).OR.(TOLD > TOL))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE POTANA_ADP
!===============================================================================================!
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
SUBROUTINE POTNUM_ADP(TOL,MMAX,XX,YY,ZZ,X,Y,Z,PS,PD)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST, January 2002                                                 !
!    Modified  : 21082012, J. Baltazar, new definition of total number of subdivisions          !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
INTEGER          :: I,M,MAUX1,MAUX2,MMAX
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4),XXX(4),YYY(4),ZZZ(4),X,Y,Z
DOUBLE PRECISION :: PS,PD,TOL,TOL1,TOLS,TOLD
DOUBLE PRECISION,ALLOCATABLE :: INFCSI(:),SUPCSI(:),INFETA(:),SUPETA(:)
DOUBLE PRECISION,ALLOCATABLE :: INFCSI1(:),SUPCSI1(:),INFETA1(:),SUPETA1(:)
DOUBLE PRECISION,ALLOCATABLE :: INTS1(:),INTS2(:),INTD1(:),INTD2(:)
DOUBLE PRECISION,ALLOCATABLE :: DELTAPS(:),DELTAPD(:)
LOGICAL          :: OKMASTER
LOGICAL,ALLOCATABLE          :: OK(:)
!-----------------------------------------------------------------------------------------------!
OKMASTER=.TRUE.
PS  =0.D0
PD  =0.D0
TOLS=0.D0
TOLD=0.D0
!-----------------------------------------------------------------------------------------------!
M=1
ALLOCATE(INFCSI(1),SUPCSI(1),INFETA(1),SUPETA(1))
INFCSI(1)=-1.D0
SUPCSI(1)= 1.D0
INFETA(1)=-1.D0
SUPETA(1)= 1.D0
!-----------------------------------------------------------------------------------------------!
DO WHILE (OKMASTER)
!-----------------------------------------------------------------------------------------------!
   ALLOCATE(INTS1(M),INTS2(M),INTD1(M),INTD2(M),OK(M),DELTAPS(M),DELTAPD(M))
   MAUX1=0
   MAUX2=0
!-----------------------------------------------------------------------------------------------!
   DO I=1,M
      CALL PANCOOR(XX,YY,ZZ,INFCSI(I),INFETA(I),XXX(1),YYY(1),ZZZ(1))
      CALL PANCOOR(XX,YY,ZZ,SUPCSI(I),INFETA(I),XXX(2),YYY(2),ZZZ(2))
      CALL PANCOOR(XX,YY,ZZ,SUPCSI(I),SUPETA(I),XXX(3),YYY(3),ZZZ(3))
      CALL PANCOOR(XX,YY,ZZ,INFCSI(I),SUPETA(I),XXX(4),YYY(4),ZZZ(4))
!-----------------------------------------------------------------------------------------------!
      CALL GAUSSPAN(4,XXX,YYY,ZZZ,X,Y,Z,INTS1(I),INTD1(I))
      CALL GAUSSPAN(8,XXX,YYY,ZZZ,X,Y,Z,INTS2(I),INTD2(I))
!-----------------------------------------------------------------------------------------------!
      DELTAPS(I)=DABS(INTS2(I)-INTS1(I))
      DELTAPD(I)=DABS(INTD2(I)-INTD1(I))
      TOL1=TOL*DABS(SUPCSI(I)-INFCSI(I))*DABS(SUPETA(I)-INFETA(I))/4.D0
!-----------------------------------------------------------------------------------------------!
!     IF ((DELTAPS(I) < TOL1).AND.(DELTAPD(I) < TOL1)) THEN
      IF (DELTAPS(I) < TOL1) THEN
         MAUX2=MAUX2+1
         OK(I)=.FALSE.
      ELSE !((DELTAPS(I) < TOL1).AND.(DELTAPD(I) < TOL1))
         MAUX1=MAUX1+4
         OK(I)=.TRUE.
      END IF !((DELTAPS(I) < TOL1).AND.(DELTAPD(I) < TOL1))
   END DO !I=1,M
!-----------------------------------------------------------------------------------------------!
   OKMASTER=.FALSE.
   ALLOCATE(INFCSI1(MAUX1),SUPCSI1(MAUX1),INFETA1(MAUX1),SUPETA1(MAUX1))
!-----------------------------------------------------------------------------------------------!
   MAUX1=0
   DO I=1,M
      IF (OK(I)) THEN
         MAUX1=MAUX1+1
         INFCSI1(MAUX1)=INFCSI(I)
         SUPCSI1(MAUX1)=(INFCSI(I)+SUPCSI(I))*0.5D0
         INFETA1(MAUX1)=INFETA(I)
         SUPETA1(MAUX1)=(INFETA(I)+SUPETA(I))*0.5D0
!-----------------------------------------------------------------------------------------------!
         MAUX1=MAUX1+1
         INFCSI1(MAUX1)=(INFCSI(I)+SUPCSI(I))*0.5D0
         SUPCSI1(MAUX1)=SUPCSI(I)
         INFETA1(MAUX1)=INFETA(I)
         SUPETA1(MAUX1)=(INFETA(I)+SUPETA(I))*0.5D0
!-----------------------------------------------------------------------------------------------!
         MAUX1=MAUX1+1
         INFCSI1(MAUX1)=(INFCSI(I)+SUPCSI(I))*0.5D0
         SUPCSI1(MAUX1)=SUPCSI(I)
         INFETA1(MAUX1)=(INFETA(I)+SUPETA(I))*0.5D0
         SUPETA1(MAUX1)=SUPETA(I)
!-----------------------------------------------------------------------------------------------!
         MAUX1=MAUX1+1
         INFCSI1(MAUX1)=INFCSI(I)
         SUPCSI1(MAUX1)=(INFCSI(I)+SUPCSI(I))*0.5D0
         INFETA1(MAUX1)=(INFETA(I)+SUPETA(I))*0.5D0
         SUPETA1(MAUX1)=SUPETA(I)
         OKMASTER=.TRUE.
      ELSE !(OK(I))
         PS  =PS  +INTS2(I)
         PD  =PD  +INTD2(I)
         TOLS=TOLS+DELTAPS(I)
         TOLD=TOLD+DELTAPD(I)
      END IF !(OK(I))
   END DO !I=1,M
!-----------------------------------------------------------------------------------------------!
   IF ((MAUX1+MAUX2) > MMAX) THEN
      OKMASTER=.FALSE.
      DO I=1,M
         IF (OK(I)) THEN
            PS  =PS  +INTS2(I)
            PD  =PD  +INTD2(I)
            TOLS=TOLS+DELTAPS(I)
            TOLD=TOLD+DELTAPD(I)
         END IF !(OK(I))
      END DO !I=1,M
   ELSE !(M > MMAX)
      DEALLOCATE(INFCSI,SUPCSI,INFETA,SUPETA)
      ALLOCATE(INFCSI(MAUX1),SUPCSI(MAUX1),INFETA(MAUX1),SUPETA(MAUX1))
      INFCSI(:)=INFCSI1(:)
      SUPCSI(:)=SUPCSI1(:)
      INFETA(:)=INFETA1(:)
      SUPETA(:)=SUPETA1(:)
      DEALLOCATE(INFCSI1,SUPCSI1,INFETA1,SUPETA1)
   END IF !(M > MMAX)
!-----------------------------------------------------------------------------------------------!
   DEALLOCATE(INTS1,INTS2,INTD1,INTD2,OK,DELTAPS,DELTAPD)
   M=MAUX1
END DO !(OKMASTER)
!-----------------------------------------------------------------------------------------------!
!IF ((TOLS > TOL).OR.(TOLD > TOL)) THEN
IF (TOLS > TOL) THEN
   WRITE(10,*) 'POTNUM_ADP=',X,Y,Z
   WRITE(10,*) 'TOL=',TOLS,TOLD
END IF !((TOLS > TOL).OR.(TOLD > TOL))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE POTNUM_ADP
!===============================================================================================!
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
SUBROUTINE FARFIELD(X0,Y0,Z0,UN0X,UN0Y,UN0Z,A0,X,Y,Z,PS,PD)
!-----------------------------------------------------------------------------------------------!
!    INPUT:                                                                                     !
!                                                                                               !
!    X0,Y0,Z0         :      CARTESIAN COORDINATES OF PANEL CENTER POINT                        !
!                                                                                               !
!    UN0X,UN0Y,UN0Z   :      COMPONENTS OF UNIT NORMAL AT THE PANEL CENTER POINT                !
!                                                                                               !
!    X,Y,Z            :      CARTESIAN COORDINATES OF FIELD POINT                               !
!                                                                                               !
!    A0               :      PANEL AREA                                                         !
!                                                                                               !
!    OUTPUT:                                                                                    !
!                                                                                               !
!    PS               :      POTENTIAL DUE TO THE SOURCE DISTRIBUTION                           !
!                                                                                               !
!    PD               :      POTENTIAL DUE TO THE DIPOLE DISTRIBUTION                           !
!-----------------------------------------------------------------------------------------------!
!    Declarations                                                                               !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
DOUBLE PRECISION :: X0,Y0,Z0,UN0X,UN0Y,UN0Z,X,Y,Z,A0,PS,PD
DOUBLE PRECISION :: RX,RY,RZ,R0,RN
!-----------------------------------------------------------------------------------------------!
!    Position Vector to the Center Point                                                        !
!-----------------------------------------------------------------------------------------------!
RX=X-X0
RY=Y-Y0
RZ=Z-Z0
R0=DSQRT(RX*RX+RY*RY+RZ*RZ)
!-----------------------------------------------------------------------------------------------!
!    Distance in the Direction of the Normal                                                    !
!-----------------------------------------------------------------------------------------------!
RN=RX*UN0X+RY*UN0Y+RZ*UN0Z
!-----------------------------------------------------------------------------------------------!
!    Source Potential                                                                           !
!-----------------------------------------------------------------------------------------------!
PS=A0/R0
!-----------------------------------------------------------------------------------------------!
!    Dipole Potential                                                                           !
!-----------------------------------------------------------------------------------------------!
PD=RN*A0/(R0**3)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE FARFIELD
!===============================================================================================!
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
SUBROUTINE GAUSSPAN(N,XX,YY,ZZ,X,Y,Z,PS,PD)
!-----------------------------------------------------------------------------------------------!
!    INPUT:                                                                                     !
!                                                                                               !
!    N                : ORDER OF GAUSSIAN RULE FOR PANEL INTEGRATION                            !
!                                                                                               !
!    XX(4),YY(4),ZZ(4): CARTESIAN COORDINATES OF PANEL CORNER POINTS                            !
!                                                                                               !
!    X,Y,Z            : CARTESIAN COORDINATES OF FIELD POINT                                    !
!                                                                                               !
!    OUTPUT:                                                                                    !
!                                                                                               !
!    PS               : POTENTIAL DUE TO THE SOURCE DISTRIBUTION                                !
!                                                                                               !
!    PD               : POTENTIAL DUE TO THE DIPOLE DISTRIBUTION                                !
!-----------------------------------------------------------------------------------------------!
INTEGER :: I,J,N
DOUBLE PRECISION :: XG(N),W(N)
DOUBLE PRECISION :: CSI,ETA
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4),X,Y,Z,PS,PD
DOUBLE PRECISION :: RM,A1X,A1Y,A1Z,A2X,A2Y,A2Z,VNX,VNY,VNZ,VNM
DOUBLE PRECISION :: SUMS,SUMD,XQ,YQ,ZQ
!-----------------------------------------------------------------------------------------------!
!    Computations                                                                               !
!-----------------------------------------------------------------------------------------------!
CALL GAUSS2M(N,XG,W)
!-----------------------------------------------------------------------------------------------!
!    Loop on the Element                                                                        !
!-----------------------------------------------------------------------------------------------!
SUMS=0.D0
SUMD=0.D0
DO I=1,N
   CSI=XG(I)
   DO J=1,N
      ETA=XG(J)
      CALL PANCOOR(XX,YY,ZZ,CSI,ETA,XQ,YQ,ZQ)
      CALL PANVECTORS(XX,YY,ZZ,CSI,ETA,A1X,A1Y,A1Z,A2X,A2Y,A2Z)
      VNX=A1Y*A2Z-A1Z*A2Y
      VNY=A1Z*A2X-A1X*A2Z
      VNZ=A1X*A2Y-A1Y*A2X
      VNM=DSQRT(VNX*VNX+VNY*VNY+VNZ*VNZ)
      RM=DSQRT((X-XQ)**2+(Y-YQ)**2+(Z-ZQ)**2)
      FS=VNM/RM
      FD=(VNX*(X-XQ)+VNY*(Y-YQ)+VNZ*(Z-ZQ))/RM**3
      SUMS=SUMS+W(I)*W(J)*FS
      SUMD=SUMD+W(I)*W(J)*FD
   END DO !J=1,N
END DO !I=1,N
!-----------------------------------------------------------------------------------------------!
PS=SUMS
PD=SUMD
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE GAUSSPAN
!===============================================================================================!
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
SUBROUTINE GAUSSPAN3(N,XX,YY,ZZ,X,Y,Z,PS,PD)
!-----------------------------------------------------------------------------------------------!
!    INPUT:                                                                                     !
!                                                                                               !
!    N                :     ORDER OF GAUSSIAN RULE FOR PANEL INTEGRATION                        !
!                                                                                               !
!    XX(3),YY(3),ZZ(3):     CARTESIAN COORDINATES OF PANEL CORNER POINTS                        !
!                                                                                               !
!    X,Y,Z            :     CARTESIAN COORDINATES OF FIELD POINT                                !
!                                                                                               !
!    OUTPUT:                                                                                    !
!                                                                                               !
!    PS               :     POTENTIAL DUE TO THE SOURCE DISTRIBUTION                            !
!                                                                                               !
!    PD               :     POTENTIAL DUE TO THE DIPOLE DISTRIBUTION                            !
!-----------------------------------------------------------------------------------------------!
!    Declarations                                                                               !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
INTEGER :: I,J,N
DOUBLE PRECISION :: CSI,ETA
DOUBLE PRECISION :: XX(3),YY(3),ZZ(3),X,Y,Z,PS,PD
DOUBLE PRECISION :: RM,A1X,A1Y,A1Z,A2X,A2Y,A2Z,VNX,VNY,VNZ,VNM
DOUBLE PRECISION :: SUMS,SUMD,XQ,YQ,ZQ,FD,FS
DOUBLE PRECISION :: CSIG(N),ETAG(N),ZETAG(N),W(N)
!-----------------------------------------------------------------------------------------------!
!    Computations                                                                               !
!-----------------------------------------------------------------------------------------------!
CALL GAUSS2M3(N,CSIG,ETAG,ZETAG,W)
!-----------------------------------------------------------------------------------------------!
!    Loop on the Element                                                                        !
!-----------------------------------------------------------------------------------------------!
SUMS=0.D0
SUMD=0.D0
DO J=1,N
   CALL PANCOOR3(XX,YY,ZZ,CSIG(J),ETAG(J),XQ,YQ,ZQ)
   A1X=XX(2)-XX(1)
   A1Y=YY(2)-YY(1)
   A1Z=ZZ(2)-ZZ(1)
   A2X=XX(3)-XX(1)
   A2Y=YY(3)-YY(1)
   A2Z=ZZ(3)-ZZ(1)
   VNX=A1Y*A2Z-A1Z*A2Y
   VNY=A1Z*A2X-A1X*A2Z
   VNZ=A1X*A2Y-A1Y*A2X
   VNM=DSQRT(VNX*VNX+VNY*VNY+VNZ*VNZ)
   RM=DSQRT((X-XQ)**2+(Y-YQ)**2+(Z-ZQ)**2)
   FS=VNM/RM
   FD=(VNX*(X-XQ)+VNY*(Y-YQ)+VNZ*(Z-ZQ))/RM**3
   SUMS=SUMS+W(J)*FS
   SUMD=SUMD+W(J)*FD
END DO !J=1,N
PS=SUMS
PD=SUMD
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE GAUSSPAN3
!-----------------------------------------------------------------------------------------------!
