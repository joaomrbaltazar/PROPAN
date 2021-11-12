!-----------------------------------------------------------------------------------------------!
!    Computation of potentials for a subdivised panel                                           !
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
SUBROUTINE POTLINSUB(TOL,MMAX,XX,YY,ZZ,X,Y,Z,PDL,PDR)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 06112013, J. Baltazar, version 1.0                                             !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
LOGICAL :: OK
INTEGER :: I,MM,MMAX
DOUBLE PRECISION :: STEP,THEP,THETA,XN(4),YN(4),ZN(4)
DOUBLE PRECISION :: XI,YI,ZI,XF,YF,ZF,LP,LI
DOUBLE PRECISION :: TOL,PRVPDL,PRVPDR,ERRPDL,ERRPDR
DOUBLE PRECISION :: X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0,PS,PD
DOUBLE PRECISION :: XX(4),YY(4),ZZ(4),X,Y,Z,PDL,PDR
!-----------------------------------------------------------------------------------------------!
!    Calculation of Mean Panel Length                                                           !
!-----------------------------------------------------------------------------------------------!
XI=(XX(1)+XX(4))*0.5D0
YI=(YY(1)+YY(4))*0.5D0
ZI=(ZZ(1)+ZZ(4))*0.5D0
!-----------------------------------------------------------------------------------------------!
XF=(XX(2)+XX(3))*0.5D0
YF=(YY(2)+YY(3))*0.5D0
ZF=(ZZ(2)+ZZ(3))*0.5D0
!-----------------------------------------------------------------------------------------------!
LP=DSQRT((XF-XI)*(XF-XI)+(YF-YI)*(YF-YI)+(ZF-ZI)*(ZF-ZI))
!-----------------------------------------------------------------------------------------------!
MM=0
OK=.TRUE.
PRVPDL=1.D10
PRVPDR=1.D10
!-----------------------------------------------------------------------------------------------!
DO WHILE (OK.AND.(MM <= MMAX))
!-----------------------------------------------------------------------------------------------!
   MM=MM+1
   PDL=0.D0
   PDR=0.D0
!-----------------------------------------------------------------------------------------------!
!    Sub-Panel Length                                                                           !
!-----------------------------------------------------------------------------------------------!
   LI=LP/DFLOAT(MM)
!-----------------------------------------------------------------------------------------------!
!    Intermediate Points                                                                        !
!-----------------------------------------------------------------------------------------------!
   STEP =1.D0/DFLOAT(MM)
   THETA=0.D0
   DO I=1,MM
      THEP =THETA+STEP       
      XN(1)=XX(1)+THETA*(XX(2)-XX(1))
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
!    Potential Calculation on each Sub-Panel                                                    !
!-----------------------------------------------------------------------------------------------!
      CALL PANEL(XN,YN,ZN,X0,Y0,Z0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX0,UNY0,UNZ0,A0)
      CALL POTPANH(XN,YN,ZN,X0,Y0,Z0,UNX0,UNY0,UNZ0,X,Y,Z,PS,PD)
!-----------------------------------------------------------------------------------------------!
      PDL=PDL+PD*(1.D0-(DFLOAT(I)-0.5D0)*LI/LP)
      PDR=PDR+PD*((DFLOAT(I)-0.5D0)*LI/LP)
!-----------------------------------------------------------------------------------------------!
   END DO !I=1,MM
!-----------------------------------------------------------------------------------------------!
   ERRPDL=DABS(PDL-PRVPDL)/(DABS(PDL)+1.D0)
   ERRPDR=DABS(PDR-PRVPDR)/(DABS(PDR)+1.D0)
!-----------------------------------------------------------------------------------------------!
   IF ((ERRPDL < TOL).AND.(ERRPDR < TOL)) OK=.FALSE.
!-----------------------------------------------------------------------------------------------!
   PRVPDL=PDL
   PRVPDR=PDR
!-----------------------------------------------------------------------------------------------!
END DO !(OK.AND.(MM <= MMAX))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE POTLINSUB
!-----------------------------------------------------------------------------------------------!
