!-----------------------------------------------------------------------------------------------!
!    Pressure-jump at nozzle trailing-edge                                                      !
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
SUBROUTINE PRESNTE(TT,DCNTE)
!-----------------------------------------------------------------------------------------------!
!    Created by: 06112014, J. Baltazar, Interpolation scheme at nozzle t.e.                     !
!    Modified  : 02102018, J. Baltazar, Different duct Kutta points at inner and outer sides    !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,IITE,TT
DOUBLE PRECISION :: CPNI,CPNO,CPNNI,CPNNO,DCNTE(NNTP)
!-----------------------------------------------------------------------------------------------!
!    Pressure-Jump at XNTE                                                                      !
!-----------------------------------------------------------------------------------------------!
IF ((XITE > TOL).OR.(XOTE > TOL)) THEN
   DO J=1,NNTP
      CALL LININT(NNX,XN0(    1:NNX : 1,J),CPN (    1:NNX : 1,J,TT),1,XITE,CPNI )
      CALL LININT(NNX,XN0(NNXT1:NNX1:-1,J),CPN (NNXT1:NNX1:-1,J,TT),1,XOTE,CPNO )
      CALL LININT(NNX,XN0(    1:NNX : 1,J),CPNN(    1:NNX : 1,J,TT),1,XITE,CPNNI)
      CALL LININT(NNX,XN0(NNXT1:NNX1:-1,J),CPNN(NNXT1:NNX1:-1,J,TT),1,XOTE,CPNNO)
      DCNTE(J)=CPNI-CPNO
!-----------------------------------------------------------------------------------------------!
!    Pressure at the Inner Side                                                                 !
!-----------------------------------------------------------------------------------------------!
      I=1
      DO WHILE ((XN0(I,J) < XITE).AND.(I <= NNX))
         I=I+1
      END DO !((XN0(I,J) < XITE).AND.(I <= NNX))
      IITE=I
      DO I=IITE,NNX
         CPN (I,J,TT)=CPNI
         CPNN(I,J,TT)=CPNNI
      END DO !I=IITE,NNX
!-----------------------------------------------------------------------------------------------!
!    Pressure at the Outer Side                                                                 !
!-----------------------------------------------------------------------------------------------!
      I=NNXT1
      DO WHILE ((XN0(I,J) < XOTE).AND.(I >= NNX1))
         I=I-1
      END DO !((XN0(I,J) < XOTE).AND.(I >= NNX1))
      IITE=I
      DO I=IITE,NNX1,-1
         CPN (I,J,TT)=CPNO
         CPNN(I,J,TT)=CPNNO
      END DO !I=IITE,NNX
   END DO !J=1,NNTP
!-----------------------------------------------------------------------------------------------!
!    Pressure-Jump at Kutta Points                                                              !
!-----------------------------------------------------------------------------------------------!
ELSE !((XITE > TOL).OR.(XOTE > TOL))
   DCNTE(:)=CPN(NNX,:,TT)-CPN(NNX1,:,TT)
END IF !((XITE > TOL).OR.(XOTE > TOL))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PRESNTE
!-----------------------------------------------------------------------------------------------!
