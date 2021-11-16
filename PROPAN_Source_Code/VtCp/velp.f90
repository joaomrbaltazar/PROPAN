!-----------------------------------------------------------------------------------------------!
!    Calculation of the tangential perturbation velocities on the blade                         !
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
SUBROUTINE VELP(TT)
!-----------------------------------------------------------------------------------------------!
!    Created by: J.A.C. Falcao de Campos, IST                                                   !
!    Modified  : 14112013, J. Baltazar, version 1.0                                             !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,TT
!-----------------------------------------------------------------------------------------------!
DO J=1,NRP
   DO I=1,NCP
!-----------------------------------------------------------------------------------------------!
!    u1 Derivative                                                                              !
!-----------------------------------------------------------------------------------------------!
      IF (I == 1) THEN
         VT1P(I,J)=A1P(I,J)*POTP(1    ,J,TT)+ &
                   B1P(I,J)*POTP(2    ,J,TT)+ &
                   C1P(I,J)*POTP(3    ,J,TT)
      END IF !(I == 1)
      IF (I == NCP) THEN
         VT1P(I,J)=A1P(I,J)*POTP(NCP-2,J,TT)+ &
                   B1P(I,J)*POTP(NCP-1,J,TT)+ &
                   C1P(I,J)*POTP(NCP  ,J,TT)
      END IF !(I == NCP)
      IF ((I /= 1).AND.(I /= NCP)) THEN
         VT1P(I,J)=A1P(I,J)*POTP(I-1  ,J,TT)+ &
                   B1P(I,J)*POTP(I    ,J,TT)+ &
                   C1P(I,J)*POTP(I+1  ,J,TT)
      END IF !((I /= 1).AND.(I /= NCP))
!-----------------------------------------------------------------------------------------------!
      IF (J == NRP.AND.ISTRIP == 2) THEN
         IF (I == 1.OR.I == NC+1) THEN
            VT1P(I,J)=A1P(I,J)*POTP(I  ,J,TT)+ &
                      B1P(I,J)*POTP(I+1,J,TT)+ &
                      C1P(I,J)*POTP(I+2,J,TT)
         ELSEIF (I == NC.OR.I == NCP) THEN
            VT1P(I,J)=A1P(I,J)*POTP(I-2,J,TT)+ &
                      B1P(I,J)*POTP(I-1,J,TT)+ &
                      C1P(I,J)*POTP(I  ,J,TT)
         ELSE
            VT1P(I,J)=A1P(I,J)*POTP(I-1,J,TT)+ &
                      B1P(I,J)*POTP(I  ,J,TT)+ &
                      C1P(I,J)*POTP(I+1,J,TT)
         END IF
      END IF !(J == NRP.AND.ISTRIP == 2)
!-----------------------------------------------------------------------------------------------!
!    u2 Derivative                                                                              !
!-----------------------------------------------------------------------------------------------!
      IF (J == 1) THEN
         IF (I == 1) THEN
            VT2P(I,J)=A2P (I,J)*POTP(1    ,1,TT)+ &
                      B2P (I,J)*POTP(1    ,2,TT)+ &
                      C2P (I,J)*POTP(1    ,3,TT)+ &
                      D2P (I,J)*POTP(2    ,1,TT)+ &
                      E2P (I,J)*POTP(3    ,1,TT)
            VT2S(I,J)=AA2P(I,J)*POTP(1    ,1,TT)+ &
                      BB2P(I,J)*POTP(1    ,2,TT)+ &
                      CC2P(I,J)*POTP(1    ,3,TT)
         END IF !(I == 1)
         IF (I == NCP) THEN
            VT2P(I,J)=A2P (I,J)*POTP(NCP  ,1,TT)+ &
                      B2P (I,J)*POTP(NCP  ,2,TT)+ &
                      C2P (I,J)*POTP(NCP  ,3,TT)+ &
                      D2P (I,J)*POTP(NCP-1,1,TT)+ &
                      E2P (I,J)*POTP(NCP-2,1,TT)
            VT2S(I,J)=AA2P(I,J)*POTP(NCP  ,1,TT)+ &
                      BB2P(I,J)*POTP(NCP  ,2,TT)+ &
                      CC2P(I,J)*POTP(NCP  ,3,TT)
         END IF !(I == NCP)
         IF ((I /= 1).AND.(I /= NCP)) THEN
            VT2P(I,J)=A2P (I,J)*POTP(I    ,1,TT)+ &
                      B2P (I,J)*POTP(I    ,2,TT)+ &
                      C2P (I,J)*POTP(I    ,3,TT)+ &
                      D2P (I,J)*POTP(I-1  ,1,TT)+ &
                      E2P (I,J)*POTP(I+1  ,1,TT)
            VT2S(I,J)=AA2P(I,J)*POTP(I    ,1,TT)+ &
                      BB2P(I,J)*POTP(I    ,2,TT)+ &
                      CC2P(I,J)*POTP(I    ,3,TT)
         END IF !((I /= 1).AND.(I /= NCP))
      END IF !(J == 1)
!-----------------------------------------------------------------------------------------------!
      IF (J == NRP) THEN
         IF (I == 1) THEN
            VT2P(I,J)=A2P (I,J)*POTP(1    ,NRP-2,TT)+ &
                      B2P (I,J)*POTP(1    ,NRP-1,TT)+ &
                      C2P (I,J)*POTP(1    ,NRP  ,TT)+ &
                      D2P (I,J)*POTP(2    ,NRP  ,TT)+ &
                      E2P (I,J)*POTP(3    ,NRP  ,TT)
            VT2S(I,J)=AA2P(I,J)*POTP(1    ,NRP-2,TT)+ &
                      BB2P(I,J)*POTP(1    ,NRP-1,TT)+ &
                      CC2P(I,J)*POTP(1    ,NRP  ,TT)
         END IF !(I == 1)
         IF (I == NCP) THEN
            VT2P(I,J)=A2P (I,J)*POTP(NCP  ,NRP-2,TT)+ &
                      B2P (I,J)*POTP(NCP  ,NRP-1,TT)+ &
                      C2P (I,J)*POTP(NCP  ,NRP  ,TT)+ &
                      D2P (I,J)*POTP(NCP-1,NRP  ,TT)+ &
                      E2P (I,J)*POTP(NCP-2,NRP  ,TT)
            VT2S(I,J)=AA2P(I,J)*POTP(NCP  ,NRP-2,TT)+ &
                      BB2P(I,J)*POTP(NCP  ,NRP-1,TT)+ &
                      CC2P(I,J)*POTP(NCP  ,NRP  ,TT)
         END IF !(I == NCP)
         IF ((I /= 1).AND.(I /= NCP)) THEN
            VT2P(I,J)=A2P (I,J)*POTP(I    ,NRP-2,TT)+ &
                      B2P (I,J)*POTP(I    ,NRP-1,TT)+ &
                      C2P (I,J)*POTP(I    ,NRP  ,TT)+ &
                      D2P (I,J)*POTP(I-1  ,NRP  ,TT)+ &
                      E2P (I,J)*POTP(I+1  ,NRP  ,TT)
            VT2S(I,J)=AA2P(I,J)*POTP(I    ,NRP-2,TT)+ &
                      BB2P(I,J)*POTP(I    ,NRP-1,TT)+ &
                      CC2P(I,J)*POTP(I    ,NRP  ,TT)
         END IF !((I /= 1).AND.(I /= NCP))
      END IF !(J == NRP)
!-----------------------------------------------------------------------------------------------!
      IF ((J > 1).AND.(J < NRP)) THEN
         IF (I == 1) THEN
            VT2P(I,J)=A2P (I,J)*POTP(1  ,J-1,TT)+ &
                      B2P (I,J)*POTP(1  ,J  ,TT)+ &
                      C2P (I,J)*POTP(1  ,J+1,TT)+ &
                      D2P (I,J)*POTP(2  ,J  ,TT)+ &
                      E2P (I,J)*POTP(3  ,J  ,TT)
            VT2S(I,J)=AA2P(I,J)*POTP(1  ,J-1,TT)+ &
                      BB2P(I,J)*POTP(1  ,J  ,TT)+ &
                      CC2P(I,J)*POTP(1  ,J+1,TT)
         END IF !(I == 1)
         IF (I == NCP) THEN
            VT2P(I,J)=A2P (I,J)*POTP(NCP,J-1,TT)+ &
                      B2P (I,J)*POTP(NCP,J  ,TT)+ &
                      C2P (I,J)*POTP(NCP,J+1,TT)+ &
                      D2P (I,J)*POTP(NCP-1,J,TT)+ &
                      E2P (I,J)*POTP(NCP-2,J,TT)
            VT2S(I,J)=AA2P(I,J)*POTP(NCP,J-1,TT)+ &
                      BB2P(I,J)*POTP(NCP,J  ,TT)+ &
                      CC2P(I,J)*POTP(NCP,J+1,TT)
         END IF !(I == NCP)
         IF ((I /= 1).AND.(I /= NCP)) THEN
            VT2P(I,J)=A2P (I,J)*POTP(I  ,J-1,TT)+ &
                      B2P (I,J)*POTP(I  ,J  ,TT)+ &
                      C2P (I,J)*POTP(I  ,J+1,TT)+ &
                      D2P (I,J)*POTP(I-1,J  ,TT)+ &
                      E2P (I,J)*POTP(I+1,J  ,TT)
            VT2S(I,J)=AA2P(I,J)*POTP(I  ,J-1,TT)+ &
                      BB2P(I,J)*POTP(I  ,J  ,TT)+ &
                      CC2P(I,J)*POTP(I  ,J+1,TT)
         END IF !((I /= 1).AND.(I /= NCP))
      END IF !((J > 1).AND.(J < NRP))
!-----------------------------------------------------------------------------------------------!
      IF (J == NRP.AND.ISTRIP == 2) THEN
         IF (I == 1.OR.I == NC+1) THEN
            VT2P(I,J)=B2P (I,J)*POTP(I      ,J,TT)+ &
                      C2P (I,J)*POTP(NCP-I+1,J,TT)+ &
                      D2P (I,J)*POTP(I+1    ,J,TT)+ &
                      E2P (I,J)*POTP(I+2    ,J,TT)
            VT2S(I,J)=BB2P(I,J)*POTP(I      ,J,TT)+ &
                      CC2P(I,J)*POTP(NCP-I+1,J,TT)
         ELSEIF (I == NC.OR.I == NCP) THEN
            VT2P(I,J)=B2P (I,J)*POTP(I      ,J,TT)+ &
                      C2P (I,J)*POTP(NCP-I+1,J,TT)+ &
                      D2P (I,J)*POTP(I-1    ,J,TT)+ &
                      E2P (I,J)*POTP(I-2    ,J,TT)
            VT2S(I,J)=BB2P(I,J)*POTP(I      ,J,TT)+ &
                      CC2P(I,J)*POTP(NCP-I+1,J,TT)
         ELSE
            VT2P(I,J)=B2P (I,J)*POTP(I      ,J,TT)+ &
                      C2P (I,J)*POTP(NCP-I+1,J,TT)+ &
                      D2P (I,J)*POTP(I-1    ,J,TT)+ &
                      E2P (I,J)*POTP(I+1    ,J,TT)
            VT2S(I,J)=BB2P(I,J)*POTP(I      ,J,TT)+ &
                      CC2P(I,J)*POTP(NCP-I+1,J,TT)
         END IF
      END IF !(J == NRP.AND.ISTRIP == 2)
!-----------------------------------------------------------------------------------------------!
   END DO !I=1,NCP
END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE VELP
!-----------------------------------------------------------------------------------------------!
