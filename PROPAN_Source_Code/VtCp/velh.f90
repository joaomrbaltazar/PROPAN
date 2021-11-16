!-----------------------------------------------------------------------------------------------!
!    Calculation of the tangential perturbation velocities on the hub                           !
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
SUBROUTINE VELH(TT)
!-----------------------------------------------------------------------------------------------!
!    Created by: J.A.C. Falcao de Campos, IST                                                   !
!    Modified  : 15112013, J. Baltazar, version 1.0                                             !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,TT
!-----------------------------------------------------------------------------------------------!
DO J=1,NHTP
   DO I=1,NHX
!-----------------------------------------------------------------------------------------------!
!    u1 Derivative                                                                              !
!-----------------------------------------------------------------------------------------------!
      IF (I == 1) THEN
         VT1H(I,J)=A1H(I,J)*POTH(1    ,J,TT)+ &
                   B1H(I,J)*POTH(2    ,J,TT)+ &
                   C1H(I,J)*POTH(3    ,J,TT)
      END IF !(I == 1)
      IF (I == NHX) THEN
         VT1H(I,J)=A1H(I,J)*POTH(NHX-2,J,TT)+ &
                   B1H(I,J)*POTH(NHX-1,J,TT)+ &
                   C1H(I,J)*POTH(NHX  ,J,TT)
      END IF !(I == NHX)
      IF ((I /= 1).AND.(I /= NHX)) THEN
         VT1H(I,J)=A1H(I,J)*POTH(I-1  ,J,TT)+ &
                   B1H(I,J)*POTH(I    ,J,TT)+ &
                   C1H(I,J)*POTH(I+1  ,J,TT)
      END IF !((I /= 1).AND.(I /= NHX))
!-----------------------------------------------------------------------------------------------!
!    u2 Derivative                                                                              !
!-----------------------------------------------------------------------------------------------!
!    Periodic Boundary Condition                                                                !
!-----------------------------------------------------------------------------------------------!
      IF (J == 1) THEN
         IF (I == 1) THEN
            VT2H(I,J)=A2H(I,J)*POTH(1    ,NHTP,TT)+ &
                      B2H(I,J)*POTH(1    ,1   ,TT)+ &
                      C2H(I,J)*POTH(1    ,2   ,TT)+ &
                      D2H(I,J)*POTH(2    ,1   ,TT)+ &
                      E2H(I,J)*POTH(3    ,1   ,TT)
         END IF !(I == 1)
         IF (I == NHX) THEN
            VT2H(I,J)=A2H(I,J)*POTH(NHX  ,NHTP,TT)+ &
                      B2H(I,J)*POTH(NHX  ,1   ,TT)+ &
                      C2H(I,J)*POTH(NHX  ,2   ,TT)+ &
                      D2H(I,J)*POTH(NHX-1,1   ,TT)+ &
                      E2H(I,J)*POTH(NHX-2,1   ,TT)
         END IF !(I == NHX)
         IF ((I /= 1).AND.(I /= NHX)) THEN
            VT2H(I,J)=A2H(I,J)*POTH(I    ,NHTP,TT)+ &
                      B2H(I,J)*POTH(I    ,1   ,TT)+ &
                      C2H(I,J)*POTH(I    ,2   ,TT)+ &
                      D2H(I,J)*POTH(I-1  ,1   ,TT)+ &
                      E2H(I,J)*POTH(I+1  ,1   ,TT)
         END IF !((I /= 1).AND.(I /= NHX))
      END IF !(J == 1)
!-----------------------------------------------------------------------------------------------!
!    Periodic Boundary Condition                                                                !
!-----------------------------------------------------------------------------------------------!
      IF (J == NHTP) THEN
         IF (I == 1) THEN
            VT2H(I,J)=A2H(I,J)*POTH(1    ,NHTP-1,TT)+ &
                      B2H(I,J)*POTH(1    ,NHTP  ,TT)+ &
                      C2H(I,J)*POTH(1    ,1     ,TT)+ &
                      D2H(I,J)*POTH(2    ,NHTP  ,TT)+ &
                      E2H(I,J)*POTH(3    ,NHTP  ,TT)
         END IF !(I == 1)
         IF (I == NHX) THEN
            VT2H(I,J)=A2H(I,J)*POTH(NHX  ,NHTP-1,TT)+ &
                      B2H(I,J)*POTH(NHX  ,NHTP  ,TT)+ &
                      C2H(I,J)*POTH(NHX  ,1     ,TT)+ &
                      D2H(I,J)*POTH(NHX-1,NHTP  ,TT)+ &
                      E2H(I,J)*POTH(NHX-2,NHTP  ,TT)
         END IF !(I == NHX)
         IF ((I /= 1).AND.(I /= NHX)) THEN
            VT2H(I,J)=A2H(I,J)*POTH(I    ,NHTP-1,TT)+ &
                      B2H(I,J)*POTH(I    ,NHTP  ,TT)+ &
                      C2H(I,J)*POTH(I    ,1     ,TT)+ &
                      D2H(I,J)*POTH(I-1  ,NHTP  ,TT)+ &
                      E2H(I,J)*POTH(I+1  ,NHTP  ,TT)
         END IF !((I /= 1).AND.(I /= NHX))
      END IF !(J == NHTP)
!-----------------------------------------------------------------------------------------------!
      IF (J == NHT) THEN
         IF (I == 1) THEN
            VT2H(I,J)=A2H(I,J)*POTH(1    ,NHT-1,TT)+ &
                      B2H(I,J)*POTH(1    ,NHT  ,TT)+ &
                      C2H(I,J)*POTH(1    ,NHT+1,TT)+ &
                      D2H(I,J)*POTH(2    ,NHT  ,TT)+ &
                      E2H(I,J)*POTH(3    ,NHT  ,TT)
         END IF !(I == 1)
         IF ((I > 1).AND.(I <= NHU)) THEN
            VT2H(I,J)=A2H(I,J)*POTH(I    ,NHT-1,TT)+ &
                      B2H(I,J)*POTH(I    ,NHT  ,TT)+ &
                      C2H(I,J)*POTH(I    ,NHT+1,TT)+ &
                      D2H(I,J)*POTH(I-1  ,NHT  ,TT)+ &
                      E2H(I,J)*POTH(I+1  ,NHT  ,TT)
         END IF !((I > 1).AND.(I <= NHU))
         IF ((I > NHU).AND.(I < NH2)) THEN
            VT2H(I,J)=A2H(I,J)*POTH(I    ,NHT-2,TT)+ &
                      B2H(I,J)*POTH(I    ,NHT-1,TT)+ &
                      C2H(I,J)*POTH(I    ,NHT  ,TT)+ &
                      D2H(I,J)*POTH(I-1  ,NHT  ,TT)+ &
                      E2H(I,J)*POTH(I+1  ,NHT  ,TT)
         END IF !((I > NHU).AND.(I < NH2))
         IF ((I >= NH2).AND.(I < NHX).AND.(NHU /= 0)) THEN
            VT2H(I,J)=A2H(I,J)*POTH(I    ,NHT-2,TT)+ &
                      B2H(I,J)*POTH(I    ,NHT-1,TT)+ &
                      C2H(I,J)*POTH(I    ,NHT  ,TT)+ &
                      D2H(I,J)*POTH(I-1  ,NHT  ,TT)+ &
                      E2H(I,J)*POTH(I+1  ,NHT  ,TT)
         END IF !((I >= NH2).AND.(I < NHX).AND.(NHU /= 0))
         IF ((I == NHX).AND.(NHU == 0)) THEN
            VT2H(I,J)=A2H(I,J)*POTH(NHX  ,NHT-1,TT)+ &
                      B2H(I,J)*POTH(NHX  ,NHT  ,TT)+ &
                      C2H(I,J)*POTH(NHX  ,NHT+1,TT)+ &
                      D2H(I,J)*POTH(NHX-1,NHT  ,TT)+ &
                      E2H(I,J)*POTH(NHX-2,NHT  ,TT)
         END IF !((I == NHX).AND.(NHU == 0))
         IF ((I == NHX).AND.(NHU /= 0)) THEN
            VT2H(I,J)=A2H(I,J)*POTH(NHX  ,NHT-2,TT)+ &
                      B2H(I,J)*POTH(NHX  ,NHT-1,TT)+ &
                      C2H(I,J)*POTH(NHX  ,NHT  ,TT)+ &
                      D2H(I,J)*POTH(NHX-1,NHT  ,TT)+ &
                      E2H(I,J)*POTH(NHX-2,NHT  ,TT)
         END IF !((I == NHX).AND.(NHU /= 0))
      END IF !(J == NHT)
!-----------------------------------------------------------------------------------------------!
      IF (J == NHT1) THEN
         IF (I == 1) THEN
            VT2H(I,J)=A2H(I,J)*POTH(1    ,NHT1-1,TT)+ &
                      B2H(I,J)*POTH(1    ,NHT1  ,TT)+ &
                      C2H(I,J)*POTH(1    ,NHT1+1,TT)+ &
                      D2H(I,J)*POTH(2    ,NHT   ,TT)+ &
                      E2H(I,J)*POTH(3    ,NHT   ,TT)
         END IF !(I == 1)
         IF ((I > 1).AND.(I <= NHU)) THEN
            VT2H(I,J)=A2H(I,J)*POTH(I    ,NHT1-1,TT)+ &
                      B2H(I,J)*POTH(I    ,NHT1  ,TT)+ &
                      C2H(I,J)*POTH(I    ,NHT1+1,TT)+ &
                      D2H(I,J)*POTH(I-1  ,NHT1  ,TT)+ &
                      E2H(I,J)*POTH(I+1  ,NHT1  ,TT)
         END IF !((I > 1).AND.(I <= NHU))
         IF ((I > NHU).AND.(I < NH2)) THEN
            VT2H(I,J)=A2H(I,J)*POTH(I    ,NHT1  ,TT)+ &
                      B2H(I,J)*POTH(I    ,NHT1+1,TT)+ &
                      C2H(I,J)*POTH(I    ,NHT1+2,TT)+ &
                      D2H(I,J)*POTH(I-1  ,NHT1  ,TT)+ &
                      E2H(I,J)*POTH(I+1  ,NHT1  ,TT)
         END IF !((I > NHU).AND.(I < NH2))
         IF ((I >= NH2).AND.(I < NHX).AND.(NHU /= 0)) THEN
            VT2H(I,J)=A2H(I,J)*POTH(I    ,NHT1  ,TT)+ &
                      B2H(I,J)*POTH(I    ,NHT1+1,TT)+ &
                      C2H(I,J)*POTH(I    ,NHT1+2,TT)+ &
                      D2H(I,J)*POTH(I-1  ,NHT1  ,TT)+ &
                      E2H(I,J)*POTH(I+1  ,NHT1  ,TT)
         END IF !((I >= NH2).AND.(I < NHX).AND.(NHU /= 0))
         IF ((I == NHX).AND.(NHU == 0)) THEN
            VT2H(I,J)=A2H(I,J)*POTH(NHX  ,NHT1-1,TT)+ &
                      B2H(I,J)*POTH(NHX  ,NHT1  ,TT)+ &
                      C2H(I,J)*POTH(NHX  ,NHT1+1,TT)+ &
                      D2H(I,J)*POTH(NHX-1,NHT1  ,TT)+ &
                      E2H(I,J)*POTH(NHX-2,NHT1  ,TT)
         END IF !((I == NHX).AND.(NHU == 0))
         IF ((I == NHX).AND.(NHU /= 0)) THEN
            VT2H(I,J)=A2H(I,J)*POTH(NHX  ,NHT1  ,TT)+ &
                      B2H(I,J)*POTH(NHX  ,NHT1+1,TT)+ &
                      C2H(I,J)*POTH(NHX  ,NHT1+2,TT)+ &
                      D2H(I,J)*POTH(NHX-1,NHT1  ,TT)+ &
                      E2H(I,J)*POTH(NHX-2,NHT1  ,TT)
         END IF !((I == NHX).AND.(NHU /= 0))
      END IF !(J == NHT1)
!-----------------------------------------------------------------------------------------------!
      IF ((J /= 1).AND.(J /= NHTP).AND.(J /= NHT).AND.(J /= NHT1)) THEN
         IF (I == 1) THEN
            VT2H(I,J)=A2H(I,J)*POTH(1    ,J-1,TT)+ &
                      B2H(I,J)*POTH(1    ,J  ,TT)+ &
                      C2H(I,J)*POTH(1    ,J+1,TT)+ &
                      D2H(I,J)*POTH(2    ,J  ,TT)+ &
                      E2H(I,J)*POTH(3    ,J  ,TT)
         END IF !(I == 1)
         IF (I == NHX) THEN
            VT2H(I,J)=A2H(I,J)*POTH(NHX  ,J-1,TT)+ &
                      B2H(I,J)*POTH(NHX  ,J  ,TT)+ &
                      C2H(I,J)*POTH(NHX  ,J+1,TT)+ &
                      D2H(I,J)*POTH(NHX-1,J  ,TT)+ &
                      E2H(I,J)*POTH(NHX-2,J  ,TT)
         END IF !(I == NHX)
         IF ((I /= 1).AND.(I /= NHX)) THEN
            VT2H(I,J)=A2H(I,J)*POTH(I    ,J-1,TT)+ &
                      B2H(I,J)*POTH(I    ,J  ,TT)+ &
                      C2H(I,J)*POTH(I    ,J+1,TT)+ &
                      D2H(I,J)*POTH(I-1  ,J  ,TT)+ &
                      E2H(I,J)*POTH(I+1  ,J  ,TT)
         END IF !((I /= 1).AND.(I /= NHX))
      END IF !((J /= 1).AND.(J /= NHTP).AND.(J /= NHT).AND.(J /= NHT1))
!-----------------------------------------------------------------------------------------------!
   END DO !I=1,NHX
END DO !J=1,NHTP
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE VELH
!-----------------------------------------------------------------------------------------------!
