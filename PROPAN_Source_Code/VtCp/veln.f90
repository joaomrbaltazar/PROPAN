!-----------------------------------------------------------------------------------------------!
!    Calculation of the tangential perturbation velocities on the nozzle                        !
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
SUBROUTINE VELN(TT)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 15112013, J. Baltazar, version 1.0                                             !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,TT
!-----------------------------------------------------------------------------------------------!
DO J=1,NNTP
   DO I=1,NNXT1
!-----------------------------------------------------------------------------------------------!
!    u1 Derivative                                                                              !
!-----------------------------------------------------------------------------------------------!
      IF (I == 1) THEN
         VT1N(I,J)=A1N(I,J)*POTN(NNXT1  ,J,TT)+ &
                   B1N(I,J)*POTN(1      ,J,TT)+ &
                   C1N(I,J)*POTN(2      ,J,TT)
      ELSEIF ((I == NNX).OR.(I == NNU).OR.(I == (NN2-1))) THEN
         VT1N(I,J)=A1N(I,J)*POTN(I-2    ,J,TT)+ &
                   B1N(I,J)*POTN(I-1    ,J,TT)+ &
                   C1N(I,J)*POTN(I      ,J,TT)
      ELSEIF ((I == NNX1).OR.(I == (NNU+2)).OR.(I == (NN2+1))) THEN
         VT1N(I,J)=A1N(I,J)*POTN(I      ,J,TT)+ &
                   B1N(I,J)*POTN(I+1    ,J,TT)+ &
                   C1N(I,J)*POTN(I+2    ,J,TT)
      ELSEIF (I == NNXT1) THEN
         VT1N(I,J)=A1N(I,J)*POTN(NNXT1-1,J,TT)+ &
                   B1N(I,J)*POTN(NNXT1  ,J,TT)+ &
                   C1N(I,J)*POTN(1      ,J,TT)
      ELSE
         VT1N(I,J)=A1N(I,J)*POTN(I-1    ,J,TT)+ &
                   B1N(I,J)*POTN(I      ,J,TT)+ &
                   C1N(I,J)*POTN(I+1    ,J,TT)
      END IF !(I)
!-----------------------------------------------------------------------------------------------!
!    u2 Derivative                                                                              !
!-----------------------------------------------------------------------------------------------!
      IF (J == 1) THEN
         IF (I == 1) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNTP,TT)+ &
                      B2N(I,J)*POTN(I    ,1   ,TT)+ &
                      C2N(I,J)*POTN(I    ,2   ,TT)+ &
                      D2N(I,J)*POTN(NNXT1,1   ,TT)+ &
                      E2N(I,J)*POTN(I+1  ,1   ,TT)
         ELSEIF ((I == NNX).OR.(I == NNU).OR.(I == (NN2-1))) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNTP,TT)+ &
                      B2N(I,J)*POTN(I    ,1   ,TT)+ &
                      C2N(I,J)*POTN(I    ,2   ,TT)+ &
                      D2N(I,J)*POTN(I-1  ,1   ,TT)+ &
                      E2N(I,J)*POTN(I-2  ,1   ,TT)
         ELSEIF ((I == NNX1).OR.(I == (NNU+2)).OR.(I == (NN2+1))) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNTP,TT)+ &
                      B2N(I,J)*POTN(I    ,1   ,TT)+ &
                      C2N(I,J)*POTN(I    ,2   ,TT)+ &
                      D2N(I,J)*POTN(I+1  ,1   ,TT)+ &
                      E2N(I,J)*POTN(I+2  ,1   ,TT)
         ELSEIF (I == NNXT1) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNTP,TT)+ &
                      B2N(I,J)*POTN(I    ,1   ,TT)+ &
                      C2N(I,J)*POTN(I    ,2   ,TT)+ &
                      D2N(I,J)*POTN(I-1  ,1   ,TT)+ &
                      E2N(I,J)*POTN(1    ,1   ,TT)
         ELSE !(I)
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNTP,TT)+ &
                      B2N(I,J)*POTN(I    ,1   ,TT)+ &
                      C2N(I,J)*POTN(I    ,2   ,TT)+ &
                      D2N(I,J)*POTN(I-1  ,1   ,TT)+ &
                      E2N(I,J)*POTN(I+1  ,1   ,TT)
         END IF !(I)
!-----------------------------------------------------------------------------------------------!
      ELSEIF (J == NNT) THEN
         IF (I == 1) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNT-1,TT)+ &
                      B2N(I,J)*POTN(I    ,NNT  ,TT)+ &
                      C2N(I,J)*POTN(I    ,NNT+1,TT)+ &
                      D2N(I,J)*POTN(NNXT1,NNT  ,TT)+ &
                      E2N(I,J)*POTN(I+1  ,NNT  ,TT)
         ELSEIF ((I == NNX).OR.(I == NNU).OR.(I == (NN2-1))) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNT-1,TT)+ &
                      B2N(I,J)*POTN(I    ,NNT  ,TT)+ &
                      C2N(I,J)*POTN(I    ,NNT+1,TT)+ &
                      D2N(I,J)*POTN(I-1  ,NNT  ,TT)+ &
                      E2N(I,J)*POTN(I-2  ,NNT  ,TT)
         ELSEIF ((I == NNX1).OR.(I == (NNU+2)).OR.(I == (NN2+1))) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNT-1,TT)+ &
                      B2N(I,J)*POTN(I    ,NNT  ,TT)+ &
                      C2N(I,J)*POTN(I    ,NNT+1,TT)+ &
                      D2N(I,J)*POTN(I+1  ,NNT  ,TT)+ &
                      E2N(I,J)*POTN(I+2  ,NNT  ,TT)
         ELSEIF (I == NNXT1) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNT-1,TT)+ &
                      B2N(I,J)*POTN(I    ,NNT  ,TT)+ &
                      C2N(I,J)*POTN(I    ,NNT+1,TT)+ &
                      D2N(I,J)*POTN(I-1  ,NNT  ,TT)+ &
                      E2N(I,J)*POTN(1    ,NNT  ,TT)
         ELSE !(I)
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNT-1,TT)+ &
                      B2N(I,J)*POTN(I    ,NNT  ,TT)+ &
                      C2N(I,J)*POTN(I    ,NNT+1,TT)+ &
                      D2N(I,J)*POTN(I-1  ,NNT  ,TT)+ &
                      E2N(I,J)*POTN(I+1  ,NNT  ,TT)
         END IF !(I)
         IF (ISTRIP == 1) THEN
            IF (((I > NNU).AND.(I < NNX)).or.((i > nnx1).and.(i <= nntp-nnu+1))) THEN
               VT2N(I,J)=A2N(I,J)*POTN(I  ,NNT-2,TT)+ &
                         B2N(I,J)*POTN(I  ,NNT-1,TT)+ &
                         C2N(I,J)*POTN(I  ,NNT  ,TT)+ &
                         D2N(I,J)*POTN(I-1,NNT  ,TT)+ &
                         E2N(I,J)*POTN(I+1,NNT  ,TT)
            END IF !((I > NNU).AND.(I < NNX))
            IF ((I == NNX).OR.(I == NNU).OR.(I == (NN2-1))) THEN
               VT2N(I,J)=A2N(I,J)*POTN(I  ,NNT-2,TT)+ &
                         B2N(I,J)*POTN(I  ,NNT-1,TT)+ &
                         C2N(I,J)*POTN(I  ,NNT  ,TT)+ &
                         D2N(I,J)*POTN(I-1,NNT  ,TT)+ &
                         E2N(I,J)*POTN(I-2,NNT  ,TT)
            END IF !((I == NNX).OR.(I == NNU).OR.(I == (NN2-1)))
            IF ((I == (NNU+2)).OR.(I == (NN2+1)).or.(i == nnx1)) THEN
               VT2N(I,J)=A2N(I,J)*POTN(I    ,NNT-2,TT)+ &
                         B2N(I,J)*POTN(I    ,NNT-1,TT)+ &
                         C2N(I,J)*POTN(I    ,NNT  ,TT)+ &
                         D2N(I,J)*POTN(I+1  ,NNT  ,TT)+ &
                         E2N(I,J)*POTN(I+2  ,NNT  ,TT)
            END IF !((I == (NNU+2)).OR.(I == (NN2+1)))
         END IF !(ISTRIP == 1)
!-----------------------------------------------------------------------------------------------!
      ELSEIF (J == NNT+1) THEN
         IF (I == 1) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNT  ,TT)+ &
                      B2N(I,J)*POTN(I    ,NNT+1,TT)+ &
                      C2N(I,J)*POTN(I    ,NNT+2,TT)+ &
                      D2N(I,J)*POTN(NNXT1,NNT+1,TT)+ &
                      E2N(I,J)*POTN(I+1  ,NNT+1,TT)
         ELSEIF ((I == NNX).OR.(I == NNU).OR.(I == (NN2-1))) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNT  ,TT)+ &
                      B2N(I,J)*POTN(I    ,NNT+1,TT)+ &
                      C2N(I,J)*POTN(I    ,NNT+2,TT)+ &
                      D2N(I,J)*POTN(I-1  ,NNT+1,TT)+ &
                      E2N(I,J)*POTN(I-2  ,NNT+1,TT)
         ELSEIF ((I == NNX1).OR.(I == (NNU+2)).OR.(I == (NN2+1))) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNT  ,TT)+ &
                      B2N(I,J)*POTN(I    ,NNT+1,TT)+ &
                      C2N(I,J)*POTN(I    ,NNT+2,TT)+ &
                      D2N(I,J)*POTN(I+1  ,NNT+1,TT)+ &
                      E2N(I,J)*POTN(I+2  ,NNT+1,TT)
         ELSEIF (I == NNXT1) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNT  ,TT)+ &
                      B2N(I,J)*POTN(I    ,NNT+1,TT)+ &
                      C2N(I,J)*POTN(I    ,NNT+2,TT)+ &
                      D2N(I,J)*POTN(I-1  ,NNT+1,TT)+ &
                      E2N(I,J)*POTN(1    ,NNT+1,TT)
         ELSE !(I)
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNT  ,TT)+ &
                      B2N(I,J)*POTN(I    ,NNT+1,TT)+ &
                      C2N(I,J)*POTN(I    ,NNT+2,TT)+ &
                      D2N(I,J)*POTN(I-1  ,NNT+1,TT)+ &
                      E2N(I,J)*POTN(I+1  ,NNT+1,TT)
         END IF !(I)
         IF (ISTRIP == 1) THEN
            IF (((I > NNU).AND.(I < NNX)).or.((i > nnx1).and.(i <= nntp-nnu+1))) THEN
               VT2N(I,J)=A2N(I,J)*POTN(I  ,NNT+1,TT)+ &
                         B2N(I,J)*POTN(I  ,NNT+2,TT)+ &
                         C2N(I,J)*POTN(I  ,NNT+3,TT)+ &
                         D2N(I,J)*POTN(I-1,NNT+1,TT)+ &
                         E2N(I,J)*POTN(I+1,NNT+1,TT)
            END IF !((I > NNU).AND.(I < NNX))
            IF ((I == NNX).OR.(I == NNU).OR.(I == (NN2-1))) THEN
               VT2N(I,J)=A2N(I,J)*POTN(I  ,NNT+1,TT)+ &
                         B2N(I,J)*POTN(I  ,NNT+2,TT)+ &
                         C2N(I,J)*POTN(I  ,NNT+3,TT)+ &
                         D2N(I,J)*POTN(I-1,NNT+1,TT)+ &
                         E2N(I,J)*POTN(I-2,NNT+1,TT)
            END IF !((I == NNX).OR.(I == NNU).OR.(I == (NN2-1)))
            IF ((I == (NNU+2)).OR.(I == (NN2+1)).or.(i == nnx1)) THEN
               VT2N(I,J)=A2N(I,J)*POTN(I    ,NNT+1,TT)+ &
                         B2N(I,J)*POTN(I    ,NNT+2,TT)+ &
                         C2N(I,J)*POTN(I    ,NNT+3,TT)+ &
                         D2N(I,J)*POTN(I+1  ,NNT+1,TT)+ &
                         E2N(I,J)*POTN(I+2  ,NNT+1,TT)
            END IF !((I == (NNU+2)).OR.(I == (NN2+1)))
         END IF !(ISTRIP == 1)
!-----------------------------------------------------------------------------------------------!
      ELSEIF (J == NNTP) THEN
         IF (I == 1) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNTP-1,TT)+ &
                      B2N(I,J)*POTN(I    ,NNTP  ,TT)+ &
                      C2N(I,J)*POTN(I    ,1     ,TT)+ &
                      D2N(I,J)*POTN(NNXT1,NNTP  ,TT)+ &
                      E2N(I,J)*POTN(I+1  ,NNTP  ,TT)
         ELSEIF ((I == NNX).OR.(I == NNU).OR.(I == (NN2-1))) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNTP-1,TT)+ &
                      B2N(I,J)*POTN(I    ,NNTP  ,TT)+ &
                      C2N(I,J)*POTN(I    ,1     ,TT)+ &
                      D2N(I,J)*POTN(I-1  ,NNTP  ,TT)+ &
                      E2N(I,J)*POTN(I-2  ,NNTP  ,TT)
         ELSEIF ((I == NNX1).OR.(I == (NNU+2)).OR.(I == (NN2+1))) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNTP-1,TT)+ &
                      B2N(I,J)*POTN(I    ,NNTP  ,TT)+ &
                      C2N(I,J)*POTN(I    ,1     ,TT)+ &
                      D2N(I,J)*POTN(I+1  ,NNTP  ,TT)+ &
                      E2N(I,J)*POTN(I+2  ,NNTP  ,TT)
         ELSEIF (I == NNXT1) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNTP-1,TT)+ &
                      B2N(I,J)*POTN(I    ,NNTP  ,TT)+ &
                      C2N(I,J)*POTN(I    ,1     ,TT)+ &
                      D2N(I,J)*POTN(I-1  ,NNTP  ,TT)+ &
                      E2N(I,J)*POTN(1    ,NNTP  ,TT)
         ELSE !(I)
            VT2N(I,J)=A2N(I,J)*POTN(I    ,NNTP-1,TT)+ &
                      B2N(I,J)*POTN(I    ,NNTP  ,TT)+ &
                      C2N(I,J)*POTN(I    ,1     ,TT)+ &
                      D2N(I,J)*POTN(I-1  ,NNTP  ,TT)+ &
                      E2N(I,J)*POTN(I+1  ,NNTP  ,TT)
         END IF !(I)
!-----------------------------------------------------------------------------------------------!
      ELSE !(J)
         IF (I == 1) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,J-1,TT)+ &
                      B2N(I,J)*POTN(I    ,J  ,TT)+ &
                      C2N(I,J)*POTN(I    ,J+1,TT)+ &
                      D2N(I,J)*POTN(NNXT1,J  ,TT)+ &
                      E2N(I,J)*POTN(I+1  ,J  ,TT)
         ELSEIF ((I == NNX).OR.(I == NNU).OR.(I == (NN2-1))) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,J-1,TT)+ &
                      B2N(I,J)*POTN(I    ,J  ,TT)+ &
                      C2N(I,J)*POTN(I    ,J+1,TT)+ &
                      D2N(I,J)*POTN(I-1  ,J  ,TT)+ &
                      E2N(I,J)*POTN(I-2  ,J  ,TT)
         ELSEIF ((I == NNX1).OR.(I == (NNU+2)).OR.(I == (NN2+1))) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,J-1,TT)+ &
                      B2N(I,J)*POTN(I    ,J  ,TT)+ &
                      C2N(I,J)*POTN(I    ,J+1,TT)+ &
                      D2N(I,J)*POTN(I+1  ,J  ,TT)+ &
                      E2N(I,J)*POTN(I+2  ,J  ,TT)
         ELSEIF (I == NNXT1) THEN
            VT2N(I,J)=A2N(I,J)*POTN(I    ,J-1,TT)+ &
                      B2N(I,J)*POTN(I    ,J  ,TT)+ &
                      C2N(I,J)*POTN(I    ,J+1,TT)+ &
                      D2N(I,J)*POTN(I-1  ,J  ,TT)+ &
                      E2N(I,J)*POTN(1    ,J  ,TT)
         ELSE !(I)
            VT2N(I,J)=A2N(I,J)*POTN(I    ,J-1,TT)+ &
                      B2N(I,J)*POTN(I    ,J  ,TT)+ &
                      C2N(I,J)*POTN(I    ,J+1,TT)+ &
                      D2N(I,J)*POTN(I-1  ,J  ,TT)+ &
                      E2N(I,J)*POTN(I+1  ,J  ,TT)
         END IF !(I)
      END IF !(J)
!-----------------------------------------------------------------------------------------------!
   END DO !I=1,NNXT1
END DO !J=1,NNTP
!-----------------------------------------------------------------------------------------------!
!    u1 Derivative Correction                                                                   !
!-----------------------------------------------------------------------------------------------!
VT1N(NNU+1,:)=0.5D0*(VT1N(NNU  ,:)+VT1N(NNU+2,:))
VT1N(NN2  ,:)=0.5D0*(VT1N(NN2-1,:)+VT1N(NN2+1,:))
!-----------------------------------------------------------------------------------------------!
!    u2 Derivative Correction                                                                   !
!-----------------------------------------------------------------------------------------------!
VT2N(NNU+1,:)=0.5D0*(VT2N(NNU  ,:)+VT2N(NNU+2,:))
VT2N(NN2  ,:)=0.5D0*(VT2N(NN2-1,:)+VT2N(NN2+1,:))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE VELN
!-----------------------------------------------------------------------------------------------!
