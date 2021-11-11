!-----------------------------------------------------------------------------------------------!
!    Compute Numerically the Jacobian Matrix. Disturb the Basic Solution                        !
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
SUBROUTINE JACOBIANCAV(JJ,TT,BETA1,PJACOB)
!-----------------------------------------------------------------------------------------------!
!    Created by: 06112014, J. Baltazar, Cavitation Model                                        !
!    Modified  : 24112014, J. Baltazar, version 3.3, Super-Cavitation Model                     !
!    Modified  : 29112014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!    Modified  : 27052016, J. Baltazar, version 1.2                                             !
!    Modified  : 02062016, J. Baltazar, 2016 version 1.3                                        !
!    Modified  : 25102016, J. Baltazar, 2016 version 1.4                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,J2,L,K,KB,JJ,TT
DOUBLE PRECISION :: BETA1
DOUBLE PRECISION :: DCPBETA(NRW),DCNBETA(NNTP)
DOUBLE PRECISION :: DELPOT(NPAN),PJACOB(NRT,NRT)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
DO L=1,NRW
!-----------------------------------------------------------------------------------------------!
   DELPOT=0.D0
   IF (NT == 0) THEN
      DO KB=1,NB
         DO I=1,NPAN
            DELPOT(I)=DELPOT(I)+WIJ(I,L,KB)*BETA1*DPOTP(L)
         END DO !I=1,NPAN
      END DO !KB=1,NB
   ELSEIF (TT == 0) THEN
      DO KB=1,NB
         DO I=1,NPAN
            DO K=((L-1)*NPW+1),(L*NPW)
               DELPOT(I)=DELPOT(I)+WIJ(I,K,KB)*BETA1*DPOTP(L)
            END DO !K=((L-1)*NPW+1),(L*NPW)
         END DO !I=1,NPAN
      END DO !KB=1,NB
   ELSE !(TT)
      DO I=1,NPAN
         DELPOT(I)=KIJ(I,L,1)*BETA1*DPOTP(L)
      END DO !I=1,NPAN
   END IF !(TT)
!-----------------------------------------------------------------------------------------------!
!    Solution of System of Equations                                                            !
!-----------------------------------------------------------------------------------------------!
   IF (ISOLVER == 0) CALL DGESL(DIJ(:,:,1),NPAN1,NPAN,IPVTC1,DELPOT,0)
!-----------------------------------------------------------------------------------------------!
!    Iterative Solver                                                                           !
!-----------------------------------------------------------------------------------------------!
   IF (ISOLVER == 1) CALL BISOF(DIJ(:,:,1),NPAN1,NPAN,DELPOT)
!-----------------------------------------------------------------------------------------------!
!    Update the Solution of the Linear Kutta Condition                                          !
!-----------------------------------------------------------------------------------------------!
   SI(1:NPAN)=RHS(1:NPAN)+DELPOT(1:NPAN)
!-----------------------------------------------------------------------------------------------!
!    Structure the Potential Solution                                                           !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
   DO J=1,NRP
      DO I=1,NCP
         K=NHPAN+(J-1)*NCP+I
         IF ((I <= IDP(J,TT)).AND.(I >= IRP(J,TT))) THEN
            SOURCEP(I,J, 1)=SI(K)
         ELSEIF ((I >= IDS(J,TT)).AND.(I <= IRS(J,TT))) THEN
            SOURCEP(I,J, 1)=SI(K)
         ELSE
            POTP   (I,J,TT)=SI(K)
         END IF
      END DO !I=1,NCP
   END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
   POTN(:,:,TT)=RESHAPE(SI(NHPAN+NPPAN+1:NHPAN+NPPAN+NNPAN),(/NNXT1,NNTP/))
!-----------------------------------------------------------------------------------------------!
!    Blade Velocities and Pressure                                                              !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      CALL VELP(TT)
      CALL PRESP(JJ,TT)
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Velocities and Pressure                                                             !
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      CALL VELN(TT)
      CALL PRESN(JJ,TT)
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Jacobian Matrix                                                                            !
!-----------------------------------------------------------------------------------------------!
   DO J=1,NRW
      J2=JI-1+J
      DCPBETA(J)=CPP(NCP,J2,TT)-CPP(1,J2,TT)
      PJACOB(J,L)=(DCPBETA(J)-DCP(J))/(BETA1*DPOTP(L))
   END DO !I=1,NRW
!-----------------------------------------------------------------------------------------------!
   CALL PRESNTE(0,DCNBETA)
   DO J=1,NNTP
!*    DCNBETA(J)=CPN(NNX,J,TT)-CPN(NNX1,J,TT)
      PJACOB(NRW+J,L)=(DCNBETA(J)-DCN(J))/(BETA1*DPOTP(L))
   END DO !I=1,NNTP
!-----------------------------------------------------------------------------------------------!
END DO !L=1,NRW
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
DO L=1,NNTP
!-----------------------------------------------------------------------------------------------!
   DELPOT=0.D0
   IF (NT == 0) THEN
      DO KB=1,NB
         DO I=1,NPAN
            DELPOT(I)=DELPOT(I)+WIJ(I,NRW+L,KB)*BETA1*DPOTN(L)
         END DO !I=1,NPAN
      END DO !KB=1,NB
   ELSEIF (TT == 0) THEN
      DO KB=1,NB
         DO I=1,NPAN
            DO K=((L-1)*NNW+1),(L*NNW)
               DELPOT(I)=DELPOT(I)+WIJ(I,NPWPAN+K,KB)*BETA1*DPOTN(L)
            END DO !K=((L-1)*NNW+1),(L*NNW)
         END DO !I=1,NPAN
      END DO !KB=1,NB
   ELSE !(TT)
      DO I=1,NPAN
         DELPOT(I)=KIJ(I,NRW+L,1)*BETA1*DPOTN(L)
      END DO !I=1,NPAN
   END IF !(TT)
!-----------------------------------------------------------------------------------------------!
!    Solution of System of Equations                                                            !
!-----------------------------------------------------------------------------------------------!
   IF (ISOLVER == 0) CALL DGESL(DIJ(:,:,1),NPAN1,NPAN,IPVTC1,DELPOT,0)
!-----------------------------------------------------------------------------------------------!
!    Iterative Solver                                                                           !
!-----------------------------------------------------------------------------------------------!
   IF (ISOLVER == 1) CALL BISOF(DIJ(:,:,1),NPAN1,NPAN,DELPOT)
!-----------------------------------------------------------------------------------------------!
!    Update the Solution of the Linear Kutta Condition                                          !
!-----------------------------------------------------------------------------------------------!
   SI(1:NPAN)=RHS(1:NPAN)+DELPOT(1:NPAN)
!-----------------------------------------------------------------------------------------------!
!    Structure the Potential Solution                                                           !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
   DO J=1,NRP
      DO I=1,NCP
         K=NHPAN+(J-1)*NCP+I
         IF ((I <= IDP(J,TT)).AND.(I >= IRP(J,TT))) THEN
            SOURCEP(I,J, 1)=SI(K)
         ELSEIF ((I >= IDS(J,TT)).AND.(I <= IRS(J,TT))) THEN
            SOURCEP(I,J, 1)=SI(K)
         ELSE
            POTP   (I,J,TT)=SI(K)
         END IF
      END DO !I=1,NCP
   END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
   POTN(:,:,TT)=RESHAPE(SI(NHPAN+NPPAN+1:NHPAN+NPPAN+NNPAN),(/NNXT1,NNTP/))
!-----------------------------------------------------------------------------------------------!
!    Blade Velocities and Pressure                                                              !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      CALL VELP(TT)
      CALL PRESP(JJ,TT)
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Velocities and Pressure                                                             !
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      CALL VELN(TT)
      CALL PRESN(JJ,TT)
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Jacobian Matrix                                                                            !
!-----------------------------------------------------------------------------------------------!
   DO J=1,NRW
      J2=JI-1+J
      DCPBETA(J)=CPP(NCP,J2,TT)-CPP(1,J2,TT)
      PJACOB(J,NRW+L)=(DCPBETA(J)-DCP(J))/(BETA1*DPOTN(L))
   END DO !I=1,NRW
!-----------------------------------------------------------------------------------------------!
   CALL PRESNTE(TT,DCNBETA)
   DO J=1,NNTP
!*    DCNBETA(J)=CPN(NNX,J,TT)-CPN(NNX1,J,TT)
      PJACOB(NRW+J,NRW+L)=(DCNBETA(J)-DCN(J))/(BETA1*DPOTN(L))
   END DO !I=1,NNTP
!-----------------------------------------------------------------------------------------------!
END DO !L=1,NNTP
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE JACOBIANCAV
!-----------------------------------------------------------------------------------------------!
