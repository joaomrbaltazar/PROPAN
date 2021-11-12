!-----------------------------------------------------------------------------------------------!
!    Compute Matrices of Influence Coefficients                                                 !
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
SUBROUTINE INFCOEFMATRIXSTD
!-----------------------------------------------------------------------------------------------!
!    Created by: 08112013, J. Baltazar, version 1.0                                             !
!    Modified  : 08112013, J. Baltazar, version 1.0                                             !
!    Modified  : 29052014, J. Baltazar, revised                                                 !
!    Modified  : 27112014, J. Baltazar, version 3.3, Super-Cavitation Model                     !
!    Modified  : 06032015, J. Baltazar, revised                                                 !
!    Modified  : 03062016, J. Baltazar, 2016 version 1.3                                        !
!    Modified  : 24102016, J. Baltazar, 2016 version 1.4                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,K,L,KB,INFO,JP,JS
!-----------------------------------------------------------------------------------------------!
SIJ= SIJ/(2.D0*PI)
DIJ=-DIJ/(2.D0*PI)
WIJ= WIJ/(2.D0*PI)
KIJ= KIJ/(2.D0*PI)
!-----------------------------------------------------------------------------------------------!
DO I=1,NHPAN+NPPAN+NNPAN
   DIJ(I,I,1)=1.D0+DIJ(I,I,1)
END DO !I=1,NHPAN+NPPAN+NNPAN
!-----------------------------------------------------------------------------------------------!
DO J=1,NRW
   DO I=1,IABS(NCPW)
      L=NHPAN+NPPAN+NNPAN+(J-1)*IABS(NCPW)+I
      DIJ(L,L,1)=2.D0
      IF (NT == 0) THEN
         IF (NCPW < 0) WIJ(L,J,1)=WIJ(L,J,1)-1.D0 !Lower side
         IF (NCPW > 0) WIJ(L,J,1)=WIJ(L,J,1)+1.D0 !Upper side
      ELSE !(NT == 0)
         K=(J-1)*NPW+I
         IF (NCPW < 0) WIJ(L,K,1)=WIJ(L,K,1)-1.D0 !Lower side
         IF (NCPW > 0) WIJ(L,K,1)=WIJ(L,K,1)+1.D0 !Upper side
      END IF !(NT == 0)
   END DO !I=1,IABS(NCPW)
END DO !J=1,NRW
!-----------------------------------------------------------------------------------------------!
!    Store Matrix of Influence Coefficients on File                                             !
!-----------------------------------------------------------------------------------------------!
REWIND(21)
WRITE (21) ((DIJ(I,J,1),I=1,NPAN),J=1,NPAN )
!-----------------------------------------------------------------------------------------------!
!    Define Steady Matrices                                                                     !
!-----------------------------------------------------------------------------------------------!
DO KB=2,NB
   DIJ(:,:,1)=DIJ(:,:,1)+DIJ(:,:,KB)
END DO !KB=2,NB
!-----------------------------------------------------------------------------------------------!
!    LU Decomposition                                                                           !
!-----------------------------------------------------------------------------------------------!
IF (ISOLVER == 0) THEN
   INFO=0
   CALL DGEFA(DIJ(:,:,1),NPAN1,NPAN,IPVT1,INFO)
   IF (INFO /= 0) THEN
      WRITE(6,*) INFO
      STOP
   END IF !(INFO /= 0)
!-----------------------------------------------------------------------------------------------!
!    Store LU Factored Matrix                                                                   !
!-----------------------------------------------------------------------------------------------!
   REWIND(22)
   WRITE (22) ((DIJ(I,J,1),I=1,NPAN),J=1,NPAN)
END IF !(ISOLVER == 0)
!-----------------------------------------------------------------------------------------------!
!    Linear Kutta Condition                                                                     !
!-----------------------------------------------------------------------------------------------!
DIJ(:,:,1)=0.d0
REWIND(21)
READ  (21) ((DIJ(I,J,1),I=1,NPAN),J=1,NPAN)
DO KB=2,NB
   DIJ(:,:,1)=DIJ(:,:,1)+DIJ(:,:,KB)
END DO !KB=2,NB
!-----------------------------------------------------------------------------------------------!
!    Kutta Condition at Blade                                                                   !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRW
   DO I=1,NPAN
      DIJ(I,NPAN+J,1)=0.D0
      DO KB=1,NB
         IF (NT == 0) THEN
            DIJ(I,NPAN+J,1)=DIJ(I,NPAN+J,1)-WIJ(I,J,KB)
         ELSE !(NT == 0)
            DO K=((J-1)*NPW+1),(J*NPW)
               DIJ(I,NPAN+J,1)=DIJ(I,NPAN+J,1)-WIJ(I,K,KB)
            END DO !K=((J-1)*NPW+1),(J*NPW)
         END IF !(NT == 0)
      END DO !KB=1,NB
   END DO !I=1,NPAN
   JP=NHPAN+(JI-1)*NCP+(J-1)*NCP+1
   JS=NHPAN+(JI-1)*NCP+J*NCP
   DIJ(NPAN+J,    JP,1)= 1.D0
   DIJ(NPAN+J,    JS,1)=-1.D0
   DIJ(NPAN+J,NPAN+J,1)= 1.D0
END DO !J=1,NRW
!-----------------------------------------------------------------------------------------------!
!    Kutta Condition at Nozzle                                                                  !
!-----------------------------------------------------------------------------------------------!
DO J=1,NNTP
   DO I=1,NPAN
      DIJ(I,NPAN+NRW+J,1)=0.D0
      DO KB=1,NB
         IF (NT == 0) THEN
            DIJ(I,NPAN+NRW+J,1)=DIJ(I,NPAN+NRW+J,1)-WIJ(I,NRW+J,KB)
         ELSE !(NT == 0)
            DO K=((J-1)*NNW+1),(J*NNW)
               DIJ(I,NPAN+NRW+J,1)=DIJ(I,NPAN+NRW+J,1)-WIJ(I,NPWPAN+K,KB)
            END DO !K=((J-1)*NNW+1),(J*NNW)
         END IF !(NT == 0)
      END DO !KB=1,NB
   END DO !I=1,NPAN
   JP=NHPAN+NPPAN+(J-1)*NNXT1+NNX1
   JS=NHPAN+NPPAN+(J-1)*NNXT1+NNX
   DIJ(NPAN+NRW+J,        JP,1)= 1.D0
   DIJ(NPAN+NRW+J,        JS,1)=-1.D0
   DIJ(NPAN+NRW+J,NPAN+NRW+J,1)= 1.D0
END DO !J=1,NNTP
!-----------------------------------------------------------------------------------------------!
!    Store Linear Kutta Matrix                                                                  !
!-----------------------------------------------------------------------------------------------!
REWIND(23)
WRITE (23) ((DIJ(I,J,1),I=1,NPAN1),J=1,NPAN1)
!-----------------------------------------------------------------------------------------------!
!    LU Decomposition                                                                           !
!-----------------------------------------------------------------------------------------------!
IF (ISOLVER == 0) THEN
   INFO=0
   CALL DGEFA(DIJ(:,:,1),NPAN1,NPAN1,IPVT2,INFO)
   IF (INFO /= 0) THEN
      WRITE(6,*) INFO
      STOP
   END IF !(INFO /= 0)
!-----------------------------------------------------------------------------------------------!
!    Store LU Factored Matrix                                                                   !
!-----------------------------------------------------------------------------------------------!
   REWIND(24)
   WRITE (24) ((DIJ(I,J,1),I=1,NPAN1),J=1,NPAN1)
END IF !(ISOLVER == 0)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE INFCOEFMATRIXSTD
!-----------------------------------------------------------------------------------------------!
