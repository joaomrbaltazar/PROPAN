!-----------------------------------------------------------------------------------------------!
!    Compute Matrices of Influence Coefficients. Unsteady Case                                  !
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
SUBROUTINE INFCOEFMATRIXUNSTD
!-----------------------------------------------------------------------------------------------!
!    Created by: 08112013, J. Baltazar, version 1.0                                             !
!    Modified  : 25112013, J. Baltazar, version 1.0                                             !
!    Modified  : 30052014, J. Baltazar, revised                                                 !
!    Modified  : 03062016, J. Baltazar, 2016 version 1.3                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER I,J,INFO,JP,JS
!-----------------------------------------------------------------------------------------------!
!    Recover Matrix of Influence Coefficients from File                                         !
!-----------------------------------------------------------------------------------------------!
REWIND(21)
READ  (21) ((DIJ(I,J,1),I=1,NPAN),J=1,NPAN)
!-----------------------------------------------------------------------------------------------!
!    LU Decomposition                                                                           !
!-----------------------------------------------------------------------------------------------!
IF (ISOLVER == 0) THEN
   INFO=0
   CALL DGEFA(DIJ(:,:,1),NPAN1,NPAN,IPVT3,INFO)
   IF (INFO /= 0) THEN
      WRITE(6,*) INFO
      STOP
   END IF !(INFO /= 0)
!-----------------------------------------------------------------------------------------------!
!    Store LU Factored Matrix                                                                   !
!-----------------------------------------------------------------------------------------------!
   REWIND(25)
   WRITE (25) ((DIJ(I,J,1),I=1,NPAN),J=1,NPAN)
END IF !(ISOLVER == 0)
!-----------------------------------------------------------------------------------------------!
!    Linear Kutta Condition                                                                     !
!-----------------------------------------------------------------------------------------------!
DIJ(:,:,1)=0.D0
REWIND(21)
READ  (21) ((DIJ(I,J,1),I=1,NPAN),J=1,NPAN)
!-----------------------------------------------------------------------------------------------!
!    Kutta Condition at Blade                                                                   !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRW
   DO I=1,NPAN
      DIJ(I,NPAN+J,1)=-KIJ(I,J,1)
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
      DIJ(I,NPAN+NRW+J,1)=-KIJ(I,NRW+J,1)
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
REWIND(26)
WRITE (26) ((DIJ(I,J,1),I=1,NPAN1),J=1,NPAN1)
!-----------------------------------------------------------------------------------------------!
!    LU Decomposition                                                                           !
!-----------------------------------------------------------------------------------------------!
IF (ISOLVER == 0) THEN
   INFO=0
   CALL DGEFA(DIJ(:,:,1),NPAN1,NPAN1,IPVT4,INFO)
   IF (INFO /= 0) THEN
      WRITE(6,*) INFO
      STOP
   END IF !(INFO /= 0)
!-----------------------------------------------------------------------------------------------!
!    Store LU Factored Matrix                                                                   !
!-----------------------------------------------------------------------------------------------!
   REWIND(27)
   WRITE (27) ((DIJ(I,J,1),I=1,NPAN1),J=1,NPAN1)
END IF !(ISOLVER == 0)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE INFCOEFMATRIXUNSTD
!-----------------------------------------------------------------------------------------------!
