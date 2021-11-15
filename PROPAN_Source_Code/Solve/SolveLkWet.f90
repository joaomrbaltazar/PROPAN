!-----------------------------------------------------------------------------------------------!
!    Solve wetted linear Kutta condition                                                        !
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
SUBROUTINE SOLVELKWET(JJ,TT)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 20112013, J. Baltazar, version 1.0                                             !
!    Modified  : 28102014, J. Baltazar, Interpolation scheme at nozzle t.e.                     !
!    Modified  : 06112014, J. Baltazar, VELH and PRESH only for IH=1                            !
!    Modified  : 27112014, J. Baltazar, version 3.3, Super-Cavitation Model                     !
!    Modified  : 26052016, J. Baltazar, version 1.2, combine Std and UnStd                      !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,J2,JJ,TT,JP,JS
INTEGER :: INFO,IPVTC(NPAN1)
!-----------------------------------------------------------------------------------------------!
!    Recover the Original Linear Kutta Influence Matrix                                         !
!-----------------------------------------------------------------------------------------------!
!    Steady Case                                                                                !
!-----------------------------------------------------------------------------------------------!
IF (TT == 0) THEN
   IF ((ISOLVER == 0).AND.(NWCAV == 0)) THEN
      REWIND(24)
      READ  (24) ((DIJ(I,J,1),I=1,NPAN1),J=1,NPAN1)
   ELSE !((ISOLVER == 0).AND.(NWCAV == 0))
      REWIND(23)
      READ  (23) ((DIJ(I,J,1),I=1,NPAN1),J=1,NPAN1)
   END IF !((ISOLVER == 0).AND.(NWCAV == 0))
END IF !(TT == 0)
!-----------------------------------------------------------------------------------------------!
!    Unsteady Case                                                                              !
!-----------------------------------------------------------------------------------------------!
IF (TT /= 0) THEN
   IF ((ISOLVER == 0).AND.(NWCAV == 0)) THEN
      REWIND(27)
      READ  (27) ((DIJ(I,J,1),I=1,NPAN1),J=1,NPAN1)
   ELSE !((ISOLVER == 0).AND.(NWCAV == 0))
      REWIND(26)
      READ  (26) ((DIJ(I,J,1),I=1,NPAN1),J=1,NPAN1)
   END IF !((ISOLVER == 0).AND.(NWCAV == 0))
END IF !(TT /= 0)
!-----------------------------------------------------------------------------------------------!
!    Correction for Super-Cavitation in the Kutta Condition                                     !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRW
   J2=JI-1+J
!-----------------------------------------------------------------------------------------------!
   IF (IRP(J2,TT) == 1) THEN
      JP=NHPAN+(JI-1)*NCP+(J-1)*NCP+1
      DIJ(NPAN+J,JP,1)=0.D0
   END IF !(IRP(J2,TT) == 1)
!-----------------------------------------------------------------------------------------------!
   IF (IRS(J2,TT) == NCP) THEN
      JS=NHPAN+(JI-1)*NCP+J*NCP
      DIJ(NPAN+J,JS,1)=0.D0
   END IF !(IRS(J2,TT) == NCP)
END DO !J=1,NRW
!-----------------------------------------------------------------------------------------------!
!    LU Decomposition for Super-Cavitation                                                      !
!-----------------------------------------------------------------------------------------------!
IF ((ISOLVER == 0).AND.(NWCAV /= 0)) THEN
   INFO=0
   CALL DGEFA(DIJ(:,:,1),NPAN1,NPAN1,IPVTC,INFO)
   IF (INFO /= 0) THEN
      WRITE(6,*) INFO
      STOP
   END IF !(INFO /= 0)
END IF !((ISOLVER == 0).AND.(NWCAV /= 0))
!-----------------------------------------------------------------------------------------------!
IF ((TT == 0).AND.(ISOLVER == 0).AND.(NWCAV == 0)) CALL DGESL(DIJ(:,:,1),NPAN1,NPAN1,IPVT2,SI,0)
IF ((TT /= 0).AND.(ISOLVER == 0).AND.(NWCAV == 0)) CALL DGESL(DIJ(:,:,1),NPAN1,NPAN1,IPVT4,SI,0)
IF ((ISOLVER == 0).AND.(NWCAV /= 0)) CALL DGESL(DIJ(:,:,1),NPAN1,NPAN1,IPVTC,SI,0)
!-----------------------------------------------------------------------------------------------!
!    Iterative Solver                                                                           !
!-----------------------------------------------------------------------------------------------!
IF (ISOLVER == 1) CALL BISOF(DIJ(:,:,1),NPAN1,NPAN1,SI)
!-----------------------------------------------------------------------------------------------!
!    Structure the Potential Solution                                                           !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
POTP(:,:,TT)=RESHAPE(SI(NHPAN+1:NHPAN+NPPAN),(/NCP,NRP/))
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
IF (NCPW < 0) POTPWP(:,:,TT)=RESHAPE(SI(NHPAN+NPPAN+NNPAN+1:NPAN),(/IABS(NCPW),NRW/))
IF (NCPW > 0) POTPWS(:,:,TT)=RESHAPE(SI(NHPAN+NPPAN+NNPAN+1:NPAN),(/IABS(NCPW),NRW/))
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
POTN(:,:,TT)=RESHAPE(SI(NHPAN+NPPAN+1:NHPAN+NPPAN+NNPAN),(/NNXT1,NNTP/))
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
POTH(:,:,TT)=RESHAPE(SI(1:NHPAN),(/NHX,NHTP/))
!-----------------------------------------------------------------------------------------------!
!    Blade Velocities and Pressure                                                              !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   CALL VELP(TT)
   CALL PRESP(JJ,TT)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake Velocities                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   CALL VELPW(TT)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Velocities and Pressure                                                             !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   CALL VELN(TT)
   CALL PRESN(JJ,TT)
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub Velocities and Pressure                                                                !
!-----------------------------------------------------------------------------------------------!
IF (IH == 1) THEN
   CALL VELH(TT)
   CALL PRESH(JJ,TT)
END IF !(IH == 1)
!-----------------------------------------------------------------------------------------------!
!    Calculate the Blade Wake Strength at the Kutta Points                                      !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRW
   DPOTP(J)=SI(NPAN+J)
END DO !J=1,NRW
!-----------------------------------------------------------------------------------------------!
!    Pressure Jump at the Blade Trailing Edge                                                   !
!-----------------------------------------------------------------------------------------------!
DCP(:)=CPP(NCP,JI:JF,TT)-CPP(1,JI:JF,TT)
!-----------------------------------------------------------------------------------------------!
!    Calculate the Nozzle Wake Strength at the Kutta Points                                     !
!-----------------------------------------------------------------------------------------------!
DO J=1,NNTP
   DPOTN(J)=SI(NPAN+NRW+J)
END DO !J=1,NNTP
!-----------------------------------------------------------------------------------------------!
!    Pressure Jump at the Nozzle Trailing Edge                                                  !
!-----------------------------------------------------------------------------------------------!
CALL PRESNTE(TT,DCN)
!-----------------------------------------------------------------------------------------------!
!    Dipole Strength in the Blade Wake                                                          !
!-----------------------------------------------------------------------------------------------!
IF (TT == 0) THEN
   DO J=1,NRW
      POTPW(:,J, 0)=DPOTP(J)
   END DO !J=1,NRW
ELSE !(TT == 0)
   DO J=1,NRW
      POTPW(1,J,TT)=DPOTP(J)
   END DO !J=1,NRW
END IF !(TT == 0)
!-----------------------------------------------------------------------------------------------!
!    Dipole Strength in the Nozzle Wake                                                         !
!-----------------------------------------------------------------------------------------------!
IF (TT == 0) THEN
   DO J=1,NNTP
      POTNW(:,J, 0)=DPOTN(J)
   END DO !J=1,NNTP
ELSE !(TT == 0)
   DO J=1,NNTP
      POTNW(1,J,TT)=DPOTN(J)
   END DO !J=1,NNTP
END IF !(TT == 0)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE SOLVELKWET
!-----------------------------------------------------------------------------------------------!
