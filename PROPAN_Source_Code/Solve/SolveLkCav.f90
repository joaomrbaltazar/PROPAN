!-----------------------------------------------------------------------------------------------!
!    Solve Cavitation Linear Kutta Condition                                                    !
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
SUBROUTINE SOLVELKCAV(JJ,TT)
!-----------------------------------------------------------------------------------------------!
!    Created by: 16102014, J. Baltazar, Cavitation Model                                        !
!    Modified  : 28102014, J. Baltazar, Interpolation scheme at nozzle t.e.                     !
!    Modified  : 06112014, J. Baltazar, version 3.0                                             !
!    Modified  : 11112014, J. Baltazar, version 3.1                                             !
!    Modified  : 27112014, J. Baltazar, version 3.3, Super-Cavitation Model                     !
!    Modified  : 12122014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!    Modified  : 05032015, J. Baltazar, correction for super-cavitation                         !
!    Modified  : 27052016, J. Baltazar, 2016 version 1.2                                        !
!    Modified  : 02062016, J. Baltazar, 2016 version 1.3                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,J2,L,K,KB,JJ,TT,JP,JS
INTEGER :: INFO
!-----------------------------------------------------------------------------------------------!
!    Rewind Matrices                                                                            !
!-----------------------------------------------------------------------------------------------!
IF (TT == 0) THEN
   REWIND(23)
   READ  (23) ((DIJ(I,J,1),I=1,NPAN1),J=1,NPAN1)
ELSE !(TT)
   REWIND(26)
   READ  (26) ((DIJ(I,J,1),I=1,NPAN1),J=1,NPAN1)
END IF !(TT)
!-----------------------------------------------------------------------------------------------!
!    Complete System - Main Matrix                                                              !
!-----------------------------------------------------------------------------------------------!
!    Blade Panels                                                                               !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRP
   DO I=1,NCP
      K=NHPAN+(J-1)*NCP+I
!-----------------------------------------------------------------------------------------------!
      DO L=1,NPAN
         IF ((I <= IDP(J,TT)).AND.(I >= IRP(J,TT))) THEN
            DIJ(L,K,1)=-SIJ(L,K,1)
            IF (TT == 0) THEN
               DO KB=2,NB
                  DIJ(L,K,1)=DIJ(L,K,1)-SIJ(L,K,KB)
               END DO !KB=2,NB
            END IF !(TT == 0)
         ELSEIF ((I >= IDS(J,TT)).AND.(I <= IRS(J,TT))) THEN
            DIJ(L,K,1)=-SIJ(L,K,1)
            IF (TT == 0) THEN
               DO KB=2,NB
                  DIJ(L,K,1)=DIJ(L,K,1)-SIJ(L,K,KB)
               END DO !KB=2,NB
            END IF !(TT == 0)
         END IF
      END DO !L=1,NPAN
!-----------------------------------------------------------------------------------------------!
   END DO !I=1,NCP
END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
!    Blade Wake Panels                                                                          !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRW
   DO I=1,IABS(NCPW)
      K=NHPAN+NPPAN+NNPAN+(J-1)*IABS(NCPW)+I
!-----------------------------------------------------------------------------------------------!
      DO L=1,NPAN
         IF ((I >= IDPWP(J,TT)).AND.(I <= IRPWP(J,TT))) THEN
            DIJ(L,K,1)=-SIJ(L,K,1)
            IF (TT == 0) THEN
               DO KB=2,NB
                  DIJ(L,K,1)=DIJ(L,K,1)-SIJ(L,K,KB)
               END DO !KB=2,NB
            END IF !(TT == 0)
         ELSEIF ((I >= IDPWS(J,TT)).AND.(I <= IRPWS(J,TT))) THEN
            DIJ(L,K,1)=-SIJ(L,K,1)
            IF (TT == 0) THEN
               DO KB=2,NB
                  DIJ(L,K,1)=DIJ(L,K,1)-SIJ(L,K,KB)
               END DO !KB=2,NB
            END IF !(TT == 0)
         END IF
      END DO !L=1,NPAN
!-----------------------------------------------------------------------------------------------!
   END DO !I=1,IABS(NCPW)
END DO !J=1,NRW
!-----------------------------------------------------------------------------------------------!
!    Correction for Super-Cavitation in the Kutta Condition                                     !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRW
   J2=JI-1+J
   IF (IRP(J2,TT) == 1) THEN
      JP=NHPAN+(JI-1)*NCP+(J-1)*NCP+1
      DIJ(NPAN+J,JP,1)=0.D0
   END IF !(IRP(J2,TT) ==   1)
   IF (IRS(J2,TT) == NCP) THEN
      JS=NHPAN+(JI-1)*NCP+J*NCP
      DIJ(NPAN+J,JS,1)=0.D0
   END IF !(IRS(J2,TT) == NCP)
END DO !J=1,NRW
!-----------------------------------------------------------------------------------------------!
!    Solution of System of Equations                                                            !
!-----------------------------------------------------------------------------------------------!
IF (ISOLVER == 0) THEN
   INFO=0
   CALL DGEFA(DIJ(:,:,1),NPAN1,NPAN1,IPVTC2,INFO)
   IF (INFO /= 0) THEN
      WRITE(6,*) INFO
      STOP
   END IF !(INFO /= 0)
   CALL DGESL(DIJ(:,:,1),NPAN1,NPAN1,IPVTC2,SI,0)
ELSEIF (ISOLVER == 1) THEN
   CALL BISOF(DIJ(:,:,1),NPAN1,NPAN1,SI)
END IF !(ISOLVER)
!-----------------------------------------------------------------------------------------------!
!    Structure the Potential Solution                                                           !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRP
   DO I=1,NCP
      K=NHPAN+(J-1)*NCP+I
      IF ((I <= IDP(J,TT)).AND.(I >= IRP(J,TT))) THEN
         SOURCEP   (I,J, 1)=SI(K)
         SOURCEPCAV(I,J,TT)=SI(K)
      ELSEIF ((I >= IDS(J,TT)).AND.(I <= IRS(J,TT))) THEN
         SOURCEP   (I,J, 1)=SI(K)
         SOURCEPCAV(I,J,TT)=SI(K)
      ELSE
         POTP      (I,J,TT)=SI(K)
      END IF
   END DO !I=1,NCP
END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
IF (TT == 0) THEN
   DO KB=2,NB
      SOURCEP(:,:,KB)=SOURCEP(:,:,1)
   END DO !KB=2,NB
END IF !(TT == 0)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
DO J=1,NRW
   DO I=1,IABS(NCPW)
      K=NHPAN+NPPAN+NNPAN+(J-1)*IABS(NCPW)+I
      IF ((I >= IDPWP(J,TT)).AND.(I <= IRPWP(J,TT))) THEN
         SOURCEPWCAV(I,J,TT)=SI(K)
      ELSEIF ((I >= IDPWS(J,TT)).AND.(I <= IRPWS(J,TT))) THEN
         SOURCEPWCAV(I,J,TT)=SI(K)
      ELSEIF (NCPW < 0) THEN
         POTPWP     (I,J,TT)=SI(K)
         SOURCEPWCAV(I,J,TT)=0.D0
      ELSEIF (NCPW > 0) THEN
         POTPWS     (I,J,TT)=SI(K)
         SOURCEPWCAV(I,J,TT)=0.D0
      END IF
   END DO !I=1,IABS(NCPW)
END DO !J=1,NRW
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
END IF !(IW == 1)
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
END SUBROUTINE SOLVELKCAV
!-----------------------------------------------------------------------------------------------!
