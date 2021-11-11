!-----------------------------------------------------------------------------------------------!
!    Output Files                                                                               !
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
SUBROUTINE OUTPUT(IFILE,JJ,TT,CC)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 03122013, J. Baltazar, version 1.0                                             !
!    Modified  : 27052014, J. Baltazar, Wake Alignment Module                                   !
!    Modified  : 03112014, J. Baltazar, version 3.0, Steady Cavitation Model                    !
!    Modified  : 12112014, J. Baltazar, version 3.1, Unsteady Cavitation Model                  !
!    Modified  : 18112014, J. Baltazar, version 3.2, Mid-chord Cavitation                       !
!    Modified  : 27112014, J. Baltazar, version 3.3, Super-Cavitation Model                     !
!    Modified  : 09122014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!    Modified  : 10022016, J. Baltazar, 2016 version 1.0                                        !
!    Modified  : 18042019, J. Baltazar, 2019 version 1.0                                        !
!    Modified  : 25052020, J. Baltazar, 2020 version 1.1, Output corrected                      !
!    Modified  : 23062020, J. Baltazar, 2020 version 1.2, Dynamic inflow                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
CHARACTER OUTFILE*50
LOGICAL :: FCHECK
INTEGER :: I,J,J2,K,KB,JJ,TT,CC,IFILE
DOUBLE PRECISION :: PHI
!-----------------------------------------------------------------------------------------------!
WRITE(ICHAR,'(I3)') TT
SELECT CASE (IFILE)
!-----------------------------------------------------------------------------------------------!
CASE(31)
!-----------------------------------------------------------------------------------------------!
!    Write Geometry in Tecplot Format                                                           !
!-----------------------------------------------------------------------------------------------!
   IF (NWA > 0) THEN
      OUTFILE='GEOMETRY_'//IDENTU(JJ)//'_T'//TRIM(ADJUSTL(ICHAR))//'.DAT'
   ELSE
      OUTFILE='GEOMETRY_'//IDENTU(JJ)//'.DAT'
   END IF
   OPEN(UNIT=31,FILE=OUTFILE,STATUS='REPLACE')
   WRITE(31,'(A)') 'TITLE = "GEOMETRY"'
   WRITE(31,'(A)') 'VARIABLES="X0" "Y0" "Z0" "UNX0" "UNY0" "UNZ0"'// &
                      ' "ET1X" "ET1Y" "ET1Z" "ET2X" "ET2Y" "ET2Z"'// &
                      ' "AT1X" "AT1Y" "AT1Z" "AT2X" "AT2Y" "AT2Z" "A0"'
!-----------------------------------------------------------------------------------------------!
!    Blade Geometry                                                                             !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      WRITE(31,200) 'ZONE T="BLADE" F=POINT, I= ',NCP,' J= ',NRP
      DO J=1,NRP
         DO I=1,NCP
            WRITE(31,210) XP0(I,J),YP0(I,J),ZP0(I,J),UNXP0(I,J),UNYP0(I,J),UNZP0(I,J), &
                    ET1XP(I,J),ET1YP(I,J),ET1ZP(I,J),ET2XP(I,J),ET2YP(I,J),ET2ZP(I,J), &
                    AT1XP(I,J),AT1YP(I,J),AT1ZP(I,J),AT2XP(I,J),AT2YP(I,J),AT2ZP(I,J),AP0(I,J)
         END DO !I=1,NCP
      END DO !J=1,NRP
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Geometry                                                                            !
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      WRITE(31,200) 'ZONE T="NOZZLE" F=POINT, I= ',NNXT1,' J= ',NNTP
      DO J=1,NNTP
         DO I=1,NNXT1
            WRITE(31,210) XN0(I,J),YN0(I,J),ZN0(I,J),UNXN0(I,J),UNYN0(I,J),UNZN0(I,J), &
                    ET1XN(I,J),ET1YN(I,J),ET1ZN(I,J),ET2XN(I,J),ET2YN(I,J),ET2ZN(I,J), &
                    AT1XN(I,J),AT1YN(I,J),AT1ZN(I,J),AT2XN(I,J),AT2YN(I,J),AT2ZN(I,J),AN0(I,J)
         END DO !I=1,NNXT1
      END DO !J=1,NNTP
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub Geometry                                                                               !
!-----------------------------------------------------------------------------------------------!
   IF (IABS(IH) == 1) THEN
      WRITE(31,200) 'ZONE T="HUB" F=POINT, I= ',NHX,' J= ',NHTP
      DO J=1,NHTP
         DO I=1,NHX
            WRITE(31,210) XH0(I,J),YH0(I,J),ZH0(I,J),UNXH0(I,J),UNYH0(I,J),UNZH0(I,J), &
                    ET1XH(I,J),ET1YH(I,J),ET1ZH(I,J),ET2XH(I,J),ET2YH(I,J),ET2ZH(I,J), &
                    AT1XH(I,J),AT1YH(I,J),AT1ZH(I,J),AT2XH(I,J),AT2YH(I,J),AT2ZH(I,J),AH0(I,J)
         END DO !I=1,NHX
      END DO !J=1,NHTP
   END IF !(IABS(IH) == 1)
   CLOSE(UNIT=31)
!-----------------------------------------------------------------------------------------------!
CASE(32)
!-----------------------------------------------------------------------------------------------!
!    Write Wake Solution in Tecplot Format                                                      !
!-----------------------------------------------------------------------------------------------!
   IF (TT == 0) THEN
      OUTFILE='CIRC_'//IDENTU(JJ)//'.DAT'
      OPEN(UNIT=32,FILE=OUTFILE,STATUS='REPLACE')
      WRITE(32,'(A)') 'TITLE = "PROPAN_'//IDENTU(JJ)//'"'
      WRITE(32,'(A)') 'VARIABLES="r/R","POTW","POT(1)","POT(N)","CP(1)","CP(N)","DELTACP"'
   END IF !(TT == 0)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      WRITE(32,205) 'ZONE T="BLADE WAKE T='//TRIM(ADJUSTL(ICHAR))//'", I= ',NRW,' F=POINT'
      DO J=1,NRW
         J2=JI-1+J
         WORK1=0.5D0*(DSQRT(YP0(1,J2)**2+ZP0(1,J2)**2)+DSQRT(YP0(NCP,J2)**2+ZP0(NCP,J2)**2))
         WRITE(32,215) WORK1,POTPW(1,J,TT),POTP(1,J2,TT),POTP(NCP,J2,TT), &
                       CPP(1,J2,TT),CPP(NCP,J2,TT),DCP(J)
      END DO !J=1,NRW
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      WRITE(32,205) 'ZONE T="NOZZLE WAKE T='//TRIM(ADJUSTL(ICHAR))//'", I= ',NNTP,' F=POINT'
      DO J=NNT+1,NNTP
         WORK1=0.5D0*(DATAN2D(ZN0(NNX1,J),YN0(NNX1,J))+DATAN2D(ZN0(NNX,J),YN0(NNX,J)))
         WORK1=WORK1-0.5D0*(DATAN2D(ZN(NNX1,NNT+1),YN(NNX1,NNT+1))+ &
                            DATAN2D(ZN(NNX ,NNT+1),YN(NNX ,NNT+1)))
         IF (WORK1 <   0.D0) WORK1=WORK1+360.D0
         IF (WORK1 > 360.D0) WORK1=WORK1-360.D0
         WRITE(32,215) WORK1,POTNW(1,J,TT),POTN(NNX1,J,TT),POTN(NNX,J,TT), &
                       CPN(NNX1,J,TT),CPN(NNX,J,TT),DCN(J)
      END DO !J=NNT+1,NNTP
      DO J=1,NNT
         WORK1=0.5D0*(DATAN2D(ZN0(NNX1,J),YN0(NNX1,J))+DATAN2D(ZN0(NNX,J),YN0(NNX,J)))
         WORK1=WORK1-0.5D0*(DATAN2D(ZN(NNX1,NNT+1),YN(NNX1,NNT+1))+ &
                            DATAN2D(ZN(NNX ,NNT+1),YN(NNX ,NNT+1)))
         WORK1=WORK1+360.D0/DFLOAT(NB)
         IF (WORK1 <   0.D0) WORK1=WORK1+360.D0
         IF (WORK1 > 360.D0) WORK1=WORK1-360.D0
         WRITE(32,215) WORK1,POTNW(1,J,TT),POTN(NNX1,J,TT),POTN(NNX,J,TT), &
                       CPN(NNX1,J,TT),CPN(NNX,J,TT),DCN(J)
      END DO !J=1,NNT
   END IF !(IN == 1)
   IF (TT == NT) CLOSE(UNIT=32)
!-----------------------------------------------------------------------------------------------!
CASE(33)
!-----------------------------------------------------------------------------------------------!
!    Write Solution in Tecplot Format                                                           !
!-----------------------------------------------------------------------------------------------!
   IF (TT == 0) THEN
      OUTFILE='SOL_'//IDENTU(JJ)//'.DAT'
      OPEN(UNIT=33,FILE=OUTFILE,STATUS='REPLACE')
      WRITE(33,'(A)') 'TITLE = "PROPAN_'//IDENTU(JJ)//'"'
      WRITE(33,'(A)') 'VARIABLES="X0" "Y0" "Z0" "POT" "SOURCE" "VT1" "VT2"'// &
                                                ' "VTT1" "VTT2" "CP" "CPN"'
   END IF !(TT == 0)
!-----------------------------------------------------------------------------------------------!
!    Blade Solution                                                                             !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      WRITE(33,200) 'ZONE T="BLADE T='//TRIM(ADJUSTL(ICHAR))//'" F=POINT, I= ',NCP,' J= ',NRP
      DO J=1,NRP
         DO I=1,NCP
            WRITE(33,220) XP0(I,J),YP0(I,J),ZP0(I,J),POTP(I,J,TT),SOURCEP(I,J,1), &
                          VT1P(I,J),VT2P(I,J),VTT1P(I,J),VTT2P(I,J),CPP(I,J,TT),CPNP(I,J,TT)
         END DO !I=1,NCP
      END DO !J=1,NRP
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake Solution                                                                        !
!-----------------------------------------------------------------------------------------------!
   IF ((IP == 1).AND.(NCPW < 0)) THEN
      WRITE(33,200) 'ZONE T="BLADE WAKE T='//TRIM(ADJUSTL(ICHAR))// &
                                                          '" F=POINT, I= ',IABS(NCPW),' J= ',NRW
      DO J=1,NRW
         DO I=1,IABS(NCPW)
            WRITE(33,220) XPW0(I,J),YPW0(I,J),ZPW0(I,J),POTPWP(I,J,TT), &
                          SOURCEPWCAV(I,J,TT),VT1PW(I,J),0.D0,0.D0,0.D0,0.D0,0.D0
         END DO !I=1,IABS(NCPW)
      END DO !J=1,NRW
   END IF !((IP == 1).AND.(NCPW < 0))
   IF ((IP == 1).AND.(NCPW > 0)) THEN
      WRITE(33,200) 'ZONE T="BLADE WAKE T='//TRIM(ADJUSTL(ICHAR))// &
                                                          '" F=POINT, I= ',IABS(NCPW),' J= ',NRW
      DO J=1,NRW
         DO I=1,IABS(NCPW)
            WRITE(33,220) XPW0(I,J),YPW0(I,J),ZPW0(I,J),POTPWS(I,J,TT), &
                          SOURCEPWCAV(I,J,TT),VT1PW(I,J),0.D0,0.D0,0.D0,0.D0,0.D0
         END DO !I=1,IABS(NCPW)
      END DO !J=1,NRW
   END IF !((IP == 1).AND.(NCPW > 0))
!-----------------------------------------------------------------------------------------------!
!    Nozzle Solution                                                                            !
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      WRITE(33,200) 'ZONE T="NOZZLE T='//TRIM(ADJUSTL(ICHAR))// &
                                                              '" F=POINT, I= ',NNXT1,' J= ',NNTP
      DO J=1,NNTP
         DO I=1,NNXT1
            WRITE(33,220) XN0(I,J),YN0(I,J),ZN0(I,J),POTN(I,J,TT),SOURCEN(I,J,1), &
                          VT1N(I,J),VT2N(I,J),VTT1N(I,J),VTT2N(I,J),CPN(I,J,TT),CPNN(I,J,TT)
         END DO !I=1,NNXT1
      END DO !J=1,NNTP
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub Solution                                                                               !
!-----------------------------------------------------------------------------------------------!
   IF (IABS(IH) == 1) THEN
      WRITE(33,200) 'ZONE T="HUB T='//TRIM(ADJUSTL(ICHAR))//'" F=POINT, I= ',NHX,' J= ',NHTP
      DO J=1,NHTP
         DO I=1,NHX
            WRITE(33,220) XH0(I,J),YH0(I,J),ZH0(I,J),POTH(I,J,TT),SOURCEH(I,J,1), &
                          VT1H(I,J),VT2H(I,J),VTT1H(I,J),VTT2H(I,J),CPH(I,J,TT),CPNH(I,J,TT)
         END DO !I=1,NHX
      END DO !J=1,NHTP
   END IF !(IABS(IH) == 1)
   IF (TT == NT) CLOSE(UNIT=33)
!-----------------------------------------------------------------------------------------------!
CASE(34)
!-----------------------------------------------------------------------------------------------!
!    Write Blade Forces in Tecplot Format                                                       !
!-----------------------------------------------------------------------------------------------!
   IF (TT == 0) THEN
      OUTFILE='FORCES_BLADE_'//IDENTU(JJ)//'.DAT'
      OPEN(UNIT=34,FILE=OUTFILE,STATUS='REPLACE')
      WRITE(34,'(A)') 'TITLE = "PROPAN_'//IDENTU(JJ)//'"'
      IF (IROTOR == -1) &
            WRITE(34,'(A)') 'VARIABLES="T" "KB" "CFX" "CMX" "CFY" "CMY" "CFZ" "CMZ"'// &
                            ' "CTX" "CQX" "CTY" "CQY" "CTZ" "CQZ"'
      IF (IROTOR == 0) &
            WRITE(34,'(A)') 'VARIABLES="T" "KB" "KTX" "KQX" "KTY" "KQY" "KTZ" "KQZ"'// &
                            '"KTTX" "KQTX" "KTTY" "KQTY" "KTTZ" "KQTZ"'
      IF (IROTOR == 1) &
            WRITE(34,'(A)') 'VARIABLES="T" "KB" "CTX" "CPX" "CTY" "CPY" "CTZ" "CPZ"'// &
                            '"CTTX" "CPTX" "CTTY" "CPTY" "CTTZ" "CPTZ"'
   END IF !(TT == 0)
   DO KB=1,NB
      WRITE(34,225) TT,KB,CFXP(TT,KB),CMXP(TT,KB),CFYP(TT,KB),CMYP(TT,KB), &
                    CFZP(TT,KB),CMZP(TT,KB),CTXP(TT),CQXP(TT),CTYP(TT), &
                    CQYP(TT),CTZP(TT),CQZP(TT)
   END DO !KB=1,NB
   IF (TT == NT) CLOSE(UNIT=34)
!-----------------------------------------------------------------------------------------------!
CASE(35)
!-----------------------------------------------------------------------------------------------!
!    Write Nozzle Forces in Tecplot Format                                                      !
!-----------------------------------------------------------------------------------------------!
   IF (TT == 0) THEN
      OUTFILE='FORCES_NOZZLE_'//IDENTU(JJ)//'.DAT'
      OPEN(UNIT=35,FILE=OUTFILE,STATUS='REPLACE')
      WRITE(35,'(A)') 'TITLE = "PROPAN_'//IDENTU(JJ)//'"'
      IF (IROTOR == -1) &
            WRITE(35,'(A)') 'VARIABLES="T" "KB" "CFX" "CMX" "CFY" "CMY" "CFZ" "CMZ"'// &
                            ' "CTX" "CQX" "CTY" "CQY" "CTZ" "CQZ"'
      IF (IROTOR == 0) &
            WRITE(35,'(A)') 'VARIABLES="T" "KB" "KTX" "KQX" "KTY" "KQY" "KTZ" "KQZ"'// &
                            '"KTTX" "KQTX" "KTTY" "KQTY" "KTTZ" "KQTZ"'
      IF (IROTOR == 1) &
            WRITE(35,'(A)') 'VARIABLES="T" "KB" "CTX" "CPX" "CTY" "CPY" "CTZ" "CPZ"'// &
                            '"CTTX" "CPTX" "CTTY" "CPTY" "CTTZ" "CPTZ"'
   END IF !(TT == 0)
   DO KB=1,NB
      WRITE(35,225) TT,KB,CFXN(TT,KB),CMXN(TT,KB),CFYN(TT,KB),CMYN(TT,KB), &
                    CFZN(TT,KB),CMZN(TT,KB),CTXN(TT),CQXN(TT),CTYN(TT), &
                    CQYN(TT),CTZN(TT),CQZN(TT)
   END DO !KB=1,NB
   IF (TT == NT) CLOSE(UNIT=35)
!-----------------------------------------------------------------------------------------------!
CASE(36)
!-----------------------------------------------------------------------------------------------!
!    Write Hub Forces in Tecplot Format                                                         !
!-----------------------------------------------------------------------------------------------!
   IF (TT == 0) THEN
      OUTFILE='FORCES_HUB_'//IDENTU(JJ)//'.DAT'
      OPEN(UNIT=36,FILE=OUTFILE,STATUS='REPLACE')
      WRITE(36,'(A)') 'TITLE = "PROPAN_'//IDENTU(JJ)//'"'
      IF (IROTOR == -1) &
            WRITE(36,'(A)') 'VARIABLES="T" "KB" "CFX" "CMX" "CFY" "CMY" "CFZ" "CMZ"'// &
                            ' "CTX" "CQX" "CTY" "CQY" "CTZ" "CQZ"'
      IF (IROTOR == 0) &
            WRITE(36,'(A)') 'VARIABLES="T" "KB" "KTX" "KQX" "KTY" "KQY" "KTZ" "KQZ"'// &
                            '"KTTX" "KQTX" "KTTY" "KQTY" "KTTZ" "KQTZ"'
      IF (IROTOR == 1) &
            WRITE(36,'(A)') 'VARIABLES="T" "KB" "CTX" "CPX" "CTY" "CPY" "CTZ" "CPZ"'// &
                            '"CTTX" "CPTX" "CTTY" "CPTY" "CTTZ" "CPTZ"'
   END IF !(TT == 0)
   DO KB=1,NB
      WRITE(36,225) TT,KB,CFXH(TT,KB),CMXH(TT,KB),CFYH(TT,KB),CMYH(TT,KB), &
                    CFZH(TT,KB),CMZH(TT,KB),CTXH(TT),CQXH(TT),CTYH(TT), &
                    CQYH(TT),CTZH(TT),CQZH(TT)
   END DO !KB=1,NB
   IF (TT == NT) CLOSE(UNIT=36)
!-----------------------------------------------------------------------------------------------!
CASE(37)
!-----------------------------------------------------------------------------------------------!
!    Write Velocity Field in Tecplot Format                                                     !
!-----------------------------------------------------------------------------------------------!
   IF (TT == 0) THEN
      OUTFILE='VELF_'//IDENTU(JJ)//'.DAT'
      OPEN(UNIT=37,FILE=OUTFILE,STATUS='REPLACE')
      WRITE(37,'(A)') 'TITLE = "PROPAN_'//IDENTU(JJ)//'"'
      WRITE(37,'(A)') 'VARIABLES="T","X","Y","Z","VX","VY","VZ","POT","DPOTDT"'
   END IF !(TT == 0)
   DO J=1,NFY
      DO I=1,NFX
         WRITE(37,230) TT,XVF(I,J),YVF(I,J),ZVF(I,J),VFX(I,J),VFY(I,J),VFZ(I,J), &
                                                                       POTF(I,J,TT),DPOTFDT(I,J)
      END DO !I=1,NFX
   END DO !J=1,NFY
   IF (TT == NT) CLOSE(UNIT=37)
!-----------------------------------------------------------------------------------------------!
CASE(38)
!-----------------------------------------------------------------------------------------------!
!    Write Panel Grid in Tecplot Format                                                         !
!-----------------------------------------------------------------------------------------------!
   OUTFILE='PROPANEL_'//IDENTU(JJ)//'_T'//TRIM(ADJUSTL(ICHAR))//'.DAT'
   OPEN(UNIT=38,FILE=OUTFILE,STATUS='REPLACE')
   WRITE(38,'(A)') ' TITLE="PROPELLER"'
   WRITE(38,'(A)') ' VARIABLES= "X" "Y" "Z" '
!-----------------------------------------------------------------------------------------------!
!    Write Blade Grid in Tecplot Format                                                         !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      K=1
      WRITE(38,206) ' ZONE T="BLADE',K,'" F=POINT, I= ',NCP1,' J= ',NRP1
      DO J=1,NRP1
         DO I=1,NCP1
            WRITE(38,235) XP(I,J),YP(I,J),ZP(I,J)
         END DO !I=1,NCP1
      END DO !J=1,NRP1
      DO K=2,NB
         PHI=DFLOAT(K-1)/DFLOAT(NB)*2.D0*PI
         WRITE(38,206) ' ZONE T="BLADE',K,'" F=POINT, I=',NCP1,' J=',NRP1
         DO J=1,NRP1
            DO I=1,NCP1
               WRITE(38,235) XP(I,J),RP(I,J)*DCOS(TP(I,J)+PHI),RP(I,J)*DSIN(TP(I,J)+PHI)
            END DO !I=1,NCP1
         END DO !J=1,NRP1
      END DO !K=2,NB
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Write Blade Wake Grid in Tecplot Format                                                    !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      K=1
      WRITE(38,206) ' ZONE T="BLADE WAKE',K,'" F=POINT, I= ',NPW1,' J= ',NRW1
      DO J=1,NRW1
         DO I=1,NPW1
            WRITE(38,235) XPW(I,J),YPW(I,J),ZPW(I,J)
         END DO !I=1,NPW1
      END DO !J=1,NRW1
      DO K=2,NB
         PHI=DFLOAT(K-1)/DFLOAT(NB)*2.D0*PI
         WRITE(38,206) ' ZONE T="BLADE WAKE',K,'" F=POINT, I= ',NPW1,' J= ',NRW1
         DO J=1,NRW1
            DO I=1,NPW1
               WRITE(38,235) XPW(I,J),RPW(I,J)*DCOS(TPW(I,J)+PHI),RPW(I,J)*DSIN(TPW(I,J)+PHI)
            END DO !I=1,NPW1
         END DO !J=1,NRW1
      END DO !K=2,NB
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Write Nozzle Grid in Tecplot Format                                                        !
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      K=1
      WRITE(38,206) ' ZONE T="NOZZLE',K,'" F=POINT, I= ',NNXT,' J= ',NNTT
      DO J=1,NNTT
         DO I=1,NNXT
            WRITE(38,235) XN(I,J),YN(I,J),ZN(I,J)
         END DO !I=1,NNXT
      END DO !J=1,NNTT
      DO K=2,NB
         PHI=DFLOAT(K-1)/DFLOAT(NB)*2.D0*PI
         WRITE(38,206) ' ZONE T="NOZZLE',K,'" F=POINT, I= ',NNXT,' J= ',NNTT
         DO J=1,NNTT
            DO I=1,NNXT
               WRITE(38,235) XN(I,J),RN(I,J)*DCOS(TN(I,J)+PHI),RN(I,J)*DSIN(TN(I,J)+PHI)
            END DO !I=1,NNXT
         END DO !J=1,NNTT
      END DO !K=2,NB
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Write Nozzle Wake Grid in Tecplot Format                                                   !
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      K=1
      WRITE(38,206) ' ZONE T="NOZZLE WAKE',K,'" F=POINT, I= ',NNW1,' J= ',NNTT
      DO J=1,NNTT
         DO I=1,NNW1
            WRITE(38,235) XNW(I,J),YNW(I,J),ZNW(I,J)
         END DO !I=1,NNW1
      END DO !J=1,NNTT
      DO K=2,NB
         PHI=DFLOAT(K-1)/DFLOAT(NB)*2.D0*PI
         WRITE(38,206) ' ZONE T="NOZZLE WAKE',K,'" F=POINT, I= ',NNW1,' J= ',NNTT
         DO J=1,NNTT
            DO I=1,NNW1
               WRITE(38,235) XNW(I,J),RNW(I,J)*DCOS(TNW(I,J)+PHI),RNW(I,J)*DSIN(TNW(I,J)+PHI)
            END DO !I=1,NNW1
         END DO !J=1,NNTT
      END DO !K=2,NB
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Write Hub Grid in Tecplot Format                                                           !
!-----------------------------------------------------------------------------------------------!
   IF (IABS(IH) == 1) THEN
      K=1
      WRITE(38,206) ' ZONE T="HUB',K,'" F=POINT, I= ',NHX1,' J= ',NHTT
      DO J=1,NHT1
         DO I=1,NHX1
            WRITE(38,235) XH(I,J),YH(I,J),ZH(I,J)
         END DO !I=1,NHX1
      END DO !J=1,NHT1
      DO J=(NHT1+1),NHTT
         DO I=1,NHX1
            WRITE(38,235) XH(I,J),YH(I,J),ZH(I,J)
         END DO !I=1,NHX1
      END DO !J=(NHT1+1),NHTT
      DO K=2,NB
         PHI=DFLOAT(K-1)/DFLOAT(NB)*2.D0*PI
         WRITE(38,206) ' ZONE T="HUB',K,'" F=POINT, I= ',NHX1,' J= ',NHTT
         DO J=1,NHT1
            DO I=1,NHX1
               WRITE(38,235) XH(I,J),RH(I,J)*DCOS(TH(I,J)+PHI),RH(I,J)*DSIN(TH(I,J)+PHI)
            END DO !I=1,NHX1
         END DO !J=1,NHT1
         DO J=(NHT1+1),NHTT
            DO I=1,NHX1
               WRITE(38,235) XH(I,J),RH(I,J)*DCOS(TH(I,J)+PHI),RH(I,J)*DSIN(TH(I,J)+PHI)
            END DO !I=1,NHX1
         END DO !J=(NHT1+1),NHTT
      END DO !K=2,NB
   END IF !(IABS(IH) == 1)
   CLOSE(UNIT=38)
!-----------------------------------------------------------------------------------------------!
CASE(39)
!-----------------------------------------------------------------------------------------------!
!    Write Cavity                                                                               !
!-----------------------------------------------------------------------------------------------!
   INQUIRE(UNIT=39,OPENED=FCHECK)
   IF (.NOT.FCHECK) THEN
      OUTFILE='ITERCAV_'//IDENTU(JJ)//'.DAT'
      OPEN(UNIT=39,FILE=OUTFILE,STATUS='REPLACE')
      WRITE(39,'(A)') ' TITLE="ITERCAV_'//IDENTU(JJ)//'"'
      IF (IROTOR == 0) WRITE(39,'(A)') ' VARIABLES="CC","LCMAXP","RACAEP","ETAP"'// &
                                                      ',"LCMAXS","RACAES","ETAS"'
      IF (IROTOR == 1) WRITE(39,'(A)') ' VARIABLES="CC","LCMAXS","RACAES","ETAS"'// &
                                                      ',"LCMAXP","RACAEP","ETAP"'
   END IF !(.NOT.FCHECK)
   IF (CC == 0) WRITE(39,'(A)') ' ZONE T="ITERCAV T='//TRIM(ADJUSTL(ICHAR))//'"'
   WRITE(39,240) CC,LCMAXP,RACAEP,ETAP,LCMAXS,RACAES,ETAS
!-----------------------------------------------------------------------------------------------!
CASE(40)
!-----------------------------------------------------------------------------------------------!
!    Write Cavity in Tecplot Format                                                             !
!-----------------------------------------------------------------------------------------------!
   IF (TT == 0) THEN
      OUTFILE='CAV_'//IDENTU(JJ)//'.DAT'
      OPEN(UNIT=40,FILE=OUTFILE,STATUS='REPLACE')
      WRITE(40,'(A)') ' TITLE="CAV"'
      WRITE(40,'(A)') ' VARIABLES= "X" "Y" "Z" "CAV"'
   END IF !(TT == 0)
!-----------------------------------------------------------------------------------------------!
!    Write Blade Cavity in Tecplot Format                                                       !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      K=1
      WRITE(40,206) ' ZONE T="BLADE',K,' T='//TRIM(ADJUSTL(ICHAR))// &
                                                                 '" F=POINT, I= ',NCP,' J= ',NRP
      DO J=1,NRP
         DO I=1,NCP
            WRITE(40,245) XP0(I,J),YP0(I,J),ZP0(I,J),THICKP(I,J,TT)
         END DO !I=1,NCP
      END DO !J=1,NRP
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Write Blade Wake Cavity in Tecplot Format                                                  !
!-----------------------------------------------------------------------------------------------!
   IF ((IP == 1).AND.(NCPW < 0)) THEN
      K=1
      WRITE(40,206) ' ZONE T="BLADE WAKE',K,' T='//TRIM(ADJUSTL(ICHAR))// &
                                                          '" F=POINT, I= ',IABS(NCPW),' J= ',NRW
      DO J=1,NRW
         DO I=1,IABS(NCPW)
            WRITE(40,245) XPW0(I,J),YPW0(I,J),ZPW0(I,J),CAMBERPWP(I,J,TT)+0.5D0*THICKPWP(I,J,TT)
         END DO !I=1,IABS(NCPW)
      END DO !J=1,NRW
   END IF !((IP == 1).AND.(NCPW < 0))
   IF ((IP == 1).AND.(NCPW > 0)) THEN
      K=1
      WRITE(40,206) ' ZONE T="BLADE WAKE',K,' T='//TRIM(ADJUSTL(ICHAR))// &
                                                          '" F=POINT, I= ',IABS(NCPW),' J= ',NRW
      DO J=1,NRW
         DO I=1,IABS(NCPW)
            WRITE(40,245) XPW0(I,J),YPW0(I,J),ZPW0(I,J),CAMBERPWS(I,J,TT)+0.5D0*THICKPWS(I,J,TT)
         END DO !I=1,IABS(NCPW)
      END DO !J=1,NRW
   END IF !((IP == 1).AND.(NCPW > 0))
   IF (TT == NT) CLOSE(UNIT=40)
!-----------------------------------------------------------------------------------------------!
CASE(41)
!-----------------------------------------------------------------------------------------------!
!    Write Cavity Detachment and Reattachment Positions                                         !
!-----------------------------------------------------------------------------------------------!
   INQUIRE(UNIT=41,OPENED=FCHECK)
   IF (.NOT.FCHECK) THEN
      OUTFILE='IPOSCAV_'//IDENTU(JJ)//'.DAT'
      OPEN(UNIT=41,FILE=OUTFILE,STATUS='REPLACE')
      WRITE(41,  205) ' TITLE="IPOSCAV_'//IDENTU(JJ)//'"'
      IF (IROTOR == 0) WRITE(41,'(A)') ' VARIABLES="J","IDP","IRP","IDPWP","IRPWP"'// &
                                                     ',"IDS","IRS","IDPWS","IRPWS"'
      IF (IROTOR == 1) WRITE(41,'(A)') ' VARIABLES="J","IDS","IRS","IDPWS","IRPWS"'// &
                                                     ',"IDP","IRP","IDPWP","IRPWP"'
   END IF !(.NOT.FCHECK)
   WRITE(41,  205) ' ZONE T="CAV ITER= ',CC,' T= '//TRIM(ADJUSTL(ICHAR))//'"'
   DO J=1,NRP
      WRITE(41,250) J,IDP(J,TT),IRP(J,TT),IDPWP(J,TT),IRPWP(J,TT), &
                      IDS(J,TT),IRS(J,TT),IDPWS(J,TT),IRPWS(J,TT)
   END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
END SELECT !(IFILE)
!-----------------------------------------------------------------------------------------------!
!    Formats                                                                                    !
!-----------------------------------------------------------------------------------------------!
200 FORMAT(A,I5,A,I5)
205 FORMAT(A,I5,A)
206 FORMAT(A,I5,A,I5,A,I5)
210 FORMAT(19(2X,E23.16))
215 FORMAT( 7(2X,E23.16))
220 FORMAT(11(2X,E23.16))
225 FORMAT( 2(2X,I5),13(2X,E23.16))
230 FORMAT( 1(2X,I5),8(2X,E23.16))
235 FORMAT( 3(2X,E23.16))
240 FORMAT( 1(2X,I5),6(2X,E23.16))
245 FORMAT( 4(2X,E23.16))
250 FORMAT( 9(2X,I5))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE OUTPUT
!-----------------------------------------------------------------------------------------------!
