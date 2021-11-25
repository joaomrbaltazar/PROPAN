!-----------------------------------------------------------------------------------------------!
!    Two-dimensional vsicous corrections                                                        !
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
SUBROUTINE VISCOR(JJ)
!-----------------------------------------------------------------------------------------------!
!    Created by: 22052014, J. Baltazar, 2014 version 1.0                                        !
!    Modified  : 23052014, J. Baltazar, 2014 version 1.0                                        !
!    Modified  : 26012015, J. Baltazar, 2015 version 1.0 Suction Force Correction               !
!    Modified  : 09032015, J. Baltazar, 2015 version 1.1 Unsteady Flow                          !
!    Modified  : 07072017, J. Baltazar, 2017 version 1.0                                        !
!    Modified  : 20032019, J. Baltazar, 2019 version 1.0                                        !
!    Modified  : 20042020, J. Baltazar, 2020 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPOST_MOD
IMPLICIT NONE
CHARACTER*50  FILEINP,OUTFILE
CHARACTER*3   ICHAR
INTEGER :: I,J,K,JJ,KK,TT
DOUBLE PRECISION :: R,R1,R2,DR,XC(4),YC(4),ZC(4)
DOUBLE PRECISION :: BETA,BETAI,PITCH,PITCH1,PITCH2,ALPHA
DOUBLE PRECISION :: XLE,XTE,TLE,TTE,YBLE,YBTE,C
DOUBLE PRECISION :: XC0,YC0,ZC0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX,UNY,UNZ,A0
DOUBLE PRECISION :: FX,FY,FZ,MX,MY,MZ,FIX,MIX,FDX,MDX,FVX,MVX,FA,FS,MS,FN,FT
DOUBLE PRECISION :: CLF,CLI,CLV,CDI,CDV,VE,FL,CIRC,RES,CFNI,CFTI,CFND,CFTD,CFNV,CFTV
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: CFIP,CMIP,CFDP,CMDP,CFVP,CMVP
DOUBLE PRECISION :: AN1,BN1,AN2,BN2,AN3,BN3,AN4,BN4,AN5,BN5,AN6,BN6
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: AMN1,BMN1,AMPLI1,PHASE1
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: AMN2,BMN2,AMPLI2,PHASE2
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: AMN3,BMN3,AMPLI3,PHASE3
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: AMN4,BMN4,AMPLI4,PHASE4
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: AMN5,BMN5,AMPLI5,PHASE5
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: AMN6,BMN6,AMPLI6,PHASE6
!-----------------------------------------------------------------------------------------------!
!    Write Forces                                                                               !
!-----------------------------------------------------------------------------------------------!
ALLOCATE(CFIP(0:NT),CMIP(0:NT),CFDP(0:NT),CMDP(0:NT),CFVP(0:NT),CMVP(0:NT))
OUTFILE='FORCES_'//IDENTU(JJ)//'.DAT'
OPEN(UNIT=30,FILE=OUTFILE,STATUS='UNKNOWN')
WRITE(30,'(A)') 'TITLE="Rotor Forces"'
WRITE(30,'(A)') 'VARIABLES="T" "CFIP" "CMIP" "CFDP" "CMDP" "CFVP" "CMVP"'
WRITE(30,'(A,I4)') 'ZONE T="BLADE" F=POINT, I=',(NT+1)
!-----------------------------------------------------------------------------------------------!
!    Write Solution Data in TecPlot Format                                                      !
!-----------------------------------------------------------------------------------------------!
OUTFILE='SECTION_'//IDENTU(JJ)//'.DAT'
OPEN(UNIT=33,FILE=OUTFILE,STATUS='UNKNOWN')
DO TT=0,NT
   WRITE(ICHAR,'(I3)') TT
   WRITE(33,*) 'TITLE = "PROPAN_'//IDENTU(JJ)//'_T'//TRIM(ADJUSTL(ICHAR))//'"'
   WRITE(33,'(A)') &
          'VARIABLES="I","R","C/R","PITCH","BETA","BETAI","ALPHA","RES","CLI","CDI",'// &
                    '"CLV","CDV","CFNI","CFTI","CFND","CFTD","CFNV","CFTV"'
!-----------------------------------------------------------------------------------------------!
!    Viscous Correction                                                                         !
!-----------------------------------------------------------------------------------------------!
   FIX=0.D0
   MIX=0.D0
   FDX=0.D0
   MDX=0.D0
   FVX=0.D0
   MVX=0.D0
   DO J=1,JF
!-----------------------------------------------------------------------------------------------!
!    Normal Direction Defined by Chord                                                          !
!-----------------------------------------------------------------------------------------------!
      XC(1)=0.5D0*(XP(1  ,J  )+XP(NCP1 ,J  ))
      XC(2)=       XP(NC1,J  )
      XC(3)=       XP(NC1,J+1)
      XC(4)=0.5D0*(XP(1  ,J+1)+XP(NCP1 ,J+1))
      YC(1)=0.5D0*(YP(1  ,J  )+YP(NCP1 ,J  ))
      YC(2)=       YP(NC1,J  )
      YC(3)=       YP(NC1,J+1)
      YC(4)=0.5D0*(YP(1  ,J+1)+YP(NCP1 ,J+1))
      ZC(1)=0.5D0*(ZP(1  ,J  )+ZP(NCP1 ,J  ))
      ZC(2)=       ZP(NC1,J  )
      ZC(3)=       ZP(NC1,J+1)
      ZC(4)=0.5D0*(ZP(1  ,J+1)+ZP(NCP1 ,J+1))
      CALL PANEL(XC,YC,ZC,XC0,YC0,ZC0,A1X,A1Y,A1Z,A2X,A2Y,A2Z,UNX,UNY,UNZ,A0)
!-----------------------------------------------------------------------------------------------!
!    Section Force                                                                              !
!-----------------------------------------------------------------------------------------------!
      FX=0.D0
      FY=0.D0
      FZ=0.D0
      MX=0.D0
      MY=0.D0
      MZ=0.D0
      FS=0.D0
      MS=0.D0
      FN=0.D0
      FT=0.D0
      DO I=1,NCP
!-----------------------------------------------------------------------------------------------!
!    Total Force                                                                                !
!-----------------------------------------------------------------------------------------------!
         FX=FX+CPNP(I,J,TT)*UNXP0(I,J)*AP0(I,J)
         FY=FY+CPNP(I,J,TT)*UNYP0(I,J)*AP0(I,J)
         FZ=FZ+CPNP(I,J,TT)*UNZP0(I,J)*AP0(I,J)
         MX=MX+CPNP(I,J,TT)*(UNYP0(I,J)*ZP0(I,J)-UNZP0(I,J)*YP0(I,J))*AP0(I,J)
         MY=MY+CPNP(I,J,TT)*(UNZP0(I,J)*XP0(I,J)-UNXP0(I,J)*ZP0(I,J))*AP0(I,J)
         MZ=MZ+CPNP(I,J,TT)*(UNXP0(I,J)*YP0(I,J)-UNYP0(I,J)*XP0(I,J))*AP0(I,J)
!-----------------------------------------------------------------------------------------------!
!    Normal Force                                                                               !
!-----------------------------------------------------------------------------------------------!
         FA=CPNP(I,J,TT)*UNXP0(I,J)*AP0(I,J)*UNX+ &
            CPNP(I,J,TT)*UNYP0(I,J)*AP0(I,J)*UNY+ &
            CPNP(I,J,TT)*UNZP0(I,J)*AP0(I,J)*UNZ
         FN=FN+FA
!-----------------------------------------------------------------------------------------------!
!    Suction Force Correction                                                                   !
!-----------------------------------------------------------------------------------------------!
         FS=FS+FA*UNX
         MS=MS+FA*UNY*ZP0(I,J)-FA*UNZ*YP0(I,J)
!-----------------------------------------------------------------------------------------------!
!    Tangential Force                                                                           !
!-----------------------------------------------------------------------------------------------!
         IF (ISF == 0) THEN
            FA=CPNP(I,J,TT)*UNXP0(I,J)*AP0(I,J)*(1.D0-UNX)+ &
               CPNP(I,J,TT)*UNYP0(I,J)*AP0(I,J)*(1.D0-UNY)+ &
               CPNP(I,J,TT)*UNZP0(I,J)*AP0(I,J)*(1.D0-UNZ)
            FT=FT+FA
         END IF !(ISF == 0)
      END DO !I=1,NCP
!-----------------------------------------------------------------------------------------------!
!    Inviscid Forces                                                                            !
!-----------------------------------------------------------------------------------------------!
      IF (IROTOR == 0) THEN
         IF (ISF == 0) THEN
            FIX=FIX+FX
            MIX=MIX+MX
         ELSEIF (ISF == 1) THEN
            FIX=FIX+FS
            MIX=MIX+MS
         END IF !(ISF)
      ELSEIF (IROTOR == 1) THEN
         IF (ISF == 0) THEN
            FIX=FIX-FX
            MIX=MIX-MX
         ELSEIF (ISF == 1) THEN
            FIX=FIX-FS
            MIX=MIX-MS
         END IF !(ISF)
         FN=-FN
      END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
!    Local Radius                                                                               !
!-----------------------------------------------------------------------------------------------!
      R =DSQRT(YC0*YC0+ZC0*ZC0)
      R1=DSQRT(YP(1,J  )*YP(1,J  )+ZP(1,J  )*ZP(1,J  ))
      R2=DSQRT(YP(1,J+1)*YP(1,J+1)+ZP(1,J+1)*ZP(1,J+1))
      DR=R2-R1
      R =0.5D0*(R1+R2)
!-----------------------------------------------------------------------------------------------!
!    Hydrodynamic Pitch                                                                         !
!-----------------------------------------------------------------------------------------------!
      IF (IROTOR == 0) THEN
         BETA =DATAND(UU(JJ)/PI/R)
         BETAI=DATAND(MX/FX/R)
      ELSEIF (IROTOR == 1) THEN
         BETA =DATAND(1.D0/UU(JJ)/R)
         BETAI=DATAND((-MX)/(-FX)/R)
      END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
!    Section Pitch                                                                              !
!-----------------------------------------------------------------------------------------------!
      PITCH1=DATAND(DABS(XP(1,J  )-XP(NC1,J  ))/R/DABS(TP(1,J  )-TP(NC1,J  )))
      PITCH2=DATAND(DABS(XP(1,J+1)-XP(NC1,J+1))/R/DABS(TP(1,J+1)-TP(NC1,J+1)))
      PITCH =PITCH1+(PITCH2-PITCH1)/(R2-R1)*(R-R1)
!-----------------------------------------------------------------------------------------------!
!    Local Angle of Attack                                                                      !
!-----------------------------------------------------------------------------------------------!
      IF (IROTOR == 0) ALPHA=PITCH-BETAI
      IF (IROTOR == 1) ALPHA=BETAI-PITCH
!-----------------------------------------------------------------------------------------------!
!    Define Leading Edge Coordinates                                                            !
!-----------------------------------------------------------------------------------------------!
      XLE =XP(NC1,J)+(XP(NC1,J+1)-XP(NC1,J))/(R2-R1)*(R-R1)
      TLE =TP(NC1,J)+(TP(NC1,J+1)-TP(NC1,J))/(R2-R1)*(R-R1)
      YBLE=R*TLE
!-----------------------------------------------------------------------------------------------!
!    Define Trailing Edge Coordinates                                                           !
!-----------------------------------------------------------------------------------------------!
      XTE =XP(1,J)+(XP(1,J+1)-XP(1,J))/(R2-R1)*(R-R1)
      TTE =TP(1,J)+(TP(1,J+1)-TP(1,J))/(R2-R1)*(R-R1)
      YBTE=R*TTE
!-----------------------------------------------------------------------------------------------!
!    Section Chord                                                                              !
!-----------------------------------------------------------------------------------------------!
      C=DSQRT((XTE-XLE)**2+(YBTE-YBLE)**2)
!-----------------------------------------------------------------------------------------------!
!    Inviscid Force                                                                             !
!-----------------------------------------------------------------------------------------------!
      IF (IROTOR == 0) THEN
         IF (ISF == 0) THEN
            FL= FX/DCOSD(BETAI)
         ELSEIF (ISF == 1) THEN
            FL= FS/DCOSD(BETAI)
         END IF !(ISF)
      END IF !(IROTOR == 0)
      IF (IROTOR == 1) THEN
         IF (ISF == 0) THEN
            FL=-FX/DCOSD(BETAI)
         ELSEIF (ISF == 1) THEN
            FL=-FS/DCOSD(BETAI)
         END IF !(ISF)
      END IF !(IROTOR == 1)
!-----------------------------------------------------------------------------------------------!
!    Section Effective Velocity                                                                 !
!-----------------------------------------------------------------------------------------------!
      IF ((J >= JI).AND.(J <= JF)) THEN
         CALL LININT(NRW,RR,POTPW(1,:,TT),1,R,CIRC)
         VE=FL/DR/PI/PI/DABS(CIRC)/2.D0
      ELSE !((J >= JI).AND.(J <= JF))
         IF (IROTOR == 0) VE=(UU(JJ)/PI)**2+ZC0**2+YC0**2
         IF (IROTOR == 1) VE=1.D0/UU(JJ)**2+ZC0**2+YC0**2
      END IF !((J >= JI).AND.(J <= JF))
!-----------------------------------------------------------------------------------------------!
!    Section Lift and Drag Coefficients                                                         !
!-----------------------------------------------------------------------------------------------!
      RES=RE(JJ)*C*VE
      FILEINP='CL_CD_INVL.DAT'
      CALL CLCD(R,RES,ALPHA,CLI,CDI,FILEINP)
      FILEINP='CL_CD_VISC.DAT'
      CALL CLCD(R,RES,ALPHA,CLV,CDV,FILEINP)
!-----------------------------------------------------------------------------------------------!
!    Inviscid                                                                                   !
!-----------------------------------------------------------------------------------------------!
      CFNI=FN/PI/PI/VE/VE/C/DR
      CFTI=FT/PI/PI/VE/VE/C/DR
!-----------------------------------------------------------------------------------------------!
!    Drag Correction                                                                            !
!-----------------------------------------------------------------------------------------------!
      CLF=FL/PI/PI/VE/VE/C/DR
      IF (IROTOR == 0) THEN
         FDX =FDX+PI*PI*(CLF*DCOSD(BETAI)-CDV*DSIND(BETAI))*VE*VE*C*DR
         MDX =MDX+PI*PI*(CLF*DSIND(BETAI)+CDV*DCOSD(BETAI))*VE*VE*C*DR*R
         CFND=CFNI+CDV*DSIND(ALPHA)
         CFTD=CFTI-CDV*DCOSD(ALPHA)
      ELSEIF (IROTOR == 1) THEN
         FDX =FDX+PI*PI*(CLF*DCOSD(BETAI)+CDV*DSIND(BETAI))*VE*VE*C*DR
         MDX =MDX+PI*PI*(CLF*DSIND(BETAI)-CDV*DCOSD(BETAI))*VE*VE*C*DR*R
         CFND=CFNI+CDV*DSIND(ALPHA)
         CFTD=CFTI-CDV*DCOSD(ALPHA)
      END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
!    Viscous Correction                                                                         !
!-----------------------------------------------------------------------------------------------!
      IF (CLI /= 0.D0) THEN
         FL  =FL*DABS(CLV/CLI) !Viscous Force
         CFNV=CFNI*DABS(CLV/CLI)*DCOSD(ALPHA)
         CFTV=CFTI*DABS(CLV/CLI)*DSIND(ALPHA)
      ELSE
         CFNV=CFNI
         CFTV=CFTI
      END IF !(CLI /= 0.D0)
      CLF=FL/PI/PI/VE/VE/C/DR
      IF (IROTOR == 0) THEN
         FVX =FVX+PI*PI*(CLF*DCOSD(BETAI)-CDV*DSIND(BETAI))*VE*VE*C*DR
         MVX =MVX+PI*PI*(CLF*DSIND(BETAI)+CDV*DCOSD(BETAI))*VE*VE*C*DR*R
         CFNV=CFNV+CDV*DSIND(ALPHA)
         CFTV=CFTV-CDV*DCOSD(ALPHA)
      ELSEIF (IROTOR == 1) THEN
         FVX =FVX+PI*PI*(CLF*DCOSD(BETAI)+CDV*DSIND(BETAI))*VE*VE*C*DR
         MVX =MVX+PI*PI*(CLF*DSIND(BETAI)-CDV*DCOSD(BETAI))*VE*VE*C*DR*R
         CFNV=CFNV+CDV*DSIND(ALPHA)
         CFTV=CFTV-CDV*DCOSD(ALPHA)
      END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
!    Write Results                                                                              !
!-----------------------------------------------------------------------------------------------!
      WRITE(33,210) J,R,C,PITCH,BETA,BETAI,ALPHA,RES,CLI,CDI,CLV,CDV, &
                    CFNI,CFTI,CFND,CFTD,CFNV,CFTV
!-----------------------------------------------------------------------------------------------!
   END DO !J=1,JF
!-----------------------------------------------------------------------------------------------!
!    Total Forces                                                                               !
!-----------------------------------------------------------------------------------------------!
   IF (IROTOR == 0) THEN
      CFIP(TT)=DFLOAT(NB)*FIX/ 8.D0
      CMIP(TT)=DFLOAT(NB)*MIX/16.D0
      CFDP(TT)=DFLOAT(NB)*FDX/ 8.D0
      CMDP(TT)=DFLOAT(NB)*MDX/16.D0
      CFVP(TT)=DFLOAT(NB)*FVX/ 8.D0
      CMVP(TT)=DFLOAT(NB)*MVX/16.D0
   ELSEIF (IROTOR == 1) THEN
      CFIP(TT)=DFLOAT(NB)*FIX*UU(JJ)**2/PI**3
      CMIP(TT)=DFLOAT(NB)*MIX*UU(JJ)**3/PI**3
      CFDP(TT)=DFLOAT(NB)*FDX*UU(JJ)**2/PI**3
      CMDP(TT)=DFLOAT(NB)*MDX*UU(JJ)**3/PI**3
      CFVP(TT)=DFLOAT(NB)*FVX*UU(JJ)**2/PI**3
      CMVP(TT)=DFLOAT(NB)*MVX*UU(JJ)**3/PI**3
   END IF !(IROTOR)
   WRITE(30,207) TT,CFIP(TT),CMIP(TT),CFDP(TT),CMDP(TT),CFVP(TT),CMVP(TT)
END DO !TT=0,NT
CLOSE(UNIT=30)
CLOSE(UNIT=33)
!-----------------------------------------------------------------------------------------------!
!    Blade Force Harmonic Analysis                                                              !
!-----------------------------------------------------------------------------------------------!
ALLOCATE(AMN1(0:NNHAR),BMN1(0:NNHAR),AMPLI1(0:NNHAR),PHASE1(0:NNHAR))
ALLOCATE(AMN2(0:NNHAR),BMN2(0:NNHAR),AMPLI2(0:NNHAR),PHASE2(0:NNHAR))
ALLOCATE(AMN3(0:NNHAR),BMN3(0:NNHAR),AMPLI3(0:NNHAR),PHASE3(0:NNHAR))
ALLOCATE(AMN4(0:NNHAR),BMN4(0:NNHAR),AMPLI4(0:NNHAR),PHASE4(0:NNHAR))
ALLOCATE(AMN5(0:NNHAR),BMN5(0:NNHAR),AMPLI5(0:NNHAR),PHASE5(0:NNHAR))
ALLOCATE(AMN6(0:NNHAR),BMN6(0:NNHAR),AMPLI6(0:NNHAR),PHASE6(0:NNHAR))
OUTFILE='FORCESHARM_'//IDENTU(JJ)//'.DAT'
OPEN(UNIT=34,FILE=OUTFILE,STATUS='UNKNOWN')
WRITE(34,'(A)') 'TITLE="Rotor Force Harmonic Analysis"'
WRITE(34,'(A)') 'VARIABLES="NNHAR" "CFIP" "CMIP" "CFDP" "CMDP" "CFVP" "CMVP"'
DO K=1,NREV
   WRITE(34,'(A,I4,A,I4)') 'ZONE T="BLADE REV=',K,'" F=POINT, I=',(NNHAR+1)
   KK=NTETA*(K-1)
!-----------------------------------------------------------------------------------------------!
!    Harmonic Analysis for CFIP                                                                 !
!-----------------------------------------------------------------------------------------------!
   DO I=0,NNHAR
      AN1=0.D0
      BN1=0.D0
      AN2=0.D0
      BN2=0.D0
      AN3=0.D0
      BN3=0.D0
      AN4=0.D0
      BN4=0.D0
      AN5=0.D0
      BN5=0.D0
      AN6=0.D0
      BN6=0.D0
      DO TT=KK+1,KK+NTETA
         AN1=AN1+CFIP(TT)*DCOS(DFLOAT(I)*DFLOAT(TT)*2.D0*PI/DFLOAT(NTETA))
         BN1=BN1+CFIP(TT)*DSIN(DFLOAT(I)*DFLOAT(TT)*2.D0*PI/DFLOAT(NTETA))
         AN2=AN2+CMIP(TT)*DCOS(DFLOAT(I)*DFLOAT(TT)*2.D0*PI/DFLOAT(NTETA))
         BN2=BN2+CMIP(TT)*DSIN(DFLOAT(I)*DFLOAT(TT)*2.D0*PI/DFLOAT(NTETA))
         AN3=AN3+CFDP(TT)*DCOS(DFLOAT(I)*DFLOAT(TT)*2.D0*PI/DFLOAT(NTETA))
         BN3=BN3+CFDP(TT)*DSIN(DFLOAT(I)*DFLOAT(TT)*2.D0*PI/DFLOAT(NTETA))
         AN4=AN4+CMDP(TT)*DCOS(DFLOAT(I)*DFLOAT(TT)*2.D0*PI/DFLOAT(NTETA))
         BN4=BN4+CMDP(TT)*DSIN(DFLOAT(I)*DFLOAT(TT)*2.D0*PI/DFLOAT(NTETA))
         AN5=AN5+CFVP(TT)*DCOS(DFLOAT(I)*DFLOAT(TT)*2.D0*PI/DFLOAT(NTETA))
         BN5=BN5+CFVP(TT)*DSIN(DFLOAT(I)*DFLOAT(TT)*2.D0*PI/DFLOAT(NTETA))
         AN6=AN6+CMVP(TT)*DCOS(DFLOAT(I)*DFLOAT(TT)*2.D0*PI/DFLOAT(NTETA))
         BN6=BN6+CMVP(TT)*DSIN(DFLOAT(I)*DFLOAT(TT)*2.D0*PI/DFLOAT(NTETA))
      END DO !TT=KK+1,KK+NTETA
      AMN1(I)=AN1*2.D0/DFLOAT(NTETA)
      BMN1(I)=BN1*2.D0/DFLOAT(NTETA)
      AMPLI1(I)=DSQRT(AMN1(I)**2+BMN1(I)**2)
      PHASE1(I)=DATAN2D(BMN1(I),AMN1(I))
      AMN2(I)=AN2*2.D0/DFLOAT(NTETA)
      BMN2(I)=BN2*2.D0/DFLOAT(NTETA)
      AMPLI2(I)=DSQRT(AMN2(I)**2+BMN2(I)**2)
      PHASE2(I)=DATAN2D(BMN2(I),AMN2(I))
      AMN3(I)=AN3*2.D0/DFLOAT(NTETA)
      BMN3(I)=BN3*2.D0/DFLOAT(NTETA)
      AMPLI3(I)=DSQRT(AMN3(I)**2+BMN3(I)**2)
      PHASE3(I)=DATAN2D(BMN3(I),AMN3(I))
      AMN4(I)=AN4*2.D0/DFLOAT(NTETA)
      BMN4(I)=BN4*2.D0/DFLOAT(NTETA)
      AMPLI4(I)=DSQRT(AMN4(I)**2+BMN4(I)**2)
      PHASE4(I)=DATAN2D(BMN4(I),AMN4(I))
      AMN5(I)=AN5*2.D0/DFLOAT(NTETA)
      BMN5(I)=BN5*2.D0/DFLOAT(NTETA)
      AMPLI5(I)=DSQRT(AMN5(I)**2+BMN5(I)**2)
      PHASE5(I)=DATAN2D(BMN5(I),AMN5(I))
      AMN6(I)=AN6*2.D0/DFLOAT(NTETA)
      BMN6(I)=BN6*2.D0/DFLOAT(NTETA)
      AMPLI6(I)=DSQRT(AMN6(I)**2+BMN6(I)**2)
      PHASE6(I)=DATAN2D(BMN6(I),AMN6(I))
!-----------------------------------------------------------------------------------------------!
      IF (I == 0) THEN
         WRITE(34,215) I,AMPLI1(I)/2.D0,PHASE1(I),AMPLI2(I)/2.D0,PHASE2(I), &
                         AMPLI3(I)/2.D0,PHASE3(I),AMPLI4(I)/2.D0,PHASE4(I), &
                         AMPLI5(I)/2.D0,PHASE5(I),AMPLI6(I)/2.D0,PHASE6(I)
      ELSE !(I == 0)
         WRITE(34,215) I,AMPLI1(I),PHASE1(I),AMPLI2(I),PHASE2(I),AMPLI3(I),PHASE3(I), &
                         AMPLI4(I),PHASE4(I),AMPLI5(I),PHASE5(I),AMPLI6(I),PHASE6(I)
      END IF !(I == 0)
   END DO !I=0,NNHAR
END DO !K=1,NREV
CLOSE(UNIT=34)
!-----------------------------------------------------------------------------------------------!
!    Formats                                                                                    !
!-----------------------------------------------------------------------------------------------!
207 FORMAT(I3,2X, 6(2X,E23.16))
210 FORMAT(I4,2X,17(2X,E23.16))
215 FORMAT(I3,2X,12(2X,E23.16))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE VISCOR
!-----------------------------------------------------------------------------------------------!
