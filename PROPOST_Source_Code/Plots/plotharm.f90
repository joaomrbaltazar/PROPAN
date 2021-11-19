!-----------------------------------------------------------------------------------------------!
!    Plot harmonic analysis                                                                     !
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
SUBROUTINE PLOTHARM(JJ)
!-----------------------------------------------------------------------------------------------!
!    Created by: 20032019, J. Baltazar, 2019 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPOST_MOD
IMPLICIT NONE
CHARACTER OUTFILE*50,ICHAR*3
INTEGER :: I,J,K,JJ,TT
DOUBLE PRECISION :: AN,BN
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: H0,H1,H2,H3
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: AMN,BMN,AMPLI,PHASE,PHASE1
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: H21,H22,H23
!-----------------------------------------------------------------------------------------------!
!    Pressure Harmonic Analysis                                                                 !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   OUTFILE='CPHARM_'//IDENTU(JJ)//'.DAT'
   OPEN(UNIT=32,FILE=OUTFILE,STATUS='UNKNOWN')
END IF !(IN == 1)
DO TT=0,NT
   IF (IN == 1) THEN
      WRITE(ICHAR,'(I3)') TT
      WRITE(32,'(A)') 'TITLE = "PROPAN_'//IDENTU(JJ)//'_T'//TRIM(ADJUSTL(ICHAR))//'"'
      WRITE(32,'(A)') 'VARIABLES="X","S/C","AN","BN","AMPLI","PHASE","PHASE1"'
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
!-----------------------------------------------------------------------------------------------!
      ALLOCATE(H0(NNX),H1(NNX),H2(NNX),H3(NNX))
      ALLOCATE(H21(NNXT1,NNTP),H22(NNXT1,NNTP),H23(NNXT1,NNTP))
      H0 =0.D0
      H1 =0.D0
      H2 =0.D0
      H3 =0.D0
      H21=0.D0
      H22=0.D0
      H23=0.D0
!-----------------------------------------------------------------------------------------------!
!    Interpolate to Constant Basic Axial Stations                                               !
!-----------------------------------------------------------------------------------------------!
      DO K=1,NNX !Inner side of nozzle
         WORK2=XN0(K,1)
         DO J=1,NNTP
            DO I=1,NNX
               H0(I)=XN0(I,J)
               H1(I)=DATAN2D(ZN0(I,J),YN0(I,J))
               H2(I)=DSQRT(ZN0(I,J)**2+YN0(I,J)**2)
               H3(I)=-0.5D0*CPNN(I,J,TT)
            END DO !I=1,NNX
!-----------------------------------------------------------------------------------------------!
            IF (INTER == 0) THEN
               CALL LININT(NNX,H0,H1,1,WORK2,H21(K,J))
               CALL LININT(NNX,H0,H2,1,WORK2,H22(K,J))
               CALL LININT(NNX,H0,H3,1,WORK2,H23(K,J))
            ELSEIF (INTER == 1) THEN
               CALL INTK1 (NNX,H0,H1,1,WORK2,H21(K,J))
               CALL INTK1 (NNX,H0,H2,1,WORK2,H22(K,J))
               CALL INTK1 (NNX,H0,H3,1,WORK2,H23(K,J))
            ELSEIF (INTER == 2) THEN
               CALL SPLINT(NNX,H0,H1,1,WORK2,H21(K,J))
               CALL SPLINT(NNX,H0,H2,1,WORK2,H22(K,J))
               CALL SPLINT(NNX,H0,H3,1,WORK2,H23(K,J))
            END IF !(INTER)
         END DO !J=1,NNTP
      END DO !K=1,NNX
!-----------------------------------------------------------------------------------------------!
      DO K=NNX+1,NNXT1 ! Outer of nozzle
         WORK2=XN0(K,1)
         DO J=1,NNTP
            DO I=1,NNX
               H0(I)=XN0(NNXT1-I+1,J)
               H1(I)=DATAN2D(ZN0(NNXT1-I+1,J),YN0(NNXT1-I+1,J))
               H2(I)=DSQRT(ZN0(NNXT1-I+1,J)**2+YN0(NNXT1-I+1,J)**2)
               H3(I)=-0.5D0*CPNN(NNXT1-I+1,J,TT)
            END DO !I=1,NNX
!-----------------------------------------------------------------------------------------------!
            IF (INTER == 0) THEN
               CALL LININT(NNX,H0,H1,1,WORK2,H21(K,J))
               CALL LININT(NNX,H0,H2,1,WORK2,H22(K,J))
               CALL LININT(NNX,H0,H3,1,WORK2,H23(K,J))
            ELSEIF (INTER == 1) THEN
               CALL INTK1 (NNX,H0,H1,1,WORK2,H21(K,J))
               CALL INTK1 (NNX,H0,H2,1,WORK2,H22(K,J))
               CALL INTK1 (NNX,H0,H3,1,WORK2,H23(K,J))
            ELSEIF (INTER == 2) THEN
               CALL SPLINT(NNX,H0,H1,1,WORK2,H21(K,J))
               CALL SPLINT(NNX,H0,H2,1,WORK2,H22(K,J))
               CALL SPLINT(NNX,H0,H3,1,WORK2,H23(K,J))
            END IF !(INTER)
         END DO !J=1,NNTP
      END DO !K=NNX+1,NNXT1
      DEALLOCATE(H0,H1,H2,H3,H22)
!-----------------------------------------------------------------------------------------------!
      ALLOCATE(H0(NNTP+1),H1(NNTP+1),H22(NNXT1,NNTP))
      ALLOCATE(AMN(NNXT1,0:NNHAR),BMN(NNXT1,0:NNHAR),AMPLI(NNXT1,0:NNHAR),PHASE(NNXT1,0:NNHAR))
      ALLOCATE(PHASE1(NNXT1,0:NNHAR))
      H0    =0.D0
      H1    =0.D0
      H22   =0.D0
      AMN   =0.D0
      BMN   =0.D0
      AMPLI =0.D0
      PHASE =0.D0
      PHASE1=0.D0
!-----------------------------------------------------------------------------------------------!
!    Interpolate to Equidistant Theta                                                           !
!-----------------------------------------------------------------------------------------------!
      DO K=1,NNXT1
         DO J=1,NNTP
            WORK2=H21(K,1)+DFLOAT(J-1)*360.D0/(DFLOAT(NB)*DFLOAT(NNTP))
            DO I=1,NNTP    
               H0(I)=H21(K,I)
               H1(I)=H23(K,I)
            END DO !I=1,NNTP
            H0(NNTP+1)=H0(1)+360.D0/(DFLOAT(NB))
            H1(NNTP+1)=H23(K,1)
!-----------------------------------------------------------------------------------------------!
            IF (INTER == 0) THEN
               CALL LININT(NNTP+1,H0,H1,1,WORK2,H22(K,J))
            ELSEIF (INTER == 1) THEN
               CALL INTK1 (NNTP+1,H0,H1,1,WORK2,H22(K,J))
            ELSEIF (INTER == 2) THEN
               CALL SPLINT(NNTP+1,H0,H1,1,WORK2,H22(K,J))
            END IF !(INTER)
         END DO !J=1,NNTP
!-----------------------------------------------------------------------------------------------!
!    Harmonic Analysis for Nozzle                                                               !
!-----------------------------------------------------------------------------------------------!
         DO I=0,NNHAR
            AN=0.D0
            BN=0.D0
            DO J=1,NNTP
               WORK2=-PI/DFLOAT(NB)+DFLOAT(J-1)*2.D0*PI/(DFLOAT(NB)*DFLOAT(NNTP))
               AN=AN+H22(K,J)*DCOS(DFLOAT(I)*DFLOAT(NB)*WORK2)
               BN=BN+H22(K,J)*DSIN(DFLOAT(I)*DFLOAT(NB)*WORK2)
            END DO !J=1,NNTP
            AMN(K,I)=AN*2.D0/DFLOAT(NNTP)
            BMN(K,I)=BN*2.D0/DFLOAT(NNTP)
!-----------------------------------------------------------------------------------------------!
!    Amplitude and Phase                                                                        !
!-----------------------------------------------------------------------------------------------!
            AMPLI(K,I)=DSQRT(AMN(K,I)**2+BMN(K,I)**2)
            PHASE(K,I)=DATAN2D(BMN(K,I),AMN(K,I))
!-----------------------------------------------------------------------------------------------!
!    Phase Relative to Theta=0                                                                  !
!-----------------------------------------------------------------------------------------------!
            IF (I >= 1) PHASE1(K,I)=PHASE(K,I)/DFLOAT(NB)+(H21(K,1)+180.D0/DFLOAT(NB))
         END DO !I=0,NNHAR
         AMPLI(K,0)=AMPLI(K,0)/2.D0
      END DO !K=1,NNXT1
      DEALLOCATE(H0,H1,H21,H22,H23)
!-----------------------------------------------------------------------------------------------!
      DO I=0,NNHAR
         WRITE(32,'(A,I3,A,I4,A)') ' ZONE T="NOZZLE N=',I,'", K=',NNXT1,' F=POINT'
         DO K=NNX,1,-1
            WRITE(32,205) XN0(K,1),SCN(K,1),AMN(K,I),BMN(K,I),AMPLI(K,I),PHASE(K,I),PHASE1(K,I)
         END DO !K=NXX,1,-1
         DO K=NNXT1,NNX+1,-1
            WRITE(32,205) XN0(K,1),SCN(K,1),AMN(K,I),BMN(K,I),AMPLI(K,I),PHASE(K,I),PHASE1(K,I)
         END DO !K=NNXT1,NNX+1,-1
      END DO !I=0,NNHAR
      DEALLOCATE(AMN,BMN,AMPLI,PHASE,PHASE1)
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Force Harmonic Analysis                                                              !
!-----------------------------------------------------------------------------------------------!
END DO !TT=0,NT
IF (IN == 1) CLOSE(UNIT=32)
!-----------------------------------------------------------------------------------------------!
!    Formats                                                                                    !
!-----------------------------------------------------------------------------------------------!
205 FORMAT(7(2X,E23.16))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PLOTHARM
!-----------------------------------------------------------------------------------------------!
