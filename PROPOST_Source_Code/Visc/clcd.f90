!-----------------------------------------------------------------------------------------------!
!    Subroutine computes angle of attack and Cd values for given blade radius, Reynolds number  !
!    and Cl values from tabulated data in a file. CL(alpha) must be single valued monotonous    !
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
!                                                                                               !
!    Input: R:     Section Radius                                                               !
!           RE:    Section Reynolds number                                                      !
!           ALPHA: Section angle of attack                                                      !
!    Output:CL:    Section Lift coefficient                                                     !
!           CD:    Section Drag coefficient                                                     !
!           CL:    Section Lift coefficient in case it exceeds the maximum                      !
!                                                                                               !
!    Modified  : 07072017, J. Baltazar, 2017 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
SUBROUTINE CLCD(R,RE,ALPHA,CL,CD,FILEINP)
!-----------------------------------------------------------------------------------------------!
CHARACTER*40 FILEINP
INTEGER :: I,J,K,J1,K1,II,JJ,KK,NR,NRE(20),NALPHA_CL(20,20),NALPHA_CD(20,20),LFILEINP
DOUBLE PRECISION :: R,RE,CL,ALPHA,CD
DOUBLE PRECISION :: R_SECTION(20),RE_TABLE(20,20),ALPHA_CL_TABLE(200,20,20),ALPHA_CD_TABLE(200,20,20)
DOUBLE PRECISION :: CL_TABLE(200,20,20),CD_TABLE(200,20,20)
DOUBLE PRECISION :: X(200),Y(200),XX,CSI,ETA,ZETA,A(2),AA(2),CDA(2),CDD(2),CLA(2),CLL(2)
!-----------------------------------------------------------------------------------------------!
!    Read Tabulated Data                                                                        !
!-----------------------------------------------------------------------------------------------!
I=1 
DO WHILE (FILEINP(I:I).NE.' ')
   I=I+1
END DO
LFILEINP=I-1
!-----------------------------------------------------------------------------------------------!
OPEN(UNIT=15,FILE=FILEINP(1:LFILEINP),STATUS='UNKNOWN')
READ(15,*) NR                            ! Number of section radii
READ(15,*) (R_SECTION(K),K=1,NR)         ! Section radii
READ(15,*) (NRE(K),K=1,NR)               ! Number of Reynolds numbers per section radius
DO K=1,NR
   READ(15,*)
   READ(15,*) (RE_TABLE(J,K),J=1,NRE(K))  ! Reynolds numbers
   READ(15,*) (NALPHA_CL(J,K),J=1,NRE(K)) ! Number of angles of attack of CL for each Reynolds number
   READ(15,*) (NALPHA_CD(J,K),J=1,NRE(K)) ! Number of angles of attack of CD for each Reynolds number
!  READ(15,*) (CL_MAX_TABLE(J,K),J=1,NRE(K))    ! Value of maximum CL for each Reynolds number
   DO J=1,NRE(K)
      READ(15,*)
      READ(15,*)
      DO I=1,NALPHA_CL(J,K)
         READ(15,*) ALPHA_CL_TABLE(I,J,K),CL_TABLE(I,J,K)
      END DO     ! I=1,NALPHA_CL(J,K)
      READ(15,*)
      DO I=1,NALPHA_CD(J,K)
         READ(15,*) ALPHA_CD_TABLE(I,J,K),CD_TABLE(I,J,K)
      END DO !I=1,NALPHA_CD(J,K)
   END DO !J=1,NRE(K)
END DO !K=1,NR
CLOSE(15)
!-----------------------------------------------------------------------------------------------!
!    Interpolation of Table Data                                                                !
!    Linear in radius and Reynolds number. Linear in angle of attack                            !
!-----------------------------------------------------------------------------------------------!
DO K=1,NR
   X(K)=R_SECTION(K)
END DO !K=1,NR
CALL INTERVAL(NR,X,R,KK,ZETA)           ! Define the interval for the section radius
IF (ZETA < 0.D0) THEN
   ZETA=0.D0
!* WRITE(30,*) 'R<R(1)  - Minimum R Limit =',R_SECTION(1)
END IF !(ZETA < 0.D0)
IF (ZETA > 1.D0) THEN
   ZETA=1.D0
!* WRITE(30,*) 'R>R(NR) - Maximum R Limit =',R_SECTION(NR)
END IF !(ZETA > 1.D0)
!-----------------------------------------------------------------------------------------------!
DO K=1,2              ! Loop on the two radii for interpolation
   K1=KK+K-1             ! Index for entering the table in radius
   DO J=1,NRE(K1)
      X(J)=RE_TABLE(J,K1)
   END DO !J=1,NRE(K1)
   CALL INTERVAL(NRE(K1),X,RE,JJ,ETA)           ! Define the interval for the Reynolds number
   IF (ETA < 0.D0) THEN
      ETA=0.D0
!*    WRITE(30,*) 'RE<RE(1,K1)        - Minimum Re Limit =',RE_TABLE(1,K1)
   END IF !(ETA < 0.D0)
   IF (ETA > 1.D0) THEN
      ETA=1.D0
!*    WRITE(30,*) 'RE>RE(NRE(K1),K1)  - Maximum Re Limit =',RE_TABLE(NRE(K1),K1)
   END IF !(ETA > 1.D0)
!-----------------------------------------------------------------------------------------------!
   DO J=1,2
      J1=JJ+J-1
!-----------------------------------------------------------------------------------------------!
      DO I=1,NALPHA_CL(J1,K1)
         X(I)=ALPHA_CL_TABLE(I,J1,K1)
         Y(I)=CL_TABLE(I,J1,K1)
      END DO !I=1,NALPHA_CL(J1,K1)
!-----------------------------------------------------------------------------------------------!
      A(J)=ALPHA
      CALL INTERVAL(NALPHA_CL(J1,K1),X,ALPHA,II,CSI)           ! Define the interval for the CL
      IF (CSI < 0.D0) THEN
!        WRITE(*,*) 'CL<CL(1,J1,K1)        - LINEAR EXTRAPOLATION J1=', J1
         CSI=0.D0
      END IF !(CSI < 0.D0)
      IF (CSI > 1.D0) THEN                                    ! Define maximum CL requires monotonous Cl(alpha)
!        WRITE(*,*) 'CL>CL(NALPHA_CL(J1,K1),J1,K1)  - Maximum CL Limit J1=', J1, ' CL=', CL_TABLE(NALPHA_CL(J1,K1),J1,K1)
         CLA(J)=Y(NALPHA_CL(J1,K1))
         A(J)=X(NALPHA_CL(J1,K1))
         CSI=1.D0
      END IF !(CSI > 1.D0)
      CLA(J)=(1.D0-CSI)*Y(II)+CSI*Y(II+1)
!     CALL INTK1(NALPHA_CL(J1,K1),X,CL,Y,A(J))    ! Interpolated angle of attack for the cl value
!-----------------------------------------------------------------------------------------------!
      DO I=1,NALPHA_CD(J1,K1)
         X(I)=ALPHA_CD_TABLE(I,J1,K1)
         Y(I)=CD_TABLE(I,J1,K1)
      END DO !I=1,NALPHA_CD(J1,K1)
!-----------------------------------------------------------------------------------------------!
      CALL INTERVAL(NALPHA_CD(J1,K1),X,A(J),II,CSI)           ! Define the interval for the CD
      IF (CSI < 0.D0) THEN
!*       WRITE(*,*) 'CD<CD(1,J1,K1)        - LINEAR EXTRAPOLATION J1=', J1
         CSI=0.D0
      END IF !(CSI < 0.D0)
      IF (CSI > 1.D0) THEN
!*       WRITE(*,*) 'CD>CD(NALPHA_CD(J1,K1),J1,K1)  - LINEAR EXTRAPOLATION J1=', J1
         CSI=1.D0
      END IF !(CSI > 1.D0)
      CDA(J)=(1.D0-CSI)*Y(II)+CSI*Y(II+1)
!     CALL INTK1(NALPHA_CD(J1,K1),X,CD,Y,CDA(J))    ! interpolated Cd for the angle of attack
   END DO !J=1,2
!-----------------------------------------------------------------------------------------------!
   AA (K)=(1.D0-ETA)*A  (1)+ETA*A  (2)                   ! Interpolated angle of attack for the Reynolds number
   CDD(K)=(1.D0-ETA)*CDA(1)+ETA*CDA(2) 
   CLL(K)=(1.D0-ETA)*CLA(1)+ETA*CLA(2) 
END DO !K=1,2 
!-----------------------------------------------------------------------------------------------!
ALPHA=(1.D0-ZETA)*AA (1)+ZETA*AA (2)
CD   =(1.D0-ZETA)*CDD(1)+ZETA*CDD(2)
CL   =(1.D0-ZETA)*CLL(1)+ZETA*CLL(2)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE CLCD
!-----------------------------------------------------------------------------------------------!
