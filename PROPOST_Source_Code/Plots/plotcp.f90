!-----------------------------------------------------------------------------------------------!
!    Plot pressure distribution                                                                 !
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
SUBROUTINE PLOTCP(JJ)
!-----------------------------------------------------------------------------------------------!
!    Created by: 21052014, J. Baltazar, 2014 version 1.0                                        !
!    Modified  : 22052014, J. Baltazar, 2014 version 1.0                                        !
!    Modified  : 02072014, J. Baltazar, 2014 version 1.2                                        !
!    Modified  : 19012015, J. Baltazar, 2015 version 1.0                                        !
!    Modified  : 10032015, J. Baltazar, 2015 version 1.1 Unsteady Flow                          !
!    Modified  : 23112015, J. Baltazar, 2015 version 1.2                                        !
!    Modified  : 07072017, J. Baltazar, 2017 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPOST_MOD
IMPLICIT NONE
CHARACTER OUTFILE*50,ICHAR*3
INTEGER :: I,J,K,L,KK,JJ,TT
DOUBLE PRECISION :: XTE,XLE,YTE,YLE,ZTE,ZLE,TTE,TLE,YBTE,YBLE,YB,PHI,SL,CL
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: R0,T0,R1,T1
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: H0,H1,H2,H3,H4
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: RP0,TP0
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: RN0,TN0
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: H21,H22,H23,H24
!-----------------------------------------------------------------------------------------------!
!    Initialise Variables                                                                       !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) ALLOCATE(SCP(NCP,NRP))
IF (IN == 1) ALLOCATE(SCN(NNXT1,NNTP*NB))
!-----------------------------------------------------------------------------------------------!
!    Write Radial Distribution in TecPlot Format                                                !
!-----------------------------------------------------------------------------------------------!
OUTFILE='CP2DP_'//IDENTU(JJ)//'.DAT'
OPEN(UNIT=31,FILE=OUTFILE,STATUS='UNKNOWN')
DO TT=0,NT
   WRITE(ICHAR,'(I3)') TT
   WRITE(31,'(A)') 'TITLE = "PROPAN_'//IDENTU(JJ)//'_T'//TRIM(ADJUSTL(ICHAR))//'"'
!-----------------------------------------------------------------------------------------------!
!    Definition of Chordwise and Radial Coordinates for Blade                                   !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
!-----------------------------------------------------------------------------------------------!
      ALLOCATE(R0(NRP1),T0(NRP1),H0(NRP1),R1(NRP1),T1(NRP1),H1(NRP1))
      ALLOCATE(RP0(NCP,NRP),TP0(NCP,NRP))
      R0 =0.D0
      T0 =0.D0
      H0 =0.D0
      R1 =0.D0
      T1 =0.D0
      H1 =0.D0
      SCP=0.D0
      RP0=0.D0
      TP0=0.D0
!-----------------------------------------------------------------------------------------------!
!    Blade Panel Control Points in Cylindrical Coordinates                                      !
!-----------------------------------------------------------------------------------------------!
      RP0=DSQRT(YP0*YP0+ZP0*ZP0)
      TP0=DATAN2(ZP0,YP0)
!-----------------------------------------------------------------------------------------------!
      I=1
      R0(:)=DSQRT(YP(I,:)*YP(I,:)+ZP(I,:)*ZP(I,:))
      T0(:)=DATAN2(ZP(I,:),YP(I,:))
      H0(:)=XP(I,:)
!-----------------------------------------------------------------------------------------------!
      I=NC1
      R1(:)=DSQRT(YP(I,:)*YP(I,:)+ZP(I,:)*ZP(I,:))
      T1(:)=DATAN2(ZP(I,:),YP(I,:))
      H1(:)=XP(I,:)
!-----------------------------------------------------------------------------------------------!
      DO J=1,NRP
         DO I=1,NCP
!-----------------------------------------------------------------------------------------------!
!    Interpolate for the Blade Trailing Edge                                                    !
!-----------------------------------------------------------------------------------------------!
            IF (ISTRIP /= 2) THEN
               IF (INTER == 0) THEN
                  CALL LININT(NRP1,R0,T0,1,RP0(I,J),TTE)
                  CALL LININT(NRP1,R0,H0,1,RP0(I,J),XTE)
               ELSEIF (INTER == 1) THEN
                  CALL INTK1 (NRP1,R0,T0,1,RP0(I,J),TTE)
                  CALL INTK1 (NRP1,R0,H0,1,RP0(I,J),XTE)
               ELSEIF (INTER == 2) THEN
                  CALL SPLINT(NRP1,R0,T0,1,RP0(I,J),TTE)
                  CALL SPLINT(NRP1,R0,H0,1,RP0(I,J),XTE)
               END IF !(INTER)
            ELSE !(ISTRIP /= 2)
               IF (INTER == 0) THEN
                  CALL LININT(NRP ,R0(1:NRP),T0(1:NRP),1,RP0(I,J),TTE)
                  CALL LININT(NRP ,R0(1:NRP),H0(1:NRP),1,RP0(I,J),XTE)
               ELSEIF (INTER == 1) THEN
                  CALL INTK1 (NRP ,R0(1:NRP),T0(1:NRP),1,RP0(I,J),TTE)
                  CALL INTK1 (NRP ,R0(1:NRP),H0(1:NRP),1,RP0(I,J),XTE)
               ELSEIF (INTER == 2) THEN
                  CALL SPLINT(NRP ,R0(1:NRP),T0(1:NRP),1,RP0(I,J),TTE)
                  CALL SPLINT(NRP ,R0(1:NRP),H0(1:NRP),1,RP0(I,J),XTE)
               END IF !(INTER)
            END IF !(ISTRIP /= 2)
!-----------------------------------------------------------------------------------------------!
            YTE=RP0(I,J)*DCOS(TTE)
            ZTE=RP0(I,J)*DSIN(TTE)
!-----------------------------------------------------------------------------------------------!
!    Interpolate for the Blade Leading Edge                                                     !
!-----------------------------------------------------------------------------------------------!
            IF (ISTRIP /= 2) THEN
               IF (INTER == 0) THEN
                  CALL LININT(NRP1,R1,T1,1,RP0(I,J),TLE)
                  CALL LININT(NRP1,R1,H1,1,RP0(I,J),XLE)
               ELSEIF (INTER == 1) THEN
                  CALL INTK1 (NRP1,R1,T1,1,RP0(I,J),TLE)
                  CALL INTK1 (NRP1,R1,H1,1,RP0(I,J),XLE)
               ELSEIF (INTER == 2) THEN
                  CALL SPLINT(NRP1,R1,T1,1,RP0(I,J),TLE)
                  CALL SPLINT(NRP1,R1,H1,1,RP0(I,J),XLE)
               END IF !(INTER)
            ELSE !(ISTRIP /= 2)
               IF (INTER == 0) THEN
                  CALL LININT(NRP ,R1(1:NRP),T1(1:NRP),1,RP0(I,J),TLE)
                  CALL LININT(NRP ,R1(1:NRP),H1(1:NRP),1,RP0(I,J),XLE)
               ELSEIF (INTER == 1) THEN
                  CALL INTK1 (NRP ,R1(1:NRP),T1(1:NRP),1,RP0(I,J),TLE)
                  CALL INTK1 (NRP ,R1(1:NRP),H1(1:NRP),1,RP0(I,J),XLE)
               ELSEIF (INTER == 2) THEN
                  CALL SPLINT(NRP ,R1(1:NRP),T1(1:NRP),1,RP0(I,J),TLE)
                  CALL SPLINT(NRP ,R1(1:NRP),H1(1:NRP),1,RP0(I,J),XLE)
               END IF !(ISTRIP /= 2)
            END IF !(ISTRIP /= 2)
!-----------------------------------------------------------------------------------------------!
            YLE=RP0(I,J)*DCOS(TLE)
            ZLE=RP0(I,J)*DSIN(TLE)
!-----------------------------------------------------------------------------------------------!
            YBLE=RP0(I,J)*TLE
            YBTE=RP0(I,J)*TTE
!-----------------------------------------------------------------------------------------------!
            PHI=DATAN2(XTE-XLE,YBTE-YBLE)
!-----------------------------------------------------------------------------------------------!
            YB=RP0(I,J)*TP0(I,J)
!-----------------------------------------------------------------------------------------------!
            SL=(YB-YBLE)*DCOS(PHI)+(XP0(I,J)-XLE)*DSIN(PHI)
            CL=DSQRT((XTE-XLE)**2+(YBTE-YBLE)**2)
!-----------------------------------------------------------------------------------------------!
            SCP(I,J)=SL/CL
         END DO !I=1,NCP
      END DO !J=1,NRP
      DEALLOCATE(R0,T0,H0,R1,T1,H1)
!-----------------------------------------------------------------------------------------------!
!    Write Radial Distribution of PHI and CP in TecPlot Format                                  !
!-----------------------------------------------------------------------------------------------!
      WRITE(31,'(A)') 'VARIABLES="S/C","r/R","POT","CP","CPN"'
!-----------------------------------------------------------------------------------------------!
!    2D Blade Pressure Distribution                                                             !
!-----------------------------------------------------------------------------------------------!
      ALLOCATE(H0(NRP),H1(NRP),H2(NRP),H3(NRP),H4(NRP))
      ALLOCATE(H21(NRS,NCP),H22(NRS,NCP),H23(NRS,NCP),H24(NRS,NCP))
      H0 =0.D0
      H1 =0.D0
      H2 =0.D0
      H3 =0.D0
      H4 =0.D0
      H21=0.D0
      H22=0.D0
      H23=0.D0
      H24=0.D0
!-----------------------------------------------------------------------------------------------!
      DO K=1,NRS 
         WRITE(31,100) 'ZONE T="BLADE r/R=',RS(K),'", I=',NCP,' F=POINT'
         DO I=1,NCP
            H0(:)=DSQRT(YP0(I,:)*YP0(I,:)+ZP0(I,:)*ZP0(I,:))
            H1(:)=POTP(I,:,TT)
            H2(:)= CPP(I,:,TT)
            H3(:)= SCP(I,:)
            H4(:)=CPNP(I,:,TT)
!-----------------------------------------------------------------------------------------------!
            IF (ISTRIP /= 2) THEN
               IF (INTER == 0) THEN
                  CALL LININT(NRP,H0,H1,1,RS(K),H21(K,I))
                  CALL LININT(NRP,H0,H2,1,RS(K),H22(K,I))
                  CALL LININT(NRP,H0,H3,1,RS(K),H23(K,I))
                  CALL LININT(NRP,H0,H4,1,RS(K),H24(K,I))
               ELSEIF (INTER == 1) THEN
                  CALL INTK1 (NRP,H0,H1,1,RS(K),H21(K,I))
                  CALL INTK1 (NRP,H0,H2,1,RS(K),H22(K,I))
                  CALL INTK1 (NRP,H0,H3,1,RS(K),H23(K,I))
                  CALL INTK1 (NRP,H0,H4,1,RS(K),H24(K,I))
               ELSEIF (INTER == 2) THEN
                  CALL SPLINT(NRP,H0,H1,1,RS(K),H21(K,I))
                  CALL SPLINT(NRP,H0,H2,1,RS(K),H22(K,I))
                  CALL SPLINT(NRP,H0,H3,1,RS(K),H23(K,I))
                  CALL SPLINT(NRP,H0,H4,1,RS(K),H24(K,I))
               END IF !(INTER)
            ELSE !(ISTRIP /= 2)
               IF (INTER == 0) THEN
                  CALL LININT(NRP-1,H0(1:NRP-1),H1(1:NRP-1),1,RS(K),H21(K,I))
                  CALL LININT(NRP-1,H0(1:NRP-1),H2(1:NRP-1),1,RS(K),H22(K,I))
                  CALL LININT(NRP-1,H0(1:NRP-1),H3(1:NRP-1),1,RS(K),H23(K,I))
                  CALL LININT(NRP-1,H0(1:NRP-1),H4(1:NRP-1),1,RS(K),H24(K,I))
               ELSEIF (INTER == 1) THEN
                  CALL INTK1 (NRP-1,H0(1:NRP-1),H1(1:NRP-1),1,RS(K),H21(K,I))
                  CALL INTK1 (NRP-1,H0(1:NRP-1),H2(1:NRP-1),1,RS(K),H22(K,I))
                  CALL INTK1 (NRP-1,H0(1:NRP-1),H3(1:NRP-1),1,RS(K),H23(K,I))
                  CALL INTK1 (NRP-1,H0(1:NRP-1),H4(1:NRP-1),1,RS(K),H24(K,I))
               ELSEIF (INTER == 2) THEN
                  CALL SPLINT(NRP-1,H0(1:NRP-1),H1(1:NRP-1),1,RS(K),H21(K,I))
                  CALL SPLINT(NRP-1,H0(1:NRP-1),H2(1:NRP-1),1,RS(K),H22(K,I))
                  CALL SPLINT(NRP-1,H0(1:NRP-1),H3(1:NRP-1),1,RS(K),H23(K,I))
                  CALL SPLINT(NRP-1,H0(1:NRP-1),H4(1:NRP-1),1,RS(K),H24(K,I))
               END IF !(INTER)
            END IF !(ISTRIP /= 2)
            WRITE(31,204) H23(K,I),RS(K),H21(K,I),H22(K,I),H24(K,I)
         END DO !I=1,NCP
      END DO !K=1,NRS
      DEALLOCATE(H0,H1,H2,H3,H4,H21,H22,H23,H24)
      DEALLOCATE(RP0,TP0)
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Definition of Chordwise and Radial Coordinates for Nozzle                                  !
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
!-----------------------------------------------------------------------------------------------!
      ALLOCATE(R0(NNTP*NB),T0(NNTP*NB),H0(NNTP*NB),R1(NNTP*NB),T1(NNTP*NB),H1(NNTP*NB))
      ALLOCATE(RN0(NNXT1,NNTP*NB),TN0(NNXT1,NNTP*NB))
      R0 =0.D0
      T0 =0.D0
      H0 =0.D0
      R1 =0.D0
      T1 =0.D0
      H1 =0.D0
      SCN=0.D0
      RN0=0.D0
      TN0=0.D0
!-----------------------------------------------------------------------------------------------!
!    Nozzle Panel Control Points in Cylindrical Coordinates                                     !
!-----------------------------------------------------------------------------------------------!
      DO K=1,NB
         PHI=DFLOAT(K-1)*360.D0/DFLOAT(NB)
         DO J=1,NNTP
            L=J+(K-1)*NNTP
            DO I=1,NNXT1
               RN0(I,L)=DSQRT(YN0(I,J)*YN0(I,J)+ZN0(I,J)*ZN0(I,J))
               TN0(I,L)=DATAN2D(ZN0(I,J),YN0(I,J))+PHI
            END DO !I=1,NNXT1
         END DO !J=1,NNTP
      END DO !K=1,NB
!-----------------------------------------------------------------------------------------------!
      I=1
      DO K=1,NB
         PHI=DFLOAT(K-1)*360.D0/DFLOAT(NB)
         DO J=1,NNTP
            L=J+(K-1)*NNTP
            R0(L)=DSQRT(YN(I,J)*YN(I,J)+ZN(I,J)*ZN(I,J))
            T0(L)=DATAN2D(ZN(I,J),YN(I,J))+PHI
            H0(L)=XN(I,J)
         END DO !J=1,NNTP
      END DO !K=1,NB
!-----------------------------------------------------------------------------------------------!
      I=NNX1
      DO K=1,NB
         PHI=DFLOAT(K-1)*360.D0/DFLOAT(NB)
         DO J=1,NNTP
            L=J+(K-1)*NNTP
            R1(L)=DSQRT(YN(I,J)*YN(I,J)+ZN(I,J)*ZN(I,J))
            T1(L)=DATAN2D(ZN(I,J),YN(I,J))+PHI
            H1(L)=XN(I,J)
         END DO !J=1,NNTP
      END DO !K=1,NB
!-----------------------------------------------------------------------------------------------!
      DO K=1,NB
         DO J=1,NNTP
            L=J+(K-1)*NNTP
            DO I=1,NNXT1
!-----------------------------------------------------------------------------------------------!
!    Interpolate for the Nozzle Leading Edge                                                    !
!-----------------------------------------------------------------------------------------------!
               WORK1=TN0(I,L)
               IF (TN0(I,L) < T0(1      )) WORK1=TN0(I,L)+360.D0
               IF (TN0(I,L) > T0(NNTP*NB)) WORK1=TN0(I,L)-360.D0
               IF (INTER == 0) THEN
                  CALL LININT(NNTP*NB,T0,H0,1,WORK1,XLE)
               ELSEIF (INTER == 1) THEN
                  CALL INTK1 (NNTP*NB,T0,H0,1,WORK1,XLE)
               ELSEIF (INTER == 2) THEN
                  CALL SPLINT(NNTP*NB,T0,H0,1,WORK1,XLE)
               END IF !(INTER)
!-----------------------------------------------------------------------------------------------!
!    Interpolate for the Nozzle Trailing Edge                                                   !
!-----------------------------------------------------------------------------------------------!
               WORK1=TN0(I,L)
               IF (TN0(I,L) < T1(1      )) WORK1=TN0(I,L)+360.D0
               IF (TN0(I,L) > T1(NNTP*NB)) WORK1=TN0(I,L)-360.D0
               IF (INTER == 0) THEN
                  CALL LININT(NNTP*NB,T1,H1,1,WORK1,XTE)
               ELSEIF (INTER == 1) THEN
                  CALL INTK1 (NNTP*NB,T1,H1,1,WORK1,XTE)
               ELSEIF (INTER == 2) THEN
                  CALL SPLINT(NNTP*NB,T1,H1,1,WORK1,XTE)
               END IF !(INTER)
!-----------------------------------------------------------------------------------------------!
               SL=XN0(I,J)-XLE
               CL=XTE-XLE
!-----------------------------------------------------------------------------------------------!
               SCN(I,L)=SL/CL
            END DO !I=1,NNXT1
         END DO !J=1,NNTP
      END DO !K=1,NB
      DEALLOCATE(R0,T0,H0,R1,T1,H1)
!-----------------------------------------------------------------------------------------------!
!    Write Radial Distribution of PHI and CP in TecPlot Format                                  !
!-----------------------------------------------------------------------------------------------!
      WRITE(31,'(A)') 'VARIABLES="S/C","THETA","POT","CP","CPN"'
!-----------------------------------------------------------------------------------------------!
!    2D Blade Pressure Distribution                                                             !
!-----------------------------------------------------------------------------------------------!
      ALLOCATE(H0(NNTP*NB),H1(NNTP*NB),H2(NNTP*NB),H3(NNTP*NB),H4(NNTP*NB))
      ALLOCATE(H21(NTS,NNXT1),H22(NTS,NNXT1),H23(NTS,NNXT1),H24(NTS,NNXT1))
      H0 =0.D0
      H1 =0.D0
      H2 =0.D0
      H3 =0.D0
      H4 =0.D0
      H21=0.D0
      H22=0.D0
      H23=0.D0
      H24=0.D0
!-----------------------------------------------------------------------------------------------!
      DO K=1,NTS
         DO I=1,NNXT1
            DO KK=1,NB
               PHI=DFLOAT(KK-1)*360.D0/DFLOAT(NB)
               DO J=1,NNTP
                  L=J+(KK-1)*NNTP
                  H0(L)=DATAN2D(ZN0(I,J),YN0(I,J))+PHI
                  H1(L)=POTN(I,J,TT)
                  WORK1=DSQRT(YN0(I,J)*YN0(I,J)+ZN0(I,J)*ZN0(I,J))
                  H2(L)= CPN(I,J,TT)
                  H3(L)= SCN(I,L)
               END DO !J=1,NNTP*NB
            END DO !KK=1,NB
            WORK2=TS(K)
            IF (TS(K) < H0(1      )) WORK2=TS(K)+360.D0
            IF (TS(K) > H0(NNTP*NB)) WORK2=TS(K)-360.D0
            IF (INTER == 0) THEN
               CALL LININT(NNTP*NB,H0,H1,1,WORK2,H21(K,I))
               CALL LININT(NNTP*NB,H0,H2,1,WORK2,H22(K,I))
               CALL LININT(NNTP*NB,H0,H3,1,WORK2,H23(K,I))
               CALL LININT(NNTP*NB,H0,H4,1,WORK2,H24(K,I))
            ELSEIF (INTER == 1) THEN
               CALL INTK1 (NNTP*NB,H0,H1,1,WORK2,H21(K,I))
               CALL INTK1 (NNTP*NB,H0,H2,1,WORK2,H22(K,I))
               CALL INTK1 (NNTP*NB,H0,H3,1,WORK2,H23(K,I))
               CALL INTK1 (NNTP*NB,H0,H4,1,WORK2,H24(K,I))
            ELSEIF (INTER == 2) THEN
               CALL SPLINT(NNTP*NB,H0,H1,1,WORK2,H21(K,I))
               CALL SPLINT(NNTP*NB,H0,H2,1,WORK2,H22(K,I))
               CALL SPLINT(NNTP*NB,H0,H3,1,WORK2,H23(K,I))
               CALL SPLINT(NNTP*NB,H0,H4,1,WORK2,H24(K,I))
            END IF !(INTER)
         END DO !I=1,NNXT1
      END DO !K=1,NTS
!-----------------------------------------------------------------------------------------------!
      DO K=1,NTS
         WRITE(31,100) 'ZONE T="NOZZLE theta=',TS(K),'", I=',NNXT1,' F=POINT'
         DO I=NNX,1,-1
            WRITE(31,204) H23(K,I),TS(K),H21(K,I),H22(K,I),H24(K,I)
         END DO !I=NXX,1,-1
         DO I=NNXT1,NNX+1,-1
            WRITE(31,204) H23(K,I),TS(K),H21(K,I),H22(K,I),H24(K,I)
         END DO !I=NNXT1,NNX+1,-1
      END DO !K=1,NTS
      DEALLOCATE(H0,H1,H2,H3,H4,H21,H22,H23,H24)
      DEALLOCATE(RN0,TN0)
   END IF !(IN == 1)
END DO !TT=0,NT
CLOSE(UNIT=31)
!-----------------------------------------------------------------------------------------------!
!    Formats                                                                                    !
!-----------------------------------------------------------------------------------------------!
100 FORMAT(A,F6.3,A,I4,A)
204 FORMAT(5(2X,E23.16))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PLOTCP
!-----------------------------------------------------------------------------------------------!
