!-----------------------------------------------------------------------------------------------!
!    Plot wake geometries                                                                       !
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
SUBROUTINE PLOTWAKE(JJ)
!-----------------------------------------------------------------------------------------------!
!    Created by: 22052014, J. Baltazar, version 1.0                                             !
!    Modified  : 22052014, J. Baltazar, version 1.0                                             !
!    Modified  : 09032015, J. Baltazar, 2015 version 1.1 Unsteady Flow                          !
!-----------------------------------------------------------------------------------------------!
USE PROPOST_MOD
IMPLICIT NONE
CHARACTER OUTFILE*50,ICHAR*3
INTEGER :: I,J,K,JJ,TT,IMIN
DOUBLE PRECISION :: XPWM
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: H0,H1,H2,H3,H4
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: H21,H22,H23,H24
!-----------------------------------------------------------------------------------------------!
!    Write Wake Geometry in TecPlot Format                                                      !
!-----------------------------------------------------------------------------------------------!
OUTFILE='WAKE_'//IDENTU(JJ)//'.DAT'
OPEN(UNIT=32,FILE=OUTFILE,STATUS='UNKNOWN')
DO TT=0,NT
   WRITE(ICHAR,'(I3)') TT
   WRITE(32,'(A)') 'TITLE = "PROPAN_'//IDENTU(JJ)//'_T'//TRIM(ADJUSTL(ICHAR))//'"'
   WRITE(32,'(A)') ' VARIABLES="x/R","y/R","z/R","r/R","P/D","BETAI"'
!-----------------------------------------------------------------------------------------------!
!    Blade Wake Geometry                                                                        !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      ALLOCATE(H0(NPW1),H1(NPW1),H2(NPW1),H3(NPW1),H4(NPW1))
      ALLOCATE(H21(NXS,NRW1),H22(NXS,NRW1),H23(NXS,NRW1),H24(NXS,NRW1))
!-----------------------------------------------------------------------------------------------!
      DO K=1,NXS
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
         WRITE(32,100) 'ZONE T="BLADE WAKE',1,'" I=',NRW1,' F=POINT'
         DO J=1,NRW1
            XPWM=XPW(1,J)
            IMIN=1
            DO I=2,NPW1
               IF (XPW(I,J) < XPWM) THEN
                  XPWM=XPW(I,J)
                  IMIN=I
               END IF !(XPW(I,J) < XPWM)
            END DO !I=2,NPW1
            H0(1:(NPW1-IMIN+1))=XPW(IMIN:NPW1,J)
            H1(1:(NPW1-IMIN+1))=YPW(IMIN:NPW1,J)
            H2(1:(NPW1-IMIN+1))=ZPW(IMIN:NPW1,J)
            DO I=IMIN,NPW1-1
               H3(I-IMIN+1)=XPW(I+1,J)-XPW(I,J)
               WORK1=DATAN2(ZPW(I  ,J),YPW(I  ,J))
               WORK2=DATAN2(ZPW(I+1,J),YPW(I+1,J))
               IF (WORK2 < WORK1) WORK2=WORK2+2.D0*PI
               H4(I-IMIN+1)=WORK2-WORK1
            END DO !I=IMIN,NPW1-1
            H3(NPW1-IMIN+1)=XPW(NPW1,J)-XPW(NPW1-1,J)
            !WORK1=DATAN2(ZPW(NPW1-1,J),YPW(NPW1-1,J))
            !WORK2=DATAN2(ZPW(NPW1  ,J),YPW(NPW1  ,J))
            !IF (WORK2 < WORK1) WORK2=WORK2+2.D0*PI
            H4(NPW1-IMIN+1)=DATAN2(ZPW(NPW1,J),YPW(NPW1,J))-DATAN2(ZPW(NPW1-1,J),YPW(NPW1-1,J))
!-----------------------------------------------------------------------------------------------!
            IF (INTER == 0) THEN
               CALL LININT((NPW1-IMIN+1),H0,H1,1,XS(K),H21(K,J))
               CALL LININT((NPW1-IMIN+1),H0,H2,1,XS(K),H22(K,J))
               CALL LININT((NPW1-IMIN+1),H0,H3,1,XS(K),H23(K,J))
               CALL LININT((NPW1-IMIN+1),H0,H4,1,XS(K),H24(K,J))
            ELSEIF (INTER == 1) THEN
               CALL INTK1 ((NPW1-IMIN+1),H0,H1,1,XS(K),H21(K,J))
               CALL INTK1 ((NPW1-IMIN+1),H0,H2,1,XS(K),H22(K,J))
               CALL INTK1 ((NPW1-IMIN+1),H0,H3,1,XS(K),H23(K,J))
               CALL INTK1 ((NPW1-IMIN+1),H0,H4,1,XS(K),H24(K,J))
            ELSEIF (INTER == 2) THEN
               CALL SPLINT((NPW1-IMIN+1),H0,H1,1,XS(K),H21(K,J))
               CALL SPLINT((NPW1-IMIN+1),H0,H2,1,XS(K),H22(K,J))
               CALL SPLINT((NPW1-IMIN+1),H0,H3,1,XS(K),H23(K,J))
               CALL SPLINT((NPW1-IMIN+1),H0,H4,1,XS(K),H24(K,J))
            END IF !(INTER)
!-----------------------------------------------------------------------------------------------!
            WORK1=DSQRT(H21(K,J)*H21(K,J)+H22(K,J)*H22(K,J))
            WRITE(32,206) XS(K),H21(K,J),H22(K,J),WORK1,PI*H23(K,J)/H24(K,J), &
                                                                DATAN2D(H23(K,J),H24(K,J)*WORK1)
         END DO !J=1,NRW1
         DO I=2,NB
            WRITE(32,100) 'ZONE T="BLADE WAKE',I,'" I=',NRW1,' F=POINT'
            DO J=1,NRW1
               WORK1=DSQRT(H21(K,J)*H21(K,J)+H22(K,J)*H22(K,J))
               WORK2=DATAN2(H22(K,J),H21(K,J))
               WORK3=DFLOAT(I-1)/DFLOAT(NB)*2.D0*PI
               WRITE(32,206) XS(K),WORK1*DCOS(WORK2+WORK3),WORK1*DSIN(WORK2+WORK3), &
                                     WORK1,PI*H23(K,J)/H24(K,J),DATAN2D(H23(K,J),H24(K,J)*WORK1)
            END DO !J=1,NRW1
         END DO !I=2,NB
      END DO !K=1,NXS
      DEALLOCATE(H0,H1,H2,H3,H4,H21,H22,H23,H24)
   END IF !(IP == 1)
END DO !TT=0,NT
CLOSE(UNIT=32)
!-----------------------------------------------------------------------------------------------!
!    Formats                                                                                    !
!-----------------------------------------------------------------------------------------------!
100 FORMAT(A,I4,A,I4,A)
206 FORMAT(6(2X,E23.16))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PLOTWAKE
!-----------------------------------------------------------------------------------------------!
