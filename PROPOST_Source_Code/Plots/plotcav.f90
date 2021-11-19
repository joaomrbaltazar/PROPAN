!-----------------------------------------------------------------------------------------------!
!    Plot cavitation results                                                                    !
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
SUBROUTINE PLOTCAV(JJ)
!-----------------------------------------------------------------------------------------------!
!    Created by: 23052016, J. Baltazar, 2016 version 1.0                                        !
!    Modified  : 06072016, J. Baltazar, 2016 version 1.0                                        !
!    Modified  : 03022017, J. Baltazar, 2017 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPOST_MOD
IMPLICIT NONE
CHARACTER INFILE*50,OUTFILE*50,ICHAR*3
LOGICAL :: EXISTS
INTEGER :: I,J,KB,JJ,TT
DOUBLE PRECISION :: XXP,YYP,ZZP
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: CAVP,CAVPW
!-----------------------------------------------------------------------------------------------!
!    Inquire File                                                                               !
!-----------------------------------------------------------------------------------------------!
INFILE='CAV_'//IDENTU(JJ)//'.DAT'
INQUIRE(FILE=INFILE,EXIST=EXISTS)
IF (EXISTS) THEN
   IF (IP == 1) ALLOCATE(CAVP (NCP,NRP,0:NT))
   IF ((IP == 1).AND.(NCPW /= 0)) ALLOCATE(CAVPW(IABS(NCPW),NRW,0:NT))
!-----------------------------------------------------------------------------------------------!
!    Read File                                                                                  !
!-----------------------------------------------------------------------------------------------!
   OPEN(UNIT=24,FILE=INFILE,STATUS='UNKNOWN')
   READ(24,*)
   READ(24,*)
   DO TT=0,NT
!-----------------------------------------------------------------------------------------------!
!    Blade Cavity                                                                               !
!-----------------------------------------------------------------------------------------------!
      IF (IP == 1) THEN
         READ(24,*)
         DO J=1,NRP
            DO I=1,NCP
               READ(24,*) WORK1,WORK2,WORK3,CAVP(I,J,TT)
            END DO !I=1,NCP
         END DO !J=NRP
      END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake Cavity                                                                          !
!-----------------------------------------------------------------------------------------------!
      IF ((IP == 1).AND.(NCPW /= 0)) THEN
         READ(24,*)
         DO J=1,NRW
            DO I=1,IABS(NCPW)
               READ(24,*) WORK1,WORK2,WORK3,CAVPW(I,J,TT)
            END DO !I=1,IABS(NCPW)
         END DO !J=NRW
      END IF !((IP == 1).AND.(NCPW /= 0))
   END DO !TT=0,NT
   CLOSE(UNIT=24)
!-----------------------------------------------------------------------------------------------!
!    Output File                                                                                !
!-----------------------------------------------------------------------------------------------!
   OUTFILE='CAVTOT_'//IDENTU(JJ)//'.DAT'
   OPEN(UNIT=35,FILE=OUTFILE,STATUS='UNKNOWN')
   DO TT=0,NT
!-----------------------------------------------------------------------------------------------!
!    Write File                                                                                 !
!-----------------------------------------------------------------------------------------------!
      WRITE(ICHAR,'(I3)') TT
      WRITE(35,'(A)') 'TITLE = "PROPAN_'//IDENTU(JJ)//'_T'//TRIM(ADJUSTL(ICHAR))//'"'
      WRITE(35,'(A)') ' VARIABLES="X","Y","Z","CAV"'
!-----------------------------------------------------------------------------------------------!
!    Blade Solution                                                                             !
!-----------------------------------------------------------------------------------------------!
      IF (IP == 1) THEN
         KB=1
         WRITE(35,100) 'ZONE T = "Blade',KB,' T='//TRIM(ADJUSTL(ICHAR))//'", I=',NCP, &
                                                                            ' J=',NRP,' F=POINT'
         DO J=1,NRP
            DO I=1,NCP
               WRITE(35,120) XP0(I,J),YP0(I,J),ZP0(I,J),CAVP(I,J,TT)
            END DO !I=1,NCP
         END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
!    Other Blades                                                                               !
!-----------------------------------------------------------------------------------------------!
         ALLOCATE(WORKA1(NCP,NRP))
         DO KB=2,NB
            CALL PERIODICFLOW(TT,KB,NCP,NRP,NT,CAVP,WORKA1)
            WRITE(35,100) 'ZONE T = "Blade',KB,' T='//TRIM(ADJUSTL(ICHAR))//'", I=',NCP,' J=' &
                                                                                 ,NRP,' F=POINT' 
            DO J=1,NRP
               DO I=1,NCP
                  WORK1=DSQRT(YP0(I,J)*YP0(I,J)+ZP0(I,J)*ZP0(I,J))
                  WORK2=DATAN2(ZP0(I,J),YP0(I,J))
	          WORK2=WORK2+DFLOAT(KB-1)*2.D0*PI/DFLOAT(NB)
	          XXP  =XP0(I,J)
                  YYP  =WORK1*DCOS(WORK2)
                  ZZP  =WORK1*DSIN(WORK2)
	          WRITE(35,120) XXP,YYP,ZZP,WORKA1(I,J)
               END DO !I=1,NCP
            END DO !J=1,NRP
         END DO !KB=2,NB
         DEALLOCATE(WORKA1)
      END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake Solution                                                                        !
!-----------------------------------------------------------------------------------------------!
      IF ((IP == 1).AND.(NCPW /= 0)) THEN
         KB=1
         WRITE(35,100) 'ZONE T = "Blade Wake',KB,' T='//TRIM(ADJUSTL(ICHAR))//'", I=' &
                                                                ,IABS(NCPW),' J=',NRW,' F=POINT'
         DO J=1,NRW
            DO I=1,IABS(NCPW)
               WRITE(35,120) XPW0(I,J),YPW0(I,J),ZPW0(I,J),CAVPW(I,J,TT)
            END DO !I=1,IABS(NCPW)
         END DO !J=1,NRW
!-----------------------------------------------------------------------------------------------!
!    Other Blade Wakes                                                                          !
!-----------------------------------------------------------------------------------------------!
         ALLOCATE(WORKA1(IABS(NCPW),NRW))
         DO KB=2,NB
            CALL PERIODICFLOW(TT,KB,IABS(NCPW),NRW,NT,CAVPW,WORKA1)
            WRITE(35,100) 'ZONE T = "Blade Wake',KB,' T='//TRIM(ADJUSTL(ICHAR))//'", I=' &
                                                                ,IABS(NCPW),' J=',NRW,' F=POINT'
            DO J=1,NRW
               DO I=1,IABS(NCPW)
                  WORK1=DSQRT(YPW0(I,J)*YPW0(I,J)+ZPW0(I,J)*ZPW0(I,J))
                  WORK2=DATAN2(ZPW0(I,J),YPW0(I,J))
	          WORK2=WORK2+DFLOAT(KB-1)*2.D0*PI/DFLOAT(NB)
	          XXP  =XPW0(I,J)
                  YYP  =WORK1*DCOS(WORK2)
                  ZZP  =WORK1*DSIN(WORK2)
	          WRITE(35,120) XXP,YYP,ZZP,WORKA1(I,J)
               END DO !I=1,IBAS(NCPW)
            END DO !J=1,NRW
         END DO !KB=2,NB
         DEALLOCATE(WORKA1)
      END IF !((IP == 1).AND.(NCPW /= 0))
!-----------------------------------------------------------------------------------------------!
   END DO !TT=0,NT
   IF (IP == 1) DEALLOCATE(CAVP)
   IF ((IP == 1).AND.(NCPW /= 0)) DEALLOCATE(CAVPW)
END IF !EXISTS
CLOSE(UNIT=35)
!-----------------------------------------------------------------------------------------------!
!    Formats                                                                                    !
!-----------------------------------------------------------------------------------------------!
100 FORMAT(A,I4,A,I4,A,I4,A,I4,A)
120 FORMAT(4(2X,E23.16))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PLOTCAV
!-----------------------------------------------------------------------------------------------!
