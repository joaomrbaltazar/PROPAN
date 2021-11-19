!-----------------------------------------------------------------------------------------------!
!    Plot solutions                                                                             !
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
SUBROUTINE PLOTSOL(JJ)
!-----------------------------------------------------------------------------------------------!
!    Created by: 02072014, J. Baltazar, 2014 version 1.2                                        !
!    Modified  : 07012015, J. Baltazar, 2015 version 1.0                                        !
!    Modified  : 09032015, J. Baltazar, 2015 version 1.1 Unsteady Flow                          !
!    Modified  : 23052016, J. Baltazar, 2016 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPOST_MOD
IMPLICIT NONE
CHARACTER OUTFILE*50,ICHAR*3
INTEGER :: I,J,KB,JJ,TT
DOUBLE PRECISION :: XXP,YYP,ZZP,RRP,VX,VY,VZ,VR,VT
!-----------------------------------------------------------------------------------------------!
OUTFILE='SOLTOT_'//IDENTU(JJ)//'.DAT'
OPEN(UNIT=34,FILE=OUTFILE,STATUS='UNKNOWN')
DO TT=0,NT
   WRITE(ICHAR,'(I3)') TT
   WRITE(34,'(A)') 'TITLE = "PROPAN_'//IDENTU(JJ)//'_T'//TRIM(ADJUSTL(ICHAR))//'"'
   WRITE(34,'(A)') ' VARIABLES="X","Y","Z","S/C","R","POT","CP","CPN","VX","VY","VZ"'
!-----------------------------------------------------------------------------------------------!
!    Blade Solution                                                                             !
!-----------------------------------------------------------------------------------------------!
   IF (IP == 1) THEN
      KB=1
      WRITE(34,100) 'ZONE T = "Blade',KB,' T='//TRIM(ADJUSTL(ICHAR))//'", I=',NCP,' J=',NRP, &
                                                                                      ' F=POINT'
      DO J=1,NRP
         DO I=1,NCP
            WORK1=DSQRT(YP0(I,J)*YP0(I,J)+ZP0(I,J)*ZP0(I,J))
            WORK2=DATAN2(ZP0(I,J),YP0(I,J))
            VX   =VTT1P(I,J,TT)*ET1XP(I,J)+VTT2P(I,J,TT)*ET2XP(I,J)
	    VY   =VTT1P(I,J,TT)*ET1YP(I,J)+VTT2P(I,J,TT)*ET2YP(I,J)
	    VZ   =VTT1P(I,J,TT)*ET1ZP(I,J)+VTT2P(I,J,TT)*ET2ZP(I,J)
	    VR   = VY*DCOS(WORK2)+VZ*DSIN(WORK2)
            VT   =-VY*DSIN(WORK2)+VZ*DCOS(WORK2) 
	    XXP  =XP0(I,J)
            YYP  =WORK1*DCOS(WORK2)
            ZZP  =WORK1*DSIN(WORK2)
	    RRP  =DSQRT(YYP*YYP+ZZP*ZZP)
            VY   =VR*DCOS(WORK2)-VT*DSIN(WORK2)
            VZ   =VR*DSIN(WORK2)+VT*DCOS(WORK2)	
	    WRITE(34,110) XXP,YYP,ZZP, &
                                     SCP(I,J),RRP,POTP(I,J,TT),CPP(I,J,TT),CPNP(I,J,TT),VX,VY,VZ
         END DO !I=1,NCP
      END DO !J=1,NRP
!-----------------------------------------------------------------------------------------------!
!    Other Blades                                                                               !
!-----------------------------------------------------------------------------------------------!
      ALLOCATE(WORKA1(NCP,NRP),WORKA2(NCP,NRP),WORKA3(NCP,NRP),WORKA4(NCP,NRP),WORKA5(NCP,NRP))
      DO KB=2,NB
         CALL PERIODICFLOW(TT,KB,NCP,NRP,NT,POTP ,WORKA1)
         CALL PERIODICFLOW(TT,KB,NCP,NRP,NT,VTT1P,WORKA2)
         CALL PERIODICFLOW(TT,KB,NCP,NRP,NT,VTT2P,WORKA3)
         CALL PERIODICFLOW(TT,KB,NCP,NRP,NT,CPP  ,WORKA4)
         CALL PERIODICFLOW(TT,KB,NCP,NRP,NT,CPNP ,WORKA5)
         WRITE(34,100) 'ZONE T = "Blade',KB,' T='//TRIM(ADJUSTL(ICHAR))//'", I=',NCP, &
                                                                            ' J=',NRP,' F=POINT'
         DO J=1,NRP
            DO I=1,NCP
               WORK1=DSQRT(YP0(I,J)*YP0(I,J)+ZP0(I,J)*ZP0(I,J))
               WORK2=DATAN2(ZP0(I,J),YP0(I,J))
               VX   =WORKA2(I,J)*ET1XP(I,J)+WORKA3(I,J)*ET2XP(I,J)
	       VY   =WORKA2(I,J)*ET1YP(I,J)+WORKA3(I,J)*ET2YP(I,J)
	       VZ   =WORKA2(I,J)*ET1ZP(I,J)+WORKA3(I,J)*ET2ZP(I,J)
	       VR   = VY*DCOS(WORK2)+VZ*DSIN(WORK2)
               VT   =-VY*DSIN(WORK2)+VZ*DCOS(WORK2) 
	       WORK2=WORK2+DFLOAT(KB-1)*2.D0*PI/DFLOAT(NB)
	       XXP  =XP0(I,J)
               YYP  =WORK1*DCOS(WORK2)
               ZZP  =WORK1*DSIN(WORK2)
	       RRP  =DSQRT(YYP*YYP+ZZP*ZZP)
               VY   =VR*DCOS(WORK2)-VT*DSIN(WORK2)
               VZ   =VR*DSIN(WORK2)+VT*DCOS(WORK2)	
	       WRITE(34,110) XXP,YYP,ZZP, &
                                       SCP(I,J),RRP,WORKA1(I,J),WORKA4(I,J),WORKA5(I,J),VX,VY,VZ
            END DO !I=1,NCP
         END DO !J=1,NRP
      END DO !KB=2,NB
      DEALLOCATE(WORKA1,WORKA2,WORKA3,WORKA4,WORKA5)
   END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Solution                                                                            !
!-----------------------------------------------------------------------------------------------!
   IF (IN == 1) THEN
      WRITE(34,100) 'ZONE T = "Nozzle T='//TRIM(ADJUSTL(ICHAR))//'", I=',(NNXT1+1),' J=', &
                                                                          (NB*NNTP+1),' F=POINT'
      KB=1
      DO J=1,NNTP
         DO I=1,NNXT1
            WORK1=DSQRT(YN0(I,J)*YN0(I,J)+ZN0(I,J)*ZN0(I,J))
            WORK2=DATAN2(ZN0(I,J),YN0(I,J))
            VX   =VTT1N(I,J,TT)*ET1XN(I,J)+VTT2N(I,J,TT)*ET2XN(I,J)
            VY   =VTT1N(I,J,TT)*ET1YN(I,J)+VTT2N(I,J,TT)*ET2YN(I,J)
            VZ   =VTT1N(I,J,TT)*ET1ZN(I,J)+VTT2N(I,J,TT)*ET2ZN(I,J)
            VR   = VY*DCOS(WORK2)+VZ*DSIN(WORK2)
            VT   =-VY*DSIN(WORK2)+VZ*DCOS(WORK2)
            WORK3=WORK2*180.D0/PI 
            XXP  =XN0(I,J)
            YYP  =WORK1*DCOS(WORK2)
            ZZP  =WORK1*DSIN(WORK2)
            VY   =VR*DCOS(WORK2)-VT*DSIN(WORK2)
            VZ   =VR*DSIN(WORK2)+VT*DCOS(WORK2)    
            WRITE(34,110) XXP,YYP,ZZP, &
                                   SCN(I,J),WORK3,POTN(I,J,TT),CPN(I,J,TT),CPNN(I,J,TT),VX,VY,VZ
         END DO !I=1,NNXT1
!-----------------------------------------------------------------------------------------------!
         I=1
         WORK1=DSQRT(YN0(I,J)*YN0(I,J)+ZN0(I,J)*ZN0(I,J))
         WORK2=DATAN2(ZN0(I,J),YN0(I,J))
         VX   =VTT1N(I,J,TT)*ET1XN(I,J)+VTT2N(I,J,TT)*ET2XN(I,J)
         VY   =VTT1N(I,J,TT)*ET1YN(I,J)+VTT2N(I,J,TT)*ET2YN(I,J)
         VZ   =VTT1N(I,J,TT)*ET1ZN(I,J)+VTT2N(I,J,TT)*ET2ZN(I,J)
         VR   = VY*DCOS(WORK2)+VZ*DSIN(WORK2)
         VT   =-VY*DSIN(WORK2)+VZ*DCOS(WORK2)
         WORK3=WORK2*180.D0/PI 
         XXP  =XN0(I,J)
         YYP  =WORK1*DCOS(WORK2)
         ZZP  =WORK1*DSIN(WORK2)
         VY   =VR*DCOS(WORK2)-VT*DSIN(WORK2)
         VZ   =VR*DSIN(WORK2)+VT*DCOS(WORK2)    
         WRITE(34,110) XXP,YYP,ZZP,SCN(I,J),WORK3,POTN(I,J,TT),CPN(I,J,TT),CPNN(I,J,TT),VX,VY,VZ
      END DO !J=1,NNTP
!-----------------------------------------------------------------------------------------------!
!    Other Blades                                                                               !
!-----------------------------------------------------------------------------------------------!
      ALLOCATE(WORKA1(NNXT1,NNTP),WORKA2(NNXT1,NNTP),WORKA3(NNXT1,NNTP))
      ALLOCATE(WORKA4(NNXT1,NNTP),WORKA5(NNXT1,NNTP))
      DO KB=2,NB
         CALL PERIODICFLOW(TT,KB,NNXT1,NNTP,NT,POTN ,WORKA1)
         CALL PERIODICFLOW(TT,KB,NNXT1,NNTP,NT,VTT1N,WORKA2)
         CALL PERIODICFLOW(TT,KB,NNXT1,NNTP,NT,VTT2N,WORKA3)
         CALL PERIODICFLOW(TT,KB,NNXT1,NNTP,NT,CPN  ,WORKA4)
         CALL PERIODICFLOW(TT,KB,NNXT1,NNTP,NT,CPNN ,WORKA5)
         DO J=1,NNTP
            DO I=1,NNXT1
               WORK1=DSQRT(YN0(I,J)*YN0(I,J)+ZN0(I,J)*ZN0(I,J))
               WORK2=DATAN2(ZN0(I,J),YN0(I,J))
               VX   =WORKA2(I,J)*ET1XN(I,J)+WORKA3(I,J)*ET2XN(I,J)
               VY   =WORKA2(I,J)*ET1YN(I,J)+WORKA3(I,J)*ET2YN(I,J)
               VZ   =WORKA2(I,J)*ET1ZN(I,J)+WORKA3(I,J)*ET2ZN(I,J)
               VR   = VY*DCOS(WORK2)+VZ*DSIN(WORK2)
               VT   =-VY*DSIN(WORK2)+VZ*DCOS(WORK2)
               WORK2=WORK2+DFLOAT(KB-1)*2.D0*PI/DFLOAT(NB)
               WORK3=WORK2*180.D0/PI 
               XXP  =XN0(I,J)
               YYP  =WORK1*DCOS(WORK2)
               ZZP  =WORK1*DSIN(WORK2)
               VY   =VR*DCOS(WORK2)-VT*DSIN(WORK2)
               VZ   =VR*DSIN(WORK2)+VT*DCOS(WORK2)    
               WRITE(34,110) XXP,YYP,ZZP, &
                                     SCN(I,J),WORK3,WORKA1(I,J),WORKA4(I,J),WORKA5(I,J),VX,VY,VZ
            END DO !I=1,NNXT1
!-----------------------------------------------------------------------------------------------!
            I=1
            WORK1=DSQRT(YN0(I,J)*YN0(I,J)+ZN0(I,J)*ZN0(I,J))
            WORK2=DATAN2(ZN0(I,J),YN0(I,J))
            VX   =WORKA2(I,J)*ET1XN(I,J)+WORKA3(I,J)*ET2XN(I,J)
            VY   =WORKA2(I,J)*ET1YN(I,J)+WORKA3(I,J)*ET2YN(I,J)
            VZ   =WORKA2(I,J)*ET1ZN(I,J)+WORKA3(I,J)*ET2ZN(I,J)
            VR   = VY*DCOS(WORK2)+VZ*DSIN(WORK2)
            VT   =-VY*DSIN(WORK2)+VZ*DCOS(WORK2)
            WORK2=WORK2+DFLOAT(KB-1)*2.D0*PI/DFLOAT(NB)
            WORK3=WORK2*180.D0/PI 
            XXP  =XN0(I,J)
            YYP  =WORK1*DCOS(WORK2)
            ZZP  =WORK1*DSIN(WORK2)
            VY   =VR*DCOS(WORK2)-VT*DSIN(WORK2)
            VZ   =VR*DSIN(WORK2)+VT*DCOS(WORK2)    
            WRITE(34,110) XXP,YYP,ZZP, &
                                     SCN(I,J),WORK3,WORKA1(I,J),WORKA4(I,J),WORKA5(I,J),VX,VY,VZ
         END DO !J=1,NNTP
      END DO !KB=2,NB
      DEALLOCATE(WORKA1,WORKA2,WORKA3,WORKA4,WORKA5)
!-----------------------------------------------------------------------------------------------!
!    Repeat First Line                                                                          !
!-----------------------------------------------------------------------------------------------!
      J=1
      DO I=1,NNXT1
         WORK1=DSQRT(YN0(I,J)*YN0(I,J)+ZN0(I,J)*ZN0(I,J))
         WORK2=DATAN2(ZN0(I,J),YN0(I,J))
         VX   =VTT1N(I,J,TT)*ET1XN(I,J)+VTT2N(I,J,TT)*ET2XN(I,J)
         VY   =VTT1N(I,J,TT)*ET1YN(I,J)+VTT2N(I,J,TT)*ET2YN(I,J)
         VZ   =VTT1N(I,J,TT)*ET1ZN(I,J)+VTT2N(I,J,TT)*ET2ZN(I,J)
         VR   = VY*DCOS(WORK2)+VZ*DSIN(WORK2)
         VT   =-VY*DSIN(WORK2)+VZ*DCOS(WORK2)
         WORK3=WORK2*180.D0/PI 
         XXP  =XN0(I,J)
         YYP  =WORK1*DCOS(WORK2)
         ZZP  =WORK1*DSIN(WORK2)
         VY   =VR*DCOS(WORK2)-VT*DSIN(WORK2)
         VZ   =VR*DSIN(WORK2)+VT*DCOS(WORK2)    
         WRITE(34,110) XXP,YYP,ZZP,SCN(I,J),WORK3,POTN(I,J,TT),CPN(I,J,TT),CPNN(I,J,TT),VX,VY,VZ
      END DO !I=1,NNXT1
!-----------------------------------------------------------------------------------------------!
      I=1
      WORK1=DSQRT(YN0(I,J)*YN0(I,J)+ZN0(I,J)*ZN0(I,J))
      WORK2=DATAN2(ZN0(I,J),YN0(I,J))
      VX   =VTT1N(I,J,TT)*ET1XN(I,J)+VTT2N(I,J,TT)*ET2XN(I,J)
      VY   =VTT1N(I,J,TT)*ET1YN(I,J)+VTT2N(I,J,TT)*ET2YN(I,J)
      VZ   =VTT1N(I,J,TT)*ET1ZN(I,J)+VTT2N(I,J,TT)*ET2ZN(I,J)
      VR   = VY*DCOS(WORK2)+VZ*DSIN(WORK2)
      VT   =-VY*DSIN(WORK2)+VZ*DCOS(WORK2)
      WORK3=WORK2*180.D0/PI 
      XXP  =XN0(I,J)
      YYP  =WORK1*DCOS(WORK2)
      ZZP  =WORK1*DSIN(WORK2)
      VY   =VR*DCOS(WORK2)-VT*DSIN(WORK2)
      VZ   =VR*DSIN(WORK2)+VT*DCOS(WORK2)    
      WRITE(34,110) XXP,YYP,ZZP,SCN(I,J),WORK3,POTN(I,J,TT),CPN(I,J,TT),CPNN(I,J,TT),VX,VY,VZ
   END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub Solution                                                                               !
!-----------------------------------------------------------------------------------------------!
   IF (IABS(IH) == 1) THEN
      WRITE(34,100) 'ZONE T = "Hub T='//TRIM(ADJUSTL(ICHAR))//'", I=',(NHX+2),' J=', &
                                                                          (NB*NHTP+1),' F=POINT'
      KB=1
      DO J=1,NHTP
         XXP  =XH0(1,J)
         YYP  =0.D0
	 ZZP  =0.D0
	 WORK2=0.D0
	 WORK3=1.D0
	 WRITE(34,110) XXP,YYP,ZZP,XXP,WORK2,WORK2,WORK3,WORK2,WORK2,WORK2,WORK2
	 DO I=1,NHX
            WORK1=DSQRT(YH0(I,J)*YH0(I,J)+ZH0(I,J)*ZH0(I,J))
            WORK2=DATAN2(ZH0(I,J),YH0(I,J))
            VX   =VTT1H(I,J,TT)*ET1XH(I,J)+VTT2H(I,J,TT)*ET2XH(I,J)
            VY   =VTT1H(I,J,TT)*ET1YH(I,J)+VTT2H(I,J,TT)*ET2YH(I,J)
            VZ   =VTT1H(I,J,TT)*ET1ZH(I,J)+VTT2H(I,J,TT)*ET2ZH(I,J)
            VR   = VY*DCOS(WORK2)+VZ*DSIN(WORK2)
            VT   =-VY*DSIN(WORK2)+VZ*DCOS(WORK2) 
            WORK3=WORK2*180.D0/PI 
            XXP  =XH0(I,J)
            YYP  =WORK1*DCOS(WORK2)
            ZZP  =WORK1*DSIN(WORK2)
            RRP  =DSQRT(YYP*YYP+ZZP*ZZP)
            VY   =VR*DCOS(WORK2)-VT*DSIN(WORK2)
            VZ   =VR*DSIN(WORK2)+VT*DCOS(WORK2)	
            WRITE(34,110) XXP,YYP,ZZP,XXP,WORK3,POTH(I,J,TT),CPH(I,J,TT),CPNH(I,J,TT),VX,VY,VZ
         END DO !I=1,NHX
         XXP  =XH0(NHX,J)
	 YYP  =0.D0
	 ZZP  =0.D0
	 WORK2=0.D0
	 WORK3=1.D0
	 WRITE(34,110) XXP,YYP,ZZP,XXP,WORK2,WORK2,WORK3,WORK2,WORK2,WORK2,WORK2
      END DO !J=1,NHTP
!-----------------------------------------------------------------------------------------------!
!    Other Blades                                                                               !
!-----------------------------------------------------------------------------------------------!
      ALLOCATE(WORKA1(NHX,NHTP),WORKA2(NHX,NHTP),WORKA3(NHX,NHTP))
      ALLOCATE(WORKA4(NHX,NHTP),WORKA5(NHX,NHTP))
      DO KB=2,NB
         CALL PERIODICFLOW(TT,KB,NHX,NHTP,NT,POTH ,WORKA1)
         CALL PERIODICFLOW(TT,KB,NHX,NHTP,NT,VTT1H,WORKA2)
         CALL PERIODICFLOW(TT,KB,NHX,NHTP,NT,VTT2H,WORKA3)
         CALL PERIODICFLOW(TT,KB,NHX,NHTP,NT,CPH  ,WORKA4)
         CALL PERIODICFLOW(TT,KB,NHX,NHTP,NT,CPNH ,WORKA5)
         DO J=1,NHTP
            XXP  =XH0(1,J)
	    YYP  =0.D0
	    ZZP  =0.D0
	    WORK2=0.D0
	    WORK3=1.D0
	    WRITE(34,110) XXP,YYP,ZZP,XXP,WORK2,WORK2,WORK3,WORK2,WORK2,WORK2,WORK2
	    DO I=1,NHX
               WORK1=DSQRT(YH0(I,J)*YH0(I,J)+ZH0(I,J)*ZH0(I,J))
               WORK2=DATAN2(ZH0(I,J),YH0(I,J))
               VX   =WORKA2(I,J)*ET1XH(I,J)+WORKA3(I,J)*ET2XH(I,J)
               VY   =WORKA2(I,J)*ET1YH(I,J)+WORKA3(I,J)*ET2YH(I,J)
               VZ   =WORKA2(I,J)*ET1ZH(I,J)+WORKA3(I,J)*ET2ZH(I,J)
               VR   = VY*DCOS(WORK2)+VZ*DSIN(WORK2)
               VT   =-VY*DSIN(WORK2)+VZ*DCOS(WORK2) 
	       WORK2=WORK2+DFLOAT(KB-1)*2.D0*PI/DFLOAT(NB)
               WORK3=WORK2*180.D0/PI 
               XXP  =XH0(I,J)
               YYP  =WORK1*DCOS(WORK2)
               ZZP  =WORK1*DSIN(WORK2)
               RRP  =DSQRT(YYP*YYP+ZZP*ZZP)
               VY   =VR*DCOS(WORK2)-VT*DSIN(WORK2)
               VZ   =VR*DSIN(WORK2)+VT*DCOS(WORK2)	
               WRITE(34,110) XXP,YYP,ZZP,XXP,WORK3,WORKA1(I,J),WORKA4(I,J),WORKA5(I,J),VX,VY,VZ
            END DO !I=1,NHX
	    XXP  =XH0(NHX,J)
	    YYP  =0.D0
	    ZZP  =0.D0
	    WORK2=0.D0
	    WORK3=1.D0
	    WRITE(34,110) XXP,YYP,ZZP,XXP,WORK2,WORK2,WORK3,WORK2,WORK2,WORK2,WORK2
         END DO !J=1,NHTP
      END DO !KB=2,NB
      DEALLOCATE(WORKA1,WORKA2,WORKA3,WORKA4,WORKA5)
!-----------------------------------------------------------------------------------------------!
!    Repeat First Line                                                                          !
!-----------------------------------------------------------------------------------------------!
      J=1
      XXP  =XH0(1,J)
      YYP  =0.D0
      ZZP  =0.D0
      WORK2=0.D0
      WORK3=1.D0
      WRITE(34,110) XXP,YYP,ZZP,XXP,WORK2,WORK2,WORK3,WORK2,WORK2,WORK2,WORK2
      DO I=1,NHX
         WORK1=DSQRT(YH0(I,J)*YH0(I,J)+ZH0(I,J)*ZH0(I,J))
         WORK2=DATAN2(ZH0(I,J),YH0(I,J))
         VX   =VTT1H(I,J,TT)*ET1XH(I,J)+VTT2H(I,J,TT)*ET2XH(I,J)
         VY   =VTT1H(I,J,TT)*ET1YH(I,J)+VTT2H(I,J,TT)*ET2YH(I,J)
         VZ   =VTT1H(I,J,TT)*ET1ZH(I,J)+VTT2H(I,J,TT)*ET2ZH(I,J)
         VR   = VY*DCOS(WORK2)+VZ*DSIN(WORK2)
         VT   =-VY*DSIN(WORK2)+VZ*DCOS(WORK2) 
         WORK3=WORK2*180.D0/PI
         XXP  =XH0(I,J)
         YYP  =WORK1*DCOS(WORK2)
         ZZP  =WORK1*DSIN(WORK2)
         RRP  =DSQRT(YYP*YYP+ZZP*ZZP)
         VY   =VR*DCOS(WORK2)-VT*DSIN(WORK2)
         VZ   =VR*DSIN(WORK2)+VT*DCOS(WORK2)	
         WRITE(34,110) XXP,YYP,ZZP,XXP,WORK3,POTH(I,J,TT),CPH(I,J,TT),CPNH(I,J,TT),VX,VY,VZ
      END DO !I=1,NHX
      XXP  =XH0(NHX,J)
      YYP  =0.D0
      ZZP  =0.D0
      WORK2=0.D0
      WORK3=1.D0
      WRITE(34,110) XXP,YYP,ZZP,XXP,WORK2,WORK2,WORK3,WORK2,WORK2,WORK2,WORK2
   END IF !(IABS(IH) == 1)
END DO !TT=0,NT
CLOSE(UNIT=34)
!-----------------------------------------------------------------------------------------------!
!    Formats                                                                                    !
!-----------------------------------------------------------------------------------------------!
100 FORMAT(A,I4,A,I4,A,I4,A,I4,A)
110 FORMAT(11(2X,E23.16))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE PLOTSOL
!-----------------------------------------------------------------------------------------------!
