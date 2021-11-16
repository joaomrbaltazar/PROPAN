!-----------------------------------------------------------------------------------------------!
!    Nozzle geometrical definition                                                              !
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
SUBROUTINE NOZZLEDEF(SIDE,XL,RL)
!-----------------------------------------------------------------------------------------------!
!    Created by:  J. Baltazar, IST, November 2010                                               !
!    Modified  : 15102013, J. Baltazar, version 1.0                                             !
!    Modified  : 21102013, J. Baltazar, differences lower than the tolerance for t.e.           !
!    Modified  : 26052014, J. Baltazar, renamed for ProPan                                      !
!    Modified  : 06072017, J. Baltazar, 2017 version 1.0                                        !
!    Modified  : 03102018, J. Baltazar, 2018 version 1.1                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
CHARACTER*5 SIDE
INTEGER :: I,LIDENTN
DOUBLE PRECISION :: XL,XD,RL,TN0,FN0,RMAX
DOUBLE PRECISION :: TTETA,TETA,STETA,CTETA,SLE,RLE,STE,RTE
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: TNN,FNN
!-----------------------------------------------------------------------------------------------!
I=1
DO WHILE (IDENTN(I:I).NE.' ')
   I=I+1
END DO !(IDENTN(I:I).NE.' ')
LIDENTN=I-1
RMAX =1.D0
!-----------------------------------------------------------------------------------------------!
IF (NRNI == NRNO) THEN
   ALLOCATE(TNN(NRNI),FNN(NRNI))
   TNN=0.D0
   FNN=0.D0
!-----------------------------------------------------------------------------------------------!
   TNN(:)=0.5D0*(YOL(:)-YIL(:))
   FNN(:)=0.5D0*(YOL(:)+YIL(:))
!-----------------------------------------------------------------------------------------------!
!    Round Trailing Edge of the Nozzle                                                          !
!-----------------------------------------------------------------------------------------------!
   STE=XIL(NRNI) !0.D0
   IF ((IDENTN(1:LIDENTN) == '19A').OR.(IDENTN(1:LIDENTN) == '19Am')) THEN
      TTETA=(TNN(NRNI-1)-TNN(NRNI))/(XIL(NRNI)-XIL(NRNI-1)) !XIL=XOL for 19A
      TETA =DATAN(TTETA)
      STETA=DSIN(TETA)
      CTETA=DCOS(TETA)
      RTE  =TNN(NRNI)/(CTETA-TTETA+STETA*TTETA)
      STE  =XIL(NRNI)-RTE+RTE*STETA
   END IF !((IDENTN(1:LIDENTN) == '19A').OR.(IDENTN(1:LIDENTN) == '19Am'))
!-----------------------------------------------------------------------------------------------!
!    Round Leading Edge of the Nozzle                                                           !
!-----------------------------------------------------------------------------------------------!
   SLE=0.D0
   IF ((IDENTN(1:LIDENTN) == '19A').OR.(IDENTN(1:LIDENTN) == '19Am')) THEN
      TTETA=(TNN(4)-TNN(3))/(XIL(4)-XIL(3)) !XIL=XOL for 19A
      TETA =DATAN(TTETA)
      STETA=DSIN(TETA)
      CTETA=DCOS(TETA)
      RLE  =(TNN(3)-TTETA*XIL(3))/(CTETA-TTETA+STETA*TTETA)
      SLE  =RLE-RLE*STETA
   END IF !((IDENTN(1:LIDENTN) == '19A').OR.(IDENTN(1:LIDENTN) == '19Am'))
END IF !(NRNI == NRNO)
!-----------------------------------------------------------------------------------------------!
!    Interpolate the Radius                                                                     !
!-----------------------------------------------------------------------------------------------!
XD=0.5D0/LD*XL+0.5D0
!-----------------------------------------------------------------------------------------------!
!    Inner Side                                                                                 !
!-----------------------------------------------------------------------------------------------!
IF (SIDE == 'INNER') THEN
   IF (XD < SLE) THEN
      TN0=DSQRT(RLE**2-(XD-RLE)**2)
      IF (DABS(XD-XIL(1)) < TOL) TN0=0.D0
      IF (INTERN == 0) THEN
         CALL LININT(NRNI,XIL,FNN,1,XD,FN0)
      ELSEIF (INTERN == 1) THEN
         CALL INTK1 (NRNI,XIL,FNN,1,XD,FN0)
      ELSEIF (INTERN == 2) THEN
         CALL SPLINT(NRNI,XIL,FNN,1,XD,FN0)
      END IF !(INTERN)
      RL=FN0-TN0
   ELSEIF (XD > STE) THEN
      TN0=DSQRT(RTE**2-(XD-XIL(NRNI)+RTE)**2)
      IF (DABS(XD-XIL(NRNI)) < TOL) TN0=0.D0
      IF (INTERN == 0) THEN
         CALL LININT(NRNI,XIL,FNN,1,XD,FN0)
      ELSEIF (INTERN == 1) THEN
         CALL INTK1 (NRNI,XIL,FNN,1,XD,FN0)
      ELSEIF (INTERN == 2) THEN
         CALL SPLINT(NRNI,XIL,FNN,1,XD,FN0)
      END IF !(INTERN)
      RL=FN0-TN0
   ELSEIF (((IDENTN(1:LIDENTN) == '19A').OR. &
            (IDENTN(1:LIDENTN) == '19Am')).AND. &
           ((XD >= 0.4D0).AND.(XD <= 0.6D0))) THEN
      CALL LININT(NRNI,XIL,YIL,1,XD,RL)
   ELSEIF ((IDENTN(1:LIDENTN) == '1399').AND. &
           ((XD >= 0.30D0).AND.(XD <= 0.55D0))) THEN
      CALL LININT(NRNI,XIL,YIL,1,XD,RL)
   ELSEIF (IDENTN(1:LIDENTN) == '37') THEN
      CALL GEODUCT37(2,1,XD,RL)
   ELSE !(XD)
      IF (INTERN == 0) THEN
         CALL LININT(NRNI,XIL,YIL,1,XD,RL)
      ELSEIF (INTERN == 1) THEN
         CALL INTK1 (NRNI,XIL,YIL,1,XD,RL)
      ELSEIF (INTERN == 2) THEN
         CALL SPLINT(NRNI,XIL,YIL,1,XD,RL)
      END IF !(INTERN)
   END IF !(XD)
END IF !(SIDE == 'INNER')
!-----------------------------------------------------------------------------------------------!
!    Outer Side                                                                                 !
!-----------------------------------------------------------------------------------------------!
IF (SIDE == 'OUTER') THEN
   IF (XD < SLE) THEN
      TN0=DSQRT(RLE**2-(XD-RLE)**2)
      IF (DABS(XD-XIL(1)) < TOL) TN0=0.D0
      IF (INTERN == 0) THEN
         CALL LININT(NRNO,XOL,FNN,1,XD,FN0)
      ELSEIF (INTERN == 1) THEN
         CALL INTK1 (NRNO,XOL,FNN,1,XD,FN0)
      ELSEIF (INTERN == 2) THEN
         CALL SPLINT(NRNO,XOL,FNN,1,XD,FN0)
      END IF !(INTERN)
      RL=FN0+TN0
   ELSEIF (XD > STE) THEN
      TN0=DSQRT(RTE**2-(XD-XOL(NRNO)+RTE)**2)
      IF (DABS(XD-XIL(NRNI)) < TOL) TN0=0.D0
      IF (INTERN == 0) THEN
         CALL LININT(NRNO,XOL,FNN,1,XD,FN0)
      ELSEIF (INTERN == 1) THEN
         CALL INTK1 (NRNO,XOL,FNN,1,XD,FN0)
      ELSEIF (INTERN == 2) THEN
         CALL SPLINT(NRNO,XOL,FNN,1,XD,FN0)
      END IF !(INTERN)
      RL=FN0+TN0
   ELSEIF (IDENTN(1:LIDENTN) == '37') THEN
      CALL GEODUCT37(1,1,(1.D0-XD),RL)
   ELSE !(XD)
      IF (INTERN == 0) THEN
         CALL LININT(NRNO,XOL,YOL,1,XD,RL)
      ELSEIF (INTERN == 1) THEN
         CALL INTK1 (NRNO,XOL,YOL,1,XD,RL)
      ELSEIF (INTERN == 2) THEN
         CALL SPLINT(NRNO,XOL,YOL,1,XD,RL)
      END IF !(INTERN)
   END IF !(XD)
END IF !(SIDE == 'OUTER')
!-----------------------------------------------------------------------------------------------!
RL=RL*LD*2.D0 !non-dimensional by propeller radius
RL=RL+RMAX+CR
!-----------------------------------------------------------------------------------------------!
DEALLOCATE(TNN,FNN)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE NOZZLEDEF
!-----------------------------------------------------------------------------------------------!
