!-----------------------------------------------------------------------------------------------!
!    Generate Nozzle Wake Grid                                                                  !
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
SUBROUTINE NOZZLEWAKEGRID
!-----------------------------------------------------------------------------------------------!
!    Created by: 21102013, J. Baltazar, 2013 version 1.0                                        !
!    Modified  : 08102014, J. Baltazar, 2014 version 1.1                                        !
!    Modified  : 24062015, J. Baltazar, 2015 version 1.2                                        !
!-----------------------------------------------------------------------------------------------!
!    Declarations                                                                               !
!-----------------------------------------------------------------------------------------------!
USE PROPANEL_MOD
IMPLICIT NONE
INTEGER :: I,J,K
DOUBLE PRECISION :: XC,TC,RR,CSI,FCSI,PHI
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DSL
!-----------------------------------------------------------------------------------------------!
!    Counters                                                                                   !
!-----------------------------------------------------------------------------------------------!
IF (ISTRIP == 1) NNW=NPW-NND
NNW1=NNW+1
!-----------------------------------------------------------------------------------------------!
ALLOCATE(XNW(NNW1,NNTT),YNW(NNW1,NNTT),ZNW(NNW1,NNTT),RNW(NNW1,NNTT),TNW(NNW1,NNTT))
XNW=0.D0
YNW=0.D0
ZNW=0.D0
RNW=0.D0
TNW=0.D0
!-----------------------------------------------------------------------------------------------!
!    Stretching along the Streamwise                                                            !
!-----------------------------------------------------------------------------------------------!
IF (ISTRIP /= 1) THEN
   ALLOCATE(DSL(NNW1))
   DSL  =0.D0
   ST1NW=0.D0
   ST2NW=0.D0
   ST3NW=0.D0
   ST4NW=0.D0
   CALL STRET_CHOICE(NNW1,DSL,ST1NW,ST2NW,ST3NW,ST4NW,ITYPENW)
END IF !(ISTRIP /= 1)
!-----------------------------------------------------------------------------------------------!
DO J=1,NNTT
   DO K=1,NB
      XNW(1,J)=XN(NNX1,J)
      RNW(1,J)=RN(NNX1,J)
      TNW(1,J)=TN(NNX1,J)
      YNW(1,J)=RNW(1,J)*DCOS(TNW(1,J))
      ZNW(1,J)=RNW(1,J)*DSIN(TNW(1,J))
   END DO !K=1,NB
!-----------------------------------------------------------------------------------------------!
!    Loop on Axial Stations                                                                     !
!-----------------------------------------------------------------------------------------------!
   DO I=2,NNW1
      IF (ISTRIP == 1) THEN
         XC=XPW(NND+I,NRW1)
         TC=TNW(1,J)+(TPW(NND+I,NRW1)-TPW(NND+1,NRW1)) !2.D0*PI*(XC-XNW(1,J))/PTN
      ELSE !(ISTRIP == 1)
         IF (PTN < TOL) THEN
            XC=XNW(1,J)+DSL(I)*XNWT
            TC=TNW(1,J)
         ELSE !(PTN < TOL)
            XC=XNW(1,J)+DSL(I)*XNWT
            TC=TNW(1,J)+2.D0*PI*DSL(I)*XNWT/PTN
         END IF !(PTN < TOL)
      END IF !(ISTRIP == 1)
!-----------------------------------------------------------------------------------------------!
      IF (ICONTRNW == 1) THEN
         IF (XC < XNWW) THEN
            CSI =(XC-LD)/(XNWW-LD)
            FCSI=DSQRT(CSI)+1.013D0*CSI-1.920D0*CSI**2+1.228D0*CSI**3-0.321D0*CSI**4
            RR  =0.5D0*(YOL(NRNI)+YIL(NRNI))
            RR  =RR*LD*2.D0+RMAX+CR
            RR  =RR+(R2-RR)*FCSI
         ELSE !(XC < XNWW)
            RR  =R2
         END IF !(XC < XNWW)
      ELSE !(ICONTRNW == 1)
         RR=0.5D0*(YOL(NRNI)+YIL(NRNI))
         RR=RR*LD*2.D0+RMAX+CR
      END IF !(ICONTRNW == 1)
!-----------------------------------------------------------------------------------------------!
      XNW(I,J)=XC
      TNW(I,J)=TC
      RNW(I,J)=RR
      YNW(I,J)=RNW(I,J)*DCOS(TNW(I,J))
      ZNW(I,J)=RNW(I,J)*DSIN(TNW(I,J))
   END DO !I=2,NNW1
END DO !J=1,NNTT
!-----------------------------------------------------------------------------------------------!
!    Write Nozzle Wake Grid in Tecplot Format                                                   !
!-----------------------------------------------------------------------------------------------!
K=1
WRITE(20,100) ' ZONE T="NOZZLE WAKE',K,'" F=POINT, I=',NNW1,' J=',NNTT
DO J=1,NNTT
   DO I=1,NNW1
      WRITE(20,110) XNW(I,J),YNW(I,J),ZNW(I,J)
   END DO !I=1,NNW1
END DO !J=1,NNTT
DO K=2,NB
   PHI=DFLOAT(K-1)/DFLOAT(NB)*2.D0*PI
   WRITE(20,100) ' ZONE T="NOZZLE WAKE',K,'" F=POINT, I=',NNW1,' J=',NNTT
   DO J=1,NNTT
      DO I=1,NNW1
         WRITE(20,110) XNW(I,J),RNW(I,J)*DCOS(TNW(I,J)+PHI),RNW(I,J)*DSIN(TNW(I,J)+PHI)
      END DO !I=1,NNW1
   END DO !J=1,NNTT
END DO !K=2,NB
!-----------------------------------------------------------------------------------------------!
!    Deallocate Variables                                                                       !
!-----------------------------------------------------------------------------------------------!
IF (ISTRIP /= 1) DEALLOCATE(DSL)
!-----------------------------------------------------------------------------------------------!
!    Formats                                                                                    !
!-----------------------------------------------------------------------------------------------!
100 FORMAT(A,I4,A,I4,A,I4)
110 FORMAT(3(2X,E23.16))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE NOZZLEWAKEGRID
!-----------------------------------------------------------------------------------------------!
