!-----------------------------------------------------------------------------------------------!
!    Generate Blade Wake Grid                                                                   !
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
!-----------------------------------------------------------------------------------------------!
SUBROUTINE BLADEWAKEGRID
!-----------------------------------------------------------------------------------------------!
!    Created by: 14102013, J. Baltazar, version 1.0                                             !
!    Modified  : 21102013, J. Baltazar, 2013 version 1.0                                        !
!    Modified  : 02122014, J. Baltazar, 2014 version 1.2                                        !
!    Modified  : 07122014, J. Baltazar, 2014 version 1.3                                        !
!    Modified  : 07012015, J. Baltazar, 2015 version 1.0                                        !
!    Modified  : 25062015, J. Baltazar, 2015 version 1.2                                        !
!    Modified  : 02072015, J. Baltazar, 2015 version 1.3                                        !
!    Modified  : 05072017, J. Baltazar, 2017 version 1.0                                        !
!    Modified  : 18042019, J. Baltazar, 2019 version 1.0                                        !
!    Modified  : 09042020, J. Baltazar, 2020 version 1.0                                        !
!    Modified  : 17062020, J. Baltazar, 2020 version 1.1                                        !
!-----------------------------------------------------------------------------------------------!
!    Declarations                                                                               !
!-----------------------------------------------------------------------------------------------!
USE PROPANEL_MOD
IMPLICIT NONE
INTEGER :: I,J,J1,J2,K,IINI,IFIN
DOUBLE PRECISION :: XW0,PP,SS,P1,PHI,RR,RU,XC,TC,S0,S1,XN2,XN3
DOUBLE PRECISION :: CSI,FCSI,ETA,FETA,CSI1,CSI2,CSI3,CSI4
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: R,PTW0,PTW,DSL,DSL2,S,XTMW
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: PTMW
!-----------------------------------------------------------------------------------------------!
!    Counters                                                                                   !
!-----------------------------------------------------------------------------------------------!
NPW1=NPW+1
NRW1=NRW+1
!-----------------------------------------------------------------------------------------------!
ALLOCATE(XPW(NPW1,NRW1),YPW(NPW1,NRW1),ZPW(NPW1,NRW1),RPW(NPW1,NRW1),TPW(NPW1,NRW1))
XPW=0.D0
YPW=0.D0
ZPW=0.D0
RPW=0.D0
TPW=0.D0
!-----------------------------------------------------------------------------------------------!
!    Interpolate Pitch Data                                                                     !
!-----------------------------------------------------------------------------------------------!
ALLOCATE(R(NRW1),PTW0(NRW1),PTW(NRW1))
R   =0.D0
PTW0=0.D0
PTW =0.D0
!-----------------------------------------------------------------------------------------------!
IF (IGEOM /= 'WINGGEOM') R(1:NRW1)=0.5D0*(RP(1,JI:(JF+1))+RP(NCP1,JI:(JF+1))) !Rotor Case
IF (IGEOM == 'WINGGEOM') R(1:NRW1)=0.5D0*(YP(1,JI:(JF+1))+YP(NCP1,JI:(JF+1))) !Wing Case
IF (INTERPW == 0) THEN
   CALL LININT(NRIW,RIW,PTW0I,NRW1,R,PTW0)
   CALL LININT(NRIW,RIW,PTWI ,NRW1,R,PTW )
ELSEIF (INTERPW == 1) THEN
   CALL INTK1 (NRIW,RIW,PTW0I,NRW1,R,PTW0)
   CALL INTK1 (NRIW,RIW,PTWI ,NRW1,R,PTW )
ELSEIF (INTERPW == 2) THEN
   CALL SPLINT(NRIW,RIW,PTW0I,NRW1,R,PTW0)
   CALL SPLINT(NRIW,RIW,PTWI ,NRW1,R,PTW )
END IF !(INTERPW)
!-----------------------------------------------------------------------------------------------!
IF (IMODELPW == 5) THEN
   ALLOCATE(XTMW(NMW+1),PTMW(NRW1,NMW+1))
   XTMW=0.D0
   PTMW=0.D0
!  XTMW(1)=XPW(1,J)
   PTMW(:,1)=PTW0(:)
   DO I=1,NMW
      J1=  1+(I-1)*NRM
      J2=NRM+(I-1)*NRM
      XTMW(I+1)=XMW(I)
      CALL LININT(NRM,RMW(J1:J2),PMW(J1:J2),NRW1,R(:),PTMW(:,I+1))
   END DO !I=1,NMW
END IF !(IMODELPW  = 5)
!-----------------------------------------------------------------------------------------------!
!    Non-Dimensional Quantities                                                                 !
!-----------------------------------------------------------------------------------------------!
PTW0=2.D0*PTW0
PTW =2.D0*PTW
IF (IGEOM == 'WINGGEOM') PTW0=0.5D0*PTW0 !Wing Case
IF (IGEOM == 'WINGGEOM') PTW =0.5D0*PTW
!-----------------------------------------------------------------------------------------------!
!    Empirical Wake Model                                                                       !
!-----------------------------------------------------------------------------------------------!
IF (IMODELPW == 0) THEN
   R0=RPH
   R1=RMAX
!-----------------------------------------------------------------------------------------------!
!    Hoshino (1989)                                                                             !
!-----------------------------------------------------------------------------------------------!
ELSEIF (IMODELPW == 1) THEN
   IF (R1 == 0.D0) THEN
      IF (INTERPW == 0) THEN
         CALL LININT(NRI,RIW,PTW0I,1,0.7D0,P1)
      ELSEIF (INTERPW == 1) THEN
         CALL INTK1 (NRI,RIW,PTW0I,1,0.7D0,P1)
      ELSEIF (INTERPW == 2) THEN
         CALL SPLINT(NRI,RIW,PTW0I,1,0.7D0,P1)
      END IF !(INTERPW)
      SS=1.D0-ADVJ/P1
      R1=0.887D0-0.125D0*SS
   END IF !(R1 == 0.D0)
!-----------------------------------------------------------------------------------------------!
!    Kawakita (1992)                                                                            !
!-----------------------------------------------------------------------------------------------!
ELSEIF (IMODELPW == 2) THEN
   IF (INTERPW == 0) THEN
      CALL LININT(NRI,RIW,PTW0I,1,0.7D0,P1)
   ELSEIF (INTERPW == 1) THEN
      CALL INTK1 (NRI,RIW,PTW0I,1,0.7D0,P1)
   ELSEIF (INTERPW == 2) THEN
      CALL SPLINT(NRI,RIW,PTW0I,1,0.7D0,P1)
   END IF !(INTERPW)
   SS=1.D0-ADVJ/P1
   IF ((R1 == 0.D0).AND.(IN == 0)) R1=0.929D0-0.162D0*SS
   IF (IN == 0) PTW(:)=PTW0(:)*(0.874D0+0.067D0*SS)
   IF (IN /= 0) PTW(:)=PTW0(:)*(0.993D0-0.251D0*SS)
!-----------------------------------------------------------------------------------------------!
!    ProPan Empirical Model                                                                     !
!-----------------------------------------------------------------------------------------------!
ELSEIF (IMODELPW == 3) THEN
   IF (IN == 0) STOP 'WAKE MODEL ONLY FOR DUCTED PROPELLERS'
   IF (INTERPW == 0) THEN
      CALL LININT(NRI,RIW,PTW0I,1,0.7D0,P1)
   ELSEIF (INTERPW == 1) THEN
      CALL INTK1 (NRI,RIW,PTW0I,1,0.7D0,P1)
   ELSEIF (INTERPW == 2) THEN
      CALL SPLINT(NRI,RIW,PTW0I,1,0.7D0,P1)
   END IF !(INTERPW)
   SS=1.D0-ADVJ/P1
   IF (IN /= 0) PTW(:)=PTW0(:)*(A0+A1*SS)
END IF !(IMODELPW)
!-----------------------------------------------------------------------------------------------!
!    Stretching along the Streamwise                                                            !
!-----------------------------------------------------------------------------------------------!
ALLOCATE(DSL(NPW1))
DSL=0.D0
CALL STRET_CHOICE(NPW1,DSL,ST1PW,ST2PW,ST3PW,ST4PW,ITYPEPW)
!-----------------------------------------------------------------------------------------------!
!     Nozzle Downstream Stretching                                                              !
!-----------------------------------------------------------------------------------------------!
IF ((IN == 1).AND.(ISTRIP == 1)) THEN
   NNX=NNU+NC+NND
   NNXT=2*NNX+1
   NNTT=2*(NNT+1)
   ALLOCATE(S(NND+1),XN(NNU+NC:NNU+NC+NND+1,NNTT/2+1))
   XN2=0.5D0*(XP(1,NRP1+1)+XP(NCP1,NRP1+1))
   XN3=LD
   S0=DABS(XN3-XN2)/DABS(0.5D0*(XN3-XN2)*(DCOS(DFLOAT(NND-1)/DFLOAT(NND)*PI)+1.D0))/DFLOAT(NND)
   S1=DABS(XN3-XN2)/DABS(0.5D0*(XP(3,NRP1+1)+XP(NCP1-2,NRP1+1))- &
                         0.5D0*(XP(2,NRP1+1)+XP(NCP1-1,NRP1+1)))/DFLOAT(NND)
   CALL STRET2(S,S1,S0,NND+1)
   XN(NNU+NC+1,NNTT/2+1)=XN2
   DO I=1,NND
      XN(NNU+NC+I+1,NNTT/2+1)=XN2+(XN3-XN2)*S(I+1)
   END DO !I=1,NND
!-----------------------------------------------------------------------------------------------!
   I=2
   DO WHILE (NNU+NC+I <= NNXT/2+1)
      DSL(I)=DSL(I-1)+DABS(XN(NNU+NC+I,NNTT/2+1)-XN(NNU+NC+I-1,NNTT/2+1))
      I=I+1
   END DO
!  ST1PW=ST1PW*PI/180.D0
   ALLOCATE(DSL2(NPW1-I+2))
   CALL STRET_CHOICE(NPW1-I+2,DSL2,ST1PW,ST2PW,ST3PW,ST4PW,ITYPEPW)
   DO J=I,NPW1
!     DSL(J)=DSL(I-1)+DABS(XPWT-DSL(I-1))*(1.D0-DCOS((PI/2.D0-ST1PW)*DFLOAT(J-I+2)/ &
!            DFLOAT(NPW1-I+2)+ST1PW))/DCOS(ST1PW)
      DSL(J)=DSL(I-1)+DABS(XPWT-DSL(I-1))*DSL2(J-I+2)
   END DO !J=I,NPW1
   DSL(:)=DSL(:)/XPWT
   DEALLOCATE(S,XN,DSL2)
END IF !((IN == 1).AND.(ISTRIP == 1))
!-----------------------------------------------------------------------------------------------!
!    Loop on Radii                                                                              !
!-----------------------------------------------------------------------------------------------!
DO J=NRW1,1,-1
   IF (IGEOM /= 'WINGGEOM') THEN
      XPW(1,J)=0.5D0*(XP(1,J+JI-1)+XP(NCP1,J+JI-1))
      TPW(1,J)=0.5D0*(TP(1,J+JI-1)+TP(NCP1,J+JI-1))
      RPW(1,J)=0.5D0*(RP(1,J+JI-1)+RP(NCP1,J+JI-1))
      YPW(1,J)=RPW(1,J)*DCOS(TPW(1,J))
      ZPW(1,J)=RPW(1,J)*DSIN(TPW(1,J))
   ELSEIF (IGEOM == 'WINGGEOM') THEN
      XPW(1,J)=0.5D0*(XP(1,J+JI-1)+XP(NCP1,J+JI-1))
      YPW(1,J)=0.5D0*(YP(1,J+JI-1)+YP(NCP1,J+JI-1))
      ZPW(1,J)=0.5D0*(ZP(1,J+JI-1)+ZP(NCP1,J+JI-1))
   END IF !(IGEOM)
!-----------------------------------------------------------------------------------------------!
!    Loop on Axial Stations                                                                     !
!-----------------------------------------------------------------------------------------------!
   XW0=XPW(1,J)
   IF (IMODELPW == 5) XTMW(1)=XPW(1,J)
   DO I=2,NPW1
!-----------------------------------------------------------------------------------------------!
!    Axial Coordinate                                                                           !
!-----------------------------------------------------------------------------------------------!
      IF (XPW(I-1,J) < XPWW) THEN
         IF (IMODELPW <= 1) THEN
            CSI =(XPW(I-1,J)-XW0)/(XPWW-XW0)
            FCSI=DSQRT(CSI)+1.013D0*CSI-1.920D0*CSI**2+1.228D0*CSI**3-0.321D0*CSI**4 !Hoshino
            PP  =PTW0(J)+(PTW(J)-PTW0(J))*FCSI
         ELSEIF (IMODELPW == 2) THEN
            CSI =(XPW(I-1,J)-XW0)/(XPWW-XW0)
            FCSI=DSQRT(CSI)+1.013D0*CSI-1.920D0*CSI**2+1.228D0*CSI**3-0.321D0*CSI**4 !Hoshino
            ETA =(RPW(I-1,J)-R0)/(RPW(I-1,NRW1)-R0) !Kawakita
            IF (IN == 0) FETA=0.579D0+0.387D0*ETA+1.073D0*ETA**2+0.393D0*ETA**3-1.435D0*ETA**4
            IF (IN == 1) FETA=0.594D0+1.076D0*ETA-0.852D0*ETA**2+1.351D0*ETA**3-1.411D0*ETA**4
            PP  =(PTW0(J)+(PTW(J)-PTW0(J))*FCSI)*FETA
         ELSEIF (IMODELPW == 3) THEN
            CSI =(XPW(I-1,J)-XW0)/(XPWW-XW0)
            FCSI=DSQRT(CSI)+1.013D0*CSI-1.920D0*CSI**2+1.228D0*CSI**3-0.321D0*CSI**4 !Hoshino
            ETA =(RPW(I-1,J)-R0)/(RMAX-R0)
            IF (IN == 1) FETA=1.1063D0+3.7456*ETA-9.6472D0*ETA**2+14.88D0*ETA**3-8.6773D0*ETA**4
            PP  =(PTW0(J)+(PTW(J)-PTW0(J))*FCSI)*FETA
         ELSEIF (IMODELPW == 4) THEN
            CSI =DFLOAT(I-1)/DFLOAT(NTETA-1)
            CSI1=2.D0*CSI*CSI*CSI-3.D0*CSI*CSI+1.D0
            CSI2=CSI*CSI*CSI-2.D0*CSI*CSI+CSI
            CSI3=-2.D0*CSI*CSI*CSI+3.D0*CSI*CSI
            CSI4=CSI*CSI*CSI-CSI*CSI
            PP  =PTW0(J)*CSI1+DPTW0*CSI2+PTW(J)*CSI3+0.D0*CSI4
         ELSEIF (IMODELPW == 5) THEN
            IF (XPW(I-1,J) < XTMW(2)) THEN
               CSI =DFLOAT(I-1)/DFLOAT(NTETA-1)
               CSI1=2.D0*CSI*CSI*CSI-3.D0*CSI*CSI+1.D0
               CSI2=CSI*CSI*CSI-2.D0*CSI*CSI+CSI
               CSI3=-2.D0*CSI*CSI*CSI+3.D0*CSI*CSI
               CSI4=CSI*CSI*CSI-CSI*CSI
               PP  =PTMW(J,1)*CSI1+DPTW0*CSI2+PTMW(J,2)*CSI3+0.D0*CSI4
            ELSE
               CALL LININT(NMW+1,XTMW(:),PTMW(J,:),1,CSI,PP)
            END IF
         END IF !(IMODELPW)
      ELSE !(XPW(I-1,J) < XPWW)
         IF (IMODELPW <= 1) THEN
            PP =PTW(J)
         ELSEIF (IMODELPW == 2) THEN
            ETA=(RPW(I-1,J)-R0)/(RPW(I-1,NRW1)-R0) !Kawakita
            IF (IN == 0) FETA=0.579D0+0.387D0*ETA+1.073D0*ETA**2+0.393D0*ETA**3-1.435D0*ETA**4
            IF (IN == 1) FETA=0.594D0+1.076D0*ETA-0.852D0*ETA**2+1.351D0*ETA**3-1.411D0*ETA**4
            PP =PTW(J)*FETA
         ELSEIF (IMODELPW == 3) THEN
            ETA=(RPW(I-1,J)-R0)/(RMAX-R0)
            IF (IN == 1) FETA=1.1063D0+3.7456*ETA-9.6472D0*ETA**2+14.88D0*ETA**3-8.6773D0*ETA**4
            PP =PTW(J)*FETA
         ELSEIF (IMODELPW == 4) THEN
            PP =PTW(J)
         ELSEIF (IMODELPW == 5) THEN
            PP =PTMW(J,NMW+1)
         END IF !(IMODELPW)
      END IF !(XPW(I-1,J) < XPWW)
!-----------------------------------------------------------------------------------------------!
      IF (ISTEADY == 0) THEN
         XC=XPW(1,J)+DSL(I)*XPWT
         IF (IGEOM /= 'WINGGEOM') TC=TPW(I-1,J)+2.D0*PI*DABS(XC-XPW(I-1,J))/PP
         IF (IGEOM == 'WINGGEOM') TC=ZPW(I-1,J)-DABS(XC-XPW(I-1,J))*DTAND(PP)
      ELSEIF (ISTEADY == 1) THEN
         XC=XPW(I-1,J)+PP/DFLOAT(NTETA)
         IF (IGEOM /= 'WINGGEOM') TC=TPW(I-1,J)+2.D0*PI/DFLOAT(NTETA)
         IF (IGEOM == 'WINGGEOM') STOP 'NOT PREPARED FOR ISTEADY=1'
      END IF !(ISTEADY)
!-----------------------------------------------------------------------------------------------!
!    Radial Coordinate With Nozzle                                                              !
!-----------------------------------------------------------------------------------------------!
      IF ((IN == 1).AND.(ISTRIP == 1)) THEN
         IF (J == NRW1) THEN
            IF (XC < LD) THEN
               CALL NOZZLEDEF('INNER',XC,RR)
            ELSE !(XC < LD)
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
                  RR  =0.5D0*(YOL(NRNI)+YIL(NRNI))
                  RR  =RR*LD*2.D0+RMAX+CR
               END IF !(ICONTRNW == 1)
            END IF !(XC < LD)
         ELSE !(J == NRW1)
            IF (XPW(I-1,J) < XPWW) THEN
               CSI =(XC-XW0)/(XPWW-XW0)
               FCSI=DSQRT(CSI)+1.013D0*CSI-1.920D0*CSI**2+1.228D0*CSI**3-0.321D0*CSI**4
               PHI=DFLOAT(J-1)*(PI-ALPHAT-ALPHAH)/DFLOAT(NRW)+ALPHAH
               RU =((R0*DCOS(ALPHAT)+RPW(I,NRW1)*DCOS(ALPHAH))- &
                    (RPW(I,NRW1)-R0)*DCOS(PHI))/(DCOS(ALPHAH)+DCOS(ALPHAT))
               RR =R(J)+(RU-R(J))*FCSI
            ELSE !(XPW(I-1,J) < XPWW)
               PHI=DFLOAT(J-1)*(PI-ALPHAT-ALPHAH)/DFLOAT(NRW)+ALPHAH
               RU =((R0*DCOS(ALPHAT)+RPW(I,NRW1)*DCOS(ALPHAH))- &
                    (RPW(I,NRW1)-R0)*DCOS(PHI))/(DCOS(ALPHAH)+DCOS(ALPHAT))
               RR =RU
            END IF !(XPW(I-1,J) < XPWW)
         END IF !(J == NRW1)
!-----------------------------------------------------------------------------------------------!
!    Radial Coordinate Without Nozzle                                                           !
!-----------------------------------------------------------------------------------------------!
      ELSE !((IN == 1).AND.(ISTRIP == 1))
         IF (XPW(I-1,J) < XPWW) THEN
            IF (IMODELPW < 4) THEN
               CSI =(XC-XW0)/(XPWW-XW0)
               FCSI=DSQRT(CSI)+1.013D0*CSI-1.920D0*CSI**2+1.228D0*CSI**3-0.321D0*CSI**4
               PHI=DFLOAT(J+JI-2)*(PI-ALPHAT-ALPHAH)/DFLOAT(NRP)+ALPHAH
               RU =((R0*DCOS(ALPHAT)+R1*DCOS(ALPHAH))- &
                    (R1-R0)*DCOS(PHI))/(DCOS(ALPHAH)+DCOS(ALPHAT))
               RR =R(J)+(RU-R(J))*FCSI
            ELSEIF (IMODELPW == 4) THEN
               PHI=DFLOAT(J+JI-2)*(PI-ALPHAT-ALPHAH)/DFLOAT(NRP)+ALPHAH
               RU =((R0*DCOS(ALPHAT)+R1*DCOS(ALPHAH))- &
                    (R1-R0)*DCOS(PHI))/(DCOS(ALPHAH)+DCOS(ALPHAT))
               RR =RU
            ELSEIF (IMODELPW == 5) THEN
               PHI=DFLOAT(J+JI-2)*(PI-ALPHAT-ALPHAH)/DFLOAT(NRP)+ALPHAH
               RU =((R0*DCOS(ALPHAT)+R1*DCOS(ALPHAH))- &
                    (R1-R0)*DCOS(PHI))/(DCOS(ALPHAH)+DCOS(ALPHAT))
               RR =RU
            END IF
         ELSE !(XPW(I-1,J) < XPWW)
            PHI=DFLOAT(J+JI-2)*(PI-ALPHAT-ALPHAH)/DFLOAT(NRP)+ALPHAH
            RU =((R0*DCOS(ALPHAT)+R1*DCOS(ALPHAH))- &
                 (R1-R0)*DCOS(PHI))/(DCOS(ALPHAH)+DCOS(ALPHAT))
            RR =RU
         END IF !(XPW(I-1,J) < XPWW)
      END IF !((IN == 1).AND.(ISTRIP == 1))
!-----------------------------------------------------------------------------------------------!
      IF (IGEOM /= 'WINGGEOM') THEN
         XPW(I,J)=XC
         TPW(I,J)=TC
         RPW(I,J)=RR
         YPW(I,J)=RPW(I,J)*DCOS(TPW(I,J))
         ZPW(I,J)=RPW(I,J)*DSIN(TPW(I,J))
      ELSEIF (IGEOM == 'WINGGEOM') THEN
         XPW(I,J)=XC
         YPW(I,J)=RR
         ZPW(I,J)=TC
      END IF !(IGEOM)
   END DO !I=2,NPW1
END DO !J=NRW1,1,-1
!-----------------------------------------------------------------------------------------------!
!    Radial Coordinate Correction due to Thick Nozzle at T.E.                                   !
!-----------------------------------------------------------------------------------------------!
IF (INTECORR == 1) THEN
   J=NRW1
   I=1
   DO WHILE ((XPW(I,J) < XI).AND.(I < NPW1))
      I=I+1
   END DO !((XPW(I,J) < XI).AND.(I < NPW1))
   I=I-1
   IINI=I
   DO WHILE ((XPW(I,J) < XF).AND.(I < NPW1))
      I=I+1
   END DO !((XPW(I,J) < XF).AND.(I < NPW1))
   IFIN=I
   DO J=NRW,2,-1
      DO I=IINI,IFIN
         CSI=(XPW(I,J)-XPW(IINI,J))/(XPW(IFIN,J)-XPW(IINI,J))
         RPW(I,J)=(1.D0-3.D0*CSI**2+2.D0*CSI**3)*RPW(IINI,J)+ &
                  (     3.D0*CSI**2-2.D0*CSI**3)*RPW(IFIN,J)+ &
                  ( CSI-2.D0*CSI**2+1.D0*CSI**3)*(RPW(IINI+1,J)-RPW(IINI,J))/ &
                       (XPW(IINI+1,J)-XPW(IINI,J))*(XPW(IFIN,J)-XPW(IINI,J))+ &
                  (    -1.D0*CSI**2+1.D0*CSI**3)*(RPW(IFIN,J)-RPW(IFIN-1,J))/ &
                       (XPW(IFIN,J)-XPW(IFIN-1,J))*(XPW(IFIN,J)-XPW(IINI,J))
         YPW(I,J)=RPW(I,J)*DCOS(TPW(I,J))
         ZPW(I,J)=RPW(I,J)*DSIN(TPW(I,J))
      END DO !I=IINI,IFIN
   END DO !J=2,NRW
END IF !(INTECORR == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub Matching With Blade Wake Root                                                          !
!-----------------------------------------------------------------------------------------------!
IF (IHCORR == 1) THEN
   DO I=2,NPW1
      IF (INTERH == 0) THEN
         CALL LININT(NHI,XHI,RHI,1,XPW(I,1),RPW(I,1))
      ELSEIF (INTERH == 1) THEN
         CALL INTK1 (NHI,XHI,RHI,1,XPW(I,1),RPW(I,1))
      ELSEIF (INTERH == 2) THEN
         CALL SPLINT(NHI,XHI,RHI,1,XPW(I,1),RPW(I,1))
      END IF !(INTERH)
      YPW(I,1)=RPW(I,1)*DCOS(TPW(I,1))
      ZPW(I,1)=RPW(I,1)*DSIN(TPW(I,1))
   END DO !I=1,NCP1
END IF !(IHCORR == 1)
!-----------------------------------------------------------------------------------------------!
!    Write Blade Wake Grid in Tecplot Format                                                    !
!-----------------------------------------------------------------------------------------------!
K=1
WRITE(20,100) ' ZONE T="BLADE WAKE',K,'" F=POINT, I=',NPW1,' J=',NRW1
DO J=1,NRW1
   DO I=1,NPW1
      WRITE(20,110) XPW(I,J),YPW(I,J),ZPW(I,J)
   END DO !I=1,NPW1
END DO !J=1,NRW1
DO K=2,NB
   PHI=DFLOAT(K-1)/DFLOAT(NB)*2.D0*PI
   WRITE(20,100) ' ZONE T="BLADE WAKE',K,'" F=POINT, I=',NPW1,' J=',NRW1
   DO J=1,NRW1
      DO I=1,NPW1
         IF (IGEOM /= 'WINGGEOM') WRITE(20,110) XPW(I,J),RPW(I,J)*DCOS(TPW(I,J)+PHI), &
                                                                     RPW(I,J)*DSIN(TPW(I,J)+PHI)
         IF (IGEOM == 'WINGGEOM') WRITE(20,110) XPW(I,J),-YPW(I,J),ZPW(I,J)
      END DO !I=1,NPW1
   END DO !J=1,NRW1
END DO !K=2,NB
!-----------------------------------------------------------------------------------------------!
!    Deallocate Variables                                                                       !
!-----------------------------------------------------------------------------------------------!
DEALLOCATE(R,PTW0,PTW,DSL)
!-----------------------------------------------------------------------------------------------!
!    Formats                                                                                    !
!-----------------------------------------------------------------------------------------------!
100 FORMAT(A,I4,A,I5,A,I4)
110 FORMAT(3(2X,E23.16))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE BLADEWAKEGRID
!-----------------------------------------------------------------------------------------------!
