!-----------------------------------------------------------------------------------------------!
!    Cavity pressure recovery                                                                   !
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
DOUBLE PRECISION FUNCTION SIGMAREC(I,J,TT)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST, April 2009                                                   !
!    Modified  : 22042016, J. Baltazar, version 1.1                                             !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,TT,I2,IM1,IP1
DOUBLE PRECISION :: SA,SB,S1(NCP)
DOUBLE PRECISION :: LAMBDA,ST,CSI
DOUBLE PRECISION :: H,DSIGMA,CPL,DCPL
!-----------------------------------------------------------------------------------------------!
LAMBDA=0.925D0
!-----------------------------------------------------------------------------------------------!
!    Pressure Side                                                                              !
!-----------------------------------------------------------------------------------------------!
IF ((I <= IDP(J,TT)).AND.(I >= IRP(J,TT))) THEN
   S1(NC)=0.D0
   DO I2=NC-1,1,-1
      IP1=I2+1
      SA=DSQRT(AT1XP(IP1,J)**2+AT1YP(IP1,J)**2+AT1ZP(IP1,J)**2)
      SB=DSQRT(AT1XP(I2 ,J)**2+AT1YP(I2 ,J)**2+AT1ZP(I2 ,J)**2)
      S1(I2)=S1(IP1)+SA+SB
   END DO !I2=NC-1,1,-1
!-----------------------------------------------------------------------------------------------!
   ST =S1(IDP(J,TT))+LAMBDA*(S1(MAX(IRP(J,TT)-1,1))-S1(IDP(J,TT)))
   CSI=(S1(I)-ST)/(S1(MAX(IRP(J,TT)-1,1))-ST)
!-----------------------------------------------------------------------------------------------!
   if (irp(j,tt) == 1) then
      sigmarec=sigma
   ELSEIF (CSI <= 0.D0) THEN
      SIGMAREC=SIGMA
   ELSE
      DSIGMA=0.D0
      CPL =CPNP(MAX(IRP(J,TT)-1,1),J,TT)
      DCPL=0.D0
      SIGMAREC=SIGMA*H(1,CSI)+DSIGMA*H(2,CSI)-CPL*H(3,CSI)-DCPL*H(4,CSI)
   END IF
!-----------------------------------------------------------------------------------------------!
!    Suction Side                                                                               !
!-----------------------------------------------------------------------------------------------!
ELSEIF ((I >= IDS(J,TT)).AND.(I <= IRS(J,TT))) THEN
   S1(NC1)=0.D0
   DO I2=NC1+1,NCP
      IM1=I2-1
      SA =DSQRT(AT1XP(IM1,J)**2+AT1YP(IM1,J)**2+AT1ZP(IM1,J)**2)
      SB =DSQRT(AT1XP(I2 ,J)**2+AT1YP(I2 ,J)**2+AT1ZP(I2 ,J)**2)
      S1(I2)=S1(IM1)+SA+SB
   END DO !I2=NC1+1,NCP
!-----------------------------------------------------------------------------------------------!
   ST =S1(IDS(J,TT))+LAMBDA*(S1(MIN(IRS(J,TT)+1,NCP))-S1(IDS(J,TT)))
   CSI=(S1(I)-ST)/(S1(MIN(IRS(J,TT)+1,NCP))-ST)
!-----------------------------------------------------------------------------------------------!
   if (irs(j,tt) == ncp) then
      sigmarec=sigma
   ELSEIF (CSI <= 0.D0) THEN
      SIGMAREC=SIGMA
   ELSE
      DSIGMA=0.D0
      CPL =CPNP(MIN(IRS(J,TT)+1,NCP),J,TT)
      DCPL=0.D0
      SIGMAREC=SIGMA*H(1,CSI)+DSIGMA*H(2,CSI)-CPL*H(3,CSI)-DCPL*H(4,CSI)
   END IF
END IF
!-----------------------------------------------------------------------------------------------!
END FUNCTION SIGMAREC
!===============================================================================================!
!    Hermite Functions                                                                          !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION FUNCTION H(I,CSI)
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
INTEGER :: I
DOUBLE PRECISION :: CSI
!-----------------------------------------------------------------------------------------------!
SELECT CASE (I)
   CASE(1)
      H=(1.D0+2.D0*CSI)*(CSI-1.D0)*(CSI-1.D0)
   CASE(2)
      H=CSI*(CSI-1.D0)*(CSI-1.D0)
   CASE(3)
      H=(1.D0-2.D0*(CSI-1.D0))*CSI*CSI
   CASE(4)
      H=(CSI-1.D0)*CSI*CSI
END SELECT !(I)
!-----------------------------------------------------------------------------------------------!
END FUNCTION H
!-----------------------------------------------------------------------------------------------!
