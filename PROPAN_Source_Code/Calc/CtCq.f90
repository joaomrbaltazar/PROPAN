!-----------------------------------------------------------------------------------------------!
!    Calculation of Rotor Forces                                                                !
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
SUBROUTINE CTCQ(JJ,TT)
!-----------------------------------------------------------------------------------------------!
!    Created by: J. Baltazar, IST                                                               !
!    Modified  : 03122013, J. Baltazar, version 1.0                                             !
!    Modified  : 26052014, J. Baltazar, Positive blade forces for turbine case                  !
!    Modified  : 06112014, J. Baltazar, Hub forces only for IH=1                                !
!    Modified  : 07012015, J. Baltazar, 2015 version 1.0, Wing Case                             !
!    Modified  : 23062020, J. Baltazar, 2020 version 1.2, Dynamic inflow                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,KB,JJ,TT
!-----------------------------------------------------------------------------------------------!
!    Integration on the Blade                                                                   !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
!-----------------------------------------------------------------------------------------------!
!    X Direction                                                                                !
!-----------------------------------------------------------------------------------------------!
   KB=1 !First Blade
   WORK1=0.D0
   WORK2=0.D0
   DO J=1,NRP-1
      DO I=1,NCP
         WORK1=WORK1+CPNP(I,J,TT)*UNXP0(I,J)*AP0(I,J)
         WORK2=WORK2+CPNP(I,J,TT)*(UNYP0(I,J)*ZP0(I,J)-UNZP0(I,J)*YP0(I,J))*AP0(I,J)
      END DO !I=1,NCP
   END DO !J=1,NRP-1
!-----------------------------------------------------------------------------------------------!
!    Include last blade strip                                                                   !
!-----------------------------------------------------------------------------------------------!
   IF (ISTRIP /= 1) THEN
      J=NRP
      DO I=1,NCP
         WORK1=WORK1+CPNP(I,J,TT)*UNXP0(I,J)*AP0(I,J)
         WORK2=WORK2+CPNP(I,J,TT)*(UNYP0(I,J)*ZP0(I,J)-UNZP0(I,J)*YP0(I,J))*AP0(I,J)
      END DO !I=1,NCP
   END IF !(ISTRIP /= 1)
!-----------------------------------------------------------------------------------------------!
!    Force Coefficients                                                                         !
!-----------------------------------------------------------------------------------------------!
   IF (IROTOR == -1) THEN
      CFXP(TT,KB)=-WORK1/4.D0
      CMXP(TT,KB)=0.D0
   ELSEIF (IROTOR == 0) THEN
      CFXP(TT,KB)=WORK1/8.D0
      CMXP(TT,KB)=WORK2/16.D0
   ELSEIF (IROTOR == 1) THEN
      CFXP(TT,KB)=-WORK1*UU(JJ)**2/PI**3
      CMXP(TT,KB)=-WORK2*UU(JJ)**3/PI**3
   END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
   DO KB=2,NB
      IF (TT < INT(P)*(KB-1)) THEN
         CFXP(TT,KB)=CFXP(0,1)
         CMXP(TT,KB)=CMXP(0,1)
      ELSE
         CFXP(TT,KB)=CFXP((TT-INT(P)*(KB-1)),1)
         CMXP(TT,KB)=CMXP((TT-INT(P)*(KB-1)),1)
      END IF !(TT < P*(KB-1))
   END DO !KB=2,NB
!-----------------------------------------------------------------------------------------------!
!    Y Direction                                                                                !
!-----------------------------------------------------------------------------------------------!
   KB=1 !First Blade
   WORK1=0.D0
   WORK2=0.D0
   DO J=1,NRP-1
      DO I=1,NCP
         WORK1=WORK1+CPNP(I,J,TT)*UNYP0(I,J)*AP0(I,J)
         WORK2=WORK2+CPNP(I,J,TT)*(UNXP0(I,J)*ZP0(I,J)-UNZP0(I,J)*XP0(I,J))*AP0(I,J)
      END DO !I=1,NCP
   END DO !J=1,NRP-1
!-----------------------------------------------------------------------------------------------!
!    Include last blade strip                                                                   !
!-----------------------------------------------------------------------------------------------!
   IF (ISTRIP /= 1) THEN
      J=NRP
      DO I=1,NCP
         WORK1=WORK1+CPNP(I,J,TT)*UNYP0(I,J)*AP0(I,J)
         WORK2=WORK2+CPNP(I,J,TT)*(UNXP0(I,J)*ZP0(I,J)-UNZP0(I,J)*XP0(I,J))*AP0(I,J)
      END DO !I=1,NCP
   END IF !(ISTRIP /= 1)
!-----------------------------------------------------------------------------------------------!
!    Force Coefficients                                                                         !
!-----------------------------------------------------------------------------------------------!
   IF (IROTOR == -1) THEN
      CFYP(TT,KB)=-WORK1/4.D0
      CMYP(TT,KB)=0.D0
   ELSEIF (IROTOR == 0) THEN
      CFYP(TT,KB)=WORK1/8.D0
      CMYP(TT,KB)=WORK2/16.D0
   ELSEIF (IROTOR == 1) THEN
      CFYP(TT,KB)=-WORK1*UU(JJ)**2/PI**3
      CMYP(TT,KB)=-WORK2*UU(JJ)**3/PI**3
   END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
   DO KB=2,NB
      IF (TT < INT(P)*(KB-1)) THEN
         CFYP(TT,KB)=CFYP(0,1)
         CMYP(TT,KB)=CMYP(0,1)
      ELSE
         CFYP(TT,KB)=CFYP((TT-INT(P)*(KB-1)),1)
         CMYP(TT,KB)=CMYP((TT-INT(P)*(KB-1)),1)
      END IF !(TT < P*(KB-1))
   END DO !KB=2,NB
!-----------------------------------------------------------------------------------------------!
!    Z Direction                                                                                !
!-----------------------------------------------------------------------------------------------!
   KB=1 !First Blade
   WORK1=0.D0
   WORK2=0.D0
   DO J=1,NRP-1
      DO I=1,NCP
         WORK1=WORK1+CPNP(I,J,TT)*UNZP0(I,J)*AP0(I,J)
         WORK2=WORK2+CPNP(I,J,TT)*(UNYP0(I,J)*XP0(I,J)-UNXP0(I,J)*YP0(I,J))*AP0(I,J)
      END DO !I=1,NCP
   END DO !J=1,NRP-1
!-----------------------------------------------------------------------------------------------!
!    Include last blade strip                                                                   !
!-----------------------------------------------------------------------------------------------!
   IF (ISTRIP /= 1) THEN
      J=NRP
      DO I=1,NCP
         WORK1=WORK1+CPNP(I,J,TT)*UNZP0(I,J)*AP0(I,J)
         WORK2=WORK2+CPNP(I,J,TT)*(UNYP0(I,J)*XP0(I,J)-UNXP0(I,J)*YP0(I,J))*AP0(I,J)
      END DO !I=1,NCP
   END IF !(ISTRIP /= 1)
!-----------------------------------------------------------------------------------------------!
!    Force Coefficients                                                                         !
!-----------------------------------------------------------------------------------------------!
   IF (IROTOR == -1) THEN
      CFZP(TT,KB)=-WORK1/4.D0
      CMZP(TT,KB)=0.D0
   ELSEIF (IROTOR == 0) THEN
      CFZP(TT,KB)=WORK1/8.D0
      CMZP(TT,KB)=WORK2/16.D0
   ELSEIF (IROTOR == 1) THEN
      CFZP(TT,KB)=-WORK1*UU(JJ)**2/PI**3
      CMZP(TT,KB)=-WORK2*UU(JJ)**3/PI**3
   END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
   DO KB=2,NB
      IF (TT < INT(P)*(KB-1)) THEN
         CFZP(TT,KB)=CFZP(0,1)
         CMZP(TT,KB)=CMZP(0,1)
      ELSE
         CFZP(TT,KB)=CFZP((TT-INT(P)*(KB-1)),1)
         CMZP(TT,KB)=CMZP((TT-INT(P)*(KB-1)),1)
      END IF !(TT < P*(KB-1))
   END DO !KB=2,NB
!-----------------------------------------------------------------------------------------------!
!    Rotor Coefficients                                                                         !
!-----------------------------------------------------------------------------------------------!
   CTXP(TT)=0.D0
   CTYP(TT)=0.D0
   CTZP(TT)=0.D0
   CQXP(TT)=0.D0
   CQYP(TT)=0.D0
   CQZP(TT)=0.D0
   DO KB=1,NB
      IF (IROTOR == -1) THEN
         CTXP(TT)=CTXP(TT)+CFZP(TT,KB)*DCOSD(UU(JJ))-CFXP(TT,KB)*DSIND(UU(JJ))
         CQXP(TT)=CQXP(TT)+CFZP(TT,KB)*DSIND(UU(JJ))+CFXP(TT,KB)*DCOSD(UU(JJ))
      ELSE !(IROTOR)
         CTXP(TT)=CTXP(TT)+CFXP(TT,KB)
         CTYP(TT)=CTXP(TT)+CFYP(TT,KB)
         CTZP(TT)=CTZP(TT)+CFZP(TT,KB)
         CQXP(TT)=CQXP(TT)+CMXP(TT,KB)
         CQYP(TT)=CQYP(TT)+CMYP(TT,KB)
         CQZP(TT)=CQZP(TT)+CMZP(TT,KB)
      END IF !(IROTOR)
   END DO !KB=1,NB
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Integration on the Nozzle                                                                  !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
!-----------------------------------------------------------------------------------------------!
!    X Direction                                                                                !
!-----------------------------------------------------------------------------------------------!
   KB=1 !First Sector
   WORK1=0.D0
   WORK2=0.D0
   DO J=1,NNTP
      DO I=1,NNXT1
         WORK1=WORK1+CPNN(I,J,TT)*UNXN0(I,J)*AN0(I,J)
         WORK2=WORK2+CPNN(I,J,TT)*(UNYN0(I,J)*ZN0(I,J)-UNZN0(I,J)*YN0(I,J))*AN0(I,J)
      END DO !I=1,NNXT1
   END DO !J=1,NNTP
!-----------------------------------------------------------------------------------------------!
!    Force Coefficients                                                                         !
!-----------------------------------------------------------------------------------------------!
   IF (IROTOR == 0) THEN
      CFXN(TT,KB)=WORK1/8.D0
      CMXN(TT,KB)=WORK2/16.D0
   ELSEIF (IROTOR == 1) THEN
      CFXN(TT,KB)=WORK1*UU(JJ)**2/PI**3
      CMXN(TT,KB)=WORK2*UU(JJ)**3/PI**3
   END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
   DO KB=2,NB
      IF (TT < INT(P)*(KB-1)) THEN
         CFXN(TT,KB)=CFXN(0,1)
         CMXN(TT,KB)=CMXN(0,1)
      ELSE
         CFXN(TT,KB)=CFXN((TT-INT(P)*(KB-1)),1)
         CMXN(TT,KB)=CMXN((TT-INT(P)*(KB-1)),1)
      END IF !(TT < P*(KB-1))
   END DO !KB=2,NB
!-----------------------------------------------------------------------------------------------!
!    Y Direction                                                                                !
!-----------------------------------------------------------------------------------------------!
   KB=1 !First Sector
   WORK1=0.D0
   WORK2=0.D0
   DO J=1,NNTP
      DO I=1,NNXT1
         WORK1=WORK1+CPNN(I,J,TT)*UNYN0(I,J)*AN0(I,J)
         WORK2=WORK2+CPNN(I,J,TT)*(UNXN0(I,J)*ZN0(I,J)-UNZN0(I,J)*XN0(I,J))*AN0(I,J)
      END DO !I=1,NNXT1
   END DO !J=1,NNTP
!-----------------------------------------------------------------------------------------------!
!    Force Coefficients                                                                         !
!-----------------------------------------------------------------------------------------------!
   IF (IROTOR == 0) THEN
      CFYN(TT,KB)=WORK1/8.D0
      CMYN(TT,KB)=WORK2/16.D0
   ELSEIF (IROTOR == 1) THEN
      CFYN(TT,KB)=WORK1*UU(JJ)**2/PI**3
      CMYN(TT,KB)=WORK2*UU(JJ)**3/PI**3
   END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
   DO KB=2,NB
      IF (TT < INT(P)*(KB-1)) THEN
         CFYN(TT,KB)=CFYN(0,1)
         CMYN(TT,KB)=CMYN(0,1)
      ELSE
         CFYN(TT,KB)=CFYN((TT-INT(P)*(KB-1)),1)
         CMYN(TT,KB)=CMYN((TT-INT(P)*(KB-1)),1)
      END IF !(TT < P*(KB-1))
   END DO !KB=2,NB
!-----------------------------------------------------------------------------------------------!
!    Z Direction                                                                                !
!-----------------------------------------------------------------------------------------------!
   KB=1 !First Sector
   WORK1=0.D0
   WORK2=0.D0
   DO J=1,NNTP
      DO I=1,NNXT1
         WORK1=WORK1+CPNN(I,J,TT)*UNZN0(I,J)*AN0(I,J)
         WORK2=WORK2+CPNN(I,J,TT)*(UNYN0(I,J)*XN0(I,J)-UNXN0(I,J)*YN0(I,J))*AN0(I,J)
      END DO !I=1,NNXT1
   END DO !J=1,NNTP
!-----------------------------------------------------------------------------------------------!
!    Force Coefficients                                                                         !
!-----------------------------------------------------------------------------------------------!
   IF (IROTOR == 0) THEN
      CFZN(TT,KB)=WORK1/8.D0
      CMZN(TT,KB)=WORK2/16.D0
   ELSEIF (IROTOR == 1) THEN
      CFZN(TT,KB)=WORK1*UU(JJ)**2/PI**3
      CMZN(TT,KB)=WORK2*UU(JJ)**3/PI**3
   END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
   DO KB=2,NB
      IF (TT < INT(P)*(KB-1)) THEN
         CFZN(TT,KB)=CFZN(0,1)
         CMZN(TT,KB)=CMZN(0,1)
      ELSE
         CFZN(TT,KB)=CFZN((TT-INT(P)*(KB-1)),1)
         CMZN(TT,KB)=CMZN((TT-INT(P)*(KB-1)),1)
      END IF !(TT < P*(KB-1))
   END DO !KB=2,NB
!-----------------------------------------------------------------------------------------------!
!    Rotor Coefficients                                                                         !
!-----------------------------------------------------------------------------------------------!
   CTXN(TT)=0.D0
   CTYN(TT)=0.D0
   CTZN(TT)=0.D0
   CQXN(TT)=0.D0
   CQYN(TT)=0.D0
   CQZN(TT)=0.D0
   DO KB=1,NB
      CTXN(TT)=CTXN(TT)+CFXN(TT,KB)
      CTYN(TT)=CTYN(TT)+CFYN(TT,KB)
      CTZN(TT)=CTZN(TT)+CFZN(TT,KB)
      CQXN(TT)=CQXN(TT)+CMXN(TT,KB)
      CQYN(TT)=CQYN(TT)+CMYN(TT,KB)
      CQZN(TT)=CQZN(TT)+CMZN(TT,KB)
   END DO !KB=1,NB
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Integration on the Hub                                                                     !
!-----------------------------------------------------------------------------------------------!
IF (IH == 1) THEN
!-----------------------------------------------------------------------------------------------!
!    X Direction                                                                                !
!-----------------------------------------------------------------------------------------------!
   KB=1 !First Sector
   WORK1=0.D0
   WORK2=0.D0
   DO J=1,NHTP
      DO I=1,NHX
         WORK1=WORK1+CPNH(I,J,TT)*UNXH0(I,J)*AH0(I,J)
         WORK2=WORK2+CPNH(I,J,TT)*(UNYH0(I,J)*ZH0(I,J)-UNZH0(I,J)*YH0(I,J))*AH0(I,J)
      END DO !I=1,NHX
   END DO !J=1,NHTP
!-----------------------------------------------------------------------------------------------!
!    Force Coefficients                                                                         !
!-----------------------------------------------------------------------------------------------!
   IF (IROTOR == 0) THEN
      CFXH(TT,KB)=WORK1/8.D0
      CMXH(TT,KB)=WORK2/16.D0
   ELSEIF (IROTOR == 1) THEN
      CFXH(TT,KB)=WORK1*UU(JJ)**2/PI**3
      CMXH(TT,KB)=WORK2*UU(JJ)**3/PI**3
   END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
   DO KB=2,NB
      IF (TT < INT(P)*(KB-1)) THEN
         CFXH(TT,KB)=CFXH(0,1)
         CMXH(TT,KB)=CMXH(0,1)
      ELSE
         CFXH(TT,KB)=CFXH((TT-INT(P)*(KB-1)),1)
         CMXH(TT,KB)=CMXH((TT-INT(P)*(KB-1)),1)
      END IF !(TT < P*(KB-1))
   END DO !KB=2,NB
!-----------------------------------------------------------------------------------------------!
!    Y Direction                                                                                !
!-----------------------------------------------------------------------------------------------!
   KB=1 !First Sector
   WORK1=0.D0
   WORK2=0.D0
   DO J=1,NHTP
      DO I=1,NHX
         WORK1=WORK1+CPNH(I,J,TT)*UNYH0(I,J)*AH0(I,J)
         WORK2=WORK2+CPNH(I,J,TT)*(UNXH0(I,J)*ZH0(I,J)-UNZH0(I,J)*XH0(I,J))*AH0(I,J)
      END DO !I=1,NHX
   END DO !J=1,NHTP
!-----------------------------------------------------------------------------------------------!
!    Force Coefficients                                                                         !
!-----------------------------------------------------------------------------------------------!
   IF (IROTOR == 0) THEN
      CFYH(TT,KB)=WORK1/8.D0
      CMYH(TT,KB)=WORK2/16.D0
   ELSEIF (IROTOR == 1) THEN
      CFYH(TT,KB)=WORK1*UU(JJ)**2/PI**3
      CMYH(TT,KB)=WORK2*UU(JJ)**3/PI**3
   END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
   DO KB=2,NB
      IF (TT < INT(P)*(KB-1)) THEN
         CFYH(TT,KB)=CFYH(0,1)
         CMYH(TT,KB)=CMYH(0,1)
      ELSE
         CFYH(TT,KB)=CFYH((TT-INT(P)*(KB-1)),1)
         CMYH(TT,KB)=CMYH((TT-INT(P)*(KB-1)),1)
      END IF !(TT < P*(KB-1))
   END DO !KB=2,NB
!-----------------------------------------------------------------------------------------------!
!    Z Direction                                                                                !
!-----------------------------------------------------------------------------------------------!
   KB=1 !First Sector
   WORK1=0.D0
   WORK2=0.D0
   DO J=1,NHTP
      DO I=1,NHX
         WORK1=WORK1+CPNH(I,J,TT)*UNZH0(I,J)*AH0(I,J)
         WORK2=WORK2+CPNH(I,J,TT)*(UNYH0(I,J)*XH0(I,J)-UNXH0(I,J)*YH0(I,J))*AH0(I,J)
      END DO !I=1,NHX
   END DO !J=1,NHTP
!-----------------------------------------------------------------------------------------------!
!    Force Coefficients                                                                         !
!-----------------------------------------------------------------------------------------------!
   IF (IROTOR == 0) THEN
      CFZH(TT,KB)=WORK1/8.D0
      CMZH(TT,KB)=WORK2/16.D0
   ELSEIF (IROTOR == 1) THEN
      CFZH(TT,KB)=WORK1*UU(JJ)**2/PI**3
      CMZH(TT,KB)=WORK2*UU(JJ)**3/PI**3
   END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
   DO KB=2,NB
      IF (TT < INT(P)*(KB-1)) THEN
         CFZH(TT,KB)=CFZH(0,1)
         CMZH(TT,KB)=CMZH(0,1)
      ELSE
         CFZH(TT,KB)=CFZH((TT-INT(P)*(KB-1)),1)
         CMZH(TT,KB)=CMZH((TT-INT(P)*(KB-1)),1)
      END IF !(TT < P*(KB-1))
   END DO !KB=2,NB
!-----------------------------------------------------------------------------------------------!
!    Rotor Coefficients                                                                         !
!-----------------------------------------------------------------------------------------------!
   CTXH(TT)=0.D0
   CTYH(TT)=0.D0
   CTZH(TT)=0.D0
   CQXH(TT)=0.D0
   CQYH(TT)=0.D0
   CQZH(TT)=0.D0
   DO KB=1,NB
      CTXH(TT)=CTXH(TT)+CFXH(TT,KB)
      CTYH(TT)=CTYH(TT)+CFYH(TT,KB)
      CTZH(TT)=CTZH(TT)+CFZH(TT,KB)
      CQXH(TT)=CQXH(TT)+CMXH(TT,KB)
      CQYH(TT)=CQYH(TT)+CMYH(TT,KB)
      CQZH(TT)=CQZH(TT)+CMZH(TT,KB)
   END DO !KB=1,NB
END IF !(IH == 1)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE CTCQ
!-----------------------------------------------------------------------------------------------!
