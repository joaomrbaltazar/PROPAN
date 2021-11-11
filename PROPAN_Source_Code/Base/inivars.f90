!-----------------------------------------------------------------------------------------------!
!    Initialise Variables                                                                       !
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
SUBROUTINE INIVARS
!-----------------------------------------------------------------------------------------------!
!    Created by: 29102013, J. Baltazar, version 1.0                                             !
!    Modified  : 02122013, J. Baltazar, version 1.0                                             !
!    Modified  : 26052014, J. Baltazar, Wake Alignment Module                                   !
!    Modified  : 06112014, J. Baltazar, version 3.0, Steady Cavitation Model                    !
!    Modified  : 11112014, J. Baltazar, version 3.1, Unsteady Cavitation Model                  !
!    Modified  : 27112014, J. Baltazar, version 3.3, Super-Cavitation                           !
!    Modified  : 29112014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!    Modified  : 10022016, J. Baltazar, 2016 version 1.0                                        !
!    Modified  : 27052016, J. Baltazar, 2016 version 1.2                                        !
!    Modified  : 24102016, J. Baltazar, 2016 version 1.4                                        !
!    Modified  : 23062020, J. Baltazar, 2020 version 1.2, Dynamic inflow                        !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------!
!    Geometry and Metrics                                                                       !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   ALLOCATE(XP0  (NCP,NRP),YP0  (NCP,NRP),ZP0  (NCP,NRP))
   ALLOCATE(UNXP0(NCP,NRP),UNYP0(NCP,NRP),UNZP0(NCP,NRP))
   ALLOCATE(ET1XP(NCP,NRP),ET1YP(NCP,NRP),ET1ZP(NCP,NRP))
   ALLOCATE(ET2XP(NCP,NRP),ET2YP(NCP,NRP),ET2ZP(NCP,NRP))
   ALLOCATE(AT1XP(NCP,NRP),AT1YP(NCP,NRP),AT1ZP(NCP,NRP))
   ALLOCATE(AT2XP(NCP,NRP),AT2YP(NCP,NRP),AT2ZP(NCP,NRP))
   ALLOCATE(UNXC0(NCP)    ,UNYC0(NCP)    ,UNZC0(NCP))
   ALLOCATE(AP0  (NCP,NRP))
!-----------------------------------------------------------------------------------------------!
   XP0  =0.D0
   YP0  =0.D0
   ZP0  =0.D0
   UNXP0=0.D0
   UNYP0=0.D0
   UNZP0=0.D0
   ET1XP=0.D0
   ET1YP=0.D0
   ET1ZP=0.D0
   ET2XP=0.D0
   ET2YP=0.D0
   ET2ZP=0.D0
   AT1XP=0.D0
   AT1YP=0.D0
   AT1ZP=0.D0
   AT2XP=0.D0
   AT2YP=0.D0
   AT2ZP=0.D0
   UNXC0=0.D0
   UNYC0=0.D0
   UNZC0=0.D0
   AP0  =0.D0
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   ALLOCATE(XPW0  (NPW,NRW),YPW0  (NPW,NRW),ZPW0  (NPW,NRW))
   ALLOCATE(UNXPW0(NPW,NRW),UNYPW0(NPW,NRW),UNZPW0(NPW,NRW))
   ALLOCATE(AT1XPW(NPW,NRW),AT1YPW(NPW,NRW),AT1ZPW(NPW,NRW))
   ALLOCATE(AT2XPW(NPW,NRW),AT2YPW(NPW,NRW),AT2ZPW(NPW,NRW))
   ALLOCATE(APW0  (NPW,NRW))
!-----------------------------------------------------------------------------------------------!
   XPW0  =0.D0
   YPW0  =0.D0
   ZPW0  =0.D0
   UNXPW0=0.D0
   UNYPW0=0.D0
   UNZPW0=0.D0
   AT1XPW=0.D0
   AT1YPW=0.D0
   AT1ZPW=0.D0
   AT2XPW=0.D0
   AT2YPW=0.D0
   AT2ZPW=0.D0
   APW0  =0.D0
   IF (NCPW /= 0) THEN
      ALLOCATE(ET1XPW(IABS(NCPW),NRW),ET1YPW(IABS(NCPW),NRW),ET1ZPW(IABS(NCPW),NRW))
      ET1XPW=0.D0
      ET1YPW=0.D0
      ET1ZPW=0.D0
   END IF !(NCPW /= 0)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   ALLOCATE(XN0  (NNXT1,NNTP),YN0  (NNXT1,NNTP),ZN0  (NNXT1,NNTP))
   ALLOCATE(UNXN0(NNXT1,NNTP),UNYN0(NNXT1,NNTP),UNZN0(NNXT1,NNTP))
   ALLOCATE(ET1XN(NNXT1,NNTP),ET1YN(NNXT1,NNTP),ET1ZN(NNXT1,NNTP))
   ALLOCATE(ET2XN(NNXT1,NNTP),ET2YN(NNXT1,NNTP),ET2ZN(NNXT1,NNTP))
   ALLOCATE(AT1XN(NNXT1,NNTP),AT1YN(NNXT1,NNTP),AT1ZN(NNXT1,NNTP))
   ALLOCATE(AT2XN(NNXT1,NNTP),AT2YN(NNXT1,NNTP),AT2ZN(NNXT1,NNTP))
   ALLOCATE(AN0  (NNXT1,NNTP))
!-----------------------------------------------------------------------------------------------!
   XN0  =0.D0
   YN0  =0.D0
   ZN0  =0.D0
   UNXN0=0.D0
   UNYN0=0.D0
   UNZN0=0.D0
   ET1XN=0.D0
   ET1YN=0.D0
   ET1ZN=0.D0
   ET2XN=0.D0
   ET2YN=0.D0
   ET2ZN=0.D0
   AT1XN=0.D0
   AT1YN=0.D0
   AT1ZN=0.D0
   AT2XN=0.D0
   AT2YN=0.D0
   AT2ZN=0.D0
   AN0  =0.D0
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   ALLOCATE(XNW0  (NNW,NNTP),YNW0  (NNW,NNTP),ZNW0  (NNW,NNTP))
   ALLOCATE(UNXNW0(NNW,NNTP),UNYNW0(NNW,NNTP),UNZNW0(NNW,NNTP))
   ALLOCATE(AT1XNW(NNW,NNTP),AT1YNW(NNW,NNTP),AT1ZNW(NNW,NNTP))
   ALLOCATE(AT2XNW(NNW,NNTP),AT2YNW(NNW,NNTP),AT2ZNW(NNW,NNTP))
   ALLOCATE(ANW0  (NNW,NNTP))
!-----------------------------------------------------------------------------------------------!
   XNW0  =0.D0
   YNW0  =0.D0
   ZNW0  =0.D0
   UNXNW0=0.D0
   UNYNW0=0.D0
   UNZNW0=0.D0
   AT1XNW=0.D0
   AT1YNW=0.D0
   AT1ZNW=0.D0
   AT2XNW=0.D0
   AT2YNW=0.D0
   AT2ZNW=0.D0
   ANW0  =0.D0
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
IF (IABS(IH) == 1) THEN
   ALLOCATE(XH0  (NHX,0:NHTP1),YH0  (NHX,0:NHTP1),ZH0  (NHX,0:NHTP1))
   ALLOCATE(UNXH0(NHX,  NHTP ),UNYH0(NHX,  NHTP ),UNZH0(NHX,  NHTP ))
   ALLOCATE(ET1XH(NHX,  NHTP ),ET1YH(NHX,  NHTP ),ET1ZH(NHX,  NHTP ))
   ALLOCATE(ET2XH(NHX,  NHTP ),ET2YH(NHX,  NHTP ),ET2ZH(NHX,  NHTP ))
   ALLOCATE(AT1XH(NHX,0:NHTP1),AT1YH(NHX,0:NHTP1),AT1ZH(NHX,0:NHTP1))
   ALLOCATE(AT2XH(NHX,0:NHTP1),AT2YH(NHX,0:NHTP1),AT2ZH(NHX,0:NHTP1))
   ALLOCATE(AH0  (NHX,  NHTP))
!-----------------------------------------------------------------------------------------------!
   XH0  =0.D0
   YH0  =0.D0
   ZH0  =0.D0
   UNXH0=0.D0
   UNYH0=0.D0
   UNZH0=0.D0
   ET1XH=0.D0
   ET1YH=0.D0
   ET1ZH=0.D0
   ET2XH=0.D0
   ET2YH=0.D0
   ET2ZH=0.D0
   AT1XH=0.D0
   AT1YH=0.D0
   AT1ZH=0.D0
   AT2XH=0.D0
   AT2YH=0.D0
   AT2ZH=0.D0
   AH0  =0.D0
END IF !(IABS(IH) == 1)
!-----------------------------------------------------------------------------------------------!
!    Numerical Differentiation                                                                  !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   ALLOCATE(A1P (NCP,NRP),B1P (NCP,NRP),C1P (NCP,NRP))
   ALLOCATE(A2P (NCP,NRP),B2P (NCP,NRP),C2P (NCP,NRP),D2P(NCP,NRP),E2P(NCP,NRP))
   ALLOCATE(AA2P(NCP,NRP),BB2P(NCP,NRP),CC2P(NCP,NRP))
   A1P =0.D0
   B1P =0.D0
   C1P =0.D0
   A2P =0.D0
   B2P =0.D0
   C2P =0.D0
   D2P =0.D0
   E2P =0.D0
   AA2P=0.D0
   BB2P=0.D0
   CC2P=0.D0
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
IF ((IP == 1).AND.(NCPW /= 0)) THEN
   ALLOCATE(A1PW(IABS(NCPW),NRW),B1PW(IABS(NCPW),NRW),C1PW(IABS(NCPW),NRW))
   A1PW=0.D0
   B1PW=0.D0
   C1PW=0.D0
END IF !((IP == 1).AND.(NCPW /= 0))
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   ALLOCATE(A1N(NNXT1,NNTP),B1N(NNXT1,NNTP),C1N(NNXT1,NNTP))
   ALLOCATE(A2N(NNXT1,NNTP),B2N(NNXT1,NNTP),C2N(NNXT1,NNTP),D2N(NNXT1,NNTP),E2N(NNXT1,NNTP))
   A1N=0.D0
   B1N=0.D0
   C1N=0.D0
   A2N=0.D0
   B2N=0.D0
   C2N=0.D0
   D2N=0.D0
   E2N=0.D0
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
IF (IABS(IH) == 1) THEN
   ALLOCATE(A1H(NHX,NHTP),B1H(NHX,NHTP),C1H(NHX,NHTP))
   ALLOCATE(A2H(NHX,NHTP),B2H(NHX,NHTP),C2H(NHX,NHTP),D2H(NHX,NHTP),E2H(NHX,NHTP))
   A1H=0.D0
   B1H=0.D0
   C1H=0.D0
   A2H=0.D0
   B2H=0.D0
   C2H=0.D0
   D2H=0.D0
   E2H=0.D0
END IF !(IABS(IH) == 1)
!-----------------------------------------------------------------------------------------------!
!    Velocity Field                                                                             !
!-----------------------------------------------------------------------------------------------!
IF (IFIELD == 1) ALLOCATE(VFX(NFX,NFY),VFY(NFX,NFY),VFZ(NFX,NFY))
!-----------------------------------------------------------------------------------------------!
!    Pressure Field                                                                             !
!-----------------------------------------------------------------------------------------------!
IF (IFIELD == 1) ALLOCATE(POTF(NFX,NFY,0:NT),DPOTFDT(NFX,NFY))
!-----------------------------------------------------------------------------------------------!
!    Influence Coefficients Matrix                                                              !
!-----------------------------------------------------------------------------------------------!
ALLOCATE(DIJ(NPAN1,NPAN1,NB),SIJ(NPAN,NPAN,NB))
DIJ=0.D0
SIJ=0.D0
IF (NT == 0) THEN
   ALLOCATE(WIJ(NPAN,NRT,NB))
   WIJ=0.D0
ELSE !(NT == 0)
   ALLOCATE(WIJ(NPAN,NWPAN,NB),KIJ(NPAN,NRT,2))
   WIJ=0.D0
   KIJ=0.D0
END IF !(NT == 0)
!-----------------------------------------------------------------------------------------------!
!    Potential                                                                                  !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   ALLOCATE(POTP(NCP,NRP,0:NT),SOURCEP(NCP,NRP,NB))
   POTP   =0.D0
   SOURCEP=0.D0
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   ALLOCATE(POTPW(NPW,NRW,0:NT),DPOTP(NRW))
   POTPW=0.D0
   DPOTP=0.D0
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   ALLOCATE(POTN(NNXT1,NNTP,0:NT),SOURCEN(NNXT1,NNTP,NB))
   POTN   =0.D0
   SOURCEN=0.D0
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   ALLOCATE(POTNW(NNW,NNTP,0:NT),DPOTN(NNTP))
   POTNW=0.D0
   DPOTN=0.D0
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
IF (IABS(IH) == 1) THEN
   ALLOCATE(POTH(NHX,NHTP,0:NT),SOURCEH(NHX,NHTP,NB))
   POTH   =0.D0
   SOURCEH=0.D0
END IF !(IABS(IH) == 1)
!-----------------------------------------------------------------------------------------------!
!    Velocities and Pressure                                                                    !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   ALLOCATE(VT1P(NCP,NRP),VT2P(NCP,NRP),VT2S(NCP,NRP),VTT1P(NCP,NRP),VTT2P(NCP,NRP))
   ALLOCATE(CPP(NCP,NRP,0:NT),CPNP(NCP,NRP,0:NT))
   VT1P =0.D0
   VT2P =0.D0
   VT2S =0.D0
   VTT1P=0.D0
   VTT2P=0.D0
   CPP  =0.D0
   CPNP =0.D0
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Gap                                                                                  !
!-----------------------------------------------------------------------------------------------!
IF (ISTRIP == 1) ALLOCATE(DCPG(NCP))
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   ALLOCATE(DCP(NRW))
   DCP=0.D0
   IF (NCPW /= 0) THEN
      ALLOCATE(VT1PW(IABS(NCPW),NRW))
      VT1PW=0.D0
   END IF !(NCPW /= 0)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   ALLOCATE(VT1N(NNXT1,NNTP),VT2N(NNXT1,NNTP),VTT1N(NNXT1,NNTP),VTT2N(NNXT1,NNTP))
   ALLOCATE(CPN(NNXT1,NNTP,0:NT),CPNN(NNXT1,NNTP,0:NT))
   VT1N =0.D0
   VT2N =0.D0
   VTT1N=0.D0
   VTT2N=0.D0
   CPN  =0.D0
   CPNP =0.D0
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   ALLOCATE(DCN(NNTP))
   DCN=0.D0
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
IF (IABS(IH) == 1) THEN
   ALLOCATE(VT1H(NHX,NHTP),VT2H(NHX,NHTP),VTT1H(NHX,NHTP),VTT2H(NHX,NHTP))
   ALLOCATE(CPH(NHX,NHTP,0:NT),CPNH(NHX,NHTP,0:NT))
   VT1H =0.D0
   VT2H =0.D0
   VTT1H=0.D0
   VTT2H=0.D0
   CPH  =0.D0
   CPNH =0.D0
END IF !(IABS(IH) == 1)
!-----------------------------------------------------------------------------------------------!
!    Cavitation                                                                                 !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   ALLOCATE(IDP(NRP,0:NT),IRP(NRP,0:NT),IDS(NRP,0:NT),IRS(NRP,0:NT))
   IDP=0
   IRP=0
   IDS=0
   IRS=0
   IF (NCAV > 0) THEN
      ALLOCATE(THICKP    (NCP,NRP,0:NT))
      ALLOCATE(SOURCEPCAV(NCP,NRP,0:NT))
      THICKP    =0.D0
      SOURCEPCAV=0.D0
   END IF !(NCAV > 0)
   IF ((NCAV > 0).AND.(IRED == 1)) THEN
      ALLOCATE(POTPWET(NCP,NRP,0:NT))
      POTPWET=0.D0
   END IF !((NCAV > 0).AND.(IRED == 1)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   ALLOCATE(IDPWP(NRW,0:NT),IRPWP(NRW,0:NT),IDPWS(NRW,0:NT),IRPWS(NRW,0:NT))
   IDPWP=0
   IRPWP=0
   IDPWS=0
   IRPWS=0
   IF ((NCAV > 0).AND.(IRED == 1)) THEN
      ALLOCATE(POTPWWET(NPW,NRW,0:NT))
      POTPWWET   =0.D0
   END IF !((NCAV > 0).AND.(IRED == 1)
   IF ((NCAV > 0).AND.(NCPW < 0)) THEN
      ALLOCATE(THICKPWP(IABS(NCPW),NRW,0:NT))
      THICKPWP   =0.D0
      ALLOCATE(CAMBERPWP(IABS(NCPW),NRW,0:NT))
      CAMBERPWP  =0.D0
      ALLOCATE(POTPWP(IABS(NCPW),NRW,0:NT),SOURCEPWCAV(IABS(NCPW),NRW,0:NT))
      POTPWP     =0.D0
      SOURCEPWCAV=0.D0
      IF (IRED == 1) THEN
         ALLOCATE(POTPWPWET(IABS(NCPW),NRW,0:NT))
         POTPWPWET=0.D0
      END IF !(IRED == 1)
   END IF !((NCAV > 0).AND.(NCPW < 0))
   IF ((NCAV > 0).AND.(NCPW > 0)) THEN
      ALLOCATE(THICKPWS(IABS(NCPW),NRW,0:NT))
      THICKPWS   =0.D0
      ALLOCATE(CAMBERPWS(IABS(NCPW),NRW,0:NT))
      CAMBERPWS  =0.D0
      ALLOCATE(POTPWS(IABS(NCPW),NRW,0:NT),SOURCEPWCAV(IABS(NCPW),NRW,0:NT))
      POTPWS     =0.D0
      SOURCEPWCAV=0.D0
      IF (IRED == 1) THEN
         ALLOCATE(POTPWSWET(IABS(NCPW),NRW,0:NT))
         POTPWSWET=0.D0
      END IF !(IRED == 1)
   END IF !((NCAV > 0).AND.(NCPW > 0))
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   IF ((NCAV > 0).AND.(IRED == 1)) THEN
      ALLOCATE(POTNWET(NNXT1,NNTP,0:NT))
      POTNWET=0.D0
   END IF !((NCAV > 0).AND.(IRED == 1)
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   IF ((NCAV > 0).AND.(IRED == 1)) THEN
      ALLOCATE(POTNWWET(NNW,NNTP,0:NT))
      POTNWWET=0.D0
   END IF !((NCAV > 0).AND.(IRED == 1)
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
IF (IABS(IH) == 1) THEN
   IF ((NCAV > 0).AND.(IRED == 1)) THEN
      ALLOCATE(POTHWET(NHX,NHTP,0:NT))
      POTHWET=0.D0
   END IF !((NCAV > 0).AND.(IRED == 1)
END IF !(IABS(IH) == 1)
!-----------------------------------------------------------------------------------------------!
!    Force Coefficients                                                                         !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   ALLOCATE(CFXP(0:NT,NB),CMXP(0:NT,NB),CFYP(0:NT,NB),CMYP(0:NT,NB),CFZP(0:NT,NB),CMZP(0:NT,NB))
   ALLOCATE(CTXP(0:NT),CQXP(0:NT),CTYP(0:NT),CQYP(0:NT),CTZP(0:NT),CQZP(0:NT))
   CFXP=0.D0
   CMXP=0.D0
   CFYP=0.D0
   CMYP=0.D0
   CFZP=0.D0
   CMZP=0.D0
   CTXP=0.D0
   CQXP=0.D0
   CTYP=0.D0
   CQYP=0.D0
   CTZP=0.D0
   CQZP=0.D0
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   ALLOCATE(CFXN(0:NT,NB),CMXN(0:NT,NB),CFYN(0:NT,NB),CMYN(0:NT,NB),CFZN(0:NT,NB),CMZN(0:NT,NB))
   ALLOCATE(CTXN(0:NT),CQXN(0:NT),CTYN(0:NT),CQYN(0:NT),CTZN(0:NT),CQZN(0:NT))
   CFXN=0.D0
   CMXN=0.D0
   CFYN=0.D0
   CMYN=0.D0
   CFZN=0.D0
   CMZN=0.D0
   CTXN=0.D0
   CQXN=0.D0
   CTYN=0.D0
   CQYN=0.D0
   CTZN=0.D0
   CQZN=0.D0
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
IF (IABS(IH) == 1) THEN
   ALLOCATE(CFXH(0:NT,NB),CMXH(0:NT,NB),CFYH(0:NT,NB),CMYH(0:NT,NB),CFZH(0:NT,NB),CMZH(0:NT,NB))
   ALLOCATE(CTXH(0:NT),CQXH(0:NT),CTYH(0:NT),CQYH(0:NT),CTZH(0:NT),CQZH(0:NT))
   CFXH=0.D0
   CMXH=0.D0
   CFYH=0.D0
   CMYH=0.D0
   CFZH=0.D0
   CMZH=0.D0
   CTXH=0.D0
   CQXH=0.D0
   CTYH=0.D0
   CQYH=0.D0
   CTZH=0.D0
   CQZH=0.D0
END IF !(IABS(IH) == 1)
!-----------------------------------------------------------------------------------------------!
!    Right Hand Side                                                                            !
!-----------------------------------------------------------------------------------------------!
ALLOCATE(SI(NPAN1),RHS(NPAN))
SI =0.D0
RHS=0.D0
!-----------------------------------------------------------------------------------------------!
!    LU Decomposition Vectors                                                                   !
!-----------------------------------------------------------------------------------------------!
ALLOCATE(IPVT1(NPAN),IPVT2(NPAN1),IPVT3(NPAN),IPVT4(NPAN1))
IPVT1=0
IPVT2=0
IPVT3=0
IPVT4=0
IF (NCAV > 0) THEN
   ALLOCATE(IPVTC1(NPAN),IPVTC2(NPAN1),IPVTC4(NPAN1))
   IPVTC1=0
   IPVTC2=0
   IPVTC4=0
END IF !(NCAV > 0)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE INIVARS
!-----------------------------------------------------------------------------------------------!
