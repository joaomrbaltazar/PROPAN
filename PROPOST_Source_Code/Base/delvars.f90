!-----------------------------------------------------------------------------------------------!
!    Delete variables                                                                           !
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
SUBROUTINE DELVARS
!-----------------------------------------------------------------------------------------------!
!    Created by: 21052014, J. Baltazar, 2014 version 1.0                                        !
!    Modified  : 19012015, J. Baltazar, 2015 version 1.0                                        !
!    Modified  : 09032015, J. Baltazar, 2015 version 1.1 Unsteady Flow                          !
!    Modified  : 07042015, J. Baltazar, 2015 version 1.2                                        !
!    Modified  : 06072016, J. Baltazar, 2016 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
USE PROPOST_MOD
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------!
!    Geometry and Metrics                                                                       !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   DEALLOCATE(XP,YP,ZP,RP,TP,SCP)
   DEALLOCATE(XP0,YP0,ZP0,AP0)
   DEALLOCATE(UNXP0,UNYP0,UNZP0)
   DEALLOCATE(ET1XP,ET1YP,ET1ZP,ET2XP,ET2YP,ET2ZP)
   DEALLOCATE(AT1XP,AT1YP,AT1ZP,AT2XP,AT2YP,AT2ZP)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   DEALLOCATE(XPW,YPW,ZPW,RPW,TPW)
   DEALLOCATE(RR)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   DEALLOCATE(XN,YN,ZN,RN,TN,SCN)
   DEALLOCATE(XN0,YN0,ZN0,AN0)
   DEALLOCATE(UNXN0,UNYN0,UNZN0)
   DEALLOCATE(ET1XN,ET1YN,ET1ZN,ET2XN,ET2YN,ET2ZN)
   DEALLOCATE(AT1XN,AT1YN,AT1ZN,AT2XN,AT2YN,AT2ZN)
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   DEALLOCATE(XNW,YNW,ZNW,RNW,TNW)
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
IF (IABS(IH) == 1) THEN
   DEALLOCATE(XH,YH,ZH,RH,TH)
   DEALLOCATE(XH0,YH0,ZH0,AH0)
   DEALLOCATE(UNXH0,UNYH0,UNZH0)
   DEALLOCATE(ET1XH,ET1YH,ET1ZH,ET2XH,ET2YH,ET2ZH)
   DEALLOCATE(AT1XH,AT1YH,AT1ZH,AT2XH,AT2YH,AT2ZH)
END IF !(IABS(IH) == 1)
!-----------------------------------------------------------------------------------------------!
!    Potential                                                                                  !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) DEALLOCATE(POTP,SOURCEP)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) DEALLOCATE(POTPW)
IF ((IP == 1).AND.(NCPW < 0)) DEALLOCATE(SOURCEPWCAV,POTPWP)
IF ((IP == 1).AND.(NCPW > 0)) DEALLOCATE(SOURCEPWCAV,POTPWS)
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) DEALLOCATE(POTN,SOURCEN)
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
IF (IABS(IH) == 1) DEALLOCATE(POTH,SOURCEH)
!-----------------------------------------------------------------------------------------------!
!    Velocities and Pressure                                                                    !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   DEALLOCATE(VT1P,VT2P,VTT1P,VTT2P)
   DEALLOCATE(CPP,CPNP)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
IF ((IP == 1).AND.(NCPW /= 0)) DEALLOCATE(VT1PW)
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   DEALLOCATE(VT1N,VT2N,VTT1N,VTT2N)
   DEALLOCATE(CPN,CPNN)
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
IF (IABS(IH) == 1) THEN
   DEALLOCATE(VT1H,VT2H,VTT1H,VTT2H)
   DEALLOCATE(CPH,CPNH)
END IF !(IABS(IH) == 1)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE DELVARS
!-----------------------------------------------------------------------------------------------!
