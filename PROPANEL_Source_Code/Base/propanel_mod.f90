!-----------------------------------------------------------------------------------------------!
!    PROPANEL Module                                                                            !
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
MODULE PROPANEL_MOD
!-----------------------------------------------------------------------------------------------!
!    Created by: 11102013, J. Baltazar, 2013 version 1.0                                        !
!    Modified  : 28102013, J. Baltazar, 2013 version 1.0                                        !
!    Modified  : 09052014, J. Baltazar, new variables ISC,ANGPITCH,ITHETA,ALPHAHT               !
!    Modified  : 01122014, J. Baltazar, 2014 version 1.2                                        !
!    Modified  : 09122014, J. Baltazar, 2014 version 1.3                                        !
!    Modified  : 06012015, J. Baltazar, 2015 version 1.0                                        !
!    Modified  : 19012015, J. Baltazar, 2015 version 1.1                                        !
!    Modified  : 24062015, J. Baltazar, 2015 version 1.2                                        !
!    Modified  : 02072015, J. Baltazar, 2015 version 1.3                                        !
!    Modified  : 05072017, J. Baltazar, 2017 version 1.0                                        !
!    Modified  : 04012018, J. Baltazar, 2018 version 1.0                                        !
!    Modified  : 09042020, J. Baltazar, 2020 version 1.0                                        !
!    Modified  : 17062020, J. Baltazar, 2020 version 1.1                                        !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
CHARACTER*8  IGEOM
DOUBLE PRECISION :: PI,TOL
PARAMETER (PI=3.141592653589793239D0,TOL=1.D-14)
!-----------------------------------------------------------------------------------------------!
!    Blade Variables                                                                            !
!-----------------------------------------------------------------------------------------------!
CHARACTER*10 IDENTP
INTEGER :: IP,NB,INTERP,ISTRIP,ISC,NRI,NCI,NRP1,NRP,NCP1,NCP,NC1,NC
DOUBLE PRECISION :: RPH,RMAX,ALPHAH,ALPHAT,ALPHALE,ALPHATE,ANGPITCH,PGAP,IHCORR
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: XP,YP,ZP,RP,TP
!-----------------------------------------------------------------------------------------------!
!    Blade Wake Variables                                                                       !
!-----------------------------------------------------------------------------------------------!
CHARACTER*10 IDENTPW
INTEGER :: INTERPW,NPW1,NPW,NRW1,NRW,NRIW,JI,JF,ITYPEPW,IMODELPW,ISTEADY,NTETA,INTECORR,NMW,NRM
DOUBLE PRECISION :: ST1PW,ST2PW,ST3PW,ST4PW,R0,R1,ADVJ,XPWW,XPWT,XI,XF,A0,A1
DOUBLE PRECISION :: RIW(50),PTW0I(50),PTWI(50),DPTW0,XMW(10),RMW(500),PMW(500)
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: XPW,YPW,ZPW,RPW,TPW
!-----------------------------------------------------------------------------------------------!
!    Nozzle Variables                                                                           !
!-----------------------------------------------------------------------------------------------!
CHARACTER*10 IDENTN
INTEGER :: IN,INTERN,IGRIDI,IGRIDO,NRNI,NRNO,NNT,NNT1,NNU,NND,NNX,NNX1,NNXT,NNTT,NNTT1,NCN
DOUBLE PRECISION :: ALPHAN,LD,CR,PTN
DOUBLE PRECISION :: XIL(50),YIL(50),XOL(50),YOL(50)
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: XN,YN,ZN,RN,TN
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake Variables                                                                      !
!-----------------------------------------------------------------------------------------------!
CHARACTER*10 IDENTNW
INTEGER :: ITYPENW,NNW1,NNW,ICONTRNW
DOUBLE PRECISION :: ST1NW,ST2NW,ST3NW,ST4NW,R2,XNWW,XNWT
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: XNW,YNW,ZNW,RNW,TNW
!-----------------------------------------------------------------------------------------------!
!    Hub Variables                                                                              !
!-----------------------------------------------------------------------------------------------!
CHARACTER*10 IDENTH
INTEGER :: IH,NHX1,NHX,NHTT,NHT1,NHT2,INTERH,NHI,NHT,NHU1,NHU,NHD,NCH,NHTT1,NHP
INTEGER :: IHR,ITHETA,ISTEP,ITERH
DOUBLE PRECISION :: XH0,XH3,ALPHAHT
DOUBLE PRECISION :: XHI(150),RHI(150),XHP(50),PTH(50)
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: XH,YH,ZH,RH,TH
!-----------------------------------------------------------------------------------------------!
END MODULE PROPANEL_MOD
!-----------------------------------------------------------------------------------------------!
