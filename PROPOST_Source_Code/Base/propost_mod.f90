!-----------------------------------------------------------------------------------------------!
!    PROPOST Module                                                                             !
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
MODULE PROPOST_MOD
!-----------------------------------------------------------------------------------------------!
!    Created by: 21052014, J. Baltazar, 2014 version 1.0                                        !
!    Modified  : 23052014, J. Baltazar, 2014 version 1.0                                        !
!    Modified  : 02072014, J. Baltazar, 2014 version 1.2                                        !
!    Modified  : 26012015, J. Baltazar, 2015 version 1.0                                        !
!    Modified  : 09032015, J. Baltazar, 2015 version 1.1 Unsteady Flow                          !
!    Modified  : 06072016, J. Baltazar, 2016 version 1.0                                        !
!    Modified  : 01022019, J. Baltazar, 2019 version 1.0 Harmonic Analysis                      !
!    Modified  : 20042020, J. Baltazar, 2020 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
DOUBLE PRECISION :: PI
PARAMETER (PI=3.141592653589793239D0)
!-----------------------------------------------------------------------------------------------!
!    Input and Output                                                                           !
!-----------------------------------------------------------------------------------------------!
CHARACTER*150 COMMENT(4)
CHARACTER*2   IDENTU(20)
!-----------------------------------------------------------------------------------------------!
!    Calculation Conditions                                                                     !
!-----------------------------------------------------------------------------------------------!
INTEGER :: IP,IN,IH,NB,NU,NWA,NREV,NTETA,INTER
INTEGER :: IROTOR,NT,ISTRIP,JI,JF,ISF
INTEGER :: NRS,NXS,NTS,NNHAR
DOUBLE PRECISION :: RE(20),P
DOUBLE PRECISION :: UU(20),RS(50),XS(50),TS(50)
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: RR
!-----------------------------------------------------------------------------------------------!
!    Dimensions                                                                                 !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
INTEGER :: NCP1,NCP,NC1,NC,NRP1,NRP
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
INTEGER :: NPW1,NPW,NRW1,NRW,NCPW
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
INTEGER :: NNXT1,NNXT,NNX1,NNX,NNTT,NNTP
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
INTEGER :: NNW1
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
INTEGER :: NHX1,NHX,NHTP,NHTT
!-----------------------------------------------------------------------------------------------!
!    Geometry and Metrics                                                                       !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XP,YP,ZP,RP,TP,SCP
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XP0,YP0,ZP0,UNXP0,UNYP0,UNZP0,AP0
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: ET1XP,ET1YP,ET1ZP,ET2XP,ET2YP,ET2ZP
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: AT1XP,AT1YP,AT1ZP,AT2XP,AT2YP,AT2ZP
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XPW,YPW,ZPW,RPW,TPW
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XPW0,YPW0,ZPW0
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XN,YN,ZN,RN,TN,SCN
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XN0,YN0,ZN0,UNXN0,UNYN0,UNZN0,AN0
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: ET1XN,ET1YN,ET1ZN,ET2XN,ET2YN,ET2ZN
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: AT1XN,AT1YN,AT1ZN,AT2XN,AT2YN,AT2ZN
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XNW,YNW,ZNW,RNW,TNW
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XH,YH,ZH,RH,TH
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XH0,YH0,ZH0,UNXH0,UNYH0,UNZH0,AH0
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: ET1XH,ET1YH,ET1ZH,ET2XH,ET2YH,ET2ZH
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: AT1XH,AT1YH,AT1ZH,AT2XH,AT2YH,AT2ZH
!-----------------------------------------------------------------------------------------------!
!    Potential                                                                                  !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: POTP,SOURCEP
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: POTPW,SOURCEPWCAV,POTPWP,POTPWS
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: POTN,SOURCEN
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: POTH,SOURCEH
!-----------------------------------------------------------------------------------------------!
!    Velocities and Pressure                                                                    !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: VT1P,VT2P,VTT1P,VTT2P
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: CPP,CPNP
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: VT1PW
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: VT1N,VT2N,VTT1N,VTT2N
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: CPN,CPNN
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: VT1H,VT2H,VTT1H,VTT2H
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: CPH,CPNH
!-----------------------------------------------------------------------------------------------!
!    Temporary Working Variables                                                                !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION :: WORK1,WORK2,WORK3
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: WORKA1,WORKA2,WORKA3,WORKA4,WORKA5
!-----------------------------------------------------------------------------------------------!
END MODULE PROPOST_MOD
!-----------------------------------------------------------------------------------------------!
