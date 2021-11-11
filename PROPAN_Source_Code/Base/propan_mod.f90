!-----------------------------------------------------------------------------------------------!
!    PROPAN Module                                                                              !
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
MODULE PROPAN_MOD
!-----------------------------------------------------------------------------------------------!
!    Created by: 28102013, J. Baltazar, version 1.0                                             !
!    Modified  : 03122013, J. Baltazar, version 1.0                                             !
!    Modified  : 29052014, J. Baltazar, Wake Alignment Model 1                                  !
!    Modified  : 29102014, J. Baltazar, Wake Alignment Model 2                                  !
!    Modified  : 03112014, J. Baltazar, Interpolation scheme at nozzle t.e.                     !
!    Modified  : 06112014, J. Baltazar, version 3.0, Steady Cavitation Model                    !
!    Modified  : 12112014, J. Baltazar, version 3.1, Unsteady Cavitation Model                  !
!    Modified  : 19112014, J. Baltazar, version 3.2, Mid-Chord Cavitation Model                 !
!    Modified  : 26112014, J. Baltazar, version 3.3, Super-Cavitation Model                     !
!    Modified  : 12112014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!    Modified  : 10022016, J. Baltazar, 2016 version 1.0                                        !
!    Modified  : 11042016, J. Baltazar, 2016 version 1.1                                        !
!    Modified  : 31052016, J. Baltazar, 2016 version 1.2                                        !
!    Modified  : 21062016, J. Baltazar, 2016 version 1.3                                        !
!    Modified  : 20102016, J. Baltazar, 2016 version 1.4, Broyden's Method for IPKC             !
!    Modified  : 04012018, J. Baltazar, 2018 version 1.0, Gap strip alignment                   !
!    Modified  : 02102018, J. Baltazar, 2018 version 1.1, Duct Kutta points at inner and outer  !
!    Modified  : 23062020, J. Baltazar, 2020 version 1.2, Dynamic inflow                        !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
DOUBLE PRECISION :: PI,TOL,FARTOL1,FARTOL2
PARAMETER (PI=3.141592653589793239D0,TOL=1.D-10)
PARAMETER (FARTOL1=3.0D0,FARTOL2=6.5D0)
!-----------------------------------------------------------------------------------------------!
!    Input and Output                                                                           !
!-----------------------------------------------------------------------------------------------!
CHARACTER*150 COMMENT(4)
CHARACTER*3   ICHAR
CHARACTER*2   IDENTU(20)
!-----------------------------------------------------------------------------------------------!
!    Calculation Conditions                                                                     !
!-----------------------------------------------------------------------------------------------!
INTEGER :: IP,IN,IH,NB,NU
INTEGER :: NKIT,MKIT,IK,NWA,NWALIGN,NGAP,NCAV,NREV,NTETA,NWR,NWT,MMAX
INTEGER :: JI,JF,IROTOR,IFREQ,IPAN,INTE,IFARP,ISTRIP,ISOLVER,IFIELD,ICP,ITE,IFARV,IFILESOL
INTEGER :: IWA,IHUB,ITIP,IBD,NN,NBASIS,CCI,CCF,CN,IRED,IMCP,IMCS
INTEGER :: NPAN1,NPAN,NWPAN,NRT,NT,NHI,NRV,JV(200),NWV
DOUBLE PRECISION :: UU(20),DTETA,P,TOLS,CQ,REXTRA,FREX,TOLG,TOLK,BETA,EPS,XV(200)
DOUBLE PRECISION :: DELTA,ALFV,SIGMA,FN,LCMAXS,RACAES,ETAS,LCMAXP,RACAEP,ETAP,CREX
DOUBLE PRECISION :: ERRLCS,ERRACS,ERRTCS,ERRLCP,ERRACP,ERRTCP,TOLC
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: WRR,WTT
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: VX,VT,VR
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: ANVX,BNVX,ANVT,BNVT,ANVR,BNVR
!-----------------------------------------------------------------------------------------------!
!    Geomtrical Definition                                                                      !
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
INTEGER :: INTERPW
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
CHARACTER*10 IDENTN
INTEGER :: INTERN,NRNI,NRNO,NNP,IGRIDI,IGRIDO
DOUBLE PRECISION :: LD,CR,PTN,PGAP,PREX,XITE,XOTE
DOUBLE PRECISION :: XIL(50),YIL(50),XOL(50),YOL(50)
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: XNN,PNN
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
INTEGER :: INTERH,NHP
DOUBLE PRECISION :: XHI(150),RHI(150),XHP(50),PTH(50)
!-----------------------------------------------------------------------------------------------!
!    Dimensions                                                                                 !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
INTEGER :: NCP1,NCP,NRP1,NRP,NC1,NC,NPPAN
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
INTEGER :: NPW1,NPW,NRW1,NRW,NPWPAN
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
INTEGER :: NNXT1,NNXT,NNX1,NNX,NNTT1,NNTT,NNTP,NNT2,NNT,NNT1,NNU,NNU1,NND,NN2,NCN,NNPAN
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
INTEGER :: NNW1,NNW,NNWPAN
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
INTEGER :: NHX1,NHX,NHTP1,NHTP,NHTT1,NHTT,NHT2,NHT1,NHT,NHU1,NHU,NHD,NH2,NCH,NHPAN
!-----------------------------------------------------------------------------------------------!
!    Velocity Field                                                                             !
!-----------------------------------------------------------------------------------------------!
INTEGER :: NFX,NFY
!-----------------------------------------------------------------------------------------------!
!    Geometry and Metrics                                                                       !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XP,YP,ZP,RP,TP
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XP0,YP0,ZP0,UNXP0,UNYP0,UNZP0,AP0
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: ET1XP,ET1YP,ET1ZP,ET2XP,ET2YP,ET2ZP
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: AT1XP,AT1YP,AT1ZP,AT2XP,AT2YP,AT2ZP
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: UNXC0,UNYC0,UNZC0
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XPW,YPW,ZPW,RPW,TPW
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XPW0,YPW0,ZPW0,UNXPW0,UNYPW0,UNZPW0,APW0
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: ET1XPW,ET1YPW,ET1ZPW
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: AT1XPW,AT1YPW,AT1ZPW,AT2XPW,AT2YPW,AT2ZPW
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XN,YN,ZN,RN,TN
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XN0,YN0,ZN0,UNXN0,UNYN0,UNZN0,AN0
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: ET1XN,ET1YN,ET1ZN,ET2XN,ET2YN,ET2ZN
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: AT1XN,AT1YN,AT1ZN,AT2XN,AT2YN,AT2ZN
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XNW,YNW,ZNW,RNW,TNW
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XNW0,YNW0,ZNW0,UNXNW0,UNYNW0,UNZNW0,ANW0
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: AT1XNW,AT1YNW,AT1ZNW,AT2XNW,AT2YNW,AT2ZNW
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XH,YH,ZH,RH,TH
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XH0,YH0,ZH0,UNXH0,UNYH0,UNZH0,AH0
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: ET1XH,ET1YH,ET1ZH,ET2XH,ET2YH,ET2ZH
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: AT1XH,AT1YH,AT1ZH,AT2XH,AT2YH,AT2ZH
!-----------------------------------------------------------------------------------------------!
!    Numerical Differentiation                                                                  !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: A1P,B1P,C1P,A2P,B2P,C2P,D2P,E2P,AA2P,BB2P,CC2P
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: A1PW,B1PW,C1PW
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: A1N,B1N,C1N,A2N,B2N,C2N,D2N,E2N
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: A1H,B1H,C1H,A2H,B2H,C2H,D2H,E2H
!-----------------------------------------------------------------------------------------------!
!    Velocity Field                                                                             !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: XVF,YVF,ZVF
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: VFX,VFY,VFZ
!-----------------------------------------------------------------------------------------------!
!    Pressure Field                                                                             !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: DPOTFDT
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: POTF
!-----------------------------------------------------------------------------------------------!
!    Influence Coefficients Matrix                                                              !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: DIJ,SIJ,WIJ,KIJ
!-----------------------------------------------------------------------------------------------!
!    Potential                                                                                  !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: POTP,SOURCEP
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: POTPW
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: DPOTP
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: POTN,SOURCEN
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: POTNW
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: DPOTN
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: POTH,SOURCEH
!-----------------------------------------------------------------------------------------------!
!    Velocities and Pressure                                                                    !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: VT1P,VT2P,VT2S,VTT1P,VTT2P
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: CPP,CPNP
!-----------------------------------------------------------------------------------------------!
!    Blade Gap                                                                                  !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: DCPG
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: VT1PW
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: DCP
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: VT1N,VT2N,VTT1N,VTT2N
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: CPN,CPNN
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: DCN
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: VT1H,VT2H,VTT1H,VTT2H
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: CPH,CPNH
!-----------------------------------------------------------------------------------------------!
!    Cavitation                                                                                 !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
INTEGER :: NPCAV,NSCAV
INTEGER,ALLOCATABLE,DIMENSION(:,:)            :: IDP,IRP,IDS,IRS
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: THICKP,SOURCEPCAV
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: POTPWET
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
INTEGER :: NPWCAV,NWCAV,NCPW
INTEGER,ALLOCATABLE,DIMENSION(:,:)            :: IDPWP,IRPWP,IDPWS,IRPWS
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: THICKPWP,THICKPWS
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: CAMBERPWP,CAMBERPWS
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: POTPWP,POTPWS,SOURCEPWCAV
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: POTPWWET,POTPWPWET,POTPWSWET
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: POTNWET
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: POTNWWET
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: POTHWET
!-----------------------------------------------------------------------------------------------!
!    Force Coefficients                                                                         !
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: CFXP,CMXP,CFYP,CMYP,CFZP,CMZP
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: CTXP,CQXP,CTYP,CQYP,CTZP,CQZP
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: CFXN,CMXN,CFYN,CMYN,CFZN,CMZN
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: CTXN,CQXN,CTYN,CQYN,CTZN,CQZN
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: CFXH,CMXH,CFYH,CMYH,CFZH,CMZH
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: CTXH,CQXH,CTYH,CQYH,CTZH,CQZH
!-----------------------------------------------------------------------------------------------!
!    Right Hand Side                                                                            !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: SI,RHS
!-----------------------------------------------------------------------------------------------!
!    LU Decomposition Vectors                                                                   !
!-----------------------------------------------------------------------------------------------!
INTEGER,ALLOCATABLE,DIMENSION(:)              :: IPVT1,IPVT2,IPVT3,IPVT4
INTEGER,ALLOCATABLE,DIMENSION(:)              :: IPVTC1,IPVTC2,IPVTC4
!-----------------------------------------------------------------------------------------------!
!    Temporary Working Variables                                                                !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION                              :: WORK1,WORK2,WORK3
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: WORKA1,WORKA2,WORKA3
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: WORKM1,WORKM2
!-----------------------------------------------------------------------------------------------!
END MODULE PROPAN_MOD
!-----------------------------------------------------------------------------------------------!
