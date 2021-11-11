!-----------------------------------------------------------------------------------------------!
!    IST SURFACE PANEL METHOD CODE FOR ROTOR ANALYSIS IN UNIFORM ONSET FLOW                     !
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
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.                     !                                                                          !
!                                                                                               !
!    Modified  : 03122013, J. Baltazar, 2014 version 1.0                                        !
!    Modified  : 29052014, J. Baltazar, 2014 version 2.0, Wake Alignment Model 1                !
!    Modified  : 29102014, J. Baltazar, 2014 version 2.1, Wake Alignment Model 2                !
!    Modified  : 03112014, J. Baltazar, 2014 version 2.2, Interpolation scheme at nozzle t.e.   !
!    Modified  : 06112014, J. Baltazar, 2014 version 3.0, Steady Cavitation Model               !
!    Modified  : 12112014, J. Baltazar, 2014 version 3.1, Unsteady Cavitation Model             !
!    Modified  : 19112014, J. Baltazar, 2014 version 3.2, Mid-Chord Cavitation                  !
!    Modified  : 28112014, J. Baltazar, 2014 version 3.3, Steady Super-Cavitation Model         !
!    Modified  : 12122014, J. Baltazar, 2014 version 3.4, Unsteady Super-Cavitation Model       !
!    Modified  : 07012015, J. Baltazar, 2015 version 1.0, Wing Case                             !
!    Modified  : 02022015, J. Baltazar, Gap sources                                             !
!    Modified  : 08032015, J. Baltazar, 2015 version 1.1, Correction for Super-Cavitation       !
!    Modified  : 06102015, J. Baltazar, 2015 version 1.2                                        !
!    Modified  : 29022016, J. Baltazar, 2016 version 1.0                                        !
!    Modified  : 04052016, J. Baltazar, 2016 version 1.1                                        !
!    Modified  : 31052016, J. Baltazar, 2016 version 1.2, Reduced System for Super-Cavitation   !
!    Modified  : 08072016, J. Baltazar, 2016 version 1.3                                        !
!    Modified  : 26102016, J. Baltazar, 2016 version 1.4, Broyden's Method for IPKC             !
!    Modified  : 06072017, J. Baltazar, 2017 version 1.0                                        !
!    Modified  : 12012018, J. Baltazar, 2018 version 1.0, Gap strip alignment                   !
!    Modified  : 03102018, J. Baltazar, 2018 version 1.1, Duct Kutta points at inner and outer  !
!    Modified  : 18042019, J. Baltazar, 2019 version 1.0, Wake alignment robustness             !
!    Modified  : 16042020, J. Baltazar, 2020 version 1.0                                        !
!    Modified  : 25052020, J. Baltazar, 2020 version 1.1, Output corrected                      !
!    Modified  : 23122020, J. Baltazar, 2020 version 1.2, Dynamic inflow                        !
!    Modified  : 05112021, J. Baltazar, 2021 version 1.0, GitHub                                !
!-----------------------------------------------------------------------------------------------!
PROGRAM PROPAN
!-----------------------------------------------------------------------------------------------!
!    Input Description                                                                          !
!                                                                                               !
!    The input is made from two files:                                                          !
!                                                                                               !
!    PROPAN.INP   : Rotor operational conditions and control variables                          !
!                                                                                               !
!    PROPANEL.DAT : Blade-nozzle-hub panel geometry                                             !
!                                                                                               !
!    WAKE.DAT     : Inflow wake                                                                 !
!                                                                                               !
!    VFIELD.DAT   : Field points in the inertial reference frame                                !
!-----------------------------------------------------------------------------------------------!
!    Input variables in PROPAN.INP                                                              !
!                                                                                               !
!    COMMENT(5)=    Comments for run description (5 lines)                                      !
!                                                                                               !
!    IP        = 1: Blade and blade wake panelling (0 no panelling)                             !
!    IN        = 1: Nozzle and nozzle wake panelling (0 no panelling)                           !
!    IH        = 1: Hub panelling (0 no panelling)                                              !
!              =-1: No-hub but closure of blade root with panels                                !
!    NB        =    Number of blades                                                            !
!    JI        =    Initial strip of the blade wake                                             !
!    JF        =    Last strip of the blade wake                                                !
!                                                                                               !
!    IROTOR    =    Choice on the rotor                                                         !
!                   IROTOR=-1: wing                                                             !
!                   IROTOR= 0: propeller                                                        !
!                   IROTOR= 1: turbine                                                          !
!    NU        =    Number of inflow conditions for the analysis                                !
!    IDENTU(NU)=    Inflow condition identification (2 characters)                              !
!    UU(NU)    =    Angle of attack in degrees   (IF IROTOR = -1)                               !
!                   Advance coefficient J=U/nD   (IF IROTOR =  0)                               !
!                   Tip-speed-ratio TSR=OmegaR/U (IF IROTOR =  1)                               !
!                                                                                               !
!    NKIT      =    Choice on the Kutta condition                                               !
!                   NKIT =0: linear Kutta condition                                             !
!                   NKIT>=1: non-linear pressure Kutta condition                                !
!    MKIT      =    Choice on the Kutta Model                                                   !
!                   MKIT=0: Newton-Raphson method (default method)                              !
!                   MKIT=1: Broyden's method                                                    !
!    IK        =    Order of the tolerance of the pressure jump                                 !
!    TOLK      =    Numerical tolerance of pressure jump at the Kutta or wake control points    !
!    BETA      =    Perturbation parameter of iterative pressure Kutta condition                !
!    XITE      =    Inner axial coordinate for pressure jump at nozzle t.e. (IF IN = 1)         !
!    XOTE      =    Outer axial coordinate for pressure jump at nozzle t.e. (IF IN = 1)         !
!                                                                                               !
!    NWA       =    Number of wake alignment iterations                                         !
!    IWA       =    Choice on the wake model                                                    !
!                   IWA=1: Model 1                                                              !
!                   IWA=2: Model 2                                                              ! 
!    NWALIGN   =    Number of panels to be aligned                                              !
!    ITE       =    Mangler-Smith condition at the TE                                           !
!    NRV       =    Number of radial points to compute the induced velocities                   !
!    JV(NRV)   =    Selection of the radial strips                                              !
!    NWV       =    Number of streamwise points to compute the induced velocities (IF IWA == 2) !
!    XV(NWV)   =    Selection of the streamwise strips (axial coordinate / Rp) (IF IWA == 2)    !
!    NBASIS    =    Number of basis functions (IF IWA == 2)                                     !
!    CCI,CCF   =    Initial and final exponents of the exponential function (IF IWA == 2)       !
!    EPS       =    Numerical tolerance for the computation of the induced velocities           !
!    IFARV     =    Control variable for the use of far-field formulas for velocities           !
!                   IFARV=0: no use of far-field formulas                                       !
!                   IFARV=1: use of far-field formulas                                          !
!                                                                                               !
!    INTERPW   =    Choice on the interpolations                                                !
!                   INTERPW=0: linear interpolation                                             !
!                   INTERPW=1: quadratic interpolation                                          !
!                   INTERPW=2: cubic interpolation                                              !
!                                                                                               !
!    IBD       = 1: Boundary layer correction model (0 no correction)                           !
!    DELTA     =    Boundary layer thickness / Rp                                               !
!    NN        =    Exponent of the power law velocity distribution                             !
!    ALFV      =    Extrapolation value to duct inner surface                                   !
!                                                                                               !
!    IHUB      = 1: Hub correction (0 no correction) (IF IWA == 2)                              !
!    ITIP      = 1: Tip correction (0 no correction) (IF IWA == 2)                              !
!    REXTRA    =    Radial coordinate correction (IF IWA == 2)                                  !
!                                                                                               !
!    IDENTN    =    Nozzle identification (10 characters)                                       !
!    LD        =    Length ratio of the nozzle                                                  !
!    CR        =    Clearance of the nozzle / Rp                                                !
!    IGRIDI    =    Grid topology in the inner side of the nozzle                               !
!                   IGRIDI=0: conventional topology                                             !
!                   IGRIDI=1: topology with tip section blade geometry                          !
!    IGRIDO    =    Grid topology in the outer side of the nozzle                               !
!                   IGRIDO=0: conventional topology                                             !
!                   IGRIDO=1: topology with tip section blade geometry                          !
!    INTERN    =    Choice on the interpolations (IF IN = 1)                                    !
!                   INTERN=0: linear interpolation                                              !
!                   INTERN=1: quadratic interpolation                                           !
!    NRNI      =    Number of input points of the inner surface of the nozzle (max 50)          !
!    XIL(NRNI) =    Axial coordinates of nozzle input points at the inner side / Ln             !
!    YIL(NRNI) =    Inner coordinates of nozzle input points / Ln                               !
!    NRNO      =    Number of input points of the outer surface of the nozzle (max 50)          !
!    XOL(NRNO) =    Axial coordinates of nozzle input points at the outer side / Ln             !
!    YOL(NRNO) =    Outer coordinates of nozzle input points / Ln                               !
!                                                                                               !
!    INTERH    =    Choice on the interpolations (IF IH = 1)                                    !
!                   INTERH=0: linear interpolation                                              !
!                   INTERH=1: quadratic interpolation                                           !
!                   INTERH=2: cubic interpolation                                               !
!    NHI       =    Number of input points of hub (max 100) (IF IH = 1)                         !
!    XHI(NHI)  =    Axial coordinates of hub input points / Rp (IF IH = 1)                      !
!    RHI(NHI)  =    Radial coordinates of hub input points / Rp (IF IH = 1)                     !
!    NHP       =    Number of input point for pitch distribution (max 50) (IF IH = 1)           !
!    XHP(NHP)  =    Axial coordinate of pitch distribution / Rp (IF IH = 1)                     !
!    PTH(NHP)  =    Pitch distribution on the hub / Diameter (IF IH = 1)                        !
!                                                                                               !
!    PGAP      =    Pitch at blade tip                                                          !
!    PREX      =    Relaxation factor of the gap strip                                          !
!                                                                                               !
!    NGAP      =    Number of gap model iterations                                              !
!    CQ        =    Discharge coefficient                                                       !
!    FREX      =    Relaxation factor of the gap model                                          !
!    TOLG      =    Numerical tolerance of transpiration velocity at the gap                    !
!                                                                                               !
!    NCAV      =    Number of cavity model iterations                                           !
!    SIGMA     =    Cavitation number (Vref=nD)                                                 !
!    FN        =    Froude number                                                               !
!    IRED      =    System of equations for cavitating flow                                     !
!                   IRED=0: complete system                                                     !
!                   IRED=1: reduced system                                                      !
!    IMCP      =    Mid-chord cavitation model on pressure side                                 !
!                   IMCP=0: no mid-chord cavitation                                             !
!                   IMCP>0: mid-chord cavitation model (thickness correction for IMCP panels)   !
!    IMCS      =    Mid-chord cavitation model on suction side                                  !
!                   IMCS=0: no mid-chord cavitation                                             !
!                   IMCS>0: mid-chord cavitation model (thickness correction for IMCS panels)   !
!    NCPW      =    Number of streamwise blade wake panels for super-cavitation                 !
!                   NCPW > 0: super-cavitation on upper side                                    !
!                   NCPW < 0: super-cavitation on lower side                                    !
!    CREX      =    Relaxation factor for unsteady terms in DBC and KBC                         !
!    TOLC      =    Numerical tolerance of cavitation model                                     !
!                                                                                               !
!    NREV      =    Number of revolutions                                                       !
!    NTETA     =    Number of theta steps                                                       !
!    IFREQ     =    Number of frequencies of the Fourier analysis                               !
!    IFILESOL  =    Starting solution file for unsteady calculations                            !
!                   IFILESOL=0: steady solution from circumferentially averaged wake            !
!                   IFILESOL=1: solution file                                                   !
!                                                                                               !
!    IPAN      =    Control variable for type of panel                                          !
!                   IPAN=0: flat panels                                                         !
!                   IPAN=1: hyperboloidal panels                                                !
!    INTE      =    Control variable for integration on the panel                               !
!                   INTE=0: analytical integration                                              !
!                   INTE=1: numerical integration                                               !
!                   INTE=2: numerical integration for wake panels                               !
!    TOLS      =    Numerical tolerance of numerical integration on the panels                  !
!    MMAX      =    Maximum number of subdivisions in the numerical integration                 !
!    IFARP     =    Control variable for the use of far-field formulas for potential            !
!                   IFARP=0: no use of far-field formulas                                       !
!                   IFARP=1: use of far-field formulas                                          !
!    ISTRIP    =    Input option to exclude last tip panel strip                                !
!                   ISTRIP=0: zero and finite tip (default value)                               !
!                   ISTRIP=1: include last panel strip for gap model                            !
!                   ISTRIP=2: include last panel strip for closed tip                           !
!    ISOLVER   =    Choice of solver for linear system of equations                             !
!                   ISOLVER=0: LU factorization                                                 !
!                   ISOLVER=1: iterative bi-conjugate gradient                                  !
!                                                                                               !
!    IFIELD    = 1: Calculation of the velocity and pressure fields (0 no calculation)          !
!    ICP       =    Cp definition                                                               !
!                   ICP=1: Vref=Vinf                                                            !
!                   ICP=2: Vref=OmegaR                                                          !
!-----------------------------------------------------------------------------------------------!
!    Declarations                                                                               !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
CHARACTER*50 OUTFILE
LOGICAL :: FCHECK
INTEGER :: JJ,TT,WW,II,CC,IWORK
REAL    :: TIME,TIME1,TIME2
DOUBLE PRECISION :: ERRG,ERRC
!-----------------------------------------------------------------------------------------------!
CALL CPU_TIME(TIME1)
PRINT*
CALL PROGRESS(0)
OUTFILE='PROPAN.OUT'
OPEN(UNIT=20,FILE=OUTFILE,STATUS='REPLACE')
WRITE(20,*)
WRITE(20,*) 'PROPAN2020 v1.2'
WRITE(20,*)
OUTFILE='PROPAN.ERR'
OPEN(UNIT=30,FILE=OUTFILE,STATUS='REPLACE')
!-----------------------------------------------------------------------------------------------!
!    Read Input                                                                                 !
!-----------------------------------------------------------------------------------------------!
CALL INPUT(10) !Operational Conditions and Control Data
CALL INPUT(11) !Read Panel Coordinates
CALL INPUT(12) !Read Input Wake
CALL INPUT(13) !Read Field Points
!-----------------------------------------------------------------------------------------------!
!    Counters                                                                                   !
!-----------------------------------------------------------------------------------------------!
CALL COUNTERS
!-----------------------------------------------------------------------------------------------!
!    Initialize Variables                                                                       !
!-----------------------------------------------------------------------------------------------!
CALL INIVARS
!-----------------------------------------------------------------------------------------------!
!    Open Files for Storing Dipole and Source Matrices                                          !
!-----------------------------------------------------------------------------------------------!
OUTFILE='.DIJUNSTD.DAT'
OPEN(UNIT=21,FILE=OUTFILE,FORM='UNFORMATTED',DISP='DELETE',STATUS='NEW')
OUTFILE='.DIJLUSTD.DAT'
OPEN(UNIT=22,FILE=OUTFILE,FORM='UNFORMATTED',DISP='DELETE',STATUS='NEW')
OUTFILE='.DLKSTD.DAT'
OPEN(UNIT=23,FILE=OUTFILE,FORM='UNFORMATTED',DISP='DELETE',STATUS='NEW')
OUTFILE='.DLKLUSTD.DAT'
OPEN(UNIT=24,FILE=OUTFILE,FORM='UNFORMATTED',DISP='DELETE',STATUS='NEW')
OUTFILE='.DIJLUUNSTD.DAT'
OPEN(UNIT=25,FILE=OUTFILE,FORM='UNFORMATTED',DISP='DELETE',STATUS='NEW')
OUTFILE='.DLKUNSTD.DAT'
OPEN(UNIT=26,FILE=OUTFILE,FORM='UNFORMATTED',DISP='DELETE',STATUS='NEW')
OUTFILE='.DLKLUUNSTD.DAT'
OPEN(UNIT=27,FILE=OUTFILE,FORM='UNFORMATTED',DISP='DELETE',STATUS='NEW')
!-----------------------------------------------------------------------------------------------!
!    Loop on Inflow Coefficient (J or TSR)                                                      !
!-----------------------------------------------------------------------------------------------!
DO JJ=1,NU
   WRITE(20,*)
   IF (IROTOR == -1) WRITE(20,*) 'AoA=',UU(JJ)
   IF (IROTOR ==  0) WRITE(20,*) '  J=',UU(JJ)
   IF (IROTOR ==  1) WRITE(20,*) 'TSR=',UU(JJ)
!-----------------------------------------------------------------------------------------------!
!    Read Panel Coordinates                                                                     !
!-----------------------------------------------------------------------------------------------!
   IF ((JJ > 1).AND.(NWA > 0)) CALL INPUT(11)
!-----------------------------------------------------------------------------------------------!
!    Loop on Revolutions and Angular Step                                                       !
!-----------------------------------------------------------------------------------------------!
   DO TT=0,NT
      WRITE(20,*)
      WRITE(20,*) 'T=',TT
!-----------------------------------------------------------------------------------------------!
!    Starting Solution for Unsteady Calculation                                                 !
!-----------------------------------------------------------------------------------------------!
      IF ((JJ < 2).AND.(TT == 0).AND.(IFILESOL == 1)) THEN
         ! Calculate Geometry
         IF (IP       == 1) CALL BLADEGEOM
         IF (IP       == 1) CALL BLADEWAKEGEOM
         IF (IN       == 1) CALL NOZZLEGEOM
         IF (IN       == 1) CALL NOZZLEWAKEGEOM
         IF (IABS(IH) == 1) CALL HUBGEOM
         ! Numerical Differentiation Coefficients
         IF (IP       == 1) CALL NUMDIFBLADE
         IF ((IP == 1).AND.(NCPW /= 0)) CALL NUMDIFBLADEWAKE
         IF (IN       == 1) CALL NUMDIFNOZZLE
         IF (IABS(IH) == 1) CALL NUMDIFHUB
         ! Influence Coefficients
         IF (IP       == 1) CALL BLADECOEF
         IF (IP       == 1) CALL BLADEWAKECOEF
         IF (IN       == 1) CALL NOZZLECOEF
         IF (IN       == 1) CALL NOZZLEWAKECOEF
         IF (IABS(IH) == 1) CALL HUBCOEF
         ! Steady Matrices of Influence Coefficients
         CALL INFCOEFMATRIXSTD
         ! Unsteady Matrices of Influence Coefficients
         IF (NT > 0) CALL INFCOEFMATRIXUNSTD
         ! Read SOL.INP
         WRITE(20,*) 'Read Solution Files => SOL.INP, SOLW.INP'
         CALL INPUT(14)
         ! Calculate Forces
         CALL CTCQ(JJ,TT)
         ! Write Forces
         IF (IROTOR == -1) THEN
            IF (IP == 1) WRITE(20,*) 'Cl=',CTXP(TT),' Cd=',CQXP(TT)
         ELSEIF (IROTOR == 0) THEN
            IF (IP == 1) WRITE(20,*) 'KTP=',CTXP(TT),' KQP=',CQXP(TT)
            IF (IN == 1) WRITE(20,*) 'KTN=',CTXN(TT),' KQN=',CQXN(TT)
            IF (IH == 1) WRITE(20,*) 'KTH=',CTXH(TT),' KQH=',CQXH(TT)
         ELSEIF (IROTOR == 1) THEN
            IF (IP == 1) WRITE(20,*) 'CTT=',CTXP(TT),' CPT=',CQXP(TT)
            IF (IN == 1) WRITE(20,*) 'CTN=',CTXN(TT),' CPN=',CQXN(TT)
            IF (IH == 1) WRITE(20,*) 'CTH=',CTXH(TT),' CPH=',CQXH(TT)
         END IF !(IROTOR)
         ! Write Results
         CALL OUTPUT(31,JJ,TT,0) !Geometry
         CALL OUTPUT(32,JJ,TT,0) !Circulation
         CALL OUTPUT(33,JJ,TT,0) !Solution
         IF (IP       == 1) CALL OUTPUT(34,JJ,TT,0)
         IF (IN       == 1) CALL OUTPUT(35,JJ,TT,0)
         IF (IABS(IH) == 1) CALL OUTPUT(36,JJ,TT,0)
         IF (IFIELD   == 1) CALL OUTPUT(37,JJ,TT,0)
         IF (NWA       > 0) CALL OUTPUT(38,JJ,TT,0)
         IF (NCAV      > 0) CALL OUTPUT(40,JJ,TT,0)
      ELSE
!-----------------------------------------------------------------------------------------------!
!    Loop on Wake Alignment                                                                     !
!-----------------------------------------------------------------------------------------------!
         DO WW=0,NWA
            IF (NWA > 0) WRITE(20,*) 'Wake Align Iter=',WW
!-----------------------------------------------------------------------------------------------!
!    Wake Model                                                                                 !
!-----------------------------------------------------------------------------------------------!
            IF (WW > 0) THEN
               IF (IWA == 1) CALL WAKEALIGN1(JJ,TT) !Model 1
               IF (IWA == 2) CALL WAKEALIGN2(JJ,TT) !Model 2
            END IF !(WW > 0)
!-----------------------------------------------------------------------------------------------!
!    Calculate Geometry                                                                         !
!-----------------------------------------------------------------------------------------------!
            IF (((JJ < 2).AND.(TT == 0)).OR.(NWA > 0)) THEN
               IF (IP       == 1) CALL BLADEGEOM
               IF (IP       == 1) CALL BLADEWAKEGEOM
               IF (IN       == 1) CALL NOZZLEGEOM
               IF (IN       == 1) CALL NOZZLEWAKEGEOM
               IF (IABS(IH) == 1) CALL HUBGEOM
!-----------------------------------------------------------------------------------------------!
!    Compute Numerical Differentiation Coefficients                                             !
!-----------------------------------------------------------------------------------------------!
               IF (IP       == 1) CALL NUMDIFBLADE
               IF ((IP == 1).AND.(NCPW /= 0)) CALL NUMDIFBLADEWAKE
               IF (IN       == 1) CALL NUMDIFNOZZLE
               IF (IABS(IH) == 1) CALL NUMDIFHUB
!-----------------------------------------------------------------------------------------------!
!    Compute of Influence Coefficients                                                          !
!-----------------------------------------------------------------------------------------------!
               IF (IP       == 1) CALL BLADECOEF
               IF (IP       == 1) CALL BLADEWAKECOEF
               IF (IN       == 1) CALL NOZZLECOEF
               IF (IN       == 1) CALL NOZZLEWAKECOEF
               IF (IABS(IH) == 1) CALL HUBCOEF
!-----------------------------------------------------------------------------------------------!
!    Compute Steady Matrices of Influence Coefficients                                          !
!-----------------------------------------------------------------------------------------------!
               CALL INFCOEFMATRIXSTD
!-----------------------------------------------------------------------------------------------!
!    Compute Unsteady Matrices of Influence Coefficients                                        !
!-----------------------------------------------------------------------------------------------!
               IF (NT > 0) CALL INFCOEFMATRIXUNSTD
            END IF !(((JJ < 2).AND.(TT == 0)).OR.(NWA > 0))
!-----------------------------------------------------------------------------------------------!
!    Loop on Gap Model                                                                          !
!-----------------------------------------------------------------------------------------------!
            DCPG=0.D0
            DO II=0,NGAP
!-----------------------------------------------------------------------------------------------!
!    Right Hand Side                                                                            !
!-----------------------------------------------------------------------------------------------!
               IF (TT > 0) CALL SOLVEWAKE(TT,0) !Wake Model
               CALL SOLVERHSWET(JJ,II,TT)
!-----------------------------------------------------------------------------------------------!
!    Linear Kutta Condition                                                                     !
!-----------------------------------------------------------------------------------------------!
               CALL SOLVELKWET(JJ,TT)
!-----------------------------------------------------------------------------------------------!
!    Iterative Pressure Kutta Condition                                                         !
!-----------------------------------------------------------------------------------------------!
               IF (NKIT >= 1) CALL SOLVEIPKCWET(JJ,TT)
!-----------------------------------------------------------------------------------------------!
!    Reduced System Variables                                                                   !
!-----------------------------------------------------------------------------------------------!
               IF (IRED == 1) THEN
                  IF (IP == 1) THEN
                     POTPWET (:,:,TT)=POTP (:,:,TT)
                     POTPWWET(:,:,TT)=POTPW(:,:,TT)
                     IF (NCPW < 0) POTPWPWET(:,:,TT)=POTPWP(:,:,TT)
                     IF (NCPW > 0) POTPWSWET(:,:,TT)=POTPWS(:,:,TT)
                  END IF !(IP == 1)
                  IF (IN == 1) THEN
                     POTNWET (:,:,TT)=POTN (:,:,TT)
                     POTNWWET(:,:,TT)=POTNW(:,:,TT)
                  END IF !(IN == 1)
                  IF (IABS(IH) == 1) THEN
                     POTHWET (:,:,TT)=POTH (:,:,TT)
                  END IF !(IABS(IH) == 1)
               END IF !(IRED == 1)
!-----------------------------------------------------------------------------------------------!
!    Cavitation Check                                                                           !
!-----------------------------------------------------------------------------------------------!
               IF (NCAV > 0) CALL CAVCHECK(TT)
!-----------------------------------------------------------------------------------------------!
!    Initialise Cavitating Flow                                                                 !
!-----------------------------------------------------------------------------------------------!
               IF ((NPCAV /= 0).OR.(NSCAV /= 0)) THEN
                  CALL CAVPROP(TT)
                  CALL CAVERRC(0,ERRC)
                  CALL OUTPUT(39,JJ,TT,0)
                  CALL OUTPUT(41,JJ,TT,0)
!-----------------------------------------------------------------------------------------------!
!    Loop on Cavitation                                                                         !
!-----------------------------------------------------------------------------------------------!
                  DO CC=1,NCAV
                     WRITE(20,*) 'Cav Model Iter=',CC
!-----------------------------------------------------------------------------------------------!
!    Dynamic Boundary Condition                                                                 !
!-----------------------------------------------------------------------------------------------!
                     IF (NPCAV /= 0) CALL CAVPOTP(JJ,TT,CC)
                     IF (NSCAV /= 0) CALL CAVPOTS(JJ,TT,CC)
!-----------------------------------------------------------------------------------------------!
!    Right Hand Side                                                                            !
!-----------------------------------------------------------------------------------------------!
                     IF (TT > 0)CALL SOLVEWAKE(TT,CC) !Wake Model
                     IF (IRED == 0) CALL SOLVERHSCAV   (JJ,II,TT)
                     IF (IRED == 1) CALL SOLVERHSCAVRED(JJ,II,TT)
!-----------------------------------------------------------------------------------------------!
!    Linear Kutta Condition                                                                     !
!-----------------------------------------------------------------------------------------------!
                     IF (IRED == 0) CALL SOLVELKCAV(JJ,TT)
                     IF (IRED == 1) CALL SOLVELKWET(JJ,TT)
!-----------------------------------------------------------------------------------------------!
!    Iterative Pressure Kutta Condition                                                         !
!-----------------------------------------------------------------------------------------------!
                     IF (NKIT >= 1) THEN
                        IF (IRED == 0) CALL SOLVEIPKCCAV(JJ,TT)
                        IF (IRED == 1) CALL SOLVEIPKCWET(JJ,TT)
                     END IF !(NKIT >= 1)
!-----------------------------------------------------------------------------------------------!
!    Kinematic Boundary Condition                                                               !
!-----------------------------------------------------------------------------------------------!
                     IF (NPCAV /= 0) CALL CAVTHICKP(JJ,TT,CC)
                     IF (NSCAV /= 0) CALL CAVTHICKS(JJ,TT,CC)
!-----------------------------------------------------------------------------------------------!
!    Cavitation Properties                                                                      !
!-----------------------------------------------------------------------------------------------!
                     CALL CAVPROP(TT)
                     CALL OUTPUT(39,JJ,TT,CC)
                     CALL OUTPUT(41,JJ,TT,CC)
!-----------------------------------------------------------------------------------------------!
!    Progress                                                                                   !
!-----------------------------------------------------------------------------------------------!
                     IWORK=NINT(FLOAT(1+(JJ-1)*(NT+1)*(NWA+1)*(NGAP+1)*(NCAV+1)+ &
                           TT*(NWA+1)*(NGAP+1)*(NCAV+1)+WW*(NGAP+1)*(NCAV+1)+II*(NCAV+1)+CC)/ &
                           FLOAT(2+NU*(NT+1)*(NWA+1)*(NGAP+1)*(NCAV+1))*100.0)
                     IF (IWORK < 100) CALL PROGRESS(IWORK)
!-----------------------------------------------------------------------------------------------!
!    Convergence Criterion of Cavitation Model                                                  !
!-----------------------------------------------------------------------------------------------!
                     CALL CAVERRC(CC,ERRC)
                     IF (ERRC  < TOLC) GOTO 5000
!-----------------------------------------------------------------------------------------------!
!    End Loop on Cavitation                                                                     !
!-----------------------------------------------------------------------------------------------!
                  END DO !CC=1,NCAV
                  5000 CONTINUE
!-----------------------------------------------------------------------------------------------!
!    End Cavitating Flow                                                                        !
!-----------------------------------------------------------------------------------------------!
               END IF !((NPCAV /= 0).OR.(NSCAV /= 0))
               NPCAV=0 !Re-Initialise
               NSCAV=0
               NWCAV=0
!-----------------------------------------------------------------------------------------------!
!    Progress                                                                                   !
!-----------------------------------------------------------------------------------------------!
               IF (NCAV == 0) THEN
                  IWORK=NINT(FLOAT(1+(JJ-1)*(NT+1)*(NWA+1)*(NGAP+1)+TT*(NWA+1)*(NGAP+1)+ &
                        WW*(NGAP+1)+II)/FLOAT(2+NU*(NT+1)*(NWA+1)*(NGAP+1))*100.0)
                  IF (IWORK < 100) CALL PROGRESS(IWORK)
               END IF !(NCAV == 0)
!-----------------------------------------------------------------------------------------------!
!    Convergence Criterion of Gap Model                                                         !
!-----------------------------------------------------------------------------------------------!
               IF (ISTRIP == 1) THEN
                  CALL GAPERRG(TT,ERRG)
!-----------------------------------------------------------------------------------------------!
!    Write Convergence of Gap Model                                                             !
!-----------------------------------------------------------------------------------------------!
                  WRITE(20,*) 'Gap Model Iter=',II,' ERRG=',ERRG
                  IF (ERRG < TOLG) GOTO 2000
               END IF !(ISTRIP == 1)
!-----------------------------------------------------------------------------------------------!
!    End Loop on Gap Model                                                                      !
!-----------------------------------------------------------------------------------------------!
            END DO !II=0,NGAP
            2000 CONTINUE
!-----------------------------------------------------------------------------------------------!
!    Calculate Forces                                                                           !
!-----------------------------------------------------------------------------------------------!
            CALL CTCQ(JJ,TT)
!-----------------------------------------------------------------------------------------------!
!    Write Forces                                                                               !
!-----------------------------------------------------------------------------------------------!
            IF (IROTOR == -1) THEN
               IF (IP == 1) WRITE(20,*) 'Cl=',CTXP(TT),' Cd=',CQXP(TT)
            ELSEIF (IROTOR == 0) THEN
               IF (IP == 1) WRITE(20,*) 'KTP=',CTXP(TT),' KQP=',CQXP(TT)
               IF (IN == 1) WRITE(20,*) 'KTN=',CTXN(TT),' KQN=',CQXN(TT)
               IF (IH == 1) WRITE(20,*) 'KTH=',CTXH(TT),' KQH=',CQXH(TT)
            ELSEIF (IROTOR == 1) THEN
               IF (IP == 1) WRITE(20,*) 'CTT=',CTXP(TT),' CPT=',CQXP(TT)
               IF (IN == 1) WRITE(20,*) 'CTN=',CTXN(TT),' CPN=',CQXN(TT)
               IF (IH == 1) WRITE(20,*) 'CTH=',CTXH(TT),' CPH=',CQXH(TT)
            END IF !(IROTOR)
!-----------------------------------------------------------------------------------------------!
!    End Loop on Wake Alignment                                                                 !
!-----------------------------------------------------------------------------------------------!
         END DO !WW=0,NWA
!-----------------------------------------------------------------------------------------------!
!    Calculation of the Velocity Field                                                          !
!-----------------------------------------------------------------------------------------------!
         IF ((IFIELD == 1).AND.(TT == 0)) CALL   VELFSTD('ROTATING',JJ,TT,NFX,NFY, &
                                                                        XVF,YVF,ZVF,VFX,VFY,VFZ)
         IF ((IFIELD == 1).AND.(TT  > 0)) CALL VELFUNSTD('INERTIAL',JJ,TT,NFX,NFY, &
                                                                        XVF,YVF,ZVF,VFX,VFY,VFZ)
!-----------------------------------------------------------------------------------------------!
!    Calculation of the Pressure Field                                                          !
!-----------------------------------------------------------------------------------------------!
         IF ((IFIELD == 1).AND.(TT == 0)) CALL   PRESFSTD('ROTATING',JJ,TT,NFX,NFY, &
                                                                            XVF,YVF,ZVF,DPOTFDT)
         IF ((IFIELD == 1).AND.(TT  > 0)) CALL PRESFUNSTD('INERTIAL',JJ,TT,NFX,NFY, &
                                                                            XVF,YVF,ZVF,DPOTFDT)
!-----------------------------------------------------------------------------------------------!
!    Write Results                                                                              !
!-----------------------------------------------------------------------------------------------!
         CALL OUTPUT(31,JJ,TT,0) !Geometry
         CALL OUTPUT(32,JJ,TT,0) !Circulation
         CALL OUTPUT(33,JJ,TT,0) !Solution
         IF (IP       == 1) CALL OUTPUT(34,JJ,TT,0)
         IF (IN       == 1) CALL OUTPUT(35,JJ,TT,0)
         IF (IABS(IH) == 1) CALL OUTPUT(36,JJ,TT,0)
         IF (IFIELD   == 1) CALL OUTPUT(37,JJ,TT,0)
         IF (NWA       > 0) CALL OUTPUT(38,JJ,TT,0)
         IF (NCAV      > 0) CALL OUTPUT(40,JJ,TT,0)
!-----------------------------------------------------------------------------------------------!
!    Starting Solution for Unsteady Calculation                                                 !
!-----------------------------------------------------------------------------------------------!
      END IF
!-----------------------------------------------------------------------------------------------!
!    End Loop on Revolutions and Angular Step                                                   !
!-----------------------------------------------------------------------------------------------!
   END DO !TT=0,NT
!-----------------------------------------------------------------------------------------------!
!    Close Files                                                                                !
!-----------------------------------------------------------------------------------------------!
   INQUIRE(UNIT=39,OPENED=FCHECK)
   IF (FCHECK) CLOSE(UNIT=39)
   INQUIRE(UNIT=41,OPENED=FCHECK)
   IF (FCHECK) CLOSE(UNIT=41)
!-----------------------------------------------------------------------------------------------!
!    End Loop on Inflow Coefficient (J or TSR)                                                  !
!-----------------------------------------------------------------------------------------------!
END DO !JJ=1,NU
!-----------------------------------------------------------------------------------------------!
!    Close Files                                                                                !
!-----------------------------------------------------------------------------------------------!
CLOSE(UNIT=21)
CLOSE(UNIT=22)
CLOSE(UNIT=23)
CLOSE(UNIT=24)
CLOSE(UNIT=25)
CLOSE(UNIT=26)
CLOSE(UNIT=27)
!-----------------------------------------------------------------------------------------------!
!    Delete Variables                                                                           !
!-----------------------------------------------------------------------------------------------!
CALL DELVARS
CALL PROGRESS(100)
!-----------------------------------------------------------------------------------------------!
CALL CPU_TIME(TIME2)
TIME=TIME2-TIME1
WRITE( *,*)
WRITE(20,*)
IF (TIME < 60.D0) THEN
   WRITE( *,'(A,F4.1,A)')' Operation time = ',TIME,' seconds'
   WRITE(20,'(A,F4.1,A)')' Operation time = ',TIME,' seconds'
ELSEIF (TIME < 3600.D0) THEN
   TIME=TIME/60.D0
   WRITE( *,'(A,F4.0,A)')' Operation time = ',TIME,' minutes'
   WRITE(20,'(A,F4.0,A)')' Operation time = ',TIME,' minutes'
ELSE
   TIME=TIME/3600.D0
   WRITE( *,'(A,F4.0,A)')' Operation time = ',TIME,' hours'
   WRITE(20,'(A,F4.0,A)')' Operation time = ',TIME,' hours'
END IF
WRITE(*,*)
!-----------------------------------------------------------------------------------------------!
!    Close Output                                                                               !
!-----------------------------------------------------------------------------------------------!
CLOSE(UNIT=20)
CLOSE(UNIT=30)
!-----------------------------------------------------------------------------------------------!
END PROGRAM PROPAN
!-----------------------------------------------------------------------------------------------!
