!-----------------------------------------------------------------------------------------------!
!    POST-PROCESSOR FOR THE IST SURFACE PANEL METHOD CODE FOR ROTOR ANALYSIS                    !
!                                                                                               !
!    Created by: J.A.C. Falcao de Campos, IST                                                   !
!    Modified  : 23052014, J. Baltazar, 2014 version 1.0                                        !
!    Modified  : 05062014, J. Baltazar, 2014 version 1.1                                        !
!    Modified  : 02072014, J. Baltazar, 2014 version 1.2                                        !
!    Modified  : 23012015, J. Baltazar, 2015 version 1.0 Suction Force Correction               !
!    Modified  : 10032015, J. Baltazar, 2015 version 1.1 Unsteady Flow                          !
!    Modified  : 23112015, J. Baltazar, 2015 version 1.2                                        !
!    Modified  : 06072016, J. Baltazar, 2016 version 1.0                                        !
!    Modified  : 07072017, J. Baltazar, 2017 version 1.0                                        !
!    Modified  : 20032019, J. Baltazar, 2019 version 1.0 Harmonic Analysis                      !
!    Modified  : 20042020, J. Baltazar, 2020 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
PROGRAM PROPOST
!-----------------------------------------------------------------------------------------------!
!    Input description:                                                                         !
!                                                                                               !
!    The input is made from two files:                                                          !
!                                                                                               !
!    PROPOST.INP   : Rotor operational conditions and control variables                         !
!                                                                                               !
!    CL_CD_INVL.DAT: Inviscid 2D CL and CD data                                                 !
!                                                                                               !
!    CL_CD_VISC.DAT: Viscous 2D CL and CD data                                                  !
!-----------------------------------------------------------------------------------------------!
!    Input variables in PROPOST.INP                                                             !
!                                                                                               !
!    COMMENT(5)=    Comments for run description (5 lines)                                      !
!                                                                                               !
!    IP        = 1: blade and blade wake panelling (0 no panelling)                             !
!    IN        = 1: nozzle and nozzle wake panelling (0 no panelling)                           !
!    IH        = 1: hub panelling (0 no panelling)                                              !
!              =-1: no-hub but closure of blade root with panels                                !
!    NB        =    Number of Blades                                                            !
!                                                                                               !
!    ISTRIP    =    Input option to exclude last tip panel strip                                !
!                   ISTRIP=0: zero and finite tip (default value)                               !
!                   ISTRIP=1: include last panel strip for gap model                            !
!                   ISTRIP=2: include last panel strip for closed tip                           !
!                                                                                               !
!    IROTOR    =    Choice on the rotor                                                         !
!                   IROTOR=0: propeller                                                         !
!                   IROTOR=1: turbine                                                           !
!    NU        =    Number of inflow conditions for the analysis                                !
!    IDENTU(NU)=    Inflow condition identification (2 characters)                              !
!    UU(NU)    =    Advance coefficient J=U/nD (IF IROTOR=0)                                    !
!                   tip-speed-ratio TSR=OmegaR/U (IF IROTOR=1)                                  !
!                                                                                               !
!    NWA       =    Number of wake alignment iterations                                         !
!                                                                                               !
!    NCPW      =    Number of Streamwise blade wake panels for super-cavitation                 !
!                   NCPW > 0: Super-cavitation on upper side                                    !
!                   NCPW < 0: Super-cavitation on lower side                                    !
!                                                                                               !
!    NREV      =    Number of revolutions                                                       !
!                                                                                               !
!    NTETA     =    Number of theta steps                                                       !
!                                                                                               !
!    NNHAR     =    Number of duct harmonics                                                    !
!                                                                                               !
!    INTER     =    Choice on the interpolations                                                !
!                   INTER=0: Linear interpolation                                               !
!                   INTER=1: Quadratic interpolation                                            !
!                   INTER=2: Cubic interpolation                                                !
!                                                                                               !
!    NRS       =    Number of selected radii                                                    !
!    RS(I)     =    Selected radii                                                              !
!                                                                                               !
!    NXS       =    Number of selected wake stations                                            !
!    XS(I)     =    Selected wake stations                                                      !
!                                                                                               !
!    NTS       =    Number of selected theta                                                    !
!    TS(I)     =    Selected theta                                                              !
!                                                                                               !
!    JI        =    Initial strip of the blade wake                                             !
!    JF        =    Last strip of the blade wake                                                !
!                                                                                               !
!    ISF       =    Suction Force Correction                                                    !
!                   ISF=0: no suction force correction                                          !
!                   ISF=1: Suction force correction                                             !
!                                                                                               !
!    RE        =    Reynolds Number (Re=rho*Omega*R**2/mu)                                      !
!-----------------------------------------------------------------------------------------------!
!    Declarations                                                                               !
!-----------------------------------------------------------------------------------------------!
USE PROPOST_MOD
IMPLICIT NONE
CHARACTER OUTFILE*50
INTEGER :: JJ,IWORK
REAL    :: TIME,TIME1,TIME2
!-----------------------------------------------------------------------------------------------!
NAMELIST /PROJ_INPUT/ COMMENT,IP,IN,IH,NB,JI,JF,ISTRIP,IROTOR,NU,IDENTU,UU,NWA,NCPW,NREV,NTETA, &
                      NNHAR,INTER,NRS,RS,NXS,XS,NTS,TS,ISF,RE
!-----------------------------------------------------------------------------------------------!
CALL CPU_TIME(TIME1)
PRINT*
CALL PROGRESS(0)
!-----------------------------------------------------------------------------------------------!
!    Read Input Operational Conditions and Control Data                                         !
!-----------------------------------------------------------------------------------------------!
OPEN(UNIT=10,FILE='PROPOST.INP',STATUS='UNKNOWN')
READ(10,PROJ_INPUT)
CLOSE(UNIT=10)
!-----------------------------------------------------------------------------------------------!
!    Counters                                                                                   !
!-----------------------------------------------------------------------------------------------!
NT=NREV*NTETA
P =DFLOAT(NTETA)/DFLOAT(NB)
!-----------------------------------------------------------------------------------------------!
!    Loop on Inflow Coefficient (J or TSR)                                                      !
!-----------------------------------------------------------------------------------------------!
DO JJ=1,NU
!-----------------------------------------------------------------------------------------------!
!    Read Panel Coordinates                                                                     !
!-----------------------------------------------------------------------------------------------!
   CALL INPUT(11,JJ)
!-----------------------------------------------------------------------------------------------!
!    Counters                                                                                   !
!-----------------------------------------------------------------------------------------------!
   NCP =IP*(NCP1-1)
   NC  =IP*(NCP/2)
   NC1 =IP*(NC+1)
   NPW =IP*(NPW1-1)
   NRW =IP*(NRW1-1)
   NNX =IN*((NNXT-1)/2)
   NNX1=IN*(NNX+1)
!-----------------------------------------------------------------------------------------------!
!    Read Input                                                                                 !
!-----------------------------------------------------------------------------------------------!
   CALL INPUT(21,JJ)
   CALL INPUT(22,JJ)
   CALL INPUT(23,JJ)
!-----------------------------------------------------------------------------------------------!
!    Plot 2D Pressure Distribution                                                              !
!-----------------------------------------------------------------------------------------------!
   CALL PLOTCP(JJ)
   IWORK=NINT((FLOAT(JJ-1)*7.0+1.0)/(FLOAT(NU)*7.0)*100.0)
   CALL PROGRESS(IWORK)
!-----------------------------------------------------------------------------------------------!
!    Plot Harmonic Analysis                                                                     !
!-----------------------------------------------------------------------------------------------!
   CALL PLOTHARM(JJ)
   IWORK=NINT((FLOAT(JJ-1)*7.0+2.0)/(FLOAT(NU)*7.0)*100.0)
   CALL PROGRESS(IWORK)
!-----------------------------------------------------------------------------------------------!
!    Plot Wake Geometry                                                                         !
!-----------------------------------------------------------------------------------------------!
   CALL PLOTWAKE(JJ)
   IWORK=NINT((FLOAT(JJ-1)*7.0+3.0)/(FLOAT(NU)*7.0)*100.0)
   CALL PROGRESS(IWORK)
!-----------------------------------------------------------------------------------------------!
!    Plot Solution                                                                              !
!-----------------------------------------------------------------------------------------------!
   CALL PLOTSOL(JJ)
   IWORK=NINT((FLOAT(JJ-1)*7.0+4.0)/(FLOAT(NU)*7.0)*100.0)
   CALL PROGRESS(IWORK)
!-----------------------------------------------------------------------------------------------!
!    Plot Cavitation                                                                            !
!-----------------------------------------------------------------------------------------------!
   CALL PLOTCAV(JJ)
   IWORK=NINT((FLOAT(JJ-1)*7.0+5.0)/(FLOAT(NU)*7.0)*100.0)
   CALL PROGRESS(IWORK)
!-----------------------------------------------------------------------------------------------!
!    Compute Viscous Corrections                                                                !
!-----------------------------------------------------------------------------------------------!
   CALL VISCOR(JJ)
   IWORK=NINT((FLOAT(JJ-1)*7.0+6.0)/(FLOAT(NU)*7.0)*100.0)
   CALL PROGRESS(IWORK)
!-----------------------------------------------------------------------------------------------!
!    Delete Variables                                                                           !
!-----------------------------------------------------------------------------------------------!
   CALL DELVARS
   IWORK=NINT((FLOAT(JJ-1)*7.0+7.0)/(FLOAT(NU)*7.0)*100.0)
   CALL PROGRESS(IWORK)
!-----------------------------------------------------------------------------------------------!
!    End Loop on Inflow Coefficient (J or TSR)                                                  !
!-----------------------------------------------------------------------------------------------!
END DO !JJ=1,NU
!-----------------------------------------------------------------------------------------------!
CALL CPU_TIME(TIME2)
TIME=TIME2-TIME1
PRINT*
IF (TIME < 60.D0) THEN
   WRITE(*,'(A,F4.1,A)')' Operation time = ',TIME,' seconds'
ELSEIF (TIME < 3600.D0) THEN
   TIME=TIME/60.D0
   WRITE(*,'(A,F4.0,A)')' Operation time = ',TIME,' minutes'
ELSE
   TIME=TIME/3600.D0
   WRITE(*,'(A,F4.0,A)')' Operation time = ',TIME,' hours'
END IF
PRINT*
!-----------------------------------------------------------------------------------------------!
END PROGRAM PROPOST
!-----------------------------------------------------------------------------------------------!
