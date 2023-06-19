!-----------------------------------------------------------------------------------------------!
SUBROUTINE GRAPE(NX,NY,XXG,YYG,ITMAX)
!-----------------------------------------------------------------------------------------------!
!    This routine generates the grid on the hub with fixed points on the boundary of a blade    !
!    sector.                                                                                    !
!                                                                                               !
!    The routine makes use of an unmodified version of the program grid2d (GRAPE) of Luís Eça.  !
!                                                                                               !
!    Adaptation by Falcao de Campos May 1999.                                                   !
!-----------------------------------------------------------------------------------------------!
!    This Program Generates a 2D Grid With an Elliptic System of Partial Differential Equations !
!                                                                                               !
!    The data in previous Tape 5 are defined as parameters internally                           !
!    TAPE 5 :                                                                                   !
!            ITMAX,ITYPE,IOUT                                                                   !
!            R(1),R(2),R(3)                                                                     !
!            SLC,ALFA,ILIM                                                                      !
!            IAK(1),IAK(2),JAK(1),JAK(2)                                                        !
!            IPLOT,IPLAN,XLOC                                                                   !
!            IPLI,IPLF,ISP,JPLI,JPLF,JSP (IPLOT NOT EQUAL TO 0)                                 !
!            BDFILE                                                                             !
!                                                                                               !
!            ITMAX  -> MAXIMUM NUMBER OF ITERATIONS                                             !
!            ITYPE  -> INITIAL APPROXIMATION OPTION                                             !
!                      0 = TRANSFINITE INTERPOLATION                                            !
!                      1 = LINEAR INTERPOLATION ON XI LINES                                     !
!                      2 = LINEAR INTERPOLATION ON ETA LINES                                    !
!                      3 = BILINEAR INTERPOLATION                                               !
!                      4 = TRANSFINITE INTERPOLATION                                            !
!            IOUT   -> GRID ANALYSIS OPTION                                                     !
!            R(1)   -> TOLERANCE FOR |X(N+1)-X(N)|                                              !
!            R(2)   -> TOLERANCE FOR |Y(N+1)-Y(N)|                                              !
!            R(3)   -> UNDERRELAXATION PARAMETER FOR THE SOLVER                                 !
!            SLC    -> 'TIME-STEP' FACTOR FOR THE MAIN DIAGONAL                                 !
!            ALFA   -> UNDERRELAXATION PARAMETER FOR THE                                        !
!                             NON-LINEAR PART OF THE CONTROL FUNCTIONS                          !
!            ILIM   -> NUMBER OF ITERATIONS SKIPPED FOR THE CALCULATION                         !
!                      OF THE CONTROL FUNCTIONS                                                 !
!            IAK(1) -> NUMBER OF GRID LINES AWAY FROM THE BOUNDARY I=1                          !
!                      WHERE THE GRID SHOULD REMAIN ORTHOGONAL TO IT                            !
!            IAK(2) -> NUMBER OF GRID LINES AWAY FROM THE BOUNDARY I=NX                         !
!                      WHERE THE GRID SHOULD REMAIN ORTHOGONAL TO IT                            !
!            JAK(1) -> NUMBER OF GRID LINES AWAY FROM THE BOUNDARY J=1                          !
!                      WHERE THE GRID SHOULD REMAIN ORTHOGONAL TO IT                            !
!            JAK(2) -> NUMBER OF GRID LINES AWAY FROM THE BOUNDARY J=NY                         !
!                      WHERE THE GRID SHOULD REMAIN ORTHOGONAL TO IT                            !
!            IPLOT  -> OUPUT FOR PLOTTING                                                       !
!                      0 = NO OUTPUT                                                            !
!                      1 = ONLY FINAL GRID                                                      !
!                      2 = INITIAL AND FINAL GRID                                               !
!                      2 = ONLY FINAL GRID                                                      !
!                                                                                               !
!            IPLAN  -> COORDINATE THAT REMAINS CONSTANT FOR A 3-D GRID                          !
!                      0 = NO 3-D GRID                                                          !
!                      1 = CONSTANT X                                                           !
!                      2 = CONSTANT Y                                                           !
!                      3 = CONSTANT Z                                                           !
!            XLOC   -> COORDINATE OF THE CONSTANT VARIABLE                                      !
!            IPLI   -> INITIAL XI LINE FOR PLOTTING                                             !
!            IPLF   -> FINAL XI LINE FOR PLOTTING                                               !
!            ISP    -> XI LINE SKIP FOR PLOTTING                                                !
!            JPLI   -> INITIAL ETA LINE FOR PLOTTING                                            !
!            JPLF   -> FINAL ETA LINE FOR PLOTTING                                              !
!            JSP    -> ETA LINE SKIP FOR PLOTTING                                               !
!            BDFILE -> NAME OF THE INPUT FILE OF THE COORDINATES OF THE BOUNDARY NODES          !
!-----------------------------------------------------------------------------------------------!
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!-----------------------------------------------------------------------------------------------!
PARAMETER(NDIM=300)
!-----------------------------------------------------------------------------------------------!
CHARACTER *40 IDGEOM
INTEGER :: I,J,NX,NY,LBOUN(4)
DOUBLE PRECISION :: XXG(NX,NY),YYG(NX,NY),R(12)
!-----------------------------------------------------------------------------------------------!
COMMON /COOR/  X(0:NDIM+1,0:NDIM+1),Y(0:NDIM+1,0:NDIM+1)
COMMON /FORCE/ RXI(NDIM,NDIM), RETA(NDIM,NDIM)
COMMON /TRANS/ FXI(NDIM,NDIM),FETA(NDIM,NDIM)
COMMON /AREA/  DSI(NDIM,4),DST(4),IDSIT(4),IAK(2),JAK(2)
COMMON /BOUN/  XL(NDIM,4),YL(NDIM,4),DXDL(NDIM,4),DYDL(NDIM,4),RL(NDIM,4)
!-----------------------------------------------------------------------------------------------!
DIMENSION DS(8*NDIM),RIS(8*NDIM),DSG(NDIM)
!-----------------------------------------------------------------------------------------------!
DATA LBOUN/4*0/ 
!-----------------------------------------------------------------------------------------------!
!    Define Parameters                                                                          !
!-----------------------------------------------------------------------------------------------!
idgeom='Generation of the grid with GRAPE'
!ITMAX =5000 !251
ITYPE =3
IOUT  =0
R(1)  =1.D-6 !0.000001
R(2)  =1.D-6 !0.000001
R(3)  =1.
IAK(1)=0
IAK(2)=0
JAK(1)=0
JAK(2)=0
IPLOT =0
IPLAN =1
XLOC  =1.0
IPLI  =0
IPLF  =0
ISP   =0
JPLI  =0
JPLF  =0
JSP   =0
ILIM  =0 !1
DO I=1,4
   LBOUN(I)=0
END DO !I=1,4
!-----------------------------------------------------------------------------------------------!
!    Write Input Parameters                                                                     !
!-----------------------------------------------------------------------------------------------!
!!WRITE(6,1000) IDGEOM
!!WRITE(6,1010) ITMAX,ITYPE,IOUT
!!WRITE(6,1020) R(1),R(2),R(3)
!!WRITE(6,1030) SLC,ALFA,ILIM
!!WRITE(6,1040) IAK(1),IAK(2),JAK(1),JAK(2)
!!WRITE(6,1050) IPLOT,IPLAN,XLOC
IF (IPLOT /= 0) WRITE(6,1060) IPLI,IPLF,ISP,JPLI,JPLF,JSP
!!WRITE(6,1075) (LBOUN(I),I=1,4)
!-----------------------------------------------------------------------------------------------!
!    Define Boundary Points                                                                     !
!-----------------------------------------------------------------------------------------------!
DO I=1,NX
   X(1 ,I)=XXG(I,1 )
   Y(1 ,I)=YYG(I,1 )
   X(NY,I)=XXG(I,NY)
   Y(NY,I)=YYG(I,NY)
END DO !I=1,NX
!-----------------------------------------------------------------------------------------------!
DO J=1,NY
   X(J,1 )=XXG(1 ,J)
   Y(J,1 )=YYG(1 ,J)
   X(J,NX)=XXG(NX,J)
   Y(J,NX)=YYG(NX,J)
END DO !J=1,NY
!-----------------------------------------------------------------------------------------------!
DO I=1,4
   IDSIT(I)=0
END DO !I=1,4
!-----------------------------------------------------------------------------------------------!
IF (IPLOT /= 0) THEN
   IPLI=MAX(0   ,IPLI)
   IPLF=MIN(NX+1,IPLF)
   JPLI=MAX(0   ,JPLI)
   JPLF=MIN(NY+1,JPLF)
   IF (ISP == 0) THEN
      IPLI=1
      IPLF=NX
      ISP =1
   END IF !(ISP == 0)
   IF (JSP == 0) THEN
      JPLI=1
      JPLF=NY
      JSP=1
   END IF !(JSP == 0)
END IF !(IPLOT /= 0)
!-----------------------------------------------------------------------------------------------!
!    Calculate Splines for Moving Nodes in Boundaries                                           !
!-----------------------------------------------------------------------------------------------!
DO I=1,4
   IF (LBOUN(I) >= 3) CALL SPLIN(NX,NY,I)
END DO !I=1,4
!-----------------------------------------------------------------------------------------------!
!    Calculate Initial Guess                                                                    !
!-----------------------------------------------------------------------------------------------!
!!PRINT*, ' Calculating initial grid ... '
!!PRINT*
CALL GUESSA(NX,NY,ITYPE)
!-----------------------------------------------------------------------------------------------!
!    Iterative Solution                                                                         !
!-----------------------------------------------------------------------------------------------!
!!PRINT*, ' Solving equations ... '
!!PRINT*
CALL SIP(R,SLC,ALFA,LBOUN,NX,NY,ITMAX,ILIM,IOUT)
!-----------------------------------------------------------------------------------------------!
!    Fill Extra Grid Lines                                                                      !
!-----------------------------------------------------------------------------------------------!
CALL BORD(NX,NY)
!-----------------------------------------------------------------------------------------------!
LT=LBOUN(1)+LBOUN(2)+LBOUN(3)+LBOUN(4)
!-----------------------------------------------------------------------------------------------!
!    Grid Coordinates                                                                           !
!-----------------------------------------------------------------------------------------------!
!XXG=X
!YYG=Y
DO J=1,NY
   DO I=1,NX
      XXG(I,J)=X(J,I)
      YYG(I,J)=Y(J,I)
   END DO !I=1,NX
END DO !J=1,NY
!-----------------------------------------------------------------------------------------------!
!    Calculate Grid Properties                                                                  !
!-----------------------------------------------------------------------------------------------!
!!PRINT*, ' Calculating grid properties ...'
!!PRINT*
!-----------------------------------------------------------------------------------------------!
CALL ANGRI(NX,NY,1)
!-----------------------------------------------------------------------------------------------!
1000 FORMAT(' IDGEOM                      :',A40)
1010 FORMAT(' ITMAX,ITYPE,IOUT            :',3I5)
1020 FORMAT(' R(1),R(2),R(3)              :',3F12.9)
1030 FORMAT(' SLC,ALFA,ILIM               :',2F7.3,I5)
1040 FORMAT(' IAK(1),IAK(2),JAK(1),JAK(2) :',4I5)
1050 FORMAT(' IPLOT,IPLAN,XLOC            :',2I5,F10.4)
1060 FORMAT(' IPLI,IPLF,ISP,JPLI,JPLF,JSP :',6I5)
1070 FORMAT(' BDFILE                      :',A32)
1075 FORMAT(' LBOUN                       :',4I5)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE GRAPE
!-----------------------------------------------------------------------------------------------!
