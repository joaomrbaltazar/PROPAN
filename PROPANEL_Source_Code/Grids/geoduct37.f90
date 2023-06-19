!-----------------------------------------------------------------------------------------------!
!    Geometry of Duct 37                                                                        !
!-----------------------------------------------------------------------------------------------!
SUBROUTINE GEODUCT37(IS,N,T,Y)
!-----------------------------------------------------------------------------------------------!
!    Created by: J.A.C. Falcao de Campos, IST, September 2007                                   !
!    Modified  : 05072017, J. Baltazar, 2017 version 1.0                                        !
!-----------------------------------------------------------------------------------------------!
!    Input Description                                                                          !
!                                                                                               !
!    IS        =1  : Outer surface                                                              !
!              =2  : Inner surface                                                              !
!    N         =     Number of points on the inner or outer surface                             !
!    T(N)      =     Parameter for stretching on x coordinate                                   !
!                                                                                               !
!    Output Description                                                                         !
!                                                                                               !
!    X(N)      =     Array of X ordinates (REMOVED)                                             !
!    Y(N)      =     Array of Y ordinates                                                       !
!-----------------------------------------------------------------------------------------------!
!    Declarations                                                                               !
!-----------------------------------------------------------------------------------------------!
!*USE IMSLF90     ! Visual FORTRAN 6
IMPLICIT NONE
INTEGER :: I,J
DOUBLE PRECISION :: PI
DOUBLE PRECISION :: ALFA,COSALFA,SINALFA,R1,R2,H1,H2,X1,X2,Y1,Y2
DOUBLE PRECISION :: XS,F,G,PPVALU !,DLEFT,DRIGHT
!-----------------------------------------------------------------------------------------------!
!    Variables                                                                                  !
!-----------------------------------------------------------------------------------------------!
INTEGER :: N,IS !,ILEFT,IRIGHT
INTEGER,PARAMETER :: N1=8,N2=10
!*INTEGER,PARAMETER :: N1=11,N2=12
!*INTEGER,PARAMETER :: N1=11,N2=11
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION,DIMENSION(N)    :: T,X,Y
!*DOUBLE PRECISION,DIMENSION(N)    :: DY,D2Y
DOUBLE PRECISION,DIMENSION(N2)   :: XIN2,YIN2,BREAK2
DOUBLE PRECISION,DIMENSION(4,N2) :: CSCOEF2
DOUBLE PRECISION,DIMENSION(N1)   :: XIN1,YIN1,BREAK1
DOUBLE PRECISION,DIMENSION(4,N1) :: CSCOEF1
XIN2 = (/0.D0,0.0125D0,0.025D0,0.05D0,0.075D0,0.1D0,0.15D0,0.2D0,0.3D0,0.4D0/)
YIN2 = (/0.1833D0,0.1461D0,0.1283D0,0.1D0,0.0792D0,0.0625D0,0.0383D0,0.0208D0,0.0033D0,0.D0/)
!*XIN2 = (/0.D0,0.0125D0,0.025D0,0.05D0,0.075D0,0.1D0,0.15D0,0.2D0,0.25D0,0.3D0,0.35D0,0.4D0/)
!*YIN2 = (/0.1833D0,0.1461D0,0.1283D0,0.1D0,0.0792D0,0.0625D0,0.0383D0,0.0208D0,0.0094D0,0.0033D0,0.0008D0,0.D0/)
!*XIN2 = (/0.D0,0.025D0,0.05D0,0.075D0,0.1D0,0.15D0,0.2D0,0.25D0,0.3D0,0.35D0,0.4D0/)
!*YIN2 = (/0.1833D0,0.1283D0,0.1D0,0.0792D0,0.0625D0,0.0383D0,0.0208D0,0.0094D0,0.0033D0,0.0008D0,0.D0/)
XIN1 = (/0.6D0,0.7D0,0.8D0,0.85D0,0.9D0,0.95D0,0.975D0,1.D0/)
YIN1 = (/0.D0,0.0008D0,0.01D0,0.0208D0,0.0383D0,0.065D0,0.0833D0,0.1242D0/)
BREAK2 =0.D0
BREAK1 =0.D0
CSCOEF2=0.D0
CSCOEF1=0.D0
!-----------------------------------------------------------------------------------------------!
PI=4.D0*DATAN(1.D0)
H1=12.42D0/100.D0          ! Ordinate of duct trailing edge / chord
H2=18.33D0/100.D0          ! Ordinate of duct leading edge  / chord
!-----------------------------------------------------------------------------------------------!
X=0.D0
Y=0.D0
!*DY=0.D0
!*D2Y=0.D0
!-----------------------------------------------------------------------------------------------!
!    Outer Surface                                                                              !
!-----------------------------------------------------------------------------------------------!
IF (IS == 1) THEN
   ALFA=3.38D0*PI/180.D0     ! Inclination angle of straight part on outer surface
   COSALFA=DCOS(ALFA)
   SINALFA=DSIN(ALFA)
   R1=3.34D0/100.D0          ! Curvature radius of duct trailing edge / chord
   R2=3.34D0/100.D0          ! Curvature radius of duct leading edge  / chord
   X1=1.D0-R1*(1.D0-SINALFA) ! X coordinate of transition point from straight to trailing edge   
   X2=R2*(1.D0+SINALFA)      ! X coordinate of transition point from straight to leading edge   
   Y1=H1+R1*COSALFA          ! Y coordinate of transition point from straight to trailing edge   
   Y2=H2+R2*COSALFA          ! Y coordinate of transition point from straight to leading edge   
!-----------------------------------------------------------------------------------------------!
!    Generate Points
!-----------------------------------------------------------------------------------------------!
   DO I=1,N  ! Do loop on the number of points
      X(I)=1.D0-T(I) ! Abcissae as parameter from te to le
!-----------------------------------------------------------------------------------------------!
!    Leading Edge                                                                               !
!-----------------------------------------------------------------------------------------------!
      IF (X(I).LE.X2) THEN
         XS=X(I)/R2
         Y(I)=H2+R2*DSQRT(1.D0-(XS-1.D0)**2)  ! Circle
!-----------------------------------------------------------------------------------------------!
!    Trailing Edge                                                                              !
!-----------------------------------------------------------------------------------------------!
      ELSEIF (X(I).GE.X1) THEN
         XS=(1.D0-X(I))/R1
         Y(I)=H1+R1*DSQRT(1.D0-(1.D0-XS)**2)  ! Circle
!-----------------------------------------------------------------------------------------------!
!    Straight Part                                                                              !
!-----------------------------------------------------------------------------------------------!
      ELSE 
         XS=(X(I)-X2)/(X1-X2)
         Y(I)=Y1*XS+Y2*(1.D0-XS)              ! Linear interpolation
      ENDIF
   END DO !I=1,N
!-----------------------------------------------------------------------------------------------!
!    Inner Surface                                                                              !
!-----------------------------------------------------------------------------------------------!
ELSEIF (IS == 2) THEN
   X1=0.6D0
   X2=0.4D0
   DO I=1,N2
      XS=XIN2(I)/X2                          ! Normalized coordinate
      G=H2*(1.D0-DSQRT(1.D0-(XS-1.D0)**2))   ! Reference ellipse
      YIN2(I)=YIN2(I)-G                      ! Deviation from reference ellipse
   END DO !I=1,N2
!* ILEFT=1                                   ! Left 1st derivative
!* IRIGHT=1                                  ! Right 1st derivative
!* DLEFT=0.D0                                ! put to zero
!* DRIGHT=0.D0                               ! put to zero
!* Call DCSDEC(N2,XIN2,YIN2,ILEFT,DLEFT,IRIGHT,DRIGHT,BREAK2,CSCOEF2)    ! Spline interpolant
   BREAK2(:)=XIN2(:)
   CSCOEF2(1,:)=YIN2(:)
   CALL CUBSPL(BREAK2,CSCOEF2,N2,1,1)          ! Spline interpolant
   DO I=1,N1
      XS=(XIN1(I)-X1)/(1.D0-X1)              ! Normalized coordinate
      G=H1*(1.D0-DSQRT(1.D0-XS**2))          ! Reference ellipse 
      YIN1(I)=YIN1(I)-G                      ! Deviation from reference ellipse
   END DO !I=1,N1
!* ILEFT=1                                   ! Left 1st derivative
!* IRIGHT=1                                  ! Right 1st derivative
!* DLEFT=0.D0                                ! put to zero
!* DRIGHT=0.D0                               ! put to zero
!* Call DCSDEC(N1,XIN1,YIN1,ILEFT,DLEFT,IRIGHT,DRIGHT,BREAK1,CSCOEF1)    ! Spline interpolant
   BREAK1(:)=XIN1(:)
   CSCOEF1(1,:)=YIN1(:)
   CALL CUBSPL(BREAK1,CSCOEF1,N1,1,1)          ! Spline interpolant
!-----------------------------------------------------------------------------------------------!
!    Generate Points
!-----------------------------------------------------------------------------------------------!
   DO I=1,N  ! Do loop on the number of points
      X(I)=T(I) ! Abcissae as parameter
!-----------------------------------------------------------------------------------------------!
!    Leading Region                                                                             !
!-----------------------------------------------------------------------------------------------!
      IF (X(I) .LE. X2 ) THEN
         XS=X(I)/X2
         G=H2*(1.D0-DSQRT(1.D0-(XS-1.D0)**2))    ! Reference ellipse
!*       F=DCSVAL(X(I),N2-1,BREAK2,CSCOEF2)      ! Evaluate the spline
         F=PPVALU(BREAK2,CSCOEF2,N2-1,4,X(I),0)  ! Evaluate the spline
!*       DY(I)=DCSDER(1,X(I),N2-1,BREAK2,CSCOEF2)! Evaluate the spline 1st derivative
!*       D2Y(I)=DCSDER(2,X(I),N2-1,BREAK2,CSCOEF2)! Evaluate the spline 2nd derivative
!-----------------------------------------------------------------------------------------------!
!    Trailing Region                                                                            !
!-----------------------------------------------------------------------------------------------!
      ELSEIF (X(I) .GE. X1 ) THEN
         XS=(X(I)-X1)/(1.D0-X1)
         G=H1*(1.D0-DSQRT(1.D0-XS**2))           ! Reference ellipse 
!*       F=DCSVAL(X(I),N1-1,BREAK1,CSCOEF1)      ! Evaluate the spline
         F=PPVALU(BREAK1,CSCOEF1,N1-1,4,X(I),0)  ! Evaluate the spline
!*       DY(I)=DCSDER(1,X(I),N2-1,BREAK1,CSCOEF1)! Evaluate the spline 1st derivative
!*       D2Y(I)=DCSDER(2,X(I),N2-1,BREAK1,CSCOEF1)! Evaluate the spline 2nd derivative
!-----------------------------------------------------------------------------------------------!
!    Straight Part                                                                              !
!-----------------------------------------------------------------------------------------------!
      ELSE 
         F=0.D0
         G=0.D0           ! Linear interpolation
!*       DY(I)=0.D0
!*       D2Y(I)=0.D0
      ENDIF
      Y(I)=F+G
   END DO !I=1,N
!-----------------------------------------------------------------------------------------------!
END IF !(IS)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE GEODUCT37
!-----------------------------------------------------------------------------------------------!
