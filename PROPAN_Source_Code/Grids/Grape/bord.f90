!-----------------------------------------------------------------------------------------------!
SUBROUTINE BORD(NX,NY)
!-----------------------------------------------------------------------------------------------!
!    Created by   : L. Eca, IST                                                                 !
!    Last Revision: Joao Baltazar, IST, April 2005                                              !
!-----------------------------------------------------------------------------------------------!
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!-----------------------------------------------------------------------------------------------!
PARAMETER(NDIM=150)
!-----------------------------------------------------------------------------------------------!
COMMON /COOR/ X(0:NDIM+1,0:NDIM+1),Y(0:NDIM+1,0:NDIM+1)
!-----------------------------------------------------------------------------------------------!
NXP1=NX+1
NYP1=NY+1
NXM1=NX-1
NYM1=NY-1
ZERO=1.D-15
!-----------------------------------------------------------------------------------------------!
!    Corners                                                                                    !
!-----------------------------------------------------------------------------------------------!
DXDKS =X(1,2)-X(1,1)
DYDKS =Y(1,2)-Y(1,1)
DXDET =X(2,1)-X(1,1)
DYDET =Y(2,1)-Y(1,1)
DSDET =DXDET**2+DYDET**2
DSDKS =DXDKS**2+DYDKS**2
DSDKS =DMAX1(DSDKS,ZERO)
FH    =2.D0*DSQRT(DSDET/DSDKS)
X(0,1)=X(2,1)+DYDKS*FH
Y(0,1)=Y(2,1)-DXDKS*FH
X(1,0)=X(1,2)-(Y(2,1)-Y(0,1))/FH
Y(1,0)=Y(1,2)+(X(2,1)-X(0,1))/FH
!-----------------------------------------------------------------------------------------------!
DXDKS    =X(1,NX)-X(1,NXM1)
DYDKS    =Y(1,NX)-Y(1,NXM1)
DXDET    =X(2,NX)-X(1,NX)
DYDET    =Y(2,NX)-Y(1,NX)
DSDET    =DXDET**2+DYDET**2
DSDKS    =DXDKS**2+DYDKS**2
DSDKS    =DMAX1(DSDKS,ZERO)
FH       =2.D0*DSQRT(DSDET/DSDKS)
X(0,NX  )=X(2,NX  )+DYDKS*FH
Y(0,NX  )=Y(2,NX  )-DXDKS*FH
X(1,NXP1)=X(1,NXM1)+(Y(2,NX)-Y(0,NX))/FH
Y(1,NXP1)=Y(1,NXM1)-(X(2,NX)-X(0,NX))/FH 
!-----------------------------------------------------------------------------------------------!
DXDKS     =X(NY,NX)-X(NY,NXM1)
DYDKS     =Y(NY,NX)-Y(NY,NXM1)
DXDET     =X(NY,NX)-X(NYM1,NX)
DYDET     =Y(NY,NX)-Y(NYM1,NX)
DSDET     =DXDET**2+DYDET**2
DSDKS     =DXDKS**2+DYDKS**2
DSDKS     =DMAX1(DSDKS,ZERO)
FH        =2.D0*DSQRT(DSDET/DSDKS)
X(NYP1,NX)=X(NYM1,NX)-DYDKS*FH
Y(NYP1,NX)=Y(NYM1,NX)+DXDKS*FH
X(NY,NXP1)=X(NY,NXM1)+(Y(NYP1,NX)-Y(NYM1,NX))/FH
Y(NY,NXP1)=Y(NY,NXM1)-(X(NYP1,NX)-X(NYM1,NX))/FH
!-----------------------------------------------------------------------------------------------!
DXDKS    =X(NY,2)-X(NY  ,1)
DYDKS    =Y(NY,2)-Y(NY  ,1)
DXDET    =X(NY,1)-X(NYM1,1)
DYDET    =Y(NY,1)-Y(NYM1,1)
DSDET    =DXDET**2+DYDET**2
DSDKS    =DXDKS**2+DYDKS**2
DSDKS    =DMAX1(DSDKS,ZERO)
FH       =2.D0*DSQRT(DSDET/DSDKS)
X(NYP1,1)=X(NYM1,1)-DYDKS*FH
Y(NYP1,1)=Y(NYM1,1)+DXDKS*FH
X(NY  ,0)=X(NY  ,2)-(Y(NYP1,1)-Y(NYM1,1))/FH
Y(NY  ,0)=Y(NY  ,2)+(X(NYP1,1)-X(NYM1,1))/FH
!-----------------------------------------------------------------------------------------------!
!    Boundary Eta=0.0 (J=1)                                                                     !
!    Boundary Eta=1.0 (J=NY)                                                                    !
!-----------------------------------------------------------------------------------------------!
DO I=2,NXM1
   IM1=I-1
   IP1=I+1
!-----------------------------------------------------------------------------------------------!
   DXDKS =X(1,IP1)-X(1,IM1)
   DYDKS =Y(1,IP1)-Y(1,IM1)
   DXDET =X(2,I)-X(1,I)
   DYDET =Y(2,I)-Y(1,I)
   DSDET =DXDET**2+DYDET**2
   DSDKS =DXDKS**2+DYDKS**2
   DSDKS =DMAX1(DSDKS,ZERO)
   FH    =2.D0*DSQRT(DSDET/DSDKS)
   X(0,I)=X(2,I)+DYDKS*FH
   Y(0,I)=Y(2,I)-DXDKS*FH
!-----------------------------------------------------------------------------------------------!
   DXDKS =X(NY,IP1)-X(NY,IM1)
   DYDKS =Y(NY,IP1)-Y(NY,IM1)
   DXDET =X(NY,I)-X(NYM1,I)
   DYDET =Y(NY,I)-Y(NYM1,I)
   DSDET =DXDET**2+DYDET**2
   DSDKS =DXDKS**2+DYDKS**2
   DSDKS =DMAX1(DSDKS,ZERO)
   FH    =2.D0*DSQRT(DSDET/DSDKS)
   X(NYP1,I)=X(NYM1,I)-DYDKS*FH
   Y(NYP1,I)=Y(NYM1,I)+DXDKS*FH
END DO !I=2,NXM1
!-----------------------------------------------------------------------------------------------!
!    Boundary XI     (I=1)                                                                      !
!    Boundary XI=1.0 (I=NX)                                                                     !
!-----------------------------------------------------------------------------------------------!
DO J=2,NYM1
   JM1=J-1
   JP1=J+1
!-----------------------------------------------------------------------------------------------!
   DXDKS =X(J,2)-X(J,1)
   DYDKS =Y(J,2)-Y(J,1)
   DXDET =X(JP1,1)-X(JM1,1)
   DYDET =Y(JP1,1)-Y(JM1,1)
   DSDET =DXDET**2+DYDET**2
   DSDKS =DXDKS**2+DYDKS**2
   DSDET =DMAX1(DSDET,ZERO)
   FH    =2.D0*DSQRT(DSDKS/DSDET)
   X(J,0)=X(J,2)-DYDET*FH
   Y(J,0)=Y(J,2)+DXDET*FH
!-----------------------------------------------------------------------------------------------!
   DXDKS    =X(J,NX)-X(J,NXM1)
   DYDKS    =Y(J,NX)-Y(J,NXM1)
   DXDET    =X(JP1,NX)-X(JM1,NX)
   DYDET    =Y(JP1,NX)-Y(JM1,NX)
   DSDET    =DXDET**2+DYDET**2
   DSDKS    =DXDKS**2+DYDKS**2
   DSDET    =DMAX1(DSDET,ZERO)
   FH       =2.D0*DSQRT(DSDKS/DSDET)
   X(J,NXP1)=X(J,NXM1)+DYDET*FH
   Y(J,NXP1)=Y(J,NXM1)-DXDET*FH
END DO !J=2,NYM1
!-----------------------------------------------------------------------------------------------!
X(0,0)=0.5D0*(X(0,1)+X(1,0))
Y(0,0)=0.5D0*(Y(0,1)+Y(1,0))
X(0,NXP1)=0.5D0*(X(0,NX)+X(1,NXP1))
Y(0,NXP1)=0.5D0*(Y(0,NX)+Y(1,NXP1))
X(NYP1,NXP1)=0.5D0*(X(NY,NXP1)+X(NYP1,NX))
Y(NYP1,NXP1)=0.5D0*(Y(NY,NXP1)+Y(NYP1,NX))
X(NYP1,0)   =0.5D0*(X(NYP1,1)+X(NY,0))
Y(NYP1,0)   =0.5D0*(Y(NYP1,1)+Y(NY,0))
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE BORD
!-----------------------------------------------------------------------------------------------!