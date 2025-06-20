!-----------------------------------------------------------------------------------------------!
SUBROUTINE CALCB(IL,IT,NX,NY)
!-----------------------------------------------------------------------------------------------!
!    Created by   : L. Eca, IST                                                                 !
!    Last Revision: Joao Baltazar, IST, April 2005                                              !
!-----------------------------------------------------------------------------------------------!
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!-----------------------------------------------------------------------------------------------!
PARAMETER(NDIM=300)
!-----------------------------------------------------------------------------------------------!
COMMON /COOR/  X(0:NDIM+1,0:NDIM+1),Y(0:NDIM+1,0:NDIM+1)
COMMON /BOUN/  XL(NDIM,4),YL(NDIM,4),DXDL(NDIM,4),DYDL(NDIM,4),RL(NDIM,4)
!-----------------------------------------------------------------------------------------------!
DIMENSION PHI(4),DPH(4),D2P(4)
!-----------------------------------------------------------------------------------------------!
NXM1=NX-1
NYM1=NY-1
NXM2=NX-2
NYM2=NY-2
!-----------------------------------------------------------------------------------------------!
IF (IL == 1) GOTO 100
IF (IL == 2) GOTO 200
IF (IL == 3) GOTO 300
IF (IL == 4) GOTO 400
!-----------------------------------------------------------------------------------------------!
100 IF(IT == 1) GOTO 110
    IF(IT == 2) GOTO 130
    IF(IT == 3) GOTO 150 
!-----------------------------------------------------------------------------------------------!
110 DO I=2,NXM1
       Y(1,I)=(4.D0*Y(2,I)-Y(3,I))/3.D0
    END DO !I=2,NXM1
RETURN
!-----------------------------------------------------------------------------------------------!
130 DO I=2,NXM1
       X(1,I)=(4.D0*X(2,I)-X(3,I))/3.D0
    END DO !I=2,NXM1
RETURN
!-----------------------------------------------------------------------------------------------!
150 DO I=2,NXM1
       DXL=2.D0*X(2,I)-0.5D0*X(3,I)
       DYL=2.D0*Y(2,I)-0.5D0*Y(3,I)
       RL0=RL(I,1)
!-----------------------------------------------------------------------------------------------!
       DO NIT=1,10
          KL=INT(RL0)
          RL1=RL0-DFLOAT(KL)
          KL1=KL+1
          CALL CALPHI(PHI,DPH,D2P,RL1)
          XLOC =XL(KL,1)*PHI(1)+XL(KL1,1)*PHI(2)+DXDL(KL,1)*PHI(3)+DXDL(KL1,1)*PHI(4)
          YLOC =YL(KL,1)*PHI(1)+YL(KL1,1)*PHI(2)+DYDL(KL,1)*PHI(3)+DYDL(KL1,1)*PHI(4)
          DXDKS=XL(KL,1)*DPH(1)+XL(KL1,1)*DPH(2)+DXDL(KL,1)*DPH(3)+DXDL(KL1,1)*DPH(4)
          DYDKS=YL(KL,1)*DPH(1)+YL(KL1,1)*DPH(2)+DYDL(KL,1)*DPH(3)+DYDL(KL1,1)*DPH(4)
          D2XKS=XL(KL,1)*D2P(1)+XL(KL1,1)*D2P(2)+DXDL(KL,1)*D2P(3)+DXDL(KL1,1)*D2P(4)
          D2YKS=YL(KL,1)*D2P(1)+YL(KL1,1)*D2P(2)+DYDL(KL,1)*D2P(3)+DYDL(KL1,1)*D2P(4)
!-----------------------------------------------------------------------------------------------!
          DXDET=DXL-1.5D0*XLOC
          DYDET=DYL-1.5D0*YLOC
          FLOC =DXDKS*DXDET+DYDKS*DYDET
          DFLOC=-1.5D0*(DXDKS**2+DYDKS**2)+D2XKS*DXDET+D2YKS*DYDET
          RL0=RL0-FLOC/DFLOC
       END DO !NIT=1,10
       X(1,I)=XLOC   
       Y(1,I)=YLOC
       RL(I,1)=RL0
    END DO !I=2,NXM1
RETURN
!-----------------------------------------------------------------------------------------------!
200 IF (IT == 1) GOTO 210
    IF (IT == 2) GOTO 230
    IF (IT == 3) GOTO 250 
!-----------------------------------------------------------------------------------------------!
210 DO J=2,NYM1
       Y(J,1)=(4.D0*Y(J,2)-Y(J,3))/3.D0
    END DO !J=2,NYM1
RETURN
!-----------------------------------------------------------------------------------------------!
230 DO J=2,NYM1
       X(J,1)=(4.D0*X(J,2)-X(J,3))/3.D0
    END DO !J=2,NYM1
RETURN
!-----------------------------------------------------------------------------------------------!
250 DO J=2,NYM1
       DXL=2.D0*X(J,2)-0.5D0*X(J,3)
       DYL=2.D0*Y(J,2)-0.5D0*Y(J,3)
       RL0=RL(J,2)
!-----------------------------------------------------------------------------------------------!
       DO NIT=1,10
          JL=INT(RL0)
          RL1=RL0-DFLOAT(JL)
          JL1=JL+1
          CALL CALPHI(PHI,DPH,D2P,RL1)
          XLOC =XL(JL,2)*PHI(1)+XL(JL1,2)*PHI(2)+DXDL(JL,2)*PHI(3)+DXDL(JL1,2)*PHI(4)
          YLOC =YL(JL,2)*PHI(1)+YL(JL1,2)*PHI(2)+DYDL(JL,2)*PHI(3)+DYDL(JL1,2)*PHI(4)
          DXDET=XL(JL,2)*DPH(1)+XL(JL1,2)*DPH(2)+DXDL(JL,2)*DPH(3)+DXDL(JL1,2)*DPH(4)
          DYDET=YL(JL,2)*DPH(1)+YL(JL1,2)*DPH(2)+DYDL(JL,2)*DPH(3)+DYDL(JL1,2)*DPH(4)
          D2XET=XL(JL,2)*D2P(1)+XL(JL1,2)*D2P(2)+DXDL(JL,2)*D2P(3)+DXDL(JL1,2)*D2P(4)
          D2YET=YL(JL,2)*D2P(1)+YL(JL1,2)*D2P(2)+DYDL(JL,2)*D2P(3)+DYDL(JL1,2)*D2P(4)
!-----------------------------------------------------------------------------------------------!
          DXDKS=DXL-1.5D0*XLOC
          DYDKS=DYL-1.5D0*YLOC
          FLOC =DXDKS*DXDET+DYDKS*DYDET
          DFLOC=-1.5D0*(DXDET**2+DYDET**2)+DXDKS*D2XET+DYDKS*D2YET
          RL0=RL0-FLOC/DFLOC
       END DO !NIT=1,10
       X(J,1)=XLOC   
       Y(J,1)=YLOC
       RL(J,2)=RL0
    END DO !J=2,NYM1
RETURN
!-----------------------------------------------------------------------------------------------!
300 IF (IT == 1) GOTO 310
    IF (IT == 2) GOTO 330
    IF (IT == 3) GOTO 350 
!-----------------------------------------------------------------------------------------------!
310 DO I=2,NXM1
       Y(NY,I)=(4.D0*Y(NYM1,I)-Y(NYM2,I))/3.D0
    END DO !I=2,NXM1
RETURN
!-----------------------------------------------------------------------------------------------!
330 DO I=2,NXM1
       X(NY,I)=(4.D0*X(NYM1,I)-X(NYM2,I))/3.D0
    END DO !I=2,NXM1
RETURN
!-----------------------------------------------------------------------------------------------!
350 DO I=2,NXM1
       DXL=0.5D0*X(NYM2,I)-2.D0*X(NYM1,I)
       DYL=0.5D0*Y(NYM2,I)-2.D0*Y(NYM1,I)
       RL0=RL(I,3)
!-----------------------------------------------------------------------------------------------!
       DO NIT=1,3
          KL=INT(RL0)
          RL1=RL0-DFLOAT(KL)
          KL1=KL+1
          CALL CALPHI(PHI,DPH,D2P,RL1)
          XLOC =XL(KL,3)*PHI(1)+XL(KL1,3)*PHI(2)+DXDL(KL,3)*PHI(3)+DXDL(KL1,3)*PHI(4)
          YLOC =YL(KL,3)*PHI(1)+YL(KL1,3)*PHI(2)+DYDL(KL,3)*PHI(3)+DYDL(KL1,3)*PHI(4)
          DXDKS=XL(KL,3)*DPH(1)+XL(KL1,3)*DPH(2)+DXDL(KL,3)*DPH(3)+DXDL(KL1,3)*DPH(4)
          DYDKS=YL(KL,3)*DPH(1)+YL(KL1,3)*DPH(2)+DYDL(KL,3)*DPH(3)+DYDL(KL1,3)*DPH(4)
          D2XKS=XL(KL,3)*D2P(1)+XL(KL1,3)*D2P(2)+DXDL(KL,3)*D2P(3)+DXDL(KL1,3)*D2P(4)
          D2YKS=YL(KL,3)*D2P(1)+YL(KL1,3)*D2P(2)+DYDL(KL,3)*D2P(3)+DYDL(KL1,3)*D2P(4)
!-----------------------------------------------------------------------------------------------!
          DXDET=1.5D0*XLOC+DXL
          DYDET=1.5D0*YLOC+DYL
          FLOC =DXDKS*DXDET+DYDKS*DYDET
          DFLOC=1.5D0*(DXDKS**2+DYDKS**2)+D2XKS*DXDET+D2YKS*DYDET
          RL0=RL0-FLOC/DFLOC
       END DO !NIT=1,3
!-----------------------------------------------------------------------------------------------!
       X(NY,I)=XLOC   
       Y(NY,I)=YLOC
       RL(I,3)=RL0
    END DO !I=2,NXM1
RETURN
!-----------------------------------------------------------------------------------------------!
400 IF (IT == 1) GOTO 410
    IF (IT == 2) GOTO 430
    IF (IT == 3) GOTO 450 
!-----------------------------------------------------------------------------------------------!
410 DO J=2,NYM1
       Y(J,NX)=(4.D0*Y(J,NXM1)-Y(J,NXM2))/3.D0
    END DO !J=2,NYM1
RETURN
!-----------------------------------------------------------------------------------------------!
430 DO J=2,NYM1
       X(J,NX)=(4.D0*X(J,NXM1)-X(J,NXM2))/3.D0
    END DO !J=2,NYM1
RETURN
!-----------------------------------------------------------------------------------------------!
450 DO J=2,NYM1
       DXL=0.5D0*X(J,NXM2)-2.D0*X(J,NXM1)
       DYL=0.5D0*Y(J,NXM2)-2.D0*Y(J,NXM1)
       RL0=RL(J,4)
!-----------------------------------------------------------------------------------------------!
       DO NIT=1,10
          JL=INT(RL0)
          RL1=RL0-DFLOAT(JL)
          JL1=JL+1
          CALL CALPHI(PHI,DPH,D2P,RL1)
          XLOC =XL(JL,4)*PHI(1)+XL(JL1,4)*PHI(2)+DXDL(JL,4)*PHI(3)+DXDL(JL1,4)*PHI(4)
          YLOC =YL(JL,4)*PHI(1)+YL(JL1,4)*PHI(2)+DYDL(JL,4)*PHI(3)+DYDL(JL1,4)*PHI(4)
          DXDET=XL(JL,4)*DPH(1)+XL(JL1,4)*DPH(2)+DXDL(JL,4)*DPH(3)+DXDL(JL1,4)*DPH(4)
          DYDET=YL(JL,4)*DPH(1)+YL(JL1,4)*DPH(2)+DYDL(JL,4)*DPH(3)+DYDL(JL1,4)*DPH(4)
          D2XET=XL(JL,4)*D2P(1)+XL(JL1,4)*D2P(2)+DXDL(JL,4)*D2P(3)+DXDL(JL1,4)*D2P(4)
          D2YET=YL(JL,4)*D2P(1)+YL(JL1,4)*D2P(2)+DYDL(JL,4)*D2P(3)+DYDL(JL1,4)*D2P(4)
!-----------------------------------------------------------------------------------------------!
          DXDKS=1.5D0*XLOC+DXL
          DYDKS=1.5D0*YLOC+DYL
          FLOC =DXDKS*DXDET+DYDKS*DYDET
          DFLOC=1.5D0*(DXDET**2+DYDET**2)+D2XET*DXDKS+D2YET*DYDKS
          RL0=RL0-FLOC/DFLOC
       END DO !NIT=1,10
!-----------------------------------------------------------------------------------------------!
       X(J,NX)=XLOC   
       Y(J,NX)=YLOC
       RL(J,4)=RL0
    END DO !J=2,NYM1
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE CALCB
!-----------------------------------------------------------------------------------------------!
