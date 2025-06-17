!-----------------------------------------------------------------------------------------------!
SUBROUTINE RHS(RXF,ALFA,LBOUN,NX,NY,ITER)
!-----------------------------------------------------------------------------------------------!
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!-----------------------------------------------------------------------------------------------!
PARAMETER(NDIM=300)
!-----------------------------------------------------------------------------------------------!
COMMON /COOR/ X(0:NDIM+1,0:NDIM+1),Y(0:NDIM+1,0:NDIM+1)
COMMON /FORCE/ RXI(NDIM,NDIM), RETA(NDIM,NDIM)
COMMON /TRANS/ FXI(NDIM,NDIM),FETA(NDIM,NDIM)
COMMON /AREA/  DSI(NDIM,4),DST(4),IDSIT(4),IAK(2),JAK(2)
!-----------------------------------------------------------------------------------------------!
DIMENSION LBOUN(4)
DIMENSION G11(NDIM,NDIM),G22(NDIM,NDIM)
DIMENSION RXIL(NDIM,NDIM),RETAL(NDIM,NDIM)
DIMENSION DXDKS(NDIM,4),DYDKS(NDIM,4),DXDET(NDIM,4),DYDET(NDIM,4)
!-----------------------------------------------------------------------------------------------!
SAVE DXDKS,DYDKS,DXDET,DYDET,G11,G22,RXIL,RETAL
!-----------------------------------------------------------------------------------------------!
ZERO=0.D0 
RXF=0.D0
!-----------------------------------------------------------------------------------------------!
NXM1=NX-1
NYM1=NY-1
IF (ITER > 0) GOTO 110
!-----------------------------------------------------------------------------------------------!
DO I=2,NXM1
   IP1=I+1
   IM1=I-1
!  ETA=0 (J=1)
   DXDKS(I,1)=0.5D0*(X(1,IP1)-X(1,IM1))
   DYDKS(I,1)=0.5D0*(Y(1,IP1)-Y(1,IM1))
   G22(1,I)  =DSQRT(DXDKS(I,1)**2+DYDKS(I,1)**2)
   FH        =DSI(I,1)/G22(1,I)
   DYDET(I,1)= FH*DXDKS(I,1)
   DXDET(I,1)=-FH*DYDKS(I,1)
!-----------------------------------------------------------------------------------------------!
   D2XKS=X(1,IP1)+X(1,IM1)-2.D0*X(1,I)
   D2YKS=Y(1,IP1)+Y(1,IM1)-2.D0*Y(1,I)
!-----------------------------------------------------------------------------------------------!
   RXIL(1,I) =-(DXDKS(I,1)*D2XKS+DYDKS(I,1)*D2YKS)
   RETAL(1,I)=-(DXDET(I,1)*D2XKS+DYDET(I,1)*D2YKS)
END DO !I=2,NXM1
!-----------------------------------------------------------------------------------------------!
DO I=2,NXM1
   IP1=I+1
   IM1=I-1
!  ETA=1 (J=NY)
   DXDKS(I,3)=0.5D0*(X(NY,IP1)-X(NY,IM1))
   DYDKS(I,3)=0.5D0*(Y(NY,IP1)-Y(NY,IM1))
   G22(NY,I) =DSQRT(DXDKS(I,3)**2+DYDKS(I,3)**2)
   FH        =DSI(I,3)/G22(NY,I)
   DYDET(I,3)= FH*DXDKS(I,3)
   DXDET(I,3)=-FH*DYDKS(I,3)
!-----------------------------------------------------------------------------------------------!
   D2XKS=X(NY,IP1)+X(NY,IM1)-2.D0*X(NY,I)
   D2YKS=Y(NY,IP1)+Y(NY,IM1)-2.D0*Y(NY,I)
!-----------------------------------------------------------------------------------------------!
   RXIL(NY,I) =-(DXDKS(I,3)*D2XKS+DYDKS(I,3)*D2YKS)
   RETAL(NY,I)=-(DXDET(I,3)*D2XKS+DYDET(I,3)*D2YKS)
END DO !I=2,NXM1
!-----------------------------------------------------------------------------------------------!
!    Xi Lines                                                                                   !
!-----------------------------------------------------------------------------------------------!
DO J=2,NYM1
   JP1=J+1
   JM1=J-1
!  XI=0 (I=1)
   DXDET(J,2)=0.5D0*(X(JP1,1)-X(JM1,1))
   DYDET(J,2)=0.5D0*(Y(JP1,1)-Y(JM1,1))
   G11(J,1)  =DSQRT(DXDET(J,2)**2+DYDET(J,2)**2)
   FH        =DSI(J,2)/G11(J,1)
   DYDKS(J,2)=-FH*DXDET(J,2)
   DXDKS(J,2)= FH*DYDET(J,2)
!-----------------------------------------------------------------------------------------------!
   D2XET=X(JP1,1)+X(JM1,1)-2.D0*X(J,1)
   D2YET=Y(JP1,1)+Y(JM1,1)-2.D0*Y(J,1)
!-----------------------------------------------------------------------------------------------!
   RXIL(J,1) =-(DXDKS(J,2)*D2XET+DYDKS(J,2)*D2YET)
   RETAL(J,1)=-(DXDET(J,2)*D2XET+DYDET(J,2)*D2YET)
END DO !J=2,NYM1
!  XI=1 (I=NX)
DO J=2,NYM1
   JP1=J+1
   JM1=J-1
   DXDET(J,4)=0.5D0*(X(JP1,NX)-X(JM1,NX))
   DYDET(J,4)=0.5D0*(Y(JP1,NX)-Y(JM1,NX))
   G11(J,NX) =DSQRT(DXDET(J,4)**2+DYDET(J,4)**2)
   FH        =DSI(J,4)/G11(J,NX)
   DYDKS(J,4)=-FH*DXDET(J,4)
   DXDKS(J,4)= FH*DYDET(J,4)
!-----------------------------------------------------------------------------------------------!
   D2XET=X(JP1,NX)+X(JM1,NX)-2.D0*X(J,NX)
   D2YET=Y(JP1,NX)+Y(JM1,NX)-2.D0*Y(J,NX)
!-----------------------------------------------------------------------------------------------!
   RXIL(J,NX) =-(DXDKS(J,4)*D2XET+DYDKS(J,4)*D2YET)
   RETAL(J,NX)=-(DXDET(J,4)*D2XET+DYDET(J,4)*D2YET)
END DO !J=2,NYM1
!-----------------------------------------------------------------------------------------------!
DO I=2,NXM1
   DO J=2,NYM1
      FAXI =1.D0-FXI(J,I)
      FBXI =     FXI(J,I)
      FAETA=1.D0-FETA(J,I)
      FBETA=     FETA(J,I)
      G11(J,I)=(G11(J,1)*FAXI + G11(J,NX)*FBXI)**2
      G22(J,I)=(G22(1,I)*FAETA+G22(NY,I)*FBETA)**2
      RXIL(J,I) =(  FAXI*RXIL(J,1)+  FBXI*RXIL(J,NX))/G11(J,I) &
                +( FAETA*RXIL(1,I)+ FBETA*RXIL(NY,I))/G22(J,I)
      RETAL(J,I)=( FAXI*RETAL(J,1)+ FBXI*RETAL(J,NX))/G11(J,I) &
                +(FAETA*RETAL(1,I)+FBETA*RETAL(NY,I))/G22(J,I)
      RXI(J,I) =RXIL(J,I)
      RETA(J,I)=RETAL(J,I)
   END DO !J=2,NYM1
END DO !I=2,NXM1
!-----------------------------------------------------------------------------------------------!
RETURN
!-----------------------------------------------------------------------------------------------!
!    Non-Linear Terms                                                                           !
!-----------------------------------------------------------------------------------------------!
110 ILOW=1+IAK(1)
    JLOW=1+JAK(1)
    ITOP=NX-IAK(2)
    JTOP=NY-JAK(2)
    ACC1=1.D0-ALFA
    ACC2=ALFA
!-----------------------------------------------------------------------------------------------!
DO I=2,NXM1
   D2XET=X(2,I)-X(1,I)-DXDET(I,1)
   D2YET=Y(2,I)-Y(1,I)-DYDET(I,1)
   RXI(1,I) =-(DXDKS(I,1)*D2XET+DYDKS(I,1)*D2YET)
   RETA(1,I)=-(DXDET(I,1)*D2XET+DYDET(I,1)*D2YET)
!-----------------------------------------------------------------------------------------------!
   D2XET     =X(NYM1,I)-X(NY,I)+DXDET(I,3)  
   D2YET     =Y(NYM1,I)-Y(NY,I)+DYDET(I,3)
   RXI(NY,I) =-(DXDKS(I,3)*D2XET+DYDKS(I,3)*D2YET)
   RETA(NY,I)=-(DXDET(I,3)*D2XET+DYDET(I,3)*D2YET)
END DO !I=2,NXM1
!-----------------------------------------------------------------------------------------------!
DO J=2,NYM1
   D2XKS=X(J,2)-X(J,1)-DXDKS(J,2)
   D2YKS=Y(J,2)-Y(J,1)-DYDKS(J,2)
   RXI(J,1) =-(DXDKS(J,2)*D2XKS+DYDKS(J,2)*D2YKS)
   RETA(J,1)=-(DXDET(J,2)*D2XKS+DYDET(J,2)*D2YKS)
!-----------------------------------------------------------------------------------------------!
   D2XKS=X(J,NXM1)-X(J,NX)+DXDKS(J,4)
   D2YKS=Y(J,NXM1)-Y(J,NX)+DYDKS(J,4)
   RXI(J,NX) =-(DXDKS(J,4)*D2XKS+DYDKS(J,4)*D2YKS)
   RETA(J,NX)=-(DXDET(J,4)*D2XKS+DYDET(J,4)*D2YKS)
END DO !J=2,NYM1
!-----------------------------------------------------------------------------------------------!
DO I=2,NXM1
   IP1 =I+1
   IM1 =I-1
   FILOW=DFLOAT(I-ILOW)
   FITOP=DFLOAT(I-ITOP)
   FILOW=-MAX(FILOW,ZERO)*0.25D0
   FITOP= MIN(FITOP,ZERO)*0.25D0
   FX1=DEXP(FILOW)
   FX2=DEXP(FITOP)
   IF (IAK(1) == 0) FX1=0.D0
   IF (IAK(2) == 0) FX2=0.D0
   IF (LBOUN(2) /= 0) FX1=0.D0
   IF (LBOUN(4) /= 0) FX2=0.D0
   DO J=2,NYM1
      JP1 =J+1
      JM1 =J-1
      FJLOW=DFLOAT(J-JLOW)
      FJTOP=DFLOAT(J-JTOP)
      FJLOW=-MAX(FJLOW,ZERO)*0.25D0
      FJTOP= MIN(FJTOP,ZERO)*0.25D0
      FY1=DEXP(FJLOW)
      FY2=DEXP(FJTOP)
      IF (JAK(1) == 0) FY1=0.D0
      IF (JAK(2) == 0) FY2=0.D0
      IF (LBOUN(1) /= 0) FY1=0.D0
      IF (LBOUN(3) /= 0) FY2=0.D0
!-----------------------------------------------------------------------------------------------!
      RADX1=FY1* RXI(1,I)/G11(J,I)
      RADY1=FY1*RETA(1,I)/G11(J,I)
!-----------------------------------------------------------------------------------------------!
      RADX3=FY2* RXI(NY,I)/G11(J,I)
      RADY3=FY2*RETA(NY,I)/G11(J,I)
!-----------------------------------------------------------------------------------------------!
      RADX2=FX1* RXI(J,1)/G22(J,I)
      RADY2=FX1*RETA(J,1)/G22(J,I)
!-----------------------------------------------------------------------------------------------!
      RADX4=FX2* RXI(J,NX)/G22(J,I)
      RADY4=FX2*RETA(J,NX)/G22(J,I)
!-----------------------------------------------------------------------------------------------!
      RNEWX=RXIL(J,I) +RADX1+RADX2+RADX3+RADX4
      RNEWY=RETAL(J,I)+RADY1+RADY2+RADY3+RADY4
      RXF=MAX(RXF,DABS(RNEWX-RXI(J,I)),DABS(RNEWY-RETA(J,I)))
      RXI(J,I) =ACC1*RXI(J,I) +ACC2*RNEWX
      RETA(J,I)=ACC1*RETA(J,I)+ACC2*RNEWY
   END DO !J=2,NYM1
END DO !I=2,NXM1
!-----------------------------------------------------------------------------------------------!
RETURN
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE RHS
!-----------------------------------------------------------------------------------------------!
