!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION FUNCTION SXX(Y)
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
DOUBLE PRECISION :: Y,UM,PI,Y2,Y3,C3,X,YM1
!-----------------------------------------------------------------------------------------------!
IF (Y < 0.26938972D0) THEN 
   UM=1.D0
   PI=4.D0*DATAN(UM)
   Y2=Y*Y
   Y3=Y2*Y
   C3=1.D0+PI*PI/6.D0
   X=1.D0-Y+Y2-C3*Y3+6.794732D0*Y2*Y2-13.205501D0*Y2*Y3+11.726095D0*Y3*Y3      
   X=PI*X
ELSE
   YM1=1.D0-Y
   Y2=YM1*YM1
   Y3=Y2*YM1
   X=1.D0+0.15D0*YM1+0.057321429D0*Y2+0.048774238D0*Y3-0.053337753D0*Y2*Y2
   X=DSQRT(6.D0*YM1)*(X+0.075845134D0*Y2*Y3)
END IF
SXX=X
!-----------------------------------------------------------------------------------------------!
END FUNCTION SXX
!-----------------------------------------------------------------------------------------------!
