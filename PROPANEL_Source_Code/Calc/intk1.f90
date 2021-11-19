!-----------------------------------------------------------------------------------------------!
!    Double Quadratic (smooth) Inter or Extrapolation                                           !
!                                                                                               !
!    INPUT     NDATA = Number of X values and number of function values                         !
!                      N must be Greater than 2                                                 !
!              X     = array of x-values                                                        !
!              Z     = array of function values                                                 !
!              N     = number of XX values and number of interpolated values                    !
!              XX    = X-value to be interpolated                                               !
!    OUTPUT    ZZ    = interpolated function value                                              !
!                                                                                               !
!    remark..the X-values need not be equidistant but must be monoyonically increasing          !
!-----------------------------------------------------------------------------------------------!
SUBROUTINE INTK1(NDATA,X,Z,N,XX,ZZ)
!-----------------------------------------------------------------------------------------------!
!    Revised by: Joao Baltazar, IST, January 2006                                               !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
INTEGER          :: I,II,IQ,IR,J,JJ,K,M,NDATA,N
DOUBLE PRECISION :: X(NDATA),Z(NDATA),XX(N),ZZ(N),FF
!-----------------------------------------------------------------------------------------------!
DO II=1,N
    I=1
    JJ=I
    ZZ(II)=0.D0
    DO 50 J=1,NDATA
    FF=XX(II)-X(J)
    IF (FF) 20,10,50
 10 ZZ(II)=Z(J)
    GOTO 100
 20 IF(J.EQ.1.OR.J.EQ.2) GOTO 60
    IF(J.NE.NDATA) JJ=2
    GOTO 56
 50 CONTINUE
    J=NDATA
 56 I=J-2
 60 M=I+2
    IR=M
    IQ=I+1
    FF=0.D0
    DO 80 K=I,M
    FF=FF+Z(K)*(XX(II)-X(IQ))*(XX(II)-X(IR))/((X(K)-X(IQ))*(X(K)-X(IR)))
    IQ=IR
    IR=K
 80 CONTINUE   
    IF(JJ.NE.1) GOTO 90
    ZZ(II)=FF
    GOTO 100
 90 IR=JJ+I
    ZZ(II)=ZZ(II)+FF*(XX(II)-X(IR))/(X(I+1)-X(IR))
    IF(JJ.EQ.0) GOTO 100
    JJ=0
    I=J-1
    GOTO 60
100 CONTINUE
END DO !II=1,N
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE INTK1
!-----------------------------------------------------------------------------------------------!
