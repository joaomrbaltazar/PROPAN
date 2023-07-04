!-----------------------------------------------------------------------------------------------!
      subroutine BIsof (A,LDA,N,B)
!-----------------------------------------------------------------------------------------------!
!    VAN DER VORST,1992 (MORTADELA)                                                             !
!    Activar instrucoes assinaladas com a), b) ou c)                                            !
!    Colocar C (coment rio) para desactivar a), b) ou c)                                        !
!    CONTADOR DE ITERACOES : ITER                                                               !
!    CONTADOR DE MULTIPLICACOES MATR.-VECT. : MVM                                               !
!c    parameter (M=2168) Leonel                                                                 !
!ccc  parameter (M=408)   Leonel                                                                !
!-----------------------------------------------------------------------------------------------!
      IMPLICIT NONE
      INTEGER LDA,N,I,MVM,ITER
      dimension a(LDA,N),b(N)
c      DIMENSION r(M),R0(M),v(M),z(M),d(M),x(M),y(M) !Leonel
      DIMENSION r(LDA),R0(LDA),v(LDA),z(LDA),d(LDA),x(LDA),y(LDA)
C      dimension dij(M)
      DOUBLE PRECISION A,B
      double precision r,v,z,d,x,y,R0
      double precision ALFA,ROI,R0V,ZZ,RZ,W,RON,BETA,RR,EPS,TOL,RES
c      DATA EPS/5.0E-5/    ! Leonel
      DATA EPS/5.D-10/
c
C        escalamento diagonal para matrizes n„o sim‚tricas
CCC   do 160 i=1, n
CCC   yi = A(I,I)
CCC   do 150 j=1, n
C        f) escalamento preliminar pela diagonal
CCC   a(i,j) = a(i,j)/yi
150      continue
CCC   b(i) = b(i)/yi
160   continue
C
c        get rhs ready
      call forwd (a,LDA,n,b)  ! aten‡„o : s¢ para b) ou c)
      do 10 i=1, N
      x(i) = b(i)
10    continue
c1      CALL MVPRD (A,X,V,LDA,N)    ! a) prec. diagonal
      call mulsgs (a,x,v,y,LDA,n)   ! c) Prec. Gauss-Seidel Sim‚trico      
      do 20 i=1, N
      r0(i) = b(i) - v(i)
      R(I) = R0(I)
      d(i) = r(i)
20    continue
C
      call dvdot (ROI,r0,1,R0,1,n)
      MVM = 2
      TOL = EPS*ROI*EPS
      do 30 iter=1, n
c1      CALL MVPRD (A,D,V,LDA,N)    ! a) prec. diagonal
      call mulsgs (a,d,v,y,LDA,n)   ! c) Prec. Gauss-Seidel Sim‚trico        
      MVM = MVM + 1
      call dvdot (r0v,r0,1,v,1,n)
      ALFA = ROI/R0V
      call dvpiv (-ALFA,v,1,r,1,R,1,n)      
c1      CALL MVPRD (A,R,Z,LDA,N)    ! a) prec. diagonal
      call mulsgs (a,R,Z,y,LDA,n)   ! c) Prec. Gauss-Seidel Sim‚trico      
      MVM = MVM + 1      
      call dvdot (ZZ,Z,1,Z,1,n)      
      call dvdot (RZ,R,1,Z,1,n)       
      w = rZ/ZZ
      call dvpiv (ALFA,D,1,x,1,x,1,n)      
      call dvpiv (w,R,1,x,1,x,1,n)
      call dvpiv (-w,Z,1,R,1,r,1,n)
      call dvdot (RON,r0,1,R,1,n) 
      BETA = (RON/ROI)*ALFA/W
      ROI = RON
      call dvpiv (-W,V,1,D,1,D,1,n)
      call dvpiv (BETA,D,1,R,1,D,1,n) 
      CALL DVDOT (RR,R,1,R,1,N)
      RES = 0.5D0*DLOG10(RR)
C        ESCREVA AGORA (ITER,RES) OU (MVM,RES)
!!    WRITE (10,*) ITER,MVM,RES  
CCC     WRITE (12,*) MVM,RES
      IF (RR .LT. TOL) THEN
c1c2      write (*,*) x    ! a) Prec. diagonal ou b) Prec. Gauss-Seidel
      CALL BACKS (A,X,Y,LDA,N)     ! c) Prec. Gauss-Seidel Sim‚trico
CCC     write (*,*) Y            ! c) Prec. Gauss-Seidel Sim‚trico
      DO 80 I=1, N
      B(I) = Y(I)
  80  ENDDO
      RETURN
      ENDIF
30    continue
      stop 001
!-----------------------------------------------------------------------------------------------!
      END SUBROUTINE BISOF
!===============================================================================================!
      SUBROUTINE FORWD (A,LDA,N,B) 
!-----------------------------------------------------------------------------------------------!
!    SOLVE (D+L)*Y=B AND STORE THE RESULT IN B                                                  !
!-----------------------------------------------------------------------------------------------!
      IMPLICIT NONE
      INTEGER LDA,K,N
      DOUBLE PRECISION A(LDA,1),B(1),T
!-----------------------------------------------------------------------------------------------!
!    SOLVE (D+L)*Y=B                                                                            !
!-----------------------------------------------------------------------------------------------!
       DO 10 K=1,N-1 !Joao Baltazar, IST, August 2007 
         T = B(K)/A(K,K)           
         B(K) = T
         CALL DVPIV (-T,A(K+1,K),1,B(K+1),1,B(K+1),1,N-K) 
10       CONTINUE 
      RETURN
!-----------------------------------------------------------------------------------------------!
      END SUBROUTINE FORWD
!===============================================================================================!
      SUBROUTINE MULSGS(A,X,Y,Z,LDA,N)
!-----------------------------------------------------------------------------------------------!
!    SYMMETRIC GAUSS-SEIDEL PRECONDITIONED MATRIX-VECTOR PRODUCT                                !
!    CAN BE USED FOR BOTH SYMMETRIC AND UNSYMMETRIC SOLVERS                                     !
!    Z IS AN AUXILIARY VECTOR                                                                   !
!-----------------------------------------------------------------------------------------------!
      IMPLICIT NONE
      INTEGER :: LDA,N,I
      DOUBLE PRECISION :: A(LDA,N),X(N),Y(N),Z(N)
!-----------------------------------------------------------------------------------------------!
!    UPPER TRIANGULAR MATRIX-VECTOR PRODUCT                                                     !
!-----------------------------------------------------------------------------------------------!
      CALL BACKS(A,X,Z,LDA,N)
      DO I=1,N
         Y(I)=X(I)-A(I,I)*Z(I)
      END DO !I=1,N
!-----------------------------------------------------------------------------------------------!
!    RESULT IN Y                                                                                !
!-----------------------------------------------------------------------------------------------!
      CALL FORWD (A,LDA,N,Y)
      DO I=1,N
         Y(I)=Y(I)+Z(I)
      END DO !I=1,N
!-----------------------------------------------------------------------------------------------!
      END SUBROUTINE MULSGS
!===============================================================================================!
      SUBROUTINE DVDOT(S,V1,INCR1,V2,INCR2,N)
!-----------------------------------------------------------------------------------------------!
      IMPLICIT NONE
      INTEGER          :: J1,J2,N,I,INCR1,INCR2
      DOUBLE PRECISION :: S, V1(1), V2(1)
      DOUBLE PRECISION :: SS,VV1,VV2
!-----------------------------------------------------------------------------------------------!
      SS=0.D0      
      IF (N > 0) THEN
         J1=1
         J2=1
         DO I=1,N
            VV1=V1(J1)
            VV2=V2(J2)
            SS=SS+VV1*VV2
            J1=J1+INCR1
            J2=J2+INCR2
         END DO !I=1,N
      END IF !(N > 0)
      S=SS
!-----------------------------------------------------------------------------------------------!
      END SUBROUTINE DVDOT
!===============================================================================================!
      SUBROUTINE DVPIV(S,V1,INCR1,V2,INCR2,V3,INCR3,N)
!-----------------------------------------------------------------------------------------------!
      IMPLICIT NONE
      INTEGER          :: J1,J2,J3,I,N,INCR1,INCR2,INCR3
      DOUBLE PRECISION :: S,V1(1),V2(1),V3(1)
      DOUBLE PRECISION :: SS,VV1,VV2
!-----------------------------------------------------------------------------------------------!
      IF (N .LT. 1) RETURN
      J1 = 1
      J2 = 1
      J3 = 1
      SS = S      
      DO 10 I=1, N
      VV1 = V1(J1)
      VV2 = V2(J2)
      V3(J3) = SS*VV1 + VV2
      J1 = J1 + INCR1
      J2 = J2 + INCR2
   10 J3 = J3 + INCR3
      RETURN
!-----------------------------------------------------------------------------------------------!
      END SUBROUTINE DVPIV
!===============================================================================================!
      SUBROUTINE BACKS(A,B,R,LDA,N) 
!-----------------------------------------------------------------------------------------------!
!    SOLVE (D+U)*R=B AND STORE THE RESULT IN R                                                  !
!-----------------------------------------------------------------------------------------------!
      IMPLICIT NONE
      INTEGER          :: LDA,I,N,KB,K
      DOUBLE PRECISION :: A(LDA,1),B(1),R(1),T
!-----------------------------------------------------------------------------------------------!
      DO I=1,N
         R(I)=B(I)
      END DO !I=1,N
!-----------------------------------------------------------------------------------------------!
!    NOW SOLVE (D+U)*R=B                                                                        !
!-----------------------------------------------------------------------------------------------!
      DO KB=1,N
         K=N+1-KB 
         R(K)=R(K)/A(K,K) 
         T=-1.D0*R(K)
         CALL DVPIV(T,A(1,K),1,R(1),1,R(1),1,K-1) 
      END DO !KB=1,N
!-----------------------------------------------------------------------------------------------!
      END SUBROUTINE BACKS
!-----------------------------------------------------------------------------------------------!