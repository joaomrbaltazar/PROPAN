C
      double precision FUNCTION dA1OT (N, SX, INCX, SY, INCY)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    N, INCX, INCY
      double precision SX(*), SY(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    I, IX, IY, M, MP1
C                                  SPECIFICATIONS FOR SPECIAL CASES
c     INTRINSIC  MOD
      INTRINSIC  MOD
      integer    MOD
C                                  SPECIFICATIONS FOR INTRINSICS
c     INTRINSIC  dabs
      INTRINSIC  dabs
      double precision dabs
C
      dA1OT = 0.0D0
      IF (N .GT. 0) THEN
         IF (INCX.NE.1 .OR. INCY.NE.1) THEN
C                                  CODE FOR UNEQUAL INCREMENTS
            IX = 1
            IY = 1
            IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
            IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
            DO 10  I=1, N
               dA1OT = dA1OT + dabs(SX(IX)*SY(IY))
               IX = IX + INCX
               IY = IY + INCY
   10       CONTINUE
         ELSE
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
            M = MOD(N,5)
C                                  CLEAN-UP LOOP SO REMAINING VECTOR
            DO 30  I=1, M
               dA1OT = dA1OT + dabs(SX(I)*SY(I))
   30       CONTINUE
            MP1 = M + 1
            DO 40  I=MP1, N, 5
               dA1OT = dA1OT + dabs(SX(I)*SY(I)) +
     &                 dabs(SX(I+1)*SY(I+1)) + dabs(SX(I+2)*SY(I+2)) +
     &                 dabs(SX(I+3)*SY(I+3)) + dabs(SX(I+4)*SY(I+4))
   40       CONTINUE
         END IF
      END IF
      RETURN
      END
