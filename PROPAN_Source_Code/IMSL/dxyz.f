C
      double precision FUNCTION DXYZ (N, DX, INCX, DY, INCY, DZ, INCZ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    N, INCX, INCY, INCZ
      double precision DX(*), DY(*), DZ(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    I, IX, IY, IZ, M, MP1
C                                  SPECIFICATIONS FOR SPECIAL CASES
c     INTRINSIC  MOD
      INTRINSIC  MOD
      integer    MOD
C
      DXYZ = 0.0D0
      IF (N .LE. 0) GO TO 9000
      IF (INCX.NE.1 .OR. INCY.NE.1 .OR. INCZ.NE.1) THEN
C
C                                  CODE FOR UNEQUAL INCREMENTS OR EQUAL
C                                    INCREMENTS NOT EQUAL TO 1
         IX = 1
         IY = 1
         IZ = 1
         IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
         IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
         IF (INCZ .LT. 0) IZ = (-N+1)*INCZ + 1
         DO 10  I=1, N
            DXYZ = DXYZ + DX(IX)*DY(IY)*DZ(IZ)
            IX = IX + INCX
            IY = IY + INCY
            IZ = IZ + INCZ
   10    CONTINUE
      ELSE
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
         M = MOD(N,3)
C                                  CLEAN-UP LOOP
         DO 30  I=1, M
            DXYZ = DXYZ + DX(I)*DY(I)*DZ(I)
   30    CONTINUE
         MP1 = M + 1
         DO 40  I=MP1, N, 3
            DXYZ = DXYZ + DX(I)*DY(I)*DZ(I) + DX(I+1)*DY(I+1)*DZ(I+1)
     &              + DX(I+2)*DY(I+2)*DZ(I+2)
   40    CONTINUE
      END IF
 9000 RETURN
      END
