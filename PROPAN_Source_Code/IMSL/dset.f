C
      SUBROUTINE DSET (N, DA, DX, INCX)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    N, INCX
      double precision DA, DX(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    I, M, MP1, NINCX
C                                  SPECIFICATIONS FOR SPECIAL CASES
c     INTRINSIC  MOD
      INTRINSIC  MOD
      integer    MOD
C
      IF (N .GT. 0) THEN
         IF (INCX .NE. 1) THEN
C                                  CODE FOR INCREMENT NOT EQUAL TO 1
            NINCX = N*INCX
            DO 10  I=1, NINCX, INCX
               DX(I) = DA
   10       CONTINUE
         ELSE
C                                  CODE FOR INCREMENT EQUAL TO 1
            M = MOD(N,8)
C                                  CLEAN-UP LOOP
            DO 30  I=1, M
               DX(I) = DA
   30       CONTINUE
            MP1 = M + 1
            DO 40  I=MP1, N, 8
               DX(I) = DA
               DX(I+1) = DA
               DX(I+2) = DA
               DX(I+3) = DA
               DX(I+4) = DA
               DX(I+5) = DA
               DX(I+6) = DA
               DX(I+7) = DA
   40       CONTINUE
         END IF
      END IF
      RETURN
      END
