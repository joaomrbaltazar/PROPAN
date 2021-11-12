C
      CHARACTER *1 FUNCTION E3PRT (I)
      integer    I
C                                  SPECIFICATIONS FOR INTRINSICS
c     INTRINSIC  CHAR
      INTRINSIC  CHAR
      CHARACTER  CHAR
C
      IF (I .LT. 0) THEN
         E3PRT = CHAR(0)
      ELSE
         E3PRT = CHAR(I)
      END IF
C
      RETURN
      END
