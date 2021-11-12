C
      SUBROUTINE M1VECH (STR1, LEN1, STR2, LEN2)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    LEN1, LEN2
      CHARACTER  STR1*(*), STR2*(*)
C
      STR2(1:LEN2) = STR1(1:LEN1)
C
      RETURN
      END
