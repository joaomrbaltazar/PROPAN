C
      integer FUNCTION I1ERIF (STR1, LEN1, STR2, LEN2)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    LEN1, LEN2
      CHARACTER  STR1(*), STR2(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    I
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   I1X
      integer    I1X
C
      IF (LEN1.LE.0 .OR. LEN2.LE.0) THEN
         I1ERIF = 1
      ELSE
         DO 10  I=1, LEN1
            IF (I1X(STR2,LEN2,STR1(I),1) .EQ. 0) THEN
               I1ERIF = I
               GO TO 9000
            END IF
   10    CONTINUE
         I1ERIF = 0
      END IF
C
 9000 RETURN
      END
