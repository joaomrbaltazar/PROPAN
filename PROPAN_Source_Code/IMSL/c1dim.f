C
      SUBROUTINE C1DIM (INT1, IARG1, NMARG1, IARG2, NMARG2, NER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    INT1, IARG1, IARG2, NER
      CHARACTER  NMARG1*(*), NMARG2*(*)
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   C12ILE, C1IARG
C
      IF (NMARG1(1:1) .EQ. '*') THEN
         NER = NER + 1
         CALL C1IARG (IARG2, NMARG2, 1, -1, NER)
         IF (IARG2 .GE. 1) THEN
            CALL C12ILE (IARG1, NMARG1(2:), IARG2, NMARG2, NER)
         ELSE
            NER = NER + 1
         END IF
      ELSE
         CALL C1IARG (IARG1, NMARG1, INT1, -1, NER)
         CALL C1IARG (IARG2, NMARG2, 1, -1, NER)
         IF (IARG2 .GE. 1) THEN
            CALL C12ILE (IARG1, NMARG1, IARG2, NMARG2, NER)
         ELSE
            NER = NER + 1
         END IF
      END IF
      RETURN
      END
