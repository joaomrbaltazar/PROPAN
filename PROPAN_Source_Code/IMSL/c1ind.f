C
      SUBROUTINE C1IND (INT1, IND, NMIND, NVAR, NMNVAR, NER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    INT1, IND, NVAR, NER
      CHARACTER  NMIND*(*), NMNVAR*(*)
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   C12ILE, C1IARG
C
      CALL C1IARG (IND, NMIND, INT1, -1, NER)
      IF (NVAR .GE. 1) THEN
         CALL C12ILE (IND, NMIND, NVAR, NMNVAR, NER)
      ELSE
         NER = NER + 1
      END IF
      RETURN
      END
