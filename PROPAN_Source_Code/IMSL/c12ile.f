C
      SUBROUTINE C12ILE (IARG1, NMARG1, IARG2, NMARG2, NER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    IARG1, IARG2, NER
      CHARACTER  NMARG1*(*), NMARG2*(*)
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1STI, E1STL
C
      IF (IARG1 .GT. IARG2) THEN
         CALL E1STI (1, IARG1)
         CALL E1STI (2, IARG2)
         CALL E1STL (1, NMARG1)
         CALL E1STL (2, NMARG2)
         CALL E1MES (5, NER, '%(L1) = %(I1) and %(L2) = %(I2).  '//
     &               '%(L1) must be less than or equal to %(L2).')
      END IF
      NER = NER + 1
      RETURN
      END
