C
      SUBROUTINE C1IARG (IARG, NMARG, INT1, INT2, NER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    IARG, INT1, INT2, NER
      CHARACTER  NMARG*(*)
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1STI, E1STL
C
      IF (INT1 .LE. INT2) THEN
         IF (IARG.LT.INT1 .OR. IARG.GT.INT2) THEN
            CALL E1STI (1, IARG)
            CALL E1STI (2, INT1)
            CALL E1STI (3, INT2)
            CALL E1STL (1, NMARG)
            CALL E1MES (5, NER, '%(L1) = %(I1).  %(L1) must be '//
     &                  'greater than or equal to %(I2) and less '//
     &                  'than or equal to %(I3).')
         END IF
      ELSE
         IF (IARG .LT. INT1) THEN
            CALL E1STI (1, IARG)
            CALL E1STI (2, INT1)
            CALL E1STL (1, NMARG)
            CALL E1MES (5, NER, '%(L1) = %(I1).  %(L1) must be '//
     &                  'greater than or equal to %(I2).')
         END IF
      END IF
      NER = NER + 1
      RETURN
      END
