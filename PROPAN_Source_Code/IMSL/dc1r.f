C
      SUBROUTINE dC1R (N, R, LDR, NER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    N, LDR, NER
      double precision R(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    I, J
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1STI, e1std
C
      DO 20  I=1, N
         IF (R(I+LDR*(I-1)) .EQ. 0.0D0) THEN
            DO 10  J=I + 1, N
               IF (R(I+LDR*(J-1)) .NE. 0.0D0) THEN
                  CALL E1STI (1, I)
                  CALL E1STI (2, J)
                  CALL e1std (1, R(I+LDR*(I-1)))
                  CALL e1std (2, R(I+LDR*(J-1)))
                  CALL E1MES (5, NER, 'R(%(I1),%(I1)) = %(d1). '//
     &                        ' Remaining elements for the row must '//
     &                        'also be zero, but R(%(I1),%(I2)) = '//
     &                        '%(d2) is not.')
                  GO TO 9000
               END IF
   10       CONTINUE
         END IF
   20 CONTINUE
      NER = NER + 1
 9000 RETURN
      END
