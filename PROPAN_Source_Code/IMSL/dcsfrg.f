c
      SUBROUTINE dCSFRG (N, A, LDA)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    N, LDA
      double precision A(LDA,*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    I
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1POP, E1PSH, E1STI, dcopy
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   N1RCD
      integer    N1RCD
C
      CALL E1PSH ('dCSFRG ')
C
      IF (N .LE. 0) THEN
         CALL E1STI (1, N)
         CALL E1MES (5, 1, 'N = %(I1).  The order of A, N, '//
     &               'must be greater than 0.')
      END IF
C
      IF (N .GT. LDA) THEN
         CALL E1STI (1, N)
         CALL E1STI (2, LDA)
         CALL E1MES (5, 2, 'N = %(I1) and LDA = %(I2).  The order of '//
     &               'A, N, must be less than or equal to the '//
     &               'leading dimension of A, LDA.')
      END IF
      IF (N1RCD(0) .NE. 0) GO TO 9000
C                                  Copy upper triangular values to lower
C                                  triangular values
      DO 10  I=1, N - 1
         CALL dcopy (N-I, A(I,I+1), LDA, A(I+1,I), 1)
   10 CONTINUE
C                                  Exit section
 9000 CALL E1POP ('dCSFRG ')
      RETURN
      END
