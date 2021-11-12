C
      SUBROUTINE dFNLSQ (F, INTCEP, NBASIS, NDATA, XDATA, FDATA, IWT,
     &                   WEIGHT, A, SSE)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    INTCEP, NBASIS, NDATA, IWT
      double precision F, SSE, XDATA(*), FDATA(*), WEIGHT(*), A(*)
      EXTERNAL   F
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    IWKX
C                                  SPECIFICATIONS FOR SPECIAL CASES
C                                  SPECIFICATIONS FOR COMMON /WORKSP/
      real       rwksp(5000)
      double precision rdwksp(2500)
      double precision dwksp(2500)
      complex    cwksp(2500)
      complex    *16 czwksp(1250)
      complex    *16 zwksp(1250)
      integer    iwksp(5000)
      logical    lwksp(5000)
      equivalence (dwksp(1), rwksp(1))
      equivalence (cwksp(1), rwksp(1)), (zwksp(1), rwksp(1))
      equivalence (iwksp(1), rwksp(1)), (lwksp(1), rwksp(1))
      equivalence (rdwksp(1), rwksp(1)), (czwksp(1), rwksp(1))
      common     /worksp/ dwksp
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1POP, E1PSH, E1STI, dF2LSQ
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   I1KGT, N1RTY
      integer    I1KGT, N1RTY
C
      CALL E1PSH ('dFNLSQ ')
C                                  Check INTCEP
      IF (INTCEP.LT.0 .OR. INTCEP.GT.1) THEN
         CALL E1STI (1, INTCEP)
         CALL E1MES (5, 1, 'The intercept option must be a zero or '//
     &               'a one while INTCEP = %(I1) is given.')
      END IF
C                                  Check NBASIS
      IF (NBASIS .LT. 1) THEN
         CALL E1STI (1, NBASIS)
         CALL E1MES (5, 2, 'The number of basis functions must be '//
     &               'at least one while NBASIS = %(I1) is given.')
      END IF
C
      IF (N1RTY(0) .NE. 0) GO TO 9000
C                                  Allocate workspace
      IWKX = I1KGT((INTCEP+NBASIS)**2+4*(INTCEP+NBASIS)+IWT+1,4)
      IF (N1RTY(0) .NE. 0) THEN
         CALL E1MES (5, 5, ' ')
         CALL E1STI (1, INTCEP)
         CALL E1STI (2, NBASIS)
         CALL E1MES (5, 5, 'The workspace requirement is '//
     &               'based on INTCEP = %(I1) and  NBASIS = %(I2).')
      ELSE
C
         CALL dF2LSQ (F, INTCEP, NBASIS, NDATA, XDATA, FDATA, IWT,
     &                WEIGHT, A, SSE, RDWKSP(IWKX))
      END IF
C                                  Exit section
 9000 CALL E1POP ('dFNLSQ ')
      RETURN
      END