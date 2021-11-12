c
      SUBROUTINE dGIRTS (N, R, LDR, NB, B, LDB, IPATH, IRANK, X, LDX,
     &                   RINV, LDRINV)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    N, LDR, NB, LDB, IPATH, IRANK, LDX, LDRINV
      double precision R(*), B(*), X(*), RINV(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    I, J, K, NER
      double precision TEMP1, TEMP2
C                                  SPECIFICATIONS FOR INTRINSICS
c     INTRINSIC  dabs
      INTRINSIC  dabs
      double precision dabs
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   C1DIM, C1IARG, E1MES, E1POP, E1PSH, E1STI, e1std,
     &           daxpy, dcopy, dger, dscal, dset, dtrsv, dC1R, dC1TRG
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   dmach, N1RTY, ddot, dA1OT
      integer    N1RTY
      double precision dmach, ddot, dA1OT
C
      CALL E1PSH ('dGIRTS ')
C                                  Check for terminal errors
      NER = 1
      CALL C1DIM (1, N, 'N', LDR, 'LDR', NER)
      CALL C1IARG (NB, 'NB', 0, -1, NER)
      CALL C1IARG (IPATH, 'IPATH', 1, 4, NER)
      IF (NB .GT. 0) THEN
         CALL C1DIM (1, N, '*N', LDB, 'LDB', NER)
         CALL C1DIM (1, N, '*N', LDX, 'LDX', NER)
      ELSE
         NER = NER + 6
      END IF
      IF (IPATH.EQ.3 .OR. IPATH.EQ.4) THEN
         CALL C1DIM (1, N, '*N', LDRINV, 'LDRINV', NER)
      END IF
      IF (N1RTY(0) .NE. 0) GO TO 9000
C                                  Linear dependent rows of R must be
C                                    represented only by rows whose
C                                    elements are zero.
      CALL dC1R (N, R, LDR, NER)
      IF (N1RTY(0) .NE. 0) GO TO 9000
C                                  Get rank
      IRANK = 0
      DO 10  I=1, N
         IF (R(I+LDR*(I-1)) .NE. 0.0D0) IRANK = IRANK + 1
   10 CONTINUE
C                                  Make A copy of B in X and work with
C                                    X
      DO 20  J=1, NB
         CALL dcopy (N, B(1+LDB*(J-1)), 1, X(1+LDX*(J-1)), 1)
   20 CONTINUE
C
      IF (IPATH.EQ.1 .OR. IPATH.EQ.3) THEN
C                                  Solve R*X = B
         IF (IRANK .LT. N) THEN
            DO 40  I=1, NB
               DO 30  J=N, 1, -1
                  IF (R(J+LDR*(J-1)) .EQ. 0.0D0) THEN
                     IF (X(J+LDX*(I-1)) .NE. 0.0D0) THEN
                        CALL E1STI (1, J)
                        CALL E1STI (2, I)
                        CALL e1std (1, X(J+LDX*(I-1)))
                        CALL E1MES (3, 1, 'The linear system of '//
     &                              'equations is inconsistent '//
     &                              'within the computed tolerance.  '//
     &                              'Elements of row %(I1) are zero, '//
     &                              'but B(%(I1),%(I2)) = %(d1).  '//
     &                              'X(%(I1),%(I2)) is set to '//
     &                              'zero.')
                     END IF
                     X(J+LDX*(I-1)) = 0.0D0
                  ELSE
                     X(J+LDX*(I-1)) = X(J+LDX*(I-1))/R(J+LDR*(J-1))
                     TEMP1 = -X(J+LDX*(I-1))
                     CALL daxpy (J-1, TEMP1, R(1+LDR*(J-1)), 1,
     &                           X(1+LDX*(I-1)), 1)
                  END IF
   30          CONTINUE
   40       CONTINUE
         ELSE
            DO 50  I=1, NB
               CALL dtrsv ('UPPER', 'NOT-TRANS', 'NOT-DIAG', N, R,
     &                     LDR, X(1+LDX*(I-1)), 1)
   50       CONTINUE
         END IF
C
      ELSE IF (IPATH.EQ.2 .OR. IPATH.EQ.4) THEN
C                                  Solve R'*X=B
         IF (IRANK .LT. N) THEN
C                                  Case of singular R
            DO 70  I=1, NB
               DO 60  J=1, N
                  TEMP1 = X(J+LDX*(I-1)) - ddot(J-1,R(1+LDR*(J-1)),1,
     &                    X(1+LDX*(I-1)),1)
                  IF (R(J+LDR*(J-1)) .EQ. 0.0D0) THEN
                     TEMP2 = dabs(X(J+LDX*(I-1))) +
     &                       dA1OT(J-1,R(1+LDR*(J-1)),1,X(1+LDX*(I-1)),
     &                       1)
                     TEMP2 = TEMP2*200.0D0*dmach(4)
                     IF (dabs(TEMP1) .GT. TEMP2) THEN
                        CALL E1STI (1, J)
                        CALL E1STI (2, I)
                        CALL E1MES (3, 1, 'The linear system of '//
     &                              'equations is inconsistent '//
     &                              'within the computed tolerance.  '//
     &                              'X(%(I1),%(I2)) is set to zero.')
                     END IF
                     X(J+LDX*(I-1)) = 0.0D0
                  ELSE
                     X(J+LDX*(I-1)) = TEMP1/R(J+LDR*(J-1))
                  END IF
   60          CONTINUE
   70       CONTINUE
         ELSE
C                                  Case of nonsingular R
            DO 80  I=1, NB
               CALL dtrsv ('UPPER', 'TRANSPOSE', 'NOT-UNIT', N, R,
     &                     LDR, X(1+LDX*(I-1)), 1)
   80       CONTINUE
         END IF
      END IF
      IF (IPATH.EQ.3 .OR. IPATH.EQ.4) THEN
C                                  Invert R
         DO 90  J=1, N
            CALL dcopy (J, R(1+LDR*(J-1)), 1, RINV(1+LDRINV*(J-1)), 1)
   90    CONTINUE
         DO 100  K=1, N
            IF (RINV(K+LDRINV*(K-1)) .EQ. 0.0D0) THEN
               CALL dset (K, 0.0D0, RINV(1+LDRINV*(K-1)), 1)
               CALL dset (N-K, 0.0D0, RINV(K+LDRINV*K), LDRINV)
            ELSE
               RINV(K+LDRINV*(K-1)) = 1.0D0/RINV(K+LDRINV*(K-1))
               TEMP1 = -RINV(K+LDRINV*(K-1))
               CALL dscal (K-1, TEMP1, RINV(1+LDRINV*(K-1)), 1)
               IF (K .LT. N) THEN
                  CALL dger (K-1, N-K, 1.0D0, RINV(1+LDRINV*(K-1)), 1,
     &                       RINV(K+LDRINV*K), LDRINV,
     &                       RINV(1+LDRINV*K), LDRINV)
                  CALL dscal (N-K, RINV(K+LDRINV*(K-1)),
     &                        RINV(K+LDRINV*K), LDRINV)
               END IF
            END IF
  100    CONTINUE
C                                  Fill lower triangle of RINV with
C                                    zeros
         CALL dC1TRG (N, RINV, LDRINV)
      END IF
C                                  Exit section
 9000 CALL E1POP ('dGIRTS ')
      RETURN
      END
