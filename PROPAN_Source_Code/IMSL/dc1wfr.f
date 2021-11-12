C
      SUBROUTINE dC1WFR (IDO, ICALL, X, LDX, IOBS, IROW, IFRQ, IWT,
     &                   XMISS, NMISS, FRQ, WT, IGO)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    IDO, ICALL, LDX, IOBS, IROW, IFRQ, IWT, NMISS, IGO
      double precision XMISS, FRQ, WT, X(*)
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1POP, E1PSH, E1STI, e1std
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   dIFNAN
      logical    dIFNAN
C
      CALL E1PSH ('dC1WFR ')
      IGO = 0
      IF (IFRQ .GT. 0) THEN
         FRQ = X((IFRQ-1)*LDX+IOBS)
         IF (dIFNAN(FRQ)) THEN
            NMISS = NMISS + IROW
            IGO = 2
         ELSE IF (FRQ .EQ. 0.0D0) THEN
            IGO = 1
            GO TO 9000
         END IF
      END IF
      IF (IWT .GT. 0) THEN
         WT = X((IWT-1)*LDX+IOBS)
         IF (dIFNAN(WT)) THEN
            IF (IGO .NE. 2) THEN
               NMISS = NMISS + IROW
               IGO = 2
            END IF
         END IF
      END IF
      IF (IFRQ .GT. 0) THEN
         IF (.NOT.dIFNAN(FRQ)) THEN
            IF (FRQ .LT. 0.0D0) THEN
               CALL E1STI (1, IOBS)
               CALL e1std (1, FRQ)
               IF (IDO .GT. 0) THEN
                  CALL E1STI (2, ICALL)
                  CALL E1MES (4, 2, 'The frequency for row '//
     &                        '%(I1) of X on invocation number %(I2) '//
     &                        'of this routine is %(d1).  '//
     &                        'Frequencies must be nonnegative.')
               ELSE
                  CALL E1MES (4, 2, 'The frequency for row '//
     &                        '%(I1) of X is %(d1).  Frequencies '//
     &                        'must be nonnegative.')
               END IF
               IGO = 3
               GO TO 9000
            END IF
         END IF
      ELSE
         FRQ = 1.0D0
      END IF
      IF (IROW .EQ. -1) FRQ = -FRQ
      IF (IWT .GT. 0) THEN
         IF (.NOT.dIFNAN(WT)) THEN
            IF (WT .LT. 0.0D0) THEN
               CALL E1STI (1, IOBS)
               CALL e1std (1, WT)
               IF (IDO .GT. 0) THEN
                  CALL E1STI (2, ICALL)
                  CALL E1MES (4, 1, 'The weight for row %(I1) of '//
     &                        'X on invocation number %(I2) of this '//
     &                        'routine was %(d1).  Weights must be '//
     &                        'nonnegative.')
               ELSE
                  CALL E1MES (4, 1, 'The weight for row %(I1) of '//
     &                        'X was %(d1).  Weights must be '//
     &                        'nonnegative.')
               END IF
               IGO = 3
               GO TO 9000
            END IF
         END IF
      ELSE
         WT = 1.0D0
      END IF
 9000 CALL E1POP ('dC1WFR ')
      RETURN
      END
