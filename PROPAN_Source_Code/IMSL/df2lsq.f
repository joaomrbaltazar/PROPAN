C
      SUBROUTINE dF2LSQ (F, INTCEP, NBASIS, NDATA, XDATA, FDATA, IWT,
     &                   WEIGHT, A, SSE, WK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    INTCEP, NBASIS, NDATA, IWT
      double precision F, SSE, XDATA(*), FDATA(*), WEIGHT(*), A(*),
     &           WK(*)
      EXTERNAL   F
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    I, ID, IDEP, IDO, IDUMMY(1), IFRQ, IIND, IMAX, IMIN,
     &           IPRNT, IR, IRANK, ISTOP, ISUB, IWT2, IX, IXX, J, LDB,
     &           LDR, LDSCPE, LDX, NCOL, NRMISS, NROW
      double precision DFE, SCPE(1), TOL
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1POP, E1POS, E1PSH, E1STI, E1USR, dR2IVN
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   dmach, N1RCD, N1RTY
      integer    N1RCD, N1RTY
      double precision dmach
C
      CALL E1PSH ('dF2LSQ ')
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
C                                  Check NDATA
      IF (NDATA .LT. 1) THEN
         CALL E1STI (1, NDATA)
         CALL E1MES (5, 3, 'The number of data points must be '//
     &               'at least one while NDATA = %(I1) is given.')
      END IF
C                                  Check IWT
      IF (IWT.LT.0 .OR. IWT.GT.1) THEN
         CALL E1STI (1, IWT)
         CALL E1MES (5, 4, 'The weighting option must be a zero or '//
     &               'a one while IWT = %(I1) is given.')
      END IF
C
      IF (N1RTY(0) .NE. 0) GO TO 9000
C
      NROW = 1
      NCOL = NBASIS + 1 + IWT
      LDX = 1
      IIND = -NBASIS
      IDEP = -1
      IFRQ = 0
      IF (IWT .EQ. 0) THEN
         IWT2 = 0
      ELSE
         IWT2 = NBASIS + 1
      END IF
      ISUB = INTCEP
      TOL = 100.0D0*dmach(4)
      LDB = INTCEP + NBASIS
      LDR = INTCEP + NBASIS
      LDSCPE = 1
C                                  Partition workspace
      IR = 1
      IX = IR + (NBASIS+INTCEP)**2
      ID = IX + INTCEP + NCOL
      IMIN = ID + NBASIS + INTCEP
      IMAX = IMIN + NBASIS + INTCEP
C
      DO 20  I=1, NDATA
         IF (I .EQ. 1) THEN
            IDO = 1
         ELSE IF (I .EQ. NDATA) THEN
            IDO = 3
         ELSE
            IDO = 2
         END IF
         IXX = IX + INTCEP
         DO 10  J=1, NBASIS
            CALL E1USR ('ON')
            WK(IXX) = F(J,XDATA(I))
            CALL E1USR ('OFF')
            IXX = IXX + 1
   10    CONTINUE
         IF (IWT .EQ. 1) THEN
            WK(IXX) = WEIGHT(I)
            IXX = IXX + 1
         END IF
         WK(IXX) = FDATA(I)
C                                  Turn off printing and stopping of
C                                  type 6 errors.
C                                  First retreive current settings.
         CALL E1POS (-6, IPRNT, ISTOP)
C                                  Then set new values
         CALL E1POS (6, 0, 0)
         CALL dR2IVN (IDO, NROW, NCOL, WK(IX+INTCEP), LDX, INTCEP,
     &                IIND, IDUMMY, IDEP, IDUMMY, IFRQ, IWT2, ISUB,
     &                TOL, A, LDB, WK(IR), LDR, WK(ID), IRANK, DFE,
     &                SCPE, LDSCPE, NRMISS, WK(IMIN), WK(IMAX), WK(IX))
C                                  Reset old values
         CALL E1POS (6, IPRNT, ISTOP)
         IF (N1RTY(0) .EQ. 4) GO TO 30
   20 CONTINUE
C                                  Clear warnings
      IF (N1RTY(1).EQ.3 .OR. N1RTY(1).EQ.6) CALL E1MES (0, 0,
     &    ' ')
      IF (IRANK .LT. INTCEP+NBASIS) THEN
         IF (INTCEP .EQ. 0) THEN
            CALL E1MES (3, 1, 'Linear dependence of the basis '//
     &                  'functions was declared.  Appropriate '//
     &                  'elements of A are set to zero.')
         ELSE
            CALL E1MES (3, 2, 'Linear dependence of the constant '//
     &                  'function (the intercept) and the basis '//
     &                  'functions was declared.  Appropriate '//
     &                  'elements of A are set to zero.')
         END IF
      END IF
   30 IF (N1RTY(1).EQ.4 .AND. N1RCD(1).EQ.1) THEN
         CALL E1MES (4, 1, 'An element of the weight vector is not '//
     &               'positive.  All elements must be positive.')
      END IF
C                                  Exit section
 9000 CALL E1POP ('dF2LSQ ')
      SSE = SCPE(1)
      RETURN
      END