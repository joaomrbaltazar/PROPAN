C
      SUBROUTINE C1TCI (CHRSTR, SLEN, NUM, IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    SLEN, NUM, IER
      CHARACTER  CHRSTR(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    COUNT, I, IMACH5, J, N, S, SIGN
      CHARACTER  ZERO
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      CHARACTER  BLANK, DIGIT*10, MINUS, PLUS
      SAVE       BLANK, DIGIT, MINUS, PLUS
C                                  SPECIFICATIONS FOR EQUIVALENCE
      EQUIVALENCE (DIGIT, ZERO)
C                                  SPECIFICATIONS FOR INTRINSICS
c     INTRINSIC  INDEX
      INTRINSIC  INDEX
      integer    INDEX
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   IMACH
      integer    IMACH
C
      DATA DIGIT/'0123456789'/
      DATA BLANK/' '/, MINUS/'-'/, PLUS/'+'/
C
C                                  CHECK SLEN
      NUM = 0
      IF (SLEN .LE. 0) THEN
         IER = -1
         GO TO 50
      END IF
C                                  HANDLE LEADING BLANKS
      SIGN = 1
      I = 1
   10 IF (I .LE. SLEN) THEN
         IF (CHRSTR(I) .EQ. BLANK) THEN
            I = I + 1
            GO TO 10
         END IF
      ELSE
         IER = 1
         GO TO 50
      END IF
C                                  CHECK FOR SIGN, IF ANY
      S = I
      IF (CHRSTR(I) .EQ. MINUS) THEN
         SIGN = -1
         I = I + 1
      ELSE IF (CHRSTR(I) .EQ. PLUS) THEN
         I = I + 1
      END IF
   20 IF (I .LE. SLEN) THEN
         IF (CHRSTR(I) .EQ. BLANK) THEN
            I = I + 1
            GO TO 20
         END IF
      ELSE
         IER = S
         GO TO 50
      END IF
C                                  SKIP LEADING ZERO
      J = I
   30 IF (I .LE. SLEN) THEN
         IF (CHRSTR(I) .EQ. ZERO) THEN
            I = I + 1
            GO TO 30
         END IF
      ELSE
         IER = 0
         GO TO 50
      END IF
C                                  CHECK FIRST NONBLANK CHARACTER
      COUNT = 0
C                                  CHECK NUMERIC CHARACTERS
      IMACH5 = IMACH(5)
   40 N = INDEX(DIGIT,CHRSTR(I))
      IF (N .NE. 0) THEN
         COUNT = COUNT + 1
         IF (NUM .GT. ((IMACH5-N)+1)/10) THEN
            IER = -2
            GO TO 50
         ELSE
            NUM = NUM*10 - 1 + N
            I = I + 1
            IF (I .LE. SLEN) GO TO 40
         END IF
      END IF
C
      IF (COUNT .EQ. 0) THEN
         IF (I .GT. J) THEN
            IER = I
         ELSE
            IER = S
         END IF
      ELSE IF (I .GT. SLEN) THEN
         NUM = SIGN*NUM
         IER = 0
      ELSE
         NUM = SIGN*NUM
         IER = I
      END IF
C
   50 CONTINUE
      RETURN
      END
