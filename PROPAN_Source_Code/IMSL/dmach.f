C
      double precision FUNCTION DMACH (N)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    N
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      double precision RMACH(8)
      SAVE       RMACH
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1POP, E1PSH, E1STI
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    IRMACH(16)
C
      EQUIVALENCE (RMACH, IRMACH)
C                                  DEFINE CONSTANTS
      DATA IRMACH(1)/1048576/
      DATA IRMACH(2)/0/
      DATA IRMACH(3)/2146435071/
      DATA IRMACH(4)/-1/
      DATA IRMACH(5)/1017118720/
      DATA IRMACH(6)/0/
      DATA IRMACH(7)/1018167296/
      DATA IRMACH(8)/0/
      DATA IRMACH(9)/1070810131/
      DATA IRMACH(10)/1352628735/
      DATA IRMACH(11)/2146959360/
      DATA IRMACH(12)/0/
      DATA IRMACH(13)/2146435072/
      DATA IRMACH(14)/0/
      DATA IRMACH(15)/-1048576/
      DATA IRMACH(16)/0/
C
      IF (N.LT.1 .OR. N.GT.8) THEN
         CALL E1PSH ('DMACH ')
         DMACH = RMACH(6)
         CALL E1STI (1, N)
         CALL E1MES (5, 5, 'The argument must be between 1 '//
     &               'and 8 inclusive. N = %(I1)')
         CALL E1POP ('DMACH ')
      ELSE
         DMACH = RMACH(N)
      END IF
C
      RETURN
      END
