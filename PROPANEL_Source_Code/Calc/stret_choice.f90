!-----------------------------------------------------------------------------------------------!
SUBROUTINE STRET_CHOICE(NT,DSL,ST1,ST2,ST3,ST4,ITYPE)
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
CHARACTER*80 FILESTRET
INTEGER :: I,NT,ITYPE
DOUBLE PRECISION :: ST1,ST2,ST3,ST4,DSL(NT)
!-----------------------------------------------------------------------------------------------!
IF (ITYPE == 0) THEN
   FILESTRET='stretf.dat'
!* PRINT *,' Reading stretching function from file:',FILESTRET
!* PRINT *
   OPEN(UNIT=60,FILE=FILESTRET,STATUS='UNKNOWN')
   READ(60,*) (DSL(I),I=1,NT)
   DSL( 1)=0.D0
   DSL(NT)=1.D0
   CLOSE(60)
END IF !(ITYPE == 0)
!-----------------------------------------------------------------------------------------------!
IF (ITYPE == 1) THEN
!* PRINT *,' Equidistant stretching function'
!* PRINT *
END IF !(ITYPE == 1)
!-----------------------------------------------------------------------------------------------!
IF (ITYPE == 2) THEN
!* PRINT *,' One-sided stretching function at I=1'
!* PRINT *
END IF !(ITYPE == 2)
!-----------------------------------------------------------------------------------------------!
IF (ITYPE == 3) THEN
!* PRINT *,' One-sided stretching function at I=NX'
!* PRINT *
END IF !(ITYPE == 3)
!-----------------------------------------------------------------------------------------------!
IF (ITYPE == 4) THEN
!* PRINT *,' Two-sided streching function'
!* PRINT *
END IF !(ITYPE == 4)
!-----------------------------------------------------------------------------------------------!
IF (ITYPE == 5) THEN
!* PRINT *,' Geometrical progression ds1 at I=1'
!* PRINT *
END IF !(ITYPE == 5)
!-----------------------------------------------------------------------------------------------!
IF (ITYPE == 6) THEN
!* PRINT *,' Geometrical progression ds1 at I=NX'
!* PRINT *
END IF !(ITYPE == 6)
!-----------------------------------------------------------------------------------------------!
IF (ITYPE == 7) THEN
!* PRINT *,' Cosine function'
!* PRINT *
END IF !(ITYPE == 7)
!-----------------------------------------------------------------------------------------------!
IF (ITYPE == 8) THEN
!* PRINT *,' Half-cosine function with stretching at I=1'
!* PRINT *
END IF !(ITYPE == 8)
!-----------------------------------------------------------------------------------------------!
IF (ITYPE == 9) THEN
!* PRINT *,' Half-cosine function with stretching at I=NX'
!* PRINT *
END IF !(ITYPE == 9)
!-----------------------------------------------------------------------------------------------!
IF (ITYPE == 10) THEN
!* PRINT *,' Double-cosine function'
!* PRINT *
END IF !(ITYPE == 10)
!-----------------------------------------------------------------------------------------------!
IF (ITYPE > 0) CALL DSCAL(DSL,ST1,ST2,ST3,ST4,NT,ITYPE)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE STRET_CHOICE
!-----------------------------------------------------------------------------------------------!
