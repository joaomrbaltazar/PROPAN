!-----------------------------------------------------------------------------------------------!
!  Subroutine SDIV performs "safe division", that is to prevent overflow, underflow, NaN, or    !
!  infinity errors.  An alternate value is returned if the division cannot be performed.        !
!                                                                                               !
!  For more information, see the discussion on:                                                 !
!  http://wiki.seas.harvard.edu/geos-chem/index.php/Floating_point_math_issues                  !
!-----------------------------------------------------------------------------------------------!
DOUBLE PRECISION FUNCTION SDIV(N,D,ALTV)
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
DOUBLE PRECISION :: N,D,ALTV
!-----------------------------------------------------------------------------------------------!
IF ( (EXPONENT(N) - EXPONENT(D) >= MAXEXPONENT(N)) .OR. (D == 0.D0) ) THEN
   SDIV = ALTV
ELSE
   SDIV = N / D
ENDIF
!-----------------------------------------------------------------------------------------------!
END FUNCTION SDIV
!-----------------------------------------------------------------------------------------------!
