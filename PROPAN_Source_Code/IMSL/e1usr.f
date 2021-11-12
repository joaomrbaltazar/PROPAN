C
      SUBROUTINE E1USR (SWITCH)
C                                  SPECIFICATIONS FOR ARGUMENTS
      CHARACTER  SWITCH*(*)
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      integer    IFINIT
      SAVE       IFINIT
C                                  SPECIFICATIONS FOR SPECIAL CASES
C                              SPECIFICATIONS FOR COMMON /ERCOM1/
      integer    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51),
     &           PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7,
     &           IALLOC(51), HDRFMT(7), TRACON(7)
      COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE,
     &           PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT,
     &           TRACON
      SAVE       /ERCOM1/
C                              SPECIFICATIONS FOR COMMON /ERCOM2/
      CHARACTER  MSGSAV(1024), PLIST(300), RNAME(51)*32
      COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
      SAVE       /ERCOM2/
C                              SPECIFICATIONS FOR COMMON /ERCOM3/
      double precision ERCKSM
      COMMON     /ERCOM3/ ERCKSM
      SAVE       /ERCOM3/
C                              SPECIFICATIONS FOR COMMON /ERCOM4/
      logical    ISUSER(51)
      COMMON     /ERCOM4/ ISUSER
      SAVE       /ERCOM4/
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1INIT, E1MES, E1STL
C
      DATA IFINIT/0/
C                                  INITIALIZE ERROR TABLE IF NECESSARY
      IF (IFINIT .EQ. 0) THEN
         CALL E1INIT
         IFINIT = 1
      END IF
      IF (SWITCH.EQ.'ON' .OR. SWITCH.EQ.'on') THEN
         ISUSER(CALLVL) = .TRUE.
      ELSE IF (SWITCH.EQ.'OFF' .OR. SWITCH.EQ.'off') THEN
         ISUSER(CALLVL) = .FALSE.
      ELSE
         CALL E1STL (1, SWITCH)
         CALL E1MES (5, 1, 'Invalid value for SWITCH in call to'//
     &               ' E1USR.  SWITCH must be set to ''ON'' or '//
     &               '''OFF''.  SWITCH = ''%(L1)'' ')
      END IF
C
      RETURN
      END
