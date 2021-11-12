C
      SUBROUTINE E1STD (ID, DVALUE)
C                                  SPECIFICATIONS FOR ARGUMENTS
      integer    ID
      double precision DVALUE
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    I, IBEG, ILEN
      CHARACTER  ARRAY(24), SAVE*25
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      integer    IFINIT
      CHARACTER  BLANK(1)
      SAVE       BLANK, IFINIT
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
      EXTERNAL   E1INIT, E1INPL
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   I1ERIF
      integer    I1ERIF
C
      DATA BLANK/' '/, IFINIT/0/
C                                  INITIALIZE IF NECESSARY
      IF (IFINIT .EQ. 0) THEN
         CALL E1INIT
         IFINIT = 1
      END IF
      IF (DVALUE .EQ. 0.0D0) THEN
         WRITE (SAVE,'(D24.15)') DVALUE
      ELSE
         WRITE (SAVE,'(1PD24.15)') DVALUE
      END IF
      DO 40  I=1, 24
   40 ARRAY(I) = SAVE(I:I)
      IBEG = I1ERIF(ARRAY,24,BLANK,1)
      IF (ID.GE.1 .AND. ID.LE.9) THEN
         ILEN = 25 - IBEG
         CALL E1INPL ('D', ID, ILEN, ARRAY(IBEG))
      END IF
C
      RETURN
      END
