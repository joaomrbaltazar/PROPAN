c
      subroutine e1pop (name)
c                                  SPECIFICATIONS FOR ARGUMENTS
      character  name*(*)
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    iertyp, ir
c                                  SPECIFICATIONS FOR SPECIAL CASES
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
c                                  SPECIFICATIONS FOR SUBROUTINES
      external   e1mes, e1prt, e1psh, e1sti, e1stl, i1krl
c                                  SPECIFICATIONS FOR FUNCTIONS
      external   i1kst, n1rgb
      integer    i1kst, n1rgb
c
      if (callvl .le. 1) then
         call e1psh ('E1POP ')
         call e1stl (1, name)
         call e1mes (5, 1, 'Error condition in E1POP.  Cannot pop '//
     &               'from %(L1) because stack is empty.')
         stop
      else if (name .ne. rname(callvl)) then
         call e1stl (1, name)
         call e1stl (2, rname(callvl))
         call e1mes (5, 2, 'Error condition in E1POP.  %(L1) does '//
     &               'not match the name %(L2) in the stack.')
         stop
      else
         iertyp = ertype(callvl)
         if (iertyp .ne. 0) then
c                                  M1VE ERROR TYPE AND ERROR CODE TO
c                                    PREVIOUS LEVEL FOR ERROR TYPES 2-7
            if (iertyp.ge.2 .and. iertyp.le.7) then
               ertype(callvl-1) = ertype(callvl)
               ercode(callvl-1) = ercode(callvl)
            end if
c                                  CHECK PRINT TABLE TO DETERMINE
c                                    WHETHER TO PRINT STORED MESSAGE
            if (iertyp .le. 4) then
               if (isuser(callvl-1) .and. printb(iertyp).eq.1)
     &             call e1prt
            else
               if (printb(iertyp) .eq. 1) call e1prt
            end if
c                                  CHECK STOP TABLE AND ERROR TYPE TO
c                                    DETERMINE WHETHER TO STOP
            if (iertyp .le. 4) then
               if (isuser(callvl-1) .and. stoptb(iertyp).eq.1) then
                  stop
               end if
            else if (iertyp .eq. 5) then
               if (stoptb(iertyp) .eq. 1) then
                  stop
               end if
            else if (hdrfmt(iertyp) .eq. 1) then
               if (isuser(callvl-1)) then
                  if (n1rgb(0) .ne. 0) then
                     stop
                  end if
               end if
            end if
         end if
c                                  SET ERROR TYPE AND CODE
         if (callvl .lt. maxlev) then
            ertype(callvl+1) = -1
            ercode(callvl+1) = -1
         end if
c                                  SET IR = AMOUNT OF WORKSPACE
c                                  ALLOCATED AT THIS LEVEL
         ir = i1kst(1) - ialloc(callvl-1)
         if (ir .gt. 0) then
c                                  RELEASE WORKSPACE
            call i1krl (ir)
            ialloc(callvl) = 0
         else if (ir .lt. 0) then
            call e1sti (1, callvl)
            call e1sti (2, ialloc(callvl-1))
            call e1sti (3, i1kst(1))
            call e1mes (5, 3, 'Error condition in E1POP. '//
     &                  ' The number of workspace allocations at '//
     &                  'level %(I1) is %(I2).  However, the total '//
     &                  'number of workspace allocations is %(I3).')
            stop
         end if
c                                  DECREASE THE STACK POINTER BY ONE
         callvl = callvl - 1
      end if
c
      return
      end