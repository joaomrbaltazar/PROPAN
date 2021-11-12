c
      subroutine e1prt
c                                  SPECIFICATIONS FOR PARAMETERS
      integer    maxrnm
      parameter  (maxrnm=32)
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, iall, ibeg, ibloc, ibloc2, iend, ier, ihdr, j,
     &           lertyp, loc, locm1, locx, maxloc, maxtmp, mloc, mod,
     &           ncbeg, nloc, nouter
      character  msgtmp(70), string(10), temp*(maxrnm)
c                                  SPECIFICATIONS FOR SAVE VARIABLES
      character  atline(9), blank(1), dbb(3), from(6), msgtyp(8,7),
     &           persla(2), qmark, unknow(8)
      save       atline, blank, dbb, from, msgtyp, persla, qmark,
     &           unknow
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
C                              SPECIFICATIONS FOR COMMON /ERCOM8/
      integer    PROLVL, XXLINE(10), XXPLEN(10), ICALOC(10), INALOC(10)
      COMMON     /ERCOM8/ PROLVL, XXLINE, XXPLEN, ICALOC, INALOC
      SAVE       /ERCOM8/
C                              SPECIFICATIONS FOR COMMON /ERCOM9/
      CHARACTER  XXPROC(10)*31
      COMMON     /ERCOM9/ XXPROC
      SAVE       /ERCOM9/
c                                  SPECIFICATIONS FOR INTRINSICS
c     intrinsic  min0
      intrinsic  min0
      integer    min0
c                                  SPECIFICATIONS FOR SUBROUTINES
      external   c1tic, m1ve, umach
c                                  SPECIFICATIONS FOR FUNCTIONS
      external   e3prt, i1dx, i1erif, iachar, icase
      integer    i1dx, i1erif, iachar, icase
      character  e3prt*1
c
      data msgtyp/'N', 'O', 'T', 'E', ' ', ' ', ' ', ' ', 'A',
     &     'L', 'E', 'R', 'T', ' ', ' ', ' ', 'W', 'A', 'R',
     &     'N', 'I', 'N', 'G', ' ', 'F', 'A', 'T', 'A', 'L',
     &     ' ', ' ', ' ', 'T', 'E', 'R', 'M', 'I', 'N', 'A',
     &     'L', 'W', 'A', 'R', 'N', 'I', 'N', 'G', ' ', 'F',
     &     'A', 'T', 'A', 'L', ' ', ' ', ' '/
      data unknow/'U', 'N', 'K', 'N', 'O', 'W', 'N', ' '/
      data atline/' ', 'a', 't', ' ', 'l', 'i', 'n', 'e', ' '/
      data blank/' '/, from/' ', 'f', 'r', 'o', 'm', ' '/
      data dbb/'.', ' ', ' '/, persla/'%', '/'/
      data qmark/'?'/
c
      if (msglen .le. 0) return
c                                  Retrieve standard error unit
c                                  number.
      call umach (3, nouter)
      maxtmp = 70
      mod = 0
      lertyp = ertype(callvl)
      ihdr = hdrfmt(lertyp)
      if (ihdr .eq. 3) then
         if (xxproc(prolvl)(1:1).eq.qmark .and. xxline(prolvl).eq.0)
     &       then
            ihdr = 1
         end if
      end if
      iend = 0
c                                  The decision was made to print
c                                  the error code for every error
c                                  type.  The if block was replaced:
c                                  if (ihdr.eq.1 .and. ertype(callvl)
c                                  .le.4) then
      if (ihdr .eq. 1) then
         msgtmp(1) = blank(1)
         iend = 1
c                                  CONVERT ERROR CODE INTO CHAR STRING
         call c1tic (ercode(callvl), string, 10, ier)
c                                  LOCATE START OF NON-BLANK CHARACTERS
         ibeg = i1erif(string,10,blank,1)
c                                  M1VE IT TO MSGTMP
         call m1ve (string, ibeg, 10, 10, msgtmp, iend+1,
     &              iend+11-ibeg, maxtmp, ier)
         iend = iend + 11 - ibeg
      end if
      if (ihdr .ne. 2) then
         call m1ve (from, 1, 6, 6, msgtmp, iend+1, iend+6, maxtmp, ier)
         iend = iend + 6
      end if
      if (ihdr .eq. 3) then
c                                  THIS IS A PROTRAN RUN TIME ERROR MSG.
c                                  RETRIEVE THE PROCEDURE NAME
         call m1ve (xxproc(prolvl), 1, xxplen(prolvl), 31, msgtmp,
     &              iend+1, iend+xxplen(prolvl), maxtmp, ier)
         mloc = iend + xxplen(prolvl) + 1
         msgtmp(mloc) = blank(1)
         iend = iend + i1dx(msgtmp(iend+1),xxplen(prolvl)+1,blank,1) -
     &          1
         if (xxline(prolvl) .gt. 0) then
c                                  INSERT ATLINE
            call m1ve (atline, 1, 9, 9, msgtmp, iend+1, iend+9,
     &                 maxtmp, ier)
            iend = iend + 9
c                                  CONVERT PROTRAN GLOBAL LINE NUMBER
            call c1tic (xxline(prolvl), string, 10, ier)
c                                  LOCATE START OF NON-BLANK CHARACTERS
            ibeg = i1erif(string,10,blank,1)
c                                  M1VE GLOBAL LINE NUMBER TO MSGTMP
            call m1ve (string, ibeg, 10, 10, msgtmp, iend+1,
     &                 iend+11-ibeg, maxtmp, ier)
            iend = iend + 11 - ibeg
         end if
      else
c                                  THIS IS EITHER A LIBRARY ERROR MSG
c                                  OR A PROTRAN PREPROCESSOR ERROR MSG
         if (ihdr .eq. 1) then
c                                  THIS IS A LIBRARY ERROR MESSAGE.
c                                  RETRIEVE ROUTINE NAME
            call m1ve (rname(callvl), 1, maxrnm, maxrnm, msgtmp,
     &                 iend+1, iend+maxrnm, maxtmp, ier)
c                                  If second letter is upper case and
c                                  first letter lower case, make the fir
c                                  letter an upper case letter.
            if (iachar(msgtmp(iend+2)).ge.65 .and.
     &          iachar(msgtmp(iend+2)).le.90) then
               if (iachar(msgtmp(iend+1)).ge.97 .and.
     &             iachar(msgtmp(iend+1)).le.122) then
                  msgtmp(iend+1) = e3prt(icase(msgtmp(iend+1)))
               end if
            end if
            msgtmp(iend+maxrnm+1) = blank(1)
            iend = iend + i1dx(msgtmp(iend+1),maxrnm+1,blank,1) - 1
         end if
c                                  ADD DOT, BLANK, BLANK IF NEEDED
         if (i1dx(msgsav,3,dbb,3) .ne. 1) then
            call m1ve (dbb, 1, 3, 3, msgtmp, iend+1, iend+3, maxtmp,
     &                 ier)
            iend = iend + 3
            mod = 3
         end if
      end if
c                                  MSGTMP AND MSGSAV NOW CONTAIN THE
c                                   ERROR MESSAGE IN FINAL FORM.
      ncbeg = 59 - iend - mod
      iall = 0
      ibloc = i1dx(msgsav,msglen,persla,2)
      if (ibloc.ne.0 .and. ibloc.lt.ncbeg) then
         locm1 = ibloc - 1
         loc = ibloc + 1
      else if (msglen .le. ncbeg) then
         locm1 = msglen
         iall = 1
      else
         loc = ncbeg
c                                  CHECK FOR APPROPRIATE PLACE TO SPLIT
   10    continue
         if (msgsav(loc) .ne. blank(1)) then
            loc = loc - 1
            if (loc .gt. 1) go to 10
            loc = ncbeg + 1
         end if
         locm1 = loc - 1
      end if
c                                  NO BLANKS FOUND IN FIRST NCBEG CHARS
      if (lertyp.ge.1 .and. lertyp.le.7) then
         write (nouter,99995) (msgtyp(i,lertyp),i=1,8),
     &                       (msgtmp(i),i=1,iend),
     &                       (msgsav(i),i=1,locm1)
      else
         write (nouter,99995) (unknow(i),i=1,8), (msgtmp(i),i=1,iend),
     &                       (msgsav(i),i=1,locm1)
      end if
      if (iall .eq. 0) then
c                                  PREPARE TO WRITE CONTINUATION OF
c                                    MESSAGE
c
c                                  FIND WHERE TO BREAK MESSAGE
c                                    LOC = NUMBER OF CHARACTERS OF
c                                          MESSAGE WRITTEN SO FAR
   20    locx = loc + 64
         nloc = loc + 1
         ibloc2 = ibloc
         maxloc = min0(msglen-loc,64)
         ibloc = i1dx(msgsav(nloc),maxloc,persla,2)
         if (msgsav(nloc).eq.blank(1) .and. ibloc2.eq.0) nloc = nloc +
     &       1
         if (ibloc .gt. 0) then
c                                  PAGE BREAK FOUND AT IBLOC
            locx = nloc + ibloc - 2
            write (nouter,99996) (msgsav(i),i=nloc,locx)
            loc = nloc + ibloc
            go to 20
c                                  DON'T BOTHER LOOKING FOR BLANK TO
c                                    BREAK AT IF LOCX .GE. MSGLEN
         else if (locx .lt. msglen) then
c                                  CHECK FOR BLANK TO BREAK THE LINE
   30       continue
            if (msgsav(locx) .eq. blank(1)) then
c                                  BLANK FOUND AT LOCX
               write (nouter,99996) (msgsav(i),i=nloc,locx)
               loc = locx
               go to 20
            end if
            locx = locx - 1
            if (locx .gt. nloc) go to 30
            locx = loc + 64
c                                  NO BLANKS FOUND IN NEXT 64 CHARS
            write (nouter,99996) (msgsav(i),i=nloc,locx)
            loc = locx
            go to 20
         else
c                                  ALL THE REST WILL FIT ON 1 LINE
            locx = msglen
            write (nouter,99996) (msgsav(i),i=nloc,locx)
         end if
      end if
c                                  SET LENGTH OF MSGSAV AND PLEN
c                                    TO SHOW THAT MESSAGE HAS
c                                    ALREADY BEEN PRINTED
      msglen = 0
      plen = 1
      if (tracon(lertyp).eq.1 .and. (callvl.gt.2.or.ihdr.gt.1)) then
c                                  SHOW TRACEBACK SINCE (FLAG IS ON) AND
c                                  ((This is LIBRARY .GT. 2 levels deep)
c                                  OR (This is protran)).
         write (nouter,99997)
         do 9005  j=callvl, 1, -1
c                                  If second letter is upper case and
c                                  first letter lower case, make the fir
c                                  letter an upper case letter.
            temp = rname(j)
            if (iachar(temp(2:2)).ge.65 .and. iachar(temp(2:2)).le.90)
     &          then
               if (iachar(temp(1:1)).ge.97 .and.
     &             iachar(temp(1:1)).le.122) then
                  temp(1:1) = e3prt(icase(temp(1:1)))
               end if
            end if
            if (j .gt. 1) then
               if (isuser(j-1)) then
                  write (nouter,99998) temp, ertype(j), ercode(j)
               else
                  write (nouter,99999) temp, ertype(j), ercode(j)
               end if
            else
               write (nouter,99998) temp, ertype(j), ercode(j)
            end if
 9005    continue
      end if
c
      return
99995 format (/, ' *** ', 8a1, ' ERROR', 59a1)
99996 format (' *** ', 9x, 64a1)
99997 format (5x, 'Here is a traceback of subprogram calls',
     &       ' in reverse order:', /, 5x, 'Routine name               '
     &       , ' ', '    Error type  Error code', /, 5x, '------------'
     &       , '  ', '              ', '    ----------  ----------')
99998 format (5x, a32, i6, 6x, i6)
99999 format (5x, a32, i6, 6x, i6, 4x, '(Called internally)')
      end
