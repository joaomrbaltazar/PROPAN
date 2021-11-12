c
      subroutine e1mes (iertyp, iercod, msgpkd)
      integer    maxmes
      parameter  (maxmes=1024)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    iertyp, iercod
      character  msgpkd*(*)
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    ertyp2, i, ier, iplen, isub, last, len2, loc, m, ms,
     &           nloc, num, pbeg
      character  msgtmp(maxmes)
c                                  SPECIFICATIONS FOR SAVE VARIABLES
      integer    ifinit, nforms
      character  blnk, dbb(3), find(4), forms(9), inref(25), lpar,
     &           ncheck(3), percnt, rpar
      save       blnk, dbb, find, forms, ifinit, inref, lpar, ncheck,
     &           nforms, percnt, rpar
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
c                                  SPECIFICATIONS FOR INTRINSICS
c     intrinsic  len,min0
      intrinsic  len, min0
      integer    len, min0
c                                  SPECIFICATIONS FOR SUBROUTINES
      external   c1tci, e1init, e1prt, e1ucs, m1ve, m1vech
c                                  SPECIFICATIONS FOR FUNCTIONS
      external   i1dx
      integer    i1dx
c
      data forms/'A', 'C', 'D', 'I', 'K', 'L', 'R', 'S', 'Z'/,
     &     nforms/9/
      data percnt/'%'/, lpar/'('/, rpar/')'/, blnk/' '/
      data inref/' ', 'i', 'n', ' ', 'r', 'e', 'f', 'e', 'r',
     &     'e', 'n', 'c', 'e', ' ', 't', 'o', ' ', 'k', 'e',
     &     'y', 'w', 'o', 'r', 'd', ' '/
      data ncheck/'N', '1', '*'/, dbb/'.', ' ', ' '/
      data find/'*', ' ', ' ', '*'/
      data ifinit/0/
c                                  INITIALIZE ERROR TABLE IF NECESSARY
      if (ifinit .eq. 0) then
         call e1init
         ifinit = 1
      end if
c                                  CHECK AND SET ERROR TYPE IF NECESSARY
      if (iertyp .ne. -1) then
         ertype(callvl) = iertyp
      else if (iertyp.lt.-1 .or. iertyp.gt.7) then
         msglen = 51
         call m1vech ('.  Error from E1MES.  Illegal error type'//
     &                ' specified. ', msglen, msgsav, msglen)
         call e1prt
         stop
      end if
c
      ertyp2 = ertype(callvl)
c                                  SET ERROR CODE IF NECESSARY
      if (iercod .gt. -1) ercode(callvl) = iercod
      len2 = len(msgpkd)
c
      if (iertyp.eq.0 .or. iercod.eq.0) then
c                                  REMOVE THE ERROR STATE
         msglen = 0
      else if (len2.eq.0 .or. (len2.eq.1.and.msgpkd(1:1).eq.blnk)) then
         if (ertyp2 .eq. 6) iferr6 = 1
         if (ertyp2 .eq. 7) iferr7 = 1
c                                  UPDATE CHECKSUM PARAMETER ERCKSM
         call e1ucs
c                                  PRINT MESSAGE IF NECESSARY
         if (ertyp2.ge.5 .and. printb(ertyp2).eq.1) call e1prt
      else
c                                  FILL UP MSGSAV WITH EXPANDED MESSAGE
         len2 = min0(len2,maxmes)
         do 10  i=1, len2
            msgtmp(i) = msgpkd(i:i)
   10    continue
         ms = 0
         m = 0
c                                  CHECK PLIST FOR KEYWORD NAME
         nloc = i1dx(plist,plen,ncheck,3)
         if (nloc.gt.0 .and. hdrfmt(ertyp2).eq.3) then
c                                  M1VE INREF INTO MSGSAV
            call m1ve (inref, 1, 25, 25, msgsav, 1, 25, 25, ier)
c                                  GET LENGTH OF KEYWORD NAME
            call c1tci (plist(nloc+3), 3, iplen, ier)
            pbeg = nloc + 3 + ier
c                                  M1VE KEYWORD NAME INTO MSGSAV
            call m1ve (plist, pbeg, pbeg+iplen-1, plen, msgsav, 26,
     &                 iplen+25, maxmes, ier)
c                                  UPDATE POINTER
            ms = iplen + 25
         end if
c                                  INSERT DOT, BLANK, BLANK
         call m1ve (dbb, 1, 3, 3, msgsav, ms+1, ms+3, maxmes, ier)
         ms = ms + 3
c                                  LOOK AT NEXT CHARACTER
   20    m = m + 1
         isub = 0
         if (m .gt. len2-4) then
            last = len2 - m + 1
            do 30  i=1, last
   30       msgsav(ms+i) = msgtmp(m+i-1)
            msglen = ms + last
            go to 40
         else if (msgtmp(m).eq.percnt .and. msgtmp(m+1).eq.lpar .and.
     &           msgtmp(m+4).eq.rpar) then
            call c1tci (msgtmp(m+3), 1, num, ier)
            if (ier.eq.0 .and. num.ne.0 .and. i1dx(forms,nforms,
     &          msgtmp(m+2),1).ne.0) then
c                                  LOCATE THE ITEM IN THE PARAMETER LIST
               call m1ve (msgtmp(m+2), 1, 2, 2, find, 2, 3, 4, ier)
               loc = i1dx(plist,plen,find,4)
               if (loc .gt. 0) then
c                                  SET IPLEN = LENGTH OF STRING
                  call c1tci (plist(loc+4), 4, iplen, ier)
                  pbeg = loc + 4 + ier
c                                  ADJUST IPLEN IF IT IS TOO BIG
                  iplen = min0(iplen,maxmes-ms)
c                                  M1VE STRING FROM PLIST INTO MSGSAV
                  call m1ve (plist, pbeg, pbeg+iplen-1, plen, msgsav,
     &                       ms+1, ms+iplen, maxmes, ier)
                  if (ier.ge.0 .and. ier.lt.iplen) then
c                                  UPDATE POINTERS
                     m = m + 4
                     ms = ms + iplen - ier
c                                  BAIL OUT IF NO MORE ROOM
                     if (ms .ge. maxmes) then
                        msglen = maxmes
                        go to 40
                     end if
c                                  SET FLAG TO SHOW SUBSTITION WAS MADE
                     isub = 1
                  end if
               end if
            end if
         end if
         if (isub .eq. 0) then
            ms = ms + 1
            msgsav(ms) = msgtmp(m)
         end if
         go to 20
   40    ertyp2 = ertype(callvl)
         if (ertyp2 .eq. 6) iferr6 = 1
         if (ertyp2 .eq. 7) iferr7 = 1
c                                  UPDATE CHECKSUM PARAMETER ERCKSM
         call e1ucs
c                                  PRINT MESSAGE IF NECESSARY
         if (ertyp2.ge.5 .and. printb(ertyp2).eq.1) call e1prt
      end if
c                                  CLEAR PARAMETER LIST
      plen = 1
c
      return
      end