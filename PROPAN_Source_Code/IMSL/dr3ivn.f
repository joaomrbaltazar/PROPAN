c
      subroutine dr3ivn (ido, nrow, ncol, ldx, intcep, iind, indind,
     &                   idep, inddep, ifrq, iwt, isub, tol, ldb, ldr,
     &                   ldscpe)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    ido, nrow, ncol, ldx, intcep, iind, idep, ifrq, iwt,
     &           isub, ldb, ldr, ldscpe, indind(*), inddep(*)
      double precision tol
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, ndep, ner, nind
c                                  SPECIFICATIONS FOR INTRINSICS
c     intrinsic  iabs
      intrinsic  iabs
      integer    iabs
c                                  SPECIFICATIONS FOR SUBROUTINES
      external   c1dim, c1iarg, c1ind, e1mes, e1sti, e1std
c                                  SPECIFICATIONS FOR FUNCTIONS
      external   n1rty
      integer    n1rty
c
      ner = 1
      call c1iarg (ido, 'IDO', 0, 3, ner)
      call c1iarg (ncol, 'NCOL', 0, -1, ner)
      call c1dim (0, iabs(nrow), 'IABS(NROW)', ldx, 'LDX', ner)
      call c1iarg (intcep, 'INTCEP', 0, 1, ner)
      call c1ind (0, ifrq, 'IFRQ', ncol, 'NCOL', ner)
      call c1ind (0, iwt, 'IWT', ncol, 'NCOL', ner)
      if (tol.lt.0.0D0 .or. tol.gt.1.0D0) then
         call e1std (1, tol)
         call e1mes (5, ner, 'TOL = %(d1).  TOL must be between 0.0 '//
     &               'and 1.0, inclusive.')
      end if
      ner = ner + 1
      if (intcep .eq. 0) then
         if (isub .ne. 0) then
            call e1sti (1, isub)
            call e1mes (5, ner, 'INTCEP = 0 and ISUB = %(I1).  '//
     &                  'When INTCEP = 0, ISUB must equal 0.')
         end if
         ner = ner + 1
      else
         call c1iarg (isub, 'ISUB', 0, 1, ner)
      end if
      nind = iabs(iind)
      ndep = iabs(idep)
      if (ndep .gt. 0) then
         call c1dim (0, intcep+nind, 'INTCEP+IABS(IIND)', ldb,
     &               'LDB', ner)
         call c1dim (0, intcep+nind, '*INTCEP+IABS(IIND)', ldr,
     &               'LDR', ner)
         call c1dim (0, ndep, 'IABS(IDEP)', ldscpe, 'LDSCPE', ner)
      else
         ner = ner + 3
      end if
      if (n1rty(0) .ne. 0) go to 9000
      do 10  i=1, idep
         if (inddep(i).le.0 .or. inddep(i).gt.ncol) then
            call e1sti (1, i)
            call e1sti (2, inddep(i))
            call e1sti (3, ncol)
            call e1mes (5, ner, 'INDDEP(%(I1)) = %(I2) and '//
     &                  'NCOL = %(I3).  INDDEP(%(I1)) must be '//
     &                  'greater than or equal 1 and less than or '//
     &                  'equal to NCOL.')
         end if
   10 continue
      ner = ner + 1
      do 20  i=1, iind
         if (indind(i).le.0 .or. indind(i).gt.ncol) then
            call e1sti (1, i)
            call e1sti (2, indind(i))
            call e1sti (3, ncol)
            call e1mes (5, ner, 'INDIND(%(I1)) = %(I2) and '//
     &                  'NCOL = %(I3).  INDIND(%(I1)) must be '//
     &                  'greater than or equal 1 and less than or '//
     &                  'equal to NCOL.')
         end if
   20 continue
 9000 return
      end
