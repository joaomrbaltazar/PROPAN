c
      subroutine dr2ivn (ido, nrow, nvar, x, ldx, intcep, iind,
     &                   indind, idep, inddep, ifrq, iwt, isub, tol,
     &                   b, ldb, r, ldr, d, irank, dfe, scpe, ldscpe,
     &                   nrmiss, xmin, xmax, wk)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    ido, nrow, nvar, ldx, intcep, iind, idep, ifrq, iwt,
     &           isub, ldb, ldr, irank, ldscpe, nrmiss, indind(*),
     &           inddep(*)
      double precision tol, dfe, x(*), b(ldb,*), r(ldr,*), d(*),
     &           scpe(ldscpe,*), xmin(*), xmax(*), wk(*)
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, i1, idepx, igo, intp1, iobs, irow, j, jdepix,
     &           jdepjx, jdepx, k, ldep, ncoef, nconst, ndep, nind,
     &           nobs
      double precision frq, sd2, sparam(5), sumwt, temp, tolsq, wt
c                                  SPECIFICATIONS FOR SAVE VARIABLES
      integer    icall
      save       icall
c                                  SPECIFICATIONS FOR INTRINSICS
c     intrinsic  iabs,dsign,dsqrt
      intrinsic  iabs, dsign, dsqrt
      integer    iabs
      double precision dsign, dsqrt
c                                  SPECIFICATIONS FOR SUBROUTINES
      external   e1mes, e1pop, e1psh, e1sti, e1std, daxpy, dcopy,
     &           drotm, drotmg, dscal, dset, dc1wfr, dcsfrg, dgirts,
     &           dr3ivn
c                                  SPECIFICATIONS FOR FUNCTIONS
      external   dmach, difnan, idanan, n1rty, dxyz
      logical    difnan
      integer    idanan, n1rty
      double precision dmach, dxyz
c
      call e1psh ('dR2IVN ')
c                                  Check for terminal errors.
c
      call dr3ivn (ido, nrow, nvar, ldx, intcep, iind, indind, idep,
     &             inddep, ifrq, iwt, isub, tol, ldb, ldr, ldscpe)
      if (n1rty(0) .ne. 0) go to 9000
      ndep = iabs(idep)
      nind = iabs(iind)
      ncoef = intcep + nind
      intp1 = intcep + 1
      idepx = ncoef + 1
      tolsq = tol**2
      if (ido .le. 1) then
c                                  Initialize ICALL, NRMISS, and DFE.
         icall = 1
         nrmiss = 0
         dfe = 0.0D0
c                                  Initialize elements of R, B, D, and
c                                    SCPE.
         do 10  i=1, ncoef
            call dset (ncoef, 0.0D0, r(1,i), 1)
   10    continue
         do 20  i=1, ndep
            call dset (ncoef, 0.0D0, b(1,i), 1)
   20    continue
         call dset (ncoef, 1.0D0, d, 1)
         do 30  i=1, ndep
            call dset (i, 0.0D0, scpe(1,i), 1)
   30    continue
c                                  Initialize XMIN and XMAX to
c                                    not-a-number.
         call dset (ncoef, dmach(6), xmin, 1)
         call dset (ncoef, dmach(6), xmax, 1)
      else
         icall = icall + 1
      end if
      if (nrow .lt. 0) then
c                                  Rows of data are to be deleted from
c                                    from analysis.
         nobs = -nrow
         irow = -1
      else
c                                  Rows of data are to be added to
c                                    analysis.
         nobs = nrow
         irow = 1
      end if
      if (isub .eq. 0) then
         i1 = 1
      else
         i1 = 2
      end if
      do 130  iobs=1, nobs
c                                  Check frequency and weight.
c
         call dc1wfr (ido, icall, x, ldx, iobs, irow, ifrq, iwt,
     &                dmach(6), nrmiss, frq, wt, igo)
         if (igo .eq. 3) go to 9000
         if (igo .eq. 2) go to 130
         if (igo .eq. 1) go to 130
c                                  Gather independent and dependent
c                                    variables in row IOBS of X into
c                                    WK. Check for missing value code
c                                    not-a-number.
         if (intcep .eq. 1) then
            wk(1) = 1.0D0
         end if
         do 40  i=1, iind
            wk(intcep+i) = x((indind(i)-1)*ldx+iobs)
   40    continue
         do 50  i=1, -iind
            wk(intcep+i) = x((i-1)*ldx+iobs)
   50    continue
         if (idanan(nind,wk(intp1),1) .gt. 0) then
            nrmiss = nrmiss + irow
            go to 130
         end if
         jdepx = idepx
         do 60  i=1, idep
            wk(jdepx) = x((inddep(i)-1)*ldx+iobs)
            jdepx = jdepx + 1
   60    continue
         do 70  i=idep + 1, 0
            wk(jdepx) = x((nvar+i-1)*ldx+iobs)
            jdepx = jdepx + 1
   70    continue
         if (idanan(ndep,wk(idepx),1) .gt. 0) then
            nrmiss = nrmiss + irow
            go to 130
         end if
c                                  Update degrees of freedom.
         dfe = dfe + frq
         if (irow .eq. 1) then
            if (ncoef .gt. 0) then
               if (difnan(xmin(1))) then
c                                  Initialize XMIN and XMAX to first
c                                    nonmissing observation.
c
                  call dcopy (ncoef, wk, 1, xmin, 1)
                  call dcopy (ncoef, wk, 1, xmax, 1)
               end if
            end if
            do 80  i=1, ncoef
c                                  Update XMIN and XMAX.
               temp = wk(i)
               if (temp .lt. xmin(i)) xmin(i) = temp
               if (temp .gt. xmax(i)) xmax(i) = temp
   80       continue
         else
            do 90  i=intp1, ncoef
               temp = wk(i)
               if (temp .eq. xmin(i)) then
                  call e1sti (1, i)
                  call e1std (1, temp)
                  call e1mes (6, 5, 'Downdating is requested, '//
     &                        'but XMIN(%(I1)) equals the current '//
     &                        'value (= %(d1)) of the associated '//
     &                        'regressor variable.  Downdating of '//
     &                        'XMIN cannot occur.')
               end if
               if (temp .eq. xmax(i)) then
                  call e1sti (1, i)
                  call e1std (1, temp)
                  call e1mes (6, 6, 'Downdating is requested, '//
     &                        'but XMAX(%(I1)) equals the current '//
     &                        'value (= %(d1)) of the associated '//
     &                        'regressor variable.  Downdating of '//
     &                        'XMAX cannot occur.')
               end if
   90       continue
         end if
         if (wt .eq. 0.0D0) go to 130
         sd2 = wt*frq
         if (isub .eq. 1) then
c                                  Update sums.
	            sumwt = r(1,1)
            r(1,1) = sumwt + sd2
            if (nind .gt. 0) then
               call daxpy (nind, sd2, wk(2), 1, r(1,2), ldr)
            end if
            call daxpy (ndep, sd2, wk(idepx), 1, b, ldb)
c
c                                  Center variables.
            if (r(1,1) .eq. 0.0D0) then
               go to 130
            end if
            d(1) = 1.0D0/r(1,1)
            if (nind .gt. 0) then
               call daxpy (nind, -d(1), r(1,2), ldr, wk(2), 1)
            end if
            call daxpy (ndep, -d(1), b(1,1), ldb, wk(idepx), 1)
            if (sumwt .eq. 0.0D0) then
               go to 130
            end if
            sd2 = sd2*r(1,1)/sumwt
         end if
c                                  Update D, R, and B via fast Givens
c                                    transformations.
         do 100  i=i1, ncoef
            call drotmg (d(i), sd2, r(i,i), wk(i), sparam)
            if (ndep .gt. 0) call drotm (ndep, b(i,1), ldb, wk(idepx),
     &          1, sparam)
            if (i .eq. ncoef) go to 100
            call drotm (ncoef-i, r(i,i+1), ldr, wk(i+1), 1, sparam)
  100    continue
c                                  Update upper triangle of SCPE.
         jdepjx = idepx
         do 120  j=1, ndep
            jdepix = idepx
            do 110  i=1, j
               scpe(i,j) = scpe(i,j) + wk(jdepix)*sd2*wk(jdepjx)
               jdepix = jdepix + 1
  110       continue
            jdepjx = jdepjx + 1
  120    continue
  130 continue
c
      if (ido.eq.0 .or. ido.eq.3) then
c                                  Remove regressors which are linearly
c                                    dependent.
         nconst = 0
         do 170  i=1, ncoef
            ldep = 0
            if (xmin(i) .eq. xmax(i)) then
               if (xmin(i) .eq. 0.0D0) then
                  ldep = 1
               else
                  nconst = nconst + 1
                  if (nconst .gt. 1) ldep = 1
               end if
            else if (i .gt. intp1) then
c                                  Scale factor for each independent
c                                    variable is SUM(X-XBAR)**2. If
c                                    INTCEP = 0, XBAR is taken to be
c                                    zero in the formula above.
c
               temp = dxyz(i-intcep,r(intp1,i),1,d(intp1),1,r(intp1,i),
     &                1)
               if (d(i)*r(i,i)**2 .le. tolsq*temp) ldep = 1
            end if
            if (ldep .eq. 1) then
               do 140  j=i + 1, ncoef
                  call drotmg (d(j), d(i), r(j,j), r(i,j), sparam)
                  if (ndep .gt. 0) call drotm (ndep, b(j,1), ldb, b(i,
     &                1), ldb, sparam)
                  if (j .eq. ncoef) go to 140
                  call drotm (ncoef-j, r(j,j+1), ldr, r(i,j+1), ldr,
     &                        sparam)
  140          continue
               do 160  j=1, ndep
                  do 150  k=1, j
                     scpe(k,j) = scpe(k,j) + b(i,k)*d(i)*b(i,j)
  150             continue
  160          continue
               call dset (ndep, 0.0D0, b(i,1), ldb)
               call dset (ncoef-i+1, 0.0D0, r(i,i), ldr)
            end if
  170    continue
c                                  Compute B by back substitution.
c
         call dgirts (ncoef, r, ldr, ndep, b, ldb, 1, irank, b, ldb,
     &                r, ldr)
         dfe = dfe - irank
         if (dfe .le. 0.0D0) then
            call e1std (1, dfe)
            call e1mes (6, 4, 'DFE = %(d1).  Statistical '//
     &                  'inference is not possible.  More '//
     &                  'observations are needed.')
         end if
c                                  Set R = SQRT(DIAG(D))*R and then put
c                                    diag(D) = 1.
         do 180  i=1, ncoef
            d(i) = dsign(dsqrt(d(i)),r(i,i))
  180    continue
         do 190  i=1, ncoef
            call dscal (ncoef-i+1, d(i), r(i,i), ldr)
  190    continue
         call dset (ncoef, 1.0D0, d, 1)
      end if
c                                  Copy upper triangle of SCPE into
c                                    lower triangle.
      if (ndep .gt. 0) then
         call dcsfrg (ndep, scpe, ldscpe)
      end if
c                                  Exit section
 9000 call e1pop ('dR2IVN ')
      return
      end
