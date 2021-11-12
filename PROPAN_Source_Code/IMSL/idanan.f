c
      integer function idanan (n, dx, incx)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    n, incx
      double precision dx(*)
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, k, ner
c                                  SPECIFICATIONS FOR SUBROUTINES
      external   e1mes, e1pop, e1psh, e1sti
c                                  SPECIFICATIONS FOR FUNCTIONS
      external   difnan
      logical    difnan
c
      ner = 1
      idanan = 0
c                                  Check N
      if (n .lt. 0) then
         call e1psh ('IDANAN')
         call e1sti (1, n)
         call e1mes (5, ner, ' The length of the input vector is'//
     &               ' %(I1).  The length must be non-negative.')
         call e1pop ('IDANAN')
         go to 9000
      end if
c                                  Check INCX
      if (incx .lt. 0) then
         call e1psh ('IDANAN')
         call e1sti (1, incx)
         call e1mes (5, ner+1, ' The length of the incrementation'//
     &               ' of the input vector is %(I1).  The '//
     &               'incrementation must be non-negative.')
         call e1pop ('IDANAN')
         go to 9000
      end if
c
      if (n.ge.1 .and. incx.ge.1) then
         i = 1
         k = 1
c                                  Repeat until NaN is found or the
c                                  vector has been checked
   10    if (difnan(dx(i))) idanan = k
         i = i + incx
         k = k + 1
         if (idanan.eq.0 .and. k.le.n) go to 10
      end if
c
 9000 return
      end
