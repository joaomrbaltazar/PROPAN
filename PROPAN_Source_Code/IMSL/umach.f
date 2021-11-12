c
      subroutine umach (n, nunit)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    n, nunit
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    nn, nout
c                                  SPECIFICATIONS FOR SAVE VARIABLES
      integer    iunit(3)
      save       iunit
c                                  SPECIFICATIONS FOR INTRINSICS
c     intrinsic  iabs
      intrinsic  iabs
      integer    iabs
c
      data iunit(1)/5/
      data iunit(2)/6/
      data iunit(3)/6/
c                                  By default the value for standard
c                                  error is the same as standard output.
c                                  This is for backward compatibilty.
c                                  The GRAPHICS library may want a
c                                  different default value (such as 0
c                                  for unix machines).
      nn = iabs(n)
      if (nn.lt.1 .or. nn.gt.3) then
c                                  Error.  Invalid range for N.
         nout = iunit(3)
         write (nout,99999) nn
99999    format (/, ' *** TERMINAL ERROR 5 from UMACH.  The absolute',
     &          /, ' ***          value of the index variable must be'
     &          , /, ' ***          1, 2, or 3.  IABS(N) = ', i6,
     &          '.', /)
         stop
c                                  Check for reset or retrieval
      else if (n .lt. 0) then
c                                  Reset
         iunit(nn) = nunit
      else
c                                  Retrieve
         nunit = iunit(n)
      end if
c
      return
      end
