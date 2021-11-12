c
      logical function difnan (x)
c                                  SPECIFICATIONS FOR ARGUMENTS
      double precision x
c
      difnan = x.ne.x
c
      return
      end
