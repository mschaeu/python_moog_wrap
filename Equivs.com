
!******************************************************************************
!     this is not really a common block, but it is a way to make
!     (in a consistent manner for several routines) definitions
!     of some temporary variables that are equivalenced (assigned,
!     or in other words, use the same memory space as) other variables.
!     This is done so that plotting of large numbers of points (<=50000)
!     for synthetic spectra may be accomplished without allocating
!     large amounts of additional memory for the code.
!******************************************************************************

      real*4      chunk(500000,5)
!      equivalence (chunk(1,1),a(1,1)),
!     .            (chunk(1,2),a(1,21)),
!     .            (chunk(1,3),a(1,41)),
!     .            (chunk(1,4),a(1,61)),
!     .            (chunk(1,5),a(1,81))
 
      real*4      xsyn(500000), dev(500000)
!      equivalence (xsyn,kapnu0(1,1)), (dev,kapnu0(1,21))

      real*4      y(500000), z(500000)
!      equivalence (y,kapnu0(1,41)), (z,kapnu0(1,61))

