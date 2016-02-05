
      real*8 function partnew (atom,k,level)
!******************************************************************************
!     This routine computes partition functions for those species that have
!     been updated from those that are in ATLAS9.  The polynomial
!     representations are in the form used by Irwin (1981, ApJS, 45, 621):
!            log10(U) = SUM{C_j*log10(T)**(j-1)}    for j = 1-->6
!     where T = temperature and C_j are the six polynomial coefficients.
!     Information on the updated partition functions are given by
!     Lawler & Sneden (2002, in preparation).  This routine will be
!     increasingly invoked as new new data are added to the *newpartdata*
!     array,  until a full replacement for the older partition functions
!     is implemented.
!******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Quants.com'


      iatom = nint(atom)
      iarray = partflag(iatom,k)

      if (level .gt. 500) then
         temp = dlog(dble(level))
      else
         temp = tlog(level)
      endif

      ulog = 0.
      do j=1,6
         ulog = ulog + newpartdata(iarray,j)*temp**(j-1)
      enddo
      partnew = dexp(ulog)

      return
      end                                  



