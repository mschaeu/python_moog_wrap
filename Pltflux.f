 
      subroutine pltflux
!******************************************************************************
!     This subroutine controls the decisions that are made around the
!     plots of flux curves
!******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Pstuff.com'


!  call up the flux plot
      if (plotopt .eq. 0) return
10    choice = 'y'
      plotroutine = 'term_land_flux'
      lscreen = 12
      call makeplot (lscreen)


!  make a hardcopy, write to a postscript file, or replot?
      array = 'WHAT TO DO NEXT ([n]/h/f/r)? '
      lscreen = 12
      nchars = 33
      call getasci (nchars,lscreen)
      choice = chinfo(1:1)
      if (choice.eq.'n' .or. nchars.le.0) then
         return
      elseif (choice .eq. 'h') then
         plotroutine = 'hard_land_flux'
         call makeplot (lscreen)
      elseif (choice .eq. 'r') then
         go to 10
      elseif (choice .eq. 'f') then
         plotroutine = 'file_land_flux'
         call makeplot (lscreen)
      endif

      return
      end




