
      subroutine pointcurs
!******************************************************************************
!     This subroutine returns the cursor position upon a click
!******************************************************************************

      include 'Pstuff.com'
      integer ichr


      call sm_graphics
      if     (whichwin .eq. '1of1') then
         call sm_window (1,1,1,1,1,1)
      elseif (whichwin .eq. '2of2') then
         call sm_defvar ('y_gutter','0.0')
         call sm_window (1,2,1,1,1,1)
      endif
      call sm_curs (xplotpos,yplotpos,ichr)
      call sm_gflush
      call sm_alpha

      
      return
      end



