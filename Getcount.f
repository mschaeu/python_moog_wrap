
      subroutine getcount (num,linechars)
!***************************************************************************
!     this routine simply counts the number of characters in the
!     ascii array 'linechars' until the first blank is seen;
!     it is most useful for discovering the length of a file name.
!***************************************************************************

      character*80 linechars
      integer num


      do i=num,1,-1
         if (linechars(i:i) .ne. ' ') go to 11
      enddo
      num = -1
      return
11    num = i


      return
      end


