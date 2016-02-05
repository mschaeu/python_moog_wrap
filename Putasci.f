
      subroutine putasci (num,line)
!******************************************************************************
!     this routine prints out the characters contained in 'array'
!******************************************************************************

      include 'Pstuff.com'

      istat = ivmove(line-1,1)
      istat = ivcleol()
      istat = ivmove(line-1,1)
      write (errmess,1001) num
      write (*,errmess) array
      return


!*****format statements
1001  format ('(a',i2,'$)')


      end







