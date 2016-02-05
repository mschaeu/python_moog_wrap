
      subroutine number (nchar,line,xnum)
!******************************************************************************
!     this routine decodes a character string into a double precision
!     floating point number
!******************************************************************************
 
      include 'Pstuff.com'
      real*8 xnum
      character form*10

 
!*****if a carriage return has been hit, return with -9999.
      if (nchar .le. 0) then
         xnum = -9999.
         return
      endif


!*****set the conversion format
      if (nchar .lt. 10) then
         write(form,1001) nchar
      else
         write(form,1002) nchar
      endif

 
!*****now do the conversiton to a number
      read (unit=chinfo,fmt=form,iostat=ierr,err=100) xnum
      return

 
!*****here an error has been detected
100   write (errmess,1004) ierr,chinfo(1:nchar)
      nchars = 65
      istat = ivwrite (line,3,errmess,nchars)
      xnum = -9999.
      return      


!*****format statements
1001  format('(f',i1,'.0)    ')
1002  format('(f',i2,'.0)   ')
1004  format ('ERROR IN NUMBER INPUT: ERROR=',i4,5x,'NUMBER=',a20) 


      end






