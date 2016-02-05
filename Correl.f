
      subroutine correl
!******************************************************************************
!     This routine cross-correlates arrays x and y using any part of the
!     arrays, as long as there are enough data points in each array.
!     The result is printed out in terms of the shift in the abscissa
!     of the y array needed to align it with the x array.
!******************************************************************************

      include 'Atmos.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Plotval.com'
      real*4 xgood(100000), ygood(100000), zgood(100000)
      real*8 spx(1000), spy(1000)
      real*4 q
      real*8 xinterp(21), yinterp(21)


      wavemin = amax1(xsyn(1),xobs(1))
      wavemax = amin1(xsyn(kount),xobs(lount))
      wavecenter = (wavemax + wavemin)/2.


!*****find the array subscripts for the minimum and maximum
!     wavelengths in the synthetic spectrum.
      do i=1,kount
         if (xsyn(i) .ge. wavemin) then
            ilosyn = i
            go to 5
         endif
      enddo
5     do i=ilosyn,kount
         if (xsyn(i) .gt. wavemax) then
            ihisyn = i - 1
            go to 10
         endif
      enddo
      ihisyn = kount


!*****dump the synthetic spectrum wavelengths and fluxes into working arrays
10    itotsyn = ihisyn - ilosyn + 1
      do i=1, itotsyn
         xgood(i) = xsyn(ilosyn-1+i)
         ygood(i) = chunk(ilosyn-1+i,1)
      enddo


!*****find the array subscript for the minimum wavelength in the
!     observed spectrum
      do j=1,lount
         if (xobs(j) .ge. xgood(1)) then
            jloobs = j
            go to 15
         endif
      enddo


!*****using a 3-point Lagrangian, make an interpolated observed
!     spectrum that matches the wavelength step size of the
!     synthetic spectrum
15    j = jloobs
      do i=1,itotsyn
         if (xgood(i) .gt. (xobs(j+1)+xobs(j))/2.) j = j + 1
         q = xgood(i)
         zgood(i) = yobs(j-1)*(q-xobs(j))*(q-xobs(j+1))/ &
                        ((xobs(j-1)-xobs(j))*(xobs(j-1)-xobs(j+1))) + &
                    yobs(j)*(q-xobs(j-1))*(q-xobs(j+1))/ &
                        ((xobs(j)-xobs(j-1))*(xobs(j)-xobs(j+1))) + &
                    yobs(j+1)*(q-xobs(j-1))*(q-xobs(j))/ &
                        ((xobs(j+1)-xobs(j-1))*(xobs(j+1)-xobs(j)))
      enddo


!*****next, EITHER do a single cross-correlation
         call crosscorr (ishift,spy(1),ygood,zgood,itotsyn)

!*****next, do the cross-correlation:
      imax = 2*maxshift + 1
      do i=1,imax
         ishift = i - 1 - maxshift
         spx(i) = dble(ishift)
         call crosscorr (ishift,spy(i),ygood,zgood,itotsyn)
      enddo

!*****interpolate to find the maximum of the correlation function
      corrmax = spy(1)
      mc = 1
      do i=2,imax
         if (spy(i) .gt. corrmax) then
            corrmax = spy(i)
            mc = i
         endif
      enddo
      if (mc.eq.1 .or. mc.eq.imax) then
         deltawave = 0.
         return
      endif
      xinterp(1) = spx(mc-1)
      yinterp(1) = spy(mc-1)
      do j=2,21
         q = spx(mc-1) + (j-1)/10.
         xinterp(j) = q
         yinterp(j) = spy(mc-1)*(q-spx(mc))*(q-spx(mc+1))/ &
                        ((spx(mc-1)-spx(mc))*(spx(mc-1)-spx(mc+1))) + &
                      spy(mc)*(q-spx(mc-1))*(q-spx(mc+1))/ &
                        ((spx(mc)-spx(mc-1))*(spx(mc)-spx(mc+1))) + &
                      spy(mc+1)*(q-spx(mc-1))*(q-spx(mc))/ &
                        ((spx(mc+1)-spx(mc-1))*(spx(mc+1)-spx(mc)))
      enddo
      jmax = 1
      corrmax = yinterp(1)
      do j=2,21
         if (yinterp(j) .gt. corrmax) then
            jmax = j
            corrmax = yinterp(j)
            shiftmax = xinterp(j)
         endif
      enddo


!*****translate this into a correction to the input velocity shift
      deltawave = -step*shiftmax
      deltavel = -3.0d5*step*shiftmax/wavecenter
      veladd = veladd + deltavel


      return
      end






      


