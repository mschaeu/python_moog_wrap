
      subroutine damping (linnumber)
!******************************************************************************
!     This subroutine computes damping 'gamma' factors
!     and then Voigt parameters 'a'
!******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'


      j = linnumber
      iwave = int(wave1(j))
      iatom10 = nint(10.*atom1(j))
      if (dampnum(j) .lt. 0.) dampnum(j) = 10.**dampnum(j)


!*****for a few lines, explicit detailed broadening terms have 
!     appeared in the literature, and so do these lines with a 
!     sepaarate subroutine
      if (itru .eq. 0) then
!     CH
         if    (iatom10 .eq. 1060) then
            if (iwave .eq. 3693) then
               call trudamp (j)
               damptype(j) = 'TRUEgam'
               return
            endif
!     Ca I
         elseif (iatom10 .eq. 200) then
            if (iwave.eq.6717 .or. iwave.eq.6318 &
                .or. iwave.eq.6343 .or. iwave.eq.6361) then
               call trudamp (j)
               damptype(j) = 'TRUEgam'
               return
            endif
!     Ca I autoionization
         elseif (iatom10 .eq. 200) then
            if (iwave.eq.6318 .or. &
                iwave.eq.6343 .or. &
                iwave.eq.6361) then
               call trudamp (j)
               damptype(j) = 'TRUEgam'
               return
            endif
         endif
      endif
  

!*****here are the calculations to set up the damping; for atomic lines 
!     there are several options:
!        dampingopt = 0 and dampnum = 0 ---> 
!                             c6 = Unsold formula
!        dampingopt = 0 and dampnum < 10^(-15) ---> 
!                             c6 = dampnum
!        dampingopt = 0 and 10^(-15) < dampnum < 10^(-5) ---> 
!                             gamma = dampnum
!        dampingopt = 0 and dampnum(i) > 10^(-5) ---> 
!                             c6 =  (Unsold formula)*dampnum
!        dampingopt = 1 --->    
!                             gammav = gamma_Barklem if possible, 
!                                      otherwise use dampingopt=0 options
!        dampingopt = 2 ---> 
!                             c6 = c6_Blackwell-group
!        dampingopt = 3 and dampnum <= 10^(-10) --->             
!                             c6 = c6_NEXTGEN for H I, He I, H2
!        dampingopt = 3 and dampnum > 10^(-10) --->             
!                             c6 = (c6_NEXTGEN for H I, He I, H2)*dampnum
!     for molecular lines (lacking a better idea) --->
!                                        c6 done as in dampingopt = 0


!*****these damping calculations are done at each atmosphere level
!      if (linprintopt .gt. 2) write (nf1out,1001) j, wave1(j)
      do i=1,ntau
         ich = nint(charge(j))
         v1 = dsqrt(2.1175d8*t(i)*(1.0/amass(j)+1.008))


!*****first calculate an Unsold approximation to gamma_VanderWaals
         if (atom1(j) .gt. 100.) then
            ebreakup = 7.0
         else 
            ebreakup = chi(j,ich)
         endif
         if (e(j,1).ge.ebreakup .or. e(j,2).ge.ebreakup) then
            unsold = 1.0e-33
         else
            unsold = dabs(1.61d-33*(13.598*charge(j)/(ebreakup - &
                     e(j,1)))**2 - 1.61d-33*(13.598*charge(j)/ &
                     (ebreakup-e(j,2)))**2)
         endif


!*****dampingopt = 0 or 
!*****dampingopt = 1 and no Barklem data
         if     (dampingopt .eq. 0 .or. &
                (dampingopt.eq.1 .and. gambark(j).lt.0)) then
            if     (dampnum(j) .eq. 0.0) then
               damptype(j) = 'UNSLDc6'
               gammav = 17.0*unsold**0.4*v1**0.6*numdens(1,1,i)
            elseif (dampnum(j) .lt. 1.0d-15) then
               damptype(j) = '   MYc6'
               gammav = 17.0*dampnum(j)**0.4*v1**0.6*numdens(1,1,i)
            elseif (dampnum(j) .lt. 1.0d-04) then
               damptype(j) = 'MYgamma'
               gammav = dampnum(j)*(t(i)/10000.)**0.3*numdens(1,1,i)
            else
               damptype(j) = 'MODUNc6'
               gammav = &
                  17.0*(unsold*dampnum(j))**0.4*v1**0.6*numdens(1,1,i)
            endif


!*****dampingopt = 1 with extant Barklem data
         elseif (dampingopt.eq.1 .and. gambark(j).gt.0.) then
               damptype(j) = 'BKgamma'
               gammav = &
                  gambark(j)*(t(i)/10000.)**alpbark(j)*numdens(1,1,i)


!*****dampingopt = 2
         elseif (dampingopt .eq. 2) then
            damptype(j) = 'BLKWLc6'
            gammav = 17.0*((1.0 + 0.67*e(j,1))*unsold)**0.4* &
                     v1**0.6*numdens(1,1,i)


!*****dampingopt = 3
         elseif (dampingopt .eq. 3) then
            damptype(j) = 'NXTGNc6'
            if (dampnum(j) .le. 1.0d-10) dampnum(j) = 1.0
                 c6h = dabs(1.01d-32*(charge(j)**2)* &
                   (13.598/(ebreakup - e(j,1)))**2 - 1.61d-33* &
                   (13.598/(ebreakup-e(j,2)))**2)
                 c6he = dabs((0.204956/0.666793)*1.01d-32* &
                        (charge(j)**2)*(13.598/(ebreakup - &
                        e(j,1)))**2 - 1.61d-33*(13.598/(ebreakup- &
                        e(j,2)))**2)
                 c6ht = dabs((0.806/0.666793)*1.01d-32* &
                        (charge(j)**2)*(13.598/(ebreakup - &
                        e(j,1)))**2 - 1.61d-33*(13.598/(ebreakup- &
                        e(j,2)))**2)
               gammav = 17.0*v1**0.6*(c6h**0.4*numdens(1,1,i) + &
                        c6he**0.4*numdens(2,1,i) + &
                        c6ht**0.4*numdens(8,1,i))*dampnum(j)**0.4
         endif


!*****compute radiative broadening either by an approximate formula or 
!*****the value in Barklem.dat) 
         if (gamrad(j).ne.0.0 .and. dampingopt .eq. 1) then
            gammar = gamrad(j)
         else
            gammar = 2.223d15/wave1(j)**2
         endif
        

!*****now Stark broadening (approximate formulae)
         excdiff = chi(j,nint(charge(j))) - e(j,2)
         if (excdiff .gt. 0.0 .and. atom1(j).lt.100.) then
            effn2 = 13.6*charge(j)**2/excdiff
         else
            effn2 = 25.
         endif
         gammas = 1.0e-8*ne(i)*effn2**2.5


!*****now finish by summing the gammas and computing the Voigt *a* values
         gammatot = gammar + gammas + gammav
         a(j,i) = gammatot*wave1(j)*1.0d-8/(12.56636*dopp(j,i))
!         if (linprintopt .gt. 2) write (nf1out,1002) i, gammar, &
!            gammas, gammav, gammatot, a(j,i)
      enddo
      return


!*****format statements
1001  format(//' LINE BROADENING PARAMETERS FOR LINE', i4, &
             ' AT WAVELENGTH',f8.2/ &
             '  i',4x,'natural',6x,'Stark',4x,'VdWaals', &
             6x,'total',5x,'a(j,i)')
1002  format (i3,1p5e11.3)


      end




