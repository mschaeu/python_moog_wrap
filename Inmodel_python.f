
      subroutine inmodel_python(aux_data, atm_data)
!******************************************************************************
!     This subroutine reads in the model 
!****************************************************************************** 

      implicit real*8 (a-h,o-z) 
      include 'Atmos.com'
      include 'Linex.com'
      include 'Mol.com'
      include 'Quants.com'
      include 'Factor.com'
      include 'Dummy.com'
      include 'Pstuff.com'
      real*8 element(95), logepsilon(95)
      real*8 kaprefmass(100)
      real*8 bmol(110)
      character list*80, list2*70
      integer append

!     added by MS to perform the internal atmosphere method
      real*4 aux_data(6)
      real*8 atm_data(100,4)

!     define some variables that we might need later
      modelnum = modelnum + 1
      modtype = "KURTYPE"

!     define the number of depth points
      ntau = int(aux_data(1))
!     and check to make sure that we're below 100 depth points
      if (ntau .gt. 100) then
         write (array,1012)
         call prinfo (10)
         stop
      endif

!     define the reference wavelength for the optical depth
      wavref = aux_data(2)
!     and now we define all the atm data
      rhox = atm_data(:,1)
      t = atm_data(:,2)
      pgas = atm_data(:,3)
      ne = atm_data(:,4)

      do i=1,ntau
         rhox(i) = atm_data(i,1)
         t(i) = atm_data(i,2)
         pgas(i) = atm_data(i,3)
         ne(i) = atm_data(i,4)
      enddo

!*****Compute other convenient forms of the temperatures
      do i=1,ntau
          theta(i) = 5040./t(i)
          tkev(i) = 8.6171d-5*t(i)
          tlog(i) = dlog(t(i))
      enddo


!*****Convert from logarithmic Pgas scales, if needed
      if (pgas(ntau)/pgas(1) .lt. 10.) then
         do i=1,ntau                                                    
            pgas(i) = 10.0**pgas(i)
         enddo
      endif


!*****Convert from logarithmic Ne scales, if needed
      if(ne(ntau)/ne(1) .lt. 20.) then
         do i=1,ntau                                                    
            ne(i) = 10.0**ne(i)
         enddo
      endif


!*****Convert from Pe to Ne, if needed
      if(ne(ntau) .lt. 1.0e7) then
         do i=1,ntau                                                    
            ne(i) = ne(i)/1.38054d-16/t(i)
         enddo
      endif


!*****compute the atomic partition functions
      do j=1,95
         elem(j) = dble(j)
         call partfn (elem(j),j)
      enddo


!*****Read the microturbulence (either a single value to apply to 
!     all layers, or a value for each of the ntau layers). 
!     Conversion to cm/sec from km/sec is done if needed
      do i=1,ntau
            vturb(i) = aux_data(3)
      enddo
      if (vturb(1) .lt. 100.) then
         write (moditle(55:62),1008) vturb(1)
         do i=1,ntau
            vturb(i) = 1.0e5*vturb(i)
         enddo
      else
         write (moditle(55:62),1008) vturb(1)/1.0e5
      endif


!*****Read in the abundance data, storing the original abundances in xabu
!*****The abundances not read in explicity are taken from the default
!*****solar ones contained in array xsolar.
      natoms = aux_data(4)
      abscale = aux_data(5)
      write (moditle(63:73),1009) abscale
      xhyd = 10.0**xsolar(1)
      xabund(1) = 1.0
      xabund(2) = 10.0**xsolar(2)/xhyd
      do i=3,95                                                      
         xabund(i) = 10.0**(xsolar(i)+abscale)/xhyd
         xabu(i) = xabund(i)
      enddo
      if (natoms .ne. 0) then
         do i=1,natoms                                                  
            xabund(idint(element(i))) = 10.0**logepsilon(i)/xhyd
            xabu(idint(element(i))) = 10.0**logepsilon(i)/xhyd
         enddo
      endif

!*****Compute the mean molecular weight, ignoring molecule formation
!     in this approximation (maybe make more general some day?)
      wtnum = 0.
      wtden = 0.
      do i=1,95
         wtnum = wtnum + xabund(i)*xam(i)
         wtden = wtden + xabund(i)
      enddo
      wtmol = wtnum/(xam(1)*wtden)
      nomolweight = 0
      if (modtype .eq. 'BEGN      ' .or. modtype .eq. 'NEXTGEN   ') then
         nomolweight = 1
      endif
      if (nomolweight .ne. 1) then
         do i=1,ntau
             molweight(i) = wtmol
         enddo
      endif

!*****Compute the density 
      if (modtype .ne. 'NEXTGEN   ') then
         do i=1,ntau                                                    
            rho(i) = pgas(i)*molweight(i)*1.6606d-24/(1.38054d-16*t(i))
         enddo
      endif


!*****Calculate the fictitious number density of hydrogen
!     Note:  ph = (-b1 + dsqrt(b1*b1 - 4.0*a1*c1))/(2.0*a1)
      iatom = 1
      call partfn (dble(iatom),iatom)
      do i=1,ntau    
         th = 5040.0/t(i)         
         ah2 = 10.0**(-(12.7422+(-5.1137+(0.1145-0.0091*th)*th)*th))
         a1 = (1.0+2.0*xabund(2))*ah2
         b1 = 1.0 + xabund(2)
         c1 = -pgas(i)         
         ph = (-b1/2.0/a1)+dsqrt((b1**2/(4.0*a1*a1))-(c1/a1))
         nhtot(i) = (ph+2.0*ph*ph*ah2)/(1.38054d-16*t(i))
      enddo


!*****Molecular equilibrium called here.
!     First, a default list of ions and molecules is considered. Then the 
!     user's list is read in. A check is performed to see if any of these 
!     species need to be added. If so, then they are appended to the 
!     default list. The molecular equilibrium routine is then called.
!     Certain species are important for continuous opacities and damping
!     calculations - these are read from the equilibrium solution and 
!     saved.


!*****Set up the default molecule list
      if (molset .eq. 0) then
         nmol = 21
      else
         nmol = 57
      endif
      if     (molset .eq. 0) then
         do i=1,110
            amol(i) = smallmollist(i)
         enddo
      elseif (molset .eq. 1) then
         do i=1,110
            amol(i) = largemollist(i)
         enddo
      else
         array = 'molset = 0 or 1 only; I quit!'
         call putasci (29,6)
         stop
      endif


!*****Read in the names of additional molecules to be used in 
!     molecular equilibrium if needed.
      moremol = aux_data(6)
      if (moremol .ne. 0) then
         read (nfmodel,*) (bmol(i),i=1,moremol)
         append = 1
         do k=1,moremol
            do l=1,nmol
               if (nint(bmol(k)) .eq. nint(amol(l))) &
               append = 0
            enddo
            if (append .eq. 1) then 
               nmol = nmol + 1
               amol(nmol) = bmol(k)
            endif
            append = 1
         enddo  
      endif


!*****do the general molecular equilibrium
101   call eqlib


!     In the number density array "numdens", the elements denoted by
!     the first subscripts are named in at the ends of the assignment
!     lines; at present these are the only ones needed for continuous 
!     opacities
!     
      do i=1,ntau
         numdens(1,1,i) = xamol(1,i)
         numdens(1,2,i) = xmol(1,i)
         numdens(2,1,i) = xamol(2,i)
         numdens(2,2,i) = xmol(2,i)
         numdens(3,1,i) = xamol(3,i)
         numdens(3,2,i) = xmol(3,i)
         numdens(4,1,i) = xamol(6,i)
         numdens(4,2,i) = xmol(6,i)
         numdens(5,1,i) = xamol(7,i)
         numdens(5,2,i) = xmol(7,i)
         numdens(6,1,i) = xamol(8,i)
         numdens(6,2,i) = xmol(8,i)
         numdens(7,1,i) = xamol(16,i)
         numdens(7,2,i) = xmol(16,i)
         numdens(8,1,i) = xmol(17,i)
      enddo



!*****SPECIAL NEEDS: for NEWMARCS models, to convert kaprefs to our units
      if (modtype .eq. 'NEWMARCS  ') then
         do i=1,ntau
            kapref(i) = kaprefmass(i)*rho(i)
         enddo
!     SPECIAL NEEDS: for KURUCZ models, to create the optical depth array,
!     and to convert kaprefs to our units
      elseif (modtype .eq. 'KURUCZ    ') then
         first = rhox(1)*kaprefmass(1)
         tottau = rinteg(rhox,kaprefmass,tauref,ntau,first) 
         tauref(1) = first
         do i=2,ntau
            tauref(i) = tauref(i-1) + tauref(i)
         enddo
         do i=1,ntau
            kapref(i) = kaprefmass(i)*rho(i)
         enddo
!     SPECIAL NEEDS: for NEXTGEN models, to convert kaprefs to our units
      elseif (modtype .eq. 'NEXTGEN   ') then
         do i=1,ntau                                                    
            kapref(i) = kaprefmass(i)*rho(i)
         enddo
!     SPECIAL NEEDS: for BEGN models, to convert kaprefs to our units
      elseif (modtype .eq. 'BEGN      ') then
         do i=1,ntau                                                    
            kapref(i) = kaprefmass(i)*rho(i)
         enddo
!     SPECIAL NEEDS: for KURTYPE models, to create internal kaprefs,
!     and to compute taurefs from the kaprefs converted to mass units
      elseif (modtype .eq. 'KURTYPE   ') then
         call opacit (1,wavref)
         do i=1,ntau                                                    
            kaprefmass(i) = kapref(i)/rho(i)
         enddo
         first = rhox(1)*kaprefmass(1)
         tottau = rinteg(rhox,kaprefmass,tauref,ntau,first) 
         tauref(1) = first
         do i=2,ntau
            tauref(i) = tauref(i-1) + tauref(i)
         enddo
!     SPECIAL NEEDS: for NEWMARCS models, to convert kaprefs to our units
      elseif (modtype .eq. 'KUR-PADOVA') then
         do i=1,ntau
            kapref(i) = kaprefmass(i)*rho(i)
         enddo
!     SPECIAL NEEDS: for generic models, to create internal kaprefs,
      elseif (modtype .eq. 'GENERIC   ' .or. &
              modtype .eq. 'WEBMARCS  ' .or. &
              modtype .eq. 'WEB2MARC  ') then
         call opacit (1,wavref)
      endif


!*****Convert from logarithmic optical depth scales, or vice versa.
!     xref will contain the log of the tauref
      if(tauref(1) .lt. 0.) then
         do i=1,ntau                                                    
            xref(i) = tauref(i)
            tauref(i) = 10.0**xref(i)
         enddo
      else
         do i=1,ntau                                                    
            xref(i) = dlog10(tauref(i))
         enddo
      endif


!*****Write information to output files
      if (modprintopt .lt. 1) return
!      write (nf1out,1002) moditle
      do i=1,ntau
         dummy1(i) = dlog10(pgas(i))
         dummy2(i) = dlog10(ne(i)*1.38054d-16*t(i))
      enddo
!      write (nf1out,1003) wavref,(i,xref(i),tauref(i),t(i),dummy1(i), &
!                          pgas(i),dummy2(i),ne(i),vturb(i),i=1,ntau)
!      write (nf1out,1004)
      do i=1,95
         dummy1(i) = dlog10(xabund(i)) + 12.0
      enddo
!      write (nf1out,1005) (names(i),i,dummy1(i),i=1,95)
!      write (nf1out,1006) modprintopt, molopt, linprintopt, fluxintopt
!      write (nf1out,1007) (kapref(i),i=1,ntau)
      return


!*****format statements
2001  format (a10)
2002  format (a80)
2003  format (6d13.0)
1001  format('permitted model types are:'/'KURUCZ, BEGN, ', &
             'KURTYPE, KUR-PADOVA, NEWMARCS, WEBMARCS, NEXTGEN, ', &
             'WEB2MARC, or GENERIC'/ 'MOOG quits!')
1002  format (/'MODEL ATMOSPHERE HEADER:'/a80/)
1003  format ('INPUT ATMOSPHERE QUANTITIES', 10x, &
              '(reference wavelength =',f10.2,')'/3x, 'i', 2x, 'xref', &
              3x, 'tauref', 7x, 'T', 6x, 'logPg', 4x, 'Pgas', &
              6x, 'logPe', 5x, 'Ne', 9x, 'Vturb'/ &
              (i4, 0pf6.2, 1pd11.4, 0pf9.1, f8.3, 1pd11.4, 0pf8.3, &
              1pd11.4, d11.2))
1004  format (/'INPUT ABUNDANCES: (log10 number densities, log H=12)'/ &
             '      Default solar abundances: Asplund et al. 2009')
1005  format (5(3x,a2, '(',i2,')=', f5.2))
1006  format (/'OPTIONS: atmosphere = ', i1, 5x, 'molecules  = ', i1/ &
              '         lines      = ', i1, 5x, 'flux/int   = ', i1)
1007  format (/'KAPREF ARRAY:'/(6(1pd12.4)))
1008  format ('vt=', f5.2)
1009  format (' M/H=', f5.2)
1010  format (13('-'),'MOOG OUTPUT FILE', 10('-'), &
              '(MOOG version from 01/28/09)', 13('-')// &
              'THE MODEL TYPE: ', a10)
1011  format ('   The Rosseland opacities and optical depths have ', &
              'been read in')
1012  format ('HOUSTON, WE HAVE MORE THAN 100 DEPTH POINTS! I QUIT!')


      end

