
      subroutine synth                   
!******************************************************************************
!     This program synthesizes a section of spectrum and compares it
!     to an observation file.
!******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Factor.com'
      include 'Mol.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Dummy.com'


!*****examine the parameter file
      call params

!*****open the files for: standard output, raw spectrum depths, smoothed
!     spectra, and (if desired) IRAF-style smoothed spectra
      nf1out = 20     
      lscreen = 4
      array = 'STANDARD OUTPUT'
      nchars = 15
      call infile ('output ',nf1out,'formatted  ',0,nchars, &
                   f1out,lscreen)
      nf2out = 21               
      lscreen = lscreen + 2
      array = 'RAW SYNTHESIS OUTPUT'
      nchars = 20
      call infile ('output ',nf2out,'formatted  ',0,nchars, &
                   f2out,lscreen)


      nf3out = 22
      lscreen = lscreen + 2
      array = 'SMOOTHED SYNTHESES OUTPUT'
      nchars = 25
      call infile ('output ',nf3out,'formatted  ',0,nchars, &
                   f3out,lscreen)

      if (f5out .ne. 'optional_output_file') then
         nf5out = 26
         lscreen = lscreen + 2
         array = 'POSTSCRIPT PLOT OUTPUT'
         nchars = 22
         call infile ('output ',nf5out,'formatted  ',0,nchars, &
                      f5out,lscreen)
      endif

      if (iraf .ne. 0) then
         nf4out = 23               
         lscreen = lscreen + 2
         array = 'IRAF ("rtext") OUTPUT'
         nchars = 24
         call infile ('output ',nf4out,'formatted  ',0,nchars, &
                       f4out,lscreen)
      endif

!*****open and read the model atmosphere file
      nfmodel = 30
      lscreen = lscreen + 2
      array = 'THE MODEL ATMOSPHERE'
      nchars = 20
      call infile ('input  ',nfmodel,'formatted  ',0,nchars, &
                   fmodel,lscreen)
      call inmodel


!*****open the line list file and the strong line list file
      nflines = 31
      lscreen = lscreen + 2
      array = 'THE LINE LIST'
      nchars = 13
      call infile ('input  ',nflines,'formatted  ',0,nchars, &
                    flines,lscreen)
      if (dostrong .gt. 0) then
         nfslines = 32
         lscreen = lscreen + 2
         array = 'THE STRONG LINE LIST'
         nchars = 20
         call infile ('input  ',nfslines,'formatted  ',0,nchars, &
                       fslines,lscreen)
      endif

!*****do the syntheses
      choice = '1'
      do i=1,100
         if (i .eq. 100) then
            write (*,1002) 
            stop
         endif
         ncall = 1
         call getsyns (lscreen,ncall)

!*****now either don't make a plot (plotopt = 0)
!                plot the synthetic spectrum, (plotopt = 1)
!                plot syntheses and observation (plotopt = 2)
!                or just smooth the syntheses (plotopt = 3)
         if (choice .eq. 'n') then
            ncall = 2
         else
            ncall = 1
         endif
         if     (plotopt .eq. 0) then
            choice = 'q'
         elseif (plotopt .eq. 1) then
!            call pltspec (lscreen,ncall)
         elseif (plotopt .eq. 2) then
            nfobs = 33               
            lscreen = lscreen + 2
            array = 'THE OBSERVED SPECTRUM'
            nchars = 21
            if (specfileopt .eq. 1) then
               call infile ('input  ',nfobs,'unformatted',2880,nchars, &
                            fobs,lscreen)
            else
               call infile ('input  ',nfobs,'formatted  ',0,nchars, &
                            fobs,lscreen)
            endif
            choice = 'q'
!            call pltspec (lscreen,ncall)
         elseif (plotopt .eq. 3) then
            nfobs = 33               
            lscreen = lscreen + 2
            array = 'THE OBSERVED SPECTRUM'
            nchars = 21
            if (specfileopt .eq. 1) then
               call infile ('input  ',nfobs,'unformatted',2880,nchars, &
                            fobs,lscreen)
            else
               call infile ('input  ',nfobs,'formatted  ',0,nchars, &
                            fobs,lscreen)
            endif

            !*****begin by reading in an observed spectrum
            call readobs (lscreen)
            if (lount .eq. -1) then
               array = 'OBSERVED SPECTRUM FILE PROBLEM!  I QUIT.'
               istat = ivwrite (line+4,3,array,40)
               stop
            endif
            call smooth (-1,ncall)
            choice = 'q'
         else
            write (*,1001)
            stop
         endif
         if (choice .eq. 'q') then
            call finish (0)
            exit
         endif
      enddo



!*****format statements
1001  format ('for syntheses, parameter "plot" must be 0, 1, 2, or 3;', &
              ' I QUIT!')
1002  format ('something wrong:  max number (99) of synthesis ', &
              'cycles exceeded; I QUIT!')



      end 





