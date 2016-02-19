
      subroutine moogsilent(driver_version, abfind_ret, synth_ret, atm_external, &
                            aux_data, atm_data)
!******************************************************************************
!     This is the main driver for the non-interactive version of MOOG.  
!     It reads the parameter file and sends MOOG to various controlling 
!     subroutines.  In this version of MOOG the parameter file must
!     be named "batch.par" (because the code cannot stop to ask the
!     user to name the parameter file)
!******************************************************************************

      include 'Atmos.com'
      include 'Pstuff.com'
      include 'Linex.com'

! this was added by MS on 02/04/16 to make moog return the correct variables
      integer*4 atm_external
      real*8 aux_data(5)
      real*8 atm_data(100,4)
      integer*4 driver_version
      real*8 abfind_ret(2500,6)
      real*4 synth_ret(500000,8)


!$$$$$$$$$$$$$$$$$$$$$$$$ USER SETUP AREA $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     in compiling MOOG, here the various machine-specific things are
!     declared.  First, define the directory where MOOG lives, in order to
!     be able to pull in auxiliary data files; executing 'make' will
!     generate a reminder of this
      moogpath = '/Applications/moog_library_complete/'

!*****What kind of machine are you using?  Possible ones are:
!     "mac" = Intel-based Apple Mac
!     "pcl" = a PC or desktop running some standard linux like Redhat
!     "uni" = a machine running Unix, specifically Sun Solaris
      machine = "mac"

!$$$$$$$$$$$$$$$$$$$$$$$ END OF USER SETUP $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!*****declare this to be the non-interactive version; variable "silent"
!     will be queried on all occasions that might call for user input;
!     DON'T CHANGE THIS VARIABLE; 
!     if silent = 'n', the normal interactive MOOG is run;
!     if silent = 'y', the non-interactive MOOG is run
      silent = 'y'


!*****invoke the overall starting routine
      control = '       '
      call begin

!*****use one of the standard driver routines ("isotop" is obsolete):
      if     (driver_version .eq. 0) then
         control = 'abfind'
         call abfind(atm_external, aux_data, atm_data)
!        now we equate the corresponding arrays
         abfind_ret(:,1) = wave1
         abfind_ret(:,2) = atom1
         abfind_ret(:,3) = e(:,1)
         abfind_ret(:,4) = dlog10(gf)
         abfind_ret(:,5) = width
         abfind_ret(:,6) = abundout
      elseif (driver_version .eq. 1) then
         control = 'synplot'
         call plotit
      elseif (driver_version .eq. 2) then
         control = 'synth  '
         call synth(atm_external, aux_data, atm_data)
         synth_ret(:,1) = xobs
         synth_ret(:,2) = yobs
         synth_ret(:,3) = xsyn
         synth_ret(:,4) = chunk(:,1)
         synth_ret(:,5) = chunk(:,2)
         synth_ret(:,6) = chunk(:,3)
         synth_ret(:,7) = chunk(:,4)
         synth_ret(:,8) = chunk(:,5)
      elseif (driver_version .eq. 3) then
         control = 'cogsyn '
         call cogsyn
      elseif (driver_version .eq. 4) then
         control = 'blends '
         call blends
      elseif (driver_version .eq. 5) then
         control = 'ewfind '
         call ewfind
      elseif (driver_version .eq. 6) then
         control = 'cog    '
         call cog
      elseif (driver_version.eq. 7) then
         control = 'calmod '
         call calmod
      elseif (driver_version .eq. 8) then
         control = 'doflux '
         call doflux
      elseif (driver_version .eq. 9) then
         control = 'weedout'
         call weedout
      elseif (driver_version .eq. 10) then
         control = 'gridsyn'
         call gridsyn
      elseif (driver_version .eq. 11) then
         control = 'gridplo'
         call gridplo
      elseif (driver_version .eq. 12) then
         control = 'binary '
         call binary
      elseif (driver_version .eq. 13) then
         control = 'abpop  '
         call abpop
      elseif (driver_version .eq. 14) then
         control = 'synpop '
         call synpop
!*****or else you are out of luck!
      else
         write(*,*) "Wrong identifier entered. Exiting."
         return
      endif

      end
