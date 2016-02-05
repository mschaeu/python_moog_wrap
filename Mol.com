!******************************************************************************
!     this common block carries the data to and from the molecular
!     equilibrium calculations; it is also used to hold isotopic
!     abundance information
!******************************************************************************

!     amol = names of the species
!     smallmolist = the small set of default molecule names
!     largemolist = the large set of default molecule names
!     xmol = number density of the species at all atmopshere layers
!     xatom = working array at a particular layer:  the number densities
!             of neutral atomic species
!     patom = working array at a particular layer:  the partial pressures
!             of neutral atomic species
!     pmol  = working array at a particular layer:  the partial pressures
!             of molecules
!     xamol = number densities of neutral atomic species at all layers
!     const = molecular constants loaded in from Bmolec.


      real*8       pmol(110), xmol(110,100), xamol(30,100), &
                  xatom(30), patom(30), &
                  amol(110), smallmollist(110), largemollist(110), &
                  datmol(7,110), const(6,110)
      integer      neq, lev, nmol, natoms, iorder(30), molopt, molset


      common /mol/ pmol, xmol, xamol, &
                   xatom, patom, &
                   amol, smallmollist, largemollist, &
                   datmol, const, &
                   neq, lev, nmol, natoms, iorder, molopt, molset

