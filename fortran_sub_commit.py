#!/usr/bin/python

import sys
import os
import re
import time
import numpy
import math
import subprocess
import scipy.stats
import matplotlib.pyplot as plt
from ctypes import *

def main():
    #now we add the 'moog' library to our list of Fortran libraries
    add_moog = cdll.LoadLibrary("/Applications/moog_library_complete/moog.so")
    moog = add_moog.moogsilent_

    #now, we declare some arrays that we will use for our interaction with MOOG
    abfind = ((c_double*2500)*6)()
    #this abfind array holds:
    #   column 0: wavelengths (in Angstroms)
    #   column 1: atomic identifier
    #   column 2: excitation potential (in eV)
    #   column 3: log(gf)
    #   column 4: measured equivalent width (in Angstroms!)
    #   column 5: abundance (in log (epsilon))
    
    synth = ((c_float*500000)*8)()
    #this synth array holds:
    #   column 0: wavelength of observed spectrum in Angstroms
    #   column 1: flux of observed spectrum
    #   column 2: wavelength of synthesis in Angstroms
    #   column 3: flux of synthesis 1
    #   column 4: flux of synthesis 2
    #   column 5: flux of synthesis 3
    #   column 6: flux of synthesis 4
    #   column 7: flux of synthesis 5
    
    #as a final step, we also declare which MOOG driver we use. The numbers
    #are defined as follows:
    #   0: abfind
    #   1: synplot
    #   2: synth
    #   3: cogsyn
    #   4: blends
    #   5: ewfind
    #   6: cog
    #   7: calmod
    #   8: doflux
    #   9: weedout
    #  10: gridsyn
    #  11: gridplot
    #  12: binary
    #  13: adpop
    #  14: synpop
    driver_version = c_int(0)
    
    #and finally call moog
    moog(byref(driver_version), abfind, synth)

    #to access these values, use the following syntax:
    #abfind[0][0] gives the wavelength


if __name__ == '__main__':
    main()
