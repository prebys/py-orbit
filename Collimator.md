# Summary #
This class constructs a collimator and calls the materials scattering routines to propagate the particles through the collimator. The stepsize of the particle in the collimator is based on the mean free path of the particle in the chosen material with a random factor. The scattering methods used are in the MaterialInteractions class. More information on the implementation of this class can be found in the CollimatorMainPage.


# Python Accessible Methods and Variables #
  1. **Collimator(int length, int ma, double density\_fac, int shape, double  a, double b, double c, double d, double angle, double pos=0, name = "my collimator")**. An instance of a collimator. Creates a collimator. Variables:
    * length: Length of the collimator in meters
    * ma: Material type index: 1=carbon, 2=aluminum, 3=iron, 4=copper, 5=tantalum, 6=tungstun, 7=platinum, 8=lead, 9 = black absorber.
    * density\_fac: A multiplier on the density of chosen material. Defaults to 1.
    * shape: Shape type index: 1=circle, 2=ellipse, 3=one sided flat, 4=two sided flat, 5=rectangular (outside is collimator), 6=rectangular (inside is collimator).
    * a: Radius for circular collimator; x-axis for elliptical; minimum horizontal extent for flat and rectangular.
    * b: Y-axis for elliptical; maximum horizontal extent for two-sided flat or rectangular.
    * c: Minimum vertical extent for rectangular.
    * d: Maximum vertical extent for rectangular.
    * angle: Tilt angle of the collimator defined in degrees.
    * pos: The optional argument of the starting z location of the collimator in the lattice.  If specified the lost particle longitudinal coordinates will be given as this pos + the relative lost position in the collimator.  If the argument is not used then the default is zero.
    * name: optional user specified name.  If the argument is left out of the list, then name is assigned to "collimator no name"
  1. **collimateBunch(Bunch bunch, Bunch lostbunch)**. The method to perform the collimation of the bunch.
    * bunch: The Bunch to be collimated.
    * lostbunch: The Bunch to be populated with particles absorbed in the collimator.
    1. **setPosition(double position)**. The method sets the position of the collimator in the lattice.
      * position: The start position of the collimator in the lattice.

# Example Scripts #

The following example demonstrates a bunch being propagated through a stand-alone circular iron collimator of length 0.5 meters and aperture 10 cm. After collimation the alive particles and the lost particles are printed. This example can also be found in: $ORBIT\_PATH/examples/Collimation as collimatebunch.py .

```

##############################################################
# This script sets up a collimator and collimates a bunch.
##############################################################

import math
import sys
from bunch import Bunch
from orbit.utils import orbitFinalize, NamedObject, ParamsDictObject
from collimator import Collimator
print "Start."

length = 0.5
ma = 3
density_fac = 1.0
shape = 1
a = 0.01
b = 0
c = 0
d = 0
angle = 0

collimator =Collimator(length, ma, density_fac, shape, a, b, c, d, angle)

#------------------------------
#Main Bunch init
#------------------------------
b = Bunch()
print "Read Bunch."
runName = "Benchmark_Collimator"

b.mass(0.93827231)
b.macroSize(1.0e+1)
energy = 1.0 #Gev
b.readBunch("parts.dat")
b.getSyncParticle().kinEnergy(energy)

#=====track bunch through Collimator Node============

lostbunch = Bunch()
lostbunch.addPartAttr("LostParticleAttributes") 

collimator.collimateBunch(b, lostbunch)
b.dumpBunch("collimatedbunch.dat")
lostbunch.dumpBunch("lostbunch.dat")
print "Stop."

}

```