# Summary #

This class injects particles into a pre-defined Bunch according to user-specified distribution functions. The user must specify the allowed transverse physical region for injection.  Particles outside of this region will be added to a lost particle Bunch. A library of distribution functions is available.

# Python Accessible Methods and Variables #
  1. **InjectParts(nparts, bunch, lostbunch, injectregion, xDistFunc, yDistFunc, lDistFunc)**. Creates an instance of an injection element.
    * nparts: Number of macro particles to be injected during each pass through the element.
    * bunch: The bunch to be injected into.
    * lostbunch: The lost particle bunch which will be populated by injected particles outside of the injection region.
    * injectregion: A list of ransverse limits of the injection region (xmin: Horizontal minimum; xmax: Horizontal maximum; ymin: Vertical minimum; ymax: Vertical maximum).
    * xDistFunc: Horizontal distribution function. One example is JohoTransverse.
    * yDistFunc: Vertical distribution function: One example is JohoTransverse.
    * lDistFunc: Longitudinal distribution function.  Examples are JohoLongitudinal and UniformLongDist.
  1. **addParticles()**. Method that gets the particle coordinates from the distributions, checks to see if they are within the injection region, and adds them to either the main bunch or the lost bunch.

# Example Scripts #

In the following example a standalone injection node with Joho style distributions is created and used to inject 10,000 macro particles. The script can be found in $ORBIT\_ROOT/examples/Injection/ORBIT\_Benchmarks/JohoXYL.  Other examples can be found in $ORBIT\_ROOT/examples/Injection/ORBIT\_Benchmarks/

```

import math
import sys
from bunch import Bunch
from foil import Foil
from collimator import Collimator
from injection import InjectParts
from injection import JohoTransverse, JohoLongitudinal
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
print "Start."

xmin = -0.050
xmax = 0.050
ymin = -0.050
ymax = 0.050

foilparams = (xmin, xmax, ymin, ymax)

#------------------------------
#Bunch init
#------------------------------
b = Bunch()
runName = "Test_Injection"

b.mass(0.93827231)
b.macroSize(1.0e+1)
energy = 1.0 #Gev
b.getSyncParticle().kinEnergy(energy)

lostfoilbunch = Bunch()
lostfoilbunch.addPartAttr("LostParticleAttributes") 

#------------------------------
#Initial Distribution Functions
#------------------------------

order = 3.
alphax = 0.063
betax = 10.209
alphay = 0.063
betay = 10.776
emitlim = 0.00152 * 2*(order + 1) * 1e-6
xcenterpos = 0.0468
xcentermom = 0.001
ycenterpos = 0.0492
ycentermom = -0.00006
tailfrac = 0
tailfac = 1
zlim = 120. * 248./360.
dElim = 0.001
nlongbunches = 100
deltazbunch = 0.688
deltaznotch = 0
ltailfrac = 0
ltailfac = 0

xFunc = JohoTransverse(order, alphax, betax, emitlim, xcenterpos, xcentermom, tailfrac, tailfac)
yFunc = JohoTransverse(order, alphay, betay, emitlim, ycenterpos, ycentermom, tailfrac, tailfac)
lFunc = JohoLongitudinal(order, zlim, dElim, nlongbunches, deltazbunch, deltaznotch, ltailfrac, ltailfac)

#------------------------------
# Inject some particles
#------------------------------

nparts = 10000
inject = InjectParts(nparts, b, lostfoilbunch, foilparams, xFunc, yFunc, lFunc)
inject.addParticles()

bunch_pyorbit_to_orbit(248.0, b, "pybunch.dat")
lostfoilbunch.dumpBunch("pylostbunch.dat")
print "Stop."
quit()

```