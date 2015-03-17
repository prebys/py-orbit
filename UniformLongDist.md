# Summary #

This class generates longitudinal particle coordinate pairs that are uniformly distributed in z and dE. The class has a method which returns a single longitudinal coordinate pair (z and dE) each time it's called.

# Python Accessible Methods and Variables #
  1. **UniformLongDist(zmin, zmax, sp, eoffset, deltaEfrac)**. Creates the class.
    * zmin: Minimum z of the distribution.
    * zmax: Maximum z of the distribution.
    * sp: SyncParticle object (representing the synchronous particle of the distribution).
    * eoffset: Mean energy offset (set to zero for none)
    * deltaEfrac:  The fractional energy 1/2 spread of the beam, i.e., the maximum dE is: dE = eoffset + sp.getEnergy()(+- deltaEfrac)
  1. **getCoordinates()**. Routine that generates and returns a single coordinate pair within specified distribution.
# Examples #

In the following example a uniform longitudinal distribution with +/-1 MeV energy spread is created for a 1 GeV beam. The complete example can be found in $ORBIT\_ROOT/examples/Injection/ORBIT\_Benchmarks/JohoXYUniformL/injectuniformlong.py

```

import math
import sys
from bunch import Bunch

from injection import UniformLongDist

#------------------------------
#Bunch init
#------------------------------
b = Bunch()
runName = "Test_Injection"

b.mass(0.93827231)
b.macroSize(1.0e+1)
energy = 1.0 #Gev
b.getSyncParticle().kinEnergy(energy)

zlim = 120. * 248./360.
zmin = -zlim
zmax = zlim
deltaEfrac = 0.001
eoffset = 0.0

sp = b.getSyncParticle()

lFunc = UniformLongDist(zmin, zmax, sp, eoffset, deltaEfrac)

```