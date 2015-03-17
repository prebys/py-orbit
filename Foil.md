# Summary #

This is the class for modeling proton passage through rectangular carbon stripping foils of user-defined thickness. The class does not create particles. Two scattering methods are available to the particles, **traverseFoilFullScatter** and **traverseFoilSimpleScatter**. The first is a collimator style scattering model which includes multiple coulomb scattering, energy loss, and nuclear scattering. This is essentially an implementation of the [Collimator](Collimator.md) scattering class. It allows for particle loss and therefore requires a Bunch for handling the lost particles. The second is a simplified scattering model which only includes multiple coulomb scattering, and does not include a mechanism particle loss.

# Python Accessible Methods and Variables #
  1. **Foil(double xmin, double xmax, double ymin, double ymax, double thick)**. Creates a foil. Variables:
    * xmin: Horizontal minimum of foil.
    * xmax: Horizontal maximum of foil.
    * ymin: Vertical minimum of foil.
    * ymax: Vertical maximum of foil.
    * thick: Thickness of foil (ug/cm^2)
  1. **traverseFoilFullScatter(Bunch bunch, Bunch lostbunch)**. Routine to scatter particles using the full scattering model ([Collimator](Collimator.md) style). Scattering methods are located in the MaterialInteractions class.
    * bunch: The Bunch object to traverse the foil.
    * lostbunch: A Bunch object for populating with lost particles.
  1. **traverseFoilSimpleScatter(Bunch bunch)**. Routine to scatter particles using a simple multiple coulomb scattering implementation.
    * bunch: The Bunch object to traverse the foil.

# Example Scripts #

In the following example a standalone foil of thickness 400 ug/cm^2 is created and then used to scatter a bunch of particles with the full scattering model.

```

import math
import sys
from bunch import Bunch
from foil import Foil

print "Start."

xmin = -0.050
xmax = 0.050
ymin = -0.050
ymax = 0.050
# Below is 1000 times the width of normal foil but will do only one turn.
thick = 400

foil = Foil(xmin, xmax, ymin, ymax, thick)

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

#=====track bunch through Foil============

lostbunch = Bunch()
lostbunch.addPartAttr("LostParticleAttributes") 

foil.traverseFoilFullScatter(b, lostbunch)
b.dumpBunch("scatteredbunch.dat")
lostbunch.dumpBunch("lostbunch.dat")

print "Stop."

```

For the simple scattering model, one could exclude the creation of the lost bunch and replace the line,

```

foil.traverseFoilFullScatter(b, lostbunch)

```

with,
```

foil.traverseFoilSimpleScatter(b)

```


