# Summary #

This class is identical to the [UniformLongDist](UniformLongDist.md) class with added functionality of time dependent zmin and zmax variables. Zmin and zmax must both be lists of no less than two pairs and functions (in the mathematical sense) of time. When getCoordinates() is called, the time of the synchronous particle is used and the given functions are interpolated for values.

The number of macro particles injected per turn is a function of the zmin and zmax values. The initial injection rate and value of zmax-zmin determine the macropaticle density that is to be used for the whole run. The number of macro particles injected per turn then becomes: **macropaticleDensity**x **(zmax-zmin)**.

Warning: If you decide to use the nmaxmacropartcles functionallity of the InjectParts class, you must take into account the paragraph above. It becomes rather burdensome.


# Python Accessible Methods and Variables #
  1. **UniformLongDistPaint(zminFunc, zmaxFunc, sp, eoffset, deltaEfrac)**. Creates the class.
    * zminFunc: Function of time (list of pairs) for the minimum z of the distribution.
    * zmaxFunc: Function of time (list of pairs) for the maximum z of the distribution.
    * sp: SyncParticle object (representing the synchronous particle of the distribution).
    * eoffset: Mean energy offset (set to zero for none)
    * deltaEfrac:  The fractional energy 1/2 spread of the beam, i.e., the maximum dE is: dE = eoffset + sp.getEnergy()(+- deltaEfrac)
  1. **getCoordinates()**. Routine that generates and returns a single coordinate pair within specified distribution.

# Example #
The following is a snippet from an example that can be found in  $ORBIT\_ROOT/examples/Injection/
The example simulates the SNS ring with the UniformLongDistPaint distribution class implemented. Each call of 'trackBunch' in the full example updates the synchronous particle time and injects particles accordingly.
```

import math
import sys
from injection import UniformLongDistPaint

#=====Main bunch parameters============
intensity = 7.8e13
turns = 1000.0
macrosperturn = 260 #this is how many particles are injected per turn
macrosize = intensity/turns/macrosperturn

b = Bunch()
b.mass(0.93827231)
b.macroSize(macrosize)
energy = 1.0 #Gev
b.getSyncParticle().kinEnergy(energy)

paramsDict = {}
lostbunch = Bunch()
paramsDict["lostbunch"]=lostbunch
paramsDict["bunch"]= b
lostbunch.addPartAttr("LostParticleAttributes")

#------------------------------
#Initial Distribution Functions
#------------------------------

sp = b.getSyncParticle()

T = 9.45424317281e-07 #peroid of the synchronous particle

#This function has zmax start at half the zlimit and linearly slope to
#its full value over 50 turns. Then it remains constant. Zmin does the same
# only it is shifted down by 2*zlim.
zmaxFunc = [[0,zlim/2],[T*49,zlim],[T*50,zlim]]
zminFunc = [[0,-(3/2)*zlim], [T*49,-zlim], [T*50,-zlim]]


lFunc = UniformLongDistPaint(zminFunc,zmaxFunc,sp,0, .001)


```