This class creates longitudinal particle coordinates that are uniform in z and Gaussian in energy.

# Python Accessible Methods and Variables #
  1. **GULongDis(zmin, zmax, sp, emean, esigma, etrunc, emin, emax)**. Creates the class.
    * zmin: Minimum z of the distribution.
    * zmax: Maximum z of the distribution.
    * sp: The SyncPart object for the distribution.
    * emean: Mean energy of the distribution.
    * esigma: Sigma spread of the energy distribution
    * etrunc: Flag for truncating the distribution.  0 is no truncation, 1 is truncation as specified by user parameters emit and emax.
    * emin: If etrunc is not 0, then this is the minimum of the energy distribution [GeV](GeV.md).
    * emax: If etrunc is not 0, then this is the maximum of the energy distribution [GeV](GeV.md).
  1. **getCoordinates()**. Routine that generates and returns a single coordinate pair within specified distribution.


# Examples #

In the following example a longitudinal distribution generator is created for a distribution that is uniform in z and Gaussian in energy.  The energy distribution is truncated at 1.003 GeV on the high end and 0.997 GeV on the low end. The complete example can be found in $ORBIT\_ROOT/examples/Injection/ORBIT\_Benchmarks/JohoXYGaussEUniformL/injectgulong.py

```

import math
import sys

from injection import JohoTransverse, GULongDist

zlim = 120. * 248./360.
zmin = -zlim
zmax = zlim
emean = 1.0
esigma = 0.001
etrunc = 1
emin = 0.997
emax = 1.003

sp = b.getSyncParticle()
lFunc = GULongDist(zmin, zmax, sp, emean, esigma, etrunc, emin, emax)

```