# Summary #

Method for creating Joho style distributions in the transverse plane.  This class has a method which returns a single transverse coordinate pair (position and momentum) each time it's called.  These distributions are based on the ACCSIM code Joho distribution models.

# Python Accessible Methods and Variables #
  1. **JohoTransverse(order, alpha, beta, emitlim, centerpos, centermom, tailfraction, tailfactor)**. Creates and instance of a transverse Joho distribution element.
    * order: The shape parameter for the Joho distribution. Options are:  0 for hollow shell in phase space;  0.5 for flat profile in real space; 1 for uniform in phase space; 1.5 for elliptical in phase space; 2 for parabolic in phase space; 3 for truncated Gaussian in phase space; infinite for gaussian in phase space.
    * alpha: Twiss alpha of the injected particle
    * beta: Twiss beta of the injected particle
    * emitlim: Twiss emittance limit for the injected particle
    * centerpos: The mean position coordinate of the distribution
    * centermom: The mean momentum coordinate of the distribution
    * tailfraction: Fraction of particles to be located in an extended tail (Range is 0 to 1. 0 is none, 1 is all).
    * tailfactor: The ratio of the tail emittance to core emittance. For instance the value 1 means the core and the tail have the same emittance.
  1. **getCoordinates()**. Routine that generates and returns a single coordinate pair within specified distribution.

# Example Scripts #

In the following example a Joho distributions are created for both the X and Y plane. The shape parameter is chosen to create a Gaussian distribution. 10% of the beam is located in an beam tail with emittance 50% larger than the core emittance. A complete example including injection using this distribution can be found in $ORBIT\_ROOT/examples/Injection/ORBIT\_Benchmarks/JohoXYL .

```

import math
import sys
from injection import JohoTransverse

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
tailfrac = 0.1
taillim = 1.5

xFunc = JohoTransverse(order, alphax, betax, emitlim, xcenterpos, xcentermom, tailfrac, taillim)
yFunc = JohoTransverse(order, alphay, betay, emitlim, ycenterpos, ycentermom, tailfrac, taillim)      

```