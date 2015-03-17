# Summary #

Method for creating Joho style distributions in the longitudinal plane.  This class has a method which returns a single longitudinal coordinate pair (z and dE) each time it's called.  These distributions are based on the ACCSIM code Joho distribution models.

# Python Accessible Methods and Variables #
  1. **JohoLongitudinal(order, zlim, dElim, nlongbunches, deltazbunch, deltaznotch, tailfraction, tailfactor)**. Creates and instance of a longitudinal Joho distribution element.
    * order: The shape parameter for the Joho distribution. Options are:  0 for hollow shell in phase space;  0.5 for flat profile in real space; 1 for uniform in phase space; 1.5 for elliptical in phase space; 2 for parabolic in phase space; 3 for truncated Gaussian in phase space; infinite for gaussian in phase space.
    * zlim:  Longitudinal coordinate limit for the distribution generator for each individual longitudinal bunch
    * dElim: dE coordinate limit for the distribution generator for each individual longitudinal bunch
    * nlongbunches: number of discrete longitudinal bunches
    * deltazbunch: bunch separation when nlongbunches > 0.
    * deltaznotch: width of an empty notch at the center of the bunch (0 = no notch)
    * tailfraction: fraction of beam to be located in an extended tail.  Range is 0 to 1, where 0 represents no beam in the tail and 1 represents all beam in the tail.
    * tailfactor: fraction of tail emittance to core emittance.
  1. **getCoordinates()**. Routine that generates and returns a single coordinate pair within specified distribution.

# Examples #

In the following example a longitudinal Joho distribution is created. The shape parameter is chosen to create a Gaussian distribution. In this example, one longitudinal bunch is generated.  The complete example including injection using this distribution can be found in $ORBIT\_ROOT/examples/Injection/ORBIT\_Benchmarks/JohoXYL/injectjoho.py .

```

import math
import sys

from injection import JohoLongitudinal
  
zlim = 120. * 248./360.
dElim = 0.001
nlongbunches = 1
deltazbunch = 0
deltaznotch = 0
ltailfrac = 0
ltailfac = 0

lFunc = JohoLongitudinal(order, zlim, dElim, nlongbunches, deltazbunch, deltaznotch, ltailfrac, ltailfac)

```

This example, 40 longitudinal bunches with 0.1 degree extent and 5 degrees apart are generated.  A complete example including injection using this distribution can be found in $ORBIT\_ROOT/examples/Injection/ORBIT\_Benchmarks/JohoXYL/injectjoho\_microbunches.py

```

import math
import sys

from injection import JohoLongitudinal
  
zlim = 120. * 248./360.
dElim = 0.001
zlim = .1 * 248./360.
dElim = 0.001
nlongbunches = 40
deltazbunch = 5*248./360
deltaznotch = 0.
ltailfrac = 0
ltailfac = 0

lFunc = JohoLongitudinal(order, zlim, dElim, nlongbunches, deltazbunch, deltaznotch, ltailfrac, ltailfac)

```