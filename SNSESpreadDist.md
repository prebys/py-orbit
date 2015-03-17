# Summary #

This class creates SNS style longitudinal particle coordinates. The distribution is uniform in z and Gaussian in energy, with user option for adding centroid energy jitter and sinusoidal energy variation. The class has a method which returns a single longitudinal coordinate pair (z and dE) each time it's called.

# Python Accessible Methods and Variables #
  1. **SNSESpreadDist(lattlength, zmin, zmax, tailfraction, sp, emean, esigma, etrunc, emin, emax, ecparams, esparams)**. Creates the class.
    * lattlength: Total lattice length.
    * zmin: Minimum z of the distribution.
    * zmax: Maximum z of the distribution.
    * tailfraction: Fractional amount of beam to be placed in an extended tail.
    * sp: SyncParticle object (representing the synchronous particle of the distribution).
    * emean: Mean energy of the distribution.
    * esigma: Sigma spread of the energy distribution
    * etrunc: Flag for truncating the distribution.  0 is no truncation, 1 is truncation as specified by user parameters emit and emax.
    * emin: If etrunc is not 0, then this is the minimum of the energy distribution.
    * emax: If etrunc is not 0, then this is the maximum of the energy distribution.
    * ecparams: Parameters of the energy centroid jitter.  List is ecparams = (ecmean, ecsigma, ectrunc, ecmin, ecmax, ecdrifti, ecdriftf, drifttime), where ecmean = mean of the centroid jitter; ecsigma = sigma of the centroid jitter; ectrunc = flag for truncation (0 for none, 1 for truncate); ecmin = the minimum energy for the jittered distribution, when ectrunc is not zero; ecmax = the minimum energy for the jittered distribution, when ectrunc is not zero; ecdrifti = initial drift energy of the centroid; ecdriftf = final drift energy of the centroid; driftime = total duration of centroid drift.
    * esparams: Parameters of the energy sinusoidal spreader.  List is esparams = (esnu, esphase, esmax, nulltime) , where esnu = frequency of the sinusoidal spreader [Hz](Hz.md); esphase = constant phase of the spreader, esmax = amplitude of the spreader; nulltime = time for decay of spreader voltage
  1. **getCoordinates()**. Routine that generates and returns a single coordinate pair within specified distribution.

# Examples #

In the following example a 1 GeV SNS style distribution with centroid energy jitter is created. There is no sinusoidal energy variation in this example.  The complete example can be found in $ORBIT\_ROOT/examples/Injection/ORBIT\_Benchmarks/JohoXYSNSESpreadL/injectSNSESpread.py


```
import math
import sys
from bunch import Bunch

from injection import JohoTransverse, SNSESpreadDist

#------------------------------
#Bunch init
#------------------------------
b = Bunch()
runName = "Test_Injection"

b.mass(0.93827231)
b.macroSize(1.0e+1)
energy = 1.0 #Gev
b.getSyncParticle().kinEnergy(energy)
sp = b.getSyncParticle()

lattlength = 248.00935;
zlim = 120. * lattlength/360.
zmin = -zlim
zmax = zlim
tailfraction = 0
emean = sp.kinEnergy()
efac = 0.784
esigma = 0.0015*efac
etrunc = 1.
emin = sp.kinEnergy() - 0.0025*efac
emax = sp.kinEnergy() + 0.0025*efac

ecmean = 0
ecsigma = 0.0015*efac
ectrunc = 1.
ecmin = -0.0035*efac
ecmax = 0.0035*efac
ecdrifti = 0
ecdriftf = 0
turns = 1000.
tturn = lattlength / (sp.beta() * 2.998e8)
drifttime= 1000.*turns*tturn

ecparams = (ecmean, ecsigma, ectrunc, ecmin, ecmax, ecdrifti, ecdriftf, drifttime)

esnu = 100.
esphase = 0.
esmax = 0
nulltime = 0
esparams = (esnu, esphase, esmax, nulltime) 

lFunc = SNSESpreadDist(lattlength, zmin, zmax, tailfraction, sp, emean, esigma, etrunc, emin, emax, ecparams, esparams)


```