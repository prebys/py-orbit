# Summary #

This class is identical to the [SNSESpreadDist](SNSESpreadDist.md) class with added functionality of time dependent zmin and zmax variables. Zmin and zmax must both be lists of no less than two pairs and functions (in the mathematical sense) of time. When getCoordinates() is called, the time of the synchronous particle is used and the given functions are interpolated for values.

The number of macro particles injected per turn is a function of the zmin and zmax values. The initial injection rate and value of zmax-zmin determine the macropaticle density that is to be used for the whole run. The number of macro particles injected per turn then becomes: **macropaticleDensity**x **(zmax-zmin)**.

Warning: If you decide to use the nmaxmacropartcles functionallity of the InjectParts class, you must take into account the paragraph above. It becomes rather burdensome.


# Python Accessible Methods and Variables #

  1. **SNSESpreadDistPaint(lattlength, zminFunc, zmaxFunc, tailfraction, sp, emean, esigma, etrunc, emin, emax, ecparams, esparams)**. Creates the class.
    * lattlength: Total lattice length.
    * zminFunc: Function of time (list of pairs) for the minimum z of the distribution.
    * zmaxFunc: Function of time (list of pairs) for the maximum z of the distribution.
    * lattlength: Total lattice length.
    * sp: [SyncParticle](SyncParticle.md) object (representing the synchronous particle of the distribution).
    * emean: Mean energy of the distribution.
    * esigma: Sigma spread of the energy distribution
    * etrunc: Flag for truncating the distribution. 0 is no truncation, 1 is truncation as specified by user parameters emit and emax.
    * emin: If etrunc is not 0, then this is the minimum of the energy distribution.
    * emax: If etrunc is not 0, then this is the maximum of the energy distribution.
    * ecparams: Parameters of the energy centroid jitter. List is ecparams = (ecmean, ecsigma, ectrunc, ecmin, ecmax, ecdrifti, ecdriftf, drifttime), where ecmean = mean of the centroid jitter; ecsigma = sigma of the centroid jitter; ectrunc = flag for truncation (0 for none, 1 for truncate); ecmin = the minimum energy for the jittered distribution, when ectrunc is not zero; ecmax = the minimum energy for the jittered distribution, when ectrunc is not zero; ecdrifti = initial drift energy of the centroid; ecdriftf = final drift energy of the centroid; driftime = total duration of centroid drift.
    * esparams: Parameters of the energy sinusoidal spreader. List is esparams = (esnu, esphase, esmax, nulltime) , where esnu = frequency of the sinusoidal spreader [Hz](Hz.md); esphase = constant phase of the spreader, esmax = amplitude of the spreader; nulltime = time for decay of spreader voltage
  1. **getCoordinates()**. Routine that generates and returns a single coordinate pair within specified distribution.

# Example #
The following is a snippet from an example that can be found in  $ORBIT\_ROOT/examples/Injection/
The example simulates the SNS ring with the SNSESpreadDistPaint distribution class implemented. Each call of 'trackBunch' in the full example updates the synchronous particle time and injects particles accordingly.


```
import math
import sys

from injection import SNSESpreadDistPaint

print "Start."

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

#=====Make a Teapot style lattice======

teapot_latt = teapot.TEAPOT_Ring()
print "Read MAD."
teapot_latt.readMAD("SNSring_pyOrbitBenchmark.LAT","RING")
print "Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()," nodes=",len(teapot_latt.getNodes())


#------------------------------
#Initial Distribution Functions
#------------------------------

lattlength = teapot_latt.getLength()
sp = b.getSyncParticle()

order = 3.
alphax = 0.063
betax = 10.209
alphay = 0.063
betay = 10.776
emitlim = 0.152 * 2*(order + 1) * 1e-6
xcenterpos = 0.0468
xcentermom = 0.00
ycenterpos = 0.0492
ycentermom = 0.00

zlim = 120. * lattlength/360.

T = 9.45424317281e-07 #peroid of the synchronous particle

#This function has zmax start at half the zlimit and linearly slope to
#its full value over 50 turns. Then it remains constant. Zmin does the same
# only it is shifted down by 2*zlim.
zmaxFunc = [[0,zlim/2],[T*49,zlim],[T*50,zlim]]
zminFunc = [[0,-(3/2)*zlim], [T*49,-zlim], [T*50,-zlim]]

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

lFunc =SNSESpreadDistPaint(lattlength, zminFunc, zmaxFunc, tailfraction, sp, emean, esigma, etrunc, emin, emax, ecparams, esparams)


```