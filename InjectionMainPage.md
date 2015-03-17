This page describes how to inject particles into a lattice according to user specified particle distributions.  A small library of particle distributions is available. The user must specify an appropriate tranverse injection region, which will often be the boundary of a foil. The injection methods do not provide foil scattering, which should be handed separately by use of the Foil class. The particles which land in the good region are added to a Bunch, and the particles that land outside of the good region are added to a lost particles Bunch.

# Associated Classes #

| **Name** | **Description** | **Address** |
|:---------|:----------------|:------------|
| InjectParts | The class which injects the particles | $ORBIT\_ROOT/py/orbit/injection/ |
| JohoTransverse | Class containing joho style transverse distributions | $ORBIT\_ROOT/py/orbit/injection// |
| JohoLongitudinal | Class containing joho style longitudinal distributions | $ORBIT\_ROOT/py/orbit/injection// |
| UniformLongDist | Class containing uniform longitudinal distributions | $ORBIT\_ROOT/py/orbit/injection// |
| UniformLongDistPaint | A variation of the UniformLongDist class that allows for time dependent zmin and zmax values | $ORBIT\_ROOT/py/orbit/injection// |
| [GULongDist](GULongDist.md) | Class containing longitudinal distributions that are Gaussian in dE and uniform in z| $ORBIT\_ROOT/py/orbit/injection// |
| [SNSESpreadDist](SNSESpreadDist.md) | Class containing SNS style longitudinal distributions | $ORBIT\_ROOT/py/orbit/injection// |
| [SNSESpreadDistPaint](SNSESpreadDistPaint.md) | A variation of the SNSESpreadDist class that allows for time dependent zmin and zmax values | $ORBIT\_ROOT/py/orbit/injection// |
| TeapotInjectionNode | The python module that allows a user to define Teapot InjectParts implementation | $ORBIT\_ROOT/py/orbit/injection// |
| injectionLatticeModifications | The python module that allows a user to add a defined TeapotInjectionNode to the teapot lattice in a drift region | $ORBIT\_ROOT/py/orbit/injection// |


# Examples #

The following example shows how to inject particles into a predefined Teapot lattice.  The script injects 10000 macro particles per pass through the lattice. The script can be found in $ORBIT\_ROOT/examples/Injection/lattice\_with\_injection\_node.py

```

##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and modifies this lattice by inserting 
# injection nodes
##############################################################

import math
import sys

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
from orbit.injection import TeapotInjectionNode
from orbit.injection import addTeapotInjectionNode
from injection import InjectParts
from injection import JohoTransverse, JohoLongitudinal, UniformLongDist

print "Start."

teapot_latt = teapot.TEAPOT_Lattice()
print "Read MAD."
teapot_latt.readMAD("MAD_Lattice/LATTICE","RING")
print "Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()," nodes=",len(teapot_latt.getNodes())

#====Injection aperature============
xmin = -0.050
xmax = 0.050
ymin = -0.050
ymax = 0.050

injectparams = (xmin, xmax, ymin, ymax)

#=====set up bunch stuff============
b = Bunch()

b.mass(0.93827231)
b.macroSize(1.0e+1)
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
tailfrac = 0.1  # 10% in tails
tailfac = 1.5   # tail emittance is 50% greater than core
zlim = 120. * 248./360.
zmin = -zlim
zmax = zlim
deltaEfrac = 0.001/2.
eoffset = 0.1
sp = b.getSyncParticle()

xFunc = JohoTransverse(order, alphax, betax, emitlim, xcenterpos, xcentermom, tailfrac, tailfac)
yFunc = JohoTransverse(order, alphay, betay, emitlim, ycenterpos, ycentermom, tailfrac, tailfac)
lFunc = UniformLongDist(zmin, zmax, sp, eoffset, deltaEfrac)

nparts = 10000.

injectnode = TeapotInjectionNode(nparts, b, lostbunch, injectparams, xFunc, yFunc, lFunc)

addTeapotInjectionNode(teapot_latt, 0., injectnode) 

print "===========Lattice modified ======================================="
print "New Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()," nodes=",len(teapot_latt.getNodes())

injectnode.track(paramsDict)

# dump ORBIT_MPI bunch to compare results
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch.dat")
bunch_pyorbit_to_orbit(teapot_latt.getLength(), lostbunch, "lostbunch.dat")
print "Stop."

```