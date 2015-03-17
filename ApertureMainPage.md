# Summary #

The aperture class provides a method for losing particles on a zero length, 100% absorbing aperture.  Three types of aperture shapes are available: 1) Circle, ellipse, and rectangle.  The user can apply horizontal and vertical offsets to the aperture. The apertures can be added singly to any drift region in a lattice, or in sets within a positional range in the lattice.  If the set is used, one aperture will be added per drift node in the accelerator.

The module requires input of two bunches from the Bunch class: one bunch for the alive particles, and one bunch to be populated with particles lost in the aperture.

# Associated Classes #

| **Name** | **Description** | **Address** |
|:---------|:----------------|:------------|
| [Aperture](Aperture.md) | The class which defines the aperture | $ORBIT\_ROOT/src/orbit/Aperture/ |
| TeapotApertureNode | The python class that allows a user to define an TeaPot style aperture  | $ORBIT\_ROOT/py/orbit/aperture/TeapotApertureNode.py |
| CircleApertureNode | Convenience class for quickly creating a  circular TeaPot aperture | $ORBIT\_ROOT/py/orbit/aperture/TeapotApertureNode.py |
| EllipseApertureNode | Convenience class for quickly creating an elliptical TeaPot aperture | $ORBIT\_ROOT/py/orbit/aperture/TeapotApertureNode.py|
| RectangleApertureNode | Convenience class for quickly creating a rectangular TeaPot aperture | $ORBIT\_ROOT/py/orbit/aperture/TeapotApertureNode.py |
| addTeapotApertureNode | The python routine that allows a user to add a defined aperture to the lattice in a drift region | $ORBIT\_ROOT/py/orbit/aperture/ApertureLatticeModifications.py|
| addCircleApertureSet | Python routine that allows a user to add a set of circular apertures to the lattice over a range of distance (one per node) | $ORBIT\_ROOT/py/orbit/aperture/ApertureRangeModifications.py|
| addEllipseApertureSet | Python routine that allows a user to add a set of ellipticalapertures to the lattice over a range of distance (one per node) | $ORBIT\_ROOT/py/orbit/aperture/ApertureRangeModifications.py|
| addRectangleApertureSet | Python routine that allows a user to add a set of rectangular apertures to the lattice over a range of distance (one per node) | $ORBIT\_ROOT/py/orbit/aperture/ApertureRangeModifications.py|

# Python Accessible Classes, Methods, and Variables #
  1. **TeapotApertureNode(int shape, double a, double b, double c, double d, name = "aperture")**.  Creates a teapot style aperture. Variables:
    * shape: Shape type index. 1=circle, 2=ellipse, 3=rectangle.
    * a: Radius for circular collimator, x-axis half length for ellipse and rectangle.
    * b: Y-axis half length for elliptical and rectangle. Doesn't apply for circle.
    * c: Horizontal offset. Default 0. Argument is optional.
    * d: Vertical offset. Default 0. Argument is optional.
    * name: The name of the aperture.
  1. **CircleApertureNode(double a, double c, double d, name = "aperture")**.  Creates a teapot style aperture. Variables:
    * a: Radius.
    * c: Horizontal offset.  Default 0. Argument is optional.
    * d: Vertical offset. Default 0. Argument is optional.
    * name: The name of the aperture.
  1. **EllipseApertureNode(double a, double b, double c, double d, name = "aperture")**.  Creates a teapot style aperture. Variables:
    * a: X-axis half length.
    * b: Y-axis half length.
    * c: Horizontal offset.  Default 0. Argument is optional.
    * d: Vertical offset. Default 0. Argument is optional.
    * name: The name of the aperture.
  1. **RectangleApertureNode(double a, double b, double c, double d, name = "aperture")**.  Creates a teapot style aperture. Variables:
    * a: X-axis half length.
    * b: Y-axis half length.
    * c: Horizontal offset.  Default 0. Argument is optional.
    * d: Vertical offset. Default 0. Argument is optional.
    * name: The name of the aperture.
  1. **addTeapotApertureNode(TEAPOT\_Lattice lattice, double position, TeapotApertureNode aperture)**. A method to add the aperture to a Teapot lattice. Variables:
    * lattice: The Teapot lattice the aperture should be added to.
    * position: The start position of the aperture in the lattice.
    * aperture: A previously defined aperture
  1. **addCircleApertureSet(double a, TEAPOT\_Lattice lattice, double s, double e, double c, double d)**. Adds a set of circular apertures in the range from variable s to e.  Variables:
    * a: radius
    * lattice: The Teapot lattice the aperture should be added to.
    * s: Start position.
    * e: End position.
    * c: Horizontal offset.  Default 0. Argument is optional.
    * d: Vertical offset. Default 0. Argument is optional.
  1. **addEllipseApertureSet(double a, double b, TEAPOT\_Lattice lattice, double s, double e, double c, double d)**. Adds a set of elliptical apertures in the range from variable s to e. Variables:
    * a: X-axis half length
    * b: Y-axis half length
    * lattice: The Teapot lattice the aperture should be added to.
    * s: Start position.
    * e: End position.
    * c: Horizontal offset.  Default 0. Argument is optional.
    * d: Vertical offset. Default 0. Argument is optional.
  1. **addRectangleApertureSet(double a, double b, TEAPOT\_Lattice lattice, double s, double e, double c, double d)**. Adds a set of rectangular apertures in the range from variable s to e. Variables:
    * a: X-axis half length
    * b: Y-axis half length
    * lattice: The Teapot lattice the aperture should be added to.
    * s: Start position.
    * e: End position.
    * c: Horizontal offset.  Default 0. Argument is optional.
    * d: Vertical offset. Default 0. Argument is optional.
# Example Scripts #

The following example demonstrates the creation of a circular aperture and it's addition to the lattice.  Also shown in comment fields are other aperture options.
```

##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and modifies this lattice by inserting 
# diagnotics nodes
##############################################################

import math
import sys

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit

from orbit.diagnostics import StatLats
from orbit.diagnostics import addTeapotDiagnosticsNode
from orbit.diagnostics import TeapotStatLatsNode
from orbit.diagnostics import addTeapotDiagnosticsNodeSet
from orbit.utils import orbitFinalize, NamedObject, ParamsDictObject
from orbit.aperture import addTeapotApertureNode
from orbit.aperture import TeapotApertureNode, CircleApertureNode, EllipseApertureNode, RectangleApertureNode
from orbit.aperture import addCircleApertureSet, addEllipseApertureSet, addRectangleApertureSet
print "Start."

teapot_latt = teapot.TEAPOT_Lattice()
print "Read MAD."
teapot_latt.readMAD("MAD_Lattice/LATTICE","RING")
print "Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()," nodes=",len(teapot_latt.getNodes())

Aperturenode = CircleApertureNode(.01)
#Aperturenode = CircleApertureNode(.01, .002, .0025)
#Aperturenode = RectangleApertureNode(.02, .015, -.002, -.0025)
addTeapotApertureNode(teapot_latt, 240, Aperturenode)
#addEllipseApertureSet(.02, .01,  teapot_latt, 150, 195)

print "===========Lattice modified ======================================="
print "New Lattice=",teapot_latt.getName()," length [m] =",teapot_latt.getLength()," nodes=",len(teapot_latt.getNodes())

print "============= nodes inside the region ==========="
# print all nodes around the specified position
for node in teapot_latt.getNodes():
	print "node=",node.getName()," type=",node.getType()," L=",node.getLength()

#------------------------------
#Main Bunch init
#------------------------------

b = Bunch()
print "Read Bunch."
runName = "Benchmark_Diagnostics"
b.mass(0.93827231)
b.macroSize(1.0e+1)
energy = 1.0 #Gev
bunch_orbit_to_pyorbit(teapot_latt.getLength(), energy, "Bm_KV_Uniform_1000",b)
b.getSyncParticle().kinEnergy(energy)
paramsDict = {}
lostbunch = Bunch()
paramsDict["lostbunch"]=lostbunch
paramsDict["bunch"]= b
lostbunch.addPartAttr("LostParticleAttributes")

#=====track bunch ============

print "Tracking..."
for turn in range(1):
	teapot_latt.trackBunch(b, paramsDict)
b.dumpBunch("bunch_final.dat")
lostbunch.dumpBunch("lostbunch_final.dat")

print "Stop."





```

Complete run examples can be found in $ORBIT\_ROOT/examples/Apertures .