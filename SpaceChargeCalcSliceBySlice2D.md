# Summary #
This class computes and applies transverse space charge kicks on 2D slices of a 3D grid.

# Description #
This class computes and applies transverse space charge kicks on 2D slices of a 3D grid. The method is to grid the particles on a rectangular 3 dimensional grid with user defined resolution. The PoissonSolverFFT2D is applied independently on each of the 2D transverse slices of a Grid3D to calculate the potential. Finally the space charge kicks are obtained from an interpolation between different 2D slices. A zero potential boundary can be taken into account optionally.

# Associated Classes #

| **Name** | **Description** | **Address** |
|:---------|:----------------|:------------|
| SpaceChargeCalcSliceBySlice2D | The class which calculates and applies the slice-by-slice 2D space charge kicks | $ORBIT\_ROOT/src/spacecharge/ |
| setSC2DSliceBySliceAccNodes | The python module that allows a user to add a space charge calculator as child nodes to elements in a lattice | $ORBIT\_ROOT/py/orbit/space\_charge/sc2dslicebyslice/scLatticeModifications.py|


setSC2DSliceBySliceAccNodes(lattice, sc\_path\_length\_min, space\_charge\_calculator, boundary = None)


# Python Accessible Methods and Variables #
  1. **SpaceChargeCalcSliceBySlice2D(xSize, ySize, zSize)**. Creates a slice-by-slice 2D space charge calculator. Variables:
    * xSize: Number of horizontal grid bins.
    * ySize: Number of vertical grid bins.
    * zSize: Number of longitudinal grid bins.
  1. **SpaceChargeCalcSliceBySlice2D.trackBunch(bunch, length)**. Computes and applies the space charge kicks to the particles in the bunch. Variables:
    * bunch: The bunch under consideration.
    * length: The length of the lattice section to be tracked.
  1. **SpaceChargeCalcSliceBySlice2D.getRhoGrid()**. Returns the 3D density grid of macro particles.
  1. **SpaceChargeCalcSliceBySlice2D.getPhiGrid()**. Returns the 3D grid of the potential.
  1. **scLatticeModifications.setSC2DSliceBySliceAccNodes(lattice, sc\_path\_length\_min, space\_charge\_calculator, boundary = None)**. Method to add space charge calculator as child nodes to the lattice elements. Returns the number of space charge nodes installed. Variables:
    * lattice: The lattice.
    * sc\_path\_length\_min: Minimum path length in node in order to do a kick.
    * space\_charge\_calculator: A predefined space charge calculator (an instance of the SpaceChargeCalcSliceBySlice2D class).
    * boundary (optional): An instance of a boundary class. If not set, the space charge kicks will be computed without boundary.

# Example Scripts #

The following example sc\_SliceBySlice2D\_boundary.py demonstrates the use of the Slice-by-slice 2D space charge module. The bunch is tracked for 2 turns and the phase advance for each particle is computed with the BunchTuneAnalysis module. All necessary files for running this example can be found in $ORBIT\_ROOT/examples/SpaceCharge/scSliceBySlice2D
```

import math
import sys
import random
import sys

from bunch import Bunch
from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
from spacecharge import Boundary2D
from orbit.space_charge.sc2dslicebyslice import scAccNodes, scLatticeModifications
from spacecharge import SpaceChargeCalcSliceBySlice2D
from orbit.diagnostics import TeapotTuneAnalysisNode
from orbit.diagnostics import addTeapotDiagnosticsNode

from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import GaussDist3D

print "Start."

#=====Make a Teapot style lattice======

lattice = teapot.TEAPOT_Lattice()
print "Read MAD."
lattice.readMAD("SNSring_pyOrbitBenchmark.LAT","RING")
print "Lattice=",lattice.getName()," length [m] =",lattice.getLength()," nodes=",len(lattice.getNodes())

#------------------------------
# Main Bunch init
#------------------------------

n_particles = 100000
twissX = TwissContainer(alpha = 0.0046902, beta = 10.207, emittance = 3.e-5)
twissY = TwissContainer(alpha = 0.056823, beta = 10.639, emittance = 3.e-5)
twissZ = TwissContainer(alpha = 0., beta = 100000., emittance = 0.008)
dist = GaussDist3D(twissX,twissY,twissZ)

b = Bunch()
for i in range(n_particles):
	(x,xp,y,yp,z,zp) = dist.getCoordinates()
	b.addParticle(x, xp, y, yp, z, zp)

total_macroSize=1.e+14
b.mass(0.93827231)
energy = 1.0 #Gev
#b.readBunch(distribution_file, n_particles)
print "Bunch Generated."
b.getSyncParticle().kinEnergy(energy)
nParticlesGlobal = b.getSizeGlobal()
b.macroSize(total_macroSize/nParticlesGlobal)

#-----------------------------------
# Add Tune Analysis node
#-----------------------------------

tunes = TeapotTuneAnalysisNode("tune_analysis")
tunes.assignTwiss(10.207, 0.0469, -0.05, 0.0061, 10.639, 0.056)
addTeapotDiagnosticsNode(lattice, 0, tunes)

#-----------------------------------
# Add Space Charge nodes
#-----------------------------------

nMacrosMin = 1
sc_path_length_min = 0.00000001
sizeX = 32   #number of grid points in horizontal direction
sizeY = 32   #number of grid points in vertical direction
sizeZ = 16   #number of longitudinal slices

# Make a boundary
bpoints = 128
bmodes = 10
xdim = 0.073
ydim = 0.073
boundary = Boundary2D(bpoints, bmodes, "Circle", xdim, ydim)

sc_calc = SpaceChargeCalcSliceBySlice2D(sizeX,sizeY,sizeZ)
scLatticeModifications.setSC2DSliceBySliceAccNodes(lattice, sc_path_length_min, sc_calc, boundary)

#-----------------------------------
# Tracking
#-----------------------------------

paramsDict = {}
paramsDict["bunch"]= b

n_turns = 2
for i in range(n_turns):
	lattice.trackBunch(b, paramsDict)

b.dumpBunch("bunch_SliceBySlice2D_boundary_turn" + str(n_turns) + ".dat")

print "Stop."

```