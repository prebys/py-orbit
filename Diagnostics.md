# Summary #
The diagnostics package includes routines for finding the statistical Twiss parameters of a bunch, the moments, and the tunes and actions.

# Associated Classes #

| **Name** | **Description** | **Address** |
|:---------|:----------------|:------------|
| BunchTwissAnalysis | The c++ class which calculates the statistical twiss parameters and moments | $ORBIT\_ROOT/src/orbit/BunchDiagnostics/ |
| BunchTuneAnalysis | The c++ class  which calculates particle tunes and actions | $ORBIT\_ROOT/src/orbit/BunchDiagnostics/ |
| StatLats | The python class which calculates and prints the statistical twiss parameters. | $ORBIT\_ROOT/py/orbit/diagnostics/diagnostics.py|
| StatLatsSetMember | The same as StatLats but for use in a  set | $ORBIT\_ROOT/py/orbit/diagnostics/diagnostics.py|
| Moments | The python class which calculates and prints the beam moments. | $ORBIT\_ROOT/py/orbit/diagnostics/diagnostics.py|
| MomentsSetMember | The same as Moments but for use in a set | $ORBIT\_ROOT/py/orbit/diagnostics/diagnostics.py|
| TeapotStatLatsNode | The teapot implementation of the StatLats class | $ORBIT\_ROOT/py/orbit/diagnostics/TeapotDiagnosticsNode.py|
| TeapotStatLatsNodeSetMember | The same as TeaPotStatLatsNode but for use in a set | $ORBIT\_ROOT/py/orbit/diagnostics/TeapotDiagnosticsNode.py|
| TeapotMomentsNode | The Teapot implementation of the moments class | $ORBIT\_ROOT/py/orbit/diagnostics/TeapotDiagnosticsNode.py|
| TeapotMomentsNodeSetMember | The same as TeapotMomentsNode but for use in a  set | $ORBIT\_ROOT/py/orbit/diagnostics/TeapotDiagnosticsNode.py|
| addTeapotDiagnosticsNode | Adds a diagnostics node to a Teapot lattice | $ORBIT\_ROOT/py/orbit/diagnostics/diagnosticsLatticeModifications.py|
| addTeapotStatLatsNodeSet | Adds a set of StatLats nodes inside a Teapot lattice | $ORBIT\_ROOT/py/orbit/diagnostics/diagnosticsLatticeModifications.py|
| addTeapotMomentsNodeSet | Adds a set of Moments nodes inside a Teapot lattice | $ORBIT\_ROOT/py/orbit/diagnostics/diagnosticsLatticeModifications.py|


# Python Accessible Methods and Variables from BunchTwissAnalysis #
  1. **analyzeBunch(Bunch bunch)**. Analyzes a bunch.  Must be run before routines that get Twiss parameters directly (e.g. not moments or statlats classes). Variables:
    * bunch: A Bunch.
  1. **getCorrelation(int ic, int jc)**. Returns the centered correlation `<(x-<x>)*(y-<y>)> = <x*y> - <x>*<y>`. Variables:
    * ic: Index of the beam phase space coordinate
    * jc: Index of the beam phase space coordinate
  1. **getAverage(int ic)**. Returns the average value for coordinate with index ic.
  1. **getGlobalCount()**. Returns the total number of analysed macroparticles.
  1. **getGlobalMacrosize()**. Returns the total macrosize.
  1. **getEmittance(int ic)**. Returns the rms emittance for index 0,1,2 - x,y,z planes (the pure betatron emittance is calculated for x and y).
  1. **getAlpha(int ic)**. Returns Twiss alpha for index 0,1,2 - x,y,z planes.
  1. **getBeta(int ic)**.  Returns Twiss beta for index 0,1,2 - x,y,z planes.
  1. **getGamma(int ic)**. Returns Twiss gamma for index 0,1,2 - x,y,z planes.
  1. **getEffectiveEmittance(int ic)**. Returns the effective rms emittance for index 0,1 - x,y planes.
  1. **getEffectiveAlpha(int ic)**. Returns effective Twiss alpha for index 0,1 - x,y planes.
  1. **getEffectiveBeta(int ic)**.  Returns effective Twiss beta for index 0,1 - x,y planes.
  1. **getEffectiveGamma(int ic)**. Returns effective Twiss gamma for index 0,1 - x,y planes.
  1. **getDispersion(int ic)**. Returns the Twiss dispersion for index 0,1 - x,y planes.
  1. **getDispersionDerivative(int ic)**. Returns the Twiss dispersion derivative for index 0,1 - x,y planes.
  1. **computeBunchMoments(Bunch bunch, int order)**. Computes the XY moments of the bunch up to a prescribed order.  Does the computation but doesn't return any values. Variables:
    * bunch: A Bunch.
    * order: The highest order of the moments.
  1. **getBunchMoment(int i, int j)**.  Returns the XY moment of the beam. Values will be nonzero only after computeBunchMoments has been run.
    * i: Moment order in X.
    * j: Moment order in Y.

# Python Accessible Methods and Variables from BunchTuneAnalysis #
  1. **assignTwiss(double bx, double ax, double dx, double dpx, double by, double ay)**.  Explicitly assigns Twiss values at the location of the tune node.  The tune calculator requires Twiss input.  Note that pyORBIT and MAD definitions for dispersion differ by a factor of relativistic beta. Variables:
    * bx: Twiss X beta.
    * ax: Twiss X alpha.
    * dx: Twiss X dispersion.
    * dpx: Twiss X dispersion derivative.
    * by: Twiss Y beta.
    * ay: Twiss Y alpha.
  1. **analyzeBunch(Bunch bunch)**.  Does the bunch tune and actions analysis, assuming the Twiss values have been assigned.  Needs to be repeated at least two times in order to arrive at a legitimate tune. This method is most conveniently accessed by adding a TeapotTuneAnalysisNode to a lattice.  The tune information is a particle attribute that gets added to the bunch, will be written out with the particle phase space coordinates when dumpBunch(Bunch bunch) is executed. Variables:
    * bunch: A Bunch.

# Methods and Variables from related diagnostics Python classes #

  1. **StatLats(String filename)**.  From diagnostics.py. The constructor for a StatLats class.
    * filename: String variable representing a filename to be written to.
  1. **writeStatLats(double position, Bunch bunch, double lattlength)**.  A method in StatLats class. Writes the statistical lattice parameters to the file. Columns are: (1) position, (2) turn (or time if lattice length not defined), (3) emitx, (4) emity, (5) betax, (6) betay, (7) alphax, (8) alphay. Variables:
    * position: The position in the lattice
    * bunch: A Bunch.
    * lattlength: The length of the lattice.  Not a required argument. Zero unless user specified.
  1. **TeapotStatLatsNode(String filename)**.  From TeapotDiagnoticsNode.py. The constructor for a Teapot style StatLats class which can be added to a Teapot lattice using the method addTeapotDiagnosticsNode (see below).
    * filename: String variable representing a filenamee to be written to.
  1. **addTeapotStatLatsNodeSet(lattice, filename)**.  This routine is inside the class diagnosticsLatticeModifications.py.  It adds a set of StatLats nodes as child elements throughout the lattice.
    * lattice: A Teapot lattice
    * filename: A string filename.
  1. **Moments(String filename, int order)**.  From diagnostics.py. The constructor for a Moments class.
    * filename: String variable representing a filename to be written to.
    * order: Highest order moment.
  1. **writeMoments(double position, Bunch bunch, double lattlength)**.  A method in the Moments class. Writes the XY moments of the beam up to a user specified order.  Columns are: (1) position, (2) turn (or time if lattice length is not defined), (3 and beyond) Bunch moment (`1, x_ave, y_ave, <x^2>, <xy>, <y^2>, etc...`)
    * position: The position in the lattice
    * bunch: A Bunch.
    * lattlength: The length of the lattice.  Not a required argument. Zero unless user specified.
  1. **TeapotMomentsNode(String filename, int order)**.  From TeapotDiagnoticsNode.py. The constructor for a Teapot style Moments class which can be added to a Teapot lattice using the method addTeapotDiagnosticsNode (see below).
    * filename: A string filename.
    * order: Highest order moment.
  1. **addTeapotMomentsNodeSet(lattice, filename, order)**.  This routine is inside the class diagnosticsLatticeModifications.py.  It adds a set of Moments nodes as child elements throughout the lattice.
    * lattice: A Teapot lattice.
    * filename: A string filename.
    * order: Highest order moment.
  1. **TeapotTuneAnalysisNode()**.  A class in TeapotDiagnostics.py.  It is a Teapot style tune calculator node which can be placed inside a Teapot lattice using the method addTeapotDiagnosticsNode (see below). The node calculates the particle tune and actions for a bunch.  The information will be present as additional particle attributes in the dump bunch files.  The horizontal and vertical  tunes will be columns 9 and 10, respectively, and the actions will be 11 and 12. The assignTwiss method should be employed once this is created.
  1. **addTeapotDiagnosticsNode(lattice, position, diagnostics\_node)**.  This routine is inside the class diagnosticsLatticeModifications.py.  It adds a singleTeapot style diagnostics node (for instance a TeaPotStatLatsNode, TeapotMomentsNode, or TeapotTuneAnalysisNode) to a Teapot Lattice at a user specified position.
    * lattice: A Teapot lattice
    * position: The position location in the lattice (meters).  Must be a drift region.
    * diagnostics\_nod: The diagnostics node.

# Example Scripts #

The following example demonstrates a bunch being propagated through a circular iron collimator of length 0.5 meters and aperture 10 cm whose front face is located 18.5 meters from the start of an existing teapot lattice. It is assumed that the main bunch, "b", is already defined, but that the lost bunch needs creation.  At the end, the main bunch and lost bunch particle parameters are printed.
```



```