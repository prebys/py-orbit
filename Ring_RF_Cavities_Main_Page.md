# Summary #

This page describes the user-accessible python and underlying C++ classes for RF cavities in rings. The models include harmonic accelerating and focusing cavities and also barrier focusing cavities.


# Associated Classes #

| **Name** | **Description** | **Address** |
|:---------|:----------------|:------------|
| [Barrier\_Cav](Barrier_Cav.md) | C++ class containing barrier cavity model - not directly accessible from python script | $ORBIT\_ROOT/src/orbit/RFCavities/ |
| [Frequency\_Cav](Frequency_Cav.md) | C++ class containing linac-type cavity model - not directly accessible from python script | $ORBIT\_ROOT/src/orbit/RFCavities/ |
| [Harmonic\_Cav](Harmonic_Cav.md) | C++ class containing harmonic cavity model - not directly accessible from python script | $ORBIT\_ROOT/src/orbit/RFCavities/ |
| [Dual\_Harmonic\_Cav](Dual_Harmonic_Cav.md) | C++ class containing dual harmonic cavity model - not directly accessible from python script | $ORBIT\_ROOT/src/orbit/RFCavities/ |
| RFNode.Barrier\_RFNode | Python class to invoke barrier cavity model from python script | $ORBIT\_ROOT/py/orbit/rf\_cavities/ |
| RFNode.Frequency\_RFNode | Python class to invoke linac-type cavity model from python script | $ORBIT\_ROOT/py/orbit/rf\_cavities/ |
| RFNode.Harmonic\_RFNode | Python class to invoke harmonic cavity model from python script | $ORBIT\_ROOT/py/orbit/rf\_cavities/ |
| RFNode.Dual\_Harmonic\_RFNode | Python class to invoke dual harmonic cavity model from python script | $ORBIT\_ROOT/py/orbit/rf\_cavities/ |
| RFNode.BRhoDep\_Harmonic\_RFNode | Python class to invoke time-dependent harmonic cavity model from python script. Time dependence is determined by user-provided parameters, with B x Rho determining the energy. | $ORBIT\_ROOT/py/orbit/rf\_cavities/ |
| RFNode.SyncPhaseDep\_Harmonic\_RFNode | Python class to invoke time-dependent harmonic cavity model from python script. Time dependence is determined by user-provided parameters, with the synchronous particle phase determining the energy. | $ORBIT\_ROOT/py/orbit/rf\_cavities/ |
| RFNode.TimeDep\_Barrier\_RFNode | Python class to invoke time-dependent barrier cavity model from python script. Time dependence is determined by user-provided parameters. | $ORBIT\_ROOT/py/orbit/rf\_cavities/ |
| RFLatticeModifications.addRFNode | Python module that allows a user to add an RF cavity node to the teapot lattice in a drift region | $ORBIT\_ROOT/py/orbit/rf\_cavities/ |


# Examples - All found in $ORBIT\_ROOT/examples/RF\_Tests #

| **Name** | **Description** |
|:---------|:----------------|
| barrier\_rf\_cavity\_test.py | Simple test of barrier cavity with parameters constant in time.|
| simple\_rf\_cavity\_test.py | Simple test of harmonic focusing cavity with parameters constant in time.|
| dual\_rf\_cavity\_test.py | Simple test of dual harmonic focusing cavity with parameters constant in time.|
| brhodep\_rf\_cavity\_test.py | Test of accelerating harmonic cavity with time-dependent parameters. Acceleration depends on B x Rho. |
| syncphasedep\_rf\_cavity\_test.py | Test of accelerating harmonic cavity with time-dependent parameters. Acceleration depends on synchronous particle phase. |
| timedep\_barrier\_rf\_cavity\_test.py | Test of barrier focusing cavity with time-dependent parameters. |


# Python Accessible Methods and Variables #
  1. **RFNode.Barrier\_RFNode(ZtoPhi, RFVoltage, RFPhasep, RFPhasem, dRFPhasep, dRFPhasem, length, name)**. Creates an instantiation of a barrier cavity.
    * ZtoPhi: unit conversion factor from longitudinal bunch coordinate in {m} to longitudinal phase {radians}, referenced to the fundamental frequency.
    * RFVoltage: Cavity voltage {GeV}.
    * RFPhasep: Center phase of positive barrier {degrees}.
    * RFPhasem: Center phase of negative barrier {degrees}.
    * dRFPhasep: Phase width of positive barrier {degrees}.
    * dRFPhasem: Phase width of negative barrier {degrees}.
    * length: Node length {0.0 m}.
    * name: Node name {string}.
  1. **RFNode.Frequency\_RFNode(RFFreq, RFE0TL, RFPhase, length, name)**. Creates an instantiation of a linac-type cavity.
    * RFFreq: Cavity frequency {Hz}.
    * RFRFE0TL: Effective cavity accelerating voltage {GeV} - E0 x T x L.
    * RFPhase: Synchronous particle phase {degrees} - dEK = q x E0TL x cos(RFPhase).
    * length: Node length {0.0 m}.
    * name: Node name {string}.
  1. **RFNode.Harmonic\_RFNode(ZtoPhi, dESync, RFHNum, RFVoltage, RFPhase, length, name)**. Creates an instantiation of a harmonic cavity.
    * ZtoPhi: unit conversion factor from longitudinal bunch coordinate in {m} to longitudinal phase {radians}, referenced to the fundamental frequency.
    * dESync: Energy gain of synchronous particle {GeV}, from which energies of all bunch particles are referenced.
    * RFVoltage: Cavity voltage {GeV}.
    * RFHNum: Cavity harmonic number in units of the fundamental frequency.
    * RFPhase: Phase of the cavity {degrees} at the zero of the bunch longitudinal coordinate {m}.
    * length: Node length {0.0 m}.
    * name: Node name {string}.
  1. **RFNode.Dual\_Harmonic\_RFNode(ZtoPhi, RFHNum, RatioRFHNum, RFVoltage, RatioVoltage, RFPhase, RFPhase2, length, name)**. Creates an instantiation of a dual harmonic cavity. See S.Y. Lee: Accelerator Physics. p. 298-305 (1999)
    * ZtoPhi: unit conversion factor from longitudinal bunch coordinate in {m} to longitudinal phase {radians}, referenced to the fundamental frequency.
    * RFHNum: Cavity harmonic number in units of the fundamental frequency for the primary rf system.
    * RatioRFHNum: Cavity harmonic ration (harmonic number of the second rf system to harmonic number of the primary rf system).
    * RFVoltage: Cavity voltage {GeV} for the primary rf system.
    * RatioVoltage: Voltage ration (Voltage of the second rf system to voltage of the primary rf system).
    * RFPhase: Phase angle {degrees} for the synchronous particle for the primary rf system.
    * RFPhase2: Phase angle {degress} for the synchronous particle for the second rf system
    * length: Node length {0.0 m}.
    * name: Node name {string}.
  1. **RFNode.BRhoDep\_Harmonic\_RFNode(ZtoPhi, accelDict, bunch, length, name)**. Creates an instantiation of a harmonic cavity.
    * ZtoPhi: unit conversion factor from longitudinal bunch coordinate in {m} to longitudinal phase {radians}, referenced to the fundamental frequency.
    * dESync: Energy gain of synchronous particle {GeV}, from which energies of all bunch particles are referenced.
    * RFVoltage: Cavity voltage {GeV}.
    * RFHNum: Cavity harmonic number in units of the fundamental frequency.
    * RFPhase: Phase of the cavity {degrees} at the zero of the bunch longitudinal coordinate {m}.
    * length: Node length {0.0 m}.
    * name: Node name {string}.

  1. **RFLatticeModifications.addRFNode(lattice, position, rf\_node)**. Places an RF node into a lattice at a specified position. Variables:
    * lattice: Instantiation of an AccLattice.
    * position: Position from the beginning of the lattice {m}. Must fall inside a drift element.
    * rf\_node: Instantiation of one of the RFNode classes listed above and described below.
  1. **track(Dict paramsDict)**. A method for doing the longitudinal space charge kick. Implemented by the lattice trackBunch method. Variables:
  1. **addLongitudinalSpaceChargeNode(lattice, position, sc1D\_node)**. Routine to add a defined SC1D\_AccNode to a drift region of a teapot lattice. Variables:
    * lattice: The teapot lattice.
    * position: Position of the SC1D\_AccNode element in meters.
    * sc1D\_node: A predefined SC1D\_AccNode





Add your content here.  Format your content with:
  * Text in **bold** or _italic_
  * Headings, paragraphs, and lists
  * Automatic links to other wiki pages