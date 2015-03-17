# Summary #

This page describes how to create a kicker with a user specified kick and custom waveform.  Currently, waveforms are available for flattop kickers (constant kick), and square-root time waveforms.

# Associated Classes #

| **Name** | **Description** | **Address** |
|:---------|:----------------|:------------|
| [XKicker](XKicker.md) | The class which creates a horizontal kicker | $ORBIT\_ROOT/py/orbit/kickernodes/ |
| [YKicker](YKicker.md) | The class which creates a vertical kicker | $ORBIT\_ROOT/py/orbit/kickernodes/ |
| TeapotXKickerNode | The python module that allows a user to define Teapot horizontal kicker node | $ORBIT\_ROOT/py/orbit/kickernodes/ |
| TeapotYKickerNode | The python module that allows a user to define Teapot vertical kicker node | $ORBIT\_ROOT/py/orbit/kickernodes/ |
| KickerLatticeModifications | The python module that allows a user to add a defined TeapotX(Y)KickerNode to a teapot lattice in a drift region | $ORBIT\_ROOT/py/orbit/kickernodes/ |
| flatTopWaveform | The python module for creating a flat top kicker waveform | $ORBIT\_ROOT/py/orbit/kickernodes/ |
| rootTWaveform | The python module for creating a square root of time kicker waveform | $ORBIT\_ROOT/py/orbit/kickernodes/ |

# Python Accessible Methods and Variables #
  1. **[SC1D\_AccNode](SC1D_AccNode.md)(nt b\_a, double length, int nMacrosMin, int useSpaceCharge, int nBins)**. Creates a longitudinal space charge node. Variables:
    * b\_a: Approximate ratio of beam pipe radius to beam radius.
    * length: Length of the lattice
    * nMacrosMin: Minimum number of macroparticles needed to do the computation
    * useSpaceCharge: Flag for turning the space charge piece on (1) or off (0).
    * nBins: Number of longitudinal slicing bins for space charge kick calculation.
  1. **assignImpedance(pyObject py\_complex\_arr)**. Convenience method for assigning an impedance array.
    * py\_complex\_arr: A python array of complex numbers which are the impedances in units of Ohms/n (n is mode number). For inductive impedances, convention is that real part is positive and imaginary part negative.
  1. **track(Dict paramsDict)**. A method for doing the longitudinal space charge kick. Implemented by the lattice trackBunch method. Variables:
  1. **addLongitudinalSpaceChargeNode(lattice, position, sc1D\_node)**. Routine to add a defined SC1D\_AccNode to a drift region of a teapot lattice. Variables:
    * lattice: The teapot lattice.
    * position: Position of the SC1D\_AccNode element in meters.
    * sc1D\_node: A predefined SC1D\_AccNode


# Example Scripts #

The following example demonstrates adding a longitudinal space charge node to a lattice.  Example can be found in ORBIT\_ROOT/py-orbit/examples/Space\_Charge/sc1d/longscbunchlattice.py
```