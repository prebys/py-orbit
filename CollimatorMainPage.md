# Summary #
The collimator class and auxiliary classes model the propagation of protons through various materials. The module only applies to protons in the energy range of 0 - 1.5 GeV. Other particle species are not handled, and there are no secondary particles generated.

There are five available collimator shapes: cirlce, elliptical, one sided flat, two sided flat, and rectangular.  The shape and aperture of the collimator are specified by the user.  0 aperture corresponds to a solid slab of material.

Six materials are available for the collimator: 1=carbon, 2=aluminum, 3=iron, 4=copper, 5=tantalum, 6=tungstun, 7=platinum, 8=lead, 9 = black absorber.  The black absorber acts like an aperture which absorbs a particle on impact. There is also a user-specified density factor for each material, useful for modeling materials that are layered with water or air.

The module requires input of two bunches from the Bunch class: one bunch for the alive particles being propagated, and one bunch to be populated with particles lost in the collimator. Particles are lost in the collimator when their energy is below 20 MeV, or when they suffer a nuclear inelastic scattering event.

A collimator must be placed in a drift section whose length is equal to or greater than the collimator length. The space in the drift will be re-allocated to the collimator and surrounding drift region(s) such that the original drift length of the section is preserved.

# Associated Classes #

| **Name** | **Description** | **Address** |
|:---------|:----------------|:------------|
| [Collimator](Collimator.md) | The class which defines the collimator | $ORBIT\_ROOT/src/orbit/MaterialInteractions/ |
| MaterialInteractions | Class containing a collection of particle scattering methods | $ORBIT\_ROOT/src/orbit/MaterialInteractions/ |
|cross\_sections | Contains the elastic and inelastic cross sections for each material (barns), and material properties (MKS) | $ORBIT\_ROOT/src/orbit/MaterialInteractions/|
|num\_recipes | Numerical recipes used in scattering models | $ORBIT\_ROOT/src/orbit/MaterialInteractions/ |
| wrap\_collimator | The python wrapper class for Collimator |  $ORBIT\_ROOT/src/orbit/MaterialInteractions/ |
| TeapotCollimatorNode | The python module that allows a user to  define collimator | $ORBIT\_ROOT/py/orbit/collimation/ |
| collimatorLatticeModifications | The python module that allows a user to  add a defined collimator to the lattice in a drift region | $ORBIT\_ROOT/py/orbit/collimation/ |

# Python Accessible Methods and Variables #
  1. **TeapotCollimatorNode(int length, int ma, double density\_fac, int shape, double  a, double b, double c, double d, double angle, string name)**. An instance of a collimator. Creates a teapot style collimator and contains methods to track through this collimator. Variables:
    * length: Length of the collimator in meters
    * ma: Material type index: 1=carbon, 2=aluminum, 3=iron, 4=copper, 5=tantalum, 6=tungstun, 7=platinum, 8=lead, 9 = black absorber.
    * density\_fac: A multiplier on the density of chosen material. Defaults to 1.
    * shape: Shape type index: 1=circle, 2=ellipse, 3=one sided flat, 4=two sided flat, 5=rectangular (outside is collimator), 6=rectangular (inside is collimator).
    * a: Radius for circular collimator, x-axis for elliptical, horizontal right side extent for flat and rectangular.
    * b: Y-axis for elliptical, left side horizontal extent for two-sided flat or rectangular.
    * c: Upper side vertical extent for rectangular.
    * d: Lower side vertical extent for rectangular.
    * angle: Tilt angle of the collimator defined in degrees.
    * name: The name of the collimator.
  1. **addTeapotCollimatorNode(TEAPOT\_Lattice lattice, double position, TeapotCollimatorNode collimator\_node)**. A method to add the collimator to a Teapot lattice. Variables:
    * lattice: The Teapot lattice the collimator should be added to.
    * position: The start position of the collimator in the lattice.
    * collimator\_node: A previously defined collimator

# Example Scripts #

The following example demonstrates a bunch being propagated through a circular iron collimator of length 0.5 meters and aperture 10 cm whose front face is located 18.5 meters from the start of an existing teapot lattice. It is assumed that the main bunch, "b", is already defined, but that the lost bunch needs creation.  At the end, the main bunch and lost bunch particle parameters are printed.
```


length = 0.5
ma = 3
density_fac = 1.0
shape = 1
a = 0.01
b = 0
c = 0
d = 0
angle = 0

collimator = TeapotCollimatorNode(length, ma, density_fac, shape, a, b, c, d, angle, "Collimator 1")

addTeapotColimatorNode(teapot_latt, 18.5, collimator) 

#------------------------------
#Lost bunch initialization
#------------------------------

lostbunch = Bunch()
lostbunch.addPartAttr("LostParticleAttributes")
paramsDict = {}
paramsDict["lostbunch"]=lostbunch
paramsDict["bunch"]= b

#=====track bunch through Collimator Node============
collimator.track(paramsDict)

# dump ORBIT_MPI bunch to compare results
bunch.dumpBunch()
lostbunch.dumpBunch()
print "Stop."


```

In the following piece of code a short 1.5 cm, one-sided flat copper collimator with aperture 5 cm is defined. The collimator is tilted at a 45 degree angle.

```

length = 0.015
ma = 4
density_fac = 1.0
shape = 3
a = 0.05
b = 0
c = 0
d = 0
angle = 45

collimator = TeapotCollimatorNode(length, ma, density_fac, shape, a, b, c, d, angle, "Collimator Scraper")

```

Complete run examples can be found in $ORBIT\_ROOT/examples/Collimation/ORBIT\_Benchmarks/ .