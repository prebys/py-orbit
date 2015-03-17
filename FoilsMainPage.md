# Summary #

This page describes how to track particles through foils. The Foil class allows the creation of rectangular carbon foils of user defined thickness and size. An arbitrary number of foils can be added to any lattice. Python convenience scripts are available for adding foils to the lattice. From the standpoint of the lattice, foils are considered as zero length elements.

Two different particle scattering methods are available in the Foil class. The default method is the full scattering method (index 0) which allows particle loss, and the alternative choice is a simple scattering method (index 1) which does not allow particle loss. More information on these routines can be found in the Foil class.

# Associated Classes #

| **Name** | **Description** | **Address** |
|:---------|:----------------|:------------|
| [Foil](Foil.md) | The class which defines the foil | $ORBIT\_ROOT/src/orbit/MaterialInteractions/ |
| MaterialInteractions | Class containing a collection of particle scattering methods | $ORBIT\_ROOT/src/orbit/MaterialInteractions/ |
|cross\_sections | Contains the elastic and inelastic cross sections for each material (barns), and material properties (MKS) | $ORBIT\_ROOT/src/orbi/MaterialInteractions/|
|num\_recipes | Numerical recipes used in scattering models | $ORBIT\_ROOT/src/orbi/MaterialInteractions/ |
| wrap\_foil| The python wrapper class for Foil |  $ORBIT\_ROOT/src/orbi/MaterialInteractions/ |
| TeapotFoilNode | The python module that allows a user to define Teapot foil | $ORBIT\_ROOT/py/orbit/foil/ |
| foilLatticeModifications | The python module that allows a user to  add a defined foil to the teapot lattice in a drift region | $ORBIT\_ROOT/py/orbit/foil/ |


# Python Accessible Methods and Variables #
  1. **TeapotFoilNode(double xmin, double xmax, double ymin, double ymax, double thick, string name)**. Creates a teapot implementation of foil. Has a track method and a setScatterChoice method. Variables:
    * xmin: Horizontal minimum of foil.
    * xmax: Horizontal maximum of foil.
    * ymin: Vertical minimum of foil.
    * ymax: Vertical maximum of foil.
    * thick: Thickness of foil (ug/cm^2)
    * name: User defined name of the foil.
  1. **addTeapotFoilNode(TEAPOT\_Lattice lattice, double position, TeapotFoilNode foil\_node)**. A method to add the collimator to a Teapot lattice. Variables:
    * lattice: The Teapot lattice the collimator should be added to.
    * position: The start position of the collimator in the lattice.
    * foil\_node: A previously defined teapot style foil.


# Example Scripts #

The following example demonstrates how a carbon 400 ug/cm^2 teapot style foil is created and then added to an existing lattice, 18.5 meters downstream of the lattice start point. The simplified scattering model is chosen and an existing bunch is transported through the foil.

```

xmin = -0.05
xmax = 0.05
ymin = -0.05
ymax = 0.05
thick = 400

foil = TeapotFoilNode(xmin, xmax, ymin, ymax, thick, "Foil 1")

addTeapotFoilNode(teapot_latt,18.5,foil) 

#=====track bunch through Collimator Node============
paramsDict = {}
paramsDict["bunch"]= b

scatterchoice = 1 //Use simple scattering model. 
foil.setScatterChoice(scatterchoice)
foil.track(paramsDict)

```

The same example as above except this time use the full scattering model and allow particle loss.

```

xmin = -0.05
xmax = 0.05
ymin = -0.05
ymax = 0.05
thick = 400

foil = TeapotFoilNode(xmin, xmax, ymin, ymax, thick, "Foil 1")

addTeapotFoilNode(teapot_latt,18.5,foil) 

lostbunch = Bunch()
lostbunch.addPartAttr("LostParticleAttributes") 
paramsDict = {}
paramsDict["bunch"]= b
paramsDict["lostbunch"]=lostbunch

scatterchoice = 0
foil.setScatterChoice(scatterchoice)
foil.track(paramsDict)

```

Complete script examples can be found in $ORBIT\_ROOT/py-orbit/examples/Foils