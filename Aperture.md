#summary Description of the Aperture class.

# Summary #

This class implements a method to check if particles are outside of a user-specified aperture.  Particles outside of the aperture are removed and added to a lost particles bunch.  The method requires the presence of a main bunch that is tracked and a second bunch for housing the lost particles.

# Python Accessible Methods and Variables #
  1. **TeapotApertureNode(int shape, double a, double pos=0, double b=0, double c=0, double d=0)**. An instance of an aperture. Variables:
    * shape: Shape type index. 1=circle, 2=ellipse, 3=rectangle.
    * a: Radius for circular collimator, x-axis half length for ellipse and rectangle.
    * b: Y-axis half length for elliptical and rectangle. Doesn't apply for circle.
    * pos: Position of aperture.  Will be overwritten if aperture is added to a lattice. Optional argument.
    * c: Horizontal offset.Optional argument.
    * d: Vertical offset. Optional argument.
  1. **setPosition(double position)**. The method sets the position of the collimator in the lattice.
    * position: The start position of the collimator in the lattice.


