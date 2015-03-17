# Summary #

This is the Teapot implementation of a collimator node.  Must be used in order to create a collimator that can be located within a Teapot lattice. Has a convenience routine for tracking the particles through the collimator.

# Python Accessible Methods #

  1. **track(_Dict_ paramsDict)**. A method of the TeapotCollimatorNode which tracks the particles through the collimator. An implementation AccNodeBunchTracker for the teapot collimator. Variables:
    * paramsDict: A dictionary which should contain two items: 1) the main bunch to be collimated, and 2) the bunch to be populated with lost particles, which should have ParticleAttribute "LostParticleAttributes"
  1. **setPosition(double position)**. The method sets the position of the collimator in the lattice.
    * position: The start position of the collimator in the lattice.