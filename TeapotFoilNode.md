# Summary #

This is the Teapot implementation of the foil class. This must be used in order to add a foil as an element in a Teapot lattice. A convenience routine is defined for switching between particle scattering algorithms.


# Python Accessible Methods #

  1. **track(_Dict_ paramsDict)**. A method for tracking the particles through the foil. Variables:
    * paramsDict: A dictionary which should contain one or two items: 1) the main bunch, and 2) the bunch to be populated with lost particles if the full scattering model is used.
  1. **setScatterChoice(int choice)**. Routine used to set the scattering method via an integer index. 0 corresponds to the full scattering model, and 1 corresponds to the simplified scattering model. Further information about the models can be found in the [Foil](Foil.md) class. The default is 0 (full scattering).
    * choice: The integer index of the scattering method choice. 0 for full scattering model and 1 for simple scattering model.