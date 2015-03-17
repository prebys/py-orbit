# Summary #

This is the Teapot implementation of the InjectParts class. This must be used in order to add an injection node (InjectParts class) as element in a Teapot lattice.

# Python Accessible Methods #

  1. **TeapotInjectionNode(nparts, bunch, lostbunch, foilparams, xDistFunc, yDistFunc, lDistFun)**. Creates the node.
    * nparts: Number of macro particles to be injected during each pass through the element.
    * bunch: The bunch to be injected into.
    * lostbunch: The lost particle bunch which will be populated by injected particles outside of the injection region.
    * injectregion: A list of ransverse limits of the injection region (xmin: Horizontal minimum; xmax: Horizontal maximum; ymin: Vertical minimum; ymax: Vertical maximum).
    * xDistFunc: Horizontal distribution function. One example is JohoTransverse.
    * yDistFunc: Vertical distribution function: One example is JohoTransverse.
    * lDistFunc: Longitudinal distribution function.  Examples are JohoLongitudinal and UniformLongDist.
  1. **track(Dict paramsDict)**. A method for injecting particles. Variables:
    * paramsDict: A dictionary which should contains two items: 1) the main bunch to be injected into, and 2) the bunch to be populated with particle that were outside of the valid injection region.