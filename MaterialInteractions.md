# Summary #

This class is a collection of materials scattering methods. The methods apply only to protons in the energy range of 0 - 1.5 GeV. No secondary particles are generated.

None of the methods in this class are directly accessible from the python level. A short description of the main scattering methods is provided below

# Methods #

  1. **momentumKick(double t, double p, double& dpx, double& dpy)**:  This routine applies a random, uniformly distributed 2D momentum kick. See Monte Carlo Methods in PDG, or the K2 code description in the PhD thesis of Nuria Catalan-Lasheras.
  1. **mcsJackson(double stepsize, double z, double a, double rho, long& idum, double beta, double pfac, double& x, double& y, double& px, double& py)**:  Applies a multiple coulomb scattering kick following JD Jackson, Chapter 13, "Collisions Between Charged Particles, Energy Loss, and Scattering", 6th edition.
  1. **ruthScattJackson(double stepsize, double z, double a, double rho, long& idum, double beta, int trackit, double pfac, double& thetax, double& thetay)**:  Applies large angle Rutherford scattering kicks following JD Jackson, Chapter 13 "Collisions Between Charged Particles, Energy Loss, and Scattering", 6th edition.
  1. **elastic\_t(double p, double a, long& idum)** Routine to generate a random momentum transfer for low energy elastic scattering (<= 0.4 GeV):  Follows K2 code described in PhD thesis of Nuria Catalan-Lasheras.
  1. **ionEnergyLoss(double beta, double z, double a)**:  Routine to calculate ionization energy loss based on Bethe-Bloch equation, see PDG Handbook.