The JAM (Jeans Anisotropic MGE) package implements the solution of the
anisotropic Jeans equations presented in Cappellari (2008, MNRAS, 390, 71).

To construct a dynamical model for a specific galaxy one first needs to
obtain a Multi-Gaussian Expansion parametrization for its surface
brightness. This can be derived with the MGE_FIT_SECTORS package,
available from http://www-astro.physics.ox.ac.uk/~mxc/idl/

The JAM package consists of three completely independent routines, plus
a set of auxiliary mathematical routines (in the sub directory /utils).

The three routines of the JAM package are:

1.  JAM_AXISYMMETRIC_RMS: Computes a prediction for
    V_RMS=sqrt(V^2+sigma^2) for an axisymmetric galaxy;

2.  JAM_AXISYMMETRIC_VEL: Computes a prediction for V;

3.  JAM_SPHERICAL_RMS: Computes a prediction for
    V_RMS=sqrt(V^2+sigma^2) for a spherical galaxy;

Also inclued in the JAM package is the following independent routine,
which turns out to be useful in many cases:

-   MGE_CIRCULAR_VELOCITY: Computes the circular velocity in the
    equatorial plane of an axisymmetric galaxy.

To learn how to use a routine you should compile it and call the
corresponding example program TEST_xxx, where xxx is one of the three
above names. The example routines can be used as starting point to run
your own models.

Michele Cappellari
Oxford, 12 August, 2008
