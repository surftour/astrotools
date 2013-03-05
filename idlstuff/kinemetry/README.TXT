This distribution contains an IDL (http://www.rsinc.com/)
implementation of the KINEMETRY method developed by
Krajnovic et al. 2006 (astro-ph/0512200).

-------------
SAURON velocity measurements and associated errors for NGC2974 are
given in file NGC2974_SAURON_kinematics.dat. There are five columns:
number of bin, X and Y bin coordinates, velocity and velocity error.

These data correspond to data published in Emsellem et al. (2004
MNRAS, 352, 271). These data are not part of the official SAURON data
release and may differ slightly from the final products.


To learn how to use KINEMETRY, run the example program
kinemetry_example.pro or read the detailed documentation at the top of
the file kinemetry.pro

To run the example KINEMETRY_EXAMPLE:
1. put all files in the directory 'kinemetry'
2. cd to that directory (within IDL: CD, kinemetry')
3. compile kinemetry_example.pro routine (.r kinemetry_example.pro)
4. type kinemetry_example at the IDL prompt

The procedure KINEMETRY needs the following extra IDL routines:
- MPFIT: by C.B. Markwardt -> http://astrog.physics.wisc.edu/~craigm/idl/

The IDL Astronomy User's Library is assumed to be installed
(http://idlastro.gfsc.nasa.gov/) and in the IDL path. In particular
the program KINEMETRY_EXAMPLE uses the routines RDFLOAT.PRO and
PLOTERROR.PRO.


The program was tested on IDL 6.0.

-------------

Davor Krajnovic
Oxford, 07 December 2005


================================================================

TJ downloaded and played with this in Feb & Mar of 2009.  We will
end up using the tj_kinemetry.pro file which has minimal changes
from the original kinemetry.pro.  And this will be copied up to
our root IDL directory.


