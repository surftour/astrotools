#!/bin/tcsh
#$ -j n
#$ -cwd
#$ -pe mpich 4
###$ -q standard
#$ -m be
#$ -M volker@mpa-garching.mpg.de
#$ -N galaxy
#
# now execute my job:

/afs/ipp-garching.mpg.de/i386_rh90/bin/mpirun -np 4  ./P-Gadget3  param.txt 



 
