#!/bin/bash 
#PBS -j oe 
#PBS -S /bin/bash 
#PBS -l select=8:ncpus=1 
#PBS -l walltime=12:00:00 
#PBS -N halo1
#PBS -M volker@mpa-garching.mpg.de
#PBS -m abe 

. /etc/profile 

cd /home/hlrb2/h1091/h1091aa/test1/P-Gadget3/

mpiexec ./P-Gadget3 param.txt 1  >& mylog.$PBS_JOBID.txt 

