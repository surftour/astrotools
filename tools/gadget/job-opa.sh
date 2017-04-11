#!/bin/tcsh
#$ -j n
#$ -cwd
#$ -pe mvapich2 8
#$ -m be
#$ -M volker@mpa-garching.mpg.de
#$ -N test
#

mpiexec -np $NSLOTS  ./P-Gadget3  param-galaxy.txt 

if (-e ../cont ) then
#    rm -f "../cont"
#    qsub job-opa.sh
endif






