#!/usr/local/bin/tcsh
#$ -j n
#$ -cwd
#$ -pe infiniband0  4
#$ -m be
#$ -M volker@mpa-garching.mpg.de
#$ -N simulation
#

setenv PYTHONHOME "/usr/local/appl/Python-2.4.3"
setenv PATH "/opt/mvapich2/bin:${PATH}"

mpiexec -machinefile $TMPDIR/machines -np $NSLOTS  ./P-Gadget3  param.txt 






