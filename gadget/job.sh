#!/usr/bin/sh
# @ output = job.$(jobid).out
# @ error = job.$(jobid).err
# @ class = chubby
# @ job_type = parallel
# @ node_usage = shared
# @ node = 1
# @ tasks_per_node = 32
# @ resources = ConsumableCpus(1)
# @ notify_user = volker@mpa-garching.mpg.de
# @ notification = always
# @ initialdir = /afs/ipp-garching.mpg.de/u/vrs/PostProcessing/P-Gadget3/
# @ network.MPI = css0,not_shared,us
# @ queue

export MP_EUILIB=us
export MP_EUIDEVICE=css0
export MP_SHARED_MEMORY=yes
export MP_SINGLE_THREAD=yes
export MEMORY_AFFINITY=MCM
export MP_EAGER_LIMIT=32K

poe ./P-Gadget3 param_800.txt 3 110 -procs 32

