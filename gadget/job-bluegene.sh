#!/bin/sh
# @ job_name = C06
# @ error = $(job_name).$(jobid).out
# @ output = $(job_name).$(jobid).out
# @ environment = COPY_ALL
# @ wall_clock_limit = 24:00:00
# @ notification = always
# @ notify_user = volker@mpa-garching.mpg.de
# @ job_type = bluegene
# @ bg_size = 512
# @ bg_connection = PREFER_TORUS
# @ queue
#


mpirun -np 1024 -verbose 1 -mode DUAL -exe /ptmp/vrs/Billennium/C06_1200/P-Gadget3/P-Gadget3 -args "param.txt 1" -cwd /ptmp/vrs/Billennium/C06_1200/P-Gadget3/


