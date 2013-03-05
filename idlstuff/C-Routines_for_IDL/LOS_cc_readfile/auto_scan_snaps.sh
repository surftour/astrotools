#!/bin/bash


mkdir ./NH_results
NH_file_header=./NH_results/NH.
snapshots_dir=.
snapshot_header=snapshot_

j=0
while [ $j -le 115 ]
do
  i=`expr $j`
  if [ $j -lt 100 ] ; then 
    i=0`expr $j`
  fi
  if [ $j -lt 10 ] ; then 
    i=00`expr $j`
  fi
  #echo $snapshots_dir/$snapshot_header$i $NH_file_header$i 
  echo $snapshots_dir/$snapshot_header$i $NH_file_header$i > temp_inputline
  sim < temp_inputline > temp_outputline
  j=`expr $j + 1`
done

rm temp_inputline
rm temp_outputline

