#!/bin/bash

#####
#
# REQUIRED AT ROOT LEVEL:
# -----
# sample.list
#
#####


task_name=spaceranger
task_number=2

if [ ! -d "results_spaceranger" ] ; then
	mkdir results_spaceranger
fi



while read line ; do
	
	qsub -v sample=$line -o "$(pwd -P)"/log."${task_name}"."${line}".output -e "$(pwd -P)"/log."${task_name}"."${line}".error sge_"${task_number}"_"${task_name}".batch

done < "sample.list"
