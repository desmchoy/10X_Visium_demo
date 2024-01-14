#!/bin/bash

task_name=10X_Visium_preliminary_analysis
task_number=3


if [ ! -d "results_preliminary_analysis" ] ; then
	mkdir results_preliminary_analysis
fi

while read line; do

	qsub -v input=$line -o "$(pwd -P)"/log."${task_name}"."${line}".output -e "$(pwd -P)"/log."${task_name}"."${line}".error sge_"${task_number}"_"${task_name}".batch

done < "sample.list"
