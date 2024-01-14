#!/bin/bash

#####
#
# REQUIRED AT ROOT LEVEL:
# -----
# spaceranger_inputs directory made
#
#####



ls ./RAW -I source | sed 's/\(^.*\)-.*/\1/g' | sort | uniq > sample.list



#if [ ! -d 'spaceranger_inputs' ] ; then
#	mkdir spaceranger_inputs
#fi



