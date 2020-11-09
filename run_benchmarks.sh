#!/bin/bash

./download_data.sh

command=$1

#for file in data_dir:
#run alignment algorithm

#time some of them - or all of them?
/usr/bin/time {command}

#valgrind some of them
valgrind {command}

#do variant calling

#do rtg_tools

