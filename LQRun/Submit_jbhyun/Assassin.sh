#!/bin/bash

declare -a TargetList=("period")
QueueOption="All" #All or fastq or longq
if [[ ${QueueOption} != "All" ]] && [[ ${QueueOption} != "fastq" ]] && [[ ${QueueOption} != "longq" ]];
then echo "Wrong Queue Option, exit"; exit 1; fi;

for Target in "${TargetList[@]}";
do
  if [[ ${QueueOption} != "All" ]]; then
    qstat | egrep ${Target} | egrep ${QueueOption} > TargetList.txt
  else qstat | egrep ${Target} > TargetList.txt
  fi;
  
  while read line
  do
    JobIdx=$( echo $line | cut -d ' ' -f1 )
    qdel ${JobIdx}
  done<TargetList.txt
done
