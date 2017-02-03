#!/bin/bash
# Author : Jihwan Bhyun, 04.Jan.2017
# Context : make a new code frame work
# Usage : ./shell original file newcyclefilename

if [[ -z $1 || -z $2 ]]; then echo "You missed either one of base of target. exited"; exit 1; fi
if [[ ! -d $LQANALYZER_SRC_PATH ]]; then echo "sourse setup.sh needed. exited"; exit 1; fi
if [[ $( pwd ) != ${LQANALYZER_SRC_PATH}Forge ]]; then echo "Running at wrong place. exited"; exit 1; fi
if [[ $1 = *"/"* ]]; then echo "Put basefile into the Forge directory. exited"; exit 1; fi

BaseCycle=$( basename $1 .cc )
if [[ $2 == *".cc" ]]; then TargetCycle=$( basename $2 .cc )
else TargetCycle=$2
fi

echo "From $BaseCycle >> To $TargetCycle"

cp ${BaseCycle}.cc ${TargetCycle}.cc
find . -name "${TargetCycle}.cc" -type f -exec sed -i s/${BaseCycle}/${TargetCycle}/g {} +

cp ${LQANALYZER_INCLUDE_PATH}/${BaseCycle}.h ${LQANALYZER_INCLUDE_PATH}/${TargetCycle}.h
find ${LQANALYZER_INCLUDE_PATH} -name "${TargetCycle}.h" -type f -exec sed -i s/${BaseCycle}/${TargetCycle}/g {} +

rm $1

sed -i "/Jihwan Bhyun Modification/a#pragma link C++ class ${TargetCycle}+;" ${LQANALYZER_INCLUDE_PATH}/LQAnalysis_LinkDef.h

echo "Fini."
