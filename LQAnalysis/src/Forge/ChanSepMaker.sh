#!/bin/bash
#Usage : ./shell <FileToCopy>
#     ex)./shell Aug.cc

if [[ -z $1 ]]; then echo "Write down the file to copy"; exit 1; fi;
if [[ -z $LQANALYZER_ANALYSIS_PATH ]]; then "source setup.sh needed"; exit 1; fi;
FileToCp=$( basename $1 .cc )

cp ${FileToCp}.cc ${FileToCp}_1e2mu.cc
find . -name "${FileToCp}_1e2mu.cc" -type f -exec sed -i s/${FileToCp}/${FileToCp}_1e2mu/g {} +
cp ${FileToCp}.cc ${FileToCp}_3mu.cc
find . -name "${FileToCp}_3mu.cc" -type f -exec sed -i s/${FileToCp}/${FileToCp}_3mu/g {} +
cp $LQANALYZER_ANALYSIS_PATH/include/${FileToCp}.h $LQANALYZER_ANALYSIS_PATH/include/${FileToCp}_1e2mu.h
find $LQANALYZER_ANALYSIS_PATH/include/ -name "${FileToCp}_1e2mu.h" -type f -exec sed -i s/${FileToCp}/${FileToCp}_1e2mu/g {} +
cp $LQANALYZER_ANALYSIS_PATH/include/${FileToCp}.h $LQANALYZER_ANALYSIS_PATH/include/${FileToCp}_3mu.h
find $LQANALYZER_ANALYSIS_PATH/include/ -name "${FileToCp}_3mu.h" -type f -exec sed -i s/${FileToCp}/${FileToCp}_3mu/g {} +



sed -i -e s/emumu_analysis=false/emumu_analysis=true/g -e s/trimu_analysis=true/trimu_analysis=false/g ${FileToCp}_1e2mu.cc
sed -i -e s/emumu_analysis=true/emumu_analysis=false/g -e s/trimu_analysis=false/trimu_analysis=true/g ${FileToCp}_3mu.cc
