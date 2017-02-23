#!/bin/bash
# 23. Dec. 2016 Author: Jihwan Bhyun
# context: When upgrade tags, you may need to copy your private function to AnalyzerCore(Those functions meant to be private not wanted to be included in central ver.)
#          Currently This is intended for copy LQAnalysis_LinkDef and AnalyzerCore
# Usaage : after type in appropriate file.(And if you add some more, check the code for usage, this is designed only for analyzer core and LQlinkdef"
#          ./shell

declare -a TargetFileList=("AnalyzerCore.cc" "AnalyzerCore.h" "LQAnalysis_LinkDef.h")

if [[ -z $CATTAG ]]; then echo "source setup needed, exit"; exit 1; fi;
if [[ ! -d ${CatBackup}/${CATTAG} ]]; then echo "made ${CatBackup}/${CATTAG}"; mkdir ${CatBackup}/${CATTAG}; fi;
if [[ ! -d ${CatBackup}/${CATTAG}/LQAnalysis ]]; then echo "made ${CatBackup}/${CATTAG}/LQAnalysis"; mkdir ${CatBackup}/${CATTAG}/LQAnalysis; fi;
if [[ ! -d ${CatBackup}/${CATTAG}/LQAnalysis/include ]]; then echo "made ${CatBackup}/${CATTAG}/LQAnalysis/include"; mkdir ${CatBackup}/${CATTAG}/LQAnalysis/include; fi;
if [[ ! -d ${CatBackup}/${CATTAG}/LQAnalysis/src ]]; then echo "made ${CatBackup}/${CATTAG}/LQAnalysis/src"; mkdir ${CatBackup}/${CATTAG}/LQAnalysis/src; fi;

for TargetFile in "${TargetFileList[@]}"
do
  if [[ ${TargetFile} == *".cc" ]]; then
    if [[ ! -e ${LQANALYZER_SRC_PATH}/${TargetFile} ]]; then echo "Target src file doesn't exist, something wrong, exited"; exit 1; fi;
    if [[ ! -e ${CatBackup}/${CATTAG}/LQAnalysis/src/${TargetFile} ]]; then echo "Source src file doesn't exist, something wrong, exited"; exit 1; fi;

    SourceFile=${CatBackup}/${CATTAG}/LQAnalysis/src/${TargetFile}
    TargetFile=${LQANALYZER_SRC_PATH}/${TargetFile}
    LineN_Start=$( nl -ba ${SourceFile} | egrep "Jihwan Bhyun Modification" | sed "s/^[\ ]*//g" | cut -d$'\t' -f1 )
    LineN_End=$( wc -l ${SourceFile} | cut -d ' ' -f1 )
    sed -n ${LineN_Start},${LineN_End}p ${SourceFile} >> ${TargetFile}

  elif [[ ${TargetFile} == *".h" ]]; then
    if [[ ! -e ${LQANALYZER_INCLUDE_PATH}/${TargetFile} ]]; then echo "Target .h file doesn't exist, something wrong, exited"; exit 1; fi;
    if [[ ! -e ${CatBackup}/${CATTAG}/LQAnalysis/include/${TargetFile} ]]; then echo "Source .h file doesn't exist, something wrong, exited"; exit 1; fi;

    SourceFile=${CatBackup}/${CATTAG}/LQAnalysis/include/${TargetFile}
    TargetFile=${LQANALYZER_INCLUDE_PATH}/${TargetFile}
    LineN_Start=$( nl -ba ${SourceFile} | egrep "Jihwan Bhyun Modification" | sed "s/^[\ ]*//g" | cut -d$'\t' -f1 )
    LineN_SrcEnd=$( wc -l ${SourceFile} | cut -d ' ' -f1 )
    LineN_TgtEnd=$( wc -l ${TargetFile} | cut -d ' ' -f1 )
    LineN_TgeEndm1=$(( ${LineN_TgtEnd} - 1 ))
    sed -i "${LineN_TgeEndm1},${LineN_TgtEnd}s/^/\/\//g" ${TargetFile}
    sed -n ${LineN_Start},${LineN_End}p ${SourceFile} >> ${TargetFile}

  else echo "not .cc, not .h files. Check. exited"; exit 1; fi;
done
