#!/bin/sh

###### SET WHAT JOBS TO RUN
runMC=true
runData=false


## RUN PARAMETERS
AnalysisCode="Sep2016_MuMubj_DataMCComp"  ##"Sep2016_3l4j_DataMCComp" ###_D3lv" ###"Jun2016_3l4j_SigBasic" ###KPS ###"ExampleAnalyzerDiMuon"
Stream="DoubleMuon" ###"SingleElectron" ### DoubleMuon DoubleEG MuonEG SingleElectron SingleMuon 
Skim="SKTree_DiLepSkim" ###"Lepton" ### "Lepton"(single lepton skim)/"DiLep"(dilepton skim)/"NoCut"(noskim) ### SKTree_NoSkim/SKTree_LeptonSkim/SKTree_Di[Tri]LepSkim/ flatcat
DataPeriod="CtoD"   ###  "C" = period C only   "ALL" or "CtoD"  = period C+D
job_logstep=1000
LogLevel="INFO"
Njob=100

#MCList="MuMubj_bkg"
MCList="QCD_All"
#MCList="Analysis_bkg_1e2mu"
#MCList="Analysis_bkg_3mu"
#MCList="Analysis_bkg"
#MCList="Analysis_bkg_TT"
#MCList="Analysis_sig_test1e2mu"
#MCList="Analysis_sig_test3mu"
###Backgound : AllSample / Analysis_bkg / Analysis_bkg_test / QCD_mu
###Signal    : Analysis_sig_All / Analysis_sig_1e2mu / Analysis_sig_3mu / tthwA_1e2mu / tthwA_3mu / Analysis_sig_test / Analysis_sig_test1

########################################################################
outputdir_cat="/data2/Users/jbhyun/cmssw/CATanalyzer/CATOut/${CATVERSION}/"

if [[ ${Stream} == 'DoubleMuon' || ${Stream} == 'SingleMuon' ]]; then outputdir_lep=${outputdir_cat}"/output/Muon/";
elif [[ ${Stream} == 'DoubleEG' || ${Stream} == 'SingleElectron' ]]; then outputdir_lep=${outputdir_cat}"/output/Electron/";
elif [[ ${Stream} == 'MuonEG' ]]; then outputdir_lep=${outputdir_cat}"/output/ElectronMuon/";
else echo "Error : job stream wrong. Check."; exit 1;
fi

if [[ ! -d "${outputdir_cat}" ]]; then echo "${outputdir_cat} doesn't exist, run OutDirSetup.sh, exit"; exit 1; fi
if [[ ! -d "${outputdir_lep}" ]]; then echo "${outputdir_lep} doesn't exist, run OutDirSetup.sh, exit"; exit 1; fi


if [[ $DataPeriod  == "C" ]]; then dir_period="periodC/";
elif [[ $DataPeriod == "D" ]]; then dir_period="periodD/";
elif [[ $DataPeriod  == "ALL" || $DataPeriod == "CtoD" ]]; then dir_period="periodCtoD/";
else echo "Error: Period Set Wrongly"; exit 1;
fi

OutputDir=${outputdir_lep}${dir_period}

if [[ ! -d "${OutputDir}" ]]; then mkdir ${OutputDir}; echo "Made ${OutputDir}"; fi

if [[ $runData == 'false' || $runData == 'False' ]]; then Stream=""; fi
if [[ $runMC == 'false' || $runMC == 'False' ]]; then MCList=""; fi

########################################################################
if [[ $runData == 'false' || $runData == 'False' ]];
  then bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -a ${AnalysisCode} -p ${DataPeriod} -s ${Skim} -list ${MCList} -n ${Njob} -o ${OutputDir} -d ${LogLevel} &> SKTree_${AnalysisCode}_${Stream}_log.txt & 
elif [[ $runMC == 'false' || $runMC == 'False' ]];
  then bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -a ${AnalysisCode} -S ${Stream} -p ${DataPeriod} -s ${Skim} -n ${Njob} -o ${OutputDir} -d ${LogLevel} &> SKTree_${AnalysisCode}_${Stream}_log.txt &
else
  bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -a ${AnalysisCode} -S ${Stream} -p ${DataPeriod} -s ${Skim} -list ${MCList} -n ${Njob} -o ${OutputDir} -d ${LogLevel} &> SKTree_${AnalysisCode}_${Stream}_log.txt &
fi
########################################################################


