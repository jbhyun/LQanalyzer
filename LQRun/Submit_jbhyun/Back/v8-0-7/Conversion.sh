#!/bin/bash

########################################################################
## MC / DATA
runMC=true
runData=false


########################################################################
## RUN PARAMETERS

AnalysisCode="Jun2017_ConversionStudy" 
Stream="MuonEG"     
#Stream="SingleMuon"     
#Stream="DoubleMuon"     
#Stream="DoubleEG"
runFake="False"
#Skim="SKTree_LeptonSkim"  ### SKTree_NoSkim/SKTree_LeptonSkim/SKTree_Di[Tri]LepSkim/ flatcat
#Skim="SKTree_DiLepSkim"  ### SKTree_NoSkim/SKTree_LeptonSkim/SKTree_Di[Tri]LepSkim/ flatcat
Skim="SKTree_TriLepSkim"  ### SKTree_NoSkim/SKTree_LeptonSkim/SKTree_Di[Tri]LepSkim/ flatcat
DataPeriod="ALL"
#job_logstep=1000
LogLevel="INFO"
#QueueOption="longq" 
QueueOption="fastq" 


MCList="Analysis_bkg"
#MCList="ZG2l"
#MCList="SignalMajor_1e2mu"#MCList="SignalMajor_3mu"
#MCList="SignalMajor_3mu"
#MCList="DY"
#MCList="CR_EMu_804"
#MCList="Analysis_bkg_test"
###Backgound : AllSample / Analysis_bkg / Analysis_bkg_test / QCD_mu
###Signal    : Analysis_sig_All / Analysis_sig_1e2mu / Analysis_sig_3mu / tthwA_1e2mu / tthwA_3mu / Analysis_sig_test / Analysis_sig_test1

########################################################################
## OUTPUT PATH CONFIG

CATOUT='/data2/Users/jbhyun/cmssw/CATanalyzer/CATOut/'
outputdir_cat="/data2/Users/jbhyun/cmssw/CATanalyzer/CATOut/${CATVERSION}"

if [[ -z ${CATVERSION} ]]; then echo "source setup.sh needed. exit"; exit 1; fi
if [[ ! -d $CATOUT ]]; then echo "Made $CATOUT"; mkdir $CATOut; fi
if [[ ! -d ${outputdir_cat} ]]; then echo "Made ${outputdir_cat}"; mkdir ${outputdir_cat}; fi
if [[ ! -d ${outputdir_cat}/${AnalysisCode} ]]; then echo "Made ${outputdir_cat}/${AnalysisCode}"; mkdir ${outputdir_cat}/${AnalysisCode}; fi;

if [[ ${Stream} == 'DoubleMuon' || ${Stream} == 'SingleMuon' ]]; then outputdir_lep="${outputdir_cat}/${AnalysisCode}/Muon";
elif [[ ${Stream} == 'DoubleEG' || ${Stream} == 'SingleElectron' ]]; then outputdir_lep="${outputdir_cat}/${AnalysisCode}/Electron";
elif [[ ${Stream} == 'MuonEG' ]]; then outputdir_lep="${outputdir_cat}/${AnalysisCode}/ElectronMuon";
else echo "Error : job stream wrong. Check."; exit 1;
fi

if [[ ${DataPeriod}   == "A" ]]; then dir_period="periodA";
elif [[ ${DataPeriod} == "B" ]]; then dir_period="periodB";
elif [[ ${DataPeriod} == "C" ]]; then dir_period="periodC";
elif [[ ${DataPeriod} == "D" ]]; then dir_period="periodD";
elif [[ ${DataPeriod} == "E" ]]; then dir_period="periodE";
elif [[ ${DataPeriod} == "F" ]]; then dir_period="periodF";
elif [[ ${DataPeriod} == "G" ]]; then dir_period="periodG";
elif [[ ${DataPeriod} == "H" ]]; then dir_period="periodH";
elif [[ ${DataPeriod}  == "ALL" || ${DataPeriod} == "BtoH" ]]; then dir_period="periodBtoH";
else echo "Error: Period Set Wrongly"; exit 1;
fi
outputdir_period="${outputdir_lep}/${dir_period}/"
outputdir_fake="${outputdir_period}Fakes/";
OutputDir=${outputdir_period}


if [[ ! -d "${outputdir_lep}" ]]; then mkdir ${outputdir_lep}; echo "Made ${outputdir_lep}"; fi
if [[ ! -d "${outputdir_period}" ]]; then mkdir ${outputdir_period}; echo "Made ${outputdir_period}"; fi
if [[ ${runFake} == "True" || ${runFake} == "true" ]];
then
  if [[ ! -d "${outputdir_fake}" ]]; then mkdir ${outputdir_fake}; echo "Made ${outputdir_fake}"; fi
  OutputDir=${outputdir_fake}
fi


if [[ $runData == 'false' || $runData == 'False' ]]; then Stream=""; fi
if [[ $runMC == 'false' || $runMC == 'False' ]]; then MCList=""; fi

########################################################################
## COMMAND

if [[ ! -e CommandHist.txt ]]; then touch CommandHist.txt; echo "Command History">> CommandHist.txt; echo >> CommandHist.txt; fi;
date >> CommandHist.txt


if [[ $runData == 'false' || $runData == 'False' ]];
  then nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -a ${AnalysisCode} -p ${DataPeriod} -s ${Skim} -list ${MCList} -o ${OutputDir} -d ${LogLevel} -q ${QueueOption} 
  echo "nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -a ${AnalysisCode} -p ${DataPeriod} -s ${Skim} -list ${MCList} -o ${OutputDir} -d ${LogLevel} -q ${QueueOption}">> CommandHist.txt
elif [[ $runMC == 'false' || $runMC == 'False' ]];
  then nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -a ${AnalysisCode} -S ${Stream} -p ${DataPeriod} -s ${Skim} -o ${OutputDir} -d ${LogLevel} -q ${QueueOption} -fake ${runFake}
  echo "nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -a ${AnalysisCode} -S ${Stream} -p ${DataPeriod} -s ${Skim} -o ${OutputDir} -d ${LogLevel} -q ${QueueOption} -fake ${runFake}" >> CommandHist.txt
else
  nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -a ${AnalysisCode} -S ${Stream} -p ${DataPeriod} -s ${Skim} -list ${MCList} -o ${OutputDir} -d ${LogLevel} -q ${QueueOption}
  echo "nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -a ${AnalysisCode} -S ${Stream} -p ${DataPeriod} -s ${Skim} -list ${MCList} -o ${OutputDir} -d ${LogLevel} -q ${QueueOption}" >> CommandHist.txt
fi

echo >> CommandHist.txt
########################################################################
