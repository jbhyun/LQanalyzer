#!/bin/sh

###### SET WHAT JOBS TO RUN
runMC=true
runData=true


## RUN PARAMETERS
job_cycle="Jul2016_3l4j_DataValid" ###KPS ###"ExampleAnalyzerDiMuon"
job_data_lumi="ALL"   ###  "C" = period C only   "ALL" or "CtoD"  = period C+D
job_stream="SingleElectron" ###"SingleElectron" ### DoubleMuon DoubleEG MuonEG SingleElectron SingleMuon 
job_skinput="True" 

job_useskim="Lepton" ###"Lepton" ### "Lepton" for single lepton skim   "DiLep" for dilepton skim
job_logstep=1000
job_loglevel="DEBUG"
job_njobs=15

version_tag=${CATVERSION}

########################################################################
outputdir_cat=$LQANALYZER_DIR"/data/output/CAT/"

if [[ ${job_stream} == 'DoubleMuon' || ${job_stream} == 'SingleMuon' ]]; then outputdir_lep=$LQANALYZER_DIR"/data/output/CAT/Muon/";
elif [[ ${job_stream} == 'DoubleEG' || ${job_stream} == 'SingleElectron' ]]; then outputdir_lep=$LQANALYZER_DIR"/data/output/CAT/Electron/";
elif [[ ${job_stream} == 'MuonEG' ]]; then outputdir_lep=$LQANALYZER_DIR"/data/output/CAT/ElectronMuon/";
else echo "Error : job stream wrong. Check."; exit 1;
fi

if [[ ! -d "${outputdir_cat}" ]]; then
    mkdir ${outputdir_cat}
    echo "Making ${outputdir_cat}"
fi

if [[ ! -d "${outputdir_lep}" ]]; then
    mkdir ${outputdir_lep}
    echo "Making ${outputdir_lep}"
fi

if [[ $job_data_lumi  == "C" ]]; then dir_tag="periodC/";
elif [[ $job_data_lumi == "D" ]]; then dir_tag="periodD/";
elif [[ $job_data_lumi  == "ALL" ]]; then dir_tag="periodCtoD/";
else echo "Error: Period Set Wrongly"; exit 1;
fi

outputdir_mc=${outputdir_lep}${dir_tag}
outputdir_data=${outputdir_lep}${dir_tag}"/Data/"
echo ${outputdir_data}
output_datafile=${outputdir_data}'/'${job_cycle}"_data_cat"${version_tag}".root"
echo $output_datafile


if [[ ! -d "${outputdir_mc}" ]]; then
    mkdir ${outputdir_mc}
    echo "Making ${outputdir_mc}"
fi

if [[ ! -d "${outputdir_data}" ]]; then
    mkdir ${outputdir_data}
    echo "Making ${outputdir_data}"
fi
 


if [[ $1  == "MC"  ]];
then
    runMC=true
    runData=false
fi

if [[ $1  == "DATA"  ]];
then
    runMC=false
    runData=true
fi

########################################################################


if [[ $runMC  == "true" ]]; 
then
    source functions.sh
    
    cycle=$job_cycle
    data_lumi=$job_data_lumi
    stream=$job_stream
    useskim=$job_useskim
    skinput=${job_skinput}
    njobs=$job_njobs
    outputdir=${outputdir_mc}

###    declare -a input_samples=( 'DY10to50_MCatNLO' 'DY50plus_MCatNLO' 'WJets_MCatNLO' "TT_MCatNLO" 'TTG_MCatNLO' 'WZ_lllnu_MCatNLO' 'WZ_llqq_MCatNLO' 'ZZ_llqq_MCatNLO' 'ZZ_llll_MCatNLO' 'singletop_tbarW_Powheg' 'singletop_tW_Powheg' 'singletop_tbar_Powheg' 'singletop_t_Powheg' 'singletop_s_MCatNLO' 'ttHtobb_Powheg' 'ttHnobb_Powheg' 'ttWJetsToLNu_MCatNLO' 'ttZToQQ_MCatNLO' 'ttZToLLNuNu_MCatNLO' )
###    declare -a input_samples=('DY10to50_MCatNLO' 'DY50plus_MCatNLO' 'DY50plus_madgraph' 'GJet_20to40_pythia8' 'GJet_40plus_pythia8' 'ggHtomm_Powheg' 'QCD_mu1000toINF_pythia8' 'QCD_em120to170_pythia8' 'QCD_mu120to170_pythia8' 'QCD_em15to20_pythia8' 'QCD_em170to300_pythia8' 'QCD_mu170to300_pythia8' 'QCD_em20to30_pythia8' 'QCD_mu20to30_pythia8' 'QCD_mu300to470_pythia8' 'QCD_em300toINF_pythia8' 'QCD_DoubleEM_30to40_pythia8' 'QCD_mu30to50_pythia8' 'QCD_DoubleEM_30toInf_pythia8' 'QCD_DoubleEM_40toInf_pythia8' 'QCD_mu470to600_pythia8' 'QCD_em50to80_pythia8' 'QCD_mu50to80_pythia8' 'QCD_mu600to800_pythia8' 'QCD_mu800to1000_pythia8' 'QCD_em80to120_pythia8' 'QCD_mu80to120_pythia8' 'QCD_170to250_bcToE_pythia8' 'QCD_20to30_bcToE_pythia8' 'QCD_250toInf_bcToE_pythia8' 'QCD_30to80_bcToE_pythia8' 'QCD_80to170_bcToE_pythia8' 'singletop_s_MCatNLO' 'singletop_tbar_Powheg' 'singletop_t_Powheg' 'singletop_tbarW_Powheg' 'singletop_tW_Powheg' 'TTG_MCatNLO' 'TT_MCatNLO' 'TT_MG5' 'ttWJetsToLNu_MCatNLO' 'ttWJetsToQQ_MCatNLO' 'ttZToQQ_MCatNLO' 'TT_powheg' 'vhf_Htomm_Powheg' 'WG_lnuG_madgraph' 'WJets_MCatNLO' 'WW_llnn_powheg' 'WW_doublescattering' 'WW_pythia8' 'WZ_lllnu_powheg' 'WZ_lllnu_MCatNLO' 'WZ_llqq_MCatNLO' 'WN_lllnu_powheg' 'WZZ_MCatNLO' 'WZ_pythia8' 'WpWp_madgraph' 'WpWp_qcd_madgraph' 'ZG_llG_MCatNLO' 'ZZ_llnunu_powheg' 'ZZ_llqq_MCatNLO' 'ZZ_llll_MCatNLO' 'ZZ_llll_powheg' 'ZZ_pythia8' 'ttHnobb_Powheg' 'ttHtobb_Powheg')
###    declare -a input_samples=("tthwA_3l4j_hc130A5_3mu" "tthwA_3l4j_hc130A5_emu")
###    declare -a input_samples=("tthwA_3l4j_hc130A30_emu_catcut")
###declare -a input_samples=("ZZ_pythia8")
#    declare -a input_samples=("QCD_mu20to30_pythia8" "QCD_mu30to50_pythia8" "QCD_mu50to80_pythia8" "QCD_mu80to120_pythia8" "QCD_mu120to170_pythia8" "QCD_mu170to300_pythia8" "QCD_mu300to470_pythia8" "QCD_mu470to600_pythia8" "QCD_mu600to800_pythia8" "QCD_mu800to1000_pythia8" "QCD_mu1000toINF_pythia8")
###    declare -a input_samples=()

    declare -a input_samples=("WJets_MCatNLO" "WG_lnuG_madgraph" "TT_MG5" "ZG_llG_MCatNLO" "WW_pythia8" "DY10to50_MCatNLO" "DY50plus_MCatNLO" 'QCD_em15to20_pythia8' 'QCD_em20to30_pythia8' 'QCD_em50to80_pythia8' 'QCD_em80to120_pythia8' 'QCD_em120to170_pythia8' 'QCD_em170to300_pythia8' 'QCD_20to30_bcToE_pythia8' 'QCD_30to80_bcToE_pythia8' 'QCD_80to170_bcToE_pythia8' 'QCD_170to250_bcToE_pythia8' 'QCD_250toInf_bcToE_pythia8')
###    declare -a input_samples=("DY50plus_madgraph")
###    declare -a input_samples=('DY10to50_MCatNLO' 'DY50plus_MCatNLO' "QCD_mu20to30_pythia8" "QCD_mu30to50_pythia8" "QCD_mu50to80_pythia8" "QCD_mu80to120_pythia8" "QCD_mu120to170_pythia8" "QCD_mu170to300_pythia8" "QCD_mu300to470_pythia8" "QCD_mu470to600_pythia8" "QCD_mu600to800_pythia8" "QCD_mu800to1000_pythia8" "QCD_mu1000toINF_pythia8")
    source submit.sh $1
fi


################ DOUBLEMUON DATA
if [[ $runData  == "true" ]];
then
    source functions.sh

    cycle=$job_cycle
    data_lumi=$job_data_lumi
    stream=$job_stream
    useskim=$job_useskim
    skinput=${job_skinput}
    njobs=$job_njobs
    outputdir=${outputdir_data}
    
    if [[ $job_data_lumi  == "ALL" ]]; then declare -a input_samples=("C" "D");    
    elif [[ $job_data_lumi  == "C" ]]; then declare -a input_samples=("C");
    elif [[ $job_data_lumi == 'D' ]]; then declare -a input_samples=('D');
    else echo "Data period set wrong. check."; exit 1;
    fi
    
    source submit.sh $1
    
fi




echo ""
echo "End of example_submit.sh script."
