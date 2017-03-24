#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################

declare -a example=("tthwA_1ta2mu_hc130A30" )
declare -a AllSample=('WJets' 'DYJets_10to50' 'DYJets' 'SingleTop_s' 'SingleTbar_t' 'SingleTop_t' 'SingleTbar_tW' 'SingleTop_tW' 'TT_powheg' 'ZZ' 'WZ' 'WW' 'WGtoLNuG' 'WGtoLNuEE' 'WGtoLNuMM' 'ZGto2LG' 'ttH_nonbb' 'ttH_bb' 'ttW' 'ttZ') 


#####################################################################################
#Background#################
#####################################################################################

#Single Process
declare -a ST=('SingleTop_s' 'SingleTbar_t' 'SingleTop_t' 'SingleTbar_tW_noHadron' 'SingleTop_tW_noHadron')
declare -a DY=('DYJets_10to50' "DYJets") 
declare -a ZZ4l=("ZZTo4L_powheg") 
declare -a ZG2l=("ZGto2LG" "WGtoLNuG") 
declare -a TT=('TT_powheg') #"TTJets_aMC") 

#Analysis Background
declare -a Analysis_bkg=('WZTo3LNu_amcatnlo' 'ZZTo4L_powheg' 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZ' 'ttH_nonbb' 'WWW' 'WWZ' 'WZZ' 'ZZZ')

#Dilepton Validation
declare -a CR_DiLep=('DYJets_10to50' 'DYJets' 'SingleTbar_tW_noHadron' 'SingleTop_tW_noHadron' 'TT_powheg' 'WZ' 'ZZ' 'WWTo2L2Nu' 'ZGto2LG')
declare -a CR_EMu=('DYJets_10to50' 'DYJets' 'SingleTbar_tW_noHadron' 'SingleTop_tW_noHadron' 'TT_powheg' 'WZ' 'ZZ' 'WWTo2L2Nu' 'ZGto2LG')
declare -a CR_MuMu_fast=('DYJets_10to50' 'DYJets')
declare -a CR_EE_fast=('DYJets_10to50' 'DYJets')

#Btag Efficiency Measurement
declare -a BtagEffSample=('TT_powheg')
#declare -a BtagEffSample=('TTJets_aMC' "DYJets")

#####################################################################################


#####################################################################################
##Signal###############
#####################################################################################
declare -a SignalMajor_All=("TTToHcToWA_1e2mu_MHc100_MA15" "TTToHcToWA_1e2mu_MHc110_MA30" "TTToHcToWA_1e2mu_MHc160_MA15" "TTToHcToWA_1e2mu_MHc160_MA30" "TTToHcToWZp_1e2mu_MHc160_MZp5" "TTToHcToWZp_1e2mu_MHc90_MZp5" "TTToHcToWA_3mu_MHc100_MA15" "TTToHcToWA_3mu_MHc110_MA30" "TTToHcToWA_3mu_MHc160_MA15" "TTToHcToWA_3mu_MHc160_MA30" "TTToHcToWZp_3mu_MHc160_MZp5" "TTToHcToWZp_3mu_MHc90_MZp5" "TTToHcToWA_1ta2mu_MHc100_MA15" "TTToHcToWA_1ta2mu_MHc110_MA30" "TTToHcToWA_1ta2mu_MHc160_MA15" "TTToHcToWA_1ta2mu_MHc160_MA30" "TTToHcToWZp_1ta2mu_MHc160_MZp5" "TTToHcToWZp_1ta2mu_MHc90_MZp5" "TTToHcToWA_2l2mu_MHc100_MA15" "TTToHcToWA_2l2mu_MHc110_MA30" "TTToHcToWA_2l2mu_MHc160_MA15" "TTToHcToWA_2l2mu_MHc160_MA30" "TTToHcToWZp_2l2mu_MHc160_MZp5" "TTToHcToWZp_2l2mu_MHc90_MZp5")


declare -a SignalMajor_1e2mu=("TTToHcToWA_1e2mu_MHc100_MA15" "TTToHcToWA_1e2mu_MHc110_MA30" "TTToHcToWA_1e2mu_MHc160_MA15" "TTToHcToWA_1e2mu_MHc160_MA30" "TTToHcToWZp_1e2mu_MHc160_MZp5" "TTToHcToWZp_1e2mu_MHc90_MZp5") 
declare -a SignalMajor_3mu=("TTToHcToWA_3mu_MHc100_MA15" "TTToHcToWA_3mu_MHc110_MA30" "TTToHcToWA_3mu_MHc160_MA15" "TTToHcToWA_3mu_MHc160_MA30" "TTToHcToWZp_3mu_MHc160_MZp5" "TTToHcToWZp_3mu_MHc90_MZp5")
declare -a SignalMajor_1tamu=("TTToHcToWA_1ta2mu_MHc100_MA15" "TTToHcToWA_1ta2mu_MHc110_MA30" "TTToHcToWA_1ta2mu_MHc160_MA15" "TTToHcToWA_1ta2mu_MHc160_MA30" "TTToHcToWZp_1ta2mu_MHc160_MZp5" "TTToHcToWZp_1ta2mu_MHc90_MZp5") 
declare -a SignalMajor_2l2mu=("TTToHcToWA_2l2mu_MHc100_MA15" "TTToHcToWA_2l2mu_MHc110_MA30" "TTToHcToWA_2l2mu_MHc160_MA15" "TTToHcToWA_2l2mu_MHc160_MA30" "TTToHcToWZp_2l2mu_MHc160_MZp5" "TTToHcToWZp_2l2mu_MHc90_MZp5")
#####################################################################################


#####################################################################################
##Just in Case####
#####################################################################################

#DataValidation-SingleE
declare -a Validation_SEle=("WJets_MCatNLO" "WG_lnuG_madgraph" "TT_MG5" "ZG_llG_MCatNLO" "WW_pythia8" "DY10to50_MCatNLO" "DY50plus_MCatNLO" 'QCD_em15to20_pythia8' 'QCD_em20to30_pythia8' 'QCD_em50to80_pythia8' 'QCD_em80to120_pythia8' 'QCD_em120to170_pythia8' 'QCD_em170to300_pythia8' 'QCD_20to30_bcToE_pythia8' 'QCD_30to80_bcToE_pythia8' 'QCD_80to170_bcToE_pythia8' 'QCD_170to250_bcToE_pythia8' 'QCD_250toInf_bcToE_pythia8') 

#MuMubj Test
declare -a MuMubj_bkg=("TT_powheg" "DY50plus_MCatNLO" "DY10to50_MCatNLO" "singletop_tbar_Powheg" "singletop_s_MCatNLO" "singletop_tW_Powheg" "singletop_tbarW_Powheg" "singletop_t_Powheg" "TTG_MCatNLO" "ttZToLLNuNu_MCatNLO" "ttWJetsToLNu_MCatNLO" "ttZToQQ_MCatNLO" "ttHnobb_Powheg" "WZ_lllnu_MCatNLO" "WZ_llqq_MCatNLO" "WW_pythia8" "ZZ_llll_MCatNLO" "ZG_llG_MCatNLO" "ZZ_llqq_MCatNLO")
#####################################################################################
