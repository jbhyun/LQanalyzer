#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################

declare -a example=('ttZ' )
declare -a AllSample=('WJets' 'DYJets_10to50' 'DYJets' 'SingleTop_s' 'SingleTbar_t' 'SingleTop_t' 'SingleTbar_tW' 'SingleTop_tW' 'TT_powheg' 'ZZ' 'WZ' 'WW' 'WGtoLNuG' 'WGtoLNuEE' 'WGtoLNuMM' 'ZGto2LG' 'ttH_nonbb' 'ttH_bb' 'ttW' 'ttZ') 



#####################################################################################
#Background#################
#####################################################################################

#Single Process
declare -a ST=('SingleTop_s' 'SingleTbar_t' 'SingleTop_t' 'SingleTbar_tW_noHadron' 'SingleTop_tW_noHadron')
declare -a DY=('DYJets_10to50' "DYJets") 
declare -a DYMG=("DYJets_MG")
declare -a ZZ4l=("ZZTo4L_powheg") 
declare -a ZG2l=("ZGto2LG" "WGtoLNuG") 
declare -a TT=('TT_powheg') #"TTJets_aMC") 

#Analysis Background
declare -a Analysis_bkg=('WZTo3LNu_powheg' 'ZZTo4L_powheg' "ggZZto2e2mu" "ggZZto2e2tau" "ggZZto2mu2tau" "ggZZto4e" "ggZZto4mu" "ggZZto4tau" "ZGto2LG" 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZ' 'tZq' 'ttH_nonbb' 'WWW' 'WWZ' 'WZZ' 'ZZZ' "vbhHtoZZ" "ggHtoZZ")


#Dilepton Validation
declare -a CR_DiLep=('DYJets_10to50' 'DYJets' 'SingleTbar_tW_noHadron' 'SingleTop_tW_noHadron' 'TT_powheg' 'WZ' 'ZZ' 'WWTo2L2Nu' 'ZGto2LG')
declare -a CR_EMu=('WJets' 'DYJets_10to50' 'DYJets' 'SingleTbar_tW_noHadron' 'SingleTop_tW_noHadron' 'WZ' 'ZZ' 'WWTo2L2Nu' 'ZGto2LG' "WGtoLNuG")
#declare -a CR_EMu=('DYJets_10to50' 'DYJets' 'SingleTbar_tW_noHadron' 'SingleTop_tW_noHadron' 'TT_powheg' 'WZ' 'ZZ' 'WWTo2L2Nu' 'ZGto2LG')
declare -a CR_MuMu_fast=('DYJets_10to50' 'DYJets')
declare -a CR_EE_fast=('DYJets_10to50' 'DYJets')

#4lepton CR
#declare -a CR_4lep=("ZZTo4L_powheg" "ggZZto2e2mu" "ggZZto2e2tau" "ggZZto2mu2tau" "ggZZto4e" "ggZZto4mu" "ggZZto4tau" "ttZToLL_M-1to10" "ttZToLL_M-10" "vbhHtoZZ" "ggHtoZZ" "ttH_nonbb" "ZZZ" "WZZ" "WWZ")
declare -a CR_4lep=("ZZTo4L_powheg" "ggZZto2e2mu" "ggZZto2e2tau" "ggZZto2mu2tau" "ggZZto4e" "ggZZto4mu" "ggZZto4tau" "ttZToLL_M-1to10" "ttZ" "vbhHtoZZ" "ggHtoZZ" "ttH_nonbb" "ZZZ" "WZZ" "WWZ" "WZG" "ZGto2LG")
declare -a CR_4lep_ZGtestAdd=("WZG" "ZGto2LG")
declare -a WWG=("WWG")
declare -a ggH=("ggHtoZZ")

#Btag Efficiency Measurement
declare -a BtagEffSample=('TT_powheg')
#declare -a BtagEffSample=('TTJets_aMC' "DYJets")

#ID, Trigger Efficiency Measurement
declare -a IDSample=('DYJets_MG' )
declare -a TrigSample=('TT_powheg' 'DYJets_MG' 'DYJets_10to50' )
declare -a TrigDiLepClosure=('TT_powheg' 'DYJets_MG' 'DYJets_10to50' )
declare -a TrigTriLepClosure=('WZTo3LNu_powheg' 'ttZ' 'ZZTo4L_powheg' )
declare -a FR_Prompt=('WJets' 'DYJets_10to50' 'DYJets' 'TT_powheg')
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
