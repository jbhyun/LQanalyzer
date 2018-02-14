#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################
declare -a tchannel_hn=('HNMoriondLL_Tchannel_EpEp_100' 'HNMoriondLL_Tchannel_EpEp_200' 'HNMoriondLL_Tchannel_EpEp_500' 'HNMoriondLL_Tchannel_EpEp_1100' 'HNMoriondLL_Tchannel_MupMup_100' 'HNMoriondLL_Tchannel_MupMup_200' 'HNMoriondLL_Tchannel_MupMup_500' 'HNMoriondLL_Tchannel_MupMup_1100' 'HNMumMum_40' 'HNMumMum_50' 'HNMumMum_200' 'HNMumMum_500' 'HNMumMum_1500'  'HNEmEm_40' 'HNEmEm_50' 'HNEmEm_200' 'HNEmEm_500' 'HNEmEm_1500'  'HNMupMup_40' 'HNMupMup_50' 'HNMupMup_200' 'HNMupMup_500' 'HNMupMup_1500' 'HNMupMup_100' 'HNMupMup_700' 'HNMupMup_1000' )
declare -a  mc_conv=('ZGto2LG' 'WGtoLNuG' 'TG' 'TTG' 'WZG' 'WWG' )

declare -a ListToProd=("DYJets_MG_10to50_TrigInfoSample" "DYJets_MG_TrigInfoSample" "TT_powheg_TrigInfoSample")

declare -a AllSample=('WJets' 'DYJets_10to50' 'DYJets' 'SingleTop_s' 'SingleTbar_t' 'SingleTop_t' 'SingleTbar_tW' 'SingleTop_tW' 'TT_powheg' 'ZZ' 'WZ' 'WW' 'WGtoLNuG' 'WGtoLNuEE' 'WGtoLNuMM' 'ZGto2LG' 'ttH_nonbb' 'ttH_bb' 'ttW' 'ttZ') 
declare -a UnRun=("TT_powheg_TrigInfoSample")


#####################################################################################
#Background#################
#####################################################################################

#Single Process
declare -a ST=('SingleTop_s' 'SingleTbar_t' 'SingleTop_t' 'SingleTbar_tW_noHadron' 'SingleTop_tW_noHadron')
declare -a DY=('DYJets_10to50' "DYJets") 
declare -a VV=('ZZ' 'WZ' 'WW')
declare -a DY50=("DYJets") 
declare -a DYMG=("DYJets_MG")
declare -a VJets=('DYJets_10to50' 'LowStat_DYJets' 'WJets')
declare -a ZZ4l=("ZZTo4L_powheg") 
declare -a ZG2l=("ZGto2LG") 
#declare -a TT=('TT_powheg' 'TTLL_powheg' 'TTLJ_powheg') #"TTJets_aMC") 
declare -a TTLJ=('TTLJ_powheg') 
declare -a TTLL=('TTLL_powheg') 
declare -a TT=('TT_powheg') #"TTJets_aMC") 
declare -a TTG=('TTG') #"TTJets_aMC") 

#Analysis Background
#declare -a Analysis_bkg=('WZTo3LNu_powheg' 'ZZTo4L_powheg' "ZGto2LG" "TTG" 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZ' 'tZq' 'ttH_nonbb' 'WWW' 'WWZ' 'WZZ' 'ZZZ' "vbhHtoZZ" "ggHtoZZ" "DYJets" "DYJets_10to50" "TTLL_powheg")
declare -a Analysis_bkg=('WZTo3LNu_powheg' 'ZZTo4L_powheg' "ZGto2LG" 'TTG' 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZ' 'tZq' 'ttH_nonbb' 'WWW' 'WWZ' 'WZZ' 'ZZZ' "vbhHtoZZ" "ggHtoZZ")
#declare -a Analysis_bkg=('WZTo3LNu_powheg' 'ZZTo4L_powheg' "ggZZto2e2mu" "ggZZto2e2tau" "ggZZto2mu2tau" "ggZZto4e" "ggZZto4mu" "ggZZto4tau" "ZGto2LG" 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZ' 'tZq' 'ttH_nonbb' 'WWW' 'WWZ' 'WZZ' 'ZZZ' "vbhHtoZZ" "ggHtoZZ")


#Dilepton Validation
declare -a CR_DiLep=('DYJets_10to50' 'DYJets' 'SingleTbar_tW_noHadron' 'SingleTop_tW_noHadron' 'TT_powheg' 'WZ' 'ZZ' 'WW')
declare -a CR_EMu=('DYJets_10to50' 'LowStat_DYJets' 'SingleTbar_tW_noHadron' 'SingleTop_tW_noHadron' 'TT_powheg' 'WZ' 'ZZ' 'WW')
declare -a CR_MuMu_fast=('DYJets_10to50' 'DYJets' 'TT_powheg')
declare -a CR_MuMu50_fast=('DYJets' 'TT_powheg')
declare -a CR_EE_fast=('DYJets_10to50' 'LowStat_DYJets' 'TT_powheg')

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
declare -a ObjEff=('TT_powheg' 'DYJets' 'ZGto2LG')
declare -a IDSample=('DYJets_MG')
declare -a TrigSample=("DYJets_MG_10to50_TrigInfoSample" "DYJets_MG_TrigInfoSample" "TT_powheg_TrigInfoSample")
declare -a TrigDiLepClosure=('TT_powheg' 'DYJets_MG' 'DYJets_10to50' )
declare -a TrigTriLepClosure=('WZTo3LNu_powheg' 'ttZ' 'ZZTo4L_powheg' )
declare -a FR_Prompt=('WJets' 'DYJets_10to50' 'DYJets' 'TT_powheg' 'SingleTop_t' 'SingleTop_s' 'SingleTbar_t' 'SingleTbar_tW_noHadron' 'SingleTop_tW_noHadron' 'WW' 'WZ' 'ZZ')
#declare -a FR_Prompt=('LowStat_WJets' 'DYJets_10to50' 'LowStat_DYJets' 'TT_powheg' 'SingleTop_t' 'SingleTop_s' 'SingleTbar_t' 'SingleTbar_tW_noHadron' 'SingleTop_tW_noHadron' 'WW' 'WZ' 'ZZ')
#declare -a MajorFakeSource=('DYJets_10to50' 'DYJets' 'TTLL_powheg' 'TT_powheg')
declare -a MajorFakeSource=('DYJets_10to50' 'DYJets' 'TT_powheg')
declare -a FakeMeasRegSample=("qcd_15to20_bctoe" "qcd_20to30_bctoe" "qcd_30to80_bctoe" "qcd_80to170_bctoe" "qcd_170to250_bctoe" "qcd_250toinf_bctoe" "QCD_Pt-20to30_EMEnriched" "QCD_Pt-30to50_EMEnriched" "QCD_Pt-50to80_EMEnriched" "QCD_Pt-80to120_EMEnriched" "QCD_Pt-120to170_EMEnriched" "QCD_Pt-170to300_EMEnriched" "QCD_Pt-300toInf_EMEnriched")
declare -a QCD_BCToE=("qcd_15to20_bctoe" "qcd_20to30_bctoe" "qcd_30to80_bctoe" "qcd_80to170_bctoe" "qcd_170to250_bctoe" "qcd_250toinf_bctoe")
declare -a QCD_EM=("QCD_Pt-20to30_EMEnriched" "QCD_Pt-30to50_EMEnriched" "QCD_Pt-50to80_EMEnriched" "QCD_Pt-80to120_EMEnriched" "QCD_Pt-120to170_EMEnriched" "QCD_Pt-170to300_EMEnriched" "QCD_Pt-300toInf_EMEnriched")
declare -a QCD_DoubleEM=("QCD_DoubleEMEnriched_30-40_mgg80toinf" "QCD_DoubleEMEnriched_30-inf_mgg40to80" "QCD_DoubleEMEnriched_40-inf_mgg80toinf") 
declare -a QCD_Mu=("QCD_Pt-15to20_MuEnriched" "QCD_Pt-20to30_MuEnriched" "QCD_Pt-30to50_MuEnriched" "QCD_Pt-50to80_MuEnriched" "QCD_Pt-80to120_MuEnriched" "QCD_Pt-120to170_MuEnriched" "QCD_Pt-170to300_MuEnriched" "QCD_Pt-300to470_MuEnriched" "QCD_Pt-470to600_MuEnriched" "QCD_Pt-600to800_MuEnriched" "QCD_Pt-800to1000_MuEnriched" "QCD_Pt-1000toInf_MuEnriched")

#GenSystSample
declare -a GenSystAnaTargetBkg=('ttZ' 'ttZToLL_M-1to10' 'ttH_nonbb' 'ttWToLNu' 'ttZToLL_M-10' 'WZTo3LNu_powheg' 'ZZTo4L_powheg')
declare -a TopPtMeasSample_SM=('TT_powheg')
declare -a TopPtMeasSample_Sig=('TTToHcToWA_incl_MHc160_MA35')
#####################################################################################


#####################################################################################
##Signal###############
#####################################################################################


declare -a SignalMajor_All=("TTToHcToWA_incl_MHc100_MA15" "TTToHcToWA_incl_MHc120_MA15" "TTToHcToWA_incl_MHc120_MA25" "TTToHcToWA_incl_MHc120_MA35" "TTToHcToWA_incl_MHc140_MA15" "TTToHcToWA_incl_MHc140_MA25" "TTToHcToWA_incl_MHc140_MA35" "TTToHcToWA_incl_MHc160_MA15" "TTToHcToWA_incl_MHc160_MA25" "TTToHcToWA_incl_MHc160_MA35")

declare -a SignalMajor_Old=("TTToHcToWA_1e2mu_MHc100_MA15" "TTToHcToWA_1e2mu_MHc120_MA15" "TTToHcToWA_1e2mu_MHc120_MA25" "TTToHcToWA_1e2mu_MHc120_MA35" "TTToHcToWA_1e2mu_MHc140_MA15" "TTToHcToWA_1e2mu_MHc140_MA25" "TTToHcToWA_1e2mu_MHc140_MA35" "TTToHcToWA_1e2mu_MHc160_MA15" "TTToHcToWA_1e2mu_MHc160_MA25" "TTToHcToWA_1e2mu_MHc160_MA35" "TTToHcToWA_3mu_MHc100_MA15" "TTToHcToWA_3mu_MHc120_MA15" "TTToHcToWA_3mu_MHc120_MA25" "TTToHcToWA_3mu_MHc120_MA35" "TTToHcToWA_3mu_MHc140_MA15" "TTToHcToWA_3mu_MHc140_MA25" "TTToHcToWA_3mu_MHc140_MA35" "TTToHcToWA_3mu_MHc160_MA15" "TTToHcToWA_3mu_MHc160_MA25" "TTToHcToWA_3mu_MHc160_MA35" "TTToHcToWA_1ta2mu_MHc100_MA15" "TTToHcToWA_1ta2mu_MHc120_MA15" "TTToHcToWA_1ta2mu_MHc120_MA25" "TTToHcToWA_1ta2mu_MHc120_MA35" "TTToHcToWA_1ta2mu_MHc140_MA15" "TTToHcToWA_1ta2mu_MHc140_MA25" "TTToHcToWA_1ta2mu_MHc140_MA35" "TTToHcToWA_1ta2mu_MHc160_MA15" "TTToHcToWA_1ta2mu_MHc160_MA25" "TTToHcToWA_1ta2mu_MHc160_MA35" "TTToHcToWA_2l2mu_MHc100_MA15" "TTToHcToWA_2l2mu_MHc120_MA15" "TTToHcToWA_2l2mu_MHc120_MA25" "TTToHcToWA_2l2mu_MHc120_MA35" "TTToHcToWA_2l2mu_MHc140_MA15" "TTToHcToWA_2l2mu_MHc140_MA25" "TTToHcToWA_2l2mu_MHc140_MA35" "TTToHcToWA_2l2mu_MHc160_MA15" "TTToHcToWA_2l2mu_MHc160_MA25" "TTToHcToWA_2l2mu_MHc160_MA35")

#declare -a SignalMajor_1e2mu=("TTToHcToWA_1e2mu_MHc100_MA15" "TTToHcToWA_1e2mu_MHc110_MA30" "TTToHcToWA_1e2mu_MHc160_MA15" "TTToHcToWA_1e2mu_MHc160_MA30") 
#declare -a SignalMajor_3mu=("TTToHcToWA_3mu_MHc100_MA15" "TTToHcToWA_3mu_MHc110_MA30" "TTToHcToWA_3mu_MHc160_MA15" "TTToHcToWA_3mu_MHc160_MA30")
#declare -a SignalMajor_1tamu=("TTToHcToWA_1ta2mu_MHc100_MA15" "TTToHcToWA_1ta2mu_MHc110_MA30" "TTToHcToWA_1ta2mu_MHc160_MA15" "TTToHcToWA_1ta2mu_MHc160_MA30") 
#declare -a SignalMajor_2l2mu=("TTToHcToWA_2l2mu_MHc100_MA15" "TTToHcToWA_2l2mu_MHc110_MA30" "TTToHcToWA_2l2mu_MHc160_MA15" "TTToHcToWA_2l2mu_MHc160_MA30")
#####################################################################################


#####################################################################################
##Just in Case####
#####################################################################################

#DataValidation-SingleE
declare -a Validation_SEle=("WJets_MCatNLO" "WG_lnuG_madgraph" "TT_MG5" "ZG_llG_MCatNLO" "WW_pythia8" "DY10to50_MCatNLO" "DY50plus_MCatNLO" 'QCD_em15to20_pythia8' 'QCD_em20to30_pythia8' 'QCD_em50to80_pythia8' 'QCD_em80to120_pythia8' 'QCD_em120to170_pythia8' 'QCD_em170to300_pythia8' 'QCD_20to30_bcToE_pythia8' 'QCD_30to80_bcToE_pythia8' 'QCD_80to170_bcToE_pythia8' 'QCD_170to250_bcToE_pythia8' 'QCD_250toInf_bcToE_pythia8') 

#MuMubj Test
declare -a MuMubj_bkg=("TT_powheg" "DY50plus_MCatNLO" "DY10to50_MCatNLO" "singletop_tbar_Powheg" "singletop_s_MCatNLO" "singletop_tW_Powheg" "singletop_tbarW_Powheg" "singletop_t_Powheg" "TTG_MCatNLO" "ttZToLLNuNu_MCatNLO" "ttWJetsToLNu_MCatNLO" "ttZToQQ_MCatNLO" "ttHnobb_Powheg" "WZ_lllnu_MCatNLO" "WZ_llqq_MCatNLO" "WW_pythia8" "ZZ_llll_MCatNLO" "ZG_llG_MCatNLO" "ZZ_llqq_MCatNLO")
#####################################################################################

#LONG LIST###########################################################################
#####################################################################################
declare -a Signal_1e2mu=("TTToHcToWA_1e2mu_MHc100_MA15" "TTToHcToWA_1e2mu_MHc100_MA20" "TTToHcToWA_1e2mu_MHc110_MA15" "TTToHcToWA_1e2mu_MHc110_MA20" "TTToHcToWA_1e2mu_MHc110_MA25" "TTToHcToWA_1e2mu_MHc110_MA30" "TTToHcToWA_1e2mu_MHc120_MA15" "TTToHcToWA_1e2mu_MHc120_MA20" "TTToHcToWA_1e2mu_MHc120_MA25" "TTToHcToWA_1e2mu_MHc120_MA30" "TTToHcToWA_1e2mu_MHc120_MA35" "TTToHcToWA_1e2mu_MHc130_MA15" "TTToHcToWA_1e2mu_MHc130_MA20" "TTToHcToWA_1e2mu_MHc130_MA25" "TTToHcToWA_1e2mu_MHc130_MA30" "TTToHcToWA_1e2mu_MHc130_MA35" "TTToHcToWA_1e2mu_MHc140_MA15" "TTToHcToWA_1e2mu_MHc140_MA20" "TTToHcToWA_1e2mu_MHc140_MA25" "TTToHcToWA_1e2mu_MHc140_MA30" "TTToHcToWA_1e2mu_MHc140_MA35" "TTToHcToWA_1e2mu_MHc150_MA15" "TTToHcToWA_1e2mu_MHc150_MA20" "TTToHcToWA_1e2mu_MHc150_MA25" "TTToHcToWA_1e2mu_MHc150_MA30" "TTToHcToWA_1e2mu_MHc150_MA35" "TTToHcToWA_1e2mu_MHc160_MA15" "TTToHcToWA_1e2mu_MHc160_MA20" "TTToHcToWA_1e2mu_MHc160_MA25" "TTToHcToWA_1e2mu_MHc160_MA30" "TTToHcToWA_1e2mu_MHc160_MA35")
declare -a Signal_2l2mu=("TTToHcToWA_2l2mu_MHc100_MA15" "TTToHcToWA_2l2mu_MHc100_MA20" "TTToHcToWA_2l2mu_MHc110_MA15" "TTToHcToWA_2l2mu_MHc110_MA20" "TTToHcToWA_2l2mu_MHc110_MA25" "TTToHcToWA_2l2mu_MHc110_MA30" "TTToHcToWA_2l2mu_MHc120_MA15" "TTToHcToWA_2l2mu_MHc120_MA20" "TTToHcToWA_2l2mu_MHc120_MA25" "TTToHcToWA_2l2mu_MHc120_MA30" "TTToHcToWA_2l2mu_MHc120_MA35" "TTToHcToWA_2l2mu_MHc130_MA15" "TTToHcToWA_2l2mu_MHc130_MA20" "TTToHcToWA_2l2mu_MHc130_MA25" "TTToHcToWA_2l2mu_MHc130_MA30" "TTToHcToWA_2l2mu_MHc130_MA35" "TTToHcToWA_2l2mu_MHc140_MA15" "TTToHcToWA_2l2mu_MHc140_MA20" "TTToHcToWA_2l2mu_MHc140_MA25" "TTToHcToWA_2l2mu_MHc140_MA30" "TTToHcToWA_2l2mu_MHc140_MA35" "TTToHcToWA_2l2mu_MHc150_MA15" "TTToHcToWA_2l2mu_MHc150_MA20" "TTToHcToWA_2l2mu_MHc150_MA25" "TTToHcToWA_2l2mu_MHc150_MA30" "TTToHcToWA_2l2mu_MHc150_MA35" "TTToHcToWA_2l2mu_MHc160_MA15" "TTToHcToWA_2l2mu_MHc160_MA20" "TTToHcToWA_2l2mu_MHc160_MA25" "TTToHcToWA_2l2mu_MHc160_MA30" "TTToHcToWA_2l2mu_MHc160_MA35")
declare -a Signal_3mu=("TTToHcToWA_3mu_MHc100_MA15" "TTToHcToWA_3mu_MHc100_MA20" "TTToHcToWA_3mu_MHc110_MA15" "TTToHcToWA_3mu_MHc110_MA20" "TTToHcToWA_3mu_MHc110_MA25" "TTToHcToWA_3mu_MHc110_MA30" "TTToHcToWA_3mu_MHc120_MA15" "TTToHcToWA_3mu_MHc120_MA20" "TTToHcToWA_3mu_MHc120_MA25" "TTToHcToWA_3mu_MHc120_MA30" "TTToHcToWA_3mu_MHc120_MA35" "TTToHcToWA_3mu_MHc130_MA15" "TTToHcToWA_3mu_MHc130_MA20" "TTToHcToWA_3mu_MHc130_MA25" "TTToHcToWA_3mu_MHc130_MA30" "TTToHcToWA_3mu_MHc130_MA35" "TTToHcToWA_3mu_MHc140_MA15" "TTToHcToWA_3mu_MHc140_MA20" "TTToHcToWA_3mu_MHc140_MA25" "TTToHcToWA_3mu_MHc140_MA30" "TTToHcToWA_3mu_MHc140_MA35" "TTToHcToWA_3mu_MHc150_MA15" "TTToHcToWA_3mu_MHc150_MA20" "TTToHcToWA_3mu_MHc150_MA25" "TTToHcToWA_3mu_MHc150_MA30" "TTToHcToWA_3mu_MHc150_MA35" "TTToHcToWA_3mu_MHc160_MA15" "TTToHcToWA_3mu_MHc160_MA20" "TTToHcToWA_3mu_MHc160_MA25" "TTToHcToWA_3mu_MHc160_MA30" "TTToHcToWA_3mu_MHc160_MA35")
declare -a Signal_1ta2mu=("TTToHcToWA_1ta2mu_MHc100_MA15" "TTToHcToWA_1ta2mu_MHc100_MA20" "TTToHcToWA_1ta2mu_MHc110_MA15" "TTToHcToWA_1ta2mu_MHc110_MA20" "TTToHcToWA_1ta2mu_MHc110_MA25" "TTToHcToWA_1ta2mu_MHc110_MA30" "TTToHcToWA_1ta2mu_MHc120_MA15" "TTToHcToWA_1ta2mu_MHc120_MA20" "TTToHcToWA_1ta2mu_MHc120_MA25" "TTToHcToWA_1ta2mu_MHc120_MA30" "TTToHcToWA_1ta2mu_MHc120_MA35" "TTToHcToWA_1ta2mu_MHc130_MA15" "TTToHcToWA_1ta2mu_MHc130_MA20" "TTToHcToWA_1ta2mu_MHc130_MA25" "TTToHcToWA_1ta2mu_MHc130_MA30" "TTToHcToWA_1ta2mu_MHc130_MA35" "TTToHcToWA_1ta2mu_MHc140_MA15" "TTToHcToWA_1ta2mu_MHc140_MA20" "TTToHcToWA_1ta2mu_MHc140_MA25" "TTToHcToWA_1ta2mu_MHc140_MA30" "TTToHcToWA_1ta2mu_MHc140_MA35" "TTToHcToWA_1ta2mu_MHc150_MA15" "TTToHcToWA_1ta2mu_MHc150_MA20" "TTToHcToWA_1ta2mu_MHc150_MA25" "TTToHcToWA_1ta2mu_MHc150_MA30" "TTToHcToWA_1ta2mu_MHc150_MA35" "TTToHcToWA_1ta2mu_MHc160_MA15" "TTToHcToWA_1ta2mu_MHc160_MA20" "TTToHcToWA_1ta2mu_MHc160_MA25" "TTToHcToWA_1ta2mu_MHc160_MA30" "TTToHcToWA_1ta2mu_MHc160_MA35")


declare -a Signal_All=("TTToHcToWA_1e2mu_MHc100_MA15" "TTToHcToWA_1e2mu_MHc100_MA20" "TTToHcToWA_1e2mu_MHc110_MA15" "TTToHcToWA_1e2mu_MHc110_MA20" "TTToHcToWA_1e2mu_MHc110_MA25" "TTToHcToWA_1e2mu_MHc110_MA30" "TTToHcToWA_1e2mu_MHc120_MA15" "TTToHcToWA_1e2mu_MHc120_MA20" "TTToHcToWA_1e2mu_MHc120_MA25" "TTToHcToWA_1e2mu_MHc120_MA30" "TTToHcToWA_1e2mu_MHc120_MA35" "TTToHcToWA_1e2mu_MHc130_MA15" "TTToHcToWA_1e2mu_MHc130_MA20" "TTToHcToWA_1e2mu_MHc130_MA25" "TTToHcToWA_1e2mu_MHc130_MA30" "TTToHcToWA_1e2mu_MHc130_MA35" "TTToHcToWA_1e2mu_MHc140_MA15" "TTToHcToWA_1e2mu_MHc140_MA20" "TTToHcToWA_1e2mu_MHc140_MA25" "TTToHcToWA_1e2mu_MHc140_MA30" "TTToHcToWA_1e2mu_MHc140_MA35" "TTToHcToWA_1e2mu_MHc150_MA15" "TTToHcToWA_1e2mu_MHc150_MA20" "TTToHcToWA_1e2mu_MHc150_MA25" "TTToHcToWA_1e2mu_MHc150_MA30" "TTToHcToWA_1e2mu_MHc150_MA35" "TTToHcToWA_1e2mu_MHc160_MA15" "TTToHcToWA_1e2mu_MHc160_MA20" "TTToHcToWA_1e2mu_MHc160_MA25" "TTToHcToWA_1e2mu_MHc160_MA30" "TTToHcToWA_1e2mu_MHc160_MA35" "TTToHcToWA_2l2mu_MHc100_MA15" "TTToHcToWA_2l2mu_MHc100_MA20" "TTToHcToWA_2l2mu_MHc110_MA15" "TTToHcToWA_2l2mu_MHc110_MA20" "TTToHcToWA_2l2mu_MHc110_MA25" "TTToHcToWA_2l2mu_MHc110_MA30" "TTToHcToWA_2l2mu_MHc120_MA15" "TTToHcToWA_2l2mu_MHc120_MA20" "TTToHcToWA_2l2mu_MHc120_MA25" "TTToHcToWA_2l2mu_MHc120_MA30" "TTToHcToWA_2l2mu_MHc120_MA35" "TTToHcToWA_2l2mu_MHc130_MA15" "TTToHcToWA_2l2mu_MHc130_MA20" "TTToHcToWA_2l2mu_MHc130_MA25" "TTToHcToWA_2l2mu_MHc130_MA30" "TTToHcToWA_2l2mu_MHc130_MA35" "TTToHcToWA_2l2mu_MHc140_MA15" "TTToHcToWA_2l2mu_MHc140_MA20" "TTToHcToWA_2l2mu_MHc140_MA25" "TTToHcToWA_2l2mu_MHc140_MA30" "TTToHcToWA_2l2mu_MHc140_MA35" "TTToHcToWA_2l2mu_MHc150_MA15" "TTToHcToWA_2l2mu_MHc150_MA20" "TTToHcToWA_2l2mu_MHc150_MA25" "TTToHcToWA_2l2mu_MHc150_MA30" "TTToHcToWA_2l2mu_MHc150_MA35" "TTToHcToWA_2l2mu_MHc160_MA15" "TTToHcToWA_2l2mu_MHc160_MA20" "TTToHcToWA_2l2mu_MHc160_MA25" "TTToHcToWA_2l2mu_MHc160_MA30" "TTToHcToWA_2l2mu_MHc160_MA35" "TTToHcToWA_3mu_MHc100_MA15" "TTToHcToWA_3mu_MHc100_MA20" "TTToHcToWA_3mu_MHc110_MA15" "TTToHcToWA_3mu_MHc110_MA20" "TTToHcToWA_3mu_MHc110_MA25" "TTToHcToWA_3mu_MHc110_MA30" "TTToHcToWA_3mu_MHc120_MA15" "TTToHcToWA_3mu_MHc120_MA20" "TTToHcToWA_3mu_MHc120_MA25" "TTToHcToWA_3mu_MHc120_MA30" "TTToHcToWA_3mu_MHc120_MA35" "TTToHcToWA_3mu_MHc130_MA15" "TTToHcToWA_3mu_MHc130_MA20" "TTToHcToWA_3mu_MHc130_MA25" "TTToHcToWA_3mu_MHc130_MA30" "TTToHcToWA_3mu_MHc130_MA35" "TTToHcToWA_3mu_MHc140_MA15" "TTToHcToWA_3mu_MHc140_MA20" "TTToHcToWA_3mu_MHc140_MA25" "TTToHcToWA_3mu_MHc140_MA30" "TTToHcToWA_3mu_MHc140_MA35" "TTToHcToWA_3mu_MHc150_MA15" "TTToHcToWA_3mu_MHc150_MA20" "TTToHcToWA_3mu_MHc150_MA25" "TTToHcToWA_3mu_MHc150_MA30" "TTToHcToWA_3mu_MHc150_MA35" "TTToHcToWA_3mu_MHc160_MA15" "TTToHcToWA_3mu_MHc160_MA20" "TTToHcToWA_3mu_MHc160_MA25" "TTToHcToWA_3mu_MHc160_MA30" "TTToHcToWA_3mu_MHc160_MA35" "TTToHcToWA_1ta2mu_MHc100_MA15" "TTToHcToWA_1ta2mu_MHc100_MA20" "TTToHcToWA_1ta2mu_MHc110_MA15" "TTToHcToWA_1ta2mu_MHc110_MA20" "TTToHcToWA_1ta2mu_MHc110_MA25" "TTToHcToWA_1ta2mu_MHc110_MA30" "TTToHcToWA_1ta2mu_MHc120_MA15" "TTToHcToWA_1ta2mu_MHc120_MA20" "TTToHcToWA_1ta2mu_MHc120_MA25" "TTToHcToWA_1ta2mu_MHc120_MA30" "TTToHcToWA_1ta2mu_MHc120_MA35" "TTToHcToWA_1ta2mu_MHc130_MA15" "TTToHcToWA_1ta2mu_MHc130_MA20" "TTToHcToWA_1ta2mu_MHc130_MA25" "TTToHcToWA_1ta2mu_MHc130_MA30" "TTToHcToWA_1ta2mu_MHc130_MA35" "TTToHcToWA_1ta2mu_MHc140_MA15" "TTToHcToWA_1ta2mu_MHc140_MA20" "TTToHcToWA_1ta2mu_MHc140_MA25" "TTToHcToWA_1ta2mu_MHc140_MA30" "TTToHcToWA_1ta2mu_MHc140_MA35" "TTToHcToWA_1ta2mu_MHc150_MA15" "TTToHcToWA_1ta2mu_MHc150_MA20" "TTToHcToWA_1ta2mu_MHc150_MA25" "TTToHcToWA_1ta2mu_MHc150_MA30" "TTToHcToWA_1ta2mu_MHc150_MA35" "TTToHcToWA_1ta2mu_MHc160_MA15" "TTToHcToWA_1ta2mu_MHc160_MA20" "TTToHcToWA_1ta2mu_MHc160_MA25" "TTToHcToWA_1ta2mu_MHc160_MA30" "TTToHcToWA_1ta2mu_MHc160_MA35")
