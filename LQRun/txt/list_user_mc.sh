#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################

#declare -a example=("tthwA_1l2ta_hc130A30" "tthwA_1ta2mu_hc130A30")  #"tthwA_4l_hc130A30")
#declare -a example=("tthwA_1ta2mu_hc130A30" "tthwA_1l2ta_hc130A30" "tthwA_4l_hc130A30")
declare -a example=("tthwA_1ta2mu_hc130A30" )
declare -a tamumu=("tthwA_1ta2mu_hc130A30" "tthwA_1ta2mu_hc155A30")
declare -a fourlep=("tthwA_4l2j_hc130A30" "tthwA_4l2j_hc155A30")
declare -a NewSamples=("tthwA_1e2mu_hc130A30_NoGenCut" "tthwA_1e2mu_hc155A30_NoGenCut") 
#declare -a NewSamples=("tthwA_3mu_hc130A30_NoGenCut") 

#declare -a DY=("DYJets" 'TT_powheg') 
declare -a DY=("DYJets") 
declare -a BtagEffSample=("TT_powheg")

declare -a tmplist=('WpWp_qcd_madgraph' 'ZG_llG_MCatNLO' 'ZZ_llnunu_powheg' 'ZZ_llqq_MCatNLO' 'ZZ_llll_MCatNLO' 'ZZ_llll_powheg' 'ZZ_pythia8' 'ttHnobb_Powheg' 'ttHtobb_Powheg')

declare -a AllSample=('WJets' 'DYJets_10to50' 'DYJets' 'SingleTop_s' 'SingleTbar_t' 'SingleTop_t' 'SingleTbar_tW' 'SingleTop_tW' 'TT_powheg' 'ZZ' 'WZ' 'WW' 'WGtoLNuG' 'WGtoLNuEE' 'WGtoLNuMM' 'ZGto2LG' 'ttH_nonbb' 'ttH_bb' 'ttW' 'ttZ') 
###Backgound
###declare -a input_samples=( 'DY10to50_MCatNLO' 'DY50plus_MCatNLO' 'WJets_MCatNLO' "TT_MCatNLO" 'TTG_MCatNLO' 'WZ_lllnu_MCatNLO' 'WZ_llqq_MCatNLO' 'ZZ_llqq_MCatNLO' 'ZZ_llll_MCatNLO' 'singletop_tbarW_Powheg' 'singletop_tW_Powheg' 'singletop_tbar_Powheg' 'singletop_t_Powheg' 'singletop_s_MCatNLO' 'ttHtobb_Powheg' 'ttHnobb_Powheg' 'ttWJetsToLNu_MCatNLO' 'ttZToQQ_MCatNLO' 'ttZToLLNuNu_MCatNLO' )
declare -a Analysis_bkg_3mu=( 'DY10to50_MCatNLO' 'DY50plus_MCatNLO' 'WJets_MCatNLO' 'TT_powheg' 'TTG_MCatNLO' 'WW_llnn_powheg' 'WZ_lllnu_MCatNLO' 'WZ_llqq_MCatNLO' 'ZZ_llqq_MCatNLO' 'ZZ_llll_MCatNLO' 'WG_lnuG_madgraph' 'ZG_llG_MCatNLO' 'singletop_tbarW_Powheg' 'singletop_tW_Powheg' 'singletop_tbar_Powheg' 'singletop_t_Powheg' 'singletop_s_MCatNLO' 'ttHtobb_Powheg' 'ttHnobb_Powheg' 'ttWJetsToLNu_MCatNLO' 'ttWJetsToQQ_MCatNLO' 'ttZToQQ_MCatNLO' 'ttZToLLNuNu_MCatNLO' )
declare -a Analysis_bkg_1e2mu=( 'DY10to50_MCatNLO' 'DY50plus_MCatNLO' 'WJets_MCatNLO' 'TT_powheg' 'TTG_MCatNLO' 'WW_llnn_powheg' 'WZ_lllnu_MCatNLO' 'WZ_llqq_MCatNLO' 'ZZ_llqq_MCatNLO' 'ZZ_llll_MCatNLO' 'WG_lnuG_madgraph' 'ZG_llG_MCatNLO' 'singletop_tbarW_Powheg' 'singletop_tW_Powheg' 'singletop_tbar_Powheg' 'singletop_t_Powheg' 'singletop_s_MCatNLO' 'ttHtobb_Powheg' 'ttHnobb_Powheg' 'ttWJetsToLNu_MCatNLO' 'ttWJetsToQQ_MCatNLO' 'ttZToQQ_MCatNLO' 'ttZToLLNuNu_MCatNLO' )

declare -a Analysis_bkg_test=( 'WZTo3LNu_powheg' )
#declare -a Analysis_bkg_test=( 'ttZ' )
declare -a Analysis_bkg_TT=('TT_powheg')

###Signal
declare -a Analysis_sig_All=("tthwA_3mu_hc90A5" "tthwA_1e2mu_hc90A5" "tthwA_3mu_hc130A5" "tthwA_1e2mu_hc130A5" "tthwA_3mu_hc130A30" "tthwA_1e2mu_hc130A30" "tthwA_3mu_hc155A30" "tthwA_1e2mu_hc155A30" "tthwzp_3mu_hc155zp5" "tthwzp_1e2mu_hc155zp5")
declare -a Analysis_sig_1e2mu=("tthwA_1e2mu_hc90A5" "tthwA_1e2mu_hc130A5" "tthwA_1e2mu_hc130A30" "tthwA_1e2mu_hc155A30" "tthwzp_1e2mu_hc155zp5")
declare -a Analysis_sig_3mu=("tthwA_3mu_hc90A5" "tthwA_3mu_hc130A5" "tthwA_3mu_hc130A30" "tthwA_3mu_hc155A30" "tthwzp_3mu_hc155zp5")
declare -a Analysis_sig_incl=("tthwA_2l6j_hc130A30" "tthwA_4l2j_hc130A30")

declare -a tthwA_1e2mu=("tthwA_1e2mu_hc90A5" "tthwA_1e2mu_hc130A5" "tthwA_1e2mu_hc130A30" "tthwA_1e2mu_hc155A30")
declare -a tthwA_3mu=("tthwA_3mu_hc90A5" "tthwA_3mu_hc130A5" "tthwA_3mu_hc130A30" "tthwA_3mu_hc155A30")

declare -a Analysis_sig_test=("tthwA_4l2j_hc155A30")
declare -a Analysis_sig_test1e2mu=("tthwA_1e2mu_hc130A30")
declare -a Analysis_sig_test3mu=("tthwA_3mu_hc130A30")

#CR##########################
#CR-DiLep+jets
declare -a CR_MuMu=('WJets' 'DYJets_10to50' 'DYJets' 'SingleTop_s' 'SingleTbar_t' 'SingleTop_t' 'SingleTbar_tW' 'SingleTop_tW' 'TT_powheg' 'WZ' 'ZZ' 'WW' 'WGtoLNuG' 'WGtoLNuEE' 'WGtoLNuMM' 'ZGto2LG' 'ttH_nonbb' 'ttH_bb' 'ttW' 'ttZ' 'QCD_Pt-15to20_MuEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_Pt-30to50_MuEnriched' 'QCD_Pt-50to80_MuEnriched' 'QCD_Pt-80to120_MuEnriched' 'QCD_Pt-120to170_MuEnriched' 'QCD_Pt-170to300_MuEnriched' 'QCD_Pt-300to470_MuEnriched' 'QCD_Pt-470to600_MuEnriched' 'QCD_Pt-600to800_MuEnriched' 'QCD_Pt-800to1000_MuEnriched' 'QCD_Pt-1000toInf_MuEnriched' 'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-30to50_EMEnriched' 'QCD_Pt-50to80_EMEnriched' 'QCD_Pt-80to120_EMEnriched' 'QCD_Pt-120to170_EMEnriched' 'QCD_Pt-170to300_EMEnriched' 'QCD_Pt-300toInf_EMEnriched' )

declare -a CR_EMu=('WJets' 'DYJets_10to50' 'DYJets' 'SingleTop_s' 'SingleTbar_t' 'SingleTop_t' 'SingleTbar_tW' 'SingleTop_tW' 'TT_powheg' 'WZ' 'ZZ' 'WW' 'WGtoLNuG' 'WGtoLNuEE' 'WGtoLNuMM' 'ZGto2LG' 'ttH_nonbb' 'ttH_bb' 'ttW' 'ttZ' 'QCD_Pt-15to20_MuEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_Pt-30to50_MuEnriched' 'QCD_Pt-50to80_MuEnriched' 'QCD_Pt-80to120_MuEnriched' 'QCD_Pt-120to170_MuEnriched' 'QCD_Pt-170to300_MuEnriched' 'QCD_Pt-300to470_MuEnriched' 'QCD_Pt-470to600_MuEnriched' 'QCD_Pt-600to800_MuEnriched' 'QCD_Pt-800to1000_MuEnriched' 'QCD_Pt-1000toInf_MuEnriched' 'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-30to50_EMEnriched' 'QCD_Pt-50to80_EMEnriched' 'QCD_Pt-80to120_EMEnriched' 'QCD_Pt-120to170_EMEnriched' 'QCD_Pt-170to300_EMEnriched' 'QCD_Pt-300toInf_EMEnriched' )


#CR-TriLep
declare -a CR_EMuMu=('DY10to50_MCatNLO' 'DY50plus_MCatNLO' 'WJets_MCatNLO' 'WZ_lllnu_MCatNLO' 'WZ_llqq_MCatNLO' 'ZZ_llll_MCatNLO' 'ZZ_llqq_MCatNLO' 'WW_llnn_powheg' 'ZG_llG_MCatNLO' 'WG_lnuG_madgraph' 'TT_powheg' 'TTG_MCatNLO' 'ttZToQQ_MCatNLO' 'ttZToLLNuNu_MCatNLO' 'ttWJetsToLNu_MCatNLO' 'ttWJetsToQQ_MCatNLO' 'ttHtobb_Powheg' 'ttHnobb_Powheg' 'singletop_tbarW_Powheg' 'singletop_tbarW_Powheg' 'singletop_tbar_Powheg' 'singletop_t_Powheg' 'singletop_s_MCatNLO')
declare -a CR_MuMuMu=('DY10to50_MCatNLO' 'DY50plus_MCatNLO' 'WJets_MCatNLO' 'WZ_lllnu_MCatNLO' 'WZ_llqq_MCatNLO' 'ZZ_llll_MCatNLO' 'ZZ_llqq_MCatNLO' 'WW_llnn_powheg' 'ZG_llG_MCatNLO' 'WG_lnuG_madgraph' 'TT_powheg' 'TTG_MCatNLO' 'ttZToQQ_MCatNLO' 'ttZToLLNuNu_MCatNLO' 'ttWJetsToLNu_MCatNLO' 'ttWJetsToQQ_MCatNLO' 'ttHtobb_Powheg' 'ttHnobb_Powheg' 'singletop_tbarW_Powheg' 'singletop_tbarW_Powheg' 'singletop_tbar_Powheg' 'singletop_t_Powheg' 'singletop_s_MCatNLO')


###QCD
declare -a QCD_mu=('QCD_Pt-15to20_MuEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_Pt-30to50_MuEnriched' 'QCD_Pt-50to80_MuEnriched' 'QCD_Pt-80to120_MuEnriched' 'QCD_Pt-120to170_MuEnriched' 'QCD_Pt-170to300_MuEnriched' 'QCD_Pt-300to470_MuEnriched' 'QCD_Pt-470to600_MuEnriched' 'QCD_Pt-600to800_MuEnriched' 'QCD_Pt-800to1000_MuEnriched' 'QCD_Pt-1000toInf_MuEnriched' )
declare -a QCD_em=('QCD_Pt-20to30_EMEnriched' 'QCD_Pt-30to50_EMEnriched' 'QCD_Pt-50to80_EMEnriched' 'QCD_Pt-80to120_EMEnriched' 'QCD_Pt-120to170_EMEnriched' 'QCD_Pt-170to300_EMEnriched' 'QCD_Pt-300toInf_EMEnriched' )

##declare -a QCD_bcToE=('QCD_20to30_bcToE_pythia8' 'QCD_30to80_bcToE_pythia8' 'QCD_80to170_bcToE_pythia8' 'QCD_170to250_bcToE_pythia8' 'QCD_250toInf_bcToE_pythia8')
##declare -a QCD_DoubleEM=('QCD_DoubleEM_30to40_pythia8' 'QCD_DoubleEM_30toInf_pythia8' 'QCD_DoubleEM_40toInf_pythia8')
declare -a QCD_All=('QCD_Pt-15to20_MuEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_Pt-30to50_MuEnriched' 'QCD_Pt-50to80_MuEnriched' 'QCD_Pt-80to120_MuEnriched' 'QCD_Pt-120to170_MuEnriched' 'QCD_Pt-170to300_MuEnriched' 'QCD_Pt-300to470_MuEnriched' 'QCD_Pt-470to600_MuEnriched' 'QCD_Pt-600to800_MuEnriched' 'QCD_Pt-800to1000_MuEnriched' 'QCD_Pt-1000toInf_MuEnriched' 'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-30to50_EMEnriched' 'QCD_Pt-50to80_EMEnriched' 'QCD_Pt-80to120_EMEnriched' 'QCD_Pt-120to170_EMEnriched' 'QCD_Pt-170to300_EMEnriched' 'QCD_Pt-300toInf_EMEnriched' )



#DataValidation-SingleE
declare -a Validation_SEle=("WJets_MCatNLO" "WG_lnuG_madgraph" "TT_MG5" "ZG_llG_MCatNLO" "WW_pythia8" "DY10to50_MCatNLO" "DY50plus_MCatNLO" 'QCD_em15to20_pythia8' 'QCD_em20to30_pythia8' 'QCD_em50to80_pythia8' 'QCD_em80to120_pythia8' 'QCD_em120to170_pythia8' 'QCD_em170to300_pythia8' 'QCD_20to30_bcToE_pythia8' 'QCD_30to80_bcToE_pythia8' 'QCD_80to170_bcToE_pythia8' 'QCD_170to250_bcToE_pythia8' 'QCD_250toInf_bcToE_pythia8') 

#MuMubj Test
declare -a MuMubj_bkg=("TT_powheg" "DY50plus_MCatNLO" "DY10to50_MCatNLO" "singletop_tbar_Powheg" "singletop_s_MCatNLO" "singletop_tW_Powheg" "singletop_tbarW_Powheg" "singletop_t_Powheg" "TTG_MCatNLO" "ttZToLLNuNu_MCatNLO" "ttWJetsToLNu_MCatNLO" "ttZToQQ_MCatNLO" "ttHnobb_Powheg" "WZ_lllnu_MCatNLO" "WZ_llqq_MCatNLO" "WW_pythia8" "ZZ_llll_MCatNLO" "ZG_llG_MCatNLO" "ZZ_llqq_MCatNLO")
