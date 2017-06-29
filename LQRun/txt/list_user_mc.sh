#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################
declare -a hnfail=('HNMumEp_50' 'HNMumEp_500' 'HNMumEp_60' 'HNMumMum_100' 'HNMumMum_1100' 'HNMumMum_1500' 'HNMumMum_200' 'HNMumMum_40' 'HNMumMum_50' 'HNMumMum_500' 'HNMumMum_60' 'HNMoriondLL_Tchannel_MumMum_100' 'HNMoriondLL_Tchannel_MumMum_1100' 'HNMoriondLL_Tchannel_MumMum_200' 'HNMoriondLL_Tchannel_MumMum_500' 'HNMumMup_100' 'HNMumMup_1100' 'HNMumMup_1500' 'HNMumMup_200' 'HNMumMup_40' 'HNMumMup_50' 'HNMumMup_500' 'HNMumMup_60' 'HNMupEm_100' 'HNMupEm_1100' 'HNMupEm_1500' 'HNMupEm_200' 'HNMupEm_40' 'HNMupEm_50' 'HNMupEm_500' 'HNMupEm_60' 'HNMupEp_100')

declare -a hn_mme=('HN_SSSF_MuMuE_1000' 'HN_SSSF_MuMuE_100' 'HN_SSSF_MuMuE_10' 'HN_SSSF_MuMuE_150' 'HN_SSSF_MuMuE_200' 'HN_SSSF_MuMuE_20' 'HN_SSSF_MuMuE_300' 'HN_SSSF_MuMuE_30' 'HN_SSSF_MuMuE_400' 'HN_SSSF_MuMuE_40' 'HN_SSSF_MuMuE_500' 'HN_SSSF_MuMuE_50' 'HN_SSSF_MuMuE_5' 'HN_SSSF_MuMuE_60' 'HN_SSSF_MuMuE_700' 'HN_SSSF_MuMuE_70' 'HN_SSSF_MuMuE_90')
declare -a qcd_ee=('qcd_15to20_bctoe' 'qcd_170to250_bctoe' 'qcd_20to30_bctoe' 'qcd_250toinf_bctoe' 'qcd_30to80_bctoe' 'qcd_80to170_bctoe' 'QCD_Pt-120to170_EMEnriched' 'QCD_Pt-170to300_EMEnriched' 'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-300toInf_EMEnriched'  'QCD_Pt-30to50_EMEnriched' 'QCD_Pt-50to80_EMEnriched' 'QCD_Pt-80to120_EMEnriched')

declare -a hn_ll_ee2=('HN_t_ch_EmEm_100_official' 'HN_t_ch_EmEm_1100_official' 'HN_t_ch_EmEm_200_official' 'HN_t_ch_EmEm_500_official' 'HN_t_ch_EpEp_100_official' 'HN_t_ch_EpEp_1100_official' 'HN_t_ch_EpEp_200_official' 'HN_t_ch_EpEp_500_official' 'HNpair_ElEl_WR5000_Zp1000_HN100_official' 'HNpair_ElEl_WR5000_Zp1000_HN200_official' 'HNpair_ElEl_WR5000_Zp1000_HN300_official' 'HNpair_ElEl_WR5000_Zp1000_HN400_official' 'HNpair_ElEl_WR5000_Zp1500_HN100_official' 'HNpair_ElEl_WR5000_Zp1500_HN200_official' 'HNpair_ElEl_WR5000_Zp1500_HN300_official' 'HNpair_ElEl_WR5000_Zp1500_HN400_official' 'HNpair_ElEl_WR5000_Zp1500_HN500_official' 'HNpair_ElEl_WR5000_Zp1500_HN600_official' 'HNpair_ElEl_WR5000_Zp1500_HN700_official' 'HNpair_ElEl_WR5000_Zp2000_HN100_official' 'HNpair_ElEl_WR5000_Zp2000_HN200_official' 'HNpair_ElEl_WR5000_Zp2000_HN300_official' 'HNpair_ElEl_WR5000_Zp2000_HN400_official' 'HNpair_ElEl_WR5000_Zp2000_HN500_official' 'HNpair_ElEl_WR5000_Zp2000_HN600_official' 'HNpair_ElEl_WR5000_Zp2000_HN700_official' 'HNpair_ElEl_WR5000_Zp2000_HN800_official' 'HNpair_ElEl_WR5000_Zp2000_HN900_official' 'HNpair_ElEl_WR5000_Zp2500_HN1000_official' 'HNpair_ElEl_WR5000_Zp2500_HN100_official' 'HNpair_ElEl_WR5000_Zp2500_HN1100_official' 'HNpair_ElEl_WR5000_Zp2500_HN1200_official' 'HNpair_ElEl_WR5000_Zp2500_HN200_official' 'HNpair_ElEl_WR5000_Zp2500_HN300_official' 'HNpair_ElEl_WR5000_Zp2500_HN400_official' 'HNpair_ElEl_WR5000_Zp2500_HN500_official' 'HNpair_ElEl_WR5000_Zp2500_HN600_official' 'HNpair_ElEl_WR5000_Zp2500_HN700_official' 'HNpair_ElEl_WR5000_Zp2500_HN800_official' 'HNpair_ElEl_WR5000_Zp2500_HN900_official' 'HNpair_ElEl_WR5000_Zp3000_HN1000_official' 'HNpair_ElEl_WR5000_Zp3000_HN100_official' 'HNpair_ElEl_WR5000_Zp3000_HN1100_official' 'HNpair_ElEl_WR5000_Zp3000_HN1200_official' 'HNpair_ElEl_WR5000_Zp3000_HN1300_official' 'HNpair_ElEl_WR5000_Zp3000_HN1400_official' 'HNpair_ElEl_WR5000_Zp3000_HN200_official' 'HNpair_ElEl_WR5000_Zp3000_HN300_official' 'HNpair_ElEl_WR5000_Zp3000_HN400_official' 'HNpair_ElEl_WR5000_Zp3000_HN500_official' 'HNpair_ElEl_WR5000_Zp3000_HN600_official' 'HNpair_ElEl_WR5000_Zp3000_HN700_official' 'HNpair_ElEl_WR5000_Zp3000_HN800_official' 'HNpair_ElEl_WR5000_Zp3000_HN900_official' 'HNpair_ElEl_WR5000_Zp4000_HN1000_official' 'HNpair_ElEl_WR5000_Zp4000_HN100_official' 'HNpair_ElEl_WR5000_Zp4000_HN1100_official' 'HNpair_ElEl_WR5000_Zp4000_HN1200_official' 'HNpair_ElEl_WR5000_Zp4000_HN1300_official' 'HNpair_ElEl_WR5000_Zp4000_HN1400_official' 'HNpair_ElEl_WR5000_Zp4000_HN1500_official' 'HNpair_ElEl_WR5000_Zp4000_HN1700_official' 'HNpair_ElEl_WR5000_Zp4000_HN1800_official' 'HNpair_ElEl_WR5000_Zp4000_HN1900_official' 'HNpair_ElEl_WR5000_Zp4000_HN200_official' 'HNpair_ElEl_WR5000_Zp4000_HN300_official' 'HNpair_ElEl_WR5000_Zp4000_HN400_official' 'HNpair_ElEl_WR5000_Zp4000_HN500_official' 'HNpair_ElEl_WR5000_Zp4000_HN600_official' 'HNpair_ElEl_WR5000_Zp4000_HN700_official' 'HNpair_ElEl_WR5000_Zp4000_HN800_official' 'HNpair_ElEl_WR5000_Zp4000_HN900_official' 'HNpair_ElEl_WR5000_Zp500_HN100_official' 'HNpair_ElEl_WR5000_Zp500_HN200_official' 'HNpair_ElEl_WR5000_Zp750_HN100_official' 'HNpair_ElEl_WR5000_Zp750_HN200_official' 'HNpair_ElEl_WR5000_Zp750_HN300_official'  'HNEmEm_100' 'HNEmEm_1100' 'HNEmEm_1500' 'HNEmEm_200' 'HNEmEm_40' 'HNEmEm_50' 'HNEmEm_500' 'HNEmEm_60' 'HNEmEp_100' 'HNEmEp_1100' 'HNEmEp_1500' 'HNEmEp_200' 'HNEmEp_40' 'HNEmEp_50' 'HNEmEp_500' 'HNEmEp_60' 'HNEpEm_100' 'HNEpEm_1100' 'HNEpEm_1500' 'HNEpEm_200' 'HNEpEm_40' 'HNEpEm_50' 'HNEpEm_500' 'HNEpEm_60' 'HNEpEp_100' 'HNEpEp_1100' 'HNEpEp_1500' 'HNEpEp_200' 'HNEpEp_40' 'HNEpEp_50' 'HNEpEp_500' 'HNEpEp_60' )

declare -a hn_ll_ee1=( 'HNEmEm_100' 'HNEmEm_1100' 'HNEmEm_1500' 'HNEmEm_200' 'HNEmEm_40' 'HNEmEm_50' 'HNEmEm_500' 'HNEmEm_60' 'HNEmEp_100' 'HNEmEp_1100' 'HNEmEp_1500' 'HNEmEp_200' 'HNEmEp_40' 'HNEmEp_50' 'HNEmEp_500' 'HNEmEp_60' 'HNEpEm_100' 'HNEpEm_1100' 'HNEpEm_1500' 'HNEpEm_200' 'HNEpEm_40' 'HNEpEm_50' 'HNEpEm_500' 'HNEpEm_60' 'HNEpEp_100' 'HNEpEp_1100' 'HNEpEp_1500' 'HNEpEp_200' 'HNEpEp_40' 'HNEpEp_50' 'HNEpEp_500' 'HNEpEp_60' )                                                                                                                                           

declare -a hn_ll_mm1=( 'HNMumMum_100' 'HNMumMum_1100' 'HNMumMum_1500' 'HNMumMum_200' 'HNMumMum_40' 'HNMumMum_50' 'HNMumMum_500' 'HNMumMum_60' 'HNMumMup_100' 'HNMumMup_1100' 'HNMumMup_1500' 'HNMumMup_200' 'HNMumMup_40' 'HNMumMup_50' 'HNMumMup_500' 'HNMumMup_60' 'HNMupMum_100' 'HNMupMum_1100' 'HNMupMum_1500' 'HNMupMum_200' 'HNMupMum_40' 'HNMupMum_50' 'HNMupMum_500' 'HNMupMum_60' 'HNMupMup_100' 'HNMupMup_1100' 'HNMupMup_1500' 'HNMupMup_200' 'HNMupMup_40' 'HNMupMup_50' 'HNMupMup_500' 'HNMupMup_60' )       


declare -a hn_pair_all=('HNpair_ElEl_WR5000_Zp1000_HN100_official' 'HNpair_ElEl_WR5000_Zp1000_HN200_official' 'HNpair_ElEl_WR5000_Zp1000_HN300_official' 'HNpair_ElEl_WR5000_Zp1000_HN400_official' 'HNpair_ElEl_WR5000_Zp1500_HN100_official' 'HNpair_ElEl_WR5000_Zp1500_HN200_official' 'HNpair_ElEl_WR5000_Zp1500_HN300_official' 'HNpair_ElEl_WR5000_Zp1500_HN400_official' 'HNpair_ElEl_WR5000_Zp1500_HN500_official' 'HNpair_ElEl_WR5000_Zp1500_HN600_official' 'HNpair_ElEl_WR5000_Zp1500_HN700_official' 'HNpair_ElEl_WR5000_Zp2000_HN100_official' 'HNpair_ElEl_WR5000_Zp2000_HN200_official' 'HNpair_ElEl_WR5000_Zp2000_HN300_official' 'HNpair_ElEl_WR5000_Zp2000_HN400_official' 'HNpair_ElEl_WR5000_Zp2000_HN500_official' 'HNpair_ElEl_WR5000_Zp2000_HN600_official' 'HNpair_ElEl_WR5000_Zp2000_HN700_official' 'HNpair_ElEl_WR5000_Zp2000_HN800_official' 'HNpair_ElEl_WR5000_Zp2000_HN900_official' 'HNpair_ElEl_WR5000_Zp2500_HN1000_official' 'HNpair_ElEl_WR5000_Zp2500_HN100_official' 'HNpair_ElEl_WR5000_Zp2500_HN1100_official' 'HNpair_ElEl_WR5000_Zp2500_HN1200_official' 'HNpair_ElEl_WR5000_Zp2500_HN200_official' 'HNpair_ElEl_WR5000_Zp2500_HN300_official' 'HNpair_ElEl_WR5000_Zp2500_HN400_official' 'HNpair_ElEl_WR5000_Zp2500_HN500_official' 'HNpair_ElEl_WR5000_Zp2500_HN600_official' 'HNpair_ElEl_WR5000_Zp2500_HN700_official' 'HNpair_ElEl_WR5000_Zp2500_HN800_official' 'HNpair_ElEl_WR5000_Zp2500_HN900_official' 'HNpair_ElEl_WR5000_Zp3000_HN1000_official' 'HNpair_ElEl_WR5000_Zp3000_HN100_official' 'HNpair_ElEl_WR5000_Zp3000_HN1100_official' 'HNpair_ElEl_WR5000_Zp3000_HN1200_official' 'HNpair_ElEl_WR5000_Zp3000_HN1300_official' 'HNpair_ElEl_WR5000_Zp3000_HN1400_official' 'HNpair_ElEl_WR5000_Zp3000_HN200_official' 'HNpair_ElEl_WR5000_Zp3000_HN300_official' 'HNpair_ElEl_WR5000_Zp3000_HN400_official' 'HNpair_ElEl_WR5000_Zp3000_HN500_official' 'HNpair_ElEl_WR5000_Zp3000_HN600_official' 'HNpair_ElEl_WR5000_Zp3000_HN700_official' 'HNpair_ElEl_WR5000_Zp3000_HN800_official' 'HNpair_ElEl_WR5000_Zp3000_HN900_official' 'HNpair_ElEl_WR5000_Zp4000_HN1000_official' 'HNpair_ElEl_WR5000_Zp4000_HN100_official' 'HNpair_ElEl_WR5000_Zp4000_HN1100_official' 'HNpair_ElEl_WR5000_Zp4000_HN1200_official' 'HNpair_ElEl_WR5000_Zp4000_HN1300_official' 'HNpair_ElEl_WR5000_Zp4000_HN1400_official' 'HNpair_ElEl_WR5000_Zp4000_HN1500_official' 'HNpair_ElEl_WR5000_Zp4000_HN1700_official' 'HNpair_ElEl_WR5000_Zp4000_HN1800_official' 'HNpair_ElEl_WR5000_Zp4000_HN1900_official' 'HNpair_ElEl_WR5000_Zp4000_HN200_official' 'HNpair_ElEl_WR5000_Zp4000_HN300_official' 'HNpair_ElEl_WR5000_Zp4000_HN400_official' 'HNpair_ElEl_WR5000_Zp4000_HN500_official' 'HNpair_ElEl_WR5000_Zp4000_HN600_official' 'HNpair_ElEl_WR5000_Zp4000_HN700_official' 'HNpair_ElEl_WR5000_Zp4000_HN800_official' 'HNpair_ElEl_WR5000_Zp4000_HN900_official' 'HNpair_ElEl_WR5000_Zp500_HN100_official' 'HNpair_ElEl_WR5000_Zp500_HN200_official' 'HNpair_ElEl_WR5000_Zp750_HN100_official' 'HNpair_ElEl_WR5000_Zp750_HN200_official' 'HNpair_ElEl_WR5000_Zp750_HN300_official' 'HNpair_MuMu_WR5000_Zp1000_HN100_official' 'HNpair_MuMu_WR5000_Zp1000_HN200_official' 'HNpair_MuMu_WR5000_Zp1000_HN300_official' 'HNpair_MuMu_WR5000_Zp1000_HN400_official' 'HNpair_MuMu_WR5000_Zp1500_HN100_official' 'HNpair_MuMu_WR5000_Zp1500_HN200_official' 'HNpair_MuMu_WR5000_Zp1500_HN300_official' 'HNpair_MuMu_WR5000_Zp1500_HN400_official' 'HNpair_MuMu_WR5000_Zp1500_HN500_official' 'HNpair_MuMu_WR5000_Zp1500_HN600_official' 'HNpair_MuMu_WR5000_Zp1500_HN700_official' 'HNpair_MuMu_WR5000_Zp2000_HN100_official' 'HNpair_MuMu_WR5000_Zp2000_HN200_official' 'HNpair_MuMu_WR5000_Zp2000_HN300_official' 'HNpair_MuMu_WR5000_Zp2000_HN400_official' 'HNpair_MuMu_WR5000_Zp2000_HN500_official' 'HNpair_MuMu_WR5000_Zp2000_HN600_official' 'HNpair_MuMu_WR5000_Zp2000_HN700_official' 'HNpair_MuMu_WR5000_Zp2000_HN800_official' 'HNpair_MuMu_WR5000_Zp2000_HN900_official' 'HNpair_MuMu_WR5000_Zp2500_HN1000_official' 'HNpair_MuMu_WR5000_Zp2500_HN100_official' 'HNpair_MuMu_WR5000_Zp2500_HN1100_official' 'HNpair_MuMu_WR5000_Zp2500_HN1200_official' 'HNpair_MuMu_WR5000_Zp2500_HN200_official' 'HNpair_MuMu_WR5000_Zp2500_HN300_official' 'HNpair_MuMu_WR5000_Zp2500_HN400_official' 'HNpair_MuMu_WR5000_Zp2500_HN500_official' 'HNpair_MuMu_WR5000_Zp2500_HN600_official' 'HNpair_MuMu_WR5000_Zp2500_HN700_official' 'HNpair_MuMu_WR5000_Zp2500_HN800_official' 'HNpair_MuMu_WR5000_Zp2500_HN900_official' 'HNpair_MuMu_WR5000_Zp3000_HN1000_official' 'HNpair_MuMu_WR5000_Zp3000_HN100_official' 'HNpair_MuMu_WR5000_Zp3000_HN1100_official' 'HNpair_MuMu_WR5000_Zp3000_HN1200_official' 'HNpair_MuMu_WR5000_Zp3000_HN1300_official' 'HNpair_MuMu_WR5000_Zp3000_HN1400_official' 'HNpair_MuMu_WR5000_Zp3000_HN200_official' 'HNpair_MuMu_WR5000_Zp3000_HN300_official' 'HNpair_MuMu_WR5000_Zp3000_HN400_official' 'HNpair_MuMu_WR5000_Zp3000_HN500_official' 'HNpair_MuMu_WR5000_Zp3000_HN600_official' 'HNpair_MuMu_WR5000_Zp3000_HN700_official' 'HNpair_MuMu_WR5000_Zp3000_HN800_official' 'HNpair_MuMu_WR5000_Zp3000_HN900_official' 'HNpair_MuMu_WR5000_Zp4000_HN1000_official' 'HNpair_MuMu_WR5000_Zp4000_HN100_official' 'HNpair_MuMu_WR5000_Zp4000_HN1100_official' 'HNpair_MuMu_WR5000_Zp4000_HN1200_official' 'HNpair_MuMu_WR5000_Zp4000_HN1300_official' 'HNpair_MuMu_WR5000_Zp4000_HN1400_official' 'HNpair_MuMu_WR5000_Zp4000_HN1500_official' 'HNpair_MuMu_WR5000_Zp4000_HN1600_official' 'HNpair_MuMu_WR5000_Zp4000_HN1700_official' 'HNpair_MuMu_WR5000_Zp4000_HN1800_official' 'HNpair_MuMu_WR5000_Zp4000_HN1900_official' 'HNpair_MuMu_WR5000_Zp4000_HN200_official' 'HNpair_MuMu_WR5000_Zp4000_HN300_official' 'HNpair_MuMu_WR5000_Zp4000_HN400_official' 'HNpair_MuMu_WR5000_Zp4000_HN500_official' 'HNpair_MuMu_WR5000_Zp4000_HN600_official' 'HNpair_MuMu_WR5000_Zp4000_HN700_official' 'HNpair_MuMu_WR5000_Zp4000_HN800_official' 'HNpair_MuMu_WR5000_Zp4000_HN900_official' 'HNpair_MuMu_WR5000_Zp500_HN100_official' 'HNpair_MuMu_WR5000_Zp500_HN200_official' 'HNpair_MuMu_WR5000_Zp750_HN100_official' 'HNpair_MuMu_WR5000_Zp750_HN200_official' 'HNpair_MuMu_WR5000_Zp750_HN300_official')


declare -a hn_pair_ee=('HNpair_ElEl_WR5000_Zp1000_HN100_official' 'HNpair_ElEl_WR5000_Zp1000_HN200_official' 'HNpair_ElEl_WR5000_Zp1000_HN300_official' 'HNpair_ElEl_WR5000_Zp1000_HN400_official' 'HNpair_ElEl_WR5000_Zp1500_HN100_official' 'HNpair_ElEl_WR5000_Zp1500_HN200_official' 'HNpair_ElEl_WR5000_Zp1500_HN300_official' 'HNpair_ElEl_WR5000_Zp1500_HN400_official' 'HNpair_ElEl_WR5000_Zp1500_HN500_official' 'HNpair_ElEl_WR5000_Zp1500_HN600_official' 'HNpair_ElEl_WR5000_Zp1500_HN700_official' 'HNpair_ElEl_WR5000_Zp2000_HN100_official' 'HNpair_ElEl_WR5000_Zp2000_HN200_official' 'HNpair_ElEl_WR5000_Zp2000_HN300_official' 'HNpair_ElEl_WR5000_Zp2000_HN400_official' 'HNpair_ElEl_WR5000_Zp2000_HN500_official' 'HNpair_ElEl_WR5000_Zp2000_HN600_official' 'HNpair_ElEl_WR5000_Zp2000_HN700_official' 'HNpair_ElEl_WR5000_Zp2000_HN800_official' 'HNpair_ElEl_WR5000_Zp2000_HN900_official' 'HNpair_ElEl_WR5000_Zp2500_HN1000_official' 'HNpair_ElEl_WR5000_Zp2500_HN100_official' 'HNpair_ElEl_WR5000_Zp2500_HN1100_official' 'HNpair_ElEl_WR5000_Zp2500_HN1200_official' 'HNpair_ElEl_WR5000_Zp2500_HN200_official' 'HNpair_ElEl_WR5000_Zp2500_HN300_official' 'HNpair_ElEl_WR5000_Zp2500_HN400_official' 'HNpair_ElEl_WR5000_Zp2500_HN500_official' 'HNpair_ElEl_WR5000_Zp2500_HN600_official' 'HNpair_ElEl_WR5000_Zp2500_HN700_official' 'HNpair_ElEl_WR5000_Zp2500_HN800_official' 'HNpair_ElEl_WR5000_Zp2500_HN900_official' 'HNpair_ElEl_WR5000_Zp3000_HN1000_official' 'HNpair_ElEl_WR5000_Zp3000_HN100_official' 'HNpair_ElEl_WR5000_Zp3000_HN1100_official' 'HNpair_ElEl_WR5000_Zp3000_HN1200_official' 'HNpair_ElEl_WR5000_Zp3000_HN1300_official' 'HNpair_ElEl_WR5000_Zp3000_HN1400_official' 'HNpair_ElEl_WR5000_Zp3000_HN200_official' 'HNpair_ElEl_WR5000_Zp3000_HN300_official' 'HNpair_ElEl_WR5000_Zp3000_HN400_official' 'HNpair_ElEl_WR5000_Zp3000_HN500_official' 'HNpair_ElEl_WR5000_Zp3000_HN600_official' 'HNpair_ElEl_WR5000_Zp3000_HN700_official' 'HNpair_ElEl_WR5000_Zp3000_HN800_official' 'HNpair_ElEl_WR5000_Zp3000_HN900_official' 'HNpair_ElEl_WR5000_Zp4000_HN1000_official' 'HNpair_ElEl_WR5000_Zp4000_HN100_official' 'HNpair_ElEl_WR5000_Zp4000_HN1100_official' 'HNpair_ElEl_WR5000_Zp4000_HN1200_official' 'HNpair_ElEl_WR5000_Zp4000_HN1300_official' 'HNpair_ElEl_WR5000_Zp4000_HN1400_official' 'HNpair_ElEl_WR5000_Zp4000_HN1500_official' 'HNpair_ElEl_WR5000_Zp4000_HN1700_official' 'HNpair_ElEl_WR5000_Zp4000_HN1800_official' 'HNpair_ElEl_WR5000_Zp4000_HN1900_official' 'HNpair_ElEl_WR5000_Zp4000_HN200_official' 'HNpair_ElEl_WR5000_Zp4000_HN300_official' 'HNpair_ElEl_WR5000_Zp4000_HN400_official' 'HNpair_ElEl_WR5000_Zp4000_HN500_official' 'HNpair_ElEl_WR5000_Zp4000_HN600_official' 'HNpair_ElEl_WR5000_Zp4000_HN700_official' 'HNpair_ElEl_WR5000_Zp4000_HN800_official' 'HNpair_ElEl_WR5000_Zp4000_HN900_official' 'HNpair_ElEl_WR5000_Zp500_HN100_official' 'HNpair_ElEl_WR5000_Zp500_HN200_official' 'HNpair_ElEl_WR5000_Zp750_HN100_official' 'HNpair_ElEl_WR5000_Zp750_HN200_official' 'HNpair_ElEl_WR5000_Zp750_HN300_official')
declare -a hn_pair_mm=('HNpair_MuMu_WR5000_Zp1000_HN100_official' 'HNpair_MuMu_WR5000_Zp1000_HN200_official' 'HNpair_MuMu_WR5000_Zp1000_HN300_official' 'HNpair_MuMu_WR5000_Zp1000_HN400_official' 'HNpair_MuMu_WR5000_Zp1500_HN100_official' 'HNpair_MuMu_WR5000_Zp1500_HN200_official' 'HNpair_MuMu_WR5000_Zp1500_HN300_official' 'HNpair_MuMu_WR5000_Zp1500_HN400_official' 'HNpair_MuMu_WR5000_Zp1500_HN500_official' 'HNpair_MuMu_WR5000_Zp1500_HN600_official' 'HNpair_MuMu_WR5000_Zp1500_HN700_official' 'HNpair_MuMu_WR5000_Zp2000_HN100_official' 'HNpair_MuMu_WR5000_Zp2000_HN200_official' 'HNpair_MuMu_WR5000_Zp2000_HN300_official' 'HNpair_MuMu_WR5000_Zp2000_HN400_official' 'HNpair_MuMu_WR5000_Zp2000_HN500_official' 'HNpair_MuMu_WR5000_Zp2000_HN600_official' 'HNpair_MuMu_WR5000_Zp2000_HN700_official' 'HNpair_MuMu_WR5000_Zp2000_HN800_official' 'HNpair_MuMu_WR5000_Zp2000_HN900_official' 'HNpair_MuMu_WR5000_Zp2500_HN1000_official' 'HNpair_MuMu_WR5000_Zp2500_HN100_official' 'HNpair_MuMu_WR5000_Zp2500_HN1100_official' 'HNpair_MuMu_WR5000_Zp2500_HN1200_official' 'HNpair_MuMu_WR5000_Zp2500_HN200_official' 'HNpair_MuMu_WR5000_Zp2500_HN300_official' 'HNpair_MuMu_WR5000_Zp2500_HN400_official' 'HNpair_MuMu_WR5000_Zp2500_HN500_official' 'HNpair_MuMu_WR5000_Zp2500_HN600_official' 'HNpair_MuMu_WR5000_Zp2500_HN700_official' 'HNpair_MuMu_WR5000_Zp2500_HN800_official' 'HNpair_MuMu_WR5000_Zp2500_HN900_official' 'HNpair_MuMu_WR5000_Zp3000_HN1000_official' 'HNpair_MuMu_WR5000_Zp3000_HN100_official' 'HNpair_MuMu_WR5000_Zp3000_HN1100_official' 'HNpair_MuMu_WR5000_Zp3000_HN1200_official' 'HNpair_MuMu_WR5000_Zp3000_HN1300_official' 'HNpair_MuMu_WR5000_Zp3000_HN1400_official' 'HNpair_MuMu_WR5000_Zp3000_HN200_official' 'HNpair_MuMu_WR5000_Zp3000_HN300_official' 'HNpair_MuMu_WR5000_Zp3000_HN400_official' 'HNpair_MuMu_WR5000_Zp3000_HN500_official' 'HNpair_MuMu_WR5000_Zp3000_HN600_official' 'HNpair_MuMu_WR5000_Zp3000_HN700_official' 'HNpair_MuMu_WR5000_Zp3000_HN800_official' 'HNpair_MuMu_WR5000_Zp3000_HN900_official' 'HNpair_MuMu_WR5000_Zp4000_HN1000_official' 'HNpair_MuMu_WR5000_Zp4000_HN100_official' 'HNpair_MuMu_WR5000_Zp4000_HN1100_official' 'HNpair_MuMu_WR5000_Zp4000_HN1200_official' 'HNpair_MuMu_WR5000_Zp4000_HN1300_official' 'HNpair_MuMu_WR5000_Zp4000_HN1400_official' 'HNpair_MuMu_WR5000_Zp4000_HN1500_official' 'HNpair_MuMu_WR5000_Zp4000_HN1600_official' 'HNpair_MuMu_WR5000_Zp4000_HN1700_official' 'HNpair_MuMu_WR5000_Zp4000_HN1800_official' 'HNpair_MuMu_WR5000_Zp4000_HN1900_official' 'HNpair_MuMu_WR5000_Zp4000_HN200_official' 'HNpair_MuMu_WR5000_Zp4000_HN300_official' 'HNpair_MuMu_WR5000_Zp4000_HN400_official' 'HNpair_MuMu_WR5000_Zp4000_HN500_official' 'HNpair_MuMu_WR5000_Zp4000_HN600_official' 'HNpair_MuMu_WR5000_Zp4000_HN700_official' 'HNpair_MuMu_WR5000_Zp4000_HN800_official' 'HNpair_MuMu_WR5000_Zp4000_HN900_official' 'HNpair_MuMu_WR5000_Zp500_HN100_official' 'HNpair_MuMu_WR5000_Zp500_HN200_official' 'HNpair_MuMu_WR5000_Zp750_HN100_official' 'HNpair_MuMu_WR5000_Zp750_HN200_official' 'HNpair_MuMu_WR5000_Zp750_HN300_official')
declare -a example=("WW" "WZ" "HNMupMup_100" "DYJets")

declare -a tmplist=('ggHtoZZ' 'ggHtoWW' 'vbfHtoWW' 'ggZZto4mu' 'ww_ds' 'WWZ')

declare -a tmpall_mc=('TTJets_aMC' 'WJets' 'WW'  'WZ' 'ZZ' 'DYJets' 'DYJets_10to50' )

declare -a hn=('DYJets_10to50'  'DYJets' 'TT_powheg' 'WJets' 'WGtoLNuG'  'ZGto2LG' 'qcd_15to20_bctoe' 'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-20to30_MuEnriched')
declare -a hn_fake=('DYJets' 'TT_powheg' 'WJets' 'qcd_15to20_bctoe' 'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_DoubleEMEnriched_30-40_mgg80toinf' 'QCD_DoubleEMEnriched_30-inf_mgg40to80' 'QCD_DoubleEMEnriched_40-inf_mgg80toinf')

declare -a pu_dilepton_list=('DYJets_10to50' 'DYJets' 'WJets' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg' 'DYJets_MG_10to50' 'DYJets_MG' 'TTJets_aMC' )

declare -a hntmp=('TTTT' 'TG' 'TTG' 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZToLL_M-10'  'tZq' 'ggHtoWW' 'ggHtoZZ' 'WWG' 'WZG' 'WZto2L2Q_amcatnlo' 'ZZTo2L2Nu_Powheg' 'ZZTo2L2Q_Powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau'  'ggZZto4e' 'ggWWto2L2Nu' 'ww_ds'  )

declare -a hn_eetmp=('DYJets_10to50' 'DYJets' 'WJets' 'WpWpEWK' 'WpWpQCD' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WWW' 'ttW' 'ttZ' 'ttH_nonbb' 'ttH_bb' 'ZZZ' 'WZZ'  'VBF_HToMuMu' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg'  'WWTo2L2Nu' 'WWToLNuQQ' 'QCD_DoubleEMEnriched_30-40_mgg80toinf' 'QCD_DoubleEMEnriched_30-inf_mgg40to80' 'QCD_DoubleEMEnriched_40-inf_mgg80toinf' 'TG' 'TTG' 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZToLL_M-10'  'tZq' 'ggHtoWW' 'ggHtoZZ' 'WWG' 'WZG' 'WZto2L2Q_amcatnlo' 'ZZTo2L2Nu_Powheg' 'ZZTo2L2Q_Powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau'  'ggZZto4e' 'ggWWto2L2Nu' 'ww_ds'  )

declare -a sktmp=('DYJets' 'WJets' 'TT_powheg')
declare -a vv=('WZTo3LNu_powheg' 'ZZTo4L_powheg' 'WGtoLNuG'  'ZGto2LG' 'ZZTo2L2Nu_Powheg' 'ZZTo2L2Q_Powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau'  'ggZZto4e' 'ggWWto2L2Nu' 'ww_ds' 'TG' 'TTG' 'ttWToLNu')

declare -a hn_ee_sig=('WpWpEWK' 'WpWpQCD'  'ZZZ' 'WZZ'  'ww_ds'  'ggZZto4e' 'WZTo3LNu_powheg' 'ZZTo4L_powheg'  'ggHtoZZ'  'WWG' 'WZG'    'ttWToLNu' 'ttZToLL_M-1to10' 'ttZToLL_M-10'  'tZq' 'ggHtoWW' )
declare -a hn_ee_sigcf=('DYJets' 'TT_powheg' 'WWTo2L2Nu' )

declare -a hn_ee_type=('DYJets' 'TT_powheg' 'WJets' 'ZGto2LG' 'QCD_Pt-30to50_EMEnriched')


declare -a tmp=(
'DYJets_MG'
'GG_HToMuMu' 
'SingleTop_t'
'TTJets_aMC'
'TTLL_powheg'
'ttZToLL_M-10'
'WWTo2L2Nu_DS'
'ZGto2LG' )

declare -a tmp2=(
'HNEpMup_100' 
'HNEpMup_1100' 
'HNEpMup_1500' 
'HNEpMup_200' 
'HNEpMup_40' 
'HNEpMup_500' 
'HNEpMup_50' 
'HNEpMup_60' 
'HNMoriondLL_Tchannel_EpMup_100' 
'HNMoriondLL_Tchannel_EpMup_1100' 
'HNMoriondLL_Tchannel_EpMup_200' 
'HNMoriondLL_Tchannel_EpMup_500' 
'HN_MuMuMu_1000'
'HN_MuMuMu_10'
'HN_MuMuMu_150'
'HN_MuMuMu_20'
'HN_MuMuMu_200'
 )
