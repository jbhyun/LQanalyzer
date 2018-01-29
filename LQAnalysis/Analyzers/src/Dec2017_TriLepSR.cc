// $Id: Dec2017_TriLepSR.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQDec2017_TriLepSR Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Dec2017_TriLepSR.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Dec2017_TriLepSR);

 Dec2017_TriLepSR::Dec2017_TriLepSR() : AnalyzerCore(), out_muons(0) {

   SetLogName("Dec2017_TriLepSR");
   Message("In Dec2017_TriLepSR constructor", INFO);
   InitialiseAnalysis();
 }


 void Dec2017_TriLepSR::InitialiseAnalysis() throw( LQError ) {
   
   /// Initialise histograms
   MakeHistograms();  
   // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
   // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
   Message("Feb2016, HwA analysis", INFO);
   return;
 }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Loop///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Dec2017_TriLepSR::ExecuteEvents()throw( LQError ){

////Basic Infos///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Apply the gen weight 
   if(!isData) weight*=MCweight;
   FillHist("GenWeight", MCweight, 1, -2., 2., 4);


   m_logger<<DEBUG<<"RunNumber/Event Number = "<<eventbase->GetEvent().RunNumber()<<" : "<<eventbase->GetEvent().EventNumber()<<LQLogger::endmsg;
   m_logger<<DEBUG<<"isData = "<<isData<<LQLogger::endmsg;


   //Pileup Reweight / Signal Weight
   float pileup_reweight=1., pileup_reweight_systdown=1., pileup_reweight_systup=1.;
   float k_factor_weight=1., geneff_weight=1., gennorm_weight=1.;
   if(!k_isdata){
      pileup_reweight          = eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central); 
      pileup_reweight_systup   = pileup_reweight>0.? eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::up)  /pileup_reweight :0. ;
      pileup_reweight_systdown = pileup_reweight>0.? eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::down)/pileup_reweight :0. ;
      k_factor_weight = GetKFactor();
      geneff_weight   = GenFilterEfficiency(k_sample_name);
      gennorm_weight  = SignalNorm(k_sample_name, 20.);
   }
   FillHist("Basic_PURW", pileup_reweight, 1., 0., 20., 200);

 
   //Total Event(MCweight+PUrw)////////////////
   FillCutFlow("NoCut", weight*pileup_reweight);


   bool EMuMu=false, TriMu=false, MultiLep=false;
   bool CutOpt=false, ObsExpComp=false, SRYield=false, MAWinOpt=false;
   bool SRDist=false;
   bool SystRun=false, TestRun=false;
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("EMuMu"))       EMuMu      = true;
     else if(k_flags.at(i).Contains("TriMu"))       TriMu      = true;
     else if(k_flags.at(i).Contains("MultiLep"))    MultiLep   = true;
     else if(k_flags.at(i).Contains("SystRun"))     SystRun    = true;
     else if(k_flags.at(i).Contains("CutOpt"))      CutOpt     = true;
     else if(k_flags.at(i).Contains("MAWinOpt"))    MAWinOpt   = true;
     else if(k_flags.at(i).Contains("ObsExpComp"))  ObsExpComp = true;
     else if(k_flags.at(i).Contains("SRYield"))     SRYield    = true;
     else if(k_flags.at(i).Contains("SRDist"))      SRDist     = true;
     else if(k_flags.at(i).Contains("TestRun"))     TestRun    = true;
   }

    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //Trigger Path of Analysis
   bool Pass_Trigger=false, Pass_EMuTrigger=false, Pass_DiMuTrigger=false;
   float trigger_ps_weight=1., trigger_period_weight=1.;
   if(EMuMu){
     int Pass_Trigger1=0, Pass_Trigger2=0;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ) Pass_Trigger1++;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger2++;

     if(isData){
       int DataPeriod=GetDataPeriod();
       if( DataPeriod>0 && DataPeriod<7 && Pass_Trigger1==1 ) Pass_Trigger=true;
       else if( DataPeriod==7 && Pass_Trigger2==1 ) Pass_Trigger=true;
     }
     else{
       if( Pass_Trigger1>0 || Pass_Trigger2>0 ) Pass_Trigger=true;
       trigger_period_weight=(Pass_Trigger1*27.257618+Pass_Trigger2*8.605696)/35.863314;
       trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
     }
   }
   if(TriMu){

     if     ( PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")   ) Pass_Trigger=true;
     else if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ) Pass_Trigger=true;

     if(!isData) trigger_ps_weight=WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi);

   }
   if(MultiLep){

     //DiMu
     if     ( PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")   ) Pass_DiMuTrigger=true;
     else if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ) Pass_DiMuTrigger=true;

     //EMu
     int Pass_Trigger1=0, Pass_Trigger2=0;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") )    Pass_Trigger1++;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger2++;

     if(isData){
       int DataPeriod=GetDataPeriod();
       if( DataPeriod>0 && DataPeriod<7 && Pass_Trigger1==1 ) Pass_EMuTrigger=true;
       else if( DataPeriod==7 && Pass_Trigger2==1 )           Pass_EMuTrigger=true;
     }
     else{
       if( Pass_Trigger1>0 || Pass_Trigger2>0 ) Pass_EMuTrigger=true;
       trigger_period_weight=(Pass_Trigger1*27.257618+Pass_Trigger2*8.605696)/35.863314;
       trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
     }
     if(Pass_DiMuTrigger || Pass_EMuTrigger) Pass_Trigger=true;

   }
   FillHist("Basic_TrigPSWeight", trigger_ps_weight, 1., 0., 1., 100);
   weight*=trigger_ps_weight;
   if(!MultiLep) weight*=trigger_period_weight;

   //Trigger Cut
   if(!Pass_Trigger) return;
   FillCutFlow("TriggerCut", weight*pileup_reweight);
   /**********************************************************************************/

   //METFilterCut : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
   if(!PassMETFilter()) return;  FillCutFlow("EventCut", weight*pileup_reweight);

   //Vertex Cut
   //(vtx.ndof>4&&maxAbsZ<=0)||abs(vtx.z)<= 24)&&((maxd0 <=0)||abs(vtx.position.rho)<=2)&&!(vtx.isFake))
   if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; FillCutFlow("VertexCut", weight*pileup_reweight);
   


   //===========================================================================//
   // Objects Selection in Analysis                                             //
   //===========================================================================//

  
   //**PreSelCut*******************************************************************************************//
   //Intended for Code speed boosting up.
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, true);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
     if(EMuMu)      { if( !(electronPreColl.size()>=1 && muonPreColl.size()>=2) ) return; }
     else if(TriMu) { if( !(muonPreColl.size()>=2) ) return; }
     FillCutFlow("PreSel", weight);
   //******************************************************************************************************//

   //Primary Object Collection
   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);

//     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
//     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
//     eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(3.);
//     eventbase->GetMuonSel()->SetBSdz(0.1);
//     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
//   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);
//
//     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
//     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
//     eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(3.);
//     eventbase->GetMuonSel()->SetBSdz(0.1);
//     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");  eventbase->GetMuonSel()->SetRelIso(0.1);
//   std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,true);
   std::vector<snu::KMuon> muonLooseColl, muonTightColl;
     for(int i=0; i<muonPreColl.size(); i++){
       if(PassIDCriteria(muonPreColl.at(i),"Test_POGLIsop4IPp5p1Chi100"))      muonLooseColl.push_back(muonPreColl.at(i));
       if(PassIDCriteria(muonPreColl.at(i),"Test_POGTIsop20IPp01p05sig4Chi4")) muonTightColl.push_back(muonPreColl.at(i));
     }
   std::vector<snu::KMuon> muonColl;  if(k_running_nonprompt){ muonColl=muonLooseColl;} else{ muonColl=muonTightColl;}


//     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_HctoWA_FAKELOOSE);
//     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
//     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
//     eventbase->GetElectronSel()->SetBETrRegIncl(false);
//     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.4, 0.4);
//     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
//     eventbase->GetElectronSel()->SetdxySigMax(4.);
//     eventbase->GetElectronSel()->SetApplyConvVeto(true);
//   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
//
//     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
//     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
//     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
//     eventbase->GetElectronSel()->SetBETrRegIncl(false);
//     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
//     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
//     eventbase->GetElectronSel()->SetdxySigMax(4.);
//     eventbase->GetElectronSel()->SetApplyConvVeto(true);
//   std::vector<snu::KElectron> electronTightColl; eventbase->GetElectronSel()->Selection(electronTightColl);

   
     std::vector<snu::KElectron> electronLooseColl, electronTightColl;
       for(int i=0; i<electronPreColl.size(); i++){
        //if(PassIDCriteria(electronPreColl.at(i), "LMVA06v1Isop4IPp5p1")) electronLooseColl.push_back(electronPreColl.at(i));
        if(PassIDCriteria(electronPreColl.at(i), "LMVA06Isop4IPp025p05sig4")) electronLooseColl.push_back(electronPreColl.at(i));
        if(PassIDCriteria(electronPreColl.at(i), "POGWP90Isop06IPp025p05sig4")) electronTightColl.push_back(electronPreColl.at(i));
       }
   std::vector<snu::KElectron> electronColl;  if(!k_running_nonprompt){ electronColl=electronTightColl;} else{ electronColl=electronLooseColl;}


     bool LeptonVeto=true;
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonLooseColl, electronLooseColl);

   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");
   std::vector<snu::KJet> ljetColl = SelLightJets(jetColl, "Medium");

   float met    = eventbase->GetEvent().MET(snu::KEvent::pfmet);
   float metphi = eventbase->GetEvent().METPhi();
   float met_x  = eventbase->GetEvent().PFMETx();
   float met_y  = eventbase->GetEvent().PFMETy();
   float Pzv, Pzv1, Pzv2;
   snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));
   //snu::KParticle v[4]; v[0].SetPx(met_x); v[0].SetPy(met_y);
   //                     v[1].SetPx(met_x); v[1].SetPy(met_y);

   int nbjets=bjetColl.size(); int njets=jetColl.size(); const int nljets=ljetColl.size();
   int Nvtx=eventbase->GetEvent().nVertices();
   //------------------------------------------------------------------------------------------------------------------//
  

   //=====================================================//
   // Correction Factors                                  //
   //=====================================================//
   float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1., btag_sf=1.;
   float trigger_sf=1.;
   float fake_weight=1.; bool EventCand=false;

   /*This part is for boosting up speed.. SF part takes rather longer time than expected*/
   if     (EMuMu){ if(electronLooseColl.size()>=1 && muonLooseColl.size()>=2) EventCand=true; }
   else if(TriMu){ if(                  muonLooseColl.size()>=3             ) EventCand=true; }


   if(EventCand & !SystRun){
     if(!isData){
       if(EMuMu || TriMu || MultiLep){
 
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electronColl);
    
         trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
         //id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muonColl);
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muonColl);
         iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);

         btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);
  
         //trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );
       }
     }
     else{
       //Perfect Prompt Ratio Approximation applied.(p=1), Applicable to generic number, combination of leptons under premise of high prompt rate.
       if(k_running_nonprompt){
//         fake_weight=GetFakeWeight(muonLooseColl, electronLooseColl, "HNTrilepFakeL2", "HNTrilepTight2", "HctoWAFakeLoose","POGMVAMIP");
       //fake_weight = GetFakeWeight(muonLooseColl, electronLooseColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA06v1Isop4IPp5p1", "POGWP90Isop06IPp025p05sig4", "TrkIsoVVLConeSUSY");

       fake_weight = GetFakeWeight(muonLooseColl, electronLooseColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA06Isop4IPp025p05sig4", "POGWP90Isop06IPp025p05sig4", "TrkIsoVVLConeSUSY");
       }
     }
   }
   weight *= k_factor_weight*geneff_weight*gennorm_weight;
   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*btag_sf*pileup_reweight*fake_weight;
   //-----------------------------------------------------------------------------------------//



   //----------------------------------------------------------------------------------//
   //==================================================================================//
   /////// MAIN ANALYSIS CODE ///////////////////////////////////////////////////////////
   //==================================================================================//
   //----------------------------------------------------------------------------------//

   if(ObsExpComp && EMuMu){

     bool SystRun=false;
     if(!SystRun){
       if( !(electronLooseColl.size()==1 && muonLooseColl.size()==2) ) return;
       if( !(electronColl.size()==1      && muonColl.size()==2) ) return;
       if( !(muonColl.at(0).Charge()!=muonColl.at(1).Charge()) ) return;
       if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>10 && muonColl.at(1).Pt()>10) ) return;
       if( fabs(electronColl.at(0).Eta())>2.5 ) return;
  
       float Mmumu=(muonColl.at(0)+muonColl.at(1)).M();
       float MTW = sqrt(2)*sqrt(met*electronColl.at(0).Pt()-met_x*electronColl.at(0).Px()-met_y*electronColl.at(0).Py());
       float M3l=(electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M();
       if(Mmumu<12) return;
 
 
       bool HasBJet=false, OnZ=false, OnZG=false, WZSel=false, ZGSel=false, ttZSel=false, ZbbSel=false, AN_Sideband=false;
       bool SRincl=false, SRoffZ=false;
       if( bjetColl.size()!=0                )       HasBJet =true;
       if( fabs(Mmumu-91.2)<10               )       OnZ     =true;
       if( fabs(M3l-91.2)<10                 )       OnZG    =true;
       if( OnZ && M3l>101.2 && met>50        )       WZSel   =true;
       if( OnZG && Mmumu<81.2 && met<50      )       ZGSel   =true;
       if( OnZ && HasBJet && jetColl.size()>2)       ttZSel  =true;
       if( OnZ && HasBJet && jetColl.size()<3)       ZbbSel  =true;
       if( Mmumu>40 && HasBJet && jetColl.size()>2 ) AN_Sideband=true;
       if( HasBJet && jetColl.size()>2         )     SRincl  =true;
       if( HasBJet && jetColl.size()>2 && !OnZ )     SRoffZ  =true;
 
       FillHist("YieldComp_CR", 0., weight, 0., 10., 10);
       
  
       //General 3l Selection
       if(!HasBJet){
         FillHist("YieldComp_CR", 1., weight, 0., 10., 10);
         FillHist("PTe_3lOS", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_3lOS", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_3lOS", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_3lOS", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_3lOS", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_3lOS", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_3lOS", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_3lOS", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_3lOS", jetColl.size(), weight, 0., 10., 10);
         FillHist("Nb_3lOS", bjetColl.size(), weight, 0., 10., 10);
         FillHist("MET_3lOS", met, weight, 0., 200., 200);
         FillHist("MTW_3lOS", MTW, weight, 0., 200., 200);
       }
       //3l OnZ ; Fake CR
       if(!HasBJet && OnZ){

         FillHist("YieldComp_CR", 2., weight, 0., 10., 10);
         FillHist("PTe_3lOSOnZ", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_3lOSOnZ", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_3lOSOnZ", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_3lOSOnZ", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_3lOSOnZ", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_3lOSOnZ", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_3lOSOnZ", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_3lOSOnZ", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_3lOSOnZ", jetColl.size(), weight, 0., 10., 10);
         FillHist("Nb_3lOSOnZ", bjetColl.size(), weight, 0., 10., 10);
         FillHist("MET_3lOSOnZ", met, weight, 0., 200., 200);
         FillHist("MTW_3lOSOnZ", MTW, weight, 0., 200., 200);
       }
       //3l OnZ & HasBJet
       if(HasBJet && OnZ){
//         FillHist("YieldComp_CR", 3., weight, 0., 10., 10);
         FillHist("PTe_3lOSOnZHasB", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_3lOSOnZHasB", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_3lOSOnZHasB", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_3lOSOnZHasB", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_3lOSOnZHasB", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_3lOSOnZHasB", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_3lOSOnZHasB", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_3lOSOnZHasB", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_3lOSOnZHasB", jetColl.size(), weight, 0., 10., 10);
         FillHist("Nb_3lOSOnZHasB", bjetColl.size(), weight, 0., 10., 10);
         FillHist("MET_3lOSOnZHasB", met, weight, 0., 200., 200);
         FillHist("MTW_3lOSOnZHasB", MTW, weight, 0., 200., 200);
       }
       if(WZSel){
         //FillHist("YieldComp_CR", 4., weight, 0., 10., 10);
         FillHist("YieldComp_CR", 3., weight, 0., 10., 10);
         FillHist("PTe_WZSel", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_WZSel", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_WZSel", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_WZSel", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_WZSel", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_WZSel", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_WZSel", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_WZSel", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_WZSel", jetColl.size(), weight, 0., 10., 10);
         FillHist("Nb_WZSel", bjetColl.size(), weight, 0., 10., 10);
         FillHist("MET_WZSel", met, weight, 0., 200., 200);
         FillHist("MTW_WZSel", MTW, weight, 0., 200., 200);
       }
       if(ZGSel){
         //FillHist("YieldComp_CR", 5., weight, 0., 10., 10);
         FillHist("YieldComp_CR", 4., weight, 0., 10., 10);
         FillHist("PTe_ZGSel", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_ZGSel", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_ZGSel", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_ZGSel", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_ZGSel", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_ZGSel", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_ZGSel", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_ZGSel", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_ZGSel", jetColl.size(), weight, 0., 10., 10);
         FillHist("Nb_ZGSel", bjetColl.size(), weight, 0., 10., 10);
         FillHist("MET_ZGSel", met, weight, 0., 200., 200);
         FillHist("MTW_ZGSel", MTW, weight, 0., 200., 200);
       }
       if(ttZSel){
         //FillHist("YieldComp_CR", 6., weight, 0., 10., 10);
         FillHist("YieldComp_CR", 5., weight, 0., 10., 10);
         FillHist("PTe_ttZSel", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_ttZSel", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_ttZSel", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_ttZSel", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_ttZSel", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_ttZSel", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_ttZSel", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_ttZSel", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_ttZSel", jetColl.size(), weight, 0., 10., 10);
         FillHist("Nb_ttZSel", bjetColl.size(), weight, 0., 10., 10);
         FillHist("MET_ttZSel", met, weight, 0., 200., 200);
         FillHist("MTW_ttZSel", MTW, weight, 0., 200., 200);
       }
       if(AN_Sideband){
         //FillHist("YieldComp_CR", 7., weight, 0., 10., 10);
         FillHist("YieldComp_CR", 6., weight, 0., 10., 10);
         FillHist("PTe_ANSideband", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_ANSideband", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_ANSideband", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_ANSideband", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_ANSideband", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_ANSideband", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_ANSideband", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_ANSideband", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_ANSideband", jetColl.size(), weight, 0., 10., 10);
         FillHist("Nb_ANSideband", bjetColl.size(), weight, 0., 10., 10);
         FillHist("MET_ANSideband", met, weight, 0., 200., 200);
         FillHist("MTW_ANSideband", MTW, weight, 0., 200., 200);
       }
       if( !isData || k_running_nonprompt ){
         if(SRincl){
           //FillHist("YieldComp_CR", 8., weight, 0., 10., 10);
           FillHist("YieldComp_CR", 7., weight, 0., 10., 10);
           FillHist("PTe_SRincl", electronColl.at(0).Pt(), weight, 0., 200., 200);
           FillHist("Etae_SRincl", electronColl.at(0).Eta(), weight, -5., 5., 100);
           FillHist("PTmu1_SRincl", muonColl.at(0).Pt(), weight, 0., 200., 200);
           FillHist("PTmu2_SRincl", muonColl.at(1).Pt(), weight, 0., 200., 200);
           FillHist("Etamu1_SRincl", muonColl.at(0).Eta(), weight, -5., 5., 100);
           FillHist("Etamu2_SRincl", muonColl.at(1).Eta(), weight, -5., 5., 100);
    
           FillHist("Mmumu_SRincl", Mmumu, weight, 0., 200., 200);
           FillHist("M3l_SRincl", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
           FillHist("Nj_SRincl", jetColl.size(), weight, 0., 10., 10);
           FillHist("Nb_SRincl", bjetColl.size(), weight, 0., 10., 10);
           FillHist("MET_SRincl", met, weight, 0., 200., 200);
           FillHist("MTW_SRincl", MTW, weight, 0., 200., 200);
         } 
         if(SRoffZ){
           //FillHist("YieldComp_CR", 9., weight, 0., 10., 10);
           FillHist("YieldComp_CR", 8., weight, 0., 10., 10);
           FillHist("PTe_SRoffZ", electronColl.at(0).Pt(), weight, 0., 200., 200);
           FillHist("Etae_SRoffZ", electronColl.at(0).Eta(), weight, -5., 5., 100);
           FillHist("PTmu1_SRoffZ", muonColl.at(0).Pt(), weight, 0., 200., 200);
           FillHist("PTmu2_SRoffZ", muonColl.at(1).Pt(), weight, 0., 200., 200);
           FillHist("Etamu1_SRoffZ", muonColl.at(0).Eta(), weight, -5., 5., 100);
           FillHist("Etamu2_SRoffZ", muonColl.at(1).Eta(), weight, -5., 5., 100);
    
           FillHist("Mmumu_SRoffZ", Mmumu, weight, 0., 200., 200);
           FillHist("M3l_SRoffZ", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
           FillHist("Nj_SRoffZ", jetColl.size(), weight, 0., 10., 10);
           FillHist("Nb_SRoffZ", bjetColl.size(), weight, 0., 10., 10);
           FillHist("MET_SRoffZ", met, weight, 0., 200., 200);
           FillHist("MTW_SRoffZ", MTW, weight, 0., 200., 200);
         }
       }
     }//End of Not SystRun
   }//End of Closure
   if(SRDist){
     if(EMuMu && !SystRun){
       CheckSRDist(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMu");
     }
     if(TriMu && !SystRun){
       CheckSRDist(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMu");
     }
   }
   if(SRYield){
     if(EMuMu && (!SystRun)){
        CheckSRYield(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMu");
     }
     if(SystRun){

       // Syst Sel. and Syst Corr. ----------------------------------------------------------------------------------------------------//
         eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
         eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
         eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
         eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
         eventbase->GetMuonSel()->SetdxySigMax(4.);
         eventbase->GetMuonSel()->SetChiNdof(4.);
       std::vector<snu::KMuon> MuTMuEnUpColl; eventbase->GetMuonSel()->Selection(MuTMuEnUpColl,true, "SystUpMuEn");

         eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
         eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
         eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
         eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
         eventbase->GetMuonSel()->SetdxySigMax(4.);
         eventbase->GetMuonSel()->SetChiNdof(4.);
       std::vector<snu::KMuon> MuTMuEnDownColl; eventbase->GetMuonSel()->Selection(MuTMuEnDownColl,true, "SystDownMuEn");

         eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
         eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
         eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
         eventbase->GetElectronSel()->SetBETrRegIncl(false);
         eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
         eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
         eventbase->GetElectronSel()->SetdxySigMax(4.);
         eventbase->GetElectronSel()->SetApplyConvVeto(true);
       std::vector<snu::KElectron> EleTElEnUpColl; eventbase->GetElectronSel()->Selection(EleTElEnUpColl, "SystUpElEn");

         eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
         eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
         eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
         eventbase->GetElectronSel()->SetBETrRegIncl(false);
         eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
         eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025);   eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
         eventbase->GetElectronSel()->SetdxySigMax(4.);
         eventbase->GetElectronSel()->SetApplyConvVeto(true);
       std::vector<snu::KElectron> EleTElEnDownColl; eventbase->GetElectronSel()->Selection(EleTElEnDownColl, "SystDownElEn");

       LeptonVeto=true;
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJESUpColl; eventbase->GetJetSel()->Selection(jetJESUpColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystUpJES");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJERUpColl; eventbase->GetJetSel()->Selection(jetJERUpColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystUpJER");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJESDownColl; eventbase->GetJetSel()->Selection(jetJESDownColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystDownJES");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJERDownColl; eventbase->GetJetSel()->Selection(jetJERDownColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystDownJER");

       std::vector<snu::KJet> bjetJESUpColl   = SelBJets(jetJESUpColl  , "Medium");
       std::vector<snu::KJet> bjetJERUpColl   = SelBJets(jetJERUpColl  , "Medium");
       std::vector<snu::KJet> bjetJESDownColl = SelBJets(jetJESDownColl, "Medium");
       std::vector<snu::KJet> bjetJERDownColl = SelBJets(jetJERDownColl, "Medium");

       float met_JESup     = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetEn       , snu::KEvent::up);
       float met_JERup     = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetRes      , snu::KEvent::up);
       float met_Unclup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::Unclustered , snu::KEvent::up);
       float met_ElEnup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::ElectronEn  , snu::KEvent::up);
       float met_MuEnup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::MuonEn      , snu::KEvent::up);
       float met_JESdown   = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetEn       , snu::KEvent::down);
       float met_JERdown   = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetRes      , snu::KEvent::down);
       float met_Uncldown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::Unclustered , snu::KEvent::down);
       float met_ElEndown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::ElectronEn  , snu::KEvent::down);
       float met_MuEndown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::MuonEn      , snu::KEvent::down);

       //Scale Factors--------------------------------------------------------------------------------//      
       float systweight=weight; //Lumi weight+Prescale weight+Period weight+PU weight applied by default;
       float id_weight_ele_ElEnup=1.  , reco_weight_ele_ElEnup=1.  , id_weight_mu_MuEnup=1.  , trk_weight_mu_MuEnup=1.;
       float id_weight_ele_ElEndown=1., reco_weight_ele_ElEndown=1., id_weight_mu_MuEndown=1., trk_weight_mu_MuEndown=1.;
       float btag_sf_JESup=1.  , btag_sf_JERup=1.  , btag_sf_LTagup=1.  , btag_sf_BCTagup=1.;
       float btag_sf_JESdown=1., btag_sf_JERdown=1., btag_sf_LTagdown=1., btag_sf_BCTagdown=1.;
       float fake_weight_FRup=1., fake_weight_FRdown=1.;
       if(!isData){
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electronColl);
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muonColl)*mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
         trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
         btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

         id_weight_ele_ElEnup   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", EleTElEnUpColl);
         reco_weight_ele_ElEnup = mcdata_correction->ElectronRecoScaleFactor(EleTElEnUpColl);
         id_weight_mu_MuEnup    = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", MuTMuEnUpColl)*mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", MuTMuEnUpColl);
         trk_weight_mu_MuEnup   = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnUpColl);

         id_weight_ele_ElEndown   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", EleTElEnDownColl);
         reco_weight_ele_ElEndown = mcdata_correction->ElectronRecoScaleFactor(EleTElEnDownColl);
         id_weight_mu_MuEndown    = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", MuTMuEnDownColl)*mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", MuTMuEnDownColl);
         trk_weight_mu_MuEndown   = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnDownColl);

         btag_sf_LTagup = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpLTag");
         btag_sf_BCTagup= BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpBCTag");
         btag_sf_JESup  = BTagScaleFactor_1a(jetJESUpColl, snu::KJet::CSVv2, snu::KJet::Medium);
         btag_sf_JERup  = BTagScaleFactor_1a(jetJERUpColl, snu::KJet::CSVv2, snu::KJet::Medium);

         btag_sf_LTagdown = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownLTag");
         btag_sf_BCTagdown= BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownBCTag");
         btag_sf_JESdown  = BTagScaleFactor_1a(jetJESDownColl, snu::KJet::CSVv2, snu::KJet::Medium);
         btag_sf_JERdown  = BTagScaleFactor_1a(jetJERDownColl, snu::KJet::CSVv2, snu::KJet::Medium);
       }
       else if(k_running_nonprompt){
         fake_weight        = GetFakeWeight(muonLooseColl, electronLooseColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA06Isop4IPp025p05sig4", "POGWP90Isop06IPp025p05sig4", "TrkIsoVVLConeSUSY");
         fake_weight_FRup   = GetFakeWeight(muonLooseColl, electronLooseColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA06Isop4IPp025p05sig4", "POGWP90Isop06IPp025p05sig4", "TrkIsoVVLConeSUSYSystUp");
         fake_weight_FRdown = GetFakeWeight(muonLooseColl, electronLooseColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA06Isop4IPp025p05sig4", "POGWP90Isop06IPp025p05sig4", "TrkIsoVVLConeSUSYSystDown");
         //fake_weight = GetFakeWeight(muonLooseColl, electronLooseColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA06v1Isop4IPp5p1", "POGWP90Isop06IPp025p05sig4", "TrkIsoVVLConeSUSY");
       }
       float systweight_central=weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight;
       float systweight_ElEnup =weight*id_weight_ele_ElEnup*reco_weight_ele_ElEnup*id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight;
       float systweight_MuEnup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu_MuEnup*trk_weight_mu_MuEnup*btag_sf        *fake_weight;
       float systweight_JESup  =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JESup  *fake_weight;
       float systweight_JERup  =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JERup  *fake_weight;
       float systweight_LTagup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_LTagup *fake_weight;
       float systweight_BCTagup=weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_BCTagup*fake_weight;
       float systweight_PUup   =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight*pileup_reweight_systup;
       float systweight_FRup   =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight_FRup;

       float systweight_ElEndown =weight*id_weight_ele_ElEndown*reco_weight_ele_ElEndown*id_weight_mu         *trk_weight_mu         *btag_sf          *fake_weight;
       float systweight_MuEndown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu_MuEndown*trk_weight_mu_MuEndown*btag_sf          *fake_weight;
       float systweight_JESdown  =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JESdown  *fake_weight;
       float systweight_JERdown  =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JERdown  *fake_weight;
       float systweight_LTagdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_LTagdown *fake_weight;
       float systweight_BCTagdown=weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_BCTagdown*fake_weight;
       float systweight_PUdown   =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *fake_weight*pileup_reweight_systdown;
       float systweight_FRdown   =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *fake_weight_FRdown;



       //----------------------------------------------------------------------------------------------------------------------//

       TString ChannelString = EMuMu? "EMuMu":"TriMu";
       DoSystRun("SRYield", ChannelString+"",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_central);
       if( !isData ){
         DoSystRun("SRYield", ChannelString+"SystUpPU",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_PUup);
         DoSystRun("SRYield", ChannelString+"SystUpUncl",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi),
                   systweight_central);
         DoSystRun("SRYield", ChannelString+"SystUpJES",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESUpColl, bjetJESUpColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi),
                   systweight_JESup);
         DoSystRun("SRYield", ChannelString+"SystUpJER",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERUpColl, bjetJERUpColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi),
                   systweight_JERup);
         DoSystRun("SRYield", ChannelString+"SystUpBTag_L",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_LTagup);
         DoSystRun("SRYield", ChannelString+"SystUpBTag_BC",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_BCTagup);
         DoSystRun("SRYield", ChannelString+"SystUpElEn",
                   EleTElEnUpColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi),
                   systweight_ElEnup);
         DoSystRun("SRYield", ChannelString+"SystUpMuEn",
                   electronColl, electronLooseColl, MuTMuEnUpColl, muonLooseColl, jetColl, bjetColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi),
                   systweight_MuEnup);


         DoSystRun("SRYield", ChannelString+"SystDownPU",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_PUdown);
         DoSystRun("SRYield", ChannelString+"SystDownUncl",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi),
                   systweight_central);
         DoSystRun("SRYield", ChannelString+"SystDownJES",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESDownColl, bjetJESDownColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi),
                   systweight_JESdown);
         DoSystRun("SRYield", ChannelString+"SystDownJER",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERDownColl, bjetJERDownColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi),
                   systweight_JERdown);
         DoSystRun("SRYield", ChannelString+"SystDownBTag_L",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_LTagdown);
         DoSystRun("SRYield", ChannelString+"SystDownBTag_BC",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_BCTagdown);
         DoSystRun("SRYield", ChannelString+"SystDownElEn",
                   EleTElEnDownColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi),
                   systweight_ElEndown);
         DoSystRun("SRYield", ChannelString+"SystDownMuEn",
                   electronColl, electronLooseColl, MuTMuEnDownColl, muonLooseColl, jetColl, bjetColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi),
                   systweight_MuEndown);

       }
       else if( isData && k_running_nonprompt ){
         DoSystRun("SRYield", ChannelString+"SystUpFR",
                  electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                  systweight_FRup);
         DoSystRun("SRYield", ChannelString+"SystDownFR",
                  electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                  systweight_FRdown);
       }
     }//End of SystRun
   }
   if(CutOpt){

     if(MAWinOpt){
       if(EMuMu){
         if( !(electronLooseColl.size()==1 && muonLooseColl.size()==2) ) return;
         if( !(electronColl.size()==1      && muonColl.size()==2) ) return;
         if( !(muonColl.at(0).Charge()!=muonColl.at(1).Charge()) ) return;
         if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>10 && muonColl.at(1).Pt()>10) ) return;
         if( fabs(electronColl.at(0).Eta())>2.5 ) return;
  
         float Mmumu=(muonColl.at(0)+muonColl.at(1)).M();
         if(Mmumu<12) return;
  
         if(bjetColl.size()==0) return;
         if(jetColl.size()<2) return;
  
         float MAhiggs=GetHiggsMass(k_sample_name, "A");
         const int Nwindow=55;
         float MwidthCuts[Nwindow]={10., 9., 8., 7.,
                                    6.0, 5.9, 5.8, 5.7, 5.6, 5.5, 5.4, 5.3, 5.2, 5.1,
                                    5.0, 4.9, 4.8, 4.7, 4.6, 4.5, 4.4, 4.3, 4.2, 4.1,
                                    4.0, 3.9, 3.8, 3.7, 3.6, 3.5, 3.4, 3.3, 3.2, 3.1,
                                    3.0, 2.9, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 
                                    2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0};
         
         FillHist("NCount_Win", 0., weight, 0., 20., 200);
         for(int it_win=0; it_win<Nwindow; it_win++){
           if(fabs(Mmumu-MAhiggs)<MwidthCuts[it_win]){
             FillHist("NCount_Win", MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
           }
         }
       }
     }
     else{

     if(EMuMu){
       if( !(electronLooseColl.size()==1 && muonLooseColl.size()==2) ) return;
       if( !(electronColl.size()==1      && muonColl.size()==2) ) return;
       if( !(muonColl.at(0).Charge()!=muonColl.at(1).Charge()) ) return;
       if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>10 && muonColl.at(1).Pt()>5) ) return;
       if( fabs(electronColl.at(0).Eta())>2.5 ) return;
  
       float Mmumu=(muonColl.at(0)+muonColl.at(1)).M();
       float MTW = sqrt(2)*sqrt(met*electronColl.at(0).Pt()-met_x*electronColl.at(0).Px()-met_y*electronColl.at(0).Py());
       float M3l=(electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M();
       if(Mmumu<12) return;

       const int NJetCut=2;     int JetCuts[NJetCut]     ={2,3};
       const int NBJetCut=1;    int BJetCuts[NBJetCut]   ={1};
       const int NZVetoCut=1;   bool ZVetoCuts[NZVetoCut]={true};

       const int NPtElCuts=2;   float ElPtCuts[NPtElCuts]    ={25., 30.};
       const int NPtMu1Cuts=3;  float Mu1PtCuts[NPtMu1Cuts]  ={10., 15., 20.};
       const int NPtMu2Cuts=5;  float Mu2PtCuts[NPtMu2Cuts]  ={5., 7., 10., 12., 15.};
       const int NPtJetCuts=3;  float JetPtCuts[NPtJetCuts]  ={25., 30., 35.};
       const int NMETCuts=3;    float METCuts[NMETCuts]      ={0., 20., 30.};

                                                     
       FillHist("Count_PreSel", 0., weight, 0., 1., 1);       
       for(int it_Z=0; it_Z<NZVetoCut; it_Z++){
         if(ZVetoCuts[it_Z] && fabs(Mmumu-91.2)<10) continue;
         FillHist("Count_ZVeto", it_Z, weight, 0., 2., 2);
       }
       for(int it_nb=0; it_nb<NBJetCut; it_nb++){
         if(bjetColl.size()<BJetCuts[it_nb]) continue;
         FillHist("Count_NBCut", it_nb, weight, 0., 2., 2);
       }
       for(int it_nj=0; it_nj<NJetCut; it_nj++){
         if(jetColl.size()<JetCuts[it_nj]) continue;
         FillHist("Count_NJCut", it_nj, weight, 0., 4., 4);
       }
       for(int it_pte=0; it_pte<NPtElCuts; it_pte++){
         if( electronColl.at(0).Pt()<ElPtCuts[it_pte] ) continue;       
         FillHist("Count_PtElCut", it_pte, weight, 0., 2., 2);
       }
       for(int it_ptmu1=0; it_ptmu1<NPtMu1Cuts; it_ptmu1++){
         if( muonColl.at(0).Pt()<Mu1PtCuts[it_ptmu1] ) continue;
         FillHist("Count_PtMu1Cut", it_ptmu1, weight, 0., 3., 3);
       }
       for(int it_ptmu2=0; it_ptmu2<NPtMu2Cuts; it_ptmu2++){
         if( muonColl.at(1).Pt()<Mu2PtCuts[it_ptmu2] ) continue;
         FillHist("Count_PtMu2Cut", it_ptmu2, weight, 0., 5., 5);
       }
       for(int it_ptj=0; it_ptj<NPtJetCuts; it_ptj++){
         if(jetColl.size()>1 && jetColl.at(1).Pt()<JetPtCuts[it_ptj]) continue;
         FillHist("Count_PtjCut", it_ptj, weight, 0., 4., 4);
       }


       int it_Sel=0;
       for(int it_Z=0 ;    it_Z<NZVetoCut;        it_Z++){//1
       for(int it_nb=0;    it_nb<NBJetCut;        it_nb++){//2
       for(int it_nj=0;    it_nj<NJetCut ;        it_nj++){//3
       for(int it_pte=0;   it_pte<NPtElCuts;      it_pte++){//4
       for(int it_ptmu1=0; it_ptmu1<NPtMu1Cuts;   it_ptmu1++){//5
       for(int it_ptmu2=0; it_ptmu2<NPtMu2Cuts;   it_ptmu2++){//6
       for(int it_ptj=0;   it_ptj<NPtJetCuts;     it_ptj++){//7
       for(int it_met=0;   it_met<NMETCuts;       it_met++){//8
         it_Sel++;

         //PassCuts----------------------------------------------//
         if( !(electronColl.at(0).Pt()>ElPtCuts[it_pte] && muonColl.at(0).Pt()>Mu1PtCuts[it_ptmu1] && muonColl.at(1).Pt()>Mu2PtCuts[it_ptmu2]) ) continue;
         if( Mu1PtCuts[it_ptmu1]<Mu2PtCuts[it_ptmu2] ) continue;
         if( ZVetoCuts[it_Z] && fabs(Mmumu-91.2)<10 ) continue;
         if( jetColl.size()<JetCuts[it_nj] || bjetColl.size()<BJetCuts[it_nb] ) continue;
         if( jetColl.at(JetCuts[it_nj]-1).Pt()<JetPtCuts[it_ptj]   ) continue;
         if( bjetColl.at(BJetCuts[it_nb]-1).Pt()<JetPtCuts[it_ptj] ) continue;
         if( met<METCuts[it_met] ) continue;
         //------------------------------------------------------//

         if(Mmumu>12 && Mmumu<40){
           FillHist("Count_AllCut", it_Sel-1, weight, 0., 10000., 10000);
         }
       }//8
       }//7
       }//6
       }//5
       }//4
       }//3
       }//2
       }//1
       
     }//End of EMuMu
     if(TriMu){
       if( !(muonLooseColl.size()==3 && electronLooseColl.size()==0) ) return;
       if( !(muonColl.size()==3) ) return;
       if( fabs(SumCharge(muonColl))!=1 ) return;
       if( !(muonColl.at(0).Pt()>20 && muonColl.at(1).Pt()>10 && muonColl.at(2).Pt()>5) ) return;

       int   IdxOS  = TriMuChargeIndex(muonColl, "OS");
       int   IdxSS1 = TriMuChargeIndex(muonColl, "SS1");
       int   IdxSS2 = TriMuChargeIndex(muonColl, "SS2");
       float MOSSS1 = (muonColl.at(IdxOS)+muonColl.at(IdxSS1)).M();
       float MOSSS2 = (muonColl.at(IdxOS)+muonColl.at(IdxSS2)).M();
       if( !(MOSSS1>12 && MOSSS2>12) ) return;


       const int NJetCut=2;     int JetCuts[NJetCut]     ={2,3};
       const int NBJetCut=1;    int BJetCuts[NBJetCut]   ={1};
       const int NZVetoCut=2;   bool ZVetoCuts[NZVetoCut]={true, false};

       const int NPtMu1Cuts=3;  float Mu1PtCuts[NPtMu1Cuts]={20., 25., 30.};
       const int NPtMu2Cuts=3;  float Mu2PtCuts[NPtMu2Cuts]={10., 12., 15.};
       const int NPtMu3Cuts=5;  float Mu3PtCuts[NPtMu3Cuts]={5., 7., 10., 12., 15.};
       const int NPtJetCuts=3;  float JetPtCuts[NPtJetCuts]={25., 30., 35.};
       const int NMETCuts=3;    float METCuts[NMETCuts]    ={0., 20., 30.};


       FillHist("Count_PreSel", 0., weight, 0., 1., 1);       

       for(int it_Z=0; it_Z<NZVetoCut; it_Z++){
         if(ZVetoCuts[it_Z] && (fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10)) continue;
         FillHist("Count_ZVeto", it_Z, weight, 0., 2., 2);
       }
       for(int it_nb=0; it_nb<NBJetCut; it_nb++){
         if(bjetColl.size()<BJetCuts[it_nb]) continue;
         FillHist("Count_NBCut", it_nb, weight, 0., 2., 2);
       }
       for(int it_nj=0; it_nj<NJetCut; it_nj++){
         if(jetColl.size()<JetCuts[it_nj]) continue;
         FillHist("Count_NJCut", it_nj, weight, 0., 4., 4);
       }
       for(int it_ptmu1=0; it_ptmu1<NPtMu1Cuts; it_ptmu1++){
         if( muonColl.at(0).Pt()<Mu1PtCuts[it_ptmu1] ) continue;
         FillHist("Count_PtMu1Cut", it_ptmu1, weight, 0., 3., 3);
       }
       for(int it_ptmu2=0; it_ptmu2<NPtMu2Cuts; it_ptmu2++){
         if( muonColl.at(1).Pt()<Mu2PtCuts[it_ptmu2] ) continue;
         FillHist("Count_PtMu2Cut", it_ptmu2, weight, 0., 5., 5);
       }
       for(int it_ptmu3=0; it_ptmu3<NPtMu3Cuts; it_ptmu3++){
         if( muonColl.at(2).Pt()<Mu3PtCuts[it_ptmu3] ) continue;       
         FillHist("Count_PtMu3Cut", it_ptmu3, weight, 0., 2., 2);
       }
       for(int it_ptj=0; it_ptj<NPtJetCuts; it_ptj++){
         if(jetColl.size()>1 && jetColl.at(1).Pt()<JetPtCuts[it_ptj]) continue;
         FillHist("Count_PtjCut", it_ptj, weight, 0., 4., 4);
       }


       int it_Sel=0;
       for(int it_Z=0 ;    it_Z<NZVetoCut;      it_Z++){//1
       for(int it_nb=0;    it_nb<NBJetCut;      it_nb++){//2
       for(int it_nj=0;    it_nj<NJetCut ;      it_nj++){//3
       for(int it_ptmu1=0; it_ptmu1<NPtMu1Cuts; it_ptmu1++){//4
       for(int it_ptmu2=0; it_ptmu2<NPtMu2Cuts; it_ptmu2++){//5
       for(int it_ptmu3=0; it_ptmu3<NPtMu3Cuts; it_ptmu3++){//6
       for(int it_ptj=0;   it_ptj<NPtJetCuts;   it_ptj++){//7
       for(int it_met=0;   it_met<NMETCuts;     it_met++){//8
         it_Sel++;


         //PassCuts----------------------------------------------//
         if( !(muonColl.at(0).Pt()>Mu1PtCuts[it_ptmu1] && muonColl.at(1).Pt()>Mu2PtCuts[it_ptmu2] && muonColl.at(2).Pt()>Mu3PtCuts[it_ptmu3]) ) continue;
         if( Mu1PtCuts[it_ptmu1]<Mu2PtCuts[it_ptmu2] || Mu2PtCuts[it_ptmu2]<Mu3PtCuts[it_ptmu3] ) continue;
         if( ZVetoCuts[it_Z] && (fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10) ) continue;
         if( jetColl.size()<JetCuts[it_nj] || bjetColl.size()<BJetCuts[it_nb] ) continue;
         if( jetColl.at(JetCuts[it_nj]-1).Pt()<JetPtCuts[it_ptj] ) continue;
         if( bjetColl.at(BJetCuts[it_nb]-1).Pt()<JetPtCuts[it_ptj] ) continue;
         if( met<METCuts[it_met] ) continue;
         //------------------------------------------------------//

         if(MOSSS1>12 && MOSSS1<40){
           FillHist("Count_AllCut1", it_Sel-1, weight, 0., 10000., 10000);
         }
         if(MOSSS2>12 && MOSSS2<40){
           FillHist("Count_AllCut2", it_Sel-1, weight, 0., 10000., 10000);
         }

       }//8
       }//7
       }//6
       }//5
       }//4
       }//3
       }//2
       }//1
 
     }//End of TriMu
     }//End of Not Mass Window Opt
   }//End of CutOpt
/////////////////////////////////////////////////////////////////////////////////// 

return;

}// End of execute event loop
  


void Dec2017_TriLepSR::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Dec2017_TriLepSR::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
//  string analysisdir = getenv("FILEDIR");  
//  if(!k_isdata) reweightPU = new Reweight((analysisdir + "SNUCAT_Pileup.root").c_str());

  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
//  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
//  DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

Dec2017_TriLepSR::~Dec2017_TriLepSR() {
  
  Message("In Dec2017_TriLepSR Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}



void Dec2017_TriLepSR::CheckSRYield(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");

  if(EMuMu){
    if( !(EleLColl.size()==1 && MuLColl.size()==2) ) return;
    if( !(EleTColl.size()==1 && MuTColl.size()==2) ) return;
    if( !(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) ) return;
    if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && MuTColl.at(1).Pt()>10) ) return;
    if( fabs(EleTColl.at(0).Eta())>2.5 ) return;

    float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
    float MTW = sqrt(2)*sqrt(MET*EleTColl.at(0).Pt()-METx*EleTColl.at(0).Px()-METy*EleTColl.at(0).Py());
    float M3l=(EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M();
    if(Mmumu<12) return;

    //if(Mmumu>40) return;
    if(fabs(Mmumu-91.2)<10) return;
    if(BJetColl.size()==0) return;

    FillHist("Nj_PreSel_1e2mu"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Mmumu_PreSel_1e2mu"+Label, Mmumu, weight, 0., 200., 40);
    //c.f. Regarding detector resolution not considering FSR tail
    //σ(m, pT) = 0.026 + 0.0065m for the barrel,
    //σ(m, pT) = 0.026 + 0.013m GeV/c2 for the endcap.
    if(JetColl.size()>1){
      FillHist("Mmumu_2jCut_1e2mu"+Label, Mmumu, weight, 0., 200., 40);

      if( Mmumu>12 && Mmumu<40 ) FillHist("Count_2jCut_1e2mu_M12to40"+Label, 0., weight, 0., 1., 1);
      if( fabs(Mmumu-15)<5     ) FillHist("Count_2jCut_1e2mu_M15"+Label, 0., weight, 0., 1., 1);
      if( fabs(Mmumu-20)<5     ) FillHist("Count_2jCut_1e2mu_M20"+Label, 0., weight, 0., 1., 1);
      if( fabs(Mmumu-25)<5     ) FillHist("Count_2jCut_1e2mu_M25"+Label, 0., weight, 0., 1., 1);
      if( fabs(Mmumu-30)<5     ) FillHist("Count_2jCut_1e2mu_M30"+Label, 0., weight, 0., 1., 1);
      if( fabs(Mmumu-35)<5     ) FillHist("Count_2jCut_1e2mu_M35"+Label, 0., weight, 0., 1., 1);
    }
  }
  if(TriMu){
    if( !(MuLColl.size()==3 && EleLColl.size()==0) ) return;
    if( !(MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 && MuTColl.at(2).Pt()>10) ) return;

    int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    if( !(MOSSS1>12 && MOSSS2>12) ) return;
    if( fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10 ) return;
    if( BJetColl.size()==0 ) return;


    if(JetColl.size()>1){
      FillHist("MOSSS1_2j"+Label, MOSSS1, weight, 0., 200., 40);
      FillHist("MOSSS2_2j"+Label, MOSSS2, weight, 0., 200., 40);
      if(MOSSS2<40)         FillHist("Count_2jCut_MOSSS2_3mu_M12to40"+Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-15)<5) FillHist("Count_2jCut_MOSSS2_3mu_M15"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-20)<5) FillHist("Count_2jCut_MOSSS2_3mu_M20"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-25)<5) FillHist("Count_2jCut_MOSSS2_3mu_M25"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-30)<5) FillHist("Count_2jCut_MOSSS2_3mu_M30"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-35)<5) FillHist("Count_2jCut_MOSSS2_3mu_M35"    +Label, 0., weight, 0., 1., 1);
    }
    if(JetColl.size()>2){
      FillHist("MOSSS1_3j"+Label, MOSSS1, weight, 0., 200., 40);
      FillHist("MOSSS2_3j"+Label, MOSSS2, weight, 0., 200., 40);
      if(MOSSS2<40)           FillHist("Count_3jCut_MOSSS2_3mu_M12to40"+Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-15)<5  ) FillHist("Count_3jCut_MOSSS2_3mu_M15"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-20)<5  ) FillHist("Count_3jCut_MOSSS2_3mu_M20"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-25)<5  ) FillHist("Count_3jCut_MOSSS2_3mu_M25"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-30)<5  ) FillHist("Count_3jCut_MOSSS2_3mu_M30"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-35)<5  ) FillHist("Count_3jCut_MOSSS2_3mu_M35"    +Label, 0., weight, 0., 1., 1);
    }
  }

}



void Dec2017_TriLepSR::CheckSRDist(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");

  if(EMuMu){
    if( !(EleLColl.size()==1 && MuLColl.size()==2) ) return;
    if( !(EleTColl.size()==1 && MuTColl.size()==2) ) return;
    if( !(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) ) return;
    if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && MuTColl.at(1).Pt()>10) ) return;
    if( fabs(EleTColl.at(0).Eta())>2.5 ) return;

    float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
    float MTW = sqrt(2)*sqrt(MET*EleTColl.at(0).Pt()-METx*EleTColl.at(0).Px()-METy*EleTColl.at(0).Py());
    float M3l=(EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M();
    if(Mmumu<12) return;

    //if(Mmumu>40) return;
    if(fabs(Mmumu-91.2)<10) return;
    if(BJetColl.size()==0) return;
    if(JetColl.size()>=1){
      FillHist("Nj_SR1j", JetColl.size(), weight, 0., 10., 10);
    }

    int Idxj1_W=-1, Idxj2_W=-1;
    snu::KParticle v; v.SetPxPyPzE(METx, METy, 0, sqrt(METx*METx+METy*METy));
    float Pzv1=GetvPz(v, EleTColl.at(0), 1), Pzv2=GetvPz(v, EleTColl.at(0), 2);
    float Pzv=fabs(Pzv1)<fabs(Pzv2)? Pzv1:Pzv2;
    v.SetXYZM(METx,METy,Pzv,0.);

    float M3lv=(EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)+v).M(), M2l2j=0.; 
    if(JetColl.size()>=2){
      Idxj1_W=GetDaughterCandIdx(JetColl, "W", -1., "Lead");
      Idxj2_W=GetDaughterCandIdx(JetColl, "W", -1., "Subl");
      M2l2j=(MuTColl.at(0)+MuTColl.at(1)+JetColl.at(Idxj1_W)+JetColl.at(Idxj2_W)).M();
    }

    if(JetColl.size()>=2){
      FillHist("Nj_SR2j", JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_SR2j", BJetColl.size(), weight, 0., 10., 10);
      FillHist("PTe_SR2j", EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("PTmu1_SR2j", MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_SR2j", MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("Mmumu_SR2j", Mmumu, weight, 0., 200., 40);
      FillHist("M2l2j_SR2j", M2l2j, weight, 0., 1000., 100);
      FillHist("M3lv_SR2j", M3lv, weight, 0., 1000., 100);
      if(Mmumu<40){
        FillHist("Mmumu_lt40_SR2j", Mmumu, weight, 0., 80., 80);
        FillHist("M2l2j_lt40_SR2j", M2l2j, weight, 0., 1000., 100);
        FillHist("M3lv_lt40_SR2j", M3lv, weight, 0., 1000., 100);
      }
      if(Mmumu<80){
        FillHist("Mmumu_lt80_SR2j", Mmumu, weight, 0., 80., 80);
        FillHist("M2l2j_lt80_SR2j", M2l2j, weight, 0., 1000., 100);
        FillHist("M3lv_lt80_SR2j", M3lv, weight, 0., 1000., 100);
      }
    }
    if(JetColl.size()>=3){
      FillHist("Nj_SR3j", JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_SR3j", BJetColl.size(), weight, 0., 10., 10);
      FillHist("PTe_SR3j", EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("PTmu1_SR3j", MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_SR3j", MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("Mmumu_SR3j", Mmumu, weight, 0., 200., 40);
      FillHist("M2l2j_SR3j", M2l2j, weight, 0., 1000., 100);
      FillHist("M3lv_SR3j", M3lv, weight, 0., 1000., 100);
      if(Mmumu<40){
        FillHist("Mmumu_lt40_SR3j", Mmumu, weight, 0., 80., 80);
        FillHist("M2l2j_lt40_SR3j", M2l2j, weight, 0., 1000., 100);
        FillHist("M3lv_lt40_SR3j", M3lv, weight, 0., 1000., 100);
      }
      if(Mmumu<80){
        FillHist("Mmumu_lt80_SR3j", Mmumu, weight, 0., 80., 80);
        FillHist("M2l2j_lt80_SR3j", M2l2j, weight, 0., 1000., 100);
        FillHist("M3lv_lt80_SR3j", M3lv, weight, 0., 1000., 100);
      }
    }

  }
  if(TriMu){
    if( !(MuLColl.size()==3 && EleLColl.size()==0) ) return;
    if( !(MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 && MuTColl.at(2).Pt()>10) ) return;

    int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    if( !(MOSSS1>12 && MOSSS2>12) ) return;
    if( fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10 ) return;
    if(BJetColl.size()==0) return;
    if(JetColl.size()>=1){
      FillHist("Nj_SR1j", JetColl.size(), weight, 0., 10., 10);
    }

    int Idxj1_W=-1, Idxj2_W=-1;
    snu::KParticle v; v.SetPxPyPzE(METx, METy, 0, sqrt(METx*METx+METy*METy));
    float Pzv1=GetvPz(v, MuTColl.at(0), 1), Pzv2=GetvPz(v, MuTColl.at(0), 2);
    float Pzv=fabs(Pzv1)<fabs(Pzv2)? Pzv1:Pzv2;
    v.SetXYZM(METx,METy,Pzv,0.);

    float M3lv=(MuTColl.at(0)+MuTColl.at(1)+MuTColl.at(2)+v).M(), M2l2j=0.; 
    if(JetColl.size()>=2){
      Idxj1_W=GetDaughterCandIdx(JetColl, "W", -1., "Lead");
      Idxj2_W=GetDaughterCandIdx(JetColl, "W", -1., "Subl");
      M2l2j=(MuTColl.at(1)+MuTColl.at(2)+JetColl.at(Idxj1_W)+JetColl.at(Idxj2_W)).M();
    }

    if(JetColl.size()>=2){
      FillHist("Nj_SR2j", JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_SR2j", BJetColl.size(), weight, 0., 10., 10);
      FillHist("PTmu1_SR2j", MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_SR2j", MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_SR2j", MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("MOSSS1_SR2j", MOSSS1, weight, 0., 200., 40);
      FillHist("MOSSS2_SR2j", MOSSS2, weight, 0., 200., 40);
      FillHist("M2l2j_SR2j", M2l2j, weight, 0., 1000., 100);
      FillHist("M3lv_SR2j", M3lv, weight, 0., 1000., 100);
      if(MOSSS2<40){
        FillHist("MOSSS2_lt40_SR2j", MOSSS2, weight, 0., 80., 80);
        FillHist("M2l2j_lt40_SR2j", M2l2j, weight, 0., 1000., 100);
        FillHist("M3lv_lt40_SR2j", M3lv, weight, 0., 1000., 100);
      }
      if(MOSSS2<80){
        FillHist("MOSSS2_lt80_SR2j", MOSSS2, weight, 0., 80., 80);
        FillHist("M2l2j_lt80_SR2j", M2l2j, weight, 0., 1000., 100);
        FillHist("M3lv_lt80_SR2j", M3lv, weight, 0., 1000., 100);
      }
    }
    if(JetColl.size()>=3){
      FillHist("Nj_SR3j", JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_SR3j", BJetColl.size(), weight, 0., 10., 10);
      FillHist("PTmu1_SR3j", MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_SR3j", MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_SR3j", MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("MOSSS1_SR3j", MOSSS1, weight, 0., 200., 40);
      FillHist("MOSSS2_SR3j", MOSSS2, weight, 0., 200., 40);
      FillHist("M2l2j_SR3j", M2l2j, weight, 0., 1000., 100);
      FillHist("M3lv_SR3j", M3lv, weight, 0., 1000., 100);
      if(MOSSS2<40){
        FillHist("MOSSS2_lt40_SR3j", MOSSS2, weight, 0., 80., 80);
        FillHist("M2l2j_lt40_SR3j", M2l2j, weight, 0., 1000., 100);
        FillHist("M3lv_lt40_SR3j", M3lv, weight, 0., 1000., 100);
      }
      if(MOSSS2<80){
        FillHist("MOSSS2_lt80_SR3j", MOSSS2, weight, 0., 80., 80);
        FillHist("M2l2j_lt80_SR3j", M2l2j, weight, 0., 1000., 100);
        FillHist("M3lv_lt80_SR3j", M3lv, weight, 0., 1000., 100);
      }
    }
  }


}


void Dec2017_TriLepSR::DoSystRun(TString Cycle, TString Mode, std::vector<snu::KElectron> EleColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight){

  int SystDir=0; TString SystKindLabel="", SystDirLabel="", ChannelLabel="";
  if(Mode.Contains("Syst")){
    if(isData&& !k_running_nonprompt) return;
    if     (Mode.Contains("Up"))      {SystDir= 1; SystDirLabel="_systup";}
    else if(Mode.Contains("Down"))    {SystDir=-1; SystDirLabel="_systdown";}
    if     (Mode.Contains("PU"))      SystKindLabel="_PU";
    else if(Mode.Contains("Nvtx"))    SystKindLabel="_Nvtx";
    else if(Mode.Contains("JES"))     SystKindLabel="_JES";
    else if(Mode.Contains("JER"))     SystKindLabel="_JER";
    else if(Mode.Contains("Uncl"))    SystKindLabel="_Uncl";
    else if(Mode.Contains("ElEn"))    SystKindLabel="_ElEn";
    else if(Mode.Contains("MuEn"))    SystKindLabel="_MuEn";
    else if(Mode.Contains("FR"))      SystKindLabel="_FR";
    else if(Mode.Contains("BTag_L"))  SystKindLabel="_BTag_L";
    else if(Mode.Contains("BTag_BC")) SystKindLabel="_BTag_BC";
  }
  if     (Mode.Contains("EMuMu"))   ChannelLabel="EMuMu";
  else if(Mode.Contains("TriMu"))   ChannelLabel="TriMu";


  if(Cycle=="SRYield"){
    CheckSRYield(MuColl, MuLColl, EleColl, EleLColl, JetColl,  BJetColl, MET, METx, METy, weight, SystDirLabel+SystKindLabel, ChannelLabel);
  }//End of Closure Cycle

  return;
}



float Dec2017_TriLepSR::ConeCorrectedPT(snu::KMuon Mu, float TightIsoCut){

  //float PTCorr=Mu.Pt();
  //float PTCorr=Mu.Pt()*(1.+RochIso(Mu,"0.4"));
  float PTCorr=Mu.Pt()*(1.+max((float) 0.,RochIso(Mu,"0.4")-TightIsoCut));
  return PTCorr;

}



float Dec2017_TriLepSR::ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut){

  //float PTCorr=Ele.Pt()*(1+Ele.PFRelIso(0.3));
  float PTCorr=Ele.Pt()*(1.+max(0.,Ele.PFRelIso(0.3)-TightIsoCut));
  return PTCorr;

}


float Dec2017_TriLepSR::FakeRateData(snu::KMuon Mu, TString Option){

  float FR=0., TightIsoCut=0.;
  int SystDir=0.; bool Syst_FR=Option.Contains("Syst");
  if(Syst_FR){ if(Option.Contains("Up")){ SystDir=1; } else if(Option.Contains("Down")){ SystDir=-1; } }
  if(Option.Contains("ConeSUSY")){
    if     (Option.Contains("TIsop15")) TightIsoCut = 0.15;
    else if(Option.Contains("TIsop20")) TightIsoCut = 0.2;
    else if(Option.Contains("TIsop25")) TightIsoCut = 0.25;
  }
  float PTCorr=ConeCorrectedPT(Mu, TightIsoCut);
    if(Option.Contains("FOPt")) PTCorr=Mu.Pt();
  float fEta=fabs(Mu.Eta());
  if(Option.Contains("POGTIsop20IPp01p1Chi3_POGLIsop6IPp5p1Chi30_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.336779 ;
      else if(PTCorr<20)  FR=0.111232 ;
      else if(PTCorr<25)  FR=0.0902   ;
      else if(PTCorr<35)  FR=0.0752774;
      else                FR=0.0670709;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.381494 ;
      else if(PTCorr<20)  FR=0.137689 ;
      else if(PTCorr<25)  FR=0.120664 ;
      else if(PTCorr<35)  FR=0.106082 ;
      else                FR=0.0923205;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.448765;
      else if(PTCorr<20)  FR=0.174516;
      else if(PTCorr<25)  FR=0.155327;
      else if(PTCorr<35)  FR=0.139882;
      else                FR=0.119441;
    }
    else{
      if     (PTCorr<15)  FR=0.464363;
      else if(PTCorr<20)  FR=0.202699;
      else if(PTCorr<25)  FR=0.17925 ;
      else if(PTCorr<35)  FR=0.167467;
      else                FR=0.138629;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p1Chi3_POGLIsop5IPp5p1Chi_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.336737;
      else if(PTCorr<20)  FR=0.138805;
      else if(PTCorr<25)  FR=0.119428;
      else if(PTCorr<35)  FR=0.102563;
      else                FR=0.095469;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.381256;
      else if(PTCorr<20)  FR=0.167567;
      else if(PTCorr<25)  FR=0.153358;
      else if(PTCorr<35)  FR=0.138543;
      else                FR=0.123659;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.448443;
      else if(PTCorr<20)  FR=0.205231;
      else if(PTCorr<25)  FR=0.191483;
      else if(PTCorr<35)  FR=0.176188;
      else                FR=0.153492;
    }
    else{
      if     (PTCorr<15)  FR=0.464362;
      else if(PTCorr<20)  FR=0.238831;
      else if(PTCorr<25)  FR=0.218972;
      else if(PTCorr<35)  FR=0.209289;
      else                FR=0.177849;
    }
  }
  else if(Option.Contains("POGTIsop15IPp01p1Chi3_POGLIsop6IPp5p1Chi30_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.260288 ;
      else if(PTCorr<20)  FR=0.0712013;
      else if(PTCorr<25)  FR=0.0583859;
      else if(PTCorr<35)  FR=0.0460673;
      else                FR=0.0396087;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.302753 ;
      else if(PTCorr<20)  FR=0.0919691;
      else if(PTCorr<25)  FR=0.0818808;
      else if(PTCorr<35)  FR=0.0667185;
      else                FR=0.0559974;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.362683 ;
      else if(PTCorr<20)  FR=0.122658 ;
      else if(PTCorr<25)  FR=0.104402 ;
      else if(PTCorr<35)  FR=0.0988028;
      else                FR=0.0768322;
    }
    else{
      if     (PTCorr<15)  FR=0.382679 ;
      else if(PTCorr<20)  FR=0.149735 ;
      else if(PTCorr<25)  FR=0.126996 ;
      else if(PTCorr<35)  FR=0.113983 ;
      else                FR=0.0974703;
    }
  }
  else if(Option.Contains("POGTIsop15IPp01p1Chi3_POGLIsop5IPp5p1Chi_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.260264 ;
      else if(PTCorr<20)  FR=0.0887165;
      else if(PTCorr<25)  FR=0.0772104;
      else if(PTCorr<35)  FR=0.0626961;
      else                FR=0.0563453;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.302525 ;
      else if(PTCorr<20)  FR=0.111957 ;
      else if(PTCorr<25)  FR=0.104062 ;
      else if(PTCorr<35)  FR=0.08743  ;
      else                FR=0.0750136;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.362481 ;
      else if(PTCorr<20)  FR=0.14414  ;
      else if(PTCorr<25)  FR=0.128693 ;
      else if(PTCorr<35)  FR=0.124427 ;
      else                FR=0.0987016;
    }
    else{
      if     (PTCorr<15)  FR=0.382678;
      else if(PTCorr<20)  FR=0.17642 ;
      else if(PTCorr<25)  FR=0.155347;
      else if(PTCorr<35)  FR=0.142424;
      else                FR=0.125045;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p1Chi3_POGLIsop6IPp5p1Chi30_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.224901 ;
      else if(PTCorr<20)  FR=0.107681 ;
      else if(PTCorr<25)  FR=0.0898675;
      else if(PTCorr<35)  FR=0.0804815;
      else                FR=0.0727239;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.267794 ;
      else if(PTCorr<20)  FR=0.136115 ;
      else if(PTCorr<25)  FR=0.131556 ;
      else if(PTCorr<35)  FR=0.106586 ;
      else                FR=0.0996845;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.319584;
      else if(PTCorr<20)  FR=0.186026;
      else if(PTCorr<25)  FR=0.163037;
      else if(PTCorr<35)  FR=0.148278;
      else                FR=0.129497;
    }
    else{
      if     (PTCorr<15)  FR=0.344918;
      else if(PTCorr<20)  FR=0.214086;
      else if(PTCorr<25)  FR=0.176514;
      else if(PTCorr<35)  FR=0.184554;
      else                FR=0.154967;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p1Chi3_POGLIsop5IPp5p1Chi_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.243216;
      else if(PTCorr<20)  FR=0.140185;
      else if(PTCorr<25)  FR=0.121813;
      else if(PTCorr<35)  FR=0.111958;
      else                FR=0.106649;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.285605;
      else if(PTCorr<20)  FR=0.171922;
      else if(PTCorr<25)  FR=0.168621;
      else if(PTCorr<35)  FR=0.141603;
      else                FR=0.137314;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.336817;
      else if(PTCorr<20)  FR=0.225045;
      else if(PTCorr<25)  FR=0.206099;
      else if(PTCorr<35)  FR=0.187497;
      else                FR=0.169763;
    }
    else{
      if     (PTCorr<15)  FR=0.363674;
      else if(PTCorr<20)  FR=0.257849;
      else if(PTCorr<25)  FR=0.217374;
      else if(PTCorr<35)  FR=0.234093;
      else                FR=0.201918;
    }
  }
  else if(Option.Contains("POGTIsop15IPp01p1Chi3_POGLIsop6IPp5p1Chi30_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.181152 ;
      else if(PTCorr<20)  FR=0.0707471;
      else if(PTCorr<25)  FR=0.0560844;
      else if(PTCorr<35)  FR=0.0483956;
      else                FR=0.0433495;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.219845 ;
      else if(PTCorr<20)  FR=0.0915236;
      else if(PTCorr<25)  FR=0.0883876;
      else if(PTCorr<35)  FR=0.0658354;
      else                FR=0.0599461;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.267231 ;
      else if(PTCorr<20)  FR=0.129102 ;
      else if(PTCorr<25)  FR=0.110895 ;
      else if(PTCorr<35)  FR=0.103245 ;
      else                FR=0.0819243;
    }
    else{
      if     (PTCorr<15)  FR=0.29299 ;
      else if(PTCorr<20)  FR=0.155378;
      else if(PTCorr<25)  FR=0.124455;
      else if(PTCorr<35)  FR=0.12888 ;
      else                FR=0.104704;
    }
  }
  else if(Option.Contains("POGTIsop15IPp01p1Chi3_POGLIsop5IPp5p1Chi_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.191885 ;
      else if(PTCorr<20)  FR=0.0915866;
      else if(PTCorr<25)  FR=0.0760472;
      else if(PTCorr<35)  FR=0.0671218;
      else                FR=0.0631913;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.23066  ;
      else if(PTCorr<20)  FR=0.115512 ;
      else if(PTCorr<25)  FR=0.113433 ;
      else if(PTCorr<35)  FR=0.0871235;
      else                FR=0.0821759;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.277334;
      else if(PTCorr<20)  FR=0.156263;
      else if(PTCorr<25)  FR=0.140098;
      else if(PTCorr<35)  FR=0.130675;
      else                FR=0.107183;
    }
    else{
      if     (PTCorr<15)  FR=0.305151;
      else if(PTCorr<20)  FR=0.187014;
      else if(PTCorr<25)  FR=0.152749;
      else if(PTCorr<35)  FR=0.164782;
      else                FR=0.136348;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi3_POGLIsop4IPp5p1_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.33842 ;
      else if(PTCorr<20)  FR=0.18468 ;
      else if(PTCorr<25)  FR=0.162466;
      else if(PTCorr<35)  FR=0.143624;
      else                FR=0.137429;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.372154;
      else if(PTCorr<20)  FR=0.210609;
      else if(PTCorr<25)  FR=0.193202;
      else if(PTCorr<35)  FR=0.181434;
      else                FR=0.163875;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.408703;
      else if(PTCorr<20)  FR=0.239243;
      else if(PTCorr<25)  FR=0.224734;
      else if(PTCorr<35)  FR=0.210624;
      else                FR=0.185866;
    }
    else{
      if     (PTCorr<15)  FR=0.397721;
      else if(PTCorr<20)  FR=0.250492;
      else if(PTCorr<25)  FR=0.241592;
      else if(PTCorr<35)  FR=0.228252;
      else                FR=0.20606 ;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi3_POGLIsop6IPp5p1_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.320068 ;
      else if(PTCorr<20)  FR=0.103925 ;
      else if(PTCorr<25)  FR=0.0828112;
      else if(PTCorr<35)  FR=0.0697201;
      else                FR=0.0599762;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.355258 ;
      else if(PTCorr<20)  FR=0.12828  ;
      else if(PTCorr<25)  FR=0.10927  ;
      else if(PTCorr<35)  FR=0.0966748;
      else                FR=0.0815638;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.394842;
      else if(PTCorr<20)  FR=0.155917;
      else if(PTCorr<25)  FR=0.137145;
      else if(PTCorr<35)  FR=0.120733;
      else                FR=0.101367;
    }
    else{
      if     (PTCorr<15)  FR=0.38534 ;
      else if(PTCorr<20)  FR=0.165267;
      else if(PTCorr<25)  FR=0.152769;
      else if(PTCorr<35)  FR=0.136345;
      else                FR=0.114394;
    }
  }
  else if(Option.Contains("HNTrilepTight2_HNTrilepFakeL2_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.314559 ;
      else if(PTCorr<20)  FR=0.115043 ;
      else if(PTCorr<25)  FR=0.0964687;
      else if(PTCorr<35)  FR=0.0742115;
      else                FR=0.0851192;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.35599  ;
      else if(PTCorr<20)  FR=0.142815 ;
      else if(PTCorr<25)  FR=0.129269 ;
      else if(PTCorr<35)  FR=0.0931603;
      else                FR=0.0976207;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.41679 ;
      else if(PTCorr<20)  FR=0.169059;
      else if(PTCorr<25)  FR=0.144527;
      else if(PTCorr<35)  FR=0.1434  ;
      else                FR=0.121598;
    }
    else{
      if     (PTCorr<15)  FR=0.428247;
      else if(PTCorr<20)  FR=0.215242;
      else if(PTCorr<25)  FR=0.188871;
      else if(PTCorr<35)  FR=0.161765;
      else                FR=0.158053;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.276234;
      else if(PTCorr<20)  FR=0.191673;
      else if(PTCorr<25)  FR=0.167833;
      else if(PTCorr<35)  FR=0.162589;
      else                FR=0.159873;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.307693;
      else if(PTCorr<20)  FR=0.220826;
      else if(PTCorr<25)  FR=0.216052;
      else if(PTCorr<35)  FR=0.18657 ;
      else                FR=0.18804 ;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.339593;
      else if(PTCorr<20)  FR=0.262275;
      else if(PTCorr<25)  FR=0.252044;
      else if(PTCorr<35)  FR=0.220044;
      else                FR=0.212143;
    }
    else{
      if     (PTCorr<15)  FR=0.339879;
      else if(PTCorr<20)  FR=0.271265;
      else if(PTCorr<25)  FR=0.249351;
      else if(PTCorr<35)  FR=0.257995;
      else                FR=0.23968 ;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop6IPp5p1_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.216822 ;
      else if(PTCorr<20)  FR=0.101333 ;
      else if(PTCorr<25)  FR=0.0833579;
      else if(PTCorr<35)  FR=0.0747546;
      else                FR=0.0662648;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.251757 ;
      else if(PTCorr<20)  FR=0.128096 ;
      else if(PTCorr<25)  FR=0.120662 ;
      else if(PTCorr<35)  FR=0.0952394;
      else                FR=0.0884953;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.28611 ;
      else if(PTCorr<20)  FR=0.166302;
      else if(PTCorr<25)  FR=0.146939;
      else if(PTCorr<35)  FR=0.126979;
      else                FR=0.110757;
    }
    else{
      if     (PTCorr<15)  FR=0.290039;
      else if(PTCorr<20)  FR=0.175547;
      else if(PTCorr<25)  FR=0.15415 ;
      else if(PTCorr<35)  FR=0.148194;
      else                FR=0.130693;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.344138;
      else if(PTCorr<20)  FR=0.187021;
      else if(PTCorr<25)  FR=0.16477 ;
      else if(PTCorr<35)  FR=0.145246;
      else                FR=0.139742;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.376491;
      else if(PTCorr<20)  FR=0.212285;
      else if(PTCorr<25)  FR=0.19565 ;
      else if(PTCorr<35)  FR=0.182581;
      else                FR=0.165813;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.413423;
      else if(PTCorr<20)  FR=0.242546;
      else if(PTCorr<25)  FR=0.227715;
      else if(PTCorr<35)  FR=0.214613;
      else                FR=0.187678;
    }
    else{
      if     (PTCorr<15)  FR=0.402425;
      else if(PTCorr<20)  FR=0.252608;
      else if(PTCorr<25)  FR=0.24362 ;
      else if(PTCorr<35)  FR=0.229344;
      else                FR=0.208459;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop6IPp5p1_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.325488 ;
      else if(PTCorr<20)  FR=0.10525  ;
      else if(PTCorr<25)  FR=0.0839633;
      else if(PTCorr<35)  FR=0.0704949;
      else                FR=0.0609975;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.35944  ;
      else if(PTCorr<20)  FR=0.129295 ;
      else if(PTCorr<25)  FR=0.110666 ;
      else if(PTCorr<35)  FR=0.0973101;
      else                FR=0.0825229;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.399404;
      else if(PTCorr<20)  FR=0.158063;
      else if(PTCorr<25)  FR=0.138922;
      else if(PTCorr<35)  FR=0.123003;
      else                FR=0.102357;
    }
    else{
      if     (PTCorr<15)  FR=0.389898;
      else if(PTCorr<20)  FR=0.166636;
      else if(PTCorr<25)  FR=0.153958;
      else if(PTCorr<35)  FR=0.137006;
      else                FR=0.115769;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.344294;
      else if(PTCorr<20)  FR=0.187095;
      else if(PTCorr<25)  FR=0.164859;
      else if(PTCorr<35)  FR=0.145231;
      else                FR=0.14    ;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.376607;
      else if(PTCorr<20)  FR=0.212338;
      else if(PTCorr<25)  FR=0.195648;
      else if(PTCorr<35)  FR=0.182624;
      else                FR=0.166085;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.413485;
      else if(PTCorr<20)  FR=0.242545;
      else if(PTCorr<25)  FR=0.227803;
      else if(PTCorr<35)  FR=0.214703;
      else                FR=0.187801;
    }
    else{
      if     (PTCorr<15)  FR=0.402582;
      else if(PTCorr<20)  FR=0.252604;
      else if(PTCorr<25)  FR=0.243616;
      else if(PTCorr<35)  FR=0.229269;
      else                FR=0.208675;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.276317;
      else if(PTCorr<20)  FR=0.191793;
      else if(PTCorr<25)  FR=0.167908;
      else if(PTCorr<35)  FR=0.162617;
      else                FR=0.160329;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.30778 ;
      else if(PTCorr<20)  FR=0.220861;
      else if(PTCorr<25)  FR=0.216049;
      else if(PTCorr<35)  FR=0.186631;
      else                FR=0.188504;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.339629;
      else if(PTCorr<20)  FR=0.262338;
      else if(PTCorr<25)  FR=0.252174;
      else if(PTCorr<35)  FR=0.220018;
      else                FR=0.212394;
    }
    else{
      if     (PTCorr<15)  FR=0.339974;
      else if(PTCorr<20)  FR=0.271258;
      else if(PTCorr<25)  FR=0.249348;
      else if(PTCorr<35)  FR=0.25788 ;
      else                FR=0.240026;
    }

    if(Syst_FR && SystDir<0){
      if(fEta<0.9){
        if     (PTCorr<15)   FR*=(1.- 0.0655089 );
        else if(PTCorr<20)   FR*=(1.- 0.304225  );
        else if(PTCorr<25)   FR*=(1.- 0.301855  );
        else if(PTCorr<35)   FR*=(1.- 0.387517  );
        else                 FR*=(1.- 0.483798  );
      }
      else if(fEta<1.6){
        if     (PTCorr<15)   FR*=(1.- 0.121366 );
        else if(PTCorr<20)   FR*=(1.- 0.365817 );
        else if(PTCorr<25)   FR*=(1.- 0.310581 );
        else if(PTCorr<35)   FR*=(1.- 0.27695  );
        else                 FR*=(1.- 0.438805 );
      }
      else if(fEta<2.1){
        if     (PTCorr<15)   FR*=(1.- 0.252924 );
        else if(PTCorr<20)   FR*=(1.- 0.368097 );
        else if(PTCorr<25)   FR*=(1.- 0.460861 );
        else if(PTCorr<35)   FR*=(1.- 0.242128 );
        else                 FR*=(1.- 0.415018 );
      }
      else{
        if     (PTCorr<15)   FR*=(1.- 0.243664 );
        else if(PTCorr<20)   FR*=(1.- 0.341107 );
        else if(PTCorr<25)   FR*=(1.- 0.442997 );
        else if(PTCorr<35)   FR*=(1.- 0.292832 );
        else                 FR*=(1.- 0.308315 );
      }
    }
    else if(Syst_FR && SystDir>0){
      if(fEta<0.9){
        if     (PTCorr<15)   FR*=(1.+ 0.0681634 );
        else if(PTCorr<20)   FR*=(1.+ 0.0553126 );
        else if(PTCorr<25)   FR*=(1.+ 0.101094  );
        else if(PTCorr<35)   FR*=(1.+ 0.07793   );
        else                 FR*=(1.+ 0.199967  );
      }
      else if(fEta<1.6){
        if     (PTCorr<15)   FR*=(1.+ 0.0625109);
        else if(PTCorr<20)   FR*=(1.+ 0.0856362);
        else if(PTCorr<25)   FR*=(1.+ 0.0966223);
        else if(PTCorr<35)   FR*=(1.+ 0.108313 );
        else                 FR*=(1.+ 0.146165 );
      }
      else if(fEta<2.1){
        if     (PTCorr<15)   FR*=(1.+ 0.0774913);
        else if(PTCorr<20)   FR*=(1.+ 0.0759516);
        else if(PTCorr<25)   FR*=(1.+ 0.147378 );
        else if(PTCorr<35)   FR*=(1.+ 0.131883 );
        else                 FR*=(1.+ 0.14489  );
      }
      else{
        if     (PTCorr<15)   FR*=(1.+ 0.0810315);
        else if(PTCorr<20)   FR*=(1.+ 0.0809193);
        else if(PTCorr<25)   FR*=(1.+ 0.0970186);
        else if(PTCorr<35)   FR*=(1.+ 0.110969 );
        else                 FR*=(1.+ 0.149472 );
      }
    }

  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop6IPp5p1_TrkIsoVVLFOPt")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.155258;
      else if(PTCorr<20)  FR=0.125076;
      else if(PTCorr<25)  FR=0.112006;
      else if(PTCorr<35)  FR=0.114503;
      else                FR=0.120845;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.189011;
      else if(PTCorr<20)  FR=0.158738;
      else if(PTCorr<25)  FR=0.159675;
      else if(PTCorr<35)  FR=0.143399;
      else                FR=0.155353;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.222263;
      else if(PTCorr<20)  FR=0.198951;
      else if(PTCorr<25)  FR=0.190684;
      else if(PTCorr<35)  FR=0.186306;
      else                FR=0.185893;
    }
    else{
      if     (PTCorr<15)  FR=0.227277;
      else if(PTCorr<20)  FR=0.205329;
      else if(PTCorr<25)  FR=0.203358;
      else if(PTCorr<35)  FR=0.223376;
      else                FR=0.222588;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLFOPt")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.247975;
      else if(PTCorr<20)  FR=0.217638;
      else if(PTCorr<25)  FR=0.202237;
      else if(PTCorr<35)  FR=0.195982;
      else                FR=0.214514;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.279385;
      else if(PTCorr<20)  FR=0.243698;
      else if(PTCorr<25)  FR=0.237864;
      else if(PTCorr<35)  FR=0.234583;
      else                FR=0.240188;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.309005;
      else if(PTCorr<20)  FR=0.284354;
      else if(PTCorr<25)  FR=0.274857;
      else if(PTCorr<35)  FR=0.272603;
      else                FR=0.281163;
    }
    else{
      if     (PTCorr<15)  FR=0.313636;
      else if(PTCorr<20)  FR=0.287708;
      else if(PTCorr<25)  FR=0.291454;
      else if(PTCorr<35)  FR=0.297233;
      else                FR=0.308716;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop6IPp5p1sig8_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.254794 ;
      else if(PTCorr<20)  FR=0.121792 ;
      else if(PTCorr<25)  FR=0.0954017;
      else if(PTCorr<35)  FR=0.0937622;
      else                FR=0.0795779;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.289057;
      else if(PTCorr<20)  FR=0.150987;
      else if(PTCorr<25)  FR=0.139373;
      else if(PTCorr<35)  FR=0.109927;
      else                FR=0.104292;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.322044;
      else if(PTCorr<20)  FR=0.18802 ;
      else if(PTCorr<25)  FR=0.163172;
      else if(PTCorr<35)  FR=0.147915;
      else                FR=0.130401;
    }
    else{
      if     (PTCorr<15)  FR=0.312404;
      else if(PTCorr<20)  FR=0.188555;
      else if(PTCorr<25)  FR=0.1747  ;
      else if(PTCorr<35)  FR=0.161571;
      else                FR=0.155695;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop6IPp5p1sig4_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.306808 ;
      else if(PTCorr<20)  FR=0.148199 ;
      else if(PTCorr<25)  FR=0.116904 ;
      else if(PTCorr<35)  FR=0.116722 ;
      else                FR=0.0994624;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.34711 ;
      else if(PTCorr<20)  FR=0.18568 ;
      else if(PTCorr<25)  FR=0.171709;
      else if(PTCorr<35)  FR=0.135723;
      else                FR=0.130152;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.382281;
      else if(PTCorr<20)  FR=0.226221;
      else if(PTCorr<25)  FR=0.195817;
      else if(PTCorr<35)  FR=0.178244;
      else                FR=0.160737;
    }
    else{
      if     (PTCorr<15)  FR=0.368437;
      else if(PTCorr<20)  FR=0.220804;
      else if(PTCorr<25)  FR=0.21116 ;
      else if(PTCorr<35)  FR=0.188548;
      else                FR=0.187302;
    }
  }

  return FR;
}



float Dec2017_TriLepSR::FakeRateData(snu::KElectron Ele, TString Option){

  float FR=0., TightIsoCut=0.;
  int SystDir=0.; bool Syst_FR=Option.Contains("Syst");
  if(Syst_FR){ if(Option.Contains("Up")){ SystDir=1; } else if(Option.Contains("Down")){ SystDir=-1; } }
  if(Option.Contains("ConeSUSY")){
    if     (Option.Contains("Isop06")) TightIsoCut = 0.06;
    else if(Option.Contains("Isop08")) TightIsoCut = 0.08;
    else if(Option.Contains("Isop1"))  TightIsoCut = 0.1;
  }
  float PTCorr=ConeCorrectedPT(Ele, TightIsoCut);
    if(Option.Contains("FOPt")) PTCorr=Ele.Pt();
  float fEta=fabs(Ele.Eta());

  if(Option.Contains("POGMVAMTIsop06IPp025p05Sig4FakeLIsop4")){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.119907 ;
      else if(PTCorr<50 ) FR= 0.065943 ;
      else                FR= 0.0682414;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.132106 ;
      else if(PTCorr<50 ) FR= 0.0682482;
      else                FR= 0.077325 ;
    }
    else{
      if     (PTCorr<35 ) FR= 0.160416 ;
      else if(PTCorr<50 ) FR= 0.0976898;
      else                FR= 0.109458 ;
    }
  }
  else if(Option.Contains("QCD_POGWP90Isop06IPp025p05sig4_LMVA06Isop4IPp025p05sig4_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.118837 ;
      else if(PTCorr<50)  FR=0.0626265;
      else                FR=0.0649617;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.129257 ;
      else if(PTCorr<50)  FR=0.0681096;
      else                FR=0.0732587;
    }
    else{
      if     (PTCorr<35)  FR=0.167549;
      else if(PTCorr<50)  FR=0.1026  ;
      else                FR=0.119577;
    }


    if(Syst_FR && SystDir<0){
      if(fEta<0.8){
        if     (PTCorr<35)   FR*=(1.- 0.077764 );
        else if(PTCorr<50)   FR*=(1.- 0.400603 );
        else                 FR*=(1.- 0.51524  );
      }
      else if(fEta<1.479){
        if     (PTCorr<35)   FR*=(1.- 0.169421 );
        else if(PTCorr<50)   FR*=(1.- 0.424976 );
        else                 FR*=(1.- 0.256287 );
      }
      else{
        if     (PTCorr<35)   FR*=(1.- 0.0561248 );
        else if(PTCorr<50)   FR*=(1.- 0.103     );
        else                 FR*=(1.- 0.111331  );
      }
    }
    else if(Syst_FR && SystDir>0){
      if(fEta<0.8){
        if     (PTCorr<35)   FR*=(1.+ 0.109443 );
        else if(PTCorr<50)   FR*=(1.+ 0.202818 );
        else                 FR*=(1.+ 0.374888 );
      }
      else if(fEta<1.479){
        if     (PTCorr<35)   FR*=(1.+ 0.0973585 );
        else if(PTCorr<50)   FR*=(1.+ 0.148029  );
        else                 FR*=(1.+ 0.246644  );
      }
      else{
        if     (PTCorr<35)   FR*=(1.+ 0.161426  );
        else if(PTCorr<50)   FR*=(1.+ 0.0707937 );
        else                 FR*=(1.+ 0.0566853 );
      }
    }

  }
  else if(Option.Contains("QCD_POGWP90Isop06IPp025p05sig4_LMVA06v1Isop4IPp025p05sig4_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.110603 ;
      else if(PTCorr<50)  FR=0.0644239;
      else                FR=0.0686986;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.136154 ;
      else if(PTCorr<50)  FR=0.0718016;
      else                FR=0.0789876;
    }
    else{
      if     (PTCorr<35)  FR=0.184978;
      else if(PTCorr<50)  FR=0.10697 ;
      else                FR=0.11768 ;
    }
  }
  else if(Option.Contains("QCD_POGWP90Isop06IPp025p05sig4_LMVA06v1Isop4IPp5p1_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.0736732;
      else if(PTCorr<50)  FR=0.0302377;
      else                FR=0.0339471;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.101205 ;
      else if(PTCorr<50)  FR=0.0460958;
      else                FR=0.0478649;
    }
    else{
      if     (PTCorr<35)  FR=0.11822  ;
      else if(PTCorr<50)  FR=0.0659004;
      else                FR=0.0723459;
    }
  }
  else if(Option.Contains("QCD_POGWP90Isop06IPp025p05sig4_LMVA06v2Isop4IPp025p05sig4_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.0974724;
      else if(PTCorr<50)  FR=0.0556398;
      else                FR=0.0571567;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.122288 ;
      else if(PTCorr<50)  FR=0.0634097;
      else                FR=0.0679921;
    }
    else{
      if     (PTCorr<35)  FR=0.166618 ;
      else if(PTCorr<50)  FR=0.0943797;
      else                FR=0.10212  ;
    }
  }
  else if(Option.Contains("QCD_POGWP90Isop08IPp025p05sig4_LMVA08v1Isop4IPp025p05sig4_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.154559 ;
      else if(PTCorr<50)  FR=0.0924185;
      else                FR=0.092448 ;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.193706;
      else if(PTCorr<50)  FR=0.107652;
      else                FR=0.11214 ;
    }
    else{
      if     (PTCorr<35)  FR=0.257529;
      else if(PTCorr<50)  FR=0.158485;
      else                FR=0.175704;
    }
  }
  else if(Option.Contains("QCD_POGWP90Isop08IPp025p05sig4_LMVA08v2Isop4IPp025p05sig4_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.137626 ;
      else if(PTCorr<50)  FR=0.0804705;
      else                FR=0.0782686;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.176499 ;
      else if(PTCorr<50)  FR=0.0959684;
      else                FR=0.0986659;
    }
    else{
      if     (PTCorr<35)  FR=0.233378;
      else if(PTCorr<50)  FR=0.14105 ;
      else                FR=0.154777;
    }
  }
  else if(Option.Contains("QCD_POGWP90Isop1IPp025p05sig4_LMVA1v1Isop4IPp025p05sig4_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.209062;
      else if(PTCorr<50)  FR=0.12362 ;
      else                FR=0.122594;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.254915;
      else if(PTCorr<50)  FR=0.145394;
      else                FR=0.149414;
    }
    else{
      if     (PTCorr<35)  FR=0.331213;
      else if(PTCorr<50)  FR=0.214573;
      else                FR=0.236407;
    }
  }
  else if(Option.Contains("QCD_POGWP90Isop1IPp025p05sig4_LMVA1v2Isop4IPp025p05sig4_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.185823;
      else if(PTCorr<50)  FR=0.107285;
      else                FR=0.103929;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.237566;
      else if(PTCorr<50)  FR=0.133319;
      else                FR=0.135416;
    }
    else{
      if     (PTCorr<35)  FR=0.306354;
      else if(PTCorr<50)  FR=0.19597 ;
      else                FR=0.213273;
    }
  }

  return FR;
} 


//float Dec2017_TriLepSR::OptimiseSelection(std::vector<snu::KMuon> MuPreColl, std::vector<snu::KElectron> ElePreColl, std::vector<snu::KJet> JetPreNoVetoColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, TString Label){
//
//  float Mmumu=
//  if(
//
//}

float Dec2017_TriLepSR::GetFakeWeight(std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleLColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, TString Option){

  float fakeweight=-1.; int NLooseNotTight=0; bool JSTrilepFR=false;
  //if(EleTID=="POGMVAMIP"      && EleLID=="HctoWAFakeLoose") EleFRLabel="POGMVAMTIsop06IPp025p05Sig4FakeLIsop4";
  if(MuTID =="HNTrilepTight2" && MuLID =="HNTrilepFakeL2" ) JSTrilepFR=true;

  TString FilterInfo="", ConeMethod="", SystOpt="";
  if     (Option.Contains("NoFilter"))  FilterInfo="NoFilter";
  else if(Option.Contains("TrkIsoVVL")) FilterInfo="TrkIsoVVL";
  if     (Option.Contains("ConeE"))     ConeMethod="ConeE";
  else if(Option.Contains("ConeSUSY"))  ConeMethod="ConeSUSY";
  else if(Option.Contains("FOPt"))      ConeMethod="FOPt";
  if     (Option.Contains("SystUp"))    SystOpt="SystUp";
  else if(Option.Contains("SystDown"))  SystOpt="SystDown";


  for(int i=0; i<MuLColl.size(); i++){
    if(!PassIDCriteria(MuLColl.at(i), MuTID)){
      float FR=0.;
      if(JSTrilepFR){
//        FR=m_datadriven_bkg->GetFakeObj()->getTrilepFakeRate_muon(false, MuLColl.at(i).Pt(), MuLColl.at(i).Eta());
      }
      else{
        FR=FakeRateData(MuLColl.at(i),MuTID.ReplaceAll("Test_","")+"_"+MuLID.ReplaceAll("Test_","")+"_"+FilterInfo+ConeMethod+SystOpt);
        //cout<<"MuFR"<<FR<<endl;
      }
      fakeweight*=-FR/(1-FR);
      NLooseNotTight++;
    }
  }
  for(int i=0; i<EleLColl.size(); i++){
    if(!PassIDCriteria(EleLColl.at(i), EleTID)){
      float FR=FakeRateData(EleLColl.at(i), "QCD_"+EleTID.ReplaceAll("Test_","")+"_"+EleLID.ReplaceAll("Test_","")+"_"+ConeMethod+SystOpt);
      fakeweight*=-FR/(1-FR);
      NLooseNotTight++;
        //cout<<"ElFR"<<FR<<endl;
    }
  }
  if(NLooseNotTight==0) return 0.;

  return fakeweight;

}


float Dec2017_TriLepSR::GetPreTrigPURW(int Nvtx){

   float weight=1.;
   if(Nvtx<2) weight=1.22087;
   else if(Nvtx<3) weight=4.79564;
   else if(Nvtx<4) weight=3.84689;
   else if(Nvtx<5) weight=2.47179;
   else if(Nvtx<6) weight=3.74121;
   else if(Nvtx<7) weight=3.16705;
   else if(Nvtx<8) weight=3.04533;
   else if(Nvtx<9) weight=2.85291;
   else if(Nvtx<10) weight=2.50907;
   else if(Nvtx<11) weight=2.29848;
   else if(Nvtx<12) weight=2.25329;
   else if(Nvtx<13) weight=1.9864;
   else if(Nvtx<14) weight=1.73597;
   else if(Nvtx<15) weight=1.57862;
   else if(Nvtx<16) weight=1.38849;
   else if(Nvtx<17) weight=1.21789;
   else if(Nvtx<18) weight=1.05361;
   else if(Nvtx<19) weight=0.909483;
   else if(Nvtx<20) weight=0.787018;
   else if(Nvtx<21) weight=0.699577;
   else if(Nvtx<22) weight=0.606739;
   else if(Nvtx<23) weight=0.533779;
   else if(Nvtx<24) weight=0.47826;
   else if(Nvtx<25) weight=0.430963;
   else if(Nvtx<26) weight=0.38291;
   else if(Nvtx<27) weight=0.368938;
   else if(Nvtx<28) weight=0.306386;
   else if(Nvtx<29) weight=0.321641;
   else if(Nvtx<30) weight=0.267767;
   else if(Nvtx<31) weight=0.303042;
   else if(Nvtx<32) weight=0.220255;
   else if(Nvtx<33) weight=0.237769;
   else if(Nvtx<34) weight=0.214341;
   else if(Nvtx<35) weight=0.254868;
   else if(Nvtx<36) weight=0.318919;
   else if(Nvtx<37) weight=0.290442;
   else if(Nvtx<38) weight=0.204838;
   else if(Nvtx<39) weight=0.389061;
   else if(Nvtx<40) weight=0.309147;
   else if(Nvtx<41) weight=0.286387;
   else if(Nvtx<42) weight=0.266532;
   else if(Nvtx<43) weight=0.327612;
   else if(Nvtx<44) weight=0.575947;
   else if(Nvtx<45) weight=0.23432;
   else if(Nvtx<46) weight=0.572198;
   else if(Nvtx<47) weight=0.541966;
   else if(Nvtx<48) weight=0.384433;
   else if(Nvtx<49) weight=1.57278;
   else weight=0.389116;

   return weight;
}



void Dec2017_TriLepSR::FillCutFlow(TString cut, float weight){
  
  if(GetHist("cutflow_W") && GetHist("cutflow_N")){
    GetHist("cutflow_W")->Fill(cut,weight);
    GetHist("cutflow_N")->Fill(cut,1);
  }
  else{
    if(!GetHist("cutflow_W")){
      AnalyzerCore::MakeHistograms("cutflow_W", 11, 0., 11.);
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(1,"NoCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(2,"TriggerCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(3,"EventCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(4,"VertexCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(5,"PreSel");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(6,"1e");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(7,"#geq1j");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(8,"1j");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(9,"met20");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(10,"ZVeto");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(11,"Nlj2");
    }
    if(!GetHist("cutflow_N")){
      AnalyzerCore::MakeHistograms("cutflow_N", 11, 0., 11.);
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(1,"NoCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(2,"TriggerCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(3,"EventCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(4,"VertexCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(5,"PreSel");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(6,"1e");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(7,"#geq1j");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(8,"1j");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(9,"met20");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(10,"ZVeto");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(11,"Nlj2");
    }
  }
}



void Dec2017_TriLepSR::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Dec2017_TriLepSR::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //Basic Hists/////////////////////////////////////////////////
  //Data properties before selection  
  AnalyzerCore::MakeHistograms("Basic_Nj_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Nb_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NmuL_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NeL_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NmuT_orig",10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NeT_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_METdist_orig", 100, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_NmuT_weT", 10, 0., 10.);

  AnalyzerCore::MakeHistograms("Basic_Ptmu1_3mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu2_3mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu1_1e2mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Pte_1e2mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j1_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j2_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j3_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j4_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_b1_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_b2_Et_orig", 200, 0., 200.);

  //After Nlcut
  AnalyzerCore::MakeHistograms("Basic_Nvtx_NoRW_wNlOScut", 50, 0., 50.);
  AnalyzerCore::MakeHistograms("Basic_Nvtx_PURW_wNlOScut", 50, 0., 50.);
  AnalyzerCore::MakeHistograms("Basic_MET_wNlOScut", 200, 0., 200.);
 


  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Dec2017_TriLepSRCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Dec2017_TriLepSR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
