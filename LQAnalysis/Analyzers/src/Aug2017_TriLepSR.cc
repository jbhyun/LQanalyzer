// $Id: Aug2017_TriLepSR.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQAug2017_TriLepSR Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Aug2017_TriLepSR.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Aug2017_TriLepSR);

 Aug2017_TriLepSR::Aug2017_TriLepSR() : AnalyzerCore(), out_muons(0) {

   SetLogName("Aug2017_TriLepSR");
   Message("In Aug2017_TriLepSR constructor", INFO);
   InitialiseAnalysis();
 }


 void Aug2017_TriLepSR::InitialiseAnalysis() throw( LQError ) {
   
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

void Aug2017_TriLepSR::ExecuteEvents()throw( LQError ){

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
      k_factor_weight          = GetKFactor();
      geneff_weight            = GenFilterEfficiency(k_sample_name);//For old samples
      gennorm_weight           = SignalNorm(k_sample_name, 20.);    //For old samples
      FillHist("Weight_PU", pileup_reweight, 1., 0., 20., 200);
   }
   weight *= k_factor_weight*geneff_weight*gennorm_weight*pileup_reweight;

 
   //Total Event(MCweight+PUrw)////////////////
   FillCutFlow("NoCut", weight*pileup_reweight);


   bool EMuMu=false, TriMu=false;
   bool CutOpt=false, ObsExpComp=false, SRYield=false, MAWinOpt=false;
   bool SRDist=false, SigBkgKin=false, CheckSystSize=false, CutFlowCheck=false, DecCompCheck=false;
   bool MACandOpt=false, MmumuShape=false;
   bool DoubleMuon=false, ElectronMuon=false;
   bool SystRun=false, TestRun=false;
   TString Cycle="";
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("EMuMu"))         EMuMu         = true;
     else if(k_flags.at(i).Contains("TriMu"))         TriMu         = true;
     else if(k_flags.at(i).Contains("SystRun"))       SystRun       = true;
     else if(k_flags.at(i).Contains("CutOpt"))        CutOpt        = true;
     else if(k_flags.at(i).Contains("MAWinOpt"))      MAWinOpt      = true;
     else if(k_flags.at(i).Contains("ObsExpComp"))   {ObsExpComp    = true; Cycle="ObsExpComp";}
     else if(k_flags.at(i).Contains("SRYield"))      {SRYield       = true; Cycle="SRYield";}
     else if(k_flags.at(i).Contains("SRDist"))       {SRDist        = true; Cycle="SRDist";}
     else if(k_flags.at(i).Contains("SigBkgKin"))    {SigBkgKin     = true; Cycle="SigBkgKin";}
     else if(k_flags.at(i).Contains("CutFlowCheck")) {CutFlowCheck  = true; }
     else if(k_flags.at(i).Contains("CheckSystSize")) CheckSystSize = true;
     else if(k_flags.at(i).Contains("DecCompCheck"))  DecCompCheck  = true;
     else if(k_flags.at(i).Contains("MACandOpt"))     MACandOpt     = true;
     else if(k_flags.at(i).Contains("MmumuShape"))    MmumuShape    = true;
     else if(k_flags.at(i).Contains("TestRun"))       TestRun       = true;
     else if(k_flags.at(i).Contains("DoubleMuon"))    DoubleMuon    = true;
     else if(k_flags.at(i).Contains("ElectronMuon"))  ElectronMuon  = true;
   }


    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //Trigger Path of Analysis
   bool Pass_Trigger=false, Pass_TriggerBG=false, Pass_TriggerH=false;
   float trigger_ps_weight=1., trigger_period_weight=1.;
   float LumiBG=27.257618, LumiH=8.605696, LumiBH=35.863314;

   if(EMuMu){
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") )    Pass_TriggerBG=true;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_TriggerH =true;

     if(isData){
       int DataPeriod=GetDataPeriod();
       if( DataPeriod>0 && DataPeriod<7 && Pass_TriggerBG ) Pass_Trigger=true;
       else if( DataPeriod==7 && Pass_TriggerH ) Pass_Trigger=true;
     }
     else{
       if( Pass_TriggerBG || Pass_TriggerH ) Pass_Trigger=true;
       trigger_period_weight=( (Pass_TriggerBG? 27.257618:0.)+(Pass_TriggerH? 8.605696:0.) )/35.863314;
       trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
     }
   }
   if(TriMu){
     if(  PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v")
        ||PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v") )    Pass_TriggerBG=true;
     if(  PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")
        ||PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ) Pass_TriggerH =true;

     if(isData){
       int DataPeriod=GetDataPeriod();
       if( DataPeriod>0 && DataPeriod<7 && Pass_TriggerBG ) Pass_Trigger=true;
       else if( DataPeriod==7 && Pass_TriggerH ) Pass_Trigger=true;
     }
     else{
       if( Pass_TriggerBG || Pass_TriggerH ) Pass_Trigger=true;
       trigger_period_weight=( (Pass_TriggerBG? 27.257618:0.)+(Pass_TriggerH? 8.605696:0.) )/35.863314;
       trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
     }
   }
   if(SigBkgKin || CheckSystSize || MmumuShape){ Pass_Trigger=true; trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi); }
   FillHist("Basic_TrigPSWeight", trigger_ps_weight, 1., 0., 1., 100);
   weight*=trigger_ps_weight*trigger_period_weight;

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
     eventbase->GetMuonSel()->SetPt(5.);                     eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, true);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
     if(EMuMu)      { if( !(electronPreColl.size()>=1 && muonPreColl.size()>=2) ) return; }
     else if(TriMu) { if( !(muonPreColl.size()>=3) ) return; }
     FillCutFlow("PreSel", weight);
   //******************************************************************************************************//

   //Primary Object Collection
   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);

     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.2);                 eventbase->GetMuonSel()->SetdxySigMax(4.);
     eventbase->GetMuonSel()->SetBSdz(0.1);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.6);
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(4.);
     eventbase->GetMuonSel()->SetBSdz(0.05);
     eventbase->GetMuonSel()->SetChiNdof(4.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.2);
   std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,true);
   std::vector<snu::KMuon> muonColl;  if(k_running_nonprompt){ muonColl=muonLooseColl;} else{ muonColl=muonTightColl;}


     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_HctoWA_FAKELOOSE);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.4, 0.4);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronTightColl; eventbase->GetElectronSel()->Selection(electronTightColl);
   std::vector<snu::KElectron> electronColl;  if(!k_running_nonprompt){ electronColl=electronTightColl;} else{ electronColl=electronLooseColl;}

   for(int i=0; i<electronColl.size(); i++){ FillHist("TestPTe", electronColl.at(i).Pt(), weight, 0., 200., 200); }


     bool LeptonVeto=false;
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
   std::vector<snu::KJet> jetNoVetoColl; eventbase->GetJetSel()->Selection(jetNoVetoColl, LeptonVeto, muonLooseColl, electronLooseColl);
   std::vector<snu::KJet> bjetNoVetoColl = SelBJets(jetNoVetoColl, "Medium");

   std::vector<snu::KJet> jetColl  = SkimJetColl(jetNoVetoColl,  electronLooseColl, muonLooseColl, "EleMuVeto");
   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");

   float met    = eventbase->GetEvent().MET(snu::KEvent::pfmet);
   float metphi = eventbase->GetEvent().METPhi();
   float met_x  = eventbase->GetEvent().PFMETx();
   float met_y  = eventbase->GetEvent().PFMETy();
   float Pzv, Pzv1, Pzv2;
   snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));
   //snu::KParticle v[4]; v[0].SetPx(met_x); v[0].SetPy(met_y);
   //                     v[1].SetPx(met_x); v[1].SetPy(met_y);

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

   if(EventCand & !SystRun && !CutFlowCheck){
     if(!isData){
       if(EMuMu || TriMu){
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl);
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
    
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl);
         trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);

         btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);
         if(TriMu){
           float trigger_sf1 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v");
           float trigger_sf2 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v");
           trigger_sf    = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
         }
         if(EMuMu){
           float trigger_sf1 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
           float trigger_sf2 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
           trigger_sf    = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
         }
       }
     }
     else{
       //Perfect Prompt Ratio Approximation applied.(p=1), Applicable to generic number, combination of leptons under premise of high prompt rate.
       if(k_running_nonprompt){
         fake_weight = GetFakeWeight_Data(muonLooseColl, electronLooseColl, "POGTIsop6IPp2p1sig4" , "POGTIsop20IPp01p05sig4Chi4", "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", "TrkIsoVVLConeSUSY");
       }
     }
   }
   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*btag_sf*trigger_sf*fake_weight;
   //-----------------------------------------------------------------------------------------//



   //----------------------------------------------------------------------------------//
   //==================================================================================//
   /////// MAIN ANALYSIS CODE ///////////////////////////////////////////////////////////
   //==================================================================================//
   //----------------------------------------------------------------------------------//

 
   if(SystRun){
     // Syst Sel. and Syst Corr. -------------------------------------------------------------------//
       eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
       eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
       eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
       eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
       eventbase->GetMuonSel()->SetdxySigMax(4.);
       eventbase->GetMuonSel()->SetChiNdof(4.);
     std::vector<snu::KMuon> MuTMuEnUpColl; eventbase->GetMuonSel()->Selection(MuTMuEnUpColl,false);
     //std::vector<snu::KMuon> MuTMuEnUpColl; eventbase->GetMuonSel()->Selection(MuTMuEnUpColl,true);
      // SetMuonResCorrection(MuTMuEnUpColl, "SystUp");

       eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
       eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
       eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
       eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
       eventbase->GetMuonSel()->SetdxySigMax(4.);
       eventbase->GetMuonSel()->SetChiNdof(4.);
     std::vector<snu::KMuon> MuTMuEnDownColl; eventbase->GetMuonSel()->Selection(MuTMuEnDownColl,false);
     //std::vector<snu::KMuon> MuTMuEnDownColl; eventbase->GetMuonSel()->Selection(MuTMuEnDownColl,true);
      // SetMuonResCorrection(MuTMuEnDownColl, "SystDown");


       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
       eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
       eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
       eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetdxySigMax(4.);
       eventbase->GetElectronSel()->SetApplyConvVeto(true);
     std::vector<snu::KElectron> EleTElEnUpColl; eventbase->GetElectronSel()->Selection(EleTElEnUpColl, "SystUpElEn");

       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
       eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
       eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
       eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
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
     float id_weight_ele_sfup=1.    , id_weight_ele_sfdown=1.    , id_weight_mu_sfup=1.    , id_weight_mu_sfdown=1.;
     float id_weight_ele_ElEnup=1.  , reco_weight_ele_ElEnup=1.  , id_weight_mu_MuEnup=1.  , trk_weight_mu_MuEnup=1.;
     float id_weight_ele_ElEndown=1., reco_weight_ele_ElEndown=1., id_weight_mu_MuEndown=1., trk_weight_mu_MuEndown=1.;
     float btag_sf_JESup=1.  , btag_sf_JERup=1.  , btag_sf_LTagup=1.  , btag_sf_BCTagup=1.;
     float btag_sf_JESdown=1., btag_sf_JERdown=1., btag_sf_LTagdown=1., btag_sf_BCTagdown=1.;
     float fake_weight_FRup=1., fake_weight_FRdown=1.;
     float ZG_SF=1., conv_up=1., conv_down=1.;
     float xsec_up=1., xsec_down=1.;
     float trigger_sf=1., trigger_sf_up=1., trigger_sf_down=1.;
     if(!isData){
       id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl);
       reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
       id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl);
       trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
       btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

       id_weight_ele_sfup       = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl,  1);
       id_weight_ele_sfdown     = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl, -1);
       id_weight_mu_sfup        = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl,  1);
       id_weight_mu_sfdown      = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl, -1);

       id_weight_ele_ElEnup     = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", EleTElEnUpColl);
       reco_weight_ele_ElEnup   = mcdata_correction->ElectronRecoScaleFactor(EleTElEnUpColl);
       id_weight_mu_MuEnup      = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuTMuEnUpColl);
       trk_weight_mu_MuEnup     = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnUpColl);

       id_weight_ele_ElEndown   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", EleTElEnDownColl);
       reco_weight_ele_ElEndown = mcdata_correction->ElectronRecoScaleFactor(EleTElEnDownColl);
       id_weight_mu_MuEndown    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuTMuEnDownColl);
       trk_weight_mu_MuEndown   = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnDownColl);

       if(TriMu){
         float trigger_sf1          = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v");
         float trigger_sf2          = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v");
         float trigger_sf1_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v"   ,"Leg1SystUp");
         float trigger_sf2_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","Leg1SystUp");
         float trigger_sf1_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v"   ,"Leg2SystUp");
         float trigger_sf2_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","Leg2SystUp");
         float trigger_sf1_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v"   ,"Leg1SystDown");
         float trigger_sf2_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","Leg1SystDown");
         float trigger_sf1_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v"   ,"Leg2SystDown");
         float trigger_sf2_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","Leg2SystDown");
         float trigger_sf2_dzup     = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","DZSystUp");
         float trigger_sf2_dzdown   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","DZSystDown");
         float trigger_sf1_up   = trigger_sf1+sqrt(pow(trigger_sf1_leg1up  -trigger_sf1,2)+pow(trigger_sf1_leg2up  -trigger_sf1,2));
         float trigger_sf1_down = trigger_sf1-sqrt(pow(trigger_sf1_leg1down-trigger_sf1,2)+pow(trigger_sf1_leg2down-trigger_sf1,2));
         float trigger_sf2_up   = trigger_sf2+sqrt(pow(trigger_sf2_leg1up  -trigger_sf2,2)+pow(trigger_sf2_leg2up  -trigger_sf2,2)+pow(trigger_sf2_dzup  -trigger_sf2,2));
         float trigger_sf2_down = trigger_sf2-sqrt(pow(trigger_sf2_leg1down-trigger_sf2,2)+pow(trigger_sf2_leg2down-trigger_sf2,2)+pow(trigger_sf2_dzdown-trigger_sf2,2));

         trigger_sf      = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
         trigger_sf_up   = ((Pass_TriggerBG ? trigger_sf1_up:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_up:0.)*LumiH)/LumiBH;
         trigger_sf_down = ((Pass_TriggerBG ? trigger_sf1_down:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_down:0.)*LumiH)/LumiBH;
       }
       else if(EMuMu){
         float trigger_sf1          = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
         float trigger_sf2          = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
         float trigger_sf1_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"   ,"Leg1SystUp");
         float trigger_sf2_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","Leg1SystUp");
         float trigger_sf1_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"   ,"Leg2SystUp");
         float trigger_sf2_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","Leg2SystUp");
         float trigger_sf1_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"   ,"Leg1SystDown");
         float trigger_sf2_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","Leg1SystDown");
         float trigger_sf1_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"   ,"Leg2SystDown");
         float trigger_sf2_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","Leg2SystDown");
         float trigger_sf2_dzup     = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","DZSystUp");
         float trigger_sf2_dzdown   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","DZSystDown");
         float trigger_sf1_up   = trigger_sf1+sqrt(pow(trigger_sf1_leg1up  -trigger_sf1,2)+pow(trigger_sf1_leg2up  -trigger_sf1,2));
         float trigger_sf1_down = trigger_sf1-sqrt(pow(trigger_sf1_leg1down-trigger_sf1,2)+pow(trigger_sf1_leg2down-trigger_sf1,2));
         float trigger_sf2_up   = trigger_sf2+sqrt(pow(trigger_sf2_leg1up  -trigger_sf2,2)+pow(trigger_sf2_leg2up  -trigger_sf2,2)+pow(trigger_sf2_dzup  -trigger_sf2,2));
         float trigger_sf2_down = trigger_sf2-sqrt(pow(trigger_sf2_leg1down-trigger_sf2,2)+pow(trigger_sf2_leg2down-trigger_sf2,2)+pow(trigger_sf2_dzdown-trigger_sf2,2));

         trigger_sf      = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
         trigger_sf_up   = ((Pass_TriggerBG ? trigger_sf1_up:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_up:0.)*LumiH)/LumiBH;
         trigger_sf_down = ((Pass_TriggerBG ? trigger_sf1_down:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_down:0.)*LumiH)/LumiBH;
       }

       btag_sf_LTagup = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpLTag");
       btag_sf_BCTagup= BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpBCTag");
       btag_sf_JESup  = BTagScaleFactor_1a(jetJESUpColl, snu::KJet::CSVv2, snu::KJet::Medium);
       btag_sf_JERup  = BTagScaleFactor_1a(jetJERUpColl, snu::KJet::CSVv2, snu::KJet::Medium);

       btag_sf_LTagdown = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownLTag");
       btag_sf_BCTagdown= BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownBCTag");
       btag_sf_JESdown  = BTagScaleFactor_1a(jetJESDownColl, snu::KJet::CSVv2, snu::KJet::Medium);
       btag_sf_JERdown  = BTagScaleFactor_1a(jetJERDownColl, snu::KJet::CSVv2, snu::KJet::Medium);

       xsec_up  +=GetXsecUncertainty(k_sample_name);
       xsec_down-=GetXsecUncertainty(k_sample_name);
       if(TriMu){
         if(k_sample_name.Contains("ZGto2LG")){ ZG_SF=0.86191; conv_up+=0.15; conv_down-=0.15; }
         else if(k_sample_name.Contains("TTG")){ conv_up+=0.5; conv_down-=0.5; }
       }
       if(EMuMu){
         if(k_sample_name.Contains("ZGto2LG")){ ZG_SF=0.963169; conv_up+=0.096; conv_down-=0.096; }
         else if(k_sample_name.Contains("TTG")){ conv_up+=0.5; conv_down-=0.5; }
       }
     }
     else if(k_running_nonprompt){
       fake_weight = GetFakeWeight_Data(muonLooseColl, electronLooseColl, "POGTIsop6IPp2p1sig4" , "POGTIsop20IPp01p05sig4Chi4", "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", "TrkIsoVVLConeSUSY");
     }
     weight*=ZG_SF;
       
     float systweight_central=weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight;
                                                                                                                                                        
     float systweight_Trigup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf_up*fake_weight;
     float systweight_ElIDup =weight*id_weight_ele_sfup  *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight;
     float systweight_MuIDup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu_sfup  *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight;
     float systweight_ElEnup =weight*id_weight_ele_ElEnup*reco_weight_ele_ElEnup*id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight;
     float systweight_MuEnup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu_MuEnup*trk_weight_mu_MuEnup*btag_sf        *trigger_sf   *fake_weight;
     float systweight_JESup  =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JESup  *trigger_sf   *fake_weight;
     float systweight_JERup  =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JERup  *trigger_sf   *fake_weight;
     float systweight_LTagup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_LTagup *trigger_sf   *fake_weight;
     float systweight_BCTagup=weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_BCTagup*trigger_sf   *fake_weight;
     float systweight_PUup   =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight*pileup_reweight_systup;
     float systweight_Xsecup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight*xsec_up;
     float systweight_Convup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight*conv_up;

     float systweight_Trigdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf_down*fake_weight;
     float systweight_ElIDdown =weight*id_weight_ele_sfdown  *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight;
     float systweight_MuIDdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu_sfdown  *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight;
     float systweight_ElEndown =weight*id_weight_ele_ElEndown*reco_weight_ele_ElEndown*id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight;
     float systweight_MuEndown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu_MuEndown*trk_weight_mu_MuEndown*btag_sf          *trigger_sf     *fake_weight;
     float systweight_JESdown  =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JESdown  *trigger_sf     *fake_weight;
     float systweight_JERdown  =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JERdown  *trigger_sf     *fake_weight;
     float systweight_LTagdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_LTagdown *trigger_sf     *fake_weight;
     float systweight_BCTagdown=weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_BCTagdown*trigger_sf     *fake_weight;
     float systweight_PUdown   =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight*pileup_reweight_systdown;
     float systweight_Xsecdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight*xsec_down;
     float systweight_Convdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight*conv_down;


     TString ChannelString = EMuMu? "EMuMu":"TriMu";
     //----------------------------------------------------------------------------------------------------------------------//


     DoSystRun(Cycle, ChannelString+"",
               electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
               systweight_central);
     if(!isData){
       DoSystRun(Cycle, ChannelString+"SystUpElID",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_ElIDup);
       DoSystRun(Cycle, ChannelString+"SystUpMuID",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_MuIDup);
       DoSystRun(Cycle, ChannelString+"SystUpTrig",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Trigup);
       DoSystRun(Cycle, ChannelString+"SystUpPU",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_PUup);
       DoSystRun(Cycle, ChannelString+"SystUpUncl",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi),
                 systweight_central);
       DoSystRun(Cycle, ChannelString+"SystUpJES",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESUpColl, bjetJESUpColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi),
                 systweight_JESup);
       DoSystRun(Cycle, ChannelString+"SystUpJER",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERUpColl, bjetJERUpColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi),
                 systweight_JERup);
       DoSystRun(Cycle, ChannelString+"SystUpBTag_L",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_LTagup);
       DoSystRun(Cycle, ChannelString+"SystUpBTag_BC",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_BCTagup);
       DoSystRun(Cycle, ChannelString+"SystUpElEn",
                 EleTElEnUpColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi),
                 systweight_ElEnup);
       DoSystRun(Cycle, ChannelString+"SystUpMuEn",
                 electronColl, electronLooseColl, MuTMuEnUpColl, muonLooseColl, jetColl, bjetColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi),
                 systweight_MuEnup);
       DoSystRun(Cycle, ChannelString+"SystUpXsec",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Xsecup);
       DoSystRun(Cycle, ChannelString+"SystUpConv",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Convup);



       DoSystRun(Cycle, ChannelString+"SystDownElID",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_ElIDdown);
       DoSystRun(Cycle, ChannelString+"SystDownMuID",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_MuIDdown);
       DoSystRun(Cycle, ChannelString+"SystDownTrig",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Trigdown);
       DoSystRun(Cycle, ChannelString+"SystDownPU",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_PUdown);
       DoSystRun(Cycle, ChannelString+"SystDownUncl",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi),
                 systweight_central);
       DoSystRun(Cycle, ChannelString+"SystDownJES",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESDownColl, bjetJESDownColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi),
                 systweight_JESdown);
       DoSystRun(Cycle, ChannelString+"SystDownJER",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERDownColl, bjetJERDownColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi),
                 systweight_JERdown);
       DoSystRun(Cycle, ChannelString+"SystDownBTag_L",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_LTagdown);
       DoSystRun(Cycle, ChannelString+"SystDownBTag_BC",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_BCTagdown);
       DoSystRun(Cycle, ChannelString+"SystDownElEn",
                 EleTElEnDownColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi),
                 systweight_ElEndown);
       DoSystRun(Cycle, ChannelString+"SystDownMuEn",
                 electronColl, electronLooseColl, MuTMuEnDownColl, muonLooseColl, jetColl, bjetColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi),
                 systweight_MuEndown);
       DoSystRun(Cycle, ChannelString+"SystDownXsec",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Xsecdown);
       DoSystRun(Cycle, ChannelString+"SystDownConv",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Convdown);
     }
   }
   if(CheckSystSize){
     
     for(int i=0; i<(int) muonColl.size(); i++){
       double dRelPt = mcdata_correction->GetRochesterMomentumWidth(muonColl.at(i));
       FillHist("Syst_MuRes", fabs(dRelPt), 1, 0., 1., 100);
     }
     for(int i=0; i<(int) electronColl.size(); i++){
       FillHist("Syst_ElRes", electronColl.at(i).PtShiftedUp()-1.  , 1, -1., 1., 200);
       FillHist("Syst_ElRes", electronColl.at(i).PtShiftedDown()-1., 1, -1., 1., 200);
     }
     for(int i=0; i<(int) jetColl.size(); i++){
       FillHist("Syst_JetES", jetColl.at(i).ScaledUpEnergy()-1.  , 1, -1., 1., 200);
       FillHist("Syst_JetES", jetColl.at(i).ScaledDownEnergy()-1., 1, -1., 1., 200);
     }
     for(int i=0; i<(int) jetColl.size(); i++){
       if(jetColl.at(i).IsMCSmeared()){
         FillHist("Syst_JetER", (jetColl.at(i).SmearedResUp()/jetColl.at(i).SmearedRes())-1.  , 1, -1., 1., 200);
         FillHist("Syst_JetER", (jetColl.at(i).SmearedResDown()/jetColl.at(i).SmearedRes())-1., 1, -1., 1., 200);
       }
       else{
         FillHist("Syst_JetER", jetColl.at(i).SmearedResUp()-1.  , 1, -1., 1., 200);
         FillHist("Syst_JetER", jetColl.at(i).SmearedResDown()-1., 1, -1., 1., 200);
       }
     }
     for(int i=0; i<(int) bjetColl.size(); i++){
       btag_sf          = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);
       float btag_sf_LTagup   = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpLTag")/btag_sf-1.;
       float btag_sf_BCTagup  = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpBCTag")/btag_sf-1.;
       float btag_sf_LTagdown = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownLTag")/btag_sf-1.;
       float btag_sf_BCTagdown= BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownBCTag")/btag_sf-1.;
       FillHist("Syst_BTag_L" , btag_sf_LTagup   , 1, -1., 1., 200);
       FillHist("Syst_BTag_L" , btag_sf_LTagdown , 1, -1., 1., 200);
       FillHist("Syst_BTag_BC", btag_sf_BCTagup  , 1, -1., 1., 200);
       FillHist("Syst_BTag_BC", btag_sf_BCTagdown, 1, -1., 1., 200);
     }

   }
   if(SRDist){
     if(EMuMu && !SystRun){
       CheckSRDist(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMu");
     }
     if(TriMu && !SystRun){
       CheckSRDist(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMu");
     }
   }
   if(SigBkgKin){
     CheckSigBkgKinematics(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "");
   }
   if(DecCompCheck){
     CheckDecComp(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMu");
   }
   if(MACandOpt){
     if(TriMu && (!SystRun)){
        OptimizeTriMuAssign(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMu");
     }
   }
   if(CutFlowCheck){
     std::vector<snu::KMuon> muonTight5Coll, muonLoose5Coll;
       for(int i=0; i<(int) muonPreColl.size(); i++){
         if(PassIDCriteria(muonPreColl.at(i),"POGTIsop20IPp01p05sig4Chi4","Roch")) muonTight5Coll.push_back(muonPreColl.at(i));
         if(PassIDCriteria(muonPreColl.at(i),       "POGTIsop6IPp2p1sig4","Roch")) muonLoose5Coll.push_back(muonPreColl.at(i));
       }

     if(EMuMu && (!SystRun)){
       CheckCutflow(muonTight5Coll, muonLoose5Coll, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMu");
     }
     if(TriMu && (!SystRun)){
       CheckCutflow(muonTight5Coll, muonLoose5Coll, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMu");

       bool ConvenerQuestion=false;

       if(ConvenerQuestion){
           eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
           eventbase->GetJetSel()->SetPt(10.);                     eventbase->GetJetSel()->SetEta(2.4);
         std::vector<snu::KJet> jetNoVeto10Coll; eventbase->GetJetSel()->Selection(jetNoVeto10Coll, false, muonLooseColl, electronLooseColl);
         std::vector<snu::KJet> bjetNoVeto10Coll = SelBJets(jetNoVeto10Coll, "Medium");
  
         CheckCutflow(muonColl, muonLooseColl, electronColl, electronLooseColl, jetNoVetoColl  , bjetNoVetoColl  , met, met*cos(metphi), met*sin(metphi), weight, "_NoVetoJet"  , "TriMu");
         CheckCutflow(muonColl, muonLooseColl, electronColl, electronLooseColl, jetNoVeto10Coll, bjetNoVeto10Coll, met, met*cos(metphi), met*sin(metphi), weight, "_NoVetoJet10", "TriMu");
       }
     }
   }
   if(SRYield){
     if(EMuMu && (!SystRun)){
        CheckSRYield(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMu");
     }
   }
   if(CutOpt){

     if(MAWinOpt){
         eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
         eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
         eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
         eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
         eventbase->GetMuonSel()->SetdxySigMax(4.);
         eventbase->GetMuonSel()->SetChiNdof(4.);
       std::vector<snu::KMuon> MuTMuEnUpColl; eventbase->GetMuonSel()->Selection(MuTMuEnUpColl, true);
       SetMuonResCorrection(MuTMuEnUpColl, "SystUp");

         eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
         eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
         eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
         eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
         eventbase->GetMuonSel()->SetdxySigMax(4.);
         eventbase->GetMuonSel()->SetChiNdof(4.);
       std::vector<snu::KMuon> MuTMuEnDownColl; eventbase->GetMuonSel()->Selection(MuTMuEnDownColl,true);
       SetMuonResCorrection(MuTMuEnDownColl, "SystDown");


       if(EMuMu){
         OptimizeMACut(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight,
                       "", "EMuMu");
         if(!isData){
           OptimizeMACut(MuTMuEnUpColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight,
                         "_systupMuEn", "EMuMu");
           OptimizeMACut(MuTMuEnDownColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight,
                         "_systdownMuEn", "EMuMu");
         }
       }
       if(TriMu){
         OptimizeMACut(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight,
                       "", "TriMu");
         if(!isData){
           OptimizeMACut(MuTMuEnUpColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight,
                         "_systupMuEn", "TriMu");
           OptimizeMACut(MuTMuEnDownColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight,
                         "_systdownMuEn", "TriMu");
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
         if((int) bjetColl.size()<BJetCuts[it_nb]) continue;
         FillHist("Count_NBCut", it_nb, weight, 0., 2., 2);
       }
       for(int it_nj=0; it_nj<NJetCut; it_nj++){
         if((int) jetColl.size()<JetCuts[it_nj]) continue;
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
         if((int)jetColl.size()>1 && jetColl.at(1).Pt()<JetPtCuts[it_ptj]) continue;
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
         if( (int) jetColl.size()<JetCuts[it_nj] || (int) bjetColl.size()<BJetCuts[it_nb] ) continue;
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
         if((int) bjetColl.size()<BJetCuts[it_nb]) continue;
         FillHist("Count_NBCut", it_nb, weight, 0., 2., 2);
       }
       for(int it_nj=0; it_nj<NJetCut; it_nj++){
         if((int) jetColl.size()<JetCuts[it_nj]) continue;
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
         if((int) jetColl.size()>1 && jetColl.at(1).Pt()<JetPtCuts[it_ptj]) continue;
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
         if( (int) jetColl.size()<JetCuts[it_nj] || (int) bjetColl.size()<BJetCuts[it_nb] ) continue;
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
   if(MmumuShape){
     TString OptionStr="";
     if(EMuMu) OptionStr="EMuMu";
     else if(TriMu) OptionStr="TriMu";
     AnalyzeMmumuShape(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", OptionStr);
   }
/////////////////////////////////////////////////////////////////////////////////// 

return;

}// End of execute event loop
  


void Aug2017_TriLepSR::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Aug2017_TriLepSR::BeginCycle() throw( LQError ){
  
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

Aug2017_TriLepSR::~Aug2017_TriLepSR() {
  
  Message("In Aug2017_TriLepSR Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}


void Aug2017_TriLepSR::AnalyzeMmumuShape(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){

  bool EMuMu=false, TriMu=false;
  //bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");
  if( EleLColl.size()==1 && MuLColl.size()==2 ) EMuMu=true;
  if( EleLColl.size()==0 && MuLColl.size()==3 ) TriMu=true;

  const int NSignal = 7;
  float Mmumu_CentSeed[NSignal]  = {15., 25., 35., 45., 55., 65., 75.};
  float Mmumu_WidthSeed[NSignal] = {1., 1., 1.3, 1.5, 1.5, 2.0, 2.5};
  float Mmumu_Step[NSignal-1]    = {1., 1., 1.2, 1.2, 1.2, 1.8};
  vector<float> Mmumu_CentFull, Mmumu_WidthFull;
  for(int it_seed=0; it_seed<NSignal; it_seed++){
    float CurrentCenter = Mmumu_CentSeed[it_seed];
    float CurrentWidth  = Mmumu_WidthSeed[it_seed];
    Mmumu_CentFull.push_back(CurrentCenter);
    Mmumu_WidthFull.push_back(CurrentWidth);
    for(int it_step=0; it_step<100; it_step++){
      if(it_seed==NSignal-1) break;
      CurrentCenter+= Mmumu_Step[it_seed];
      CurrentWidth  = (CurrentCenter-Mmumu_CentSeed[it_seed])/10.*Mmumu_WidthSeed[it_seed+1]
                     +(Mmumu_CentSeed[it_seed+1]-CurrentCenter)/10.*Mmumu_WidthSeed[it_seed];
      if(CurrentCenter<Mmumu_CentSeed[it_seed]+10.){
        Mmumu_CentFull.push_back(CurrentCenter);
        Mmumu_WidthFull.push_back(CurrentWidth);
      }
      else break;
    }
  }


  if(EMuMu){
    if( !(EleLColl.size()==1 && MuLColl.size()==2) ) return;
    if( !(EleTColl.size()==1 && MuTColl.size()==2) ) return;
    if( !(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) ) return;
    if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && MuTColl.at(1).Pt()>10) ) return;
    if( fabs(EleTColl.at(0).Eta())>2.5 ) return;
  
    float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
    if( Mmumu<12 ) return;
    //if( Mmumu<12 || fabs(Mmumu-91.2)<10 ) return;
  
    if(BJetColl.size()==0) return;
    if(JetColl.size()<2) return;

    float MA=GetHiggsMass(k_sample_name, "A");
    int IdxSeed=-1;
    for(int i=0; i<NSignal; i++){ if(fabs(MA-Mmumu_CentSeed[i])<0.01){ IdxSeed=i; break; } }
  
    FillHist("NCount_SR_1e2mu"+Label, 0., weight, 0., 1., 1);

    if(!isData){
      std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
      std::vector<snu::KMuon> MuAColl;
      int MuA=0, MuW=0, MutaW=0;
        for(int i=0; i<(int) MuTColl.size(); i++){
          int LepType=GetLeptonType(MuTColl.at(i),truthColl);
          if(LepType==2){ MuAColl.push_back(MuTColl.at(i)); MuA++; }
          else if(LepType==1) MuW++;
          else if(LepType==3) MutaW++;
        }
      float IdxFillBin=0.; 
      if     (MuA==2            ) IdxFillBin=1.;
      else if(MuA==1 && MuW==1  ) IdxFillBin=2.;
      else if(MuA==1 && MutaW==1) IdxFillBin=3.;
      else if(MuW==2            ) IdxFillBin=4.;
      else if(MuW==1 && MutaW==1) IdxFillBin=5.;
      FillHist("SRmumu_GenComp_1e2mu"+Label, 0., weight, 0., 6., 6);
      FillHist("SRmumu_GenComp_1e2mu"+Label, IdxFillBin+0.01, weight, 0., 6., 6);
      if(IdxSeed==-1) FillHist("HistError", 0., weight, 0., 1., 1);
      else{
        if(fabs(Mmumu-Mmumu_CentSeed[IdxSeed])<Mmumu_WidthSeed[IdxSeed]){
          FillHist("SRWinmumu_GenComp_1e2mu"+Label, 0., weight, 0., 6., 6);
          FillHist("SRWinmumu_GenComp_1e2mu"+Label, IdxFillBin+0.01, weight, 0., 6., 6);
        }
      }

      if(MuAColl.size()>1){
        FillHist("NCount_muAAX_1e2mu"+Label, 0., weight, 0., 1., 1);
        float Mmumu_A=(MuAColl.at(0)+MuAColl.at(1)).M();
        FillHist("Mmumu_A"+Label, Mmumu_A, weight, 10., 110., 2000);
      }
    }
  }
  if(TriMu){
    if( !(EleLColl.size()==0 && MuLColl.size()==3) ) return;
    if( !(EleTColl.size()==0 && MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 && MuTColl.at(2).Pt()>10) ) return;
  
    int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    int   IdxSSW = TriMuChargeIndex(MuTColl, MET, METx, METy, "SSW");
    int   IdxSSA = TriMuChargeIndex(MuTColl, MET, METx, METy, "SSA");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    float MmumuWA= (MuTColl.at(IdxOS)+MuTColl.at(IdxSSW)).M();
    float Mmumu  = (MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)).M();
    if( !(MOSSS1>12 && MOSSS2>12) ) return;
    //if( fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10 ) return;

    if(BJetColl.size()==0) return;
    if(JetColl.size()<2) return;
                                 
    float MA=GetHiggsMass(k_sample_name, "A");
    int   IdxSeed=-1;
    for(int i=0; i<NSignal; i++){ if(fabs(MA-Mmumu_CentSeed[i])<0.01){ IdxSeed=i; break; } }
  
    FillHist("NCount_SR_3mu"+Label, 0., weight, 0., 1., 1);

    if(!isData){
      std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
      std::vector<snu::KMuon> MuAColl;
      int MuA=0, MuW=0, MutaW=0;
      float PairType=0.;
      bool OSisA=false, OSisW=false, OSista=false, SSisA=false, SSisW=false, SSista=false;
        for(int i=0; i<(int) MuTColl.size(); i++){
          int LepType=GetLeptonType(MuTColl.at(i),truthColl);
          if(LepType==2){ MuAColl.push_back(MuTColl.at(i)); MuA++; }
          else if(LepType==1) MuW++;
          else if(LepType==3) MutaW++;

          if     (i==IdxOS  && LepType==1) OSisW =true;
          else if(i==IdxOS  && LepType==2) OSisA =true;
          else if(i==IdxOS  && LepType==3) OSista=true;
          else if(i==IdxSSA && LepType==1) SSisW =true;
          else if(i==IdxSSA && LepType==2) SSisA =true;
          else if(i==IdxSSA && LepType==3) SSista=true;
        }
      if     (  OSisA  && SSisA                        ) PairType = 1.;
      else if( (OSisA  && SSisW)  || (OSisW && SSisA)  ) PairType = 2.;
      else if( (OSisA  && SSista) || (OSisW && SSista) ) PairType = 3.;
      else if(  OSisW  && SSisW                        ) PairType = 4.;
      else if( (OSisW  && SSista) || (OSista && SSisW) ) PairType = 5.;
      else if(  OSista && SSista                       ) PairType = 6.;
      FillHist("TriMu_PairType"+Label, PairType+0.01, weight, 0., 7., 7);

      
      float IdxFillBin=0.; 
      if     (MuA==2            ) IdxFillBin=1.;
      else if(MuA==1 && MuW==1  ) IdxFillBin=2.;
      else if(MuA==1 && MutaW==1) IdxFillBin=3.;
      else if(MuW==2            ) IdxFillBin=4.;
      else if(MuW==1 && MutaW==1) IdxFillBin=5.;
      FillHist("SRmumu_GenComp_3mu"+Label, 0., weight, 0., 6., 6);
      FillHist("SRmumu_GenComp_3mu"+Label, IdxFillBin+0.01, weight, 0., 6., 6);
      if(IdxSeed==-1) FillHist("HistError", 0., weight, 0., 1., 1);
      else{
        if(fabs(Mmumu-Mmumu_CentSeed[IdxSeed])<Mmumu_WidthSeed[IdxSeed]){
          FillHist("SRWinmumu_GenComp_3mu"+Label, 0., weight, 0., 6., 6);
          FillHist("SRWinmumu_GenComp_3mu"+Label, IdxFillBin+0.01, weight, 0., 6., 6);
        }
      }

      if(MuAColl.size()>1){
        FillHist("NCount_muAAX_3mu"+Label, 0., weight, 0., 1., 1);
        float Mmumu_A=(MuAColl.at(0)+MuAColl.at(1)).M();
        FillHist("Mmumu_A"+Label, Mmumu_A, weight, 10., 110., 2000);
      }
    }
  }
}


void Aug2017_TriLepSR::OptimizeMACut(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){

  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");
  if(EMuMu){
    if( !(EleLColl.size()==1 && MuLColl.size()==2) ) return;
    if( !(EleTColl.size()==1 && MuTColl.size()==2) ) return;
    if( !(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) ) return;
    if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && MuTColl.at(1).Pt()>10) ) return;
    if( fabs(EleTColl.at(0).Eta())>2.5 ) return;
  
    float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
    if( Mmumu<12 || fabs(Mmumu-91.2)<10 ) return;
  
    if(BJetColl.size()==0) return;
    if(JetColl.size()<2) return;
  
    const int Nwindow=100;
    float MwidthCuts[Nwindow]={10.0, 9.9, 9.8, 9.7, 9.6, 9.5, 9.4, 9.3, 9.2, 9.1,
                                9.0, 8.9, 8.8, 8.7, 8.6, 8.5, 8.4, 8.3, 8.2, 8.1,
                                8.0, 7.9, 7.8, 7.7, 7.6, 7.5, 7.4, 7.3, 7.2, 7.1,
                                7.0, 6.9, 6.8, 6.7, 6.6, 6.5, 6.4, 6.3, 6.2, 6.1,
                                6.0, 5.9, 5.8, 5.7, 5.6, 5.5, 5.4, 5.3, 5.2, 5.1,
                                5.0, 4.9, 4.8, 4.7, 4.6, 4.5, 4.4, 4.3, 4.2, 4.1,
                                4.0, 3.9, 3.8, 3.7, 3.6, 3.5, 3.4, 3.3, 3.2, 3.1,
                                3.0, 2.9, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 
                                2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1,
                                1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};
    

    FillHist("NCount_Win"+Label, 0., weight, 0., 1., 1);
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-15.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM15"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-25.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM25"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-35.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM35"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-45.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM45"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-55.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM55"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-65.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM65"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-75.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM75"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }



    if(!isData){
      std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
      std::vector<snu::KMuon> MuAColl;
        for(int i=0; i<(int) MuTColl.size(); i++){
          int LepType=GetLeptonType(MuTColl.at(i),truthColl);
          if(LepType==2) MuAColl.push_back(MuTColl.at(i));
        }

      if(MuAColl.size()>1){
        FillHist("NCount_TrueWin"+Label, 0., weight, 0., 1., 1);
        float Mmumu_A=(MuAColl.at(0)+MuAColl.at(1)).M();
        FillHist("Mmumu_A"+Label, Mmumu_A, weight, 10., 90., 8000);
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-15.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM15"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-25.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM25"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-35.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM35"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-45.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM45"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-55.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM55"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-65.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM65"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-75.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM75"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
      }
    }
  }
  if(TriMu){
    if( !(EleLColl.size()==0 && MuLColl.size()==3) ) return;
    if( !(EleTColl.size()==0 && MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 && MuTColl.at(2).Pt()>10) ) return;
  
    int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    int   IdxSSA = TriMuChargeIndex(MuTColl, MET, METx, METy, "SSA");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    float Mmumu  = (MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)).M();
    if( !(MOSSS1>12 && MOSSS2>12) ) return;
    if( fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10 ) return;

    if(BJetColl.size()==0) return;
    if(JetColl.size()<2) return;
  
    const int Nwindow=100;
    float MwidthCuts[Nwindow]={10.0, 9.9, 9.8, 9.7, 9.6, 9.5, 9.4, 9.3, 9.2, 9.1,
                                9.0, 8.9, 8.8, 8.7, 8.6, 8.5, 8.4, 8.3, 8.2, 8.1,
                                8.0, 7.9, 7.8, 7.7, 7.6, 7.5, 7.4, 7.3, 7.2, 7.1,
                                7.0, 6.9, 6.8, 6.7, 6.6, 6.5, 6.4, 6.3, 6.2, 6.1,
                                6.0, 5.9, 5.8, 5.7, 5.6, 5.5, 5.4, 5.3, 5.2, 5.1,
                                5.0, 4.9, 4.8, 4.7, 4.6, 4.5, 4.4, 4.3, 4.2, 4.1,
                                4.0, 3.9, 3.8, 3.7, 3.6, 3.5, 3.4, 3.3, 3.2, 3.1,
                                3.0, 2.9, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 
                                2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1,
                                1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};
                               
    
    FillHist("NCount_Win"+Label, 0., weight, 0., 1., 1);
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-15.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM15"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-25.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM25"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-35.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM35"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-45.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM45"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-55.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM55"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-65.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM65"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-75.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM75"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }



    if(!isData){
      std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
      std::vector<snu::KMuon> MuAColl;
        bool OSisA=false, SSAisA=false;
        for(int i=0; i<(int) MuTColl.size(); i++){
          int LepType=GetLeptonType(MuTColl.at(i),truthColl);
          if(LepType==2) MuAColl.push_back(MuTColl.at(i));
          if(i==IdxOS && LepType==2) OSisA=true;
          if(i==IdxSSA && LepType==2) SSAisA=true;
        }

      if(MuAColl.size()>1){
        FillHist("NCount_TrueWin"+Label, 0., weight, 0., 1., 1);
        float Mmumu_A=(MuAColl.at(0)+MuAColl.at(1)).M();
        FillHist("Mmumu_A"+Label, Mmumu_A, weight, 10., 90., 8000);
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-15.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM15"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
            if(OSisA && SSAisA) FillHist("NCount_TrueWinM15_Cor"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-25.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM25"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
            if(OSisA && SSAisA) FillHist("NCount_TrueWinM25_Cor"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-35.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM35"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
            if(OSisA && SSAisA) FillHist("NCount_TrueWinM35_Cor"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-45.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM45"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
            if(OSisA && SSAisA) FillHist("NCount_TrueWinM45_Cor"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-55.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM55"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
            if(OSisA && SSAisA) FillHist("NCount_TrueWinM55_Cor"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-65.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM65"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
            if(OSisA && SSAisA) FillHist("NCount_TrueWinM65_Cor"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-75.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM75"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
            if(OSisA && SSAisA) FillHist("NCount_TrueWinM75_Cor"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
      }
    }
  }
}


void Aug2017_TriLepSR::OptimizeTriMuAssign(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){

  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");
  if(TriMu){
    if( !(MuLColl.size()==3 && EleLColl.size()==0) ) return;
    if( !(MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 && MuTColl.at(2).Pt()>10) ) return;
      float CurrentIdx=0.;
      int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
      int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
      int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
      float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
      float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    if( !(MOSSS1>12 && MOSSS2>12) ) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//3lep within simulation reg. & Triggered.
    if( fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10 ) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//3lep within simulation reg. & Triggered.
    if(BJetColl.size()==0) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//HasB
    if(JetColl.size()<2) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//geq2j

    //Now optimize Assign
    //1st check MT Dist.
    std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
    int LepTypeSS1=GetLeptonType(MuTColl.at(IdxSS1), truthColl);
    int LepTypeSS2=GetLeptonType(MuTColl.at(IdxSS2), truthColl);
    if( !( (LepTypeSS1==1 && LepTypeSS2==2) || (LepTypeSS1==2 && LepTypeSS2==1) ) ) return;
    int IdxSSA=LepTypeSS1==2? IdxSS1:IdxSS2;
    int IdxSSW=LepTypeSS1==1? IdxSS1:IdxSS2;
    float METphi  = eventbase->GetEvent().METPhi();
      snu::KParticle v; v.SetPxPyPzE(METx, METy, 0, sqrt(METx*METx+METy*METy));
    int Idxj1_W=-1, Idxj2_W=-1, Idxlj1_W=-1, Idxlj2_W=-1;
      Idxj1_W=GetDaughterCandIdx(JetColl, "W", -1., "Lead");
      Idxj2_W=GetDaughterCandIdx(JetColl, "W", -1., "Subl");
    std::vector<snu::KJet> LJetColl = SelLightJets(JetColl, "Medium");
      if(LJetColl.size()>1){
        Idxlj1_W=GetDaughterCandIdx(LJetColl, "W", -1., "Lead");
        Idxlj2_W=GetDaughterCandIdx(LJetColl, "W", -1., "Subl");
      }

    float MT1       = sqrt(2)*sqrt(MET*MuTColl.at(IdxSS1).Pt()-METx*MuTColl.at(IdxSS1).Px()-METy*MuTColl.at(IdxSS1).Py());
    float MT2       = sqrt(2)*sqrt(MET*MuTColl.at(IdxSS2).Pt()-METx*MuTColl.at(IdxSS2).Px()-METy*MuTColl.at(IdxSS2).Py());
    float MTW       = sqrt(2)*sqrt(MET*MuTColl.at(IdxSSW).Pt()-METx*MuTColl.at(IdxSSW).Px()-METy*MuTColl.at(IdxSSW).Py());
    float MTA       = sqrt(2)*sqrt(MET*MuTColl.at(IdxSSA).Pt()-METx*MuTColl.at(IdxSSA).Px()-METy*MuTColl.at(IdxSSA).Py());
    float dPhivlW   = TVector2::Phi_mpi_pi(MuTColl.at(IdxSSW).Phi()- METphi);
    float dPhivlA   = TVector2::Phi_mpi_pi(MuTColl.at(IdxSSA).Phi()- METphi);
    float fdPhivlW  = fabs(dPhivlW); 
    float fdPhivlA  = fabs(dPhivlA);
    float dPt       = MuTColl.at(IdxSSW).Pt()-MuTColl.at(IdxSSA).Pt();
    float dMTWA     = MTW-MTA;
    float dMT_W70   = fabs(MTW-70);
    float dMT_A70   = fabs(MTA-70);

    float fdPhi_AWlv = fabs(TVector2::Phi_mpi_pi((MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)).Phi() - (v+MuTColl.at(IdxSSW)).Phi()));
    float fdPhi_AWjj = fabs(TVector2::Phi_mpi_pi((MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)).Phi() - (JetColl.at(Idxj1_W)+JetColl.at(Idxj2_W)).Phi()));
    float fdPhi_AWljlj = LJetColl.size()>1? fabs(TVector2::Phi_mpi_pi((MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)).Phi() - (LJetColl.at(Idxlj1_W)+LJetColl.at(Idxlj2_W)).Phi())):-1;
    float fdPhi_AWlv_F = fabs(TVector2::Phi_mpi_pi((MuTColl.at(IdxOS)+MuTColl.at(IdxSSW)).Phi() - (v+MuTColl.at(IdxSSA)).Phi()));
    float fdPhi_AWjj_F = fabs(TVector2::Phi_mpi_pi((MuTColl.at(IdxOS)+MuTColl.at(IdxSSW)).Phi() - (JetColl.at(Idxj1_W)+JetColl.at(Idxj2_W)).Phi()));
    float fdPhi_AWljlj_F = LJetColl.size()>1? fabs(TVector2::Phi_mpi_pi((MuTColl.at(IdxOS)+MuTColl.at(IdxSSW)).Phi() - (LJetColl.at(Idxlj1_W)+LJetColl.at(Idxlj2_W)).Phi())):-1;

    float dR_Ab1    = BJetColl.at(0).DeltaR( (MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)) );
    float dR_Ab1_F  = BJetColl.at(0).DeltaR( (MuTColl.at(IdxOS)+MuTColl.at(IdxSSW)) );
    float dR_AWjj   = (MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)).DeltaR(JetColl.at(Idxj1_W)+JetColl.at(Idxj2_W));
    float dR_AWjj_F = (MuTColl.at(IdxOS)+MuTColl.at(IdxSSW)).DeltaR(JetColl.at(Idxj1_W)+JetColl.at(Idxj2_W));
    float dR_AWljlj   = LJetColl.size()>1? (MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)).DeltaR(LJetColl.at(Idxlj1_W)+LJetColl.at(Idxlj2_W)):-1;
    float dR_AWljlj_F = LJetColl.size()>1? (MuTColl.at(IdxOS)+MuTColl.at(IdxSSW)).DeltaR(LJetColl.at(Idxlj1_W)+LJetColl.at(Idxlj2_W)):-1;


    FillHist("dR_b1lA", BJetColl.at(0).DeltaR(MuTColl.at(IdxSSA)), weight, 0., 5., 50);
    FillHist("dR_b1lW", BJetColl.at(0).DeltaR(MuTColl.at(IdxSSW)), weight, 0., 5., 50);
    FillHist("dR_j1lW", JetColl.at(0).DeltaR(MuTColl.at(IdxSSW)), weight, 0., 5., 50);
    FillHist("dR_j2lW", JetColl.at(1).DeltaR(MuTColl.at(IdxSSW)), weight, 0., 5., 50);
    FillHist("dR_j1lA", JetColl.at(0).DeltaR(MuTColl.at(IdxSSA)), weight, 0., 5., 50);
    FillHist("dR_j2lA", JetColl.at(1).DeltaR(MuTColl.at(IdxSSA)), weight, 0., 5., 50);

    FillHist("dR_Ab1"      ,       dR_Ab1, weight, 0.,   5., 50);
    FillHist("dR_Ab1_F"    ,     dR_Ab1_F, weight, 0.,   5., 50);
    FillHist("dR_AWjj"     ,      dR_AWjj, weight, 0.,   5., 50);
    FillHist("dR_AWjj_F"   ,    dR_AWjj_F, weight, 0.,   5., 50);
    FillHist("fdPhi_AWlv"  ,   fdPhi_AWlv, weight, 0., 3.15, 10);
    FillHist("fdPhi_AWjj"  ,   fdPhi_AWjj, weight, 0., 3.15, 10);
    FillHist("fdPhi_AWlv_F", fdPhi_AWlv_F, weight, 0., 3.15, 10);
    FillHist("fdPhi_AWjj_F", fdPhi_AWjj_F, weight, 0., 3.15, 10);

    if(LJetColl.size()>1){
      FillHist("fdPhi_AWjj_lj2"  ,   fdPhi_AWljlj, weight, 0., 3.15, 10);
      FillHist("fdPhi_AWjj_lj2_F", fdPhi_AWljlj_F, weight, 0., 3.15, 10);
      FillHist("dR_AWljlj"       ,      dR_AWljlj, weight, 0.,   5., 50);
      FillHist("dR_AWljlj_F"     ,    dR_AWljlj_F, weight, 0.,   5., 50);
    }
    if(MuTColl.at(IdxSSW).Charge()>0){
      FillHist("dR_b1lW_HcTommlv", BJetColl.at(0).DeltaR(MuTColl.at(IdxSSW)), weight, 0., 5., 50);
      FillHist("dR_j1lW_HcTommlv",  JetColl.at(0).DeltaR(MuTColl.at(IdxSSW)), weight, 0., 5., 50);
      FillHist("dR_j2lW_HcTommlv",  JetColl.at(1).DeltaR(MuTColl.at(IdxSSW)), weight, 0., 5., 50);

      FillHist("fdPhi_AWlv_HcTommlv", fdPhi_AWlv, weight, 0., 3.15, 10);
      FillHist("fdPhi_AWjj_HcTommlv", fdPhi_AWjj, weight, 0., 3.15, 10);
      if(JetColl.size()-BJetColl.size()>1){
        FillHist("fdPhi_AWjj_lj2_HcTommlv",  fdPhi_AWjj, weight, 0., 3.15, 10);
        FillHist("dR_AWljlj_HcTommlv"     ,   dR_AWljlj, weight, 0.,   5., 50);
        FillHist("dR_AWljlj_F_HcTommlv"   , dR_AWljlj_F, weight, 0.,   5., 50);
      }
    }
    else{
      FillHist("dR_b1lW_HcTommqq", BJetColl.at(0).DeltaR(MuTColl.at(IdxSSW)), weight, 0., 5., 50);
      FillHist("dR_j1lW_HcTommqq",  JetColl.at(0).DeltaR(MuTColl.at(IdxSSW)), weight, 0., 5., 50);
      FillHist("dR_j2lW_HcTommqq",  JetColl.at(1).DeltaR(MuTColl.at(IdxSSW)), weight, 0., 5., 50);

      FillHist("fdPhi_AWlv_HcTommqq", fdPhi_AWlv, weight, 0., 3.15, 10);
      FillHist("fdPhi_AWjj_HcTommqq", fdPhi_AWjj, weight, 0., 3.15, 10);
      if(JetColl.size()-BJetColl.size()>1){
        FillHist("fdPhi_AWjj_lj2_HcTommqq",  fdPhi_AWjj, weight, 0., 3.15, 10);
        FillHist("dR_AWljlj_HcTommqq"     ,   dR_AWljlj, weight, 0.,   5., 50);
        FillHist("dR_AWljlj_F_HcTommqq"   , dR_AWljlj_F, weight, 0.,   5., 50);
      }
    }


    //Single Property
    FillHist("MT_muA", MTA, weight, 0., 200., 40);
    FillHist("MT_muW", MTW, weight, 0., 200., 40);
    FillHist("fdPhi_muA", fdPhivlA, weight, 0., 3.15, 10);
    FillHist("fdPhi_muW", fdPhivlW, weight, 0., 3.15, 10);
    FillHist("dMTW70", dMT_W70, weight, 0., 100., 20);   
    FillHist("dMTA70", dMT_A70, weight, 0., 100., 20);   
//    FillHist("dPtmuAW_MTW_2D"   , dPt, MTW    , weight, -250., 250., 50,    0., 200., 20);   
//    FillHist("dPtmuAW_MTA_2D"   , dPt, MTA    , weight, -250., 250., 50,    0., 200., 20);   
//    FillHist("PtlA_PtlW_2D"     , MuTColl.at(IdxSSA).Pt(), MuTColl.at(IdxSSW).Pt(), weight, 0., 150., 30, 0., 150., 30);
//    FillHist("MTA_MTW_2D"       , MTA, MTW    , weight, -  0., 200., 20,    0., 200., 20);   


    //Relative Property
    FillHist("dPt_muAW" , dPt, weight, 0., 300, 60);
    FillHist("dPhi_muAW", fdPhivlW-fdPhivlA, weight, -3.15, 3.15, 20);
    FillHist("dMT", MTW-MTA, weight, -200., 200., 40);   
    FillHist("ddMT70", dMT_W70-dMT_A70, weight, -100., 100., 20);   
//    FillHist("dPtmuAW_dMTW70_2D", dPt, dMTWA            , weight, -250., 250., 50, -200., 200., 40);   
//    FillHist("dPtmuAW_ddMT70_2D", dPt, dMT_W70-dMT_A70  , weight, -250., 250., 50, -100., 100., 40);   
//    FillHist("dPtmuAW_dPhi_2D"  , dPt, fdPhivlW-fdPhivlA, weight, -250., 250., 50, -3.15, 3.15, 40);   

    FillHist("CompTest_Pt", -1., weight, -1., 2., 3);
    if(dPt>0) FillHist("CompTest_Pt", 0., weight, -1., 2., 3);
    else      FillHist("CompTest_Pt", 1., weight, -1., 2., 3);

    FillHist("CompTest_MT", -1., weight, -1., 2., 3);
    if(dMTWA>0) FillHist("CompTest_MT", 0., weight, -1., 2., 3);
    else        FillHist("CompTest_MT", 1., weight, -1., 2., 3);

    FillHist("CompTest_MT70", -1., weight, -1., 2., 3);
    if(dMT_W70<dMT_A70) FillHist("CompTest_MT70", 0., weight, -1., 2., 3);
    else                FillHist("CompTest_MT70", 1., weight, -1., 2., 3);

    FillHist("CompTest_Phi", -1., weight, -1., 2., 3);
    if(fdPhivlW>fdPhivlA) FillHist("CompTest_Phi", 0., weight, -1., 2., 3);
    else                  FillHist("CompTest_Phi", 1., weight, -1., 2., 3);

    FillHist("CompTest_dRlb", -1., weight, -1., 2., 3);
    if(BJetColl.at(0).DeltaR(MuTColl.at(IdxSSA))>BJetColl.at(0).DeltaR(MuTColl.at(IdxSSW))){
       FillHist("CompTest_dRlb", 0., weight, -1., 2., 3);
    }
    else FillHist("CompTest_dRlb", 1., weight, -1., 2., 3);


    //Ver2
    bool PhiAssign=true;
    if(PhiAssign){
      const int NCut_W1=7;  float WCut1[NCut_W1]={20., 30., 40., 50., 60., 70., 80.};
      const int NCut_W2=14; float WCut2[NCut_W2]={70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200.};
      const int NCut_dPt=9; float dPtCut[NCut_dPt]={0., 5., 10., 15., 20., 25., 30., 35., 40.};
      const int NCut_dRlb=8;float dRlbCut[NCut_dRlb]={0., 0.5, 1., 1.5, 2., 2.5, 3., 5.} ;
      //const int NCut_W1=7;  float WCut1[NCut_W1]={30., 35., 40., 45., 50., 55., 60.};
      //const int NCut_W2=9;  float WCut2[NCut_W2]={100., 105., 110., 115., 120., 125., 130., 135., 140.};
      //const int NCut_dPt=11; float dPtCut[NCut_dPt]={15., 17., 19., 21., 24., 25., 27., 29., 31., 33., 35.};
  
      int CutIdx=0;
      FillHist("Norm_MTWdPtdR", 0., weight, 0., 1., 1);
      for(int it_dpt=0; it_dpt<NCut_dPt; it_dpt++){
      for(int it_w1=0; it_w1<NCut_W1; it_w1++){
      for(int it_w2=0; it_w2<NCut_W2; it_w2++){
      for(int it_dRlb=0; it_dRlb<NCut_dRlb; it_dRlb++){
        bool CorrectAssign=false;
        if(fabs(MuTColl.at(IdxSS1).Pt()-MuTColl.at(IdxSS2).Pt())<dPtCut[it_dpt]){
          bool SS1_W= MT1>WCut1[it_w1] && MT1<WCut2[it_w2];
          bool SS2_W= MT2>WCut1[it_w1] && MT2<WCut2[it_w2];
          if     (SS2_W && (!SS1_W)){ if(IdxSS2==IdxSSW) CorrectAssign=true; }
          else if(SS1_W && (!SS2_W)){ if(IdxSS1==IdxSSW) CorrectAssign=true; }
          else{
           float dRbl1=BJetColl.at(0).DeltaR(MuTColl.at(IdxSS1));
           float dRbl2=BJetColl.at(0).DeltaR(MuTColl.at(IdxSS2));
           if     (dRbl1<dRlbCut[it_dRlb] && dRbl2>dRlbCut[it_dRlb]){ if(IdxSS1==IdxSSW) CorrectAssign=true; }
           else if(dRbl1>dRlbCut[it_dRlb] && dRbl2<dRlbCut[it_dRlb]){ if(IdxSS2==IdxSSW) CorrectAssign=true; }
           else  { if(IdxSS1==IdxSSW) CorrectAssign=true; }
//           if     (dRbl1<dRbl2){ if(IdxSS1==IdxSSW) CorrectAssign=true; }
//           else                { if(IdxSS2==IdxSSW) CorrectAssign=true; }
          }
        }
        else{ if(IdxSS1==IdxSSW) CorrectAssign=true; }
  
        if(CorrectAssign) FillHist("MTWdPtdRCutOpt", CutIdx, weight, 0., 8000., 8000);
        CutIdx++;
      }
      }
      }
      }
    }



    //Ver1
    bool dPTPtMTOpt=true;
    if(dPTPtMTOpt){
      //Ver1.1
      //const int NCut_W1=7;  float WCut1[NCut_W1]={30., 35., 40., 45., 50., 55., 60.};
      //const int NCut_W2=9;  float WCut2[NCut_W2]={100., 105., 110., 115., 120., 125., 130., 135., 140.};
      //const int NCut_dPt=11; float dPtCut[NCut_dPt]={15., 17., 19., 21., 24., 25., 27., 29., 31., 33., 35.};
      //Ver1.2
      const int NCut_W1=13;  float WCut1[NCut_W1]={20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80.};
      const int NCut_W2=14;  float WCut2[NCut_W2]={70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 1000.};
      const int NCut_dPt=13; float dPtCut[NCut_dPt]={0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60};
      //Ver1.Test
      //const int NCut_W1=7;  float WCut1[NCut_W1]={20., 30., 40., 50., 60., 70., 80.};
      //const int NCut_W2=14; float WCut2[NCut_W2]={70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200.};
      //const int NCut_dPt=9; float dPtCut[NCut_dPt]={0., 5., 10., 15., 20., 25., 30., 35., 40.};


      int CutIdx=0;
      FillHist("Norm_MTWdPt", 0., weight, 0., 1., 1);
      for(int it_dpt=0; it_dpt<NCut_dPt; it_dpt++){
      for(int it_w1=0; it_w1<NCut_W1; it_w1++){
      for(int it_w2=0; it_w2<NCut_W2; it_w2++){
        bool CorrectAssign=false;
        if(fabs(MuTColl.at(IdxSS1).Pt()-MuTColl.at(IdxSS2).Pt())<dPtCut[it_dpt]){
          bool SS1_W= MT1>WCut1[it_w1] && MT1<WCut2[it_w2];
          bool SS2_W= MT2>WCut1[it_w1] && MT2<WCut2[it_w2];
          if(SS2_W && (!SS1_W)){ if(IdxSS2==IdxSSW) CorrectAssign=true; }
          else{ if(IdxSS1==IdxSSW) CorrectAssign=true; }
        }
        else{ if(IdxSS1==IdxSSW) CorrectAssign=true; }
  
        if(CorrectAssign) FillHist("MTWdPtCutOpt", CutIdx, weight, 0., 3000., 3000);
        CutIdx++;
      }
      }
      }
    }

  

  }

}



void Aug2017_TriLepSR::CheckCutflow(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){

  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");
  std::vector<snu::KMuon> MuL10Coll;
    for(int i=0; i<MuLColl.size(); i++){ if(MuLColl.at(i).Pt()>10.){ MuL10Coll.push_back(MuLColl.at(i)); } }
  std::vector<snu::KMuon> MuL5Coll;
    for(int i=0; i<MuLColl.size(); i++){ if(MuLColl.at(i).Pt()>5.){   MuL5Coll.push_back(MuLColl.at(i)); } }
  std::vector<snu::KMuon> MuT10Coll;
    for(int i=0; i<MuTColl.size(); i++){ if(MuTColl.at(i).Pt()>10.){ MuT10Coll.push_back(MuTColl.at(i)); } }
  std::vector<snu::KMuon> MuT5Coll;
    for(int i=0; i<MuTColl.size(); i++){ if(MuTColl.at(i).Pt()>5.){   MuT5Coll.push_back(MuTColl.at(i)); } }

  if(EMuMu){
    if( !(EleLColl.size()==1 && EleTColl.size()==1)) return;
    if( !( MuLColl.size()>=2 &&  MuTColl.size()>=2)) return;
    if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10) ) return;
    if( fabs(EleTColl.at(0).Eta())>2.5 ) return;
      if(k_sample_name.Contains("DY") || k_sample_name.Contains("TT_powheg")){
        std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
        int LeptonType=0;
        LeptonType=GetLeptonType(MuTColl.at(1),truthColl);
        if(LeptonType>=-4 && LeptonType<0) goto flowstart_1e2mu;
        LeptonType=GetLeptonType(EleTColl.at(0),truthColl);
        if(LeptonType>=-4 && LeptonType<0) goto flowstart_1e2mu;
        LeptonType=GetLeptonType(MuTColl.at(0),truthColl);
        if(LeptonType>=-4 && LeptonType<0) goto flowstart_1e2mu;
        return;
      }
      flowstart_1e2mu:

    float CurrentIdx=0.;
    //Scenario that what if used 5GeV cut on 3rd lepton.
    if( MuL5Coll.size()==2 && MuT5Coll.size()==2 ){
      if( MuT5Coll.at(0).Charge()!=MuT5Coll.at(1).Charge() ){
        float Mmumu=(MuT5Coll.at(0)+MuT5Coll.at(1)).M();
        if(Mmumu>12) FillHist("CutFlow_1e2mu"+Label, CurrentIdx, weight, 0., 10., 10);
      }
    } CurrentIdx++;

    //Original scenario
    if( !(MuL10Coll.size()==2 && MuT10Coll.size()==2) ) return;
    if( !(MuT10Coll.at(0).Charge()!=MuT10Coll.at(1).Charge()) ) return;
      float Mmumu=(MuT10Coll.at(0)+MuT10Coll.at(1)).M();
    if(Mmumu<12) return;
    if( !(EleTColl.at(0).Pt()>25 && MuT10Coll.at(0).Pt()>10 && MuT10Coll.at(1).Pt()>10) ) return;
      FillHist("CutFlow_1e2mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//Analysis PT Cut
    if(fabs(Mmumu-91.2)<10) return;
      FillHist("CutFlow_1e2mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//OffZ
    if(BJetColl.size()==0) return;
      FillHist("CutFlow_1e2mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//HasB
    if(JetColl.size()<2) return;
      FillHist("CutFlow_1e2mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//geq2j
  }
  if(TriMu){
    if( !(EleLColl.size()==0) ) return;
    if( !(MuTColl.size()>=3 ) ) return;
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10) ) return;
      if(k_sample_name.Contains("DY") || k_sample_name.Contains("TT_powheg")){
        std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
        bool HasFake=false;
        for(int i=0; i<MuTColl.size(); i++){
          int LeptonType=GetLeptonType(MuTColl.at(i),truthColl);
          if(LeptonType>=-4 && LeptonType<0){ HasFake=true; break; }
        }
        if(!HasFake) return;
      }

    float CurrentIdx=0.;
    //Scenario that what if used 5GeV cut on 3rd lepton.
    if( MuL5Coll.size()==3 && MuT5Coll.size()==3 ){
      if( fabs(SumCharge(MuT5Coll))==1 ){
        int   IdxOS  = TriMuChargeIndex(MuT5Coll, "OS");
        int   IdxSS1 = TriMuChargeIndex(MuT5Coll, "SS1");
        int   IdxSS2 = TriMuChargeIndex(MuT5Coll, "SS2");
        float MOSSS1 = (MuT5Coll.at(IdxOS)+MuT5Coll.at(IdxSS1)).M();
        float MOSSS2 = (MuT5Coll.at(IdxOS)+MuT5Coll.at(IdxSS2)).M();

        if(MOSSS1>12 && MOSSS2>12){
          FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10);
        }
      }
    } CurrentIdx++;


    if( !(MuL10Coll.size()==3 && MuT10Coll.size()==3) ) return;
    if( !(MuT10Coll.at(0).Pt()>20 && MuT10Coll.at(1).Pt()>10 && MuT10Coll.at(2).Pt()>10) ) return;
    if( fabs(SumCharge(MuT10Coll))!=1 ) return;
      int   IdxOS  = TriMuChargeIndex(MuT10Coll, "OS");
      int   IdxSS1 = TriMuChargeIndex(MuT10Coll, "SS1");
      int   IdxSS2 = TriMuChargeIndex(MuT10Coll, "SS2");
      float MOSSS1 = (MuT10Coll.at(IdxOS)+MuT10Coll.at(IdxSS1)).M();
      float MOSSS2 = (MuT10Coll.at(IdxOS)+MuT10Coll.at(IdxSS2)).M();
    if( !(MOSSS1>12 && MOSSS2>12) ) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//3lep within simulation reg. & Triggered.
    if( fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10 ) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//3lep within simulation reg. & Triggered.
    if(BJetColl.size()==0) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//HasB
    if(JetColl.size()<2) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//geq2j
  }

}


void Aug2017_TriLepSR::CheckDecComp(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){

  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");
  if(EMuMu){
    if( !(EleLColl.size()==1 && MuLColl.size()==2) ) return;
    if( !(EleTColl.size()==1 && MuTColl.size()==2) ) return;
    if( !(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) ) return;
    if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10) ) return;
    if( fabs(EleTColl.at(0).Eta())>2.5 ) return;
      float CurrentIdx=0.;
      float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
    if(Mmumu<12) return;
      FillHist("CutFlow_1e2mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//3lep within simulation reg. & Triggered.
    if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && MuTColl.at(1).Pt()>10) ) return;
      FillHist("CutFlow_1e2mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//Analysis PT Cut
    if(fabs(Mmumu-91.2)<10) return;
      FillHist("CutFlow_1e2mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//OffZ
    if(BJetColl.size()==0) return;
      FillHist("CutFlow_1e2mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//HasB
    if(JetColl.size()<2) return;
      FillHist("CutFlow_1e2mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//geq2j



    std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
    int NtaW=0, NmuW=0, NeW=0, NlW=0;
    int NmuA=0;
    for(int i=2; i<truthColl.size(); i++){
      int fpid=fabs(truthColl.at(i).PdgId());
      if( !(fpid==11 || fpid==13 || fpid==15) ) continue;
      int fmpid=fabs(truthColl.at(truthColl.at(i).IndexMother()).PdgId());
      if( fmpid==24 ){
        if     (fpid==11) NeW++;
        else if(fpid==13) NmuW++;
        else if(fpid==15) NtaW++;
      }
      if( fpid==13 && fmpid==36 ){ NmuW++; }
    }
    NlW = NeW + NmuW + NtaW;
    //if(NlWlta==0) PrintTruth();

    if( fabs(Mmumu-15.)<0.7 ) FillHist("NlWltaW_M15"+Label, NlW, weight, 0., 10., 10);
    if( fabs(Mmumu-25.)<1.  ) FillHist("NlWltaW_M25"+Label, NlW, weight, 0., 10., 10);
    if( fabs(Mmumu-35.)<1.3 ) FillHist("NlWltaW_M35"+Label, NlW, weight, 0., 10., 10);
    if( fabs(Mmumu-45.)<1.5 ) FillHist("NlWltaW_M45"+Label, NlW, weight, 0., 10., 10);
    if( fabs(Mmumu-55.)<1.5 ) FillHist("NlWltaW_M55"+Label, NlW, weight, 0., 10., 10);
    if( fabs(Mmumu-65.)<2.  ) FillHist("NlWltaW_M65"+Label, NlW, weight, 0., 10., 10);
    if( fabs(Mmumu-75.)<2.5 ) FillHist("NlWltaW_M75"+Label, NlW, weight, 0., 10., 10);
  }
  if(TriMu){
    if( !(MuLColl.size()==3 && EleLColl.size()==0) ) return;
    if( !(MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 && MuTColl.at(2).Pt()>5) ) return;
      if(k_sample_name.Contains("DY") || k_sample_name.Contains("TT_powheg")){
        std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
        int LeptonType=0;
        LeptonType=GetLeptonType(MuTColl.at(2),truthColl);
        if(LeptonType>=-4 && LeptonType<0) goto flowstart_3mu;
        LeptonType=GetLeptonType(MuTColl.at(1),truthColl);
        if(LeptonType>=-4 && LeptonType<0) goto flowstart_3mu;
        LeptonType=GetLeptonType(MuTColl.at(0),truthColl);
        if(LeptonType>=-4 && LeptonType<0) goto flowstart_3mu;
        return;
      }
      flowstart_3mu:
      float CurrentIdx=0.;
      int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
      int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
      int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
      float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
      float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    if( !(MOSSS1>12 && MOSSS2>12) ) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//3lep within simulation reg. & Triggered.
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 && MuTColl.at(2).Pt()>10) ) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//3lep within simulation reg. & Triggered.
    if( fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10 ) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//3lep within simulation reg. & Triggered.
    if(BJetColl.size()==0) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//HasB
    if(JetColl.size()<2) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//geq2j
  }

}



void Aug2017_TriLepSR::CheckSRYield(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");

  const int Ncent=7;
  float Center[Ncent]      = { 15., 25., 35., 45., 55., 65., 75.};
  float Width[Ncent]       = { 1., 1. , 1.3, 1.5, 1.5,  2., 2.5};
  TString CentStr[Ncent]   = { "15", "25", "35", "45", "55", "65", "75"};


  const int NSignal = 7;
  float Mmumu_CentSeed[NSignal]  = {15., 25., 35., 45., 55., 65., 75.};
  float Mmumu_WidthSeed[NSignal] = {1., 1., 1.3, 1.5, 1.5, 2.0, 2.5};
  float Mmumu_Step[NSignal-1]    = {1., 1., 1.2, 1.2, 1.2, 1.8};
  vector<float> Mmumu_CentFull, Mmumu_WidthFull;
  for(int it_seed=0; it_seed<NSignal; it_seed++){
    float CurrentCenter = Mmumu_CentSeed[it_seed];
    float CurrentWidth  = Mmumu_WidthSeed[it_seed];
    Mmumu_CentFull.push_back(CurrentCenter);
    Mmumu_WidthFull.push_back(CurrentWidth);
    for(int it_step=0; it_step<100; it_step++){
      if(it_seed==NSignal-1) break;
      CurrentCenter+= Mmumu_Step[it_seed];
      CurrentWidth  = (CurrentCenter-Mmumu_CentSeed[it_seed])/10.*Mmumu_WidthSeed[it_seed+1]
                     +(Mmumu_CentSeed[it_seed+1]-CurrentCenter)/10.*Mmumu_WidthSeed[it_seed];
      if(CurrentCenter<Mmumu_CentSeed[it_seed]+10.){
        Mmumu_CentFull.push_back(CurrentCenter);
        Mmumu_WidthFull.push_back(CurrentWidth);
      }
      else break;
    }
  }

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

    if(fabs(Mmumu-91.2)<10) return;

    float MAWidth=2.5;
    //V1:  1./  1./ 1.5/ 1.5/ 1.5/ 2./ 2.5
    //V2: 0.7/  1./ 1.3/ 1.5/ 1.5/ 2./ 2.5
    //V3: 0.7/ 0.9/ 0.9/ 1.5/ 1.5/ 2./ 2.5
    //V4:  1./  1./ 1.3/ 1.5/ 1.5/ 2./ 2.5
    if(JetColl.size()>1){
      if(Label=="") FillHist("Nb_2jCutPreB_1e2mu"+Label, BJetColl.size(), weight, 0., 10., 10);

      FillHist("Count_2jCutPreB_1e2mu_PreSel"+Label, 0., weight, 0., 1., 1);
      if( Mmumu>12. && Mmumu<80 ) FillHist("Count_2jCutPreB_1e2mu_M12to80"+Label, 0., weight, 0., 1., 1);

      for( int it_m=0; it_m<Ncent; it_m++){
        float Shift  = 0.;
        if( fabs(Mmumu-Center[it_m]-Shift)<Width[it_m] ){
          FillHist("Count_2jCutPreB_1e2mu_M"+CentStr[it_m]+Label, 0., weight, 0., 1., 1);
        }
      }
    }

    if(BJetColl.size()==0) return;


    if(Label==""){
      FillHist("Nj_PreSel_1e2mu"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Mmumu_PreSel_1e2mu"+Label, Mmumu, weight, 0., 200., 40);
    }
    //c.f. Regarding detector resolution not considering FSR tail
    //σ(m, pT) = 0.026 + 0.0065m for the barrel,
    //σ(m, pT) = 0.026 + 0.013m GeV/c2 for the endcap.
    if(JetColl.size()>1){
      if(Label=="") FillHist("Mmumu_2jCut_1e2mu"+Label, Mmumu, weight, 0., 200., 40);

      FillHist("Count_2jCut_1e2mu_PreSel"+Label, 0., weight, 0., 1., 1);
      if( Mmumu>12. && Mmumu<80.  ) FillHist("Count_2jCut_1e2mu_M12to80"+Label, 0., weight, 0., 1., 1);

      for( int it_m=0; it_m<Ncent; it_m++){
        float Shift = 0.;
        if( fabs(Mmumu-Center[it_m]-Shift)<Width[it_m] ){
          FillHist("Count_2jCut_M"+CentStr[it_m]+Label, 0., weight, 0., 1., 1);
          //FillHist("Count_2jCut_1e2mu_M"+CentStr[it_m]+Label, 0., weight, 0., 1., 1);
        }
      }

      for(int it_cent=0; it_cent<Mmumu_CentFull.size(); it_cent++){
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<Mmumu_WidthFull.at(it_cent)){
          FillHist("Mmumu_CntBin_SR2j"+Label, it_cent+0.01, weight, 0., 54., 54);
        }
  
        float LargeWidth=5.;
        if((Mmumu_CentFull.at(it_cent)-LargeWidth)<12.) LargeWidth=Mmumu_CentFull.at(it_cent)-12.;
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidth){
          FillHist("Mmumu_Width10_SR2j"+Label, it_cent+0.01, weight, 0., 54., 54);
        }
        if(Label==""){
          float Variation=0.2;
          float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
          if((Mmumu_CentFull.at(it_cent)-LargeWidthUp)<12.) LargeWidthUp=Mmumu_CentFull.at(it_cent)-12.;
          if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthUp){
            FillHist("Mmumu_Width10_SR2j_systup_Mrange", it_cent+0.01, weight, 0., 54., 54);
          }
          if((Mmumu_CentFull.at(it_cent)-LargeWidthDown)<12.) LargeWidthDown=Mmumu_CentFull.at(it_cent)-12.;
          if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthDown){
            FillHist("Mmumu_Width10_SR2j_systdown_Mrange", it_cent+0.01, weight, 0., 54., 54);
          }
        }
      }
    }
  }
  if(TriMu){
    if( !(MuLColl.size()==3 && EleLColl.size()==0) ) return;
    if( !(MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>20. && MuTColl.at(1).Pt()>10. && MuTColl.at(2).Pt()>10.) ) return;

    int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    int   IdxSSA = TriMuChargeIndex(MuTColl, MET, METx, METy, "SSA");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    //float Mmumu  = MOSSS2;
    float Mmumu  = (MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)).M();
    float MAWidth=2.5;
    if( !(MOSSS1>12. && MOSSS2>12.) ) return;
    if( fabs(MOSSS1-91.2)<10. || fabs(MOSSS2-91.2)<10. ) return;

    if(JetColl.size()>1){
      if(Label=="") FillHist("Nb_2jCutPreB_3mu"+Label, BJetColl.size(), weight, 0., 10., 10);

      FillHist("Count_2jCutPreB_Mmumu_3mu_PreSel"+Label, 0., weight, 0., 1., 1);
      if(Mmumu<80.) FillHist("Count_2jCutPreB_Mmumu_3mu_M12to80"+Label, 0., weight, 0., 1., 1);

      for( int it_m=0; it_m<Ncent; it_m++){
        float Shift = 0.;
        if( fabs(Mmumu-Center[it_m]-Shift)<Width[it_m] ){
          FillHist("Count_2jCutPreB_Mmumu_3mu_M"+CentStr[it_m]+Label, 0., weight, 0., 1., 1);
        }
      }
    }

    if( BJetColl.size()==0 ) return;

    if(JetColl.size()>1){
      if(Label==""){
        FillHist("MOSSS1_2j"+Label, MOSSS1, weight, 0., 200., 10);
        FillHist("MOSSS2_2j"+Label, MOSSS2, weight, 0., 200., 20);
        FillHist("Mmumu_2j" +Label, Mmumu , weight, 0., 200., 20);
      }
      FillHist("Count_2jCut_Mmumu_3mu_PreSel"+Label, 0., weight, 0., 1., 1);
      if(Mmumu<80) FillHist("Count_2jCut_Mmumu_3mu_M12to80"+Label, 0., weight, 0., 1., 1);

      for( int it_m=0; it_m<Ncent; it_m++){
        float Shift = 0.;
        if( fabs(Mmumu-Center[it_m]-Shift)<Width[it_m] ){
          FillHist("Count_2jCut_M"+CentStr[it_m]+Label, 0., weight, 0., 1., 1);
          //FillHist("Count_2jCut_Mmumu_3mu_M"+CentStr[it_m]+Label, 0., weight, 0., 1., 1);
        }
      }

      for(int it_cent=0; it_cent<Mmumu_CentFull.size(); it_cent++){
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<Mmumu_WidthFull.at(it_cent)){
          FillHist("Mmumu_CntBin_SR2j"+Label, it_cent+0.01, weight, 0., 54., 54);
        }
  
        float LargeWidth=5.;
        if((Mmumu_CentFull.at(it_cent)-LargeWidth)<12.) LargeWidth=Mmumu_CentFull.at(it_cent)-12.;
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidth){
          FillHist("Mmumu_Width10_SR2j"+Label, it_cent+0.01, weight, 0., 54., 54);
        }
        if(Label==""){
          float Variation=0.2;
          float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
          if((Mmumu_CentFull.at(it_cent)-LargeWidthUp)<12.) LargeWidthUp=Mmumu_CentFull.at(it_cent)-12.;
          if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthUp){
            FillHist("Mmumu_Width10_SR2j_systup_Mrange", it_cent+0.01, weight, 0., 54., 54);
          }
          if((Mmumu_CentFull.at(it_cent)-LargeWidthDown)<12.) LargeWidthDown=Mmumu_CentFull.at(it_cent)-12.;
          if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthDown){
            FillHist("Mmumu_Width10_SR2j_systdown_Mrange", it_cent+0.01, weight, 0., 54., 54);
          }
        }
      }
    }
  }
}


void Aug2017_TriLepSR::CheckSigBkgKinematics(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool TriMu=false, EMuMu=false;
  int IdxOS=-1, IdxSS1=-1, IdxSS2=-1;
  float MOSSS1=0., MOSSS2=0., Mmumu=0.;
  std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
  int IdxA1=-1, IdxA2=-1, CountA=0;
  for(int i=2; i<truthColl.size(); i++){
    if(truthColl.at(i).GenStatus()!=1) continue;
    int fpid=fabs(truthColl.at(i).PdgId());
    if( fpid!=11 && fpid!=13 ) continue;
    if( (fpid==11 && fabs(truthColl.at(i).Eta())>2.5) || (fpid==13 && fabs(truthColl.at(i).Eta())>2.4) ) continue;

    int LepType=GetLeptonType(i,truthColl);
    if(LepType==1) FillHist("PTlW_Gen", truthColl.at(i).Pt(), weight, 0., 200., 40);
    if(LepType==2){
      if(CountA==0){ IdxA1=i; CountA++; }
      else{ if(truthColl.at(i).Pt()>truthColl.at(IdxA1).Pt()){ IdxA2=IdxA1; IdxA1=i; CountA++; }
            else{ IdxA2=i; CountA++; }
      }
      FillHist("PTlA_Gen", truthColl.at(i).Pt(), weight, 0., 200., 40);
    }
  }
  if(IdxA1!=-1) FillHist("PTlA1_Gen", truthColl.at(IdxA1).Pt(), weight, 0., 200., 40);
  if(IdxA2!=-1) FillHist("PTlA2_Gen", truthColl.at(IdxA2).Pt(), weight, 0., 200., 100);

  if( MuTColl.size()==3 && MuLColl.size()==3 && EleTColl.size()==0 && EleLColl.size()==0 ) TriMu=true;
  if( MuTColl.size()==2 && MuLColl.size()==2 && EleTColl.size()==1 && EleLColl.size()==1 ) EMuMu=true;
  if( TriMu && fabs(SumCharge(MuTColl))!=1 ) TriMu=false;
  if( EMuMu && SumCharge(MuTColl)!=0       ) EMuMu=false;
  if( TriMu ){
    IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    Mmumu  = MOSSS2;
  }
  if( EMuMu ){ Mmumu = (MuTColl.at(0)+MuTColl.at(1)).M(); }
  if( EMuMu && !(Mmumu>12 && fabs(Mmumu-91.2)>10) ) EMuMu=false;
  if( TriMu && !(MOSSS1>12 && MOSSS2>12 && fabs(MOSSS1-91.2)>10 && fabs(MOSSS2-91.2)>10) ) TriMu=false;
  if( !(TriMu || EMuMu) ) return;

  std::vector<snu::KMuon>     MuAColl, MuWColl,  MuFkColl; 
  std::vector<snu::KElectron>         EleWColl, EleFkColl; 
    for(int i=0; i<(int) MuTColl.size(); i++){
      int LepType=GetLeptonType(MuTColl.at(i), truthColl);
      if     (LepType==2) MuAColl.push_back(MuTColl.at(i));
      else if(LepType==1) MuWColl.push_back(MuTColl.at(i));
      else if(LepType<0 ) MuFkColl.push_back(MuTColl.at(i));
    }
    for(int i=0; i<(int) EleTColl.size(); i++){
      int LepType=GetLeptonType(EleTColl.at(i), truthColl);
      if     (LepType==1) EleWColl.push_back(EleTColl.at(i));
      else if(LepType<0 ) EleFkColl.push_back(EleTColl.at(i));
    }


  for(int i=0; i<(int) EleWColl.size(); i++){
    FillHist("PT_lW", EleWColl.at(i).Pt(), weight, 0., 200., 40);
    float MTW = sqrt(2)*sqrt(MET*EleWColl.at(i).Pt()-METx*EleWColl.at(i).Px()-METy*EleWColl.at(i).Py());
    FillHist("MTW_lW", MTW, weight, 0., 200., 40);
  }
  for(int i=0; i<(int) MuWColl.size(); i++){
    FillHist("PT_lW", MuWColl.at(i).Pt(), weight, 0., 200., 40);
    float MTW = sqrt(2)*sqrt(MET*MuWColl.at(i).Pt()-METx*MuWColl.at(i).Px()-METy*MuWColl.at(i).Py());
    FillHist("MTW_lW", MTW, weight, 0., 200., 40);
  }
  for(int i=0; i<(int) MuAColl.size(); i++){
    FillHist("PT_lA", MuAColl.at(i).Pt(), weight, 0., 200., 40);
    float MTW = sqrt(2)*sqrt(MET*MuAColl.at(i).Pt()-METx*MuAColl.at(i).Px()-METy*MuAColl.at(i).Py());
    FillHist("MTW_lA", MTW, weight, 0., 200., 40);
    if(i==0) FillHist("PT_lA1", MuAColl.at(i).Pt(), weight, 0., 200., 40);
    if(i==1) FillHist("PT_lA2", MuAColl.at(i).Pt(), weight, 0., 200., 100);
  }
  for(int i=0; i<(int) MuFkColl.size(); i++){
    FillHist("PT_MuFk", MuFkColl.at(i).Pt(), weight, 0., 200., 40);
  }
  for(int i=0; i<(int) EleFkColl.size(); i++){
    FillHist("PT_EleFk", EleFkColl.at(i).Pt(), weight, 0., 200., 40);
  }
  if(EMuMu){
    FillHist("PT_Ele_1e2mu", EleTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("PT_Mu1_1e2mu", MuTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("PT_Mu2_1e2mu", MuTColl.at(1).Pt(), weight, 0., 200., 100);
  }
  if(TriMu){
    FillHist("PT_Mu1_3mu", MuTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("PT_Mu2_3mu", MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("PT_Mu3_3mu", MuTColl.at(2).Pt(), weight, 0., 200., 100);
  }

}


void Aug2017_TriLepSR::CheckSRDist(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");
  if(k_sample_name.Contains("TTToHcToWA") && Label!="") return;//No systrun for signals

  vector<float> MmumuEdges;
  float CurrentEdge=0.;
   for(int i=0; i<19; i++){
     if     (i==0)  CurrentEdge=12.;
     else if(i==1)  CurrentEdge=20.;
     else if(i<7)   CurrentEdge+=10.;
     else if(i==7)  CurrentEdge=81.2;
     else if(i==8)  CurrentEdge=101.2;     
     else if(i==9)  CurrentEdge=110.;
     else if(i<18)  CurrentEdge+=10;
     else CurrentEdge=200.;
     MmumuEdges.push_back(CurrentEdge);
   }

  const int NSignal = 7;
  float Mmumu_CentSeed[NSignal]  = {15., 25., 35., 45., 55., 65., 75.};
  float Mmumu_WidthSeed[NSignal] = {1., 1., 1.3, 1.5, 1.5, 2.0, 2.5};
  float Mmumu_Step[NSignal-1]    = {1., 1., 1.2, 1.2, 1.2, 1.8};
  vector<float> Mmumu_CentFull, Mmumu_WidthFull;
  for(int it_seed=0; it_seed<NSignal; it_seed++){
    float CurrentCenter = Mmumu_CentSeed[it_seed];
    float CurrentWidth  = Mmumu_WidthSeed[it_seed];
    Mmumu_CentFull.push_back(CurrentCenter);
    Mmumu_WidthFull.push_back(CurrentWidth);
    for(int it_step=0; it_step<100; it_step++){
      if(it_seed==NSignal-1) break;
      CurrentCenter+= Mmumu_Step[it_seed];
      CurrentWidth  = (CurrentCenter-Mmumu_CentSeed[it_seed])/10.*Mmumu_WidthSeed[it_seed+1]
                     +(Mmumu_CentSeed[it_seed+1]-CurrentCenter)/10.*Mmumu_WidthSeed[it_seed];
      if(CurrentCenter<Mmumu_CentSeed[it_seed]+10.){
        Mmumu_CentFull.push_back(CurrentCenter);
        Mmumu_WidthFull.push_back(CurrentWidth);
      }
      else break;
    }
  }


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

      FillHist("CutFlow_1e2mu"+Label, 0., weight, 0., 20., 20);
      if(Label=="") FillHist("Mmumu_3lOS"+Label, Mmumu, weight, 0., 200., 200);
      if(Label==""){ if(JetColl.size()>1) FillHist("Mmumu_3lOS2j"+Label, Mmumu, weight, 0., 200., 200); }
    if(fabs(Mmumu-91.2)<10) return;
      FillHist("CutFlow_1e2mu"+Label, 1., weight, 0., 20., 20);
      FillHist("Nb_3lOSOffZ"+Label, BJetColl.size(), weight, 0., 5., 5);
      //Only for shape study
      if(Label=="") FillHist("Mmumu_3lOSOffZ"+Label, Mmumu, weight, 0., 200., 200);
      for(int it_cent=0; it_cent<Mmumu_CentFull.size(); it_cent++){
        if(Label!="") continue;
        float LargeWidth=5.;
        if((Mmumu_CentFull.at(it_cent)-LargeWidth)<12.) LargeWidth=Mmumu_CentFull.at(it_cent)-12.;
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidth){
          FillHist("Mmumu_CntBin_3lOSOffZ_Width10"+Label, it_cent+0.01, weight, 0., 54., 54);
        }
      }//----------------

    if(BJetColl.size()==0) return;
      FillHist("CutFlow_1e2mu"+Label, 2., weight, 0., 20., 20);
      FillHist("Nj_3lOSOffZHasB"+Label, JetColl.size(), weight, 0., 10., 10);
    if(JetColl.size()<2) return;
      FillHist("CutFlow_1e2mu"+Label, 3., weight, 0., 20., 20);

    //if(isData && !k_running_nonprompt && Mmumu<81.2) return;
    if(!isData && k_sample_name.Contains("TT_powheg")){
      std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
      int NFk=0;
      for(int it_lep=0; it_lep<EleTColl.size(); it_lep++){
        int LepType = GetLeptonType(EleTColl.at(it_lep),truthColl);
        if(LepType<0 && LepType>-5) NFk++;
      }
      for(int it_lep=0; it_lep<MuTColl.size(); it_lep++){
        int LepType = GetLeptonType(MuTColl.at(it_lep),truthColl);
        if(LepType<0 && LepType>-5) NFk++;
      }
      if(NFk<3) return;
    }

    int Idxj1_W=-1, Idxj2_W=-1;
    snu::KParticle v; v.SetPxPyPzE(METx, METy, 0, sqrt(METx*METx+METy*METy));
    float Pzv1=GetvPz(v, EleTColl.at(0), 1), Pzv2=GetvPz(v, EleTColl.at(0), 2);
    float Pzv=fabs(Pzv1)<fabs(Pzv2)? Pzv1:Pzv2;
    v.SetXYZM(METx,METy,Pzv,0.);

    float M3lv=(EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)+v).M(), M2l2j=0.; 
    Idxj1_W=GetDaughterCandIdx(JetColl, "W", -1., "Lead");
    Idxj2_W=GetDaughterCandIdx(JetColl, "W", -1., "Subl");
    M2l2j=(MuTColl.at(0)+MuTColl.at(1)+JetColl.at(Idxj1_W)+JetColl.at(Idxj2_W)).M();

    FillHist("Nj_SR2j"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_SR2j"+Label, BJetColl.size(), weight, 0., 5., 5);
    FillHist("PTe_SR2j"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("PTmu1_SR2j"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
    FillHist("PTmu2_SR2j"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("Mmumu_BinW10_SR2j"+Label, Mmumu, weight, &MmumuEdges[0], MmumuEdges.size()-1);
    FillHist("Mmumu_BinW2_SR2j" +Label, Mmumu, weight, 12., 80., 34);
    for(int it_cent=0; it_cent<34; it_cent++){
      float BinCenter=13.+2.*((float) it_cent);
      float LargeWidth=5.;
      if     ((BinCenter-LargeWidth)<12.)  LargeWidth=BinCenter-12.;
      else if((BinCenter+LargeWidth)>81.2) LargeWidth=81.2-BinCenter;
      if(fabs(Mmumu-BinCenter)<LargeWidth){
        FillHist("Mmumu_BinW2_SR2j_Width10"+Label, BinCenter, weight, 12., 80., 34);
      }
      if(Label==""){
        float Variation=0.2;
        float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
        if     ((BinCenter-LargeWidthUp)<12. ) LargeWidthUp=BinCenter-12.;
        else if((BinCenter+LargeWidthUp)>81.2) LargeWidthUp=81.2-BinCenter;
        if(fabs(Mmumu-BinCenter)<LargeWidthUp){
          FillHist("Mmumu_BinW2_SR2j_Width10_systup_Mrange", BinCenter, weight, 12., 80., 34);
        }
        if((BinCenter-LargeWidthDown)<12.) LargeWidthDown=BinCenter-12.;
        if(fabs(Mmumu-BinCenter)<LargeWidthDown){
          FillHist("Mmumu_BinW2_SR2j_Width10_systdown_Mrange", BinCenter, weight, 12., 80., 34);
        }
      }
    }
    for(int it_cent=0; it_cent<Mmumu_CentFull.size(); it_cent++){
      if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<Mmumu_WidthFull.at(it_cent)){
        FillHist("Mmumu_CntBin_SR2j"+Label, it_cent+0.01, weight, 0., 54., 54);
      }
      float LargeWidth=5.;
      if((Mmumu_CentFull.at(it_cent)-LargeWidth)<12.) LargeWidth=Mmumu_CentFull.at(it_cent)-12.;
      if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidth){
        FillHist("Mmumu_CntBin_SR2j_Width10"+Label, it_cent+0.01, weight, 0., 54., 54);
      }
      if(Label==""){
        float Variation=0.2;
        float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
        if((Mmumu_CentFull.at(it_cent)-LargeWidthUp)<12.) LargeWidthUp=Mmumu_CentFull.at(it_cent)-12.;
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthUp){
          FillHist("Mmumu_CntBin_SR2j_Width10_systup_Mrange", it_cent+0.01, weight, 0., 54., 54);
        }
        if((Mmumu_CentFull.at(it_cent)-LargeWidthDown)<12.) LargeWidthDown=Mmumu_CentFull.at(it_cent)-12.;
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthDown){
          FillHist("Mmumu_CntBin_SR2j_Width10_systdown_Mrange", it_cent+0.01, weight, 0., 54., 54);
        }
      }
    }
    FillHist("Mmumu_SR2j"+Label, Mmumu, weight, 0., 200., 200);
    FillHist("M2l2j_SR2j"+Label, M2l2j, weight, 0., 1000., 100);
    FillHist("M3lv_SR2j"+Label, M3lv, weight, 0., 1000., 100);
    if(Mmumu<80){
      FillHist("M2l2j_lt80_SR2j"+Label, M2l2j, weight, 0., 1000., 100);
      FillHist("M3lv_lt80_SR2j"+Label, M3lv, weight, 0., 1000., 100);
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
    int   IdxSSW = TriMuChargeIndex(MuTColl, MET, METx, METy, "SSW");
    int   IdxSSA = TriMuChargeIndex(MuTColl, MET, METx, METy, "SSA");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    float Mmumu  = (MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)).M();

    if( !(MOSSS1>12 && MOSSS2>12) ) return;

      FillHist("CutFlow_3mu"+Label, 0., weight, 0., 20., 20);
      if(Label=="") FillHist("Mmumu_3lOS"+Label, Mmumu, weight, 0., 200., 200);
    if( fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10 ) return;
      FillHist("CutFlow_3mu"+Label, 1., weight, 0., 20., 20);
      FillHist("Nb_3lOSOffZ"+Label, BJetColl.size(), weight, 0., 5., 5);
      //Only for shape studies
      if(Label=="") FillHist("Mmumu_3lOSOffZ"+Label, Mmumu, weight, 0., 200., 200);
      for(int it_cent=0; it_cent<Mmumu_CentFull.size(); it_cent++){
        if(Label!="") continue;
        float LargeWidth=5.;
        if((Mmumu_CentFull.at(it_cent)-LargeWidth)<12.) LargeWidth=Mmumu_CentFull.at(it_cent)-12.;
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidth){
          FillHist("Mmumu_CntBin_3lOSOffZ_Width10"+Label, it_cent+0.01, weight, 0., 54., 54);
        }
      }//----------------------

    if(BJetColl.size()==0) return;
      FillHist("CutFlow_3mu"+Label, 2., weight, 0., 20., 20);
      FillHist("Nj_3lOSOffZHasB"+Label, JetColl.size(), weight, 0., 10., 10);
    if(JetColl.size()<2) return;
      FillHist("CutFlow_3mu"+Label, 3., weight, 0., 20., 20);

    //if(isData && !k_running_nonprompt && Mmumu<81.2) return;
    if(!isData && k_sample_name.Contains("TT_powheg")){
      std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
      int NFk=0;
      for(int it_lep=0; it_lep<MuTColl.size(); it_lep++){
        int LepType = GetLeptonType(MuTColl.at(it_lep),truthColl);
        if(LepType<0 && LepType>-5) NFk++;
      }
      if(NFk<3) return;
    }


    int Idxj1_W=-1, Idxj2_W=-1;
    snu::KParticle v; v.SetPxPyPzE(METx, METy, 0, sqrt(METx*METx+METy*METy));
    float Pzv1=GetvPz(v, MuTColl.at(IdxSSW), 1), Pzv2=GetvPz(v, MuTColl.at(IdxSSW), 2);
    float Pzv=fabs(Pzv1)<fabs(Pzv2)? Pzv1:Pzv2;
    v.SetXYZM(METx,METy,Pzv,0.);

    float M3lv=(MuTColl.at(0)+MuTColl.at(1)+MuTColl.at(2)+v).M(), M2l2j=0.; 
    Idxj1_W=GetDaughterCandIdx(JetColl, "W", -1., "Lead");
    Idxj2_W=GetDaughterCandIdx(JetColl, "W", -1., "Subl");
    M2l2j=(MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)+JetColl.at(Idxj1_W)+JetColl.at(Idxj2_W)).M();

    FillHist("Nj_SR2j"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_SR2j"+Label, BJetColl.size(), weight, 0., 5., 5);
    FillHist("PTmu1_SR2j"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
    FillHist("PTmu2_SR2j"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("PTmu3_SR2j"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
    FillHist("MOSSS1_SR2j"+Label, MOSSS1, weight, 0., 200., 80);
    FillHist("Mmumu_BinW10_SR2j"+Label, Mmumu, weight, &MmumuEdges[0], MmumuEdges.size()-1);
    FillHist("Mmumu_BinW2_SR2j" +Label, Mmumu, weight, 12., 80., 34);
    for(int it_cent=0; it_cent<34; it_cent++){
      float BinCenter=13.+2.*((float) it_cent);
      float LargeWidth=5.;
      if     ((BinCenter-LargeWidth)<12. ) LargeWidth=BinCenter-12.;
      else if((BinCenter+LargeWidth)>81.2) LargeWidth=81.2-BinCenter;
      if(fabs(Mmumu-BinCenter)<LargeWidth){
        FillHist("Mmumu_BinW2_SR2j_Width10"+Label, BinCenter, weight, 12., 80., 34);
      }
      if(Label==""){
        float Variation=0.2;
        float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
        if     ((BinCenter-LargeWidthUp)<12. ) LargeWidthUp=BinCenter-12.;
        else if((BinCenter+LargeWidthUp)>81.2) LargeWidthUp=81.2-BinCenter;
        if(fabs(Mmumu-BinCenter)<LargeWidthUp){
          FillHist("Mmumu_BinW2_SR2j_Width10_systup_Mrange", BinCenter, weight, 12., 80., 34);
        }
        if((BinCenter-LargeWidthDown)<12.) LargeWidthDown=BinCenter-12.;
        if(fabs(Mmumu-BinCenter)<LargeWidthDown){
          FillHist("Mmumu_BinW2_SR2j_Width10_systdown_Mrange", BinCenter, weight, 12., 80., 34);
        }
      }
    }
    for(int it_cent=0; it_cent<Mmumu_CentFull.size(); it_cent++){
      if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<Mmumu_WidthFull.at(it_cent)){
        FillHist("Mmumu_CntBin_SR2j"+Label, it_cent+0.01, weight, 0., 54., 54);
      }
      float LargeWidth=5.;
      if((Mmumu_CentFull.at(it_cent)-LargeWidth)<12.) LargeWidth=Mmumu_CentFull.at(it_cent)-12.;
      if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidth){
        FillHist("Mmumu_CntBin_SR2j_Width10"+Label, it_cent+0.01, weight, 0., 54., 54);
      }
      if(Label==""){
        float Variation=0.2;
        float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
        if((Mmumu_CentFull.at(it_cent)-LargeWidthUp)<12.) LargeWidthUp=Mmumu_CentFull.at(it_cent)-12.;
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthUp){
          FillHist("Mmumu_CntBin_SR2j_Width10_systup_Mrange", it_cent+0.01, weight, 0., 54., 54);
        }
        if((Mmumu_CentFull.at(it_cent)-LargeWidthDown)<12.) LargeWidthDown=Mmumu_CentFull.at(it_cent)-12.;
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthDown){
          FillHist("Mmumu_CntBin_SR2j_Width10_systdown_Mrange", it_cent+0.01, weight, 0., 54., 54);
        }
      }
    }
    FillHist("Mmumu_SR2j"+Label, Mmumu, weight, 0., 200., 200);
    FillHist("M2l2j_SR2j"+Label, M2l2j, weight, 0., 1000., 100);
    FillHist("M3lv_SR2j"+Label, M3lv, weight, 0., 1000., 100);
    if(Mmumu<80){
      FillHist("M2l2j_lt80_SR2j"+Label, M2l2j, weight, 0., 1000., 100);
      FillHist("M3lv_lt80_SR2j"+Label, M3lv, weight, 0., 1000., 100);
    }
  }
}


void Aug2017_TriLepSR::DoSystRun(TString Cycle, TString Mode, std::vector<snu::KElectron> EleColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight){

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
    else if(Mode.Contains("ElID"))    SystKindLabel="_ElID";
    else if(Mode.Contains("MuID"))    SystKindLabel="_MuID";
    else if(Mode.Contains("ElEn"))    SystKindLabel="_ElEn";
    else if(Mode.Contains("MuEn"))    SystKindLabel="_MuEn";
    else if(Mode.Contains("FR"))      SystKindLabel="_FR";
    else if(Mode.Contains("BTag_L"))  SystKindLabel="_BTag_L";
    else if(Mode.Contains("BTag_BC")) SystKindLabel="_BTag_BC";
    else if(Mode.Contains("Trig"))    SystKindLabel="_Trig";
    else if(Mode.Contains("Xsec"))    SystKindLabel="_Xsec";
    else if(Mode.Contains("Conv"))    SystKindLabel="_Conv";
  }
  if     (Mode.Contains("EMuMu"))   ChannelLabel="EMuMu";
  else if(Mode.Contains("TriMu"))   ChannelLabel="TriMu";


  if     (Cycle=="SRYield"){
    CheckSRYield(MuColl, MuLColl, EleColl, EleLColl, JetColl,  BJetColl, MET, METx, METy, weight, SystDirLabel+SystKindLabel, ChannelLabel);
  }//End of Closure Cycle
  else if(Cycle=="SRDist"){
    CheckSRDist(MuColl, MuLColl, EleColl, EleLColl, JetColl, BJetColl, MET, METx, METy, weight, SystDirLabel+SystKindLabel, ChannelLabel);
  }


  return;
}



float Aug2017_TriLepSR::ConeCorrectedPT(snu::KMuon Mu, float TightIsoCut){

  //float PTCorr=Mu.Pt();
  //float PTCorr=Mu.Pt()*(1.+RochIso(Mu,"0.4"));
  float PTCorr=Mu.Pt()*(1.+max((float) 0.,RochIso(Mu,"0.4")-TightIsoCut));
  return PTCorr;

}



float Aug2017_TriLepSR::ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut){

  //float PTCorr=Ele.Pt()*(1+Ele.PFRelIso(0.3));
  float PTCorr=Ele.Pt()*(1.+max(0.,Ele.PFRelIso(0.3)-TightIsoCut));
  return PTCorr;

}



void Aug2017_TriLepSR::FillCutFlow(TString cut, float weight){
  
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



void Aug2017_TriLepSR::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Aug2017_TriLepSR::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //AnalyzerCore::MakeHistograms("Basic_Nj_orig", 10, 0., 10.);
  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Aug2017_TriLepSRCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Aug2017_TriLepSR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
