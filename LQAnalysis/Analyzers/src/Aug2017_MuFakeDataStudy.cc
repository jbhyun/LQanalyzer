/***************************************************************************
 * @Project: Aug2017_MuFakeDataStudy 
 * @Package: LQanalyzer, ROOT-based analysis framework for Korea SNU / Main package author: John Almond
 *
 * @Code Author: Jihwan Bhyun 
 *
 ***************************************************************************/

/// Local includes
 #include "Aug2017_MuFakeDataStudy.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Aug2017_MuFakeDataStudy);

 Aug2017_MuFakeDataStudy::Aug2017_MuFakeDataStudy() : AnalyzerCore(), out_muons(0) {

   SetLogName("Aug2017_MuFakeDataStudy");
   Message("In Aug2017_MuFakeDataStudy constructor", INFO);
   InitialiseAnalysis();
 }


 void Aug2017_MuFakeDataStudy::InitialiseAnalysis() throw( LQError ) {
   
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

void Aug2017_MuFakeDataStudy::ExecuteEvents()throw( LQError ){

////Basic Infos///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Apply the gen weight 
   float geneff_weight=1., gennorm_weight=1., k_factor_weight=1.;

   if(!isData){ weight*=MCweight;
              geneff_weight   = GenFilterEfficiency(k_sample_name);
              gennorm_weight  = SignalNorm(k_sample_name, 200.);
              k_factor_weight = GetKFactor();
              weight *= geneff_weight*gennorm_weight*k_factor_weight;
   }
   FillHist("GenWeight", MCweight, 1, -2., 2., 4);

   /// Acts on data to remove bad reconstructed event 
   //if(isData&& (! eventbase->GetEvent().LumiMask(lumimask))) return;

   m_logger<<DEBUG<<"RunNumber/Event Number = "<<eventbase->GetEvent().RunNumber()<<" : "<<eventbase->GetEvent().EventNumber()<<LQLogger::endmsg;
   m_logger<<DEBUG<<"isData = "<<isData<<LQLogger::endmsg;


   //Numbet of Vertex NoPUreweight
   FillHist("Nvtx_nocut_NoRW", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);

   //Pileup Reweight
   float pileup_reweight=1., pileup_reweight_systdown=1., pileup_reweight_systup=1.;
   if(!k_isdata){ pileup_reweight          = eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central); 
                  pileup_reweight_systup   = pileup_reweight>0.? eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::up)  /pileup_reweight :0. ;
                  pileup_reweight_systdown = pileup_reweight>0.? eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::down)/pileup_reweight :0. ;
   }

   //Numbet of Vertex PUreweight
   FillHist("Nvtx_nocut_PURW", eventbase->GetEvent().nVertices(), weight*pileup_reweight, 0., 50., 50);


   //Total Event(MCweight+PUrw)////////////////
   FillCutFlow("NoCut", weight*pileup_reweight);


  
   bool IPEffScan=false, CheckIDVarDist=false;
   bool NormCheck=false, FRMeasure=false, SiglWP=false, FRScan=false, IDValidation=false;
   bool Closure=false, SSDilepCR=false, SystRun=false;
   bool DoubleMuon=false, ElectronMuon=false;
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("IPEffScan"))      IPEffScan      = true;
     else if(k_flags.at(i).Contains("CheckIDVarDist")) CheckIDVarDist = true;
     else if(k_flags.at(i).Contains("NormCheck"))      NormCheck      = true;
     else if(k_flags.at(i).Contains("FRMeasure"))      FRMeasure      = true;
     else if(k_flags.at(i).Contains("SiglWP"))         SiglWP         = true;
     else if(k_flags.at(i).Contains("FRScan"))         FRScan         = true;
     else if(k_flags.at(i).Contains("IDValidation"))   IDValidation   = true;
     else if(k_flags.at(i).Contains("Closure"))        Closure        = true;
     else if(k_flags.at(i).Contains("SystRun"))        SystRun        = true;
     else if(k_flags.at(i).Contains("DoubleMuon"))     DoubleMuon     = true;
     else if(k_flags.at(i).Contains("ElectronMuon"))   ElectronMuon   = true;
     else if(k_flags.at(i).Contains("SSDilepCR"))      SSDilepCR      = true;
   }

    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //ListTriggersAvailable();

   //Normalisation Lumi Setting
   bool TrigMuIso8=false, TrigMuIso17=false, TrigMu8=false, TrigMu17=false;
   bool Pass_Trigger=false, Pass_TriggerBG=false, Pass_TriggerH=false;
   float trigger_ps_weight=1., trigger_period_weight=1.;
   float LumiBG=27.257618, LumiH=8.605696, LumiBH=35.863314;

   if(IPEffScan || CheckIDVarDist || IDValidation || Closure || SSDilepCR){
     if(DoubleMuon){
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
     if(ElectronMuon){
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
   }
   if(NormCheck || FRMeasure){
     if(PassTrigger("HLT_Mu17_TrkIsoVVL_v")) TrigMuIso17 =true;
     if(PassTrigger("HLT_Mu8_TrkIsoVVL_v"))  TrigMuIso8  =true;
     if(PassTrigger("HLT_Mu17_v"))           TrigMu17    =true;
     if(PassTrigger("HLT_Mu8_v"))            TrigMu8     =true;

     if(TrigMu17 || TrigMu8 || TrigMuIso17 || TrigMuIso8) Pass_Trigger=true;
   }
   FillHist("TriggerPSWeight", trigger_ps_weight, 1., 0., 1., 100);
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
   


///////Objects in Analysis/////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

   //Primary Object Collection
   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
  
   /**PreSelCut***********************************************************************************************/
   //Intended for Code speed boosting up.
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, true);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
   if(IPEffScan || CheckIDVarDist || IDValidation){ if(!(muonPreColl.size()>=2)) return; }
   else if(Closure){ if(!(muonPreColl.size()>=3)) return; }
   else if(SSDilepCR){ bool Pass=false;
     int Nlep=muonPreColl.size()+electronPreColl.size();
     if(Nlep>=3) Pass=true;
     else if(Nlep==2){ int SumCharge=0; for(int i=0; i<(int) muonPreColl.size(); i++){ SumCharge+=muonPreColl.at(i).Charge();}
                                        for(int i=0; i<(int) electronPreColl.size(); i++){ SumCharge+=electronPreColl.at(i).Charge(); }  
                       if(SumCharge!=0) Pass=true; }
     if(!Pass) return;
   }
   /**********************************************************************************************************/

   //For Fake Study
   //Muon ID's to Test
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.2);                 eventbase->GetMuonSel()->SetdxySigMax(4.);
     eventbase->GetMuonSel()->SetBSdz(0.1);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.6);
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);

//     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
//     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
//     eventbase->GetMuonSel()->SetBSdxy(0.5);         
//     eventbase->GetMuonSel()->SetBSdz(0.1);
//     eventbase->GetMuonSel()->SetChiNdof(100.);
//     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
//   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);

     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(4.);
     eventbase->GetMuonSel()->SetBSdz(0.05);
     eventbase->GetMuonSel()->SetChiNdof(4.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.2);
   std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,true);


   std::vector<snu::KMuon> muonColl;
     if(!k_running_nonprompt){ muonColl=muonTightColl; }else{ muonColl=muonLooseColl; }


     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_HctoWA_FAKELOOSE);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
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
   std::vector<snu::KElectron> electronColl;
     if(!k_running_nonprompt){ electronColl=electronTightColl; }else{ electronColl=electronLooseColl; }

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     bool LeptonVeto=false;
   std::vector<snu::KJet> jetNoVetoColl; eventbase->GetJetSel()->Selection(jetNoVetoColl, LeptonVeto, muonLooseColl, electronLooseColl);
   std::vector<snu::KJet> bjetNoVetoColl = SelBJets(jetNoVetoColl, "Medium");

   std::vector<snu::KJet> jetColl  = SkimJetColl(jetNoVetoColl,  electronLooseColl, muonLooseColl, "EleMuVeto");
   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");


   float  met    = eventbase->GetEvent().MET(snu::KEvent::pfmet);
   float  metphi = eventbase->GetEvent().METPhi();
   double met_x  = eventbase->GetEvent().PFMETx();
   double met_y  = eventbase->GetEvent().PFMETy();
   double Pzv,Pzv1, Pzv2;
   snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));

   int Nvtx=eventbase->GetEvent().nVertices();


   /*****************************************************
   **Scale Factors
   *****************************************************/
   float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1., btag_sf=1.;
   float trigger_sf=1.;
   float fake_weight=1.; bool EventCand=false;
   float nvtx_reweight=1.; 

   /*This part is for boosting up speed.. SF part takes rather longer time than expected*/
   if(!SystRun){
     if     (NormCheck   ){ EventCand=false; }
     else if(IDValidation){ if(muonLooseColl.size()>=2) EventCand=true; }
     else if(Closure     ){ if(muonLooseColl.size()>=3) EventCand=true; }
     else if(SSDilepCR   ){ if(muonLooseColl.size()>=2) EventCand=true; }

     if(!isData){
       if(EventCand && (NormCheck || IDValidation || Closure || SSDilepCR)){
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl);
      
         trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl);

         //if(UnPreTrig) nvtx_reweight = mcdata_correction->UserPileupWeight(eventbase->GetEvent(), jetColl.size());
         //if(PreTrig)   nvtx_reweight = GetPreTrigPURW(Nvtx);
         if(Closure){
           btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);
         }
         if(DoubleMuon){
           float trigger_sf1 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v");
           float trigger_sf2 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v");
           trigger_sf    = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
         }
         if(ElectronMuon){
           float trigger_sf1 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
           float trigger_sf2 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
           trigger_sf    = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
         }
       }
     }
     else if(k_running_nonprompt && EventCand && Closure ){

       //fake_weight = GetFakeWeight(muonLooseColl, "POGTIsop6IPp2p1NoChi" , "POGTIsop20IPp01p05sig4Chi4", bjetNoVetoColl, "TrkIsoVVLConeSUSY");
       //fake_weight = GetFakeWeight(muonLooseColl, "Test_POGLIsop6IPp5p1", "Test_POGTIsop20IPp01p05sig4Chi4", "TrkIsoVVLConeSUSY");
     }
   }
   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight*fake_weight*btag_sf;
   /***************************************************************************************************/



   //=============================================================================//
   // Main Analysis Code     -----------------------------------------------------//
   //=============================================================================//


   if(IPEffScan){

     ScanIPEfficiency(muonPreColl, electronLooseColl, 0.15, weight, "", "");

     if(!isData){
       std::vector<snu::KMuon> FakeColl = SkimLepColl(muonPreColl, truthColl, "HFake");
       std::vector<snu::KMuon> FakeTypem1Coll;
       std::vector<snu::KMuon> FakeTypem23Coll;
        for(int i=0; i<muonPreColl.size(); i++){
          int LepType=GetLeptonType(muonPreColl.at(i),truthColl);
          if(LepType==-1) FakeTypem1Coll.push_back(muonPreColl.at(i));
          else if(LepType==-2 || LepType==-3) FakeTypem23Coll.push_back(muonPreColl.at(i));
        }
       ScanIPEfficiency(FakeColl, electronLooseColl, 0.15, weight, "_HFake", "");
       ScanIPEfficiency(FakeTypem1Coll, electronLooseColl, 0.15, weight, "_Typem1", "");
       ScanIPEfficiency(FakeTypem23Coll, electronLooseColl, 0.15, weight, "_Typem23", "");
     }
   }
   if(CheckIDVarDist){

     CheckIDVarDistribution(muonPreColl, electronLooseColl, weight, "", "");

     if(!isData){
       std::vector<snu::KMuon> FakeColl = SkimLepColl(muonPreColl, truthColl, "HFake");
       std::vector<snu::KMuon> FakeTypem1Coll;
       std::vector<snu::KMuon> FakeTypem23Coll;
        for(int i=0; i<muonPreColl.size(); i++){
          int LepType=GetLeptonType(muonPreColl.at(i),truthColl);
          if(LepType==-1) FakeTypem1Coll.push_back(muonPreColl.at(i));
          else if(LepType==-2 || LepType==-3) FakeTypem23Coll.push_back(muonPreColl.at(i));
        }
       CheckIDVarDistribution(FakeColl, electronLooseColl, weight, "_HFake", "");
       CheckIDVarDistribution(FakeTypem1Coll, electronLooseColl, weight, "_Typem1", "");
       CheckIDVarDistribution(FakeTypem23Coll, electronLooseColl, weight, "_Typem23", "");
     }

   }
   if(NormCheck){

     float trigger_ps_weight_MuIso17=1., trigger_ps_weight_MuIso8=1., nvtx_reweight_Iso17=1., nvtx_reweight_Iso8=1.;
     float normsf_MuIso17=1., normsf_MuIso8=1.;
     if(!isData){
       trigger_ps_weight_MuIso17 = WeightByTrigger("HLT_Mu17_TrkIsoVVL_v", TargetLumi);
       trigger_ps_weight_MuIso8  = WeightByTrigger("HLT_Mu8_TrkIsoVVL_v" , TargetLumi);
       nvtx_reweight_Iso17       = NvtxWeight(Nvtx, "HLT_Mu17_TrkIsoVVL_v"); 
       nvtx_reweight_Iso8        = NvtxWeight(Nvtx, "HLT_Mu8_TrkIsoVVL_v"); 
       trk_weight_mu             = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
       normsf_MuIso17 = 0.887126*trigger_ps_weight_MuIso17;
       normsf_MuIso8  = 0.884726*trigger_ps_weight_MuIso8;
     }

     if(!SystRun){
       if(TrigMuIso17){
         CheckNormCR(muonColl, muonLooseColl, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight*trigger_ps_weight_MuIso17*trk_weight_mu, "_MuIso17", "TrigMu17");
       }
       if(TrigMuIso8){
         CheckNormCR(muonColl, muonLooseColl, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight*trigger_ps_weight_MuIso8*trk_weight_mu, "_MuIso8" , "TrigMu8");
       }
     }
     else{
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
       float trk_weight_mu_MuEnup=1.;
       float trk_weight_mu_MuEndown=1.;
       if(!isData){
         trk_weight_mu           = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
         trk_weight_mu_MuEnup    = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnUpColl);
         trk_weight_mu_MuEndown  = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnDownColl);
       }
       float systweightIso17_central  =weight*nvtx_reweight_Iso17*normsf_MuIso17*trk_weight_mu       ;
       float systweightIso17_Nvtxup   =weight*                    normsf_MuIso17*trk_weight_mu       ;
       float systweightIso17_MuEnup   =weight*nvtx_reweight_Iso17*normsf_MuIso17*trk_weight_mu_MuEnup;
       float systweightIso17_PUup     =weight*nvtx_reweight_Iso17*normsf_MuIso17*trk_weight_mu         *pileup_reweight_systup;
       float systweightIso17_Nvtxdown =weight*                    normsf_MuIso17*trk_weight_mu       ;
       float systweightIso17_MuEndown =weight*nvtx_reweight_Iso17*normsf_MuIso17*trk_weight_mu_MuEndown;
       float systweightIso17_PUdown   =weight*nvtx_reweight_Iso17*normsf_MuIso17*trk_weight_mu         *pileup_reweight_systdown;

       float systweightIso8_central  =weight*nvtx_reweight_Iso8*normsf_MuIso8*trk_weight_mu       ;
       float systweightIso8_Nvtxup   =weight*                   normsf_MuIso8*trk_weight_mu       ;
       float systweightIso8_MuEnup   =weight*nvtx_reweight_Iso8*normsf_MuIso8*trk_weight_mu_MuEnup;
       float systweightIso8_PUup     =weight*nvtx_reweight_Iso8*normsf_MuIso8*trk_weight_mu       *pileup_reweight_systup;
       float systweightIso8_Nvtxdown =weight*                   normsf_MuIso8*trk_weight_mu       ;
       float systweightIso8_MuEndown =weight*nvtx_reweight_Iso8*normsf_MuIso8*trk_weight_mu_MuEndown;
       float systweightIso8_PUdown   =weight*nvtx_reweight_Iso8*normsf_MuIso8*trk_weight_mu       *pileup_reweight_systdown;


       if(TrigMuIso17){

         DoSystRun("NormCheck", "Mu17",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweightIso17_central);
         if(!isData){
           DoSystRun("NormCheck", "Mu17SystUpPU",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweightIso17_PUup);
           DoSystRun("NormCheck", "Mu17SystUpNvtx",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweightIso17_Nvtxup);
           DoSystRun("NormCheck", "Mu17SystUpMuEn",
                   electronColl, electronLooseColl, MuTMuEnUpColl, muonLooseColl, jetColl, bjetColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi),
                   systweightIso17_MuEnup);
           DoSystRun("NormCheck", "Mu17SystUpElEn",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi),
                   systweightIso17_central);
           DoSystRun("NormCheck", "Mu17SystUpUncl",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi),
                   systweightIso17_central);
           DoSystRun("NormCheck", "Mu17SystUpJES",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESUpColl, bjetColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi),
                   systweightIso17_central);
           DoSystRun("NormCheck", "Mu17SystUpJER",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERUpColl, bjetColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi),
                   systweightIso17_central);

           DoSystRun("NormCheck", "Mu17SystDownPU",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweightIso17_PUdown);
           DoSystRun("NormCheck", "Mu17SystDownNvtx",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweightIso17_Nvtxdown);
           DoSystRun("NormCheck", "Mu17SystDownMuEn",
                   electronColl, electronLooseColl, MuTMuEnDownColl, muonLooseColl, jetColl, bjetColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi),
                   systweightIso17_MuEndown);
           DoSystRun("NormCheck", "Mu17SystDownElEn",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi),
                   systweightIso17_central);
           DoSystRun("NormCheck", "Mu17SystDownUncl",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi),
                   systweightIso17_central);
           DoSystRun("NormCheck", "Mu17SystDownJES",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESDownColl, bjetColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi),
                   systweightIso17_central);
           DoSystRun("NormCheck", "Mu17SystDownJER",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERDownColl, bjetColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi),
                   systweightIso17_central);
         }
       }//End of Mu17
       if(TrigMuIso8){
         DoSystRun("NormCheck", "Mu8",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweightIso8_central);
         if(!isData){
           DoSystRun("NormCheck", "Mu8SystUpPU",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweightIso8_PUup);
           DoSystRun("NormCheck", "Mu8SystUpNvtx",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweightIso8_Nvtxup);
           DoSystRun("NormCheck", "Mu8SystUpMuEn",
                   electronColl, electronLooseColl, MuTMuEnUpColl, muonLooseColl, jetColl, bjetColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi),
                   systweightIso8_MuEnup);
           DoSystRun("NormCheck", "Mu8SystUpElEn",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi),
                   systweightIso8_central);
           DoSystRun("NormCheck", "Mu8SystUpUncl",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi),
                   systweightIso8_central);
           DoSystRun("NormCheck", "Mu8SystUpJES",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESUpColl, bjetColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi),
                   systweightIso8_central);
           DoSystRun("NormCheck", "Mu8SystUpJER",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERUpColl, bjetColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi),
                   systweightIso8_central);

           DoSystRun("NormCheck", "Mu8SystDownPU",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweightIso8_PUdown);
           DoSystRun("NormCheck", "Mu8SystDownNvtx",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweightIso8_Nvtxdown);
           DoSystRun("NormCheck", "Mu8SystDownMuEn",
                   electronColl, electronLooseColl, MuTMuEnDownColl, muonLooseColl, jetColl, bjetColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi),
                   systweightIso8_MuEndown);
           DoSystRun("NormCheck", "Mu8SystDownElEn",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi),
                   systweightIso8_central);
           DoSystRun("NormCheck", "Mu8SystDownUncl",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi),
                   systweightIso8_central);
           DoSystRun("NormCheck", "Mu8SystDownJES",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESDownColl, bjetColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi),
                   systweightIso8_central);
           DoSystRun("NormCheck", "Mu8SystDownJER",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERDownColl, bjetColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi),
                   systweightIso8_central);
         }//End of Not isData
       }//End of Mu8

     }//End of SystRun
   }//End of NormCheck
   if(FRMeasure){
     if(SiglWP){
       std::vector<snu::KMuon> muonPOGTIsop6IPp01p05sig4Chi4Coll, muonPOGTIsop6IPp2p1sig4Coll, muonPOGTIsop6IPp2p1Coll;
       for(int i=0; i<muonPreColl.size(); i++){
         if(PassIDCriteria(muonPreColl.at(i),"POGTIsop6IPp01p05sig4Chi4","Roch")) muonPOGTIsop6IPp01p05sig4Chi4Coll.push_back(muonPreColl.at(i));
         if(PassIDCriteria(muonPreColl.at(i),"POGTIsop6IPp2p1sig4","Roch"))       muonPOGTIsop6IPp2p1sig4Coll.push_back(muonPreColl.at(i));
         if(PassIDCriteria(muonPreColl.at(i),"POGTIsop6IPp2p1","Roch"))           muonPOGTIsop6IPp2p1Coll.push_back(muonPreColl.at(i));
       }

       float trigger_ps_weight_MuIso17=1., trigger_ps_weight_MuIso8=1., nvtx_reweight_Iso17=1., nvtx_reweight_Iso8=1.;
       float trigger_ps_weight_Mu17=1., trigger_ps_weight_Mu8=1., nvtx_reweight_17=1., nvtx_reweight_8=1.;
       float normsf_MuIso8=1., normsf_MuIso17=1.;
       if(!isData){
         trigger_ps_weight_MuIso17 = WeightByTrigger("HLT_Mu17_TrkIsoVVL_v", TargetLumi);
         trigger_ps_weight_MuIso8  = WeightByTrigger("HLT_Mu8_TrkIsoVVL_v" , TargetLumi);
         trigger_ps_weight_Mu17    = WeightByTrigger("HLT_Mu17_v", TargetLumi);
         trigger_ps_weight_Mu8     = WeightByTrigger("HLT_Mu8_v" , TargetLumi);
         nvtx_reweight_Iso17       = NvtxWeight(Nvtx, "HLT_Mu17_TrkIsoVVL_v"); 
         nvtx_reweight_Iso8        = NvtxWeight(Nvtx, "HLT_Mu8_TrkIsoVVL_v"); 
         nvtx_reweight_17          = NvtxWeight(Nvtx, "HLT_Mu17_v"); 
         nvtx_reweight_8           = NvtxWeight(Nvtx, "HLT_Mu8_v"); 
         normsf_MuIso17            = 0.887126*trigger_ps_weight_MuIso17;
         normsf_MuIso8             = 0.884726*trigger_ps_weight_MuIso8;
       }

       MeasureFakeRate(muonPOGTIsop6IPp01p05sig4Chi4Coll,
         electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
         weight*normsf_MuIso8*nvtx_reweight_8,
         "POGTIsop20IPp01p05sig4Chi4", "_POGTIsop20IPp01p05sig4Chi4LIsop6", "TrigMu8");
       MeasureFakeRate(muonPOGTIsop6IPp01p05sig4Chi4Coll,
         electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
         weight*normsf_MuIso17*nvtx_reweight_17,
         "POGTIsop20IPp01p05sig4Chi4", "_POGTIsop20IPp01p05sig4Chi4LIsop6", "TrigMu17");

       MeasureFakeRate(muonPOGTIsop6IPp2p1sig4Coll,
         electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
         weight*normsf_MuIso8*nvtx_reweight_8,
         "POGTIsop20IPp01p05sig4Chi4", "_POGTIsop20IPp01p05sig4Chi4POGTIsop6IPp2p1sig4", "TrigMu8");
       MeasureFakeRate(muonPOGTIsop6IPp2p1sig4Coll,
         electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
         weight*normsf_MuIso17*nvtx_reweight_17,
         "POGTIsop20IPp01p05sig4Chi4", "_POGTIsop20IPp01p05sig4Chi4POGTIsop6IPp2p1sig4", "TrigMu17");

       MeasureFakeRate(muonPOGTIsop6IPp2p1Coll,
         electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
         weight*normsf_MuIso8*nvtx_reweight_8,
         "POGTIsop20IPp01p05sig4Chi4", "_POGTIsop20IPp01p05sig4Chi4POGTIsop6IPp2p1", "TrigMu8");
       MeasureFakeRate(muonPOGTIsop6IPp2p1Coll,
         electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
         weight*normsf_MuIso17*nvtx_reweight_17,
         "POGTIsop20IPp01p05sig4Chi4", "_POGTIsop20IPp01p05sig4Chi4POGTIsop6IPp2p1", "TrigMu17");

     }
     if(FRScan){
     }
   }
   if(IDValidation){

     ValidateID(muonTightColl, weight, "");
   }
   if(SSDilepCR){

     float fake_weight1=1., fake_weight2=1.;
     if(k_running_nonprompt){
       fake_weight2 = GetFakeWeight_Data(muonLooseColl, electronLooseColl, "POGTIsop6IPp2p1sig4" , "POGTIsop20IPp01p05sig4Chi4", "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", "TrkIsoVVLConeSUSY");
       //fake_weight2 = GetFakeWeight(muonLooseColl, "POGTIsop6IPp2p1sig4", "POGTIsop20IPp01p05sig4Chi4", "TrkIsoVVLConeSUSY");
       //fake_weight2 = GetFakeWeight(muonLooseColl, "POGLIsop4IPp5p1Chi100", "POGTIsop20IPp01p05sig4Chi4", "TrkIsoVVLConeSUSY");
     }
     if(!isData && (k_sample_name.Contains("TT_powheg") || k_sample_name.Contains("DYJets") || k_sample_name.Contains("WJets"))){
       int NLepCand=0;
       for(int i=0; i<(int) muonColl.size(); i++){
         int LepType=GetLeptonType(muonColl.at(i), truthColl);
         if(LepType>0 || LepType<-4) NLepCand++;
       }
       for(int i=0; i<(int) electronColl.size(); i++){
         int LepType=GetLeptonType(electronColl.at(i), truthColl);
         if(LepType>0 || LepType<-4) NLepCand++;
       }
       if(NLepCand!=2) return;
     }
   
     //CheckSSDilepCRs(muonColl, muonLooseColl, electronColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight*fake_weight2, "_POGLIsop4IPp5p1Chi100_IsoFilter", "");
     CheckSSDilepCRs(muonColl, muonLooseColl, electronColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight*fake_weight2, "_POGTIsop6IPp2p1sig4_IsoFilter", "");
   }
   if(Closure){

     if(!SystRun){
       std::vector<snu::KMuon> muonPOGLIsop4IPp5p1Chi100Coll,
                               muonPOGTIsop4IPp2p1NoChiColl     , muonPOGTIsop6IPp2p1Coll,
                               muonPOGTIsop4IPp2p1sig4NoChiColl , muonPOGTIsop6IPp2p1sig4Coll, 
                               muonPOGTIsop4IPp01p05sig4Chi4Coll, muonPOGTIsop6IPp01p05sig4Chi4Coll;


       for(int i=0; i<(int) muonPreColl.size(); i++){
         //Round2
         if(PassIDCriteria(muonPreColl.at(i),"POGLIsop4IPp5p1Chi100","Roch"))      muonPOGLIsop4IPp5p1Chi100Coll.push_back(muonPreColl.at(i));
         if(PassIDCriteria(muonPreColl.at(i),"POGTIsop4IPp2p1NoChi","Roch"))       muonPOGTIsop4IPp2p1NoChiColl.push_back(muonPreColl.at(i));
         if(PassIDCriteria(muonPreColl.at(i),"POGTIsop4IPp2p1sig4NoChi","Roch"))   muonPOGTIsop4IPp2p1sig4NoChiColl.push_back(muonPreColl.at(i));
         if(PassIDCriteria(muonPreColl.at(i),"POGTIsop6IPp2p1","Roch"))            muonPOGTIsop6IPp2p1Coll.push_back(muonPreColl.at(i));
         if(PassIDCriteria(muonPreColl.at(i),"POGTIsop6IPp2p1sig4","Roch"))        muonPOGTIsop6IPp2p1sig4Coll.push_back(muonPreColl.at(i));
         if(PassIDCriteria(muonPreColl.at(i),"POGTIsop6IPp01p05sig4Chi4","Roch"))  muonPOGTIsop6IPp01p05sig4Chi4Coll.push_back(muonPreColl.at(i));
       }

       float fake_weight1=1., fake_weight2=1., fake_weight3=1., fake_weight4=1., fake_weight5=1., fake_weight6=1.; 
       float fake_weight11=1., fake_weight12=1., fake_weight13=1., fake_weight14=1., fake_weight15=1.; 
       if(k_running_nonprompt){
        fake_weight1 = GetFakeWeight(muonPOGLIsop4IPp5p1Chi100Coll    , "POGLIsop4IPp5p1Chi100"    , "POGTIsop20IPp01p05sig4Chi4", "NoFilterConeSUSY");
        fake_weight2 = GetFakeWeight(muonPOGTIsop4IPp2p1NoChiColl     , "POGTIsop4IPp2p1NoChi"     , "POGTIsop20IPp01p05sig4Chi4", "NoFilterConeSUSY");
        fake_weight3 = GetFakeWeight(muonPOGTIsop4IPp2p1sig4NoChiColl , "POGTIsop4IPp2p1sig4NoChi" , "POGTIsop20IPp01p05sig4Chi4", "NoFilterConeSUSY");
        fake_weight4 = GetFakeWeight(muonPOGTIsop6IPp2p1sig4Coll      , "POGTIsop6IPp2p1sig4"      , "POGTIsop20IPp01p05sig4Chi4", "NoFilterConeSUSY");
        fake_weight5 = GetFakeWeight(muonPOGTIsop6IPp01p05sig4Chi4Coll, "POGTIsop6IPp01p05sig4Chi4", "POGTIsop20IPp01p05sig4Chi4", "NoFilterConeSUSY");
        fake_weight6 = GetFakeWeight(muonPOGTIsop6IPp2p1Coll          , "POGTIsop6IPp2p1"          , "POGTIsop20IPp01p05sig4Chi4", "NoFilterConeSUSY");
       }

       if(!k_running_nonprompt){
         CheckTrilepCRs(muonColl, muonPOGLIsop4IPp5p1Chi100Coll    , electronColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight*fake_weight1, "_POGLIsop4IPp5p1Chi100", "");
         CheckTrilepCRs(muonColl, muonPOGTIsop6IPp2p1sig4Coll      , electronColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight*fake_weight4, "_POGTIsop6IPp2p1sig4", "");
         CheckTrilepCRs(muonColl, muonPOGTIsop6IPp01p05sig4Chi4Coll, electronColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight*fake_weight5, "_POGTIsop6IPp01p05sig4Chi4", "");

       }
       else{
         CheckTrilepCRs(muonPOGLIsop4IPp5p1Chi100Coll    , muonPOGLIsop4IPp5p1Chi100Coll    , electronColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight*fake_weight1, "_POGLIsop4IPp5p1Chi100", "");
        CheckTrilepCRs(muonPOGTIsop6IPp2p1sig4Coll      , muonPOGTIsop6IPp2p1sig4Coll      , electronColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight*fake_weight4, "_POGTIsop6IPp2p1sig4", "");
         CheckTrilepCRs(muonPOGTIsop6IPp01p05sig4Chi4Coll, muonPOGTIsop6IPp01p05sig4Chi4Coll, electronColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight*fake_weight5, "_POGTIsop6IPp01p05sig4Chi4", "");

       }
       //CheckTrilepCRs(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
        // weight, "", "");
     }
     else{

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
       float systweight=weight; //Lumi weight+Prescale weight+Period weight+PU weight applied by default;
       float id_weight_ele_ElEnup=1.  , reco_weight_ele_ElEnup=1.  , id_weight_mu_MuEnup=1.  , trk_weight_mu_MuEnup=1.;
       float id_weight_ele_ElEndown=1., reco_weight_ele_ElEndown=1., id_weight_mu_MuEndown=1., trk_weight_mu_MuEndown=1.;
       float btag_sf_JESup=1.  , btag_sf_JERup=1.  , btag_sf_LTagup=1.  , btag_sf_BCTagup=1.;
       float btag_sf_JESdown=1., btag_sf_JERdown=1., btag_sf_LTagdown=1., btag_sf_BCTagdown=1.;
       float fake_weight_FRup=1., fake_weight_FRdown=1.;
       float trigger_sf=1., trigger_sf_up=1., trigger_sf_down=1.;
       if(!isData){
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl);
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl);
         trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
         btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

         id_weight_ele_ElEnup   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", EleTElEnUpColl);
         reco_weight_ele_ElEnup = mcdata_correction->ElectronRecoScaleFactor(EleTElEnUpColl);
         id_weight_mu_MuEnup    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuTMuEnUpColl);
         trk_weight_mu_MuEnup   = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnUpColl);

         id_weight_ele_ElEndown   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", EleTElEnDownColl);
         reco_weight_ele_ElEndown = mcdata_correction->ElectronRecoScaleFactor(EleTElEnDownColl);
         id_weight_mu_MuEndown    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuTMuEnDownColl);
         trk_weight_mu_MuEndown   = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnDownColl);

         if(ElectronMuon){
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
         fake_weight = GetFakeWeight(muonLooseColl, "POGTIsop6IPp2p1sig4", "POGTIsop20IPp01p05sig4Chi4", bjetNoVetoColl, "TrkIsoVVLConeSUSY");
       }
       float systweight_central=weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight *trigger_sf;

       float systweight_Trigup =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight *trigger_sf_up;
       float systweight_ElEnup =weight*k_factor_weight*id_weight_ele_ElEnup*reco_weight_ele_ElEnup*id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight *trigger_sf;
       float systweight_MuEnup =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu_MuEnup*trk_weight_mu_MuEnup*btag_sf        *fake_weight *trigger_sf;
       float systweight_JESup  =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JESup  *fake_weight *trigger_sf;
       float systweight_JERup  =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JERup  *fake_weight *trigger_sf;
       float systweight_LTagup =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_LTagup *fake_weight *trigger_sf;
       float systweight_BCTagup=weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_BCTagup*fake_weight *trigger_sf;
       float systweight_PUup   =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight*pileup_reweight_systup *trigger_sf;

       float systweight_Trigdown =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *fake_weight *trigger_sf_down;
       float systweight_ElEndown =weight*k_factor_weight*id_weight_ele_ElEndown*reco_weight_ele_ElEndown*id_weight_mu         *trk_weight_mu         *btag_sf          *fake_weight *trigger_sf;
       float systweight_MuEndown =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu_MuEndown*trk_weight_mu_MuEndown*btag_sf          *fake_weight *trigger_sf;
       float systweight_JESdown  =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JESdown  *fake_weight *trigger_sf;
       float systweight_JERdown  =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JERdown  *fake_weight *trigger_sf;
       float systweight_LTagdown =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_LTagdown *fake_weight *trigger_sf;
       float systweight_BCTagdown=weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_BCTagdown*fake_weight *trigger_sf;
       float systweight_PUdown   =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *fake_weight*pileup_reweight_systdown *trigger_sf;



       //----------------------------------------------------------------------------------------------------------------------//


       DoSystRun("Closure", "",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_central);
       if(!isData){
         DoSystRun("Closure", "SystUpTrig",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_Trigup);
         DoSystRun("Closure", "SystUpPU",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_PUup);
         DoSystRun("Closure", "SystUpUncl",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi),
                   systweight_central);
         DoSystRun("Closure", "SystUpJES",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESUpColl, bjetJESUpColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi),
                   systweight_JESup);
         DoSystRun("Closure", "SystUpJER",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERUpColl, bjetJERUpColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi),
                   systweight_JERup);
         DoSystRun("Closure", "SystUpBTag_L",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_LTagup);
         DoSystRun("Closure", "SystUpBTag_BC",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_BCTagup);
         DoSystRun("Closure", "SystUpElEn", EleTElEnUpColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi),
                   systweight_ElEnup);
         DoSystRun("Closure", "SystUpMuEn", electronColl, electronLooseColl, MuTMuEnUpColl, muonLooseColl, jetColl, bjetColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi),
                   systweight_MuEnup);


         DoSystRun("Closure", "SystDownTrig",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_Trigdown);
         DoSystRun("Closure", "SystDownPU",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_PUdown);
         DoSystRun("Closure", "SystDownUncl",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi),
                   systweight_central);
         DoSystRun("Closure", "SystDownJES",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESDownColl, bjetJESDownColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi),
                   systweight_JESdown);
         DoSystRun("Closure", "SystDownJER",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERDownColl, bjetJERDownColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi),
                   systweight_JERdown);
         DoSystRun("Closure", "SystDownBTag_L",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_LTagdown);
         DoSystRun("Closure", "SystDownBTag_BC",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_BCTagdown);
         DoSystRun("Closure", "SystDownElEn", EleTElEnDownColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi),
                   systweight_ElEndown);
         DoSystRun("Closure", "SystDownMuEn", electronColl, electronLooseColl, MuTMuEnDownColl, muonLooseColl, jetColl, bjetColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi),
                   systweight_MuEndown);

       }
       else if( isData && k_running_nonprompt ){

//         DoSystRun("Closure", "", electronLooseColl, electronLooseColl, muonLooseColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), systweight);

       }

     }//End of SystRun

  }//End of Closure


return;
}// End of execute event loop
  

void Aug2017_MuFakeDataStudy::ScanIPEfficiency(std::vector<snu::KMuon> muonColl, std::vector<snu::KElectron> electronLooseColl, float IsoCut, float weight, TString Label, TString Option){
  //Input : 
  //Purpose : Check Data IP Efficiency Map
  //         - If loosen iso, tight IP must compromise the loss of fake rejection, i.e. IP must be very tight.


  if(isData){
    if( !(muonColl.size()==2 && electronLooseColl.size()==0) ) return;
    if( !(muonColl.at(0).Charge()!= muonColl.at(1).Charge()) ) return;
    if( !(muonColl.at(0).Pt()>20 && muonColl.at(1).Pt()>10 ) ) return;
    if( !(fabs((muonColl.at(0)+muonColl.at(1)).M()-91.2)<15) ) return;
  }

  std::ostringstream s1; s1<<IsoCut;
  TString Str_IsoWP=s1.str();
  Str_IsoWP.ReplaceAll(".","p");

  for(int i=0; i<muonColl.size(); i++){
    if(!muonColl.at(i).IsTight()) return;
    if(!(muonColl.at(i).RelIso04()<IsoCut)) return;

    float d0=fabs(muonColl.at(i).dXY()), dz=fabs(muonColl.at(i).dZ()), d0Sig=fabs(muonColl.at(i).dXYSig());
    FillHist("Data_d0_iso"+Str_IsoWP+Label, d0, weight, 0., 0.05, 20);
    FillHist("Data_dz_iso"+Str_IsoWP+Label, dz, weight, 0., 0.1, 40);
    FillHist("Data_d0Sig_iso"+Str_IsoWP+Label, d0Sig, weight, 0., 10., 10);

    //IP ROC Plane
    int Nd0WP=20, NdzWP=40;
    for(int it_d0=1; it_d0<=Nd0WP; it_d0++){
      for(int it_dz=1; it_dz<=NdzWP; it_dz++){
        FillHist("NEleSumW_d0dz_iso"+Str_IsoWP+Label, (it_d0-1)*0.0025001, (it_dz-1)*0.0025001, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
        if(d0<it_d0*0.0025 && dz< it_dz*0.0025){
          FillHist("NEleIDSumW_d0dz_iso"+Str_IsoWP+Label, (it_d0-1)*0.0025001, (it_dz-1)*0.0025001, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
        }
      }
    }

    if(muonColl.at(i).Pt()<15){
      FillHist("Data_d0_iso"+Str_IsoWP+"Pt10to15"+Label, d0, weight, 0., 0.05, 20);
      FillHist("Data_dz_iso"+Str_IsoWP+"Pt10to15"+Label, dz, weight, 0., 0.1, 40);
      FillHist("Data_d0Sig_iso"+Str_IsoWP+"Pt10to15"+Label, d0Sig, weight, 0., 10., 10);
  
      //IP ROC Plane
      int Nd0WP=20, NdzWP=40;
      for(int it_d0=1; it_d0<=Nd0WP; it_d0++){
        for(int it_dz=1; it_dz<=NdzWP; it_dz++){
          FillHist("NEleSumW_d0dz_iso"+Str_IsoWP+"Pt10to15"+Label, (it_d0-1)*0.0025001, (it_dz-1)*0.0025001, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
          if(d0<it_d0*0.0025 && dz< it_dz*0.0025){
            FillHist("NEleIDSumW_d0dz_iso"+Str_IsoWP+"Pt10to15"+Label, (it_d0-1)*0.0025001, (it_dz-1)*0.0025001, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
          }
        }
      }
    }//End of Pt10to15
    else if(muonColl.at(i).Pt()<20){
      FillHist("Data_d0_iso"+Str_IsoWP+"Pt15to20"+Label, d0, weight, 0., 0.05, 20);
      FillHist("Data_dz_iso"+Str_IsoWP+"Pt15to20"+Label, dz, weight, 0., 0.1, 40);
      FillHist("Data_d0Sig_iso"+Str_IsoWP+"Pt15to20"+Label, d0Sig, weight, 0., 10., 10);
  
      //IP ROC Plane
      int Nd0WP=20, NdzWP=40;
      for(int it_d0=1; it_d0<=Nd0WP; it_d0++){
        for(int it_dz=1; it_dz<=NdzWP; it_dz++){
          FillHist("NEleSumW_d0dz_iso"+Str_IsoWP+"Pt15to20"+Label, (it_d0-1)*0.0025001, (it_dz-1)*0.0025001, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
          if(d0<it_d0*0.0025 && dz< it_dz*0.0025){
            FillHist("NEleIDSumW_d0dz_iso"+Str_IsoWP+"Pt15to20"+Label, (it_d0-1)*0.0025001, (it_dz-1)*0.0025001, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
          }
        }
      }
    }//End of Pt15to20

  }//End of MuLoop


}



void Aug2017_MuFakeDataStudy::CheckIDVarDistribution(std::vector<snu::KMuon> muonColl, std::vector<snu::KElectron> electronLooseColl, float weight, TString Label, TString Option){
  //Input : 
  //Purpose : Check Data IP Efficiency Map
  //         - If loosen iso, tight IP must compromise the loss of fake rejection, i.e. IP must be very tight.


  if(isData){
    if( !(muonColl.size()==2 && electronLooseColl.size()==0) ) return;
    if( !(muonColl.at(0).Charge()!= muonColl.at(1).Charge()) ) return;
    if( !(muonColl.at(0).Pt()>20 && muonColl.at(1).Pt()>10 ) ) return;
    if( !(fabs((muonColl.at(0)+muonColl.at(1)).M()-91.2)<15) ) return;
  }

  for(int i=0; i<muonColl.size(); i++){

    float PT=muonColl.at(i).Pt(), fEta=fabs(muonColl.at(i).Eta());
    bool PassCut = muonColl.at(i).IsPF() && muonColl.at(i).IsGlobal()
                  && muonColl.at(i).validPixHits()>0 && muonColl.at(i).ActiveLayer()>5
                  && muonColl.at(i).validStations()>1;
                  //&& muonColl.at(i).validHits()>0 && muonColl.at(i).validStations()>1

    FillHist("IsGlobal"+Label, 0., weight, 0., 2., 2);
    if(muonColl.at(i).IsGlobal()) FillHist("IsGlobal"+Label, 1., weight, 0., 2., 2); 
    FillHist("NValidHits"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
    FillHist("NValidPixHits"+Label, muonColl.at(i).validPixHits(), weight, 0., 20., 20);
    FillHist("NValidStations"+Label, muonColl.at(i).validStations(), weight, 0., 20., 20);
    FillHist("NActiveLayer"+Label, muonColl.at(i).ActiveLayer(), weight, 0., 20., 20);
    FillHist("GlobalChi2"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 20., 20);
    FillHist("D0"+Label, muonColl.at(i).dXY(), weight, 0., 0.05, 50);
    FillHist("DZ"+Label, muonColl.at(i).dZ(), weight, 0., 0.2, 200); 

    if(fEta<1.2){
      if(PT<20){
        FillHist("IsGlobal_MB10to20"+Label, 0., weight, 0., 2., 2);
        if(muonColl.at(i).IsGlobal()) FillHist("IsGlobal_MB10to20"+Label, 1., weight, 0., 2., 2); 
        FillHist("NValidHits_MB10to20"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
        FillHist("NValidPixHits_MB10to20"+Label, muonColl.at(i).validPixHits(), weight, 0., 20., 20);
        FillHist("NValidStations_MB10to20"+Label, muonColl.at(i).validStations(), weight, 0., 20., 20);
        FillHist("NActiveLayer_MB10to20"+Label, muonColl.at(i).ActiveLayer(), weight, 0., 20., 20);
        FillHist("GlobalChi2_MB10to20"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 20., 20);
        FillHist("D0_MB10to20"+Label, muonColl.at(i).dXY(), weight, 0., 0.05, 50);
        FillHist("DZ_MB10to20"+Label, muonColl.at(i).dZ(), weight, 0., 0.2, 200); 
      }
      else if(PT<50){
        FillHist("IsGlobal_MB20to50"+Label, 0., weight, 0., 2., 2);
        if(muonColl.at(i).IsGlobal()) FillHist("IsGlobal_MB20to50"+Label, 1., weight, 0., 2., 2); 
        FillHist("NValidHits_MB20to50"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
        FillHist("NValidPixHits_MB20to50"+Label, muonColl.at(i).validPixHits(), weight, 0., 20., 20);
        FillHist("NValidStations_MB20to50"+Label, muonColl.at(i).validStations(), weight, 0., 20., 20);
        FillHist("NActiveLayer_MB20to50"+Label, muonColl.at(i).ActiveLayer(), weight, 0., 20., 20);
        FillHist("GlobalChi2_MB20to50"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 20., 20);
        FillHist("D0_MB20to50"+Label, muonColl.at(i).dXY(), weight, 0., 0.05, 50);
        FillHist("DZ_MB20to50"+Label, muonColl.at(i).dZ(), weight, 0., 0.2, 200); 
      }
      else if(PT<100){
        FillHist("IsGlobal_MB50to100"+Label, 0., weight, 0., 2., 2);
        if(muonColl.at(i).IsGlobal()) FillHist("IsGlobal_MB50to100"+Label, 1., weight, 0., 2., 2); 
        FillHist("NValidHits_MB50to100"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
        FillHist("NValidPixHits_MB50to100"+Label, muonColl.at(i).validPixHits(), weight, 0., 20., 20);
        FillHist("NValidStations_MB50to100"+Label, muonColl.at(i).validStations(), weight, 0., 20., 20);
        FillHist("NActiveLayer_MB50to100"+Label, muonColl.at(i).ActiveLayer(), weight, 0., 20., 20);
        FillHist("GlobalChi2_MB50to100"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 20., 20);
        FillHist("D0_MB50to100"+Label, muonColl.at(i).dXY(), weight, 0., 0.05, 50);
        FillHist("DZ_MB50to100"+Label, muonColl.at(i).dZ(), weight, 0., 0.2, 200); 
      }
    }
    else if(fEta<2.4){
      if(PT<20){
        FillHist("IsGlobal_ME10to20"+Label, 0., weight, 0., 2., 2);
        if(muonColl.at(i).IsGlobal()) FillHist("IsGlobal_ME10to20"+Label, 1., weight, 0., 2., 2); 
        FillHist("NValidHits_ME10to20"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
        FillHist("NValidPixHits_ME10to20"+Label, muonColl.at(i).validPixHits(), weight, 0., 20., 20);
        FillHist("NValidStations_ME10to20"+Label, muonColl.at(i).validStations(), weight, 0., 20., 20);
        FillHist("NActiveLayer_ME10to20"+Label, muonColl.at(i).ActiveLayer(), weight, 0., 20., 20);
        FillHist("GlobalChi2_ME10to20"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 20., 20);
        FillHist("D0_ME10to20"+Label, muonColl.at(i).dXY(), weight, 0., 0.05, 50);
        FillHist("DZ_ME10to20"+Label, muonColl.at(i).dZ(), weight, 0., 0.2, 200); 
      }
      else if(PT<50){
        FillHist("IsGlobal_ME20to50"+Label, 0., weight, 0., 2., 2);
        if(muonColl.at(i).IsGlobal()) FillHist("IsGlobal_ME20to50"+Label, 1., weight, 0., 2., 2); 
        FillHist("NValidHits_ME20to50"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
        FillHist("NValidPixHits_ME20to50"+Label, muonColl.at(i).validPixHits(), weight, 0., 20., 20);
        FillHist("NValidStations_ME20to50"+Label, muonColl.at(i).validStations(), weight, 0., 20., 20);
        FillHist("NActiveLayer_ME20to50"+Label, muonColl.at(i).ActiveLayer(), weight, 0., 20., 20);
        FillHist("GlobalChi2_ME20to50"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 20., 20);
        FillHist("D0_ME20to50"+Label, muonColl.at(i).dXY(), weight, 0., 0.05, 50);
        FillHist("DZ_ME20to50"+Label, muonColl.at(i).dZ(), weight, 0., 0.2, 200); 
      }
      else if(PT<100){
        FillHist("IsGlobal_ME50to100"+Label, 0., weight, 0., 2., 2);
        if(muonColl.at(i).IsGlobal()) FillHist("IsGlobal_ME50to100"+Label, 1., weight, 0., 2., 2); 
        FillHist("NValidHits_ME50to100"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
        FillHist("NValidPixHits_ME50to100"+Label, muonColl.at(i).validPixHits(), weight, 0., 20., 20);
        FillHist("NValidStations_ME50to100"+Label, muonColl.at(i).validStations(), weight, 0., 20., 20);
        FillHist("NActiveLayer_ME50to100"+Label, muonColl.at(i).ActiveLayer(), weight, 0., 20., 20);
        FillHist("GlobalChi2_ME50to100"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 20., 20);
        FillHist("D0_ME50to100"+Label, muonColl.at(i).dXY(), weight, 0., 0.05, 50);
        FillHist("DZ_ME50to100"+Label, muonColl.at(i).dZ(), weight, 0., 0.2, 200); 
      }
    }

    if(PassCut){
      FillHist("NValidHits_PassCut"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
      FillHist("GlobalChi2_PassCut"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 20., 20);
  
      if(fEta<1.2){
        if(PT<20){
          FillHist("NValidHits_MB10to20_PassCut"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
          FillHist("GlobalChi2_MB10to20_PassCut"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 20., 20);
        }
        else if(PT<50){
          FillHist("NValidHits_MB20to50_PassCut"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
          FillHist("GlobalChi2_MB20to50_PassCut"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 20., 20);
        }
        else if(PT<100){
          FillHist("NValidHits_MB50to100_PassCut"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
          FillHist("GlobalChi2_MB50to100_PassCut"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 20., 20);
        }
      }
      else if(fEta<2.4){
        if(PT<20){
          FillHist("NValidHits_ME10to20_PassCut"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
          FillHist("GlobalChi2_ME10to20_PassCut"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 20., 20);
        }
        else if(PT<50){
          FillHist("NValidHits_ME20to50_PassCut"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
          FillHist("GlobalChi2_ME20to50_PassCut"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 20., 20);
        }
        else if(PT<100){
          FillHist("NValidHits_ME50to100_PassCut"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
          FillHist("GlobalChi2_ME50to100_PassCut"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 20., 20);
        }
      }
    }//End of PassCut

  }//End of MuonLoop
}



void Aug2017_MuFakeDataStudy::CheckNormCR(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){

  bool UnPreTrig=false, TrigMu8=false, TrigMu17=false;
  if     (Option.Contains("UnPreTrig")) UnPreTrig=true;
  else if(Option.Contains("TrigMu8"))   TrigMu8  =true;
  else if(Option.Contains("TrigMu17"))  TrigMu17 =true;

  std::vector<snu::KJet> JetVetoColl = SkimJetColl(JetColl, EleLColl, MuLColl, "EleMuVeto");
  int Nvtx=eventbase->GetEvent().nVertices();

  if(MuTColl.size()==1 && MuLColl.size()==1 && EleLColl.size()==0){
    if( UnPreTrig && MuTColl.at(0).Pt()<27.) return;
    if( TrigMu8   && MuTColl.at(0).Pt()<10.) return;
    if( TrigMu17  && MuTColl.at(0).Pt()<20.) return;

    float MTW = sqrt(2)*sqrt(MET*MuTColl.at(0).Pt()-METx*MuTColl.at(0).Px()-METy*MuTColl.at(0).Py());
   
    FillHist("MET"+Label, MET, weight, 0., 200., 40);
    FillHist("MTW"+Label, MTW, weight, 0., 200., 40);
    FillHist("PTmu"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);

    if(MET>50.){
      FillHist("MTW_met50"+Label, MTW, weight, 0., 200., 40);

      //Inclusive Selection CR plot - Incl W+jet test
      if(MTW>70.){
        FillHist("Count_NormCRIncl"+Label, 0., weight, 0., 2., 2);
        FillHist("PTmu_met50mtw70"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
        FillHist("Etamu_met50mtw70"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
        FillHist("Nj_met50mtw70"+Label, JetVetoColl.size(), weight, 0., 10., 10);
        FillHist("MET_met50mtw70"+Label, MET, weight, 0., 200., 40);
      }
    }
    if(JetVetoColl.size()>0 && JetVetoColl.at(0).Pt()>40.){
      FillHist("MET_gt1j40"+Label, MET, weight, 0., 200., 40);
      FillHist("MTW_gt1j40"+Label, MTW, weight, 0., 200., 40);
      FillHist("Nvtx_gt1j40"+Label, Nvtx, weight, 0., 50., 50);

      if(MET>50.){
        FillHist("MTW_gt1j40met50"+Label, MTW, weight, 0., 200., 40);
        FillHist("Nvtx_gt1j40met50"+Label, Nvtx, weight, 0., 50., 50);

        //e+geq1j Selection CR plot - W+geq1j test
        if(MTW>70.){
          FillHist("Count_NormCRgt1j"+Label,   0., weight, 0., 2., 2); 

          FillHist("PTmu_gt1j40met50mtw70"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
          FillHist("Etamu_gt1j40met50mtw70"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
          FillHist("PTj1_gt1j40met50mtw70"+Label, JetVetoColl.at(0).Pt(), weight, 0., 200., 40);
          FillHist("Etaj1_gt1j40met50mtw70"+Label, JetVetoColl.at(0).Eta(), weight, -5., 5., 20);
          FillHist("dRmuj1_gt1j40met50mtw70"+Label, MuTColl.at(0).DeltaR(JetVetoColl.at(0)), weight, 0., 5., 50);
          FillHist("Nj_gt1j40met50mtw70"+Label, JetVetoColl.size(), weight, 0., 10., 10);
          FillHist("MET_gt1j40met50mtw70"+Label, MET, weight, 0., 200., 40);
          FillHist("Nvtx_gt1j40met50mtw70"+Label, Nvtx, weight, 0., 50., 50);
        }
      }
    }
  }//End of 1l

  if(MuTColl.size()==2 && MuLColl.size()==2 && EleLColl.size()==0){
    if( UnPreTrig && !(MuTColl.at(0).Pt()>27. && MuTColl.at(1).Pt()>10.) ) return;
    if( TrigMu8   && !(MuTColl.at(0).Pt()>20. && MuTColl.at(1).Pt()>10.) ) return;
    if( TrigMu17  && !(MuTColl.at(0).Pt()>20. && MuTColl.at(1).Pt()>10.) ) return;
    if( MuTColl.at(0).Charge()==MuTColl.at(1).Charge()  ) return;
    if( fabs((MuTColl.at(0)+MuTColl.at(1)).M()-91.2)>15 ) return;

    FillHist("Mmumu_OS2mu"+Label, (MuTColl.at(0)+MuTColl.at(1)).M(), weight, 60., 120., 60);
    FillHist("Count_NormCRIncl"+Label, 1., weight, 0., 2., 2);
    if(JetVetoColl.size()>0 && JetVetoColl.at(0).Pt()>40.){
      FillHist("Mmumu_OS2mugt1j"+Label, (MuTColl.at(0)+MuTColl.at(1)).M(), weight, 60., 120., 60);
      FillHist("Count_NormCRgt1j"+Label, 1., weight, 0., 2., 2);

      FillHist("Nvtx_OS2mugt1j"+Label, Nvtx, weight, 0., 50., 50);
      FillHist("MET_OS2mugt1j"+Label, MET, weight, 0., 200., 40);
    }
  }//End of 2l


}



void Aug2017_MuFakeDataStudy::MeasureFakeRate(std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString TightID, TString Label, TString Option){

  bool TrigMu8=false, TrigMu17=false;
  TString AddFilterLabel="";
  if(Option.Contains("TrkIsoVVL")) AddFilterLabel="_TrkIsoVVL";
  if(Option.Contains("TrigMu8"))   TrigMu8 =true;
  if(Option.Contains("TrigMu17"))  TrigMu17=true;

  const int NPtEdges=9;
  float PtEdges[NPtEdges] = {0., 10., 15., 20., 25., 35., 50., 70., 100.};
  float PT=0., PTCorr=0., fEta=0., MTW=0., TightIsoCut=0.;
  float IDSF=1., IsoSF=1., TkSF=1.;
  if     (TightID.Contains("TIsop15")) TightIsoCut = 0.15;
  else if(TightID.Contains("TIsop20")) TightIsoCut = 0.2;
  else if(TightID.Contains("TIsop25")) TightIsoCut = 0.25;

  if(MuLColl.size()>0){
    PT     = MuLColl.at(0).Pt();
    PTCorr = ConeCorrectedPT(MuLColl.at(0), TightIsoCut);
    fEta   = fabs(MuLColl.at(0).Eta());
    MTW    = sqrt(2)*sqrt(MET*MuLColl.at(0).Pt()-METx*MuLColl.at(0).Px()-METy*MuLColl.at(0).Py());
    if(!isData){
//      IDSF   = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuLColl);
//      IsoSF  = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", MuLColl);
      TkSF   = mcdata_correction->MuonTrackingEffScaleFactor(MuLColl);
    }
  }
  weight*=TkSF;

  //Selection------------------------------------------------//
  bool PassJetReq=false, TrigSel=false;
  if( !(MuLColl.size()==1 && EleLColl.size()==0) ) return;

  if     (TrigMu17 && PassTrigger("HLT_Mu17"+AddFilterLabel+"_v") && PT>20 && PTCorr>35) TrigSel=true;
  else if(TrigMu8  && PassTrigger("HLT_Mu8" +AddFilterLabel+"_v") && PT>10 && PTCorr<35) TrigSel=true;
  else return;
  //if     (TrigMu17 && PassTrigger("HLT_Mu17"+AddFilterLabel+"_v") && PT>20 ) TrigSel=true;
  //else if(TrigMu8  && PassTrigger("HLT_Mu8" +AddFilterLabel+"_v") && PT>10 && PT<20) TrigSel=true;
  //else return;
  if(!TrigSel) return;  

  for(int i=0; i<JetColl.size(); i++){
    if(JetColl.at(i).Pt()<40) continue;
    if(MuLColl.at(0).DeltaR(JetColl.at(i))<1.0) continue;
    PassJetReq=true; break;
  }
  if(!PassJetReq) return;

  //if( !(MET<30 && MTW<25) ) return;
  if( !(MET<25 && MTW<25) ) return;
  //---------------------------------------------------------//
  

  bool IsNearB = IsNearBJet(MuLColl.at(0), BJetColl);

  if(fEta<0.9){
    FillHist("MuBSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
    if(IsNearB) FillHist("MuBSumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    else        FillHist("MuBSumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
    if(PassIDCriteria(MuLColl.at(0), TightID)){
      FillHist("MuBIDSumW_PT1D"+Label, PTCorr, weight*IDSF*IsoSF, PtEdges, NPtEdges-1);  
      if(IsNearB) FillHist("MuBIDSumW_NearB_PT1D"+Label, PTCorr, weight*IDSF*IsoSF, PtEdges, NPtEdges-1);
      else        FillHist("MuBIDSumW_AwayB_PT1D"+Label, PTCorr, weight*IDSF*IsoSF, PtEdges, NPtEdges-1);
    }
  }
  else if(fEta<1.6){
    FillHist("MuBESumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
    if(IsNearB) FillHist("MuBESumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    else        FillHist("MuBESumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
    if(PassIDCriteria(MuLColl.at(0), TightID)){
      FillHist("MuBEIDSumW_PT1D"+Label, PTCorr, weight*IDSF*IsoSF, PtEdges, NPtEdges-1);  
      if(IsNearB) FillHist("MuBEIDSumW_NearB_PT1D"+Label, PTCorr, weight*IDSF*IsoSF, PtEdges, NPtEdges-1);
      else        FillHist("MuBEIDSumW_AwayB_PT1D"+Label, PTCorr, weight*IDSF*IsoSF, PtEdges, NPtEdges-1);
    }
  }
  else{
    FillHist("MuESumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
    if(IsNearB)  FillHist("MuESumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    else         FillHist("MuESumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
    if(PassIDCriteria(MuLColl.at(0), TightID)){
      FillHist("MuEIDSumW_PT1D"+Label, PTCorr, weight*IDSF*IsoSF, PtEdges, NPtEdges-1);  
      if(IsNearB)  FillHist("MuEIDSumW_NearB_PT1D"+Label, PTCorr, weight*IDSF*IsoSF, PtEdges, NPtEdges-1);
      else         FillHist("MuEIDSumW_AwayB_PT1D"+Label, PTCorr, weight*IDSF*IsoSF, PtEdges, NPtEdges-1);
    }
  }

}


void Aug2017_MuFakeDataStudy::ValidateID(std::vector<snu::KMuon> MuTColl, float weight, TString Label){
  //Purpose: 1) agreement btw data & MC - bias check

  //Check agreement
  bool PassSel=true; float Mmumu=0.;
  if( MuTColl.size()!=2 ) PassSel=false;
  if( PassSel && !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10) ) PassSel=false;
  if( PassSel &&   MuTColl.at(0).Charge()== MuTColl.at(1).Charge() ) PassSel=false;
  if( PassSel ) Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
  if( PassSel && fabs(Mmumu-91.2)>15 ) PassSel=false;
  
  if(PassSel){
    FillHist("PTmu1"+Label, MuTColl.at(0).Pt(),  weight, 0., 200., 40);
    FillHist("PTmu2"+Label, MuTColl.at(1).Pt(),  weight, 0., 200., 40);
    FillHist("Etamu1"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("Etamu2"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
    FillHist("Mmumu"+Label, Mmumu, weight, 60., 120., 60);
  }
  
}



void Aug2017_MuFakeDataStudy::CheckSSDilepCRs(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetNoVetoColl, std::vector<snu::KJet> BJetNoVetoColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool EMu=Option.Contains("EMu"), MuMu=!EMu;//By default trimu.
  std::vector<snu::KJet> JetColl  = SkimJetColl(JetNoVetoColl,  EleLColl, MuLColl, "EleMuVeto");
  std::vector<snu::KJet> BJetColl = SelBJets(JetColl, "Medium");

  if(MuMu){

  if( !(MuLColl.size()==2 && EleLColl.size()==0) ) return;
  if( !(MuTColl.size()==2) ) return;
  if(  BJetColl.size()==0  ) return;
  if( !(JetColl.size()>=3) ) return;
  
  if( !(MuTColl.at(0).Charge()==MuTColl.at(1).Charge()) ) return;
  float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();

  if( fabs(Mmumu-91.2)<10 || Mmumu<12) return;

  FillHist("PTMu1" +Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
  FillHist("PTMu2" +Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
  FillHist("EtaMu1"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
  FillHist("EtaMu2"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
  FillHist("Mmumu" +Label, Mmumu, weight, 0., 200., 40);
  FillHist("Nj"    +Label, JetColl.size(), weight, 0., 10., 10);
  FillHist("Nb"    +Label, BJetColl.size(), weight, 0., 10., 10);
  FillHist("MET"   +Label, MET, weight, 0., 200., 40);

  }//End of SS MuMu
  else if(EMu){

  if( !(MuLColl.size()==1 && EleLColl.size()==1) ) return;
  if( !(MuTColl.size()==1 && EleTColl.size()==1) ) return;
  if(  BJetColl.size()==0  ) return;
  if( !(JetColl.size()>=3) ) return;
  
  if( !(MuTColl.at(0).Charge()==EleTColl.at(0).Charge()) ) return;

  FillHist("PTMu1" +Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
  FillHist("PTEl1" +Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
  FillHist("EtaMu1"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
  FillHist("EtaEl1"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
  FillHist("Nj"    +Label, JetColl.size(), weight, 0., 10., 10);
  FillHist("Nb"    +Label, BJetColl.size(), weight, 0., 10., 10);
  FillHist("MET"   +Label, MET, weight, 0., 200., 40);

  }//End of SS EMu
}


void Aug2017_MuFakeDataStudy::CheckTrilepCRs(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetNoVetoColl, std::vector<snu::KJet> BJetNoVetoColl, float MET, float METx, float METy, float weight, TString Label, TString Option){

  bool EMuMu=Option.Contains("EMuMu"), TriMu=!EMuMu;//By default trimu.
   std::vector<snu::KJet> JetColl  = SkimJetColl(JetNoVetoColl,  EleLColl, MuLColl, "EleMuVeto");
   std::vector<snu::KJet> BJetColl = SelBJets(JetColl, "Medium");

  if(TriMu){
    FillHist("CutFlow", 0., weight, 0., 10., 10);
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
    FillHist("CutFlow", 1., weight, 0., 10., 10);
  
    int   IdxSSZ    = fabs(MOSSS1-91.2)<fabs(MOSSS2-91.2)? IdxSS1:IdxSS2;
    int   IdxSSNonZ = 3-IdxOS-IdxSSZ;
    float MOSSSZ = (MuTColl.at(IdxOS)+MuTColl.at(IdxSSZ)).M();
    float M3l    = (MuTColl.at(0)+MuTColl.at(1)+MuTColl.at(2)).M();
    float MTW    = sqrt(2)*sqrt(MET*MuTColl.at(IdxSSNonZ).Pt()-METx*MuTColl.at(IdxSSNonZ).Px()-METy*MuTColl.at(IdxSSNonZ).Py());
  
    bool HasBJet=false, OnZ=false, OffZ=false, OnZG=false, WZSel=false, ZGSel=false, ttZSel=false, ZbbSel=false;
    if( BJetColl.size()!=0                 ) HasBJet     = true;
    if( fabs(MOSSSZ-91.2)<10               ) OnZ         = true;
    if( fabs(MOSSS1-91.2)>10 && fabs(MOSSS2-91.2)>10)OffZ= true;
    if( fabs(M3l-91.2)<10                  ) OnZG        = true;
    if( OnZ && M3l>101.2 && MET>50         ) WZSel       = true;
    if( OnZG && OffZ && MET<50      ) ZGSel       = true;
    if( OnZ && HasBJet && JetColl.size()>2 ) ttZSel      = true;
  
    //General 3l Selection
    if(!HasBJet){
      FillHist("CutFlow", 1., weight, 0., 10., 10);

      FillHist("PTmu1_3lOS"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("PTmu2_3lOS"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_3lOS"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_3lOS"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_3lOS"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_3lOS"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSS1_3lOS"+Label, MOSSS1, weight, 0., 200., 40);
      FillHist("MOSSS2_3lOS"+Label, MOSSS2, weight, 0., 200., 40);
      FillHist("M3l_3lOS"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_3lOS"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_3lOS"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_3lOS"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_3lOS"+Label, MTW, weight, 0., 200., 20);
    }
    //3l OnZ ; Fake CR
    if(OnZ){
      FillHist("CutFlow", 2., weight, 0., 10., 10);

      FillHist("PTmu1_3lOSOnZ"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("PTmu2_3lOSOnZ"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_3lOSOnZ"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_3lOSOnZ"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_3lOSOnZ"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_3lOSOnZ"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSSZ_3lOSOnZ"+Label, MOSSSZ, weight, 60., 120., 30);
      FillHist("M3l_3lOSOnZ"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_3lOSOnZ"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_3lOSOnZ"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_3lOSOnZ"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_3lOSOnZ"+Label, MTW, weight, 0., 200., 20);
    }
    //3l OnZ & HasBJet
    if(HasBJet && OnZ){
      FillHist("CutFlow", 3., weight, 0., 10., 10);

      FillHist("PTmu1_3lOSOnZHasB"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_3lOSOnZHasB"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_3lOSOnZHasB"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_3lOSOnZHasB"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_3lOSOnZHasB"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_3lOSOnZHasB"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSSZ_3lOSOnZHasB"+Label, MOSSSZ, weight, 60., 120., 30);
      FillHist("M3l_3lOSOnZHasB"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_3lOSOnZHasB"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_3lOSOnZHasB"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_3lOSOnZHasB"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_3lOSOnZHasB"+Label, MTW, weight, 0., 200., 20);
    }
    if(WZSel){
      FillHist("CutFlow", 4., weight, 0., 10., 10);

      FillHist("PTmu1_WZSel"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("PTmu2_WZSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_WZSel"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_WZSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_WZSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_WZSel"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSSZ_WZSel"+Label, MOSSSZ, weight, 60., 120., 30);
      FillHist("M3l_WZSel"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_WZSel"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_WZSel"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_WZSel"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_WZSel"+Label, MTW, weight, 0., 200., 20);
    }
    if(ZGSel){
      FillHist("CutFlow", 5., weight, 0., 10., 10);

      FillHist("PTmu1_ZGSel"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("PTmu2_ZGSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_ZGSel"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_ZGSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_ZGSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_ZGSel"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSSZ_ZGSel"+Label, MOSSSZ, weight, 0., 200., 40);
      FillHist("M3l_ZGSel"+Label, M3l, weight, 60., 120., 30);
      FillHist("Nj_ZGSel"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_ZGSel"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_ZGSel"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_ZGSel"+Label, MTW, weight, 0., 200., 20);
    }
    if(ttZSel){
      FillHist("CutFlow", 6., weight, 0., 10., 10);

      FillHist("PTmu1_ttZSel"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_ttZSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_ttZSel"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_ttZSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_ttZSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_ttZSel"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSSZ_ttZSel"+Label, MOSSSZ, weight, 60., 120., 30);
      FillHist("M3l_ttZSel"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_ttZSel"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_ttZSel"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_ttZSel"+Label, MET, weight, 0., 300., 30);
      FillHist("MTW_ttZSel"+Label, MTW, weight, 0., 200., 20);
    }
    if(OffZ && BJetColl.size()==1 && JetColl.size()==1){
      FillHist("PTmu1_TopFakeCR"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_TopFakeCR"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_TopFakeCR"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_TopFakeCR"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_TopFakeCR"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_TopFakeCR"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSS1_TopFakeCR"+Label, MOSSS1, weight, 0., 300., 60);
      FillHist("MOSSS2_TopFakeCR"+Label, MOSSS2, weight, 0., 200., 40);
      FillHist("M3l_TopFakeCR"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_TopFakeCR"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_TopFakeCR"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_TopFakeCR"+Label, MET, weight, 0., 300., 30);
      FillHist("MTW_TopFakeCR"+Label, MTW, weight, 0., 200., 20);
    }
    if( ( (!isData) || k_running_nonprompt )
        && OffZ && HasBJet && JetColl.size()>1 ){
      
      FillHist("PTmu1_SR2j"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_SR2j"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_SR2j"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_SR2j"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_SR2j"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_SR2j"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSS1_SR2j"+Label, MOSSS1, weight, 0., 300., 60);
      FillHist("MOSSS2_SR2j"+Label, MOSSS2, weight, 0., 200., 40);
      FillHist("M3l_SR2j"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_SR2j"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_SR2j"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_SR2j"+Label, MET, weight, 0., 300., 30);
      FillHist("MTW_SR2j"+Label, MTW, weight, 0., 200., 20);
    }
    if( ( (!isData) || k_running_nonprompt )
        && OffZ && HasBJet && JetColl.size()>2 ){
      
      FillHist("PTmu1_SR3j"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_SR3j"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_SR3j"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_SR3j"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_SR3j"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_SR3j"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSS1_SR3j"+Label, MOSSS1, weight, 0., 300., 60);
      FillHist("MOSSS2_SR3j"+Label, MOSSS2, weight, 0., 200., 40);
      FillHist("M3l_SR3j"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_SR3j"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_SR3j"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_SR3j"+Label, MET, weight, 0., 300., 30);
      FillHist("MTW_SR3j"+Label, MTW, weight, 0., 200., 20);
    }
  }
  else if(EMuMu){
    if( !(EleLColl.size()==1 && MuLColl.size()==2) ) return;
    if( !(EleTColl.size()==1 && MuTColl.size()==2) ) return;
    if( !(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) ) return;
    if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && MuTColl.at(1).Pt()>10) ) return;
  
    float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
    float MTW = sqrt(2)*sqrt(MET*EleTColl.at(0).Pt()-METx*EleTColl.at(0).Px()-METy*EleTColl.at(0).Py());
    float M3l=(EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M();
    if(Mmumu<12) return;
    
    
    bool HasBJet=false, OnZ=false, OnZG=false, WZSel=false, ZGSel=false, ttZSel=false, ZbbSel=false;
    if( BJetColl.size()!=0 )                      HasBJet = true;
    if( fabs(Mmumu-91.2)<10     )                 OnZ     = true;
    if( fabs(M3l-91.2)<10       )                 OnZG    = true;
    if( OnZ && M3l>101.2 && MET>50 )              WZSel   = true;
    if( OnZG && Mmumu<81.2 && MET<50 )            ZGSel   = true;
    if( OnZ && HasBJet && JetColl.size()>2)       ttZSel  = true;
    if( OnZ && HasBJet && JetColl.size()<3)       ZbbSel  = true;
     
    
    //General 3l Selection
    if(!HasBJet){
      FillHist("PTe_3lOS"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("Etae_3lOS"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("PTmu1_3lOS"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("PTmu2_3lOS"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 200);
      FillHist("Etamu1_3lOS"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("Etamu2_3lOS"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 100);
    
      FillHist("Mmumu_3lOS"+Label, Mmumu, weight, 0., 200., 200);
      FillHist("M3l_3lOS"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 500);
      FillHist("Nj_3lOS"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_3lOS"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_3lOS"+Label, MET, weight, 0., 200., 200);
      FillHist("MTW_3lOS"+Label, MTW, weight, 0., 200., 200);
    }
    //3l OnZ ; Fake CR
    if(OnZ){
      FillHist("PTe_3lOSOnZ"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("Etae_3lOSOnZ"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("PTmu1_3lOSOnZ"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("PTmu2_3lOSOnZ"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 200);
      FillHist("Etamu1_3lOSOnZ"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("Etamu2_3lOSOnZ"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 100);
    
      FillHist("Mmumu_3lOSOnZ"+Label, Mmumu, weight, 0., 200., 200);
      FillHist("M3l_3lOSOnZ"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 500);
      FillHist("Nj_3lOSOnZ"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_3lOSOnZ"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_3lOSOnZ"+Label, MET, weight, 0., 200., 200);
      FillHist("MTW_3lOSOnZ"+Label, MTW, weight, 0., 200., 200);
      if(!OnZG){
        FillHist("PTe_3lOSOnZOffZG"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 200);
        FillHist("Etae_3lOSOnZOffZG"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 100);
        FillHist("PTmu1_3lOSOnZOffZG"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 200);
        FillHist("PTmu2_3lOSOnZOffZG"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 200);
        FillHist("Etamu1_3lOSOnZOffZG"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 100);
        FillHist("Etamu2_3lOSOnZOffZG"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 100);
    
        FillHist("Mmumu_3lOSOnZOffZG"+Label, Mmumu, weight, 0., 200., 200);
        FillHist("M3l_3lOSOnZOffZG"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 500);
        FillHist("Nj_3lOSOnZOffZG"+Label, JetColl.size(), weight, 0., 10., 10);
        FillHist("Nb_3lOSOnZOffZG"+Label, BJetColl.size(), weight, 0., 10., 10);
        FillHist("MET_3lOSOnZOffZG"+Label, MET, weight, 0., 200., 200);
        FillHist("MTW_3lOSOnZOffZG"+Label, MTW, weight, 0., 200., 200);
      }
    }
    //3l OnZ & HasBJet
    if(HasBJet && OnZ){
      FillHist("PTe_3lOSOnZHasB"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("Etae_3lOSOnZHasB"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("PTmu1_3lOSOnZHasB"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("PTmu2_3lOSOnZHasB"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 200);
      FillHist("Etamu1_3lOSOnZHasB"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("Etamu2_3lOSOnZHasB"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 100);
    
      FillHist("Mmumu_3lOSOnZHasB"+Label, Mmumu, weight, 0., 200., 200);
      FillHist("M3l_3lOSOnZHasB"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 500);
      FillHist("Nj_3lOSOnZHasB"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_3lOSOnZHasB"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_3lOSOnZHasB"+Label, MET, weight, 0., 200., 200);
      FillHist("MTW_3lOSOnZHasB"+Label, MTW, weight, 0., 200., 200);
    }
    if(WZSel){
      FillHist("PTe_WZSel"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("Etae_WZSel"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("PTmu1_WZSel"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("PTmu2_WZSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 200);
      FillHist("Etamu1_WZSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("Etamu2_WZSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 100);
    
      FillHist("Mmumu_WZSel"+Label, Mmumu, weight, 0., 200., 200);
      FillHist("M3l_WZSel"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 500);
      FillHist("Nj_WZSel"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_WZSel"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_WZSel"+Label, MET, weight, 0., 200., 200);
      FillHist("MTW_WZSel"+Label, MTW, weight, 0., 200., 200);
    }
    if(ZGSel){
      FillHist("PTe_ZGSel"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("Etae_ZGSel"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("PTmu1_ZGSel"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("PTmu2_ZGSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 200);
      FillHist("Etamu1_ZGSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("Etamu2_ZGSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 100);
    
      FillHist("Mmumu_ZGSel"+Label, Mmumu, weight, 0., 200., 200);
      FillHist("M3l_ZGSel"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 500);
      FillHist("Nj_ZGSel"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_ZGSel"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_ZGSel"+Label, MET, weight, 0., 200., 200);
      FillHist("MTW_ZGSel"+Label, MTW, weight, 0., 200., 200);
    }
    if(ttZSel){
      FillHist("PTe_ttZSel"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("Etae_ttZSel"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("PTmu1_ttZSel"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("PTmu2_ttZSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 200);
      FillHist("Etamu1_ttZSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("Etamu2_ttZSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 100);
    
      FillHist("Mmumu_ttZSel"+Label, Mmumu, weight, 0., 200., 200);
      FillHist("M3l_ttZSel"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 500);
      FillHist("Nj_ttZSel"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_ttZSel"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_ttZSel"+Label, MET, weight, 0., 200., 200);
      FillHist("MTW_ttZSel"+Label, MTW, weight, 0., 200., 200);
    }
    if(ZbbSel){
      FillHist("PTe_ZbbSel"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("Etae_ZbbSel"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("PTmu1_ZbbSel"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("PTmu2_ZbbSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 200);
      FillHist("Etamu1_ZbbSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("Etamu2_ZbbSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 100);
    
      FillHist("Mmumu_ZbbSel"+Label, Mmumu, weight, 0., 200., 200);
      FillHist("M3l_ZbbSel"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 500);
      FillHist("Nj_ZbbSel"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_ZbbSel"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_ZbbSel"+Label, MET, weight, 0., 200., 200);
      FillHist("MTW_ZbbSel"+Label, MTW, weight, 0., 200., 200);
    }
  }

}


void Aug2017_MuFakeDataStudy::DoSystRun(TString Cycle, TString Mode, std::vector<snu::KElectron> EleColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight){

  int SystDir=0; TString SystKindLabel="", SystDirLabel="", TrigOption="", TrigLabel="";
  if(Mode.Contains("Syst")){
    if     (isData)                   return;
    if     (Mode.Contains("Up"))      {SystDir= 1; SystDirLabel="_systup";}
    else if(Mode.Contains("Down"))    {SystDir=-1; SystDirLabel="_systdown";}
    if     (Mode.Contains("PU"))      SystKindLabel="_PU";
    else if(Mode.Contains("Trig"))    SystKindLabel="_Trig";
    else if(Mode.Contains("Nvtx"))    SystKindLabel="_Nvtx";
    else if(Mode.Contains("JES"))     SystKindLabel="_JES";
    else if(Mode.Contains("JER"))     SystKindLabel="_JER";
    else if(Mode.Contains("Uncl"))    SystKindLabel="_Uncl";
    else if(Mode.Contains("ElEn"))    SystKindLabel="_ElEn";
    else if(Mode.Contains("MuEn"))    SystKindLabel="_MuEn";
    else if(Mode.Contains("BTag_L"))  SystKindLabel="_BTag_L";
    else if(Mode.Contains("BTag_BC")) SystKindLabel="_BTag_BC";
  }
  if     (Mode.Contains("Mu8"))     {TrigOption  ="TrigMu8",  TrigLabel="_MuIso8";}
  else if(Mode.Contains("Mu17"))    {TrigOption  ="TrigMu17", TrigLabel="_MuIso17";}


  if(Cycle=="NormCheck"){
    CheckNormCR(MuColl, MuLColl, EleLColl, JetColl, MET, METx, METy, weight, TrigLabel+SystDirLabel+SystKindLabel, TrigOption);
  }//End of NormCheck
  else if(Cycle=="Closure"){
    CheckTrilepCRs(MuColl, MuLColl, EleColl, EleLColl, JetColl, BJetColl, MET, METx, METy, weight, SystDirLabel+SystKindLabel, "");
  }//End of Closure Cycle

  return;
}


void Aug2017_MuFakeDataStudy::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Aug2017_MuFakeDataStudy::BeginCycle() throw( LQError ){
  
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

Aug2017_MuFakeDataStudy::~Aug2017_MuFakeDataStudy() {
  
  Message("In Aug2017_MuFakeDataStudy Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}



float Aug2017_MuFakeDataStudy::NvtxWeight(int Nvtx, TString Option){

  float weight=1.;

  if(Option.Contains("HLT_Mu17_TrkIsoVVL_v")){
    if     (Nvtx<2)  weight=4.01869;
    else if(Nvtx<3)  weight=5.8951;
    else if(Nvtx<4)  weight=5.57468;
    else if(Nvtx<5)  weight=5.44934;
    else if(Nvtx<6)  weight=4.63174;
    else if(Nvtx<7)  weight=4.22565;
    else if(Nvtx<8)  weight=3.63892;
    else if(Nvtx<9)  weight=3.0894;
    else if(Nvtx<10) weight=2.70132;
    else if(Nvtx<11) weight=2.37994;
    else if(Nvtx<12) weight=2.08104;
    else if(Nvtx<13) weight=1.84174;
    else if(Nvtx<14) weight=1.62352;
    else if(Nvtx<15) weight=1.44319;
    else if(Nvtx<16) weight=1.26092;
    else if(Nvtx<17) weight=1.12051;
    else if(Nvtx<18) weight=0.983225;
    else if(Nvtx<19) weight=0.8638;
    else if(Nvtx<20) weight=0.776799;
    else if(Nvtx<21) weight=0.68109;
    else if(Nvtx<22) weight=0.597061;
    else if(Nvtx<23) weight=0.536559;
    else if(Nvtx<24) weight=0.482552;
    else if(Nvtx<25) weight=0.433793;
    else if(Nvtx<26) weight=0.398808;
    else if(Nvtx<27) weight=0.365049;
    else if(Nvtx<28) weight=0.332766;
    else if(Nvtx<29) weight=0.32254;
    else if(Nvtx<30) weight=0.318423;
    else if(Nvtx<31) weight=0.291816;
    else if(Nvtx<32) weight=0.286975;
    else if(Nvtx<33) weight=0.271214;
    else if(Nvtx<34) weight=0.28137;
    else if(Nvtx<35) weight=0.293271;
    else if(Nvtx<36) weight=0.284085;
    else if(Nvtx<37) weight=0.283364;
    else if(Nvtx<38) weight=0.285382;
    else if(Nvtx<39) weight=0.310926;
    else if(Nvtx<40) weight=0.367502;
    else if(Nvtx<41) weight=0.389082;
    else if(Nvtx<42) weight=0.421021;
    else if(Nvtx<43) weight=0.488235;
    else if(Nvtx<44) weight=0.48878;
    else if(Nvtx<45) weight=0.506944;
    else if(Nvtx<46) weight=0.603992;
    else if(Nvtx<47) weight=0.429591;
    else if(Nvtx<48) weight=0.645653;
    else if(Nvtx<49) weight=0.774392;
    else if(Nvtx<50) weight=0.831004;
    else             weight=1.38403;    
  }
  else if(Option.Contains("HLT_Mu8_TrkIsoVVL_v")){
    if     (Nvtx<2)  weight=1.76887;
    else if(Nvtx<3)  weight=3.25332;
    else if(Nvtx<4)  weight=4.25007;
    else if(Nvtx<5)  weight=4.57121;
    else if(Nvtx<6)  weight=4.2404;
    else if(Nvtx<7)  weight=3.59236;
    else if(Nvtx<8)  weight=3.44674;
    else if(Nvtx<9)  weight=2.98316;
    else if(Nvtx<10) weight=2.68859;
    else if(Nvtx<11) weight=2.39977;
    else if(Nvtx<12) weight=2.06791;
    else if(Nvtx<13) weight=1.90538;
    else if(Nvtx<14) weight=1.70796;
    else if(Nvtx<15) weight=1.49451;
    else if(Nvtx<16) weight=1.30118;
    else if(Nvtx<17) weight=1.13331;
    else if(Nvtx<18) weight=1.01193;
    else if(Nvtx<19) weight=0.893801;
    else if(Nvtx<20) weight=0.779264;
    else if(Nvtx<21) weight=0.673095;
    else if(Nvtx<22) weight=0.573211;
    else if(Nvtx<23) weight=0.507485;
    else if(Nvtx<24) weight=0.464434;
    else if(Nvtx<25) weight=0.409781;
    else if(Nvtx<26) weight=0.346145;
    else if(Nvtx<27) weight=0.351931;
    else if(Nvtx<28) weight=0.314329;
    else if(Nvtx<29) weight=0.282073;
    else if(Nvtx<30) weight=0.27587;
    else if(Nvtx<31) weight=0.271491;
    else if(Nvtx<32) weight=0.243343;
    else if(Nvtx<33) weight=0.27028;
    else if(Nvtx<34) weight=0.240851;
    else if(Nvtx<35) weight=0.279104;
    else if(Nvtx<36) weight=0.235309;
    else if(Nvtx<37) weight=0.316236;
    else if(Nvtx<38) weight=0.291917;
    else if(Nvtx<39) weight=0.245509;
    else if(Nvtx<40) weight=0.369302;
    else if(Nvtx<41) weight=0.351038;
    else if(Nvtx<42) weight=0.628974;
    else if(Nvtx<43) weight=0.605933;
    else if(Nvtx<44) weight=0.795525;
    else if(Nvtx<45) weight=0.464135;
    else if(Nvtx<46) weight=0.432921;
    else if(Nvtx<47) weight=0.66136;
    else if(Nvtx<48) weight=0.47839;
    else if(Nvtx<49) weight=1.07925;
    else if(Nvtx<50) weight=0.2382;
    else             weight=1.52761;
  }
  else if(Option.Contains("HLT_Mu17_v")){
    if(Nvtx<2) weight=2.24611;
    else if(Nvtx<3) weight=4.16115;
    else if(Nvtx<4) weight=4.71952;
    else if(Nvtx<5) weight=5.03327;
    else if(Nvtx<6) weight=4.48741;
    else if(Nvtx<7) weight=3.95478;
    else if(Nvtx<8) weight=3.64015;
    else if(Nvtx<9) weight=3.12166;
    else if(Nvtx<10) weight=2.76316;
    else if(Nvtx<11) weight=2.55775;
    else if(Nvtx<12) weight=2.22586;
    else if(Nvtx<13) weight=1.9445;
    else if(Nvtx<14) weight=1.67536;
    else if(Nvtx<15) weight=1.52012;
    else if(Nvtx<16) weight=1.34184;
    else if(Nvtx<17) weight=1.13934;
    else if(Nvtx<18) weight=0.995811;
    else if(Nvtx<19) weight=0.864904;
    else if(Nvtx<20) weight=0.746838;
    else if(Nvtx<21) weight=0.65137;
    else if(Nvtx<22) weight=0.550073;
    else if(Nvtx<23) weight=0.483515;
    else if(Nvtx<24) weight=0.41918;
    else if(Nvtx<25) weight=0.376139;
    else if(Nvtx<26) weight=0.335294;
    else if(Nvtx<27) weight=0.293373;
    else if(Nvtx<28) weight=0.270018;
    else if(Nvtx<29) weight=0.257578;
    else if(Nvtx<30) weight=0.251079;
    else if(Nvtx<31) weight=0.235611;
    else if(Nvtx<32) weight=0.234609;
    else if(Nvtx<33) weight=0.217548;
    else if(Nvtx<34) weight=0.229281;
    else if(Nvtx<35) weight=0.223943;
    else if(Nvtx<36) weight=0.222331;
    else if(Nvtx<37) weight=0.231875;
    else if(Nvtx<38) weight=0.24248;
    else if(Nvtx<39) weight=0.247389;
    else if(Nvtx<40) weight=0.321184;
    else if(Nvtx<41) weight=0.307271;
    else if(Nvtx<42) weight=0.384262;
    else if(Nvtx<43) weight=0.363736;
    else if(Nvtx<44) weight=0.439185;
    else if(Nvtx<45) weight=0.426933;
    else if(Nvtx<46) weight=0.406609;
    else if(Nvtx<47) weight=0.367442;
    else if(Nvtx<48) weight=0.740154;
    else if(Nvtx<49) weight=0.600453;
    else if(Nvtx<50) weight=0.463006;
    else weight=1.2347;
  }
  else if(Option.Contains("HLT_Mu8_v")){
    if(Nvtx<2) weight=1.13637;
    else if(Nvtx<3) weight=1.49427;
    else if(Nvtx<4) weight=4.24496;
    else if(Nvtx<5) weight=4.09766;
    else if(Nvtx<6) weight=4.03281;
    else if(Nvtx<7) weight=3.24854;
    else if(Nvtx<8) weight=3.18818;
    else if(Nvtx<9) weight=2.7531;
    else if(Nvtx<10) weight=2.5599;
    else if(Nvtx<11) weight=2.34081;
    else if(Nvtx<12) weight=2.09831;
    else if(Nvtx<13) weight=1.8428;
    else if(Nvtx<14) weight=1.62385;
    else if(Nvtx<15) weight=1.50268;
    else if(Nvtx<16) weight=1.33354;
    else if(Nvtx<17) weight=1.14302;
    else if(Nvtx<18) weight=1.03362;
    else if(Nvtx<19) weight=0.937813;
    else if(Nvtx<20) weight=0.805245;
    else if(Nvtx<21) weight=0.723832;
    else if(Nvtx<22) weight=0.577172;
    else if(Nvtx<23) weight=0.527943;
    else if(Nvtx<24) weight=0.481102;
    else if(Nvtx<25) weight=0.430203;
    else if(Nvtx<26) weight=0.348326;
    else if(Nvtx<27) weight=0.351447;
    else if(Nvtx<28) weight=0.316442;
    else if(Nvtx<29) weight=0.255616;
    else if(Nvtx<30) weight=0.265576;
    else if(Nvtx<31) weight=0.273567;
    else if(Nvtx<32) weight=0.296904;
    else if(Nvtx<33) weight=0.243986;
    else if(Nvtx<34) weight=0.291674;
    else if(Nvtx<35) weight=0.316465;
    else if(Nvtx<36) weight=0.186769;
    else if(Nvtx<37) weight=0.241667;
    else if(Nvtx<38) weight=0.362795;
    else if(Nvtx<39) weight=0.302989;
    else if(Nvtx<40) weight=0.435794;
    else if(Nvtx<41) weight=0.423558;
    else if(Nvtx<42) weight=0.448057;
    else if(Nvtx<43) weight=0.563595;
    else if(Nvtx<44) weight=0.655235;
    else if(Nvtx<45) weight=0.675944;
    else if(Nvtx<46) weight=0.564465;
    else if(Nvtx<47) weight=0.821991;
    else if(Nvtx<48) weight=0.529438;
    else if(Nvtx<49) weight=0.612416;
    else if(Nvtx<50) weight=1.10675;
    else weight=1.85859;
  }

  return weight;

}

float Aug2017_MuFakeDataStudy::FakeRateData(snu::KMuon Mu, TString Option){

  float FR=0., TightIsoCut=0.;
  if(Option.Contains("ConeSUSY")){
    if     (Option.Contains("TIsop15")) TightIsoCut = 0.15;
    else if(Option.Contains("TIsop20")) TightIsoCut = 0.2;
    else if(Option.Contains("TIsop25")) TightIsoCut = 0.25;
  }
  float PTCorr=ConeCorrectedPT(Mu, TightIsoCut);
    if(Option.Contains("FOPt")) PTCorr=Mu.Pt();
  float fEta=fabs(Mu.Eta());

  if(Option=="POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.103166 ;
      else if(PTCorr<20)  FR=0.0800024;
      else if(PTCorr<25)  FR=0.0847176;
      else if(PTCorr<35)  FR=0.0680578;
      else                FR=0.0658314;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.095482 ;
      else if(PTCorr<20)  FR=0.0713299;
      else if(PTCorr<25)  FR=0.102061 ;
      else if(PTCorr<35)  FR=0.114623 ;
      else                FR=0.0850804;
    }
    else{
      if     (PTCorr<15)  FR=0.0603521;
      else if(PTCorr<20)  FR=0.0918451;
      else if(PTCorr<25)  FR=0.0671922;
      else if(PTCorr<35)  FR=0.128504 ;
      else                FR=0.0964251;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.282211;
      else if(PTCorr<20)  FR=0.2111  ;
      else if(PTCorr<25)  FR=0.187283;
      else if(PTCorr<35)  FR=0.195845;
      else                FR=0.193727;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.321604;
      else if(PTCorr<20)  FR=0.255115;
      else if(PTCorr<25)  FR=0.253692;
      else if(PTCorr<35)  FR=0.197637;
      else                FR=0.238313;
    }
    else{
      if     (PTCorr<15)  FR=0.35875 ;
      else if(PTCorr<20)  FR=0.303843;
      else if(PTCorr<25)  FR=0.298489;
      else if(PTCorr<35)  FR=0.255767;
      else                FR=0.266638;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.276247;
      else if(PTCorr<20)  FR=0.191163;
      else if(PTCorr<25)  FR=0.158219;
      else if(PTCorr<35)  FR=0.155122;
      else                FR=0.153071;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.309303;
      else if(PTCorr<20)  FR=0.219685;
      else if(PTCorr<25)  FR=0.211428;
      else if(PTCorr<35)  FR=0.170437;
      else                FR=0.19083 ;
    }
    else{
      if     (PTCorr<15)  FR=0.339582;
      else if(PTCorr<20)  FR=0.263034;
      else if(PTCorr<25)  FR=0.244504;
      else if(PTCorr<35)  FR=0.222135;
      else                FR=0.219002;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1NoChi_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.104934 ;
      else if(PTCorr<20)  FR=0.0810565;
      else if(PTCorr<25)  FR=0.085597 ;
      else if(PTCorr<35)  FR=0.0684988;
      else                FR=0.0665389;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.0977776;
      else if(PTCorr<20)  FR=0.0729936;
      else if(PTCorr<25)  FR=0.10299  ;
      else if(PTCorr<35)  FR=0.116561 ;
      else                FR=0.0859202;
    }
    else{
      if     (PTCorr<15)  FR=0.0610895;
      else if(PTCorr<20)  FR=0.0923621;
      else if(PTCorr<25)  FR=0.0671893;
      else if(PTCorr<35)  FR=0.128802 ;
      else                FR=0.0967269;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1NoChi_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.308976;
      else if(PTCorr<20)  FR=0.229748;
      else if(PTCorr<25)  FR=0.202024;
      else if(PTCorr<35)  FR=0.211457;
      else                FR=0.215925;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.351577;
      else if(PTCorr<20)  FR=0.275117;
      else if(PTCorr<25)  FR=0.266206;
      else if(PTCorr<35)  FR=0.207133;
      else                FR=0.252437;
    }
    else{
      if     (PTCorr<15)  FR=0.364771;
      else if(PTCorr<20)  FR=0.309155;
      else if(PTCorr<25)  FR=0.302773;
      else if(PTCorr<35)  FR=0.265689;
      else                FR=0.275625;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1NoChi_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.30168 ;
      else if(PTCorr<20)  FR=0.205724;
      else if(PTCorr<25)  FR=0.167485;
      else if(PTCorr<35)  FR=0.163753;
      else                FR=0.16684 ;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.336895;
      else if(PTCorr<20)  FR=0.234467;
      else if(PTCorr<25)  FR=0.219452;
      else if(PTCorr<35)  FR=0.176959;
      else                FR=0.200026;
    }
    else{
      if     (PTCorr<15)  FR=0.34518 ;
      else if(PTCorr<20)  FR=0.267018;
      else if(PTCorr<25)  FR=0.247139;
      else if(PTCorr<35)  FR=0.228795;
      else                FR=0.225253;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1NoChi_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.0730267;
      else if(PTCorr<20)  FR=0.0422678;
      else if(PTCorr<25)  FR=0.0394787;
      else if(PTCorr<35)  FR=0.0292311;
      else                FR=0.0247411;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.075509 ;
      else if(PTCorr<20)  FR=0.0421607;
      else if(PTCorr<25)  FR=0.0535583;
      else if(PTCorr<35)  FR=0.0579744;
      else                FR=0.0374653;
    }
    else{
      if     (PTCorr<15)  FR=0.0488572;
      else if(PTCorr<20)  FR=0.0598638;
      else if(PTCorr<25)  FR=0.036773 ;
      else if(PTCorr<35)  FR=0.0710471;
      else                FR=0.0479654;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1NoChi_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.243067 ;
      else if(PTCorr<20)  FR=0.120063 ;
      else if(PTCorr<25)  FR=0.100958 ;
      else if(PTCorr<35)  FR=0.0993144;
      else                FR=0.0935684;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.288401;
      else if(PTCorr<20)  FR=0.159155;
      else if(PTCorr<25)  FR=0.150465;
      else if(PTCorr<35)  FR=0.106834;
      else                FR=0.122566;
    }
    else{
      if     (PTCorr<15)  FR=0.308748;
      else if(PTCorr<20)  FR=0.195025;
      else if(PTCorr<25)  FR=0.182063;
      else if(PTCorr<35)  FR=0.153384;
      else                FR=0.148075;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1NoChi_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.23621  ;
      else if(PTCorr<20)  FR=0.1076   ;
      else if(PTCorr<25)  FR=0.0816986;
      else if(PTCorr<35)  FR=0.0747883;
      else                FR=0.068686 ;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.275343 ;
      else if(PTCorr<20)  FR=0.135587 ;
      else if(PTCorr<25)  FR=0.121004 ;
      else if(PTCorr<35)  FR=0.0901332;
      else                FR=0.0936505;
    }
    else{
      if     (PTCorr<15)  FR=0.291057;
      else if(PTCorr<20)  FR=0.169294;
      else if(PTCorr<25)  FR=0.145199;
      else if(PTCorr<35)  FR=0.130402;
      else                FR=0.118149;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.220591 ;
      else if(PTCorr<20)  FR=0.0957646;
      else if(PTCorr<25)  FR=0.0756069;
      else if(PTCorr<35)  FR=0.0664369;
      else                FR=0.0596936;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.251618 ;
      else if(PTCorr<20)  FR=0.115786 ;
      else if(PTCorr<25)  FR=0.102281 ;
      else if(PTCorr<35)  FR=0.0847535;
      else                FR=0.0776838;
    }
    else{
      if     (PTCorr<15)  FR=0.278727;
      else if(PTCorr<20)  FR=0.147944;
      else if(PTCorr<25)  FR=0.125158;
      else if(PTCorr<35)  FR=0.116875;
      else                FR=0.10205 ;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1sig4NoChi_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.333358;
      else if(PTCorr<20)  FR=0.236995;
      else if(PTCorr<25)  FR=0.23349 ;
      else if(PTCorr<35)  FR=0.17589 ;
      else                FR=0.168586;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.323098;
      else if(PTCorr<20)  FR=0.207417;
      else if(PTCorr<25)  FR=0.261432;
      else if(PTCorr<35)  FR=0.288036;
      else                FR=0.21438 ;
    }
    else{
      if     (PTCorr<15)  FR=0.191651;
      else if(PTCorr<20)  FR=0.237084;
      else if(PTCorr<25)  FR=0.170018;
      else if(PTCorr<35)  FR=0.283506;
      else                FR=0.221901;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1sig4NoChi_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.434692;
      else if(PTCorr<20)  FR=0.303991;
      else if(PTCorr<25)  FR=0.246176;
      else if(PTCorr<35)  FR=0.265866;
      else                FR=0.267146;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.472631;
      else if(PTCorr<20)  FR=0.360498;
      else if(PTCorr<25)  FR=0.33543 ;
      else if(PTCorr<35)  FR=0.252456;
      else                FR=0.315846;
    }
    else{
      if     (PTCorr<15)  FR=0.463157;
      else if(PTCorr<20)  FR=0.373889;
      else if(PTCorr<25)  FR=0.363318;
      else if(PTCorr<35)  FR=0.321106;
      else                FR=0.339224;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1sig4NoChi_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.432992;
      else if(PTCorr<20)  FR=0.2985  ;
      else if(PTCorr<25)  FR=0.244148;
      else if(PTCorr<35)  FR=0.248196;
      else                FR=0.247938;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.468875;
      else if(PTCorr<20)  FR=0.344572;
      else if(PTCorr<25)  FR=0.323177;
      else if(PTCorr<35)  FR=0.259591;
      else                FR=0.29687 ;
    }
    else{
      if     (PTCorr<15)  FR=0.455616;
      else if(PTCorr<20)  FR=0.359923;
      else if(PTCorr<25)  FR=0.338585;
      else if(PTCorr<35)  FR=0.314778;
      else                FR=0.318862;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.235023 ;
      else if(PTCorr<20)  FR=0.126737 ;
      else if(PTCorr<25)  FR=0.10877  ;
      else if(PTCorr<35)  FR=0.0810774;
      else                FR=0.0657765;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.23597  ;
      else if(PTCorr<20)  FR=0.121983 ;
      else if(PTCorr<25)  FR=0.135113 ;
      else if(PTCorr<35)  FR=0.146711 ;
      else                FR=0.0987709;
    }
    else{
      if     (PTCorr<15)  FR=0.150243 ;
      else if(PTCorr<20)  FR=0.151052 ;
      else if(PTCorr<25)  FR=0.0953203;
      else if(PTCorr<35)  FR=0.159848 ;
      else                FR=0.114909 ;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.345149;
      else if(PTCorr<20)  FR=0.164104;
      else if(PTCorr<25)  FR=0.124837;
      else if(PTCorr<35)  FR=0.126013;
      else                FR=0.118599;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.389373;
      else if(PTCorr<20)  FR=0.208879;
      else if(PTCorr<25)  FR=0.189129;
      else if(PTCorr<35)  FR=0.131404;
      else                FR=0.155852;
    }
    else{
      if     (PTCorr<15)  FR=0.39351 ;
      else if(PTCorr<20)  FR=0.237086;
      else if(PTCorr<25)  FR=0.217314;
      else if(PTCorr<35)  FR=0.183934;
      else                FR=0.185743;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.342129;
      else if(PTCorr<20)  FR=0.161717;
      else if(PTCorr<25)  FR=0.124337;
      else if(PTCorr<35)  FR=0.124102;
      else                FR=0.112075;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.3853  ;
      else if(PTCorr<20)  FR=0.200014;
      else if(PTCorr<25)  FR=0.181963;
      else if(PTCorr<35)  FR=0.143144;
      else                FR=0.149275;
    }
    else{
      if     (PTCorr<15)  FR=0.386585;
      else if(PTCorr<20)  FR=0.228739;
      else if(PTCorr<25)  FR=0.203977;
      else if(PTCorr<35)  FR=0.190597;
      else                FR=0.180573;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.323899 ;
      else if(PTCorr<20)  FR=0.145267 ;
      else if(PTCorr<25)  FR=0.114264 ;
      else if(PTCorr<35)  FR=0.102888 ;
      else                FR=0.0929075;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.355195;
      else if(PTCorr<20)  FR=0.170462;
      else if(PTCorr<25)  FR=0.15116 ;
      else if(PTCorr<35)  FR=0.127486;
      else                FR=0.118954;
    }
    else{
      if     (PTCorr<15)  FR=0.370076;
      else if(PTCorr<20)  FR=0.199035;
      else if(PTCorr<25)  FR=0.170545;
      else if(PTCorr<35)  FR=0.161002;
      else                FR=0.146339;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.276577 ;
      else if(PTCorr<20)  FR=0.142112 ;
      else if(PTCorr<25)  FR=0.117096 ;
      else if(PTCorr<35)  FR=0.0858086;
      else                FR=0.0684839;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.351973;
      else if(PTCorr<20)  FR=0.154119;
      else if(PTCorr<25)  FR=0.154107;
      else if(PTCorr<35)  FR=0.162744;
      else                FR=0.109306;
    }
    else{
      if     (PTCorr<15)  FR=0.325502;
      else if(PTCorr<20)  FR=0.252336;
      else if(PTCorr<25)  FR=0.144032;
      else if(PTCorr<35)  FR=0.212277;
      else                FR=0.151429;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.363658;
      else if(PTCorr<20)  FR=0.169303;
      else if(PTCorr<25)  FR=0.127615;
      else if(PTCorr<35)  FR=0.128469;
      else                FR=0.121185;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.423617;
      else if(PTCorr<20)  FR=0.219131;
      else if(PTCorr<25)  FR=0.196251;
      else if(PTCorr<35)  FR=0.134896;
      else                FR=0.160779;
    }
    else{
      if     (PTCorr<15)  FR=0.475586;
      else if(PTCorr<20)  FR=0.269348;
      else if(PTCorr<25)  FR=0.242709;
      else if(PTCorr<35)  FR=0.19999 ;
      else                FR=0.201519;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.362155;
      else if(PTCorr<20)  FR=0.167222;
      else if(PTCorr<25)  FR=0.125889;
      else if(PTCorr<35)  FR=0.120206;
      else                FR=0.110192;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.422128;
      else if(PTCorr<20)  FR=0.213486;
      else if(PTCorr<25)  FR=0.189325;
      else if(PTCorr<35)  FR=0.140276;
      else                FR=0.151332;
    }
    else{
      if     (PTCorr<15)  FR=0.473036;
      else if(PTCorr<20)  FR=0.268132;
      else if(PTCorr<25)  FR=0.232456;
      else if(PTCorr<35)  FR=0.201764;
      else                FR=0.193886;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.341034 ;
      else if(PTCorr<20)  FR=0.151259 ;
      else if(PTCorr<25)  FR=0.118491 ;
      else if(PTCorr<35)  FR=0.106629 ;
      else                FR=0.0954289;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.38962 ;
      else if(PTCorr<20)  FR=0.183194;
      else if(PTCorr<25)  FR=0.160408;
      else if(PTCorr<35)  FR=0.133272;
      else                FR=0.124299;
    }
    else{
      if     (PTCorr<15)  FR=0.455303;
      else if(PTCorr<20)  FR=0.232687;
      else if(PTCorr<25)  FR=0.197538;
      else if(PTCorr<35)  FR=0.182767;
      else                FR=0.163027;
    }
  }
  else {cout<<"No Such FR!"<<endl;}



  return FR;
}


float  Aug2017_MuFakeDataStudy::GetFakeWeight(std::vector<snu::KMuon>& MuLColl, TString MuLID, TString MuTID, vector<snu::KJet>& BJetNoVetoColl, TString Option){

  float fakeweight=-1.; int NLooseNotTight=0;

  TString FilterInfo="", ConeMethod="";
  if     (Option.Contains("NoFilter"))  FilterInfo="NoFilter";
  else if(Option.Contains("TrkIsoVVL")) FilterInfo="TrkIsoVVL";
  if     (Option.Contains("ConeE"))     ConeMethod="ConeE";
  else if(Option.Contains("ConeSUSY"))  ConeMethod="ConeSUSY";
  else if(Option.Contains("FOPt"))      ConeMethod="FOPt";


  for(int i=0; i<(int) MuLColl.size(); i++){
    if(!PassIDCriteria(MuLColl.at(i), MuTID, "Roch")){
      float FR=0.;
      bool IsNearB = IsNearBJet(MuLColl.at(i), BJetNoVetoColl);
      TString FlavOpt = IsNearB? "_NearB":"_AwayB";

      FR=FakeRateData(MuLColl.at(i),MuTID+"_"+MuLID+FlavOpt+"_"+FilterInfo+ConeMethod);
      //cout<<i<<FlavOpt<<" MuFR"<<FR<<endl;
      fakeweight*=-FR/(1.-FR);
      NLooseNotTight++;
    }
  }
  if(NLooseNotTight==0) return 0.;

  return fakeweight;


}

float Aug2017_MuFakeDataStudy::GetFakeWeight(std::vector<snu::KMuon>& MuLColl, TString MuLID, TString MuTID, TString Option){

  TString FilterInfo="", ConeMethod="";
  if     (Option.Contains("NoFilter"))  FilterInfo="NoFilter";
  else if(Option.Contains("TrkIsoVVL")) FilterInfo="TrkIsoVVL";
  if     (Option.Contains("ConeE"))     ConeMethod="ConeE";
  else if(Option.Contains("ConeSUSY"))  ConeMethod="ConeSUSY";
  else if(Option.Contains("FOPt"))      ConeMethod="FOPt";

  float fakeweight=-1.; int NLooseNotTight=0;
  for(int i=0; i<MuLColl.size(); i++){
    if(!PassIDCriteria(MuLColl.at(i), MuTID,"Roch")){
      float FR=0.;
      FR=FakeRateData(MuLColl.at(i),MuTID+"_"+MuLID+"_"+FilterInfo+ConeMethod);
      fakeweight*=-FR/(1.-FR);
      NLooseNotTight++;
    }
  }
  if(NLooseNotTight==0) return 0.;

  return fakeweight;
}



int Aug2017_MuFakeDataStudy::GetFakeLepSrcIdx(snu::KMuon Mu, std::vector<snu::KTruth> TruthColl){

  int MatchedIdx=-1;
  int LepType=GetLeptonType(Mu,TruthColl);
  if(LepType==-1){
    float mindR=999., maxdR=0.4; int IdxMindR=-1;
    for(int j=2; j<TruthColl.size(); j++){
      if( TruthColl.at(j).GenStatus()!=1 ) continue;
      if( fabs(TruthColl.at(j).PdgId())<50 ) continue;
      if( !(TruthColl.at(j).Pt()>5 && fabs(TruthColl.at(j).Eta())<2.5) ) continue;
      if( Mu.DeltaR(TruthColl.at(j))<mindR ){
        mindR=Mu.DeltaR(TruthColl.at(j));
        IdxMindR=j;
      }
    }
    if(mindR<maxdR){ MatchedIdx=IdxMindR; }
  }
  if(LepType==-2 || LepType==-3){
    int MatchedTruthIdx = GenMatchedIdx(Mu,TruthColl);
    int MotherIdx       = FirstNonSelfMotherIdx(MatchedTruthIdx,TruthColl);
    int GrMotherIdx     = FirstNonSelfMotherIdx(MotherIdx,TruthColl);
    
    if(LepType==-2) MatchedIdx=MotherIdx;
    if(LepType==-3) MatchedIdx=GrMotherIdx;    
  }

  return MatchedIdx;
}



int Aug2017_MuFakeDataStudy::GetFakeLepJetSrcType(snu::KMuon Mu, std::vector<snu::KJet> JetColl){
  //Type: -1: Unmatched, 1:L, 2:C, 3:B
  int SrcType=-1;
  bool NearB=false, NearC=false, NearL=false;
  for(int i=0; i<JetColl.size(); i++){
    if(Mu.DeltaR(JetColl.at(i))<0.4){
      if(JetColl.at(i).HadronFlavour()==5){ NearB=true; break; }//1)
      else if(JetColl.at(i).HadronFlavour()==4){ NearC=true; }
      else if(JetColl.at(i).HadronFlavour()==0){ NearL=true; }
    }
  }

  if     (NearB) SrcType=3;
  else if(NearC) SrcType=2;
  else if(NearL) SrcType=1;

  return SrcType;
//1) Higher Priority to B. if there's multiple near jets, then b-jet has higher priority
}


float Aug2017_MuFakeDataStudy::ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut){

  //float PTCorr=Ele.Pt()*(1+Ele.PFRelIso(0.3));
  float PTCorr=Ele.Pt()*(1+min((float) 0,(float) Ele.PFRelIso(0.3)-TightIsoCut));
  return PTCorr;

}


float Aug2017_MuFakeDataStudy::ConeCorrectedPT(snu::KMuon Mu, float TightIsoCut){

  //float PTCorr=Mu.Pt();
  //float PTCorr=Mu.Pt()*(1.+RochIso(Mu,"0.4"));
  float PTCorr=Mu.Pt()*(1.+max((float) 0.,RochIso(Mu,"0.4")-TightIsoCut));
  return PTCorr;

}

bool Aug2017_MuFakeDataStudy::IsNearBJet(snu::KMuon Mu, std::vector<snu::KJet> bjetNoVetoColl){

  bool IsNearB=false;
  for(int i=0; i<bjetNoVetoColl.size(); i++){
    if(Mu.DeltaR(bjetNoVetoColl.at(i))<0.4){
      IsNearB=true; break;
    }
  }

  return IsNearB;
}




int Aug2017_MuFakeDataStudy::StepPassed(std::vector<snu::KMuon> MuColl, std::vector<snu::KElectron> EleColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, TString Option){

   int PassedSteps=0;
   bool EMuMu=false, TriMu=false;
   if     (Option.Contains("EMuMu")) EMuMu=true;
   else if(Option.Contains("TriMu")) TriMu=true;

   if(EMuMu){
     //Step1
     if(EleColl.size()==1 && MuColl.size()==2){
        if(EleColl.at(0).Pt()>25 && MuColl.at(0).Pt()>15 && MuColl.at(1).Pt()>10){PassedSteps++;} else return PassedSteps;
     }else return PassedSteps;
     

     //Step2
     if(MuColl.at(0).Charge()!=MuColl.at(1).Charge()) {PassedSteps++;} else return PassedSteps;
     //Step3
     if((MuColl.at(0)+MuColl.at(1)).M()>12)           {PassedSteps++;} else return PassedSteps;
     //Step4
     if(BJetColl.size()>=1)                           {PassedSteps++;} else return PassedSteps;
     //Step5
     if(JetColl.size()>=3)                            {PassedSteps++;} else return PassedSteps;
     //Step6
     if((MuColl.at(0)+MuColl.at(1)).M()>12 && (MuColl.at(0)+MuColl.at(1)).M()<40) {PassedSteps++;} else return PassedSteps;
   }
   else if(TriMu){
     //Step1
     if(EleColl.size()==0 && MuColl.size()==3){
       if(MuColl.at(0).Pt()>20 && MuColl.at(1).Pt()>10 && MuColl.at(2).Pt()>10 ){PassedSteps++;} else return PassedSteps;
     }else return PassedSteps;

     //Step2
     if(fabs(SumCharge(MuColl))==1){ PassedSteps++;} else return PassedSteps;

     int Idx_OS=TriMuChargeIndex(MuColl,"OS");
     int Idx_SS1=TriMuChargeIndex(MuColl,"SS1");
     int Idx_SS2=TriMuChargeIndex(MuColl,"SS2");
     float MOSSS1=(MuColl.at(Idx_OS)+MuColl.at(Idx_SS1)).M();
     float MOSSS2=(MuColl.at(Idx_OS)+MuColl.at(Idx_SS2)).M();
     //Step3
     if(MOSSS1>12 && MOSSS2>12){PassedSteps++;} else return PassedSteps;
     //Step4
     if(BJetColl.size()>=1){PassedSteps++;} else return PassedSteps;
     //Step5
     if(JetColl.size()>=3){PassedSteps++;} else return PassedSteps;
    
   }

   return PassedSteps;
}


int Aug2017_MuFakeDataStudy::NPromptFake_Mu(std::vector<snu::KMuon> MuColl, std::vector<snu::KTruth> truthColl, TString Option){

   int Nprompt=0, Nfake=0;

   bool ReturnPrompt=false, ReturnFake=false;
   if(Option.Contains("EWPrompt")) ReturnPrompt=true;
   if(Option.Contains("HFake"))    ReturnFake=true;

   for(int i=0; i<MuColl.size(); i++){
     int MuType=GetLeptonType(MuColl.at(i), truthColl);
     //if     (MuType>0 ) Nprompt++;

     if     (MuType>0 && MuType<4) Nprompt++;
     else if(MuType<0 && MuType>-5 ) Nfake++;
//     cout<<i<<" MuType "<<MuType<<endl;
   }

//   cout<<"Tot "<<Nprompt+Nfake<<" NPr "<<Nprompt<<" Nfake "<<Nfake<<endl;
   if(ReturnPrompt && ReturnFake) return Nprompt+Nfake;
   else if(ReturnPrompt) return Nprompt;
   else if(ReturnFake)   return Nfake;

   return 0.;
  
}

int Aug2017_MuFakeDataStudy::NPromptFake_Ele(std::vector<snu::KElectron> EleColl, std::vector<snu::KTruth> truthColl, TString Option){

   int Nprompt=0, Nfake=0;

   bool ReturnPrompt=false, ReturnFake=false, ReturnConv=false;
   if(Option.Contains("EWPrompt")) ReturnPrompt=true;
   if(Option.Contains("HFake"))    ReturnFake=true;
   if(Option.Contains("Conv"))     ReturnConv=true;

   for(int i=0; i<EleColl.size(); i++){
     int EleType=GetLeptonType(EleColl.at(i), truthColl);
     //if     (EleType>0 )  Nprompt++;

     if     (EleType>0 && EleType<4)  Nprompt++;
     else if(EleType<0 && EleType>-5) Nfake++;
     //else if(EleType<-4)              Nfake++;

//     cout<<i<<" ElType "<<EleType<<endl;
   }


//   cout<<"Tot "<<Nprompt+Nfake<<" NPr "<<Nprompt<<" Nfake "<<Nfake<<endl;
   if(ReturnPrompt && ReturnFake) return Nprompt+Nfake;
   else if(ReturnPrompt) return Nprompt;
   else if(ReturnFake)   return Nfake;

   return 0.;
  
}



bool Aug2017_MuFakeDataStudy::IsConvCand(snu::KElectron Ele, std::vector<snu::KTruth> TruthColl, TString Option){

  bool IsConversionCandidate=false;
  int EleType=GetLeptonType(Ele, TruthColl);
  int HardPhotonIdx=-1;
  if(EleType==1 || fabs(EleType)==2 || fabs(EleType)==3 ) return false; 

  for(int i=2; i<TruthColl.size(); i++){

    int pid=TruthColl.at(i).PdgId();
    int midx=TruthColl.at(i).IndexMother();
    int mpid=TruthColl.at(midx).PdgId();
    int GenSt=TruthColl.at(i).GenStatus();

    if( fabs(pid)==22 && ((GenSt>20 && GenSt<30) || GenSt==1 ) ) HardPhotonIdx=i;
  }
  if(HardPhotonIdx<2) return false;
  
  if(TruthColl.at(HardPhotonIdx).DeltaR(Ele)<0.4) IsConversionCandidate=true;

  return IsConversionCandidate;

}



void Aug2017_MuFakeDataStudy::FillCutFlow(TString cut, float weight){
  
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
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(5,"3lCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(6,"M(#mu#mu)>4");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(7,"#geq1b");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(8,"2j");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(9,"3j");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(10,"ZVeto");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(11,"Nlj2");
    }
    if(!GetHist("cutflow_N")){
      AnalyzerCore::MakeHistograms("cutflow_N", 11, 0., 11.);
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(1,"NoCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(2,"TriggerCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(3,"EventCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(4,"VertexCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(5,"3lCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(6,"M(#mu#mu)>4");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(7,"#geq1b");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(8,"2j");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(9,"3j");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(10,"ZVeto");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(11,"Nlj2");
    }
  }
}



void Aug2017_MuFakeDataStudy::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Aug2017_MuFakeDataStudy::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //Basic Hists/////////////////////////////////////////////////
  //Data properties before selection  
  AnalyzerCore::MakeHistograms("Basic_b2_Et_orig", 200, 0., 200.);



  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Aug2017_MuFakeDataStudyCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Aug2017_MuFakeDataStudy::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
