// $Id: Feb2018_TrigAccept.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQFeb2018_TrigAccept Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Feb2018_TrigAccept.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Feb2018_TrigAccept);

 Feb2018_TrigAccept::Feb2018_TrigAccept() : AnalyzerCore(), out_muons(0) {

   SetLogName("Feb2018_TrigAccept");
   Message("In Feb2018_TrigAccept constructor", INFO);
   InitialiseAnalysis();
 }


 void Feb2018_TrigAccept::InitialiseAnalysis() throw( LQError ) {
   
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

void Feb2018_TrigAccept::ExecuteEvents()throw( LQError ){

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
   float top_pt_reweight=1.;
   if(!k_isdata){
      pileup_reweight          = eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central); 
      k_factor_weight = GetKFactor();
      geneff_weight   = GenFilterEfficiency(k_sample_name);
      gennorm_weight  = SignalNorm(k_sample_name, 20.);
      if(k_sample_name.Contains("TT_powheg")){
       std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
       top_pt_reweight = TopPTReweight(truthColl);
      }
   }
   weight*=pileup_reweight;
   FillHist("Basic_PURW", pileup_reweight, 1., 0., 20., 200);

 
   //Total Event(MCweight+PUrw)////////////////
   FillCutFlow("NoCut", weight*pileup_reweight);


   bool DoubleMuon=false, ElectronMuon=false, CheckTrigAccept=false;
   bool EMuMuEComp=false;
   bool SystRun=false;
   for(int i=0; i<(int) k_flags.size(); i++){
     if     (k_flags.at(i).Contains("DoubleMuon"))      DoubleMuon      = true;
     else if(k_flags.at(i).Contains("ElectronMuon"))    ElectronMuon    = true;
     else if(k_flags.at(i).Contains("SystRun"))         SystRun         = true;
     else if(k_flags.at(i).Contains("CheckTrigAccept")) CheckTrigAccept = true;
     else if(k_flags.at(i).Contains("EMuMuEComp"))      EMuMuEComp      = true;
   }


    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //Trigger Path of Analysis
   bool Pass_Trigger=false, Pass_TriggerBG=false, Pass_TriggerH=false;
   float trigger_ps_weight=1., trigger_period_weight=1.;
   float LumiBG=27.257618, LumiH=8.605696, LumiBH=35.863314;
   
   if(!isData) trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
   weight*=trigger_ps_weight*trigger_period_weight;

//   if(!Pass_Trigger) return;
   /**********************************************************************************/

   //METFilterCut : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
   if(!PassMETFilter()) return;  FillCutFlow("EventCut", weight*pileup_reweight);

   //Vertex Cut
   //(vtx.ndof>4&&maxAbsZ<=0)||abs(vtx.z)<= 24)&&((maxd0 <=0)||abs(vtx.position.rho)<=2)&&!(vtx.isFake))
   if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; FillCutFlow("VertexCut", weight*pileup_reweight);
   


   //===========================================================================//
   // Objects Selection in Analysis                                             //
   //===========================================================================//

   //Primary Object Collection
   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);

   //weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*btag_sf*pileup_reweight*fake_weight;
   //-----------------------------------------------------------------------------------------//



   //----------------------------------------------------------------------------------//
   //==================================================================================//
   /////// MAIN ANALYSIS CODE ///////////////////////////////////////////////////////////
   //==================================================================================//
   //----------------------------------------------------------------------------------//

   if(EMuMuEComp){
     bool Pass_TriggerMuE=false, Pass_TriggerEMu=false;
     bool Pass_TriggerBG =false, Pass_TriggerH  =false;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") )    {Pass_TriggerEMu=true; Pass_TriggerBG=true;}
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) {Pass_TriggerEMu=true; Pass_TriggerH =true;}
     if( PassTrigger("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v") )    {Pass_TriggerMuE=true; Pass_TriggerBG=true;}
     if( PassTrigger("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v") ) {Pass_TriggerMuE=true; Pass_TriggerH =true;}
     if( !(Pass_TriggerMuE || Pass_TriggerEMu) ) return;

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
       eventbase->GetElectronSel()->SetPt(10.,15);             eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.4, 0.4);
       eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetdxySigMax(4.);
       eventbase->GetElectronSel()->SetApplyConvVeto(true);
     std::vector<snu::KElectron> electronTempLColl; eventbase->GetElectronSel()->Selection(electronTempLColl);
       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_HctoWA_FAKELOOSE);
       eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
       eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.4, 0.4);
       eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetdxySigMax(4.);
       eventbase->GetElectronSel()->SetApplyConvVeto(true);
     std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
       for(int i=0; i<(int) electronTempLColl.size(); i++){ electronLooseColl.push_back(electronTempLColl.at(i)); }
       

       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
       eventbase->GetElectronSel()->SetPt(10.,15.);            eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
       eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetdxySigMax(4.);
       eventbase->GetElectronSel()->SetApplyConvVeto(true);
     std::vector<snu::KElectron> electronTempTColl; eventbase->GetElectronSel()->Selection(electronTempTColl);
       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
       eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
       eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
       eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetdxySigMax(4.);
       eventbase->GetElectronSel()->SetApplyConvVeto(true);
     std::vector<snu::KElectron> electronTightColl; eventbase->GetElectronSel()->Selection(electronTightColl);
       for(int i=0; i<(int) electronTempTColl.size(); i++){ electronTightColl.push_back(electronTempTColl.at(i)); }
     std::vector<snu::KElectron> electronColl;  if(!k_running_nonprompt){ electronColl=electronTightColl;} else{ electronColl=electronLooseColl;}

       bool LeptonVeto=false;
       eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
       eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     std::vector<snu::KJet> jetNoVetoColl; eventbase->GetJetSel()->Selection(jetNoVetoColl, LeptonVeto, muonLooseColl, electronLooseColl);
     std::vector<snu::KJet> bjetNoVetoColl = SelBJets(jetNoVetoColl, "Medium");
     std::vector<snu::KJet> jetColl  = SkimJetColl(jetNoVetoColl,  electronLooseColl, muonLooseColl, "EleMuVeto");
     std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");

     float met    = eventbase->GetEvent().MET(snu::KEvent::pfmet);
     float metphi = eventbase->GetEvent().METPhi();

     float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1., btag_sf=1.;
     float trigger_sf=1.;
     float fake_weight=1.; bool EventCand=false;
  
     /*This part is for boosting up speed.. SF part takes rather longer time than expected*/
     if(!isData){
       //id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl);
       reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
       //id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl);
       trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
       btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);
       //float trigger_sf1 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
       //float trigger_sf2 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
       //trigger_sf        = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
     }
     weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*btag_sf*trigger_sf;

     if(Pass_TriggerEMu){
       if( !(electronLooseColl.size()==1 && muonLooseColl.size()==2) ) return;
       if( !(electronColl.size()==1      && muonColl.size()==2     ) ) return;
       if( !(muonColl.at(0).Charge()     != muonColl.at(1).Charge()) ) return;
       if( !(electronColl.at(0).Pt()>25  && muonColl.at(0).Pt()>10 && muonColl.at(1).Pt()>10) ) return;
       if( fabs(electronColl.at(0).Eta())>2.5 ) return;

       float Mmumu=(muonColl.at(0)+muonColl.at(1)).M();
        
       if(Mmumu<12) return;
         FillHist("CutFlow_EMu_1e2mu", 0., weight, 0., 20., 20);
       if(fabs(Mmumu-91.2)<10) return;
         FillHist("CutFlow_EMu_1e2mu", 1., weight, 0., 20., 20);
       if(bjetColl.size()==0) return;
         FillHist("CutFlow_EMu_1e2mu", 2., weight, 0., 20., 20);
       if(jetColl.size()<2) return;
         FillHist("CutFlow_EMu_1e2mu", 3., weight, 0., 20., 20);
       if(Mmumu>40) return;
         FillHist("CutFlow_EMu_1e2mu", 4., weight, 0., 20., 20);
         if( fabs(Mmumu-15)<5 ) FillHist("Count_EMu_1e2mu_M15", 0., weight, 0., 1., 1);
         if( fabs(Mmumu-20)<5 ) FillHist("Count_EMu_1e2mu_M20", 0., weight, 0., 1., 1);
         if( fabs(Mmumu-25)<5 ) FillHist("Count_EMu_1e2mu_M25", 0., weight, 0., 1., 1);
         if( fabs(Mmumu-30)<5 ) FillHist("Count_EMu_1e2mu_M30", 0., weight, 0., 1., 1);
         if( fabs(Mmumu-35)<5 ) FillHist("Count_EMu_1e2mu_M35", 0., weight, 0., 1., 1);
     }
     if( (!Pass_TriggerEMu) &&  Pass_TriggerMuE ){
       if( !(electronLooseColl.size()==1 && muonLooseColl.size()==2) ) return;
       if( !(electronColl.size()==1      && muonColl.size()==2     ) ) return;
       if( !(muonColl.at(0).Charge()     != muonColl.at(1).Charge()) ) return;
       if( !(electronColl.at(0).Pt()>10  && muonColl.at(0).Pt()>25 && muonColl.at(1).Pt()>10) ) return;
       if( fabs(electronColl.at(0).Eta())>2.5 ) return;

       float Mmumu=(muonColl.at(0)+muonColl.at(1)).M();
       if(Mmumu<12) return;
         FillHist("CutFlow_MuE_1e2mu", 0., weight, 0., 20., 20);
       if(fabs(Mmumu-91.2)<10) return;
         FillHist("CutFlow_MuE_1e2mu", 1., weight, 0., 20., 20);
       if(bjetColl.size()==0) return;
         FillHist("CutFlow_MuE_1e2mu", 2., weight, 0., 20., 20);
       if(jetColl.size()<2) return;
         FillHist("CutFlow_MuE_1e2mu", 3., weight, 0., 20., 20);
       if(Mmumu>40) return;
         FillHist("CutFlow_MuE_1e2mu", 4., weight, 0., 20., 20);
         if( fabs(Mmumu-15)<5 ) FillHist("Count_MuE_1e2mu_M15", 0., weight, 0., 1., 1);
         if( fabs(Mmumu-20)<5 ) FillHist("Count_MuE_1e2mu_M20", 0., weight, 0., 1., 1);
         if( fabs(Mmumu-25)<5 ) FillHist("Count_MuE_1e2mu_M25", 0., weight, 0., 1., 1);
         if( fabs(Mmumu-30)<5 ) FillHist("Count_MuE_1e2mu_M30", 0., weight, 0., 1., 1);
         if( fabs(Mmumu-35)<5 ) FillHist("Count_MuE_1e2mu_M35", 0., weight, 0., 1., 1);
     }
   }
   if(CheckTrigAccept){

    int NEle=0, NMu=0, NEle_25=0, NMu_20=0, NMu_10=0;
    bool GenAccept=false, Gen1e2mu=false, Gen3mu=false, PassMmumu=true;
    std::vector<snu::KTruth> TrueMuColl;
    for(int i=2; i<(int) truthColl.size(); i++){
     
      int LeptonType=GetLeptonType(i,truthColl);
      if( !(LeptonType==1 || LeptonType==2 || LeptonType==3) ) continue;

      int fpid=fabs(truthColl.at(i).PdgId());
      if( !(fpid==11 || fpid==13) ) continue;

      float PT=truthColl.at(i).Pt(), fEta=fabs(truthColl.at(i).Eta());
      if(fpid==11){
        NEle++;
        FillHist("PTEl", truthColl.at(i).Pt(), weight, 0., 200., 40); 
        if(PT>25 && fEta<2.5) NEle_25++;
      }
      if(fpid==13){
        NMu++;
        FillHist("PTMu", truthColl.at(i).Pt(), weight, 0., 200., 40); 
        if(LeptonType==2) FillHist("PTMu_A", truthColl.at(i).Pt(), weight, 0., 200., 40);
        else if(LeptonType==1) FillHist("PTMu_W", truthColl.at(i).Pt(), weight, 0., 200., 40);
        if(PT>20 && fEta<2.4)  NMu_20++;
        if(PT>10 && fEta<2.4){ NMu_10++; TrueMuColl.push_back(truthColl.at(i));}
      }
    }//End of truthcoll
    FillHist("NEl" ,NEle, weight, 0., 10., 10);
    FillHist("NMu", NMu, weight, 0., 10., 10);

    if(NEle>0 && NMu>1){
      FillHist("AccEffFlow_1e2mu", 0., weight, 0., 10., 10);
      if(NEle_25>=1 && NMu_10>=1) FillHist("AccEffFlow_1e2mu", 1., weight, 0., 10., 10);
      if(NEle_25>=1 && NMu_10>=2) FillHist("AccEffFlow_1e2mu", 2., weight, 0., 10., 10);
      if(NEle_25>=1 && NMu_10>=2){ GenAccept=true; Gen1e2mu=true;}
    }
    if(NMu>2 && !Gen1e2mu){
      FillHist("AccEffFlow_3mu", 0., weight, 0., 10., 10);
      if(NMu_20>=1 && NMu_10>=2) FillHist("AccEffFlow_3mu", 1., weight, 0., 10., 10);
      if(NMu_20>=1 && NMu_10>=3 && NEle_25==0) FillHist("AccEffFlow_3mu", 2., weight, 0., 10., 10);
      if(NMu_20>=1 && NMu_10>=3 && NEle_25==0){ GenAccept=true; Gen3mu=true; }
    }


   //Now Reco
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, true);
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);


     std::vector<snu::KMuon> muonLooseColl, muonTightColl;
     for(int i=0; i<(int) muonPreColl.size(); i++){
       if(PassIDCriteria(muonPreColl.at(i),"POGTIsop6IPp2p1sig4NoChi"))   muonLooseColl.push_back(muonPreColl.at(i));
       if(PassIDCriteria(muonPreColl.at(i),"POGTIsop20IPp01p05sig4Chi4")) muonTightColl.push_back(muonPreColl.at(i));
     }

     std::vector<snu::KElectron> electronLooseColl, electronTightColl;
     for(int i=0; i<(int) electronPreColl.size(); i++){
       if(PassIDCriteria(electronPreColl.at(i), "HctoWAFakeLoose"))           electronLooseColl.push_back(electronPreColl.at(i));
       if(PassIDCriteria(electronPreColl.at(i), "POGWP90Isop06IPp025p1sig4")) electronTightColl.push_back(electronPreColl.at(i));
     }

     bool LeptonVeto=true;
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonLooseColl, electronLooseColl);
     std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");

     float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1., btag_sf=1.;
     float trigger_sf=1.;
     bool EventCand=false;
     int NLepTight=electronTightColl.size()+muonTightColl.size();
     if(NLepTight>=3 && muonTightColl.size()>=2) EventCand=true;
     if(EventCand){
       id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronTightColl);
       reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronTightColl);
       id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonTightColl);
       trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonTightColl);

       if(electronTightColl.size()>0){
         float trigger_sf1 = mcdata_correction->GetTriggerSF(electronTightColl, muonTightColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
         float trigger_sf2 = mcdata_correction->GetTriggerSF(electronTightColl, muonTightColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
         Pass_TriggerBG = PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
         Pass_TriggerH  = PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"); 
         if(Pass_TriggerBG || Pass_TriggerH) Pass_Trigger=true;
         trigger_sf    = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
       }
       else{
         float trigger_sf1 = mcdata_correction->GetTriggerSF(electronTightColl, muonTightColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v");
         float trigger_sf2 = mcdata_correction->GetTriggerSF(electronTightColl, muonTightColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v");
         Pass_TriggerBG = PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
         Pass_TriggerH  = PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
         if(Pass_TriggerBG || Pass_TriggerH) Pass_Trigger=true;

         trigger_sf    = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
       }
       btag_sf       = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);
     }
     weight*=id_weight_ele*reco_weight_ele*id_weight_mu*trk_weight_mu*trigger_sf*btag_sf;


     //4LepOpt
     if(electronTightColl.size()>0){
       bool Pass=true;
       if( !(electronTightColl.size()>=1 && muonTightColl.size()>=2) ) Pass=false;
       //if( !(electronTightColl.size()==1 && muonTightColl.size()==2 && 
       //      electronLooseColl.size()==1 && muonLooseColl.size()==2) ) Pass=false;
       if( Pass && !(electronTightColl.at(0).Pt()>25 && muonTightColl.at(1).Pt()>10) ) Pass=false;
       if(Pass) FillHist("AccEffFlow_1e2mu", 3., weight, 0., 10., 10);
       if( Pass && !Pass_Trigger ) Pass=false;
       if(Pass) FillHist("AccEffFlow_1e2mu", 4., weight, 0., 10., 10);
         std::vector<snu::KElectron> ElepColl, ElemColl;
         for(int i=0; i<electronTightColl.size(); i++){
           if(electronTightColl.at(i).Charge()>0) ElepColl.push_back(electronTightColl.at(i));
           else                                   ElemColl.push_back(electronTightColl.at(i));
         }
         std::vector<snu::KMuon> MupColl, MumColl;
         for(int i=0; i<muonTightColl.size(); i++){
           if(muonTightColl.at(i).Charge()>0) MupColl.push_back(muonTightColl.at(i));
           else                               MumColl.push_back(muonTightColl.at(i));
         }
         bool HasOS=true, OffZQCD=true;
         if(MupColl.size()==0 || MumColl.size()==0) HasOS=false;
         for(int i=0; i<(int) MupColl.size(); i++){
           for(int j=0; j<(int) MumColl.size(); j++){
             float Mmumu=(MupColl.at(i)+MumColl.at(j)).M();
             if( Mmumu<12 || fabs(Mmumu-91.2)<10 ) OffZQCD=false;
           }
         }
         for(int i=0; i<(int) ElepColl.size(); i++){
           for(int j=0; j<(int) ElemColl.size(); j++){
             float Melel=(ElepColl.at(i)+ElemColl.at(j)).M();
             if( Melel<12 || fabs(Melel-91.2)<10 ) OffZQCD=false;
           }
         }
       if( Pass && !HasOS                                            ) Pass=false;
       if(Pass) FillHist("AccEffFlow_1e2mu", 5., weight, 0., 10., 10);
       if( Pass && !OffZQCD                                          ) Pass=false;
       if(Pass) FillHist("AccEffFlow_1e2mu", 6., weight, 0., 10., 10);
       if( Pass && !(jetColl.size()>=2)                              ) Pass=false;
       if(Pass) FillHist("AccEffFlow_1e2mu", 7., weight, 0., 10., 10);
       if( Pass && bjetColl.size()==0                                ) Pass=false;
       if(Pass) FillHist("AccEffFlow_1e2mu", 8., weight, 0., 10., 10);
         float MA=GetHiggsMass(k_sample_name, "A");
         float Mmumu=Pass? (MupColl.at(MupColl.size()-1)+MumColl.at(MumColl.size()-1)).M():0.;
       if( Pass && fabs(Mmumu-MA)>5                                  ) Pass=false;
       if(Pass) FillHist("AccEffFlow_1e2mu", 9., weight, 0., 10., 10);
     }
     else{
       bool Pass=true;
       if( !(muonTightColl.size()>=3 && electronTightColl.size()==0) ) Pass=false;
       //if( !(muonTightColl.size()==3 && muonLooseColl.size()==3 
       //      && electronLooseColl.size()==0) ) Pass=false;       
       if( Pass && !(muonTightColl.at(0).Pt()>20)                    ) Pass=false;
       if(Pass) FillHist("AccEffFlow_3mu", 3., weight, 0., 10., 10);
       if( Pass && !Pass_Trigger                                     ) Pass=false;
       if(Pass) FillHist("AccEffFlow_3mu", 4., weight, 0., 10., 10);
         std::vector<snu::KMuon> MupColl, MumColl;
         for(int i=0; i<muonTightColl.size(); i++){
           if(muonTightColl.at(i).Charge()>0) MupColl.push_back(muonTightColl.at(i));
           else                               MumColl.push_back(muonTightColl.at(i));
         }
         bool HasOS=true, OffZQCD=true;
         if(MupColl.size()==0 || MumColl.size()==0) HasOS=false;
         for(int i=0; i<(int) MupColl.size(); i++){
           for(int j=0; j<(int) MumColl.size(); j++){
             float Mmumu=(MupColl.at(i)+MumColl.at(j)).M();
             if( Mmumu<12 || fabs(Mmumu-91.2)<10 ) OffZQCD=false;
           }
         }
       if( Pass && !HasOS                                            ) Pass=false;
       if(Pass) FillHist("AccEffFlow_3mu", 5., weight, 0., 10., 10);
       if( Pass && !OffZQCD                                          ) Pass=false;
       if(Pass) FillHist("AccEffFlow_3mu", 6., weight, 0., 10., 10);
         bool DoOptimization=true;
         if( Pass && DoOptimization ){
           for(int i=0; i<(int) MupColl.size(); i++){
             for(int j=0; j<(int) MumColl.size(); j++){
               if(GetLeptonType(MupColl.at(i),truthColl)==2 && GetLeptonType(MumColl.at(j),truthColl)==2){
                 if(muonTightColl.size()==3){ if(MupColl.size()==2) FillHist("PairTruth3_3mu", i+0.0001, weight, -1., 2., 3); 
                                              else                  FillHist("PairTruth3_3mu", j+0.0001, weight, -1., 2., 3);
                                              FillHist("PairTruth3_3mu", -1., weight, -1., 2., 3);
                                            }
                 if(muonTightColl.size()==4){ FillHist("PairTruth4_3mu", i*2.+j+0.0001,weight, -1., 4., 5);
                                              FillHist("PairTruth4_3mu", -1., weight, -1., 4., 5);
                                            }
               }
             }
           }
         }//End of Optimization Check
       if( Pass && !(jetColl.size()>=2)                              ) Pass=false;
       if(Pass) FillHist("AccEffFlow_3mu", 7., weight, 0., 10., 10);
       if( Pass && bjetColl.size()==0                                ) Pass=false;
       if(Pass) FillHist("AccEffFlow_3mu", 8., weight, 0., 10., 10);
         float MA=GetHiggsMass(k_sample_name, "A");
         float Mmumu=Pass? (MupColl.at(MupColl.size()-1)+MumColl.at(MumColl.size()-1)).M():0.;
       if( Pass && fabs(Mmumu-MA)>5                                  ) Pass=false;
       if(Pass) FillHist("AccEffFlow_3mu", 9., weight, 0., 10., 10);
     }


  }

/////////////////////////////////////////////////////////////////////////////////// 

return;

}// End of execute event loop
  


void Feb2018_TrigAccept::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Feb2018_TrigAccept::BeginCycle() throw( LQError ){
  
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

Feb2018_TrigAccept::~Feb2018_TrigAccept() {
  
  Message("In Feb2018_TrigAccept Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}




void Feb2018_TrigAccept::CheckDiMuValidation(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  if( !(MuTColl.size()==2 && MuLColl.size()==2 && EleTColl.size()==0 && EleLColl.size()==0) ) return;
  if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10) ) return;
  float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();

  if(Mmumu<50) return;

  FillHist("PTMu1"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
  FillHist("PTMu2"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
  FillHist("EtaMu1"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
  FillHist("EtaMu2"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
  FillHist("Nj"+Label, JetColl.size(), weight, 0., 10., 10);
  if(JetColl.size()>0){
    FillHist("PTj1"+Label, JetColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("Etaj1"+Label, JetColl.at(0).Eta(), weight, -5., 5., 20);
  }
  FillHist("Mmumu"+Label, (MuTColl.at(0)+MuTColl.at(1)).M(), weight, 60., 120., 60);

  float UnUsedSum= Option.Contains("NO")? BJetColl.at(0).Pt()+MET+METx+METy:0.;
  

}


void Feb2018_TrigAccept::CheckEMuValidation(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METphi, float weight, TString Label, TString Option){

  if( !(MuTColl.size()==1 && MuLColl.size()==1 && EleTColl.size()==1 && EleLColl.size()==1) ) return;
  if( !(MuTColl.at(0).Charge()!=EleTColl.at(0).Charge()) ) return;
  if( !(MuTColl.at(0).Pt()>10 && EleTColl.at(0).Pt()>25 && fabs(EleTColl.at(0).Eta())<2.5)  ) return;
  if( !(MuTColl.at(0).DeltaR(EleTColl.at(0))>0.4) ) return;
  if( JetColl.size()<2 ) return;

  FillHist("PTMu1"  +Label, MuTColl.at(0).Pt(),   weight, 0., 300., 30);
  FillHist("EtaMu1" +Label, MuTColl.at(0).Eta(),  weight, -5., 5., 20);
  FillHist("PTEle1" +Label, EleTColl.at(0).Pt(),  weight, 0., 300., 30);
  FillHist("EtaEle1"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
  FillHist("Nj"     +Label, JetColl.size(), weight, 0., 10., 10);
  FillHist("Nb"     +Label, BJetColl.size(), weight, 0., 5., 5);
  FillHist("PTj1"   +Label, JetColl.at(0).Pt(), weight, 0., 500., 50);
  FillHist("PTj2"   +Label, JetColl.at(1).Pt(), weight, 0., 400., 40);
  FillHist("Etaj1"  +Label, JetColl.at(0).Eta(), weight, -5., 5., 20);
  FillHist("Etaj2"  +Label, JetColl.at(1).Eta(), weight, -5., 5., 20);
  FillHist("MET"    +Label, MET, weight, 0., 400., 40);
  FillHist("METphi" +Label, METphi, weight, -3.15, 3.15, 10); 
   
}

void Feb2018_TrigAccept::DoSystRun(TString Cycle, TString Mode, std::vector<snu::KElectron> EleColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight){

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
    //CheckSRYield(MuColl, MuLColl, EleColl, EleLColl, JetColl,  BJetColl, MET, METx, METy, weight, SystDirLabel+SystKindLabel, ChannelLabel);
  }//End of Closure Cycle

  return;
}



float Feb2018_TrigAccept::ConeCorrectedPT(snu::KMuon Mu, float TightIsoCut){

  //float PTCorr=Mu.Pt();
  //float PTCorr=Mu.Pt()*(1.+RochIso(Mu,"0.4"));
  float PTCorr=Mu.Pt()*(1.+max((float) 0.,RochIso(Mu,"0.4")-TightIsoCut));
  return PTCorr;

}



float Feb2018_TrigAccept::ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut){

  //float PTCorr=Ele.Pt()*(1+Ele.PFRelIso(0.3));
  float PTCorr=Ele.Pt()*(1.+max(0.,Ele.PFRelIso(0.3)-TightIsoCut));
  return PTCorr;

}


float Feb2018_TrigAccept::FakeRateData(snu::KMuon Mu, TString Option){

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



float Feb2018_TrigAccept::FakeRateData(snu::KElectron Ele, TString Option){

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


//float Feb2018_TrigAccept::OptimiseSelection(std::vector<snu::KMuon> MuPreColl, std::vector<snu::KElectron> ElePreColl, std::vector<snu::KJet> JetPreNoVetoColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, TString Label){
//
//  float Mmumu=
//  if(
//
//}

float Feb2018_TrigAccept::GetFakeWeight(std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleLColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, TString Option){

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


  for(int i=0; i<(int) MuLColl.size(); i++){
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
  for(int i=0; i<(int) EleLColl.size(); i++){
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



void Feb2018_TrigAccept::FillCutFlow(TString cut, float weight){
  
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



void Feb2018_TrigAccept::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Feb2018_TrigAccept::MakeHistograms(){
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
  // *  Remove//Overide this Feb2018_TrigAcceptCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Feb2018_TrigAccept::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
