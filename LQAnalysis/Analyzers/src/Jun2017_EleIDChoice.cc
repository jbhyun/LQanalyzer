// $Id: Jun2017_EleIDChoice.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQJun2017_EleIDChoice Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Jun2017_EleIDChoice.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Jun2017_EleIDChoice);

 Jun2017_EleIDChoice::Jun2017_EleIDChoice() : AnalyzerCore(), out_muons(0) {

   SetLogName("Jun2017_EleIDChoice");
   Message("In Jun2017_EleIDChoice constructor", INFO);
   InitialiseAnalysis();
 }


 void Jun2017_EleIDChoice::InitialiseAnalysis() throw( LQError ) {
   
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

void Jun2017_EleIDChoice::ExecuteEvents()throw( LQError ){

////Basic Infos///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Apply the gen weight 
   if(!isData) weight*=MCweight;
   FillHist("GenWeight", MCweight, 1, -2., 2., 4);

   /// Acts on data to remove bad reconstructed event 
   //if(isData&& (! eventbase->GetEvent().LumiMask(lumimask))) return;

   m_logger<<DEBUG<<"RunNumber/Event Number = "<<eventbase->GetEvent().RunNumber()<<" : "<<eventbase->GetEvent().EventNumber()<<LQLogger::endmsg;
   m_logger<<DEBUG<<"isData = "<<isData<<LQLogger::endmsg;


   //Numbet of Vertex NoPUreweight
   FillHist("Nvtx_nocut_NoRW", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);

   //Pileup Reweight
   float pileup_reweight=1.;
   if (!k_isdata) { pileup_reweight=eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central); }

   //Numbet of Vertex PUreweight
   FillHist("Nvtx_nocut_PURW", eventbase->GetEvent().nVertices(), weight*pileup_reweight, 0., 50., 50);


   //Total Event(MCweight+PUrw)////////////////
   FillCutFlow("NoCut", weight*pileup_reweight);


   bool EleFake=false, BasicCompositionStudy=false, SigRegFake=false, EMuMu=false, TriMu=false, MatchingTest=false;
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("EleFake"))    EleFake=true;
     else if(k_flags.at(i).Contains("BasicCompositionStudy")) BasicCompositionStudy=true;
     else if(k_flags.at(i).Contains("SigRegFake")) SigRegFake=true;
     else if(k_flags.at(i).Contains("EMuMu"))      EMuMu=true;
     else if(k_flags.at(i).Contains("TriMu"))      TriMu=true;
     else if(k_flags.at(i).Contains("MatchingTest")) MatchingTest=true;
   }

    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //ListTriggersAvailable();

   //Normalisation Lumi Setting
   float trigger_ps_weight=1.;
   if(!isData) trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
   weight*=trigger_ps_weight;
   FillHist("TriggerPSWeight", trigger_ps_weight, 1., 0., 1., 100);


   //Trigger Cut
   //if(!Pass_Trigger) return;
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
   //std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, false);
    eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
   //  if     (MuEff ) { if( !(muonPreColl.size()>=1)     ) return; }
   //  else if(EleEff) { if( !(electronPreColl.size()>=1) ) return; }
     if     (SigRegFake && EMuMu){ if( !(muonPreColl.size()>=2 && electronPreColl.size()>=1) ) return; }
     else if(SigRegFake && TriMu){ if( !(muonPreColl.size()>=3) ) return; }
   /**********************************************************************************************************/

   //For Fake Study
   //Muon ID's to Test
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
   std::vector<snu::KMuon> muonVetoColl; eventbase->GetMuonSel()->Selection(muonVetoColl, true);


   //Electron ID's to Test
   //
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.5, 0.5);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.1);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetdxySigMax(3.);
   std::vector<snu::KElectron> electronVetoColl; eventbase->GetElectronSel()->Selection(electronVetoColl);

     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.5, 0.5);
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(3.);
   std::vector<snu::KElectron> electronFakeLColl; eventbase->GetElectronSel()->Selection(electronFakeLColl);

     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.1);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetdxySigMax(3.);
   std::vector<snu::KElectron> electronPOGTIPColl; eventbase->GetElectronSel()->Selection(electronPOGTIPColl);
//     eventbase->GetElectronSel()->SetCheckCharge(true);


   //------------------------------//
   //For SignalRegion Fake Study   //
   //------------------------------//
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
     eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetBSdz(0.1);
   std::vector<snu::KMuon> muonHN2FakeLColl; eventbase->GetMuonSel()->Selection(muonHN2FakeLColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.1);
     eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetBSdz(0.1);
   std::vector<snu::KMuon> muonHN2FakeTColl; eventbase->GetMuonSel()->Selection(muonHN2FakeTColl, true);


     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MEDIUM);
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0994, 0.107);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.2);//Not in ID, but additional safe WP
   std::vector<snu::KElectron> electronPOGMColl; eventbase->GetElectronSel()->Selection(electronPOGMColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.1);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetCheckCharge(true);
   std::vector<snu::KElectron> electronPOGTIPQColl; eventbase->GetElectronSel()->Selection(electronPOGTIPQColl);



     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronMVAMColl; eventbase->GetElectronSel()->Selection(electronMVAMColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP80);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronMVATColl; eventbase->GetElectronSel()->Selection(electronMVATColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);    eventbase->GetElectronSel()->SetCheckCharge(true);
   std::vector<snu::KElectron> electronMVAMQColl; eventbase->GetElectronSel()->Selection(electronMVAMQColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP80);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);    eventbase->GetElectronSel()->SetCheckCharge(true);
   std::vector<snu::KElectron> electronMVATQColl; eventbase->GetElectronSel()->Selection(electronMVATQColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_HN_MVA_TIGHT);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.05, 0.05);
     eventbase->GetElectronSel()->SetdxyBEMax(0.01, 0.01);   eventbase->GetElectronSel()->SetdzBEMax(0.04, 0.04);
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetCheckCharge(true);      eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronHNMVATColl; eventbase->GetElectronSel()->Selection(electronHNMVATColl);



   std::vector<snu::KElectron> electronNull;

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
   //  //eventbase->GetJetSel()->SetPileUpJetID(true,"Loose");
     bool LeptonVeto=true;
     //bool LeptonVeto=false;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonHN2FakeLColl, electronVetoColl);


   //Method to apply 2a method SF
   //std::vector<int> bIdxColl  = GetSFBJetIdx(jetColl,"Medium");
   //std::vector<int> ljIdxColl = GetSFLJetIdx(jetColl, bIdxColl, "Medium");
   //std::vector<snu::KJet> bjetColl2a; for(int i=0; i<bIdxColl.size(); i++) {bjetColl2a.push_back(jetColl.at(bIdxColl.at(i)));}
   //std::vector<snu::KJet> ljetColl2a; for(int i=0; i<ljIdxColl.size(); i++){ljetColl2a.push_back(jetColl.at(ljIdxColl.at(i)));}

   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");
   //std::vector<snu::KJet> ljetColl = SelLightJets(jetColl, "Medium");

   double met = eventbase->GetEvent().PFMETType1();
   double met_x = eventbase->GetEvent().PFMETType1x();
   double met_y = eventbase->GetEvent().PFMETType1y();
   //double met = eventbase->GetEvent().MET();
   double Pzv,Pzv1, Pzv2;
   snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));
   //snu::KParticle v[4]; v[0].SetPx(met_x); v[0].SetPy(met_y);
   //                     v[1].SetPx(met_x); v[1].SetPy(met_y);

   //int nbjets=bjetColl.size(); int njets=jetColl.size(); const int nljets=ljetColl.size();//njets-nbjets;//number of light jets
   int Nvtx=eventbase->GetEvent().nVertices();


   /*****************************************************
   **Scale Factors
   *****************************************************/
   float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1., btag_sf=1.;
   float trigger_sf=1.;
   float fake_weight=1.; bool EventCand=false;
   float geneff_weight=1., gennorm_weight=1.;

   /*This part is for boosting up speed.. SF part takes rather longer time than expected*/
//   if     (SigRegFake && TriMu) { if(muonHN2FakeLColl.size()==3)     EventCand=true; }
//   else if(SigRegFake && EMuMu) { if(muonHN2FakeLColl.size()==2 && electronLooseColl.size()==1) EventCand=true; }
//   else if(DiMuon_analysis){ if(muonLooseColl.size()==2)     EventCand=true; }
//   else if(DiEle_analysis) { if(electronLooseColl.size()==2) EventCand=true; }


   if(true){
   //if(SigRegFake){
     if(!isData){
      //trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );
  
      //id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronColl);
      //reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronFakeLColl);
  
      //id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muonColl);
      //iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
      //trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
  
      //btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

       geneff_weight   = GenFilterEfficiency(k_sample_name);
       gennorm_weight  = SignalNorm(k_sample_name, 200.);
     }
     else{
       //Perfect Prompt Ratio Approximation applied.(p=1), Applicable to generic number, combination of leptons under premise of high prompt rate.
       //if(k_running_nonprompt){
       //  fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", muonLooseColl.size(), electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", electronLooseColl.size());
       //}
     }
   }

   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight*fake_weight*btag_sf*geneff_weight*gennorm_weight;
   /***************************************************************************************************/

   //////Basic Objects Check//////////////////////////////////////////////////////////
/*   FillHist("Basic_Nj_orig", njets, weight, 0., 10., 10);
   FillHist("Basic_Nb_orig", nbjets, weight, 0., 10., 10);//Bjet def CVSincV2>0.89;pfKJet::CSVv2IVF_discriminator 
   FillHist("Basic_NmuT_orig", muonColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NmuL_orig", muonLooseColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeT_orig", electronColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeL_orig", electronLooseColl.size(), weight, 0., 10., 10);
   if(electronColl.size()==1) FillHist("Basic_NmuT_weT", muonColl.size(), weight, 0., 10., 10);
   if(TriMu_analysis){
     if(muonColl.size()>0) FillHist("Basic_Ptmu1_3mu_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>1) FillHist("Basic_Ptmu2_3mu_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>2) FillHist("Basic_Ptmu3_3mu_orig", muonColl.at(2).Pt(), weight, 0., 200., 200);
   }
   if(EMuMu_analysis){
     if(muonColl.size()>0) FillHist("Basic_Ptmu1_1e2mu_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>1) FillHist("Basic_Ptmu2_1e2mu_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
     if(electronColl.size()>0) FillHist("Basic_Pte_1e2mu_orig", electronColl.at(0).Pt(), weight, 0., 200., 200);
   }
   if(njets>0)  FillHist("Basic_j1_Et_orig", jetColl.at(0).Et(), weight, 0., 200., 200);
   if(njets>1)  FillHist("Basic_j2_Et_orig", jetColl.at(1).Et(), weight, 0., 200., 200);
   if(njets>2)  FillHist("Basic_j3_Et_orig", jetColl.at(2).Et(), weight, 0., 200., 200);
   if(njets>3)  FillHist("Basic_j4_Et_orig", jetColl.at(3).Et(), weight, 0., 200., 200);
   if(nbjets>0) FillHist("Basic_b1_Et_orig", bjetColl.at(0).Et(), weight, 0, 200., 200);
   if(nbjets>1) FillHist("Basic_b2_Et_orig", bjetColl.at(1).Et(), weight, 0., 200., 200);
   FillHist("Basic_METdist_orig", met, weight, 0., 200., 100);*/
   ///////////////////////////////////////////////////////////////////////////////////


/************************************************************************************/
//////////////////////////////////////////////////////////////////////////////////////
/////// MAIN ANALYSIS CODE ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
/************************************************************************************/

   if(SigRegFake){

     int NPr=NLeptonicBosonDecay(truthColl);
     if(k_sample_name.Contains("TT_powheg") && NPr>=1 ) return;


     //How Many Fake Events of TT DY pass analysis selection?
     //Which ID gives best significance?
     //And fake composition after that selection?
   
     int N_ElePOGCBM_MuHN2FakeT  = StepPassed(muonHN2FakeTColl, electronPOGMColl, jetColl, bjetColl, met, "EMuMu");
     int N_ElePOGCBT_MuHN2FakeT  = StepPassed(muonHN2FakeTColl, electronPOGTIPColl, jetColl, bjetColl, met, "EMuMu");
     int N_ElePOGCBTQ_MuHN2FakeT  = StepPassed(muonHN2FakeTColl, electronPOGTIPQColl, jetColl, bjetColl, met, "EMuMu");
     int N_ElePOGMVAM_MuHN2FakeT  = StepPassed(muonHN2FakeTColl, electronMVAMColl, jetColl, bjetColl, met, "EMuMu");
     int N_ElePOGMVAMQ_MuHN2FakeT = StepPassed(muonHN2FakeTColl, electronMVAMQColl, jetColl, bjetColl, met, "EMuMu");
     int N_ElePOGMVAT_MuHN2FakeT  = StepPassed(muonHN2FakeTColl, electronMVATColl, jetColl, bjetColl, met, "EMuMu");
     int N_ElePOGMVATQ_MuHN2FakeT = StepPassed(muonHN2FakeTColl, electronMVATQColl, jetColl, bjetColl, met, "EMuMu");
     int N_EleHNMVAT_MuHN2FakeT   = StepPassed(muonHN2FakeTColl, electronHNMVATColl, jetColl, bjetColl, met, "EMuMu");



     //Step6 : 3lep+>=1b+>=3j+Mmumu12-40Cut
     if( N_ElePOGCBM_MuHN2FakeT   >=6 ) FillHist("Yield_Step6", 0., weight, 0., 10., 10);
     if( N_ElePOGCBT_MuHN2FakeT   >=6 ) FillHist("Yield_Step6", 1., weight, 0., 10., 10);
     if( N_ElePOGCBTQ_MuHN2FakeT  >=6 ) FillHist("Yield_Step6", 2., weight, 0., 10., 10);
     if( N_ElePOGMVAM_MuHN2FakeT  >=6 ) FillHist("Yield_Step6", 3., weight, 0., 10., 10);
     if( N_ElePOGMVAMQ_MuHN2FakeT >=6 ) FillHist("Yield_Step6", 4., weight, 0., 10., 10);
     if( N_ElePOGMVAT_MuHN2FakeT  >=6 ) FillHist("Yield_Step6", 5., weight, 0., 10., 10);
     if( N_ElePOGMVATQ_MuHN2FakeT >=6 ) FillHist("Yield_Step6", 6., weight, 0., 10., 10);
     if( N_EleHNMVAT_MuHN2FakeT   >=6 ) FillHist("Yield_Step6", 7., weight, 0., 10., 10);

     //Step5 : 3lep(Mmumu>12)+>=1b+>=3j
     if( N_ElePOGCBM_MuHN2FakeT   >=5 ) FillHist("Yield_Step5", 0., weight, 0., 10., 10);
     if( N_ElePOGCBT_MuHN2FakeT   >=5 ) FillHist("Yield_Step5", 1., weight, 0., 10., 10);
     if( N_ElePOGCBTQ_MuHN2FakeT  >=5 ) FillHist("Yield_Step5", 2., weight, 0., 10., 10);
     if( N_ElePOGMVAM_MuHN2FakeT  >=5 ) FillHist("Yield_Step5", 3., weight, 0., 10., 10);
     if( N_ElePOGMVAMQ_MuHN2FakeT >=5 ) FillHist("Yield_Step5", 4., weight, 0., 10., 10);
     if( N_ElePOGMVAT_MuHN2FakeT  >=5 ) FillHist("Yield_Step5", 5., weight, 0., 10., 10);
     if( N_ElePOGMVATQ_MuHN2FakeT >=5 ) FillHist("Yield_Step5", 6., weight, 0., 10., 10);
     if( N_EleHNMVAT_MuHN2FakeT   >=5 ) FillHist("Yield_Step5", 7., weight, 0., 10., 10);

     //Step4 : 3lep(Mmumu>12)+>=1b
     if( N_ElePOGCBM_MuHN2FakeT   >=4 ) FillHist("Yield_Step4", 0., weight, 0., 10., 10);
     if( N_ElePOGCBT_MuHN2FakeT   >=4 ) FillHist("Yield_Step4", 1., weight, 0., 10., 10);
     if( N_ElePOGCBTQ_MuHN2FakeT  >=4 ) FillHist("Yield_Step4", 2., weight, 0., 10., 10);
     if( N_ElePOGMVAM_MuHN2FakeT  >=4 ) FillHist("Yield_Step4", 3., weight, 0., 10., 10);
     if( N_ElePOGMVAMQ_MuHN2FakeT >=4 ) FillHist("Yield_Step4", 4., weight, 0., 10., 10);
     if( N_ElePOGMVAT_MuHN2FakeT  >=4 ) FillHist("Yield_Step4", 5., weight, 0., 10., 10);
     if( N_ElePOGMVATQ_MuHN2FakeT >=4 ) FillHist("Yield_Step4", 6., weight, 0., 10., 10);
     if( N_EleHNMVAT_MuHN2FakeT   >=4 ) FillHist("Yield_Step4", 7., weight, 0., 10., 10);

     //Step3 : 3lep(Mmumu>12)
     if( N_ElePOGCBM_MuHN2FakeT   >=3 ) FillHist("Yield_Step3", 0., weight, 0., 10., 10);
     if( N_ElePOGCBT_MuHN2FakeT   >=3 ) FillHist("Yield_Step3", 1., weight, 0., 10., 10);
     if( N_ElePOGCBTQ_MuHN2FakeT  >=3 ) FillHist("Yield_Step3", 2., weight, 0., 10., 10);
     if( N_ElePOGMVAM_MuHN2FakeT  >=3 ) FillHist("Yield_Step3", 3., weight, 0., 10., 10);
     if( N_ElePOGMVAMQ_MuHN2FakeT >=3 ) FillHist("Yield_Step3", 4., weight, 0., 10., 10);
     if( N_ElePOGMVAT_MuHN2FakeT  >=3 ) FillHist("Yield_Step3", 5., weight, 0., 10., 10);
     if( N_ElePOGMVATQ_MuHN2FakeT >=3 ) FillHist("Yield_Step3", 6., weight, 0., 10., 10);
     if( N_EleHNMVAT_MuHN2FakeT   >=3 ) FillHist("Yield_Step3", 7., weight, 0., 10., 10);

     //Step2 : 3lep
     if( N_ElePOGCBM_MuHN2FakeT   >=2 ) FillHist("Yield_Step2", 0., weight, 0., 10., 10);
     if( N_ElePOGCBT_MuHN2FakeT   >=2 ) FillHist("Yield_Step2", 1., weight, 0., 10., 10);
     if( N_ElePOGCBTQ_MuHN2FakeT  >=2 ) FillHist("Yield_Step2", 2., weight, 0., 10., 10);
     if( N_ElePOGMVAM_MuHN2FakeT  >=2 ) FillHist("Yield_Step2", 3., weight, 0., 10., 10);
     if( N_ElePOGMVAMQ_MuHN2FakeT >=2 ) FillHist("Yield_Step2", 4., weight, 0., 10., 10);
     if( N_ElePOGMVAT_MuHN2FakeT  >=2 ) FillHist("Yield_Step2", 5., weight, 0., 10., 10);
     if( N_ElePOGMVATQ_MuHN2FakeT >=2 ) FillHist("Yield_Step2", 6., weight, 0., 10., 10);
     if( N_EleHNMVAT_MuHN2FakeT   >=2 ) FillHist("Yield_Step2", 7., weight, 0., 10., 10);


     
     //Fake Composition
     int NPrompt_Ele = NPromptFake_Ele(electronMVAMColl, truthColl, "Prompt");
     int NFake_Ele   = NPromptFake_Ele(electronMVAMColl, truthColl, "Fake"  );
     int NPrompt_Mu  = NPromptFake_Mu(muonHN2FakeTColl, truthColl, "Prompt");
     int NFake_Mu    = NPromptFake_Mu(muonHN2FakeTColl, truthColl, "Fake"  );

     if( N_ElePOGMVAM_MuHN2FakeT  >=6 ){
       if     (NPrompt_Ele==1 && NPrompt_Mu==2) FillHist("OriginComposition_Step6", 0., weight, 0., 10., 10);
       else if(NPrompt_Ele==1 && NPrompt_Mu==1) FillHist("OriginComposition_Step6", 1., weight, 0., 10., 10);
       else if(NPrompt_Ele==1 && NPrompt_Mu==0) FillHist("OriginComposition_Step6", 2., weight, 0., 10., 10);
       else if(NPrompt_Ele==0 && NPrompt_Mu==2) FillHist("OriginComposition_Step6", 3., weight, 0., 10., 10);
       else if(NPrompt_Ele==0 && NPrompt_Mu==1) FillHist("OriginComposition_Step6", 4., weight, 0., 10., 10);
       else if(NPrompt_Ele==0 && NPrompt_Mu==0) FillHist("OriginComposition_Step6", 5., weight, 0., 10., 10);
       else FillHist("OriginComposition_Step6", 9., weight, 0., 10., 10);

       if(NFake_Ele==0 && NFake_Mu==0) FillHist("NPrAndNFk_Step6", 0., weight, 0., 10., 10);
       else                            FillHist("NPrAndNFk_Step6", 1., weight, 0., 10., 10);
     } 
     if( N_ElePOGMVAM_MuHN2FakeT  >=5 ){
       if     (NPrompt_Ele==1 && NPrompt_Mu==2) FillHist("OriginComposition_Step5", 0., weight, 0., 10., 10);
       else if(NPrompt_Ele==1 && NPrompt_Mu==1) FillHist("OriginComposition_Step5", 1., weight, 0., 10., 10);
       else if(NPrompt_Ele==1 && NPrompt_Mu==0) FillHist("OriginComposition_Step5", 2., weight, 0., 10., 10);
       else if(NPrompt_Ele==0 && NPrompt_Mu==2) FillHist("OriginComposition_Step5", 3., weight, 0., 10., 10);
       else if(NPrompt_Ele==0 && NPrompt_Mu==1) FillHist("OriginComposition_Step5", 4., weight, 0., 10., 10);
       else if(NPrompt_Ele==0 && NPrompt_Mu==0) FillHist("OriginComposition_Step5", 5., weight, 0., 10., 10);
       else FillHist("OriginComposition_Step5", 9., weight, 0., 10., 10);

       if(NFake_Ele==0 && NFake_Mu==0) FillHist("NPrAndNFk_Step5", 0., weight, 0., 10., 10);
       else                            FillHist("NPrAndNFk_Step5", 1., weight, 0., 10., 10);
     } 
     if( N_ElePOGMVAM_MuHN2FakeT  >=4 ){
       if(NFake_Ele==0 && NFake_Mu==0) FillHist("NPrAndNFk_Step4", 0., weight, 0., 10., 10);
       else                            FillHist("NPrAndNFk_Step4", 1., weight, 0., 10., 10);
     } 
     if( N_ElePOGMVAM_MuHN2FakeT  >=3 ){
       if(NFake_Ele==0 && NFake_Mu==0) FillHist("NPrAndNFk_Step3", 0., weight, 0., 10., 10);
       else                            FillHist("NPrAndNFk_Step3", 1., weight, 0., 10., 10);
     } 



   }
   if(MatchingTest){


     for(int i=0; i<electronPreColl.size(); i++){
       int LepType=GetLeptonType(electronPreColl.at(i), truthColl);
       //if(fabs(LepType)>=6){
       //if(LepType==6 || LepType==-11 ){
       if(false){
         int MatchedIdx=GenMatchedIdx(electronPreColl.at(i),truthColl);
         cout<<"LepType : "<<LepType<<" MatchedIdx "<<MatchedIdx<<endl;
         PrintTruth();
       }
       FillHist("LepType", LepType, weight, -20., 20., 40);


       int MatchedIdx=GenMatchedIdx(electronPreColl.at(i),truthColl);
       bool FromHad=HasHadronicAncestor(MatchedIdx, truthColl);
       //if(fabs(LepType)>=4){
       //if( true ){
       if( (LepType<=-2 && !FromHad) || (LepType>=-1 && FromHad) ){
         cout<<"LepType : "<<LepType<<" MatchedIdx : "<<MatchedIdx<<" FromHad : "<<FromHad<<endl;
         PrintTruth();
       }
       if(LepType==0 && FromHad) FillHist("ErrorType", 0., weight, -20., 20., 40);
       if(LepType==1 && FromHad) FillHist("ErrorType", 1., weight, -20., 20., 40);
       if(LepType==2 && FromHad) FillHist("ErrorType", 2., weight, -20., 20., 40);
       if(LepType==3 && FromHad) FillHist("ErrorType", 3., weight, -20., 20., 40);
       if(LepType==4 && FromHad) FillHist("ErrorType", 4., weight, -20., 20., 40);
       if(LepType==5 && FromHad) FillHist("ErrorType", 5., weight, -20., 20., 40);

       if(LepType==-1 && FromHad)  FillHist("ErrorType", -1., weight, -20., 20., 40);
       if(LepType==-2 && !FromHad) FillHist("ErrorType", -2., weight, -20., 20., 40);
       if(LepType==-3 && !FromHad) FillHist("ErrorType", -3., weight, -20., 20., 40);
       if(LepType==-4 && !FromHad) FillHist("ErrorType", -4., weight, -20., 20., 40);
       if(LepType==-5 && !FromHad) FillHist("ErrorType", -5., weight, -20., 20., 40);
       if(LepType==-6 && !FromHad) FillHist("ErrorType", -6., weight, -20., 20., 40);
       if(LepType==-7 && !FromHad) FillHist("ErrorType", -7., weight, -20., 20., 40);
       if(LepType==-8 && !FromHad) FillHist("ErrorType", -8., weight, -20., 20., 40);
       if(LepType==-9 && !FromHad) FillHist("ErrorType", -9., weight, -20., 20., 40);
       if(LepType==-10 && !FromHad) FillHist("ErrorType", -10., weight, -20., 20., 40);
       if(LepType==-11 && !FromHad) FillHist("ErrorType", -11., weight, -20., 20., 40);
      

     }
     for(int i=0; i<muonPreColl.size(); i++){
       int LepType=GetLeptonType(muonPreColl.at(i), truthColl);
       //if( fabs(LepType)>=6 ){
       //if(LepType==6 || LepType==-11 ){
       if(false){
         int MatchedIdx=GenMatchedIdx(muonPreColl.at(i),truthColl);
         cout<<"LepType : "<<LepType<<" MatchedIdx "<<MatchedIdx<<endl;
         PrintTruth();
       }
       FillHist("LepType", LepType, weight, -20., 20., 40);


       int MatchedIdx=GenMatchedIdx(muonPreColl.at(i),truthColl);
       bool FromHad=HasHadronicAncestor(MatchedIdx, truthColl);
       //if(fabs(LepType)>=4){
       if( (LepType<=-2 && !FromHad) || (LepType>=-1 && FromHad) ){
       //if(true){
         cout<<"LepType : "<<LepType<<" MatchedIdx : "<<MatchedIdx<<" FromHad : "<<FromHad<<endl;
         PrintTruth();
       }

       if(LepType==0 && FromHad) FillHist("ErrorType", 0., weight, -20., 20., 40);
       if(LepType==1 && FromHad) FillHist("ErrorType", 1., weight, -20., 20., 40);
       if(LepType==2 && FromHad) FillHist("ErrorType", 2., weight, -20., 20., 40);
       if(LepType==3 && FromHad) FillHist("ErrorType", 3., weight, -20., 20., 40);
       if(LepType==4 && FromHad) FillHist("ErrorType", 4., weight, -20., 20., 40);
       if(LepType==5 && FromHad) FillHist("ErrorType", 5., weight, -20., 20., 40);

       if(LepType==-1 && FromHad)  FillHist("ErrorType", -1., weight, -20., 20., 40);
       if(LepType==-2 && !FromHad) FillHist("ErrorType", -2., weight, -20., 20., 40);
       if(LepType==-3 && !FromHad) FillHist("ErrorType", -3., weight, -20., 20., 40);
       if(LepType==-4 && !FromHad) FillHist("ErrorType", -4., weight, -20., 20., 40);
       if(LepType==-5 && !FromHad) FillHist("ErrorType", -5., weight, -20., 20., 40);
       if(LepType==-6 && !FromHad) FillHist("ErrorType", -6., weight, -20., 20., 40);
       if(LepType==-7 && !FromHad) FillHist("ErrorType", -7., weight, -20., 20., 40);
       if(LepType==-8 && !FromHad) FillHist("ErrorType", -8., weight, -20., 20., 40);
       if(LepType==-9 && !FromHad) FillHist("ErrorType", -9., weight, -20., 20., 40);
       if(LepType==-10 && !FromHad) FillHist("ErrorType", -10., weight, -20., 20., 40);
       if(LepType==-11 && !FromHad) FillHist("ErrorType", -11., weight, -20., 20., 40);

     }


     std::vector<snu::KElectron> FakeEleColl     = SkimLepColl(electronPreColl, truthColl, "Fake");
     std::vector<snu::KElectron> PromptEleColl   = SkimLepColl(electronPreColl, truthColl, "Prompt");
     std::vector<snu::KMuon>     FakeMuColl      = SkimLepColl(muonPreColl, truthColl, "Fake");
     std::vector<snu::KMuon>     PromptMuColl    = SkimLepColl(muonPreColl, truthColl, "Prompt");


     if(PromptEleColl.size()==2) FillHist("Mll_Prompt", (PromptEleColl.at(0)+PromptEleColl.at(1)).M(), weight, 60., 120., 60);
     if(PromptMuColl.size()==2) FillHist("Mll_Prompt", (PromptMuColl.at(0)+PromptMuColl.at(1)).M(), weight, 60., 120., 60);
     if(FakeEleColl.size()>=2) FillHist("Mll_Fake", (FakeEleColl.at(0)+FakeEleColl.at(1)).M(), weight, 60., 120., 60);
     if(FakeMuColl.size()>=2) FillHist("Mll_Fake", (FakeMuColl.at(0)+FakeMuColl.at(1)).M(), weight, 60., 120., 60);

     
   }

   //Dummy Area
   //NotMuch Interest since reliso is much more evidently different
   //FillHist("FakeElAbsIso_0l", FakeColl.at(i).PFAbsIso(0.3), weight, 0., 100., 100);
   //FillHist("FakeElChIso_0l", FakeColl.at(i).PFChargedHadronIso(0.3), weight, 0., 100., 100);
   //FillHist("FakeElPhoIso_0l", FakeColl.at(i).PFPhotonIso(0.3), weight, 0., 100., 100);
   //FillHist("FakeElNeuIso_0l", FakeColl.at(i).PFNeutralHadronIso(0.3), weight, 0., 100., 100);
   //FillHist("PromptElAbsIso_0l", PromptColl.at(i).PFAbsIso(0.3), weight, 0., 100., 100);
   //FillHist("PromptElChIso_0l", PromptColl.at(i).PFChargedHadronIso(0.3), weight, 0., 100., 100);
   //FillHist("PromptElPhoIso_0l", PromptColl.at(i).PFPhotonIso(0.3), weight, 0., 100., 100);
   //FillHist("PromptElNeuIso_0l", PromptColl.at(i).PFNeutralHadronIso(0.3), weight, 0., 100., 100);

//     if(NPromptGenLepAll>2){ //{PrintTruth(); cout<<"NPrL "<<NPromptGenLepAll<<endl;}
//       int NnewEle=0;
//       std::vector<int> PromptIdxColl;
//       std::vector<int> PromptLepTypeColl;
//       int NTauHard=0;
//       for(int i=2; i<truthColl.size(); i++){
//         if(truthColl.at(i).IndexMother()<2) continue;
//         if(truthColl.at(i).GenStatus()!=1)  continue;
//         int pid=truthColl.at(i).PdgId();
//         if( fabs(pid)==15 ){
//            int mpid=fabs(truthColl.at(truthColl.at(i).IndexMother()).PdgId());
//            if((mpid==23 || mpid==24) && truthColl.at(i).GenStatus()==2)              NTauHard++;
//            else if(truthColl.at(i).GenStatus()>20 && truthColl.at(i).GenStatus()<30) NTauHard++;
//         }
//         if( !(fabs(pid)==11 || fabs(pid)==13) ) continue;
//
//         int LepType=GetLeptonType(i, truthColl);
//         if(LepType==1 || LepType==3) {PromptIdxColl.push_back(i); PromptLepTypeColl.push_back(LepType);}
//
//       }
//
//       std::vector<int> ConvBranchIdxColl;
//       std::vector<int> ConvResultColl;
//       for(int i=0; i<PromptIdxColl.size(); i++){
//         for(int j=i+1; j<PromptIdxColl.size(); j++){
//           int pid1 =truthColl.at(PromptIdxColl.at(i)).PdgId();
//           int midx1=truthColl.at(PromptIdxColl.at(i)).IndexMother(), mid1=truthColl.at(midx1).PdgId();
//           while(mid1==pid1){
//             int pid2 =truthColl.at(PromptIdxColl.at(j)).PdgId();
//             int midx2=truthColl.at(PromptIdxColl.at(j)).IndexMother(), mid2=truthColl.at(midx2).PdgId();
//             while(mid2==pid2){
//               if(midx1==midx2) break;
//               else {midx2=truthColl.at(midx2).IndexMother(); mid2=truthColl.at(midx2).PdgId();}
//             }
//             if(midx1==midx2){ ConvBranchIdxColl.push_back(midx1); ConvResultColl.push_back(10*i+j); break; }
//             else {midx1=truthColl.at(midx1).IndexMother(); mid1=truthColl.at(midx1).PdgId();}
//           }
//         }
//       }
//       int NConv=ConvBranchIdxColl.size(); int NConvEE=0, NConvMuMu=0, NConvEMu=0;
//       int NMultiConv=0, SameBranchCount=0;
//       for(int i=0; i<ConvBranchIdxColl.size(); i++){
//         for(int j=i+1; j<ConvBranchIdxColl.size(); j++){
//           if(ConvBranchIdxColl.at(i)==ConvBranchIdxColl.at(j)) SameBranchCount++;
//         }
//       }
//       if(SameBranchCount==3) NMultiConv++;//Basically we should also count all the nC2 (n=3,4...) but higher order terms need not be considered till now yet.
//       NConv-=NMultiConv;
//       if(ConvResultColl.size()>0){
//         for(int i=0; i<ConvBranchIdxColl.size(); i++){
////           if(truthColl.at(ConvBranchIdxColl.at(i)).Pt()<0.1)      continue;
////           else if(truthColl.at(ConvResultColl.at(i)/10).Pt()<0.1) continue;
////           else if(truthColl.at(ConvResultColl.at(i)%10).Pt()<0.1) continue;
//
//           int PrIdx     =fabs(truthColl.at(ConvBranchIdxColl.at(i)).Pt()-truthColl.at(PromptIdxColl.at(ConvResultColl.at(i)/10)).Pt())<fabs(truthColl.at(ConvBranchIdxColl.at(i)).Pt()-truthColl.at(PromptIdxColl.at(ConvResultColl.at(i)%10)).Pt()) ? PromptIdxColl.at(ConvResultColl.at(i)/10) : PromptIdxColl.at(ConvResultColl.at(i)%10);
//           int ConvSonIdx=fabs(truthColl.at(ConvBranchIdxColl.at(i)).Pt()-truthColl.at(PromptIdxColl.at(ConvResultColl.at(i)/10)).Pt())<fabs(truthColl.at(ConvBranchIdxColl.at(i)).Pt()-truthColl.at(PromptIdxColl.at(ConvResultColl.at(i)%10)).Pt()) ? PromptIdxColl.at(ConvResultColl.at(i)%10) : PromptIdxColl.at(ConvResultColl.at(i)/10);
//
//           FillHist("Pt_PrCand", truthColl.at(PrIdx).Pt(), weight, 0., 200., 400);
//           FillHist("Pt_CvCand", truthColl.at(ConvSonIdx).Pt(), weight, 0., 200., 400);
//           FillHist("Pt_MotherCand", truthColl.at(ConvBranchIdxColl.at(i)).Pt(), weight, 0., 200., 400);
//           FillHist("dPt_MPr", fabs(truthColl.at(ConvBranchIdxColl.at(i)).Pt()-truthColl.at(PrIdx).Pt()), weight, 0., 200., 400);
//           FillHist("dPt_MCv", fabs(truthColl.at(ConvBranchIdxColl.at(i)).Pt()-truthColl.at(ConvSonIdx).Pt()), weight, 0., 200., 400);
//
//           if(truthColl.at(ConvBranchIdxColl.at(i)).Pt()>0 && truthColl.at(PrIdx).Pt()>0 && truthColl.at(ConvSonIdx).Pt()>0){
//
//           FillHist("dRPrMa", truthColl.at(ConvBranchIdxColl.at(i)).DeltaR(truthColl.at(PrIdx)), weight, 0., 1., 100);
//           FillHist("dRCvMa", truthColl.at(ConvBranchIdxColl.at(i)).DeltaR(truthColl.at(ConvSonIdx)), weight, 0., 1., 100);
//           FillHist("PtRelPrMa", truthColl.at(PrIdx).Pt()/truthColl.at(ConvBranchIdxColl.at(i)).Pt(), weight, 0., 2., 200);
//           FillHist("PtRelCvMa", truthColl.at(ConvSonIdx).Pt()/truthColl.at(ConvBranchIdxColl.at(i)).Pt(), weight, 0., 2., 200);
//           //}
//           }
//           if(ConvBranchIdxColl.at(i)!=truthColl.at(PromptIdxColl.at(ConvResultColl.at(i)/10)).IndexMother() || NMultiConv>0 ) PrintTruth();
//         }
//       }
//       FillHist("NPrL_New", NPromptGenLepAll-NConv, 1., -10., 10., 20);
//       FillHist("NTauHard", NTauHard, 1., 0., 10., 10);
//       //if((NPromptGenLepAll-NConv)!=2) PrintTruth();
//       //if(ConvBranchIdxColl.size()>0) PrintTruth();
//     }


 
/////////////////////////////////////////////////////////////////////////////////// 

return;
}// End of execute event loop
  


void Jun2017_EleIDChoice::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Jun2017_EleIDChoice::BeginCycle() throw( LQError ){
  
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

Jun2017_EleIDChoice::~Jun2017_EleIDChoice() {
  
  Message("In Jun2017_EleIDChoice Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}


int Jun2017_EleIDChoice::StepPassed(std::vector<snu::KMuon> MuColl, std::vector<snu::KElectron> EleColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, TString Option){

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

   return PassedSteps;
}


int Jun2017_EleIDChoice::NPromptFake_Mu(std::vector<snu::KMuon> MuColl, std::vector<snu::KTruth> truthColl, TString Option){

   int Nprompt=0, Nfake=0;

   bool ReturnPrompt=false, ReturnFake=false;
   bool DrawHist=false;
   if     (Option.Contains("Prompt")) ReturnPrompt=true;
   else if(Option.Contains("Fake"))   ReturnFake=true;

   for(int i=0; i<MuColl.size(); i++){
     int MuType=GetLeptonType(MuColl.at(i), truthColl);

     if     (MuType==1) Nprompt++;
     else if(MuType==2) Nprompt++;
     else if(MuType==3) Nprompt++;
     else if(MuType<0 ) Nfake++;
     else if(MuType>=4) Nfake++;
   }

   if     (ReturnPrompt) return Nprompt;
   else if(ReturnFake)   return Nfake;

   return 0.;
  
}

int Jun2017_EleIDChoice::NPromptFake_Ele(std::vector<snu::KElectron> EleColl, std::vector<snu::KTruth> truthColl, TString Option){

   int Nprompt=0, Nfake=0;

   bool ReturnPrompt=false, ReturnFake=false;
   bool DrawHist=false;
   if     (Option.Contains("Prompt")) ReturnPrompt=true;
   else if(Option.Contains("Fake"))   ReturnFake=true;

   for(int i=0; i<EleColl.size(); i++){
     int EleType=GetLeptonType(EleColl.at(i), truthColl);

     if     (EleType==1) Nprompt++;
     else if(EleType==2) Nprompt++;
     else if(EleType==3) Nprompt++;
     else if(EleType<0 ) Nfake++;
     else if(EleType>=4) Nfake++;
   }

   if     (ReturnPrompt) return Nprompt;
   else if(ReturnFake)   return Nfake;

   return 0.;
  
}



bool Jun2017_EleIDChoice::IsConvCand(snu::KElectron Ele, std::vector<snu::KTruth> TruthColl, TString Option){

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



void Jun2017_EleIDChoice::FillCutFlow(TString cut, float weight){
  
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



void Jun2017_EleIDChoice::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Jun2017_EleIDChoice::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //Basic Hists/////////////////////////////////////////////////
  //Data properties before selection  
  AnalyzerCore::MakeHistograms("Basic_b2_Et_orig", 200, 0., 200.);



  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Jun2017_EleIDChoiceCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Jun2017_EleIDChoice::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
