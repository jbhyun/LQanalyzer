// $Id: May2017_ObjectEff.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQMay2017_ObjectEff Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "May2017_ObjectEff.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (May2017_ObjectEff);

 May2017_ObjectEff::May2017_ObjectEff() : AnalyzerCore(), out_muons(0) {

   SetLogName("May2017_ObjectEff");
   Message("In May2017_ObjectEff constructor", INFO);
   InitialiseAnalysis();
 }


 void May2017_ObjectEff::InitialiseAnalysis() throw( LQError ) {
   
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

void May2017_ObjectEff::ExecuteEvents()throw( LQError ){

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


   bool EleEff=false, MuEff=false, ConversionEff=false, MatchingEff=false;
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("EleEff"))        EleEff=true;
     else if(k_flags.at(i).Contains("MuEff"))         MuEff =true;
     else if(k_flags.at(i).Contains("ConversionEff")) ConversionEff=true;
     else if(k_flags.at(i).Contains("MatchingEff"))   MatchingEff=true;
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
    eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
   //  if     (MuEff ) { if( !(muonPreColl.size()>=1)     ) return; }
   //  else if(EleEff) { if( !(electronPreColl.size()>=1) ) return; }
   /**********************************************************************************************************/

   //Muon ID's to Test
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
   std::vector<snu::KMuon> muonVetoColl; eventbase->GetMuonSel()->Selection(muonVetoColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.15);
   std::vector<snu::KMuon> muonPOGTColl; eventbase->GetMuonSel()->Selection(muonPOGTColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
     eventbase->GetMuonSel()->SetBSdxy(0.05);                eventbase->GetMuonSel()->SetdxySigMax(3.);
   std::vector<snu::KMuon> muonHN1FakeLColl; eventbase->GetMuonSel()->Selection(muonHN1FakeLColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.1);
     eventbase->GetMuonSel()->SetBSdxy(0.05);                eventbase->GetMuonSel()->SetdxySigMax(3.);
   std::vector<snu::KMuon> muonHN1FakeTColl; eventbase->GetMuonSel()->Selection(muonHN1FakeTColl, true);
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


   //Electron ID's to Test
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
   std::vector<snu::KElectron> electronPOGTColl; eventbase->GetElectronSel()->Selection(electronPOGTColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.1);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetdxySigMax(3.);
   std::vector<snu::KElectron> electronPOGTIPColl; eventbase->GetElectronSel()->Selection(electronPOGTIPColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.1);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronPOGTIPConvColl; eventbase->GetElectronSel()->Selection(electronPOGTIPConvColl);

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

   //  eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
   //  eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
   //  //eventbase->GetJetSel()->SetPileUpJetID(true,"Loose");
   //  bool LeptonVeto=true;
   //std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonLooseColl, electronPOGMColl);



   //Method to apply 2a method SF
   //std::vector<int> bIdxColl  = GetSFBJetIdx(jetColl,"Medium");
   //std::vector<int> ljIdxColl = GetSFLJetIdx(jetColl, bIdxColl, "Medium");
   //std::vector<snu::KJet> bjetColl2a; for(int i=0; i<bIdxColl.size(); i++) {bjetColl2a.push_back(jetColl.at(bIdxColl.at(i)));}
   //std::vector<snu::KJet> ljetColl2a; for(int i=0; i<ljIdxColl.size(); i++){ljetColl2a.push_back(jetColl.at(ljIdxColl.at(i)));}

   //std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");
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
//   if     (TriMu_analysis) { if(muonLooseColl.size()==3)     EventCand=true; }
//   else if(EMuMu_analysis) { if(muonLooseColl.size()==2 && electronLooseColl.size()==1) EventCand=true; }
//   else if(DiMuon_analysis){ if(muonLooseColl.size()==2)     EventCand=true; }
//   else if(DiEle_analysis) { if(electronLooseColl.size()==2) EventCand=true; }


   if(true){
     if(!isData){
      //trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );
  
      //id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronColl);
      //reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
  
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


   if(EleEff){

     //Np/f All - Denominator, 0:Prompt, 1:Fake
     int Nprompt_All     = NPromptFake_Ele(electronPreColl,truthColl,"Prompt");
     int Nfake_All       = NPromptFake_Ele(electronPreColl,truthColl,"Fake");

     FillHist("IDeff_NEle", 0., Nprompt_All*weight, 0., 2., 2);
     FillHist("IDeff_NEle", 1., Nfake_All*weight, 0., 2., 2);


     
     //N_IDpass - Numerator, 
     int Nprompt_POGM    = NPromptFake_Ele(electronPOGMColl,truthColl,"Prompt");
     int Nprompt_POGT    = NPromptFake_Ele(electronPOGTColl,truthColl,"Prompt");
     int Nprompt_POGTIP  = NPromptFake_Ele(electronPOGTIPColl,truthColl,"Prompt");
     int Nprompt_POGTIPQ = NPromptFake_Ele(electronPOGTIPQColl,truthColl,"Prompt");
     int Nprompt_MVAM    = NPromptFake_Ele(electronMVAMColl,truthColl,"Prompt");
     int Nprompt_MVAMQ   = NPromptFake_Ele(electronMVAMQColl,truthColl,"Prompt");
     int Nprompt_MVAT    = NPromptFake_Ele(electronMVATColl,truthColl,"Prompt");
     int Nprompt_MVATQ   = NPromptFake_Ele(electronMVATQColl,truthColl,"Prompt");
     int Nprompt_HNMVAT  = NPromptFake_Ele(electronHNMVATColl,truthColl,"Prompt");

     int Nfake_POGM    = NPromptFake_Ele(electronPOGMColl,truthColl,"Fake");
     int Nfake_POGT    = NPromptFake_Ele(electronPOGTColl,truthColl,"Fake");
     int Nfake_POGTIP  = NPromptFake_Ele(electronPOGTIPColl,truthColl,"Fake");
     int Nfake_POGTIPQ = NPromptFake_Ele(electronPOGTIPQColl,truthColl,"Fake");
     int Nfake_MVAM    = NPromptFake_Ele(electronMVAMColl,truthColl,"Fake");
     int Nfake_MVAMQ   = NPromptFake_Ele(electronMVAMQColl,truthColl,"Fake");
     int Nfake_MVAT    = NPromptFake_Ele(electronMVATColl,truthColl,"Fake");
     int Nfake_MVATQ   = NPromptFake_Ele(electronMVATQColl,truthColl,"Fake");
     int Nfake_HNMVAT  = NPromptFake_Ele(electronHNMVATColl,truthColl,"Fake");


     FillHist("IDeff_NEleIDp", 0., Nprompt_POGM*weight, 0., 10., 10);
     FillHist("IDeff_NEleIDp", 1., Nprompt_POGT*weight, 0., 10., 10);
     FillHist("IDeff_NEleIDp", 2., Nprompt_POGTIP*weight, 0., 10., 10);
     FillHist("IDeff_NEleIDp", 3., Nprompt_POGTIPQ*weight, 0., 10., 10);
     FillHist("IDeff_NEleIDp", 4., Nprompt_MVAM*weight, 0., 10., 10);
     FillHist("IDeff_NEleIDp", 5., Nprompt_MVAMQ*weight, 0., 10., 10);
     FillHist("IDeff_NEleIDp", 6., Nprompt_MVAT*weight, 0., 10., 10);
     FillHist("IDeff_NEleIDp", 7., Nprompt_MVATQ*weight, 0., 10., 10);
     FillHist("IDeff_NEleIDp", 8., Nprompt_HNMVAT*weight, 0., 10., 10);

     FillHist("IDeff_NEleIDf", 0., Nfake_POGM*weight, 0., 10., 10);
     FillHist("IDeff_NEleIDf", 1., Nfake_POGT*weight, 0., 10., 10);
     FillHist("IDeff_NEleIDf", 2., Nfake_POGTIP*weight, 0., 10., 10);
     FillHist("IDeff_NEleIDf", 3., Nfake_POGTIPQ*weight, 0., 10., 10);
     FillHist("IDeff_NEleIDf", 4., Nfake_MVAM*weight, 0., 10., 10);
     FillHist("IDeff_NEleIDf", 5., Nfake_MVAMQ*weight, 0., 10., 10);
     FillHist("IDeff_NEleIDf", 6., Nfake_MVAT*weight, 0., 10., 10);
     FillHist("IDeff_NEleIDf", 7., Nfake_MVATQ*weight, 0., 10., 10);
     FillHist("IDeff_NEleIDf", 8., Nfake_HNMVAT*weight, 0., 10., 10);


     for(int i=0; i<electronPreColl.size(); i++){
       int LeptonType=GetLeptonType(electronPreColl.at(i), truthColl);
       if(LeptonType==1 || LeptonType==3){
         FillHist("IDeff_NElep_PtEta", electronPreColl.at(i).Pt(), fabs(electronPreColl.at(i).Eta()), weight, 0., 200., 20, 0., 2.5, 5);
         if(PassIDCriteria(electronPreColl.at(i), "HNMVATIPQ")){
           FillHist("IDeff_NEleIDp_HN_PtEta", electronPreColl.at(i).Pt(), fabs(electronPreColl.at(i).Eta()), weight, 0., 200., 20, 0., 2.5, 5);
         }
         if(PassIDCriteria(electronPreColl.at(i), "POGMVAMIP")){
           FillHist("IDeff_NEleIDp_POGMVA_PtEta", electronPreColl.at(i).Pt(), fabs(electronPreColl.at(i).Eta()), weight, 0., 200., 20, 0., 2.5, 5);
         }
       }
       else if(LeptonType<0 || LeptonType>=4){
          FillHist("IDeff_NElef_PtEta", electronPreColl.at(i).Pt(), fabs(electronPreColl.at(i).Eta()), weight, 0., 200., 20, 0., 2.5, 5);
         if(PassIDCriteria(electronPreColl.at(i), "HNMVATIPQ")){
           FillHist("IDeff_NEleIDf_HN_PtEta", electronPreColl.at(i).Pt(), fabs(electronPreColl.at(i).Eta()), weight, 0., 200., 20, 0., 2.5, 5);
         }
         if(PassIDCriteria(electronPreColl.at(i), "POGMVAMIP")){
           FillHist("IDeff_NEleIDf_POGMVA_PtEta", electronPreColl.at(i).Pt(), fabs(electronPreColl.at(i).Eta()), weight, 0., 200., 20, 0., 2.5, 5);
         }
       }
     }
   }
   if(MuEff){

     //Np/f All - Denominator, 0:Prompt, 1:Fake
     int Nprompt_All     = NPromptFake_Mu(muonPreColl,truthColl,"BSM");
     int Nfake_All       = NPromptFake_Mu(muonPreColl,truthColl,"Fake");


     FillHist("IDeff_NMu", 0., Nprompt_All*weight, 0., 2., 2);
     FillHist("IDeff_NMu", 1., Nfake_All*weight, 0., 2., 2);

     
     //N_IDpass - Numerator, 
     int Nprompt_Veto     = NPromptFake_Mu(muonVetoColl,truthColl,"BSM");
     int Nprompt_POGT     = NPromptFake_Mu(muonPOGTColl,truthColl,"BSM");
     int Nprompt_HN1FakeL = NPromptFake_Mu(muonHN1FakeLColl,truthColl,"BSM");
     int Nprompt_HN1FakeT = NPromptFake_Mu(muonHN1FakeTColl,truthColl,"BSM");
     int Nprompt_HN2FakeL = NPromptFake_Mu(muonHN2FakeLColl,truthColl,"BSM");
     int Nprompt_HN2FakeT = NPromptFake_Mu(muonHN2FakeTColl,truthColl,"BSM");

     int Nfake_Veto     = NPromptFake_Mu(muonVetoColl,truthColl,"Fake");
     int Nfake_POGT     = NPromptFake_Mu(muonPOGTColl,truthColl,"Fake");
     int Nfake_HN1FakeL = NPromptFake_Mu(muonHN1FakeLColl,truthColl,"Fake");
     int Nfake_HN1FakeT = NPromptFake_Mu(muonHN1FakeTColl,truthColl,"Fake");
     int Nfake_HN2FakeL = NPromptFake_Mu(muonHN2FakeLColl,truthColl,"Fake");
     int Nfake_HN2FakeT = NPromptFake_Mu(muonHN2FakeTColl,truthColl,"Fake");


     FillHist("IDeff_NMuIDp", 0., Nprompt_Veto*weight, 0., 10., 10);
     FillHist("IDeff_NMuIDp", 1., Nprompt_POGT*weight, 0., 10., 10);
     FillHist("IDeff_NMuIDp", 2., Nprompt_HN1FakeL*weight, 0., 10., 10);
     FillHist("IDeff_NMuIDp", 3., Nprompt_HN1FakeT*weight, 0., 10., 10);
     FillHist("IDeff_NMuIDp", 4., Nprompt_HN2FakeL*weight, 0., 10., 10);
     FillHist("IDeff_NMuIDp", 5., Nprompt_HN2FakeT*weight, 0., 10., 10);

     FillHist("IDeff_NMuIDf", 0., Nfake_Veto*weight, 0., 10., 10);
     FillHist("IDeff_NMuIDf", 1., Nfake_POGT*weight, 0., 10., 10);
     FillHist("IDeff_NMuIDf", 2., Nfake_HN1FakeL*weight, 0., 10., 10);
     FillHist("IDeff_NMuIDf", 3., Nfake_HN1FakeT*weight, 0., 10., 10);
     FillHist("IDeff_NMuIDf", 4., Nfake_HN2FakeL*weight, 0., 10., 10);
     FillHist("IDeff_NMuIDf", 5., Nfake_HN2FakeT*weight, 0., 10., 10);

   }
   if(ConversionEff){

     if(k_sample_name.Contains("ZG")){

       //PrintTruth();
       //cout<<HardPhotonIdx<<endl;
  
       int Nconv_All        = 0;
       int Nconv_POGM       = 0;
       int Nconv_POGT       = 0;
       int Nconv_POGTIP     = 0;
       int Nconv_POGTIPConv = 0;
       int Nconv_POGTIPQ    = 0;
       int Nconv_MVAM       = 0;
       int Nconv_MVAMQ      = 0;
       int Nconv_MVAT       = 0;
       int Nconv_MVATQ      = 0;
       int Nconv_HNMVAT     = 0;
  
       for(int i=0; i<electronPreColl.size(); i++)       { if(IsConvCand(electronPreColl.at(i), truthColl))        Nconv_All++; }
       for(int i=0; i<electronPOGMColl.size(); i++)      { if(IsConvCand(electronPOGMColl.at(i), truthColl))       Nconv_POGM++; }
       for(int i=0; i<electronPOGTColl.size(); i++)      { if(IsConvCand(electronPOGTColl.at(i), truthColl))       Nconv_POGT++; }
       for(int i=0; i<electronPOGTIPColl.size(); i++)    { if(IsConvCand(electronPOGTIPColl.at(i), truthColl))     Nconv_POGTIP++; }
       for(int i=0; i<electronPOGTIPConvColl.size(); i++){ if(IsConvCand(electronPOGTIPConvColl.at(i), truthColl)) Nconv_POGTIPConv++; }
       for(int i=0; i<electronPOGTIPQColl.size(); i++)   { if(IsConvCand(electronPOGTIPQColl.at(i), truthColl))    Nconv_POGTIPQ++; }
       for(int i=0; i<electronMVAMColl.size(); i++)      { if(IsConvCand(electronMVAMColl.at(i), truthColl))       Nconv_MVAM++; }
       for(int i=0; i<electronMVAMQColl.size(); i++)     { if(IsConvCand(electronMVAMQColl.at(i), truthColl))      Nconv_MVAMQ++; }
       for(int i=0; i<electronMVATColl.size(); i++)      { if(IsConvCand(electronMVATColl.at(i), truthColl))       Nconv_MVAT++; }
       for(int i=0; i<electronMVATQColl.size(); i++)     { if(IsConvCand(electronMVAMQColl.at(i), truthColl))      Nconv_MVATQ++; }
       for(int i=0; i<electronHNMVATColl.size(); i++)    { if(IsConvCand(electronHNMVATColl.at(i), truthColl))     Nconv_HNMVAT++; }
  
       //N all conversion ele
       FillHist("IDeff_NEleConv", 0., Nconv_All*weight, 0., 1., 1);
  
       //After ID
       FillHist("IDeff_NEleIDConv", 0., Nconv_POGM*weight, 0., 10., 10);
       FillHist("IDeff_NEleIDConv", 1., Nconv_POGT*weight, 0., 10., 10);
       FillHist("IDeff_NEleIDConv", 2., Nconv_POGTIP*weight, 0., 10., 10);
       FillHist("IDeff_NEleIDConv", 3., Nconv_POGTIPConv*weight, 0., 10., 10);
       FillHist("IDeff_NEleIDConv", 4., Nconv_POGTIPQ*weight, 0., 10., 10);
       FillHist("IDeff_NEleIDConv", 5., Nconv_MVAM*weight, 0., 10., 10);
       FillHist("IDeff_NEleIDConv", 6., Nconv_MVAMQ*weight, 0., 10., 10);
       FillHist("IDeff_NEleIDConv", 7., Nconv_MVAT*weight, 0., 10., 10);
       FillHist("IDeff_NEleIDConv", 8., Nconv_MVATQ*weight, 0., 10., 10);
       FillHist("IDeff_NEleIDConv", 9., Nconv_HNMVAT*weight, 0., 10., 10);
  
     }//Endof Zgamma
   }//End of Conversion Study
   if(MatchingEff){
     

     /**************************************************
     **ELECTRON MATCHING EFF STUDY
     **************************************************/
     for(int j=0; j<electronPreColl.size(); j++){
   
       int LeptonType=0;
       int MatchedTruthIdx   = GenMatchedIdx(electronPreColl.at(j),truthColl);
       int LastSelfIdx       = LastSelfMotherIdx(MatchedTruthIdx,truthColl);
       int MotherIdx         = FirstNonSelfMotherIdx(MatchedTruthIdx, truthColl);
       int LastSelfMIdx      = LastSelfMotherIdx(MotherIdx,truthColl);
       int GrMotherIdx       = FirstNonSelfMotherIdx(MotherIdx, truthColl);
       int MPID=-1, GrMPID=-1;
         if(MotherIdx!=-1)   MPID=truthColl.at(MotherIdx).PdgId();
         if(GrMotherIdx!=-1) GrMPID=truthColl.at(GrMotherIdx).PdgId();
     
       if     (MatchedTruthIdx==-1) LeptonType=-1;
       else if(fabs(MPID)>50)       LeptonType=-2;
       else if(fabs(MPID)==15){
                if(truthColl.at(MotherIdx).GenStatus()==2){
                  if(fabs(GrMPID)==23 || fabs(GrMPID)==24) LeptonType=3;
                  else if(truthColl.at(LastSelfMIdx).GenStatus()>20 && truthColl.at(LastSelfMIdx).GenStatus()<30) LeptonType=3;
                  else if(fabs(GrMPID)>50)                 LeptonType=-3;
                }
                else                                       LeptonType=3;//Assigned to remove margin but do we need this?
              }
       else if(fabs(MPID)==23 || fabs(MPID)==24 ) LeptonType=1;
       else if(fabs(MPID)==36 || fabs(MPID)==32 ) LeptonType=2;
       else if(truthColl.at(LastSelfIdx).GenStatus()>20 && truthColl.at(LastSelfIdx).GenStatus()<30) LeptonType=5;
       else LeptonType=4;
     
       if(LeptonType==0){
     //    PrintTruth();
     //    cout<<"Idx "<<MatchedTruthIdx<<" MIdx_ "<<MotherIdx<<" _MIdx "<<LastSelfMIdx<<" GrMIdx "<<GrMotherIdx<<" MPID "<<MPID<<" GPID "<<GrMPID<<" LepType "<<LeptonType<<endl;
       }
       FillHist("LepType", LeptonType, weight, -10., 10., 20);

       if(MatchedTruthIdx!=-1) FillHist("dRgenTele", truthColl.at(MatchedTruthIdx).DeltaR(electronPreColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==1) FillHist("dRgenTele_Type1", truthColl.at(MatchedTruthIdx).DeltaR(electronPreColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==3) FillHist("dRgenTele_Type3", truthColl.at(MatchedTruthIdx).DeltaR(electronPreColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==5) FillHist("dRgenTele_Type5", truthColl.at(MatchedTruthIdx).DeltaR(electronPreColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==-2) FillHist("dRgenTele_TypeM2", truthColl.at(MatchedTruthIdx).DeltaR(electronPreColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==-3) FillHist("dRgenTele_TypeM3", truthColl.at(MatchedTruthIdx).DeltaR(electronPreColl.at(j)), weight, 0., 0.02, 2000);
   
       float dR=999.;
       for(int i=2; i<truthColl.size(); i++){
         if(truthColl.at(i).IndexMother()<0 )  continue;
         if(truthColl.at(i).GenStatus()!=1)    continue;
         if(fabs(truthColl.at(i).PdgId())!=11) continue;
      
         if(truthColl.at(i).DeltaR(electronPreColl.at(j))<dR){ dR=truthColl.at(i).DeltaR(electronPreColl.at(j));}
       }
       if(LeptonType==-1){
         FillHist("dRgenTele_TypeM1_bW", dR, weight, 0.1, 5.1, 500);
         if(dR<0.4){
           FillHist("IDeff_NEle_TypeM1dR04", 0., weight, 0., 2., 2);
           if(electronPreColl.at(j).PassTight() && electronPreColl.at(j).PFRelIso(0.3)<0.05 && electronPreColl.at(j).dxy()<0.05 && electronPreColl.at(j).dz()<0.1 && electronPreColl.at(j).dxySig()<3. ) FillHist("IDeff_NEleID_TypeM1dR04", 0., weight, 0., 2., 2);
         }
         else{
           FillHist("IDeff_NEle_TypeM1dR04", 1., weight, 0., 2., 2);
           if(electronPreColl.at(j).PassTight() && electronPreColl.at(j).PFRelIso(0.3)<0.05 && electronPreColl.at(j).dxy()<0.05 && electronPreColl.at(j).dz()<0.1 && electronPreColl.at(j).dxySig()<3. ) FillHist("IDeff_NEleID_TypeM1dR04", 1., weight, 0., 2., 2);
         }
       }
       if(MatchedTruthIdx!=-1) FillHist("dRgenTele_bW", truthColl.at(MatchedTruthIdx).DeltaR(electronPreColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==1) FillHist("dRgenTele_Type1_bW", truthColl.at(MatchedTruthIdx).DeltaR(electronPreColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==3) FillHist("dRgenTele_Type3_bW", truthColl.at(MatchedTruthIdx).DeltaR(electronPreColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==5) FillHist("dRgenTele_Type5_bW", truthColl.at(MatchedTruthIdx).DeltaR(electronPreColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==-2) FillHist("dRgenTele_TypeM2_bW", truthColl.at(MatchedTruthIdx).DeltaR(electronPreColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==-3) FillHist("dRgenTele_TypeM3_bW", truthColl.at(MatchedTruthIdx).DeltaR(electronPreColl.at(j)), weight, 0., 0.1, 10);

     }//End of ElePrecoll loop(EleMatchingStudy)

     for(int j=0; j<electronPOGTIPColl.size(); j++){
   
       int LeptonType=0;
       int MatchedTruthIdx   = GenMatchedIdx(electronPOGTIPColl.at(j),truthColl);
       int LastSelfIdx       = LastSelfMotherIdx(MatchedTruthIdx,truthColl);
       int MotherIdx         = FirstNonSelfMotherIdx(MatchedTruthIdx, truthColl);
       int LastSelfMIdx      = LastSelfMotherIdx(MotherIdx,truthColl);
       int GrMotherIdx       = FirstNonSelfMotherIdx(MotherIdx, truthColl);
       int MPID=-1, GrMPID=-1;
         if(MotherIdx!=-1)   MPID=truthColl.at(MotherIdx).PdgId();
         if(GrMotherIdx!=-1) GrMPID=truthColl.at(GrMotherIdx).PdgId();
     
       if     (MatchedTruthIdx==-1) LeptonType=-1;
       else if(fabs(MPID)>50)       LeptonType=-2;
       else if(fabs(MPID)==15){
                if(truthColl.at(MotherIdx).GenStatus()==2){
                  if(fabs(GrMPID)==23 || fabs(GrMPID)==24) LeptonType=3;
                  else if(truthColl.at(LastSelfMIdx).GenStatus()>20 && truthColl.at(LastSelfMIdx).GenStatus()<30) LeptonType=3;
                  else if(fabs(GrMPID)>50)                 LeptonType=-3;
                }
                else                                       LeptonType=3;//Assigned to remove margin but do we need this?
              }
       else if(fabs(MPID)==23 || fabs(MPID)==24 ) LeptonType=1;
       else if(fabs(MPID)==36 || fabs(MPID)==32 ) LeptonType=2;
       else if(truthColl.at(LastSelfIdx).GenStatus()>20 && truthColl.at(LastSelfIdx).GenStatus()<30) LeptonType=5;
       else LeptonType=4;
     
       if(LeptonType==0){
     //    PrintTruth();
     //    cout<<"Idx "<<MatchedTruthIdx<<" MIdx_ "<<MotherIdx<<" _MIdx "<<LastSelfMIdx<<" GrMIdx "<<GrMotherIdx<<" MPID "<<MPID<<" GPID "<<GrMPID<<" LepType "<<LeptonType<<endl;
       }
       FillHist("LepType_POGT", LeptonType, weight, -10., 10., 20);

       if(MatchedTruthIdx!=-1) FillHist("dRgenTele_POGT", truthColl.at(MatchedTruthIdx).DeltaR(electronPOGTIPColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==1) FillHist("dRgenTele_Type1_POGT", truthColl.at(MatchedTruthIdx).DeltaR(electronPOGTIPColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==3) FillHist("dRgenTele_Type3_POGT", truthColl.at(MatchedTruthIdx).DeltaR(electronPOGTIPColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==5) FillHist("dRgenTele_Type5_POGT", truthColl.at(MatchedTruthIdx).DeltaR(electronPOGTIPColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==-2) FillHist("dRgenTele_TypeM2_POGT", truthColl.at(MatchedTruthIdx).DeltaR(electronPOGTIPColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==-3) FillHist("dRgenTele_TypeM3_POGT", truthColl.at(MatchedTruthIdx).DeltaR(electronPOGTIPColl.at(j)), weight, 0., 0.02, 2000);
   
       float dR=999.;
       for(int i=2; i<truthColl.size(); i++){
         if(truthColl.at(i).IndexMother()<0 )  continue;
         if(truthColl.at(i).GenStatus()!=1)    continue;
         if(fabs(truthColl.at(i).PdgId())!=11) continue;
      
         if(truthColl.at(i).DeltaR(electronPOGTIPColl.at(j))<dR){ dR=truthColl.at(i).DeltaR(electronPOGTIPColl.at(j));}
       }
       if(LeptonType==-1){
         FillHist("dRgenTele_TypeM1_bW_POGT", dR, weight, 0.1, 5.1, 500);
         if(dR<0.4){
           FillHist("IDeff_NEle_TypeM1dR04_POGT", 0., weight, 0., 2., 2);
         }
         else{
           FillHist("IDeff_NEle_TypeM1dR04_POGT", 1., weight, 0., 2., 2);
         }
       }
       if(MatchedTruthIdx!=-1) FillHist("dRgenTele_bW_POGT", truthColl.at(MatchedTruthIdx).DeltaR(electronPOGTIPColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==1) FillHist("dRgenTele_Type1_bW_POGT", truthColl.at(MatchedTruthIdx).DeltaR(electronPOGTIPColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==3) FillHist("dRgenTele_Type3_bW_POGT", truthColl.at(MatchedTruthIdx).DeltaR(electronPOGTIPColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==5) FillHist("dRgenTele_Type5_bW_POGT", truthColl.at(MatchedTruthIdx).DeltaR(electronPOGTIPColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==-2) FillHist("dRgenTele_TypeM2_bW_POGT", truthColl.at(MatchedTruthIdx).DeltaR(electronPOGTIPColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==-3) FillHist("dRgenTele_TypeM3_bW_POGT", truthColl.at(MatchedTruthIdx).DeltaR(electronPOGTIPColl.at(j)), weight, 0., 0.1, 10);

     }//End of ElePrecoll loop(EleMatchingStudy)



     /**************************************************
     **MUON MATCHING EFF STUDY
     **************************************************/
     for(int j=0; j<muonPreColl.size(); j++){

   
       int LeptonType=0;
       int MatchedTruthIdx   = GenMatchedIdx(muonPreColl.at(j),truthColl);
       int LastSelfIdx       = LastSelfMotherIdx(MatchedTruthIdx,truthColl);
       int MotherIdx         = FirstNonSelfMotherIdx(MatchedTruthIdx, truthColl);
       int LastSelfMIdx      = LastSelfMotherIdx(MotherIdx,truthColl);
       int GrMotherIdx       = FirstNonSelfMotherIdx(MotherIdx, truthColl);
       int MPID=-1, GrMPID=-1;
         if(MotherIdx!=-1)   MPID=truthColl.at(MotherIdx).PdgId();
         if(GrMotherIdx!=-1) GrMPID=truthColl.at(GrMotherIdx).PdgId();
     
       if     (MatchedTruthIdx==-1) LeptonType=-1;
       else if(fabs(MPID)>50)       LeptonType=-2;
       else if(fabs(MPID)==15){
                if(truthColl.at(MotherIdx).GenStatus()==2){
                  if(fabs(GrMPID)==23 || fabs(GrMPID)==24) LeptonType=3;
                  else if(truthColl.at(LastSelfMIdx).GenStatus()>20 && truthColl.at(LastSelfMIdx).GenStatus()<30) LeptonType=3;
                  else if(fabs(GrMPID)>50)                 LeptonType=-3;
                }
                else                                       LeptonType=3;//Assigned to remove margin but do we need this?
              }
       else if(fabs(MPID)==23 || fabs(MPID)==24 ) LeptonType=1;
       else if(fabs(MPID)==36 || fabs(MPID)==32 ) LeptonType=2;
       else if(truthColl.at(LastSelfIdx).GenStatus()>20 && truthColl.at(LastSelfIdx).GenStatus()<30) LeptonType=5;
       else LeptonType=4;
     
       if(LeptonType==0){
     //    PrintTruth();
     //    cout<<"Idx "<<MatchedTruthIdx<<" MIdx_ "<<MotherIdx<<" _MIdx "<<LastSelfMIdx<<" GrMIdx "<<GrMotherIdx<<" MPID "<<MPID<<" GPID "<<GrMPID<<" LepType "<<LeptonType<<endl;
       }
       FillHist("LepType", LeptonType, weight, -10., 10., 20);

       if(MatchedTruthIdx!=-1) FillHist("dRgenTmu", truthColl.at(MatchedTruthIdx).DeltaR(muonPreColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==1) FillHist("dRgenTmu_Type1", truthColl.at(MatchedTruthIdx).DeltaR(muonPreColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==3) FillHist("dRgenTmu_Type3", truthColl.at(MatchedTruthIdx).DeltaR(muonPreColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==5) FillHist("dRgenTmu_Type5", truthColl.at(MatchedTruthIdx).DeltaR(muonPreColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==-2) FillHist("dRgenTmu_TypeM2", truthColl.at(MatchedTruthIdx).DeltaR(muonPreColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==-3) FillHist("dRgenTmu_TypeM3", truthColl.at(MatchedTruthIdx).DeltaR(muonPreColl.at(j)), weight, 0., 0.02, 2000);
   
       float dR=999.;
       for(int i=2; i<truthColl.size(); i++){
         if(truthColl.at(i).IndexMother()<0 )  continue;
         if(truthColl.at(i).GenStatus()!=1)    continue;
         if(fabs(truthColl.at(i).PdgId())!=13) continue;
      
         if(truthColl.at(i).DeltaR(muonPreColl.at(j))<dR){ dR=truthColl.at(i).DeltaR(muonPreColl.at(j));}
       }
       if(LeptonType==-1){
         FillHist("dRgenTmu_TypeM1_bW", dR, weight, 0.1, 5.1, 500);
         if(dR<0.4){
           FillHist("IDeff_NMu_TypeM1dR04", 0., weight, 0., 2., 2);
           if(muonPreColl.at(j).IsTight() && muonPreColl.at(j).RelIso04()<0.15) FillHist("IDeff_NMuID_TypeM1dR04", 0., weight, 0., 2., 2);
         }
         else{
           FillHist("IDeff_NMu_TypeM1dR04", 1., weight, 0., 2., 2);
           if(muonPreColl.at(j).IsTight() && muonPreColl.at(j).RelIso04()<0.15) FillHist("IDeff_NMuID_TypeM1dR04", 1., weight, 0., 2., 2);
         }
       }
       if(MatchedTruthIdx!=-1) FillHist("dRgenTmu_bW", truthColl.at(MatchedTruthIdx).DeltaR(muonPreColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==1) FillHist("dRgenTmu_Type1_bW", truthColl.at(MatchedTruthIdx).DeltaR(muonPreColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==3) FillHist("dRgenTmu_Type3_bW", truthColl.at(MatchedTruthIdx).DeltaR(muonPreColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==5) FillHist("dRgenTmu_Type5_bW", truthColl.at(MatchedTruthIdx).DeltaR(muonPreColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==-2) FillHist("dRgenTmu_TypeM2_bW", truthColl.at(MatchedTruthIdx).DeltaR(muonPreColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==-3) FillHist("dRgenTmu_TypeM3_bW", truthColl.at(MatchedTruthIdx).DeltaR(muonPreColl.at(j)), weight, 0., 0.1, 10);



     }//End of MuonPreColl Loop(Muon Matching Study)


     //dRdist due to ID(POGT - PFMUON) Bias Test
     for(int j=0; j<muonPOGTColl.size(); j++){
   
       int LeptonType=0;
       int MatchedTruthIdx   = GenMatchedIdx(muonPOGTColl.at(j),truthColl);
       int LastSelfIdx       = LastSelfMotherIdx(MatchedTruthIdx,truthColl);
       int MotherIdx         = FirstNonSelfMotherIdx(MatchedTruthIdx, truthColl);
       int LastSelfMIdx      = LastSelfMotherIdx(MotherIdx,truthColl);
       int GrMotherIdx       = FirstNonSelfMotherIdx(MotherIdx, truthColl);
       int MPID=-1, GrMPID=-1;
         if(MotherIdx!=-1)   MPID=truthColl.at(MotherIdx).PdgId();
         if(GrMotherIdx!=-1) GrMPID=truthColl.at(GrMotherIdx).PdgId();
     
       if     (MatchedTruthIdx==-1) LeptonType=-1;
       else if(fabs(MPID)>50)       LeptonType=-2;
       else if(fabs(MPID)==15){
                if(truthColl.at(MotherIdx).GenStatus()==2){
                  if(fabs(GrMPID)==23 || fabs(GrMPID)==24) LeptonType=3;
                  else if(truthColl.at(LastSelfMIdx).GenStatus()>20 && truthColl.at(LastSelfMIdx).GenStatus()<30) LeptonType=3;
                  else if(fabs(GrMPID)>50)                 LeptonType=-3;
                }
                else                                       LeptonType=3;//Assigned to remove margin but do we need this?
              }
       else if(fabs(MPID)==23 || fabs(MPID)==24 ) LeptonType=1;
       else if(fabs(MPID)==36 || fabs(MPID)==32 ) LeptonType=2;
       else if(truthColl.at(LastSelfIdx).GenStatus()>20 && truthColl.at(LastSelfIdx).GenStatus()<30) LeptonType=5;
       else LeptonType=4;
     
       if(LeptonType==0){
     //    PrintTruth();
     //    cout<<"Idx "<<MatchedTruthIdx<<" MIdx_ "<<MotherIdx<<" _MIdx "<<LastSelfMIdx<<" GrMIdx "<<GrMotherIdx<<" MPID "<<MPID<<" GPID "<<GrMPID<<" LepType "<<LeptonType<<endl;
       }
       FillHist("LepType_POGT", LeptonType, weight, -10., 10., 20);

       if(MatchedTruthIdx!=-1) FillHist("dRgenTmu_POGT", truthColl.at(MatchedTruthIdx).DeltaR(muonPOGTColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==1) FillHist("dRgenTmu_Type1_POGT", truthColl.at(MatchedTruthIdx).DeltaR(muonPOGTColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==3) FillHist("dRgenTmu_Type3_POGT", truthColl.at(MatchedTruthIdx).DeltaR(muonPOGTColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==5) FillHist("dRgenTmu_Type5_POGT", truthColl.at(MatchedTruthIdx).DeltaR(muonPOGTColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==-2) FillHist("dRgenTmu_TypeM2_POGT", truthColl.at(MatchedTruthIdx).DeltaR(muonPOGTColl.at(j)), weight, 0., 0.02, 2000);
       if(LeptonType==-3) FillHist("dRgenTmu_TypeM3_POGT", truthColl.at(MatchedTruthIdx).DeltaR(muonPOGTColl.at(j)), weight, 0., 0.02, 2000);
   
       float dR=999.;
       for(int i=2; i<truthColl.size(); i++){
         if(truthColl.at(i).IndexMother()<0 )  continue;
         if(truthColl.at(i).GenStatus()!=1)    continue;
         if(fabs(truthColl.at(i).PdgId())!=13) continue;
      
         if(truthColl.at(i).DeltaR(muonPOGTColl.at(j))<dR){ dR=truthColl.at(i).DeltaR(muonPOGTColl.at(j));}
       }
       if(MatchedTruthIdx!=-1) FillHist("dRgenTmu_bW_POGT", truthColl.at(MatchedTruthIdx).DeltaR(muonPOGTColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==-1) FillHist("dRgenTmu_TypeM1_bW_POGT", dR, weight, 0.1, 5.1, 50);
       if(LeptonType==1) FillHist("dRgenTmu_Type1_bW_POGT", truthColl.at(MatchedTruthIdx).DeltaR(muonPOGTColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==3) FillHist("dRgenTmu_Type3_bW_POGT", truthColl.at(MatchedTruthIdx).DeltaR(muonPOGTColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==5) FillHist("dRgenTmu_Type5_bW_POGT", truthColl.at(MatchedTruthIdx).DeltaR(muonPOGTColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==-2) FillHist("dRgenTmu_TypeM2_bW_POGT", truthColl.at(MatchedTruthIdx).DeltaR(muonPOGTColl.at(j)), weight, 0., 0.1, 10);
       if(LeptonType==-3) FillHist("dRgenTmu_TypeM3_bW_POGT", truthColl.at(MatchedTruthIdx).DeltaR(muonPOGTColl.at(j)), weight, 0., 0.1, 10);

     }//End of MuonTightColl Loop(MuonMatching EffStudy)

   }//End of MatchingEff Study 
 
/////////////////////////////////////////////////////////////////////////////////// 

return;
}// End of execute event loop
  


void May2017_ObjectEff::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void May2017_ObjectEff::BeginCycle() throw( LQError ){
  
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

May2017_ObjectEff::~May2017_ObjectEff() {
  
  Message("In May2017_ObjectEff Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

int May2017_ObjectEff::NPromptFake_Ele(std::vector<snu::KElectron> EleColl, TString Option){

   int Nprompt=0, Nfake=0;

   bool ReturnPrompt=false, ReturnFake=false;
   if     (Option.Contains("Prompt")) ReturnPrompt=true;
   else if(Option.Contains("Fake"))   ReturnFake=true;

   for(int i=0; i<EleColl.size(); i++){ if(EleColl.at(i).MCIsPrompt()){Nprompt++;} else{Nfake++;} }

   if     (ReturnPrompt) return Nprompt;
   else if(ReturnFake)   return Nfake;

   return 0.;
  
}


int May2017_ObjectEff::NPromptFake_Mu(std::vector<snu::KMuon> MuColl, TString Option){

   int Nprompt=0, Nfake=0;

   bool ReturnPrompt=false, ReturnFake=false;
   if     (Option.Contains("Prompt")) ReturnPrompt=true;
   else if(Option.Contains("Fake"))   ReturnFake=true;

   for(int i=0; i<MuColl.size(); i++){ if(MuColl.at(i).MCIsPrompt()){Nprompt++;} else{Nfake++;} }

   if     (ReturnPrompt) return Nprompt;
   else if(ReturnFake)   return Nfake;

   return 0.;
  
}

int May2017_ObjectEff::NPromptFake_Mu(std::vector<snu::KMuon> MuColl, std::vector<snu::KTruth> truthColl, TString Option){

   int Nprompt=0, Nfake=0;

   bool ReturnPrompt=false, ReturnFake=false, ReturnBSM=false;
   if     (Option.Contains("Prompt")) ReturnPrompt=true;
   else if(Option.Contains("Fake"))   ReturnFake=true;
   else if(Option.Contains("BSM"))    ReturnBSM=true;

   for(int i=0; i<MuColl.size(); i++){
     int MuType=GetLeptonType(MuColl.at(i), truthColl);

     if(!ReturnBSM){
       if     (MuType==1) Nprompt++;
       else if(MuType<0) Nfake++;
       else if(MuType==2) Nprompt++;
     }
     else if(ReturnBSM && MuType==2) Nprompt++;
     
   }

   if     (ReturnPrompt || ReturnBSM) return Nprompt;
   else if(ReturnFake)   return Nfake;

   return 0.;
  
}

int May2017_ObjectEff::NPromptFake_Ele(std::vector<snu::KElectron> EleColl, std::vector<snu::KTruth> truthColl, TString Option){

   int Nprompt=0, Nfake=0;

   bool ReturnPrompt=false, ReturnFake=false;
   if     (Option.Contains("Prompt")) ReturnPrompt=true;
   else if(Option.Contains("Fake"))   ReturnFake=true;

   for(int i=0; i<EleColl.size(); i++){
     int EleType=GetLeptonType(EleColl.at(i), truthColl);

     if     (EleType==1) Nprompt++;
     else if(EleType<0) Nfake++;
     else if(EleType==2) Nprompt++;
   }

   if     (ReturnPrompt) return Nprompt;
   else if(ReturnFake)   return Nfake;

   return 0.;
  
}



bool May2017_ObjectEff::IsConvCand(snu::KElectron Ele, std::vector<snu::KTruth> TruthColl, TString Option){

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



void May2017_ObjectEff::FillCutFlow(TString cut, float weight){
  
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



void May2017_ObjectEff::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void May2017_ObjectEff::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //Basic Hists/////////////////////////////////////////////////
  //Data properties before selection  
  AnalyzerCore::MakeHistograms("Basic_b2_Et_orig", 200, 0., 200.);



  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this May2017_ObjectEffCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void May2017_ObjectEff::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
