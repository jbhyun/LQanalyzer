// $Id: Mar2017_3l4j_TriLepComp.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQMar2017_3l4j_TriLepComp Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Mar2017_3l4j_TriLepComp.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Mar2017_3l4j_TriLepComp);

 Mar2017_3l4j_TriLepComp::Mar2017_3l4j_TriLepComp() : AnalyzerCore(), out_muons(0) {

   SetLogName("Mar2017_3l4j_TriLepComp");
   Message("In Mar2017_3l4j_TriLepComp constructor", INFO);
   InitialiseAnalysis();
 }


 void Mar2017_3l4j_TriLepComp::InitialiseAnalysis() throw( LQError ) {
   
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

void Mar2017_3l4j_TriLepComp::ExecuteEvents()throw( LQError ){

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


   //bool TriMu_analysis=true, EMuMu_analysis=false, DiMuon_analysis=false, DiEle_analysis=false;
   bool TriMu_analysis=false, EMuMu_analysis=true, DiMuon_analysis=false, DiEle_analysis=false;
   //bool TriMu_analysis=false, EMuMu_analysis=false, DiMuon_analysis=true, DiEle_analysis=false;
   //bool TriMu_analysis=false, EMuMu_analysis=false, DiMuon_analysis=false, DiEle_analysis=true;

   //bool FakeEstimation=false;

    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //ListTriggersAvailable();

   //Normalisation Lumi Setting
   float trigger_ps_weight=1.;
   if(!isData) trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
   weight*=trigger_ps_weight;
   FillHist("TriggerPSWeight", trigger_ps_weight, 1., 0., 1., 100);

   //Trigger Path of Analysis
   bool Pass_Trigger=false;
   //if(false){
   if(EMuMu_analysis){
     //if( PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v") ) Pass_Trigger=true;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
          || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
          || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger=true;
   }
   else if(TriMu_analysis){
   //else if(EMuMu_analysis || TriMu_analysis){
     //if( PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v") ) Pass_Trigger=true;
     //if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ) Pass_Trigger=true;
     if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")
          || PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") ) Pass_Trigger=true;
   }
   else if(DiMuon_analysis){
     //if( PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v") ) Pass_Trigger=true;
     if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ) Pass_Trigger=true;
   }
   else if(DiEle_analysis){
     if( PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger=true;
   }

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
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, true);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
     if     (TriMu_analysis) { if( !(muonPreColl.size()>=3)     ) return; }
     else if(EMuMu_analysis) { if( !(muonPreColl.size()>=2 && electronPreColl.size()>=1) ) return; }
   /**********************************************************************************************************/

     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.25);
   std::vector<snu::KMuon> muonPOGLColl; eventbase->GetMuonSel()->Selection(muonPOGLColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.6);
   std::vector<snu::KMuon> muonVetoColl; eventbase->GetMuonSel()->Selection(muonVetoColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.05);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.05);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.1);
   std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,true);
   std::vector<snu::KMuon> muonColl;      if(k_running_nonprompt){ muonColl=muonLooseColl;} else{ muonColl=muonTightColl;}
   //std::vector<snu::KMuon> muonPromptColl; muonPromptColl=GetTruePrompt(muonColl, false);

     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0994, 0.107);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.2);//Not in ID, but additional safe WP
   std::vector<snu::KElectron> electronVetoColl; eventbase->GetElectronSel()->Selection(electronVetoColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(20.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.5, 0.5);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.2);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetCheckCharge(true);
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(20.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.2);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetCheckCharge(true);
   std::vector<snu::KElectron> electronTightColl; eventbase->GetElectronSel()->Selection(electronTightColl);
   std::vector<snu::KElectron> electronColl;      if(k_running_nonprompt){ electronColl=electronLooseColl; } else{ electronColl=electronTightColl; }
   std::vector<snu::KElectron> electronNull;

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     //eventbase->GetJetSel()->SetPileUpJetID(true,"Loose");
     bool LeptonVeto=true;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonLooseColl, electronLooseColl);
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(20.);                     eventbase->GetJetSel()->SetEta(2.4);
     LeptonVeto=true;
   std::vector<snu::KJet> jet20POGLVetoColl; eventbase->GetJetSel()->Selection(jet20POGLVetoColl, LeptonVeto, muonPOGLColl, electronVetoColl);
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(30.);                     eventbase->GetJetSel()->SetEta(2.4);
     LeptonVeto=true;
   std::vector<snu::KJet> jet30POGLVetoColl; eventbase->GetJetSel()->Selection(jet30POGLVetoColl, LeptonVeto, muonPOGLColl, electronVetoColl);
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(20.);                     eventbase->GetJetSel()->SetEta(2.4);
     LeptonVeto=true;
   std::vector<snu::KJet> jet20HNVVetoColl; eventbase->GetJetSel()->Selection(jet20HNVVetoColl, LeptonVeto, muonVetoColl, electronVetoColl);
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(30.);                     eventbase->GetJetSel()->SetEta(2.4);
     LeptonVeto=true;
   std::vector<snu::KJet> jet30HNVVetoColl; eventbase->GetJetSel()->Selection(jet30HNVVetoColl, LeptonVeto, muonVetoColl, electronVetoColl);
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(20.);                     eventbase->GetJetSel()->SetEta(2.4);
     LeptonVeto=true;
   std::vector<snu::KJet> jet20FakeLVetoColl; eventbase->GetJetSel()->Selection(jet20FakeLVetoColl, LeptonVeto, muonLooseColl, electronVetoColl);
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(30.);                     eventbase->GetJetSel()->SetEta(2.4);
     LeptonVeto=true;
   std::vector<snu::KJet> jet30FakeLVetoColl; eventbase->GetJetSel()->Selection(jet30FakeLVetoColl, LeptonVeto, muonLooseColl, electronVetoColl);
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(20.);                     eventbase->GetJetSel()->SetEta(2.4);
     LeptonVeto=true;
   std::vector<snu::KJet> jet20FakeTVetoColl; eventbase->GetJetSel()->Selection(jet20FakeTVetoColl, LeptonVeto, muonTightColl, electronVetoColl);
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(30.);                     eventbase->GetJetSel()->SetEta(2.4);
     LeptonVeto=true;
   std::vector<snu::KJet> jet30FakeTVetoColl; eventbase->GetJetSel()->Selection(jet30FakeTVetoColl, LeptonVeto, muonTightColl, electronVetoColl);



   //Method to apply 2a method SF
   //std::vector<int> bIdxColl  = GetSFBJetIdx(jetColl,"Medium");
   //std::vector<int> ljIdxColl = GetSFLJetIdx(jetColl, bIdxColl, "Medium");
   //std::vector<snu::KJet> bjetColl2a; for(int i=0; i<bIdxColl.size(); i++) {bjetColl2a.push_back(jetColl.at(bIdxColl.at(i)));}
   //std::vector<snu::KJet> ljetColl2a; for(int i=0; i<ljIdxColl.size(); i++){ljetColl2a.push_back(jetColl.at(ljIdxColl.at(i)));}

   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");
   std::vector<snu::KJet> ljetColl = SelLightJets(jetColl, "Medium");

   double met = eventbase->GetEvent().PFMETType1();
   double met_x = eventbase->GetEvent().PFMETType1x();
   double met_y = eventbase->GetEvent().PFMETType1y();
   //double met = eventbase->GetEvent().MET();
   double Pzv,Pzv1, Pzv2;
   snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));
   //snu::KParticle v[4]; v[0].SetPx(met_x); v[0].SetPy(met_y);
   //                     v[1].SetPx(met_x); v[1].SetPy(met_y);

   int nbjets=bjetColl.size(); int njets=jetColl.size(); const int nljets=ljetColl.size();//njets-nbjets;//number of light jets
   int Nvtx=eventbase->GetEvent().nVertices();


   /*****************************************************
   **Scale Factors
   *****************************************************/
   float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1., btag_sf=1.;
   float trigger_sf=1.;
   float fake_weight=1.; bool EventCand=false;

   /*This part is for boosting up speed.. SF part takes rather longer time than expected*/
   if     (TriMu_analysis) { if(muonLooseColl.size()==3)     EventCand=true; }
   else if(EMuMu_analysis) { if(muonLooseColl.size()==2 && electronLooseColl.size()==1) EventCand=true; }
   else if(DiMuon_analysis){ if(muonLooseColl.size()==2)     EventCand=true; }
   else if(DiEle_analysis) { if(electronLooseColl.size()==2) EventCand=true; }


   if(EventCand){
     if(!isData){
      //trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );
  
       id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronColl);
       reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
  
       id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muonColl);
       //iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
       trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
  
       btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);
     }
     else{
       //Perfect Prompt Ratio Approximation applied.(p=1), Applicable to generic number, combination of leptons under premise of high prompt rate.
       if(k_running_nonprompt){
         fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", muonLooseColl.size(), electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", electronLooseColl.size());
       }
     }
   }

   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight*fake_weight*btag_sf;
   /***************************************************************************************************/

   //////Basic Objects Check//////////////////////////////////////////////////////////
   FillHist("Basic_Nj_orig", njets, weight, 0., 10., 10);
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
   FillHist("Basic_METdist_orig", met, weight, 0., 200., 100);
   ///////////////////////////////////////////////////////////////////////////////////


/************************************************************************************/
//////////////////////////////////////////////////////////////////////////////////////
/////// MAIN ANALYSIS CODE ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
/************************************************************************************/


   if(EMuMu_analysis){

     //Step1 : 1e+2mu +PTCut
     //if( !(muonColl.at(0).Pt()>27) ) return; //SingleMuon Trig Case
     if( !(electronLooseColl.size()==1 && muonLooseColl.size()==2) ) return;
     if( !(electronColl.size()==1 && muonColl.size()==2) ) return;
     //if( !(electronColl.at(0).Pt()>20 && muonColl.at(0).Pt()>20 && muonColl.at(1).Pt()>10) ) return;
     if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>15 && muonColl.at(1).Pt()>10) ) return;
     if( SumCharge(muonColl)!=0 ) return;
     FillCutFlow("3lCut", weight);
     FillHist("WeightDist", weight, 1, -10., 10., 2000);


     int IdxZCandLead = GetDaughterCandIdx(muonColl, "Z", 10., "Lead");
     int IdxZCandSubl = GetDaughterCandIdx(muonColl, "Z", 10., "Subl");

     float Mmumu=(muonColl.at(0)+muonColl.at(1)).M();
     float MTW = sqrt( 2*(met*electronColl.at(0).Pt()-met_x*electronColl.at(0).Px()-met_y*electronColl.at(0).Py()) );


     //if(bjetColl.size()!=0) return;;

     if(IdxZCandLead!=-1 && IdxZCandSubl!=-1){
       FillHist("Mmumu_Zwin10_1e2mu", (muonColl.at(IdxZCandLead)+muonColl.at(IdxZCandSubl)).M(), weight, 0., 200., 200);
       FillHist("PTe_Zwin10_1e2mu", electronColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTmu1_Zwin10_1e2mu", muonColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTmu2_Zwin10_1e2mu", muonColl.at(1).Pt(), weight, 0., 200., 200);
       FillHist("Etae_Zwin10_1e2mu", electronColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Etamu1_Zwin10_1e2mu", muonColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Etamu2_Zwin10_1e2mu", muonColl.at(1).Eta(), weight, -5., 5., 100);
     }


     //Generic 3lep agreement
     FillHist("PTe_3lOSCut_1e2mu", electronColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("PTmu1_3lOSCut_1e2mu", muonColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("PTmu2_3lOSCut_1e2mu", muonColl.at(1).Pt(), weight, 0., 200., 200);
     FillHist("Etae_3lOSCut_1e2mu", electronColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etamu1_3lOSCut_1e2mu", muonColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etamu2_3lOSCut_1e2mu", muonColl.at(1).Eta(), weight, -5., 5., 100);

     FillHist("Nj_3lOSCut_1e2mu", njets, weight, 0., 10., 10);
     FillHist("Nb_3lOSCut_1e2mu", nbjets, weight, 0., 10., 10);
     FillHist("Nlj_3lOSCut_1e2mu", nljets, weight, 0., 10., 10);
     FillHist("NlVeto_3lOSCut_1e2mu", muonVetoColl.size()+electronVetoColl.size(), weight, 0., 10., 10);
     FillHist("NeVeto_3lOSCut_1e2mu", electronVetoColl.size(), weight, 0., 10., 10);
     FillHist("NmuVeto_3lOSCut_1e2mu", muonVetoColl.size(), weight, 0., 10., 10);
     FillHist("NmuL_3lOSCut_1e2mu", muonLooseColl.size(), weight, 0., 10., 10);
     FillHist("NeL_3lOSCut_1e2mu", electronLooseColl.size(), weight, 0., 10., 10);

     FillHist("MET_3lOSCut_1e2mu", met, weight, 0., 200., 200); 
     FillHist("Mmumu_3lOSCut_1e2mu", Mmumu, weight, 0., 200., 200);
     FillHist("M3l_3lOSCut_1e2mu", (muonColl.at(0)+muonColl.at(1)+electronColl.at(0)).M(), weight, 0., 500., 500);


     
     if(Mmumu>4){
       FillHist("PTe_3lOSM12Cut_1e2mu", electronColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTmu1_3lOSM12Cut_1e2mu", muonColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTmu2_3lOSM12Cut_1e2mu", muonColl.at(1).Pt(), weight, 0., 200., 200);
       FillHist("Etae_3lOSM12Cut_1e2mu", electronColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Etamu1_3lOSM12Cut_1e2mu", muonColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Etamu2_3lOSM12Cut_1e2mu", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
       FillHist("Nj_3lOSM12Cut_1e2mu", njets, weight, 0., 10., 10);
       FillHist("Nb_3lOSM12Cut_1e2mu", nbjets, weight, 0., 10., 10);
       FillHist("Nlj_3lOSM12Cut_1e2mu", nljets, weight, 0., 10., 10);
       FillHist("NlVeto_3lOSM12Cut_1e2mu", muonVetoColl.size()+electronVetoColl.size(), weight, 0., 10., 10);
       FillHist("NeVeto_3lOSM12Cut_1e2mu", electronVetoColl.size(), weight, 0., 10., 10);
       FillHist("NmuVeto_3lOSM12Cut_1e2mu", muonVetoColl.size(), weight, 0., 10., 10);
       FillHist("NmuL_3lOSM12Cut_1e2mu", muonLooseColl.size(), weight, 0., 10., 10);
       FillHist("NeL_3lOSCut_1e2mu", electronLooseColl.size(), weight, 0., 10., 10);

       FillHist("MET_3lOSM12Cut_1e2mu", met, weight, 0., 200., 200); 
       FillHist("dRmumu_3lOSM12Cut_1e2mu", muonColl.at(0).DeltaR(muonColl.at(1)), weight, 0., 5., 500);
       FillHist("dRemu1_3lOSM12Cut_1e2mu", electronColl.at(0).DeltaR(muonColl.at(0)), weight, 0., 5., 500);
       FillHist("dRemu2_3lOSM12Cut_1e2mu", electronColl.at(0).DeltaR(muonColl.at(1)), weight, 0., 5., 500);
       FillHist("Mmumu_3lOSCut_1e2mu", Mmumu, weight, 0., 200., 200);
       FillHist("M3l_3lOSCut_1e2mu", (muonColl.at(0)+muonColl.at(1)+electronColl.at(0)).M(), weight, 0., 500., 500);
       FillHist("MTW_3lOSM12Cut_1e2mu", MTW, weight, 0., 300., 300);
     }


   }
   if(TriMu_analysis){

     //Step1 : 3mu +PTCut+ SSS veto
     //if( !(muonVetoColl.size()==3 && electronVetoColl.size()==0) ) return;//No additional leptons - keep 4th mu ? ~ 5% These guys have impact of ~25%
     if( !(muonLooseColl.size()==3) ) return;//For valid usage of fakeable object method. FakeLooseSel should be exactly 3.
     int NLmF=0;
     for(int i=0; i<muonVetoColl.size(); i++){
       if(k_running_nonprompt){
         if(!(muonVetoColl.at(i).IsTight() && muonVetoColl.at(i).dXY()<0.05 && muonVetoColl.at(i).dXYSig()<3. && muonVetoColl.at(i).RelIso04()*muonVetoColl.at(i).Pt()/muonVetoColl.at(i).RochPt()<0.4)) NLmF++; }
       else{
         if(!(muonVetoColl.at(i).IsTight() && muonVetoColl.at(i).dXY()<0.05 && muonVetoColl.at(i).dXYSig()<3. && muonVetoColl.at(i).RelIso04()*muonVetoColl.at(i).Pt()/muonVetoColl.at(i).RochPt()<0.1)) NLmF++; }
     }
     FillHist("NMuPOGL_NoFakeL_NFakeL3Cut", NLmF, weight, 0., 10., 10);
     //if(NLmF!=0) return;//No additional isolated loose muon

     if( !(muonColl.size()==3     && electronColl.size()==0) ) return;
     if( !(muonColl.at(0).Pt()>20. && muonColl.at(2).Pt()>10.) ) return;
     if( fabs(SumCharge(muonColl))!=1 ) return;
     FillHist("NMuPOGL_NoFakeL_NFakeT3Cut", NLmF, weight, 0., 10., 10);
     //if( k_sample_name.Contains("ZGto2LG") ){}
     FillCutFlow("3lCut", weight);

     FillHist("WeightDist", weight, 1, -10., 10., 2000);

     //General variables
     int IdxZCandLead = GetDaughterCandIdx(muonColl, "Z", 10., "Lead");
     int IdxZCandSubl = GetDaughterCandIdx(muonColl, "Z", 10., "Subl");

     int IdxOS  = TriMuChargeIndex(muonColl,"OS");
     int IdxSS1 = TriMuChargeIndex(muonColl,"SS1");
     int IdxSS2 = TriMuChargeIndex(muonColl,"SS2");

     float Mmumu_OSSS1=(muonColl.at(IdxOS)+muonColl.at(IdxSS1)).M();
     float Mmumu_OSSS2=(muonColl.at(IdxOS)+muonColl.at(IdxSS2)).M();

     int IdxZmuOS  = IdxOS;
     int IdxZmuSS  = fabs(Mmumu_OSSS1-91.2)<fabs(Mmumu_OSSS2-91.2) ? IdxSS1 : IdxSS2;
     int IdxNonZSS = fabs(Mmumu_OSSS1-91.2)<fabs(Mmumu_OSSS2-91.2) ? IdxSS2 : IdxSS1;
     
     float M3l = (muonColl.at(0)+muonColl.at(1)+muonColl.at(2)).M();

     int IdxZCandLead_NoLim = GetDaughterCandIdx(muonColl, "Z", -1., "Lead");
     int IdxZCandSubl_NoLim = GetDaughterCandIdx(muonColl, "Z", -1., "Subl");
     int IdxWCand_NoLim = 3-IdxZCandLead_NoLim-IdxZCandSubl_NoLim; 
     FillHist("Mmumu_ZCand", (muonColl.at(IdxZCandLead_NoLim)+muonColl.at(IdxZCandSubl_NoLim)).M(), weight, 0., 200., 200);
     //float MTW = (muonColl.at(IdxWCand_NoLim)+v).Mt();
     float MTW = sqrt(met*muonColl.at(IdxWCand_NoLim).Pt()-met_x*muonColl.at(IdxWCand_NoLim).Px()-met_y*muonColl.at(IdxWCand_NoLim).Py());


     //CR - TriLep + b veto
     if(bjetColl.size()!=0) return;;

//     FillCutFlow("bVeto", weight);
     //3lep Z peak test : Fake Control Check
     if(IdxZCandLead!=-1 && IdxZCandSubl!=-1){
       FillHist("Mmumu_Zwin10_3mu", (muonColl.at(IdxZCandLead)+muonColl.at(IdxZCandSubl)).M(), weight, 0., 200., 200);
       FillHist("PTmu1_Zwin10_3mu", muonColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTmu2_Zwin10_3mu", muonColl.at(1).Pt(), weight, 0., 200., 200);
       FillHist("PTmu3_Zwin10_3mu", muonColl.at(2).Pt(), weight, 0., 200., 200);
       FillHist("Etamu1_Zwin10_3mu", muonColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Etamu2_Zwin10_3mu", muonColl.at(1).Eta(), weight, -5., 5., 100);
       FillHist("Etamu3_Zwin10_3mu", muonColl.at(2).Eta(), weight, -5., 5., 100);
     }


     //Generic 3lep agreement
     FillHist("PTmu1_3lOSCut_3mu", muonColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("PTmu2_3lOSCut_3mu", muonColl.at(1).Pt(), weight, 0., 200., 200);
     FillHist("PTmu3_3lOSCut_3mu", muonColl.at(2).Pt(), weight, 0., 200., 200);
     FillHist("Etamu1_3lOSCut_3mu", muonColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etamu2_3lOSCut_3mu", muonColl.at(1).Eta(), weight, -5., 5., 100);
     FillHist("Etamu3_3lOSCut_3mu", muonColl.at(2).Eta(), weight, -5., 5., 100);

     FillHist("Nj_3lOSCut_3mu", jetColl.size(), weight, 0., 10., 10);
     FillHist("NlVeto_3lOSCut_3mu", muonVetoColl.size()+electronVetoColl.size(), weight, 0., 10., 10);
     FillHist("NeVeto_3lOSCut_3mu", electronVetoColl.size(), weight, 0., 10., 10);
     FillHist("NmuVeto_3lOSCut_3mu", muonVetoColl.size(), weight, 0., 10., 10);
     FillHist("NmuL_3lOSCut_3mu", muonLooseColl.size(), weight, 0., 10., 10);

     FillHist("MET_3lOSCut_3mu", met, weight, 0., 200., 200); 
     FillHist("Nb_3lOSCut_3mu", bjetColl.size(), weight, 0., 10., 10);
     FillHist("Mmumu_OSSS1_3lOSCut_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS1)).M(), weight, 0., 200., 200);
     FillHist("Mmumu_OSSS2_3lOSCut_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS2)).M(), weight, 0., 200., 200);
     FillHist("M3mu_3lOSCut_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS1)+muonColl.at(IdxSS2)).M(), weight, 0., 500., 500);

     //Large Difference at Low Mass, but we don't have low prompt estimate from MC. 
     //So Remove low mass part
     if(Mmumu_OSSS2>4 && Mmumu_OSSS1>4){
       FillHist("PTmu1_3lOSM12Cut_3mu", muonColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTmu2_3lOSM12Cut_3mu", muonColl.at(1).Pt(), weight, 0., 200., 200);
       FillHist("PTmu3_3lOSM12Cut_3mu", muonColl.at(2).Pt(), weight, 0., 200., 200);
       FillHist("Etamu1_3lOSM12Cut_3mu", muonColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Etamu2_3lOSM12Cut_3mu", muonColl.at(1).Eta(), weight, -5., 5., 100);
       FillHist("Etamu3_3lOSM12Cut_3mu", muonColl.at(2).Eta(), weight, -5., 5., 100);
  
       FillHist("Nj_3lOSM12Cut_3mu", jetColl.size(), weight, 0., 10., 10);
       FillHist("NlVeto_3lOSM12Cut_3mu", muonVetoColl.size()+electronVetoColl.size(), weight, 0., 10., 10);
       FillHist("NeVeto_3lOSM12Cut_3mu", electronVetoColl.size(), weight, 0., 10., 10);
       FillHist("NmuVeto_3lOSM12Cut_3mu", muonVetoColl.size(), weight, 0., 10., 10);
       FillHist("NmuL_3lOSM12Cut_3mu", muonLooseColl.size(), weight, 0., 10., 10);

       FillHist("MET_3lOSM12Cut_3mu", met, weight, 0., 200., 200); 
       FillHist("Nb_3lOSM12Cut_3mu", bjetColl.size(), weight, 0., 10., 10);
       FillHist("dRmumu_OSSS1_3lOSM12Cut_3mu", muonColl.at(IdxOS).DeltaR(muonColl.at(IdxSS1)), weight, 0., 5., 500);
       FillHist("dRmumu_OSSS2_3lOSM12Cut_3mu", muonColl.at(IdxOS).DeltaR(muonColl.at(IdxSS2)), weight, 0., 5., 500);
       FillHist("Mmumu_OSSS1_3lOSM12Cut_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS1)).M(), weight, 0., 200., 200);
       FillHist("Mmumu_OSSS2_3lOSM12Cut_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS2)).M(), weight, 0., 200., 200);
       FillHist("M3mu_3lOSM12Cut_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS1)+muonColl.at(IdxSS2)).M(), weight, 0., 500., 500);
       FillHist("MTW_3lOSM12Cut_3mu", MTW, weight, 0., 300., 300);
      

       //Njet related Study
       FillHist("Nj20_POGLVeto_3lOSM12Cut_3mu", jet20POGLVetoColl.size(), weight, 0., 10., 10);
       FillHist("Nj30_POGLVeto_3lOSM12Cut_3mu", jet30POGLVetoColl.size(), weight, 0., 10., 10);
       FillHist("Nj20_HNVVeto_3lOSM12Cut_3mu",  jet20HNVVetoColl.size(), weight, 0., 10., 10);
       FillHist("Nj30_HNVVeto_3lOSM12Cut_3mu",  jet30HNVVetoColl.size(), weight, 0., 10., 10);
       FillHist("Nj20_FakeLVeto_3lOSM12Cut_3mu", jet20FakeLVetoColl.size(), weight, 0., 10., 10);
       FillHist("Nj30_FakeLVeto_3lOSM12Cut_3mu", jet30FakeLVetoColl.size(), weight, 0., 10., 10);
       FillHist("Nj20_FakeTVeto_3lOSM12Cut_3mu", jet20FakeTVetoColl.size(), weight, 0., 10., 10);
       FillHist("Nj30_FakeTVeto_3lOSM12Cut_3mu", jet30FakeTVetoColl.size(), weight, 0., 10., 10);

       if(njets==0){
         FillHist("PTmu1_3lOSM12Cut_0j_3mu", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_3lOSM12Cut_0j_3mu", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("PTmu3_3lOSM12Cut_0j_3mu", muonColl.at(2).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_3lOSM12Cut_0j_3mu", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_3lOSM12Cut_0j_3mu", muonColl.at(1).Eta(), weight, -5., 5., 100);
         FillHist("Etamu3_3lOSM12Cut_0j_3mu", muonColl.at(2).Eta(), weight, -5., 5., 100);
         FillHist("MET_3lOSM12Cut_0j_3mu", met, weight, 0., 200., 200); 
         FillHist("Mmumu_OSSS1_3lOSM12Cut_0j_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS1)).M(), weight, 0., 200., 200);
         FillHist("Mmumu_OSSS2_3lOSM12Cut_0j_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS2)).M(), weight, 0., 200., 200);
         FillHist("M3mu_3lOSM12Cut_0j_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS1)+muonColl.at(IdxSS2)).M(), weight, 0., 500., 500);
         FillHist("MTW_3lOSM12Cut_0j_3mu", MTW, weight, 0., 300., 300);
       }
       else if(njets==1){
         FillHist("PTmu1_3lOSM12Cut_1j_3mu", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_3lOSM12Cut_1j_3mu", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("PTmu3_3lOSM12Cut_1j_3mu", muonColl.at(2).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_3lOSM12Cut_1j_3mu", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_3lOSM12Cut_1j_3mu", muonColl.at(1).Eta(), weight, -5., 5., 100);
         FillHist("Etamu3_3lOSM12Cut_1j_3mu", muonColl.at(2).Eta(), weight, -5., 5., 100);
         FillHist("MET_3lOSM12Cut_1j_3mu", met, weight, 0., 200., 200); 
         FillHist("Mmumu_OSSS1_3lOSM12Cut_1j_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS1)).M(), weight, 0., 200., 200);
         FillHist("Mmumu_OSSS2_3lOSM12Cut_1j_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS2)).M(), weight, 0., 200., 200);
         FillHist("M3mu_3lOSM12Cut_1j_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS1)+muonColl.at(IdxSS2)).M(), weight, 0., 500., 500);
         FillHist("MTW_3lOSM12Cut_1j_3mu", MTW, weight, 0., 300., 300);
       }
       else if(njets==2){
         FillHist("PTmu1_3lOSM12Cut_2j_3mu", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_3lOSM12Cut_2j_3mu", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("PTmu3_3lOSM12Cut_2j_3mu", muonColl.at(2).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_3lOSM12Cut_2j_3mu", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_3lOSM12Cut_2j_3mu", muonColl.at(1).Eta(), weight, -5., 5., 100);
         FillHist("Etamu3_3lOSM12Cut_2j_3mu", muonColl.at(2).Eta(), weight, -5., 5., 100);
         FillHist("MET_3lOSM12Cut_2j_3mu", met, weight, 0., 200., 200); 
         FillHist("Mmumu_OSSS1_3lOSM12Cut_2j_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS1)).M(), weight, 0., 200., 200);
         FillHist("Mmumu_OSSS2_3lOSM12Cut_2j_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS2)).M(), weight, 0., 200., 200);
         FillHist("M3mu_3lOSM12Cut_2j_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS1)+muonColl.at(IdxSS2)).M(), weight, 0., 500., 500);
         FillHist("MTW_3lOSM12Cut_2j_3mu", MTW, weight, 0., 300., 300);
       }

      
       

       /*-----------------------------------------------------------------------------------------------------------------------*/
       //Most difference removed by 4GeV cut. But still have M(3l)~90 difference
       //Need to test whether this is from ZZTo4L & out of 4GeV gen cut and one of them out of acceptance
       FillHist("PTmuNonZ_3lOSM12Cut_3mu", muonColl.at(IdxNonZSS).Pt(), weight, 0., 200., 200);
       FillHist("EtamuNonZ_3lOSM12Cut_3mu", muonColl.at(IdxNonZSS).Eta(), weight, -5., 5., 100);
       float dRS=min(muonColl.at(IdxNonZSS).DeltaR(muonColl.at(IdxZmuOS)),muonColl.at(IdxNonZSS).DeltaR(muonColl.at(IdxZmuSS)));
       float dRL=max(muonColl.at(IdxNonZSS).DeltaR(muonColl.at(IdxZmuOS)),muonColl.at(IdxNonZSS).DeltaR(muonColl.at(IdxZmuSS)));
       float dRZmu=muonColl.at(IdxNonZSS).DeltaR(muonColl.at(IdxZmuOS)+muonColl.at(IdxZmuSS));
       float dPhiS=min(fabs(muonColl.at(IdxNonZSS).DeltaPhi(muonColl.at(IdxZmuOS))),fabs(muonColl.at(IdxNonZSS).DeltaPhi(muonColl.at(IdxZmuSS))));
       float dPhiL=max(fabs(muonColl.at(IdxNonZSS).DeltaPhi(muonColl.at(IdxZmuOS))),fabs(muonColl.at(IdxNonZSS).DeltaPhi(muonColl.at(IdxZmuSS))));
       float dPhiZmu=fabs(muonColl.at(IdxNonZSS).DeltaPhi(muonColl.at(IdxZmuOS)+muonColl.at(IdxZmuSS)));

       FillHist("dRmuZmuNonZS_3lOSM12Cut_3mu", dRS, weight, 0., 5., 500);
       FillHist("dRmuZmuNonZL_3lOSM12Cut_3mu", dRL, weight, 0., 5., 500);
       FillHist("dRmuZ_3lOSM12Cut_3mu", dRZmu, weight, 0., 5., 500);
       FillHist("AbsdPhimuZmuNonZS_3lOSM12Cut_3mu", dPhiS, weight, 0., 3.15, 200);
       FillHist("AbsdPhimuZmuNonZL_3lOSM12Cut_3mu", dPhiL, weight, 0., 3.15, 200);
       FillHist("AbsdPhimuZ_3lOSM12Cut_3mu", dPhiZmu, weight, 0., 3.15, 200);

       //M(3l)~90 events inspection
       if(M3l>70 && M3l<100){
         FillHist("dRmuZmuNonZS_3lOSM12M3lCut_3mu", dRS, weight, 0., 5., 500);
         FillHist("dRmuZmuNonZL_3lOSM12M3lCut_3mu", dRL, weight, 0., 5., 500);
         FillHist("dRmuZ_3lOSM12M3lCut_3mu", dRZmu, weight, 0., 5., 500);
         FillHist("AbsdPhimuZmuNonZS_3lOSM12M3lCut_3mu", dPhiS, weight, 0., 3.15, 200);
         FillHist("AbsdPhimuZmuNonZL_3lOSM12M3lCut_3mu", dPhiL, weight, 0., 3.15, 200);
         FillHist("AbsdPhimuZ_3lOSM12M3lCut_3mu", dPhiZmu, weight, 0., 3.15, 200);
       }
       /*-----------------------------------------------------------------------------------------------------------------------*/
      
       //After removing 3l~90GeV peak
       if(fabs(M3l-90)>10){
         FillHist("PTmu1_3lOSM12M3l90Cut_3mu", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_3lOSM12M3l90Cut_3mu", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("PTmu3_3lOSM12M3l90Cut_3mu", muonColl.at(2).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_3lOSM12M3l90Cut_3mu", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_3lOSM12M3l90Cut_3mu", muonColl.at(1).Eta(), weight, -5., 5., 100);
         FillHist("Etamu3_3lOSM12M3l90Cut_3mu", muonColl.at(2).Eta(), weight, -5., 5., 100);
    
         FillHist("Nj_3lOSM12M3l90Cut_3mu", jetColl.size(), weight, 0., 10., 10);
         FillHist("NlVeto_3lOSM12M3l90Cut_3mu", muonVetoColl.size()+electronVetoColl.size(), weight, 0., 10., 10);
         FillHist("NeVeto_3lOSM12M3l90Cut_3mu", electronVetoColl.size(), weight, 0., 10., 10);
         FillHist("NmuVeto_3lOSM12M3l90Cut_3mu", muonVetoColl.size(), weight, 0., 10., 10);
         FillHist("MET_3lOSM12M3l90Cut_3mu", met, weight, 0., 200., 200); 
         FillHist("Nb_3lOSM12M3l90Cut_3mu", bjetColl.size(), weight, 0., 10., 10);
         FillHist("dRmumu_OSSS1_3lOSM12M3l90Cut_3mu", muonColl.at(IdxOS).DeltaR(muonColl.at(IdxSS1)), weight, 0., 5., 500);
         FillHist("dRmumu_OSSS2_3lOSM12M3l90Cut_3mu", muonColl.at(IdxOS).DeltaR(muonColl.at(IdxSS2)), weight, 0., 5., 500);
         FillHist("Mmumu_OSSS1_3lOSM12M3l90Cut_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS1)).M(), weight, 0., 200., 200);
         FillHist("Mmumu_OSSS2_3lOSM12M3l90Cut_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS2)).M(), weight, 0., 200., 200);
         FillHist("M3mu_3lOSM12M3l90Cut_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS1)+muonColl.at(IdxSS2)).M(), weight, 0., 500., 500);
       }

     }

     //ttZ peak test : ttZ cross section && Fake Control with b present
//     Reg_3ljb: ;

   }
   if(DiMuon_analysis){

     //Step1 : 2OSmu +PTCut +Zwindow
     //if( !(muonColl.at(0).Pt()>27) ) return; //SingleMuon Trig Case
     if( !(muonColl.size()==2) ) return;
     if( !(muonColl.at(0).Pt()>20 && muonColl.at(1).Pt()>10) ) return;
     if( fabs(SumCharge(muonColl))!=0 ) return;
     if( fabs((muonColl.at(0)+muonColl.at(1)).M()-91.2)>15 )   return;

   }
   if(DiEle_analysis){

     //Step1 : 1e+2mu +PTCut
     //if( !(muonColl.at(0).Pt()>27) ) return; //SingleMuon Trig Case
     if( !(electronColl.size()==2) ) return;
     if( !(electronColl.at(0).Pt()>25 && electronColl.at(1).Pt()>15) ) return;
     if( electronColl.at(0).Charge() == electronColl.at(1).Charge() )  return;
     if( fabs((electronColl.at(0)+electronColl.at(1)).M()-91.2)>15 )   return;


   }



/////////////////////////////////////////////////////////////////////////////////// 


return;
}// End of execute event loop
  


void Mar2017_3l4j_TriLepComp::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Mar2017_3l4j_TriLepComp::BeginCycle() throw( LQError ){
  
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

Mar2017_3l4j_TriLepComp::~Mar2017_3l4j_TriLepComp() {
  
  Message("In Mar2017_3l4j_TriLepComp Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Mar2017_3l4j_TriLepComp::FillCutFlow(TString cut, float weight){
  
  if(GetHist("cutflow_W") && GetHist("cutflow_N")){
    GetHist("cutflow_W")->Fill(cut,weight);
    GetHist("cutflow_N")->Fill(cut,1);
  }
  else{
    if(!GetHist("cutflow_W")){
      AnalyzerCore::MakeHistograms("cutflow_W", 6, 0., 6.);
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(1,"NoCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(2,"TriggerCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(3,"EventCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(4,"VertexCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(5,"3lCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(6,"bVeto");
    }
    if(!GetHist("cutflow_N")){
      AnalyzerCore::MakeHistograms("cutflow_N", 6, 0., 6.);
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(1,"NoCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(2,"TriggerCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(3,"EventCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(4,"VertexCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(5,"3lCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(6,"bVeto");
    }
  }
}



void Mar2017_3l4j_TriLepComp::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Mar2017_3l4j_TriLepComp::MakeHistograms(){
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


  AnalyzerCore::MakeHistograms("Basic_Nj_wNlOScut", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Nb_wNlOScut", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Ne_wNlOScut", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Nmu_wNlOScut", 10, 0., 10.);

  AnalyzerCore::MakeHistograms("Basic_Pte_wNlOScut", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu1_wNlOScut", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu2_wNlOScut", 200, 0., 200.);

  AnalyzerCore::MakeHistograms("Basic_Mmumu_wNlOScut", 200, 0., 200.);


  //After Nljcut
  AnalyzerCore::MakeHistograms("Basic_Nb_wNljcut", 10, 0., 10.);


  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Mar2017_3l4j_TriLepCompCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Mar2017_3l4j_TriLepComp::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
