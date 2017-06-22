// $Id: Mar2017_Validation.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQMar2017_Validation Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Mar2017_Validation.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Mar2017_Validation);

 Mar2017_Validation::Mar2017_Validation() : AnalyzerCore(), out_muons(0) {

   SetLogName("Mar2017_Validation");
   Message("In Mar2017_Validation constructor", INFO);
   InitialiseAnalysis();
 }


 void Mar2017_Validation::InitialiseAnalysis() throw( LQError ) {
   
   /// Initialise histograms
   MakeHistograms();  
   // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
   // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
   Message("HwA analysis", INFO);
   return;
 }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Loop///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Mar2017_Validation::ExecuteEvents()throw( LQError ){

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
   //if (!k_isdata) { pileup_reweight=eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::down); }
   if (!k_isdata) { pileup_reweight=eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central); }

   //Numbet of Vertex PUreweight
   FillHist("Nvtx_nocut_PURW", eventbase->GetEvent().nVertices(), weight*pileup_reweight, 0., 50., 50);


   //Total Event(MCweight+PUrw)////////////////
   FillCutFlow("NoCut", weight*pileup_reweight);


   bool DoubleEle_analysis=false, EMu_analysis=false, SingleMu_analysis=false, DoubleMu_analysis=false, SingleEle_analysis=false;
   bool SiglMuTrig=false, SiglEleTrig=false, EMuTrig=false, DiMuTrig=false, DiEleTrig=false;
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("DoubleEle_analysis")) DoubleEle_analysis=true;
     else if(k_flags.at(i).Contains("EMu_analysis"))       EMu_analysis      =true;
     else if(k_flags.at(i).Contains("SingleMu_analysis"))  SingleMu_analysis =true;
     else if(k_flags.at(i).Contains("SingleEle_analysis")) SingleEle_analysis=true;
     else if(k_flags.at(i).Contains("DoubleMu_analysis"))  DoubleMu_analysis =true;
     else if(k_flags.at(i).Contains("SiglMuTrig"))         SiglMuTrig        =true;
     else if(k_flags.at(i).Contains("SiglEleTrig"))        SiglEleTrig       =true;
     else if(k_flags.at(i).Contains("EMuTrig"))            EMuTrig           =true;
     else if(k_flags.at(i).Contains("DiMuTrig"))           DiMuTrig          =true;
     else if(k_flags.at(i).Contains("DiEleTrig"))          DiEleTrig         =true;
   }


    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //ListTriggersAvailable();


   //Trigger Path of Analysis & Lumi Normalisation Setting
   bool Pass_Trigger=false;
   float trigger_ps_weight=1.;

   if(DoubleMu_analysis){
     if(DiMuTrig){
       if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")
          || PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") ) Pass_Trigger=true;
     }
     if(SiglMuTrig){
       if( PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v") ) Pass_Trigger=true;
     }
     
     if(!isData) trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
   }
   else if(SingleMu_analysis){
     if( PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v") ) Pass_Trigger=true;
     if(!isData) trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
   }
   else if(SingleEle_analysis){
     if(SiglEleTrig){
       if( PassTrigger("HLT_Ele27_WPTight_Gsf_v") ) Pass_Trigger=true;
       if(!isData) trigger_ps_weight=WeightByTrigger("HLT_Ele27_WPTight_Gsf_v", TargetLumi);
     }
   }
   else if(DoubleEle_analysis){
     if(DiEleTrig){
       if( PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger=true;
       if(!isData) trigger_ps_weight=WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);
     }
     if(SiglEleTrig){
       if( PassTrigger("HLT_Ele27_WPTight_Gsf_v") ) Pass_Trigger=true;
       if(!isData) trigger_ps_weight=WeightByTrigger("HLT_Ele27_WPTight_Gsf_v", TargetLumi);
     }
   }
   else if(EMu_analysis){
     if(SiglMuTrig){
       if( PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v") ) Pass_Trigger=true;
     }
     if(EMuTrig){
       if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
          || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger=true;
     }
     if(!isData) trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
   }
   weight*=trigger_ps_weight;
   FillHist("TriggerPSWeight", trigger_ps_weight, 1., 0., 1., 100);

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
  
   
   //POG IDs
     //eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     //eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     //eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.25);//POG WP L
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.05);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);//(muonColl, bool RochCorr, bool debug)
     //eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     //eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     //eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.15);//POG WP L
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(5.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.05);               eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");  eventbase->GetMuonSel()->SetRelIso(0.1);
   std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,false);

   std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->Selection(muonColl, true);//(muonColl, bool RochCorr, bool debug)

   
   //HNIDs
   /*  eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.25);//POG WP L
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);//(muonColl, bool RochCorr, bool debug)
   
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.05);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);//POG WP L
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);//(muonColl, bool RochCorr, bool debug)
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.05);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.1);//POG WP T
   std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,true);//(muonTightColl, bool RochCorr, bool debug)
   std::vector<snu::KMuon> muonColl;      if(k_running_nonprompt){ muonColl=muonLooseColl;} else{ muonColl=muonTightColl;}*/
   //std::vector<snu::KMuon> muonPromptColl; muonPromptColl=GetTruePrompt(muonColl, false);

     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.5, 0.5);//2016 80X tuned WP
     //eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0994, 0.107);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.2);//Not in ID, but additional safe WP
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.1);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetdxySigMax(3.);
   std::vector<snu::KElectron> electronColl;      eventbase->GetElectronSel()->Selection(electronColl);

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     //eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     eventbase->GetJetSel()->SetPt(30.);                     eventbase->GetJetSel()->SetEta(2.4);//ForFakeStudyValidation
     //eventbase->GetJetSel()->SetPileUpJetID(true,"Loose");
     bool LeptonVeto=true;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonLooseColl, electronLooseColl);



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
   double metphi =eventbase->GetEvent().METPhi();
   double Pzv,Pzv1, Pzv2;
   //cout<<"METallcorr "<<eventbase->GetEvent().MET()<<" METtype1 "<<eventbase->GetEvent().PFMETType1()<<endl;
   //cout<<"METtype1x "<<eventbase->GetEvent().PFMETType1x()<<" METtype1y "<<eventbase->GetEvent().PFMETType1y()<<endl;
   //snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));
   //snu::KParticle v[4]; v[0].SetPx(met_x); v[0].SetPy(met_y);
   //                     v[1].SetPx(met_x); v[1].SetPy(met_y);

   int   nbjets=bjetColl.size(); int njets=jetColl.size(); const int nljets=ljetColl.size();
   int   Nvtx=eventbase->GetEvent().nVertices();
   float Rho=eventbase->GetEvent().Rho();


   /*****************************************************
   **Scale Factors
   *****************************************************/
   float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1., btag_sf=1.;
   float trigger_sf=1.;
   float top_pt_reweight=1.;
   float fake_weight=1.; bool EventCand=false;

   /*This part is for boosting up speed.. SF part takes rather longer time than expected*/
   if     (DoubleMu_analysis) { if(muonLooseColl.size()>=2)                                EventCand=true; }
   else if(SingleMu_analysis) { if(muonLooseColl.size()>=2)                                EventCand=true; }
   else if(SingleEle_analysis){ if(electronLooseColl.size()>=1)                            EventCand=true; }
   else if(DoubleEle_analysis){ if(electronLooseColl.size()>=2)                            EventCand=true; }
   else if(EMu_analysis)      { if(electronLooseColl.size()>=1 && muonLooseColl.size()>=1) EventCand=true; }

   if(EventCand){
     if(!isData){
       if(SiglMuTrig)  trigger_sf = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );
       if(SiglEleTrig) trigger_sf = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_Ele27_WPTight_Gsf" );
       //double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electronColl, muonColl, 0, 0, 0);
       //double trigger_eff_MC   = mcdata_correction->TriggerEfficiencyLegByLeg(electronColl, muonColl, 0, 1, 0);
       //trigger_sf = trigger_eff_Data/trigger_eff_MC;
  
       id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronColl);
       reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
  
       //id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muonColl);
       id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muonColl);
       iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
       trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
  
       btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

       if(k_sample_name.Contains("TT_powheg")){
         top_pt_reweight = TopPTReweight(truthColl);
       }
     }
     else{
       //Perfect Prompt Ratio Approximation applied.(p=1), Applicable to generic number, combination of leptons under premise of high prompt rate.
       if(k_running_nonprompt){
         //fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", muonLooseColl.size(), electronNull, "ELECTRON_HN_TIGHT", electronNull.size());
       }
     }
   }

   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight*fake_weight*btag_sf*top_pt_reweight*trigger_sf;
   /***************************************************************************************************/

   //////Basic Objects Check//////////////////////////////////////////////////////////
   FillHist("Basic_Nj_orig", njets, weight, 0., 10., 10);
   FillHist("Basic_Nb_orig", nbjets, weight, 0., 10., 10);//Bjet def CVSincV2>0.89;pfKJet::CSVv2IVF_discriminator 
   FillHist("Basic_NmuT_orig", muonColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NmuL_orig", muonLooseColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeT_orig", electronColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeL_orig", electronLooseColl.size(), weight, 0., 10., 10);
   if(electronColl.size()==1) FillHist("Basic_NmuT_weT", muonColl.size(), weight, 0., 10., 10);
   if(muonColl.size()>0) FillHist("Basic_Ptmu1_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
   if(muonColl.size()>1) FillHist("Basic_Ptmu2_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
   if(muonColl.size()>2) FillHist("Basic_Ptmu3_orig", muonColl.at(2).Pt(), weight, 0., 200., 200);
   if(electronColl.size()>0) FillHist("Basic_Pte_orig", electronColl.at(0).Pt(), weight, 0., 200., 200);
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


   if(EMu_analysis){


     //Step1 : 1e+1mu
     if( !(electronColl.size()==1 && muonColl.size()==1) ) return;
     if( EMuTrig    && !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>27) ) return;
     if( SiglMuTrig && !(muonColl.at(0).Pt()>27 && electronColl.at(0).Pt()>25) ) return;
     FillHist("Cutflow_TrkRecoIDIsoTrigPUBW", 0., weight, 0., 10., 10);


     //Step2 : OS
     if(electronColl.at(0).Charge()==muonColl.at(0).Charge()) return;
     FillHist("Cutflow_TrkRecoIDIsoTrigPUBW", 1., weight, 0., 10., 10);

     FillHist("Nj_NlOSCut_TrkRecoIDIsoTrigPUBW", njets, weight, 0., 10., 10);
     FillHist("MET_NlOSCut_TrkRecoIDIsoTrigPUBW", met, weight, 0., 500., 500);



     //Step3 : Nj>=1
     if(jetColl.size()<1) return;
     FillHist("Cutflow_TrkRecoIDIsoTrigPUBW", 2., weight, 0., 10., 10);



     //Step4 : Nj>=2
     if(jetColl.size()<2) return;
     FillHist("Cutflow_TrkRecoIDIsoTrigPUBW", 3., weight, 0., 10., 10);

     FillHist("MET_NljOSCut_TrkRecoIDIsoTrigPUBW", met, weight, 0., 500., 500);
     FillHist("Nb_NljOSCut_TrkRecoIDIsoTrigPUBW", bjetColl.size(), weight, 0., 10., 10);


     FillHist("PTe_NljOSCut_TrkRecoIDIsoTrigPUBW", electronColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("Etae_NljOSCut_TrkRecoIDIsoTrigPUBW", electronColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("PTmu1_NljOSCut_TrkRecoIDIsoTrigPUBW", muonColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("Etamu1_NljOSCut_TrkRecoIDIsoTrigPUBW", muonColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("PTj1_NljOSCut_TrkRecoIDIsoTrigPUBW", jetColl.at(0).Pt(), weight, 0., 500., 500);
     FillHist("Etaj1_NljOSCut_TrkRecoIDIsoTrigPUBW", jetColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("PTj2_NljOSCut_TrkRecoIDIsoTrigPUBW", jetColl.at(1).Pt(), weight, 0., 500., 500);
     FillHist("Etaj2_NljOSCut_TrkRecoIDIsoTrigPUBW", jetColl.at(1).Eta(), weight, -5., 5., 100);
     FillHist("Nvtx_NljOSCut_TrkRecoIDIsoTrigPUBW", Nvtx, weight, 0., 50., 50);
     FillHist("dRemu_NljOSCut_TrkRecoIDIsoTrigPUBW", electronColl.at(0).DeltaR(muonColl.at(0)), weight, 0., 5., 100);
     FillHist("dPhiemu_NljOSCut_TrkRecoIDIsoTrigPUBW", electronColl.at(0).DeltaPhi(muonColl.at(0)), weight, -3.15, 3.15, 200);



     //Step5 : Nb>=1
     if(bjetColl.size()==0) return;
     FillHist("Cutflow_TrkRecoIDIsoTrigPUBW", 4., weight, 0., 10., 10);


     FillHist("PTe_NljbOSCut_TrkRecoIDIsoTrigPUBW", electronColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("Etae_NljbOSCut_TrkRecoIDIsoTrigPUBW", electronColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("PTmu1_NljbOSCut_TrkRecoIDIsoTrigPUBW", muonColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("Etamu1_NljbOSCut_TrkRecoIDIsoTrigPUBW", muonColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("PTj1_NljbOSCut_TrkRecoIDIsoTrigPUBW", jetColl.at(0).Pt(), weight, 0., 500., 500);
     FillHist("Etaj1_NljbOSCut_TrkRecoIDIsoTrigPUBW", jetColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("PTj2_NljbOSCut_TrkRecoIDIsoTrigPUBW", jetColl.at(1).Pt(), weight, 0., 500., 500);
     FillHist("Etaj2_NljbOSCut_TrkRecoIDIsoTrigPUBW", jetColl.at(1).Eta(), weight, -5., 5., 100);
     FillHist("PTb1_NljbOSCut_TrkRecoIDIsoTrigPUBW", bjetColl.at(0).Pt(), weight, 0., 500., 500);
     FillHist("Etab1_NljbOSCut_TrkRecoIDIsoTrigPUBW", bjetColl.at(0).Eta(), weight, -5., 5., 100);

     FillHist("MET_NljbOSCut_TrkRecoIDIsoTrigPUBW", met, weight, 0., 500., 500);
     FillHist("METphi_NljbOSCut_TrkRecoIDIsoTrigPUBW", metphi, weight, -3.15, 3.15, 100);
     FillHist("Nvtx_NljbOSCut_TrkRecoIDIsoTrigPUBW", Nvtx, weight, 0., 50., 50);
     FillHist("Rho_NljbOSCut_TrkRecoIDIsoTrigPUBW", Rho, weight, 0., 50., 50);
     FillHist("Nj_NljbOSCut_TrkRecoIDIsoTrigPUBW", njets, weight, 0., 10., 10);
     FillHist("Nb_NljbOSCut_TrkRecoIDIsoTrigPUBW", bjetColl.size(), weight, 0., 10., 10);
     FillHist("dRemu_NljbOSCut_TrkRecoIDIsoTrigPUBW", electronColl.at(0).DeltaR(muonColl.at(0)), weight, 0., 5., 100);
     FillHist("dPhiemu_NljbOSCut_TrkRecoIDIsoTrigPUBW", electronColl.at(0).DeltaPhi(muonColl.at(0)), weight, -3.15, 3.15, 200);


     //Step6 : MET>40
     if(met<40) return;
     FillHist("Cutflow_TrkRecoIDIsoTrigPUBW", 5., weight, 0., 10., 10);

     //Step7 : Nj>=3
     if(jetColl.size()<3) return;
     FillHist("Cutflow_TrkRecoIDIsoTrigPUBW", 6., weight, 0., 10., 10);


   }
   if(SingleMu_analysis){

     if( !(PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v")) ) return;
     if( !(muonColl.size()==2) )     return;
     if( !(electronColl.size()==0) )  return;
     if( SumCharge(muonColl)!=0 ) return;
     if( !(muonColl.at(0).Pt()>27 && muonColl.at(1).Pt()>27) )  return;
     if( !(fabs((muonColl.at(0)+muonColl.at(1)).M()-91.2)<15) ) return;
    
     float weightSFapplied = weight*trigger_sf*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight;


   
   }
   if(DoubleMu_analysis){


     if( muonColl.size()!=2 ) return;
     if( SumCharge(muonColl)!=0 ) return;
     if( DiMuTrig   && !(muonColl.at(0).Pt()>20 && muonColl.at(1).Pt()>20) ) return; //DoubleMu Trig
     if( SiglMuTrig && !(muonColl.at(0).Pt()>27 && muonColl.at(1).Pt()>20) ) return; //SignleMu Trig
     if( fabs((muonColl.at(0)+muonColl.at(1)).M()-91.2)>15 ) return; //Only Z peak
     //if( !((muonColl.at(0)+muonColl.at(1)).M()>12) ) return; //Whole DiMuon Region


     //Inclusive Plots
     FillHist("YieldComparison_TkIDIsoTrigPUW", 0., weight, 0., 10., 10);
     FillHist("PTmu1_TkIDIsoTrigPUW", muonColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("PTmu2_TkIDIsoTrigPUW", muonColl.at(1).Pt(), weight, 0., 200., 200);
     FillHist("Etamu1_TkIDIsoTrigPUW", muonColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etamu2_TkIDIsoTrigPUW", muonColl.at(1).Eta(), weight, -5., 5., 100);
     FillHist("Mmumu_TkIDIsoTrigPUW", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 200., 200);
     FillHist("Nj_TkIDIsoTrigPUW", jetColl.size(), weight, 0., 10., 10);
     FillHist("MET_TkIDIsoTrigPUW", met, weight, 0., 500., 500);
     FillHist("Nvtx_TkIDIsoTrigPUW", Nvtx, weight, 0., 50., 50);
     FillHist("Rho_TkIDIsoTrigPUW", Rho, weight, 0., 50., 50);
     if(pileup_reweight!=0) FillHist("Rho_TkIDIsoTrigW", Rho, weight/pileup_reweight, 0., 50., 50);

     for(int i=0; i<muonColl.size(); i++){     
       FillHist("RelIso04_TkIDIsoTrigPUW", muonColl.at(i).RelIso04()*muonColl.at(i).Pt()/muonColl.at(i).RochPt(), weight, 0., 0.2, 200);
       FillHist("Absd0_TkIDIsoTrigPUW", fabs(muonColl.at(i).dXY()), weight, 0., 0.2, 2000);
       FillHist("Absd0sig_TkIDIsoTrigPUW", fabs(muonColl.at(i).dXYSig()), weight, 0., 10., 100);
       FillHist("Absdz_TkIDIsoTrigPUW", fabs(muonColl.at(i).dZ()), weight, 0., 0.5, 5000);
       FillHist("NhitTkLayer_TkIDIsoTrigPUW", muonColl.at(i).ActiveLayer(), weight, 0., 20., 20);//Nhit on trackLayer
       FillHist("Nhitpixel_TkIDIsoTrigPUW", muonColl.at(i).validPixHits(), weight, 0., 20., 20);//Nhit on pixel
       FillHist("NhitMuChamber_TkIDIsoTrigPUW", muonColl.at(i).validHits(), weight, 0., 50., 50);//Nhit on muon chamber
       FillHist("NmatchedStation_TkIDIsoTrigPUW", muonColl.at(i).validStations(), weight, 0., 20., 20);//N muon station having muon segment
     }
     
     
     //Per Jet Bins
     if(njets==0){
       FillHist("YieldComparison_TkIDIsoTrigPUW", 1., weight, 0., 10., 10);
       FillHist("PTmu1_0j_TkIDIsoTrigPUW", muonColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTmu2_0j_TkIDIsoTrigPUW", muonColl.at(1).Pt(), weight, 0., 200., 200);
       FillHist("Etamu1_0j_TkIDIsoTrigPUW", muonColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Etamu2_0j_TkIDIsoTrigPUW", muonColl.at(1).Eta(), weight, -5., 5., 100);
       FillHist("MET_0j_TkIDIsoTrigPUW", met, weight, 0., 500., 500);
     }
     if(njets==1){
       FillHist("YieldComparison_TkIDIsoTrigPUW", 2., weight, 0., 10., 10);
       FillHist("PTmu1_1j_TkIDIsoTrigPUW", muonColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTmu2_1j_TkIDIsoTrigPUW", muonColl.at(1).Pt(), weight, 0., 200., 200);
       FillHist("Etamu1_1j_TkIDIsoTrigPUW", muonColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Etamu2_1j_TkIDIsoTrigPUW", muonColl.at(1).Eta(), weight, -5., 5., 100);
       FillHist("MET_1j_TkIDIsoTrigPUW", met, weight, 0., 500., 500);
       FillHist("PTj1_1j_TkIDIsoTrigPUW", jetColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("Etaj1_1j_TkIDIsoTrigPUW", jetColl.at(0).Eta(), weight, -5., 5., 100);
     }
     if(njets==2){
       FillHist("YieldComparison_TkIDIsoTrigPUW", 3., weight, 0., 10., 10);
       FillHist("PTmu1_2j_TkIDIsoTrigPUW", muonColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTmu2_2j_TkIDIsoTrigPUW", muonColl.at(1).Pt(), weight, 0., 200., 200);
       FillHist("Etamu1_2j_TkIDIsoTrigPUW", muonColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Etamu2_2j_TkIDIsoTrigPUW", muonColl.at(1).Eta(), weight, -5., 5., 100);
       FillHist("MET_2j_TkIDIsoTrigPUW", met, weight, 0., 500., 500);
       FillHist("PTj1_2j_TkIDIsoTrigPUW", jetColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTj2_2j_TkIDIsoTrigPUW", jetColl.at(1).Pt(), weight, 0., 200., 200);
       FillHist("Etaj1_2j_TkIDIsoTrigPUW", jetColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Etaj2_2j_TkIDIsoTrigPUW", jetColl.at(1).Eta(), weight, -5., 5., 100);
     }


   }
   if(DoubleEle_analysis){

     if( electronColl.size()!=2 ) return;
     if( electronColl.at(0).Charge()==electronColl.at(1).Charge() )    return;
     if( SiglEleTrig && !(electronColl.at(0).Pt()>30. && electronColl.at(1).Pt()>30.) ) return;
     if( DiEleTrig   && !(electronColl.at(0).Pt()>25. && electronColl.at(1).Pt()>20.) ) return;
     if( SiglEleTrig && !(electronColl.at(0).TriggerMatched("HLT_Ele27_WPTight_Gsf_v")
                         && electronColl.at(1).TriggerMatched("HLT_Ele27_WPTight_Gsf_v")) ) return;
     if( fabs((electronColl.at(0)+electronColl.at(1)).M()-91.2)>15 )   return;
     
      
         
     //Inclusive Plots
     FillHist("YieldComparison_RecoIDTrigPUW", 0., weight, 0., 10., 10);
     FillHist("PTele1_RecoIDTrigPUW", electronColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("PTele2_RecoIDTrigPUW", electronColl.at(1).Pt(), weight, 0., 200., 200);
     FillHist("Etaele1_RecoIDTrigPUW", electronColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etaele2_RecoIDTrigPUW", electronColl.at(1).Eta(), weight, -5., 5., 100);
     FillHist("Mee_RecoIDTrigPUW", (electronColl.at(0)+electronColl.at(1)).M(), weight, 0., 200., 200);
     FillHist("Nj_RecoIDTrigPUW", jetColl.size(), weight, 0., 10., 10);
     FillHist("MET_RecoIDTrigPUW", met, weight, 0., 500., 500);
     FillHist("Nvtx_RecoIDTrigPUW", Nvtx, weight, 0., 50., 50);

     //Inclusive Per Barrel & Endcap
     for(int i=0; i<electronColl.size(); i++){     
       if(fabs(electronColl.at(i).Eta())<1.479){
         FillHist("Absd0_eB_RecoIDTrigPUW", fabs(electronColl.at(i).dxy()), weight, 0., 0.05, 500);
         FillHist("Absdz_eB_RecoIDTrigPUW", fabs(electronColl.at(i).dz()), weight, 0., 0.1, 1000);
         FillHist("RelIso03_eB_RecoIDTrigPUW", electronColl.at(i).PFRelIso(0.3), weight, 0., 0.1, 100);
       }
       else{
         FillHist("Absd0_eE_RecoIDTrigPUW", fabs(electronColl.at(i).dxy()), weight, 0., 0.05, 500);
         FillHist("Absdz_eE_RecoIDTrigPUW", fabs(electronColl.at(i).dz()), weight, 0., 0.1, 1000);
         FillHist("RelIso03_eE_RecoIDTrigPUW", electronColl.at(i).PFRelIso(0.3), weight, 0., 0.1, 100);
       }
     }
     if     ( fabs(electronColl.at(0).Eta())<1.479 && fabs(electronColl.at(1).Eta())<1.479 ){
       FillHist("MeBeB_RecoIDTrigPUW", (electronColl.at(0)+electronColl.at(1)).M(), weight, 0., 200., 200);
     }
     else if( fabs(electronColl.at(0).Eta())>1.479 && fabs(electronColl.at(1).Eta())>1.479 ){
       FillHist("MeEeE_RecoIDTrigPUW", (electronColl.at(0)+electronColl.at(1)).M(), weight, 0., 200., 200);
     }
     

     //Per Jet Bins
     if(njets==0){
       FillHist("YieldComparison_RecoIDTrigPUW", 1., weight, 0., 10., 10);
       FillHist("PTele1_0j_RecoIDTrigPUW", electronColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTele2_0j_RecoIDTrigPUW", electronColl.at(1).Pt(), weight, 0., 200., 200);
       FillHist("Etaele1_0j_RecoIDTrigPUW", electronColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Etaele2_0j_RecoIDTrigPUW", electronColl.at(1).Eta(), weight, -5., 5., 100);
       FillHist("MET_0j_RecoIDTrigPUW", met, weight, 0., 500., 500);
     }
     if(njets==1){
       FillHist("YieldComparison_RecoIDTrigPUW", 2., weight, 0., 10., 10);
       FillHist("PTele1_1j_RecoIDTrigPUW", electronColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTele2_1j_RecoIDTrigPUW", electronColl.at(1).Pt(), weight, 0., 200., 200);
       FillHist("Etaele1_1j_RecoIDTrigPUW", electronColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Etaele2_1j_RecoIDTrigPUW", electronColl.at(1).Eta(), weight, -5., 5., 100);
       FillHist("MET_1j_RecoIDTrigPUW", met, weight, 0., 500., 500);
       FillHist("PTj1_1j_RecoIDTrigPUW", jetColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("Etaj1_1j_RecoIDTrigPUW", jetColl.at(0).Eta(), weight, -5., 5., 100);
     }
     if(njets==2){
       FillHist("YieldComparison_RecoIDTrigPUW", 3., weight, 0., 10., 10);
       FillHist("PTele1_2j_RecoIDTrigPUW", electronColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTele2_2j_RecoIDTrigPUW", electronColl.at(1).Pt(), weight, 0., 200., 200);
       FillHist("Etaele1_2j_RecoIDTrigPUW", electronColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Etaele2_2j_RecoIDTrigPUW", electronColl.at(1).Eta(), weight, -5., 5., 100);
       FillHist("MET_2j_RecoIDTrigPUW", met, weight, 0., 500., 500);
       FillHist("PTj1_2j_RecoIDTrigPUW", jetColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTj2_2j_RecoIDTrigPUW", jetColl.at(1).Pt(), weight, 0., 200., 200);
       FillHist("Etaj1_2j_RecoIDTrigPUW", jetColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Etaj2_2j_RecoIDTrigPUW", jetColl.at(1).Eta(), weight, -5., 5., 100);
     }



   }//DoubleEle End
   if(SingleEle_analysis){
     if( electronLooseColl.size()!=1 ) return;
     if( electronColl.size()!=1 ) return;
     if( SiglEleTrig && !(electronColl.at(0).Pt()>30.) ) return;

     float MTW = sqrt(2)*sqrt(met*electronColl.at(0).Pt()-met_x*electronColl.at(0).Px()-met_y*electronColl.at(0).Py());

     FillHist("MET_e27_e30", met, weight, 0., 300., 300);
     if(met>30) FillHist("MTW_e27_e30met30", MTW, weight, 0., 200., 200);
     if(jetColl.size()>0){
       FillHist("MET_e27_e30gt1j", met, weight, 0., 300., 300);
       if(met>30) FillHist("MTW_e27_e30met30gt1j", MTW, weight, 0., 200., 200);
     }
     if(met>30 && MTW>70){
       FillHist("Count_NormCR", 0., weight, 0., 5., 5);
       FillHist("MTW_e27_e30met30mtw70", MTW, weight, 50., 200., 150);
       if(jetColl.size()==1){
          if(jetColl.at(0).Pt()>30 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0){
            FillHist("Count_NormCR", 1., weight, 0., 5., 5);
            FillHist("MTW_e27_e30met30mtw701jdR01", MTW, weight, 50., 200., 150);
          }
       }
       if(jetColl.size()>0){
          FillHist("Count_NormCR", 2., weight, 0., 5., 5);
       }
     }

   }



/////////////////////////////////////////////////////////////////////////////////// 


return;
}// End of execute event loop
  


void Mar2017_Validation::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Mar2017_Validation::BeginCycle() throw( LQError ){
  
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

Mar2017_Validation::~Mar2017_Validation() {
  
  Message("In Mar2017_Validation Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Mar2017_Validation::FillCutFlow(TString cut, float weight){
  
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



void Mar2017_Validation::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Mar2017_Validation::MakeHistograms(){
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


  AnalyzerCore::MakeHistograms("YieldComparison_RecoIDPUW", 10, 0., 10.);
    GetHist("YieldComparison_RecoIDPUW")->GetXaxis()->SetBinLabel(1,"Incl");
    GetHist("YieldComparison_RecoIDPUW")->GetXaxis()->SetBinLabel(2,"0j");
    GetHist("YieldComparison_RecoIDPUW")->GetXaxis()->SetBinLabel(3,"1j");
    GetHist("YieldComparison_RecoIDPUW")->GetXaxis()->SetBinLabel(4,"2j");
  AnalyzerCore::MakeHistograms("YieldComparison_TkIDIsoPUW", 10, 0., 10.);
    GetHist("YieldComparison_TkIDIsoPUW")->GetXaxis()->SetBinLabel(1,"Incl");
    GetHist("YieldComparison_TkIDIsoPUW")->GetXaxis()->SetBinLabel(2,"0j");
    GetHist("YieldComparison_TkIDIsoPUW")->GetXaxis()->SetBinLabel(3,"1j");
    GetHist("YieldComparison_TkIDIsoPUW")->GetXaxis()->SetBinLabel(4,"2j");


  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Mar2017_ValidationCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Mar2017_Validation::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
