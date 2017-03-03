// $Id: Feb2017_3l4j_BTagClosure.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQFeb2017_3l4j_BTagClosure Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Feb2017_3l4j_BTagClosure.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Feb2017_3l4j_BTagClosure);

 Feb2017_3l4j_BTagClosure::Feb2017_3l4j_BTagClosure() : AnalyzerCore(), out_muons(0) {

   SetLogName("Feb2017_3l4j_BTagClosure");
   Message("In Feb2017_3l4j_BTagClosure constructor", INFO);
   InitialiseAnalysis();
 }


 void Feb2017_3l4j_BTagClosure::InitialiseAnalysis() throw( LQError ) {
   
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

void Feb2017_3l4j_BTagClosure::ExecuteEvents()throw( LQError ){

////Basic Infos///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Apply the gen weight 
   float weight_nopu=1; 
   if(!isData) weight*=MCweight;
   if(!isData) weight_nopu*=MCweight;
   FillHist("GenWeight", MCweight, 1, -10, 10, 1000);

   FillHist("test", weight, 1, -1., 1., 2);
   //Total Event  
   FillCutFlow("NoCut", weight);

   /// Acts on data to remove bad reconstructed event 
   //if(isData&& (! eventbase->GetEvent().LumiMask(lumimask))) return;

   m_logger<<DEBUG<<"RunNumber/Event Number = "<<eventbase->GetEvent().RunNumber()<<" : "<<eventbase->GetEvent().EventNumber()<<LQLogger::endmsg;
   m_logger<<DEBUG<<"isData = "<<isData<<LQLogger::endmsg;


   //Numbet of Vertex NoPUreweight
   FillHist("Nvtx_nocut_NoRW", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);

   //Pileup Reweight
   float pileup_reweight=(1.0);
   if (!k_isdata) { pileup_reweight=eventbase->GetEvent().PileUpWeight(); }

   //Numbet of Vertex PUreweight
   FillHist("Nvtx_nocut_PURW", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);


//   bool EMu_analysis=true, DoubleMu_analysis=false;
   bool DoubleEle_analysis=false, SingleMu_analysis=false, DoubleMu_analysis=false, EMu_analysis=true;


   //Trigers
   std::vector<TString> triggerslist, triggerlist1, triggerlist2;//  ListTriggersAvailable();

   TString analysis_trigger;
   //if(EMu_analysis) {analysis_trigger="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";}
   if(EMu_analysis) {analysis_trigger="HLT_IsoMu24_v";}
   else if(DoubleMu_analysis) analysis_trigger="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v";
   else if(SingleMu_analysis) analysis_trigger="HLT_IsoMu24_v";
   else if(DoubleEle_analysis) analysis_trigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
   //if(!PassTrigger(analysis_trigger)) return; //Not for now..

   //PassTrigger Uses BeginWith so don't stick to exact name

   float trigger_ps_weight=1, weight_trigger_sf=1, trigger_eff=1;
/*   std::vector<snu::KMuon> trigmuColl; std::vector<snu::KElectron> trigeColl;
     eventbase->GetMuonSel()->SetPt(20.);eventbase->GetMuonSel()->SetEta(2.4);eventbase->GetMuonSel()->SetRelIso(0.15);eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);  eventbase->GetMuonSel()->Selection(trigmuColl);
     eventbase->GetElectronSel()->SetPt(20.);eventbase->GetElectronSel()->SetEta(2.4);eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MEDIUM); eventbase->GetElectronSel()->Selection(trigeColl);*/

   trigger_ps_weight=WeightByTrigger(analysis_trigger, TargetLumi);
   //weight_trigger_sf=1.;//TriggerScaleFactor(trigeColl, trigmuColl, analysis_trigger);



   FillHist("TriggerSFWeight" , weight_trigger_sf, 1., 0. , 2., 200); FillHist("TriggerPSWeight", trigger_ps_weight, 1., 0., 500., 500);
   weight*=weight_trigger_sf*trigger_ps_weight;
   weight_nopu*=weight_trigger_sf*trigger_ps_weight;
   //FillCutFlow("TriggerCut", weight); //;in trigger Cut mode
   //FillHist("Basic_prescale", prescale, 1., 0., 2000., 2000);

   //Initial Event Cut : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
   //if(!PassMETFilter()) return;  FillCutFlow("EventCut", weight);

   //Vertex Cut
   //if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; FillCutFlow("VertexCut", weight);
   //Good Primary vtx def:(vtx.ndof()>4&&maxAbsZ<=0)||std::abs(vtx.z())<= 24)&&((maxd0 <=0)||std::abs(vtx.position().rho())<=2)&&!(vtx.isFake()))  



///////Objects in Analysis/////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

   //Primary Object Collection
   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
  
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.25);//POG WP L
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);//(muonColl, bool RochCorr, bool debug)
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.15);//POG WP T
   std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->Selection(muonColl,true);//(muonColl, bool RochCorr, bool debug)
 //std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->SelectMuons(muonColl, BaseSelection::MUON_POG_TIGHT, 10., 2.4);
   std::vector<snu::KMuon> muonPromptColl; muonPromptColl=GetTruePrompt(muonColl, false);

     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
     eventbase->GetElectronSel()->SetPt(20.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0994, 0.107);//2016 80X tuned WP
   //eventbase->GetElectronSel()->SetdxyBEMax(0.005, 0.01);  eventbase->GetElectronSel()->SetdzBEMax(0.005, 0.01);//Not in ID, but additional safe WP
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(20.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.2);//Not in ID, but additional safe WP
   std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->Selection(electronColl);
 //std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->SelectElectrons(electronColl, "ELECTRON_POG_TIGHT", 20., 2.4);

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     //eventbase->GetJetSel()->SetPt(30.);                     eventbase->GetJetSel()->SetEta(2.4);
     eventbase->GetJetSel()->SetPt(20.);                     eventbase->GetJetSel()->SetEta(2.4);
   //eventbase->GetJetSel()->SetUseJetPileUp(true);
     bool LeptonVeto=true;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonColl, electronColl);
   std::vector<snu::KJet> jetLooseColl; eventbase->GetJetSel()->SelectJets(jetLooseColl, muonColl, electronColl, "PFJET_LOOSE", 20., 2.4);

   //Method to apply 2a method SF
   std::vector<int> bIdxColl  = GetSFBJetIdx(jetColl,"Medium");
   std::vector<int> ljIdxColl = GetSFLJetIdx(jetColl, bIdxColl, "Medium");
   std::vector<snu::KJet> bjetColl2a; for(int i=0; i<bIdxColl.size(); i++) {bjetColl2a.push_back(jetColl.at(bIdxColl.at(i)));}
   std::vector<snu::KJet> ljetColl2a; for(int i=0; i<ljIdxColl.size(); i++){ljetColl2a.push_back(jetColl.at(ljIdxColl.at(i)));}

   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");
   std::vector<snu::KJet> ljetColl = SelLightJets(jetColl, "Medium");

   //Temporarily Placed Cuts(Activate ones at original place if trig simul. fully available
   //METFilter
   if(!PassMETFilter()) return;  FillCutFlow("EventCut", weight);

   //Vertex Cut
   if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; FillCutFlow("VertexCut", weight);

   //Trigger eff. emul. and cuts
//   FillCutFlow("TriggerCut", weight); //;in trigger Cut mode
   /////////////////////////////////////////////////////////////////////////////////////////

   bool emu=false, mumu=false;
   int mu1_Ai=-1, mu2_Ai=-1, mu_Wi=-1;
   int nbjets=bjetColl.size(); int njets=jetColl.size(); const int nljets=ljetColl.size();//njets-nbjets;//number of light jets
   double Pzv1, Pzv2;

   double met = eventbase->GetEvent().MET();
   double met_x = eventbase->GetEvent().MET()*TMath::Cos(eventbase->GetEvent().METPhi());
   double met_y = eventbase->GetEvent().MET()*TMath::Sin(eventbase->GetEvent().METPhi());
   int Nvtx=eventbase->GetEvent().nVertices();
   double Pzv;
   snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));
//   snu::KParticle v[4]; v[0].SetPx(met_x); v[0].SetPy(met_y);
//                        v[1].SetPx(met_x); v[1].SetPy(met_y);

   //Scale Factors
   float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1.;
   float btag_sf=1.;
   float trigger_sf=1.;

   if(!isData){
    // trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );

     id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronColl);
     reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);

     id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muonColl);
     iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
     trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);

     btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);
   }


//////Basic Objects Distribution////////////////////////////////////////////////////////////////////////
   FillHist("Basic_Nj_orig", njets, weight, 0., 10., 10);
   FillHist("Basic_Nb_orig", nbjets, weight, 0., 10., 10);//Bjet def CVSincV2>0.89;pfKJet::CSVv2IVF_discriminator 
   FillHist("Basic_NmuT_orig", muonColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NmuL_orig", muonLooseColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeT_orig", electronColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeL_orig", electronLooseColl.size(), weight, 0., 10., 10);
   if(electronColl.size()==1) FillHist("Basic_NmuT_weT", muonColl.size(), weight, 0., 10., 10);
   if(DoubleMu_analysis){
     if(muonColl.size()>0) FillHist("Basic_Ptmu1_3mu_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>1) FillHist("Basic_Ptmu2_3mu_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>2) FillHist("Basic_Ptmu3_3mu_orig", muonColl.at(2).Pt(), weight, 0., 200., 200);
   }
   if(EMu_analysis){
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
   FillHist("Basic_METdist_orig", met, weight, 0., 200., 100);//It is orig because even though I havent put METcut below, if there is correlation between Nl ,Nj, Nb and MET then MET distrubution will also be affected by selection. The distribution remains the same only when the correlation is absent.

///////Event Selection&Histograms//////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////   

   if(EMu_analysis){

     //if( !(PassTrigger("HLT_IsoMu24_v")||PassTrigger("HLT_IsoTkMu24_v")) ) return;
     if( !(PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
          || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
          || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v")) ) return;


     //Step1 : 1e+1mu
     if( !(electronColl.size()==1 && muonColl.size()==1) ) return;
     if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>15) ) return; //EMu Trig Case
     //if( !(muonColl.at(0).Pt()>27) ) return; //SingleMuon Trig Case
     FillHist("Cutflow_NoW", 0., weight, 0., 10., 10);
     FillHist("Cutflow_TrkRecoW", 0., weight*trk_weight_mu*reco_weight_ele, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoW", 0., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoTrigW", 0., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoTrigPUW", 0., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, 0., 10., 10);

     //Step2 : OS
     if(electronColl.at(0).Charge()==muonColl.at(0).Charge()) return;
     FillHist("Cutflow_NoW", 1., weight, 0., 10., 10);
     FillHist("Cutflow_TrkRecoW", 1., weight*trk_weight_mu*reco_weight_ele, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoW", 1., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoTrigW", 1., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoTrigPUW", 1., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, 0., 10., 10);


     FillHist("Nj_NlOSCut_TrkRecoIDIsoTrigPUW", njets, weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, 0., 10., 10);
     FillHist("MET_NlOSCut_TrkRecoIDIsoTrigPUW", met, weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, 0., 500., 500);


     if(jetColl.size()==0){
       FillHist("Pte_0j_NlOSCut_TrkRecoIDIsoTrigPUW", electronColl.at(0).Pt(), weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, 0., 500., 500);
       FillHist("Etae_0j_NlOSCut_TrkRecoIDIsoTrigPUW", electronColl.at(0).Eta(), weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, -5., 5., 100);
       FillHist("Ptmu_0j_NlOSCut_TrkRecoIDIsoTrigPUW", muonColl.at(0).Pt(), weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, 0., 500., 500);
       FillHist("Etamu_0j_NlOSCut_TrkRecoIDIsoTrigPUW", muonColl.at(0).Eta(), weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, -5., 5., 100);
     }
     FillHist("dRemu_NlOSCut_TrkRecoIDIsoTrigPUW", electronColl.at(0).DeltaR(muonColl.at(0)), weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, 0., 5., 100);
     FillHist("dPhiemu_NlOSCut_TrkRecoIDIsoTrigPUW", electronColl.at(0).DeltaPhi(muonColl.at(0)), weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, -3.15, 3.15, 200);


     //Step3 : Nj>=1
     if(jetColl.size()<1) return;
     FillHist("Cutflow_NoW", 2., weight, 0., 10., 10);
     FillHist("Cutflow_TrkRecoW", 2., weight*trk_weight_mu*reco_weight_ele, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoW", 2., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoTrigW", 2., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoTrigPUW", 2., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, 0., 10., 10);


     //Step4 : Nj>=2
     if(jetColl.size()<2) return;
     FillHist("Cutflow_NoW", 3., weight, 0., 10., 10);
     FillHist("Cutflow_TrkRecoW", 3., weight*trk_weight_mu*reco_weight_ele, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoW", 3., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoTrigW", 3., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoTrigPUW", 3., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, 0., 10., 10);


     FillHist("MET_NljOSCut_TrkRecoIDIsoTrigPUW", met, weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, 0., 500., 500);
     FillHist("dRemu_NljOSCut_TrkRecoIDIsoTrigPUW", electronColl.at(0).DeltaR(muonColl.at(0)), weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, 0., 5., 100);
     FillHist("dPhiemu_NljOSCut_TrkRecoIDIsoTrigPUW", electronColl.at(0).DeltaPhi(muonColl.at(0)), weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, -3.15, 3.15, 200);


     //Step5 : MET>40
     if(met<40) return;
     FillHist("Cutflow_NoW", 4., weight, 0., 10., 10);
     FillHist("Cutflow_TrkRecoW", 4., weight*trk_weight_mu*reco_weight_ele, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoW", 4., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoTrigW", 4., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoTrigPUW", 4., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, 0., 10., 10);



     float SFappliedWeight=weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight;
     FillHist("PTe", electronColl.at(0).Pt(), SFappliedWeight, 0., 200., 200);
     FillHist("Etae", electronColl.at(0).Eta(), SFappliedWeight, -5., 5., 100);
     FillHist("PTmu1", muonColl.at(0).Pt(), SFappliedWeight, 0., 200., 200);
     FillHist("Etamu1", muonColl.at(0).Eta(), SFappliedWeight, -5., 5., 100);
     FillHist("PTj1", jetColl.at(0).Pt(), SFappliedWeight, 0., 500., 500);
     FillHist("Etaj1", jetColl.at(0).Eta(), SFappliedWeight, -5., 5., 100);
     FillHist("PTj2", jetColl.at(1).Pt(), SFappliedWeight, 0., 500., 500);
     FillHist("Etaj2", jetColl.at(1).Eta(), SFappliedWeight, -5., 5., 100);
//     FillHist("PTj3", jetColl.at(2).Pt(), SFappliedWeight, 0., 200., 200);
     FillHist("MET", met, SFappliedWeight, 0., 500., 500);
     FillHist("Nvtx", Nvtx, SFappliedWeight, 0., 50., 50);
     FillHist("dRemu", electronColl.at(0).DeltaR(muonColl.at(0)), SFappliedWeight, 0., 5., 100);
     FillHist("dPhiemu", electronColl.at(0).DeltaPhi(muonColl.at(0)), SFappliedWeight, -3.15, 3.15, 200);

     FillHist("Nb_raw", nbjets, SFappliedWeight, 0., 10., 10);
     FillHist("Nb_1aSFed", nbjets, SFappliedWeight*btag_sf, 0., 10., 10);
     FillHist("Nb_2aSFed", bjetColl2a.size(), SFappliedWeight, 0., 10., 10);
     if(nbjets!=0){
       FillHist("PTb1_raw", bjetColl.at(0).Pt(), SFappliedWeight, 0., 200., 200);
       FillHist("PTb1_1aSFed", bjetColl.at(0).Pt(), SFappliedWeight*btag_sf, 0., 200., 200);

       FillHist("Etab1_raw", bjetColl.at(0).Eta(), SFappliedWeight, -5., 5., 100);
       FillHist("Etab1_1aSFed", bjetColl.at(0).Eta(), SFappliedWeight*btag_sf, -5., 5., 100);
     }
     if(bjetColl2a.size()!=0){
       FillHist("PTb1_2aSFed", bjetColl2a.at(0).Pt(), SFappliedWeight, 0., 200., 200);
       FillHist("Etab1_2aSFed", bjetColl2a.at(0).Eta(), SFappliedWeight, -5., 5., 100);
     }



     //Step6 : Nj>=3
     if(jetColl.size()<3) return;
     FillHist("Cutflow_NoW", 5., weight, 0., 10., 10);
     FillHist("Cutflow_TrkRecoW", 5., weight*trk_weight_mu*reco_weight_ele, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoW", 5., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoTrigW", 5., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoTrigPUW", 5., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, 0., 10., 10);


     FillHist("Nb_Nl3jMETCut_raw", nbjets, SFappliedWeight, 0., 10., 10);
     FillHist("Nb_Nl3jMETCut_1aSFed", nbjets, SFappliedWeight*btag_sf, 0., 10., 10);
     FillHist("Nb_Nl3jMETCut_2aSFed", bjetColl2a.size(), SFappliedWeight, 0., 10., 10);
     if(nbjets!=0){
       FillHist("PTb1_Nl3jMETCut_raw", bjetColl.at(0).Pt(), SFappliedWeight, 0., 200., 200);
       FillHist("PTb1_Nl3jMETCut_1aSFed", bjetColl.at(0).Pt(), SFappliedWeight*btag_sf, 0., 200., 200);

       FillHist("Etab1_Nl3jMETCut_raw", bjetColl.at(0).Eta(), SFappliedWeight, -5., 5., 100);
       FillHist("Etab1_Nl3jMETCut_1aSFed", bjetColl.at(0).Eta(), SFappliedWeight*btag_sf, -5., 5., 100);
     }
     if(bjetColl2a.size()!=0){
       FillHist("PTb1_Nl3jMETCut_2aSFed", bjetColl2a.at(0).Pt(), SFappliedWeight, 0., 200., 200);
       FillHist("Etab1_Nl3jMETCut_2aSFed", bjetColl2a.at(0).Eta(), SFappliedWeight, -5., 5., 100);
     }




     //Step7 : Nb>=1
     if(bjetColl.size()==0) return;
     FillHist("Cutflow_NoW", 6., weight, 0., 10., 10);
     FillHist("Cutflow_TrkRecoW", 6., weight*trk_weight_mu*reco_weight_ele, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoW", 6., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoTrigW", 6., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf, 0., 10., 10);
     FillHist("Cutflow_TrkRecoIDIsoTrigPUW", 6., weight*trk_weight_mu*reco_weight_ele*id_weight_ele*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, 0., 10., 10);



   }
   if(SingleMu_analysis){
    bool cycle1=false, cycle2=true;
    if(cycle1){

    if(muonLooseColl.size()==2){
     if(SumCharge(muonLooseColl)==0){
      if(muonLooseColl.at(0).Pt()<27){
     
        if(fabs((muonLooseColl.at(0)+muonLooseColl.at(1)).M()-91.2)>10) return;
   
        FillHist("PTmu1_Zwin", muonLooseColl.at(0).Pt(), weight, 0., 200., 200);
        FillHist("PTmu2_Zwin", muonLooseColl.at(1).Pt(), weight, 0., 200., 200);
   
        float PTmu2=muonLooseColl.at(1).Pt();
        float PTmu1=muonLooseColl.at(0).Pt();
   
        FillHist("dxy_mu", fabs(muonLooseColl.at(0).dXY()), weight, 0., 1., 1000);
        FillHist("dz_mu", fabs(muonLooseColl.at(0).dZ()), weight, 0., 1., 1000);
        FillHist("dxysig_mu", fabs(muonLooseColl.at(0).dXYSig()), weight, 0., 10., 1000);
   
        if(PTmu2<10){
          FillHist("dxy_mu_PT5_10", fabs(muonLooseColl.at(1).dXY()), weight, 0., 1., 1000);
          FillHist("dz_mu_PT5_10", fabs(muonLooseColl.at(1).dZ()), weight, 0., 1., 1000);
          FillHist("dxysig_mu_PT5_10", fabs(muonLooseColl.at(1).dXYSig()), weight, 0., 10., 1000);
        }
        if(PTmu2<20){//default cut muPT>10
          FillHist("dxy_mu_PT10_20", fabs(muonLooseColl.at(1).dXY()), weight, 0., 1., 1000);
          FillHist("dz_mu_PT10_20", fabs(muonLooseColl.at(1).dZ()), weight, 0., 1., 1000);
          FillHist("dxysig_mu_PT10_20", fabs(muonLooseColl.at(1).dXYSig()), weight, 0., 10., 1000);
        }
   
        if(PTmu1<30){//default cut muPT>20
          FillHist("dxy_mu_PT20_30", fabs(muonLooseColl.at(0).dXY()), weight, 0., 1., 1000);
          FillHist("dz_mu_PT20_30", fabs(muonLooseColl.at(0).dZ()), weight, 0., 1., 1000);
          FillHist("dxysig_mu_PT20_30", fabs(muonLooseColl.at(0).dXYSig()), weight, 0., 10., 1000);
        }
        else if(PTmu1<40){
          FillHist("dxy_mu_PT30_40", fabs(muonLooseColl.at(0).dXY()), weight, 0., 1., 1000);
          FillHist("dz_mu_PT30_40", fabs(muonLooseColl.at(0).dZ()), weight, 0., 1., 1000);
          FillHist("dxysig_mu_PT30_40", fabs(muonLooseColl.at(0).dXYSig()), weight, 0., 10., 1000);
        }
        else if(PTmu1<50){
          FillHist("dxy_mu_PT40_50", fabs(muonLooseColl.at(0).dXY()), weight, 0., 1., 1000);
          FillHist("dz_mu_PT40_50", fabs(muonLooseColl.at(0).dZ()), weight, 0., 1., 1000);
          FillHist("dxysig_mu_PT40_50", fabs(muonLooseColl.at(0).dXYSig()), weight, 0., 10., 1000);
        }
        else if(PTmu1<60){
          FillHist("dxy_mu_PT50_60", fabs(muonLooseColl.at(0).dXY()), weight, 0., 1., 1000);
          FillHist("dz_mu_PT50_60", fabs(muonLooseColl.at(0).dZ()), weight, 0., 1., 1000);
          FillHist("dxysig_mu_PT50_60", fabs(muonLooseColl.at(0).dXYSig()), weight, 0., 10., 1000);
        }
        else if(PTmu1<70){
          FillHist("dxy_mu_PT60_70", fabs(muonLooseColl.at(0).dXY()), weight, 0., 1., 1000);
          FillHist("dz_mu_PT60_70", fabs(muonLooseColl.at(0).dZ()), weight, 0., 1., 1000);
          FillHist("dxysig_mu_PT60_70", fabs(muonLooseColl.at(0).dXYSig()), weight, 0., 10., 1000);
        }
        else if(PTmu1<80){
          FillHist("dxy_mu_PT70_80", fabs(muonLooseColl.at(0).dXY()), weight, 0., 1., 1000);
          FillHist("dz_mu_PT70_80", fabs(muonLooseColl.at(0).dZ()), weight, 0., 1., 1000);
          FillHist("dxysig_mu_PT70_80", fabs(muonLooseColl.at(0).dXYSig()), weight, 0., 10., 1000);
        }
        else if(PTmu1<90){
          FillHist("dxy_mu_PT80_90", fabs(muonLooseColl.at(0).dXY()), weight, 0., 1., 1000);
          FillHist("dz_mu_PT80_90", fabs(muonLooseColl.at(0).dZ()), weight, 0., 1., 1000);
          FillHist("dxysig_mu_PT80_90", fabs(muonLooseColl.at(0).dXYSig()), weight, 0., 10., 1000);
        }
        else if(PTmu1<100){
          FillHist("dxy_mu_PT90_100", fabs(muonLooseColl.at(0).dXY()), weight, 0., 1., 1000);
          FillHist("dz_mu_PT90_100", fabs(muonLooseColl.at(0).dZ()), weight, 0., 1., 1000);
          FillHist("dxysig_mu_PT90_100", fabs(muonLooseColl.at(0).dXYSig()), weight, 0., 10., 1000);
        }
        else if(PTmu1<130){
          FillHist("dxy_mu_PT100_130", fabs(muonLooseColl.at(0).dXY()), weight, 0., 1., 1000);
          FillHist("dz_mu_PT100_130", fabs(muonLooseColl.at(0).dZ()), weight, 0., 1., 1000);
          FillHist("dxysig_mu_PT100_130", fabs(muonLooseColl.at(0).dXYSig()), weight, 0., 10., 1000);
        }
        else if(PTmu1<200){
          FillHist("dxy_mu_PT130_200", fabs(muonLooseColl.at(0).dXY()), weight, 0., 1., 1000);
          FillHist("dz_mu_PT130_200", fabs(muonLooseColl.at(0).dZ()), weight, 0., 1., 1000);
          FillHist("dxysig_mu_PT130_200", fabs(muonLooseColl.at(0).dXYSig()), weight, 0., 10., 1000);
        }

      }
     }
    }
    
    }//End of cycle1
    if(cycle2){

      if( !(PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v")) ) return;
      if( !(muonColl.size()==2) )     return;
      if( !(electronColl.size()==0) )  return;
      if( SumCharge(muonColl)!=0 ) return;
      if( !(muonColl.at(0).Pt()>27 && muonColl.at(1).Pt()>27) )  return;
      if( !(fabs((muonColl.at(0)+muonColl.at(1)).M()-91.2)<15) ) return;
     
      float weightSFapplied = weight*trigger_sf*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight;


      FillHist("PTmu1_TrkIDIsoTrigPUW", muonColl.at(0).Pt(), weightSFapplied, 0., 200., 200);
      FillHist("PTmu2_TrkIDIsoTrigPUW", muonColl.at(1).Pt(), weightSFapplied, 0., 200., 200);
      FillHist("Etamu1_TrkIDIsoTrigPUW", muonColl.at(0).Eta(), weightSFapplied, -5., 5., 100);
      FillHist("Etamu2_TrkIDIsoTrigPUW", muonColl.at(1).Eta(), weightSFapplied, -5., 5., 100);
      FillHist("Nj_TrkIDIsoTrigPUW", jetColl.size(), weightSFapplied, 0., 10., 10);
      if(njets>0){
        FillHist("PTj1_TrkIDIsoTrigPUW", jetColl.at(0).Pt(), weightSFapplied, 0., 200., 200);
        FillHist("Etaj1_TrkIDIsoTrigPUW", jetColl.at(0).Eta(), weightSFapplied, -5., 5., 100);
      }
      if(njets>1){
        FillHist("PTj2_TrkIDIsoTrigPUW", jetColl.at(1).Pt(), weightSFapplied, 0., 200., 200);
        FillHist("Etaj2_TrkIDIsoTrigPUW", jetColl.at(1).Eta(), weightSFapplied, -5., 5., 100);
      }
      FillHist("MET_TrkIDIsoTrigPUW", met, weightSFapplied, 0., 500., 500);
      FillHist("Nvtx_TrkIDIsoTrigPUW", Nvtx, weightSFapplied, 0., 50., 50);

      FillHist("Iso_TrkIDIsoTrigPUW", muonColl.at(0).RelIso04(), weightSFapplied, 0., 0.2, 1000);
      FillHist("dxy_TrkIDIsoTrigPUW", muonColl.at(0).dXY(), weightSFapplied, 0., 0.05, 1000);
      FillHist("dz_TrkIDIsoTrigPUW", muonColl.at(0).dZ(), weightSFapplied, 0., 0.1, 1000);
      FillHist("NhitTkLayer_TrkIDIsoTrigPUW", muonColl.at(0).ActiveLayer(), weightSFapplied, 0., 20., 20);//Nhit on trackLayer
      FillHist("Nhitpixel_TrkIDIsoTrigPUW", muonColl.at(0).validPixHits(), weightSFapplied, 0., 20., 20);//Nhit on pixel
      FillHist("NhitMuChamber_TrkIDIsoTrigPUW", muonColl.at(0).validHits(), weightSFapplied, 0., 50., 50);//Nhit on muon chamber
      FillHist("NmatchedStation_TrkIDIsoTrigPUW", muonColl.at(0).validStations(), weightSFapplied, 0., 20., 20);//N muon station having muon segment

      if(njets==0){
        FillHist("PTmu1_0j_TrkIDIsoTrigPUW", muonColl.at(0).Pt(), weightSFapplied, 0., 200., 200);
        FillHist("PTmu2_0j_TrkIDIsoTrigPUW", muonColl.at(1).Pt(), weightSFapplied, 0., 200., 200);
        FillHist("Etamu1_0j_TrkIDIsoTrigPUW", muonColl.at(0).Eta(), weightSFapplied, -5., 5., 100);
        FillHist("Etamu2_0j_TrkIDIsoTrigPUW", muonColl.at(1).Eta(), weightSFapplied, -5., 5., 100);
        FillHist("Nj_0j_TrkIDIsoTrigPUW", jetColl.size(), weightSFapplied, 0., 10., 10);
        FillHist("MET_0j_TrkIDIsoTrigPUW", met, weightSFapplied, 0., 500., 500);
        FillHist("Nvtx_0j_TrkIDIsoTrigPUW", Nvtx, weightSFapplied, 0., 50., 50);
      }
      if(njets==1){
        FillHist("PTmu1_1j_TrkIDIsoTrigPUW", muonColl.at(0).Pt(), weightSFapplied, 0., 200., 200);
        FillHist("PTmu2_1j_TrkIDIsoTrigPUW", muonColl.at(1).Pt(), weightSFapplied, 0., 200., 200);
        FillHist("Etamu1_1j_TrkIDIsoTrigPUW", muonColl.at(0).Eta(), weightSFapplied, -5., 5., 100);
        FillHist("Etamu2_1j_TrkIDIsoTrigPUW", muonColl.at(1).Eta(), weightSFapplied, -5., 5., 100);
        FillHist("Nj_1j_TrkIDIsoTrigPUW", jetColl.size(), weightSFapplied, 0., 10., 10);
        FillHist("PTj1_1j_TrkIDIsoTrigPUW", jetColl.at(0).Pt(), weightSFapplied, 0., 200., 200);
        FillHist("Etaj1_1j_TrkIDIsoTrigPUW", jetColl.at(0).Eta(), weightSFapplied, -5., 5., 100);
        FillHist("MET_1j_TrkIDIsoTrigPUW", met, weightSFapplied, 0., 500., 500);
        FillHist("Nvtx_1j_TrkIDIsoTrigPUW", Nvtx, weightSFapplied, 0., 50., 50);
      }
      if(njets==2){
        FillHist("PTmu1_2j_TrkIDIsoTrigPUW", muonColl.at(0).Pt(), weightSFapplied, 0., 200., 200);
        FillHist("PTmu2_2j_TrkIDIsoTrigPUW", muonColl.at(1).Pt(), weightSFapplied, 0., 200., 200);
        FillHist("Etamu1_2j_TrkIDIsoTrigPUW", muonColl.at(0).Eta(), weightSFapplied, -5., 5., 100);
        FillHist("Etamu2_2j_TrkIDIsoTrigPUW", muonColl.at(1).Eta(), weightSFapplied, -5., 5., 100);
        FillHist("Nj_2j_TrkIDIsoTrigPUW", jetColl.size(), weightSFapplied, 0., 10., 10);
        FillHist("PTj1_2j_TrkIDIsoTrigPUW", jetColl.at(0).Pt(), weightSFapplied, 0., 200., 200);
        FillHist("Etaj1_2j_TrkIDIsoTrigPUW", jetColl.at(0).Eta(), weightSFapplied, -5., 5., 100);
        FillHist("PTj2_2j_TrkIDIsoTrigPUW", jetColl.at(1).Pt(), weightSFapplied, 0., 200., 200);
        FillHist("Etaj2_2j_TrkIDIsoTrigPUW", jetColl.at(1).Eta(), weightSFapplied, -5., 5., 100);
        FillHist("MET_2j_TrkIDIsoTrigPUW", met, weightSFapplied, 0., 500., 500);
        FillHist("Nvtx_2j_TrkIDIsoTrigPUW", Nvtx, weightSFapplied, 0., 50., 50);
      }
    
   
    }//End of cycle 2
   }
   if(DoubleMu_analysis){
     if(PassTrigger(analysis_trigger)){
       FillHist("Counter_Nobj_NoW", 0., muonLooseColl.size()*weight, 0., 10., 10);
       FillHist("Counter_Nobj_NoW", 1., muonColl.size()*weight, 0., 10., 10);
       FillHist("Counter_Nobj_NoW", 2., jetColl.size()*weight, 0., 10., 10);
       FillHist("Counter_Nobj_NoW", 3., bjetColl.size()*weight, 0., 10., 10);

       if(muonColl.size()!=2) return;
       if(SumCharge(muonColl)!=0) return;
       if(muonColl.at(0).Pt()<20) return;
       if(muonColl.at(1).Pt()<20) return;
       if(fabs((muonColl.at(0)+muonColl.at(1)).M()-91.2)>15) return;

       FillHist("PTmu1_NoW", muonColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTmu2_NoW", muonColl.at(1).Pt(), weight, 0., 200., 200);
       FillHist("Etamu1_NoW", muonColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Etamu2_NoW", muonColl.at(1).Eta(), weight, -5., 5., 100);
       FillHist("Nj_NoW", jetColl.size(), weight, 0., 10., 10);
       FillHist("NjL_NoW", jetLooseColl.size(), weight, 0., 10., 10);
       FillHist("MET_NoW", met, weight, 0., 500., 500);
       FillHist("Nvtx_NoW", Nvtx, weight, 0., 50., 50);
       
       FillHist("Iso_NoW", muonColl.at(0).RelIso04(), weight, 0., 0.2, 1000);
       FillHist("dxy_NoW", muonColl.at(0).dXY(), weight, 0., 0.4, 1000);
       FillHist("dz_NoW", muonColl.at(0).dZ(), weight, 0., 0.7, 1000);
       FillHist("NhitTkLayer_NoW", muonColl.at(0).ActiveLayer(), weight, 0., 20., 20);//Nhit on trackLayer
       FillHist("Nhitpixel_NoW", muonColl.at(0).validPixHits(), weight, 0., 20., 20);//Nhit on pixel
       FillHist("NhitMuChamber_NoW", muonColl.at(0).validHits(), weight, 0., 50., 50);//Nhit on muon chamber
       FillHist("NmatchedStation_NoW", muonColl.at(0).validStations(), weight, 0., 20., 20);//N muon station having muon segment
       
       FillHist("PTmu1_TkIDIsoPUW", muonColl.at(0).Pt(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);
       FillHist("PTmu2_TkIDIsoPUW", muonColl.at(1).Pt(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);
       FillHist("Etamu1_TkIDIsoPUW", muonColl.at(0).Eta(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, -5., 5., 100);
       FillHist("Etamu2_TkIDIsoPUW", muonColl.at(1).Eta(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, -5., 5., 100);
       FillHist("Nj_TkIDIsoPUW", jetColl.size(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 10., 10);
       FillHist("NjL_TkIDIsoPUW", jetLooseColl.size(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 10., 10);
       FillHist("MET_TkIDIsoPUW", met, weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 500., 500);
       FillHist("Nvtx_TkIDIsoPUW", Nvtx, weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 50., 50);
       
       FillHist("Iso_TkIDIsoPUW", muonColl.at(0).RelIso04(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 0.2, 1000);
       FillHist("dxy_TkIDIsoPUW", muonColl.at(0).dXY(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 0.4, 1000);
       FillHist("dz_TkIDIsoPUW", muonColl.at(0).dZ(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 0.7, 1000);
       FillHist("NhitTkLayer_TkIDIsoPUW", muonColl.at(0).ActiveLayer(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 20., 20);//Nhit on trackLayer
       FillHist("Nhitpixel_TkIDIsoPUW", muonColl.at(0).validPixHits(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 20., 20);//Nhit on pixel
       FillHist("NhitMuChamber_TkIDIsoPUW", muonColl.at(0).validHits(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 50., 50);//Nhit on muon chamber
       FillHist("NmatchedStation_TkIDIsoPUW", muonColl.at(0).validStations(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 20., 20);//N muon station having muon segment
     
       if(njets==0){
         FillHist("PTmu1_0j", muonColl.at(0).Pt(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);
         FillHist("PTmu2_0j", muonColl.at(1).Pt(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);
         FillHist("Etamu1_0j", muonColl.at(0).Eta(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, -5., 5., 100);
         FillHist("Etamu2_0j", muonColl.at(1).Eta(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, -5., 5., 100);
         FillHist("MET_0j", met, weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 500., 500);
       }
       if(njets==1){
         FillHist("PTmu1_1j", muonColl.at(0).Pt(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);
         FillHist("PTmu2_1j", muonColl.at(1).Pt(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);
         FillHist("Etamu1_1j", muonColl.at(0).Eta(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, -5., 5., 100);
         FillHist("Etamu2_1j", muonColl.at(1).Eta(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, -5., 5., 100);
         FillHist("MET_1j", met, weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 500., 500);
         FillHist("PTj1_1j", jetColl.at(0).Pt(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);
         FillHist("Etaj1_1j", jetColl.at(0).Eta(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, -5., 5., 100);
       }
       if(njets==2){
         FillHist("PTmu1_2j", muonColl.at(0).Pt(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);
         FillHist("PTmu2_2j", muonColl.at(1).Pt(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);
         FillHist("Etamu1_2j", muonColl.at(0).Eta(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, -5., 5., 100);
         FillHist("Etamu2_2j", muonColl.at(1).Eta(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, -5., 5., 100);
         FillHist("MET_2j", met, weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 500., 500);
         FillHist("PTj1_2j", jetColl.at(0).Pt(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);
         FillHist("PTj2_2j", jetColl.at(1).Pt(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);
         FillHist("Etaj1_2j", jetColl.at(0).Eta(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, -5., 5., 100);
         FillHist("Etaj2_2j", jetColl.at(1).Eta(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, -5., 5., 100);
       }
     }
   }
   if(DoubleEle_analysis){
     bool cycle1=false, cycle2=true;
     if(cycle1){
     ////////////////////////////////////////////////////////////////////////////

     if(electronLooseColl.size()!=2) return;
     FillCutFlow("NlCut(2l)", weight);
     if(electronLooseColl.at(0).Charge()==electronLooseColl.at(1).Charge()) return;
     FillCutFlow("OSlep", weight);
     if(electronLooseColl.at(0).Pt()<30) return;
     
     if(fabs((electronLooseColl.at(0)+electronLooseColl.at(1)).M()-91.2)>10) return;

     FillHist("PTe1_Zwin", electronLooseColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("PTe2_Zwin", electronLooseColl.at(1).Pt(), weight, 0., 200., 200);

     if(fabs(electronLooseColl.at(0).Eta())<1.479){
       FillHist("PTeB1_Zwin", electronLooseColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTeB2_Zwin", electronLooseColl.at(1).Pt(), weight, 0., 200., 200);

       FillHist("dxy_eleB", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
       FillHist("dz_eleB", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
     }
     else{
       FillHist("PTeE1_Zwin", electronLooseColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("PTeE2_Zwin", electronLooseColl.at(1).Pt(), weight, 0., 200., 200);

       FillHist("dxy_eleE", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
       FillHist("dz_eleE", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
     }


     float PTele2=electronLooseColl.at(1).Pt();
     float PTele1=electronLooseColl.at(0).Pt();
     if(fabs(electronLooseColl.at(1).Eta())<1.479){
       if(PTele2<20){//default cut elePT>15
         FillHist("dxy_eleB_PT15_20", fabs(electronLooseColl.at(1).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleB_PT15_20", fabs(electronLooseColl.at(1).dz()), weight, 0., 1., 1000);
       }
       else if(PTele2<30){
         FillHist("dxy_eleB_PT20_30", fabs(electronLooseColl.at(1).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleB_PT20_30", fabs(electronLooseColl.at(1).dz()), weight, 0., 1., 1000);
       }
     }
     else{
       if(PTele2<20){//default cut elePT>15
         FillHist("dxy_eleE_PT15_20", fabs(electronLooseColl.at(1).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleE_PT15_20", fabs(electronLooseColl.at(1).dz()), weight, 0., 1., 1000);
       }
       else if(PTele2<30){
         FillHist("dxy_eleE_PT20_30", fabs(electronLooseColl.at(1).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleE_PT20_30", fabs(electronLooseColl.at(1).dz()), weight, 0., 1., 1000);
       }
     }


     if(fabs(electronLooseColl.at(0).Eta())<1.479){
       if(PTele1<40){//defaultcut elePT>30
         FillHist("dxy_eleB_PT30_40", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleB_PT30_40", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
       else if(PTele1<50){
         FillHist("dxy_eleB_PT40_50", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleB_PT40_50", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
       else if(PTele1<60){
         FillHist("dxy_eleB_PT50_60", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleB_PT50_60", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
       else if(PTele1<70){
         FillHist("dxy_eleB_PT60_70", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleB_PT60_70", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
       else if(PTele1<80){
         FillHist("dxy_eleB_PT70_80", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleB_PT70_80", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
       else if(PTele1<90){
         FillHist("dxy_eleB_PT80_90", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleB_PT80_90", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
       else if(PTele1<100){
         FillHist("dxy_eleB_PT90_100", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleB_PT90_100", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
       else if(PTele1<130){
         FillHist("dxy_eleB_PT100_130", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleB_PT100_130", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
       else if(PTele1<200){
         FillHist("dxy_eleB_PT130_200", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleB_PT130_200", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
     }
     else{
       if(PTele1<40){//default cut elePT>30
         FillHist("dxy_eleE_PT30_40", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleE_PT30_40", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
       else if(PTele1<50){
         FillHist("dxy_eleE_PT40_50", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleE_PT40_50", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
       else if(PTele1<60){
         FillHist("dxy_eleE_PT50_60", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleE_PT50_60", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
       else if(PTele1<70){
         FillHist("dxy_eleE_PT60_70", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleE_PT60_70", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
       else if(PTele1<80){
         FillHist("dxy_eleE_PT70_80", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleE_PT70_80", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
       else if(PTele1<90){
         FillHist("dxy_eleE_PT80_90", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleE_PT80_90", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
       else if(PTele1<100){
         FillHist("dxy_eleE_PT90_100", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleE_PT90_100", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
       else if(PTele1<130){
         FillHist("dxy_eleE_PT100_130", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleE_PT100_130", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }
       else if(PTele1<200){
         FillHist("dxy_eleE_PT130_200", fabs(electronLooseColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleE_PT130_200", fabs(electronLooseColl.at(0).dz()), weight, 0., 1., 1000);
       }

     }//EndCap End
     /////////////////////////////////
     }//Cycle1 End
     if(cycle2){

     if(isData){
       if(!PassTrigger(analysis_trigger)) return;
     }
     if(electronColl.size()!=2) return;
     if(electronColl.at(0).Charge()==electronColl.at(1).Charge()) return;
     if(electronColl.at(0).Pt()<25) return;
     if(electronColl.at(1).Pt()<15) return;
     if(fabs((electronColl.at(0)+electronColl.at(1)).M()-91.2)<15){ 
     ///////////////////////////////////////////////////////////////////////////
     
       if(!isData){
         if(PassTrigger(analysis_trigger)){
           FillHist("PTele1_NoW_TrigCut", electronColl.at(0).Pt(), weight, 0., 200., 200);
           FillHist("Etaele1_NoW_TrigCut", electronColl.at(0).Eta(), weight, -5., 5., 200);
           FillHist("PTele2_NoW_TrigCut", electronColl.at(1).Pt(), weight, 0., 200., 200);
           FillHist("Etaele2_NoW_TrigCut", electronColl.at(1).Eta(), weight, -5., 5., 200);
           FillHist("Nj_NoW_TrigCut", jetColl.size(), weight, 0., 10., 10);
           FillHist("NjL_NoW_TrigCut", jetLooseColl.size(), weight, 0., 10., 10);
           FillHist("Nvtx_NoW_TrigCut", Nvtx, weight, 0., 50., 50);
      
           FillHist("PTele1_IDSF_TrigCut", electronColl.at(0).Pt(), weight*id_weight_ele, 0., 200., 200);
           FillHist("Etaele1_IDSF_TrigCut", electronColl.at(0).Eta(), weight*id_weight_ele, -5., 5., 200);
           FillHist("PTele2_IDSF_TrigCut", electronColl.at(1).Pt(), weight*id_weight_ele, 0., 200., 200);
           FillHist("Etaele2_IDSF_TrigCut", electronColl.at(1).Eta(), weight*id_weight_ele, -5., 5., 200);
           FillHist("Nj_IDSF_TrigCut", jetColl.size(), weight*id_weight_ele, 0., 10., 10);
           FillHist("NjL_IDSF_TrigCut", jetLooseColl.size(), weight*id_weight_ele, 0., 10., 10);
           FillHist("Nvtx_IDSF_TrigCut", Nvtx, weight*id_weight_ele, 0., 50., 50);
      
           FillHist("PTele1_IDPUSF_TrigCut", electronColl.at(0).Pt(), weight*id_weight_ele*pileup_reweight, 0., 200., 200);
           FillHist("Etaele1_IDPUSF_TrigCut", electronColl.at(0).Eta(), weight*id_weight_ele*pileup_reweight, -5., 5., 200);
           FillHist("PTele2_IDPUSF_TrigCut", electronColl.at(1).Pt(), weight*id_weight_ele*pileup_reweight, 0., 200., 200);
           FillHist("Etaele2_IDPUSF_TrigCut", electronColl.at(1).Eta(), weight*id_weight_ele*pileup_reweight, -5., 5., 200);
           FillHist("Nj_IDPUSF_TrigCut", jetColl.size(), weight*id_weight_ele*pileup_reweight, 0., 10., 10);
           FillHist("NjL_IDPUSF_TrigCut", jetLooseColl.size(), weight*id_weight_ele*pileup_reweight, 0., 10., 10);
           FillHist("Nvtx_IDPUSF_TrigCut", Nvtx, weight*id_weight_ele*pileup_reweight, 0., 50., 50);
         }
           FillHist("PTele1_NoW_TrigW", electronColl.at(0).Pt(), weight*trigger_eff, 0., 200., 200);
           FillHist("Etaele1_NoW_TrigW", electronColl.at(0).Eta(), weight*trigger_eff, -5., 5., 200);
           FillHist("PTele2_NoW_TrigW", electronColl.at(1).Pt(), weight*trigger_eff, 0., 200., 200);
           FillHist("Etaele2_NoW_TrigW", electronColl.at(1).Eta(), weight*trigger_eff, -5., 5., 200);
           FillHist("Nj_NoW_TrigW", jetColl.size(), weight*trigger_eff, 0., 10., 10);
           FillHist("NjL_NoW_TrigW", jetLooseColl.size(), weight*trigger_eff, 0., 10., 10);
           FillHist("Nvtx_NoW_TrigW", Nvtx, weight*trigger_eff, 0., 50., 50);
      
           FillHist("PTele1_IDSF_TrigW", electronColl.at(0).Pt(), weight*trigger_eff*id_weight_ele, 0., 200., 200);
           FillHist("Etaele1_IDSF_TrigW", electronColl.at(0).Eta(), weight*trigger_eff*id_weight_ele, -5., 5., 200);
           FillHist("PTele2_IDSF_TrigW", electronColl.at(1).Pt(), weight*trigger_eff*id_weight_ele, 0., 200., 200);
           FillHist("Etaele2_IDSF_TrigW", electronColl.at(1).Eta(), weight*trigger_eff*id_weight_ele, -5., 5., 200);
           FillHist("Nj_IDSF_TrigW", jetColl.size(), weight*trigger_eff*id_weight_ele, 0., 10., 10);
           FillHist("NjL_IDSF_TrigW", jetLooseColl.size(), weight*trigger_eff*id_weight_ele, 0., 10., 10);
           FillHist("Nvtx_IDSF_TrigW", Nvtx, weight*trigger_eff*id_weight_ele, 0., 50., 50);
      
           FillHist("PTele1_IDPUSF_TrigW", electronColl.at(0).Pt(), weight*trigger_eff*id_weight_ele*pileup_reweight, 0., 200., 200);
           FillHist("Etaele1_IDPUSF_TrigW", electronColl.at(0).Eta(), weight*trigger_eff*id_weight_ele*pileup_reweight, -5., 5., 200);
           FillHist("PTele2_IDPUSF_TrigW", electronColl.at(1).Pt(), weight*trigger_eff*id_weight_ele*pileup_reweight, 0., 200., 200);
           FillHist("Etaele2_IDPUSF_TrigW", electronColl.at(1).Eta(), weight*trigger_eff*id_weight_ele*pileup_reweight, -5., 5., 200);
           FillHist("Nj_IDPUSF_TrigW", jetColl.size(), weight*trigger_eff*id_weight_ele*pileup_reweight, 0., 10., 10);
           FillHist("NjL_IDPUSF_TrigW", jetLooseColl.size(), weight*trigger_eff*id_weight_ele*pileup_reweight, 0., 10., 10);
           FillHist("Nvtx_IDPUSF_TrigW", Nvtx, weight*trigger_eff*id_weight_ele*pileup_reweight, 0., 50., 50);
       }
       else{
           FillHist("PTele1_NoW_TrigCut", electronColl.at(0).Pt(), weight, 0., 200., 200);
           FillHist("Etaele1_NoW_TrigCut", electronColl.at(0).Eta(), weight, -5., 5., 200);
           FillHist("PTele2_NoW_TrigCut", electronColl.at(1).Pt(), weight, 0., 200., 200);
           FillHist("Etaele2_NoW_TrigCut", electronColl.at(1).Eta(), weight, -5., 5., 200);
           FillHist("Nj_NoW_TrigCut", jetColl.size(), weight, 0., 10., 10);
           FillHist("NjL_NoW_TrigCut", jetLooseColl.size(), weight, 0., 10., 10);
           FillHist("Nvtx_NoW_TrigCut", Nvtx, weight, 0., 50., 50);
      
           FillHist("PTele1_IDSF_TrigCut", electronColl.at(0).Pt(), weight*id_weight_ele, 0., 200., 200);
           FillHist("Etaele1_IDSF_TrigCut", electronColl.at(0).Eta(), weight*id_weight_ele, -5., 5., 200);
           FillHist("PTele2_IDSF_TrigCut", electronColl.at(1).Pt(), weight*id_weight_ele, 0., 200., 200);
           FillHist("Etaele2_IDSF_TrigCut", electronColl.at(1).Eta(), weight*id_weight_ele, -5., 5., 200);
           FillHist("Nj_IDSF_TrigCut", jetColl.size(), weight*id_weight_ele, 0., 10., 10);
           FillHist("NjL_IDSF_TrigCut", jetLooseColl.size(), weight*id_weight_ele, 0., 10., 10);
           FillHist("Nvtx_IDSF_TrigCut", Nvtx, weight*id_weight_ele, 0., 50., 50);
      
           FillHist("PTele1_IDPUSF_TrigCut", electronColl.at(0).Pt(), weight*id_weight_ele*pileup_reweight, 0., 200., 200);
           FillHist("Etaele1_IDPUSF_TrigCut", electronColl.at(0).Eta(), weight*id_weight_ele*pileup_reweight, -5., 5., 200);
           FillHist("PTele2_IDPUSF_TrigCut", electronColl.at(1).Pt(), weight*id_weight_ele*pileup_reweight, 0., 200., 200);
           FillHist("Etaele2_IDPUSF_TrigCut", electronColl.at(1).Eta(), weight*id_weight_ele*pileup_reweight, -5., 5., 200);
           FillHist("Nj_IDPUSF_TrigCut", jetColl.size(), weight*id_weight_ele*pileup_reweight, 0., 10., 10);
           FillHist("NjL_IDPUSF_TrigCut", jetLooseColl.size(), weight*id_weight_ele*pileup_reweight, 0., 10., 10);
           FillHist("Nvtx_IDPUSF_TrigCut", Nvtx, weight*id_weight_ele*pileup_reweight, 0., 50., 50);
  
           FillHist("PTele1_NoW_TrigW", electronColl.at(0).Pt(), weight*trigger_eff, 0., 200., 200);
           FillHist("Etaele1_NoW_TrigW", electronColl.at(0).Eta(), weight*trigger_eff, -5., 5., 200);
           FillHist("PTele2_NoW_TrigW", electronColl.at(1).Pt(), weight*trigger_eff, 0., 200., 200);
           FillHist("Etaele2_NoW_TrigW", electronColl.at(1).Eta(), weight*trigger_eff, -5., 5., 200);
           FillHist("Nj_NoW_TrigW", jetColl.size(), weight*trigger_eff, 0., 10., 10);
           FillHist("NjL_NoW_TrigW", jetLooseColl.size(), weight*trigger_eff, 0., 10., 10);
           FillHist("Nvtx_NoW_TrigW", Nvtx, weight*trigger_eff, 0., 50., 50);
      
           FillHist("PTele1_IDSF_TrigW", electronColl.at(0).Pt(), weight*trigger_eff*id_weight_ele, 0., 200., 200);
           FillHist("Etaele1_IDSF_TrigW", electronColl.at(0).Eta(), weight*trigger_eff*id_weight_ele, -5., 5., 200);
           FillHist("PTele2_IDSF_TrigW", electronColl.at(1).Pt(), weight*trigger_eff*id_weight_ele, 0., 200., 200);
           FillHist("Etaele2_IDSF_TrigW", electronColl.at(1).Eta(), weight*trigger_eff*id_weight_ele, -5., 5., 200);
           FillHist("Nj_IDSF_TrigW", jetColl.size(), weight*trigger_eff*id_weight_ele, 0., 10., 10);
           FillHist("NjL_IDSF_TrigW", jetLooseColl.size(), weight*trigger_eff*id_weight_ele, 0., 10., 10);
           FillHist("Nvtx_IDSF_TrigW", Nvtx, weight*trigger_eff*id_weight_ele, 0., 50., 50);
      
           FillHist("PTele1_IDPUSF_TrigW", electronColl.at(0).Pt(), weight*trigger_eff*id_weight_ele*pileup_reweight, 0., 200., 200);
           FillHist("Etaele1_IDPUSF_TrigW", electronColl.at(0).Eta(), weight*trigger_eff*id_weight_ele*pileup_reweight, -5., 5., 200);
           FillHist("PTele2_IDPUSF_TrigW", electronColl.at(1).Pt(), weight*trigger_eff*id_weight_ele*pileup_reweight, 0., 200., 200);
           FillHist("Etaele2_IDPUSF_TrigW", electronColl.at(1).Eta(), weight*trigger_eff*id_weight_ele*pileup_reweight, -5., 5., 200);
           FillHist("Nj_IDPUSF_TrigW", jetColl.size(), weight*trigger_eff*id_weight_ele*pileup_reweight, 0., 10., 10);
           FillHist("NjL_IDPUSF_TrigW", jetLooseColl.size(), weight*trigger_eff*id_weight_ele*pileup_reweight, 0., 10., 10);
           FillHist("Nvtx_IDPUSF_TrigW", Nvtx, weight*trigger_eff*id_weight_ele*pileup_reweight, 0., 50., 50);
       }
         
  /*
       if(fabs(electronColl.at(0).Eta())<1.479){
         FillHist("PTeleB_NoW", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("dxy_eleB_NoW", fabs(electronColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleB_NoW", fabs(electronColl.at(0).dz()), weight, 0., 1., 1000);
         FillHist("Iso_eleB_NoW", electronColl.at(0).PFRelIso(0.3), weight, 0., 0.2, 1000);
  
         FillHist("PTeleB_IDSF", electronColl.at(0).Pt(), weight*id_weight_ele, 0., 200., 200);
         FillHist("dxy_eleB_IDSF", fabs(electronColl.at(0).dxy()), weight*id_weight_ele, 0., 1., 1000);
         FillHist("dz_eleB_IDSF", fabs(electronColl.at(0).dz()), weight*id_weight_ele, 0., 1., 1000);
         FillHist("Iso_eleB_IDSF", electronColl.at(0).PFRelIso(0.3), weight*id_weight_ele, 0., 0.2, 1000);
  
         FillHist("PTeleB_IDPUSF", electronColl.at(0).Pt(), weight*id_weight_ele*pileup_reweight, 0., 200., 200);
         FillHist("dxy_eleB_IDPUSF", fabs(electronColl.at(0).dxy()), weight*id_weight_ele*pileup_reweight, 0., 1., 1000);
         FillHist("dz_eleB_IDPUSF", fabs(electronColl.at(0).dz()), weight*id_weight_ele*pileup_reweight, 0., 1., 1000);
         FillHist("Iso_eleB_IDPUSF", electronColl.at(0).PFRelIso(0.3), weight*id_weight_ele*pileup_reweight, 0., 0.2, 1000);
  
       }
       else{
         FillHist("PTeleE_NoW", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("dxy_eleE_NoW", fabs(electronColl.at(0).dxy()), weight, 0., 1., 1000);
         FillHist("dz_eleE_NoW", fabs(electronColl.at(0).dz()), weight, 0., 1., 1000);
         FillHist("Iso_eleE_NoW", electronColl.at(0).PFRelIso(0.3), weight, 0., 0.2, 1000);
  
         FillHist("PTeleE_IDSF", electronColl.at(0).Pt(), weight*id_weight_ele, 0., 200., 200);
         FillHist("dxy_eleE_IDSF", fabs(electronColl.at(0).dxy()), weight*id_weight_ele, 0., 1., 1000);
         FillHist("dz_eleE_IDSF", fabs(electronColl.at(0).dz()), weight*id_weight_ele, 0., 1., 1000);
         FillHist("Iso_eleE_IDSF", electronColl.at(0).PFRelIso(0.3), weight*id_weight_ele, 0., 0.2, 1000);
  
         FillHist("PTeleE_IDPUSF", electronColl.at(0).Pt(), weight*id_weight_ele*pileup_reweight, 0., 200., 200);
         FillHist("dxy_eleE_IDPUSF", fabs(electronColl.at(0).dxy()), weight*id_weight_ele*pileup_reweight, 0., 1., 1000);
         FillHist("dz_eleE_IDPUSF", fabs(electronColl.at(0).dz()), weight*id_weight_ele*pileup_reweight, 0., 1., 1000);
         FillHist("Iso_eleE_IDPUSF", electronColl.at(0).PFRelIso(0.3), weight*id_weight_ele*pileup_reweight, 0., 0.2, 1000);
  
       }       
  */
       }//fabs(Mee-MZ)<15 End
       if((electronColl.at(0)+electronColl.at(1)).M()>50){
         if(PassTrigger(analysis_trigger)){
           FillHist("DY50Count_NoW_TrigCut", 0., weight, 0., 1., 1);
           FillHist("DY50Count_IDW_TrigCut", 0., weight*id_weight_ele, 0., 1., 1);
           FillHist("DY50Count_IDPUW_TrigCut", 0., weight*id_weight_ele*pileup_reweight, 0., 1., 1);

           FillHist("PTe1_NoW_TrigCut_DY50", electronColl.at(0).Pt(), weight, 0., 200., 200);
           FillHist("PTe2_NoW_TrigCut_DY50", electronColl.at(1).Pt(), weight, 0., 200., 200);
           FillHist("Etae1_NoW_TrigCut_DY50", electronColl.at(0).Eta(), weight, -5., 5., 200);
           FillHist("Etae2_NoW_TrigCut_DY50", electronColl.at(1).Eta(), weight, -5., 5., 200);
           FillHist("Mee_NoW_TrigCut_DY50", (electronColl.at(0)+electronColl.at(1)).M(), weight, 0., 200., 200);

           FillHist("PTe1_IDW_TrigCut_DY50", electronColl.at(0).Pt(), weight*id_weight_ele, 0., 200., 200);
           FillHist("PTe2_IDW_TrigCut_DY50", electronColl.at(1).Pt(), weight*id_weight_ele, 0., 200., 200);
           FillHist("Etae1_IDW_TrigCut_DY50", electronColl.at(0).Eta(), weight*id_weight_ele, -5., 5., 200);
           FillHist("Etae2_IDW_TrigCut_DY50", electronColl.at(1).Eta(), weight*id_weight_ele, -5., 5., 200);
           FillHist("Mee_IDW_TrigCut_DY50", (electronColl.at(0)+electronColl.at(1)).M(), weight*id_weight_ele, 0., 200., 200);

           FillHist("PTe1_IDPUW_TrigCut_DY50", electronColl.at(0).Pt(), weight*id_weight_ele*pileup_reweight, 0., 200., 200);
           FillHist("PTe2_IDPUW_TrigCut_DY50", electronColl.at(1).Pt(), weight*id_weight_ele*pileup_reweight, 0., 200., 200);
           FillHist("Etae1_IDPUW_TrigCut_DY50", electronColl.at(0).Eta(), weight*id_weight_ele*pileup_reweight, -5., 5., 200);
           FillHist("Etae2_IDPUW_TrigCut_DY50", electronColl.at(1).Eta(), weight*id_weight_ele*pileup_reweight, -5., 5., 200);
           FillHist("Mee_IDPUW_TrigCut_DY50", (electronColl.at(0)+electronColl.at(1)).M(), weight*id_weight_ele*pileup_reweight, 0., 200., 200);

         }

       }
     }//Cycle2 End


   }//DoubleEle End




/////////////////////////////////////////////////////////////////////////////////// 


return;
}// End of execute event loop
  


void Feb2017_3l4j_BTagClosure::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Feb2017_3l4j_BTagClosure::BeginCycle() throw( LQError ){
  
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

Feb2017_3l4j_BTagClosure::~Feb2017_3l4j_BTagClosure() {
  
  Message("In Feb2017_3l4j_BTagClosure Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Feb2017_3l4j_BTagClosure::FillCutFlow(TString cut, float weight){
  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 12, 0., 12.);
    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"NlCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"OS(2mu)");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"NjCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(8,"NbCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(9,"NljCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(10,"lPtCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(11,"MjjCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(12,"MmumuCut");
    
  }
}



void Feb2017_3l4j_BTagClosure::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Feb2017_3l4j_BTagClosure::MakeHistograms(){
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
/*
  //After all preselections
  AnalyzerCore::MakeHistograms("Basic_MET_wPreSel", 100, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Nj_wPreSel", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Nb_wPreSel", 5, 0., 5.);
  AnalyzerCore::MakeHistograms("Basic_Jet1_Et", 50, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Jet2_Et", 50, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Jet3_Et", 50, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Jet4_Et", 50, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Jet1_Eta", 20, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Jet2_Eta", 20, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Jet3_Eta", 20, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Jet4_Eta", 20, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_bJet1_Et", 20, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_bJet2_Et", 20, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_bJet1_Eta", 20, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_bJet2_Eta", 20, -5., 5.);

  AnalyzerCore::MakeHistograms("Basic_Pte_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu1_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu2_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Etae_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etamu1_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etamu2_wPreSel", 100, -5., 5.);

  AnalyzerCore::MakeHistograms("Basic_Ptmu3_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Etamu3_wPreSel", 100, -5., 5.);

  AnalyzerCore::MakeHistograms("Basic_Ptj1_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptj2_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptj3_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptj4_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptb1_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptb2_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptlj1_wPreSel", 200, 0., 200.);

  AnalyzerCore::MakeHistograms("Basic_Etaj1_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etaj2_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etaj3_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etaj4_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etab1_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etab2_wPreSel", 100, -5., 5.);
*/

  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Feb2017_3l4j_BTagClosureCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Feb2017_3l4j_BTagClosure::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
