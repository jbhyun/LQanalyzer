// $Id: Jan2017_3l4j_DYCheck.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQJan2017_3l4j_DYCheck Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Jan2017_3l4j_DYCheck.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Jan2017_3l4j_DYCheck);

 Jan2017_3l4j_DYCheck::Jan2017_3l4j_DYCheck() : AnalyzerCore(), out_muons(0) {

   SetLogName("Jan2017_3l4j_DYCheck");
   Message("In Jan2017_3l4j_DYCheck constructor", INFO);
   InitialiseAnalysis();
 }


 void Jan2017_3l4j_DYCheck::InitialiseAnalysis() throw( LQError ) {
   
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

void Jan2017_3l4j_DYCheck::ExecuteEvents()throw( LQError ){

////Basic Infos///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Apply the gen weight 
   float weight_nopu=1; 
   if(!isData) weight*=MCweight;
   if(!isData) weight_nopu*=MCweight;
   FillHist("GenWeight", MCweight, 1, -10, 10, 1000);

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
   bool DoubleEle_analysis=false, SingleMu_analysis=true, DoubleMu_analysis=false, EMu_analysis=false;


   //Trigers
   std::vector<TString> triggerslist, triggerlist1, triggerlist2;//  ListTriggersAvailable();

   TString analysis_trigger;
   if(EMu_analysis) {analysis_trigger="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";}
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
     eventbase->GetMuonSel()->SetPt(20.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.15);//POG WP T
   std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->Selection(muonColl,true);//(muonColl, bool RochCorr, bool debug)
 //std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->SelectMuons(muonColl, BaseSelection::MUON_POG_TIGHT, 10., 2.4);
   std::vector<snu::KMuon> muonPromptColl; muonPromptColl=GetTruePrompt(muonColl, false);

     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0994, 0.107);//2016 80X tuned WP
   //eventbase->GetElectronSel()->SetdxyBEMax(0.005, 0.01);  eventbase->GetElectronSel()->SetdzBEMax(0.005, 0.01);//Not in ID, but additional safe WP
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(15.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.2);//Not in ID, but additional safe WP
   std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->Selection(electronColl);
 //std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->SelectElectrons(electronColl, "ELECTRON_POG_TIGHT", 20., 2.4);

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(20.);                     eventbase->GetJetSel()->SetEta(2.4);
   //eventbase->GetJetSel()->SetUseJetPileUp(true);
     bool LeptonVeto=true;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonColl, electronColl);
 //std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->SelectJets(jetColl, muonColl, electronColl, "PFJET_TIGHT", 20., 2.4);


//   std::vector<int> bIdxColl=GetSFBJetIdx(jetColl,"Medium");
//   std::vector<int> ljIdxColl=GetSFLJetIdx(jetColl, bIdxColl, "Medium");

//   std::vector<snu::KJet> bjetColl; for(int i=0; i<bIdxColl.size(); i++){bjetColl.push_back(jetColl.at(bIdxColl.at(i)));}
//   std::vector<snu::KJet> ljetColl; for(int i=0; i<ljIdxColl.size(); i++){ljetColl.push_back(jetColl.at(ljIdxColl.at(i)));}

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
   float trigger_sf=1.;

   if(!isData){
     trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );

     id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronColl);
     reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);

     id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muonColl);
     iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
     trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
   }


//////Basic Objects Distribution////////////////////////////////////////////////////////////////////////
//   FillHist("Basic_Nj_orig", njets, weight, 0., 10., 10);
//   FillHist("Basic_Nb_orig", nbjets, weight, 0., 10., 10);//Bjet def CVSincV2>0.89;pfKJet::CSVv2IVF_discriminator 
//   FillHist("Basic_NmuT_orig", muonColl.size(), weight, 0., 10., 10);
//   FillHist("Basic_NmuL_orig", muonLooseColl.size(), weight, 0., 10., 10);
//   FillHist("Basic_NeT_orig", electronColl.size(), weight, 0., 10., 10);
//   FillHist("Basic_NeL_orig", electronLooseColl.size(), weight, 0., 10., 10);
//   if(electronColl.size()==1) FillHist("Basic_NmuT_weT", muonColl.size(), weight, 0., 10., 10);
//   if(DoubleMu_analysis){
//     if(muonColl.size()>0) FillHist("Basic_Ptmu1_3mu_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
//     if(muonColl.size()>1) FillHist("Basic_Ptmu2_3mu_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
//     if(muonColl.size()>2) FillHist("Basic_Ptmu3_3mu_orig", muonColl.at(2).Pt(), weight, 0., 200., 200);
//   }
//   if(EMu_analysis){
//     if(muonColl.size()>0) FillHist("Basic_Ptmu1_1e2mu_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
//     if(muonColl.size()>1) FillHist("Basic_Ptmu2_1e2mu_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
//     if(electronColl.size()>0) FillHist("Basic_Pte_1e2mu_orig", electronColl.at(0).Pt(), weight, 0., 200., 200);
//   }
//   if(njets>0) FillHist("Basic_j1_Et_orig", jetColl.at(0).Et(), weight, 0., 200., 200);
//   if(njets>1) FillHist("Basic_j2_Et_orig", jetColl.at(1).Et(), weight, 0., 200., 200);
//   if(njets>2) FillHist("Basic_j3_Et_orig", jetColl.at(2).Et(), weight, 0., 200., 200);
//   if(njets>3) FillHist("Basic_j4_Et_orig", jetColl.at(3).Et(), weight, 0., 200., 200);
//   if(nbjets>0) FillHist("Basic_b1_Et_orig", bjetColl.at(0).Et(), weight, 0, 200., 200);
//   if(nbjets>1) FillHist("Basic_b2_Et_orig", bjetColl.at(1).Et(), weight, 0., 200., 200);
//   FillHist("Basic_METdist_orig", met, weight, 0., 200., 100);//It is orig because even though I havent put METcut below, if there is correlation between Nl ,Nj, Nb and MET then MET distrubution will also be affected by selection. The distribution remains the same only when the correlation is absent.

///////Event Selection&Histograms//////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////   
   if(SingleMu_analysis){

     //Check Difference of triggerCut and trigger eff
     TString analysis_trigger1="HLT_IsoMu24_v";
     TString analysis_trigger2="HLT_IsoTkMu24_v";

     bool Pass_OS2Mu=false, Pass_Mdimu=false, Pass_Trig=false;
     if(muonColl.size()==2){
       if(muonColl.at(0).Pt()>27 && muonColl.at(1).Pt()>27){
         if(SumCharge(muonColl)==0) Pass_OS2Mu=true;
       }
       if(fabs((muonColl.at(0)+muonColl.at(1)).M()-91.2)<15) Pass_Mdimu=true;
     }
     if(PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v")) Pass_Trig=true;

     if(Pass_OS2Mu && Pass_Mdimu && Pass_Trig){
       float Mmumu=(muonColl.at(0)+muonColl.at(1)).M();
       FillHist("Mmumu_NoW", Mmumu, weight, 0., 200., 200);
       FillHist("Mmumu_TrkSF", Mmumu, weight*trk_weight_mu, 0., 200., 200);
       FillHist("Mmumu_TrkIDSF", Mmumu, weight*trk_weight_mu*id_weight_mu, 0., 200., 200);
       FillHist("Mmumu_TrkIDIsoSF", Mmumu, weight*trk_weight_mu*id_weight_mu*iso_weight_mu, 0., 200., 200);
       FillHist("Mmumu_TrkIDIsoTrigSF", Mmumu, weight*trk_weight_mu*id_weight_mu*iso_weight_mu*trigger_sf, 0., 200., 200);
       FillHist("Mmumu_TrkIDIsoTrigPUSF", Mmumu, weight*trk_weight_mu*id_weight_mu*iso_weight_mu*trigger_sf*pileup_reweight, 0., 200., 200);
     }

   }
   if(DoubleMu_analysis){
     //At this point weight=MCweight*trigger_ps_weight
     //Checking Nevt Passing 2 Trigger
     if(isData){
       if(PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")) FillHist("TriggerCheck", 0., weight, 0., 2., 2);
       if(PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")) FillHist("TriggerCheck", 1., weight, 0., 2., 2);
     }
     else{
       if(PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")) FillHist("TriggerCheck", 0., weight, 0., 2., 2);
       if(PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")) FillHist("TriggerCheck", 1., weight, 0., 2., 2);
     }

     //Check Difference of triggerCut and trigger eff
     TString analysis_trigger2="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";
     if(muonColl.size()==2 && SumCharge(muonColl)==0){
       if((muonColl.at(0)+muonColl.at(1)).M()>50){
         if(isData){
           if(PassTrigger(analysis_trigger)) FillHist("Mmumu_Trig1_Trig1W", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 200., 200);
         }
         else FillHist("Mmumu_Trig1_Trig1W", (muonColl.at(0)+muonColl.at(1)).M(), weight*trigger_eff, 0., 200., 200); 

         if(isData){
           if(PassTrigger(analysis_trigger2)) FillHist("Mmumu_Trig2_Trig1W", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 200., 200);
         }
         else FillHist("Mmumu_Trig2_Trig1W", (muonColl.at(0)+muonColl.at(1)).M(), weight*trigger_eff, 0., 200., 200); 

       }
     }


     if(PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")){
       if(muonColl.size()!=2) return;
       if(SumCharge(muonColl)!=0) return;
       if((muonColl.at(0)+muonColl.at(1)).M()<50) return;
  
       //Check effect of keep true prompt
       FillHist("NmuNPmNmuP", fabs(muonColl.size()-muonPromptColl.size()), weight, 0., 50., 50);
  
       //Only Trigger Applied
       FillHist("Mmumu_Trig1_NoW", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 200., 200);
  
       //MuIDSF
       FillHist("Mmumu_Trig1_IDSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*id_weight_mu, 0., 200., 200);
  
       //MuIsoSF
       FillHist("Mmumu_Trig1_IsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*iso_weight_mu, 0., 200., 200);
  
       //PUrw
       FillHist("Mmumu_Trig1_PUSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*pileup_reweight, 0., 200., 200);
  
  
       //Combination
       //TrkSF
       FillHist("Mmumu_Trig1_TrkSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu, 0., 200., 200);
  
       //MuIDSF
       FillHist("Mmumu_Trig1_TrkIDSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu, 0., 200., 200);
  
       //MuIsoSF
       FillHist("Mmumu_Trig1_TrkIDIsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu, 0., 200., 200);
  
       //PUrw
       FillHist("Mmumu_Trig1_PUTrkIDIsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);

       //DYCounter
       if(fabs((muonColl.at(0)+muonColl.at(1)).M()-91.2)<15){
         FillHist("DYCounter_NoW", 0., weight, 0., 3., 3);
         FillHist("DYCounter_PUTrkIDIsoSF", 0., weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 3., 3);
       }
     }
     if(PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")){
       if(muonColl.size()!=2) return;
       if(SumCharge(muonColl)!=0) return;
       if((muonColl.at(0)+muonColl.at(1)).M()<50) return;
  
       //Only Trigger Applied
       FillHist("Mmumu_Trig2_NoW", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 200., 200);
  
       //MuIDSF
       FillHist("Mmumu_Trig2_IDSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*id_weight_mu, 0., 200., 200);
  
       //MuIsoSF
       FillHist("Mmumu_Trig2_IsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*iso_weight_mu, 0., 200., 200);
  
       //PUrw
       FillHist("Mmumu_Trig2_PUSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*pileup_reweight, 0., 200., 200);
  
  

       //Combination
       //TrkSF
       FillHist("Mmumu_Trig2_TrkSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu, 0., 200., 200);
  
       //MuIDSF
       FillHist("Mmumu_Trig2_TrkIDSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu, 0., 200., 200);
  
       //MuIsoSF
       FillHist("Mmumu_Trig2_TrkIDIsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu, 0., 200., 200);
  
       //PUrw
       FillHist("Mmumu_Trig2_PUTrkIDIsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);

       //DYCounter
       if(fabs((muonColl.at(0)+muonColl.at(1)).M()-91.2)<15){
         FillHist("DYCounter_NoW", 1., weight, 0., 3., 3);
         FillHist("DYCounter_PUTrkIDIsoSF", 1., weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 3., 3);
       }

     }
     if(PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") || PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")){
       //DYCounter
       if(fabs((muonColl.at(0)+muonColl.at(1)).M()-91.2)<15){
         FillHist("DYCounter_NoW", 2., weight, 0., 3., 3);
         FillHist("DYCounter_PUTrkIDIsoSF", 2., weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 3., 3);
       }
     }
     if(isData){
       if(PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")){
         if(muonColl.size()!=2) return;
         if(SumCharge(muonColl)!=0) return;
         if((muonColl.at(0)+muonColl.at(1)).M()<50) return;

         //Only Trigger Applied
         FillHist("Mmumu_Trig1W_NoW", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 200., 200);
    
         //MuIDSF
         FillHist("Mmumu_Trig1W_IDSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*id_weight_mu, 0., 200., 200);
    
         //MuIsoSF
         FillHist("Mmumu_Trig1W_IsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*iso_weight_mu, 0., 200., 200);
    
         //PUrw
         FillHist("Mmumu_Trig1W_PUSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*pileup_reweight, 0., 200., 200);
    
    
         //Combination
         //TrkSF
         FillHist("Mmumu_Trig1W_TrkSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu, 0., 200., 200);
    
         //MuIDSF
         FillHist("Mmumu_Trig1W_TrkIDSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu, 0., 200., 200);
    
         //MuIsoSF
         FillHist("Mmumu_Trig1W_TrkIDIsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu, 0., 200., 200);
    
         //PUrw
         FillHist("Mmumu_Trig1W_PUTrkIDIsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);
       }
       if(PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")){
         if(muonColl.size()!=2) return;
         if(SumCharge(muonColl)!=0) return;
         if((muonColl.at(0)+muonColl.at(1)).M()<50) return;

         //Only Trigger Applied
         FillHist("Mmumu_Trig2W_NoW", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 200., 200);
    
         //MuIDSF
         FillHist("Mmumu_Trig2W_IDSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*id_weight_mu, 0., 200., 200);
    
         //MuIsoSF
         FillHist("Mmumu_Trig2W_IsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*iso_weight_mu, 0., 200., 200);
    
         //PUrw
         FillHist("Mmumu_Trig2W_PUSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*pileup_reweight, 0., 200., 200);
    
    
         //Combination
         //TrkSF
         FillHist("Mmumu_Trig2W_TrkSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu, 0., 200., 200);
    
         //MuIDSF
         FillHist("Mmumu_Trig2W_TrkIDSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu, 0., 200., 200);
    
         //MuIsoSF
         FillHist("Mmumu_Trig2W_TrkIDIsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu, 0., 200., 200);
    
         //PUrw
         FillHist("Mmumu_Trig2W_PUTrkIDIsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);
       }
     }
     else{
         if(muonColl.size()!=2) return;
         if(SumCharge(muonColl)!=0) return;
         if((muonColl.at(0)+muonColl.at(1)).M()<50) return;

         weight*=trigger_eff;
         //Only Trigger Applied
         FillHist("Mmumu_Trig1W_NoW", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 200., 200);
    
         //MuIDSF
         FillHist("Mmumu_Trig1W_IDSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*id_weight_mu, 0., 200., 200);
    
         //MuIsoSF
         FillHist("Mmumu_Trig1W_IsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*iso_weight_mu, 0., 200., 200);
    
         //PUrw
         FillHist("Mmumu_Trig1W_PUSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*pileup_reweight, 0., 200., 200);
    
    
         //Combination
         //TrkSF
         FillHist("Mmumu_Trig1W_TrkSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu, 0., 200., 200);
    
         //MuIDSF
         FillHist("Mmumu_Trig1W_TrkIDSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu, 0., 200., 200);
    
         //MuIsoSF
         FillHist("Mmumu_Trig1W_TrkIDIsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu, 0., 200., 200);
    
         //PUrw
         FillHist("Mmumu_Trig1W_PUTrkIDIsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);

         //Only Trigger Applied
         FillHist("Mmumu_Trig2W_NoW", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 200., 200);
    
         //MuIDSF
         FillHist("Mmumu_Trig2W_IDSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*id_weight_mu, 0., 200., 200);
    
         //MuIsoSF
         FillHist("Mmumu_Trig2W_IsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*iso_weight_mu, 0., 200., 200);
    
         //PUrw
         FillHist("Mmumu_Trig2W_PUSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*pileup_reweight, 0., 200., 200);
    
    
         //Combination
         //TrkSF
         FillHist("Mmumu_Trig2W_TrkSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu, 0., 200., 200);
    
         //MuIDSF
         FillHist("Mmumu_Trig2W_TrkIDSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu, 0., 200., 200);
    
         //MuIsoSF
         FillHist("Mmumu_Trig2W_TrkIDIsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu, 0., 200., 200);
    
         //PUrw
         FillHist("Mmumu_Trig2W_PUTrkIDIsoSF", (muonColl.at(0)+muonColl.at(1)).M(), weight*trk_weight_mu*id_weight_mu*iso_weight_mu*pileup_reweight, 0., 200., 200);
     }
    
   }
   if(DoubleEle_analysis){
     if(isData){ if(!PassTrigger(analysis_trigger)) return; }
     if(electronColl.size()!=2) return;
     if(electronColl.at(0).Charge()==electronColl.at(1).Charge()) return;
     if(electronColl.at(0).Pt()<25) return;
     if(electronColl.at(1).Pt()<15) return;
     if((electronColl.at(0)+electronColl.at(1)).M()<50) return;

     if(fabs((electronColl.at(0)+electronColl.at(1)).M()-91.2)<15){
       FillHist("DYCounter_NoW_NoSF", 0., weight, 0., 1., 1);
       FillHist("DYCounter_NoW_IDPUSF", 0., weight*id_weight_ele*pileup_reweight, 0., 1., 1);
       FillHist("DYCounter_TrigW_NoSF", 0., weight*trigger_eff, 0., 1., 1);
       FillHist("DYCounter_TrigW_IDSF", 0., weight*trigger_eff*id_weight_ele, 0., 1., 1);
       FillHist("DYCounter_TrigW_IDPUSF", 0., weight*trigger_eff*id_weight_ele*pileup_reweight, 0., 1., 1);
     }
     FillHist("Mee_TrigW_NoSF", (electronColl.at(0)+electronColl.at(1)).M(), weight*trigger_eff, 0., 200., 200);
     FillHist("Mee_TrigW_IDSF", (electronColl.at(0)+electronColl.at(1)).M(), weight*trigger_eff*id_weight_ele, 0., 200., 200);
     FillHist("Mee_TrigW_IDPUSF", (electronColl.at(0)+electronColl.at(1)).M(), weight*trigger_eff*id_weight_ele*pileup_reweight, 0., 200., 200);

   }


/////////////////////////////////////////////////////////////////////////////////// 


return;
}// End of execute event loop
  


void Jan2017_3l4j_DYCheck::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Jan2017_3l4j_DYCheck::BeginCycle() throw( LQError ){
  
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

Jan2017_3l4j_DYCheck::~Jan2017_3l4j_DYCheck() {
  
  Message("In Jan2017_3l4j_DYCheck Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Jan2017_3l4j_DYCheck::FillCutFlow(TString cut, float weight){
  
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



void Jan2017_3l4j_DYCheck::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Jan2017_3l4j_DYCheck::MakeHistograms(){
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
  // *  Remove//Overide this Jan2017_3l4j_DYCheckCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Jan2017_3l4j_DYCheck::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
