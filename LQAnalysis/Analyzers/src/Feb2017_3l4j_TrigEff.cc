// $Id: Feb2017_3l4j_TrigEff.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQFeb2017_3l4j_TrigEff Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

 // Local includes
 #include "Feb2017_3l4j_TrigEff.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 // Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Feb2017_3l4j_TrigEff);

 Feb2017_3l4j_TrigEff::Feb2017_3l4j_TrigEff() : AnalyzerCore(), out_muons(0) {

   SetLogName("Feb2017_3l4j_TrigEff");
   Message("In Feb2017_3l4j_TrigEff constructor", INFO);
   InitialiseAnalysis();
 }


 void Feb2017_3l4j_TrigEff::InitialiseAnalysis() throw( LQError ) {
   
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

void Feb2017_3l4j_TrigEff::ExecuteEvents()throw( LQError ){

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


   //bool TriMu_analysis=true, EMuMu_analysis=false;
   bool TriMu_analysis=false, EMuMu_analysis=true;


   //Trigers
   std::vector<TString> triggerslist, triggerlist1, triggerlist2;//  ListTriggersAvailable();

   TString analysis_trigger;
   if(EMuMu_analysis) {analysis_trigger="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";}
   else if(TriMu_analysis) analysis_trigger="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";
   //if(!PassTrigger(analysis_trigger)) return; //Not for now..

   //PassTrigger Uses BeginWith so don't stick to exact name

   float trigger_ps_weight=1, weight_trigger_sf=1, trigger_eff=1;
   //trigger_ps_weight=WeightByTrigger(analysis_trigger, TargetLumi);
   trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);



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
     eventbase->GetMuonSel()->SetPt(5.);                    eventbase->GetMuonSel()->SetEta(2.4);
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
     eventbase->GetElectronSel()->SetPt(20.);                eventbase->GetElectronSel()->SetEta(2.5);
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
     //trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );

     id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronColl);
     reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);

     id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muonColl);
     iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
     trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
   }


//////Basic Objects Distribution////////////////////////////////////////////////////////////////////////
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
   if(njets>0) FillHist("Basic_j1_Et_orig", jetColl.at(0).Et(), weight, 0., 200., 200);
   if(njets>1) FillHist("Basic_j2_Et_orig", jetColl.at(1).Et(), weight, 0., 200., 200);
   if(njets>2) FillHist("Basic_j3_Et_orig", jetColl.at(2).Et(), weight, 0., 200., 200);
   if(njets>3) FillHist("Basic_j4_Et_orig", jetColl.at(3).Et(), weight, 0., 200., 200);
   if(nbjets>0) FillHist("Basic_b1_Et_orig", bjetColl.at(0).Et(), weight, 0, 200., 200);
   if(nbjets>1) FillHist("Basic_b2_Et_orig", bjetColl.at(1).Et(), weight, 0., 200., 200);
   FillHist("Basic_METdist_orig", met, weight, 0., 200., 100);//It is orig because even though I havent put METcut below, if there is correlation between Nl ,Nj, Nb and MET then MET distrubution will also be affected by selection. The distribution remains the same only when the correlation is absent.

///////Event Selection&Histograms//////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////   

     bool Pass_Mu17TkMu8Dz  = PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
     bool Pass_Mu17Mu8Dz    = PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
     bool Pass_IsoMu24      = PassTrigger("HLT_IsoMu24_v");
     bool Pass_IsoTkMu24    = PassTrigger("HLT_IsoTkMu24_v");
     bool Pass_Mu12Mu8Mu5   = PassTrigger("HLT_TripleMu_12_10_5_v");

     bool Pass_Ele27        = PassTrigger("HLT_Ele27_WPTight_Gsf_v");
     //bool Pass_Ele27        = PassTrigger("HLT_Ele32_eta2p1_WPTight_Gsf_v");
//     bool Pass_Ele23Mu8or12 = (  PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
//                              || PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v")
//                              || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
//                              || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") );
     bool Pass_Ele23Mu8or12 = (  PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
                              || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") );

     bool Pass_DiMu9Ele9    = PassTrigger("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v");


     weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight;
     FillHist("Count_NoCut", 0., weight, 0., 1., 1);

   if(TriMu_analysis){

     bool Pass_LepMinSel=false, Pass_LepMeanSel=false, Pass_SiglLepSel=false;
     
     if(muonColl.size()==3){
       if(fabs(SumCharge(muonColl))==1){
         if(muonColl.at(0).Pt()>20 && muonColl.at(1).Pt()>10) Pass_LepMinSel=true;
         if(muonColl.at(0).Pt()>20 && muonColl.at(2).Pt()>10) Pass_LepMeanSel=true;
         if(muonColl.at(0).Pt()>26)                           Pass_SiglLepSel=true;
       }
     }

     //Only Trigger Skim
     //Minimum lepton requirements needed, because events with leptons out of acceptance have no meaning to count
     //And there are also cuts that I will apply in any case, determined from Genlevel study
     if(Pass_LepMinSel){
       if(Pass_Mu17TkMu8Dz) FillHist("Count_EachTrigMinSel_3mu", 0., weight, 0., 10., 10);
       if(Pass_Mu17Mu8Dz)   FillHist("Count_EachTrigMinSel_3mu", 1., weight, 0., 10., 10);
       if(Pass_Mu17TkMu8Dz || Pass_Mu17Mu8Dz) FillHist("Count_EachTrigMinSel_3mu", 2., weight, 0., 10., 10);
     }
     if(Pass_SiglLepSel){
       if(Pass_IsoMu24)     FillHist("Count_EachTrigMinSel_3mu", 3., weight, 0., 10., 10);
      // if(Pass_IsoTkMu24)   FillHist("Count_EachTrigMinSel_3mu", 4., weight, 0., 10., 10);
      // if(Pass_IsoMu24 || Pass_IsoTkMu24) FillHist("Count_EachTrigMinSel_3mu", 5., weight, 0., 10., 10);
     }
     if(Pass_LepMinSel){
      // if(Pass_Mu12Mu8Mu5)  FillHist("Count_EachTrigMinSel_3mu", 6., weight, 0., 10., 10);
       if(Pass_Mu12Mu8Mu5)  FillHist("Count_EachTrigMinSel_3mu", 4., weight, 0., 10., 10);
     }

     if(Pass_LepMeanSel){
       if(Pass_Mu17TkMu8Dz) FillHist("Count_EachTrigMeanSel_3mu", 0., weight, 0., 10., 10);
       if(Pass_Mu17Mu8Dz)   FillHist("Count_EachTrigMeanSel_3mu", 1., weight, 0., 10., 10);
       if(Pass_Mu17TkMu8Dz || Pass_Mu17Mu8Dz) FillHist("Count_EachTrigMeanSel_3mu", 2., weight, 0., 10., 10);
       if(Pass_SiglLepSel && Pass_IsoMu24)     FillHist("Count_EachTrigMeanSel_3mu", 3., weight, 0., 10., 10);
       //if(Pass_SiglLepSel && Pass_IsoTkMu24)   FillHist("Count_EachTrigMeanSel_3mu", 4., weight, 0., 10., 10);
       //if(Pass_SiglLepSel && (Pass_IsoMu24 || Pass_IsoTkMu24)) FillHist("Count_EachTrigMeanSel_3mu", 5., weight, 0., 10., 10);
       //if(Pass_Mu12Mu8Mu5)  FillHist("Count_EachTrigMeanSel_3mu", 6., weight, 0., 10., 10);
       if(Pass_Mu12Mu8Mu5)  FillHist("Count_EachTrigMeanSel_3mu", 4., weight, 0., 10., 10);
     }


     //Combinations
     if(Pass_LepMinSel){
       if(Pass_Mu17TkMu8Dz || Pass_Mu17Mu8Dz ) FillHist("Count_CombTrigMinSel_3mu", 0., weight, 0., 10., 10);
       if(Pass_Mu17TkMu8Dz || Pass_Mu17Mu8Dz || Pass_IsoMu24 || Pass_IsoTkMu24){
         FillHist("Count_CombTrigMinSel_3mu", 1., weight, 0., 10., 10);
       }
       if(Pass_Mu17TkMu8Dz || Pass_Mu17Mu8Dz || Pass_IsoMu24 || Pass_IsoTkMu24 || Pass_Mu12Mu8Mu5 ){
         FillHist("Count_CombTrigMinSel_3mu", 2., weight, 0., 10., 10);
       }
     }
     if(Pass_LepMeanSel){
       if(Pass_Mu17TkMu8Dz || Pass_Mu17Mu8Dz ) FillHist("Count_CombTrigMeanSel_3mu", 0., weight, 0., 10., 10);
       if(Pass_Mu17TkMu8Dz || Pass_Mu17Mu8Dz || Pass_IsoMu24 || Pass_IsoTkMu24){
         FillHist("Count_CombTrigMeanSel_3mu", 1., weight, 0., 10., 10);
       }
       if(Pass_Mu17TkMu8Dz || Pass_Mu17Mu8Dz || Pass_IsoMu24 || Pass_IsoTkMu24 || Pass_Mu12Mu8Mu5 ){
         FillHist("Count_CombTrigMeanSel_3mu", 2., weight, 0., 10., 10);
       }
     }


   }//End of TriMu Loop

   if(EMuMu_analysis){

     bool Pass_EMuMinSel=false, Pass_DiMuMinSel=false, Pass_EMuMeanSel=false, Pass_EMuMuSel=false, Pass_SiglESel=false, Pass_SiglMuSel=false, Pass_MuLegMeanSel=false;
     
     if(muonColl.size()==2 && electronColl.size()==1){
       if(fabs(SumCharge(muonColl))==0){
         if(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>14) Pass_EMuMinSel=true;//EMu Trig Base approach
         if(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>14 && muonColl.at(1).Pt()>10) Pass_EMuMeanSel=true;//EMu Trig Base approach
         if(electronColl.at(0).Pt()>20 && muonColl.at(1).Pt()>11) Pass_EMuMuSel=true;//EMuMu Trig Base approach
         if(electronColl.at(0).Pt()>30 ) Pass_SiglESel=true;//EMuMu Trig Base approach
         if(muonColl.at(0).Pt()>26 )     Pass_SiglMuSel=true;//EMuMu Trig Base approach
         if(muonColl.at(0).Pt()>20 && muonColl.at(1).Pt()>10) Pass_DiMuMinSel=true;//EMuMu Trig Base approach
         if(muonColl.at(1).Pt()>10) Pass_MuLegMeanSel=true;//EMuMu Trig Base approach
       }
     }


     //Only Trigger Skim
     //Minimum Cuts
     if(Pass_EMuMinSel   && Pass_Ele23Mu8or12)                    FillHist("Count_EachTrigMinSel_1e2mu", 0., weight, 0., 10., 10);
     if(Pass_EMuMuSel    && Pass_DiMu9Ele9)                       FillHist("Count_EachTrigMinSel_1e2mu", 1., weight, 0., 10., 10);
     if(Pass_SiglESel    && Pass_Ele27)                           FillHist("Count_EachTrigMinSel_1e2mu", 2., weight, 0., 10., 10);
     if(Pass_SiglMuSel   && Pass_IsoMu24)                         FillHist("Count_EachTrigMinSel_1e2mu", 3., weight, 0., 10., 10);
     if(Pass_DiMuMinSel  && (Pass_Mu17TkMu8Dz || Pass_Mu17Mu8Dz)) FillHist("Count_EachTrigMinSel_1e2mu", 4., weight, 0., 10., 10);


     //Reasonable Cuts
     if(Pass_EMuMeanSel && Pass_Ele23Mu8or12)                   FillHist("Count_EachTrigMeanSel_1e2mu", 0., weight, 0., 10., 10);
     if(Pass_EMuMuSel   && Pass_DiMu9Ele9)                      FillHist("Count_EachTrigMeanSel_1e2mu", 1., weight, 0., 10., 10);
     if(Pass_SiglESel && Pass_MuLegMeanSel && Pass_Ele27)       FillHist("Count_EachTrigMeanSel_1e2mu", 2., weight, 0., 10., 10);
     if(Pass_SiglMuSel && Pass_MuLegMeanSel && (Pass_IsoMu24 || Pass_IsoTkMu24))     FillHist("Count_EachTrigMeanSel_1e2mu", 3., weight, 0., 10., 10);
     if(Pass_DiMuMinSel && (Pass_Mu17TkMu8Dz || Pass_Mu17Mu8Dz))FillHist("Count_EachTrigMeanSel_1e2mu", 4., weight, 0., 10., 10);
     



     //Combinations Approach 1
     if(Pass_EMuMinSel){
       if(Pass_Ele23Mu8or12)                 FillHist("Count_CombTrigMinSel1_1e2mu", 0., weight, 0., 10., 10);
       if(Pass_Ele23Mu8or12 || Pass_Ele27)   FillHist("Count_CombTrigMinSel1_1e2mu", 1., weight, 0., 10., 10);
       if(Pass_Ele23Mu8or12 || Pass_IsoMu24) FillHist("Count_CombTrigMinSel1_1e2mu", 2., weight, 0., 10., 10);
       if(Pass_Ele23Mu8or12 || Pass_Ele27 || Pass_IsoMu24){
         FillHist("Count_CombTrigMinSel1_1e2mu", 3., weight, 0., 10., 10);
       }
     }
     if(Pass_EMuMeanSel){
       if(Pass_Ele23Mu8or12)                 FillHist("Count_CombTrigMeanSel1_1e2mu", 0., weight, 0., 10., 10);
       if(Pass_Ele23Mu8or12 || Pass_Ele27)   FillHist("Count_CombTrigMeanSel1_1e2mu", 1., weight, 0., 10., 10);
       if(Pass_Ele23Mu8or12 || Pass_IsoMu24) FillHist("Count_CombTrigMeanSel1_1e2mu", 2., weight, 0., 10., 10);
       if(Pass_Ele23Mu8or12 || Pass_Ele27 || Pass_IsoMu24){
         FillHist("Count_CombTrigMeanSel1_1e2mu", 3., weight, 0., 10., 10);
       }
     }


     //Combinations Approach 2
     if(Pass_EMuMuSel){
       if(Pass_DiMu9Ele9)                    FillHist("Count_CombTrigOnSel2_1e2mu", 0., weight, 0., 10., 10);
       if(Pass_DiMu9Ele9 || Pass_Ele27)      FillHist("Count_CombTrigOnSel2_1e2mu", 1., weight, 0., 10., 10);
       if(Pass_DiMu9Ele9 || Pass_IsoMu24)    FillHist("Count_CombTrigOnSel2_1e2mu", 2., weight, 0., 10., 10);
       if(Pass_DiMu9Ele9 || Pass_Ele27 || Pass_IsoMu24){
         FillHist("Count_CombTrigOnSel2_1e2mu", 3., weight, 0., 10., 10);
       }
     }


   }//End of EMuMu Loop



/////////////////////////////////////////////////////////////////////////////////// 


return;
}// End of execute event loop
  


void Feb2017_3l4j_TrigEff::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Feb2017_3l4j_TrigEff::BeginCycle() throw( LQError ){
  
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

Feb2017_3l4j_TrigEff::~Feb2017_3l4j_TrigEff() {
  
  Message("In Feb2017_3l4j_TrigEff Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Feb2017_3l4j_TrigEff::FillCutFlow(TString cut, float weight){
  
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



void Feb2017_3l4j_TrigEff::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Feb2017_3l4j_TrigEff::MakeHistograms(){
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

  //Analysis Histograms
  AnalyzerCore::MakeHistograms("Count_EachTrigMinSel_3mu", 10, 0., 10.);
   GetHist("Count_EachTrigMinSel_3mu")->GetXaxis()->SetBinLabel(1,"Mu17TkMu8Dz");
   GetHist("Count_EachTrigMinSel_3mu")->GetXaxis()->SetBinLabel(2,"Mu17Mu8Dz");
   GetHist("Count_EachTrigMinSel_3mu")->GetXaxis()->SetBinLabel(3,"Mu17(Tk)Mu8Dz");
   GetHist("Count_EachTrigMinSel_3mu")->GetXaxis()->SetBinLabel(4,"IsoMu24");
   //GetHist("Count_EachTrigMinSel_3mu")->GetXaxis()->SetBinLabel(5,"IsoTkMu24");
   //GetHist("Count_EachTrigMinSel_3mu")->GetXaxis()->SetBinLabel(6,"Iso(Tk)Mu24");
   //GetHist("Count_EachTrigMinSel_3mu")->GetXaxis()->SetBinLabel(7,"Mu12_8_5");
   GetHist("Count_EachTrigMinSel_3mu")->GetXaxis()->SetBinLabel(5,"Mu12_8_5");

  AnalyzerCore::MakeHistograms("Count_EachTrigMeanSel_3mu", 10, 0., 10.);
   GetHist("Count_EachTrigMeanSel_3mu")->GetXaxis()->SetBinLabel(1,"Mu17TkMu8Dz");
   GetHist("Count_EachTrigMeanSel_3mu")->GetXaxis()->SetBinLabel(2,"Mu17Mu8Dz");
   GetHist("Count_EachTrigMeanSel_3mu")->GetXaxis()->SetBinLabel(3,"Mu17(Tk)Mu8Dz");
   GetHist("Count_EachTrigMeanSel_3mu")->GetXaxis()->SetBinLabel(4,"IsoMu24");
   //GetHist("Count_EachTrigMeanSel_3mu")->GetXaxis()->SetBinLabel(5,"IsoTkMu24");
   //GetHist("Count_EachTrigMeanSel_3mu")->GetXaxis()->SetBinLabel(6,"Iso(Tk)Mu24");
   //GetHist("Count_EachTrigMeanSel_3mu")->GetXaxis()->SetBinLabel(7,"Mu12_8_5");
   GetHist("Count_EachTrigMeanSel_3mu")->GetXaxis()->SetBinLabel(5,"Mu12_8_5");

  AnalyzerCore::MakeHistograms("Count_CombTrigMinSel_3mu", 10, 0., 10.);
   GetHist("Count_CombTrigMinSel_3mu")->GetXaxis()->SetBinLabel(1,"Mu17(Tk)Mu8Dz");
   GetHist("Count_CombTrigMinSel_3mu")->GetXaxis()->SetBinLabel(2,"OR_Iso(Tk)Mu24");
   GetHist("Count_CombTrigMinSel_3mu")->GetXaxis()->SetBinLabel(3,"OR_Mu12_8_5");

  AnalyzerCore::MakeHistograms("Count_CombTrigMeanSel_3mu", 10, 0., 10.);
   GetHist("Count_CombTrigMeanSel_3mu")->GetXaxis()->SetBinLabel(1,"Mu17(Tk)Mu8Dz");
   GetHist("Count_CombTrigMeanSel_3mu")->GetXaxis()->SetBinLabel(2,"OR_Iso(Tk)Mu24");
   GetHist("Count_CombTrigMeanSel_3mu")->GetXaxis()->SetBinLabel(3,"OR_Mu12_8_5");


  //EMuMu
  AnalyzerCore::MakeHistograms("Count_EachTrigMinSel_1e2mu", 10, 0., 10.);
   GetHist("Count_EachTrigMinSel_1e2mu")->GetXaxis()->SetBinLabel(1,"Ele23Mu8or12");
   GetHist("Count_EachTrigMinSel_1e2mu")->GetXaxis()->SetBinLabel(2,"DiMu9Ele9");
   GetHist("Count_EachTrigMinSel_1e2mu")->GetXaxis()->SetBinLabel(3,"Ele27");
   GetHist("Count_EachTrigMinSel_1e2mu")->GetXaxis()->SetBinLabel(4,"IsoMu24");
   GetHist("Count_EachTrigMinSel_1e2mu")->GetXaxis()->SetBinLabel(5,"Mu17(Tk)Mu8Dz");


  AnalyzerCore::MakeHistograms("Count_EachTrigMeanSel_1e2mu", 10, 0., 10.);
   GetHist("Count_EachTrigMeanSel_1e2mu")->GetXaxis()->SetBinLabel(1,"Ele23Mu8or12");
   GetHist("Count_EachTrigMeanSel_1e2mu")->GetXaxis()->SetBinLabel(2,"DiMu9Ele9");
   GetHist("Count_EachTrigMeanSel_1e2mu")->GetXaxis()->SetBinLabel(3,"Ele27");
   GetHist("Count_EachTrigMeanSel_1e2mu")->GetXaxis()->SetBinLabel(4,"IsoMu24");
   GetHist("Count_EachTrigMeanSel_1e2mu")->GetXaxis()->SetBinLabel(5,"Mu17(Tk)Mu8Dz");



  AnalyzerCore::MakeHistograms("Count_CombTrigMinSel1_1e2mu", 10, 0., 10.);
   GetHist("Count_CombTrigMinSel1_1e2mu")->GetXaxis()->SetBinLabel(1,"Ele23Mu8or12");
   GetHist("Count_CombTrigMinSel1_1e2mu")->GetXaxis()->SetBinLabel(2,"OR_Ele27");
   GetHist("Count_CombTrigMinSel1_1e2mu")->GetXaxis()->SetBinLabel(3,"OR_IsoMu24");
   GetHist("Count_CombTrigMinSel1_1e2mu")->GetXaxis()->SetBinLabel(4,"OR_Ele27||IsoMu24");

  AnalyzerCore::MakeHistograms("Count_CombTrigMeanSel1_1e2mu", 10, 0., 10.);
   GetHist("Count_CombTrigMeanSel1_1e2mu")->GetXaxis()->SetBinLabel(1,"Ele23Mu8or12");
   GetHist("Count_CombTrigMeanSel1_1e2mu")->GetXaxis()->SetBinLabel(2,"OR_Ele27");
   GetHist("Count_CombTrigMeanSel1_1e2mu")->GetXaxis()->SetBinLabel(3,"OR_IsoMu24");
   GetHist("Count_CombTrigMeanSel1_1e2mu")->GetXaxis()->SetBinLabel(4,"OR_Ele27||IsoMu24");


  AnalyzerCore::MakeHistograms("Count_CombTrigOnSel2_1e2mu", 10, 0., 10.);
   GetHist("Count_CombTrigOnSel2_1e2mu")->GetXaxis()->SetBinLabel(1,"DiMu9Ele9");
   GetHist("Count_CombTrigOnSel2_1e2mu")->GetXaxis()->SetBinLabel(2,"OR_Ele27");
   GetHist("Count_CombTrigOnSel2_1e2mu")->GetXaxis()->SetBinLabel(3,"OR_IsoMu24");
   GetHist("Count_CombTrigOnSel2_1e2mu")->GetXaxis()->SetBinLabel(4,"OR_Ele27||IsoMu24");


  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Feb2017_3l4j_TrigEffCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Feb2017_3l4j_TrigEff::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
