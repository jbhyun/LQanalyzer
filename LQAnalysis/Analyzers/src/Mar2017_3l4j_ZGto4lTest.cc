// $Id: Mar2017_3l4j_ZGto4lTest.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQMar2017_3l4j_ZGto4lTest Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Mar2017_3l4j_ZGto4lTest.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Mar2017_3l4j_ZGto4lTest);

 Mar2017_3l4j_ZGto4lTest::Mar2017_3l4j_ZGto4lTest() : AnalyzerCore(), out_muons(0) {

   SetLogName("Mar2017_3l4j_ZGto4lTest");
   Message("In Mar2017_3l4j_ZGto4lTest constructor", INFO);
   InitialiseAnalysis();
 }


 void Mar2017_3l4j_ZGto4lTest::InitialiseAnalysis() throw( LQError ) {
   
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

void Mar2017_3l4j_ZGto4lTest::ExecuteEvents()throw( LQError ){

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
   if (!k_isdata) { pileup_reweight=eventbase->GetEvent().PileUpWeight(); }

   //Numbet of Vertex PUreweight
   FillHist("Nvtx_nocut_PURW", eventbase->GetEvent().nVertices(), weight*pileup_reweight, 0., 50., 50);


   //Total Event(MCweight+PUrw)////////////////
   FillCutFlow("NoCut", weight*pileup_reweight);


   bool TriMu_analysis=false, EMuMu_analysis=false, DiMuon_analysis=false, DiEle_analysis=false;
   //bool TriMu_analysis=true, EMuMu_analysis=false, DiMuon_analysis=false, DiEle_analysis=false;
   //bool TriMu_analysis=false, EMuMu_analysis=true, DiMuon_analysis=false, DiEle_analysis=false;
   //bool TriMu_analysis=false, EMuMu_analysis=false, DiMuon_analysis=true, DiEle_analysis=false;
   //bool TriMu_analysis=false, EMuMu_analysis=false, DiMuon_analysis=false, DiEle_analysis=true;
   bool ZGtest=true;

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
   if(EMuMu_analysis){
     //if( PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v") ) Pass_Trigger=true;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
          || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
          || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger=true;
   }
   else if(TriMu_analysis){
     //if( PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v") ) Pass_Trigger=true;
     if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ) Pass_Trigger=true;
   }
   else if(DiMuon_analysis){
     //if( PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v") ) Pass_Trigger=true;
     if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ) Pass_Trigger=true;
   }
   else if(DiEle_analysis){
     if( PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger=true;
   }

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
  
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.25);//POG WP L
   std::vector<snu::KMuon> muonVetoColl; eventbase->GetMuonSel()->Selection(muonVetoColl, true);//(muonColl, bool RochCorr, bool debug)
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
   std::vector<snu::KMuon> muonColl;      if(k_running_nonprompt){ muonColl=muonLooseColl;} else{ muonColl=muonTightColl;}
   //std::vector<snu::KMuon> muonHNColl; eventbase->GetMuonSel()->SelectMuons(muonHNColl, "MUON_HN_TRI_TIGHT", 10., 2.4);
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
     eventbase->GetElectronSel()->SetCheckCharge(true);
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(20.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.2);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetCheckCharge(true);
   std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->Selection(electronColl);
   std::vector<snu::KElectron> electronNull;
   //std::vector<snu::KElectron> electronHNColl; eventbase->GetElectronSel()->SelectElectrons(electronHNColl, "ELECTRON_HN_TIGHT", 20., 2.4);

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(30.);                     eventbase->GetJetSel()->SetEta(2.4);
     //eventbase->GetJetSel()->SetPt(20.);                     eventbase->GetJetSel()->SetEta(2.4);
     eventbase->GetJetSel()->SetPileUpJetID(true,"Loose");
     bool LeptonVeto=true;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonTightColl, electronColl);
   //std::vector<snu::KJet> jetLooseColl; eventbase->GetJetSel()->SelectJets(jetLooseColl, muonColl, electronColl, "PFJET_LOOSE", 20., 2.4);

   //Method to apply 2a method SF
   //std::vector<int> bIdxColl  = GetSFBJetIdx(jetColl,"Medium");
   //std::vector<int> ljIdxColl = GetSFLJetIdx(jetColl, bIdxColl, "Medium");
   //std::vector<snu::KJet> bjetColl2a; for(int i=0; i<bIdxColl.size(); i++) {bjetColl2a.push_back(jetColl.at(bIdxColl.at(i)));}
   //std::vector<snu::KJet> ljetColl2a; for(int i=0; i<ljIdxColl.size(); i++){ljetColl2a.push_back(jetColl.at(ljIdxColl.at(i)));}

   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");
   std::vector<snu::KJet> ljetColl = SelLightJets(jetColl, "Medium");

   double met = eventbase->GetEvent().MET();
   double met_x = eventbase->GetEvent().MET()*TMath::Cos(eventbase->GetEvent().METPhi());
   double met_y = eventbase->GetEvent().MET()*TMath::Sin(eventbase->GetEvent().METPhi());
   double Pzv,Pzv1, Pzv2;
   snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));
//   snu::KParticle v[4]; v[0].SetPx(met_x); v[0].SetPy(met_y);
//                        v[1].SetPx(met_x); v[1].SetPy(met_y);

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
       if(k_running_nonprompt && TriMu_analysis){
         fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", muonLooseColl.size(), electronNull, "ELECTRON_HN_TIGHT", electronNull.size());
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
     if( !(electronColl.size()==1 && muonColl.size()==2) ) return;
     if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>15) ) return;
     if( SumCharge(muonColl)!=0 ) return;

   }
   if(TriMu_analysis){

     //Step1 : 3mu +PTCut+ SSS veto
     //if( !(muonColl.at(0).Pt()>27) ) return; //SingleMuon Trig Case
     if( !(muonColl.size()==3 && electronColl.size()==0) ) return;
     if( !(muonColl.at(0).Pt()>20 && muonColl.at(2).Pt()>10) ) return;
     if( fabs(SumCharge(muonColl))!=1 ) return;
     FillCutFlow("3lCut", weight);

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

     //CR - TriLep + b veto
//     if(bjetColl.size()!=0) goto Reg_3ljb;

//     FillCutFlow("bVeto", weight);
     //3lep Z peak test : Fake Control Check
     if(IdxZCandLead!=-1 && IdxZCandSubl!=-1){
       FillHist("Mmumu_Zwin10_3mu", (muonColl.at(IdxZCandLead)+muonColl.at(IdxZCandSubl)).M(), weight, 0., 200., 200);
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
       FillHist("MET_3lOSM12Cut_3mu", met, weight, 0., 200., 200); 
       FillHist("Nb_3lOSM12Cut_3mu", bjetColl.size(), weight, 0., 10., 10);
       FillHist("dRmumu_OSSS1_3lOSM12Cut_3mu", muonColl.at(IdxOS).DeltaR(muonColl.at(IdxSS1)), weight, 0., 5., 500);
       FillHist("dRmumu_OSSS2_3lOSM12Cut_3mu", muonColl.at(IdxOS).DeltaR(muonColl.at(IdxSS2)), weight, 0., 5., 500);
       FillHist("Mmumu_OSSS1_3lOSM12Cut_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS1)).M(), weight, 0., 200., 200);
       FillHist("Mmumu_OSSS2_3lOSM12Cut_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS2)).M(), weight, 0., 200., 200);
       FillHist("M3mu_3lOSM12Cut_3mu", (muonColl.at(IdxOS)+muonColl.at(IdxSS1)+muonColl.at(IdxSS2)).M(), weight, 0., 500., 500);

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

       if(M3l>70 && M3l<100){
         FillHist("dRmuZmuNonZS_3lOSM12M3lCut_3mu", dRS, weight, 0., 5., 500);
         FillHist("dRmuZmuNonZL_3lOSM12M3lCut_3mu", dRL, weight, 0., 5., 500);
         FillHist("dRmuZ_3lOSM12M3lCut_3mu", dRZmu, weight, 0., 5., 500);
         FillHist("AbsdPhimuZmuNonZS_3lOSM12M3lCut_3mu", dPhiS, weight, 0., 3.15, 200);
         FillHist("AbsdPhimuZmuNonZL_3lOSM12M3lCut_3mu", dPhiL, weight, 0., 3.15, 200);
         FillHist("AbsdPhimuZ_3lOSM12M3lCut_3mu", dPhiZmu, weight, 0., 3.15, 200);
       }

     }

     //ttZ peak test : ttZ cross section && Fake Control with b present
//     Reg_3ljb: ;

   }
   if(ZGtest){
     if(k_sample_name.Contains("ZZTo4L")){
       //Step1 : 3mu +PTCut+ SSS veto
       //if( !(muonColl.at(0).Pt()>27) ) return; //SingleMuon Trig Case
       
       int IdxZdecay1=-1, IdxZdecay2=-1;
       int Idxl1Z1=-1, Idxl2Z1=-1, Idxl1Z2=-1, Idxl2Z2=-1;
       for(int i=2; i<truthColl.size(); i++){
         int pid=truthColl.at(i).PdgId();
         int midx=truthColl.at(i).IndexMother();
         int mpid=truthColl.at(midx).PdgId();
         if(fabs(pid)>10 && fabs(pid)<20 && mpid==23 ){
           if(IdxZdecay1==-1)        IdxZdecay1=midx;
           else if(midx!=IdxZdecay1) IdxZdecay2=midx;
         }       
       }
       for(int i=2; i<truthColl.size(); i++){
         int pid=truthColl.at(i).PdgId();
         int midx=truthColl.at(i).IndexMother();
         int mpid=truthColl.at(midx).PdgId();
         if(fabs(pid)>10 && fabs(pid)<20 && mpid==23 ){
           if(midx==IdxZdecay1){
             if(Idxl1Z1==-1) Idxl1Z1=i;
             else Idxl2Z1=i;
           }
           else{
             if(Idxl1Z2==-1) Idxl1Z2=i;
             else Idxl2Z2=i;
           }
         }
       }
  
       FillHist("MZdist", truthColl.at(IdxZdecay1).M(), weight, 0., 200., 200);
       FillHist("MZdist", truthColl.at(IdxZdecay2).M(), weight, 0., 200., 200);
       FillHist("MZZdist", (truthColl.at(IdxZdecay1)+truthColl.at(IdxZdecay2)).M(), weight, 0., 500., 500);
       float MZ1=truthColl.at(IdxZdecay1).M();
       float MZ2=truthColl.at(IdxZdecay2).M();
       int IdxlZlead=-1, IdxlZsubl=-1, IdxlGlead=-1, IdxlGsubl=-1;
       //if( true ){//NoCut
       //if( fabs(MZ1-91.2)<15 && MZ2<5 ){
       //if( fabs(MZ1-91.2)<15 && MZ2<6 ){
       //if( fabs(MZ1-91.2)<15 && MZ2<10 ){
       //if( fabs(MZ1-91.2)<15 && MZ2<30 ){
       //if( fabs(MZ1-91.2)<15 && fabs(MZ2-91.2)<15 ){
       //if( MZ1<50 && MZ2<50 ){
       //if( MZ1>30 && MZ2<6 ){
       if( MZ1>4 && MZ2<6 ){
         IdxlZlead = truthColl.at(Idxl1Z1).Pt()>truthColl.at(Idxl2Z1).Pt()? Idxl1Z1:Idxl2Z1;
         IdxlZsubl = truthColl.at(Idxl1Z1).Pt()>truthColl.at(Idxl2Z1).Pt()? Idxl2Z1:Idxl1Z1;
         IdxlGlead = truthColl.at(Idxl1Z2).Pt()>truthColl.at(Idxl2Z2).Pt()? Idxl1Z2:Idxl2Z2;
         IdxlGsubl = truthColl.at(Idxl1Z2).Pt()>truthColl.at(Idxl2Z2).Pt()? Idxl2Z2:Idxl1Z2;
       }
       //else if( fabs(MZ2-91.2)<15 && MZ1<10 ){
       //else if( fabs(MZ2-91.2)<15 && MZ1<30 ){
       else if( false ){
       //else if( MZ1<7 && MZ2>30 ){
         IdxlZlead = truthColl.at(Idxl1Z2).Pt()>truthColl.at(Idxl2Z2).Pt()? Idxl1Z2:Idxl2Z2;
         IdxlZsubl = truthColl.at(Idxl1Z2).Pt()>truthColl.at(Idxl2Z2).Pt()? Idxl2Z2:Idxl1Z2;
         IdxlGlead = truthColl.at(Idxl1Z1).Pt()>truthColl.at(Idxl2Z1).Pt()? Idxl1Z1:Idxl2Z1;
         IdxlGsubl = truthColl.at(Idxl1Z1).Pt()>truthColl.at(Idxl2Z1).Pt()? Idxl2Z1:Idxl1Z1;
       }
  
  
       int Naccept=0;
       bool AcclZlead=false, AcclZsubl=false, AcclGlead=false, AcclGsubl=false;
       if(IdxlZlead!=-1 && IdxlZsubl!=-1 && IdxlGlead!=-1 && IdxlGsubl!=-1){
         if( truthColl.at(IdxlZlead).Pt()>10 && fabs(truthColl.at(IdxlZlead).Eta())<2.4 ) {Naccept++; AcclZlead=true;}
         if( truthColl.at(IdxlZsubl).Pt()>10 && fabs(truthColl.at(IdxlZsubl).Eta())<2.4 ) {Naccept++; AcclZsubl=true;}
         if( truthColl.at(IdxlGlead).Pt()>10 && fabs(truthColl.at(IdxlGlead).Eta())<2.4 ) {Naccept++; AcclGlead=true;}
         if( truthColl.at(IdxlGsubl).Pt()>10 && fabs(truthColl.at(IdxlGsubl).Eta())<2.4 ) {Naccept++; AcclGsubl=true;}
         if( Naccept!= 3) return;
         if( !(truthColl.at(IdxlZlead).Pt()>20 || truthColl.at(IdxlGlead).Pt()>20) ) return;

         FillHist("Count_NAcc3l", 0., MCweight, 0., 1., 1);
  
         if( !(AcclZlead && AcclZsubl ) ){
           if(AcclZlead){
             FillHist("M3l_ZF", (truthColl.at(IdxlZlead)+truthColl.at(IdxlGlead)+truthColl.at(IdxlGsubl)).M(), weight, 0., 500., 500);
           }
           else{
             FillHist("M3l_ZF", (truthColl.at(IdxlZsubl)+truthColl.at(IdxlGlead)+truthColl.at(IdxlGsubl)).M(), weight, 0., 500., 500);
           }
         }
         else if(AcclGlead){
           FillHist("M3l_ZT", (truthColl.at(IdxlZlead)+truthColl.at(IdxlZsubl)+truthColl.at(IdxlGlead)).M(), weight, 0., 500., 500);
           FillHist("AbsdPhi_ZmuG_ZT", fabs(truthColl.at(IdxlGlead).DeltaPhi(truthColl.at(IdxlZlead)+truthColl.at(IdxlZsubl))), weight, 0., 3.15, 200);
           FillHist("dR_ZmuG_ZT", truthColl.at(IdxlGlead).DeltaR(truthColl.at(IdxlZlead)+truthColl.at(IdxlZsubl)), weight, 0., 5., 500);
           FillHist("PTmuG_ZT", truthColl.at(IdxlGlead).Pt(), weight, 0., 100., 100);
         }
         else if(AcclGsubl){
           FillHist("M3l_ZT", (truthColl.at(IdxlZlead)+truthColl.at(IdxlZsubl)+truthColl.at(IdxlGsubl)).M(), weight, 0., 500., 500);
           FillHist("AbsdPhi_ZmuG_ZT", fabs(truthColl.at(IdxlGsubl).DeltaPhi(truthColl.at(IdxlZlead)+truthColl.at(IdxlZsubl))), weight, 0., 3.15, 200);
           FillHist("dR_ZmuG_ZT", truthColl.at(IdxlGsubl).DeltaR(truthColl.at(IdxlZlead)+truthColl.at(IdxlZsubl)), weight, 0., 5., 500);
           FillHist("PTmuG_ZT", truthColl.at(IdxlGsubl).Pt(), weight, 0., 100., 100);
         }
         
       }
     }

//     PrintTruth();     
//////////////////////////////////////////////////////////////////////////////////////////////////

     //ZG conversion test
     else if(k_sample_name.Contains("ZGto2LG") || k_sample_name.Contains("WGtoLNuG") ){

       int NZGee=0, NZGmumu=0;
       int NZGhardee=0, NZGsoftee=0, NZGhardmumu=0, NZGsoftmumu=0;
       int Idxmu1G=-1, Idxmu2G=-1, Idxe1G=-1, Idxe2G=-1;
       int Idxl1Else=-1, Idxl2Else=-1;
       int IdxGhard=-1;
       for(int i=2; i<truthColl.size(); i++){
         int pid=truthColl.at(i).PdgId();
         int midx=truthColl.at(i).IndexMother();
         int mpid=truthColl.at(midx).PdgId();
         //Decay counting
         if(pid==11 && mpid==22 ){
           NZGee++;
           int finalmpid=mpid, finalmidx=midx, finalSt=-1;
           while(finalmpid==22){
             finalSt   = truthColl.at(finalmidx).GenStatus();
             finalmidx = truthColl.at(finalmidx).IndexMother();
             finalmpid = truthColl.at(finalmidx).PdgId();
           }
           if(finalSt==23) NZGhardee++;
           else NZGsoftee++;

           FillHist("Gdecee_St", finalSt, 1., 0., 100., 100);
          
         }
         else if(pid==13 && mpid==22 ){
           NZGmumu++;
           int finalmpid=mpid, finalmidx=midx, finalSt=-1;
           while(finalmpid==22){
             finalSt   = truthColl.at(finalmidx).GenStatus();
             finalmidx = truthColl.at(finalmidx).IndexMother();
             finalmpid = truthColl.at(finalmidx).PdgId();
           }
           if(finalSt==23) NZGhardmumu++;
           else NZGsoftmumu++;
           FillHist("Gdecmumu_St", finalSt, 1., 0., 100., 100);
         }

         //For kinematic dist of gen leptons from gamma
         if( fabs(pid)==13 && truthColl.at(i).GenStatus()==1 ){

           int finalmpid=mpid, finalmidx=midx;
           while(fabs(finalmpid)==13){
             finalmidx = truthColl.at(finalmidx).IndexMother();
             finalmpid = truthColl.at(finalmidx).PdgId();
           }
           if(finalmpid==22){ 
             if(pid>0) Idxmu1G=i;
             if(pid<0) Idxmu2G=i;
           }
         }
         if( fabs(pid)==11 && truthColl.at(i).GenStatus()==1 ){

           int finalmpid=mpid, finalmidx=midx;
           while(fabs(finalmpid)==11){
             finalmidx = truthColl.at(finalmidx).IndexMother();
             finalmpid = truthColl.at(finalmidx).PdgId();
           }
           if(finalmpid==22){ 
             if(pid>0) Idxe1G=i;
             if(pid<0) Idxe2G=i;
           }
         }
         if( fabs(pid)>10 && fabs(pid)<20 && (mpid==23 || fabs(mpid)==24) ){

           if(pid>0) Idxl1Else=i;
           if(pid<0) Idxl2Else=i;
         }
         if( pid==22 && (truthColl.at(i).GenStatus()==23 || truthColl.at(i).GenStatus()==1) ) IdxGhard=i;

       }//Endof TruthFor loops
       FillHist("ZG4lComposition", 0., NZGee, 0., 10., 10); //NZG
       FillHist("ZG4lComposition", 1., NZGhardee, 0., 10., 10); //NZG
       FillHist("ZG4lComposition", 2., NZGsoftee, 0., 10., 10); //NZG
       FillHist("ZG4lComposition", 3., NZGmumu, 0., 10., 10); //NZG
       FillHist("ZG4lComposition", 4., NZGhardmumu, 0., 10., 10); //NZG
       FillHist("ZG4lComposition", 5., NZGsoftmumu, 0., 10., 10); //NZG

       if(IdxGhard!=-1) FillHist("MGhard", truthColl.at(IdxGhard).M(), weight, 0., 100., 100);

       if(Idxmu1G!=-1 && Idxmu2G!=-1 && Idxl1Else!=-1 && Idxl2Else!=-1){
         int IdxmuGlead=truthColl.at(Idxmu1G).Pt()>truthColl.at(Idxmu2G).Pt()? Idxmu1G:Idxmu2G;
         int IdxmuGsubl=truthColl.at(Idxmu1G).Pt()>truthColl.at(Idxmu2G).Pt()? Idxmu2G:Idxmu1G;
         int IdxlElselead=truthColl.at(Idxl1Else).Pt()>truthColl.at(Idxl2Else).Pt()? Idxl1Else:Idxl2Else;
         int IdxlElsesubl=truthColl.at(Idxl1Else).Pt()>truthColl.at(Idxl2Else).Pt()? Idxl2Else:Idxl1Else;

         bool AcclElselead=false, AcclElsesubl=false, AccmuGlead=false, AccmuGsubl=false;
         int NAccept=0;
         FillHist("TestM2lG", (truthColl.at(Idxmu1G)+truthColl.at(Idxmu2G)).M(), weight, 0., 10., 200);
         FillHist("TestM2lZ", (truthColl.at(Idxl1Else)+truthColl.at(Idxl2Else)).M(), weight, 0., 200., 200);
         FillHist("PTmuGlead", truthColl.at(IdxmuGlead).Pt(), weight, 0., 200., 200);
         FillHist("PTmuGsubl", truthColl.at(IdxmuGsubl).Pt(), weight, 0., 200., 200);
         if( truthColl.at(IdxmuGlead).Pt()>10 && fabs(truthColl.at(IdxmuGlead).Eta())<2.4 ) {AccmuGlead=true; NAccept++;}
         if( truthColl.at(IdxmuGsubl).Pt()>10 && fabs(truthColl.at(IdxmuGsubl).Eta())<2.4 ) {AccmuGsubl=true; NAccept++;}
         if( truthColl.at(IdxlElselead).Pt()>10 && fabs(truthColl.at(IdxlElselead).Eta())<2.4 ) {AcclElselead=true; NAccept++;}
         if( truthColl.at(IdxlElsesubl).Pt()>10 && fabs(truthColl.at(IdxlElsesubl).Eta())<2.4 ) {AcclElsesubl=true; NAccept++;}

         if( NAccept!=3 ) return;
         if( !(truthColl.at(IdxmuGlead).Pt()>20 || truthColl.at(IdxlElselead).Pt()>20) ) return;
         if( IdxGhard!=-1 ){ if( truthColl.at(IdxGhard).M()>4 ) return;}


         FillHist("dRmumuG", truthColl.at(IdxmuGlead).DeltaR(truthColl.at(IdxmuGsubl)), weight, 0., 5., 500);
         FillHist("dRmuleadG", truthColl.at(IdxmuGlead).DeltaR(truthColl.at(IdxmuGlead)+truthColl.at(IdxmuGsubl)), weight, 0., 5., 500);
         FillHist("dRmusublG", truthColl.at(IdxmuGsubl).DeltaR(truthColl.at(IdxmuGsubl)+truthColl.at(IdxmuGlead)), weight, 0., 5., 500);
         if((truthColl.at(IdxmuGlead)+truthColl.at(IdxmuGsubl)).Pt()>0){
           FillHist("PTrelmuleadG", truthColl.at(IdxmuGlead).Pt()/(truthColl.at(IdxmuGlead)+truthColl.at(IdxmuGsubl)).Pt(), weight, 0., 5., 500);
         }

         if(AccmuGlead) FillHist("AcceptComp", 0., weight, 0., 10., 10);
         else if(AccmuGsubl) FillHist("AcceptComp", 1., weight, 0., 10., 10);

       }
     }//ZG, WG sample if ends
   }//ZGtest if ends
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
  


void Mar2017_3l4j_ZGto4lTest::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Mar2017_3l4j_ZGto4lTest::BeginCycle() throw( LQError ){
  
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

Mar2017_3l4j_ZGto4lTest::~Mar2017_3l4j_ZGto4lTest() {
  
  Message("In Mar2017_3l4j_ZGto4lTest Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Mar2017_3l4j_ZGto4lTest::FillCutFlow(TString cut, float weight){
  
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



void Mar2017_3l4j_ZGto4lTest::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Mar2017_3l4j_ZGto4lTest::MakeHistograms(){
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


  AnalyzerCore::MakeHistograms("ZG4lComposition", 10, 0., 10.);
    GetHist("ZG4lComposition")->GetXaxis()->SetBinLabel(1,"G,ee");
    GetHist("ZG4lComposition")->GetXaxis()->SetBinLabel(2,"HardG,ee");
    GetHist("ZG4lComposition")->GetXaxis()->SetBinLabel(3,"SoftG,ee");
    GetHist("ZG4lComposition")->GetXaxis()->SetBinLabel(4,"G,mumu");
    GetHist("ZG4lComposition")->GetXaxis()->SetBinLabel(5,"HardG,mumu");
    GetHist("ZG4lComposition")->GetXaxis()->SetBinLabel(6,"SoftG,mumu");


  Message("Made histograms", INFO);

  // **
  // *  Remove//Overide this Mar2017_3l4j_ZGto4lTestCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Mar2017_3l4j_ZGto4lTest::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
