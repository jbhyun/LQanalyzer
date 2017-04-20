// $Id: Apr2017_SignalRegion.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQApr2017_SignalRegion Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Apr2017_SignalRegion.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Apr2017_SignalRegion);

 Apr2017_SignalRegion::Apr2017_SignalRegion() : AnalyzerCore(), out_muons(0) {

   SetLogName("Apr2017_SignalRegion");
   Message("In Apr2017_SignalRegion constructor", INFO);
   InitialiseAnalysis();
 }


 void Apr2017_SignalRegion::InitialiseAnalysis() throw( LQError ) {
   
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

void Apr2017_SignalRegion::ExecuteEvents()throw( LQError ){

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


   bool TriMu_analysis=true, EMuMu_analysis=false, DiMuon_analysis=false, DiEle_analysis=false;
   //bool TriMu_analysis=false, EMuMu_analysis=true, DiMuon_analysis=false, DiEle_analysis=false;
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
   if(EMuMu_analysis){
     //if( PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v") ) Pass_Trigger=true;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
          || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
          || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger=true;
   }
   else if(TriMu_analysis){
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
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.5, 0.5);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.2);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetCheckCharge(true);
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.2);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetCheckCharge(true);
   std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->Selection(electronColl);
   std::vector<snu::KElectron> electronNull;

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     //eventbase->GetJetSel()->SetPileUpJetID(true,"Loose");
     bool LeptonVeto=true;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonLooseColl, electronVetoColl);



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
   float geneff_weight=1., gennorm_weight=1.;

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

       geneff_weight   = GenFilterEfficiency(k_sample_name);
       gennorm_weight  = SignalNorm(k_sample_name, 200.);
     }
     else{
       //Perfect Prompt Ratio Approximation applied.(p=1), Applicable to generic number, combination of leptons under premise of high prompt rate.
       if(k_running_nonprompt){
         fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", muonLooseColl.size(), electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", electronLooseColl.size());
       }
     }
   }

   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight*fake_weight*btag_sf*geneff_weight*gennorm_weight;
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
     if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>15) ) return;
     if( SumCharge(muonColl)!=0 ) return;
     FillCutFlow("3lCut", weight);
     FillHist("WeightDist", weight, 1, -10., 10., 2000);

     //Variables
     int IdxZCandLead = GetDaughterCandIdx(muonColl, "Z", 10., "Lead");
     int IdxZCandSubl = GetDaughterCandIdx(muonColl, "Z", 10., "Subl");
     float Pzv1=GetvPz(v,electronColl.at(0),1), Pzv2=GetvPz(v,electronColl.at(0),2);
     float Pzv_absS=fabs(Pzv1)<fabs(Pzv2) ? Pzv1 : Pzv2;
     snu::KParticle v_absS;
     v_absS.SetXYZM(met_x,met_y,Pzv_absS,0);

     float Mmumu=(muonColl.at(0)+muonColl.at(1)).M();

     float MTW = sqrt(2)*sqrt(met*electronColl.at(0).Pt()-met_x*electronColl.at(0).Px()-met_y*electronColl.at(0).Py());

     //Dimuon 4Cut - MC Valid Region
     if( Mmumu<4.) return;
     FillCutFlow("M(#mu#mu)>4", weight);

     FillHist("Nb_3lM4Cut_1e2mu", bjetColl.size(), weight, 0., 10., 10);
     FillHist("Nlj_3lM4Cut_1e2mu", ljetColl.size(), weight, 0., 10., 10);
     FillHist("Nj_3lM4Cut_1e2mu", jetColl.size(), weight, 0., 10., 10);
     FillHist("MET_3lM4Cut_1e2mu", met, weight, 0., 200., 200);
     FillHist("MTW_3lM4Cut_1e2mu", MTW, weight, 0., 300., 300); 
 
     if(njets>=1) FillHist("Nb_3lM41jCut_1e2mu", nbjets, weight, 0., 10., 10);
     if(njets>=2) FillHist("Nb_3lM42jCut_1e2mu", nbjets, weight, 0., 10., 10);
     if(njets>=3) FillHist("Nb_3lM42jCut_1e2mu", nbjets, weight, 0., 10., 10);


     //1b Cut (Implicit 1jet cut)
     if(nbjets==0) return;
     FillCutFlow("#geq1b", weight);

     FillHist("Nb_3lM41bCut_1e2mu", bjetColl.size(), weight, 0., 10., 10);
     FillHist("Nlj_3lM41bCut_1e2mu", ljetColl.size(), weight, 0., 10., 10);
     FillHist("Nj_3lM41bCut_1e2mu", jetColl.size(), weight, 0., 10., 10);
     FillHist("MET_3lM41bCut_1e2mu", met, weight, 0., 200., 200);
     FillHist("MTW_3lM41bCut_1e2mu", MTW, weight, 0., 300., 300); 
     FillHist("Mmumu_3lM41bCut_1e2mu", Mmumu, weight, 0., 300., 300);


     //2j Cut
     if(njets<2) return;
     FillCutFlow("2j", weight);

     FillHist("Mmumu_3lM41b2jCut_1e2mu", Mmumu, weight, 0., 300., 300);

     //Signal Region
     //3j Cut
     if(njets<3) return;
     FillCutFlow("3j", weight);

     int IdxjWCandLead_NoLim=GetDaughterCandIdx(jetColl, "W", -1., "Lead");
     int IdxjWCandSubl_NoLim=GetDaughterCandIdx(jetColl, "W", -1., "Subl");
     int IdxljWCandLead_NoLim=-1, IdxljWCandSubl_NoLim=-1; 
     if(nljets>1){
       IdxljWCandLead_NoLim=GetDaughterCandIdx(ljetColl, "W", -1., "Lead");
       IdxljWCandSubl_NoLim=GetDaughterCandIdx(ljetColl, "W", -1., "Subl");
     }
     float M2l2j=(muonColl.at(0)+muonColl.at(1)+jetColl.at(IdxjWCandLead_NoLim)+jetColl.at(IdxjWCandSubl_NoLim)).M();
     float M2l2lj=(muonColl.at(0)+muonColl.at(1)+jetColl.at(IdxjWCandLead_NoLim)+jetColl.at(IdxjWCandSubl_NoLim)).M();
     float M3lv=(muonColl.at(0)+muonColl.at(1)+muonColl.at(2)+v_absS).M();


     //For Sake of expected ttZ without 40Cut
     if(IdxZCandLead!=-1 && IdxZCandSubl!=-1){
       FillHist("Mmumu_ttZ_3lM41b3jCut_1e2mu", (muonColl.at(IdxZCandLead)+muonColl.at(IdxZCandSubl)).M(), weight, 0., 200., 200);
       FillHist("Nb_ttZ_3lM41b3jCut_1e2mu", nbjets, weight, 0., 10., 10);
       FillHist("Nlj_ttZ_3lM41b3jCut_1e2mu", nljets, weight, 0., 10., 10);
       FillHist("Nj_ttZ_3lM41b3jCut_1e2mu", njets, weight, 0., 10., 10);
       FillHist("PTmu2Z_ttZ_3lM41b3jCut_1e2mu", muonColl.at(IdxZCandSubl).Pt(), weight, 0., 200., 200);
     }

     FillHist("Nb_3lM41b3jCut_1e2mu", nbjets, weight, 0., 10., 10);
     FillHist("Nlj_3lM41b3jCut_1e2mu", nljets, weight, 0., 10., 10);
     FillHist("Nj_3lM41b3jCut_1e2mu", njets, weight, 0., 10., 10);
     FillHist("MET_3lM41b3jCut_1e2mu", met, weight, 0., 200., 200);
     FillHist("MTW_3lM41b3jCut_1e2mu", MTW, weight, 0., 300., 300); 
     FillHist("PTmu1_3lM41b3jCut_1e2mu", muonColl.at(0).Pt(), weight, 0., 300., 300);
     FillHist("PTmu2_3lM41b3jCut_1e2mu", muonColl.at(1).Pt(), weight, 0., 200., 200);
     FillHist("PTmu3_3lM41b3jCut_1e2mu", muonColl.at(2).Pt(), weight, 0., 200., 200);
     FillHist("Etamu1_3lM41b3jCut_1e2mu", muonColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etamu2_3lM41b3jCut_1e2mu", muonColl.at(1).Eta(), weight, -5., 5., 100);
     FillHist("Etamu3_3lM41b3jCut_1e2mu", muonColl.at(2).Eta(), weight, -5., 5., 100);


     //Signatures
     FillHist("Mmumu_3lM41b3jCut_1e2mu", Mmumu, weight, 0., 200., 200);
     FillHist("M2l2j_3lM41b3jCut_1e2mu", M2l2j, weight, 0., 1000., 1000);
     FillHist("M3lv_3lM41b3jCut_1e2mu", M3lv, weight, 0., 1000., 1000);


     //Related to 
     //FillHist("Mmumu_M3lv_MdimuBase_3lM41b3jCut_3mu", (muonColl.at(IdxmuAOS_Mdimu)+muonColl.at(IdxmuASS_Mdimu)).M(), (muonColl.at(0)+muonColl.at(1)+muonColl.at(2)+v_absS_MdimuBase).M(), weight, 2., 202., 40, 0., 1000., 40);
     //FillHist("Mmumu_M2l2j_MdimuBase_3lM41b3jCut_3mu", (muonColl.at(IdxmuAOS_Mdimu)+muonColl.at(IdxmuASS_Mdimu)).M(), (muonColl.at(IdxmuAOS_Mdimu)+muonColl.at(IdxmuASS_Mdimu)+jetColl.at(IdxjWCandLead_NoLim)+jetColl.at(IdxjWCandSubl_NoLim)).M(), weight, 2., 202., 40, 0., 1000., 40);
     //FillHist("M3lv_M2l2j_MdimuBase_3lM41b3jCut_3mu", (muonColl.at(0)+muonColl.at(1)+muonColl.at(2)+v_absS_MdimuBase).M(), (muonColl.at(IdxmuAOS_Mdimu)+muonColl.at(IdxmuASS_Mdimu)+jetColl.at(IdxjWCandLead_NoLim)+jetColl.at(IdxjWCandSubl_NoLim)).M(), weight, 0., 1000., 40, 0., 1000., 40);
    
   

     if(IdxZCandLead!=-1 && IdxZCandSubl!=-1) return;
     FillCutFlow("ZVeto", weight);


     FillHist("Nb_3lM41b3jMZCut_3mu", nbjets, weight, 0., 10., 10);
     FillHist("Nlj_3lM41b3jMZCut_3mu", nljets, weight, 0., 10., 10);
     FillHist("Nj_3lM41b3jMZCut_3mu", njets, weight, 0., 10., 10);
     FillHist("MET_3lM41b3jMZCut_3mu", met, weight, 0., 200., 200);
     FillHist("MTW_3lM41b3jMZCut_3mu", MTW, weight, 0., 300., 300); 
     FillHist("PTmu1_3lM41b3jMZCut_3mu", muonColl.at(0).Pt(), weight, 0., 300., 300);
     FillHist("PTmu2_3lM41b3jMZCut_3mu", muonColl.at(1).Pt(), weight, 0., 200., 200);
     FillHist("PTmu3_3lM41b3jMZCut_3mu", muonColl.at(2).Pt(), weight, 0., 200., 200);
     FillHist("Etamu1_3lM41b3jMZCut_3mu", muonColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etamu2_3lM41b3jMZCut_3mu", muonColl.at(1).Eta(), weight, -5., 5., 100);
     FillHist("Etamu3_3lM41b3jMZCut_3mu", muonColl.at(2).Eta(), weight, -5., 5., 100);
 


     if(nljets<2) return;
     FillCutFlow("Nlj2", weight);






   }
   if(TriMu_analysis){

     //Step1 : 3mu +PTCut+ SSS veto
     //if( !(muonVetoColl.size()==3 && electronVetoColl.size()==0) ) return;//No additional leptons - keep 4th mu ? ~ 5% These guys have impact of ~25%
     if( !(muonLooseColl.size()==3 && electronLooseColl.size()==0) ) return;//For valid usage of fakeable object method. FakeLooseSel should be exactly 3.
     if( !(muonColl.size()==3     && electronColl.size()==0) ) return;
     if( !(muonColl.at(0).Pt()>20. && muonColl.at(2).Pt()>10.) ) return;
     if( fabs(SumCharge(muonColl))!=1 ) return;
     FillCutFlow("3lCut", weight);

     FillHist("WeightDist", weight, 1, -10., 10., 2000);

     //General variables
     int IdxZCandLead = GetDaughterCandIdx(muonColl, "Z", 10., "Lead");
     int IdxZCandSubl = GetDaughterCandIdx(muonColl, "Z", 10., "Subl");

     int IdxOS  = TriMuChargeIndex(muonColl,"OS");
     int IdxSS1 = TriMuChargeIndex(muonColl,"SS1");
     int IdxSS2 = TriMuChargeIndex(muonColl,"SS2");

     float Pzv1_1=GetvPz(v,muonColl.at(IdxSS1),1), Pzv1_2=GetvPz(v,muonColl.at(IdxSS1),2);
     float Pzv2_1=GetvPz(v,muonColl.at(IdxSS2),1), Pzv2_2=GetvPz(v,muonColl.at(IdxSS2),2);
     float Pzv_absS1=fabs(Pzv1_1)<fabs(Pzv1_2) ? Pzv1_1 : Pzv1_2;
     float Pzv_absS2=fabs(Pzv2_1)<fabs(Pzv2_2) ? Pzv2_1 : Pzv2_2;
     snu::KParticle v_absS1, v_absS2;
     v_absS1.SetXYZM(met_x,met_y,Pzv_absS1,0);
     v_absS2.SetXYZM(met_x,met_y,Pzv_absS2,0);


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
     float MTW = sqrt(2)*sqrt(met*muonColl.at(IdxWCand_NoLim).Pt()-met_x*muonColl.at(IdxWCand_NoLim).Px()-met_y*muonColl.at(IdxWCand_NoLim).Py());


     //Dimuon 4Cut - MC Valid Region
     if( !(Mmumu_OSSS1>4. && Mmumu_OSSS2>4.) ) return;
     FillCutFlow("M(#mu#mu)>4", weight);

     FillHist("Nb_3lM4Cut_3mu", bjetColl.size(), weight, 0., 10., 10);
     FillHist("Nlj_3lM4Cut_3mu", ljetColl.size(), weight, 0., 10., 10);
     FillHist("Nj_3lM4Cut_3mu", jetColl.size(), weight, 0., 10., 10);
     FillHist("MET_3lM4Cut_3mu", met, weight, 0., 200., 200);
     FillHist("MTW_3lM4Cut_3mu", MTW, weight, 0., 300., 300); 
 
     if(njets>=1) FillHist("Nb_3lM41jCut_3mu", nbjets, weight, 0., 10., 10);
     if(njets>=2) FillHist("Nb_3lM42jCut_3mu", nbjets, weight, 0., 10., 10);
     if(njets>=3) FillHist("Nb_3lM42jCut_3mu", nbjets, weight, 0., 10., 10);


     FillHist("LeakageCount", 0., weight, 0., 10., 10);
     if( IdxZCandLead!=-1 && IdxZCandSubl!=-1 ) FillHist("LeakageCount", 1., weight, 0., 10., 10);//Zwin
     if( nbjets==0 && IdxZCandLead!=-1 && IdxZCandSubl!=-1 ) FillHist("LeakageCount", 2., weight, 0., 10., 10);//Zwin no bjet
     if( nbjets==1 && IdxZCandLead!=-1 && IdxZCandSubl!=-1 ) FillHist("LeakageCount", 3., weight, 0., 10., 10);//Zwin 1b
     if( nbjets==1 && njets<3 && IdxZCandLead!=-1 && IdxZCandSubl!=-1 ) FillHist("LeakageCount", 4., weight, 0., 10., 10);//Zwin 1b, 1,2j
     if( nbjets==1 && njets<3 && IdxZCandLead!=-1 && IdxZCandSubl!=-1 && Mmumu_OSSS1>40. && Mmumu_OSSS2>40. ){
       FillHist("LeakageCount", 5., weight, 0., 10., 10);//Zwin 1b, 1,2j, Mmumu>40
     }
     if( nbjets>0 && njets>2 && IdxZCandLead!=-1 && IdxZCandSubl!=-1 ) FillHist("LeakageCount", 6., weight, 0., 10., 10);//SigReg ttZ
     if( nbjets>0 && njets>2 && Mmumu_OSSS1>40. && Mmumu_OSSS2>40. ) FillHist("LeakageCount", 7., weight, 0., 10., 10);//SigReg M>40
     if( nbjets>0 && njets>2 && Mmumu_OSSS1>40. && Mmumu_OSSS2>40. && IdxZCandLead!=-1 && IdxZCandSubl!=-1 ) FillHist("LeakageCount", 8., weight, 0., 10., 10);//SigReg ttZ, M>40
     
     FillHist("MOSSS1_MOSSS2_3mu", Mmumu_OSSS1, Mmumu_OSSS2, weight, 0., 200., 40, 0., 200., 40);


     //1b Cut (Implicit 1jet cut)
     if(nbjets==0) return;
     FillCutFlow("#geq1b", weight);

     FillHist("Nb_3lM41bCut_3mu", bjetColl.size(), weight, 0., 10., 10);
     FillHist("Nlj_3lM41bCut_3mu", ljetColl.size(), weight, 0., 10., 10);
     FillHist("Nj_3lM41bCut_3mu", jetColl.size(), weight, 0., 10., 10);
     FillHist("MET_3lM41bCut_3mu", met, weight, 0., 200., 200);
     FillHist("MTW_3lM41bCut_3mu", MTW, weight, 0., 300., 300); 
     FillHist("PTb_3lM41bCut_3mu", bjetColl.at(0).Pt(), weight, 0., 300., 300);


     //2j Cut
     if(njets<2) return;
     FillCutFlow("2j", weight);


     //Signal Region
     //3j Cut
     if(njets<3) return;
     FillCutFlow("3j", weight);


     int IdxjWCandLead_NoLim=GetDaughterCandIdx(jetColl, "W", -1., "Lead");
     int IdxjWCandSubl_NoLim=GetDaughterCandIdx(jetColl, "W", -1., "Subl");
     int IdxljWCandLead_NoLim=-1, IdxljWCandSubl_NoLim=-1; 
     if(nljets>1){
       IdxljWCandLead_NoLim=GetDaughterCandIdx(ljetColl, "W", -1., "Lead");
       IdxljWCandSubl_NoLim=GetDaughterCandIdx(ljetColl, "W", -1., "Subl");
     }
     //PT based
     float M3lv1=(muonColl.at(0)+muonColl.at(1)+muonColl.at(2)+v_absS1).M();
     float M3lv2=(muonColl.at(0)+muonColl.at(1)+muonColl.at(2)+v_absS2).M();
     float M2l2j1W=(muonColl.at(IdxOS)+muonColl.at(IdxSS1)+jetColl.at(IdxjWCandLead_NoLim)+jetColl.at(IdxjWCandSubl_NoLim)).M();
     float M2l2j2W=(muonColl.at(IdxOS)+muonColl.at(IdxSS2)+jetColl.at(IdxjWCandLead_NoLim)+jetColl.at(IdxjWCandSubl_NoLim)).M();

     //Mmumu based
     float IdxmuAOS_Mdimu=-1, IdxmuASS_Mdimu=-1, IdxmuW_Mdimu=-1;
     snu::KParticle v_absS_MdimuBase;
     if(Mmumu_OSSS1>Mmumu_OSSS2){
       IdxmuW_Mdimu=IdxSS1; IdxmuAOS_Mdimu=IdxOS; IdxmuASS_Mdimu=IdxSS2;
       v_absS_MdimuBase=v_absS1;
     }
     else{
       IdxmuW_Mdimu=IdxSS2; IdxmuAOS_Mdimu=IdxOS; IdxmuASS_Mdimu=IdxSS1;
       v_absS_MdimuBase=v_absS2;
     }

     //For Sake of expected ttZ without 40Cut
     if(IdxZCandLead!=-1 && IdxZCandSubl!=-1){
       FillHist("Mmumu_ttZ_3lM41b3jCut_3mu", (muonColl.at(IdxZCandLead)+muonColl.at(IdxZCandSubl)).M(), weight, 0., 200., 200);
       FillHist("Nb_ttZ_3lM41b3jCut_3mu", nbjets, weight, 0., 10., 10);
       FillHist("Nlj_ttZ_3lM41b3jCut_3mu", nljets, weight, 0., 10., 10);
       FillHist("Nj_ttZ_3lM41b3jCut_3mu", njets, weight, 0., 10., 10);
       FillHist("PTmu2Z_ttZ_3lM41b3jCut_3mu", muonColl.at(IdxZCandSubl).Pt(), weight, 0., 200., 200);
     }

     FillHist("Nb_3lM41b3jCut_3mu", nbjets, weight, 0., 10., 10);
     FillHist("Nlj_3lM41b3jCut_3mu", nljets, weight, 0., 10., 10);
     FillHist("Nj_3lM41b3jCut_3mu", njets, weight, 0., 10., 10);
     FillHist("MET_3lM41b3jCut_3mu", met, weight, 0., 200., 200);
     FillHist("MTW_3lM41b3jCut_3mu", MTW, weight, 0., 300., 300); 
     FillHist("PTmu1_3lM41b3jCut_3mu", muonColl.at(0).Pt(), weight, 0., 300., 300);
     FillHist("PTmu2_3lM41b3jCut_3mu", muonColl.at(1).Pt(), weight, 0., 200., 200);
     FillHist("PTmu3_3lM41b3jCut_3mu", muonColl.at(2).Pt(), weight, 0., 200., 200);
     FillHist("Etamu1_3lM41b3jCut_3mu", muonColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etamu2_3lM41b3jCut_3mu", muonColl.at(1).Eta(), weight, -5., 5., 100);
     FillHist("Etamu3_3lM41b3jCut_3mu", muonColl.at(2).Eta(), weight, -5., 5., 100);
 
     //PT based
     FillHist("Mmumu_OSSS1_3lM41b3jCut_3mu", Mmumu_OSSS1, weight, 2., 202., 200);
     FillHist("Mmumu_OSSS2_3lM41b3jCut_3mu", Mmumu_OSSS2, weight, 2., 202., 200);
     FillHist("M3lv1_3lM41b3jCut_3mu", M3lv1, weight, 0., 1000., 1000);
     FillHist("M3lv2_3lM41b3jCut_3mu", M3lv2, weight, 0., 1000., 1000);
     FillHist("M2l2j1W_3lM41b3jCut_3mu", M2l2j1W, weight, 0., 1000., 1000);
     FillHist("M2l2j2W_3lM41b3jCut_3mu", M2l2j2W, weight, 0., 1000., 1000);

     if(nljets>1){
       float M2l2lj1=(muonColl.at(IdxOS)+muonColl.at(IdxSS1)+ljetColl.at(0)+ljetColl.at(1)).M();
       float M2l2lj2=(muonColl.at(IdxOS)+muonColl.at(IdxSS2)+ljetColl.at(0)+ljetColl.at(1)).M();
       float M2l2ljW1=(muonColl.at(IdxOS)+muonColl.at(IdxSS1)+ljetColl.at(IdxljWCandLead_NoLim)+ljetColl.at(IdxljWCandSubl_NoLim)).M();
       float M2l2ljW2=(muonColl.at(IdxOS)+muonColl.at(IdxSS2)+ljetColl.at(IdxljWCandLead_NoLim)+ljetColl.at(IdxljWCandSubl_NoLim)).M();
       FillHist("M2l2ljW1_3lM41b3jCut_3mu", M2l2ljW1, weight, 0., 1000., 1000);
       FillHist("M2l2ljW2_3lM41b3jCut_3mu", M2l2ljW2, weight, 0., 1000., 1000);
     }


     //Mmumu based
     FillHist("Mmumu_OSSSL_3lM41b3jCut_3mu", max(Mmumu_OSSS1,Mmumu_OSSS2), weight, 2., 202., 200);
     FillHist("Mmumu_OSSSS_3lM41b3jCut_3mu", min(Mmumu_OSSS1,Mmumu_OSSS2), weight, 2., 202., 200);
     FillHist("M3lv_MdimuBase_3lM41b3jCut_3mu", (muonColl.at(0)+muonColl.at(1)+muonColl.at(2)+v_absS_MdimuBase).M(), weight, 0., 1000., 1000);
     FillHist("M2l2j_MdimuBase_3lM41b3jCut_3mu", (muonColl.at(IdxmuAOS_Mdimu)+muonColl.at(IdxmuASS_Mdimu)+jetColl.at(IdxjWCandLead_NoLim)+jetColl.at(IdxjWCandSubl_NoLim)).M(), weight, 0., 1000., 1000);
     FillHist("Mmumu_M3lv_MdimuBase_3lM41b3jCut_3mu", (muonColl.at(IdxmuAOS_Mdimu)+muonColl.at(IdxmuASS_Mdimu)).M(), (muonColl.at(0)+muonColl.at(1)+muonColl.at(2)+v_absS_MdimuBase).M(), weight, 2., 202., 40, 0., 1000., 40);
     FillHist("Mmumu_M2l2j_MdimuBase_3lM41b3jCut_3mu", (muonColl.at(IdxmuAOS_Mdimu)+muonColl.at(IdxmuASS_Mdimu)).M(), (muonColl.at(IdxmuAOS_Mdimu)+muonColl.at(IdxmuASS_Mdimu)+jetColl.at(IdxjWCandLead_NoLim)+jetColl.at(IdxjWCandSubl_NoLim)).M(), weight, 2., 202., 40, 0., 1000., 40);
     FillHist("M3lv_M2l2j_MdimuBase_3lM41b3jCut_3mu", (muonColl.at(0)+muonColl.at(1)+muonColl.at(2)+v_absS_MdimuBase).M(), (muonColl.at(IdxmuAOS_Mdimu)+muonColl.at(IdxmuASS_Mdimu)+jetColl.at(IdxjWCandLead_NoLim)+jetColl.at(IdxjWCandSubl_NoLim)).M(), weight, 0., 1000., 40, 0., 1000., 40);
    
   

     if(IdxZCandLead!=-1 && IdxZCandSubl!=-1) return;
     FillCutFlow("ZVeto", weight);


     FillHist("Nb_3lM41b3jMZCut_3mu", nbjets, weight, 0., 10., 10);
     FillHist("Nlj_3lM41b3jMZCut_3mu", nljets, weight, 0., 10., 10);
     FillHist("Nj_3lM41b3jMZCut_3mu", njets, weight, 0., 10., 10);
     FillHist("MET_3lM41b3jMZCut_3mu", met, weight, 0., 200., 200);
     FillHist("MTW_3lM41b3jMZCut_3mu", MTW, weight, 0., 300., 300); 
     FillHist("Mmumu_OSSS1_3lM41b3jMZCut_3mu", Mmumu_OSSS1, weight, 0., 200., 200);
     FillHist("Mmumu_OSSS2_3lM41b3jMZCut_3mu", Mmumu_OSSS2, weight, 0., 200., 200);
     FillHist("M3lv1_3lM41b3jMZCut_3mu", M3lv1, weight, 0., 1000., 1000);
     FillHist("M3lv2_3lM41b3jMZCut_3mu", M3lv2, weight, 0., 1000., 1000);
     FillHist("M2l2j1W_3lM41b3jMZCut_3mu", M2l2j1W, weight, 0., 1000., 1000);
     FillHist("M2l2j2W_3lM41b3jMZCut_3mu", M2l2j2W, weight, 0., 1000., 1000);

     FillHist("PTmu1_3lM41b3jMZCut_3mu", muonColl.at(0).Pt(), weight, 0., 300., 300);
     FillHist("PTmu2_3lM41b3jMZCut_3mu", muonColl.at(1).Pt(), weight, 0., 200., 200);
     FillHist("PTmu3_3lM41b3jMZCut_3mu", muonColl.at(2).Pt(), weight, 0., 200., 200);
     FillHist("Etamu1_3lM41b3jMZCut_3mu", muonColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etamu2_3lM41b3jMZCut_3mu", muonColl.at(1).Eta(), weight, -5., 5., 100);
     FillHist("Etamu3_3lM41b3jMZCut_3mu", muonColl.at(2).Eta(), weight, -5., 5., 100);
 


     if(nljets<2) return;
     FillCutFlow("Nlj2", weight);



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
  


void Apr2017_SignalRegion::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Apr2017_SignalRegion::BeginCycle() throw( LQError ){
  
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

Apr2017_SignalRegion::~Apr2017_SignalRegion() {
  
  Message("In Apr2017_SignalRegion Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Apr2017_SignalRegion::FillCutFlow(TString cut, float weight){
  
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



void Apr2017_SignalRegion::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Apr2017_SignalRegion::MakeHistograms(){
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
 
  AnalyzerCore::MakeHistograms("LeakageCount", 10, 0., 10.);
    GetHist("LeakageCount")->GetXaxis()->SetBinLabel(1,"Total,M>4");
    GetHist("LeakageCount")->GetXaxis()->SetBinLabel(2,"Z");
    GetHist("LeakageCount")->GetXaxis()->SetBinLabel(3,"Z,0b");
    GetHist("LeakageCount")->GetXaxis()->SetBinLabel(4,"Z,1b");
    GetHist("LeakageCount")->GetXaxis()->SetBinLabel(5,"Z,1b,1-2j");
    GetHist("LeakageCount")->GetXaxis()->SetBinLabel(6,"Z,1b,1-2j,M>40");
    GetHist("LeakageCount")->GetXaxis()->SetBinLabel(7,"SigttZ");
    GetHist("LeakageCount")->GetXaxis()->SetBinLabel(8,"SigM40");
    GetHist("LeakageCount")->GetXaxis()->SetBinLabel(9,"SigttZM40");


  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Apr2017_SignalRegionCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Apr2017_SignalRegion::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
