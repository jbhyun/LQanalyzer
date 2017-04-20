// $Id: Mar2017_ZZControl.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQMar2017_ZZControl Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Mar2017_ZZControl.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Mar2017_ZZControl);

 Mar2017_ZZControl::Mar2017_ZZControl() : AnalyzerCore(), out_muons(0) {

   SetLogName("Mar2017_ZZControl");
   Message("In Mar2017_ZZControl constructor", INFO);
   InitialiseAnalysis();
 }


 void Mar2017_ZZControl::InitialiseAnalysis() throw( LQError ) {
   
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

void Mar2017_ZZControl::ExecuteEvents()throw( LQError ){

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
     if     (TriMu_analysis) { if( !(muonPreColl.size()>=4)     ) return; }
     else if(EMuMu_analysis) { if( !(muonPreColl.size()>=2 && electronPreColl.size()>=2) ) return; }
   /**********************************************************************************************************/

     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.25);
   std::vector<snu::KMuon> muonPOGLColl; eventbase->GetMuonSel()->Selection(muonPOGLColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.05);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
   std::vector<snu::KMuon> muonFakeLColl; eventbase->GetMuonSel()->Selection(muonFakeLColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.05);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.1);
   std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,true);
   std::vector<snu::KMuon> muonColl;      if(k_running_nonprompt){ muonColl=muonFakeLColl;} else{ muonColl=muonTightColl;}
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

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     eventbase->GetJetSel()->SetPileUpJetID(true,"Loose");
     bool LeptonVeto=true;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonFakeLColl, electronVetoColl);



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
   if     (TriMu_analysis) { if(muonFakeLColl.size()==4)     EventCand=true; }
   else if(EMuMu_analysis) { if(muonFakeLColl.size()==2 && electronLooseColl.size()==1) EventCand=true; }
   else if(DiMuon_analysis){ if(muonFakeLColl.size()==2)     EventCand=true; }
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
         fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonFakeLColl, "MUON_HN_TRI_TIGHT", muonFakeLColl.size(), electronNull, "ELECTRON_HN_TIGHT", electronNull.size());
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

     //4Lepton
     //Step1 : 3mu +PTCut+ SSS veto
     //if( !(muonVetoColl.size()==3 && electronVetoColl.size()==0) ) return;//No additional leptons - keep 4th mu ? ~ 5% These guys have impact of ~25%
     if( !(muonFakeLColl.size()==4 && electronLooseColl.size()==0) ) return;//For valid usage of fakeable object method. FakeLooseSel should be exactly 3.
     if( !(muonColl.size()==4     && electronColl.size()==0) ) return;
     if( !(muonColl.at(0).Pt()>20. && muonColl.at(3).Pt()>10.) ) return;
     if( fabs(SumCharge(muonColl))!=0 ) return;
     FillHist("WeightDist", weight, 1, -10., 10., 2000);

//     if(bjetColl.size()!=0) return;

     //General variables
     std::vector<int> IdxmupColl; for(int i=0;i<muonColl.size();i++){ if(muonColl.at(i).Charge()>0) IdxmupColl.push_back(i); }
     std::vector<int> IdxmumColl; for(int i=0;i<muonColl.size();i++){ if(muonColl.at(i).Charge()<0) IdxmumColl.push_back(i); }

     float Mmumu1_1=(muonColl.at(IdxmupColl.at(0))+muonColl.at(IdxmumColl.at(0))).M();
     float Mmumu1_2=(muonColl.at(IdxmupColl.at(1))+muonColl.at(IdxmumColl.at(1))).M();
     float Mmumu2_1=(muonColl.at(IdxmupColl.at(0))+muonColl.at(IdxmumColl.at(1))).M();
     float Mmumu2_2=(muonColl.at(IdxmupColl.at(1))+muonColl.at(IdxmumColl.at(0))).M();
     float dM1_1=fabs(Mmumu1_1-91.2);
     float dM1_2=fabs(Mmumu1_2-91.2);
     float dM2_1=fabs(Mmumu2_1-91.2);
     float dM2_2=fabs(Mmumu2_2-91.2);
     float M4mu=(muonColl.at(0)+muonColl.at(1)+muonColl.at(2)+muonColl.at(3)).M();

      
     if( Mmumu1_1>4 && Mmumu1_2>4 && Mmumu2_1>4 && Mmumu2_2>4 ){

       bool onZonZ=false, onZoffZ_FSR=false, elseZZ=false; 
       int Idxmu1_Z1=-1, Idxmu2_Z1=-1, Idxmu1_Z2=-1, Idxmu2_Z2=-1;
       if( dM1_1<10. && dM1_2<10. ){
         Idxmu1_Z1=muonColl.at(IdxmupColl.at(0)).Pt()>muonColl.at(IdxmumColl.at(0)).Pt()? IdxmupColl.at(0):IdxmumColl.at(0);
         Idxmu2_Z1=muonColl.at(IdxmupColl.at(0)).Pt()>muonColl.at(IdxmumColl.at(0)).Pt()? IdxmumColl.at(0):IdxmupColl.at(0);
         Idxmu1_Z2=muonColl.at(IdxmupColl.at(1)).Pt()>muonColl.at(IdxmumColl.at(1)).Pt()? IdxmupColl.at(1):IdxmumColl.at(1);
         Idxmu2_Z2=muonColl.at(IdxmupColl.at(1)).Pt()>muonColl.at(IdxmumColl.at(1)).Pt()? IdxmumColl.at(1):IdxmupColl.at(1);
         if(muonColl.at(Idxmu1_Z1).Pt()>20. && muonColl.at(Idxmu1_Z2).Pt()>20.){
           onZonZ=true; FillHist("Category", 0., weight, 0., 10., 10); 
         }
       }
       if( dM2_1<10. && dM2_2<10. && ( (dM2_1+dM2_2)<(dM1_1+dM1_2) ) ){
         Idxmu1_Z1=muonColl.at(IdxmupColl.at(0)).Pt()>muonColl.at(IdxmumColl.at(1)).Pt()? IdxmupColl.at(0):IdxmumColl.at(1);
         Idxmu2_Z1=muonColl.at(IdxmupColl.at(0)).Pt()>muonColl.at(IdxmumColl.at(1)).Pt()? IdxmumColl.at(1):IdxmupColl.at(0);
         Idxmu1_Z2=muonColl.at(IdxmupColl.at(1)).Pt()>muonColl.at(IdxmumColl.at(0)).Pt()? IdxmupColl.at(1):IdxmumColl.at(0);
         Idxmu2_Z2=muonColl.at(IdxmupColl.at(1)).Pt()>muonColl.at(IdxmumColl.at(0)).Pt()? IdxmumColl.at(0):IdxmupColl.at(1);
         if(muonColl.at(Idxmu1_Z1).Pt()>20. && muonColl.at(Idxmu1_Z2).Pt()>20.){
           onZonZ=true; FillHist("Category", 0., weight, 0., 10., 10);
         }
       }
       if( (!onZonZ) && fabs(M4mu-91.2)<10. ){
         onZoffZ_FSR=true; FillHist("Category", 1., weight, 0., 10., 10);
       }
       if( !(onZonZ || onZoffZ_FSR) ){ elseZZ=true; FillHist("Category", 2., weight, 0., 10., 10); }
  
       if(onZonZ || onZoffZ_FSR) FillHist("Yield_ZZ", 0., weight, 0., 1., 1);
  
       FillHist("Mmumu1_1", Mmumu1_1, weight, 0., 200., 200);
       FillHist("Mmumu1_2", Mmumu1_2, weight, 0., 200., 200);
       FillHist("Mmumu2_1", Mmumu2_1, weight, 0., 200., 200);
       FillHist("Mmumu2_2", Mmumu2_2, weight, 0., 200., 200);
       FillHist("M4mu", M4mu, weight, 0., 700., 700);

       if(onZonZ){
         FillHist("MZ1_onZonZ", (muonColl.at(Idxmu1_Z1)+muonColl.at(Idxmu2_Z1)).M(), weight, 0., 200., 200);
         FillHist("MZ2_onZonZ", (muonColl.at(Idxmu1_Z2)+muonColl.at(Idxmu2_Z2)).M(), weight, 0., 200., 200);
         FillHist("M4mu_onZonZ", M4mu, weight, 0., 700., 700);
       }
     }//End of ZZ test

     //ZG test
     if( fabs(M4mu-91.2)<10 ){
       //if( !(fabs(Mmumu1_1-3.1)>0.1 && fabs(Mmumu1_2-3.1)>0.1 && fabs(Mmumu2_1-3.1)>0.1 && fabs(Mmumu2_2-3.1)>0.1) ) return;
       int Idxmu1_Z=-1, Idxmu2_Z=-1, Idxmu1_G=-1, Idxmu2_G=-1;
       if( Mmumu1_1>4 && Mmumu1_2<4 ){
         Idxmu1_Z=muonColl.at(IdxmupColl.at(0)).Pt()>muonColl.at(IdxmumColl.at(0)).Pt()? IdxmupColl.at(0):IdxmumColl.at(0);
         Idxmu2_Z=muonColl.at(IdxmupColl.at(0)).Pt()>muonColl.at(IdxmumColl.at(0)).Pt()? IdxmumColl.at(0):IdxmupColl.at(0);
         Idxmu1_G=muonColl.at(IdxmupColl.at(1)).Pt()>muonColl.at(IdxmumColl.at(1)).Pt()? IdxmupColl.at(1):IdxmumColl.at(1);
         Idxmu2_G=muonColl.at(IdxmupColl.at(1)).Pt()>muonColl.at(IdxmumColl.at(1)).Pt()? IdxmumColl.at(1):IdxmupColl.at(1);
       }
       else if( Mmumu1_2>4 && Mmumu1_1<4 ){
         Idxmu1_Z=muonColl.at(IdxmupColl.at(1)).Pt()>muonColl.at(IdxmumColl.at(1)).Pt()? IdxmupColl.at(1):IdxmumColl.at(1);
         Idxmu2_Z=muonColl.at(IdxmupColl.at(1)).Pt()>muonColl.at(IdxmumColl.at(1)).Pt()? IdxmumColl.at(1):IdxmupColl.at(1);
         Idxmu1_G=muonColl.at(IdxmupColl.at(0)).Pt()>muonColl.at(IdxmumColl.at(0)).Pt()? IdxmupColl.at(0):IdxmumColl.at(0);
         Idxmu2_G=muonColl.at(IdxmupColl.at(0)).Pt()>muonColl.at(IdxmumColl.at(0)).Pt()? IdxmumColl.at(0):IdxmupColl.at(0);
       }
       else if( Mmumu2_1>4 && Mmumu2_2<4 ){
         Idxmu1_Z=muonColl.at(IdxmupColl.at(0)).Pt()>muonColl.at(IdxmumColl.at(1)).Pt()? IdxmupColl.at(0):IdxmumColl.at(1);
         Idxmu2_Z=muonColl.at(IdxmupColl.at(0)).Pt()>muonColl.at(IdxmumColl.at(1)).Pt()? IdxmumColl.at(1):IdxmupColl.at(0);
         Idxmu1_G=muonColl.at(IdxmupColl.at(1)).Pt()>muonColl.at(IdxmumColl.at(0)).Pt()? IdxmupColl.at(1):IdxmumColl.at(0);
         Idxmu2_G=muonColl.at(IdxmupColl.at(1)).Pt()>muonColl.at(IdxmumColl.at(0)).Pt()? IdxmumColl.at(0):IdxmupColl.at(1);
       }
       else if( Mmumu2_2>4 && Mmumu2_1<4 ){
         Idxmu1_G=muonColl.at(IdxmupColl.at(0)).Pt()>muonColl.at(IdxmumColl.at(1)).Pt()? IdxmupColl.at(0):IdxmumColl.at(1);
         Idxmu2_G=muonColl.at(IdxmupColl.at(0)).Pt()>muonColl.at(IdxmumColl.at(1)).Pt()? IdxmumColl.at(1):IdxmupColl.at(0);
         Idxmu1_Z=muonColl.at(IdxmupColl.at(1)).Pt()>muonColl.at(IdxmumColl.at(0)).Pt()? IdxmupColl.at(1):IdxmumColl.at(0);
         Idxmu2_Z=muonColl.at(IdxmupColl.at(1)).Pt()>muonColl.at(IdxmumColl.at(0)).Pt()? IdxmumColl.at(0):IdxmupColl.at(1);
       }
       else return;

       FillHist("M4mu_ZG", M4mu, weight, 0., 200., 200);
       FillHist("MZ_ZG", (muonColl.at(Idxmu1_Z)+muonColl.at(Idxmu2_Z)).M(), weight, 0., 300., 300);
       FillHist("MG_ZG", (muonColl.at(Idxmu1_G)+muonColl.at(Idxmu2_G)).M(), weight, 0., 10., 100);
     }//End of ZG Test

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
  


void Mar2017_ZZControl::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Mar2017_ZZControl::BeginCycle() throw( LQError ){
  
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

Mar2017_ZZControl::~Mar2017_ZZControl() {
  
  Message("In Mar2017_ZZControl Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Mar2017_ZZControl::FillCutFlow(TString cut, float weight){
  
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



void Mar2017_ZZControl::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Mar2017_ZZControl::MakeHistograms(){
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
  // *  Remove//Overide this Mar2017_ZZControlCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Mar2017_ZZControl::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
