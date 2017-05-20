// $Id: May2017_FakeRateProd.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQMay2017_FakeRateProd Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "May2017_FakeRateProd.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (May2017_FakeRateProd);

 May2017_FakeRateProd::May2017_FakeRateProd() : AnalyzerCore(), out_muons(0) {

   SetLogName("May2017_FakeRateProd");
   Message("In May2017_FakeRateProd constructor", INFO);
   InitialiseAnalysis();
 }


 void May2017_FakeRateProd::InitialiseAnalysis() throw( LQError ) {
   
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

void May2017_FakeRateProd::ExecuteEvents()throw( LQError ){

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


   bool EleFR=false;
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("EleFR"))   EleFR=true;
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

   //Trigger Path of Analysis
   bool Pass_Trigger=false;
   if(EleFR){
     //if( PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v") )  FillHist("TrCount", 0., weight, 0., 5., 5);//0.1M
     //if( PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v") )         FillHist("TrCount", 1., weight, 0., 5., 5);//0.2M
     //if( PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ) FillHist("TrCount", 2., weight, 0., 5., 5);//0.25M
     //if( PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v") )         FillHist("TrCount", 3., weight, 0., 5., 5);//1M
     //if( PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ) FillHist("TrCount", 4., weight, 0., 5., 5);//0.9M
     //HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v   6.992 0.1M
     //HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v          6.162 0.2M
     //HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v 14.888 0.25M
     //HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v         30.397 1M
     //HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v 58.896 0.9M
     //HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v         16.43
     //HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v 63.046

     if( PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v")
        || PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v")     ) Pass_Trigger=true;
     


     //if( PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v")
     //    || PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v") ) Pass_Trigger=true;
   }

//   if(EMuMu_analysis){
//     //if( PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v") ) Pass_Trigger=true;
//     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
//          || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
//          || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger=true;
//   }
//   else if(TriMu_analysis){
//     //if( PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v") ) Pass_Trigger=true;
//     //if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ) Pass_Trigger=true;
//     if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")
//          || PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") ) Pass_Trigger=true;
//   }
//   else if(DiMuon_analysis){
//     //if( PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v") ) Pass_Trigger=true;
//     if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ) Pass_Trigger=true;
//   }
//   else if(DiEle_analysis){
//     if( PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger=true;
//   }

   //Trigger Cut
   if(!Pass_Trigger) return;
   FillCutFlow("TriggerCut", weight*pileup_reweight);
   /**********************************************************************************/

   //METFilterCut : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
   if(!PassMETFilter()) return;  FillCutFlow("EventCut", weight*pileup_reweight);

   //Vertex Cut
   //(vtx.ndof>4&&maxAbsZ<=0)||abs(vtx.z)<= 24)&&((maxd0 <=0)||abs(vtx.position.rho)<=2)&&!(vtx.isFake))
   if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; FillCutFlow("VertexCut", weight*pileup_reweight);
   
   bool EleIDSF=false, EMuTrigSF=false, EleLeg=false, MuonLeg=false, DZeff=false;

///////Objects in Analysis/////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

   //Primary Object Collection
   //std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
  
   /**PreSelCut***********************************************************************************************/
   //Intended for Code speed boosting up.
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(5.);                    eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, false);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
     if(EleFR)    { if( !(electronPreColl.size()>=1) ) return; }
     FillCutFlow("PreSel", weight);
   /**********************************************************************************************************/

     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(5.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.25);
   std::vector<snu::KMuon> muonPOGLColl; eventbase->GetMuonSel()->Selection(muonPOGLColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(5.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.6);
   std::vector<snu::KMuon> muonVetoColl; eventbase->GetMuonSel()->Selection(muonVetoColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(5.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.05);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(5.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.05);               eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");  eventbase->GetMuonSel()->SetRelIso(0.1);
   std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,false);
   std::vector<snu::KMuon> muonColl;      if(k_running_nonprompt){ muonColl=muonLooseColl;} else{ muonColl=muonTightColl;}


     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.5, 0.5);//2016 80X tuned WP
     //eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0994, 0.107);//2016 80X tuned WP
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
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.5, 0.5);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.1);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetdxySigMax(3.);
   std::vector<snu::KElectron> electronFakeLColl; eventbase->GetElectronSel()->Selection(electronFakeLColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.1);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetdxySigMax(3.);
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

/*   if(EventCand){
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
   }*/

   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight*fake_weight*btag_sf*geneff_weight*gennorm_weight;
   /***************************************************************************************************/

   //////Basic Objects Check//////////////////////////////////////////////////////////
   FillHist("Basic_Nj_orig", njets, weight, 0., 10., 10);
   FillHist("Basic_Nb_orig", nbjets, weight, 0., 10., 10);//Bjet def CVSincV2>0.89;pfKJet::CSVv2IVF_discriminator 
   FillHist("Basic_NmuT_orig", muonColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NmuL_orig", muonLooseColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeT_orig", electronColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeL_orig", electronLooseColl.size(), weight, 0., 10., 10);
   FillHist("Basic_METdist_orig", met, weight, 0., 200., 100);
   ///////////////////////////////////////////////////////////////////////////////////


/************************************************************************************/
//////////////////////////////////////////////////////////////////////////////////////
/////// MAIN ANALYSIS CODE ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
/************************************************************************************/


   if(EleFR){
     //Normalisation Check
     if( electronColl.size()==2){
       if( electronColl.at(0).Charge()!=electronColl.at(1).Charge() ){
         if( fabs((electronColl.at(0)+electronColl.at(1)).M()-91.2)<15 ){
           if( PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v") ){
             FillHist("Mee_Ele17", (electronColl.at(0)+electronColl.at(1)).M(), weight, 60., 120., 60);
           }
           if( PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ){
             FillHist("Mee_Ele17J30", (electronColl.at(0)+electronColl.at(1)).M(), weight, 60., 120., 60);
           }
         }
       }
     } 

     if( electronVetoColl.size()!=1 ) return;
     if( electronFakeLColl.size()!=1 ) return;
     if( electronFakeLColl.at(0).Pt()<20 ) return;
       FillCutFlow("1e", weight);
     if( jetColl.size()==0 ) return;
       FillCutFlow("#geq1j", weight);
     if( jetColl.size()!=1 ) return;
       FillCutFlow("1j", weight);

     float MTW = sqrt(2)*sqrt(met*electronFakeLColl.at(0).Pt()-met_x*electronFakeLColl.at(0).Px()-met_y*electronFakeLColl.at(0).Py());
     FillHist("MTW", MTW, weight, 0., 200., 200);
     if(electronColl.size()==1){
       if(met>30 && MTW>70 && MTW<120){
         if( PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v") ){
           FillHist("MTW_Tight_Ele17", MTW, weight, 0., 200., 200);
         }
         if( PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ){
           FillHist("MTW_Tight_Ele17J30", MTW, weight, 0., 200., 200);
         }

       }
     }
 

     if( met>25 ) return;
       FillCutFlow("met20", weight);
 
     if( MTW>25 ) return;
     

     float EleXbinEdges[6]={25.,30.,40.,50.,70.,100.};
     float EleYbinEdges[5]={0., 1.0, 1.5, 2.0, 2.5};
     int EleNbinsX=5, EleNbinsY=4;

     if( PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v") ){
    
       if(jetColl.at(0).DeltaR(electronFakeLColl.at(0))>1.0){
         if(jetColl.at(0).Pt()>30){
           FillHist("FR_Ne_POGTIPIsop5Ele17_Pt30_Pt", electronFakeLColl.at(0).Pt(), weight, EleXbinEdges, EleNbinsX);
           FillHist("FR_Ne_POGTIPIsop5Ele17_Pt30_PtEta", electronFakeLColl.at(0).Pt(), fabs(electronFakeLColl.at(0).Eta()), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);
           if( ( fabs(electronFakeLColl.at(0).Eta())<1.447 && electronFakeLColl.at(0).PFRelIso(0.3)<0.0588 )
               || (fabs(electronFakeLColl.at(0).Eta())>1.447 && electronFakeLColl.at(0).PFRelIso(0.3)<0.0571 ) ){
             FillHist("FR_Ne_POGTIPIsop5Ele17Isop05_Pt30_Pt", electronFakeLColl.at(0).Pt(), weight, EleXbinEdges, EleNbinsX);
             FillHist("FR_Ne_POGTIPIsop5Ele17Isop05_Pt30_PtEta", electronFakeLColl.at(0).Pt(), fabs(electronFakeLColl.at(0).Eta()), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);
           }
         }
         if(jetColl.at(0).Pt()>40){
           FillHist("FR_Ne_POGTIPIsop5Ele17_Pt40_Pt", electronFakeLColl.at(0).Pt(), weight, EleXbinEdges, EleNbinsX);
           FillHist("FR_Ne_POGTIPIsop5Ele17_Pt40_PtEta", electronFakeLColl.at(0).Pt(), fabs(electronFakeLColl.at(0).Eta()), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);
           if( ( fabs(electronFakeLColl.at(0).Eta())<1.447 && electronFakeLColl.at(0).PFRelIso(0.3)<0.0588 )
               || (fabs(electronFakeLColl.at(0).Eta())>1.447 && electronFakeLColl.at(0).PFRelIso(0.3)<0.0571 ) ){
             FillHist("FR_Ne_POGTIPIsop5Ele17Isop05_Pt40_Pt", electronFakeLColl.at(0).Pt(), weight, EleXbinEdges, EleNbinsX);
             FillHist("FR_Ne_POGTIPIsop5Ele17Isop05_Pt40_PtEta", electronFakeLColl.at(0).Pt(), fabs(electronFakeLColl.at(0).Eta()), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);
           }
         }
         if(jetColl.at(0).Pt()>50){
           FillHist("FR_Ne_POGTIPIsop5Ele17_Pt50_Pt", electronFakeLColl.at(0).Pt(), weight, EleXbinEdges, EleNbinsX);
           FillHist("FR_Ne_POGTIPIsop5Ele17_Pt50_PtEta", electronFakeLColl.at(0).Pt(), fabs(electronFakeLColl.at(0).Eta()), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);
           if( ( fabs(electronFakeLColl.at(0).Eta())<1.447 && electronFakeLColl.at(0).PFRelIso(0.3)<0.0588 )
               || (fabs(electronFakeLColl.at(0).Eta())>1.447 && electronFakeLColl.at(0).PFRelIso(0.3)<0.0571 ) ){
             FillHist("FR_Ne_POGTIPIsop5Ele17Isop05_Pt50_Pt", electronFakeLColl.at(0).Pt(), weight, EleXbinEdges, EleNbinsX);
             FillHist("FR_Ne_POGTIPIsop5Ele17Isop05_Pt50_PtEta", electronFakeLColl.at(0).Pt(), fabs(electronFakeLColl.at(0).Eta()), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);
           }
         }

       }//End of dR1.0
     }//End of Ele17NoJet Trig

     if( PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ){
    
       if(jetColl.at(0).DeltaR(electronFakeLColl.at(0))>1.0){
         if(jetColl.at(0).Pt()>30){
           FillHist("FR_Ne_POGTIPIsop5Ele17J30_Pt30_Pt", electronFakeLColl.at(0).Pt(), weight, EleXbinEdges, EleNbinsX);
           FillHist("FR_Ne_POGTIPIsop5Ele17J30_Pt30_PtEta", electronFakeLColl.at(0).Pt(), fabs(electronFakeLColl.at(0).Eta()), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);
           if( ( fabs(electronFakeLColl.at(0).Eta())<1.447 && electronFakeLColl.at(0).PFRelIso(0.3)<0.0588 )
               || (fabs(electronFakeLColl.at(0).Eta())>1.447 && electronFakeLColl.at(0).PFRelIso(0.3)<0.0571 ) ){
             FillHist("FR_Ne_POGTIPIsop5Ele17J30Isop05_Pt30_Pt", electronFakeLColl.at(0).Pt(), weight, EleXbinEdges, EleNbinsX);
             FillHist("FR_Ne_POGTIPIsop5Ele17J30Isop05_Pt30_PtEta", electronFakeLColl.at(0).Pt(), fabs(electronFakeLColl.at(0).Eta()), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);
           }
         }
         if(jetColl.at(0).Pt()>40){
           FillHist("FR_Ne_POGTIPIsop5Ele17J30_Pt40_Pt", electronFakeLColl.at(0).Pt(), weight, EleXbinEdges, EleNbinsX);
           FillHist("FR_Ne_POGTIPIsop5Ele17J30_Pt40_PtEta", electronFakeLColl.at(0).Pt(), fabs(electronFakeLColl.at(0).Eta()), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);
           if( ( fabs(electronFakeLColl.at(0).Eta())<1.447 && electronFakeLColl.at(0).PFRelIso(0.3)<0.0588 )
               || (fabs(electronFakeLColl.at(0).Eta())>1.447 && electronFakeLColl.at(0).PFRelIso(0.3)<0.0571 ) ){
             FillHist("FR_Ne_POGTIPIsop5Ele17J30Isop05_Pt40_Pt", electronFakeLColl.at(0).Pt(), weight, EleXbinEdges, EleNbinsX);
             FillHist("FR_Ne_POGTIPIsop5Ele17J30Isop05_Pt40_PtEta", electronFakeLColl.at(0).Pt(), fabs(electronFakeLColl.at(0).Eta()), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);

           }
         }
         if(jetColl.at(0).Pt()>50){
           FillHist("FR_Ne_POGTIPIsop5Ele17J30_Pt50_Pt", electronFakeLColl.at(0).Pt(), weight, EleXbinEdges, EleNbinsX);
           FillHist("FR_Ne_POGTIPIsop5Ele17J30_Pt50_PtEta", electronFakeLColl.at(0).Pt(), fabs(electronFakeLColl.at(0).Eta()), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);
           if( ( fabs(electronFakeLColl.at(0).Eta())<1.447 && electronFakeLColl.at(0).PFRelIso(0.3)<0.0588 )
               || (fabs(electronFakeLColl.at(0).Eta())>1.447 && electronFakeLColl.at(0).PFRelIso(0.3)<0.0571 ) ){
             FillHist("FR_Ne_POGTIPIsop5Ele17J30Isop05_Pt50_Pt", electronFakeLColl.at(0).Pt(), weight, EleXbinEdges, EleNbinsX);
             FillHist("FR_Ne_POGTIPIsop5Ele17J30Isop05_Pt50_PtEta", electronFakeLColl.at(0).Pt(), fabs(electronFakeLColl.at(0).Eta()), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);

           }
         }
       }//End of dRjl1.0

     }//End of Ele17Jet30 Trig


   }//EleFR Ends


/////////////////////////////////////////////////////////////////////////////////// 


return;
}// End of execute event loop
  


void May2017_FakeRateProd::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void May2017_FakeRateProd::BeginCycle() throw( LQError ){
  
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

May2017_FakeRateProd::~May2017_FakeRateProd() {
  
  Message("In May2017_FakeRateProd Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void May2017_FakeRateProd::FillCutFlow(TString cut, float weight){
  
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
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(5,"PreSel");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(6,"1e");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(7,"#geq1j");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(8,"1j");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(9,"met20");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(10,"ZVeto");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(11,"Nlj2");
    }
    if(!GetHist("cutflow_N")){
      AnalyzerCore::MakeHistograms("cutflow_N", 11, 0., 11.);
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(1,"NoCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(2,"TriggerCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(3,"EventCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(4,"VertexCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(5,"PreSel");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(6,"1e");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(7,"#geq1j");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(8,"1j");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(9,"met20");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(10,"ZVeto");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(11,"Nlj2");
    }
  }
}



void May2017_FakeRateProd::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void May2017_FakeRateProd::MakeHistograms(){
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
 


  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this May2017_FakeRateProdCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void May2017_FakeRateProd::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
