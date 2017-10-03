// $Id: Jul2017_EMuMuComp.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQJul2017_EMuMuComp Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Jul2017_EMuMuComp.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Jul2017_EMuMuComp);

 Jul2017_EMuMuComp::Jul2017_EMuMuComp() : AnalyzerCore(), out_muons(0) {

   SetLogName("Jul2017_EMuMuComp");
   Message("In Jul2017_EMuMuComp constructor", INFO);
   InitialiseAnalysis();
 }


 void Jul2017_EMuMuComp::InitialiseAnalysis() throw( LQError ) {
   
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

void Jul2017_EMuMuComp::ExecuteEvents()throw( LQError ){

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


   bool EleFR=false, TrigSel=false, NormCheck=false, UnPreTrig=false, PreTrig=false, Closure=false;
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("EleFR"))     EleFR     =true;
     else if(k_flags.at(i).Contains("TrigSel"))   TrigSel   =true;
     else if(k_flags.at(i).Contains("NormCheck")) NormCheck =true;
     else if(k_flags.at(i).Contains("UnPreTrig")) UnPreTrig =true;
     else if(k_flags.at(i).Contains("PreTrig"))   PreTrig   =true;
     else if(k_flags.at(i).Contains("Closure"))   Closure   =true;
   }

    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //ListTriggersAvailable();

   //Trigger Path of Analysis
   bool Pass_Trigger=false;
   float trigger_ps_weight=1., trigger_period_weight=1.;
   if(EleFR || NormCheck){
     if(UnPreTrig){
       if( PassTrigger("HLT_Ele27_WPTight_Gsf_v") ) Pass_Trigger=true;
       if(!isData) trigger_ps_weight=WeightByTrigger("HLT_Ele27_WPTight_Gsf_v", TargetLumi);
     }
     if(PreTrig){
       if( PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ) Pass_Trigger=true;
       if(!isData) trigger_ps_weight=WeightByTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v", TargetLumi);
     }
   }
   if(Closure){
     int Pass_Trigger1=0, Pass_Trigger2=0;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ) Pass_Trigger1++;
     if( PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger2++;

     if(isData){
       int DataPeriod=GetDataPeriod();
       if( DataPeriod>0 && DataPeriod<7 && Pass_Trigger1==1 ) Pass_Trigger=true;
       else if( DataPeriod==7 && Pass_Trigger2==1 ) Pass_Trigger=true;
     }
     else{
       if( Pass_Trigger1>0 || Pass_Trigger2>0 ) Pass_Trigger=true;
       trigger_period_weight=(Pass_Trigger1*27.257618+Pass_Trigger2*8.605696)/35.863314;
     }
   }
   FillHist("TriggerPSWeight", trigger_ps_weight, 1., 0., 1., 100);
   weight*=trigger_ps_weight*trigger_period_weight;


     //HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v   6.992 0.1M
     //HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v          6.162 0.2M
     //HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v 14.888 0.25M
     //HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v         30.397 1M
     //HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v 58.896 0.9M
     //HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v         16.43
     //HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v 63.046

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
   //std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
  
   /**PreSelCut***********************************************************************************************/
   //Intended for Code speed boosting up.
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, true);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
     if(EleFR || TrigSel || NormCheck)    { if( !(electronPreColl.size()>=1) ) return; }
     else if(Closure)                     { if( !(electronPreColl.size()>=1 && muonPreColl.size()>=2) ) return; }
     FillCutFlow("PreSel", weight);
   /**********************************************************************************************************/

     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetBSdz(0.1);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetBSdz(0.1);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");  eventbase->GetMuonSel()->SetRelIso(0.1);
   std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,false);
   std::vector<snu::KMuon> muonColl;      if(k_running_nonprompt){ muonColl=muonLooseColl;} else{ muonColl=muonTightColl;}


     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronFakeLPreOptColl; eventbase->GetElectronSel()->Selection(electronFakeLPreOptColl);

     std::vector<snu::KElectron> electronFakeL1Coll;
       for(int i=0; i<electronFakeLPreOptColl.size(); i++){
         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "POGMVAMFakeLIso04Opt1")) electronFakeL1Coll.push_back(electronFakeLPreOptColl.at(i));
       }
     std::vector<snu::KElectron> electronFakeL2Coll;
       for(int i=0; i<electronFakeLPreOptColl.size(); i++){
         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "POGMVAMFakeLIso04Opt2")) electronFakeL2Coll.push_back(electronFakeLPreOptColl.at(i));
       }
     std::vector<snu::KElectron> electronFakeL2Iso05Coll;
       for(int i=0; i<electronFakeLPreOptColl.size(); i++){
         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "POGMVAMFakeLIso05Opt2")) electronFakeL2Iso05Coll.push_back(electronFakeLPreOptColl.at(i));
       }

     std::vector<snu::KElectron> electronFakeVetoColl;
       for(int i=0; i<electronFakeLPreOptColl.size(); i++){
         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "POGMVAMFakeLIso04Opt2NoTrig")) electronFakeVetoColl.push_back(electronFakeLPreOptColl.at(i));
       }


     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronTightColl; eventbase->GetElectronSel()->Selection(electronTightColl);
   std::vector<snu::KElectron> electronColl;  if(!k_running_nonprompt) electronColl=electronTightColl;


     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     //eventbase->GetJetSel()->SetPileUpJetID(true,"Loose");
     bool LeptonVeto=false;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonLooseColl, electronFakeLPreOptColl);



   //Method to apply 2a method SF
   //std::vector<int> bIdxColl  = GetSFBJetIdx(jetColl,"Medium");
   //std::vector<int> ljIdxColl = GetSFLJetIdx(jetColl, bIdxColl, "Medium");
   //std::vector<snu::KJet> bjetColl2a; for(int i=0; i<bIdxColl.size(); i++) {bjetColl2a.push_back(jetColl.at(bIdxColl.at(i)));}
   //std::vector<snu::KJet> ljetColl2a; for(int i=0; i<ljIdxColl.size(); i++){ljetColl2a.push_back(jetColl.at(ljIdxColl.at(i)));}

   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");
   std::vector<snu::KJet> ljetColl = SelLightJets(jetColl, "Medium");

   double PFMETType1  = eventbase->GetEvent().PFMETType1();
   double met = eventbase->GetEvent().MET(snu::KEvent::pfmet);
   double met_x = eventbase->GetEvent().PFMETx();
   double met_y = eventbase->GetEvent().PFMETy();
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
//
   /*This part is for boosting up speed.. SF part takes rather longer time than expected*/

   if(true){
     if(!isData){
       if(NormCheck || Closure){
       //trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );
  
       id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electronColl);
       //id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronColl);
       reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
  
       id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muonColl);
       //iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
       trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
       //btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

       //geneff_weight   = GenFilterEfficiency(k_sample_name);
       //gennorm_weight  = SignalNorm(k_sample_name, 200.);
       }
     }
     else{
       //Perfect Prompt Ratio Approximation applied.(p=1), Applicable to generic number, combination of leptons under premise of high prompt rate.
   //    if(k_running_nonprompt){
   //      fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", muonLooseColl.size(), electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", electronLooseColl.size());
//       }
     }
   }

   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight*fake_weight*btag_sf*geneff_weight*gennorm_weight;
   //-----------------------------------------------------------------------------------------//


//----------------------------------------------------------------------------------//
//////////////////////////////////////////////////////////////////////////////////////
/////// MAIN ANALYSIS CODE ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//----------------------------------------------------------------------------------//

   if(TrigSel){
     if(electronColl.size()!=1) return;
     if(electronColl.at(0).Pt()<25) return;
     float MTW = sqrt(2)*sqrt(met*electronColl.at(0).Pt()-met_x*electronColl.at(0).Px()-met_y*electronColl.at(0).Py());

     if(PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v")){
       FillHist("Count_e17", 0., weight, 0., 10., 10);//1: Ele25

       if(jetColl.size()==1){
         if( jetColl.at(0).Pt()>30 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e17", 1., weight, 0., 10., 10);//2: 1j away 30
           if( met<30 && MTW<30 ) FillHist("Count_e17", 2., weight, 0., 10., 10);//FakeEnrich 1j30
         }
         if( jetColl.at(0).Pt()>40 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e17", 3., weight, 0., 10., 10);//2: 1j away 40
           if( met<30 && MTW<30 ) FillHist("Count_e17", 4., weight, 0., 10., 10);//FakeEnrich 1j40
         }
         if( jetColl.at(0).Pt()>60 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e17", 5., weight, 0., 10., 10);//2: 1j away 60
           if( met<30 && MTW<30 ) FillHist("Count_e17", 6., weight, 0., 10., 10);//FakeEnrich 1j60
         }
       }//1j 

     }//End of Ele17 
     if(PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v")){
       FillHist("Count_e17j30", 0., weight, 0., 10., 10);//1: Ele25

       if(jetColl.size()==1){
         if( jetColl.at(0).Pt()>30 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e17j30", 1., weight, 0., 10., 10);//2: 1j away 30
           if( met<30 && MTW<30 ) FillHist("Count_e17j30", 2., weight, 0., 10., 10);//FakeEnrich 1j30
         }
         if( jetColl.at(0).Pt()>40 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e17j30", 3., weight, 0., 10., 10);//2: 1j away 40
           if( met<30 && MTW<30 ) FillHist("Count_e17j30", 4., weight, 0., 10., 10);//FakeEnrich 1j40
         }
         if( jetColl.at(0).Pt()>60 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e17j30", 5., weight, 0., 10., 10);//2: 1j away 60
           if( met<30 && MTW<30 ) FillHist("Count_e17j30", 6., weight, 0., 10., 10);//FakeEnrich 1j60
         }
       }//1j 

     }//End of Ele17Jet30 
     if(PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v")){
       FillHist("Count_e23", 0., weight, 0., 10., 10);//1: Ele25

       if(jetColl.size()==1){
         if( jetColl.at(0).Pt()>30 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e23", 1., weight, 0., 10., 10);//2: 1j away 30
           if( met<30 && MTW<30 ) FillHist("Count_e23", 2., weight, 0., 10., 10);//FakeEnrich 1j30
         }
         if( jetColl.at(0).Pt()>40 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e23", 3., weight, 0., 10., 10);//2: 1j away 40
           if( met<30 && MTW<30 ) FillHist("Count_e23", 4., weight, 0., 10., 10);//FakeEnrich 1j40
         }
         if( jetColl.at(0).Pt()>60 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e23", 5., weight, 0., 10., 10);//2: 1j away 60
           if( met<30 && MTW<30 ) FillHist("Count_e23", 6., weight, 0., 10., 10);//FakeEnrich 1j60
         }
       }//1j 

     }//End of Ele23 
     if(PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v")){
       FillHist("Count_e23j30", 0., weight, 0., 10., 10);//1: Ele25

       if(jetColl.size()==1){
         if( jetColl.at(0).Pt()>30 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e23j30", 1., weight, 0., 10., 10);//2: 1j away 30
           if( met<30 && MTW<30 ) FillHist("Count_e23j30", 2., weight, 0., 10., 10);//FakeEnrich 1j30
         }
         if( jetColl.at(0).Pt()>40 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e23j30", 3., weight, 0., 10., 10);//2: 1j away 40
           if( met<30 && MTW<30 ) FillHist("Count_e23j30", 4., weight, 0., 10., 10);//FakeEnrich 1j40
         }
         if( jetColl.at(0).Pt()>60 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e23j30", 5., weight, 0., 10., 10);//2: 1j away 60
           if( met<30 && MTW<30 ) FillHist("Count_e23j30", 6., weight, 0., 10., 10);//FakeEnrich 1j60
         }
       }//1j 

     }//End of Ele23Jet30

   //170602 - Conclusion -Ele23Jet30 is the best.
   //though for single lepton events other triggers as Ele17 may give better events, but for dijet topology lepton&jet events
   //Ele23 Jet30 gives the most events. So I chose this trigger.
    
   }//End of TriggerChoice

   if(NormCheck){
     //1) High MTW region(MTW>70, MET>50) mostly includes W related events, no fakes
     //   -a) For Norm check with 2)
     //   -b) For NLO xsec check.
     //2) DY+Jets events ; very clear for prescale check
     
     std::vector<snu::KJet> jetVeto2Coll = SkimJetColl(jetColl, electronFakeL2Coll, muonLooseColl, "EleMuVeto");
     float nvtx_reweight=1., purw_down=1.;      
     if(!isData){
        nvtx_reweight=mcdata_correction->UserPileupWeight(eventbase->GetEvent(), jetVeto2Coll.size());
        if(pileup_reweight!=0) purw_down=eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::down)/pileup_reweight;
     }

     if(electronColl.size()==1 && electronFakeL2Coll.size()==1){
       if( UnPreTrig && electronColl.at(0).Pt()<30) return;
       if( PreTrig   && electronColl.at(0).Pt()<25) return;
       float MTW = sqrt(2)*sqrt(met*electronColl.at(0).Pt()-met_x*electronColl.at(0).Px()-met_y*electronColl.at(0).Py());

       //Nvtx
       if(met>30){
         FillHist("Nvtx_met30", Nvtx, weight, 0., 50., 50);
         FillHist("Nvtx_met30_PUdown", Nvtx, weight*purw_down, 0., 50., 50);
         FillHist("Nvtx_met30_NvtxRW", Nvtx, weight*nvtx_reweight, 0., 50., 50);
         if(jetVeto2Coll.size()==0) {
           FillHist("Nvtx_0jmet30", Nvtx, weight, 0., 50., 50);
           FillHist("Nvtx_0jmet30_NvtxRW", Nvtx, weight*nvtx_reweight, 0., 50., 50);
           FillHist("NvtxRW_0jmet30", Nvtx, nvtx_reweight, 0., 50., 50);
         }
         if(jetVeto2Coll.size()>0 && jetVeto2Coll.at(0).Pt()>40){
           if(jetVeto2Coll.size()==1) {
             FillHist("Nvtx_1jmet30", Nvtx, weight, 0., 50., 50);
             FillHist("Nvtx_1jmet30_NvtxRW", Nvtx, weight*nvtx_reweight, 0., 50., 50);
             FillHist("NvtxRW_1jmet30", Nvtx, nvtx_reweight, 0., 50., 50);
           }
           FillHist("Nvtx_gt1jmet30", Nvtx, weight, 0., 50., 50);
           FillHist("Nvtx_gt1jmet30_PUdown", Nvtx, weight*purw_down, 0., 50., 50);
           FillHist("Nvtx_gt1jmet30_NvtxRW", Nvtx, weight*nvtx_reweight, 0., 50., 50);
         }         
       }

       //For Validation in unprescaled path
       FillHist("PTe", electronColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("Etae", electronColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Nj", jetVeto2Coll.size(), weight, 0., 10., 10);
       if(jetVeto2Coll.size()>0){
         FillHist("PTj1", jetVeto2Coll.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etaj1", jetVeto2Coll.at(0).Eta(), weight, -5., 5., 100);
       }
       FillHist("MET", met, weight, 0., 200., 200);
       FillHist("MET_noxy", PFMETType1, weight, 0., 200., 200);
       FillHist("MTW", MTW, weight, 0., 200., 200);
       FillHist("MTW_noxy", MTW, weight, 0., 200., 200);
       FillHist("MET_NvtxRW", met, weight*nvtx_reweight, 0., 200., 200);
       FillHist("MTW_NvtxRW", MTW, weight*nvtx_reweight, 0., 200., 200);
       if(met>30){
         FillHist("PTe_met30", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_met30", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Nj_met30", jetVeto2Coll.size(), weight, 0., 10., 10);
         if(jetVeto2Coll.size()>0){
           FillHist("PTj1_met30", jetVeto2Coll.at(0).Pt(), weight, 0., 200., 200);
           FillHist("Etaj1_met30", jetVeto2Coll.at(0).Eta(), weight, -5., 5., 100);
         }
         FillHist("MET_met30", met, weight, 0., 200., 200);
         FillHist("MTW_met30", MTW, weight, 0., 200., 200);

         //NvtxRWeighted
         FillHist("PTe_met30_NvtxRW", electronColl.at(0).Pt(), weight*nvtx_reweight, 0., 200., 200);
         FillHist("Etae_met30_NvtxRW", electronColl.at(0).Eta(), weight*nvtx_reweight, -5., 5., 100);
         FillHist("Nj_met30_NvtxRW", jetVeto2Coll.size(), weight*nvtx_reweight, 0., 10., 10);
         if(jetVeto2Coll.size()>0){
           FillHist("PTj1_met30_NvtxRW", jetVeto2Coll.at(0).Pt(), weight*nvtx_reweight, 0., 200., 200);
           FillHist("Etaj1_met30_NvtxRW", jetVeto2Coll.at(0).Eta(), weight*nvtx_reweight, -5., 5., 100);
         }
         FillHist("MET_met30_NvtxRW", met, weight*nvtx_reweight, 0., 200., 200);
         FillHist("MTW_met30_NvtxRW", MTW, weight*nvtx_reweight, 0., 200., 200);
       }
       //For both prescaled & unprescaled
       if(jetVeto2Coll.size()>0 && jetVeto2Coll.at(0).Pt()>40){
         FillHist("PTe_gt1j40", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_gt1j40", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Nj_gt1j40", jetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("PTj1_gt1j40", jetVeto2Coll.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etaj1_gt1j40", jetVeto2Coll.at(0).Eta(), weight, -5., 5., 100);
         FillHist("MET_gt1j40", met, weight, 0., 200., 200);
         FillHist("MTW_gt1j40", MTW, weight, 0., 200., 200);
         if(met>30){
           FillHist("PTe_gt1j40met30", electronColl.at(0).Pt(), weight, 0., 200., 200);
           FillHist("Etae_gt1j40met30", electronColl.at(0).Eta(), weight, -5., 5., 100);
           FillHist("Nj_gt1j40met30", jetVeto2Coll.size(), weight, 0., 10., 10);
           FillHist("PTj1_gt1j40met30", jetVeto2Coll.at(0).Pt(), weight, 0., 200., 200);
           FillHist("Etaj1_gt1j40met30", jetVeto2Coll.at(0).Eta(), weight, -5., 5., 100);
           FillHist("MET_gt1j40met30", met, weight, 0., 200., 200);
           FillHist("MTW_gt1j40met30", MTW, weight, 0., 200., 200);
           FillHist("MTW_gt1j40met30", MTW, weight, 0., 200., 200);
         }

         //NvtxRW
         FillHist("PTe_gt1j40_NvtxRW", electronColl.at(0).Pt(), weight*nvtx_reweight, 0., 200., 200);
         FillHist("Etae_gt1j40_NvtxRW", electronColl.at(0).Eta(), weight*nvtx_reweight, -5., 5., 100);
         FillHist("Nj_gt1j40_NvtxRW", jetVeto2Coll.size(), weight*nvtx_reweight, 0., 10., 10);
         FillHist("PTj1_gt1j40_NvtxRW", jetVeto2Coll.at(0).Pt(), weight*nvtx_reweight, 0., 200., 200);
         FillHist("Etaj1_gt1j40_NvtxRW", jetVeto2Coll.at(0).Eta(), weight*nvtx_reweight, -5., 5., 100);
         FillHist("MET_gt1j40_NvtxRW", met, weight*nvtx_reweight, 0., 200., 200);
         FillHist("MTW_gt1j40_NvtxRW", MTW, weight*nvtx_reweight, 0., 200., 200);
         if(met>30){
           FillHist("MTW_gt1j40met30_NvtxRW", MTW, weight*nvtx_reweight, 0., 200., 200);
         }
       }

       //NormCheck
       if(met>30 && MTW>50){
         FillHist("Count_NormCR", 0., weight, 0., 5., 5);
         if(jetVeto2Coll.size()>0){
            if(jetVeto2Coll.at(0).Pt()>40){
              FillHist("Count_NormCR", 1., weight, 0., 5., 5);
              FillHist("Count_NormCR", 2., weight*nvtx_reweight, 0., 5., 5);
            }
         }
       }
     }//End of 1 electron 

     if(electronColl.size()==2 && electronFakeL2Coll.size()==2){
       bool Execute=false;
       if( PreTrig && electronColl.at(0).Charge()!=electronColl.at(1).Charge()
                   && electronColl.at(0).Pt()>25. && electronColl.at(1).Pt()>10. 
                   && fabs((electronColl.at(0)+electronColl.at(1)).M()-91.2)<15.   ){ Execute=true; }
       if( UnPreTrig && electronColl.at(0).Charge()!=electronColl.at(1).Charge()
                     && electronColl.at(0).Pt()>30. && electronColl.at(1).Pt()>10. 
                     && fabs((electronColl.at(0)+electronColl.at(1)).M()-91.2)<15. ){ Execute=true;}

       if(Execute){
         FillHist("Mee_e25e10", (electronColl.at(0)+electronColl.at(1)).M(), weight, 60., 120., 60);
         FillHist("Count_NormCR", 3., weight, 0., 5., 5);
         if(jetVeto2Coll.size()>0 && jetVeto2Coll.at(0).Pt()>40){
           FillHist("Mee_e25e10gt1j", (electronColl.at(0)+electronColl.at(1)).M(), weight, 60., 120., 60);
           FillHist("Count_NormCR", 4., weight, 0., 5., 5);
           FillHist("MET_Zwin_e25e10gt1j", met, weight, 0., 200., 200);
           FillHist("MET_Zwin_e25e10gt1j_NvtxRW", met, weight*nvtx_reweight, 0., 200., 200);
           
           FillHist("Nvtx_2OSEle_PURW", Nvtx, weight, 0., 50., 50);
           if(isData){
             FillHist("Nvtx_nocut_data", Nvtx, weight, 0., 50., 50);
             if(njets==0) FillHist("Nvtx_0j_nocut_data", Nvtx, weight, 0., 50., 50);
             if(njets==1) FillHist("Nvtx_1j_nocut_data", Nvtx, weight, 0., 50., 50);
           }
           else{
             FillHist("Nvtx_nocut_mc", Nvtx, weight, 0., 50., 50);
             if(njets==0) FillHist("Nvtx_0j_nocut_mc", Nvtx, weight, 0., 50., 50);
             if(njets==1) FillHist("Nvtx_1j_nocut_mc", Nvtx, weight, 0., 50., 50);
           }

         }
       }
     }//End of DY Count 
     
   }
   if(EleFR){

     std::vector<snu::KJet> jetVeto2Coll = SkimJetColl(jetColl, electronFakeL2Coll, muonLooseColl, "EleMuVeto");
     std::vector<snu::KJet> jetVeto2Iso05Coll = SkimJetColl(jetColl, electronFakeL2Iso05Coll, muonLooseColl, "EleMuVeto");
     //float nvtx_reweight=mcdata_correction->UserPileupWeight(eventbase->GetEvent(), jetVeto2Coll.size());
     //float purw_down=0.; if(pileup_reweight!=0) purw_down=eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::down)/pileup_reweight; 
     //pileup_reweight electronFakeVetoColl

     if(electronFakeL2Coll.size()==1 ) FillHist("IDFlow", 0., weight, 0., 5., 5);
     if(electronFakeL2Coll.size()==1 && electronFakeVetoColl.size()==1) FillHist("IDFlow", 1., weight, 0., 5., 5);
     if(electronFakeL2Coll.size()==1 && electronPreColl.size()==1) FillHist("IDFlow", 2., weight, 0., 5., 5);

     if(electronFakeL2Coll.size()==1 && jetVeto2Coll.size()>0){
       bool PassSel=true; float MTW=0.;
       if( electronFakeL2Coll.at(0).Pt()<25 ) PassSel=false;
       if( PassSel && !(jetVeto2Coll.at(0).Pt()>40) ) PassSel=false;
       if( PassSel && !(electronFakeL2Coll.at(0).DeltaR(jetVeto2Coll.at(0))>1.0) ) PassSel=false;
       if( PassSel ){
         MTW = sqrt(2)*sqrt(met*electronFakeL2Coll.at(0).Pt()-met_x*electronFakeL2Coll.at(0).Px()-met_y*electronFakeL2Coll.at(0).Py());
         
         for(int i=1; i<=20; i++){
           if(met<i*5.){
             for(int j=1; j<=20; j++){
               if(MTW<j*5.){
                 FillHist("NEvt_ltMETltMTW_2D", (i-1)*5., (j-1)*5., weight, 0., 100., 20, 0., 100., 20); 
                 if(PassIDCriteria(electronFakeL2Coll.at(0),"POGMVAMIP")){
                   FillHist("NIDEvt_ltMETltMTW_2D", (i-1)*5., (j-1)*5., weight, 0., 100., 20, 0., 100., 20); 
                 }
               }
             }
           }
         }

         if( !(met<25 && MTW<35) ) PassSel=false;
       }
       if( PassSel ){
         float PTCorr=electronFakeL2Coll.at(0).Pt()*(1.+electronFakeL2Coll.at(0).PFRelIso(0.3));
         float fEta=fabs(electronFakeL2Coll.at(0).Eta());
         const int NPtEdges=7, NEtaEdges=4;
         float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100., 200.};
         float EtaEdges[NEtaEdges]={0., 0.8, 1.479, 2.5};
  
         if( fEta<1.479 ) FillHist("EleEBOpt2SumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
         else             FillHist("EleEEOpt2SumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
         FillHist("EleOpt2SumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
         FillHist("EleOpt2SumW_PTEta_2D", PTCorr, fEta, weight, PtEdges, NPtEdges-1, EtaEdges, NEtaEdges-1);
         if(PassIDCriteria(electronFakeL2Coll.at(0),"POGMVAMIP")){
           if( fEta<1.479 ) FillHist("EleEBOpt2IDSumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
           else             FillHist("EleEEOpt2IDSumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
           FillHist("EleOpt2IDSumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
           FillHist("EleOpt2IDSumW_PTEta_2D", PTCorr, fEta, weight, PtEdges, NPtEdges-1, EtaEdges, NEtaEdges-1);
         }

         //Validation
         FillHist("dRej1", electronFakeL2Coll.at(0).DeltaR(jetVeto2Coll.at(0)),weight, 0., 5., 50);
         FillHist("PTe1", electronFakeL2Coll.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae1", electronFakeL2Coll.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTj1", jetVeto2Coll.at(0).Pt(), weight, 0., 500., 500);
         FillHist("MET", met, weight, 0., 200., 200);
         FillHist("MTW", MTW, weight, 0., 200., 200);
         FillHist("Nj", jetVeto2Coll.size(), weight, 0., 10., 10);
         
       }
     }//Endof Opt2


     if(electronFakeL2Iso05Coll.size()==1 && jetVeto2Iso05Coll.size()>0){
       bool PassSel=true; float MTW=0.;
       if( electronFakeL2Iso05Coll.at(0).Pt()<25 ) PassSel=false;
       if( PassSel && !(jetVeto2Iso05Coll.at(0).Pt()>40) ) PassSel=false;
       if( PassSel && !(electronFakeL2Iso05Coll.at(0).DeltaR(jetVeto2Iso05Coll.at(0))>1.0) ) PassSel=false;
       if( PassSel ){
         MTW = sqrt(2)*sqrt(met*electronFakeL2Iso05Coll.at(0).Pt()-met_x*electronFakeL2Iso05Coll.at(0).Px()-met_y*electronFakeL2Iso05Coll.at(0).Py());
         
         for(int i=1; i<=20; i++){
           if(met<i*5.){
             for(int j=1; j<=20; j++){
               if(MTW<j*5.){
                 FillHist("NEvtIso05_ltMETltMTW_2D", (i-1)*5., (j-1)*5., weight, 0., 100., 20, 0., 100., 20); 
                 if(PassIDCriteria(electronFakeL2Iso05Coll.at(0),"POGMVAMIP")){
                   FillHist("NIDEvtIso05_ltMETltMTW_2D", (i-1)*5., (j-1)*5., weight, 0., 100., 20, 0., 100., 20); 
                 }
               }
             }
           }
         }

         if( !(met<25 && MTW<35) ) PassSel=false;
       }
       if( PassSel ){
         float PTCorr=electronFakeL2Iso05Coll.at(0).Pt()*(1.+electronFakeL2Iso05Coll.at(0).PFRelIso(0.3));
         float fEta=fabs(electronFakeL2Iso05Coll.at(0).Eta());
         const int NPtEdges=7, NEtaEdges=3;
         float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100., 200.};
         float EtaEdges[NEtaEdges]={0., 1.479, 2.5};
  
         if( fEta<1.479 ) FillHist("EleEBOpt2Iso05SumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
         else             FillHist("EleEEOpt2Iso05SumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
         FillHist("EleOpt2Iso05SumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
         FillHist("EleOpt2Iso05SumW_PTEta_2D", PTCorr, fEta, weight, PtEdges, NPtEdges-1, EtaEdges, NEtaEdges-1);
         if(PassIDCriteria(electronFakeL2Iso05Coll.at(0),"POGMVAMIP")){
           if( fEta<1.479 ) FillHist("EleEBOpt2Iso05IDSumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
           else             FillHist("EleEEOpt2Iso05IDSumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
           FillHist("EleOpt2Iso05IDSumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
           FillHist("EleOpt2Iso05IDSumW_PTEta_2D", PTCorr, fEta, weight, PtEdges, NPtEdges-1, EtaEdges, NEtaEdges-1);
         }

         //Validation
         FillHist("dRej1", electronFakeL2Iso05Coll.at(0).DeltaR(jetVeto2Iso05Coll.at(0)),weight, 0., 5., 50);
         FillHist("PTe1", electronFakeL2Iso05Coll.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae1", electronFakeL2Iso05Coll.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTj1", jetVeto2Iso05Coll.at(0).Pt(), weight, 0., 500., 500);
         FillHist("MET", met, weight, 0., 200., 200);
         FillHist("MTW", MTW, weight, 0., 200., 200);
         FillHist("Nj", jetVeto2Iso05Coll.size(), weight, 0., 10., 10);
         
       }
     }//Endof Opt2 Iso05



   }
   if(Closure){
      

//     if(k_running_nonprompt) electronColl=electronFakeL2Coll;
//     if( !(electronFakeL2Coll.size()==1 && muonLooseColl.size()==2) ) return;
//     if( !(electronColl.size()==1 && muonColl.size()==2) ) return;
//     if( !(muonColl.at(0).Charge()!=muonColl.at(1).Charge()) ) return;
//     if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>15 && muonColl.at(1).Pt()>10) ) return;
//     if( fabs(electronColl.at(0).Eta())>2.5 ) return;
//
//     float Mmumu=(muonColl.at(0)+muonColl.at(1)).M();
//     float MTW = sqrt(2)*sqrt(met*electronColl.at(0).Pt()-met_x*electronColl.at(0).Px()-met_y*electronColl.at(0).Py());
//     if(Mmumu<12) return;
    
     
     
   }


/////////////////////////////////////////////////////////////////////////////////// 

return;

}// End of execute event loop
  


void Jul2017_EMuMuComp::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Jul2017_EMuMuComp::BeginCycle() throw( LQError ){
  
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

Jul2017_EMuMuComp::~Jul2017_EMuMuComp() {
  
  Message("In Jul2017_EMuMuComp Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Jul2017_EMuMuComp::FillCutFlow(TString cut, float weight){
  
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



void Jul2017_EMuMuComp::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Jul2017_EMuMuComp::MakeHistograms(){
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
  // *  Remove//Overide this Jul2017_EMuMuCompCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Jul2017_EMuMuComp::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
