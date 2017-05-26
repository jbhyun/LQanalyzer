// $Id: May2017_TmpIDTrigSF.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQMay2017_TmpIDTrigSF Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "May2017_TmpIDTrigSF.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (May2017_TmpIDTrigSF);

 May2017_TmpIDTrigSF::May2017_TmpIDTrigSF() : AnalyzerCore(), out_muons(0) {

   SetLogName("May2017_TmpIDTrigSF");
   Message("In May2017_TmpIDTrigSF constructor", INFO);
   InitialiseAnalysis();
 }


 void May2017_TmpIDTrigSF::InitialiseAnalysis() throw( LQError ) {
   
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

void May2017_TmpIDTrigSF::ExecuteEvents()throw( LQError ){

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


   bool EleIDSF=false, EMuTrigSF=false, MuIDSF=false, MuTrigSF=false;
   bool MuonLeg=false, EleLeg=false, DZeff=false, MCClosure_2l=false, MCClosure_3l=false;
   bool CorrelationStudy=false, EleEffbyMu=false, MuEffbyEle=false, MuEffbyMu=false;
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("EleIDSF"))          EleIDSF=true;
     else if(k_flags.at(i).Contains("EMuTrigSF"))        EMuTrigSF=true;
     else if(k_flags.at(i).Contains("MuonLeg"))          MuonLeg=true;
     else if(k_flags.at(i).Contains("EleLeg"))           EleLeg=true;
     else if(k_flags.at(i).Contains("DZeff"))            DZeff=true;
     else if(k_flags.at(i).Contains("CorrelationStudy")) CorrelationStudy=true;
     else if(k_flags.at(i).Contains("EleEffbyMu"))       EleEffbyMu=true;
     else if(k_flags.at(i).Contains("MuEffbyEle"))       MuEffbyEle=true;
     else if(k_flags.at(i).Contains("MuEffbyMu"))        MuEffbyMu=true;
     else if(k_flags.at(i).Contains("MCClosure_2l"))     MCClosure_2l=true;
     else if(k_flags.at(i).Contains("MCClosure_3l"))     MCClosure_3l=true;
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
   if(EleIDSF || MuonLeg){
     if( PassTrigger("HLT_Ele27_WPTight_Gsf_v") ) Pass_Trigger=true;
   }
   if(EleLeg){
     if( PassTrigger("HLT_IsoMu24_v") ) Pass_Trigger=true;
   }
   if(DZeff){
     if(EMuTrigSF){
       if( PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ) Pass_Trigger=true;
     }
   }
   if(MCClosure_2l || MCClosure_3l || CorrelationStudy) Pass_Trigger=true;

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
     if(EleIDSF)       { if( !(electronPreColl.size()>=2) ) return; }
     if(EMuTrigSF)     { if( !(muonPreColl.size()>=1 && electronPreColl.size()>=1) ) return;}
     if(EMuTrigSF && MCClosure_3l){
                         if( !(muonPreColl.size()>=2 && electronPreColl.size()>=1) ) return;
     }
     if(EMuTrigSF && CorrelationStudy){
                         if( !( (muonPreColl.size()+electronPreColl.size())>=3 ) ) return;
     }

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
     eventbase->GetMuonSel()->SetBSdxy(0.05);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.1);
   std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,false);
   std::vector<snu::KMuon> muonColl;      if(k_running_nonprompt){ muonColl=muonLooseColl;} else{ muonColl=muonTightColl;}


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
     //eventbase->GetElectronSel()->SetCheckCharge(true);
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
   std::vector<snu::KElectron> electronTightColl; eventbase->GetElectronSel()->Selection(electronTightColl);
   std::vector<snu::KElectron> electronColl;      if(k_running_nonprompt){ electronColl=electronLooseColl;} else{ electronColl=electronTightColl;}
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
   int NTighte=0, NTightmu=0;
   if(true){
     if(!isData){
      //trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );
  
//       id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronColl);
//       reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
  
//       id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muonColl);
       //iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
//       trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
  
//       btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

       geneff_weight   = GenFilterEfficiency(k_sample_name);
       gennorm_weight  = SignalNorm(k_sample_name, 200.);
     }
     else{
       //Perfect Prompt Ratio Approximation applied.(p=1), Applicable to generic number, combination of leptons under premise of high prompt rate.
       if(k_running_nonprompt){
         //fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", muonLooseColl.size(), electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", 0);
         for(int i=0; i<electronLooseColl.size(); i++){
            if( ( fabs(electronLooseColl.at(i).Eta())<1.447 && electronLooseColl.at(i).PFRelIso(0.3)<0.0588 )
                 || (fabs(electronLooseColl.at(i).Eta())>1.447 && electronLooseColl.at(i).PFRelIso(0.3)<0.0571 ) ){
            NTighte++;
            }
         }
         for(int i=0; i<muonLooseColl.size(); i++){
            if( muonLooseColl.at(i).RelIso04()<0.1 ){
            NTightmu++;
            }
         }

         if(NTighte==0)      {fake_weight=-1.; }
         else if(NTighte==1) {fake_weight=1.; }
         else if(NTighte==2) {fake_weight=0.; }
         else return;

         FillHist("NTighte", NTighte, 1, 0., 3., 3); 
         FillHist("NTightmu", NTightmu, 1, 0., 3., 3); 
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
   FillHist("Basic_METdist_orig", met, weight, 0., 200., 100);
   ///////////////////////////////////////////////////////////////////////////////////


/************************************************************************************/
//////////////////////////////////////////////////////////////////////////////////////
/////// MAIN ANALYSIS CODE ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
/************************************************************************************/


   if(EleIDSF){

     if( electronColl.size()!=2 ) return;
     if( electronColl.at(0).Charge()==electronColl.at(1).Charge() ) return;
     if( fabs((electronColl.at(0)+electronColl.at(1)).M()-90)>30 ) return;
     if( !(electronColl.at(0).Pt()>30 && electronColl.at(1).Pt()>20) ) return;
     FillCutFlow("BeforeTrigMatch",weight);
     int IdxTag=-1, IdxProbe=-1;
     if( electronColl.at(0).TriggerMatched("HLT_Ele27_WPTight_Gsf_v") ) IdxTag=0;
     else if( electronColl.at(1).TriggerMatched("HLT_Ele27_WPTight_Gsf_v") 
              && electronColl.at(1).Pt()>30 ) IdxTag=1;
     else return; FillCutFlow("AfterTrigMatch", weight);

     IdxProbe=1-IdxTag;


     FillHist("NTighteSel", NTighte, 1, 0., 3., 3); 
     
     //ID efficiency
     float EleXbinEdges[9]={20.,30.,40.,50.,60.,80.,120.,200.,500.};
     float EleYbinEdges[9]={-2.5, -2.0, -1.5, -1.0, 0., 1.0, 1.5, 2.0, 2.5};
     int EleNbinsX=8, EleNbinsY=8;

     FillHist("ID_Ne_POGT_PtEta", electronColl.at(IdxProbe).Pt(), electronColl.at(IdxProbe).Eta(), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);
     FillHist("ID_Ne_POGT_Pt", electronColl.at(IdxProbe).Pt(), weight, EleXbinEdges, EleNbinsX);
     FillHist("ID_Ne_POGT_Eta", electronColl.at(IdxProbe).Eta(), weight, EleYbinEdges, EleNbinsY);
     if( electronColl.at(IdxProbe).dxy()<0.05 && electronColl.at(IdxProbe).dz()<0.1
        && electronColl.at(IdxProbe).dxySig()<3. ){
        FillHist("ID_Ne_POGTIP_PtEta", electronColl.at(IdxProbe).Pt(), electronColl.at(IdxProbe).Eta(), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);
        FillHist("ID_Ne_POGTIP_Pt", electronColl.at(IdxProbe).Pt(), weight, EleXbinEdges, EleNbinsX);
        FillHist("ID_Ne_POGTIP_Eta", electronColl.at(IdxProbe).Eta(), weight, EleYbinEdges, EleNbinsY);
     }


   }
   if(EMuTrigSF){
     if(EleLeg){

       if( !(electronColl.size()==1 && muonColl.size()==1) ) return;
       if( !(muonColl.at(0).Pt()>27) ) return;
       if( !(electronColl.at(0).DeltaR(muonColl.at(0))>0.3) ) return;
  
       if( !muonColl.at(0).TriggerMatched("HLT_IsoMu24_v") ) return;
  
       //Trigger Leg efficiency       
       float EleXbinEdges[13]={0.,10.,20.,23.,25.,30.,40.,50.,60.,80.,120.,200.,500.};
       float EleYbinEdges[5]={0., 1.0, 1.5, 2.0, 2.5};
       int EleNbinsX=12, EleNbinsY=4;
       int IdxProbe=0;
       if( electronColl.at(IdxProbe).dxy()<0.05 && electronColl.at(IdxProbe).dz()<0.1
          && electronColl.at(IdxProbe).dxySig()<3. ){
          FillHist("Tr_Ne_POGTIP_PtEta", electronColl.at(IdxProbe).Pt(), fabs(electronColl.at(IdxProbe).Eta()), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);
          FillHist("Tr_Ne_POGTIP_Pt", electronColl.at(IdxProbe).Pt(), weight, EleXbinEdges, EleNbinsX);
          FillHist("Tr_Ne_POGTIP_Eta", fabs(electronColl.at(IdxProbe).Eta()), weight, EleYbinEdges, EleNbinsY);
          if( electronColl.at(0).TriggerMatched("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
             || electronColl.at(0).TriggerMatched("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){
            FillHist("Tr_Ne_POGTIPEle23_PtEta", electronColl.at(IdxProbe).Pt(), fabs(electronColl.at(IdxProbe).Eta()), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);
            FillHist("Tr_Ne_POGTIPEle23_Pt", electronColl.at(IdxProbe).Pt(), weight, EleXbinEdges, EleNbinsX);
            FillHist("Tr_Ne_POGTIPEle23_Eta", fabs(electronColl.at(IdxProbe).Eta()), weight, EleYbinEdges, EleNbinsY);
          }
       }

     }//End of EleLeg
     if(MuonLeg){

       if( !(electronVetoColl.size()==1 && muonPOGLColl.size()==1) ) return;
       if( !(electronColl.size()==1 && muonColl.size()==1) ) return;
       if( !(electronColl.at(0).Pt()>30) ) return;
       if( !(electronColl.at(0).DeltaR(muonColl.at(0))>0.3) ) return;
       if( !(electronColl.at(0).dxy()<0.05 && electronColl.at(0).dz()<0.1 && electronColl.at(0).dxySig()<3. ) ) return;
       if( !electronColl.at(0).TriggerMatched("HLT_Ele27_WPTight_Gsf_v") ) return; 

  
  
       //Trigger Leg efficiency
       int DataPeriodCheck = -1;//-1: MC, 1:BtoF, 2:GtoH
       if(isData){
         DataPeriodCheck=GetDataPeriod();
       }


       float MuXbinEdges[15]={0.,5.,8.,10.,12.,15.,20.,30.,40.,50.,60.,80.,120.,200.,500.};
       float Mu1XbinEdges[14]={0.,5.,8.,10.,15.,20.,30.,40.,50.,60.,80.,120.,200.,500.};
       float Mu2XbinEdges[14]={0.,8.,10.,12.,15.,20.,30.,40.,50.,60.,80.,120.,200.,500.};
       float MuYbinEdges[5]={0., 1.0, 1.5, 2.0, 2.4};
       int Mu1NbinsX=13, Mu2NbinsX=13, MuNbinsX=14, MuNbinsY=4;
       int IdxProbe=0;

       FillHist("Ptmu", muonColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("Etamu", muonColl.at(0).Eta(), weight, -2.5, 2.5, 100);

       //Only Mu8 Leg (<=PeriodG)
       if(DataPeriodCheck<0 || DataPeriodCheck<=6){
         FillHist("Tr_Nmu8_HNID1_PtEta", muonColl.at(IdxProbe).Pt(), fabs(muonColl.at(IdxProbe).Eta()), weight, Mu1XbinEdges, Mu1NbinsX, MuYbinEdges, MuNbinsY);
         FillHist("Tr_Nmu8_HNID1_Pt", muonColl.at(IdxProbe).Pt(), weight, Mu1XbinEdges, Mu1NbinsX);
         FillHist("Tr_Nmu8_HNID1_Eta", fabs(muonColl.at(IdxProbe).Eta()), weight, MuYbinEdges, MuNbinsY);

         //if( muonColl.at(IdxProbe).TriggerMatched("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){
         if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){
           FillHist("Tr_Nmu8_HNID1Mu8_PtEta", muonColl.at(IdxProbe).Pt(), fabs(muonColl.at(IdxProbe).Eta()), weight, Mu1XbinEdges, Mu1NbinsX, MuYbinEdges, MuNbinsY);
           FillHist("Tr_Nmu8_HNID1Mu8_Pt", muonColl.at(IdxProbe).Pt(), weight, Mu1XbinEdges, Mu1NbinsX);
           FillHist("Tr_Nmu8_HNID1Mu8_Eta", fabs(muonColl.at(IdxProbe).Eta()), weight, MuYbinEdges, MuNbinsY);
         }
       }
       //Only Mu12 Leg(>=PeriodG)
       if(DataPeriodCheck<0 || DataPeriodCheck==6){

         //Only Mu12 Leg(>=PeriodG)
         FillHist("Tr_Nmu12_HNID1_PtEta", muonColl.at(IdxProbe).Pt(), fabs(muonColl.at(IdxProbe).Eta()), weight, Mu2XbinEdges, Mu2NbinsX, MuYbinEdges, MuNbinsY);
         FillHist("Tr_Nmu12_HNID1_Pt", muonColl.at(IdxProbe).Pt(), weight, Mu2XbinEdges, Mu2NbinsX);
         FillHist("Tr_Nmu12_HNID1_Eta", fabs(muonColl.at(IdxProbe).Eta()), weight, MuYbinEdges, MuNbinsY);

         //if( muonColl.at(IdxProbe).TriggerMatched("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){
         if( PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){
           FillHist("Tr_Nmu12_HNID1Mu12_PtEta", muonColl.at(IdxProbe).Pt(), fabs(muonColl.at(IdxProbe).Eta()), weight, Mu2XbinEdges, Mu2NbinsX, MuYbinEdges, MuNbinsY);
           FillHist("Tr_Nmu12_HNID1Mu12_Pt", muonColl.at(IdxProbe).Pt(), weight, Mu2XbinEdges, Mu2NbinsX);
           FillHist("Tr_Nmu12_HNID1Mu12_Eta", fabs(muonColl.at(IdxProbe).Eta()), weight, MuYbinEdges, MuNbinsY);
         }

         FillHist("Tr_NmuG_HNID1_PtEta", muonColl.at(IdxProbe).Pt(), fabs(muonColl.at(IdxProbe).Eta()), weight, MuXbinEdges, MuNbinsX, MuYbinEdges, MuNbinsY);
         FillHist("Tr_NmuG_HNID1_Pt", muonColl.at(IdxProbe).Pt(), weight, MuXbinEdges, MuNbinsX);
         FillHist("Tr_NmuG_HNID1_Eta", fabs(muonColl.at(IdxProbe).Eta()), weight, MuYbinEdges, MuNbinsY);
         if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
            || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){

         //if( muonColl.at(IdxProbe).TriggerMatched("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
         //   || muonColl.at(IdxProbe).TriggerMatched("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){
            FillHist("Tr_NmuG_HNID1Mu8orMu12_PtEta", muonColl.at(IdxProbe).Pt(), fabs(muonColl.at(IdxProbe).Eta()), weight, MuXbinEdges, MuNbinsX, MuYbinEdges, MuNbinsY);
            FillHist("Tr_NmuG_HNID1Mu8orMu12_Pt", muonColl.at(IdxProbe).Pt(), weight, MuXbinEdges, MuNbinsX);
            FillHist("Tr_NmuG_HNID1Mu8orMu12_Eta", fabs(muonColl.at(IdxProbe).Eta()), weight, MuYbinEdges, MuNbinsY);
         }
       }


       FillHist("Tr_Nmu_HNID1_PtEta", muonColl.at(IdxProbe).Pt(), fabs(muonColl.at(IdxProbe).Eta()), weight, MuXbinEdges, MuNbinsX, MuYbinEdges, MuNbinsY);
       FillHist("Tr_Nmu_HNID1_Pt", muonColl.at(IdxProbe).Pt(), weight, MuXbinEdges, MuNbinsX);
       FillHist("Tr_Nmu_HNID1_Eta", fabs(muonColl.at(IdxProbe).Eta()), weight, MuYbinEdges, MuNbinsY);
       if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
          || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){

       //if( muonColl.at(IdxProbe).TriggerMatched("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
       //   || muonColl.at(IdxProbe).TriggerMatched("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){
          FillHist("Tr_Nmu_HNID1MuLeg_PtEta", muonColl.at(IdxProbe).Pt(), fabs(muonColl.at(IdxProbe).Eta()), weight, MuXbinEdges, MuNbinsX, MuYbinEdges, MuNbinsY);
          FillHist("Tr_Nmu_HNID1MuLeg_Pt", muonColl.at(IdxProbe).Pt(), weight, MuXbinEdges, MuNbinsX);
          FillHist("Tr_Nmu_HNID1MuLeg_Eta", fabs(muonColl.at(IdxProbe).Eta()), weight, MuYbinEdges, MuNbinsY);
       }
       
     }//End of Muon Leg
     if(DZeff){

       float EleXbinEdges[9]={25.,30.,40.,50.,60.,80.,120.,200.,500.};
       float EleYbinEdges[5]={0., 1.0, 1.5, 2.0, 2.5};

       float MuXbinEdges[10]={15.,20.,30.,40.,50.,60.,80.,120.,200.,500.};
       float MuYbinEdges[5]={0., 1.0, 1.5, 2.0, 2.4};
       int MuNbinsX=9, MuNbinsY=4;
       int EleNbinsX=8, EleNbinsY=4;

       int DataPeriodCheck = -1;//-1: MC, 1:BtoF, 2:GtoH
       if(isData){
         if(eventbase->GetEvent().RunNumber()<278809) DataPeriodCheck=1;//BtoF
         else if(eventbase->GetEvent().RunNumber()<280919) DataPeriodCheck=2;//G
         else DataPeriodCheck=3;//H
       }

       if( !(DataPeriodCheck<0 || DataPeriodCheck==2) )              return;
       if( !(electronColl.size()==1 && muonColl.size()==1) )         return;
       if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>15) ) return;
       if( !(electronColl.at(0).dxy()<0.05 && electronColl.at(0).dz()<0.1 && electronColl.at(0).dxySig()<3. ) ) return;
       if( !(electronColl.at(0).DeltaR(muonColl.at(0))>0.3) )        return;
       if( !(PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")) ) return;

       float dz=fabs(electronColl.at(0).dz()-muonColl.at(0).dZ());
       FillHist("Tr_Nemu_POGTIPEleHNID1MuTrig_Tot", 0., weight, 0., 1., 1);
       FillHist("Tr_Nemu_POGTIPEleHNID1MuTrig_DZ", dz, weight, 0., 0.3, 300);
       FillHist("Tr_Nemu_POGTIPEleHNID1MuTrig_dRemu", electronColl.at(0).DeltaR(muonColl.at(0)), weight, 0., 5., 50);
       FillHist("Tr_Nemu_POGTIPEleHNID1MuTrig_ElePt", electronColl.at(0).Pt(), weight, EleXbinEdges, EleNbinsX);
       FillHist("Tr_Nemu_POGTIPEleHNID1MuTrig_EleEta", fabs(electronColl.at(0).Eta()), weight, EleYbinEdges, EleNbinsY);
       FillHist("Tr_Nemu_POGTIPEleHNID1MuTrig_MuPt", muonColl.at(0).Pt(), weight, MuXbinEdges, MuNbinsX);
       FillHist("Tr_Nemu_POGTIPEleHNID1MuTrig_MuEta", fabs(muonColl.at(0).Eta()), weight, MuYbinEdges, MuNbinsY);
       FillHist("Tr_Nemu_POGTIPEleHNID1MuTrig_MuEtaEleEta", fabs(muonColl.at(0).Eta()), fabs(electronColl.at(0).Eta()), weight, MuYbinEdges, MuNbinsY, EleYbinEdges, EleNbinsY);
       if( PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ){
         FillHist("Tr_Nemu_POGTIPEleHNID1MuTrigdz_Tot", 0., weight, 0., 1., 1);
         FillHist("Tr_Nemu_POGTIPEleHNID1MuTrigdz_DZ", dz, weight, 0., 0.3, 300);
         FillHist("Tr_Nemu_POGTIPEleHNID1MuTrigdz_dRemu", electronColl.at(0).DeltaR(muonColl.at(0)), weight, 0., 5., 50);
         FillHist("Tr_Nemu_POGTIPEleHNID1MuTrigdz_ElePt", electronColl.at(0).Pt(), weight, EleXbinEdges, EleNbinsX);
         FillHist("Tr_Nemu_POGTIPEleHNID1MuTrigdz_EleEta", fabs(electronColl.at(0).Eta()), weight, EleYbinEdges, EleNbinsY);
         FillHist("Tr_Nemu_POGTIPEleHNID1MuTrigdz_MuPt", muonColl.at(0).Pt(), weight, MuXbinEdges, MuNbinsX);
         FillHist("Tr_Nemu_POGTIPEleHNID1MuTrigdz_MuEta", fabs(muonColl.at(0).Eta()), weight, MuYbinEdges, MuNbinsY);
         FillHist("Tr_Nemu_POGTIPEleHNID1MuTrigdz_MuEtaEleEta", fabs(muonColl.at(0).Eta()), fabs(electronColl.at(0).Eta()), weight, MuYbinEdges, MuNbinsY, EleYbinEdges, EleNbinsY);
       }
       
     }//End of DZeff
     if(CorrelationStudy){
       if(EleEffbyMu){

         FillHist("Cutflow", 0., weight, 0., 10., 10);
         if( !(electronColl.size()==1 && muonColl.size()==2) ) return;
         FillHist("Cutflow", 1., weight, 0., 10., 10);
         if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>15 && muonColl.at(1).Pt()>10) ) return;
         FillHist("Cutflow", 2., weight, 0., 10., 10);
        
         int IdxTag= electronColl.at(0).DeltaR(muonColl.at(0))>electronColl.at(0).DeltaR(muonColl.at(1)) ? 0:1;
         int IdxProbe=1-IdxTag;
         bool Pass=true;
         if( muonColl.at(0).DeltaR(muonColl.at(1))<0.4 )            Pass=false;
         if(Pass) FillHist("Cutflow", 3., weight, 0., 10., 10);
         if( electronColl.at(0).DeltaR(muonColl.at(IdxTag))<0.4 )   Pass=false;
         if(Pass) FillHist("Cutflow", 4., weight, 0., 10., 10);
         if( !muonColl.at(IdxTag).TriggerMatched("HLT_IsoMu24_v") ) Pass=false;
         if(Pass) FillHist("Cutflow", 5., weight, 0., 10., 10);
         if( muonColl.at(IdxTag).Pt()<27. )                         Pass=false;
         if(Pass) FillHist("Cutflow", 6., weight, 0., 10., 10);
         if( electronColl.at(0).DeltaR(muonColl.at(IdxProbe))>0.5 ) Pass=false;
         if(Pass) FillHist("Cutflow", 7., weight, 0., 10., 10);
         if( !Pass ) return;
         if(Pass) FillHist("Cutflow", 8., weight, 0., 10., 10);


         float EleXbinEdges[9]={25.,30.,40.,50.,60.,80.,120.,200.,500.};
         float EleYbinEdges[5]={0., 1.0, 1.5, 2.0, 2.5};
         int EleNbinsX=8, EleNbinsY=4;

         float dRemu=electronColl.at(0).DeltaR(muonColl.at(IdxProbe)); 

         FillHist("Tr_Ne_POGTIP_dRemu", dRemu, weight, 0., 0.5, 50);
         if( dRemu<0.3 ){
           FillHist("Tr_Ne_POGTIPdRemult03_Pt", electronColl.at(0).Pt(), weight, EleXbinEdges, EleNbinsX);
           FillHist("Tr_Ne_POGTIPdRemult03_PtEta", electronColl.at(0).Pt(), electronColl.at(0).Eta(), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);
         }
         if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){
           FillHist("Tr_Ne_POGTIPEle23_dRemu", dRemu, weight, 0., 0.5, 50);
           if( dRemu<0.3 ){
             FillHist("Tr_Ne_POGTIPdRemult03Ele23_Pt", electronColl.at(0).Pt(), weight, EleXbinEdges, EleNbinsX);
             FillHist("Tr_Ne_POGTIPdRemult03Ele23_PtEta", electronColl.at(0).Pt(), electronColl.at(0).Eta(), weight, EleXbinEdges, EleNbinsX, EleYbinEdges, EleNbinsY);
           }
         }
       }//End of EleLeg By Mu
       if(MuEffbyEle){
         if( !(electronColl.size()==2 && muonColl.size()==1) ) return;
         if( !(electronColl.at(0).Pt()>25 && electronColl.at(1).Pt()>25 && muonColl.at(0).Pt()>10) ) return;
        
         int IdxTag= muonColl.at(0).DeltaR(electronColl.at(0))>muonColl.at(0).DeltaR(electronColl.at(1)) ? 0:1;
         int IdxProbe=1-IdxTag;
         bool Pass=true;
         if( electronColl.at(0).DeltaR(electronColl.at(1))<0.3 )    Pass=false;
         if( muonColl.at(0).DeltaR(electronColl.at(IdxTag))<0.3 )   Pass=false;
         if( !electronColl.at(IdxTag).TriggerMatched("HLT_Ele27_WPTight_Gsf_v") ) Pass=false;
         if( electronColl.at(IdxTag).Pt()<30. )                     Pass=false;
         if( muonColl.at(0).DeltaR(electronColl.at(IdxProbe))>0.5 ) Pass=false;
         if( !Pass ) return;


         float MuXbinEdges[12]={10.,12.,15.,20.,30.,40.,50.,60.,80.,120.,200.,500.};
         float Mu1XbinEdges[14]={0.,5.,8.,10.,15.,20.,30.,40.,50.,60.,80.,120.,200.,500.};
         float Mu2XbinEdges[14]={0.,8.,10.,12.,15.,20.,30.,40.,50.,60.,80.,120.,200.,500.};
         float MuYbinEdges[5]={0., 1.0, 1.5, 2.0, 2.4};
         int Mu1NbinsX=13, Mu2NbinsX=13, MuNbinsX=12, MuNbinsY=4;

         float dRemu=muonColl.at(0).DeltaR(electronColl.at(IdxProbe)); 

         FillHist("Tr_Nmu8_HNID1_dRemu", dRemu, weight, 0., 0.5, 50);
         FillHist("Tr_Nmu12_HNID1_dRemu", dRemu, weight, 0., 0.5, 50);
         if( dRemu<0.3 ){
           FillHist("Tr_Nmu8_HNID1dRemult03_Pt", muonColl.at(0).Pt(), weight, Mu1XbinEdges, Mu1NbinsX);
           FillHist("Tr_Nmu8_HNID1dRemult03_PtEta", muonColl.at(0).Pt(), muonColl.at(0).Eta(), weight, MuXbinEdges, MuNbinsX, MuYbinEdges, MuNbinsY);
           FillHist("Tr_Nmu12_HNID1dRemult03_Pt", muonColl.at(0).Pt(), weight, Mu2XbinEdges, Mu2NbinsX);
           FillHist("Tr_Nmu12_HNID1dRemult03_PtEta", muonColl.at(0).Pt(), muonColl.at(0).Eta(), weight, MuXbinEdges, MuNbinsX, MuYbinEdges, MuNbinsY);
         }
         if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){
           FillHist("Tr_Nmu8_HNID1Mu8_dRemu", dRemu, weight, 0., 0.5, 50);
           if( dRemu<0.3 ){
             FillHist("Tr_Nmu8_HNID1dRemult03Mu8_Pt", muonColl.at(0).Pt(), weight, Mu1XbinEdges, Mu1NbinsX);
             FillHist("Tr_Nmu8_HNID1dRemult03Mu8_PtEta", muonColl.at(0).Pt(), muonColl.at(0).Eta(), weight, MuXbinEdges, MuNbinsX, MuYbinEdges, MuNbinsY);
           }
         }
         if( PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){
           FillHist("Tr_Nmu12_HNID1Mu12_dRemu", dRemu, weight, 0., 0.5, 50);
           if( dRemu<0.3 ){
             FillHist("Tr_Nmu12_HNID1dRemult03Mu12_Pt", muonColl.at(0).Pt(), weight, Mu2XbinEdges, Mu2NbinsX);
             FillHist("Tr_Nmu12_HNID1dRemult03Mu12_PtEta", muonColl.at(0).Pt(), muonColl.at(0).Eta(), weight, MuXbinEdges, MuNbinsX, MuYbinEdges, MuNbinsY);
           }
         }

       }//End of EleLeg By Mu
       if(MuEffbyMu){
         
         if( !(electronColl.size()==1 && muonColl.size()==2) ) return;
         if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>15 && muonColl.at(1).Pt()>10) ) return;
        
         int IdxTag=0;
         int IdxProbe=1;
         bool Pass=true;
         if( electronColl.at(0).DeltaR(muonColl.at(IdxTag))<0.3 )            Pass=false;
         if( electronColl.at(0).DeltaR(muonColl.at(IdxProbe))<0.3 )          Pass=false;
         if( !electronColl.at(0).TriggerMatched("HLT_Ele27_WPTight_Gsf_v") ) Pass=false;
         if( electronColl.at(0).Pt()<30. )                                   Pass=false;
         if( muonColl.at(0).DeltaR(muonColl.at(1))>0.5 )                     Pass=false;

         if( !Pass ) return;


         float MuXbinEdges[12]={10.,12.,15.,20.,30.,40.,50.,60.,80.,120.,200.,500.};
         float Mu1XbinEdges[14]={0.,5.,8.,10.,15.,20.,30.,40.,50.,60.,80.,120.,200.,500.};
         float Mu2XbinEdges[14]={0.,8.,10.,12.,15.,20.,30.,40.,50.,60.,80.,120.,200.,500.};
         float MuYbinEdges[5]={0., 1.0, 1.5, 2.0, 2.4};
         int Mu1NbinsX=13, Mu2NbinsX=13, MuNbinsX=12, MuNbinsY=4;

         float dRmumu=muonColl.at(0).DeltaR(muonColl.at(1)); 


         FillHist("Tr_Nmu8_HNID1_dRmumu", dRmumu, weight, 0., 0.5, 50);
         FillHist("Tr_Nmu12_HNID1_dRmumu", dRmumu, weight, 0., 0.5, 50);
         if( dRmumu<0.3 ){
           FillHist("Tr_Nmu8_HNID1dRmumult03_Pt", muonColl.at(0).Pt(), weight, Mu1XbinEdges, Mu1NbinsX);
           FillHist("Tr_Nmu8_HNID1dRmumult03_PtEta", muonColl.at(0).Pt(), muonColl.at(0).Eta(), weight, MuXbinEdges, MuNbinsX, MuYbinEdges, MuNbinsY);
           FillHist("Tr_Nmu12_HNID1dRmumult03_Pt", muonColl.at(0).Pt(), weight, Mu2XbinEdges, Mu2NbinsX);
           FillHist("Tr_Nmu12_HNID1dRmumult03_PtEta", muonColl.at(0).Pt(), muonColl.at(0).Eta(), weight, MuXbinEdges, MuNbinsX, MuYbinEdges, MuNbinsY);
         }
         if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){
           FillHist("Tr_Nmu8_HNID1Mu8_dRmumu", dRmumu, weight, 0., 0.5, 50);
           if( dRmumu<0.3 ){
             FillHist("Tr_Nmu8_HNID1dRmumult03Mu8_Pt", muonColl.at(0).Pt(), weight, Mu1XbinEdges, Mu1NbinsX);
             FillHist("Tr_Nmu8_HNID1dRmumult03Mu8_PtEta", muonColl.at(0).Pt(), muonColl.at(0).Eta(), weight, MuXbinEdges, MuNbinsX, MuYbinEdges, MuNbinsY);
           }
         }
         if( PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){
           FillHist("Tr_Nmu12_HNID1Mu12_dRmumu", dRmumu, weight, 0., 0.5, 50);
           if( dRmumu<0.3 ){
             FillHist("Tr_Nmu12_HNID1dRmumult03Mu12_Pt", muonColl.at(0).Pt(), weight, Mu2XbinEdges, Mu2NbinsX);
             FillHist("Tr_Nmu12_HNID1dRmumult03Mu12_PtEta", muonColl.at(0).Pt(), muonColl.at(0).Eta(), weight, MuXbinEdges, MuNbinsX, MuYbinEdges, MuNbinsY);
           }
         }



       }

     }//End of Correlation Study
     if(MCClosure_2l){
     
       if( !(electronColl.size()==1 && muonColl.size()==1) )         return;
       if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>10) ) return;
       if( !(electronColl.at(0).dxy()<0.05 && electronColl.at(0).dz()<0.1 && electronColl.at(0).dxySig()<3. ) ) return;
       //For justification of methodology for dRemu>0.3 where correlation is limited
       //if( !(electronColl.at(0).DeltaR(muonColl.at(0))>0.3) )        return;
     
       //Expectation
       float EleLegEff=1., MuLegEff=1., DzEff=0.99;
       float ElePt=electronColl.at(0).Pt(), MuPt=muonColl.at(0).Pt();
       if     (ElePt<30.)  EleLegEff=0.941;
       else if(ElePt<40.)  EleLegEff=0.957;
       else if(ElePt<50.)  EleLegEff=0.973;
       else if(ElePt<60.)  EleLegEff=0.98;
       else if(ElePt<80.)  EleLegEff=0.985;
       else if(ElePt<120.) EleLegEff=0.988;
       else                EleLegEff=0.99;

       if     (MuPt<15.)   MuLegEff=0.922;
       else if(MuPt<20.)   MuLegEff=0.938;
       else if(MuPt<30.)   MuLegEff=0.946;
       else if(MuPt<60.)   MuLegEff=0.95;
       else                MuLegEff=0.948;

       float TrigEffMu8Ele23=1., TrigEffMu12Ele23=1., TrigEffMu12Ele23Dz=1.;
       TrigEffMu8Ele23=EleLegEff*MuLegEff;
       TrigEffMu12Ele23=EleLegEff*MuLegEff;
       TrigEffMu12Ele23Dz=EleLegEff*MuLegEff*DzEff;

       //PT10 - Mu8Ele23 Closure
       FillHist("MCClosure_TotPt10", 0., weight, 0., 2., 2);
       FillHist("MCClosure_TotPt10", 1., weight, 0., 2., 2);

       FillHist("MCClosure_ExpObs_Mu8Ele23Pt10", 0., weight*TrigEffMu8Ele23, 0., 2., 2);
       if(  PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){
         FillHist("MCClosure_ExpObs_Mu8Ele23Pt10", 1., weight, 0., 2., 2);
       }

       //PT15 - Mu12Ele23 Closure, Mu8Ele23 Closure
       if(MuPt>15){
         FillHist("MCClosure_TotPt15", 0., weight, 0., 2., 2);
         FillHist("MCClosure_TotPt15", 1., weight, 0., 2., 2);

         FillHist("MCClosure_ExpObs_Mu8Ele23Pt15", 0., weight*TrigEffMu8Ele23, 0., 2., 2);
         FillHist("MCClosure_ExpObs_Mu12Ele23DzPt15", 0., weight*TrigEffMu12Ele23Dz, 0., 2., 2);
         if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){
           FillHist("MCClosure_ExpObs_Mu8Ele23Pt15", 1., weight, 0., 2., 2);
         }
         if( PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ){
           FillHist("MCClosure_ExpObs_Mu12Ele23DzPt15", 1., weight, 0., 2., 2);
         }

       }

//       Probably Not Going to Use. But This was used when considering oring of all 3 triggers.
//       But considering non-zero prescale, this will cause additional systematic. So I will not do.
//       if(  PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
//          || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
//          || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") )
//       {
//         FillHist("MCClosure_ExpObs", 1., weight, 0., 2., 2);
//       }

      
     }//End of MCClosure(DiLep)

     if(MCClosure_3l){

       if( !(electronColl.size()==1 && muonColl.size()==2) )                                                    return;
       if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>15 && muonColl.at(1).Pt()>10) )                  return;
       if( !(electronColl.at(0).dxy()<0.05 && electronColl.at(0).dz()<0.1 && electronColl.at(0).dxySig()<3. ) ) return;

       float dRemuS=min(electronColl.at(0).DeltaR(muonColl.at(0)), electronColl.at(0).DeltaR(muonColl.at(1)));
       float dRemuL=max(electronColl.at(0).DeltaR(muonColl.at(0)), electronColl.at(0).DeltaR(muonColl.at(1)));
       float dRmumu=muonColl.at(0).DeltaR(muonColl.at(1));
  
       FillHist("dRemuS", dRemuS, weight, 0., 5., 500);
       FillHist("dRemuL", dRemuL, weight, 0., 5., 500);
       FillHist("dRmumu", dRmumu, weight, 0., 5., 500);
       //For justification of methodology for dRemu>0.3 where correlation is limited
       if( !(electronColl.at(0).DeltaR(muonColl.at(0))>0.3 && electronColl.at(0).DeltaR(muonColl.at(1))>0.3 ) ) return;
       //if( !(muonColl.at(0).DeltaR(muonColl.at(1))>0.3) ) return;
     
       //Expectation -Mu8Ele23
       float EleLegEff=1., Mu1LegEff8=1., Mu2LegEff8=1., DzEff=0.99;
       float ElePt=electronColl.at(0).Pt(), Mu1Pt=muonColl.at(0).Pt(), Mu2Pt=muonColl.at(1).Pt();
       if     (ElePt<30.)  EleLegEff=0.941;
       else if(ElePt<40.)  EleLegEff=0.957;
       else if(ElePt<50.)  EleLegEff=0.973;
       else if(ElePt<60.)  EleLegEff=0.98;
       else if(ElePt<80.)  EleLegEff=0.985;
       else if(ElePt<120.) EleLegEff=0.988;
       else                EleLegEff=0.99;

       if     (Mu1Pt<15.)   Mu1LegEff8=0.922;
       else if(Mu1Pt<20.)   Mu1LegEff8=0.938;
       else if(Mu1Pt<30.)   Mu1LegEff8=0.946;
       else if(Mu1Pt<60.)   Mu1LegEff8=0.95;
       else                 Mu1LegEff8=0.948;

       if     (Mu2Pt<15.)   Mu2LegEff8=0.922;
       else if(Mu2Pt<20.)   Mu2LegEff8=0.938;
       else if(Mu2Pt<30.)   Mu2LegEff8=0.946;
       else if(Mu2Pt<60.)   Mu2LegEff8=0.95;
       else                 Mu2LegEff8=0.948;

       //Expectation - Mu12Ele23Dz
       float Mu2LegEff12=1.;
       if     (Mu2Pt<12.) Mu2LegEff12=0.04;
       else if(Mu2Pt<15.) Mu2LegEff12=0.914;
       else               Mu2LegEff12=Mu2LegEff8;


       //Expectation from Data efficiency(Check SF impact)
       //Expectation -Mu8Ele23
       float EleLegEffData=1., Mu1LegEff8Data=1., Mu2LegEff8Data=1., DzEffData=0.964;
       if     (ElePt<30.)  EleLegEffData=0.922;
       else if(ElePt<40.)  EleLegEffData=0.952;
       else if(ElePt<50.)  EleLegEffData=0.962;
       else if(ElePt<60.)  EleLegEffData=0.971;
       else if(ElePt<80.)  EleLegEffData=0.974;
       else if(ElePt<120.) EleLegEffData=0.977;
       else                EleLegEffData=0.975;

       if     (Mu1Pt<15.)   Mu1LegEff8Data=0.875;
       else if(Mu1Pt<20.)   Mu1LegEff8Data=0.909;
       else if(Mu1Pt<30.)   Mu1LegEff8Data=0.922;
       else if(Mu1Pt<40.)   Mu1LegEff8Data=0.927;
       else if(Mu1Pt<50.)   Mu1LegEff8Data=0.925;
       else if(Mu1Pt<60.)   Mu1LegEff8Data=0.932;
       else if(Mu1Pt<80.)   Mu1LegEff8Data=0.924;
       else                 Mu1LegEff8Data=0.915;

       if     (Mu2Pt<15.)   Mu2LegEff8Data=0.875;
       else if(Mu2Pt<20.)   Mu2LegEff8Data=0.909;
       else if(Mu2Pt<30.)   Mu2LegEff8Data=0.922;
       else if(Mu2Pt<40.)   Mu2LegEff8Data=0.927;
       else if(Mu2Pt<50.)   Mu2LegEff8Data=0.925;
       else if(Mu2Pt<60.)   Mu2LegEff8Data=0.932;
       else if(Mu2Pt<80.)   Mu2LegEff8Data=0.924;
       else                 Mu2LegEff8Data=0.915;


       //Expectation - Mu12Ele23Dz
       float Mu1LegEff12Data=1., Mu2LegEff12Data=1.;
       if     (Mu1Pt<12.) Mu1LegEff12Data=0.028;
       else if(Mu1Pt<15.) Mu1LegEff12Data=0.88;
       else if(Mu1Pt<20.) Mu1LegEff12Data=0.914;
       else if(Mu1Pt<50.) Mu1LegEff12Data=0.93;
       else if(Mu1Pt<60.) Mu1LegEff12Data=0.938;
       else               Mu1LegEff12Data=0.93;

       if     (Mu2Pt<12.) Mu2LegEff12Data=0.028;
       else if(Mu2Pt<15.) Mu2LegEff12Data=0.88;
       else if(Mu2Pt<20.) Mu2LegEff12Data=0.914;
       else if(Mu2Pt<50.) Mu2LegEff12Data=0.93;
       else if(Mu2Pt<60.) Mu2LegEff12Data=0.938;
       else               Mu2LegEff12Data=0.93;


       //MC eff
       float TrigEffMu8Ele23MC=1., TrigEffMu12Ele23DzMC=1.;
       float ProbMu1_NoDz=Mu1LegEff8,        ProbMu2_NoDz=Mu2LegEff8;
       float ProbMu1_Dz  =Mu1LegEff8*DzEff,  ProbMu2_Dz  =Mu2LegEff12*DzEff;
       TrigEffMu8Ele23MC    = EleLegEff*(1-(1-ProbMu1_NoDz)*(1-ProbMu2_NoDz));
       TrigEffMu12Ele23DzMC = EleLegEff*(1-(1-ProbMu1_Dz)*(1-ProbMu2_Dz));

       //DataEff
       float TrigEffMu8Ele23Data=1., TrigEffMu12Ele23DzData=1.;
       ProbMu1_NoDz=Mu1LegEff8Data,            ProbMu2_NoDz=Mu2LegEff8Data;
       ProbMu1_Dz  =Mu1LegEff12Data*DzEffData, ProbMu2_Dz  =Mu2LegEff12Data*DzEffData;
       TrigEffMu8Ele23Data    = EleLegEffData*(1-(1-ProbMu1_NoDz)*(1-ProbMu2_NoDz));
       TrigEffMu12Ele23DzData = EleLegEffData*(1-(1-ProbMu1_Dz)*(1-ProbMu2_Dz));


       FillHist("MCClosure_Tot", 0., weight, 0., 2., 2);
       FillHist("MCClosure_Tot", 1., weight, 0., 2., 2);
       FillHist("MCClosure_ExpObs_Mu8Ele23", 0., weight*TrigEffMu8Ele23MC, 0., 2., 2);
       FillHist("MCClosure_ExpObs_Mu12Ele23Dz", 0., weight*TrigEffMu12Ele23DzMC, 0., 2., 2);
       FillHist("MCClosure_ExpMC_MuE_MuEdz", 0., weight*TrigEffMu8Ele23MC, 0., 2., 2);
       FillHist("MCClosure_ExpMC_MuE_MuEdz", 1., weight*TrigEffMu12Ele23DzMC, 0., 2., 2);
       FillHist("MCClosure_ExpData_MuE_MuEdz", 0., weight*TrigEffMu8Ele23Data, 0., 2., 2);
       FillHist("MCClosure_ExpData_MuE_MuEdz", 1., weight*TrigEffMu12Ele23DzData, 0., 2., 2);

       if(  PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ){
         FillHist("MCClosure_ExpObs_Mu8Ele23", 1., weight, 0., 2., 2);
       }
       if( PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ){
         FillHist("MCClosure_ExpObs_Mu12Ele23Dz", 1., weight, 0., 2., 2);
       }


     }//End of MCClosure(TriLep)


   }//End of EMu Trigger SF



/////////////////////////////////////////////////////////////////////////////////// 


return;
}// End of execute event loop
  


void May2017_TmpIDTrigSF::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void May2017_TmpIDTrigSF::BeginCycle() throw( LQError ){
  
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

May2017_TmpIDTrigSF::~May2017_TmpIDTrigSF() {
  
  Message("In May2017_TmpIDTrigSF Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void May2017_TmpIDTrigSF::FillCutFlow(TString cut, float weight){
  
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
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(6,"BeforeTrigMatch");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(7,"AfterTrigMatch");
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
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(5,"PreSel");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(6,"BeforeTrigMatch");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(7,"AfterTrigMatch");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(8,"2j");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(9,"3j");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(10,"ZVeto");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(11,"Nlj2");
    }
  }
}



void May2017_TmpIDTrigSF::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void May2017_TmpIDTrigSF::MakeHistograms(){
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
  // *  Remove//Overide this May2017_TmpIDTrigSFCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void May2017_TmpIDTrigSF::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
