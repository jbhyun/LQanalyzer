// $Id: Oct2017_EleFRSyst.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQOct2017_EleFRSyst Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Oct2017_EleFRSyst.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Oct2017_EleFRSyst);

 Oct2017_EleFRSyst::Oct2017_EleFRSyst() : AnalyzerCore(), out_muons(0) {

   SetLogName("Oct2017_EleFRSyst");
   Message("In Oct2017_EleFRSyst constructor", INFO);
   InitialiseAnalysis();
 }


 void Oct2017_EleFRSyst::InitialiseAnalysis() throw( LQError ) {
   
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

void Oct2017_EleFRSyst::ExecuteEvents()throw( LQError ){

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
   float pileup_reweight=1., pileup_reweight_systdown=1., pileup_reweight_systup=1.;
   if(!k_isdata){ pileup_reweight          = eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central); 
                  pileup_reweight_systup   = pileup_reweight>0.? eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::up)  /pileup_reweight :0. ;
                  pileup_reweight_systdown = pileup_reweight>0.? eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::down)/pileup_reweight :0. ;
   }

   //Numbet of Vertex PUreweight
   FillHist("Nvtx_nocut_PURW", eventbase->GetEvent().nVertices(), weight*pileup_reweight, 0., 50., 50);


   //Total Event(MCweight+PUrw)////////////////
   FillCutFlow("NoCut", weight*pileup_reweight);


   bool NormCheck=false, FRMeasure=false, SiglWP=false, FRScan=false, OptMETMTWCuts=false, Closure=false, IDValidation=false;
   bool UnPreTrig=false, PreTrig=false, SiglPreTrig=false, MultPreTrig=false, SystRun=false;
   bool HighdXYFRCheck=false, IsoIPOpt=false;
   bool EleFRSyst=false, NvtxProf=false;
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("NormCheck"))      NormCheck      = true;
     else if(k_flags.at(i).Contains("FRMeasure"))      FRMeasure      = true;
     else if(k_flags.at(i).Contains("SiglWP"))         SiglWP         = true;
     else if(k_flags.at(i).Contains("FRScan"))         FRScan         = true;
     else if(k_flags.at(i).Contains("UnPreTrig"))      UnPreTrig      = true;
     else if(k_flags.at(i).Contains("SiglPreTrig"))    SiglPreTrig    = true;
     else if(k_flags.at(i).Contains("MultPreTrig"))    MultPreTrig    = true;
     else if(k_flags.at(i).Contains("OptMETMTWCuts"))  OptMETMTWCuts  = true;
     else if(k_flags.at(i).Contains("Closure"))        Closure        = true;
     else if(k_flags.at(i).Contains("HighdXYFRCheck")) HighdXYFRCheck = true;
     else if(k_flags.at(i).Contains("IsoIPOpt"))       IsoIPOpt       = true;
     else if(k_flags.at(i).Contains("IDValidation"))   IDValidation   = true;
     else if(k_flags.at(i).Contains("EleFRSyst"))      EleFRSyst      = true;
     else if(k_flags.at(i).Contains("NvtxProf"))       NvtxProf       = true;
     else if(k_flags.at(i).Contains("SystRun"))        SystRun        = true;
   }

   if(SiglPreTrig || MultPreTrig) PreTrig=true;
    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //ListTriggersAvailable();

   //Trigger Path of Analysis
   bool Pass_Trigger=false;
   float trigger_ps_weight=1., trigger_period_weight=1.;
   if(FRMeasure || NormCheck || OptMETMTWCuts || HighdXYFRCheck || EleFRSyst){
     if(UnPreTrig){
       if( PassTrigger("HLT_Ele27_WPTight_Gsf_v") ) Pass_Trigger=true;
       if(!isData) trigger_ps_weight=WeightByTrigger("HLT_Ele27_WPTight_Gsf_v", TargetLumi);
     }
     if(SiglPreTrig){
       if( PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ) Pass_Trigger=true;
       if(!isData) trigger_ps_weight=WeightByTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v", TargetLumi);
     }
     if(MultPreTrig){
         eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
         eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
         eventbase->GetElectronSel()->SetBETrRegIncl(false);
         eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
         eventbase->GetElectronSel()->SetdxySigMax(3.);
         eventbase->GetElectronSel()->SetApplyConvVeto(true);
       std::vector<snu::KElectron> electronFakeLColl; eventbase->GetElectronSel()->Selection(electronFakeLColl);

       int Case=0;
       if( electronFakeLColl.size()==1 && electronFakeLColl.at(0).Pt()>25
           && PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ) {Pass_Trigger=true; Case=1;}
       if( electronFakeLColl.size()==1 && electronFakeLColl.at(0).Pt()<25 && electronFakeLColl.at(0).Pt()>20 
           && PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ) {Pass_Trigger=true; Case=2;}
       if( electronFakeLColl.size()==1 && electronFakeLColl.at(0).Pt()<20 && electronFakeLColl.at(0).Pt()>15 
           && PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ) {Pass_Trigger=true; Case=3;}

       if     (!isData && Case==1) trigger_ps_weight=WeightByTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v", TargetLumi);
       else if(!isData && Case==2) trigger_ps_weight=WeightByTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v", TargetLumi);
       else if(!isData && Case==3) trigger_ps_weight=WeightByTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v", TargetLumi);
     }

   }
   if(Closure){
     int Pass_Trigger1=0, Pass_Trigger2=0;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ) Pass_Trigger1++;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger2++;

     if(isData){
       int DataPeriod=GetDataPeriod();
       if( DataPeriod>0 && DataPeriod<7 && Pass_Trigger1==1 ) Pass_Trigger=true;
       else if( DataPeriod==7 && Pass_Trigger2==1 ) Pass_Trigger=true;
     }
     else{
       if( Pass_Trigger1>0 || Pass_Trigger2>0 ) Pass_Trigger=true;
       trigger_period_weight=(Pass_Trigger1*27.257618+Pass_Trigger2*8.605696)/35.863314;
       trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
     }
   }
   if(IsoIPOpt || IDValidation){
     if( PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger=true;
     if(!isData) trigger_ps_weight=WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);
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
   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
  
   /**PreSelCut***********************************************************************************************/
   //Intended for Code speed boosting up.
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, true);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
     if(FRMeasure || NormCheck || HighdXYFRCheck || EleFRSyst ){ if( !(electronPreColl.size()>=1) ) return; }
     else if(Closure)                   { if( !(electronPreColl.size()>=1 && muonPreColl.size()>=2) ) return; }
     else if(IsoIPOpt)                  { if( !(electronPreColl.size()>=2) ) return; }
     FillCutFlow("PreSel", weight);
   /**********************************************************************************************************/

   std::vector<snu::KMuon> muonLooseColl, muonTightColl;
     for(int i=0; i<muonPreColl.size(); i++){
       if(PassIDCriteria(muonPreColl.at(i),"Test_POGLIsop4IPp5p1Chi100"))      muonLooseColl.push_back(muonPreColl.at(i));
       if(PassIDCriteria(muonPreColl.at(i),"Test_POGTIsop20IPp01p05sig4Chi4")) muonTightColl.push_back(muonPreColl.at(i));
     }
   std::vector<snu::KMuon> muonColl;      if(k_running_nonprompt){ muonColl=muonLooseColl;} else{ muonColl=muonTightColl;}


     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025);   eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronFakeLPreOptColl; eventbase->GetElectronSel()->Selection(electronFakeLPreOptColl);


     std::vector<snu::KElectron> electronFakeLColl, electronTightColl;
       for(int i=0; i<electronPreColl.size(); i++){
         if(PassIDCriteria(electronPreColl.at(i), "LMVA06Isop4IPp025p05sig4")) electronFakeLColl.push_back(electronPreColl.at(i));
       //Tight
         if(PassIDCriteria(electronPreColl.at(i), "POGWP90Isop06IPp025p05sig4")) electronTightColl.push_back(electronPreColl.at(i));
       }

//     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
//     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
//     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
//     eventbase->GetElectronSel()->SetBETrRegIncl(false);
//     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
//     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
//     eventbase->GetElectronSel()->SetdxySigMax(4.);
//     eventbase->GetElectronSel()->SetApplyConvVeto(true);
//   std::vector<snu::KElectron> electronTightColl; eventbase->GetElectronSel()->Selection(electronTightColl);
   std::vector<snu::KElectron> electronColl;
     if(!k_running_nonprompt){ electronColl=electronTightColl;} else{electronColl=electronFakeLColl;}

     bool LeptonVeto=false;
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
   std::vector<snu::KJet> jetNoVetoColl; eventbase->GetJetSel()->Selection(jetNoVetoColl, LeptonVeto, muonLooseColl, electronFakeLColl);
   std::vector<snu::KJet> bjetNoVetoColl = SelBJets(jetNoVetoColl, "Medium");

   std::vector<snu::KJet> jetColl  = SkimJetColl(jetNoVetoColl,  electronFakeLColl, muonLooseColl, "EleMuVeto");
   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");

   float met           = eventbase->GetEvent().MET(snu::KEvent::pfmet);
   float metphi = eventbase->GetEvent().METPhi();
   double met_x = eventbase->GetEvent().PFMETx();
   double met_y = eventbase->GetEvent().PFMETy();
   double Pzv,Pzv1, Pzv2;
   snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));
   //snu::KParticle v[4]; v[0].SetPx(met_x); v[0].SetPy(met_y);
   //                     v[1].SetPx(met_x); v[1].SetPy(met_y);

   int nbjets=bjetColl.size(); int njets=jetColl.size();
   int Nvtx=eventbase->GetEvent().nVertices();


   /*****************************************************
   **Scale Factors
   *****************************************************/
   float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1., btag_sf=1.;
   float trigger_sf=1.;
   float fake_weight=1.; bool EventCand=false;
   float geneff_weight=1., gennorm_weight=1., k_factor_weight=1.;
   float nvtx_reweight=1.; 


   //This part is for boosting up speed.. SF part takes rather longer time than expected
   //Systematic study case should be treated in other part, since there are so many diff. kinds of variations affecting each other.
   if(!SystRun){
     if     (Closure     ){ if(electronFakeLColl.size()>=1 && muonLooseColl.size()>=2) EventCand=true; }
     else if(IDValidation){ if(electronFakeLColl.size()>=2) EventCand=true;}
     else if(NormCheck || OptMETMTWCuts || FRMeasure ){ if(electronFakeLColl.size()>=1) EventCand=true;}

     if(!isData){
       if(EventCand && (NormCheck || Closure || OptMETMTWCuts || IDValidation || FRMeasure || EleFRSyst)){
         //geneff_weight   = GenFilterEfficiency(k_sample_name);
         //gennorm_weight  = SignalNorm(k_sample_name, 200.);
         k_factor_weight = GetKFactor();
  
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electronColl);
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
    
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muonColl);
         //id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muonColl);
         iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
         trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);

         if(UnPreTrig) nvtx_reweight = mcdata_correction->UserPileupWeight(eventbase->GetEvent(), jetColl.size());
         if(PreTrig)   nvtx_reweight = GetPreTrigPURW(Nvtx);

         if(Closure){
           btag_sf       = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);}
       }
       //trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );
       if(HighdXYFRCheck && (k_sample_name.Contains("qcd") ||k_sample_name.Contains("QCD"))) pileup_reweight=1.;
     }
     else if(k_running_nonprompt && EventCand && Closure ){
       //fake_weight = GetFakeWeight(muonLooseColl, electronFakeLColl, "HNTrilepFakeL2", "HNTrilepTight2", "HctoWAFakeLoose", "POGMVAMIP", "");
       fake_weight = GetFakeWeight(muonLooseColl, electronFakeLColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA06Isop4IPp025p05sig4", "POGWP90Isop06IPp025p05sig4", "TrkIsoVVLConeSUSY");
       //fake_weight = GetFakeWeight(muonLooseColl, electronFakeLColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA06v1Isop4IPp025p05sig4", "POGWP90Isop06IPp025p05sig4", "TrkIsoVVLConeSUSY");
       //fake_weight = GetFakeWeight(muonLooseColl, electronFakeLColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA06v2Isop4IPp025p05sig4", "POGWP90Isop06IPp025p05sig4", "TrkIsoVVLConeSUSY");
       //fake_weight = GetFakeWeight(muonLooseColl, electronFakeLColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA08v1Isop4IPp025p05sig4", "POGWP90Isop08IPp025p05sig4", "TrkIsoVVLConeSUSY");
       //fake_weight = GetFakeWeight(muonLooseColl, electronFakeLColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA08v2Isop4IPp025p05sig4", "POGWP90Isop08IPp025p05sig4", "TrkIsoVVLConeSUSY");
       //fake_weight = GetFakeWeight(muonLooseColl, electronFakeLColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA1v1Isop4IPp025p05sig4", "POGWP90Isop1IPp025p05sig4", "TrkIsoVVLConeSUSY");
       //fake_weight = GetFakeWeight(muonLooseColl, electronFakeLColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA1v2Isop4IPp025p05sig4", "POGWP90Isop1IPp025p05sig4", "TrkIsoVVLConeSUSY");


     }
   }

   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight*fake_weight*btag_sf;
   weight *= nvtx_reweight;
   weight *= geneff_weight*gennorm_weight*k_factor_weight;
   //-----------------------------------------------------------------------------------------//


//----------------------------------------------------------------------------------//
//////////////////////////////////////////////////////////////////////////////////////
/////// MAIN ANALYSIS CODE ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//----------------------------------------------------------------------------------//

   if(EleFRSyst){
     // FR Syst.
     //1. Prompt Uncertainty
     //   a) Norm Unc. ~10% conservatively; incl Trig+ID+Norm Unc.
     //   b) Pile-Up + Prescale Unc ; Nvtx Variation
     //2. Impact from cut value : MET, MTW cut variation
     //3. Balance Jet Systematic : dR, jet Pt. Njet
     //4. Variance when B-jet present
     //5. Prompt Rate: What happens if prompt is included 

     //How
     //1.   By weight variation
     //2-4. Better with in code.
     //5.   Measure Prompt Rate in DY events, and apply it.


     const int NPtEdges=6;
     float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100.};

     CheckEleFRSyst(electronColl, electronPreColl, muonColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight, NPtEdges, PtEdges, "LMVA06Isop4IPp025p05sig4", "POGWP90Isop06IPp025p05sig4", "", "SiglPreTrig");

   }
   if(NvtxProf){

     if(isData && PreTrig){
       CheckNvtxProfile(electronColl, electronFakeLColl, muonLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight, "PreTrigData", "");
     }
     if(isData && UnPreTrig){
       CheckNvtxProfile(electronColl, electronFakeLColl, muonLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight, "UnPreTrigData", "");
     }
     if(!isData && SiglPreTrig){
       CheckNvtxProfile(electronColl, electronFakeLColl, muonLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight, "PreTrigMCPUcent", "");
     }
     if(!isData && UnPreTrig){
       CheckNvtxProfile(electronColl, electronFakeLColl, muonLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight, "UnPreTrigMCPUcent", "");
     }

   }
   if(NormCheck){
     
     if(!SystRun){
       if(PreTrig){
        float nvtx_reweight_UnPreTrig=1., nvtx_reweight_PreTrig=1.;
        if(!isData){
           nvtx_reweight_UnPreTrig = mcdata_correction->UserPileupWeight(eventbase->GetEvent(), jetColl.size()); //Measured in DY UnPreTrig
           nvtx_reweight_PreTrig = GetPreTrigPURW(Nvtx);
        }

        CheckNormCR(electronColl, electronFakeLColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "PreTrig");
        CheckNormCR(electronColl, electronFakeLColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight*nvtx_reweight_UnPreTrig, "_NvtxRWUnPreTrig", "PreTrig");
        CheckNormCR(electronColl, electronFakeLColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight*nvtx_reweight_PreTrig, "_NvtxRWPreTrig", "PreTrig");


       }
       if(UnPreTrig){ CheckNormCR(electronColl, electronFakeLColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "UnPreTrig");
         
       }                
     }//End of not SystRun
     if(SystRun){

       if(isData){
         if(PreTrig){
           DoSystRun("NormCheck", "PreTrig",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     1.);
         }
         if(UnPreTrig){
           DoSystRun("NormCheck", "UnPreTrig",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     1.);
         }
       }
       else{
           eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
           eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
           eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
           eventbase->GetElectronSel()->SetBETrRegIncl(false);
           eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
           eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
           eventbase->GetElectronSel()->SetdxySigMax(4.);
           eventbase->GetElectronSel()->SetApplyConvVeto(true);
         std::vector<snu::KElectron> EleTElEnUpColl; eventbase->GetElectronSel()->Selection(EleTElEnUpColl, "SystUpElEn");
  
           eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
           eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
           eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
           eventbase->GetElectronSel()->SetBETrRegIncl(false);
           eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
           eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025);   eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
           eventbase->GetElectronSel()->SetdxySigMax(4.);
           eventbase->GetElectronSel()->SetApplyConvVeto(true);
         std::vector<snu::KElectron> EleTElEnDownColl; eventbase->GetElectronSel()->Selection(EleTElEnDownColl, "SystDownElEn");
  
  
           eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
           eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
         std::vector<snu::KJet> jetJESUpColl; eventbase->GetJetSel()->Selection(jetJESUpColl, LeptonVeto, muonLooseColl, electronFakeLColl, "SystUpJES");
           eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
           eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
         std::vector<snu::KJet> jetJERUpColl; eventbase->GetJetSel()->Selection(jetJERUpColl, LeptonVeto, muonLooseColl, electronFakeLColl, "SystUpJER");
           eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
           eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
         std::vector<snu::KJet> jetJESDownColl; eventbase->GetJetSel()->Selection(jetJESDownColl, LeptonVeto, muonLooseColl, electronFakeLColl, "SystDownJES");
           eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
           eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
         std::vector<snu::KJet> jetJERDownColl; eventbase->GetJetSel()->Selection(jetJERDownColl, LeptonVeto, muonLooseColl, electronFakeLColl, "SystDownJER");
  
         float met_JESup     = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetEn       , snu::KEvent::up);
         float met_JERup     = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetRes      , snu::KEvent::up);
         float met_Unclup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::Unclustered , snu::KEvent::up);
         float met_ElEnup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::ElectronEn  , snu::KEvent::up);
         float met_MuEnup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::MuonEn      , snu::KEvent::up);
         float met_JESdown   = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetEn       , snu::KEvent::down);
         float met_JERdown   = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetRes      , snu::KEvent::down);
         float met_Uncldown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::Unclustered , snu::KEvent::down);
         float met_ElEndown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::ElectronEn  , snu::KEvent::down);
         float met_MuEndown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::MuonEn      , snu::KEvent::down);
  
  
         float id_weight_ele_ElEnup=1.  , reco_weight_ele_ElEnup=1.  ;
         float id_weight_ele_ElEndown=1., reco_weight_ele_ElEndown=1.;
  
         if(!isData){
           id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electronColl);
           reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
  
           id_weight_ele_ElEnup   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", EleTElEnUpColl);
           reco_weight_ele_ElEnup = mcdata_correction->ElectronRecoScaleFactor(EleTElEnUpColl);
  
           id_weight_ele_ElEndown   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", EleTElEnDownColl);
           reco_weight_ele_ElEndown = mcdata_correction->ElectronRecoScaleFactor(EleTElEnDownColl);

           if(UnPreTrig) nvtx_reweight = mcdata_correction->UserPileupWeight(eventbase->GetEvent(), jetColl.size()); //Measured in DY UnPreTrig
           if(PreTrig) nvtx_reweight = GetPreTrigPURW(Nvtx);
         }
         float systweight_central  =weight*nvtx_reweight*id_weight_ele       *reco_weight_ele       ;

         float systweight_Nvtxup   =weight              *id_weight_ele       *reco_weight_ele       ;
         float systweight_ElEnup   =weight*nvtx_reweight*id_weight_ele_ElEnup*reco_weight_ele_ElEnup;
         float systweight_PUup     =weight*nvtx_reweight*id_weight_ele       *reco_weight_ele       *pileup_reweight_systup;
  
         float systweight_Nvtxdown =weight              *id_weight_ele         *reco_weight_ele         ;
         float systweight_ElEndown =weight*nvtx_reweight*id_weight_ele_ElEndown*reco_weight_ele_ElEndown;
         float systweight_PUdown   =weight*nvtx_reweight*id_weight_ele         *reco_weight_ele         *pileup_reweight_systdown;


         if(PreTrig){
           //Central
           DoSystRun("NormCheck", "PreTrig",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     systweight_central);

           //Variation Up
           DoSystRun("NormCheck", "PreTrigSystUpPU",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     systweight_PUup);
           DoSystRun("NormCheck", "PreTrigSystUpNvtx",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     systweight_Nvtxup);
           DoSystRun("NormCheck", "PreTrigSystUpElEn",
                     EleTElEnUpColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi),
                     systweight_ElEnup);
           DoSystRun("NormCheck", "PreTrigSystUpMuEn",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi),
                     systweight_central);//No Muon is used, so fine to put anything there.
           DoSystRun("NormCheck", "PreTrigSystUpUncl",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi),
                     systweight_central);
           DoSystRun("NormCheck", "PreTrigSystUpJES",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetJESUpColl, bjetColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi),
                     systweight_central);
           DoSystRun("NormCheck", "PreTrigSystUpJER",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetJERUpColl, bjetColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi),
                     systweight_central);
           //Variation Down
           DoSystRun("NormCheck", "PreTrigSystDownPU",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     systweight_PUdown);
           DoSystRun("NormCheck", "PreTrigSystDownNvtx",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     systweight_Nvtxdown);
           DoSystRun("NormCheck", "PreTrigSystDownElEn",
                     EleTElEnDownColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi),
                     systweight_ElEndown);
           DoSystRun("NormCheck", "PreTrigSystDownMuEn",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi),
                     systweight_central);//No Muon is used, so fine to put anything there.
           DoSystRun("NormCheck", "PreTrigSystDownUncl",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi),
                     systweight_central);
           DoSystRun("NormCheck", "PreTrigSystDownJES",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetJESDownColl, bjetColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi),
                     systweight_central);
           DoSystRun("NormCheck", "PreTrigSystDownJER",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetJERDownColl, bjetColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi),
                     systweight_central);
         }//End of PreTrig
         if(UnPreTrig){
           //Central
           DoSystRun("NormCheck", "UnPreTrig",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     systweight_central);

           //Variation Up
           DoSystRun("NormCheck", "UnPreTrigSystUpPU",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     systweight_PUup);
           DoSystRun("NormCheck", "UnPreTrigSystUpElEn",
                     EleTElEnUpColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi),
                     systweight_ElEnup);
           DoSystRun("NormCheck", "UnPreTrigSystUpMuEn",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi),
                     systweight_central);//No Muon is used, so fine to put anything there.
           DoSystRun("NormCheck", "UnPreTrigSystUpUncl",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi),
                     systweight_central);
           DoSystRun("NormCheck", "UnPreTrigSystUpJES",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetJESUpColl, bjetColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi),
                     systweight_central);
           DoSystRun("NormCheck", "UnPreTrigSystUpJER",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetJERUpColl, bjetColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi),
                     systweight_central);
           //Variation Down
           DoSystRun("NormCheck", "UnPreTrigSystDownPU",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     systweight_PUdown);
           DoSystRun("NormCheck", "UnPreTrigSystDownElEn",
                     EleTElEnDownColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi),
                     systweight_ElEndown);
           DoSystRun("NormCheck", "UnPreTrigSystDownMuEn",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi),
                     systweight_central);//No Muon is used, so fine to put anything there.
           DoSystRun("NormCheck", "UnPreTrigSystDownUncl",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi),
                     systweight_central);
           DoSystRun("NormCheck", "UnPreTrigSystDownJES",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetJESDownColl, bjetColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi),
                     systweight_central);
           DoSystRun("NormCheck", "UnPreTrigSystDownJER",
                     electronColl, electronFakeLColl, muonColl, muonLooseColl, jetJERDownColl, bjetColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi),
                     systweight_central);
         }//End of UnPreTrig
       }//End of !isData
     }//End of SystRun
     
   }//End of NormCheck
   if(FRMeasure){

     const int NPtEdges=6;
     float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100.};

     if(!SystRun){
       if(SiglWP){
         MeasureFakeRate(electronPreColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight, NPtEdges, PtEdges,
              "LMVA06Isop4IPp025p05sig4", "POGWP90Isop06IPp025p05sig4", "_WP90Isop06_L928676Isop4", "SiglPreTrig");
  
       }
       if(FRScan){
  
         //Full Scan
         const int NIsoCuts=7, NMVACuts=21;
         float IsoCuts[NIsoCuts]={0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
         float MVACuts[NMVACuts]={-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
    
         ScanFakeRate(electronFakeLPreOptColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met_x, met_y, weight, NMVACuts, MVACuts, NIsoCuts, IsoCuts, NPtEdges, PtEdges, "SiglPreTrig2D1D");
  
       }
     }
     if(SiglWP && SystRun){
         eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
         eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
         eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
         eventbase->GetElectronSel()->SetBETrRegIncl(false);
         eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
         eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
         eventbase->GetElectronSel()->SetdxySigMax(4.);
         eventbase->GetElectronSel()->SetApplyConvVeto(true);
       std::vector<snu::KElectron> EleTElEnUpColl; eventbase->GetElectronSel()->Selection(EleTElEnUpColl, "SystUpElEn");
  
         eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
         eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
         eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
         eventbase->GetElectronSel()->SetBETrRegIncl(false);
         eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
         eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025);   eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
         eventbase->GetElectronSel()->SetdxySigMax(4.);
         eventbase->GetElectronSel()->SetApplyConvVeto(true);
       std::vector<snu::KElectron> EleTElEnDownColl; eventbase->GetElectronSel()->Selection(EleTElEnDownColl, "SystDownElEn");
  
  
       LeptonVeto=false;
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJESUpColl; eventbase->GetJetSel()->Selection(jetJESUpColl, LeptonVeto, muonLooseColl, electronFakeLColl, "SystUpJES");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJERUpColl; eventbase->GetJetSel()->Selection(jetJERUpColl, LeptonVeto, muonLooseColl, electronFakeLColl, "SystUpJER");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJESDownColl; eventbase->GetJetSel()->Selection(jetJESDownColl, LeptonVeto, muonLooseColl, electronFakeLColl, "SystDownJES");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJERDownColl; eventbase->GetJetSel()->Selection(jetJERDownColl, LeptonVeto, muonLooseColl, electronFakeLColl, "SystDownJER");
  
       float met_JESup     = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetEn       , snu::KEvent::up);
       float met_JERup     = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetRes      , snu::KEvent::up);
       float met_Unclup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::Unclustered , snu::KEvent::up);
       float met_ElEnup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::ElectronEn  , snu::KEvent::up);
       float met_MuEnup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::MuonEn      , snu::KEvent::up);
       float met_JESdown   = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetEn       , snu::KEvent::down);
       float met_JERdown   = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetRes      , snu::KEvent::down);
       float met_Uncldown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::Unclustered , snu::KEvent::down);
       float met_ElEndown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::ElectronEn  , snu::KEvent::down);
       float met_MuEndown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::MuonEn      , snu::KEvent::down);
  
  
       float id_weight_ele_ElEnup=1.  , reco_weight_ele_ElEnup=1.  ;
       float id_weight_ele_ElEndown=1., reco_weight_ele_ElEndown=1.;
  
       if(!isData){
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electronColl);
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
  
         id_weight_ele_ElEnup   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", EleTElEnUpColl);
         reco_weight_ele_ElEnup = mcdata_correction->ElectronRecoScaleFactor(EleTElEnUpColl);
  
         id_weight_ele_ElEndown   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", EleTElEnDownColl);
         reco_weight_ele_ElEndown = mcdata_correction->ElectronRecoScaleFactor(EleTElEnDownColl);

         if(PreTrig) nvtx_reweight = GetPreTrigPURW(Nvtx);
       }
       float systweight_central  =weight*nvtx_reweight*id_weight_ele       *reco_weight_ele       ;

       float systweight_Nvtxup   =weight              *id_weight_ele       *reco_weight_ele       ;
       float systweight_ElEnup   =weight*nvtx_reweight*id_weight_ele_ElEnup*reco_weight_ele_ElEnup;
       float systweight_PUup     =weight*nvtx_reweight*id_weight_ele       *reco_weight_ele       *pileup_reweight_systup;
  
       float systweight_Nvtxdown =weight              *id_weight_ele         *reco_weight_ele         ;
       float systweight_ElEndown =weight*nvtx_reweight*id_weight_ele_ElEndown*reco_weight_ele_ElEndown;
       float systweight_PUdown   =weight*nvtx_reweight*id_weight_ele         *reco_weight_ele         *pileup_reweight_systdown;


       if(PreTrig){
         //Central
         DoSystRun("PromptVarCheck", "PreTrig",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_central);

         //Variation Up
         DoSystRun("PromptVarCheck", "PreTrigSystUpPU",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_PUup);
         DoSystRun("PromptVarCheck", "PreTrigSystUpNvtx",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_Nvtxup);
         DoSystRun("PromptVarCheck", "PreTrigSystUpElEn",
                   EleTElEnUpColl, electronFakeLColl, muonColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi),
                   systweight_ElEnup);
         DoSystRun("PromptVarCheck", "PreTrigSystUpMuEn",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi),
                   systweight_central);//No Muon is used, so fine to put anything there.
         DoSystRun("PromptVarCheck", "PreTrigSystUpUncl",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi),
                   systweight_central);
         DoSystRun("PromptVarCheck", "PreTrigSystUpJES",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetJESUpColl, bjetNoVetoColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi),
                   systweight_central);
         DoSystRun("PromptVarCheck", "PreTrigSystUpJER",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetJERUpColl, bjetNoVetoColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi),
                   systweight_central);
         //Variation Down
         DoSystRun("PromptVarCheck", "PreTrigSystDownPU",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_PUdown);
         DoSystRun("PromptVarCheck", "PreTrigSystDownNvtx",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_Nvtxdown);
         DoSystRun("PromptVarCheck", "PreTrigSystDownElEn",
                   EleTElEnDownColl, electronFakeLColl, muonColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi),
                   systweight_ElEndown);
         DoSystRun("PromptVarCheck", "PreTrigSystDownMuEn",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi),
                   systweight_central);//No Muon is used, so fine to put anything there.
         DoSystRun("PromptVarCheck", "PreTrigSystDownUncl",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi),
                   systweight_central);
         DoSystRun("PromptVarCheck", "PreTrigSystDownJES",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetJESDownColl, bjetNoVetoColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi),
                   systweight_central);
         DoSystRun("PromptVarCheck", "PreTrigSystDownJER",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetJERDownColl, bjetNoVetoColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi),
                   systweight_central);
       }//End of PreTrig

     }
   }
   if(Closure){

     if(!SystRun){

       if(k_sample_name.Contains("TT_powheg")){ if(NLeptonicBosonDecay(truthColl)==2) return; }
       CheckTrilepCRs(electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "");
     }//End of Not SystRun
     if(SystRun){

       // Syst Sel. and Syst Corr. ----------------------------------------------------------------------------------------------------//
         eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
         eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
         eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
         eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
         eventbase->GetMuonSel()->SetdxySigMax(4.);
         eventbase->GetMuonSel()->SetChiNdof(4.);
       std::vector<snu::KMuon> MuTMuEnUpColl; eventbase->GetMuonSel()->Selection(MuTMuEnUpColl,true, "SystUpMuEn");

         eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
         eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
         eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
         eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
         eventbase->GetMuonSel()->SetdxySigMax(4.);
         eventbase->GetMuonSel()->SetChiNdof(4.);
       std::vector<snu::KMuon> MuTMuEnDownColl; eventbase->GetMuonSel()->Selection(MuTMuEnDownColl,true, "SystDownMuEn");

         eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
         eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
         eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
         eventbase->GetElectronSel()->SetBETrRegIncl(false);
         eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
         eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
         eventbase->GetElectronSel()->SetdxySigMax(4.);
         eventbase->GetElectronSel()->SetApplyConvVeto(true);
       std::vector<snu::KElectron> EleTElEnUpColl; eventbase->GetElectronSel()->Selection(EleTElEnUpColl, "SystUpElEn");

         eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
         eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
         eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
         eventbase->GetElectronSel()->SetBETrRegIncl(false);
         eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
         eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025);   eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
         eventbase->GetElectronSel()->SetdxySigMax(4.);
         eventbase->GetElectronSel()->SetApplyConvVeto(true);
       std::vector<snu::KElectron> EleTElEnDownColl; eventbase->GetElectronSel()->Selection(EleTElEnDownColl, "SystDownElEn");

       LeptonVeto=true;
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJESUpColl; eventbase->GetJetSel()->Selection(jetJESUpColl, LeptonVeto, muonLooseColl, electronFakeLColl, "SystUpJES");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJERUpColl; eventbase->GetJetSel()->Selection(jetJERUpColl, LeptonVeto, muonLooseColl, electronFakeLColl, "SystUpJER");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJESDownColl; eventbase->GetJetSel()->Selection(jetJESDownColl, LeptonVeto, muonLooseColl, electronFakeLColl, "SystDownJES");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJERDownColl; eventbase->GetJetSel()->Selection(jetJERDownColl, LeptonVeto, muonLooseColl, electronFakeLColl, "SystDownJER");

       std::vector<snu::KJet> bjetJESUpColl   = SelBJets(jetJESUpColl  , "Medium");
       std::vector<snu::KJet> bjetJERUpColl   = SelBJets(jetJERUpColl  , "Medium");
       std::vector<snu::KJet> bjetJESDownColl = SelBJets(jetJESDownColl, "Medium");
       std::vector<snu::KJet> bjetJERDownColl = SelBJets(jetJERDownColl, "Medium");

       float met_JESup     = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetEn       , snu::KEvent::up);
       float met_JERup     = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetRes      , snu::KEvent::up);
       float met_Unclup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::Unclustered , snu::KEvent::up);
       float met_ElEnup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::ElectronEn  , snu::KEvent::up);
       float met_MuEnup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::MuonEn      , snu::KEvent::up);
       float met_JESdown   = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetEn       , snu::KEvent::down);
       float met_JERdown   = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetRes      , snu::KEvent::down);
       float met_Uncldown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::Unclustered , snu::KEvent::down);
       float met_ElEndown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::ElectronEn  , snu::KEvent::down);
       float met_MuEndown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::MuonEn      , snu::KEvent::down);

       //Scale Factors--------------------------------------------------------------------------------//      
       float systweight=weight; //Lumi weight+Prescale weight+Period weight+PU weight applied by default;
       float id_weight_ele_ElEnup=1.  , reco_weight_ele_ElEnup=1.  , id_weight_mu_MuEnup=1.  , trk_weight_mu_MuEnup=1.;
       float id_weight_ele_ElEndown=1., reco_weight_ele_ElEndown=1., id_weight_mu_MuEndown=1., trk_weight_mu_MuEndown=1.;
       float btag_sf_JESup=1.  , btag_sf_JERup=1.  , btag_sf_LTagup=1.  , btag_sf_BCTagup=1.;
       float btag_sf_JESdown=1., btag_sf_JERdown=1., btag_sf_LTagdown=1., btag_sf_BCTagdown=1.;
       float fake_weight_FRup=1., fake_weight_FRdown=1.;
       if(!isData){
         //geneff_weight   = GenFilterEfficiency(k_sample_name);
         //gennorm_weight  = SignalNorm(k_sample_name, 200.);
         k_factor_weight = GetKFactor();
    
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electronColl);
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muonColl)*mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
         trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
         btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

         id_weight_ele_ElEnup   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", EleTElEnUpColl);
         reco_weight_ele_ElEnup = mcdata_correction->ElectronRecoScaleFactor(EleTElEnUpColl);
         id_weight_mu_MuEnup    = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", MuTMuEnUpColl)*mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", MuTMuEnUpColl);
         trk_weight_mu_MuEnup   = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnUpColl);

         id_weight_ele_ElEndown   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", EleTElEnDownColl);
         reco_weight_ele_ElEndown = mcdata_correction->ElectronRecoScaleFactor(EleTElEnDownColl);
         id_weight_mu_MuEndown    = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", MuTMuEnDownColl)*mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", MuTMuEnDownColl);
         trk_weight_mu_MuEndown   = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnDownColl);

         btag_sf_LTagup = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpLTag");
         btag_sf_BCTagup= BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpBCTag");
         btag_sf_JESup  = BTagScaleFactor_1a(jetJESUpColl, snu::KJet::CSVv2, snu::KJet::Medium);
         btag_sf_JERup  = BTagScaleFactor_1a(jetJERUpColl, snu::KJet::CSVv2, snu::KJet::Medium);

         btag_sf_LTagdown = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownLTag");
         btag_sf_BCTagdown= BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownBCTag");
         btag_sf_JESdown  = BTagScaleFactor_1a(jetJESDownColl, snu::KJet::CSVv2, snu::KJet::Medium);
         btag_sf_JERdown  = BTagScaleFactor_1a(jetJERDownColl, snu::KJet::CSVv2, snu::KJet::Medium);
       }
       else if(k_running_nonprompt){
         fake_weight        = GetFakeWeight(muonLooseColl, electronFakeLColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA06Isop4IPp025p05sig4", "POGWP90Isop06IPp025p05sig4", "TrkIsoVVLConeSUSY");
         fake_weight_FRup   = GetFakeWeight(muonLooseColl, electronFakeLColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA06Isop4IPp025p05sig4", "POGWP90Isop06IPp025p05sig4", "TrkIsoVVLConeSUSYSystUp");
         fake_weight_FRdown = GetFakeWeight(muonLooseColl, electronFakeLColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA06Isop4IPp025p05sig4", "POGWP90Isop06IPp025p05sig4", "TrkIsoVVLConeSUSYSystDown");
         //fake_weight = GetFakeWeight(muonLooseColl, electronFakeLColl, "HNTrilepFakeL2", "HNTrilepTight2", "HctoWAFakeLoose", "POGMVAMIP", "");
       }
       float systweight_central=weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight;
       float systweight_ElEnup =weight*k_factor_weight*id_weight_ele_ElEnup*reco_weight_ele_ElEnup*id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight;
       float systweight_MuEnup =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu_MuEnup*trk_weight_mu_MuEnup*btag_sf        *fake_weight;
       float systweight_JESup  =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JESup  *fake_weight;
       float systweight_JERup  =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JERup  *fake_weight;
       float systweight_LTagup =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_LTagup *fake_weight;
       float systweight_BCTagup=weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_BCTagup*fake_weight;
       float systweight_PUup   =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight*pileup_reweight_systup;
       float systweight_FRup   =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight_FRup;

       float systweight_ElEndown =weight*k_factor_weight*id_weight_ele_ElEndown*reco_weight_ele_ElEndown*id_weight_mu         *trk_weight_mu         *btag_sf          *fake_weight;
       float systweight_MuEndown =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu_MuEndown*trk_weight_mu_MuEndown*btag_sf          *fake_weight;
       float systweight_JESdown  =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JESdown  *fake_weight;
       float systweight_JERdown  =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JERdown  *fake_weight;
       float systweight_LTagdown =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_LTagdown *fake_weight;
       float systweight_BCTagdown=weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_BCTagdown*fake_weight;
       float systweight_PUdown   =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *fake_weight*pileup_reweight_systdown;
       float systweight_FRdown   =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *fake_weight_FRdown;

       //----------------------------------------------------------------------------------------------------------------------//


       DoSystRun("Closure", "",
                 electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_central);
       if(!isData){
         DoSystRun("Closure", "SystUpPU",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_PUup);
         DoSystRun("Closure", "SystUpUncl",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi),
                   systweight_central);
         DoSystRun("Closure", "SystUpJES",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetJESUpColl, bjetJESUpColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi),
                   systweight_JESup);
         DoSystRun("Closure", "SystUpJER",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetJERUpColl, bjetJERUpColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi),
                   systweight_JERup);
         DoSystRun("Closure", "SystUpBTag_L",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_LTagup);
         DoSystRun("Closure", "SystUpBTag_BC",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_BCTagup);
         DoSystRun("Closure", "SystUpElEn", EleTElEnUpColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi),
                   systweight_ElEnup);
         DoSystRun("Closure", "SystUpMuEn", electronColl, electronFakeLColl, MuTMuEnUpColl, muonLooseColl, jetColl, bjetColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi),
                   systweight_MuEnup);


         DoSystRun("Closure", "SystDownPU",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_PUdown);
         DoSystRun("Closure", "SystDownUncl",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi),
                   systweight_central);
         DoSystRun("Closure", "SystDownJES",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetJESDownColl, bjetJESDownColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi),
                   systweight_JESdown);
         DoSystRun("Closure", "SystDownJER",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetJERDownColl, bjetJERDownColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi),
                   systweight_JERdown);
         DoSystRun("Closure", "SystDownBTag_L",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_LTagdown);
         DoSystRun("Closure", "SystDownBTag_BC",
                   electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_BCTagdown);
         DoSystRun("Closure", "SystDownElEn", EleTElEnDownColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi),
                   systweight_ElEndown);
         DoSystRun("Closure", "SystDownMuEn", electronColl, electronFakeLColl, MuTMuEnDownColl, muonLooseColl, jetColl, bjetColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi),
                   systweight_MuEndown);

       }
       else if( isData && k_running_nonprompt ){

         DoSystRun("Closure", "SystUpFR", electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_FRup);
         DoSystRun("Closure", "SystDownFR", electronColl, electronFakeLColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_FRdown);
       }
     }//End of SystRun
   }//End of Closure


/////////////////////////////////////////////////////////////////////////////////// 

return;

}// End of execute event loop
  


void Oct2017_EleFRSyst::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Oct2017_EleFRSyst::BeginCycle() throw( LQError ){
  
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

Oct2017_EleFRSyst::~Oct2017_EleFRSyst() {
  
  Message("In Oct2017_EleFRSyst Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}


/*
void Oct2017_EleFRSyst::CheckNvtxProfile(std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  std::vector<snu::KJet> JetVetoColl = SkimJetColl(JetColl, EleLColl, MuLColl, "EleMuVeto");
  int Nvtx=eventbase->GetEvent().nVertices();

  if(EleTColl.size()==1 && EleLColl.size()==1 && MuLColl.size()==0 && JetVetoColl.size()>0){
    if( !(EleTColl.at(0).Pt()>30 && JetVetoColl.at(0).Pt()>40) ) return;

    float MTW = sqrt(2)*sqrt(MET*EleTColl.at(0).Pt()-METx*EleTColl.at(0).Px()-METy*EleTColl.at(0).Py());
    FillHist("Nvtx_gt1j40"+Label, Nvtx, weight, 0., 50., 50);
  }

}
*/
void Oct2017_EleFRSyst::CheckNormCR(std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool UnPreTrig=false, PreTrig=false;
  if     (Option.Contains("UnPreTrig")) UnPreTrig=true;
  else if(Option.Contains("PreTrig"))   PreTrig  =true;

  std::vector<snu::KJet> JetVetoColl = SkimJetColl(JetColl, EleLColl, MuLColl, "EleMuVeto");
  int Nvtx=eventbase->GetEvent().nVertices();

  if(EleTColl.size()==1 && EleLColl.size()==1 && MuLColl.size()==0){
    if( UnPreTrig && EleTColl.at(0).Pt()<30.) return;
    if( PreTrig   && EleTColl.at(0).Pt()<25.) return;

    float MTW = sqrt(2)*sqrt(MET*EleTColl.at(0).Pt()-METx*EleTColl.at(0).Px()-METy*EleTColl.at(0).Py());
   
    if( UnPreTrig ){
      FillHist("MET"+Label, MET, weight, 0., 200., 40);
      FillHist("MTW"+Label, MTW, weight, 0., 200., 40);
      if(MET>50.){
        FillHist("MTW_met50"+Label, MTW, weight, 0., 200., 40);

        //Inclusive Selection CR plot - Incl W+jet test
        if(MTW>50.){
          FillHist("Count_UnPreNormCR"+Label, 0., weight, 0., 4., 4);
          FillHist("PTe_met50mtw50"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
          FillHist("Etae_met50mtw50"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
          FillHist("Nj_met50mtw50"+Label, JetVetoColl.size(), weight, 0., 10., 10);
          FillHist("MET_met50mtw50"+Label, MET, weight, 0., 200., 40);
        }
      }
    }
    if( UnPreTrig || PreTrig ){
      if(JetVetoColl.size()>0 && JetVetoColl.at(0).Pt()>40.){
        FillHist("MET_gt1j40"+Label, MET, weight, 0., 200., 40);
        FillHist("MTW_gt1j40"+Label, MTW, weight, 0., 200., 40);
        FillHist("Nvtx_gt1j40"+Label, Nvtx, weight, 0., 50., 50);
        FillHist("PTe_gt1j40"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);

        if(MET>50.){
          FillHist("MTW_gt1j40met50"+Label, MTW, weight, 0., 200., 40);
          FillHist("Nvtx_gt1j40met50"+Label, Nvtx, weight, 0., 50., 50);

          //e+geq1j Selection CR plot - W+geq1j test
          if(MTW>50.){
            FillHist("Count_UnPreNormCR"+Label, 1., weight, 0., 4., 4); 
            FillHist("Count_PreNormCR"+Label,   0., weight, 0., 2., 2); 

            FillHist("PTe_gt1j40met50mtw50"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
            FillHist("Etae_gt1j40met50mtw50"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
            FillHist("PTj1_gt1j40met50mtw50"+Label, JetVetoColl.at(0).Pt(), weight, 0., 200., 40);
            FillHist("Etaj1_gt1j40met50mtw50"+Label, JetVetoColl.at(0).Eta(), weight, -5., 5., 20);
            FillHist("dRej1_gt1j40met50mtw50"+Label, EleTColl.at(0).DeltaR(JetVetoColl.at(0)), weight, 0., 5., 50);
            FillHist("Nj_gt1j40met50mtw50"+Label, JetVetoColl.size(), weight, 0., 10., 10);
            FillHist("MET_gt1j40met50mtw50"+Label, MET, weight, 0., 200., 40);

            FillHist("Nvtx_gt1j40met50mtw50"+Label, Nvtx, weight, 0., 50., 50);
          }
        }
      }
    }//End of UnPreTrig & PreTrig
  }//End of 1l

  if(EleTColl.size()==2 && EleLColl.size()==2 && MuLColl.size()==0){
    if( UnPreTrig && !(EleTColl.at(0).Pt()>30. && EleTColl.at(1).Pt()>10.) ) return;
    if( PreTrig   && !(EleTColl.at(0).Pt()>25. && EleTColl.at(1).Pt()>10.) ) return;
    if( EleTColl.at(0).Charge()==EleTColl.at(1).Charge()  ) return;
    if( fabs((EleTColl.at(0)+EleTColl.at(1)).M()-91.2)>15 ) return;

    if(UnPreTrig){
      FillHist("Mee_e25e10"+Label, (EleTColl.at(0)+EleTColl.at(1)).M(), weight, 60., 120., 60);
      FillHist("Count_UnPreNormCR"+Label, 2., weight, 0., 4., 4);
    }
    if(UnPreTrig || PreTrig){
      if(JetVetoColl.size()>0 && JetVetoColl.at(0).Pt()>40.){
        FillHist("Mee_e25e10gt1j"+Label, (EleTColl.at(0)+EleTColl.at(1)).M(), weight, 60., 120., 60);
        FillHist("Count_UnPreNormCR"+Label, 3., weight, 0., 4., 4);
        FillHist("Count_PreNormCR"+Label, 1., weight, 0., 2., 2);

        FillHist("Nvtx_e25e10gt1j"+Label, Nvtx, weight, 0., 50., 50);
        FillHist("MET_e25e10gt1j"+Label, MET, weight, 0., 200., 40);
      }
    }
  }//End of 2l

}


void Oct2017_EleFRSyst::OptimiseMETMTWCuts(std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, float MET, float METx, float METy, float weight, TString Option){

  bool SiglPreTrig=false, MultPreTrig=false;
  float MinElePt=0.;
  if     (Option.Contains("SiglPreTrig")) {SiglPreTrig = true; MinElePt=25.;}
  else if(Option.Contains("MultPreTrig")) {MultPreTrig = true; MinElePt=15.;}

  std::vector<snu::KJet> JetVetoColl = SkimJetColl(JetColl, EleLColl, MuLColl, "EleMuVeto");


  if(EleLColl.size()==1 && MuLColl.size()==0 && JetVetoColl.size()>0){
    bool PassJetReq=false;
    float PTCorr=EleLColl.at(0).Pt()*(1.+EleLColl.at(0).PFRelIso(0.3)), fEta=fabs(EleLColl.at(0).Eta());

    //Selection Requirement--------------------------------------------------------------//
    if( EleLColl.at(0).Pt()<MinElePt ) return;
    if( PTCorr<25. ) return;
    for(int j=0; j<JetVetoColl.size(); j++){
      if(JetVetoColl.at(j).Pt()<40) continue;
      if(JetVetoColl.at(j).DeltaR(EleLColl.at(0))<1.0) continue;
      PassJetReq=true; break;
    }
    if(!PassJetReq) return;
    //-----------------------------------------------------------------------------------//

    float MTW = sqrt(2)*sqrt(MET*EleLColl.at(0).Pt()-METx*EleLColl.at(0).Px()-METy*EleLColl.at(0).Py());
    if( EleLColl.size()==1 ){
      FillHist("NEvt_MET_1D", MET, weight, 0., 200., 40);
      FillHist("NEvt_MTW_1D", MTW, weight, 0., 200., 40);
      FillHist("NEvt_METMTW_2D", MET, MTW, weight, 0., 200., 40, 0., 200., 40); 
      if(MET<25) FillHist("NEvt_MTW_met25_1D", MTW, weight, 0., 200., 40);
      if(MTW<35) FillHist("NEvt_MET_mtw35_1D", MET, weight, 0., 200., 40);

      for(int i=1; i<=20; i++){
        if(MET<i*5.){
          for(int j=1; j<=20; j++){
            if(MTW<j*5.){
              FillHist("NEvt_ltMETltMTW_2D", (i-1)*5., (j-1)*5., weight, 0., 100., 20, 0., 100., 20); 
            }
          }
        }
      }
    }
    if( EleTColl.size()==1 ){
      FillHist("NIDEvt_MET_1D", MET, weight, 0., 200., 40);
      FillHist("NIDEvt_MTW_1D", MTW, weight, 0., 200., 40);
      FillHist("NIDEvt_METMTW_2D", MET, MTW, weight, 0., 200., 40, 0., 200., 40); 
      if(MET<25) FillHist("NIDEvt_MTW_met25_1D", MTW, weight, 0., 200., 40);
      if(MTW<35) FillHist("NIDEvt_MET_mtw35_1D", MET, weight, 0., 200., 40);

      for(int i=1; i<=20; i++){
        if(MET<i*5.){
          for(int j=1; j<=20; j++){
            if(MTW<j*5.){
              FillHist("NIDEvt_ltMETltMTW_2D", (i-1)*5., (j-1)*5., weight, 0., 100., 20, 0., 100., 20); 
            }
          }
        }
      }
    }
  }//End of 1e+0mu+geq1j
 
}



void Oct2017_EleFRSyst::OptimiseIsoIPWP(std::vector<snu::KElectron> EleTColl, float weight, TString Option){
  //Purpose : Using DY Data evts, check ROC curve of Iso IP on top of specific MVA WP(without iso req.)
  //EleTColl should be without Iso,IP cuts

  if(                   EleTColl.size()!=2                  ) return;
  if( !(EleTColl.at(0).Charge()!=EleTColl.at(1).Charge())   ) return;
  if( !(fabs((EleTColl.at(0)+EleTColl.at(1)).M()-91.2)<15)  ) return;
  if( !(EleTColl.at(0).Pt()>25. && EleTColl.at(1).Pt()>15.) ) return;

  float IsoWP1=0.06, IsoWP2=0.08, IsoWP3=0.1;
  std::ostringstream s1,s2,s3; s1<<IsoWP1; s2<<IsoWP2; s3<<IsoWP3;
  TString Str_IsoWP1=s1.str(), Str_IsoWP2=s2.str(), Str_IsoWP3=s3.str();
  Str_IsoWP1.ReplaceAll(".","p"); Str_IsoWP2.ReplaceAll(".","p"); Str_IsoWP3.ReplaceAll(".","p");


  for(int i=0; i<EleTColl.size(); i++){
    float Iso=EleTColl.at(i).PFRelIso(0.3);
    float d0=fabs(EleTColl.at(i).dxy()), dz=fabs(EleTColl.at(i).dz());

    FillHist("Data_Iso", EleTColl.at(i).PFRelIso(0.3), weight, 0., 0.3, 300);
    FillHist("Data_d0", d0, weight, 0., 0.1, 20);
    FillHist("Data_dz", dz, weight, 0., 0.2, 40);

    //IsoWP Study
    int NIsoWP=40;
    for(int it_iso=1; it_iso<=NIsoWP; it_iso++){
      FillHist("NEleSumW_Iso", (it_iso-1.)*0.00501, weight, 0., 0.2, 40);
      if(EleTColl.at(i).PFRelIso(0.3)<it_iso*0.005){
        FillHist("NEleIDSumW_Iso", (it_iso-1.)*0.00501, weight, 0., 0.2, 40);
      }
    }


    //IP ROC Plane
    int Nd0WP=10, NdzWP=20;
    for(int it_d0=1; it_d0<=Nd0WP; it_d0++){
      for(int it_dz=1; it_dz<=NdzWP; it_dz++){
        if(EleTColl.at(i).PFRelIso(0.3)<IsoWP1){
          FillHist("NEleSumW_d0dz_iso"+Str_IsoWP1, (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
          if(d0<it_d0*0.005 && dz< it_dz*0.005){
            FillHist("NEleIDSumW_d0dz_iso"+Str_IsoWP1, (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
          }
        }
        if(EleTColl.at(i).PFRelIso(0.3)<IsoWP2){
          FillHist("NEleSumW_d0dz_iso"+Str_IsoWP2, (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
          if(d0<it_d0*0.005 && dz< it_dz*0.005){
            FillHist("NEleIDSumW_d0dz_iso"+Str_IsoWP2, (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
          }
        }
        if(EleTColl.at(i).PFRelIso(0.3)<IsoWP3){
          FillHist("NEleSumW_d0dz_iso"+Str_IsoWP3, (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
          if(d0<it_d0*0.005 && dz< it_dz*0.005){
            FillHist("NEleIDSumW_d0dz_iso"+Str_IsoWP3, (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
          }
        }
      }
    }
  }//End of EleLoop


}

void Oct2017_EleFRSyst::CheckPromptVariation(std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){
  //BJet : NoVetoColl,  Jet: No matter whether lepton vetoed or not, it will veto in code any way.

  bool SiglPreTrig=false, MultPreTrig=false;
  float MinElePt=0.;
  if     (Option.Contains("SiglPreTrig")) {SiglPreTrig = true; MinElePt=25.;}
  else if(Option.Contains("MultPreTrig")) {MultPreTrig = true; MinElePt=15.;}

  const int NPtEdges=6;
  float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100.};

  std::vector<snu::KJet> JetVetoColl = SkimJetColl(JetColl, EleLColl, MuLColl, "EleMuVeto");
  std::vector<snu::KJet> BJetVetoColl = SelBJets(JetVetoColl, "Medium");
  int NBJets = BJetVetoColl.size();


  if(EleTColl.size()==1 && EleLColl.size()==1 && MuLColl.size()==0 && JetVetoColl.size()>0){
    bool PassJetReq=false;
    float PTCorr=ConeCorrectedPT(EleLColl.at(0), 0.06), fEta=fabs(EleLColl.at(0).Eta());
  
    //Selection Requirement--------------------------------------------------------------//
    if( EleLColl.at(0).Pt()<MinElePt ) return;
    if( PTCorr<25. ) return;
    for(int j=0; j<JetVetoColl.size(); j++){
      if(JetVetoColl.at(j).Pt()<40) continue;
      if(JetVetoColl.at(j).DeltaR(EleLColl.at(0))<1.0) continue;
      PassJetReq=true; break;
    }
    if(!PassJetReq) return;
    float MTW = sqrt(2)*sqrt(MET*EleLColl.at(0).Pt()-METx*EleLColl.at(0).Px()-METy*EleLColl.at(0).Py());

    if( !(MET<25 && MTW<25) ) return;
    //-----------------------------------------------------------------------------------//
    //bool IsNearB = IsNearBJet(EleLColl.at(0), BJetColl);
    if     (fEta<0.8)   FillHist("EleB1IDSumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    else if(fEta<1.479) FillHist("EleB2IDSumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    else                FillHist("EleEIDSumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  }//End of 1e+0mu+geq1j
 
}


void Oct2017_EleFRSyst::MeasureFakeRate(std::vector<snu::KElectron> ElePreColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, int NPtEdges, float PtEdges[], TString LooseID, TString TightID, TString Label, TString Option){
  //BJet : NoVetoColl,  Jet: No matter whether lepton vetoed or not, it will veto in code any way.

  std::vector<snu::KElectron> EleLColl;
    for(int i=0; i<ElePreColl.size(); i++){ if(PassIDCriteria(ElePreColl.at(i),LooseID)) EleLColl.push_back(ElePreColl.at(i)); }

  bool SiglPreTrig=false, MultPreTrig=false;
  float MinElePt=0.;
  if     (Option.Contains("SiglPreTrig")) {SiglPreTrig = true; MinElePt=25.;}
  else if(Option.Contains("MultPreTrig")) {MultPreTrig = true; MinElePt=15.;}

  std::vector<snu::KJet> JetVetoColl = SkimJetColl(JetColl, EleLColl, MuLColl, "EleMuVeto");
  std::vector<snu::KJet> BJetVetoColl = SelBJets(JetVetoColl, "Medium");
  int NBJets = BJetVetoColl.size();


  if(EleLColl.size()==1 && MuLColl.size()==0 && JetVetoColl.size()>0){
    bool PassJetReq=false;
    float PTCorr=ConeCorrectedPT(EleLColl.at(0), 0.06), fEta=fabs(EleLColl.at(0).Eta());

    //Selection Requirement--------------------------------------------------------------//
    if( EleLColl.at(0).Pt()<MinElePt ) return;
    if( PTCorr<25. ) return;
    for(int j=0; j<JetVetoColl.size(); j++){
      if(JetVetoColl.at(j).Pt()<40) continue;
      if(JetVetoColl.at(j).DeltaR(EleLColl.at(0))<1.0) continue;
      PassJetReq=true; break;
    }
    if(!PassJetReq) return;
    float MTW = sqrt(2)*sqrt(MET*EleLColl.at(0).Pt()-METx*EleLColl.at(0).Px()-METy*EleLColl.at(0).Py());

    if( !(MET<25 && MTW<25) ) return;
    //if( MET<30 && MTW<25 ) return;//Final Optimised
    //-----------------------------------------------------------------------------------//

    //bool IsNearB = IsNearBJet(EleLColl.at(0), BJetColl);
    if(fEta<0.8){
      FillHist("EleB1SumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
   
      if(NBJets>0){ FillHist("EleB1SumW_BjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
                   if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_BjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1); } 
    }
    else if(fEta<1.479){
      FillHist("EleB2SumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
   
      if(NBJets>0){ FillHist("EleB2SumW_BjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
                   if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_BjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1); } 
    }
    else{
      FillHist("EleESumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
   
      if(NBJets>0){ FillHist("EleESumW_BjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
                   if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_BjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1); } 
    }
  }//End of 1e+0mu+geq1j
 
}


void Oct2017_EleFRSyst::ScanFakeRate(std::vector<snu::KElectron> EleColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, int NMVACuts, float MVACuts[], int NIsoCuts, float IsoCuts[], int NPtEdges, float PtEdges[], TString Option){
  //EleColl should be without iso, mva requirement, since those are objectives of the optimisation
  //We see near bjet, so bjet should be noveto, jet should be noveto.(optimising proc.->Veto cond is ambiguous, and 1.0 cut is already there)

  float MinElePt=0.; TString TightID="POGMVAMIP";
  bool Scan1D=false, Scan2D=false;
  if(Option.Contains("SiglPreTrig")) MinElePt=25.;
  if(Option.Contains("MultPreTrig")) MinElePt=15.;
  if(Option.Contains("2D"))          Scan2D  =true;
  if(Option.Contains("1D"))          Scan1D  =true;

      
  for(int it_iso=0; it_iso<NIsoCuts; it_iso++){
    for(int it_mva=0; it_mva<NMVACuts; it_mva++){

      //Selection Requirement--------------------------------------------------------------//
      int NLooseEle=0;
      if( IsoCuts[it_iso]==0. || MVACuts[it_mva]==1. ) continue;
      for(int i_ele=0; i_ele<EleColl.size(); i_ele++){
        if(EleColl.at(i_ele).PFRelIso(0.3)<IsoCuts[it_iso] && EleColl.at(i_ele).MVA()>MVACuts[it_mva]) NLooseEle++;
      }
      if(NLooseEle!=1) continue;
      if(MuLColl.size()!=0) continue;
      if(EleColl.at(0).Pt()<MinElePt) continue;
      bool PassJetReq=false; 
      for(int j=0; j<JetColl.size(); j++){
        if(JetColl.at(j).Pt()<40) continue;
        if(JetColl.at(j).DeltaR(EleColl.at(0))<1.0) continue;
        PassJetReq=true;
      }
      if(!PassJetReq) continue;
      float MTW = sqrt(2)*sqrt(MET*EleColl.at(0).Pt()-METx*EleColl.at(0).Px()-METy*EleColl.at(0).Py());
      if( !(MET<30 && MTW<25) ) continue;
      //-----------------------------------------------------------------------------------//

      float PTCorr       = ConeCorrectedPT(EleColl.at(0), 0.1);
      bool  IsNearB      = IsNearBJet(EleColl.at(0), BJetColl);
      bool  TrigUnbiased = PTCorr>MinElePt*(1.+IsoCuts[it_iso]);
      //TrigUnbiased Cut is used only for 2D scan, since 2D scan is avg FR of each cut value, but we are interested plateau value of 1D dist.

      std::ostringstream s2; s2<<IsoCuts[it_iso];     std::ostringstream s1; s1<<MVACuts[it_mva];
      TString Str_IsoCut=s2.str();                    TString Str_MVACut=s1.str();
      Str_IsoCut.ReplaceAll(".","p");                 Str_MVACut.ReplaceAll(".","p");  Str_MVACut.ReplaceAll("-","m");



      if(EleColl.at(0).PFRelIso(0.3)<IsoCuts[it_iso] && EleColl.at(0).MVA()>MVACuts[it_mva]){

        if(fabs(EleColl.at(0).Eta())<0.8){//Barrel1
          if(Scan2D && TrigUnbiased){
            FillHist("EleB1SumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
            if(PassIDCriteria(EleColl.at(0), TightID)) FillHist("EleB1IDSumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
  
            if(IsNearB){
              FillHist("EleB1SumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
              if(PassIDCriteria(EleColl.at(0), TightID)) FillHist("EleB1IDSumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
            }
          }
          if(Scan1D){
            FillHist("EleB1SumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
            if(PassIDCriteria(EleColl.at(0), TightID)) FillHist("EleB1IDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
  
            if(IsNearB){
              FillHist("EleB1SumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
              if(PassIDCriteria(EleColl.at(0), TightID)) FillHist("EleB1IDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
            }
          }
        }
        else if(fabs(EleColl.at(0).Eta())<1.479){//Barrel2
          if(Scan2D && TrigUnbiased){
            FillHist("EleB2SumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
            if(PassIDCriteria(EleColl.at(0), TightID)) FillHist("EleB2IDSumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
  
            if(IsNearB){
              FillHist("EleB2SumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
              if(PassIDCriteria(EleColl.at(0), TightID)) FillHist("EleB2IDSumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
            }
          }
          if(Scan1D){
            FillHist("EleB2SumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
            if(PassIDCriteria(EleColl.at(0), TightID)) FillHist("EleB2IDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
  
            if(IsNearB){
              FillHist("EleB2SumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
              if(PassIDCriteria(EleColl.at(0), TightID)) FillHist("EleB2IDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
            }
          }
        }
        else if(fabs(EleColl.at(0).Eta())<2.5){//EndCap
          if(Scan2D && TrigUnbiased){
            FillHist("EleESumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
            if(PassIDCriteria(EleColl.at(0), TightID)) FillHist("EleEIDSumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
  
            if(IsNearB){
              FillHist("EleESumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
              if(PassIDCriteria(EleColl.at(0), TightID)) FillHist("EleEIDSumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
            }
          }
          if(Scan1D){
            FillHist("EleESumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
            if(PassIDCriteria(EleColl.at(0), TightID)) FillHist("EleEIDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
  
            if(IsNearB){
              FillHist("EleESumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
              if(PassIDCriteria(EleColl.at(0), TightID)) FillHist("EleEIDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
            }
          }
        }
      }//End of if Iso<Cut && MVA<Cut
    }//End of MVACut loop
  }//End of IsoCut Loop


}


void Oct2017_EleFRSyst::MeasureHighdXYFakeRate(std::vector<snu::KElectron> EleLowd0TColl, std::vector<snu::KElectron> EleLowd0LColl, std::vector<snu::KElectron> EleHighd0TColl, std::vector<snu::KElectron> EleHighd0LColl, std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, std::vector<snu::KTruth> TruthColl, float MET, float METx, float METy, float weight, TString Label, TString Option){
     //jet:LeptonNoVetoedColl both for ambuguity and near jet checking

  //Series of Cuts for FR measurement--------------------------------------------------------------------------------//
  bool PassSelHighd0=true, PassSelLowd0=true; float MTW=0.;
  if( EleHighd0LColl.size()!=1 ) PassSelHighd0=false;
  if( PassSelHighd0 && !(EleHighd0LColl.at(0).Pt()>25. && fabs(EleHighd0LColl.at(0).dxySig())>5.) ) PassSelHighd0=false;
  if( PassSelHighd0 && !(JetColl.size()>0 && JetColl.at(0).Pt()>40) ) PassSelHighd0=false;
  if( PassSelHighd0 ){
    bool HasAway40Jet=false;
    for( int i=0; i<JetColl.size(); i++){
      if(JetColl.at(i).Pt()<40) continue;
      if(JetColl.at(i).DeltaR(EleHighd0LColl.at(0))>1.0) {HasAway40Jet=true; break;}
    }
    if(!HasAway40Jet) PassSelHighd0=false;
  }
  if( !isData && PassSelHighd0 && !(k_sample_name.Contains("qcd") || k_sample_name.Contains("QCD")) ){
    int EleType=isData? 0 : GetLeptonType(EleHighd0LColl.at(0),TruthColl);
    if( EleType<0 && EleType>-5 ) PassSelHighd0=false;// for MC, we only subtract that originate from prompts. Fakes are estimated in data.
  }


  if( EleLowd0LColl.size()!=1 ) PassSelLowd0=false;
  if( PassSelLowd0 && !(EleLowd0LColl.at(0).Pt()>25.) ) PassSelLowd0=false;
  if( PassSelLowd0 && !(JetColl.size()>0 && JetColl.at(0).Pt()>40) ) PassSelLowd0=false;
  if( PassSelLowd0 ){
    bool HasAway40Jet=false;
    for( int i=0; i<JetColl.size(); i++){
      if(JetColl.at(i).Pt()<40) continue;
      if(JetColl.at(i).DeltaR(EleLowd0LColl.at(0))>1.0) {HasAway40Jet=true; break;}
    }
    if(!HasAway40Jet) PassSelLowd0=false;
  }
  if( !isData && PassSelLowd0 && !(k_sample_name.Contains("qcd") || k_sample_name.Contains("QCD")) ){
    int EleType=isData? 0 : GetLeptonType(EleLowd0LColl.at(0),TruthColl);
    if( EleType<0 && EleType>-5 ) PassSelLowd0=false;// for MC, we only subtract that originate from prompts. Fakes are estimated in data.
  }


  //-----------------------------------------------------------------------------------------------------------------//
  
  if( PassSelHighd0 ){
    
    float PTCorr=EleHighd0LColl.at(0).Pt()*(1.+EleHighd0LColl.at(0).PFRelIso(0.3));
    float fEta=fabs(EleHighd0LColl.at(0).Eta());
    const int NPtEdges=8; float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100., 150., 200.};
    bool  IsNearB = IsNearBJet(EleHighd0LColl.at(0), BJetColl);
    float MTW = sqrt(2)*sqrt(MET*EleHighd0LColl.at(0).Pt()-METx*EleHighd0LColl.at(0).Px()-METy*EleHighd0LColl.at(0).Py());
    bool PassMETMTW=(MET<25 && MTW<35);


    //Pass Loose
    if(fEta<0.8){
      FillHist("EleHighd0B1SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      if(IsNearB) FillHist("EleHighd0B1SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
    }
    else if(fEta<1.479){
      FillHist("EleHighd0B2SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      if(IsNearB) FillHist("EleHighd0B2SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
    }
    else{
      FillHist("EleHighd0ESumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      if(IsNearB) FillHist("EleHighd0ESumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
    }

    if(PassMETMTW){
      if(fEta<0.8){
        FillHist("EleHighd0METMTWB1SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        if(IsNearB) FillHist("EleHighd0METMTWB1SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else if(fEta<1.479){
        FillHist("EleHighd0METMTWB2SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        if(IsNearB) FillHist("EleHighd0METMTWB2SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else{
        FillHist("EleHighd0METMTWESumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        if(IsNearB) FillHist("EleHighd0METMTWESumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }

    //Pass IDCriteria
    if(EleHighd0TColl.size()==1){

      if(fEta<0.8){
        FillHist("EleHighd0B1IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        if(IsNearB) FillHist("EleHighd0B1IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else if(fEta<1.479){
        FillHist("EleHighd0B2IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        if(IsNearB) FillHist("EleHighd0B2IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else{
        FillHist("EleHighd0EIDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        if(IsNearB) FillHist("EleHighd0EIDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      }
  
      if(PassMETMTW){
        if(fEta<0.8){
          FillHist("EleHighd0METMTWB1IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
          if(IsNearB) FillHist("EleHighd0METMTWB1IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        }
        else if(fEta<1.479){
          FillHist("EleHighd0METMTWB2IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
          if(IsNearB) FillHist("EleHighd0METMTWB2IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        }
        else{
          FillHist("EleHighd0METMTWEIDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
          if(IsNearB) FillHist("EleHighd0METMTWEIDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
    }//End of highd0 Tight
  }//End of High d0 MeasReg

  //QCDMC Highd0 to Lowd0 SF Measurement-----------------------------------------------------------------------------//
  if( !isData && (k_sample_name.Contains("qcd") || k_sample_name.Contains("QCD") || k_sample_name.Contains("TT_powheg")) && PassSelLowd0 ){

    float PTCorr=EleLowd0LColl.at(0).Pt()*(1.+EleLowd0LColl.at(0).PFRelIso(0.3));
    float fEta=fabs(EleLowd0LColl.at(0).Eta());
    const int NPtEdges=8; float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100., 150., 200.};
    bool  IsNearB = IsNearBJet(EleLowd0LColl.at(0), BJetColl);
    float MTW = sqrt(2)*sqrt(MET*EleLowd0LColl.at(0).Pt()-METx*EleLowd0LColl.at(0).Px()-METy*EleLowd0LColl.at(0).Py());
    bool PassMETMTW=(MET<25 && MTW<35);


    //Pass Loose
    if(fEta<0.8){
      FillHist("EleLowd0B1SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      if(IsNearB) FillHist("EleLowd0B1SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
    }
    else if(fEta<1.479){
      FillHist("EleLowd0B2SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      if(IsNearB) FillHist("EleLowd0B2SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
    }
    else{
      FillHist("EleLowd0ESumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      if(IsNearB) FillHist("EleLowd0ESumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
    }

    if(PassMETMTW){
      if(fEta<0.8){
        FillHist("EleLowd0METMTWB1SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        if(IsNearB) FillHist("EleLowd0METMTWB1SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else if(fEta<1.479){
        FillHist("EleLowd0METMTWB2SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        if(IsNearB) FillHist("EleLowd0METMTWB2SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else{
        FillHist("EleLowd0METMTWESumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        if(IsNearB) FillHist("EleLowd0METMTWESumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }

    //Pass IDCriteria
    if(EleLowd0TColl.size()==1){

      if(fEta<0.8){
        FillHist("EleLowd0B1IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        if(IsNearB) FillHist("EleLowd0B1IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else if(fEta<1.479){
        FillHist("EleLowd0B2IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        if(IsNearB) FillHist("EleLowd0B2IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else{
        FillHist("EleLowd0EIDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        if(IsNearB) FillHist("EleLowd0EIDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
      }
  
      if(PassMETMTW){
        if(fEta<0.8){
          FillHist("EleLowd0METMTWB1IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
          if(IsNearB) FillHist("EleLowd0METMTWB1IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        }
        else if(fEta<1.479){
          FillHist("EleLowd0METMTWB2IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
          if(IsNearB) FillHist("EleLowd0METMTWB2IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        }
        else{
          FillHist("EleLowd0METMTWEIDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
          if(IsNearB) FillHist("EleLowd0METMTWEIDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
    }//End of highd0 Tight
  }//End of QCDMC Lowd0 MeasReg



}



void Oct2017_EleFRSyst::CheckEleFRSyst(std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> ElePreColl, std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetNoVetoColl, std::vector<snu::KJet> BJetNoVetoColl, float MET, float METx, float METy, float weight, int NPtEdges, float PtEdges[], TString LooseID, TString TightID, TString Label, TString Option){

  // FR Syst.
  //1. Prompt Uncertainty
  //   - Norm Unc. ~15% conservatively;
  //    In syst variation for tight sel, variations are included in 15%(Nvtx, Xsec, JES, JER, ElEn, Uncl)
  //2. Impact from cut value : MET, MTW cut variation
  //3. Balance Jet Systematic : dR, jet Pt. Njet
  //4. Variance when B-jet present
  //5. Prompt Rate: What happens if prompt is included 

  //How
  //1.   By weight variation
  //2-4. Better with in code.
  //5.   Measure Prompt Rate in DY events, and apply it.

  std::vector<snu::KElectron> EleLColl;
    for(int i=0; i<ElePreColl.size(); i++){ if(PassIDCriteria(ElePreColl.at(i),LooseID)) EleLColl.push_back(ElePreColl.at(i)); }

  float MinElePt=25.;
  std::vector<snu::KJet> JetColl  = SkimJetColl(JetNoVetoColl, EleLColl, MuLColl, "EleMuVeto");
  std::vector<snu::KJet> BJetColl = SelBJets(JetColl, "Medium");
  int NBJets = BJetNoVetoColl.size();


  if(EleLColl.size()==1 && MuLColl.size()==0 && JetColl.size()>0){
    float PTCorr=ConeCorrectedPT(EleLColl.at(0), 0.06), fEta=fabs(EleLColl.at(0).Eta());

    //Selection Requirement--------------------------------------------------------------//
    bool Selection_cent=true;
    bool PassJetReq=false, PassJetReq_30=false, PassJetReq_40=false, PassJetReq_60=false;
    int Naway_30=0, Naway_40=0, Naway_60=0;
    

    if( !(EleLColl.at(0).Pt()>MinElePt && PTCorr>25) ) return;//Applies to All
    for(int j=0; j<JetColl.size(); j++){
      if(JetColl.at(j).Pt()>30 && JetColl.at(j).DeltaR(EleLColl.at(0))>1.0) {PassJetReq_30 = true; Naway_30++;}
      if(JetColl.at(j).Pt()>40 && JetColl.at(j).DeltaR(EleLColl.at(0))>1.0) {PassJetReq_40 = true; Naway_40++;}
      if(JetColl.at(j).Pt()>60 && JetColl.at(j).DeltaR(EleLColl.at(0))>1.0) {PassJetReq_60 = true; Naway_60++;}
    }   
    FillHist("Nj_1eCut", JetColl.size(), weight, 0., 10., 10);
    FillHist("Nawayj30_1eCut", Naway_30, weight, 0., 10., 10);
    FillHist("Nawayj40_1eCut", Naway_40, weight, 0., 10., 10);
    FillHist("Nawayj60_1eCut", Naway_60, weight, 0., 10., 10);
    if( !PassJetReq_30 ) return;//If !30 then that means !40, !60 as well.

    float MTW = sqrt(2)*sqrt(MET*EleLColl.at(0).Pt()-METx*EleLColl.at(0).Px()-METy*EleLColl.at(0).Py());
    //-----------------------------------------------------------------------------------//

    //bool IsNearB = IsNearBJet(EleLColl.at(0), BJetColl);

    //Central
    if(PassJetReq_40 && MET<25 && MTW<25){
      FillHist("EleAllSumW_PT_FR1D_Cent"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleAllIDSumW_PT_FR1D_Cent"+Label, PTCorr, weight, PtEdges, NPtEdges-1);

      if(fEta<0.8){
        FillHist("EleB1SumW_PT_FR1D_Cent"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_PT_FR1D_Cent"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else if(fEta<1.479){
        FillHist("EleB2SumW_PT_FR1D_Cent"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_PT_FR1D_Cent"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else{
        FillHist("EleESumW_PT_FR1D_Cent"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_PT_FR1D_Cent"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }

    //BJet Presence ; Src Dep. Stability Test
    if(PassJetReq_40 && NBJets>0 && MET<25 && MTW<25){
      FillHist("EleAllSumW_PT_FR1D_HasB"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleAllIDSumW_PT_FR1D_HasB"+Label, PTCorr, weight, PtEdges, NPtEdges-1);

      if(fEta<0.8){
        FillHist("EleB1SumW_PT_FR1D_HasB"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_PT_FR1D_HasB"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else if(fEta<1.479){
        FillHist("EleB2SumW_PT_FR1D_HasB"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_PT_FR1D_HasB"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else{
        FillHist("EleESumW_PT_FR1D_HasB"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_PT_FR1D_HasB"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }

    //Away Jet PT Variation; Mother Energy Scale Stability Test
    if(PassJetReq_30 && MET<25 && MTW<25){
      FillHist("EleAllSumW_PT_FR1D_JetPt30"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleAllIDSumW_PT_FR1D_JetPt30"+Label, PTCorr, weight, PtEdges, NPtEdges-1);

      if(fEta<0.8){
        FillHist("EleB1SumW_PT_FR1D_JetPt30"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_PT_FR1D_JetPt30"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else if(fEta<1.479){
        FillHist("EleB2SumW_PT_FR1D_JetPt30"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_PT_FR1D_JetPt30"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else{
        FillHist("EleESumW_PT_FR1D_JetPt30"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_PT_FR1D_JetPt30"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }
    if(PassJetReq_60 && MET<25 && MTW<25){
      FillHist("EleAllSumW_PT_FR1D_JetPt60"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleAllIDSumW_PT_FR1D_JetPt60"+Label, PTCorr, weight, PtEdges, NPtEdges-1);

      if(fEta<0.8){
        FillHist("EleB1SumW_PT_FR1D_JetPt60"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_PT_FR1D_JetPt60"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else if(fEta<1.479){
        FillHist("EleB2SumW_PT_FR1D_JetPt60"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_PT_FR1D_JetPt60"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else{
        FillHist("EleESumW_PT_FR1D_JetPt60"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_PT_FR1D_JetPt60"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }

    //MET & MTW Cut Bias Test
    if(PassJetReq_40 && MET<25){
      FillHist("EleAllSumW_PT_FR1D_NoMTW"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleAllIDSumW_PT_FR1D_NoMTW"+Label, PTCorr, weight, PtEdges, NPtEdges-1);

      if(fEta<0.8){
        FillHist("EleB1SumW_PT_FR1D_NoMTW"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_PT_FR1D_NoMTW"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else if(fEta<1.479){
        FillHist("EleB2SumW_PT_FR1D_NoMTW"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_PT_FR1D_NoMTW"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else{
        FillHist("EleESumW_PT_FR1D_NoMTW"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_PT_FR1D_NoMTW"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }
    if(PassJetReq_40 && MTW<25){
      FillHist("EleAllSumW_PT_FR1D_NoMET"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleAllIDSumW_PT_FR1D_NoMET"+Label, PTCorr, weight, PtEdges, NPtEdges-1);

      if(fEta<0.8){
        FillHist("EleB1SumW_PT_FR1D_NoMET"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_PT_FR1D_NoMET"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else if(fEta<1.479){
        FillHist("EleB2SumW_PT_FR1D_NoMET"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_PT_FR1D_NoMET"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else{
        FillHist("EleESumW_PT_FR1D_NoMET"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_PT_FR1D_NoMET"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }

    //Prompt Variation: Prompt Variation ~<15%(Xsec+JES+JER+Uncl+Nvtx+ElEn)
    if(PassJetReq_40 && MET<25 && MTW<25){
      float weight_PrUp   = isData? weight:weight*1.15;
      float weight_PrDown = isData? weight:weight*0.85;

      FillHist("EleAllSumW_PT_FR1D_PrUp"+Label, PTCorr, weight_PrUp, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleAllIDSumW_PT_FR1D_PrUp"+Label, PTCorr, weight_PrUp, PtEdges, NPtEdges-1);
      FillHist("EleAllSumW_PT_FR1D_PrDown"+Label, PTCorr, weight_PrDown, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleAllIDSumW_PT_FR1D_PrDown"+Label, PTCorr, weight_PrDown, PtEdges, NPtEdges-1);

      if(fEta<0.8){
        FillHist("EleB1SumW_PT_FR1D_PrUp"+Label, PTCorr, weight_PrUp, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_PT_FR1D_PrUp"+Label, PTCorr, weight_PrUp, PtEdges, NPtEdges-1);
        FillHist("EleB1SumW_PT_FR1D_PrDown"+Label, PTCorr, weight_PrDown, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_PT_FR1D_PrDown"+Label, PTCorr, weight_PrDown, PtEdges, NPtEdges-1);

      }
      else if(fEta<1.479){
        FillHist("EleB2SumW_PT_FR1D_PrUp"+Label, PTCorr, weight_PrUp, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_PT_FR1D_PrUp"+Label, PTCorr, weight_PrUp, PtEdges, NPtEdges-1);
        FillHist("EleB2SumW_PT_FR1D_PrDown"+Label, PTCorr, weight_PrDown, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_PT_FR1D_PrDown"+Label, PTCorr, weight_PrDown, PtEdges, NPtEdges-1);

      }
      else{
        FillHist("EleESumW_PT_FR1D_PrUp"+Label, PTCorr, weight_PrUp, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_PT_FR1D_PrUp"+Label, PTCorr, weight_PrUp, PtEdges, NPtEdges-1);
        FillHist("EleESumW_PT_FR1D_PrDown"+Label, PTCorr, weight_PrDown, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_PT_FR1D_PrDown"+Label, PTCorr, weight_PrDown, PtEdges, NPtEdges-1);

      }
    }

    //1away jet
    if(PassJetReq_40 && MET<25 && MTW<25 && Naway_40>1){
      FillHist("EleAllSumW_PT_FR1D_away1j"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleAllIDSumW_PT_FR1D_away1j"+Label, PTCorr, weight, PtEdges, NPtEdges-1);

      if(fEta<0.8){
        FillHist("EleB1SumW_PT_FR1D_away1j"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_PT_FR1D_away1j"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else if(fEta<1.479){
        FillHist("EleB2SumW_PT_FR1D_away1j"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_PT_FR1D_away1j"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else{
        FillHist("EleESumW_PT_FR1D_away1j"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_PT_FR1D_away1j"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }

  }//End of 1e+0mu+geq1j
  
}

void Oct2017_EleFRSyst::CheckTrilepCRs(std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  if( !(EleLColl.size()==1 && MuLColl.size()==2) ) return;
  if( !(EleTColl.size()==1 && MuTColl.size()==2) ) return;
  if( !(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) ) return;
  if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && MuTColl.at(1).Pt()>10) ) return;

  float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
  float MTW = sqrt(2)*sqrt(MET*EleTColl.at(0).Pt()-METx*EleTColl.at(0).Px()-METy*EleTColl.at(0).Py());
  float M3l=(EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M();
  if(Mmumu<12) return;
 
  
  bool HasBJet=false, OnZ=false, OnZG=false, WZSel=false, ZGSel=false, ttZSel=false, ZbbSel=false, AN_Sideband=false;
  if( BJetColl.size()!=0 )                      HasBJet = true;
  if( fabs(Mmumu-91.2)<10     )                 OnZ     = true;
  if( fabs(M3l-91.2)<10       )                 OnZG    = true;
  if( OnZ && M3l>101.2 && MET>50 )              WZSel   = true;
  if( OnZG && Mmumu<81.2 && MET<50 )            ZGSel   = true;
  if( OnZ && HasBJet && JetColl.size()>2)       ttZSel  = true;
  if( OnZ && HasBJet && JetColl.size()<3)       ZbbSel  = true;
  if( Mmumu>45. && HasBJet && JetColl.size()>2 ) AN_Sideband=true;
   


  
  //General 3l Selection
  if(!HasBJet){
    FillHist("PTe_3lOS"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("Etae_3lOS"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("PTmu1_3lOS"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("PTmu2_3lOS"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("Etamu1_3lOS"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("Etamu2_3lOS"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
    FillHist("dRemu1_3lOS"+Label, MuTColl.at(0).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
    FillHist("dRemu2_3lOS"+Label, MuTColl.at(1).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
    FillHist("dRmu1mu2_3lOS"+Label, MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
  
    FillHist("Mmumu_3lOS"+Label, Mmumu, weight, 0., 200., 40);
    FillHist("M3l_3lOS"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 100);
    FillHist("Nj_3lOS"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_3lOS"+Label, BJetColl.size(), weight, 0., 10., 10);
    FillHist("MET_3lOS"+Label, MET, weight, 0., 200., 20);
    FillHist("MTW_3lOS"+Label, MTW, weight, 0., 200., 20);
  }
  //3l OnZ ; Fake CR
  if(!HasBJet && OnZ){
    FillHist("PTe_3lOSOnZ"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("Etae_3lOSOnZ"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("PTmu1_3lOSOnZ"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("PTmu2_3lOSOnZ"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("Etamu1_3lOSOnZ"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("Etamu2_3lOSOnZ"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
  
    FillHist("Mmumu_3lOSOnZ"+Label, Mmumu, weight, 60., 120., 30);
    FillHist("M3l_3lOSOnZ"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 100);
    FillHist("Nj_3lOSOnZ"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_3lOSOnZ"+Label, BJetColl.size(), weight, 0., 10., 10);
    FillHist("MET_3lOSOnZ"+Label, MET, weight, 0., 300., 30);
    FillHist("MTW_3lOSOnZ"+Label, MTW, weight, 0., 200., 20);
  }
  //3l OnZ & HasBJet
  if(HasBJet && OnZ){
    FillHist("PTe_3lOSOnZHasB"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("Etae_3lOSOnZHasB"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("PTmu1_3lOSOnZHasB"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
    FillHist("PTmu2_3lOSOnZHasB"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("Etamu1_3lOSOnZHasB"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("Etamu2_3lOSOnZHasB"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
  
    FillHist("Mmumu_3lOSOnZHasB"+Label, Mmumu, weight, 60., 120., 30);
    FillHist("M3l_3lOSOnZHasB"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 100);
    FillHist("Nj_3lOSOnZHasB"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_3lOSOnZHasB"+Label, BJetColl.size(), weight, 0., 10., 10);
    FillHist("MET_3lOSOnZHasB"+Label, MET, weight, 0., 200., 40);
    FillHist("MTW_3lOSOnZHasB"+Label, MTW, weight, 0., 200., 40);
  }
  if(WZSel){
    FillHist("PTe_WZSel"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("Etae_WZSel"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("PTmu1_WZSel"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("PTmu2_WZSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("Etamu1_WZSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("Etamu2_WZSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
  
    FillHist("Mmumu_WZSel"+Label, Mmumu, weight, 60., 120., 30);
    FillHist("M3l_WZSel"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 100);
    FillHist("Nj_WZSel"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_WZSel"+Label, BJetColl.size(), weight, 0., 10., 10);
    FillHist("MET_WZSel"+Label, MET, weight, 0., 200., 20);
    FillHist("MTW_WZSel"+Label, MTW, weight, 0., 200., 20);
  }
  if(ZGSel){
    FillHist("PTe_ZGSel"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("Etae_ZGSel"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("PTmu1_ZGSel"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("PTmu2_ZGSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("Etamu1_ZGSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("Etamu2_ZGSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
  
    FillHist("Mmumu_ZGSel"+Label, Mmumu, weight, 0., 200., 40);
    FillHist("M3l_ZGSel"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 60., 120., 30);
    FillHist("Nj_ZGSel"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_ZGSel"+Label, BJetColl.size(), weight, 0., 10., 10);
    FillHist("MET_ZGSel"+Label, MET, weight, 0., 200., 20);
    FillHist("MTW_ZGSel"+Label, MTW, weight, 0., 200., 20);
  }
  if(ttZSel){
    FillHist("PTe_ttZSel"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("Etae_ttZSel"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("PTmu1_ttZSel"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
    FillHist("PTmu2_ttZSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("Etamu1_ttZSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("Etamu2_ttZSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
    FillHist("dRemu1_ttZSel"+Label, MuTColl.at(0).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
    FillHist("dRemu2_ttZSel"+Label, MuTColl.at(1).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
    FillHist("dRmu1mu2_ttZSel"+Label, MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);

  
    FillHist("Mmumu_ttZSel"+Label, Mmumu, weight, 60., 120., 30);
    FillHist("M3l_ttZSel"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 50);
    FillHist("Nj_ttZSel"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_ttZSel"+Label, BJetColl.size(), weight, 0., 10., 10);
    FillHist("MET_ttZSel"+Label, MET, weight, 0., 300., 30);
    FillHist("MTW_ttZSel"+Label, MTW, weight, 0., 200., 20);
  }
  if(AN_Sideband){
    FillHist("PTe_ANSideband"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("Etae_ANSideband"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("PTmu1_ANSideband"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
    FillHist("PTmu2_ANSideband"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("Etamu1_ANSideband"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("Etamu2_ANSideband"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
  
    FillHist("Mmumu_ANSideband"+Label, Mmumu, weight, 0., 200., 40);
    FillHist("M3l_ANSideband"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 50);
    FillHist("Nj_ANSideband"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_ANSideband"+Label, BJetColl.size(), weight, 0., 10., 10);
    FillHist("MET_ANSideband"+Label, MET, weight, 0., 300., 30);
    FillHist("MTW_ANSideband"+Label, MTW, weight, 0., 200., 20);
  }
  if(!OnZ){
    FillHist("PTe_OffZ"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("Etae_OffZ"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("PTmu1_OffZ"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
    FillHist("PTmu2_OffZ"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("Etamu1_OffZ"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("Etamu2_OffZ"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
    FillHist("dRemu1_OffZ"+Label, MuTColl.at(0).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
    FillHist("dRemu2_OffZ"+Label, MuTColl.at(1).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
    FillHist("dRmu1mu2_OffZ"+Label, MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
  
    FillHist("Mmumu_OffZ"+Label, Mmumu, weight, 0., 200., 40);
    FillHist("M3l_OffZ"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 100);
    FillHist("Nj_OffZ"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_OffZ"+Label, BJetColl.size(), weight, 0., 10., 10);
    FillHist("MET_OffZ"+Label, MET, weight, 0., 200., 20);
    FillHist("MTW_OffZ"+Label, MTW, weight, 0., 200., 20);
  }
  if(!OnZ && BJetColl.size()==0 && JetColl.size()<3){
    FillHist("PTe_NoBlt2JOffZ"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("Etae_NoBlt2JOffZ"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("PTmu1_NoBlt2JOffZ"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
    FillHist("PTmu2_NoBlt2JOffZ"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("Etamu1_NoBlt2JOffZ"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("Etamu2_NoBlt2JOffZ"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
    FillHist("dRemu1_NoBlt2JOffZ"+Label, MuTColl.at(0).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
    FillHist("dRemu2_NoBlt2JOffZ"+Label, MuTColl.at(1).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
    FillHist("dRmu1mu2_NoBlt2JOffZ"+Label, MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
  
    FillHist("Mmumu_NoBlt2JOffZ"+Label, Mmumu, weight, 0., 200., 40);
    FillHist("M3l_NoBlt2JOffZ"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 100);
    FillHist("Nj_NoBlt2JOffZ"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_NoBlt2JOffZ"+Label, BJetColl.size(), weight, 0., 10., 10);
    FillHist("MET_NoBlt2JOffZ"+Label, MET, weight, 0., 300., 30);
    FillHist("MTW_NoBlt2JOffZ"+Label, MTW, weight, 0., 200., 20);
  }
  if(!OnZ && BJetColl.size()==1 && JetColl.size()<2){
    FillHist("PTe_1Blt2JOffZ"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("Etae_1Blt2JOffZ"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("PTmu1_1Blt2JOffZ"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
    FillHist("PTmu2_1Blt2JOffZ"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("Etamu1_1Blt2JOffZ"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("Etamu2_1Blt2JOffZ"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
    FillHist("dRemu1_1Blt2JOffZ"+Label, MuTColl.at(0).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
    FillHist("dRemu2_1Blt2JOffZ"+Label, MuTColl.at(1).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
    FillHist("dRmu1mu2_1Blt2JOffZ"+Label, MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
  
    FillHist("Mmumu_1Blt2JOffZ"+Label, Mmumu, weight, 0., 200., 40);
    FillHist("M3l_1Blt2JOffZ"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 100);
    FillHist("Nj_1Blt2JOffZ"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_1Blt2JOffZ"+Label, BJetColl.size(), weight, 0., 10., 10);
    FillHist("MET_1Blt2JOffZ"+Label, MET, weight, 0., 300., 30);
    FillHist("MTW_1Blt2JOffZ"+Label, MTW, weight, 0., 200., 20);

  }
  if(!OnZ && BJetColl.size()>0 && JetColl.size()>2 && Mmumu>45){
    FillHist("PTe_ANSBOffZ"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("Etae_ANSBOffZ"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("PTmu1_ANSBOffZ"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
    FillHist("PTmu2_ANSBOffZ"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("Etamu1_ANSBOffZ"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("Etamu2_ANSBOffZ"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
    FillHist("dRemu1_ANSBOffZ"+Label, MuTColl.at(0).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
    FillHist("dRemu2_ANSBOffZ"+Label, MuTColl.at(1).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
    FillHist("dRmu1mu2_ANSBOffZ"+Label, MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
  
    FillHist("Mmumu_ANSBOffZ"+Label, Mmumu, weight, 0., 200., 40);
    FillHist("M3l_ANSBOffZ"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 100);
    FillHist("Nj_ANSBOffZ"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_ANSBOffZ"+Label, BJetColl.size(), weight, 0., 10., 10);
    FillHist("MET_ANSBOffZ"+Label, MET, weight, 0., 300., 30);
    FillHist("MTW_ANSBOffZ"+Label, MTW, weight, 0., 200., 20);
  }
  if( !isData || k_running_nonprompt){
    if(BJetColl.size()>0 && JetColl.size()>2){
      FillHist("PTe_SRincl"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("Etae_SRincl"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("PTmu1_SRincl"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_SRincl"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_SRincl"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_SRincl"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("dRemu1_SRincl"+Label, MuTColl.at(0).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
      FillHist("dRemu2_SRincl"+Label, MuTColl.at(1).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
      FillHist("dRmu1mu2_SRincl"+Label, MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
    
      FillHist("Mmumu_SRincl"+Label, Mmumu, weight, 0., 200., 40);
      FillHist("M3l_SRincl"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 100);
      FillHist("Nj_SRincl"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_SRincl"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_SRincl"+Label, MET, weight, 0., 300., 30);
      FillHist("MTW_SRincl"+Label, MTW, weight, 0., 200., 20);
    }
    if( !OnZ && BJetColl.size()>0 && JetColl.size()>1){
      FillHist("PTe_SR2j"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("Etae_SR2j"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("PTmu1_SR2j"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_SR2j"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_SR2j"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_SR2j"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("dRemu1_SR2j"+Label, MuTColl.at(0).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
      FillHist("dRemu2_SR2j"+Label, MuTColl.at(1).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
      FillHist("dRmu1mu2_SR2j"+Label, MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
    
      FillHist("Mmumu_SR2j"+Label, Mmumu, weight, 0., 200., 40);
      FillHist("M3l_SR2j"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 100);
      FillHist("Nj_SR2j"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_SR2j"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_SR2j"+Label, MET, weight, 0., 300., 30);
      FillHist("MTW_SR2j"+Label, MTW, weight, 0., 200., 20);

      FillHist("Count_SR2j"+Label, 0., weight, 0., 1., 1);
    }
    if( !OnZ && BJetColl.size()>0 && JetColl.size()>2){
      FillHist("PTe_SR3j"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("Etae_SR3j"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("PTmu1_SR3j"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_SR3j"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_SR3j"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_SR3j"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("dRemu1_SR3j"+Label, MuTColl.at(0).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
      FillHist("dRemu2_SR3j"+Label, MuTColl.at(1).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
      FillHist("dRmu1mu2_SR3j"+Label, MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
    
      FillHist("Mmumu_SR3j"+Label, Mmumu, weight, 0., 200., 40);
      FillHist("M3l_SR3j"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 100);
      FillHist("Nj_SR3j"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_SR3j"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_SR3j"+Label, MET, weight, 0., 300., 30);
      FillHist("MTW_SR3j"+Label, MTW, weight, 0., 200., 20);
    }
  }


}


void Oct2017_EleFRSyst::ValidateID(std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, float weight){
  //Purpose: 1) agreement btw data & MC - bias check
  //         2) check efficiency of ID w.r.t. GsfGedEle.
  //*EleLshould be without any ID cuts.


  //Check agreement
  bool PassSel=true; float Mee=0.;
  if( EleTColl.size()!=2 ) PassSel=false;
  if( PassSel && !(EleTColl.at(0).Pt()>25 && EleTColl.at(1).Pt()>15) ) PassSel=false;
  if( PassSel &&   EleTColl.at(0).Charge()== EleTColl.at(1).Charge() ) PassSel=false;
  if( PassSel ) Mee=(EleTColl.at(0)+EleTColl.at(1)).M();
  if( PassSel && fabs(Mee-91.2)>15 ) PassSel=false;
  
  if(PassSel){
    FillHist("PTe1" , EleTColl.at(0).Pt(),  weight, 0., 200., 200);
    FillHist("PTe2" , EleTColl.at(1).Pt(),  weight, 0., 200., 200);
    FillHist("Etae1", EleTColl.at(0).Eta(), weight, -5., 5., 100);
    FillHist("Etae2", EleTColl.at(1).Eta(), weight, -5., 5., 100);
    FillHist("Mee"  , Mee, weight, 60., 120., 60);
  }
  
  //Check Efficiency
  if(isData || k_sample_name.Contains("DY") || k_sample_name.Contains("TT")){
    bool PassLoose=true;
    if( EleLColl.size()!=2 ) PassLoose=false;
    if( PassLoose && EleLColl.at(0).Charge()==EleLColl.at(1).Charge() ) PassLoose=false;
    if( PassLoose && !(EleLColl.at(0).Pt()>25 && EleLColl.at(1).Pt()>15) ) PassLoose=false;
    
    if(PassLoose){
      for(int i=0; i<EleLColl.size(); i++){
        FillHist("NEleSumW_PT" , EleLColl.at(i).Pt(),  weight, 0., 200., 20);
        FillHist("NEleSumW_Eta", EleLColl.at(i).Eta(), weight, -5., 5., 20);
        FillHist("NEleSumW_PTEta", EleLColl.at(i).Pt(), fabs(EleLColl.at(i).Eta()), weight, 0., 200., 20, 0, 2.5, 5);
        if(PassIDCriteria(EleLColl.at(i), "POGMVATIP")){
          FillHist("NEleIDSumW_PT" , EleLColl.at(i).Pt(),  weight, 0., 200., 20);
          FillHist("NEleIDSumW_Eta", EleLColl.at(i).Eta(), weight, -5., 5., 20);
          FillHist("NEleIDSumW_PTEta", EleLColl.at(i).Pt(), fabs(EleLColl.at(i).Eta()), weight, 0., 200., 20, 0., 2.5, 5);
        } 
      }
    }
  }
  

}


bool Oct2017_EleFRSyst::IsNearBJet(snu::KElectron Ele, std::vector<snu::KJet> bjetNoVetoColl){

  bool IsNearB=false;
  for(int i=0; i<bjetNoVetoColl.size(); i++){
    if(Ele.DeltaR(bjetNoVetoColl.at(i))<0.4){
      IsNearB=true; break;
    }
  }

  return IsNearB;
}


float Oct2017_EleFRSyst::ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut){

  //float PTCorr=Ele.Pt()*(1+Ele.PFRelIso(0.3));
  float PTCorr=Ele.Pt()*(1.+max(0.,Ele.PFRelIso(0.3)-TightIsoCut));
  return PTCorr;

}

float Oct2017_EleFRSyst::ConeCorrectedPT(snu::KMuon Mu, float TightIsoCut){

  //float PTCorr=Mu.Pt();
  //float PTCorr=Mu.Pt()*(1.+RochIso(Mu,"0.4"));
  float PTCorr=Mu.Pt()*(1.+max((float) 0.,RochIso(Mu,"0.4")-TightIsoCut));
  return PTCorr;

}



void Oct2017_EleFRSyst::DoSystRun(TString Cycle, TString Mode, std::vector<snu::KElectron> EleColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight){

  int SystDir=0; TString SystKindLabel="", SystDirLabel="";
  if(Mode.Contains("Syst")){
    if(isData&& !k_running_nonprompt) return;
    if     (Mode.Contains("Up"))      {SystDir= 1; SystDirLabel="_systup";}
    else if(Mode.Contains("Down"))    {SystDir=-1; SystDirLabel="_systdown";}
    if     (Mode.Contains("PU"))      SystKindLabel="_PU";
    else if(Mode.Contains("Nvtx"))    SystKindLabel="_Nvtx";
    else if(Mode.Contains("JES"))     SystKindLabel="_JES";
    else if(Mode.Contains("JER"))     SystKindLabel="_JER";
    else if(Mode.Contains("Uncl"))    SystKindLabel="_Uncl";
    else if(Mode.Contains("ElEn"))    SystKindLabel="_ElEn";
    else if(Mode.Contains("MuEn"))    SystKindLabel="_MuEn";
    else if(Mode.Contains("FR"))      SystKindLabel="_FR";
    else if(Mode.Contains("BTag_L"))  SystKindLabel="_BTag_L";
    else if(Mode.Contains("BTag_BC")) SystKindLabel="_BTag_BC";
  }

  if(Cycle=="NormCheck"){

    bool UnPreTrig=false, PreTrig=false;
    if     (Mode.Contains("UnPreTrig")) UnPreTrig=true;
    else if(Mode.Contains("PreTrig"))   PreTrig  =true;

    if     (PreTrig  ) CheckNormCR(EleColl, EleLColl, MuLColl, JetColl, MET, METx, METy, weight, SystDirLabel+SystKindLabel, "PreTrig");
    else if(UnPreTrig) CheckNormCR(EleColl, EleLColl, MuLColl, JetColl, MET, METx, METy, weight, SystDirLabel+SystKindLabel, "UnPreTrig"); 

  }//End of NormCheck
  else if(Cycle=="PromptVarCheck"){
    bool UnPreTrig=false, PreTrig=false;
    if     (Mode.Contains("UnPreTrig")) UnPreTrig=true;
    else if(Mode.Contains("PreTrig"))   PreTrig  =true;

    if(PreTrig) CheckPromptVariation(EleColl, EleLColl, MuLColl, JetColl, BJetColl, MET, METx, METy, weight, SystDirLabel+SystKindLabel, "SiglPreTrig");
  }
  else if(Cycle=="Closure"){
    CheckTrilepCRs(EleColl, EleLColl, MuColl, MuLColl, JetColl, BJetColl, MET, METx, METy, weight, SystDirLabel+SystKindLabel);
  }//End of Closure Cycle

  return;
}


float Oct2017_EleFRSyst::GetFakeWeight(std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleLColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, TString Option){

  float fakeweight=-1.; int NLooseNotTight=0; bool JSTrilepFR=false;
  //if(EleTID=="POGMVAMIP"      && EleLID=="HctoWAFakeLoose") EleFRLabel="POGMVAMTIsop06IPp025p05Sig4FakeLIsop4";
  if(MuTID =="HNTrilepTight2" && MuLID =="HNTrilepFakeL2" ) JSTrilepFR=true;

  TString FilterInfo="", ConeMethod="", SystOpt="";
  if     (Option.Contains("NoFilter"))  FilterInfo="NoFilter";
  else if(Option.Contains("TrkIsoVVL")) FilterInfo="TrkIsoVVL";
  if     (Option.Contains("ConeE"))     ConeMethod="ConeE";
  else if(Option.Contains("ConeSUSY"))  ConeMethod="ConeSUSY";
  else if(Option.Contains("FOPt"))      ConeMethod="FOPt";
  if     (Option.Contains("SystUp"))    SystOpt="SystUp";
  else if(Option.Contains("SystDown"))  SystOpt="SystDown";


  for(int i=0; i<MuLColl.size(); i++){
    if(!PassIDCriteria(MuLColl.at(i), MuTID)){
      float FR=0.;
      if(JSTrilepFR){
//        FR=m_datadriven_bkg->GetFakeObj()->getTrilepFakeRate_muon(false, MuLColl.at(i).Pt(), MuLColl.at(i).Eta());
      }
      else{
        FR=FakeRateData(MuLColl.at(i),MuTID.ReplaceAll("Test_","")+"_"+MuLID.ReplaceAll("Test_","")+"_"+FilterInfo+ConeMethod+SystOpt);
        //cout<<"MuFR"<<FR<<endl;
      }
      fakeweight*=-FR/(1-FR);
      NLooseNotTight++;
    }
  }
  for(int i=0; i<EleLColl.size(); i++){
    if(!PassIDCriteria(EleLColl.at(i), EleTID)){
      float FR=FakeRateData(EleLColl.at(i), "QCD_"+EleTID.ReplaceAll("Test_","")+"_"+EleLID.ReplaceAll("Test_","")+"_"+ConeMethod+SystOpt);
      fakeweight*=-FR/(1-FR);
      NLooseNotTight++;
        //cout<<"ElFR"<<FR<<endl;
    }
  }
  if(NLooseNotTight==0) return 0.;

  return fakeweight;
}


float Oct2017_EleFRSyst::FakeRateData(snu::KMuon Mu, TString Option){

  float FR=0., TightIsoCut=0.;
  int SystDir=0.; bool Syst_FR=Option.Contains("Syst");
  if(Syst_FR){ if(Option.Contains("Up")){ SystDir=1; } else if(Option.Contains("Down")){ SystDir=-1; } }
  if(Option.Contains("ConeSUSY")){
    if     (Option.Contains("TIsop15")) TightIsoCut = 0.15;
    else if(Option.Contains("TIsop20")) TightIsoCut = 0.2;
    else if(Option.Contains("TIsop25")) TightIsoCut = 0.25;
  }
  float PTCorr=ConeCorrectedPT(Mu, TightIsoCut);
    if(Option.Contains("FOPt")) PTCorr=Mu.Pt();
  float fEta=fabs(Mu.Eta());
  if(Option.Contains("POGTIsop20IPp01p1Chi3_POGLIsop6IPp5p1Chi30_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.336779 ;
      else if(PTCorr<20)  FR=0.111232 ;
      else if(PTCorr<25)  FR=0.0902   ;
      else if(PTCorr<35)  FR=0.0752774;
      else                FR=0.0670709;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.381494 ;
      else if(PTCorr<20)  FR=0.137689 ;
      else if(PTCorr<25)  FR=0.120664 ;
      else if(PTCorr<35)  FR=0.106082 ;
      else                FR=0.0923205;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.448765;
      else if(PTCorr<20)  FR=0.174516;
      else if(PTCorr<25)  FR=0.155327;
      else if(PTCorr<35)  FR=0.139882;
      else                FR=0.119441;
    }
    else{
      if     (PTCorr<15)  FR=0.464363;
      else if(PTCorr<20)  FR=0.202699;
      else if(PTCorr<25)  FR=0.17925 ;
      else if(PTCorr<35)  FR=0.167467;
      else                FR=0.138629;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p1Chi3_POGLIsop5IPp5p1Chi_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.336737;
      else if(PTCorr<20)  FR=0.138805;
      else if(PTCorr<25)  FR=0.119428;
      else if(PTCorr<35)  FR=0.102563;
      else                FR=0.095469;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.381256;
      else if(PTCorr<20)  FR=0.167567;
      else if(PTCorr<25)  FR=0.153358;
      else if(PTCorr<35)  FR=0.138543;
      else                FR=0.123659;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.448443;
      else if(PTCorr<20)  FR=0.205231;
      else if(PTCorr<25)  FR=0.191483;
      else if(PTCorr<35)  FR=0.176188;
      else                FR=0.153492;
    }
    else{
      if     (PTCorr<15)  FR=0.464362;
      else if(PTCorr<20)  FR=0.238831;
      else if(PTCorr<25)  FR=0.218972;
      else if(PTCorr<35)  FR=0.209289;
      else                FR=0.177849;
    }
  }
  else if(Option.Contains("POGTIsop15IPp01p1Chi3_POGLIsop6IPp5p1Chi30_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.260288 ;
      else if(PTCorr<20)  FR=0.0712013;
      else if(PTCorr<25)  FR=0.0583859;
      else if(PTCorr<35)  FR=0.0460673;
      else                FR=0.0396087;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.302753 ;
      else if(PTCorr<20)  FR=0.0919691;
      else if(PTCorr<25)  FR=0.0818808;
      else if(PTCorr<35)  FR=0.0667185;
      else                FR=0.0559974;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.362683 ;
      else if(PTCorr<20)  FR=0.122658 ;
      else if(PTCorr<25)  FR=0.104402 ;
      else if(PTCorr<35)  FR=0.0988028;
      else                FR=0.0768322;
    }
    else{
      if     (PTCorr<15)  FR=0.382679 ;
      else if(PTCorr<20)  FR=0.149735 ;
      else if(PTCorr<25)  FR=0.126996 ;
      else if(PTCorr<35)  FR=0.113983 ;
      else                FR=0.0974703;
    }
  }
  else if(Option.Contains("POGTIsop15IPp01p1Chi3_POGLIsop5IPp5p1Chi_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.260264 ;
      else if(PTCorr<20)  FR=0.0887165;
      else if(PTCorr<25)  FR=0.0772104;
      else if(PTCorr<35)  FR=0.0626961;
      else                FR=0.0563453;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.302525 ;
      else if(PTCorr<20)  FR=0.111957 ;
      else if(PTCorr<25)  FR=0.104062 ;
      else if(PTCorr<35)  FR=0.08743  ;
      else                FR=0.0750136;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.362481 ;
      else if(PTCorr<20)  FR=0.14414  ;
      else if(PTCorr<25)  FR=0.128693 ;
      else if(PTCorr<35)  FR=0.124427 ;
      else                FR=0.0987016;
    }
    else{
      if     (PTCorr<15)  FR=0.382678;
      else if(PTCorr<20)  FR=0.17642 ;
      else if(PTCorr<25)  FR=0.155347;
      else if(PTCorr<35)  FR=0.142424;
      else                FR=0.125045;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p1Chi3_POGLIsop6IPp5p1Chi30_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.224901 ;
      else if(PTCorr<20)  FR=0.107681 ;
      else if(PTCorr<25)  FR=0.0898675;
      else if(PTCorr<35)  FR=0.0804815;
      else                FR=0.0727239;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.267794 ;
      else if(PTCorr<20)  FR=0.136115 ;
      else if(PTCorr<25)  FR=0.131556 ;
      else if(PTCorr<35)  FR=0.106586 ;
      else                FR=0.0996845;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.319584;
      else if(PTCorr<20)  FR=0.186026;
      else if(PTCorr<25)  FR=0.163037;
      else if(PTCorr<35)  FR=0.148278;
      else                FR=0.129497;
    }
    else{
      if     (PTCorr<15)  FR=0.344918;
      else if(PTCorr<20)  FR=0.214086;
      else if(PTCorr<25)  FR=0.176514;
      else if(PTCorr<35)  FR=0.184554;
      else                FR=0.154967;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p1Chi3_POGLIsop5IPp5p1Chi_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.243216;
      else if(PTCorr<20)  FR=0.140185;
      else if(PTCorr<25)  FR=0.121813;
      else if(PTCorr<35)  FR=0.111958;
      else                FR=0.106649;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.285605;
      else if(PTCorr<20)  FR=0.171922;
      else if(PTCorr<25)  FR=0.168621;
      else if(PTCorr<35)  FR=0.141603;
      else                FR=0.137314;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.336817;
      else if(PTCorr<20)  FR=0.225045;
      else if(PTCorr<25)  FR=0.206099;
      else if(PTCorr<35)  FR=0.187497;
      else                FR=0.169763;
    }
    else{
      if     (PTCorr<15)  FR=0.363674;
      else if(PTCorr<20)  FR=0.257849;
      else if(PTCorr<25)  FR=0.217374;
      else if(PTCorr<35)  FR=0.234093;
      else                FR=0.201918;
    }
  }
  else if(Option.Contains("POGTIsop15IPp01p1Chi3_POGLIsop6IPp5p1Chi30_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.181152 ;
      else if(PTCorr<20)  FR=0.0707471;
      else if(PTCorr<25)  FR=0.0560844;
      else if(PTCorr<35)  FR=0.0483956;
      else                FR=0.0433495;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.219845 ;
      else if(PTCorr<20)  FR=0.0915236;
      else if(PTCorr<25)  FR=0.0883876;
      else if(PTCorr<35)  FR=0.0658354;
      else                FR=0.0599461;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.267231 ;
      else if(PTCorr<20)  FR=0.129102 ;
      else if(PTCorr<25)  FR=0.110895 ;
      else if(PTCorr<35)  FR=0.103245 ;
      else                FR=0.0819243;
    }
    else{
      if     (PTCorr<15)  FR=0.29299 ;
      else if(PTCorr<20)  FR=0.155378;
      else if(PTCorr<25)  FR=0.124455;
      else if(PTCorr<35)  FR=0.12888 ;
      else                FR=0.104704;
    }
  }
  else if(Option.Contains("POGTIsop15IPp01p1Chi3_POGLIsop5IPp5p1Chi_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.191885 ;
      else if(PTCorr<20)  FR=0.0915866;
      else if(PTCorr<25)  FR=0.0760472;
      else if(PTCorr<35)  FR=0.0671218;
      else                FR=0.0631913;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.23066  ;
      else if(PTCorr<20)  FR=0.115512 ;
      else if(PTCorr<25)  FR=0.113433 ;
      else if(PTCorr<35)  FR=0.0871235;
      else                FR=0.0821759;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.277334;
      else if(PTCorr<20)  FR=0.156263;
      else if(PTCorr<25)  FR=0.140098;
      else if(PTCorr<35)  FR=0.130675;
      else                FR=0.107183;
    }
    else{
      if     (PTCorr<15)  FR=0.305151;
      else if(PTCorr<20)  FR=0.187014;
      else if(PTCorr<25)  FR=0.152749;
      else if(PTCorr<35)  FR=0.164782;
      else                FR=0.136348;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi3_POGLIsop4IPp5p1_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.33842 ;
      else if(PTCorr<20)  FR=0.18468 ;
      else if(PTCorr<25)  FR=0.162466;
      else if(PTCorr<35)  FR=0.143624;
      else                FR=0.137429;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.372154;
      else if(PTCorr<20)  FR=0.210609;
      else if(PTCorr<25)  FR=0.193202;
      else if(PTCorr<35)  FR=0.181434;
      else                FR=0.163875;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.408703;
      else if(PTCorr<20)  FR=0.239243;
      else if(PTCorr<25)  FR=0.224734;
      else if(PTCorr<35)  FR=0.210624;
      else                FR=0.185866;
    }
    else{
      if     (PTCorr<15)  FR=0.397721;
      else if(PTCorr<20)  FR=0.250492;
      else if(PTCorr<25)  FR=0.241592;
      else if(PTCorr<35)  FR=0.228252;
      else                FR=0.20606 ;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi3_POGLIsop6IPp5p1_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.320068 ;
      else if(PTCorr<20)  FR=0.103925 ;
      else if(PTCorr<25)  FR=0.0828112;
      else if(PTCorr<35)  FR=0.0697201;
      else                FR=0.0599762;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.355258 ;
      else if(PTCorr<20)  FR=0.12828  ;
      else if(PTCorr<25)  FR=0.10927  ;
      else if(PTCorr<35)  FR=0.0966748;
      else                FR=0.0815638;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.394842;
      else if(PTCorr<20)  FR=0.155917;
      else if(PTCorr<25)  FR=0.137145;
      else if(PTCorr<35)  FR=0.120733;
      else                FR=0.101367;
    }
    else{
      if     (PTCorr<15)  FR=0.38534 ;
      else if(PTCorr<20)  FR=0.165267;
      else if(PTCorr<25)  FR=0.152769;
      else if(PTCorr<35)  FR=0.136345;
      else                FR=0.114394;
    }
  }
  else if(Option.Contains("HNTrilepTight2_HNTrilepFakeL2_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.314559 ;
      else if(PTCorr<20)  FR=0.115043 ;
      else if(PTCorr<25)  FR=0.0964687;
      else if(PTCorr<35)  FR=0.0742115;
      else                FR=0.0851192;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.35599  ;
      else if(PTCorr<20)  FR=0.142815 ;
      else if(PTCorr<25)  FR=0.129269 ;
      else if(PTCorr<35)  FR=0.0931603;
      else                FR=0.0976207;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.41679 ;
      else if(PTCorr<20)  FR=0.169059;
      else if(PTCorr<25)  FR=0.144527;
      else if(PTCorr<35)  FR=0.1434  ;
      else                FR=0.121598;
    }
    else{
      if     (PTCorr<15)  FR=0.428247;
      else if(PTCorr<20)  FR=0.215242;
      else if(PTCorr<25)  FR=0.188871;
      else if(PTCorr<35)  FR=0.161765;
      else                FR=0.158053;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.276234;
      else if(PTCorr<20)  FR=0.191673;
      else if(PTCorr<25)  FR=0.167833;
      else if(PTCorr<35)  FR=0.162589;
      else                FR=0.159873;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.307693;
      else if(PTCorr<20)  FR=0.220826;
      else if(PTCorr<25)  FR=0.216052;
      else if(PTCorr<35)  FR=0.18657 ;
      else                FR=0.18804 ;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.339593;
      else if(PTCorr<20)  FR=0.262275;
      else if(PTCorr<25)  FR=0.252044;
      else if(PTCorr<35)  FR=0.220044;
      else                FR=0.212143;
    }
    else{
      if     (PTCorr<15)  FR=0.339879;
      else if(PTCorr<20)  FR=0.271265;
      else if(PTCorr<25)  FR=0.249351;
      else if(PTCorr<35)  FR=0.257995;
      else                FR=0.23968 ;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop6IPp5p1_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.216822 ;
      else if(PTCorr<20)  FR=0.101333 ;
      else if(PTCorr<25)  FR=0.0833579;
      else if(PTCorr<35)  FR=0.0747546;
      else                FR=0.0662648;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.251757 ;
      else if(PTCorr<20)  FR=0.128096 ;
      else if(PTCorr<25)  FR=0.120662 ;
      else if(PTCorr<35)  FR=0.0952394;
      else                FR=0.0884953;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.28611 ;
      else if(PTCorr<20)  FR=0.166302;
      else if(PTCorr<25)  FR=0.146939;
      else if(PTCorr<35)  FR=0.126979;
      else                FR=0.110757;
    }
    else{
      if     (PTCorr<15)  FR=0.290039;
      else if(PTCorr<20)  FR=0.175547;
      else if(PTCorr<25)  FR=0.15415 ;
      else if(PTCorr<35)  FR=0.148194;
      else                FR=0.130693;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.344138;
      else if(PTCorr<20)  FR=0.187021;
      else if(PTCorr<25)  FR=0.16477 ;
      else if(PTCorr<35)  FR=0.145246;
      else                FR=0.139742;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.376491;
      else if(PTCorr<20)  FR=0.212285;
      else if(PTCorr<25)  FR=0.19565 ;
      else if(PTCorr<35)  FR=0.182581;
      else                FR=0.165813;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.413423;
      else if(PTCorr<20)  FR=0.242546;
      else if(PTCorr<25)  FR=0.227715;
      else if(PTCorr<35)  FR=0.214613;
      else                FR=0.187678;
    }
    else{
      if     (PTCorr<15)  FR=0.402425;
      else if(PTCorr<20)  FR=0.252608;
      else if(PTCorr<25)  FR=0.24362 ;
      else if(PTCorr<35)  FR=0.229344;
      else                FR=0.208459;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop6IPp5p1_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.325488 ;
      else if(PTCorr<20)  FR=0.10525  ;
      else if(PTCorr<25)  FR=0.0839633;
      else if(PTCorr<35)  FR=0.0704949;
      else                FR=0.0609975;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.35944  ;
      else if(PTCorr<20)  FR=0.129295 ;
      else if(PTCorr<25)  FR=0.110666 ;
      else if(PTCorr<35)  FR=0.0973101;
      else                FR=0.0825229;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.399404;
      else if(PTCorr<20)  FR=0.158063;
      else if(PTCorr<25)  FR=0.138922;
      else if(PTCorr<35)  FR=0.123003;
      else                FR=0.102357;
    }
    else{
      if     (PTCorr<15)  FR=0.389898;
      else if(PTCorr<20)  FR=0.166636;
      else if(PTCorr<25)  FR=0.153958;
      else if(PTCorr<35)  FR=0.137006;
      else                FR=0.115769;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeE")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.344294;
      else if(PTCorr<20)  FR=0.187095;
      else if(PTCorr<25)  FR=0.164859;
      else if(PTCorr<35)  FR=0.145231;
      else                FR=0.14    ;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.376607;
      else if(PTCorr<20)  FR=0.212338;
      else if(PTCorr<25)  FR=0.195648;
      else if(PTCorr<35)  FR=0.182624;
      else                FR=0.166085;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.413485;
      else if(PTCorr<20)  FR=0.242545;
      else if(PTCorr<25)  FR=0.227803;
      else if(PTCorr<35)  FR=0.214703;
      else                FR=0.187801;
    }
    else{
      if     (PTCorr<15)  FR=0.402582;
      else if(PTCorr<20)  FR=0.252604;
      else if(PTCorr<25)  FR=0.243616;
      else if(PTCorr<35)  FR=0.229269;
      else                FR=0.208675;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.2771  ;
      else if(PTCorr<20)  FR=0.193227;
      else if(PTCorr<25)  FR=0.162945;
      else if(PTCorr<35)  FR=0.165255;
      else                FR=0.158028;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.310046;
      else if(PTCorr<20)  FR=0.221337;
      else if(PTCorr<25)  FR=0.214164;
      else if(PTCorr<35)  FR=0.179193;
      else                FR=0.184973;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.340931;
      else if(PTCorr<20)  FR=0.263206;
      else if(PTCorr<25)  FR=0.245025;
      else if(PTCorr<35)  FR=0.222435;
      else                FR=0.214271;
    }
    else{
      if     (PTCorr<15)  FR=0.338846;
      else if(PTCorr<20)  FR=0.269771;
      else if(PTCorr<25)  FR=0.255547;
      else if(PTCorr<35)  FR=0.25313 ;
      else                FR=0.250546;
    }

    if(Syst_FR && SystDir<0){
      if(fEta<0.9){
        if     (PTCorr<15)   FR*=(1.- 0.0655089 );
        else if(PTCorr<20)   FR*=(1.- 0.304225  );
        else if(PTCorr<25)   FR*=(1.- 0.301855  );
        else if(PTCorr<35)   FR*=(1.- 0.387517  );
        else                 FR*=(1.- 0.483798  );
      }
      else if(fEta<1.6){
        if     (PTCorr<15)   FR*=(1.- 0.121366 );
        else if(PTCorr<20)   FR*=(1.- 0.365817 );
        else if(PTCorr<25)   FR*=(1.- 0.310581 );
        else if(PTCorr<35)   FR*=(1.- 0.27695  );
        else                 FR*=(1.- 0.438805 );
      }
      else if(fEta<2.1){
        if     (PTCorr<15)   FR*=(1.- 0.252924 );
        else if(PTCorr<20)   FR*=(1.- 0.368097 );
        else if(PTCorr<25)   FR*=(1.- 0.460861 );
        else if(PTCorr<35)   FR*=(1.- 0.242128 );
        else                 FR*=(1.- 0.415018 );
      }
      else{
        if     (PTCorr<15)   FR*=(1.- 0.243664 );
        else if(PTCorr<20)   FR*=(1.- 0.341107 );
        else if(PTCorr<25)   FR*=(1.- 0.442997 );
        else if(PTCorr<35)   FR*=(1.- 0.292832 );
        else                 FR*=(1.- 0.308315 );
      }
    }
    else if(Syst_FR && SystDir>0){
      if(fEta<0.9){
        if     (PTCorr<15)   FR*=(1.+ 0.0681634 );
        else if(PTCorr<20)   FR*=(1.+ 0.0553126 );
        else if(PTCorr<25)   FR*=(1.+ 0.101094  );
        else if(PTCorr<35)   FR*=(1.+ 0.07793   );
        else                 FR*=(1.+ 0.199967  );
      }
      else if(fEta<1.6){
        if     (PTCorr<15)   FR*=(1.+ 0.0625109);
        else if(PTCorr<20)   FR*=(1.+ 0.0856362);
        else if(PTCorr<25)   FR*=(1.+ 0.0966223);
        else if(PTCorr<35)   FR*=(1.+ 0.108313 );
        else                 FR*=(1.+ 0.146165 );
      }
      else if(fEta<2.1){
        if     (PTCorr<15)   FR*=(1.+ 0.0774913);
        else if(PTCorr<20)   FR*=(1.+ 0.0759516);
        else if(PTCorr<25)   FR*=(1.+ 0.147378 );
        else if(PTCorr<35)   FR*=(1.+ 0.131883 );
        else                 FR*=(1.+ 0.14489  );
      }
      else{
        if     (PTCorr<15)   FR*=(1.+ 0.0810315);
        else if(PTCorr<20)   FR*=(1.+ 0.0809193);
        else if(PTCorr<25)   FR*=(1.+ 0.0970186);
        else if(PTCorr<35)   FR*=(1.+ 0.110969 );
        else                 FR*=(1.+ 0.149472 );
      }
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop6IPp5p1_TrkIsoVVLFOPt")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.155258;
      else if(PTCorr<20)  FR=0.125076;
      else if(PTCorr<25)  FR=0.112006;
      else if(PTCorr<35)  FR=0.114503;
      else                FR=0.120845;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.189011;
      else if(PTCorr<20)  FR=0.158738;
      else if(PTCorr<25)  FR=0.159675;
      else if(PTCorr<35)  FR=0.143399;
      else                FR=0.155353;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.222263;
      else if(PTCorr<20)  FR=0.198951;
      else if(PTCorr<25)  FR=0.190684;
      else if(PTCorr<35)  FR=0.186306;
      else                FR=0.185893;
    }
    else{
      if     (PTCorr<15)  FR=0.227277;
      else if(PTCorr<20)  FR=0.205329;
      else if(PTCorr<25)  FR=0.203358;
      else if(PTCorr<35)  FR=0.223376;
      else                FR=0.222588;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLFOPt")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.247975;
      else if(PTCorr<20)  FR=0.217638;
      else if(PTCorr<25)  FR=0.202237;
      else if(PTCorr<35)  FR=0.195982;
      else                FR=0.214514;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.279385;
      else if(PTCorr<20)  FR=0.243698;
      else if(PTCorr<25)  FR=0.237864;
      else if(PTCorr<35)  FR=0.234583;
      else                FR=0.240188;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.309005;
      else if(PTCorr<20)  FR=0.284354;
      else if(PTCorr<25)  FR=0.274857;
      else if(PTCorr<35)  FR=0.272603;
      else                FR=0.281163;
    }
    else{
      if     (PTCorr<15)  FR=0.313636;
      else if(PTCorr<20)  FR=0.287708;
      else if(PTCorr<25)  FR=0.291454;
      else if(PTCorr<35)  FR=0.297233;
      else                FR=0.308716;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop6IPp5p1sig8_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.254794 ;
      else if(PTCorr<20)  FR=0.121792 ;
      else if(PTCorr<25)  FR=0.0954017;
      else if(PTCorr<35)  FR=0.0937622;
      else                FR=0.0795779;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.289057;
      else if(PTCorr<20)  FR=0.150987;
      else if(PTCorr<25)  FR=0.139373;
      else if(PTCorr<35)  FR=0.109927;
      else                FR=0.104292;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.322044;
      else if(PTCorr<20)  FR=0.18802 ;
      else if(PTCorr<25)  FR=0.163172;
      else if(PTCorr<35)  FR=0.147915;
      else                FR=0.130401;
    }
    else{
      if     (PTCorr<15)  FR=0.312404;
      else if(PTCorr<20)  FR=0.188555;
      else if(PTCorr<25)  FR=0.1747  ;
      else if(PTCorr<35)  FR=0.161571;
      else                FR=0.155695;
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop6IPp5p1sig4_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.306808 ;
      else if(PTCorr<20)  FR=0.148199 ;
      else if(PTCorr<25)  FR=0.116904 ;
      else if(PTCorr<35)  FR=0.116722 ;
      else                FR=0.0994624;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.34711 ;
      else if(PTCorr<20)  FR=0.18568 ;
      else if(PTCorr<25)  FR=0.171709;
      else if(PTCorr<35)  FR=0.135723;
      else                FR=0.130152;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.382281;
      else if(PTCorr<20)  FR=0.226221;
      else if(PTCorr<25)  FR=0.195817;
      else if(PTCorr<35)  FR=0.178244;
      else                FR=0.160737;
    }
    else{
      if     (PTCorr<15)  FR=0.368437;
      else if(PTCorr<20)  FR=0.220804;
      else if(PTCorr<25)  FR=0.21116 ;
      else if(PTCorr<35)  FR=0.188548;
      else                FR=0.187302;
    }
  }

  return FR;
}



float Oct2017_EleFRSyst::FakeRateData(snu::KElectron Ele, TString Option){

  float FR=0., TightIsoCut=0.;
  int SystDir=0.; bool Syst_FR=Option.Contains("Syst");
  if(Syst_FR){ if(Option.Contains("Up")){ SystDir=1; } else if(Option.Contains("Down")){ SystDir=-1; } }
  if(Option.Contains("ConeSUSY")){
    if     (Option.Contains("Isop06")) TightIsoCut = 0.06;
    else if(Option.Contains("Isop08")) TightIsoCut = 0.08;
    else if(Option.Contains("Isop1"))  TightIsoCut = 0.1;
  }
  float PTCorr=ConeCorrectedPT(Ele, TightIsoCut);
    if(Option.Contains("FOPt")) PTCorr=Ele.Pt();
  float fEta=fabs(Ele.Eta());


  if(Option.Contains("QCD_POGWP90Isop06IPp025p05sig4_LMVA06v1Isop4IPp025p05sig4_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.110603 ;
      else if(PTCorr<50)  FR=0.0644239;
      else                FR=0.0686986;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.136154 ;
      else if(PTCorr<50)  FR=0.0718016;
      else                FR=0.0789876;
    }
    else{
      if     (PTCorr<35)  FR=0.184978;
      else if(PTCorr<50)  FR=0.10697 ;
      else                FR=0.11768 ;
    }
  }
  else if(Option.Contains("QCD_POGWP90Isop06IPp025p05sig4_LMVA06Isop4IPp025p05sig4_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.118837 ;
      else if(PTCorr<50)  FR=0.0626265;
      else                FR=0.0649617;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.129257 ;
      else if(PTCorr<50)  FR=0.0681096;
      else                FR=0.0732587;
    }
    else{
      if     (PTCorr<35)  FR=0.167549;
      else if(PTCorr<50)  FR=0.1026  ;
      else                FR=0.119577;
    }

    if(Syst_FR && SystDir<0){
      if(fEta<0.8){
        if     (PTCorr<35)   FR*=(1.- 0.0791571 );
        else if(PTCorr<50)   FR*=(1.- 0.391781  );
        else                 FR*=(1.- 0.477299  );
      }
      else if(fEta<1.479){
        if     (PTCorr<35)   FR*=(1.- 0.160179  );
        else if(PTCorr<50)   FR*=(1.- 0.510462  );
        else                 FR*=(1.- 0.23345   );
      }
      else{
        if     (PTCorr<35)   FR*=(1.- 0.0610817 );
        else if(PTCorr<50)   FR*=(1.- 0.112503  );
        else                 FR*=(1.- 0.114225  );
      }
    }
    else if(Syst_FR && SystDir>0){
      if(fEta<0.8){
        if     (PTCorr<35)   FR*=(1.+ 0.122614 );
        else if(PTCorr<50)   FR*=(1.+ 0.190562 );
        else                 FR*=(1.+ 0.282318 );
      }
      else if(fEta<1.479){
        if     (PTCorr<35)   FR*=(1.+ 0.10469  );
        else if(PTCorr<50)   FR*=(1.+ 0.137677 );
        else                 FR*=(1.+ 0.208225 );
      }
      else{
        if     (PTCorr<35)   FR*=(1.+ 0.156983  );
        else if(PTCorr<50)   FR*=(1.+ 0.0675837 );
        else                 FR*=(1.+ 0.05555   );
      }
    }
  }
  else if(Option.Contains("QCD_POGWP90Isop06IPp025p05sig4_LMVA06v1Isop4IPp5p1_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.0736732;
      else if(PTCorr<50)  FR=0.0302377;
      else                FR=0.0339471;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.101205 ;
      else if(PTCorr<50)  FR=0.0460958;
      else                FR=0.0478649;
    }
    else{
      if     (PTCorr<35)  FR=0.11822  ;
      else if(PTCorr<50)  FR=0.0659004;
      else                FR=0.0723459;
    }
  }
  else if(Option.Contains("QCD_POGWP90Isop06IPp025p05sig4_LMVA06v2Isop4IPp025p05sig4_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.0974724;
      else if(PTCorr<50)  FR=0.0556398;
      else                FR=0.0571567;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.122288 ;
      else if(PTCorr<50)  FR=0.0634097;
      else                FR=0.0679921;
    }
    else{
      if     (PTCorr<35)  FR=0.166618 ;
      else if(PTCorr<50)  FR=0.0943797;
      else                FR=0.10212  ;
    }
  }
  else if(Option.Contains("QCD_POGWP90Isop08IPp025p05sig4_LMVA08v1Isop4IPp025p05sig4_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.154559 ;
      else if(PTCorr<50)  FR=0.0924185;
      else                FR=0.092448 ;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.193706;
      else if(PTCorr<50)  FR=0.107652;
      else                FR=0.11214 ;
    }
    else{
      if     (PTCorr<35)  FR=0.257529;
      else if(PTCorr<50)  FR=0.158485;
      else                FR=0.175704;
    }
  }
  else if(Option.Contains("QCD_POGWP90Isop08IPp025p05sig4_LMVA08v2Isop4IPp025p05sig4_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.137626 ;
      else if(PTCorr<50)  FR=0.0804705;
      else                FR=0.0782686;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.176499 ;
      else if(PTCorr<50)  FR=0.0959684;
      else                FR=0.0986659;
    }
    else{
      if     (PTCorr<35)  FR=0.233378;
      else if(PTCorr<50)  FR=0.14105 ;
      else                FR=0.154777;
    }
  }
  else if(Option.Contains("QCD_POGWP90Isop1IPp025p05sig4_LMVA1v1Isop4IPp025p05sig4_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.209062;
      else if(PTCorr<50)  FR=0.12362 ;
      else                FR=0.122594;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.254915;
      else if(PTCorr<50)  FR=0.145394;
      else                FR=0.149414;
    }
    else{
      if     (PTCorr<35)  FR=0.331213;
      else if(PTCorr<50)  FR=0.214573;
      else                FR=0.236407;
    }
  }
  else if(Option.Contains("QCD_POGWP90Isop1IPp025p05sig4_LMVA1v2Isop4IPp025p05sig4_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.185823;
      else if(PTCorr<50)  FR=0.107285;
      else                FR=0.103929;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.237566;
      else if(PTCorr<50)  FR=0.133319;
      else                FR=0.135416;
    }
    else{
      if     (PTCorr<35)  FR=0.306354;
      else if(PTCorr<50)  FR=0.19597 ;
      else                FR=0.213273;
    }
  }
  else if(Option.Contains("POGMVAMIPOpt2Iso04")){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.193686 ;
      else if(PTCorr<50 ) FR= 0.115564 ;
      else                FR= 0.11773  ;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.233713;
      else if(PTCorr<50 ) FR= 0.127104;
      else                FR= 0.139068;
    }
    else{
      if     (PTCorr<35 ) FR= 0.261925;
      else if(PTCorr<50 ) FR= 0.168822;
      else                FR= 0.193664;
    }
  }
  else if(Option.Contains("POGMVAMTIsop06IPp025p05Sig4FakeLIsop4")){
    //OLD(MVA:(-0.9,-0.9,-0.8) ; Coarse Optimisation
//    if(fEta<0.8){
//      if     (PTCorr<35 ) FR= 0.122623 ;
//      else if(PTCorr<50 ) FR= 0.0679368;
//      else                FR= 0.0707609;
//    }
//    else if(fEta<1.479){
//      if     (PTCorr<35 ) FR= 0.12563  ;
//      else if(PTCorr<50 ) FR= 0.0646067;
//      else                FR= 0.0715022;
//    }
//    else{
//      if     (PTCorr<35 ) FR= 0.156385 ;
//      else if(PTCorr<50 ) FR= 0.0942703;
//      else                FR= 0.104913 ;
//    }
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.119907 ;
      else if(PTCorr<50 ) FR= 0.065943 ;
      else                FR= 0.0682414;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.132106 ;
      else if(PTCorr<50 ) FR= 0.0682482;
      else                FR= 0.077325 ;
    }
    else{
      if     (PTCorr<35 ) FR= 0.160416 ;
      else if(PTCorr<50 ) FR= 0.0976898;
      else                FR= 0.109458 ;
    }
  }
  else if(Option.Contains("POGMVAMTIsop1IPp025p05Sig4FakeLIsop4")){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.204863;
      else if(PTCorr<50 ) FR= 0.120178;
      else                FR= 0.123658;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.236692;
      else if(PTCorr<50 ) FR= 0.128829;
      else                FR= 0.141241;
    }
    else{
      if     (PTCorr<35 ) FR= 0.294857;
      else if(PTCorr<50 ) FR= 0.193047;
      else                FR= 0.219773;
    }
  }
  else if(Option.Contains("POGMVATTIsop06IPp025p05Sig4FakeLIsop4")){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.0813384;
      else if(PTCorr<50 ) FR= 0.0473129;
      else                FR= 0.0490471;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.0949602;
      else if(PTCorr<50 ) FR= 0.0517796;
      else                FR= 0.0630693;
    }
    else{
      if     (PTCorr<35 ) FR= 0.109435 ;
      else if(PTCorr<50 ) FR= 0.0670226;
      else                FR= 0.079485 ;
    }
  }
  else if(Option.Contains("POGMVATTIsop1IPp025p05Sig4FakeLIsop4")){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.146427 ;
      else if(PTCorr<50 ) FR= 0.0908473;
      else                FR= 0.0960944;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.18927 ;
      else if(PTCorr<50 ) FR= 0.108929;
      else                FR= 0.126564;
    }
    else{
      if     (PTCorr<35 ) FR= 0.217664;
      else if(PTCorr<50 ) FR= 0.141156;
      else                FR= 0.173934;
    }
  }
  else if(Option.Contains("POGMVAMIPOpt2Iso06")){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.193733;
      else if(PTCorr<50 ) FR= 0.103705;
      else                FR= 0.104228;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.233731;
      else if(PTCorr<50 ) FR= 0.115073;
      else                FR= 0.12269 ;
    }
    else{
      if     (PTCorr<35 ) FR= 0.261928;
      else if(PTCorr<50 ) FR= 0.15727 ;
      else                FR= 0.179709;
    }
  }
  else if(Option.Contains("POGMVAMFakeLIso04NoMVA")){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.134663 ;
      else if(PTCorr<50 ) FR= 0.0726663;
      else                FR= 0.0687263;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.139805 ;
      else if(PTCorr<50 ) FR= 0.0696455;
      else                FR= 0.0688009;
    }
    else{
      if     (PTCorr<35 ) FR= 0.0977045;
      else if(PTCorr<50 ) FR= 0.0474763;
      else                FR= 0.0413024;
    }
  }
  else if(Option.Contains("POGMVAMFakeLIso04Opt2Data")){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.169891 ;
      else if(PTCorr<50 ) FR= 0.0980297;
      else                FR= 0.0972565;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.196989;
      else if(PTCorr<50 ) FR= 0.103387;
      else                FR= 0.109215;
    }
    else{
      if     (PTCorr<35 ) FR= 0.1734  ;
      else if(PTCorr<50 ) FR= 0.102403;
      else                FR= 0.108606;
    }
  }
  else if(Option.Contains("POGMVAMFakeLIso04Opt1")){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.232943;
      else if(PTCorr<50 ) FR= 0.143792;
      else                FR= 0.154709;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.277706;
      else if(PTCorr<50 ) FR= 0.156332;
      else                FR= 0.177498;
    }
    else{
      if     (PTCorr<35 ) FR= 0.336682;
      else if(PTCorr<50 ) FR= 0.226226;
      else                FR= 0.266821;
    }
  }
  return FR;
} 


float Oct2017_EleFRSyst::GetPreTrigPURW(int Nvtx){

   float weight=1.;
   if(Nvtx<2) weight=2.44008;
   else if(Nvtx<3) weight=3.766;
   else if(Nvtx<4) weight=4.23426;
   else if(Nvtx<5) weight=3.07165;
   else if(Nvtx<6) weight=3.80078;
   else if(Nvtx<7) weight=3.53227;
   else if(Nvtx<8) weight=3.11202;
   else if(Nvtx<9) weight=2.95288;
   else if(Nvtx<10) weight=2.49038;
   else if(Nvtx<11) weight=2.32367;
   else if(Nvtx<12) weight=2.15409;
   else if(Nvtx<13) weight=1.88366;
   else if(Nvtx<14) weight=1.68599;
   else if(Nvtx<15) weight=1.50876;
   else if(Nvtx<16) weight=1.31791;
   else if(Nvtx<17) weight=1.15175;
   else if(Nvtx<18) weight=0.995099;
   else if(Nvtx<19) weight=0.869139;
   else if(Nvtx<20) weight=0.746565;
   else if(Nvtx<21) weight=0.642937;
   else if(Nvtx<22) weight=0.575486;
   else if(Nvtx<23) weight=0.509482;
   else if(Nvtx<24) weight=0.451994;
   else if(Nvtx<25) weight=0.399878;
   else if(Nvtx<26) weight=0.373236;
   else if(Nvtx<27) weight=0.339509;
   else if(Nvtx<28) weight=0.310082;
   else if(Nvtx<29) weight=0.326138;
   else if(Nvtx<30) weight=0.262004;
   else if(Nvtx<31) weight=0.272648;
   else if(Nvtx<32) weight=0.224301;
   else if(Nvtx<33) weight=0.239823;
   else if(Nvtx<34) weight=0.225261;
   else if(Nvtx<35) weight=0.226247;
   else if(Nvtx<36) weight=0.316148;
   else if(Nvtx<37) weight=0.272005;
   else if(Nvtx<38) weight=0.209038;
   else if(Nvtx<39) weight=0.296724;
   else if(Nvtx<40) weight=0.323122;
   else if(Nvtx<41) weight=0.243641;
   else if(Nvtx<42) weight=0.242447;
   else if(Nvtx<43) weight=0.305836;
   else if(Nvtx<44) weight=0.708821;
   else if(Nvtx<45) weight=0.29469;
   else if(Nvtx<46) weight=0.502953;
   else if(Nvtx<47) weight=0.666761;
   else if(Nvtx<48) weight=0.423836;
   else if(Nvtx<49) weight=1.06903;
   else if(Nvtx<50) weight=1.74698;
   else weight=0.890914;

//normed at gt1jmet30
//   if(Nvtx<2) weight=1.22087;
//   else if(Nvtx<3) weight=4.79564;
//   else if(Nvtx<4) weight=3.84689;
//   else if(Nvtx<5) weight=2.47179;
//   else if(Nvtx<6) weight=3.74121;
//   else if(Nvtx<7) weight=3.16705;
//   else if(Nvtx<8) weight=3.04533;
//   else if(Nvtx<9) weight=2.85291;
//   else if(Nvtx<10) weight=2.50907;
//   else if(Nvtx<11) weight=2.29848;
//   else if(Nvtx<12) weight=2.25329;
//   else if(Nvtx<13) weight=1.9864;
//   else if(Nvtx<14) weight=1.73597;
//   else if(Nvtx<15) weight=1.57862;
//   else if(Nvtx<16) weight=1.38849;
//   else if(Nvtx<17) weight=1.21789;
//   else if(Nvtx<18) weight=1.05361;
//   else if(Nvtx<19) weight=0.909483;
//   else if(Nvtx<20) weight=0.787018;
//   else if(Nvtx<21) weight=0.699577;
//   else if(Nvtx<22) weight=0.606739;
//   else if(Nvtx<23) weight=0.533779;
//   else if(Nvtx<24) weight=0.47826;
//   else if(Nvtx<25) weight=0.430963;
//   else if(Nvtx<26) weight=0.38291;
//   else if(Nvtx<27) weight=0.368938;
//   else if(Nvtx<28) weight=0.306386;
//   else if(Nvtx<29) weight=0.321641;
//   else if(Nvtx<30) weight=0.267767;
//   else if(Nvtx<31) weight=0.303042;
//   else if(Nvtx<32) weight=0.220255;
//   else if(Nvtx<33) weight=0.237769;
//   else if(Nvtx<34) weight=0.214341;
//   else if(Nvtx<35) weight=0.254868;
//   else if(Nvtx<36) weight=0.318919;
//   else if(Nvtx<37) weight=0.290442;
//   else if(Nvtx<38) weight=0.204838;
//   else if(Nvtx<39) weight=0.389061;
//   else if(Nvtx<40) weight=0.309147;
//   else if(Nvtx<41) weight=0.286387;
//   else if(Nvtx<42) weight=0.266532;
//   else if(Nvtx<43) weight=0.327612;
//   else if(Nvtx<44) weight=0.575947;
//   else if(Nvtx<45) weight=0.23432;
//   else if(Nvtx<46) weight=0.572198;
//   else if(Nvtx<47) weight=0.541966;
//   else if(Nvtx<48) weight=0.384433;
//   else if(Nvtx<49) weight=1.57278;
//   else weight=0.389116;

   return weight;
}



void Oct2017_EleFRSyst::FillCutFlow(TString cut, float weight){
  
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



void Oct2017_EleFRSyst::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Oct2017_EleFRSyst::MakeHistograms(){
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
  // *  Remove//Overide this Oct2017_EleFRSystCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Oct2017_EleFRSyst::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
/*
  if(MuTColl.at(0).Pt()>130. && MuTColl.at(0).Pt()<140.){
    if(!HasBJet){
      FillHist("IsoEl_Pt130140_3lOS", EleTColl.at(0).PFRelIso(0.3), weight, 0., 0.4, 40);
      FillHist("d0SigEl_Pt130140_3lOS", fabs(EleTColl.at(0).dxySig()), weight, 0., 8., 40);
      FillHist("dzEl_Pt130140_3lOS", fabs(EleTColl.at(0).dz()), weight, 0., 0.1, 100);

      FillHist("Isomu1_Pt130140_3lOS", MuTColl.at(0).RelIso04(), weight, 0., 0.4, 40);
      FillHist("d0Sigmu1_Pt130140_3lOS", fabs(MuTColl.at(0).dXYSig()), weight, 0., 8., 40);
      FillHist("dzmu1_Pt130140_3lOS", fabs(MuTColl.at(0).dZ()), weight, 0., 0.1, 100);

      FillHist("Isomu2_Pt130140_3lOS", MuTColl.at(1).RelIso04(), weight, 0., 0.4, 40);
      FillHist("d0Sigmu2_Pt130140_3lOS", fabs(MuTColl.at(1).dXYSig()), weight, 0., 8., 40);
      FillHist("dzmu2_Pt130140_3lOS", fabs(MuTColl.at(1).dZ()), weight, 0., 0.1, 100);


      FillHist("dRemu1_Pt130140_3lOS", EleTColl.at(0).DeltaR(MuTColl.at(0)), weight, 0., 5., 50);
      FillHist("dRemu2_Pt130140_3lOS", EleTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
      FillHist("dRmu1mu2_Pt130140_3lOS", MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
  
      FillHist("Pte_Pt130140_3lOS", EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("Ptmu2_Pt130140_3lOS", MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("Etae_Pt130140_3lOS", EleTColl.at(0).Eta(), weight, -2.5, 2.5, 10);
      FillHist("Etamu1_Pt130140_3lOS", MuTColl.at(0).Eta(), weight, -2.5, 2.5, 10);
      FillHist("Etamu2_Pt130140_3lOS", MuTColl.at(1).Eta(), weight, -2.5, 2.5, 10);
      FillHist("Mmumu_Pt130140_3lOS", (MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 200., 40);
      FillHist("Memu1_Pt130140_3lOS", (EleTColl.at(0)+MuTColl.at(0)).M(), weight, 0., 200., 40);
      FillHist("Memu2_Pt130140_3lOS", (EleTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 200., 40);
      FillHist("MET_Pt130140_3lOS", MET, weight, 0., 200., 20);
      FillHist("MTW_Pt130140_3lOS", MTW, weight, 0., 200., 20);
  
      FillHist("Nj_Pt130140_3lOS", JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_Pt130140_3lOS", BJetColl.size(), weight, 0., 10., 10);
    }
    if(ttZSel){
      FillHist("IsoEl_Pt130140_ttZSel", EleTColl.at(0).PFRelIso(0.3), weight, 0., 0.4, 40);
      FillHist("d0SigEl_Pt130140_ttZSel", fabs(EleTColl.at(0).dxySig()), weight, 0., 8., 40);
      FillHist("dzEl_Pt130140_ttZSel", fabs(EleTColl.at(0).dz()), weight, 0., 0.1, 100);

      FillHist("Isomu1_Pt130140_ttZSel", MuTColl.at(0).RelIso04(), weight, 0., 0.4, 40);
      FillHist("d0Sigmu1_Pt130140_ttZSel", fabs(MuTColl.at(0).dXYSig()), weight, 0., 8., 40);
      FillHist("dzmu1_Pt130140_ttZSel", fabs(MuTColl.at(0).dZ()), weight, 0., 0.1, 100);

      FillHist("Isomu2_Pt130140_ttZSel", MuTColl.at(1).RelIso04(), weight, 0., 0.4, 40);
      FillHist("d0Sigmu2_Pt130140_ttZSel", fabs(MuTColl.at(1).dXYSig()), weight, 0., 8., 40);
      FillHist("dzmu2_Pt130140_ttZSel", fabs(MuTColl.at(1).dZ()), weight, 0., 0.1, 100);

      FillHist("dRemu1_Pt130140_ttZSel", EleTColl.at(0).DeltaR(MuTColl.at(0)), weight, 0., 5., 50);
      FillHist("dRemu2_Pt130140_ttZSel", EleTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
      FillHist("dRmu1mu2_Pt130140_ttZSel", MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
  
      FillHist("Pte_Pt130140_ttZSel", EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("Ptmu2_Pt130140_ttZSel", MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("Etae_Pt130140_ttZSel", EleTColl.at(0).Eta(), weight, -2.5, 2.5, 10);
      FillHist("Etamu1_Pt130140_ttZSel", MuTColl.at(0).Eta(), weight, -2.5, 2.5, 10);
      FillHist("Etamu2_Pt130140_ttZSel", MuTColl.at(1).Eta(), weight, -2.5, 2.5, 10);
      FillHist("MET_Pt130140_ttZSel", MET, weight, 0., 200., 20);
      FillHist("MTW_Pt130140_ttZSel", MTW, weight, 0., 200., 20);
  
      FillHist("Nj_Pt130140_ttZSel", JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_Pt130140_ttZSel", BJetColl.size(), weight, 0., 10., 10);
      FillHist("Mmumu_Pt130140_ttZSel", (MuTColl.at(0)+MuTColl.at(1)).M(), weight, 60., 120., 30);
    }
  }
*/
