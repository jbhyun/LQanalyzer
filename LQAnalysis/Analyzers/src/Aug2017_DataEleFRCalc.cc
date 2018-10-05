// $Id: Aug2017_DataEleFRCalc.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQAug2017_DataEleFRCalc Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Aug2017_DataEleFRCalc.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Aug2017_DataEleFRCalc);

 Aug2017_DataEleFRCalc::Aug2017_DataEleFRCalc() : AnalyzerCore(), out_muons(0) {

   SetLogName("Aug2017_DataEleFRCalc");
   Message("In Aug2017_DataEleFRCalc constructor", INFO);
   InitialiseAnalysis();
 }


 void Aug2017_DataEleFRCalc::InitialiseAnalysis() throw( LQError ) {
   
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

void Aug2017_DataEleFRCalc::ExecuteEvents()throw( LQError ){

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
   bool SSDilepCR=false;
   bool DoubleMuon=false, ElectronMuon=false;
   bool HighdXYFRCheck=false, IsoIPOpt=false;
   for(int i=0; i<(int) k_flags.size(); i++){
     if     (k_flags.at(i).Contains("NormCheck"))      NormCheck      = true;
     else if(k_flags.at(i).Contains("FRMeasure"))      FRMeasure      = true;
     else if(k_flags.at(i).Contains("SiglWP"))         SiglWP         = true;
     else if(k_flags.at(i).Contains("FRScan"))         FRScan         = true;
     else if(k_flags.at(i).Contains("UnPreTrig"))      UnPreTrig      = true;
     else if(k_flags.at(i).Contains("SiglPreTrig"))    SiglPreTrig    = true;
     else if(k_flags.at(i).Contains("MultPreTrig"))    MultPreTrig    = true;
     else if(k_flags.at(i).Contains("OptMETMTWCuts"))  OptMETMTWCuts  = true;
     else if(k_flags.at(i).Contains("Closure"))        Closure        = true;
     else if(k_flags.at(i).Contains("DoubleMuon"))     DoubleMuon     = true;
     else if(k_flags.at(i).Contains("ElectronMuon"))   ElectronMuon   = true;
     else if(k_flags.at(i).Contains("HighdXYFRCheck")) HighdXYFRCheck = true;
     else if(k_flags.at(i).Contains("IsoIPOpt"))       IsoIPOpt       = true;
     else if(k_flags.at(i).Contains("IDValidation"))   IDValidation   = true;
     else if(k_flags.at(i).Contains("SystRun"))        SystRun        = true;
     else if(k_flags.at(i).Contains("SSDilepCR"))      SSDilepCR      = true;
   }

   if(SiglPreTrig || MultPreTrig) PreTrig=true;
    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //ListTriggersAvailable();

   //Trigger Path of Analysis
   bool  Pass_Trigger=false, Pass_TriggerBG=false, Pass_TriggerH=false;
   float trigger_ps_weight=1., trigger_period_weight=1.;
   float LumiBG=27.257618, LumiH=8.605696, LumiBH=35.863314;
   if(FRMeasure || NormCheck || OptMETMTWCuts || HighdXYFRCheck){
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
       std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);

       int Case=0;
       if( electronLooseColl.size()==1 && electronLooseColl.at(0).Pt()>25
           && PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ) {Pass_Trigger=true; Case=1;}
       if( electronLooseColl.size()==1 && electronLooseColl.at(0).Pt()<25 && electronLooseColl.at(0).Pt()>20 
           && PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ) {Pass_Trigger=true; Case=2;}
       if( electronLooseColl.size()==1 && electronLooseColl.at(0).Pt()<20 && electronLooseColl.at(0).Pt()>15 
           && PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ) {Pass_Trigger=true; Case=3;}

       if     (!isData && Case==1) trigger_ps_weight=WeightByTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v", TargetLumi);
       else if(!isData && Case==2) trigger_ps_weight=WeightByTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v", TargetLumi);
       else if(!isData && Case==3) trigger_ps_weight=WeightByTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v", TargetLumi);
     }

   }
   if(Closure || SSDilepCR){
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") )    Pass_TriggerBG=true;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_TriggerH =true;

     if(isData){
       int DataPeriod=GetDataPeriod();
       if( DataPeriod>0 && DataPeriod<7 && Pass_TriggerBG ) Pass_Trigger=true;
       else if( DataPeriod==7 && Pass_TriggerH ) Pass_Trigger=true;
     }
     else{
       if( Pass_TriggerBG || Pass_TriggerH ) Pass_Trigger=true;
       trigger_period_weight=( (Pass_TriggerBG? 27.257618:0.)+(Pass_TriggerH? 8.605696:0.) )/35.863314;
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
     eventbase->GetMuonSel()->SetPt(5.);                     eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, true);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
     if(FRMeasure || NormCheck || HighdXYFRCheck ){ if( !(electronPreColl.size()>=1) ) return; }
     else if(Closure)                   { if( !(electronPreColl.size()>=1 && muonPreColl.size()>=2) ) return; }
     else if(IsoIPOpt)                  { if( !(electronPreColl.size()>=2) ) return; }
     else if(SSDilepCR){ bool Pass=false;
       int Nlep=muonPreColl.size()+electronPreColl.size();
       if(Nlep>=3) Pass=true;
       else if(Nlep==2){ int SumCharge=0; for(int i=0; i<(int) muonPreColl.size(); i++){ SumCharge+=muonPreColl.at(i).Charge();}
                                          for(int i=0; i<(int) electronPreColl.size(); i++){ SumCharge+=electronPreColl.at(i).Charge(); }  
                         if(SumCharge!=0) Pass=true; }
       if(!Pass) return;
     }
     FillCutFlow("PreSel", weight);
   /**********************************************************************************************************/

     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.2);                 eventbase->GetMuonSel()->SetdxySigMax(4.);
     eventbase->GetMuonSel()->SetBSdz(0.1);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.6);
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);

     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(4.);
     eventbase->GetMuonSel()->SetBSdz(0.05);
     eventbase->GetMuonSel()->SetChiNdof(4.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.2);
   std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,true);
   std::vector<snu::KMuon> muonColl;      if(k_running_nonprompt){ muonColl=muonLooseColl;} else{ muonColl=muonTightColl;}

     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronFakeLPreOptColl; eventbase->GetElectronSel()->Selection(electronFakeLPreOptColl);


     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_HctoWA_FAKELOOSE);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);     eventbase->GetElectronSel()->SetSCEta(2.5);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.4, 0.4);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);     eventbase->GetElectronSel()->SetSCEta(2.5);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronTightColl; eventbase->GetElectronSel()->Selection(electronTightColl);
   std::vector<snu::KElectron> electronColl;
     if(!k_running_nonprompt){ electronColl=electronTightColl; }else{ electronColl=electronLooseColl; }

     bool LeptonVeto=false;
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
   std::vector<snu::KJet> jetNoVetoColl; eventbase->GetJetSel()->Selection(jetNoVetoColl, LeptonVeto, muonLooseColl, electronLooseColl);
   std::vector<snu::KJet> bjetNoVetoColl = SelBJets(jetNoVetoColl, "Medium");

   std::vector<snu::KJet> jetColl  = SkimJetColl(jetNoVetoColl,  electronLooseColl, muonLooseColl, "EleMuVeto");
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
   float trignorm_sf=1., nvtx_reweight=1.; 

   //This part is for boosting up speed.. SF part takes rather longer time than expected
   //Systematic study case should be treated in other part, since there are so many diff. kinds of variations affecting each other.
   if(!SystRun){
     if     (Closure     ){ if(electronLooseColl.size()>=1 && muonLooseColl.size()>=2) EventCand=true; }
     else if(IDValidation){ if(electronLooseColl.size()>=2) EventCand=true;}
     else if(NormCheck || OptMETMTWCuts || FRMeasure ){ if(electronLooseColl.size()>=1) EventCand=true;}
     else if(SSDilepCR   ){ if(muonLooseColl.size()>=1 && electronLooseColl.size()>=1) EventCand=true; }

     if(!isData){
       if(EventCand && (Closure || OptMETMTWCuts || IDValidation || SSDilepCR )){
         //geneff_weight   = GenFilterEfficiency(k_sample_name);
         //gennorm_weight  = SignalNorm(k_sample_name, 200.);
         k_factor_weight = GetKFactor();
  
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl);
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
    
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl);
         trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);

         if(NormCheck && PreTrig){
           trigger_sf = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
         }
         if(DoubleMuon){
           float trigger_sf1 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v");
           float trigger_sf2 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v");
           trigger_sf    = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
         }
         if(ElectronMuon){
           float trigger_sf1 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
           float trigger_sf2 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
           trigger_sf    = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
         }

         if(PreTrig && !NormCheck)   nvtx_reweight = GetPreTrigPURW(Nvtx);

         if(Closure || SSDilepCR){
           btag_sf       = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);}
       }
       if(HighdXYFRCheck && (k_sample_name.Contains("qcd") ||k_sample_name.Contains("QCD"))) pileup_reweight=1.;
     }
     else if(k_running_nonprompt && EventCand && Closure ){
       //fake_weight = GetFakeWeight(muonLooseColl, electronFakeLColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", "LMVA06Isop4IPp025p05sig4", "POGWP90Isop06IPp025p05sig4", "TrkIsoVVLConeSUSY");
     }
   }

   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight*fake_weight*btag_sf*trigger_sf;
   weight *= nvtx_reweight;
   weight *= geneff_weight*gennorm_weight*k_factor_weight;
   //-----------------------------------------------------------------------------------------//


//----------------------------------------------------------------------------------//
//////////////////////////////////////////////////////////////////////////////////////
/////// MAIN ANALYSIS CODE ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//----------------------------------------------------------------------------------//


   if(NormCheck){
     
     if(!SystRun){
       if(PreTrig){
        float nvtx_reweight_UnPreTrig=1., nvtx_reweight_PreTrig=1.;
        if(!isData){
           nvtx_reweight_PreTrig = GetPreTrigPURW(Nvtx);
        }

        CheckNormCR(electronColl, electronLooseColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "PreTrig");
        CheckNormCR(electronColl, electronLooseColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight*nvtx_reweight_PreTrig, "_NvtxRWPreTrig", "PreTrig");
       }
     }//End of not SystRun
     if(SystRun){

       if(isData){
         if(PreTrig){
           DoSystRun("NormCheck", "PreTrig",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     1.);
         }
         if(UnPreTrig){
           DoSystRun("NormCheck", "UnPreTrig",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     1.);
         }
       }
       else{
           eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
           eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
           eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
           eventbase->GetElectronSel()->SetBETrRegIncl(false);     eventbase->GetElectronSel()->SetSCEta(2.5);
           eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
           eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
           eventbase->GetElectronSel()->SetdxySigMax(4.);
           eventbase->GetElectronSel()->SetApplyConvVeto(true);
         std::vector<snu::KElectron> EleTElEnUpColl; eventbase->GetElectronSel()->Selection(EleTElEnUpColl, "SystUpElEn");
  
           eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
           eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
           eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
           eventbase->GetElectronSel()->SetBETrRegIncl(false);     eventbase->GetElectronSel()->SetSCEta(2.5);
           eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
           eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
           eventbase->GetElectronSel()->SetdxySigMax(4.);
           eventbase->GetElectronSel()->SetApplyConvVeto(true);
         std::vector<snu::KElectron> EleTElEnDownColl; eventbase->GetElectronSel()->Selection(EleTElEnDownColl, "SystDownElEn");
  
  
           eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
           eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
         std::vector<snu::KJet> jetJESUpColl; eventbase->GetJetSel()->Selection(jetJESUpColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystUpJES");
           eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
           eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
         std::vector<snu::KJet> jetJERUpColl; eventbase->GetJetSel()->Selection(jetJERUpColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystUpJER");
           eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
           eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
         std::vector<snu::KJet> jetJESDownColl; eventbase->GetJetSel()->Selection(jetJESDownColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystDownJES");
           eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
           eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
         std::vector<snu::KJet> jetJERDownColl; eventbase->GetJetSel()->Selection(jetJERDownColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystDownJER");
  
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
           reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
           reco_weight_ele_ElEnup = mcdata_correction->ElectronRecoScaleFactor(EleTElEnUpColl);
           reco_weight_ele_ElEndown = mcdata_correction->ElectronRecoScaleFactor(EleTElEnDownColl);
           if(UnPreTrig) nvtx_reweight = mcdata_correction->UserPileupWeight(eventbase->GetEvent(), jetColl.size()); //Measured in DY UnPreTrig
           if(PreTrig) nvtx_reweight = GetPreTrigPURW(Nvtx);
           trignorm_sf=0.947182;
         }
         float systweight_central  =weight*trignorm_sf*trigger_sf*nvtx_reweight*id_weight_ele       *reco_weight_ele       ;

         float systweight_Nvtxup   =weight*trignorm_sf*trigger_sf              *id_weight_ele       *reco_weight_ele       ;
         float systweight_ElEnup   =weight*trignorm_sf*trigger_sf*nvtx_reweight*id_weight_ele_ElEnup*reco_weight_ele_ElEnup;
         float systweight_PUup     =weight*trignorm_sf*trigger_sf*nvtx_reweight*id_weight_ele       *reco_weight_ele       *pileup_reweight_systup;
  
         float systweight_Nvtxdown =weight*trignorm_sf*trigger_sf              *id_weight_ele         *reco_weight_ele         ;
         float systweight_ElEndown =weight*trignorm_sf*trigger_sf*nvtx_reweight*id_weight_ele_ElEndown*reco_weight_ele_ElEndown;
         float systweight_PUdown   =weight*trignorm_sf*trigger_sf*nvtx_reweight*id_weight_ele         *reco_weight_ele         *pileup_reweight_systdown;


         if(PreTrig){
           //Central
           DoSystRun("NormCheck", "PreTrig",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     systweight_central);

           //Variation Up
           DoSystRun("NormCheck", "PreTrigSystUpPU",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     systweight_PUup);
           DoSystRun("NormCheck", "PreTrigSystUpNvtx",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     systweight_Nvtxup);
           DoSystRun("NormCheck", "PreTrigSystUpElEn",
                     EleTElEnUpColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi),
                     systweight_ElEnup);
           DoSystRun("NormCheck", "PreTrigSystUpMuEn",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi),
                     systweight_central);//No Muon is used, so fine to put anything there.
           DoSystRun("NormCheck", "PreTrigSystUpUncl",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi),
                     systweight_central);
           DoSystRun("NormCheck", "PreTrigSystUpJES",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESUpColl, bjetColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi),
                     systweight_central);
           DoSystRun("NormCheck", "PreTrigSystUpJER",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERUpColl, bjetColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi),
                     systweight_central);
           //Variation Down
           DoSystRun("NormCheck", "PreTrigSystDownPU",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     systweight_PUdown);
           DoSystRun("NormCheck", "PreTrigSystDownNvtx",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     systweight_Nvtxdown);
           DoSystRun("NormCheck", "PreTrigSystDownElEn",
                     EleTElEnDownColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi),
                     systweight_ElEndown);
           DoSystRun("NormCheck", "PreTrigSystDownMuEn",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi),
                     systweight_central);//No Muon is used, so fine to put anything there.
           DoSystRun("NormCheck", "PreTrigSystDownUncl",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi),
                     systweight_central);
           DoSystRun("NormCheck", "PreTrigSystDownJES",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESDownColl, bjetColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi),
                     systweight_central);
           DoSystRun("NormCheck", "PreTrigSystDownJER",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERDownColl, bjetColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi),
                     systweight_central);
         }//End of PreTrig
         if(UnPreTrig){
           //Central
           DoSystRun("NormCheck", "UnPreTrig",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     systweight_central);

           //Variation Up
           DoSystRun("NormCheck", "UnPreTrigSystUpPU",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     systweight_PUup);
           DoSystRun("NormCheck", "UnPreTrigSystUpElEn",
                     EleTElEnUpColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi),
                     systweight_ElEnup);
           DoSystRun("NormCheck", "UnPreTrigSystUpMuEn",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi),
                     systweight_central);//No Muon is used, so fine to put anything there.
           DoSystRun("NormCheck", "UnPreTrigSystUpUncl",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi),
                     systweight_central);
           DoSystRun("NormCheck", "UnPreTrigSystUpJES",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESUpColl, bjetColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi),
                     systweight_central);
           DoSystRun("NormCheck", "UnPreTrigSystUpJER",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERUpColl, bjetColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi),
                     systweight_central);
           //Variation Down
           DoSystRun("NormCheck", "UnPreTrigSystDownPU",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                     systweight_PUdown);
           DoSystRun("NormCheck", "UnPreTrigSystDownElEn",
                     EleTElEnDownColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi),
                     systweight_ElEndown);
           DoSystRun("NormCheck", "UnPreTrigSystDownMuEn",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi),
                     systweight_central);//No Muon is used, so fine to put anything there.
           DoSystRun("NormCheck", "UnPreTrigSystDownUncl",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi),
                     systweight_central);
           DoSystRun("NormCheck", "UnPreTrigSystDownJES",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESDownColl, bjetColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi),
                     systweight_central);
           DoSystRun("NormCheck", "UnPreTrigSystDownJER",
                     electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERDownColl, bjetColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi),
                     systweight_central);
         }//End of UnPreTrig
       }//End of !isData
     }//End of SystRun
     
   }//End of NormCheck
   if(OptMETMTWCuts){

     OptimiseMETMTWCuts(electronColl, electronLooseColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight, "SiglPreTrig");

   }
   if(IsoIPOpt){
     std::vector<snu::KElectron> electronPOGMVAMColl;
       for(int i=0; i<(int) electronPreColl.size(); i++){
         if(electronPreColl.at(i).Pt()<25) continue;
         if(PassIDCriteria(electronPreColl.at(i), "POGMVAM")) electronPOGMVAMColl.push_back(electronPreColl.at(i));
       }

     OptimiseIsoIPWP(electronPOGMVAMColl, weight, "");

   }//End of IsoIPOpt
   if(FRMeasure){

     const int NPtEdges=6;
     float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100.};
     if(!isData) trignorm_sf=0.947182;
     weight*=trignorm_sf;

     if(SiglWP){

       MeasureFakeRate(electronPreColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight, NPtEdges, PtEdges,
            "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", "_WP90Isop06IPp025p1sig4_L928878Isop4", "SiglPreTrig");

     }
     if(FRScan){

       //Full Scan
       const int NIsoCuts=7, NMVACuts=21;
       float IsoCuts[NIsoCuts]={0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
       float MVACuts[NMVACuts]={-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
  
       ScanFakeRate(electronFakeLPreOptColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met_x, met_y, weight, NMVACuts, MVACuts, NIsoCuts, IsoCuts, NPtEdges, PtEdges, "SiglPreTrig2D1D");

     }
   }
   if(SSDilepCR){

     float fake_weight1=1., fake_weight2=1.;
     if(k_running_nonprompt){
       fake_weight1 = GetFakeWeight(muonLooseColl , electronLooseColl, "POGTIsop6IPp2p1sig4", "POGTIsop20IPp01p05sig4Chi4", "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", "NoFilterConeSUSY");
       fake_weight2 = GetFakeWeight(muonLooseColl , electronLooseColl, "POGTIsop6IPp2p1sig4", "POGTIsop20IPp01p05sig4Chi4", "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", "TrkIsoVVLConeSUSY");

     }
     if(!isData && (k_sample_name.Contains("TT_powheg") || k_sample_name.Contains("DYJets") || k_sample_name.Contains("WJets"))){
       int NLepCand=0;
       for(int i=0; i<(int) muonColl.size(); i++){
         int LepType=GetLeptonType(muonColl.at(i), truthColl);
         if(LepType>0 || LepType<-4) NLepCand++;
       }
       for(int i=0; i<(int) electronColl.size(); i++){
         int LepType=GetLeptonType(electronColl.at(i), truthColl);
         if(LepType>0 || LepType<-4) NLepCand++;
       }
       if(NLepCand!=2) return;
     }
   
     CheckSSDilepCRs(muonColl, muonLooseColl, electronColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight*fake_weight1, "_POGTIsop6IPp2p1sig4_NoFilter", "");
     CheckSSDilepCRs(muonColl, muonLooseColl, electronColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight*fake_weight2, "_POGTIsop6IPp2p1sig4_IsoFilter", "");
     
   }
   if(Closure){

     if(!SystRun){

       if(k_sample_name.Contains("TT_powheg")){ if(NLeptonicBosonDecay(truthColl)==2) return; }

       std::vector<snu::KMuon> muonPOGLIsop4IPp5p1Chi100Coll,
                               muonPOGTIsop4IPp2p1NoChiColl     , muonPOGTIsop6IPp2p1Coll,
                               muonPOGTIsop4IPp2p1sig4NoChiColl , muonPOGTIsop6IPp2p1sig4NoChiColl,
                                                                  muonPOGTIsop6IPp2p1sig4Coll,
                               muonPOGTIsop4IPp01p05sig4Chi4Coll, muonPOGTIsop6IPp01p05sig4Chi4Coll;


       for(int i=0; i<(int) muonPreColl.size(); i++){
         //Round2
         if(PassIDCriteria(muonPreColl.at(i),"POGLIsop4IPp5p1Chi100","Roch"))      muonPOGLIsop4IPp5p1Chi100Coll.push_back(muonPreColl.at(i));
         if(PassIDCriteria(muonPreColl.at(i),"POGTIsop6IPp2p1","Roch"))            muonPOGTIsop6IPp2p1Coll.push_back(muonPreColl.at(i));
         if(PassIDCriteria(muonPreColl.at(i),"POGTIsop6IPp2p1sig4","Roch"))        muonPOGTIsop6IPp2p1sig4Coll.push_back(muonPreColl.at(i));
         if(PassIDCriteria(muonPreColl.at(i),"POGTIsop6IPp01p05sig4Chi4","Roch"))  muonPOGTIsop6IPp01p05sig4Chi4Coll.push_back(muonPreColl.at(i));
       }

       float fake_weight0=1., fake_weight1=1., fake_weight2=1., fake_weight3=1., fake_weight4=1., fake_weight5=1.; 
       float fake_weight10=1., fake_weight11=1., fake_weight12=1., fake_weight13=1., fake_weight14=1., fake_weight15=1.; 
       if(k_running_nonprompt){
        fake_weight3 = GetFakeWeight(muonPOGTIsop6IPp2p1Coll          , electronLooseColl, "POGTIsop6IPp2p1"          , "POGTIsop20IPp01p05sig4Chi4", "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", "NoFilterConeSUSY");
        fake_weight4 = GetFakeWeight(muonPOGTIsop6IPp01p05sig4Chi4Coll, electronLooseColl, "POGTIsop6IPp01p05sig4Chi4", "POGTIsop20IPp01p05sig4Chi4", "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", "NoFilterConeSUSY");
        fake_weight5 = GetFakeWeight(muonPOGTIsop6IPp2p1sig4Coll      , electronLooseColl, "POGTIsop6IPp2p1sig4"      , "POGTIsop20IPp01p05sig4Chi4", "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", "NoFilterConeSUSY");

       }


       if(!k_running_nonprompt){
         CheckTrilepCRs(electronColl, electronLooseColl, muonColl, muonPOGTIsop6IPp2p1Coll , jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
           weight*fake_weight3, "_POGTIsop6IPp2p1", "");
         CheckTrilepCRs(electronColl, electronLooseColl, muonColl, muonPOGTIsop6IPp01p05sig4Chi4Coll, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
           weight*fake_weight4, "_POGTIsop6IPp01p05sig4Chi4", "");
         CheckTrilepCRs(electronColl, electronLooseColl, muonColl, muonPOGTIsop6IPp2p1sig4Coll      , jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
           weight*fake_weight5, "_POGTIsop6IPp2p1sig4", "");

       }
       else{
         CheckTrilepCRs(electronColl, electronLooseColl, muonPOGTIsop6IPp2p1Coll , muonPOGTIsop6IPp2p1Coll , jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
           weight*fake_weight3, "_POGTIsop6IPp2p1", "");
         CheckTrilepCRs(electronColl, electronLooseColl, muonPOGTIsop6IPp01p05sig4Chi4Coll, muonPOGTIsop6IPp01p05sig4Chi4Coll, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
           weight*fake_weight4, "_POGTIsop6IPp01p05sig4Chi4", "");
         CheckTrilepCRs(electronColl, electronLooseColl, muonPOGTIsop6IPp2p1sig4Coll      , muonPOGTIsop6IPp2p1sig4Coll      , jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
           weight*fake_weight5, "_POGTIsop6IPp2p1sig4", "");

       }
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
         eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
         eventbase->GetElectronSel()->SetBETrRegIncl(false);
         eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
         eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
         eventbase->GetElectronSel()->SetdxySigMax(4.);
         eventbase->GetElectronSel()->SetApplyConvVeto(true);
       std::vector<snu::KElectron> EleTElEnUpColl; eventbase->GetElectronSel()->Selection(EleTElEnUpColl, "SystUpElEn");

         eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
         eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
         eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
         eventbase->GetElectronSel()->SetBETrRegIncl(false);
         eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
         eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
         eventbase->GetElectronSel()->SetdxySigMax(4.);
         eventbase->GetElectronSel()->SetApplyConvVeto(true);
       std::vector<snu::KElectron> EleTElEnDownColl; eventbase->GetElectronSel()->Selection(EleTElEnDownColl, "SystDownElEn");

       LeptonVeto=true;
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJESUpColl; eventbase->GetJetSel()->Selection(jetJESUpColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystUpJES");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJERUpColl; eventbase->GetJetSel()->Selection(jetJERUpColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystUpJER");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJESDownColl; eventbase->GetJetSel()->Selection(jetJESDownColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystDownJES");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJERDownColl; eventbase->GetJetSel()->Selection(jetJERDownColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystDownJER");

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
       float trigger_sf=1., trigger_sf_up=1., trigger_sf_down=1.;
       if(!isData){
         //geneff_weight   = GenFilterEfficiency(k_sample_name);
         //gennorm_weight  = SignalNorm(k_sample_name, 200.);
         k_factor_weight = GetKFactor();
    
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl);
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl);
         trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
         btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

         id_weight_ele_ElEnup   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", EleTElEnUpColl);
         reco_weight_ele_ElEnup = mcdata_correction->ElectronRecoScaleFactor(EleTElEnUpColl);
         id_weight_mu_MuEnup    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuTMuEnUpColl);
         trk_weight_mu_MuEnup   = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnUpColl);

         id_weight_ele_ElEndown   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", EleTElEnDownColl);
         reco_weight_ele_ElEndown = mcdata_correction->ElectronRecoScaleFactor(EleTElEnDownColl);
         id_weight_mu_MuEndown    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuTMuEnDownColl);
         trk_weight_mu_MuEndown   = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnDownColl);

         if(ElectronMuon){
           float trigger_sf1          = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
           float trigger_sf2          = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
           float trigger_sf1_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"   ,"Leg1SystUp");
           float trigger_sf2_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","Leg1SystUp");
           float trigger_sf1_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"   ,"Leg2SystUp");
           float trigger_sf2_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","Leg2SystUp");
           float trigger_sf1_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"   ,"Leg1SystDown");
           float trigger_sf2_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","Leg1SystDown");
           float trigger_sf1_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"   ,"Leg2SystDown");
           float trigger_sf2_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","Leg2SystDown");
           float trigger_sf2_dzup     = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","DZSystUp");
           float trigger_sf2_dzdown   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","DZSystDown");
           float trigger_sf1_up   = trigger_sf1+sqrt(pow(trigger_sf1_leg1up  -trigger_sf1,2)+pow(trigger_sf1_leg2up  -trigger_sf1,2));
           float trigger_sf1_down = trigger_sf1-sqrt(pow(trigger_sf1_leg1down-trigger_sf1,2)+pow(trigger_sf1_leg2down-trigger_sf1,2));
           float trigger_sf2_up   = trigger_sf2+sqrt(pow(trigger_sf2_leg1up  -trigger_sf2,2)+pow(trigger_sf2_leg2up  -trigger_sf2,2)+pow(trigger_sf2_dzup  -trigger_sf2,2));
           float trigger_sf2_down = trigger_sf2-sqrt(pow(trigger_sf2_leg1down-trigger_sf2,2)+pow(trigger_sf2_leg2down-trigger_sf2,2)+pow(trigger_sf2_dzdown-trigger_sf2,2));

           trigger_sf      = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
           trigger_sf_up   = ((Pass_TriggerBG ? trigger_sf1_up:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_up:0.)*LumiH)/LumiBH;
           trigger_sf_down = ((Pass_TriggerBG ? trigger_sf1_down:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_down:0.)*LumiH)/LumiBH;
         }

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
         fake_weight = GetFakeWeight(muonLooseColl, electronLooseColl, "POGTIsop6IPp2p1sig4" , "POGTIsop20IPp01p05sig4Chi4", "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", bjetNoVetoColl, "TrkIsoVVLConeSUSY");
       }
       float systweight_central=weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight *trigger_sf;

       float systweight_Trigup =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight *trigger_sf_up;
       float systweight_ElEnup =weight*k_factor_weight*id_weight_ele_ElEnup*reco_weight_ele_ElEnup*id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight *trigger_sf;
       float systweight_MuEnup =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu_MuEnup*trk_weight_mu_MuEnup*btag_sf        *fake_weight *trigger_sf;
       float systweight_JESup  =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JESup  *fake_weight *trigger_sf;
       float systweight_JERup  =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JERup  *fake_weight *trigger_sf;
       float systweight_LTagup =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_LTagup *fake_weight *trigger_sf;
       float systweight_BCTagup=weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_BCTagup*fake_weight *trigger_sf;
       float systweight_PUup   =weight*k_factor_weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight*pileup_reweight_systup*trigger_sf;

       float systweight_Trigdown =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *fake_weight *trigger_sf_down;
       float systweight_ElEndown =weight*k_factor_weight*id_weight_ele_ElEndown*reco_weight_ele_ElEndown*id_weight_mu         *trk_weight_mu         *btag_sf          *fake_weight *trigger_sf;
       float systweight_MuEndown =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu_MuEndown*trk_weight_mu_MuEndown*btag_sf          *fake_weight *trigger_sf;
       float systweight_JESdown  =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JESdown  *fake_weight *trigger_sf;
       float systweight_JERdown  =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JERdown  *fake_weight *trigger_sf;
       float systweight_LTagdown =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_LTagdown *fake_weight *trigger_sf;
       float systweight_BCTagdown=weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_BCTagdown*fake_weight *trigger_sf;
       float systweight_PUdown   =weight*k_factor_weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *fake_weight*pileup_reweight_systdown*trigger_sf;


       //----------------------------------------------------------------------------------------------------------------------//


       DoSystRun("Closure", "",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_central);
       if(!isData){
         DoSystRun("Closure", "SystUpTrig",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_Trigup);
         DoSystRun("Closure", "SystUpPU",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_PUup);
         DoSystRun("Closure", "SystUpUncl",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi),
                   systweight_central);
         DoSystRun("Closure", "SystUpJES",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESUpColl, bjetJESUpColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi),
                   systweight_JESup);
         DoSystRun("Closure", "SystUpJER",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERUpColl, bjetJERUpColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi),
                   systweight_JERup);
         DoSystRun("Closure", "SystUpBTag_L",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_LTagup);
         DoSystRun("Closure", "SystUpBTag_BC",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_BCTagup);
         DoSystRun("Closure", "SystUpElEn", EleTElEnUpColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi),
                   systweight_ElEnup);
         DoSystRun("Closure", "SystUpMuEn", electronColl, electronLooseColl, MuTMuEnUpColl, muonLooseColl, jetColl, bjetColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi),
                   systweight_MuEnup);


         DoSystRun("Closure", "SystDownTrig",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_Trigdown);
         DoSystRun("Closure", "SystDownPU",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_PUdown);
         DoSystRun("Closure", "SystDownUncl",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi),
                   systweight_central);
         DoSystRun("Closure", "SystDownJES",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESDownColl, bjetJESDownColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi),
                   systweight_JESdown);
         DoSystRun("Closure", "SystDownJER",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERDownColl, bjetJERDownColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi),
                   systweight_JERdown);
         DoSystRun("Closure", "SystDownBTag_L",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_LTagdown);
         DoSystRun("Closure", "SystDownBTag_BC",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_BCTagdown);
         DoSystRun("Closure", "SystDownElEn", EleTElEnDownColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi),
                   systweight_ElEndown);
         DoSystRun("Closure", "SystDownMuEn", electronColl, electronLooseColl, MuTMuEnDownColl, muonLooseColl, jetColl, bjetColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi),
                   systweight_MuEndown);

       }
       else if( isData && k_running_nonprompt ){

//         DoSystRun("Closure", "", electronFakeLColl, electronFakeLColl, muonLooseColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), systweight);

       }
     }//End of SystRun
   }//End of Closure

   if(HighdXYFRCheck){

     //Additional Selection for highdxy measurement-----------------------------------------------------------------//
     std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);

       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
       eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
       eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetApplyConvVeto(true);
     std::vector<snu::KElectron> EleHighdxyTColl; eventbase->GetElectronSel()->Selection(EleHighdxyTColl);

       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_HctoWA_FAKELOOSE);
       eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
       eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.4, 0.4);
       eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetApplyConvVeto(true);
     std::vector<snu::KElectron> EleHighdxyLColl; eventbase->GetElectronSel()->Selection(EleHighdxyLColl);
     //-------------------------------------------------------------------------------------------------------------//


     MeasureHighdXYFakeRate(electronColl, electronLooseColl, EleHighdxyTColl, EleHighdxyLColl, muonColl, muonLooseColl, jetNoVetoColl, bjetNoVetoColl, truthColl, met, met*cos(metphi), met*sin(metphi), weight, "", "");
   }//End of HighdXYCheck

   if(IDValidation){
     
     ValidateID(electronColl, electronPreColl, weight);
     
   }//End of validation

/////////////////////////////////////////////////////////////////////////////////// 

return;

}// End of execute event loop
  


void Aug2017_DataEleFRCalc::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Aug2017_DataEleFRCalc::BeginCycle() throw( LQError ){
  
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

Aug2017_DataEleFRCalc::~Aug2017_DataEleFRCalc() {
  
  Message("In Aug2017_DataEleFRCalc Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Aug2017_DataEleFRCalc::CheckNormCR(std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


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
        if(MTW>70.){
          FillHist("Count_UnPreNormCR"+Label, 0., weight, 0., 4., 4);
          FillHist("PTe_met50mtw70"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
          FillHist("Etae_met50mtw70"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
          FillHist("Nj_met50mtw70"+Label, JetVetoColl.size(), weight, 0., 10., 10);
          FillHist("MET_met50mtw70"+Label, MET, weight, 0., 200., 40);
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
          if(MTW>70.){
            FillHist("Count_UnPreNormCR"+Label, 1., weight, 0., 4., 4); 
            FillHist("Count_PreNormCR"+Label,   0., weight, 0., 2., 2); 

            FillHist("PTe_gt1j40met50mtw70"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
            FillHist("Etae_gt1j40met50mtw70"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
            FillHist("PTj1_gt1j40met50mtw70"+Label, JetVetoColl.at(0).Pt(), weight, 0., 200., 40);
            FillHist("Etaj1_gt1j40met50mtw70"+Label, JetVetoColl.at(0).Eta(), weight, -5., 5., 20);
            FillHist("dRej1_gt1j40met50mtw70"+Label, EleTColl.at(0).DeltaR(JetVetoColl.at(0)), weight, 0., 5., 50);
            FillHist("Nj_gt1j40met50mtw70"+Label, JetVetoColl.size(), weight, 0., 10., 10);
            FillHist("MET_gt1j40met50mtw70"+Label, MET, weight, 0., 200., 40);

            FillHist("Nvtx_gt1j40met50mtw70"+Label, Nvtx, weight, 0., 50., 50);
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


void Aug2017_DataEleFRCalc::OptimiseMETMTWCuts(std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, float MET, float METx, float METy, float weight, TString Option){

  bool SiglPreTrig=false, MultPreTrig=false;
  float MinElePt=0.;
  if     (Option.Contains("SiglPreTrig")) {SiglPreTrig = true; MinElePt=25.;}
  else if(Option.Contains("MultPreTrig")) {MultPreTrig = true; MinElePt=15.;}

  std::vector<snu::KJet> JetVetoColl = SkimJetColl(JetColl, EleLColl, MuLColl, "EleMuVeto");


  if(EleLColl.size()==1 && MuLColl.size()==0 && JetVetoColl.size()>0){
    bool PassJetReq=false;
    float PTCorr=ConeCorrectedPT(EleLColl.at(0), 0.06), fEta=fabs(EleLColl.at(0).Eta());
    //float PTCorr=EleLColl.at(0).Pt()*(1.+EleLColl.at(0).PFRelIso(0.3)), fEta=fabs(EleLColl.at(0).Eta());

    //Selection Requirement--------------------------------------------------------------//
    if( EleLColl.at(0).Pt()<MinElePt ) return;
    if( PTCorr<25. ) return;
    for(int j=0; j<(int) JetVetoColl.size(); j++){
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



void Aug2017_DataEleFRCalc::OptimiseIsoIPWP(std::vector<snu::KElectron> EleTColl, float weight, TString Option){
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


  for(int i=0; i<(int) EleTColl.size(); i++){
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

void Aug2017_DataEleFRCalc::MeasureFakeRate(std::vector<snu::KElectron> ElePreColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, int NPtEdges, float PtEdges[], TString LooseID, TString TightID, TString Label, TString Option){
  //BJet : NoVetoColl,  Jet: No matter whether lepton vetoed or not, it will veto in code any way.

  std::vector<snu::KElectron> EleLColl;
    for(int i=0; i<(int) ElePreColl.size(); i++){ if(PassIDCriteria(ElePreColl.at(i),LooseID)) EleLColl.push_back(ElePreColl.at(i)); }

  bool SiglPreTrig=false, MultPreTrig=false;
  float MinElePt=0.;
  if     (Option.Contains("SiglPreTrig")) {SiglPreTrig = true; MinElePt=25.;}
  else if(Option.Contains("MultPreTrig")) {MultPreTrig = true; MinElePt=15.;}

  std::vector<snu::KJet> JetVetoColl = SkimJetColl(JetColl, EleLColl, MuLColl, "EleMuVeto");
  std::vector<snu::KJet> BJetVetoColl = SelBJets(JetVetoColl, "Medium");
  int NBJets = BJetVetoColl.size();


  if(EleLColl.size()==1 && MuLColl.size()==0 && JetVetoColl.size()>0){
    bool  PassJetReq=false;
    float PTCorr=ConeCorrectedPT(EleLColl.at(0), 0.06), fEta=fabs(EleLColl.at(0).Eta());

    //Selection Requirement--------------------------------------------------------------//
    if( EleLColl.at(0).Pt()<MinElePt ) return;
    if( PTCorr<25. ) return;
    for(int j=0; j<(int) JetVetoColl.size(); j++){
      if(JetVetoColl.at(j).Pt()<40) continue;
      if(JetVetoColl.at(j).DeltaR(EleLColl.at(0))<1.0) continue;
      PassJetReq=true; break;
    }
    if(!PassJetReq) return;
    float MTW = sqrt(2)*sqrt(MET*EleLColl.at(0).Pt()-METx*EleLColl.at(0).Px()-METy*EleLColl.at(0).Py());

    if( !(MET<25 && MTW<25) ) return;
    //if( MET<30 && MTW<25 ) return;//Final Optimised
    //-----------------------------------------------------------------------------------//

    bool IsNearB = IsNearBJet(EleLColl.at(0), BJetColl);
    if(fEta<0.8){
      FillHist("EleB1SumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
   
      if(NBJets>0){ FillHist("EleB1SumW_HasB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
                   if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_HasB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1); } 
      if(IsNearB){ FillHist("EleB1SumW_NearB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
                   if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_NearB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1); } 
      else{ FillHist("EleB1SumW_AwayB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_AwayB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1); } 
    }
    else if(fEta<1.479){
      FillHist("EleB2SumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
   
      if(NBJets>0){ FillHist("EleB2SumW_HasB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
                   if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_HasB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1); } 
      if(IsNearB){ FillHist("EleB2SumW_NearB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
                   if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_NearB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1); } 
      else { FillHist("EleB2SumW_AwayB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
             if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_AwayB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1); } 
    }
    else{
      FillHist("EleESumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
   
      if(NBJets>0){ FillHist("EleESumW_HasB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
                   if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_HasB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1); } 
      if(IsNearB){ FillHist("EleESumW_NearB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
                   if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_NearB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1); } 
      else       { FillHist("EleESumW_AwayB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
                   if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_AwayB_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1); } 
    }
  }//End of 1e+0mu+geq1j
 
}


void Aug2017_DataEleFRCalc::ScanFakeRate(std::vector<snu::KElectron> EleColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, int NMVACuts, float MVACuts[], int NIsoCuts, float IsoCuts[], int NPtEdges, float PtEdges[], TString Option){
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
      for(int i_ele=0; i_ele<(int) EleColl.size(); i_ele++){
        if(EleColl.at(i_ele).PFRelIso(0.3)<IsoCuts[it_iso] && EleColl.at(i_ele).MVA()>MVACuts[it_mva]) NLooseEle++;
      }
      if(NLooseEle!=1) continue;
      if(MuLColl.size()!=0) continue;
      if(EleColl.at(0).Pt()<MinElePt) continue;
      bool PassJetReq=false; 
      for(int j=0; j<(int) JetColl.size(); j++){
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


void Aug2017_DataEleFRCalc::MeasureHighdXYFakeRate(std::vector<snu::KElectron> EleLowd0TColl, std::vector<snu::KElectron> EleLowd0LColl, std::vector<snu::KElectron> EleHighd0TColl, std::vector<snu::KElectron> EleHighd0LColl, std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, std::vector<snu::KTruth> TruthColl, float MET, float METx, float METy, float weight, TString Label, TString Option){
     //jet:LeptonNoVetoedColl both for ambuguity and near jet checking

  //Series of Cuts for FR measurement--------------------------------------------------------------------------------//
  bool PassSelHighd0=true, PassSelLowd0=true; float MTW=0.;
  if( EleHighd0LColl.size()!=1 ) PassSelHighd0=false;
  if( PassSelHighd0 && !(EleHighd0LColl.at(0).Pt()>25. && fabs(EleHighd0LColl.at(0).dxySig())>5.) ) PassSelHighd0=false;
  if( PassSelHighd0 && !(JetColl.size()>0 && JetColl.at(0).Pt()>40) ) PassSelHighd0=false;
  if( PassSelHighd0 ){
    bool HasAway40Jet=false;
    for( int i=0; i<(int) JetColl.size(); i++){
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
    for( int i=0; i<(int) JetColl.size(); i++){
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



void Aug2017_DataEleFRCalc::CheckSSDilepCRs(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetNoVetoColl, std::vector<snu::KJet> BJetNoVetoColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool MuMu=Option.Contains("MuMu"), EMu=!MuMu;//By default trimu.
  std::vector<snu::KJet> JetColl  = SkimJetColl(JetNoVetoColl,  EleLColl, MuLColl, "EleMuVeto");
  std::vector<snu::KJet> BJetColl = SelBJets(JetColl, "Medium");

  if(MuMu){

  if( !(MuLColl.size()==2 && EleLColl.size()==0) ) return;
  if( !(MuTColl.size()==2) ) return;
  if(  BJetColl.size()==0  ) return;
  if( !(JetColl.size()>=3) ) return;
  
  if( !(MuTColl.at(0).Charge()==MuTColl.at(1).Charge()) ) return;
  float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();

  if( fabs(Mmumu-91.2)<10 || Mmumu<12) return;

  FillHist("PTMu1" +Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
  FillHist("PTMu2" +Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
  FillHist("EtaMu1"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
  FillHist("EtaMu2"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
  FillHist("Mmumu" +Label, Mmumu, weight, 0., 200., 40);
  FillHist("Nj"    +Label, JetColl.size(), weight, 0., 10., 10);
  FillHist("Nb"    +Label, BJetColl.size(), weight, 0., 10., 10);
  FillHist("MET"   +Label, MET, weight, 0., 200., 40);

  }//End of SS MuMu
  else if(EMu){

  if( !(MuLColl.size()==1 && EleLColl.size()==1) ) return;
  if( !(MuTColl.size()==1 && EleTColl.size()==1) ) return;
  if(  BJetColl.size()==0  ) return;
  if( !(JetColl.size()>=3) ) return;
  
  if( !(MuTColl.at(0).Charge()==EleTColl.at(0).Charge()) ) return;

  FillHist("PTMu1" +Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
  FillHist("PTEl1" +Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
  FillHist("EtaMu1"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
  FillHist("EtaEl1"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
  FillHist("Nj"    +Label, JetColl.size(), weight, 0., 10., 10);
  FillHist("Nb"    +Label, BJetColl.size(), weight, 0., 10., 10);
  FillHist("MET"   +Label, MET, weight, 0., 200., 40);

  }//End of SS EMu

}


void Aug2017_DataEleFRCalc::CheckTrilepCRs(std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetNoVetoColl, std::vector<snu::KJet> BJetNoVetoColl, float MET, float METx, float METy, float weight, TString Label, TString Option){

  std::vector<snu::KJet> JetColl  = SkimJetColl(JetNoVetoColl, EleLColl, MuLColl, "EleMuVeto");
  std::vector<snu::KJet> BJetColl = SelBJets(JetColl, "Medium");

  if( !(EleLColl.size()==1 && MuLColl.size()==2) ) return;
  if( !(EleTColl.size()==1 && MuTColl.size()==2) ) return;
  if( !(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) ) return;
  if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && MuTColl.at(1).Pt()>10) ) return;
  if( !(fabs(EleTColl.at(0).Eta())<2.5) ) return;

  float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
  float MTW = sqrt(2)*sqrt(MET*EleTColl.at(0).Pt()-METx*EleTColl.at(0).Px()-METy*EleTColl.at(0).Py());
  float M3l=(EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M();
  if(Mmumu<12) return;
 

  bool HasBJet=false, OnZ=false, OnZG=false, WZSel=false, ZGSel=false, ttZSel=false;
  if( BJetColl.size()!=0 )                      HasBJet = true;
  if( fabs(Mmumu-91.2)<10     )                 OnZ     = true;
  if( fabs(M3l-91.2)<10       )                 OnZG    = true;
  if( OnZ && M3l>101.2 && MET>50 )              WZSel   = true;
  if( OnZG && Mmumu<81.2 && MET<50 )            ZGSel   = true;
  if( OnZ && HasBJet && JetColl.size()>2)       ttZSel  = true;


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
  if(OnZ){
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
  if(!HasBJet & OnZ){
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
  if(!OnZ && BJetColl.size()==1 && JetColl.size()==1){
    FillHist("PTe_TopFakeCR"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("Etae_TopFakeCR"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("PTmu1_TopFakeCR"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
    FillHist("PTmu2_TopFakeCR"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("Etamu1_TopFakeCR"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
    FillHist("Etamu2_TopFakeCR"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
  
    FillHist("Mmumu_TopFakeCR"+Label, Mmumu, weight, 0., 200., 40);
    FillHist("M3l_TopFakeCR"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 100);
    FillHist("Nj_TopFakeCR"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_TopFakeCR"+Label, BJetColl.size(), weight, 0., 10., 10);
    FillHist("MET_TopFakeCR"+Label, MET, weight, 0., 200., 40);
    FillHist("MTW_TopFakeCR"+Label, MTW, weight, 0., 200., 40);
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


void Aug2017_DataEleFRCalc::ValidateID(std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, float weight){
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
      for(int i=0; i<(int) EleLColl.size(); i++){
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


bool Aug2017_DataEleFRCalc::IsNearBJet(snu::KElectron Ele, std::vector<snu::KJet> bjetNoVetoColl){

  bool IsNearB=false;
  for(int i=0; i<(int) bjetNoVetoColl.size(); i++){
    if(Ele.DeltaR(bjetNoVetoColl.at(i))<0.4){
      IsNearB=true; break;
    }
  }

  return IsNearB;
}


float Aug2017_DataEleFRCalc::ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut){

  float PTCorr=Ele.Pt()*(1.+max(0.,Ele.PFRelIso(0.3)-TightIsoCut));
  return PTCorr;

}

float Aug2017_DataEleFRCalc::ConeCorrectedPT(snu::KMuon Mu, float TightIsoCut){

  float PTCorr=Mu.Pt()*(1.+max((float) 0.,RochIso(Mu,"0.4")-TightIsoCut));
  return PTCorr;

}



void Aug2017_DataEleFRCalc::DoSystRun(TString Cycle, TString Mode, std::vector<snu::KElectron> EleColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight){

  int SystDir=0; TString SystKindLabel="", SystDirLabel="";
  if(Mode.Contains("Syst")){
    if     (isData)                   return;
    if     (Mode.Contains("Up"))      {SystDir= 1; SystDirLabel="_systup";}
    else if(Mode.Contains("Down"))    {SystDir=-1; SystDirLabel="_systdown";}
    if     (Mode.Contains("PU"))      SystKindLabel="_PU";
    else if(Mode.Contains("Nvtx"))    SystKindLabel="_Nvtx";
    else if(Mode.Contains("JES"))     SystKindLabel="_JES";
    else if(Mode.Contains("JER"))     SystKindLabel="_JER";
    else if(Mode.Contains("Uncl"))    SystKindLabel="_Uncl";
    else if(Mode.Contains("ElEn"))    SystKindLabel="_ElEn";
    else if(Mode.Contains("MuEn"))    SystKindLabel="_MuEn";
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
  else if(Cycle=="Closure"){
    CheckTrilepCRs(EleColl, EleLColl, MuColl, MuLColl, JetColl, BJetColl, MET, METx, METy, weight, SystDirLabel+SystKindLabel);
  }//End of Closure Cycle

  return;
}


float Aug2017_DataEleFRCalc::GetFakeWeight(std::vector<snu::KMuon>& MuLColl, std::vector<snu::KElectron>& EleLColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, TString Option){

  float fakeweight=-1.; int NLooseNotTight=0;

  TString FilterInfo="", ConeMethod="";
  if     (Option.Contains("NoFilter"))  FilterInfo="NoFilter";
  else if(Option.Contains("TrkIsoVVL")) FilterInfo="TrkIsoVVL";
  if     (Option.Contains("ConeE"))     ConeMethod="ConeE";
  else if(Option.Contains("ConeSUSY"))  ConeMethod="ConeSUSY";
  else if(Option.Contains("FOPt"))      ConeMethod="FOPt";


  for(int i=0; i<(int) MuLColl.size(); i++){
    if(!PassIDCriteria(MuLColl.at(i), MuTID)){
      float FR=0.;
      FR=FakeRateData(MuLColl.at(i),MuTID+"_"+MuLID+"_"+FilterInfo+ConeMethod);
      //cout<<"MuFR"<<FR<<endl;
      fakeweight*=-FR/(1-FR);
      NLooseNotTight++;
    }
  }
  for(int i=0; i<(int) EleLColl.size(); i++){
    if(!PassIDCriteria(EleLColl.at(i), EleTID)){
      float FR=FakeRateData(EleLColl.at(i), "QCD_"+EleTID+"_"+EleLID+"_"+ConeMethod);
      fakeweight*=-FR/(1-FR);
      NLooseNotTight++;
        //cout<<"ElFR"<<FR<<endl;
    }
  }
  if(NLooseNotTight==0) return 0.;

  return fakeweight;
}


float Aug2017_DataEleFRCalc::GetFakeWeight(std::vector<snu::KMuon>& MuLColl, TString MuLID, TString MuTID, vector<snu::KJet>& BJetNoVetoColl, TString Option){

  TString FilterInfo="", ConeMethod="";
  if     (Option.Contains("NoFilter"))  FilterInfo="NoFilter";
  else if(Option.Contains("TrkIsoVVL")) FilterInfo="TrkIsoVVL";
  if     (Option.Contains("ConeE"))     ConeMethod="ConeE";
  else if(Option.Contains("ConeSUSY"))  ConeMethod="ConeSUSY";
  else if(Option.Contains("FOPt"))      ConeMethod="FOPt";

  int NLooseNotTight=0; float fakeweight=-1.;
  for(int i=0; i<(int) MuLColl.size(); i++){
    if(!PassIDCriteria(MuLColl.at(i), MuTID,"Roch")){
      float FR=0.;
      bool IsNearB = IsNearBJet(MuLColl.at(i), BJetNoVetoColl);
      TString FlavOpt = IsNearB? "_NearB":"_AwayB";
      
      FR=FakeRateData(MuLColl.at(i),MuTID+"_"+MuLID+FlavOpt+"_"+FilterInfo+ConeMethod);
      fakeweight*=-FR/(1-FR);
      NLooseNotTight++;
    }
  }
  if(NLooseNotTight==0) return 0.;

  return fakeweight;
}



float Aug2017_DataEleFRCalc::GetFakeWeight(std::vector<snu::KMuon>& MuLColl, std::vector<snu::KElectron>& EleLColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, vector<snu::KJet>& BJetNoVetoColl, TString Option){

  float fakeweight=-1.; int NLooseNotTight=0;

  TString FilterInfo="", ConeMethod="";
  if     (Option.Contains("NoFilter"))  FilterInfo="NoFilter";
  else if(Option.Contains("TrkIsoVVL")) FilterInfo="TrkIsoVVL";
  if     (Option.Contains("ConeE"))     ConeMethod="ConeE";
  else if(Option.Contains("ConeSUSY"))  ConeMethod="ConeSUSY";
  else if(Option.Contains("FOPt"))      ConeMethod="FOPt";


  for(int i=0; i<(int) MuLColl.size(); i++){
    if(!PassIDCriteria(MuLColl.at(i), MuTID, "Roch")){
      float FR=0.;
      bool IsNearB = IsNearBJet(MuLColl.at(i), BJetNoVetoColl);
      TString FlavOpt = IsNearB? "_NearB":"_AwayB";

      FR=FakeRateData(MuLColl.at(i),MuTID+"_"+MuLID+FlavOpt+"_"+FilterInfo+ConeMethod);
      //cout<<"MuFR"<<FR<<endl;
      fakeweight*=-FR/(1-FR);
      NLooseNotTight++;
    }
  }
  for(int i=0; i<(int) EleLColl.size(); i++){
    if(!PassIDCriteria(EleLColl.at(i), EleTID)){
      float FR=FakeRateData(EleLColl.at(i), "QCD_"+EleTID+"_"+EleLID+"_"+ConeMethod);
      fakeweight*=-FR/(1-FR);
      NLooseNotTight++;
        //cout<<"ElFR"<<FR<<endl;
    }
  }
  if(NLooseNotTight==0) return 0.;

  return fakeweight;
}


bool Aug2017_DataEleFRCalc::IsNearBJet(snu::KMuon Mu, std::vector<snu::KJet>& bjetNoVetoColl){

  bool IsNearB=false;
  for(int i=0; i<(int) bjetNoVetoColl.size(); i++){
    if(Mu.DeltaR(bjetNoVetoColl.at(i))<0.4){
      IsNearB=true; break;
    }
  }

  return IsNearB;
}


float Aug2017_DataEleFRCalc::FakeRateData(snu::KMuon Mu, TString Option){

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

  if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.342129;
      else if(PTCorr<20)  FR=0.161717;
      else if(PTCorr<25)  FR=0.124337;
      else if(PTCorr<35)  FR=0.124102;
      else                FR=0.112075;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.3853  ;
      else if(PTCorr<20)  FR=0.200014;
      else if(PTCorr<25)  FR=0.181963;
      else if(PTCorr<35)  FR=0.143144;
      else                FR=0.149275;
    }
    else{
      if     (PTCorr<15)  FR=0.386585;
      else if(PTCorr<20)  FR=0.228739;
      else if(PTCorr<25)  FR=0.203977;
      else if(PTCorr<35)  FR=0.190597;
      else                FR=0.180573;
    }

    if(Syst_FR && SystDir<0){
      if(fEta<0.9){
        if     (PTCorr<15)   FR*=(1.- 0.00570765 );
        else if(PTCorr<20)   FR*=(1.- 0.146032   );
        else if(PTCorr<25)   FR*=(1.- 0.079018   );
        else if(PTCorr<35)   FR*=(1.- 0.211878   );
        else                 FR*=(1.- 0.364666   );
      }
      else if(fEta<1.6){
        if     (PTCorr<15)   FR*=(1.- 0.0378566 );
        else if(PTCorr<20)   FR*=(1.- 0.173204  );
        else if(PTCorr<25)   FR*=(1.- 0.13663   );
        else if(PTCorr<35)   FR*=(1.- 0.0838593 );
        else                 FR*=(1.- 0.27241   );
      }
      else{
        if     (PTCorr<15)   FR*=(1.- 0.0761964 );
        else if(PTCorr<20)   FR*=(1.- 0.137059  );
        else if(PTCorr<25)   FR*=(1.- 0.246004  );
        else if(PTCorr<35)   FR*=(1.- 0.0928262 );
        else                 FR*=(1.- 0.214074  );
      }
    }
    else if(Syst_FR && SystDir>0){
      if(fEta<0.9){
        if     (PTCorr<15)   FR*=(1.+ 0.0467043  );
        else if(PTCorr<20)   FR*=(1.+ 0.00598351 );
        else if(PTCorr<25)   FR*=(1.+ 0.0564788  );
        else if(PTCorr<35)   FR*=(1.+ 0.0280901  );
        else                 FR*=(1.+ 0.0889982  );
      }
      else if(fEta<1.6){
        if     (PTCorr<15)   FR*=(1.+ 0.0339258 );
        else if(PTCorr<20)   FR*=(1.+ 0.0216653 );
        else if(PTCorr<25)   FR*=(1.+ 0.0344182 );
        else if(PTCorr<35)   FR*=(1.+ 0.072051  );
        else                 FR*=(1.+ 0.0674714 );
      }
      else{
        if     (PTCorr<15)   FR*=(1.+ 0.043051 );
        else if(PTCorr<20)   FR*=(1.+ 0.0239283);
        else if(PTCorr<25)   FR*=(1.+ 0.108454 );
        else if(PTCorr<35)   FR*=(1.+ 0.0785984);
        else                 FR*=(1.+ 0.0663562);
      }
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY"){
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
  else {cout<<"No Such FR!"<<endl;}

  return FR;
}



float Aug2017_DataEleFRCalc::FakeRateData(snu::KElectron Ele, TString Option){

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


  if(Option=="QCD_POGWP90Isop06IPp025p05sig4_LMVA06Isop4IPp025p05sig4_ConeSUSY"){
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
  }
  else if(Option=="QCD_POGWP90Isop06IPp025p1sig4_HctoWAFakeLoose_ConeSUSY"){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.10963  ;
      else if(PTCorr<50)  FR=0.072741 ;
      else                FR=0.0709206;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.110035 ;
      else if(PTCorr<50)  FR=0.0752341;
      else                FR=0.0718434;
    }
    else{
      if     (PTCorr<35)  FR=0.144356;
      else if(PTCorr<50)  FR=0.106974;
      else                FR=0.128996;
    }

    if(Syst_FR && SystDir<0){
      if(fEta<0.8){
        if     (PTCorr<35)   FR*=(1.- 0.0677277 );
        else if(PTCorr<50)   FR*=(1.- 0.404283  );
        else                 FR*=(1.- 0.442664  );
      }
      else if(fEta<1.479){
        if     (PTCorr<35)   FR*=(1.- 0.144574 );
        else if(PTCorr<50)   FR*=(1.- 0.43086  );
        else                 FR*=(1.- 0.224811 );
      }
      else{
        if     (PTCorr<35)   FR*=(1.- 0.055181  );
        else if(PTCorr<50)   FR*=(1.- 0.0939935 );
        else                 FR*=(1.- 0.135551  );
      }
    }
    else if(Syst_FR && SystDir>0){
      if(fEta<0.8){
        if     (PTCorr<35)   FR*=(1.+ 0.0956376 );
        else if(PTCorr<50)   FR*=(1.+ 0.169789  );
        else                 FR*=(1.+ 0.287828  );
      }
      else if(fEta<1.479){
        if     (PTCorr<35)   FR*=(1.+ 0.0953948 );
        else if(PTCorr<50)   FR*=(1.+ 0.127346  );
        else                 FR*=(1.+ 0.221431  );
      }
      else{
        if     (PTCorr<35)   FR*=(1.+ 0.136633  );
        else if(PTCorr<50)   FR*=(1.+ 0.0603159 );
        else                 FR*=(1.+ 0.047715  );
      }
    }
  }

  return FR;
} 


float Aug2017_DataEleFRCalc::GetPreTrigPURW(int Nvtx){

   float weight=1.;
   if(Nvtx<2) weight=2.42633;
   else if(Nvtx<3) weight=4.40376;
   else if(Nvtx<4) weight=4.1189;
   else if(Nvtx<5) weight=3.87844;
   else if(Nvtx<6) weight=3.84634;
   else if(Nvtx<7) weight=3.47177;
   else if(Nvtx<8) weight=3.16573;
   else if(Nvtx<9) weight=2.89114;
   else if(Nvtx<10) weight=2.57004;
   else if(Nvtx<11) weight=2.36693;
   else if(Nvtx<12) weight=2.13589;
   else if(Nvtx<13) weight=1.88736;
   else if(Nvtx<14) weight=1.67763;
   else if(Nvtx<15) weight=1.48384;
   else if(Nvtx<16) weight=1.29959;
   else if(Nvtx<17) weight=1.16943;
   else if(Nvtx<18) weight=1.00198;
   else if(Nvtx<19) weight=0.855362;
   else if(Nvtx<20) weight=0.759321;
   else if(Nvtx<21) weight=0.646064;
   else if(Nvtx<22) weight=0.581407;
   else if(Nvtx<23) weight=0.504966;
   else if(Nvtx<24) weight=0.4465;
   else if(Nvtx<25) weight=0.406343;
   else if(Nvtx<26) weight=0.377922;
   else if(Nvtx<27) weight=0.330518;
   else if(Nvtx<28) weight=0.311081;
   else if(Nvtx<29) weight=0.323852;
   else if(Nvtx<30) weight=0.275319;
   else if(Nvtx<31) weight=0.268598;
   else if(Nvtx<32) weight=0.222533;
   else if(Nvtx<33) weight=0.233084;
   else if(Nvtx<34) weight=0.223157;
   else if(Nvtx<35) weight=0.229524;
   else if(Nvtx<36) weight=0.324319;
   else if(Nvtx<37) weight=0.277683;
   else if(Nvtx<38) weight=0.200669;
   else if(Nvtx<39) weight=0.300899;
   else if(Nvtx<40) weight=0.289072;
   else if(Nvtx<41) weight=0.267605;
   else if(Nvtx<42) weight=0.226886;
   else if(Nvtx<43) weight=0.336191;
   else if(Nvtx<44) weight=0.57703;
   else if(Nvtx<45) weight=0.338796;
   else if(Nvtx<46) weight=0.451115;
   else if(Nvtx<47) weight=0.52386;
   else if(Nvtx<48) weight=0.475035;
   else if(Nvtx<49) weight=0.71763;
   else if(Nvtx<50) weight=1.65286;
   else             weight=0.910066;

   return weight;
}



void Aug2017_DataEleFRCalc::FillCutFlow(TString cut, float weight){
  
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



void Aug2017_DataEleFRCalc::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Aug2017_DataEleFRCalc::MakeHistograms(){
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
  // *  Remove//Overide this Aug2017_DataEleFRCalcCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Aug2017_DataEleFRCalc::ClearOutputVectors() throw(LQError) {

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
