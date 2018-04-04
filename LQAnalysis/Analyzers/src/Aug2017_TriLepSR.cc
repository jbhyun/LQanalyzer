// $Id: Aug2017_TriLepSR.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQAug2017_TriLepSR Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Aug2017_TriLepSR.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Aug2017_TriLepSR);

 Aug2017_TriLepSR::Aug2017_TriLepSR() : AnalyzerCore(), out_muons(0) {

   SetLogName("Aug2017_TriLepSR");
   Message("In Aug2017_TriLepSR constructor", INFO);
   InitialiseAnalysis();
 }


 void Aug2017_TriLepSR::InitialiseAnalysis() throw( LQError ) {
   
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

void Aug2017_TriLepSR::ExecuteEvents()throw( LQError ){

////Basic Infos///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Apply the gen weight 
   if(!isData) weight*=MCweight;
   FillHist("GenWeight", MCweight, 1, -2., 2., 4);


   m_logger<<DEBUG<<"RunNumber/Event Number = "<<eventbase->GetEvent().RunNumber()<<" : "<<eventbase->GetEvent().EventNumber()<<LQLogger::endmsg;
   m_logger<<DEBUG<<"isData = "<<isData<<LQLogger::endmsg;


   //Pileup Reweight / Signal Weight
   float pileup_reweight=1., pileup_reweight_systdown=1., pileup_reweight_systup=1.;
   float k_factor_weight=1., geneff_weight=1., gennorm_weight=1.;
   if(!k_isdata){
      pileup_reweight          = eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central); 
      pileup_reweight_systup   = pileup_reweight>0.? eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::up)  /pileup_reweight :0. ;
      pileup_reweight_systdown = pileup_reweight>0.? eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::down)/pileup_reweight :0. ;
      k_factor_weight          = GetKFactor();
      geneff_weight            = GenFilterEfficiency(k_sample_name);
      gennorm_weight           = SignalNorm(k_sample_name, 20.);
      FillHist("Weight_PU", pileup_reweight, 1., 0., 20., 200);
   }
   weight *= k_factor_weight*geneff_weight*gennorm_weight*pileup_reweight;

 
   //Total Event(MCweight+PUrw)////////////////
   FillCutFlow("NoCut", weight*pileup_reweight);


   bool EMuMu=false, TriMu=false;
   bool CutOpt=false, ObsExpComp=false, SRYield=false, MAWinOpt=false;
   bool SRDist=false, SigBkgKin=false, GenNormCheck=false, CheckSystSize=false, CutFlowCheck=false;
   bool DoubleMuon=false, ElectronMuon=false;
   bool SystRun=false, TestRun=false;
   TString Cycle="";
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("EMuMu"))         EMuMu         = true;
     else if(k_flags.at(i).Contains("TriMu"))         TriMu         = true;
     else if(k_flags.at(i).Contains("SystRun"))       SystRun       = true;
     else if(k_flags.at(i).Contains("CutOpt"))        CutOpt        = true;
     else if(k_flags.at(i).Contains("MAWinOpt"))      MAWinOpt      = true;
     else if(k_flags.at(i).Contains("ObsExpComp"))   {ObsExpComp    = true; Cycle="ObsExpComp";}
     else if(k_flags.at(i).Contains("SRYield"))      {SRYield       = true; Cycle="SRYield";}
     else if(k_flags.at(i).Contains("SRDist"))       {SRDist        = true; Cycle="SRDist";}
     else if(k_flags.at(i).Contains("SigBkgKin"))    {SigBkgKin     = true; Cycle="SigBkgKin";}
     else if(k_flags.at(i).Contains("CutFlowCheck")) {CutFlowCheck  = true; }
     else if(k_flags.at(i).Contains("GenNormCheck"))  GenNormCheck  = true;
     else if(k_flags.at(i).Contains("CheckSystSize")) CheckSystSize = true;
     else if(k_flags.at(i).Contains("TestRun"))       TestRun       = true;
     else if(k_flags.at(i).Contains("DoubleMuon"))    DoubleMuon    = true;
     else if(k_flags.at(i).Contains("ElectronMuon"))  ElectronMuon  = true;
   }

   if(GenNormCheck){
     std::vector<float> PDFWVec  = eventbase->GetEvent().PdfWeights();
     std::vector<float> ScaleWVec= eventbase->GetEvent().ScaleWeights();

     std::vector<float> LOPDFWVec;
     std::vector<float> NNPDF30LOAS0130WVec;
     LOPDFWVec.push_back(PDFWVec.at(0));
     LOPDFWVec.push_back(PDFWVec.at(332));
     LOPDFWVec.push_back(PDFWVec.at(435));
     LOPDFWVec.push_back(PDFWVec.at(1067));
     LOPDFWVec.push_back(PDFWVec.at(1068));
     LOPDFWVec.push_back(PDFWVec.at(1069));
     LOPDFWVec.push_back(PDFWVec.at(1070));
     LOPDFWVec.push_back(PDFWVec.at(1172)/PDFWVec.at(0)*PDFWVec.at(1070));
     LOPDFWVec.push_back(PDFWVec.at(1173)/PDFWVec.at(0)*PDFWVec.at(1070));
     for(int i=1071; i<=1171; i++){
       NNPDF30LOAS0130WVec.push_back(PDFWVec.at(i)/PDFWVec.at(0)*PDFWVec.at(1070));
     }

     for(int it_pdf=0; it_pdf<101; it_pdf++){
     //for(int it_pdf=0; it_pdf<102; it_pdf++){
       FillHist("Norm_PDF_M12to40", it_pdf+0.0001, weight*PDFWVec.at(it_pdf), 0., 102., 102);
       FillHist("Norm_PDF_M15"    , it_pdf+0.0001, weight*PDFWVec.at(it_pdf), 0., 102., 102);
       FillHist("Norm_PDF_M20"    , it_pdf+0.0001, weight*PDFWVec.at(it_pdf), 0., 102., 102);
       FillHist("Norm_PDF_M25"    , it_pdf+0.0001, weight*PDFWVec.at(it_pdf), 0., 102., 102);
       FillHist("Norm_PDF_M30"    , it_pdf+0.0001, weight*PDFWVec.at(it_pdf), 0., 102., 102);
       FillHist("Norm_PDF_M35"    , it_pdf+0.0001, weight*PDFWVec.at(it_pdf), 0., 102., 102);
     }
     for(int it_pdf=0; it_pdf<(int) NNPDF30LOAS0130WVec.size(); it_pdf++){
       FillHist("Norm_NNPDFLO_M12to40", it_pdf+0.0001, weight*NNPDF30LOAS0130WVec.at(it_pdf), 0., 102., 102);
       FillHist("Norm_NNPDFLO_M15"    , it_pdf+0.0001, weight*NNPDF30LOAS0130WVec.at(it_pdf), 0., 102., 102);
       FillHist("Norm_NNPDFLO_M20"    , it_pdf+0.0001, weight*NNPDF30LOAS0130WVec.at(it_pdf), 0., 102., 102);
       FillHist("Norm_NNPDFLO_M25"    , it_pdf+0.0001, weight*NNPDF30LOAS0130WVec.at(it_pdf), 0., 102., 102);
       FillHist("Norm_NNPDFLO_M30"    , it_pdf+0.0001, weight*NNPDF30LOAS0130WVec.at(it_pdf), 0., 102., 102);
       FillHist("Norm_NNPDFLO_M35"    , it_pdf+0.0001, weight*NNPDF30LOAS0130WVec.at(it_pdf), 0., 102., 102);
     }
     for(int it_pdf=0; it_pdf<(int) LOPDFWVec.size(); it_pdf++){
       FillHist("Norm_LOPDF_M12to40", it_pdf+0.0001, weight*LOPDFWVec.at(it_pdf), 0., 10., 10);
       FillHist("Norm_LOPDF_M15"    , it_pdf+0.0001, weight*LOPDFWVec.at(it_pdf), 0., 10., 10);
       FillHist("Norm_LOPDF_M20"    , it_pdf+0.0001, weight*LOPDFWVec.at(it_pdf), 0., 10., 10);
       FillHist("Norm_LOPDF_M25"    , it_pdf+0.0001, weight*LOPDFWVec.at(it_pdf), 0., 10., 10);
       FillHist("Norm_LOPDF_M30"    , it_pdf+0.0001, weight*LOPDFWVec.at(it_pdf), 0., 10., 10);
       FillHist("Norm_LOPDF_M35"    , it_pdf+0.0001, weight*LOPDFWVec.at(it_pdf), 0., 10., 10);
     }
   }

    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //Trigger Path of Analysis
   bool Pass_Trigger=false, Pass_TriggerBG=false, Pass_TriggerH=false;
   float trigger_ps_weight=1., trigger_period_weight=1.;
   float LumiBG=27.257618, LumiH=8.605696, LumiBH=35.863314;

   if(EMuMu){
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
   if(TriMu){
     if(  PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v")
        ||PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v") )    Pass_TriggerBG=true;
     if(  PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")
        ||PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ) Pass_TriggerH =true;

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
   if(SigBkgKin || CheckSystSize){ Pass_Trigger=true; trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi); }
   FillHist("Basic_TrigPSWeight", trigger_ps_weight, 1., 0., 1., 100);
   weight*=trigger_ps_weight*trigger_period_weight;

   //Trigger Cut
   if(!Pass_Trigger) return;
   FillCutFlow("TriggerCut", weight*pileup_reweight);
   /**********************************************************************************/

   //METFilterCut : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
   if(!PassMETFilter()) return;  FillCutFlow("EventCut", weight*pileup_reweight);

   //Vertex Cut
   //(vtx.ndof>4&&maxAbsZ<=0)||abs(vtx.z)<= 24)&&((maxd0 <=0)||abs(vtx.position.rho)<=2)&&!(vtx.isFake))
   if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; FillCutFlow("VertexCut", weight*pileup_reweight);
   


   //===========================================================================//
   // Objects Selection in Analysis                                             //
   //===========================================================================//

  
   //**PreSelCut*******************************************************************************************//
   //Intended for Code speed boosting up.
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(5.);                     eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, true);
     eventbase->GetElectronSel()->SetPt(5.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
     if(EMuMu)      { if( !(electronPreColl.size()>=1 && muonPreColl.size()>=2) ) return; }
     else if(TriMu) { if( !(muonPreColl.size()>=3) ) return; }
     FillCutFlow("PreSel", weight);
   //******************************************************************************************************//

   //Primary Object Collection
   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);

     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(5.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.2);                 eventbase->GetMuonSel()->SetdxySigMax(4.);
     eventbase->GetMuonSel()->SetBSdz(0.1);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.6);
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(5.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(4.);
     eventbase->GetMuonSel()->SetBSdz(0.05);
     eventbase->GetMuonSel()->SetChiNdof(4.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.2);
   std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,true);
   std::vector<snu::KMuon> muonColl;  if(k_running_nonprompt){ muonColl=muonLooseColl;} else{ muonColl=muonTightColl;}


     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_HctoWA_FAKELOOSE);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(5.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.4, 0.4);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(5.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronTightColl; eventbase->GetElectronSel()->Selection(electronTightColl);
   std::vector<snu::KElectron> electronColl;  if(!k_running_nonprompt){ electronColl=electronTightColl;} else{ electronColl=electronLooseColl;}


     bool LeptonVeto=false;
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
   std::vector<snu::KJet> jetNoVetoColl; eventbase->GetJetSel()->Selection(jetNoVetoColl, LeptonVeto, muonLooseColl, electronLooseColl);
   std::vector<snu::KJet> bjetNoVetoColl = SelBJets(jetNoVetoColl, "Medium");

   std::vector<snu::KJet> jetColl  = SkimJetColl(jetNoVetoColl,  electronLooseColl, muonLooseColl, "EleMuVeto");
   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");

   float met    = eventbase->GetEvent().MET(snu::KEvent::pfmet);
   float metphi = eventbase->GetEvent().METPhi();
   float met_x  = eventbase->GetEvent().PFMETx();
   float met_y  = eventbase->GetEvent().PFMETy();
   float Pzv, Pzv1, Pzv2;
   snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));
   //snu::KParticle v[4]; v[0].SetPx(met_x); v[0].SetPy(met_y);
   //                     v[1].SetPx(met_x); v[1].SetPy(met_y);

   int Nvtx=eventbase->GetEvent().nVertices();
   //------------------------------------------------------------------------------------------------------------------//
  

   //=====================================================//
   // Correction Factors                                  //
   //=====================================================//
   float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1., btag_sf=1.;
   float trigger_sf=1.;
   float fake_weight=1.; bool EventCand=false;

   /*This part is for boosting up speed.. SF part takes rather longer time than expected*/
   if     (EMuMu){ if(electronLooseColl.size()>=1 && muonLooseColl.size()>=2) EventCand=true; }
   else if(TriMu){ if(                  muonLooseColl.size()>=3             ) EventCand=true; }

   if(EventCand & !SystRun && !CutFlowCheck){
     if(!isData){
       if(EMuMu || TriMu){
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl);
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
    
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl);
         trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);

         btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);
         if(TriMu){
           float trigger_sf1 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v");
           float trigger_sf2 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v");
           trigger_sf    = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
         }
         if(EMuMu){
           float trigger_sf1 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
           float trigger_sf2 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
           trigger_sf    = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
         }
       }
     }
     else{
       //Perfect Prompt Ratio Approximation applied.(p=1), Applicable to generic number, combination of leptons under premise of high prompt rate.
       if(k_running_nonprompt){
         fake_weight = GetFakeWeight_Data(muonLooseColl, electronLooseColl, "POGTIsop6IPp2p1sig4" , "POGTIsop20IPp01p05sig4Chi4", "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", "TrkIsoVVLConeSUSY");
       }
     }
   }
   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*btag_sf*trigger_sf*fake_weight;
   //-----------------------------------------------------------------------------------------//



   //----------------------------------------------------------------------------------//
   //==================================================================================//
   /////// MAIN ANALYSIS CODE ///////////////////////////////////////////////////////////
   //==================================================================================//
   //----------------------------------------------------------------------------------//

 
   if(SystRun){
     // Syst Sel. and Syst Corr. -------------------------------------------------------------------//
       eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
       eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
       eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
       eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
       eventbase->GetMuonSel()->SetdxySigMax(4.);
       eventbase->GetMuonSel()->SetChiNdof(4.);
     std::vector<snu::KMuon> MuTMuEnUpColl; eventbase->GetMuonSel()->Selection(MuTMuEnUpColl,true);
       SetMuonResCorrection(MuTMuEnUpColl, "SystUp");

       eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
       eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
       eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
       eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
       eventbase->GetMuonSel()->SetdxySigMax(4.);
       eventbase->GetMuonSel()->SetChiNdof(4.);
     std::vector<snu::KMuon> MuTMuEnDownColl; eventbase->GetMuonSel()->Selection(MuTMuEnDownColl,true);
       SetMuonResCorrection(MuTMuEnDownColl, "SystDown");


       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
       eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
       eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
       eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetdxySigMax(4.);
       eventbase->GetElectronSel()->SetApplyConvVeto(true);
     std::vector<snu::KElectron> EleTElEnUpColl; eventbase->GetElectronSel()->Selection(EleTElEnUpColl, "SystUpElEn");

       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
       eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
       eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
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
     float id_weight_ele_sfup=1.    , id_weight_ele_sfdown=1.    , id_weight_mu_sfup=1.    , id_weight_mu_sfdown=1.;
     float id_weight_ele_ElEnup=1.  , reco_weight_ele_ElEnup=1.  , id_weight_mu_MuEnup=1.  , trk_weight_mu_MuEnup=1.;
     float id_weight_ele_ElEndown=1., reco_weight_ele_ElEndown=1., id_weight_mu_MuEndown=1., trk_weight_mu_MuEndown=1.;
     float btag_sf_JESup=1.  , btag_sf_JERup=1.  , btag_sf_LTagup=1.  , btag_sf_BCTagup=1.;
     float btag_sf_JESdown=1., btag_sf_JERdown=1., btag_sf_LTagdown=1., btag_sf_BCTagdown=1.;
     float fake_weight_FRup=1., fake_weight_FRdown=1.;
     float ZG_SF=1., conv_up=1., conv_down=1.;
     float xsec_up=1., xsec_down=1.;
     float trigger_sf=1., trigger_sf_up=1., trigger_sf_down=1.;
     if(!isData){
       id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl);
       reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
       id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl);
       trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
       btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

       id_weight_ele_sfup       = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl,  1);
       id_weight_ele_sfdown     = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl, -1);
       id_weight_mu_sfup        = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl,  1);
       id_weight_mu_sfdown      = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl, -1);

       id_weight_ele_ElEnup     = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", EleTElEnUpColl);
       reco_weight_ele_ElEnup   = mcdata_correction->ElectronRecoScaleFactor(EleTElEnUpColl);
       id_weight_mu_MuEnup      = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuTMuEnUpColl);
       trk_weight_mu_MuEnup     = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnUpColl);

       id_weight_ele_ElEndown   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", EleTElEnDownColl);
       reco_weight_ele_ElEndown = mcdata_correction->ElectronRecoScaleFactor(EleTElEnDownColl);
       id_weight_mu_MuEndown    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuTMuEnDownColl);
       trk_weight_mu_MuEndown   = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnDownColl);

       if(TriMu){
         float trigger_sf1          = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v");
         float trigger_sf2          = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v");
         float trigger_sf1_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v"   ,"Leg1SystUp");
         float trigger_sf2_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","Leg1SystUp");
         float trigger_sf1_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v"   ,"Leg2SystUp");
         float trigger_sf2_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","Leg2SystUp");
         float trigger_sf1_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v"   ,"Leg1SystDown");
         float trigger_sf2_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","Leg1SystDown");
         float trigger_sf1_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v"   ,"Leg2SystDown");
         float trigger_sf2_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","Leg2SystDown");
         float trigger_sf2_dzup     = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","DZSystUp");
         float trigger_sf2_dzdown   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","DZSystDown");
         float trigger_sf1_up   = trigger_sf1+sqrt(pow(trigger_sf1_leg1up  -trigger_sf1,2)+pow(trigger_sf1_leg2up  -trigger_sf1,2));
         float trigger_sf1_down = trigger_sf1-sqrt(pow(trigger_sf1_leg1down-trigger_sf1,2)+pow(trigger_sf1_leg2down-trigger_sf1,2));
         float trigger_sf2_up   = trigger_sf2+sqrt(pow(trigger_sf2_leg1up  -trigger_sf2,2)+pow(trigger_sf2_leg2up  -trigger_sf2,2)+pow(trigger_sf2_dzup  -trigger_sf2,2));
         float trigger_sf2_down = trigger_sf2-sqrt(pow(trigger_sf2_leg1down-trigger_sf2,2)+pow(trigger_sf2_leg2down-trigger_sf2,2)+pow(trigger_sf2_dzdown-trigger_sf2,2));

         trigger_sf      = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
         trigger_sf_up   = ((Pass_TriggerBG ? trigger_sf1_up:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_up:0.)*LumiH)/LumiBH;
         trigger_sf_down = ((Pass_TriggerBG ? trigger_sf1_down:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_down:0.)*LumiH)/LumiBH;
       }
       else if(EMuMu){
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

       xsec_up  +=GetXsecUncertainty(k_sample_name);
       xsec_down-=GetXsecUncertainty(k_sample_name);
       if(TriMu){
         if(k_sample_name.Contains("ZGto2LG")){ ZG_SF=0.86191; conv_up+=0.15; conv_down-=0.15; }
         else if(k_sample_name.Contains("TTG")){ conv_up+=0.5; conv_down-=0.5; }
       }
       if(EMuMu){
         if(k_sample_name.Contains("ZGto2LG")){ ZG_SF=0.963169; conv_up+=0.096; conv_down-=0.096; }
         else if(k_sample_name.Contains("TTG")){ conv_up+=0.5; conv_down-=0.5; }
       }
     }
     else if(k_running_nonprompt){
       fake_weight = GetFakeWeight_Data(muonLooseColl, electronLooseColl, "POGTIsop6IPp2p1sig4" , "POGTIsop20IPp01p05sig4Chi4", "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", "TrkIsoVVLConeSUSY");
     }
     weight*=ZG_SF;
       
     float systweight_central=weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight;
                                                                                                                                                        
     float systweight_Trigup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf_up*fake_weight;
     float systweight_ElIDup =weight*id_weight_ele_sfup  *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight;
     float systweight_MuIDup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu_sfup  *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight;
     float systweight_ElEnup =weight*id_weight_ele_ElEnup*reco_weight_ele_ElEnup*id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight;
     float systweight_MuEnup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu_MuEnup*trk_weight_mu_MuEnup*btag_sf        *trigger_sf   *fake_weight;
     float systweight_JESup  =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JESup  *trigger_sf   *fake_weight;
     float systweight_JERup  =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JERup  *trigger_sf   *fake_weight;
     float systweight_LTagup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_LTagup *trigger_sf   *fake_weight;
     float systweight_BCTagup=weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_BCTagup*trigger_sf   *fake_weight;
     float systweight_PUup   =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight*pileup_reweight_systup;
     float systweight_Xsecup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight*xsec_up;
     float systweight_Convup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight*conv_up;

     float systweight_Trigdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf_down*fake_weight;
     float systweight_ElIDdown =weight*id_weight_ele_sfdown  *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight;
     float systweight_MuIDdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu_sfdown  *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight;
     float systweight_ElEndown =weight*id_weight_ele_ElEndown*reco_weight_ele_ElEndown*id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight;
     float systweight_MuEndown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu_MuEndown*trk_weight_mu_MuEndown*btag_sf          *trigger_sf     *fake_weight;
     float systweight_JESdown  =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JESdown  *trigger_sf     *fake_weight;
     float systweight_JERdown  =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JERdown  *trigger_sf     *fake_weight;
     float systweight_LTagdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_LTagdown *trigger_sf     *fake_weight;
     float systweight_BCTagdown=weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_BCTagdown*trigger_sf     *fake_weight;
     float systweight_PUdown   =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight*pileup_reweight_systdown;
     float systweight_Xsecdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight*xsec_down;
     float systweight_Convdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight*conv_down;


     TString ChannelString = EMuMu? "EMuMu":"TriMu";
     //----------------------------------------------------------------------------------------------------------------------//


     DoSystRun(Cycle, ChannelString+"",
               electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
               systweight_central);
     if(!isData){
       DoSystRun(Cycle, ChannelString+"SystUpElID",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_ElIDup);
       DoSystRun(Cycle, ChannelString+"SystUpMuID",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_MuIDup);
       DoSystRun(Cycle, ChannelString+"SystUpTrig",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Trigup);
       DoSystRun(Cycle, ChannelString+"SystUpPU",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_PUup);
       DoSystRun(Cycle, ChannelString+"SystUpUncl",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi),
                 systweight_central);
       DoSystRun(Cycle, ChannelString+"SystUpJES",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESUpColl, bjetJESUpColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi),
                 systweight_JESup);
       DoSystRun(Cycle, ChannelString+"SystUpJER",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERUpColl, bjetJERUpColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi),
                 systweight_JERup);
       DoSystRun(Cycle, ChannelString+"SystUpBTag_L",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_LTagup);
       DoSystRun(Cycle, ChannelString+"SystUpBTag_BC",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_BCTagup);
       DoSystRun(Cycle, ChannelString+"SystUpElEn",
                 EleTElEnUpColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi),
                 systweight_ElEnup);
       DoSystRun(Cycle, ChannelString+"SystUpMuEn",
                 electronColl, electronLooseColl, MuTMuEnUpColl, muonLooseColl, jetColl, bjetColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi),
                 systweight_MuEnup);
       DoSystRun(Cycle, ChannelString+"SystUpXsec",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Xsecup);
       DoSystRun(Cycle, ChannelString+"SystUpConv",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Convup);



       DoSystRun(Cycle, ChannelString+"SystDownElID",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_ElIDdown);
       DoSystRun(Cycle, ChannelString+"SystDownMuID",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_MuIDdown);
       DoSystRun(Cycle, ChannelString+"SystDownTrig",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Trigdown);
       DoSystRun(Cycle, ChannelString+"SystDownPU",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_PUdown);
       DoSystRun(Cycle, ChannelString+"SystDownUncl",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi),
                 systweight_central);
       DoSystRun(Cycle, ChannelString+"SystDownJES",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESDownColl, bjetJESDownColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi),
                 systweight_JESdown);
       DoSystRun(Cycle, ChannelString+"SystDownJER",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERDownColl, bjetJERDownColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi),
                 systweight_JERdown);
       DoSystRun(Cycle, ChannelString+"SystDownBTag_L",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_LTagdown);
       DoSystRun(Cycle, ChannelString+"SystDownBTag_BC",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_BCTagdown);
       DoSystRun(Cycle, ChannelString+"SystDownElEn",
                 EleTElEnDownColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi),
                 systweight_ElEndown);
       DoSystRun(Cycle, ChannelString+"SystDownMuEn",
                 electronColl, electronLooseColl, MuTMuEnDownColl, muonLooseColl, jetColl, bjetColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi),
                 systweight_MuEndown);
       DoSystRun(Cycle, ChannelString+"SystDownXsec",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Xsecdown);
       DoSystRun(Cycle, ChannelString+"SystDownConv",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Convdown);
     }
   }
   if(CheckSystSize){
     
     for(int i=0; i<(int) muonColl.size(); i++){
       double dRelPt = mcdata_correction->GetRochesterMomentumWidth(muonColl.at(i));
       FillHist("Syst_MuRes", fabs(dRelPt), 1, 0., 1., 100);
     }
     for(int i=0; i<(int) electronColl.size(); i++){
       FillHist("Syst_ElRes", electronColl.at(i).PtShiftedUp()-1.  , 1, -1., 1., 200);
       FillHist("Syst_ElRes", electronColl.at(i).PtShiftedDown()-1., 1, -1., 1., 200);
     }
     for(int i=0; i<(int) jetColl.size(); i++){
       FillHist("Syst_JetES", jetColl.at(i).ScaledUpEnergy()-1.  , 1, -1., 1., 200);
       FillHist("Syst_JetES", jetColl.at(i).ScaledDownEnergy()-1., 1, -1., 1., 200);
     }
     for(int i=0; i<(int) jetColl.size(); i++){
       if(jetColl.at(i).IsMCSmeared()){
         FillHist("Syst_JetER", (jetColl.at(i).SmearedResUp()/jetColl.at(i).SmearedRes())-1.  , 1, -1., 1., 200);
         FillHist("Syst_JetER", (jetColl.at(i).SmearedResDown()/jetColl.at(i).SmearedRes())-1., 1, -1., 1., 200);
       }
       else{
         FillHist("Syst_JetER", jetColl.at(i).SmearedResUp()-1.  , 1, -1., 1., 200);
         FillHist("Syst_JetER", jetColl.at(i).SmearedResDown()-1., 1, -1., 1., 200);
       }
     }
     for(int i=0; i<(int) bjetColl.size(); i++){
       btag_sf          = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);
       float btag_sf_LTagup   = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpLTag")/btag_sf-1.;
       float btag_sf_BCTagup  = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpBCTag")/btag_sf-1.;
       float btag_sf_LTagdown = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownLTag")/btag_sf-1.;
       float btag_sf_BCTagdown= BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownBCTag")/btag_sf-1.;
       FillHist("Syst_BTag_L" , btag_sf_LTagup   , 1, -1., 1., 200);
       FillHist("Syst_BTag_L" , btag_sf_LTagdown , 1, -1., 1., 200);
       FillHist("Syst_BTag_BC", btag_sf_BCTagup  , 1, -1., 1., 200);
       FillHist("Syst_BTag_BC", btag_sf_BCTagdown, 1, -1., 1., 200);
     }

   }
   if(SRDist){
     if(EMuMu && !SystRun){
       CheckSRDist(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMu");
     }
     if(TriMu && !SystRun){
       CheckSRDist(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMu");
     }
   }
   if(SigBkgKin){
     CheckSigBkgKinematics(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "");
   }
   if(CutFlowCheck){
     if(EMuMu && (!SystRun)){
        CheckCutflow(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMu");
     }
     if(TriMu && (!SystRun)){
        CheckCutflow(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMu");
     }
   }
   if(SRYield){
     if(EMuMu && (!SystRun)){
        CheckSRYield(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMu");
     }
//     if(SystRun){
//
//       // Syst Sel. and Syst Corr. ----------------------------------------------------------------------------------------------------//
//         eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
//         eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
//         eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
//         eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
//         eventbase->GetMuonSel()->SetdxySigMax(4.);
//         eventbase->GetMuonSel()->SetChiNdof(4.);
//       std::vector<snu::KMuon> MuTMuEnUpColl; eventbase->GetMuonSel()->Selection(MuTMuEnUpColl,true, "SystUpMuEn");
//
//         eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
//         eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
//         eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
//         eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
//         eventbase->GetMuonSel()->SetdxySigMax(4.);
//         eventbase->GetMuonSel()->SetChiNdof(4.);
//       std::vector<snu::KMuon> MuTMuEnDownColl; eventbase->GetMuonSel()->Selection(MuTMuEnDownColl,true, "SystDownMuEn");
//
//         eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
//         eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
//         eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
//         eventbase->GetElectronSel()->SetBETrRegIncl(false);
//         eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
//         eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
//         eventbase->GetElectronSel()->SetdxySigMax(4.);
//         eventbase->GetElectronSel()->SetApplyConvVeto(true);
//       std::vector<snu::KElectron> EleTElEnUpColl; eventbase->GetElectronSel()->Selection(EleTElEnUpColl, "SystUpElEn");
//
//         eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
//         eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
//         eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
//         eventbase->GetElectronSel()->SetBETrRegIncl(false);
//         eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
//         eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
//         eventbase->GetElectronSel()->SetdxySigMax(4.);
//         eventbase->GetElectronSel()->SetApplyConvVeto(true);
//       std::vector<snu::KElectron> EleTElEnDownColl; eventbase->GetElectronSel()->Selection(EleTElEnDownColl, "SystDownElEn");
//
//       LeptonVeto=true;
//         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
//         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
//       std::vector<snu::KJet> jetJESUpColl; eventbase->GetJetSel()->Selection(jetJESUpColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystUpJES");
//         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
//         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
//       std::vector<snu::KJet> jetJERUpColl; eventbase->GetJetSel()->Selection(jetJERUpColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystUpJER");
//         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
//         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
//       std::vector<snu::KJet> jetJESDownColl; eventbase->GetJetSel()->Selection(jetJESDownColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystDownJES");
//         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
//         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
//       std::vector<snu::KJet> jetJERDownColl; eventbase->GetJetSel()->Selection(jetJERDownColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystDownJER");
//
//       std::vector<snu::KJet> bjetJESUpColl   = SelBJets(jetJESUpColl  , "Medium");
//       std::vector<snu::KJet> bjetJERUpColl   = SelBJets(jetJERUpColl  , "Medium");
//       std::vector<snu::KJet> bjetJESDownColl = SelBJets(jetJESDownColl, "Medium");
//       std::vector<snu::KJet> bjetJERDownColl = SelBJets(jetJERDownColl, "Medium");
//
//       float met_JESup     = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetEn       , snu::KEvent::up);
//       float met_JERup     = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetRes      , snu::KEvent::up);
//       float met_Unclup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::Unclustered , snu::KEvent::up);
//       float met_ElEnup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::ElectronEn  , snu::KEvent::up);
//       float met_MuEnup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::MuonEn      , snu::KEvent::up);
//       float met_JESdown   = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetEn       , snu::KEvent::down);
//       float met_JERdown   = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetRes      , snu::KEvent::down);
//       float met_Uncldown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::Unclustered , snu::KEvent::down);
//       float met_ElEndown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::ElectronEn  , snu::KEvent::down);
//       float met_MuEndown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::MuonEn      , snu::KEvent::down);
//
//       //Scale Factors--------------------------------------------------------------------------------//      
//       float systweight=weight; //Lumi weight+Prescale weight+Period weight+PU weight applied by default;
//       float id_weight_ele_ElEnup=1.  , reco_weight_ele_ElEnup=1.  , id_weight_mu_MuEnup=1.  , trk_weight_mu_MuEnup=1.;
//       float id_weight_ele_ElEndown=1., reco_weight_ele_ElEndown=1., id_weight_mu_MuEndown=1., trk_weight_mu_MuEndown=1.;
//       float btag_sf_JESup=1.  , btag_sf_JERup=1.  , btag_sf_LTagup=1.  , btag_sf_BCTagup=1.;
//       float btag_sf_JESdown=1., btag_sf_JERdown=1., btag_sf_LTagdown=1., btag_sf_BCTagdown=1.;
//       float fake_weight_FRup=1., fake_weight_FRdown=1.;
//       float trigger_sf=1., trigger_sf_up=1., trigger_sf_down=1.;
//       TString Str_TrigBase="";
//       if(!isData){
//         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl);
//         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
//         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl);
//         trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
//         btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);
//
//         id_weight_ele_ElEnup   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", EleTElEnUpColl);
//         reco_weight_ele_ElEnup = mcdata_correction->ElectronRecoScaleFactor(EleTElEnUpColl);
//         id_weight_mu_MuEnup    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuTMuEnUpColl);
//         trk_weight_mu_MuEnup   = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnUpColl);
//
//         id_weight_ele_ElEndown   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", EleTElEnDownColl);
//         reco_weight_ele_ElEndown = mcdata_correction->ElectronRecoScaleFactor(EleTElEnDownColl);
//         id_weight_mu_MuEndown    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuTMuEnDownColl);
//         trk_weight_mu_MuEndown   = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnDownColl);
//
//         if     (EMuMu) Str_TrigBase="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL";
//         else if(TriMu) Str_TrigBase="HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL";
//        
//         float trigger_sf1          = mcdata_correction->GetTriggerSF(electronColl, muonColl, Str_TrigBase+"_v");
//         float trigger_sf2          = mcdata_correction->GetTriggerSF(electronColl, muonColl, Str_TrigBase+"_DZ_v");
//         float trigger_sf1_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, Str_TrigBase+"_v"   ,"Leg1SystUp");
//         float trigger_sf2_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, Str_TrigBase+"_DZ_v","Leg1SystUp");
//         float trigger_sf1_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, Str_TrigBase+"_v"   ,"Leg2SystUp");
//         float trigger_sf2_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, Str_TrigBase+"_DZ_v","Leg2SystUp");
//         float trigger_sf1_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, Str_TrigBase+"_v"   ,"Leg1SystDown");
//         float trigger_sf2_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, Str_TrigBase+"_DZ_v","Leg1SystDown");
//         float trigger_sf1_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, Str_TrigBase+"_v"   ,"Leg2SystDown");
//         float trigger_sf2_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, Str_TrigBase+"_DZ_v","Leg2SystDown");
//         float trigger_sf2_dzup     = mcdata_correction->GetTriggerSF(electronColl, muonColl, Str_TrigBase+"_DZ_v","DZSystUp");
//         float trigger_sf2_dzdown   = mcdata_correction->GetTriggerSF(electronColl, muonColl, Str_TrigBase+"_DZ_v","DZSystDown");
//         float trigger_sf1_up   = trigger_sf1+sqrt(pow(trigger_sf1_leg1up  -trigger_sf1,2)+pow(trigger_sf1_leg2up  -trigger_sf1,2));
//         float trigger_sf1_down = trigger_sf1-sqrt(pow(trigger_sf1_leg1down-trigger_sf1,2)+pow(trigger_sf1_leg2down-trigger_sf1,2));
//         float trigger_sf2_up   = trigger_sf2+sqrt(pow(trigger_sf2_leg1up  -trigger_sf2,2)+pow(trigger_sf2_leg2up  -trigger_sf2,2)+pow(trigger_sf2_dzup  -trigger_sf2,2));
//         float trigger_sf2_down = trigger_sf2-sqrt(pow(trigger_sf2_leg1down-trigger_sf2,2)+pow(trigger_sf2_leg2down-trigger_sf2,2)+pow(trigger_sf2_dzdown-trigger_sf2,2));
//
//         trigger_sf      = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
//         trigger_sf_up   = ((Pass_TriggerBG ? trigger_sf1_up:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_up:0.)*LumiH)/LumiBH;
//         trigger_sf_down = ((Pass_TriggerBG ? trigger_sf1_down:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_down:0.)*LumiH)/LumiBH;
//         
//         btag_sf_LTagup = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpLTag");
//         btag_sf_BCTagup= BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpBCTag");
//         btag_sf_JESup  = BTagScaleFactor_1a(jetJESUpColl, snu::KJet::CSVv2, snu::KJet::Medium);
//         btag_sf_JERup  = BTagScaleFactor_1a(jetJERUpColl, snu::KJet::CSVv2, snu::KJet::Medium);
//
//         btag_sf_LTagdown = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownLTag");
//         btag_sf_BCTagdown= BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownBCTag");
//         btag_sf_JESdown  = BTagScaleFactor_1a(jetJESDownColl, snu::KJet::CSVv2, snu::KJet::Medium);
//         btag_sf_JERdown  = BTagScaleFactor_1a(jetJERDownColl, snu::KJet::CSVv2, snu::KJet::Medium);
//       }
//       else if(k_running_nonprompt){
//         fake_weight = GetFakeWeight(muonLooseColl, electronLooseColl, "POGTIsop6IPp2p1sig4" , "POGTIsop20IPp01p05sig4Chi4", "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", bjetNoVetoColl, "TrkIsoVVLConeSUSY");
//       }
//
//       float systweight_central=weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight *trigger_sf;
//
//       float systweight_Trigup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight *trigger_sf_up;
//       float systweight_ElEnup =weight*id_weight_ele_ElEnup*reco_weight_ele_ElEnup*id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight *trigger_sf;
//       float systweight_MuEnup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu_MuEnup*trk_weight_mu_MuEnup*btag_sf        *fake_weight *trigger_sf;
//       float systweight_JESup  =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JESup  *fake_weight *trigger_sf;
//       float systweight_JERup  =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JERup  *fake_weight *trigger_sf;
//       float systweight_LTagup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_LTagup *fake_weight *trigger_sf;
//       float systweight_BCTagup=weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_BCTagup*fake_weight *trigger_sf;
//       float systweight_PUup   =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight*pileup_reweight_systup *trigger_sf;
//       float systweight_FRup   =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *fake_weight_FRup *trigger_sf;
//
//       float systweight_Trigdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *fake_weight *trigger_sf_down;
//       float systweight_ElEndown =weight*id_weight_ele_ElEndown*reco_weight_ele_ElEndown*id_weight_mu         *trk_weight_mu         *btag_sf          *fake_weight *trigger_sf;
//       float systweight_MuEndown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu_MuEndown*trk_weight_mu_MuEndown*btag_sf          *fake_weight *trigger_sf;
//       float systweight_JESdown  =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JESdown  *fake_weight *trigger_sf;
//       float systweight_JERdown  =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JERdown  *fake_weight *trigger_sf;
//       float systweight_LTagdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_LTagdown *fake_weight *trigger_sf;
//       float systweight_BCTagdown=weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_BCTagdown*fake_weight *trigger_sf;
//       float systweight_PUdown   =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *fake_weight*pileup_reweight_systdown *trigger_sf;
//       float systweight_FRdown   =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *fake_weight_FRdown *trigger_sf;
//
//       //----------------------------------------------------------------------------------------------------------------------//
//
//       TString ChannelString = EMuMu? "EMuMu":"TriMu";
//       DoSystRun("SRYield", ChannelString+"",
//                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
//                 systweight_central);
//       if( !isData ){
//         DoSystRun("SRYield", ChannelString+"SystUpTrig",
//                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
//                   systweight_Trigup);
//         DoSystRun("SRYield", ChannelString+"SystUpPU",
//                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
//                   systweight_PUup);
//         DoSystRun("SRYield", ChannelString+"SystUpUncl",
//                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi),
//                   systweight_central);
//         DoSystRun("SRYield", ChannelString+"SystUpJES",
//                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESUpColl, bjetJESUpColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi),
//                   systweight_JESup);
//         DoSystRun("SRYield", ChannelString+"SystUpJER",
//                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERUpColl, bjetJERUpColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi),
//                   systweight_JERup);
//         DoSystRun("SRYield", ChannelString+"SystUpBTag_L",
//                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
//                   systweight_LTagup);
//         DoSystRun("SRYield", ChannelString+"SystUpBTag_BC",
//                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
//                   systweight_BCTagup);
//         DoSystRun("SRYield", ChannelString+"SystUpElEn",
//                   EleTElEnUpColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi),
//                   systweight_ElEnup);
//         DoSystRun("SRYield", ChannelString+"SystUpMuEn",
//                   electronColl, electronLooseColl, MuTMuEnUpColl, muonLooseColl, jetColl, bjetColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi),
//                   systweight_MuEnup);
//
//
//         DoSystRun("SRYield", ChannelString+"SystDownTrig",
//                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
//                   systweight_Trigdown);
//         DoSystRun("SRYield", ChannelString+"SystDownPU",
//                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
//                   systweight_PUdown);
//         DoSystRun("SRYield", ChannelString+"SystDownUncl",
//                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi),
//                   systweight_central);
//         DoSystRun("SRYield", ChannelString+"SystDownJES",
//                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESDownColl, bjetJESDownColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi),
//                   systweight_JESdown);
//         DoSystRun("SRYield", ChannelString+"SystDownJER",
//                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERDownColl, bjetJERDownColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi),
//                   systweight_JERdown);
//         DoSystRun("SRYield", ChannelString+"SystDownBTag_L",
//                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
//                   systweight_LTagdown);
//         DoSystRun("SRYield", ChannelString+"SystDownBTag_BC",
//                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
//                   systweight_BCTagdown);
//         DoSystRun("SRYield", ChannelString+"SystDownElEn",
//                   EleTElEnDownColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi),
//                   systweight_ElEndown);
//         DoSystRun("SRYield", ChannelString+"SystDownMuEn",
//                   electronColl, electronLooseColl, MuTMuEnDownColl, muonLooseColl, jetColl, bjetColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi),
//                   systweight_MuEndown);
//
//       }
//       else if( isData && k_running_nonprompt ){
////         DoSystRun("SRYield", ChannelString+"SystUpFR",
////                  electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
////                  systweight_FRup);
////         DoSystRun("SRYield", ChannelString+"SystDownFR",
////                  electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
////                  systweight_FRdown);
//       }
//     }//End of SystRun
   }
   if(CutOpt){

     if(MAWinOpt){
         eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
         eventbase->GetMuonSel()->SetPt(5.);                    eventbase->GetMuonSel()->SetEta(2.4);
         eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
         eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
         eventbase->GetMuonSel()->SetdxySigMax(4.);
         eventbase->GetMuonSel()->SetChiNdof(4.);
       std::vector<snu::KMuon> MuTMuEnUpColl; eventbase->GetMuonSel()->Selection(MuTMuEnUpColl,true, "SystUpMuEn");

         eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
         eventbase->GetMuonSel()->SetPt(5.);                    eventbase->GetMuonSel()->SetEta(2.4);
         eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
         eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
         eventbase->GetMuonSel()->SetdxySigMax(4.);
         eventbase->GetMuonSel()->SetChiNdof(4.);
       std::vector<snu::KMuon> MuTMuEnDownColl; eventbase->GetMuonSel()->Selection(MuTMuEnDownColl,true, "SystDownMuEn");


       if(EMuMu){
         OptimizeMACut(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight,
                       "", "EMuMu");
         if(!isData){
           OptimizeMACut(MuTMuEnUpColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight,
                         "_systupMuEn", "EMuMu");
           OptimizeMACut(MuTMuEnDownColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight,
                         "_systdownMuEn", "EMuMu");
         }
       }
       if(TriMu){
         OptimizeMACut(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight,
                       "", "TriMu");
         if(!isData){
           OptimizeMACut(MuTMuEnUpColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight,
                         "_systupMuEn", "TriMu");
           OptimizeMACut(MuTMuEnDownColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight,
                         "_systdownMuEn", "TriMu");
         }
       }
     }
     else{

     if(EMuMu){
       if( !(electronLooseColl.size()==1 && muonLooseColl.size()==2) ) return;
       if( !(electronColl.size()==1      && muonColl.size()==2) ) return;
       if( !(muonColl.at(0).Charge()!=muonColl.at(1).Charge()) ) return;
       if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>10 && muonColl.at(1).Pt()>5) ) return;
       if( fabs(electronColl.at(0).Eta())>2.5 ) return;
  
       float Mmumu=(muonColl.at(0)+muonColl.at(1)).M();
       float MTW = sqrt(2)*sqrt(met*electronColl.at(0).Pt()-met_x*electronColl.at(0).Px()-met_y*electronColl.at(0).Py());
       float M3l=(electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M();
       if(Mmumu<12) return;

       const int NJetCut=2;     int JetCuts[NJetCut]     ={2,3};
       const int NBJetCut=1;    int BJetCuts[NBJetCut]   ={1};
       const int NZVetoCut=1;   bool ZVetoCuts[NZVetoCut]={true};

       const int NPtElCuts=2;   float ElPtCuts[NPtElCuts]    ={25., 30.};
       const int NPtMu1Cuts=3;  float Mu1PtCuts[NPtMu1Cuts]  ={10., 15., 20.};
       const int NPtMu2Cuts=5;  float Mu2PtCuts[NPtMu2Cuts]  ={5., 7., 10., 12., 15.};
       const int NPtJetCuts=3;  float JetPtCuts[NPtJetCuts]  ={25., 30., 35.};
       const int NMETCuts=3;    float METCuts[NMETCuts]      ={0., 20., 30.};

                                                     
       FillHist("Count_PreSel", 0., weight, 0., 1., 1);       
       for(int it_Z=0; it_Z<NZVetoCut; it_Z++){
         if(ZVetoCuts[it_Z] && fabs(Mmumu-91.2)<10) continue;
         FillHist("Count_ZVeto", it_Z, weight, 0., 2., 2);
       }
       for(int it_nb=0; it_nb<NBJetCut; it_nb++){
         if((int) bjetColl.size()<BJetCuts[it_nb]) continue;
         FillHist("Count_NBCut", it_nb, weight, 0., 2., 2);
       }
       for(int it_nj=0; it_nj<NJetCut; it_nj++){
         if((int) jetColl.size()<JetCuts[it_nj]) continue;
         FillHist("Count_NJCut", it_nj, weight, 0., 4., 4);
       }
       for(int it_pte=0; it_pte<NPtElCuts; it_pte++){
         if( electronColl.at(0).Pt()<ElPtCuts[it_pte] ) continue;       
         FillHist("Count_PtElCut", it_pte, weight, 0., 2., 2);
       }
       for(int it_ptmu1=0; it_ptmu1<NPtMu1Cuts; it_ptmu1++){
         if( muonColl.at(0).Pt()<Mu1PtCuts[it_ptmu1] ) continue;
         FillHist("Count_PtMu1Cut", it_ptmu1, weight, 0., 3., 3);
       }
       for(int it_ptmu2=0; it_ptmu2<NPtMu2Cuts; it_ptmu2++){
         if( muonColl.at(1).Pt()<Mu2PtCuts[it_ptmu2] ) continue;
         FillHist("Count_PtMu2Cut", it_ptmu2, weight, 0., 5., 5);
       }
       for(int it_ptj=0; it_ptj<NPtJetCuts; it_ptj++){
         if((int)jetColl.size()>1 && jetColl.at(1).Pt()<JetPtCuts[it_ptj]) continue;
         FillHist("Count_PtjCut", it_ptj, weight, 0., 4., 4);
       }


       int it_Sel=0;
       for(int it_Z=0 ;    it_Z<NZVetoCut;        it_Z++){//1
       for(int it_nb=0;    it_nb<NBJetCut;        it_nb++){//2
       for(int it_nj=0;    it_nj<NJetCut ;        it_nj++){//3
       for(int it_pte=0;   it_pte<NPtElCuts;      it_pte++){//4
       for(int it_ptmu1=0; it_ptmu1<NPtMu1Cuts;   it_ptmu1++){//5
       for(int it_ptmu2=0; it_ptmu2<NPtMu2Cuts;   it_ptmu2++){//6
       for(int it_ptj=0;   it_ptj<NPtJetCuts;     it_ptj++){//7
       for(int it_met=0;   it_met<NMETCuts;       it_met++){//8
         it_Sel++;

         //PassCuts----------------------------------------------//
         if( !(electronColl.at(0).Pt()>ElPtCuts[it_pte] && muonColl.at(0).Pt()>Mu1PtCuts[it_ptmu1] && muonColl.at(1).Pt()>Mu2PtCuts[it_ptmu2]) ) continue;
         if( Mu1PtCuts[it_ptmu1]<Mu2PtCuts[it_ptmu2] ) continue;
         if( ZVetoCuts[it_Z] && fabs(Mmumu-91.2)<10 ) continue;
         if( (int) jetColl.size()<JetCuts[it_nj] || (int) bjetColl.size()<BJetCuts[it_nb] ) continue;
         if( jetColl.at(JetCuts[it_nj]-1).Pt()<JetPtCuts[it_ptj]   ) continue;
         if( bjetColl.at(BJetCuts[it_nb]-1).Pt()<JetPtCuts[it_ptj] ) continue;
         if( met<METCuts[it_met] ) continue;
         //------------------------------------------------------//

         if(Mmumu>12 && Mmumu<40){
           FillHist("Count_AllCut", it_Sel-1, weight, 0., 10000., 10000);
         }
       }//8
       }//7
       }//6
       }//5
       }//4
       }//3
       }//2
       }//1
       
     }//End of EMuMu
     if(TriMu){
       if( !(muonLooseColl.size()==3 && electronLooseColl.size()==0) ) return;
       if( !(muonColl.size()==3) ) return;
       if( fabs(SumCharge(muonColl))!=1 ) return;
       if( !(muonColl.at(0).Pt()>20 && muonColl.at(1).Pt()>10 && muonColl.at(2).Pt()>5) ) return;

       int   IdxOS  = TriMuChargeIndex(muonColl, "OS");
       int   IdxSS1 = TriMuChargeIndex(muonColl, "SS1");
       int   IdxSS2 = TriMuChargeIndex(muonColl, "SS2");
       float MOSSS1 = (muonColl.at(IdxOS)+muonColl.at(IdxSS1)).M();
       float MOSSS2 = (muonColl.at(IdxOS)+muonColl.at(IdxSS2)).M();
       if( !(MOSSS1>12 && MOSSS2>12) ) return;


       const int NJetCut=2;     int JetCuts[NJetCut]     ={2,3};
       const int NBJetCut=1;    int BJetCuts[NBJetCut]   ={1};
       const int NZVetoCut=2;   bool ZVetoCuts[NZVetoCut]={true, false};

       const int NPtMu1Cuts=3;  float Mu1PtCuts[NPtMu1Cuts]={20., 25., 30.};
       const int NPtMu2Cuts=3;  float Mu2PtCuts[NPtMu2Cuts]={10., 12., 15.};
       const int NPtMu3Cuts=5;  float Mu3PtCuts[NPtMu3Cuts]={5., 7., 10., 12., 15.};
       const int NPtJetCuts=3;  float JetPtCuts[NPtJetCuts]={25., 30., 35.};
       const int NMETCuts=3;    float METCuts[NMETCuts]    ={0., 20., 30.};


       FillHist("Count_PreSel", 0., weight, 0., 1., 1);       

       for(int it_Z=0; it_Z<NZVetoCut; it_Z++){
         if(ZVetoCuts[it_Z] && (fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10)) continue;
         FillHist("Count_ZVeto", it_Z, weight, 0., 2., 2);
       }
       for(int it_nb=0; it_nb<NBJetCut; it_nb++){
         if((int) bjetColl.size()<BJetCuts[it_nb]) continue;
         FillHist("Count_NBCut", it_nb, weight, 0., 2., 2);
       }
       for(int it_nj=0; it_nj<NJetCut; it_nj++){
         if((int) jetColl.size()<JetCuts[it_nj]) continue;
         FillHist("Count_NJCut", it_nj, weight, 0., 4., 4);
       }
       for(int it_ptmu1=0; it_ptmu1<NPtMu1Cuts; it_ptmu1++){
         if( muonColl.at(0).Pt()<Mu1PtCuts[it_ptmu1] ) continue;
         FillHist("Count_PtMu1Cut", it_ptmu1, weight, 0., 3., 3);
       }
       for(int it_ptmu2=0; it_ptmu2<NPtMu2Cuts; it_ptmu2++){
         if( muonColl.at(1).Pt()<Mu2PtCuts[it_ptmu2] ) continue;
         FillHist("Count_PtMu2Cut", it_ptmu2, weight, 0., 5., 5);
       }
       for(int it_ptmu3=0; it_ptmu3<NPtMu3Cuts; it_ptmu3++){
         if( muonColl.at(2).Pt()<Mu3PtCuts[it_ptmu3] ) continue;       
         FillHist("Count_PtMu3Cut", it_ptmu3, weight, 0., 2., 2);
       }
       for(int it_ptj=0; it_ptj<NPtJetCuts; it_ptj++){
         if((int) jetColl.size()>1 && jetColl.at(1).Pt()<JetPtCuts[it_ptj]) continue;
         FillHist("Count_PtjCut", it_ptj, weight, 0., 4., 4);
       }


       int it_Sel=0;
       for(int it_Z=0 ;    it_Z<NZVetoCut;      it_Z++){//1
       for(int it_nb=0;    it_nb<NBJetCut;      it_nb++){//2
       for(int it_nj=0;    it_nj<NJetCut ;      it_nj++){//3
       for(int it_ptmu1=0; it_ptmu1<NPtMu1Cuts; it_ptmu1++){//4
       for(int it_ptmu2=0; it_ptmu2<NPtMu2Cuts; it_ptmu2++){//5
       for(int it_ptmu3=0; it_ptmu3<NPtMu3Cuts; it_ptmu3++){//6
       for(int it_ptj=0;   it_ptj<NPtJetCuts;   it_ptj++){//7
       for(int it_met=0;   it_met<NMETCuts;     it_met++){//8
         it_Sel++;


         //PassCuts----------------------------------------------//
         if( !(muonColl.at(0).Pt()>Mu1PtCuts[it_ptmu1] && muonColl.at(1).Pt()>Mu2PtCuts[it_ptmu2] && muonColl.at(2).Pt()>Mu3PtCuts[it_ptmu3]) ) continue;
         if( Mu1PtCuts[it_ptmu1]<Mu2PtCuts[it_ptmu2] || Mu2PtCuts[it_ptmu2]<Mu3PtCuts[it_ptmu3] ) continue;
         if( ZVetoCuts[it_Z] && (fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10) ) continue;
         if( (int) jetColl.size()<JetCuts[it_nj] || (int) bjetColl.size()<BJetCuts[it_nb] ) continue;
         if( jetColl.at(JetCuts[it_nj]-1).Pt()<JetPtCuts[it_ptj] ) continue;
         if( bjetColl.at(BJetCuts[it_nb]-1).Pt()<JetPtCuts[it_ptj] ) continue;
         if( met<METCuts[it_met] ) continue;
         //------------------------------------------------------//

         if(MOSSS1>12 && MOSSS1<40){
           FillHist("Count_AllCut1", it_Sel-1, weight, 0., 10000., 10000);
         }
         if(MOSSS2>12 && MOSSS2<40){
           FillHist("Count_AllCut2", it_Sel-1, weight, 0., 10000., 10000);
         }

       }//8
       }//7
       }//6
       }//5
       }//4
       }//3
       }//2
       }//1
 
     }//End of TriMu
     }//End of Not Mass Window Opt
   }//End of CutOpt
/////////////////////////////////////////////////////////////////////////////////// 

return;

}// End of execute event loop
  


void Aug2017_TriLepSR::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Aug2017_TriLepSR::BeginCycle() throw( LQError ){
  
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

Aug2017_TriLepSR::~Aug2017_TriLepSR() {
  
  Message("In Aug2017_TriLepSR Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}


void Aug2017_TriLepSR::OptimizeMACut(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){

  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");
  if(EMuMu){
    if( !(EleLColl.size()==1 && MuLColl.size()==2) ) return;
    if( !(EleTColl.size()==1 && MuTColl.size()==2) ) return;
    if( !(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) ) return;
    if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && MuTColl.at(1).Pt()>10) ) return;
    if( fabs(EleTColl.at(0).Eta())>2.5 ) return;
  
    float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
    if( Mmumu<12 || fabs(Mmumu-91.2)<10 ) return;
  
    if(BJetColl.size()==0) return;
    if(JetColl.size()<2) return;
  
    const int Nwindow=64;
    float MwidthCuts[Nwindow]={10., 9., 8., 7.,
                               6.0, 5.9, 5.8, 5.7, 5.6, 5.5, 5.4, 5.3, 5.2, 5.1,
                               5.0, 4.9, 4.8, 4.7, 4.6, 4.5, 4.4, 4.3, 4.2, 4.1,
                               4.0, 3.9, 3.8, 3.7, 3.6, 3.5, 3.4, 3.3, 3.2, 3.1,
                               3.0, 2.9, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 
                               2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1,
                               1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};
    
    FillHist("NCount_Win"+Label, 0., weight, 0., 1., 1);
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-15.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM15"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-25.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM25"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-35.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM35"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }

    if(!isData){
      std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
      std::vector<snu::KMuon> MuAColl;
        for(int i=0; i<(int) MuTColl.size(); i++){
          int LepType=GetLeptonType(MuTColl.at(i),truthColl);
          if(LepType==2) MuAColl.push_back(MuTColl.at(i));
        }

      if(MuAColl.size()>1){
        FillHist("NCount_TrueWin"+Label, 0., weight, 0., 1., 1);
        float Mmumu_A=(MuAColl.at(0)+MuAColl.at(1)).M();
        FillHist("Mmumu_A"+Label, Mmumu_A, weight, 10., 40., 300);
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-15.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM15"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-25.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM25"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-35.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM35"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
      }
    }

    if(JetColl.size()<3) return;

    FillHist("NCount_3j_Win"+Label, 0., weight, 0., 1., 1);
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-15.)<MwidthCuts[it_win]){
        FillHist("NCount_3j_WinM15"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-25.)<MwidthCuts[it_win]){
        FillHist("NCount_3j_WinM25"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-35.)<MwidthCuts[it_win]){
        FillHist("NCount_3j_WinM35"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }

  }
  if(TriMu){
    if( !(EleLColl.size()==0 && MuLColl.size()==3) ) return;
    if( !(EleTColl.size()==0 && MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 && MuTColl.at(2).Pt()>10) ) return;
  
    int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    float Mmumu  = MOSSS2;
    if( !(MOSSS1>12 && MOSSS2>12) ) return;
    if( fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10 ) return;

    if(BJetColl.size()==0) return;
    if(JetColl.size()<2) return;
  
    const int Nwindow=64;
    float MwidthCuts[Nwindow]={10., 9., 8., 7.,
                               6.0, 5.9, 5.8, 5.7, 5.6, 5.5, 5.4, 5.3, 5.2, 5.1,
                               5.0, 4.9, 4.8, 4.7, 4.6, 4.5, 4.4, 4.3, 4.2, 4.1,
                               4.0, 3.9, 3.8, 3.7, 3.6, 3.5, 3.4, 3.3, 3.2, 3.1,
                               3.0, 2.9, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 
                               2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1,
                               1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};
    
    FillHist("NCount_Win"+Label, 0., weight, 0., 1., 1);
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-15.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM15"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-25.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM25"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-35.)<MwidthCuts[it_win]){
        FillHist("NCount_WinM35"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }

    if(!isData){
      std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
      std::vector<snu::KMuon> MuAColl;
        for(int i=0; i<(int) MuTColl.size(); i++){
          int LepType=GetLeptonType(MuTColl.at(i),truthColl);
          if(LepType==2) MuAColl.push_back(MuTColl.at(i));
        }

      if(MuAColl.size()>1){
        FillHist("NCount_TrueWin"+Label, 0., weight, 0., 1., 1);
        float Mmumu_A=(MuAColl.at(0)+MuAColl.at(1)).M();
        FillHist("Mmumu_A"+Label, Mmumu_A, weight, 10., 40., 300);
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-15.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM15"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-25.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM25"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
        for(int it_win=0; it_win<Nwindow; it_win++){
          if(fabs(Mmumu_A-35.)<MwidthCuts[it_win]){
            FillHist("NCount_TrueWinM35"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
          }
        }
      }
    }


    if(JetColl.size()<3) return;

    FillHist("NCount_3j_Win"+Label, 0., weight, 0., 1., 1);
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-15.)<MwidthCuts[it_win]){
        FillHist("NCount_3j_WinM15"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-25.)<MwidthCuts[it_win]){
        FillHist("NCount_3j_WinM25"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
    for(int it_win=0; it_win<Nwindow; it_win++){
      if(fabs(Mmumu-35.)<MwidthCuts[it_win]){
        FillHist("NCount_3j_WinM35"+Label, MwidthCuts[it_win]+0.001, weight, 0., 20., 200);
      }
    }
  }
}


void Aug2017_TriLepSR::CheckCutflow(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){

  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");
  if(EMuMu){
    if( !(EleLColl.size()==1 && MuLColl.size()==2) ) return;
    if( !(EleTColl.size()==1 && MuTColl.size()==2) ) return;
    if( !(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) ) return;
    if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10) ) return;
    if( fabs(EleTColl.at(0).Eta())>2.5 ) return;
      if(k_sample_name.Contains("DY") || k_sample_name.Contains("TT_powheg")){
        std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
        int LeptonType=0;
        LeptonType=GetLeptonType(MuTColl.at(1),truthColl);
        if(LeptonType>=-4 && LeptonType<0) goto flowstart_1e2mu;
        LeptonType=GetLeptonType(EleTColl.at(0),truthColl);
        if(LeptonType>=-4 && LeptonType<0) goto flowstart_1e2mu;
        LeptonType=GetLeptonType(MuTColl.at(0),truthColl);
        if(LeptonType>=-4 && LeptonType<0) goto flowstart_1e2mu;
        return;
      }
      flowstart_1e2mu:
      float CurrentIdx=0.;
      float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
    if(Mmumu<12) return;
      FillHist("CutFlow_1e2mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//3lep within simulation reg. & Triggered.
    if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && MuTColl.at(1).Pt()>10) ) return;
      FillHist("CutFlow_1e2mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//Analysis PT Cut
    if(fabs(Mmumu-91.2)<10) return;
      FillHist("CutFlow_1e2mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//OffZ
    if(BJetColl.size()==0) return;
      FillHist("CutFlow_1e2mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//HasB
    if(JetColl.size()<2) return;
      FillHist("CutFlow_1e2mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//geq2j
  }
  if(TriMu){
    if( !(MuLColl.size()==3 && EleLColl.size()==0) ) return;
    if( !(MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 && MuTColl.at(2).Pt()>5) ) return;
      if(k_sample_name.Contains("DY") || k_sample_name.Contains("TT_powheg")){
        std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
        int LeptonType=0;
        LeptonType=GetLeptonType(MuTColl.at(2),truthColl);
        if(LeptonType>=-4 && LeptonType<0) goto flowstart_3mu;
        LeptonType=GetLeptonType(MuTColl.at(1),truthColl);
        if(LeptonType>=-4 && LeptonType<0) goto flowstart_3mu;
        LeptonType=GetLeptonType(MuTColl.at(0),truthColl);
        if(LeptonType>=-4 && LeptonType<0) goto flowstart_3mu;
        return;
      }
      flowstart_3mu:
      float CurrentIdx=0.;
      int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
      int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
      int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
      float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
      float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    if( !(MOSSS1>12 && MOSSS2>12) ) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//3lep within simulation reg. & Triggered.
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 && MuTColl.at(2).Pt()>10) ) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//3lep within simulation reg. & Triggered.
    if( fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10 ) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//3lep within simulation reg. & Triggered.
    if(BJetColl.size()==0) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//HasB
    if(JetColl.size()<2) return;
      FillHist("CutFlow_3mu"+Label, CurrentIdx, weight, 0., 10., 10); CurrentIdx++;//geq2j
  }

}


void Aug2017_TriLepSR::CheckSRYield(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");
//  if(k_sample_name.Contains("TTG")){
//    std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
//    if(EleTColl.size()==0) return;
//    int LepType=GetLeptonType(EleTColl.at(0),truthColl);
//    if(LepType>-5) return;
//  }

  if(EMuMu){
    if( !(EleLColl.size()==1 && MuLColl.size()==2) ) return;
    if( !(EleTColl.size()==1 && MuTColl.size()==2) ) return;
    if( !(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) ) return;
    if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && MuTColl.at(1).Pt()>10) ) return;
    if( fabs(EleTColl.at(0).Eta())>2.5 ) return;

    float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
    float MTW = sqrt(2)*sqrt(MET*EleTColl.at(0).Pt()-METx*EleTColl.at(0).Px()-METy*EleTColl.at(0).Py());
    float M3l=(EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M();
    if(Mmumu<12) return;

    if(fabs(Mmumu-91.2)<10) return;

    float MAWidth=2.5;
    if(JetColl.size()>1){
      FillHist("Nb_2jCutPreB_1e2mu"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("Count_2jCutPreB_1e2mu_PreSel"+Label, 0., weight, 0., 1., 1);
      if( Mmumu>12 && Mmumu<80   ) FillHist("Count_2jCutPreB_1e2mu_M12to80"+Label, 0., weight, 0., 1., 1);
      if( fabs(Mmumu-15)<MAWidth ) FillHist("Count_2jCutPreB_1e2mu_M15"+Label, 0., weight, 0., 1., 1);
      if( fabs(Mmumu-20)<MAWidth ) FillHist("Count_2jCutPreB_1e2mu_M20"+Label, 0., weight, 0., 1., 1);
      if( fabs(Mmumu-25)<MAWidth ) FillHist("Count_2jCutPreB_1e2mu_M25"+Label, 0., weight, 0., 1., 1);
      if( fabs(Mmumu-30)<MAWidth ) FillHist("Count_2jCutPreB_1e2mu_M30"+Label, 0., weight, 0., 1., 1);
      if( fabs(Mmumu-35)<MAWidth ) FillHist("Count_2jCutPreB_1e2mu_M35"+Label, 0., weight, 0., 1., 1);
    }

    if(BJetColl.size()==0) return;

    FillHist("Nj_PreSel_1e2mu"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Mmumu_PreSel_1e2mu"+Label, Mmumu, weight, 0., 200., 40);
    //c.f. Regarding detector resolution not considering FSR tail
    //σ(m, pT) = 0.026 + 0.0065m for the barrel,
    //σ(m, pT) = 0.026 + 0.013m GeV/c2 for the endcap.
    if(JetColl.size()>1){
      FillHist("Mmumu_2jCut_1e2mu"+Label, Mmumu, weight, 0., 200., 40);
      FillHist("Count_2jCut_1e2mu_PreSel"+Label, 0., weight, 0., 1., 1);
      if( Mmumu>12 && Mmumu<80   ) FillHist("Count_2jCut_1e2mu_M12to80"+Label, 0., weight, 0., 1., 1);
      if( fabs(Mmumu-15)<MAWidth ) FillHist("Count_2jCut_1e2mu_M15"+Label, 0., weight, 0., 1., 1);
      if( fabs(Mmumu-20)<MAWidth ) FillHist("Count_2jCut_1e2mu_M20"+Label, 0., weight, 0., 1., 1);
      if( fabs(Mmumu-25)<MAWidth ) FillHist("Count_2jCut_1e2mu_M25"+Label, 0., weight, 0., 1., 1);
      if( fabs(Mmumu-30)<MAWidth ) FillHist("Count_2jCut_1e2mu_M30"+Label, 0., weight, 0., 1., 1);
      if( fabs(Mmumu-35)<MAWidth ) FillHist("Count_2jCut_1e2mu_M35"+Label, 0., weight, 0., 1., 1);
    }
  }
  if(TriMu){
    if( !(MuLColl.size()==3 && EleLColl.size()==0) ) return;
    if( !(MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 && MuTColl.at(2).Pt()>10) ) return;

    int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    float MAWidth=2.5;
    if( !(MOSSS1>12 && MOSSS2>12) ) return;
    if( fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10 ) return;

    if(JetColl.size()>1){
      FillHist("Nb_2jCutPreB_3mu"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("Count_2jCutPreB_MOSSS2_3mu_PreSel"+Label, 0., weight, 0., 1., 1);
      if(MOSSS2<80)               FillHist("Count_2jCutPreB_MOSSS2_3mu_M12to80"+Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-15)<MAWidth) FillHist("Count_2jCutPreB_MOSSS2_3mu_M15"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-20)<MAWidth) FillHist("Count_2jCutPreB_MOSSS2_3mu_M20"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-25)<MAWidth) FillHist("Count_2jCutPreB_MOSSS2_3mu_M25"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-30)<MAWidth) FillHist("Count_2jCutPreB_MOSSS2_3mu_M30"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-35)<MAWidth) FillHist("Count_2jCutPreB_MOSSS2_3mu_M35"    +Label, 0., weight, 0., 1., 1);
    }

    if( BJetColl.size()==0 ) return;

    if(JetColl.size()>1){
      FillHist("MOSSS1_2j"+Label, MOSSS1, weight, 0., 200., 40);
      FillHist("MOSSS2_2j"+Label, MOSSS2, weight, 0., 200., 40);
      FillHist("Count_2jCut_MOSSS2_3mu_PreSel"+Label, 0., weight, 0., 1., 1);
      if(MOSSS2<80)               FillHist("Count_2jCut_MOSSS2_3mu_M12to80"+Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-15)<MAWidth) FillHist("Count_2jCut_MOSSS2_3mu_M15"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-20)<MAWidth) FillHist("Count_2jCut_MOSSS2_3mu_M20"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-25)<MAWidth) FillHist("Count_2jCut_MOSSS2_3mu_M25"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-30)<MAWidth) FillHist("Count_2jCut_MOSSS2_3mu_M30"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-35)<MAWidth) FillHist("Count_2jCut_MOSSS2_3mu_M35"    +Label, 0., weight, 0., 1., 1);
    }
    if(JetColl.size()>2){
      FillHist("MOSSS1_3j"+Label, MOSSS1, weight, 0., 200., 40);
      FillHist("MOSSS2_3j"+Label, MOSSS2, weight, 0., 200., 40);
      FillHist("Count_3jCut_MOSSS2_3mu_PreSel"+Label, 0., weight, 0., 1., 1);
      if(MOSSS2<80)                 FillHist("Count_3jCut_MOSSS2_3mu_M12to80"+Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-15)<MAWidth  ) FillHist("Count_3jCut_MOSSS2_3mu_M15"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-20)<MAWidth  ) FillHist("Count_3jCut_MOSSS2_3mu_M20"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-25)<MAWidth  ) FillHist("Count_3jCut_MOSSS2_3mu_M25"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-30)<MAWidth  ) FillHist("Count_3jCut_MOSSS2_3mu_M30"    +Label, 0., weight, 0., 1., 1);
      if(fabs(MOSSS2-35)<MAWidth  ) FillHist("Count_3jCut_MOSSS2_3mu_M35"    +Label, 0., weight, 0., 1., 1);
    }
  }

}


void Aug2017_TriLepSR::CheckSigBkgKinematics(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool TriMu=false, EMuMu=false;
  int IdxOS=-1, IdxSS1=-1, IdxSS2=-1;
  float MOSSS1=0., MOSSS2=0., Mmumu=0.;
  std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
  int IdxA1=-1, IdxA2=-1, CountA=0;
  for(int i=2; i<truthColl.size(); i++){
    if(truthColl.at(i).GenStatus()!=1) continue;
    int fpid=fabs(truthColl.at(i).PdgId());
    if( fpid!=11 && fpid!=13 ) continue;
    if( (fpid==11 && fabs(truthColl.at(i).Eta())>2.5) || (fpid==13 && fabs(truthColl.at(i).Eta())>2.4) ) continue;

    int LepType=GetLeptonType(i,truthColl);
    if(LepType==1) FillHist("PTlW_Gen", truthColl.at(i).Pt(), weight, 0., 200., 40);
    if(LepType==2){
      if(CountA==0){ IdxA1=i; CountA++; }
      else{ if(truthColl.at(i).Pt()>truthColl.at(IdxA1).Pt()){ IdxA2=IdxA1; IdxA1=i; CountA++; }
            else{ IdxA2=i; CountA++; }
      }
      FillHist("PTlA_Gen", truthColl.at(i).Pt(), weight, 0., 200., 40);
    }
  }
  if(IdxA1!=-1) FillHist("PTlA1_Gen", truthColl.at(IdxA1).Pt(), weight, 0., 200., 40);
  if(IdxA2!=-1) FillHist("PTlA2_Gen", truthColl.at(IdxA2).Pt(), weight, 0., 200., 40);

  if( MuTColl.size()==3 && MuLColl.size()==3 && EleTColl.size()==0 && EleLColl.size()==0 ) TriMu=true;
  if( MuTColl.size()==2 && MuLColl.size()==2 && EleTColl.size()==1 && EleLColl.size()==1 ) EMuMu=true;
  if( TriMu && fabs(SumCharge(MuTColl))!=1 ) TriMu=false;
  if( EMuMu && SumCharge(MuTColl)!=0       ) EMuMu=false;
  if( TriMu ){
    IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    Mmumu  = MOSSS2;
  }
  if( EMuMu ){ Mmumu = (MuTColl.at(0)+MuTColl.at(1)).M(); }
  if( EMuMu && !(Mmumu>12 && fabs(Mmumu-91.2)>10) ) EMuMu=false;
  if( TriMu && !(MOSSS1>12 && MOSSS2>12 && fabs(MOSSS1-91.2)>10 && fabs(MOSSS2-91.2)>10) ) TriMu=false;
  if( !(TriMu || EMuMu) ) return;

  std::vector<snu::KMuon>     MuAColl, MuWColl,  MuFkColl; 
  std::vector<snu::KElectron>         EleWColl, EleFkColl; 
    for(int i=0; i<(int) MuTColl.size(); i++){
      int LepType=GetLeptonType(MuTColl.at(i), truthColl);
      if     (LepType==2) MuAColl.push_back(MuTColl.at(i));
      else if(LepType==1) MuWColl.push_back(MuTColl.at(i));
      else if(LepType<0 ) MuFkColl.push_back(MuTColl.at(i));
    }
    for(int i=0; i<(int) EleTColl.size(); i++){
      int LepType=GetLeptonType(EleTColl.at(i), truthColl);
      if     (LepType==1) EleWColl.push_back(EleTColl.at(i));
      else if(LepType<0 ) EleFkColl.push_back(EleTColl.at(i));
    }


  for(int i=0; i<(int) EleWColl.size(); i++){
    FillHist("PT_lW", EleWColl.at(i).Pt(), weight, 0., 200., 40);
  }
  for(int i=0; i<(int) MuWColl.size(); i++){
    FillHist("PT_lW", MuWColl.at(i).Pt(), weight, 0., 200., 40);
  }
  for(int i=0; i<(int) MuAColl.size(); i++){
    FillHist("PT_lA", MuAColl.at(i).Pt(), weight, 0., 200., 40);
    if(i==0) FillHist("PT_lA1", MuAColl.at(i).Pt(), weight, 0., 200., 40);
    if(i==1) FillHist("PT_lA2", MuAColl.at(i).Pt(), weight, 0., 200., 40);
  }
  for(int i=0; i<(int) MuFkColl.size(); i++){
    FillHist("PT_MuFk", MuFkColl.at(i).Pt(), weight, 0., 200., 40);
  }
  for(int i=0; i<(int) EleFkColl.size(); i++){
    FillHist("PT_EleFk", EleFkColl.at(i).Pt(), weight, 0., 200., 40);
  }
  if(EMuMu){
    FillHist("PT_Ele_1e2mu", EleTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("PT_Mu1_1e2mu", MuTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("PT_Mu2_1e2mu", MuTColl.at(1).Pt(), weight, 0., 200., 40);
  }
  if(TriMu){
    FillHist("PT_Mu1_3mu", MuTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("PT_Mu2_3mu", MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("PT_Mu3_3mu", MuTColl.at(2).Pt(), weight, 0., 200., 40);
  }

}


void Aug2017_TriLepSR::CheckSRDist(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");
  vector<float> MmumuEdges;
  float CurrentEdge=0.;
   for(int i=0; i<16; i++){
     if     (i==0) CurrentEdge=10.;
     else if(i==1) CurrentEdge=12.5;
     else if(i!=15) CurrentEdge+=5;
     else CurrentEdge=80;
     MmumuEdges.push_back(CurrentEdge);
   }
     
  if(EMuMu){
    if( !(EleLColl.size()==1 && MuLColl.size()==2) ) return;
    if( !(EleTColl.size()==1 && MuTColl.size()==2) ) return;
    if( !(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) ) return;
    if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && MuTColl.at(1).Pt()>10) ) return;
    if( fabs(EleTColl.at(0).Eta())>2.5 ) return;

    float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
    float MTW = sqrt(2)*sqrt(MET*EleTColl.at(0).Pt()-METx*EleTColl.at(0).Px()-METy*EleTColl.at(0).Py());
    float M3l=(EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M();
    if(Mmumu<12) return;

      FillHist("CutFlow_1e2mu"+Label, 0., weight, 0., 20., 20);
    if(fabs(Mmumu-91.2)<10) return;
      FillHist("CutFlow_1e2mu"+Label, 1., weight, 0., 20., 20);
      FillHist("Nb_3lOSOffZ"+Label, BJetColl.size(), weight, 0., 5., 5);
    if(BJetColl.size()==0) return;
      FillHist("CutFlow_1e2mu"+Label, 2., weight, 0., 20., 20);
      FillHist("Nj_3lOSOffZHasB"+Label, JetColl.size(), weight, 0., 10., 10);
    if(JetColl.size()<2) return;
      FillHist("CutFlow_1e2mu"+Label, 3., weight, 0., 20., 20);

    int Idxj1_W=-1, Idxj2_W=-1;
    snu::KParticle v; v.SetPxPyPzE(METx, METy, 0, sqrt(METx*METx+METy*METy));
    float Pzv1=GetvPz(v, EleTColl.at(0), 1), Pzv2=GetvPz(v, EleTColl.at(0), 2);
    float Pzv=fabs(Pzv1)<fabs(Pzv2)? Pzv1:Pzv2;
    v.SetXYZM(METx,METy,Pzv,0.);

    float M3lv=(EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)+v).M(), M2l2j=0.; 
    Idxj1_W=GetDaughterCandIdx(JetColl, "W", -1., "Lead");
    Idxj2_W=GetDaughterCandIdx(JetColl, "W", -1., "Subl");
    M2l2j=(MuTColl.at(0)+MuTColl.at(1)+JetColl.at(Idxj1_W)+JetColl.at(Idxj2_W)).M();

    FillHist("Nj_SR2j"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_SR2j"+Label, BJetColl.size(), weight, 0., 5., 5);
    FillHist("PTe_SR2j"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
    FillHist("PTmu1_SR2j"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
    FillHist("PTmu2_SR2j"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("Mmumu_CntBin_SR2j"+Label, Mmumu, weight, &MmumuEdges[0], MmumuEdges.size()-1);
    FillHist("Mmumu_SR2j"+Label, Mmumu, weight, 0., 200., 200);
    FillHist("M2l2j_SR2j"+Label, M2l2j, weight, 0., 1000., 100);
    FillHist("M3lv_SR2j"+Label, M3lv, weight, 0., 1000., 100);
    if(Mmumu<80){
      FillHist("M2l2j_lt80_SR2j"+Label, M2l2j, weight, 0., 1000., 100);
      FillHist("M3lv_lt80_SR2j"+Label, M3lv, weight, 0., 1000., 100);
    }
  }
  if(TriMu){
    if( !(MuLColl.size()==3 && EleLColl.size()==0) ) return;
    if( !(MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 && MuTColl.at(2).Pt()>10) ) return;

    int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    float Mmumu  = MOSSS2;
    if( !(MOSSS1>12 && MOSSS2>12) ) return;

      FillHist("CutFlow_3mu"+Label, 0., weight, 0., 20., 20);
    if( fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10 ) return;
      FillHist("CutFlow_3mu"+Label, 1., weight, 0., 20., 20);
      FillHist("Nb_3lOSOffZ"+Label, BJetColl.size(), weight, 0., 5., 5);
    if(BJetColl.size()==0) return;
      FillHist("CutFlow_3mu"+Label, 2., weight, 0., 20., 20);
      FillHist("Nj_3lOSOffZHasB"+Label, JetColl.size(), weight, 0., 10., 10);
    if(JetColl.size()<2) return;
      FillHist("CutFlow_3mu"+Label, 3., weight, 0., 20., 20);


    int Idxj1_W=-1, Idxj2_W=-1;
    snu::KParticle v; v.SetPxPyPzE(METx, METy, 0, sqrt(METx*METx+METy*METy));
    float Pzv1=GetvPz(v, MuTColl.at(0), 1), Pzv2=GetvPz(v, MuTColl.at(0), 2);
    float Pzv=fabs(Pzv1)<fabs(Pzv2)? Pzv1:Pzv2;
    v.SetXYZM(METx,METy,Pzv,0.);

    float M3lv=(MuTColl.at(0)+MuTColl.at(1)+MuTColl.at(2)+v).M(), M2l2j=0.; 
    Idxj1_W=GetDaughterCandIdx(JetColl, "W", -1., "Lead");
    Idxj2_W=GetDaughterCandIdx(JetColl, "W", -1., "Subl");
    M2l2j=(MuTColl.at(1)+MuTColl.at(2)+JetColl.at(Idxj1_W)+JetColl.at(Idxj2_W)).M();

    FillHist("Nj_SR2j"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_SR2j"+Label, BJetColl.size(), weight, 0., 5., 5);
    FillHist("PTmu1_SR2j"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
    FillHist("PTmu2_SR2j"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("PTmu3_SR2j"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
    FillHist("MOSSS1_SR2j"+Label, MOSSS1, weight, 0., 200., 80);
    FillHist("Mmumu_CntBin_SR2j"+Label, Mmumu, weight, &MmumuEdges[0], MmumuEdges.size()-1);
    FillHist("Mmumu_SR2j"+Label, Mmumu, weight, 0., 200., 200);
    FillHist("M2l2j_SR2j"+Label, M2l2j, weight, 0., 1000., 100);
    FillHist("M3lv_SR2j"+Label, M3lv, weight, 0., 1000., 100);
    if(MOSSS2<80){
      FillHist("M2l2j_lt80_SR2j"+Label, M2l2j, weight, 0., 1000., 100);
      FillHist("M3lv_lt80_SR2j"+Label, M3lv, weight, 0., 1000., 100);
    }
  }
}


void Aug2017_TriLepSR::DoSystRun(TString Cycle, TString Mode, std::vector<snu::KElectron> EleColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight){

  int SystDir=0; TString SystKindLabel="", SystDirLabel="", ChannelLabel="";
  if(Mode.Contains("Syst")){
    if(isData&& !k_running_nonprompt) return;
    if     (Mode.Contains("Up"))      {SystDir= 1; SystDirLabel="_systup";}
    else if(Mode.Contains("Down"))    {SystDir=-1; SystDirLabel="_systdown";}
    if     (Mode.Contains("PU"))      SystKindLabel="_PU";
    else if(Mode.Contains("Nvtx"))    SystKindLabel="_Nvtx";
    else if(Mode.Contains("JES"))     SystKindLabel="_JES";
    else if(Mode.Contains("JER"))     SystKindLabel="_JER";
    else if(Mode.Contains("Uncl"))    SystKindLabel="_Uncl";
    else if(Mode.Contains("ElID"))    SystKindLabel="_ElID";
    else if(Mode.Contains("MuID"))    SystKindLabel="_MuID";
    else if(Mode.Contains("ElEn"))    SystKindLabel="_ElEn";
    else if(Mode.Contains("MuEn"))    SystKindLabel="_MuEn";
    else if(Mode.Contains("FR"))      SystKindLabel="_FR";
    else if(Mode.Contains("BTag_L"))  SystKindLabel="_BTag_L";
    else if(Mode.Contains("BTag_BC")) SystKindLabel="_BTag_BC";
    else if(Mode.Contains("Trig"))    SystKindLabel="_Trig";
    else if(Mode.Contains("Xsec"))    SystKindLabel="_Xsec";
    else if(Mode.Contains("Conv"))    SystKindLabel="_Conv";
  }
  if     (Mode.Contains("EMuMu"))   ChannelLabel="EMuMu";
  else if(Mode.Contains("TriMu"))   ChannelLabel="TriMu";


  if     (Cycle=="SRYield"){
    CheckSRYield(MuColl, MuLColl, EleColl, EleLColl, JetColl,  BJetColl, MET, METx, METy, weight, SystDirLabel+SystKindLabel, ChannelLabel);
  }//End of Closure Cycle
  else if(Cycle=="SRDist"){
    CheckSRDist(MuColl, MuLColl, EleColl, EleLColl, JetColl, BJetColl, MET, METx, METy, weight, SystDirLabel+SystKindLabel, ChannelLabel);
  }


  return;
}



float Aug2017_TriLepSR::ConeCorrectedPT(snu::KMuon Mu, float TightIsoCut){

  //float PTCorr=Mu.Pt();
  //float PTCorr=Mu.Pt()*(1.+RochIso(Mu,"0.4"));
  float PTCorr=Mu.Pt()*(1.+max((float) 0.,RochIso(Mu,"0.4")-TightIsoCut));
  return PTCorr;

}



float Aug2017_TriLepSR::ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut){

  //float PTCorr=Ele.Pt()*(1+Ele.PFRelIso(0.3));
  float PTCorr=Ele.Pt()*(1.+max(0.,Ele.PFRelIso(0.3)-TightIsoCut));
  return PTCorr;

}


float Aug2017_TriLepSR::FakeRateData(snu::KMuon Mu, TString Option){

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
  if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.276317;
      else if(PTCorr<20)  FR=0.191793;
      else if(PTCorr<25)  FR=0.167908;
      else if(PTCorr<35)  FR=0.162617;
      else                FR=0.160329;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.30778 ;
      else if(PTCorr<20)  FR=0.220861;
      else if(PTCorr<25)  FR=0.216049;
      else if(PTCorr<35)  FR=0.186631;
      else                FR=0.188504;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.339629;
      else if(PTCorr<20)  FR=0.262338;
      else if(PTCorr<25)  FR=0.252174;
      else if(PTCorr<35)  FR=0.220018;
      else                FR=0.212394;
    }
    else{
      if     (PTCorr<15)  FR=0.339974;
      else if(PTCorr<20)  FR=0.271258;
      else if(PTCorr<25)  FR=0.249348;
      else if(PTCorr<35)  FR=0.25788 ;
      else                FR=0.240026;
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
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.103166 ;
      else if(PTCorr<20)  FR=0.0800024;
      else if(PTCorr<25)  FR=0.0847176;
      else if(PTCorr<35)  FR=0.0680578;
      else                FR=0.0658314;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.095482 ;
      else if(PTCorr<20)  FR=0.0713299;
      else if(PTCorr<25)  FR=0.102061 ;
      else if(PTCorr<35)  FR=0.114623 ;
      else                FR=0.0850804;
    }
    else{
      if     (PTCorr<15)  FR=0.0603521;
      else if(PTCorr<20)  FR=0.0918451;
      else if(PTCorr<25)  FR=0.0671922;
      else if(PTCorr<35)  FR=0.128504 ;
      else                FR=0.0964251;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.282211;
      else if(PTCorr<20)  FR=0.2111  ;
      else if(PTCorr<25)  FR=0.187283;
      else if(PTCorr<35)  FR=0.195845;
      else                FR=0.193727;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.321604;
      else if(PTCorr<20)  FR=0.255115;
      else if(PTCorr<25)  FR=0.253692;
      else if(PTCorr<35)  FR=0.197637;
      else                FR=0.238313;
    }
    else{
      if     (PTCorr<15)  FR=0.35875 ;
      else if(PTCorr<20)  FR=0.303843;
      else if(PTCorr<25)  FR=0.298489;
      else if(PTCorr<35)  FR=0.255767;
      else                FR=0.266638;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_NoReq_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.276247;
      else if(PTCorr<20)  FR=0.191163;
      else if(PTCorr<25)  FR=0.158219;
      else if(PTCorr<35)  FR=0.155122;
      else                FR=0.153071;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.309303;
      else if(PTCorr<20)  FR=0.219685;
      else if(PTCorr<25)  FR=0.211428;
      else if(PTCorr<35)  FR=0.170437;
      else                FR=0.19083 ;
    }
    else{
      if     (PTCorr<15)  FR=0.339582;
      else if(PTCorr<20)  FR=0.263034;
      else if(PTCorr<25)  FR=0.244504;
      else if(PTCorr<35)  FR=0.222135;
      else                FR=0.219002;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1NoChi_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.104934 ;
      else if(PTCorr<20)  FR=0.0810565;
      else if(PTCorr<25)  FR=0.085597 ;
      else if(PTCorr<35)  FR=0.0684988;
      else                FR=0.0665389;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.0977776;
      else if(PTCorr<20)  FR=0.0729936;
      else if(PTCorr<25)  FR=0.10299  ;
      else if(PTCorr<35)  FR=0.116561 ;
      else                FR=0.0859202;
    }
    else{
      if     (PTCorr<15)  FR=0.0610895;
      else if(PTCorr<20)  FR=0.0923621;
      else if(PTCorr<25)  FR=0.0671893;
      else if(PTCorr<35)  FR=0.128802 ;
      else                FR=0.0967269;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1NoChi_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.308976;
      else if(PTCorr<20)  FR=0.229748;
      else if(PTCorr<25)  FR=0.202024;
      else if(PTCorr<35)  FR=0.211457;
      else                FR=0.215925;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.351577;
      else if(PTCorr<20)  FR=0.275117;
      else if(PTCorr<25)  FR=0.266206;
      else if(PTCorr<35)  FR=0.207133;
      else                FR=0.252437;
    }
    else{
      if     (PTCorr<15)  FR=0.364771;
      else if(PTCorr<20)  FR=0.309155;
      else if(PTCorr<25)  FR=0.302773;
      else if(PTCorr<35)  FR=0.265689;
      else                FR=0.275625;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1NoChi_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.0730267;
      else if(PTCorr<20)  FR=0.0422678;
      else if(PTCorr<25)  FR=0.0394787;
      else if(PTCorr<35)  FR=0.0292311;
      else                FR=0.0247411;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.075509 ;
      else if(PTCorr<20)  FR=0.0421607;
      else if(PTCorr<25)  FR=0.0535583;
      else if(PTCorr<35)  FR=0.0579744;
      else                FR=0.0374653;
    }
    else{
      if     (PTCorr<15)  FR=0.0488572;
      else if(PTCorr<20)  FR=0.0598638;
      else if(PTCorr<25)  FR=0.036773 ;
      else if(PTCorr<35)  FR=0.0710471;
      else                FR=0.0479654;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1NoChi_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.243067 ;
      else if(PTCorr<20)  FR=0.120063 ;
      else if(PTCorr<25)  FR=0.100958 ;
      else if(PTCorr<35)  FR=0.0993144;
      else                FR=0.0935684;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.288401;
      else if(PTCorr<20)  FR=0.159155;
      else if(PTCorr<25)  FR=0.150465;
      else if(PTCorr<35)  FR=0.106834;
      else                FR=0.122566;
    }
    else{
      if     (PTCorr<15)  FR=0.308748;
      else if(PTCorr<20)  FR=0.195025;
      else if(PTCorr<25)  FR=0.182063;
      else if(PTCorr<35)  FR=0.153384;
      else                FR=0.148075;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1sig4NoChi_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.333358;
      else if(PTCorr<20)  FR=0.236995;
      else if(PTCorr<25)  FR=0.23349 ;
      else if(PTCorr<35)  FR=0.17589 ;
      else                FR=0.168586;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.323098;
      else if(PTCorr<20)  FR=0.207417;
      else if(PTCorr<25)  FR=0.261432;
      else if(PTCorr<35)  FR=0.288036;
      else                FR=0.21438 ;
    }
    else{
      if     (PTCorr<15)  FR=0.191651;
      else if(PTCorr<20)  FR=0.237084;
      else if(PTCorr<25)  FR=0.170018;
      else if(PTCorr<35)  FR=0.283506;
      else                FR=0.221901;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1sig4NoChi_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.434692;
      else if(PTCorr<20)  FR=0.303991;
      else if(PTCorr<25)  FR=0.246176;
      else if(PTCorr<35)  FR=0.265866;
      else                FR=0.267146;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.472631;
      else if(PTCorr<20)  FR=0.360498;
      else if(PTCorr<25)  FR=0.33543 ;
      else if(PTCorr<35)  FR=0.252456;
      else                FR=0.315846;
    }
    else{
      if     (PTCorr<15)  FR=0.463157;
      else if(PTCorr<20)  FR=0.373889;
      else if(PTCorr<25)  FR=0.363318;
      else if(PTCorr<35)  FR=0.321106;
      else                FR=0.339224;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4NoChi_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.235023 ;
      else if(PTCorr<20)  FR=0.126541 ;
      else if(PTCorr<25)  FR=0.108612 ;
      else if(PTCorr<35)  FR=0.0807077;
      else                FR=0.0657107;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.23597  ;
      else if(PTCorr<20)  FR=0.121365 ;
      else if(PTCorr<25)  FR=0.134578 ;
      else if(PTCorr<35)  FR=0.146404 ;
      else                FR=0.0986332;
    }
    else{
      if     (PTCorr<15)  FR=0.149694 ;
      else if(PTCorr<20)  FR=0.150774 ;
      else if(PTCorr<25)  FR=0.0950449;
      else if(PTCorr<35)  FR=0.158923 ;
      else                FR=0.114237 ;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4NoChi_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.344249;
      else if(PTCorr<20)  FR=0.163312;
      else if(PTCorr<25)  FR=0.124449;
      else if(PTCorr<35)  FR=0.125493;
      else                FR=0.116852;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.388531;
      else if(PTCorr<20)  FR=0.208789;
      else if(PTCorr<25)  FR=0.188631;
      else if(PTCorr<35)  FR=0.131483;
      else                FR=0.153284;
    }
    else{
      if     (PTCorr<15)  FR=0.391918;
      else if(PTCorr<20)  FR=0.235671;
      else if(PTCorr<25)  FR=0.216606;
      else if(PTCorr<35)  FR=0.182999;
      else                FR=0.181455;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.235023 ;
      else if(PTCorr<20)  FR=0.126737 ;
      else if(PTCorr<25)  FR=0.10877  ;
      else if(PTCorr<35)  FR=0.0810774;
      else                FR=0.0657765;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.23597  ;
      else if(PTCorr<20)  FR=0.121983 ;
      else if(PTCorr<25)  FR=0.135113 ;
      else if(PTCorr<35)  FR=0.146711 ;
      else                FR=0.0987709;
    }
    else{
      if     (PTCorr<15)  FR=0.150243 ;
      else if(PTCorr<20)  FR=0.151052 ;
      else if(PTCorr<25)  FR=0.0953203;
      else if(PTCorr<35)  FR=0.159848 ;
      else                FR=0.114909 ;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.345149;
      else if(PTCorr<20)  FR=0.164104;
      else if(PTCorr<25)  FR=0.124837;
      else if(PTCorr<35)  FR=0.126013;
      else                FR=0.118599;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.389373;
      else if(PTCorr<20)  FR=0.208879;
      else if(PTCorr<25)  FR=0.189129;
      else if(PTCorr<35)  FR=0.131404;
      else                FR=0.155852;
    }
    else{
      if     (PTCorr<15)  FR=0.39351 ;
      else if(PTCorr<20)  FR=0.237086;
      else if(PTCorr<25)  FR=0.217314;
      else if(PTCorr<35)  FR=0.183934;
      else                FR=0.185743;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.276577 ;
      else if(PTCorr<20)  FR=0.142112 ;
      else if(PTCorr<25)  FR=0.117096 ;
      else if(PTCorr<35)  FR=0.0858086;
      else                FR=0.0684839;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.351973;
      else if(PTCorr<20)  FR=0.154119;
      else if(PTCorr<25)  FR=0.154107;
      else if(PTCorr<35)  FR=0.162744;
      else                FR=0.109306;
    }
    else{
      if     (PTCorr<15)  FR=0.325502;
      else if(PTCorr<20)  FR=0.252336;
      else if(PTCorr<25)  FR=0.144032;
      else if(PTCorr<35)  FR=0.212277;
      else                FR=0.151429;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.363658;
      else if(PTCorr<20)  FR=0.169303;
      else if(PTCorr<25)  FR=0.127615;
      else if(PTCorr<35)  FR=0.128469;
      else                FR=0.121185;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.423617;
      else if(PTCorr<20)  FR=0.219131;
      else if(PTCorr<25)  FR=0.196251;
      else if(PTCorr<35)  FR=0.134896;
      else                FR=0.160779;
    }
    else{
      if     (PTCorr<15)  FR=0.475586;
      else if(PTCorr<20)  FR=0.269348;
      else if(PTCorr<25)  FR=0.242709;
      else if(PTCorr<35)  FR=0.19999 ;
      else                FR=0.201519;
    }
  }



  return FR;
}



float Aug2017_TriLepSR::FakeRateData(snu::KElectron Ele, TString Option){

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

  if(Option=="QCD_POGWP90Isop06IPp025p1sig4_HctoWAFakeLoose_ConeSUSY"){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.106151 ;
      else if(PTCorr<50)  FR=0.0663149;
      else                FR=0.0648656;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.1088   ;
      else if(PTCorr<50)  FR=0.0704312;
      else                FR=0.0653196;
    }
    else{
      if     (PTCorr<35)  FR=0.143974;
      else if(PTCorr<50)  FR=0.104829;
      else                FR=0.126112;
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
        if     (PTCorr<35)   FR*=(1.- 0.077764 );
        else if(PTCorr<50)   FR*=(1.- 0.400603 );
        else                 FR*=(1.- 0.51524  );
      }
      else if(fEta<1.479){
        if     (PTCorr<35)   FR*=(1.- 0.169421 );
        else if(PTCorr<50)   FR*=(1.- 0.424976 );
        else                 FR*=(1.- 0.256287 );
      }
      else{
        if     (PTCorr<35)   FR*=(1.- 0.0561248 );
        else if(PTCorr<50)   FR*=(1.- 0.103     );
        else                 FR*=(1.- 0.111331  );
      }
    }
    else if(Syst_FR && SystDir>0){
      if(fEta<0.8){
        if     (PTCorr<35)   FR*=(1.+ 0.109443 );
        else if(PTCorr<50)   FR*=(1.+ 0.202818 );
        else                 FR*=(1.+ 0.374888 );
      }
      else if(fEta<1.479){
        if     (PTCorr<35)   FR*=(1.+ 0.0973585 );
        else if(PTCorr<50)   FR*=(1.+ 0.148029  );
        else                 FR*=(1.+ 0.246644  );
      }
      else{
        if     (PTCorr<35)   FR*=(1.+ 0.161426  );
        else if(PTCorr<50)   FR*=(1.+ 0.0707937 );
        else                 FR*=(1.+ 0.0566853 );
      }
    }
  }

  return FR;
} 


//float Aug2017_TriLepSR::OptimiseSelection(std::vector<snu::KMuon> MuPreColl, std::vector<snu::KElectron> ElePreColl, std::vector<snu::KJet> JetPreNoVetoColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, TString Label){
//
//  float Mmumu=
//  if(
//
//}

float Aug2017_TriLepSR::GetFakeWeight(std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleLColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, TString Option){

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


  for(int i=0; i<(int) MuLColl.size(); i++){
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
  for(int i=0; i<(int) EleLColl.size(); i++){
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


float Aug2017_TriLepSR::GetFakeWeight(std::vector<snu::KMuon>& MuLColl, std::vector<snu::KElectron>& EleLColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, vector<snu::KJet>& BJetNoVetoColl, TString Option){

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


bool Aug2017_TriLepSR::IsNearBJet(snu::KMuon Mu, std::vector<snu::KJet>& bjetNoVetoColl){

  bool IsNearB=false;
  for(int i=0; i<(int) bjetNoVetoColl.size(); i++){
    if(Mu.DeltaR(bjetNoVetoColl.at(i))<0.4){
      IsNearB=true; break;
    }
  }

  return IsNearB;
}



float Aug2017_TriLepSR::GetPreTrigPURW(int Nvtx){

   float weight=1.;
   if(Nvtx<2) weight=1.22087;
   else if(Nvtx<3) weight=4.79564;
   else if(Nvtx<4) weight=3.84689;
   else if(Nvtx<5) weight=2.47179;
   else if(Nvtx<6) weight=3.74121;
   else if(Nvtx<7) weight=3.16705;
   else if(Nvtx<8) weight=3.04533;
   else if(Nvtx<9) weight=2.85291;
   else if(Nvtx<10) weight=2.50907;
   else if(Nvtx<11) weight=2.29848;
   else if(Nvtx<12) weight=2.25329;
   else if(Nvtx<13) weight=1.9864;
   else if(Nvtx<14) weight=1.73597;
   else if(Nvtx<15) weight=1.57862;
   else if(Nvtx<16) weight=1.38849;
   else if(Nvtx<17) weight=1.21789;
   else if(Nvtx<18) weight=1.05361;
   else if(Nvtx<19) weight=0.909483;
   else if(Nvtx<20) weight=0.787018;
   else if(Nvtx<21) weight=0.699577;
   else if(Nvtx<22) weight=0.606739;
   else if(Nvtx<23) weight=0.533779;
   else if(Nvtx<24) weight=0.47826;
   else if(Nvtx<25) weight=0.430963;
   else if(Nvtx<26) weight=0.38291;
   else if(Nvtx<27) weight=0.368938;
   else if(Nvtx<28) weight=0.306386;
   else if(Nvtx<29) weight=0.321641;
   else if(Nvtx<30) weight=0.267767;
   else if(Nvtx<31) weight=0.303042;
   else if(Nvtx<32) weight=0.220255;
   else if(Nvtx<33) weight=0.237769;
   else if(Nvtx<34) weight=0.214341;
   else if(Nvtx<35) weight=0.254868;
   else if(Nvtx<36) weight=0.318919;
   else if(Nvtx<37) weight=0.290442;
   else if(Nvtx<38) weight=0.204838;
   else if(Nvtx<39) weight=0.389061;
   else if(Nvtx<40) weight=0.309147;
   else if(Nvtx<41) weight=0.286387;
   else if(Nvtx<42) weight=0.266532;
   else if(Nvtx<43) weight=0.327612;
   else if(Nvtx<44) weight=0.575947;
   else if(Nvtx<45) weight=0.23432;
   else if(Nvtx<46) weight=0.572198;
   else if(Nvtx<47) weight=0.541966;
   else if(Nvtx<48) weight=0.384433;
   else if(Nvtx<49) weight=1.57278;
   else weight=0.389116;

   return weight;
}



void Aug2017_TriLepSR::FillCutFlow(TString cut, float weight){
  
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



void Aug2017_TriLepSR::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Aug2017_TriLepSR::MakeHistograms(){
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
  // *  Remove//Overide this Aug2017_TriLepSRCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Aug2017_TriLepSR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
