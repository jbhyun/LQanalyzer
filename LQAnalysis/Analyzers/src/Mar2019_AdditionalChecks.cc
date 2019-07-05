// $Id: Mar2019_AdditionalChecks.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQMar2019_AdditionalChecks Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Mar2019_AdditionalChecks.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Mar2019_AdditionalChecks);

 Mar2019_AdditionalChecks::Mar2019_AdditionalChecks() : AnalyzerCore(), out_muons(0) {

   SetLogName("Mar2019_AdditionalChecks");
   Message("In Mar2019_AdditionalChecks constructor", INFO);
   InitialiseAnalysis();
 }


 void Mar2019_AdditionalChecks::InitialiseAnalysis() throw( LQError ) {
   
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

void Mar2019_AdditionalChecks::ExecuteEvents()throw( LQError ){

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
      geneff_weight            = GenFilterEfficiency(k_sample_name);//For old samples
      gennorm_weight           = SignalNorm(k_sample_name, 20.);    //For old samples
      FillHist("Weight_PU", pileup_reweight, 1., 0., 20., 200);
   }
   weight *= k_factor_weight*geneff_weight*gennorm_weight*pileup_reweight;

 
   //Total Event(MCweight+PUrw)////////////////
   FillCutFlow("NoCut", weight*pileup_reweight);


   bool EMuMu=false, TriMu=false;
   bool SRYield=false, SRDist=false, SigEffValid=false, JetVetoCheck=false, FkSrcCheck=false;
   bool SystRun=false, TestRun=false;
   TString Cycle="";
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("EMuMu"))         EMuMu         = true;
     else if(k_flags.at(i).Contains("TriMu"))         TriMu         = true;
     else if(k_flags.at(i).Contains("SystRun"))       SystRun       = true;
     else if(k_flags.at(i).Contains("SRYield"))      {SRYield       = true; Cycle="SRYield";}
     else if(k_flags.at(i).Contains("SRDist"))       {SRDist        = true; Cycle="SRDist";}
     else if(k_flags.at(i).Contains("TestRun"))       TestRun       = true;
     else if(k_flags.at(i).Contains("SigEffValid"))   SigEffValid   = true;
     else if(k_flags.at(i).Contains("JetVetoCheck"))  JetVetoCheck  = true;
     else if(k_flags.at(i).Contains("FkSrcCheck"))    FkSrcCheck    = true;
   }

   bool Dec_lvqq=false, Dec_llvv=false, Dec_else=false;
   bool Dec_evqq=false, Dec_mvqq=false, Dec_tvqq=false;
   bool Dec_evqq_hard = false, Dec_mvqq_hard = false;
   bool Dec_eevv=false, Dec_emvv=false, Dec_etvv=false, Dec_mmvv=false, Dec_mtvv=false, Dec_ttvv=false;

   if(SigEffValid){

     std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
//       PrintTruth();

     int Nb_t=0, Nq_W=0, Nl_W=0, Ne_W=0, Nm_W=0, Nt_W=0, Nv_W=0, Nm_A=0;
     int Nb_t_hard=0, Nq_W_hard=0, Nl_W_hard=0, Ne_W_hard=0, Nm_W_hard=0, Nm_A_hard=0;

     for(int it=2; it<truthColl.size(); it++){
       int   apid   = abs(truthColl.at(it).PdgId());
       int   GenSt  = truthColl.at(it).GenStatus();
       int   midx   = truthColl.at(it).IndexMother();
       int   ampid  = abs(truthColl.at(midx).PdgId());
       float pt     = truthColl.at(it).Pt();
       float feta   = pt>1E-8? fabs(truthColl.at(it).Eta()): 0.;
       bool  HardScattered = GenSt>20 && GenSt<30;
       bool  PassQGenCut = apid<10 && pt>10. && feta<5.;
       bool  PassLGenCut = (apid==11 or apid==13) && pt>5. && feta<2.6;

       if( apid==5 && ampid==6  && HardScattered ){Nb_t++; if(PassQGenCut) Nb_t_hard++;}
       if( apid<10 && ampid==24 && HardScattered ){Nq_W++; if(PassQGenCut) Nq_W_hard++;}
       if( apid==13 && ampid==36                 ){Nm_A++; if(PassLGenCut) Nm_A_hard++;}
       if( apid>10 && apid<20   && ampid==24     ){
         Nl_W++; if(PassLGenCut) Nl_W_hard++;
         if( apid==11 ){Ne_W++; if(PassLGenCut) Ne_W_hard++;}
         if( apid==13 ){Nm_W++; if(PassLGenCut) Nm_W_hard++;} 
         if( apid==15 ){Nt_W++;}
         if( apid==12 or apid==14 or apid==16 ){Nv_W++;}
       }
     }
     if(Nq_W==2){
       Dec_lvqq=true;
       if     (Ne_W==1) Dec_evqq=true;
       else if(Nm_W==1) Dec_mvqq=true;
       else if(Nt_W==1) Dec_tvqq=true;
       else             Dec_else=true;
     }
     else if(Nq_W==0){
       Dec_llvv=true;
       if     (Ne_W==2            ) Dec_eevv=true;
       else if(Ne_W==1 and Nm_W==1) Dec_emvv=true;
       else if(Ne_W==1 and Nt_W==1) Dec_etvv=true;
       else if(Nm_W==2            ) Dec_mmvv=true;
       else if(Nm_W==1 and Nt_W==1) Dec_mtvv=true;
       else if(Nt_W==2            ) Dec_ttvv=true;
       else                         Dec_else=true;
     }
     else                           Dec_else=true;

     Dec_evqq_hard = (Nm_A_hard==2 && Ne_W_hard==1 && Nq_W_hard==2 && Nb_t_hard==2);
     Dec_mvqq_hard = (Nm_A_hard==2 && Nm_W_hard==1 && Nq_W_hard==2 && Nb_t_hard==2);


     FillHist("Nb_t", Nb_t, weight, 0., 10., 10); 
     FillHist("Nq_W", Nq_W, weight, 0., 10., 10); 
     FillHist("Nl_W", Nl_W, weight, 0., 10., 10); 
     FillHist("Nm_A", Nm_A, weight, 0., 10., 10); 
     if(Dec_lvqq){
       FillHist("DecMode", 0., weight, 0., 3., 3);
       if     (Dec_evqq) FillHist("DecMode_lvqq", 0., weight, 0., 6., 6);
       else if(Dec_mvqq) FillHist("DecMode_lvqq", 2., weight, 0., 6., 6);
       else if(Dec_tvqq) FillHist("DecMode_lvqq", 4., weight, 0., 6., 6);

       if     (Dec_evqq_hard) FillHist("DecMode_lvqq", 1., weight, 0., 6., 6);
       else if(Dec_mvqq_hard) FillHist("DecMode_lvqq", 3., weight, 0., 6., 6);
       else if(Dec_tvqq     ) FillHist("DecMode_lvqq", 5., weight, 0., 6., 6);
     }
     else if(Dec_llvv) FillHist("DecMode", 1., weight, 0., 3., 3);
     else if(Dec_else) FillHist("DecMode", 2., weight, 0., 3., 3);

     float LumiWeight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
     if(Dec_evqq     ) FillHist("NCutDecMode_Raw", 0., weight*LumiWeight, 0., 6., 6);
     if(Dec_evqq_hard) FillHist("NCutDecMode_Raw", 1., weight*LumiWeight, 0., 6., 6);
     if(Dec_mvqq     ) FillHist("NCutDecMode_Raw", 2., weight*LumiWeight, 0., 6., 6);
     if(Dec_mvqq_hard) FillHist("NCutDecMode_Raw", 3., weight*LumiWeight, 0., 6., 6);
     if(Dec_tvqq     ) FillHist("NCutDecMode_Raw", 4., weight*LumiWeight, 0., 6., 6);
     if(Dec_llvv     ) FillHist("NCutDecMode_Raw", 5., weight*LumiWeight, 0., 6., 6);

     if(Dec_evqq) FillHist("NDecMode_Raw", 0., weight*LumiWeight, 0., 9., 9);
     if(Dec_mvqq) FillHist("NDecMode_Raw", 1., weight*LumiWeight, 0., 9., 9);
     if(Dec_tvqq) FillHist("NDecMode_Raw", 2., weight*LumiWeight, 0., 9., 9);
     if(Dec_eevv) FillHist("NDecMode_Raw", 3., weight*LumiWeight, 0., 9., 9);
     if(Dec_emvv) FillHist("NDecMode_Raw", 4., weight*LumiWeight, 0., 9., 9);
     if(Dec_etvv) FillHist("NDecMode_Raw", 5., weight*LumiWeight, 0., 9., 9);
     if(Dec_mmvv) FillHist("NDecMode_Raw", 6., weight*LumiWeight, 0., 9., 9);
     if(Dec_mtvv) FillHist("NDecMode_Raw", 7., weight*LumiWeight, 0., 9., 9);
     if(Dec_ttvv) FillHist("NDecMode_Raw", 8., weight*LumiWeight, 0., 9., 9);
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
   if(false){ Pass_Trigger=true; trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi); }
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
     eventbase->GetMuonSel()->SetPt(5.);                     eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, true);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
     if(EMuMu)      { if( !(electronPreColl.size()>=1 && muonPreColl.size()>=2) ) return; }
     else if(TriMu) { if( !(muonPreColl.size()>=3) ) return; }
     FillCutFlow("PreSel", weight);
   //******************************************************************************************************//

   //Primary Object Collection
   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);

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
   std::vector<snu::KMuon> muonColl;  if(k_running_nonprompt){ muonColl=muonLooseColl;} else{ muonColl=muonTightColl;}


     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_HctoWA_FAKELOOSE);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.4, 0.4);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
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

   if(EventCand & !SystRun){
     if(!isData){
       if(EMuMu || TriMu){
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl);
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
    
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl);
//         trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);//No longer recommendation

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
     std::vector<snu::KMuon> MuTMuEnUpColl; eventbase->GetMuonSel()->Selection(MuTMuEnUpColl,false);
     //std::vector<snu::KMuon> MuTMuEnUpColl; eventbase->GetMuonSel()->Selection(MuTMuEnUpColl,true);
      // SetMuonResCorrection(MuTMuEnUpColl, "SystUp");

       eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
       eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
       eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
       eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
       eventbase->GetMuonSel()->SetdxySigMax(4.);
       eventbase->GetMuonSel()->SetChiNdof(4.);
     std::vector<snu::KMuon> MuTMuEnDownColl; eventbase->GetMuonSel()->Selection(MuTMuEnDownColl,false);
     //std::vector<snu::KMuon> MuTMuEnDownColl; eventbase->GetMuonSel()->Selection(MuTMuEnDownColl,true);
      // SetMuonResCorrection(MuTMuEnDownColl, "SystDown");


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
       //trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);//No longer recommendation
       btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

       id_weight_ele_sfup       = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl,  1);
       id_weight_ele_sfdown     = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl, -1);
       id_weight_mu_sfup        = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl,  1);
       id_weight_mu_sfdown      = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl, -1);

       id_weight_ele_ElEnup     = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", EleTElEnUpColl);
       reco_weight_ele_ElEnup   = mcdata_correction->ElectronRecoScaleFactor(EleTElEnUpColl);
       id_weight_mu_MuEnup      = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuTMuEnUpColl);
       //trk_weight_mu_MuEnup     = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnUpColl);//No longer recommendation

       id_weight_ele_ElEndown   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", EleTElEnDownColl);
       reco_weight_ele_ElEndown = mcdata_correction->ElectronRecoScaleFactor(EleTElEnDownColl);
       id_weight_mu_MuEndown    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuTMuEnDownColl);
       //trk_weight_mu_MuEndown   = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnDownColl);//No longer recommendation

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
         if(k_sample_name.Contains("ZGto2LG")){ ZG_SF=0.845; conv_up+=0.15; conv_down-=0.15; }
         else if(k_sample_name.Contains("TTG")){ conv_up+=0.5; conv_down-=0.5; }
       }
       if(EMuMu){
         if(k_sample_name.Contains("ZGto2LG")){ ZG_SF=0.961; conv_up+=0.092; conv_down-=0.092; }
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
   if(SRDist){
     if(EMuMu && !SystRun){
       CheckSRDist(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMu");
     }
     if(TriMu && !SystRun){
       CheckSRDist(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMu");
     }
   }
   if(SigEffValid){
     if(EMuMu && !SystRun){
       if(Dec_evqq     ) CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMuevqq");
       if(Dec_evqq_hard) CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMuevqq_hard");
       if(Dec_mvqq     ) CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMumvqq");
       if(Dec_mvqq_hard) CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMumvqq_hard");
       if(Dec_tvqq)      CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMutvqq");

       if(Dec_eevv)      CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMueevv");
       if(Dec_emvv)      CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMuemvv");
       if(Dec_etvv)      CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMuetvv");
       if(Dec_mmvv)      CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMummvv");
       if(Dec_mtvv)      CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMumtvv");
       if(Dec_ttvv)      CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMuttvv");
     }
     if(TriMu && !SystRun){
       if(Dec_evqq     ) CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMuevqq");
       if(Dec_evqq_hard) CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMuevqq_hard");
       if(Dec_mvqq     ) CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMumvqq");
       if(Dec_mvqq_hard) CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMumvqq_hard");
       if(Dec_tvqq)      CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMutvqq");

       if(Dec_eevv)      CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMueevv");
       if(Dec_emvv)      CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMuemvv");
       if(Dec_etvv)      CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMuetvv");
       if(Dec_mmvv)      CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMummvv");
       if(Dec_mtvv)      CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMumtvv");
       if(Dec_ttvv)      CompareSignalEfficiency(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMuttvv");
     }
   }
   if(JetVetoCheck){
     if(EMuMu){
       CheckJetVetoImpact(muonColl, muonLooseColl, electronColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMu");
     }
     if(TriMu){
       CheckJetVetoImpact(muonColl, muonLooseColl, electronColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMu");
     }
   }
   if(FkSrcCheck){
     if(EMuMu){
       CheckFakeSource(muonColl, muonLooseColl, electronColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight, "", "EMuMu");
     }
     if(TriMu){
       CheckFakeSource(muonColl, muonLooseColl, electronColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), weight, "", "TriMu");
     }
   }





     

/////////////////////////////////////////////////////////////////////////////////// 

return;

}// End of execute event loop
  


void Mar2019_AdditionalChecks::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Mar2019_AdditionalChecks::BeginCycle() throw( LQError ){
  
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

Mar2019_AdditionalChecks::~Mar2019_AdditionalChecks() {
  
  Message("In Mar2019_AdditionalChecks Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}


void Mar2019_AdditionalChecks::CompareSignalEfficiency(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");
  float FillCutBinIndex=-1., FillRawBinIndex=-1.;
  if(Option.Contains("evqq_hard")) FillCutBinIndex=1.;
  else if(Option.Contains("evqq")){FillCutBinIndex=0.; FillRawBinIndex=0.;}
  if(Option.Contains("mvqq_hard")) FillCutBinIndex=3.;
  else if(Option.Contains("mvqq")){FillCutBinIndex=2.; FillRawBinIndex=1.;}
  if(Option.Contains("tvqq"))     {FillCutBinIndex=4.; FillRawBinIndex=2.;}

  if(Option.Contains("eevv"))     {FillCutBinIndex=5.; FillRawBinIndex=3.;}
  if(Option.Contains("emvv"))     {FillCutBinIndex=5.; FillRawBinIndex=4.;}
  if(Option.Contains("etvv"))     {FillCutBinIndex=5.; FillRawBinIndex=5.;}
  if(Option.Contains("mmvv"))     {FillCutBinIndex=5.; FillRawBinIndex=6.;}
  if(Option.Contains("mtvv"))     {FillCutBinIndex=5.; FillRawBinIndex=7.;}
  if(Option.Contains("ttvv"))     {FillCutBinIndex=5.; FillRawBinIndex=8.;}
  FillCutBinIndex+=1E-3; FillRawBinIndex+=1E-3;
  


  const int NSignal = 7;
  float Mmumu_CentSeed[NSignal]  = {15., 25., 35., 45., 55., 65., 75.};
  //float Mmumu_WidthSeed[NSignal] = {1., 1., 1.3, 1.5, 1.5, 2.0, 2.5}; - Old - OffPreAppVer
  //float Mmumu_Step[NSignal-1]    = {1., 1., 1.2, 1.2, 1.2, 1.8};      - Old - OffPreAppVer
  float Mmumu_WidthSeed[NSignal] = {0.5,  0.7, 0.8,  1., 1.2, 1.5, 1.8};
  float Mmumu_Step[NSignal-1]    = {  0.45, 0.55, 0.6, 0.75, 0.9, 1.15};
  TString CentStr[NSignal]       = { "15", "25", "35", "45", "55", "65", "75"};

  vector<float> Mmumu_CentFull, Mmumu_WidthFull;
  for(int it_seed=0; it_seed<NSignal; it_seed++){
    float CurrentCenter = Mmumu_CentSeed[it_seed];
    float CurrentWidth  = Mmumu_WidthSeed[it_seed];
    Mmumu_CentFull.push_back(CurrentCenter);
    Mmumu_WidthFull.push_back(CurrentWidth);
    for(int it_step=0; it_step<100; it_step++){
      if(it_seed==NSignal-1) break;
      CurrentCenter+= Mmumu_Step[it_seed];
      CurrentWidth  = (CurrentCenter-Mmumu_CentSeed[it_seed])/10.*Mmumu_WidthSeed[it_seed+1]
                     +(Mmumu_CentSeed[it_seed+1]-CurrentCenter)/10.*Mmumu_WidthSeed[it_seed];
      if(CurrentCenter<Mmumu_CentSeed[it_seed]+10.){
        Mmumu_CentFull.push_back(CurrentCenter);
        Mmumu_WidthFull.push_back(CurrentWidth);
      }
      else break;
    }
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

    if(fabs(Mmumu-91.2)<10) return;

    float MAWidth=2.5;
    if(BJetColl.size()==0) return;


    if(JetColl.size()>1){

      FillHist("Count_2jCut_1e2mu_PreSel"+Label, 0., weight, 0., 1., 1);
      if( Mmumu>12. && Mmumu<80.  ) FillHist("Count_2jCut_1e2mu_M12to80"+Label, 0., weight, 0., 1., 1);

      for( int it_m=0; it_m<NSignal; it_m++){
        float MA=GetHiggsMass(k_sample_name, "A");
        if(fabs(MA-Mmumu_CentSeed[it_m])>1.) continue;
        float Shift = 0.;
        if( fabs(Mmumu-Mmumu_CentSeed[it_m])<Mmumu_WidthSeed[it_m] ){
          if(FillCutBinIndex>0.) FillHist("NCutDecMode_Final", FillCutBinIndex, weight, 0., 6., 6);
          if(FillRawBinIndex>0.) FillHist("NDecMode_Final", FillRawBinIndex, weight, 0., 9., 9);
        }
      }


      for(int it_cent=0; it_cent<Mmumu_CentFull.size(); it_cent++){
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<Mmumu_WidthFull.at(it_cent)){
          FillHist("Mmumu_CntBin_SR2j"+Label, it_cent+0.01, weight, 0., 95., 95);
        }
      }
    }
  }
  if(TriMu){
    if( !(MuLColl.size()==3 && EleLColl.size()==0) ) return;
    if( !(MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>20. && MuTColl.at(1).Pt()>10. && MuTColl.at(2).Pt()>10.) ) return;

    int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    int   IdxSSA = TriMuChargeIndex(MuTColl, MET, METx, METy, "SSA");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    float Mmumu  = (MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)).M();
    float MAWidth=2.5;
    if( !(MOSSS1>12. && MOSSS2>12.) ) return;
    if( fabs(MOSSS1-91.2)<10. || fabs(MOSSS2-91.2)<10. ) return;


    if( BJetColl.size()==0 ) return;

    if(JetColl.size()>1){
      FillHist("Count_2jCut_Mmumu_3mu_PreSel"+Label, 0., weight, 0., 1., 1);
      if(Mmumu<80) FillHist("Count_2jCut_Mmumu_3mu_M12to80"+Label, 0., weight, 0., 1., 1);

      for( int it_m=0; it_m<NSignal; it_m++){
        float MA=GetHiggsMass(k_sample_name, "A");
        if(fabs(MA-Mmumu_CentSeed[it_m])>1.) continue;
        if( fabs(Mmumu-Mmumu_CentSeed[it_m])<Mmumu_WidthSeed[it_m] ){
          if(FillCutBinIndex>0.) FillHist("NCutDecMode_Final", FillCutBinIndex, weight, 0., 6., 6);
          if(FillRawBinIndex>0.) FillHist("NDecMode_Final", FillRawBinIndex, weight, 0., 9., 9);
        }
      }
    }
  }
}


void Mar2019_AdditionalChecks::CheckJetVetoImpact(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetNoVetoColl, std::vector<snu::KJet> BJetNoVetoColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");

  std::vector<snu::KJet> JetColl  = SkimJetColl(JetNoVetoColl,  EleLColl, MuLColl, "EleMuVeto");
  std::vector<snu::KJet> BJetColl = SelBJets(JetColl, "Medium");


  const int NSignal = 7;
  float Mmumu_CentSeed[NSignal]  = {15., 25., 35., 45., 55., 65., 75.};
  float Mmumu_WidthSeed[NSignal] = {0.5,  0.7, 0.8,  1., 1.2, 1.5, 1.8};
  float Mmumu_Step[NSignal-1]    = {0.45, 0.55, 0.6, 0.75, 0.9, 1.15};
  TString CentStr[NSignal]       = {"15", "25", "35", "45", "55", "65", "75"};

  vector<float> Mmumu_CentFull, Mmumu_WidthFull;
  for(int it_seed=0; it_seed<NSignal; it_seed++){
    float CurrentCenter = Mmumu_CentSeed[it_seed];
    float CurrentWidth  = Mmumu_WidthSeed[it_seed];
    Mmumu_CentFull.push_back(CurrentCenter);
    Mmumu_WidthFull.push_back(CurrentWidth);
    for(int it_step=0; it_step<100; it_step++){
      if(it_seed==NSignal-1) break;
      CurrentCenter+= Mmumu_Step[it_seed];
      CurrentWidth  = (CurrentCenter-Mmumu_CentSeed[it_seed])/10.*Mmumu_WidthSeed[it_seed+1]
                     +(Mmumu_CentSeed[it_seed+1]-CurrentCenter)/10.*Mmumu_WidthSeed[it_seed];
      if(CurrentCenter<Mmumu_CentSeed[it_seed]+10.){
        Mmumu_CentFull.push_back(CurrentCenter);
        Mmumu_WidthFull.push_back(CurrentWidth);
      }
      else break;
    }
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

    if(fabs(Mmumu-91.2)<10) return;

    FillHist("Njet" ,  JetColl.size(), weight, 0., 10., 10);
    FillHist("Nbjet", BJetColl.size(), weight, 0., 10., 10);
    FillHist("NjetNoVeto" ,  JetNoVetoColl.size(), weight, 0., 10., 10);
    FillHist("NbjetNoVeto", BJetNoVetoColl.size(), weight, 0., 10., 10);


    float MAWidth=2.5;
    if(BJetColl.size()>0 && JetColl.size()>1){
      for( int it_m=0; it_m<NSignal; it_m++){
        float MA=GetHiggsMass(k_sample_name, "A");
        if(fabs(MA-Mmumu_CentSeed[it_m])>1.) continue;
        if( fabs(Mmumu-Mmumu_CentSeed[it_m])<Mmumu_WidthSeed[it_m] ){
          FillHist("NFinal_VetoTF", 0., weight, 0., 3., 3);
        }
      }
      if(!k_sample_name.Contains("HcToWA")){
        FillHist("NFinal_VetoTF", 0., weight, 0., 3., 3);
      }
    }
    if(BJetNoVetoColl.size()>0 && JetColl.size()>1){
      for( int it_m=0; it_m<NSignal; it_m++){
        float MA=GetHiggsMass(k_sample_name, "A");
        if(fabs(MA-Mmumu_CentSeed[it_m])>1.) continue;
        if( fabs(Mmumu-Mmumu_CentSeed[it_m])<Mmumu_WidthSeed[it_m] ){
          FillHist("NFinal_VetoTF", 1., weight, 0., 3., 3);
        }
      }
      if(!k_sample_name.Contains("HcToWA")){
        FillHist("NFinal_VetoTF", 1., weight, 0., 3., 3);
      }
    }
    if(BJetNoVetoColl.size()>0 && JetNoVetoColl.size()>1){
      for( int it_m=0; it_m<NSignal; it_m++){
        float MA=GetHiggsMass(k_sample_name, "A");
        if(fabs(MA-Mmumu_CentSeed[it_m])>1.) continue;
        if( fabs(Mmumu-Mmumu_CentSeed[it_m])<Mmumu_WidthSeed[it_m] ){
          FillHist("NFinal_VetoTF", 2., weight, 0., 3., 3);
        }
      }
      if(!k_sample_name.Contains("HcToWA")){
        FillHist("NFinal_VetoTF", 2., weight, 0., 3., 3);
      }
    }
  }
  if(TriMu){
    if( !(MuLColl.size()==3 && EleLColl.size()==0) ) return;
    if( !(MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>20. && MuTColl.at(1).Pt()>10. && MuTColl.at(2).Pt()>10.) ) return;

    int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    int   IdxSSA = TriMuChargeIndex(MuTColl, MET, METx, METy, "SSA");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    float Mmumu  = (MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)).M();
    float MAWidth=2.5;
    if( !(MOSSS1>12. && MOSSS2>12.) ) return;
    if( fabs(MOSSS1-91.2)<10. || fabs(MOSSS2-91.2)<10. ) return;


    FillHist("Njet" ,  JetColl.size(), weight, 0., 10., 10);
    FillHist("Nbjet", BJetColl.size(), weight, 0., 10., 10);
    FillHist("NjetNoVeto" ,  JetNoVetoColl.size(), weight, 0., 10., 10);
    FillHist("NbjetNoVeto", BJetNoVetoColl.size(), weight, 0., 10., 10);


    if(BJetColl.size()>0 && JetColl.size()>1){
      for( int it_m=0; it_m<NSignal; it_m++){
        float MA=GetHiggsMass(k_sample_name, "A");
        if(fabs(MA-Mmumu_CentSeed[it_m])>1.) continue;
        if( fabs(Mmumu-Mmumu_CentSeed[it_m])<Mmumu_WidthSeed[it_m] ){
          FillHist("NFinal_VetoTF", 0., weight, 0., 3., 3);
        }
      }
      if(!k_sample_name.Contains("HcToWA")){
        FillHist("NFinal_VetoTF", 0., weight, 0., 3., 3);
      }
    }
    if(BJetNoVetoColl.size()>0 && JetColl.size()>1){
      for( int it_m=0; it_m<NSignal; it_m++){
        float MA=GetHiggsMass(k_sample_name, "A");
        if(fabs(MA-Mmumu_CentSeed[it_m])>1.) continue;
        if( fabs(Mmumu-Mmumu_CentSeed[it_m])<Mmumu_WidthSeed[it_m] ){
          FillHist("NFinal_VetoTF", 1., weight, 0., 3., 3);
        }
      }
      if(!k_sample_name.Contains("HcToWA")){
        FillHist("NFinal_VetoTF", 1., weight, 0., 3., 3);
      }
    }
    if(BJetNoVetoColl.size()>0 && JetNoVetoColl.size()>1){
      for( int it_m=0; it_m<NSignal; it_m++){
        float MA=GetHiggsMass(k_sample_name, "A");
        if(fabs(MA-Mmumu_CentSeed[it_m])>1.) continue;
        if( fabs(Mmumu-Mmumu_CentSeed[it_m])<Mmumu_WidthSeed[it_m] ){
          FillHist("NFinal_VetoTF", 2., weight, 0., 3., 3);
        }
      }
      if(!k_sample_name.Contains("HcToWA")){
        FillHist("NFinal_VetoTF", 2., weight, 0., 3., 3);
      }
    }
  }
}



void Mar2019_AdditionalChecks::CheckFakeSource(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetNoVetoColl, std::vector<snu::KJet> BJetNoVetoColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");
  std::vector<snu::KJet> JetColl  = SkimJetColl(JetNoVetoColl,  EleLColl, MuLColl, "EleMuVeto");
  std::vector<snu::KJet> BJetColl = SelBJets(JetColl, "Medium");
  if(isData) return;


  const int NSignal = 7;
  float Mmumu_CentSeed[NSignal]  = {15., 25., 35., 45., 55., 65., 75.};
  float Mmumu_WidthSeed[NSignal] = {0.5,  0.7, 0.8,  1., 1.2, 1.5, 1.8};
  float Mmumu_Step[NSignal-1]    = {0.45, 0.55, 0.6, 0.75, 0.9, 1.15};
  TString CentStr[NSignal]       = {"15", "25", "35", "45", "55", "65", "75"};

  vector<float> Mmumu_CentFull, Mmumu_WidthFull;
  for(int it_seed=0; it_seed<NSignal; it_seed++){
    float CurrentCenter = Mmumu_CentSeed[it_seed];
    float CurrentWidth  = Mmumu_WidthSeed[it_seed];
    Mmumu_CentFull.push_back(CurrentCenter);
    Mmumu_WidthFull.push_back(CurrentWidth);
    for(int it_step=0; it_step<100; it_step++){
      if(it_seed==NSignal-1) break;
      CurrentCenter+= Mmumu_Step[it_seed];
      CurrentWidth  = (CurrentCenter-Mmumu_CentSeed[it_seed])/10.*Mmumu_WidthSeed[it_seed+1]
                     +(Mmumu_CentSeed[it_seed+1]-CurrentCenter)/10.*Mmumu_WidthSeed[it_seed];
      if(CurrentCenter<Mmumu_CentSeed[it_seed]+10.){
        Mmumu_CentFull.push_back(CurrentCenter);
        Mmumu_WidthFull.push_back(CurrentWidth);
      }
      else break;
    }
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

    if(fabs(Mmumu-91.2)<10) return;

    std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);

    int LeptonType_el  = GetLeptonType(EleTColl.at(0), truthColl);
    int LeptonType_mu1 = GetLeptonType(MuTColl.at(0), truthColl);
    int LeptonType_mu2 = GetLeptonType(MuTColl.at(1), truthColl);

    int JetType_el         = GetFakeLepSrcJetType(EleTColl.at(0), JetNoVetoColl);
    int MatchedTruthIdx_el = GenMatchedIdx(EleTColl.at(0),truthColl);
    int LastSelfIdx_el     = LastSelfMotherIdx(MatchedTruthIdx_el,truthColl);
    int MotherIdx_el       = FirstNonSelfMotherIdx(MatchedTruthIdx_el,truthColl);
    int LastSelfMIdx_el    = LastSelfMotherIdx(MotherIdx_el,truthColl);
    int GrMotherIdx_el     = FirstNonSelfMotherIdx(MotherIdx_el,truthColl);
    int LastSelfGrMIdx_el  = LastSelfMotherIdx(GrMotherIdx_el,truthColl);
    int MPID_el=0, GrMPID_el=0;
    if(  MotherIdx_el!=-1  ) MPID_el   = truthColl.at(MotherIdx_el).PdgId();
    if( GrMotherIdx_el!=-1 ) GrMPID_el = truthColl.at(GrMotherIdx_el).PdgId();

    int JetType_mu1         = GetFakeLepSrcJetType(MuTColl.at(0), JetNoVetoColl);
    int MatchedTruthIdx_mu1 = GenMatchedIdx(MuTColl.at(0),truthColl);
    int LastSelfIdx_mu1     = LastSelfMotherIdx(MatchedTruthIdx_mu1,truthColl);
    int MotherIdx_mu1       = FirstNonSelfMotherIdx(MatchedTruthIdx_mu1,truthColl);
    int LastSelfMIdx_mu1    = LastSelfMotherIdx(MotherIdx_mu1,truthColl);
    int GrMotherIdx_mu1     = FirstNonSelfMotherIdx(MotherIdx_mu1,truthColl);
    int LastSelfGrMIdx_mu1  = LastSelfMotherIdx(GrMotherIdx_mu1,truthColl);
    int MPID_mu1=0, GrMPID_mu1=0;
    if(  MotherIdx_mu1!=-1  ) MPID_mu1   = truthColl.at(MotherIdx_mu1).PdgId();
    if( GrMotherIdx_mu1!=-1 ) GrMPID_mu1 = truthColl.at(GrMotherIdx_mu1).PdgId();

    int JetType_mu2         = GetFakeLepSrcJetType(MuTColl.at(1), JetNoVetoColl);
    int MatchedTruthIdx_mu2 = GenMatchedIdx(MuTColl.at(1),truthColl);
    int LastSelfIdx_mu2     = LastSelfMotherIdx(MatchedTruthIdx_mu2,truthColl);
    int MotherIdx_mu2       = FirstNonSelfMotherIdx(MatchedTruthIdx_mu2,truthColl);
    int LastSelfMIdx_mu2    = LastSelfMotherIdx(MotherIdx_mu2,truthColl);
    int GrMotherIdx_mu2     = FirstNonSelfMotherIdx(MotherIdx_mu2,truthColl);
    int LastSelfGrMIdx_mu2  = LastSelfMotherIdx(GrMotherIdx_mu2,truthColl);
    int MPID_mu2=0, GrMPID_mu2=0;
    if(  MotherIdx_mu2!=-1  ) MPID_mu2   = truthColl.at(MotherIdx_mu2).PdgId();
    if( GrMotherIdx_mu2!=-1 ) GrMPID_mu2 = truthColl.at(GrMotherIdx_mu2).PdgId();


    if(LeptonType_el<0 && LeptonType_el>-5){
      FillHist("LepType_FkEl_NoJCut", LeptonType_el, weight, -10., 10., 20);
      FillHist("JetType_FkEl_NoJCut", JetType_el, weight, -5., 5., 10);
      if(LeptonType_el==-2){
        FillHist("MHadIdx_FkEl_NoJCut", abs(MPID_el), weight, 0., 10000000., 1E7);
      }
      if(LeptonType_el==-3 && GrMPID_el>0){
        FillHist("MHadIdx_FkEl_NoJCut", abs(GrMPID_el), weight, 0., 10000000., 1E7);
      }
    }
    if(LeptonType_mu1<0 && LeptonType_mu1>-5){
      FillHist("LepType_FkMu_NoJCut", LeptonType_mu1, weight, -10., 10., 20);
      FillHist("JetType_FkMu_NoJCut", JetType_mu1, weight, -5., 5., 10);
      if(LeptonType_mu1==-2){
        FillHist("MHadIdx_FkMu_NoJCut", abs(MPID_mu1), weight, 0., 10000000., 1E7);
      }
      if(LeptonType_mu1==-3 && GrMPID_mu1>0){
        FillHist("MHadIdx_FkMu_NoJCut", abs(GrMPID_mu1), weight, 0., 10000000., 1E7);
      }
    }
    if(LeptonType_mu2<0 && LeptonType_mu2>-5){
      FillHist("LepType_FkMu_NoJCut", LeptonType_mu2, weight, -10., 10., 20);
      FillHist("JetType_FkMu_NoJCut", JetType_mu2, weight, -5., 5., 10);
      if(LeptonType_mu2==-2){
        FillHist("MHadIdx_FkMu_NoJCut", abs(MPID_mu2), weight, 0., 10000000., 1E7);
      }
      if(LeptonType_mu2==-3 && GrMPID_mu2>0){
        FillHist("MHadIdx_FkMu_NoJCut", abs(GrMPID_mu2), weight, 0., 10000000., 1E7);
      }
    }


    if(JetColl.size()<2) return;
    if(BJetColl.size()==0) return;


    if(LeptonType_el<0 && LeptonType_el>-5){
      FillHist("LepType_FkEl_JCut", LeptonType_el, weight, -10., 10., 20);
      FillHist("JetType_FkEl_JCut", JetType_el, weight, -5., 5., 10);
      if(LeptonType_el==-2){
        FillHist("MHadIdx_FkEl_JCut", abs(MPID_el), weight, 0., 10000000., 1E7);
      }
      if(LeptonType_el==-3 && GrMPID_el>0){
        FillHist("MHadIdx_FkEl_JCut", abs(GrMPID_el), weight, 0., 10000000., 1E7);
      }
    }
    if(LeptonType_mu1<0 && LeptonType_mu1>-5){
      FillHist("LepType_FkMu_JCut", LeptonType_mu1, weight, -10., 10., 20);
      FillHist("JetType_FkMu_JCut", JetType_mu1, weight, -5., 5., 10);
      if(LeptonType_mu1==-2){
        FillHist("MHadIdx_FkMu_JCut", abs(MPID_mu1), weight, 0., 10000000., 1E7);
      }
      if(LeptonType_mu1==-3 && GrMPID_mu1>0){
        FillHist("MHadIdx_FkMu_JCut", abs(GrMPID_mu1), weight, 0., 10000000., 1E7);
      }
    }
    if(LeptonType_mu2<0 && LeptonType_mu2>-5){
      FillHist("LepType_FkMu_JCut", LeptonType_mu2, weight, -10., 10., 20);
      FillHist("JetType_FkMu_JCut", JetType_mu2, weight, -5., 5., 10);
      if(LeptonType_mu2==-2){
        FillHist("MHadIdx_FkMu_JCut", abs(MPID_mu2), weight, 0., 10000000., 1E7);
      }
      if(LeptonType_mu2==-3 && GrMPID_mu2>0){
        FillHist("MHadIdx_FkMu_JCut", abs(GrMPID_mu2), weight, 0., 10000000., 1E7);
      }
    }
  }
  if(TriMu){
    if( !(MuLColl.size()==3 && EleLColl.size()==0) ) return;
    if( !(MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>20. && MuTColl.at(1).Pt()>10. && MuTColl.at(2).Pt()>10.) ) return;

    int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    int   IdxSSA = TriMuChargeIndex(MuTColl, MET, METx, METy, "SSA");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    float Mmumu  = (MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)).M();
    float MAWidth=2.5;
    if( !(MOSSS1>12. && MOSSS2>12.) ) return;
    if( fabs(MOSSS1-91.2)<10. || fabs(MOSSS2-91.2)<10. ) return;


    std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);

    int LeptonType_mu1 = GetLeptonType(MuTColl.at(0), truthColl);
    int LeptonType_mu2 = GetLeptonType(MuTColl.at(1), truthColl);
    int LeptonType_mu3 = GetLeptonType(MuTColl.at(2), truthColl);

    int JetType_mu1         = GetFakeLepSrcJetType(MuTColl.at(0), JetNoVetoColl);
    int MatchedTruthIdx_mu1 = GenMatchedIdx(MuTColl.at(0),truthColl);
    int LastSelfIdx_mu1     = LastSelfMotherIdx(MatchedTruthIdx_mu1,truthColl);
    int MotherIdx_mu1       = FirstNonSelfMotherIdx(MatchedTruthIdx_mu1,truthColl);
    int LastSelfMIdx_mu1    = LastSelfMotherIdx(MotherIdx_mu1,truthColl);
    int GrMotherIdx_mu1     = FirstNonSelfMotherIdx(MotherIdx_mu1,truthColl);
    int LastSelfGrMIdx_mu1  = LastSelfMotherIdx(GrMotherIdx_mu1,truthColl);
    int MPID_mu1=0, GrMPID_mu1=0;
    if(  MotherIdx_mu1!=-1  ) MPID_mu1   = truthColl.at(MotherIdx_mu1).PdgId();
    if( GrMotherIdx_mu1!=-1 ) GrMPID_mu1 = truthColl.at(GrMotherIdx_mu1).PdgId();

    int JetType_mu2         = GetFakeLepSrcJetType(MuTColl.at(1), JetNoVetoColl);
    int MatchedTruthIdx_mu2 = GenMatchedIdx(MuTColl.at(1),truthColl);
    int LastSelfIdx_mu2     = LastSelfMotherIdx(MatchedTruthIdx_mu2,truthColl);
    int MotherIdx_mu2       = FirstNonSelfMotherIdx(MatchedTruthIdx_mu2,truthColl);
    int LastSelfMIdx_mu2    = LastSelfMotherIdx(MotherIdx_mu2,truthColl);
    int GrMotherIdx_mu2     = FirstNonSelfMotherIdx(MotherIdx_mu2,truthColl);
    int LastSelfGrMIdx_mu2  = LastSelfMotherIdx(GrMotherIdx_mu2,truthColl);
    int MPID_mu2=0, GrMPID_mu2=0;
    if(  MotherIdx_mu2!=-1  ) MPID_mu2   = truthColl.at(MotherIdx_mu2).PdgId();
    if( GrMotherIdx_mu2!=-1 ) GrMPID_mu2 = truthColl.at(GrMotherIdx_mu2).PdgId();

    int JetType_mu3         = GetFakeLepSrcJetType(MuTColl.at(2), JetNoVetoColl);
    int MatchedTruthIdx_mu3 = GenMatchedIdx(MuTColl.at(2),truthColl);
    int LastSelfIdx_mu3     = LastSelfMotherIdx(MatchedTruthIdx_mu3,truthColl);
    int MotherIdx_mu3       = FirstNonSelfMotherIdx(MatchedTruthIdx_mu3,truthColl);
    int LastSelfMIdx_mu3    = LastSelfMotherIdx(MotherIdx_mu3,truthColl);
    int GrMotherIdx_mu3     = FirstNonSelfMotherIdx(MotherIdx_mu3,truthColl);
    int LastSelfGrMIdx_mu3  = LastSelfMotherIdx(GrMotherIdx_mu3,truthColl);
    int MPID_mu3=0, GrMPID_mu3=0;
    if(  MotherIdx_mu3!=-1  ) MPID_mu3   = truthColl.at(MotherIdx_mu3).PdgId();
    if( GrMotherIdx_mu3!=-1 ) GrMPID_mu3 = truthColl.at(GrMotherIdx_mu3).PdgId();


    if(LeptonType_mu1<0 && LeptonType_mu1>-5){
      FillHist("LepType_FkMu_NoJCut", LeptonType_mu1, weight, -10., 10., 20);
      FillHist("JetType_FkMu_NoJCut", JetType_mu1, weight, -5., 5., 10);
      if(LeptonType_mu1==-2){
        FillHist("MHadIdx_FkMu_NoJCut", abs(MPID_mu1), weight, 0., 10000000., 1E7);
      }
      if(LeptonType_mu1==-3 && GrMPID_mu1>0){
        FillHist("MHadIdx_FkMu_NoJCut", abs(GrMPID_mu1), weight, 0., 10000000., 1E7);
      }
    }
    if(LeptonType_mu2<0 && LeptonType_mu2>-5){
      FillHist("LepType_FkMu_NoJCut", LeptonType_mu2, weight, -10., 10., 20);
      FillHist("JetType_FkMu_NoJCut", JetType_mu2, weight, -5., 5., 10);
      if(LeptonType_mu2==-2){
        FillHist("MHadIdx_FkMu_NoJCut", abs(MPID_mu2), weight, 0., 10000000., 1E7);
      }
      if(LeptonType_mu2==-3 && GrMPID_mu2>0){
        FillHist("MHadIdx_FkMu_NoJCut", abs(GrMPID_mu2), weight, 0., 10000000., 1E7);
      }
    }
    if(LeptonType_mu3<0 && LeptonType_mu3>-5){
      FillHist("LepType_FkMu_NoJCut", LeptonType_mu3, weight, -10., 10., 20);
      FillHist("JetType_FkMu_NoJCut", JetType_mu3, weight, -5., 5., 10);
      if(LeptonType_mu3==-2){
        FillHist("MHadIdx_FkMu_NoJCut", abs(MPID_mu3), weight, 0., 10000000., 1E7);
      }
      if(LeptonType_mu3==-3 && GrMPID_mu3>0){
        FillHist("MHadIdx_FkMu_NoJCut", abs(GrMPID_mu3), weight, 0., 10000000., 1E7);
      }
    }


    if(JetColl.size()<2) return;
    if(BJetColl.size()==0) return;


    if(LeptonType_mu1<0 && LeptonType_mu1>-5){
      FillHist("LepType_FkMu_JCut", LeptonType_mu1, weight, -10., 10., 20);
      FillHist("JetType_FkMu_JCut", JetType_mu1, weight, -5., 5., 10);
      if(LeptonType_mu1==-2){
        FillHist("MHadIdx_FkMu_JCut", abs(MPID_mu1), weight, 0., 10000000., 1E7);
      }
      if(LeptonType_mu1==-3 && GrMPID_mu1>0){
        FillHist("MHadIdx_FkMu_JCut", abs(GrMPID_mu1), weight, 0., 10000000., 1E7);
      }
    }
    if(LeptonType_mu2<0 && LeptonType_mu2>-5){
      FillHist("LepType_FkMu_JCut", LeptonType_mu2, weight, -10., 10., 20);
      FillHist("JetType_FkMu_JCut", JetType_mu2, weight, -5., 5., 10);
      if(LeptonType_mu2==-2){
        FillHist("MHadIdx_FkMu_JCut", abs(MPID_mu2), weight, 0., 10000000., 1E7);
      }
      if(LeptonType_mu2==-3 && GrMPID_mu2>0){
        FillHist("MHadIdx_FkMu_JCut", abs(GrMPID_mu2), weight, 0., 10000000., 1E7);
      }
    }
    if(LeptonType_mu3<0 && LeptonType_mu3>-5){
      FillHist("LepType_FkMu_JCut", LeptonType_mu3, weight, -10., 10., 20);
      FillHist("JetType_FkMu_JCut", JetType_mu3, weight, -5., 5., 10);
      if(LeptonType_mu3==-2){
        FillHist("MHadIdx_FkMu_JCut", abs(MPID_mu3), weight, 0., 10000000., 1E7);
      }
      if(LeptonType_mu3==-3 && GrMPID_mu3>0){
        FillHist("MHadIdx_FkMu_JCut", abs(GrMPID_mu3), weight, 0., 10000000., 1E7);
      }
    }
  }
}






void Mar2019_AdditionalChecks::CheckSRYield(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");

  const int NSignal = 7;
  float Mmumu_CentSeed[NSignal]  = {15., 25., 35., 45., 55., 65., 75.};
  //float Mmumu_WidthSeed[NSignal] = {1., 1., 1.3, 1.5, 1.5, 2.0, 2.5}; - Old - OffPreAppVer
  //float Mmumu_Step[NSignal-1]    = {1., 1., 1.2, 1.2, 1.2, 1.8};      - Old - OffPreAppVer
  float Mmumu_WidthSeed[NSignal] = {0.5,  0.7, 0.8,  1., 1.2, 1.5, 1.8};
  float Mmumu_Step[NSignal-1]    = {  0.45, 0.55, 0.6, 0.75, 0.9, 1.15};
  TString CentStr[NSignal]       = { "15", "25", "35", "45", "55", "65", "75"};

  vector<float> Mmumu_CentFull, Mmumu_WidthFull;
  for(int it_seed=0; it_seed<NSignal; it_seed++){
    float CurrentCenter = Mmumu_CentSeed[it_seed];
    float CurrentWidth  = Mmumu_WidthSeed[it_seed];
    Mmumu_CentFull.push_back(CurrentCenter);
    Mmumu_WidthFull.push_back(CurrentWidth);
    for(int it_step=0; it_step<100; it_step++){
      if(it_seed==NSignal-1) break;
      CurrentCenter+= Mmumu_Step[it_seed];
      CurrentWidth  = (CurrentCenter-Mmumu_CentSeed[it_seed])/10.*Mmumu_WidthSeed[it_seed+1]
                     +(Mmumu_CentSeed[it_seed+1]-CurrentCenter)/10.*Mmumu_WidthSeed[it_seed];
      if(CurrentCenter<Mmumu_CentSeed[it_seed]+10.){
        Mmumu_CentFull.push_back(CurrentCenter);
        Mmumu_WidthFull.push_back(CurrentWidth);
      }
      else break;
    }
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

    if(fabs(Mmumu-91.2)<10) return;

    float MAWidth=2.5;
    //V1:  1./  1./ 1.5/ 1.5/ 1.5/ 2./ 2.5
    //V2: 0.7/  1./ 1.3/ 1.5/ 1.5/ 2./ 2.5
    //V3: 0.7/ 0.9/ 0.9/ 1.5/ 1.5/ 2./ 2.5
    //V4:  1./  1./ 1.3/ 1.5/ 1.5/ 2./ 2.5
    if(JetColl.size()>1){
      if(Label=="") FillHist("Nb_2jCutPreB_1e2mu"+Label, BJetColl.size(), weight, 0., 10., 10);

      FillHist("Count_2jCutPreB_1e2mu_PreSel"+Label, 0., weight, 0., 1., 1);
      if( Mmumu>12. && Mmumu<80 ) FillHist("Count_2jCutPreB_1e2mu_M12to80"+Label, 0., weight, 0., 1., 1);

      for( int it_m=0; it_m<NSignal; it_m++){
        float Shift  = 0.;
        if( fabs(Mmumu-Mmumu_CentSeed[it_m]-Shift)<Mmumu_WidthSeed[it_m] ){
          FillHist("Count_2jCutPreB_1e2mu_M"+CentStr[it_m]+Label, 0., weight, 0., 1., 1);
        }
      }
    }

    if(BJetColl.size()==0) return;


    if(Label==""){
      FillHist("Nj_PreSel_1e2mu"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Mmumu_PreSel_1e2mu"+Label, Mmumu, weight, 0., 200., 40);
    }
    //c.f. Regarding detector resolution not considering FSR tail
    //σ(m, pT) = 0.026 + 0.0065m for the barrel,
    //σ(m, pT) = 0.026 + 0.013m GeV/c2 for the endcap.
    if(JetColl.size()>1){
      if(Label=="") FillHist("Mmumu_2jCut_1e2mu"+Label, Mmumu, weight, 0., 200., 40);

      FillHist("Count_2jCut_1e2mu_PreSel"+Label, 0., weight, 0., 1., 1);
      if( Mmumu>12. && Mmumu<80.  ) FillHist("Count_2jCut_1e2mu_M12to80"+Label, 0., weight, 0., 1., 1);

      for( int it_m=0; it_m<NSignal; it_m++){
        float Shift = 0.;
        if( fabs(Mmumu-Mmumu_CentSeed[it_m]-Shift)<Mmumu_WidthSeed[it_m] ){
          FillHist("Count_2jCut_M"+CentStr[it_m]+Label, 0., weight, 0., 1., 1);
          //FillHist("Count_2jCut_1e2mu_M"+CentStr[it_m]+Label, 0., weight, 0., 1., 1);
        }
      }

      for(int it_cent=0; it_cent<Mmumu_CentFull.size(); it_cent++){
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<Mmumu_WidthFull.at(it_cent)){
          FillHist("Mmumu_CntBin_SR2j"+Label, it_cent+0.01, weight, 0., 95., 95);
        }
  
        float LargeWidth=5.;
        if((Mmumu_CentFull.at(it_cent)-LargeWidth)<12.) LargeWidth=Mmumu_CentFull.at(it_cent)-12.;
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidth){
          FillHist("Mmumu_Width10_SR2j"+Label, it_cent+0.01, weight, 0., 95., 95);
        }
        if(Label==""){
          float Variation=0.2;
          float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
          if((Mmumu_CentFull.at(it_cent)-LargeWidthUp)<12.) LargeWidthUp=Mmumu_CentFull.at(it_cent)-12.;
          if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthUp){
            FillHist("Mmumu_Width10_SR2j_systup_Mrange", it_cent+0.01, weight, 0., 95., 95);
          }
          if((Mmumu_CentFull.at(it_cent)-LargeWidthDown)<12.) LargeWidthDown=Mmumu_CentFull.at(it_cent)-12.;
          if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthDown){
            FillHist("Mmumu_Width10_SR2j_systdown_Mrange", it_cent+0.01, weight, 0., 95., 95);
          }
        }
      }
    }
  }
  if(TriMu){
    if( !(MuLColl.size()==3 && EleLColl.size()==0) ) return;
    if( !(MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>20. && MuTColl.at(1).Pt()>10. && MuTColl.at(2).Pt()>10.) ) return;

    int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    int   IdxSSA = TriMuChargeIndex(MuTColl, MET, METx, METy, "SSA");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    //float Mmumu  = MOSSS2;
    float Mmumu  = (MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)).M();
    float MAWidth=2.5;
    if( !(MOSSS1>12. && MOSSS2>12.) ) return;
    if( fabs(MOSSS1-91.2)<10. || fabs(MOSSS2-91.2)<10. ) return;

    if(JetColl.size()>1){
      if(Label=="") FillHist("Nb_2jCutPreB_3mu"+Label, BJetColl.size(), weight, 0., 10., 10);

      FillHist("Count_2jCutPreB_Mmumu_3mu_PreSel"+Label, 0., weight, 0., 1., 1);
      if(Mmumu<80.) FillHist("Count_2jCutPreB_Mmumu_3mu_M12to80"+Label, 0., weight, 0., 1., 1);

      for( int it_m=0; it_m<NSignal; it_m++){
        float Shift = 0.;
        if( fabs(Mmumu-Mmumu_CentSeed[it_m]-Shift)<Mmumu_WidthSeed[it_m] ){
          FillHist("Count_2jCutPreB_Mmumu_3mu_M"+CentStr[it_m]+Label, 0., weight, 0., 1., 1);
        }
      }
    }

    if( BJetColl.size()==0 ) return;

    if(JetColl.size()>1){
      if(Label==""){
        FillHist("MOSSS1_2j"+Label, MOSSS1, weight, 0., 200., 10);
        FillHist("MOSSS2_2j"+Label, MOSSS2, weight, 0., 200., 20);
        FillHist("Mmumu_2j" +Label, Mmumu , weight, 0., 200., 20);
      }
      FillHist("Count_2jCut_Mmumu_3mu_PreSel"+Label, 0., weight, 0., 1., 1);
      if(Mmumu<80) FillHist("Count_2jCut_Mmumu_3mu_M12to80"+Label, 0., weight, 0., 1., 1);

      for( int it_m=0; it_m<NSignal; it_m++){
        float Shift = 0.;
        if( fabs(Mmumu-Mmumu_CentSeed[it_m]-Shift)<Mmumu_WidthSeed[it_m] ){
          FillHist("Count_2jCut_M"+CentStr[it_m]+Label, 0., weight, 0., 1., 1);
          //FillHist("Count_2jCut_Mmumu_3mu_M"+CentStr[it_m]+Label, 0., weight, 0., 1., 1);
        }
      }

      for(int it_cent=0; it_cent<Mmumu_CentFull.size(); it_cent++){
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<Mmumu_WidthFull.at(it_cent)){
          FillHist("Mmumu_CntBin_SR2j"+Label, it_cent+0.01, weight, 0., 95., 95);
        }
  
        float LargeWidth=5.;
        if((Mmumu_CentFull.at(it_cent)-LargeWidth)<12.) LargeWidth=Mmumu_CentFull.at(it_cent)-12.;
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidth){
          FillHist("Mmumu_Width10_SR2j"+Label, it_cent+0.01, weight, 0., 95., 95);
        }
        if(Label==""){
          float Variation=0.2;
          float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
          if((Mmumu_CentFull.at(it_cent)-LargeWidthUp)<12.) LargeWidthUp=Mmumu_CentFull.at(it_cent)-12.;
          if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthUp){
            FillHist("Mmumu_Width10_SR2j_systup_Mrange", it_cent+0.01, weight, 0., 95., 95);
          }
          if((Mmumu_CentFull.at(it_cent)-LargeWidthDown)<12.) LargeWidthDown=Mmumu_CentFull.at(it_cent)-12.;
          if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthDown){
            FillHist("Mmumu_Width10_SR2j_systdown_Mrange", it_cent+0.01, weight, 0., 95., 95);
          }
        }
      }
    }
  }
}




void Mar2019_AdditionalChecks::CheckSRDist(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool EMuMu=Option.Contains("EMuMu"), TriMu=Option.Contains("TriMu");
  if(k_sample_name.Contains("TTToHcToWA") && Label!="") return;//No systrun for signals

  vector<float> MmumuEdges, MmumuEdges2, MmumuEdges3, MmumuEdges4;
  float CurrentEdge=0.;
   for(int i=0; i<19; i++){
     if     (i==0)  CurrentEdge=12.;
     else if(i==1)  CurrentEdge=20.;
     else if(i<7)   CurrentEdge+=10.;
     else if(i==7)  CurrentEdge=81.2;
     else if(i==8)  CurrentEdge=101.2;     
     else if(i==9)  CurrentEdge=110.;
     else if(i<18)  CurrentEdge+=10;
     else CurrentEdge=200.;
     MmumuEdges.push_back(CurrentEdge);
   }
   for(int i=0; i<36; i++){
     if     (i==0)  CurrentEdge=12.;
     else if(i<35)  CurrentEdge+=2.;
     else if(i==35) CurrentEdge=81.2;
     MmumuEdges2.push_back(CurrentEdge);
   }
   for(int i=0; i<70; i++){
     if     (i==0)               CurrentEdge =12.;
     else if(CurrentEdge<80.  )  CurrentEdge+=1.;
     else if(CurrentEdge<81.2 )  CurrentEdge =81.2;
     else continue;
     MmumuEdges3.push_back(CurrentEdge);
   }
   for(int i=0; i<114; i++){
     if     (i==0)                 CurrentEdge=12.;
     else if(CurrentEdge<20.-1E-3) CurrentEdge+=0.25;
     else if(CurrentEdge<40.-1E-3) CurrentEdge+=0.5;
     else if(CurrentEdge<80.-1E-3) CurrentEdge+= 1.;
     else if(CurrentEdge<81.2 )    CurrentEdge=81.2;
     else continue;
     MmumuEdges4.push_back(CurrentEdge);
   }


  const int NSignal = 7;
  float Mmumu_CentSeed[NSignal]  = {15., 25., 35., 45., 55., 65., 75.};
  //float Mmumu_WidthSeed[NSignal] = {1., 1., 1.3, 1.5, 1.5, 2.0, 2.5};//Old - OffPreAppVer
  //float Mmumu_Step[NSignal-1]    = {1., 1., 1.2, 1.2, 1.2, 1.8};//Old - OffPreAppVer
  float Mmumu_WidthSeed[NSignal] = {0.5,  0.7, 0.8,  1., 1.2, 1.5, 1.8};
  float Mmumu_Step[NSignal-1]    = {  0.45, 0.55, 0.6, 0.75, 0.9, 1.15};

  vector<float> Mmumu_CentFull, Mmumu_WidthFull;
  for(int it_seed=0; it_seed<NSignal; it_seed++){
    float CurrentCenter = Mmumu_CentSeed[it_seed];
    float CurrentWidth  = Mmumu_WidthSeed[it_seed];
    Mmumu_CentFull.push_back(CurrentCenter);
    Mmumu_WidthFull.push_back(CurrentWidth);
    for(int it_step=0; it_step<100; it_step++){
      if(it_seed==NSignal-1) break;
      CurrentCenter+= Mmumu_Step[it_seed];
      CurrentWidth  = (CurrentCenter-Mmumu_CentSeed[it_seed])/10.*Mmumu_WidthSeed[it_seed+1]
                     +(Mmumu_CentSeed[it_seed+1]-CurrentCenter)/10.*Mmumu_WidthSeed[it_seed];
      if(CurrentCenter<Mmumu_CentSeed[it_seed]+10.){
        Mmumu_CentFull.push_back(CurrentCenter);
        Mmumu_WidthFull.push_back(CurrentWidth);
      }
      else break;
    }
  }


  if(EMuMu){
    if( !(EleLColl.size()==1 && MuLColl.size()==2) ) return;
    if( !(EleTColl.size()==1 && MuTColl.size()==2) ) return;
    if( !(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) ) return;
    if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && MuTColl.at(1).Pt()>5) ) return;
    //if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && MuTColl.at(1).Pt()>10) ) return;
    if( fabs(EleTColl.at(0).Eta())>2.5 ) return;

    float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
    float MTW = sqrt(2)*sqrt(MET*EleTColl.at(0).Pt()-METx*EleTColl.at(0).Px()-METy*EleTColl.at(0).Py());
    float M3l=(EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M();
    //if(Mmumu<12) return;

      FillHist("CutFlow_1e2mu"+Label, 0., weight, 0., 20., 20);
      if(Label=="") FillHist("Mmumu_3lOS"+Label, Mmumu, weight, 0., 200., 200);
      if(Label==""){ if(JetColl.size()>1) FillHist("Mmumu_3lOS2j"+Label, Mmumu, weight, 0., 200., 200); }
    if(fabs(Mmumu-91.2)<10) return;
      FillHist("CutFlow_1e2mu"+Label, 1., weight, 0., 20., 20);
      FillHist("Nb_3lOSOffZ"+Label, BJetColl.size(), weight, 0., 5., 5);
      //Only for shape study
      if(Label=="") FillHist("Mmumu_3lOSOffZ"+Label, Mmumu, weight, 0., 200., 200);
      for(int it_cent=0; it_cent<Mmumu_CentFull.size(); it_cent++){
        if(Label!="") continue;
        float LargeWidth=5.;
        if((Mmumu_CentFull.at(it_cent)-LargeWidth)<12.) LargeWidth=Mmumu_CentFull.at(it_cent)-12.;
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidth){
          FillHist("Mmumu_CntBin_3lOSOffZ_Width10"+Label, it_cent+0.01, weight, 0., 95., 95);
        }
      }//----------------

    if(BJetColl.size()==0) return;
      FillHist("CutFlow_1e2mu"+Label, 2., weight, 0., 20., 20);
      FillHist("Nj_3lOSOffZHasB"+Label, JetColl.size(), weight, 0., 10., 10);
    if(JetColl.size()<2) return;
      FillHist("CutFlow_1e2mu"+Label, 3., weight, 0., 20., 20);

    //if(isData && !k_running_nonprompt && Mmumu<81.2) return;
    if(!isData && k_sample_name.Contains("TT_powheg")){
      std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
      int NFk=0;
      for(int it_lep=0; it_lep<EleTColl.size(); it_lep++){
        int LepType = GetLeptonType(EleTColl.at(it_lep),truthColl);
        if(LepType<0 && LepType>-5) NFk++;
      }
      for(int it_lep=0; it_lep<MuTColl.size(); it_lep++){
        int LepType = GetLeptonType(MuTColl.at(it_lep),truthColl);
        if(LepType<0 && LepType>-5) NFk++;
      }
      if(NFk<3) return;
    }

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
    if(isData or (k_sample_name.Contains("TTToHcToWA") and Label=="")) FillHist("Mmumu_BinWp1_SR2j"    +Label, Mmumu, weight, 0., 80., 800);
    //if(isData or (k_sample_name.Contains("TTToHcToWA") and Label=="")) FillHist("Mmumu_BinWp1_SR2j"    +Label, Mmumu, weight, 12., 80., 680);
    FillHist("Mmumu_BinW10_SR2j"   +Label, Mmumu, weight, &MmumuEdges[0], MmumuEdges.size()-1);
    FillHist("Mmumu_BinW2_SR2j"    +Label, Mmumu, weight, 12., 80., 34);
    FillHist("Mmumu_BinW2Full_SR2j"+Label, Mmumu, weight, &MmumuEdges2[0], MmumuEdges2.size()-1);
    FillHist("Mmumu_BinW1_SR2j"    +Label, Mmumu, weight, &MmumuEdges3[0], MmumuEdges3.size()-1);
    FillHist("Mmumu_BinRes_SR2j"   +Label, Mmumu, weight, &MmumuEdges4[0], MmumuEdges4.size()-1);
    for(int it_cent=0; it_cent<113; it_cent++){
      float BinCenter=(MmumuEdges4.at(it_cent)+MmumuEdges4.at(it_cent+1))/2.;
      float LargeWidth=5.;
      if     ((BinCenter-LargeWidth)<12.)  LargeWidth=BinCenter-12.;
      else if((BinCenter+LargeWidth)>81.2) LargeWidth=81.2-BinCenter;
      if(fabs(Mmumu-BinCenter)<LargeWidth){
        FillHist("Mmumu_BinRes_SR2j_Width10"+Label, BinCenter, weight, &MmumuEdges4[0], MmumuEdges4.size()-1);
      }
      if(Label==""){
        float Variation=LargeWidth>3.? 0.2:0.;
        float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
        if     ((BinCenter-LargeWidthUp)<12. ) LargeWidthUp=BinCenter-12.;
        else if((BinCenter+LargeWidthUp)>81.2) LargeWidthUp=81.2-BinCenter;
        if(fabs(Mmumu-BinCenter)<LargeWidthUp){
          FillHist("Mmumu_BinRes_SR2j_Width10_systup_Mrange", BinCenter, weight, &MmumuEdges4[0], MmumuEdges4.size()-1);
        }
        if(fabs(Mmumu-BinCenter)<LargeWidthDown){
          FillHist("Mmumu_BinRes_SR2j_Width10_systdown_Mrange", BinCenter, weight, &MmumuEdges4[0], MmumuEdges4.size()-1);
        }
      }
    }
    for(int it_cent=0; it_cent<69; it_cent++){
      float BinCenter=it_cent!=68? 12.5+it_cent:80.6;
      float LargeWidth=5.;
      if     ((BinCenter-LargeWidth)<12.)  LargeWidth=BinCenter-12.;
      else if((BinCenter+LargeWidth)>81.2) LargeWidth=81.2-BinCenter;
      if(fabs(Mmumu-BinCenter)<LargeWidth){
        FillHist("Mmumu_BinW1_SR2j_Width10"+Label, BinCenter, weight, &MmumuEdges3[0], MmumuEdges3.size()-1);
      }
      if(Label==""){
        float Variation=it_cent!=68? 0.2:0.;
        float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
        if     ((BinCenter-LargeWidthUp)<12. ) LargeWidthUp=BinCenter-12.;
        else if((BinCenter+LargeWidthUp)>81.2) LargeWidthUp=81.2-BinCenter;
        if(fabs(Mmumu-BinCenter)<LargeWidthUp){
          FillHist("Mmumu_BinW1_SR2j_Width10_systup_Mrange", BinCenter, weight, &MmumuEdges3[0], MmumuEdges3.size()-1);
        }
        if(fabs(Mmumu-BinCenter)<LargeWidthDown){
          FillHist("Mmumu_BinW1_SR2j_Width10_systdown_Mrange", BinCenter, weight, &MmumuEdges3[0], MmumuEdges3.size()-1);
        }
      }
    }
    for(int it_cent=0; it_cent<35; it_cent++){
      float BinCenter=it_cent!=34? 13.+2.*((float) it_cent):80.6;
      float LargeWidth=5.;
      if     ((BinCenter-LargeWidth)<12.)  LargeWidth=BinCenter-12.;
      else if((BinCenter+LargeWidth)>81.2) LargeWidth=81.2-BinCenter;
      if(fabs(Mmumu-BinCenter)<LargeWidth){
        FillHist("Mmumu_BinW2Full_SR2j_Width10"+Label, BinCenter, weight, &MmumuEdges2[0], MmumuEdges2.size()-1);
      }
      if(Label==""){
        float Variation=it_cent!=34? 0.2:0.;
        float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
        if     ((BinCenter-LargeWidthUp)<12. ) LargeWidthUp=BinCenter-12.;
        else if((BinCenter+LargeWidthUp)>81.2) LargeWidthUp=81.2-BinCenter;
        if(fabs(Mmumu-BinCenter)<LargeWidthUp){
          FillHist("Mmumu_BinW2Full_SR2j_Width10_systup_Mrange", BinCenter, weight, &MmumuEdges2[0], MmumuEdges2.size()-1);
        }
        if(fabs(Mmumu-BinCenter)<LargeWidthDown){
          FillHist("Mmumu_BinW2Full_SR2j_Width10_systdown_Mrange", BinCenter, weight, &MmumuEdges2[0], MmumuEdges2.size()-1);
        }
      }
    }
    for(int it_cent=0; it_cent<34; it_cent++){
      float BinCenter=13.+2.*((float) it_cent);
      float LargeWidth=5.;
      if     ((BinCenter-LargeWidth)<12.)  LargeWidth=BinCenter-12.;
      else if((BinCenter+LargeWidth)>81.2) LargeWidth=81.2-BinCenter;
      if(fabs(Mmumu-BinCenter)<LargeWidth){
        FillHist("Mmumu_BinW2_SR2j_Width10"+Label, BinCenter, weight, 12., 80., 34);
      }
      if(Label==""){
        float Variation=0.2;
        float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
        if     ((BinCenter-LargeWidthUp)<12. ) LargeWidthUp=BinCenter-12.;
        else if((BinCenter+LargeWidthUp)>81.2) LargeWidthUp=81.2-BinCenter;
        if(fabs(Mmumu-BinCenter)<LargeWidthUp){
          FillHist("Mmumu_BinW2_SR2j_Width10_systup_Mrange", BinCenter, weight, 12., 80., 34);
        }
        if((BinCenter-LargeWidthDown)<12.) LargeWidthDown=BinCenter-12.;
        if(fabs(Mmumu-BinCenter)<LargeWidthDown){
          FillHist("Mmumu_BinW2_SR2j_Width10_systdown_Mrange", BinCenter, weight, 12., 80., 34);
        }
      }
    }
    for(int it_cent=0; it_cent<Mmumu_CentFull.size(); it_cent++){
      if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<Mmumu_WidthFull.at(it_cent)){
        FillHist("Mmumu_CntBin_SR2j"+Label, it_cent+0.01, weight, 0., 95., 95);
      }
      float LargeWidth=5.;
      if((Mmumu_CentFull.at(it_cent)-LargeWidth)<12.) LargeWidth=Mmumu_CentFull.at(it_cent)-12.;
      if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidth){
        FillHist("Mmumu_CntBin_SR2j_Width10"+Label, it_cent+0.01, weight, 0., 95., 95);
      }
      if(Label==""){
        float Variation=0.2;
        float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
        if((Mmumu_CentFull.at(it_cent)-LargeWidthUp)<12.) LargeWidthUp=Mmumu_CentFull.at(it_cent)-12.;
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthUp){
          FillHist("Mmumu_CntBin_SR2j_Width10_systup_Mrange", it_cent+0.01, weight, 0., 95., 95);
        }
        if((Mmumu_CentFull.at(it_cent)-LargeWidthDown)<12.) LargeWidthDown=Mmumu_CentFull.at(it_cent)-12.;
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthDown){
          FillHist("Mmumu_CntBin_SR2j_Width10_systdown_Mrange", it_cent+0.01, weight, 0., 95., 95);
        }
      }
    }
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
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 && MuTColl.at(2).Pt()>5) ) return;
    //if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 && MuTColl.at(2).Pt()>10) ) return;

    int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    int   IdxSSW = TriMuChargeIndex(MuTColl, MET, METx, METy, "SSW");
    int   IdxSSA = TriMuChargeIndex(MuTColl, MET, METx, METy, "SSA");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    float Mmumu  = (MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)).M();

    //if( !(MOSSS1>12 && MOSSS2>12) ) return;

      FillHist("CutFlow_3mu"+Label, 0., weight, 0., 20., 20);
      if(Label=="") FillHist("Mmumu_3lOS"+Label, Mmumu, weight, 0., 200., 200);
    if( fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10 ) return;
      FillHist("CutFlow_3mu"+Label, 1., weight, 0., 20., 20);
      FillHist("Nb_3lOSOffZ"+Label, BJetColl.size(), weight, 0., 5., 5);
      //Only for shape studies
      if(Label=="") FillHist("Mmumu_3lOSOffZ"+Label, Mmumu, weight, 0., 200., 200);
      for(int it_cent=0; it_cent<Mmumu_CentFull.size(); it_cent++){
        if(Label!="") continue;
        float LargeWidth=5.;
        if((Mmumu_CentFull.at(it_cent)-LargeWidth)<12.) LargeWidth=Mmumu_CentFull.at(it_cent)-12.;
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidth){
          FillHist("Mmumu_CntBin_3lOSOffZ_Width10"+Label, it_cent+0.01, weight, 0., 95., 95);
        }
      }//----------------------

    if(BJetColl.size()==0) return;
      FillHist("CutFlow_3mu"+Label, 2., weight, 0., 20., 20);
      FillHist("Nj_3lOSOffZHasB"+Label, JetColl.size(), weight, 0., 10., 10);
    if(JetColl.size()<2) return;
      FillHist("CutFlow_3mu"+Label, 3., weight, 0., 20., 20);

    //if(isData && !k_running_nonprompt && Mmumu<81.2) return;
    if(!isData && k_sample_name.Contains("TT_powheg")){
      std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
      int NFk=0;
      for(int it_lep=0; it_lep<MuTColl.size(); it_lep++){
        int LepType = GetLeptonType(MuTColl.at(it_lep),truthColl);
        if(LepType<0 && LepType>-5) NFk++;
      }
      if(NFk<3) return;
    }


    int Idxj1_W=-1, Idxj2_W=-1;
    snu::KParticle v; v.SetPxPyPzE(METx, METy, 0, sqrt(METx*METx+METy*METy));
    float Pzv1=GetvPz(v, MuTColl.at(IdxSSW), 1), Pzv2=GetvPz(v, MuTColl.at(IdxSSW), 2);
    float Pzv=fabs(Pzv1)<fabs(Pzv2)? Pzv1:Pzv2;
    v.SetXYZM(METx,METy,Pzv,0.);

    float M3lv=(MuTColl.at(0)+MuTColl.at(1)+MuTColl.at(2)+v).M(), M2l2j=0.; 
    Idxj1_W=GetDaughterCandIdx(JetColl, "W", -1., "Lead");
    Idxj2_W=GetDaughterCandIdx(JetColl, "W", -1., "Subl");
    M2l2j=(MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)+JetColl.at(Idxj1_W)+JetColl.at(Idxj2_W)).M();

    FillHist("Nj_SR2j"+Label, JetColl.size(), weight, 0., 10., 10);
    FillHist("Nb_SR2j"+Label, BJetColl.size(), weight, 0., 5., 5);
    FillHist("PTmu1_SR2j"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
    FillHist("PTmu2_SR2j"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
    FillHist("PTmu3_SR2j"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
    FillHist("MOSSS1_SR2j"+Label, MOSSS1, weight, 0., 200., 80);
    if(isData or (k_sample_name.Contains("TTToHcToWA") and Label=="")) FillHist("Mmumu_BinWp1_SR2j"    +Label, Mmumu, weight, 0., 80., 800);
    //if(isData or (k_sample_name.Contains("TTToHcToWA") and Label=="")) FillHist("Mmumu_BinWp1_SR2j"    +Label, Mmumu, weight, 12., 80., 680);
    FillHist("Mmumu_BinW10_SR2j"+Label, Mmumu, weight, &MmumuEdges[0], MmumuEdges.size()-1);
    FillHist("Mmumu_BinW2_SR2j" +Label, Mmumu, weight, 12., 80., 34);
    FillHist("Mmumu_BinW2Full_SR2j"+Label, Mmumu, weight, &MmumuEdges2[0], MmumuEdges2.size()-1);
    FillHist("Mmumu_BinW1_SR2j"    +Label, Mmumu, weight, &MmumuEdges3[0], MmumuEdges3.size()-1);
    FillHist("Mmumu_BinRes_SR2j"   +Label, Mmumu, weight, &MmumuEdges4[0], MmumuEdges4.size()-1);
    for(int it_cent=0; it_cent<113; it_cent++){
      float BinCenter=(MmumuEdges4.at(it_cent)+MmumuEdges4.at(it_cent+1))/2.;
      float LargeWidth=5.;
      if     ((BinCenter-LargeWidth)<12.)  LargeWidth=BinCenter-12.;
      else if((BinCenter+LargeWidth)>81.2) LargeWidth=81.2-BinCenter;
      if(fabs(Mmumu-BinCenter)<LargeWidth){
        FillHist("Mmumu_BinRes_SR2j_Width10"+Label, BinCenter, weight, &MmumuEdges4[0], MmumuEdges4.size()-1);
      }
      if(Label==""){
        float Variation=LargeWidth>3.? 0.2:0.;
        float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
        if     ((BinCenter-LargeWidthUp)<12. ) LargeWidthUp=BinCenter-12.;
        else if((BinCenter+LargeWidthUp)>81.2) LargeWidthUp=81.2-BinCenter;
        if(fabs(Mmumu-BinCenter)<LargeWidthUp){
          FillHist("Mmumu_BinRes_SR2j_Width10_systup_Mrange", BinCenter, weight, &MmumuEdges4[0], MmumuEdges4.size()-1);
        }
        if(fabs(Mmumu-BinCenter)<LargeWidthDown){
          FillHist("Mmumu_BinRes_SR2j_Width10_systdown_Mrange", BinCenter, weight, &MmumuEdges4[0], MmumuEdges4.size()-1);
        }
      }
    }
    for(int it_cent=0; it_cent<69; it_cent++){
      float BinCenter=it_cent!=68? 12.5+it_cent:80.6;
      float LargeWidth=5.;
      if     ((BinCenter-LargeWidth)<12.)  LargeWidth=BinCenter-12.;
      else if((BinCenter+LargeWidth)>81.2) LargeWidth=81.2-BinCenter;
      if(fabs(Mmumu-BinCenter)<LargeWidth){
        FillHist("Mmumu_BinW1_SR2j_Width10"+Label, BinCenter, weight, &MmumuEdges3[0], MmumuEdges3.size()-1);
      }
      if(Label==""){
        float Variation=it_cent!=68? 0.2:0.;
        float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
        if     ((BinCenter-LargeWidthUp)<12. ) LargeWidthUp=BinCenter-12.;
        else if((BinCenter+LargeWidthUp)>81.2) LargeWidthUp=81.2-BinCenter;
        if(fabs(Mmumu-BinCenter)<LargeWidthUp){
          FillHist("Mmumu_BinW1_SR2j_Width10_systup_Mrange", BinCenter, weight, &MmumuEdges3[0], MmumuEdges3.size()-1);
        }
        if(fabs(Mmumu-BinCenter)<LargeWidthDown){
          FillHist("Mmumu_BinW1_SR2j_Width10_systdown_Mrange", BinCenter, weight, &MmumuEdges3[0], MmumuEdges3.size()-1);
        }
      }
    }
    for(int it_cent=0; it_cent<35; it_cent++){
      float BinCenter=it_cent!=34? 13.+2.*((float) it_cent):80.6;
      float LargeWidth=5.;
      if     ((BinCenter-LargeWidth)<12.)  LargeWidth=BinCenter-12.;
      else if((BinCenter+LargeWidth)>81.2) LargeWidth=81.2-BinCenter;
      if(fabs(Mmumu-BinCenter)<LargeWidth){
        FillHist("Mmumu_BinW2Full_SR2j_Width10"+Label, BinCenter, weight, &MmumuEdges2[0], MmumuEdges2.size()-1);
      }
      if(Label==""){
        float Variation=it_cent!=34? 0.2:0.;
        float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
        if     ((BinCenter-LargeWidthUp)<12. ) LargeWidthUp=BinCenter-12.;
        else if((BinCenter+LargeWidthUp)>81.2) LargeWidthUp=81.2-BinCenter;
        if(fabs(Mmumu-BinCenter)<LargeWidthUp){
          FillHist("Mmumu_BinW2Full_SR2j_Width10_systup_Mrange", BinCenter, weight, &MmumuEdges2[0], MmumuEdges2.size()-1);
        }
        if(fabs(Mmumu-BinCenter)<LargeWidthDown){
          FillHist("Mmumu_BinW2Full_SR2j_Width10_systdown_Mrange", BinCenter, weight, &MmumuEdges2[0], MmumuEdges2.size()-1);
        }
      }
    }
    for(int it_cent=0; it_cent<34; it_cent++){
      float BinCenter=13.+2.*((float) it_cent);
      float LargeWidth=5.;
      if     ((BinCenter-LargeWidth)<12. ) LargeWidth=BinCenter-12.;
      else if((BinCenter+LargeWidth)>81.2) LargeWidth=81.2-BinCenter;
      if(fabs(Mmumu-BinCenter)<LargeWidth){
        FillHist("Mmumu_BinW2_SR2j_Width10"+Label, BinCenter, weight, 12., 80., 34);
      }
      if(Label==""){
        float Variation=0.2;
        float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
        if     ((BinCenter-LargeWidthUp)<12. ) LargeWidthUp=BinCenter-12.;
        else if((BinCenter+LargeWidthUp)>81.2) LargeWidthUp=81.2-BinCenter;
        if(fabs(Mmumu-BinCenter)<LargeWidthUp){
          FillHist("Mmumu_BinW2_SR2j_Width10_systup_Mrange", BinCenter, weight, 12., 80., 34);
        }
        if((BinCenter-LargeWidthDown)<12.) LargeWidthDown=BinCenter-12.;
        if(fabs(Mmumu-BinCenter)<LargeWidthDown){
          FillHist("Mmumu_BinW2_SR2j_Width10_systdown_Mrange", BinCenter, weight, 12., 80., 34);
        }
      }
    }
    for(int it_cent=0; it_cent<Mmumu_CentFull.size(); it_cent++){
      if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<Mmumu_WidthFull.at(it_cent)){
        FillHist("Mmumu_CntBin_SR2j"+Label, it_cent+0.01, weight, 0., 95., 95);
      }
      float LargeWidth=5.;
      if((Mmumu_CentFull.at(it_cent)-LargeWidth)<12.) LargeWidth=Mmumu_CentFull.at(it_cent)-12.;
      if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidth){
        FillHist("Mmumu_CntBin_SR2j_Width10"+Label, it_cent+0.01, weight, 0., 95., 95);
      }
      if(Label==""){
        float Variation=0.2;
        float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
        if((Mmumu_CentFull.at(it_cent)-LargeWidthUp)<12.) LargeWidthUp=Mmumu_CentFull.at(it_cent)-12.;
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthUp){
          FillHist("Mmumu_CntBin_SR2j_Width10_systup_Mrange", it_cent+0.01, weight, 0., 95., 95);
        }
        if((Mmumu_CentFull.at(it_cent)-LargeWidthDown)<12.) LargeWidthDown=Mmumu_CentFull.at(it_cent)-12.;
        if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthDown){
          FillHist("Mmumu_CntBin_SR2j_Width10_systdown_Mrange", it_cent+0.01, weight, 0., 95., 95);
        }
      }
    }
    FillHist("Mmumu_SR2j"+Label, Mmumu, weight, 0., 200., 200);
    FillHist("M2l2j_SR2j"+Label, M2l2j, weight, 0., 1000., 100);
    FillHist("M3lv_SR2j"+Label, M3lv, weight, 0., 1000., 100);
    if(Mmumu<80){
      FillHist("M2l2j_lt80_SR2j"+Label, M2l2j, weight, 0., 1000., 100);
      FillHist("M3lv_lt80_SR2j"+Label, M3lv, weight, 0., 1000., 100);
    }
  }
}




void Mar2019_AdditionalChecks::DoSystRun(TString Cycle, TString Mode, std::vector<snu::KElectron> EleColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight){

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



float Mar2019_AdditionalChecks::ConeCorrectedPT(snu::KMuon Mu, float TightIsoCut){

  //float PTCorr=Mu.Pt();
  //float PTCorr=Mu.Pt()*(1.+RochIso(Mu,"0.4"));
  float PTCorr=Mu.Pt()*(1.+max((float) 0.,RochIso(Mu,"0.4")-TightIsoCut));
  return PTCorr;

}



float Mar2019_AdditionalChecks::ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut){

  //float PTCorr=Ele.Pt()*(1+Ele.PFRelIso(0.3));
  float PTCorr=Ele.Pt()*(1.+max(0.,Ele.PFRelIso(0.3)-TightIsoCut));
  return PTCorr;

}



void Mar2019_AdditionalChecks::FillCutFlow(TString cut, float weight){
  
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



void Mar2019_AdditionalChecks::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Mar2019_AdditionalChecks::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //AnalyzerCore::MakeHistograms("Basic_Nj_orig", 10, 0., 10.);
  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Mar2019_AdditionalChecksCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Mar2019_AdditionalChecks::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


int Mar2019_AdditionalChecks::GetFakeLepSrcJetType(snu::KElectron Ele, std::vector<snu::KJet> JetColl){
  //Type: -1: Unmatched, 1:L, 2:C, 3:B
  int SrcType=-1;
  bool NearB=false, NearC=false, NearL=false;
  for(int i=0; i<(int) JetColl.size(); i++){
    if(Ele.DeltaR(JetColl.at(i))<0.4){
      if(JetColl.at(i).HadronFlavour()==5){ NearB=true; break; }//1)
      else if(JetColl.at(i).HadronFlavour()==4){ NearC=true; }
      else if(JetColl.at(i).HadronFlavour()==0){ NearL=true; }
    }
  }

  if     (NearB) SrcType=3;
  else if(NearC) SrcType=2;
  else if(NearL) SrcType=1;

  return SrcType;
//1) Higher Priority to B. if there's multiple near jets, then b-jet has higher priority
}


int Mar2019_AdditionalChecks::GetFakeLepSrcJetType(snu::KMuon Mu, std::vector<snu::KJet> JetColl){
  //Type: -1: Unmatched, 1:L, 2:C, 3:B
  int SrcType=-1;
  bool NearB=false, NearC=false, NearL=false;
  for(int i=0; i<(int) JetColl.size(); i++){
    if(Mu.DeltaR(JetColl.at(i))<0.4){
      if(JetColl.at(i).HadronFlavour()==5){ NearB=true; break; }//1)
      else if(JetColl.at(i).HadronFlavour()==4){ NearC=true; }
      else if(JetColl.at(i).HadronFlavour()==0){ NearL=true; }
    }
  }

  if     (NearB) SrcType=3;
  else if(NearC) SrcType=2;
  else if(NearL) SrcType=1;

  return SrcType;
//1) Higher Priority to B. if there's multiple near jets, then b-jet has higher priority
}
