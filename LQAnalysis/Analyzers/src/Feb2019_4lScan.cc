// $Id: Feb2019_4lScan.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQFeb2019_4lScan Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Feb2019_4lScan.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Feb2019_4lScan);

 Feb2019_4lScan::Feb2019_4lScan() : AnalyzerCore(), out_muons(0) {

   SetLogName("Feb2019_4lScan");
   Message("In Feb2019_4lScan constructor", INFO);
   InitialiseAnalysis();
 }


 void Feb2019_4lScan::InitialiseAnalysis() throw( LQError ) {
   
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

void Feb2019_4lScan::ExecuteEvents()throw( LQError ){

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


   bool DoubleMuon=false, DoubleEG=false, MuonEG=false;
   bool HighMassZZ=false;
   bool SystRun=false;
   TString Cycle="";
   for(int i=0; i<(int) k_flags.size(); i++){
     if     (k_flags.at(i).Contains("SystRun"))       SystRun       = true;
     else if(k_flags.at(i).Contains("DoubleMuon"))    DoubleMuon    = true;
     else if(k_flags.at(i).Contains("MuonEG"))        MuonEG        = true;
     else if(k_flags.at(i).Contains("DoubleEG"))      DoubleEG      = true;
     else if(k_flags.at(i).Contains("HighMassZZ"))    HighMassZZ    = true;
   }



    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //Trigger Path of Analysis
   bool Pass_Trigger=false, Pass_TriggerBG=false, Pass_TriggerH=false;
   float trigger_ps_weight=1., trigger_period_weight=1.;
   float LumiBG=27.257618, LumiH=8.605696, LumiBH=35.863314;

   if(MuonEG){
     if(  PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") 
       or PassTrigger("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v") )  Pass_TriggerBG=true;
     if(  PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v")
       or PassTrigger("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_TriggerH =true;

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
   if(DoubleMuon){
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
   if(DoubleEG){
     if( PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger=true;
     if(!isData){trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);}
   }
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
     if(DoubleMuon)   { if( !(muonPreColl.size()>=4) ) return; }
     else if(MuonEG)  { if( !(muonPreColl.size()>=2 && electronPreColl.size()>=2) ) return; }
     else if(DoubleEG){ if( !(electronPreColl.size()>=4) ) return; }
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
//     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.4, 0.4);
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
//     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.15, 0.15);
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
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

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
   std::vector<snu::KJet> jetNoVeto10Coll; eventbase->GetJetSel()->Selection(jetNoVeto10Coll, LeptonVeto, muonLooseColl, electronLooseColl);


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
   if     (MuonEG)    { if(electronColl.size()>=2 && muonColl.size()>=2) EventCand=true; }
   else if(DoubleMuon){ if(               muonColl.size()>=4           ) EventCand=true; }
   else if(DoubleEG)  { if(             electronColl.size()>=4         ) EventCand=true; }

   if(EventCand & !SystRun){
     if(!isData){
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ID_ELECTRON_MVA_90", electronColl);
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
    
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl);

         btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);
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

 
   if(HighMassZZ){
     if     (DoubleMuon){
       CheckMZZSpectrum(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                        weight, "", "DiMu");
     }
     else if(MuonEG){
       CheckMZZSpectrum(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                        weight, "", "EMu");
     }
     else if(DoubleEG){
       CheckMZZSpectrum(muonColl, muonLooseColl, electronColl, electronLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                        weight, "", "DiEl");
     }
   }

return;

}// End of execute event loop
  


void Feb2019_4lScan::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Feb2019_4lScan::BeginCycle() throw( LQError ){
  
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

Feb2019_4lScan::~Feb2019_4lScan() {
  
  Message("In Feb2019_4lScan Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}


void Feb2019_4lScan::CheckMZZSpectrum(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option){


  bool DoubleMuon=Option.Contains("DiMu"), MuonEG=Option.Contains("EMu"), DoubleEG=Option.Contains("DiEl");

  float M4l=0., Nj=0., Nb=0., Mjj=0.;
  bool  OnShellZZ =false;

  if(DoubleMuon){
    if( !(MuTColl.size()==4) ) return;
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(3).Pt()>10) ) return;
    int Q = 0; for(int it_m=0; it_m<MuTColl.size(); it_m++){ Q+=MuTColl.at(it_m).Charge(); }
    if(Q!=0) return;

    int NZOnShell=0;
    for(unsigned int i=0; i<MuTColl.size(); i++){
      for(unsigned int j=i+1; j<MuTColl.size(); j++){
        if(MuTColl.at(i).Charge()==MuTColl.at(j).Charge()) continue;
        if(fabs( (MuTColl.at(i)+MuTColl.at(j)).M()-91.2 )<15) NZOnShell++;
        if(      (MuTColl.at(i)+MuTColl.at(j)).M()<4.       ) return;
      }
    }
    if(NZOnShell>1) OnShellZZ=true;

    Nj = JetColl.size();
    Nb = BJetColl.size(); 
    Mjj = Nj>1.5? (JetColl.at(0)+JetColl.at(1)).M():-1.;
    M4l = (MuTColl.at(0)+MuTColl.at(1)+MuTColl.at(2)+MuTColl.at(3)).M();

    FillHist("M4l_incl", M4l, weight, 0., 500., 500);
    FillHist("Nj_incl",  Nj, weight, 0., 10., 10);
    FillHist("Nb_incl",  Nb, weight, 0., 5., 5);
    if(Nj>1.5) FillHist("Mjj_incl", Mjj, weight, 0., 500., 50);
    if(OnShellZZ){
      FillHist("M4l_ZZ", M4l, weight, 0., 500., 500);
      FillHist("Nj_ZZ",  Nj, weight, 0., 10., 10);
      FillHist("Nb_ZZ",  Nb, weight, 0., 5., 5);
      if(Nj>1.5){
        FillHist("Mjj_ZZ", Mjj, weight, 0., 500., 50);
        if(Nb>=1) FillHist("Mjj_ZZbj", Mjj, weight, 0., 500., 50);
      }
    }
  }
  if(MuonEG){
    if( !(MuTColl.size()==2 && EleTColl.size()==2) ) return;
    if( !(MuTColl.at(1).Pt()>10 && EleTColl.at(1).Pt()>10) ) return;
    if( !(MuTColl.at(0).Pt()>25 || EleTColl.at(0).Pt()>25) ) return;
    if(   fabs( ( MuTColl.at(0)+ MuTColl.at(1)).M()-91.2 )<15 
       && fabs( (EleTColl.at(0)+EleTColl.at(1)).M()-91.2 )<15 ) OnShellZZ=true;
    if(   ( MuTColl.at(0)+ MuTColl.at(1)).M()<4. 
       || (EleTColl.at(0)+EleTColl.at(1)).M()<4. ) return;
    if( MuTColl.at(0).Charge()==MuTColl.at(1).Charge() || EleTColl.at(0).Charge()==EleTColl.at(1).Charge() ) return;


    Nj = JetColl.size();
    Nb = BJetColl.size(); 
    Mjj = Nj>1.5? (JetColl.at(0)+JetColl.at(1)).M():-1.;
    M4l = (MuTColl.at(0)+MuTColl.at(1)+EleTColl.at(0)+EleTColl.at(1)).M();

    FillHist("M4l_incl", M4l, weight, 0., 500., 500);
    FillHist("Nj_incl",  Nj, weight, 0., 10., 10);
    FillHist("Nb_incl",  Nb, weight, 0., 5., 5);
    if(Nj>1.5) FillHist("Mjj_incl", Mjj, weight, 0., 500., 50);
    if(OnShellZZ){
      FillHist("M4l_ZZ", M4l, weight, 0., 500., 500);
      FillHist("Nj_ZZ",  Nj, weight, 0., 10., 10);
      FillHist("Nb_ZZ",  Nb, weight, 0., 5., 5);
      if(Nj>1.5){
        FillHist("Mjj_ZZ", Mjj, weight, 0., 500., 50);
        if(Nb>=1) FillHist("Mjj_ZZbj", Mjj, weight, 0., 500., 50);
      }
    }
  }
  if(DoubleEG){
    if( !(EleTColl.size()==4) ) return;
    if( !(EleTColl.at(0).Pt()>25 && EleTColl.at(1).Pt()>15 && EleTColl.at(3).Pt()>10) ) return;
    int Q = 0; for(int it_el=0; it_el<EleTColl.size(); it_el++){ Q+=EleTColl.at(it_el).Charge(); }
    if(Q!=0) return;

    int NZOnShell=0;
    for(unsigned int i=0; i<EleTColl.size(); i++){
      for(unsigned int j=i+1; j<EleTColl.size(); j++){
        if(EleTColl.at(i).Charge()==EleTColl.at(j).Charge()) continue;
        if(fabs( (EleTColl.at(i)+EleTColl.at(j)).M()-91.2 )<15) NZOnShell++;
        if(      (EleTColl.at(i)+EleTColl.at(j)).M()<4.       ) return;
      }
    }
    if(NZOnShell>1) OnShellZZ=true;

    Nj = JetColl.size();
    Nb = BJetColl.size(); 
    Mjj = Nj>1.5? (JetColl.at(0)+JetColl.at(1)).M():-1.;
    M4l = (EleTColl.at(0)+EleTColl.at(1)+EleTColl.at(2)+EleTColl.at(3)).M();

    FillHist("M4l_incl", M4l, weight, 0., 500., 500);
    FillHist("Nj_incl",  Nj, weight, 0., 10., 10);
    FillHist("Nb_incl",  Nb, weight, 0., 5., 5);
    if(Nj>1.5) FillHist("Mjj_incl", Mjj, weight, 0., 500., 50);
    if(OnShellZZ){
      FillHist("M4l_ZZ", M4l, weight, 0., 500., 500);
      FillHist("Nj_ZZ",  Nj, weight, 0., 10., 10);
      FillHist("Nb_ZZ",  Nb, weight, 0., 5., 5);
      if(Nj>1.5){
        FillHist("Mjj_ZZ", Mjj, weight, 0., 500., 50);
        if(Nb>=1) FillHist("Mjj_ZZbj", Mjj, weight, 0., 500., 50);
      }
    }
  }

}




void Feb2019_4lScan::DoSystRun(TString Cycle, TString Mode, std::vector<snu::KElectron> EleColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight){

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

  return;
}



float Feb2019_4lScan::ConeCorrectedPT(snu::KMuon Mu, float TightIsoCut){

  //float PTCorr=Mu.Pt();
  //float PTCorr=Mu.Pt()*(1.+RochIso(Mu,"0.4"));
  float PTCorr=Mu.Pt()*(1.+max((float) 0.,RochIso(Mu,"0.4")-TightIsoCut));
  return PTCorr;

}



float Feb2019_4lScan::ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut){

  //float PTCorr=Ele.Pt()*(1+Ele.PFRelIso(0.3));
  float PTCorr=Ele.Pt()*(1.+max(0.,Ele.PFRelIso(0.3)-TightIsoCut));
  return PTCorr;

}



void Feb2019_4lScan::FillCutFlow(TString cut, float weight){
  
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



void Feb2019_4lScan::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Feb2019_4lScan::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //AnalyzerCore::MakeHistograms("Basic_Nj_orig", 10, 0., 10.);
  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Feb2019_4lScanCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Feb2019_4lScan::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
