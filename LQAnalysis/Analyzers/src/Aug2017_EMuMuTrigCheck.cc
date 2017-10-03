// $Id: Aug2017_EMuMuTrigCheck.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQAug2017_EMuMuTrigCheck Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Aug2017_EMuMuTrigCheck.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Aug2017_EMuMuTrigCheck);

 Aug2017_EMuMuTrigCheck::Aug2017_EMuMuTrigCheck() : AnalyzerCore(), out_muons(0) {

   SetLogName("Aug2017_EMuMuTrigCheck");
   Message("In Aug2017_EMuMuTrigCheck constructor", INFO);
   InitialiseAnalysis();
 }


 void Aug2017_EMuMuTrigCheck::InitialiseAnalysis() throw( LQError ) {
   
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

void Aug2017_EMuMuTrigCheck::ExecuteEvents()throw( LQError ){

////Basic Infos///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Apply the gen weight 
   if(!isData) weight*=MCweight;
   FillHist("GenWeight", MCweight, 1, -2., 2., 4);


   m_logger<<DEBUG<<"RunNumber/Event Number = "<<eventbase->GetEvent().RunNumber()<<" : "<<eventbase->GetEvent().EventNumber()<<LQLogger::endmsg;
   m_logger<<DEBUG<<"isData = "<<isData<<LQLogger::endmsg;


   //Pileup Reweight
   float pileup_reweight=1., pileup_reweight_systdown=1., pileup_reweight_systup=1.;
   if(!k_isdata){ pileup_reweight          = eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central); 
                  pileup_reweight_systup   = pileup_reweight>0.? eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::up)  /pileup_reweight :0. ;
                  pileup_reweight_systdown = pileup_reweight>0.? eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::down)/pileup_reweight :0. ;
   }
   FillHist("Basic_PURW", pileup_reweight, 1., 0., 20., 200);

 
   //Total Event(MCweight+PUrw)////////////////
   FillCutFlow("NoCut", weight*pileup_reweight);


   bool EMuMu=false, EMu=false;
   bool Mu12Ele23=false, Mu8Ele23=false;
   bool SystRun=false;
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("EMuMu"))     EMuMu     = true;
     else if(k_flags.at(i).Contains("EMu"))       EMu       = true;
     else if(k_flags.at(i).Contains("SystRun"))   SystRun   = true;
     else if(k_flags.at(i).Contains("Mu12Ele23")) Mu12Ele23 = true;
     else if(k_flags.at(i).Contains("Mu8Ele23"))  Mu8Ele23  = true;
   }

    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //Trigger Path of Analysis
   bool Pass_Trigger=false;
   float trigger_ps_weight=1., trigger_period_weight=1.;
   if(EMuMu || EMu){
     int Pass_Trigger1=0, Pass_Trigger2=0, Pass_Trigger3=0;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ) Pass_Trigger1++;
     if( PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger2++;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v")  ) Pass_Trigger3++;

     if(isData){
       int DataPeriod=GetDataPeriod();
       if     ( DataPeriod>0 && DataPeriod<7 && Pass_Trigger1==1 ) Pass_Trigger=true;
       else if( DataPeriod==7 && Mu12Ele23 && Pass_Trigger2==1 ) Pass_Trigger=true;
       else if( DataPeriod==7 && Mu8Ele23  && Pass_Trigger3==1 ) Pass_Trigger=true;
     }
     else{
       if( Pass_Trigger1>0 || (Mu12Ele23 && Pass_Trigger2>0) || (Mu8Ele23 && Pass_Trigger3>0) ) Pass_Trigger=true;
       if( Mu12Ele23 ) trigger_period_weight=(Pass_Trigger1*27.257618+Pass_Trigger2*8.605696)/35.863314;
       else if( Mu8Ele23 ) trigger_period_weight=(Pass_Trigger1*27.257618+Pass_Trigger3*8.605696)/35.863314;
       trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
     }
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
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, true);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
     if     (EMuMu) { if( !(electronPreColl.size()>=1 && muonPreColl.size()>=2) ) return; }
     else if(EMu  ) { if( !(electronPreColl.size()>=1 && muonPreColl.size()>=1) ) return; }
     FillCutFlow("PreSel", weight);
   //******************************************************************************************************//

   //Primary Object Collection
   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);

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
   std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,true);
   std::vector<snu::KMuon> muonColl;  if(k_running_nonprompt){ muonColl=muonLooseColl;} else{ muonColl=muonTightColl;}


     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_HctoWA_FAKELOOSE);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.4, 0.4);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);

     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronTightColl; eventbase->GetElectronSel()->Selection(electronTightColl);
   std::vector<snu::KElectron> electronColl;  if(!k_running_nonprompt){ electronColl=electronTightColl;} else{ electronColl=electronLooseColl;}

     bool LeptonVeto=true;
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonLooseColl, electronLooseColl);

   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");
   std::vector<snu::KJet> ljetColl = SelLightJets(jetColl, "Medium");

   float met    = eventbase->GetEvent().MET(snu::KEvent::pfmet);
   float metphi = eventbase->GetEvent().METPhi();
   float met_x  = eventbase->GetEvent().PFMETx();
   float met_y  = eventbase->GetEvent().PFMETy();
   float Pzv, Pzv1, Pzv2;
   snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));
   //snu::KParticle v[4]; v[0].SetPx(met_x); v[0].SetPy(met_y);
   //                     v[1].SetPx(met_x); v[1].SetPy(met_y);

   int nbjets=bjetColl.size(); int njets=jetColl.size(); const int nljets=ljetColl.size();
   int Nvtx=eventbase->GetEvent().nVertices();
   //------------------------------------------------------------------------------------------------------------------//
  

   //=====================================================//
   // Correction Factors                                  //
   //=====================================================//
   float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1., btag_sf=1.;
   float trigger_sf=1.;
   float fake_weight=1.; bool EventCand=false;
   float k_factor_weight=1., geneff_weight=1., gennorm_weight=1.;

   /*This part is for boosting up speed.. SF part takes rather longer time than expected*/
   if     (EMuMu){ if(electronLooseColl.size()>=1 && muonLooseColl.size()>=2) EventCand=true; }
   else if(EMu  ){ if(electronLooseColl.size()>=1 && muonLooseColl.size()>=1) EventCand=true; }


   if(EventCand & !SystRun){
     if(!isData){
       if(EMuMu || EMu){
         k_factor_weight = GetKFactor();
         geneff_weight   = GenFilterEfficiency(k_sample_name);
         gennorm_weight  = SignalNorm(k_sample_name, 200.);
 
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electronColl);
    
         trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muonColl);
         //iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);

         btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);
  
         //trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );
       }
     }
     else{
       //Perfect Prompt Ratio Approximation applied.(p=1), Applicable to generic number, combination of leptons under premise of high prompt rate.
       if(k_running_nonprompt){
         fake_weight=GetFakeWeight(muonLooseColl, electronLooseColl, "HNTrilepFakeL2", "HNTrilepTight2", "HctoWAFakeLoose","POGMVAMIP");
       }
     }
   }
   weight *= k_factor_weight*geneff_weight*gennorm_weight;
   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*btag_sf*pileup_reweight*fake_weight;
   //-----------------------------------------------------------------------------------------//



   //----------------------------------------------------------------------------------//
   //==================================================================================//
   /////// MAIN ANALYSIS CODE ///////////////////////////////////////////////////////////
   //==================================================================================//
   //----------------------------------------------------------------------------------//


   if(EMuMu){

     bool SystRun=false;
     if(!SystRun){
       if( !(electronLooseColl.size()==1 && muonLooseColl.size()==2) ) return;
       if( !(electronColl.size()==1      && muonColl.size()==2) ) return;
       if( !(muonColl.at(0).Charge()!=muonColl.at(1).Charge()) ) return;
       if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>15 && muonColl.at(1).Pt()>10) ) return;
       if( fabs(electronColl.at(0).Eta())>2.5 ) return;
  
       float Mmumu=(muonColl.at(0)+muonColl.at(1)).M();
       float MTW = sqrt(2)*sqrt(met*electronColl.at(0).Pt()-met_x*electronColl.at(0).Px()-met_y*electronColl.at(0).Py());
       float M3l=(electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M();
       if(Mmumu<12) return;
 
 
       bool HasBJet=false, OnZ=false, OnZG=false, WZSel=false, ZGSel=false, ttZSel=false, ZbbSel=false, AN_Sideband=false;
       bool SRincl=false, SRoffZ=false;
       if( bjetColl.size()!=0                )       HasBJet =true;
       if( fabs(Mmumu-91.2)<10               )       OnZ     =true;
       if( fabs(M3l-91.2)<10                 )       OnZG    =true;
       if( OnZ && M3l>101.2 && met>50        )       WZSel   =true;
       if( OnZG && Mmumu<81.2 && met<50      )       ZGSel   =true;
       if( OnZ && HasBJet && jetColl.size()>2)       ttZSel  =true;
       if( OnZ && HasBJet && jetColl.size()<3)       ZbbSel  =true;
       if( Mmumu>40 && HasBJet && jetColl.size()>2 ) AN_Sideband=true;
       if( HasBJet && jetColl.size()>2         )     SRincl  =true;
       if( HasBJet && jetColl.size()>2 && !OnZ )     SRoffZ  =true;
 
       FillHist("YieldComp_CR", 0., weight, 0., 10., 10);
       
  
       //General 3l Selection
       if(!HasBJet){
         FillHist("YieldComp_CR", 1., weight, 0., 10., 10);
         FillHist("PTe_3lOS", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_3lOS", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_3lOS", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_3lOS", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_3lOS", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_3lOS", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_3lOS", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_3lOS", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_3lOS", jetColl.size(), weight, 0., 10., 10);
         FillHist("Nb_3lOS", bjetColl.size(), weight, 0., 10., 10);
         FillHist("MET_3lOS", met, weight, 0., 200., 200);
         FillHist("MTW_3lOS", MTW, weight, 0., 200., 200);
       }
       //3l OnZ ; Fake CR
       if(!HasBJet && OnZ){

         FillHist("YieldComp_CR", 2., weight, 0., 10., 10);
         FillHist("PTe_3lOSOnZ", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_3lOSOnZ", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_3lOSOnZ", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_3lOSOnZ", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_3lOSOnZ", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_3lOSOnZ", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_3lOSOnZ", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_3lOSOnZ", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_3lOSOnZ", jetColl.size(), weight, 0., 10., 10);
         FillHist("Nb_3lOSOnZ", bjetColl.size(), weight, 0., 10., 10);
         FillHist("MET_3lOSOnZ", met, weight, 0., 200., 200);
         FillHist("MTW_3lOSOnZ", MTW, weight, 0., 200., 200);
       }
       //3l OnZ & HasBJet
       if(HasBJet && OnZ){
         FillHist("YieldComp_CR", 3., weight, 0., 10., 10);
         FillHist("PTe_3lOSOnZHasB", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_3lOSOnZHasB", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_3lOSOnZHasB", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_3lOSOnZHasB", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_3lOSOnZHasB", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_3lOSOnZHasB", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_3lOSOnZHasB", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_3lOSOnZHasB", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_3lOSOnZHasB", jetColl.size(), weight, 0., 10., 10);
         FillHist("Nb_3lOSOnZHasB", bjetColl.size(), weight, 0., 10., 10);
         FillHist("MET_3lOSOnZHasB", met, weight, 0., 200., 200);
         FillHist("MTW_3lOSOnZHasB", MTW, weight, 0., 200., 200);
       }
       if(WZSel){
         FillHist("YieldComp_CR", 4., weight, 0., 10., 10);
         FillHist("PTe_WZSel", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_WZSel", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_WZSel", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_WZSel", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_WZSel", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_WZSel", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_WZSel", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_WZSel", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_WZSel", jetColl.size(), weight, 0., 10., 10);
         FillHist("Nb_WZSel", bjetColl.size(), weight, 0., 10., 10);
         FillHist("MET_WZSel", met, weight, 0., 200., 200);
         FillHist("MTW_WZSel", MTW, weight, 0., 200., 200);
       }
       if(ZGSel){
         FillHist("YieldComp_CR", 5., weight, 0., 10., 10);
         FillHist("PTe_ZGSel", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_ZGSel", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_ZGSel", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_ZGSel", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_ZGSel", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_ZGSel", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_ZGSel", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_ZGSel", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_ZGSel", jetColl.size(), weight, 0., 10., 10);
         FillHist("Nb_ZGSel", bjetColl.size(), weight, 0., 10., 10);
         FillHist("MET_ZGSel", met, weight, 0., 200., 200);
         FillHist("MTW_ZGSel", MTW, weight, 0., 200., 200);
       }
       if(ttZSel){
         FillHist("YieldComp_CR", 6., weight, 0., 10., 10);
         FillHist("PTe_ttZSel", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_ttZSel", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_ttZSel", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_ttZSel", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_ttZSel", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_ttZSel", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_ttZSel", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_ttZSel", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_ttZSel", jetColl.size(), weight, 0., 10., 10);
         FillHist("Nb_ttZSel", bjetColl.size(), weight, 0., 10., 10);
         FillHist("MET_ttZSel", met, weight, 0., 200., 200);
         FillHist("MTW_ttZSel", MTW, weight, 0., 200., 200);
       }
       if(AN_Sideband){
         FillHist("YieldComp_CR", 7., weight, 0., 10., 10);
         FillHist("PTe_ANSideband", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_ANSideband", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_ANSideband", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_ANSideband", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_ANSideband", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_ANSideband", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_ANSideband", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_ANSideband", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_ANSideband", jetColl.size(), weight, 0., 10., 10);
         FillHist("Nb_ANSideband", bjetColl.size(), weight, 0., 10., 10);
         FillHist("MET_ANSideband", met, weight, 0., 200., 200);
         FillHist("MTW_ANSideband", MTW, weight, 0., 200., 200);
       }
     }//End of Not SystRun
   }//End of Closure
   if(EMu){
     if( !(electronColl.size()==1 && muonColl.size()==1) ) return;
     if( !(electronColl.at(0).Charge()!=muonColl.at(0).Charge()) ) return;
     if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>15) ) return;
     
     //OS emu above trig turnon
      FillHist("Cutflow_EMu", 0., weight, 0., 10., 10);

     //1j cut
     if( jetColl.size()<1 ) return;
      FillHist("Cutflow_EMu", 1., weight, 0., 10., 10);

     //2jcut
     if( jetColl.size()<2 ) return;
      FillHist("Cutflow_EMu", 2., weight, 0., 10., 10);

     //1bcut
     if( bjetColl.size()==0 ) return;
      FillHist("Cutflow_EMu", 3., weight, 0., 10., 10);

     //MET40
     if( met<40 ) return;
      FillHist("Cutflow_EMu", 4., weight, 0., 10., 10);

 
   }
/////////////////////////////////////////////////////////////////////////////////// 

return;

}// End of execute event loop
  


void Aug2017_EMuMuTrigCheck::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Aug2017_EMuMuTrigCheck::BeginCycle() throw( LQError ){
  
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

Aug2017_EMuMuTrigCheck::~Aug2017_EMuMuTrigCheck() {
  
  Message("In Aug2017_EMuMuTrigCheck Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}




float Aug2017_EMuMuTrigCheck::ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut){

  float PTCorr=Ele.Pt()*(1+Ele.PFRelIso(0.3));
  //float PTCorr=Ele.Pt()*(1+min(0,Ele.PFRelIso(0.3)-TightIsoCut));
  return PTCorr;

}




float Aug2017_EMuMuTrigCheck::FakeRateData(snu::KElectron Ele, TString ID){

  float FR=0., PTCorr=Ele.Pt()*(1.+Ele.PFRelIso(0.3)), fEta=fabs(Ele.Eta());
  if(ID=="POGMVAMTIsop06IPp025p05Sig4FakeLIsop4"){
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
  return FR;
} 


float Aug2017_EMuMuTrigCheck::GetFakeWeight(std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleLColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID){

  float fakeweight=-1.; int NLooseNotTight=0;
  TString EleFRLabel=""; bool JSTrilepFR=false;
  if(EleTID=="POGMVAMIP"      && EleLID=="HctoWAFakeLoose") EleFRLabel="POGMVAMTIsop06IPp025p05Sig4FakeLIsop4";
  if(MuTID =="HNTrilepTight2" && MuLID =="HNTrilepFakeL2" ) JSTrilepFR=true;

  for(int i=0; i<MuLColl.size(); i++){
    if(!PassIDCriteria(MuLColl.at(i), MuTID)){
      float FR=0.;
      if(JSTrilepFR){
//        FR=m_datadriven_bkg->GetFakeObj()->getTrilepFakeRate_muon(false, MuLColl.at(i).Pt(), MuLColl.at(i).Eta());
      }
      fakeweight*=-FR/(1-FR);
      NLooseNotTight++;
    }
  }
  for(int i=0; i<EleLColl.size(); i++){
    if(!PassIDCriteria(EleLColl.at(i), EleTID)){
      float FR=FakeRateData(EleLColl.at(i), EleFRLabel);
     
      fakeweight*=-FR/(1-FR);
      NLooseNotTight++;
    }
  }
  if(NLooseNotTight==0) return 0.;

  return fakeweight;
}


float Aug2017_EMuMuTrigCheck::GetPreTrigPURW(int Nvtx){

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



void Aug2017_EMuMuTrigCheck::FillCutFlow(TString cut, float weight){
  
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



void Aug2017_EMuMuTrigCheck::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Aug2017_EMuMuTrigCheck::MakeHistograms(){
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
  // *  Remove//Overide this Aug2017_EMuMuTrigCheckCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Aug2017_EMuMuTrigCheck::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
