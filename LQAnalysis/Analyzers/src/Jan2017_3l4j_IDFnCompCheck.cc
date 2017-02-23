// $Id: Jan2017_3l4j_IDFnCompCheck.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQJan2017_3l4j_IDFnCompCheck Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Jan2017_3l4j_IDFnCompCheck.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Jan2017_3l4j_IDFnCompCheck);

 Jan2017_3l4j_IDFnCompCheck::Jan2017_3l4j_IDFnCompCheck() : AnalyzerCore(), out_muons(0) {

   SetLogName("Jan2017_3l4j_IDFnCompCheck");
   Message("In Jan2017_3l4j_IDFnCompCheck constructor", INFO);
   InitialiseAnalysis();
 }


 void Jan2017_3l4j_IDFnCompCheck::InitialiseAnalysis() throw( LQError ) {
   
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

void Jan2017_3l4j_IDFnCompCheck::ExecuteEvents()throw( LQError ){

////Basic Infos///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Apply the gen weight 
   float weight_nopu=1; 
   if(!isData) weight*=MCweight;
   if(!isData) weight_nopu*=MCweight;
   FillHist("GenWeight", MCweight, 1, -10, 10, 1000);

   //Total Event  
   FillCutFlow("NoCut", weight);

   /// Acts on data to remove bad reconstructed event 
   //if(isData&& (! eventbase->GetEvent().LumiMask(lumimask))) return;

   m_logger<<DEBUG<<"RunNumber/Event Number = "<<eventbase->GetEvent().RunNumber()<<" : "<<eventbase->GetEvent().EventNumber()<<LQLogger::endmsg;
   m_logger<<DEBUG<<"isData = "<<isData<<LQLogger::endmsg;


   //Numbet of Vertex NoPUreweight
   FillHist("Nvtx_nocut_NoRW", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);

   //Pileup Reweight
   float pileup_reweight=(1.0);
   //if (!k_isdata) { pileup_reweight = TempPileupWeight(); weight*=pileup_reweight;}

   //Numbet of Vertex PUreweight
   FillHist("Nvtx_nocut_PURW", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);


   bool emumu_analysis=true, trimu_analysis=false;
//   bool emumu_analysis=false, trimu_analysis=true;


   //Trigers
   std::vector<TString> triggerslist, triggerlist1, triggerlist2;//  ListTriggersAvailable();

   TString analysis_trigger;
   if(emumu_analysis) {analysis_trigger="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";}
   else if(trimu_analysis) analysis_trigger="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v";
   //if(!PassTrigger(analysis_trigger)) return; //Not for now..

   //PassTrigger Uses BeginWith so don't stick to exact name

   float trigger_ps_weight=1, weight_trigger_sf=1;
/*   std::vector<snu::KMuon> trigmuColl; std::vector<snu::KElectron> trigeColl;
     eventbase->GetMuonSel()->SetPt(20.);eventbase->GetMuonSel()->SetEta(2.4);eventbase->GetMuonSel()->SetRelIso(0.15);eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);  eventbase->GetMuonSel()->Selection(trigmuColl);
     eventbase->GetElectronSel()->SetPt(20.);eventbase->GetElectronSel()->SetEta(2.4);eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MEDIUM); eventbase->GetElectronSel()->Selection(trigeColl);*/

   trigger_ps_weight=WeightByTrigger(analysis_trigger, TargetLumi);
   //weight_trigger_sf=1.;//TriggerScaleFactor(trigeColl, trigmuColl, analysis_trigger);



   FillHist("TriggerSFWeight" , weight_trigger_sf, 1., 0. , 2., 200); FillHist("TriggerPSWeight", trigger_ps_weight, 1., 0., 500., 500);
   weight*=weight_trigger_sf*trigger_ps_weight;
   weight_nopu*=weight_trigger_sf*trigger_ps_weight;
   //FillCutFlow("TriggerCut", weight); //;in trigger Cut mode
   //FillHist("Basic_prescale", prescale, 1., 0., 2000., 2000);

   //Initial Event Cut : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
   //if(!PassMETFilter()) return;  FillCutFlow("EventCut", weight);

   //Vertex Cut
   //if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; FillCutFlow("VertexCut", weight);
   //Good Primary vtx def:(vtx.ndof()>4&&maxAbsZ<=0)||std::abs(vtx.z())<= 24)&&((maxd0 <=0)||std::abs(vtx.position().rho())<=2)&&!(vtx.isFake()))  



///////Objects in Analysis/////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

   //Primary Object Collection
   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
  
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.25);
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.15);
//     eventbase->GetMuonSel()->SetBSdxy(999.);                eventbase->GetMuonSel()->SetBSdz(999.);
//     eventbase->GetMuonSel()->SetChiNdof(999.);
   std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->Selection(muonColl);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
 
   std::vector<snu::KMuon> muonLQTightColl;  eventbase->GetMuonSel()->SelectMuons(muonLQTightColl, "MUON_POG_TIGHT_TEST", 10., 2.4);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
   std::vector<snu::KMuon> muonLQTightColl2; eventbase->GetMuonSel()->SelectMuons(muonLQTightColl2, "MUON_POG_TIGHT", 10., 2.4);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
   std::vector<snu::KMuon> muonLQLooseColl;  eventbase->GetMuonSel()->SelectMuons(muonLQLooseColl, "MUON_POG_LOOSE", 10., 2.4);
   
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0994, 0.107);//ISO Ichep16 Loose WP
   //eventbase->GetElectronSel()->SetdxyBEMax(0.005, 0.01);  eventbase->GetElectronSel()->SetdzBEMax(0.41, 0.822);
//cout << eventbase->GetElectronSel()->pt_cut_min<<endl;
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);//pt>10/eta<2.4(hole region ~1.4 excluded by default)/EGAMMA_LOOSE/NoExtIso*/
//cout << eventbase->GetElectronSel()->pt_cut_min<<endl;
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//Iso Ichep16 Loose WP
   /*eventbase->GetElectronSel()->SetdxyBEMax(0.005, 0.01);  eventbase->GetElectronSel()->SetdzBEMax(0.0466, 0.417);*/
   std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->Selection(electronColl);//pt>10/eta<2.4(hole region ~1.4 excluded by default)/EGAMMA_LOOSE/NoExtIso*/

  std::vector<snu::KElectron> electronLQLooseColl; eventbase->GetElectronSel()->SelectElectrons(electronLQLooseColl, "ELECTRON_POG_16LOOSE", 25., 2.5);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
  std::vector<snu::KElectron> electronLQTightColl; eventbase->GetElectronSel()->SelectElectrons(electronLQTightColl, "ELECTRON_POG_16TIGHT", 25., 2.5);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_TIGHT);
     eventbase->GetJetSel()->SetPt(20.);                     eventbase->GetJetSel()->SetEta(2.4);
//   eventbase->GetJetSel()->SetUseJetPileUp(true);
     bool LeptonVeto=true;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonColl, electronColl);
// std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->SelectJets(isData, jetColl, muonColl, electronColl, "PFJET_TIGHT", 20., 2.4, true);


//   std::vector<int> bIdxColl=GetSFBJetIdx(jetColl,"Medium");
//   std::vector<int> ljIdxColl=GetSFLJetIdx(jetColl, bIdxColl, "Medium");

//   std::vector<snu::KJet> bjetColl; for(int i=0; i<bIdxColl.size(); i++){bjetColl.push_back(jetColl.at(bIdxColl.at(i)));}
//   std::vector<snu::KJet> ljetColl; for(int i=0; i<ljIdxColl.size(); i++){ljetColl.push_back(jetColl.at(ljIdxColl.at(i)));}

   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");
   std::vector<snu::KJet> ljetColl = SelLightJets(jetColl, "Medium");

   //Temporarily Placed Cuts(Activate ones at original place if trig simul. fully available
   //METFilter
   if(!PassMETFilter()) return;  FillCutFlow("EventCut", weight);

   //Vertex Cut
   if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; FillCutFlow("VertexCut", weight);

   //Trigger eff. emul. and cuts
   float trigger_eff=1.;
   if(isData){
     if(!PassTrigger(analysis_trigger)) return;
   }
   else{
     weight*=trigger_eff;
     weight_nopu*=trigger_eff;
   }
   FillCutFlow("TriggerCut", weight); //;in trigger Cut mode
   /////////////////////////////////////////////////////////////////////////////////////////


   bool emu=false, mumu=false;
   int mu1_Ai=-1, mu2_Ai=-1, mu_Wi=-1;
   int nbjets=bjetColl.size(); int njets=jetColl.size(); const int nljets=ljetColl.size();//njets-nbjets;//number of light jets
   double Pzv1, Pzv2;

   double met = eventbase->GetEvent().MET();
   double met_x = eventbase->GetEvent().MET()*TMath::Cos(eventbase->GetEvent().METPhi());
   double met_y = eventbase->GetEvent().MET()*TMath::Sin(eventbase->GetEvent().METPhi());
   int Nvtx=eventbase->GetEvent().nVertices();
   double Pzv;
   snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));
//   snu::KParticle v[4]; v[0].SetPx(met_x); v[0].SetPy(met_y);
//                        v[1].SetPx(met_x); v[1].SetPy(met_y);

   //Scale Factors
   float id_weight=1., reco_weight=1., iso_weight=1.;
   if(!isData){
//     id_weight *= ElectronScaleFactor(BaseSelection::ELECTRON_POG_TIGHT, electronColl);
//     id_weight *= MuonScaleFactor(BaseSelection::MUON_POG_TIGHT, muonColl,0);
//     iso_weight *= MuonISOScaleFactor(BaseSelection::MUON_POG_TIGHT, muonColl,0);
//     reco_weight *= ElectronRecoScaleFactor(electronColl);
//     weight*=id_weight*reco_weight*iso_weight;
   }





//////Basic Objects Distribution////////////////////////////////////////////////////////////////////////
   FillHist("Basic_Nj_orig", njets, weight, 0., 10., 10);
   FillHist("Basic_Nb_orig", nbjets, weight, 0., 10., 10);//Bjet def CVSincV2>0.89;pfKJet::CSVv2IVF_discriminator 
   FillHist("Basic_NmuT_orig", muonColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NmuL_orig", muonLooseColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeT_orig", electronColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeL_orig", electronLooseColl.size(), weight, 0., 10., 10);
   if(muonLooseColl.size()!=muonLQLooseColl.size()) FillHist("DiffSizeMu_LQL_L", 0.5, 1, 0., 1., 1);
   if(muonColl.size()!=muonLQTightColl.size()) FillHist("DiffSizeMu_LQT_T", 0.5, 1, 0., 1., 1);
   if(muonLooseColl.size()>0 && muonLQLooseColl.size()>0) FillHist("RatioSizeMu_LQL_L", muonLQLooseColl.size()/muonLooseColl.size(), 1, 0., 10., 1000);
   if(muonColl.size()>0 && muonLQTightColl.size()>0) FillHist("RatioSizeMu_LQT_T", muonLQTightColl.size()/muonColl.size(), 1, 0., 10., 1000);
   if(muonLQTightColl.size()!=muonLQTightColl2.size()) FillHist("DiffSizeMu_LQT1_LQT2", 0.5, 1, 0., 1., 1);

   if(electronLooseColl.size()!=electronLQLooseColl.size()) FillHist("DiffSizeEle_LQL_L", 0.5, 1, 0., 1., 1);
   if(electronColl.size()!=electronLQTightColl.size()) FillHist("DiffSizeEle_LQT_T", 0.5, 1, 0., 1., 1);
   if(electronLooseColl.size()>0 && electronLQLooseColl.size()>0) FillHist("RatioSizeEle_LQL_L", electronLQLooseColl.size()/electronLooseColl.size(), 1, 0., 10., 1000);
   if(electronColl.size()>0 && electronLQTightColl.size()>0) FillHist("RatioSizeEle_LQT_T", electronLQTightColl.size()/electronColl.size(), 1, 0., 10., 1000);

/*
   if(electronColl.size()==1) FillHist("Basic_NmuT_weT", muonColl.size(), weight, 0., 10., 10);
   if(trimu_analysis){
     if(muonColl.size()>0) FillHist("Basic_Ptmu1_3mu_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>1) FillHist("Basic_Ptmu2_3mu_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>2) FillHist("Basic_Ptmu3_3mu_orig", muonColl.at(2).Pt(), weight, 0., 200., 200);
   }
   if(emumu_analysis){
     if(muonColl.size()>0) FillHist("Basic_Ptmu1_1e2mu_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>1) FillHist("Basic_Ptmu2_1e2mu_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
     if(electronColl.size()>0) FillHist("Basic_Pte_1e2mu_orig", electronColl.at(0).Pt(), weight, 0., 200., 200);
   }
   if(njets>0) FillHist("Basic_j1_Et_orig", jetColl.at(0).Et(), weight, 0., 200., 200);
   if(njets>1) FillHist("Basic_j2_Et_orig", jetColl.at(1).Et(), weight, 0., 200., 200);
   if(njets>2) FillHist("Basic_j3_Et_orig", jetColl.at(2).Et(), weight, 0., 200., 200);
   if(njets>3) FillHist("Basic_j4_Et_orig", jetColl.at(3).Et(), weight, 0., 200., 200);
   if(nbjets>0) FillHist("Basic_b1_Et_orig", bjetColl.at(0).Et(), weight, 0, 200., 200);
   if(nbjets>1) FillHist("Basic_b2_Et_orig", bjetColl.at(1).Et(), weight, 0., 200., 200);
   FillHist("Basic_METdist_orig", met, weight, 0., 200., 100);//It is orig because even though I havent put METcut below, if there is correlation between Nl ,Nj, Nb and MET then MET distrubution will also be affected by selection. The distribution remains the same only when the correlation is absent.

///////Event Selection&Histograms//////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////   
   if((emumu_analysis)&&(electronColl.size()==1)&&(muonColl.size()==1)) emu=true;
   else if((trimu_analysis)&&(muonColl.size()==2)) mumu=true;
   else return; //veto if not emu or mumu
   FillCutFlow("NlCut", weight);

   if((trimu_analysis)&&(SumCharge(muonColl)!=0)) return;
   FillCutFlow("OS(2mu)", weight);

   FillHist("Basic_Nvtx_NoRW_wNlOScut", Nvtx, weight_nopu, 0., 50., 50);
   FillHist("Basic_Nvtx_PURW_wNlOScut", Nvtx, weight, 0., 50., 50);

   FillHist("Basic_Nj_wNlOScut", njets, weight, 0., 10., 10);
   FillHist("Basic_Nb_wNlOScut", nbjets, weight, 0., 10., 10);
   FillHist("Basic_MET_wNlOScut", met, weight, 0., 200., 200);
   if(emumu_analysis){
     FillHist("Basic_Ne_wNlOScut", electronColl.size(), weight, 0., 10., 10);
     FillHist("Basic_Nmu_wNlOScut", muonColl.size(), weight, 0., 10., 10);
     FillHist("Basic_Pte_wNlOScut", electronColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("Basic_Ptmu1_wNlOScut", muonColl.at(0).Pt(), weight, 0., 200., 200);
   }
   else if(trimu_analysis){
     FillHist("Basic_Ne_wNlOScut", electronColl.size(), weight, 0., 10., 10);
     FillHist("Basic_Nmu_wNlOScut", muonColl.size(), weight, 0., 10., 10);
     FillHist("Basic_Ptmu1_wNlOScut", muonColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("Basic_Ptmu2_wNlOScut", muonColl.at(1).Pt(), weight, 0., 200., 200);
     FillHist("Basic_Mmumu_wNlOScut", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 200., 200);
   }

   if(njets<3) return;
     FillCutFlow("NjCut", weight);
     FillHist("Basic_Nb_wNljcut", nbjets, weight, 0., 10., 10);

   if(bjetColl.size()==0) return;
     FillCutFlow("NbCut", weight);

   if(nljets<2) return;
     FillCutFlow("NljCut", weight);

   if(trimu_analysis){ if(muonColl.at(0).Pt()<20) return;}
   else if(emumu_analysis){ if(electronColl.at(0).Pt()<20) return;}
   FillCutFlow("lPtCut", weight); //leading lepton Pt>20GeV

/////////////////////////////////////////////////////////////////////////////////// 

///Validation/////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

   //Basic Kinematic Variables//////////////////////////////////////////////////////
   FillHist("Nvtx", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);
   FillHist("MET", met, weight, 0., 200., 100);
   FillHist("Nj", njets, weight, 0., 10., 10);//Njets after all cut
   FillHist("Nb", nbjets, weight, 0., 5., 5);//Nbjet after all cuts

   if(emumu_analysis){
     FillHist("Pte", electronColl.at(0).Pt(), weight, 0., 500., 500);
     FillHist("Ptmu1", muonColl.at(0).Pt(), weight, 0., 500., 500);
     FillHist("Etae", electronColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etamu1", muonColl.at(0).Eta(), weight, -5., 5., 100);

     FillHist("Ptj1", jetColl.at(0).Pt(), weight, 0., 500., 500);
     FillHist("Ptj2", jetColl.at(1).Pt(), weight, 0., 500., 500);
     FillHist("Ptj3", jetColl.at(2).Pt(), weight, 0., 500., 500);
     if(njets>3) FillHist("Ptj4", jetColl.at(3).Pt(), weight, 0., 500., 500);
     FillHist("Ptb1", bjetColl.at(0).Pt(), weight, 0., 500., 500);
     if(nbjets>1) FillHist("Ptb2", bjetColl.at(1).Pt(), weight, 0., 500., 500);
     FillHist("Ptlj1", ljetColl.at(0).Pt(), weight, 0., 500., 500);
     FillHist("Ptlj2", ljetColl.at(1).Pt(), weight, 0., 500., 500);
  
     FillHist("Etaj1", jetColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etaj2", jetColl.at(1).Eta(), weight, -5., 5., 100);
     FillHist("Etaj3", jetColl.at(2).Eta(), weight, -5., 5., 100);
     if(njets>3) FillHist("Etaj4", jetColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etab1", bjetColl.at(0).Eta(), weight, -5., 5., 100);
     if(nbjets>1) FillHist("Etab2", bjetColl.at(1).Eta(), weight, -5., 5., 100);
     FillHist("Etalj1", ljetColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etalj2", ljetColl.at(1).Eta(), weight, -5., 5., 100);
   }
   else if(trimu_analysis){
     FillHist("Ptmu1", muonColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("Ptmu2", muonColl.at(1).Pt(), weight, 0., 200., 200);
     FillHist("Etamu1", muonColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etamu2", muonColl.at(1).Eta(), weight, -5., 5., 100);

     FillHist("Ptj1", jetColl.at(0).Pt(), weight, 0., 500., 500);
     FillHist("Ptj2", jetColl.at(1).Pt(), weight, 0., 500., 500);
     FillHist("Ptj3", jetColl.at(2).Pt(), weight, 0., 500., 500);
     if(njets>3) FillHist("Ptj4", jetColl.at(0).Pt(), weight, 0., 500., 500);
     FillHist("Ptb1", bjetColl.at(0).Pt(), weight, 0., 500., 500);
     if(nbjets>1) FillHist("Ptb2", bjetColl.at(1).Pt(), weight, 0., 500., 500);
     FillHist("Ptlj1", ljetColl.at(0).Pt(), weight, 0., 500., 500);
     FillHist("Ptlj2", ljetColl.at(1).Pt(), weight, 0., 500., 500);
  
     FillHist("Etaj1", jetColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etaj2", jetColl.at(1).Eta(), weight, -5., 5., 100);
     FillHist("Etaj3", jetColl.at(2).Eta(), weight, -5., 5., 100);
     if(njets>3) FillHist("Etaj4", jetColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etab1", bjetColl.at(0).Eta(), weight, -5., 5., 100);
     if(nbjets>1) FillHist("Etab2", bjetColl.at(1).Eta(), weight, -5., 5., 100);
     FillHist("Etalj1", ljetColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etalj2", ljetColl.at(1).Eta(), weight, -5., 5., 100);

     FillHist("Mmumu", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
     FillHist("Mmumu_20S", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 20, 400);

     //Jet Assignment//////////////////////////
     //j_W candidate decision
     float bestM=0;
     int j1_Wi=-1, j2_Wi=-1;
     for(int i=0; i<nljets; i++) {
        for(int j=i+1; j<nljets; j++) {
           float tmpM=(ljetColl.at(i)+ljetColl.at(j)).M();
           if(i==0) {bestM=tmpM; j1_Wi=i; j2_Wi=j;}
           if(fabs(tmpM-80.4)<fabs(bestM-80.4)) {bestM=tmpM; j1_Wi=i; j2_Wi=j;}
        }
     }
     FillHist("Mmumujj", (muonColl.at(0)+muonColl.at(1)+ljetColl.at(j1_Wi)+ljetColl.at(j2_Wi)).M(), weight, 0., 1000., 1000);   

     //After RochCorr
     CorrectMuonMomentum(muonColl);
     FillHist("Ptmu1_RochCorr", muonColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("Ptmu2_RochCorr", muonColl.at(1).Pt(), weight, 0., 200., 200);
     FillHist("Etamu1_RochCorr", muonColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Etamu2_RochCorr", muonColl.at(1).Eta(), weight, -5., 5., 100);
     FillHist("Mmumu_RochCorr", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
     FillHist("Mmumu_10S_RochCorr", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 20., 400);
     FillHist("Mmumujj_RochCorr", (muonColl.at(0)+muonColl.at(1)+ljetColl.at(j1_Wi)+ljetColl.at(j2_Wi)).M(), weight, 0., 1000., 1000);
   } 

*/
return;
}// End of execute event loop
  


void Jan2017_3l4j_IDFnCompCheck::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Jan2017_3l4j_IDFnCompCheck::BeginCycle() throw( LQError ){
  
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

Jan2017_3l4j_IDFnCompCheck::~Jan2017_3l4j_IDFnCompCheck() {
  
  Message("In Jan2017_3l4j_IDFnCompCheck Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Jan2017_3l4j_IDFnCompCheck::FillCutFlow(TString cut, float weight){
  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 12, 0., 12.);
    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"NlCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"OS(2mu)");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"NjCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(8,"NbCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(9,"NljCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(10,"lPtCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(11,"MjjCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(12,"MmumuCut");
    
  }
}



void Jan2017_3l4j_IDFnCompCheck::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Jan2017_3l4j_IDFnCompCheck::MakeHistograms(){
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
/*
  //After all preselections
  AnalyzerCore::MakeHistograms("Basic_MET_wPreSel", 100, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Nj_wPreSel", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Nb_wPreSel", 5, 0., 5.);
  AnalyzerCore::MakeHistograms("Basic_Jet1_Et", 50, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Jet2_Et", 50, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Jet3_Et", 50, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Jet4_Et", 50, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Jet1_Eta", 20, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Jet2_Eta", 20, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Jet3_Eta", 20, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Jet4_Eta", 20, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_bJet1_Et", 20, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_bJet2_Et", 20, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_bJet1_Eta", 20, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_bJet2_Eta", 20, -5., 5.);

  AnalyzerCore::MakeHistograms("Basic_Pte_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu1_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu2_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Etae_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etamu1_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etamu2_wPreSel", 100, -5., 5.);

  AnalyzerCore::MakeHistograms("Basic_Ptmu3_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Etamu3_wPreSel", 100, -5., 5.);

  AnalyzerCore::MakeHistograms("Basic_Ptj1_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptj2_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptj3_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptj4_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptb1_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptb2_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptlj1_wPreSel", 200, 0., 200.);

  AnalyzerCore::MakeHistograms("Basic_Etaj1_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etaj2_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etaj3_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etaj4_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etab1_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etab2_wPreSel", 100, -5., 5.);
*/

  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Jan2017_3l4j_IDFnCompCheckCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Jan2017_3l4j_IDFnCompCheck::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
