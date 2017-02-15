// $Id: Sep2016_3l4j_DataMCComp.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQSep2016_3l4j_DataMCComp Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Sep2016_3l4j_DataMCComp.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Sep2016_3l4j_DataMCComp);

 Sep2016_3l4j_DataMCComp::Sep2016_3l4j_DataMCComp() : AnalyzerCore(), out_muons(0) {

   SetLogName("Sep2016_3l4j_DataMCComp");
   Message("In Sep2016_3l4j_DataMCComp constructor", INFO);
   InitialiseAnalysis();
 }


 void Sep2016_3l4j_DataMCComp::InitialiseAnalysis() throw( LQError ) {
   
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

void Sep2016_3l4j_DataMCComp::ExecuteEvents()throw( LQError ){

////Basic Infos///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Apply the gen weight  
   if(!isData) weight*=MCweight;
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
   if (!k_isdata) { pileup_reweight = TempPileupWeight(); weight*=pileup_reweight;}

   //Numbet of Vertex NoPUreweight
   FillHist("Nvtx_nocut_PURW", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);


   bool emumu_analysis=false, trimu_analysis=true;

   //Trigers
   std::vector<TString> triggerslist, triggerlist1, triggerlist2;//  ListTriggersAvailable();

   TString analysis_trigger;
   if(emumu_analysis) {analysis_trigger="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";}
   else if(trimu_analysis) analysis_trigger="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v";
   //if(!PassTrigger(analysis_trigger)) return; //Not for now..

   //PassTrigger Uses BeginWith so don't stick to exact name

   float trigger_ps_weight, weight_trigger_sf;
/*   std::vector<snu::KMuon> trigmuColl; std::vector<snu::KElectron> trigeColl;
     eventbase->GetMuonSel()->SetPt(20.);eventbase->GetMuonSel()->SetEta(2.4);eventbase->GetMuonSel()->SetRelIso(0.15);eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);  eventbase->GetMuonSel()->Selection(trigmuColl);
     eventbase->GetElectronSel()->SetPt(20.);eventbase->GetElectronSel()->SetEta(2.4);eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MEDIUM); eventbase->GetElectronSel()->Selection(trigeColl);*/

   trigger_ps_weight=WeightByTrigger(analysis_trigger, TargetLumi);
   weight_trigger_sf=1.;//TriggerScaleFactor(trigeColl, trigmuColl, analysis_trigger);



   FillHist("TriggerSFWeight" , weight_trigger_sf, 1., 0. , 2., 200); FillHist("TriggerPSWeight", trigger_ps_weight, 1., 0., 500., 500);
   weight*=weight_trigger_sf*trigger_ps_weight;
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
//   eventbase->GetMuonSel()->SetBSdxy(0.005);               eventbase->GetMuonSel()->SetBSdz(0.03);
//   eventbase->GetMuonSel()->SetChiNdof(3);
   std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->Selection(muonColl);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
// std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->SelectMuons(muonColl, BaseSelection::MUON_POG_LOOSE, 10., 2.4);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
   
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
     eventbase->GetElectronSel()->SetPt(20.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0994, 0.107);//ISO Ichep16 Loose WP
   //eventbase->GetElectronSel()->SetdxyBEMax(0.005, 0.01);  eventbase->GetElectronSel()->SetdzBEMax(0.005, 0.01);
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);//pt>10/eta<2.4(hole region ~1.4 excluded by default)/EGAMMA_LOOSE/NoExtIso*/

     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(20.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//Iso Ichep16 Loose WP
//   eventbase->GetElectronSel()->SetdxyBEMax(0.005, 0.01);  eventbase->GetElectronSel()->SetdzBEMax(0.005, 0.01);
   std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->Selection(electronColl);//pt>10/eta<2.4(hole region ~1.4 excluded by default)/EGAMMA_LOOSE/NoExtIso*/
// std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->SelectElectrons(electronColl, BaseSelection::ELECTRON_POG_TIGHT, 20., 2.4);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_TIGHT);
     eventbase->GetJetSel()->SetPt(20.);                     eventbase->GetJetSel()->SetEta(2.4);
//   eventbase->GetJetSel()->SetUseJetPileUp(true);
     bool LeptonVeto=true;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonColl, electronColl, isData, true);
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
     if(emumu_analysis) trigger_eff = TriggerEff(analysis_trigger, muonColl, electronColl);
     else if (trimu_analysis) trigger_eff = TriggerEff(analysis_trigger, muonColl);
     weight*=trigger_eff;
   }
   FillCutFlow("TriggerCut", weight); //;in trigger Cut mode
   /////////////////////////////////////////////////////////////////////////////////////////

   bool emumu=false, trimu=false;
   int mu1_Ai=-1, mu2_Ai=-1, mu_Wi=-1, j1_Wi=-1, j2_Wi=-1;
   int nbjets=bjetColl.size(); int njets=jetColl.size(); const int nljets=njets-nbjets;//number of light jets

   double met = eventbase->GetEvent().MET();
   double met_x = eventbase->GetEvent().MET()*TMath::Cos(eventbase->GetEvent().METPhi());
   double met_y = eventbase->GetEvent().MET()*TMath::Sin(eventbase->GetEvent().METPhi());
   double Pzv, Pzv_1=0, Pzv_2=0, Pzv_absS=0, Pzv_absL=0, Pzv_dRS=0, Pzv_dRL=0,Pzv_truth=0;
   double Pzv_1_truth=0, Pzv_2_truth=0, Pzv_absS_truth=0;
   snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));
   snu::KParticle v_absS, v_truth;

/*
   //Scale Factors
   float id_weight=1., reco_weight=1., iso_weight=1.;
   if(!isData){
     id_weight *= ElectronScaleFactor(BaseSelection::ELECTRON_POG_TIGHT, electronColl);
     id_weight *= MuonScaleFactor(BaseSelection::MUON_POG_TIGHT, muonColl,0);
     iso_weight *= MuonISOScaleFactor(BaseSelection::MUON_POG_TIGHT, muonColl,0);
     reco_weight *= ElectronRecoScaleFactor(electronColl);
     weight*=id_weight*reco_weight*iso_weight;
   }
*/



   ///Basic Objects Distribution// 
   FillHist("Basic_Nj_orig", njets, weight, 0., 10., 10);
   FillHist("Basic_Nb_orig", nbjets, weight, 0., 10., 10);//Bjet def CVSincV2>0.89;pfKJet::CSVv2IVF_discriminator 
   FillHist("Basic_NmuT_orig", muonColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NmuL_orig", muonLooseColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeT_orig", electronColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeL_orig", electronLooseColl.size(), weight, 0., 10., 10);
   if(electronColl.size()==1) FillHist("Basic_NmuT_weT", muonColl.size(), weight, 0., 10., 10);
   if(trimu_analysis){
     if(muonColl.size()>0) FillHist("Basic_Ptmu1_3mu_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>1) FillHist("Basic_Ptmu2_3mu_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>2) FillHist("Basic_Ptmu3_3mu_orig", muonColl.at(2).Pt(), weight, 0., 200., 200);
   }
   if(emumu_analysis){
     if(muonColl.size()>0) FillHist("Basic_Ptmu1_1e2mu_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>1) FillHist("Basic_Ptmu2_1e2mu_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
     if(electronColl.size()>0) FillHist("Basic_Pte1_1e2mu_orig", electronColl.at(0).Pt(), weight, 0., 200., 200);
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
   if((emumu_analysis)&&(electronColl.size()==1)&&(muonColl.size()==2)) emumu=true;
   else if ((trimu_analysis)&&(electronColl.size()==0)&&(muonColl.size()==3)) trimu=true;
   else return; //veto if not emumu or mumumu
   FillCutFlow("NlCut", weight);

   if((emumu_analysis)&&(SumCharge(muonColl)==0)) FillHist("Basic_NmuT2eT1_wOS", 0., weight, 0., 1., 1);
   else if((trimu_analysis)&&(fabs(SumCharge(muonColl))==1)) FillHist("Basic_NmuT3eT0_wOS", 0., weight, 0., 1., 1);
   else return;
   FillCutFlow("OSmuon", weight);

   if(trimu_analysis){ if(muonColl.at(0).Pt()<20) return; }
   FillHist("Basic_Nvtx_wNlOScut", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);
   FillHist("Basic_Nj_wNlOScut", njets, weight, 0., 10., 10);
   FillHist("Basic_Nb_wNlOScut", nbjets, weight, 0., 10., 10);
//Data-MC Comparison at Nl,OScut sample/////
   FillHist("Basic_mu1_Pt_wNlOScut", muonColl.at(0).Pt(), weight, 0., 200., 200);
   FillHist("Basic_mu1_Eta_wNlOScut", muonColl.at(0).Eta(), weight, -5., 5., 100);
   FillHist("Basic_mu2_Pt_wNlOScut", muonColl.at(1).Pt(), weight, 0., 200., 200);
   FillHist("Basic_mu2_Eta_wNlOScut", muonColl.at(1).Eta(), weight, -5., 5., 100);
   if(emumu_analysis){
     FillHist("Basic_e_Pt_wNlOScut", electronColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("Basic_e_Eta_wNlOScut", electronColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Basic_Mmumu_wNlOScut_1e2mu", (muonColl.at(0)+muonColl.at(1)).M(), 0., 1000., 1000);
     FillHist("Basic_MTev_wNlOScut_1e2mu", (electronColl.at(0)+v).Mt(), 0., 1000., 1000);
   }
   else if(trimu_analysis){
     FillHist("Basic_mu3_Pt_wNlOScut", muonColl.at(2).Pt(), weight, 0., 200., 200);
     FillHist("Basic_mu3_Eta_wNlOScut", muonColl.at(2).Eta(), weight, -5., 5., 100);
     int mu1_Zi, mu2_Zi; float tmp=1000;
     for(int i=0; i<3; i++){
       for(int j=i+1; j<3; j++){
         if(fabs((muonColl.at(i)+muonColl.at(j)).M()-90)<tmp) {tmp=fabs((muonColl.at(i)+muonColl.at(j)).M()-90); mu1_Zi=i; mu2_Zi=j;}
       }
     }
     int mu_Witmp=3-mu1_Zi-mu2_Zi;
     FillHist("Basic_Mmumu_wNlOScut_3mu", (muonColl.at(mu1_Zi)+muonColl.at(mu2_Zi)).M(), weight, 0., 200., 200);
     FillHist("Basic_MTmuv_wNlOScut_3mu", (muonColl.at(mu_Witmp)+v).Mt(), weight, 0., 1000., 1000);
   }
   
   
/////////////////////////////////////////////                  
   if(njets<3) return;
     FillCutFlow("NjCut", weight);
     FillHist("Basic_Nb_wNljcut", nbjets, weight, 0., 10., 10);

   if(bjetColl.size()==0) return;
     FillCutFlow("NbCut", weight);

   if(nljets<2) return;
     FillCutFlow("NljCut", weight);

   if(trimu_analysis){ if(muonColl.at(0).Pt()<20) return;}
   else if(emumu_analysis){
      if(electronColl.at(0).Pt()<20) return;
   }
   FillCutFlow("lPtCut", weight); //leading lepton Pt>20GeV


   //Kinematic Variables after all cuts//////////////////////////////////////////////////////
   FillHist("Basic_MET_wPreSel", met, weight, 0., 200., 100);
   FillHist("Basic_Nj_wPreSel", njets, weight, 0., 10., 10);//Njets after all cut
   FillHist("Basic_Nb_wPreSel", nbjets, weight, 0., 5., 5);//Nbjet after all cuts

   if(emumu_analysis){
     FillHist("Basic_Pte_wPreSel", electronColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("Basic_Ptmu1_wPreSel", muonColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("Basic_Ptmu2_wPreSel", muonColl.at(1).Pt(), weight, 0., 200., 200);
     FillHist("Basic_Etae_wPreSel", electronColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Basic_Etamu1_wPreSel", muonColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Basic_Etamu2_wPreSel", muonColl.at(1).Eta(), weight, -5., 5., 100);
   }
   else if(trimu_analysis){
     FillHist("Basic_Ptmu1_wPreSel", muonColl.at(0).Pt(), weight, 0., 200., 200);
     FillHist("Basic_Ptmu2_wPreSel", muonColl.at(1).Pt(), weight, 0., 200., 200);
     FillHist("Basic_Ptmu3_wPreSel", muonColl.at(2).Pt(), weight, 0., 200., 200);
     FillHist("Basic_Etamu1_wPreSel", muonColl.at(0).Eta(), weight, -5., 5., 100);
     FillHist("Basic_Etamu2_wPreSel", muonColl.at(1).Eta(), weight, -5., 5., 100);
     FillHist("Basic_Etamu3_wPreSel", muonColl.at(2).Eta(), weight, -5., 5., 100);
   } 
   FillHist("Basic_Ptj1_wPreSel", jetColl.at(0).Pt(), weight, 0., 200., 200);
   FillHist("Basic_Ptj2_wPreSel", jetColl.at(1).Pt(), weight, 0., 200., 200);
   FillHist("Basic_Ptj3_wPreSel", jetColl.at(2).Pt(), weight, 0., 200., 200);
   if(njets>3) FillHist("Basic_Ptj4_wPreSel", jetColl.at(0).Pt(), weight, 0., 200., 200);
   FillHist("Basic_Ptb1_wPreSel", bjetColl.at(0).Pt(), weight, 0., 200., 200);
   if(nbjets>1) FillHist("Basic_Ptb2_wPreSel", bjetColl.at(1).Pt(), weight, 0., 200., 200);
   FillHist("Basic_Ptlj1_wPreSel", ljetColl.at(0).Pt(), weight, 0., 200., 200);
   FillHist("Basic_Ptlj2_wPreSel", ljetColl.at(1).Pt(), weight, 0., 200., 200);

   FillHist("Basic_Etaj1_wPreSel", jetColl.at(0).Eta(), weight, -5., 5., 100);
   FillHist("Basic_Etaj2_wPreSel", jetColl.at(1).Eta(), weight, -5., 5., 100);
   FillHist("Basic_Etaj3_wPreSel", jetColl.at(2).Eta(), weight, -5., 5., 100);
   if(njets>3) FillHist("Basic_Etaj4_wPreSel", jetColl.at(0).Eta(), weight, -5., 5., 100);
   FillHist("Basic_Etab1_wPreSel", bjetColl.at(0).Eta(), weight, -5., 5., 100);
   if(nbjets>1) FillHist("Basic_Etab2_wPreSel", bjetColl.at(1).Eta(), weight, -5., 5., 100);
   FillHist("Basic_Etalj1_wPreSel", ljetColl.at(0).Eta(), weight, -5., 5., 100);
   FillHist("Basic_Etalj2_wPreSel", ljetColl.at(1).Eta(), weight, -5., 5., 100);

//////////////////////////////////////////////////////////////////////////////////////////////   
/////////////Analysis/////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////////////////////////////////////
   /////Hc->AW,A->mumu/W->lv////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////////

   //Jet Assignment//////////////////////////
   float tmp;
   for(int i=0; i<nljets; i++){
     for(int j=i+1; j<nljets; j++){
       if((i==0)&&(j==1)){ tmp=fabs((ljetColl.at(i)+ljetColl.at(j)).M()-80); j1_Wi=i; j2_Wi=j;}
       else if(fabs((ljetColl.at(i)+ljetColl.at(j)).M()-80)<tmp) {tmp=fabs((ljetColl.at(i)+ljetColl.at(j)).M()-80); j1_Wi=i; j2_Wi=j;}
     }
   }

   
   //emumu case////////////////////////////////////////////////////////////////////////////

   if(emumu_analysis){

     //v assignment
     if(emumu_analysis) {
       Pzv_1=GetvPz(v,electronColl.at(0),1); Pzv_2=GetvPz(v,electronColl.at(0),2);
       Pzv_absS=fabs(Pzv_1)<fabs(Pzv_2) ? Pzv_1 : Pzv_2;
       v_absS.SetXYZM(met_x,met_y,Pzv_absS,0);
     }

     /////////////////////////////////////////////////////


      FillHist("Mmumu_1e2mu", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 200., 200);          
      FillHist("MW_jj_1e2mu", (ljetColl.at(j1_Wi)+ljetColl.at(j2_Wi)).M(), weight, 0., 200., 50);
      FillHist("M3lv_1e2mu", (muonColl.at(0)+muonColl.at(1)+electronColl.at(0)+v_absS).M(), weight, 0., 1000., 1000);
      FillHist("M2l2j_1e2mu", (muonColl.at(0)+muonColl.at(1)+ljetColl.at(j1_Wi)+ljetColl.at(j2_Wi)).M(), weight, 0., 1000., 1000);
      FillHist("MW_lv_1e2mu", (electronColl.at(0)+v_absS).M(), weight, 0., 200., 50);

     float Mmumu=(muonColl.at(0)+muonColl.at(1)).M();
     float MAWindow=5;
     if(fabs(Mmumu-5)<MAWindow){
       FillHist("M3lv_atMdimu5", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)+v_absS).M(), weight, 0., 1000., 1000);
       FillHist("M2l2j_atMdimu5", (muonColl.at(0)+muonColl.at(1)+ljetColl.at(j1_Wi)+ljetColl.at(j2_Wi)).M(), weight, 0., 1000., 1000);
     }

     if(fabs(Mmumu-30)<MAWindow){
       FillHist("M3lv_atMdimu30", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)+v_absS).M(), weight, 0., 1000., 1000);
       FillHist("M2l2j_atMdimu30", (muonColl.at(0)+muonColl.at(1)+ljetColl.at(j1_Wi)+ljetColl.at(j2_Wi)).M(), weight, 0., 1000., 1000);
     }

   }
   if(trimu_analysis){
//Select leading SS muon as W muon
     int OSindex=TriMuChargeIndex(muonColl,"OS");
     int SSindex1=TriMuChargeIndex(muonColl,"SS1");
     int SSindex2=TriMuChargeIndex(muonColl, "SS2");

     //mu_Wi assignment by PT method      
     mu1_Ai=OSindex;
     mu_Wi=muonColl.at(SSindex1)>muonColl.at(SSindex2) ? SSindex1 : SSindex2;
     mu2_Ai=muonColl.at(SSindex1)<muonColl.at(SSindex2) ? SSindex1 : SSindex2;

     //v assignment by absS
     Pzv_1=GetvPz(v,muonColl.at(mu_Wi),1); Pzv_2=GetvPz(v,muonColl.at(mu_Wi),2);
     Pzv_absS=fabs(Pzv_1)<fabs(Pzv_2) ? Pzv_1 : Pzv_2;
     v_absS.SetXYZM(met_x,met_y,Pzv_absS,0);


     FillHist("Mmumu_3mu_PT", (muonColl.at(mu1_Ai)+muonColl.at(mu2_Ai)).M(), weight, 0., 1000., 1000);
     FillHist("M3lv_3mu_PT", (muonColl.at(0)+muonColl.at(1)+muonColl.at(2)+v_absS).M(), weight, 0., 1000., 1000);
     FillHist("M2l2j_3mu_PT", (muonColl.at(mu1_Ai)+muonColl.at(mu2_Ai)+ljetColl.at(j1_Wi)+ljetColl.at(j2_Wi)).M(), weight, 0., 1000., 1000);


     float Mmumu=(muonColl.at(mu1_Ai)+muonColl.at(mu2_Ai)).M();
     float MAWindow=5;
     if(fabs(Mmumu-5)<MAWindow){
       FillHist("M3lv_atMdimu5_PT", (muonColl.at(mu1_Ai)+muonColl.at(mu2_Ai)+muonColl.at(mu_Wi)+v_absS).M(), weight, 0., 1000., 1000);
       FillHist("M2l2j_atMdimu5_PT", (muonColl.at(mu1_Ai)+muonColl.at(mu2_Ai)+ljetColl.at(j1_Wi)+ljetColl.at(j2_Wi)).M(), weight, 0., 1000., 1000);
     }

     if(fabs(Mmumu-30)<MAWindow){
       FillHist("M3lv_atMdimu30_PT", (muonColl.at(0)+muonColl.at(1)+muonColl.at(2)+v_absS).M(), weight, 0., 1000., 1000);
       FillHist("M2l2j_atMdimu30_PT", (muonColl.at(mu1_Ai)+muonColl.at(mu2_Ai)+ljetColl.at(j1_Wi)+ljetColl.at(j2_Wi)).M(), weight, 0., 1000., 1000);
     }


     //Mmumu based assignment
     float M_OSSS1=(muonColl.at(OSindex)+muonColl.at(SSindex1)).M();
     float M_OSSS2=(muonColl.at(OSindex)+muonColl.at(SSindex2)).M();
     float M_OS1  =M_OSSS1 > M_OSSS2 ? M_OSSS1 : M_OSSS2;
     float M_OS2  =M_OSSS1 < M_OSSS2 ? M_OSSS1 : M_OSSS2;
     int mu1_Ai_Mmumu=OSindex;
     int mu_Wi_Mmumu =M_OSSS1>M_OSSS2 ? SSindex1 : SSindex2;
     int mu2_Ai_Mmumu=M_OSSS1>M_OSSS2 ? SSindex2 : SSindex1;
 
     FillHist("Mmumu_3mu_Mmumu", (muonColl.at(mu1_Ai_Mmumu)+muonColl.at(mu2_Ai_Mmumu)).M(), weight, 0., 1000., 1000);
     FillHist("M3lv_3mu_Mmumu", (muonColl.at(0)+muonColl.at(1)+muonColl.at(2)+v_absS).M(), weight, 0., 1000., 1000);
     FillHist("M2l2j_3mu_Mmumu", (muonColl.at(mu1_Ai_Mmumu)+muonColl.at(mu2_Ai_Mmumu)+ljetColl.at(j1_Wi)+ljetColl.at(j2_Wi)).M(), weight, 0., 1000., 1000);


     Mmumu=(muonColl.at(mu1_Ai_Mmumu)+muonColl.at(mu2_Ai_Mmumu)).M();
     MAWindow=5;
     if(fabs(Mmumu-5)<MAWindow){
       FillHist("M3lv_atMdimu5_Mmumu", (muonColl.at(mu1_Ai_Mmumu)+muonColl.at(mu2_Ai_Mmumu)+muonColl.at(mu_Wi_Mmumu)+v_absS).M(), weight, 0., 1000., 1000);
       FillHist("M2l2j_atMdimu5_Mmumu", (muonColl.at(mu1_Ai_Mmumu)+muonColl.at(mu2_Ai_Mmumu)+ljetColl.at(j1_Wi)+ljetColl.at(j2_Wi)).M(), weight, 0., 1000., 1000);
     }

     if(fabs(Mmumu-30)<MAWindow){
       FillHist("M3lv_atMdimu30_Mmumu", (muonColl.at(0)+muonColl.at(1)+muonColl.at(2)+v_absS).M(), weight, 0., 1000., 1000);
       FillHist("M2l2j_atMdimu30_Mmumu", (muonColl.at(mu1_Ai_Mmumu)+muonColl.at(mu2_Ai_Mmumu)+ljetColl.at(j1_Wi)+ljetColl.at(j2_Wi)).M(), weight, 0., 1000., 1000);
     }


   }

/*   
   Pzv1 = GetvPz(v[0], electronColl.at(0), 1);
   Pzv2 = GetvPz(v[1], electronColl.at(0), 2);
   v[0].SetPz(Pzv1);v[0].SetE(sqrt(Pzv1*Pzv1+met*met)); FillHist("MW1", (v[0]+electronColl.at(0)).M(), weight, 0, 200, 200);
   v[1].SetPz(Pzv2);v[1].SetE(sqrt(Pzv2*Pzv2+met*met)); FillHist("MW2", (v[1]+electronColl.at(0)).M(), weight, 0, 200, 200);

     double tmpMt1, tmpMt2, tmpMt, tmpMtot, tmpMtot1, tmpMtot2;
     double tmpMtx1, tmpMtx2; 
     int vti=-1, vtxi=-1, bti=-1, bti2=-1, btxi=-1, btxi2=-1;
     int Decaytlv=0, Decaytxlv=0;
     double deltatlv=0, deltatxlv=0;
     snu::KParticle tmpJ;  

     //THc side lv decay, tx side qqdecay scenario

     for(int i=0; i<4; i++){
        tmpMt1=(v[0]+electronColl.at(0)+muonColl.at(0)+muonColl.at(1)+jetColl.at(i)).M();
        tmpMt2=(v[1]+electronColl.at(0)+muonColl.at(0)+muonColl.at(1)+jetColl.at(i)).M();

        tmpJ.SetPxPyPzE(0,0,0,0);
        for(int j=0;j<4;j++){ if(j!=i) tmpJ+=jetColl.at(j); } //Mass from jets from hadronically decaying tx
      
        tmpMtot1=fabs(tmpMt1-172)+fabs(tmpJ.M()-172);
        tmpMtot2=fabs(tmpMt2-172)+fabs(tmpJ.M()-172);
        if(tmpMtot1<tmpMtot2) {tmpMtot=tmpMtot1; vti=0;}
        else {tmpMtot=tmpMtot2; vti=1;}

        if(i==0) {deltatlv=tmpMtot; bti=0;}
        else if(tmpMtot<deltatlv) {deltatlv=tmpMtot; bti=i;}
     } 

     //tHc side qq decay, tx side lv decay scenario
     for(int i=0; i<4; i++){
        tmpMtx1=(v[0]+electronColl.at(0)+jetColl.at(i)).M();
        tmpMtx2=(v[1]+electronColl.at(0)+jetColl.at(i)).M();

        tmpJ.SetPxPyPzE(0,0,0,0);
        for(int j=0;j<4;j++){ if(j!=i) tmpJ+=jetColl.at(j); }
        tmpMt=(muonColl.at(0)+muonColl.at(1)+tmpJ).M();

        tmpMtot1=fabs(tmpMtx1-172)+fabs(tmpMt-172);
        tmpMtot2=fabs(tmpMtx2-172)+fabs(tmpMt-172);
        if(tmpMtot1<tmpMtot2) {tmpMtot=tmpMtot1; vtxi=0;}
        else {tmpMtot=tmpMtot2; vtxi=1;}

        if(i==0) {deltatxlv=tmpMtot; btxi=0;}
        else if(tmpMtot<deltatxlv) {deltatxlv=tmpMtot; btxi=i;}
     }
     if(deltatlv<deltatxlv) {Decaytlv=1; Decaytxlv=0;}
     else {Decaytlv=0; Decaytxlv=1;}

     
     int k=0, nb=0;
     int j1_Wi=-1, j2_Wi=-1;
     if(Decaytlv) {btxi=-1; btxi2=-1;} if(Decaytxlv) {bti=-1; bti2=-1;}//by initializing to 4, we know if it is less than 4, if st. exec.
     if(NBJet(jetColl)<3){
       for(int i=0; i<4; i++){//If the remaining 3 jets have btagged jet, then the other 2 must be originated from W
          if((Decaytlv)&&(i!=bti)){
             if((jetColl.at(i).BJetTaggerValue(KJet::CSVv2)>=0.89)&&(nb==1)) {btxi2=i; nb+=1;}
             if((jetColl.at(i).BJetTaggerValue(KJet::CSVv2)>=0.89)&&(nb==0)) {btxi=i; nb+=1;}
             if((jetColl.at(i).BJetTaggerValue(KJet::CSVv2)<0.89)&&(k==1)) {j2_Wi=i; k+=1;}
             if((jetColl.at(i).BJetTaggerValue(KJet::CSVv2)<0.89)&&(k==0)) {j1_Wi=i; k+=1;}
          }       
          if((Decaytxlv)&&(i!=btxi)){
             if((jetColl.at(i).BJetTaggerValue(KJet::CSVv2)>=0.89)&&(nb==1)) {bti2=i; nb+=1;}
             if((jetColl.at(i).BJetTaggerValue(KJet::CSVv2)>=0.89)&&(nb==0)) {bti=i; nb+=1;}
             if((jetColl.at(i).BJetTaggerValue(KJet::CSVv2)<0.89)&&(k==1)) {j2_Wi=i; k+=1;}
             if((jetColl.at(i).BJetTaggerValue(KJet::CSVv2)<0.89)&&(k==0)) {j1_Wi=i; k+=1;}
          }
       }    
    }
     
     if((btxi==-1)||(bti==-1)){//In case none of the 3 jets doesn't include b tagged j at certain level
        double tmpMjj_W; double delta=1000;
        for(int i=0; i<4; i++){
           if((Decaytlv)&&(i==bti)) continue; if((Decaytxlv)&&(i==btxi)) continue;
           for(int j=i+1; j<4; j++){
               if((Decaytlv)&&(j==bti)) continue; if((Decaytxlv)&&(j==btxi)) continue;
               tmpMjj_W=(jetColl.at(i)+jetColl.at(j)).M();
               if(fabs(tmpMjj_W-80.4)<delta) {delta=fabs(tmpMjj_W-80.4); j1_Wi=i; j2_Wi=j; if(Decaytlv) {btxi=6-bti-i-j;} if(Decaytxlv) {bti=6-btxi-i-j;}}
           }
       }
     }
     if(NBJet(jetColl)<3){
        double tmpMjj1=0, tmpMjj2=0;
        if((Decaytlv)&&(nb==2)){
           tmpMjj1=(jetColl.at(j1_Wi)+jetColl.at(btxi)).M(); tmpMjj2=(jetColl.at(j1_Wi)+jetColl.at(btxi2)).M();
           if(fabs(tmpMjj1-80.4)<fabs(tmpMjj2-80.4)) {j2_Wi=btxi;btxi=btxi2; btxi2=-1;}
           else {j2_Wi=btxi2; btxi2=-1;}
        }
        if((Decaytxlv)&&(nb==2)){
           tmpMjj1=(jetColl.at(j1_Wi)+jetColl.at(bti)).M(); tmpMjj2=(jetColl.at(j1_Wi)+jetColl.at(bti2)).M();
           if(fabs(tmpMjj1-80.4)<fabs(tmpMjj2-80.4)) {j2_Wi=bti; bti=bti2; bti2=-1;}
           else {j2_Wi=bti2; bti2=-1;}
        }
     }
        

     if(Decaytlv){
     //emu channel
     FillHist("RecoMhc_emu", (v[vti]+electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0, 200, 200);
     FillHist("RecoMt_emu",  (v[vti]+electronColl.at(0)+muonColl.at(0)+muonColl.at(1)+jetColl.at(bti)).M(), weight, 0, 300, 300);

     FillHist("RecoMW_jj_emu",(jetColl.at(j1_Wi)+jetColl.at(j2_Wi)).M(), weight, 0, 150, 150);
     FillHist("RecoMtx_emu", (jetColl.at(j1_Wi)+jetColl.at(j2_Wi)+jetColl.at(btxi)).M(), weight, 0, 300, 300);
     FillHist("RecoMA_emu", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0, 150, 150);

     //All channel
     FillHist("RecoMhc", (v[vti]+electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0, 200, 200);
     FillHist("RecoMt",  (v[vti]+electronColl.at(0)+muonColl.at(0)+muonColl.at(1)+jetColl.at(bti)).M(), weight, 0, 300, 300);

     FillHist("RecoMW_jj",(jetColl.at(j1_Wi)+jetColl.at(j2_Wi)).M(), weight, 0, 150, 150);
     FillHist("RecoMtx", (jetColl.at(j1_Wi)+jetColl.at(j2_Wi)+jetColl.at(btxi)).M(), weight, 0, 300, 300);
     FillHist("RecoMA", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0, 150, 150);
     }
     if(Decaytxlv){
     //emu channel
     FillHist("RecoMW_jj_emu", (jetColl.at(j1_Wi)+jetColl.at(j2_Wi)).M(), weight, 0, 150, 150);
     FillHist("RecoMhc_emu", (jetColl.at(j1_Wi)+jetColl.at(j2_Wi)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0, 200, 200);
     FillHist("RecoMt_emu", (jetColl.at(j1_Wi)+jetColl.at(j2_Wi)+jetColl.at(bti)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0, 300, 300);

     FillHist("RecoMtx_emu", (jetColl.at(btxi)+electronColl.at(0)+v[vtxi]).M(), weight, 0, 300, 300);
     FillHist("RecoMA_emu", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0, 150, 150);
     
     //All channel
     FillHist("RecoMW_jj", (jetColl.at(j1_Wi)+jetColl.at(j2_Wi)).M(), weight, 0, 150, 150);
     FillHist("RecoMhc", (jetColl.at(j1_Wi)+jetColl.at(j2_Wi)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0, 200, 200);
     FillHist("RecoMt", (jetColl.at(j1_Wi)+jetColl.at(j2_Wi)+jetColl.at(bti)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0, 300, 300);

     FillHist("RecoMtx", (jetColl.at(btxi)+electronColl.at(0)+v[vtxi]).M(), weight, 0, 300, 300);
     FillHist("RecoMA", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0, 150, 150);
     }
   }


   ////////////////////////////////////////////////////////////////////////////////////////////////////
   ////TriMu Channel///////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////////////////////
   if((electronColl.size()==0)&&(muonColl.size()==3)&&(fabs(SumCharge(muonColl))==1)){
     //For muons, we don't know which muon from which. for that we need to tag charge to make it easy.
     //This part we tag which index corresponds to minus of plus charge
     int OSindex=TriMuChargeIndex(muonColl,"OS");
     int SSindex1=TriMuChargeIndex(muonColl,"SS1");
     int SSindex2=TriMuChargeIndex(muonColl, "SS2");

     //Neutrino Filling
     double Pzv11 = GetvPz(v[0], muonColl.at(SSindex1), 1);
     double Pzv12 = GetvPz(v[0], muonColl.at(SSindex1), 2);
     double Pzv21 = GetvPz(v[0], muonColl.at(SSindex2), 1);
     double Pzv22 = GetvPz(v[0], muonColl.at(SSindex2), 2);
     v[0].SetPz(Pzv11);v[0].SetE(sqrt(Pzv11*Pzv11+met*met)); //FillHist("MW", (v[0]+muonColl.at(SSindex1)).M(), weight, 0, 200, 200);
     v[1].SetPxPyPzE(met_x, met_y, Pzv12, sqrt(Pzv12*Pzv12+met*met)); //FillHist("MW", (v[1]+muonColl.at(SSindex1)).M(), weight, 0, 200, 200);
     v[2].SetPxPyPzE(met_x, met_y, Pzv21, sqrt(Pzv21*Pzv21+met*met)); //FillHist("MW", (v[2]+muonColl.at(SSindex2)).M(), weight, 0, 200, 200);
     v[3].SetPxPyPzE(met_x, met_y, Pzv22, sqrt(Pzv22*Pzv22+met*met)); //FillHist("MW", (v[3]+muonColl.at(SSindex2)).M(), weight, 0, 200, 200);


     double tmpMt[4], tmpMtx[4], tmpMtot, tmpMtot1;
     int n=0, tmpindex=-1; 
     int vti=-1, vtxi=-1, bti=-1, bti2=-1, btxi=-1, btxi2=-1, muWindex=-1, muAindex1=-1, muAindex2=-1;
     int Decaytlv=0, Decaytxlv=0;
     double deltatlv=0, deltatxlv=0;
     snu::KParticle tmpJ;  

     //THc side lv decay, tx side qqdecay scenario
     for(int i=0; i<4; i++){

        tmpJ.SetPxPyPzE(0,0,0,0);
        for(int j=0;j<4;j++) {if(j!=i) tmpJ+=jetColl.at(j);} //Mass from jets from hadronically decaying tx

        //Sorting out which direction Pzv and which muon from W is likely to be
        for(int j=0;j<4; j++){
           tmpMt[j]=(v[j]+muonColl.at(0)+muonColl.at(1)+muonColl.at(2)+jetColl.at(i)).M();
           tmpMtot1=fabs(tmpMt[j]-172)+fabs(tmpJ.M()-172);
           if(j==0) {tmpMtot=tmpMtot1; n=j;}
           if(tmpMtot1<tmpMtot) {tmpMtot=tmpMtot1; n=j;}
        }
        if((n==0)||(n==1)) {muWindex=SSindex1; muAindex1=SSindex2; muAindex2=OSindex; vti=n;}
        else if((n==2)||(n==3)) {muWindex=SSindex2; muAindex1=SSindex1; muAindex2=OSindex; vti=n;}
        
        //Sort out which jet is likely from b ftom t
        if(i==0) {deltatlv=tmpMtot; bti=0;}
        else if(tmpMtot<deltatlv) {deltatlv=tmpMtot; bti=i;}
     } 

     //tHc side qq decay, tx side lv decay scenario
     for(int i=0; i<4; i++){

        tmpJ.SetPxPyPzE(0,0,0,0);
        for(int j=0;j<4;j++){ if(j!=i) tmpJ+=jetColl.at(j); }

        for(int j=0; j<4; j++){
            if((j==0)||(j==1)) tmpindex=SSindex1;
            if((j==2)||(j==3)) tmpindex=SSindex2;
            tmpMtx[j]=(v[j]+muonColl.at(tmpindex)+jetColl.at(i)).M();
            tmpMtot1=fabs(tmpMtx[j]-172)+fabs(tmpJ.M()-172);
            if(j==0) {tmpMtot=tmpMtot1; n=j;}
            if(tmpMtot1<tmpMtot) {tmpMtot=tmpMtot1; n=j;}
        }
        if((n==0)||(n==1)) {muWindex=SSindex1; muAindex1=SSindex2; muAindex2=OSindex; vtxi=n;}
        else if((n==2)||(n==3)) {muWindex=SSindex2; muAindex1=SSindex1; muAindex2=OSindex; vtxi=n;}

        if(i==0) {deltatxlv=tmpMtot; btxi=0;}
        else if(tmpMtot<deltatxlv) {deltatxlv=tmpMtot; btxi=i;}
     }
     if(deltatlv<deltatxlv) {Decaytlv=1; Decaytxlv=0;}
     else {Decaytlv=0; Decaytxlv=1;}

     
     int k=0, nb=0;
     int j1_Wi=-1, j2_Wi=-1;
     if(Decaytlv) {btxi=-1; btxi2=-1;} if(Decaytxlv) {bti=-1; bti2=-1;}//by initializing to 4, we know if it is less than 4, if st. exec.
     //If the remaining 3 jets have btagged jet, then the other 2 must be originated from W
     //Purpose : check how many btaggedj is in jetColl & which jet index is b or light
     if(NBJet(jetColl)<3){
       for(int i=0; i<4; i++){
          if((Decaytlv)&&(i!=bti)){
             if((jetColl.at(i).BJetTaggerValue(KJet::CSVv2)>=0.89)&&(nb==1)) {btxi2=i; nb+=1;}
             if((jetColl.at(i).BJetTaggerValue(KJet::CSVv2)>=0.89)&&(nb==0)) {btxi=i; nb+=1;}
             if((jetColl.at(i).BJetTaggerValue(KJet::CSVv2)<0.89)&&(k==1)) {j2_Wi=i; k+=1;}
             if((jetColl.at(i).BJetTaggerValue(KJet::CSVv2)<0.89)&&(k==0)) {j1_Wi=i; k+=1;}
          }       
          if((Decaytxlv)&&(i!=btxi)){
             if((jetColl.at(i).BJetTaggerValue(KJet::CSVv2)>=0.89)&&(nb==1)) {bti2=i; nb+=1;}
             if((jetColl.at(i).BJetTaggerValue(KJet::CSVv2)>=0.89)&&(nb==0)) {bti=i; nb+=1;}
             if((jetColl.at(i).BJetTaggerValue(KJet::CSVv2)<0.89)&&(k==1)) {j2_Wi=i; k+=1;}
             if((jetColl.at(i).BJetTaggerValue(KJet::CSVv2)<0.89)&&(k==0)) {j1_Wi=i; k+=1;}
          }
       }    
    }

     //In case none of the 3 jets doesn't include b tagged j at certain level
     //Purpose : If no b tag, Wjj is chosen to Mjj~Mw, the other should be B.
     if((btxi==-1)||(bti==-1)){
        double tmpMjj_W; double delta=1000;
        for(int i=0; i<4; i++){
           if((Decaytlv)&&(i==bti)) continue; if((Decaytxlv)&&(i==btxi)) continue;
           for(int j=i+1; j<4; j++){
               if((Decaytlv)&&(j==bti)) continue; if((Decaytxlv)&&(j==btxi)) continue;
               tmpMjj_W=(jetColl.at(i)+jetColl.at(j)).M();
               if(fabs(tmpMjj_W-80.4)<delta) {delta=fabs(tmpMjj_W-80.4); j1_Wi=i; j2_Wi=j; if(Decaytlv) {btxi=6-bti-i-j;} if(Decaytxlv) {bti=6-btxi-i-j;}}
           }
       }
     }

     //Pupose : if 2 b tag on 1 branch, then one of M(lj)(bj)~mw is light jet. the other remains as b jet
     if(NBJet(jetColl)<3){
        double tmpMjj1=0, tmpMjj2=0;
        if((Decaytlv)&&(nb==2)){
           tmpMjj1=(jetColl.at(j1_Wi)+jetColl.at(btxi)).M(); tmpMjj2=(jetColl.at(j1_Wi)+jetColl.at(btxi2)).M();
           if(fabs(tmpMjj1-80.4)<fabs(tmpMjj2-80.4)) {j2_Wi=btxi;btxi=btxi2; btxi2=-1;}
           else {j2_Wi=btxi2; btxi2=-1;}
        }
        if((Decaytxlv)&&(nb==2)){
           tmpMjj1=(jetColl.at(j1_Wi)+jetColl.at(bti)).M(); tmpMjj2=(jetColl.at(j1_Wi)+jetColl.at(bti2)).M();
           if(fabs(tmpMjj1-80.4)<fabs(tmpMjj2-80.4)) {j2_Wi=bti; bti=bti2; bti2=-1;}
           else {j2_Wi=bti2; bti2=-1;}
        }
     }
        
//     cout<<"Decaytlv: "<<Decaytlv<<" Decaytxlv: "<<Decaytxlv<<" vti: "<<vti<<" bti: "<<bti<<" j1_Wi: "<<j1_Wi<<" j2_Wi: "<<j2_Wi<<" muAindex1: "<<muAindex1<<" muAindex2: "<<muAindex2<<" muWindex: "<<muWindex<<" btxi: "<<btxi<<" OSindex: "<<OSindex<<" SSindex1: "<<SSindex1<<" SSindex2: "<<SSindex2<<endl;
     if(Decaytlv){
     //3mu channel
     FillHist("RecoMhc_3mu", (v[vti]+muonColl.at(0)+muonColl.at(1)+muonColl.at(2)).M(), weight, 0, 200, 200);
     FillHist("RecoMt_3mu",  (v[vti]+muonColl.at(0)+muonColl.at(1)+muonColl.at(2)+jetColl.at(bti)).M(), weight, 0, 300, 300);

     FillHist("RecoMW_jj_3mu",(jetColl.at(j1_Wi)+jetColl.at(j2_Wi)).M(), weight, 0, 150, 150);
     FillHist("RecoMtx_3mu", (jetColl.at(j1_Wi)+jetColl.at(j2_Wi)+jetColl.at(btxi)).M(), weight, 0, 300, 300);
     FillHist("RecoMA_3mu", (muonColl.at(muAindex1)+muonColl.at(muAindex2)).M(), weight, 0, 150, 150);
     
     //All channel
     FillHist("RecoMhc", (v[vti]+muonColl.at(0)+muonColl.at(1)+muonColl.at(2)).M(), weight, 0, 200, 200);
     FillHist("RecoMt",  (v[vti]+muonColl.at(0)+muonColl.at(1)+muonColl.at(2)+jetColl.at(bti)).M(), weight, 0, 300, 300);

     FillHist("RecoMW_jj",(jetColl.at(j1_Wi)+jetColl.at(j2_Wi)).M(), weight, 0, 150, 150);
     FillHist("RecoMtx", (jetColl.at(j1_Wi)+jetColl.at(j2_Wi)+jetColl.at(btxi)).M(), weight, 0, 300, 300);
     FillHist("RecoMA", (muonColl.at(muAindex1)+muonColl.at(muAindex2)).M(), weight, 0, 150, 150);
     }
     if(Decaytxlv){
     //3mu channel
     FillHist("RecoMW_jj_3mu", (jetColl.at(j1_Wi)+jetColl.at(j2_Wi)).M(), weight, 0, 150, 150);
     FillHist("RecoMhc_3mu", (jetColl.at(j1_Wi)+jetColl.at(j2_Wi)+muonColl.at(muAindex1)+muonColl.at(muAindex2)).M(), weight, 0, 200, 200);
     FillHist("RecoMt_3mu", (jetColl.at(j1_Wi)+jetColl.at(j2_Wi)+jetColl.at(bti)+muonColl.at(muAindex1)+muonColl.at(muAindex2)).M(), weight, 0, 300, 300);

     FillHist("RecoMtx_3mu", (jetColl.at(btxi)+muonColl.at(muWindex)+v[vtxi]).M(), weight, 0, 300, 300);
     FillHist("RecoMA_3mu", (muonColl.at(muAindex1)+muonColl.at(muAindex2)).M(), weight, 0, 150, 150);

     //All channel
     FillHist("RecoMW_jj", (jetColl.at(j1_Wi)+jetColl.at(j2_Wi)).M(), weight, 0, 150, 150);
     FillHist("RecoMhc", (jetColl.at(j1_Wi)+jetColl.at(j2_Wi)+muonColl.at(muAindex1)+muonColl.at(muAindex2)).M(), weight, 0, 200, 200);
     FillHist("RecoMt", (jetColl.at(j1_Wi)+jetColl.at(j2_Wi)+jetColl.at(bti)+muonColl.at(muAindex1)+muonColl.at(muAindex2)).M(), weight, 0, 300, 300);

     FillHist("RecoMtx", (jetColl.at(btxi)+muonColl.at(muWindex)+v[vtxi]).M(), weight, 0, 300, 300);
     FillHist("RecoMA", (muonColl.at(muAindex1)+muonColl.at(muAindex2)).M(), weight, 0, 150, 150);
     }
   }
*/
return;
}// End of execute event loop
  


void Sep2016_3l4j_DataMCComp::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Sep2016_3l4j_DataMCComp::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  string analysisdir = getenv("FILEDIR");  
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

Sep2016_3l4j_DataMCComp::~Sep2016_3l4j_DataMCComp() {
  
  Message("In Sep2016_3l4j_DataMCComp Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}


void Sep2016_3l4j_DataMCComp::FillCutFlow(TString cut, float weight){
  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 12, 0., 12.);
    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"NlCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"OSmuon");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"NjCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(8,"NbCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(9,"NljCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(10,"lPtCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(11,"MjjCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(12,"MmumuCut");
    
  }
}


void Sep2016_3l4j_DataMCComp::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void Sep2016_3l4j_DataMCComp::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //Basic Hists/////////////////////////////////////////////////
  //Data properties before selection  
  AnalyzerCore::MakeHistograms("Basic_Nj_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Nb_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NmuL_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NeL_orig", 10, 0., 10.);
//  AnalyzerCore::MakeHistograms("Basic_NmuT_orig",10, 0., 10.);
//  AnalyzerCore::MakeHistograms("Basic_NeT_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_METdist_orig", 100, 0., 200.);

  //After Nlcut
  AnalyzerCore::MakeHistograms("Basic_NmuL2eL1", 1, 0., 1.);
  AnalyzerCore::MakeHistograms("Basic_NmuL3eL0", 1, 0., 1.);
  AnalyzerCore::MakeHistograms("Basic_NmuL2eL1_wOS", 1, 0., 1.);
  AnalyzerCore::MakeHistograms("Basic_NmuL3eL0_wOS", 1, 0., 1.);

  AnalyzerCore::MakeHistograms("Basic_Nj_wNlcut", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Nb_wNlcut", 10, 0., 10.);
/*
  AnalyzerCore::MakeHistograms("Basic_NmuT2eL1", 1, 0., 1.);
  AnalyzerCore::MakeHistograms("Basic_NmuL2eT1", 1, 0., 1.);
  AnalyzerCore::MakeHistograms("Basic_NmuT2eT1", 1, 0., 1.);
  AnalyzerCore::MakeHistograms("Basic_NmuL2muT1eL0", 1, 0., 1.);
  AnalyzerCore::MakeHistograms("Basic_NmuL1muT2eL0", 1, 0., 1.);
  AnalyzerCore::MakeHistograms("Basic_NmuT3eL0", 1, 0., 1.);
*/
  //After Nljcut
  AnalyzerCore::MakeHistograms("Basic_Nb_wNljcut", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Pt_LeadLep", 50, 0., 200.);

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

  AnalyzerCore::MakeHistograms("Basic_Etaj1_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etaj2_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etaj3_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etaj4_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etab1_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etab2_wPreSel", 100, -5., 5.);


  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Sep2016_3l4j_DataMCCompCore::MakeHistograms() to make new hists for your analysis
  // **
  
}


void Sep2016_3l4j_DataMCComp::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
