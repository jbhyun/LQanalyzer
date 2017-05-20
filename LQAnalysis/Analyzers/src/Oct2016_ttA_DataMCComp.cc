// $Id: Oct2016_ttA_DataMCComp.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQOct2016_ttA_DataMCComp Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Oct2016_ttA_DataMCComp.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Oct2016_ttA_DataMCComp);

 Oct2016_ttA_DataMCComp::Oct2016_ttA_DataMCComp() : AnalyzerCore(), out_muons(0) {

   SetLogName("Oct2016_ttA_DataMCComp");
   Message("In Oct2016_ttA_DataMCComp constructor", INFO);
   InitialiseAnalysis();
 }


 void Oct2016_ttA_DataMCComp::InitialiseAnalysis() throw( LQError ) {
   
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

void Oct2016_ttA_DataMCComp::ExecuteEvents()throw( LQError ){

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
   if (!k_isdata) { weight*=pileup_reweight;}

   //Numbet of Vertex NoPUreweight
   FillHist("Nvtx_nocut_PURW", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);


   bool emumu_analysis=false, trimu_analysis=true;

   //Trigers
   std::vector<TString> triggerslist, triggerlist1, triggerlist2;//  ListTriggersAvailable();

   TString analysis_trigger;
   if(emumu_analysis) {analysis_trigger="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";}
   //if(emumu_analysis) {analysis_trigger="HLT_Ele32_eta2p1_WPTight_Gsf";}
   //if(emumu_analysis) {analysis_trigger="HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL";}//"HLT_Ele27_eta2p1_WPLoose_Gsf";}
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
   //if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; FillCutFlow("VertexCut", weight);
   //Good Primary vtx def:(vtx.ndof()>4&&maxAbsZ<=0)||std::abs(vtx.z())<= 24)&&((maxd0 <=0)||std::abs(vtx.position().rho())<=2)&&!(vtx.isFake()))  



///////Objects in Analysis/////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

   //Primary Object Collection
   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
   
   //MuonID////////////////////////////////////////////////////////////////////////////////////////
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.25);
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.25);
//   eventbase->GetMuonSel()->SetBSdxy(0.005);               eventbase->GetMuonSel()->SetBSdz(0.03);
//   eventbase->GetMuonSel()->SetChiNdof(3);
   std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->Selection(muonColl);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
// std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->SelectMuons(muonColl, BaseSelection::MUON_POG_LOOSE, 10., 2.4);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE

   //ElectronID//////////////////////////////////////////////////////////////////////////////////////
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0994, 0.107);//ISO Ichep16 Loose WP
   //eventbase->GetElectronSel()->SetdxyBEMax(0.005, 0.01);  eventbase->GetElectronSel()->SetdzBEMax(0.005, 0.01);
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);//pt>10/eta<2.4(hole region ~1.4 excluded by default)/EGAMMA_LOOSE/NoExtIso*/
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_VETO);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.175, 0.159);//ISO Ichep16 Loose WP
   //eventbase->GetElectronSel()->SetdxyBEMax(0.005, 0.01);  eventbase->GetElectronSel()->SetdzBEMax(0.005, 0.01);
   std::vector<snu::KElectron> electronVetoColl; eventbase->GetElectronSel()->Selection(electronVetoColl);//pt>10/eta<2.4(hole region ~1.4 excluded by default)/EGAMMA_LOOSE/NoExtIso*/

     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MEDIUM);
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.4);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//Iso Ichep16 Loose WP
//   eventbase->GetElectronSel()->SetdxyBEMax(0.005, 0.01);  eventbase->GetElectronSel()->SetdzBEMax(0.005, 0.01);
   std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->Selection(electronColl);//pt>10/eta<2.4(hole region ~1.4 excluded by default)/EGAMMA_LOOSE/NoExtIso*/
// std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->SelectElectrons(electronColl, BaseSelection::ELECTRON_POG_TIGHT, 20., 2.4);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE


    //JetID///////////////////////////////////////////////////////////////////////////////////////////////
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(30.);                     eventbase->GetJetSel()->SetEta(2.4);
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
//
   if((emumu_analysis)&&(electronColl.size()==1)) FillCutFlow("1eCut", weight);
   else return; //veto if not emumu or mumumu

   if(njets>0) FillCutFlow("1jCut", weight);
   else return;

   if(njets>1) FillCutFlow("2jCut", weight);
   else return;

   if(njets>2) FillCutFlow("3jCut", weight);
   else return;

   if(njets>3) FillCutFlow("4jCut", weight);
   else return;

   if(nbjets>0) FillCutFlow("1bCut", weight);
   else return;

   FillHist("Basic_Pte_wl4j1b", electronColl.at(0).Pt(), weight, 0., 200., 200);
   FillHist("Basic_Etae_wl4j1b", electronColl.at(0).Eta(), weight, -5., 5., 100);
   FillHist("Basic_Ptj1_wl4j1b", jetColl.at(0).Pt(), weight, 0., 500., 500);
   FillHist("Basic_Etaj1_wl4j1b", jetColl.at(0).Eta(), weight, -5., 5., 100);
   FillHist("Basic_Ptj2_wl4j1b", jetColl.at(1).Pt(), weight, 0., 500., 500);
   FillHist("Basic_Etaj2_wl4j1b", jetColl.at(1).Eta(), weight, -5., 5., 100);
   FillHist("Basic_MET_wl4j1b", met, weight, 0., 200., 200);


   if(muonColl.size()!=2) return;

   if(SumCharge(muonColl)==0) FillCutFlow("OS2mu", weight);
   else return;

   if(muonColl.at(0).Pt()>20) FillCutFlow("l#muPT",weight);
   else return;

   
/////////////////////////////////////////////                  

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
   FillHist("Basic_Ptj4_wPreSel", jetColl.at(0).Pt(), weight, 0., 200., 200);
   FillHist("Basic_Ptb1_wPreSel", bjetColl.at(0).Pt(), weight, 0., 200., 200);
   if(nbjets>1) FillHist("Basic_Ptb2_wPreSel", bjetColl.at(1).Pt(), weight, 0., 200., 200);

   FillHist("Basic_Etaj1_wPreSel", jetColl.at(0).Eta(), weight, -5., 5., 100);
   FillHist("Basic_Etaj2_wPreSel", jetColl.at(1).Eta(), weight, -5., 5., 100);
   FillHist("Basic_Etaj3_wPreSel", jetColl.at(2).Eta(), weight, -5., 5., 100);
   FillHist("Basic_Etaj4_wPreSel", jetColl.at(0).Eta(), weight, -5., 5., 100);
   FillHist("Basic_Etab1_wPreSel", bjetColl.at(0).Eta(), weight, -5., 5., 100);
   if(nbjets>1) FillHist("Basic_Etab2_wPreSel", bjetColl.at(1).Eta(), weight, -5., 5., 100);

//////////////////////////////////////////////////////////////////////////////////////////////   
/////////////Analysis/////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

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
     
     if(muonColl.at(0).RelIso04()<0.15 && muonColl.at(1).RelIso04()<0.15){
       FillHist("Mmumu_IsoT_1e2mu", (muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 200., 200);
     }
   }


return;
}// End of execute event loop
  


void Oct2016_ttA_DataMCComp::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Oct2016_ttA_DataMCComp::BeginCycle() throw( LQError ){
  
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

Oct2016_ttA_DataMCComp::~Oct2016_ttA_DataMCComp() {
  
  Message("In Oct2016_ttA_DataMCComp Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Oct2016_ttA_DataMCComp::FillCutFlow(TString cut, float weight){
  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 12, 0., 12.);
    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"1eCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"1jCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"2jCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(8,"3jCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(9,"4jCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(10,"1bCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(11,"OS2mu");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(12,"l#muPT");
    
  }
}



void Oct2016_ttA_DataMCComp::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}

void Oct2016_ttA_DataMCComp::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //Basic Hists/////////////////////////////////////////////////
  //Data properties before selection  
  AnalyzerCore::MakeHistograms("Basic_prescale", 2000, 0., 2000.);
  AnalyzerCore::MakeHistograms("Basic_Nj_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Nb_orig", 10, 0., 10.);//Bjet def CVSincV2>0.89;pfKJet::CSVv2IVF_discriminator 
  AnalyzerCore::MakeHistograms("Basic_NmuT_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NmuL_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NeT_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NeL_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NmuT_weT", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu1_3mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu2_3mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu3_3mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu1_1e2mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu2_1e2mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Pte1_1e2mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j1_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j2_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j3_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j4_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_b1_Et_orig", 200, 0, 200.);
  AnalyzerCore::MakeHistograms("Basic_b2_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_METdist_orig", 100, 0., 200.);//It is orig because even though I havent put METcut below, if there is correlation between Nl ,Nj, Nb and MET then MET distrubution will also be affected by selection. The distribution remains the same only when the correlation is absent.
  AnalyzerCore::MakeHistograms("Basic_Pte_wl4j1b", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Etae_wl4j1b", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Ptj1_wl4j1b", 500, 0., 500.);
  AnalyzerCore::MakeHistograms("Basic_Etaj1_wl4j1b", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Ptj2_wl4j1b", 500, 0., 500.);
  AnalyzerCore::MakeHistograms("Basic_Etaj2_wl4j1b", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_MET_wl4j1b", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_MET_wPreSel", 100, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Nj_wPreSel", 10, 0., 10.);//Njets after all cut
  AnalyzerCore::MakeHistograms("Basic_Nb_wPreSel", 5, 0., 5.);//Nbjet after all cuts
  AnalyzerCore::MakeHistograms("Basic_Pte_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu1_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu2_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Etae_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etamu1_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etamu2_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu1_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu2_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu3_wPreSel", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Etamu1_wPreSel", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etamu2_wPreSel", 100, -5., 5.);
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
  // *  Remove//Overide this Oct2016_ttA_DataMCCompCore::MakeHistograms() to make new hists for your analysis
  // **
  
}





void Oct2016_ttA_DataMCComp::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}

