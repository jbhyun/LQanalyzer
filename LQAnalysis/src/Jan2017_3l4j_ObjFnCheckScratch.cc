// $Id: Jan2017_3l4j_ObjFnCheckScratch.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQJan2017_3l4j_ObjFnCheckScratch Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Jan2017_3l4j_ObjFnCheckScratch.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Jan2017_3l4j_ObjFnCheckScratch);

 Jan2017_3l4j_ObjFnCheckScratch::Jan2017_3l4j_ObjFnCheckScratch() : AnalyzerCore(), out_muons(0) {

   SetLogName("Jan2017_3l4j_ObjFnCheckScratch");
   Message("In Jan2017_3l4j_ObjFnCheckScratch constructor", INFO);
   InitialiseAnalysis();
 }


 void Jan2017_3l4j_ObjFnCheckScratch::InitialiseAnalysis() throw( LQError ) {
   
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

void Jan2017_3l4j_ObjFnCheckScratch::ExecuteEvents()throw( LQError ){

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
   if (!k_isdata) { pileup_reweight = TempPileupWeight(); weight*=pileup_reweight;}

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
//   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
  
  //MY FUNCTION
//     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
//     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
//     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.25);
//   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl);
//     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
//     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
//     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.15);
//     eventbase->GetMuonSel()->SetBSdxy(999.);                eventbase->GetMuonSel()->SetBSdz(999.);
//     eventbase->GetMuonSel()->SetChiNdof(999.);
//   std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->Selection(muonColl);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
 
   std::vector<snu::KMuon> muonLQTightColl;  eventbase->GetMuonSel()->SelectMuons(muonLQTightColl, "MUON_POG_TIGHT_TEST", 10., 2.4);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
   std::vector<snu::KMuon> muonLQTightColl2; eventbase->GetMuonSel()->SelectMuons(muonLQTightColl2, "MUON_POG_TIGHT", 10., 2.4);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE

  //MY FUNCTION   
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//Iso Ichep16 Loose WP
   /*eventbase->GetElectronSel()->SetdxyBEMax(0.005, 0.01);  eventbase->GetElectronSel()->SetdzBEMax(0.0466, 0.417);*/
   std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->Selection(electronColl);

  //LQ SELECTION FUNCTION
  std::vector<snu::KElectron> electronLQColl; eventbase->GetElectronSel()->SelectElectrons(electronLQColl, "ELECTRON_POG_16TIGHT", 25., 2.5);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
   

  //Counter//////////
  if(electronColl.size()!=electronLQColl.size()) FillHist("DiffCounter_Ele", 0.5, 1, 0., 1., 1);
  if(muonLQTightColl.size()!=muonLQTightColl2.size()) FillHist("DiffCounter_LQMu", 0.5, 1, 0., 1., 1);






return;
}// End of execute event loop
  


void Jan2017_3l4j_ObjFnCheckScratch::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Jan2017_3l4j_ObjFnCheckScratch::BeginCycle() throw( LQError ){
  
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

Jan2017_3l4j_ObjFnCheckScratch::~Jan2017_3l4j_ObjFnCheckScratch() {
  
  Message("In Jan2017_3l4j_ObjFnCheckScratch Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Jan2017_3l4j_ObjFnCheckScratch::FillCutFlow(TString cut, float weight){
  
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



void Jan2017_3l4j_ObjFnCheckScratch::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Jan2017_3l4j_ObjFnCheckScratch::MakeHistograms(){
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
  // *  Remove//Overide this Jan2017_3l4j_ObjFnCheckScratchCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Jan2017_3l4j_ObjFnCheckScratch::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
