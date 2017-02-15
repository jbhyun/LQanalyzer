// $Id: Sep2016_MuMubj_DataMCComp.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQSep2016_MuMubj_DataMCComp Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Sep2016_MuMubj_DataMCComp.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Sep2016_MuMubj_DataMCComp);

 Sep2016_MuMubj_DataMCComp::Sep2016_MuMubj_DataMCComp() : AnalyzerCore(), out_muons(0) {

   SetLogName("Sep2016_MuMubj_DataMCComp");
   Message("In Sep2016_MuMubj_DataMCComp constructor", INFO);
   InitialiseAnalysis();
 }


 void Sep2016_MuMubj_DataMCComp::InitialiseAnalysis() throw( LQError ) {
   
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

void Sep2016_MuMubj_DataMCComp::ExecuteEvents()throw( LQError ){

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



   //Trigers
   std::vector<TString> triggerslist, triggerlist1, triggerlist2;//  ListTriggersAvailable();

   TString analysis_trigger;
   analysis_trigger="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v";
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

     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(20.);                    eventbase->GetMuonSel()->SetEta(2.1);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.25);
   std::vector<snu::KMuon> muonVetoColl; eventbase->GetMuonSel()->Selection(muonVetoColl);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(25.);                    eventbase->GetMuonSel()->SetEta(2.1);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.25);
   //eventbase->GetMuonSel()->SetBSdxy(0.005);               eventbase->GetMuonSel()->SetBSdz(0.03);
   //eventbase->GetMuonSel()->SetChiNdof(3);
   std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->Selection(muonColl);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
 //std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->SelectMuons(muonColl, BaseSelection::MUON_POG_LOOSE, 10., 2.4);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
   
   //eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
   //eventbase->GetElectronSel()->SetPt(20.);                eventbase->GetElectronSel()->SetEta(2.5);
   //eventbase->GetElectronSel()->SetBETrRegIncl(false);
   //eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//Iso Ichep16 Loose WP
   //eventbase->GetElectronSel()->SetdxyBEMax(0.005, 0.01);  eventbase->GetElectronSel()->SetdzBEMax(0.005, 0.01);
   std::vector<snu::KElectron> electronColl; //eventbase->GetElectronSel()->Selection(electronColl);//pt>10/eta<2.4(hole region ~1.4 excluded by default)/EGAMMA_LOOSE/NoExtIso*/

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_TIGHT);
     eventbase->GetJetSel()->SetPt(30.);                     eventbase->GetJetSel()->SetEta(4.7);
//   eventbase->GetJetSel()->SetUseJetPileUp(true);
     bool LeptonVeto=true;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonVetoColl, electronColl, isData, false);
// std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->SelectJets(isData, jetColl, muonColl, electronColl, "PFJET_TIGHT", 20., 2.4, true);

//   std::vector<int> bIdxMColl=GetSFBJetIdx(jetColl,"Medium");
//   std::vector<int> bIdxTColl=GetSFBJetIdx(jetColl,"Tight");

   //RegionalJets/////
   std::vector<snu::KJet> bjetBarrelMColl, bjetBarrelTColl;
   std::vector<snu::KJet> jetBarrelColl, jetEndCapColl;
     for(unsigned int i=0; i<jetColl.size(); i++){
       if(fabs(jetColl.at(i).Eta())<2.4){
         jetBarrelColl.push_back(jetColl.at(i));
         if(jetColl.at(i).BJetTaggerValue(snu::KJet::CSVv2)>0.935) bjetBarrelTColl.push_back(jetColl.at(i));
         if(jetColl.at(i).BJetTaggerValue(snu::KJet::CSVv2)>0.800) bjetBarrelMColl.push_back(jetColl.at(i));
       }
       else{
         jetEndCapColl.push_back(jetColl.at(i));
       }
     }



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
     trigger_eff = TriggerEff(analysis_trigger, muonColl);
     weight*=trigger_eff;
   }
   FillCutFlow("TriggerCut", weight); //;in trigger Cut mode
   /////////////////////////////////////////////////////////////////////////////////////////


   double met = eventbase->GetEvent().MET();

   //Scale Factors
   float id_weight=1., reco_weight=1., iso_weight=1.;
   if(!isData){
//     id_weight *= ElectronScaleFactor(BaseSelection::ELECTRON_POG_TIGHT, electronColl);
//     id_weight *= MuonScaleFactor(BaseSelection::MUON_POG_TIGHT, muonColl,0);
//     iso_weight *=MuonISOScaleFactor(BaseSelection::MUON_POG_TIGHT, muonColl,0);
//     reco_weight *= ElectronRecoScaleFactor(electronColl);
//     weight*=id_weight*reco_weight*iso_weight;
   }



   //Variables//////////////////////
   int NbBarrelM=bjetBarrelMColl.size();
   int NbBarrelT=bjetBarrelTColl.size();
   int NjBarrel=jetBarrelColl.size();
   int NjEndcap=jetEndCapColl.size();
   int njets=jetColl.size();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///Basic Objects Distribution// 
   FillHist("Basic_Nj_orig", njets, weight, 0., 10., 10);
   FillHist("Basic_NbM_orig", NbBarrelM, weight, 0., 10., 10);//Bjet def CVSincV2>0.89;pfKJet::CSVv2IVF_discriminator 
   FillHist("Basic_NmuT_orig", muonColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeT_orig", electronColl.size(), weight, 0., 10., 10);

   if(muonColl.size()>0) FillHist("Basic_Ptmu1_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
   if(muonColl.size()>1) FillHist("Basic_Ptmu2_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
   if(electronColl.size()>0) FillHist("Basic_Pte1_orig", electronColl.at(0).Pt(), weight, 0., 200., 200);
   if(electronColl.size()>1) FillHist("Basic_Pte2_orig", electronColl.at(1).Pt(), weight, 0., 200., 200);

   if(njets>0) FillHist("Basic_j1_Et_orig", jetColl.at(0).Et(), weight, 0., 200., 200);
   if(njets>1) FillHist("Basic_j2_Et_orig", jetColl.at(1).Et(), weight, 0., 200., 200);
   if(njets>2) FillHist("Basic_j3_Et_orig", jetColl.at(2).Et(), weight, 0., 200., 200);
   if(njets>3) FillHist("Basic_j4_Et_orig", jetColl.at(3).Et(), weight, 0., 200., 200);
   if(NbBarrelM>0) FillHist("Basic_b1_Et_orig", bjetBarrelMColl.at(0).Et(), weight, 0, 200., 200);
   if(NbBarrelM>1) FillHist("Basic_b2_Et_orig", bjetBarrelMColl.at(1).Et(), weight, 0., 200., 200);
   FillHist("Basic_METdist_orig", met, weight, 0., 200., 100);//It is orig because even though I havent put METcut below, if there is correlation between Nl ,Nj, Nb and MET then MET distrubution will also be affected by selection. The distribution remains the same only when the correlation is absent.

///////Event Selection&Histograms//////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////   

   if(muonColl.size()!=2) return;
   FillCutFlow("NlCut", weight);

   if(SumCharge(muonColl)!=0) return;
   FillCutFlow("OSmuon", weight);

   if(NbBarrelM<1) return;
   FillCutFlow("NbM", weight);

   if(njets<2) return; 
   FillCutFlow("Nj", weight);


   FillHist("Basic_Nvtx_wNlOSbjcut", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);
   FillHist("Basic_MET_wNlOSbjcut", met, weight, 0., 200., 100);
   FillHist("Basic_Nj_wNlOSbjcut", njets, weight, 0., 10., 10);//Njets after all cut
   FillHist("Basic_Nb_wNlOSbjcut", NbBarrelM, weight, 0., 5., 5);//Nbjet after all cuts

   FillHist("Basic_mu1_Pt_wNlOSbjcut", muonColl.at(0).Pt(), weight, 0., 200., 200);
   FillHist("Basic_mu1_Eta_wNlOSbjcut", muonColl.at(0).Eta(), weight, -5., 5., 100);
   FillHist("Basic_mu2_Pt_wNlOSbjcut", muonColl.at(1).Pt(), weight, 0., 200., 200);
   FillHist("Basic_mu2_Eta_wNlOSbjcut", muonColl.at(1).Eta(), weight, -5., 5., 100);
   FillHist("Basic_Ptj1_wNlOSbjcut", jetColl.at(0).Pt(), weight, 0., 200., 200);
   FillHist("Basic_Ptj2_wNlOSbjcut", jetColl.at(1).Pt(), weight, 0., 200., 200);
   FillHist("Basic_Ptb1_wNlOSbjcut", bjetBarrelMColl.at(0).Pt(), weight, 0., 200., 200);

   FillHist("Basic_Etaj1_wNlOSbjcut", jetColl.at(0).Eta(), weight, -5., 5., 100);
   FillHist("Basic_Etaj2_wNlOSbjcut", jetColl.at(1).Eta(), weight, -5., 5., 100);
   FillHist("Basic_Etab1_wNlOSbjcut", bjetBarrelMColl.at(0).Eta(), weight, -5., 5., 100);

   //CSV Medium WP
   //Region 1
   if((NbBarrelM==1)&&(NjBarrel==1)&&(NjEndcap>0)){
     FillCutFlow("Reg1_CSVM",weight);
     FillHist("Mmumu_Reg1_CSVM",(muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 100., 200);
   }
   //Region2
   else if((NbBarrelM>0)&&(NjBarrel==2)&&(NjEndcap==0)&&(met<40)){
     if((muonColl.at(0)+muonColl.at(1)).DeltaPhi(jetBarrelColl.at(0)+jetBarrelColl.at(1))>2.5){
       FillCutFlow("Reg2_CSVM",weight);
       FillHist("Mmumu_Reg2_CSVM",(muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 100., 200);
     }
   }

   //CSV Tight WP
   //Region 1
   if((NbBarrelT==1)&&(NjBarrel==1)&&(NjEndcap>0)){
     FillCutFlow("Reg1_CSVT",weight);
     FillHist("Mmumu_Reg1_CSVT",(muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 100., 200);
   }
   //Region2
   else if((NbBarrelT>0)&&(NjBarrel==2)&&(NjEndcap==0)&&(met<40)){
     if((muonColl.at(0)+muonColl.at(1)).DeltaPhi(jetBarrelColl.at(0)+jetBarrelColl.at(1))>2.5){
       FillCutFlow("Reg2_CSVT",weight);
       FillHist("Mmumu_Reg2_CSVT",(muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 100., 200);
     }
   }


//////////////////////////////////////////////////////////////////////////////////////////////


return;
}// End of execute event loop
  


void Sep2016_MuMubj_DataMCComp::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Sep2016_MuMubj_DataMCComp::BeginCycle() throw( LQError ){
  
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

Sep2016_MuMubj_DataMCComp::~Sep2016_MuMubj_DataMCComp() {
  
  Message("In Sep2016_MuMubj_DataMCComp Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}


void Sep2016_MuMubj_DataMCComp::FillCutFlow(TString cut, float weight){
  
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
    GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"NbM");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(8,"Nj");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(9,"Reg1_CSVM");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(10,"Reg2_CSVM");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(11,"Reg1_CSVT");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(12,"Reg2_CSVT");
    
  }
}


void Sep2016_MuMubj_DataMCComp::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Sep2016_MuMubj_DataMCComp::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //Basic Hists/////////////////////////////////////////////////
  //Data properties before selection  
  AnalyzerCore::MakeHistograms("Basic_Nvtx_raw", 50, 0., 50.);
  AnalyzerCore::MakeHistograms("Basic_Nvtx", 50, 0., 50.);
  AnalyzerCore::MakeHistograms("Basic_Nj_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Nb_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NmuT_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NeT_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu1_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu2_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Pte1_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Pte2_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j1_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j2_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j3_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j4_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_b1_Et_orig", 200, 0, 200.);
  AnalyzerCore::MakeHistograms("Basic_b2_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_METdist_orig", 100, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Nvtx_wNlOSbjcut", 50, 0., 50.);
  AnalyzerCore::MakeHistograms("Basic_MET_wNlOSbjcut", 100, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Nj_wNlOSbjcut", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Nb_wNlOSbjcut", 5, 0., 5.);
  AnalyzerCore::MakeHistograms("Basic_mu1_Pt_wNlOSbjcut", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_mu1_Eta_wNlOSbjcut", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_mu2_Pt_wNlOSbjcut", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_mu2_Eta_wNlOSbjcut", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Ptj1_wNlOSbjcut", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptj2_wNlOSbjcut", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptb1_wNlOSbjcut", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Etaj1_wNlOSbjcut", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etaj2_wNlOSbjcut", 100, -5., 5.);
  AnalyzerCore::MakeHistograms("Basic_Etab1_wNlOSbjcut", 100, -5., 5.);

  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Sep2016_MuMubj_DataMCCompCore::MakeHistograms() to make new hists for your analysis
  // **
  
}




void Sep2016_MuMubj_DataMCComp::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
