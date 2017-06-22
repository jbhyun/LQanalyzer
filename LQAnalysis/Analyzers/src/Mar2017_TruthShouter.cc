// $Id: Mar2017_TruthShouter.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQMar2017_TruthShouter Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Mar2017_TruthShouter.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Mar2017_TruthShouter);

 Mar2017_TruthShouter::Mar2017_TruthShouter() : AnalyzerCore(), out_muons(0) {

   SetLogName("Mar2017_TruthShouter");
   Message("In Mar2017_TruthShouter constructor", INFO);
   InitialiseAnalysis();
 }


 void Mar2017_TruthShouter::InitialiseAnalysis() throw( LQError ) {
   
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

void Mar2017_TruthShouter::ExecuteEvents()throw( LQError ){


   //PrintTruth(); return;


   if(!PassMETFilter()) return;
   if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return;
   if(!isData) weight*=MCweight;
   if(!k_isdata) { weight*=eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central);}


   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     bool LeptonVeto=false;
   std::vector<snu::KMuon> NullMuon; std::vector<snu::KElectron> NullEle;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, NullMuon, NullEle);
   std::vector<snu::KJet> JetCleanColl = SkimJetColl(jetColl, truthColl, "NoPrNoTau");
   //How many jets are unmatched to a parton?
   //How many matched jets have different matching from parton?
   for(int i=0; i<JetCleanColl.size(); i++){
     int JetHadFlav=JetCleanColl.at(i).HadronFlavour();
     int MatchedIdx=GenMatchedIdx(JetCleanColl.at(i),truthColl);
     FillHist("PartonMatchingFr", 0., weight, 0., 3., 3);
     if(MatchedIdx==-1) FillHist("PartonMatchingFr", 1., weight, 0., 3., 3);
     if(MatchedIdx==-2) FillHist("PartonMatchingFr", 2., weight, 0., 3., 3);

     if(MatchedIdx<0) continue;
     //Now only parton matched jets are remained

     FillHist("PartHadMatchConsist", 0., weight, 0., 4., 4);
     if(JetHadFlav==5) FillHist("PartHadMatchConsist_B", 0., weight, 0., 4., 4);
     if(JetHadFlav==4) FillHist("PartHadMatchConsist_C", 0., weight, 0., 4., 4);
     if(JetHadFlav==0) FillHist("PartHadMatchConsist_L", 0., weight, 0., 4., 4);

     if(IsJetConsistentPartonHadronMatch(JetCleanColl.at(i), truthColl, "AllFlav")){
       FillHist("PartHadMatchConsist", 1., weight, 0., 4., 4);
       if(JetHadFlav==5) FillHist("PartHadMatchConsist_B", 1., weight, 0., 4., 4);
       if(JetHadFlav==4) FillHist("PartHadMatchConsist_C", 1., weight, 0., 4., 4);
       if(JetHadFlav==0) FillHist("PartHadMatchConsist_L", 1., weight, 0., 4., 4);
     }
     if(IsJetConsistentPartonHadronMatch(JetCleanColl.at(i), truthColl, "BFlav")){
       FillHist("PartHadMatchConsist", 2., weight, 0., 4., 4);
       if(JetHadFlav==5) FillHist("PartHadMatchConsist_B", 2., weight, 0., 4., 4);
       if(JetHadFlav==4) FillHist("PartHadMatchConsist_C", 2., weight, 0., 4., 4);
       if(JetHadFlav==0) FillHist("PartHadMatchConsist_L", 2., weight, 0., 4., 4);
     }
     if(IsJetConsistentPartonHadronMatch(JetCleanColl.at(i), truthColl, "Heavy")){
       FillHist("PartHadMatchConsist", 3., weight, 0., 4., 4);
       if(JetHadFlav==5) FillHist("PartHadMatchConsist_B", 3., weight, 0., 4., 4);
       if(JetHadFlav==4) FillHist("PartHadMatchConsist_C", 3., weight, 0., 4., 4);
       if(JetHadFlav==0) FillHist("PartHadMatchConsist_L", 3., weight, 0., 4., 4);
     }

   }

/////////////////////////////////////////////////////////////////////////////////// 


return;
}// End of execute event loop
  


void Mar2017_TruthShouter::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Mar2017_TruthShouter::BeginCycle() throw( LQError ){
  
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

Mar2017_TruthShouter::~Mar2017_TruthShouter() {
  
  Message("In Mar2017_TruthShouter Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Mar2017_TruthShouter::FillCutFlow(TString cut, float weight){
  
  if(GetHist("cutflow_W") && GetHist("cutflow_N")){
    GetHist("cutflow_W")->Fill(cut,weight);
    GetHist("cutflow_N")->Fill(cut,1);
  }
  else{
    if(!GetHist("cutflow_W")){
      AnalyzerCore::MakeHistograms("cutflow_W", 6, 0., 6.);
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(1,"NoCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(2,"TriggerCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(3,"EventCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(4,"VertexCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(5,"3lCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(6,"bVeto");
    }
    if(!GetHist("cutflow_N")){
      AnalyzerCore::MakeHistograms("cutflow_N", 6, 0., 6.);
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(1,"NoCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(2,"TriggerCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(3,"EventCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(4,"VertexCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(5,"3lCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(6,"bVeto");
    }
  }
}



void Mar2017_TruthShouter::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Mar2017_TruthShouter::MakeHistograms(){
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


  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Mar2017_TruthShouterCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Mar2017_TruthShouter::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
