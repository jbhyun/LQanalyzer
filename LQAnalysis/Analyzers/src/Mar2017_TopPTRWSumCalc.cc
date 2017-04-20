// $Id: Mar2017_TopPTRWSumCalc.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQMar2017_TopPTRWSumCalc Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Mar2017_TopPTRWSumCalc.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Mar2017_TopPTRWSumCalc);

 Mar2017_TopPTRWSumCalc::Mar2017_TopPTRWSumCalc() : AnalyzerCore(), out_muons(0) {

   SetLogName("Mar2017_TopPTRWSumCalc");
   Message("In Mar2017_TopPTRWSumCalc constructor", INFO);
   InitialiseAnalysis();
 }


 void Mar2017_TopPTRWSumCalc::InitialiseAnalysis() throw( LQError ) {
   
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

void Mar2017_TopPTRWSumCalc::ExecuteEvents()throw( LQError ){


////Basic Infos///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Apply the gen weight 
   FillHist("GenWeight", MCweight, 1, -2., 2., 4);

   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
   double TopPTreweight=TopPTReweight(truthColl);

 
   FillHist("SumGenW", 0., MCweight, 0., 1., 1);
   FillHist("SumTopPTRW", 0., MCweight*TopPTreweight, 0., 1., 1);

   int Nhardt=0;
   for(std::vector<snu::KTruth>::iterator it_truth= truthColl.begin(); it_truth!=truthColl.end(); it_truth++){
     if( fabs(it_truth->PdgId())==6 && it_truth->GenStatus()<30 && it_truth->GenStatus()>20 ){
       Nhardt++;
       FillHist("NhardTopGenStDist", it_truth->GenStatus(), 1., 0., 100., 100);
     }
     if( fabs(it_truth->PdgId())==6 ){
       FillHist("NTopGenStDist", it_truth->GenStatus(), 1., 0., 100., 100);
     }
   }
   FillHist("NHardScatterTopCount", Nhardt, 1., 0., 10., 10);



/////////////////////////////////////////////////////////////////////////////////// 


return;
}// End of execute event loop
  


void Mar2017_TopPTRWSumCalc::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Mar2017_TopPTRWSumCalc::BeginCycle() throw( LQError ){
  
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

Mar2017_TopPTRWSumCalc::~Mar2017_TopPTRWSumCalc() {
  
  Message("In Mar2017_TopPTRWSumCalc Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Mar2017_TopPTRWSumCalc::FillCutFlow(TString cut, float weight){
  
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



void Mar2017_TopPTRWSumCalc::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Mar2017_TopPTRWSumCalc::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //Basic Hists/////////////////////////////////////////////////
  //Data properties before selection  
  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Mar2017_TopPTRWSumCalcCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Mar2017_TopPTRWSumCalc::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
