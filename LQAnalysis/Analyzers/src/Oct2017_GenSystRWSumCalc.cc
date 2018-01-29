// $Id: Oct2017_GenSystRWSumCalc.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQOct2017_GenSystRWSumCalc Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Oct2017_GenSystRWSumCalc.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Oct2017_GenSystRWSumCalc);

 Oct2017_GenSystRWSumCalc::Oct2017_GenSystRWSumCalc() : AnalyzerCore(), out_muons(0) {

   SetLogName("Oct2017_GenSystRWSumCalc");
   Message("In Oct2017_GenSystRWSumCalc constructor", INFO);
   InitialiseAnalysis();
 }


 void Oct2017_GenSystRWSumCalc::InitialiseAnalysis() throw( LQError ) {
   
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

void Oct2017_GenSystRWSumCalc::ExecuteEvents()throw( LQError ){


////Basic Infos///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Apply the gen weight 
   bool EMuMu=false, TriMu=false;
   bool GenWSum=false;
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("GenWSum"))       GenWSum    = true;
     //else if(k_flags.at(i).Contains("TriMu"))       TriMu      = true;
   }

   FillHist("GenWeight", MCweight, 1, -2., 2., 4);

   if(GenWSum){
     
     std::vector<float> PdfWVec  = eventbase->GetEvent().PdfWeights();
     std::vector<float> ScaleWVec= eventbase->GetEvent().ScaleWeights();


     FillHist("TotalW_muRmuF", 0.0001, 1, 0., 7., 7);
     FillHist("TotalW_PDF", 0.0001, 1, 0., 103., 103);
     for(int i=0; i<ScaleWVec.size(); i++){
       //cout<<i+0.0001<<" "<<ScaleWVec.at(i)<<endl;
       FillHist("TotalW_muRmuF", i+1.0001, ScaleWVec.at(i)/PdfWVec.at(0), 0., 7., 7);
     }
     for(int i=0; i<PdfWVec.size(); i++){
       //cout<<i+0.0001<<" "<<PdfWVec.at(i)<<endl;
       FillHist("TotalW_PDF", i+1.0001, PdfWVec.at(i)/PdfWVec.at(0), 0., 103., 103);
     }
     //cout<<"LHEWeight "<<eventbase->GetEvent().LHEWeight()<<" ";
     //cout<<"Id1 "<<eventbase->GetEvent().Id1()<<" ";
     //cout<<"Id2 "<<eventbase->GetEvent().Id2()<<" ";
     //cout<<"Q "<<eventbase->GetEvent().Q()<<" ";
     //cout<<"x1 "<<eventbase->GetEvent().x1()<<" ";
     //cout<<"x2 "<<eventbase->GetEvent().x2()<<" ";
     //cout<<"Scale1 "<<eventbase->GetEvent().ScaleWeights().size()<<" ";
     //cout<<"Pdf1 "<<eventbase->GetEvent().PdfWeights().size()<<" ";
     //if(eventbase->GetEvent().ScaleWeights().size()>0) cout<<"Scale1 "<<eventbase->GetEvent().ScaleWeights().at(0)<<" ";
     //if(eventbase->GetEvent().PdfWeights().size()>0) cout<<"Pdf1 "<<eventbase->GetEvent().PdfWeights().at(0)<<" ";
     //cout<<endl;


     return;
   }



/////////////////////////////////////////////////////////////////////////////////// 


return;
}// End of execute event loop
  


void Oct2017_GenSystRWSumCalc::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Oct2017_GenSystRWSumCalc::BeginCycle() throw( LQError ){
  
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

Oct2017_GenSystRWSumCalc::~Oct2017_GenSystRWSumCalc() {
  
  Message("In Oct2017_GenSystRWSumCalc Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Oct2017_GenSystRWSumCalc::FillCutFlow(TString cut, float weight){
  
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



void Oct2017_GenSystRWSumCalc::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Oct2017_GenSystRWSumCalc::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //Basic Hists/////////////////////////////////////////////////
  //Data properties before selection  
  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Oct2017_GenSystRWSumCalcCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Oct2017_GenSystRWSumCalc::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
