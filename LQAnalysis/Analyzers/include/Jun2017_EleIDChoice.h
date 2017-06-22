#ifndef Jun2017_EleIDChoice_h
#define Jun2017_EleIDChoice_h

#include "AnalyzerCore.h"


class Jun2017_EleIDChoice : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  Jun2017_EleIDChoice();
  ~Jun2017_EleIDChoice();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);

  int NPromptFake_Ele(std::vector<snu::KElectron> EleColl, std::vector<snu::KTruth> TruthColl, TString Option);
  int NPromptFake_Mu(std::vector<snu::KMuon> MuColl, std::vector<snu::KTruth> TruthColl, TString Option);
  bool IsConvCand(snu::KElectron Ele, std::vector<snu::KTruth> TruthColl, TString Option="");
  int StepPassed(std::vector<snu::KMuon> MuColl, std::vector<snu::KElectron> EleColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, TString Option="");


 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( Jun2017_EleIDChoice, 1);
};
#endif
