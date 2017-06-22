#ifndef May2017_MCFakeStudy_h
#define May2017_MCFakeStudy_h

#include "AnalyzerCore.h"


class May2017_MCFakeStudy : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  May2017_MCFakeStudy();
  ~May2017_MCFakeStudy();

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


  ClassDef ( May2017_MCFakeStudy, 1);
};
#endif
