#ifndef Aug2017_ObjEffCalc_h
#define Aug2017_ObjEffCalc_h

#include "AnalyzerCore.h"


class Aug2017_ObjEffCalc : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  Aug2017_ObjEffCalc();
  ~Aug2017_ObjEffCalc();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);

  int NPromptFake_Ele(std::vector<snu::KElectron> EleColl, TString Option);
  int NPromptFake_Ele(std::vector<snu::KElectron> EleColl, std::vector<snu::KTruth> TruthColl, TString Option);
  int NPromptFake_Mu(std::vector<snu::KMuon> MuColl, TString Option);
  int NPromptFake_Mu(std::vector<snu::KMuon> MuColl, std::vector<snu::KTruth> TruthColl, TString Option);
  bool IsConvCand(snu::KElectron Ele, std::vector<snu::KTruth> TruthColl, TString Option="");


 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( Aug2017_ObjEffCalc, 1);
};
#endif
