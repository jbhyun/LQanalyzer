#ifndef Sep2016_3l4j_Validation_DiLep_h
#define Sep2016_3l4j_Validation_DiLep_h

#include "AnalyzerCore.h"


class Sep2016_3l4j_Validation_DiLep : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  Sep2016_3l4j_Validation_DiLep();
  ~Sep2016_3l4j_Validation_DiLep();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);
  void FillTrigDist(TString cut, float w);
  void GetGenMatchedSigIndex(std::vector<snu::KTruth>& truthColl, std::vector<snu::KMuon>& muonColl, std::vector<snu::KElectron>& electronColl, std::vector<snu::KJet>& jetColl, int mum_Ai, int mup_Ai, int e_Wi, int mu_Wi, int j1_Wi, int j2_Wi, int b_ti, int bx_txi, float weight);
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( Sep2016_3l4j_Validation_DiLep, 1);
};
#endif
