#ifndef Jul2017_NvtxWeightCheck_h
#define Jul2017_NvtxWeightCheck_h

#include "AnalyzerCore.h"


class Jul2017_NvtxWeightCheck : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  Jul2017_NvtxWeightCheck();
  ~Jul2017_NvtxWeightCheck();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( Jul2017_NvtxWeightCheck, 1);
};
#endif
