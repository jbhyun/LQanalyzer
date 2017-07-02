#ifndef Jun2017_ConversionStudy_h
#define Jun2017_ConversionStudy_h

#include "AnalyzerCore.h"


class Jun2017_ConversionStudy : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  Jun2017_ConversionStudy();
  ~Jun2017_ConversionStudy();

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


//  bool IsFinalPhotonSt23(std::vector<snu::KTruth> TruthColl);
//  int GetPhotonType(int PhotonIdx, std::vector<snu::KTruth> TruthColl);
//  bool IsHardPhotonConverted(std::vector<snu::KTruth> TruthColl);
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( Jun2017_ConversionStudy, 1);
};
#endif
