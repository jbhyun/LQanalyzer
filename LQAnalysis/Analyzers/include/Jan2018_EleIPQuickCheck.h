#ifndef Jan2018_EleIPQuickCheck_h
#define Jan2018_EleIPQuickCheck_h

#include "AnalyzerCore.h"


class Jan2018_EleIPQuickCheck : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  Jan2018_EleIPQuickCheck();
  ~Jan2018_EleIPQuickCheck();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);
  float FakeRateData(snu::KElectron Ele, TString Option);
  float FakeRateData(snu::KMuon      Mu, TString Option);
  float GetPreTrigPURW(int Nvtx);
  void DoSystRun(TString Cycle, TString Mode, std::vector<snu::KElectron> EleColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight);
  bool IsNearBJet(snu::KElectron Ele, std::vector<snu::KJet> bjetNoVetoColl);
  float ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut);
  float ConeCorrectedPT(snu::KMuon     Mu , float TightIsoCut);

  float GetFakeWeight(std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleLColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, TString Option);


  void CheckIPIsoComp(std::vector<snu::KElectron> EleColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, std::vector<snu::KTruth> TruthColl, float MET, float METx, float METy, float weight, TString Label, TString Option);

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( Jan2018_EleIPQuickCheck, 1);
};
#endif
