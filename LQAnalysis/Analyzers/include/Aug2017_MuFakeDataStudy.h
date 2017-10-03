#ifndef Aug2017_MuFakeDataStudy_h
#define Aug2017_MuFakeDataStudy_h

#include "AnalyzerCore.h"


class Aug2017_MuFakeDataStudy : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  Aug2017_MuFakeDataStudy();
  ~Aug2017_MuFakeDataStudy();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);

  int   NPromptFake_Ele(std::vector<snu::KElectron> EleColl, std::vector<snu::KTruth> TruthColl, TString Option);
  int   NPromptFake_Mu(std::vector<snu::KMuon> MuColl, std::vector<snu::KTruth> TruthColl, TString Option);
  bool  IsConvCand(snu::KElectron Ele, std::vector<snu::KTruth> TruthColl, TString Option="");
  int   StepPassed(std::vector<snu::KMuon> MuColl, std::vector<snu::KElectron> EleColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, TString Option="");
  float ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut);
  float ConeCorrectedPT(snu::KMuon Mu, float TightIsoCut);
  int   GetFakeLepJetSrcType(snu::KMuon Mu, std::vector<snu::KJet> JetColl);
  float FakeRateData(snu::KMuon Mu, TString Option);
  float GetFakeWeight(std::vector<snu::KMuon> MuLColl, TString MuLID, TString MuTID, TString Option);

  float NvtxWeight(int Nvtx, TString Option);
  bool IsNearBJet(snu::KMuon Mu, std::vector<snu::KJet> bjetNoVetoColl);


  void DoSystRun(TString Cycle, TString Mode, std::vector<snu::KElectron> EleColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight);

  void ScanIPEfficiency(std::vector<snu::KMuon> muonColl, std::vector<snu::KElectron> electronLooseColl, float IsoCut, float weight, TString Label, TString Option);
  void CheckIDVarDistribution(std::vector<snu::KMuon> muonColl, std::vector<snu::KElectron> electronLooseColl, float weight, TString Label, TString Option);
  void CheckNormCR(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, float MET, float MET_x, float MET_y, float weight, TString Label, TString Option);
  void MeasureFakeRate(std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString TightID, TString Label, TString Option);
  void ValidateID(std::vector<snu::KMuon> MuTColl, float weight, TString Label);
  void CheckTrilepCRs(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option);

  int GetFakeLepSrcIdx(snu::KMuon Mu, std::vector<snu::KTruth> TruthColl);
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( Aug2017_MuFakeDataStudy, 1);
};
#endif
