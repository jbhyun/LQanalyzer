#ifndef Aug2017_EMuMuTrigCheck_h
#define Aug2017_EMuMuTrigCheck_h

#include "AnalyzerCore.h"


class Aug2017_EMuMuTrigCheck : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  Aug2017_EMuMuTrigCheck();
  ~Aug2017_EMuMuTrigCheck();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);
  float FakeRateData(snu::KElectron Ele, TString ID);
  float GetPreTrigPURW(int Nvtx);
  void DoSystRun(TString Cycle, TString Mode, std::vector<snu::KElectron> EleColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, float MET, float METx, float METy, float weight=1.);
  bool IsNearBJet(snu::KElectron Ele, std::vector<snu::KJet> bjetNoVetoColl);
  void ScanFakeRate(std::vector<snu::KElectron> EleColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, int NMVACuts, float MVACuts[], int NIsoCuts, float IsoCuts[], int NPtEdges, float PtEdges[], TString Option);
  float ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut);

  float GetFakeWeight(std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleLColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID);

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( Aug2017_EMuMuTrigCheck, 1);
};
#endif
