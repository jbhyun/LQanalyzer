#ifndef Aug2017_TriLepSR_h
#define Aug2017_TriLepSR_h

#include "AnalyzerCore.h"


class Aug2017_TriLepSR : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  Aug2017_TriLepSR();
  ~Aug2017_TriLepSR();

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
  void CheckSRYield(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option);

  void CheckSRDist(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option);
  bool IsNearBJet(snu::KElectron Ele, std::vector<snu::KJet> bjetNoVetoColl);
  bool IsNearBJet(snu::KMuon Mu, std::vector<snu::KJet>& bjetNoVetoColl);
  void ScanFakeRate(std::vector<snu::KElectron> EleColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, int NMVACuts, float MVACuts[], int NIsoCuts, float IsoCuts[], int NPtEdges, float PtEdges[], TString Option);
  float ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut);
  float ConeCorrectedPT(snu::KMuon     Mu , float TightIsoCut);

  float GetFakeWeight(std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleLColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, TString Option);
  float GetFakeWeight(std::vector<snu::KMuon>& MuLColl, std::vector<snu::KElectron>& EleLColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, vector<snu::KJet>& BJetNoVetoColl, TString Option);


  void OptimizeMACut(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option);
  void CheckSigBkgKinematics(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option);
  void CheckCutflow(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option);

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( Aug2017_TriLepSR, 1);
};
#endif
