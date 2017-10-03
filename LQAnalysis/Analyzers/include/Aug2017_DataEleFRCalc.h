#ifndef Aug2017_DataEleFRCalc_h
#define Aug2017_DataEleFRCalc_h

#include "AnalyzerCore.h"


class Aug2017_DataEleFRCalc : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  Aug2017_DataEleFRCalc();
  ~Aug2017_DataEleFRCalc();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);
  float FakeRateData(snu::KMuon Mu, TString Option);
  float FakeRateData(snu::KElectron Ele, TString ID);
  float GetPreTrigPURW(int Nvtx);
  void DoSystRun(TString Cycle, TString Mode, std::vector<snu::KElectron> EleColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight=1.);
  bool IsNearBJet(snu::KElectron Ele, std::vector<snu::KJet> bjetNoVetoColl);
  float ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut);
  float ConeCorrectedPT(snu::KMuon Mu, float TightIsoCut);

  float GetFakeWeight(std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleLColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, TString Option);

  void CheckNormCR(std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, float MET, float METx, float METy, float weight, TString Label, TString Option);
  void OptimiseMETMTWCuts(std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, float MET, float METx, float METy, float weight, TString Option);
  void OptimiseIsoIPWP(std::vector<snu::KElectron> EleTColl, float weight, TString Option);
  void MeasureFakeRate(std::vector<snu::KElectron> ElePreColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, int NPtEdges, float PtEdges[], TString LooseID, TString TightID, TString Label, TString Option);

  void ScanFakeRate(std::vector<snu::KElectron> EleColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, int NMVACuts, float MVACuts[], int NIsoCuts, float IsoCuts[], int NPtEdges, float PtEdges[], TString Option);
  void MeasureHighdXYFakeRate(std::vector<snu::KElectron> EleLowd0TColl, std::vector<snu::KElectron> EleLowd0LColl, std::vector<snu::KElectron> EleHighd0TColl, std::vector<snu::KElectron> EleHighd0LColl, std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, std::vector<snu::KTruth> TruthColl, float MET, float METx, float METy, float weight, TString Label, TString Option);
  void CheckTrilepCRs(std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label="", TString Option="");
  void ValidateID(std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, float weight);

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( Aug2017_DataEleFRCalc, 1);
};
#endif
