#ifndef Nov2017_GenPDFQ2Syst_h
#define Nov2017_GenPDFQ2Syst_h

#include "AnalyzerCore.h"


class Nov2017_GenPDFQ2Syst : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  Nov2017_GenPDFQ2Syst();
  ~Nov2017_GenPDFQ2Syst();

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
  void ScanFakeRate(std::vector<snu::KElectron> EleColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, int NMVACuts, float MVACuts[], int NIsoCuts, float IsoCuts[], int NPtEdges, float PtEdges[], TString Option);
  float ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut);
  float ConeCorrectedPT(snu::KMuon     Mu , float TightIsoCut);

  float GetFakeWeight(std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleLColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, TString Option);


//  void CheckGenPDFQ2Syst(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option);

  void CheckGenPDFQ2Syst(std::vector<snu::KMuon>& MuTColl, std::vector<snu::KMuon>& MuLColl, std::vector<snu::KElectron>& EleTColl, std::vector<snu::KElectron>& EleLColl, std::vector<snu::KJet>& JetColl, std::vector<snu::KJet>& BJetColl, float MET, float METx, float METy, float weight, TString Label, TString Option);

  void CheckTopPtSpectrum(std::vector<snu::KTruth>& TruthColl, float weight);
  float GetTopLOtoNLOReweight(std::vector<snu::KTruth>& TruthColl);
  void CheckLOvsNLOTopPtBias(std::vector<snu::KMuon>& MuTColl, std::vector<snu::KMuon>& MuLColl, std::vector<snu::KElectron>& EleTColl, std::vector<snu::KElectron>& EleLColl, std::vector<snu::KJet>& JetColl, std::vector<snu::KJet>& BJetColl, float MET, float METx, float METy, std::vector<snu::KTruth>& TruthColl, float weight, TString Label, TString Option);
  void CheckGenFilterBias(std::vector<snu::KMuon>& MuTColl, std::vector<snu::KMuon>& MuLColl, std::vector<snu::KElectron>& EleTColl, std::vector<snu::KElectron>& EleLColl, std::vector<snu::KJet>& JetColl, std::vector<snu::KJet>& BJetColl, float MET, float METx, float METy, std::vector<snu::KTruth>& TruthColl, float weight, TString Label, TString Option);

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( Nov2017_GenPDFQ2Syst, 1);
};
#endif
