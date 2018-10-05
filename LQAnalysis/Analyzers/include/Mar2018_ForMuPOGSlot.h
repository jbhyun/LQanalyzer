#ifndef Mar2018_ForMuPOGSlot_h
#define Mar2018_ForMuPOGSlot_h

#include "AnalyzerCore.h"


class Mar2018_ForMuPOGSlot : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  Mar2018_ForMuPOGSlot();
  ~Mar2018_ForMuPOGSlot();

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
  float ConeCorrectedPT(snu::KMuon     Mu,  double TightIsoCut);
  int   GetFakeLepJetSrcType(snu::KMuon Mu, std::vector<snu::KJet> JetColl);
  float FakeRateMC(snu::KElectron Ele, TString Option);
  float FakeRateMC(snu::KMuon Mu, TString Option);
  void  Draw1DFakePlot(std::vector<snu::KElectron> FakeColl, std::vector<snu::KJet> JetColl, float MVACut, float IsoCut, int NPtEdges, float PtEdges[], TString Option="");
  void  Draw1DFakePlot(std::vector<snu::KElectron> FakeColl, std::vector<snu::KJet> JetColl, float MVAB1Cut, float MVAB2Cut, float MVAECut, float IsoCut, int NPtEdges, float PtEdges[], TString Option="");
  void Draw1DClosurePlot(std::vector<snu::KMuon> MuPreColl, std::vector<snu::KElectron> ElePreColl, std::vector<snu::KJet> jetColl, std::vector<snu::KJet> bjetColl, std::vector<snu::KTruth> truthColl, float met, float MVAB1Cut, float MVAB2Cut, float MVAECut, float IsoCut, TString Option);


  void CheckFakeSources(std::vector<snu::KMuon> muonColl, std::vector<snu::KMuon> muonLooseColl, std::vector<snu::KElectron> electronLooseColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, std::vector<snu::KTruth> truthColl, float weight, TString Label, TString Option);
  void CheckIDEfficiency(std::vector<snu::KMuon> muonColl, std::vector<snu::KTruth> truthColl, float weight, TString Label, TString Option);
  void CheckIDVarSensitivity(std::vector<snu::KMuon> muonColl, std::vector<snu::KJet> JetColl, std::vector<snu::KTruth> truthColl, float weight, TString Label, TString Option);
  void InspectFakeRate(std::vector<snu::KMuon> muonLooseColl, std::vector<snu::KJet> JetColl, std::vector<snu::KTruth> truthColl, TString LooseID, TString TightID, float weight, TString Label, TString Option);
  void ScanFakeRate(std::vector<snu::KMuon> muonLooseColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KTruth> truthColl, TString LooseID, TString TightID, int Nd0Cuts, float d0Cuts[], int NChi2Cuts, float Chi2Cuts[], float weight, TString Label, TString Option);
  void CheckTriggerBias(std::vector<snu::KMuon> muonLooseColl, std::vector<snu::KJet> JetColl, std::vector<snu::KTruth> truthColl, TString LooseID, TString TightID, float weight, TString Label, TString Option);
  void EmulateFRMeasurement(std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, float MET, float METx, float METy, std::vector<snu::KTruth> truthColl, TString LooseID, TString TightID, float weight, TString Label, TString Option);
  void CheckMCClosure(std::vector<snu::KMuon> MuPreColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetNoVetoColl, std::vector<snu::KTruth> TruthColl, TString LooseID, TString TightID, float weight, TString Label, TString Option);


  void DrawPlotsForPOGSlot(std::vector<snu::KMuon> muonColl, std::vector<snu::KMuon> muonLooseColl, std::vector<snu::KElectron> electronLooseColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, std::vector<snu::KTruth> truthColl, float weight, TString Label, TString Option);

  int GetFakeLepSrcIdx(snu::KMuon Mu, std::vector<snu::KTruth> TruthColl);

  void PerformanceComp(std::vector<snu::KMuon> muonColl, std::vector<snu::KMuon> muonLooseColl, std::vector<snu::KElectron> electronTightColl, std::vector<snu::KElectron> electronLooseColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, std::vector<snu::KTruth> truthColl, float weight, TString Label, TString Option);
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( Mar2018_ForMuPOGSlot, 1);
};
#endif
