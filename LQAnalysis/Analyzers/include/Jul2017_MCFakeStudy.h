#ifndef Jul2017_MCFakeStudy_h
#define Jul2017_MCFakeStudy_h

#include "AnalyzerCore.h"


class Jul2017_MCFakeStudy : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  Jul2017_MCFakeStudy();
  ~Jul2017_MCFakeStudy();

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
  bool  IsNearBJet(snu::KElectron Ele, std::vector<snu::KJet> bjetNoVetoColl);
  int   StepPassed(std::vector<snu::KMuon> MuColl, std::vector<snu::KElectron> EleColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, TString Option="");
  float ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut);
  float ConeCorrectedPT(snu::KMuon Mu, float TightIsoCut);
  int   GetFakeLepSrcType(snu::KElectron Ele, std::vector<snu::KJet> JetColl);
  float FakeRateMC(snu::KElectron Ele, TString Option);
  float FakeRateMC(snu::KMuon Mu, TString Option);
  void  ScanFakeRate(std::vector<snu::KElectron> FakePreColl, std::vector<snu::KJet> JetColl, float MVACut, float IsoCut, int NPtEdges, float PtEdges[], TString PreID, TString TightID, TString Label, TString Option);
  void  Draw1DFakePlot(std::vector<snu::KElectron> FakeColl, std::vector<snu::KJet> JetColl, float MVAB1Cut, float MVAB2Cut, float MVAECut, float IsoCut, int NPtEdges, float PtEdges[], TString Option="");
  void EmulateFRMeasurement(std::vector<snu::KElectron> ElePreColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, int NPtEdges, float PtEdges[], TString LooseID, TString TightID, TString Label, TString Option);
  void Draw1DClosurePlot(std::vector<snu::KMuon> MuPreColl, std::vector<snu::KElectron> ElePreColl, std::vector<snu::KJet> jetColl, std::vector<snu::KJet> bjetColl, std::vector<snu::KTruth> truthColl, float met, float MVAB1Cut, float MVAB2Cut, float MVAECut, float IsoCut, TString Option);
  void CheckMCClosure(std::vector<snu::KMuon> MuPreColl, std::vector<snu::KElectron> ElePreColl, std::vector<snu::KJet> JetNoVetoColl, std::vector<snu::KJet> BJetNoVetoColl, std::vector<snu::KTruth> TruthColl, TString EleLID, TString EleTID, TString MuLID, TString MuTID, TString Label, TString Option);

  void Draw2DFakePlot(std::vector<snu::KElectron> FakeColl, std::vector<snu::KJet> JetColl, int NMVACuts, float MVACuts[], int NIsoCuts, float IsoCuts[], int NPtEdges, float PtEdges[], TString Option="ConePt");

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( Jul2017_MCFakeStudy, 1);
};
#endif
