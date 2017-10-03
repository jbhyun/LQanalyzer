#ifndef Aug2017_MuIDChoice_h
#define Aug2017_MuIDChoice_h

#include "AnalyzerCore.h"


class Aug2017_MuIDChoice : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  Aug2017_MuIDChoice();
  ~Aug2017_MuIDChoice();

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
  int   GetFakeLepJetSrcType(snu::KMuon Mu, std::vector<snu::KJet> JetColl);
  float FakeRateMC(snu::KElectron Ele, TString Option);
  void  Draw1DFakePlot(std::vector<snu::KElectron> FakeColl, std::vector<snu::KJet> JetColl, float MVACut, float IsoCut, int NPtEdges, float PtEdges[], TString Option="");
  void  Draw1DFakePlot(std::vector<snu::KElectron> FakeColl, std::vector<snu::KJet> JetColl, float MVAB1Cut, float MVAB2Cut, float MVAECut, float IsoCut, int NPtEdges, float PtEdges[], TString Option="");
  void Draw1DClosurePlot(std::vector<snu::KMuon> MuPreColl, std::vector<snu::KElectron> ElePreColl, std::vector<snu::KJet> jetColl, std::vector<snu::KJet> bjetColl, std::vector<snu::KTruth> truthColl, float met, float MVAB1Cut, float MVAB2Cut, float MVAECut, float IsoCut, TString Option);



  void ComparePunziDisc(std::vector<snu::KMuon> muonColl, std::vector<snu::KMuon> muonLooseColl, std::vector<snu::KElectron> electronColl, std::vector<snu::KElectron> electronLooseColl, std::vector<snu::KJet> jetColl, std::vector<snu::KJet> bjetColl, float weight, TString Label, TString Option);

  int GetFakeLepSrcIdx(snu::KMuon Mu, std::vector<snu::KTruth> TruthColl);
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( Aug2017_MuIDChoice, 1);
};
#endif
