// $Id: Feb2018_SihyunsProblem.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQFeb2018_SihyunsProblem Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Feb2018_SihyunsProblem.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Feb2018_SihyunsProblem);

 Feb2018_SihyunsProblem::Feb2018_SihyunsProblem() : AnalyzerCore(), out_muons(0) {

   SetLogName("Feb2018_SihyunsProblem");
   Message("In Feb2018_SihyunsProblem constructor", INFO);
   InitialiseAnalysis();
 }


 void Feb2018_SihyunsProblem::InitialiseAnalysis() throw( LQError ) {
   
   /// Initialise histograms
   MakeHistograms();  
   // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
   // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
   Message("Feb2016, HwA analysis", INFO);
   return;
 }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Loop///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Feb2018_SihyunsProblem::ExecuteEvents()throw( LQError ){


   bool SihyunsProblem=true;
   if(SihyunsProblem){
     if(!PassMETFilter()) return;
     if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return;
     if(!isData) weight*=MCweight;
     if(!k_isdata) { weight*=eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central);}

     if(!PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) return;
     if(!isData) weight*=WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);

       eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
     std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
  
       std::vector<snu::KElectron> electronLooseColl, electronTightColl;
         for(int i=0; i<(int)electronPreColl.size(); i++){
          if(PassIDCriteria(electronPreColl.at(i), "LMVA06Isop4IPp025p1sig4"))   electronLooseColl.push_back(electronPreColl.at(i));
          if(PassIDCriteria(electronPreColl.at(i), "POGWP90Isop06IPp025p1sig4")) electronTightColl.push_back(electronPreColl.at(i));
         }
     std::vector<snu::KElectron> electronColl;  if(!k_running_nonprompt){ electronColl=electronTightColl;} else{ electronColl=electronLooseColl;}
     std::vector<snu::KTruth>    truthColl;    eventbase->GetTruthSel()->Selection(truthColl);

     if( electronLooseColl.size()!=2 ) return;
     if( electronTightColl.size()!=2 ) return;
     if( !(electronColl.at(0).Pt()>25 && fabs(electronColl.at(0).Eta())<2.5 && electronColl.at(1).Pt()>15 && fabs(electronColl.at(1).Eta())<2.5) ) return;
     if( electronColl.at(0).Charge() != electronColl.at(1).Charge() ) return;
     if( !(fabs((electronColl.at(0)+electronColl.at(1)).M()-91.2)<30) ) return;

     int El1Type = GetLeptonType(electronColl.at(0), truthColl);
     int El2Type = GetLeptonType(electronColl.at(1), truthColl);

     int El1Idx  = GenMatchedIdx(electronColl.at(0), truthColl);
     int El2Idx  = GenMatchedIdx(electronColl.at(1), truthColl);

     int PhoIdx1=-1, PhoIdx2=-1;
     if(El1Idx==-1) PhoIdx1 = GetNearPhotonIdx(electronColl.at(0), truthColl); 
     if(El2Idx==-1) PhoIdx2 = GetNearPhotonIdx(electronColl.at(1), truthColl); 

     float dPtRel_l1gen=El1Idx!=-1? fabs(electronColl.at(0).Pt()-truthColl.at(El1Idx).Pt())/electronColl.at(0).Pt():-1.;
     float dPtRel_l2gen=El2Idx!=-1? fabs(electronColl.at(1).Pt()-truthColl.at(El2Idx).Pt())/electronColl.at(1).Pt():-1.;
 

     FillHist("Mee_NoTypeCut", (electronColl.at(0)+electronColl.at(1)).M(), weight, 60., 120., 60);
     FillHist("LepType_NoTypeCut", El1Type, weight, -10., 10., 20);
     FillHist("LepType_NoTypeCut", El2Type, weight, -10., 10., 20);
     FillHist("PTe1_NoTypeCut", electronColl.at(0).Pt(), weight, 0., 100., 20);
     FillHist("PTe2_NoTypeCut", electronColl.at(1).Pt(), weight, 0., 100., 20);
     FillHist("Etae1_NoTypeCut", electronColl.at(0).Eta(), weight, -2.5, 2.5, 10);
     FillHist("Etae2_NoTypeCut", electronColl.at(1).Eta(), weight, -2.5, 2.5, 10);
     FillHist("dPtRel_l1gen_NoTypeCut", dPtRel_l1gen, weight, -1., 2., 30);
     FillHist("dPtRel_l2gen_NoTypeCut", dPtRel_l2gen, weight, -1., 2., 30);
     FillHist("dPtRel_lgen_NoTypeCut", dPtRel_l1gen, weight, -1., 2., 30);
     FillHist("dPtRel_lgen_NoTypeCut", dPtRel_l2gen, weight, -1., 2., 30);

     
     bool SihyunTypeReq=false;
     if(El1Type==4 || El1Type==5 || El1Type==-5 || El1Type==-6) SihyunTypeReq=true;
     if(El2Type==4 || El2Type==5 || El2Type==-5 || El2Type==-6) SihyunTypeReq=true;
     if(!SihyunTypeReq) return;

     FillHist("Mee_TypeCut", (electronColl.at(0)+electronColl.at(1)).M(), weight, 60., 120., 60);
     FillHist("LepType_TypeCut", El1Type, weight, -10., 10., 20);
     FillHist("LepType_TypeCut", El2Type, weight, -10., 10., 20);
     FillHist("PTe1_TypeCut", electronColl.at(0).Pt(), weight, 0., 100., 20);
     FillHist("PTe2_TypeCut", electronColl.at(1).Pt(), weight, 0., 100., 20);
     FillHist("Etae1_TypeCut", electronColl.at(0).Eta(), weight, -2.5, 2.5, 10);
     FillHist("Etae2_TypeCut", electronColl.at(1).Eta(), weight, -2.5, 2.5, 10);
     FillHist("dPtRel_l1gen_TypeCut", dPtRel_l1gen, weight, -1., 2., 30);
     FillHist("dPtRel_l2gen_TypeCut", dPtRel_l2gen, weight, -1., 2., 30);
     FillHist("dPtRel_lgen_TypeCut", dPtRel_l1gen, weight, -1., 2., 30);
     FillHist("dPtRel_lgen_TypeCut", dPtRel_l2gen, weight, -1., 2., 30);
      
     
     bool PrintHistory=false;
     if(PrintHistory){
       PrintTruth();
       cout<<"REle1Pt: "<<electronColl.at(0).Pt()<<" REle1Eta: "<<electronColl.at(0).Eta()<<" REle1Phi: "<<electronColl.at(0).Phi()<<" REle2Pt: "<<electronColl.at(1).Pt()<<" REle2Eta: "<<electronColl.at(1).Eta()<<" REle2Phi: "<<electronColl.at(1).Phi()<<endl;
       cout<<"El1Idx: "<<El1Idx<<" El1Type: "<<El1Type<<" El2Idx: "<<El2Idx<<" El2Type: "<<El2Type<<endl;
       cout<<"PhoIdx1: "<<PhoIdx1<<" PhoIdx2: "<<PhoIdx2<<endl;
       cout<<"==============================================================="<<endl;
     }

   }

   bool ClosureTest=false;
   if(ClosureTest){
     if(!PassMETFilter()) return;
     if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return;
     if(!isData) weight*=MCweight;
     if(!k_isdata) { weight*=eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central);}


       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
       eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
       eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.1);//Not in ID, but additional safe WP
       eventbase->GetElectronSel()->SetdxySigMax(3.);
     std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->Selection(electronColl);
     std::vector<snu::KTruth>    truthColl;    eventbase->GetTruthSel()->Selection(truthColl);

     std::vector<snu::KElectron> PromptColl=SkimLepColl(electronColl, truthColl, "Prompt");
     std::vector<snu::KElectron> FakeColl=SkimLepColl(electronColl, truthColl, "HFakeNHConv");

     FillHist("NPr", PromptColl.size(), weight, 0., 10., 10);
     FillHist("NFake", FakeColl.size(), weight, 0., 10., 10);
     if(PromptColl.size()!=2) return;
     if(FakeColl.size()==0) return;
     //if( !(PromptColl.at(0).Pt()>25 || FakeColl.at(0).Pt()>25) ) return;

     for(int i=0; i<(int) FakeColl.size(); i++){
       int LepType=GetLeptonType(FakeColl.at(i), truthColl);

       if(LepType==-1) FillHist("M3e_TypeM1", (PromptColl.at(0)+PromptColl.at(1)+FakeColl.at(i)).M(), weight, 0., 200., 200);
       if(LepType==-2) FillHist("M3e_TypeM2", (PromptColl.at(0)+PromptColl.at(1)+FakeColl.at(i)).M(), weight, 0., 200., 200);
       if(LepType==-3) FillHist("M3e_TypeM3", (PromptColl.at(0)+PromptColl.at(1)+FakeColl.at(i)).M(), weight, 0., 200., 200);
       if(LepType==-4) FillHist("M3e_TypeM4", (PromptColl.at(0)+PromptColl.at(1)+FakeColl.at(i)).M(), weight, 0., 200., 200);
       if(LepType==-5) FillHist("M3e_TypeM5", (PromptColl.at(0)+PromptColl.at(1)+FakeColl.at(i)).M(), weight, 0., 200., 200);
       if(LepType==-6) FillHist("M3e_TypeM6", (PromptColl.at(0)+PromptColl.at(1)+FakeColl.at(i)).M(), weight, 0., 200., 200);
       if(LepType==-6) FillHist("M3e_TypeM6", (PromptColl.at(0)+PromptColl.at(1)+FakeColl.at(i)).M(), weight, 0., 200., 200);
       if(LepType==4) FillHist("M3e_Type4", (PromptColl.at(0)+PromptColl.at(1)+FakeColl.at(i)).M(), weight, 0., 200., 200);
       if(LepType==5) FillHist("M3e_Type5", (PromptColl.at(0)+PromptColl.at(1)+FakeColl.at(i)).M(), weight, 0., 200., 200);

     }
   }


/////////////////////////////////////////////////////////////////////////////////// 


return;
}// End of execute event loop
  


void Feb2018_SihyunsProblem::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}

/*
bool Feb2018_SihyunsProblem::IsFinalPhotonSt23(std::vector<snu::KTruth> TruthColl){
//In Some XG proc events, it seems there is status 23 photon, yet no status 1 photon and no other genparticle is daughter of this photon.
//This is to check whether this is the case for the event.
//And this is designed only for 1 hard photon case as W+G or Z+G or TT+G

  bool IsFinalGammaStatus23 = false;
  bool HasStatus23Photon    = false;
  for(int i=2; i<TruthColl.size(); i++){
    int fpid  = fabs(TruthColl.at(i).PdgId());
    int GenSt = TruthColl.at(i).GenStatus();
    int MPID_direct= TruthColl.at(TruthColl.at(i).IndexMother()).PdgId();
    if( !((fpid!=22 && MPID_direct==22) || (fpid==22 && (GenSt==23||GenSt==1))) ) continue;

    int LastSelfIdx  = LastSelfMotherIdx(i,TruthColl);
    int LastSelfSt   = TruthColl.at(LastSelfIdx).GenStatus();
    int MotherIdx    = FirstNonSelfMotherIdx(i,TruthColl);
    int LastSelfMIdx=-1, MStatus_orig=-1;
    if(MotherIdx!=-1){
      LastSelfMIdx = LastSelfMotherIdx(MotherIdx,TruthColl);
      MStatus_orig = TruthColl.at(LastSelfMIdx).GenStatus();
    }

    if(fpid==22){
      if(GenSt==23) {HasStatus23Photon=true; IsFinalGammaStatus23=true;}
      else if(GenSt==1 && LastSelfSt==23) {IsFinalGammaStatus23=false; break;}//a)
    }
    else if( MPID_direct==22 && MStatus_orig==23 ){ IsFinalGammaStatus23=false; break;}//b)
  }

  if(!HasStatus23Photon) return false;
  
  return IsFinalGammaStatus23;

//footnotes
//a) The status 1 photon is end of the history of status 23 photon.
//b) Some particle is daughter of status 23 photon.
}

int Feb2018_SihyunsProblem::GetPhotonType(int PhotonIdx, std::vector<snu::KTruth> TruthColl){
//Type : 
// 0: Invalid input or Error or HardScatter is input when hardscatter is not final state
// 1: HardScatter / 2: Else prompt daughter(l,V,t)
//-1: Reserved for unmatched(Not used now) / -2: Hadronic origin

  if( PhotonIdx<2 ) return 0;
  if( !(TruthColl.at(PhotonIdx).PdgId()==22 && (TruthColl.at(PhotonIdx).GenStatus()==1 || TruthColl.at(PhotonIdx).GenStatus()==23)) ) return 0;

  if(TruthColl.at(PhotonIdx).GenStatus()==23){
    if(IsFinalPhotonSt23(TruthColl)) return 1;
    else                             return 0;
  }//From this pt, only St1 Photon is treated.

  int PhotonType=0;
  int LastSelfIdx    = LastSelfMotherIdx(PhotonIdx,TruthColl);
  int MotherIdx      = FirstNonSelfMotherIdx(PhotonIdx,TruthColl);
  int fMPID=0, Status_orig=0;
  bool HadronicOrigin = false;
    if( PhotonIdx!=-1 ){ Status_orig    = TruthColl.at(LastSelfIdx).GenStatus();
                         HadronicOrigin = HasHadronicAncestor(PhotonIdx, TruthColl);
                       }                           
    if( MotherIdx!=-1 ){ fMPID          = fabs(TruthColl.at(MotherIdx).PdgId()); }


  if     (  Status_orig>20 && Status_orig<30   ) PhotonType= 1;//1)
  else if(       fMPID==23 || fMPID==25        ) PhotonType= 1;//2)
  else if( fMPID==24 || fMPID==6  || fMPID==37 ) PhotonType= 2;//3)
  else if(           HadronicOrigin            ) PhotonType=-2;//4)
  else if( fMPID==11 || fMPID==13 || fMPID==15 ) PhotonType= 2;//5)
  else                                           PhotonType= 0;
  
  return PhotonType;
//footnote
//1) In case of hard scattered photon, they may have history; GenSt=23>...>GenSt1, And depending on generator, their mother can be explicitly
//  written in history as Z>GG St2 but sometimes their field mother is not designated for avoiding confusion of gauge symmetry.
//  e.g. qq>llG instead of qq>Z>llG
//  To cover all the case, first thing to check is original state of photon is hard scattered or not
//2) Sometimes, if there is no correction on gamma is applied on PS step, photon's final state is 1 before any history.
//   e.g. G;St=1, Mother=Z ; cannot find hard scatter Z / in such case only possibility is to check mother.
//   But in some case, it is not obvious. because in PS step, charged ptls can radiate photons including bosons.
//   So it is artificial to distinguish hard scattered photon and photon from PS step. And fraction of this case is not negilgible; very frequently observed.
//3) top and charged bosons radiate photons, and some case the photon is very energetic.
//4) This category does not include tops. Photons from hadrons and quarks. But predominantly, in most of the cases they are daughter of pi0.
//   But rarely other mesons as eta, B, or even some quarks can also radiate energetic photons.
//5) Photons radiated from lepton FSR, but sometimes they radiate quite energetic photons.
}

bool Feb2018_SihyunsProblem::IsHardPhotonConverted(std::vector<snu::KTruth> TruthColl){

  bool IsConverted=false;

  for(int i=2; i<TruthColl.size(); i++){
    int fpid  = fabs(TruthColl.at(i).PdgId());
    int MPID_direct= TruthColl.at(TruthColl.at(i).IndexMother()).PdgId();
    if( !(fpid!=22 && MPID_direct==22) ) continue;

    int MotherIdx  = FirstNonSelfMotherIdx(i,TruthColl);
    int LastSelfMIdx=-1, MStatus_orig=-1;
    if(MotherIdx!=-1){
      LastSelfMIdx = LastSelfMotherIdx(MotherIdx,TruthColl);
      MStatus_orig = TruthColl.at(LastSelfMIdx).GenStatus();
    }

    if( MPID_direct==22 && MStatus_orig==23 ){ IsConverted=true; break; }
  }

  return IsConverted;

}
*/

void Feb2018_SihyunsProblem::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
//  string analysisdir = getenv("FILEDIR");  
//  if(!k_isdata) reweightPU = new Reweight((analysisdir + "SNUCAT_Pileup.root").c_str());

  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
//  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
//  DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

Feb2018_SihyunsProblem::~Feb2018_SihyunsProblem() {
  
  Message("In Feb2018_SihyunsProblem Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Feb2018_SihyunsProblem::FillCutFlow(TString cut, float weight){
  
  if(GetHist("cutflow_W") && GetHist("cutflow_N")){
    GetHist("cutflow_W")->Fill(cut,weight);
    GetHist("cutflow_N")->Fill(cut,1);
  }
  else{
    if(!GetHist("cutflow_W")){
      AnalyzerCore::MakeHistograms("cutflow_W", 10, 0., 10.);
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(1,"NoCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(2,"TriggerCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(3,"EventCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(4,"VertexCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(5,"3lOSCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(6,"M(#mu#mu)>12");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(7,"#geq3j");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(8,"#geq1b");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(9,"ZVeto");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(10,"M(A)Range");

    }
    if(!GetHist("cutflow_N")){
      AnalyzerCore::MakeHistograms("cutflow_N", 10, 0., 10.);
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(1,"NoCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(2,"TriggerCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(3,"EventCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(4,"VertexCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(5,"3lOSCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(6,"M(#mu#mu)>12");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(7,"#geq3j");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(8,"#geq1b");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(9,"ZVeto");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(10,"M(A)Range");
    }
  }
}



void Feb2018_SihyunsProblem::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Feb2018_SihyunsProblem::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //Basic Hists/////////////////////////////////////////////////
  //Data properties before selection  
  AnalyzerCore::MakeHistograms("Basic_Nj_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Nb_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NmuL_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NeL_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NmuT_orig",10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NeT_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_METdist_orig", 100, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_NmuT_weT", 10, 0., 10.);

  AnalyzerCore::MakeHistograms("Basic_Ptmu1_3mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu2_3mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu1_1e2mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Pte_1e2mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j1_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j2_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j3_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j4_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_b1_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_b2_Et_orig", 200, 0., 200.);

  //After Nlcut
  AnalyzerCore::MakeHistograms("Basic_Nvtx_NoRW_wNlOScut", 50, 0., 50.);
  AnalyzerCore::MakeHistograms("Basic_Nvtx_PURW_wNlOScut", 50, 0., 50.);
  AnalyzerCore::MakeHistograms("Basic_MET_wNlOScut", 200, 0., 200.);


  AnalyzerCore::MakeHistograms("Basic_Nj_wNlOScut", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Nb_wNlOScut", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Ne_wNlOScut", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Nmu_wNlOScut", 10, 0., 10.);

  AnalyzerCore::MakeHistograms("Basic_Pte_wNlOScut", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu1_wNlOScut", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu2_wNlOScut", 200, 0., 200.);

  AnalyzerCore::MakeHistograms("Basic_Mmumu_wNlOScut", 200, 0., 200.);


  //After Nljcut
  AnalyzerCore::MakeHistograms("Basic_Nb_wNljcut", 10, 0., 10.);


  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Feb2018_SihyunsProblemCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Feb2018_SihyunsProblem::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
