// $Id: Jun2017_ConversionStudy.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQJun2017_ConversionStudy Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Jun2017_ConversionStudy.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Jun2017_ConversionStudy);

 Jun2017_ConversionStudy::Jun2017_ConversionStudy() : AnalyzerCore(), out_muons(0) {

   SetLogName("Jun2017_ConversionStudy");
   Message("In Jun2017_ConversionStudy constructor", INFO);
   InitialiseAnalysis();
 }


 void Jun2017_ConversionStudy::InitialiseAnalysis() throw( LQError ) {
   
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

void Jun2017_ConversionStudy::ExecuteEvents()throw( LQError ){


   bool ExternalConversionGenTest=false;
   if(ExternalConversionGenTest){
     if(!PassMETFilter()) return;
     if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return;
     if(!isData) weight*=MCweight;
     if(!k_isdata) { weight*=eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central);}


       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
       eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
       eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.1);//Not in ID, but additional safe WP
       eventbase->GetElectronSel()->SetdxySigMax(3.);
     std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->Selection(electronColl);
       eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
       eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
       eventbase->GetMuonSel()->SetBSdxy(0.05);                eventbase->GetMuonSel()->SetdxySigMax(3.);
       eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.1);
     std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->Selection(muonColl,true);
     std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
    
     //electron External Conversion
     for(int i=0; i<electronColl.size(); i++){

       int LeptonType=GetLeptonType(electronColl.at(i), truthColl);

       //External Conversion Properties
       if(LeptonType==-1){
         int NMatched=0; float dRmin=999.; int MatchedIdx=-1;
         for(int j=2; j<truthColl.size(); j++){
           if( truthColl.at(j).IndexMother()<0 ) continue;
           if( !(fabs(truthColl.at(j).PdgId())==22 && (truthColl.at(j).GenStatus()==1 || truthColl.at(j).GenStatus()==23)) ) continue;
           if( GetPhotonType(j, truthColl)!=1  )  continue;
           //if( truthColl.at(i).GenStatus()==23 && !IsFinalPhotonSt23(truthColl) ) continue;
    
           if( electronColl.at(i).DeltaR(truthColl.at(j))>0.4 ) continue;
           NMatched++;
           if( electronColl.at(i).DeltaR(truthColl.at(j))<dRmin ){ dRmin=electronColl.at(i).DeltaR(truthColl.at(j)); MatchedIdx=j; }
         }//End of TruthColl
         FillHist("NMatched_ext", NMatched, weight, 0., 10., 10);
  
         if(MatchedIdx!=-1){
           FillHist("dRmin_eg_ext", electronColl.at(i).DeltaR(truthColl.at(MatchedIdx)), weight, 0., 1., 100);
           if(truthColl.at(MatchedIdx).Pt()>0){
             FillHist("RelPT_eog_ext", electronColl.at(i).Pt()/truthColl.at(MatchedIdx).Pt(), weight, 0., 2., 200);
           }
         }
       }//End of Unmatched
     }//End of electron Coll

     //muon External Conversion
     for(int i=0; i<muonColl.size(); i++){

       int LeptonType=GetLeptonType(muonColl.at(i), truthColl);

       //External Conversion Properties
       if(LeptonType==-1){
         int NMatched=0; float dRmin=999.; int MatchedIdx=-1;
         for(int j=2; j<truthColl.size(); j++){
           if( truthColl.at(j).IndexMother()<0 ) continue;
           if( !(fabs(truthColl.at(j).PdgId())==22 && (truthColl.at(j).GenStatus()==1 || truthColl.at(j).GenStatus()==23)) ) continue;
           if( GetPhotonType(j, truthColl)!=1  )  continue;
           //if( truthColl.at(i).GenStatus()==23 && !IsFinalPhotonSt23(truthColl) ) continue;
    
           if( muonColl.at(i).DeltaR(truthColl.at(j))>0.4 ) continue;
           NMatched++;
           if( muonColl.at(i).DeltaR(truthColl.at(j))<dRmin ){ dRmin=muonColl.at(i).DeltaR(truthColl.at(j)); MatchedIdx=j; }
         }//End of TruthColl
         FillHist("NMatched_ext", NMatched, weight, 0., 10., 10);
  
         if(MatchedIdx!=-1){
           FillHist("dRmin_mug_ext", muonColl.at(i).DeltaR(truthColl.at(MatchedIdx)), weight, 0., 1., 100);
           if(truthColl.at(MatchedIdx).Pt()>0){
             FillHist("RelPT_muog_ext", muonColl.at(i).Pt()/truthColl.at(MatchedIdx).Pt(), weight, 0., 2., 200);
           }
         }
       }//End of Unmatched
     }//End of muon Coll


     //Internal Conversion
     int Mu1Idx=-1, Mu2Idx=-1, TmpMuIdx1=-1, TmpMuIdx2=-1;
     int Ele1Idx=-1, Ele2Idx=-1, TmpEleIdx1=-1, TmpEleIdx2=-1;
     int NCount=0;
     for(int i=2; i<truthColl.size(); i++){
       if( truthColl.at(i).IndexMother()<0 ) continue;
       if( truthColl.at(i).GenStatus()!=1 )  continue; 
       int fpid=fabs(truthColl.at(i).PdgId());
       if( !(fpid==11 || fpid==13) ) continue;

       int LeptonType=GetLeptonType(i, truthColl);
       if(LeptonType==4){
         NCount++;
         if(NCount==1){
           if(fpid==11) TmpEleIdx1=i;
           else if(fpid==13) TmpMuIdx1=i;
         }
         else if(NCount==2){
           if(fpid==11) TmpEleIdx2=i;
           else if(fpid==13) TmpMuIdx2=i;
         }
       }
     }
     if(TmpEleIdx1!=-1 && TmpEleIdx2!=-1){
       Ele1Idx = truthColl.at(TmpEleIdx1).Pt()>truthColl.at(TmpEleIdx2).Pt()? TmpEleIdx1:TmpEleIdx2;
       Ele2Idx = truthColl.at(TmpEleIdx1).Pt()>truthColl.at(TmpEleIdx2).Pt()? TmpEleIdx2:TmpEleIdx1;
       int PhoIdx  = FirstNonSelfMotherIdx(Ele1Idx,truthColl);
       bool Ele1AccPass=false, Ele2AccPass=false, Ele1AccMinPass=false, Ele2AccMinPass=false;
       if(truthColl.at(Ele1Idx).Pt()>25 && fabs(truthColl.at(Ele1Idx).Eta())<2.5) Ele1AccPass=true;
       if(truthColl.at(Ele2Idx).Pt()>25 && fabs(truthColl.at(Ele2Idx).Eta())<2.5) Ele2AccPass=true;
       if(truthColl.at(Ele1Idx).Pt()>10 && fabs(truthColl.at(Ele1Idx).Eta())<2.5) Ele1AccMinPass=true;
       if(truthColl.at(Ele2Idx).Pt()>10 && fabs(truthColl.at(Ele2Idx).Eta())<2.5) Ele2AccMinPass=true;
       if(Ele1AccPass) FillHist("AccComp", 0., weight, 0., 10., 10);
       if(Ele1AccPass && Ele2AccPass)  FillHist("AccComp", 1., weight, 0., 10., 10);
       if(Ele1AccPass && !Ele2AccMinPass) FillHist("AccComp", 2., weight, 0., 10., 10);
       if( (Ele1AccPass || Ele2AccPass) && ( !Ele1AccMinPass || !Ele2AccMinPass) ) FillHist("AccComp", 3., weight, 0., 10., 10);

       FillHist("Pte1_int", truthColl.at(Ele1Idx).Pt(), weight, 0., 200., 200);
       FillHist("Pte2_int", truthColl.at(Ele2Idx).Pt(), weight, 0., 200., 200);
       FillHist("Mee_int", (truthColl.at(Ele1Idx)+truthColl.at(Ele2Idx)).M(), weight, 0., 10., 200);
       FillHist("dRe1g_int", truthColl.at(Ele1Idx).DeltaR(truthColl.at(PhoIdx)), weight, 0., 1., 100);
       if(truthColl.at(PhoIdx).Pt()>0){
         FillHist("RelPt_e1og_int", truthColl.at(Ele1Idx).Pt()/truthColl.at(PhoIdx).Pt(), weight, 0., 2., 200);
       }
       if(Ele1AccPass && !Ele2AccPass){
         FillHist("Pte1_int_Case1", truthColl.at(Ele1Idx).Pt(), weight, 0., 200., 200);
         FillHist("Pte2_int_Case1", truthColl.at(Ele2Idx).Pt(), weight, 0., 200., 200);
         FillHist("Mee_int_Case1", (truthColl.at(Ele1Idx)+truthColl.at(Ele2Idx)).M(), weight, 0., 10., 200);
         FillHist("dRe1g_int_Case1", truthColl.at(Ele1Idx).DeltaR(truthColl.at(PhoIdx)), weight, 0., 1., 100);
         if(truthColl.at(PhoIdx).Pt()>0){
           FillHist("RelPt_e1og_int_Case1", truthColl.at(Ele1Idx).Pt()/truthColl.at(PhoIdx).Pt(), weight, 0., 2., 200);
         }
       }
       if(Ele1AccPass && !Ele2AccMinPass){
         FillHist("Pte1_int_Case2", truthColl.at(Ele1Idx).Pt(), weight, 0., 200., 200);
         FillHist("Pte2_int_Case2", truthColl.at(Ele2Idx).Pt(), weight, 0., 200., 200);
         FillHist("Mee_int_Case2", (truthColl.at(Ele1Idx)+truthColl.at(Ele2Idx)).M(), weight, 0., 10., 200);
         FillHist("dRe1g_int_Case2", truthColl.at(Ele1Idx).DeltaR(truthColl.at(PhoIdx)), weight, 0., 1., 100);
         if(truthColl.at(PhoIdx).Pt()>0){
           FillHist("RelPt_e1og_int_Case2", truthColl.at(Ele1Idx).Pt()/truthColl.at(PhoIdx).Pt(), weight, 0., 2., 200);
         }
       }
       if( (Ele1AccPass || Ele2AccPass) && ( !Ele1AccMinPass || !Ele2AccMinPass) ){
         FillHist("Pte1_int_Case3", truthColl.at(Ele1Idx).Pt(), weight, 0., 200., 200);
         FillHist("Pte2_int_Case3", truthColl.at(Ele2Idx).Pt(), weight, 0., 200., 200);
         FillHist("Mee_int_Case3", (truthColl.at(Ele1Idx)+truthColl.at(Ele2Idx)).M(), weight, 0., 10., 200);
         FillHist("dRe1g_int_Case3", truthColl.at(Ele1Idx).DeltaR(truthColl.at(PhoIdx)), weight, 0., 1., 100);
         if(truthColl.at(PhoIdx).Pt()>0){
           FillHist("RelPt_e1og_int_Case3", truthColl.at(Ele1Idx).Pt()/truthColl.at(PhoIdx).Pt(), weight, 0., 2., 200);
         }

       }

     }
     else if(TmpMuIdx1!=-1 && TmpMuIdx2!=-1){
       Mu1Idx = truthColl.at(TmpMuIdx1).Pt()>truthColl.at(TmpMuIdx2).Pt()? TmpMuIdx1:TmpMuIdx2;
       Mu2Idx = truthColl.at(TmpMuIdx1).Pt()>truthColl.at(TmpMuIdx2).Pt()? TmpMuIdx2:TmpMuIdx1;

       int PhoIdx  = FirstNonSelfMotherIdx(Mu1Idx,truthColl);
       bool Mu1AccPass=false, Mu2AccPass=false, Mu1AccMinPass=false, Mu2AccMinPass=false;
       if(truthColl.at(Mu1Idx).Pt()>10 && fabs(truthColl.at(Mu1Idx).Eta())<2.5) Mu1AccPass=true;
       if(truthColl.at(Mu2Idx).Pt()>10 && fabs(truthColl.at(Mu2Idx).Eta())<2.5) Mu2AccPass=true;
       if(Mu1AccPass) FillHist("AccComp_mu", 0., weight, 0., 10., 10);
       if(Mu1AccPass && Mu2AccPass)  FillHist("AccComp_mu", 1., weight, 0., 10., 10);
       if(Mu1AccPass && !Mu2AccPass) FillHist("AccComp_mu", 2., weight, 0., 10., 10);
       if( (Mu1AccPass || Mu2AccPass) && ( !Mu1AccPass || !Mu2AccPass) ) FillHist("AccComp_mu", 3., weight, 0., 10., 10);

       FillHist("Ptmu1_int", truthColl.at(Mu1Idx).Pt(), weight, 0., 200., 200);
       FillHist("Ptmu2_int", truthColl.at(Mu2Idx).Pt(), weight, 0., 200., 200);
       FillHist("Mmumu_int", (truthColl.at(Mu1Idx)+truthColl.at(Mu2Idx)).M(), weight, 0., 10., 200);
       FillHist("dRmu1g_int", truthColl.at(Mu1Idx).DeltaR(truthColl.at(PhoIdx)), weight, 0., 1., 100);
       if(truthColl.at(PhoIdx).Pt()>0){
         FillHist("RelPt_mu1og_int", truthColl.at(Mu1Idx).Pt()/truthColl.at(PhoIdx).Pt(), weight, 0., 2., 200);
       }
       if(Mu1AccPass && !Mu2AccPass){
         FillHist("Ptmu1_int_Case1", truthColl.at(Mu1Idx).Pt(), weight, 0., 200., 200);
         FillHist("Ptmu2_int_Case1", truthColl.at(Mu2Idx).Pt(), weight, 0., 200., 200);
         FillHist("Mmumu_int_Case1", (truthColl.at(Mu1Idx)+truthColl.at(Mu2Idx)).M(), weight, 0., 10., 200);
         FillHist("dRmu1g_int_Case1", truthColl.at(Mu1Idx).DeltaR(truthColl.at(PhoIdx)), weight, 0., 1., 100);
         if(truthColl.at(PhoIdx).Pt()>0){
           FillHist("RelPt_mu1og_int_Case1", truthColl.at(Mu1Idx).Pt()/truthColl.at(PhoIdx).Pt(), weight, 0., 2., 200);
         }
       }
       if( (Mu1AccPass || Mu2AccPass) && ( !Mu1AccPass || !Mu2AccPass) ){
         FillHist("Ptmu1_int_Case3", truthColl.at(Mu1Idx).Pt(), weight, 0., 200., 200);
         FillHist("Ptmu2_int_Case3", truthColl.at(Mu2Idx).Pt(), weight, 0., 200., 200);
         FillHist("Mmumu_int_Case3", (truthColl.at(Mu1Idx)+truthColl.at(Mu2Idx)).M(), weight, 0., 10., 200);
         FillHist("dRmu1g_int_Case3", truthColl.at(Mu1Idx).DeltaR(truthColl.at(PhoIdx)), weight, 0., 1., 100);
         if(truthColl.at(PhoIdx).Pt()>0){
           FillHist("RelPt_mu1og_int_Case3", truthColl.at(Mu1Idx).Pt()/truthColl.at(PhoIdx).Pt(), weight, 0., 2., 200);
         }
       }

     }

     return;


//       for(int i=2; i<truthColl.size(); i++){
//         if( truthColl.at(i).IndexMother()<0 ) continue;
//         if( !(fabs(truthColl.at(i).PdgId())==22 && (truthColl.at(i).GenStatus()==1 || truthColl.at(i).GenStatus()==23)) ) continue;
//         if( truthColl.at(i).GenStatus()==23 && !IsFinalPhotonSt23(truthColl) ) continue;
//  
//         int PhotonType=GetPhotonType(i, truthColl);
//         FillHist("PhotonType_nocut", PhotonType, weight, -10., 10., 20);
//  
//         if(truthColl.at(i).Pt()>25){ FillHist("PhotonType_Ptgt25", PhotonType, weight, -10., 10., 20);
//           //if(PhotonType==6) {cout<<"Idx "<<i<<endl; PrintTruth();}
//         }

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

     for(int i=0; i<FakeColl.size(); i++){
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


   bool ConversionIPCheck=false;
   if(ConversionIPCheck){

       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
       eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
     //  eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);//Not in ID, but additional safe WP
     std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->Selection(electronColl);
       eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
       eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
       eventbase->GetMuonSel()->SetBSdxy(0.05);                eventbase->GetMuonSel()->SetdxySigMax(3.);
       eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.1);
     std::vector<snu::KMuon>   muonColl; eventbase->GetMuonSel()->Selection(muonColl,true);
     std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);

     std::vector<snu::KElectron> EleIntConvColl;
       for(int i=0; i<electronColl.size(); i++){
         int LepType=GetLeptonType(electronColl.at(i), truthColl);
         if(LepType>3) EleIntConvColl.push_back(electronColl.at(i));
       }
     std::vector<snu::KElectron> EleExtConvColl;
       for(int i=0; i<electronColl.size(); i++){
         int LepType=GetLeptonType(electronColl.at(i), truthColl);
         if(LepType>3) EleIntConvColl.push_back(electronColl.at(i));
       }

     for(int i=0; i<electronColl.size(); i++){

       int LepType=GetLeptonType(electronColl.at(i), truthColl);
       float d0=fabs(electronColl.at(i).dxy()), dz=fabs(electronColl.at(i).dz()), d0Sig=fabs(electronColl.at(i).dxySig());

       //d0
       if     (LepType==1)                 FillHist("d0_Type1", d0, weight, 0., 0.2, 200);
       else if(LepType==3)                 FillHist("d0_Type3", d0, weight, 0., 0.2, 200);
       else if(LepType>3)                  FillHist("d0_Type45_IntConv", d0, weight, 0., 0.2, 200);
       else if(LepType==-1)                FillHist("d0_Typem1", d0, weight, 0., 0.2, 200);
       else if(LepType<=-2 && LepType>=-4) FillHist("d0_Typem234", d0, weight, 0., 0.2, 200);
       else if(LepType<-4)                 FillHist("d0_Typem56_ExtConv", d0, weight, 0., 0.2, 200);

       //dz
       if     (LepType==1)                 FillHist("dz_Type1", dz, weight, 0., 0.2, 200);
       else if(LepType==3)                 FillHist("dz_Type3", dz, weight, 0., 0.2, 200);
       else if(LepType>3)                  FillHist("dz_Type45_IntConv", dz, weight, 0., 0.2, 200);
       else if(LepType==-1)                FillHist("dz_Typem1", dz, weight, 0., 0.2, 200);
       else if(LepType<=-2 && LepType>=-4) FillHist("dz_Typem234", dz, weight, 0., 0.2, 200);
       else if(LepType<-4)                 FillHist("dz_Typem56_ExtConv", dz, weight, 0., 0.2, 200);

       //d0Sig
       if     (LepType==1)                 FillHist("d0Sig_Type1", d0Sig, weight, 0., 10, 100);
       else if(LepType==3)                 FillHist("d0Sig_Type3", d0Sig, weight, 0., 10, 100);
       else if(LepType>3)                  FillHist("d0Sig_Type45_IntConv", d0Sig, weight, 0., 10, 100);
       else if(LepType==-1)                FillHist("d0Sig_Typem1", d0Sig, weight, 0., 10, 100);
       else if(LepType<=-2 && LepType>=-4) FillHist("d0Sig_Typem234", d0Sig, weight, 0., 10, 100);
       else if(LepType<-4)                 FillHist("d0Sig_Typem56_ExtConv", d0Sig, weight, 0., 10, 100);
     }//End of ele loop
 
   
   }//End of IPCheck


   bool ZGDataMCComp=false;
   if(ZGDataMCComp){

     if(!isData) weight*=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
     if(!isData) weight*=MCweight;
     float pileup_reweight=1.;
     if(!k_isdata) { pileup_reweight*=eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central);}


     int Pass_Trigger1=0, Pass_Trigger2=0;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ) Pass_Trigger1++;
     if( PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger2++;

     float trigger_period_weight=1.;
     bool Pass_Trigger=false;
     if(isData){
       int DataPeriod=GetDataPeriod();
       if( DataPeriod>0 && DataPeriod<7 && Pass_Trigger1==1 ) Pass_Trigger=true;
       else if( DataPeriod==7 && Pass_Trigger2==1 ) Pass_Trigger=true;
     }
     else{
       if( Pass_Trigger1>0 || Pass_Trigger2>0 ) Pass_Trigger=true;
       trigger_period_weight=(Pass_Trigger1*27.257618+Pass_Trigger2*8.605696)/35.863314;
     }
     weight*=trigger_period_weight;
     if(!Pass_Trigger) return;
     if(!PassMETFilter()) return;  FillCutFlow("EventCut", weight*pileup_reweight);
     if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; FillCutFlow("VertexCut", weight*pileup_reweight);

     /**PreSelCut***********************************************************************************************/
     //Intended for Code speed boosting up.
       eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
       eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, true);
       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
       eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
     std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
       if( !(muonPreColl.size()>=2 && electronPreColl.size()>=1) ) return;
     /**********************************************************************************************************/

     std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);

       eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
       eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
       eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(3.);
       eventbase->GetMuonSel()->SetBSdz(0.1);
       eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
     std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);
       eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
       eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
       eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(3.);
       eventbase->GetMuonSel()->SetBSdz(0.1);
       eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.1);
     std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,true);
     std::vector<snu::KMuon> muonColl;      
       if     ( isData){ if(k_running_nonprompt){ muonColl=muonLooseColl;} else{ muonColl=muonTightColl;} }
       else if(!isData) muonColl=muonTightColl;
//       else if(!isData) muonColl=SkimLepColl(muonTightColl, truthColl, "PromptEWtauNHConv");
  
       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
       eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
       eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.5, 0.5);
       eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.2);
       eventbase->GetElectronSel()->SetApplyConvVeto(true);    eventbase->GetElectronSel()->SetCheckCharge(true);
     std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP80);
       eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
       eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.05, 0.05);
       eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.2);
       eventbase->GetElectronSel()->SetApplyConvVeto(true);    eventbase->GetElectronSel()->SetCheckCharge(true);
     std::vector<snu::KElectron> electronTightColl; eventbase->GetElectronSel()->Selection(electronTightColl);
     std::vector<snu::KElectron> electronColl;
       if     ( isData){ if(k_running_nonprompt){ electronColl=electronLooseColl;} else{ electronColl=electronTightColl;} }
       else if(!isData) electronColl=electronTightColl;
//       else if(!isData) electronColl=SkimLepColl(electronTightColl, truthColl, "PromptEWtauNHConv");

  
       eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
       eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       //eventbase->GetJetSel()->SetPileUpJetID(true,"Loose");
       bool LeptonVeto=true;
     std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonLooseColl, electronLooseColl);
  
     std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");
     std::vector<snu::KJet> ljetColl = SelLightJets(jetColl, "Medium");
  
  
     double met = eventbase->GetEvent().PFMETType1();
     double met_x = eventbase->GetEvent().PFMETType1x();
     double met_y = eventbase->GetEvent().PFMETType1y();
     int nbjets=bjetColl.size(); int njets=jetColl.size(); const int nljets=ljetColl.size();
  
     /*****************************************************
     **Scale Factors
     *****************************************************/
     float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1., btag_sf=1.;
     float trigger_sf=1.;
     float fake_weight=1.; bool EventCand=false;
  
     /*This part is for boosting up speed.. SF part takes rather longer time than expected*/
     if(muonLooseColl.size()==2 && electronLooseColl.size()==1) EventCand=true; 
  
     if(EventCand){
  
       if(!isData){
        //trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );
    
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_80", electronColl);
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
    
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muonColl);
         //iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
         trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
    
         btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

         if(k_sample_name.Contains("ZZTo4L")) weight*=GetKFactor();
       }
       else{
         //Perfect Prompt Ratio Approximation applied.(p=1), Applicable to generic number, combination of leptons under premise of high prompt rate.
         if(k_running_nonprompt){
           fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", muonLooseColl.size(), electronLooseColl, "ELECTRON_MVA_TIGHT", electronLooseColl.size(), "ELECTRON_MVA_FAKELOOSE", "dijet_ajet40");
         }
       }
     }
  
     weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight*fake_weight*btag_sf;
  
     if( !(electronLooseColl.size()==1 && muonLooseColl.size()==2) ) return;
     if( !(electronColl.size()==1 && muonColl.size()==2) ) return;
     if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>15 && muonColl.at(1).Pt()>10) ) return;
     if( SumCharge(muonColl)!=0 ) return;
  
     FillCutFlow("3lCut", weight);

     float Mmumu=(muonColl.at(0)+muonColl.at(1)).M();
     float M3l=(electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M();
     float MTW = sqrt( 2*(met*electronColl.at(0).Pt()-met_x*electronColl.at(0).Px()-met_y*electronColl.at(0).Py()) );

     if( Mmumu<12 ) return;
     if( fabs(Mmumu-91.2)<10 && M3l>101.2 && met>50 ){
       if(MTW>50){
         FillHist("WZ_Yield", 0., weight, 0., 1., 1);
       }
       FillHist("WZ_Mmumu", Mmumu, weight, 60., 120., 60);
       FillHist("WZ_MET", met, weight, 0., 200., 200);
       FillHist("WZ_MTW", MTW, weight, 0., 200., 200);
     }

     if( fabs(M3l-91.2)<10 && Mmumu<81.2 && met<50 ){
       FillHist("ConvSel_M3l", M3l, weight, 60., 120., 60);
       if(k_sample_name.Contains("ZGto2LG")){
         int EleType=GetLeptonType(electronColl.at(0), truthColl);
         if     (EleType== 4) FillHist("ConvSel_M3l_Type4", M3l, weight, 60., 120., 60);
         else if(EleType== 5) FillHist("ConvSel_M3l_Type5", M3l, weight, 60., 120., 60);
         else if(EleType<=-5) FillHist("ConvSel_M3l_TypeM56", M3l, weight, 60., 120., 60);
       }
     }

     if(k_sample_name.Contains("ZGto2LG")){

       int EleType=GetLeptonType(electronColl.at(0), truthColl);
       int Mu1Type=GetLeptonType(muonColl.at(0), truthColl);
       int Mu2Type=GetLeptonType(muonColl.at(1), truthColl);
       int NEWAllTau=(EleType>0 && EleType<4) + (Mu1Type>0 && Mu1Type<4) + (Mu2Type>0 && Mu2Type<4);
       int NEWlep=(EleType>0 && EleType<3) + (Mu1Type>0 && Mu1Type<3) + (Mu2Type>0 && Mu2Type<3);
       int NEWtau=NEWAllTau-NEWlep;

//       cout<<"NEWAllTau "<<NEWAllTau<<" NEWlep "<<NEWlep<<endl;
//       cout<<"ElType "<<EleType<<" Mu1Type "<<Mu1Type<<" Mu2Type "<<Mu2Type<<endl;

       //yield comparison
       FillHist("ZGLepComposition", 0., weight, 0., 10., 10);//All
       if(NEWAllTau==2)   FillHist("ZGLepComposition", 1., weight, 0., 10., 10);//2 prompt or tau
       if(NEWlep==2   )   FillHist("ZGLepComposition", 2., weight, 0., 10., 10);//2 prompt
       if(EleType<0 || EleType>3) FillHist("ZGLepComposition", 3., weight, 0., 10., 10);//Ele fake
       if( (EleType<0 && EleType>=-4) || (Mu1Type<0 && Mu1Type>=-4) || (Mu2Type<0 && Mu2Type>=-4) ){//at least 1 hadron fake
         FillHist("ZGLepComposition", 4., weight, 0., 10., 10);
       }
       if( EleType<=-5 || EleType>3 ) FillHist("ZGLepComposition", 5., weight, 0., 10., 10);//ele is conversion fake
       if( EleType<=-5 ) FillHist("ZGLepComposition", 6., weight, 0., 10., 10);//Ele is external conversion
       if( EleType>=4 )  FillHist("ZGLepComposition", 7., weight, 0., 10., 10);//Ele is internal conversion
     }

   }



   bool SRContributionCheck=true;
   if(SRContributionCheck){

     if(!isData) weight*=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
     if(!isData) weight*=MCweight;
     float pileup_reweight=1.;
     if(!k_isdata) { pileup_reweight*=eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central);}
     FillCutFlow("NoCut", weight*pileup_reweight);


     int Pass_Trigger1=0, Pass_Trigger2=0;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ) Pass_Trigger1++;
     if( PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger2++;

     float trigger_period_weight=1.;
     bool Pass_Trigger=false;
     if(isData){
       int DataPeriod=GetDataPeriod();
       if( DataPeriod>0 && DataPeriod<7 && Pass_Trigger1==1 ) Pass_Trigger=true;
       else if( DataPeriod==7 && Pass_Trigger2==1 ) Pass_Trigger=true;
     }
     else{
       if( Pass_Trigger1>0 || Pass_Trigger2>0 ) Pass_Trigger=true;
       trigger_period_weight=(Pass_Trigger1*27.257618+Pass_Trigger2*8.605696)/35.863314;
     }
     weight*=trigger_period_weight;
     if(!Pass_Trigger) return;     FillCutFlow("TriggerCut", weight*pileup_reweight);
     if(!PassMETFilter()) return;  FillCutFlow("EventCut", weight*pileup_reweight);
     if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; FillCutFlow("VertexCut", weight*pileup_reweight);

     /**PreSelCut***********************************************************************************************/
     //Intended for Code speed boosting up.
       eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, true);
       eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
     std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
       if( !(muonPreColl.size()>=2 && electronPreColl.size()>=1) ) return;
     /**********************************************************************************************************/

     std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);

       eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
       eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
       eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(3.);
       eventbase->GetMuonSel()->SetBSdz(0.1);
       eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
     std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);
       eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
       eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
       eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(3.);
       eventbase->GetMuonSel()->SetBSdz(0.1);
       eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.1);
     std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,true);
     std::vector<snu::KMuon> muonColl;      
       if     ( isData){ if(k_running_nonprompt){ muonColl=muonLooseColl;} else{ muonColl=muonTightColl;} }
       else if(!isData) muonColl=muonTightColl;
//       else if(!isData) muonColl=SkimLepColl(muonTightColl, truthColl, "PromptEWtauNHConv");
  

       eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     std::vector<snu::KElectron> electronPreLooseColl; eventbase->GetElectronSel()->Selection(electronPreLooseColl);
       
     std::vector<snu::KElectron> electronLooseColl; 
       for(int i=0; i<electronPreLooseColl.size(); i++){
         if(PassIDCriteria(electronPreLooseColl.at(i), "POGMVAMFakeLIso04Opt1")) electronLooseColl.push_back(electronPreLooseColl.at(i));
         //if(PassIDCriteria(electronPreLooseColl.at(i), "POGMVAMFakeLIso04Opt2")) electronLooseColl.push_back(electronPreLooseColl.at(i));
       }
     
       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
       eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
       eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetdxySigMax(3.);
       eventbase->GetElectronSel()->SetApplyConvVeto(true);
     std::vector<snu::KElectron> electronTightColl; eventbase->GetElectronSel()->Selection(electronTightColl);
     std::vector<snu::KElectron> electronColl=electronTightColl;
//     std::vector<snu::KElectron> electronColl;
//       if     ( isData){ if(k_running_nonprompt){ electronColl=electronLooseColl;} else{ electronColl=electronTightColl;} }
//       else if(!isData) electronColl=electronTightColl;
//       else if(!isData) electronColl=SkimLepColl(electronTightColl, truthColl, "PromptEWtauNHConv");

  
       eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
       eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       //eventbase->GetJetSel()->SetPileUpJetID(true,"Loose");
       bool LeptonVeto=true;
     std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonTightColl, electronTightColl);
  
     std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");
     std::vector<snu::KJet> ljetColl = SelLightJets(jetColl, "Medium");
  
  
     double met = eventbase->GetEvent().PFMETType1();
     double met_x = eventbase->GetEvent().PFMETType1x();
     double met_y = eventbase->GetEvent().PFMETType1y();
     int nbjets=bjetColl.size(); int njets=jetColl.size(); const int nljets=ljetColl.size();
  
     /*****************************************************
     **Scale Factors
     *****************************************************/
     float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1., btag_sf=1.;
     float trigger_sf=1.;
     float fake_weight=1.; bool EventCand=false;
  
     /*This part is for boosting up speed.. SF part takes rather longer time than expected*/
     if(muonLooseColl.size()==2 && electronLooseColl.size()==1) EventCand=true; 
  
     if(EventCand){
  
       if(!isData){
        //trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );
    
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electronColl);
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
    
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muonColl);
         //iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
         trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
    
         btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

         if(k_sample_name.Contains("ZZTo4L")) weight*=GetKFactor();
       }
       else{
         //Perfect Prompt Ratio Approximation applied.(p=1), Applicable to generic number, combination of leptons under premise of high prompt rate.
         if(k_running_nonprompt){
           fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", muonLooseColl.size(), electronLooseColl, "ELECTRON_MVA_TIGHT", electronLooseColl.size(), "ELECTRON_MVA_FAKELOOSE", "dijet_ajet40");
         }
       }
     }
  
     weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight*fake_weight*btag_sf;
  
     if( !(electronLooseColl.size()==1 && muonLooseColl.size()==2) ) return;
     if( !(electronColl.size()==1 && muonColl.size()==2) ) return;
     if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>15 && muonColl.at(1).Pt()>10) ) return;
     if( SumCharge(muonColl)!=0 ) return;
     FillCutFlow("3lOSCut", weight);

     float Mmumu=(muonColl.at(0)+muonColl.at(1)).M();
     float M3l=(electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M();
     float MTW = sqrt( 2*(met*electronColl.at(0).Pt()-met_x*electronColl.at(0).Px()-met_y*electronColl.at(0).Py()) );

     if( Mmumu<12 ) return;
     FillCutFlow("M(#mu#mu)>12", weight);

     if(njets<3) return;
     FillCutFlow("#geq3j", weight);

     if(nbjets==0) return;
     FillCutFlow("#geq1b", weight);

     if( fabs(Mmumu-91.2)<10 ) return;
     FillCutFlow("ZVeto", weight);

     if( Mmumu>40 ) return;
     FillCutFlow("M(A)Range", weight);

   }


/////////////////////////////////////////////////////////////////////////////////// 


return;
}// End of execute event loop
  


void Jun2017_ConversionStudy::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}

/*
bool Jun2017_ConversionStudy::IsFinalPhotonSt23(std::vector<snu::KTruth> TruthColl){
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

int Jun2017_ConversionStudy::GetPhotonType(int PhotonIdx, std::vector<snu::KTruth> TruthColl){
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

bool Jun2017_ConversionStudy::IsHardPhotonConverted(std::vector<snu::KTruth> TruthColl){

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

void Jun2017_ConversionStudy::BeginCycle() throw( LQError ){
  
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

Jun2017_ConversionStudy::~Jun2017_ConversionStudy() {
  
  Message("In Jun2017_ConversionStudy Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Jun2017_ConversionStudy::FillCutFlow(TString cut, float weight){
  
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



void Jun2017_ConversionStudy::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Jun2017_ConversionStudy::MakeHistograms(){
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
  // *  Remove//Overide this Jun2017_ConversionStudyCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Jun2017_ConversionStudy::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
