#include "MCDataCorrections.h"

float MCDataCorrections::GetTriggerSF(vector<snu::KElectron>& EleColl, vector<snu::KMuon>& MuColl, TString TrigName, TString Option){

  if(corr_isdata) return 1.;

  float TriggerEff_Data = TriggerEfficiency(EleColl, MuColl, TrigName, true,  Option);
  float TriggerEff_MC   = TriggerEfficiency(EleColl, MuColl, TrigName, false, Option);
  
  float TriggerScaleFactor = TriggerEff_MC!=0.? TriggerEff_Data/TriggerEff_MC:0.;
   if(TriggerScaleFactor<0) TriggerScaleFactor=0.;

  return TriggerScaleFactor;

}


float MCDataCorrections::TriggerEfficiency(vector<snu::KElectron>& EleColl, vector<snu::KMuon>& MuColl, TString TrigName, bool ReturnDataEff, TString Option){
  //DataorMC : T: Return DataEff, F: Return MCEff
  //Option : NoCorr  : Each Leg Perfect Independent(No Correlation)(One lep failing a leg can fire same flav leg without any Bias) 
  //                   P(leg2|lep1 && !leg1|lep1)=P(leg2|lep1)=Measured leg eff.
  //         FullCorr: Each Leg Perfect Correlation(One lep failing a leg cannot fire any same flav leg)
  //                   P(leg2|lep1 && !leg1|lep1)=0

  if(corr_isdata) return 1.;

  bool NoCorr=Option.Contains("NoCorr"), FullCorr=Option.Contains("FullCorr");
  TString StrMCorData = ReturnDataEff? "_DATA":"_MC";
  float TriggerEff=-1.;
  float SystDir=0;
  bool Syst_Leg1=false, Syst_Leg2=false, Syst_DZ=false;
  if(Option.Contains("Syst")){
    if     (Option.Contains("Up"))   SystDir= 1.;
    else if(Option.Contains("Down")) SystDir=-1.;
    if     (Option.Contains("Leg1")) Syst_Leg1=true;
    if     (Option.Contains("Leg2")) Syst_Leg2=true;
    if     (Option.Contains("DZ"))   Syst_DZ  =true;
  }


  if(TrigName.Contains("HLT_Mu17_TrkIsoVVL_v")){
    float MinPt=10., MaxPt=120.;
    if(MuColl.size()==1){
      float mu1pt = MuColl.at(0).Pt(), mu1feta = fabs(MuColl.at(0).Eta());
      if     (mu1pt<MinPt) mu1pt=MinPt;
      else if(mu1pt>MaxPt) mu1pt=MaxPt-1;

      TH2F* HistEff_LegMu = GetCorrectionHist("TRIGEFF_DiMu_Mu17_HctoWA"+StrMCorData);

      float Eff_Mu1LegMu = HistEff_LegMu->GetBinContent(HistEff_LegMu->FindBin(mu1feta, mu1pt));

      TriggerEff = Eff_Mu1LegMu;

    }//End of 2Mu
    if(MuColl.size()==2){
      float mu1pt = MuColl.at(0).Pt(), mu1feta = fabs(MuColl.at(0).Eta());
      float mu2pt = MuColl.at(1).Pt(), mu2feta = fabs(MuColl.at(1).Eta());
      if     (mu1pt<MinPt) mu1pt=MinPt;
      else if(mu1pt>MaxPt) mu1pt=MaxPt-1;
      if     (mu2pt<MinPt) mu2pt=MinPt;
      else if(mu2pt>MaxPt) mu2pt=MaxPt-1;

      TH2F* HistEff_LegMu = GetCorrectionHist("TRIGEFF_DiMu_Mu17_HctoWA"+StrMCorData);

      float Eff_Mu1LegMu = HistEff_LegMu->GetBinContent(HistEff_LegMu->FindBin(mu1feta, mu1pt));
      float Eff_Mu2LegMu = HistEff_LegMu->GetBinContent(HistEff_LegMu->FindBin(mu2feta, mu2pt));

      TriggerEff = 1.-(1.-Eff_Mu1LegMu)*(1.-Eff_Mu2LegMu);

    }//End of 2Mu
  }
  else if(TrigName.Contains("HLT_Mu8_TrkIsoVVL_v")){
    float MinPt=10., MaxPt=120.;
    if(MuColl.size()==1){
      float mu1pt = MuColl.at(0).Pt(), mu1feta = fabs(MuColl.at(0).Eta());
      if     (mu1pt<MinPt) mu1pt=MinPt;
      else if(mu1pt>MaxPt) mu1pt=MaxPt-1;

      TH2F* HistEff_LegMu = GetCorrectionHist("TRIGEFF_EMu_Mu8_HctoWA"+StrMCorData);

      float Eff_Mu1LegMu = HistEff_LegMu->GetBinContent(HistEff_LegMu->FindBin(mu1pt, mu1feta));

      TriggerEff = Eff_Mu1LegMu;

    }//End of 2Mu
    if(MuColl.size()==2){
      float mu1pt = MuColl.at(0).Pt(), mu1feta = fabs(MuColl.at(0).Eta());
      float mu2pt = MuColl.at(1).Pt(), mu2feta = fabs(MuColl.at(1).Eta());
      if     (mu1pt<MinPt) mu1pt=MinPt;
      else if(mu1pt>MaxPt) mu1pt=MaxPt-1;
      if     (mu2pt<MinPt) mu2pt=MinPt;
      else if(mu2pt>MaxPt) mu2pt=MaxPt-1;

      TH2F* HistEff_LegMu = GetCorrectionHist("TRIGEFF_EMu_Mu8_HctoWA"+StrMCorData);

      float Eff_Mu1LegMu = HistEff_LegMu->GetBinContent(HistEff_LegMu->FindBin(mu1pt, mu1feta));
      float Eff_Mu2LegMu = HistEff_LegMu->GetBinContent(HistEff_LegMu->FindBin(mu2pt, mu2feta));

      TriggerEff = 1.-(1.-Eff_Mu1LegMu)*(1.-Eff_Mu2LegMu);

    }//End of 2Mu
  }
  else if(TrigName.Contains("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v")){
    float MinPt=25., MaxPt=200.;
    if(EleColl.size()==1){
      float elpt = EleColl.at(0).Pt(), elfeta = fabs(EleColl.at(0).Eta());
      if     (elpt<MinPt) elpt=MinPt;
      else if(elpt>MaxPt) elpt=MaxPt-1;

      TH2F* HistEff_LegEl = GetCorrectionHist("TRIGEFF_EMu_Ele23_HctoWA"+StrMCorData);

      float Eff_El1LegEl = HistEff_LegEl->GetBinContent(HistEff_LegEl->FindBin(elpt, elfeta));

      TriggerEff = Eff_El1LegEl;
    }//End of 2Mu
    if(EleColl.size()==2){
      float el1pt = EleColl.at(0).Pt(), el1feta = fabs(EleColl.at(0).Eta());
      float el2pt = EleColl.at(1).Pt(), el2feta = fabs(EleColl.at(1).Eta());
      if     (el1pt<MinPt) el1pt=MinPt;
      else if(el1pt>MaxPt) el1pt=MaxPt-1;
      if     (el2pt<MinPt) el2pt=MinPt;
      else if(el2pt>MaxPt) el2pt=MaxPt-1;

      TH2F* HistEff_LegEl = GetCorrectionHist("TRIGEFF_EMu_Ele23_HctoWA"+StrMCorData);

      float Eff_El1LegEl = HistEff_LegEl->GetBinContent(HistEff_LegEl->FindBin(el1pt, el1feta));
      float Eff_El2LegEl = HistEff_LegEl->GetBinContent(HistEff_LegEl->FindBin(el2pt, el2feta));

      TriggerEff = 1.-(1.-Eff_El1LegEl)*(1.-Eff_El2LegEl);
    }//End of 2Mu
  }
  else if(TrigName.Contains("HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL")){
    float MinPt=10., MaxPt=120.;
    if(MuColl.size()==2){
      float mu1pt = MuColl.at(0).Pt(), mu1feta = fabs(MuColl.at(0).Eta());
      float mu2pt = MuColl.at(1).Pt(), mu2feta = fabs(MuColl.at(1).Eta());
      if     (mu1pt<MinPt) mu1pt=MinPt;
      else if(mu1pt>MaxPt) mu1pt=MaxPt-1;
      if     (mu2pt<MinPt) mu2pt=MinPt;
      else if(mu2pt>MaxPt) mu2pt=MaxPt-1;

      TH2F* HistEff_Leg1 = GetCorrectionHist("TRIGEFF_DiMu_Mu17_HctoWA"+StrMCorData);
      TH2F* HistEff_Leg2 = GetCorrectionHist("TRIGEFF_DiMu_Mu8ORTkMu8_HctoWA"+StrMCorData);

      float Eff_Mu1Leg1 = HistEff_Leg1->GetBinContent(HistEff_Leg1->FindBin(mu1feta, mu1pt));
      float Eff_Mu1Leg2 = HistEff_Leg2->GetBinContent(HistEff_Leg2->FindBin(mu1feta, mu1pt));
      float Eff_Mu2Leg1 = HistEff_Leg1->GetBinContent(HistEff_Leg1->FindBin(mu2feta, mu2pt));
      float Eff_Mu2Leg2 = HistEff_Leg2->GetBinContent(HistEff_Leg2->FindBin(mu2feta, mu2pt));
      float Eff_dz      = ReturnDataEff? 0.979389:0.989453;

      if(Syst_Leg1){ Eff_Mu1Leg1+=SystDir*HistEff_Leg1->GetBinError(HistEff_Leg1->FindBin(mu1feta, mu1pt));
                     Eff_Mu2Leg1+=SystDir*HistEff_Leg1->GetBinError(HistEff_Leg1->FindBin(mu2feta, mu2pt)); }
      if(Syst_Leg2){ Eff_Mu1Leg2+=SystDir*HistEff_Leg2->GetBinError(HistEff_Leg2->FindBin(mu1feta, mu1pt));
                     Eff_Mu2Leg2+=SystDir*HistEff_Leg2->GetBinError(HistEff_Leg2->FindBin(mu2feta, mu2pt)); }
      if(Syst_DZ  ){ Eff_dz     +=SystDir*( ReturnDataEff? 0.00269428:0.00397434); }

      if(!TrigName.Contains("DZ")) Eff_dz=1.;

      if(NoCorr){
        TriggerEff = (1.-(1.-Eff_Mu1Leg1*Eff_Mu2Leg2)*(1.-Eff_Mu2Leg1*Eff_Mu1Leg2))*Eff_dz;
      }
      else if(FullCorr){
        TriggerEff = Eff_Mu1Leg1*Eff_Mu2Leg2*Eff_dz;
      }

    }//End of 2Mu
    else if(MuColl.size()==3){
      float mu1pt = MuColl.at(0).Pt(), mu1feta = fabs(MuColl.at(0).Eta());
      float mu2pt = MuColl.at(1).Pt(), mu2feta = fabs(MuColl.at(1).Eta());
      float mu3pt = MuColl.at(2).Pt(), mu3feta = fabs(MuColl.at(2).Eta());
      if     (mu1pt<MinPt) mu1pt=MinPt;
      else if(mu1pt>MaxPt) mu1pt=MaxPt-1;
      if     (mu2pt<MinPt) mu2pt=MinPt;
      else if(mu2pt>MaxPt) mu2pt=MaxPt-1;
      if     (mu3pt<MinPt) mu3pt=MinPt;
      else if(mu3pt>MaxPt) mu3pt=MaxPt-1;

      TH2F* HistEff_Leg1 = GetCorrectionHist("TRIGEFF_DiMu_Mu17_HctoWA"+StrMCorData);
      TH2F* HistEff_Leg2 = GetCorrectionHist("TRIGEFF_DiMu_Mu8ORTkMu8_HctoWA"+StrMCorData);

      float Eff_Mu1Leg1 = HistEff_Leg1->GetBinContent(HistEff_Leg1->FindBin(mu1feta, mu1pt));
      float Eff_Mu1Leg2 = HistEff_Leg2->GetBinContent(HistEff_Leg2->FindBin(mu1feta, mu1pt));
      float Eff_Mu2Leg1 = HistEff_Leg1->GetBinContent(HistEff_Leg1->FindBin(mu2feta, mu2pt));
      float Eff_Mu2Leg2 = HistEff_Leg2->GetBinContent(HistEff_Leg2->FindBin(mu2feta, mu2pt));
      float Eff_Mu3Leg1 = HistEff_Leg1->GetBinContent(HistEff_Leg1->FindBin(mu3feta, mu3pt));
      float Eff_Mu3Leg2 = HistEff_Leg2->GetBinContent(HistEff_Leg2->FindBin(mu3feta, mu3pt));
      float Eff_dz      = ReturnDataEff? 0.979389:0.989453;

      if(Syst_Leg1){ Eff_Mu1Leg1+=SystDir*HistEff_Leg1->GetBinError(HistEff_Leg1->FindBin(mu1feta, mu1pt));
                     Eff_Mu2Leg1+=SystDir*HistEff_Leg1->GetBinError(HistEff_Leg1->FindBin(mu2feta, mu2pt)); 
                     Eff_Mu3Leg1+=SystDir*HistEff_Leg1->GetBinError(HistEff_Leg1->FindBin(mu3feta, mu3pt)); }
      if(Syst_Leg2){ Eff_Mu1Leg2+=SystDir*HistEff_Leg2->GetBinError(HistEff_Leg2->FindBin(mu1feta, mu1pt));
                     Eff_Mu2Leg2+=SystDir*HistEff_Leg2->GetBinError(HistEff_Leg2->FindBin(mu2feta, mu2pt)); 
                     Eff_Mu3Leg2+=SystDir*HistEff_Leg2->GetBinError(HistEff_Leg2->FindBin(mu3feta, mu3pt)); }
      if(Syst_DZ  ){ Eff_dz     +=SystDir*( ReturnDataEff? 0.00269428:0.00397434); }

      if(!TrigName.Contains("DZ")) Eff_dz=1.;


      if(NoCorr){
        TriggerEff=( Eff_Mu1Leg1*(1.-(1.-Eff_Mu2Leg2)*(1.-Eff_Mu3Leg2))
                    +(1.-Eff_Mu1Leg1)*Eff_Mu1Leg2*(1.-(1.-Eff_Mu2Leg1)*(1.-Eff_Mu3Leg1))
                    +(1.-Eff_Mu1Leg1)*(1.-Eff_Mu1Leg2)*(1.-(1.-Eff_Mu2Leg1*Eff_Mu3Leg2)*(1.-Eff_Mu3Leg1*Eff_Mu2Leg2)) )*Eff_dz;
      }
      else if(FullCorr){
        TriggerEff= Eff_Mu1Leg1*(1.-(1.-Eff_Mu2Leg2*Eff_dz)*(1.-Eff_Mu3Leg2*Eff_dz))
                   +(1.-Eff_Mu1Leg1)*Eff_Mu2Leg1*Eff_Mu3Leg2*Eff_dz;
      }
    }//End of 3Mu
    else if(MuColl.size()==4){
      float mu1pt = MuColl.at(0).Pt(), mu1feta = fabs(MuColl.at(0).Eta());
      float mu2pt = MuColl.at(1).Pt(), mu2feta = fabs(MuColl.at(1).Eta());
      float mu3pt = MuColl.at(2).Pt(), mu3feta = fabs(MuColl.at(2).Eta());
      float mu4pt = MuColl.at(3).Pt(), mu4feta = fabs(MuColl.at(3).Eta());
      if     (mu1pt<MinPt) mu1pt=MinPt;
      else if(mu1pt>MaxPt) mu1pt=MaxPt-1;
      if     (mu2pt<MinPt) mu2pt=MinPt;
      else if(mu2pt>MaxPt) mu2pt=MaxPt-1;
      if     (mu3pt<MinPt) mu3pt=MinPt;
      else if(mu3pt>MaxPt) mu3pt=MaxPt-1;
      if     (mu4pt<MinPt) mu4pt=MinPt;
      else if(mu4pt>MaxPt) mu4pt=MaxPt-1;

      TH2F* HistEff_Leg1 = GetCorrectionHist("TRIGEFF_DiMu_Mu17_HctoWA"+StrMCorData);
      TH2F* HistEff_Leg2 = GetCorrectionHist("TRIGEFF_DiMu_Mu8ORTkMu8_HctoWA"+StrMCorData);

      float Eff_Mu1Leg1 = HistEff_Leg1->GetBinContent(HistEff_Leg1->FindBin(mu1feta, mu1pt));
      float Eff_Mu2Leg1 = HistEff_Leg1->GetBinContent(HistEff_Leg1->FindBin(mu2feta, mu2pt));
      float Eff_Mu2Leg2 = HistEff_Leg2->GetBinContent(HistEff_Leg2->FindBin(mu2feta, mu2pt));
      float Eff_Mu3Leg1 = HistEff_Leg1->GetBinContent(HistEff_Leg1->FindBin(mu3feta, mu3pt));
      float Eff_Mu3Leg2 = HistEff_Leg2->GetBinContent(HistEff_Leg2->FindBin(mu3feta, mu3pt));
      float Eff_Mu4Leg2 = HistEff_Leg2->GetBinContent(HistEff_Leg2->FindBin(mu4feta, mu4pt));
      float Eff_dz      = ReturnDataEff? 0.979389:0.989453;

      if(Syst_Leg1){ Eff_Mu1Leg1+=SystDir*HistEff_Leg1->GetBinError(HistEff_Leg1->FindBin(mu1feta, mu1pt));
                     Eff_Mu2Leg1+=SystDir*HistEff_Leg1->GetBinError(HistEff_Leg1->FindBin(mu2feta, mu2pt));
                     Eff_Mu3Leg1+=SystDir*HistEff_Leg1->GetBinError(HistEff_Leg1->FindBin(mu3feta, mu3pt)); }
      if(Syst_Leg2){ Eff_Mu2Leg2+=SystDir*HistEff_Leg2->GetBinError(HistEff_Leg2->FindBin(mu2feta, mu2pt));
                     Eff_Mu3Leg2+=SystDir*HistEff_Leg2->GetBinError(HistEff_Leg2->FindBin(mu3feta, mu3pt));
                     Eff_Mu4Leg2+=SystDir*HistEff_Leg2->GetBinError(HistEff_Leg2->FindBin(mu4feta, mu4pt)); }
      if(Syst_DZ  ){ Eff_dz     +=SystDir*( ReturnDataEff? 0.00269428:0.00397434); }

      if(!TrigName.Contains("DZ")) Eff_dz=1.;


      //Full Correlation Formulae: MC closure only full correlation agrees with obs.
      TriggerEff= Eff_Mu1Leg1*(1.-(1.-Eff_Mu2Leg2*Eff_dz)*(1.-Eff_Mu3Leg2*Eff_dz)*(1.-Eff_Mu4Leg2*Eff_dz))
                 +(1.-Eff_Mu1Leg1)*Eff_Mu2Leg1*(1.-(1.-Eff_Mu3Leg2*Eff_dz)*(1.-Eff_Mu4Leg2*Eff_dz))
                 +(1.-Eff_Mu1Leg1)*(1.-Eff_Mu2Leg1)*Eff_Mu3Leg1*Eff_Mu4Leg2*Eff_dz;
    }//End of 4Mu
  }//End of DiMu Trigger
  else if(TrigName.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL")){
    float MinPt_Ele=25., MaxPt_Ele=200.;
    float MinPt_Mu =10., MaxPt_Mu =120.;
    if(MuColl.size()==1 && EleColl.size()==1){
      float mupt = MuColl.at(0).Pt(),  mufeta = fabs(MuColl.at(0).Eta());
      float elpt = EleColl.at(0).Pt(), elfeta = fabs(EleColl.at(0).Eta());
      if     (mupt<MinPt_Mu)  mupt=MinPt_Mu;
      else if(mupt>MaxPt_Mu)  mupt=MaxPt_Mu-1;
      if     (elpt<MinPt_Ele) elpt=MinPt_Ele;
      else if(elpt>MaxPt_Ele) elpt=MaxPt_Ele-1;


      TH2F* HistEff_LegMu = GetCorrectionHist("TRIGEFF_EMu_Mu8_HctoWA"+StrMCorData);
      TH2F* HistEff_LegEl = GetCorrectionHist("TRIGEFF_EMu_Ele23_HctoWA"+StrMCorData);

      float Eff_LegMu = HistEff_LegMu->GetBinContent(HistEff_LegMu->FindBin(mupt, mufeta));
      float Eff_LegEl = HistEff_LegEl->GetBinContent(HistEff_LegEl->FindBin(elpt, elfeta));
      float Eff_dz    = ReturnDataEff? 0.962337:0.991454;

      if(Syst_Leg1){ Eff_LegMu+=SystDir*HistEff_LegMu->GetBinError(HistEff_LegMu->FindBin(mupt, mufeta));  }
      if(Syst_Leg2){ Eff_LegEl+=SystDir*HistEff_LegEl->GetBinError(HistEff_LegEl->FindBin(elpt , elfeta)); }
      if(Syst_DZ  ){ Eff_dz   +=SystDir*( ReturnDataEff? 0.00500431:0.00239173 ); }

      if(!TrigName.Contains("DZ")) Eff_dz=1.;

      TriggerEff = Eff_LegMu*Eff_LegEl*Eff_dz;

    }
    else if(MuColl.size()==2 && EleColl.size()==1){
      float mu1pt = MuColl.at(0).Pt(),  mu1feta = fabs(MuColl.at(0).Eta());
      float mu2pt = MuColl.at(1).Pt(),  mu2feta = fabs(MuColl.at(1).Eta());
      float elpt  = EleColl.at(0).Pt(), elfeta  = fabs(EleColl.at(0).Eta());

      if     (mu1pt<MinPt_Mu) mu1pt=MinPt_Mu;
      else if(mu1pt>MaxPt_Mu) mu1pt=MaxPt_Mu-1;
      if     (mu2pt<MinPt_Mu) mu2pt=MinPt_Mu;
      else if(mu2pt>MaxPt_Mu) mu2pt=MaxPt_Mu-1;
      if     (elpt<MinPt_Ele) elpt =MinPt_Ele;
      else if(elpt>MaxPt_Ele) elpt =MaxPt_Ele-1;

      TH2F* HistEff_LegMu = GetCorrectionHist("TRIGEFF_EMu_Mu8_HctoWA"+StrMCorData);
      TH2F* HistEff_LegEl = GetCorrectionHist("TRIGEFF_EMu_Ele23_HctoWA"+StrMCorData);

      float Eff_Mu1LegMu = HistEff_LegMu->GetBinContent(HistEff_LegMu->FindBin(mu1pt, mu1feta));
      float Eff_Mu2LegMu = HistEff_LegMu->GetBinContent(HistEff_LegMu->FindBin(mu2pt, mu2feta));
      float Eff_LegEl    = HistEff_LegEl->GetBinContent(HistEff_LegEl->FindBin(elpt , elfeta ));
      float Eff_dz       = ReturnDataEff? 0.962337:0.991454;

      if(Syst_Leg1){ Eff_Mu1LegMu+=SystDir*HistEff_LegMu->GetBinError(HistEff_LegMu->FindBin(mu1pt, mu1feta));
                     Eff_Mu2LegMu+=SystDir*HistEff_LegMu->GetBinError(HistEff_LegMu->FindBin(mu2pt, mu2feta)); }
      if(Syst_Leg2){ Eff_LegEl   +=SystDir*HistEff_LegEl->GetBinError(HistEff_LegEl->FindBin(elpt , elfeta )); }
      if(Syst_DZ  ){ Eff_dz      +=SystDir*( ReturnDataEff? 0.00500431:0.00239173 ); }

      if(!TrigName.Contains("DZ")) Eff_dz=1.;

      TriggerEff = Eff_LegEl*(1.-(1.-Eff_Mu1LegMu*Eff_dz)*(1.-Eff_Mu2LegMu*Eff_dz));

    }
    else if(MuColl.size()==3 && EleColl.size()==1){
      float mu1pt = MuColl.at(0).Pt(),  mu1feta = fabs(MuColl.at(0).Eta());
      float mu2pt = MuColl.at(1).Pt(),  mu2feta = fabs(MuColl.at(1).Eta());
      float mu3pt = MuColl.at(2).Pt(),  mu3feta = fabs(MuColl.at(2).Eta());
      float elpt  = EleColl.at(0).Pt(), elfeta  = fabs(EleColl.at(0).Eta());

      if     (mu1pt<MinPt_Mu) mu1pt=MinPt_Mu;
      else if(mu1pt>MaxPt_Mu) mu1pt=MaxPt_Mu-1;
      if     (mu2pt<MinPt_Mu) mu2pt=MinPt_Mu;
      else if(mu2pt>MaxPt_Mu) mu2pt=MaxPt_Mu-1;
      if     (mu3pt<MinPt_Mu) mu3pt=MinPt_Mu;
      else if(mu3pt>MaxPt_Mu) mu3pt=MaxPt_Mu-1;
      if     (elpt<MinPt_Ele) elpt =MinPt_Ele;
      else if(elpt>MaxPt_Ele) elpt =MaxPt_Ele-1;

      TH2F* HistEff_LegMu = GetCorrectionHist("TRIGEFF_EMu_Mu8_HctoWA"+StrMCorData);
      TH2F* HistEff_LegEl = GetCorrectionHist("TRIGEFF_EMu_Ele23_HctoWA"+StrMCorData);

      float Eff_Mu1LegMu = HistEff_LegMu->GetBinContent(HistEff_LegMu->FindBin(mu1pt, mu1feta));
      float Eff_Mu2LegMu = HistEff_LegMu->GetBinContent(HistEff_LegMu->FindBin(mu2pt, mu2feta));
      float Eff_Mu3LegMu = HistEff_LegMu->GetBinContent(HistEff_LegMu->FindBin(mu3pt, mu3feta));
      float Eff_LegEl    = HistEff_LegEl->GetBinContent(HistEff_LegEl->FindBin(elpt , elfeta ));
      float Eff_dz       = ReturnDataEff? 0.962337:0.991454;

      if(Syst_Leg1){ Eff_Mu1LegMu+=SystDir*HistEff_LegMu->GetBinError(HistEff_LegMu->FindBin(mu1pt, mu1feta));
                     Eff_Mu2LegMu+=SystDir*HistEff_LegMu->GetBinError(HistEff_LegMu->FindBin(mu2pt, mu2feta));
                     Eff_Mu3LegMu+=SystDir*HistEff_LegMu->GetBinError(HistEff_LegMu->FindBin(mu3pt, mu3feta)); }
      if(Syst_Leg2){ Eff_LegEl   +=SystDir*HistEff_LegEl->GetBinError(HistEff_LegEl->FindBin(elpt , elfeta )); }
      if(Syst_DZ  ){ Eff_dz      +=SystDir*( ReturnDataEff? 0.00500431:0.00239173 ); }

      if(!TrigName.Contains("DZ")) Eff_dz=1.;

      TriggerEff = Eff_LegEl*(1.-(1.-Eff_Mu1LegMu*Eff_dz)*(1.-Eff_Mu2LegMu*Eff_dz)*(1.-Eff_Mu3LegMu*Eff_dz));

    }
    else if(MuColl.size()==2 && EleColl.size()==2){
      float mu1pt = MuColl.at(0).Pt(),  mu1feta = fabs(MuColl.at(0).Eta());
      float mu2pt = MuColl.at(1).Pt(),  mu2feta = fabs(MuColl.at(1).Eta());
      float el1pt = EleColl.at(0).Pt(), el1feta  = fabs(EleColl.at(0).Eta());
      float el2pt = EleColl.at(1).Pt(), el2feta  = fabs(EleColl.at(1).Eta());

      if     (mu1pt<MinPt_Mu)  mu1pt=MinPt_Mu;
      else if(mu1pt>MaxPt_Mu)  mu1pt=MaxPt_Mu-1;
      if     (mu2pt<MinPt_Mu)  mu2pt=MinPt_Mu;
      else if(mu2pt>MaxPt_Mu)  mu2pt=MaxPt_Mu-1;
      if     (el1pt<MinPt_Ele) el1pt=MinPt_Ele;
      else if(el1pt>MaxPt_Ele) el1pt=MaxPt_Ele-1;
      if     (el2pt<MinPt_Ele) el2pt=MinPt_Ele;
      else if(el2pt>MaxPt_Ele) el2pt=MaxPt_Ele-1;

      TH2F* HistEff_LegMu = GetCorrectionHist("TRIGEFF_EMu_Mu8_HctoWA"+StrMCorData);
      TH2F* HistEff_LegEl = GetCorrectionHist("TRIGEFF_EMu_Ele23_HctoWA"+StrMCorData);

      float Eff_Mu1LegMu = HistEff_LegMu->GetBinContent(HistEff_LegMu->FindBin(mu1pt, mu1feta));
      float Eff_Mu2LegMu = HistEff_LegMu->GetBinContent(HistEff_LegMu->FindBin(mu2pt, mu2feta));
      float Eff_El1LegEl = HistEff_LegEl->GetBinContent(HistEff_LegEl->FindBin(el1pt, el1feta));
      float Eff_El2LegEl = HistEff_LegEl->GetBinContent(HistEff_LegEl->FindBin(el2pt, el2feta));
      float Eff_dz       = ReturnDataEff? 0.962337:0.991454;

      if(Syst_Leg1){ Eff_Mu1LegMu+=SystDir*HistEff_LegMu->GetBinError(HistEff_LegMu->FindBin(mu1pt, mu1feta));
                     Eff_Mu2LegMu+=SystDir*HistEff_LegMu->GetBinError(HistEff_LegMu->FindBin(mu2pt, mu2feta)); }
      if(Syst_Leg1){ Eff_El1LegEl+=SystDir*HistEff_LegEl->GetBinError(HistEff_LegEl->FindBin(el1pt, el1feta));
                     Eff_El2LegEl+=SystDir*HistEff_LegEl->GetBinError(HistEff_LegEl->FindBin(el2pt, el2feta)); }
      if(Syst_DZ  ){ Eff_dz      +=SystDir*( ReturnDataEff? 0.00500431:0.00239173 ); }

      if(!TrigName.Contains("DZ")) Eff_dz=1.;

      TriggerEff = (1.-(1.-Eff_El1LegEl)*(1.-Eff_El2LegEl))*(1.-(1.-Eff_Mu1LegMu*Eff_dz)*(1.-Eff_Mu2LegMu*Eff_dz));

    }
  }

  return TriggerEff;
}
