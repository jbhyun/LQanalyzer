#include "AnalyzerCore.h"

float AnalyzerCore::GetFakeWeight_Data(std::vector<snu::KMuon>& MuLColl, std::vector<snu::KElectron>& EleLColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, TString Option){

  float fakeweight=-1.; int NLooseNotTight=0;

  TString FilterInfo="", ConeMethod="", SystOpt="", SystKind="";
  if     (Option.Contains("NoFilter"))  FilterInfo="NoFilter";
  else if(Option.Contains("TrkIsoVVL")) FilterInfo="TrkIsoVVL";
  if     (Option.Contains("ConeE"))     ConeMethod="ConeE";
  else if(Option.Contains("ConeSUSY"))  ConeMethod="ConeSUSY";
  else if(Option.Contains("FOPt"))      ConeMethod="FOPt";
  if     (Option.Contains("SystUp"))    SystOpt="SystUp";
  else if(Option.Contains("SystDown"))  SystOpt="SystDown";
  if     (Option.Contains("PrNorm"))    SystKind="PrN";
  else if(Option.Contains("MES"))       SystKind="MES";
  else if(Option.Contains("JFlav"))     SystKind="JFlav";
  

  for(int i=0; i<(int) MuLColl.size(); i++){
    if(!PassIDCriteria(MuLColl.at(i), MuTID)){
      float FR=GetFakeRate_Data(MuLColl.at(i),MuTID+"_"+MuLID+"_"+FilterInfo+ConeMethod+SystOpt+SystKind);
      fakeweight*=-FR/(1-FR);
      NLooseNotTight++;
    }
  }
  for(int i=0; i<(int) EleLColl.size(); i++){
    if(!PassIDCriteria(EleLColl.at(i), EleTID)){
      float FR=GetFakeRate_Data(EleLColl.at(i), "QCD_"+EleTID+"_"+EleLID+"_"+ConeMethod+SystOpt+SystKind);
      fakeweight*=-FR/(1-FR);
      NLooseNotTight++;
    }
  }
  if(NLooseNotTight==0) return 0.;

  return fakeweight;

}


float AnalyzerCore::GetFakeWeight_MC(std::vector<snu::KMuon>& MuLColl, std::vector<snu::KElectron>& EleLColl, TString MuLID, TString MuTID, TString EleLID, TString EleTID, TString Option){

  float fakeweight=-1.; int NLooseNotTight=0;

  TString FilterInfo="", ConeMethod="";
  if     (Option.Contains("NoFilter"))  FilterInfo="NoFilter";
  else if(Option.Contains("TrkIsoVVL")) FilterInfo="TrkIsoVVL";
  if     (Option.Contains("ConeE"))     ConeMethod="ConeE";
  else if(Option.Contains("ConeSUSY"))  ConeMethod="ConeSUSY";
  else if(Option.Contains("FOPt"))      ConeMethod="FOPt";


  for(int i=0; i<(int) MuLColl.size(); i++){
    if(!PassIDCriteria(MuLColl.at(i), MuTID)){
      float FR=GetFakeRate_MC(MuLColl.at(i), "QCD_"+MuTID+"_"+MuLID+"_"+FilterInfo+ConeMethod);
      fakeweight*=-FR/(1-FR);
      NLooseNotTight++;
    }
  }
  for(int i=0; i<(int) EleLColl.size(); i++){
    if(!PassIDCriteria(EleLColl.at(i), EleTID)){
      float FR=GetFakeRate_MC(EleLColl.at(i), "QCD_"+EleTID+"_"+EleLID+"_"+ConeMethod);
      fakeweight*=-FR/(1-FR);
      NLooseNotTight++;
    }
  }
  if(NLooseNotTight==0) return 0.;

  return fakeweight;

}


float AnalyzerCore::GetFakeRate_Data(snu::KElectron Ele, TString Option){

  float FR=0., TightIsoCut=0.;
  int SystDir=0.; bool Syst_FR=Option.Contains("Syst");
  bool bool_PrN=false, bool_MES=false, bool_JFlav=false, bool_All=true;
  if(Syst_FR){
    if(Option.Contains("Up")){ SystDir=1; } else if(Option.Contains("Down")){ SystDir=-1; } 
    if     (Option.Contains("PrN")  ){ bool_PrN=true;   bool_All=false;}
    else if(Option.Contains("MES")  ){ bool_MES=true;   bool_All=false;}
    else if(Option.Contains("JFlav")){ bool_JFlav=true; bool_All=false;}
    
  }
  if(Option.Contains("ConeSUSY")){
    if     (Option.Contains("Isop06")) TightIsoCut = 0.06;
    else if(Option.Contains("Isop08")) TightIsoCut = 0.08;
    else if(Option.Contains("Isop1"))  TightIsoCut = 0.1;
  }
  float PTCorr=PTConeCorr(Ele, TightIsoCut);
    if(Option.Contains("FOPt")) PTCorr=Ele.Pt();
  float fEta=fabs(Ele.Eta());

  if(Option.Contains("QCD_POGWP90Isop06IPp025p1sig4_HctoWAFakeLoose_ConeSUSY")){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.10963  ;
      else if(PTCorr<50)  FR=0.072741 ;
      else                FR=0.0709206;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.110035 ;
      else if(PTCorr<50)  FR=0.0752341;
      else                FR=0.0718434;
    }
    else{
      if     (PTCorr<35)  FR=0.144356;
      else if(PTCorr<50)  FR=0.106974;
      else                FR=0.128996;
    }

    if(Syst_FR && SystDir<0){
      if(fEta<0.8){
        if(bool_All){
          if     (PTCorr<35)   FR*=(1.- 0.0677277 );
          else if(PTCorr<50)   FR*=(1.- 0.404283  );
          else                 FR*=(1.- 0.442664  );
        }
        else if(bool_PrN){
          if     (PTCorr<35)   FR*=(1.- 0.0677277 );
          else if(PTCorr<50)   FR*=(1.- 0.172508 );
          else                 FR*=(1.- 0.279723 );
        }
        else if(bool_MES){
          if     (PTCorr<35)   FR*=(1.- 0.000000 );
          else if(PTCorr<50)   FR*=(1.- 0.262032 );
          else                 FR*=(1.- 0.343083 );
        }
        else if(bool_JFlav){
          if     (PTCorr<35)   FR*=(1.- 0.000000 );
          else if(PTCorr<50)   FR*=(1.- 0.255001 );
          else                 FR*=(1.- 0.000000 );
        }
      }
      else if(fEta<1.479){
        if(bool_All){
          if     (PTCorr<35)   FR*=(1.- 0.144574 );
          else if(PTCorr<50)   FR*=(1.- 0.43086  );
          else                 FR*=(1.- 0.224811 );
        }
        else if(bool_PrN){
          if     (PTCorr<35)   FR*=(1.- 0.040827 );
          else if(PTCorr<50)   FR*=(1.- 0.122534 );
          else                 FR*=(1.- 0.198652 );
        }
        else if(bool_MES){
          if     (PTCorr<35)   FR*=(1.- 0.12475 );
          else if(PTCorr<50)   FR*=(1.- 0.300166 );
          else                 FR*=(1.- 0.105249 );
        }
        else if(bool_JFlav){
          if     (PTCorr<35)   FR*=(1.- 0.0605977 );
          else if(PTCorr<50)   FR*=(1.- 0.283771 );
          else                 FR*=(1.- 0.000000 );
        }
      }
      else{
        if(bool_All){
          if     (PTCorr<35)   FR*=(1.- 0.055181  );
          else if(PTCorr<50)   FR*=(1.- 0.0939935 );
          else                 FR*=(1.- 0.135551  );
        }
        else if(bool_PrN){
          if     (PTCorr<35)   FR*=(1.- 0.0132059 );
          else if(PTCorr<50)   FR*=(1.- 0.0401004 );
          else                 FR*=(1.- 0.0465298 );
        }
        else if(bool_MES){
          if     (PTCorr<35)   FR*=(1.- 0.0535775 );
          else if(PTCorr<50)   FR*=(1.- 0.0569569 );
          else                 FR*=(1.- 0.0715485 );
        }
        else if(bool_JFlav){
          if     (PTCorr<35)   FR*=(1.- 0.000000 );
          else if(PTCorr<50)   FR*=(1.- 0.0631082 );
          else                 FR*=(1.- 0.105308 );
        }
      }
    }
    else if(Syst_FR && SystDir>0){
      if(fEta<0.8){
        if(bool_All){
          if     (PTCorr<35)   FR*=(1.+ 0.0956376 );
          else if(PTCorr<50)   FR*=(1.+ 0.169789  );
          else                 FR*=(1.+ 0.287828  );
        }
        else if(bool_PrN){
          if     (PTCorr<35)   FR*=(1.+ 0.0663645 );
          else if(PTCorr<50)   FR*=(1.+ 0.167165 );
          else                 FR*=(1.+ 0.266838 );
        }
        else if(bool_MES){
          if     (PTCorr<35)   FR*=(1.+ 0.0512289 );
          else if(PTCorr<50)   FR*=(1.+ 0.0297342 );
          else                 FR*=(1.+ 0.0000000 );
        }
        else if(bool_JFlav){
          if     (PTCorr<35)   FR*=(1.+ 0.0460207 );
          else if(PTCorr<50)   FR*=(1.+ 0.0000000 );
          else                 FR*=(1.+ 0.1079    );
        }
      }
      else if(fEta<1.479){
        if(bool_All){
          if     (PTCorr<35)   FR*=(1.+ 0.0953948 );
          else if(PTCorr<50)   FR*=(1.+ 0.127346  );
          else                 FR*=(1.+ 0.221431  );
        }
        else if(bool_PrN){
          if     (PTCorr<35)   FR*=(1.+ 0.0403147 );
          else if(PTCorr<50)   FR*=(1.+ 0.119671 );
          else                 FR*=(1.+ 0.191853 );
        }
        else if(bool_MES){
          if     (PTCorr<35)   FR*=(1.+ 0.0864575 );
          else if(PTCorr<50)   FR*=(1.+ 0.0435422 );
          else                 FR*=(1.+ 0.0650656 );
        }
        else if(bool_JFlav){
          if     (PTCorr<35)   FR*=(1.+ 0.0 );
          else if(PTCorr<50)   FR*=(1.+ 0.0 );
          else                 FR*=(1.+ 0.0893912 );
        }
      }
      else{
        if(bool_All){
          if     (PTCorr<35)   FR*=(1.+ 0.136633  );
          else if(PTCorr<50)   FR*=(1.+ 0.0603159 );
          else                 FR*=(1.+ 0.047715  );
        }
        else if(bool_PrN){
          if     (PTCorr<35)   FR*=(1.+ 0.0131328 );
          else if(PTCorr<50)   FR*=(1.+ 0.0396541 );
          else                 FR*=(1.+ 0.0458165 );
        }
        else if(bool_MES){
          if     (PTCorr<35)   FR*=(1.+ 0.0835295 );
          else if(PTCorr<50)   FR*=(1.+ 0.0454484 );
          else                 FR*=(1.+ 0.0133254 );
        }
        else if(bool_JFlav){
          if     (PTCorr<35)   FR*=(1.+ 0.107326 );
          else if(PTCorr<50)   FR*=(1.+ 0.0 );
          else                 FR*=(1.+ 0.0 );
        }
      }
    }
  }
  else{ cout<<"No such Ele FR!"<<endl; }

  return FR;
}

float AnalyzerCore::GetFakeRate_Data(snu::KMuon Mu, TString Option){

  float FR=0., TightIsoCut=0.;
  int SystDir=0.; bool Syst_FR=Option.Contains("Syst");
  bool bool_PrN=false, bool_MES=false, bool_JFlav=false, bool_All=true;
  if(Syst_FR){
    if(Option.Contains("Up")){ SystDir=1; } else if(Option.Contains("Down")){ SystDir=-1; } 
    if     (Option.Contains("PrN")  ){ bool_PrN=true;   bool_All=false;}
    else if(Option.Contains("MES")  ){ bool_MES=true;   bool_All=false;}
    else if(Option.Contains("JFlav")){ bool_JFlav=true; bool_All=false;}
  }
  if(Option.Contains("ConeSUSY")){
    if     (Option.Contains("TIsop15")) TightIsoCut = 0.15;
    else if(Option.Contains("TIsop20")) TightIsoCut = 0.2;
    else if(Option.Contains("TIsop25")) TightIsoCut = 0.25;
  }

  float PTCorr=PTConeCorr(Mu, TightIsoCut);
    if(Option.Contains("FOPt")) PTCorr=Mu.Pt();
  float fEta=fabs(Mu.Eta());

  if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.342129;
      else if(PTCorr<20)  FR=0.161717;
      else if(PTCorr<25)  FR=0.124337;
      else if(PTCorr<35)  FR=0.124102;
      else                FR=0.112075;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.3853  ;
      else if(PTCorr<20)  FR=0.200014;
      else if(PTCorr<25)  FR=0.181963;
      else if(PTCorr<35)  FR=0.143144;
      else                FR=0.149275;
    }
    else{
      if     (PTCorr<15)  FR=0.386585;
      else if(PTCorr<20)  FR=0.228739;
      else if(PTCorr<25)  FR=0.203977;
      else if(PTCorr<35)  FR=0.190597;
      else                FR=0.180573;
    }

    if(Syst_FR && SystDir<0){
      if(fEta<0.9){
        if(bool_All){
          if     (PTCorr<15)   FR*=(1.- 0.00570765 );
          else if(PTCorr<20)   FR*=(1.- 0.146032   );
          else if(PTCorr<25)   FR*=(1.- 0.079018   );
          else if(PTCorr<35)   FR*=(1.- 0.211878   );
          else                 FR*=(1.- 0.364666   );
        }
        else if(bool_PrN){
          if     (PTCorr<15)   FR*=(1.- 0.000860614 );
          else if(PTCorr<20)   FR*=(1.- 0.00358884 );
          else if(PTCorr<25)   FR*=(1.- 0.00991981 );
          else if(PTCorr<35)   FR*=(1.- 0.0237898 );
          else                 FR*=(1.- 0.0911435 );
        }
        else if(bool_MES){
          if     (PTCorr<15)   FR*=(1.- 0.0011784 );
          else if(PTCorr<20)   FR*=(1.- 0.123918 );
          else if(PTCorr<25)   FR*=(1.- 0.0 );
          else if(PTCorr<35)   FR*=(1.- 0.0196743 );
          else                 FR*=(1.- 0.0976778 );
        }
        else if(bool_JFlav){
          if     (PTCorr<15)   FR*=(1.- 0.00551797 );
          else if(PTCorr<20)   FR*=(1.- 0.0771795 );
          else if(PTCorr<25)   FR*=(1.- 0.0783928 );
          else if(PTCorr<35)   FR*=(1.- 0.209617 );
          else                 FR*=(1.- 0.339313 );
        }
      }
      else if(fEta<1.6){
        if(bool_All){
          if     (PTCorr<15)   FR*=(1.- 0.0378566 );
          else if(PTCorr<20)   FR*=(1.- 0.173204  );
          else if(PTCorr<25)   FR*=(1.- 0.13663   );
          else if(PTCorr<35)   FR*=(1.- 0.0838593 );
          else                 FR*=(1.- 0.27241   );
        }
        else if(bool_PrN){
          if     (PTCorr<15)   FR*=(1.- 0.000628793 );
          else if(PTCorr<20)   FR*=(1.- 0.00280486 );
          else if(PTCorr<25)   FR*=(1.- 0.00629162 );
          else if(PTCorr<35)   FR*=(1.- 0.0202333 );
          else                 FR*=(1.- 0.0691845 );
        }
        else if(bool_MES){
          if     (PTCorr<15)   FR*=(1.- 0.0378514 );
          else if(PTCorr<20)   FR*=(1.- 0.093347 );
          else if(PTCorr<25)   FR*=(1.- 0.0988228 );
          else if(PTCorr<35)   FR*=(1.- 0.0 );
          else                 FR*=(1.- 0.0863021 );
        }
        else if(bool_JFlav){
          if     (PTCorr<15)   FR*=(1.- 0.0 );
          else if(PTCorr<20)   FR*=(1.- 0.14587 );
          else if(PTCorr<25)   FR*=(1.- 0.0941392 );
          else if(PTCorr<35)   FR*=(1.- 0.0813817 );
          else                 FR*=(1.- 0.248943 );
        }
      }
      else{
        if(bool_All){
          if     (PTCorr<15)   FR*=(1.- 0.0761964 );
          else if(PTCorr<20)   FR*=(1.- 0.137059  );
          else if(PTCorr<25)   FR*=(1.- 0.246004  );
          else if(PTCorr<35)   FR*=(1.- 0.0928262 );
          else                 FR*=(1.- 0.214074  );
        }
        else if(bool_PrN){
          if     (PTCorr<15)   FR*=(1.- 0.00059424 );
          else if(PTCorr<20)   FR*=(1.- 0.00233966 );
          else if(PTCorr<25)   FR*=(1.- 0.00553085 );
          else if(PTCorr<35)   FR*=(1.- 0.0146442 );
          else                 FR*=(1.- 0.0657493 );
        }
        else if(bool_MES){
          if     (PTCorr<15)   FR*=(1.- 0.0384606 );
          else if(PTCorr<20)   FR*=(1.- 0.0123603 );
          else if(PTCorr<25)   FR*=(1.- 0.119505 );
          else if(PTCorr<35)   FR*=(1.- 0.071105 );
          else                 FR*=(1.- 0.0360618 );
        }
        else if(bool_JFlav){
          if     (PTCorr<15)   FR*=(1.- 0.0657747 );
          else if(PTCorr<20)   FR*=(1.- 0.136481 );
          else if(PTCorr<25)   FR*=(1.- 0.214955 );
          else if(PTCorr<35)   FR*=(1.- 0.0578475 );
          else                 FR*=(1.- 0.20051 );
        }
      }
    }
    else if(Syst_FR && SystDir>0){
      if(fEta<0.9){
        if(bool_All){
          if     (PTCorr<15)   FR*=(1.+ 0.0467043  );
          else if(PTCorr<20)   FR*=(1.+ 0.00598351 );
          else if(PTCorr<25)   FR*=(1.+ 0.0564788  );
          else if(PTCorr<35)   FR*=(1.+ 0.0280901  );
          else                 FR*=(1.+ 0.0889982  );
        }
        else if(bool_PrN){
          if     (PTCorr<15)   FR*=(1.+ 0.000859589 );
          else if(PTCorr<20)   FR*=(1.+ 0.00358265 );
          else if(PTCorr<25)   FR*=(1.+ 0.00988772 );
          else if(PTCorr<35)   FR*=(1.+ 0.023615 );
          else                 FR*=(1.+ 0.0889982 );
        }
        else if(bool_MES){
          if     (PTCorr<15)   FR*=(1.+ 0.0466964 );
          else if(PTCorr<20)   FR*=(1.+ 0.00479239 );
          else if(PTCorr<25)   FR*=(1.+ 0.0556066 );
          else if(PTCorr<35)   FR*=(1.+ 0.0152114 );
          else                 FR*=(1.+ 0.0 );
        }
        else if(bool_JFlav){
          if     (PTCorr<15)   FR*=(1.+ 0.0 );
          else if(PTCorr<20)   FR*=(1.+ 0.0 );
          else if(PTCorr<25)   FR*=(1.+ 0.0 );
          else if(PTCorr<35)   FR*=(1.+ 0.0 );
          else                 FR*=(1.+ 0.0 );
        }
      }
      else if(fEta<1.6){
        if(bool_All){
          if     (PTCorr<15)   FR*=(1.+ 0.0339258 );
          else if(PTCorr<20)   FR*=(1.+ 0.0216653 );
          else if(PTCorr<25)   FR*=(1.+ 0.0344182 );
          else if(PTCorr<35)   FR*=(1.+ 0.072051  );
          else                 FR*=(1.+ 0.0674714 );
        }
        else if(bool_PrN){
          if     (PTCorr<15)   FR*=(1.+ 0.000628106 );
          else if(PTCorr<20)   FR*=(1.+ 0.00280003 );
          else if(PTCorr<25)   FR*=(1.+ 0.00627158 );
          else if(PTCorr<35)   FR*=(1.+ 0.0200861 );
          else                 FR*=(1.+ 0.0674714 );
        }
        else if(bool_MES){
          if     (PTCorr<15)   FR*=(1.+ 0.0338968 );
          else if(PTCorr<20)   FR*=(1.+ 0.0214836 );
          else if(PTCorr<25)   FR*=(1.+ 0.0338419 );
          else if(PTCorr<35)   FR*=(1.+ 0.0691946 );
          else                 FR*=(1.+ 0. );
        }
        else if(bool_JFlav){
          if     (PTCorr<15)   FR*=(1.+ 0.00125572 );
          else if(PTCorr<20)   FR*=(1.+ 0. );
          else if(PTCorr<25)   FR*=(1.+ 0. );
          else if(PTCorr<35)   FR*=(1.+ 0. );
          else                 FR*=(1.+ 0. );
        }
      }
      else{
        if(bool_All){
          if     (PTCorr<15)   FR*=(1.+ 0.043051 );
          else if(PTCorr<20)   FR*=(1.+ 0.0239283);
          else if(PTCorr<25)   FR*=(1.+ 0.108454 );
          else if(PTCorr<35)   FR*=(1.+ 0.0785984);
          else                 FR*=(1.+ 0.0663562);
        }
        else if(bool_PrN){
          if     (PTCorr<15)   FR*=(1.+ 0.000593604 );
          else if(PTCorr<20)   FR*=(1.+ 0.00233561 );
          else if(PTCorr<25)   FR*=(1.+ 0.00551262 );
          else if(PTCorr<35)   FR*=(1.+ 0.0145331 );
          else                 FR*=(1.+ 0.0637743 );
        }
        else if(bool_MES){
          if     (PTCorr<15)   FR*=(1.+ 0.0430469 );
          else if(PTCorr<20)   FR*=(1.+ 0.023814 );
          else if(PTCorr<25)   FR*=(1.+ 0.108314 );
          else if(PTCorr<35)   FR*=(1.+ 0.0772431 );
          else                 FR*=(1.+ 0.0183298 );
        }
        else if(bool_JFlav){
          if     (PTCorr<15)   FR*=(1.+ 0. );
          else if(PTCorr<20)   FR*=(1.+ 0. );
          else if(PTCorr<25)   FR*=(1.+ 0. );
          else if(PTCorr<35)   FR*=(1.+ 0. );
          else                 FR*=(1.+ 0. );
        }
      }
    }
  }
  else if(Option.Contains("POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY")){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.2771  ;
      else if(PTCorr<20)  FR=0.193227;
      else if(PTCorr<25)  FR=0.162945;
      else if(PTCorr<35)  FR=0.165255;
      else                FR=0.158028;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.310046;
      else if(PTCorr<20)  FR=0.221337;
      else if(PTCorr<25)  FR=0.214164;
      else if(PTCorr<35)  FR=0.179193;
      else                FR=0.184973;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.340931;
      else if(PTCorr<20)  FR=0.263206;
      else if(PTCorr<25)  FR=0.245025;
      else if(PTCorr<35)  FR=0.222435;
      else                FR=0.214271;
    }
    else{
      if     (PTCorr<15)  FR=0.338846;
      else if(PTCorr<20)  FR=0.269771;
      else if(PTCorr<25)  FR=0.255547;
      else if(PTCorr<35)  FR=0.25313 ;
      else                FR=0.250546;
    }

    if(Syst_FR && SystDir<0){
      if(fEta<0.9){
        if     (PTCorr<15)   FR*=(1.- 0.0655089 );
        else if(PTCorr<20)   FR*=(1.- 0.304225  );
        else if(PTCorr<25)   FR*=(1.- 0.301855  );
        else if(PTCorr<35)   FR*=(1.- 0.387517  );
        else                 FR*=(1.- 0.483798  );
      }
      else if(fEta<1.6){
        if     (PTCorr<15)   FR*=(1.- 0.121366 );
        else if(PTCorr<20)   FR*=(1.- 0.365817 );
        else if(PTCorr<25)   FR*=(1.- 0.310581 );
        else if(PTCorr<35)   FR*=(1.- 0.27695  );
        else                 FR*=(1.- 0.438805 );
      }
      else if(fEta<2.1){
        if     (PTCorr<15)   FR*=(1.- 0.252924 );
        else if(PTCorr<20)   FR*=(1.- 0.368097 );
        else if(PTCorr<25)   FR*=(1.- 0.460861 );
        else if(PTCorr<35)   FR*=(1.- 0.242128 );
        else                 FR*=(1.- 0.415018 );
      }
      else{
        if     (PTCorr<15)   FR*=(1.- 0.243664 );
        else if(PTCorr<20)   FR*=(1.- 0.341107 );
        else if(PTCorr<25)   FR*=(1.- 0.442997 );
        else if(PTCorr<35)   FR*=(1.- 0.292832 );
        else                 FR*=(1.- 0.308315 );
      }
    }
    else if(Syst_FR && SystDir>0){
      if(fEta<0.9){
        if     (PTCorr<15)   FR*=(1.+ 0.0681634 );
        else if(PTCorr<20)   FR*=(1.+ 0.0553126 );
        else if(PTCorr<25)   FR*=(1.+ 0.101094  );
        else if(PTCorr<35)   FR*=(1.+ 0.07793   );
        else                 FR*=(1.+ 0.199967  );
      }
      else if(fEta<1.6){
        if     (PTCorr<15)   FR*=(1.+ 0.0625109);
        else if(PTCorr<20)   FR*=(1.+ 0.0856362);
        else if(PTCorr<25)   FR*=(1.+ 0.0966223);
        else if(PTCorr<35)   FR*=(1.+ 0.108313 );
        else                 FR*=(1.+ 0.146165 );
      }
      else if(fEta<2.1){
        if     (PTCorr<15)   FR*=(1.+ 0.0774913);
        else if(PTCorr<20)   FR*=(1.+ 0.0759516);
        else if(PTCorr<25)   FR*=(1.+ 0.147378 );
        else if(PTCorr<35)   FR*=(1.+ 0.131883 );
        else                 FR*=(1.+ 0.14489  );
      }
      else{
        if     (PTCorr<15)   FR*=(1.+ 0.0810315);
        else if(PTCorr<20)   FR*=(1.+ 0.0809193);
        else if(PTCorr<25)   FR*=(1.+ 0.0970186);
        else if(PTCorr<35)   FR*=(1.+ 0.110969 );
        else                 FR*=(1.+ 0.149472 );
      }
    }
  }
  else {cout<<"No Such Mu FR!"<<endl;}

  return FR;
}


float AnalyzerCore::GetFakeRate_MC(snu::KMuon Mu, TString Option){

  float FR=0., TightIsoCut=0.;
  if(Option.Contains("ConeSUSY")){
    if     (Option.Contains("TIsop15")) TightIsoCut = 0.15;
    else if(Option.Contains("TIsop20")) TightIsoCut = 0.2;
    else if(Option.Contains("TIsop25")) TightIsoCut = 0.25;
  }
  float PTCorr=PTConeCorr(Mu, TightIsoCut);
    if(Option.Contains("FOPt")) PTCorr=Mu.Pt();
  float fEta=fabs(Mu.Eta());

  if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.395195 ;
      else if(PTCorr<20)  FR=0.202967 ;
      else if(PTCorr<25)  FR=0.169647 ;
      else if(PTCorr<35)  FR=0.143029 ;
      else if(PTCorr<50)  FR=0.111431 ;
      else                FR=0.0951533;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.418279;
      else if(PTCorr<20)  FR=0.238722;
      else if(PTCorr<25)  FR=0.217297;
      else if(PTCorr<35)  FR=0.177343;
      else if(PTCorr<50)  FR=0.161958;
      else                FR=0.136925;
    }
    else{
      if     (PTCorr<15)  FR=0.435522;
      else if(PTCorr<20)  FR=0.276537;
      else if(PTCorr<25)  FR=0.259684;
      else if(PTCorr<35)  FR=0.245707;
      else if(PTCorr<50)  FR=0.207092;
      else                FR=0.187308;
    }
  }
  else if(Option=="QCD_AwayB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.396487;
      else if(PTCorr<20)  FR=0.207776;
      else if(PTCorr<25)  FR=0.183522;
      else if(PTCorr<35)  FR=0.155841;
      else if(PTCorr<50)  FR=0.123879;
      else                FR=0.107774;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.421205;
      else if(PTCorr<20)  FR=0.245923;
      else if(PTCorr<25)  FR=0.230051;
      else if(PTCorr<35)  FR=0.192691;
      else if(PTCorr<50)  FR=0.170713;
      else                FR=0.150102;
    }
    else{
      if     (PTCorr<15)  FR=0.440847;
      else if(PTCorr<20)  FR=0.289948;
      else if(PTCorr<25)  FR=0.272952;
      else if(PTCorr<35)  FR=0.261306;
      else if(PTCorr<50)  FR=0.224457;
      else                FR=0.20334 ;
    }
  }
  else if(Option=="QCD_HasB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.300225 ;
      else if(PTCorr<20)  FR=0.143649 ;
      else if(PTCorr<25)  FR=0.117231 ;
      else if(PTCorr<35)  FR=0.103437 ;
      else if(PTCorr<50)  FR=0.0767057;
      else                FR=0.0631149;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.304083 ;
      else if(PTCorr<20)  FR=0.176359 ;
      else if(PTCorr<25)  FR=0.164183 ;
      else if(PTCorr<35)  FR=0.12428  ;
      else if(PTCorr<50)  FR=0.13231  ;
      else                FR=0.0932271;
    }
    else{
      if     (PTCorr<15)  FR=0.239239;
      else if(PTCorr<20)  FR=0.160474;
      else if(PTCorr<25)  FR=0.192764;
      else if(PTCorr<35)  FR=0.176604;
      else if(PTCorr<50)  FR=0.137952;
      else                FR=0.131124;
    }
  }
  else if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.376251 ;
      else if(PTCorr<20)  FR=0.180005 ;
      else if(PTCorr<25)  FR=0.148464 ;
      else if(PTCorr<35)  FR=0.122854 ;
      else if(PTCorr<50)  FR=0.0952236;
      else                FR=0.081313 ;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.394537;
      else if(PTCorr<20)  FR=0.206124;
      else if(PTCorr<25)  FR=0.18436 ;
      else if(PTCorr<35)  FR=0.151379;
      else if(PTCorr<50)  FR=0.135037;
      else                FR=0.117474;
    }
    else{
      if     (PTCorr<15)  FR=0.409673;
      else if(PTCorr<20)  FR=0.23768 ;
      else if(PTCorr<25)  FR=0.21949 ;
      else if(PTCorr<35)  FR=0.206409;
      else if(PTCorr<50)  FR=0.172383;
      else                FR=0.158458;
    }
  }
  else if(Option=="QCD_AwayB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.378176 ;
      else if(PTCorr<20)  FR=0.184572 ;
      else if(PTCorr<25)  FR=0.160799 ;
      else if(PTCorr<35)  FR=0.135949 ;
      else if(PTCorr<50)  FR=0.107135 ;
      else                FR=0.0933122;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.398919;
      else if(PTCorr<20)  FR=0.213567;
      else if(PTCorr<25)  FR=0.196568;
      else if(PTCorr<35)  FR=0.164957;
      else if(PTCorr<50)  FR=0.145089;
      else                FR=0.130676;
    }
    else{
      if     (PTCorr<15)  FR=0.415817;
      else if(PTCorr<20)  FR=0.250131;
      else if(PTCorr<25)  FR=0.231834;
      else if(PTCorr<35)  FR=0.221139;
      else if(PTCorr<50)  FR=0.187845;
      else                FR=0.173545;
    }
  }
  else if(Option=="QCD_HasB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.28293  ;
      else if(PTCorr<20)  FR=0.13529  ;
      else if(PTCorr<25)  FR=0.106368 ;
      else if(PTCorr<35)  FR=0.0874014;
      else if(PTCorr<50)  FR=0.0660323;
      else                FR=0.0541946;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.276023 ;
      else if(PTCorr<20)  FR=0.154409 ;
      else if(PTCorr<25)  FR=0.140863 ;
      else if(PTCorr<35)  FR=0.110194 ;
      else if(PTCorr<50)  FR=0.105941 ;
      else                FR=0.0796622;
    }
    else{
      if     (PTCorr<15)  FR=0.255897;
      else if(PTCorr<20)  FR=0.14699 ;
      else if(PTCorr<25)  FR=0.163857;
      else if(PTCorr<35)  FR=0.149301;
      else if(PTCorr<50)  FR=0.117876;
      else                FR=0.111363;
    }
  }//For Here Non Main ID Studies
  else if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.31422 ;
      else if(PTCorr<20)  FR=0.2225  ;
      else if(PTCorr<25)  FR=0.201252;
      else if(PTCorr<35)  FR=0.174787;
      else if(PTCorr<50)  FR=0.144774;
      else                FR=0.137056;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.328501;
      else if(PTCorr<20)  FR=0.250262;
      else if(PTCorr<25)  FR=0.239327;
      else if(PTCorr<35)  FR=0.204919;
      else if(PTCorr<50)  FR=0.200699;
      else                FR=0.188323;
    }
    else{
      if     (PTCorr<15)  FR=0.369477;
      else if(PTCorr<20)  FR=0.29852 ;
      else if(PTCorr<25)  FR=0.28201 ;
      else if(PTCorr<35)  FR=0.271146;
      else if(PTCorr<50)  FR=0.238203;
      else                FR=0.234892;
    }
  }
  else if(Option=="QCD_AwayB_POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.319891;
      else if(PTCorr<20)  FR=0.246482;
      else if(PTCorr<25)  FR=0.249877;
      else if(PTCorr<35)  FR=0.223699;
      else if(PTCorr<50)  FR=0.195308;
      else                FR=0.18359 ;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.338585;
      else if(PTCorr<20)  FR=0.286107;
      else if(PTCorr<25)  FR=0.297202;
      else if(PTCorr<35)  FR=0.262352;
      else if(PTCorr<50)  FR=0.257737;
      else                FR=0.252453;
    }
    else{
      if     (PTCorr<15)  FR=0.385327;
      else if(PTCorr<20)  FR=0.345415;
      else if(PTCorr<25)  FR=0.336095;
      else if(PTCorr<35)  FR=0.333724;
      else if(PTCorr<50)  FR=0.298049;
      else                FR=0.288971;
    }
  }
  else if(Option=="QCD_HasB_POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){//NearB
      if     (PTCorr<15)  FR=0.082816 ;
      else if(PTCorr<20)  FR=0.0835979;
      else if(PTCorr<25)  FR=0.096017 ;
      else if(PTCorr<35)  FR=0.0864538;
      else if(PTCorr<50)  FR=0.0610422;
      else                FR=0.0662293;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.0954808;
      else if(PTCorr<20)  FR=0.100042 ;
      else if(PTCorr<25)  FR=0.108263 ;
      else if(PTCorr<35)  FR=0.0848722;
      else if(PTCorr<50)  FR=0.120748 ;
      else                FR=0.0643254;
    }
    else{
      if     (PTCorr<15)  FR=0.0803636;
      else if(PTCorr<20)  FR=0.0874016;
      else if(PTCorr<25)  FR=0.113497 ;
      else if(PTCorr<35)  FR=0.114078 ;
      else if(PTCorr<50)  FR=0.103653 ;
      else                FR=0.167258 ;
    }
  }
  else if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.420458;
      else if(PTCorr<20)  FR=0.210949;
      else if(PTCorr<25)  FR=0.174416;
      else if(PTCorr<35)  FR=0.146076;
      else if(PTCorr<50)  FR=0.113188;
      else                FR=0.097384;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.462609;
      else if(PTCorr<20)  FR=0.256771;
      else if(PTCorr<25)  FR=0.228671;
      else if(PTCorr<35)  FR=0.185761;
      else if(PTCorr<50)  FR=0.167516;
      else                FR=0.142433;
    }
    else{
      if     (PTCorr<15)  FR=0.542713;
      else if(PTCorr<20)  FR=0.330838;
      else if(PTCorr<25)  FR=0.300315;
      else if(PTCorr<35)  FR=0.282436;
      else if(PTCorr<50)  FR=0.234402;
      else                FR=0.210434;
    }
  }
  else if(Option=="QCD_AwayB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.421107;
      else if(PTCorr<20)  FR=0.214535;
      else if(PTCorr<25)  FR=0.186623;
      else if(PTCorr<35)  FR=0.157792;
      else if(PTCorr<50)  FR=0.124924;
      else                FR=0.110091;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.463734;
      else if(PTCorr<20)  FR=0.260781;
      else if(PTCorr<25)  FR=0.236924;
      else if(PTCorr<35)  FR=0.198074;
      else if(PTCorr<50)  FR=0.173928;
      else                FR=0.15409 ;
    }
    else{
      if     (PTCorr<15)  FR=0.543834;
      else if(PTCorr<20)  FR=0.335507;
      else if(PTCorr<25)  FR=0.303724;
      else if(PTCorr<35)  FR=0.288788;
      else if(PTCorr<50)  FR=0.245314;
      else                FR=0.22206 ;
    }
  }
  else if(Option=="QCD_HasB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){//NearB
      if     (PTCorr<15)  FR=0.297745 ;
      else if(PTCorr<20)  FR=0.15785  ;
      else if(PTCorr<25)  FR=0.123827 ;
      else if(PTCorr<35)  FR=0.100511 ;
      else if(PTCorr<50)  FR=0.06543  ;
      else                FR=0.0693239;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.324575 ;
      else if(PTCorr<20)  FR=0.207702 ;
      else if(PTCorr<25)  FR=0.176821 ;
      else if(PTCorr<35)  FR=0.119094 ;
      else if(PTCorr<50)  FR=0.160515 ;
      else                FR=0.0776659;
    }
    else{
      if     (PTCorr<15)  FR=0.397288;
      else if(PTCorr<20)  FR=0.221381;
      else if(PTCorr<25)  FR=0.232421;
      else if(PTCorr<35)  FR=0.208876;
      else if(PTCorr<50)  FR=0.167898;
      else                FR=0.214372;
    }
  }
  else if(Option=="QCD_BjMatch_POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.212887;
      else if(PTCorr<20)  FR=0.170592;
      else if(PTCorr<25)  FR=0.16553 ;
      else if(PTCorr<35)  FR=0.153397;
      else if(PTCorr<50)  FR=0.123731;
      else                FR=0.104981;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.239631;
      else if(PTCorr<20)  FR=0.206327;
      else if(PTCorr<25)  FR=0.216909;
      else if(PTCorr<35)  FR=0.182379;
      else if(PTCorr<50)  FR=0.170383;
      else                FR=0.156924;
    }
    else{
      if     (PTCorr<15)  FR=0.22021 ;
      else if(PTCorr<20)  FR=0.228938;
      else if(PTCorr<25)  FR=0.24682 ;
      else if(PTCorr<35)  FR=0.241839;
      else if(PTCorr<50)  FR=0.204816;
      else                FR=0.19224 ;
    }
  }
  else if(Option=="QCD_CjMatch_POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.280125;
      else if(PTCorr<20)  FR=0.206693;
      else if(PTCorr<25)  FR=0.237403;
      else if(PTCorr<35)  FR=0.232245;
      else if(PTCorr<50)  FR=0.222692;
      else                FR=0.233209;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.285224;
      else if(PTCorr<20)  FR=0.263604;
      else if(PTCorr<25)  FR=0.263886;
      else if(PTCorr<35)  FR=0.261374;
      else if(PTCorr<50)  FR=0.276096;
      else                FR=0.250051;
    }
    else{
      if     (PTCorr<15)  FR=0.316214;
      else if(PTCorr<20)  FR=0.275994;
      else if(PTCorr<25)  FR=0.315642;
      else if(PTCorr<35)  FR=0.333266;
      else if(PTCorr<50)  FR=0.308791;
      else                FR=0.331183;
    }
  }
  else if(Option=="QCD_LjMatch_POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.235959;
      else if(PTCorr<20)  FR=0.163328;
      else if(PTCorr<25)  FR=0.255066;
      else if(PTCorr<35)  FR=0.220489;
      else if(PTCorr<50)  FR=0.150643;
      else                FR=0.217929;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.290873;
      else if(PTCorr<20)  FR=0.200495;
      else if(PTCorr<25)  FR=0.300621;
      else if(PTCorr<35)  FR=0.262121;
      else if(PTCorr<50)  FR=0.316673;
      else                FR=0.371495;
    }
    else{
      if     (PTCorr<15)  FR=0.389803;
      else if(PTCorr<20)  FR=0.364795;
      else if(PTCorr<25)  FR=0.398313;
      else if(PTCorr<35)  FR=0.382796;
      else if(PTCorr<50)  FR=0.383164;
      else                FR=0.421589;
    }
  }
  else { cout<<"RateMC Error"<<endl; }

  return FR;
}


float AnalyzerCore::GetFakeRate_MC(snu::KElectron Ele, TString Option){

  float FR=0., TightIsoCut=0.;
  if(Option.Contains("ConeSUSY")){
    if     (Option.Contains("Isop06")) TightIsoCut = 0.06;
    else if(Option.Contains("Isop08")) TightIsoCut = 0.08;
    else if(Option.Contains("Isop1"))  TightIsoCut = 0.1;
  }
  float PTCorr=PTConeCorr(Ele, TightIsoCut);
    if(Option.Contains("FOPt")) PTCorr=Ele.Pt();
  float fEta=fabs(Ele.Eta());

  if(Option=="QCD_POGWP90Isop06IPp025p1sig4_HctoWAFakeLoose_ConeSUSY"){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.0927136;
      else if(PTCorr<50)  FR=0.0720359;
      else if(PTCorr<70)  FR=0.0711901;
      else                FR=0.0624766;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.149674 ;
      else if(PTCorr<50)  FR=0.0613648;
      else if(PTCorr<70)  FR=0.0846761;
      else                FR=0.0534003;
    }
    else{
      if     (PTCorr<35)  FR=0.159323;
      else if(PTCorr<50)  FR=0.104854;
      else if(PTCorr<70)  FR=0.109523;
      else                FR=0.120316;
    }
  }
  else if(Option=="QCD_POGWP90Isop06IPp025p05sig4_LMVA06Isop4IPp025p05sig4_ConeSUSY"){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.0933608;
      else if(PTCorr<50)  FR=0.0721879;
      else if(PTCorr<70)  FR=0.0714611;
      else                FR=0.0627619;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.156215 ;
      else if(PTCorr<50)  FR=0.0664036;
      else if(PTCorr<70)  FR=0.0919452;
      else                FR=0.0562619;
    }
    else{
      if     (PTCorr<35)  FR=0.161256;
      else if(PTCorr<50)  FR=0.109329;
      else if(PTCorr<70)  FR=0.114383;
      else                FR=0.126938;
    }
  }

  return FR;
}


float AnalyzerCore::PTConeCorr(snu::KElectron Ele, float TightIsoCut){

  float PTCorr=Ele.Pt()*(1.+max(0.,Ele.PFRelIso(0.3)-TightIsoCut));
  return PTCorr;

}


float AnalyzerCore::PTConeCorr(snu::KMuon Mu, float TightIsoCut){

  float PTCorr=Mu.Pt()*(1.+max((float) 0.,RochIso(Mu,"0.4")-TightIsoCut));
  return PTCorr;

}
