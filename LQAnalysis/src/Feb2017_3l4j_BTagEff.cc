// $Id: Feb2017_3l4j_BTagEff.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQFeb2017_3l4j_BTagEff Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Feb2017_3l4j_BTagEff.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Feb2017_3l4j_BTagEff);

 Feb2017_3l4j_BTagEff::Feb2017_3l4j_BTagEff() : AnalyzerCore(), out_muons(0) {

   SetLogName("Feb2017_3l4j_BTagEff");
   Message("In Feb2017_3l4j_BTagEff constructor", INFO);
   InitialiseAnalysis();
 }


 void Feb2017_3l4j_BTagEff::InitialiseAnalysis() throw( LQError ) {
   
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

void Feb2017_3l4j_BTagEff::ExecuteEvents()throw( LQError ){

////Basic Infos///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Apply the gen weight 
   float weight_nopu=1; 
   if(!isData) weight*=MCweight;
   if(!isData) weight_nopu*=MCweight;
   FillHist("GenWeight", MCweight, 1, -10, 10, 1000);

   //Total Event  
   FillCutFlow("NoCut", weight);

   /// Acts on data to remove bad reconstructed event 
   //if(isData&& (! eventbase->GetEvent().LumiMask(lumimask))) return;

   m_logger<<DEBUG<<"RunNumber/Event Number = "<<eventbase->GetEvent().RunNumber()<<" : "<<eventbase->GetEvent().EventNumber()<<LQLogger::endmsg;
   m_logger<<DEBUG<<"isData = "<<isData<<LQLogger::endmsg;


   //Numbet of Vertex NoPUreweight
   FillHist("Nvtx_nocut_NoRW", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);

   //Pileup Reweight
   float pileup_reweight=(1.0);
   if(!k_isdata) { 
     //pileup_reweight=eventbase->GetEvent().PileUpWeight(snu::KEvent::down);
     pileup_reweight=eventbase->GetEvent().PileUpWeight();
     //pileup_reweight = TempPileupWeight();
     //weight*=pileup_reweight;
   } 
   //Numbet of Vertex PUreweight
   FillHist("Nvtx_nocut_PURW", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);


   bool EMu_analysis=false, SingleMu_analysis=false, DoubleMu_analysis=false, DoubleEle_analysis=true;


   //Trigers
   std::vector<TString> triggerslist, triggerlist1, triggerlist2;//  ListTriggersAvailable();

   TString analysis_trigger;
   if(EMu_analysis) {analysis_trigger="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";}
   else if(DoubleMu_analysis) analysis_trigger="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";
   else if(SingleMu_analysis) analysis_trigger="HLT_IsoMu24_v";
   else if(DoubleEle_analysis) analysis_trigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
   //if(!PassTrigger(analysis_trigger)) return; //Not for now..

   //PassTrigger Uses BeginWith so don't stick to exact name

   float trigger_ps_weight=1, weight_trigger_sf=1;
/*   std::vector<snu::KMuon> trigmuColl; std::vector<snu::KElectron> trigeColl;
     eventbase->GetMuonSel()->SetPt(20.);eventbase->GetMuonSel()->SetEta(2.4);eventbase->GetMuonSel()->SetRelIso(0.15);eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);  eventbase->GetMuonSel()->Selection(trigmuColl);
     eventbase->GetElectronSel()->SetPt(20.);eventbase->GetElectronSel()->SetEta(2.4);eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MEDIUM); eventbase->GetElectronSel()->Selection(trigeColl);*/

   trigger_ps_weight=WeightByTrigger(analysis_trigger, TargetLumi);
   //weight_trigger_sf=1.;//TriggerScaleFactor(trigeColl, trigmuColl, analysis_trigger);



   FillHist("TriggerSFWeight" , weight_trigger_sf, 1., 0. , 2., 200); FillHist("TriggerPSWeight", trigger_ps_weight, 1., 0., 500., 500);
   weight*=weight_trigger_sf*trigger_ps_weight;
   weight_nopu*=weight_trigger_sf*trigger_ps_weight;
   //FillCutFlow("TriggerCut", weight); //;in trigger Cut mode
   //FillHist("Basic_prescale", prescale, 1., 0., 2000., 2000);

   //Initial Event Cut : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
   //if(!PassMETFilter()) return;  FillCutFlow("EventCut", weight);

   //Vertex Cut
   //if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; FillCutFlow("VertexCut", weight);
   //Good Primary vtx def:(vtx.ndof()>4&&maxAbsZ<=0)||std::abs(vtx.z())<= 24)&&((maxd0 <=0)||std::abs(vtx.position().rho())<=2)&&!(vtx.isFake()))  



///////Objects in Analysis/////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

   //Primary Object Collection
   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
 
   //SingleMuAnalysis 
   //  eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
   //  eventbase->GetMuonSel()->SetPt(5.);                    eventbase->GetMuonSel()->SetEta(2.4);
   //  eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.25);
   //std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.25);
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl);
   CorrectMuonMomentum(muonLooseColl);

     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.15);
   std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->Selection(muonColl);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
   //std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->SelectMuons(muonColl, BaseSelection::MUON_POG_TIGHT, 10., 2.4);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
   CorrectMuonMomentum(muonColl);
   
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0994, 0.107);//ISO Ichep16 Loose WP
   //eventbase->GetElectronSel()->SetdxyBEMax(0.005, 0.01);  eventbase->GetElectronSel()->SetdzBEMax(0.005, 0.01);
//cout << eventbase->GetElectronSel()->pt_cut_min<<endl;
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);//pt>10/eta<2.4(hole region ~1.4 excluded by default)/EGAMMA_LOOSE/NoExtIso*/
//cout << eventbase->GetElectronSel()->pt_cut_min<<endl;
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//Iso Ichep16 Loose WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.2);
   //std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->Selection(electronColl);//pt>10/eta<2.4(hole region ~1.4 excluded by default)/EGAMMA_LOOSE/NoExtIso*/
   std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->SelectElectrons(electronColl, "ELECTRON_POG_TIGHT", 10., 2.5);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
//cout<<electronLooseColl.size()<<" "<<electronColl.size()<<endl;
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_TIGHT);
     eventbase->GetJetSel()->SetPt(20.);                     eventbase->GetJetSel()->SetEta(2.4);
//   eventbase->GetJetSel()->SetUseJetPileUp(true);
     bool LeptonVeto=true;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonColl, electronColl, isData, true);
   std::vector<snu::KJet> jetLooseColl; eventbase->GetJetSel()->SelectJets(isData, jetLooseColl, muonColl, electronColl, "PFJET_LOOSE", 20., 2.4, true);


//   std::vector<int> bIdxColl=GetSFBJetIdx(jetColl,"Medium");
//   std::vector<int> ljIdxColl=GetSFLJetIdx(jetColl, bIdxColl, "Medium");

//   std::vector<snu::KJet> bjetColl; for(int i=0; i<bIdxColl.size(); i++){bjetColl.push_back(jetColl.at(bIdxColl.at(i)));}
//   std::vector<snu::KJet> ljetColl; for(int i=0; i<ljIdxColl.size(); i++){ljetColl.push_back(jetColl.at(ljIdxColl.at(i)));}

   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");
   std::vector<snu::KJet> ljetColl = SelLightJets(jetColl, "Medium");

   //Temporarily Placed Cuts(Activate ones at original place if trig simul. fully available
   //METFilter
   if(!PassMETFilter()) return;  FillCutFlow("EventCut", weight);

   //Vertex Cut
   if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; FillCutFlow("VertexCut", weight);

   //Trigger eff. emul. and cuts
   float trigger_eff=1.;
   if(isData){
     //if(!PassTrigger(analysis_trigger)) return;
   }
   else{
     if(EMu_analysis) trigger_eff = TriggerEff(analysis_trigger, muonColl, electronColl);
     else if (DoubleMu_analysis) trigger_eff = TriggerEff(analysis_trigger, muonColl);
     else if (SingleMu_analysis) trigger_eff = TriggerEff(analysis_trigger, muonLooseColl);
     else if (DoubleEle_analysis) trigger_eff = TriggerEff(analysis_trigger, electronColl);
     //weight*=trigger_eff;
     //weight_nopu*=trigger_eff;
   }
   FillCutFlow("TriggerCut", weight); //;in trigger Cut mode
   /////////////////////////////////////////////////////////////////////////////////////////


   bool emu=false, mumu=false;
   int mu1_Ai=-1, mu2_Ai=-1, mu_Wi=-1;
   int nbjets=bjetColl.size(); int njets=jetColl.size(); const int nljets=ljetColl.size();//njets-nbjets;//number of light jets
   double Pzv1, Pzv2;

   double met = eventbase->GetEvent().MET();
   double met_x = eventbase->GetEvent().MET()*TMath::Cos(eventbase->GetEvent().METPhi());
   double met_y = eventbase->GetEvent().MET()*TMath::Sin(eventbase->GetEvent().METPhi());
   int Nvtx=eventbase->GetEvent().nVertices();
   double Pzv;
   snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));
//   snu::KParticle v[4]; v[0].SetPx(met_x); v[0].SetPy(met_y);
//                        v[1].SetPx(met_x); v[1].SetPy(met_y);

   //Scale Factors
   float id_weight_mu=1., reco_weight_mu=1., iso_weight_mu=1., trk_weight_mu=1., id_weight_ele=1.;
   if(!isData){
     id_weight_mu  *= MuonScaleFactor("MUON_POG_TIGHT", muonColl,0);
     iso_weight_mu *= MuonISOScaleFactor("MUON_POG_TIGHT", muonColl,0);
     trk_weight_mu *= MuonTrackingEffScaleFactor(muonColl);
     id_weight_ele *= ElectronScaleFactor(BaseSelection::ELECTRON_POG_TIGHT, electronColl);
//     id_weight *= MuonScaleFactor(BaseSelection::MUON_POG_TIGHT, muonColl,0);
//     iso_weight *= MuonISOScaleFactor(BaseSelection::MUON_POG_TIGHT, muonColl,0);
//     reco_weight *= ElectronRecoScaleFactor(electronColl);
//     weight*=id_weight*reco_weight*iso_weight;
   }





//////Basic Objects Distribution////////////////////////////////////////////////////////////////////////
   FillHist("Basic_Nj_orig", njets, weight, 0., 10., 10);
   FillHist("Basic_NjL_orig", jetLooseColl.size(), weight, 0., 10., 10);
   FillHist("Basic_Nb_orig", nbjets, weight, 0., 10., 10);//Bjet def CVSincV2>0.89;pfKJet::CSVv2IVF_discriminator 
   FillHist("Basic_NmuT_orig", muonColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NmuL_orig", muonLooseColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeT_orig", electronColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeL_orig", electronLooseColl.size(), weight, 0., 10., 10);
   if(electronColl.size()==1) FillHist("Basic_NmuT_weT", muonColl.size(), weight, 0., 10., 10);
   if(DoubleMu_analysis || SingleMu_analysis){
     if(muonColl.size()>0) FillHist("Basic_Ptmu1_3mu_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>1) FillHist("Basic_Ptmu2_3mu_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>2) FillHist("Basic_Ptmu3_3mu_orig", muonColl.at(2).Pt(), weight, 0., 200., 200);
   }
   if(EMu_analysis){
     if(muonColl.size()>0) FillHist("Basic_Ptmu1_1e2mu_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>1) FillHist("Basic_Ptmu2_1e2mu_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
     if(electronColl.size()>0) FillHist("Basic_Pte_1e2mu_orig", electronColl.at(0).Pt(), weight, 0., 200., 200);
   }
   if(njets>0) FillHist("Basic_j1_Et_orig", jetColl.at(0).Et(), weight, 0., 200., 200);
   if(njets>1) FillHist("Basic_j2_Et_orig", jetColl.at(1).Et(), weight, 0., 200., 200);
   if(njets>2) FillHist("Basic_j3_Et_orig", jetColl.at(2).Et(), weight, 0., 200., 200);
   if(njets>3) FillHist("Basic_j4_Et_orig", jetColl.at(3).Et(), weight, 0., 200., 200);
   if(nbjets>0) FillHist("Basic_b1_Et_orig", bjetColl.at(0).Et(), weight, 0, 200., 200);
   if(nbjets>1) FillHist("Basic_b2_Et_orig", bjetColl.at(1).Et(), weight, 0., 200., 200);
   FillHist("Basic_METdist_orig", met, weight, 0., 200., 100);//It is orig because even though I havent put METcut below, if there is correlation between Nl ,Nj, Nb and MET then MET distrubution will also be affected by selection. The distribution remains the same only when the correlation is absent.

///////Event Selection&Histograms//////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////   
 
  int Nlep=electronColl.size()+muonColl.size();
  float BTagWP=0.8484;

  //Curious 1: lepVeto needed? -> Yes, we don't want lepton to enter selection
  //Curious 2: lepton cut needed? -> Well, Let's check it.
  if(jetColl.size()>0){
    for( int i=0; i<jetColl.size(); i++){
      if(jetColl.at(i).HadronFlavour()==5){
        FillHist("NB_True_PT_Eta", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
        if(jetColl.at(i).BJetTaggerValue(snu::KJet::CSVv2)>BTagWP){
          FillHist("NBtag_TrueB_PT_Eta", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
        }
      }
      else if(jetColl.at(i).HadronFlavour()==4){
        FillHist("NC_True_PT_Eta", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
        if(jetColl.at(i).BJetTaggerValue(snu::KJet::CSVv2)>BTagWP){
          FillHist("NBtag_TrueC_PT_Eta", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
        }
      }
      else if(jetColl.at(i).HadronFlavour()==0){
        FillHist("NL_True_PT_Eta", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
        if(jetColl.at(i).BJetTaggerValue(snu::KJet::CSVv2)>BTagWP){
          FillHist("NBtag_TrueL_PT_Eta", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
        }
      }
    }
  }//End of Inclusive

  if(Nlep==1){
    if(jetColl.size()>0){
      for( int i=0; i<jetColl.size(); i++){
        if(jetColl.at(i).HadronFlavour()==5){
          FillHist("NB_True_PT_Eta_1l", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
          if(jetColl.at(i).BJetTaggerValue(snu::KJet::CSVv2)>BTagWP){
            FillHist("NBtag_TrueB_PT_Eta_1l", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
          }
        }
        else if(jetColl.at(i).HadronFlavour()==4){
          FillHist("NC_True_PT_Eta_1l", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
          if(jetColl.at(i).BJetTaggerValue(snu::KJet::CSVv2)>BTagWP){
            FillHist("NBtag_TrueC_PT_Eta_1l", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
          }
        }
        else if(jetColl.at(i).HadronFlavour()==0){
          FillHist("NL_True_PT_Eta_1l", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
          if(jetColl.at(i).BJetTaggerValue(snu::KJet::CSVv2)>BTagWP){
            FillHist("NBtag_TrueL_PT_Eta_1l", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
          }
        }
      }
    }
  }//End of lep 1 loop
  else if(Nlep==2){
    if(jetColl.size()>0){
      for( int i=0; i<jetColl.size(); i++){
        if(jetColl.at(i).HadronFlavour()==5){
          FillHist("NB_True_PT_Eta_2l", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
          if(jetColl.at(i).BJetTaggerValue(snu::KJet::CSVv2)>BTagWP){
            FillHist("NBtag_TrueB_PT_Eta_2l", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
          }
        }
        else if(jetColl.at(i).HadronFlavour()==4){
          FillHist("NC_True_PT_Eta_2l", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
          if(jetColl.at(i).BJetTaggerValue(snu::KJet::CSVv2)>BTagWP){
            FillHist("NBtag_TrueC_PT_Eta_2l", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
          }
        }
        else if(jetColl.at(i).HadronFlavour()==0){
          FillHist("NL_True_PT_Eta_2l", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
          if(jetColl.at(i).BJetTaggerValue(snu::KJet::CSVv2)>BTagWP){
            FillHist("NBtag_TrueL_PT_Eta_2l", jetColl.at(i).Pt(), jetColl.at(i).Eta(), 1, 20, 3000, 149, -2.4, 2.4, 8);
          }
        }
      }
    }
  }//End of lep 2 loop

 
return;
}// End of execute event loop
  


void Feb2017_3l4j_BTagEff::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Feb2017_3l4j_BTagEff::BeginCycle() throw( LQError ){
  
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

Feb2017_3l4j_BTagEff::~Feb2017_3l4j_BTagEff() {
  
  Message("In Feb2017_3l4j_BTagEff Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Feb2017_3l4j_BTagEff::FillCutFlow(TString cut, float weight){
  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 12, 0., 12.);
    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"NlCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"OS(2mu)");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"NjCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(8,"NbCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(9,"NljCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(10,"lPtCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(11,"MjjCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(12,"MmumuCut");
    
  }
}



void Feb2017_3l4j_BTagEff::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Feb2017_3l4j_BTagEff::MakeHistograms(){
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


  //2D Hist
  AnalyzerCore::MakeHistograms2D("NB_True_PT_Eta", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NB_True_PT_Eta")->SetTitle("N(j)(TruthB);PT(GeV);#eta");
    GetHist2D("NB_True_PT_Eta")->SetOption("textcolz");
  AnalyzerCore::MakeHistograms2D("NBtag_TrueB_PT_Eta", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NBtag_TrueB_PT_Eta")->SetTitle("N(j)(BtaggedTruthB);PT(GeV);#eta");
    GetHist2D("NBtag_TrueB_PT_Eta")->SetOption("textcolz");
  AnalyzerCore::MakeHistograms2D("NC_True_PT_Eta", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NC_True_PT_Eta")->SetTitle("N(j)(TruthC);PT(GeV);#eta");
    GetHist2D("NC_True_PT_Eta")->SetOption("textcolz");
  AnalyzerCore::MakeHistograms2D("NBtag_TrueC_PT_Eta", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NBtag_TrueC_PT_Eta")->SetTitle("N(j)(BtaggedTruthC);PT(GeV);#eta");
    GetHist2D("NBtag_TrueC_PT_Eta")->SetOption("textcolz");
  AnalyzerCore::MakeHistograms2D("NL_True_PT_Eta", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NL_True_PT_Eta")->SetTitle("N(j)(TruthL);PT(GeV);#eta");
    GetHist2D("NL_True_PT_Eta")->SetOption("textcolz");
  AnalyzerCore::MakeHistograms2D("NBtag_TrueL_PT_Eta", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NBtag_TrueL_PT_Eta")->SetTitle("N(j)(BtaggedTruthL);PT(GeV);#eta");
    GetHist2D("NBtag_TrueL_PT_Eta")->SetOption("textcolz");
 
  //1l category
  AnalyzerCore::MakeHistograms2D("NB_True_PT_Eta_1l", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NB_True_PT_Eta_1l")->SetTitle("N(j)(TruthB);PT(GeV);#eta");
    GetHist2D("NB_True_PT_Eta_1l")->SetOption("textcolz");
  AnalyzerCore::MakeHistograms2D("NBtag_TrueB_PT_Eta_1l", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NBtag_TrueB_PT_Eta_1l")->SetTitle("N(j)(BtaggedTruthB);PT(GeV);#eta");
    GetHist2D("NBtag_TrueB_PT_Eta_1l")->SetOption("textcolz");
  AnalyzerCore::MakeHistograms2D("NC_True_PT_Eta_1l", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NC_True_PT_Eta_1l")->SetTitle("N(j)(TruthC);PT(GeV);#eta");
    GetHist2D("NC_True_PT_Eta_1l")->SetOption("textcolz");
  AnalyzerCore::MakeHistograms2D("NBtag_TrueC_PT_Eta_1l", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NBtag_TrueC_PT_Eta_1l")->SetTitle("N(j)(BtaggedTruthC);PT(GeV);#eta");
    GetHist2D("NBtag_TrueC_PT_Eta_1l")->SetOption("textcolz");
  AnalyzerCore::MakeHistograms2D("NL_True_PT_Eta_1l", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NL_True_PT_Eta_1l")->SetTitle("N(j)(TruthL);PT(GeV);#eta");
    GetHist2D("NL_True_PT_Eta_1l")->SetOption("textcolz");
  AnalyzerCore::MakeHistograms2D("NBtag_TrueL_PT_Eta_1l", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NBtag_TrueL_PT_Eta_1l")->SetTitle("N(j)(BtaggedTruthL);PT(GeV);#eta");
    GetHist2D("NBtag_TrueL_PT_Eta_1l")->SetOption("textcolz");


  //2l category
  AnalyzerCore::MakeHistograms2D("NB_True_PT_Eta_2l", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NB_True_PT_Eta_2l")->SetTitle("N(j)(TruthB);PT(GeV);#eta");
    GetHist2D("NB_True_PT_Eta_2l")->SetOption("textcolz");
  AnalyzerCore::MakeHistograms2D("NBtag_TrueB_PT_Eta_2l", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NBtag_TrueB_PT_Eta_2l")->SetTitle("N(j)(BtaggedTruthB);PT(GeV);#eta");
    GetHist2D("NBtag_TrueB_PT_Eta_2l")->SetOption("textcolz");
  AnalyzerCore::MakeHistograms2D("NC_True_PT_Eta_2l", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NC_True_PT_Eta_2l")->SetTitle("N(j)(TruthC);PT(GeV);#eta");
    GetHist2D("NC_True_PT_Eta_2l")->SetOption("textcolz");
  AnalyzerCore::MakeHistograms2D("NBtag_TrueC_PT_Eta_2l", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NBtag_TrueC_PT_Eta_2l")->SetTitle("N(j)(BtaggedTruthC);PT(GeV);#eta");
    GetHist2D("NBtag_TrueC_PT_Eta_2l")->SetOption("textcolz");
  AnalyzerCore::MakeHistograms2D("NL_True_PT_Eta_2l", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NL_True_PT_Eta_2l")->SetTitle("N(j)(TruthL);PT(GeV);#eta");
    GetHist2D("NL_True_PT_Eta_2l")->SetOption("textcolz");
  AnalyzerCore::MakeHistograms2D("NBtag_TrueL_PT_Eta_2l", 149, 20, 3000, 8, -2.4, 2.4);
    GetHist2D("NBtag_TrueL_PT_Eta_2l")->SetTitle("N(j)(BtaggedTruthL);PT(GeV);#eta");
    GetHist2D("NBtag_TrueL_PT_Eta_2l")->SetOption("textcolz");


  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Feb2017_3l4j_BTagEffCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Feb2017_3l4j_BTagEff::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
