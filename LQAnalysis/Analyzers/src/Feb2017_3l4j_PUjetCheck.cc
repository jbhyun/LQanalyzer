// $Id: Feb2017_3l4j_PUjetCheck.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQFeb2017_3l4j_PUjetCheck Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Feb2017_3l4j_PUjetCheck.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Feb2017_3l4j_PUjetCheck);

 Feb2017_3l4j_PUjetCheck::Feb2017_3l4j_PUjetCheck() : AnalyzerCore(), out_muons(0) {

   SetLogName("Feb2017_3l4j_PUjetCheck");
   Message("In Feb2017_3l4j_PUjetCheck constructor", INFO);
   InitialiseAnalysis();
 }


 void Feb2017_3l4j_PUjetCheck::InitialiseAnalysis() throw( LQError ) {
   
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

void Feb2017_3l4j_PUjetCheck::ExecuteEvents()throw( LQError ){

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
     //weight*=pileup_reweight;
   } 
   //Numbet of Vertex PUreweight
   FillHist("Nvtx_nocut_PURW", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);


   bool EMu_analysis=false, SingleMu_analysis=false, DoubleMu_analysis=true, DoubleEle_analysis=false;


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
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);

     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(20.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.15);
   std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->Selection(muonColl, true);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
   //std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->SelectMuons(muonColl, BaseSelection::MUON_POG_TIGHT, 10., 2.4);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
   
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
     eventbase->GetElectronSel()->SetPt(20.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0994, 0.107);//ISO Ichep16 Loose WP
   //eventbase->GetElectronSel()->SetdxyBEMax(0.005, 0.01);  eventbase->GetElectronSel()->SetdzBEMax(0.005, 0.01);
//cout << eventbase->GetElectronSel()->pt_cut_min<<endl;
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);//pt>10/eta<2.4(hole region ~1.4 excluded by default)/EGAMMA_LOOSE/NoExtIso*/
//cout << eventbase->GetElectronSel()->pt_cut_min<<endl;
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(20.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//Iso Ichep16 Loose WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.2);
   std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->Selection(electronColl);//pt>10/eta<2.4(hole region ~1.4 excluded by default)/EGAMMA_LOOSE/NoExtIso*/
   //std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->SelectElectrons(electronColl, "ELECTRON_POG_TIGHT", 10., 2.5);//pt>10/eta<2.4/RelIso03<0.2/MUON_LOOSE
//cout<<electronLooseColl.size()<<" "<<electronColl.size()<<endl;
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(20.);                     eventbase->GetJetSel()->SetEta(2.4);
//   eventbase->GetJetSel()->SetUseJetPileUp(true);
     bool LeptonVeto=true;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonColl, electronColl);
   std::vector<snu::KJet> jetLooseColl; eventbase->GetJetSel()->SelectJets(jetLooseColl, muonColl, electronColl, "PFJET_LOOSE", 20., 2.4);


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
   float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1.;
   float trigger_sf=1.;

   if(!isData){
     //trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );

     id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronColl);
     reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);

     id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muonColl);
     iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
     trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
   }
   weight *= id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight;




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

//  TruthPrintOut();
   if(!isData){
   //if(false){
     std::vector<snu::KGenJet> genjetColl;
     eventbase->GetGenJetSel()->Selection(genjetColl);//PT>8 at ntuple level skim

     //PU Check
     bool SelectionPass=false;
     bool DY_analysis=true, EMuMu_analysis=false, TriMu_analysis=false;
     int PassCount=0;

     //For DY 2mu
     if(DY_analysis){
       if(muonColl.size()==2 && electronColl.size()==0 ){
         PassCount++;
         if(muonColl.at(0).Charge()!=muonColl.at(1).Charge()) PassCount++;
         if(muonColl.at(0).Pt()>20 && muonColl.at(1).Pt()>20) PassCount++;
       }
       if(PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")) PassCount++;     
       if(PassCount==4) SelectionPass=true;
     }

     //For Signal
     if(EMuMu_analysis){
       if(muonColl.size()==2 && electronColl.size()==1 ){
         PassCount++;
         if(muonColl.at(0).Charge()!=muonColl.at(1).Charge()) PassCount++;
         if(muonColl.at(0).Pt()>15 && muonColl.at(1).Pt()>10) PassCount++;
         if(electronColl.at(0).Pt()>25) PassCount++;
       }
       if(PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
         || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) PassCount++;
      
       if(PassCount==5) SelectionPass=true;
     }
     if(TriMu_analysis){
       if(muonColl.size()==3){
         PassCount++;
         if(muonColl.at(0).Charge()!=muonColl.at(1).Charge()) PassCount++;
         if(muonColl.at(0).Pt()>20 && muonColl.at(2).Pt()>10) PassCount++;
       }
       if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ) PassCount++;
       if(PassCount==4) SelectionPass=true;
     }

   //  weight=MCweight;

     if(!SelectionPass) return;

     //jet-genjet Matching
     std::vector<int> MatchedGenJetIdxLooseColl;
     std::vector<int> MatchedGenJetIdxTightColl;
     std::vector<int> MatchedPartonIdxColl;
     std::vector<int> JetCMSFlavColl;
     for(int i=0; i<jetColl.size(); i++){
       float dRminjgenj=10000, dRminjq=10000, dR=10000, dPtRel1=10000, dPtRel2=10000;
       int MatchedIdx_genj=-1, MatchedIdx_gq=-1;
       //genjet matching
       for(int j=0; j<genjetColl.size(); j++){         
         dR=jetColl.at(i).DeltaR(genjetColl.at(j));
         if(dR<dRminjgenj){
           dRminjgenj=dR; MatchedIdx_genj=j;
         }
       }
       //Parton matching
       for(int j=0; j<truthColl.size(); j++){
         if( !(truthColl.at(j).GenStatus()==23 || truthColl.at(j).GenStatus()==3) ) continue;
         if( !(fabs(truthColl.at(j).PdgId())<10 || fabs(truthColl.at(j).PdgId())==21 ) ) continue;
         if( truthColl.at(j).Pt()==0 ) continue;
 
         dR=jetColl.at(i).DeltaR(truthColl.at(j));
         dPtRel1=fabs(jetColl.at(i).Pt()-truthColl.at(j).Pt())/jetColl.at(i).Pt();
         dPtRel2=fabs(jetColl.at(i).Pt()-truthColl.at(j).Pt())/truthColl.at(j).Pt();
         if( dPtRel1>3 || dPtRel2>3 ) continue;
         if( dR>0.4 ) continue;

         if(dR<dRminjq){
           dRminjq=dR; MatchedIdx_gq=j;
         }         
       }
       if(dRminjgenj<0.2){
         //Maximum tolerable level: 0.2 : from 0.4dRcone of jet && PFjet dR uncertainty ~0.03
         //And also prevents duplicate matching
         MatchedGenJetIdxLooseColl.push_back(MatchedIdx_genj);
       }
       else{
         MatchedGenJetIdxLooseColl.push_back(-1);
         FillHist("dRminjPULjgen", dRminjgenj, weight, 0., 5., 500);
       }

       //if(dRminjq<0.4 && dRminjgenj<0.2 && jetColl.at(i).PartonFlavour()!=0 ){
       if(dRminjgenj<0.2 && (jetColl.at(i).PartonFlavour()!=0 || jetColl.at(i).HadronFlavour()!=0)){
       //if(dRminjgenj<0.2){
         //Identified Physical Jet
         MatchedGenJetIdxTightColl.push_back(MatchedIdx_genj);
         MatchedPartonIdxColl.push_back(MatchedIdx_gq);
         JetCMSFlavColl.push_back(jetColl.at(i).PartonFlavour());
       }
       //else if (dRminjq>0.4 && dRminjgenj>0.2 && jetColl.at(i).PartonFlavour()==0){
       else if (dRminjgenj>0.2 && jetColl.at(i).PartonFlavour()==0 && jetColl.at(i).HadronFlavour()==0 && MatchedIdx_gq==-1){
         //Identified PU Jet
         MatchedGenJetIdxTightColl.push_back(-1);
         MatchedPartonIdxColl.push_back(-1);
         JetCMSFlavColl.push_back(jetColl.at(i).PartonFlavour());        

         FillHist("dRminjPUTjgen", dRminjgenj, weight, 0., 5., 500);
       }
       else{
         //Else case : either PU jet but close to physical jet OR physical jet severely affected by PU activity
         MatchedGenJetIdxTightColl.push_back(-2);
         MatchedPartonIdxColl.push_back(-2);
         JetCMSFlavColl.push_back(jetColl.at(i).PartonFlavour());        


         FillHist("dRminjElseTjgen", dRminjgenj, weight, 0., 5., 500);
       }

     }

     //Property Check
     for(int i=0; i<jetColl.size(); i++){

       //Loops with Loose classification Coll
       if(MatchedGenJetIdxLooseColl.at(i) != -1){
         FillHist("dRjjgen_L", jetColl.at(i).DeltaR(genjetColl.at(MatchedGenJetIdxLooseColl.at(i))), weight, 0., 0.2, 20);

         FillHist("RealJet_PT_Eta_L", jetColl.at(i).Pt(), jetColl.at(i).Eta(), weight, 20., 500., 48, -2.4, 2.4, 12);
         FillHist("RealJet_PT_L", jetColl.at(i).Pt(), weight, 0., 300., 300);
         FillHist("RealJet_Eta_L", jetColl.at(i).Eta(), weight, -2.4, 2.4, 48);

         if(k_sample_name.Contains("DYJets")){
           if((muonColl.at(0)+muonColl.at(1)).Pt()>0){
             FillHist("RealJet_dPhiZj_L", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), weight, -3.2, 3.2, 128);
             FillHist("RealJet_dPhiZj_PTrelZ_L", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), jetColl.at(i).Pt()/(muonColl.at(0)+muonColl.at(1)).Pt(), weight, -3.2, 3.2, 64, 0., 10., 20);
           }
           if(fabs((muonColl.at(0)+muonColl.at(1)).DeltaPhi(jetColl.at(i)))>2.5){
             FillHist("RealJet_dRjgenj_RealEnriched_L", jetColl.at(i).DeltaR(genjetColl.at(MatchedGenJetIdxLooseColl.at(i))), weight, 0., 0.2, 20); 
           }
           else if(fabs((muonColl.at(0)+muonColl.at(1)).DeltaPhi(jetColl.at(i)))<1.5){
             FillHist("RealJet_dRjgenj_PUEnriched_L", jetColl.at(i).DeltaR(genjetColl.at(MatchedGenJetIdxLooseColl.at(i))), weight, 0., 0.2, 20); 
           }
         }//DYLoop Ends

       }
       else{
         FillHist("PUJet_PT_Eta_L", jetColl.at(i).Pt(), jetColl.at(i).Eta(), weight, 20., 500., 48, -2.4, 2.4, 12);
         FillHist("PUJet_PT_L", jetColl.at(i).Pt(), weight, 0., 300., 300);
         FillHist("PUJet_Eta_L", jetColl.at(i).Eta(), weight, -2.4, 2.4, 48);

         if(k_sample_name.Contains("DYJets")){
           if((muonColl.at(0)+muonColl.at(1)).Pt()>0){
             FillHist("PUJet_dPhiZj_L", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), weight, -3.2, 3.2, 128);
             FillHist("PUJet_dPhiZj_PTrelZ_L", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), jetColl.at(i).Pt()/(muonColl.at(0)+muonColl.at(1)).Pt(), weight, -3.2, 3.2, 64, 0., 10., 20);     }
         }

       }

       //Loops with Tight classification Coll
       if(MatchedGenJetIdxTightColl.at(i) > 0){
         float jetpt=jetColl.at(i).Pt();
         float jetmva=jetColl.at(i).PileupJetIDMVA();

         //Quite Confidently Physical Jets
         //Gen matching quality 
         FillHist("dRjjgen_T", jetColl.at(i).DeltaR(genjetColl.at(MatchedGenJetIdxTightColl.at(i))), weight, 0., 0.2, 20);
         if(MatchedPartonIdxColl.at(i)>0){
           FillHist("dRjq_T", jetColl.at(i).DeltaR(truthColl.at(MatchedPartonIdxColl.at(i))), weight, 0., 0.4, 40);
           FillHist("PtqoPtj_T", truthColl.at(MatchedPartonIdxColl.at(i)).Pt()/jetColl.at(i).Pt(), weight, 0., 4., 400);
         }

        //Real Jet Basic Prop
         FillHist("RealJet_PT_Eta_T", jetColl.at(i).Pt(), jetColl.at(i).Eta(), weight, 20., 500., 48, -2.4, 2.4, 12);
         FillHist("RealJet_PT_T", jetColl.at(i).Pt(), weight, 0., 300., 300);
         FillHist("RealJet_Eta_T", jetColl.at(i).Eta(), weight, -2.4, 2.4, 48);
         FillHist("RealJet_Count_T", 0., weight, 0., 1., 1);

         //Real Jet Prop in separation region(only DY)
         if(k_sample_name.Contains("DYJets")){
           if((muonColl.at(0)+muonColl.at(1)).Pt()>0){
             FillHist("RealJet_dPhiZj_T", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), weight, -3.2, 3.2, 128);
             FillHist("RealJet_dPhiZj_PTrelZ_T", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), jetColl.at(i).Pt()/(muonColl.at(0)+muonColl.at(1)).Pt(), weight, -3.2, 3.2, 64, 0., 10., 20);
             if(fabs((muonColl.at(0)+muonColl.at(1)).DeltaPhi(jetColl.at(i)))>2.5){
               FillHist("RealJet_dRjgenj_RealEnriched_T", jetColl.at(i).DeltaR(genjetColl.at(MatchedGenJetIdxTightColl.at(i))), weight, 0., 0.2, 20); 
             }
             else if(fabs((muonColl.at(0)+muonColl.at(1)).DeltaPhi(jetColl.at(i)))<1.5){
               FillHist("RealJet_dRjgenj_PUEnriched_T", jetColl.at(i).DeltaR(genjetColl.at(MatchedGenJetIdxTightColl.at(i))), weight, 0., 0.2, 20); 
             }
           }
         }//DYLoop Ends

         //Realjets passing MVA WP
         if((jetpt<30 && jetmva>-0.97) || (jetpt>=30 && jetmva>-0.89)){
           FillHist("RealJet_PT_TmvaL", jetColl.at(i).Pt(), weight, 0., 300., 300);
           FillHist("RealJet_Eta_TmvaL", jetColl.at(i).Eta(), weight, -2.4, 2.4, 48);
           FillHist("RealJet_Count_TmvaL", 0., weight, 0., 1., 1);
           if(k_sample_name.Contains("DYJets") && (muonColl.at(0)+muonColl.at(1)).Pt()>0){
             FillHist("RealJet_dPhiZj_TmvaL", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), weight, -3.2, 3.2, 128);
           }
         }
         if((jetpt<30 && jetmva>0.18) || (jetpt>=30 && jetmva>0.61)){
           FillHist("RealJet_PT_TmvaM", jetColl.at(i).Pt(), weight, 0., 300., 300);
           FillHist("RealJet_Eta_TmvaM", jetColl.at(i).Eta(), weight, -2.4, 2.4, 48);
           FillHist("RealJet_Count_TmvaM", 0., weight, 0., 1., 1);
           if(k_sample_name.Contains("DYJets") && (muonColl.at(0)+muonColl.at(1)).Pt()>0){
             FillHist("RealJet_dPhiZj_TmvaM", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), weight, -3.2, 3.2, 128);
           }
         }
         if((jetpt<30 && jetmva>0.69) || (jetpt>=30 && jetmva>0.86)){
           FillHist("RealJet_PT_TmvaT", jetColl.at(i).Pt(), weight, 0., 300., 300);
           FillHist("RealJet_Eta_TmvaT", jetColl.at(i).Eta(), weight, -2.4, 2.4, 48);
           FillHist("RealJet_Count_TmvaT", 0., weight, 0., 1., 1);
           if(k_sample_name.Contains("DYJets") && (muonColl.at(0)+muonColl.at(1)).Pt()>0){
             FillHist("RealJet_dPhiZj_TmvaT", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), weight, -3.2, 3.2, 128);
           }
         }

       }
       else if(MatchedGenJetIdxTightColl.at(i)==-1){
         //Quite confidently PU Jets
         //PUJet Basic Prop
         FillHist("PUJet_PT_Eta_T", jetColl.at(i).Pt(), jetColl.at(i).Eta(), weight, 20., 500., 48, -2.4, 2.4, 12);
         FillHist("PUJet_PT_T", jetColl.at(i).Pt(), weight, 0., 300., 300);
         FillHist("PUJet_Eta_T", jetColl.at(i).Eta(), weight, -2.4, 2.4, 48);
         FillHist("PUJet_Count_T", 0., weight, 0., 1., 1);

         //PUJet Prop in separation region(only DY)
         if(k_sample_name.Contains("DYJets")){
           if((muonColl.at(0)+muonColl.at(1)).Pt()>0){
             FillHist("PUJet_dPhiZj_T", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), weight, -3.2, 3.2, 128);
             FillHist("PUJet_dPhiZj_PTrelZ_T", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), jetColl.at(i).Pt()/(muonColl.at(0)+muonColl.at(1)).Pt(), weight, -3.2, 3.2, 64, 0., 10., 20);
           }
         }

         //PUJets passing MVA WP
         float jetpt=jetColl.at(i).Pt();
         float jetmva=jetColl.at(i).PileupJetIDMVA();
         if((jetpt<30 && jetmva>-0.97) || (jetpt>=30 && jetmva>-0.89)){
           FillHist("PUJet_PT_TmvaL", jetColl.at(i).Pt(), weight, 0., 300., 300);
           FillHist("PUJet_Eta_TmvaL", jetColl.at(i).Eta(), weight, -2.4, 2.4, 48);
           FillHist("PUJet_Count_TmvaL", 0., weight, 0., 1., 1);
           if(k_sample_name.Contains("DYJets") && (muonColl.at(0)+muonColl.at(1)).Pt()>0){
             FillHist("PUJet_dPhiZj_TmvaL", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), weight, -3.2, 3.2, 128);
           }

         }
         if((jetpt<30 && jetmva>0.18) || (jetpt>=30 && jetmva>0.61)){
           FillHist("PUJet_PT_TmvaM", jetColl.at(i).Pt(), weight, 0., 300., 300);
           FillHist("PUJet_Eta_TmvaM", jetColl.at(i).Eta(), weight, -2.4, 2.4, 48);
           FillHist("PUJet_Count_TmvaM", 0., weight, 0., 1., 1);
           if(k_sample_name.Contains("DYJets") && (muonColl.at(0)+muonColl.at(1)).Pt()>0){
             FillHist("PUJet_dPhiZj_TmvaM", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), weight, -3.2, 3.2, 128);
           }
         }
         if((jetpt<30 && jetmva>0.69) || (jetpt>=30 && jetmva>0.86)){
           FillHist("PUJet_PT_TmvaT", jetColl.at(i).Pt(), weight, 0., 300., 300);
           FillHist("PUJet_Eta_TmvaT", jetColl.at(i).Eta(), weight, -2.4, 2.4, 48);
           FillHist("PUJet_Count_TmvaT", 0., weight, 0., 1., 1);
           if(k_sample_name.Contains("DYJets") && (muonColl.at(0)+muonColl.at(1)).Pt()>0){
             FillHist("PUJet_dPhiZj_TmvaT", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), weight, -3.2, 3.2, 128);
           }
         }

       }
       else{
         //Not confident. could be either physical jets heavily affected by PU activity or PU jets close to physical jets
         //Else Jet Basic Prop
         FillHist("ElseJet_PT_Eta_T", jetColl.at(i).Pt(), jetColl.at(i).Eta(), weight, 20., 500., 48, -2.4, 2.4, 12);
         FillHist("ElseJet_PT_T", jetColl.at(i).Pt(), weight, 0., 300., 300);
         FillHist("ElseJet_Eta_T", jetColl.at(i).Eta(), weight, -2.4, 2.4, 48);
         FillHist("ElseJet_Count_T", 0., weight, 0., 1., 1);

         //ElseJet Prop in separation region ( only DY)
         if(k_sample_name.Contains("DYJets")){
           if((muonColl.at(0)+muonColl.at(1)).Pt()>0){
             FillHist("ElseJet_dPhiZj_T", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), weight, -3.2, 3.2, 128);
             FillHist("ElseJet_dPhiZj_PTrelZ_T", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), jetColl.at(i).Pt()/(muonColl.at(0)+muonColl.at(1)).Pt(), weight, -3.2, 3.2, 64, 0., 10., 20);
           }
         }

         //Else jets passing MVA WP
         float jetpt=jetColl.at(i).Pt();
         float jetmva=jetColl.at(i).PileupJetIDMVA();
         if((jetpt<30 && jetmva>-0.97) || (jetpt>=30 && jetmva>-0.89)){
           FillHist("ElseJet_PT_TmvaL", jetColl.at(i).Pt(), weight, 0., 300., 300);
           FillHist("ElseJet_Eta_TmvaL", jetColl.at(i).Eta(), weight, -2.4, 2.4, 48);
           FillHist("ElseJet_Count_TmvaL", 0., weight, 0., 1., 1);
           if(k_sample_name.Contains("DYJets") && (muonColl.at(0)+muonColl.at(1)).Pt()>0){
             FillHist("ElseJet_dPhiZj_TmvaL", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), weight, -3.2, 3.2, 128);
           }

         }
         if((jetpt<30 && jetmva>0.18) || (jetpt>=30 && jetmva>0.61)){
           FillHist("ElseJet_PT_TmvaM", jetColl.at(i).Pt(), weight, 0., 300., 300);
           FillHist("ElseJet_Eta_TmvaM", jetColl.at(i).Eta(), weight, -2.4, 2.4, 48);
           FillHist("ElseJet_Count_TmvaM", 0., weight, 0., 1., 1);
           if(k_sample_name.Contains("DYJets") && (muonColl.at(0)+muonColl.at(1)).Pt()>0){
             FillHist("ElseJet_dPhiZj_TmvaM", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), weight, -3.2, 3.2, 128);
           }
         }
         if((jetpt<30 && jetmva>0.69) || (jetpt>=30 && jetmva>0.86)){
           FillHist("ElseJet_PT_TmvaT", jetColl.at(i).Pt(), weight, 0., 300., 300);
           FillHist("ElseJet_Eta_TmvaT", jetColl.at(i).Eta(), weight, -2.4, 2.4, 48);
           FillHist("ElseJet_Count_TmvaT", 0., weight, 0., 1., 1);
           if(k_sample_name.Contains("DYJets") && (muonColl.at(0)+muonColl.at(1)).Pt()>0){
             FillHist("ElseJet_dPhiZj_TmvaT", jetColl.at(i).DeltaPhi((muonColl.at(0)+muonColl.at(1))), weight, -3.2, 3.2, 128);
           }
         }


       }

     }//Property Loop over jets Ends


   }//MC Lines End


 
return;
}// End of execute event loop
  


void Feb2017_3l4j_PUjetCheck::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Feb2017_3l4j_PUjetCheck::BeginCycle() throw( LQError ){
  
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

Feb2017_3l4j_PUjetCheck::~Feb2017_3l4j_PUjetCheck() {
  
  Message("In Feb2017_3l4j_PUjetCheck Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Feb2017_3l4j_PUjetCheck::FillCutFlow(TString cut, float weight){
  
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



void Feb2017_3l4j_PUjetCheck::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Feb2017_3l4j_PUjetCheck::MakeHistograms(){
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


  //Analysis Histograms
  AnalyzerCore::MakeHistograms2D("RealJet_PT_Eta_L"       , 48, 20., 500., 12, -2.4, 2.4);
    GetHist2D("RealJet_PT_Eta_L")->SetOption("colz");
  AnalyzerCore::MakeHistograms2D("RealJet_dPhiZj_PTrelZ_L", 64, -3.2, 3.2, 20, 0., 10.);
    GetHist2D("RealJet_dPhiZj_PTrelZ_L")->SetOption("colz");
  AnalyzerCore::MakeHistograms2D("PUJet_PT_Eta_L"         , 48, 20., 500., 12, -2.4, 2.4);
    GetHist2D("PUJet_PT_Eta_L")->SetOption("colz");
  AnalyzerCore::MakeHistograms2D("PUJet_dPhiZj_PTrelZ_L"  , 64, -3.2, 3.2, 20, 0., 10.);
    GetHist2D("PUJet_dPhiZj_PTrelZ_L")->SetOption("colz");
  AnalyzerCore::MakeHistograms2D("RealJet_PT_Eta_T"       , 48, 20., 500., 12, -2.4, 2.4);
    GetHist2D("RealJet_PT_Eta_T")->SetOption("colz");
  AnalyzerCore::MakeHistograms2D("RealJet_dPhiZj_PTrelZ_T", 64, -3.2, 3.2, 20, 0., 10.);
    GetHist2D("RealJet_dPhiZj_PTrelZ_T")->SetOption("colz");
  AnalyzerCore::MakeHistograms2D("PUJet_PT_Eta_T"         , 48, 20., 500., 12, -2.4, 2.4);
    GetHist2D("PUJet_PT_Eta_T")->SetOption("colz");
  AnalyzerCore::MakeHistograms2D("PUJet_dPhiZj_PTrelZ_T"  , 64, -3.2, 3.2, 20, 0., 10.);
    GetHist2D("PUJet_dPhiZj_PTrelZ_T")->SetOption("colz");
  AnalyzerCore::MakeHistograms2D("ElseJet_PT_Eta_T"       , 48, 20., 500., 12, -2.4, 2.4);
    GetHist2D("ElseJet_PT_Eta_T")->SetOption("colz");
  AnalyzerCore::MakeHistograms2D("ElseJet_dPhiZj_PTrelZ_T", 64, -3.2, 3.2, 20, 0., 10.);
    GetHist2D("ElseJet_dPhiZj_PTrelZ_T")->SetOption("colz");




  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Feb2017_3l4j_PUjetCheckCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Feb2017_3l4j_PUjetCheck::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
