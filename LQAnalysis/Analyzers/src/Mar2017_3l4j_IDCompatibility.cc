// $Id: Mar2017_3l4j_IDCompatibility.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQMar2017_3l4j_IDCompatibility Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Mar2017_3l4j_IDCompatibility.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Mar2017_3l4j_IDCompatibility);

 Mar2017_3l4j_IDCompatibility::Mar2017_3l4j_IDCompatibility() : AnalyzerCore(), out_muons(0) {

   SetLogName("Mar2017_3l4j_IDCompatibility");
   Message("In Mar2017_3l4j_IDCompatibility constructor", INFO);
   InitialiseAnalysis();
 }


 void Mar2017_3l4j_IDCompatibility::InitialiseAnalysis() throw( LQError ) {
   
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

void Mar2017_3l4j_IDCompatibility::ExecuteEvents()throw( LQError ){

////Basic Infos///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Apply the gen weight 
   float weight_nopu=1; 
   if(!isData) weight*=MCweight;
   FillHist("GenWeight", MCweight, 1, -2., 2., 4);

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
   if (!k_isdata) { pileup_reweight=eventbase->GetEvent().PileUpWeight(); }

   //Numbet of Vertex PUreweight
   FillHist("Nvtx_nocut_PURW", eventbase->GetEvent().nVertices(), weight*pileup_reweight, 0., 50., 50);


   //bool TriMu_analysis=true, EMuMu_analysis=false, DiMuon_analysis=false, DiEle_analysis=false;
   //bool TriMu_analysis=false, EMuMu_analysis=true, DiMuon_analysis=false, DiEle_analysis=false;
   bool TriMu_analysis=false, EMuMu_analysis=false, DiMuon_analysis=true, DiEle_analysis=false;
   //bool TriMu_analysis=false, EMuMu_analysis=false, DiMuon_analysis=false, DiEle_analysis=true;

   //Trigers
   std::vector<TString> triggerslist, triggerlist1, triggerlist2;//  ListTriggersAvailable();

   TString analysis_trigger;
   if(EMuMu_analysis) {analysis_trigger="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";}
   else if(TriMu_analysis) analysis_trigger="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";
   //if(!PassTrigger(analysis_trigger)) return; //Not for now..
   //PassTrigger Uses BeginWith so don't stick to exact name


   float trigger_ps_weight=1;
   trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
    weight*=trigger_ps_weight;
    weight_nopu*=trigger_ps_weight;
   FillHist("TriggerPSWeight", trigger_ps_weight, 1., 0., 1., 100);

   //Initial Event Cut : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
   //METFilter
   if(!PassMETFilter()) return;  FillCutFlow("EventCut", weight);

   //Vertex Cut
   //Good Primary vtx def:(vtx.ndof()>4&&maxAbsZ<=0)||std::abs(vtx.z())<= 24)
   //                     &&((maxd0 <=0) || std::abs(vtx.position().rho())<=2)&&!(vtx.isFake()))  
   if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; FillCutFlow("VertexCut", weight);
   


///////Objects in Analysis/////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

   //Primary Object Collection
   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
  
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.25);//POG WP L
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);//(muonColl, bool RochCorr, bool debug)
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(5.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.15);//POG WP T
   std::vector<snu::KMuon> muonColl; eventbase->GetMuonSel()->Selection(muonColl,true);//(muonColl, bool RochCorr, bool debug)
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(5.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.05);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.15);//POG WP T
   std::vector<snu::KMuon> muonHNColl; eventbase->GetMuonSel()->Selection(muonHNColl,true);//(muonColl, bool RochCorr, bool debug)
 //std::vector<snu::KMuon> muonPromptColl; muonPromptColl=GetTruePrompt(muonColl, false);

     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
     eventbase->GetElectronSel()->SetPt(20.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0994, 0.107);//2016 80X tuned WP
   //eventbase->GetElectronSel()->SetdxyBEMax(0.005, 0.01);  eventbase->GetElectronSel()->SetdzBEMax(0.005, 0.01);//Not in ID, but additional safe WP
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(20.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.2);//Not in ID, but additional safe WP
   std::vector<snu::KElectron> electronColl; eventbase->GetElectronSel()->Selection(electronColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(20.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.0588, 0.0571);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.2);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetCheckCharge(true);
   std::vector<snu::KElectron> electronHNColl; eventbase->GetElectronSel()->Selection(electronHNColl);

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(20.);                     eventbase->GetJetSel()->SetEta(2.4);
     //eventbase->GetJetSel()->SetPileUpJetID(true,"Loose");
     bool LeptonVeto=true;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonColl, electronColl);
   //std::vector<snu::KJet> jetLooseColl; eventbase->GetJetSel()->SelectJets(jetLooseColl, muonColl, electronColl, "PFJET_LOOSE", 20., 2.4);

   //Method to apply 2a method SF
   //std::vector<int> bIdxColl  = GetSFBJetIdx(jetColl,"Medium");
   //std::vector<int> ljIdxColl = GetSFLJetIdx(jetColl, bIdxColl, "Medium");
   //std::vector<snu::KJet> bjetColl2a; for(int i=0; i<bIdxColl.size(); i++) {bjetColl2a.push_back(jetColl.at(bIdxColl.at(i)));}
   //std::vector<snu::KJet> ljetColl2a; for(int i=0; i<ljIdxColl.size(); i++){ljetColl2a.push_back(jetColl.at(ljIdxColl.at(i)));}

   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");
   std::vector<snu::KJet> ljetColl = SelLightJets(jetColl, "Medium");

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
   float btag_sf=1.;
   float trigger_sf=1.;

   if(!isData){
     trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );

     id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronColl);
     reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);

     id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muonColl);
     iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
     trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);

     //btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);
   }


//////Basic Objects Distribution////////////////////////////////////////////////////////////////////////
   FillHist("Basic_Nj_orig", njets, weight, 0., 10., 10);
   FillHist("Basic_Nb_orig", nbjets, weight, 0., 10., 10);//Bjet def CVSincV2>0.89;pfKJet::CSVv2IVF_discriminator 
   FillHist("Basic_NmuT_orig", muonColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NmuL_orig", muonLooseColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeT_orig", electronColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeL_orig", electronLooseColl.size(), weight, 0., 10., 10);
   if(electronColl.size()==1) FillHist("Basic_NmuT_weT", muonColl.size(), weight, 0., 10., 10);
   if(TriMu_analysis){
     if(muonColl.size()>0) FillHist("Basic_Ptmu1_3mu_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>1) FillHist("Basic_Ptmu2_3mu_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>2) FillHist("Basic_Ptmu3_3mu_orig", muonColl.at(2).Pt(), weight, 0., 200., 200);
   }
   if(EMuMu_analysis){
     if(muonColl.size()>0) FillHist("Basic_Ptmu1_1e2mu_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>1) FillHist("Basic_Ptmu2_1e2mu_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
     if(electronColl.size()>0) FillHist("Basic_Pte_1e2mu_orig", electronColl.at(0).Pt(), weight, 0., 200., 200);
   }
   if(njets>0)  FillHist("Basic_j1_Et_orig", jetColl.at(0).Et(), weight, 0., 200., 200);
   if(njets>1)  FillHist("Basic_j2_Et_orig", jetColl.at(1).Et(), weight, 0., 200., 200);
   if(njets>2)  FillHist("Basic_j3_Et_orig", jetColl.at(2).Et(), weight, 0., 200., 200);
   if(njets>3)  FillHist("Basic_j4_Et_orig", jetColl.at(3).Et(), weight, 0., 200., 200);
   if(nbjets>0) FillHist("Basic_b1_Et_orig", bjetColl.at(0).Et(), weight, 0, 200., 200);
   if(nbjets>1) FillHist("Basic_b2_Et_orig", bjetColl.at(1).Et(), weight, 0., 200., 200);
   FillHist("Basic_METdist_orig", met, weight, 0., 200., 100);//It is orig because even though I havent put METcut below, if there is correlation between Nl ,Nj, Nb and MET then MET distrubution will also be affected by selection. The distribution remains the same only when the correlation is absent.

///////Event Selection&Histograms//////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////   


   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight;
   if(EMuMu_analysis){
     //if( !(PassTrigger("HLT_IsoMu24_v")||PassTrigger("HLT_IsoTkMu24_v")) ) return;
     if( !(PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
          || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")
          || PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v")) ) return;


     //Step1 : 1e+2mu +PTCut
     //if( !(muonColl.at(0).Pt()>27) ) return; //SingleMuon Trig Case
     if( !(electronColl.size()==1 && muonColl.size()==2) ) return;
     if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>15) ) return;
     if( SumCharge(muonColl)!=0 ) return;

     //Overall Difference
     FillHist("NevtDiff_1e2mu", 0., weight, 0., 10., 10);
     if(electronHNColl.size()==1 && muonHNColl.size()==2){
       if(electronHNColl.at(0).Pt()>25 && muonHNColl.at(0).Pt()>15 && SumCharge(muonColl)==0){
         FillHist("NevtDiff_1e2mu", 1., weight, 0., 10., 10);
       }
     }

     //Ele ID Property
     FillHist("IDflow_e_1e2mu", 0., weight, 0., 10., 10);
     if(fabs(electronColl.at(0).Eta())<1.479){
       FillHist("Absd0_eB_1e2mu", fabs(electronColl.at(0).dxy()), weight, 0., 0.05, 500);
       FillHist("Absdz_eB_1e2mu", fabs(electronColl.at(0).dz()), weight, 0., 0.1, 1000);
       
       if(fabs(electronColl.at(0).dxy())<0.0111){
         FillHist("IDflow_e_1e2mu", 1., weight, 0., 10., 10);
         if(fabs(electronColl.at(0).dz())<0.0466){
           FillHist("IDflow_e_1e2mu", 2., weight, 0., 10., 10);
           if(fabs(electronColl.at(0).dxySig())<3.){
             FillHist("IDflow_e_1e2mu", 3., weight, 0., 10., 10);
             if(electronColl.at(0).GsfCtfScPixChargeConsistency()){
               FillHist("IDflow_e_1e2mu", 4., weight, 0., 10., 10);
             }
           }
         }
       }
     }
     else{
       FillHist("Absd0_eE_1e2mu", fabs(electronColl.at(0).dxy()), weight, 0., 0.1, 1000);
       FillHist("Absdz_eE_1e2mu", fabs(electronColl.at(0).dz()), weight, 0., 0.2, 2000);

       if(fabs(electronColl.at(0).dxy())<0.0351){
         FillHist("IDflow_e_1e2mu", 1., weight, 0., 10., 10);
         if(fabs(electronColl.at(0).dz())<0.417){
           FillHist("IDflow_e_1e2mu", 2., weight, 0., 10., 10);
           if(fabs(electronColl.at(0).dxySig())<3.){
             FillHist("IDflow_e_1e2mu", 3., weight, 0., 10., 10);
             if(electronColl.at(0).GsfCtfScPixChargeConsistency()){
               FillHist("IDflow_e_1e2mu", 4., weight, 0., 10., 10);
             }
           }
         }
       }
     }
     FillHist("Absd0sig_e_1e2mu", fabs(electronColl.at(0).dxySig()), weight, 0., 10., 100);


     //Muon ID Property
     float RochIso04=999.;
     for(unsigned int i=0; i<muonColl.size(); i++){
       RochIso04=muonColl.at(i).RelIso04()/muonColl.at(i).RochPt()*muonColl.at(i).Pt();

       FillHist("IDflow_mu_1e2mu", 0., weight, 0., 10., 10);
       if(fabs(muonColl.at(i).dXY())<0.05){
         FillHist("IDflow_mu_1e2mu", 1., weight, 0., 10., 10);
         if(fabs(muonColl.at(i).dXYSig())<3.){
           FillHist("IDflow_mu_1e2mu", 2., weight, 0., 10., 10);
           if(RochIso04<0.1) FillHist("IDflow_mu_1e2mu", 3., weight, 0., 10., 10); 
         }
       }
       FillHist("Absd0_mu_1e2mu", fabs(muonColl.at(i).dXY()), weight, 0., 0.2, 2000);
       FillHist("Absd0sig_mu_1e2mu", fabs(muonColl.at(i).dXYSig()), weight, 0., 10., 100);
       FillHist("RelIso04_mu_1e2mu", RochIso04, weight, 0., 0.2, 200);

     }
   }
   if(TriMu_analysis){
     //if( !(PassTrigger("HLT_IsoMu24_v")||PassTrigger("HLT_IsoTkMu24_v")) ) return;
     if( !(PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")) ) return;


     //Step1 : 3mu +PTCut
     //if( !(muonColl.at(0).Pt()>27) ) return; //SingleMuon Trig Case
     if( !(muonColl.size()==3) ) return;
     if( !(muonColl.at(0).Pt()>20 && muonColl.at(2).Pt()>10) ) return;
     if( fabs(SumCharge(muonColl))!=1 ) return;

     //Overall Difference
     FillHist("NevtDiff_3mu", 0., weight, 0., 10., 10);
     if( muonHNColl.size()==3 ){
       if( muonHNColl.at(0).Pt()>20 && muonHNColl.at(2).Pt()>10 && fabs(SumCharge(muonHNColl))==1){
         FillHist("NevtDiff_3mu", 1., weight, 0., 10., 10);
       }
     }
         

     //Muon ID Property
     float RochIso04=999.;
     for(unsigned int i=0; i<muonColl.size(); i++){
       RochIso04=muonColl.at(i).RelIso04()/muonColl.at(i).RochPt()*muonColl.at(i).Pt();

       FillHist("IDflow_mu_3mu", 0., weight, 0., 10., 10);
       if(fabs(muonColl.at(i).dXY())<0.05){
         FillHist("IDflow_mu_3mu", 1., weight, 0., 10., 10);
         if(fabs(muonColl.at(i).dXYSig())<3.){
           FillHist("IDflow_mu_3mu", 2., weight, 0., 10., 10);
           if(RochIso04<0.1) FillHist("IDflow_mu_3mu", 3., weight, 0., 10., 10); 
         }
       }
       FillHist("Absd0_mu_3mu", fabs(muonColl.at(i).dXY()), weight, 0., 0.2, 2000);
       FillHist("Absd0sig_mu_3mu", fabs(muonColl.at(i).dXYSig()), weight, 0., 10., 100);
       FillHist("RelIso04_mu_3mu", RochIso04, weight, 0., 0.2, 200);

     }
   }

   if(DiMuon_analysis){
     if( !(PassTrigger("HLT_IsoMu24_v")||PassTrigger("HLT_IsoTkMu24_v")) ) return;
     weight *= trigger_sf;
     //if( !(PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")) ) return;


     //Step1 : 2OSmu +PTCut +Zwindow
     //if( !(muonColl.at(0).Pt()>27) ) return; //SingleMuon Trig Case
     if( !(muonColl.size()==2) ) return;
     //if( !(muonColl.at(0).Pt()>20 && muonColl.at(1).Pt()>10) ) return;
     if( !(muonColl.at(0).Pt()>27) ) return;
     if( fabs(SumCharge(muonColl))!=0 ) return;
     if( fabs((muonColl.at(0)+muonColl.at(1)).M()-91.2)>15 )   return;


     //Muon ID Property
     float RochIso04=999.;
     for(unsigned int i=0; i<muonColl.size(); i++){
       RochIso04=muonColl.at(i).RelIso04()/muonColl.at(i).RochPt()*muonColl.at(i).Pt();

       FillHist("Absd0_mu", fabs(muonColl.at(i).dXY()), weight, 0., 0.2, 2000);
       FillHist("Absdz_mu", fabs(muonColl.at(i).dZ()), weight, 0., 0.2, 2000);
       FillHist("Absd0sig_mu", fabs(muonColl.at(i).dXYSig()), weight, 0., 10., 100);
       FillHist("RelIso04_mu", RochIso04, weight, 0., 0.15, 200);

       float mupt=muonColl.at(i).Pt();
       if     (mupt<5 )  continue;//For safety
       else if(mupt<10){
         FillHist("Absd0_mu_5_10", fabs(muonColl.at(i).dXY()), weight, 0., 0.2, 2000);
         FillHist("Absdz_mu_5_10", fabs(muonColl.at(i).dZ()), weight, 0., 0.2, 2000);
         FillHist("Absd0sig_mu_5_10", fabs(muonColl.at(i).dXYSig()), weight, 0., 10., 100);
         FillHist("RelIso04_mu_5_10", RochIso04, weight, 0., 0.15, 200);
       }
       else if(mupt<20){
         FillHist("Absd0_mu_10_20", fabs(muonColl.at(i).dXY()), weight, 0., 0.2, 2000);
         FillHist("Absdz_mu_10_20", fabs(muonColl.at(i).dZ()), weight, 0., 0.2, 2000);
         FillHist("Absd0sig_mu_10_20", fabs(muonColl.at(i).dXYSig()), weight, 0., 10., 100);
         FillHist("RelIso04_mu_10_20", RochIso04, weight, 0., 0.15, 200);
       }
       else if(mupt<30){
         FillHist("Absd0_mu_20_30", fabs(muonColl.at(i).dXY()), weight, 0., 0.2, 2000);
         FillHist("Absdz_mu_20_30", fabs(muonColl.at(i).dZ()), weight, 0., 0.2, 2000);
         FillHist("Absd0sig_mu_20_30", fabs(muonColl.at(i).dXYSig()), weight, 0., 10., 100);
         FillHist("RelIso04_mu_20_30", RochIso04, weight, 0., 0.15, 200);
       }
       else if(mupt<50){
         FillHist("Absd0_mu_30_50", fabs(muonColl.at(i).dXY()), weight, 0., 0.2, 2000);
         FillHist("Absdz_mu_30_50", fabs(muonColl.at(i).dZ()), weight, 0., 0.2, 2000);
         FillHist("Absd0sig_mu_30_50", fabs(muonColl.at(i).dXYSig()), weight, 0., 10., 100);
         FillHist("RelIso04_mu_30_50", RochIso04, weight, 0., 0.15, 200);
       }
       else if(mupt<70){
         FillHist("Absd0_mu_50_70", fabs(muonColl.at(i).dXY()), weight, 0., 0.2, 2000);
         FillHist("Absdz_mu_50_70", fabs(muonColl.at(i).dZ()), weight, 0., 0.2, 2000);
         FillHist("Absd0sig_mu_50_70", fabs(muonColl.at(i).dXYSig()), weight, 0., 10., 100);
         FillHist("RelIso04_mu_50_70", RochIso04, weight, 0., 0.15, 200);
       }
       else if(mupt<100){
         FillHist("Absd0_mu_70_100", fabs(muonColl.at(i).dXY()), weight, 0., 0.2, 2000);
         FillHist("Absdz_mu_70_100", fabs(muonColl.at(i).dZ()), weight, 0., 0.2, 2000);
         FillHist("Absd0sig_mu_70_100", fabs(muonColl.at(i).dXYSig()), weight, 0., 10., 100);
         FillHist("RelIso04_mu_70_100", RochIso04, weight, 0., 0.15, 200);
       }
       else{
         FillHist("Absd0_mu_100_inf", fabs(muonColl.at(i).dXY()), weight, 0., 0.2, 2000);
         FillHist("Absdz_mu_100_inf", fabs(muonColl.at(i).dZ()), weight, 0., 0.2, 2000);
         FillHist("Absd0sig_mu_100_inf", fabs(muonColl.at(i).dXYSig()), weight, 0., 10., 100);
         FillHist("RelIso04_mu_100_inf", RochIso04, weight, 0., 0.15, 200);
       }

     }

   }
   if(DiEle_analysis){
     if( !(PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) ) return;


     //Step1 : 1e+2mu +PTCut
     //if( !(muonColl.at(0).Pt()>27) ) return; //SingleMuon Trig Case
     if( !(electronColl.size()==2) ) return;
     if( !(electronColl.at(0).Pt()>25 && electronColl.at(1).Pt()>15) ) return;
     if( electronColl.at(0).Charge() == electronColl.at(1).Charge() )  return;
     if( fabs((electronColl.at(0)+electronColl.at(1)).M()-91.2)>15 )   return;


     for(unsigned int i=0; i<electronColl.size(); i++){

       if(fabs(electronColl.at(i).Eta())<1.479){
         FillHist("Absd0_eB", fabs(electronColl.at(i).dxy()), weight, 0., 0.2, 2000);
         FillHist("Absdz_eB", fabs(electronColl.at(i).dz()), weight, 0., 0.2, 2000);
         FillHist("Absd0sig_eB", fabs(electronColl.at(i).dxySig()), weight, 0., 10., 100);
  
         float elpt=electronColl.at(i).Pt();
         if     (elpt<10 )  continue;//For safety
         else if(elpt<20){
           FillHist("Absd0_eB_10_20", fabs(electronColl.at(i).dxy()), weight, 0., 0.2, 2000);
           FillHist("Absdz_eB_10_20", fabs(electronColl.at(i).dz()), weight, 0., 0.2, 2000);
           FillHist("Absd0sig_eB_10_20", fabs(electronColl.at(i).dxySig()), weight, 0., 10., 100);
         }
         else if(elpt<30){
           FillHist("Absd0_eB_20_30", fabs(electronColl.at(i).dxy()), weight, 0., 0.2, 2000);
           FillHist("Absdz_eB_20_30", fabs(electronColl.at(i).dz()), weight, 0., 0.2, 2000);
           FillHist("Absd0sig_eB_20_30", fabs(electronColl.at(i).dxySig()), weight, 0., 10., 100);
         }
         else if(elpt<50){
           FillHist("Absd0_eB_30_50", fabs(electronColl.at(i).dxy()), weight, 0., 0.2, 2000);
           FillHist("Absdz_eB_30_50", fabs(electronColl.at(i).dz()), weight, 0., 0.2, 2000);
           FillHist("Absd0sig_eB_30_50", fabs(electronColl.at(i).dxySig()), weight, 0., 10., 100);
         }
         else if(elpt<70){
           FillHist("Absd0_eB_50_70", fabs(electronColl.at(i).dxy()), weight, 0., 0.2, 2000);
           FillHist("Absdz_eB_50_70", fabs(electronColl.at(i).dz()), weight, 0., 0.2, 2000);
           FillHist("Absd0sig_eB_50_70", fabs(electronColl.at(i).dxySig()), weight, 0., 10., 100);
         }
         else if(elpt<100){
           FillHist("Absd0_eB_70_100", fabs(electronColl.at(i).dxy()), weight, 0., 0.2, 2000);
           FillHist("Absdz_eB_70_100", fabs(electronColl.at(i).dz()), weight, 0., 0.2, 2000);
           FillHist("Absd0sig_eB_70_100", fabs(electronColl.at(i).dxySig()), weight, 0., 10., 100);
         }
         else{
           FillHist("Absd0_eB_100_inf", fabs(electronColl.at(i).dxy()), weight, 0., 0.2, 2000);
           FillHist("Absdz_eB_100_inf", fabs(electronColl.at(i).dz()), weight, 0., 0.2, 2000);
           FillHist("Absd0sig_eB_100_inf", fabs(electronColl.at(i).dxySig()), weight, 0., 10., 100);
         }
       }//BarrelEnds
       else{
         FillHist("Absd0_eE", fabs(electronColl.at(i).dxy()), weight, 0., 0.2, 2000);
         FillHist("Absdz_eE", fabs(electronColl.at(i).dz()), weight, 0., 0.2, 2000);
         FillHist("Absd0sig_eE", fabs(electronColl.at(i).dxySig()), weight, 0., 10., 100);
  
         float elpt=electronColl.at(i).Pt();
         if     (elpt<10 )  continue;//For safety
         else if(elpt<20){
           FillHist("Absd0_eE_10_20", fabs(electronColl.at(i).dxy()), weight, 0., 0.2, 2000);
           FillHist("Absdz_eE_10_20", fabs(electronColl.at(i).dz()), weight, 0., 0.2, 2000);
           FillHist("Absd0sig_eE_10_20", fabs(electronColl.at(i).dxySig()), weight, 0., 10., 100);
         }
         else if(elpt<30){
           FillHist("Absd0_eE_20_30", fabs(electronColl.at(i).dxy()), weight, 0., 0.2, 2000);
           FillHist("Absdz_eE_20_30", fabs(electronColl.at(i).dz()), weight, 0., 0.2, 2000);
           FillHist("Absd0sig_eE_20_30", fabs(electronColl.at(i).dxySig()), weight, 0., 10., 100);
         }
         else if(elpt<50){
           FillHist("Absd0_eE_30_50", fabs(electronColl.at(i).dxy()), weight, 0., 0.2, 2000);
           FillHist("Absdz_eE_30_50", fabs(electronColl.at(i).dz()), weight, 0., 0.2, 2000);
           FillHist("Absd0sig_eE_30_50", fabs(electronColl.at(i).dxySig()), weight, 0., 10., 100);
         }
         else if(elpt<70){
           FillHist("Absd0_eE_50_70", fabs(electronColl.at(i).dxy()), weight, 0., 0.2, 2000);
           FillHist("Absdz_eE_50_70", fabs(electronColl.at(i).dz()), weight, 0., 0.2, 2000);
           FillHist("Absd0sig_eE_50_70", fabs(electronColl.at(i).dxySig()), weight, 0., 10., 100);
         }
         else if(elpt<100){
           FillHist("Absd0_eE_70_100", fabs(electronColl.at(i).dxy()), weight, 0., 0.2, 2000);
           FillHist("Absdz_eE_70_100", fabs(electronColl.at(i).dz()), weight, 0., 0.2, 2000);
           FillHist("Absd0sig_eE_70_100", fabs(electronColl.at(i).dxySig()), weight, 0., 10., 100);
         }
         else{
           FillHist("Absd0_eE_100_inf", fabs(electronColl.at(i).dxy()), weight, 0., 0.2, 2000);
           FillHist("Absdz_eE_100_inf", fabs(electronColl.at(i).dz()), weight, 0., 0.2, 2000);
           FillHist("Absd0sig_eE_100_inf", fabs(electronColl.at(i).dxySig()), weight, 0., 10., 100);
         }
       }//Endcap Ends

     }//Ele Loop Ends


   }



/////////////////////////////////////////////////////////////////////////////////// 


return;
}// End of execute event loop
  


void Mar2017_3l4j_IDCompatibility::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Mar2017_3l4j_IDCompatibility::BeginCycle() throw( LQError ){
  
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

Mar2017_3l4j_IDCompatibility::~Mar2017_3l4j_IDCompatibility() {
  
  Message("In Mar2017_3l4j_IDCompatibility Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}

void Mar2017_3l4j_IDCompatibility::FillCutFlow(TString cut, float weight){
  
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



void Mar2017_3l4j_IDCompatibility::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Mar2017_3l4j_IDCompatibility::MakeHistograms(){
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
  // *  Remove//Overide this Mar2017_3l4j_IDCompatibilityCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Mar2017_3l4j_IDCompatibility::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
