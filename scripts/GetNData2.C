#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include <iomanip>
#include <sstream>
#include "TSystem.h"

#include <map>

float GetEventsProcessed(std::string file);
float GetEventsPassed(std::string file);
float GetSumWeights(std::string filename);

map<TString, Double_t> map_lumi;
map<TString, Double_t> neventmap;
map<TString, Double_t> n_w_eventmap;
float GetNData2(TString SampleName) {
  

  TString def_version = TString(getenv("CATVERSION"));
  
  vector<int> v763; 
  vector<int> v765; 
    
  
  //TString dir = "ls /data7/DATA/CatNtuples/v8-0-7/SKTrees/DataTriLep/DoubleMuon/"+SampleName+"/*.root > inputlist_efflumi.txt";
  TString dir = "ls /data7/DATA/CatNtuples/v8-0-7/SKTrees/DataTriLep/MuonEG/"+SampleName+"/*.root > inputlist_efflumi.txt";
  //TString dir = "ls /data2/CatNtuples/v8-0-7/SKTrees/MCTriLep/"+SampleName+"/*.root > inputlist_efflumi.txt";
  
  bool use_sum_genweight(false);
  
  system(dir.Data());
  
  
  std::ifstream fin("inputlist_efflumi.txt");

  std::string word;
  
  float number_events_processed(0.);
  float number_events_processed2(0.);
  float number_events_passed(0.);
  float number_events_passed2(0.);
  float sum_of_weights(0.);

  while ( fin >> word ) {
    number_events_processed+= GetEventsProcessed(word);
    number_events_processed+= 1;
  }
  //cout<<"Total Nevents: "<<number_events_processed<<endl;
  return number_events_processed;
  
}


float GetEventsProcessed(std::string filename){
  TFile* file = TFile::Open(filename.c_str());
  //  cout << file << endl;
  TH1F*  EventCounter = (TH1F*) (file ->Get("cutflow"));
  
  float value = EventCounter->GetBinContent(4);
  file->Close();
  return value;
}

float GetEventsPassed(std::string filename){
  TFile* file = TFile::Open(filename.c_str());

  TH1F*  EventCounter = (TH1F*) (file ->Get("cutflow"));

  float value = EventCounter->GetBinContent(1);
  file->Close();
  return value;
}

float GetSumWeights(std::string filename){
  TFile* file = TFile::Open(filename.c_str());

  TDirectory * dir = (TDirectory*)file->Get(TString(filename) + ":/ntuple");
  TTree * tree;
  dir->GetObject("event",tree);
  
  float genWeight=0;
  TBranch        *b_genWeight; 
  
  tree->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
  
  float sum_weight=0.;
  Long64_t nentries = tree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    tree->LoadTree(jentry);
    nb = tree->GetEntry(jentry);
    sum_weight += genWeight;
  }
  cout << "sum_weight = " << sum_weight << endl;
  return sum_weight;
}
