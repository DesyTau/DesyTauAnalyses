/////////////////////////////////////////////////////////////
// Analysis Macro for Synch Phys14 Ntuple for h->tau tau
// Author: Francesco Costanza <francesco.costanza@desy.de>
//
// Wed Jul 15 11:50:56 2015 by ROOT version 5.34/18
/////////////////////////////////////////////////////////////

#include <iostream>
#include <utility>
#include <algorithm>

#include "TFile.h"
#include "TChain.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Phys14Tree.h"

struct sort_second {
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
        return left.second < right.second;
    }
};

void load_entries( Phys14Tree& t, std::vector< std::pair< int, int> >& entry){
  while(t.GetEntry(-1) > 0)
    entry.push_back(std::make_pair(t.LoadedEntryId(), t.evt));

  std::sort(entry.begin(), entry.end(), sort_second());
}


// main

int main(int argc, char * argv[]) {

  // first argument - ref file
  // second argument - test file
  
  using namespace std;

  TChain* tref = new TChain("TauCheck");
  TChain* ttest = new TChain("TauCheck");
  
  if(argc < 2) return -1;

  std::cout<<"ref: "<< argv[1]<<std::endl;
  std::cout<<"test: "<< argv[2]<<std::endl;  
  
  TFile* fref = new TFile(argv[1], "read");
  if(!fref->FindObjectAny("TauCheck")) return -2;
  fref->Close();
  
  tref->Add(argv[1]);
  if(!tref) return -3;

  TFile* ftest = new TFile(argv[2], "read");
  if(!ftest->FindObjectAny("TauCheck")) return -2;
  ftest->Close();
  
  ttest->Add(argv[2]);
  if(!ttest) return -3;

  Phys14Tree* ref = new Phys14Tree(tref);
  std::vector<std::pair<int, int> > ref_entry;  
  load_entries(*ref, ref_entry);
  
  Phys14Tree* test = new Phys14Tree(ttest);
  std::vector<std::pair<int, int> > test_entry;  
  load_entries(*test, test_entry);

  if(test_entry.size() != ref_entry.size())
    std::cout<<"Ntry mismatch: ref = "<<ref_entry.size()<<"; test = "<<test_entry.size()<<";"<<std::endl;

  UInt_t iref = 0;
  UInt_t itest = 0;

  int ref_id = 0;
  int test_id = 0;

  std::vector<std::pair<int, int> > ref_missing; 
  std::vector<std::pair<int, int> > test_missing; 

  std::vector<std::pair< std::pair< int, int>, std::pair< int, int> > > bad_events; 
  //std::vector<std::pair<int, std::pair< int, int> > > good;  

  
  for( ; iref < ref_entry.size(); iref++, itest++){
    ref_id = ref_entry.at(iref).second;
    
    while(test_entry.at(itest).second < ref_id){
      ref_missing.push_back(test_entry.at(itest));
      itest++;
    }
    test_id = test_entry.at(itest).second;
    
    if (ref_id != test_id){
      test_missing.push_back(ref_entry.at(iref));
      itest--;
      continue;
    }

    // lets compare!
    ref->GetEntry(ref_entry.at(iref).first);
    test->GetEntry(test_entry.at(itest).first);
      
    if( ref->pt_1 != test->pt_1){
      std::cout<<ref->pt_1<<" "<<test->pt_1<<std::endl;

      bad_events.push_back(std::make_pair( ref_entry.at(iref), ref_entry.at(itest)));
    }
  }
  std::cout<<std::endl;

  std::cout<<test_entry.size()-ref_missing.size()-bad_events.size()<<" events with perfect match."<<std::endl;
  std::cout<<ref_entry.size()-test_missing.size()-bad_events.size()<<" events with perfect match."<<std::endl;  
  std::cout<<std::endl;
  
  std::cout<<test_missing.size()<<" events missing in test tree:"<<std::endl;
  for( UInt_t i = 0; i < test_missing.size(); i++)
    std::cout<<test_missing.at(i).second<<", ";
  std::cout<<std::endl<<std::endl;

  std::cout<<ref_missing.size()<<" events missing in ref tree:"<<std::endl;
  for( UInt_t i = 0; i < ref_missing.size(); i++)
    std::cout<<ref_missing.at(i).second<<", ";
  std::cout<<std::endl<<std::endl;

  std::cout<<bad_events.size()<<" events show differences:"<<std::endl;
  for( UInt_t i = 0; i < bad_events.size(); i++)
    std::cout<<bad_events.at(i).second.second<<", ";
  std::cout<<std::endl<<std::endl;
    
  return 0;
}
