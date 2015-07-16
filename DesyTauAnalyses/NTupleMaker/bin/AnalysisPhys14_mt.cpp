/////////////////////////////////////////////////////////////
// Analysis Macro for Synch Phys14 Ntuple for h->tau tau
// Author: Francesco Costanza <francesco.costanza@desy.de>
//
// Wed Jul 15 11:50:56 2015 by ROOT version 5.34/18
/////////////////////////////////////////////////////////////

#include <iostream>

#include "TFile.h"
#include "TChain.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Phys14Tree.h"

const float MuMass = 0.105658367;

int main(int argc, char * argv[]) {

  // first argument - filename

  using namespace std;

  TChain* t = new TChain("TauCheck");

  if(argc < 1) return -1;

  TFile* f = new TFile(argv[1], "read");
  if(!f->FindObjectAny("TauCheck")) return -2;
  f->Close();
  
  t->Add(argv[1]);
  if(!t) return -3;
  
  Phys14Tree* phys14 = new Phys14Tree(t);

  //Initialiazing leafs for convenience:
  Int_t           *run = &phys14->run;
  Int_t           *lumi = &phys14->lumi;
  Int_t           *evt = &phys14->evt;
  Int_t           *npv = &phys14->npv;
  Int_t           *npu = &phys14->npu;
  Float_t         *rho = &phys14->rho;
  Float_t         *mcweight = &phys14->mcweight;
  Float_t         *puweight = &phys14->puweight;
  Float_t         *trigweight_1 = &phys14->trigweight_1;
  Float_t         *trigweight_2 = &phys14->trigweight_2;
  Float_t         *idweight_1 = &phys14->idweight_1;
  Float_t         *idweight_2 = &phys14->idweight_2;
  Float_t         *isoweight_1 = &phys14->isoweight_1;
  Float_t         *isoweight_2 = &phys14->isoweight_2;
  Float_t         *effweight = &phys14->effweight;
  Float_t         *fakeweight = &phys14->fakeweight;
  Float_t         *embeddedWeight = &phys14->embeddedWeight;
  Float_t         *signalWeight = &phys14->signalWeight;
  Float_t         *weight = &phys14->weight;
  Float_t         *m_vis = &phys14->m_vis;
  Float_t         *m_sv = &phys14->m_sv;
  Float_t         *pt_sv = &phys14->pt_sv;
  Float_t         *eta_sv = &phys14->eta_sv;
  Float_t         *phi_sv = &phys14->phi_sv;
  Float_t         *pt_1 = &phys14->pt_1;
  Float_t         *phi_1 = &phys14->phi_1;
  Float_t         *eta_1 = &phys14->eta_1;
  Float_t         *m_1 = &phys14->m_1;
  Int_t           *q_1 = &phys14->q_1;
  Float_t         *iso_1 = &phys14->iso_1;
  Float_t         *mva_1 = &phys14->mva_1;
  Float_t         *d0_1 = &phys14->d0_1;
  Float_t         *dZ_1 = &phys14->dZ_1;
  Float_t         *mt_1 = &phys14->mt_1;
  Float_t         *pt_2 = &phys14->pt_2;
  Float_t         *phi_2 = &phys14->phi_2;
  Float_t         *eta_2 = &phys14->eta_2;
  Float_t         *m_2 = &phys14->m_2;
  Int_t           *q_2 = &phys14->q_2;
  Float_t         *iso_2 = &phys14->iso_2;
  Float_t         *d0_2 = &phys14->d0_2;
  Float_t         *dZ_2 = &phys14->dZ_2;
  Float_t         *mva_2 = &phys14->mva_2;
  Float_t         *mt_2 = &phys14->mt_2;
  Char_t          *os = &phys14->os;
  Char_t          *dilepton_veto = &phys14->dilepton_veto;
  Char_t          *extraelec_veto = &phys14->extraelec_veto;
  Char_t          *extramuon_veto = &phys14->extramuon_veto;
  Float_t         *byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = &phys14->byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
  Float_t         *againstElectronLooseMVA5_1 = &phys14->againstElectronLooseMVA5_1;
  Float_t         *againstElectronMediumMVA5_1 = &phys14->againstElectronMediumMVA5_1;
  Float_t         *againstElectronTightMVA5_1 = &phys14->againstElectronTightMVA5_1;
  Float_t         *againstElectronVLooseMVA5_1 = &phys14->againstElectronVLooseMVA5_1;
  Float_t         *againstElectronVTightMVA5_1 = &phys14->againstElectronVTightMVA5_1;
  Float_t         *againstMuonLoose3_1 = &phys14->againstMuonLoose3_1;
  Float_t         *againstMuonTight3_1 = &phys14->againstMuonTight3_1;
  Float_t         *byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = &phys14->byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
  Float_t         *againstElectronLooseMVA5_2 = &phys14->againstElectronLooseMVA5_2;
  Float_t         *againstElectronMediumMVA5_2 = &phys14->againstElectronMediumMVA5_2;
  Float_t         *againstElectronTightMVA5_2 = &phys14->againstElectronTightMVA5_2;
  Float_t         *againstElectronVLooseMVA5_2 = &phys14->againstElectronVLooseMVA5_2;
  Float_t         *againstElectronVTightMVA5_2 = &phys14->againstElectronVTightMVA5_2;
  Float_t         *againstMuonLoose3_2 = &phys14->againstMuonLoose3_2;
  Float_t         *againstMuonTight3_2 = &phys14->againstMuonTight3_2;
  Float_t         *met = &phys14->met;
  Float_t         *metphi = &phys14->metphi;
  Float_t         *metcov00 = &phys14->metcov00;
  Float_t         *metcov01 = &phys14->metcov01;
  Float_t         *metcov10 = &phys14->metcov10;
  Float_t         *metcov11 = &phys14->metcov11;
  Float_t         *mvamet = &phys14->mvamet;
  Float_t         *mvametphi = &phys14->mvametphi;
  Float_t         *mvacov00 = &phys14->mvacov00;
  Float_t         *mvacov01 = &phys14->mvacov01;
  Float_t         *mvacov10 = &phys14->mvacov10;
  Float_t         *mvacov11 = &phys14->mvacov11;
  Float_t         *pt_tt = &phys14->pt_tt;
  Float_t         *pzetavis = &phys14->pzetavis;
  Float_t         *pzetamiss = &phys14->pzetamiss;
  Float_t         *mva_gf = &phys14->mva_gf;
  Int_t           *njets = &phys14->njets;
  Int_t           *njetspt20 = &phys14->njetspt20;
  Float_t         *jpt_1 = &phys14->jpt_1;
  Float_t         *jeta_1 = &phys14->jeta_1;
  Float_t         *jphi_1 = &phys14->jphi_1;
  Float_t         *jptraw_1 = &phys14->jptraw_1;
  Float_t         *jptunc_1 = &phys14->jptunc_1;
  Float_t         *jmva_1 = &phys14->jmva_1;
  Float_t         *jlrm_1 = &phys14->jlrm_1;
  Int_t           *jctm_1 = &phys14->jctm_1;
  Float_t         *jpt_2 = &phys14->jpt_2;
  Float_t         *jeta_2 = &phys14->jeta_2;
  Float_t         *jphi_2 = &phys14->jphi_2;
  Float_t         *jptraw_2 = &phys14->jptraw_2;
  Float_t         *jptunc_2 = &phys14->jptunc_2;
  Float_t         *jmva_2 = &phys14->jmva_2;
  Float_t         *jlrm_2 = &phys14->jlrm_2;
  Int_t           *jctm_2 = &phys14->jctm_2;
  Float_t         *mjj = &phys14->mjj;
  Float_t         *jdeta = &phys14->jdeta;
  Int_t           *njetingap = &phys14->njetingap;
  Int_t           *nbtag = &phys14->nbtag;
  Float_t         *bpt = &phys14->bpt;
  Float_t         *beta = &phys14->beta;
  Float_t         *bphi = &phys14->bphi;  

  std::cout<<"Starting loop.."<<std::endl;
  std::cout<<"Ntries = "<<phys14->GetEntries()<<std::endl;  
  
  while(phys14->GetEntry(-1) > 0)
    cout<<"Muon pt = "<<*pt_1<<"GeV"<<std::endl;

  return 0;
}
