#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TProfile.h>
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include "TLorentzVector.h"         // 4-vector class
#include "TCanvas.h"
#include "TH1F.h"
#endif

void test_sendToIlya(){

  TFile *f = new TFile("electron_ntuple.root");
  TTree *electronTree_ = (TTree*) f->FindObjectAny("ElectronTree");

  // Vars for PVs
  Int_t pvNTracks_;

  // Vars for pile-up
  Int_t nPUTrue_;    // true pile-up
  Int_t nPU_;        // generated pile-up
  Int_t nPV_;        // number of reconsrtucted primary vertices
  Float_t rho_;      // the rho variable

  // all electron variables
  Int_t nElectrons_;

  Float_t pt_;
  Float_t etaSC_;
  Float_t phiSC_;
  Float_t dEtaIn_;
  Float_t dPhiIn_;
  Float_t hOverE_;
  Float_t full5x5_sigmaIetaIeta_;
  Float_t isoChargedHadrons_;
  Float_t isoNeutralHadrons_;
  Float_t isoPhotons_;
  Float_t isoChargedFromPU_;
  Float_t ooEmooP_;
  Float_t d0_;
  Float_t dz_;
  Int_t expectedMissingInnerHits_;
  Int_t passConversionVeto_;     
  Int_t passVetoId_;
  Int_t passLooseId_;
  Int_t passMediumId_;
  Int_t passTightId_;
  Int_t isTrue_;
  
  electronTree_->SetBranchAddress("pvNTracks"    ,  &pvNTracks_);

  electronTree_->SetBranchAddress("nPV"        ,  &nPV_);
  electronTree_->SetBranchAddress("nPU"        ,  &nPU_);
  electronTree_->SetBranchAddress("nPUTrue"    ,  &nPUTrue_);
  electronTree_->SetBranchAddress("rho"        ,  &rho_);

  //electronTree_->SetBranchAddress("nEle"    ,  &nElectrons_);
  electronTree_->SetBranchAddress("pt"    ,  &pt_    );
  electronTree_->SetBranchAddress("eta" ,    &etaSC_ );
  electronTree_->SetBranchAddress("phi" ,    &phiSC_ );
  electronTree_->SetBranchAddress("dEtaIn",  &dEtaIn_);
  electronTree_->SetBranchAddress("dPhiIn",  &dPhiIn_);
  electronTree_->SetBranchAddress("hOverE",  &hOverE_);
  electronTree_->SetBranchAddress("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta_);
  electronTree_->SetBranchAddress("isoChargedHadrons"      , &isoChargedHadrons_);
  electronTree_->SetBranchAddress("isoNeutralHadrons"      , &isoNeutralHadrons_);
  electronTree_->SetBranchAddress("isoPhotons"             , &isoPhotons_);
  electronTree_->SetBranchAddress("isoChargedFromPU"       , &isoChargedFromPU_);
  electronTree_->SetBranchAddress("ooEmooP", &ooEmooP_);
  electronTree_->SetBranchAddress("d0"     , &d0_);
  electronTree_->SetBranchAddress("dz"     , &dz_);
  electronTree_->SetBranchAddress("expectedMissingInnerHits", &expectedMissingInnerHits_);
  electronTree_->SetBranchAddress("passConversionVeto", &passConversionVeto_);
  electronTree_->SetBranchAddress("passVetoId"  ,  &passVetoId_ );
  electronTree_->SetBranchAddress("passLooseId"  ,  &passLooseId_ );
  electronTree_->SetBranchAddress("passMediumId" ,  &passMediumId_ );
  electronTree_->SetBranchAddress("passTightId"  ,  &passTightId_ );

  electronTree_->SetBranchAddress("isTrue"    , &isTrue_);

  Int_t p=0;
  TProfile *profID_pt = new TProfile("profID_pt", "profID_pt", 9,10,100);
  TProfile *profCuts_pt = new TProfile("profCuts_pt", "profCuts_pt", 9,10,100);
  TProfile *profID_eta = new TProfile("profID_eta", "profID_eta", 10,-2.5,2.5);
  TProfile *profCuts_eta = new TProfile("profCuts_eta", "profCuts_eta", 10,-2.5,2.5);

  int IDisone = 0;
  int IDiszero = 0;
  for (Int_t i=0; i<electronTree_->GetEntries(); i++) {
    //for (Int_t i=0; i<10; i++) {
    electronTree_->GetEntry(i);

    int passVetoId_cuts = 1;

    double ea;
    if(0.0000 <= fabs(etaSC_) && fabs(etaSC_) < 0.8000) ea = 0.1013;
    else if(0.8000 <= fabs(etaSC_) && fabs(etaSC_) < 1.3000) ea = 0.0988;
    else if(1.3000 <= fabs(etaSC_) && fabs(etaSC_) < 2.0000) ea = 0.0572;
    else if(2.0000 <= fabs(etaSC_) && fabs(etaSC_) < 2.2000) ea = 0.0842;
    else if(2.2000 <= fabs(etaSC_) && fabs(etaSC_) < 5.0000) ea = 0.1530;

    double relIsoWithEA = 1/pt_*(isoChargedHadrons_+TMath::Max(0.,isoNeutralHadrons_+isoPhotons_-rho_*ea));

    // barrel
    if(fabs(etaSC_)<=1.479) {
      if(dEtaIn_ >= 0.013625) passVetoId_cuts = 0;
      if(dPhiIn_ >= 0.230374) passVetoId_cuts = 0;
      if(full5x5_sigmaIetaIeta_ >= 0.011586) passVetoId_cuts = 0;
      if(hOverE_ >= 0.181130) passVetoId_cuts = 0;
      if(d0_ >= 0.094095) passVetoId_cuts = 0;
      if(dz_ >= 0.713070) passVetoId_cuts = 0;
      if(ooEmooP_ >= 0.295751) passVetoId_cuts = 0;
      if(relIsoWithEA >= 0.158721) passVetoId_cuts = 0;
      if(expectedMissingInnerHits_>2) passVetoId_cuts = 0;
      if(!passConversionVeto_) passVetoId_cuts = 0;
    } else if(1.479 < fabs(etaSC_) && fabs(etaSC_) < 2.5) { // endcap
      if(dEtaIn_ >= 0.011932) passVetoId_cuts = 0;
      if(dPhiIn_ >= 0.255450) passVetoId_cuts = 0;
      if(full5x5_sigmaIetaIeta_ >= 0.031849) passVetoId_cuts = 0;
      if(hOverE_ >= 0.223870) passVetoId_cuts = 0;
      if(d0_ >= 0.342293) passVetoId_cuts = 0;
      if(dz_ >= 0.953461) passVetoId_cuts = 0;
      if(ooEmooP_ >= 0.155501) passVetoId_cuts = 0;
      if(relIsoWithEA >= 0.177032) passVetoId_cuts = 0;
      if(expectedMissingInnerHits_>3) passVetoId_cuts = 0;
      if(!passConversionVeto_) passVetoId_cuts = 0;
    }

    if(passVetoId_==0 && passVetoId_cuts==1) IDiszero++;
    if(passVetoId_==1 && passVetoId_cuts==0) IDisone++;

    profID_pt->Fill(pt_, passVetoId_);
    profID_eta->Fill(etaSC_, passVetoId_);
    profCuts_pt->Fill(pt_,passVetoId_cuts);
    profCuts_eta->Fill(etaSC_,passVetoId_cuts);

    if (passVetoId_) p++;
  }

  cout << ((double)p)/electronTree_->GetEntries() << endl;
  cout << "total electrons = " << electronTree_->GetEntries() << endl;
  cout << "IDiszero = " << IDiszero << " ---> " << IDiszero/((double)electronTree_->GetEntries()) << endl;
  cout << "IDisone = " << IDisone << " ---> " << IDisone/((double)electronTree_->GetEntries()) << endl;

  TCanvas* cprofID_pt = new TCanvas("cprofID_pt","cprofID_pt");
  profID_pt->SetLineColor(kRed); profCuts_pt->SetLineColor(kBlue);
  profID_pt->GetXaxis()->SetTitle("pT");
  profID_pt->GetYaxis()->SetTitle("Efficiency");
  profID_pt->Draw(); profCuts_pt->Draw("same"); gPad->BuildLegend();
  cprofID_pt->Print("cprofID_pt.png");

  TCanvas* cprofID_eta = new TCanvas("cprofID_eta","cprofID_eta");
  profID_eta->SetLineColor(kRed); profCuts_eta->SetLineColor(kBlue);
  profID_eta->GetXaxis()->SetTitle("eta");
  profID_eta->GetYaxis()->SetTitle("Efficiency");
  profID_eta->Draw(); profCuts_eta->Draw("same"); gPad->BuildLegend();
  cprofID_eta->Print("cprofID_eta.png");

}
