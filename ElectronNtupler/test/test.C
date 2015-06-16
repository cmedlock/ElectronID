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

void test(){

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
  TProfile *prof = new TProfile("prof", "prof", 9,10,100);
  TProfile *prof_check = new TProfile("prof_check", "prof_check", 9,10,100);

  TH1F* hrelIsoWithEA = new TH1F("hrelIsoWithEA","hrelIsoWithEA",50,0,0.2);
  TH1F* hrelIsoWithEA_agree = new TH1F("hrelIsoWithEA_agree","hrelIsoWithEA_agree",50,0,0.2);
  TH1F* hdEtaIn = new TH1F("hdEtaIn","hdEtaIn",50,-0.2,0.2);
  TH1F* hdEtaIn_agree = new TH1F("hdEtaIn_agree","hdEtaIn_agree",50,-0.2,0.2);
  TH1F* hdPhiIn = new TH1F("hdPhiIn","hdPhiIn",50,-0.3,0.3);
  TH1F* hdPhiIn_agree = new TH1F("hdPhiIn_agree","hdPhiIn_agree",50,-0.3,0.3);
  TH1F* hsieie = new TH1F("hsieie","hsieie",50,0,0.1);
  TH1F* hsieie_agree = new TH1F("hsieie_agree","hsieie_agree",50,0,0.1);
  TH1F* hhOverE = new TH1F("hhOverE","hhOverE",50,0,0.4);
  TH1F* hhOverE_agree = new TH1F("hhOverE_agree","hhOverE_agree",50,0,0.4);
  TH1F* hd0 = new TH1F("hd0","hd0",50,-0.4,0.4);
  TH1F* hd0_agree = new TH1F("hd0_agree","hd0_agree",50,-0.4,0.4);
  TH1F* hdz = new TH1F("hdz","hdz",50,-0.4,0.4);
  TH1F* hdz_agree = new TH1F("hdz_agree","hdz_agree",50,-0.4,0.4);
  TH1F* hooEmooP = new TH1F("hooEmooP","hooEmooP",50,0,0.3);
  TH1F* hooEmooP_agree = new TH1F("hooEmooP_agree","hooEmooP_agree",50,0,0.3);
  TH1F* hmissingHits = new TH1F("hmissingHits","hmissingHits",3,0,3);
  TH1F* hmissingHits_agree = new TH1F("hmissingHits_agree","hmissingHits_agree",3,0,3);
  TH1F* hconversion = new TH1F("hconversion","hconversion",3,0,3);
  TH1F* hconversion_agree = new TH1F("hconversion_agree","hconversion_agree",3,0,3);

  int trueisone = 0;
  int trueiszero = 0;
  for (Int_t i=0; i<electronTree_->GetEntries(); i++) {
    //for (Int_t i=0; i<10; i++) {
    electronTree_->GetEntry(i);

    int passVetoId_check = 1;

    double ea;
    if(0.0000 <= fabs(etaSC_) && fabs(etaSC_) < 0.8000) ea = 0.1013;
    else if(0.8000 <= fabs(etaSC_) && fabs(etaSC_) < 1.3000) ea = 0.0988;
    else if(1.3000 <= fabs(etaSC_) && fabs(etaSC_) < 2.0000) ea = 0.0572;
    else if(2.0000 <= fabs(etaSC_) && fabs(etaSC_) < 2.2000) ea = 0.0842;
    else if(2.2000 <= fabs(etaSC_) && fabs(etaSC_) < 5.0000) ea = 0.1530;

    double relIsoWithEA = 1/pt_*(isoChargedHadrons_+TMath::Max(0.,isoNeutralHadrons_+isoPhotons_-rho_*ea));

    // barrel
    if(fabs(etaSC_)<=1.479) {
      if(dEtaIn_ >= 0.013625) passVetoId_check = 0;
      if(dPhiIn_ >= 0.230374) passVetoId_check = 0;
      if(full5x5_sigmaIetaIeta_ >= 0.011586) passVetoId_check = 0;
      if(hOverE_ >= 0.181130) passVetoId_check = 0;
      if(d0_ >= 0.094095) passVetoId_check = 0;
      if(dz_ >= 0.713070) passVetoId_check = 0;
      if(ooEmooP_ >= 0.295751) passVetoId_check = 0;
      if(relIsoWithEA >= 0.158721) passVetoId_check = 0;
      if(expectedMissingInnerHits_>2) passVetoId_check = 0;
      if(!passConversionVeto_) passVetoId_check = 0;
    } else if(1.479 < fabs(etaSC_) && fabs(etaSC_) < 2.5) { // endcap
      if(dEtaIn_ >= 0.011932) passVetoId_check = 0;
      if(dPhiIn_ >= 0.255450) passVetoId_check = 0;
      if(full5x5_sigmaIetaIeta_ >= 0.031849) passVetoId_check = 0;
      if(hOverE_ >= 0.223870) passVetoId_check = 0;
      if(d0_ >= 0.342293) passVetoId_check = 0;
      if(dz_ >= 0.953461) passVetoId_check = 0;
      if(ooEmooP_ >= 0.155501) passVetoId_check = 0;
      if(relIsoWithEA >= 0.177032) passVetoId_check = 0;
      if(expectedMissingInnerHits_>3) passVetoId_check = 0;
      if(!passConversionVeto_) passVetoId_check = 0;
    }

    if(passVetoId_==0 && passVetoId_check==1) trueiszero++;
    if(passVetoId_==1 && passVetoId_check==0) trueisone++;

    if(passVetoId_==passVetoId_check) {
      hrelIsoWithEA_agree->Fill(relIsoWithEA);
      hdEtaIn_agree->Fill(dEtaIn_);
      hdPhiIn_agree->Fill(dPhiIn_);
      hsieie_agree->Fill(full5x5_sigmaIetaIeta_);
      hhOverE_agree->Fill(hOverE_);
      hd0_agree->Fill(d0_);
      hdz_agree->Fill(dz_);
      hooEmooP_agree->Fill(ooEmooP_);
      hmissingHits_agree->Fill(expectedMissingInnerHits_);
      hconversion_agree->Fill(passConversionVeto_);
    } else if(passVetoId_==0 && passVetoId_check==1) {
      hrelIsoWithEA->Fill(relIsoWithEA);
      hdEtaIn->Fill(dEtaIn_);
      hdPhiIn->Fill(dPhiIn_);
      hsieie->Fill(full5x5_sigmaIetaIeta_);
      hhOverE->Fill(hOverE_);
      hd0->Fill(d0_);
      hdz->Fill(dz_);
      hooEmooP->Fill(ooEmooP_);
      hmissingHits->Fill(expectedMissingInnerHits_);
      hconversion->Fill(passConversionVeto_);
    }

    prof->Fill(pt_, passVetoId_);
    prof_check->Fill(pt_,passVetoId_check);

    if (passVetoId_) p++;
  }

  // scale
  hrelIsoWithEA->Scale(1/hrelIsoWithEA->Integral());
  hrelIsoWithEA_agree->Scale(1/hrelIsoWithEA_agree->Integral());
  hdEtaIn->Scale(1/hdEtaIn->Integral());
  hdEtaIn_agree->Scale(1/hdEtaIn_agree->Integral());
  hdPhiIn->Scale(1/hdPhiIn->Integral());
  hdPhiIn_agree->Scale(1/hdPhiIn_agree->Integral());
  hsieie->Scale(1/hsieie->Integral());
  hsieie_agree->Scale(1/hsieie_agree->Integral());
  hhOverE->Scale(1/hhOverE->Integral());
  hhOverE_agree->Scale(1/hhOverE_agree->Integral());
  hd0->Scale(1/hd0->Integral());
  hd0_agree->Scale(1/hd0_agree->Integral());
  hdz->Scale(1/hdz->Integral());
  hdz_agree->Scale(1/hdz_agree->Integral());
  hooEmooP->Scale(1/hooEmooP->Integral());
  hooEmooP_agree->Scale(1/hooEmooP_agree->Integral());
  hmissingHits->Scale(1/hmissingHits->Integral());
  hmissingHits_agree->Scale(1/hmissingHits_agree->Integral());
  hconversion->Scale(1/hconversion->Integral());
  hconversion_agree->Scale(1/hconversion_agree->Integral());
  
  cout << ((double)p)/electronTree_->GetEntries() << endl;
  cout << "total electrons = " << electronTree_->GetEntries() << endl;
  cout << "trueiszero = " << trueiszero << " ---> " << trueiszero/((double)electronTree_->GetEntries()) << endl;
  cout << "trueisone = " << trueisone << " ---> " << trueisone/((double)electronTree_->GetEntries()) << endl;
/*
  TCanvas* cprof = new TCanvas();
  prof->Draw(); prof_check->Draw("same"); gPad->BuildLegend();
  cprof->Print("cprof.png");

  TCanvas* crelIsoWithEA = new TCanvas();
  hrelIsoWithEA->Draw();  hrelIsoWithEA->SetLineColor(kRed); hrelIsoWithEA_agree->Draw("same");
  crelIsoWithEA->SetLogy(); crelIsoWithEA->Print("crelIsoWithEA.png");

  TCanvas* cdEtaIn = new TCanvas();
  hdEtaIn_agree->Draw(); hdEtaIn->Draw("same"); gPad->BuildLegend();
  cdEtaIn->SetLogy(); cdEtaIn->Print("cdEtaIn.png");

  TCanvas* cdPhiIn = new TCanvas();
  hdPhiIn_agree->Draw(); hdPhiIn->Draw("same"); gPad->BuildLegend();
  cdPhiIn->SetLogy(); cdPhiIn->Print("cdPhiIn.png");

  TCanvas* csieie = new TCanvas();
  hsieie_agree->Draw(); hsieie->Draw("same"); gPad->BuildLegend();
  csieie->SetLogy(); csieie->Print("csieie.png");

  TCanvas* chOverE = new TCanvas();
  hhOverE_agree->Draw(); hhOverE->Draw("same"); gPad->BuildLegend();
  chOverE->SetLogy(); chOverE->Print("chOverE.png");

  TCanvas* cd0 = new TCanvas();
  hd0_agree->Draw(); hd0->Draw("same"); gPad->BuildLegend();
  cd0->SetLogy(); cd0->Print("cd0.png");

  TCanvas* cdz = new TCanvas();
  hdz_agree->Draw(); hdz->Draw("same"); gPad->BuildLegend();
  cdz->SetLogy(); cdz->Print("cdz.png");

  TCanvas* cooEmooP = new TCanvas();
  hooEmooP_agree->Draw(); hooEmooP->Draw("same"); gPad->BuildLegend();
  cooEmooP->SetLogy(); cooEmooP->Print("cooEmooP.png");

  TCanvas* cmissingHits = new TCanvas();
  hmissingHits_agree->Draw(); hmissingHits->Draw("same"); gPad->BuildLegend();
  cmissingHits->SetLogy(); cmissingHits->Print("cmissingHits.png");

  TCanvas* cconversion = new TCanvas();
  hconversion_agree->Draw(); hconversion->Draw("same"); gPad->BuildLegend();
  cconversion->SetLogy(); cconversion->Print("cconversion.png");
*/
}
