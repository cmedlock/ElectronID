// -*- C++ -*-
//
// Package:    EgammaWork/ElectronNtupler
// Class:      ElectronNtuplerVIDDemo
// 
/**\class ElectronNtuplerVIDDemo ElectronNtuplerVIDDemo.cc EgammaWork/ElectronNtupler/plugins/ElectronNtuplerVIDDemo.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ilya Kravchenko
//         Created:  Thu, 10 Jul 2014 09:54:13 GMT
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "Math/VectorUtil.h"

//
// class declaration
//

class ElectronNtuplerVIDDemo : public edm::EDAnalyzer {
   public:
      explicit ElectronNtuplerVIDDemo(const edm::ParameterSet&);
      ~ElectronNtuplerVIDDemo();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  enum ElectronMatchType {UNMATCHED = 0, 
			  TRUE_PROMPT_ELECTRON, 
			  TRUE_ELECTRON_FROM_TAU,
			  TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      int matchToTruth(const edm::Ptr<reco::GsfElectron> el, 
		       const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);

      void findFirstNonElectronMother(const reco::Candidate *particle,
				    int &ancestorPID, int &ancestorStatus);

      // ----------member data ---------------------------

      // Data members that are the same for AOD and miniAOD
      edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupToken_;
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;

      // AOD case data members
      edm::EDGetToken electronsToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;

      // MiniAOD case data members
      edm::EDGetToken electronsMiniAODToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxMiniAODToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsMiniAODToken_;

      // ID decisions objects
      edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;

  TTree *electronTree_;

  // Vars for PVs
  Int_t pvNTracks_;

  // Vars for pile-up
  Int_t nPUTrue_;    // true pile-up
  Int_t nPU_;        // generated pile-up
  Int_t nPV_;        // number of reconsrtucted primary vertices
  Float_t rho_;      // the rho variable

  // all variables for the output tree
  //Int_t nElectrons_;

  Float_t pt_;
  Float_t eta_;
  Float_t phi_;
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
  Int_t   expectedMissingInnerHits_;
  Int_t   passConversionVeto_;

  Int_t passVetoId_;
  Int_t passLooseId_;
  Int_t passMediumId_;
  Int_t passTightId_;

  Int_t isTrue_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ElectronNtuplerVIDDemo::ElectronNtuplerVIDDemo(const edm::ParameterSet& iConfig):
  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap")))
{

  //
  // Prepare tokens for all input collections and objects
  //

  pileupToken_ = consumes<edm::View<PileupSummaryInfo> >
    (iConfig.getParameter <edm::InputTag>
     ("pileup"));

  rhoToken_    = consumes<double> 
    (iConfig.getParameter <edm::InputTag>
     ("rho"));

  beamSpotToken_    = consumes<reco::BeamSpot> 
    (iConfig.getParameter <edm::InputTag>
     ("beamSpot"));

  // AOD tokens
  electronsToken_    = mayConsume<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
     ("electrons"));

  genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticles"));

  vtxToken_          = mayConsume<reco::VertexCollection>
    (iConfig.getParameter<edm::InputTag>
     ("vertices"));

  conversionsToken_ = mayConsume< reco::ConversionCollection >
    (iConfig.getParameter<edm::InputTag>
     ("conversions"));

  // MiniAOD tokens
  // For electrons, use the fact that pat::Electron can be cast into 
  // GsfElectron
  electronsMiniAODToken_    = mayConsume<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
     ("electronsMiniAOD"));

  genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticlesMiniAOD"));

  vtxMiniAODToken_          = mayConsume<reco::VertexCollection>
    (iConfig.getParameter<edm::InputTag>
     ("verticesMiniAOD"));

  conversionsMiniAODToken_ = mayConsume< reco::ConversionCollection >
    (iConfig.getParameter<edm::InputTag>
     ("conversionsMiniAOD"));

  edm::Service<TFileService> fs;
  electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");

  electronTree_->Branch("pvNTracks"    ,  &pvNTracks_);

  electronTree_->Branch("nPV"        ,  &nPV_);
  electronTree_->Branch("nPU"        ,  &nPU_);
  electronTree_->Branch("nPUTrue"    ,  &nPUTrue_);
  electronTree_->Branch("rho"        ,  &rho_);

  //electronTree_->Branch("nEle",  &nElectrons_ , "nEle/I");
  electronTree_->Branch("pt"  ,  &pt_    );
  electronTree_->Branch("eta" ,  &eta_ );
  electronTree_->Branch("phi" ,  &phi_ );

  electronTree_->Branch("dEtaIn",  &dEtaIn_);
  electronTree_->Branch("dPhiIn",  &dPhiIn_);
  electronTree_->Branch("hOverE",  &hOverE_);
  electronTree_->Branch("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta_);
  electronTree_->Branch("isoChargedHadrons"      , &isoChargedHadrons_);
  electronTree_->Branch("isoNeutralHadrons"      , &isoNeutralHadrons_);
  electronTree_->Branch("isoPhotons"             , &isoPhotons_);
  electronTree_->Branch("isoChargedFromPU"       , &isoChargedFromPU_);
  electronTree_->Branch("ooEmooP", &ooEmooP_);
  electronTree_->Branch("d0"     , &d0_);
  electronTree_->Branch("dz"     , &dz_);
  electronTree_->Branch("expectedMissingInnerHits", &expectedMissingInnerHits_);
  electronTree_->Branch("passConversionVeto", &passConversionVeto_);

  electronTree_->Branch("passVetoId"  ,  &passVetoId_ );
  electronTree_->Branch("passLooseId"  ,  &passLooseId_ );
  electronTree_->Branch("passMediumId" ,  &passMediumId_ );
  electronTree_->Branch("passTightId"  ,  &passTightId_ );

  electronTree_->Branch("isTrue"             , &isTrue_);

}


ElectronNtuplerVIDDemo::~ElectronNtuplerVIDDemo()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElectronNtuplerVIDDemo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;

  Int_t nPU=0, nPUTrue=0;
  Float_t rho=0;

  // Get Pileup info
  Handle<edm::View<PileupSummaryInfo> > pileupHandle;
  iEvent.getByToken(pileupToken_, pileupHandle);
  for( auto & puInfoElement : *pileupHandle){
    if( puInfoElement.getBunchCrossing() == 0 ){
      nPU     = puInfoElement.getPU_NumInteractions();
      nPUTrue = puInfoElement.getTrueNumInteractions();
    }
  }

  // Get rho value
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho  = *rhoH;

  // Get the beam spot
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_,theBeamSpot);  
  
  // Retrieve the collection of electrons from the event.
  // If we fail to retrieve the collection with the standard AOD
  // name, we next look for the one with the stndard miniAOD name.
  //   We use exactly the same handle for AOD and miniAOD formats
  // since pat::Electron objects can be recast as reco::GsfElectron objects.
  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  bool isAOD = true;
  iEvent.getByToken(electronsToken_, electrons);
  if( !electrons.isValid() ){
    isAOD = false;
    iEvent.getByToken(electronsMiniAODToken_,electrons);
  }
  
  // Get the MC collection
  Handle<edm::View<reco::GenParticle> > genParticles;
  if( isAOD )
    iEvent.getByToken(genParticlesToken_,genParticles);
  else
    iEvent.getByToken(genParticlesMiniAODToken_,genParticles);

  // Get PV
  edm::Handle<reco::VertexCollection> vertices;
  if( isAOD )
    iEvent.getByToken(vtxToken_, vertices);
  else
    iEvent.getByToken(vtxMiniAODToken_, vertices);
  
  if (vertices->empty()) return; // skip the event if no PV found
  
  // Find the first vertex in the collection that passes
  // good quality criteria
  VertexCollection::const_iterator firstGoodVertex = vertices->end();
  int firstGoodVertexIdx = 0;
  for (VertexCollection::const_iterator vtx = vertices->begin(); 
       vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    bool isFake = vtx->isFake();
    if( !isAOD )
      isFake = (vtx->chi2()==0 && vtx->ndof()==0);
    // Check the goodness
    if ( !isFake
	 &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
	 && fabs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
  }

  if ( firstGoodVertex==vertices->end() )
    return; // skip event if there are no good PVs

  // Get the conversions collection
  edm::Handle<reco::ConversionCollection> conversions;
  if(isAOD)
    iEvent.getByToken(conversionsToken_, conversions);
  else
    iEvent.getByToken(conversionsMiniAODToken_, conversions);

  // Get the electron ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  // If you need more info, check with the EGM group.
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  iEvent.getByToken(eleVetoIdMapToken_ ,veto_id_decisions);
  iEvent.getByToken(eleLooseIdMapToken_ ,loose_id_decisions);
  iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(eleTightIdMapToken_ ,tight_id_decisions);

  // Clear vectors
  //nElectrons_ = 0;
  /*  nPU_.clear();
  nPUTrue_.clear();
  rho_.clear();
  pt_.clear();
  eta_.clear();
  phi_.clear();
  dEtaIn_.clear();
  dPhiIn_.clear();
  hOverE_.clear();
  full5x5_sigmaIetaIeta_.clear();
  isoChargedHadrons_.clear();
  isoNeutralHadrons_.clear();
  isoPhotons_.clear();
  isoChargedFromPU_.clear();
  ooEmooP_.clear();
  d0_.clear();
  dz_.clear();
  expectedMissingInnerHits_.clear();
  passConversionVeto_.clear();     
  //
  passVetoId_ .clear();
  passLooseId_ .clear();
  passMediumId_.clear();
  passTightId_ .clear();
  //
  isTrue_.clear();*/

  // Loop over electrons
  for (size_t i = 0; i < electrons->size(); ++i){
    const auto el = electrons->ptrAt(i);

    // Kinematics
    if( el->pt() < 10 ) // keep only electrons above 10 GeV
      continue;
    
    //nElectrons_++;

    //
    // Save electron kinematics
    //
    pt_  = el->pt() ;
    eta_ = el->superCluster()->eta() ;
    phi_ = el->superCluster()->phi() ;

    nPU_=nPU;
    nPUTrue_=nPUTrue;
    rho_=rho;

    // ID and matching
    dEtaIn_= el->deltaEtaSuperClusterTrackAtVtx() ;
    dPhiIn_= el->deltaPhiSuperClusterTrackAtVtx() ;
    hOverE_= el->hcalOverEcal() ;
    full5x5_sigmaIetaIeta_= el->full5x5_sigmaIetaIeta() ;
    // |1/E-1/p| = |1/E - EoverPinner/E| is computed below
    // The if protects against ecalEnergy == inf or zero
    // (always the case for miniAOD for electrons <5 GeV)
    if( el->ecalEnergy() == 0 ){
      printf("Electron energy is zero!\n");
      ooEmooP_= 1e30 ;
    }else if( !std::isfinite(el->ecalEnergy())){
      printf("Electron energy is not finite!\n");
      ooEmooP_= 1e30 ;
    }else{
      ooEmooP_= fabs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy() ) ;
    }
    
    // Isolation
    GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();
    // Compute isolation with delta beta correction for PU
    isoChargedHadrons_= pfIso.sumChargedHadronPt ;
    isoNeutralHadrons_= pfIso.sumNeutralHadronEt ;
    isoPhotons_= pfIso.sumPhotonEt ;
    isoChargedFromPU_= pfIso.sumPUPt ;

    // Impact parameter
    reco::GsfTrackRef theTrack = el->gsfTrack();
    d0_= (-1) * theTrack->dxy(firstGoodVertex->position() ) ;
    dz_= theTrack->dz( firstGoodVertex->position() ) ;

    // Conversion rejection
    expectedMissingInnerHits_=el->gsfTrack()->hitPattern()
					.numberOfHits(reco::HitPattern::MISSING_INNER_HITS) ;

    bool passConvVeto = !ConversionTools::hasMatchedConversion(*el, 
							       conversions,
							       theBeamSpot->position());
    passConversionVeto_= (int) passConvVeto ;

    //
    // Look up and save the ID decisions
    // 
    bool isPassVeto  = (*veto_id_decisions)[el];
    bool isPassLoose  = (*loose_id_decisions)[el];
    bool isPassMedium = (*medium_id_decisions)[el];
    bool isPassTight  = (*tight_id_decisions)[el];
    passVetoId_= (int)isPassVeto  ;
    passLooseId_ = (int)isPassLoose ;
    passMediumId_ = (int)isPassMedium;
    passTightId_ = (int)isPassTight ;

    // Save MC truth match
    isTrue_= matchToTruth( el, genParticles) ;

    // Save the info
    electronTree_->Fill();


   }
   

}


// ------------ method called once each job just before starting event loop  ------------
void 
ElectronNtuplerVIDDemo::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronNtuplerVIDDemo::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
ElectronNtuplerVIDDemo::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
ElectronNtuplerVIDDemo::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
ElectronNtuplerVIDDemo::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
ElectronNtuplerVIDDemo::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronNtuplerVIDDemo::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int ElectronNtuplerVIDDemo::matchToTruth(const edm::Ptr<reco::GsfElectron> el, 
				  const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

  // 
  // Explicit loop and geometric matching method (advised by Josh Bendavid)
  //

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  const reco::Candidate *closestElectron = 0;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestElectron = particle;
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if( !(closestElectron != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }

  // 
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);

  if( ancestorPID == -999 && ancestorStatus == -999 ){
    // No non-electron parent??? This should never happen.
    // Complain.
    printf("ElectronNtupler: ERROR! Electron does not apper to have a non-electron parent\n");
    return UNMATCHED;
  }
  
  if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
    return TRUE_NON_PROMPT_ELECTRON;

  if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
    return TRUE_ELECTRON_FROM_TAU;

  // What remains is true prompt electrons
  return TRUE_PROMPT_ELECTRON;
}

void ElectronNtuplerVIDDemo::findFirstNonElectronMother(const reco::Candidate *particle,
						 int &ancestorPID, int &ancestorStatus){

  if( particle == 0 ){
    printf("ElectronNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-electron parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 11 ){
    findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }

  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronNtuplerVIDDemo);
