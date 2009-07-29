#include "TauAnalysis/RecoTools/interface/PATElecTauPairZeeHypothesisAlgorithm.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "PhysicsTools/Utilities/interface/deltaR.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"

#include <TMath.h>

typedef edm::Ptr<PATElecTauPair> PATElecTauPairPtr;
typedef edm::Ptr<reco::GenParticle> GenParticlePtr;
typedef edm::Ptr<reco::CaloJet> CaloJetPtr;
typedef edm::Ptr<reco::PFJet> PFJetPtr;
typedef edm::Ptr<reco::Track> TrackPtr;
typedef edm::Ptr<reco::GsfElectron> GsfElectronPtr;
typedef edm::Ptr<reco::GsfTrack> GsfTrackPtr;

PATElecTauPairZeeHypothesisAlgorithm::PATElecTauPairZeeHypothesisAlgorithm(const edm::ParameterSet& cfg)
{
  srcGenElectrons_ = cfg.getParameter<edm::InputTag>("genElectronsFromZsSource");
  srcCaloJets_ = cfg.getParameter<edm::InputTag>("caloJetSource");
  srcPFJets_ = cfg.getParameter<edm::InputTag>("pfJetSource");
  srcTracks_ = cfg.getParameter<edm::InputTag>("trackSource");
  srcGsfElectrons_ = cfg.getParameter<edm::InputTag>("gsfElectronSource");
  srcGsfTracks_ = cfg.getParameter<edm::InputTag>("gsfTrackSource");

  tkminPixelHits_ = cfg.getParameter<int>("tkminPixelHits");
  tkminTrackerHits_ = cfg.getParameter<int>("tkminTrackerHits");
  tkmaxChi2_ = cfg.getParameter<double>("tkmaxChi2");

  dRmatch_ = cfg.getParameter<double>("dRmatch");

  verbosity_ = cfg.getUntrackedParameter<int>("verbosity", 0);
}

PATElecTauPairZeeHypothesisAlgorithm::~PATElecTauPairZeeHypothesisAlgorithm() 
{
//--- nothing to be done yet...
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

template <typename T>
edm::Ptr<T> getBestMatchParticle(const edm::Event& evt, const edm::InputTag& src,
				 const reco::Particle::LorentzVector& direction, double dRmatch)
{
  edm::Handle<edm::View<T> > particles;
  evt.getByLabel(src, particles);

  edm::Ptr<T> particle_bestMatch;

  double dRmin = dRmatch;
  for ( unsigned idx = 0, numParticles = particles->size(); 
	idx < numParticles; ++idx ) {
    edm::Ptr<T> particle = particles->ptrAt(idx);
    if ( reco::deltaR(particle->p4(), direction) < dRmin ) {
      particle_bestMatch = particle;
    }
  }
  
  return particle_bestMatch;
}

template <typename T>
edm::Ptr<T> getBestMatchTrack(const edm::Event& evt, const edm::InputTag& src,
			      const reco::Particle::LorentzVector& direction, double dRmatch,
			      int tkminPixelHits, int tkminTrackerHits, double tkmaxChi2)
{
  edm::Handle<edm::View<T> > tracks;
  evt.getByLabel(src, tracks);

  edm::Ptr<T> track_bestMatch;

  double ptMax = -1.;
  for ( unsigned idx = 0, numTracks = tracks->size(); 
	idx < numTracks; ++idx ) {
    edm::Ptr<T> track = tracks->ptrAt(idx);

    if ( track->hitPattern().numberOfValidPixelHits() >= tkminPixelHits &&
	 track->numberOfValidHits() >= tkminTrackerHits &&
	 track->normalizedChi2() < tkmaxChi2 ) {
      double pt = track->pt();
      if ( reco::deltaR(track->eta(), track->phi(), direction.eta(), direction.phi()) < dRmatch && pt > ptMax ) {
	track_bestMatch = track;
	ptMax = pt;
      }
    }
  }
  
  return track_bestMatch;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

reco::Particle::LorentzVector getP4(const CaloJetPtr& caloJet, const PFJetPtr& pfJet, const TrackPtr& track,
				    const GsfElectronPtr& gsfElectron, const GsfTrackPtr& gsfTrack,
				    int type, bool& isValid)
{
  assert(type >= PATElecTauPairZeeHypothesis::kCaloJet && type <= PATElecTauPairZeeHypothesis::kGsfTrack);

  isValid = true;
  
  if ( type == PATElecTauPairZeeHypothesis::kCaloJet && caloJet.isAvailable() ) {
    return caloJet->p4();
  } else if ( type == PATElecTauPairZeeHypothesis::kPFJet && pfJet.isAvailable() ) {
    return pfJet->p4();
  } else if ( type == PATElecTauPairZeeHypothesis::kRecoTrack && track.isAvailable() ) {
    return reco::Particle::LorentzVector(track->px(), track->py(), track->pz(), track->p());
  } else if ( type == PATElecTauPairZeeHypothesis::kGsfElecron && gsfElectron.isAvailable() ) {
    return gsfElectron->p4();
  } else if ( type == PATElecTauPairZeeHypothesis::kGsfTrack && gsfTrack.isAvailable() ) {
    return reco::Particle::LorentzVector(gsfTrack->px(), gsfTrack->py(), gsfTrack->pz(), gsfTrack->p());
  }

  isValid = false;

  return reco::Particle::LorentzVector(0,0,0,0);
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

PATElecTauPairZeeHypothesis PATElecTauPairZeeHypothesisAlgorithm::buildZeeHypothesis(edm::Ptr<PATElecTauPair> elecTauPair,
										     const edm::Event& evt, const edm::EventSetup& es)
{
  PATElecTauPairZeeHypothesis ZeeHypothesis(elecTauPair);

  const reco::Particle::LorentzVector& direction1 = elecTauPair->leg1()->p4();
  const reco::Particle::LorentzVector& direction2 = elecTauPair->leg2()->p4();

  ZeeHypothesis.genElec1_ = getBestMatchParticle<reco::GenParticle>(evt, srcGenElectrons_, direction1, dRmatch_);
  ZeeHypothesis.genElec2_ = getBestMatchParticle<reco::GenParticle>(evt, srcGenElectrons_, direction2, dRmatch_);

  ZeeHypothesis.elec1matchedCaloJet_ = getBestMatchParticle<reco::CaloJet>(evt, srcCaloJets_, direction1, dRmatch_);
  ZeeHypothesis.elec1matchedPFJet_ = getBestMatchParticle<reco::PFJet>(evt, srcPFJets_, direction1, dRmatch_);
  ZeeHypothesis.elec1matchedTrack_ = 
    getBestMatchTrack<reco::Track>(evt, srcTracks_, direction1, dRmatch_, tkminPixelHits_, tkminTrackerHits_, tkmaxChi2_);
  ZeeHypothesis.elec1matchedGsfElectron_ = getBestMatchParticle<reco::GsfElectron>(evt, srcGsfElectrons_, direction1, dRmatch_);
  ZeeHypothesis.elec1matchedGsfTrack_ = 
    getBestMatchTrack<reco::GsfTrack>(evt, srcGsfTracks_, direction1, dRmatch_, tkminPixelHits_, tkminTrackerHits_, tkmaxChi2_);

  ZeeHypothesis.elec2matchedCaloJet_ = getBestMatchParticle<reco::CaloJet>(evt, srcCaloJets_, direction2, dRmatch_);
  ZeeHypothesis.elec2matchedPFJet_ = getBestMatchParticle<reco::PFJet>(evt, srcPFJets_, direction2, dRmatch_);
  ZeeHypothesis.elec2matchedTrack_ = 
    getBestMatchTrack<reco::Track>(evt, srcTracks_, direction2, dRmatch_, tkminPixelHits_, tkminTrackerHits_, tkmaxChi2_);
  ZeeHypothesis.elec2matchedGsfElectron_ = getBestMatchParticle<reco::GsfElectron>(evt, srcGsfElectrons_, direction2, dRmatch_);
  ZeeHypothesis.elec2matchedGsfTrack_ = 
    getBestMatchTrack<reco::GsfTrack>(evt, srcGsfTracks_, direction2, dRmatch_, tkminPixelHits_, tkminTrackerHits_, tkmaxChi2_);

  //edm::ESHandle<HepPDT::ParticleDataTable> pdt;
  //es.getData(pdt);
  //double nominalZ0mass = pdt->particle(23)->mass();
  double nominalZ0mass = 91.188;

  double visMassBestMatch = 1.e+3;
  reco::Particle::LorentzVector p4Elec1bestMatch;
  int elec1bestMatchType = PATElecTauPairZeeHypothesis::kUndefined;
  reco::Particle::LorentzVector p4Elec2bestMatch;
  int elec2bestMatchType = PATElecTauPairZeeHypothesis::kUndefined;
  for ( int elec1Type = PATElecTauPairZeeHypothesis::kCaloJet;
	elec1Type <= PATElecTauPairZeeHypothesis::kGsfTrack; ++elec1Type ) {
    
    bool isValid1;
    reco::Particle::LorentzVector elec1Momentum 
      = getP4(ZeeHypothesis.elec1matchedCaloJet_, ZeeHypothesis.elec1matchedPFJet_, ZeeHypothesis.elec1matchedTrack_,
	      ZeeHypothesis.elec1matchedGsfElectron_, ZeeHypothesis.elec1matchedGsfTrack_, elec1Type, isValid1);
    if ( !isValid1 ) continue;
    
    for ( int elec2Type = PATElecTauPairZeeHypothesis::kCaloJet;
	elec2Type <= PATElecTauPairZeeHypothesis::kGsfTrack; ++elec2Type ) {

      bool isValid2;
      reco::Particle::LorentzVector elec2Momentum = 
	getP4(ZeeHypothesis.elec2matchedCaloJet_, ZeeHypothesis.elec2matchedPFJet_, ZeeHypothesis.elec2matchedTrack_,
	      ZeeHypothesis.elec2matchedGsfElectron_, ZeeHypothesis.elec2matchedGsfTrack_, elec2Type, isValid2);
      if ( !isValid2 ) continue;

      double visMass = (elec1Momentum + elec2Momentum).mass();
      if ( TMath::Abs(visMass - nominalZ0mass) < TMath::Abs(visMassBestMatch - nominalZ0mass) ) {
        p4Elec1bestMatch = elec1Momentum;
	elec1bestMatchType = elec1Type;
        p4Elec2bestMatch = elec2Momentum;
	elec2bestMatchType = elec2Type;
	visMassBestMatch = visMass;
      }
    }
  }

  ZeeHypothesis.p4Elec1bestMatch_ = p4Elec1bestMatch;
  ZeeHypothesis.typeElec1bestMatch_ = elec1bestMatchType;
  ZeeHypothesis.p4Elec2bestMatch_ = p4Elec2bestMatch;
  ZeeHypothesis.typeElec2bestMatch_ = elec2bestMatchType;

  if ( verbosity_ ) {
    if ( ZeeHypothesis.genElec1_.isAvailable() &&
	 ZeeHypothesis.genElec2_.isAvailable() ) {
      std::cout << "<PATElecTauPairZeeHypothesisAlgorithm::buildZeeHypothesis>:" << std::endl;
      std::cout << " PATElecTauPair matches generated Z --> e+ e- decay," 
		<< " visMassBestMatch = " << visMassBestMatch << std::endl;
    }
  }

  return ZeeHypothesis;
}
