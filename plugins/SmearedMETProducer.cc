#include "TauAnalysis/RecoTools/plugins/SmearedMETProducer.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

typedef edm::View<reco::Candidate> recoCandidateView;

SmearedMETProducer::SmearedMETProducer(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("src");

  edm::ParameterSet cfgSmearedParticleCollections = cfg.getParameter<edm::ParameterSet>("smearedParticles");
  typedef std::vector<std::string> vstring;
  vstring smearedParticleCollectionNames = cfgSmearedParticleCollections.getParameterNamesForType<edm::ParameterSet>(); 
  for ( vstring::const_iterator smearedParticleCollectionName = smearedParticleCollectionNames.begin();
	smearedParticleCollectionName != smearedParticleCollectionNames.end(); ++smearedParticleCollectionName ) {
    edm::ParameterSet cfgSmearedParticleCollection = 
      cfgSmearedParticleCollections.getParameter<edm::ParameterSet>(*smearedParticleCollectionName);
    smearedParticleCollections_.push_back(smearedParticleType(cfgSmearedParticleCollection));
  }

  produces<pat::METCollection>();
}

SmearedMETProducer::~SmearedMETProducer()
{
// nothing to be done yet...
}
    
reco::Particle::LorentzVector compSumP4(const recoCandidateView& particles)
{
  reco::Particle::LorentzVector sumP4;
  
  for ( recoCandidateView::const_iterator particle = particles.begin();
	particle != particles.end(); ++particle ) {
    sumP4 += particle->p4();
  }

  return sumP4;
}

void updateMETcorrection(reco::Particle::LorentzVector& metCorr,
			 const reco::Particle::LorentzVector& projDirection_parallel, double corr_parallel, double corr_perpendicular)
{
  double projParallelCosPhi = TMath::Cos(projDirection_parallel.phi());
  double projParallelSinPhi = TMath::Sin(projDirection_parallel.phi());

  double projPerpendicularCosPhi = TMath::Cos(projDirection_parallel.phi() - 0.5*TMath::Pi());
  double projPerpendicularSinPhi = TMath::Sin(projDirection_parallel.phi() - 0.5*TMath::Pi());

  double metCorrDeltaPx = corr_parallel*projParallelCosPhi + corr_perpendicular*projPerpendicularCosPhi;
  double metCorrDeltaPy = corr_parallel*projParallelSinPhi + corr_perpendicular*projPerpendicularSinPhi;

  double metCorrPx = metCorr.px() + metCorrDeltaPx;
  double metCorrPy = metCorr.py() + metCorrDeltaPy;

  metCorr = math::XYZTLorentzVector(metCorrPx, metCorrPy, 0., 0.);
}

void SmearedMETProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  reco::Particle::LorentzVector metCorrection;

  for ( std::vector<smearedParticleType>::const_iterator smearedParticleCollection = smearedParticleCollections_.begin();
	smearedParticleCollection != smearedParticleCollections_.end(); ++smearedParticleCollection ) {
    edm::Handle<recoCandidateView> originalParticles;
    evt.getByLabel(smearedParticleCollection->srcOriginal_, originalParticles);
    metCorrection += compSumP4(*originalParticles);
    edm::Handle<recoCandidateView> smearedParticles;
    evt.getByLabel(smearedParticleCollection->srcOriginal_, smearedParticles);
    metCorrection -= compSumP4(*smearedParticles);
  }

  std::auto_ptr<pat::METCollection> smearedMETs(new pat::METCollection);

  edm::Handle<pat::METCollection> originalMETs;
  evt.getByLabel(src_, originalMETs);

  for ( pat::METCollection::const_iterator originalMET = originalMETs->begin();
	originalMET != originalMETs->end(); ++originalMET ) {
    pat::MET smearedMET(*originalMET);

    double smearedMETpx = originalMET->px() + metCorrection.px();
    double smearedMETpy = originalMET->py() + metCorrection.py();
    double smearedMETpt = TMath::Sqrt(smearedMETpx*smearedMETpx + smearedMETpy*smearedMETpy);

    smearedMET.setP4(math::XYZTLorentzVector(smearedMETpx, smearedMETpy, 0., smearedMETpt));
    
    smearedMETs->push_back(smearedMET);
  }

  evt.put(smearedMETs);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(SmearedMETProducer);
