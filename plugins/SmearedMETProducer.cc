#include "TauAnalysis/RecoTools/plugins/SmearedMETProducer.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

typedef edm::View<reco::Candidate> recoCandidateView;

SmearedMETProducer::SmearedMETProducer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
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
    //std::cout << " ptParticle = " << particle->pt() << std::endl;
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
  //std::cout << "<SmearedMETProducer::produce>:" << std::endl;
  //std::cout << " moduleLabel = " << moduleLabel_ << std::endl;

  reco::Particle::LorentzVector metCorrection;

  for ( std::vector<smearedParticleType>::const_iterator smearedParticleCollection = smearedParticleCollections_.begin();
	smearedParticleCollection != smearedParticleCollections_.end(); ++smearedParticleCollection ) {
    edm::Handle<recoCandidateView> originalParticles;
    //std::cout << " srcOriginal = " << smearedParticleCollection->srcOriginal_.label() << std::endl;
    evt.getByLabel(smearedParticleCollection->srcOriginal_, originalParticles);
    reco::Particle::LorentzVector p4original = compSumP4(*originalParticles);
    //std::cout << "--> ptOriginal = " << p4original.pt() << std::endl;
    metCorrection += p4original;

    edm::Handle<recoCandidateView> smearedParticles;
    //std::cout << " srcSmeared = " << smearedParticleCollection->srcSmeared_.label() << std::endl;
    evt.getByLabel(smearedParticleCollection->srcSmeared_, smearedParticles);
    reco::Particle::LorentzVector p4smeared = compSumP4(*smearedParticles);
    //std::cout << "--> ptSmeared = " << p4smeared.pt() << std::endl;
    metCorrection -= p4smeared;

    if ( smearedParticleCollection->smearByResolutionUncertainty_ ) {
      double sumPt = 0.;
      for ( recoCandidateView::const_iterator particle = originalParticles->begin();
	    particle != originalParticles->end(); ++particle ) {
	//std::cout << " ptParticle = " << particle->pt() << std::endl;
	sumPt += particle->pt();
      }
      
      //std::cout << "--> sumPt = " << sumPt << std::endl;

      if ( sumPt > 0. && p4original.pt() > 0. ) {
        double sigma = smearedParticleCollection->smearByResolutionUncertainty_*TMath::Sqrt(sumPt);

        double ptSmearResUncertainty = rnd_.Gaus(0., sigma);
        //std::cout << " ptSmearResUncertainty = " << ptSmearResUncertainty << std::endl;

        double pxSmearResUncertainty = p4original.px()*(ptSmearResUncertainty/p4original.pt());
        double pySmearResUncertainty = p4original.py()*(ptSmearResUncertainty/p4original.pt());
 
        metCorrection += math::XYZTLorentzVector(pxSmearResUncertainty, pySmearResUncertainty, 0., 0.);
      }
    }
  }

  std::auto_ptr<pat::METCollection> smearedMETs(new pat::METCollection);

  edm::Handle<pat::METCollection> originalMETs;
  evt.getByLabel(src_, originalMETs);

  for ( pat::METCollection::const_iterator originalMET = originalMETs->begin();
	originalMET != originalMETs->end(); ++originalMET ) {
    //std::cout << "originalMET: px = " << originalMET->px() << ", py = " << originalMET->py() << std::endl;
    
    pat::MET smearedMET(*originalMET);
    
    double smearedMETpx = originalMET->px() + metCorrection.px();
    double smearedMETpy = originalMET->py() + metCorrection.py();
    double smearedMETpt = TMath::Sqrt(smearedMETpx*smearedMETpx + smearedMETpy*smearedMETpy);

    smearedMET.setP4(math::XYZTLorentzVector(smearedMETpx, smearedMETpy, 0., smearedMETpt));
    //std::cout << "smearedMET: px = " << smearedMET.px() << ", py = " << smearedMET.py() << std::endl;
    
    smearedMETs->push_back(smearedMET);
  }
  
  evt.put(smearedMETs);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(SmearedMETProducer);
