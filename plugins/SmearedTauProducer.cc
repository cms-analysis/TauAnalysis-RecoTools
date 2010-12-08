#include "TauAnalysis/RecoTools/plugins/SmearedTauProducer.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Candidate/interface/Particle.h"

#include <TMath.h>

SmearedTauProducer::SmearedTauProducer(const edm::ParameterSet& cfg)
  : src_(cfg.getParameter<edm::InputTag>("src")),
    shiftByJECuncertainty_(0.),
    smearingModule_(0) 
{
  if ( cfg.exists("shiftByJECuncertainty") ) {
    shiftByJECuncertainty_ = cfg.getParameter<double>("shiftByJECuncertainty");

    jetCorrPayloadName_ = cfg.getParameter<std::string>("jetCorrPayloadName");
    jetCorrUncertaintyTag_ = cfg.getParameter<std::string>("jetCorrUncertaintyTag");

    jecFlavorUncertainty_ = cfg.getParameter<double>("jecFlavorUncertainty");
  }

  if ( cfg.exists("fileName") ) {
    smearingModule_ = new SmearedParticleMaker<pat::Tau, GenJetRetriever<pat::Tau> >(cfg);
  }

  produces<pat::TauCollection>();
}

SmearedTauProducer::~SmearedTauProducer()
{
  delete smearingModule_;
}

void SmearedTauProducer::produce(edm::Event& evt, const edm::EventSetup& es) 
{
  std::auto_ptr<pat::TauCollection> smearedTaus(new pat::TauCollection);
  
  edm::Handle<pat::TauCollection> patTaus;
  evt.getByLabel(src_, patTaus);
  
  JetCorrectionUncertainty* jecUncertainty = 0;
  if ( shiftByJECuncertainty_ != 0. ) {
    edm::ESHandle<JetCorrectorParametersCollection> jetCorrParameterSet;
    es.get<JetCorrectionsRecord>().get(jetCorrPayloadName_, jetCorrParameterSet); 
    const JetCorrectorParameters& jetCorrParameters = (*jetCorrParameterSet)[jetCorrUncertaintyTag_];
    jecUncertainty = new JetCorrectionUncertainty(jetCorrParameters);
  }

  for ( pat::TauCollection::const_iterator patTau = patTaus->begin();
	patTau != patTaus->end(); ++patTau ) {
    pat::Tau smearedTau = (*patTau);
    
    if ( shiftByJECuncertainty_ != 0. ) {
      jecUncertainty->setJetEta(patTau->eta());
      jecUncertainty->setJetPt(patTau->pt());
      double shift = shiftByJECuncertainty_*jecUncertainty->getUncertainty(true);

//--- add flaor uncertainty
//   (to account for difference in uncertainty on jet response
//    between tau-jets and quark/gluon jets)
      shift = TMath::Sqrt(shift*shift + jecFlavorUncertainty_*jecFlavorUncertainty_);
      
      reco::Particle::LorentzVector patTauP4 = patTau->p4();
      smearedTau.setP4((1. + shift)*patTauP4);

      std::cout << "patTau: Pt = " << patTau->pt() << "," 
		<< " eta = " << patTau->eta() << ", phi = " << patTau->phi() << std::endl;
      std::cout << "smearedTau: Pt = " << smearedTau.pt() << "," 
		<< " eta = " << smearedTau.eta() << ", phi = " << smearedTau.phi() << std::endl;
    }

    smearingModule_->smear(smearedTau);
    
    smearedTaus->push_back(smearedTau);
  }
  
  evt.put(smearedTaus);
  
  delete jecUncertainty;
} 

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(SmearedTauProducer);
