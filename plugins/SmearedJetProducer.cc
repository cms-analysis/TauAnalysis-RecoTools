#include "TauAnalysis/RecoTools/plugins/SmearedJetProducer.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Candidate/interface/Particle.h"

SmearedJetProducer::SmearedJetProducer(const edm::ParameterSet& cfg)
  : src_(cfg.getParameter<edm::InputTag>("src")),
    shiftByJECuncertainty_(0.),
    smearingModule_(0)
{
  if ( cfg.exists("shiftByJECuncertainty") ) {
    shiftByJECuncertainty_ = cfg.getParameter<double>("shiftByJECuncertainty");

    jetCorrPayloadName_ = cfg.getParameter<std::string>("jetCorrPayloadName");
    jetCorrUncertaintyTag_ = cfg.getParameter<std::string>("jetCorrUncertaintyTag");
  }

  if ( cfg.exists("fileName") ) {
    smearingModule_ = new SmearedParticleMaker<pat::Jet, GenJetRetriever<pat::Jet> >(cfg);
  }

  produces<pat::JetCollection>();
}

SmearedJetProducer::~SmearedJetProducer()
{
  delete smearingModule_;
}

void SmearedJetProducer::produce(edm::Event& evt, const edm::EventSetup& es) 
{
  std::auto_ptr<pat::JetCollection> smearedJets(new pat::JetCollection);
  
  edm::Handle<pat::JetCollection> patJets;
  evt.getByLabel(src_, patJets);
  
  JetCorrectionUncertainty* jecUncertainty = 0;
  if ( shiftByJECuncertainty_ != 0. ) {
    edm::ESHandle<JetCorrectorParametersCollection> jetCorrParameterSet;
    es.get<JetCorrectionsRecord>().get(jetCorrPayloadName_, jetCorrParameterSet); 
    const JetCorrectorParameters& jetCorrParameters = (*jetCorrParameterSet)[jetCorrUncertaintyTag_];
    jecUncertainty = new JetCorrectionUncertainty(jetCorrParameters);
  }

  for ( pat::JetCollection::const_iterator patJet = patJets->begin();
	patJet != patJets->end(); ++patJet ) {
    pat::Jet smearedJet = (*patJet);
    
    if ( shiftByJECuncertainty_ != 0. ) {
      jecUncertainty->setJetEta(patJet->eta());
      jecUncertainty->setJetPt(patJet->pt());
      double shift = shiftByJECuncertainty_*jecUncertainty->getUncertainty(true);
      
      reco::Particle::LorentzVector patJetP4 = patJet->p4();
      smearedJet.setP4((1. + shift)*patJetP4);

      std::cout << "patJet: Pt = " << patJet->pt() << "," 
		<< " eta = " << patJet->eta() << ", phi = " << patJet->phi() << std::endl;
      std::cout << "smearedJet: Pt = " << smearedJet.pt() << "," 
		<< " eta = " << smearedJet.eta() << ", phi = " << smearedJet.phi() << std::endl;
    }

    smearingModule_->smear(smearedJet);
    
    smearedJets->push_back(smearedJet);
  }
  
  evt.put(smearedJets);
  
  delete jecUncertainty;
} 

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(SmearedJetProducer);
