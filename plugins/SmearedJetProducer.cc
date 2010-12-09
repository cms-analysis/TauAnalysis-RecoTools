#include "TauAnalysis/RecoTools/plugins/SmearedJetProducer.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Candidate/interface/Particle.h"

SmearedJetProducer::SmearedJetProducer(const edm::ParameterSet& cfg)
  : src_(cfg.getParameter<edm::InputTag>("src")),
    jecUncertainty_(0),
    isJECuncertaintyFromFile_(false),
    shiftByJECuncertainty_(0.),
    smearingModule_(0)
{
  if ( cfg.exists("shiftByJECuncertainty") ) {
    if ( cfg.exists("jecUncertaintyInputFileName") ) {
      edm::FileInPath jecUncertaintyInputFileName = cfg.getParameter<edm::FileInPath>("jecUncertaintyInputFileName");
      if ( jecUncertaintyInputFileName.isLocal() ) {  
        jecUncertainty_ = new JetCorrectionUncertainty(jecUncertaintyInputFileName.fullPath().data());
        isJECuncertaintyFromFile_ = true;
      } else {
        throw cms::Exception("SmearedJetProducer") 
          << " Failed to find file = " << jecUncertaintyInputFileName << " !!";
      }  
    }

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
  delete jecUncertainty_;
  delete smearingModule_;
}

void SmearedJetProducer::produce(edm::Event& evt, const edm::EventSetup& es) 
{
  std::auto_ptr<pat::JetCollection> smearedJets(new pat::JetCollection);

  edm::Handle<pat::JetCollection> patJets;
  evt.getByLabel(src_, patJets);

  if ( shiftByJECuncertainty_ != 0. && !isJECuncertaintyFromFile_ ) {
    edm::ESHandle<JetCorrectorParametersCollection> jetCorrParameterSet;
    es.get<JetCorrectionsRecord>().get(jetCorrPayloadName_, jetCorrParameterSet); 
    const JetCorrectorParameters& jetCorrParameters = (*jetCorrParameterSet)[jetCorrUncertaintyTag_];
    delete jecUncertainty_;
    jecUncertainty_ = new JetCorrectionUncertainty(jetCorrParameters);
  }

  for ( pat::JetCollection::const_iterator patJet = patJets->begin();
	patJet != patJets->end(); ++patJet ) {
    pat::Jet smearedJet = (*patJet);
    if ( shiftByJECuncertainty_ != 0. ) {
      jecUncertainty_->setJetEta(patJet->eta());
      jecUncertainty_->setJetPt(patJet->pt());

      double shift = shiftByJECuncertainty_*jecUncertainty_->getUncertainty(true);

      reco::Particle::LorentzVector patJetP4 = patJet->p4();
      smearedJet.setP4((1. + shift)*patJetP4);

      //std::cout << "patJet: Pt = " << patJet->pt() << "," 
      //	  << " eta = " << patJet->eta() << ", phi = " << patJet->phi() << std::endl;
      //std::cout << "smearedJet: Pt = " << smearedJet.pt() << "," 
      //	  << " eta = " << smearedJet.eta() << ", phi = " << smearedJet.phi() << std::endl;
    }

    if ( smearingModule_ ) smearingModule_->smear(smearedJet);
    
    smearedJets->push_back(smearedJet);
  }
  
  evt.put(smearedJets);
} 

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(SmearedJetProducer);
