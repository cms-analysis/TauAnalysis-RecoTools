#include "TauAnalysis/RecoTools/plugins/SmearedTauProducer.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Candidate/interface/Particle.h"

#include <TMath.h>

SmearedTauProducer::SmearedTauProducer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label")),
    src_(cfg.getParameter<edm::InputTag>("src")),
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
    varyByNsigmas_ = ( cfg.exists("varyByNsigmas") ) ?
      cfg.getParameter<double>("varyByNsigmas") : 1.0;

    jetCorrPayloadName_ = cfg.getParameter<std::string>("jetCorrPayloadName");
    jetCorrUncertaintyTag_ = cfg.getParameter<std::string>("jetCorrUncertaintyTag");

    jecFlavorUncertainty_ = cfg.getParameter<double>("jecFlavorUncertainty");
  }

  produces<pat::TauCollection>();
}

SmearedTauProducer::~SmearedTauProducer()
{
  delete jecUncertainty_;
  delete smearingModule_;
}

void SmearedTauProducer::produce(edm::Event& evt, const edm::EventSetup& es) 
{
  //std::cout << "<SmearedTauProducer::produce>:" << std::endl;
  //std::cout << " moduleLabel = " << moduleLabel_ << std::endl;
  //std::cout << " shiftByJECuncertainty = " << shiftByJECuncertainty_ << std::endl;

  std::auto_ptr<pat::TauCollection> smearedTaus(new pat::TauCollection);
  
  edm::Handle<pat::TauCollection> patTaus;
  evt.getByLabel(src_, patTaus);
  
  if ( shiftByJECuncertainty_ != 0. && !isJECuncertaintyFromFile_ ) {
    edm::ESHandle<JetCorrectorParametersCollection> jetCorrParameterSet;
    es.get<JetCorrectionsRecord>().get("AK5PF", jetCorrParameterSet); 
    const JetCorrectorParameters& jetCorrParameters = (*jetCorrParameterSet)[jetCorrUncertaintyTag_];
    delete jecUncertainty_;
    jecUncertainty_ = new JetCorrectionUncertainty(jetCorrParameters);
  }

  for ( pat::TauCollection::const_iterator patTau = patTaus->begin();
	patTau != patTaus->end(); ++patTau ) {
    pat::Tau smearedTau = (*patTau);
    
    if ( shiftByJECuncertainty_ != 0. ) {
      jecUncertainty_->setJetEta(patTau->eta());
      jecUncertainty_->setJetPt(patTau->pt());

      double jecUncertainty = jecUncertainty_->getUncertainty(true);
      //std::cout << " jecUncertainty = " << jecUncertainty << std::endl;

//--- add flavor uncertainty
//   (to account for difference in uncertainty on jet response
//    between tau-jets and quark/gluon jets)
      double shift = shiftByJECuncertainty_*TMath::Sqrt(jecUncertainty*jecUncertainty + jecFlavorUncertainty_*jecFlavorUncertainty_);
      //std::cout << " shift = " << shift << std::endl;

      shift *= varyByNsigmas_;

      reco::Particle::LorentzVector patTauP4 = patTau->p4();
      smearedTau.setP4((1. + shift)*patTauP4);

      smearedTau.addUserFloat("shiftByJECuncertainty", shift);

      //std::cout << "patTau: Pt = " << patTau->pt() << "," 
      //          << " eta = " << patTau->eta() << ", phi = " << patTau->phi() << std::endl;
      //std::cout << "smearedTau: Pt = " << smearedTau.pt() << "," 
      //	  << " eta = " << smearedTau.eta() << ", phi = " << smearedTau.phi() << std::endl;
    }

    if ( smearingModule_ ) smearingModule_->smear(smearedTau);
    
    smearedTaus->push_back(smearedTau);
  }
  
  evt.put(smearedTaus);
} 

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(SmearedTauProducer);
