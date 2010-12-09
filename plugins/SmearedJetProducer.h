/**
   \class   SmearedJetProducer

   \brief   Class for estimating systematic uncertainties on jets

   Author: Mike Bachtis, Wisconsin
	   Christian Veelken, UC Davis
*/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "TauAnalysis/RecoTools/interface/SmearedParticleMaker.h"

#include <string>

class SmearedJetProducer : public edm::EDProducer  
{
 public:
  explicit SmearedJetProducer(const edm::ParameterSet&);
  ~SmearedJetProducer();
    
 private:
  void produce(edm::Event&, const edm::EventSetup&);

// ----------member data ---------------------------
  edm::InputTag src_; // input collection

  JetCorrectionUncertainty* jecUncertainty_;

  bool isJECuncertaintyFromFile_;

  double shiftByJECuncertainty_;

  std::string jetCorrPayloadName_;
  std::string jetCorrUncertaintyTag_;

  SmearedParticleMaker<pat::Jet, GenJetRetriever<pat::Jet> >* smearingModule_;
};



 

