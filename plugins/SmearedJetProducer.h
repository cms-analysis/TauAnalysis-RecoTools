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

//---------------------- member data ---------------------------
  std::string moduleLabel_;

  edm::InputTag src_; // input collection

  JetCorrectionUncertainty* jecUncertainty_;

  double jecUncertaintyValue_;

  bool isJECuncertaintyFromFile_;

  double shiftByJECuncertainty_;

  std::string jetCorrPayloadName_;
  std::string jetCorrUncertaintyTag_;

  double varyByNsigmas_; // option to shift energy by N standard-deviations, useful for template morphing

  SmearedParticleMaker<pat::Jet, GenJetRetriever<pat::Jet> >* smearingModule_;
};



 

