#include "TauAnalysis/RecoTools/interface/JetMETsmearInputProducerT.h"

namespace JetMETsmearInputProducer_namespace
{
  template <>
  class GenJetMatcherT<pat::Jet>
  {
    public:

     GenJetMatcherT(const edm::ParameterSet&) {}
     ~GenJetMatcherT() {}

     const reco::GenJet* operator()(const pat::Jet& jet, edm::Event* evt = 0) const
     {
       return jet.genJet();
     }
  };
}

typedef JetMETsmearInputProducerT<pat::Jet> PATJetMETsmearInputProducer;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATJetMETsmearInputProducer);
