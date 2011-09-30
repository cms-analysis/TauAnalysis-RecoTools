#include "TauAnalysis/RecoTools/interface/JetMETsmearInputProducerT.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"

typedef JetMETsmearInputProducerT<reco::CaloJet> CaloJetMETsmearInputProducer;
typedef JetMETsmearInputProducerT<reco::PFJet> PFJetMETsmearInputProducer;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(CaloJetMETsmearInputProducer);
DEFINE_FWK_MODULE(PFJetMETsmearInputProducer);
