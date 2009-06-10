#include "TauAnalysis/CandidateTools/interface/CompositePtrCandidateT1T2MEtProducer.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

typedef CompositePtrCandidateT1T2MEtProducer<pat::Electron, pat::Tau> PATElecTauPairProducer;
typedef CompositePtrCandidateT1T2MEtProducer<pat::Muon, pat::Tau> PATMuTauPairProducer;
typedef CompositePtrCandidateT1T2MEtProducer<pat::Tau, pat::Tau> PATDiTauPairProducer;
typedef CompositePtrCandidateT1T2MEtProducer<pat::Electron, pat::Muon> PATElecMuPairProducer;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_ANOTHER_FWK_MODULE(PATElecTauPairProducer);
DEFINE_ANOTHER_FWK_MODULE(PATMuTauPairProducer);
DEFINE_ANOTHER_FWK_MODULE(PATDiTauPairProducer);
DEFINE_ANOTHER_FWK_MODULE(PATElecMuPairProducer);
