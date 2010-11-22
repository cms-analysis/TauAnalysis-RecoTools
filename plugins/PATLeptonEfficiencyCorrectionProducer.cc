#include "TauAnalysis/RecoTools/plugins/PATLeptonEfficiencyCorrectionProducer.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

typedef PATLeptonEfficiencyCorrectionProducer<pat::Electron> PATElectronEfficiencyCorrectionProducer;
typedef PATLeptonEfficiencyCorrectionProducer<pat::Muon> PATMuonEfficiencyCorrectionProducer;
typedef PATLeptonEfficiencyCorrectionProducer<pat::Tau> PATTauEfficiencyCorrectionProducer;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATElectronEfficiencyCorrectionProducer);
DEFINE_FWK_MODULE(PATMuonEfficiencyCorrectionProducer);
DEFINE_FWK_MODULE(PATTauEfficiencyCorrectionProducer);

