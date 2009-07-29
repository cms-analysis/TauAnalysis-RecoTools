#include "TauAnalysis/RecoTools/plugins/PATElecTauPairZeeHypothesisEventSelector.h"

#include "PhysicsTools/UtilAlgos/interface/EventSelectorBase.h"

#include "FWCore/Framework/interface/MakerMacros.h"

using namespace reco;

DEFINE_EDM_PLUGIN(EventSelectorPluginFactory, PATElecTauPairZeeHypothesisMinEventSelector, "PATElecTauPairZeeHypothesisMinEventSelector");

DEFINE_EDM_PLUGIN(EventSelectorPluginFactory, PATElecTauPairZeeHypothesisMaxEventSelector, "PATElecTauPairZeeHypothesisMaxEventSelector");
