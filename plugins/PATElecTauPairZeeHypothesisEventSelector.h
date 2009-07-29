//
// $Id: PATElecTauPairZeeHypothesisEventSelector.h,v 1.1 2009/02/11 17:10:16 veelken Exp $
//

#ifndef TauAnalysis_RecoTools_PATElecTauPairZeeHypothesisEventSelector_h
#define TauAnalysis_RecoTools_PATElecTauPairZeeHypothesisEventSelector_h

#include "PhysicsTools/UtilAlgos/interface/AnySelector.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectCountEventSelector.h"
#include "PhysicsTools/UtilAlgos/interface/MinNumberSelector.h"
#include "PhysicsTools/PatUtils/interface/MaxNumberSelector.h"

#include "AnalysisDataFormats/TauAnalysis/interface/PATElecTauPairZeeHypothesis.h"

#include <vector>

namespace reco {

  typedef ObjectCountEventSelector<std::vector<PATElecTauPairZeeHypothesis>, AnySelector, MinNumberSelector> PATElecTauPairZeeHypothesisMinEventSelector;

  typedef ObjectCountEventSelector<std::vector<PATElecTauPairZeeHypothesis>, AnySelector, MaxNumberSelector> PATElecTauPairZeeHypothesisMaxEventSelector;

}

#endif
