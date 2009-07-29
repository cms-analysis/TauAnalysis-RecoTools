//
// $Id: PATElecTauPairZeeHypothesisSelector.h,v 1.1 2009/06/10 09:33:09 veelken Exp $
//

#ifndef TauAnalysis_RecoTools_PATElecTauPairZeeHypothesisSelector_h
#define TauAnalysis_RecoTools_PATElecTauPairZeeHypothesisSelector_h

#include "DataFormats/Common/interface/RefVector.h"

#include "PhysicsTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "PhysicsTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"
#include "PhysicsTools/UtilAlgos/interface/SingleElementCollectionSelector.h"

#include "AnalysisDataFormats/TauAnalysis/interface/PATElecTauPairZeeHypothesis.h"

#include <vector>

typedef SingleObjectSelector<
            std::vector<PATElecTauPairZeeHypothesis>,
            StringCutObjectSelector<PATElecTauPairZeeHypothesis>
        > PATElecTauPairZeeHypothesisSelector;

typedef SingleObjectSelector<
            std::vector<PATElecTauPairZeeHypothesis>,
            StringCutObjectSelector<PATElecTauPairZeeHypothesis>,
            edm::RefVector<std::vector<PATElecTauPairZeeHypothesis> >
        > PATElecTauPairZeeHypothesisRefSelector;

#endif
