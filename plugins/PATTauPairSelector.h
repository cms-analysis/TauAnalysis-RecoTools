//
// $Id: PATTauPairSelector.h,v 1.1 2009/02/04 17:32:15 veelken Exp $
//

#ifndef TauAnalysis_RecoTools_PATTauPairSelector_h
#define TauAnalysis_RecoTools_PATTauPairSelector_h

#include "DataFormats/Common/interface/RefVector.h"

#include "PhysicsTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "PhysicsTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"
#include "PhysicsTools/UtilAlgos/interface/SingleElementCollectionSelector.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"

#include <vector>

typedef SingleObjectSelector<
            std::vector<PATElecTauPair>,
            StringCutObjectSelector<PATElecTauPair>
        > PATElecTauPairSelector;
typedef SingleObjectSelector<
            std::vector<PATMuTauPair>,
            StringCutObjectSelector<PATMuTauPair>
        > PATMuTauPairSelector;
typedef SingleObjectSelector<
            std::vector<PATDiTauPair>,
            StringCutObjectSelector<PATDiTauPair>
        > PATDiTauPairSelector;
typedef SingleObjectSelector<
            std::vector<PATElecMuPair>,
            StringCutObjectSelector<PATElecMuPair>
        > PATElecMuPairSelector;

typedef SingleObjectSelector<
            std::vector<PATElecTauPair>,
            StringCutObjectSelector<PATElecTauPair>,
            edm::RefVector<std::vector<PATElecTauPair> >
        > PATElecTauPairRefSelector;
typedef SingleObjectSelector<
            std::vector<PATMuTauPair>,
            StringCutObjectSelector<PATMuTauPair>,
            edm::RefVector<std::vector<PATMuTauPair> >
        > PATMuTauPairRefSelector;
typedef SingleObjectSelector<
            std::vector<PATDiTauPair>,
            StringCutObjectSelector<PATDiTauPair>,
            edm::RefVector<std::vector<PATDiTauPair> >
        > PATDiTauPairRefSelector;
typedef SingleObjectSelector<
            std::vector<PATElecMuPair>,
            StringCutObjectSelector<PATElecMuPair>,
            edm::RefVector<std::vector<PATElecMuPair> >
        > PATElecMuPairRefSelector;

#endif
