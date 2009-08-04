//
// $Id: VertexEventSelector.h,v 1.1 2009/02/11 17:10:16 veelken Exp $
//

#ifndef TauAnalysis_RecoTools_VertexEventSelector_h
#define TauAnalysis_RecoTools_VertexEventSelector_h

#include "PhysicsTools/UtilAlgos/interface/AnySelector.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectCountEventSelector.h"
#include "PhysicsTools/UtilAlgos/interface/MinNumberSelector.h"
#include "PhysicsTools/UtilAlgos/interface/MaxNumberSelector.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include <vector>

namespace reco {

  typedef ObjectCountEventSelector<std::vector<Vertex>, AnySelector, MinNumberSelector> VertexMinEventSelector;

  typedef ObjectCountEventSelector<std::vector<Vertex>, AnySelector, MaxNumberSelector> VertexMaxEventSelector;

}

#endif
