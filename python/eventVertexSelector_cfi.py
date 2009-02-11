import FWCore.ParameterSet.Config as cms
import copy

# module to select candidates for "the primary event vertex".
# Selection based on PhysicsTools/PatAlogos/plugins/PATSingleVertexSelector;
# see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePhysicsCutParser
# on how to use the cut-string

# Note: offlinePrimaryVerticesWithBS collection is sorted
#       in order of decreasing sum of Pt of tracks fitted to each vertex

selectedPrimaryVertexHighestPtTrackSum = cms.EDFilter("PATSingleVertexSelector",
    mode = cms.string('firstVertex'),
    vertices = cms.InputTag('offlinePrimaryVerticesWithBS'),
    filter = cms.bool(False)                                                    
)

selectedPrimaryVertexQuality = cms.EDFilter("VertexSelector",
    src = cms.InputTag('selectedPrimaryVertexHighestPtTrackSum'),
    cut = cms.string("isValid & (chi2prob(chi2,ndof) > 0.01) & (tracksSize >= 2)"),
    filter = cms.bool(False)                                          
)

selectedPrimaryVertexPosition = cms.EDFilter("VertexSelector",
    src = cms.InputTag('selectedPrimaryVertexQuality'),
    cut = cms.string("z > -25 & z < +25"),
    filter = cms.bool(False)                                           
)

selectPrimaryVertexForTauAnalyses = cms.Sequence( selectedPrimaryVertexHighestPtTrackSum
                                                 *selectedPrimaryVertexQuality
                                                 *selectedPrimaryVertexPosition )
