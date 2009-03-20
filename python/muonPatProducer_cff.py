from TauAnalysis.RecoTools.muonIsolation_cff import *
from TauAnalysis.RecoTools.muonPatConfig_cfi import *
from TauAnalysis.RecoTools.muonPatSelector_cfi import *

from PhysicsTools.PatAlgos.triggerLayer0.trigMatchSequences_cff import *
muonTrigMatchHLT1MuonIso.src = "allLayer0MuonsForTauAnalyses"
muonTrigMatchHLT1MuonNonIso.src = "allLayer0MuonsForTauAnalyses"
patTrigMatchMuonsForTauAnalyses = cms.Sequence( (patHLT1MuonIso * muonTrigMatchHLT1MuonIso)
	                                       +(patHLT1MuonNonIso * muonTrigMatchHLT1MuonNonIso) )

from PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi import *
muonMatch.src = "allLayer0MuonsForTauAnalyses" # matching of recontructed muons to generator level information

patLayer0MuonsForTauAnalyses = cms.Sequence( patAODMuonIsolation
                                            +allLayer0MuonsForTauAnalyses
                                            +patLayer0MuonIsolation
                                            +muonMatch
                                            +patTrigMatchMuonsForTauAnalyses )
patLayer1MuonsForTauAnalyses = cms.Sequence(allLayer1MuonsForTauAnalyses)

allLayer1MuonsSelForTauAnalyses = cms.EDProducer("PATMuonSelProducer",

  leptonSource = cms.InputTag("allLayer1MuonsForTauAnalyses"),

  selFlags = cms.PSet(
    tauAnalysisSelMuonGlobal = cms.PSet(
      src = cms.InputTag('selectedLayer1MuonsGlobal')
    ),
    tauAnalysisSelMuonEta21 = cms.PSet(
      src = cms.InputTag('selectedLayer1MuonsEta21Individual')
    ),
    tauAnalysisSelMuonPt15 = cms.PSet(
      src = cms.InputTag('selectedLayer1MuonsPt15Individual')
    ),
    tauAnalysisSelMuonTrkIso = cms.PSet(
      src = cms.InputTag('selectedLayer1MuonsTrkIsoIndividual')
    ),
    tauAnalysisSelMuonEcalIso = cms.PSet(
      src = cms.InputTag('selectedLayer1MuonsEcalIsoIndividual')
    ),
    tauAnalysisSelMuonPionVeto = cms.PSet(
      src = cms.InputTag('selectedLayer1MuonsPionVetoIndividual')
    ),
    tauAnalysisSelMuonTrkIP = cms.PSet(
      src = cms.InputTag('selectedLayer1MuonsTrkIPindividual')
    )
  )
)

produceMuonsForTauAnalyses = cms.Sequence( patLayer0MuonsForTauAnalyses
                                          +patLayer1MuonsForTauAnalyses
                                          +selectMuonsForTauAnalyses
                                          +allLayer1MuonsSelForTauAnalyses )
