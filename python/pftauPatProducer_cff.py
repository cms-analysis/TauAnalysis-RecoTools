import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.RecoTools.pftauReconstruction_cff import *
from TauAnalysis.RecoTools.pftauPatConfig_cfi import *
from TauAnalysis.RecoTools.pftauPatSelector_cfi import *
from TauAnalysis.RecoTools.pftauPatSelectorForElecTau_cfi import *
from TauAnalysis.RecoTools.pftauPatSelectorForMuTau_cfi import *

allLayer1PFTausSelForTauAnalyses = cms.EDProducer("PATTauSelProducer",

  leptonSource = cms.InputTag("allLayer1PFTausForTauAnalyses"),

  selFlags = cms.PSet(
    tauAnalysisSelTauEta21 = cms.PSet(
      src = cms.InputTag('selectedLayer1TausEta21Individual')
    ),
    tauAnalysisSelTauPt20 = cms.PSet(
      src = cms.InputTag('selectedLayer1TausPt20Individual')
    ),
    tauAnalysisSelTauLeadTrk = cms.PSet(
      src = cms.InputTag('selectedLayer1TausLeadTrkIndividual')
    ),
    tauAnalysisSelTauLeadTrkPt = cms.PSet(
      src = cms.InputTag('selectedLayer1TausLeadTrkPtIndividual')
    ),
    tauAnalysisSelTauTrkIso = cms.PSet(
      src = cms.InputTag('selectedLayer1TausTrkIsoIndividual')
    ),
    tauAnalysisSelTauEcalIso = cms.PSet(
      src = cms.InputTag('selectedLayer1TausEcalIsoIndividual')
    ),
    tauAnalysisSelTauElecVeto = cms.PSet(
      src = cms.InputTag('selectedLayer1TausElectronVetoIndividual')
    ),
    tauAnalysisSelTauMuonVeto = cms.PSet(
      src = cms.InputTag('selectedLayer1TausMuonVetoIndividual')
    ),
    tauAnalysisSelTauProng = cms.PSet(
      src = cms.InputTag('selectedLayer1TausProngIndividual')
    )
  )
)

producePFTausForTauAnalyses = cms.Sequence( pfRecoTauProducerForTauAnalyses
                                           *pfRecoTauIsoDiscrForTauAnalyses
                                           *allLayer0PFTausForTauAnalyses
                                           *pfRecoTauLdgTrkFindForTauAnalyses
                                           *pfRecoTauLdgTrkPtCutForTauAnalyses
                                           *pfRecoTauTrkIsoDiscrForTauAnalyses
                                           *pfRecoTauEcalIsoDiscrForTauAnalyses
                                           *pfRecoTauElecRejDiscrForTauAnalyses
                                           *pfRecoTauMuonRejDiscrForTauAnalyses
                                           *allLayer1PFTausForTauAnalyses
                                           *selectPFTausForTauAnalyses
                                           *selectPFTausForElecTau
                                           *selectPFTausForMuTau
                                           *allLayer1PFTausSelForTauAnalyses )




