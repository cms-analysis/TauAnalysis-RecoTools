import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.RecoTools.pftauHighEffReconstruction_cff import *
from TauAnalysis.RecoTools.pftauHighEffPatConfig_cfi import *
from TauAnalysis.RecoTools.pftauHighEffPatSelector_cfi import *

allLayer1PFTausSelForTauAnalysesEff = cms.EDProducer("PATTauSelProducer",

  leptonSource = cms.InputTag("allLayer1PFTausForTauAnalysesEff"),

  selFlags = cms.PSet(
    tauAnalysisSelTauEta21 = cms.PSet(
      src = cms.InputTag('selectedLayer1TausEta21IndividualEff')
    ),
    tauAnalysisSelTauEt20 = cms.PSet(
      src = cms.InputTag('selectedLayer1TausEt20IndividualEff')
    ),
    tauAnalysisSelTauLeadTrk = cms.PSet(
      src = cms.InputTag('selectedLayer1TausLeadTrkIndividualEff')
    ),
    tauAnalysisSelTauLeadTrkPt = cms.PSet(
      src = cms.InputTag('selectedLayer1TausLeadTrkPtIndividualEff')
    ),
    tauAnalysisSelTauTrkIso = cms.PSet(
      src = cms.InputTag('selectedLayer1TausTrkIsoIndividualEff')
    ),
    tauAnalysisSelTauEcalIso = cms.PSet(
      src = cms.InputTag('selectedLayer1TausEcalIsoIndividualEff')
    ),
    tauAnalysisSelTauElecVeto = cms.PSet(
      src = cms.InputTag('selectedLayer1TausElecVetoIndividualEff')
    ),
    tauAnalysisSelTauMuonVeto = cms.PSet(
      src = cms.InputTag('selectedLayer1TausMuonVetoIndividualEff')
    ),
    tauAnalysisSelTauProng = cms.PSet(
      src = cms.InputTag('selectedLayer1TausProngIndividualEff')
    )
  )
)

produceTausForTauAnalysesEff = cms.Sequence( pfRecoTauIsoDiscrForTauAnalysesEff
                                            *allLayer0PFTausForTauAnalysesEff
                                            *pfRecoTauLdgTrkFindForTauAnalysesEff
                                            *pfRecoTauLdgTrkPtCutForTauAnalysesEff
                                            *pfRecoTauTrkIsoDiscrForTauAnalysesEff
                                            *pfRecoTauEcalIsoDiscrForTauAnalysesEff
                                            *pfRecoTauElecRejDiscrForTauAnalysesEff
                                            *pfRecoTauMuonRejDiscrForTauAnalysesEff
                                            *allLayer1PFTausForTauAnalysesEff
                                            *selectPFTausForTauAnalysesEff
                                            *allLayer1PFTausSelForTauAnalysesEff )
                                   


