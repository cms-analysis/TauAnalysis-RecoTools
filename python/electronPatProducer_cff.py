from TauAnalysis.RecoTools.electronPatConfig_cfi import *
from TauAnalysis.RecoTools.electronPatSelector_cfi import *
from TauAnalysis.RecoTools.electronAntiIsoPatSelector_cfi import *
from TauAnalysis.RecoTools.electronIdentification_cff import *
from TauAnalysis.RecoTools.electronIsolation_cff import *

from PhysicsTools.PatAlgos.triggerLayer0.trigMatchSequences_cff import *
electronTrigMatchCandHLT1ElectronStartup.src = "allLayer0ElectronsForTauAnalyses"
electronTrigMatchHLT1Electron.src = "allLayer0ElectronsForTauAnalyses"
electronTrigMatchHLT1ElectronRelaxed.src = "allLayer0ElectronsForTauAnalyses"

patTrigMatchElectronsForTauAnalyses = cms.Sequence( (patCandHLT1ElectronStartup * electronTrigMatchCandHLT1ElectronStartup)
	                                           +(patHLT1ElectronRelaxed * electronTrigMatchHLT1ElectronRelaxed)
                                                   +(patHLT1Electron * electronTrigMatchHLT1Electron) )

from PhysicsTools.PatAlgos.mcMatchLayer0.electronMatch_cfi import *
electronMatch.src = "allLayer0ElectronsForTauAnalyses" # matching of recontructed electrons to generator level information

patLayer0ElectronsForTauAnalyses = cms.Sequence( patAODElectronIsolation
                                                +allLayer0ElectronsForTauAnalyses
                                                +patLayer0ElectronId
                                                +patLayer0ElectronIsolation
                                                +electronMatch
                                                +patTrigMatchElectronsForTauAnalyses )
patLayer1ElectronsForTauAnalyses = cms.Sequence(allLayer1ElectronsForTauAnalyses)

allLayer1ElectronsSelForTauAnalyses = cms.EDProducer("PATElectronSelProducer",

  leptonSource = cms.InputTag("allLayer1ElectronsForTauAnalyses"),

  selFlags = cms.PSet(
    tauAnalysisSelElectronTightIdGlobal = cms.PSet(
      src = cms.InputTag('selectedLayer1ElectronsTightId')
    ),
    tauAnalysisSelElectronAntiCrackCut = cms.PSet(
      src = cms.InputTag('selectedLayer1ElectronsAntiCrackCutIndividual')
    ),
    tauAnalysisSelElectronEta21 = cms.PSet(
      src = cms.InputTag('selectedLayer1ElectronsEta21Individual')
    ),
    tauAnalysisSelElectronPt15 = cms.PSet(
      src = cms.InputTag('selectedLayer1ElectronsPt15Individual')
    ),
    tauAnalysisSelElectronHLTmatch = cms.PSet(
      src = cms.InputTag('selectedLayer1ElectronsHLTmatchIndividual')
    ),
    tauAnalysisSelElectronTrkIso = cms.PSet(
      src = cms.InputTag('selectedLayer1ElectronsTrkIsoIndividual')
    ),
    tauAnalysisSelElectronEcalIso = cms.PSet(
      src = cms.InputTag('selectedLayer1ElectronsEcalIsoIndividual')
    ),
    tauAnalysisSelElectronHcalIso = cms.PSet(
      src = cms.InputTag('selectedLayer1ElectronsHcalIsoIndividual')
    ),
    tauAnalysisSelElectronTrk = cms.PSet(
      src = cms.InputTag('selectedLayer1ElectronsTrkIndividual')
    ),
    tauAnalysisSelElectronTrkIP = cms.PSet(
      src = cms.InputTag('selectedLayer1ElectronsTrkIPindividual')
    )
  )
)

produceElectronsForTauAnalyses = cms.Sequence( patLayer0ElectronsForTauAnalyses
                                              +patLayer1ElectronsForTauAnalyses
                                              +selectElectronsForTauAnalyses
                                              +allLayer1ElectronsSelForTauAnalyses )





                                         
