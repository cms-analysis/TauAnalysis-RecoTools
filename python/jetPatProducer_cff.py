from PhysicsTools.PatAlgos.cleaningLayer0.caloJetCleaner_cfi import *
from PhysicsTools.PatAlgos.recoLayer0.jetFlavourId_cff import *
from PhysicsTools.PatAlgos.recoLayer0.bTagging_cff import *
from PhysicsTools.PatAlgos.recoLayer0.jetTracksCharge_cff import *
from PhysicsTools.PatAlgos.recoLayer0.jetMETCorrections_cff import *
from PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi import *
from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff import *
allLayer0Jets.removeOverlaps.electrons.collection = cms.InputTag("allLayer0ElectronsForTauAnalyses")
allLayer1Jets.addTrigMatch = cms.bool(False)
allLayer1Jets.jetCorrFactorsSource = cms.InputTag("layer0JetCorrFactors")

patLayer0JetsForTauAnalyses = cms.Sequence( allLayer0Jets
                                           *patAODBTagging
                                           *patJetFlavourId
                                           *patLayer0BTagging
                                           #*patAODJetMETCorrections --> no need to do this here, as it is already done in met
                                           *patLayer0JetMETCorrections
                                           *patLayer0JetTracksCharge
                                           *(jetPartonMatch + jetGenJetMatch) )

patLayer1JetsForTauAnalyses = cms.Sequence(allLayer1Jets)

produceJetsForTauAnalyses = cms.Sequence(patLayer0JetsForTauAnalyses + patLayer1JetsForTauAnalyses)
