import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.RecoTools.electronIsolation_cff import *

inputTagEcalRecHitEB = cms.InputTag("reducedEcalRecHitsEB")
inputTagEcalRecHitEE = cms.InputTag("reducedEcalRecHitsEE")

EleIsoTrackExtractorBlock.barrelEcalHits = inputTagEcalRecHitEB
EleIsoTrackExtractorBlock.endcapEcalHits = inputTagEcalRecHitEE
eleIsoDepositForTauAnalysesTk.ExtractorPSet = cms.PSet(EleIsoTrackExtractorBlock)

EleIsoEcalFromHitsExtractorBlock.barrelRecHits = inputTagEcalRecHitEB
EleIsoEcalFromHitsExtractorBlock.endcapRecHits = inputTagEcalRecHitEE
EleIsoEcalFromHitsExtractorBlock.barrelEcalHits = inputTagEcalRecHitEB
EleIsoEcalFromHitsExtractorBlock.endcapEcalHits = inputTagEcalRecHitEE
eleIsoDepositForTauAnalysesEcalFromHits.ExtractorPSet = cms.PSet(EleIsoEcalFromHitsExtractorBlock)

EleIsoEcalFromClustsExtractorBlock.barrelEcalHits = inputTagEcalRecHitEB
EleIsoEcalFromClustsExtractorBlock.endcapEcalHits = inputTagEcalRecHitEE
eleIsoDepositForTauAnalysesEcalFromClusts.ExtractorPSet = cms.PSet(EleIsoEcalFromClustsExtractorBlock)

EleIsoEcalSCVetoFromClustsExtractorBlock.barrelEcalHits = inputTagEcalRecHitEB
EleIsoEcalSCVetoFromClustsExtractorBlock.endcapEcalHits = inputTagEcalRecHitEE
eleIsoDepositForTauAnalysesEcalSCVetoFromClusts.ExtractorPSet = cms.PSet(EleIsoEcalSCVetoFromClustsExtractorBlock)

EleIsoHcalFromHitsExtractorBlock.barrelEcalHits = inputTagEcalRecHitEB
EleIsoHcalFromHitsExtractorBlock.endcapEcalHits = inputTagEcalRecHitEE
eleIsoDepositForTauAnalysesHcalFromHits.ExtractorPSet = cms.PSet(EleIsoHcalFromHitsExtractorBlock)

EleIsoHcalFromTowersExtractorBlock.barrelEcalHits = inputTagEcalRecHitEB
EleIsoHcalFromTowersExtractorBlock.endcapEcalHits = inputTagEcalRecHitEE
eleIsoDepositForTauAnalysesHcalFromTowers.ExtractorPSet = cms.PSet(EleIsoHcalFromTowersExtractorBlock)

# HCAL recHits don't exist on AOD (and hence in the output files produced by the FastSimulation);
# as a consequence, eleIsoDepositForTauAnalysesHcalFromHits needs to be removed from the sequence
patAODElectronIsolationLabels = cms.VInputTag(
  cms.InputTag("eleIsoDepositForTauAnalysesTk"),
  cms.InputTag("eleIsoDepositForTauAnalysesEcalFromHits"),
  cms.InputTag("eleIsoDepositForTauAnalysesEcalFromClusts"),
  cms.InputTag("eleIsoDepositForTauAnalysesHcalFromTowers"),
  cms.InputTag("eleIsoDepositForTauAnalysesEcalSCVetoFromClusts")
)

patAODElectronIsolations = cms.EDFilter("MultipleIsoDepositsToValueMaps",
  collection   = cms.InputTag("pixelMatchGsfElectrons"),
  associations = patAODElectronIsolationLabels,
)

layer0ElectronIsolations = cms.EDFilter("CandManyValueMapsSkimmerIsoDeposits",
  collection   = cms.InputTag("allLayer0ElectronsForTauAnalyses"),
  backrefs     = cms.InputTag("allLayer0ElectronsForTauAnalyses"),
  commonLabel  = cms.InputTag("patAODElectronIsolations"),
  associations = patAODElectronIsolationLabels
)

eleIsoDepositAOD = cms.Sequence(
   eleIsoDepositForTauAnalysesTk * eleIsoDepositForTauAnalysesEcalFromClusts * eleIsoDepositForTauAnalysesEcalFromHits 
  *eleIsoDepositForTauAnalysesHcalFromTowers 
  *eleIsoDepositForTauAnalysesEcalSCVetoFromClusts
)

patAODElectronIsolation = cms.Sequence(
   egammaSuperClusterMerger 
  *egammaBasicClusterMerger
  *eleIsoDepositAOD
  *patAODElectronIsolations
)

patLayer0ElectronIsolation = cms.Sequence(layer0ElectronIsolations)
