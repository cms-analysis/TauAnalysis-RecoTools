import FWCore.ParameterSet.Config as cms
import copy

from PhysicsTools.PatAlgos.cleaningLayer0.electronCleaner_cfi import *
from PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi import *

allLayer0ElectronsForTauAnalyses = copy.deepcopy(allLayer0Electrons)
allLayer0ElectronsForTauAnalyses.electronSource = cms.InputTag("pixelMatchGsfElectrons")
allLayer0ElectronsForTauAnalyses.removeDuplicates = cms.bool(True)
allLayer0ElectronsForTauAnalyses.selection = cms.PSet(type = cms.string('none'))
allLayer0ElectronsForTauAnalyses.isolation = cms.PSet()
allLayer0ElectronsForTauAnalyses.markItems = cms.bool(True)
allLayer0ElectronsForTauAnalyses.bitsToIgnore = cms.vstring('')
allLayer0ElectronsForTauAnalyses.saveRejected = cms.string('')
allLayer0ElectronsForTauAnalyses.saveAll = cms.string('')
     
allLayer1ElectronsForTauAnalyses = copy.deepcopy(allLayer1Electrons)
allLayer1ElectronsForTauAnalyses.electronSource = cms.InputTag("allLayer0ElectronsForTauAnalyses")
allLayer1ElectronsForTauAnalyses.embedTrack = cms.bool(False)
allLayer1ElectronsForTauAnalyses.embedGsfTrack = cms.bool(False)
allLayer1ElectronsForTauAnalyses.embedSuperCluster = cms.bool(False)
allLayer1ElectronsForTauAnalyses.isolation = cms.PSet(
   tracker = cms.PSet(
      src = cms.InputTag("layer0ElectronIsolations","eleIsoDepositForTauAnalysesTk"),
      deltaR = cms.double(0.6),
      vetos = cms.vstring(
         '0.02', # inner radius veto cone (was 0.015)
         'Threshold(1.0)' # threshold on individual track pt
      ),       
      skipDefaultVeto = cms.bool(True),
   ),
   ecal = cms.PSet(
      src = cms.InputTag("layer0ElectronIsolations","eleIsoDepositForTauAnalysesEcalFromHits"),
      deltaR = cms.double(0.6),
      vetos = cms.vstring(
         'EcalBarrel:0.040', 
         'EcalBarrel:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)',
         'EcalEndcaps:0.070',
         'EcalEndcaps:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)',
         'EcalBarrel:ThresholdFromTransverse(0.08)',
         'EcalEndcaps:ThresholdFromTransverse(0.3)'
      ),
      skipDefaultVeto = cms.bool(True),
   ),
   hcal = cms.PSet(
      src = cms.InputTag("layer0ElectronIsolations","eleIsoDepositForTauAnalysesHcalFromTowers"),
      deltaR = cms.double(0.6),
      skipDefaultVeto = cms.bool(True),
   ),
)
allLayer1ElectronsForTauAnalyses.isoDeposits = cms.PSet(
   tracker = cms.InputTag("layer0ElectronIsolations","eleIsoDepositForTauAnalysesTk"),
   ecal    = cms.InputTag("layer0ElectronIsolations","eleIsoDepositForTauAnalysesEcalFromHits"),
   hcal    = cms.InputTag("layer0ElectronIsolations","eleIsoDepositForTauAnalysesHcalFromTowers"),
)
allLayer1ElectronsForTauAnalyses.addElectronID = cms.bool(True)
allLayer1ElectronsForTauAnalyses.electronIDSources = cms.PSet(
   robust     = cms.InputTag("elecIdForTauAnalysesCutBasedRobust"),
   loose      = cms.InputTag("elecIdForTauAnalysesCutBasedLoose"),
   tight      = cms.InputTag("elecIdForTauAnalysesCutBasedTight"),        
)
allLayer1ElectronsForTauAnalyses.addTrigMatch = cms.bool(True)
allLayer1ElectronsForTauAnalyses.trigPrimMatch = cms.VInputTag(
   cms.InputTag("electronTrigMatchHLT1Electron"),
   #cms.InputTag("electronTrigMatchHLT1ElectronRelaxed"), 
)
allLayer1ElectronsForTauAnalyses.addGenMatch = cms.bool(False)
allLayer1ElectronsForTauAnalyses.embedGenMatch    = cms.bool(False)
allLayer1ElectronsForTauAnalyses.genParticleMatch = cms.InputTag("electronMatch")

