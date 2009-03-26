import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------  
# reconstruct IsoDeposits for electrons in an enlarged isolation cone 
#--------------------------------------------------------------------------------

from RecoEgamma.EgammaIsolationAlgos.eleTrackExtractorBlocks_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleEcalExtractorBlocks_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleHcalExtractorBlocks_cff import *

from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositTk_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositEcalFromHits_cff import *

#--------------------------------------------------------------------------------
#
# CV: eleIsoDepositHcalFromTowers module already included by
#      PhysicsTools/PatAlgos/python/recoLayer0/electronIsolation_cff.py
#     in CMSSW_2_2_x
#    --> causes run-time error
#
#     %MSG-s ConfigFileReadError:  26-Mar-2009 09:29:00 CET pre-events
#     Problem with configuration file runZtoMuTau_cfg.py
#     ---- Configuration BEGIN
#     python encountered the error: exceptions.AttributeError 'EDProducer' object has no attribute '_Labelable__label'
#     ---- Configuration END
#
#        if module gets included twice.
#
# NOTE: include temporarily disabled for CMSSW_2_2_x only !!
#
#from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositHcalFromTowers_cff import *
#--------------------------------------------------------------------------------

#
# Track based isolation
#
EleIsoTrackExtractorBlock.DR_Max = cms.double(1.)

eleIsoDepositTk.ExtractorPSet = cms.PSet(EleIsoTrackExtractorBlock)

#
# ECAL recHit based isolation
#
EleIsoEcalFromHitsExtractorBlock.extRadius = cms.double(1.)

eleIsoDepositEcalFromHits.ExtractorPSet = cms.PSet(EleIsoEcalFromHitsExtractorBlock)

#
# HCAL isolation based on CaloTowers
#
EleIsoHcalFromTowersExtractorBlock.extRadius = cms.double(1.)

#eleIsoDepositHcalFromTowers.ExtractorPSet = cms.PSet(EleIsoHcalFromTowersExtractorBlock)

electronIsoDeposits = cms.Sequence( eleIsoDepositEcalFromHits )
#--------------------------------------------------------------------------------
#
# CV: eleIsoDepositTk module already executed by sequence eleIsoDepositAOD
#     defined by
#      PhysicsTools/PatAlgos/python/recoLayer0/electronIsolation_cff.py
#     in CMSSW_2_2_x
#    --> causes run-time error
#
#     <ValueMap::rawIndexOf>:
#      id = 2018
#      idx = 0
#     %MSG-w ScheduleExecutionFailure:  PostModule 26-Mar-2009 10:15:47 CET  Run: 1 Event: 799
#     an exception occurred and all paths for the event are being skipped: 
#     ---- ScheduleExecutionFailure BEGIN
#     ProcessingStopped
#     ---- InvalidReference BEGIN
#     ValueMap: no associated value for given product and index
#     cms::Exception going through module PATElectronProducer/allLayer1Electrons run: 1 event: 799
#     ---- InvalidReference END
#     Exception going through path p
#     ---- ScheduleExecutionFailure END
#
#        if module gets executed twice.
#electronIsoDeposits = cms.Sequence( eleIsoDepositTk
#                                   *eleIsoDepositEcalFromHits 
#                                   *eleIsoDepositHcalFromTowers )
#--------------------------------------------------------------------------------
