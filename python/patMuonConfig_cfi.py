import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.recoLayer0.muonIsolation_cff import *
from PhysicsTools.PatAlgos.triggerLayer0.trigMatchSequences_cff import *
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import *
from PhysicsTools.PatAlgos.cleaningLayer1.muonCleaner_cfi import *

#--------------------------------------------------------------------------------  
# PAT layer 0 muon configuration parameters
#--------------------------------------------------------------------------------

# add isolated HLT muon trigger to PAT trigger match sequence
patTrigMatch._seq = patTrigMatch._seq * patTrigMatchHLT1MuonIso

#--------------------------------------------------------------------------------  
# PAT layer 1 muon configuration parameters
#--------------------------------------------------------------------------------

# increase size of muon isolation cone from default of deltaR = 0.3 to 0.6
allLayer1Muons.isolation.tracker.deltaR = cms.double(0.6)

allLayer1Muons.isolation.ecal.deltaR = cms.double(0.6)

allLayer1Muons.isolation.hcal.deltaR = cms.double(0.6)

allLayer1Muons.isolation.user.deltaR = cms.double(0.6)
#
# add IsoDeposit objects for Track, ECAL and HCAL based isolation
#
allLayer1Muons.isoDeposits = cms.PSet(
   tracker         = allLayer1Muons.isolation.tracker.src,
   ecal            = allLayer1Muons.isolation.ecal.src,
   hcal            = allLayer1Muons.isolation.hcal.src,
   particle        = cms.InputTag("muonIsolationValueMap", "pfmuIsoDepositPFCandidates"),
   chargedparticle = cms.InputTag("muonIsolationValueMap", "pfmuIsoChDepositPFCandidates"),
   neutralparticle = cms.InputTag("muonIsolationValueMap", "pfmuIsoNeDepositPFCandidates"),
   gammaparticle   = cms.InputTag("muonIsolationValueMap", "pfmuIsoGaDepositPFCandidates")
)
#
# enable matching to HLT trigger information;
# match offline reconstructed muons to isolated and non-isolated HLT muon paths
#
allLayer1Muons.addTrigMatch = cms.bool(True)
allLayer1Muons.trigPrimMatch = cms.VInputTag(
    cms.InputTag("muonTrigMatchHLT1MuonNonIso"),
    cms.InputTag("muonTrigMatchHLT1MuonIso")
)
#
# enable matching to generator level information
#
allLayer1Muons.addGenMatch = cms.bool(True)
