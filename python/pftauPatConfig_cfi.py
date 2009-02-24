import FWCore.ParameterSet.Config as cms
import copy

from PhysicsTools.PatAlgos.producersLayer1.tauProducer_cfi import * 

from PhysicsTools.PatAlgos.mcMatchLayer0.tauMatch_cfi import *
tauMatch.src = "allLayer0PFTausForTauAnalyses" # matching of recontructed tau-jets to generator level information
tauGenJetMatch.src = "allLayer0PFTausForTauAnalyses"

# computation of particle flow based IsoDeposits
from PhysicsTools.PatAlgos.recoLayer0.tauIsolation_cff import *
tauIsoDepositPFCandidatesForTauAnalyses = copy.deepcopy(tauIsoDepositPFCandidates)
tauIsoDepositPFCandidatesForTauAnalyses.src = cms.InputTag("pfRecoTauProducerForTauAnalyses")
#tauIsoFromDepsPFCandidatesForTauAnalyses = copy.deepcopy(tauIsoFromDepsPFCandidates)
#tauIsoFromDepsPFCandidatesForTauAnalyses.deposits[0].src = cms.InputTag("tauIsoDepositPFCandidatesForTauAnalyses")

tauIsoDepositPFChargedHadronsForTauAnalyses = copy.deepcopy(tauIsoDepositPFChargedHadrons)
tauIsoDepositPFChargedHadronsForTauAnalyses.src = cms.InputTag("pfRecoTauProducerForTauAnalyses")
#tauIsoFromDepsPFChargedHadronsForTauAnalyses = copy.deepcopy(tauIsoFromDepsPFChargedHadrons)
#tauIsoFromDepsPFChargedHadronsForTauAnalyses.deposits[0].src = cms.InputTag("tauIsoDepositPFChargedHadronsForTauAnalyses")

tauIsoDepositPFNeutralHadronsForTauAnalyses = copy.deepcopy(tauIsoDepositPFNeutralHadrons)
tauIsoDepositPFNeutralHadronsForTauAnalyses.src = cms.InputTag("pfRecoTauProducerForTauAnalyses")
#tauIsoFromDepsPFNeutralHadronsForTauAnalyses = copy.deepcopy(tauIsoFromDepsPFNeutralHadrons)
#tauIsoFromDepsPFNeutralHadronsForTauAnalyses.deposits[0].src = cms.InputTag("tauIsoDepositPFNeutralHadronsForTauAnalyses")

tauIsoDepositPFGammasForTauAnalyses = copy.deepcopy(tauIsoDepositPFGammas)
tauIsoDepositPFGammasForTauAnalyses.src = cms.InputTag("pfRecoTauProducerForTauAnalyses")
#tauIsoFromDepsPFGammasForTauAnalyses = copy.deepcopy(tauIsoFromDepsPFGammas)
#tauIsoFromDepsPFGammasForTauAnalyses.deposits[0].src = cms.InputTag("tauIsoDepositPFGammasForTauAnalyses")

patPFTauIsolationForTauAnalyses = cms.Sequence( tauIsoDepositPFCandidatesForTauAnalyses
                                               #*tauIsoFromDepsPFCandidatesForTauAnalyses
                                               *tauIsoDepositPFChargedHadronsForTauAnalyses
                                               #*tauIsoFromDepsPFChargedHadronsForTauAnalyses
                                               *tauIsoDepositPFNeutralHadronsForTauAnalyses
                                               #*tauIsoFromDepsPFNeutralHadronsForTauAnalyses
                                               *tauIsoDepositPFGammasForTauAnalyses
                                               #*tauIsoFromDepsPFGammasForTauAnalyses
                                              )
                                 
patAODPFTauIsolationLabelsForTauAnalyses = cms.VInputTag(
    cms.InputTag("tauIsoDepositPFCandidatesForTauAnalyses"),
    cms.InputTag("tauIsoDepositPFChargedHadronsForTauAnalyses"),
    cms.InputTag("tauIsoDepositPFNeutralHadronsForTauAnalyses"),
    cms.InputTag("tauIsoDepositPFGammasForTauAnalyses")
)

patAODPFTauIsolationForTauAnalyses = copy.deepcopy(patAODPFTauIsolation)
patAODPFTauIsolationForTauAnalyses.collection = cms.InputTag("pfRecoTauProducerForTauAnalyses")
patAODPFTauIsolationForTauAnalyses.associations = patAODPFTauIsolationLabelsForTauAnalyses

patLayer0PFTauIsolationForTauAnalyses = copy.deepcopy(patLayer0PFTauIsolation)
patLayer0PFTauIsolationForTauAnalyses.collection = cms.InputTag("allLayer0PFTausForTauAnalyses")
patLayer0PFTauIsolationForTauAnalyses.backrefs = cms.InputTag("allLayer0PFTausForTauAnalyses")
patLayer0PFTauIsolationForTauAnalyses.commonLabel = cms.InputTag("patAODPFTauIsolationForTauAnalyses")
patLayer0PFTauIsolationForTauAnalyses.associations = patAODPFTauIsolationLabelsForTauAnalyses

#####################  PAT LAYER 0  #########################

allLayer0PFTausForTauAnalyses = cms.EDFilter("PATPFTauCleaner",
     # EWK tau analysis specific tau reconstruction                                        
     tauSource              = cms.InputTag("pfRecoTauProducerForTauAnalyses"),
     tauDiscriminatorSource = cms.InputTag("pfRecoTauIsoDiscrForTauAnalyses"),
     # CMS "standard" tau reconstruction                                       
     #tauSource              = cms.InputTag("pfRecoTauProducer"),
     #tauDiscriminatorSource = cms.InputTag("pfRecoTauDiscriminationByIsolation"),
 
     removeOverlaps = cms.PSet(
       # Flag or discard taus that match with clean electrons
       #electrons = cms.PSet(     
       #    collection = cms.InputTag("allLayer0ElectronsForTauAnalyses"),
       #    deltaR     = cms.double(0.3)
       #),
       # Flag or discard taus that match with clean muons
       #muons = cms.PSet(     
       #    collection = cms.InputTag("allLayer0MuonsForTauAnalyses"),
       #    deltaR     = cms.double(0.3)
       #)
     ),

     markItems    = cms.bool(True),
     bitsToIgnore = cms.vstring('Core/Preselection'),
     saveRejected = cms.string(''),
     saveAll      = cms.string('')
)

#################### PAT LAYER 1 ##########################

allLayer1PFTausForTauAnalyses = copy.deepcopy(allLayer1Taus)
allLayer1PFTausForTauAnalyses.tauSource = "allLayer0PFTausForTauAnalyses"
allLayer1PFTausForTauAnalyses.tauIDSources = cms.PSet(
    leadingTrackFinding = cms.InputTag("pfRecoTauLdgTrkFindForTauAnalyses"),
    leadingTrackPtCut = cms.InputTag("pfRecoTauLdgTrkPtCutForTauAnalyses"),
    trackIsolation = cms.InputTag("pfRecoTauTrkIsoDiscrForTauAnalyses"),
    ecalIsolation = cms.InputTag("pfRecoTauEcalIsoDiscrForTauAnalyses"),
    #byIsolation = cms.InputTag("pfRecoTauIsoDiscrForTauAnalyses"),
    againstElectron = cms.InputTag("pfRecoTauElecRejDiscrForTauAnalyses"),
    againstMuon = cms.InputTag("pfRecoTauMuonRejDiscrForTauAnalyses")
)
allLayer1PFTausForTauAnalyses.isoDeposits.pfAllParticles = cms.InputTag("patLayer0PFTauIsolationForTauAnalyses", "tauIsoDepositPFCandidatesForTauAnalyses")
allLayer1PFTausForTauAnalyses.isoDeposits.pfChargedHadron = cms.InputTag("patLayer0PFTauIsolationForTauAnalyses", "tauIsoDepositPFChargedHadronsForTauAnalyses")
allLayer1PFTausForTauAnalyses.isoDeposits.pfNeutralHadron = cms.InputTag("patLayer0PFTauIsolationForTauAnalyses", "tauIsoDepositPFNeutralHadronsForTauAnalyses")
allLayer1PFTausForTauAnalyses.isoDeposits.pfGamma = cms.InputTag("patLayer0PFTauIsolationForTauAnalyses", "tauIsoDepositPFGammasForTauAnalyses")
allLayer1PFTausForTauAnalyses.isolation.pfAllParticles.src = allLayer1PFTausForTauAnalyses.isoDeposits.pfAllParticles
allLayer1PFTausForTauAnalyses.isolation.pfChargedHadron.src = allLayer1PFTausForTauAnalyses.isoDeposits.pfChargedHadron
allLayer1PFTausForTauAnalyses.isolation.pfNeutralHadron.src = allLayer1PFTausForTauAnalyses.isoDeposits.pfNeutralHadron
allLayer1PFTausForTauAnalyses.isolation.pfGamma.src = allLayer1PFTausForTauAnalyses.isoDeposits.pfGamma
allLayer1PFTausForTauAnalyses.addTrigMatch = False
allLayer1PFTausForTauAnalyses.addGenMatch = True
allLayer1PFTausForTauAnalyses.addGenJetMatch = True
