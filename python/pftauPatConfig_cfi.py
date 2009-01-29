import FWCore.ParameterSet.Config as cms
import copy

from PhysicsTools.PatAlgos.producersLayer1.tauProducer_cfi import * 

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
allLayer1PFTausForTauAnalyses.addTrigMatch = False
allLayer1PFTausForTauAnalyses.addGenMatch = False
allLayer1PFTausForTauAnalyses.addGenJetMatch = False
