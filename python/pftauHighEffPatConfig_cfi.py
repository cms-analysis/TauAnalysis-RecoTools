import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.RecoTools.pftauPatConfig_cfi import *

#####################  PAT LAYER 0  #########################

allLayer0PFTausForTauAnalysesEff = copy.deepcopy(allLayer0PFTausForTauAnalyses)
allLayer0PFTausForTauAnalysesEff.tauSource              = "pfRecoTauProducerForTauAnalysesEff"
allLayer0PFTausForTauAnalysesEff.tauDiscriminatorSource = "pfRecoTauIsoDiscrForTauAnalysesEff"
allLayer0PFTausForTauAnalysesEff.bitsToIgnore = 'Core/Preselection'

#################### PAT LAYER 1 ##########################

allLayer1PFTausForTauAnalysesEff = copy.deepcopy(allLayer1PFTausForTauAnalyses)
allLayer1PFTausForTauAnalysesEff.tauSource = "allLayer0PFTausForTauAnalysesEff"
allLayer1PFTausForTauAnalysesEff.tauIDSources = cms.PSet(
    leadingTrackFinding = cms.InputTag("pfRecoTauLdgTrkFindForTauAnalysesEff"),
    leadingTrackPtCut = cms.InputTag("pfRecoTauLdgTrkPtCutForTauAnalysesEff"),
    trackIsolation = cms.InputTag("pfRecoTauTrkIsoDiscrForTauAnalysesEff"),
    ecalIsolation = cms.InputTag("pfRecoTauEcataulIsoDiscrForTauAnalysesEff"),
    #byIsolation = cms.InputTag("pfRecoTauIsoDiscrForTauAnalysesEff"),
    againstElectron = cms.InputTag("pfRecoTauElecRejDiscrForTauAnalysesEff"),
    againstMuon = cms.InputTag("pfRecoTauMuonRejDiscrForTauAnalysesEff")
)
