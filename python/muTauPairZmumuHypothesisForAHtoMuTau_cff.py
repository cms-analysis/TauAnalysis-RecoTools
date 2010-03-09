import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.RecoTools.muTauPairZmumuHypothesis_cff import *

#--------------------------------------------------------------------------------
# produce data-formats providing information about compatibility of muon + tau-jet pairs
# with hypothesis of being a pair of muons resulting from a Z --> mu+ mu- decay
#--------------------------------------------------------------------------------

muTauPairZmumuHypothesesForAHtoMuTauCentralJetVeto = muTauPairZmumuHypotheses.clone(
    diCandidatePairSource = cms.InputTag('selectedMuTauPairsForAHtoMuTauPzetaDiffCumulative')
)

muTauPairVisMassHypothesesForAHtoMuTauCentralJetVeto = muTauPairVisMassHypotheses.clone(
    diCandidatePairSource = muTauPairZmumuHypothesesForAHtoMuTauCentralJetVeto.diCandidatePairSource
)
muTauPairVisMassHypothesesForAHtoMuTauCentralJetVeto.ZllHypotheses[0].src = \
  cms.InputTag('muTauPairZmumuHypothesesForAHtoMuTauCentralJetVeto')

muTauPairZmumuHypothesesForAHtoMuTauCentralJetVetoLooseMuonIsolation = muTauPairZmumuHypothesesForAHtoMuTauCentralJetVeto.clone(
    diCandidatePairSource = cms.InputTag('selectedMuTauPairsForAHtoMuTauPzetaDiffLooseMuonIsolationCumulative')
)

muTauPairVisMassHypothesesForAHtoMuTauCentralJetVetoLooseMuonIsolation = muTauPairVisMassHypothesesForAHtoMuTauCentralJetVeto.clone(
    diCandidatePairSource = muTauPairZmumuHypothesesForAHtoMuTauCentralJetVetoLooseMuonIsolation.diCandidatePairSource
)
muTauPairVisMassHypothesesForAHtoMuTauCentralJetVetoLooseMuonIsolation.ZllHypotheses[0].src = \
  cms.InputTag('muTauPairZmumuHypothesesForAHtoMuTauCentralJetVetoLooseMuonIsolation')

muTauPairZmumuHypothesesForAHtoMuTauCentralJetBtag = muTauPairZmumuHypotheses.clone(
    diCandidatePairSource = cms.InputTag('selectedMuTauPairsForAHtoMuTauValidCollinearApproxCumulative')
)

muTauPairVisMassHypothesesForAHtoMuTauCentralJetBtag = muTauPairVisMassHypotheses.clone(
    diCandidatePairSource = muTauPairZmumuHypothesesForAHtoMuTauCentralJetBtag.diCandidatePairSource
)
muTauPairVisMassHypothesesForAHtoMuTauCentralJetBtag.ZllHypotheses[0].src = \
  cms.InputTag('muTauPairZmumuHypothesesForAHtoMuTauCentralJetBtag')

muTauPairZmumuHypothesesForAHtoMuTauCentralJetBtagLooseMuonIsolation = muTauPairZmumuHypothesesForAHtoMuTauCentralJetBtag.clone(
    diCandidatePairSource = cms.InputTag('selectedMuTauPairsForAHtoMuTauValidCollinearApproxLooseMuonIsolationCumulative')
)

muTauPairVisMassHypothesesForAHtoMuTauCentralJetBtagLooseMuonIsolation = muTauPairVisMassHypothesesForAHtoMuTauCentralJetBtag.clone(
    diCandidatePairSource = muTauPairZmumuHypothesesForAHtoMuTauCentralJetBtagLooseMuonIsolation.diCandidatePairSource
)
muTauPairVisMassHypothesesForAHtoMuTauCentralJetBtagLooseMuonIsolation.ZllHypotheses[0].src = \
  cms.InputTag('muTauPairZmumuHypothesesForAHtoMuTauCentralJetBtagLooseMuonIsolation')

produceMuTauPairZmumuHypothesesForAHtoMuTau = cms.Sequence(
   muTauPairZmumuHypothesesForAHtoMuTauCentralJetVeto
  * muTauPairVisMassHypothesesForAHtoMuTauCentralJetVeto
  * muTauPairZmumuHypothesesForAHtoMuTauCentralJetVetoLooseMuonIsolation
  * muTauPairVisMassHypothesesForAHtoMuTauCentralJetVetoLooseMuonIsolation
  * muTauPairZmumuHypothesesForAHtoMuTauCentralJetBtag
  * muTauPairVisMassHypothesesForAHtoMuTauCentralJetBtag
  * muTauPairZmumuHypothesesForAHtoMuTauCentralJetBtagLooseMuonIsolation
  * muTauPairVisMassHypothesesForAHtoMuTauCentralJetBtagLooseMuonIsolation
)
