import FWCore.ParameterSet.Config as cms
import copy

from RecoTauTag.RecoTau.PFRecoTauProducer_cfi import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByIsolation_cfi import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByLeadingTrackPtCut_cfi import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByLeadingTrackFinding_cfi import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByTrackIsolation_cfi import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByECALIsolation_cfi import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationAgainstElectron_cfi import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationAgainstMuon_cfi import *

pfRecoTauProducerForTauAnalyses = copy.deepcopy(pfRecoTauProducer)
### Match ###
pfRecoTauProducerForTauAnalyses.MatchingConeMetric = "DR"
pfRecoTauProducerForTauAnalyses.MatchingConeSizeFormula = "0.1"
pfRecoTauProducerForTauAnalyses.MatchingConeSize_max = 0.1
pfRecoTauProducerForTauAnalyses.MatchingConeSize_min = 0.0
### Trk Sig ###
pfRecoTauProducerForTauAnalyses.TrackerSignalConeMetric = "DR"
pfRecoTauProducerForTauAnalyses.TrackerSignalConeSizeFormula = "5.0/ET"
pfRecoTauProducerForTauAnalyses.TrackerSignalConeSize_min = 0.07
pfRecoTauProducerForTauAnalyses.TrackerSignalConeSize_max = 0.15
pfRecoTauProducerForTauAnalyses.LeadTrack_minPt = 5.0
pfRecoTauProducerForTauAnalyses.LeadChargedHadrCand_minPt = 5.0
### Trk Iso ###
pfRecoTauProducerForTauAnalyses.TrackerIsolConeMetric = "DR"
pfRecoTauProducerForTauAnalyses.TrackerIsolConeSizeFormula = "0.5"
pfRecoTauProducerForTauAnalyses.TrackerIsolConeSize_max = "1.0"
pfRecoTauProducerForTauAnalyses.TrackerIsolConeSize_min = "0.0"
pfRecoTauProducerForTauAnalyses.Track_IsolAnnulus_minNhits = 8
pfRecoTauProducerForTauAnalyses.ChargedHadrCand_IsolAnnulus_minNhits = 8
pfRecoTauProducerForTauAnalyses.ChargedHadrCand_minPt = 1.0
pfRecoTauProducerForTauAnalyses.Track_minPt = 1.0
### Ecl Sig ###
pfRecoTauProducerForTauAnalyses.ECALSignalConeMetric = "DR"
pfRecoTauProducerForTauAnalyses.ECALSignalConeSizeFormula = "0.15"
pfRecoTauProducerForTauAnalyses.ECALSignalConeSize_max = "0.15"
pfRecoTauProducerForTauAnalyses.ECALSignalConeSize_min = "0.0"
### Ecl Iso ###
pfRecoTauProducerForTauAnalyses.ECALIsolConeMetric = "DR"   
pfRecoTauProducerForTauAnalyses.ECALIsolConeSizeFormula = "0.50"
pfRecoTauProducerForTauAnalyses.ECALIsolConeSize_max = "1.0"
pfRecoTauProducerForTauAnalyses.ECALIsolConeSize_min = "0.0"
pfRecoTauProducerForTauAnalyses.GammaCand_minPt = 1.5
pfRecoTauProducerForTauAnalyses.DataType = cms.string("AOD")

### Discriminators ###
pfRecoTauLdgTrkFindForTauAnalyses = copy.deepcopy(pfRecoTauDiscriminationByLeadingTrackFinding)
pfRecoTauLdgTrkFindForTauAnalyses.PFTauProducer = 'allLayer0PFTausForTauAnalyses'

pfRecoTauLdgTrkPtCutForTauAnalyses = copy.deepcopy(pfRecoTauDiscriminationByLeadingTrackPtCut)
pfRecoTauLdgTrkPtCutForTauAnalyses.PFTauProducer = 'allLayer0PFTausForTauAnalyses'
pfRecoTauLdgTrkPtCutForTauAnalyses.MinPtLeadingTrack = 5.0

pfRecoTauTrkIsoDiscrForTauAnalyses = copy.deepcopy(pfRecoTauDiscriminationByTrackIsolation)
pfRecoTauTrkIsoDiscrForTauAnalyses.PFTauProducer = 'allLayer0PFTausForTauAnalyses'

pfRecoTauEcalIsoDiscrForTauAnalyses = copy.deepcopy(pfRecoTauDiscriminationByECALIsolation)          
pfRecoTauEcalIsoDiscrForTauAnalyses.PFTauProducer = 'allLayer0PFTausForTauAnalyses'

pfRecoTauElecRejDiscrForTauAnalyses = copy.deepcopy(pfRecoTauDiscriminationAgainstElectron)
pfRecoTauElecRejDiscrForTauAnalyses.PFTauProducer = 'allLayer0PFTausForTauAnalyses'

pfRecoTauMuonRejDiscrForTauAnalyses = copy.deepcopy(pfRecoTauDiscriminationAgainstMuon)
pfRecoTauMuonRejDiscrForTauAnalyses.PFTauProducer = 'allLayer0PFTausForTauAnalyses' 

pfRecoTauIsoDiscrForTauAnalyses = copy.deepcopy(pfRecoTauDiscriminationByIsolation)
pfRecoTauIsoDiscrForTauAnalyses.PFTauProducer = 'pfRecoTauProducerForTauAnalyses'


