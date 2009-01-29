import FWCore.ParameterSet.Config as cms
import copy

from RecoTauTag.RecoTau.PFRecoTauProducer_cfi import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByIsolation_cfi import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByLeadingTrackPtCut_cfi import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByLeadingTrackFinding_cfi import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationAgainstElectron_cfi import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationAgainstMuon_cfi import *

pfRecoTauProducerForTauAnalysesEff = copy.deepcopy(pfRecoTauProducer)
pfRecoTauProducerForTauAnalysesEff.PFTauTagInfoProducer = "pfRecoTauTagInfoProducerForTauAnalysesEff"
### Match ###    
pfRecoTauProducerForTauAnalysesEff.MatchingConeMetric = "DR"       
pfRecoTauProducerForTauAnalysesEff.MatchingConeSizeFormula = "0.1"
pfRecoTauProducerForTauAnalysesEff.MatchingConeSize_max = 0.1
pfRecoTauProducerForTauAnalysesEff.MatchingConeSize_min = 0.0
### Trk Sig ###
pfRecoTauProducerForTauAnalysesEff.TrackerSignalConeMetric = "DR"
pfRecoTauProducerForTauAnalysesEff.TrackerSignalConeSizeFormula = "0.15"
pfRecoTauProducerForTauAnalysesEff.TrackerSignalConeSize_min = 0.0
pfRecoTauProducerForTauAnalysesEff.TrackerSignalConeSize_max = 0.15
pfRecoTauProducerForTauAnalysesEff.LeadTrack_minPt = 5.0
pfRecoTauProducerForTauAnalysesEff.LeadChargedHadrCand_minPt = 5.0
### Trk Iso ###
pfRecoTauProducerForTauAnalysesEff.TrackerIsolConeMetric = "DR"
pfRecoTauProducerForTauAnalysesEff.TrackerIsolConeSizeFormula = "0.7" 
pfRecoTauProducerForTauAnalysesEff.TrackerIsolConeSize_max = "1.0"
pfRecoTauProducerForTauAnalysesEff.TrackerIsolConeSize_min = "0.0"
pfRecoTauProducerForTauAnalysesEff.Track_IsolAnnulus_minNhits = 8
pfRecoTauProducerForTauAnalysesEff.ChargedHadrCand_minPt = 0.3
### Ecl Sig ###
pfRecoTauProducerForTauAnalysesEff.ECALSignalConeMetric = "DR"   
pfRecoTauProducerForTauAnalysesEff.ECALSignalConeSizeFormula = "0.15"
pfRecoTauProducerForTauAnalysesEff.ECALSignalConeSize_max = "0.15"
pfRecoTauProducerForTauAnalysesEff.ECALSignalConeSize_min = "0.0"
### Ecl Iso ###
pfRecoTauProducerForTauAnalysesEff.ECALIsolConeMetric = "DR"     
pfRecoTauProducerForTauAnalysesEff.ECALIsolConeSizeFormula = "0.7"
pfRecoTauProducerForTauAnalysesEff.ECALIsolConeSize_max = "1.0"
pfRecoTauProducerForTauAnalysesEff.ECALIsolConeSize_min = "0.0"
pfRecoTauProducerForTauAnalysesEff.GammaCand_minPt = 1.5

### Discriminators ###
pfRecoTauLdgTrkFindForTauAnalysesEff = copy.deepcopy(pfRecoTauDiscriminationByLeadingTrackFinding)
pfRecoTauLdgTrkFindForTauAnalysesEff.PFTauProducer = 'allLayer0PFTausForTauAnalysesEff'

pfRecoTauLdgTrkPtCutForTauAnalysesEff = copy.deepcopy(pfRecoTauDiscriminationByLeadingTrackPtCut)
pfRecoTauLdgTrkPtCutForTauAnalysesEff.PFTauProducer = 'allLayer0PFTausForTauAnalysesEff'
pfRecoTauLdgTrkPtCutForTauAnalysesEff.MinPtLeadingTrack = 5.0

pfRecoTauTrkIsoDiscrForTauAnalysesEff = copy.deepcopy(pfRecoTauDiscriminationByIsolation)          
pfRecoTauTrkIsoDiscrForTauAnalysesEff.ApplyDiscriminationByECALIsolation = False
pfRecoTauTrkIsoDiscrForTauAnalysesEff.PFTauProducer = 'allLayer0PFTausForTauAnalysesEff'

pfRecoTauEcalIsoDiscrForTauAnalysesEff = copy.deepcopy(pfRecoTauDiscriminationByIsolation)          
pfRecoTauEcalIsoDiscrForTauAnalysesEff.ApplyDiscriminationByTrackerIsolation = False
pfRecoTauEcalIsoDiscrForTauAnalysesEff.ApplyDiscriminationByECALIsolation = True 
pfRecoTauEcalIsoDiscrForTauAnalysesEff.PFTauProducer = 'allLayer0PFTausForTauAnalysesEff'

pfRecoTauElecRejDiscrForTauAnalysesEff = copy.deepcopy(pfRecoTauDiscriminationAgainstElectron)
pfRecoTauElecRejDiscrForTauAnalysesEff.PFTauProducer = 'allLayer0PFTausForTauAnalysesEff'

pfRecoTauMuonRejDiscrForTauAnalysesEff = copy.deepcopy(pfRecoTauDiscriminationAgainstElectron)
pfRecoTauMuonRejDiscrForTauAnalysesEff.PFTauProducer = 'allLayer0PFTausForTauAnalysesEff' 

pfRecoTauIsoDiscrForTauAnalysesEff = copy.deepcopy(pfRecoTauDiscriminationByIsolation)
pfRecoTauIsoDiscrForTauAnalysesEff.PFTauProducer = 'pfRecoTauProducerForTauAnalysesEff'

PFTauForTauAnalysesEff = cms.Sequence(pfRecoTauProducerForTauAnalysesEff)
