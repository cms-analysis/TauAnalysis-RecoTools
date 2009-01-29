import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.RecoTools.pftauPatSelectorForMuTau_cfi import *

selectedLayer1TausForMuTauAntiOverlapWithMuonsVetoEff = copy.deepcopy(selectedLayer1TausForMuTauAntiOverlapWithMuonsVeto)
selectedLayer1TausForMuTauAntiOverlapWithMuonsVetoEff.src = "allLayer1TausForTauAnalysesEff" 

selectedLayer1TausForMuTauEta21CumulativeEff = copy.deepcopy(selectedLayer1TausForMuTauEta21Cumulative)
selectedLayer1TausForMuTauEta21CumulativeEff.src = cms.InputTag("selectedLayer1TausForMuTauAntiOverlapWithMuonsVetoEff")

selectedLayer1TausForMuTauPt20CumulativeEff = copy.deepcopy(selectedLayer1TausForMuTauPt20Cumulative)
selectedLayer1TausForMuTauPt20CumulativeEff.src = cms.InputTag("selectedLayer1TausForMuTauEta21CumulativeEff")

selectedLayer1TausForMuTauLeadTrkCumulativeEff = copy.deepcopy(selectedLayer1TausForMuTauLeadTrkCumulative)
selectedLayer1TausForMuTauLeadTrkCumulativeEff.src = cms.InputTag("selectedLayer1TausForMuTauPt20CumulativeEff")

selectedLayer1TausForMuTauLeadTrkPtCumulativeEff = copy.deepcopy(selectedLayer1TausForMuTauLeadTrkPtCumulative)
selectedLayer1TausForMuTauLeadTrkPtCumulativeEff.src = cms.InputTag("selectedLayer1TausForMuTauLeadTrkCumulativeEff")

selectedLayer1TausForMuTauTrkIsoCumulativeEff = copy.deepcopy(selectedLayer1TausForMuTauTrkIsoCumulative)
selectedLayer1TausForMuTauTrkIsoCumulativeEff.src = cms.InputTag("selectedLayer1TausForMuTauLeadTrkPtCumulativeEff")

selectedLayer1TausForMuTauEcalIsoCumulativeEff = copy.deepcopy(selectedLayer1TausForMuTauEcalIsoCumulative)
selectedLayer1TausForMuTauEcalIsoCumulativeEff.src = cms.InputTag("selectedLayer1TausForMuTauTrkIsoCumulativeEff")

selectedLayer1TausForMuTauProngCumulativeEff = copy.deepcopy(selectedLayer1TausForMuTauProngCumulative)
selectedLayer1TausForMuTauProngCumulativeEff.src = cms.InputTag("selectedLayer1TausForMuTauEcalIsoCumulativeEff")

selectedLayer1TausForMuTauMuonVetoCumulativeEff = copy.deepcopy(selectedLayer1TausForMuTauMuonVetoCumulative)
selectedLayer1TausForMuTauMuonVetoCumulativeEff.src = cms.InputTag("selectedLayer1TausForMuTauProngCumulativeEff")

selectPFTausForMuTauEff = cms.Sequence( selectedLayer1TausForMuTauAntiOverlapWithMuonsVetoEff
                                       *selectedLayer1TausForMuTauEta21CumulativeEff 
                                       *selectedLayer1TausForMuTauPt20CumulativeEff
                                       *selectedLayer1TausForMuTauLeadTrkCumulativeEff
                                       *selectedLayer1TausForMuTauLeadTrkPtCumulativeEff
                                       *selectedLayer1TausForMuTauTrkIsoCumulativeEff
                                       *selectedLayer1TausForMuTauEcalIsoCumulativeEff
                                       *selectedLayer1TausForMuTauProngCumulativeEff
                                       *selectedLayer1TausForMuTauMuonVetoCumulativeEff )

