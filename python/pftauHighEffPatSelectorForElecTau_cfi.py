import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.RecoTools.pftauPatSelectorForElecTau_cfi import *

selectedLayer1TausForElecTauAntiOverlapWithElectronsVetoEff = copy.deepcopy(selectedLayer1TausForElecTauAntiOverlapWithElectronsVeto)
selectedLayer1TausForElecTauAntiOverlapWithElectronsVetoEff.src = "allLayer1TausForTauAnalysesEff" 

selectedLayer1TausForElecTauEta21CumulativeEff = copy.deepcopy(selectedLayer1TausForElecTauEta21Cumulative)
selectedLayer1TausForElecTauEta21CumulativeEff.src = cms.InputTag("selectedLayer1TausForElecTauAntiOverlapWithElectronsVetoEff")

selectedLayer1TausForElecTauPt20CumulativeEff = copy.deepcopy(selectedLayer1TausForElecTauPt20Cumulative)
selectedLayer1TausForElecTauPt20CumulativeEff.src = cms.InputTag("selectedLayer1TausForElecTauEta21CumulativeEff")

selectedLayer1TausForElecTauLeadTrkCumulativeEff = copy.deepcopy(selectedLayer1TausForElecTauLeadTrkCumulative)
selectedLayer1TausForElecTauLeadTrkCumulativeEff.src = cms.InputTag("selectedLayer1TausForElecTauPt20CumulativeEff")

selectedLayer1TausForElecTauLeadTrkPtCumulativeEff = copy.deepcopy(selectedLayer1TausForElecTauLeadTrkPtCumulative)
selectedLayer1TausForElecTauLeadTrkPtCumulativeEff.src = cms.InputTag("selectedLayer1TausForElecTauLeadTrkCumulativeEff")

selectedLayer1TausForElecTauTrkIsoCumulativeEff = copy.deepcopy(selectedLayer1TausForElecTauTrkIsoCumulative)
selectedLayer1TausForElecTauTrkIsoCumulativeEff.src = cms.InputTag("selectedLayer1TausForElecTauLeadTrkPtCumulativeEff")

selectedLayer1TausForElecTauEcalIsoCumulativeEff = copy.deepcopy(selectedLayer1TausForElecTauEcalIsoCumulative)
selectedLayer1TausForElecTauEcalIsoCumulativeEff.src = cms.InputTag("selectedLayer1TausForElecTauTrkIsoCumulativeEff")

selectedLayer1TausForElecTauProngCumulativeEff = copy.deepcopy(selectedLayer1TausForElecTauProngCumulative)
selectedLayer1TausForElecTauProngCumulativeEff.src = cms.InputTag("selectedLayer1TausForElecTauEcalIsoCumulativeEff")

selectedLayer1TausForElecTauElectronVetoCumulativeEff = copy.deepcopy(selectedLayer1TausForElecTauElectronVetoCumulative)
selectedLayer1TausForElecTauElectronVetoCumulativeEff.src = cms.InputTag("selectedLayer1TausForElecTauProngCumulativeEff")

selectPFTausForElecTauEff = cms.Sequence( selectedLayer1TausForElecTauAntiOverlapWithElectronsVetoEff
                                         *selectedLayer1TausForElecTauEta21CumulativeEff 
                                         *selectedLayer1TausForElecTauPt20CumulativeEff
                                         *selectedLayer1TausForElecTauLeadTrkCumulativeEff
                                         *selectedLayer1TausForElecTauLeadTrkPtCumulativeEff
                                         *selectedLayer1TausForElecTauTrkIsoCumulativeEff
                                         *selectedLayer1TausForElecTauEcalIsoCumulativeEff
                                         *selectedLayer1TausForElecTauProngCumulativeEff
                                         *selectedLayer1TausForElecTauElectronVetoCumulativeEff )

