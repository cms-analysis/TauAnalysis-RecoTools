import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.RecoTools.pftauPatSelector_cfi import *

selectedLayer1TausEta21IndividualEff = copy.deepcopy(selectedLayer1TausEta21Individual)
selectedLayer1TausEta21IndividualEff.src = "allLayer1TausForTauAnalysesEff" 

selectedLayer1TausPt20IndividualEff = copy.deepcopy(selectedLayer1TausPt20Individual)
selectedLayer1TausPt20IndividualEff.src = selectedLayer1TausEta21IndividualEff.src

selectedLayer1TausLeadTrkIndividualEff = copy.deepcopy(selectedLayer1TausLeadTrkIndividual)
selectedLayer1TausLeadTrkIndividualEff.src = selectedLayer1TausEta21IndividualEff.src

selectedLayer1TausLeadTrkPtIndividualEff = copy.deepcopy(selectedLayer1TausLeadTrkPtIndividual)
selectedLayer1TausLeadTrkPtIndividualEff.src = selectedLayer1TausEta21IndividualEff.src

selectedLayer1TausTrkIsoIndividualEff = copy.deepcopy(selectedLayer1TausTrkIsoIndividual)
selectedLayer1TausTrkIsoIndividualEff.src = selectedLayer1TausEta21IndividualEff.src

selectedLayer1TausEcalIsoIndividualEff = copy.deepcopy(selectedLayer1TausEcalIsoIndividual)
selectedLayer1TausEcalIsoIndividualEff.src = selectedLayer1TausEta21IndividualEff.src

selectedLayer1TausProngIndividualEff = copy.deepcopy(selectedLayer1TausProngIndividual)
selectedLayer1TausProngIndividualEff.src = selectedLayer1TausEta21IndividualEff.src

selectedLayer1TausElectronVetoIndividualEff = copy.deepcopy(selectedLayer1TausElectronVetoIndividual)
selectedLayer1TausElectronVetoIndividualEff.src = selectedLayer1TausEta21IndividualEff.src

selectedLayer1TausMuonVetoIndividualEff = copy.deepcopy(selectedLayer1TausMuonVetoIndividual)
selectedLayer1TausMuonVetoIndividualEff.src = selectedLayer1TausEta21IndividualEff.src

selectPFTausForTauAnalysesEff = cms.Sequence( selectedLayer1TausEta21IndividualEff 
                                             *selectedLayer1TausPt20IndividualEff
                                             *selectedLayer1TausLeadTrkIndividualEff
                                             *selectedLayer1TausLeadTrkPtIndividualEff
                                             *selectedLayer1TausTrkIsoIndividualEff
                                             *selectedLayer1TausEcalIsoIndividualEff
                                             *selectedLayer1TausProngIndividualEff
                                             *selectedLayer1TausElectronVetoIndividualEff
                                             *selectedLayer1TausMuonVetoIndividualEff )

