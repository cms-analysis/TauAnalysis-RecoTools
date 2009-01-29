import FWCore.ParameterSet.Config as cms
import copy
  
# module to select Electrons
# See https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePhysicsCutParser
# on how to use the cut-string

# require electron candidate to pass the tight electron id. criteria
selectedLayer1ElectronsTightId = cms.EDFilter("PATElectronSelector",
     src = cms.InputTag("allLayer1ElectronsForTauAnalyses"),
     cut = cms.string('electronID("tight") > 0.'),
     filter = cms.bool(False)
)

# require electron candidate to not be within eta-crack
# between Barrel and Encap ECAL calorimeter
selectedLayer1ElectronsAntiCrackCutCumulative = cms.EDFilter("PATElectronSelector",
     src = cms.InputTag("selectedLayer1ElectronsTightId"),
     cut = cms.string('abs(superCluster.eta) < 1.442 || abs(superCluster.eta) > 1.560'),
     filter = cms.bool(False)
)

selectedLayer1ElectronsAntiCrackCutIndividual = copy.deepcopy(selectedLayer1ElectronsAntiCrackCutCumulative)
selectedLayer1ElectronsAntiCrackCutIndividual.src = selectedLayer1ElectronsTightId.src

# require electron candidate to be within geometric acceptance of electron trigger
selectedLayer1ElectronsEta21Cumulative = cms.EDFilter("PATElectronSelector",
     src = cms.InputTag("selectedLayer1ElectronsAntiCrackCutCumulative"),
     cut = cms.string('abs(eta) < 2.1'),
     filter = cms.bool(False)
)

selectedLayer1ElectronsEta21Individual = copy.deepcopy(selectedLayer1ElectronsEta21Cumulative)
selectedLayer1ElectronsEta21Individual.src = selectedLayer1ElectronsTightId.src

# require electron candidate to have transverse momentum above threshold
selectedLayer1ElectronsPt15Cumulative = cms.EDFilter("PATElectronSelector",
     src = cms.InputTag("selectedLayer1ElectronsEta21Cumulative"),
     cut = cms.string('pt > 15.'),
     filter = cms.bool(False)
)

selectedLayer1ElectronsPt15Individual = copy.deepcopy(selectedLayer1ElectronsPt15Cumulative)
selectedLayer1ElectronsPt15Individual.src = selectedLayer1ElectronsTightId.src

# require electron candidate to match L1 && HLT electron trigger 
selectedLayer1ElectronsHLTmatchCumulative = cms.EDFilter("PATElectronSelector",
     src = cms.InputTag("selectedLayer1ElectronsPt15Cumulative"),
     cut = cms.string('triggerMatches.size()'),
     filter = cms.bool(False)
)

selectedLayer1ElectronsHLTmatchIndividual = copy.deepcopy(selectedLayer1ElectronsHLTmatchCumulative)
selectedLayer1ElectronsHLTmatchIndividual.src = selectedLayer1ElectronsTightId.src

# require electron candidate to be isolated
# with respect to tracks (of Pt > 1. GeV)
selectedLayer1ElectronsTrkIsoCumulative = cms.EDFilter("PATElectronSelector",
     src = cms.InputTag("selectedLayer1ElectronsHLTmatchCumulative"),
     cut = cms.string('trackIso = 0.'),
     filter = cms.bool(False)
)

selectedLayer1ElectronsTrkIsoIndividual = copy.deepcopy(selectedLayer1ElectronsTrkIsoCumulative)
selectedLayer1ElectronsTrkIsoIndividual.src = selectedLayer1ElectronsTightId.src

# require electron candidate to be isolated
# with respect to energy deposits in ECAL
# (not associated to electron candidate)
selectedLayer1ElectronsEcalIsoCumulative = cms.EDFilter("PATElectronSelector",
     src = cms.InputTag("selectedLayer1ElectronsTrkIsoCumulative"),
     cut = cms.string('ecalIso < 3.8'),
     filter = cms.bool(False)
)

selectedLayer1ElectronsEcalIsoIndividual = copy.deepcopy(selectedLayer1ElectronsEcalIsoCumulative)
selectedLayer1ElectronsEcalIsoIndividual.src = selectedLayer1ElectronsTightId.src

# require electron candidate to be isolated
# with respect to energy deposits in HCAL
selectedLayer1ElectronsHcalIsoCumulative = cms.EDFilter("PATElectronSelector",
     src = cms.InputTag("selectedLayer1ElectronsEcalIsoCumulative"),
     cut = cms.string('hcalIso < 1.5'),
     filter = cms.bool(False)
)

selectedLayer1ElectronsHcalIsoIndividual = copy.deepcopy(selectedLayer1ElectronsHcalIsoCumulative)
selectedLayer1ElectronsHcalIsoIndividual.src = selectedLayer1ElectronsTightId.src

# require electron candidate to be linked to (GSF) track
selectedLayer1ElectronsTrkCumulative = cms.EDFilter("PATElectronSelector",
     src = cms.InputTag("selectedLayer1ElectronsHcalIsoCumulative"),
     cut = cms.string('gsfTrack.isNonnull'),
     filter = cms.bool(False)
)

selectedLayer1ElectronsTrkIndividual = copy.deepcopy(selectedLayer1ElectronsTrkCumulative)
selectedLayer1ElectronsTrkIndividual.src = selectedLayer1ElectronsTightId.src

# require track of electron candidate to have small transverse impact parameter
# (in order to veto electrons resulting from b-quark decays)
selectedLayer1ElectronsTrkIPcumulative = cms.EDProducer("PATElectronSelector",
     src = cms.InputTag("selectedLayer1ElectronsTrkCumulative"),
     cut = cms.string('abs(gsfTrack.d0) < 0.05'),
     filter = cms.bool(False)
)

selectedLayer1ElectronsTrkIPindividual = copy.deepcopy(selectedLayer1ElectronsTrkIPcumulative)
selectedLayer1ElectronsTrkIPindividual.src = selectedLayer1ElectronsTightId.src

selectElectronsForTauAnalyses = cms.Sequence( selectedLayer1ElectronsTightId
                                             *selectedLayer1ElectronsAntiCrackCutCumulative * selectedLayer1ElectronsAntiCrackCutIndividual
                                             *selectedLayer1ElectronsEta21Cumulative * selectedLayer1ElectronsEta21Individual
                                             *selectedLayer1ElectronsPt15Cumulative * selectedLayer1ElectronsPt15Individual
                                             *selectedLayer1ElectronsHLTmatchCumulative * selectedLayer1ElectronsHLTmatchIndividual
                                             *selectedLayer1ElectronsTrkIsoCumulative * selectedLayer1ElectronsTrkIsoIndividual
                                             *selectedLayer1ElectronsEcalIsoCumulative * selectedLayer1ElectronsEcalIsoIndividual
                                             *selectedLayer1ElectronsHcalIsoCumulative * selectedLayer1ElectronsHcalIsoIndividual
                                             *selectedLayer1ElectronsTrkCumulative * selectedLayer1ElectronsTrkIndividual
                                             *selectedLayer1ElectronsTrkIPcumulative * selectedLayer1ElectronsTrkIPindividual )
