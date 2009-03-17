import FWCore.ParameterSet.Config as cms
import copy

# module to select Muons
# See https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePhysicsCutParser
# on how to use the cut-string

# require muon candidate to be a global muon
# (track in muon system linked to track in Pixel + SiTracker detectors)
selectedLayer1MuonsGlobal = cms.EDFilter("PATMuonSelector",
     src = cms.InputTag("allLayer1MuonsForTauAnalyses"),
     cut = cms.string('isGlobalMuon()'),
     filter = cms.bool(False)
)

# require muon candidate to be within geometric acceptance of muon trigger
selectedLayer1MuonsEta21Cumulative = cms.EDFilter("PATMuonSelector",
     src = cms.InputTag("selectedLayer1MuonsGlobal"),
     cut = cms.string('abs(eta) < 2.1'),
     filter = cms.bool(False)
)

selectedLayer1MuonsEta21Individual = copy.deepcopy(selectedLayer1MuonsEta21Cumulative)
selectedLayer1MuonsEta21Individual.src = selectedLayer1MuonsGlobal.src

# require muon candidate to have transverse momentum above threshold
selectedLayer1MuonsPt15Cumulative = cms.EDFilter("PATMuonSelector",
     src = cms.InputTag("selectedLayer1MuonsEta21Cumulative"),
     cut = cms.string('pt > 15.'),
     filter = cms.bool(False)
)

selectedLayer1MuonsPt15Individual = copy.deepcopy(selectedLayer1MuonsPt15Cumulative)
selectedLayer1MuonsPt15Individual.src = selectedLayer1MuonsGlobal.src

# require muon candidate to match L1 && HLT muon trigger 
selectedLayer1MuonsHLTmatchCumulative = cms.EDFilter("PATMuonSelector",
     src = cms.InputTag("selectedLayer1MuonsPt15Cumulative"),
     cut = cms.string('triggerMatches.size()'),
     filter = cms.bool(False)
)

selectedLayer1MuonsHLTmatchIndividual = copy.deepcopy(selectedLayer1MuonsHLTmatchCumulative)
selectedLayer1MuonsHLTmatchIndividual.src = selectedLayer1MuonsGlobal.src

# require muon candidate to be isolated
# with respect to tracks (of Pt > 1. GeV)
#selectedLayer1MuonsTrkIsoCumulative = cms.EDFilter("PATMuonSelector",
#     src = cms.InputTag("selectedLayer1MuonsHLTmatchCumulative"),
#     cut = cms.string('trackIso = 0.'),
#     filter = cms.bool(False)
#)

selectedLayer1MuonsTrkIsoCumulative = cms.EDFilter("PATMuonIsoDepositSelector",
     src = cms.InputTag("selectedLayer1MuonsHLTmatchCumulative"),
     type = cms.string('tracker'),
     vetos = cms.vstring("0.01", "Threshold(0.9)"),
     #vetos = cms.vstring("0.01"),                          
     dRisoCone = cms.double(0.6),
     numMax = cms.int32(0),
     #sumPtMax = cms.double(0.9),
     filter = cms.bool(False)
)

selectedLayer1MuonsTrkIsoIndividual = copy.deepcopy(selectedLayer1MuonsTrkIsoCumulative)
selectedLayer1MuonsTrkIsoIndividual.src = selectedLayer1MuonsGlobal.src

# require muon candidate to be isolated
# with respect to energy deposits in ECAL
# (not associated to muon candidate)
selectedLayer1MuonsEcalIsoCumulative = cms.EDFilter("PATMuonSelector",
     src = cms.InputTag("selectedLayer1MuonsTrkIsoCumulative"),
     cut = cms.string('ecalIso < 1.'),
     filter = cms.bool(False)
)

selectedLayer1MuonsEcalIsoIndividual = copy.deepcopy(selectedLayer1MuonsEcalIsoCumulative)
selectedLayer1MuonsEcalIsoIndividual.src = selectedLayer1MuonsGlobal.src

# require muon candidate to pass pion veto
selectedLayer1MuonsPionVetoCumulative = cms.EDProducer("PATMuonAntiPionSelector",
     src = cms.InputTag("selectedLayer1MuonsEcalIsoCumulative"),
     CaloCompCoefficient = cms.double(0.8),
     SegmCompCoefficient = cms.double(1.2),
     AntiPionCut = cms.double(1.0),
     filter = cms.bool(False)
)

selectedLayer1MuonsPionVetoIndividual = copy.deepcopy(selectedLayer1MuonsPionVetoCumulative)
selectedLayer1MuonsPionVetoIndividual.src = selectedLayer1MuonsGlobal.src

# require muon candidate to be linked to track in silicon strip + pixel detectors
# (all global muons should be linked to tracks in the "inner" tracking detectors;
#  in case the muon is not linked to an "inner" track,
#  the track impact parameter selection will cause processing of the entire event to be skipped !!)
selectedLayer1MuonsTrkCumulative = cms.EDFilter("PATMuonSelector",
     src = cms.InputTag("selectedLayer1MuonsPionVetoCumulative"),
     cut = cms.string('innerTrack.isNonnull'),
     filter = cms.bool(False)
)

selectedLayer1MuonsTrkIndividual = copy.deepcopy(selectedLayer1MuonsTrkCumulative)
selectedLayer1MuonsTrkIndividual.src = selectedLayer1MuonsGlobal.src

# require track of muon candidate to have small transverse impact parameter
# (in order to veto muons resulting from b-quark decays)
selectedLayer1MuonsTrkIPcumulative = cms.EDProducer("PATMuonIpSelector",
     src = cms.InputTag("selectedLayer1MuonsTrkCumulative"),
     vertexSource = cms.InputTag("selectedPrimaryVertexPosition"),
     IpMax = cms.double(0.05),
     filter = cms.bool(False)                                               
)

selectedLayer1MuonsTrkIPindividual = copy.deepcopy(selectedLayer1MuonsTrkIPcumulative)
selectedLayer1MuonsTrkIPindividual.src = selectedLayer1MuonsGlobal.src

selectMuonsForTauAnalyses = cms.Sequence( selectedLayer1MuonsGlobal
                                         *selectedLayer1MuonsEta21Cumulative * selectedLayer1MuonsEta21Individual
                                         *selectedLayer1MuonsPt15Cumulative * selectedLayer1MuonsPt15Individual 
                                         *selectedLayer1MuonsHLTmatchCumulative * selectedLayer1MuonsHLTmatchIndividual
                                         *selectedLayer1MuonsTrkIsoCumulative * selectedLayer1MuonsTrkIsoIndividual
                                         *selectedLayer1MuonsEcalIsoCumulative * selectedLayer1MuonsEcalIsoIndividual
                                         *selectedLayer1MuonsPionVetoCumulative * selectedLayer1MuonsPionVetoIndividual
                                         *selectedLayer1MuonsTrkCumulative * selectedLayer1MuonsTrkIndividual
                                         *selectedLayer1MuonsTrkIPcumulative * selectedLayer1MuonsTrkIPindividual )

# define additional collections of muon candidates
# with loose track and ECAL isolation applied
# (NOTE: to be used for the purpose of factorizing efficiencies
#        of muon isolation from other event selection criteria,
#        in order to avoid problems with limited Monte Carlo statistics)

selectedLayer1MuonsTrkIsoCumulativeLooseMuonIsolation = copy.deepcopy(selectedLayer1MuonsTrkIsoCumulative)
selectedLayer1MuonsTrkIsoCumulativeLooseMuonIsolation.vetos = cms.vstring("0.01")
selectedLayer1MuonsTrkIsoCumulativeLooseMuonIsolation.numMax = cms.int32(-1)
selectedLayer1MuonsTrkIsoCumulativeLooseMuonIsolation.sumPtMax = cms.double(8.)

selectedLayer1MuonsTrkIsoIndividualLooseMuonIsolation = copy.deepcopy(selectedLayer1MuonsTrkIsoCumulativeLooseMuonIsolation)
selectedLayer1MuonsTrkIsoIndividualLooseMuonIsolation.src = selectedLayer1MuonsGlobal.src

selectedLayer1MuonsEcalIsoCumulativeLooseMuonIsolation = copy.deepcopy(selectedLayer1MuonsEcalIsoCumulative)
selectedLayer1MuonsEcalIsoCumulativeLooseMuonIsolation.src = cms.InputTag("selectedLayer1MuonsTrkIsoCumulativeLooseMuonIsolation")
selectedLayer1MuonsEcalIsoCumulativeLooseMuonIsolation.cut = cms.string('ecalIso < 8.')

selectedLayer1MuonsEcalIsoIndividualLooseMuonIsolation = copy.deepcopy(selectedLayer1MuonsEcalIsoCumulativeLooseMuonIsolation)
selectedLayer1MuonsEcalIsoIndividualLooseMuonIsolation.src = selectedLayer1MuonsGlobal.src

selectedLayer1MuonsPionVetoCumulativeLooseMuonIsolation = copy.deepcopy(selectedLayer1MuonsPionVetoCumulative)
selectedLayer1MuonsPionVetoCumulativeLooseMuonIsolation.src = cms.InputTag("selectedLayer1MuonsEcalIsoCumulativeLooseMuonIsolation")

selectedLayer1MuonsPionVetoIndividualLooseMuonIsolation = copy.deepcopy(selectedLayer1MuonsPionVetoCumulativeLooseMuonIsolation)
selectedLayer1MuonsPionVetoIndividualLooseMuonIsolation.src = selectedLayer1MuonsGlobal.src

selectedLayer1MuonsTrkCumulativeLooseMuonIsolation = copy.deepcopy(selectedLayer1MuonsTrkCumulative)
selectedLayer1MuonsTrkCumulativeLooseMuonIsolation.src = cms.InputTag("selectedLayer1MuonsPionVetoCumulativeLooseMuonIsolation")

selectedLayer1MuonsTrkIndividualLooseMuonIsolation = copy.deepcopy(selectedLayer1MuonsTrkCumulativeLooseMuonIsolation)
selectedLayer1MuonsTrkIndividualLooseMuonIsolation.src = selectedLayer1MuonsGlobal.src

selectedLayer1MuonsTrkIPcumulativeLooseMuonIsolation = copy.deepcopy(selectedLayer1MuonsTrkIPcumulative)
selectedLayer1MuonsTrkIPcumulativeLooseMuonIsolation.src = src = cms.InputTag("selectedLayer1MuonsTrkCumulativeLooseMuonIsolation")

selectedLayer1MuonsTrkIPindividualLooseMuonIsolation = copy.deepcopy(selectedLayer1MuonsTrkIPcumulativeLooseMuonIsolation)
selectedLayer1MuonsTrkIPindividualLooseMuonIsolation.src = selectedLayer1MuonsGlobal.src

selectMuonsForTauAnalysesLooseMuonIsolation = cms.Sequence( selectedLayer1MuonsTrkIsoCumulativeLooseMuonIsolation
                                                           *selectedLayer1MuonsTrkIsoIndividualLooseMuonIsolation
                                                           *selectedLayer1MuonsEcalIsoCumulativeLooseMuonIsolation
                                                           *selectedLayer1MuonsEcalIsoIndividualLooseMuonIsolation
                                                           *selectedLayer1MuonsPionVetoCumulativeLooseMuonIsolation
                                                           *selectedLayer1MuonsPionVetoIndividualLooseMuonIsolation
                                                           *selectedLayer1MuonsTrkCumulativeLooseMuonIsolation
                                                           *selectedLayer1MuonsTrkIndividualLooseMuonIsolation
                                                           *selectedLayer1MuonsTrkIPcumulativeLooseMuonIsolation
                                                           *selectedLayer1MuonsTrkIPindividualLooseMuonIsolation )
