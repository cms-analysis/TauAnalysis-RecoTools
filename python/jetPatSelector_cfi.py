import FWCore.ParameterSet.Config as cms
import copy
  
# module to select Jets
# See https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePhysicsCutParser
# on how to use the cut-string

# select jets not identified as electron, muon or tau-jets
selectedLayer1JetsAntiOverlapWithLeptonsVeto = cms.EDFilter("PATJetAntiOverlapSelector",
     src = cms.InputTag("allLayer1Jets"),
     srcNotToBeFiltered = cms.VInputTag("selectedLayer1ElectronsTrkIPcumulative",
                                        "selectedLayer1MuonsTrkIPcumulative",
                                        "selectedLayer1TausForDiTauProngCumulative"),
     dRmin = cms.double(0.7),
     filter = cms.bool(False)                                           
)

# select central jets
selectedLayer1JetsEta21 = cms.EDFilter("PATJetSelector",
     src = cms.InputTag("selectedLayer1JetsAntiOverlapWithLeptonsVeto"),
     cut = cms.string('abs(eta) < 2.1'),
     filter = cms.bool(False)
)

# select jets with Et > 20 GeV
selectedLayer1JetsEt20 = cms.EDFilter("PATJetSelector",
     src = cms.InputTag("selectedLayer1JetsEta21"),
     cut = cms.string('et > 20.'),
     filter = cms.bool(False)
)

selectJetsForTauAnalyses = cms.Sequence( selectedLayer1JetsAntiOverlapWithLeptonsVeto
                                        *selectedLayer1JetsEta21
                                        *selectedLayer1JetsEt20 )
