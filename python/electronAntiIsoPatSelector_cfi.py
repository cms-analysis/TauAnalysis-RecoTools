import FWCore.ParameterSet.Config as cms
  
# module to select Electrons
# See https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePhysicsCutParser
# on how to use the cut-string

# selection of electron sample passing
# **inverted** isolation criteria
# of electronPatSelector_cfi.py

selectedLayer1ElectronsTrkAntiIso = cms.EDFilter("PATElectronSelector",
     src = cms.InputTag("selectedLayer1ElectronsHLTmatchCumulative"),
     cut = cms.string('trackIso > 0.'),
     filter = cms.bool(False)
)

selectedLayer1ElectronsEcalAntiIso = cms.EDFilter("PATElectronSelector",
     src = cms.InputTag("selectedLayer1ElectronsTrkAntiIso"),
     cut = cms.string('ecalIso > 3.8'),
     filter = cms.bool(False)
)

selectedLayer1ElectronsHcalAntiIso = cms.EDFilter("PATElectronSelector",
     src = cms.InputTag("selectedLayer1ElectronsEcalAntiIso"),
     cut = cms.string('hcalIso > 1.5'),
     filter = cms.bool(False)
)

selectedLayer1ElectronsTrkForAntiIso = cms.EDFilter("PATElectronSelector",
     src = cms.InputTag("selectedLayer1ElectronsHcalAntiIso"),
     cut = cms.string('gsfTrack != 0'),
     filter = cms.bool(False)
)
selectedLayer1ElectronsTrkIPforAntiIso = cms.EDProducer("PATElectronSelector",
     src = cms.InputTag("selectedLayer1ElectronsTrkForAntiIso"),
     cut = cms.string('abs(gsfTrack.d0) < 0.02'),
     filter = cms.bool(False)
)

selectAntiIsoElectrons = cms.Sequence( selectedLayer1ElectronsTrkAntiIso
                                      *selectedLayer1ElectronsEcalAntiIso
                                      *selectedLayer1ElectronsHcalAntiIso
                                      *selectedLayer1ElectronsTrkForAntiIso
                                      *selectedLayer1ElectronsTrkIPforAntiIso )
