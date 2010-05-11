import FWCore.ParameterSet.Config as cms

from TauAnalysis.RecoTools.patElectronSelection_cfi import *
from TauAnalysis.RecoTools.patElectronSelectionForElecTau_cfi import *
from TauAnalysis.RecoTools.patElectronSelectionForElecMu_cfi import *
from TauAnalysis.RecoTools.patMuonSelection_cfi import *
from TauAnalysis.RecoTools.patPFTauSelection_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForElecTau_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForMuTau_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForDiTau_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForWTauNu_cfi import *

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *

#--------------------------------------------------------------------------------
# define selection criteria for pat::Electrons
# (settings made here overwrite values defined in electronPatSelector_cfi)
#--------------------------------------------------------------------------------

selectedLayer1ElectronsTightId.cut = cms.string('(abs(superCluster.eta) < 1.479 & electronID("eidRobustTight") > 0 & eSuperClusterOverP < 1.05 & eSuperClusterOverP > 0.95) | (abs(superCluster.eta) > 1.479 & electronID("eidRobustTight") > 0 & eSuperClusterOverP < 1.12 & eSuperClusterOverP > 0.95)')
selectedLayer1ElectronsLooseId.cut = cms.string('(abs(superCluster.eta) < 1.479 & electronID("eidRobustLoose") > 0 & eSuperClusterOverP < 1.05 & eSuperClusterOverP > 0.95) | (abs(superCluster.eta) > 1.479 & electronID("eidRobustLoose") > 0 & eSuperClusterOverP < 1.12 & eSuperClusterOverP > 0.95)')
selectedLayer1ElectronsAntiCrackCut.cut = cms.string('abs(superCluster.eta) < 1.442 | abs(superCluster.eta) > 1.560')
selectedLayer1ElectronsEta21.cut = cms.string('abs(eta) < 2.1')
selectedLayer1ElectronsPt15.cut = cms.string('pt > 15.')
selectedLayer1ElectronsTrkIso.cut = cms.string('userIsolation("pat::TrackIso") < 1.')
selectedLayer1ElectronsEcalIso.cut = cms.string('(abs(superCluster.eta) < 1.479 & userIsolation("pat::EcalIso") < 2.5) | (abs(superCluster.eta) > 1.479 & userIsolation("pat::EcalIso") < 3.5)')
selectedLayer1ElectronsTrk.cut = cms.string('gsfTrack.isNonnull')
selectedLayer1ElectronsTrkIP.vertexSource = cms.InputTag("selectedPrimaryVertexPosition")
selectedLayer1ElectronsTrkIP.IpMax = cms.double(0.05)

patElectronSelConfigurator = objSelConfigurator(
    [ selectedLayer1ElectronsTightId,
      selectedLayer1ElectronsAntiCrackCut,
      selectedLayer1ElectronsEta21,
      selectedLayer1ElectronsPt15,
      selectedLayer1ElectronsTrkIso,
      selectedLayer1ElectronsEcalIso,
      selectedLayer1ElectronsTrk,
      selectedLayer1ElectronsTrkIP ],
    src = "cleanPatElectrons",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectLayer1Electrons = patElectronSelConfigurator.configure(pyNameSpace = locals())

selectedLayer1ElectronsTrkIsoLooseIsolation.cut = cms.string('userIsolation("pat::TrackIso") < 8.')
selectedLayer1ElectronsEcalIsoLooseIsolation.cut = cms.string('userIsolation("pat::EcalIso") < 8.')
selectedLayer1ElectronsTrkLooseIsolation.cut = selectedLayer1ElectronsTrk.cut
selectedLayer1ElectronsTrkIPlooseIsolation.vertexSource = selectedLayer1ElectronsTrkIP.vertexSource
selectedLayer1ElectronsTrkIPlooseIsolation.IpMax = selectedLayer1ElectronsTrkIP.IpMax

patElectronSelConfiguratorLooseIsolation = objSelConfigurator(
    [ selectedLayer1ElectronsTightId,
      selectedLayer1ElectronsAntiCrackCut,
      selectedLayer1ElectronsEta21,
      selectedLayer1ElectronsPt15,
      selectedLayer1ElectronsTrkIsoLooseIsolation,
      selectedLayer1ElectronsEcalIsoLooseIsolation,
      selectedLayer1ElectronsTrkLooseIsolation,
      selectedLayer1ElectronsTrkIPlooseIsolation ],
    src = "cleanPatElectrons",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectLayer1ElectronsLooseIsolation = patElectronSelConfiguratorLooseIsolation.configure(pyNameSpace = locals())

#
# select electrons for Z->electron + tau-jet analysis
#

selectedLayer1ElectronsForElecTauTightId.cut = cms.string('(abs(superCluster.eta) < 1.479 & electronID("eidRobustTight") > 0 & eSuperClusterOverP < 1.05 & eSuperClusterOverP > 0.95) | (abs(superCluster.eta) > 1.479 & electronID("eidRobustTight") > 0 & eSuperClusterOverP < 1.12 & eSuperClusterOverP > 0.95)')
selectedLayer1ElectronsForElecTauLooseId.cut = cms.string('(abs(superCluster.eta) < 1.479 & electronID("eidRobustLoose") > 0 & eSuperClusterOverP < 1.4 & eSuperClusterOverP > 0.8) | (abs(superCluster.eta) > 1.479 & electronID("eidRobustLoose") > 0 & eSuperClusterOverP < 1.6 & eSuperClusterOverP > 0.8)')
selectedLayer1ElectronsForElecTauAntiCrackCut.cut = cms.string('abs(superCluster.eta) < 1.442 | abs(superCluster.eta) > 1.560')
selectedLayer1ElectronsForElecTauEta21.cut = cms.string('abs(eta) < 2.1')
selectedLayer1ElectronsForElecTauPt15.cut = cms.string('pt > 15.')
selectedLayer1ElectronsForElecTauTrkIso.cut = cms.string('userIsolation("pat::TrackIso") < 1.')
selectedLayer1ElectronsForElecTauEcalIso.cut = cms.string('(abs(superCluster.eta) < 1.479 & userIsolation("pat::EcalIso") < 2.5) | (abs(superCluster.eta) > 1.479 & userIsolation("pat::EcalIso") < 3.5)')
selectedLayer1ElectronsForElecTauTrk.cut = cms.string('gsfTrack.isNonnull')
selectedLayer1ElectronsForElecTauTrkIP.vertexSource = cms.InputTag("selectedPrimaryVertexPosition")
selectedLayer1ElectronsForElecTauTrkIP.IpMax = cms.double(0.05)
selectedLayer1ElectronsForElecTauConversionVeto.cotThetaCut = cms.double(0.05)
selectedLayer1ElectronsForElecTauConversionVeto.docaElecTrack = cms.double(0)
selectedLayer1ElectronsForElecTauConversionVeto.dRElecTrack = cms.double(0.1)
selectedLayer1ElectronsForElecTauConversionVeto.doPixCut = cms.bool(True)
selectedLayer1ElectronsForElecTauConversionVeto.nTrkMax = cms.double(1)
selectedLayer1ElectronsForElecTauConversionVeto.useConversionColl = cms.bool(True)

patElectronSelConfiguratorForElecTau = objSelConfigurator(
    [ selectedLayer1ElectronsForElecTauLooseId,
      selectedLayer1ElectronsForElecTauAntiCrackCut,
      selectedLayer1ElectronsForElecTauEta21,
      selectedLayer1ElectronsForElecTauPt15,
      selectedLayer1ElectronsForElecTauTrkIso,
      selectedLayer1ElectronsForElecTauEcalIso,
      selectedLayer1ElectronsForElecTauTrk,
      selectedLayer1ElectronsForElecTauTrkIP,
	  selectedLayer1ElectronsForElecTauConversionVeto ],
    src = "cleanPatElectrons",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectLayer1ElectronsForElecTau = patElectronSelConfiguratorForElecTau.configure(pyNameSpace = locals())

patElectronSelConfiguratorForElecTauLooseIsolation = objSelConfigurator(
    [ selectedLayer1ElectronsForElecTauLooseId,
      selectedLayer1ElectronsForElecTauAntiCrackCut,
      selectedLayer1ElectronsForElecTauEta21,
      selectedLayer1ElectronsForElecTauPt15,
      selectedLayer1ElectronsForElecTauTrkIsoLooseIsolation,
      selectedLayer1ElectronsForElecTauEcalIsoLooseIsolation,
      selectedLayer1ElectronsForElecTauTrkLooseIsolation,
      selectedLayer1ElectronsForElecTauTrkIPlooseIsolation,
	  selectedLayer1ElectronsForElecTauConversionVetoLooseIsolation ],
    src = "cleanPatElectrons",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectLayer1ElectronsForElecTauLooseIsolation = patElectronSelConfiguratorForElecTauLooseIsolation.configure(pyNameSpace = locals())

#
# select electrons for Z->electron + muon analysis
#

selectedLayer1ElectronsForElecMuAntiOverlapWithMuonsVeto.dRmin = cms.double(0.3)
selectedLayer1ElectronsForElecMuTightId.cut = cms.string('electronID("eidRobustTight") > 0')
selectedLayer1ElectronsForElecMuAntiCrackCut.cut = selectedLayer1ElectronsAntiCrackCut.cut
selectedLayer1ElectronsForElecMuEta21.cut = selectedLayer1ElectronsEta21.cut 
selectedLayer1ElectronsForElecMuPt15.cut = selectedLayer1ElectronsPt15.cut 
selectedLayer1ElectronsForElecMuTrkIso.cut = cms.string('userIsolation("pat::TrackIso") < 2.')
selectedLayer1ElectronsForElecMuEcalIso.cut = selectedLayer1ElectronsEcalIso.cut
selectedLayer1ElectronsForElecMuTrk.cut = selectedLayer1ElectronsTrk.cut 
selectedLayer1ElectronsForElecMuTrkIP.vertexSource = selectedLayer1ElectronsTrkIP.vertexSource 
selectedLayer1ElectronsForElecMuTrkIP.IpMax = selectedLayer1ElectronsTrkIP.IpMax 

patElectronSelConfiguratorForElecMu = objSelConfigurator(
    [ selectedLayer1ElectronsForElecMuAntiOverlapWithMuonsVeto,
      selectedLayer1ElectronsForElecMuTightId,
      selectedLayer1ElectronsForElecMuAntiCrackCut,
      selectedLayer1ElectronsForElecMuEta21,
      selectedLayer1ElectronsForElecMuPt15,
      selectedLayer1ElectronsForElecMuTrkIso,
      selectedLayer1ElectronsForElecMuEcalIso,
      selectedLayer1ElectronsForElecMuTrk,
      selectedLayer1ElectronsForElecMuTrkIP ],
    src = "cleanPatElectrons",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectLayer1ElectronsForElecMu = patElectronSelConfiguratorForElecMu.configure(pyNameSpace = locals())

selectedLayer1ElectronsForElecMuTrkIsoLooseIsolation.cut = cms.string('userIsolation("pat::TrackIso") < 8.')
selectedLayer1ElectronsForElecMuEcalIsoLooseIsolation.cut = cms.string('userIsolation("pat::EcalIso") < 8.')
selectedLayer1ElectronsForElecMuTrkLooseIsolation.cut = selectedLayer1ElectronsForElecMuTrk.cut
selectedLayer1ElectronsForElecMuTrkIPlooseIsolation.vertexSource = selectedLayer1ElectronsForElecMuTrkIP.vertexSource
selectedLayer1ElectronsForElecMuTrkIPlooseIsolation.IpMax = selectedLayer1ElectronsForElecMuTrkIP.IpMax

patElectronSelConfiguratorForElecMuLooseIsolation = objSelConfigurator(
    [ selectedLayer1ElectronsForElecMuAntiOverlapWithMuonsVeto,
      selectedLayer1ElectronsForElecMuTightId,
      selectedLayer1ElectronsForElecMuAntiCrackCut,
      selectedLayer1ElectronsForElecMuEta21,
      selectedLayer1ElectronsForElecMuPt15,
      selectedLayer1ElectronsForElecMuTrkIsoLooseIsolation,
      selectedLayer1ElectronsForElecMuEcalIsoLooseIsolation,
      selectedLayer1ElectronsForElecMuTrkLooseIsolation,
      selectedLayer1ElectronsForElecMuTrkIPlooseIsolation ],
    src = "cleanPatElectrons",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectLayer1ElectronsForElecMuLooseIsolation = patElectronSelConfiguratorForElecMuLooseIsolation.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------
# define selection criteria for pat::Muons
# (settings made here overwrite values defined in muonPatSelector_cfi)
#--------------------------------------------------------------------------------

selectedLayer1MuonsGlobal.cut = cms.string('isGlobalMuon()')
selectedLayer1MuonsEta21.cut = cms.string('abs(eta) < 2.1')
selectedLayer1MuonsPt15.cut = cms.string('pt > 15.')
selectedLayer1MuonsTrkIso.vetos = vetos = cms.vstring("0.01")
selectedLayer1MuonsTrkIso.dRisoCone = cms.double(0.6)
selectedLayer1MuonsTrkIso.sumPtMax = cms.double(1.)
selectedLayer1MuonsTrkIso.sumPtMethod = cms.string("absolute")
selectedLayer1MuonsEcalIso.cut = cms.string('userIsolation("pat::EcalIso") < 1.')
selectedLayer1MuonsPionVeto.CaloCompCoefficient = cms.double(0.8)
selectedLayer1MuonsPionVeto.SegmCompCoefficient = cms.double(1.2)
selectedLayer1MuonsPionVeto.AntiPionCut = cms.double(1.0)
selectedLayer1MuonsTrk.cut = cms.string('innerTrack.isNonnull')
selectedLayer1MuonsTrkIP.vertexSource = cms.InputTag("selectedPrimaryVertexPosition")
selectedLayer1MuonsTrkIP.IpMax = cms.double(0.05)

patMuonSelConfigurator = objSelConfigurator(
    [ selectedLayer1MuonsGlobal,
      selectedLayer1MuonsEta21,
      selectedLayer1MuonsPt15,
      selectedLayer1MuonsTrkIso,
      selectedLayer1MuonsEcalIso,
      selectedLayer1MuonsPionVeto,
      selectedLayer1MuonsTrk,
      selectedLayer1MuonsTrkIP ],
    src = "cleanPatMuons",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectLayer1Muons = patMuonSelConfigurator.configure(pyNameSpace = locals())

selectedLayer1MuonsTrkIsoLooseIsolation.vetos = cms.vstring("0.01")
selectedLayer1MuonsTrkIsoLooseIsolation.numMax = cms.int32(-1)
selectedLayer1MuonsTrkIsoLooseIsolation.sumPtMax = cms.double(8.)
selectedLayer1MuonsTrkIsoLooseIsolation.sumPtMethod = cms.string("absolute")
selectedLayer1MuonsEcalIsoLooseIsolation.cut = cms.string('userIsolation("pat::EcalIso") < 8.')
selectedLayer1MuonsPionVetoLooseIsolation.CaloCompCoefficient = selectedLayer1MuonsPionVeto.CaloCompCoefficient
selectedLayer1MuonsPionVetoLooseIsolation.SegmCompCoefficient = selectedLayer1MuonsPionVeto.SegmCompCoefficient
selectedLayer1MuonsPionVetoLooseIsolation.AntiPionCut = selectedLayer1MuonsPionVeto.AntiPionCut
selectedLayer1MuonsTrkLooseIsolation.cut = selectedLayer1MuonsTrk.cut
selectedLayer1MuonsTrkIPlooseIsolation.vertexSource = selectedLayer1MuonsTrkIP.vertexSource
selectedLayer1MuonsTrkIPlooseIsolation.IpMax = selectedLayer1MuonsTrkIP.IpMax

patMuonSelConfiguratorLooseIsolation = objSelConfigurator(
    [ selectedLayer1MuonsGlobal,
      selectedLayer1MuonsEta21,
      selectedLayer1MuonsPt15,
      selectedLayer1MuonsTrkIsoLooseIsolation,
      selectedLayer1MuonsEcalIsoLooseIsolation,
      selectedLayer1MuonsPionVetoLooseIsolation,
      selectedLayer1MuonsTrkLooseIsolation,
      selectedLayer1MuonsTrkIPlooseIsolation ],
    src = "cleanPatMuons",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectLayer1MuonsLooseIsolation = patMuonSelConfiguratorLooseIsolation.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------
# define selection criteria for pat::(PF)Taus
# (settings made here overwrite values defined in pftauPatSelector_cfi)
#--------------------------------------------------------------------------------

selectedLayer1TausEta21.cut = cms.string("abs(eta) < 2.1")
selectedLayer1TausPt20.cut = cut = cms.string("pt > 20.")
selectedLayer1TausLeadTrk.cut = cms.string('tauID("leadingTrackFinding") > 0.5')
selectedLayer1TausLeadTrkPt.cut = cms.string('tauID("leadingTrackPtCut") > 0.5')
selectedLayer1TausTaNCdiscr.cut = cms.string('tauID("byTaNCfrQuarterPercent") > -1.e+3') # cut on TaNC output disabled per default
selectedLayer1TausTrkIso.cut = cms.string('tauID("trackIsolation") > 0.5')
selectedLayer1TausEcalIso.cut = cms.string('tauID("ecalIsolation") > 0.5')
selectedLayer1TausProng.cut = cms.string("signalPFChargedHadrCands.size() = 1 | signalPFChargedHadrCands.size() = 3")
selectedLayer1TausCharge.cut = cms.string('abs(charge) > 0.5 & abs(charge) < 1.5')
selectedLayer1TausMuonVeto.cut = cms.string('tauID("againstMuon") > 0.5')
selectedLayer1TausElectronVeto.cut = cms.string('tauID("againstElectron") > 0.5')
selectedLayer1TausEcalCrackVeto.cut = cms.string("abs(eta) < 1.460 | abs(eta) > 1.558")

patTauSelConfigurator = objSelConfigurator(
    [ selectedLayer1TausEta21,
      selectedLayer1TausPt20,
      selectedLayer1TausLeadTrk,
      selectedLayer1TausLeadTrkPt,
      selectedLayer1TausTaNCdiscr,
      selectedLayer1TausTrkIso,
      selectedLayer1TausEcalIso,
      selectedLayer1TausProng,
      selectedLayer1TausCharge,
      selectedLayer1TausMuonVeto,
      selectedLayer1TausElectronVeto,
      selectedLayer1TausEcalCrackVeto ],
    src = "cleanPatTaus",
    pyModuleName = __name__,
    doSelCumulative = True,
    doSelIndividual = True
)

selectLayer1Taus = patTauSelConfigurator.configure(pyNameSpace = locals())
#
# define collections of pat::(PF)Taus used in semi-leptonic e + tau-jet channel
# (require electron and tau-jet candidates to be separated in eta-phi,
#  in order to avoid double-counting the same particle both as electron and as tau-jet;
#  apply anti-electron veto only; no need to apply anti-muon veto)
#
selectedLayer1TausForElecTauAntiOverlapWithElectronsVeto.dRmin = cms.double(0.3)
selectedLayer1TausForElecTauEta21.cut = selectedLayer1TausEta21.cut
selectedLayer1TausForElecTauPt20.cut = selectedLayer1TausPt20.cut
selectedLayer1TausForElecTauLeadTrk.cut = selectedLayer1TausLeadTrk.cut
selectedLayer1TausForElecTauLeadTrkPt.cut = selectedLayer1TausLeadTrkPt.cut
selectedLayer1TausForElecTauTaNCdiscr.cut = cms.string('tauID("byTaNCfrQuarterPercent") > -1.e+3') # cut on TaNC output disabled per default
selectedLayer1TausForElecTauTrkIso.cut = selectedLayer1TausTrkIso.cut
selectedLayer1TausForElecTauEcalIso.cut = selectedLayer1TausEcalIso.cut
selectedLayer1TausForElecTauProng.cut = selectedLayer1TausProng.cut
selectedLayer1TausForElecTauCharge.cut = selectedLayer1TausCharge.cut
selectedLayer1TausForElecTauElectronVeto.cut = selectedLayer1TausElectronVeto.cut
selectedLayer1TausForElecTauEcalCrackVeto.cut =  selectedLayer1TausEcalCrackVeto.cut
selectedLayer1TausForElecTauMuonVeto.cut = selectedLayer1TausMuonVeto.cut

patTauSelConfiguratorForElecTau = objSelConfigurator(
    [ selectedLayer1TausForElecTauAntiOverlapWithElectronsVeto,
      selectedLayer1TausForElecTauEta21,
      selectedLayer1TausForElecTauPt20,
      selectedLayer1TausForElecTauLeadTrk,
      selectedLayer1TausForElecTauLeadTrkPt,
      selectedLayer1TausForElecTauTaNCdiscr,
      selectedLayer1TausForElecTauTrkIso,
      selectedLayer1TausForElecTauEcalIso,
      selectedLayer1TausForElecTauProng,
      selectedLayer1TausForElecTauCharge,
      selectedLayer1TausForElecTauElectronVeto,
      selectedLayer1TausForElecTauEcalCrackVeto,
      selectedLayer1TausForElecTauMuonVeto ],
    src = "cleanPatTaus",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectLayer1TausForElecTau = patTauSelConfiguratorForElecTau.configure(pyNameSpace = locals())
#
# define collections of pat::(PF)Taus used in semi-leptonic mu + tau-jet channel
# (require muon and tau-jet candidates to be separated in eta-phi,
#  in order to avoid double-counting the same particle both as muon and as tau-jet;
#  apply anti-muon veto only; no need to apply anti-electron veto)
#
selectedLayer1TausForMuTauAntiOverlapWithMuonsVeto.dRmin = cms.double(0.3)
selectedLayer1TausForMuTauEta21.cut = selectedLayer1TausEta21.cut
selectedLayer1TausForMuTauPt20.cut = selectedLayer1TausPt20.cut
selectedLayer1TausForMuTauLeadTrk.cut = selectedLayer1TausLeadTrk.cut
selectedLayer1TausForMuTauLeadTrkPt.cut = selectedLayer1TausLeadTrkPt.cut
selectedLayer1TausForMuTauTaNCdiscr.cut = cms.string('tauID("byTaNCfrQuarterPercent") > -1.e+3') # cut on TaNC output disabled per default
selectedLayer1TausForMuTauTrkIso.cut = selectedLayer1TausTrkIso.cut
selectedLayer1TausForMuTauEcalIso.cut = selectedLayer1TausEcalIso.cut
selectedLayer1TausForMuTauProng.cut = selectedLayer1TausProng.cut
selectedLayer1TausForMuTauCharge.cut = selectedLayer1TausCharge.cut
selectedLayer1TausForMuTauMuonVeto.cut = selectedLayer1TausMuonVeto.cut
selectedLayer1TausForMuTauElectronVeto.cut = selectedLayer1TausElectronVeto.cut

patTauSelConfiguratorForMuTau = objSelConfigurator(
    [ selectedLayer1TausForMuTauAntiOverlapWithMuonsVeto,
      selectedLayer1TausForMuTauEta21,
      selectedLayer1TausForMuTauPt20,
      selectedLayer1TausForMuTauLeadTrk,
      selectedLayer1TausForMuTauLeadTrkPt,
      selectedLayer1TausForMuTauTaNCdiscr,
      selectedLayer1TausForMuTauTrkIso,
      selectedLayer1TausForMuTauEcalIso,
      selectedLayer1TausForMuTauProng,
      selectedLayer1TausForMuTauCharge,
      selectedLayer1TausForMuTauMuonVeto,
      selectedLayer1TausForMuTauElectronVeto ],
    src = "cleanPatTaus",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectLayer1TausForMuTau = patTauSelConfiguratorForMuTau.configure(pyNameSpace = locals())
#
# define collections of pat::(PF)Taus used in pure hadronic tau-jet + tau-jet channel
# (no need to apply anti-electron or anti-muon vetos)
#
selectedLayer1TausForDiTauEta21.cut = selectedLayer1TausEta21.cut
selectedLayer1TausForDiTauPt20.cut = selectedLayer1TausPt20.cut
selectedLayer1TausForDiTauLeadTrk.cut = selectedLayer1TausLeadTrk.cut
selectedLayer1TausForDiTauLeadTrkPt.cut = selectedLayer1TausLeadTrkPt.cut
selectedLayer1TausForDiTauTaNCdiscr.cut = cms.string('tauID("byTaNCfrQuarterPercent") > -1.e+3') # cut on TaNC output disabled per default
selectedLayer1TausForDiTauTrkIso.cut = selectedLayer1TausTrkIso.cut
selectedLayer1TausForDiTauEcalIso.cut = selectedLayer1TausEcalIso.cut
selectedLayer1TausForDiTauProng.cut = selectedLayer1TausProng.cut
selectedLayer1TausForDiTauCharge.cut = selectedLayer1TausCharge.cut

patTauSelConfiguratorForDiTau = objSelConfigurator(
    [ selectedLayer1TausForDiTauEta21,
      selectedLayer1TausForDiTauPt20,
      selectedLayer1TausForDiTauLeadTrk,
      selectedLayer1TausForDiTauLeadTrkPt,
      selectedLayer1TausForDiTauTaNCdiscr,
      selectedLayer1TausForDiTauTrkIso,
      selectedLayer1TausForDiTauEcalIso,
      selectedLayer1TausForDiTauProng,
      selectedLayer1TausForDiTauCharge ],
    src = "cleanPatTaus",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectLayer1TausForDiTau = patTauSelConfiguratorForDiTau.configure(pyNameSpace = locals())

# define collections of pat::(PF)Taus used in W->tau-jet + nu channel

selectedLayer1TausForWTauNuEta21.cut = selectedLayer1TausEta21.cut
selectedLayer1TausForWTauNuPt20.cut = selectedLayer1TausPt20.cut 
selectedLayer1TausForWTauNuLeadTrk.cut = selectedLayer1TausLeadTrk.cut
selectedLayer1TausForWTauNuLeadTrkPt.cut = cms.string("leadPFChargedHadrCand().isNonnull() & leadPFChargedHadrCand().pt() > 15.")
selectedLayer1TausForWTauNuTaNCdiscr.cut = cms.string('tauID("byTaNCfrQuarterPercent") > -1.e+3') # cut on TaNC output disabled per default
selectedLayer1TausForWTauNuEcalIso.cut = cms.string('tauID("byIsolation") > 0.5')
selectedLayer1TausForWTauNuTrkIso.cut = selectedLayer1TausTrkIso.cut
selectedLayer1TausForWTauNuProng.cut = selectedLayer1TausProng.cut
selectedLayer1TausForWTauNuCharge.cut = selectedLayer1TausCharge.cut
selectedLayer1TausForWTauNuMuonVeto.cut = selectedLayer1TausMuonVeto.cut
selectedLayer1TausForWTauNuElectronVeto.cut = selectedLayer1TausElectronVeto.cut
selectedLayer1TausForWTauNuEcalCrackVeto.cut = selectedLayer1TausEcalCrackVeto.cut

patTauSelConfiguratorForWTauNu = objSelConfigurator(
    [ selectedLayer1TausForWTauNuEta21,
      selectedLayer1TausForWTauNuPt20,
      selectedLayer1TausForWTauNuLeadTrk,
      selectedLayer1TausForWTauNuLeadTrkPt,
      selectedLayer1TausForWTauNuTaNCdiscr,
      selectedLayer1TausForWTauNuEcalIso,
      selectedLayer1TausForWTauNuTrkIso,
      selectedLayer1TausForWTauNuProng,
      selectedLayer1TausForWTauNuCharge, 
      selectedLayer1TausForWTauNuMuonVeto,
      selectedLayer1TausForWTauNuElectronVeto,
      selectedLayer1TausForWTauNuEcalCrackVeto ],
    src = "cleanPatTaus",
    pyModuleName = __name__,
    doSelIndividual = True
)
selectLayer1TausForWTauNu = patTauSelConfiguratorForWTauNu.configure(pyNameSpace = locals())

# loose isolation selection
selectedLayer1TausForWTauNuTaNCdiscrLooseIsolation.cut = cms.string('tauID("byTaNCfrQuarterPercent") > -1.e+3') # cut on TaNC output disabled per default
selectedLayer1TausForWTauNuEcalIsoLooseIsolation.cut = cms.string("isolationPFChargedHadrCandsPtSum()<10")
selectedLayer1TausForWTauNuTrkIsoLooseIsolation.cut = cms.string("isolationPFChargedHadrCandsPtSum()<10")
selectedLayer1TausForWTauNuProngLooseIsolation.cut = selectedLayer1TausForWTauNuTrkIsoLooseIsolation.cut
selectedLayer1TausForWTauNuChargeLooseIsolation.cut = selectedLayer1TausForWTauNuTrkIsoLooseIsolation.cut
selectedLayer1TausForWTauNuMuonVetoLooseIsolation.cut = selectedLayer1TausForWTauNuMuonVeto.cut
selectedLayer1TausForWTauNuElectronVetoLooseIsolation.cut = selectedLayer1TausForWTauNuElectronVeto.cut
selectedLayer1TausForWTauNuEcalCrackVetoLooseIsolation.cut = selectedLayer1TausForWTauNuEcalCrackVeto.cut

patTauSelConfiguratorForWTauNuLooseIsolation = objSelConfigurator(
    [ selectedLayer1TausForWTauNuEta21,
      selectedLayer1TausForWTauNuPt20,
      selectedLayer1TausForWTauNuLeadTrk,
      selectedLayer1TausForWTauNuLeadTrkPt,
      selectedLayer1TausForWTauNuTaNCdiscrLooseIsolation,
      selectedLayer1TausForWTauNuEcalIsoLooseIsolation,
      selectedLayer1TausForWTauNuTrkIsoLooseIsolation,
      selectedLayer1TausForWTauNuProngLooseIsolation,
      selectedLayer1TausForWTauNuChargeLooseIsolation, 
      selectedLayer1TausForWTauNuMuonVetoLooseIsolation,
      selectedLayer1TausForWTauNuElectronVetoLooseIsolation,
      selectedLayer1TausForWTauNuEcalCrackVetoLooseIsolation ],
    src = "cleanPatTaus",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectLayer1TausForWTauNuLooseIsolation = patTauSelConfiguratorForWTauNuLooseIsolation.configure(pyNameSpace = locals())

produceLayer1SelLeptons = cms.Sequence (
    selectLayer1Electrons + selectLayer1ElectronsLooseIsolation
   + selectLayer1Muons + selectLayer1MuonsLooseIsolation
   + selectLayer1Taus
   + selectLayer1ElectronsForElecTau + selectLayer1ElectronsForElecTauLooseIsolation
   + selectLayer1ElectronsForElecMu + selectLayer1ElectronsForElecMuLooseIsolation
   + selectLayer1TausForElecTau + selectLayer1TausForMuTau + selectLayer1TausForDiTau
   + selectLayer1TausForWTauNu +selectLayer1TausForWTauNuLooseIsolation
)
