import FWCore.ParameterSet.Config as cms

process = cms.Process("produceUnclEnCalibrationNtuple")

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('Configuration/Geometry/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#--------------------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        ##'file:/data1/veelken/CMSSW_5_3_x/skims/goldenZmumuEvents_ZplusJets_madgraph_2013Feb14_AOD_99_0_MvW.root'
        ##'file:/data1/veelken/CMSSW_5_3_x/skims/goldenZmumuEvents_Data_runs203894to208686_2013Feb14_AOD_12_2_KCO.root'
        'file:/data1/veelken/CMSSW_5_3_x/skims/goldenZmumuEvents_ZplusJets_madgraph_qTgt80_AOD.root'
        ##'file:/data1/veelken/CMSSW_5_3_x/skims/goldenZmumuEvents_Data_runs203894to208686_qTgt80_AOD.root'                        
    ),
    ##eventsToProcess = cms.untracked.VEventRange(
    ##    '1:48430:19356511'
    ##),
    skipEvents = cms.untracked.uint32(0)            
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#--------------------------------------------------------------------------------
# define configuration parameter default values

isMC = True # use for MC
##isMC = False # use for Data
##HLTprocessName = "HLT" # use for 2012 Data
HLTprocessName = "HLT" # use for Spring'12 MC
#type1JetPtThreshold = 20.0 # increased jet Pt threshold to reduce sensitivity of Type 1 corrected MET to pile-up
type1JetPtThreshold = 10.0 # current default value recommended by JetMET POG
#jecUncertaintyTag = "Total"
#jecUncertaintyTag = "SubTotalDataMC"
jecUncertaintyTag = "SubTotalMC"
runPeriod = "2012RunABCD" # use for MET sys. shift correction vs. Nvtx
applyRochesterMuonCorr = True
##applyUnclEnergyResidualCorr = False
applyUnclEnergyResidualCorr = True
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define GlobalTag to be used for event reconstruction
if isMC:
    process.GlobalTag.globaltag = cms.string('START53_V22::All')
else:
    process.GlobalTag.globaltag = cms.string('FT_53_V21_AN4::All')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system/grid
#
#__isMC = #isMC#
#__HLTprocessName = #HLTprocessName#
#__runPeriod = #runPeriod#
#__process.GlobalTag.globaltag = cms.string(#globaltag#)
#
#--------------------------------------------------------------------------------

process.reconstructionSequence = cms.Sequence()

#--------------------------------------------------------------------------------
# load definition of VBTF Z --> mu+ mu- event selection
# (with no isolation cuts applied on one of the two muons)
# load definitions of data-quality filters
process.load("TauAnalysis.TauIdEfficiency.filterDataQuality_cfi")
if isMC:
    process.dataQualityFilters.remove(process.hltPhysicsDeclared)
    process.dataQualityFilters.remove(process.dcsstatus)
    ##process.dataQualityFilters.remove(process.HBHENoiseFilter)
    process.dataQualityFilters.remove(process.hcalLaserEventFilter)
    process.dataQualityFilters.remove(process.hcallasereventfilter2012)
    process.dataQualityFilters.remove(process.ecalLaserCorrFilter)
else:
    # CV: disable DCS status check
    #    (already applied already in skimming stage;
    #     causes many error messages, due to dropped collections)
    process.dataQualityFilters.remove(process.dcsstatus)
process.reconstructionSequence += process.dataQualityFilters

process.load("TauAnalysis.Skimming.goldenZmmSelectionVBTFnoMuonIsolation_cfi")
process.reconstructionSequence += process.goldenZmumuSelectionSequence

if applyRochesterMuonCorr:
    print "Enabling Rochester muon momentum corrections"
    process.patMuonsForGoldenZmmSelectionRochesterMomentumCorr = cms.EDProducer("RochesterCorrPATMuonProducer",
        src = cms.InputTag('patMuonsForGoldenZmmSelection'),
        isMC = cms.bool(isMC)
    )
    process.goldenZmumuSelectionSequence.replace(process.patMuonsForGoldenZmmSelection, process.patMuonsForGoldenZmmSelection*process.patMuonsForGoldenZmmSelectionRochesterMomentumCorr)
    process.goodMuons.src = cms.InputTag('patMuonsForGoldenZmmSelectionRochesterMomentumCorr')
else:
    print "Rochester muon momentum corrections disabled"
    
process.uniqueGoodMuonPairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("goodMuons"),
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(2)                          
)
process.reconstructionSequence += process.uniqueGoodMuonPairFilter # CV: disable only for testing !!!
process.uniqueGoodIsoMuonPairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("goodIsoMuons"),
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(2)                          
)
process.reconstructionSequence += process.uniqueGoodIsoMuonPairFilter # CV: disable only for testing !!!
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce rerun kt6PFJets to produce energy density rho needed for L1FastJet corrections
# (neccessary because kt6PFJet collection was dropped during skimming)
process.load("RecoJets.Configuration.RecoPFJets_cff")
process.kt6PFJets.doAreaFastjet = cms.bool(True)
process.kt6PFJets.doRhoFastjet = cms.bool(True)
process.reconstructionSequence += process.kt6PFJets

# produce rerun kt6CaloJets to produce energy density rho needed for L1FastJet corrections
# (neccessary because kt6CaloJet collection was dropped during skimming)
process.load("RecoJets.Configuration.RecoJets_cff")
process.kt6CaloJets.doAreaFastjet = cms.bool(True)
process.kt6CaloJets.doRhoFastjet = cms.bool(True)
process.reconstructionSequence += process.kt6CaloJets
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# select reco::Vertex collections
process.load("TauAnalysis/RecoTools/recoVertexSelection_cff")

process.reconstructionSequence += process.selectPrimaryVertex
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce pile-up reweighting factors

if isMC:
    process.load("TauAnalysis/RecoTools/vertexMultiplicityReweight_cfi")
    process.vertexMultiplicityReweight3d2012RunABCDruns190456to208686 = process.vertexMultiplicityReweight.clone(
        inputFileName = cms.FileInPath("TauAnalysis/RecoTools/data/expPUpoissonMean_runs190456to208686_Mu17_Mu8.root"),
        type = cms.string("gen3d"),
        mcPeriod = cms.string("Summer12_S10")
    )
    process.reconstructionSequence += process.vertexMultiplicityReweight3d2012RunABCDruns190456to208686
    process.vertexMultiplicityReweight1d2012RunABCDruns190456to208686 = process.vertexMultiplicityReweight.clone(
        inputFileName = cms.FileInPath("TauAnalysis/RecoTools/data/expPUpoissonDist_runs190456to208686_Mu17_Mu8.root"),
        type = cms.string("gen"),
        mcPeriod = cms.string("Summer12_S10")
    )
    process.reconstructionSequence += process.vertexMultiplicityReweight1d2012RunABCDruns190456to208686
    process.vertexMultiplicityReweight3d2012RunDruns203894to208686 = process.vertexMultiplicityReweight3d2012RunABCDruns190456to208686.clone(
        inputFileName = cms.FileInPath("TauAnalysis/RecoTools/data/expPUpoissonMean_runs203894to208686_IsoMu8_eta2p1_LooseIsoPFTau20_L1ETM26.root")
    )
    process.reconstructionSequence += process.vertexMultiplicityReweight3d2012RunDruns203894to208686
    process.vertexMultiplicityReweight1d2012RunDruns203894to208686 = process.vertexMultiplicityReweight1d2012RunABCDruns190456to208686.clone(
        inputFileName = cms.FileInPath("TauAnalysis/RecoTools/data/expPUpoissonDist_runs203894to208686_IsoMu8_eta2p1_LooseIsoPFTau20_L1ETM26.root")
    )
    process.reconstructionSequence += process.vertexMultiplicityReweight1d2012RunDruns203894to208686
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# reconstruct PFMET and CaloMET and produce collections of PFCandidates and CaloTowers not in jets

process.load("JetMETCorrections/Configuration/DefaultJEC_cff")

process.load("JetMETCorrections/Type1MET/pfMETCorrections_cff")
if isMC:
    process.pfJetMETcorr.offsetCorrLabel = cms.string("ak5PFL1Fastjet")
    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")
else:
    process.pfJetMETcorr.offsetCorrLabel = cms.string("ak5PFL1Fastjet")
    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
process.pfCandsNotInJet.name = cms.untracked.string("pfCandsNotInJet")
process.pfCandsNotInJet.bottomCollection = cms.InputTag('particleFlow')
process.pfCandMETcorr.verbosity = cms.int32(0)
process.pfJetMETcorr.verbosity = cms.int32(0)
process.pfType1CorrectedMet.verbosity = cms.int32(0)
process.reconstructionSequence += process.producePFMETCorrections

process.load("JetMETCorrections/Type1MET/caloMETCorrections_cff")
if isMC:    
    ##process.caloJetMETcorr.offsetCorrLabel = cms.string("ak5CaloL1Offset")
    ##process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL1L2L3")
    process.caloJetMETcorr.offsetCorrLabel = cms.string("ak5CaloL1Fastjet")
    process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL1FastL2L3")
else:
    ##process.caloJetMETcorr.offsetCorrLabel = cms.string("ak5CaloL1Offset")
    ##process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL1L2L3Residual")
    process.caloJetMETcorr.offsetCorrLabel = cms.string("ak5CaloL1Fastjet")
    process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL1FastL2L3Residual")
process.reconstructionSequence += process.produceCaloMETCorrections

process.caloTowersNotInJet = cms.EDProducer("TPCaloJetsOnCaloTowers",
    enable = cms.bool(True),
    name = cms.untracked.string("caloTowersNotInJet"),
    topCollection = cms.InputTag('ak5CaloJets'),
    bottomCollection = cms.InputTag('towerMaker'),
    verbose = cms.untracked.bool(False)
)
process.reconstructionSequence += process.caloTowersNotInJet
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce collections of CaloJets/CaloTowers and PFJets/PFCandidates not overlapping with reconstructed muons
##
##process.ak5CaloJetsAntiOverlapWithLeptonsVeto = cms.EDFilter("CaloJetAntiOverlapSelector",
##    src = cms.InputTag('ak5CaloJets'),                                                
##    srcNotToBeFiltered = cms.VInputTag('goodMuons'),
##    dRmin = cms.double(0.20),
##    filter = cms.bool(False)                                           
##)
##process.reconstructionSequence += process.ak5CaloJetsAntiOverlapWithLeptonsVeto
##process.caloTowersAntiOverlapWithLeptonsVeto = cms.EDFilter("CaloTowerAntiOverlapSelector",
##    src = cms.InputTag('towerMaker'),  
##    srcNotToBeFiltered = cms.VInputTag('goodMuons'),
##    dRmin = cms.double(0.10),
##    filter = cms.bool(False)                                           
##)
##process.reconstructionSequence += process.caloTowersAntiOverlapWithLeptonsVeto
##process.caloTowersNotInJetAntiOverlapWithLeptonsVeto = cms.EDFilter("CaloTowerAntiOverlapSelector",
##    src = cms.InputTag('caloTowersNotInJet'),  
##    srcNotToBeFiltered = cms.VInputTag('goodMuons'),
##    dRmin = cms.double(0.10),
##    filter = cms.bool(False)                                           
##)
##process.reconstructionSequence += process.caloTowersNotInJetAntiOverlapWithLeptonsVeto
##
##process.ak5PFJetsAntiOverlapWithLeptonsVeto = cms.EDFilter("PFJetAntiOverlapSelector",
##    src = cms.InputTag('ak5PFJets'),                                                
##    srcNotToBeFiltered = cms.VInputTag('goodMuons'),
##    dRmin = cms.double(0.20),
##    filter = cms.bool(False)                                           
##)
##process.reconstructionSequence += process.ak5PFJetsAntiOverlapWithLeptonsVeto
##process.pfCandsAntiOverlapWithLeptonsVeto = cms.EDFilter("PFCandidateAntiOverlapSelector",
##    src = cms.InputTag('particleFlow'),  
##    srcNotToBeFiltered = cms.VInputTag('goodMuons'),
##    dRmin = cms.double(0.02),
##    filter = cms.bool(False)                                           
##)
##process.reconstructionSequence += process.pfCandsAntiOverlapWithLeptonsVeto
##process.pfCandsNotInJetAntiOverlapWithLeptonsVeto = cms.EDFilter("PFCandidateAntiOverlapSelector",
##    src = cms.InputTag('pfCandsNotInJet'),  
##    srcNotToBeFiltered = cms.VInputTag('goodMuons'),
##    dRmin = cms.double(0.02),
##    filter = cms.bool(False)                                           
##)
##process.reconstructionSequence += process.pfCandsNotInJetAntiOverlapWithLeptonsVeto
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# sum PFCandidates, CaloTowers and reco::Tracks in slices of eta,
# using 2 differents binnings:
#  o 28 bins = binning used for older version of 'Residual' jet energy corrections
#  o 32 bins = binning used for latest 'Residual' corrections
#
process.load("LLRAnalysis/TauTauStudies/sumCaloTowersInEtaSlices_cfi")
process.reconstructionSequence += process.sumCaloTowersInEtaSlicesNoHF28bins
process.reconstructionSequence += process.sumCaloTowersInEtaSlicesNoHF32bins

process.sumCaloTowersInEtaSlices28bins = process.sumCaloTowersInEtaSlicesNoHF28bins.clone(
    noHF = cms.bool(False)
)
process.reconstructionSequence += process.sumCaloTowersInEtaSlices28bins
process.sumCaloTowersInEtaSlices32bins = process.sumCaloTowersInEtaSlicesNoHF32bins.clone(
    noHF = cms.bool(False)
)
process.reconstructionSequence += process.sumCaloTowersInEtaSlices32bins

process.load("JetMETCorrections/Type1MET/sumTracksInEtaSlices_cfi")
process.reconstructionSequence += process.sumTracksInEtaSlices28bins
process.reconstructionSequence += process.sumTracksInEtaSlices32bins
#--------------------------------------------------------------------------------    

process.analysisSequence = cms.Sequence()

from LLRAnalysis.TauTauStudies.sumCaloTowersInEtaSlices_cfi import getBinning, etaBinsForResidualCorr28bins, etaBinsForResidualCorr32bins
binning28bins = getBinning(etaBinsForResidualCorr28bins)
binning32bins = getBinning(etaBinsForResidualCorr32bins)

srcGenParticles    = None 
caloJetCorrections = None
pfJetCorrections   = None
srcWeights         = None
if isMC:
    srcGenParticles = 'genParticles'
    ##caloJetCorrections = [ 'ak5CaloL1Offset', 'ak5CaloL1Fastjet', 'ak5CaloL1L2L3', 'ak5CaloL1FastL2L3' ]
    ##pfJetCorrections   = [ 'ak5PFL1Offset', 'ak5PFL1Fastjet', 'ak5PFL1L2L3', 'ak5PFL1FastL2L3' ]
    caloJetCorrections = [ 'ak5CaloL1Fastjet', 'ak5CaloL1FastL2L3' ]
    pfJetCorrections   = [ 'ak5PFL1Fastjet', 'ak5PFL1FastL2L3' ]
    srcWeights = [
        'vertexMultiplicityReweight3d2012RunABCDruns190456to208686',
        'vertexMultiplicityReweight3d2012RunDruns203894to208686'
    ]
else:
    srcGenParticles = ''
    ##caloJetCorrections = [ 'ak5CaloL1Offset', 'ak5CaloL1Fastjet', 'ak5CaloL1L2L3Residual', 'ak5CaloL1FastL2L3Residual' ]
    ##pfJetCorrections   = [ 'ak5PFL1Offset', 'ak5PFL1Fastjet', 'ak5PFL1L2L3Residual', 'ak5PFL1FastL2L3Residual' ]
    caloJetCorrections = [ 'ak5CaloL1Fastjet', 'ak5CaloL1FastL2L3Residual' ]
    pfJetCorrections   = [ 'ak5PFL1Fastjet', 'ak5PFL1FastL2L3Residual' ]
    srcWeights = []
    
process.unclEnCalibrationNtupleProducer28bins = cms.EDAnalyzer("UnclEnCalibrationNtupleProducer",
    treeName = cms.string("unclEnCalibrationNtuple28"),                                                           
    srcZllCandidates = cms.InputTag('goldenZmumuCandidatesGe1IsoMuons'),
    srcGenParticles = cms.InputTag(srcGenParticles),
    #--------------------------------------------------------                                                       
    srcCaloMEt = cms.InputTag('caloType1CorrectedMet'),
    ##srcCaloJets = cms.InputTag('ak5CaloJetsAntiOverlapWithLeptonsVeto'),
    srcCaloJets = cms.InputTag('ak5CaloJets'), 
    caloJetCorrections = cms.vstring(caloJetCorrections),
    ##caloJetCorrectionRef = cms.string(caloJetCorrections[2]), # "ak5CaloL1L2L3" for MC / "ak5CaloL1L2L3Residual" for Data
    caloJetCorrectionRef = cms.string(caloJetCorrections[1]), # "ak5CaloL1FastL2L3" for MC / "ak5CaloL1FastL2L3Residual" for Data
    srcCaloTowers = cms.InputTag('towerMaker'),
    srcCaloTowersNotInJets = cms.InputTag('caloTowersNotInJet'),                                                           
    ##srcCaloTowers = cms.InputTag('caloTowersAntiOverlapWithLeptonsVeto'),
    ##srcCaloTowersNotInJets = cms.InputTag('caloTowersNotInJetAntiOverlapWithLeptonsVeto'),
    type1CaloJetPtThresholds = cms.vdouble(20., 30.),                                                             
    #--------------------------------------------------------                                
    srcPFMEt = cms.InputTag('pfType1CorrectedMet'),
    ##srcPFJets = cms.InputTag('ak5PFJetsAntiOverlapWithLeptonsVeto'),                                                               
    srcPFJets = cms.InputTag('ak5PFJets'),
    skipMuons = cms.bool(True),
    skipMuonSelection = cms.string("isGlobalMuon | isStandAloneMuon"),                                                           
    pfJetCorrections = cms.vstring(pfJetCorrections),
    ##pfJetCorrectionRef = cms.string(pfJetCorrections[3]), # "ak5PFL1FastL2L3" for MC / "ak5PFL1FastL2L3Residual" for Data
    pfJetCorrectionRef = cms.string(pfJetCorrections[1]), # "ak5PFL1FastL2L3" for MC / "ak5PFL1FastL2L3Residual" for Data
    ##srcPFCandidates = cms.InputTag('pfCandsAntiOverlapWithLeptonsVeto'),
    ##srcPFCandidatesNotInJets = cms.InputTag('pfCandsNotInJetAntiOverlapWithLeptonsVeto'),
    srcPFCandidates = cms.InputTag('particleFlow'),
    srcPFCandidatesNotInJets = cms.InputTag('pfCandsNotInJet'),
    type1PFJetPtThresholds = cms.vdouble(10., 15., 20.),  
    #-------------------------------------------------------- 
    sumsInEtaSlices = cms.PSet(
        caloTowersNoHF = cms.PSet(                                                  
            src = cms.InputTag('sumCaloTowersInEtaSlicesNoHF28bins')
        ),
        caloTowers = cms.PSet(                                                  
            src = cms.InputTag('sumCaloTowersInEtaSlices28bins')
        ),
        tracks = cms.PSet(                                                  
            src = cms.InputTag('sumTracksInEtaSlices28bins')
        )
    ),
    srcVertices = cms.InputTag('selectedPrimaryVertexPosition'),
    #--------------------------------------------------------
    # CV: pile-up information for Monte Carlo and data                                                               
    srcGenPileUpSummary = cms.InputTag('addPileupInfo'),
    inputFileNameLumiCalc = cms.FileInPath('TauAnalysis/RecoTools/data_nocrab/lumiCalc_2012RunABCD_byLS.out'),
    isMC = cms.bool(isMC),
    #--------------------------------------------------------                                                      
    srcWeights = cms.VInputTag(srcWeights),
    binning = cms.VPSet(binning28bins),
    ##binning = cms.VPSet(
    ##    cms.PSet(
    ##        binLabel = cms.string("allEta"),
    ##        binSelection = cms.string("eta >= -1.e+3 & eta < +1.e+3")
    ##    )
    ##),                                                           
    verbosity = cms.int32(0)
)
process.analysisSequence += process.unclEnCalibrationNtupleProducer28bins

process.unclEnCalibrationNtupleProducer32bins = process.unclEnCalibrationNtupleProducer28bins.clone(
    treeName = cms.string("unclEnCalibrationNtuple32"),
    sumsInEtaSlices = cms.PSet(
        caloTowersNoHF = cms.PSet(                                                  
            src = cms.InputTag('sumCaloTowersInEtaSlicesNoHF32bins')
        ),
        caloTowers = cms.PSet(                                                  
            src = cms.InputTag('sumCaloTowersInEtaSlices32bins')
        ),
        tracks = cms.PSet(                                                  
            src = cms.InputTag('sumTracksInEtaSlices32bins')
        )
    ),
    binning = cms.VPSet(binning32bins),
    verbosity = cms.int32(0)
)    
process.analysisSequence += process.unclEnCalibrationNtupleProducer32bins

process.TFileService = cms.Service("TFileService",
    ##fileName = cms.string("/data2/veelken/CMSSW_5_3_x/Ntuples/unclEnCalibration/unclEnCalibrationNtuple_Data_runs203894to208686_qTgt80.root")
    fileName = cms.string("/data2/veelken/CMSSW_5_3_x/Ntuples/unclEnCalibration/unclEnCalibrationNtuple_ZplusJets_madgraph_qTgt80.root")                               
)

process.p = cms.Path(
    process.reconstructionSequence
   + process.analysisSequence
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

processDumpFile = open('produceUnclEnCalibrationNtuple.dump' , 'w')
print >> processDumpFile, process.dumpPython()
