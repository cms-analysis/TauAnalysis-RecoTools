import FWCore.ParameterSet.Config as cms

import os
import re

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
# define configuration parameter default values

sample = 'ZplusJets_madgraph'
##sample = 'Data_runs203894to208686'

version = 'v9_07'
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system/grid
#
#__sample = #sample#
#__version = #version#
#
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# customize configuration depending on whether we run on Monte Carlo or Data

isMC = None
if sample.find("Data") != -1:
    isMC = False
else:
    isMC = True

srcGenParticles = None
srcWeights = None
if isMC:
    srcGenParticles = 'genParticles'
    srcWeights = [
        'vertexMultiplicityReweight3d2012RunABCDruns190456to208686',
        'vertexMultiplicityReweight3d2012RunDruns203894to208686'
    ]
else:
    srcGenParticles = ''
    srcWeights = []
#--------------------------------------------------------------------------------    

#--------------------------------------------------------------------------------
# set input files

inputFilePath = '/data2/veelken/CMSSW_5_3_x/PATtuples/ZllRecoilCorrection/%s/' % version
inputFile_regex = r"[a-zA-Z0-9_/:.]*ZllRecoilCorrectionPATtuple_%s_%s_(?P<jobId>\d*).root" % (sample, version)

inputFileNames = []
files = None
if inputFilePath.startswith('/castor/'):
    files = [ "".join([ "rfio:", file_info['path'] ]) for file_info in castor.nslsl(inputFilePath) ]
elif inputFilePath.find("/store") != -1:
    files = [ file_info['path'] for file_info in eos.lsl(inputFilePath) ]
else:    
    files = [ "".join([ "file:", inputFilePath, file ]) for file in os.listdir(inputFilePath) ]
for file in files:
    #print "file = %s" % file 
    inputFile_matcher = re.compile(inputFile_regex)
    if inputFile_matcher.match(file):
        inputFileNames.append(file)
#print "inputFileNames = %s" % inputFileNames

##inputFileNames = [ "file:/data1/veelken/CMSSW_5_3_x/skims/ZllRecoilCorrectionPATtuple_ZplusJets_madgraph_qTgt80_v9_07.root" ]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputFileNames),
    ##eventsToProcess = cms.untracked.VEventRange(
    ##    '1:1:115392'
    ##),                        
    skipEvents = cms.untracked.uint32(0),    
)
#--------------------------------------------------------------------------------

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#--------------------------------------------------------------------------------
# define GlobalTag to be used for event reconstruction
if isMC:
    process.GlobalTag.globaltag = cms.string('START53_V15::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_P_V42_AN3::All')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce collections of objects not stored in PAT-tuple

process.reconstructionSequence = cms.Sequence()

process.load("JetMETCorrections/Configuration/DefaultJEC_cff")

process.load("JetMETCorrections/Type1MET/pfMETCorrections_cff")
process.pfCandsNotInJet.name = cms.untracked.string("pfCandsNotInJet")
process.reconstructionSequence += process.pfCandsNotInJet

process.caloTowersNotInJet = cms.EDProducer("TPCaloJetsOnCaloTowers",
    enable = cms.bool(True),
    name = cms.untracked.string("caloTowersNotInJet"),
    topCollection = cms.InputTag('ak5CaloJets'),
    bottomCollection = cms.InputTag('towerMaker'),
    verbose = cms.untracked.bool(False)
)
process.reconstructionSequence += process.caloTowersNotInJet
#--------------------------------------------------------------------------------

process.analysisSequence = cms.Sequence()

def getBinning(etaBinsForResidualCorr):
    binning = []
    numEtaBinsForResidualCorr = len(etaBinsForResidualCorr) - 1
    for idxEtaBinForResidualCorr in range(numEtaBinsForResidualCorr):
        binEdgeLow = etaBinsForResidualCorr[idxEtaBinForResidualCorr]
        binEdgeHigh = etaBinsForResidualCorr[idxEtaBinForResidualCorr + 1]
        binLabel = "eta%1.3ftoeta%1.3f" % (binEdgeLow, binEdgeHigh)
        binLabel = binLabel.replace(".", "p")
        binLabel = binLabel.replace("+", "P")
        binLabel = binLabel.replace("-", "M")
        binSelection = "eta >= %f & eta < %f" % (binEdgeLow, binEdgeHigh)
        cfgBinnining = cms.PSet(
            binLabel = cms.string(binLabel),
            binSelection = cms.string(binSelection)
        )
        binning.append(cfgBinnining)
    return binning

etaBinsForResidualCorr28bins = [
    -5.191, -2.964, -2.853, -2.5, -2.411, -2.322, -1.93, -1.479, -1.305, -1.131, -0.957, -0.783, -0.522, -0.261, 0.,
    0.261, 0.522, 0.783, 0.957, 1.131, 1.305, 1.479, 1.93, 2.322, 2.411, 2.5, 2.853, 2.964, 5.191
]
numEtaBinsForResidualCorr28bins = len(etaBinsForResidualCorr28bins) - 1
if numEtaBinsForResidualCorr28bins != 28:
    raise ValueError("Invalid binning (28 bins expected) = ", etaBinsForResidualCorr28bins)
binning28bins = getBinning(etaBinsForResidualCorr28bins)

etaBinsForResidualCorr32bins = [
    -5.191, -3.489, -3.139, -2.964, -2.853, -2.5, -2.411, -2.322, -1.93, -1.479, -1.305, -1.131, -0.957, -0.783, -0.522, -0.261, 0.,
    0.261, 0.522, 0.783, 0.957, 1.131, 1.305, 1.479, 1.93, 2.322, 2.411, 2.5, 2.853, 2.964, 3.139, 3.489, 5.191
]
numEtaBinsForResidualCorr32bins = len(etaBinsForResidualCorr32bins) - 1
if numEtaBinsForResidualCorr32bins != 32:
    raise ValueError("Invalid binning (32 bins expected) = ", etaBinsForResidualCorr32bins)
binning32bins = getBinning(etaBinsForResidualCorr32bins)

caloJetCorrections = None
pfJetCorrections   = None
if isMC:
    caloJetCorrections = [ 'ak5CaloL1Offset', 'ak5CaloL1Fastjet', 'ak5CaloL1L2L3', 'ak5CaloL1FastL2L3' ]
    pfJetCorrections   = [ 'ak5PFL1Offset', 'ak5PFL1Fastjet', 'ak5PFL1L2L3', 'ak5PFL1FastL2L3' ]
else:
    caloJetCorrections = [ 'ak5CaloL1Offset', 'ak5CaloL1Fastjet', 'ak5CaloL1L2L3Residual', 'ak5CaloL1FastL2L3Residual' ]
    pfJetCorrections   = [ 'ak5PFL1Offset', 'ak5PFL1Fastjet', 'ak5PFL1L2L3Residual', 'ak5PFL1FastL2L3Residual' ]
    
process.unclEnCalibrationNtupleProducer28bins = cms.EDAnalyzer("UnclEnCalibrationNtupleProducer",
    treeName = cms.string("unclEnCalibrationNtuple28"),                                                           
    srcZllCandidates = cms.InputTag('goldenZmumuCandidatesGe1IsoMuons'),
    srcGenParticles = cms.InputTag(srcGenParticles),
    #--------------------------------------------------------                                                       
    srcCaloMEt = cms.InputTag('patType1CorrectedCaloMet'),
    srcCaloJets = cms.InputTag('ak5CaloJets'),
    caloJetCorrections = cms.vstring(caloJetCorrections),
    caloJetCorrectionRef = cms.string(caloJetCorrections[2]), # "ak5CaloL1L2L3" for MC / "ak5CaloL1L2L3Residual" for Data
    srcCaloTowers = cms.InputTag('towerMaker'),
    srcCaloTowersNotInJets = cms.InputTag('caloTowersNotInJet'),
    #--------------------------------------------------------                                
    srcPFMEt = cms.InputTag('patType1CorrectedPFMet'),
    srcPFJets = cms.InputTag('ak5PFJets'),
    pfJetCorrections = cms.vstring(pfJetCorrections),
    pfJetCorrectionRef = cms.string(pfJetCorrections[2]), # "ak5PFL1FastL2L3" for MC / "ak5PFL1FastL2L3Residual" for Data
    srcPFCandidates = cms.InputTag('particleFlow'),
    srcPFCandidatesNotInJets = cms.InputTag('pfCandsNotInJet'),                                                       
    #-------------------------------------------------------- 
    type1JetPtThresholds = cms.vdouble(10., 15., 20., 25., 30.),                                                           
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
    fileName = cms.string("/data2/veelken/CMSSW_5_3_x/Ntuples/unclEnCalibration/unclEnCalibrationNtuple_%s_%s_2013May08.root" % (sample, version))
)

process.p = cms.Path(
    process.reconstructionSequence
   + process.analysisSequence
)

processDumpFile = open('produceUnclEnCalibrationNtuple.dump' , 'w')
print >> processDumpFile, process.dumpPython()
