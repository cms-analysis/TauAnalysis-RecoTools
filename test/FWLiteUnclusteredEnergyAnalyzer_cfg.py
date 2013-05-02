import FWCore.ParameterSet.Config as cms

import TauAnalysis.Configuration.tools.castor as castor
import TauAnalysis.Configuration.tools.eos as eos

import os
import re

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring(),
    maxEvents = cms.int32(-1)
)

#--------------------------------------------------------------------------------
# define configuration parameter default values

sample = 'ZplusJets_madgraph'
#shiftBy = -0.10
shiftBy = 0.
#shiftBy = +0.10
branchNames_weights = [ 'vertexMultiplicityReweight3d2012RunDruns203894to208686_' ] 
applyZvtxReweight = True
applyQtReweight = True

version = 'v9_07_2013May01'

# PFCandidates
metType = "pfCands"
branchName_met = 'pfMEt'
branchName_jetSum = ''
branchName_unclEnSum = 'sumPFCands'
subtract_qT = True
etaMin = -9.9
etaMax = +9.9

##numEtaBinsForResidualCorr = 28
numEtaBinsForResidualCorr = 32

outputFileName = None
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system/grid
#
#__sample = #sample#
#__metType = #metType#
#__version = #version#
#__shiftBy = #shiftBy#
#__type1JetPtThreshold = #type1JetPtThreshold#
#__applyResidualCorr = #applyResidualCorr#
#__outputFileName = #outputFileName#
#
if sample == 'ZplusJets_madgraph':
    ##branchNames_weights = [ 'vertexMultiplicityReweight3d2012RunDruns203894to208686_' ]
    branchNames_weights = []
    applyZvtxReweight = True
    ##applyZvtxReweight = False
    applyQtReweight = True
    numPileUpBinning = [ 0., 12.5, 15., 17.5, 20., 22.5, 1.e+3 ]
elif sample == 'Data_runs190456to193621' or \
     sample == 'Data_runs193834to196531' or \
     sample == 'Data_runs190782to190949_recover' or \
     sample == 'Data_runs198022to198523' or \
     sample == 'Data_runs198934to202016' or \
     sample == 'Data_runs202044to203002' or \
     sample == 'Data_runs203894to208686':
    branchNames_weights = []
    applyZvtxReweight = False
    applyQtReweight = False
    numPileUpBinning = [ 0., 12.5, 15., 17.5, 20., 22.5, 1.e+3 ]
else:
    raise ValueError("Invalid Configuration Parameter 'sample' = %s !!" % sample)
if metType == "caloTowersNoHF":
    branchName_met = 'caloMEtNoHF'
    branchName_jetSum = ''
    branchName_unclEnSum = 'caloTowersNoHF'
    subtract_qT = False
    etaMin = -3.0
    etaMax = +3.0
    residualCorrFileName = "JetMETCorrections/Type1MET/data/unclEnResidualCorr_2012RunABCDruns190456to208686_caloTowersNoHF.txt"
elif metType == "caloTowers":
    branchName_met = 'caloMEt'
    branchName_jetSum = ''
    branchName_unclEnSum = 'caloTowers'
    subtract_qT = False
    etaMin = -9.9
    etaMax = +9.9
    residualCorrFileName = "JetMETCorrections/Type1MET/data/unclEnResidualCorr_2012RunABCDruns190456to208686_caloTowers.txt"
elif metType == "pfCands":
    branchName_met = 'pfMEt'
    branchName_jetSum = ''
    branchName_unclEnSum = 'sumPFCands'
    subtract_qT = True
    etaMin = -9.9
    etaMax = +9.9
    residualCorrFileName = "JetMETCorrections/Type1MET/data/unclEnResidualCorr_2012RunABCDruns190456to208686_pfCands.txt"
elif metType == "pfJetsPlusCandsNotInJet":
    branchName_met = 'pfMEt'
    branchName_jetSum = 'sumPFJetsPtGt%1.0f' % type1JetPtThreshold
    branchName_unclEnSum = 'sumPFJetsPtLt%1.0fPlusCandsNotInJet' % type1JetPtThreshold
    subtract_qT = True
    etaMin = -9.9
    etaMax = +9.9
    residualCorrFileName = "JetMETCorrections/Type1MET/data/unclEnResidualCorr_2012RunABCDruns190456to208686_pfCands.txt"
elif metType == "tracks":
    branchName_met = 'trackMEt'
    branchName_jetSum = ''
    branchName_unclEnSum = 'tracks'
    subtract_qT = True
    etaMin = -9.9
    etaMax = +9.9
    residualCorrFileName = "JetMETCorrections/Type1MET/data/unclEnResidualCorr_2012RunABCDruns190456to208686_tracks.txt"
elif metType == "genParticles":
    branchName_met = 'genMEt'
    branchName_jetSum = ''
    branchName_unclEnSum = 'sumGenParticles'
    subtract_qT = True
    etaMin = -9.9
    etaMax = +9.9    
else:    
    raise ValueError("Invalid Configuration Parameter 'metType' = %s !!" % metType)  
#--------------------------------------------------------------------------------

etaBinsForResidualCorr = None
if numEtaBinsForResidualCorr == 28:
    etaBinsForResidualCorr = [
        -5.191, -2.964, -2.853, -2.5, -2.411, -2.322, -1.93, -1.479, -1.305, -1.131, -0.957, -0.783, -0.522, -0.261, 0.,
        0.261, 0.522, 0.783, 0.957, 1.131, 1.305, 1.479, 1.93, 2.322, 2.411, 2.5, 2.853, 2.964, 5.191
    ]
elif numEtaBinsForResidualCorr == 32:
    etaBinsForResidualCorr = [
        -5.191, -3.489, -3.139, -2.964, -2.853, -2.5, -2.411, -2.322, -1.93, -1.479, -1.305, -1.131, -0.957, -0.783, -0.522, -0.261, 0.,
        0.261, 0.522, 0.783, 0.957, 1.131, 1.305, 1.479, 1.93, 2.322, 2.411, 2.5, 2.853, 2.964, 3.139, 3.489, 5.191
    ]
else:
    raise ValueError("Invalid Configuration Parameter 'numEtaBinsForResidualCorr' = %i" % numEtaBinsForResidualCorr)
if not (len(etaBinsForResidualCorr) - 1) == numEtaBinsForResidualCorr:
    raise ValueError("Invalid binning (%i bins expected) = ", etaBinsForResidualCorr)

etaBinning = []
for idxEtaBinForResidualCorr in range(numEtaBinsForResidualCorr):
    binEdgeLow = etaBinsForResidualCorr[idxEtaBinForResidualCorr]
    binEdgeHigh = etaBinsForResidualCorr[idxEtaBinForResidualCorr + 1]
    binLabel = "eta%1.3ftoeta%1.3f" % (binEdgeLow, binEdgeHigh)
    binLabel = binLabel.replace(".", "p")
    binLabel = binLabel.replace("+", "P")
    binLabel = binLabel.replace("-", "M")
    binCenter = 0.5*(binEdgeLow + binEdgeHigh)
    if not (binCenter > etaMin and binCenter < etaMax):
        continue
    cfgBinnining = cms.PSet(
        binLabel = cms.string(binLabel),
        binCenter = cms.double(binCenter),
        binEdgeLow = cms.double(binEdgeLow),
        binEdgeHigh = cms.double(binEdgeHigh)
    )
    etaBinning.append(cfgBinnining)

inputFilePath = '/data2/veelken/CMSSW_5_3_x/Ntuples/unclEnCalibration/'
inputFile_regex = r"[a-zA-Z0-9_/:.]*unclEnCalibrationNtuple_%s_%s.root" % (sample, version)

#--------------------------------------------------------------------------------
inputFileNames = []
if inputFilePath.find('/castor/') != -1:
    inputFileNames = [ 'rfio:%s' % file_info['path'] for file_info in castor.nslsl(inputFilePath) ]
elif inputFilePath.find("/store") != -1:
    inputFileNames = [ file_info['path'] for file_info in eos.lsl(inputFilePath) ]    
else:
    inputFileNames = [ 'file:%s' % os.path.join(inputFilePath, file_name) for file_name in os.listdir(inputFilePath) ]
#print "inputFileNames = %s" % inputFileNames

inputFileNames_matched = []
for inputFileName in inputFileNames:
    inputFile_matcher = re.compile(inputFile_regex)
    if inputFile_matcher.match(inputFileName):
	inputFileNames_matched.append(inputFileName)
#print "inputFileNames_matched = %s" % inputFileNames_matched

if len(inputFileNames_matched) == 0:
    raise ValueError("Found no input files !!")

setattr(process.fwliteInput, "fileNames", cms.vstring(inputFileNames_matched))
#--------------------------------------------------------------------------------

directory = None
if sample.find("Data") != -1:
    directory = 'Data'
else:
    directory = sample
print "directory = %s" % directory
    
label = ""
if applyResidualCorr:
    label += "_wCalibration"
else:    
    if abs(shiftBy) < 1.e-3:
        label += "_central"
    elif shiftBy > +1.e-3:
        label += "_shiftUp"
    elif shiftBy < -1.e-3:
        label += "_shiftDown"
print "label = %s" % label
outputFileName = 'UnclusteredEnergyAnalyzer_%s_%s_%s.root' % (sample, metType, label)
process.fwliteOutput = cms.PSet(
    fileName = cms.string(outputFileName)
)

process.UnclusteredEnergyAnalyzer = cms.PSet(
    metType = cms.string(metType),
    treeName = cms.string("unclEnCalibrationNtupleProducer/unclEnCalibrationNtuple"),
    branchName_recZ = cms.string('recZ'),
    branchName_recMuPlus = cms.string('recMuPlus'),
    branchName_recMuMinus = cms.string('recMuMinus'),
    branchName_met = cms.string(branchName_met),
    branchName_jetSum = cms.string(branchName_jetSum),
    branchName_unclEnSum = cms.string(branchName_unclEnSum),
    branchName_numVertices = cms.string('numVertices'),
    branchName_zVtx = cms.string('theVertexZ'),
    branchName_numPileUp = cms.string('numPileUp'),
    branchNames_weights = cms.vstring(branchNames_weights),
    applyZvtxReweight = cms.bool(applyZvtxReweight),
    zVtxReweightInputFileName = cms.FileInPath("TauAnalysis/RecoTools/data/unclEnCalibration_zVtxReweight_runs203894to208686_vs_Summer12mc.root"),
    zVtxReweightHistogramName = cms.string('zVtxReweight'),
    applyQtReweight = cms.bool(applyQtReweight),
    qTreweightInputFileName = cms.FileInPath("TauAnalysis/RecoTools/data/unclEnCalibration_qTreweight.root"),
    qTreweightHistogramName = cms.string('qTreweight'),
    subtract_qT = cms.bool(subtract_qT),
    shiftBy = cms.double(shiftBy),
    etaBinning = cms.VPSet(etaBinning),
    numPileUpBinning = cms.vdouble(numPileUpBinning),
    applyResidualCorr = cms.bool(applyResidualCorr),
    isMC = cms.bool(False), # CV: only used to decide whether to apply "unclustered energy" calibration to MC or Data
    directory = cms.string("%s_%s_%s" % (directory, metType, label))
)

if applyResidualCorr:
    process.UnclusteredEnergyAnalyzer.residualCorrFileName = cms.FileInPath(residualCorrFileName)
