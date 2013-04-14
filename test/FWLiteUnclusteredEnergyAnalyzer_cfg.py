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

##sample = 'Data_runs203894to208686'
##shiftBy = 0.
##branchNames_weights = []
##applyZvtxReweight = False

version = 'v9_04'

# CaloTowers
##type = "caloTowers"
##branchName_met = 'caloMEt'
##branchName_unclEnSum = 'sumCaloTowers'
##subtract_qT = False
##etaMin = -3.0
##etaMax = +3.0

# PFCandidates
type = "pfCands"
branchName_met = 'pfMEt'
branchName_unclEnSum = 'sumPFCands'
subtract_qT = True
etaMin = -9.9
etaMax = +9.9
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system/grid
#
#__sample = 'SAMPLE_HOOK'
#__type = 'TYPE_HOOK'
#__shiftBy = SHIFTBY_HOOK
#
if sample == 'ZplusJets_madgraph':
    branchNames_weights = [ 'vertexMultiplicityReweight3d2012RunDruns203894to208686_' ]
    applyZvtxReweight = True
    ##applyZvtxReweight = False
elif sample == 'Data_runs203894to208686':
    branchNames_weights = []
    applyZvtxReweight = False
else:
    raise ValueError("Invalid Configuration Parameter 'sample' = %s !!" % sample)
if type == "caloTowersNoHF":
    branchName_met = 'caloMEtNoHF'
    branchName_unclEnSum = 'sumCaloTowersNoHF'
    subtract_qT = False
    etaMin = -3.0
    etaMax = +3.0
elif type == "caloTowers":
    branchName_met = 'caloMEt'
    branchName_unclEnSum = 'sumCaloTowers'
    subtract_qT = False
    etaMin = -9.9
    etaMax = +9.9
elif type == "pfCands":
    branchName_met = 'pfMEt'
    branchName_unclEnSum = 'sumPFCands'
    subtract_qT = True
    etaMin = -9.9
    etaMax = +9.9
elif type == "tracks":
    branchName_met = 'trackMEt'
    branchName_unclEnSum = 'sumTracks'
    subtract_qT = True
    etaMin = -9.9
    etaMax = +9.9
elif type == "genParticles":
    branchName_met = 'genMEt'
    branchName_unclEnSum = 'sumGenParticles'
    subtract_qT = True
    etaMin = -9.9
    etaMax = +9.9    
else:    
    raise ValueError("Invalid Configuration Parameter 'type' = %s !!" % type)  
#--------------------------------------------------------------------------------

etaBinsForResidualCorr = [
    -5.191, -3.489, -3.139, -2.964, -2.853, -2.5, -2.411, -2.322, -1.93, -1.479, -1.305, -1.131, -0.957, -0.783, -0.522, -0.261, 0.,
    0.261, 0.522, 0.783, 0.957, 1.131, 1.305, 1.479, 1.93, 2.322, 2.411, 2.5, 2.853, 2.964, 3.139, 3.489, 5.191
]
numEtaBinsForResidualCorr = len(etaBinsForResidualCorr) - 1
if numEtaBinsForResidualCorr != 32:
    raise ValueError("Invalid binning = ", etaBinsForResidualCorr)

binning = []
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
    binning.append(cfgBinnining)

inputFilePath = '/data2/veelken/CMSSW_5_3_x/Ntuples/unclEnCalibration/'
inputFile_regex = r"[a-zA-Z0-9_/:.]*unclEnCalibrationNtuple_%s_%s_2012Mar31.root" % (sample, version)

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

label = "%s_%s" % (sample, type)
if abs(shiftBy) < 1.e-3:
    label += "_central"
elif shiftBy > +1.e-3:
    label += "_shiftUp"
elif shiftBy < -1.e-3:
    label += "_shiftDown"
print "label = %s" % label
outputFileName = 'UnclusteredEnergyAnalyzer_%s.root' % label
process.fwliteOutput = cms.PSet(
    fileName = cms.string(outputFileName)
)

process.UnclusteredEnergyAnalyzer = cms.PSet(
    type = cms.string(type),
    treeName = cms.string("unclEnCalibrationNtupleProducer/unclEnCalibrationNtuple"),
    branchName_recZ = cms.string('recZ'),
    branchName_recMuPlus = cms.string('recMuPlus'),
    branchName_recMuMinus = cms.string('recMuMinus'),
    branchName_met = cms.string(branchName_met),
    branchName_unclEnSum = cms.string(branchName_unclEnSum),
    branchName_numVertices = cms.string('numVertices'),
    branchName_zVtx = cms.string('theVertexZ'),
    branchNames_weights = cms.vstring(branchNames_weights),
    applyZvtxReweight = cms.bool(applyZvtxReweight),
    zVtxReweightInputFileName = cms.FileInPath("TauAnalysis/RecoTools/data/zVtxReweight_runs203894to208686_vs_Summer12mc.root"),
    zVtxReweightHistogramName = cms.string('zVtxReweight'),
    subtract_qT = cms.bool(subtract_qT),
    shiftBy = cms.double(shiftBy),
    binning = cms.VPSet(binning),
    directory = cms.string(label)
)
