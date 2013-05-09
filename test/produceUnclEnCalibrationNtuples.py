#!/usr/bin/env python

from TauAnalysis.RecoTools.tools.buildConfigFilesZllRecoilCorrectionAnalysis import getPATtupleFileNames

import os

samplesToAnalyze = [
    'Data_runs190456to193621',
    'Data_runs193834to196531',
    'Data_runs190782to190949_recover',
    'Data_runs198022to198523',
    'Data_runs198934to202016',
    'Data_runs202044to203002',
    'Data_runs203894to208686',
    'ZplusJets_madgraph'
]

version = 'v9_07'

numInputFilesPerJob = 5

configFile = 'produceUnclEnCalibrationNtuple_cfg.py'
executable_cmsRun = 'cmsRun'
executable_rm = 'rm -f'

nice = 'nice '

inputFilePath = "/data2/veelken/CMSSW_5_3_x/PATtuples/ZllRecoilCorrection/%s/" % version
outputFilePath = "/data2/veelken/CMSSW_5_3_x/Ntuples/unclEnCalibration/Ntuples/%s/" % version

# create outputFilePath in case it does not yet exist
def createFilePath_recursively(filePath):
    filePath_items = filePath.split('/')
    currentFilePath = "/"
    for filePath_item in filePath_items:
        currentFilePath = os.path.join(currentFilePath, filePath_item)
        if len(currentFilePath) <= 1:
            continue
        if not os.path.exists(currentFilePath):
            os.mkdir(currentFilePath)

if not os.path.isdir(outputFilePath):
    print "outputFilePath does not yet exist, creating it."
    createFilePath_recursively(outputFilePath)

if not os.path.isdir("tmp"):
    os.mkdir('tmp')

configFilePath = os.path.join(os.getcwd(), "tmp")

# Function to prepare customized config files
def customizeConfigFile(sample, version, inputFileNames, outputFileName, cfgFileName_original, cfgFileName_modified):
    cfgFile_original = open(cfgFileName_original, "r")
    cfg_original = cfgFile_original.read()
    cfgFile_original.close()

    cfg_modified = cfg_original.replace("#__", "")
    cfg_modified = cfg_modified.replace("#sample#", "'%s'" % sample)
    cfg_modified = cfg_modified.replace("#version#", "'%s'" % version)
    cfg_modified = cfg_modified.replace("#outputFileName#", "'%s'" % outputFileName)

    cfg_modified += "\n"
    cfg_modified += "process.source.fileNames = cms.untracked.vstring()\n"
    for inputFileName in inputFileNames:
        cfg_modified += "process.source.fileNames.append('file:%s')\n" % inputFileName
    cfg_modified += "\n"
    cfg_modified += "process.TFileService.fileName = cms.string('%s')\n" % outputFileName
    cfg_modified += "\n"

    cfgFile_modified = open(cfgFileName_modified, "w")
    cfgFile_modified.write(cfg_modified)
    cfgFile_modified.close()

# Split list of inputFileNames into groups of size N
def chunks(inputFileNames, N):
    for idx in range(0, len(inputFileNames), N):
        yield inputFileNames[idx:idx + N]

#--------------------------------------------------------------------------------
# Build config files

cfgFileNames = {} # key = outputFileName
logFileNames = {} # key = outputFileName
for sampleToAnalyze in samplesToAnalyze:

    inputFileNames_matched = getPATtupleFileNames( [ sampleToAnalyze ], inputFilePath)[0]
    inputFileNames_matched = sorted(inputFileNames_matched)
    if len(inputFileNames_matched) == 0:
        print("Sample %s has no input files --> skipping !!" % sampleToAnalyze)
        continue
    
    for jobId, inputFileNames_chunk in enumerate(chunks(inputFileNames_matched, numInputFilesPerJob)):
        print "building config file for: sample = %s, part = %i" % (sampleToAnalyze, jobId)
                        
        cfgFileName_original = configFile
        cfgFileName_modified = os.path.join(configFilePath, cfgFileName_original.replace("_cfg.py", "_%s_%s_part%i_cfg.py" % (sampleToAnalyze, version, jobId)))
        print " cfgFileName_modified = '%s'" % cfgFileName_modified

        outputFileName = os.path.join(outputFilePath, "unclEnCalibrationNtuple_%s_%s_2013May08_part%i.root" % (sampleToAnalyze, version, jobId))
        print " outputFileName = '%s'" % outputFileName

        customizeConfigFile(sampleToAnalyze, version, inputFileNames_chunk, outputFileName, cfgFileName_original, cfgFileName_modified)
                
        logFileName = outputFileName
        logFileName = logFileName.replace(".root", ".out")
                
        cfgFileNames[outputFileName] = cfgFileName_modified
        logFileNames[outputFileName] = logFileName

def make_MakeFile_vstring(list_of_strings):
    retVal = ""
    for i, string_i in enumerate(list_of_strings):
        if i > 0:
            retVal += " "
        retVal += string_i
    return retVal

# done building config files, now build Makefile...
makeFileName = os.path.join(configFilePath, "Makefile_produceUnclEnCalibrationNtuples_%s" % version)
makeFile = open(makeFileName, "w")
makeFile.write("\n")
outputFileNames = []
for outputFileName in cfgFileNames.keys():
    outputFileNames.append(outputFileName)
makeFile.write("all: %s\n" % make_MakeFile_vstring(outputFileNames))
makeFile.write("\techo 'Finished producing Ntuples for unclustered energy Calibration.'\n")
makeFile.write("\n")
for outputFileName, cfgFileName in cfgFileNames.items():
    makeFile.write("%s: %s\n" %
      (outputFileName,
       ""))
    makeFile.write("\t%s%s %s &> %s\n" %
      (nice, executable_cmsRun,
       cfgFileName,
       logFileNames[outputFileName]))
makeFile.write("\n")
makeFile.write(".PHONY: clean\n")
makeFile.write("clean:\n")
makeFile.write("\t%s %s\n" % (executable_rm, make_MakeFile_vstring(outputFileNames)))
makeFile.write("\techo 'Finished deleting old files.'\n")
makeFile.write("\n")
makeFile.close()

print("Finished building Makefile. Now execute 'make -j 12 -f %s'." % makeFileName)
