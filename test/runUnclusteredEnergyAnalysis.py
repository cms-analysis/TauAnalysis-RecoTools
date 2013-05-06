#!/usr/bin/env python

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

version = 'v9_07_2013May01'

metTypes_Data = {
    'caloTowers' : {
        ##'jobs'                 : [ 'central', 'wCalibration' ],
        'jobs'                 : [ 'central' ],
        'type1JetPtThresholds' : [ 0. ]
    },
    'caloTowersNoHF' : {
        ##'jobs'                 : [ 'central', 'wCalibration' ],
        'jobs'                 : [ 'central' ],
        'type1JetPtThresholds' : [ 0. ]
    },
    'pfCands' : {
        ##'jobs'                 : [ 'central', 'wCalibration' ],
        'jobs'                 : [ 'central' ],
        'type1JetPtThresholds' : [ 0. ]
    },
    'pfJetsPlusCandsNotInJet' : {
        ##'jobs'                 : [ 'central', 'wCalibration' ],
        'jobs'                 : [ 'central' ],
        'type1JetPtThresholds' : [ 10., 15., 20. ]
    },
    'tracks' : {
        ##'jobs'                 : [ 'central', 'wCalibration' ],
        'jobs'                 : [ 'central' ],
        'type1JetPtThresholds' : [ 0. ]
    }
}

metTypes_MC = {
    'caloTowers' : {
        'jobs'                 : [ 'central', 'shiftUp', 'shiftDown' ],
        'type1JetPtThresholds' : [ 0. ]
    },        
    'caloTowersNoHF' : {
        'jobs'                 : [ 'central', 'shiftUp', 'shiftDown' ],
        'type1JetPtThresholds' : [ 0. ]
    },
    'pfCands' : {
        'jobs'                 : [ 'central', 'shiftUp', 'shiftDown' ],
        'type1JetPtThresholds' : [ 0. ]
    },
    'pfJetsPlusCandsNotInJet' : {
        'jobs'                 : [ 'central', 'shiftUp', 'shiftDown' ],
        'type1JetPtThresholds' : [ 10., 15., 20. ]
    },
    'tracks' : {
        'jobs'                 : [ 'central' ],
        'type1JetPtThresholds' : [ 0. ]
    },
    'genParticles' : {
        'jobs'                 : [ 'central' ],
        'type1JetPtThresholds' : [ 0. ]
    }
}

configFile = 'FWLiteUnclusteredEnergyAnalyzer_cfg.py'
executable_FWLiteUnclusteredEnergyAnalyzer = 'FWLiteUnclusteredEnergyAnalyzer'
executable_hadd = 'hadd'
executable_rm = 'rm'

outputFilePath = "/data2/veelken/CMSSW_5_3_x/Ntuples/unclEnCalibration/plots/%s/" % version

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
def customizeConfigFile(sample, metType, version, shiftBy, type1JetPtThreshold, applyResidualCorr, outputFileName, cfgFileName_original, cfgFileName_modified):
    cfgFile_original = open(cfgFileName_original, "r")
    cfg_original = cfgFile_original.read()
    cfgFile_original.close()

    cfg_modified = cfg_original.replace("#__", "")
    cfg_modified = cfg_modified.replace("#sample#", "'%s'" % sample)
    cfg_modified = cfg_modified.replace("#metType#", "'%s'" % metType)
    cfg_modified = cfg_modified.replace("#version#", "'%s'" % version)
    cfg_modified = cfg_modified.replace("#shiftBy#", "%1.2f" % shiftBy)
    cfg_modified = cfg_modified.replace("#type1JetPtThreshold#", "%1.1f" % type1JetPtThreshold)
    applyResidualCorr_string = None
    if applyResidualCorr:
        applyResidualCorr_string = "True"
    else:
        applyResidualCorr_string = "False"
    cfg_modified = cfg_modified.replace("#applyResidualCorr#", "%s" % applyResidualCorr_string)
    cfg_modified = cfg_modified.replace("#outputFileName#", "'%s'" % outputFileName)
       
    cfgFile_modified = open(cfgFileName_modified, "w")
    cfgFile_modified.write(cfg_modified)
    cfgFile_modified.close()

#--------------------------------------------------------------------------------
# Build config files

commandsToRun = {} # key = metType
outputFileNames = {} # key = metType
for sampleToAnalyze in samplesToAnalyze:
    metTypes = None
    if sampleToAnalyze.find("Data") != -1:
        metTypes = metTypes_Data
    else:
        metTypes = metTypes_MC
    for metType, jobConfig in metTypes.items():
        for job in jobConfig['jobs']:
            for type1JetPtThreshold in jobConfig['type1JetPtThresholds']:
                print "building '%s' config file for: sample = %s, metType = %s, job = %s, type1JetPtThreshold = %1.1f" % \
                  (executable_FWLiteUnclusteredEnergyAnalyzer, sampleToAnalyze, metType, job, type1JetPtThreshold)
                shiftBy = None
                if job.find("shiftUp") != -1:
                    shiftBy = +0.10
                elif job.find("shiftDown") != -1:
                    shiftBy = -0.10
                else:
                    shiftBy = 0.
                
                applyResidualCorr = None
                if job.find("wCalibration") != -1:
                    applyResidualCorr = True
                else:
                    applyResidualCorr = False
                    
                label = "%s_%s_jetPtThreshold%1.1f" % (sampleToAnalyze, metType, type1JetPtThreshold)
                label = label.replace(".", "_")
                if applyResidualCorr:
                    label += "_wCalibration"
                else:    
                    if abs(shiftBy) < 1.e-3:
                        label += "_central"
                    elif shiftBy > +1.e-3:
                        label += "_shiftUp"
                    elif shiftBy < -1.e-3:
                        label += "_shiftDown"
                        
                cfgFileName_original = configFile
                cfgFileName_modified = os.path.join(configFilePath, cfgFileName_original.replace("_cfg.py", "_%s_cfg.py" % label))
                print " cfgFileName_modified = '%s'" % cfgFileName_modified
                
                outputFileName = cfgFileName_modified
                outputFileName = outputFileName.replace(configFilePath, "")
                if outputFileName.find("/") == 0:
                    outputFileName = outputFileName[1:]
                outputFileName = outputFileName.replace("_cfg.py", ".root")
                outputFileName = outputFileName.replace("FWLiteUnclusteredEnergyAnalyzer", "UnclusteredEnergyAnalyzer")
                outputFileName = os.path.join(outputFilePath, outputFileName)
                print " outputFileName = '%s'" % outputFileName

                customizeConfigFile(sampleToAnalyze, metType, version, shiftBy, type1JetPtThreshold, applyResidualCorr, outputFileName, cfgFileName_original, cfgFileName_modified)
                
                logFileName = outputFileName
                logFileName = logFileName.replace(".root", ".out")
                
                if not metType in commandsToRun.keys():
                    commandsToRun[metType] = []
                command = 'nice %s %s >&! %s' % (executable_FWLiteUnclusteredEnergyAnalyzer, cfgFileName_modified, logFileName)
                print "--> command = '%s'" % command
                commandsToRun[metType].append(command)

                if not metType in outputFileNames.keys():
                    outputFileNames[metType] = []
                outputFileNames[metType].append(outputFileName)

def make_vstring(list_of_strings):
    retVal = ""
    for i, string_i in enumerate(list_of_strings):
        if i > 0:
            retVal += " "
        retVal += string_i
    return retVal

for metType in outputFileNames.keys():
    print "configuring new '%s' for: metType = %s" % (executable_hadd, metType)
    outputFileName_hadd = os.path.join(outputFilePath, "UnclusteredEnergyAnalyzer_all_%s.root" % metType)
    command = '%s %s' % (executable_rm, outputFileName_hadd)
    print "--> command = '%s'" % command
    commandsToRun[metType].append(command)
    command = '%s %s %s' % (executable_hadd, outputFileName_hadd, make_vstring(outputFileNames[metType]))
    print "--> command = '%s'" % command
    commandsToRun[metType].append(command)
   
shellScriptName = os.path.join(configFilePath, "runUnclusteredEnergyAnalysis.sh")
shellScript = open(shellScriptName, "w")
shellScript.write("#!/bin/csh -f\n")
shellScript.write("\n")
for metType in commandsToRun.keys():
    for command in commandsToRun[metType]:
        shellScript.write("%s\n" % command)
shellScript.write("echo 'Finished running UnclusteredEnergyAnalysis.'\n")
shellScript.write("\n")
shellScript.close()

print("Finished building shell script. Now execute 'source %s'." % shellScriptName)
