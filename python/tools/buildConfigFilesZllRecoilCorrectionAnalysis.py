import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.tools.buildConfigFilesTauIdEffAnalysis import make_inputFileNames_vstring

import os
import re

#--------------------------------------------------------------------------------
#
# define auxiliary functions
#
def getPATtupleFileNames(sampleNames, inputFilePath):
    
    inputFileNames = os.listdir(inputFilePath)
    #print(inputFileNames)

    # check if inputFile is PAT-tuple and
    # matches sample to be analyzed
    inputFileNames_matched = []
    fwliteInput_fileNames = ""
    for inputFileName in inputFileNames:
        isMatched = False
        for sampleName in sampleNames:
            inputFile_regex = \
              r"ZllRecoilCorrectionPATtuple_%s_[a-zA-Z0-9_]*.root" % sampleName
            inputFile_matcher = re.compile(inputFile_regex)
            if inputFile_matcher.match(inputFileName):
                isMatched = True
        if isMatched:
            inputFileNames_matched.append(os.path.join(inputFilePath, inputFileName))
            fwliteInput_fileNames += "process.fwliteInput.fileNames.append('%s')\n" % os.path.join(inputFilePath, inputFileName)

    print " found %i input files." % len(inputFileNames_matched)

    return (inputFileNames_matched, fwliteInput_fileNames) 
#--------------------------------------------------------------------------------

def buildConfigFile_produceZllRecoilNtuples(maxEvents,
                                            sampleName, metOptionName, inputFilePath, outputFilePath, samplesToAnalyze,
                                            central_or_shift, srcMEt, srcJets, hltPaths, srcWeights):

    """Build cfg.py file to run FWLiteZllRecoilCorrectionNtupleProducer macro on PAT-tuples,
       and produce 'plain' ROOT Ntuple needed for fitting Z-recoil correction parameters"""

    print "<buildConfigFile_produceZllRecoilNtuples>:"
    print " processing sample %s" % sampleName

    inputFileNames = getPATtupleFileNames(samplesToAnalyze[sampleName]['samples'], inputFilePath)
    if len(inputFileNames[0]) == 0:
        print("Sample %s has no input files --> skipping !!" % sampleName)
        return
    fwliteInput_fileNames = inputFileNames[1]

    print(" building config file...")

    outputFileName = 'ZllRecoilCorrectionNtuple_%s_%s_%s.root' % (sampleName, metOptionName, central_or_shift)
    outputFileName_full = os.path.join(outputFilePath, outputFileName)
    
    directory = sampleName

    hltPaths_string = make_inputFileNames_vstring(hltPaths)
    srcWeights_string = make_inputFileNames_vstring(srcWeights)

    addPUreweight_string = ""
# CV: do not apply rho_neutral reweighting to Monte Carlo samples other than Zmumu 
#     until difference in rho_neutral distribution between Zmummu, 
#     TTbar and WW/WZ/ZZ is understood and corrected for
    if samplesToAnalyze[sampleName]['isMC'] and samplesToAnalyze[sampleName]['applyRhoNeutralReweighting']:
        addPUreweight_string = \
"""
    addPUreweight = cms.PSet(
        inputFileName = cms.FileInPath('TauAnalysis/RecoTools/data/vertexMultiplicityVsRhoPFNeutralReweight.root'),
        meName = cms.string('histoReweight_fitted'),
        minPUreweight = cms.double(1.e-1),
        maxPUreweight = cms.double(1.e+1)
    )
"""

    config = \
"""
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring(),
    
    maxEvents = cms.int32(%i),
    
    outputEvery = cms.uint32(1000)
)

%s
    
process.fwliteOutput = cms.PSet(
    fileName = cms.string('%s')
)

process.ZllRecoilCorrectionNtupleProducer = cms.PSet(

    directory = cms.string('%s'),

    srcZllCandidates = cms.InputTag('goldenZmumuCandidatesGe1IsoMuons'),
    srcMEt = cms.InputTag('%s'),
    srcJets = cms.InputTag('%s'),
    srcUnclPFCands = cms.InputTag('pfCandsNotInJet'),

    srcTrigger = cms.InputTag('TriggerResults::HLT'),
    hltPaths = cms.vstring(%s),

    srcWeights = cms.VInputTag(%s),

    srcVertices = cms.InputTag('selectedPrimaryVertexPosition'),
    srcRhoNeutral = cms.InputTag('kt6PFNeutralJetsForVtxMultReweighting', 'rho'),
%s
)
""" % (maxEvents, fwliteInput_fileNames, outputFileName_full, directory, 
       srcMEt, srcJets, hltPaths_string, srcWeights_string, addPUreweight_string)

    configFileName = "produceZllRecoilCorrectionNtuple_%s_%s_%s_cfg.py" % (sampleName, metOptionName, central_or_shift)
    configFileName_full = os.path.join(outputFilePath, configFileName)    
    configFile = open(configFileName_full, "w")
    configFile.write(config)
    configFile.close()

    logFileName = configFileName.replace('_cfg.py', '.log')
    logFileName_full = os.path.join(outputFilePath, logFileName)

    retVal = {}
    retVal['configFileName'] = configFileName_full
    retVal['outputFileName'] = outputFileName_full
    retVal['logFileName']    = logFileName_full

    return retVal

def buildConfigFile_fitZllRecoilNtuples(sampleName, metOptionName, inputFileName, outputFilePath, samplesToAnalyze, central_or_shift, 
                                        refBranchName = "qT", projParlBranchName = "uParl", projPerpBranchName = "uPerp"):

    """Build cfg.py file to run fitZllRecoilCorrection macro to run on 'plain' ROOT Ntuples
       and fit Z-recoil correction parameters"""

    print "<buildConfigFile_fitZllRecoilNtuples>:"
    print " processing sample %s" % sampleName


    print(" building config file...")

    directory = sampleName
    
    processType = None
    if samplesToAnalyze[sampleName]['isMC']:
        processType = 'MC'
    else:
        processType = 'Data'

    outputFileName = "fittedZllRecoilCorrectionParameters_%s_%s_%s_%s_cfi.py" % (sampleName, metOptionName, refBranchName, central_or_shift)
    outputFileName_full = os.path.join(outputFilePath, outputFileName)

    config = \
"""
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring('%s'),
    
    maxEvents = cms.int32(-1),
    
    outputEvery = cms.uint32(1000)
)

process.fitZllRecoilCorrection = cms.PSet(

    directory = cms.string('%s'),

    type = cms.string('%s'),

    refBranchName = cms.string('%s'),
    projParlBranchName = cms.string('%s'),
    projPerpBranchName = cms.string('%s'),

    outputFileName = cms.string('%s')
)
""" % (inputFileName, directory, 
       processType, refBranchName, projParlBranchName, projPerpBranchName, outputFileName_full)

    configFileName = "fitZllRecoilCorrectionNtuple_%s_%s_%s_%s_cfg.py" % (sampleName, metOptionName, refBranchName, central_or_shift)
    configFileName_full = os.path.join(outputFilePath, configFileName)    
    configFile = open(configFileName_full, "w")
    configFile.write(config)
    configFile.close()

    logFileName = configFileName.replace('_cfg.py', '.log')
    logFileName_full = os.path.join(outputFilePath, logFileName)

    retVal = {}
    retVal['configFileName'] = configFileName_full
    retVal['outputFileName'] = outputFileName_full
    retVal['logFileName']    = logFileName_full

    return retVal

def buildConfigFile_FWLiteZllRecoilCorrectionAnalyzer(maxEvents,
                                                      sampleName, metOptionName, inputFilePath, outputFilePath, samplesToAnalyze, 
                                                      central_or_shift, srcMEt, srcJets, hltPaths, srcWeights, 
                                                      ZllRecoilCorrectionParameterFileNames, intLumiData):

    """Build cfg.py file to run FWLiteZllRecoilCorrectionAnalyzer macro on PAT-tuples,
       and fill control plots of MET in Data compared to Monte Carlo simulation with Z-recoil corrections applied"""

    print "<buildConfigFile_FWLiteZllRecoilCorrectionAnalyzer>:"
    print " processing sample %s" % sampleName

    inputFileNames = getPATtupleFileNames(samplesToAnalyze[sampleName]['samples'], inputFilePath)
    if len(inputFileNames[0]) == 0:
        print("Sample %s has no input files --> skipping !!" % sampleName)
        return
    fwliteInput_fileNames = inputFileNames[1]

    print(" building config file...")

    outputFileName = 'analyzeZllRecoilCorrectionHistograms_%s_%s_%s.root' % (sampleName, metOptionName, central_or_shift)
    outputFileName_full = os.path.join(outputFilePath, outputFileName)
    
    directory = sampleName
    if central_or_shift != 'central':
        directory += '/%s' % central_or_shift

    processType = None
    if samplesToAnalyze[sampleName]['isMC']:
        if samplesToAnalyze[sampleName]['Type'] == 'Signal':
            processType = 'MC_signal'
        elif samplesToAnalyze[sampleName]['Type'] == 'Background':
            processType = 'MC_background'
        else:
            raise ValueError("Invalid MC type = %s !!" % samplesToAnalyze[sampleName]['Type'])
    else:
        processType = 'Data'

    recoZllRecoilCorrectionParameters_string = ""
    if ZllRecoilCorrectionParameterFileNames is not None:
        recoZllRecoilCorrectionParameters_string = "    algorithm = cms.PSet(\n"
        recoZllRecoilCorrectionParameters_string += "        parameter = cms.PSet(\n"
        for parameterSetName in [ 'data', 'mc' ]: 
            ZllRecoilCorrectionParameterFile = open(ZllRecoilCorrectionParameterFileNames[parameterSetName], 'r')
            lines = ZllRecoilCorrectionParameterFile.readlines()
            for lineNumber, line in enumerate(lines):
                #print "%i: %s" % (lineNumber, line)
                if lineNumber == 0:
                    recoZllRecoilCorrectionParameters_string += "            %s = cms.PSet(\n" % parameterSetName
                elif lineNumber < (len(lines) - 1):
                    recoZllRecoilCorrectionParameters_string += "            %s" % line
            recoZllRecoilCorrectionParameters_string += "            ),\n"
            ZllRecoilCorrectionParameterFile.close()
        recoZllRecoilCorrectionParameters_string += "        )\n"
        recoZllRecoilCorrectionParameters_string += "    ),\n"

    hltPaths_string = make_inputFileNames_vstring(hltPaths)
    srcWeights_string = make_inputFileNames_vstring(srcWeights)
    
    allEvents_DBS = 0
    xSection = 0.
    if samplesToAnalyze[sampleName]['isMC']:
        allEvents_DBS = samplesToAnalyze[sampleName]['allEvents_DBS']
        xSection = samplesToAnalyze[sampleName]['xSection']

    addPUreweight_string = ""
# CV: do not apply rho_neutral reweighting to Monte Carlo samples other than Zmumu
#     until difference in rho_neutral distribution between Zmummu,
#     TTbar and WW/WZ/ZZ is understood and corrected for
    if samplesToAnalyze[sampleName]['isMC'] and samplesToAnalyze[sampleName]['applyRhoNeutralReweighting']:
        addPUreweight_string = \
"""
    addPUreweight = cms.PSet(
        inputFileName = cms.FileInPath('TauAnalysis/RecoTools/data/vertexMultiplicityVsRhoPFNeutralReweight.root'),
        meName = cms.string('histoReweight_fitted'),
        minPUreweight = cms.double(1.e-1),
        maxPUreweight = cms.double(1.e+1)
    ),
"""

    selEventsFileName = None 
    if central_or_shift == 'central':
        selEventsFileName = 'selEvents_%s_%s.txt' % (sampleName, metOptionName)
    else:
        selEventsFileName = ''
        
    config = \
"""
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring(),
    
    maxEvents = cms.int32(%i),
    
    outputEvery = cms.uint32(1000)
)

%s
    
process.fwliteOutput = cms.PSet(
    fileName = cms.string('%s')
)

process.ZllRecoilCorrectionAnalyzer = cms.PSet(

    directory = cms.string('%s'),

    type = cms.string('%s'),

%s

    srcZllCandidates = cms.InputTag('goldenZmumuCandidatesGe1IsoMuons'),
    srcMuons = cms.InputTag('patMuons'), # CV: pat::Muon collection contains 'goodMuons' only
    srcMEt = cms.InputTag('%s'),
    srcJets = cms.InputTag('%s'),

    srcTrigger = cms.InputTag('TriggerResults::HLT'),
    hltPaths = cms.vstring(%s),

    srcWeights = cms.VInputTag(%s),

    srcVertices = cms.InputTag('selectedPrimaryVertexPosition'),
    srcRhoNeutral = cms.InputTag('kt6PFNeutralJetsForVtxMultReweighting', 'rho'),
%s

    selEventsFileName = cms.string('%s'),

    # CV: 'srcEventCounter' is defined in TauAnalysis/Skimming/test/skimTauIdEffSample_cfg.py
    srcEventCounter = cms.InputTag('processedEventsSkimming'),
    allEvents_DBS = cms.int32(%i),
    
    xSection = cms.double(%f),
    
    intLumiData = cms.double(%f)
)
""" % (maxEvents, fwliteInput_fileNames, outputFileName_full, directory,
       processType, recoZllRecoilCorrectionParameters_string,
       srcMEt, srcJets, hltPaths_string, srcWeights_string, addPUreweight_string,
       os.path.join(outputFilePath, selEventsFileName), allEvents_DBS, xSection, intLumiData)

    configFileName = "analyzeZllRecoilCorrectionPATtuple_%s_%s_%s_cfg.py" % (sampleName, metOptionName, central_or_shift)
    configFileName_full = os.path.join(outputFilePath, configFileName)    
    configFile = open(configFileName_full, "w")
    configFile.write(config)
    configFile.close()

    logFileName = configFileName.replace('_cfg.py', '.log')
    logFileName_full = os.path.join(outputFilePath, logFileName)

    retVal = {}
    retVal['configFileName'] = configFileName_full
    retVal['outputFileName'] = outputFileName_full
    retVal['logFileName']    = logFileName_full

    return retVal

def buildConfigFile_makeZllRecoilCorrectionFinalPlots(sampleNameData, sampleNameMC_signal, sampleNameMCs_background, runPeriod, 
                                                      metOptionName, inputFileName, outputFilePath, samplesToAnalyze, corrLevelMC,
                                                      central_or_shift):

    """Build cfg.py file to run makeZllRecoilCorrectionFinalPlots macro
       and make final control plots of MET in Data compared to Monte Carlo simulation with Z-recoil corrections applied"""

    print "<buildConfigFile_makeZllRecoilCorrectionFinalPlots>:"
    print " processing combination of sampleData %s, sampleMC (signal) %s" % (sampleNameData, sampleNameMC_signal)

    print(" building config file...")

    #corrLevelData = "beforeGenPUreweight"
    corrLevelData = corrLevelMC
    
    directoryData           =   "/".join([ sampleNameData,      corrLevelData ])
    directoryMC_signal      =   "/".join([ sampleNameMC_signal, corrLevelMC   ])
    directoryMCs_background = [ "/".join([ sampleNameMC_bgr,    corrLevelData ]) for sampleNameMC_bgr in sampleNameMCs_background ]

    sampleNameMCs = []
    sampleNameMCs.append(sampleNameMC_signal)
    sampleNameMCs.extend(sampleNameMCs_background)
    mcScaleFactors_string = ""
    mcScaleFactors_string += "    mcScaleFactors = cms.PSet(\n"
    for sampleNameMC in sampleNameMCs:
        if 'scaleFactor' in samplesToAnalyze[sampleNameMC].keys():
            mcScaleFactors_string += "        %s = cms.double(%f),\n" % (sampleNameMC, samplesToAnalyze[sampleNameMC]['scaleFactor'])
    mcScaleFactors_string += "    ),\n"

    runPeriod_string = ""
    if runPeriod == "2011RunA":
        runPeriod_string = "Run A"
    elif runPeriod == "2011RunB":
        runPeriod_string = "Run B"

    sysShiftsUp   = []
    sysShiftsDown = []
    for sysShift in central_or_shift:
        print "sysShift = %s" % sysShift
        if sysShift.find('Up') != -1:
            sysShiftsUp.append(sysShift)
            sysShiftsDown.append(sysShift.replace("Up", "Down"))

    outputFileName = "plotZllRecoilCorrection_%s_%s.png" % (metOptionName, corrLevelMC)
    outputFilePath_plots = os.path.join(outputFilePath, "plots")
    if not os.path.exists(outputFilePath_plots):
        os.mkdir(outputFilePath_plots)    
    outputFileName_full = os.path.join(outputFilePath_plots, outputFileName)

    config = \
"""
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring('%s'),
    
    maxEvents = cms.int32(-1),
    
    outputEvery = cms.uint32(1000)
)

process.makeZllRecoilCorrectionFinalPlots = cms.PSet(

    runPeriod = cms.string('%s'),

    directoryData           = cms.string('%s'),
    directoryMC_signal      = cms.string('%s'),
    directoryMCs_background = cms.vstring(%s),

%s

    sysShiftsUp   = cms.vstring(%s),
    sysShiftsDown = cms.vstring(%s),

    variables = cms.VPSet(
        cms.PSet(
            meName = cms.string('lPlusPt'),
            xAxisTitle = cms.string('P_{T}^{#mu+} / GeV')
        ),
        cms.PSet(
            meName = cms.string('lMinusPt'),
            xAxisTitle = cms.string('P_{T}^{#mu-} / GeV')
        ),
        cms.PSet(
            meName = cms.string('ZllCandMass'),
            xAxisTitle = cms.string('M_{#mu #mu}GeV')
        ),
        cms.PSet(
            meName = cms.string('metS'),
            xAxisTitle = cms.string('E_{T}^{miss} / GeV')
        ),
        cms.PSet(
            meName = cms.string('metL'),
            xAxisTitle = cms.string('E_{T}^{miss} / GeV')
        ),
        cms.PSet(
            meName = cms.string('metProjParlZ'),
            xAxisTitle = cms.string('u_{#parallel} / GeV')
        ),
        cms.PSet(
            meName = cms.string('metProjPerpZ'),
            xAxisTitle = cms.string('u_{#perp}  / GeV')
        ),
        cms.PSet(
            meName = cms.string('numJetsRawPtGt10'),
            xAxisTitle = cms.string('Num. Jets, P_{T}^{raw} > 10 GeV')
        ),
        cms.PSet(
            meName = cms.string('numJetsRawPtGt15'),
            xAxisTitle = cms.string('Num. Jets, P_{T}^{raw} > 15 GeV')
        ),
        cms.PSet(
            meName = cms.string('numJetsRawPtGt20'),
            xAxisTitle = cms.string('Num. Jets, P_{T}^{raw} > 20 GeV')
        ),
        cms.PSet(
            meName = cms.string('numJetsCorrPtGt10'),
            xAxisTitle = cms.string('Num. Jets, P_{T}^{corr} > 10 GeV')
        ),
        cms.PSet(
            meName = cms.string('numJetsCorrPtGt15'),
            xAxisTitle = cms.string('Num. Jets, P_{T}^{corr} > 15 GeV')
        ),
        cms.PSet(
            meName = cms.string('numJetsCorrPtGt20'),
            xAxisTitle = cms.string('Num. Jets, P_{T}^{corr} > 20 GeV')
        ),
        cms.PSet(
            meName = cms.string('numVertices'),
            xAxisTitle = cms.string('rec. Vertex Multiplicity')
        ),
        cms.PSet(
            meName = cms.string('vertexZ'),
            xAxisTitle = cms.string('z_{vtx} / cm')
        ),
        cms.PSet(
            meName = cms.string('rhoNeutral'),
            xAxisTitle = cms.string('#rho_{h0} / GeV')
        )
    ),

    outputFileName = cms.string('%s')
)
""" % (inputFileName, runPeriod_string, directoryData, directoryMC_signal, make_inputFileNames_vstring(directoryMCs_background),
       mcScaleFactors_string,
       sysShiftsUp, sysShiftsDown,
       outputFileName_full)

    configFileName = "makeZllRecoilCorrectionFinalPlots_%s_%s_cfg.py" % (metOptionName, corrLevelMC)
    configFileName_full = os.path.join(outputFilePath, configFileName)    
    configFile = open(configFileName_full, "w")
    configFile.write(config)
    configFile.close()

    logFileName = configFileName.replace('_cfg.py', '.log')
    logFileName_full = os.path.join(outputFilePath, logFileName)

    retVal = {}
    retVal['configFileName'] = configFileName_full
    retVal['outputFileName'] = outputFileName_full
    retVal['logFileName']    = logFileName_full

    return retVal


