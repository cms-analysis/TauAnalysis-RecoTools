import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.tools.buildConfigFilesTauIdEffAnalysis import make_inputFileNames_vstring, getStringRep_bool

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
    if len(inputFileNames_matched) == 0:
        raise ValueError("Sample = %s has no input files !!" % sampleNames)

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
    ##srcUnclPFCands = cms.InputTag('pfCandsNotInJet'),
    srcUnclPFCands = cms.InputTag('particleFlow'), # WARNING: logic to skip PFCandidates which are PFJet constituents not implemented yet !!

    srcTrigger = cms.InputTag('TriggerResults::HLT'),
    hltPaths = cms.vstring(%s),

    srcWeights = cms.VInputTag(%s),

    srcVertices = cms.InputTag('selectedPrimaryVertexPosition'),
    srcRhoNeutral = cms.InputTag('kt6PFNeutralJetsForVtxMultReweighting', 'rho'),
    srcRhoCharged = cms.InputTag('kt6PFChargedHadronJetsForVtxMultReweighting', 'rho')
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
                                        isCaloMEt,
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

    isCaloMEt = cms.bool(%s),

    refBranchName = cms.string('%s'),
    projParlBranchName = cms.string('%s'),
    projPerpBranchName = cms.string('%s'),

    outputFileName = cms.string('%s')
)
""" % (inputFileName, directory, 
       processType,
       getStringRep_bool(isCaloMEt),
       refBranchName, projParlBranchName, projPerpBranchName, outputFileName_full)

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
                                                      sampleName, runPeriod,
                                                      metOptionName, inputFilePath, outputFilePath, samplesToAnalyze, 
                                                      central_or_shift,
                                                      srcMEt, applyMEtShiftCorrection, srcMEtCov, sfMEtCov, srcJets,
                                                      srcNoPileUpMEtInputs, plotNoPileUpMEtInputs,
                                                      hltPaths, srcWeights,                                                      
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


    shiftedMEtCorrX_string = None
    shiftedMEtCorrY_string = None
    if metOptionName == "caloMEtNoHF":
        if runPeriod == "2012RunABCD": # CV: 2012RunD only
            if processType == 'Data':
                shiftedMEtCorrX_string = "-6.61856e-01 + 5.93230e-02*y" # CV: x = sumEt, y = numVertices
                shiftedMEtCorrY_string = "-2.26588e-01 - 1.24721e-02*y"
            else:
                shiftedMEtCorrX_string = "+1.52239e-01 - 1.14200e-01*y"
                shiftedMEtCorrY_string = "+8.57729e-01 - 5.68892e-02*y"
    elif metOptionName == "caloMEtTypeIcorrected":
        if runPeriod == "2012RunABCD": # CV: 2012RunD only
            if processType == 'Data':
                shiftedMEtCorrX_string = "-8.24198e-01 + 2.63584e-01*y" # CV: x = sumEt, y = numVertices
                shiftedMEtCorrY_string = "+7.00992e-02 - 1.34908e-01*y"
            else:
                shiftedMEtCorrX_string = "+3.76453e-01 - 2.83174e-01*y"
                shiftedMEtCorrY_string = "+1.23929e+00 - 3.26232e-01*y"            
    elif metOptionName == "pfMEt" or metOptionName == "pfMEtTypeIcorrected" or metOptionName == "pfMEtTypeIcorrectedSmeared":
        if runPeriod == "2012RunABCD":
            if processType == 'Data':
                shiftedMEtCorrX_string = "+4.83642e-02 + 2.48870e-01*y" # CV: x = sumEt, y = numVertices
                shiftedMEtCorrY_string = "-1.50135e-01 - 8.27917e-02*y"
            else:
                shiftedMEtCorrX_string = "+1.62861e-01 - 2.38517e-02*y"
                shiftedMEtCorrY_string = "+3.60860e-01 - 1.30335e-01*y"
    elif metOptionName == "pfMEtNoPileUp" or metOptionName == "pfMEtNoPileUpSmeared" or metOptionName == "pfchsMEtNoPileUpSmeared":
        if runPeriod == "2012RunABCD":
            if processType == 'Data':
                shiftedMEtCorrX_string = "+9.91346e-02 + 1.19915e-01*y"
                shiftedMEtCorrY_string = "-3.31467e-01 - 5.96704e-02*y"
            else:
                shiftedMEtCorrX_string = "+1.88526e-01 - 3.44137e-03*y"
                shiftedMEtCorrY_string = "-9.92166e-02 - 8.85794e-02*y"
    elif metOptionName == "pfMEtMVA" or metOptionName == "pfMEtMVASmeared":
        if runPeriod == "2012RunABCD":
            if processType == 'Data':
                shiftedMEtCorrX_string = "+9.20426e-02 + 5.24015e-02*y"
                shiftedMEtCorrY_string = "-1.07972e-01 - 5.67229e-03*y"
            else:
                shiftedMEtCorrX_string = "-2.35531e-01 - 2.16895e-02*y"
                shiftedMEtCorrY_string = "+6.28718e-02 - 2.74277e-02*y"
    elif metOptionName == "pfMEtMVAunityResponse" or metOptionName == "pfMEtMVAunityResponseSmeared":
        if runPeriod == "2012RunABCD":
            if processType == 'Data':
                shiftedMEtCorrX_string = "+7.45955e-02 + 7.61148e-02*y"
                shiftedMEtCorrY_string = "-7.40447e-02 - 2.41926e-02*y"
            else:
                shiftedMEtCorrX_string = "-1.82175e-01 - 4.26010e-02*y"
                shiftedMEtCorrY_string = "+2.10292e-01 - 6.09202e-02*y"
    elif metOptionName == "caloMEtTypeIcorrected":
        if runPeriod == "2012RunABCD":
            if processType == 'Data':
                shiftedMEtCorrX_string = "-7.62178e-01 + 2.30814e-01*y"
                shiftedMEtCorrY_string = "+3.36434e-01 - 1.22073e-01*y"
            else:
                shiftedMEtCorrX_string = "+3.12039e-01 - 2.75528e-01*y"
                shiftedMEtCorrY_string = "+6.39895e-01 - 2.38409e-01*y"
    shiftedMEtCorr_string = ""
    if shiftedMEtCorrX_string and shiftedMEtCorrY_string:
        shiftedMEtCorr_string  = "cms.PSet(\n"
        shiftedMEtCorr_string += "    numJetsMin = cms.int32(-1),\n"
        shiftedMEtCorr_string += "    numJetsMax = cms.int32(-1),\n"
        shiftedMEtCorr_string += "    px = cms.string('%s'),\n" % shiftedMEtCorrX_string
        shiftedMEtCorr_string += "    py = cms.string('%s')\n" % shiftedMEtCorrY_string
        shiftedMEtCorr_string += ")\n"        
    if applyMEtShiftCorrection and len(shiftedMEtCorr_string) == 0:
        raise ValueError("No sys. MET shifts defined for MET type = %s, run-period = %s !!" % (metOptionName, runPeriod))

    hltPaths_string = make_inputFileNames_vstring(hltPaths)
    srcWeights_string = make_inputFileNames_vstring(srcWeights)

    srcGenPileUpInfo = ""
    if processType == "MC_signal" or processType == "MC_background":
        srcGenPileUpInfo = "    srcGenPileUpInfo = cms.InputTag('addPileupInfo'),\n"
    
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

    plotNoPileUpMEtInputs_string = make_inputFileNames_vstring(plotNoPileUpMEtInputs)
    srcNoPileUpMEtSF = None
    if len(plotNoPileUpMEtInputs) >= 1:
        srcNoPileUpMEtSF = '%s:sfNoPU' % srcNoPileUpMEtInputs
    else:
        srcNoPileUpMEtSF = ''
    srcNoPileUpMEtSF = '' # CV: temporarily disabled, as PAT-tuples do not contain SF yet

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
    srcGenParticles = cms.InputTag('genParticles'),
    qT = cms.string("rec"),
    srcMuons = cms.InputTag('patMuons'), # CV: pat::Muon collection contains 'goodMuons' only
    srcMEt = cms.InputTag('%s'),
    srcMEtSignCovMatrix = cms.InputTag('%s'),
    sfMEtSignCovMatrix = cms.double(%f),
    srcJets = cms.InputTag('%s'),
    srcPFCandidates = cms.InputTag('particleFlow'),
    ##srcPFCandidates = cms.InputTag('pfCandsNotInJet'),

    shiftedMEtCorr = cms.PSet(
        jetPtThreshold = cms.double(30.0),
        parameter = cms.VPSet(
%s
        )
    ),
    applyMEtShiftCorr = cms.bool(%s),

    srcTrigger = cms.InputTag('TriggerResults::HLT'),
    hltPaths = cms.vstring(%s),

    srcWeights = cms.VInputTag(%s),

%s    

    srcVertices = cms.InputTag('selectedPrimaryVertexPosition'),
    srcRhoNeutral = cms.InputTag('kt6PFNeutralJetsForVtxMultReweighting', 'rho'),
    srcRhoCharged = cms.InputTag('kt6PFChargedHadronJetsForVtxMultReweighting', 'rho'),
%s

    srcNoPileUpMEtInputs = cms.InputTag('%s'),
    srcNoPileUpMEtSF = cms.InputTag('%s'),
    plotNoPileUpMEtInputs = cms.vstring(%s),

    selEventsFileName = cms.string('%s'),

    # CV: 'srcEventCounter' is defined in TauAnalysis/Skimming/test/skimTauIdEffSample_cfg.py
    srcEventCounter = cms.InputTag('processedEventsSkimming'),
    allEvents_DBS = cms.int32(%i),
    
    xSection = cms.double(%f),
    
    intLumiData = cms.double(%f)
)
""" % (maxEvents, fwliteInput_fileNames, outputFileName_full, directory,
       processType, recoZllRecoilCorrectionParameters_string,
       srcMEt, srcMEtCov, sfMEtCov, srcJets, shiftedMEtCorr_string, getStringRep_bool(applyMEtShiftCorrection),
       hltPaths_string, srcWeights_string, srcGenPileUpInfo, addPUreweight_string,
       srcNoPileUpMEtInputs, srcNoPileUpMEtSF, plotNoPileUpMEtInputs_string,
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
                                                      central_or_shift,
                                                      isCaloMEt,
                                                      plotNoPileUpMEtInputs):

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
    if runPeriod == "2012RunA":
        runPeriod_string = "Run A"
    elif runPeriod == "2012RunABC":
        runPeriod_string = "Run ABC"
    elif runPeriod == "2012RunABCD":
        runPeriod_string = "Run ABCD"

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

    plotNoPileUpMEtInputs_string = ""
    for plotNoPileUpMEtInput in plotNoPileUpMEtInputs:
        for histogramName_and_xAxisTitle in [ [ 'Pparl', 'P_{#parallel} / GeV' ],
                                              [ 'Pperp', 'P_{#perp} / GeV'     ],
                                              [ 'SumEt', '#Sigma E_{T} / GeV'  ] ] :
            plotNoPileUpMEtInputs_string += \
"""        
        cms.PSet(
            meName = cms.string('%s/%s'),
            xAxisTitle = cms.string('%s')
        ),
""" % (plotNoPileUpMEtInput, histogramName_and_xAxisTitle[0], histogramName_and_xAxisTitle[1])

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

    isCaloMEt = cms.bool(%s),

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
            xAxisTitle = cms.string('M_{#mu #mu} / GeV')
        ),
        cms.PSet(
            meName = cms.string('ZllCandPt'),
            xAxisTitle = cms.string('q_{T} / GeV')
        ),
        cms.PSet(
            meName = cms.string('metL'),
            xAxisTitle = cms.string('E_{T}^{miss} / GeV')
        ),
        cms.PSet(
            meName = cms.string('metPhi'),
            xAxisTitle = cms.string('#phi^{miss}')
        ),
        cms.PSet(
            meName = cms.string('metProjParlZ'),
            xAxisTitle = cms.string('-(u_{#parallel} + q_{T}) / GeV')
        ),
        cms.PSet(
            meName = cms.string('metProjPerpZ'),
            xAxisTitle = cms.string('u_{#perp}  / GeV')
        ),
        cms.PSet(
            meName = cms.string('metSigmaParlZ'),
            xAxisTitle = cms.string('#sigmaE_{#parallel}^{miss}')
        ),
        cms.PSet(
            meName = cms.string('metPullParlZ'),
            xAxisTitle = cms.string('E_{#parallel}^{miss} / #sigmaE_{#parallel}^{miss}')
        ),
        cms.PSet(
            meName = cms.string('metPullParlZNumJetsPtGt30Eq0'),
            xAxisTitle = cms.string('E_{#parallel}^{miss} / #sigmaE_{#parallel}^{miss}')
        ),
        cms.PSet(
            meName = cms.string('metPullParlZNumJetsPtGt30Eq1'),
            xAxisTitle = cms.string('E_{#parallel}^{miss} / #sigmaE_{#parallel}^{miss}')
        ),
        cms.PSet(
            meName = cms.string('metPullParlZNumJetsPtGt30Eq2'),
            xAxisTitle = cms.string('E_{#parallel}^{miss} / #sigmaE_{#parallel}^{miss}')
        ),
        cms.PSet(
            meName = cms.string('metPullParlZNumJetsPtGt30Ge3'),
            xAxisTitle = cms.string('E_{#parallel}^{miss} / #sigmaE_{#parallel}^{miss}')
        ),
        cms.PSet(
            meName = cms.string('metSigmaPerpZ'),
            xAxisTitle = cms.string('#sigmaE_{#perp  }^{miss}')
        ),        
        cms.PSet(
            meName = cms.string('metPullPerpZ'),
            xAxisTitle = cms.string('E_{#perp}^{miss}  / #sigmaE_{#perp}^{miss}')
        ),
        cms.PSet(
            meName = cms.string('metPullPerpZNumJetsPtGt30Eq0'),
            xAxisTitle = cms.string('E_{#perp}^{miss}  / #sigmaE_{#perp}^{miss}')
        ),
        cms.PSet(
            meName = cms.string('metPullPerpZNumJetsPtGt30Eq1'),
            xAxisTitle = cms.string('E_{#perp}^{miss}  / #sigmaE_{#perp}^{miss}')
        ),
        cms.PSet(
            meName = cms.string('metPullPerpZNumJetsPtGt30Eq2'),
            xAxisTitle = cms.string('E_{#perp}^{miss}  / #sigmaE_{#perp}^{miss}')
        ),
        cms.PSet(
            meName = cms.string('metPullPerpZNumJetsPtGt30Ge3'),
            xAxisTitle = cms.string('E_{#perp}^{miss}  / #sigmaE_{#perp}^{miss}')
        ),        
        cms.PSet(
            meName = cms.string('uParl'),
            xAxisTitle = cms.string('u_{#parallel} / GeV')
        ),
        cms.PSet(
            meName = cms.string('uPerp'),
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
            meName = cms.string('numJetsRawPtGt30'),
            xAxisTitle = cms.string('Num. Jets, P_{T}^{raw} > 30 GeV')
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
            meName = cms.string('numJetsCorrPtGt30'),
            xAxisTitle = cms.string('Num. Jets, P_{T}^{corr} > 30 GeV')
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
        ),
%s        
    ),

    outputFileName = cms.string('%s')
)
""" % (inputFileName, runPeriod_string, directoryData, directoryMC_signal, make_inputFileNames_vstring(directoryMCs_background),
       mcScaleFactors_string,
       sysShiftsUp, sysShiftsDown,
       getStringRep_bool(isCaloMEt),
       plotNoPileUpMEtInputs_string,
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


