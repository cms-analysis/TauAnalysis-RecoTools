import FWCore.ParameterSet.Config as cms

from RecoEcal.EgammaClusterProducers.geometryForClustering_cff import *
from RecoEgamma.ElectronIdentification.likelihoodPdfsDB_cfi import *
from RecoEgamma.ElectronIdentification.likelihoodESetup_cfi import *
from RecoEgamma.ElectronIdentification.neuralNetElectronId_cfi import *

from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi        import eidCutBasedExt
from RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi import eidCutBasedClassesExt as eidPtdrExt ## PTDR is how people know it
from RecoEgamma.ElectronIdentification.electronIdNeuralNetExt_cfi       import eidNeuralNetExt
from RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi      import eidLikelihoodExt 
 
elecIdForTauAnalysesCutBasedExtPAT = eidCutBasedExt.copy();  
elecIdForTauAnalysesCutBasedExtPAT.src = cms.InputTag("allLayer0ElectronsForTauAnalyses")

elecIdForTauAnalysesCutBasedRobust = elecIdForTauAnalysesCutBasedExtPAT.copy();
elecIdForTauAnalysesCutBasedRobust.electronQuality = 'robust'

elecIdForTauAnalysesCutBasedLoose  = elecIdForTauAnalysesCutBasedExtPAT.copy();
elecIdForTauAnalysesCutBasedLoose.electronQuality  = 'loose'

elecIdForTauAnalysesCutBasedTight  = elecIdForTauAnalysesCutBasedExtPAT.copy();
elecIdForTauAnalysesCutBasedTight.electronQuality  = 'tight'

patLayer0ElectronId = cms.Sequence( elecIdForTauAnalysesCutBasedRobust
                                   *elecIdForTauAnalysesCutBasedLoose
                                   *elecIdForTauAnalysesCutBasedTight )
