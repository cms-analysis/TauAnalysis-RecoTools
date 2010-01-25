import FWCore.ParameterSet.Config as cms
import copy
import sys

from TauAnalysis.CandidateTools.tools.getInstanceName import getInstanceName
from TauAnalysis.CandidateTools.tools.composeModuleName import composeModuleName

#--------------------------------------------------------------------------------
# utility function for configuring BoolEventSelFlagProducer modules
#
# Author: Christian Veelken, UC Davis
#
#--------------------------------------------------------------------------------

class eventSelFlagProdConfigurator(cms._ParameterTypeBase):

    def __init__(self, objSelList, boolEventSelFlagProducer = "BoolEventSelFlagProducer", pyModuleName = None):
        self.objSelList = objSelList
        self.boolEventSelFlagProducer = boolEventSelFlagProducer
        self.pyModuleName = pyModuleName

    def _addModule(self, objSelItem, srcParam, sysName = None):
        # create module
        module = cms.EDFilter(self.boolEventSelFlagProducer)

        # check if multiple plugins need to be created for
        # "cumulative" and "individual" collections
        if len(srcParam) >= 2:
            # multiple plugin case;
            # prepare configuration parameters for plugins
            # derrived from EventSelectorBase class
            exclAttributeNames = [ "src_cumulative", "src_individual" ]
            selector = cms.PSet()
            for objSelAttrName in dir(objSelItem):
                objSelAttr = getattr(objSelItem, objSelAttrName)
                if isinstance(objSelAttr, cms._ParameterTypeBase) and not objSelAttrName in exclAttributeNames:
                    setattr(selector, objSelAttrName, objSelAttr)
                
            selectors = []        
            for src_i in srcParam:
                assert len(src_i) == 2, "'src_i' Parameter must have exactly 2 elements !!"
                selector_i = copy.deepcopy(selector)
                if sysName is None:
                    setattr(selector_i, "src", src_i[0])
                else:
                    setattr(selector_i, "src", composeModuleName(src_i[0], sysName))
                if len(src_i) >= 2 : setattr(selector_i, "instanceName", cms.string(src_i[1]))
                selectors.append(selector_i)

            setattr(module, "selectors", cms.VPSet(*selectors))
        else:
            # single plugin case;
            # add configuration parameters directly to module
            exclAttributeNames = [ "src", ]
            for objSelAttrName in dir(objSelItem):
                objSelAttr = getattr(objSelItem, objSelAttrName)
                if isinstance(objSelAttr, cms._ParameterTypeBase) and not objSelAttrName in exclAttributeNames:
                    setattr(module, objSelAttrName, objSelAttr)

            src = srcParam[0]
            if sysName is None:
                setattr(module, "src", src[0])
            else:
                setattr(module, "src", cms.InputTag(composeModuleName(src[0].value(), sysName)))
            if len(src) >= 2 : setattr(module, "instanceName", cms.string(src[1]))

        # set module attributes
        moduleName = None
        if sysName is None:
            moduleName = objSelItem.pluginName.value()
        else:
            moduleName = composeModuleName(objSelItem.pluginName.value(), sysName)
        module.setLabel(moduleName)
               
        # register module in global python name-space
        pyModule = sys.modules[self.pyModuleName]
        if pyModule is None:
            raise ValueError("'pyModuleName' Parameter invalid !!")
        setattr(pyModule, moduleName, module)

        # add module to sequence
        if self.sequence == None:
            self.sequence = module
        else:
            self.sequence *= module

    def configure(self):
        # configure 'BoolEventSelFlagProducer' modules
        # storing results of event selection as boolean flags in the event

        if self.boolEventSelFlagProducer is None:
            raise ValueError("'boolEventSelFlagProducer' Parameter must not be empty !!")

        self.sequence = None

        for objSelItem in self.objSelList:
            # check if configuration parameter set includes
            # either "src" or "src_cumulative" and "src_individual" parameters;
            # in the latter case, add two plugins with different instanceName parameters to the module
            src = None
            src_cumulative = None
            src_individual = None
            systematics = None
            for objSelAttrName in dir(objSelItem):
                if objSelAttrName == "src":
                    src = getattr(objSelItem, objSelAttrName)
                if objSelAttrName == "src_cumulative":
                    src_cumulative = getattr(objSelItem, objSelAttrName)
                if objSelAttrName == "src_individual":
                    src_individual = getattr(objSelItem, objSelAttrName)
                if objSelAttrName == "systematics":
                    systematics = getattr(objSelItem, objSelAttrName)

            if ( src is not None and src_cumulative is not None and src_individual is not None ) or \
               ( src is None and (src_cumulative is None or src_individual is None) ):
                raise ValueError("must specify either 'src' or 'src_cumulative' and 'src_individual' Parameters !!")

            if src is not None:
                self._addModule(objSelItem, [ [ src, ] ] )
                if systematics is not None:
                    for sysName in systematics:
                        self._addModule(objSelItem, [ src, ], sysName )
            if src_cumulative is not None and src_individual is not None:
                self._addModule(objSelItem, [ [ src_cumulative, "cumulative" ], [ src_individual, "individual" ] ])
                if systematics is not None:
                    for sysName in systematics:
                        self._addModule(objSelItem, [ [ src_cumulative, "cumulative" ] ], sysName)

        return cms.Sequence(self.sequence)