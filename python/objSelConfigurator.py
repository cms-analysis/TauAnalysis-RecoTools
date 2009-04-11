import FWCore.ParameterSet.Config as cms
import sys

class objSelConfigurator(cms._ParameterTypeBase):

    def __init__(self, objSelList, src = None, pyModuleName = None, doSelCumulative = True, doSelIndividual = False):
        self.objSelList = objSelList
        self.src = src
        self.pyModuleName = pyModuleName,
        self.doSelCumulative = doSelCumulative
        self.doSelIndividual = doSelIndividual

    @staticmethod
    def _composeModuleName(part_1, part_2):
        # auxiliary function for concatenating two strings;
        # if the last character of part_1 is lower-case (upper-case),
        # capitalize (lowercase) the first character of part_2
        if part_1[-1].islower():
            return part_1 + part_2.capitalize()
        else:
            return part_1 + part_2[0].lower() + part_2[1:]

    class _getterCumulative:
        # auxiliary class for composing name of module selecting "cumulative" collection
        @staticmethod
        def get_src(src, lastModuleName):
            if lastModuleName == None:
                return src
            else:
                return lastModuleName
        @staticmethod
        def get_moduleName(name):
            return objSelConfigurator._composeModuleName(name, "Cumulative")

    class _getterIndividual:
        # auxiliary class for composing name of module selecting "individual" collection
        @staticmethod
        def get_src(src, lastModuleName):
            return src
        @staticmethod
        def get_moduleName(name):
            return objSelConfigurator._composeModuleName(name, "Individual")

    def _addModule(self, objSelItem, getter):
        # create module
        moduleType = objSelItem.type.value()
        module = cms.EDFilter(moduleType)

        # set module attributes
        for objSelAttrName in dir(objSelItem):
            objSelAttr = getattr(objSelItem, objSelAttrName)
            if isinstance(objSelAttr, cms._ParameterTypeBase) and not objSelAttrName in ["name", "type"]:
                setattr(module, objSelAttrName, objSelAttr)
        src = getter.get_src(self.src, self.lastModuleName)
        setattr(module, "src", cms.InputTag(src))

        moduleName = getter.get_moduleName(objSelItem.name.value())
        
        # register module in global python name-space
        pyModule = sys.modules[self.pyModuleName[0]]
        if pyModule is None:
            raise ValueError("pyModuleName Parameter invalid !!")
        setattr(pyModule, moduleName, module)
        
        self.lastModuleName = moduleName

        # add module to sequence
        if self.sequence == None:
            self.sequence = module
        else:
            self.sequence *= module

    def configure(self):
        # configure modules for "cumulative" and "individual" collections
        # of objects passing selection

        if self.src is None:
            raise ValueError("src Parameter must not be empty !!")

        self.sequence = None
        self.lastModuleName = None

        if self.doSelCumulative:
            getter = objSelConfigurator._getterCumulative()
            for objSelItem in self.objSelList:
                self._addModule(objSelItem, getter)

        if self.doSelIndividual:
            getter = objSelConfigurator._getterIndividual()
            for objSelItem in self.objSelList:
                self._addModule(objSelItem, getter)

        return cms.Sequence(self.sequence)
