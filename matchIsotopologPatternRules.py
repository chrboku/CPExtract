from formulaTools import formulaTools
from utils import Bunch

import abc

_formulaTools=formulaTools()

class Rule:
    def getAllRequiredIsotopologs(self):
        return []

    def check(self, isotopologDict, log=False):
        return False

    def getChromatographicPeaks(self):
        return []

    @abc.abstractmethod
    def getMessage(self):
        pass




class PresenceRule(Rule):
    def __init__(self, otherIsotopolog="[13C]2", minIntensity=10000, mustBePresent=True, verifyChromPeakSimilarity=True, ratioWindows={'X': [1/2., 2]}):
        self.otherIsotopolog=otherIsotopolog
        self.minIntensity=minIntensity
        self.mustBePresent=mustBePresent
        self.verifyChromPeakSimilarity=verifyChromPeakSimilarity
        self.ratioWindows=ratioWindows


    def getAllRequiredIsotopologs(self):
        return [self.otherIsotopolog] + self.ratioWindows.keys()

    def check(self, isotopologDict, log=False):
        if self.otherIsotopolog not in isotopologDict.keys():
            if self.mustBePresent:
                return False
            if not self.mustBePresent:
                return True

        if self.mustBePresent and isotopologDict[self.otherIsotopolog].intensity<self.minIntensity:
            return False

        ratWinPassed=True
        for isotopolog, window in self.ratioWindows.items():
            if isotopolog in isotopologDict.keys():
                ratio=isotopologDict[self.otherIsotopolog].intensity/isotopologDict[isotopolog].intensity

                if window[0] <= ratio <= window[1]:
                    pass
                else:
                    ratWinPassed=False

        return ratWinPassed

    def getChromatographicPeaks(self):
        if self.verifyChromPeakSimilarity:
            return [self.otherIsotopolog]
        else:
            return []

    def getMessage(self):
        return "Isotopolog %s is present and its ratios check out as well"%(self.otherIsotopolog)




class AbsenceRule(Rule):
    def __init__(self, otherIsotopolog="[13C]-1", maxRatio=0.1):
        self.otherIsotopolog=otherIsotopolog
        self.maxRatio=maxRatio

    def getAllRequiredIsotopologs(self):
        return [self.otherIsotopolog]

    def check(self, isotopologDict, log=False):
        if self.otherIsotopolog not in isotopologDict.keys():
            return True

        if self.otherIsotopolog in isotopologDict.keys():
            if (isotopologDict[self.otherIsotopolog].intensity/isotopologDict["X"].intensity)<self.maxRatio:
                return True
            else:
                return False

    def getMessage(self):
        return "Isotopolog %s is not present or below the limit (%.2f)"%(self.otherIsotopolog, self.maxRatio)




class AnyIntensityRule(Rule):
    def __init__(self, anyIsotopolog=["X"], minimumIntensity=1E5):
        self.anyIsotopolog=anyIsotopolog
        self.minimumIntensity=minimumIntensity

    def getAllRequiredIsotopologs(self):
        return self.anyIsotopolog

    def check(self, isotopologDict, log=False):
        for iso in self.anyIsotopolog:
            if iso in isotopologDict.keys():
                if isotopologDict[iso].intensity>=self.minimumIntensity:
                    return True

        return False

    def getMessage(self):
        return "No isotopolog %s is above the set minimum intensity threshold (%.2f)"%(self.anyIsotopolog, self.minimumIntensity)




class AllIntensityRule(Rule):
    def __init__(self, allIsotopolog=["X"], minimumIntensity=1E5):
        self.allIsotopolog=allIsotopolog
        self.minimumIntensity=minimumIntensity

    def getAllRequiredIsotopologs(self):
        return self.allIsotopolog

    def check(self, isotopologDict, log=False):
        for iso in self.allIsotopolog:
            if iso in isotopologDict.keys():
                if isotopologDict[iso].intensity<self.minimumIntensity:
                    return False
            else:
                return False

        return True

    def getMessage(self):
        return "All isotopologs %s are below the set minimum intensity threshold (%.2f)"%(self.allIsotopolog, self.minimumIntensity)



















class RuleMatcher:

    def __init__(self, rules, ppm=5., log=False):
        self.rules=rules
        self.ppm=ppm
        self.log=log

        self.isos={"X": 0}
        for rule in self.rules:
            for iso in rule.getAllRequiredIsotopologs():
                if iso not in self.isos.keys():
                    self.isos[iso]=_formulaTools.calcIsotopologOffsetWeight(_formulaTools.parseFormula(iso))


    def matchIsoPatternRules(self, msScan, signalInd, charge=1):
        mz=msScan.mz_list[signalInd]

        isotopologDict={}
        noNeedToCheckFurther=[]

        for iso in self.isos.keys():
            if iso not in isotopologDict:
                b=Bunch(intensity=None)

                if iso=="X":
                    b.mz=mz
                else:
                    b.mz=mz+self.isos[iso]/charge

                bound=msScan.findMZ(b.mz, self.ppm)
                ind=msScan.getMostIntensePeak(bound[0], bound[1])

                if ind!=-1:
                    b.intensity=msScan.intensity_list[ind]

                    isotopologDict[iso]=b
                    noNeedToCheckFurther.append(ind)

        rulesValid=True
        for rule in self.rules:
            if rulesValid:
                ruleValid=rule.check(isotopologDict, log=self.log)
                rulesValid=rulesValid and ruleValid

        return rulesValid, noNeedToCheckFurther if rulesValid else []

    def checkChromPeaks(self, peakAreas):
        rulesValid=True
        for rule in self.rules:
            rulesValid=rulesValid and rule.check(peakAreas, log=self.log)

        return rulesValid


    def getChromatographicPeaks(self):
        chromPeaks={}
        for rule in self.rules:
            for iso in rule.getChromatographicPeaks():
                chromPeaks[iso]= Bunch(mzInc=_formulaTools.calcIsotopologOffsetWeight(_formulaTools.parseFormula(iso)), requiredChromPeak=True)

        return chromPeaks