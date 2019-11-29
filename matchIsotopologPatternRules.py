from formulaTools import formulaTools
from utils import Bunch

import abc

_formulaTools=formulaTools()

class Rule:
    def getAllRequiredIsotopologs(self):
        return {}

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
            return True

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
                b=Bunch(isoInd=None, intensity=None, mz=None)

                if iso=="X":
                    b.mz=mz
                else:
                    b.mz=mz+self.isos[iso]/charge

                bound=msScan.findMZ(b.mz, self.ppm)
                ind=msScan.getMostIntensePeak(bound[0], bound[1])

                if ind!=-1:
                    b.isoInd=ind
                    b.intensity=msScan.intensity_list[b.isoInd]

                    isotopologDict[iso]=b
                    noNeedToCheckFurther.append(ind)

        passedRules=True
        for rule in self.rules:
            if passedRules:
                passedRules=passedRules and rule.check(isotopologDict, log=self.log)

        return passedRules, noNeedToCheckFurther if passedRules else []

    def getChromatographicPeaks(self):
        chromPeaks={}
        for rule in self.rules:
            for iso in rule.getChromatographicPeaks():
                chromPeaks[iso]= Bunch(mzInc=_formulaTools.calcIsotopologOffsetWeight(_formulaTools.parseFormula(iso)), requiredChromPeak=True)

        return chromPeaks