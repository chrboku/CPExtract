from formulaTools import formulaTools
from utils import Bunch

import abc

_formulaTools=formulaTools()

class Rule:
    def isStatic(self):
        return True

    def getAllRequiredIsotopologs(self):
        return []

    def check(self, isotopologDict, log=False):
        return False  ## or return False, str for dynamic information

    def getChromatographicPeaks(self):
        return []  ## or return False, str for dynamic information

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

    def isStatic(self):
        return True

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



class RatioRule(Rule):
    def __init__(self, numeratorA="[13C]2+X", denominatorA="[13C]1", numeratorB=None, denominatorB=None, ratio = 1.0, ratioWindowMultiplier=[0.666, 1.5], mustBePresent=True):
        self.numeratorA = numeratorA
        self.denominatorA = denominatorA
        self.numeratorB = numeratorB
        self.denominatorB = denominatorB
        self.ratio = ratio
        self.ratioWindowMultiplier = ratioWindowMultiplier
        self.mustBePresent = mustBePresent

    def isStatic(self):
        return True

    def getAllRequiredIsotopologs(self):
        a = []
        a.extend(self.numeratorA.split("+"))
        a.extend(self.denominatorA.split("+"))
        if self.numeratorB is not None:
            a.extend(self.numeratorB.split("+"))
        if self.denominatorB is not None:
            a.extend(self.denominatorB.split("+"))
        return a

    def check(self, isotopologDict, log=False):
        intNA = 0
        intDA = 0

        for iso in self.numeratorA.split("+"):
            if iso in isotopologDict.keys():
                intNA = intNA + isotopologDict[iso].intensity
            elif self.mustBePresent: return False
        for iso in self.denominatorA.split("+"):
            if iso in isotopologDict.keys():
                intDA = intDA + isotopologDict[iso].intensity
            elif self.mustBePresent: return False

        if intNA == 0 or intDA == 0:
            return not self.mustBePresent

        ratA = intNA / intDA

        if self.numeratorB is not None and self.denominatorB is not None:
            intNB = 0
            intDB = 0

            for iso in self.numeratorB.split("+"):
                if iso in isotopologDict.keys():
                    intNB = intNB + isotopologDict[iso].intensity
                elif self.mustBePresent: return False
            for iso in self.denominatorB.split("+"):
                if iso in isotopologDict.keys():
                    intDB = intDB + isotopologDict[iso].intensity
                elif self.mustBePresent: return False

            if intNB == 0 or intNB == 0:
                return not self.mustBePresent

            ratB = intNB / intDB
        else:
            ratB = self.ratio

        return self.ratioWindowMultiplier[0] <= ratA/ratB <= self.ratioWindowMultiplier[1]

    def getChromatographicPeaks(self):
        return []

    def getMessage(self):
        if self.numeratorB is not None and self.denominatorB is not None:
            return "Ratio test for (%s)/(%s) is (%s)/(%s) with allowed deviation of %.3f - %.3f"%(self.numeratorA, self.denominatorA, self.numeratorB, self.denominatorB, self.ratioWindowMultiplier[0], self.ratioWindowMultiplier[1])
        return "Ratio test for (%s)/(%s) is %.3f with allowed deviation of %.3f - %.3f"%(self.numeratorA, self.denominatorA, self.ratio, self.ratioWindowMultiplier[0], self.ratioWindowMultiplier[1])








class AbsenceRule(Rule):
    def __init__(self, otherIsotopolog="[13C]-1", maxRatio=0.1):
        self.otherIsotopolog=otherIsotopolog
        self.maxRatio=maxRatio

    def isStatic(self):
        return True

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

    def isStatic(self):
        return True

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

    def isStatic(self):
        return True

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
        dynamicInfo = []
        for rule in self.rules:
            if rulesValid:
                if rule.isStatic():
                    temp = rulesValid and rule.check(isotopologDict, log=self.log)
                else:
                    temp, dynInfo = rule.check(isotopologDict, log=self.log)
                    dynamicInfo.append(dynInfo)
                rulesValid=rulesValid and temp

        return rulesValid, noNeedToCheckFurther if rulesValid else [], "Static" if len(dynamicInfo) == 0 else ";".join(dynamicInfo)

    def checkChromPeaks(self, peakAreas):
        rulesValid=True
        dynamicInfo = []
        for rule in self.rules:
            if rulesValid:
                if rule.isStatic():
                    temp = rule.check(peakAreas, log=self.log)
                else:
                    temp, dynInfo = rule.check(peakAreas, log=self.log)
                    dynamicInfo.append(dynInfo)
                rulesValid = rulesValid and temp

        return rulesValid, "Static" if len(dynamicInfo) == 0 else ";".join(dynamicInfo)


    def getChromatographicPeaks(self):
        chromPeaks={}
        for rule in self.rules:
            for iso in rule.getChromatographicPeaks():
                chromPeaks[iso] = Bunch(mzInc=_formulaTools.calcIsotopologOffsetWeight(_formulaTools.parseFormula(iso)), requiredChromPeak=True)

        return chromPeaks