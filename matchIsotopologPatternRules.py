from formulaTools import formulaTools
from utils import Bunch


ruleChecker=formulaTools()

class Rule:
    def getAllRequiredIsotopologs(self):
        return {}

    def check(self, isotopologDict, log=False):
        return False



class PresenceRule(Rule):
    def __init__(self, otherIsotopolog="[13C]2", minIntensity=10000, mustBePresent=True, ratioWindows={'X': [1/2., 2]}):
        self.otherIsotopolog=otherIsotopolog
        self.minIntensity=minIntensity
        self.mustBePresent=mustBePresent
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




class AbsenceRule(Rule):
    def __init__(self, otherIsotopolog="[13C]-1", maxIntensity=10000):
        self.otherIsotopolog=otherIsotopolog
        self.maxIntensity=maxIntensity

    def getAllRequiredIsotopologs(self):
        return [self.otherIsotopolog]

    def check(self, isotopologDict, log=False):
        if self.otherIsotopolog not in isotopologDict.keys():
            return True

        if self.otherIsotopolog in isotopologDict.keys():
            if isotopologDict[self.otherIsotopolog].intensity<self.maxIntensity:
                return True
            else:
                return False



def matchIsoPatternRules(msScan, signalInd, rules, charge=1, ppm=5, log=False):
    mz=msScan.mz_list[signalInd]

    isotopologDict={}
    noNeedToCheckFurther=[]

    for rule in rules:
        isos=rule.getAllRequiredIsotopologs()
        for iso in isos:
            if iso not in isotopologDict:
                b=Bunch(isoInd=None, intensity=None, mz=None)

                if iso=="X":
                    b.mz=mz
                else:
                    b.mz=mz+ruleChecker.calcIsotopologOffsetWeight(ruleChecker.parseFormula(iso))/charge

                bound=msScan.findMZ(b.mz, ppm)
                ind=msScan.getMostIntensePeak(bound[0], bound[1])

                if ind!=-1:
                    b.isoInd=ind
                    b.intensity=msScan.intensity_list[b.isoInd]

                    isotopologDict[iso]=b
                    noNeedToCheckFurther.append(ind)

    passedRules=True
    for rule in rules:
        if passedRules:
            passedRules=passedRules and rule.check(isotopologDict, log=log)

    return passedRules, noNeedToCheckFurther if passedRules else []