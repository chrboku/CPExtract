#    MetExtract II
#    Copyright (C) 2015
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


from utils import getNormRatio
from utils import Bunch
from copy import deepcopy
from matchIsotopologPatternRules import RuleMatcher


def getSubstitutionArray(purity, xMax, maxSub):
    ret = []
    ret.append([-1 for n in range(0, maxSub + 1)])
    for i in range(1, xMax + 1):
        cur = []
        for j in range(0, maxSub + 1):
            cur.append(getNormRatio(purity, i, j))
        ret.append(cur)
    return ret


# results object (ion signal pairs)
class mzFeature:
    def __init__(self, mz, similarityString, scanIndex, loading, nIntensity, ionMode, type="pair", otherIsotopologs=None):
        self.mz = mz
        self.similarityString = similarityString
        self.scanIndex = scanIndex
        self.loading = loading
        self.nIntensity = nIntensity
        self.ionMode = ionMode
        self.type=type
        if otherIsotopologs==None:
            otherIsotopologs={}
        self.otherIsotopologs=otherIsotopologs

    def __str__(self):
        return "mz: %.4f, lmz: %.4f, xCount: %d, charge: %d, NInt: %.1f, LInt: %.1f, ionMode: %s, scanIndex: %d"%(
                    self.mz, self.lmz, self.xCount, self.loading, self.nIntensity, self.lIntensity, self.ionMode, self.scanIndex)

## calculates all possible combinations of labeling elements
## required for double labeling
def getCombinationsOfLabel(useElems, labelingElements, minLabelingAtoms, maxLabelingAtoms, used=None, startAt=0,
                           ind=2):
    combs = []
    if used is None:
        used = {}

    x = sum(used.values())
    if minLabelingAtoms <= x <= maxLabelingAtoms:
        b = Bunch(atoms=deepcopy(used), atomsCount=sum([used[e] for e in used.keys()]), mzdelta=sum(
            [(labelingElements[e].massLabeled - labelingElements[e].massNative) * used[e] for e in used.keys()]))
        combs.append(b)

    if startAt < len(useElems):
        for cStartAt in range(startAt, len(useElems)):
            e = useElems[cStartAt]
            for i in range(labelingElements[e].minXn, labelingElements[e].maxXn + 1):
                c = deepcopy(used)
                c[e] = i

                combs.extend(
                    getCombinationsOfLabel(useElems, labelingElements, minLabelingAtoms, maxLabelingAtoms, c,
                                           cStartAt + 1, ind=ind + 2))

    return combs


# detects in each recorded MS scan (lvl. 1) isotope patterns originating from a native and a (partially) labelled
# metabolite / biotransformation product. It is important, that all atoms that can be labelled are actually
# labelled and have the same chance to be labelled. Thus, fluxomics applications or such isotope patterns
# are in general not supported. They can, however, be detected if not isotopolog verification step is
# used (not recommended)
def matchPartners(mzXMLData, rules,
                  labellingIsotopeB, useCIsotopePatternValidation, intensityThres, maxLoading, xCounts, xOffset, ppm,
                  purityN, purityL, startTime, stopTime, filterLine, ionMode,
                  metabolisationExperiment,
                  reportFunction=None, writeExtendedDiagnostics=True, searchDirection="asc"):
    scanRTRange = stopTime - startTime

    cValidationOffset = 1.00335484   # mass difference between 12C and 13C

    detectedIonPairs = []

    oriXOffset = xOffset
    oriCValidationOffset = cValidationOffset


    # substitution arrays for checking the number of carbon atoms
    maxSub=5
    if len(xCounts)==0:
        xCounts=[1,2,3]
    purityNArray = getSubstitutionArray(purityN, max(xCounts)+1, maxSub)   # native metabolite
    purityLArray = getSubstitutionArray(purityL, max(xCounts)+1, maxSub)   # labelled metabolite


    ruleMatcher=RuleMatcher(rules, ppm=ppm, log=False)

    # iterate over all MS scans (lvl. 1)
    curScanIndex=0
    for j in range(0, len(mzXMLData.MS1_list)):

        try:
            curScan = mzXMLData.MS1_list[j]
            curScanDetectedIonPairs = []

            # check for correct filterline and scan time
            if curScan.filter_line == filterLine:
                if startTime <= (curScan.retention_time / 60.) <= stopTime:

                    if reportFunction is not None:
                        reportFunction((curScan.retention_time / 60. - startTime) / scanRTRange,
                                       "FilterLine: %s, RT %.2f" % (filterLine, curScan.retention_time / 60.))

                    dontUsePeakIndices = []

                    # assume each peak to be a valid M (monoisotopic, native metabolite ion)
                    # and verify this assumption (search for (partially) labelled pendant)

                    indices=range(0, len(curScan.mz_list))
                    if searchDirection == "asc":
                        pass
                    else:
                        indices = range(len(curScan.mz_list)-1, 0-1, -1)

                    for currentPeakIndex in indices:
                        if not (currentPeakIndex in dontUsePeakIndices):
                            curPeakmz = curScan.mz_list[currentPeakIndex]
                            curPeakIntensity = curScan.intensity_list[currentPeakIndex]

                            curPeakDetectedIonPairs = []

                            # only consider peaks above the threshold
                            if curPeakIntensity >= intensityThres:
                                skipOtherLoadings = False

                                ## do not process peaks that are likely isotopologs
                                backIsos=[]
                                for l in range(maxLoading, 0, -1):
                                    iso = curScan.findMZ(curPeakmz - oriCValidationOffset / l, ppm)
                                    iso = curScan.getMostIntensePeak(iso[0], iso[1])

                                    if iso != -1 and curScan.intensity_list[iso]>curPeakIntensity:
                                        backIsos.append(l)
                                if len(backIsos)>0:
                                    continue


                                possibleLoadings=[]
                                ## figure out possible loadings
                                for l in range(maxLoading, 0, -1):
                                    iso = curScan.findMZ(curPeakmz + oriCValidationOffset / l, ppm, start=currentPeakIndex)
                                    iso = curScan.getMostIntensePeak(iso[0], iso[1])

                                    if iso != -1:
                                        possibleLoadings.append(l)
                                        break ## skip other loadings



                                if len(possibleLoadings)==0:
                                    possibleLoadings=[1]


                                for curLoading in possibleLoadings:
                                    if not skipOtherLoadings:
                                        xOffset = oriXOffset / float(curLoading)
                                        cValidationOffset = oriCValidationOffset / float(curLoading)


                                        # labeling patters derived from small tracer compounds (e.g. 13C)
                                        # Idea for the detection has been developed by Bernhard Seidl for the detecton of polyketides
                                        #
                                        # E.g. tracer is [13C]2
                                        # |
                                        # |
                                        # |   |
                                        # |   |   |
                                        # | | |   |   |
                                        # | | | | |   |   |
                                        # | | | | | | | | |
                                        #
                                        # E.g. [13C]3D3
                                        # |  <-- [13C]3D3 -->  |
                                        # ||                  ||
                                        # |||                ||||
                                        #
                                        # NOTE: these isotope patterns are strange
                                        if useCIsotopePatternValidation == 3:

                                            rulesValid, noNeedToCheckFurther=ruleMatcher.matchIsoPatternRules(curScan, currentPeakIndex, curLoading)

                                            if rulesValid:
                                                curPeakDetectedIonPairs.append(
                                                    mzFeature(mz=curPeakmz,
                                                              similarityString="CP",
                                                              scanIndex=curScanIndex,
                                                              loading=curLoading,
                                                              nIntensity=curPeakIntensity,
                                                              ionMode=ionMode,

                                                              type="CP",
                                                              otherIsotopologs=ruleMatcher.getChromatographicPeaks()))

                                                dontUsePeakIndices.extend(noNeedToCheckFurther)

                            curScanDetectedIonPairs.extend(curPeakDetectedIonPairs)

                if len(curScanDetectedIonPairs)>0 and False:
                    from utils import printObjectsAsTable
                    print "\n"
                    print curScan.retention_time/60.
                    printObjectsAsTable(curScanDetectedIonPairs, attrs=["mz", "xCount", "loading", "nIntensity", "lIntensity", "ionMode"])

                detectedIonPairs.extend(curScanDetectedIonPairs)
                curScanIndex += 1

        except Exception as e:
            import traceback
            traceback.print_exc()

    return detectedIonPairs


















