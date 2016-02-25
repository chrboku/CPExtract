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

from formulaTools import formulaTools


maxSub = 5


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
    def __init__(self, mz, lmz, deltamzTheoretical, xCount, ratioNative, ratioLabeled, scanIndex, loading, nIntensity, lIntensity, ionMode):
        self.mz = mz
        self.lmz = lmz
        self.deltamzTheoretical = deltamzTheoretical
        self.xCount = xCount
        self.ratioNative = ratioNative
        self.ratioLabeled = ratioLabeled
        self.scanIndex = scanIndex
        self.loading = loading
        self.nIntensity = nIntensity
        self.lIntensity = lIntensity
        self.ionMode = ionMode

    def __str__(self):
        return "mz: %.4f, lmz: %.4f, xCount: %d, charge: %d, NInt: %.1f, LInt: %.1f, ionMode: %s, scanIndex: %d"%(
                    self.mz, self.lmz, self.xCount, self.loading, self.nIntensity, self.lIntensity, self.ionMode, self.scanIndex)






## calculates all possible combinations of labeling elements
## required for double labeling
def getCombinationsOfLabel(useElems, labelingElements, minLabelingAtoms, maxLabelingAtoms, used=None, startAt=0, ind=2):
        combs=[]
        if used is None:
            used = {}

        x=sum(used.values())
        if minLabelingAtoms <= x <= maxLabelingAtoms:
            b=Bunch(atoms=deepcopy(used), atomsCount=sum([used[e] for e in used.keys()]), mzdelta=sum([(labelingElements[e].massLabeled-labelingElements[e].massNative)*used[e] for e in used.keys()]))
            combs.append(b)

        if startAt < len(useElems):
            for cStartAt in range(startAt, len(useElems)):
                e=useElems[cStartAt]
                for i in range(labelingElements[e].minXn, labelingElements[e].maxXn+1):
                    c=deepcopy(used)
                    c[e]=i

                    combs.extend(getCombinationsOfLabel(useElems, labelingElements, minLabelingAtoms, maxLabelingAtoms, c, cStartAt+1, ind=ind+2))

        return combs



# detects in each recorded MS scan (lvl. 1) isotope patterns originating from a native and a (partially) labelled
# metabolite / biotransformation product. It is important, that all atoms that can be labelled are actually
# labelled and have the same chance to be labelled. Thus, fluxomics applications or such isotope patterns
# are in general not supported. They can, however, be detected if not isotopolog verification step is
# used (not recommended)
def matchPartners(mzXMLData, labellingElement, useCValidation, intensityThres, isotopologIntensityThres, maxLoading, xMin, xMax, xOffset, ppm,
                  intensityErrorN, intensityErrorL, purityN, purityL, startTime, stopTime, filterLine, ionMode,
                  peakCountLeft, peakCountRight, lowAbundanceIsotopeCutoff, metabolisationExperiment, checkRatio,
                  minRatio, maxRatio, reportFunction=None):

    fT=formulaTools()


    scanRTRange = stopTime - startTime

    carbonMassOffset = 1.00335   # mass difference between 12C and 13C

    detectedIonPairs = []


    #### for testing only. labellingElement is overwritten with a new object (Bunch)
    labelingElements={}
    #labelingElements[labellingElement]=Bunch(nativeIsotope="12C", labelingIsotope="13C", massNative=0      , massLabeled=0+xOffset , isotopicEnrichmentNative=purityN, isotopicEnrichmentLabeled=purityL, minXn=xMin, maxXn=xMax)
    labelingElements["C"]=Bunch(nativeIsotope="12C", labelingIsotope="13C", massNative=12.      , massLabeled=13.00335 , isotopicEnrichmentNative=0.9893, isotopicEnrichmentLabeled=0.995, minXn=1, maxXn=2)
    labelingElements["H"]=Bunch(nativeIsotope= "1H", labelingIsotope= "2H", massNative=1.0078250, massLabeled=2.0141018, isotopicEnrichmentNative=0.9999, isotopicEnrichmentLabeled=0.96 , minXn=1, maxXn=3)

    minLabelingAtoms=3
    maxLabelingAtoms=5


    ## substitution arrays for checking the number of carbon atoms
    for elem in labelingElements.keys():
        b=labelingElements[elem]
        b.purityNArray = getSubstitutionArray(b.isotopicEnrichmentNative , b.maxXn + 3, maxSub)   # native metabolite
        b.purityLArray = getSubstitutionArray(b.isotopicEnrichmentLabeled, b.maxXn + 3, maxSub)   # labelled metabolite

    ## combinations of labeling elements
    tempCombs=getCombinationsOfLabel(["C", "H"], labelingElements, minLabelingAtoms, maxLabelingAtoms)
    combs=[]
    for comb in tempCombs:
        if comb.atomsCount==3 and comb.atoms["H"]==3:
            pass
        else:
            combs.append(comb)

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
                                       "RT %.2f" % (curScan.retention_time / 60.))

                    dontUsePeakIndices = []

                    # assume each peak to be a valid M (monoisotopic, native metabolite ion)
                    # and verify this assumption (search for (partially) labelled pendant)
                    for currentPeakIndex in range(0, len(curScan.mz_list)):
                        if not (currentPeakIndex in dontUsePeakIndices):
                            curPeakmz = curScan.mz_list[currentPeakIndex]
                            curPeakIntensity = curScan.intensity_list[currentPeakIndex]

                            curPeakDetectedIonPairs = []

                            # only consider peaks above the threshold
                            if curPeakIntensity >= intensityThres:
                                skipOtherLoadings = False


                                possibleLoadings=[]
                                ## figure out possible loadings
                                for l in range(1, maxLoading+1, 1):
                                    iso = curScan.findMZ(curPeakmz + carbonMassOffset / l, ppm, start=currentPeakIndex)
                                    iso = curScan.getMostIntensePeak(iso[0], iso[1])

                                    if iso != -1:
                                        possibleLoadings.append(l)

                                if len(possibleLoadings)==0:
                                    possibleLoadings=[1]


                                for curLoading in possibleLoadings:
                                    if not skipOtherLoadings:


                                        # C-isotope distribution validation for labelling with N, S, ... (useCValidation == 2)
                                        # region
                                        # The carbon distribution of both isotopologs is checked for equality
                                        # checks if the isotope patterns of M, M+1.. and M', M'+1.. are approximately the same
                                        # E.g. 15N-labelling
                                        # |    <--   Nn   -->    |
                                        # ||                     ||
                                        # |||                   ||||
                                        # Required for some labelling applications (e.g. S, N, Cl)
                                        # Requires: - a high resolution and separation of different isotoplogs (especially carbon)
                                        # EXPERIMENTAL: has not been tested with real data (not N or S labelled sample material
                                        #               was available)
                                        if useCValidation == 2:

                                            # search for M+1 isotoplog (if required it needs to be present; slight speed gain)
                                            isoM_p1 = curScan.findMZ(curPeakmz + carbonMassOffset/curLoading, ppm, start=currentPeakIndex)
                                            isoM_p1 = curScan.getMostIntensePeak(isoM_p1[0], isoM_p1[1])
                                            if isoM_p1 != -1:
                                                dontUseXCount = []

                                                # test certain number of labelled atoms
                                                for xCount in range(xMax, xMin - 1, -1):
                                                    if not (xCount in dontUseXCount):

                                                        # search for M'
                                                        isoM_pX = curScan.findMZ(curPeakmz + xCount * xOffset / curLoading, ppm,
                                                                                 start=currentPeakIndex)
                                                        isoM_pX = curScan.getMostIntensePeak(isoM_pX[0], isoM_pX[1],
                                                                                             intensityThres)
                                                        if isoM_pX != -1:
                                                            # check if ratio is okay
                                                            rat=curPeakIntensity/curScan.intensity_list[isoM_pX]
                                                            if not(checkRatio) or minRatio<=rat<=maxRatio:

                                                                # search for M'+1
                                                                isoM_pXp1 = curScan.findMZ(
                                                                    curPeakmz + (xCount * xOffset + carbonMassOffset) / curLoading, ppm,
                                                                    start=currentPeakIndex)
                                                                isoM_pXp1 = curScan.getMostIntensePeak(isoM_pXp1[0],
                                                                                                       isoM_pXp1[1])
                                                                if isoM_pXp1 != -1:
                                                                    acceptIsoPeak = []

                                                                    if not (isoM_p1 in dontUsePeakIndices):
                                                                        acceptIsoPeak.append(
                                                                            curScan.intensity_list[isoM_p1] / curPeakIntensity)
                                                                    acceptIsoPeakO = []

                                                                    # check if ratios of M+1/M and M'+1/M' are approximately
                                                                    # the same (<=intensityErrorL)
                                                                    if len(acceptIsoPeak) > 0:

                                                                        if not (isoM_pX in dontUsePeakIndices):
                                                                            if not (isoM_pXp1 in dontUsePeakIndices):
                                                                                for isoRatio in acceptIsoPeak:
                                                                                    if abs(curScan.intensity_list[isoM_pXp1] /
                                                                                                   curScan.intensity_list[
                                                                                                       isoM_pX] - isoRatio) < intensityErrorL:
                                                                                        acceptIsoPeakO.append(
                                                                                            (isoM_pX, isoM_pXp1))
                                                                    if len(acceptIsoPeak) > 0 and len(acceptIsoPeakO) > 0:
                                                                        #mzs.append(
                                                                        #    [mz, xCount, curScanNum, curLoading, intensity])
                                                                        curScanDetectedIonPairs.append(
                                                                            mzFeature(mz=curPeakmz,
                                                                                      lmz=curScan.mz_list[isoM_pX],
                                                                                      xCount=xCount,
                                                                                      scanIndex=curScanIndex,
                                                                                      loading=curLoading,
                                                                                      nIntensity=curPeakIntensity,
                                                                                      lIntensity=curScan.intensity_list[acceptIsoPeakO[0][0]],
                                                                                      ionMode=ionMode))
                                                                        skipOtherLoadings = True
                                        # endregion
                                        # Mixed isotope validation (useCValidation == 1)
                                        # region
                                        # First, the distinct isotope patterns of the native and
                                        # x-labelled metabolite forms are tested. If these do not successfully validate
                                        # a ion signal pairing, C-validation (similar carbon isotope pattern) is tested
                                        # Required for some labelling applications (e.g. S, N, Cl)
                                        # Requires: - a high resolution and separation of different isotoplogs (especially carbon)
                                        # EXPERIMENTAL: has not been tested with real data (not N or S labelled sample material
                                        #               was available)
                                        # CAUTION: Not tested, partially implemented. Use with caution
                                        if useCValidation == 1:
                                            #Mixed Isotope Validation
                                            isoM_p1 = curScan.findMZ(curPeakmz + xOffset / curLoading, ppm, start=currentPeakIndex)
                                            isoM_p1 = curScan.getMostIntensePeak(isoM_p1[0], isoM_p1[1])
                                            if peakCountLeft == 1 or isoM_p1 != -1:
                                                dontUseXCount = []
                                                for xCount in range(xMax, xMin - 1, -1):
                                                    if not (xCount in dontUseXCount):
                                                        isoM_pX = curScan.findMZ(curPeakmz + xCount * xOffset / curLoading, ppm,
                                                                                 start=currentPeakIndex)
                                                        isoM_pX = curScan.getMostIntensePeak(isoM_pX[0], isoM_pX[1],
                                                                                             intensityThres)
                                                        if isoM_pX != -1:
                                                            isoM_pXm1 = curScan.findMZ(curPeakmz + (xCount - 1) * xOffset / curLoading, ppm,
                                                                                       start=currentPeakIndex)
                                                            isoM_pXm1 = curScan.getMostIntensePeak(isoM_pXm1[0],
                                                                                                   isoM_pXm1[1])
                                                            if peakCountRight == 1 or isoM_pXm1[0] != -1:
                                                                acceptIsoPeak = []
                                                                normRatio = purityNArray[xCount][1]
                                                                if peakCountLeft > 1:
                                                                    if metabolisationExperiment:
                                                                        isoM_pXp1 = curScan.findMZ(
                                                                            curPeakmz + (xCount + 1) * xOffset / curLoading, ppm,
                                                                            start=currentPeakIndex)
                                                                        isoM_pXp1 = curScan.getMostIntensePeak(
                                                                            isoM_pXp1[0], isoM_pXp1[1])


                                                                        if not (isoM_p1 in dontUsePeakIndices):
                                                                            if abs(curScan.intensity_list[
                                                                                       isoM_p1] / curPeakIntensity - normRatio -
                                                                                           curScan.intensity_list[
                                                                                               isoM_pXp1]) < intensityErrorN:
                                                                                acceptIsoPeak.append(isoM_p1)
                                                                    else:

                                                                        if not (isoM_p1 in dontUsePeakIndices):
                                                                            if abs(curScan.intensity_list[
                                                                                       isoM_p1] / curPeakIntensity - normRatio) < intensityErrorN:
                                                                                acceptIsoPeak.append(isoM_p1)
                                                                else:
                                                                    acceptIsoPeak.append(1)
                                                                acceptIsoPeakO = []
                                                                normRatio = purityLArray[xCount][1]

                                                                if peakCountRight > 1:

                                                                    if not (isoM_pX in dontUsePeakIndices):
                                                                        if not (isoM_pXm1 in dontUsePeakIndices):
                                                                            if abs(curScan.intensity_list[isoM_pXm1] /
                                                                                           curScan.intensity_list[
                                                                                               isoM_pX] - normRatio) < intensityErrorL:
                                                                                acceptIsoPeakO.append((isoM_pX, isoM_pXm1))
                                                                else:
                                                                    acceptIsoPeakO.append(1)
                                                                if (peakCountLeft == 1 or len(acceptIsoPeak) > 0) and (
                                                                                peakCountRight == 1 or len(
                                                                            acceptIsoPeakO) > 0):
                                                                    #mzs.append(
                                                                    #    [mz, xCount, curScanNum, curLoading, intensity])

                                                                    curScanDetectedIonPairs.append(
                                                                        mzFeature(mz=curPeakmz,
                                                                                  xCount=xCount,
                                                                                  lmz=curScan.mz_list[isoM_pX],
                                                                                  scanIndex=curScanIndex,
                                                                                  loading=curLoading,
                                                                                  nIntensity=curPeakIntensity,
                                                                                  lIntensity=curScan.intensity_list[acceptIsoPeakO[0][0]],
                                                                                  ionMode=ionMode))
                                                                    skipOtherLoadings = True
                                                                    #dontUseXCount.append(xCount-1)

                                                            else:
                                                                isoM_pXm1 = curScan.findMZ(
                                                                    curPeakmz + (xCount + 1) * carbonMassOffset / curLoading, ppm,
                                                                    start=currentPeakIndex)
                                                                isoM_pXm1 = curScan.getMostIntensePeak(isoM_pXm1[0],
                                                                                                       isoM_pXm1[1])
                                                                if isoM_pXm1 != -1:
                                                                    acceptIsoPeak = []

                                                                    if not (isoM_p1 in dontUsePeakIndices):
                                                                        acceptIsoPeak.append(
                                                                            curScan.intensity_list[isoM_p1] / curPeakIntensity)
                                                                    acceptIsoPeakO = []

                                                                    if len(acceptIsoPeak) > 0:
                                                                        if not (isoM_pX in dontUsePeakIndices):
                                                                            if not (isoM_pXm1 in dontUsePeakIndices):
                                                                                for isoRatio in acceptIsoPeak:
                                                                                    if abs(curScan.intensity_list[
                                                                                               isoM_pXm1] /
                                                                                                   curScan.intensity_list[
                                                                                                       isoM_pX] - isoRatio) < intensityErrorL:
                                                                                        acceptIsoPeakO.append(
                                                                                            (isoM_pX, isoM_pXm1))
                                                                    if len(acceptIsoPeak) > 0 and len(
                                                                            acceptIsoPeakO) > 0:
                                                                        #mzs.append([mz, xCount, curScanNum, curLoading,
                                                                        #            intensity])

                                                                        curScanDetectedIonPairs.append(
                                                                            mzFeature(mz=curPeakmz,
                                                                                      lmz=curScan.mz_list[isoM_pX],
                                                                                      xCount=xCount,
                                                                                      scanIndex=curScanIndex,
                                                                                      loading=curLoading,
                                                                                      nIntensity=curPeakIntensity,
                                                                                      lIntensity=curScan.intensity_list[acceptIsoPeakO[0][0]],
                                                                                      ionMode=ionMode))
                                                                        skipOtherLoadings = True
                                        # endregion
                                        # Isotope pattern validation (useCValidation == 0)
                                        # region
                                        # It is tested, if the expected isotope patterns
                                        # separate for the native and labelled metabolite follow a theoretical pattern
                                        # E.g. native 12C and uniformly / partially 13C-labelled metabolites
                                        # |    <--   Cn   -->    |
                                        # ||                    ||
                                        # |||                  ||||
                                        # Necessary mainly for 13C-labelling with mirror-symmetric isotope patterns
                                        # NOTE: - Approach is mainly used for 13C-labelling
                                        if useCValidation==0:
                                            # find M+1 peak
                                            isoM_p1 = curScan.findMZ(curPeakmz + carbonMassOffset / curLoading, ppm, start=currentPeakIndex)
                                            isoM_p1 = curScan.getMostIntensePeak(isoM_p1[0], isoM_p1[1])
                                            if isoM_p1 != -1 or peakCountLeft == 1 or lowAbundanceIsotopeCutoff:
                                                # test certain number of labelled carbon atoms

                                                for comb in combs:
                                                    # find corresponding M' peak
                                                    isoM_pX = curScan.findMZ(curPeakmz + comb.mzdelta / curLoading, ppm, start=currentPeakIndex)
                                                    isoM_pX = curScan.getMostIntensePeak(isoM_pX[0], isoM_pX[1], intensityThres)
                                                    if isoM_pX != -1:

                                                        labPeakmz=curScan.mz_list[isoM_pX]
                                                        labPeakIntensity=curScan.intensity_list[isoM_pX]

                                                        # (1.) check if M' and M ratio are as expected
                                                        if checkRatio:
                                                            rat = curPeakIntensity / labPeakIntensity
                                                            if minRatio <= rat <= maxRatio:
                                                                pass     ## ratio check passed
                                                            else:
                                                                continue ## ratio check not passed

                                                        ## no isotopolog verification needs to be performed
                                                        if peakCountLeft == 1 and peakCountRight == 1:
                                                            curPeakDetectedIonPairs.append(
                                                                mzFeature(mz=curPeakmz,
                                                                          lmz=curScan.mz_list[isoM_pX],
                                                                          deltamzTheoretical=comb.mzdelta / curLoading,
                                                                          xCount=fT.flatToString(comb.atoms),
                                                                          ratioNative=0,
                                                                          ratioLabeled=0,
                                                                          scanIndex=curScanIndex,
                                                                          loading=curLoading,
                                                                          nIntensity=curPeakIntensity,
                                                                          lIntensity=labPeakIntensity,
                                                                          ionMode=ionMode))

                                                            skipOtherLoadings = True
                                                            continue

                                                        # find M'-1 peak
                                                        isoM_pXm1 = curScan.findMZ(curPeakmz + (comb.mzdelta - carbonMassOffset) / curLoading, ppm, start=currentPeakIndex)
                                                        isoM_pXm1 = curScan.getMostIntensePeak(isoM_pXm1[0], isoM_pXm1[1])
                                                        normRatioL = labelingElements["C"].purityLArray[comb.atomsCount][1]
                                                        normRatioN = labelingElements["C"].purityNArray[comb.atomsCount][1]


                                                        # 2. check if the observed M'-1/M' ratio fits the theoretical one
                                                        if peakCountRight > 1 or (lowAbundanceIsotopeCutoff and labPeakIntensity*normRatioL <= isotopologIntensityThres):
                                                            if isoM_pXm1 != -1:
                                                                observedRatioMp=curScan.intensity_list[isoM_pXm1] / labPeakIntensity
                                                                if abs(normRatioL-observedRatioMp) <= intensityErrorL:
                                                                    pass     ## accept peak
                                                                else:
                                                                    continue ## discard current peak
                                                            elif lowAbundanceIsotopeCutoff:
                                                                if labPeakIntensity*normRatioL <= isotopologIntensityThres:
                                                                    pass     ## accept peak
                                                                else:
                                                                    continue ## discard current peak
                                                            else:
                                                                continue     ## discard current peak
                                                        else:
                                                            pass ## accept peak

                                                        # 3. check if the observed M+1/M ratio fits the theoretical one
                                                        if peakCountLeft > 1 or (lowAbundanceIsotopeCutoff and curPeakIntensity*normRatioN <= isotopologIntensityThres):
                                                            observedRatioM = curScan.intensity_list[isoM_p1] / curPeakIntensity
                                                            adjRatio=0

                                                            if metabolisationExperiment:
                                                                # if experiment is a tracer-fate study, assume a conjugated moiety and correct M+1/M ratio for it
                                                                isoM_pXp1 = curScan.findMZ( curPeakmz + (comb.mzdelta + carbonMassOffset) / curLoading, ppm, start=currentPeakIndex)
                                                                isoM_pXp1 = curScan.getMostIntensePeak( isoM_pXp1[0], isoM_pXp1[1])
                                                                if isoM_pXp1 != -1:
                                                                    adjRatio=curScan.intensity_list[isoM_pXp1] / labPeakIntensity

                                                            if abs(abs(observedRatioM-adjRatio)-normRatioN) <= intensityErrorN:
                                                                pass         ## acceptPeak
                                                            elif lowAbundanceIsotopeCutoff:
                                                                if curPeakIntensity*normRatioN <= isotopologIntensityThres:
                                                                    pass     ## accept peak
                                                                else:
                                                                    continue ## discard peak
                                                            else:
                                                                continue     ## discard current peak




                                                        # All verification criteria are passed, store the ion signal pair
                                                        # for further processing
                                                        curPeakDetectedIonPairs.append(
                                                            mzFeature(mz=curPeakmz,
                                                                      lmz=curScan.mz_list[isoM_pX],
                                                                      deltamzTheoretical=comb.mzdelta / curLoading,
                                                                      xCount=fT.flatToString(comb.atoms),
                                                                      ratioNative=normRatioN,
                                                                      ratioLabeled=normRatioL,
                                                                      scanIndex=curScanIndex,
                                                                      loading=curLoading,
                                                                      nIntensity=curPeakIntensity,
                                                                      lIntensity=labPeakIntensity,
                                                                      ionMode=ionMode))

                                                        skipOtherLoadings = True
                                                        # endregion

                            if False:  ## select best fit
                                if len(curPeakDetectedIonPairs)>0:
                                    bestFit=None
                                    bestFitPPMDiff=1000000000

                                    ## TODO select best fit based on isotopic pattern (e.g. intensity)


                                    for mt in curPeakDetectedIonPairs:
                                        if abs(mt.lmz-mt.mz-mt.deltamzTheoretical)*1000000./mt.mz < bestFitPPMDiff:
                                            bestFit=mt
                                            bestFitPPMDiff=abs(mt.lmz-mt.mz-mt.deltamzTheoretical)*1000000./mt.mz

                                    curScanDetectedIonPairs.append(bestFit)
                            else:     ## use all peak pairs
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


















