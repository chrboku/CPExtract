import logging

from sqlite3 import *
from operator import itemgetter
import os
from math import floor
import ast
from collections import defaultdict

from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.lib import colors
from reportlab.pdfgen import canvas
from reportlab.graphics import renderPDF

from XICAlignment import XICAlignment

from utils import ChromPeakFeature, getSubGraphs, Bunch, SQLInsert, natSort, getSubGraphsFromDictDict, CallBackMethod
from runIdentification import ChromPeakFeature, getDBSuffix
from TableUtils import TableUtils
from MZHCA import HierarchicalClustering, cutTreeSized

import base64
from pickle import dumps, loads

from math import isnan


import time

import HCA_general
from Chromatogram import Chromatogram
from utils import mean, sd, corr


import exportAsFeatureML

def get_main_dir():
    from utils import get_main_dir
    return os.path.join(get_main_dir(), '')


# HELPER METHOD for writing first page of PDF (unused)
def _writeFirstPage(pdf, groupSizePPM, maxTimeDeviation, align, nPolynom):
    curHeight = 780

    pdf.drawString(60, curHeight, "Max. Group Size")
    pdf.drawString(200, curHeight, "%.1f" % groupSizePPM)
    curHeight -= 20

    pdf.drawString(60, curHeight, "Max. Time Deviation")
    pdf.drawString(200, curHeight, "%.2f seconds" % maxTimeDeviation)
    curHeight -= 20

    pdf.drawString(60, curHeight, "Align")
    pdf.drawString(200, curHeight, "%s" % ("yes" if align else "no"))
    curHeight -= 20

    if align:
        pdf.drawString(60, curHeight, "Polynom Order")
        pdf.drawString(200, curHeight, "%d" % nPolynom)
        curHeight -= 20

    pdf.showPage()





# store used configuration to DB file
def writeConfigToDB(curs, align, file, groupSizePPM, maxLoading, maxTimeDeviation, xCounts, nPolynom,
                    negativeScanEvent, positiveScanEvent, rVersion, meVersion):
    SQLInsert(curs, "config", key="MEVersion", value=str(meVersion))
    SQLInsert(curs, "config", key="RVersion", value=str(rVersion))

    SQLInsert(curs, "config", key="FPBRACK_xCounts", value=str(xCounts))
    SQLInsert(curs, "config", key="FPBRACK_groupSizePPM", value=str(groupSizePPM))
    SQLInsert(curs, "config", key="FPBRACK_positiveScanEvent", value=str(positiveScanEvent))
    SQLInsert(curs, "config", key="FPBRACK_negativeScanEvent", value=str(negativeScanEvent))
    SQLInsert(curs, "config", key="FPBRACK_maxTimeDeviation", value=str(maxTimeDeviation))
    SQLInsert(curs, "config", key="FPBRACK_maxLoading", value=str(maxLoading))
    SQLInsert(curs, "config", key="FPBRACK_file", value=str(file))
    SQLInsert(curs, "config", key="FPBRACK_align", value=str(align))
    SQLInsert(curs, "config", key="FPBRACK_nPolynom", value=str(nPolynom))

# bracket results
def bracketResults(indGroups, xCounts, groupSizePPM, positiveScanEvent=None, negativeScanEvent=None,
                 maxTimeDeviation=0.36 * 60, maxLoading=1, file="./results.tsv", align=True, nPolynom=1,
                 pwMaxSet=None, pwValSet=None, pwTextSet=None, rVersion="", meVersion="", generalProcessingParams=Bunch(), start=0):


    # create results SQLite tables
    resDB=Bunch(conn=None, curs=None)
    if os.path.exists(file+getDBSuffix()) and os.path.isfile(file+getDBSuffix()):
        os.remove(file+getDBSuffix())
    resDB.conn=connect(file+getDBSuffix())
    resDB.curs=resDB.conn.cursor()

    try:

        cpf = get_main_dir() + "./XICAlignment.r"       # initialise chromatographic alignment script
        xicAlign = XICAlignment(cpf)


        resDB.curs.execute("DROP TABLE IF EXISTS GroupResults")
        resDB.curs.execute("CREATE TABLE GroupResults(id INTEGER PRIMARY KEY, mz FLOAT, rt FLOAT, similarityString TEXT, charge INTEGER, scanEvent TEXT, ionisationMode TEXT, OGroup INTEGER, Ion TEXT, LOSS TEXT, M TEXT)")

        resDB.curs.execute("DROP TABLE IF EXISTS FoundFeaturePairs")
        resDB.curs.execute("CREATE TABLE FoundFeaturePairs(resID INTEGER, file TEXT, featurePairID INTEGER, featureGroupID INTEGER, PeakScale FLOAT, " \
                            "Area FLOAT, featureType TEXT)")

        resDB.curs.execute("DROP TABLE IF EXISTS FileMapping")
        resDB.curs.execute("CREATE TABLE FileMapping(fileName TEXT, filePath TEXT, groupID INTEGER)")

        resDB.curs.execute("DROP TABLE IF EXISTS FileGroups")
        resDB.curs.execute("CREATE TABLE FileGroups(groupName TEXT, id INTEGER PRIMARY KEY)")

        resDB.curs.execute("DROP TABLE IF EXISTS config")
        resDB.curs.execute("CREATE TABLE config(id INTEGER PRIMARY KEY AUTOINCREMENT, key TEXT, value TEXT)")

        writeConfigToDB(resDB.curs, align, file, groupSizePPM, maxLoading, maxTimeDeviation, xCounts, nPolynom,
                        negativeScanEvent, positiveScanEvent, rVersion, meVersion)


        results = []

        # generate results for each input file
        i=1
        for key in natSort(indGroups.keys()):
            files=indGroups[key]
            SQLInsert(resDB.curs, "FileGroups", groupName=key, id=i)

            for ident in files:
                conn = connect(ident + getDBSuffix())
                curs = conn.cursor()
                fname=ident
                if ".mzxml" in ident.lower():
                    fname=ident[max(ident.rfind("/") + 1, ident.rfind("\\") + 1):ident.lower().rfind(".mzxml")]
                if ".mzml" in ident.lower():
                    fname=ident[max(ident.rfind("/") + 1, ident.rfind("\\") + 1):ident.lower().rfind(".mzml")]
                b=Bunch(filePath=ident, fileName=fname,featurePairs=[], conn=conn, curs=curs)

                results.append(b)

                SQLInsert(resDB.curs, "FileMapping", fileName=b.fileName, filePath=b.filePath, groupID=i)

            i+=1

        with open(file, "wb") as f:
            # initialise TSV results file

            f.write("Num\tComment\tMZ\tMZ_Range\tRT\tRT_Range\tPeakScales\tXn\tCharge\tScanEvent\tIonMode\tOGroup\tIon\tLoss\tM")
            for res in results:
                f.write("\t" + res.fileName + "_Area")

            f.write("\tdoublePeak")
            f.write("\n")

            # get ionisation modes scan events
            ionModes = {}
            if positiveScanEvent != "None":
                ionModes["+"] = positiveScanEvent
            if negativeScanEvent != "None":
                ionModes["-"] = negativeScanEvent
            assert len(ionModes) > 0

            totalChromPeaks = 0
            totalChromPeaksProc = 0
            tracersDeltaMZ = {}


            ## TODO this used the last open file. May be invalid. improve
            rows = []
            for row in res.curs.execute("SELECT key, value FROM config WHERE key='metabolisationExperiment'"):
                rows.append(str(row[1]))
            assert len(rows) == 1
            isMetabolisationExperiment = rows[0].lower() == "true"

            ## TODO this used the last open file. May be invalid. improve
            if isMetabolisationExperiment:
                for row in res.curs.execute("SELECT id, name, deltaMZ FROM tracerConfiguration"):
                    tracerName = str(row[1])
                    dmz = float(row[2])
                    if tracerName not in tracersDeltaMZ:
                        tracersDeltaMZ[tracerName] = dmz

                    if tracersDeltaMZ[tracerName] != dmz:
                        logging.warning("Warning: Tracers have not been configured identical in all measurement files")
            else:
                rows = []
                ## TODO this used the last open file. May be invalid. improve
                for row in res.curs.execute("SELECT key, value FROM config WHERE key='xOffset'"):
                    rows.append(str(row[1]))
                assert len(rows) == 1

                dmz = float(rows[0])
                if "FLE" not in tracersDeltaMZ:
                    tracersDeltaMZ["FLE"] = dmz

                if tracersDeltaMZ["FLE"] != dmz:
                    logging.warning("Warning: Tracers have not been configured identical in all measurement files")


            grp = 1

            curNum = 0

            mzXMLs={}
            xCounts=set()
            if pwMaxSet is not None: pwMaxSet.put(Bunch(mes="max", val=len(results)))
            if pwValSet is not None: pwValSet.put(Bunch(mes="value", val=0))
            for resi, res in enumerate(results):
                for row in res.curs.execute("SELECT DISTINCT similarityString FROM chromPeaks"):
                    xCounts.add(str(row[0]))
                if pwTextSet is not None: pwTextSet.put(Bunch(mes="text", val="Importing raw-file %s"%res.filePath))
                mzXML = Chromatogram()
                mzXML.parse_file(res.filePath, intensityCutoff=10000)

                mzXMLs[res.filePath]=mzXML

                if pwValSet is not None: pwValSet.put(Bunch(mes="value", val=resi))

            xCounts=sorted(list(xCounts))

            # xCountshelp = []
            # a=xCounts.replace(" ", "").replace(";", ",").split(",")
            # for j in a:
            #     if "-" in j:
            #         xCountshelp.extend(range(int(j[0:j.find("-")]), int(j[j.find("-")+1:len(j)])+1))
            #     else:
            #         xCountshelp.append(int(j))
            # xCounts=sorted(list(set(xCountshelp)))

            # prefetch number of bracketing steps for the progress bar
            totalThingsToDo=0
            for ionMode in ionModes:
                for tracer in tracersDeltaMZ:
                    for xCount in xCounts:
                        for cLoading in range(maxLoading, 0, -1):
                            totalThingsToDo+=1

            doneSoFar=0
            if pwMaxSet is not None: pwMaxSet.put(Bunch(mes="max", val=totalThingsToDo))

            # bracket data
            # process ionModes, tracers, xCount and loading separately
            for ionMode in ionModes:
                scanEvent = ionModes[ionMode]
                for xCount in xCounts:
                    for cLoading in range(maxLoading, 0, -1):

                        # check each processed LC-HRMS file if the used data processing parameters match
                        for res in results:
                            res.featurePairs=[]

                            if pwTextSet is not None:
                                # Log time used for bracketing
                                elapsed = (time.time() - start) / 60.
                                hours = ""
                                if elapsed >= 60.:
                                    hours = "%d hours " % (elapsed // 60)
                                mins = "%.2f mins" % (elapsed % 60.)
                                pwTextSet.put(Bunch(mes="text", val="<p align='right' >%s%s elapsed</p>\n\n\n\nProcessing \n  (Ionmode: %s, XCount: %s, Charge: %d) \n  File: %s" % (
                                                    hours, mins, ionMode, xCount, cLoading, res.fileName)))

                            for row in res.curs.execute(
                                    "SELECT id, fGroupID, assignedName, mz, loading, " \
                                           "ionMode, similarityString, PeakCenter, PeakCenterMin, PeakScale, " \
                                           "SNR, PeakArea, PeakAbundance, heteroIsotoplogues, assignedMZs, " \
                                           "comments, foundMatches " \
                                    "FROM chromPeaks WHERE ionMode='%s' AND similarityString='%s' and loading=%d"%(ionMode, xCount, cLoading)):
                                try:
                                    cp = ChromPeakFeature(id=row[0], PeakCenter=row[7], PeakCenterMin=row[8],
                                                          PeakScale=row[9], SNR=row[10], PeakArea=row[11],
                                                          mz=row[3], similarityString=str(row[6]),
                                                          loading=row[4], fGroupID=row[1],
                                                          ionMode=str(row[5]), file=str(res.fileName))

                                    assert cp.ionMode in ionModes.keys()
                                    res.featurePairs.append(cp)
                                except TypeError as err:
                                    print "  TypeError in file %s, id %s, (Ionmode: %s, XCount: %s, Charge: %d)"%(str(res.fileName), str(row[0]), ionMode, xCount, cLoading), err.message
                                except:
                                    print "  some general error in file %s, id %s, (Ionmode: %s, XCount: %s, Charge: %d)"%(str(res.fileName), str(row[0]), ionMode, xCount, cLoading)

                            totalChromPeaks = totalChromPeaks + len(res.featurePairs)



                        if pwTextSet is not None:
                            # Log time used for bracketing
                            elapsed = (time.time() - start) / 60.
                            hours = ""
                            if elapsed >= 60.:
                                hours = "%d hours " % (elapsed // 60)
                            mins = "%.2f mins" % (elapsed % 60.)
                            pwTextSet.put(Bunch(mes="text", val="<p align='right' >%s%s elapsed</p>\n\n\n\nClustering \n  (Ionmode: %s, XCount: %s, Charge: %d)" % (hours, mins, ionMode, xCount, cLoading)))


                            totalChromPeaks = totalChromPeaks + len(res.featurePairs)


                        for tracer in tracersDeltaMZ:
                            ## TODO this does not work correctly. improve
                            if pwValSet is not None: pwValSet.put(Bunch(mes="value", val=doneSoFar))
                            doneSoFar+=1

                            # get all results that match current bracketing criteria
                            chromPeaksAct = 0
                            allmz = []
                            for res in results:
                                chromPeaksAct = chromPeaksAct + len([chromPeak.mz for chromPeak in res.featurePairs if
                                                                     chromPeak.similarityString == xCount and chromPeak.loading == cLoading  and chromPeak.ionMode == ionMode])

                            allMZs=[]
                            for res in results:
                                allMZs.extend([chromPeak.mz for chromPeak in res.featurePairs if
                                                               chromPeak.similarityString == xCount and chromPeak.loading == cLoading and chromPeak.ionMode == ionMode])

                            # cluster all current results with HCA
                            allMZs=sorted(allMZs)
                            lastMZ=None
                            u=0
                            curAllMZs=[]

                            ## pre-separate detected feature pairs to improve speed of HCA
                            while len(allMZs)>0:
                                procCurSet=False
                                if u<len(allMZs) and (lastMZ is None or (allMZs[u]-lastMZ)<3*(groupSizePPM*allMZs[u]/1E6)):
                                    curAllMZs.append(allMZs[u])
                                    lastMZ=allMZs[u]
                                    u+=1
                                else:
                                    lastMZ=None
                                    allMZs=allMZs[u:]
                                    u=0
                                    procCurSet=True

                                if procCurSet or len(allMZs)==u:
                                    if len(curAllMZs)>0:
                                        hc = HierarchicalClustering(curAllMZs,
                                                                dist=lambda x, y: x.getValue() - y.getValue(),
                                                                val=lambda x: x, mean=lambda x, y: x / y,
                                                                add=lambda x, y: x + y)

                                        # cut HCA tree in subclusters
                                        for n in cutTreeSized(hc.getTree(), groupSizePPM):
                                            try:
                                                lowMZ=min([kid.getValue() for kid in n.getKids()])
                                                minMZ=max([kid.getValue() for kid in n.getKids()])

                                                maxMZInGroup = lowMZ

                                                partChromPeaks = {}
                                                partXICs = {}
                                                toDel = {}
                                                allChromPeaks=[]
                                                for res in results:
                                                    chromPeaks = []
                                                    toDel[res.filePath] = set()
                                                    mzs=[]
                                                    for i in range(len(res.featurePairs)):
                                                        chromPeak = res.featurePairs[i]
                                                        if chromPeak.similarityString == xCount and chromPeak.loading == cLoading  and chromPeak.ionMode == ionMode and chromPeak.mz <= minMZ:
                                                            mzs.append(chromPeak.mz)
                                                            toDel[res.filePath].add(i)
                                                            chromPeaks.append(chromPeak)
                                                            allChromPeaks.append(chromPeak)
                                                            maxMZInGroup = max(maxMZInGroup, chromPeak.mz)

                                                    xics = []
                                                    # get EICs of current subCluster
                                                    for chromPeak in chromPeaks:
                                                        eic, times, scanIds, mzs=mzXMLs[res.filePath].getEIC(mean(mzs), 5., filterLine=positiveScanEvent if chromPeak.ionMode=="+" else negativeScanEvent,
                                                                                    removeSingles=True, intThreshold=1000, useMS1=True, useMS2=False, startTime=0, endTime=1000000)
                                                        xics.append([len(eic), times, eic, 1])
                                                    partChromPeaks[res.filePath] = chromPeaks
                                                    partXICs[res.filePath] = xics

                                                    partChromPeaks[res.filePath] = chromPeaks
                                                eics = []
                                                peaks = []
                                                scantimes = []
                                                for k in partChromPeaks.keys():
                                                    for i in range(len(partChromPeaks[k])):
                                                        eics.append(partXICs[k][i][2])
                                                        peaks.append([partChromPeaks[k][i]])
                                                        scantimes.append(partXICs[k][i][1])

                                                # optional: align EICs; get bracketed chromatographic peaks
                                                aligned = xicAlign.alignXIC(eics, peaks, scantimes, align=align, maxTimeDiff=maxTimeDeviation, nPolynom=nPolynom)

                                                aligned = [(x[0][0], int(x[0][1])) for x in aligned]

                                                maxGroup = max(aligned, key=itemgetter(1))[1]
                                                minGroup = min(aligned, key=itemgetter(1))[1]

                                                groupedChromPeaks = []
                                                groupedChromPeaksAVGMz = []
                                                groupedChromPeaksAVGTimes = []
                                                groupedChromPeaksAVGPeakScale = []

                                                for i in range(maxGroup + 1):
                                                    groupedChromPeaks.append({})
                                                    groupedChromPeaksAVGMz.append([])
                                                    groupedChromPeaksAVGTimes.append([])
                                                    groupedChromPeaksAVGPeakScale.append([])

                                                # calculate average values (RT and mz) for feature pairs in the sub-subclusters
                                                j = 0
                                                for k in partChromPeaks.keys():
                                                    for i in range(len(partChromPeaks[k])):
                                                        if not (groupedChromPeaks[aligned[j][1]].has_key(k)):
                                                            groupedChromPeaks[aligned[j][1]][k] = []
                                                        groupedChromPeaks[aligned[j][1]][k].append((aligned[j], partChromPeaks[k][i]))
                                                        groupedChromPeaksAVGMz[aligned[j][1]].append(partChromPeaks[k][i].mz)
                                                        groupedChromPeaksAVGTimes[aligned[j][1]].append(partChromPeaks[k][i].PeakCenterMin)
                                                        groupedChromPeaksAVGPeakScale[aligned[j][1]].append(partChromPeaks[k][i].PeakScale)

                                                        j = j + 1

                                                assert (j == len(aligned))
                                                assert (len(groupedChromPeaks) == len(groupedChromPeaksAVGMz) == len(groupedChromPeaksAVGTimes))

                                                # write results to data matrix and SQLite DB
                                                for i in range(minGroup, maxGroup + 1):
                                                    if len(groupedChromPeaks[i]) > 0:
                                                        curNum = curNum + 1
                                                        avgmz = sum(groupedChromPeaksAVGMz[i]) / len(groupedChromPeaksAVGMz[i])
                                                        avgtime = sum(groupedChromPeaksAVGTimes[i]) / len(groupedChromPeaksAVGTimes[i])
                                                        avgPeakScale = sum(groupedChromPeaksAVGPeakScale[i]) / len(groupedChromPeaksAVGPeakScale[i])

                                                        f.write("\t".join([str(a) for a in [curNum, "", avgmz, "%f - %f"%(min(groupedChromPeaksAVGMz[i]), max(groupedChromPeaksAVGMz[i])),
                                                                                            avgtime , "%.2f - %.2f"%(min(groupedChromPeaksAVGTimes[i]), max(groupedChromPeaksAVGTimes[i])),
                                                                                            avgPeakScale, xCount, cLoading, scanEvent, ionMode,
                                                                                            "", "", "", ""]]))

                                                        SQLInsert(resDB.curs, "GroupResults", id=curNum, OGroup=curNum, mz=avgmz, rt=avgtime, similarityString=xCount, charge=cLoading, ionisationMode=ionMode, scanEvent=scanEvent)

                                                        doublePeak = 0
                                                        for j in range(len(results)):
                                                            res = results[j].filePath
                                                            f.write("\t")
                                                            if groupedChromPeaks[i].has_key(res) and len(groupedChromPeaks[i][res]) > 0:

                                                                for peak in groupedChromPeaks[i][res]:
                                                                    SQLInsert(resDB.curs, "FoundFeaturePairs", file=peak[1].file, featurePairID=peak[1].id,
                                                                              resID=curNum,
                                                                              PeakScale=peak[1].PeakScale,
                                                                              Area=peak[1].PeakArea, featureType="foundPattern")

                                                                doublePeak = doublePeak + len(groupedChromPeaks[i][res]) -1

                                                                ## Peak area
                                                                f.write(";".join([str(peak[1].PeakArea) for peak in groupedChromPeaks[i][res]]))
                                                        f.write("\t%d" % doublePeak)
                                                        f.write("\n")


                                            except Exception as e:
                                                import traceback
                                                traceback.print_exc()
                                                logging.error(str(traceback))

                                                logging.error("Error", str(e))
                                            grp = grp + 1

                                            for res in results:
                                                if toDel.has_key(res.filePath):
                                                    toDel[res.filePath] = [a for a in toDel[res.filePath]]
                                                    toDel[res.filePath].sort()
                                                    toDel[res.filePath].reverse()
                                                    for i in toDel[res.filePath]:
                                                        res.featurePairs.pop(i)
                                                        chromPeaksAct = chromPeaksAct - 1
                                                        totalChromPeaksProc = totalChromPeaksProc + 1

                                            if pwTextSet is not None:

                                                # Log time used for bracketing
                                                elapsed = (time.time() - start) / 60.
                                                hours = ""
                                                if elapsed >= 60.:
                                                    hours = "%d hours " % (elapsed // 60)
                                                mins = "%.2f mins" % (elapsed % 60.)
                                                pwTextSet.put(Bunch(mes="text",val="<p align='right' >%s%s elapsed</p>\n\n\n\nBracketing results\n%d feature pairs (%d individual pairs processed).. (Ionmode: %s, XCount: %s, Charge: %d)"%(hours, mins, curNum, totalChromPeaksProc, ionMode, xCount, cLoading)))
                                    curAllMZs=[]

            import uuid
            import platform
            import datetime

            identifier="%s_%s_%s"%(str(uuid.uuid1()), str(platform.node()), str(datetime.datetime.now()))

            f.write("## MetExtract II %s\n"%(Bunch(MetExtractVersion=meVersion, RVersion=rVersion, UUID_ext=identifier).dumpAsJSon().replace("\"", "'")))
            f.write("## Individual files processing parameters %s\n"%(generalProcessingParams.dumpAsJSon().replace("\"", "'")))
            processingParams=Bunch()

            xCountsFmt = []
            try:
                for xn in sorted(xCounts):
                    if xn-1 in xCounts and xn+1 in xCounts:
                        xCountsFmt.append("-")
                    else:
                        if xn-1 not in xCounts:
                            xCountsFmt.append(", ")
                        xCountsFmt.append(xn)
                xCountsFmt.pop(0)
                seen = set()
                seen_add = seen.add
                xCountsFmt = [x for x in xCountsFmt if not (x in seen or seen_add(x)) or x == ", "]
            except Exception:
                xCountsFmt=xCounts


            processingParams.FPBracketing_xCounts="".join([str(t) for t in xCountsFmt])
            processingParams.FPBracketing_groupSizePPM=groupSizePPM
            processingParams.FPBracketing_positiveScanEvent=positiveScanEvent
            processingParams.FPBracketing_negativeScanEvent=negativeScanEvent
            processingParams.FPBracketing_maxTimeDeviation=maxTimeDeviation
            processingParams.FPBracketing_maxLoading=maxLoading
            processingParams.FPBracketing_resultsFile=file
            processingParams.FPBracketing_align=align
            if align: processingParams.FPBracketing_nPolynom=nPolynom
            f.write("## Bracketing files processing parameters %s\n"%(processingParams.dumpAsJSon().replace("\"", "'")))

        exportAsFeatureML.convertMEMatrixToFeatureMLSepPolarities(file, postfix="_ab")

    except Exception as ex:
        import traceback
        traceback.print_exc()
        print(str(traceback))

        logging.error("Error during bracketing of files")

    finally:
        resDB.conn.commit()
        resDB.curs.close()
        resDB.conn.close()






























########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################














# store used configuration to DB file
def writeMetaboliteGroupingConfigToDB(curs, minConnectionsInFiles, minConnectionRate, groups):
    SQLInsert(curs, "config", key="MEGROUP_groups", value=str(groups))
    SQLInsert(curs, "config", key="MEGROUP_minConnectionsInFiles", value=str(minConnectionsInFiles))
    SQLInsert(curs, "config", key="MEGROUP_minConnectionRate", value=str(minConnectionRate))


def getBordersFor(curs, fID, file):
    for row in curs.execute("SELECT NBorderLeft, NBorderRight, LBorderLeft, LBorderRight FROM foundfeaturepairs WHERE file=? and resID=?", (file, fID)):
        return (float(row[0]), float(row[1]), float(row[2]), float(row[3]))
    return None

def getPeak(times, rt, borderleft, borderright):
    ## find best matching rt
    ind, tim=min(enumerate(times), key=lambda x:abs(x[1]-rt))
    minind=int(max(0, ind-int(borderleft)))
    maxind=int(min(len(times)-1, ind+borderright))

    return (minind, maxind)


def calculateMetaboliteGroups(file="./results.tsv", groups=[], eicPPM=10., maxAnnotationTimeWindow=0.05,
                              minConnectionsInFiles=1, minConnectionRate=0.4, minPeakCorrelation=0.85, useRatio=False,
                              runIdentificationInstance=None, pwMaxSet=None, pwValSet=None, pwTextSet=None, cpus=1, toFile=None):

    if toFile is None:
        toFile=file

    logging.info("Starting convoluting feature pairs..")

    resDB=Bunch(conn=None, curs=None)
    resDB.conn=connect(file+getDBSuffix())
    resDB.curs=resDB.conn.cursor()

    try:
        # Select only those groups that shall be used for metabolite group grouping
        ## load mzxml files for improved metabolite convolution
        useGroups=[]
        useGroupsForConfig=[]
        numFilesForConvolution=0
        for group in groups:
            if group.useForMetaboliteGrouping:
                useGroups.append(group)
                useGroupsForConfig.append(str(group.name)+":"+str(group.files))

                for fi in group.files:
                    numFilesForConvolution+=1

        if pwMaxSet is not None: pwMaxSet.put(Bunch(mes="max", val=numFilesForConvolution))
        if pwValSet is not None: pwValSet.put(Bunch(mes="value", val=0))

        # connect to results db
        writeMetaboliteGroupingConfigToDB(resDB.curs, minConnectionsInFiles, minConnectionRate, str(useGroupsForConfig).replace("'", "").replace("\"", ""))


        # read results table
        table = TableUtils.readFile(file, delim="\t")
        cols = [col.getName() for col in table.getColumns()]

        # check if necessary columns are present
        ## overall columns
        assert "Num" in cols
        assert "OGroup" in cols
        assert "Ion" in cols
        assert "Loss" in cols
        assert "M" in cols
        assert "doublePeak" in cols

        doublePeaks=0
        for row in table.getData(cols=["COUNT(*)"], where="doublePeak>0"):
            doublePeaks=int(row)

        if doublePeaks>0:
            logging.info("  found double peaks: %d"%doublePeaks)
        writeCols=[]
        if doublePeaks>0:
            TableUtils.saveFile(table, toFile.replace(".tsv", ".doublePeaks.tsv"), cols=writeCols, order="OGroup, MZ, Xn", where="doublePeak>0")

        table.executeSQL("DELETE FROM :table: WHERE doublePeak>0")
        ## file specific columns

        for group in useGroups:
            for f in group.files:
                fShort=f[f.rfind("/")+1 : f.rfind(".")]
                if fShort[0].isdigit():
                    fShort="_"+fShort

                #assert "%s_FID"%fShort in cols
                #assert "%s_GroupID"%fShort in cols


        # fetch all feature pairs from the results file
        nodes={}
        for fpNum, mz, rt, charge, scanEvent, ionMode in table.getData(cols=["Num", "MZ", "RT", "Charge", "ScanEvent", "IonMode"]):
            nodes[fpNum]=Bunch(fpNum=fpNum, mz=mz, rt=rt, charge=charge, scanEvent=scanEvent, ionMode=ionMode, oGroup=None, correlationToOtherFPs={})


        ## generate all correlations and sil ratios
        fileCorrelations = {}
        fileSILRatios = {}

        borders={}
        for row in resDB.curs.execute("SELECT PeakScale, file, resID FROM foundfeaturepairs"):
            fil=str(row[1])
            resID=int(row[2])
            if fil not in borders.keys():
                borders[fil]={}
            borders[fil][resID]=float(row[0])


        ## Iterate all files
        processedFiles=0
        procObjects=[]
        for group in useGroups:
            for fi in group.files:

                if fi not in fileCorrelations.keys():
                    #if pwTextSet is not None: pwTextSet.put(Bunch(mes="text", val="Convoluting feature pairs in file %s" % (fi)))
                    #if pwValSet is not None: pwValSet.put(Bunch(mes="value", val=processedFiles))

                    s=ConvoluteFPsInFile(eicPPM, fi, maxAnnotationTimeWindow, nodes, borders)
                    procObjects.append(s)

                    #processedFiles += 1

        from multiprocessing import Pool, Manager
        p = Pool(processes=cpus, maxtasksperchild=1)  # only in python >=2.7; experimental
        manager = Manager()
        queue = manager.Queue()

        start=time.time()
        # start the multiprocessing pool
        res = p.imap_unordered(processConvoluteFPsInFile, procObjects)

        # wait until all subprocesses have finished re-integrating their respective LC-HRMS data file
        loop = True
        freeSlots = range(min(len(procObjects), cpus))
        assignedThreads = {}
        while loop:
            completed = res._index
            if completed == len(procObjects):
                loop = False
            else:

                if pwValSet is not None: pwValSet.put(Bunch(mes="value", val=completed))

                elapsed = (time.time() - start) / 60.
                hours = ""
                if elapsed >= 60.:
                    hours = "%d hours " % (elapsed // 60)
                mins = "%.2f mins" % (elapsed % 60.)

                if pwTextSet is not None: pwTextSet.put(Bunch(mes="text", val="<p align='right' >%s%s elapsed</p>\n\n%d / %d files done (%d parallel)" % (hours, mins, completed , len(procObjects), cpus)))

                time.sleep(.5)

        p.close()
        p.terminate()
        p.join()

        ## fetch multiprocessing results
        for rei, re in enumerate(res):
            fileCorrelations[procObjects[rei].fi] = re[0]
            fileSILRatios[procObjects[rei].fi] = re[1]

        logging.info("Testing SIL ratios")
        ## test all feature pair pairs for co-elution and similar SIL ratio
        connections={}
        for fpNumA in nodes.keys():
            nodeA=nodes[fpNumA]
            if fpNumA not in connections.keys():
                connections[fpNumA]={}

            for fpNumB in nodes.keys():
                nodeB=nodes[fpNumB]
                if fpNumB not in connections.keys():
                    connections[fpNumB]={}

                if nodeA.fpNum!=nodeB.fpNum and nodeA.fpNum<nodeB.fpNum and abs(nodeB.rt-nodeA.rt)<=maxAnnotationTimeWindow:

                    ## test feature pair pair convolution in all samples
                    allCorrelations=[]
                    allSILRatios=[]
                    for group in useGroups:
                        for fi in group.files:

                            if fi in fileCorrelations.keys() and fpNumA in fileCorrelations[fi].keys() and fpNumB in fileCorrelations[fi][fpNumA].keys():
                                co = fileCorrelations[fi][fpNumA][fpNumB]
                                allCorrelations.append(co)
                                allSILRatios.append(1)
                            #if fi in fileSILRatios.keys() and fpNumA in fileSILRatios[fi].keys() and fpNumB in fileSILRatios[fi][fpNumA].keys():
                            #    silRatio = fileSILRatios[fi][fpNumA][fpNumB]
                            #    allSILRatios.append(silRatio)

                    if len(allCorrelations)>0 and len(allSILRatios)>0:
                        if sum([co>=minPeakCorrelation for co in allCorrelations])>=minConnectionRate*len(allCorrelations) and \
                                (not useRatio or sum([rat for rat in allSILRatios])>=minConnectionRate*len(allSILRatios)):
                            nodes[nodeA.fpNum].correlationToOtherFPs[nodeB.fpNum] = True
                            nodes[nodeB.fpNum].correlationToOtherFPs[nodeA.fpNum] = True

                        connections[fpNumA][fpNumB]=mean(allCorrelations)
                        connections[fpNumB][fpNumA]=mean(allCorrelations)
                    else:
                        connections[fpNumA][fpNumB] = 0
                        connections[fpNumB][fpNumA] = 0



        nodes2={}
        for fpNumA in nodes.keys():
            if fpNumA not in nodes2.keys():
                nodes2[fpNumA]=[]
            for fpNumB in nodes[fpNumA].correlationToOtherFPs.keys():
                nodes2[fpNumA].append(fpNumB)


        for k in nodes2.keys():
            uniq = []
            for u in nodes2[k]:
                if u not in uniq:
                    uniq.append(u)
            nodes2[k] = uniq

        # get subgraphs from the feature pair graph. Each subgraph represents one convoluted
        # feature group
        tGroups = getSubGraphs(nodes2)



        def splitGroupWithHCA(tGroup, correlations):
            # if 1 or two feature pairs are in a group, automatically use them (no further splitting possible)
            if len(tGroup) <= 2:
                return [tGroup]

            # construct 2-dimensional data matrix
            data = []
            for tG1 in tGroup:
                t = []
                for tG2 in tGroup:
                    if tG1 in correlations.keys() and tG2 in correlations[tG1].keys():
                        t.append(correlations[tG1][tG2])
                    elif tG1 == tG2:
                        t.append(1.)
                    else:
                        t.append(0.)
                data.append(t)

            # calculate HCA tree using the correlations between all feature pairs
            hc = HCA_general.HCA_generic()
            tree = hc.generateTree(objs=data, ids=tGroup)

            # split the HCA tree in sub-clusters
            def checkSubCluster(tree, hca, corrThreshold, cutOffMinRatio):
                if isinstance(tree, HCA_general.HCALeaf):
                    return False
                elif isinstance(tree, HCA_general.HCAComposite):
                    corrsLKid = hca.link(tree.getLeftKid())
                    corrs = hca.link(tree)
                    corrsRKid = hca.link(tree.getRightKid())

                    aboveThresholdLKid = sum([corr > corrThreshold for corr in corrsLKid])
                    aboveThreshold = sum([corr > corrThreshold for corr in corrs])
                    aboveThresholdRKid = sum([corr > corrThreshold for corr in corrsRKid])

                if (aboveThreshold * 1. / len(corrs)) >= cutOffMinRatio:
                    return False
                else:
                    return True

            subClusts = hc.splitTreeWithCallback(tree,CallBackMethod(_target=checkSubCluster, corrThreshold=minPeakCorrelation,cutOffMinRatio=minConnectionRate).getRunMethod(),recursive=False)

            # convert the subclusters into arrays of feature pairs belonging to the same metabolite
            return [[leaf.getID() for leaf in subClust.getLeaves()] for subClust in subClusts]


        groups = []
        done = 0
        for tGroup in tGroups:
            print "Group", done, "with elements", len(tGroup), "of", len(tGroups)
            cGroups = [tGroup]
            while len(cGroups) > 0:
                gGroup = cGroups.pop(0)

                sGroups = splitGroupWithHCA(gGroup, connections)

                if len(sGroups) == 1:
                    groups.append(sGroups[0])
                else:
                    cGroups.extend(sGroups)

            done = done + 1


        # Separate feature pair clusters; softer
        curGroup=1
        for tGroup in groups:
            if len(tGroup)==1:
                # if only one feature pair is in the group, do nothing
                table.setData(cols=["OGroup"], vals=[curGroup], where="Num=%d"%tGroup[0])
                resDB.curs.execute("UPDATE GroupResults SET OGroup=%d WHERE id = %d"%(curGroup, tGroup[0]))
            else:
                table.setData(cols=["OGroup"], vals=[curGroup], where="Num IN (%s)"%(",".join([str(t) for t in tGroup])))
                resDB.curs.execute("UPDATE GroupResults SET OGroup=%d WHERE id IN (%s)"%(curGroup, ",".join([str(t) for t in tGroup])))

            curGroup+=1

        logging.info("Annotating ions in feature groups")
        # Annotate groups with common adducts and in-source fragments
        if runIdentificationInstance is not None:
            groups=defaultdict(list)
            for row in table.getData(cols=["Num", "OGroup", "MZ", "IonMode", "Charge", "Ion", "Loss", "M"]):
                num, ogrp, mz, ionMode, loading, adducts, fDesc, ms=row
                groups[ogrp].append(ChromPeakFeature(id=num, fGroupID=ogrp, mz=mz, ionMode=ionMode, loading=loading, adducts=[], heteroAtomsFeaturePairs=[], fDesc=[]))

            for ogrp in groups.keys():
                chromPeaks={}
                for fp in groups[ogrp]:
                    chromPeaks[fp.id]=fp

                runIdentificationInstance.annotateChromPeaks(chromPeaks.keys(), chromPeaks)

            for ogrp in groups.keys():
                for fp in groups[ogrp]:
                    table.setData(cols=["Ion", "Loss", "M"], vals=[",".join([str(a) for a in fp.adducts]),
                                                              ",".join([str(a) for a in fp.fDesc]),
                                                              ",".join([str(a) for a in fp.Ms])],
                                  where="Num = %d"%(fp.id))
                    resDB.curs.execute("UPDATE GroupResults SET Ion='%s', Loss='%s', M='%s' WHERE id = %d"%(",".join([str(a) for a in fp.adducts]), ",".join([str(a) for a in fp.fDesc]), ",".join([str(a) for a in fp.Ms]), fp.id))

        resDB.curs.execute("DELETE FROM GroupResults WHERE id NOT IN (%s)"%(",".join([str(num) for num in table.getData(cols=["Num"])])))


        ## reassign feature group ids
        resDB.curs.execute("UPDATE GroupResults SET OGroup='X'||OGroup")
        table.executeSQL("UPDATE :table: SET OGroup='X'||OGroup")
        oGrps=[]
        curGrp=1
        curFP=1
        for row in resDB.curs.execute("SELECT OGroup, AVG(rt) FROM GroupResults GROUP BY OGroup ORDER BY AVG(rt)"):
            oGrps.append(str(row[0]))
        for tgrp in oGrps:
            resDB.curs.execute("UPDATE GroupResults SET OGroup='%d' WHERE OGroup='%s'"%(curGrp, tgrp))
            table.setData(cols=["OGroup"], vals=[curGrp], where="OGroup='%s'"%(tgrp))
            curGrp=curGrp+1




        processingParams=Bunch()
        processingParams.MEConvoluting_groups=str(useGroupsForConfig).replace("'","").replace("\"","")
        processingParams.MEConvoluting_connThreshold=minConnectionRate
        processingParams.MEConvoluting_minPeakCorrelation=minPeakCorrelation
        table.addComment("## Convolution FPs processing parameters %s"%(processingParams.dumpAsJSon().replace("\"", "'")))

        writeCols=[]
        for column in table.getColumns():
            if not column.name.endswith("_CorrelatesTo"):
                writeCols.append(column.name)

        TableUtils.saveFile(table, toFile, cols=writeCols, order="OGroup, MZ, Xn")

        exportAsFeatureML.convertMEMatrixToFeatureMLSepPolarities(file, postfix="_ac")

    except Exception as ex:
        import traceback
        traceback.print_exc()

        raise ex

    finally:
        resDB.conn.commit()
        resDB.curs.close()
        resDB.conn.close()




class ConvoluteFPsInFile:
    def __init__(self, eicPPM, fi, maxAnnotationTimeWindow, nodes, borders):
        self.eicPPM=eicPPM
        self.fi=fi
        self.maxAnnotationTimeWindow=maxAnnotationTimeWindow
        self.nodes=nodes
        self.borders=borders

        self.fileCorrelations=None
        self.fileSILRatios=None

    def getConvolutionInFile(self):

        startProc=time.time()

        eicPPM=self.eicPPM
        fi=self.fi
        maxAnnotationTimeWindow=self.maxAnnotationTimeWindow
        nodes=self.nodes
        borders=self.borders

        fileCorrelations = {}
        fileSILRatios = {}

        logging.info("  Convoluting feature pairs in file %s" % (fi))
        b = fi.replace("\\", "/")
        fiName = b[(b.rfind("/") + 1):b.rfind(".mzXML")]
        mzXML = Chromatogram()
        mzXML.parse_file(fi)
        for fpNumA in nodes.keys():
            nodeA = nodes[fpNumA]
            if fpNumA not in fileCorrelations.keys():
                fileCorrelations[fpNumA] = {}
                fileSILRatios[fpNumA] = {}

            for fpNumB in nodes.keys():
                nodeB = nodes[fpNumB]
                if fpNumB not in fileCorrelations.keys():
                    fileCorrelations[fpNumB] = {}
                    fileSILRatios[fpNumB] = {}

                if nodeA.fpNum != nodeB.fpNum and nodeA.fpNum < nodeB.fpNum and abs(nodeB.rt - nodeA.rt) <= maxAnnotationTimeWindow:

                    ra=None
                    if fiName in borders.keys() and nodeA.fpNum in borders[fiName].keys():
                        ra=borders[fiName][nodeA.fpNum]
                    rb=None
                    if fiName in borders.keys() and nodeB.fpNum in borders[fiName].keys():
                        rb=borders[fiName][nodeB.fpNum]

                    if ra != None and rb != None:

                        if nodeA.scanEvent in mzXML.getFilterLines(includeMS1=True, includeMS2=False,
                                                                   includePosPolarity=True, includeNegPolarity=True) and \
                                        nodeB.scanEvent in mzXML.getFilterLines(includeMS1=True, includeMS2=False,
                                                                                includePosPolarity=True,
                                                                                includeNegPolarity=True):

                            meanRT = mean([nodeA.rt, nodeB.rt])

                            eicAL, timesA, scanIdsA, mzs = mzXML.getEIC(nodeA.mz, ppm=eicPPM, filterLine=nodeA.scanEvent,
                                                                      removeSingles=True, intThreshold=0, useMS1=True,
                                                                      useMS2=False, startTime=meanRT * 60 - 120,
                                                                      endTime=meanRT * 60 + 120)
                            eicBL, timesB, scanIdsB, mzs = mzXML.getEIC(nodeB.mz, ppm=eicPPM, filterLine=nodeB.scanEvent,
                                                                      removeSingles=True, intThreshold=0, useMS1=True,
                                                                      useMS2=False, startTime=meanRT * 60 - 120,
                                                                      endTime=meanRT * 60 + 120)

                            eicA, timesA, scanIdsA, mzs = mzXML.getEIC(nodeA.mz, ppm=eicPPM, filterLine=nodeA.scanEvent,
                                                                     removeSingles=True, intThreshold=0, useMS1=True,
                                                                     useMS2=False, startTime=meanRT * 60 - 120,
                                                                     endTime=meanRT * 60 + 120)
                            eicB, timesB, scanIdsB, mzs = mzXML.getEIC(nodeB.mz, ppm=eicPPM, filterLine=nodeB.scanEvent,
                                                                     removeSingles=True, intThreshold=0, useMS1=True,
                                                                     useMS2=False, startTime=meanRT * 60 - 120,
                                                                     endTime=meanRT * 60 + 120)



                            timesMin = [t for t in timesA]
                            timesMin.extend([t for t in timesB])
                            timesMin=sorted(timesMin)

                            from utils import mapArrayToRefTimes

                            eicAL=mapArrayToRefTimes(eicAL, timesA, timesMin)
                            eicBL=mapArrayToRefTimes(eicBL, timesB, timesMin)
                            eicA=mapArrayToRefTimes(eicA, timesA, timesMin)
                            eicB=mapArrayToRefTimes(eicB, timesB, timesMin)

                            timesMin=[t/60. for t in timesMin]


                            ## A) Test correlation of different feature pairs
                            try:
                                lI, rI = getPeak(timesMin, meanRT, min(ra, rb)*2, min(ra, rb)*2)

                                eicAC = eicAL[lI:rI]
                                eicBC = eicBL[lI:rI]

                                co = corr(eicAC, eicBC)

                                if not (isnan(co)) and co != None:
                                    fileCorrelations[fpNumA][fpNumB] = co
                            except Exception as err:
                                logging.error(
                                    "  Error during convolution of feature pairs (Peak-correlation, Nums: %s and %s, message: %s).." % (
                                    fpNumA, fpNumB, err.message))

                            ## B) Test similarity of native:labeled ratio
                            try:
                                lISIL, rISIL = getPeak(timesMin, nodeA.rt, min(ra, ra) * .8, min(ra, ra) * .8)

                                eicANCSIL = eicA[lISIL:rISIL]
                                eicALCSIL = eicAL[lISIL:rISIL]

                                folds = [eicANCSIL[i] / eicALCSIL[i] for i in range(len(eicANCSIL)) if eicALCSIL[i] > 0]
                                ma = mean(folds)
                                sa = sd(folds)

                                lISIL, rISIL = getPeak(timesMin, nodeB.rt, min(ra, ra) * .8, min(ra, ra) * .8)

                                eicBNCSIL = eicB[lISIL:rISIL]
                                eicBLCSIL = eicBL[lISIL:rISIL]

                                folds = [eicBNCSIL[i] / eicBLCSIL[i] for i in range(len(eicBNCSIL)) if eicBLCSIL[i] > 0]
                                mb = mean(folds)
                                sb = sd(folds)

                                if ma != None and mb != None and sa != None and sb != None:
                                    silRatioFold = (max([ma, mb]) / min([ma, mb]))
                                    #print ma, mb, sa, sb, max([ma, mb]), min([ma, mb])
                                    fileSILRatios[fpNumA][fpNumB] = silRatioFold <= 1 + max(0.5, 3 * mean([sa, sb]))

                                    #print fiName, nodeA.mz, nodeB.mz, meanRT, silRatioFold, co
                            except Exception as err:
                                logging.error(
                                    "  Error during convolution of feature pairs (SIL-ratio, Nums: %s and %s, message: %s).." % (
                                    fpNumA, fpNumB, err.message))

        logging.info("  Finished convoluting feature pairs in file %s (%.1f minutes)" % (fi, (time.time()-startProc)/60.))
        return (fileCorrelations, fileSILRatios)

def processConvoluteFPsInFile(cif):
    se=cif.getConvolutionInFile()
    return se




