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

import base64
import logging
import os
import platform
import pprint
import re
from copy import copy
from math import floor
from pickle import loads, dumps
from sqlite3 import *
from collections import OrderedDict

import numpy as np

import HCA_general
from utils import USEGRADIENTDESCENDPEAKPICKING
from utils import getNormRatio, CallBackMethod

pp = pprint.PrettyPrinter(indent=1)

from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics.charts.lineplots import ScatterPlot
from reportlab.lib import pagesizes
from reportlab.lib import colors
from reportlab.lib.colors import Color
from reportlab.lib.units import mm
from reportlab.platypus import Paragraph, Table
from reportlab.platypus.flowables import Flowable
from reportlab.pdfgen import canvas
from reportlab.graphics import renderPDF
from reportlab.lib.styles import getSampleStyleSheet

# make behavior os independent
def get_main_dir():
    from utils import get_main_dir
    return os.path.join(get_main_dir(), '')

# helper function for ReportLab
def coord(x, y):
    return x, y

# helper function for ReportLab
def noLabel(a):
    return ""

# helper class for ReportLab
class TTR(Flowable):  #TableTextRotate
    '''Rotates a tex in a table cell.'''

    def __init__(self, text):
        Flowable.__init__(self)
        self.text = text

    def draw(self):
        canvas = self.canv
        canvas.rotate(90)
        canvas.drawString(0, -1, self.text)


from Chromatogram import Chromatogram
from runIdentification_matchPartners import matchPartners
from formulaTools import formulaTools, getIsotopeMass
from utils import getAtomAdd, mean, weightedMean, sd, weightedSd, smoothDataSeries, SQLInsert, SQLSelectAsObject
from MZHCA import HierarchicalClustering, cutTreeSized
from chromPeakPicking.MassSpecWavelet import MassSpecWavelet
import Baseline
from utils import corr, getSubGraphs, ChromPeakFeature, Bunch, natSort
from SGR import SGRGenerator
from mePyGuis.TracerEdit import ConfiguredTracer
import exportAsFeatureML

from matchIsotopologPatternRules import *


def getDBSuffix():
    return ".identified.sqlite"


# returns the abbreviation for a given element
# e.g. element="Carbon" --> return="C"
def getShortName(element):
    fT = formulaTools()
    for i in fT.elemDetails:
        d = fT.elemDetails[i]
        if d[0] == element:
            return d[1]
    return ""


# counts how ofter an entry i is listed in the vector x and returns it as a dictionary
# e.g. x=[1,2,1,3,4,4,5,4] --> return={1:2, 2:1, 3:1, 4:3, 5:1}
def countEntries(x):
    ret = {}
    for i in x:
        if not (ret.has_key(i)):
            ret[i] = 0
        ret[i] = ret[i] + 1
    return ret


peakAbundanceUseSignals = 5
peakAbundanceUseSignalsSides = int((peakAbundanceUseSignals - 1) / 2)
# This class is used as a Command for each LC-HRMS file and is called by the multiprocessing module in MExtract.py
# during the processing of each individual LC-HRMS data files
class RunIdentification:

    # Constructor of the class, which stores the processing parameters
    # writeMZXML: 0001: 12C  0010: 12C-Iso  0100: 13C-Iso  1000: 13C
    def __init__(self, file, exOperator="", exExperimentID="", exComments="", exExperimentName="",
                 writePDF=False, writeFeatureML=False, writeTSV=False, writeMZXML=9,
                 metabolisationExperiment=False,
                 labellingisotopeA='12C', labellingisotopeB='13C', useCIsotopePatternValidation=0, xOffset=1.00335,
                 configuredTracer=None, intensityThreshold=0, intensityCutoff=0, maxLoading=1, xCounts="", ppm=2.,
                 purityN=0.99, purityL=0.99, minSpectraCount=1, clustPPM=8.,
                 chromPeakPPM=5., snrTh=1., scales=[1, 35], peakCenterError=5, peakScaleError=3, minPeakCorr=0.85,
                 checkPeaksRatio=False,
                 calcIsoRatioNative=1, calcIsoRatioLabelled=-1, calcIsoRatioMoiety=1,
                 startTime=2, stopTime=37, positiveScanEvent="None", negativeScanEvent="None",
                 eicSmoothingWindow="None", eicSmoothingWindowSize=0, eicSmoothingPolynom=0, artificialMPshift_start=0, artificialMPshift_stop=0,
                 correctCCount=True, minCorrelation=0.85, minCorrelationConnections=0.4, hAIntensityError=5., hAMinScans=3, adducts=[],
                 elements=[], heteroAtoms=[], simplifyInSourceFragments=True, chromPeakFile=None, lock=None, queue=None, pID=-1, rVersion="NA",
                 meVersion="NA", rules=""):

        self.labellingIsotopeA=labellingisotopeA
        self.labellingIsotopeB=labellingisotopeB
        ma, ea = getIsotopeMass(labellingisotopeA)
        mb, eb = getIsotopeMass(labellingisotopeB)

        if (not metabolisationExperiment and len(ea) > 0 and ea == eb) or metabolisationExperiment or (
                        ea == "Hydrogen" and eb == "Deuterium"):
            pass
        else:
            self.printMessage("Labelling specifications are invalid..", type="warning")
            raise Exception("Labelling specifications are invalid..")

        if positiveScanEvent == "None" and negativeScanEvent == "None":
            self.printMessage("No scan event was specified..", type="warning")
            raise Exception("No scan event was specified..")

        self.file = file
        self.writePDF = writePDF
        self.writeFeatureML = writeFeatureML
        self.writeTSV = writeTSV
        self.writeMZXML = writeMZXML

        #E. Experiment
        self.experimentOperator=exOperator
        self.experimentID=exExperimentID
        self.experimentComments=exComments
        self.experimentName=exExperimentName

        #0. General
        self.startTime = startTime
        self.stopTime = stopTime
        self.positiveScanEvent = positiveScanEvent
        self.negativeScanEvent = negativeScanEvent

        self.metabolisationExperiment = metabolisationExperiment

        self.labellingElement = getShortName(ea)
        self.isotopeA = ma
        self.isotopeB = mb
        self.useCIsotopePatternValidation = useCIsotopePatternValidation

        self.configuredTracer = configuredTracer
        if not self.metabolisationExperiment:
            self.configuredTracer = ConfiguredTracer(name="FML", id=0)

        #1. Mass picking
        self.intensityThreshold = intensityThreshold
        self.intensityCutoff = intensityCutoff
        self.maxLoading = maxLoading

        self.xCounts = []
        self.xCountsString=xCounts
        a=xCounts.replace(" ", "").replace(";",",").split(",")
        for j in a:
            if "-" in j:
                self.xCounts.extend(range(int(j[0:j.find("-")]), int(j[j.find("-")+1:len(j)])+1))
            elif j!="":
                self.xCounts.append(int(j))
        self.xCounts=sorted(list(set(self.xCounts)))

        self.xOffset = xOffset
        self.ppm = ppm
        self.purityN = purityN
        self.purityL = purityL

        #2. Results clustering
        self.minSpectraCount = max(1, minSpectraCount)
        self.clustPPM = clustPPM

        #3. Peak detection
        self.chromPeakPPM = chromPeakPPM

        if eicSmoothingWindow == None or eicSmoothingWindow == "" or eicSmoothingWindow == "none":
            eicSmoothingWindow = "None"
        self.eicSmoothingWindow = eicSmoothingWindow
        self.eicSmoothingWindowSize = eicSmoothingWindowSize
        self.eicSmoothingPolynom = eicSmoothingPolynom
        self.artificialMPshift_start = artificialMPshift_start
        self.artificialMPshift_stop  = artificialMPshift_stop
        self.snrTh = snrTh
        self.scales = scales
        self.peakCenterError = peakCenterError
        self.peakScaleError = peakScaleError
        self.minPeakCorr = minPeakCorr
        self.checkPeaksRatio=checkPeaksRatio

        self.calcIsoRatioNative=calcIsoRatioNative
        self.calcIsoRatioLabelled=calcIsoRatioLabelled
        self.calcIsoRatioMoiety=calcIsoRatioMoiety

        self.chromPeakFile = chromPeakFile
        if self.chromPeakFile is None:
            cpf = get_main_dir()+ "./chromPeakPicking/MassSpecWaveletIdentification.r"
            self.chromPeakFile = cpf

        self.performCorrectCCount = correctCCount

        self.minCorrelation = minCorrelation
        self.minCorrelationConnections = minCorrelationConnections

        self.hAIntensityError = hAIntensityError
        self.hAMinScans = hAMinScans
        self.adducts = adducts
        self.elements = {}
        for elem in elements:
            self.elements[elem.name] = elem
        self.heteroAtoms = {}
        for heteroAtom in heteroAtoms:
            self.heteroAtoms[heteroAtom.name] = heteroAtom
        self.simplifyInSourceFragments=simplifyInSourceFragments

        #Rules
        self.rules=rules

        #System
        self.lock = lock
        self.queue = queue
        self.pID = pID
        self.rVersion = rVersion
        self.meVersion = meVersion

    # Thread safe printing function
    def printMessage(self, message, type="info"):
        if self.lock is not None:
            self.lock.acquire()
            if type.lower()=="info":
                logging.info("   %d: %s" % (self.pID, message))
            elif type.lower()=="warning":
                logging.warning("   %d: %s" % (self.pID, message))
            elif type.lower()=="error":
                logging.error("   %d: %s" % (self.pID, message))
            else:
                logging.debug("   %d: %s" % (self.pID, message))
            self.lock.release()

    # helper function used to update the status of the current processing in the Process Dialog
    def postMessageToProgressWrapper(self, mes, val=""):
        if self.pID != -1 and self.queue is not None:
            if mes.lower() == "text":
                self.queue.put(Bunch(pid=self.pID, mes="text", val="%d: %s" % (self.pID, val)))
            elif mes == "value" or mes == "max":
                self.queue.put(Bunch(pid=self.pID, mes=mes, val=val))
            elif mes == "start" or mes == "end" or mes == "failed":
                self.queue.put(Bunch(pid=self.pID, mes=mes))

    def getMostLikelyHeteroIsotope(self, foundIsotopes):

        if len(foundIsotopes) == 0:
            return

        for iso in foundIsotopes:
            isoD = foundIsotopes[iso]
            useCn = -1
            useCnRatio = 0.
            useCnRatioTheo = 0.

            for cn in isoD:
                cnF = isoD[cn]

                cnRatio = mean([x[1] for x in cnF])
                cnRatioTheo = mean([x[2] for x in cnF])
                cnScans = len(cnF)

                if useCn == -1 or (abs(cnRatioTheo - cnRatio) < abs(useCnRatio - useCnRatioTheo)):
                    useCn = cn
                    useCnRatio = cnRatio
                    useCnRatioTheo = cnRatioTheo

            foundIsotopes[iso] = {useCn: isoD[useCn]}

    # store configuration used for processing the LC-HRMS file into the database
    def writeConfigurationToDB(self, conn, curs):
        curs.execute("DROP TABLE IF EXISTS config")
        curs.execute("CREATE TABLE config (id INTEGER PRIMARY KEY AUTOINCREMENT, key TEXT, value TEXT)")

        curs.execute("DROP TABLE IF EXISTS tracerConfiguration")
        curs.execute("create table tracerConfiguration(id INTEGER PRIMARY KEY, name TEXT, elementCount INTEGER, natural TEXT, labelling TEXT, "
                     "deltaMZ REAL, purityN REAL, purityL REAL, amountN REAL, amountL REAL, monoisotopicRatio REAL, lowerError REAL, "
                     "higherError REAL, tracertype TEXT)")

        curs.execute("DROP TABLE IF EXISTS MZs")
        curs.execute("create table MZs(id INTEGER PRIMARY KEY, tracer INTEGER, mz REAL, similarityObject TEXT, scanid INTEGER, scantime REAL, "
                     "loading INTEGER, intensity FLOAT, ionMode TEXT, type TEXT, otherIsotopologs TEXT)")

        curs.execute("DROP TABLE IF EXISTS MZBins")
        curs.execute("CREATE TABLE MZBins(id INTEGER PRIMARY KEY, mz REAL, ionMode TEXT)")

        curs.execute("DROP TABLE IF EXISTS MZBinsKids")
        curs.execute("CREATE TABLE MZBinsKids(mzbinID INTEGER, mzID INTEGER)")

        curs.execute("DROP TABLE IF EXISTS chromPeaks")
        curs.execute("CREATE TABLE chromPeaks(id INTEGER PRIMARY KEY, fGroupID INTEGER, assignedName TEXT, mz FLOAT, loading INTEGER, "
                     "ionMode TEXT, similarityString TEXT, PeakCenter INTEGER, PeakCenterMin FLOAT, PeakScale FLOAT, SNR FLOAT, "
                     "PeakArea FLOAT, PeakAbundance FLOAT, heteroIsotoplogues TEXT, assignedMZs TEXT, comments TEXT, foundMatches TEXT)")

        curs.execute("DROP TABLE IF EXISTS allChromPeaks")
        curs.execute("CREATE TABLE allChromPeaks(id INTEGER PRIMARY KEY, fGroupID INTEGER, assignedName TEXT, mz FLOAT, loading INTEGER, "
                     "ionMode TEXT, similarityString TEXT, PeakCenter INTEGER, PeakCenterMin FLOAT, PeakScale FLOAT, SNR FLOAT, "
                     "PeakArea FLOAT, PeakAbundance FLOAT, heteroIsotoplogues TEXT, assignedMZs TEXT, comments TEXT, foundMatches TEXT)")

        curs.execute("DROP TABLE IF EXISTS stats")
        curs.execute("CREATE TABLE stats (key TEXT, value TEXT)")

        SQLInsert(curs, "config", key="MetExtractVersion", value=self.meVersion)
        SQLInsert(curs, "config", key="RVersion", value=self.rVersion)

        SQLInsert(curs, "config", key="ExperimentName", value=self.experimentName)
        SQLInsert(curs, "config", key="ExperimentOperator", value=self.experimentOperator)
        SQLInsert(curs, "config", key="ExperimentID", value=self.experimentID)
        SQLInsert(curs, "config", key="ExperimentComments", value=self.experimentComments)

        SQLInsert(curs, "config", key="labellingElement", value=self.labellingElement)
        SQLInsert(curs, "config", key="isotopeA", value=self.isotopeA)
        SQLInsert(curs, "config", key="isotopeB", value=self.isotopeB)
        SQLInsert(curs, "config", key="useCValidation", value=self.useCIsotopePatternValidation)
        SQLInsert(curs, "config", key="metabolisationExperiment", value=str(self.metabolisationExperiment))
        SQLInsert(curs, "config", key="configuredTracer", value=base64.b64encode(dumps(self.configuredTracer)))
        SQLInsert(curs, "config", key="startTime", value=self.startTime)
        SQLInsert(curs, "config", key="stopTime", value=self.stopTime)
        SQLInsert(curs, "config", key="positiveScanEvent", value=self.positiveScanEvent)
        SQLInsert(curs, "config", key="negativeScanEvent", value=self.negativeScanEvent)
        SQLInsert(curs, "config", key="intensityThreshold", value=self.intensityThreshold)
        SQLInsert(curs, "config", key="intensityCutoff", value=self.intensityCutoff)
        SQLInsert(curs, "config", key="maxLoading", value=self.maxLoading)
        SQLInsert(curs, "config", key="xCounts", value=self.xCountsString)
        SQLInsert(curs, "config", key="xOffset", value=self.xOffset)
        SQLInsert(curs, "config", key="ppm", value=self.ppm)
        SQLInsert(curs, "config", key="purityN", value=self.purityN)
        SQLInsert(curs, "config", key="purityL", value=self.purityL)
        SQLInsert(curs, "config", key="clustPPM", value=self.clustPPM)
        SQLInsert(curs, "config", key="chromPeakPPM", value=self.chromPeakPPM)
        SQLInsert(curs, "config", key="eicSmoothing", value=self.eicSmoothingWindow)
        SQLInsert(curs, "config", key="eicSmoothingWindowSize", value=self.eicSmoothingWindowSize)
        SQLInsert(curs, "config", key="eicSmoothingPolynom", value=self.eicSmoothingPolynom)
        SQLInsert(curs, "config", key="artificialMPshift_start", value=self.artificialMPshift_start)
        SQLInsert(curs, "config", key="artificialMPshift_stop", value=self.artificialMPshift_stop)
        SQLInsert(curs, "config", key="snrTh", value=self.snrTh)
        SQLInsert(curs, "config", key="scales", value=base64.b64encode(dumps(self.scales)))
        SQLInsert(curs, "config", key="peakAbundanceCriteria", value="Center +- %d signals (%d total)"%(peakAbundanceUseSignalsSides, peakAbundanceUseSignals))
        SQLInsert(curs, "config", key="peakCenterError", value=self.peakCenterError)
        SQLInsert(curs, "config", key="peakScaleError", value=self.peakScaleError)
        SQLInsert(curs, "config", key="minPeakCorr", value=self.minPeakCorr)
        SQLInsert(curs, "config", key="checkPeaksRatio", value=str(self.checkPeaksRatio))
        SQLInsert(curs, "config", key="calcIsoRatioNative", value=self.calcIsoRatioNative)
        SQLInsert(curs, "config", key="calcIsoRatioLabelled", value=self.calcIsoRatioLabelled)
        SQLInsert(curs, "config", key="calcIsoRatioMoiety", value=self.calcIsoRatioMoiety)
        SQLInsert(curs, "config", key="minSpectraCount", value=self.minSpectraCount)
        SQLInsert(curs, "config", key="configuredHeteroAtoms", value=base64.b64encode(dumps(self.heteroAtoms)))
        SQLInsert(curs, "config", key="haIntensityError", value=self.hAIntensityError)
        SQLInsert(curs, "config", key="haMinScans", value=self.hAMinScans)
        SQLInsert(curs, "config", key="minCorrelation", value=self.minCorrelation)
        SQLInsert(curs, "config", key="minCorrelationConnections", value=self.minCorrelationConnections)
        SQLInsert(curs, "config", key="adducts", value=base64.b64encode(dumps(self.adducts)))
        SQLInsert(curs, "config", key="elements", value=base64.b64encode(dumps(self.elements)))
        SQLInsert(curs, "config", key="simplifyInSourceFragments", value=str(self.simplifyInSourceFragments))

        import uuid
        import platform
        import datetime

        self.processingUUID="%s_%s_%s"%(str(uuid.uuid1()), str(platform.node()), str(datetime.datetime.now()))
        SQLInsert(curs, "config", key="processingUUID_ext", value=base64.b64encode(str(self.processingUUID)))

        if self.metabolisationExperiment:
            self.configuredTracer.id = 1
            SQLInsert(curs, "tracerConfiguration", id=tracer.id, name=tracer.name, elementCount=tracer.elementCount, natural=tracer.isotopeA, labelling=tracer.isotopeB,
                  deltaMZ=getIsotopeMass(tracer.isotopeB)[0] - getIsotopeMass(tracer.isotopeA)[0], purityN=tracer.enrichmentA, purityL=tracer.enrichmentB,
                  amountN=tracer.amountA, amountL=tracer.amountB, monoisotopicRatio=tracer.monoisotopicRatio, lowerError=tracer.maxRelNegBias,
                  higherError=tracer.maxRelPosBias, tracertype=tracer.tracerType)

        else:
            #ConfiguredTracer(name="Full metabolome labeling experiment", id=0)
            SQLInsert(curs, "tracerConfiguration", id=0, name="FLE")
        conn.commit()

    def parseMzXMLFile(self):
        mzxml = Chromatogram()
        mzxml.parse_file(self.file, intensityCutoff=self.intensityCutoff)
        return mzxml

    # creates a new PDF page which contains the used data processing parameters
    def writeSettingsToPDF(self, pdf):
        currentHeight = 800

        pdf.drawString(50, currentHeight, "Experiment name")
        pdf.drawString(240, currentHeight, self.experimentName)
        pdf.line(50, currentHeight - 2, 540, currentHeight - 2)
        currentHeight -= 20

        pdf.drawString(70, currentHeight, "Operator")
        pdf.drawString(240, currentHeight, self.experimentOperator)
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "ID")
        pdf.drawString(240, currentHeight, self.experimentID)
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Comments")
        currentHeight -= 15


        p = Paragraph(self.experimentComments, style=getSampleStyleSheet()["Normal"])
        w, h = p.wrap(450, 120)
        p.wrapOn(pdf, 450, 120)
        p.drawOn(pdf, 60, currentHeight - h + 10);

        pdf.showPage()

        currentHeight = 800

        pdf.drawString(50, currentHeight, "Parameters");
        pdf.line(50, currentHeight - 2, 540, currentHeight - 2);
        currentHeight -= 20

        if self.metabolisationExperiment:
            pdf.drawString(70, currentHeight, "Metabolisation Experiment")
            currentHeight -= 25

        else:
            pdf.drawString(50, currentHeight, "Full Metabolome Experiment");
            currentHeight -= 15
            pdf.drawString(50, currentHeight, "Labelling")

            pdf.drawString(240, currentHeight, "%s (%s %s)" % (self.labellingElement, self.isotopeA, self.isotopeB));
            currentHeight -= 15

            pdf.drawString(50, currentHeight, "Abundance %d%s" % (self.isotopeA, str(self.labellingElement)))
            pdf.drawString(240, currentHeight, "%.2f%%" % self.purityN);
            currentHeight -= 15

            pdf.drawString(70, currentHeight, "Abundance %d%s" % (self.isotopeB, str(self.labellingElement)))
            pdf.drawString(240, currentHeight, "%.2f%%" % self.purityL);
            currentHeight -= 35

        pdf.drawString(50, currentHeight, "M/Z Picking");
        pdf.line(50, currentHeight - 2, 540, currentHeight - 2);
        currentHeight -= 15
        pdf.drawString(70, currentHeight, "Scan event(s)");
        currentHeight -= 15

        if self.positiveScanEvent != "None":
            pdf.drawString(90, currentHeight, "Positive mode")
            pdf.drawString(240, currentHeight, "%s" % self.positiveScanEvent);
            currentHeight -= 15
        if self.negativeScanEvent != "None":
            pdf.drawString(90, currentHeight, "Negative mode")
            pdf.drawString(220, currentHeight, "%s" % self.negativeScanEvent);
            currentHeight -= 15

        pdf.drawString(70, currentHeight, "Intensity threshold")
        pdf.drawString(240, currentHeight, "%d%s" % (self.intensityThreshold,
                                                     " (Low abundance cutoff for isotopologues: %.0f)"%self.intensityThresholdIsotopologs if self.lowAbundanceIsotopeCutoff else ""));
        currentHeight -= 15
        pdf.drawString(70, currentHeight, "Intensity cutoff")
        pdf.drawString(240, currentHeight, "%d" % (self.intensityCutoff));
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Time")
        pdf.drawString(240, currentHeight, "%.1f-%.1f min" % (self.startTime, self.stopTime));
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Atom count")
        pdf.drawString(240, currentHeight, "%s" % (",".join(str(x) for x in self.xCountsString)));
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Max. charge")
        pdf.drawString(240, currentHeight, "%d" % self.maxLoading);
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Max. mass deviation ")
        pdf.drawString(240, currentHeight, "%.2f" % self.ppm);
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Isotopic pattern count")
        pdf.drawString(240, currentHeight, "Native: %d Labelled: %d" % (self.isotopicPatternCountLeft, self.isotopicPatternCountRight));
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Intensity abundance error")
        pdf.drawString(240, currentHeight, "Native: %.1f%% Labelled: %.1f%%" % (self.intensityErrorN * 100., self.intensityErrorL * 100.));
        currentHeight -= 35

        pdf.drawString(50, currentHeight, "Post Processing");
        pdf.line(50, currentHeight - 2, 540, currentHeight - 2);
        currentHeight -= 20

        pdf.drawString(70, currentHeight, "Clustering PPM")
        pdf.drawString(240, currentHeight, "%.2f" % self.clustPPM);
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "EIC ppm")
        pdf.drawString(240, currentHeight, "%.2f" % self.chromPeakPPM);
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "EIC Smoothing Window")
        pdf.drawString(240, currentHeight, self.eicSmoothingWindow)
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "M' artificial shift")
        pdf.drawString(240, currentHeight, "%d - %d"%(self.artificialMPshift_start, self.artificialMPshift_stop))
        currentHeight -=15

        if self.eicSmoothingWindow.lower() != "none":
            pdf.drawString(70, currentHeight, "EIC Smoothing Window Size")
            pdf.drawString(240, currentHeight, str(self.eicSmoothingWindowSize))
            currentHeight -= 15

        if self.eicSmoothingWindow.lower() == "SavitzkyGolay":
            pdf.drawString(70, currentHeight, "EIC Smoothing Polynom")
            pdf.drawString(240, currentHeight, str(self.eicSmoothingPolynom))
            currentHeight -= 15

        pdf.drawString(70, currentHeight, "Scales")
        pdf.drawString(240, currentHeight, "%d - %d" % (self.scales[0], self.scales[1]));
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Peak matching")
        #pdf.drawString(240, currentHeight, "Center: %d Scales: %d Min. Corr: %.2f" % (
        #    self.peakCenterError, self.peakScaleError, self.minPeakCorr));
        pdf.drawString(240, currentHeight, "Center: %d, Min. Corr: %.2f" % (
            self.peakCenterError, self.minPeakCorr));
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Min. Corr:")
        pdf.drawString(240, currentHeight, "%.2f" % self.minPeakCorr);
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Min. scans")
        pdf.drawString(240, currentHeight, "%d" % self.minSpectraCount);
        currentHeight -= 35

        pdf.drawString(50, currentHeight, "Hetero atom annotation");
        pdf.line(50, currentHeight - 2, 540, currentHeight - 2);
        currentHeight -= 20

        pdf.drawString(70, currentHeight, "Max. intensity error")
        pdf.drawString(240, currentHeight, "%.1f%%" % (self.hAIntensityError * 100.));
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Min. scans")
        pdf.drawString(240, currentHeight, "%d" % self.hAMinScans);
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Hetero atom isotopes")
        p = Paragraph(str(", ".join(["%s (m: %.4f, rel. ab. %.1f%%, min: %d, max: %d)" % (
            pIso, self.heteroAtoms[pIso].mzOffset, self.heteroAtoms[pIso].relativeAbundance * 100.,
            self.heteroAtoms[pIso].minCount, self.heteroAtoms[pIso].maxCount) for pIso in self.heteroAtoms])), style=getSampleStyleSheet()["Normal"])
        w, h = p.wrap(300, 60)
        p.wrapOn(pdf, 300, 60)
        p.drawOn(pdf, 240, currentHeight - h + 10);
        currentHeight -= max(35, h + 10)

        pdf.drawString(50, currentHeight, "Non-targeted feature grouping");
        pdf.line(50, currentHeight - 2, 540, currentHeight - 2);
        currentHeight -= 20

        pdf.drawString(70, currentHeight, "Min. correlation")
        pdf.drawString(240, currentHeight, "%.2f" % self.minCorrelation);
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Min. correlation connections")
        pdf.drawString(240, currentHeight, "%.2f" % self.minCorrelationConnections);
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Adducts")
        currentHeight -=15
        p = Paragraph(str(", ".join(
            ["%s (m/z: %.4f, z:%d%s)" % (ad.name, ad.mzoffset, ad.charge, ad.polarity) for ad in self.adducts])), style=getSampleStyleSheet()["Normal"])
        w, h = p.wrap(460, 60)
        p.wrapOn(pdf, 460, 60)
        p.drawOn(pdf, 80, currentHeight - h + 10);
        currentHeight -= max(20, h + 10)

        pdf.drawString(70, currentHeight, "Elements")
        currentHeight -=15
        p = Paragraph(
            str(", ".join(["%s (m: %.4f)" % (el, self.elements[el].weight) for el in self.elements.keys()])), style=getSampleStyleSheet()["Normal"])
        w, h = p.wrap(460, 60)
        p.wrapOn(pdf, 460, 60)
        p.drawOn(pdf, 80, currentHeight - h + 10);
        currentHeight -= max(35, h + 10)

        pdf.drawString(70, currentHeight, "Simplify in-source fragments: "+str(self.simplifyInSourceFragments))
        currentHeight -= 15

        pdf.drawString(50, currentHeight, "Software");
        pdf.line(50, currentHeight - 2, 540, currentHeight - 2);
        currentHeight -= 20

        pdf.drawString(70, currentHeight, "MetExtract version")
        pdf.drawString(240, currentHeight, str(self.meVersion));
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "R-Version")
        pdf.drawString(240, currentHeight, str(self.rVersion));
        currentHeight -= 15


        p = Paragraph("UUID_ext: "+self.processingUUID, style=getSampleStyleSheet()["Normal"])
        w, h = p.wrap(450, 120)
        p.wrapOn(pdf, 450, 120)
        p.drawOn(pdf, 60, currentHeight - h + 10);
        currentHeight -= 15

        pdf.showPage()

    # creates a new PDF page which contains the tracer used in this experiment
    def writeCurrentTracerToPDF(self, pdf, tracer):
        currentHeight = 800
        pdf.drawString(50, currentHeight, "Tracer: %s" % tracer.name);
        pdf.line(50, currentHeight - 2, 540, currentHeight - 1);
        currentHeight -= 20

        pdf.drawString(70, currentHeight, "Labelling")
        pdf.drawString(240, currentHeight, "%s, %s" % (tracer.isotopeA, tracer.isotopeB));
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Delta m/z")
        pdf.drawString(240, currentHeight, "%.5f" % self.xOffset);
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Element count")
        pdf.drawString(240, currentHeight, "%d" % tracer.elementCount);
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Purity %s" % tracer.isotopeA)
        pdf.drawString(240, currentHeight, "%.2f%%" % (tracer.enrichmentA * 100.));
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Purity %s" % tracer.isotopeB)
        pdf.drawString(240, currentHeight, "%.2f%%" % (tracer.enrichmentB * 100.));
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Monoisotopic ratio (%s:%s)" % (tracer.isotopeA, tracer.isotopeB))
        pdf.drawString(240, currentHeight,
                       "%.5f (%.2f:%.2f [v/v])" % (tracer.monoisotopicRatio, tracer.amountA, tracer.amountB));
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Max. intensity error")
        pdf.drawString(240, currentHeight, "%.1f%%, %.1f%% (%.3f - %.3f)" % (
            tracer.maxRelNegBias * 100., tracer.maxRelPosBias * 100., tracer.monoisotopicRatio *tracer.maxRelNegBias,
            tracer.monoisotopicRatio *tracer.maxRelPosBias));
        currentHeight -= 15

        pdf.showPage()

    # data processing step 1: searches each mass spectrum for isotope patterns of native and highly isotope enriched
    # metabolite ions. The actual calculation and processing of the data is performed in the file runIdentification_matchPartners.py.
    # The positive and negative ionisation modes are processed separately.
    def findSignalPairs(self, mzxml, tracer, rules, reportFunction=None):
        mzs = []
        posFound = 0
        negFound = 0


        def reportFunctionHelper(curVal, text):
            reportFunction(curVal, text)

        if self.positiveScanEvent != "None":
            if self.negativeScanEvent != "None":
                def reportFunctionHelper(curVal, text):
                    reportFunction(curVal / 2, text)

            p = matchPartners(mzXMLData=mzxml, rules=rules,
                              labellingIsotopeB=self.labellingIsotopeB,
                              useCIsotopePatternValidation=self.useCIsotopePatternValidation,
                              intensityThres=self.intensityThreshold,
                              maxLoading=self.maxLoading,
                              xCounts=self.xCounts,
                              xOffset=self.xOffset,
                              ppm=self.ppm,
                              purityN=tracer.enrichmentA if self.metabolisationExperiment > 1 else self.purityN,
                              purityL=tracer.enrichmentB if self.metabolisationExperiment > 1 else self.purityL,
                              startTime=self.startTime, stopTime=self.stopTime,
                              filterLine=self.positiveScanEvent,
                              ionMode="+",
                              metabolisationExperiment=self.metabolisationExperiment,
                              reportFunction=reportFunctionHelper)
            posFound = len(p)
            mzs.extend(p)

        def reportFunctionHelper(curVal, text):
            reportFunction(curVal, text)

        if self.negativeScanEvent != "None":
            if self.positiveScanEvent != "None":
                def reportFunctionHelper(curVal, text):
                    reportFunction(.5 + curVal / 2, text)

            n = matchPartners(mzXMLData=mzxml, rules=rules,
                              labellingIsotopeB=self.labellingIsotopeB,
                              useCIsotopePatternValidation=self.useCIsotopePatternValidation,
                              intensityThres=self.intensityThreshold,
                              maxLoading=self.maxLoading,
                              xCounts=self.xCounts,
                              xOffset=self.xOffset,
                              ppm=self.ppm,
                              purityN=tracer.enrichmentA if self.metabolisationExperiment > 1 else self.purityN,
                              purityL=tracer.enrichmentB if self.metabolisationExperiment > 1 else self.purityL,
                              startTime=self.startTime, stopTime=self.stopTime,
                              filterLine=self.negativeScanEvent,
                              ionMode="-",
                              metabolisationExperiment=self.metabolisationExperiment,
                              reportFunction=reportFunctionHelper)
            negFound = len(n)
            mzs.extend(n)

        return mzs, negFound, posFound

    # store detected signal pairs (1st data processing step) in the database
    def writeSignalPairsToDB(self, mzs, mzxml):
        conn = connect(self.file + getDBSuffix())
        curs = conn.cursor()

        for mz in mzs:
            mz.id = self.curMZId

            scanEvent = ""
            if mz.ionMode == "+":
                scanEvent = self.positiveScanEvent
            elif mz.ionMode == "-":
                scanEvent = self.negativeScanEvent

            SQLInsert(curs, "MZs", id=mz.id, mz=mz.mz, similarityObject=mz.similarityString, scanid=mzxml.getIthMS1Scan(mz.scanIndex, scanEvent).id, scanTime=mzxml.getIthMS1Scan(mz.scanIndex, scanEvent).retention_time,
                      loading=mz.loading, intensity=mz.nIntensity, ionMode=mz.ionMode, type=mz.type, otherIsotopologs=base64.b64encode(dumps(mz.otherIsotopologs)))

            self.curMZId = self.curMZId + 1

        conn.commit()
        curs.close()
        conn.close()

    # data processing step 2: cluster detected signal pairs with HCA
    def clusterFeaturePairs(self, mzs, reportFunction=None):
        mzbins = {}
        mzbins['+'] = []
        mzbins['-'] = []

        similarityStrings=list(set([mz.similarityString for mz in mzs]))

        # cluster each detected number of carbon atoms separately
        donei=0
        for similarityString in similarityStrings:
            if reportFunction is not None:
                reportFunction(donei / len(similarityStrings), "Current: %s" % similarityString)
            # cluster each detected number of loadings separately
            for loading in range(self.maxLoading, 0, -1):
                for ionMode in ['+', '-']:

                    xAct = sorted(
                        [mz for mz in mzs if mz.ionMode == ionMode and mz.similarityString == similarityString and mz.loading == loading],
                        key=lambda x: x.mz)

                    doClusterings=[]

                    if len(xAct) > 0:
                        lastUnused = 0
                        usedCount = 0

                        # non-consecutive mz regions are prematurely separated without HCA to reduce the number of
                        # signal pairs for HCA (faster processing)
                        for i in range(1, len(xAct)):
                            ppmDiff = (xAct[i].mz - xAct[i - 1].mz) * 1000000 / xAct[i - 1].mz
                            if ppmDiff > self.clustPPM:
                                doClusterings.append(xAct[lastUnused:i])

                                usedCount = usedCount + i - lastUnused
                                lastUnused = i

                        # cluster remainign signal pairs
                        if lastUnused < len(xAct):
                            doClusterings.append(xAct[lastUnused:len(xAct)])
                            usedCount = usedCount + len(xAct) - lastUnused

                        assert usedCount == len(xAct)

                        temp=doClusterings
                        doClusterings=[]
                        for t in temp:
                            xAct=sorted([mz for mz in t], key=lambda x: x.scanIndex)

                            lastUnused = 0
                            usedCount = 0

                            for i in range(1, len(xAct)):
                                scanDiff = xAct[i].scanIndex - xAct[i - 1].scanIndex
                                if scanDiff > 20:
                                    doClusterings.append(xAct[lastUnused:i])

                                    usedCount = usedCount + i - lastUnused
                                    lastUnused = i

                            # cluster remainign signal pairs
                            if lastUnused < len(xAct):
                                doClusterings.append(xAct[lastUnused:len(xAct)])
                                usedCount = usedCount + len(xAct) - lastUnused

                            assert usedCount == len(xAct)


                        for doClustering in doClusterings:
                            hc = HierarchicalClustering(doClustering,
                                                        dist=lambda x, y: x.getValue() - y.getValue(),
                                                        val=lambda x: x.mz,
                                                        mean=lambda x, y: x / y,
                                                        add=lambda x, y: x + y)

                            for n in cutTreeSized(hc.getTree(), self.clustPPM):
                                mzbins[ionMode].append(n)

        return mzbins

    # store signal pair clusters (2nd data processing step) in the database
    def writeFeaturePairClustersToDB(self, mzbins):
        conn = connect(self.file + getDBSuffix())
        curs = conn.cursor()

        for ionMode in ['+', '-']:
            for mzbin in mzbins[ionMode]:
                SQLInsert(curs, "MZBins", id=self.curMZBinId, mz=mzbin.getValue(), ionMode=ionMode)
                for kid in mzbin.getKids():
                    SQLInsert(curs, "MZBinsKids", mzbinID=self.curMZBinId, mzID=kid.getObject().id)
                self.curMZBinId = self.curMZBinId + 1

        conn.commit()
        curs.close()
        conn.close()

    def removeImpossibleFeaturePairClusters(self, mzbins):
        mzBinsNew={}
        for ionMode in mzbins.keys():
            mzBinsNew[ionMode]=[mzbin for mzbin in mzbins[ionMode] if len(mzbin.getKids())>=self.minSpectraCount]
        return mzBinsNew





    # EXPERIMENTAL: calulcate the mean intensity ratio of two chromatographic peaks at the same retention time
    # in two different EICs. Thus, the internal standardisation is not calculated using the peak area
    # but rather the ratio of each MS peak that contributes to the chromatographic peak.
    def getMeanRatioOfScans(self, eicA, eicB, lib, rib, perfWeighted=True, minInt=1000, minRatiosNecessary=3):
        try:
            if perfWeighted:
                sumEIC = sum(eicA[lib:rib])
                sumEICL = sum(eicB[lib:rib])

                normEIC = eicA
                if sumEICL > sumEIC:
                    normEIC = eicB

                os=[o for o in range(lib, rib) if eicA[o] >= minInt and eicB[o] >= minInt and eicB[o]>0]
                normSum = sum([normEIC[o] for o in os])
                weights = [1. * normEIC[o] / normSum for o in os]
                ratios = [1. * eicA[o] / eicB[o] for o in os if eicB[o]>minInt and eicB[o]>0]


                if len(ratios) >= minRatiosNecessary:
                    assert len(ratios) == len(weights)
                    ratio = sum([ratios[o] * weights[o] for o in range(len(ratios))])
                    deviation = sum([ratios[o] * weights[o] for o in range(len(ratios))])
                    return ratio
                else:
                    return -1
            else:
                ratios = [(eicA[o] / eicB[o]) for o in range(lib, rib) if eicA[o] >= minInt and eicB[o] >= minInt]
                if len(ratios) >= minRatiosNecessary:
                    return mean(ratios)
                else:
                    return -1

        except IndexError:
            return -1

    def __getEICFor(self, mz, mzxml, scanEvent):
        eic, times, scanIds, mzs = mzxml.getEIC(mz, self.chromPeakPPM, filterLine=scanEvent)
        eicBaseline = self.BL.getBaseline(copy(eic), times)
        eicSmoothed = smoothDataSeries(times, copy(eic), windowLen=self.eicSmoothingWindowSize,
                                       window=self.eicSmoothingWindow, polynom=self.eicSmoothingPolynom)
        return eic, eicSmoothed, times

    def __getChromPeaksFor(self, mz, mzxml, scanEvent):
        # extract the EIC of the native ion and detect its chromatographic peaks
        eic, eicSmoothed, times = self.__getEICFor(mz, mzxml, scanEvent)
        peaks = []
        try:
            peaks = self.CP.getPeaksFor(times, eicSmoothed, scales=self.scales, snrTh=self.snrTh, startIndex=0, endIndex=len(eic) - 1)
        except Exception as ex:
            self.printMessage("Error getting chromatographic peaks for mz %.5f, ScanEvent %s: %s" % (mz, scanEvent, str(ex)), type="error")

        return eicSmoothed, peaks, times

    def __findBestArtificialShift(self, eicN, eicL, lb, rb, shiftFrom=0, shiftTo=0):
        correlations = []

        for artShift in range(shiftFrom, shiftTo + 1):
            peakN = eicN[lb:rb]
            peakL = eicL[(lb + artShift):(rb + artShift)]
            silRatios = [peakN[i] / peakL[i] for i in range(int(len(peakN) * .25), int(len(peakN) * .75) + 1) if peakL[i] > 0 and peakN[i] > 0]
            correlations.append(Bunch(correlation=corr(peakN, peakL), artificialShift=artShift, silRatios=silRatios,
                                      peakNInts=[peakN[i] for i in range(int(len(peakN) * .25), int(len(peakN) * .75) + 1) if peakL[i] > 0 and peakN[i] > 0],
                                      peakLInts=[peakL[i] for i in range(int(len(peakN) * .25), int(len(peakN) * .75) + 1) if peakL[i] > 0 and peakN[i] > 0]))
        bestFit = max(correlations, key=lambda x: x.correlation)

        return bestFit

    # data processing step 3: for each signal pair cluster extract the EICs, detect chromatographic peaks
    # present in both EICs at approximately the same retention time and very their chromatographic peak shapes.
    # if all criteria are passed, write this detected feature pair to the database
    def findChromatographicPeaksAndWriteToDB(self, mzbins, mzxml, rules, reportFunction=None):
        conn = connect(self.file + getDBSuffix())
        curs = conn.cursor()
        chromPeaks = []
        totalBins = (len(mzbins['+']) + len(mzbins['-']))

        # process each signal pair cluster
        for ionMode in ['+', '-']:
            scanEvent = self.positiveScanEvent if ionMode == "+" else self.negativeScanEvent

            if scanEvent != "None":
                for mzbin, i in zip(mzbins[ionMode], range(0, len(mzbins[ionMode]))):
                    if reportFunction is not None:
                        doneBins = i
                        if ionMode == '-':
                            doneBins += len(mzbins['+'])

                        reportFunction(1. * doneBins / totalBins, "%d/%d mzbins done / %d features found so far" % (doneBins, totalBins, len(chromPeaks)))

                    kids = mzbin.getKids()
                    if len(kids) < self.minSpectraCount:
                        continue

                    otherIsotopologs = kids[0].getObject().otherIsotopologs
                    similarityString = kids[0].getObject().similarityString
                    loading=kids[0].getObject().loading
                    assert all([kid.getObject().similarityString == similarityString for kid in kids])


                    # calulcate mean mz value for this signal pair cluster
                    mz = weightedMean([kid.getObject().mz for kid in kids],[kid.getObject().nIntensity for kid in kids])
                    eic, peaks, times = self.__getChromPeaksFor(mz, mzxml, scanEvent)
                    for peak in peaks:

                        peakAreas={}
                        peakAreas["X"]=Bunch(intensity=peak.peakArea)

                        ## Test if pattern was present in at least n scans
                        assignedScans=[]
                        for kid in kids:
                            kido = kid.getObject()
                            if abs(kido.scanIndex - peak.peakIndex) < peak.peakScale * 2:
                                assignedScans.append(kido.scanIndex)
                        if len(assignedScans)<self.minSpectraCount:
                            continue

                        ## Test other isotopologs that are required
                        peak.foundMatches = {}
                        allIsotopologsFound=True
                        for otherIsotopolog in otherIsotopologs:

                            mzl=mz + otherIsotopologs[otherIsotopolog].mzInc / loading
                            eicL, peaksL, times=self.__getChromPeaksFor(mzl, mzxml, scanEvent)

                            # match detected chromatographic peaks from both EICs
                            closestMatch=-1
                            closestOffset=10000000000

                            # search closest peak pair
                            for li, peakL in enumerate(peaksL):
                                if abs(peak.peakIndex - peakL.peakIndex) <= self.peakCenterError:
                                    if closestMatch==-1 or closestOffset > abs(peak.peakIndex - peakL.peakIndex):
                                        closestMatch=peakL
                                        closestOffset=abs(peak.peakIndex - peakL.peakIndex)

                            if closestMatch==-1:
                                if otherIsotopologs[otherIsotopolog].requiredChromPeak:
                                    allIsotopologsFound=False
                                    break
                                continue


                            # Test correlation and find artificial shift
                            lb = int(max(0, min(peak.peakIndex - peak.peakLeftFlank, peakL.peakIndex - peakL.peakLeftFlank)))
                            rb = int(min(max(peak.peakIndex + peak.peakRightFlank, peakL.peakIndex + peakL.peakRightFlank) + 1, len(eic) - 1))
                            co = self.__findBestArtificialShift(eic, eicL, lb, rb, self.artificialMPshift_start, self.artificialMPshift_stop)
                            # check peak shape (Pearson correlation)
                            if co.correlation <= self.minPeakCorr:
                                if otherIsotopologs[otherIsotopolog].requiredChromPeak:
                                    allIsotopologsFound = False
                                    break
                                continue

                            peak.foundMatches[otherIsotopolog] = closestMatch
                            closestMatch.peaksCorr=co.correlation
                            closestMatch.artificialShift=co.artificialShift
                            closestMatch.peaksRatio=self.getMeanRatioOfScans(eicL, eic, lb, rb)
                            peakAreas[otherIsotopolog]=Bunch(intensity=closestMatch.peakArea)



                        if allIsotopologsFound and (not(self.checkPeaksRatio) or RuleMatcher(rules, ppm=self.ppm).checkChromPeaks(peakAreas)):
                            feat=ChromPeakFeature(id=len(chromPeaks)+1, fGroupID=len(chromPeaks)+1, assignedName="",
                                 mz=mz, loading=loading, ionMode=ionMode, similarityString=similarityString,
                                 PeakCenter=peak.peakIndex, PeakCenterMin=times[peak.peakIndex]/60., PeakScale=peak.peakScale,
                                 SNR=peak.peakSNR, PeakArea=peak.peakArea, PeakAbundance=eic[peak.peakIndex],
                                 heteroIsotoplogs={}, assignedMZs=[], comments=[],
                                 foundMatches=peak.foundMatches)
                            chromPeaks.append(feat)

                            curs.execute("INSERT INTO chromPeaks (id, fGroupID, assignedName, mz, loading, ionMode, similarityString, PeakCenter, PeakCenterMin, PeakScale, SNR, PeakArea, PeakAbundance, "
                                         "heteroIsotoplogues, assignedMZs, comments, foundMatches) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                                         (feat.id, feat.fGroupID, feat.assignedName, feat.mz, feat.loading, feat.ionMode, feat.similarityString, feat.PeakCenter,
                                          feat.PeakCenterMin, feat.PeakScale, feat.SNR, feat.PeakArea, feat.PeakAbundance,
                                          base64.b64encode(dumps(feat.heteroIsotoplogs)), base64.b64encode(dumps(feat.assignedMZs)), base64.b64encode(dumps(feat.comments)),
                                          base64.b64encode(dumps(feat.foundMatches)))  )

                            if True:
                                self.printMessage("\n\n%d: Peak similarities found: mz %.5f, RT %.2f min, PeakArea %s, ScanEvent %s"%(len(chromPeaks), mz, times[peak.peakIndex]/60., peak.peakArea, scanEvent), type="warning")
                                for iso in natSort(peak.foundMatches.keys()):
                                    peakB=peak.foundMatches[iso]
                                    self.printMessage("    --> %s, RT %5.2f min, Area %12.1f (AreaRatio %8.3f%%, IntensityRatio %8.3f%%), Pearson correlation %6.3f, Artificial shift %4d" %
                                                      (iso, times[peakB.peakIndex] / 60., peakB.peakArea,100*peakB.peakArea/peak.peakArea, 100*peakB.peaksRatio, peakB.peaksCorr, peakB.artificialShift), type="warning")
                                self.printMessage("Rules are "+str(RuleMatcher(rules, ppm=self.ppm).checkChromPeaks(peakAreas))+" using the peak areas")

        conn.commit()
        curs.close()
        conn.close()
        return chromPeaks


    # data processing step 6: convolute different feature pairs into feature groups using the chromatographic
    # profiles of different metabolite ions
    def groupFeaturePairsUntargetedAndWriteToDB(self, chromPeaks, mzxml, reportFunction=None):
        try:
            conn = connect(self.file + getDBSuffix())
            curs = conn.cursor()

            nodes = {}
            correlations = {}

            allPeaks = {}
            for peak in chromPeaks:
                nodes[peak.id] = []
                allPeaks[peak.id] = peak
                peak.correlationsToOthers=[]

            # compare all detected feature pairs at approximately the same retention time
            for piA in range(len(chromPeaks)):
                peakA = chromPeaks[piA]
                eicA, eicSmoothedA, timesA = self.__getEICFor(peakA.mz, mzxml, self.positiveScanEvent if peakA.ionMode == "+" else self.negativeScanEvent)

                if reportFunction is not None:
                    reportFunction(.7 * piA / len(chromPeaks), "Matching features (%d remaining)" % (len(chromPeaks) - piA))

                if peakA.id not in correlations.keys():
                    correlations[peakA.id]={}

                for peakB in chromPeaks:
                    eicB, eicSmoothedB, timesB = self.__getEICFor(peakB.mz, mzxml,self.positiveScanEvent if peakB.ionMode == "+" else self.negativeScanEvent)

                    if peakB.id not in correlations.keys():
                        correlations[peakB.id]={}

                    if peakA.mz < peakB.mz:
                        if abs(peakA.PeakCenter - peakB.PeakCenter) < self.peakCenterError:
                            bmin = int(max(0, mean([peakA.PeakCenter - 1 * peakA.BorderLeft,
                                                    peakB.PeakCenter - 1 * peakB.BorderLeft])))
                            bmax = int(min(len(peakB.NXICSmoothed) - 1, mean([peakB.PeakCenter + 1 * peak.NBorderRight,
                                                                      peakA.PeakCenter + 1 * peakA.BorderRight])))

                            pb = corr(eicSmoothedA[bmin:bmax], eicSmoothedB[bmin:bmax])
                            if str(pb)=="nan":
                                pb=-1

                            correlations[peakA.id][peakB.id]=pb
                            correlations[peakB.id][peakA.id]=pb

                            ## TODO
                            silRatiosA=1#peakA.silRatios.silRatios
                            silRatiosB=1#peakB.silRatios.silRatios

                            meanSilRatioA = 1#weightedMean(silRatiosA, peakA.silRatios.peakNInts)
                            meanSilRatioB = 1#weightedMean(silRatiosB, peakB.silRatios.peakNInts)

                            silRatiosFold=0
                            silRatiosSD=1

                            try:

                                silRatiosFold=1#max(meanSilRatioA, meanSilRatioB)/min(meanSilRatioA, meanSilRatioB)
                                silRatiosSD=0.1#weightedSd([abs(r-meanSilRatioA)/meanSilRatioA for r in silRatiosA if min(r, meanSilRatioA)>0]+
                                               #        [abs(r-meanSilRatioB)/meanSilRatioB for r in silRatiosB if min(r, meanSilRatioB)>0],
                                               #        peakA.silRatios.peakNInts+peakB.silRatios.peakNInts)


                                # check for similar chromatographic peak profile and similar native to labeled ratio
                                if pb >= self.minCorrelation and silRatiosFold <= 1+max(0.25, 3*silRatiosSD):
                                    nodes[peakA.id].append(peakB.id)
                                    nodes[peakB.id].append(peakA.id)

                                    peakB.correlationsToOthers.append(peakA.id)
                                    peakA.correlationsToOthers.append(peakB.id)

                            except Exception as e:
                                logging.error("Error while convoluting two feature pairs, skipping.. (%s)"%str(e))

                            try:
                                SQLInsert(curs, "featurefeatures", fID1=peakA.id, fID2=peakB.id, corr=pb, silRatioValue=silRatiosFold)
                            except Exception as e:
                                logging.error("Error while convoluting two feature pairs, skipping.. (%s)"%str(e))
                                SQLInsert(curs, "featurefeatures", fID1=peakA.id, fID2=peakB.id, corr=0, silRatioValue=0)

            self.postMessageToProgressWrapper("text", "%s: Convoluting feature groups" % tracer.name)


            for k in nodes.keys():
                uniq = []
                for u in nodes[k]:
                    if u not in uniq:
                        uniq.append(u)
                nodes[k] = uniq

            # get subgraphs from the feature pair graph. Each subgraph represents one convoluted
            # feature group
            tGroups = getSubGraphs(nodes)

            def splitGroupWithHCA(tGroup, correlations):
                # if 1 or two feature pairs are in a group, automatically use them (no further splitting possible)
                if len(tGroup)<=2:
                    return [tGroup]

                # construct 2-dimensional data matrix
                data=[]
                for tG1 in tGroup:
                    t=[]
                    for tG2 in tGroup:
                        if tG1 in correlations.keys() and tG2 in correlations[tG1].keys():
                            t.append(correlations[tG1][tG2])
                        elif tG1==tG2:
                            t.append(1.)
                        else:
                            t.append(0.)
                    data.append(t)

                # calculate HCA tree using the correlations between all feature pairs
                hc=HCA_general.HCA_generic()
                tree=hc.generateTree(objs=data, ids=tGroup)
                #hc.plotTree(tree)
                #print "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"

                # split the HCA tree in sub-clusters
                def checkSubCluster(tree, hca, corrThreshold, cutOffMinRatio):
                    if isinstance(tree, HCA_general.HCALeaf):
                        return False
                    elif isinstance(tree, HCA_general.HCAComposite):
                        corrsLKid=hca.link(tree.getLeftKid())
                        corrs=hca.link(tree)
                        corrsRKid=hca.link(tree.getRightKid())

                        aboveThresholdLKid=sum([corr>corrThreshold for corr in corrsLKid])
                        aboveThreshold=sum([corr>corrThreshold for corr in corrs])
                        aboveThresholdRKid=sum([corr>corrThreshold for corr in corrsRKid])

                        #print aboveThresholdLKid, aboveThreshold, aboveThresholdRKid

                    if (aboveThresholdLKid*1./len(corrs))>=cutOffMinRatio and (aboveThreshold*1./len(corrs))>=cutOffMinRatio and (aboveThresholdRKid*1./len(corrs))>=cutOffMinRatio:
                        return False
                    else:
                        return True
                #subClusts=hc.splitTreeWithCallbackBottomUp(tree,
                #                                   CallBackMethod(_target=checkSubCluster, corrThreshold=self.minCorrelation, cutOffMinRatio=self.minCorrelationConnections).getRunMethod())
                subClusts=hc.splitTreeWithCallback(tree,
                                                   CallBackMethod(_target=checkSubCluster, corrThreshold=self.minCorrelation, cutOffMinRatio=self.minCorrelationConnections).getRunMethod(),
                                                   recursive=False)

                # convert the subclusters into arrays of feature pairs belonging to the same metabolite
                return [[leaf.getID() for leaf in subClust.getLeaves()] for subClust in subClusts]


            groups=[]
            done=0
            for tGroup in tGroups:
                #print "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
                cGroups=[tGroup]
                while len(cGroups)>0:
                    #print "HCA with", len(cGroups)
                    gGroup=cGroups.pop(0)
                    sGroups=splitGroupWithHCA(gGroup, correlations)
                    #print "first group", len(gGroup), gGroup, "split into", sGroups

                    if len(sGroups)==1:
                        groups.append(sGroups[0])
                    else:
                        cGroups.extend(sGroups)

                done=done+1
                self.postMessageToProgressWrapper("text", "%s: Convoluting feature groups (%d/%d done)" % (tracer.name, done, len(tGroups)))

            if True:
                done=0
                for gi in range(len(groups)):
                    group = groups[gi]
                    if reportFunction is not None:
                        reportFunction(0.7 + .3 * piA / len(chromPeaks),"Searching for relationships (%d remaining)" % (len(groups) - gi))

                    # first, search for adduct relationships (same number of carbon atoms) in each convoluted
                    # feature group

                    peaksInGroup={}
                    for a in group:
                        peaksInGroup[a]=allPeaks[a]

                        temp=[]
                        for j in range(len(peaksInGroup[a].correlationsToOthers)):
                            if peaksInGroup[a].correlationsToOthers[j] in group:
                                temp.append(peaksInGroup[a].correlationsToOthers[j])

                        peaksInGroup[a].correlationsToOthers=temp

                    done=done+1
                    self.postMessageToProgressWrapper("text", "%s: Annotating feature groups (%d/%d done)" % (tracer.name, done, len(groups)))


            for peak in chromPeaks:
                adds = countEntries(peak.adducts)
                peak.adducts = adds.keys()
                curs.execute("UPDATE chromPeaks SET adducts='%s', fDesc='%s', correlationsToOthers='%s' WHERE id=%d" % (
                              base64.b64encode(dumps(peak.adducts)), base64.b64encode(dumps(peak.fDesc)), base64.b64encode(dumps(peak.correlationsToOthers)), peak.id))
                curs.execute("UPDATE allChromPeaks SET adducts='%s', fDesc='%s' WHERE id=%d" % (
                              base64.b64encode(dumps(peak.adducts)), base64.b64encode(dumps(peak.fDesc)), peak.id))

                curs.execute("UPDATE chromPeaks SET heteroAtomsFeaturePairs=? WHERE id=?",
                             (base64.b64encode(dumps(peak.heteroAtomsFeaturePairs)), peak.id))
                curs.execute("UPDATE allChromPeaks SET heteroAtomsFeaturePairs=? WHERE id=?",
                             (base64.b64encode(dumps(peak.heteroAtomsFeaturePairs)), peak.id))

            # store feature group in the database
            for groupi, group in enumerate(sorted(groups, key=lambda x: sum([allPeaks[p].NPeakCenterMin / 60. for p in x]) / len(x))):
                SQLInsert(curs, "featureGroups", id=self.curFeatureGroupID, featureName="fg_%d"%self.curFeatureGroupID, tracer=tracerID)
                groupMeanElutionIndex = 0

                hasPos = False
                hasNeg = False

                for p in sorted(group, key=lambda x: allPeaks[x].mz):
                    hasPos = hasPos or allPeaks[p].ionMode == "+"
                    hasNeg = hasNeg or allPeaks[p].ionMode == "-"

                    groupMeanElutionIndex += allPeaks[p].NPeakCenter

                    allPeaks[p].fGroupID = self.curFeatureGroupID
                    SQLInsert(curs, "featureGroupFeatures", fID=p, fDesc="", fGroupID=self.curFeatureGroupID)

                groupMeanElutionIndex = groupMeanElutionIndex / len(group)

                # store one positve and one negative ionisation mode MS scan in the database for
                # later visualisation (one for each convoluted feature group)
                if hasPos:
                    scan = mzxml.getIthMS1Scan(int(groupMeanElutionIndex), self.positiveScanEvent)

                    SQLInsert(curs, "massspectrum", mID=self.curMassSpectrumID, fgID=self.curFeatureGroupID, time=scan.retention_time,
                              mzs=";".join([str(u) for u in scan.mz_list]), intensities=";".join([str(u) for u in scan.intensity_list]), ionMode="+")
                    self.curMassSpectrumID += 1
                if hasNeg:
                    scan = mzxml.getIthMS1Scan(int(groupMeanElutionIndex), self.negativeScanEvent)

                    SQLInsert(curs, "massspectrum", mID=self.curMassSpectrumID, fgID=self.curFeatureGroupID, time=scan.retention_time,
                              mzs=";".join([str(u) for u in scan.mz_list]), intensities=";".join([str(u) for u in scan.intensity_list]), ionMode="-")
                    self.curMassSpectrumID += 1

                self.curFeatureGroupID += 1

            conn.commit()
            self.printMessage("%s: Feature grouping done. " % tracer.name + str(len(groups)) + " feature groups", type="info")

            conn.commit()
            curs.close()
            conn.close()

        except Exception as ex:
            import traceback
            traceback.print_exc()

            self.printMessage("Error in %s: %s" % (self.file, str(ex)), type="error")
            self.postMessageToProgressWrapper("failed", self.pID)

    # store one MS scan for each detected feature pair in the database
    def writeMassSpectraToDB(self, chromPeaks, mzxml, reportFunction=None):
        conn = connect(self.file + getDBSuffix())
        curs = conn.cursor()
        massSpectraWrittenPos = {}
        massSpectraWrittenNeg = {}
        for pi in range(len(chromPeaks)):
            peak = chromPeaks[pi]
            if reportFunction is not None:
                reportFunction(1. * pi / len(chromPeaks), "%d feature pairs remaining" % (len(chromPeaks) - pi))

            scanEvent = ""
            iMode = ''
            massSpectraWritten = None
            if peak.ionMode == "+":
                iMode = '+'
                scanEvent = self.positiveScanEvent
                massSpectraWritten = massSpectraWrittenPos
            elif peak.ionMode == "-":
                iMode = '-'
                scanEvent = self.negativeScanEvent
                massSpectraWritten = massSpectraWrittenNeg

            if peak.NPeakCenter not in massSpectraWritten.keys():
                scan = mzxml.getIthMS1Scan(peak.NPeakCenter, scanEvent)

                SQLInsert(curs, "massspectrum", mID=self.curMassSpectrumID, fgID=-1, time=scan.retention_time, mzs=";".join([str(u) for u in scan.mz_list]),
                          intensities=";".join([str(u) for u in scan.intensity_list]), ionMode=iMode)

                massSpectraWritten[peak.NPeakCenter] = self.curMassSpectrumID
                self.curMassSpectrumID = self.curMassSpectrumID + 1

            curs.execute(
                "UPDATE chromPeaks SET massSpectrumID=%d WHERE id=%d" % (massSpectraWritten[peak.NPeakCenter], peak.id))

        conn.commit()
        curs.close()
        conn.close()

    ## write a new featureML file
    def writeResultsToFeatureML(self, forFile):
        conn = connect(forFile + getDBSuffix())
        curs = conn.cursor()

        features=[]
        cols = OrderedDict([("id", "id"),
                            ("ogroup", "fGroupID"),
                            ("mz", "mz"),
                            ("rt", "PeakCenterMin"),
                            ("lmz", "mz"),
                            ("charge", "loading"),
                            ("name", "id")
                            ])
        for chromPeak in SQLSelectAsObject(curs, "SELECT %s FROM chromPeaks" % (", ".join(["%s AS %s" % (cols[col], col) for col in cols.keys()])), newObject=Bunch):
            chromPeak.ogroup=-1
            chromPeak.Xn=0
            chromPeak.rt=chromPeak.rt*60.
            features.append(chromPeak)

        exportAsFeatureML.writeFeatureListToFeatureML(features, forFile+".featureML", ppmPM=self.ppm, rtPM=0.25*60)

    # write feature pairs detected in this LC-HRMS data file into a new TSV file.
    # Each row represent one feature pair
    def writeResultsToTSVFile(self, forFile):
        conn = connect(forFile + getDBSuffix())
        curs = conn.cursor()

        chromPeaks = []
        configTracers = {}

        cols=OrderedDict([("Num",     "id"),
                          #("OGroup",  "fGroupID"),
                          ("MZ",      "mz"),
                          ("RT",      "PeakCenterMin"),
                          ("Label",   "similarityString"),
                          ("Z",       "loading"),
                          ("Pol",     "ionMode"),
                          ("Scale",   "PeakScale"),
                          ("Area",    "PeakArea")
        ])

        with open(forFile + ".tsv", "w") as tsvFile:
            tsvFile.write("\t".join(cols.keys()))
            tsvFile.write("\n")


            for chromPeak in SQLSelectAsObject(curs, "SELECT %s FROM chromPeaks"%(", ".join(["%s AS %s"%(cols[col], col) for col in cols.keys()])), newObject=Bunch):
                tsvFile.write("\t".join(str(chromPeak.__dict__.get(k)) for k in cols.keys()))
                tsvFile.write("\n")



            tsvFile.write("## MetExtract II %s\n"%(Bunch(MetExtractVersion=self.meVersion, RVersion=self.rVersion, UUID_ext=self.processingUUID).dumpAsJSon().replace("\"", "'")))

            processingParams=Bunch()
            processingParams.experimentOperator=self.experimentOperator
            processingParams.experimentID=self.experimentID
            processingParams.experimentComments=self.experimentComments
            processingParams.experimentName=self.experimentName
            tsvFile.write("## Experiment parameters %s\n"%(processingParams.dumpAsJSon().replace("\"", "'")))

            processingParams=Bunch()
            processingParams.startTime=self.startTime
            processingParams.stopTime=self.stopTime
            processingParams.positiveScanEvent=self.positiveScanEvent
            processingParams.negativeScanEvent=self.negativeScanEvent
            processingParams.metabolisationExperiment=self.metabolisationExperiment
            processingParams.intensityThreshold=self.intensityThreshold
            processingParams.intensityCutoff=self.intensityCutoff
            processingParams.maxLoading=self.maxLoading
            processingParams.xCounts=self.xCountsString
            processingParams.xoffset=self.xOffset
            processingParams.ppm=self.ppm
            processingParams.purityN=self.purityN
            processingParams.purityL=self.purityL

            #2. Results clustering
            processingParams.minSpectraCount=self.minSpectraCount
            processingParams.clustPPM=self.clustPPM

            #3. Peak detection
            processingParams.chromPeakPPM=self.chromPeakPPM

            processingParams.eicSmoothingWindow=self.eicSmoothingWindow
            processingParams.eicSmoothingWindowSize=self.eicSmoothingWindowSize
            processingParams.eicSmoothingPolynom=self.eicSmoothingPolynom
            processingParams.scales=self.scales
            processingParams.minCorr=self.minPeakCorr
            processingParams.minCorrelationConvolution=self.minCorrelation
            processingParams.minCorrelationConnections=self.minCorrelationConnections
            tsvFile.write("## Data processing parameters %s\n"%(processingParams.dumpAsJSon().replace("\"", "'")))



    # method, which is called by the multiprocessing module to actually process the LC-HRMS data
    def identify(self):
        try:

            start = time.time()

            # region Initialize data processing pipeline
            ######################################################################################

            self.printMessage("File: %s" % self.file, type="info")

            exec(self.rules) in globals()
            self.rules = rules

            self.postMessageToProgressWrapper("start")
            self.postMessageToProgressWrapper("max", 100.)
            self.postMessageToProgressWrapper("value", 0.)
            self.postMessageToProgressWrapper("text", "Initialising")

            if not USEGRADIENTDESCENDPEAKPICKING:
                self.CP = MassSpecWavelet(self.chromPeakFile)
            else:
                from chromPeakPicking.GradientPeaks import GradientPeaks
                self.CP=GradientPeaks()                                                                       ## generic gradient descend peak picking - do not use. Parameters need to be optimized
                self.CP=GradientPeaks(minInt=1000, minIntFlanks=1, minIncreaseRatio=.01)                                                                       ## LTQ Orbitrap XL
                self.CP=GradientPeaks(minInt=10000, minIntFlanks=10, minIncreaseRatio=.15, expTime=[10, 250]) ## Swiss Orbitrap HF data
                self.CP=GradientPeaks(minInt=1000, minIntFlanks=10, minIncreaseRatio=.05, minDelta=10000, expTime=[5, 150]) ## Bernhard HSS
                self.CP=GradientPeaks(minInt=1000,minIntFlanks=100,minIncreaseRatio=.5,minDelta=100,expTime=[5,150])    ##Lin
                #self.CP=GradientPeaks(minInt=5, minIntFlanks=2, minIncreaseRatio=.05, expTime=[15, 150], minDelta=1, minInflectionDelta=2) ## Roitinger
                #self.CP=GradientPeaks(minInt=10000, minIntFlanks=10, minIncreaseRatio=.05, expTime=[5, 45])       ## RNA

            self.BL = Baseline.Baseline()


            self.curPeakId = 1
            self.curEICId = 1
            self.curMZId = 1
            self.curMZBinId = 1
            self.curFeatureGroupID = 1
            self.curMassSpectrumID = 1
            # endregion


            # region Create results database
            ######################################################################################

            self.postMessageToProgressWrapper("text", "Creating results DB")

            if os.path.exists(self.file + getDBSuffix()) and os.path.isfile(self.file + getDBSuffix()):
                os.remove(self.file + getDBSuffix())
            conn = connect(self.file + getDBSuffix())
            curs = conn.cursor()

            self.writeConfigurationToDB(conn, curs)

            curs.close()
            conn.close()
            # endregion

            # region Parse mzXML file
            ######################################################################################

            self.postMessageToProgressWrapper("text", "Parsing chromatogram file")

            mzxml = self.parseMzXMLFile()


            # endregion

            # Start calculation
            ######################################################################################

            self.postMessageToProgressWrapper("text", "Starting data processing")

            if self.metabolisationExperiment:
                self.printMessage("Tracer: %s" % self.configuredTracer.name, type="info")
                self.xOffset = getIsotopeMass(self.configuredTracer.isotopeB)[0] - getIsotopeMass(self.configuredTracer.isotopeA)[0]

            else:
                # Full metabolome labelling experiment
                pass
            # endregion






            # region 1. Find 12C 13C partners in the mz dimension (0-25%)
            ######################################################################################

            self.postMessageToProgressWrapper("value", 0)
            self.postMessageToProgressWrapper("text", "Extracting signal pairs")

            def reportFunction(curVal, text):
                self.postMessageToProgressWrapper("value",0.25 * curVal)
                self.postMessageToProgressWrapper("text", "Extracting signal pairs (%s)" % (text))

            mzs, negFound, posFound = self.findSignalPairs(mzxml, self.configuredTracer, self.rules, reportFunction)
            self.writeSignalPairsToDB(mzs, mzxml)

            self.printMessage("Extracting signal pairs done. pos: %d neg: %d mzs (including mismatches)" % (posFound, negFound), type="info")
            # endregion

            # region 2. Cluster found mz values according to mz value and number of x atoms (25-35%)
            ######################################################################################

            self.postMessageToProgressWrapper("value", 0.25)
            self.postMessageToProgressWrapper("text", "Clustering found signal pairs")

            def reportFunction(curVal, text):
                self.postMessageToProgressWrapper("value", 0.25 + 0.1 * curVal)
                self.postMessageToProgressWrapper("text", "Clustering found signal pairs (%s)" % (text))

            mzbins = self.clusterFeaturePairs(mzs, reportFunction)
            self.writeFeaturePairClustersToDB(mzbins)
            mzbins = self.removeImpossibleFeaturePairClusters(mzbins)

            self.printMessage("Clustering found signal pairs done. pos: %d neg: %d mz bins (including mismatches)" % (len(mzbins['+']), len(mzbins['-'])), type="info")
            # endregion

            # region 3. Extract chromatographic peaks (35-65%)
            ######################################################################################

            self.postMessageToProgressWrapper("value", 0.35)
            self.postMessageToProgressWrapper("text", "Separating feature pairs")

            def reportFunction(curVal, text):
                self.postMessageToProgressWrapper("value",0.35 + 0.3 * curVal)
                self.postMessageToProgressWrapper("text", "Separating feature pairs (%s)" % (text))

            chromPeaks = self.findChromatographicPeaksAndWriteToDB(mzbins, mzxml, self.rules, reportFunction)

            self.printMessage("Separating feature pairs done. pos: %d neg: %d chromatographic peaks (including mismatches)" % (len([c for c in chromPeaks if c.ionMode == "+"]),len([c for c in chromPeaks if c.ionMode == "-"])), type="info")
            # endregion

            # region 4. Group feature pairs untargeted using chromatographic peak shape into feature groups (75-80%)
            ######################################################################################

            self.postMessageToProgressWrapper("value", 0.75)
            self.postMessageToProgressWrapper("text", "Grouping feature pairs")

            def reportFunction(curVal, text):
                self.postMessageToProgressWrapper("value", 0.75  + 0.05 * curVal)
                self.postMessageToProgressWrapper("text", "Grouping feature pairs (%s)" % (text))

            #self.groupFeaturePairsUntargetedAndWriteToDB(chromPeaks, mzxml, reportFunction)
            # endregion

            # Log time used for processing of individual files
            elapsed = (time.time() - start) / 60.
            hours = ""
            if elapsed >= 60.:
                if elapsed < 120.:
                    hours = "1 hour "
                else:
                    hours = "%d hours " % (elapsed // 60)
            mins = "%.2f min(s)" % (elapsed % 60.)

            self.printMessage("Calculations finished (%s%s).." % (hours, mins), type="info")
            # endregion

            # region 8. Write results to files (95-100%, without progress indicator)
            ######################################################################################

            # W.1 Write results to TSV File
            ##########################################################################################
            if self.writeTSV:
                self.postMessageToProgressWrapper("text", "Writing results to TSV..")
                self.writeResultsToTSVFile(self.file)

            # W.2 Write results to featureML File
            ##########################################################################################
            if self.writeFeatureML:
                self.postMessageToProgressWrapper("text", "Writing results to featureML..")

                self.writeResultsToFeatureML(self.file)

            # W.3 Write resutls to PDF
            ##########################################################################################
            if False:
                self.postMessageToProgressWrapper("text", "Writing results to PDF")

                def reportFunction(curVal, text):
                    self.postMessageToProgressWrapper("value", 90. + 5. * curVal)
                    self.postMessageToProgressWrapper("text", "Writing results to PDF (%s)" % text)

                self.writeResultsToPDF(mzxml, reportFunction=reportFunction)
            # endregion

            mzxml.freeMe();


            self.printMessage("%s done.." % self.file, type="info")
            self.postMessageToProgressWrapper("end")


        except Exception as ex:
            import traceback
            traceback.print_exc()

            self.printMessage("Error in %s: %s" % (self.file, str(ex)), type="error")
            self.postMessageToProgressWrapper("failed")


























    # annotate a metabolite group (consisting of ChromPeakPair instances) with the defined
    # hetero atoms based on detected feature pairs
    def annotateFeaturePairsWithHeteroAtoms(self, group, peaksInGroup):

        # iterate all peaks pairwise to find different adducts of the metabolite
        for pa in group:
            peakA = peaksInGroup[pa]

            for pb in group:
                peakB = peaksInGroup[pb]

                # peakA shall always have the lower mass
                if peakA.mz < peakB.mz:

                    ## check, if it could be a hetero atom
                    for pIso in self.heteroAtoms:
                        pIsotope = self.heteroAtoms[pIso]

                        bestFit = None
                        bestFitRatio = 100000
                        bestFitID = None

                        if abs(peakB.mz - peakA.mz - pIsotope.mzOffset) <= (
                                    max(peakA.mz, peakB.mz) * 2.5 * self.ppm / 1000000.):

                            ratio = peakB.PeakArea / peakA.PeakArea

                            for nHAtoms in range(pIsotope.minCount, pIsotope.maxCount + 1):
                                if abs(ratio - nHAtoms * pIsotope.relativeAbundance) <= self.hAIntensityError:
                                    if abs(ratio - nHAtoms * pIsotope.relativeAbundance) < bestFitRatio:
                                        bestFitRatio = abs(ratio - nHAtoms * pIsotope.relativeAbundance)
                                        bestFit = Bunch(pIso=pIso, nHAtoms=nHAtoms)
                                        bestFitID = pb

                        if bestFit is not None:
                            peakB = peaksInGroup[bestFitID]
                            # peakA.heteroAtomsFeaturePairs.append(Bunch(name="M_%s%d"%(pIso, bestFit),  partnerAdd="_%s%d"%(pIso, bestFit),  toPeak=pb))
                            peakB.heteroAtomsFeaturePairs.append(
                                Bunch(name="_%s%d" % (bestFit.pIso, bestFit.nHAtoms),
                                      partnerAdd="M_%s%d" % (bestFit.pIso, bestFit.nHAtoms), toPeak=pa))







    # annotate a metabolite group (consisting of ChromPeakPair instances) with the defined
    # adducts and generate possible in-source fragments. Remove inconsistencies
    # in the form of impossible adduct combinations (e.g. [M+H]+ and [M+Na]+ for the same ion)
    def annotateChromPeaks(self, group, peaksInGroup):


        for pe in group:
            peak = peaksInGroup[pe]

            if not hasattr(peak, "fDesc"):
                setattr(peak, "fDesc", [])
            if not hasattr(peak, "adducts"):
                setattr(peak, "adducts", [])
            if not hasattr(peak, "Ms"):
                setattr(peak, "Ms", [])


        if len(group)<=40:
            self.annotateFeaturePairsWithHeteroAtoms(group, peaksInGroup)
            for pa in group:
                peakA = peaksInGroup[pa]
                peakA.adducts.extend(peakA.heteroAtomsFeaturePairs)

            fT = formulaTools()

            # prepare adducts list
            addPol = {}
            addPol['+'] = []
            addPol['-'] = []
            adductsDict={}
            for adduct in self.adducts:
                adductsDict[adduct.name]=adduct
                if adduct.polarity != "":
                    addPol[adduct.polarity].append(adduct)
            adducts = addPol
            adducts['+']=sorted(adducts['+'], key=lambda x:x.mzoffset)
            adducts['-']=sorted(adducts['-'], key=lambda x:x.mzoffset)

            # 1. iterate all peaks pairwise to find different adducts of the metabolite
            for pa in group:
                peakA = peaksInGroup[pa]
                for pb in group:
                    peakB = peaksInGroup[pb]

                    # peakA shall always have the lower mass
                    if peakA.mz < peakB.mz and (len(peakB.heteroAtomsFeaturePairs)==0 or not all([s.name.startswith("_") for s in peakB.heteroAtomsFeaturePairs])):

                        # check, if it could be some kind of adduct combination
                        for adA in adducts[peakA.ionMode]:
                            for adB in adducts[peakB.ionMode]:
                                if adA.mCount == 1 and adB.mCount == 1 and adA.mzoffset < adB.mzoffset:

                                    if peakA.ionMode == '-' and peakB.ionMode == '+' and \
                                            peakA.loading == peakB.loading and \
                                            abs(abs(adB.mzoffset - adA.mzoffset) - 2 * 1.007276) <= (max(peakA.mz, peakB.mz) * 2.5 * self.ppm / 1000000.) and \
                                            abs(peakB.mz - peakA.mz - 2 * 1.007276 / peakA.loading) <= (max(peakA.mz, peakB.mz) * 2.5 * self.ppm / 1000000.):
                                        peakA.adducts.append(Bunch(name="[M-H]-", partnerAdd="[M+H]+", toPeak=pb))
                                        peakB.adducts.append(Bunch(name="[M+H]+", partnerAdd="[M-H]-", toPeak=pa))

                                    elif peakA.loading == adA.charge and peakB.loading == adB.charge and \
                                            abs((peakA.mz - adA.mzoffset) / adA.mCount * peakA.loading - (peakB.mz - adB.mzoffset) / adB.mCount * peakB.loading) <= (max(peakA.mz, peakB.mz) * 2.5 * self.ppm / 1000000.):
                                        peakA.adducts.append(Bunch(name=adA.name, partnerAdd=adB.name, toPeak=pb))
                                        peakB.adducts.append(Bunch(name=adB.name, partnerAdd=adA.name, toPeak=pa))

                        for adA in adducts[peakA.ionMode]:
                            for adB in adducts[peakB.ionMode]:

                                if adA.mCount == 1 and adB.mCount == 2:
                                    if peakA.loading == adA.charge and peakB.loading == adB.charge and \
                                            abs((peakA.mz - adA.mzoffset) / adA.mCount * peakA.loading - (peakB.mz - adB.mzoffset) / adB.mCount * peakB.loading) <= (max(peakA.mz, peakB.mz) * 2.5 * self.ppm / 1000000.):
                                        peakA.adducts.append(Bunch(name=adA.name, partnerAdd=adB.name, toPeak=pb))
                                        peakB.adducts.append(Bunch(name=adB.name, partnerAdd=adA.name, toPeak=pa))


            def isAdductPairPresent(pA, adductAName, adductsA, pB, adductBName, adductsB, checkPartners=True):
                aFound=False
                for add in adductsA:
                    if add.name == adductAName and (not checkPartners or (add.partnerAdd == adductBName and add.toPeak == pB)):
                        aFound=True
                bFound=False
                for add in adductsB:
                    if add.name == adductBName and (not checkPartners or (add.partnerAdd == adductAName and add.toPeak == pA)):
                        bFound=True
                return aFound and bFound

            def removeAdductFromFeaturePair(pA, adductAName, adductsA):
                aFound=[]
                for i, add in enumerate(adductsA):
                    if add.name == adductAName:
                        aFound.append(i)

                aFound=reversed(sorted(list(set(aFound))))
                for i in aFound:
                    del adductsA[i]

            def removeAdductPair(pA, adductAName, adductsA, pB, adductBName, adductsB):
                aFound=[]
                for i, add in enumerate(adductsA):
                    if add.name == adductAName and add.partnerAdd == adductBName and add.toPeak == pB:
                        aFound.append(i)
                bFound=[]
                for i, add in enumerate(adductsB):
                    if add.name == adductBName and add.partnerAdd == adductAName and add.toPeak == pA:
                        bFound.append(i)

                aFound=reversed(sorted(list(set(aFound))))
                for i in aFound:
                    del adductsA[i]
                bFound=reversed(sorted(list(set(bFound))))
                for i in bFound:
                    del adductsB[i]

            # 2. remove incorrect annotations from feature pairs e.g. A: [M+H]+ and [M+Na]+ with B: [M+Na]+ and [M+2Na-H]+
            for p in group:
                peak = peaksInGroup[p]
                peak.adducts = list(set(peak.adducts))

            for pa in group:
                peakA = peaksInGroup[pa]
                if len(peakA.adducts)>0:

                    for pb in group:
                        peakB = peaksInGroup[pb]
                        if len(peakB.adducts)>0:
                            if peakA.mz < peakB.mz:

                                for pc in group:
                                    peakC = peaksInGroup[pc]
                                    if len(peakC.adducts)>0:

                                        if      isAdductPairPresent(pa, "[M+H]+",     peakA.adducts,     pb, "[M+Na]+",        peakB.adducts) and \
                                                isAdductPairPresent(pa, "[M+Na]+",    peakA.adducts,     pb, "[M+2Na-H]+",     peakB.adducts) and \
                                                isAdductPairPresent(pa, "[M+H]+",     peakA.adducts,     pc, "[M+2Na-H]+",     peakC.adducts) and \
                                                isAdductPairPresent(pb, "[M+H]+",     peakB.adducts,     pc, "[M+Na]+",        peakC.adducts) and \
                                                isAdductPairPresent(pb, "[M+Na]+",    peakB.adducts,     pc, "[M+2Na-H]+",     peakC.adducts) and \
                                                peakA.mz < peakB.mz < peakC.mz:
                                            removeAdductPair(pa, "[M+Na]+",    peakA.adducts,     pb, "[M+2Na-H]+",     peakB.adducts)
                                            removeAdductPair(pb, "[M+H]+",     peakB.adducts,     pc, "[M+Na]+",        peakC.adducts)

                                        if      isAdductPairPresent(pa, "[M+H]+",     peakA.adducts,     pb, "[M+K]+",         peakB.adducts) and \
                                                isAdductPairPresent(pa, "[M+K]+",     peakA.adducts,     pb, "[M+2K-H]+",      peakB.adducts) and \
                                                isAdductPairPresent(pa, "[M+H]+",     peakA.adducts,     pc, "[M+2K-H]+",      peakC.adducts) and \
                                                isAdductPairPresent(pb, "[M+H]+",     peakB.adducts,     pc, "[M+K]+",         peakC.adducts) and \
                                                isAdductPairPresent(pb, "[M+K]+",     peakB.adducts,     pc, "[M+2K-H]+",      peakC.adducts) and \
                                                peakA.mz < peakB.mz < peakC.mz:
                                            removeAdductPair(pa, "[M+K]+",    peakA.adducts,     pb, "[M+2K-H]+",     peakB.adducts)
                                            removeAdductPair(pb, "[M+H]+",    peakB.adducts,     pc, "[M+K]+",        peakC.adducts)

                                        if      isAdductPairPresent(pa, "[M+2H]++",      peakA.adducts,     pb, "[M+H+Na]++",         peakB.adducts) and \
                                                isAdductPairPresent(pa, "[M+H+Na]++",    peakA.adducts,     pb, "[M+2Na]++",          peakB.adducts) and \
                                                isAdductPairPresent(pa, "[M+H]++",       peakA.adducts,     pc, "[M+2Na]++",          peakC.adducts) and \
                                                isAdductPairPresent(pb, "[M+2H]++",      peakB.adducts,     pc, "[M+H+Na]++",         peakC.adducts) and \
                                                isAdductPairPresent(pb, "[M+H+Na]++",    peakB.adducts,     pc, "[M+2Na]++",          peakC.adducts) and \
                                                peakA.mz < peakB.mz < peakC.mz:
                                            removeAdductPair(pa, "[M+H+Na]++",    peakA.adducts,     pb, "[M+2Na]++",      peakB.adducts)
                                            removeAdductPair(pb, "[M+2H]++",      peakB.adducts,     pc, "[M+H+Na]++",     peakC.adducts)

                                        if      isAdductPairPresent(pa, "[M+2H]++",      peakA.adducts,     pb, "[M+H+K]++",         peakB.adducts) and \
                                                isAdductPairPresent(pa, "[M+H+K]++",     peakA.adducts,     pb, "[M+2K]++",          peakB.adducts) and \
                                                isAdductPairPresent(pa, "[M+H]++",       peakA.adducts,     pc, "[M+2K]++",          peakC.adducts) and \
                                                isAdductPairPresent(pb, "[M+2H]++",      peakB.adducts,     pc, "[M+H+K]++",         peakC.adducts) and \
                                                isAdductPairPresent(pb, "[M+H+K]++",     peakB.adducts,     pc, "[M+2K]++",          peakC.adducts) and \
                                                peakA.mz < peakB.mz < peakC.mz:
                                            removeAdductPair(pa, "[M+H+K]++",    peakA.adducts,     pb, "[M+2K]++",      peakB.adducts)
                                            removeAdductPair(pb, "[M+2H]++",     peakB.adducts,     pc, "[M+H+K]++",     peakC.adducts)

                                        if      isAdductPairPresent(pb, "[M+H]+",        peakB.adducts,     pc, "[M+Na]+",           peakC.adducts) and \
                                                isAdductPairPresent(pa, "[M+2Na-H]+",    peakA.adducts,     pc, "[M+CH3FeN]+",       peakC.adducts) and \
                                                peakA.mz < peakB.mz < peakC.mz:
                                            removeAdductPair(pa, "[M+2Na-H]+",    peakA.adducts,     pc, "[M+CH3FeN]+",      peakC.adducts)

                                        if      isAdductPairPresent(pa, "[M-H]-",        peakA.adducts,     pb, "[M+H]+",        peakB.adducts) and \
                                                isAdductPairPresent(pa, "[M+Na-2H]-",    peakA.adducts,     pc, "[M+Na]+",       peakC.adducts) and \
                                                peakA.mz < peakB.mz < peakC.mz:
                                            aFound=[]
                                            for i, add in enumerate(peakA.adducts):
                                                if add.name == "[M+Na-2H]-":
                                                    aFound.append(i)
                                            aFound=reversed(sorted(list(set(aFound))))
                                            for i in aFound:
                                                del peakA.adducts[i]

                                        if      isAdductPairPresent(pa, "[M-H]-",     peakA.adducts,     pb, "[M+H]+",        peakB.adducts) and \
                                                isAdductPairPresent(pa, "[M+Na]+",    peakA.adducts,     pc, "[M+Na]+",       peakC.adducts) and \
                                                peakA.mz < peakB.mz < peakC.mz:
                                            bFound=[]
                                            for i, add in enumerate(peakB.adducts):
                                                if add.name == "[M+Na-2H]-":
                                                    bFound.append(i)
                                            bFound=reversed(sorted(list(set(bFound))))
                                            for i in bFound:
                                                del peakB.adducts[i]


            for pa in group:
                peakA = peaksInGroup[pa]
                if len(peakA.adducts)>0:

                    for pb in group:
                        peakB = peaksInGroup[pb]
                        if len(peakB.adducts)>0:
                            if peakA.mz < peakB.mz:

                                if      isAdductPairPresent(pa, "[M+H]+",       peakA.adducts,     pb, "[M+Na]+",         peakB.adducts) and \
                                        isAdductPairPresent(pa, "[M+Na]+",      peakA.adducts,     pb, "[M+2Na-H]+",      peakB.adducts) and \
                                        peakA.mz < peakB.mz:
                                    removeAdductPair(pa, "[M+Na]+",    peakA.adducts,     pb, "[M+2Na-H]+",      peakB.adducts)

                                if      isAdductPairPresent(pa, "[M+H]+",      peakA.adducts,     pb, "[M+K]+",         peakB.adducts) and \
                                        isAdductPairPresent(pa, "[M+K]+",      peakA.adducts,     pb, "[M+2K-H]+",      peakB.adducts) and \
                                        peakA.mz < peakB.mz:
                                    removeAdductPair(pa, "[M+K]+",    peakA.adducts,     pb, "[M+2K-H]+",      peakB.adducts)

                                if      isAdductPairPresent(pa, "[M-H]-",          peakA.adducts,     pb, "[M+Na]+",         peakB.adducts) and \
                                        isAdductPairPresent(pa, "[M+Na-2H]-",      peakA.adducts,     pb, "[M+2Na-H]+",      peakB.adducts) and \
                                        peakA.mz < peakB.mz:
                                    removeAdductPair(pa, "[M+Na-2H]-",    peakA.adducts,     pb, "[M+2Na-H]+",      peakB.adducts)

                                if      isAdductPairPresent(pa, "[M-H]-",         peakA.adducts,     pb, "[M+K]+",         peakB.adducts) and \
                                        isAdductPairPresent(pa, "[M+K-2H]-",      peakA.adducts,     pb, "[M+2K-H]+",      peakB.adducts) and \
                                        peakA.mz < peakB.mz:
                                    removeAdductPair(pa, "[M+K-2H]-",    peakA.adducts,     pb, "[M+2K-H]+",      peakB.adducts)

                                if      isAdductPairPresent(pa, "[M-H]-",          peakA.adducts,     pb, "[M+Na]+",     peakB.adducts) and \
                                        isAdductPairPresent(pa, "[M+Na-2H]-",      peakA.adducts,     pb, "[M+H]+",      peakB.adducts) and \
                                        peakA.mz < peakB.mz:
                                    removeAdductPair(pa, "[M+Na-2H]-",    peakA.adducts,     pb, "[M+Na]+",      peakB.adducts)

                                if      isAdductPairPresent(pa, "[M-H]-",         peakA.adducts,     pb, "[M+K]+",     peakB.adducts) and \
                                        isAdductPairPresent(pa, "[M+K-2H]-",      peakA.adducts,     pb, "[M+H]+",      peakB.adducts) and \
                                        peakA.mz < peakB.mz:
                                    removeAdductPair(pa, "[M+K-2H]-",    peakA.adducts,     pb, "[M+K]+",      peakB.adducts)

                                if      isAdductPairPresent(pa, "[M+2H]++",         peakA.adducts,     pb, "[M+H+Na]++",     peakB.adducts) and \
                                        isAdductPairPresent(pa, "[M+H+Na]++",       peakA.adducts,     pb, "[M+2Na]++",      peakB.adducts) and \
                                        peakA.mz < peakB.mz:
                                    removeAdductPair(pa, "[M+H+Na]++",    peakA.adducts,     pb, "[M+2Na]++",      peakB.adducts)

                                if      isAdductPairPresent(pa, "[M+2H]++",        peakA.adducts,     pb, "[M+H+K]++",     peakB.adducts) and \
                                        isAdductPairPresent(pa, "[M+H+K]++",       peakA.adducts,     pb, "[M+2K]++",      peakB.adducts) and \
                                        peakA.mz < peakB.mz:
                                    removeAdductPair(pa, "[M+H+Na]++",    peakA.adducts,     pb, "[M+2Na]++",      peakB.adducts)



            for pa in group:
                peakA = peaksInGroup[pa]
                if len(peakA.adducts)>0:
                    if      isAdductPairPresent(pa, "[M+H]+",            peakA.adducts,     pa, "[2M+H]+",         peakA.adducts, checkPartners=False):
                        removeAdductFromFeaturePair(pa, "[M+H]+",   peakA.adducts)
                    if      isAdductPairPresent(pa, "[M+NH4]+",          peakA.adducts,     pa, "[2M+NH4]+",       peakA.adducts, checkPartners=False):
                        removeAdductFromFeaturePair(pa, "[M+NH4]+", peakA.adducts)
                    if      isAdductPairPresent(pa, "[M+Na]+",           peakA.adducts,     pa, "[2M+Na]+",        peakA.adducts, checkPartners=False):
                        removeAdductFromFeaturePair(pa, "[M+Na]+",  peakA.adducts)
                    if      isAdductPairPresent(pa, "[M+K]+",            peakA.adducts,     pa, "[2M+K]+",         peakA.adducts, checkPartners=False):
                        removeAdductFromFeaturePair(pa, "[M+K]+",   peakA.adducts)


            inSourceFragments={}

            # 3. iterate all peaks pairwise to find common in-source fragments
            for pa in group:
                peakA = peaksInGroup[pa]

                for pb in group:
                    peakB = peaksInGroup[pb]

                    # peakA shall always have the lower mass
                    if peakA.mz < peakB.mz and abs(peakA.mz - peakB.mz) <= 101:
                        mzDif = peakB.mz - peakA.mz
                        done = False

                        # generate putative in-source fragments (using the labelled carbon atoms)
                        if len(self.elements) > 0:

                            elemDictReq = {}
                            for elem in self.elements:
                                elemDictReq[elem] = [self.elements[elem].weight, self.elements[elem].numberValenzElectrons]
                            t = SGRGenerator(atoms=elemDictReq)

                            useAtoms = self.elements.keys()


                            atomsRange = []
                            for elem in self.elements.keys():
                                elem = self.elements[elem]
                                atomsRange.append([elem.minCount, elem.maxCount])

                            corrFact = abs(peakB.mz - peakA.mz)
                            if corrFact <= 1:
                                corrFact = 1.
                            pFs = t.findFormulas(mzDif, ppm=self.ppm * 2. * peakA.mz / corrFact, useAtoms=useAtoms,
                                                 atomsRange=atomsRange,
                                                 useSevenGoldenRules=False)

                            for pF in pFs:
                                if pa not in inSourceFragments.keys():
                                    inSourceFragments[pa]={}
                                if pb not in inSourceFragments[pa].keys():
                                    inSourceFragments[pa][pb]=[]


                                c = fT.parseFormula(pF)
                                mw = fT.calcMolWeight(c)
                                dif = abs(abs(peakB.mz - peakA.mz) - mw)
                                sf = fT.flatToString(c, prettyPrintWithHTMLTags=False)
                                inSourceFragments[pa][pb].append("%.4f-%s"%(peakB.mz, sf))

            if self.simplifyInSourceFragments:
                for pa in group:
                    peakA = peaksInGroup[pa]
                    peakA.fDesc = []

                    for pb in group:
                        peakB = peaksInGroup[pb]

                        for pc in group:
                            peakC = peaksInGroup[pc]

                            if peakA.mz < peakC.mz < peakB.mz:

                                if pa in inSourceFragments.keys() and pb in inSourceFragments[pa].keys() and \
                                    pc in inSourceFragments[pa].keys() and \
                                    pc in inSourceFragments.keys() and pb in inSourceFragments[pc].keys():
                                    del inSourceFragments[pa][pc]

            for pa in group:
                peakA = peaksInGroup[pa]

                if pa in inSourceFragments.keys():
                    for pb in inSourceFragments[pa].keys():
                        for inFrag in inSourceFragments[pa][pb]:
                            peakA.fDesc.append(inFrag)





            for pe in group:
                peak = peaksInGroup[pe]
                peak.fDesc=list(set(peak.fDesc))
                peak.adducts=list(set([a.name for a in peak.adducts]))

                if not hasattr(peak, "Ms"):
                    setattr(peak, "Ms", [])

                peak.Ms=[]
                for assignedAdduct in peak.adducts:
                    if assignedAdduct in adductsDict.keys():
                        peak.Ms.append((peak.mz-adductsDict[assignedAdduct].mzoffset)/adductsDict[assignedAdduct].mCount/peak.loading)
