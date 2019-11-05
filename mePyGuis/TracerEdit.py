import sys
import os
import pickle
import base64
from copy import copy, deepcopy

from PyQt4 import QtGui, QtCore

from mePyGuis.TracerEditor import Ui_Dialog
from utils import getRatio, getXCombinations
from formulaTools import formulaTools, getIsotopeMass


fT = formulaTools()

class ConfiguredTracer():
    def __init__(self, name="", labeling="[13C]15", enrichmentA=.9893, enrichmentB=.995,
                 amountA=.9, amountB=1., monoisotopicRatio=1, maxRelNegBias=70, maxRelPosBias=130,
                 id=-1, mzDelta=None):
        self.id = id

        self.name = name
        self.labeling = labeling
        self.enrichmentA = enrichmentA
        self.enrichmentB = enrichmentB
        self.amountA = amountA
        self.amountB = amountB
        self.monoisotopicRatio = monoisotopicRatio
        self.maxRelNegBias = maxRelNegBias
        self.maxRelPosBias = maxRelPosBias

        if mzDelta is None:
            self.mzDelta = fT.calcIsotopologOffsetWeight(fT.parseFormula(self.labeling))
        else:
            self.mzDelta = mzDelta

    def __str__(self):
        return "ConfiguredTracer: %s %s enrichment: %.3f %.3f amount %.1f %.1f monoisotopicRatio %.3f bias: %.1f %.1f" % (
            self.name, self.labeling, self.enrichmentA, self.enrichmentB,
            self.amountA, self.amountB, self.monoisotopicRatio, self.maxRelNegBias, self.maxRelPosBias)

class tracerEdit(QtGui.QDialog, Ui_Dialog):
    def __init__(self, parent=None, initDir=None):
        QtGui.QDialog.__init__(self, parent)
        self.setupUi(self)
        self.setWindowTitle("Tracer editor")

        self.discardButton.clicked.connect(self.dialogCan)
        self.acceptButton.clicked.connect(self.dialogFin)

        self.tNativeAmount.valueChanged.connect(self.updateRatio)
        self.tLabeledAmount.valueChanged.connect(self.updateRatio)

        self.acceptButton.setFocus(True)

    def updateRatio(self):

        elems=fT.parseFormula(str(self.tLabelingIsotopes.text()), keepOnlyMainIsotopes=True)
        if len(elems)==1:
            t=None
            for k in elems.keys():
                t=k
            natR=getRatio(self.tNativeIsotopicPurity.value()/100, elems[t], 0)*self.tNativeAmount.value()
            labR=getRatio(self.tLabeledIsotopicPurity.value()/100, elems[t], 0)*self.tLabeledAmount.value()

            rat=natR/labR
            self.tRatio.setText("Ratio is %.1f%% (Ratio is based on amount and isotopic purity)" % (100. * rat))
        else:
            rat = self.tNativeAmount.value() / self.tLabeledAmount.value()
            self.tRatio.setText("Ratio is %.1f%% (Ratio is not based on isotopic purity. Please adjust it)"%(100.*rat))
        return rat

    def dialogCan(self):
        self.reject()

    def dialogFin(self):
        self.accept()

    def setTracer(self, tracer):
        if tracer is not None:
            self.tName.setText(tracer.name)
            self.tLabelingIsotopes.setText(tracer.labeling)
            self.tNativeIsotopicPurity.setValue(tracer.enrichmentA*100)
            self.tLabeledIsotopicPurity.setValue(tracer.enrichmentB*100)
            self.tNativeAmount.setValue(tracer.amountA)
            self.tLabeledAmount.setValue(tracer.amountB)
            self.monoisotopicRatio=self.updateRatio()
            self.tMinRatio.setValue(tracer.maxRelNegBias*100)
            self.tMaxRatio.setValue(tracer.maxRelPosBias*100)

            print tracer

        self.updateRatio()

    def getTracer(self):

        x=ConfiguredTracer(str(self.tName.text()),
                           str(self.tLabelingIsotopes.text()),
                           self.tNativeIsotopicPurity.value()/100,
                           self.tLabeledIsotopicPurity.value()/100,
                           self.tNativeAmount.value(),
                           self.tLabeledAmount.value(),
                           self.updateRatio(),
                           self.tMinRatio.value()/100,
                           self.tMaxRatio.value()/100)

        return x

    def executeDialog(self):
        x = self.exec_()
        return x


if __name__ == "__main__":
    import sys


    app = QtGui.QApplication(sys.argv)
    Dialog = tracerEdit()

    Dialog.show()
    x = app.exec_()

    sys.exit(x)
    