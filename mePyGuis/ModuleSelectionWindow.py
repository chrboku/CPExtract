# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\mePyGuis\guis\ModuleSelectionWindow.ui'
#
# Created: Wed Mar 02 10:42:37 2016
#      by: PyQt4 UI code generator 4.10
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(794, 503)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/MEIcon/ressources/MEIcon.ico")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.setWindowIcon(icon)
        MainWindow.setStyleSheet(_fromUtf8("background-color: rgb(255, 255, 255);"))
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.gridLayout_2 = QtGui.QGridLayout(self.centralwidget)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.line_2 = QtGui.QFrame(self.centralwidget)
        self.line_2.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_2.setFrameShape(QtGui.QFrame.HLine)
        self.line_2.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_2.setObjectName(_fromUtf8("line_2"))
        self.gridLayout.addWidget(self.line_2, 3, 0, 1, 2)
        self.tracExtractLabel = QtGui.QLabel(self.centralwidget)
        self.tracExtractLabel.setObjectName(_fromUtf8("tracExtractLabel"))
        self.gridLayout.addWidget(self.tracExtractLabel, 2, 1, 1, 1)
        self.fragExtractLabel = QtGui.QLabel(self.centralwidget)
        self.fragExtractLabel.setObjectName(_fromUtf8("fragExtractLabel"))
        self.gridLayout.addWidget(self.fragExtractLabel, 4, 1, 1, 1)
        self.allExtractLabel = QtGui.QLabel(self.centralwidget)
        self.allExtractLabel.setStyleSheet(_fromUtf8(""))
        self.allExtractLabel.setObjectName(_fromUtf8("allExtractLabel"))
        self.gridLayout.addWidget(self.allExtractLabel, 0, 1, 1, 1)
        self.line = QtGui.QFrame(self.centralwidget)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.gridLayout.addWidget(self.line, 1, 0, 1, 2)
        self.documentationLabel = QtGui.QLabel(self.centralwidget)
        self.documentationLabel.setObjectName(_fromUtf8("documentationLabel"))
        self.gridLayout.addWidget(self.documentationLabel, 7, 1, 1, 1)
        self.line_3 = QtGui.QFrame(self.centralwidget)
        self.line_3.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_3.setFrameShape(QtGui.QFrame.HLine)
        self.line_3.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_3.setObjectName(_fromUtf8("line_3"))
        self.gridLayout.addWidget(self.line_3, 5, 0, 1, 2)
        self.allExtractIcon = QtGui.QPushButton(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(90)
        sizePolicy.setVerticalStretch(90)
        sizePolicy.setHeightForWidth(self.allExtractIcon.sizePolicy().hasHeightForWidth())
        self.allExtractIcon.setSizePolicy(sizePolicy)
        self.allExtractIcon.setMinimumSize(QtCore.QSize(90, 90))
        self.allExtractIcon.setStyleSheet(_fromUtf8("background-image: url(:/AllExtract/ressources/AllExtract.png);"))
        self.allExtractIcon.setText(_fromUtf8(""))
        self.allExtractIcon.setFlat(True)
        self.allExtractIcon.setObjectName(_fromUtf8("allExtractIcon"))
        self.gridLayout.addWidget(self.allExtractIcon, 0, 0, 1, 1)
        self.tracExtractIcon = QtGui.QPushButton(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(90)
        sizePolicy.setVerticalStretch(90)
        sizePolicy.setHeightForWidth(self.tracExtractIcon.sizePolicy().hasHeightForWidth())
        self.tracExtractIcon.setSizePolicy(sizePolicy)
        self.tracExtractIcon.setMinimumSize(QtCore.QSize(90, 90))
        self.tracExtractIcon.setStyleSheet(_fromUtf8("background-image: url(:/TracExtract/ressources/TracExtract.png);"))
        self.tracExtractIcon.setText(_fromUtf8(""))
        self.tracExtractIcon.setFlat(True)
        self.tracExtractIcon.setObjectName(_fromUtf8("tracExtractIcon"))
        self.gridLayout.addWidget(self.tracExtractIcon, 2, 0, 1, 1)
        self.fragExtractIcon = QtGui.QPushButton(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(90)
        sizePolicy.setVerticalStretch(90)
        sizePolicy.setHeightForWidth(self.fragExtractIcon.sizePolicy().hasHeightForWidth())
        self.fragExtractIcon.setSizePolicy(sizePolicy)
        self.fragExtractIcon.setMinimumSize(QtCore.QSize(90, 90))
        self.fragExtractIcon.setStyleSheet(_fromUtf8("background-image: url(:/FragExtract/ressources/FragExtract.png);"))
        self.fragExtractIcon.setText(_fromUtf8(""))
        self.fragExtractIcon.setFlat(True)
        self.fragExtractIcon.setObjectName(_fromUtf8("fragExtractIcon"))
        self.gridLayout.addWidget(self.fragExtractIcon, 4, 0, 1, 1)
        self.documentationIcon = QtGui.QPushButton(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(90)
        sizePolicy.setVerticalStretch(90)
        sizePolicy.setHeightForWidth(self.documentationIcon.sizePolicy().hasHeightForWidth())
        self.documentationIcon.setSizePolicy(sizePolicy)
        self.documentationIcon.setMinimumSize(QtCore.QSize(90, 90))
        self.documentationIcon.setStyleSheet(_fromUtf8("background-image: url(:/Documentation/ressources/Documentation.png);"))
        self.documentationIcon.setText(_fromUtf8(""))
        self.documentationIcon.setFlat(True)
        self.documentationIcon.setObjectName(_fromUtf8("documentationIcon"))
        self.gridLayout.addWidget(self.documentationIcon, 7, 0, 1, 1)
        self.line_4 = QtGui.QFrame(self.centralwidget)
        self.line_4.setFrameShape(QtGui.QFrame.HLine)
        self.line_4.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_4.setObjectName(_fromUtf8("line_4"))
        self.gridLayout.addWidget(self.line_4, 6, 0, 1, 2)
        self.gridLayout_2.addLayout(self.gridLayout, 1, 1, 1, 1)
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem, 0, 1, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem1, 2, 1, 1, 1)
        spacerItem2 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem2, 1, 0, 1, 1)
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        spacerItem3 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem3)
        self.version = QtGui.QLabel(self.centralwidget)
        self.version.setStyleSheet(_fromUtf8("color: slategrey;\n"
""))
        self.version.setText(_fromUtf8(""))
        self.version.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.version.setObjectName(_fromUtf8("version"))
        self.verticalLayout.addWidget(self.version)
        self.gridLayout_2.addLayout(self.verticalLayout, 2, 4, 1, 1)
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        spacerItem4 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem4)
        self.label = QtGui.QLabel(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setMinimumSize(QtCore.QSize(80, 80))
        self.label.setStyleSheet(_fromUtf8("image: url(:/MEIcon_Large/ressources/MEIcon_Large.png);"))
        self.label.setText(_fromUtf8(""))
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout.addWidget(self.label)
        self.verticalLayout_3.addLayout(self.horizontalLayout)
        spacerItem5 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_3.addItem(spacerItem5)
        spacerItem6 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.verticalLayout_3.addItem(spacerItem6)
        self.gridLayout_2.addLayout(self.verticalLayout_3, 0, 4, 2, 1)
        MainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow", None))
        self.tracExtractLabel.setText(_translate("MainWindow", "TracExtract", None))
        self.fragExtractLabel.setText(_translate("MainWindow", "FragExtract", None))
        self.allExtractLabel.setText(_translate("MainWindow", "AllExtract", None))
        self.documentationLabel.setText(_translate("MainWindow", "Documentation", None))

import ressources_rc

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

