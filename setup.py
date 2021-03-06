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



from distutils.core import setup
import py2exe
from shutil import rmtree, copy, move
import os
import sys
import zipfile

from MetExtractII_Main import MetExtractVersion

import matplotlib

import sys
sys.path.append("../PyMassBankSearchTool")
sys.setrecursionlimit(5000)
from TableUtils import TableUtils








def replaceInFile(textToReplace, newText, filein, fileout=None):
    f = open(filein,'r')
    filedata = f.read()
    f.close()

    newdata = filedata.replace(textToReplace, newText)

    if fileout is None:
        fileout=filein
    f = open(fileout,'w')
    f.write(newdata)
    f.close()


########################################################################################################################
########################################################################################################################
########################################################################################################################



#<editor-fold desc="### check if R is installed and accessible">
def checkR():
    try:
        import rpy2.robjects as ro              # import RPy2 module
        r = ro.r                                # make R globally accessible

        v = r("R.Version()$version.string")     # if this R-commando is executed, the RPy2 connection to the
                                                # R subprocess has been established
        return True
    except:
        # The R subprocess could not be started / accessed successfully
        return False

def loadRConfFile(path):
    import os
    if os.path.isfile(path+"/RPATH.conf"):
        with open(path+"/RPATH.conf", "rb") as rconf:
            line=rconf.readline()
            os.environ["R_HOME"]=line
            return True
    else:
        return False


__RHOMEENVVAR=""
import os
def get_main_dir():
    from utils import get_main_dir
    return os.path.join(get_main_dir(), '')
if "R_HOME" in os.environ.keys():
    __RHOMEENVVAR=os.environ["R_HOME"]

os.environ["R_USER"]=get_main_dir()+"/Ruser"
# try to load r configuration file (does not require any environment variables or registry keys)
if not loadRConfFile(path=get_main_dir()) or not checkR():
    os.environ["R_HOME"]=get_main_dir()+"/R"

    if checkR():
        with open("RPATH.conf", "wb") as rconf:
            rconf.write(get_main_dir()+"/R")
            tryLoad=False

            # Show a dialog box to the user that R could not be started
            from os import sys
            from PyQt4 import QtGui, QtCore

            app = QtGui.QApplication(sys.argv)

            QtGui.QMessageBox.information(None, "MetExtract",
                      "R successfully configured\nUsing MetExtract R-Installation\nPlease restart",
                      QtGui.QMessageBox.Ok)
            sys.exit(0)
    else:

        os.environ["R_HOME"]=__RHOMEENVVAR
        os.environ["R_HOME_FROM"]="RPATH environment variable"
        if not checkR():

            print "Error: R could not be loaded correctly (No RPATH.conf file or R_HOME environment variable found)\nPlease make sure it is installed and accessible"

            # Show a dialog box to the user that R could not be started
            from os import sys
            from PyQt4 import QtGui, QtCore

            app = QtGui.QApplication(sys.argv)

            if QtGui.QMessageBox.warning(None, "MetExtract",
                                      "Error: R could not be loaded\nPlease make sure it is installed and accessible\n"
                                      "The default installation path is C:\\Program Files\\R\n"
                                      "Do you want to specify the folder?",
                                      QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)==QtGui.QMessageBox.Yes:
                tryLoad=True
                from utils import get_main_dir
                lastDir=get_main_dir()
                while tryLoad:
                    folder = str(QtGui.QFileDialog.getExistingDirectory(None, "Select R-directory (not bin folder)", directory=lastDir))
                    if folder=="":
                        sys.exit(1)
                    else:
                        lastDir=folder
                        os.environ["R_HOME"]=folder
                        if checkR():
                            with open("RPATH.conf", "wb") as rconf:
                                rconf.write(folder)
                                tryLoad=False

                                QtGui.QMessageBox.information(None, "MetExtract",
                                          "R successfully configured\nPlease restart",
                                          QtGui.QMessageBox.Ok)
                                sys.exit(0)
                        else:
                            if QtGui.QMessageBox.warning(None, "MetExtract",
                                          "Error: R could not be loaded from the specified location\n"
                                          "%s\n\n"
                                          "Please make sure it is installed and accessible\n"
                                          "The default installation path is C:\\Program Files\\R\n"
                                          "Do you want to specify the folder?"%folder,
                                          QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)==QtGui.QMessageBox.Yes:
                                pass
                            else:
                                sys.exit(1)
            else:
                sys.exit(1)
else:
    os.environ["R_HOME_FROM"]="RPATH.conf of MetExtract II"
#</editor-fold>
#Used to locate R.dll (no idea why this is necessary)
import rpy2.robjects as ro
r = ro.r




# remove previous setup files
try:
    rmtree("./dist/")
except:
    pass
try:
    rmtree("./build/")
except:
    pass
try:
    rmtree("./distribute/")
except:
    print "Error: could not clean up py2exe environment prior to compilation\n==============================\n"

# get local files (images, R-Scripts, ...)
data_files = matplotlib.get_py2exe_datafiles()
data_files.append("./chromPeakPicking/MassSpecWaveletIdentification.r")
data_files.append("./XICAlignment.r")

err = False

import os
import shutil

def mergeDirs (root_src_dir, root_dst_dir):
    for src_dir, dirs, files in os.walk(root_src_dir):
        dst_dir = src_dir.replace(root_src_dir, root_dst_dir)
        if not os.path.exists(dst_dir):
            os.mkdir(dst_dir)
        for file_ in files:
            src_file = os.path.join(src_dir, file_)
            dst_file = os.path.join(dst_dir, file_)
            if os.path.exists(dst_file):
                os.remove(dst_file)
            shutil.move(src_file, dst_dir)
import openpyxl


class Target:
    def __init__(self, **kw):
        self.__dict__.update(kw)

c=Target(script = "CPExtract.py")

print "###################################################"
print "########## Packing MetExtractII_Main"
print "###################################################"
setup(console=[c],
      options={"py2exe": {
                 "includes": ["sip", "matplotlib.backends.backend_tkagg", 'scipy', 'scipy.integrate', 'scipy.special.*','scipy.linalg.*', 'scipy.sparse.csgraph._validation', 'scipy._lib.messagestream'],  # use this line if above does not work
                 "dll_excludes": ["MSVCP90.dll", "api-ms-win-core-string-l1-1-0.dll","api-ms-win-core-registry-l1-1-0.dll","api-ms-win-core-errorhandling-l1-1-0.dll","api-ms-win-core-string-l2-1-0.dll",
                                  "api-ms-win-core-profile-l1-1-0.dll","api-ms-win*.dll","api-ms-win-core-processthreads-l1-1-2.dll","api-ms-win-core-libraryloader-l1-2-1.dll","api-ms-win-core-file-l1-2-1.dll",
                                  "api-ms-win-security-base-l1-2-0.dll","api-ms-win-eventing-provider-l1-1-0.dll","api-ms-win-core-heap-l2-1-0.dll","api-ms-win-core-libraryloader-l1-2-0.dll","api-ms-win-core-localization-l1-2-1.dll",
                                  "api-ms-win-core-sysinfo-l1-1-0.dll","api-ms-win-core-synch-l1-2-0.dll","api-ms-win-core-heap-l1-2-0.dll","api-ms-win-core-handle-l1-1-0.dll","api-ms-win-core-io-l1-1-1.dll","api-ms-win-core-com-l1-1-1.dll",
                                  "api-ms-win-core-memory-l1-1-2.dll","api-ms-win-core-version-l1-1-1.dll","api-ms-win-core-version-l1-1-0.dll","api-ms-win-core-processthreads-l1-1-0.dll"],
                 "excludes": ["_gtkagg", "_tkagg", 'jinja2.asyncsupport','jinja2.asyncfilters'],
                 "packages": ["FileDialog", "openpyxl", 'reportlab','reportlab.graphics.charts','reportlab.graphics.samples','reportlab.graphics.widgets','reportlab.graphics.barcode','reportlab.graphics','reportlab.lib','reportlab.pdfbase','reportlab.pdfgen','reportlab.platypus', 'zeep', 'lxml', 'scipy'],
                 'dist_dir': "./dist"
      }},
      data_files=data_files,
      requires=['matplotlib'])
#rmtree("./build/")


print "forced waiting..."
import time
time.sleep(3)


print "Setup finished\n==============================\n"

# copy Settings and necessary files to distribution folder
try:
    os.makedirs("./dist/Settings/")
    copy("./Settings/defaultSettings.ini", "./dist/Settings/defaultSettings.ini")

    os.makedirs("./dist/chromPeakPicking/")
    copy("./chromPeakPicking/MassSpecWaveletIdentification.r", "./dist/chromPeakPicking/MassSpecWaveletIdentification.r")
    copy("./XICAlignment.r", "./dist/XICAlignment.r")
    copy("./LICENSE.txt", "./dist/LICENSE.txt")
    copy("./calculateIsotopeEnrichment.R", "./dist/calculateIsotopeEnrichment.R")
    print "Additional resources copied\n==============================\n"
except:
    print "Error: Could not copy all required files"
    sys.exit(1)
    err = True


# copy documentation
def copyAllFilesInFolder(src, dest):
    src_files = os.listdir(src)
    for file_name in src_files:
        full_file_name = os.path.join(src, file_name)
        if (os.path.isfile(full_file_name)):
            shutil.copy(full_file_name, dest)


# rename dist folder to PyMetExtract and current version of the software
meDistFolder="./dist"
try:
    from time import sleep  # sometimes, the re-naming does not work (probably some kind of lock from the OS)
    sleep(3)                # this short waiting time decreases the number of times the renaming does not work
    os.rename("./dist", "./CPExtract_%s"%MetExtractVersion)
    meDistFolder="./CPExtract_%s"%MetExtractVersion
    print "Distribution renamed\n==============================\n"

except:
    print "Error: Could not rename dist folder"
    err = True

try:
    os.makedirs('./distribute')
except:
    pass

import os

## update MetExtract Version in NSIS setup file
replaceInFile("$$CPEXTRACTVERSION$$", MetExtractVersion, filein="setup.nsi", fileout="setup_curVersion.nsi")

os.system("\"c:\\Program Files (x86)\\NSIS\\makensis.exe\" setup_curVersion.nsi")
os.remove("./setup_curVersion.nsi")
os.rename("./Setup.exe", "./Setup_CPExtract_%s.exe"%MetExtractVersion)
move("./Setup_CPExtract_%s.exe"%MetExtractVersion, "./distribute/Setup_CPExtractII_%s.exe"%MetExtractVersion)
#shutil.copy("../../distribution/vcredist_x86.exe", "./distribute/vcredist_x86.exe")


def zipdir(path, zip):
    for root, dirs, files in os.walk(path):
        for file in files:
            zip.write(os.path.join(root, file))

# create zip archive from the executables and documentation
if not err:

    print "Zipping MetExtract (%s)" % MetExtractVersion
    zipFileName = "PyCPExtract_%s.zip" % MetExtractVersion

    zipF = zipfile.ZipFile(zipFileName, 'w')
    zipdir(meDistFolder, zipF)
    zipF.close()

    move('./%s' % zipFileName, './distribute/%s' % zipFileName)

    print "CPExtract (%s) created\n see %s\n==============================\n" % (
        MetExtractVersion, './distribute/%s' % zipFileName)

    try:
        rmtree(meDistFolder)
    except:
        print "Cleanup failed. dist and/or build directories still there\n==============================\n"

