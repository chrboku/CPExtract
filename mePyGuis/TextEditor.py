
from TextEditorGUI import Ui_Dialog

from PyQt4 import QtGui
from PyQt4 import QtCore



class TextEditor(QtGui.QDialog, Ui_Dialog):
    def __init__(self, parent=None, labelText=None, text=None):
        QtGui.QDialog.__init__(self, parent)
        self.setWindowTitle("Text editor")
        self.setupUi(self)

        if labelText is not None:
            self.label.setText(labelText)
        if text is not None:
            self.plainTextEdit.setPlainText(text)


    def getText(self):
        return str(self.plainTextEdit.toPlainText())





if __name__=="__main__":
    import sys

    app = QtGui.QApplication(sys.argv)
    Dialog = TextEditor()

    Dialog.show()
    x = app.exec_()

    print Dialog.getText()