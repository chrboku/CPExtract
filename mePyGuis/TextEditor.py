
from TextEditorGUI import Ui_Dialog

from PyQt4 import QtGui
from PyQt4 import QtCore



class TextEditor(QtGui.QDialog, Ui_Dialog):
    def __init__(self, parent=None, labelText=None, text=None, textChangeCallback=None):
        QtGui.QDialog.__init__(self, parent)
        self.setWindowTitle("Text editor")
        self.setupUi(self)

        if labelText is not None:
            self.label.setText(labelText)
        if text is not None:
            self.plainTextEdit.setPlainText(text)

        self.textChangeCallback=textChangeCallback
        if textChangeCallback is not None:
            self.plainTextEdit.textChanged.connect(self.textChanged)


    def getText(self):
        return str(self.plainTextEdit.toPlainText())


    def textChanged(self):
        x=self.textChangeCallback(self.getText())
        self.buttonBox.setEnabled(x[0])

if __name__=="__main__":


    def testForPythonCode(text):
        try:
            exec(text)
            return (True, "")
        except Exception as ex:
            return (False, ex.message)


    import sys

    app = QtGui.QApplication(sys.argv)
    Dialog = TextEditor(labelText="Please enter some valid Python code", text="2+2", textChangeCallback=testForPythonCode)

    Dialog.show()
    x = app.exec_()

    print Dialog.getText()