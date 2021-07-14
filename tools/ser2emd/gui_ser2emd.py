'''
GUI tool to convert SER file to EMD file.
'''


import sys
import os
from PyQt5 import QtGui, QtCore, QtWidgets

import ncempy.io.ser

class Converter(QtWidgets.QWidget):

    def __init__(self):
        super().__init__()
        
        self.initUI()
        
    def initUI(self):
        #self.statusBar().showMessage('Welcome')
        
        self.resize(600,250)
        qr = self.frameGeometry()
        cp = QtWidgets.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
        
        ser_lbl = QtWidgets.QLabel('input SER file:', self) 
        
        self.ser_txt = QtWidgets.QLineEdit(self)
        self.ser_txt.setReadOnly(True)
        self.serButton = QtWidgets.QPushButton('Open', self)
        
        hbox_ser = QtWidgets.QHBoxLayout()
        hbox_ser.addWidget(self.ser_txt)
        hbox_ser.addWidget(self.serButton)       
        
        ser_sep = QtWidgets.QFrame()
        ser_sep.setFrameStyle(QtWidgets.QFrame.HLine | QtWidgets.QFrame.Sunken)
        
        
        emd_lbl = QtWidgets.QLabel('output EMD file:', self) 
        
        self.emd_txt = QtWidgets.QLineEdit(self)
        self.emd_txt.setReadOnly(True)
        self.emdButton = QtWidgets.QPushButton('Save', self)
        
        hbox_emd = QtWidgets.QHBoxLayout()
        hbox_emd.addWidget(self.emd_txt)
        hbox_emd.addWidget(self.emdButton)
        
        emd_sep = QtWidgets.QFrame()
        emd_sep.setFrameStyle(QtWidgets.QFrame.HLine | QtWidgets.QFrame.Sunken)
        
        
        self.msg = QtWidgets.QLabel('Ready', self)
        self.convButton = QtWidgets.QPushButton('Convert', self)
        self.exitButton = QtWidgets.QPushButton('Exit', self)
        
        hbox_buttons = QtWidgets.QHBoxLayout()
        hbox_buttons.addWidget(self.msg)
        hbox_buttons.addStretch(1)
        hbox_buttons.addWidget(self.convButton)
        hbox_buttons.addWidget(self.exitButton)
  
        
        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(ser_lbl)
        vbox.addLayout(hbox_ser)
        vbox.addWidget(ser_sep)
        vbox.addWidget(emd_lbl)
        vbox.addLayout(hbox_emd)
        vbox.addWidget(emd_sep)
        vbox.addStretch(1)
        vbox.addLayout(hbox_buttons)
        
        self.setLayout(vbox)
        
        self.setWindowTitle('Convert SER to EMD')
        self.show()
        
        
        self.serButton.clicked.connect(self.clicked_serButton)
        self.emdButton.clicked.connect(self.clicked_emdButton)
        
        self.convButton.clicked.connect(self.convert)
        self.exitButton.clicked.connect(QtCore.QCoreApplication.instance().quit)
        
        
        
    def keyPressEvent(self, e):
        
        # esc to close
        if e.key() == QtCore.Qt.Key_Escape:
            self.close()
            
    def clicked_serButton(self):
        self.msg.setText('Ready')
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open SER file', filter='SER files (*.ser);;All files (*.*)')
        if type(fname)==list:
            fname = fname[0]
        self.ser_txt.setText(fname)
       
    
    def clicked_emdButton(self):
        self.msg.setText('Ready')
        fname = QtWidgets.QFileDialog.getSaveFileName(self, 'Save EMD file')
        if type(fname)==list:
            fname = fname[0]
        self.emd_txt.setText(fname)
        
    def convert(self):
        ser_name = self.ser_txt.text()
            
        fser = ncempy.io.ser.fileSER(ser_name)
        
        if os.path.isfile(self.emd_txt.text()):
            os.remove(self.emd_txt.text())
        
        fser.writeEMD(self.emd_txt.text())
        
        self.msg.setText('Done')
            
    
    
        
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    con = Converter()
    sys.exit(app.exec_())
