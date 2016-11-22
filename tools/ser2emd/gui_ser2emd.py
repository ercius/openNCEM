'''
GUI tool to convert SER file to EMD file.
'''


import sys
import os
from PySide import QtGui, QtCore

import ncempy.fio.ser

class Converter(QtGui.QWidget):

    def __init__(self):
        super().__init__()
        
        self.initUI()
        
    def initUI(self):
        #self.statusBar().showMessage('Welcome')
        
        self.resize(600,250)
        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
        
        ser_lbl = QtGui.QLabel('input SER file:', self) 
        
        self.ser_txt = QtGui.QLineEdit(self)
        self.ser_txt.setReadOnly(True)
        self.serButton = QtGui.QPushButton('Open', self)
        
        hbox_ser = QtGui.QHBoxLayout()
        hbox_ser.addWidget(self.ser_txt)
        hbox_ser.addWidget(self.serButton)       
        
        ser_sep = QtGui.QFrame()
        ser_sep.setFrameStyle(QtGui.QFrame.HLine | QtGui.QFrame.Sunken)
        
        emi_lbl = QtGui.QLabel('input EMI file: (optional)', self) 
        
        self.emi_txt = QtGui.QLineEdit(self)
        self.emi_txt.setReadOnly(True)
        self.emiButton = QtGui.QPushButton('Open', self)
        
        hbox_emi = QtGui.QHBoxLayout()
        hbox_emi.addWidget(self.emi_txt)
        hbox_emi.addWidget(self.emiButton)
        
        emi_sep = QtGui.QFrame()
        emi_sep.setFrameStyle(QtGui.QFrame.HLine | QtGui.QFrame.Sunken)
        
        emd_lbl = QtGui.QLabel('output SER file:', self) 
        
        self.emd_txt = QtGui.QLineEdit(self)
        self.emd_txt.setReadOnly(True)
        self.emdButton = QtGui.QPushButton('Save', self)
        
        hbox_emd = QtGui.QHBoxLayout()
        hbox_emd.addWidget(self.emd_txt)
        hbox_emd.addWidget(self.emdButton)
        
        emd_sep = QtGui.QFrame()
        emd_sep.setFrameStyle(QtGui.QFrame.HLine | QtGui.QFrame.Sunken)
        
        
        self.msg = QtGui.QLabel('Ready', self)
        self.convButton = QtGui.QPushButton('Convert', self)
        self.exitButton = QtGui.QPushButton('Exit', self)
        
        hbox_buttons = QtGui.QHBoxLayout()
        hbox_buttons.addWidget(self.msg)
        hbox_buttons.addStretch(1)
        hbox_buttons.addWidget(self.convButton)
        hbox_buttons.addWidget(self.exitButton)
  
        
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(ser_lbl)
        vbox.addLayout(hbox_ser)
        vbox.addWidget(ser_sep)
        vbox.addWidget(emi_lbl)
        vbox.addLayout(hbox_emi)
        vbox.addWidget(emi_sep)       
        vbox.addWidget(emd_lbl)
        vbox.addLayout(hbox_emd)
        vbox.addWidget(emd_sep)
        vbox.addStretch(1)
        vbox.addLayout(hbox_buttons)
        
        self.setLayout(vbox)
        
        self.setWindowTitle('Convert SER to EMD')
        self.show()
        
        
        self.serButton.clicked.connect(self.clicked_serButton)
        self.emiButton.clicked.connect(self.clicked_emiButton)
        self.emdButton.clicked.connect(self.clicked_emdButton)
        
        self.convButton.clicked.connect(self.convert)
        self.exitButton.clicked.connect(QtCore.QCoreApplication.instance().quit)
        
        
        
    def keyPressEvent(self, e):
        
        # esc to close
        if e.key() == QtCore.Qt.Key_Escape:
            self.close()
            
    def clicked_serButton(self):
        self.msg.setText('Ready')
        fname, _ = QtGui.QFileDialog.getOpenFileName(self, 'Open SER file', filter='SER files (*.ser);;All files (*.*)')
        self.ser_txt.setText(fname)
        
    def clicked_emiButton(self):
        self.msg.setText('Ready')
        fname, _ = QtGui.QFileDialog.getOpenFileName(self, 'Open EMI file', filter='EMI files (*.emi);;All files (*.*)')    
        self.emi_txt.setText(fname)
    
    def clicked_emdButton(self):
        self.msg.setText('Ready')
        fname, _ = QtGui.QFileDialog.getSaveFileName(self, 'Save EMD file')
        self.emd_txt.setText(fname)
        
    def convert(self):
        ser_name = self.ser_txt.text()
        if not self.emi_txt.text() == '':
            emi_name = self.emi_txt.text() 
        else:
            emi_name = None
            
        fser = ncempy.fio.ser.fileSER(ser_name, emifile=emi_name)
        
        if os.path.isfile(self.emd_txt.text()):
            os.remove(self.emd_txt.text())
        
        fser.writeEMD(self.emd_txt.text())
        
        self.msg.setText('Done')
            
    
    
        
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    con = Converter()
    sys.exit(app.exec_())
