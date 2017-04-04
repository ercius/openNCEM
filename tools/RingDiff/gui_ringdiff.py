'''
GUI tool to evaluate ring diffraction patterns.
'''


import sys
import numpy as np
import os
import copy
import h5py

import ncempy.io.emd
import ncempy.algo.local_max
import ncempy.algo.distortion
import ncempy.algo.radial_profile
import ncempy.eval.ring_diff

from PyQt4 import QtGui, QtCore

import pyqtgraph as pg

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')


class MaskDialog(QtGui.QDialog):

    def __init__(self, parent):
        super().__init__(parent=parent)
        
        self.parent=parent
        
        self.setModal(True)
        
        self.imgview = pg.ImageView()
        self.imgview.setImage(parent.data[parent.idx,:,:].transpose())
        self.imgview.ui.roiBtn.hide()
        self.imgview.ui.menuBtn.hide()
        
        self.roi = pg.PolyLineROI([ [int(0.25*parent.data.shape[2]), int(0.25*parent.data.shape[1])], [int(0.75*parent.data.shape[2]), int(0.25*parent.data.shape[1])], [int(0.75*parent.data.shape[2]), int(0.75*parent.data.shape[1])], [int(0.25*parent.data.shape[2]), int(0.75*parent.data.shape[1])]], closed=True)
        self.imgview.getView().addItem(self.roi)
        
        okbtn = QtGui.QPushButton('Use new mask', self)
        okbtn.clicked.connect(self.accept)
        okbtn.setDefault(True)
        okbtn.setToolTip('Create mask to remove selected ROI from affecting the radial profile.')
        
        cclbtn = QtGui.QPushButton('Remove mask', self)
        cclbtn.clicked.connect(self.reject)
        cclbtn.setToolTip('Cancel procedure and remove any mask present in the settings.')

        vbox = QtGui.QVBoxLayout(self)
        vbox.addWidget(self.imgview)
        hbox_btns = QtGui.QHBoxLayout()
        hbox_btns.addStretch(1)
        hbox_btns.addWidget(okbtn)
        hbox_btns.addWidget(cclbtn)
        vbox.addLayout(hbox_btns)
        
        self.setGeometry(100,100,768,768+64)
        self.setWindowTitle('Create Mask')
        
        self.setToolTip('This dialog allows to create a mask.\nUse the ROI to select the area, which should be removed from affecting the radial profile.')
        
        self.show()
        
        
    def getResult(self):
        
        im = QtGui.QImage(self.parent.data[self.parent.idx,:,:].shape[0], self.parent.data[self.parent.idx,:,:,].shape[1], QtGui.QImage.Format_ARGB32) 
        im.fill(1)
                
        p = QtGui.QPainter(im)
        p.setPen(pg.functions.mkPen(None))
        p.setBrush(pg.functions.mkBrush('k'))
        p.setTransform(self.roi.itemTransform(self.imgview.getImageItem())[0])
        p.drawPath(self.roi.shape())
        p.end()
        
        return pg.functions.imageToArray(im)[:,:,0].astype(int).transpose()
        
        
    @staticmethod
    def getMask(parent):
        dialog = MaskDialog(parent)
        result = dialog.exec_()
        mask = dialog.getResult()
        return (mask, result==QtGui.QDialog.Accepted) 



class Main(QtGui.QMainWindow):

    def __init__(self):
        super().__init__()
        
        self.gui_file = {}
        self.femd_in = None
        self.gui_localmax = {}
        self.gui_polar = {}
        self.gui_radprof = {}
        self.gui_run = {}
        self.log_wdg = None
        
        self.plt_localmax = None
        self.plt_localmax_img = None
        self.plt_polar = None
        self.plt_radprof = None
        
        self.reset()
        
        self.initUI()
        
        self.log('This is the "Ring Diffraction Pattern Evaluation" tool, part of the OpenNCEM package.')
        self.log('----------')
        self.log('Ready to go ..')
        
        
        
    def reset(self):
        self.data = None
        self.dims = None
        self.settings = {}
        self.points = None
        self.center = None
        self.dists = None
        self.radprof = None
        self.back_params = None
        self.res = None
        self.mask = None
        self.idx = None
    
    
    def initUI(self):
    
        self.statusBar().showMessage('')
        
        self.mnwid = QtGui.QWidget()
        
        self.setCentralWidget(self.mnwid)
        
        self.log_wdg = QtGui.QWidget()
        self.log_txt = QtGui.QPlainTextEdit(self.log_wdg)
        self.log_txt.setReadOnly(True)
        hbox_log = QtGui.QHBoxLayout(self.log_wdg)
        hbox_log.addWidget(self.log_txt)
        self.log_wdg.setGeometry(300,300,768,200)
        self.log_wdg.setWindowTitle('Log')
        self.log_wdg.show()
        
        ## file informations
        frame_files = QtGui.QGroupBox('Input', self.mnwid)
        layout_files = QtGui.QVBoxLayout(frame_files)
        
        self.gui_file['re_btn'] = QtGui.QPushButton('Reopen Evaluation File')
        self.gui_file['re_btn'].clicked.connect(self.on_reopen)
        layout_files.addWidget(self.gui_file['re_btn'])

        self.gui_file['in_btn'] = QtGui.QPushButton('Open', frame_files)
        self.gui_file['in_btn'].clicked.connect(self.on_open_infile)
        self.gui_file['in_lbl'] = QtGui.QLabel( 'input emd file: ', frame_files )
        self.gui_file['in_txt'] = QtGui.QLineEdit( '', frame_files )
        self.gui_file['in_txt'].setReadOnly(True)
        hbox_infile = QtGui.QHBoxLayout()
        hbox_infile.addWidget(self.gui_file['in_lbl'])        
        hbox_infile.addWidget(self.gui_file['in_txt'])
        hbox_infile.addWidget(self.gui_file['in_btn'])
        layout_files.addLayout(hbox_infile)
        
        self.gui_file['idx_lbl'] = QtGui.QLabel('index: ', frame_files)
        self.gui_file['idx_slider'] = QtGui.QSlider(QtCore.Qt.Horizontal, frame_files)
        self.gui_file['idx_slider'].setMinimum(0)
        self.gui_file['idx_slider'].setMaximum(0)
        self.gui_file['idx_slider'].valueChanged.connect(self.on_indexSlider)
        self.gui_file['idx_value'] = QtGui.QLabel('0', frame_files)
        hbox_idx = QtGui.QHBoxLayout()
        hbox_idx.addWidget(self.gui_file['idx_lbl'])
        hbox_idx.addWidget(self.gui_file['idx_slider'])
        hbox_idx.addWidget(self.gui_file['idx_value'])
        layout_files.addLayout(hbox_idx)
        

        ## local maxima stuff
        frame_localmax = QtGui.QGroupBox('Local Maxima', self.mnwid)
        layout_localmax = QtGui.QVBoxLayout(frame_localmax)
        
        self.gui_localmax['lbl_lmax_r'] = QtGui.QLabel('local radius: ', frame_localmax)
        self.gui_localmax['txt_lmax_r'] = QtGui.QLineEdit( '', frame_localmax)
        hbox_lmax_r = QtGui.QHBoxLayout()
        hbox_lmax_r.addWidget(self.gui_localmax['lbl_lmax_r'])
        hbox_lmax_r.addWidget(self.gui_localmax['txt_lmax_r'])
        layout_localmax.addLayout(hbox_lmax_r)
        
        self.gui_localmax['lbl_lmax_thresh'] = QtGui.QLabel('threshold: ', frame_localmax)
        self.gui_localmax['txt_lmax_thresh'] = QtGui.QLineEdit( '', frame_localmax)
        hbox_lmax_thresh = QtGui.QHBoxLayout()
        hbox_lmax_thresh.addWidget(self.gui_localmax['lbl_lmax_thresh'])
        hbox_lmax_thresh.addWidget(self.gui_localmax['txt_lmax_thresh'])
        layout_localmax.addLayout(hbox_lmax_thresh)
        
        self.gui_localmax['lbl_lmax_cinit'] = QtGui.QLabel('init center: ')
        self.gui_localmax['txt_lmax_cinit'] = QtGui.QLineEdit( '', frame_localmax)
        self.gui_localmax['btn_lmax_cinit'] = QtGui.QPushButton( 'Select', frame_localmax)
        self.gui_localmax['btn_lmax_cinit'].clicked.connect(self.update_localmax)
        self.gui_localmax['btn_lmax_cinit'].setCheckable(True)
        hbox_lmax_cinit = QtGui.QHBoxLayout()
        hbox_lmax_cinit.addWidget(self.gui_localmax['lbl_lmax_cinit'])
        hbox_lmax_cinit.addWidget(self.gui_localmax['txt_lmax_cinit'])
        hbox_lmax_cinit.addWidget(self.gui_localmax['btn_lmax_cinit'])
        layout_localmax.addLayout(hbox_lmax_cinit)
        
        self.gui_localmax['lbl_lmax_range'] = QtGui.QLabel('radial range: ', frame_localmax)
        self.gui_localmax['txt_lmax_range'] = QtGui.QLineEdit( '', frame_localmax)
        hbox_lmax_range = QtGui.QHBoxLayout()
        hbox_lmax_range.addWidget(self.gui_localmax['lbl_lmax_range'])
        hbox_lmax_range.addWidget(self.gui_localmax['txt_lmax_range'])
        layout_localmax.addLayout(hbox_lmax_range)
        
        self.gui_localmax['lmax_btn'] = QtGui.QPushButton('Find Local Maxima', frame_localmax)
        self.gui_localmax['lmax_btn'].clicked.connect(self.on_localmax)
        layout_localmax.addWidget(self.gui_localmax['lmax_btn'])
        
        self.gui_localmax['min_lbl'] = QtGui.QLabel('min: ', frame_localmax)
        self.gui_localmax['min_slider'] = QtGui.QSlider(QtCore.Qt.Horizontal, frame_localmax)
        self.gui_localmax['min_slider'].setMinimum(0)
        self.gui_localmax['min_slider'].setMaximum(0)
        self.gui_localmax['min_slider'].valueChanged.connect(self.on_intensitySlider)
        self.gui_localmax['min_value'] = QtGui.QLabel('0', frame_localmax)
        hbox_lmax_min = QtGui.QHBoxLayout()
        hbox_lmax_min.addWidget(self.gui_localmax['min_lbl'])
        hbox_lmax_min.addWidget(self.gui_localmax['min_slider'])
        hbox_lmax_min.addWidget(self.gui_localmax['min_value'])
        layout_localmax.addLayout(hbox_lmax_min)
        
        self.gui_localmax['max_lbl'] = QtGui.QLabel('max: ', frame_localmax)
        self.gui_localmax['max_slider'] = QtGui.QSlider(QtCore.Qt.Horizontal, frame_localmax)
        self.gui_localmax['max_slider'].setMinimum(0)
        self.gui_localmax['max_slider'].setMaximum(0)
        self.gui_localmax['max_slider'].valueChanged.connect(self.on_intensitySlider)
        self.gui_localmax['max_value'] = QtGui.QLabel('0', frame_localmax)
        hbox_lmax_max = QtGui.QHBoxLayout()
        hbox_lmax_max.addWidget(self.gui_localmax['max_lbl'])
        hbox_lmax_max.addWidget(self.gui_localmax['max_slider'])
        hbox_lmax_max.addWidget(self.gui_localmax['max_value'])
        layout_localmax.addLayout(hbox_lmax_max)
        
        #self.gui_localmax['upd_btn'] = QtGui.QPushButton('Update', frame_localmax)
        #self.gui_localmax['upd_btn'].clicked.connect(self.update_localmax)
        #layout_localmax.addWidget(self.gui_localmax['upd_btn'])


        ## polar plot stuff
        frame_polar = QtGui.QGroupBox('Polar Plot',self.mnwid)
        layout_polar = QtGui.QVBoxLayout(frame_polar)
        
        self.gui_polar['center_lbl'] = QtGui.QLabel('center: ', frame_polar)
        self.gui_polar['cencpy_btn'] = QtGui.QPushButton('Copy Init', frame_polar)
        self.gui_polar['cencpy_btn'].clicked.connect(self.on_copyCenter)
        hbox_polar_center = QtGui.QHBoxLayout()
        hbox_polar_center.addWidget(self.gui_polar['center_lbl'])
        hbox_polar_center.addWidget(self.gui_polar['cencpy_btn'])
        layout_polar.addLayout(hbox_polar_center)
        
        self.gui_polar['cenopt_btn'] = QtGui.QPushButton('Optimize Center', frame_polar)
        self.gui_polar['cenopt_btn'].clicked.connect(self.on_optimizeCenter)
        layout_polar.addWidget(self.gui_polar['cenopt_btn'])
        
        self.gui_polar['dist_lbl'] = QtGui.QLabel('distortion orders: ', frame_polar)
        self.gui_polar['dist_txt'] = QtGui.QLineEdit('', frame_polar)
        hbox_polar_dists = QtGui.QHBoxLayout()
        hbox_polar_dists.addWidget(self.gui_polar['dist_lbl'])
        hbox_polar_dists.addWidget(self.gui_polar['dist_txt'])
        layout_polar.addLayout(hbox_polar_dists)
        
        self.gui_polar['dists_btn'] = QtGui.QPushButton('Fit Distortions', frame_polar)
        self.gui_polar['dists_btn'].clicked.connect(self.on_fitDist)
        layout_polar.addWidget(self.gui_polar['dists_btn'])
        
        #self.gui_polar['upd_btn'] = QtGui.QPushButton('Update', frame_polar)
        #self.gui_polar['upd_btn'].clicked.connect(self.update_polar)
        #layout_polar.addWidget(self.gui_polar['upd_btn'])
        
        
        ## radial profile stuff
        frame_radprof = QtGui.QGroupBox('Radial Profile', self.mnwid)
        layout_radprof = QtGui.QVBoxLayout(frame_radprof)
        
        self.gui_radprof['rad_lbl'] = QtGui.QLabel('r_max, dr, sigma: (opt.)', frame_radprof)
        self.gui_radprof['rad_txt'] = QtGui.QLineEdit('', frame_radprof)
        hbox_radprof_rad = QtGui.QHBoxLayout()
        hbox_radprof_rad.addWidget(self.gui_radprof['rad_lbl'])
        hbox_radprof_rad.addWidget(self.gui_radprof['rad_txt'])
        layout_radprof.addLayout(hbox_radprof_rad)
        
        self.gui_radprof['crct_check'] = QtGui.QCheckBox('correct distortions', frame_radprof)
        self.gui_radprof['mask_btn'] = QtGui.QPushButton('Mask', frame_radprof)
        self.gui_radprof['mask_btn'].clicked.connect(self.on_mask)
        hbox_radprof_dist = QtGui.QHBoxLayout()
        hbox_radprof_dist.addWidget(self.gui_radprof['crct_check'])
        hbox_radprof_dist.addWidget(self.gui_radprof['mask_btn'])
        layout_radprof.addLayout(hbox_radprof_dist)
        
        self.gui_radprof['ext_btn'] = QtGui.QPushButton('Extract', frame_radprof)
        self.gui_radprof['ext_btn'].clicked.connect(self.on_extractRadProf)
        layout_radprof.addWidget(self.gui_radprof['ext_btn'])
        
        self.gui_radprof['fitxs_lbl'] = QtGui.QLabel('supp. x: ', frame_radprof)
        self.gui_radprof['fitxs_txt'] = QtGui.QLineEdit('', frame_radprof)
        self.gui_radprof['fitxs_btn'] = QtGui.QPushButton('Select', frame_radprof)
        self.gui_radprof['fitxs_btn'].setCheckable(True)
        self.gui_radprof['fitxs_btn'].clicked.connect(self.update_RadProf)
        hbox_radprof_fitxs = QtGui.QHBoxLayout()
        hbox_radprof_fitxs.addWidget(self.gui_radprof['fitxs_lbl'])
        hbox_radprof_fitxs.addWidget(self.gui_radprof['fitxs_txt'])
        hbox_radprof_fitxs.addWidget(self.gui_radprof['fitxs_btn'])
        layout_radprof.addLayout(hbox_radprof_fitxs)
        
        self.gui_radprof['fitxsw_lbl'] = QtGui.QLabel('width: ', frame_radprof)
        self.gui_radprof['fitxsw_txt'] = QtGui.QLineEdit('', frame_radprof)
        hbox_radprof_fitxsw = QtGui.QHBoxLayout()
        hbox_radprof_fitxsw.addWidget(self.gui_radprof['fitxsw_lbl'])
        hbox_radprof_fitxsw.addWidget(self.gui_radprof['fitxsw_txt'])
        layout_radprof.addLayout(hbox_radprof_fitxsw)
        
        self.gui_radprof['back_init_label'] = QtGui.QLabel('init offset, ampl., exp.: ')
        self.gui_radprof['back_init_txt'] = QtGui.QLineEdit('', frame_radprof)
        hbox_radprof_backinit = QtGui.QHBoxLayout()
        hbox_radprof_backinit.addWidget(self.gui_radprof['back_init_label'])
        hbox_radprof_backinit.addWidget(self.gui_radprof['back_init_txt'])
        layout_radprof.addLayout(hbox_radprof_backinit)
        
        self.gui_radprof['fitback_btn'] = QtGui.QPushButton('Fit Background', frame_radprof)
        self.gui_radprof['fitback_btn'].clicked.connect(self.on_subtractBackground)
        self.gui_radprof['back_check'] = QtGui.QCheckBox('Subtract Background', frame_radprof)
        self.gui_radprof['back_check'].stateChanged.connect(self.update_RadProf)
        #self.gui_radprof['back_btn'].setCheckable(True)
        #self.gui_radprof['back_btn'].clicked.connect(self.on_subtractBackground)
        hbox_radprof_backbtn = QtGui.QHBoxLayout()
        hbox_radprof_backbtn.addWidget(self.gui_radprof['fitback_btn'])
        hbox_radprof_backbtn.addWidget(self.gui_radprof['back_check'])
        layout_radprof.addLayout(hbox_radprof_backbtn)
        
        self.gui_radprof['fit_tbl'] = QtGui.QTableWidget(0,2, frame_radprof)
        self.gui_radprof['fit_tbl'].setHorizontalHeaderLabels(['function', 'initial parameters'])
        self.gui_radprof['fit_tbl'].horizontalHeader().setStretchLastSection(True)
        layout_radprof.addWidget(self.gui_radprof['fit_tbl'])
        
        self.gui_radprof['fit_add_btn'] = QtGui.QPushButton('Add', frame_radprof)
        self.gui_radprof['fit_add_btn'].clicked.connect(self.on_addFit)
        self.gui_radprof['fit_del_btn'] = QtGui.QPushButton('Delete', frame_radprof)
        self.gui_radprof['fit_del_btn'].clicked.connect(self.on_delFit)
        hbox_radprof_fitbtns = QtGui.QHBoxLayout()
        hbox_radprof_fitbtns.addWidget(self.gui_radprof['fit_add_btn'])
        hbox_radprof_fitbtns.addWidget(self.gui_radprof['fit_del_btn'])
        layout_radprof.addLayout(hbox_radprof_fitbtns)
        
        self.gui_radprof['fit_range_lbl'] = QtGui.QLabel('fit range: ', frame_radprof)
        self.gui_radprof['fit_range_txt'] = QtGui.QLineEdit('', frame_radprof)
        hbox_radprof_fitrange = QtGui.QHBoxLayout()
        hbox_radprof_fitrange.addWidget(self.gui_radprof['fit_range_lbl'])
        hbox_radprof_fitrange.addWidget(self.gui_radprof['fit_range_txt'])
        layout_radprof.addLayout(hbox_radprof_fitrange)
        
        self.gui_radprof['fit_btn'] = QtGui.QPushButton('Fit Radial Profile', frame_radprof)
        self.gui_radprof['fit_btn'].clicked.connect(self.on_fitRadProf)
        self.gui_radprof['fit_check'] = QtGui.QRadioButton('Plot Fit', frame_radprof)
        self.gui_radprof['fit_check'].toggled.connect(self.update_RadProf)
        self.gui_radprof['init_check'] = QtGui.QRadioButton('Plot Init', frame_radprof)
        self.gui_radprof['init_check'].toggled.connect(self.update_RadProf)
        self.gui_radprof['init_check'].toggle()
        hbox_radprof_fitbtns2 = QtGui.QHBoxLayout()
        hbox_radprof_fitbtns2.addWidget(self.gui_radprof['fit_btn'])
        vbox_radprof_tglbtns = QtGui.QVBoxLayout()
        vbox_radprof_tglbtns.addWidget(self.gui_radprof['fit_check'])
        vbox_radprof_tglbtns.addWidget(self.gui_radprof['init_check'])
        hbox_radprof_fitbtns2.addLayout(vbox_radprof_tglbtns)
        layout_radprof.addLayout(hbox_radprof_fitbtns2)
        
       
        ## Output stuff
        frame_run = QtGui.QGroupBox('Run', self.mnwid)
        layout_run = QtGui.QVBoxLayout(frame_run)
        
        self.gui_run['run_btn'] = QtGui.QPushButton('Run Current', frame_run)
        self.gui_run['run_btn'].clicked.connect(self.on_runsgl)
        self.gui_run['all_btn'] = QtGui.QPushButton('Run All', frame_run)
        self.gui_run['all_btn'].clicked.connect(self.on_runall)
        hbox_run = QtGui.QHBoxLayout()
        hbox_run.addWidget(self.gui_run['run_btn'])
        hbox_run.addWidget(self.gui_run['all_btn'])
        layout_run.addLayout(hbox_run)
        
        self.gui_run['out_btn'] = QtGui.QPushButton('Save evaluation to EMD', frame_run)
        self.gui_run['out_btn'].clicked.connect(self.on_saveEMDFile)
        layout_run.addWidget(self.gui_run['out_btn'])
        
        self.gui_run['savedist_btn'] = QtGui.QPushButton('Save corrected diffraction pattern', frame_run)
        self.gui_run['savedist_btn'].clicked.connect(self.on_saveCorrPattern)
        layout_run.addWidget(self.gui_run['savedist_btn'])
        
        
        vbox_left = QtGui.QVBoxLayout()
        vbox_left.addWidget(frame_files)
        
        sep1 = QtGui.QFrame(self.mnwid)
        sep1.setFrameStyle(QtGui.QFrame.HLine | QtGui.QFrame.Sunken)
        vbox_left.addWidget(sep1)
        
        vbox_left.addWidget(frame_localmax)
        
        sep2 = QtGui.QFrame(self.mnwid)
        sep2.setFrameStyle(QtGui.QFrame.HLine | QtGui.QFrame.Sunken)
        vbox_left.addWidget(sep2)
        
        vbox_left.addWidget(frame_polar)
        
        sep3 = QtGui.QFrame(self.mnwid)
        sep3.setFrameStyle(QtGui.QFrame.HLine | QtGui.QFrame.Sunken)
        vbox_left.addWidget(sep3)
        
        vbox_left.addWidget(frame_radprof)
        
        sep4 = QtGui.QFrame(self.mnwid)
        sep4.setFrameStyle(QtGui.QFrame.HLine | QtGui.QFrame.Sunken)
        vbox_left.addWidget(sep4)
        
        vbox_left.addWidget(frame_run)
        
        vbox_left.addStretch(1)
        
        
        # localmax tab
        self.plt_localmax = pg.PlotWidget()
        self.plt_localmax.setAspectLocked(True)
        self.plt_localmax.invertY(True)
        
        self.ch_localmax = {}
        self.ch_localmax['vLine'] = pg.InfiniteLine(angle=90, movable=False)
        self.ch_localmax['vLine'].setPen((0,255,0))
        self.ch_localmax['hLine'] = pg.InfiniteLine(angle=0, movable=False)
        self.ch_localmax['hLine'].setPen((0,255,0))
        
        proxy = pg.SignalProxy(self.plt_localmax.scene().sigMouseMoved, rateLimit=60, slot=self.localmax_mouseMoved)
        self.plt_localmax.scene().sigMouseClicked.connect(self.localmax_mouseClicked)
        
        # polar plot tab
        self.plt_polar = pg.PlotWidget()
        self.plt_polar.setMouseEnabled(x=False, y=True)
        self.plt_polar.setXRange(-np.pi, np.pi)
        axis1 = self.plt_polar.getAxis('bottom')
        axis1.setLabel('theta','rad')
        axis2 = self.plt_polar.getAxis('left')
        axis2.setLabel('r')
        
        # radial profile plot tab
        self.plt_radprof = pg.PlotWidget()
        
        self.ch_radprof = {}
        self.ch_radprof['vLine'] = pg.InfiniteLine(angle=90, movable=False)
        self.ch_radprof['vLine'].setPen((0,255,0))
        self.ch_radprof['hLine'] = pg.InfiniteLine(angle=0, movable=False)
        self.ch_radprof['hLine'].setPen((0,255,0))
        
        proxy2 = pg.SignalProxy(self.plt_radprof.scene().sigMouseMoved, rateLimit=60, slot=self.radprof_mouseMoved)
        self.plt_radprof.scene().sigMouseClicked.connect(self.radprof_mouseClicked)
        
        
        left = QtGui.QWidget(self.mnwid)
        hbox_left = QtGui.QHBoxLayout(left)
        hbox_left.addLayout(vbox_left)
        
        left_scroll = QtGui.QScrollArea()
        left_scroll.setWidget(left)
        
        #hbox.addWidget(self.imv_localmax)
        #hbox.addStretch(1)
        
        self.right = QtGui.QTabWidget(self.mnwid)
        self.right_tabs = {}
        self.right_tabs['localmax'] = self.right.addTab(self.plt_localmax, 'Local Maxima')
        self.right_tabs['polar'] = self.right.addTab(self.plt_polar, 'Polar Plot')
        self.right_tabs['radprof'] = self.right.addTab(self.plt_radprof, 'Radial Profile')
        #hbox_right = QtGui.QHBoxLayout(self.right)
        #hbox_right.addWidget(self.imv_localmax)
        
        splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        splitter.addWidget(left_scroll)
        splitter.addWidget(self.right)
        splitter.setStretchFactor(0, 0.2)
        splitter.setStretchFactor(1, 1.0)   
        
        
        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(splitter)
        
        self.mnwid.setLayout(hbox)
        
        self.setGeometry(100,100,1500,1100)
        self.setWindowTitle('Evaluation of Ring Diffraction Patterns')
        
        
        ## tooltips
        self.gui_file['re_btn'].setToolTip('Reopen an EMD evaluation file to conveniently load all the settings or resume previous work.\nNote that you have to rerun the evaluation, as no results are loaded, even if they are present in the EMD file.')
        self.gui_file['in_txt'].setToolTip('Use the Open button to load a file.')
        self.gui_file['in_btn'].setToolTip('Open an EMD file with diffraction patterns.\nOnly the first EMD group is read in. Dataset can be a single diffraction patterns or a series.')
        self.gui_file['idx_slider'].setToolTip('Control the index of the diffraction patterns to work on, if a series is loaded.\nEach diffraction patterns in the series has its own set of results, the settings are shared between all.')
        
        self.gui_localmax['txt_lmax_r'].setToolTip('Radius to use for local maxima determination.\nA local maxmimum will suppress all other local maxima within this radius.\nShould be a single integer.')
        self.gui_localmax['txt_lmax_thresh'].setToolTip('Threshold for local maxima determination.\nA local maximum is considered a local maximum, if its value at least this threshold above the minimum value found in the area described by the local radius.\nShould be a single value in intensity units.')
        self.gui_localmax['txt_lmax_cinit'].setToolTip('Initial guess for the center of the diffraction pattern.\nShould be two comma separated values in the displayed coordinate system.')
        self.gui_localmax['btn_lmax_cinit'].setToolTip('Select the initial guess for the center of the diffraction pattern from the image plot.\nMove the green crosshair to the center position and left-click.')
        self.gui_localmax['txt_lmax_range'].setToolTip('Filter the local maxima to only select those within a certain radial range to the initially guessed center.\nShould be two comma separated values in displayed coordinate system units.\nIf not provided, no filtering will take place.')
        self.gui_localmax['lmax_btn'].setToolTip('Run the local maxima detection.\nUpon completion, these are drawn as red circles in the image plot.\nThe goal is to find and select only points on a single ring, which will be used to determine the polar coordinate system of the diffraction patterns.')
        self.gui_localmax['min_slider'].setToolTip('Black level used to plot the image.')
        self.gui_localmax['max_slider'].setToolTip('White level used to plot the image.')
        
        self.gui_polar['cencpy_btn'].setToolTip('Copy initial center guess from above to use for polar plot.')
        self.gui_polar['cenopt_btn'].setToolTip('Optimize the current center to minimize the deviation from the mean radial distance of the found positions on one ring.')
        self.gui_polar['dist_txt'].setToolTip('Distortion orders to correct for.\nShould be a comma separated list of integers.')
        self.gui_polar['dists_btn'].setToolTip('Fit the distortions and plot the results to the polar plot.\nThe goal is to minimize the radial deviations.')
        
        self.gui_radprof['rad_txt'].setToolTip('Parameters for extracting the radial profile:\nr_max is the maximum radial distance considered\ndr is the stepsize on the radial axis\nsigma the width of the Gaussian kernel density estimatior applied\nShould be a comma separated list of three values. Can be left blank to use default values.')
        self.gui_radprof['crct_check'].setToolTip('Toogles whether the distortions fitted above are corrected during extraction of radial profile or not.')
        self.gui_radprof['mask_btn'].setToolTip('Launches a dialog to mask the input image. This allows to select areas, not to be used during extraction of the radial profile.')
        self.gui_radprof['ext_btn'].setToolTip('Extract the radial profile and plot it in the radial profile tab.')
        self.gui_radprof['fitxs_txt'].setToolTip('r values of supporting points in the radial profile, which will be used to fit a power law background.\nShould be a comma separated list of values.')
        self.gui_radprof['fitxs_btn'].setToolTip('Conveniently select the r values for supporting points in the radial profile plot.\nMove the green crosshair and left-click to select a point.\nPress the select button again to stop the selection process.')
        self.gui_radprof['fitxsw_txt'].setToolTip('Range in which points in the radial profile around the supporting points are considered for fitting the background.')
        self.gui_radprof['back_init_txt'].setToolTip('Initial parameters for fitting the power law background (offset + amplitude*r^exponent).\nShould be a list of three comma separated values.\nGood values usually are 1, 1000, -1.')
        self.gui_radprof['fitback_btn'].setToolTip('Fit power law background to the supporting points selected on the radial profile.')
        self.gui_radprof['back_check'].setToolTip('Check to subtract the fitted background.\nMake sure first, that the fitted background fits nicely to the radial profile.\nThis can be judged from the plot, when subtract background is not checked.')
        self.gui_radprof['fit_tbl'].setToolTip('List of peak functions fitted to the radial profile.\nFirst column is the name of the fit function, second column contains comma separated lists of values used as initial guess for the parameters.\nUse the option to plot the initial guess to get the initial parameters right.\nFor now just use \'voigt\' functions, which require lists of four comma separated values as initial parameters.')
        self.gui_radprof['fit_add_btn'].setToolTip('Add a row to the list of fit functions.')
        self.gui_radprof['fit_del_btn'].setToolTip('Delete the selected row(s) of the list of fit functions.')
        self.gui_radprof['fit_range_txt'].setToolTip('Range of the radial profile used during fit.\nGenerally used to cut of lower and upper end of the radial profile.\nShould be a list of two comma separated values.')
        self.gui_radprof['fit_btn'].setToolTip('Fit the radial profile with the given peak functions.')
        self.gui_radprof['fit_check'].setToolTip('Show the fitted functions in the radial profile plot.')
        self.gui_radprof['init_check'].setToolTip('Show the functions with initial parameters in the radial profile plot.')
        
        self.gui_run['run_btn'].setToolTip('Run the evaluation on the currently selected image.\nEquivalent to subsequently hitting the single buttons in order.')
        self.gui_run['all_btn'].setToolTip('Run the evaluation on all images in the series.\nEquivalent to hitting \'Run Current\' for each index in the series.')
        self.gui_run['out_btn'].setToolTip('Save the settings and the results in an EMD file.')
        self.gui_run['savedist_btn'].setToolTip('Save the diffraction pattern corrected for the fitted distortions in an EMD file.')
        
        
        self.plt_localmax.setToolTip('Image plot used to find local maxima.\nThe diffraction pattern is plotted according to the intensity sliders on the left.\nThe detected and optionally filtered local maxima are plotted as red circles.')
        self.plt_polar.setToolTip('Plot of the selected points on a ring in polar coordinates.\nThis plot is used to optimize the center and fit the distortions.\nRed crosses are points before correction, green after correction.')
        self.plt_radprof.setToolTip('Plot of the radial profile (red line).\nIf subtract background is not checked, the fitted power law is plotted as blue line with black crosses marking the supporting points used for the fit.\nOtherwise the fit functions are either plotted using the fitted parameters or the initial guess.')
        
        self.show()
        


    def log(self, msg):
        self.log_txt.appendPlainText(msg)

        
        
    def keyPressEvent(self, e):
        
        if e.key() == QtCore.Qt.Key_Escape:
            self.close()
            


    def closeEvent(self, event):
        self.log('Goodbye!')
        del self.femd_in
        self.log_wdg.close()
        QtGui.QMainWindow.closeEvent(self, event)
        
         
            
    def on_open_infile(self):
        '''
        Open an EMD file with input images.
        '''
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open EMD file', filter='EMD files (*.emd);;All files (*.*)')
        try:
            self.femd_in = ncempy.io.emd.fileEMD(fname, readonly=True)
            
            self.gui_file['in_txt'].setText(fname)
            
            self.log('Opening file "{}" ..'.format(fname))
            
            # get and save data
            data, dims = self.femd_in.get_emdgroup(self.femd_in.list_emds[0])
            self.log('.. found EMD type group at "{}".'.format(self.femd_in.list_emds[0].name))
            
            if len(data.shape)==2:
                data = data[np.newaxis,:,:]
                self.log('.. single image.')
            elif len(data.shape)==3:
                self.log('.. series of {:d} images.'.format(data.shape[2])) 
            
            self.data = np.copy(data)
            self.dims = copy.deepcopy(dims)
            
            # configure index slider
            self.idx = 0
            
            if len(data.shape)==3:
                self.gui_file['idx_slider'].setMinimum(0)
                self.gui_file['idx_slider'].setMaximum(self.data.shape[0]-1)
            
            # reset parameters
            self.settings = {}
            self.points = [None]*data.shape[0]
            self.center = [None]*data.shape[0]
            self.dists = [None]*data.shape[0]
            self.radprof = [None]*data.shape[0]
            self.back_params = [None]*data.shape[0]
            self.res = [None]*data.shape[0]
            self.mask = None
            
            # configure intensity sliders
            min_data = np.min(self.data)
            max_data = np.max(self.data)
            
            self.gui_localmax['min_slider'].setMinimum(min_data)
            self.gui_localmax['min_slider'].setMaximum(max_data)
            self.gui_localmax['min_slider'].setValue(min_data)
            
            self.gui_localmax['max_slider'].setMinimum(min_data)
            self.gui_localmax['max_slider'].setMaximum(max_data)
            self.gui_localmax['max_slider'].setValue(max_data)
            
            self.update_localmax()
            
        except: 
            self.gui_file['in_txt'].setText('')
            self.log('Failed to open EMD file "{}".'.format(fname))
            raise TypeError()
    
    
    def on_reopen(self):
    
        fname = QtGui.QFileDialog.getOpenFileName( self, 'Open EMD evaluation file', filter='EMD files (*.emd);;All files (*.*)')
        
        try:
            femd = ncempy.io.emd.fileEMD(fname, readonly=True)
            
            self.log('Reopening evaluation file "{}" ..'.format(fname))
            
            # search todo group
            todo = []
        
            # recursive function
            def proc_group(grp,todo):
                for item in grp:
                    if grp.get(item, getclass=True) == h5py._hl.group.Group:
                        item = grp.get(item)
                        if 'type' in item.attrs:
                            if item.attrs['type'] == np.string_(ncempy.eval.ring_diff.cur_eva_vers):
                                todo.append(item)
                        proc_group(item, todo)
            
            proc_group(femd.file_hdl, todo)
            
            grp_eva = todo[0]
            self.log('.. found evaluation group at "{}".'.format(grp_eva.name))
            
            self.femd_in = ncempy.io.emd.fileEMD(grp_eva.attrs['filename'].decode('utf-8'), readonly=True)
            self.gui_file['in_txt'].setText(grp_eva.attrs['filename'].decode('utf-8'))
            
            # get and save data
            data, dims = self.femd_in.get_emdgroup(self.femd_in.list_emds[0])
            
            self.log('.. loading data from "{}"-"{}".'.format(self.femd_in.file_hdl.filename, self.femd_in.list_emds[0].name))
            
            if len(data.shape)==2:
                data = data[np.newaxis,:,:]
                self.log('.. single image.')
            elif len(data.shape)==3:
                self.log('.. series of {:d} images.'.format(data.shape[2])) 
            
            self.data = np.copy(data)
            self.dims = copy.deepcopy(dims)
                
            # configure index slider
            self.idx = 0
                
            if len(data.shape)==3:
                self.gui_file['idx_slider'].setMinimum(0)
                self.gui_file['idx_slider'].setMaximum(self.data.shape[0]-1)
            
            # reset parameters
            self.settings = ncempy.eval.ring_diff.get_settings(grp_eva['settings_ringdiffraction'])
            
            self.gui_localmax['txt_lmax_r'].setText('{}'.format(self.settings['lmax_r']))
            self.gui_localmax['txt_lmax_thresh'].setText('{}'.format(self.settings['lmax_thresh']))
            self.gui_localmax['txt_lmax_cinit'].setText('{:g}, {:g}'.format(self.settings['lmax_cinit'][1], self.settings['lmax_cinit'][0]))
            self.gui_localmax['txt_lmax_range'].setText('{:g}, {:g}'.format(self.settings['lmax_range'][0], self.settings['lmax_range'][1]))
            
            self.gui_polar['dist_txt'].setText( ', '.join(map(str, self.settings['ns']) ) )
            
            self.gui_radprof['rad_txt'].setText( '{:g}, {:g}, {:g}'.format(self.settings['rad_rmax'], self.settings['rad_dr'], self.settings['rad_sigma']) )
            self.gui_radprof['fitxs_txt'].setText( ', '.join(map('{:g}'.format, self.settings['back_xs'])) )
            self.gui_radprof['fitxsw_txt'].setText( '{:g}'.format(self.settings['back_xswidth']) )
            self.gui_radprof['back_init_txt'].setText( ', '.join(map('{:g}'.format, self.settings['back_init'])) )
               
            n = 0
            for i in range(len(self.settings['fit_funcs'])):
                self.gui_radprof['fit_tbl'].insertRow(self.gui_radprof['fit_tbl'].rowCount())
                new_func = QtGui.QTableWidgetItem( self.settings['fit_funcs'][i] )
                new_init = QtGui.QTableWidgetItem( ', '.join(map('{:g}'.format, self.settings['fit_init'][n:n+ncempy.algo.math.lkp_funcs[self.settings['fit_funcs'][i]][1]])) )
                n += ncempy.algo.math.lkp_funcs[self.settings['fit_funcs'][i]][1]
                self.gui_radprof['fit_tbl'].setItem(i,0,new_func)
                self.gui_radprof['fit_tbl'].setItem(i,1,new_init)
            
            self.gui_radprof['fit_range_txt'].setText( ', '.join(map('{:g}'.format, self.settings['fit_rrange'])) )
            
            self.points = [None]*data.shape[0]
            self.center = [None]*data.shape[0]
            self.dists = [None]*data.shape[0]
            self.radprof = [None]*data.shape[0]
            self.back_params = [None]*data.shape[0]
            self.res = [None]*data.shape[0]
            
            # get mask from settings but move it to self.
            if self.settings['mask'] is None:
                self.mask = None
            else:
                self.mask = self.settings['mask']
            del self.settings['mask']
                
            self.log('.. reusing saved parameter settings.')
                
                
            # configure intensity sliders
            min_data = np.min(self.data)
            max_data = np.max(self.data)
            
            new_min = int(min_data + self.settings['plt_imgminmax'][0]*(max_data - min_data) )
            new_max = int(min_data + self.settings['plt_imgminmax'][1]*(max_data - min_data) )
               
            self.gui_localmax['min_slider'].setMinimum(min_data)
            self.gui_localmax['min_slider'].setMaximum(max_data)
            self.gui_localmax['min_slider'].setValue( new_min )
               
            self.gui_localmax['max_slider'].setMinimum(min_data)
            self.gui_localmax['max_slider'].setMaximum(max_data)
            self.gui_localmax['max_slider'].setValue( new_max )
            
            self.update_localmax()
      
            del femd
            
        except:
            self.log('Reopening EMD evaluation file "{}" failed.'.format(fname))
            raise RuntimeError()
            
    
    def on_saveEMDFile(self):
        '''
        Open or create an EMD file to hold the evaluation.
        '''
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Save to EMD file', filter='EMD files (*.emd);;All files (*.*)')
        
        try:
            self.log('Saving evaluation to file "{}".'.format(fname))
            
            # create/overwrite outfile
            if os.path.isfile(fname):
                os.remove(fname)
            femd = ncempy.io.emd.fileEMD(fname)            
                
            if 'evaluation' in femd.file_hdl:
                grp_eva = femd.file_hdl['evaluation']
            else:
                grp_eva = femd.file_hdl.create_group('evaluation')
            
            hdl = ncempy.eval.ring_diff.put_sglgroup(grp_eva, self.femd_in.file_hdl.filename.split('/')[-1], self.femd_in.list_emds[0])
            
            # insert the result datasets        
            
            n = None
            if len(self.dims) == 3:
                n = None
                for i in range(len(self.radprof)):
                    if not self.radprof[i] is None:
                        n = i
                        radprof = np.zeros( (self.radprof[i].shape[0], len(self.radprof)) )
                        for ii in range(len(self.radprof)):
                            if not self.radprof[ii] is None:
                                radprof[:,ii] = self.radprof[ii][:,1]
                        femd.put_emdgroup('radial_profile', radprof, ( (self.radprof[n][:,0], 'radial distance', self.dims[2][2]), self.dims[0] ), parent=hdl, overwrite=True)
                        break
                
                for i in range(len(self.res)):
                    if not self.res[i] is None:
                        res = np.zeros( (self.res[i].shape[0], len(self.res)) )
                        for ii in range(len(self.res)):
                            if not self.res[ii] is None:
                                res[:,ii] = self.res[ii][:]
                        femd.put_emdgroup('fit_results', res, ( ( np.array(range(res.shape[0])), 'parameters', '[]'), self.dims[0] ), parent=hdl, overwrite=True)
                        break
                
                for i in range(len(self.center)):
                    if not self.center[i] is None:
                        centers = np.zeros( (2, len(self.center)) )
                        for ii in range(len(self.center)):
                            if not self.center[ii] is None:
                                centers[:,ii] = self.center[ii][:]
                        femd.put_emdgroup('centers', centers, ( ( np.array(range(2)), 'dimension', self.dims[2][2]), self.dims[0] ), parent=hdl, overwrite=True) 
                        break
                
                for i in range(len(self.dists)):
                    if not self.dists[i] is None:
                        dists = np.zeros( (self.dists[i].shape[0], len(self.dists)) )
                        for ii in range(len(self.dists)):
                            if not self.dists[ii] is None:
                                dists[:,ii] = self.dists[ii][:]
                        femd.put_emdgroup('distortions', dists, ( ( np.array(range(dists.shape[0])), 'parameters', '[]'), self.dims[0] ), parent=hdl, overwrite=True)
                        break
                
                for i in range(len(self.back_params)):
                    if not self.back_params[i] is None:
                        back_params = np.zeros( (self.back_params[i].shape[0], len(self.back_params)) )
                        for ii in range(len(self.back_params)):
                            if not self.back_params[ii] is None:
                                back_params[:,ii] = self.back_params[ii][:]
                        femd.put_emdgroup('back_results', back_params, ( ( np.array(range(back_params.shape[0])), 'background parameters', '[]'), self.dims[0] ), parent=hdl, overwrite=True)
                        break
                
            else:
                if not self.radprof[self.idx] is None:
                    femd.put_emdgroup('radial_profile', self.radprof[0][:,1], ( (self.radprof[0][:,0], 'radial_distance', self.dims[0][2]), ), parent=hdl, overwrite=True)
                if not self.res[self.idx] is None:
                    femd.put_emdgroup('fit_results', self.res[0], ( ( np.array(range(self.res[0].shape[0])), 'parameters', '[]'), ), parent=hdl, overwrite=True)
                if not self.center[self.idx] is None:
                    femd.put_emdgroup('centers', self.center[0], ( ( np.array(range(2)), 'dimension', self.dims[0][2]), ), parent=hdl, overwrite=True) 
                if not self.dists[self.idx] is None:
                    femd.put_emdgroup('distortions', self.dists[0], ( ( np.array(range(self.dists[0].shape[0])), 'parameters', '[]'), ), parent=hdl, overwrite=True)
                if not self.back_params[self.idx] is None:
                    femd.put_emdgroup('back_results', self.back_params[0], ( ( np.array(range(self.back_params[0].shape[0])), 'background parameters', '[]'), ), parent=hdl, overwrite=True)

            # save settings 
            # make them valid
            # create a copy of settings to decouple
            mysettings = copy.deepcopy(self.settings)
            # reinsert mask
            mysettings['mask'] = self.mask
            for key in ncempy.eval.ring_diff.min_dummie_settings:
                if not key in mysettings:
                    mysettings[key] = ncempy.eval.ring_diff.min_dummie_settings[key]
            
            ncempy.eval.ring_diff.put_settings( hdl, mysettings )
            
            # save a comment
            femd.put_comment( 'Ring diffraction evaluation saved from gui_ringdiff.' )

            # close file
            del femd
            
        except:
            self.log('Saving evaluation to file "{}" failed.'.format(fname))
            raise RuntimeError()
            
    
    def on_intensitySlider(self):
        '''
        Intensity sliders to change the img plot in local maxima.
        '''
        
        # get values
        min_val = self.gui_localmax['min_slider'].value()
        max_val = self.gui_localmax['max_slider'].value()
        
        # update labels
        self.gui_localmax['min_value'].setText('{:d}'.format(min_val))
        self.gui_localmax['max_value'].setText('{:d}'.format(max_val))     
         
        # update plot
        if not self.plt_localmax_img is None:
            self.plt_localmax_img.setLevels( (min_val, max_val) )
        
        # update settings
        min_data = np.min(self.data)
        max_data = np.max(self.data)
        self.settings['plt_imgminmax'] = ( (min_val-min_data)/(max_data - min_data) , (max_val-min_data)/(max_data - min_data) )



    def on_indexSlider(self):
        '''
        Index slider to run through image series.
        '''
        
        # get value
        val = self.gui_file['idx_slider'].value()
        
        # update label
        self.gui_file['idx_value'].setText('{:d}'.format(val))
        
        # save in main
        self.idx = val
        
        # update plot
        temp = self.right.currentWidget()
        self.update_RadProf()
        self.update_polar()
        self.update_localmax()
        self.right.setCurrentWidget(temp)


    def update_localmax(self):
        '''
        Update the localmax plot view.
        '''
        
        # update plot in local maxima
        self.plt_localmax.clear()
        
        if not self.data is None:
        
            # sgl image versus series
            if len(self.dims) == 2:
            
                # set axis
                axis1 = self.plt_localmax.getAxis('bottom')
                axis1.setLabel(self.dims[1][1], self.dims[1][2])
                axis2 = self.plt_localmax.getAxis('left')
                axis2.setLabel(self.dims[0][1], self.dims[0][2])
               
                # plot image
                self.plt_localmax_img = pg.ImageItem(self.data[self.idx,:,:].transpose().astype('float64'), levels=(self.gui_localmax['min_slider'].value(), self.gui_localmax['max_slider'].value()))
                self.plt_localmax_img.setZValue(-100)
                self.plt_localmax_img.setRect(pg.QtCore.QRectF( self.dims[1][0][0],self.dims[0][0][0],self.dims[1][0][-1]-self.dims[1][0][0],self.dims[0][0][-1]-self.dims[0][0][0]))
                self.plt_localmax.addItem(self.plt_localmax_img)
                
            else:
            
                # set axis
                axis1 = self.plt_localmax.getAxis('bottom')
                axis1.setLabel(self.dims[2][1], self.dims[2][2])
                axis2 = self.plt_localmax.getAxis('left')
                axis2.setLabel(self.dims[1][1], self.dims[1][2])
               
                # plot image
                self.plt_localmax_img = pg.ImageItem(self.data[self.idx,:,:].transpose().astype('float64'), levels=(self.gui_localmax['min_slider'].value(), self.gui_localmax['max_slider'].value()))
                self.plt_localmax_img.setZValue(-100)
                self.plt_localmax_img.setRect(pg.QtCore.QRectF( self.dims[2][0][0],self.dims[1][0][0],self.dims[2][0][-1]-self.dims[2][0][0],self.dims[1][0][-1]-self.dims[1][0][0]))
                self.plt_localmax.addItem(self.plt_localmax_img)
                
            # plot points
            if not self.points[self.idx] is None:
            
                # draw points
                self.plt_localmax.plot(self.points[self.idx][:,1], self.points[self.idx][:,0], pen=None, symbol='o', symbolPen=(255,0,0), symbolBrush=None)

        if self.gui_localmax['btn_lmax_cinit'].isChecked():
            self.plt_localmax.getPlotItem().addItem(self.ch_localmax['vLine'], ignoreBounds=True)
            self.plt_localmax.getPlotItem().addItem(self.ch_localmax['hLine'], ignoreBounds=True)


        self.right.setCurrentWidget( self.right.widget(self.right_tabs['localmax']))

        

    def on_localmax(self):
        '''
        Calculate local maxima.
        '''
    
        self.log('Calculating local maxima ..')
    
        # parse the input
        try:
            max_r = float(self.gui_localmax['txt_lmax_r'].text())
            thresh = float(self.gui_localmax['txt_lmax_thresh'].text())
                
            cinit = self.gui_localmax['txt_lmax_cinit'].text().strip()
            rrange = self.gui_localmax['txt_lmax_range'].text().strip()
                
            if cinit == '' and rrange == '':
                cinit = []
                rrange = []
            else:
                cinit = [float(item.strip()) for item in cinit.split(',')]
                cinit = cinit[::-1]
                rrange = [float(item.strip()) for item in rrange.split(',')]
                assert(len(cinit)==2 and len(rrange)==2)
                
        except:
            self.log('.. could not parse the given input.')
            raise TypeError('Bad input to local maxima')

        # update settings
        self.settings['lmax_r'] = max_r
        self.settings['lmax_thresh'] = thresh
        self.settings['lmax_cinit'] = cinit
        self.settings['lmax_range'] = rrange
        
        if len(self.dims) == 2:
            c_dims = self.dims
        else:
            c_dims = self.dims[1:3]
        
        
        # find local max
        points = ncempy.algo.local_max.local_max(self.data[self.idx,:,:], self.settings['lmax_r'], self.settings['lmax_thresh'])
        points = ncempy.algo.local_max.points_todim(points, c_dims)
                
        self.log('.. found {:d} candidate points.'.format(points.shape[0]))
                
        # filter to single ring if input provided
        if (len(self.settings['lmax_cinit'])==2 and len(self.settings['lmax_range'])==2):
            points = ncempy.algo.distortion.filter_ring(points, self.settings['lmax_cinit'], self.settings['lmax_range'])
            self.log('.. filtered to {:d} points.'.format(points.shape[0]))
        
        # save points in main
        self.points[self.idx] = points
        
        # update localmax view
        self.update_localmax() 
            

    
    def localmax_mouseMoved(self, evt):
        pos = evt[0]
        
        mousePoint = self.plt_localmax.getPlotItem().vb.mapSceneToView(pos)
        
        self.ch_localmax['vLine'].setPos(mousePoint.x())
        self.ch_localmax['hLine'].setPos(mousePoint.y())
        
            
    def localmax_mouseClicked(self, evt):
    
        if self.gui_localmax['btn_lmax_cinit'].isChecked():
    
            mousePoint = self.plt_localmax.getPlotItem().vb.mapSceneToView(evt.scenePos())
        
            cinit = (mousePoint.y(), mousePoint.x())
            self.settings['lmax_cinit'] = cinit
            self.log('Selected initial center guess: ({:g}, {:g}).'.format(self.settings['lmax_cinit'][0], self.settings['lmax_cinit'][1]))
            self.gui_localmax['txt_lmax_cinit'].setText( ', '.join(map('{:g}'.format, self.settings['lmax_cinit'][::-1])) )

            self.gui_localmax['btn_lmax_cinit'].setChecked(False)
            
            self.update_localmax()
            
    
    
    def update_polar(self):
        '''
        Update the polar plot depending on which data is present.
        '''
        
        self.plt_polar.clear()
        
        if not self.center[self.idx] is None:
            # update center in left column
            self.gui_polar['center_lbl'].setText('center: ({:.3f}, {:.3f})'.format(self.center[self.idx][0], self.center[self.idx][1]))
        
            if not self.points[self.idx] is None:
                points_plr = ncempy.algo.distortion.points_topolar(self.points[self.idx], self.center[self.idx])
        
                # horizontal mean line
                self.plt_polar.plot( [-np.pi, np.pi], [np.mean(points_plr[:,0]), np.mean(points_plr[:,0])], pen=pg.mkPen('k', style=QtCore.Qt.DashLine))
                
                # uncorrected radial positions
                self.plt_polar.plot(points_plr[:,1], points_plr[:,0], pen=None, symbol='x', symbolPen=(255,0,0), symbolBrush=None)
        
                if not self.dists[self.idx] is None:
                    
                    xpl_ell = np.linspace(-np.pi, np.pi, 200)
                    # single distortions
                    for i in range(len(self.settings['ns'])):
                        self.plt_polar.plot( xpl_ell, self.dists[self.idx][0]*ncempy.algo.distortion.rad_dis(xpl_ell, self.dists[self.idx][i*2+1], self.dists[self.idx][i*2+2], self.settings['ns'][i]), pen=pg.mkPen('m', style=QtCore.Qt.DashLine) )
                    
                    # sum of distortions
                    sum_dists = self.dists[self.idx][0]*np.ones(xpl_ell.shape)
                    for i in range(len(self.settings['ns'])):
                        sum_dists *= ncempy.algo.distortion.rad_dis( xpl_ell, self.dists[self.idx][i*2+1], self.dists[self.idx][i*2+2], self.settings['ns'][i])
                    self.plt_polar.plot( xpl_ell, sum_dists, pen=(0,0,255))

                    # corrected radial positions
                    points_plr_corr = np.copy(points_plr)
                    for i in range(len(self.settings['ns'])):
                        points_plr_corr[:,0] /= ncempy.algo.distortion.rad_dis(points_plr_corr[:,1], self.dists[self.idx][i*2+1], self.dists[self.idx][i*2+2], self.settings['ns'][i])
                    
                    self.plt_polar.plot( points_plr_corr[:,1], points_plr_corr[:,0], pen=None, symbol='x', symbolPen=(0,180,0), symbolBrush=None) 

        self.right.setCurrentWidget( self.right.widget(self.right_tabs['polar']))



    def on_copyCenter(self):
        
        try:
            cinit = self.gui_localmax['txt_lmax_cinit'].text().strip()   
            if cinit == '':
                cinit = []
            else:
                cinit = [float(item.strip()) for item in cinit.split(',')]
                assert(len(cinit)==2)
        except:
            self.log('Failed to parse initial center guess from input.')
            raise TypeError('Bad input to local maxima')
            
        if len(cinit) == 0:
            self.log('Failed to copy initial center guess from input.')
        else:
            self.settings['lmax_cinit'] = cinit
            self.center[self.idx]=np.array(cinit)
            self.log('Copied initial center guess: ({:g}, {:g}).'.format(self.center[self.idx][0], self.center[self.idx][1]))

        self.update_polar()
        


    def on_optimizeCenter(self):
    
        self.log('Optimizing center ..')
        center = ncempy.algo.distortion.optimize_center(self.points[self.idx], self.center[self.idx])
        self.center[self.idx]= center
        self.log('.. optimized to ({:g}, {:g}).'.format(self.center[self.idx][0], self.center[self.idx][1]))
        
        self.update_polar()
        
        
        
    def on_fitDist(self):
        
        try:
            ns = self.gui_polar['dist_txt'].text().strip()
            if ns == '':
                ns = []
            else:
                ns = [int(item.strip()) for item in ns.split(',')]
        except:
            self.log('Failed to parse input to fit distortions.')
            raise TypeError('Bad input in ')
            
        # update setting
        self.settings['ns'] = ns 
       
        self.log('Fitting distortions of order(s) {}'.format(self.settings['ns']))
       
        # run optimization
        if len(ns) >= 1:
            points_plr = ncempy.algo.distortion.points_topolar(self.points[self.idx], self.center[self.idx])
            dists = ncempy.algo.distortion.optimize_distortion(points_plr, self.settings['ns'])
            
            self.dists[self.idx] = dists
            
            self.log('.. fitted to:')
            self.log('.. .. radius: {:g}'.format(self.dists[self.idx][0]))
            
            for i in range(len(self.settings['ns'])):
                self.log('.. .. order: {:d}, alpha: {:g}, beta: {:g}'.format(self.settings['ns'][i], self.dists[self.idx][1+i*2], self.dists[self.idx][2+i*2]))
            
        self.update_polar()
        
        
    
    def update_RadProf(self):
        
        self.plt_radprof.clear()
        
        if len(self.dims) == 2:
            c_dims = self.dims
        else:
            c_dims = self.dims[1:3]
        
        if (not self.radprof is None) and (not self.radprof[self.idx] is None):
        
            axis1 = self.plt_radprof.getAxis('bottom')
            axis1.setLabel('radial distance', c_dims[0][2])
            axis2 = self.plt_radprof.getAxis('left')
            axis2.setLabel('intensity')
        
            # rad profile has been calculated
            R = np.copy(self.radprof[self.idx][:,0])
            I = np.copy(self.radprof[self.idx][:,1])
            
            # cut to fitrange
            if (not len(self.settings['fit_rrange']) == 0) and (self.gui_radprof['back_check'].isChecked()):
                sel = (R>=self.settings['fit_rrange'][0])*(R<=self.settings['fit_rrange'][1])
                I = I[sel]
                R = R[sel]
        
            back = None
        
            if not self.back_params[self.idx] is None:
            
                if 'back_xs' in self.settings:
                    # calulate background
                    fit_R = np.array([])
                    fit_I = np.array([])
                    for xpoint in self.settings['back_xs']:
                        ix = np.where( np.abs(R-xpoint) <= self.settings['back_xswidth'])
                        fit_R = np.append(fit_R, R[ix])
                        fit_I = np.append(fit_I, I[ix])
                       
                    back = ncempy.algo.math.sum_functions( R, ('const', 'powlaw'), self.back_params[self.idx] )
     
                    if not self.gui_radprof['back_check'].isChecked():
                        self.plt_radprof.plot(fit_R, fit_I, pen=None, symbol='x', symbolPen=(0,0,0))
                        self.plt_radprof.plot(R, back, pen=(0,0,255))
     
     
            if self.gui_radprof['back_check'].isChecked() and not back is None:
                # subtract background
                I -= back
                
            
            if (not self.res[self.idx] is None) and (self.gui_radprof['fit_check'].isChecked()):
            
                # draw fit results
                fitsum = ncempy.algo.math.sum_functions( R, self.settings['fit_funcs'], self.res[self.idx] )
                self.plt_radprof.plot(R, fitsum, pen=(0,0,255))
                
                i = 0
                for n in range(len(self.settings['fit_funcs'])):
                    self.plt_radprof.plot(R, ncempy.algo.math.lkp_funcs[self.settings['fit_funcs'][n]][0](R, self.res[self.idx][i:i+ncempy.algo.math.lkp_funcs[self.settings['fit_funcs'][n]][1]]), pen=(0,180,0)) 
                    i += ncempy.algo.math.lkp_funcs[self.settings['fit_funcs'][n]][1]
                    
            elif ('fit_init' in self.settings) and (self.gui_radprof['init_check'].isChecked()):
            
                # draw fit results
                initsum = ncempy.algo.math.sum_functions( R, self.settings['fit_funcs'], self.settings['fit_init'] )
                self.plt_radprof.plot(R, initsum, pen=(0,0,255))
                
                i = 0
                for n in range(len(self.settings['fit_funcs'])):
                    self.plt_radprof.plot(R, ncempy.algo.math.lkp_funcs[self.settings['fit_funcs'][n]][0](R, self.settings['fit_init'][i:i+ncempy.algo.math.lkp_funcs[self.settings['fit_funcs'][n]][1]]), pen=(0,180,0)) 
                    i += ncempy.algo.math.lkp_funcs[self.settings['fit_funcs'][n]][1]
                    
    
            # plot radial profile
            self.plt_radprof.plot(R, I, pen=(255,0,0))
        
        if self.gui_radprof['fitxs_btn'].isChecked():
            self.plt_radprof.getPlotItem().addItem(self.ch_radprof['vLine'], ignoreBounds=True)
            self.plt_radprof.getPlotItem().addItem(self.ch_radprof['hLine'], ignoreBounds=True)
        
        self.right.setCurrentWidget( self.right.widget(self.right_tabs['radprof']))
        
        
    
    def on_extractRadProf(self):
        '''
        Extract the radial profile from the diffraction pattern.
        '''
        
        self.log('Extracting radial profile ..')
        
        try:
            pars = self.gui_radprof['rad_txt'].text().strip()
            if pars == '':
                pars = []
            else:
                pars = [float(item.strip()) for item in pars.split(',')]
                assert(len(pars)==3)
                
            fitrange = self.gui_radprof['fit_range_txt'].text().strip()
            if fitrange == '':
                fitrange = []
            else:
                fitrange = [float(item.strip()) for item in fitrange.split(',')]
                assert(len(fitrange)==2)
                   
        except:
            self.log('Failed to parse given input to extract radial profile.')
            raise TypeError('Bad input in ')
        
        if len(self.dims) == 2:
            c_dims = self.dims
        else:
            c_dims = self.dims[1:3]
        
        # save settings
        if len(pars) == 0:
            self.settings['rad_rmax'] = np.abs(c_dims[1][0][0]-c_dims[1][0][1])*np.min(self.data[self.idx,:,:].shape)/2.0
            self.settings['rad_dr'] = np.abs(c_dims[1][0][0]-c_dims[1][0][1])/10.
            self.settings['rad_sigma'] = np.abs(c_dims[1][0][0]-c_dims[1][0][1])
            self.log('.. calculating adaptive defaults: r_max: {:g}, dr: {:g}, sigma: {:g}'.format(self.settings['rad_rmax'], self.settings['rad_dr'], self.settings['rad_sigma']))
        else:
            self.settings['rad_rmax'] = pars[0]
            self.settings['rad_dr'] = pars[1]
            self.settings['rad_sigma'] = pars[2]
            self.log('.. using given r_max, dr and sigma.')
        
        self.settings['fit_rrange'] = fitrange
        
        # get the polar coordinate system
        if self.gui_radprof['crct_check'].isChecked():
            rs, thes = ncempy.algo.radial_profile.calc_polarcoords( self.center[self.idx], c_dims, self.settings['ns'], self.dists[self.idx] )
            self.log('.. calculating coordinate system, correcting for distortions.')
        else:
            rs, thes = ncempy.algo.radial_profile.calc_polarcoords( self.center[self.idx], c_dims )
            self.log('.. calculating coordinate system, not correcting distortions.')
        
        # get the radial profile
        R, I = ncempy.algo.radial_profile.calc_radialprofile( self.data[self.idx,:,:], rs, self.settings['rad_rmax'], self.settings['rad_dr'], self.settings['rad_sigma'], self.mask )
        
        # save in main
        self.radprof[self.idx] = np.array([R,I]).transpose()
        self.log('.. extracted radial profile.')
    
        # update the plot
        self.update_RadProf()



    def on_mask(self):
        
        if not self.data is None:
            mask, flag = MaskDialog.getMask(self)
            
            if flag:
                self.mask = mask
                self.log('Created new mask for extracting radial profile.')
            else:
                self.mask = None
                self.log('Removed any mask for extracting radial profile.')


        
    def on_subtractBackground(self):
        '''
        Fit power law background and subtract.
        '''
        
        self.log('Fitting background to radial profile.')
        
        # parse the input
        try:
            fitxs_width = float(self.gui_radprof['fitxsw_txt'].text())
                
            fitxs = self.gui_radprof['fitxs_txt'].text().strip()
            back_init = self.gui_radprof['back_init_txt'].text().strip()    
                
            if fitxs == '':
                fitxs = []
            else:
                fitxs = [float(item.strip()) for item in fitxs.split(',')]
                assert(len(fitxs) >= 1)
            
            if back_init == '':
                back_init = [1,1,1]
            else:
                back_init = [float(item.strip()) for item in back_init.split(',')]
                assert(len(back_init) == 3)
        except:
            self.log('Failed to parse input given to fit background.')
            raise TypeError('Bad input to subtract background.')
            
        # save settings
        self.settings['back_xs'] = fitxs
        self.settings['back_xswidth'] = fitxs_width
        self.settings['back_init'] = back_init
       
        # get background points
        fit_R = np.array([])
        fit_I = np.array([])
        for xpoint in self.settings['back_xs']:
            ix = np.where( np.abs(self.radprof[self.idx][:,0]-xpoint) <= self.settings['back_xswidth'])
            fit_R = np.append(fit_R, self.radprof[self.idx][ix,0])
            fit_I = np.append(fit_I, self.radprof[self.idx][ix,1])
        
        self.log('.. found {:d} points to support background.'.format(fit_I.shape[0]))
        
        # fit power law
        funcs_back = [ 'const', 'powlaw' ]
        res_back = ncempy.algo.radial_profile.fit_radialprofile( fit_R, fit_I, funcs_back, self.settings['back_init'], maxfev=1000 )
        
        # save output
        self.back_params[self.idx] = res_back
        
        self.log('.. background fitted to: offset: {:g}, amplitude: {:g}, exponent: {:g}'.format(self.back_params[self.idx][0], self.back_params[self.idx][1], self.back_params[self.idx][2]))
        
        self.update_RadProf()
       
        
   
    def on_fitRadProf(self):
        '''
        Fit the radial profile with peak functions.
        '''
        
        self.log('Fitting radial profile ..')
        
        # parse the input
        try:
            fitrange = self.gui_radprof['fit_range_txt'].text().strip()
            
            if fitrange == '':
                fitrange = []
            else:
                fitrange = [float(item.strip()) for item in fitrange.split(',')]
                assert(len(fitrange)==2)
                
            funcs = []
            init_guess = []    
               
            for n in range(self.gui_radprof['fit_tbl'].rowCount()):
                label = self.gui_radprof['fit_tbl'].item(n,0).text()
                funcs.append(label)

                initparams = self.gui_radprof['fit_tbl'].item(n,1).text()
                initparams = [float(item.strip()) for item in initparams.split(',')]
                assert(len(initparams)==ncempy.algo.math.lkp_funcs[label][1])  
                for param in initparams:
                    init_guess.append(param)        
        
        except:
            self.log('Failed to parse input given to fit radial profile.')
            raise TypeError('Bad input to fit radial profile.')
        
        
        # save settings
        self.settings['fit_rrange'] = fitrange
        self.settings['fit_funcs'] = funcs
        self.settings['fit_init'] = init_guess
        
        
        if not self.radprof is None:
        
            # rad profile has been calculated
            R = np.copy(self.radprof[self.idx][:,0])
            I = np.copy(self.radprof[self.idx][:,1])
            
            # cut to fitrange
            if not len(self.settings['fit_rrange']) == 0:
                sel = (R>=self.settings['fit_rrange'][0])*(R<=self.settings['fit_rrange'][1])
                I = I[sel]
                R = R[sel]
                self.log('.. fitting in the range of {:g} to {:g}.'.format(self.settings['fit_rrange'][0], self.settings['fit_rrange'][1]))
        
            back = None
        
            if (not self.back_params[self.idx] is None) and (self.gui_radprof['back_check'].isChecked()):
            
                # calulate background
                fit_R = np.array([])
                fit_I = np.array([])
                for xpoint in self.settings['back_xs']:
                    ix = np.where( np.abs(R-xpoint) <= self.settings['back_xswidth'])
                    fit_R = np.append(fit_R, R[ix])
                    fit_I = np.append(fit_I, I[ix])
                       
                back = ncempy.algo.math.sum_functions( R, ('const', 'powlaw'), self.back_params[self.idx] )
         
                I -= back
                
                self.log('.. subtracting background.')
 
 
            res = ncempy.algo.radial_profile.fit_radialprofile( R, I, self.settings['fit_funcs'], self.settings['fit_init'], 10000 )
            
            self.res[self.idx] = res
            self.log('.. peaks fitted to:')
            n = 0
            for i in range(len(self.settings['fit_funcs'])):
                self.log('.. .. {}: {}'.format(self.settings['fit_funcs'][i],  ', '.join(map('{:g}'.format, self.res[self.idx][n:n+ncempy.algo.math.lkp_funcs[self.settings['fit_funcs'][i]][1]]))))
                n += ncempy.algo.math.lkp_funcs[self.settings['fit_funcs'][i]][1]
                        
        self.update_RadProf()
        


    def on_addFit(self):
        self.gui_radprof['fit_tbl'].insertRow(self.gui_radprof['fit_tbl'].rowCount())
    
    
    def on_delFit(self):
    
        row_to_delete = []
    
        for item in self.gui_radprof['fit_tbl'].selectedIndexes():
            this_row = item.row()
            if not this_row in row_to_delete:
                row_to_delete.append(this_row)
        
        for row in sorted(row_to_delete)[::-1]:
            self.gui_radprof['fit_tbl'].removeRow(row)
    
    
    def radprof_mouseMoved(self, evt):
        pos = evt[0]
        
        mousePoint = self.plt_radprof.getPlotItem().vb.mapSceneToView(pos)
        
        self.ch_radprof['vLine'].setPos(mousePoint.x())
        self.ch_radprof['hLine'].setPos(mousePoint.y())
        
            
    def radprof_mouseClicked(self, evt):
    
        if self.gui_radprof['fitxs_btn'].isChecked():
    
            mousePoint = self.plt_radprof.getPlotItem().vb.mapSceneToView(evt.scenePos())
            
            fitxs = self.gui_radprof['fitxs_txt'].text().strip()
            
            if fitxs == '':
                fitxs = []
            else:
                fitxs = [float(item.strip()) for item in fitxs.split(',')]
            
            fitxs.append(mousePoint.x())
            
            self.settings['back_xs'] = fitxs
            
            self.gui_radprof['fitxs_txt'].setText( ', '.join(map('{:g}'.format, self.settings['back_xs'])) )
            
            self.update_RadProf()
            

    
    def on_runsgl(self):
    
        self.on_localmax()
        
        self.on_copyCenter()
        
        self.on_optimizeCenter()
        
        self.on_fitDist()
        
        self.on_extractRadProf()
        
        self.on_subtractBackground()
        
        self.on_fitRadProf()
        
    
    def on_runall(self):
        
        if len(self.dims) == 3:
            for i in range(self.data.shape[0]):
            
                self.idx = i
            
                self.log('Running evaluation of {:d}/{:d}.'.format(self.idx+1, self.data.shape[0]))
                
                self.on_runsgl()
        else:
            self.on_runsgl()
    
    
    def on_saveCorrPattern(self):
    
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Save to EMD file', filter='EMD files (*.emd);;All files (*.*)')
        
        try:
            self.log('Saving corrected diffraction pattern to file "{}".'.format(fname))
            
            # create/overwrite outfile
            if os.path.isfile(fname):
                os.remove(fname)
            femd = ncempy.io.emd.fileEMD(fname)            
                
            grp = femd.file_hdl['data']
            
            data_corr = np.zeros(self.data.shape)
            
            if len(self.dims) == 2:
                c_dims = self.dims
            else:
                c_dims = self.dims[1:3]
            
            for i in range(self.data.shape[0]):
                
                self.idx = i
                
                if (not self.center[self.idx] is None) and (not self.dists[self.idx] is None):
                 
                    self.log('.. correcting image {:d}/{:d}.'.format(self.idx+1, self.data.shape[0]))
                        
                    data_corr[self.idx,:,:] = ncempy.algo.radial_profile.correct_distortion( self.data[self.idx,:,:], c_dims, self.center[self.idx], self.settings['ns'], self.dists[self.idx])
                else:
                    self.log('.. skipping image {:d}/{:d} for missing results.'.format(self.idx+1, self.data.shape[0]))
            
            if len(self.dims) == 3:
                femd.put_emdgroup(self.femd_in.file_hdl.filename.split('/')[-1], data_corr, self.dims, parent=grp)
            else:
                femd.put_emdgroup(self.femd_in.file_hdl.filename.split('/')[-1], data_corr[self.idx,:,:], self.dims, parent=grp)
            
            del femd
        
        except:
            self.log('Saving corrected diffraction pattern failed.')
            raise
    
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    main = Main()
    sys.exit(app.exec_())
