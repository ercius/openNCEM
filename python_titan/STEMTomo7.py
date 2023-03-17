""" Python script for acquiring STEM tomography tilt series
manually. This script saves all metdata and images
in an EMD Berkeley file. It allows the user to take
multiple images per tilt angle or a STEM scan
rotation series for dose fractionation and drift
correction. It works with both the TEAM Stage
and the FEI Compustage which can be selected
using input arguments when the script is initially run.

author: Peter Ercius, percius@lbl.gov, Wolfgang Theis, Birmingham UK
"""

import uuid

import numpy as np
import scipy as sci
import socket
from sys import path as path2
import os
import wx  # wxPython GUI package
import h5py
import threading
import time
import math
from copy import copy

# Import only necessary modules from matplotlib
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.figure import Figure

# For connections to microscope and TIA
from comtypes.client import CreateObject
from comtypes.safearray import safearray_as_ndarray  # get data across COM barrier fast

# Connect to TEAM Stage
try:
    import TEAMstageclass
except ModuleNotFoundError:
    print('TEAM stage module not found.')

import argparse
version = 7.2  # version number for this program


class TEAMFrame(wx.Frame):
    def __init__(self, parent, title, stagetype='compustage'):

        try:
            from matplotlib.pyplot import set_cmap
        except:
            pass
        # Setup the stagetype as compustage or teamstage
        self.stagetype = stagetype

        # Set the update stage position to either
        # compustage or teamstage function
        if stagetype == 'teamstage':
            self.getPosition = self.getPosition_teamstage  # point to this function
            self.microscopeName = 'TEAM 0.5'
            self.stageType = 'NCEM TEAM Stage'
        elif stagetype == 'compustage':
            self.getPosition = self.getPosition_compustage  # point to this function
            self.microscopeName = 'FEI'
            self.stageType = 'compustage'

        self.TIA = None
        self._microscope = None
        self.Stage = None
        self.Acq = None
        self.Proj = None
        self. Ill = None
        self.detector0 = None

        # For review tab
        self.showindex = 0
        self.showsubindex = 0
        self.showset = 0
        self.curtiltnr = 0

        # Initialize the base Frame
        wx.Frame.__init__(self, parent, title=title, size=(1300, 1000))

        # menu 
        #menubar = wx.MenuBar()
        #filem = wx.Menu()
        #editm = wx.Menu()
        #helpm = wx.Menu()
        #menubar.Append(filem, '&File')
        #menubar.Append(editm, '&Edit')
        #menubar.Append(helpm, '&Help')
        #self.SetMenuBar(menubar)

        # define sizers     
        mainbox = wx.BoxSizer(wx.HORIZONTAL)
        imagecolumnbox = wx.BoxSizer(wx.VERTICAL)
        infobox = wx.BoxSizer(wx.VERTICAL)
        imagebox = wx.BoxSizer(wx.HORIZONTAL)
        zoomfftbox = wx.BoxSizer(wx.VERTICAL)
        displayvalbox = wx.BoxSizer(wx.VERTICAL)
        imagebuttonbox = wx.BoxSizer(wx.HORIZONTAL)
        currentvalbox = wx.BoxSizer(wx.VERTICAL)
        aquirebox = wx.BoxSizer(wx.VERTICAL)
        tiltbox = wx.BoxSizer(wx.VERTICAL)

        notebook = wx.Notebook(self)

        # create the page windows as children of the notebook
        searchpage = wx.Panel(notebook)
        tomosetuppage = wx.Panel(notebook)
        aquirepage = wx.Panel(notebook)
        reviewpage = wx.Panel(notebook)

        # Figures
        # Create a figure in a canvas
        self.TEMfigure = Figure()
        self.canvas = FigureCanvasWxAgg(self, -1, self.TEMfigure)
        self.axes = self.TEMfigure.add_subplot(111)
        self.axes.xaxis.set_visible(False)
        self.axes.yaxis.set_visible(False)
        try:
            set_cmap('gray')
        except:
            pass

        # Create zoom figure in a canvas
        # self.TEMfigureZoom = Figure()
        # self.Zoomcanvas = FigureCanvasWxAgg(self, -1, self.TEMfigureZoom)
        # self.Zoomaxes = self.TEMfigureZoom.add_subplot(111)
        # self.Zoomaxes.xaxis.set_visible(False)
        # self.Zoomaxes.yaxis.set_visible(False)

        # Create fft figure in a canvas
        # self.TEMfigurefft = Figure()
        # self.fftcanvas = FigureCanvasWxAgg(self, -1, self.TEMfigurefft)
        # self.fftaxes = self.TEMfigurefft.add_subplot(111)
        # self.fftaxes.xaxis.set_visible(False)
        # self.fftaxes.yaxis.set_visible(False)
        
        # set_cmap does not work in older matplotlib versions
        try:
            set_cmap('gray')
        except:
            pass

        # Image Acquisition page
        self._lLabel1 = wx.StaticText(aquirepage, wx.ID_ANY, 'File Parameters')
        self._lFprefix = wx.StaticText(aquirepage, wx.ID_ANY, 'File name')
        self._iFprefix = wx.TextCtrl(aquirepage, wx.ID_ANY, value='TiltSeries1')
        self._lDir = wx.StaticText(aquirepage, wx.ID_ANY, 'Output directory')
        try:
            self._iDir = wx.TextCtrl(aquirepage, wx.ID_ANY, value='G:\UserData',
                                     size=(300, 20))  # correct for TEAM 0.5
        except:
            self._iDir = wx.TextCtrl(aquirepage, wx.ID_ANY, value=os.getcwd(), size=(300, 20))
        self._btnNewFile = wx.Button(aquirepage, label='Update File Name')

        self._lLabel3 = wx.StaticText(aquirepage, wx.ID_ANY, 'Acquire parameters')
        self._lDwell = wx.StaticText(aquirepage, wx.ID_ANY, 'Dwell time (usec)')
        self._iDwell = wx.TextCtrl(aquirepage, wx.ID_ANY, value='12')
        self._lBin = wx.StaticText(aquirepage, wx.ID_ANY, 'Binning')
        self._iBin = wx.TextCtrl(aquirepage, wx.ID_ANY, value='4')
        self._lRep = wx.StaticText(aquirepage, wx.ID_ANY, 'Number of images per aquire')
        self._iRep = wx.TextCtrl(aquirepage, wx.ID_ANY, value='1')
        self._lRot = wx.StaticText(aquirepage, wx.ID_ANY, 'Num scan rotations (2 or 4 override repeats)')
        self._iRot = wx.RadioBox(aquirepage, wx.ID_ANY, label='Rotation Series', choices=('0', '2', '4'))
        self._iRot.SetSelection(0)
        self._lDel = wx.StaticText(aquirepage, wx.ID_ANY, 'Delay between images (sec)')
        self._iDel = wx.TextCtrl(aquirepage, wx.ID_ANY, value='0')
        self._btnAcq = wx.Button(aquirepage, label='Acquire')
        aquirebox.Add(self._lLabel1)
        aquirebox.Add(self._lFprefix)
        aquirebox.Add(self._iFprefix)
        aquirebox.Add(self._btnNewFile)
        aquirebox.Add(self._lDir)
        aquirebox.Add(self._iDir)
        aquirebox.Add(self._lLabel3)
        aquirebox.Add(self._lDwell)
        aquirebox.Add(self._iDwell)
        aquirebox.Add(self._lBin)
        aquirebox.Add(self._iBin)
        aquirebox.Add(self._lRep)
        aquirebox.Add(self._iRep)
        aquirebox.Add(self._lRot)
        aquirebox.Add(self._iRot)
        aquirebox.Add(self._lDel)
        aquirebox.Add(self._iDel)
        aquirebox.Add(self._btnAcq)
        
        # Add Go to button to simplify moving the stage.
        # Use a dictionary instead of many variables to try a new way
        self._gotoGUI = {}
        self._gotoGUI['_lLabel31'] = wx.StaticText(aquirepage, wx.ID_ANY, 'Set stage angles (TS only)')
        self._gotoGUI['_lGotoA'] = wx.StaticText(aquirepage, wx.ID_ANY, 'Alpha (deg)')
        self._gotoGUI['_lGotoG'] = wx.StaticText(aquirepage, wx.ID_ANY, 'Gamma (deg)')
        self._gotoGUI['_iGotoA'] = wx.TextCtrl(aquirepage, wx.ID_ANY, value='0')
        self._gotoGUI['_iGotoG'] = wx.TextCtrl(aquirepage, wx.ID_ANY, value='0')
        self._gotoGUI['_btnGotoAngles'] = wx.Button(aquirepage, label='Go to angles')
        self._gotoGUI['_cGotoConfirm'] = wx.CheckBox(aquirepage, label='Confirm go to', pos=(20, 20))
        self._gotoGUI['_cGotoConfirm'].SetValue(True)
        for ii in ('_lLabel31', '_lGotoA','_iGotoA', '_lGotoG', '_iGotoG','_cGotoConfirm','_btnGotoAngles'):
            aquirebox.Add(self._gotoGUI[ii])

        # Review page
        self._lDispPos = wx.StaticText(reviewpage, wx.ID_ANY, 'Displ pos')
        self._lDispAlphaGamma = wx.StaticText(reviewpage, wx.ID_ANY, 'Displ alpha gamma')
        self._lDispFinePos = wx.StaticText(reviewpage, wx.ID_ANY, 'Displ fine pos')
        self._btnDispGotoPos = wx.Button(reviewpage, label='Goto display position')
        self._lDispIndex = wx.StaticText(reviewpage, wx.ID_ANY, 'Index')
        displayvalbox.Add(self._lDispPos)
        displayvalbox.Add(self._lDispAlphaGamma)
        displayvalbox.Add(self._lDispFinePos)
        displayvalbox.Add(self._lDispIndex)
        displayvalbox.Add(self._btnDispGotoPos, proportion=0, flag=wx.ALIGN_RIGHT | wx.RIGHT, border=10)
        self._btnprev = wx.Button(reviewpage, label='Show prev')
        self._btnnext = wx.Button(reviewpage, label='Show next')
        imagebuttonbox.Add(self._btnprev, proportion=0, flag=wx.ALIGN_RIGHT | wx.RIGHT, border=10)
        imagebuttonbox.Add(self._btnnext, proportion=0, flag=wx.ALIGN_RIGHT | wx.RIGHT, border=10)
        displayvalbox.Add(imagebuttonbox, proportion=0, flag=wx.LEFT | wx.RIGHT | wx.EXPAND, border=20)

        # search page
        self._lPos = wx.StaticText(searchpage, wx.ID_ANY, 'positions')
        self._lAlphaGamma = wx.StaticText(searchpage, wx.ID_ANY, 'alpha gamma')
        self._lFinePos = wx.StaticText(searchpage, wx.ID_ANY, 'fine positions')
        # self._lEucentric = wx.StaticText(searchpage, wx.ID_ANY, 'eucentric')
        # self._cbfft = wx.CheckBox(searchpage, label='display fft', pos=(20, 20))
        # self._cbfft.SetValue(False)
        self._btnGetPos = wx.Button(searchpage, label='Get Position')
        self._btnSearch = wx.Button(searchpage, label='Search')
        currentvalbox.Add(self._lPos, proportion=0, flag=wx.LEFT | wx.RIGHT | wx.EXPAND, border=10)
        currentvalbox.Add(self._lAlphaGamma, proportion=0, flag=wx.LEFT | wx.RIGHT | wx.EXPAND, border=10)
        currentvalbox.Add(self._lFinePos, proportion=0, flag=wx.LEFT | wx.RIGHT | wx.EXPAND, border=10)
        currentvalbox.Add(self._btnGetPos, proportion=0, flag=wx.ALIGN_RIGHT | wx.RIGHT, border=10)
        currentvalbox.Add(self._btnSearch, proportion=0, flag=wx.ALIGN_RIGHT | wx.RIGHT, border=10)
        # currentvalbox.Add(self._cbfft, proportion=0, flag=wx.ALIGN_RIGHT | wx.RIGHT, border=10)

        # Set label fonts X
        # note: wx.Font(pointSize, family, style, weight, underline=False, faceName="", encoding=wx.FONTENCODING_DEFAULT)
        self.setFonts(True)
        # self.fontTopLabels = wx.Font(14, wx.DEFAULT, wx.NORMAL,wx.BOLD)
        # self.fontLabels = wx.Font(12, wx.DEFAULT, wx.NORMAL,wx.BOLD)
        # self.fontDisable = wx.Font(12,wx.DEFAULT,wx.ITALIC,wx.NORMAL) #not sure why this does not work, but its not used yet
        # self._lLabel1.SetFont(fontTopLabels)
        # self._lLabel3.SetFont(fontTopLabels)
        # self._lFprefix.SetFont(fontLabels)
        # self._lDir.SetFont(fontLabels)
        # self._lDwell.SetFont(fontLabels)
        # self._lBin.SetFont(fontLabels)
        # self._lRep.SetFont(fontLabels)
        # self._lRot.SetFont(fontLabels)
        # self._lDel.SetFont(fontLabels)
        # self._lPos.SetFont(fontLabels)
        # self._lAlphaGamma.SetFont(fontLabels)
        # self._lFinePos.SetFont(fontLabels)

        # tomo setup page
        self._lslope = wx.StaticText(tomosetuppage, wx.ID_ANY, 'slope')
        self._islope = wx.TextCtrl(tomosetuppage, wx.ID_ANY, value='16')
        self._lA1 = wx.StaticText(tomosetuppage, wx.ID_ANY, 'alpha 1')
        self._lA2 = wx.StaticText(tomosetuppage, wx.ID_ANY, 'alpha 2')
        self._lG1 = wx.StaticText(tomosetuppage, wx.ID_ANY, 'gamma 1')
        self._lG2 = wx.StaticText(tomosetuppage, wx.ID_ANY, 'gamma 2')
        self._iA1 = wx.TextCtrl(tomosetuppage, wx.ID_ANY, value='63.4')
        self._iA2 = wx.TextCtrl(tomosetuppage, wx.ID_ANY, value='91.1')
        self._iG1 = wx.TextCtrl(tomosetuppage, wx.ID_ANY, value='117.6')
        self._iG2 = wx.TextCtrl(tomosetuppage, wx.ID_ANY, value='-98.8')
        self._btnCalcTilt = wx.Button(tomosetuppage, label='Calc tilt sequence')
        self._btnShowNextTilt = wx.Button(tomosetuppage, label='Show next')
        self._btnShowPrevTilt = wx.Button(tomosetuppage, label='Show prev')
        self._lselectedtilt = wx.StaticText(tomosetuppage, wx.ID_ANY, 'tilt data')
        self._btnCalcTilt.Bind(wx.EVT_BUTTON, self.onCalcTilt)  # bind the button event
        self._btnShowPrevTilt.Bind(wx.EVT_BUTTON, self.onPrevTilt)  # bind the button event
        self._btnShowNextTilt.Bind(wx.EVT_BUTTON, self.onNextTilt)  # bind the button event
        self._btnTilMove = wx.Button(tomosetuppage, label='Move to orientation - not implemented')
        self._cb = wx.CheckBox(tomosetuppage, label='gamma tilt', pos=(20, 20))
        self._cb.SetValue(True)
        tiltgrid = wx.GridSizer(4, 2, 0, 0)
        tiltgrid.AddMany(
            [self._lslope, self._islope, self._lA1, self._iA1, self._lG1, self._iG1, self._lA2, self._iA2, self._lG2,
             self._iG2])
        tiltbox.Add(self._cb)
        tiltbox.Add(tiltgrid)
        tiltbox.Add(self._btnCalcTilt)
        tiltbox.Add(self._btnShowNextTilt)
        tiltbox.Add(self._btnShowPrevTilt)
        tiltbox.Add(self._lselectedtilt)
        tiltbox.Add(self._btnTilMove)

        # Set color
        self.SetBackgroundColour("0000FF")

        # set resonable axis for start (will be 17 micro meters at 5000x, don't want to start at 1m)
        # for now remember none loaded yet
        self.imagevalid = 0

        # Bind the events
        self._btnAcq.Bind(wx.EVT_BUTTON, self.onAcquire)  # bind the button event
        self._btnSearch.Bind(wx.EVT_BUTTON, self.onSearch)  # bind the button event
        self._btnGetPos.Bind(wx.EVT_BUTTON, self.onGetPos)  # bind the button event
        ## Disabled mouse events
        # self.canvas.mpl_connect('button_press_event', self.onClick) #bind the left click event
        self._btnNewFile.Bind(wx.EVT_BUTTON, self.onNewFileName)  # bind the button event
        self._btnprev.Bind(wx.EVT_BUTTON, self.onShowPrev)  # bind the button event
        self._btnnext.Bind(wx.EVT_BUTTON, self.onShowNext)  # bind the button event
        self._btnDispGotoPos.Bind(wx.EVT_BUTTON, self.onGotoDispPos)
        self._gotoGUI['_btnGotoAngles'].Bind(wx.EVT_BUTTON, self.onGotoAngles)
        #zoomfftbox.Add(self.Zoomcanvas, proportion=1, flag=wx.LEFT | wx.RIGHT | wx.EXPAND, border=0)
        #zoomfftbox.Add(self.fftcanvas, proportion=1, flag=wx.LEFT | wx.RIGHT | wx.EXPAND, border=0)

        imagebox.Add(self.canvas, proportion=3, flag=wx.LEFT | wx.RIGHT | wx.EXPAND, border=0)
        #imagebox.Add(zoomfftbox, proportion=1, flag=wx.LEFT | wx.RIGHT | wx.EXPAND, border=0)

        # imagecolumnbox.Add(displayvalbox, proportion=0, flag=wx.LEFT|wx.RIGHT|wx.EXPAND, border=20)
        imagecolumnbox.Add(imagebox, proportion=1, flag=wx.LEFT | wx.RIGHT | wx.EXPAND, border=0)

        infobox.Add(currentvalbox, proportion=0, flag=wx.LEFT | wx.RIGHT | wx.EXPAND, border=0)
        # infobox.Add(aquirebox, proportion=0, flag=wx.LEFT|wx.RIGHT|wx.EXPAND, border=0)

        searchpage.SetSizer(infobox);
        tomosetuppage.SetSizer(tiltbox);
        aquirepage.SetSizer(aquirebox);
        reviewpage.SetSizer(displayvalbox);

        # add the pages to the notebook with the label to show on the tab
        notebook.AddPage(searchpage, "search")
        notebook.AddPage(tomosetuppage, "setup")
        notebook.AddPage(aquirepage, "aquire")
        notebook.AddPage(reviewpage, "review")

        mainbox.Add(imagecolumnbox, proportion=1, flag=wx.LEFT | wx.RIGHT | wx.EXPAND, border=0)
        mainbox.Add(notebook, proportion=0, flag=wx.LEFT | wx.RIGHT | wx.EXPAND, border=20)

        self.SetSizer(mainbox)  # Initialize the parent sizer

        self.sb = self.CreateStatusBar()
        self.sb.SetStatusText('Idle...')

        self.onConnect()

        self.Show(True)  # show the frame
        # self.NewFileName() #crashes - still busy with some init process?

    def __del__(self):
        try:
            self.TS.TS_Disconnect()
        except:
            pass
        self.f.close()
        self.TS.TS_Disconnect()

    def drawPlot(self, imArray):
        self.axes.clear()
        self.axes.imshow(imArray)  # add the image to the figure
        self.canvas.draw()  # draw the panel
        self.Zoomcanvas.draw()  # draw the panel
        self.fftcanvas.draw()  # draw the panel

    # Connect to the microscope and setup STEM detectors
    def onConnect(self):
        # connect using comtypes
        try:
            self._microscope = CreateObject('TEMScripting.Instrument')
            print("Connected to microscope")
        except:
            print("Microscope connection failed")
            raise
        try:
            self.TIA = CreateObject('ESVision.Application')
            print('Connected to TIA')
        except:
            print('Connection to TIA failed')
            raise

        # Connect to TEAM stage
        print('Connecting to TEAM Stage...')
        self.TS = TEAMstageclass.TEAMstage()
        try:
            self.TS.TS_Connect('localhost', 5557)
        except Exception as e:
            print('Error with TEAM Stage connection. Exception type is {}'.format(e))
        print('Success')
        
        try:
            # Get microscope interfaces
            self.Acq = self._microscope.Acquisition
            self.Proj = self._microscope.Projection
            self.Ill = self._microscope.Illumination
            self.Stage = self._microscope.Stage  # FEI stage. Not needed for TEAM Stage

            self.detector0 = self.Acq.Detectors(0)  # older pythoncom versions might require square brackets []
        except:
            print('Connections to microscope interfaces failed.')
            raise

        self._microscope.Acquisition.AddAcqDevice(self.detector0);
        
        # Determine the maximum binning value
        self.determineMaxBinning()
        
        # Initialize user acquisition values
        self.getUsrAcqVals()

    def setDatasetAttributesEditable(self, ed):
        self._iBin.SetEditable(ed)
        self._iRep.SetEditable(ed)
        self._iRot.EnableItem(0, ed)  # disable/enable radio buttons
        self._iRot.EnableItem(1, ed)
        self._iRot.EnableItem(2, ed)
        self._iDwell.SetEditable(ed)
        self._iDel.SetEditable(ed)

        self.setFonts(ed)

    def setFonts(self, ed):
        if ed:
            fontTopLabels = wx.Font(14, wx.DEFAULT, wx.NORMAL, wx.BOLD)
            fontLabels = wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.BOLD)
        else:
            fontTopLabels = wx.Font(14, wx.DEFAULT, wx.ITALIC, wx.BOLD)
            fontLabels = wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD)

        self._lLabel1.SetFont(fontTopLabels)
        self._lLabel3.SetFont(fontTopLabels)
        self._lFprefix.SetFont(fontLabels)
        self._lDir.SetFont(fontLabels)
        self._lDwell.SetFont(fontLabels)
        self._lBin.SetFont(fontLabels)
        self._lRep.SetFont(fontLabels)
        self._lRot.SetFont(fontLabels)
        self._lDel.SetFont(fontLabels)
        self._lPos.SetFont(fontLabels)
        self._lAlphaGamma.SetFont(fontLabels)
        self._lFinePos.SetFont(fontLabels)
        
        for ii in self._gotoGUI.values():
            ii.SetFont(fontLabels)
        
        self._gotoGUI['_lLabel31'].SetFont(fontTopLabels)
        
    def initializeEMDFile(self):
        """ Initialize the file and tags for HDF5 / EMD file format

        """
        print('init emd file: ' + self.fullName)

        with h5py.File(self.fullName + '.emd', 'w') as f:
            f['/'].attrs['version'] = version

            dataroot = f.create_group('data')
            dataroot.attrs['filename'] = str(self.Fprefix)
            dataroot.attrs['stemtomo version'] = version

            dataTop = dataroot.create_group('raw')
            dataTop.attrs['nextindex'] = 0
            dataTop.attrs['nextset'] = 0
            dataTop.attrs['emd_group_type'] = 1
            dataTop.attrs['uuid'] = str(uuid.uuid1())
            
            dataTop.attrs['rotSetting'] = self.rotSetting
            if self.rotSetting == 1:
                dataTop.attrs['rotations'] = [0, 90]
                self.Rep = 2
            elif self.rotSetting == 2:
                dataTop.attrs['rotations'] = [0, 90, 45, 135]
                self.Rep = 4
            else:
                dataTop.attrs['rotations'] = 'none'

            dataTop.attrs['images per set'] = self.Rep
            dataTop.attrs['binning'] = self.Bin
            dataTop.attrs['max binning'] = self.maxBin
            dataTop.attrs['dwell time'] = self.Dwell
            dataTop.attrs['delay time'] = self.Del

            pix = self.maxBin / self.Bin
            # Create dataset with shape [alpha, repeats, Y, X]
            _ = dataTop.create_dataset('data', (1, self.Rep, pix, pix,), '<u2', maxshape=(None, self.Rep, pix, pix),
                                       chunks=(1, 1, pix, pix))

            # Create the EMD dimension datasets
            dim4 = dataTop.create_dataset('dim4', (1,), 'f', maxshape=(pix,))
            dim4.attrs['name'] = np.string_('X')
            dim4.attrs['units'] = np.string_('n_m')
            dim3 = dataTop.create_dataset('dim3', (1,), 'f', maxshape=(pix,))
            dim3.attrs['name'] = np.string_('Y')
            dim3.attrs['units'] = np.string_('n_m')
            dim1 = dataTop.create_dataset('dim1', (1,), 'f', maxshape=(None,))
            dim1.attrs['name'] = np.string_('Alpha')
            dim1.attrs['units'] = np.string_('degrees')

            dim2 = dataTop.create_dataset('dim2', (1,), 'f', maxshape=(self.Rep,))
            if self.rotSetting == 1 | self.rotSetting == 2:
                dim2.attrs['name'] = np.string_('Scan rotation')  # Need to make this STEM rotation angle if applicable
                dim2.attrs['units'] = np.string_('degrees')
            else:
                dim2.attrs['name'] = np.string_('image in set')  # Need to make this STEM rotation angle if applicable
                dim2.attrs['units'] = np.string_('num')

            # Create other datasets to hold other metadata during acquisition
            tiltangles = dataTop.create_dataset('tiltangles', (1,), 'f', maxshape=(None,))
            tiltangles.attrs['values'] = ['tilt angle']
            tiltangles.attrs['index'] = ['alpha tilt (dim1)']

            imagesetparameters = dataTop.create_dataset('imageSetParameters', (6, 1), 'f', maxshape=(6, None))
            imagesetparameters.attrs['values'] = ['x', 'y', 'z', 'alpha', 'gamma', 'timestamp']
            imagesetparameters.attrs['index'] = ['tilt number (dim1)']

            imageparameters = dataTop.create_dataset('imageParameters', (6, 1), 'f', maxshape=(6, None))
            imageparameters.attrs['values'] = ['image set', 'rep number', 'scanRotation', 'pixelSizeX', 'pixelSizeY',
                                               'timestamp']
            imageparameters.attrs['index'] = ['image number']

            microscope = f.create_group('microscope')
            microscope.attrs['name'] = self.microscopeName
            microscope.attrs['high tension'] = self._microscope.Gun.HTValue
            microscope.attrs['spot size'] = self.Ill.SpotsizeIndex
            microscope.attrs['magnification'] = self.Ill.StemMagnification
            microscope.attrs['convergence angle'] = self.Ill.ConvergenceAngle
            microscope.attrs['camera length'] = self.Proj.CameraLength
            microscope.attrs['stage type'] = self.stageType
            
            stage = f.create_group('user')

    def setStemAcqVals(self):
        self.myStemAcqParams = self.Acq.Detectors.AcqParams
        self.myStemAcqParams.Binning = self.Bin
        self.myStemAcqParams.ImageSize = 0  # self._microscope.ACQIMAGESIZE_FULL
        self.myStemAcqParams.DwellTime = self.Dwell
        self.Acq.Detectors.AcqParams = self.myStemAcqParams

    def setStemSearchVals(self):
        self.myStemSearchParams = self.Acq.Detectors.AcqParams
        self.myStemSearchParams.Binning = 8
        self.myStemSearchParams.ImageSize = 0  # self._microscope.ACQIMAGESIZE_FULL
        self.myStemSearchParams.DwellTime = 6e-6
        self.Acq.Detectors.AcqParams = self.myStemSearchParams

    def getUsrAcqVals(self):
        self.Fprefix = self._iFprefix.GetValue()
        self.Fdir = self._iDir.GetValue()
        self.fullName = self.Fdir + os.sep + self.Fprefix
        # self.Angle = int(self._iAngle.GetValue())
        self.Dwell = int(self._iDwell.GetValue()) * 1e-6  # Change to microseconds
        self.Bin = int(self._iBin.GetValue())
        self.Rep = int(self._iRep.GetValue())
        self.rotSetting = self._iRot.GetSelection()
        self.Del = int(self._iDel.GetValue())

    def onClick(self, event):
        # Determine mouse position and move stage. xdata = x coord of mouse in data coords
        if (self.imagevalid == 0) or (event.xdata == None) or (event.ydata == None):
            pass
        else:
            # get slider values
            self.updateFinePosition()
            self._lFinePos.SetLabel('Fine position = ({0[0]}, {0[1]}, {0[2]}'.format(self.TS.finecoords))

            # image center position
            XC = self.calX * self.imageData.shape[0] / 2.0
            YC = self.calY * self.imageData.shape[1] / 2.0

            print("XC/YC")
            print(XC)
            print(YC)

            # required motion
            XM = event.xdata - XC
            YM = event.ydata - YC

            print("event data")
            print(event.xdata)
            print(event.ydata)

            print("XM/YM")
            print(XM)
            print(YM)

            # Project onto TEAM Stage axes
            rot = self.TS.STEMRot / 180.0 * np.pi  # change to radians
            XMr = XM * np.cos(rot) + YM * np.sin(rot)
            YMr = -XM * np.sin(rot) + YM * np.cos(rot)

            print("XMr/YMr")
            print(XMr)
            print(YMr)

            # Determine X&ZY slider values (-100 to 100)
            newSliderX = self.TS.finecoords[0] + XMr / self.TS.sliderX_cal
            newSliderY = self.TS.finecoords[1] + YMr / self.TS.sliderY_cal

            # Test if in range else show error message and don't move (later move coarse?)
            newSliders = [newSliderX, newSliderY, self.TS.finecoords[2]]

            if ((newSliderX > -100) and (newSliderX < 100) and (newSliderY > -100) and (newSliderY < 100)):
                # Make a tuple of the new values
                self.setFinePosition(newSliders)
                self.sb.SetStatusText(
                    'New sliders: ' + str(newSliders[0]) + ", " + str(newSliders[1]) + ", " + str(newSliders[2]))
                self._lFinePos.SetLabel('Xf = ' + "{0:.2f}".format(newSliders[0]) + ", " + 'Yf = ' + "{0:.2f}".format(
                    newSliders[1]) + ", " + 'Zf = ' + "{0:.2f}".format(newSliders[2]))
                # wait 1 second
                time.sleep(1)
                self.quickimage()
            else:
                steps = [1, 2, 4]
                steps[0] = int(XM / self.TS.step10X_cal)
                steps[1] = int(YM / self.TS.step10Y_cal)
                steps[2] = 0
                size = [1, 1, 1]
                self.sb.SetStatusText('Steps: ' + str(steps[0]) + ', ' + str(steps[1]) + ', ' + str(steps[2]))
                # self.step(steps, size)
                time.sleep(1)
                self.quickimage()

            pass
        pass

    def step(self, steps, size):
        self.TS.TS_Connect('localhost', 5557)
        xsteps = steps[0]
        if xsteps < 0:
            self.TS.TS_Step([-xsteps, 0, 0], -size)
        elif xsteps > 0:
            self.TS.TS_Step([xsteps, 0, 0], size)
        ysteps = steps[1]
        if ysteps < 0:
            self.TS.TS_Step([-ysteps, 0, 0], -size)
        elif ysteps > 0:
            self.TS.TS_Step([ysteps, 0, 0], size)
        zsteps = steps[1]
        if zsteps < 0:
            self.TS.TS_Step([-zsteps, 0, 0], -size)
        elif zsteps > 0:
            self.TS.TS_Step([zsteps, 0, 0], size)
        self.TS.TS_Disconnect()
        return 1

    def updateFinePosition(self):
        # Get the TEAM stage fine position coordinates
        self.TS.TS_Connect('localhost', 5557)
        self.TS.TS_GetFinePosition()
        self.TS.TS_Disconnect()
        return self.TS.finecoords

    def onCalcTilt(self, event):
        self.calcgammatiltangles(int(self._islope.GetValue()),
                                 (float(self._iA1.GetValue()), float(self._iG1.GetValue())),
                                 (float(self._iA2.GetValue()), float(self._iG2.GetValue())))
        self.curtiltnr = 0

    def onNextTilt(self, event):
        if self.curtiltnr < self.tiltanglelist.shape[1] - 1:
            self.curtiltnr = self.curtiltnr + 1
        self.showtiltdetails()

    def onPrevTilt(self, event):
        if self.curtiltnr > 0:
            self.curtiltnr = self.curtiltnr - 1
        self.showtiltdetails()

    def showtiltdetails(self):
        #mystr = "nr= " + str(self.curtiltnr) + '  phi= ' + "{0:.2f}".format(
        #    self.tiltanglelist[0, self.curtiltnr]) + '  al= ' + "{0:.2f}".format(
        #    self.tiltanglelist[1, self.curtiltnr]) + '  ga= ' + "{0:.2f}".format(
        #    self.tiltanglelist[2, self.curtiltnr]) + '  rot= ' + "{0:.2f}".format(self.tiltanglelist[3, self.curtiltnr])
        mystr = 'nr= {} phi= {:0.2f} al= {:0.2f} ga= {:0.2f} rot= {:0.2f}'.format(self.curtiltnr,
                                                                                 self.tiltanglelist[0, self.curtiltnr],
                                                                                 self.tiltanglelist[1, self.curtiltnr],
                                                                                 self.tiltanglelist[2, self.curtiltnr],
                                                                                 self.tiltanglelist[3, self.curtiltnr])
        self._lselectedtilt.SetLabel(mystr)

    def calcgammatiltangles(self, slope, algam1, algam2):
        # For gamma tilting in TEAM Stage uses EST angles
        alpha1 = algam1[0]
        gamma1 = algam1[1]
        alpha2 = algam2[0]
        gamma2 = algam2[1]
        dir = 1.0
        if gamma2 < gamma1:
            dir = -1.0
        angles = range(slope * 4 + 1)
        for i in xrange(0, slope):
            angles[i] = math.atan(i / (slope * 1.0))
            angles[slope] = math.pi / 4
        for i in xrange(0, slope):
            angles[4 * slope - i - 1] = math.pi - angles[i + 1]
        for i in xrange(0, slope):
            angles[2 * slope + i] = math.pi / 2 + angles[i]
        for i in xrange(0, slope):
            angles[2 * slope - i] = math.pi - angles[i + 2 * slope]
        angles[4 * slope] = math.pi

        alphaoffsetdegrees = (alpha1 + alpha2 - 180) / 2
        axistiltdegrees = (alpha1 - alpha2) / 2

        gammafactor = (gamma1 - gamma2) / 180
        alphaoffset = alphaoffsetdegrees * math.pi / 180
        axistilt = axistiltdegrees * math.pi / 180

        print('alphaoffset = {}'.format(alphaoffsetdegrees))
        print('axistilt = {}'.format(axistiltdegrees))
        print('gammafactor = {}'.format(gammafactor))
        print('1/slope second angle = {}'.format(slope))
        print(' ')

        print('n \t phi \t alpha \t gamma')

        self.tiltanglelist = np.empty((4, slope * 4 + 1))
        self.tiltanglelist.fill(0)
        for n in xrange(slope * 4 + 1):
            phi = angles[n]
            phidegrees = phi / math.pi * 180
            alpha = math.acos(-math.sin(axistilt) * math.cos(phi))
            gamma = math.acos(math.cos(phi) * math.cos(axistilt) / math.sin(alpha))
            gammadegrees = gamma * 180 / math.pi
            alphadegrees = alpha * 180 / math.pi + alphaoffsetdegrees
            expgammadegrees = dir * gammadegrees * gammafactor + gamma1
            imagerotdegrees = phidegrees * 5  # insert correct function!!!!
            #print(str(n) + '\t' "{0:.2f}".format(phidegrees) + '\t' + "{0:.2f}".format(
            #    alphadegrees) + '\t' + "{0:.2f}".format(expgammadegrees))
            print('{} \t {:.2f} \t {:0.2f} \t {:0.2f}'.format(n, phidegrees, alphadegrees, expgammadegrees))
            self.tiltanglelist[0, n] = phidegrees
            self.tiltanglelist[1, n] = alphadegrees
            self.tiltanglelist[2, n] = expgammadegrees
            self.tiltanglelist[3, n] = imagerotdegrees

        self.curtiltnr = 0
        self.showtiltdetails()

    def getPosition_compustage(self):
        """ Get compustage position

        Returns
        -------
            : tuple
                The (X, Y, Z, alpha, beta) position and tilt angles.

        """
        stageObj = self.Stage.Position
        posArray = (stageObj.X, stageObj.Y, stageObj.Z,
                    stageObj.A, stageObj.B)
        return posArray

    def getPosition_teamstage(self):
        """ Return the current TEAM Stage position.

        Returns
        -------
            : tuple
                The (X, Y, Z, alpha, gamma) position and tilt angles.

        """
        self.TS.TS_Connect('localhost', 5557)
        self.TS.TS_GetPosition()
        self.TS.TS_Disconnect()
        return self.TS.coords

    def getFinePosition(self):
        """ Get the TEAM stage fine (slider) position coordinates

        """
        self.TS.TS_Connect('localhost', 5557)
        self.TS.TS_GetFinePosition()
        self.TS.TS_Disconnect()
        return self.TS.finecoords

    def setFinePosition(self, newSliders):
        """ Set the stage coordinates.

        """
        self.TS.TS_Connect('localhost', 5557)
        self.TS.TS_SetFinePosition(newSliders)
        self.TS.TS_Disconnect()
        return self.TS.finecoords

    def updatePositionSIM(self):
        """ For testing offline

        """
        print("using simulated TEAM Stage Positions")

        return (0, 1, 2, 3, 4)

    def updateEucentric(self):
        """ Get the TEAM Stage eucentric parameters

        """
        self.TS.TS_Connect('localhost', 5557)
        self.TS.TS_GetEucentricAxis()
        self.TS.TS_Disconnect()

    def updatePositionandFinePosition(self):
        self.TS.TS_Connect('localhost', 5557)
        self.TS.TS_GetPosition()
        self.TS.TS_GetFinePosition()
        self.TS.TS_Disconnect()

    def getPixelSize(self):
        """ Get the pixel calibration from TIA

        """
        window1 = self.TIA.ActiveDisplayWindow()
        Im1 = window1.FindDisplay(window1.DisplayNames[0])  # returns an image display object
        unit1 = Im1.SpatialUnit  # returns SpatialUnit object
        self.unitName = unit1.unitstring  # returns a string (such as nm)
        self.calX = Im1.image.calibration.deltaX  # returns the x calibration in meters
        self.calY = Im1.image.calibration.deltaY

    def onSearch(self, event):
        self.quickimage()

    def onGetPos(self, event):
        """ get position and slider values (if available) and
        show them in the GUI.

        """
        pos = self.getPosition()
        try:
            self._lFinePos.SetLabel('Xf = {}, Yf = {}, Zf = {}'.format(self.TS.finecoords[0], self.TS.finecoords[1]),
                                    self.TS.finecoords[2])
        except:
            pass
        self._lPos.SetLabel('(X, Y, Z) = ({0[0]}, {0[1]}, {0[2]})'.format(pos))
        self._lAlphaGamma.SetLabel('Alpha = {0[3]}, Gamma or Beta = {0[4]}'.format(pos))

    def quickimage(self):
        """ Acquire an image quickly.

        """

        # Stop an ongoing acquisition
        self.stopAcqusition()

        # Acquire a quick image
        self.setStemSearchVals()
        acquiredImageSet = self.Acq.AcquireImages()

        with safearray_as_ndarray:
            self.imageData = acquiredImageSet(0).AsSafeArray  # get data as ndarray

        self.getPixelSize()  # get the pixel calibration
        imageShape = self.imageData.shape
        self.axes.clear()
        self.axes.imshow(np.fliplr(np.rot90(self.imageData, 3)), extent=(0, self.calX * imageShape[0],
                                                                         0, self.calY * imageShape[1]))
        self.canvas.draw()  # draw the panel
        self.sb.SetStatusText('Search image')

        # self.ZoomimageData = self.imageData[int(imageShape[0] * 0.4):int(imageShape[0] * 0.6),
                             # int(imageShape[1] * 0.4):int(imageShape[1] * 0.6)]
        # ZoomimageShape = self.ZoomimageData.shape
        # self.Zoomaxes.clear()
        # self.Zoomaxes.imshow(np.fliplr(np.rot90(self.ZoomimageData, 3)), extent=(0, self.calX * ZoomimageShape[0],
                                                                                 # 0, self.calY * ZoomimageShape[1]))
        # self.Zoomcanvas.draw()  # draw the panel
        self._lDispFinePos.SetLabel(self._lFinePos.GetLabel())
        self._lDispPos.SetLabel(self._lPos.GetLabel())
        self._lDispAlphaGamma.SetLabel(self._lAlphaGamma.GetLabel())
        self._lDispIndex.SetLabel('Image shown: search (not saved) ')

        # Show FFT
        # if self._cbfft.GetValue():
            # fftimageShape = self.imageData.shape
            # cft = np.fft.fft2(self.imageData, self.imageData.shape)
            # vecfunc = np.vectorize(abs)
            # vecfunc2 = np.vectorize(math.log10)
            # vecfunc3 = np.vectorize(int)
            # a = vecfunc3(vecfunc2(vecfunc(cft)))
            # n = imageShape[0]
            # self.fftimageData = copy(a)
            # self.fftimageData[0:n / 2. - 1, 0:n / 2. - 1] = a[n / 2.:n - 1, n / 2.:n - 1]
            # self.fftimageData[n / 2.:n - 1, n / 2.:n - 1] = a[0:n / 2. - 1, 0:n / 2. - 1]
            # self.fftimageData[0:n / 2. - 1, n / 2.:n - 1] = a[n / 2.:n - 1, 0:n / 2. - 1]
            # self.fftimageData[n / 2.:n - 1, 0:n / 2. - 1] = a[0:n / 2. - 1, n / 2.:n - 1]
            # self.fftaxes.clear()
            # self.fftaxes.imshow(np.fliplr(np.rot90(self.fftimageData, 3)), extent=(0, self.calX * fftimageShape[0],
                                                                                   # 0, self.calY * fftimageShape[1]))
            # self.fftcanvas.draw()  # draw the panel
        # else:
            # self.fftaxes.clear()
            # self.fftcanvas.draw()  # draw the panel

        self.imagevalid = 1

    def determineMaxBinning(self):
        print('Determine maximum STEM binning')
        
        # Stop and ongiong acquisition
        self.stopAcqusition()
        self.Ill.BeamBlanked = True
        
        # Acquire an image
        self.setStemSearchVals()
        acquiredImageSet = self.Acq.AcquireImages()

        with safearray_as_ndarray:
                imageData = acquiredImageSet(0).AsSafeArray # get data as ndarray
        
        sh = imageData.shape
        
        self.maxBin = sh[0] * self.myStemSearchParams.Binning

    def onShowPrev(self, event):
        self.showindex = self.showindex - 1
        self.showimage()

    def onShowNext(self, event):
        self.showindex = self.showindex + 1
        self.showimage()

    def onGotoDispPos(self, event):
        self.stage_goto(self.DispPosangles)

    def onNewFileName(self, event):
        self.NewFileName()

    def NewFileName(self):
        self.Fprefix = self._iFprefix.GetValue()
        self.Fdir = self._iDir.GetValue()
        self.fullName = self.Fdir + os.sep + self.Fprefix
        if os.path.exists(self.fullName + '.emd'):
            self.setDatasetAttributesEditable(False)
            # copy all parameters from file
            f = h5py.File(self.fullName + '.emd', 'r')
            dataroot = f['data']
            dataTop = dataroot['raw']
            self.Rep = int(dataTop.attrs.get('images per set'))
            self.rotSetting = dataTop.attrs['rotSetting']
            self.Bin = int(dataTop.attrs.get('binning'))
            self.maxBin = int(dataTop.attrs.get('max binning'))
            self.Dwell = float(dataTop.attrs.get('dwell time'))
            self.Del = float(dataTop.attrs.get('delay time'))
            self._iRep.SetValue(str(self.Rep))
            self._iRot.SetSelection(self.rotSetting)
            self._iBin.SetValue(str(self.Bin))
            self._iDwell.SetValue(str(int(self.Dwell / 1e-6)))
            self._iDel.SetValue(str(int(self.Del)))
            self.setDatasetAttributesEditable(False)
            f.close()
        else:
            self.setDatasetAttributesEditable(True)

    def stage_goto(self, p):
        self.TS.TS_Connect('localhost', 5557)
        self.TS.TS_GoTo(p)
        self.TS.TS_Disconnect()

    def showimage(self):
        # Show image from data set
        with h5py.File(self.fullName + '.emd', 'a') as f:  # open and expand if file exists
            dataroot = f['data']
            dataTop = dataroot['raw']
            dset = dataTop['data']
            nindex = dataTop.attrs.get('nextindex')

            # check that in range
            if (self.showindex < dataTop.attrs.get('nextindex', 0)) and (self.showindex >= 0):
                self.showset = dataTop['imageParameters'][0, self.showindex]
                self.showsubindex = dataTop['imageParameters'][1, self.showindex]
                # self.imageData = dset[self.showsubindex,self.showset,:,:]
                self.imageData = dset[self.showset, self.showsubindex, :, :]
                self.calX = dataTop['dim4'][1] - dataTop['dim4'][0]
                self.calY = dataTop['dim3'][1] - dataTop['dim3'][0]
                imageShape = self.imageData.shape
                self.axes.clear()
                self.axes.imshow(np.fliplr(np.rot90(self.imageData, 3)), extent=(0, self.calX * imageShape[0],
                                                                                 0, self.calY * imageShape[1]))
                self.canvas.draw()  # draw the panel
                self._lDispIndex.SetLabel('Image shown: ' + str(self.showindex))
                # stage = f['stage']
                Position = dataTop['imageSetParameters'][0:5, self.showset]
                self._lDispFinePos.SetLabel('Fine position not recorded')
                self._lDispPos.SetLabel('(X, Y , Z) = ({0[0]}, {0[1]}, {0[2]})'.format(Position))
                self._lDispAlphaGamma.SetLabel('Alpha = {0[3]}, Gamma or Beta = {0[4]}'.format(Position))
                self.DispPosangles = Position[:]

                # self.ZoomimageData = self.imageData[int(imageShape[0] * 0.4):int(imageShape[0] * 0.6),
                                     # int(imageShape[1] * 0.4):int(imageShape[1] * 0.6)]
                # ZoomimageShape = self.ZoomimageData.shape
                # self.Zoomaxes.clear()
                # self.Zoomaxes.imshow(np.fliplr(np.rot90(self.ZoomimageData, 3)), extent=(0, self.calX * ZoomimageShape[0],
                                                                                         # 0, self.calY * ZoomimageShape[1]))
                # self.Zoomcanvas.draw()  # draw the panel

                # if self._cbfft.GetValue():
                    # fftimageShape = self.imageData.shape
                    # cft = np.fft.fft2(self.imageData, self.imageData.shape)
                    # vecfunc = np.vectorize(abs)
                    # vecfunc2 = np.vectorize(math.log10)
                    # vecfunc3 = np.vectorize(int)
                    # a = vecfunc3(vecfunc2(vecfunc(cft)))
                    # n = imageShape[0]
                    # self.fftimageData = copy(a)
                    # self.fftimageData[0:n / 2 - 1, 0:n / 2 - 1] = a[n / 2:n - 1, n / 2:n - 1]
                    # self.fftimageData[n / 2:n - 1, n / 2:n - 1] = a[0:n / 2 - 1, 0:n / 2 - 1]
                    # self.fftimageData[0:n / 2 - 1, n / 2:n - 1] = a[n / 2:n - 1, 0:n / 2 - 1]
                    # self.fftimageData[n / 2:n - 1, 0:n / 2 - 1] = a[0:n / 2 - 1, n / 2:n - 1]
                    # self.fftaxes.clear()
                    # self.fftaxes.imshow(np.fliplr(np.rot90(self.fftimageData, 3)), extent=(0, self.calX * fftimageShape[0],
                                                                                           # 0, self.calY * fftimageShape[1]))
                    # self.fftcanvas.draw()  # draw the panel
                # else:
                    # self.fftaxes.clear()
                    # self.fftaxes.imshow(np.fliplr(np.rot90(self.ZoomimageData, 3)), extent=(0, self.calX * ZoomimageShape[0],
                                                                                            # 0, self.calY * ZoomimageShape[1]))
                    # self.fftcanvas.draw()  # draw the panel

                self.imagevalid = 0

    def stopAcqusition(self):
        if self.TIA.AcquisitionManager().isAcquiring:
            self.TIA.AcquisitionManager().Stop()
    
    def onGotoAngles(self, event):
        """ Go to the angles that the user requested while keeping the x,y,z
        consistent.
        
        """
        alpha = float(self._gotoGUI['_iGotoA'].GetValue())
        gamma = float(self._gotoGUI['_iGotoG'].GetValue())
        cur_p = self.getPosition()
        p = (cur_p[0], cur_p[1], cur_p[2], alpha, gamma)
        if self._gotoGUI['_cGotoConfirm'].GetValue():
            answer = wx.MessageBox('Goto: X = {0[0]}, Y= {0[1]}, Z = {0[2]}, Alpha = {0[3]}, Gamma = {0[4]}?'.format(p),
                                   'Stage movement confirm', wx.YES_NO | wx.ICON_EXCLAMATION)
        else:
            answer = True
        
        if answer:
            self.stage_goto(p)
        
    def onAcquire(self, event):
        """ Acquire the image(s) for the tilt series stack. Save
        them to the EMD file after all acquisitions are complete
        to speed up the acquisition.

        """
        # Stop an ongoing acquisition
        self.stopAcqusition()

        # Save the data and acquisition information to HDF5/EMD
        # Expand the HDF5 dataset if exists. Otherwise create a new file
        newfile = 0
        self.getUsrAcqVals()
        if not (os.path.exists(self.fullName + '.emd')):
            self.initializeEMDFile()
            self.setDatasetAttributesEditable(False)
            newfile = 1
        else:
            self.NewFileName()

        self.setStemAcqVals()

        myImageSetStack = []
        mytimestack = []

        # Create a STEM rotation array
        rot0 = self.Ill.StemRotation
        rotNum = self.rotSetting
        if rotNum == 1:  # 2 rotations is selected
            rotArray = np.array((rot0, rot0 + np.pi / 2.0))  # 0, 90 pair
            self.Rep = 2
        elif rotNum == 2:  # 4 rotations is selected
            rotArray = np.array((rot0, rot0 + np.pi / 2.0, rot0 + np.pi / 4.0,
                                 rot0 + 3.0 * np.pi / 4.0))  # 0, 90 , 45, 135 degrees in that order
            self.Rep = 4
        else:
            rotArray = rot0 * np.ones(self.Rep)

        self._iRep.SetValue(str(self.Rep))

        # Acquire the images
        tempData = np.zeros((self.Rep, self.maxBin / self.Bin, self.maxBin / self.Bin))
        for x in range(0, self.Rep):
            if (x != 0):  # and (self.Del<10)):
                self.sb.SetStatusText('Waiting ' + str(self.Del) + ' seconds ...')
                time.sleep(self.Del)
            self.sb.SetStatusText('Acquiring image ' + str(x + 1) + ' ...')

            self.Ill.StemRotation = rotArray[x]

            # Acquire an image
            self.Ill.BeamBlanked = False
            mytimestack.append(time.time())
            acquiredImageSet = self.Acq.AcquireImages()
            self.Ill.BeamBlanked = True

            with safearray_as_ndarray:
                self.imageData = acquiredImageSet(0).AsSafeArray  # get data as ndarray

            tempData[x, :, :] = self.imageData  # store data to write to disk later
            imageShape = self.imageData.shape
            self.getPixelSize()

            # Show the image in the main window
            self.axes.clear()
            self.axes.imshow(np.fliplr(np.rot90(self.imageData, 3)), extent=(
            0, self.calX * imageShape[0], 0, self.calY * imageShape[1]))  # add the image to the figure
            self.canvas.draw()  # draw the panel
            self.imagevalid = 1

            # Show the zoom image
            # self.ZoomimageData = self.imageData[int(imageShape[0] * 0.4):int(imageShape[0] * 0.6),
                                 # int(imageShape[1] * 0.4):int(imageShape[1] * 0.6)]
            # ZoomimageShape = self.ZoomimageData.shape
            # self.Zoomaxes.clear()
            # self.Zoomaxes.imshow(np.fliplr(np.rot90(self.ZoomimageData, 3)), extent=(
            # 0, self.calX * ZoomimageShape[0], 0, self.calY * ZoomimageShape[1]))  # add the image to the figure
            # self.Zoomcanvas.draw()  # draw the panel
            # self.fftcanvas.draw()  # draw the panel
            
            self._lDispIndex.SetLabel('Image shown: from current aquired set nr. ' + str(x))
            self.sb.SetStatusText('Finished Acquisition image ' + str(x + 1) + ' ...')

        # Set STEM rotation back to the original value
        self.Ill.StemRotation = rot0

        # Save the data to the currently open file name
        with h5py.File(self.fullName + '.emd', 'a') as f:
            dataroot = f['data']
            dataTop = dataroot['raw']
            dset = dataTop['data']
            # stage = f['stage']
            user = f['user']

            if newfile == 1:
                print('Check EMD initialization is correct. This note replaces an untested removal of deprecated code.')
                xdim = np.linspace(0, (imageShape[0] - 1) * self.calX * 1e9,
                                   imageShape[0])  # multiply by 1e9 for nanometers
                ydim = np.linspace(0, (imageShape[1] - 1) * self.calY * 1e9, imageShape[1])

                dataTop['dim4'].resize(imageShape[0], 0)
                dataTop['dim3'].resize(imageShape[1], 0)
                dataTop['dim4'][:] = xdim
                dataTop['dim3'][:] = ydim

                dataTop['dim2'].resize(self.Rep, 0)
                if (self.rotSetting == 1) | (
                        self.rotSetting == 2):  # use rotation series angles if taking a rotation series
                    dataTop['dim2'][:] = rotArray * 180.0 / np.pi
                    dataTop['dim2'].attrs['name'] = np.string_('scan rotation')
                    dataTop['dim2'].attrs['units'] = np.string_('degrees')
                else:
                    dataTop['dim2'][:] = range(0, self.Rep)

            # Data common to set
            pos = self.getPosition()
            nset = dataTop.attrs.get('nextset')
            dataTop.attrs['nextset'] = nset + 1
            dataTop['dim1'].resize(nset + 1, 0)
            dataTop['dim1'][nset] = pos[3]  # Use Alpha tilt angle in dim1
            dataTop['tiltangles'].resize((nset + 1,))
            dataTop['tiltangles'][-1] = pos[3]  # Save the current tilt angle in an easily readable array
            dataTop['imageSetParameters'].resize(nset + 1, 1)  # resize the second dimension
            dataTop['imageSetParameters'][0:5, nset] = pos  # save the current position
            dataTop['imageSetParameters'][5, nset] = mytimestack[0]  # save the acquisition time

            # Write the data to disk
            for x in range(0, self.Rep):
                self.sb.SetStatusText('Writing image(s) ' + str(x + 1) + ' ...')

                self.imageData = np.squeeze(tempData[x, :, :])  # write temporary data to disk

                nindex = dataTop.attrs.get('nextindex')
                dataTop.attrs['nextindex'] = nindex + 1
                dset.resize(nset + 1, 0)  # increase the alpha axis data shape
                dset[nset, x, :, :] = self.imageData

                dataTop['imageParameters'].resize(nindex + 1, 1)  # resize the second dimension
                dataTop['imageParameters'][0, nindex] = nset
                dataTop['imageParameters'][1, nindex] = x
                dataTop['imageParameters'][2, nindex] = rotArray[x]
                dataTop['imageParameters'][3, nindex] = self.calX
                dataTop['imageParameters'][4, nindex] = self.calY
                dataTop['imageParameters'][5, nindex] = mytimestack[x]

            f['microscope'].attrs['STEM rotations'] = rotArray
        # f.close()

        self.sb.SetStatusText('Idle...')

        self._lDispIndex.SetLabel('Image shown: ' + str(nindex))
        self.showindex = nindex


# Parse the arguments
parser = argparse.ArgumentParser(description='Acquire atomic electron tomography tilt series.')

# Allow user to set TeamStage or Compustage
parser.add_argument('--teamstage', '-t', action='store_const', const=True, default=False)
parser.add_argument('--compustage', '-c', action='store_const', const=True, default=True)

parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(version))

args = parser.parse_args()

if args.teamstage:
    stagetype0 = 'teamstage'
else:
    stagetype0 = 'compustage'

print('Stage type = {}'.format(stagetype0))

app = wx.App(False)
frame = TEAMFrame(None, 'NCEM: STEM Tomography v{}'.format(version), stagetype=stagetype0)
app.MainLoop()
