''' Python script for acquiring STEM data experiments
such as Focal series, Scanning rotation series, time series and
single images. All data is saved as a Berkeley EMD file. Use
ncempy.io.emd to read these files.

author: Peter Ercius, percius@lbl.gov
'''

import uuid
import argparse
version = 1.1 #version number for this program

import numpy as np
import socket
import os
import wx # wxPython GUI package
import h5py
import time

# For connections to FEI TEMScripting and TIA
from comtypes.client import CreateObject
from comtypes.safearray import safearray_as_ndarray # get data across COM barrier fast

# For connecting to TEAM Stage
try:
    import TEAMstageclass
except:
    print('No TEAM Stage functions available.')

# Main window
class TEAMFrame(wx.Frame):
    def __init__(self, parent, title, stagetype = 'compustage'):
        #Import only necessary modules from matplotlib

        from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
        from matplotlib.figure import Figure
        try:
            from matplotlib.pyplot import set_cmap
        except:
            pass
        # Setup the stagetype as compustage or teamstage
        self.stagetype = stagetype
        
        self.showindex = 0
        self.showsubindex = 0
        self.showset = 0
        self.curtiltnr = 0
        self.maxBin = 0
        
        self.II = 0 # integer to put after file names
        self.times = None # acquisition time(s)
        
        self.calX = self.calY = 0
        
        # Set the update stage position to either
        # compustage or teamstage function
        if stagetype == 'teamstage':
            self.getPosition = self.getPosition_teamstage
            self.microscopeName = 'TEAM 0.5'
            self.stageName = 'NCEM TEAM Stage'
        elif stagetype == 'compustage':
            self.getPosition = self.getPosition_compustage
            self.microscopeName = 'FEI'
            self.stageName = 'compustage'
        
        #Initialize the base Frame
        wx.Frame.__init__(self, parent, title=title, size=(1300,1000))
        
        # menu 
        # menubar = wx.MenuBar()
        # filem = wx.Menu()
        # editm = wx.Menu()
        # helpm = wx.Menu()
        # menubar.Append(filem, '&File')
        # menubar.Append(editm, '&Edit')
        # menubar.Append(helpm, '&Help')
        # self.SetMenuBar(menubar)
        
        # Define sizers   
        mainbox = wx.BoxSizer(wx.HORIZONTAL)
        imagecolumnbox = wx.BoxSizer(wx.VERTICAL)
        infobox = wx.BoxSizer(wx.VERTICAL)
        imagebox = wx.BoxSizer(wx.HORIZONTAL)
        displayvalbox = wx.BoxSizer(wx.VERTICAL)
        imagebuttonbox = wx.BoxSizer(wx.HORIZONTAL)
        currentvalbox = wx.BoxSizer(wx.VERTICAL)
        setupbox = wx.BoxSizer(wx.VERTICAL)
        acquirebox = wx.BoxSizer(wx.VERTICAL)
        focalbox = wx.BoxSizer(wx.VERTICAL)
        driftbox = wx.BoxSizer(wx.VERTICAL)
        seriesbox = wx.BoxSizer(wx.VERTICAL)

        notebook = wx.Notebook(self)
        
        # create the page windows as children of the notebook
        setupPage = wx.Panel(notebook)
        acquirePage = wx.Panel(notebook) # single acquires
        focalPage = wx.Panel(notebook) # focal series
        driftPage = wx.Panel(notebook) # 0-90 degree drift compensation
        timeSeriesPage = wx.Panel(notebook) # time series
        reviewpage = wx.Panel(notebook) # review previous data
        
        # Figure
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
            
        self.axesIm = self.axes.imshow(np.zeros((10,10)))
        self.TEMfigure.tight_layout()
        
        # Setup page
        self._lDir = wx.StaticText(setupPage, wx.ID_ANY, 'Output directory')
        try:
            self._iDir = wx.TextCtrl(setupPage,wx.ID_ANY, value=r'G:\UserData\Ercius\temp', size=(300,20)) #correct for TEAM 0.5
        except:
            self._iDir = wx.TextCtrl(setupPage,wx.ID_ANY, value=os.getcwd(), size=(300,20))
        
        self._lUser = wx.StaticText(setupPage, wx.ID_ANY, 'User name')
        self._iUser = wx.TextCtrl(setupPage,wx.ID_ANY, value='', size=(300,20))
        self._lSample = wx.StaticText(setupPage, wx.ID_ANY, 'Sample name')
        self._iSample = wx.TextCtrl(setupPage, wx.ID_ANY, value='', size=(300,20))
        setupbox.Add(self._lDir)
        setupbox.Add(self._iDir)
        setupbox.Add(self._lUser)
        setupbox.Add(self._iUser)
        setupbox.Add(self._lSample)
        setupbox.Add(self._iSample)
        
        # Single image acquisition page
        self._lLabel00 = wx.StaticText(acquirePage, wx.ID_ANY, 'File Parameters')
        self._lFprefix0 = wx.StaticText(acquirePage, wx.ID_ANY, 'File name prefix')
        self._iFprefix0 = wx.TextCtrl(acquirePage,wx.ID_ANY,value='Acquisition')
        
        self._lLabel01 = wx.StaticText(acquirePage, wx.ID_ANY, 'Acquire parameters')
        self._lDwell0 = wx.StaticText(acquirePage, wx.ID_ANY, 'Dwell time (usec)')
        self._iDwell0 = wx.TextCtrl(acquirePage,wx.ID_ANY,value='1')
        self._lBin0 = wx.StaticText(acquirePage, wx.ID_ANY, 'Binning')
        self._iBin0 = wx.TextCtrl(acquirePage,wx.ID_ANY,value='8')
        self._btnAcqSingle = wx.Button(acquirePage, label='Acquire')
        acquirebox.Add(self._lLabel00)
        acquirebox.Add(self._lFprefix0)
        acquirebox.Add(self._iFprefix0)
        acquirebox.Add(self._lLabel01)
        acquirebox.Add(self._lDwell0)
        acquirebox.Add(self._iDwell0)
        acquirebox.Add(self._lBin0)
        acquirebox.Add(self._iBin0)
        acquirebox.Add(self._btnAcqSingle)
        
        # Focal series acquisition page
        self._lFprefix1 = wx.StaticText(focalPage, wx.ID_ANY, 'File name prefix')
        self._iFprefix1 = wx.TextCtrl(focalPage,wx.ID_ANY,value='FocalSeries')
        
        self._lLabel10 = wx.StaticText(focalPage, wx.ID_ANY, 'Acquire parameters')
        self._lDwell1 = wx.StaticText(focalPage, wx.ID_ANY, 'Dwell time (usec)')
        self._iDwell1 = wx.TextCtrl(focalPage,wx.ID_ANY,value='1')
        self._lBin1 = wx.StaticText(focalPage, wx.ID_ANY, 'Binning')
        self._iBin1 = wx.TextCtrl(focalPage,wx.ID_ANY,value='8')
        self._lRep1 = wx.StaticText(focalPage, wx.ID_ANY, 'Number of images per aquire')
        self._iRep1 = wx.TextCtrl(focalPage,wx.ID_ANY,value='1')
        self._lStart1 = wx.StaticText(focalPage, wx.ID_ANY, 'Start defocus (nm)')
        self._iStart1 = wx.TextCtrl(focalPage, wx.ID_ANY,value='-40')
        self._lStep1 = wx.StaticText(focalPage, wx.ID_ANY, 'Focal step (nm)')
        self._iStep1 = wx.TextCtrl(focalPage,wx.ID_ANY,value='10')
        self._lNum1 = wx.StaticText(focalPage, wx.ID_ANY, 'Number of defoci')
        self._iNum1 = wx.TextCtrl(focalPage,wx.ID_ANY,value='9')
        self._btnAcqFocal = wx.Button(focalPage, label='Acquire')
        focalbox.Add(self._lFprefix1)
        focalbox.Add(self._iFprefix1)
        focalbox.Add(self._lLabel10)
        focalbox.Add(self._lDwell1)
        focalbox.Add(self._iDwell1)
        focalbox.Add(self._lBin1)
        focalbox.Add(self._iBin1)
        focalbox.Add(self._lStart1)
        focalbox.Add(self._iStart1)
        focalbox.Add(self._lStep1)
        focalbox.Add(self._iStep1)
        focalbox.Add(self._lNum1)
        focalbox.Add(self._iNum1)
        focalbox.Add(self._lRep1)
        focalbox.Add(self._iRep1)
        focalbox.Add(self._btnAcqFocal)
        
        # Drift series acquisition page
        self._lFprefix2 = wx.StaticText(driftPage, wx.ID_ANY, 'File name prefix')
        self._iFprefix2 = wx.TextCtrl(driftPage,wx.ID_ANY,value='Rotation')
        
        self._lLabel20 = wx.StaticText(driftPage, wx.ID_ANY, 'Acquire parameters')
        self._lDwell2 = wx.StaticText(driftPage, wx.ID_ANY, 'Dwell time (usec)')
        self._iDwell2 = wx.TextCtrl(driftPage,wx.ID_ANY,value='1')
        self._lBin2 = wx.StaticText(driftPage, wx.ID_ANY, 'Binning')
        self._iBin2 = wx.TextCtrl(driftPage,wx.ID_ANY,value='8')
        self._lRot2 = wx.StaticText(driftPage,wx.ID_ANY, 'Num scan rotations')
        self._iRot2 = wx.RadioBox(driftPage,wx.ID_ANY,label='Rotation Series (degrees)', choices=('0, 90','0, 90, 45, 135'))
        self._iRot2.SetSelection(0)
        self._btnAcqDrift = wx.Button(driftPage, label = 'Acquire')
        
        driftbox.Add(self._lLabel20)
        driftbox.Add(self._lFprefix2)
        driftbox.Add(self._iFprefix2)
        driftbox.Add(self._lDwell2)
        driftbox.Add(self._iDwell2)
        driftbox.Add(self._lBin2)
        driftbox.Add(self._iBin2)
        driftbox.Add(self._lRot2)
        driftbox.Add(self._iRot2)
        driftbox.Add(self._btnAcqDrift)
        
        # Time series acquisition page
        self._lFprefix3 = wx.StaticText(timeSeriesPage, wx.ID_ANY, 'File name prefix')
        self._iFprefix3 = wx.TextCtrl(timeSeriesPage,wx.ID_ANY,value='TimeSeries')
        
        self._lLabel30 = wx.StaticText(timeSeriesPage, wx.ID_ANY, 'Acquire parameters')
        self._lDwell3 = wx.StaticText(timeSeriesPage, wx.ID_ANY, 'Dwell time (usec)')
        self._iDwell3 = wx.TextCtrl(timeSeriesPage,wx.ID_ANY,value='1')
        self._lBin3 = wx.StaticText(timeSeriesPage, wx.ID_ANY, 'Binning')
        self._iBin3 = wx.TextCtrl(timeSeriesPage,wx.ID_ANY,value='8')
        self._lRep3 = wx.StaticText(timeSeriesPage, wx.ID_ANY, 'Number of images per aquire')
        self._iRep3 = wx.TextCtrl(timeSeriesPage,wx.ID_ANY,value='3')
        self._lDel3 = wx.StaticText(timeSeriesPage, wx.ID_ANY, 'Delay between images (sec)')
        self._iDel3 = wx.TextCtrl(timeSeriesPage,wx.ID_ANY,value='0')
        
        self._btnAcqTimeSeries = wx.Button(timeSeriesPage, label='Acquire')
        seriesbox.Add(self._lFprefix3)
        seriesbox.Add(self._iFprefix3)
        seriesbox.Add(self._lDwell3)
        seriesbox.Add(self._iDwell3)
        seriesbox.Add(self._lBin3)
        seriesbox.Add(self._iBin3)
        seriesbox.Add(self._lRep3)
        seriesbox.Add(self._iRep3)
        seriesbox.Add(self._lDel3)
        seriesbox.Add(self._iDel3)
        seriesbox.Add(self._btnAcqTimeSeries)
        
        # Review page
        '''
        self._lDispPos = wx.StaticText(reviewpage, wx.ID_ANY, 'Displ pos')
        self._lDispAlphaGamma = wx.StaticText(reviewpage, wx.ID_ANY, 'Displ alpha gamma')
        self._lDispFinePos = wx.StaticText(reviewpage, wx.ID_ANY, 'Displ fine pos')
        self._btnDispGotoPos = wx.Button(reviewpage, label='Goto display position')     
        self._lDispIndex  = wx.StaticText(reviewpage, wx.ID_ANY, 'Index')
        displayvalbox.Add(self._lDispPos)
        displayvalbox.Add(self._lDispAlphaGamma)
        displayvalbox.Add(self._lDispFinePos)
        displayvalbox.Add(self._lDispIndex)
        displayvalbox.Add(self._btnDispGotoPos, proportion=0, flag=wx.ALIGN_RIGHT|wx.RIGHT, border=10)
        self._btnprev = wx.Button(reviewpage, label='Show prev')
        self._btnnext = wx.Button(reviewpage, label='Show next')
        imagebuttonbox.Add(self._btnprev, proportion=0, flag=wx.ALIGN_RIGHT|wx.RIGHT, border=10)
        imagebuttonbox.Add(self._btnnext, proportion=0, flag=wx.ALIGN_RIGHT|wx.RIGHT, border=10)
        displayvalbox.Add(imagebuttonbox, proportion=0, flag=wx.LEFT|wx.RIGHT|wx.EXPAND, border=20)
        '''
        
        #Set color
        self.SetBackgroundColour("0000FF")
        
        self.imagevalid = 0
        
        #Bind the events
        self._btnAcqSingle.Bind(wx.EVT_BUTTON, self.onAcquireSingle) #bind the button event
        self._btnAcqFocal.Bind(wx.EVT_BUTTON, self.onAcquireFocal) #bind the button event
        self._btnAcqTimeSeries.Bind(wx.EVT_BUTTON, self.onAcquireTimeSeries) #bind the button event
        self._btnAcqDrift.Bind(wx.EVT_BUTTON, self.onAcquireDrift) #bind the button event

        # Eanble Mouse events
        #self.canvas.mpl_connect('button_press_event', self.onClick) #bind the left click event
                
        #self._btnprev.Bind(wx.EVT_BUTTON, self.onShowPrev) #bind the button event
        #self._btnnext.Bind(wx.EVT_BUTTON, self.onShowNext) #bind the button event
        
        #self._btnDispGotoPos.Bind(wx.EVT_BUTTON, self.onGotoDispPos)
        
        imagebox.Add(self.canvas, proportion=3, flag=wx.LEFT|wx.RIGHT|wx.EXPAND, border=0)
                
        imagecolumnbox.Add(imagebox, proportion=1, flag=wx.LEFT|wx.RIGHT|wx.EXPAND, border=0)
        
        infobox.Add(currentvalbox, proportion=0, flag=wx.LEFT|wx.RIGHT|wx.EXPAND, border=0)
        
        setupPage.SetSizer(setupbox)
        acquirePage.SetSizer(acquirebox)
        focalPage.SetSizer(focalbox)
        driftPage.SetSizer(driftbox)
        timeSeriesPage.SetSizer(seriesbox)
        #reviewpage.SetSizer(displayvalbox)
        
        # add the pages to the notebook with the label to show on the tab
        notebook.AddPage(setupPage, 'Setup')
        notebook.AddPage(acquirePage, 'Single')
        notebook.AddPage(focalPage, 'Focal')
        notebook.AddPage(driftPage, 'Rotation')
        notebook.AddPage(timeSeriesPage, 'Time series')
        #notebook.AddPage(reviewpage, 'Review')
        
        mainbox.Add(imagecolumnbox, proportion=1, flag=wx.LEFT|wx.RIGHT|wx.EXPAND, border=0)
        mainbox.Add(notebook, proportion=0, flag=wx.LEFT|wx.RIGHT|wx.EXPAND, border=20)
        
        self.SetSizer(mainbox) #Initialize the parent sizer
        
        self.sb = self.CreateStatusBar()
        self.sb.SetStatusText('Idle...')
        
        self.onConnect()
        
        self.Show(True) #show the frame
    
    def __del__(self):
        self.f.close()
        self.TS.TS_Disconnect()
        
    def drawPlot(self, imArray):
        #print('Draw plot')
        
        # Need to add updates to the color limits. Not sure how to do that.
        #self.axesIm.set_array(imArray) # This should be faster...but it does not work. Need to update 

        self.axes.clear()
        if imArray.shape[0] <= 2048:
            axIm = self.axes.imshow(imArray) 
        else:
            axIm = self.axes.imshow(imArray[::2,::2]) # 4kx4k cant be shown. reduce the size of the image
        self.canvas.draw() #draw the panel
    
    #Connect to the microscope and setup STEM detectors
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
        
        if self.stagetype == 'teamstage':
            #Connect to TEAM stage
            try:
                self.TS = TEAMstageclass.TEAMstage()
                self.TS.TS_Connect('localhost',5557)
            except Exception as e:
                print('Error with TEAM Stage connection. Exception type is %s' % e)
            
        try:
            # Get microscope interfaces
            self.Acq = self._microscope.Acquisition
            self.Proj = self._microscope.Projection
            self.Ill = self._microscope.Illumination
            self.Stage = self._microscope.Stage #FEI stage. Not needed for TEAM Stage
            
            self.detector0 = self.Acq.Detectors(0) #older pythoncom versions might require square brackets []
        except:
            print('Connections to microscope interfaces failed.')
            raise

        # Add the detector
        self._microscope.Acquisition.AddAcqDevice(self.detector0)
        
        self.determineMaxBinning()
        print('Max binning = {}'.format(self.maxBin))

    def setFonts(self,ed):
        if ed:
            fontTopLabels = wx.Font(14, wx.DEFAULT, wx.NORMAL,wx.BOLD)
            fontLabels = wx.Font(12, wx.DEFAULT, wx.NORMAL,wx.BOLD)
        else:
            fontTopLabels = wx.Font(14, wx.DEFAULT, wx.ITALIC,wx.BOLD)
            fontLabels = wx.Font(12, wx.DEFAULT, wx.ITALIC,wx.BOLD)

        self._lLabel1.SetFont(fontTopLabels)
        self._lLabel3.SetFont(fontTopLabels)
        self._lFprefix.SetFont(fontLabels)
        self._lDir.SetFont(fontLabels)
        self._lDwell.SetFont(fontLabels)
        self._lBin.SetFont(fontLabels)
        self._lRep.SetFont(fontLabels)
        self._lRot.SetFont(fontLabels)
        self._lDel.SetFont(fontLabels)
    
    ''' OLD
    def createDims(self, dataTop, type, pix):
        if type == 'single':
            dim2 = dataTop.create_dataset('dim2',(pix,),'f')
            dim2.attrs['name'] = np.string_('X')
            dim2.attrs['units'] = np.string_('n_m')
            dim1 = dataTop.create_dataset('dim1',(pix,),'f')
            dim1.attrs['name'] = np.string_('Y')
            dim1.attrs['units'] = np.string_('n_m')
            dims = (dim1,dim2)
        else:
            dims = None
        return dims
    '''
    
    def setStemAcqVals(self):
        self.myStemAcqParams = self.Acq.Detectors.AcqParams
        self.myStemAcqParams.Binning = self.Bin
        self.myStemAcqParams.ImageSize = 0 #self._microscope.ACQIMAGESIZE_FULL
        self.myStemAcqParams.DwellTime = self.Dwell
        self.Acq.Detectors.AcqParams = self.myStemAcqParams
    
    def setStemSearchVals(self):
        self.myStemSearchParams = self.Acq.Detectors.AcqParams
        self.myStemSearchParams.Binning = 8
        self.myStemSearchParams.ImageSize = 0 # self._microscope.ACQIMAGESIZE_FULL
        self.myStemSearchParams.DwellTime = 2e-6
        self.Acq.Detectors.AcqParams = self.myStemSearchParams
        
    def getUsrAcqVals(self, exType):
        self.Fdir = self._iDir.GetValue()
        
        if exType == 'drift':
            self.rotSetting = self._iRot2.GetSelection()
            self.Fprefix = self._iFprefix2.GetValue()
            self.fullName = self.Fdir + os.sep + self.Fprefix
            self.Dwell = float(self._iDwell2.GetValue()) * 1e-6 #Change to microseconds
            self.Bin = int(self._iBin2.GetValue())
        elif exType == 'timeseries':
            self.Fprefix = self._iFprefix3.GetValue()
            self.fullName = self.Fdir + os.sep + self.Fprefix
            self.Rep = int(self._iRep3.GetValue())
            self.Del = int(self._iDel3.GetValue())
            self.Dwell = float(self._iDwell3.GetValue()) * 1e-6 #Change to microseconds
            self.Bin = int(self._iBin3.GetValue())
            print('timeseries repeats = {}'.format(self.Rep))
        elif exType == 'single':
            self.Fprefix = self._iFprefix0.GetValue()
            self.fullName = self.Fdir + os.sep + self.Fprefix
            self.Dwell = float(self._iDwell0.GetValue()) * 1e-6 #Change to microseconds
            self.Bin = int(self._iBin0.GetValue())
        elif exType == 'focalseries':
            self.Fprefix = self._iFprefix1.GetValue()
            self.fullName = self.Fdir + os.sep + self.Fprefix
            self.Dwell = float(self._iDwell1.GetValue()) * 1e-6 #Change to microseconds
            self.Bin = int(self._iBin1.GetValue())
            self.numDF = int(self._iNum1.GetValue())
            self.numPerDF = int(self._iRep1.GetValue())
            self.stepDF = int(self._iStep1.GetValue()) * 1e-9 # change to meters
            
            self.Rep = int(self.numDF * self.numPerDF)
            
        self.expectedImageShape = (self.maxBin/self.Bin,self.maxBin/self.Bin)
        #print('self.expectedImageShape = {}'.format(self.expectedImageShape))
        
    def onClick(self,event):
        # Determine mouse position and move stage. xdata = x coord of mouse in data coords
        if ((self.imagevalid ==0) or (event.xdata == None) or (event.ydata == None)):
            pass
        else:
            #get slider values
            self.updateFinePosition()
            #self._lFinePos.SetLabel('XF = ' + str(self.TS.finecoords[0]) + ", " + 'YF = ' + str(self.TS.finecoords[1]) + ", " + 'ZF = ' + str(self.TS.finecoords[2]))
        
            #get position and slider values (might be needed in later versions)
            #self.updatePositionandFinePosition()
            
            # image center position
            XC = self.calX*self.imageData.shape[0]/2.0
            YC = self.calY*self.imageData.shape[1]/2.0
            
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
            
            #Project onto TEAM Stage axes
            rot = self.TS.STEMRot/180.0*np.pi #change to radians
            XMr = XM*np.cos(rot) + YM*np.sin(rot)
            YMr = -XM*np.sin(rot) + YM*np.cos(rot)
            
            print("XMr/YMr")
            print(XMr)
            print(YMr)
            
            #Determine X&ZY slider values (-100 to 100)
            newSliderX = self.TS.finecoords[0] + XMr / self.TS.sliderX_cal
            newSliderY = self.TS.finecoords[1] + YMr / self.TS.sliderY_cal
        
            #Test if in range else show error message and don't move (later move coarse?)
            newSliders = [newSliderX, newSliderY, self.TS.finecoords[2]]
            
            if ((newSliderX > -100) and (newSliderX < 100) and (newSliderY > -100) and (newSliderY < 100)):
                #Make a tuple of the new values
                self.setFinePosition(newSliders)
                self.sb.SetStatusText('New sliders: ' + str(newSliders[0]) + ", " + str(newSliders[1]) + ", " + str(newSliders[2]))
                self._lFinePos.SetLabel('Xf = ' + "{0:.2f}".format(newSliders[0]) + ", " + 'Yf = ' + "{0:.2f}".format(newSliders[1]) + ", " + 'Zf = ' +"{0:.2f}".format(newSliders[2]))
                #wait 1 second
                time.sleep(1)
                self.quickimage()
            else:
                steps = [1,2,4]
                steps[0] = int(XM / self.TS.step10X_cal)
                steps[1] = int(YM / self.TS.step10Y_cal)
                steps[2] = 0
                size = [1,1,1]
                self.sb.SetStatusText('Steps: ' + str(steps[0]) + ', ' + str(steps[1])+ ', ' + str(steps[2]))
                #self.step(steps, size)
                time.sleep(1)
                self.quickimage()
            
            pass
        pass
            
    def step(self, steps, size):
        self.TS.TS_Connect('localhost',5557)
        xsteps = steps[0]
        if xsteps < 0:
            self.TS.TS_Step([-xsteps,0,0],-size)
        elif xsteps>0:
            self.TS.TS_Step([xsteps,0,0],size)
        ysteps = steps[1]
        if ysteps < 0:
            self.TS.TS_Step([-ysteps,0,0],-size)
        elif ysteps>0:
            self.TS.TS_Step([ysteps,0,0],size)
        zsteps = steps[1]
        if zsteps < 0:
            self.TS.TS_Step([-zsteps,0,0],-size)
        elif zsteps>0:
            self.TS.TS_Step([zsteps,0,0],size)
        self.TS.TS_Disconnect()
        return 1
    
    def getPosition_compustage(self):
        ''' Get compustage position
        
        
        '''
        stageObj = self.Stage.Position
        posArray = (stageObj.X,stageObj.Y,stageObj.Z,stageObj.A,stageObj.B)
        return posArray
        
    def getPosition_teamstage(self):
        ''' Get the TEAM Stage position
        
        '''
        self.TS.TS_Connect('localhost',5557)
        self.TS.TS_GetPosition()
        self.TS.TS_Disconnect()
        return self.TS.coords
    
    def getFinePosition(self):   
        # Get the TEAM stage fine position coordinates
        self.TS.TS_Connect('localhost',5557)
        self.TS.TS_GetFinePosition()
        self.TS.TS_Disconnect()
        return self.TS.finecoords

    def setFinePosition(self, newSliders):  
        #Set the TEAM Stage fine coordinates
        self.TS.TS_Connect('localhost',5557)
        self.TS.TS_SetFinePosition(newSliders)
        self.TS.TS_Disconnect()
        return self.TS.finecoords
    
    def updatePositionSIM(self):
        # For testing offline
        print("using simulated Stage Positions")
        return (0,1,2,3,4)
        
    #def updatePosition(self):
    #    ''' Update the TEAM Stage position. Class holds the position.
    #
    #    '''
    #    
    #    self.TS.TS_Connect('localhost',5557)
    #    self.TS.TS_GetPosition()
    #    self.TS.TS_GetFinePosition()
    #    self.TS.TS_Disconnect()
    #    return 1
    
    def getPixelSize(self):
        # Get the pixel calibration from TIA
        window1 = self.TIA.ActiveDisplayWindow()
        Im1 = window1.FindDisplay(window1.DisplayNames[0]) #returns an image display object
        unit1 = Im1.SpatialUnit #returns SpatialUnit object
        self.unitName = unit1.unitstring #returns a string (such as nm)
        print('Image calibration = {}, {}'.format(Im1.image.calibration.deltaX,Im1.image.calibration.deltay))
        #self.calX = Im1.image.calibration.deltaX #returns the x calibration in meters
        #self.calY = Im1.image.calibration.deltaY
        print('Scan resolution = {}'.format(self.TIA.ScanningServer().ScanResolution))
        self.calX = self.TIA.ScanningServer().ScanResolution
        self.calY = self.TIA.ScanningServer().ScanResolution
        
    def quickimage(self):
        #Stop and ongiong acquisition
        self.stopAcqusition()
        
        #Acquire a quick image
        self.setStemSearchVals()
        acquiredImageSet = self.Acq.AcquireImages()

        with safearray_as_ndarray:
            self.imageData = acquiredImageSet(0).AsSafeArray # get data as ndarray
        
        self.getPixelSize() #get the pixel calibration
        imageShape = self.imageData.shape
        self.axes.clear()
        self.axes.imshow(np.fliplr(np.rot90(self.imageData,3)),extent=(0,self.calX*imageShape[0],0,self.calY*imageShape[1])) #add the image to the figure
        self.canvas.draw() #draw the panel
        self.sb.SetStatusText('Search image')
        
        self.imagevalid =1
    
    def determineMaxBinning(self):
        print('Determine maximum STEM binning')
        
        #Stop and ongiong acquisition
        self.stopAcqusition()
        
        #Acquire a quick image
        self.setStemSearchVals()
        acquiredImageSet = self.Acq.AcquireImages()

        with safearray_as_ndarray:
                imageData = acquiredImageSet(0).AsSafeArray # get data as ndarray
        
        sh = imageData.shape
        
        self.maxBin = sh[0] * self.myStemSearchParams.Binning
        
    def onShowPrev(self, event):
        self.showindex = self.showindex -1
        self.showimage()

    def onShowNext(self, event):
        self.showindex = self.showindex +1
        self.showimage()

    def showimage(self):
        #Show image from data set
        f = h5py.File(self.fullName + '.emd','a') #open and expand if file exists
        dataroot = f['data']
        dataTop = dataroot['raw']
        dset=dataTop['data']
        nindex = dataTop.attrs.get('nextindex')
        
        # check that in range
        if (self.showindex < dataTop.attrs.get('nextindex',0)) and (self.showindex>=0):
            self.showset = dataTop['imageParameters'][0,self.showindex] 
            self.showsubindex = dataTop['imageParameters'][1,self.showindex] 
            #self.imageData = dset[self.showsubindex,self.showset,:,:]
            self.imageData = dset[self.showset,self.showsubindex,:,:]
            self.calX = dataTop['dim4'][1]-dataTop['dim4'][0]
            self.calY = dataTop['dim3'][1]-dataTop['dim3'][0]
            imageShape = self.imageData.shape
            self.axes.clear()
            self.axes.imshow(np.fliplr(np.rot90(self.imageData,3)),extent=(0,self.calX*imageShape[0],0,self.calY*imageShape[1])) #add the image to the figure
            self.canvas.draw() #draw the panel
            self._lDispIndex.SetLabel('Image shown: ' + str(self.showindex))
            stage = f['stage']
            Position = dataTop['imageSetParameters'][0:5,self.showset]
            self._lDispFinePos.SetLabel('Fine position not recorded')
            self._lDispPos.SetLabel('X = ' + str(Position[0]) + ", " + 'Y = ' + str(Position[1]) + ", " + 'Z = ' + str(Position[2]))
            self._lDispAlphaGamma.SetLabel('Alpha = ' + str(Position[3]) + ", " + 'Gamma = ' + str(Position[4]))
            self.DispPosangles = Position[:]
            
            self.imagevalid =0
            # display related properties
        f.close()
        
    def stopAcqusition(self):
        if self.TIA.AcquisitionManager().isAcquiring:
            self.TIA.AcquisitionManager().Stop()
    
    def onAcquireFocal(self, event):

        self.getUsrAcqVals('focalseries')
        self.setStemAcqVals()
        
        # Setup focal series list
        initialDF = self.Proj.Defocus
        self.dfList = (np.arange(self.numDF)-(self.numDF-1)/2)*self.stepDF + initialDF
        print(self.dfList)
        self.initializeEMDFile(exType = 'focalseries')

        dataArray = np.zeros((self.numDF,self.numPerDF,self.expectedImageShape[0],self.expectedImageShape[1]),dtype='<u2')
        acqTimes = np.zeros((self.Rep,))
        
        for ii,df in enumerate(self.dfList):
                               
            # Set the defocus value
            self.Proj.Defocus = df
            print('Acquire at focus: {} nm'.format(np.round(df*1e9)))

            for jj in range(self.numPerDF):
                # Update the GUI
                self.sb.SetStatusText('Image #{} of {} at defocus {} nm'.format(ii+1, int(self.numDF), np.round(1e9*df)))
                self.acquireImages()
                dataArray[ii,jj,:,:] = self.imageData
            
            acqTimes[ii] = self.acqTime
        
        self.Proj.Defocus = initialDF
        self.imageData = dataArray
        self.times = acqTimes
    
        #Save the data to the file
        self.writeEMDdata('focalseries')
        
        del dataArray
        
        self.sb.SetStatusText('Finished focal series acquisition...')

    def onAcquireTimeSeries(self, event):
        ''' Acquire a STEM time series
        
        '''
        self.getUsrAcqVals('timeseries')
        self.setStemAcqVals()
        
        self.initializeEMDFile(exType = 'timeseries')
        
        dataArray = np.zeros((self.Rep,self.expectedImageShape[0],self.expectedImageShape[1]),dtype='<u2')
        acqTimes = np.zeros((self.Rep,))
        for ii in range(self.Rep):
            self.acquireImages()
            dataArray[ii,:,:] = self.imageData
            acqTimes[ii] = self.acqTime
            time.sleep(self.Del)
            
        self.imageData = dataArray
        self.times = acqTimes
    
        #Save the data to the file
        self.writeEMDdata('timeseries')
        
        del dataArray
        
        
        self.sb.SetStatusText('Finished time series acquisition...')
        
    def onAcquireDrift(self, event):
        ''' Acquire a set of stem rotation angles. Either
        0,90 or 0,90,45,135 degrees
        
        '''
        
        self.getUsrAcqVals('drift')
        self.setStemAcqVals()
        
        #Create a STEM rotation array
        rot0 = self.Ill.StemRotation # Current STEM rotation
        rotNum = self.rotSetting
        
        if rotNum == 0: #2 rotations is selected
            rotArray = np.array((rot0,rot0+np.pi/2.0)) #0, 90 pair
            self.Rep = 2
        elif rotNum == 1: #4 rotations is selected
            rotArray = np.array((rot0,rot0+np.pi/2.0,rot0+np.pi/4.0,rot0+3.0*np.pi/4.0)) # 0, 90 , 45, 135 degrees in that order
            self.Rep = 4
        else:
            self.Rep = 1
            rotArray = rot0*np.ones((1,))
        
        self.initializeEMDFile(exType = 'drift')
        
        dataArray = np.zeros((self.Rep,self.expectedImageShape[0],self.expectedImageShape[1]),dtype='<u2')
        acqTimes = np.zeros((self.Rep,))
        for ii,r in enumerate(rotArray):
            self.Ill.StemRotation = r
            self.acquireImages()
            acqTimes[ii] = self.acqTime
            dataArray[ii,:,:] = self.imageData
            
        self.imageData = dataArray
        self.times = acqTimes
        
        # Reset the stem rotation angle
        self.Ill.StemRotation = rot0
        
        #Save the data to the file
        self.writeEMDdata('drift')
        del dataArray
        
        self.sb.SetStatusText('Finished drift acquisition...')
    
    def onAcquireSingle(self, event):

        self.getUsrAcqVals('single')
        self.setStemAcqVals()
        self.Rep = 1 # only 1 image will be acquired
        
        self.initializeEMDFile(exType = 'single')
        
        #Acquire a single image
        self.acquireImages()
        self.times = time.time()
        #Save the data to the file
        self.writeEMDdata('single')
        
        self.sb.SetStatusText('Finished single acquisition...')
    
    def acquireImages(self):
        ''' Acquire a single image with the currently
        setup image acquisition settings.
        
        Will be called by onAcquireSingle, onAcquireFocal, etc.
        
        '''
        #print('pre pixel size = {} (m)'.format(self.calX))
        #Stop an ongoing acquisition
        self.stopAcqusition()
        
        #Update the TEAM Stage position
        #self.updatePosition()
        
        self.sb.SetStatusText('Acquiring image ...')
        
        # Get the acquisition time
        self.acqTime = time.time()

        self.Ill.BeamBlanked = False
        
        #Acquire the image
        acquiredImageSet = self.Acq.AcquireImages()

        with safearray_as_ndarray:
            self.imageData = acquiredImageSet(0).AsSafeArray.T # get data as ndarray oriented properly

        self.Ill.BeamBlanked = True
        
        # Get the current pixel size
        self.getPixelSize()
        
        # Show the image in the main window
        self.drawPlot(self.imageData) #add the image to the figure
        
        self.sb.SetStatusText('Finished Acquisition image ...')
        #print('post pixel size = {} (m)'.format(self.calX))

    def initializeEMDFile(self, exType):
        ''' Initialize the file and tags for HDF5 / EMD file format
        
        '''
        while os.path.exists(self.fullName + '_{:04}.emd'.format(self.II)):
            self.II += 1
        self.fullName += '_{:04}'.format(self.II)
        
        print('init emd file: {}'.format(self.fullName))
        print('experiment type: {}'.format(exType))
        
        with h5py.File(self.fullName + '.emd','w') as f:
            dataroot = f.create_group('data')
            dataroot.attrs['filename'] = str(self.Fprefix)
            dataTop = dataroot.create_group(exType)
            dataTop.attrs['emd_group_type'] = 1
            dataTop.attrs['uuid'] = str(uuid.uuid1())
            
            if exType == 'drift':
                dataTop.attrs['rotSetting'] = self.rotSetting
                if self.rotSetting == 0:
                    dataTop.attrs['rotations'] = [0,90]
                elif self.rotSetting == 1:
                    dataTop.attrs['rotations'] = [0,90,45,135]
                else:
                    dataTop.attrs['rotations'] = 'none'
            elif exType == 'timeseries':
                dataTop.attrs['number images'] = self.Rep
                dataTop.attrs['delay time'] = self.Del
            elif exType == 'focalseries':
                dataTop.attrs['number images'] = self.Rep
                dataTop.attrs['number defocus'] = self.numDF
                dataTop.attrs['number per df'] = self.numPerDF
            
            pix = self.maxBin/self.Bin
            print('num pixesl = {}'.format(pix))
            
            # Create dataset with correct shape
            if exType == 'single':
                sh = (pix,pix)
                maxsh = (pix,pix)
                chunks = (pix,pix)
            elif exType == 'focalseries':
                sh = (self.numDF,self.numPerDF, pix, pix)
                maxsh = (self.numDF,self.numPerDF, pix, pix)
                chunks = (1,1, pix, pix)
            else:
                sh = (self.Rep, pix, pix)
                maxsh = (self.Rep, pix, pix)
                chunks = (1, pix, pix)

            # Initialize the data set
            dset = dataTop.create_dataset('data',sh,'<u2',maxshape=maxsh,chunks=chunks)
            #print(sh)
            # Create the EMD dimension datasets
            _ = self.createDims(dataTop, exType, pix)
             
            microscope = f.create_group('microscope')
            microscope.attrs['name'] = self.microscopeName
            microscope.attrs['high tension'] = self._microscope.Gun.HTValue
            microscope.attrs['spot size'] = self.Ill.SpotsizeIndex
            microscope.attrs['magnification'] = self.Ill.StemMagnification
            microscope.attrs['defocus'] = self.Ill.ProbeDefocus
            microscope.attrs['convergence angle'] = self.Ill.ConvergenceAngle
            microscope.attrs['camera length'] = self.Proj.CameraLength
            microscope.attrs['binning'] = self.Bin
            microscope.attrs['max binning'] = self.maxBin
            microscope.attrs['dwell time'] = self.Dwell
            
            stage = f.create_group('stage')
            stage.attrs['type'] = self.stageName
            stage.attrs['position'] = self.getPosition()

            user = f.create_group('user')
            user.attrs['user name'] = self._iUser.GetValue()
            
            sample = f.create_group('sample')
            sample.attrs['sample name'] = self._iSample.GetValue()
            
        #Close the file
    
    def createDims(self, dataTop, exType, pix):
        
        if exType == 'single':
            dim2 = dataTop.create_dataset('dim2',(pix,),'f')
            dim2.attrs['name'] = np.string_('X')
            dim2.attrs['units'] = np.string_('n_m')
            dim1 = dataTop.create_dataset('dim1',(pix,),'f')
            dim1.attrs['name'] = np.string_('Y')
            dim1.attrs['units'] = np.string_('n_m')
            #dims = (dim1,dim2)
        elif exType == 'focalseries':
            dim4 = dataTop.create_dataset('dim4',(pix,),'f')
            dim4.attrs['name'] = np.string_('X')
            dim4.attrs['units'] = np.string_('n_m')
            dim3 = dataTop.create_dataset('dim3',(pix,),'f')
            dim3.attrs['name'] = np.string_('Y')
            dim3.attrs['units'] = np.string_('n_m')
        else:
            dim3 = dataTop.create_dataset('dim3',(pix,),'f')
            dim3.attrs['name'] = np.string_('X')
            dim3.attrs['units'] = np.string_('n_m')
            dim2 = dataTop.create_dataset('dim2',(pix,),'f')
            dim2.attrs['name'] = np.string_('Y')
            dim2.attrs['units'] = np.string_('n_m')
            
        if exType == 'drift':
            dim1 = dataTop.create_dataset('dim1',(self.Rep,),'f')
            dim1.attrs['name'] = np.string_('STEM rotation')
            dim1.attrs['units'] = np.string_('deg')
        elif exType == 'timeseries':
            dim1 = dataTop.create_dataset('dim1',(self.Rep,),'f')
            dim1.attrs['name'] = np.string_('time')
            dim1.attrs['units'] = np.string_('sec')
        elif exType == 'focalseries':
            dim1 = dataTop.create_dataset('dim1',(self.numDF,),'f')
            dim1.attrs['name'] = np.string_('defocus')
            dim1.attrs['units'] = np.string_('n_m')
            dim2 = dataTop.create_dataset('dim2',(self.numPerDF,),'f')
            dim2.attrs['name'] = np.string_('')
            dim2.attrs['units'] = np.string_('')

        return 1
    
    def writeEMDdata(self, exType):
        with h5py.File(self.fullName + '.emd','a') as f:
            dataroot = f['data']
            dataTop = dataroot[exType]
            dims = [dataTop['dim1'],dataTop['dim2']]
            try:
                # Get third dim if it exists
                dims.append(dataTop['dim3'])
            except:
                pass
            try:
                # Get fourth dim if it exists
                dims.append(dataTop['dim4'])
            except:
                pass
                
            dset = dataTop['data']
            stage = f['stage']
            user = f['user']
            sample = f['sample']
            
            imageShape = self.imageData.shape[-2:]
            xdim = np.linspace(0,(imageShape[0]-1)*self.calX*1e9,imageShape[0]) #multiply by 1e9 for nanometers
            ydim = np.linspace(0,(imageShape[1]-1)*self.calY*1e9,imageShape[1])
            dims[-1][:] = xdim
            dims[-2][:] = ydim
            dims[-2][:] = ydim
            
            if len(dims) > 2:
                if exType == 'drift':
                    if self.rotSetting == 0:
                        dims[0][:] = (0,90)
                    elif self.rotSetting == 1:
                        dims[0][:] = (0,90,45,135)
                elif exType == 'timeseries':
                    dims[0][:] = self.times - self.times[0]
                elif exType == 'focalseries':
                    dims[0][:] = self.dfList
                    dims[1][:] = range(self.numPerDF)

            dataTop.create_dataset('acquisition times', data = self.times)
            
            self.sb.SetStatusText('Writing image(s)...')
            dset[:] = self.imageData
            
            self.imageData = 0 # set the image data to 0 to free memory
# Parse the arguments
parser = argparse.ArgumentParser(description='Acquire sets of STEM images.')

# Allow user to set TeamStage or Compustage
parser.add_argument('--teamstage', '-t', action='store_const', const=True, default=False)
parser.add_argument('--compustage', '-c', action='store_const', const=True, default=True)

parser.add_argument('-v','--version', action='version', version='%(prog)s {}'.format(version))

args = parser.parse_args()

if args.teamstage:
    stagetype = 'teamstage'
else:
    stagetype = 'compustage'

print('Stage type = {}'.format(stagetype))
 
app = wx.App(False)
frame = TEAMFrame(None, "NCEM: STEM Experiments", stagetype = stagetype)
app.MainLoop()
