""" Python script for determining the tilt required
to bring a crystal on axis by clicking on a position
in a CBED pattern.

Currently, only works for the NCEM TEAM Stage and the
FEI Flucam. The flucam data is gotten from TIA after acquisition.

300 kV: similar to 200 kV
200kV: alpha orientation = -18, polarity = -1

author: Peter Ercius, percius@lbl.gov
"""

import argparse

import numpy as np
import socket
import os
import wx  # wxPython GUI package
import h5py
import time

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.figure import Figure
from matplotlib import patches

# For connections to FEI TEMScripting and TIA
from comtypes.client import CreateObject
from comtypes.safearray import safearray_as_ndarray

import quaternions

# For connecting to TEAM Stage
try:
    import TEAMstageclass
except ModuleNotFoundError:
    print('No TEAM Stage functions available.')

version = 0.1

# Main window
class TilterFrame(wx.Frame):
    def __init__(self, parent, title, stagetype='compustage', verbose=False):

        try:
            # This does not exist in some versions
            from matplotlib.pyplot import set_cmap
        except:
            pass
        
        self.v = verbose
        
        # Setup the stagetype as compustage or teamstage
        self.stagetype = stagetype
        
        self.stage_alpha = 0
        self.stage_gamma = 0
        self.stage_beta = 0
        
        # Where to save the diffraction patterns and clicks
        self.filename = 'g:/UserData/Ercius/flucam_data/KStilter.emd'
        print('Saving data to: {}'.format(self.filename))
        
        self.cbedXY = (255, 255)
        self.axisXY = (0, 0)
        self.cbed = None
        self.pix = (512, 512) # size of a flucam image
        self.alpha = 0
        self.gamma = 0
        self.new_alpha = 0
        self.new_gamma = 0
        
        # TODO Add in configuration option for cbed pixel sizes and change
        # to radians (see calc_tilts also.)
        self.cbed_pixel_size = 0.041  # in degrees
        
        # Todo: Add configuration option for alpha axis orientation
        self.alpha_orientation = 18
        self.alpha_polarity = -1
        print('Alpha polarity = {} degrees for 200 and 300 kV'.format(self.alpha_polarity))
        print('Alpha axis orientation = {} degrees for 200 and 300 kV'.format(self.alpha_orientation))

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
        
        # Initialize the base Frame
        wx.Frame.__init__(self, parent, title=title, size=(1000, 700))
        
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
        #displayvalbox = wx.BoxSizer(wx.VERTICAL)
        #imagebuttonbox = wx.BoxSizer(wx.HORIZONTAL)
        currentvalbox = wx.BoxSizer(wx.VERTICAL)
        setupbox = wx.BoxSizer(wx.VERTICAL)
        tiltbox = wx.BoxSizer(wx.VERTICAL)
        alignbox = wx.BoxSizer(wx.VERTICAL)

        notebook = wx.Notebook(self)
        
        # create the page windows as children of the notebook
        setupPage = wx.Panel(notebook)
        tiltPage = wx.Panel(notebook)
        alignPage = wx.Panel(notebook)
        
        # Create a figure in a canvas
        self.TEMfigure = Figure()
        self.canvas = FigureCanvasWxAgg(self, -1, self.TEMfigure)
        self.ax0 = self.TEMfigure.add_subplot(111)
        self.ax0.xaxis.set_visible(False)
        self.ax0.yaxis.set_visible(False)
        try:
            set_cmap('viridis')
        except NameError:
            pass
        except ValueError:
            pass
            
        self.ax0im = self.ax0.imshow(np.zeros(self.pix,dtype='<u2'))
        self.TEMfigure.tight_layout()
        
        # Add central cbed circle
        self.zero_beam = patches.Circle(self.cbedXY,
                                        10, fill=False, lw=2, edgecolor='r')
        self.ax0.add_patch(self.zero_beam)
        
        # Add a circle to show up when clicked
        self.axis_center = patches.Circle((128, 127),
                                          10, fill=False, lw=2, edgecolor='w')
        self.ax0.add_patch(self.axis_center)
        
        # Setup page        
        self._lUser = wx.StaticText(setupPage, wx.ID_ANY, 'User name')
        self._iUser = wx.TextCtrl(setupPage,wx.ID_ANY, value='', 
                                  size=(300,20))
        self._lSample = wx.StaticText(setupPage, wx.ID_ANY, 'Sample name')
        self._iSample = wx.TextCtrl(setupPage, wx.ID_ANY, value='', 
                                    size=(300,20))
        self._label_alpha_orientation = wx.StaticText(setupPage, wx.ID_ANY, 'Alpha orientation (deg)')
        self._val_alpha_orientation = wx.TextCtrl(setupPage, wx.ID_ANY, value='18', 
                                    size=(300,20))
        self._label_alpha_polarity = wx.StaticText(setupPage, wx.ID_ANY, 'Alpha polarity')
        self._val_alpha_polarity = wx.TextCtrl(setupPage, wx.ID_ANY, value='-1', 
                                    size=(300,20))
        self._label_cbed_pixel_size = wx.StaticText(setupPage, wx.ID_ANY, 'CBED pixel size (deg)')
        self._val_cbed_pixel_size = wx.TextCtrl(setupPage, wx.ID_ANY, value='0.041', 
                                    size=(300,20))
        self._lStage = wx.StaticText(setupPage, wx.ID_ANY, 'Stage type')
        self._iStage = wx.RadioBox(setupPage, wx.ID_ANY,
                                   label='stage type', choices=('alpha-gamma',
                                                                'alpha-beta'))

        setupbox.Add(self._lUser)
        setupbox.Add(self._iUser)
        setupbox.Add(self._lSample)
        setupbox.Add(self._iSample)
        setupbox.Add(self._lStage)
        setupbox.Add(self._iStage)
        setupbox.Add(self._label_alpha_orientation)
        setupbox.Add(self._val_alpha_orientation)
        setupbox.Add(self._label_alpha_polarity)
        setupbox.Add(self._val_alpha_polarity)
        setupbox.Add(self._label_cbed_pixel_size)
        setupbox.Add(self._val_cbed_pixel_size)
        
        # Tilt page
        self._label_alpha = wx.StaticText(tiltPage, wx.ID_ANY, 'Alpha')
        self._val_alpha = wx.StaticText(tiltPage, wx.ID_ANY, '0')
        self._label_gamma = wx.StaticText(tiltPage, wx.ID_ANY, 'Gamma')
        self._val_gamma = wx.StaticText(tiltPage, wx.ID_ANY, '0')
        self._btn_get_cbed = wx.Button(tiltPage, label='Get CBED')
        self._label_new_alpha = wx.StaticText(tiltPage, wx.ID_ANY, 
                                              'Target Alpha')
        self._val_new_alpha = wx.StaticText(tiltPage, wx.ID_ANY, '0')
        self._label_new_gamma = wx.StaticText(tiltPage, wx.ID_ANY, 
                                              'Target Gamma')
        self._val_new_gamma = wx.StaticText(tiltPage, wx.ID_ANY, '0')
        self._btn_accept_tilt = wx.Button(tiltPage, label='Accept tilt')
        
        tiltbox.Add(self._label_alpha)
        tiltbox.Add(self._val_alpha)
        tiltbox.Add(self._label_gamma)
        tiltbox.Add(self._val_gamma)
        tiltbox.Add(self._btn_get_cbed)
        tiltbox.Add(self._label_new_alpha)
        tiltbox.Add(self._val_new_alpha)
        tiltbox.Add(self._label_new_gamma)
        tiltbox.Add(self._val_new_gamma)
        tiltbox.Add(self._btn_accept_tilt)

        # Align page
        self._label_cbed_center = wx.StaticText(alignPage, wx.ID_ANY, 'CBED center (pixels)')
        self._val_cbed_center = wx.StaticText(alignPage, wx.ID_ANY, str(self.cbedXY))
        self._label_cbed_pixel_size = wx.StaticText(alignPage, wx.ID_ANY, 'CBED pixel size (deg)')
        self._val_align_cbed_pixel_size = wx.StaticText(alignPage, wx.ID_ANY, str(self.cbed_pixel_size))
        self._btn_get_cbed_align = wx.Button(alignPage, label='Get CBED')
        # allow manual convergence angle
        self._label_convergence = wx.StaticText(alignPage, wx.ID_ANY, 
                                              'Convergence angle (mrad)')
        self._val_convergence = wx.TextCtrl(alignPage, wx.ID_ANY, value='30.0', 
                                            size=(300,20))
        self._btn_auto_cbed_align = wx.Button(alignPage, label='Auto align')
        
        alignbox.Add(self._label_cbed_center)
        alignbox.Add(self._val_cbed_center)
        alignbox.Add(self._label_cbed_pixel_size)
        alignbox.Add(self._val_align_cbed_pixel_size)
        alignbox.Add(self._label_convergence)
        alignbox.Add(self._val_convergence)
        alignbox.Add(self._btn_get_cbed_align)
        alignbox.Add(self._btn_auto_cbed_align)

        # Bind the events
        self._btn_get_cbed.Bind(wx.EVT_BUTTON, self.on_get_cbed)
        self._btn_accept_tilt.Bind(wx.EVT_BUTTON, self.on_accept_tilt)
        self._btn_get_cbed_align.Bind(wx.EVT_BUTTON, self.on_get_cbed)
        self._btn_auto_cbed_align.Bind(wx.EVT_BUTTON, self.on_auto_align)

        # Enable Mouse events
        self.canvas.mpl_connect('button_press_event', self.on_click) 

        imagebox.Add(self.canvas, proportion=3, 
                     flag=wx.LEFT|wx.RIGHT|wx.EXPAND, border=0)
                
        imagecolumnbox.Add(imagebox, proportion=1, 
                           flag=wx.LEFT|wx.RIGHT|wx.EXPAND, border=0)
        
        infobox.Add(currentvalbox, proportion=0, 
                    flag=wx.LEFT|wx.RIGHT|wx.EXPAND, border=0)
        
        setupPage.SetSizer(setupbox)
        tiltPage.SetSizer(tiltbox)
        alignPage.SetSizer(alignbox)
        #focalPage.SetSizer(focalbox)
        #driftPage.SetSizer(driftbox)
        #timeSeriesPage.SetSizer(seriesbox)
        #reviewpage.SetSizer(displayvalbox)
        
        # add the pages to the notebook with the label to show on the tab
        notebook.AddPage(setupPage, 'Setup')
        notebook.AddPage(tiltPage, 'Tilt')
        notebook.AddPage(alignPage, 'Align')
        
        mainbox.Add(imagecolumnbox, proportion=1, 
                    flag=wx.LEFT|wx.RIGHT|wx.EXPAND, border=0)
        mainbox.Add(notebook, proportion=0, 
                    flag=wx.LEFT|wx.RIGHT|wx.EXPAND, border=20)
        
        self.SetSizer(mainbox)  # Initialize the parent sizer
        
        # Set color
        self.SetBackgroundColour("0000FF")
        
        self.sb = self.CreateStatusBar()
        self.sb.SetStatusText('Idle...')
        
        self.on_connect() # Connect to the microscope and stage
        
        self.Show(True) # show the GUI
    
    def __del__(self):
        self.f.close()
        self.TS.TS_Disconnect()
        
    def drawPlot(self, imArray):
        #print('Draw plot')
        
        # Need to add updates to the color limits. Not sure how to do that.
        # self.ax0Im.set_data(imArray) <-- This should be faster.
        # but it does not work. Need to update 

        self.ax0.clear()
        if imArray.shape[0] <= 2048:
            axIm = self.ax0.imshow(imArray)
        else:
            # 4kx4k cant be shown
            axIm = self.ax0.imshow(imArray[::2,::2])
        
        # Add the zero beam marker
        self.ax0.add_patch(self.zero_beam)
        self.canvas.draw()
    
    def on_connect(self):
        """ Connect to the microscope and the stage. 
        Setup the detectors
        
        """
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
            #  Connect to TEAM stage
            try:
                print('Trying to connect to TEAM Stage')
                self.TS = TEAMstageclass.TEAMstage()
                self.TS.TS_Connect('localhost', 5557)
            except:
                print('Error with TEAM Stage connection.')
                raise
            
        try:
            # Get microscope interfaces
            self.Acq = self._microscope.Acquisition
            self.Proj = self._microscope.Projection
            self.Ill = self._microscope.Illumination
            self.Stage = self._microscope.Stage # FEI stage
            
            # Older pythoncom versions might require square brackets [] below
            self.detector0 = self.Acq.Detectors(0)
        except:
            print('Connections to microscope interfaces failed.')
            raise

        # Add the first detector
        self._microscope.Acquisition.AddAcqDevice(self.detector0)
        
        # Initialize the EMD file
        print('init emd')
        self.initializeEMDFile()

    def setFonts(self, ed):
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
    
    def on_get_cbed(self, event):
        """ Gets the image in the currently active display window.
        Also, gest the cbed and alpha axis input by user
        
        """
        self.axis_center.center = (None, None) # remove circle from plot
        #self.TEMfigure.canvas.draw() # probably not needed
        
        # Get user input for alpha orientation and cbed pixel size
        # Todo: get this from a config file instead
        self.alpha_orientation = float(self._val_alpha_orientation.GetValue())
        self.cbed_pixel_size = float(self._val_cbed_pixel_size.GetValue())
        self.alpha_polarity = float(self._val_alpha_polarity.GetValue())
        pos = self.getPosition()
        self._val_alpha.SetLabel('{:0.3}'.format(pos[3]))
        self._val_gamma.SetLabel('{:0.3}'.format(pos[4]))
        self._val_new_alpha.SetLabel('{:0.3}'.format(pos[3]))
        self._val_new_gamma.SetLabel('{:0.3}'.format(pos[4]))
        
        # Show convergence angle
        self._val_convergence.SetLabel('{:0.3}'.format(self.Ill.ConvergenceAngle*1e3))
        
        window1 = self.TIA.ActiveDisplayWindow()
        disp1 = window1.FindDisplay(window1.DisplayNames[0])
        
        self.cbed = np.rot90(np.asarray(disp1.Image.data.array),1)
        self.drawPlot(self.cbed)
    
    def on_accept_tilt(self, event):
        # Tell TS to rotate
        # save to EMD
        print('Rotate to: {}, {} degrees'.format(self.new_alpha * 180 / np.pi,
                                                 self.new_gamma * 180 / np.pi))
        self.writeEMDdata()

    def on_auto_align(self, event):
        """ Determine center and radius of vacuum probe
        Threshold cbed pattern
        Use area to estimate radius
        Find center of mass

        """
        im = self.cbed
        imr = im.ravel()
        
        th = np.mean(imr)
        for ii in range(100):
            # an auotomated threshold
            foreg = imr >= th
            backg = imr < th

            thf = np.mean(imr[foreg])
            thg = np.mean(imr[backg])

            th2 = (thf + thg) / 2.0

            err = np.abs(th - th2)
            
            if err < 0.1:
                break
            else:
                th = th2

        im_th = im.copy()
        im_th[im >= th2] = 1
        im_th[im < th2] = 0

        # Find the center of mass
        mn = im_th.sum()

        XX, YY = np.mgrid[0:im.shape[0], 0:im.shape[1]]
        cx = np.sum(XX * im_th)/mn
        cy = np.sum(YY * im_th)/mn
        if self.v:
            print('cx,cy = {}, {}'.format(cx, cy))
        self.cbedXY = (cy, cx)

        # Estimate the radius using area
        # Convert to degrees for now
        # Todo: Change cbed pixel size to radians
        radius0 = np.sqrt(im_th.sum() / np.pi) # pixels
        
        # Get the convergence angle
        #a = self.Ill.ConvergenceAngle
        a = float(self._val_convergence.GetLabel())
        a = a * 1.0e-3 # convert to radians
        
        self.cbed_pixel_size = (a / radius0) * 180 / np.pi
        if self.v:
            print('New cbed pixel size = {} deg'.format(self.cbed_pixel_size))
            print('CBED radius = {} pixels'.format(radius0))
            print('CBED radius = {} rad'.format(radius0*self.cbed_pixel_size*np.pi/180.))
        # Update the GUI
        self._val_cbed_pixel_size.SetLabel(str(self.cbed_pixel_size))
        self._val_align_cbed_pixel_size.SetLabel(str(self.cbed_pixel_size))
        self._val_cbed_center.SetLabel(str(self.cbedXY))
        
        # Update the zero beam circle center and radius
        self.zero_beam.center = self.cbedXY
        self.zero_beam.set_radius(radius0)
        self.TEMfigure.canvas.draw()
        
        if self.v:
            print('self.cbedXY = {}'.format(self.cbedXY))
            print('self.zero_beam_center = {}'.format(self.zero_beam.center))
        
    def on_click(self, event):
        if self.v:
            print('{} click: button={}, x={}, y={}, xdata={}, ydata={}'.format(
                 'double' if event.dblclick else 'single', event.button,
                  event.x, event.y, event.xdata, event.ydata))
        if event.button == 1 and event.xdata and event.ydata:
            # left click calculates new angles
            self.new_alpha, self.new_gamma = self.calc_tilts(event.xdata, 
                                                             event.ydata)
            if self.v:
                print('New tilts = {0:0.3}, {1:0.3}'.format(self.new_alpha*180/np.pi, 
                                                            self.new_gamma*180/np.pi))
            # Update the position of the axis
            self.axisXY = (event.xdata, event.ydata)
            
            self._val_new_alpha.SetLabel('{:0.3}'.format(self.new_alpha * 180/np.pi))
            self._val_new_gamma.SetLabel('{:0.3}'.format(self.new_gamma * 180/np.pi))
            
            self.axis_center.center = (event.xdata, event.ydata)
            self.TEMfigure.canvas.draw()
            return
        elif event.button == 3 and event.xdata and event.ydata:
            # right click sets cbed center
            self.zero_beam.center = (event.xdata, event.ydata)
            
            self.TEMfigure.canvas.draw()
            self.cbedXY = [event.xdata, event.ydata] # update the cbed center
            if self.v:
                print('new cbedXY = {}'.format(self.cbedXY))
            
            self._val_align_cbed_pixel_size.SetLabel(str(self.cbed_pixel_size))
            self._val_cbed_center.SetLabel(str(self.cbedXY))
            self.zero_beam.center = self.cbedXY
            return
            
    
    def calc_tilts(self, x, y):
        """ Pass through function which is based on alpha-gamma or 
        alpha-beta radio button.
        
        Parameters
        ----------
        x, y : float, pixels
            The location clicked on the diffraction pattern image.
        
        """
        print(self._iStage.GetSelection())
        if self._iStage.GetSelection() == 0:
            out = self.calc_tilts_ag(x, y)
        else:
            out = self.calc_tilts_ab(x, y)
        return out
        
        
    def calc_tilts_ag(self, x, y):
        """ Calculate tilts for alpha gamma TEAM stage. Values are
        returned are the new stage angles to rotate to in
        radians.
        
        Parameters
        ----------
        x, y : float, pixels
            The location clicked on the diffraction pattern image
            
        Returns
        -------
        : 2-tuple
            The alpha and gamm tilt to bring the axis to the 
            beam orientation. In radians.
        """
        alpha_orientation = self.alpha_orientation * np.pi / 180.
        
        pos = self.getPosition()
        tsAlpha = pos[3] * np.pi / 180
        tsGamma = pos[4] * np.pi / 180
        beam = [0, 0, 1] # beam points along z axis in world-beam frame
        
        # both in pixels
        axisXY = (x, y) # where clicked
        cbedXY = self.cbedXY # where the zero beam is
        
        if self.v:
            print('cbedXY = {}'.format(cbedXY))
            print('axisXY = {}'.format(axisXY))
        
        # Find x and y distance on the plot
        deltaXY = [x1 - x2 for (x1, x2) in zip(axisXY, cbedXY)]
        
        # Reverse X to get gamma direction correct
        deltaXY[0] *= -1
        
        # Rotate the delta vector to accommodate the alpha axis rotation
        rot = np.array(((np.cos(alpha_orientation), -np.sin(alpha_orientation)), 
                        (np.sin(alpha_orientation), np.cos(alpha_orientation))))
        deltaXY = np.dot(deltaXY, rot)
        
        if self.v:
            print("deltaXY = {0[0]:0.3}, {0[1]:0.3}".format(deltaXY))
        
        # Determine if positive or negative alpha
        if deltaXY[1] <= 0:
            alpha_sign = self.alpha_polarity
        else:
            alpha_sign = -self.alpha_polarity
        if self.v:
            print('alpha sign = {}'.format(alpha_sign))
        
        # Measure theta and phi using polar coordinates
        theta = alpha_sign * np.sqrt(deltaXY[0]**2 + deltaXY[1]**2) * self.cbed_pixel_size * np.pi/180 #theta can be calculated directly from the user click point and diff pattern calibration
        phi = np.arctan(deltaXY[0]/(deltaXY[1] + 1e-5)) # use arctan2 instead?
        phi_arctan2 = np.arctan2(deltaXY[0],(deltaXY[1] + 1e-5)) 
        
        if self.v:
            print("theta = {} rad".format(theta))
            print("theta = {} deg".format(theta * 180/np.pi))
            print("phi = {} rad, {} deg".format(phi, phi*180/np.pi))
            print("phi2 = {} rad, {} deg".format(phi_arctan2, phi_arctan2*180/np.pi))

        # Find the quaternion of the stage rotation
        # this rotates the beam [0, 0, 1] into stage coordinates
        qA = quaternions.fromAngleAxis(tsAlpha, [-1, 0, 0]) #alpha axis is X-axis
        qG = quaternions.fromAngleAxis(tsGamma, [0, 0, 1])
        world2tsQ = quaternions.normalize(quaternions.multiply(qA, qG)) #first gamma then A; gamma is attached to Alpha stage

        # Rotate the beam to the axis in world coordinates
        axis_qA = quaternions.fromAngleAxis(theta, [-1, 0, 0]) #alpha axis is X-axis
        axis_qG = quaternions.fromAngleAxis(phi, [0, 0, 1])
        beam2axisQ = quaternions.normalize(quaternions.multiply(axis_qA, axis_qG)) #first gamma then A; gamma is attached to Alpha stage
        axisV = quaternions.applyToVector(beam2axisQ, beam)
        
        # Find the axis location in the tilted TS coordinates
        axisV_TSframe = quaternions.applyToVector(world2tsQ, axisV)
        
        if self.v:
            print("axis in world frame = {}".format(axisV))
            print("axis in TS frame = {}".format(axisV_TSframe))
        
        # The alpha and gamma values to tilt to (radians)
        new_alpha_quaternion = np.arccos(axisV_TSframe[2]) # this loses the sign
        new_alpha = tsAlpha + theta # directly from the calculated theta value
        new_gamma = np.arctan(axisV_TSframe[0]/(axisV_TSframe[1] + 1e-5))
        print('TRY ARCTAN2 HERE for TEAM STAGE')
        
        if self.v:
            print("axis in world frame = {}".format(axisV))
            print("new alpha quaternion = {:0.3}".format(new_alpha_quaternion * 180/np.pi))
            print("new gamma quaternion = {:0.3}".format(new_gamma * 180/np.pi))
            print("new alpha from theta = {:0.3}".format((tsAlpha + theta)*180/np.pi))
        return new_alpha, new_gamma
    
    def calc_tilts_ab(self, x, y):
        """ Calculate tilts for alpha beta double tilt holder. Values are
        returned are the new stage angles to rotate to in
        radians.
        
            Parameters
            ----------
            x, y : float, pixels
                The location clicked on the image diffraction pattern
            
            cbedXY : 2-tuple
                The location of the center of the pattern to calculate
                the relative tilts
                
            verbose : bool, default False
                Print out information about the calculation
                
            Returns
            -------
            : 2-tuple
                The alpha and beta tilt to bring the axis to the 
                beam orientation. In radians.
        """
        alpha_orientation = self.alpha_orientation * np.pi / 180. # beta orientation is assumed to be orthogonal
        
        pos = self.getPosition()
        stage_alpha = pos[3] * np.pi / 180
        stage_gamma = pos[4] * np.pi / 180
        
        # both in pixels
        axisXY = (x, y) # where clicked
        cbedXY = self.cbedXY # where the zero beam is
        
        if self.v:
            print('cbedXY = {}'.format(cbedXY))
            print('axisXY = {}'.format(axisXY))
        
        # Find x and y distance on the plot
        deltaXY = [x1 - x2 for (x1, x2) in zip(axisXY, cbedXY)]
        
        # Reverse X to get beta direction correct
        print('Check beta direction reversal')
        deltaXY[0] *= -1
        
        # Rotate the delta vector to accommodate the alpha axis rotation
        rot = np.array(((np.cos(alpha_orientation), -np.sin(alpha_orientation)), 
                        (np.sin(alpha_orientation), np.cos(alpha_orientation))))
        deltaXY = np.dot(deltaXY, rot)
        
        if self.v:
            print("deltaXY = {0[0]:0.3}, {0[1]:0.3}".format(deltaXY))
        
        # Determine if positive or negative alpha
        if deltaXY[1] <= 0:
            alpha_sign = self.alpha_polarity
        else:
            alpha_sign = -self.alpha_polarity
        if self.v:
            print('alpha sign = {}'.format(alpha_sign))
        
        # Determine positive or negative beta
        self.beta_polarity = 1
        print('beta_polarity set to +1')
        if deltaXY[1] <= 0:
            beta_sign = self.beta_polarity
        else:
            beta_sign = -self.beta_polarity
        if self.v:
            print('alpha sign = {}'.format(alpha_sign))
        
        alpha = alpha_sign * deltaXY[0] * self.cbed_pixel_size * np.pi/180
        beta  = beta_sign * deltaXY[1] * self.cbed_pixel_size * np.pi/180
        
        new_alpha = stage_alpha + alpha
        new_beta = stage_alpha + beta
        
        if self.v:
            print('alpha, beta = {}, {}'.format(alpha, beta))
            print('new_alpha, new_beta = {}, {}'.format(new_alpha, new_beta))
        return new_alpha, new_beta
        
    def getPosition_compustage(self):
        """ Get compustage position
        
        
        """
        stageObj = self.Stage.Position
        posArray = (stageObj.X, stageObj.Y, stageObj.Z, stageObj.A, stageObj.B)
        return posArray
        
    def getPosition_teamstage(self):
        """ Get the TEAM Stage position
        
        """
        self.TS.TS_Connect('localhost',5557)
        self.TS.TS_GetPosition()
        self.TS.TS_Disconnect()
        return self.TS.coords
    
    def getFinePosition(self):   
        """ Get the TEAM stage fine position coordinates
        
        """
        self.TS.TS_Connect('localhost',5557)
        self.TS.TS_GetFinePosition()
        self.TS.TS_Disconnect()
        return self.TS.finecoords

    def setFinePosition(self, newSliders):  
        """ Set the TEAM Stage fine coordinates
        """
        self.TS.TS_Connect('localhost',5557)
        self.TS.TS_SetFinePosition(newSliders)
        self.TS.TS_Disconnect()
        return self.TS.finecoords
    
    def updatePositionSIM(self):
        # For testing offline
        print("using simulated Stage Positions")
        return (0, 1, 2, 3, 4)
             
    def initializeEMDFile(self):
        """ Initialize the file and tags for HDF5 / EMD file format only
        if it does not already exist.
        
        """
        
        if not os.path.exists(self.filename):
            print('init emd file: {}'.format(self.filename))
        
            with h5py.File(self.filename, 'a') as f:
                data_root = f.create_group('data')
                data_root.attrs['filename'] = str(self.filename)
                g_cbed = data_root.create_group('cbed')
                g_cbed.attrs['emd_group_type'] = 1
                
                # Initialize the data set
                dset = g_cbed.create_dataset('data', dtype='<u2', shape=(1, self.pix[0], self.pix[1]),
                                             maxshape=(None, self.pix[0], self.pix[1]),
                                             chunks=(1, self.pix[0], self.pix[1]),compression='lzf')
                # Create the EMD dimension datasets
                _ = self.createDims(g_cbed, self.pix)
                
                # click_x, click_y, cbed_x, cbed_y, HT, CL, cbed_pixel_size, 
                # convergence angle, time
                dset_params = g_cbed.create_dataset('params', dtype='f', shape=(1, 9),
                                                    maxshape=(None, 9))
                dset_params.attrs['columns'] = ('click_x', 'click_y', 
                                                'cbed_x', 'cbed_y',
                                                'HT', 'camera_length',
                                                'cbed_pixel_size',
                                                'convergence_angle',
                                                'time')
                
                # User name, sample name
                self.names_dt = h5py.special_dtype(vlen=str)
                dset_names = g_cbed.create_dataset('names', shape=(1, 2),
                                                   dtype=self.names_dt,
                                                   maxshape=(None, 2))
                dset_names.attrs['columns'] = ('user_name', 'sample_name')
                
                # Sample positions
                dset_position = g_cbed.create_dataset('position', dtype='f', shape=(1, 5),
                                                      maxshape=(None, 5))
                
                microscope = f.create_group('microscope')
                microscope.attrs['name'] = self.microscopeName
    #            microscope.attrs['high tension'] = self._microscope.Gun.HTValue
    #            microscope.attrs['spot size'] = self.Ill.SpotsizeIndex
    #            microscope.attrs['magnification'] = self.Ill.StemMagnification
    #            microscope.attrs['defocus'] = self.Ill.ProbeDefocus
    #            microscope.attrs['convergence angle'] = self.Ill.ConvergenceAngle
    #            microscope.attrs['camera length'] = self.Proj.CameraLength
    #            microscope.attrs['binning'] = self.Bin
    #            microscope.attrs['max binning'] = self.maxBin
    #            microscope.attrs['dwell time'] = self.Dwell
    #            
                stage = f.create_group('stage')
                stage.attrs['type'] = self.stageName
    #            stage.attrs['position'] = self.getPosition()
    
                user = f.create_group('user')
    #            user.attrs['user name'] = self._iUser.GetValue()
                
                sample = f.create_group('sample')
    #            sample.attrs['sample name'] = self._iSample.GetValue()
            
        #Close the file
    
    def createDims(self, dataTop, pix):
        """ Create simple dim vectors. No pixel size is added because
        each image might have a different diffraction pixel size.
        
        """
        dim3 = dataTop.create_dataset('dim3', data=(0, 1))
        dim3.attrs['name'] = np.string_('kx')
        dim3.attrs['units'] = np.string_('pixels')
        
        dim2 = dataTop.create_dataset('dim2', data=(0, 1))
        dim2.attrs['name'] = np.string_('ky')
        dim2.attrs['units'] = np.string_('pixels')
        
        dim1 = dataTop.create_dataset('dim1', data=(0, 1))
        dim1.attrs['name'] = np.string_('')
        dim1.attrs['units'] = np.string_('')

        return 1
    
    def writeEMDdata(self):
        """ Write the cbed and meta data. This can be used to track
        how the user uses the software, data for ML learning, and 
        tracks the experiment for the user.
        
        """
        self.sb.SetStatusText('Writing cbed...')
        with h5py.File(self.filename, 'a') as f:
            data_root = f['data']
            g_cbed = data_root['cbed']
                
            dset_cbed = g_cbed['data']
            dset_params = g_cbed['params']
            dset_names = g_cbed['names']
            dset_position = g_cbed['position']
            
            # Resize the data sets
            cur_number = dset_cbed.shape[0]
            dset_cbed.resize(cur_number + 1, 0) #number, axis
            dset_params.resize(cur_number + 1, 0)
            dset_names.resize(cur_number + 1, 0)
            dset_position.resize(cur_number + 1, 0)
            
            # Save the cbed image
            dset_cbed[-1, :, :] = self.cbed
            
            # save the metadata
            print('WARNING: Check the params and data in the EMD.')
            dset_params[-1, :] = (self.axisXY[0], self.axisXY[1],
                                 self.cbedXY[0], self.cbedXY[1],
                                 self._microscope.Gun.HTValue,
                                 self.Proj.CameraLength,
                                 self.cbed_pixel_size,
                                 self.Ill.ConvergenceAngle,
                                 time.time())
            dset_names[-1, :] = (self._iUser.GetValue(),
                                self._iSample.GetValue())
            
            dset_position[-1, :] = self.getPosition()
        self.sb.SetStatusText('Idle...')
            

# Parse the arguments
parser = argparse.ArgumentParser(description='Measure tilt needed to bring a zone onto axis.')

# Allow user to set TeamStage or Compustage
parser.add_argument('--teamstage', '-t', action='store_const', const=True, default=False)
parser.add_argument('--compustage', '-c', action='store_const', const=True, default=True)

parser.add_argument('-v','--version', action='version', version='%(prog)s {}'.format(version))

parser.add_argument('-vb','--verbose', action='store_const', const=True, default=False)

args = parser.parse_args()

if args.teamstage:
    stagetype = 'teamstage'
else:
    stagetype = 'compustage'
    
print('Stage type = {}'.format(stagetype))
 
app = wx.App(False)
frame = TilterFrame(None, "NCEM: KSPace Tilter", stagetype = stagetype, verbose=args.verbose)
app.MainLoop()
