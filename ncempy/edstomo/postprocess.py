import sys, os, shutil
from collections import OrderedDict
import numpy as np
import genfire

CropBounds = {'xmin': 0.10, 'xmax': 0.90,
              'ymin': 0.10, 'ymax': 0.90,
              'zmin': 0.10, 'zmax': 0.90}

def ReadGENFIRESignals(GENFIREDirectory, SignalNames):
    SignalDict = OrderedDict()

    for SigName in SignalNames:
        Sig = genfire.fileio.readMRC(SigName+'_reconstruction.mrc')
        SignalDict[SigName] = Sig

def SquareCropSignal(Signal, CropBounds = None):
    # x axis
    Signal[:,:, :int(CropBounds['xmin']*Signal.shape[2])] = 0 # Min for axis
    Signal[:,:, int(CropBounds['xmax']*Signal.shape[2]):] = 0 # Max for axis.

    # y axis
    Signal[:, :int(CropBounds['ymin']*Signal.shape[1]), :] = 0 # Min for axis
    Signal[:, int(CropBounds['ymax']*Signal.shape[1]):, :] = 0 # Max for axis.
    
    # z axis
    Signal[:int(CropBounds['zmin']*Signal.shape[0]), :, :] = 0 # Min for axis
    Signal[int(CropBounds['zmax']*Signal.shape[0]):, :, :] = 0 # Max for axis.

    return Signal
