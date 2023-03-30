# A minimal script to test connecting to the microsocpe and TIA in STEM mode.
# An image is acquired, but not saved.

# Connect to FEI TEMScripting and TIA using COM 
from comtypes.client import CreateObject
from comtypes.safearray import safearray_as_ndarray # get data across COM barrier fast

def getPixelSize():
    """Get the pixel calibration from TIA
    """
    window1 = TIA.ActiveDisplayWindow()
    Im1 = window1.FindDisplay(window1.DisplayNames[0]) # returns an image display object
    unit1 = Im1.SpatialUnit # returns SpatialUnit object
    unitName = unit1.unitstring # returns a string (such as nm)
    calX = Im1.image.calibration.deltaX # returns the x calibration in meters
    calY = Im1.image.calibration.deltaY
    print('pixel size = {:03}, {:03}'.format(calX, calY))
    return((calX, calY))

_microscope = CreateObject('TEMScripting.Instrument')
print("Connected to microscope")
TIA = CreateObject('ESVision.Application')
print('Connected to TIA')

# Get microscope interfaces
Acq = _microscope.Acquisition
Proj = _microscope.Projection
Ill = _microscope.Illumination
Stage = _microscope.Stage 

detector0 = Acq.Detectors(0) #old microscopes might require square brackets [] instead of ()
# Add the detector
Acq.AddAcqDevice(detector0)

myStemAcqParams = Acq.Detectors.AcqParams
myStemAcqParams.Binning = 8
myStemAcqParams.ImageSize = 0 #_microscope.ACQIMAGESIZE_FULL
myStemAcqParams.DwellTime = 1e-6
Acq.Detectors.AcqParams = myStemAcqParams

print(detector0.Info.Name)

print('Acquire image')
acquiredImageSet = Acq.AcquireImages() 

# This quickly retrieves the data across the COM barrier
# and returns the data as a numpy ndarray
with safearray_as_ndarray:
    image = acquiredImageSet(0).AsSafeArray

getPixelSize()

print('convergence angle = {}'.format(Ill.ConvergenceAngle))
